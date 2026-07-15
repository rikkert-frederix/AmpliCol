use libc::{c_char, c_double, c_int, size_t};
use rusticol::{NativeRuntime, RusticolError};
use std::any::Any;
use std::cell::RefCell;
use std::collections::BTreeMap;
use std::ffi::{CStr, CString};
use std::panic::{AssertUnwindSafe, catch_unwind};
use std::path::Path;
use std::ptr;
use std::slice;

pub const RUSTICOL_STATUS_OK: c_int = 0;
pub const RUSTICOL_STATUS_INVALID_ARGUMENT: c_int = 1;
pub const RUSTICOL_STATUS_BUFFER_TOO_SMALL: c_int = 2;
pub const RUSTICOL_STATUS_RUNTIME_ERROR: c_int = 3;
pub const RUSTICOL_STATUS_PANIC: c_int = 4;

pub struct RusticolRuntimeHandle {
    runtime: NativeRuntime,
}

struct AbiError {
    status: c_int,
    message: String,
}

type AbiResult<T> = Result<T, AbiError>;

impl From<RusticolError> for AbiError {
    fn from(error: RusticolError) -> Self {
        Self {
            status: RUSTICOL_STATUS_RUNTIME_ERROR,
            message: error.to_string(),
        }
    }
}

thread_local! {
    static LAST_ERROR: RefCell<CString> = RefCell::new(CString::new("").expect("empty CString"));
}

fn sanitize_c_string(message: impl AsRef<str>) -> CString {
    let sanitized = message.as_ref().replace('\0', "\\0");
    CString::new(sanitized)
        .unwrap_or_else(|_| CString::new("Rusticol error").expect("literal CString"))
}

fn set_last_error(message: impl AsRef<str>) {
    LAST_ERROR.with(|slot| *slot.borrow_mut() = sanitize_c_string(message));
}

fn guard(operation: impl FnOnce() -> AbiResult<()>) -> c_int {
    finish_guard(catch_unwind(AssertUnwindSafe(operation)))
}

fn finish_guard(result: Result<AbiResult<()>, Box<dyn Any + Send>>) -> c_int {
    match result {
        Ok(Ok(())) => RUSTICOL_STATUS_OK,
        Ok(Err(error)) => {
            set_last_error(&error.message);
            error.status
        }
        Err(payload) => {
            let message = payload
                .downcast_ref::<&str>()
                .copied()
                .or_else(|| payload.downcast_ref::<String>().map(String::as_str))
                .unwrap_or("unknown Rust panic");
            set_last_error(format!("Rusticol panic: {message}"));
            RUSTICOL_STATUS_PANIC
        }
    }
}

fn abi_error(status: c_int, message: impl Into<String>) -> AbiError {
    AbiError {
        status,
        message: message.into(),
    }
}

fn invalid(message: impl Into<String>) -> AbiError {
    abi_error(RUSTICOL_STATUS_INVALID_ARGUMENT, message)
}

fn buffer_too_small(message: impl Into<String>) -> AbiError {
    abi_error(RUSTICOL_STATUS_BUFFER_TOO_SMALL, message)
}

unsafe fn required_handle<'a>(
    handle: *const RusticolRuntimeHandle,
) -> AbiResult<&'a RusticolRuntimeHandle> {
    if handle.is_null() {
        return Err(invalid("Rusticol runtime handle is null"));
    }
    // SAFETY: The caller promises that a non-null handle was returned by rusticol_runtime_load
    // and remains alive for the duration of this call.
    Ok(unsafe { &*handle })
}

unsafe fn required_handle_mut<'a>(
    handle: *mut RusticolRuntimeHandle,
) -> AbiResult<&'a mut RusticolRuntimeHandle> {
    if handle.is_null() {
        return Err(invalid("Rusticol runtime handle is null"));
    }
    // SAFETY: The ABI documents handles as mutable and non-concurrently callable.
    Ok(unsafe { &mut *handle })
}

unsafe fn optional_c_string<'a>(value: *const c_char) -> AbiResult<Option<&'a str>> {
    if value.is_null() {
        return Ok(None);
    }
    // SAFETY: The caller supplies a NUL-terminated string for non-null pointers.
    let value = unsafe { CStr::from_ptr(value) }
        .to_str()
        .map_err(|error| invalid(format!("C string is not valid UTF-8: {error}")))?;
    Ok(Some(value))
}

unsafe fn required_c_string<'a>(value: *const c_char, name: &str) -> AbiResult<&'a str> {
    unsafe { optional_c_string(value) }?.ok_or_else(|| invalid(format!("{name} is null")))
}

unsafe fn write_size(value: usize, output: *mut size_t, name: &str) -> AbiResult<()> {
    if output.is_null() {
        return Err(invalid(format!("{name} is null")));
    }
    // SAFETY: The caller supplies writable storage for one size_t.
    unsafe { *output = value as size_t };
    Ok(())
}

unsafe fn write_i32(value: i32, output: *mut i32, name: &str) -> AbiResult<()> {
    if output.is_null() {
        return Err(invalid(format!("{name} is null")));
    }
    // SAFETY: The caller supplies writable storage for one i32.
    unsafe { *output = value };
    Ok(())
}

unsafe fn write_string(
    value: &str,
    buffer: *mut c_char,
    capacity: size_t,
    required: *mut size_t,
) -> AbiResult<()> {
    let bytes = value.as_bytes();
    let required_capacity = bytes
        .len()
        .checked_add(1)
        .ok_or_else(|| invalid("string length overflow"))?;
    unsafe { write_size(required_capacity, required, "required string capacity")? };
    if buffer.is_null() {
        if capacity == 0 {
            return Ok(());
        }
        return Err(invalid("string output buffer is null"));
    }
    if capacity < required_capacity {
        return Err(buffer_too_small(format!(
            "string output buffer has capacity {capacity}, requires {required_capacity}"
        )));
    }
    // SAFETY: Capacity was checked and the source does not overlap the caller's buffer.
    unsafe {
        ptr::copy_nonoverlapping(bytes.as_ptr(), buffer.cast::<u8>(), bytes.len());
        *buffer.add(bytes.len()) = 0;
    }
    Ok(())
}

unsafe fn read_selector_ids(
    values: *const *const c_char,
    count: size_t,
    name: &str,
) -> AbiResult<Option<Vec<String>>> {
    if count == 0 {
        return Ok(None);
    }
    if values.is_null() {
        return Err(invalid(format!("{name} array is null")));
    }
    // SAFETY: The caller supplies count string pointers.
    let values = unsafe { slice::from_raw_parts(values, count) };
    let mut result = Vec::with_capacity(count);
    for (index, value) in values.iter().enumerate() {
        result.push(unsafe { required_c_string(*value, &format!("{name}[{index}]")) }?.to_string());
    }
    Ok(Some(result))
}

unsafe fn read_f64_slice<'a>(
    values: *const c_double,
    count: size_t,
    name: &str,
) -> AbiResult<&'a [f64]> {
    if count == 0 {
        return Err(invalid(format!("{name} must not be empty")));
    }
    if values.is_null() {
        return Err(invalid(format!("{name} is null")));
    }
    // SAFETY: The caller supplies count readable f64 values.
    Ok(unsafe { slice::from_raw_parts(values, count) })
}

unsafe fn write_f64_slice(
    values: &[f64],
    output: *mut c_double,
    capacity: size_t,
    name: &str,
) -> AbiResult<()> {
    if output.is_null() {
        return Err(invalid(format!("{name} is null")));
    }
    if capacity < values.len() {
        return Err(buffer_too_small(format!(
            "{name} has capacity {capacity}, requires {}",
            values.len()
        )));
    }
    // SAFETY: Capacity was checked and the source cannot overlap the caller-owned output.
    unsafe { ptr::copy_nonoverlapping(values.as_ptr(), output, values.len()) };
    Ok(())
}

#[unsafe(no_mangle)]
pub extern "C" fn rusticol_abi_version() -> u32 {
    catch_unwind(|| NativeRuntime::ABI_VERSION).unwrap_or(0)
}

#[unsafe(no_mangle)]
pub extern "C" fn rusticol_last_error_message(
    buffer: *mut c_char,
    capacity: size_t,
    required: *mut size_t,
) -> c_int {
    match catch_unwind(AssertUnwindSafe(|| {
        LAST_ERROR.with(|slot| {
            let value = slot.borrow();
            // SAFETY: Pointer validation and copying are performed by write_string.
            unsafe {
                write_string(
                    value.to_str().unwrap_or("Rusticol error"),
                    buffer,
                    capacity,
                    required,
                )
            }
        })
    })) {
        Ok(Ok(())) => RUSTICOL_STATUS_OK,
        Ok(Err(error)) => {
            set_last_error(&error.message);
            RUSTICOL_STATUS_BUFFER_TOO_SMALL
        }
        Err(_) => RUSTICOL_STATUS_PANIC,
    }
}

#[unsafe(no_mangle)]
pub extern "C" fn rusticol_runtime_load(
    process_dir: *const c_char,
    process_key: *const c_char,
    model_parameters_path: *const c_char,
    output: *mut *mut RusticolRuntimeHandle,
) -> c_int {
    guard(|| {
        if output.is_null() {
            return Err(invalid("runtime output handle is null"));
        }
        // SAFETY: String pointers follow the ABI contract; output is checked above.
        let process_dir = unsafe { required_c_string(process_dir, "process_dir") }?;
        let process_key = unsafe { optional_c_string(process_key) }?;
        let model_parameters = unsafe { optional_c_string(model_parameters_path) }?;
        let runtime =
            NativeRuntime::load(process_dir, process_key, model_parameters.map(Path::new))?;
        let boxed = Box::new(RusticolRuntimeHandle { runtime });
        // SAFETY: output points to writable handle storage supplied by the caller.
        unsafe { *output = Box::into_raw(boxed) };
        Ok(())
    })
}

#[unsafe(no_mangle)]
pub extern "C" fn rusticol_runtime_free(handle: *mut RusticolRuntimeHandle) -> c_int {
    guard(|| {
        if !handle.is_null() {
            // SAFETY: The caller must free a live handle exactly once.
            unsafe { drop(Box::from_raw(handle)) };
        }
        Ok(())
    })
}

fn runtime_string(
    handle: *const RusticolRuntimeHandle,
    buffer: *mut c_char,
    capacity: size_t,
    required: *mut size_t,
    get: impl FnOnce(&NativeRuntime) -> AbiResult<String>,
) -> c_int {
    guard(|| {
        // SAFETY: Handle and output pointers are validated by helpers.
        let handle = unsafe { required_handle(handle) }?;
        let value = get(&handle.runtime)?;
        unsafe { write_string(&value, buffer, capacity, required) }
    })
}

#[unsafe(no_mangle)]
pub extern "C" fn rusticol_runtime_metadata_json(
    handle: *const RusticolRuntimeHandle,
    buffer: *mut c_char,
    capacity: size_t,
    required: *mut size_t,
) -> c_int {
    runtime_string(handle, buffer, capacity, required, |runtime| {
        Ok(runtime.metadata_json()?)
    })
}

#[unsafe(no_mangle)]
pub extern "C" fn rusticol_runtime_physics_json(
    handle: *const RusticolRuntimeHandle,
    buffer: *mut c_char,
    capacity: size_t,
    required: *mut size_t,
) -> c_int {
    runtime_string(handle, buffer, capacity, required, |runtime| {
        Ok(runtime.physics_json()?)
    })
}

#[unsafe(no_mangle)]
pub extern "C" fn rusticol_runtime_process(
    handle: *const RusticolRuntimeHandle,
    buffer: *mut c_char,
    capacity: size_t,
    required: *mut size_t,
) -> c_int {
    runtime_string(handle, buffer, capacity, required, |runtime| {
        Ok(runtime.metadata().process)
    })
}

#[unsafe(no_mangle)]
pub extern "C" fn rusticol_runtime_process_key(
    handle: *const RusticolRuntimeHandle,
    buffer: *mut c_char,
    capacity: size_t,
    required: *mut size_t,
) -> c_int {
    runtime_string(handle, buffer, capacity, required, |runtime| {
        Ok(runtime.metadata().process_key)
    })
}

#[unsafe(no_mangle)]
pub extern "C" fn rusticol_runtime_color_accuracy(
    handle: *const RusticolRuntimeHandle,
    buffer: *mut c_char,
    capacity: size_t,
    required: *mut size_t,
) -> c_int {
    runtime_string(handle, buffer, capacity, required, |runtime| {
        Ok(runtime.metadata().color_accuracy)
    })
}

#[unsafe(no_mangle)]
pub extern "C" fn rusticol_runtime_external_count(
    handle: *const RusticolRuntimeHandle,
    output: *mut size_t,
) -> c_int {
    guard(|| {
        // SAFETY: Helpers validate pointers.
        let handle = unsafe { required_handle(handle) }?;
        unsafe { write_size(handle.runtime.external_count(), output, "external count") }
    })
}

#[unsafe(no_mangle)]
pub extern "C" fn rusticol_runtime_external_pdg(
    handle: *const RusticolRuntimeHandle,
    index: size_t,
    output: *mut i32,
) -> c_int {
    guard(|| {
        // SAFETY: Helpers validate pointers.
        let handle = unsafe { required_handle(handle) }?;
        let particles = handle.runtime.external_particles()?;
        let particle = particles
            .get(index)
            .ok_or_else(|| invalid(format!("external particle index {index} is out of range")))?;
        unsafe { write_i32(particle.pdg, output, "external PDG output") }
    })
}

#[unsafe(no_mangle)]
pub extern "C" fn rusticol_runtime_helicity_count(
    handle: *const RusticolRuntimeHandle,
    output: *mut size_t,
) -> c_int {
    guard(|| {
        // SAFETY: Helpers validate pointers.
        let handle = unsafe { required_handle(handle) }?;
        unsafe { write_size(handle.runtime.helicities()?.len(), output, "helicity count") }
    })
}

#[unsafe(no_mangle)]
pub extern "C" fn rusticol_runtime_helicity_id(
    handle: *const RusticolRuntimeHandle,
    index: size_t,
    buffer: *mut c_char,
    capacity: size_t,
    required: *mut size_t,
) -> c_int {
    runtime_string(handle, buffer, capacity, required, |runtime| {
        runtime
            .helicities()?
            .get(index)
            .map(|item| item.id.clone())
            .ok_or_else(|| invalid(format!("helicity index {index} is out of range")))
    })
}

#[unsafe(no_mangle)]
pub extern "C" fn rusticol_runtime_helicity_vector(
    handle: *const RusticolRuntimeHandle,
    index: size_t,
    output: *mut i32,
    capacity: size_t,
    required: *mut size_t,
) -> c_int {
    guard(|| {
        // SAFETY: Helpers validate pointers.
        let handle = unsafe { required_handle(handle) }?;
        let helicities = handle.runtime.helicities()?;
        let item = helicities
            .get(index)
            .ok_or_else(|| invalid(format!("helicity index {index} is out of range")))?;
        unsafe { write_size(item.helicities.len(), required, "helicity vector length")? };
        if output.is_null() {
            if capacity == 0 {
                return Ok(());
            }
            return Err(invalid("helicity vector output is null"));
        }
        if capacity < item.helicities.len() {
            return Err(buffer_too_small(format!(
                "helicity vector capacity {capacity} is smaller than {}",
                item.helicities.len()
            )));
        }
        // SAFETY: Capacity was checked.
        unsafe {
            ptr::copy_nonoverlapping(item.helicities.as_ptr(), output, item.helicities.len())
        };
        Ok(())
    })
}

#[unsafe(no_mangle)]
pub extern "C" fn rusticol_runtime_color_count(
    handle: *const RusticolRuntimeHandle,
    output: *mut size_t,
) -> c_int {
    guard(|| {
        // SAFETY: Helpers validate pointers.
        let handle = unsafe { required_handle(handle) }?;
        unsafe {
            write_size(
                handle.runtime.color_components()?.len(),
                output,
                "color count",
            )
        }
    })
}

#[unsafe(no_mangle)]
pub extern "C" fn rusticol_runtime_color_id(
    handle: *const RusticolRuntimeHandle,
    index: size_t,
    buffer: *mut c_char,
    capacity: size_t,
    required: *mut size_t,
) -> c_int {
    runtime_string(handle, buffer, capacity, required, |runtime| {
        runtime
            .color_components()?
            .get(index)
            .map(|item| item.id.clone())
            .ok_or_else(|| invalid(format!("color index {index} is out of range")))
    })
}

#[unsafe(no_mangle)]
pub extern "C" fn rusticol_runtime_color_kind(
    handle: *const RusticolRuntimeHandle,
    index: size_t,
    buffer: *mut c_char,
    capacity: size_t,
    required: *mut size_t,
) -> c_int {
    runtime_string(handle, buffer, capacity, required, |runtime| {
        runtime
            .color_components()?
            .get(index)
            .map(|item| item.kind.clone())
            .ok_or_else(|| invalid(format!("color index {index} is out of range")))
    })
}

#[unsafe(no_mangle)]
pub extern "C" fn rusticol_runtime_color_word(
    handle: *const RusticolRuntimeHandle,
    index: size_t,
    output: *mut size_t,
    capacity: size_t,
    required: *mut size_t,
) -> c_int {
    guard(|| {
        // SAFETY: Helpers validate pointers.
        let handle = unsafe { required_handle(handle) }?;
        let colors = handle.runtime.color_components()?;
        let item = colors
            .get(index)
            .ok_or_else(|| invalid(format!("color index {index} is out of range")))?;
        unsafe { write_size(item.word.len(), required, "color word length")? };
        if output.is_null() {
            if capacity == 0 {
                return Ok(());
            }
            return Err(invalid("color word output is null"));
        }
        if capacity < item.word.len() {
            return Err(buffer_too_small(format!(
                "color word capacity {capacity} is smaller than {}",
                item.word.len()
            )));
        }
        // SAFETY: Capacity was checked and usize matches C size_t.
        unsafe { ptr::copy_nonoverlapping(item.word.as_ptr(), output, item.word.len()) };
        Ok(())
    })
}

#[unsafe(no_mangle)]
pub extern "C" fn rusticol_runtime_model_parameter_count(
    handle: *const RusticolRuntimeHandle,
    output: *mut size_t,
) -> c_int {
    guard(|| {
        // SAFETY: Helpers validate pointers.
        let handle = unsafe { required_handle(handle) }?;
        unsafe {
            write_size(
                handle.runtime.model_parameters()?.len(),
                output,
                "model parameter count",
            )
        }
    })
}

#[unsafe(no_mangle)]
pub extern "C" fn rusticol_runtime_model_parameter_name(
    handle: *const RusticolRuntimeHandle,
    index: size_t,
    buffer: *mut c_char,
    capacity: size_t,
    required: *mut size_t,
) -> c_int {
    runtime_string(handle, buffer, capacity, required, |runtime| {
        runtime
            .model_parameters()?
            .get(index)
            .map(|item| item.name.clone())
            .ok_or_else(|| invalid(format!("model parameter index {index} is out of range")))
    })
}

#[unsafe(no_mangle)]
pub extern "C" fn rusticol_runtime_resolved_shape(
    handle: *const RusticolRuntimeHandle,
    helicity_ids: *const *const c_char,
    helicity_count: size_t,
    color_ids: *const *const c_char,
    color_count: size_t,
    output_helicity_count: *mut size_t,
    output_color_count: *mut size_t,
) -> c_int {
    guard(|| {
        // SAFETY: Helpers validate pointers and arrays.
        let handle = unsafe { required_handle(handle) }?;
        let helicities =
            unsafe { read_selector_ids(helicity_ids, helicity_count, "helicity ids") }?;
        let colors = unsafe { read_selector_ids(color_ids, color_count, "color ids") }?;
        let (helicity_count, color_count) = handle
            .runtime
            .resolved_shape(helicities.as_deref(), colors.as_deref())?;
        unsafe {
            write_size(
                helicity_count,
                output_helicity_count,
                "resolved helicity count",
            )?;
            write_size(color_count, output_color_count, "resolved color count")
        }
    })
}

#[unsafe(no_mangle)]
pub extern "C" fn rusticol_runtime_evaluate_f64(
    handle: *mut RusticolRuntimeHandle,
    momenta: *const c_double,
    momentum_count: size_t,
    point_count: size_t,
    output: *mut c_double,
    output_capacity: size_t,
) -> c_int {
    guard(|| {
        // SAFETY: Helpers validate all pointers.
        let handle = unsafe { required_handle_mut(handle) }?;
        let momenta = unsafe { read_f64_slice(momenta, momentum_count, "momenta") }?;
        let values = handle.runtime.evaluate_f64(momenta, point_count)?;
        unsafe { write_f64_slice(&values, output, output_capacity, "total output") }
    })
}

#[unsafe(no_mangle)]
pub extern "C" fn rusticol_runtime_evaluate_resolved_f64(
    handle: *mut RusticolRuntimeHandle,
    momenta: *const c_double,
    momentum_count: size_t,
    point_count: size_t,
    helicity_ids: *const *const c_char,
    helicity_count: size_t,
    color_ids: *const *const c_char,
    color_count: size_t,
    output: *mut c_double,
    output_capacity: size_t,
    output_helicity_count: *mut size_t,
    output_color_count: *mut size_t,
) -> c_int {
    guard(|| {
        // SAFETY: Helpers validate all pointers and arrays.
        let handle = unsafe { required_handle_mut(handle) }?;
        let momenta = unsafe { read_f64_slice(momenta, momentum_count, "momenta") }?;
        let helicities =
            unsafe { read_selector_ids(helicity_ids, helicity_count, "helicity ids") }?;
        let colors = unsafe { read_selector_ids(color_ids, color_count, "color ids") }?;
        let resolved = handle.runtime.evaluate_resolved_f64(
            momenta,
            point_count,
            helicities.as_deref(),
            colors.as_deref(),
        )?;
        let (_, helicity_count, color_count) = resolved.shape();
        unsafe {
            write_size(
                helicity_count,
                output_helicity_count,
                "resolved helicity count",
            )?;
            write_size(color_count, output_color_count, "resolved color count")?;
            write_f64_slice(&resolved.values, output, output_capacity, "resolved output")
        }
    })
}

#[unsafe(no_mangle)]
pub extern "C" fn rusticol_runtime_set_model_parameters(
    handle: *mut RusticolRuntimeHandle,
    names: *const *const c_char,
    real: *const c_double,
    imaginary: *const c_double,
    count: size_t,
) -> c_int {
    guard(|| {
        if count == 0 {
            return Err(invalid("model parameter update must not be empty"));
        }
        // SAFETY: Helpers and explicit checks validate all pointers.
        let handle = unsafe { required_handle_mut(handle) }?;
        let names = unsafe { read_selector_ids(names, count, "model parameter names") }?
            .expect("positive count returns names");
        let real = unsafe { read_f64_slice(real, count, "model parameter real values") }?;
        let imaginary =
            unsafe { read_f64_slice(imaginary, count, "model parameter imaginary values") }?;
        let mut values = BTreeMap::new();
        for index in 0..count {
            if values
                .insert(names[index].clone(), (real[index], imaginary[index]))
                .is_some()
            {
                return Err(invalid(format!(
                    "duplicate model parameter update {:?}",
                    names[index]
                )));
            }
        }
        Ok(handle.runtime.set_model_parameters(&values)?)
    })
}

#[unsafe(no_mangle)]
pub extern "C" fn rusticol_runtime_set_model_parameter(
    handle: *mut RusticolRuntimeHandle,
    name: *const c_char,
    real: c_double,
    imaginary: c_double,
) -> c_int {
    guard(|| {
        // SAFETY: Helpers validate pointers.
        let handle = unsafe { required_handle_mut(handle) }?;
        let name = unsafe { required_c_string(name, "model parameter name") }?;
        Ok(handle.runtime.set_model_parameter(name, real, imaginary)?)
    })
}

#[unsafe(no_mangle)]
pub extern "C" fn rusticol_runtime_set_model_parameters_json(
    handle: *mut RusticolRuntimeHandle,
    path: *const c_char,
) -> c_int {
    guard(|| {
        // SAFETY: Helpers validate pointers.
        let handle = unsafe { required_handle_mut(handle) }?;
        let path = unsafe { required_c_string(path, "model parameter JSON path") }?;
        Ok(handle.runtime.set_model_parameters_json(Path::new(path))?)
    })
}

#[unsafe(no_mangle)]
pub extern "C" fn rusticol_runtime_mute_warnings(
    handle: *mut RusticolRuntimeHandle,
    muted: c_int,
) -> c_int {
    guard(|| {
        // SAFETY: Helper validates the handle.
        let handle = unsafe { required_handle_mut(handle) }?;
        if muted == 0 {
            handle.runtime.unmute_warnings();
        } else {
            handle.runtime.mute_warnings();
        }
        Ok(())
    })
}

#[unsafe(no_mangle)]
pub extern "C" fn rusticol_runtime_take_warnings_json(
    handle: *mut RusticolRuntimeHandle,
    buffer: *mut c_char,
    capacity: size_t,
    required: *mut size_t,
) -> c_int {
    guard(|| {
        // SAFETY: Helpers validate pointers.
        let handle = unsafe { required_handle_mut(handle) }?;
        let warnings = handle.runtime.take_warnings();
        let json = serde_json::to_string(&warnings)
            .map_err(|error| invalid(format!("could not serialize warnings: {error}")))?;
        unsafe { write_string(&json, buffer, capacity, required) }
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn every_abi_error_is_thread_local_and_nul_terminated() {
        set_last_error("one\0two");
        LAST_ERROR.with(|slot| {
            assert_eq!(slot.borrow().to_str().unwrap(), "one\\0two");
        });
    }

    #[test]
    fn null_handle_is_reported_without_panicking() {
        let mut count = 0;
        let status = rusticol_runtime_external_count(ptr::null(), &mut count);
        assert_eq!(status, RUSTICOL_STATUS_INVALID_ARGUMENT);
    }

    #[test]
    fn short_string_buffer_has_a_distinct_status() {
        set_last_error("a message larger than one byte");
        let mut byte = 0_i8;
        let mut required = 0;
        let status = rusticol_last_error_message(&mut byte, 1, &mut required);
        assert_eq!(status, RUSTICOL_STATUS_BUFFER_TOO_SMALL);
        assert!(required > 1);
    }

    #[test]
    fn panic_payloads_map_to_the_panic_status() {
        let payload: Box<dyn Any + Send> = Box::new(String::from("contained test panic"));
        let status = finish_guard(Err(payload));
        assert_eq!(status, RUSTICOL_STATUS_PANIC);
        LAST_ERROR.with(|slot| {
            assert!(
                slot.borrow()
                    .to_str()
                    .unwrap()
                    .contains("contained test panic")
            );
        });
    }
}
