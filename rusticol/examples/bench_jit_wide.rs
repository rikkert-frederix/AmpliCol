use std::env;
use std::fs;
use std::path::Path;
use std::time::Instant;

use symbolica::evaluate::JITCompiledEvaluator;
use symbolica::prelude::{
    BatchEvaluator, Complex, ExpressionEvaluator, JITCompilationSettings, Rational,
};

type SavedEvaluator = (
    bool,
    JITCompilationSettings,
    ExpressionEvaluator<Complex<Rational>>,
    Option<JITCompiledEvaluator<f64>>,
    Option<JITCompiledEvaluator<Complex<f64>>>,
);

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() != 6 {
        eprintln!("usage: bench_jit_wide STATE INPUT_LEN OUTPUT_LEN BATCH REPEATS");
        std::process::exit(2);
    }
    let path = Path::new(&args[1]);
    let input_len: usize = args[2].parse().expect("invalid input length");
    let output_len: usize = args[3].parse().expect("invalid output length");
    let batch: usize = args[4].parse().expect("invalid batch size");
    let repeats: usize = args[5].parse().expect("invalid repeat count");

    let bytes = fs::read(path).expect("could not read evaluator state");
    let ((_, settings, exact, _, saved), _) =
        bincode::decode_from_slice::<SavedEvaluator, _>(&bytes, bincode::config::standard())
            .expect("could not decode evaluator state");
    let mut scalar = saved.expect("state has no saved complex JIT evaluator");

    let scalar_compile_started = Instant::now();
    let mut fresh_scalar = exact
        .jit_compile::<Complex<f64>>(settings.clone())
        .expect("could not compile fresh scalar complex evaluator");
    let scalar_compile_s = scalar_compile_started.elapsed().as_secs_f64();

    let compile_started = Instant::now();
    let mut wide = exact
        .jit_compile::<Complex<wide::f64x2>>(settings)
        .expect("could not compile wide complex evaluator");
    let wide_compile_s = compile_started.elapsed().as_secs_f64();

    let params = vec![Complex::new(0.75, 0.125); batch * input_len];
    let mut scalar_out = vec![Complex::new(0.0, 0.0); batch * output_len];
    let mut fresh_scalar_out = scalar_out.clone();
    let mut wide_out = scalar_out.clone();
    let native_rows = batch.div_ceil(2);
    let mut native_params =
        vec![Complex::new(wide::f64x2::ZERO, wide::f64x2::ZERO); native_rows * input_len];
    for native_row in 0..native_rows {
        let first_row = native_row * 2;
        let second_row = usize::min(first_row + 1, batch - 1);
        for parameter in 0..input_len {
            let first = params[first_row * input_len + parameter];
            let second = params[second_row * input_len + parameter];
            native_params[native_row * input_len + parameter] = Complex::new(
                wide::f64x2::new([first.re, second.re]),
                wide::f64x2::new([first.im, second.im]),
            );
        }
    }
    let mut native_out =
        vec![Complex::new(wide::f64x2::ZERO, wide::f64x2::ZERO); native_rows * output_len];

    scalar
        .evaluate_batch(batch, &params, &mut scalar_out)
        .expect("scalar JIT warmup failed");
    fresh_scalar
        .evaluate_batch(batch, &params, &mut fresh_scalar_out)
        .expect("fresh scalar JIT warmup failed");
    wide.evaluate_batch(batch, &params, &mut wide_out)
        .expect("wide JIT warmup failed");
    for (output, input) in native_out
        .chunks_mut(output_len)
        .zip(native_params.chunks(input_len))
    {
        wide.evaluate(input, output);
    }

    let mut max_abs_diff = 0.0_f64;
    let mut max_fresh_scalar_diff = 0.0_f64;
    let mut max_rel_diff = 0.0_f64;
    let mut max_output_abs = 0.0_f64;
    let mut max_diff_index = 0usize;
    for (index, (scalar_value, wide_value)) in scalar_out.iter().zip(&wide_out).enumerate() {
        let fresh_scalar_value = fresh_scalar_out[index];
        max_fresh_scalar_diff = max_fresh_scalar_diff.max(
            (scalar_value.re - fresh_scalar_value.re)
                .hypot(scalar_value.im - fresh_scalar_value.im),
        );
        let abs_diff = ((scalar_value.re - wide_value.re).powi(2)
            + (scalar_value.im - wide_value.im).powi(2))
        .sqrt();
        let output_abs = scalar_value.re.hypot(scalar_value.im);
        if abs_diff > max_abs_diff {
            max_abs_diff = abs_diff;
            max_diff_index = index;
        }
        max_output_abs = max_output_abs.max(output_abs);
        max_rel_diff = max_rel_diff.max(abs_diff / output_abs.max(1.0e-300));
    }
    let reference_row = &wide_out[..output_len];
    let max_wide_lane_spread = wide_out
        .chunks(output_len)
        .skip(1)
        .flat_map(|row| row.iter().zip(reference_row))
        .map(|(a, b)| (a.re - b.re).hypot(a.im - b.im))
        .fold(0.0_f64, f64::max);
    let scalar_at_max = scalar_out[max_diff_index];
    let wide_at_max = wide_out[max_diff_index];
    let mut max_diff_by_lane = [0.0_f64; 4];
    for (row, (scalar_row, wide_row)) in scalar_out
        .chunks(output_len)
        .zip(wide_out.chunks(output_len))
        .enumerate()
    {
        for (scalar_value, wide_value) in scalar_row.iter().zip(wide_row) {
            max_diff_by_lane[row % 4] = max_diff_by_lane[row % 4]
                .max((scalar_value.re - wide_value.re).hypot(scalar_value.im - wide_value.im));
        }
    }

    let scalar_started = Instant::now();
    for _ in 0..repeats {
        scalar
            .evaluate_batch(batch, &params, &mut scalar_out)
            .expect("scalar JIT timing failed");
    }
    let scalar_s = scalar_started.elapsed().as_secs_f64();

    let wide_started = Instant::now();
    for _ in 0..repeats {
        wide.evaluate_batch(batch, &params, &mut wide_out)
            .expect("wide JIT timing failed");
    }
    let wide_s = wide_started.elapsed().as_secs_f64();

    let native_started = Instant::now();
    for _ in 0..repeats {
        for (output, input) in native_out
            .chunks_mut(output_len)
            .zip(native_params.chunks(input_len))
        {
            wide.evaluate(input, output);
        }
    }
    let native_s = native_started.elapsed().as_secs_f64();

    let mut max_native_diff = 0.0_f64;
    for (native_row, output) in native_out.chunks(output_len).enumerate() {
        for (column, value) in output.iter().enumerate() {
            for lane in 0..2 {
                let row = usize::min(native_row * 2 + lane, batch - 1);
                let scalar_value = scalar_out[row * output_len + column];
                max_native_diff = max_native_diff.max(
                    (scalar_value.re - value.re.as_array()[lane])
                        .hypot(scalar_value.im - value.im.as_array()[lane]),
                );
            }
        }
    }

    println!(
        "{{\"batch\":{batch},\"repeats\":{repeats},\"scalar_compile_s\":{scalar_compile_s},\"wide_compile_s\":{wide_compile_s},\"scalar_us_per_point\":{},\"wide_us_per_point\":{},\"native_us_per_point\":{},\"speedup\":{},\"native_speedup\":{},\"max_abs_diff\":{max_abs_diff},\"max_fresh_scalar_diff\":{max_fresh_scalar_diff},\"max_native_diff\":{max_native_diff},\"max_rel_diff\":{max_rel_diff},\"max_output_abs\":{max_output_abs},\"max_diff_index\":{max_diff_index},\"scalar_at_max\":[{},{}],\"wide_at_max\":[{},{}],\"max_wide_lane_spread\":{max_wide_lane_spread},\"max_diff_by_lane\":{:?}}}",
        scalar_s * 1.0e6 / (batch * repeats) as f64,
        wide_s * 1.0e6 / (batch * repeats) as f64,
        native_s * 1.0e6 / (batch * repeats) as f64,
        scalar_s / wide_s,
        scalar_s / native_s,
        scalar_at_max.re,
        scalar_at_max.im,
        wide_at_max.re,
        wide_at_max.im,
        max_diff_by_lane,
    );
}
