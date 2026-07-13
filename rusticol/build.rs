use std::path::Path;
use std::process::Command;

fn main() {
    let profile = std::env::var("PROFILE").unwrap_or_else(|_| "unknown".to_string());
    let target = std::env::var("TARGET").unwrap_or_else(|_| "unknown".to_string());
    println!("cargo:rustc-env=RUSTICOL_BUILD_PROFILE={profile}");
    println!("cargo:rustc-env=RUSTICOL_BUILD_TARGET={target}");
    if target.contains("apple-darwin") {
        add_gcc_runtime_rpath();
    }
}

fn add_gcc_runtime_rpath() {
    let Ok(output) = Command::new("gfortran")
        .arg("-print-file-name=libgcc_s.1.1.dylib")
        .output()
    else {
        return;
    };
    if !output.status.success() {
        return;
    }
    let path_text = String::from_utf8_lossy(&output.stdout);
    let path = Path::new(path_text.trim());
    let Some(parent) = path.parent() else {
        return;
    };
    let Ok(runtime_dir) = parent.canonicalize() else {
        return;
    };
    println!("cargo:rustc-link-arg=-Wl,-rpath,{}", runtime_dir.display());
}
