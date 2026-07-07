fn main() {
    let profile = std::env::var("PROFILE").unwrap_or_else(|_| "unknown".to_string());
    let target = std::env::var("TARGET").unwrap_or_else(|_| "unknown".to_string());
    println!("cargo:rustc-env=RUSTICOL_BUILD_PROFILE={profile}");
    println!("cargo:rustc-env=RUSTICOL_BUILD_TARGET={target}");
}
