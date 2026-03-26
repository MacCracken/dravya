/// Initialize tracing subscriber for applications using dravya.
///
/// Uses the `DRAVYA_LOG` environment variable for filtering (defaults to `warn`).
/// Safe to call multiple times — subsequent calls are no-ops.
///
/// Note: this is a convenience for binaries. Library consumers should
/// configure their own subscriber; dravya emits standard `tracing` events.
pub fn init() {
    use tracing_subscriber::EnvFilter;
    let filter = EnvFilter::try_from_env("DRAVYA_LOG").unwrap_or_else(|_| EnvFilter::new("warn"));
    let _ = tracing_subscriber::fmt().with_env_filter(filter).try_init();
}
