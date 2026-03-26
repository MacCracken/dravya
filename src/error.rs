use thiserror::Error;

#[derive(Debug, Error)]
#[non_exhaustive]
pub enum DravyaError {
    #[error("invalid material: {0}")]
    InvalidMaterial(String),
    #[error("invalid stress: {0}")]
    InvalidStress(String),
    #[error("invalid strain: {0}")]
    InvalidStrain(String),
    #[error("yield exceeded: {0}")]
    YieldExceeded(String),
    #[error("computation error: {0}")]
    ComputationError(String),
}

pub type Result<T> = std::result::Result<T, DravyaError>;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn error_display() {
        let e = DravyaError::YieldExceeded("steel failed".into());
        assert!(e.to_string().contains("steel failed"));
    }
}
