use thiserror::Error;

#[derive(Debug, Clone, Error)]
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
    #[error("division by zero: {0}")]
    DivisionByZero(&'static str),
    #[error("invalid parameter: {name} = {value} ({reason})")]
    InvalidParameter {
        name: &'static str,
        value: f64,
        reason: &'static str,
    },
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

    #[test]
    fn division_by_zero_display() {
        let e = DravyaError::DivisionByZero("zero modulus");
        assert_eq!(e.to_string(), "division by zero: zero modulus");
    }

    #[test]
    fn invalid_parameter_display() {
        let e = DravyaError::InvalidParameter {
            name: "poisson_ratio",
            value: -2.0,
            reason: "must be between -1 and 0.5",
        };
        assert!(e.to_string().contains("poisson_ratio"));
        assert!(e.to_string().contains("-2"));
    }

    #[test]
    fn error_is_clone() {
        let e = DravyaError::ComputationError("test".into());
        let e2 = e.clone();
        assert_eq!(e.to_string(), e2.to_string());
    }
}
