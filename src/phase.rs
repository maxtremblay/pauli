use num_complex::Complex;
use std::fmt;
use std::ops::{Mul, MulAssign};

type Scalar = i8;

/// A phase for Pauli operators.
///
/// Phases are limited to 1, -1, i and -i which
/// all have their own constructors.
///
/// # Example
///
/// ```
/// use pauli::Phase;
///
/// assert_eq!(Phase::minus_one() * Phase::i(), Phase::minus_i());
/// assert_eq!(Phase::i() * Phase::minus_i(), Phase::one());
/// ```
#[derive(Clone, Copy, PartialEq, Eq, Hash)]
pub struct Phase(Complex<Scalar>);

impl Phase {
    /// Create a new phase of 1.
    pub fn one() -> Self {
        Self(Complex::new(1, 0))
    }

    /// Create a new phase of -1.
    pub fn minus_one() -> Self {
        Self(Complex::new(-1, 0))
    }

    /// Create a new phase of i.
    pub fn i() -> Self {
        Self(Complex::new(0, 1))
    }

    /// Create a new phase of -i.
    pub fn minus_i() -> Self {
        Self(Complex::new(0, -1))
    }
}

impl Mul<Phase> for Phase {
    type Output = Phase;

    fn mul(self, other: Phase) -> Phase {
        Self(self.0 * other.0)
    }
}

impl MulAssign<Phase> for Phase {
    fn mul_assign(&mut self, other: Phase) {
        *self = *self * other;
    }
}

impl fmt::Debug for Phase {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        if *self == Phase::one() {
            write!(f, "1")
        } else if *self == Phase::minus_one() {
            write!(f, "-1")
        } else if *self == Phase::i() {
            write!(f, "i")
        } else if *self == Phase::minus_i() {
            write!(f, "-i")
        } else {
            // This should never happen.
            write!(f, "--")
        }
    }
}

impl fmt::Display for Phase {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        fmt::Debug::fmt(self, f)
    }
}
