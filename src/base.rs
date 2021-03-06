use std::fmt::{self, Display};
use std::ops::Mul;
use Pauli::{I, X, Y, Z};
use serde::{Serialize, Deserialize};

/// A single qubit Pauli operator without a phase.
///
/// These operators form a multiplicative group
/// and follow the usual commutation and anti-commutation relations.
///
/// # Example
///
/// ```
/// use pauli::{I, X, Y, Z};
///
/// assert_eq!(X * Y, Z);
/// assert!(X.commutes_with(I));
/// assert!(Y.anticommutes_with(Z));
/// ```
#[derive(Debug, PartialEq, Eq, Clone, Copy, Hash, Serialize, Deserialize)]
pub enum Pauli {
    I,
    X,
    Y,
    Z,
}

impl Pauli {
    /// Checks if the operator anti-commutes with the other operator.
    ///
    /// # Example
    ///
    /// ```
    /// use pauli::{I, X, Y, Z};
    ///
    /// assert!(X.anticommutes_with(Y));
    /// assert!(Y.anticommutes_with(Z));
    /// assert!(Z.anticommutes_with(X));
    /// ```
    pub fn anticommutes_with(self, other: Self) -> bool {
        !self.commutes_with(other)
    }

    /// Checks if the operator commutes with the other operator.
    ///
    /// # Example
    ///
    /// ```
    /// use pauli::{I, X, Y, Z};
    ///
    /// assert!(I.commutes_with(X));
    /// assert!(Y.commutes_with(Y));
    /// assert!(Z.commutes_with(I));
    /// ```
    pub fn commutes_with(self, other: Self) -> bool {
        self == I || other == I || self == other
    }

    /// Checks if the operator is not the identity.
    pub fn is_non_trivial(self) -> bool {
        self != I
    }

    /// Checks if the operator is the identity.
    pub fn is_trivial(self) -> bool {
        self == I
    }
}

impl Mul<Pauli> for Pauli {
    type Output = Pauli;

    fn mul(self, other: Pauli) -> Pauli {
        match (self, other) {
            (I, p) => p,
            (p, q) if p == q => I,
            (X, Y) => Z,
            (Y, Z) => X,
            (Z, X) => Y,
            (p, q) => q * p,
        }
    }
}

impl Display for Pauli {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            I => write!(f, "I"),
            X => write!(f, "X"),
            Y => write!(f, "Y"),
            Z => write!(f, "Z"),
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn commutations() {
        assert!(I.commutes_with(I));
        assert!(I.commutes_with(X));
        assert!(I.commutes_with(Y));
        assert!(I.commutes_with(Z));

        assert!(X.commutes_with(I));
        assert!(X.commutes_with(X));
        assert!(X.anticommutes_with(Y));
        assert!(X.anticommutes_with(Z));

        assert!(Y.commutes_with(I));
        assert!(Y.anticommutes_with(X));
        assert!(Y.commutes_with(Y));
        assert!(Y.anticommutes_with(Z));

        assert!(Z.commutes_with(I));
        assert!(Z.anticommutes_with(X));
        assert!(Z.anticommutes_with(Y));
        assert!(Z.commutes_with(Z));
    }

    #[test]
    fn multiplications() {
        assert_eq!(I * I, I);
        assert_eq!(I * X, X);
        assert_eq!(I * Y, Y);
        assert_eq!(I * Z, Z);

        assert_eq!(X * I, X);
        assert_eq!(X * X, I);
        assert_eq!(X * Y, Z);
        assert_eq!(X * Z, Y);

        assert_eq!(Y * I, Y);
        assert_eq!(Y * X, Z);
        assert_eq!(Y * Y, I);
        assert_eq!(Y * Z, X);

        assert_eq!(Z * I, Z);
        assert_eq!(Z * X, Y);
        assert_eq!(Z * Y, X);
        assert_eq!(Z * Z, I);
    }
}
