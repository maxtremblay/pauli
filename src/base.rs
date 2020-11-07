pub use Pauli::{I, X, Y, Z};

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
#[derive(Debug, PartialEq, Eq, Clone, Copy, Hash)]
pub enum Pauli {
    I,
    X,
    Y,
    Z,
}

impl Pauli {
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
}

impl std::ops::Mul for Pauli {
    type Output = Self;

    fn mul(self, other: Self) -> Self {
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
