use crate::Pauli;
use sprs::vec::{NnzEither, SparseIterTools, VectorIterator};
use sprs::CsVec;
use std::error::Error;
use std::fmt;
use std::ops::Mul;
use Pauli::{X, Z};

/// A Pauli operator optimized for sparse operations.
///
/// A Pauli operator is a string of Paulis
/// such as `IXIX` or `XIYIZ`.
/// However,
/// we usually only care about the non-identity positions
/// and we refer to the previous as operators as `X_1 X_3`
/// and `X_0 Y_2 Z_4`.
#[derive(Debug, PartialEq, Eq, Clone, Hash)]
pub struct PauliOperator {
    paulis: CsVec<Pauli>,
}

impl PauliOperator {
    /// Builds a new Pauli Operator.
    ///
    /// To build an operator,
    /// we specify the length,
    /// the position of non-identity elements
    /// and their values.
    ///
    /// # Exemple
    ///
    /// This creates the `XIYIZ` operator.
    ///
    /// ```
    /// # use pauli::PauliOperator;
    /// # use pauli::{X, Y, Z};
    /// let operator = PauliOperator::new(5, vec![0, 2, 4], vec![X, Y, Z]);
    /// ```
    ///
    /// # Panic
    ///
    /// Panics if a position is greater or equal to the length or if
    /// the number of positions and Paulis are different.
    pub fn new(length: usize, positions: Vec<usize>, paulis: Vec<Pauli>) -> Self {
        Self::try_new(length, positions, paulis).expect("invalid operator")
    }

    /// Builds a new Pauli Operator or returns an error
    /// if either a position is greater or equal to the length or if
    /// the numbers of positions and Paulis are different.
    ///
    /// # Exemple
    /// This creates the `XIYIZ` operator.  
    /// ```
    /// # use pauli::PauliOperator;
    /// # use pauli::{X, Y, Z};
    /// let operator = PauliOperator::try_new(5, vec![0, 2, 4], vec![X, Y, Z]);
    /// assert!(operator.is_ok());
    /// ```
    pub fn try_new(
        length: usize,
        positions: Vec<usize>,
        paulis: Vec<Pauli>,
    ) -> Result<Self, PauliError> {
        if positions.len() != paulis.len() {
            Err(PauliError::IncompatibleLength(
                positions.len(),
                paulis.len(),
            ))
        } else if let Some(pos) = positions.iter().find(|pos| **pos >= length) {
            Err(PauliError::OutOfBound(*pos, length))
        } else {
            Ok(Self {
                paulis: CsVec::new(length, positions, paulis),
            })
        }
    }

    /// Creates a Pauli operator of zero length.
    pub fn empty() -> Self {
        Self {
            paulis: CsVec::empty(0),
        }
    }

    /// Checks if two operators commute.
    ///
    /// If an operator is smaller than the other,
    /// it is padded with identities.
    ///
    /// # Example
    ///
    /// ```
    /// # use pauli::PauliOperator;
    /// # use pauli::{X, Y, Z};
    /// let op1 = PauliOperator::new(5, vec![1, 2, 3], vec![X, Y, Z]);
    /// let op2 = PauliOperator::new(5, vec![2, 3, 4], vec![X, X, X]);
    /// let op3 = PauliOperator::new(5, vec![0, 1], vec![Z, Z]);
    ///
    /// assert!(op1.commutes_with(&op2));
    /// assert!(!op1.commutes_with(&op3));
    /// ```
    pub fn commutes_with(&self, other: &Self) -> bool {
        self.iter()
            .nnz_zip(other.iter())
            .filter(|(_, pauli, other_pauli)| pauli.anticommutes_with(**other_pauli))
            .count()
            % 2
            == 0
    }

    /// Checks if two operators anticommute.
    ///
    /// If an operator is smaller than the other,
    /// it is padded with identities.
    ///
    /// # Example
    ///
    /// ```
    /// # use pauli::PauliOperator;
    /// # use pauli::{X, Y, Z};
    /// let op1 = PauliOperator::new(5, vec![1, 2, 3], vec![X, Y, Z]);
    /// let op2 = PauliOperator::new(5, vec![2, 3, 4], vec![X, X, X]);
    /// let op3 = PauliOperator::new(5, vec![0, 1], vec![Z, Z]);
    ///
    /// assert!(!op1.anticommutes_with(&op2));
    /// assert!(op1.anticommutes_with(&op3));
    /// ```
    pub fn anticommutes_with(&self, other: &Self) -> bool {
        !self.commutes_with(other)
    }

    /// Returns an iterator over pairs of positions and Paulis.
    ///
    /// # Example
    ///
    /// ```
    /// # use pauli::{PauliOperator, X, Y, Z};
    /// let operator = PauliOperator::new(5, vec![0, 2, 4], vec![X, Y, Z]);
    /// let mut iter = operator.iter();
    ///
    /// assert_eq!(iter.next(), Some((0, &X)));
    /// assert_eq!(iter.next(), Some((2, &Y)));
    /// assert_eq!(iter.next(), Some((4, &Z)));
    /// assert_eq!(iter.next(), None);
    /// ```
    pub fn iter(&self) -> VectorIterator<Pauli, usize> {
        self.paulis.iter()
    }

    /// Returns the Pauli at the given position
    /// or None if the position is out of bound.
    ///
    /// # Example
    ///
    /// ```
    /// # use pauli::{PauliOperator, I, X, Y, Z};
    /// let operator = PauliOperator::new(5, vec![0, 2, 4], vec![X, Y, Z]);
    ///
    /// assert_eq!(operator.get(0), Some(X));
    /// assert_eq!(operator.get(1), Some(I));
    /// assert_eq!(operator.get(2), Some(Y));
    /// assert_eq!(operator.get(10), None);
    /// ```
    pub fn get(&self, position: usize) -> Option<Pauli> {
        self.paulis.get(position).cloned().or_else(|| {
            if position < self.len() {
                Some(Pauli::I)
            } else {
                None
            }
        })
    }

    /// Returns the length of the operator.
    pub fn len(&self) -> usize {
        self.paulis.dim()
    }

    /// Returns the number of non identity elements.
    pub fn weight(&self) -> usize {
        self.paulis.nnz()
    }

    /// Returns a slice of the positions
    /// where the element is not identity.
    ///
    /// # Example
    ///
    /// ```
    /// # use pauli::{PauliOperator, X, Y, Z};
    /// let operator = PauliOperator::new(5, vec![0, 2, 4], vec![X, Y, Z]);
    ///
    /// assert_eq!(operator.non_trivial_positions(), &[0, 2, 4]);
    /// ```
    pub fn non_trivial_positions(&self) -> &[usize] {
        self.paulis.indices()
    }

    /// Returns two operators such that there product is the
    /// original operator and the first contains only Xs and
    /// the second only Zs.
    ///
    /// # Example
    ///
    /// ```
    /// # use pauli::{PauliOperator, X, Y, Z};
    /// let operator = PauliOperator::new(5, vec![0, 2, 4], vec![X, Y, Z]);
    /// let (x_operator, z_operator) = operator.partition_x_and_z();
    ///
    /// assert_eq!(x_operator, PauliOperator::new(5, vec![0, 2], vec![X, X]));
    /// assert_eq!(z_operator, PauliOperator::new(5, vec![2, 4], vec![Z, Z]));
    /// ```
    pub fn partition_x_and_z(&self) -> (Self, Self) {
        (self.x_part(), self.z_part())
    }

    /// Returns the X part of the operator.
    ///
    /// # Example
    ///
    /// ```
    /// # use pauli::{PauliOperator, X, Y, Z};
    /// let operator = PauliOperator::new(5, vec![0, 2, 4], vec![X, Y, Z]);
    /// let x_operator = operator.x_part();
    ///
    /// assert_eq!(x_operator, PauliOperator::new(5, vec![0, 2], vec![X, X]));
    /// ```
    pub fn x_part(&self) -> Self {
        let (positions, paulis) = self
            .iter()
            .filter_map(|(position, pauli)| {
                if *pauli != Z {
                    Some((position, X))
                } else {
                    None
                }
            })
            .unzip();
        Self::new(self.len(), positions, paulis)
    }

    /// Returns the Z part of the operator.
    ///
    /// # Example
    ///
    /// ```
    /// # use pauli::{PauliOperator, X, Y, Z};
    /// let operator = PauliOperator::new(5, vec![0, 2, 4], vec![X, Y, Z]);
    /// let z_operator = operator.z_part();
    ///
    /// assert_eq!(z_operator, PauliOperator::new(5, vec![2, 4], vec![Z, Z]));
    /// ```
    pub fn z_part(&self) -> Self {
        let (positions, paulis) = self
            .iter()
            .filter_map(|(position, pauli)| {
                if *pauli != X {
                    Some((position, Z))
                } else {
                    None
                }
            })
            .unzip();
        Self::new(self.len(), positions, paulis)
    }

    /// Returns the element-wise product of two operators
    /// or an Error if they have different lengths.
    ///
    /// For a panicking version, use the `*` operator.
    ///
    /// # Example
    ///
    /// ```
    /// # use pauli::PauliOperator;
    /// # use pauli::{X, Y, Z};
    /// let op1 = PauliOperator::new(5, vec![1, 2, 3], vec![X, Y, Z]);
    /// let op2 = PauliOperator::new(5, vec![2, 3, 4], vec![Y, X, Z]);
    ///
    /// let product = PauliOperator::new(5, vec![1, 3, 4], vec![X, Y, Z]);
    ///
    /// assert_eq!(op1.multiply_with(&op2), Ok(product))
    /// ```
    pub fn multiply_with(&self, other: &Self) -> Result<Self, PauliError> {
        if self.len() != other.len() {
            Err(PauliError::IncompatibleLength(self.len(), other.len()))
        } else {
            let (positions, paulis) = self
                .iter()
                .nnz_or_zip(other.iter())
                .map(|values| match values {
                    NnzEither::Left((position, &pauli)) => (position, pauli),
                    NnzEither::Right((position, &pauli)) => (position, pauli),
                    NnzEither::Both((position, &p0, &p1)) => (position, p0 * p1),
                })
                .filter(|(_, pauli)| pauli.is_non_trivial())
                .unzip();
            Ok(PauliOperator::new(self.len(), positions, paulis))
        }
    }

    /// Converts a PauliOperator to a Vec of its non trivial positions
    /// consumming the operator.
    ///
    /// # Example
    ///
    /// ```
    /// # use pauli::PauliOperator;
    /// # use pauli::{X, Y, Z};
    /// let operator = PauliOperator::new(5, vec![1, 2, 3], vec![X, Y, Z]);
    /// let positions = operator.into_raw_positions();
    /// assert_eq!(positions, vec![1, 2, 3]);
    /// ```
    pub fn into_raw_positions(self) -> Vec<usize> {
        self.paulis.into_raw_storage().0
    }

    /// Converts a PauliOperator to a Vec of the Paulis
    /// consumming the operator.
    ///
    /// # Example
    ///
    /// ```
    /// # use pauli::PauliOperator;
    /// # use pauli::{X, Y, Z};
    /// let operator = PauliOperator::new(5, vec![1, 2, 3], vec![X, Y, Z]);
    /// let paulis = operator.into_raw_paulis();
    /// assert_eq!(paulis, vec![X, Y, Z]);
    /// ```
    pub fn into_raw_paulis(self) -> Vec<Pauli> {
        self.paulis.into_raw_storage().1
    }

    /// Converts a PauliOperator to a Vec of the positions and
    /// a Vec of Paulis consumming the operator.
    ///
    /// # Example
    ///
    /// ```
    /// # use pauli::PauliOperator;
    /// # use pauli::{X, Y, Z};
    /// let operator = PauliOperator::new(5, vec![1, 2, 3], vec![X, Y, Z]);
    /// let (positions, paulis) = operator.into_raw();
    /// assert_eq!(positions, vec![1, 2, 3]);
    /// assert_eq!(paulis, vec![X, Y, Z]);
    /// ```
    pub fn into_raw(self) -> (Vec<usize>, Vec<Pauli>) {
        self.paulis.into_raw_storage()
    }
}

impl<'a> Mul<&'a PauliOperator> for &'a PauliOperator {
    type Output = PauliOperator;

    fn mul(self, other: Self) -> Self::Output {
        self.multiply_with(other).unwrap()
    }
}

impl fmt::Display for PauliOperator {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "[")?;
        for (idx, (position, pauli)) in self.iter().enumerate() {
            write!(f, "({}, {})", position, pauli)?;
            if idx + 1 < self.weight() {
                write!(f, ", ")?;
            }
        }
        write!(f, "]")
    }
}

#[derive(Debug, PartialEq, Eq, Clone, Copy, Hash)]
pub enum PauliError {
    IncompatibleLength(usize, usize),
    OutOfBound(usize, usize),
}

impl fmt::Display for PauliError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            Self::IncompatibleLength(l0, l1) => {
                write!(f, "incompatible length {} and {}", l0, l1)
            }
            Self::OutOfBound(pos, len) => {
                write!(f, "position {} is out of bound for length {}", pos, len)
            }
        }
    }
}

impl Error for PauliError {}

#[cfg(test)]
mod test {
    use super::*;
    use Pauli::{X, Z};

    #[test]
    fn commutes_with_different_lengths() {
        let long_operator = PauliOperator::new(10, vec![0, 2, 7, 9], vec![X, X, X, X]);
        let short_operator = PauliOperator::new(4, vec![0, 1, 2], vec![Z, Z, Z]);
        assert!(long_operator.commutes_with(&short_operator));
    }

    #[test]
    fn anticommutes_with_different_lengths() {
        let long_operator = PauliOperator::new(10, vec![0, 2, 7, 9], vec![X, Z, X, Z]);
        let short_operator = PauliOperator::new(4, vec![0, 1, 2], vec![Z, Z, Z]);
        assert!(long_operator.anticommutes_with(&short_operator));
    }
}
