use crate::Pauli;
use crate::Phase;
use itertools::Itertools;

/// A dense Pauli operator is a global phase
/// and a list of single qubit Paulis.
///
/// # Example
///
/// ```
/// # use pauli::DensePauliOperator;
/// use pauli::{I, X, Y, Z};
/// use pauli::Phase;
///
/// // This instanciate a new operator with global phase -i
/// // and X, Y and Z acting on the first, second and third
/// // qubits respectively.
/// let operator = DensePauliOperator::with_phase_and_paulis(Phase::i(), vec![X, Y, Z]);
///
/// // You don't need to specified the phase. This will create an
/// // operator with phase 1 and X acting on each qubits.
/// let other_operator = DensePauliOperator::with_paulis(vec![X, X, X]);
///
/// // But it is easy to create the same operator with a different phase afterwhile.
/// let other_operator_with_phase = Phase::minus_i() * &other_operator;
/// assert_eq!(
///     other_operator_with_phase,
///     DensePauliOperator::with_phase_and_paulis(Phase::minus_i(), vec![X, X, X])
/// );
///
/// // Or to apply the same Pauli to all qubits.
/// let other_operator_with_pauli = Y * &other_operator;
/// assert_eq!(other_operator_with_pauli, DensePauliOperator::with_phase_and_paulis(Phase::minus_i(), vec![Z, Z, Z]));
///
/// // These two operators commute.
/// assert!(operator.commutes_with(&other_operator));
///
/// // Finally, it is easy to compute their product.
/// let product = &operator * &other_operator;
/// assert_eq!(
///     product,
///     DensePauliOperator::with_phase_and_paulis(Phase::one(), vec![I, Z, Y])
/// );
/// ```
#[derive(Debug, PartialEq, Eq, Clone, Hash)]
pub struct DensePauliOperator {
    paulis: Vec<Pauli>,
    phase: Phase,
}

impl DensePauliOperator {
    /// Creates a new empty operator.
    pub fn new() -> Self {
        Self {
            paulis: Vec::new(),
            phase: Phase::one(),
        }
    }

    /// Creates a new operator from the given Paulis with
    /// a phase of 1.
    pub fn with_paulis(paulis: Vec<Pauli>) -> Self {
        Self {
            paulis,
            phase: Phase::one(),
        }
    }

    /// Creates a new operator from the given phase and Paulis.
    pub fn with_phase_and_paulis(phase: Phase, paulis: Vec<Pauli>) -> Self {
        Self { paulis, phase }
    }

    /// Returns the length of the operator.
    pub fn len(&self) -> usize {
        self.paulis.len()
    }

    /// Returns true if the operator contains no Paulis.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Returns an iterator over all single qubit Pauli
    /// in the operators.
    pub fn paulis(&self) -> impl Iterator<Item = Pauli> + '_ {
        self.paulis.iter().cloned()
    }

    /// Returns the phase of the operator.
    pub fn phase(&self) -> Phase {
        self.phase
    }

    /// Returns an iterator over the positions where single qubit
    /// Paulis are not identity.
    ///
    /// # Example
    ///
    /// ```
    /// # use pauli::DensePauliOperator;
    /// # use pauli::{I, X, Y, Z};
    /// let operator = DensePauliOperator::with_paulis(vec![X, I, Y, I, Z, I]);
    /// let mut iter = operator.non_trivial_positions();
    ///
    /// assert_eq!(iter.next(), Some(0));
    /// assert_eq!(iter.next(), Some(2));
    /// assert_eq!(iter.next(), Some(4));
    /// assert_eq!(iter.next(), None);
    /// ```
    pub fn non_trivial_positions(&self) -> impl Iterator<Item = usize> + '_ {
        self.paulis().positions(Pauli::is_non_trivial)
    }

    /// Returns an iterator yielding (position, Pauli) pairs where
    /// single qubit Paulis are not identity.
    ///
    /// # Example
    ///
    /// ```
    /// # use pauli::DensePauliOperator;
    /// # use pauli::{I, X, Y, Z};
    /// let operator = DensePauliOperator::with_paulis(vec![X, I, Y, I, Z, I]);
    /// let mut iter = operator.non_trivial_paulis();
    ///
    /// assert_eq!(iter.next(), Some((0, X)));
    /// assert_eq!(iter.next(), Some((2, Y)));
    /// assert_eq!(iter.next(), Some((4, Z)));
    /// assert_eq!(iter.next(), None);
    /// ```
    pub fn non_trivial_paulis(&self) -> impl Iterator<Item = (usize, Pauli)> + '_ {
        self.paulis().enumerate().filter_map(|(position, pauli)| {
            if pauli.is_non_trivial() {
                Some((position, pauli.clone()))
            } else {
                None
            }
        })
    }

    /// Checks if two operators commute.
    ///
    /// if operators have different lengths, the shortest one is padded with I.
    ///
    /// # Example
    ///
    /// ```
    /// # use pauli::DensePauliOperator;
    /// # use pauli::{I, X, Y, Z};
    /// let first_operator = DensePauliOperator::with_paulis(vec![X, Y, Z]);
    /// let second_operator = DensePauliOperator::with_paulis(vec![Y, Y, Y]);
    /// let third_operator = DensePauliOperator::with_paulis(vec![I, X, I]);
    ///
    /// assert_eq!(first_operator.commutes_with(&second_operator), true);
    /// assert_eq!(first_operator.commutes_with(&third_operator), false);
    /// assert_eq!(second_operator.commutes_with(&third_operator), false);
    /// ```
    ///
    /// # Panic
    ///
    /// Panics if operator have different lengths.
    pub fn commutes_with(&self, other: &Self) -> bool {
        assert_same_length(self, other);
        self.paulis()
            .zip(other.paulis())
            .filter(|(p, q)| p.anticommutes_with(*q))
            .count()
            % 2
            == 0
    }

    /// Checks if two operators anticommute.
    ///
    /// # Example
    ///
    /// ```
    /// # use pauli::DensePauliOperator;
    /// # use pauli::{I, X, Y, Z};
    /// let first_operator = DensePauliOperator::with_paulis(vec![X, Y, Z]);
    /// let second_operator = DensePauliOperator::with_paulis(vec![Y, Y, Y]);
    /// let third_operator = DensePauliOperator::with_paulis(vec![I, X, I]);
    ///
    /// assert_eq!(first_operator.anticommutes_with(&second_operator), false);
    /// assert_eq!(first_operator.anticommutes_with(&third_operator), true);
    /// assert_eq!(second_operator.anticommutes_with(&third_operator), true);
    /// ```
    ///
    /// # Panic
    ///
    /// Panics if operator have different lengths.
    pub fn anticommutes_with(&self, other: &Self) -> bool {
        assert_same_length(self, other);
        self.paulis()
            .zip(other.paulis())
            .filter(|(p, q)| p.anticommutes_with(*q))
            .count()
            % 2
            == 1
    }
}

impl std::ops::Mul<&DensePauliOperator> for &DensePauliOperator {
    type Output = DensePauliOperator;

    fn mul(self, other: &DensePauliOperator) -> DensePauliOperator {
        assert_same_length(self, other);
        let (phase, paulis) = self.paulis().zip(other.paulis()).fold(
            (Phase::one(), Vec::with_capacity(self.len())),
            |(total_phase, mut paulis), (p, q)| {
                let (phase, pauli) = p.multiply_with_phase(q);
                paulis.push(pauli);
                (total_phase * phase, paulis)
            },
        );
        DensePauliOperator {
            phase: phase * self.phase * other.phase,
            paulis,
        }
    }
}

impl std::ops::Mul<Phase> for &DensePauliOperator {
    type Output = DensePauliOperator;

    fn mul(self, phase: Phase) -> DensePauliOperator {
        DensePauliOperator {
            paulis: self.paulis.clone(),
            phase: self.phase * phase,
        }
    }
}

impl std::ops::Mul<Pauli> for &DensePauliOperator {
    type Output = DensePauliOperator;

    fn mul(self, pauli: Pauli) -> DensePauliOperator {
        let (phase, paulis) = self.paulis().fold(
            (Phase::one(), Vec::with_capacity(self.len())),
            |(total_phase, mut paulis), p| {
                let (phase, product) = p.multiply_with_phase(pauli);
                paulis.push(product);
                (total_phase * phase, paulis)
            },
        );
        DensePauliOperator {
            paulis,
            phase: self.phase * phase,
        }
    }
}

macro_rules! impl_commutative_mul {
    ($t:ty) => {
        impl std::ops::Mul<&DensePauliOperator> for $t {
            type Output = DensePauliOperator;

            fn mul(self, operator: &DensePauliOperator) -> DensePauliOperator {
                operator * self
            }
        }
    };
}

impl_commutative_mul!(Phase);
impl_commutative_mul!(Pauli);

fn assert_same_length(first: &DensePauliOperator, second: &DensePauliOperator) {
    if first.len() != second.len() {
        panic!(
            "operators have different length: {} and {}",
            first.len(),
            second.len()
        );
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use Pauli::{I, X, Y, Z};

    #[test]
    fn operator_operator_multiplication() {
        let first = DensePauliOperator::with_phase_and_paulis(Phase::i(), vec![I, X, Y, Z]);
        let second = DensePauliOperator::with_phase_and_paulis(Phase::one(), vec![X, Z, X, Z]);
        let product = DensePauliOperator::with_phase_and_paulis(Phase::minus_i(), vec![X, Y, Z, I]);
        assert_eq!(&first * &second, product);
        assert_eq!(&second * &first, product);
    }

    #[test]
    fn operator_phase_multiplication() {
        let operator =
            DensePauliOperator::with_phase_and_paulis(Phase::minus_one(), vec![X, Z, I, Y]);
        let phase = Phase::i();
        let product = DensePauliOperator::with_phase_and_paulis(Phase::minus_i(), vec![X, Z, I, Y]);
        assert_eq!(&operator * phase, product);
        assert_eq!(phase * &operator, product);
    }

    #[test]
    fn operator_pauli_multiplication() {
        let operator =
            DensePauliOperator::with_phase_and_paulis(Phase::minus_one(), vec![Y, X, Z, I]);
        let pauli = Z;
        let product =
            DensePauliOperator::with_phase_and_paulis(Phase::minus_one(), vec![X, Y, I, Z]);
        assert_eq!(&operator * pauli, product);
        assert_eq!(pauli * &operator, product);
    }
}
