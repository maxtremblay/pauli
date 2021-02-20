//! A library for quantum Pauli operators.
//!
//! This contains a definition of the single-qubit
//! Pauli operators `I`, `X`, `Y` and `Z` and
//! of general multi-qubit Pauli operators such as `XYZ`.
//!
//! This library is built with error correction in mind
//! thus the phases are ignored.

mod base;
pub use base::Pauli;
pub use Pauli::{I, X, Y, Z};

mod operator;
pub use operator::{PauliError, PauliOperator};
