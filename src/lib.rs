/*!
 pauli is a library for manipulating quantum Pauli operators.
*/

pub mod base;
pub mod sparse;

pub use base::Pauli;
pub use Pauli::{I, X, Y, Z};
