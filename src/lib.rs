/*!
 pauli is a library for manipulating quantum Pauli operators.
*/

#[macro_use]
extern crate impl_ops;

pub mod base;

pub use base::Pauli;
pub use Pauli::{I, X, Y, Z};
