#![feature(test)]
extern crate test;
mod benches;
mod common;
mod gadget;
mod poseidon;
mod rescue;
mod rescue_prime;
mod sponge;
mod transcript;
#[cfg(test)]
mod tests;

pub use gadget::{
    poseidon::{poseidon_gadget, poseidon_gadget_var_length},
    rescue::{rescue_gadget, rescue_gadget_var_length},
};
pub use poseidon::{poseidon_hash, poseidon_hash_var_length};
pub use rescue::{rescue_hash, rescue_hash_var_length};
