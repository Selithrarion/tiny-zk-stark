use anyhow::Result;
use ark_ff::{One, PrimeField, Zero};
use ark_poly::EvaluationDomain;
use crate::air::{Air, FibonacciAir};
use crate::domain::StarkDomain;
use crate::fft::{intt, ntt};
use crate::field_babybear::FpBabyBear;
use crate::merkle_tree_fp::{Hasher, MerkleTree};
use crate::poseidon2::Poseidon2;
use crate::poseidon2_instance_babybear::POSEIDON2_BABYBEAR_16_PARAMS;

pub mod air;
pub mod fri;
pub mod poly;
pub mod fft;
pub mod poseidon2;
pub mod poseidon2_instance_babybear;
pub mod merkle_tree_fp;
pub mod field_babybear;
pub mod poseidon2_params;
pub mod utils;
pub mod domain;

fn main() -> Result<()> {
    println!("start main");
    const TRACE_LEN: usize = 16;
    const BLOWUP_FACTOR: usize = 4;
    let fib_result = FpBabyBear::from(987u64); // fib(16)

    let air = FibonacciAir::new(TRACE_LEN, fib_result);
    let domain = StarkDomain::new(TRACE_LEN, BLOWUP_FACTOR)?;
    let hasher = Poseidon2::new(&POSEIDON2_BABYBEAR_16_PARAMS);

    println!("main: trace generation");
    let trace_columns = air.get_execution_trace()?;
    let trace_evals = &trace_columns[0];

    let mut trace_coeffs = trace_evals.clone();
    intt(&mut trace_coeffs)?;

    let mut trace_lde = trace_coeffs.clone();
    println!("main: trace_lde len before ntt: {}", trace_lde.len());
    trace_lde.resize(domain.lde_domain.size(), FpBabyBear::zero());
    ntt(&mut trace_lde)?;
    let trace_tree = MerkleTree::new(&hasher, &trace_lde)?;
    let trace_root = trace_tree.root()?;
    println!("trace root: {}", trace_root);

    let constraint_evals_on_lde = air.evaluate_constraints_on_lde(&trace_lde, &domain)?;
    println!("main: got {} constraint evaluations", constraint_evals_on_lde.len());

    let hash_bytes = hasher.hash(&[&trace_root]);
    let base_alpha = FpBabyBear::from_le_bytes_mod_order(&hash_bytes);

    let mut composition_lde = vec![FpBabyBear::zero(); domain.lde_domain.size()];
    let mut current_alpha = FpBabyBear::one();
    for constraint_evals in constraint_evals_on_lde.iter() {
        current_alpha *= base_alpha;
        for (i, eval) in constraint_evals.iter().enumerate() {
            composition_lde[i] += current_alpha * *eval;
        }
    }

    let max_degree_composition = TRACE_LEN;
    println!("main: max_degree_composition: {}", max_degree_composition);

    let fri_proof = fri::fri_prove(composition_lde, &domain, &hasher)?;
    let is_valid = fri::fri_verify(
        &fri_proof,
        max_degree_composition,
        domain.lde_domain.size(),
        3,
        &hasher,
    )?;

    println!("proof is valid: {}", is_valid);

    Ok(())
}