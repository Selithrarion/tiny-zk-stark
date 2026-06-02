use crate::domain::StarkDomain;
use crate::merkle_tree_fp::{MerkleProof, MerkleTree, Hasher};
use crate::poly::{eval, Poly};
use crate::fft::intt;
use anyhow::{bail, Context, Result};
use ark_ff::{FftField, PrimeField};
use ark_poly::{EvaluationDomain, GeneralEvaluationDomain};

// todo: [u8; 32]?
pub type FriHash<F> = F;

pub struct FriQueryProof<F: FftField + PrimeField> {
	pub layer_proofs: Vec<(MerkleProof<F>, (F, F), MerkleProof<F>)>,
}

pub struct FriCommitment<F: FftField + PrimeField> {
	pub roots: Vec<FriHash<F>>,
	pub final_poly: Poly<F>,
}

pub struct FriProof<F: FftField + PrimeField> {
	pub commitments: Vec<FriCommitment<F>>,
	pub queries: Vec<FriQueryProof<F>>,
}

fn compute_alpha<F, H>(hasher: &H, root: &FriHash<F>, layer_index: usize) -> F
where
	F: FftField + PrimeField,
	H: Hasher<F>,
{
	let hash_bytes = hasher.hash(&[root, &F::from(layer_index as u64)]);
	F::from_le_bytes_mod_order(hash_bytes.as_ref())
}

fn generate_query_indices<F, H>(hasher: &H, roots: &[FriHash<F>], final_poly: &Poly<F>, domain_size: usize, num_queries: usize) -> Vec<usize>
where
	F: FftField + PrimeField,
	H: Hasher<F>,
{
	let final_poly_hash = hasher.hash(&final_poly.iter().collect::<Vec<_>>());

	(0..num_queries).map(|i| {
		let i_as_field_element = F::from(i as u64);
		let final_poly_hash_as_field_element = F::from_le_bytes_mod_order(final_poly_hash.as_ref());
		let mut seed: Vec<&F> = Vec::with_capacity(roots.len() + 2);
		seed.extend(roots.iter());
		seed.push(&final_poly_hash_as_field_element);
		seed.push(&i_as_field_element);
		let hash_bytes = hasher.hash(&seed);
		let mut seed_bytes = [0u8; 8];
		seed_bytes.copy_from_slice(&hash_bytes.as_ref()[0..8]);
		(u64::from_le_bytes(seed_bytes) as usize) % domain_size
	}).collect()
}

pub fn fri_prove<F, H>(
	mut poly: Poly<F>,
	domain: &StarkDomain,
	hasher: &H,
) -> Result<FriProof<F>>
where
	F: FftField + PrimeField,
	H: Hasher<F>,
{
	// commit phase
	println!("start fri_prove");
	const FRI_STOP_LEN: usize = 1;

	let mut layers: Vec<(Poly<F>, MerkleTree<F, H>)> = Vec::new();
	let mut roots: Vec<FriHash<F>> = Vec::new();
	let mut domains: Vec<StarkDomain> = Vec::new();
	let initial_domain_size = domain.lde_domain.size();
	let mut current_domain = domain.clone();

	if !poly.len().is_power_of_two() {
		bail!("!poly.len().is_power_of_two() err");
	}

	while poly.len() > FRI_STOP_LEN {
		domains.push(current_domain.clone());

		let tree = MerkleTree::new(hasher, &poly)?;
		let root = tree.root()?;
		roots.push(root);
		layers.push((poly.clone(), tree));

		let alpha = compute_alpha(hasher, &root, layers.len());

		poly = poly
			.chunks_exact(2)
			.map(|chunk| chunk[0] + alpha * chunk[1])
			.collect();

		current_domain = StarkDomain::new(current_domain.trace_domain.size() / 2, 1)?;
	}
	let mut final_poly = poly;
	intt(&mut final_poly)?;

	// query phase
	let num_queries = 3; // todo: configurable?
	let query_indices = generate_query_indices(hasher, &roots, &final_poly, initial_domain_size, num_queries);
	let mut queries: Vec<FriQueryProof<F>> = Vec::new();

	for &initial_idx in &query_indices {
		let mut layer_proofs = Vec::new();
		let mut current_idx = initial_idx;
		for (poly, tree) in &layers {
			let sibling_idx = current_idx ^ 1;
			let val = poly[current_idx];
			let sibling_val = poly[sibling_idx];
			let val_proof = tree.prove(current_idx)?;
			let sibling_proof = tree.prove(sibling_idx)?;

			layer_proofs.push((val_proof, (val, sibling_val), sibling_proof));

			current_idx /= 2;
		}
		queries.push(FriQueryProof { layer_proofs });
	}

	let commitment = FriCommitment { roots, final_poly };

	Ok(FriProof {
		commitments: vec![commitment], // todo: handle multiple commitments
		queries,
	})
}

pub fn fri_verify<F, H>(
	proof: &FriProof<F>,
	max_degree: usize,
	initial_domain_size: usize,
	num_queries: usize,
	hasher: &H,
) -> Result<bool>
where
	F: FftField + PrimeField,
	H: Hasher<F>,
{
	// todo: handle multiple commitments
	println!("start fri_verify");
	let commitment = proof.commitments.get(0).context("no commitment found")?;
	let query_indices = generate_query_indices(hasher, &commitment.roots, &commitment.final_poly, initial_domain_size, num_queries);

	for (i, query) in proof.queries.iter().enumerate() {
		let mut current_idx = *query_indices.get(i).context("query_indices.get(i) err")?;
		let mut expected_value_for_next_layer: Option<F> = None;

		for (j, (val_proof, (val, sibling_val), sibling_proof)) in query.layer_proofs.iter().enumerate() {
			let root = commitment.roots.get(j).context("root not found")?;
			let alpha = compute_alpha(hasher, root, j + 1);
			let sibling_idx = current_idx ^ 1;

			if !MerkleTree::<F, H>::verify(hasher, root, val, current_idx, val_proof)? {
				bail!("invalid merkle proof for val at layer {j} query {i}");
			}
			if !MerkleTree::<F, H>::verify(hasher, root, sibling_val, sibling_idx, sibling_proof)? {
				bail!("invalid merkle proof for sibling_val at layer {j} query {i}");
			}

			if let Some(expected) = expected_value_for_next_layer {
				if *val != expected {
					bail!("folding consistency check failed at layer {j} query {i}");
				}
			}

			let (left_val, right_val) = if current_idx % 2 == 0 { (val, sibling_val) } else { (sibling_val, val) };
			expected_value_for_next_layer = Some(*left_val + alpha * right_val);

			current_idx /= 2;
		}

		if let Some(expected) = expected_value_for_next_layer {
			let final_domain_size = commitment.final_poly.len();
			let final_domain = GeneralEvaluationDomain::<F>::new(final_domain_size).context("failed to create final domain")?;
			let final_domain_point = final_domain.element(current_idx);
			let final_poly_eval = eval(&commitment.final_poly, final_domain_point);
			if final_poly_eval != expected {
				bail!("final polynomial evaluation mismatch at query {i}");
			}
		}
	}

	let num_foldings = commitment.roots.len();
	let max_allowed_degree = max_degree >> num_foldings;
	let final_degree = degree(&commitment.final_poly);
	if final_degree > max_allowed_degree {
		println!("fri_verify: final_degree: {}, max_allowed_degree: {}", final_degree, max_allowed_degree);
		println!("fri_verify: final_poly (coeffs): {:?}", commitment.final_poly.iter().map(|f| f.to_string()).collect::<Vec<_>>());
		bail!("final polynomial degree is too high");
	}

	Ok(true)
}

fn degree<F: PrimeField>(poly: &Poly<F>) -> usize {
	poly.iter().rposition(|&c| !c.is_zero()).unwrap_or(0)
}