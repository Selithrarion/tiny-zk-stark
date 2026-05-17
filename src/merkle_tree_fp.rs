use ark_ff::PrimeField;
use std::marker::PhantomData;
use anyhow::{Result, Context};

pub trait Hasher<F: PrimeField> {
	type Hash: AsRef<[u8]>;
	fn hash(&self, input: &[&F]) -> Self::Hash;
}

#[derive(Debug)]
pub struct MerkleProof<F: PrimeField> {
	pub siblings: Vec<F>,
}

#[derive(Debug)]
pub struct MerkleTree<F: PrimeField, H: Hasher<F>> {
	layers: Vec<Vec<F>>,
	_phantom: PhantomData<(F, H)>,
}

impl<F: PrimeField, H: Hasher<F>> MerkleTree<F, H> {
	pub fn new(hasher: &H, leaves: &[F]) -> Result<Self> {
		if leaves.is_empty() {
			return Ok(MerkleTree::<F, H> {
				layers: vec![vec![F::zero()]],
				_phantom: PhantomData,
			});
		}

		let mut layers = Vec::new();
		let size = leaves.len().next_power_of_two();
		let mut current_layer = leaves.to_vec();
		current_layer.resize(size, F::zero());
		layers.push(current_layer);

		while layers.last().context("layers should not be empty")?.len() > 1 {
			let prev_layer = layers.last().context("layers should not be empty")?;
			let new_layer = prev_layer
				.chunks_exact(2)
				.map(|chunk| F::from_le_bytes_mod_order(hasher.hash(&[&chunk[0], &chunk[1]]).as_ref()))
				.collect();
			layers.push(new_layer);
		}

		Ok(MerkleTree::<F, H> {
			layers,
			_phantom: PhantomData,
		})
	}

	pub fn root(&self) -> Result<F> {
		self.layers
			.last()
			.and_then(|l| l.get(0))
			.copied()
			.context("tree is empty, no root")
	}

	pub fn prove(&self, leaf_index: usize) -> Result<MerkleProof<F>> {
		let mut siblings = Vec::new();
		let mut current_index = leaf_index;

		for layer in self.layers.iter().take(self.layers.len() - 1) {
			if current_index >= layer.len() {
				return Err(anyhow::anyhow!("leaf_index out of bounds"));
			}
			let sibling_index = current_index ^ 1;
			siblings.push(layer[sibling_index]);
			current_index /= 2;
		}

		Ok(MerkleProof { siblings })
	}

	pub fn verify(hasher: &H, root: &F, leaf: &F, leaf_index: usize, proof: &MerkleProof<F>) -> Result<bool> where F: PrimeField {
		let mut current_hash = *leaf;
		let mut current_index = leaf_index;

		for sibling in proof.siblings.iter() {
			let (left, right) = if current_index % 2 == 0 {
				(&current_hash, sibling)
			} else {
				(sibling, &current_hash)
			};
			current_hash = F::from_le_bytes_mod_order(hasher.hash(&[left, right]).as_ref());
			current_index /= 2;
		}

		Ok(current_hash == *root)
	}
}
impl<'a, F: PrimeField, H: Hasher<F>> Hasher<F> for &'a H {
	type Hash = H::Hash;
	fn hash(&self, input: &[&F]) -> Self::Hash {
		(*self).hash(input)
	}
}