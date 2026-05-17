// https://github.com/HorizenLabs/poseidon2/tree/main

use super::poseidon2_params::Poseidon2Params;
use ark_ff::{BigInteger, PrimeField};
use std::sync::Arc;
use crate::merkle_tree_fp::Hasher;

#[derive(Clone, Debug)]
pub struct Poseidon2<F: PrimeField> {
	pub(crate) params: Arc<Poseidon2Params<F>>,
}

impl<F: PrimeField> Poseidon2<F> {
	pub fn new(params: &Arc<Poseidon2Params<F>>) -> Self {
		Poseidon2 {
			params: Arc::clone(params),
		}
	}

	pub fn get_t(&self) -> usize {
		self.params.t
	}

	pub fn permutation(&self, input: &[F]) -> Vec<F> {
		let t = self.params.t;
		assert_eq!(input.len(), t);

		let mut current_state = input.to_owned();

		// Linear layer at beginning
		self.matmul_external(&mut current_state);

		for r in 0..self.params.rounds_f_beginning {
			current_state = self.add_rc(&current_state, &self.params.round_constants[r]);
			current_state = self.sbox(&current_state);
			self.matmul_external(&mut current_state);
		}

		let p_end = self.params.rounds_f_beginning + self.params.rounds_p;
		for r in self.params.rounds_f_beginning..p_end {
			current_state[0].add_assign(&self.params.round_constants[r][0]);
			current_state[0] = self.sbox_p(&current_state[0]);
			self.matmul_internal(&mut current_state, &self.params.mat_internal_diag_m_1);
		}

		for r in p_end..self.params.rounds {
			current_state = self.add_rc(&current_state, &self.params.round_constants[r]);
			current_state = self.sbox(&current_state);
			self.matmul_external(&mut current_state);
		}
		current_state
	}

	fn sbox(&self, input: &[F]) -> Vec<F> {
		input.iter().map(|el| self.sbox_p(el)).collect()
	}

	fn sbox_p(&self, input: &F) -> F {
		let mut input2 = *input;
		input2.square_in_place();

		match self.params.d {
			3 => {
				let mut out = input2;
				out.mul_assign(input);
				out
			}
			5 => {
				let mut out = input2;
				out.square_in_place();
				out.mul_assign(input);
				out
			}
			7 => {
				let mut out = input2;
				out.square_in_place();
				out.mul_assign(&input2);
				out.mul_assign(input);
				out
			}
			_ => {
				panic!()
			}
		}
	}

	fn matmul_m4(&self, input: &mut[F]) {
		let t = self.params.t;
		let t4 = t / 4;
		for i in 0..t4 {
			let start_index = i * 4;
			let mut t_0 = input[start_index];
			t_0.add_assign(&input[start_index + 1]);
			let mut t_1 = input[start_index + 2];
			t_1.add_assign(&input[start_index + 3]);
			let mut t_2 = input[start_index + 1];
			t_2.double_in_place();
			t_2.add_assign(&t_1);
			let mut t_3 = input[start_index + 3];
			t_3.double_in_place();
			t_3.add_assign(&t_0);
			let mut t_4 = t_1;
			t_4.double_in_place();
			t_4.double_in_place();
			t_4.add_assign(&t_3);
			let mut t_5 = t_0;
			t_5.double_in_place();
			t_5.double_in_place();
			t_5.add_assign(&t_2);
			let mut t_6 = t_3;
			t_6.add_assign(&t_5);
			let mut t_7 = t_2;
			t_7.add_assign(&t_4);
			input[start_index] = t_6;
			input[start_index + 1] = t_5;
			input[start_index + 2] = t_7;
			input[start_index + 3] = t_4;
		}
	}

	fn matmul_external(&self, input: &mut[F]) {
		let t = self.params.t;
		match t {
			2 => {
				// Matrix circ(2, 1)
				let mut sum = input[0];
				sum.add_assign(&input[1]);
				input[0].add_assign(&sum);
				input[1].add_assign(&sum);
			}
			3 => {
				// Matrix circ(2, 1, 1)
				let mut sum = input[0];
				sum.add_assign(&input[1]);
				sum.add_assign(&input[2]);
				input[0].add_assign(&sum);
				input[1].add_assign(&sum);
				input[2].add_assign(&sum);
			}
			4 => {
				// Applying cheap 4x4 MDS matrix to each 4-element part of the state
				self.matmul_m4(input);
			}
			8 | 12 | 16 | 20 | 24 => {
				// Applying cheap 4x4 MDS matrix to each 4-element part of the state
				self.matmul_m4(input);

				// Applying second cheap matrix for t > 4
				let t4 = t / 4;
				let mut stored = [F::zero(); 4];
				for l in 0..4 {
					stored[l] = input[l];
					for j in 1..t4 {
						stored[l].add_assign(&input[4 * j + l]);
					}
				}
				for i in 0..input.len() {
					input[i].add_assign(&stored[i % 4]);
				}
			}
			_ => {
				panic!()
			}
		}
	}

	fn matmul_internal(&self, input: &mut[F], mat_internal_diag_m_1: &[F]) {
		let t = self.params.t;

		match t {
			2 => {
				// [2, 1]
				// [1, 3]
				let mut sum = input[0];
				sum.add_assign(&input[1]);
				input[0].add_assign(&sum);
				input[1].double_in_place();
				input[1].add_assign(&sum);
			}
			3 => {
				// [2, 1, 1]
				// [1, 2, 1]
				// [1, 1, 3]
				let mut sum = input[0];
				sum.add_assign(&input[1]);
				sum.add_assign(&input[2]);
				input[0].add_assign(&sum);
				input[1].add_assign(&sum);
				input[2].double_in_place();
				input[2].add_assign(&sum);
			}
			4 | 8 | 12 | 16 | 20 | 24 => {
				// Compute input sum
				let mut sum = input[0];
				input
					.iter()
					.skip(1)
					.take(t-1)
					.for_each(|el| sum.add_assign(el));
				// Add sum + diag entry * element to each element
				for i in 0..input.len() {
					input[i].mul_assign(&mat_internal_diag_m_1[i]);
					input[i].add_assign(&sum);
				}
			}
			_ => {
				panic!()
			}
		}
	}

	fn add_rc(&self, input: &[F], rc: &[F]) -> Vec<F> {
		input
			.iter()
			.zip(rc.iter())
			.map(|(a, b)| {
				let mut r = *a;
				r.add_assign(b);
				r
			})
			.collect()
	}
}

impl<F: PrimeField> Hasher<F> for Poseidon2<F> {
	type Hash = [u8; 32];

	fn hash(&self, input: &[&F]) -> Self::Hash {
		let t = self.get_t();
		let mut state = vec![F::zero(); t];
		for (i, &val) in input.iter().enumerate() {
			if i >= t { break; }
			state[i] = *val;
		}

		let result_state = self.permutation(&state);

		let mut output = [0u8; 32];
		for (i, val) in result_state.iter().take(4).enumerate() {
			let repr = val.into_bigint();

			let bytes = repr.to_bytes_le();

			let start = i * 8;
			output[start..start + 8].copy_from_slice(&bytes[..8]);
		}
		output
	}
}

// #[allow(unused_imports)]
// #[cfg(test)]
// mod poseidon2_tests_goldilocks {
// 	use super::*;
// 	use crate::{fields::{goldilocks::FpGoldiLocks}};
// 	use crate::poseidon2::poseidon2_instance_goldilocks::{
// 		POSEIDON2_GOLDILOCKS_8_PARAMS,
// 		POSEIDON2_GOLDILOCKS_12_PARAMS,
// 		POSEIDON2_GOLDILOCKS_16_PARAMS,
// 		POSEIDON2_GOLDILOCKS_20_PARAMS,
// 	};
// 	use std::convert::TryFrom;
// 	use crate::utils::from_hex;
//
// 	type Scalar = FpGoldiLocks;
//
// 	static TESTRUNS: usize = 5;
//
// 	#[test]
// 	fn consistent_perm() {
// 		let instances = vec![
// 			Poseidon2::new(&POSEIDON2_GOLDILOCKS_8_PARAMS),
// 			Poseidon2::new(&POSEIDON2_GOLDILOCKS_12_PARAMS),
// 			Poseidon2::new(&POSEIDON2_GOLDILOCKS_16_PARAMS),
// 			Poseidon2::new(&POSEIDON2_GOLDILOCKS_20_PARAMS),
// 		];
// 		for instance in instances {
// 			let t = instance.params.t;
// 			for _ in 0..TESTRUNS {
// 				let input1: Vec<Scalar> = (0..t).map(|_| random_scalar()).collect();
//
// 				let mut input2: Vec<Scalar>;
// 				loop {
// 					input2 = (0..t).map(|_| random_scalar()).collect();
// 					if input1 != input2 {
// 						break;
// 					}
// 				}
//
// 				let perm1 = instance.permutation(&input1);
// 				let perm2 = instance.permutation(&input1);
// 				let perm3 = instance.permutation(&input2);
// 				assert_eq!(perm1, perm2);
// 				assert_ne!(perm1, perm3);
// 			}
// 		}
// 	}
//
// 	#[test]
// 	fn kats() {
// 		let poseidon2 = Poseidon2::new(&POSEIDON2_GOLDILOCKS_12_PARAMS);
// 		let mut input: Vec<Scalar> = vec![];
// 		for i in 0..poseidon2.params.t {
// 			input.push(Scalar::from(i as u64));
// 		}
// 		let perm = poseidon2.permutation(&input);
// 		assert_eq!(perm[0], from_hex("0x01eaef96bdf1c0c1"));
// 		assert_eq!(perm[1], from_hex("0x1f0d2cc525b2540c"));
// 		assert_eq!(perm[2], from_hex("0x6282c1dfe1e0358d"));
// 		assert_eq!(perm[3], from_hex("0xe780d721f698e1e6"));
// 		assert_eq!(perm[4], from_hex("0x280c0b6f753d833b"));
// 		assert_eq!(perm[5], from_hex("0x1b942dd5023156ab"));
// 		assert_eq!(perm[6], from_hex("0x43f0df3fcccb8398"));
// 		assert_eq!(perm[7], from_hex("0xe8e8190585489025"));
// 		assert_eq!(perm[8], from_hex("0x56bdbf72f77ada22"));
// 		assert_eq!(perm[9], from_hex("0x7911c32bf9dcd705"));
// 		assert_eq!(perm[10], from_hex("0xec467926508fbe67"));
// 		assert_eq!(perm[11], from_hex("0x6a50450ddf85a6ed"));
// 	}
// }

#[allow(unused_imports)]
#[cfg(test)]
mod poseidon2_tests_babybear {
	use super::*;
	use crate::field_babybear::{FpBabyBear};
	use crate::poseidon2_instance_babybear::{
		POSEIDON2_BABYBEAR_16_PARAMS,
		POSEIDON2_BABYBEAR_24_PARAMS,
	};
	use std::convert::TryFrom;
	use crate::utils::{from_hex, random_scalar};

	type Scalar = FpBabyBear;

	static TESTRUNS: usize = 5;

	#[test]
	fn consistent_perm() {
		let instances = vec![
			Poseidon2::new(&POSEIDON2_BABYBEAR_16_PARAMS),
			Poseidon2::new(&POSEIDON2_BABYBEAR_24_PARAMS)
		];
		for instance in instances {
			let t = instance.params.t;
			for _ in 0..TESTRUNS {
				let input1: Vec<Scalar> = (0..t).map(|_| random_scalar()).collect();

				let mut input2: Vec<Scalar>;
				loop {
					input2 = (0..t).map(|_| random_scalar()).collect();
					if input1 != input2 {
						break;
					}
				}

				let perm1 = instance.permutation(&input1);
				let perm2 = instance.permutation(&input1);
				let perm3 = instance.permutation(&input2);
				assert_eq!(perm1, perm2);
				assert_ne!(perm1, perm3);
			}
		}
	}

	#[test]
	fn kats() {
		let poseidon2 = Poseidon2::new(&POSEIDON2_BABYBEAR_24_PARAMS);
		let mut input: Vec<Scalar> = vec![];
		for i in 0..poseidon2.params.t {
			input.push(Scalar::from(i as u64));
		}
		let perm = poseidon2.permutation(&input);
		assert_eq!(perm[0], from_hex("0x2ed3e23d"));
		assert_eq!(perm[1], from_hex("0x12921fb0"));
		assert_eq!(perm[2], from_hex("0x0e659e79"));
		assert_eq!(perm[3], from_hex("0x61d81dc9"));
		assert_eq!(perm[4], from_hex("0x32bae33b"));
		assert_eq!(perm[5], from_hex("0x62486ae3"));
		assert_eq!(perm[6], from_hex("0x1e681b60"));
		assert_eq!(perm[7], from_hex("0x24b91325"));
		assert_eq!(perm[8], from_hex("0x2a2ef5b9"));
		assert_eq!(perm[9], from_hex("0x50e8593e"));
		assert_eq!(perm[10], from_hex("0x5bc818ec"));
		assert_eq!(perm[11], from_hex("0x10691997"));
		assert_eq!(perm[12], from_hex("0x35a14520"));
		assert_eq!(perm[13], from_hex("0x2ba6a3c5"));
		assert_eq!(perm[14], from_hex("0x279d47ec"));
		assert_eq!(perm[15], from_hex("0x55014e81"));
		assert_eq!(perm[16], from_hex("0x5953a67f"));
		assert_eq!(perm[17], from_hex("0x2f403111"));
		assert_eq!(perm[18], from_hex("0x6b8828ff"));
		assert_eq!(perm[19], from_hex("0x1801301f"));
		assert_eq!(perm[20], from_hex("0x2749207a"));
		assert_eq!(perm[21], from_hex("0x3dc9cf21"));
		assert_eq!(perm[22], from_hex("0x3c985ba2"));
		assert_eq!(perm[23], from_hex("0x57a99864"));
	}
}