use crate::field_babybear::FpBabyBear;
use anyhow::{anyhow, Result};
use ark_ff::FftField;
use ark_poly::{EvaluationDomain, GeneralEvaluationDomain};

#[derive(Debug, Clone)]
pub struct StarkDomain {
	pub trace_domain: GeneralEvaluationDomain<FpBabyBear>,
	pub lde_domain: GeneralEvaluationDomain<FpBabyBear>,
	pub blowup_factor: usize,
}

impl StarkDomain {
	pub fn new(trace_len: usize, blowup_factor: usize) -> Result<Self> {
		// if !trace_len.is_power_of_two() {
		// 	return Err(anyhow!("trace_len must be a power of two"));
		// }
		if !blowup_factor.is_power_of_two() {
			return Err(anyhow!("blowup_factor must be a power of two"));
		}

		let trace_domain = GeneralEvaluationDomain::new(trace_len)
			.ok_or_else(|| anyhow!("failed to create trace domain of size {}", trace_len))?;

		let lde_len = trace_len * blowup_factor;
		let lde_domain = GeneralEvaluationDomain::new_coset(lde_len, FpBabyBear::GENERATOR)
			.ok_or_else(|| anyhow!("failed to create lde coset domain of size {}", lde_len))?;

		Ok(Self {
			trace_domain,
			lde_domain,
			blowup_factor,
		})
	}

	pub fn trace_generator(&self) -> FpBabyBear {
		self.trace_domain.group_gen()
	}

	pub fn lde_generator(&self) -> FpBabyBear {
		self.lde_domain.group_gen()
	}
}