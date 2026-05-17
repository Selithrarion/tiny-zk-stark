use crate::domain::StarkDomain;
use crate::field_babybear::FpBabyBear;
use crate::poly::{poly_div, poly_mul, poly_shift, poly_sub, Poly};
use anyhow::{anyhow, bail, Context, Result};
use ark_ff::{Field, One, Zero};
use ark_poly::EvaluationDomain;

pub trait Air {
	fn trace_len(&self) -> usize;
	fn public_inputs(&self) -> Vec<FpBabyBear>;
	fn get_execution_trace(&self) -> Result<Vec<Vec<FpBabyBear>>>;
	fn get_constraint_polynomials(
		&self,
		trace_polys: &[Poly<FpBabyBear>],
		domain: &StarkDomain,
	) -> Result<Vec<Poly<FpBabyBear>>>;
	fn evaluate_constraints_on_lde(
		&self,
		trace_lde: &[FpBabyBear],
		domain: &StarkDomain,
	) -> Result<Vec<Vec<FpBabyBear>>>;
	fn evaluate_constraints(
		&self,
		z: FpBabyBear,
		trace_evals: &[FpBabyBear],
		domain: &StarkDomain,
	) -> Result<Vec<FpBabyBear>>;
}

pub struct FibonacciAir {
	trace_len: usize,
	result: FpBabyBear,
}

impl FibonacciAir {
	pub fn new(trace_len: usize, result: FpBabyBear) -> Self {
		Self { trace_len, result }
	}
}

impl Air for FibonacciAir {
	fn trace_len(&self) -> usize {
		self.trace_len
	}

	fn public_inputs(&self) -> Vec<FpBabyBear> {
		vec![FpBabyBear::one(), self.result]
	}

	fn get_execution_trace(&self) -> Result<Vec<Vec<FpBabyBear>>> {
		let mut trace = vec![FpBabyBear::from(0u64); self.trace_len];
		trace[0] = FpBabyBear::from(1u64);
		trace[1] = FpBabyBear::from(1u64);

		for i in 2..self.trace_len {
			trace[i] = trace[i - 1] + trace[i - 2];
		}

		let last_element = *trace.last().context("trace is empty")?;
		if last_element != self.result {
			return Err(anyhow::anyhow!(
				"fibonacci result mismatch: expected {}, got {}",
				self.result,
				last_element
			));
		}

		Ok(vec![trace])
	}

	fn get_constraint_polynomials(
		&self,
		trace_polys: &[Poly<FpBabyBear>],
		domain: &StarkDomain,
	) -> Result<Vec<Poly<FpBabyBear>>> {
		if trace_polys.len() != 1 {
			return Err(anyhow!("expected 1 trace polynomial for fibonacci"));
		}
		let p = &trace_polys[0];

		let g = domain.trace_domain.group_gen();
		let g_inv = g.inverse().context("g must be invertible")?;

		// P(x) - 1 / (x - g^0)
		let p_minus_1 = poly_sub(p, &[FpBabyBear::one()]);
		let (boundary_constraint_0, rem) =
			poly_div(&p_minus_1, &[-FpBabyBear::one(), FpBabyBear::one()])?;
		if !rem.iter().all(|c| c.is_zero()) {
			return Err(anyhow!("boundary constraint a_0=1 does not hold"));
		}

		// P(x) / (x - g^{n-1})
		let last_point = g.pow([(self.trace_len - 1) as u64]);
		let p_minus_res = poly_sub(p, &[self.result]);
		let (boundary_constraint_1, rem) = poly_div(&p_minus_res, &[-last_point, FpBabyBear::one()])?;
		if !rem.iter().all(|c| c.is_zero()) {
			return Err(anyhow!("boundary constraint a_n-1=result does not hold"));
		}

		// P(x*g^2) - P(x*g) - P(x) = 0
		let p_next = poly_shift(p, g);
		let p_next_next = poly_shift(&p_next, g);
		let transition_numerator = poly_sub(&poly_sub(&p_next_next, &p_next), p);

		// the zerofier
		// Z(x) = (x^n - 1) / ( (x - g^{n-1}) * (x - g^{n-2}) )
		let mut z_num = vec![FpBabyBear::zero(); self.trace_len + 1];
		z_num[self.trace_len] = FpBabyBear::one();
		z_num[0] = -FpBabyBear::one();
		let z_den_0 = &[-last_point, FpBabyBear::one()];
		let z_den_1 = &[-last_point * g_inv, FpBabyBear::one()];
		let z_den = poly_mul(z_den_0, z_den_1)?;
		let (z, rem) = poly_div(&z_num, &z_den)?;
		if !rem.iter().all(|c| c.is_zero()) {
			return Err(anyhow!("zerofier polynomial division has a remainder"));
		}

		let (transition_constraint, rem) = poly_div(&transition_numerator, &z)?;
		if !rem.iter().all(|c| c.is_zero()) {
			return Err(anyhow!("transition constraint division has a remainder"));
		}

		Ok(vec![
			boundary_constraint_0,
			boundary_constraint_1,
			transition_constraint,
		])
	}

	fn evaluate_constraints_on_lde(
		&self,
		trace_lde: &[FpBabyBear],
		domain: &StarkDomain,
	) -> Result<Vec<Vec<FpBabyBear>>> {
		// todo: precompute evaluations of inverse
		let lde_size = domain.lde_domain.size();
		let trace_gen = domain.trace_generator();
		let lde_domain_elements: Vec<_> = domain.lde_domain.elements().collect();

		let mut constraint_evals = Vec::new();

		// boundary constraint for a_0 = 1
		// (P(z) - 1) / (z - 1)
		let boundary_constraint_0: Vec<_> = trace_lde
			.iter()
			.zip(lde_domain_elements.iter())
			.map(|(&p_val, &x)| (p_val - FpBabyBear::one()) * (x - FpBabyBear::one()).inverse().unwrap())
			.collect();
		constraint_evals.push(boundary_constraint_0);

		// boundary constraint for a_{n-1} = result
		// (P(z) - result) / (z - g^{n-1})
		let last_point = trace_gen.pow([(self.trace_len - 1) as u64]);
		let boundary_constraint_1: Vec<_> = trace_lde
			.iter()
			.zip(lde_domain_elements.iter())
			.map(|(&p_val, &x)| (p_val - self.result) * (x - last_point).inverse().unwrap())
			.collect();
		constraint_evals.push(boundary_constraint_1);

		// transition constraint
		let mut transition_numerators = vec![FpBabyBear::zero(); lde_size];
		for i in 0..lde_size {
			let p_i = trace_lde[i];
			let p_i_next = trace_lde[(i + domain.blowup_factor) % lde_size];
			let p_i_next_next = trace_lde[(i + 2 * domain.blowup_factor) % lde_size];
			transition_numerators[i] = p_i_next_next - p_i_next - p_i;
		}

		let n = self.trace_len as u64;
		let z_den_p1 = trace_gen.pow([n - 1]);
		let z_den_p2 = trace_gen.pow([n - 2]);
		let z_den_evals: Vec<_> = lde_domain_elements.iter().map(|x| (*x - z_den_p1) * (*x - z_den_p2)).collect();
		let z_num_evals: Vec<_> = lde_domain_elements.iter().map(|x| x.pow([n]) - FpBabyBear::one()).collect();
		let z_inv_evals: Vec<_> = z_num_evals.iter().zip(z_den_evals.iter()).map(|(num, den)| *den * num.inverse().unwrap()).collect();

		let transition_constraint: Vec<_> = transition_numerators.iter().zip(z_inv_evals.iter()).map(|(num, z_inv)| *num * *z_inv).collect();
		constraint_evals.push(transition_constraint);

		Ok(constraint_evals)
	}

	fn evaluate_constraints(
		&self,
		z: FpBabyBear,
		trace_evals: &[FpBabyBear],
		domain: &StarkDomain,
	) -> Result<Vec<FpBabyBear>> {
		if trace_evals.len() < 1 {
			bail!("at least one trace evaluation is required");
		}
		let p_z = trace_evals[0];

		// boundary constraint for a_0 = 1
		// (P(z) - 1) / (z - 1)
		let boundary_0_num = p_z - FpBabyBear::one();
		let boundary_0_den = z - FpBabyBear::one();
		let boundary_constraint_0 =
			boundary_0_num * boundary_0_den.inverse().context("z cannot be 1")?;

		// boundary constraint for a_{n-1} = result
		// (P(z) - result) / (z - g^{n-1})
		let g = domain.trace_domain.group_gen();
		let last_point = g.pow([(self.trace_len - 1) as u64]);
		let boundary_1_num = p_z - self.result;
		let boundary_1_den = z - last_point;
		let boundary_constraint_1 =
			boundary_1_num * boundary_1_den.inverse().context("z cannot be the last point")?;

		// transition constraint
		if trace_evals.len() < 3 {
			bail!("at least 3 trace evaluations are required for transition constraints");
		}
		let p_z_g = trace_evals[1];
		let p_z_g2 = trace_evals[2];
		let transition_constraint_val = p_z_g2 - p_z_g - p_z;

		Ok(vec![
			boundary_constraint_0,
			boundary_constraint_1,
			transition_constraint_val,
		])
	}
}