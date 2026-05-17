use ark_ff::{FftField, Field, PrimeField};
use std::cmp::max;
use std::iter::Sum;
use anyhow::Context;
use crate::fft::{intt, ntt};

pub type Poly<F> = Vec<F>;

pub fn eval<F: Field>(poly: &[F], t: F) -> F {
	poly.iter()
		.enumerate()
		.map(|(i, &c)| c * t.pow([i as u64]))
		.sum()
}

pub fn poly_mul<F: FftField>(a: &[F], b: &[F]) -> anyhow::Result<Poly<F>> {
	if a.is_empty() || b.is_empty() {
		return Ok(vec![]);
	}

	let res_len = a.len() + b.len() - 1;
	let n = res_len.next_power_of_two();

	let mut a_padded = a.to_vec();
	a_padded.resize(n, F::zero());

	let mut b_padded = b.to_vec();
	b_padded.resize(n, F::zero());

	ntt(&mut a_padded)?;
	ntt(&mut b_padded)?;

	let mut res_padded: Vec<F> = a_padded
		.iter()
		.zip(b_padded.iter())
		.map(|(&x, &y)| x * y)
		.collect();

	intt(&mut res_padded)?;
	res_padded.truncate(res_len);
	Ok(res_padded)
}

pub fn poly_add<F: Field>(a: &[F], b: &[F]) -> Poly<F> {
	let max_len = max(a.len(), b.len());
	let mut result = vec![F::zero(); max_len];

	for i in 0..max_len {
		let val_a = a.get(i).copied().unwrap_or_else(F::zero);
		let val_b = b.get(i).copied().unwrap_or_else(F::zero);
		result[i] = val_a + val_b;
	}
	result
}

pub fn poly_sub<F: Field>(a: &[F], b: &[F]) -> Poly<F> {
	let max_len = max(a.len(), b.len());
	let mut result = vec![F::zero(); max_len];

	for i in 0..max_len {
		let val_a = a.get(i).copied().unwrap_or_else(F::zero);
		let val_b = b.get(i).copied().unwrap_or_else(F::zero);
		result[i] = val_a - val_b;
	}
	while result.last().map_or(false, |c| c.is_zero()) && result.len() > 1 {
		result.pop();
	}
	result
}


pub fn poly_div<F: PrimeField>(a: &[F], b: &[F]) -> anyhow::Result<(Poly<F>, Poly<F>)> {
	if b.iter().all(|c| c.is_zero()) {
		return Err(anyhow::anyhow!("poly_div zero div err"));
	}

	let b_deg = b.len() - 1;
	let a_deg = a.len() - 1;
	if b_deg > a_deg {
		return Ok((vec![], a.to_vec()));
	}

	if a.is_empty() {
		return Ok((vec![], vec![]));
	}

	let mut quotient = vec![F::zero(); a_deg - b_deg + 1];
	let mut remainder = a.to_vec();

	for i in (0..quotient.len()).rev() {
		let lead_b = *b.last().context("lead_b err")?;
		let lead_b_inv = lead_b.inverse().context("lead_b_inv err")?;
		let coeff = *remainder.last().context("remainder.last() err")? * lead_b_inv;
		quotient[i] = coeff;

		for (j, &b_coeff) in b.iter().enumerate() {
			remainder[i + j] -= coeff * b_coeff;
		}
		remainder.pop();
	}

	while remainder.last().map_or(false, |c| c.is_zero()) && remainder.len() > 1 {
		remainder.pop();
	}

	Ok((quotient, remainder))
}

pub fn poly_shift<F: Field>(p: &[F], shift: F) -> Poly<F> {
	let mut current_shift_power = F::one();
	p.iter().map(|&c| {
		let res = c * current_shift_power;
		current_shift_power *= shift;
		res
	}).collect()
}

pub fn poly_scale<F: Field>(p: &[F], scalar: F) -> Poly<F> {
	p.iter().map(|&x| x * scalar).collect()
}

pub fn lagrange_interpolate<F: FftField + Sum>(points: &[(F, F)]) -> anyhow::Result<Poly<F>> {
	let n = points.len();
	let mut res = vec![F::zero()];

	for i in 0..n {
		let (xi, yi) = points[i];
		let mut li = vec![F::one()];
		let mut denominator = F::one();

		for j in 0..n {
			if i == j {
				continue;
			}
			let (xj, _) = points[j];
			li = poly_mul(&li, &[-xj, F::one()])?;
			denominator = denominator * (xi - xj);
		}

		let term = poly_scale(&li, yi * denominator.inverse().unwrap()); // TODO: unwrap
		res = poly_add(&res, &term);
	}

	Ok(res)
}