use anyhow::{anyhow, Context, Result};
use ark_ff::{FftField};

fn ntt_in_place<F: FftField>(a: &mut [F], inverse: bool) -> Result<()> {
	let n = a.len();
	if n <= 1 {
		return Ok(());
	}
	if !n.is_power_of_two() {
		return Err(anyhow!("ntt input size must be a power of two"));
	}

	let n_log = n.trailing_zeros();
	for i in 1..n {
		let j = i.reverse_bits() >> (usize::BITS - n_log);
		if i < j {
			a.swap(i, j);
		}
	}

	let mut len = 2;
	while len <= n {
		let omega = F::get_root_of_unity(len as u64).context("failed to get root of unity")?;
		let omega = if inverse {
			omega.inverse().context("failed to invert root of unity")?
		} else {
			omega
		};

		for i in (0..n).step_by(len) {
			let mut current_omega = F::one();
			for j in 0..(len / 2) {
				let u = a[i + j];
				let v = a[i + j + len / 2] * current_omega;
				a[i + j] = u + v;
				a[i + j + len / 2] = u - v;
				current_omega *= omega;
			}
		}
		len *= 2;
	}

	if inverse {
		let n_inv = F::from(n as u64)
			.inverse()
			.context("failed to invert n")?;
		for x in a.iter_mut() {
			*x *= n_inv;
		}
	}

	Ok(())
}

pub fn ntt<F: FftField>(coeffs: &mut [F]) -> Result<()> {
	ntt_in_place(coeffs, false)
}

pub fn intt<F: FftField>(evals: &mut [F]) -> Result<()> {
	ntt_in_place(evals, true)
}