#pragma once

#include <Eigen/Sparse>
#include <Eigen/Dense>

#include "PhysConstants.h"


Eigen::SparseMatrix<double> T_sparse(int N, double dx);
Eigen::SparseMatrix<double> H_nl(int l, int N, double dx, const std::vector<double>& r_std);


inline Eigen::SparseMatrix<double> centrifugalTerm(int l, const Eigen::VectorXd& r)   // grid points (size N)
{
	/* centrifugal_term.h
   -----------------------------------------------------------------
   Builds the diagonal matrix  ħ² l(l+1) / (2 m r²)
   in sparse form (one non‑zero per row, still additive with T, V).
   Requires: Eigen ≥ 3.4 and PhysConst header.
   -----------------------------------------------------------------*/

	const int    N = r.size();
	const double pref = PhysConst::hbar2_over_two * l * (l + 1); // ħ²/2m ⋅ l(l+1)

	Eigen::SparseMatrix<double> C(N, N);
	C.reserve(Eigen::VectorXi::Constant(N, 1));	// 1 nz per row

	for (int i = 0; i < N; ++i)
		C.insert(i, i) = pref / (r[i] * r[i]);	// r must \neq 0

	return C;
}
