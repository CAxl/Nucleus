#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <vector>

#include "hamil.h"
#include "PhysConstants.h"
#include "potentials.h"
#include "matrix.h"



Eigen::SparseMatrix<double> T_sparse(int N, double dx)
{
	const double hc = PhysConst::hbarc;
	const double m = PhysConst::mN;
	const double c = - hc * hc / (2.0 * m * dx * dx);

	Eigen::SparseMatrix<double> T(N, N);
	T.reserve(Eigen::VectorXi::Constant(N, 3));

	std::vector<Eigen::Triplet<double>> triplets;
	triplets.reserve(3 * N - 2);

	for (int i = 0; i < N; ++i)
	{
		triplets.emplace_back(i, i, -2 * c);
		if (i > 0)     triplets.emplace_back(i, i - 1, c);
		if (i < N - 1)   triplets.emplace_back(i, i + 1, c);
	}
	
	T.setFromTriplets(triplets.begin(), triplets.end());
	
	return T;
}


Eigen::SparseMatrix<double> H_nl(int l, int N, double dx, const std::vector<double>& r_std)
{
	Eigen::Map<const Eigen::VectorXd> r(r_std.data(), r_std.size());

	// potential
	std::vector<double> V_std = V_HO(r_std);
	Eigen::SparseMatrix<double> V = diagSparse(V_std);

	// kinetic
	auto T = T_sparse(N, dx);

	// cent
	Eigen::SparseMatrix<double> L2 = centrifugalTerm(l, r);

	
	return T + V + L2;
}


