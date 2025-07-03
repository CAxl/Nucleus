#include <Eigen/Sparse>
#include <Eigen/Dense>

#include "hamil.h"
#include "PhysConstants.h"




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


