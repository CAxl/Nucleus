#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Spectra/SymEigsSolver.h>
#include <Spectra/MatOp/SparseSymMatProd.h>
#include <numeric>
#include <cmath>


#include "matrix.h"
#include "util.h"
#include "numerics.h"
#include "shells/potentials.h"
#include "shells/hamil.h"


int main() {


	// discretization of x-range
	int N = 100;
	double dx = 0.1;
	std::vector<double> r_std = range(dx, 10, N);
	Eigen::Map<const Eigen::VectorXd> r(r_std.data(), r_std.size());


	int l = 0;

	/*---------- Hamiltonian -------------------------------*/

	for (int l = 0; l <= 2; ++l) { // s,p,d

		Eigen::SparseMatrix<double> H = H_nl(l, N, dx, r_std);

		// convert to Spectra:
		Spectra::SparseSymMatProd<double> op(H);

		int k = 5; // number of eigenvals to collect
		int m = 2 * k + 1; // size of Krylov subspace (~2k)

		Spectra::SymEigsSolver<Spectra::SparseSymMatProd<double>> eigs(op, k, m);

		eigs.init();
		int nconv = eigs.compute(Spectra::SortRule::SmallestAlge);

		Eigen::VectorXd eigvals = eigs.eigenvalues();

		std::cout << "l = " << l << " eigenvalues (MeV):\n" << eigvals.transpose() << "\n\n";
	}


}




