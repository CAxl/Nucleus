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


	/*----------- grid ------------------------------------ */
	int N = 100;
	double dx = 0.01;
	std::vector<double> r_std = range(dx, 10, N);
	Eigen::Map<const Eigen::VectorXd> r(r_std.data(), r_std.size());


	/* ---- operators --------------------------------------*/

	std::vector<double> V_std = V_HO(r_std);
	Eigen::SparseMatrix<double> V = diagSparse(V_std);

	int l = 2;
	Eigen::SparseMatrix<double> L2 = centrifugalTerm(l, r);

	auto T = T_sparse(N, dx);

	/*---------- Hamiltonian -------------------------------*/

	Eigen::SparseMatrix<double> H = T + V + L2;

	


	


}




