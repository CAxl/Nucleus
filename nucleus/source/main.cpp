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
#include "observables/BE.h"


int main() {


	// discretization of x-range
	int N = 1000;
	double dx = 0.01;

	// build r-grid
	std::vector<double> r(N);
	for (int i = 0; i < N; ++i) r[i] = (i + 1) * dx; // start at dx to avoid r=0 singularity


	// for now (no SO)
	int l = 0;


	// oxygen
	int A = 16;
	int Z = 8;
	

	/*---------- Hamiltonian -------------------------------*/

	// choose potentials
	std::vector<PotentialTerm> pots;
	pots.push_back(WS_potential(A, Z));
	//pots.push_back(Coulomb_potential(A, Z)); // only include for protons

	// example: compute Hamiltonian for s-wave (l=0)
	Eigen::SparseMatrix<double> H = H_nl(0, r, dx, pots);


	// Spectra eigensolver
	int k = 5;				// eigenvalues wanted
	int m = 10 * k;			// Krylov subspace
	Spectra::SparseSymMatProd<double> op(H);
	Spectra::SymEigsSolver<Spectra::SparseSymMatProd<double>> eigs(op, k, m);

	eigs.init();
	int nconv = eigs.compute(Spectra::SortRule::SmallestAlge);


	Eigen::VectorXd eigvals = eigs.eigenvalues();
	std::cout << "Eigenvalues for O16 (s-wave, l=0):\n" << eigvals << "\n";



	/* BE test */


	double BE = BindingEnergy(A, Z)[0];
	double BE_per_A = BindingEnergy(A, Z)[1];
	std::cout << "Binding energy 16O = " << BE << std::endl;
	std::cout << "Binding energy per nucleon = " << BE_per_A << std::endl;


}




