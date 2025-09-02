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
#include "shells/mfSolver.h"


int main() {


	// discretization of x-range
	int N = 1000;
	double dx = 0.01;

	// build r-grid
	std::vector<double> r(N);
	for (int i = 0; i < N; ++i) r[i] = (i + 1) * dx; // start at dx to avoid r=0 singularity

	// oxygen
	int A = 16;
	int Z = 8;
	

	// proton orbital: 0s_1/2
	{
		int l = 0;
		double j = 0.5;
		Eigen::VectorXd eigvals = solve_mean_field(A, Z, l, j, r, dx, true, 5);
		std::cout << "Proton eigenvalues (O16, 0s1/2):\n" << eigvals << "\n";
	}


	// neutron 0s1/2
	{
		int l = 0;
		double j = 0.5;
		auto wfns = solve_mean_field_wavefuncs(A, Z, l, j, r, dx, false, 5);

		write_xy_to_file(r, wfns[0], "../../../../Output/n_wf_s1_2_1.txt");
	}

	// neutron 0p_3/2
	{
		int l = 1;
		double j = 1.5;
		auto wfns = solve_mean_field_wavefuncs(A, Z, l, j, r, dx, false, 5);

		write_xy_to_file(r, wfns[0], "../../../../Output/n_wf_p3_2_1.txt");
	}

	// neutron 0p_1/2
	{
		int l = 1;
		double j = 0.5;
		auto wfns = solve_mean_field_wavefuncs(A, Z, l, j, r, dx, false, 5);

		write_xy_to_file(r, wfns[0], "../../../../Output/n_wf_p1_2_1.txt");
	}
	




	/* BE test */


	double BE = BindingEnergy(A, Z)[0];
	double BE_per_A = BindingEnergy(A, Z)[1];
	std::cout << "Binding energy 16O = " << BE << std::endl;
	std::cout << "Binding energy per nucleon = " << BE_per_A << std::endl;




	/* (x,v(r)) data */
	std::vector<double> V_tot(r.size());  
	std::vector<double> VC = V_C(r, A, Z);
	std::vector<double> VWS = V_WS(r, A, Z);

	for (size_t i = 0; i < r.size(); ++i) {
		V_tot[i] = VC[i] + VWS[i];
	}

	write_xy_to_file(r, V_tot, "../../../../Output/V_WS+V_C");
	write_xy_to_file(r, VC, "../../../../Output/V_C");
	write_xy_to_file(r, VWS, "../../../../Output/V_WS");


}








