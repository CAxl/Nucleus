#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Spectra/SymEigsSolver.h>
#include <Spectra/MatOp/SparseSymMatProd.h>
#include <numeric>
#include <cmath>

#include <matrix.h>
#include <util.h>
#include <numerics.h>
#include <shells/potentials.h>
#include <shells/hamil.h>
#include <shells/mfSolver.h>
#include <observables/BE.h>



int main() {


	// discretization of x-range
	int N = 1000;
	double dx = 0.1;

	// build r-grid
	std::vector<double> r(N);
	for (int i = 0; i < N; ++i) r[i] = (i + 1) * dx; // start at dx to avoid r=0 singularity

	// oxygen
	int A = 16;
	int Z = 8;

	// print r test
	for (int i = 0; i < r.size(); ++i) std::cout << r[i] << std::endl;
	

	// proton orbital: 0s_1/2
	{
		int l = 0;
		double j = 0.5;
		Eigen::VectorXd eigvals = solve_mean_field(A, Z, l, j, r, dx, true, 5);
		std::cout << "Proton eigenvalues (O16, 0s1/2):\n" << eigvals << "\n";
	}


	// proton 0s1/2
	{
		int l = 0;
		double j = 0.5;
		auto wfns = solve_mean_field_wavefuncs(A, Z, l, j, r, dx, true, 5, false);// last true flag is for R = u/r

		Eigen::MatrixXd u_norm = normalize_wavefunction(r, wfns[0]);

		for (int i = 0; i < r.size(); ++i) {
			std::cout << wfns[0][i] << std::endl;
		}

		/* 
		note 7/9: 
		- norm of R not correct! \int r^2 R(r)\: dr, whereas \int u(r)\: dr
		- If I actually write the u_norm() to xy-file and plot (u_norm)^2 in python I get 1 (FOR u(r) ONLY)
		
		
		*/ 
		write_xy_to_file(r, u_norm, "../../../../Output/p_norm_test.txt");
		//write_xy_to_file(r, wfns[0], "../../../../Output/p_Rwf_s1_2_1.txt");
		//write_xy_to_file(r, u_norm, "../../../../Output/p_NormRwf_s1_2_1.txt");
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

	//write_xy_to_file(r, V_tot, "../../../../Output/V_WS+V_C");
	//write_xy_to_file(r, VC, "../../../../Output/V_C");
	//write_xy_to_file(r, VWS, "../../../../Output/V_WS");

	


}








