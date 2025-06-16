#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <numeric>
#include <cmath>


#include "matrix.h"
#include "util.h"
#include "mathematics.h"
#include "potentials.h"


int main() {


	// testing funtion r_range
	int N = 100;
	double dx = 0.01;
	std::vector<double> r_std = range(dx,10, N);


	std::vector<double> V_std = V_HO(r_std);
	
	
	Eigen::Map<const Eigen::VectorXd> r(r_std.data(), r_std.size());
	Eigen::Map<const Eigen::VectorXd> V(V_std.data(), V_std.size());
	Eigen::MatrixXd V_mat = V.asDiagonal();

	std::cout << V_mat << std::endl;


	/*Eigen::MatrixXd Lap = Laplacian(N, dx);
	std::cout << Lap << std::endl;*/


	


}




