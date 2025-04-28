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


	//std::vector<double> mat = {
	//	4, 1, 1,
	//	1, 3, 0,
	//	1, 0, 2

	//};
	//std::vector<double> vec = { 1.0,2.0,3.0 };

	//vec.resize(5);	// vec = [1.0, 2.0, 3.0, 0, 0]
	//vec[3] = 4.5;

	//vec.push_back(6.0);	// appends 6.0 as new last item in vec

	//vec.pop_back();

	//for (int i = 0; i < vec.size(); i++)
	//	std::cout << vec[i] << std::endl;

	//Eigen::MatrixXd A(3, 3);

	//A << 1, 2, 3, 
	//	4, 5, 6,
	//	7, 8, 9;

	//Eigen::VectorXd x(3);

	//x << 1, 0, -1;

	//Eigen::VectorXd y = 2 * x;
	//Eigen::MatrixXd At = A.transpose();


	//// create 100 x 100 Id matrix
	//int N1 = 100;
	//Eigen::VectorXd ones(N1);
	//
	//for (int i = 0; i < N1; i++)
	//	ones[i] = 1;
	//
	//Eigen::MatrixXd I(N1, N1);
	//I = ones.asDiagonal();
	////std::cout << I << std::endl;
	//
	//Eigen::MatrixXd y = 2 * I;



	// testing funtion of r_range
	//int N2 = 1000;
	//std::vector<double> r = range(0.01,10,N2);

	//print_matrix(r);

	//std::vector<double> V = V_HO(r);
	//print_matrix(V);


	// testing off diags for (eventually) central difference in matrix form
	//int dim = 5;
	//Eigen::MatrixXd m = Eigen::MatrixXd::Zero(dim, dim);

	//m.diagonal(1).setConstant(1);
	//m.diagonal(-1).setConstant(1);
	//m.diagonal(0).setConstant(-2);

	//std::cout << m << std::endl;


	int N = 100;
	double start = 0.0;
	double stop = 10;
	double dx = (stop - start) / (N - 1);
	
	std::vector<double> x = range(start, stop, N);
	std::vector<double> f;
	f.reserve(N);


	for (int i = 0; i < N; ++i) {	// f(x) = x^2
		f.push_back(x[i] * x[i]);
	}
	std::vector<double> d2f = ddiff(f, dx);
	
	write_xy_to_file(x, d2f, "./plot_data/x_d2f.txt");

	std::vector<double> V = V_HO(x);
	write_xy_to_file(x, V, "./plot_data/HO.txt");

	std::vector<double> W = V_WS(x, 10);
	write_xy_to_file(x, W, "./plot_data/WS.txt");

	std::vector<double> Wprime = ddiff(W, dx);
	write_xy_to_file(x, Wprime, "./plot_data/WSderivative.txt");

}




