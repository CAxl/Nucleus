#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <numeric>


#include "matrix.h"
#include "util.h"
#include "mathematics.h"




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

	for (int i = 0; i < N; ++i) {
		std::cout << d2f[i] << std::endl;
	}








}



std::vector<double> V_HO(std::vector<double> r)
{
	double m = 1;
	double E0 = -55;
	double hw = 8.6;
	double mc2 = 939.57;
	double hc = 197.326;

	double C = 2 * mc2 / (hc * hc);

	std::vector<double> V;

	V.reserve(r.size());	// allocate memory size(r) = N

	for (auto ri : r) {
		double e = (0.5 * (mc2) * (hw * hw) * (ri * ri) / (hc * hc)) + E0;	// calculate V(r_i)

		V.push_back(C * e);
	}
	return V;
}
