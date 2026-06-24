#include <vector>
#include <Eigen/Dense>
#include <Eigen/Sparse>


#include <matrix.h>


std::vector<double> range(double start, double stop, int points)
{
	/*
	Function r_range:
	:param: float start: starting value of vector
	:param: flaot stop: stopping point of vector (included)
	:param: int points: gives number of elements in the vector -> dr

	functions like np.linspace
	*/


	std::vector<double> r_range;
	r_range.reserve(points);

	double dr = (stop - start) / (points - 1);

	for (int i = 0; i < points; i++)
		r_range.push_back(start + i * dr);

	return r_range;
}

Eigen::MatrixXd IdMat(int N) {

	Eigen::MatrixXd Id = Eigen::MatrixXd::Zero(N, N);
	Id.diagonal(0).setConstant(1);

	return Id;
}


Eigen::MatrixXd Laplacian(int N, double dx) {

	Eigen::MatrixXd Lapl = Eigen::MatrixXd::Zero(N, N);
	
	Lapl.diagonal(0).setConstant(-2);
	Lapl.diagonal(1).setConstant(1);
	Lapl.diagonal(-1).setConstant(1);

	Lapl = std::pow(1/dx,2) * Lapl;

	return Lapl;
}

