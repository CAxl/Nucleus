#include <iostream>
#include <vector>
#include <Eigen/Dense>

#include "numerics.h"




std::vector<double> ddiff(std::vector<double>& f, double dx) {
	/*
	param: std::vector<double>& f -> function (vector) to be differentiated
	param: dx -> mesh size

	returns 2nd derivative of function (vector) f.
	Central finite differences are used in the interior 
	Forward and backward finte differences at the boundaries (i = 0, i = N-1)

	*/

	std::vector<double> d2f;

	int N = f.size();
	d2f.reserve(N);

	double forward = (f[2] - 2 * f[1] + f[0]) / (dx * dx);
	d2f.push_back(forward);
	
	for (int i = 1; i < N-1; ++i) {
		double cdiff = (f[i + 1] - 2 * f[i] + f[i - 1]) / (dx * dx);
		d2f.push_back(cdiff);
	}

	double backward = (f[N - 1] - 2 * f[N - 2] + f[N - 3]) / (dx * dx);
	d2f.push_back(backward);


	return d2f;
}


double simpson(const std::vector<double>& x, const Eigen::VectorXd& y)
{
	int N = x.size();
	if (N < 2 || N != y.size()) {
		throw std::runtime_error("Invalid input size for Simpson integration");
	}

	double h = x[1] - x[0]; // assume uniform spacing
	double sum = y[0] + y[N - 1];

	for (int i = 1; i < N - 1; i++)
	{
		if (i % 2 == 0) {
			sum += 2.0 * y[i];
		}
		else {
			sum += 4.0 * y[i];
		}
	}

	return sum * h / 3.0;
}




// Euler integrator

