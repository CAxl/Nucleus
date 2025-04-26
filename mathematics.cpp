#include <iostream>
#include <vector>

#include "mathematics.h"




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
