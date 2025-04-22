#include <vector>
#include "matrix.h"



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
