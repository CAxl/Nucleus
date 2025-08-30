#include <vector>
#include <cmath>

#include "potentials.h"


std::vector<double> V_HO(const std::vector<double>& r)
{
	double m = 1;
	double E0 = -55;
	double hw = 8.6;
	double mc2 = 939.57;
	double hc = 197.326;

	//double C = 2 * mc2 / (hc * hc);

	std::vector<double> V;

	V.reserve(r.size());	// allocate memory size(r) = N

	for (auto ri : r) {
		double e = (0.5 * (mc2) * (hw * hw) * (ri * ri) / (hc * hc)) + E0;	// calculate V(r_i)

		V.push_back(e);
	}
	return V;
}

std::vector<double> V_WS(std::vector<double>& r, int A) {

	std::vector<double> V(r.size());

	double V_0 = 50;
	double r_0 = 1.5;
	double R = r_0 * std::pow(A, 1 / 3);
	double a = 0.5;



	for (size_t i = 0; i < r.size(); ++i) {
		V[i] = -V_0 / (1 + std::exp((r[i] - R) / a));

	}

	return V;
}