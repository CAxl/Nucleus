#include <vector>
#include <cmath>
#include "potentials.h"
#include "matrix.h" // for diagSparse()
#include "PhysConstants.h"



// raw radial potential functions:

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

std::vector<double> V_WS(const std::vector<double>& r, int A, int Z) {	// A-dependence of V_0, we use Suhonen's example parametrization on pp. 44

	std::vector<double> V(r.size());

	int N = A - Z;
	double V_0 = 51 - 33 * (N - Z) / A;	// for now: we can only deal with neutrons ('+' for protons, suhonen)

	double r_0 = 1.27;
	double R = r_0 * std::pow(A, 1.0 / 3.0);
	double a = 0.67;



	for (size_t i = 0; i < r.size(); ++i) {
		V[i] = -V_0 / (1 + std::exp((r[i] - R) / a));

	}

	return V;
}

std::vector<double> V_C(const std::vector<double>& r, int A, int Z)
{
	std::vector<double> V;
	V.reserve(r.size());

	double R = 1.27 * pow(A, 1.0 / 3.0); // nuclear radius (fm, same as in WS)

	for (double ri : r)
	{
		double val;
		if (ri <= R)
		{
			val = PhysConst::e2 * Z * (3.0 - pow(ri / R, 2.0)) / (2.0 * R);
		}
		else
		{
			val = PhysConst::e2 * Z / ri;
		}

		V.push_back(val);
	}

	return V;
}

//std::vector<double> V_SO(const std::vector<double>& r, int l, double j, int A, int Z)
//{
//	
//}



// wrappers into PotentialTerm:

PotentialTerm HO_potential()
{
	return [](const std::vector<double>& r) {
		return diagSparse(V_HO(r));
		};
}

PotentialTerm WS_potential(int A, int Z)
{
	return [=](const std::vector<double>& r) {
		return diagSparse(V_WS(r, A, Z));
		};
}

PotentialTerm Coulomb_potential(int A, int Z)
{
	return [=](const std::vector<double>& r) {
		return diagSparse(V_C(r, A, Z));
		};
}

//PotentialTerm SO_potential(int l, double j, int A, int Z)
//{
//	return [=](const std::vector<double>& r) {
//		return diagSparse(V_SO(r, l, j, A, Z));
//		};
//}



