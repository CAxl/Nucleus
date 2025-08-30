#include <cmath>
#include <array>





int delta(int a, int b) { // Kroenecker delta
	return a == b;
}



std::array<double,2> BindingEnergy(int A, int Z) { // Weizsäcker mass formula

	/*
	Calculates the total Binding energy (B, BE) of nucleus with A, Z. 
	
	Returns a size-2 array with indexing:
	[0] -> total binding energy
	[1] -> binding energy per nucleon
	*/

	int delt = (1 - (A % 2)) * (1 - (2 * (Z % 2)));

	double a_V = 15.40;
	double a_s = 16.71;
	double a_C = 0.701;
	double a_a = 22.56;
	double a_p = 11.88;

	double t1 = a_V * A;
	double t2 = a_s * pow(A, 2.0 / 3.0);
	double t3 = a_C * Z * (Z - 1) * pow(A, - 1.0 / 3.0);
	double t4 = a_a * pow(A - 2.0 * Z, 2) / A;
	double t5 = a_p * delt * pow(A, - 1.0 / 2.0);


	double B_tot = t1 - t2 - t3 - t4 + t5;
	double BpN = B_tot / A;

	return { B_tot, BpN };
}





