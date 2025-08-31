#pragma once

#include <functional>
#include <vector>
#include <Eigen/Sparse>


// "Potential term" type: takes r-grid -> sparse diagonal matrix
using PotentialTerm = std::function<Eigen::SparseMatrix<double>(const std::vector<double>& r)>;


// Radial potential functions (return: raw vectors)
std::vector<double> V_HO(const std::vector<double>& r);
std::vector<double> V_WS(const std::vector<double>& r, int A, int Z);
std::vector<double> V_C(const std::vector<double>& r, int A, int Z);
//std::vector<double> V_SO(const std::vector<double>& r, int l, double j, int A, int Z);


// Wrappers that turn raw functions into PotentialTerms
PotentialTerm HO_potential();
PotentialTerm WS_potential(int A, int Z);
PotentialTerm Coulomb_potential(int A, int Z);
//PotentialTerm SO_potential(int l, double j, int A, int Z);



