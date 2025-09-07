#pragma once

#include <vector>
#include "numerics.h"


std::vector<double> ddiff(std::vector<double>& f, double dx);
double simpson(const std::vector<double>& x, const Eigen::VectorXd& y);
Eigen::VectorXd normalize_wavefunction(const std::vector<double>& x, const Eigen::VectorXd& psi);



