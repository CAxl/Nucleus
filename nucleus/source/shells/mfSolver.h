#pragma once

#include <vector>
#include <Eigen/Sparse>
#include <Eigen/Dense>

// Forward declare type alias for potentials
#include "potentials.h"

// Solve mean-field Hamiltonian for given l, j, nucleus (A,Z).
Eigen::VectorXd solve_mean_field(
    int A, int Z,
    int l, double j,
    const std::vector<double>& r, double dx,
    bool is_proton,
    int n_eigs = 5
);

// Expand eigenfunctions in HO basis
//std::vector<double> expand_on_HO(
//    const Eigen::MatrixXd& eigenvecs,
//    const std::vector<double>& r,
//    int l,
//    double hw, int nmax
//);



