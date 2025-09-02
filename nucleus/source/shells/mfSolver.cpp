#include "mfSolver.h"
#include "hamil.h"
#include "potentials.h"
#include <Spectra/SymEigsSolver.h>
#include <Spectra/MatOp/SparseSymMatProd.h>
#include <iostream>

//--------------------------------------
// Solve MF Hamiltonian for given orbital
//--------------------------------------
Eigen::VectorXd solve_mean_field(
    int A, int Z,
    int l, double j,
    const std::vector<double>& r, double dx,
    bool is_proton,
    int n_eigs
)
{
    // build potentials
    std::vector<PotentialTerm> V;
    V.push_back(WS_potential(A, Z));
    if (is_proton) V.push_back(Coulomb_potential(A, Z));    // only for protons
    V.push_back(SO_potential(l, j, A, Z));

    // build Hamil
    Eigen::SparseMatrix<double> H = H_nl(l, r, dx, V);

    // Spectra eigensolver
    int m = 10 * n_eigs;    // krylov subspace dimension
    
    Spectra::SparseSymMatProd<double> op(H);
    Spectra::SymEigsSolver<Spectra::SparseSymMatProd<double>> eigs(op, n_eigs, m);

    eigs.init();
    int nconv = eigs.compute(Spectra::SortRule::SmallestAlge);

    if (eigs.info() != Spectra::CompInfo::Successful)
        std::cerr << "Eigensolver failed!" << std::endl;


    return eigs.eigenvalues();

}