#pragma once
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <type_traits>

std::vector<double> range(double start, double stop, int points);
Eigen::MatrixXd IdMat(int N);
Eigen::MatrixXd Laplacian(int N, double dx);


// template: works for both std::vector<double> and Eigen::VectorXd
template <typename VectorType>
Eigen::SparseMatrix<double> diagSparse(const VectorType& vec)
{
    const int N = static_cast<int>(vec.size());
    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(N);

    for (int i = 0; i < N; ++i)
    {
        if constexpr (std::is_same_v<VectorType, Eigen::VectorXd>)
        {
            triplets.emplace_back(i, i, vec(i));  // Eigen indexing
        }
        else
        {
            triplets.emplace_back(i, i, vec[i]);  // std::vector indexing
        }
    }

    Eigen::SparseMatrix<double> D(N, N);
    D.setFromTriplets(triplets.begin(), triplets.end());

    return D;
}
