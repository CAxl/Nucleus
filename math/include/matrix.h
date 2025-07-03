#pragma once
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Sparse>

std::vector<double> range(double start, double stop, int points);
Eigen::MatrixXd IdMat(int N);
Eigen::MatrixXd Laplacian(int N, double dx);



template <typename VectorType> Eigen::SparseMatrix<double> diagSparse(const VectorType& vec)
{
	/*
	Takes an arbitrary vector (std::vector or Eigen::VectorXd) and puts it on a diagonal Sparse matrix

	*/

	const int N = static_cast<int>(vec.size());
	std::vector<Eigen::Triplet<double>> triplets;
	triplets.reserve(N);

	for (int i = 0; i < N; ++i)
		triplets.emplace_back(i, i, vec[i]);

	Eigen::SparseMatrix<double> D(N, N);
	D.setFromTriplets(triplets.begin(), triplets.end());

	return D;
}
