#pragma once
#include <vector>
#include <Eigen/Dense>

std::vector<double> range(double start, double stop, int points);
Eigen::MatrixXd IdMat(int N);
Eigen::MatrixXd Laplacian(int N, double dx);
