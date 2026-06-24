#pragma once
#include <vector>
#include <algorithm>

void print_matrix(const std::vector<double>& matrix);
void write_xy_to_file(const std::vector<double>& x, const std::vector<double>& y, const std::string& filename);
// Overload for Eigen
void write_xy_to_file(const std::vector<double>& x, const Eigen::VectorXd& y, const std::string& filename);