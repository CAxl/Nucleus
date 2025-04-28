#pragma once
#include <vector>
#include <algorithm>

void print_matrix(const std::vector<double>& matrix);
void quick_plot(const std::vector<double>& y, int height = 20);
void write_xy_to_file(const std::vector<double>& x, const std::vector<double>& y, const std::string& filename);
