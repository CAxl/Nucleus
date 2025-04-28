#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>  // for std::minmax_element
#include <fstream>


#include "util.h"



void print_matrix(const std::vector<double>& matrix) {

    /*
    param: matrix: takes input matrix/vector and prints it in a nicely formatted way
    */
    
    const int digits_per_row = 4;
    const int precision = 8;
    const int width = 12;

    std::cout << std::fixed << std::setprecision(precision);
    std::cout << "[";

    for (size_t i = 0; i < matrix.size(); ++i) {
        std::cout << std::setw(width) << matrix[i];

        if (i != matrix.size() - 1) {
            std::cout << ", ";
            if ((i + 1) % digits_per_row == 0) {
                std::cout << "\n ";  // newline and indentation
            }
        }
    }

    std::cout << " ]\n";
}


void quick_plot(const std::vector<double>& data, int height) {
    if (data.empty()) {
        std::cout << "Empty data." << std::endl;
        return;
    }

    double max_val = *std::max_element(data.begin(), data.end());
    double min_val = *std::min_element(data.begin(), data.end());

    if (max_val == min_val) {
        // Special case: all values are constant
        for (int h = height; h >= 0; --h) {
            if (h == height / 2) {
                for (size_t i = 0; i < data.size(); ++i) {
                    std::cout << "*";
                }
            }
            else {
                for (size_t i = 0; i < data.size(); ++i) {
                    std::cout << " ";
                }
            }
            std::cout << "\n";
        }
        return;
    }

    // Regular case: non-constant data
    for (int h = height; h >= 0; --h) {
        double threshold = min_val + (max_val - min_val) * h / height;
        for (size_t i = 0; i < data.size(); ++i) {
            if (data[i] >= threshold) {
                std::cout << "*";
            }
            else {
                std::cout << " ";
            }
        }
        std::cout << "\n";
    }
}


void write_xy_to_file(const std::vector<double>& x, const std::vector<double>& y, const std::string& filename) {
    if (x.size() != y.size()) {
        std::cerr << "Vectors x and y must have the same size!" << std::endl;
        return;
    }

    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }

    for (size_t i = 0; i < x.size(); ++i) {
        file << x[i] << " " << y[i] << "\n";
    }

    file.close();
}





