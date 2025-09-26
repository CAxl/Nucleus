#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>  // for std::minmax_element
#include <fstream>
#include <Eigen/Dense>


#include <util.h>



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

// overload for Eigen::VectorXd
void write_xy_to_file(const std::vector<double>& x, const Eigen::VectorXd& y, const std::string& filename) {
    if (x.size() != static_cast<size_t>(y.size())) {
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





