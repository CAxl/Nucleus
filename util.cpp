#include <iostream>
#include <iomanip>
#include <vector>


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


