#include "FEM/src/interpolate.h"

#include <iostream>

using namespace std;
using namespace Eigen;
using namespace marvel;

using FLOAT_TYPE = float;

int main() {

    Matrix<FLOAT_TYPE, -1, -1> vertex(3, 4);
    Matrix<int, -1, -1> element(4, 1);
    Matrix<FLOAT_TYPE, -1, -1> point(3, 4);
    Eigen::SparseMatrix<FLOAT_TYPE> coeffient;

    vertex.col(0) << 0, 0, 0;
    vertex.col(1) << 1, 0, 0;
    vertex.col(2) << 0, 1, 0;
    vertex.col(3) << 0, 0, 1;

    element << 0, 1, 2, 3;

    point.col(0) << 0.5, 0, 0;
    point.col(1) << 0, 0.5, 0;
    point.col(2) << 0.01, 0.01, 0.01;
    point.col(3) << 0, 0, 0;

    point = vertex;
    for(int i = 0; i < point.cols(); i++)
        point(0, i) += 1.0;

    interp_pts_in_tets<FLOAT_TYPE, 3>(vertex, element, point, coeffient);

    std::cout << "Vertex: " << std::endl << vertex << std::endl;

    std::cout << "Face: " << std::endl << element << std::endl;

    std::cout << "Point: " << std::endl << point << std::endl;

    std::cout << "Coefficient: " << std::endl << coeffient << std::endl;
    
    for(int i = 0; i < point.cols(); i++)
        vertex(0, i) += 2.0;
    
    
    point = coeffient * vertex;
    std::cout << point << std::endl;
    
}
