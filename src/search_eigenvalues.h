#ifndef MARVEL_SEARCH_EIGENVALUES
#define MARVEL_SEARCH_EIGENVALUES

#include <Eigen/Core>
#include <Eigen/SparseCore>

#include <iostream>

namespace marvel{
int find_max_min_eigenvalues(const Eigen::SparseMatrix<double>& A, double& max_eigvalue, double& min_eigvalue);
}
#endif
