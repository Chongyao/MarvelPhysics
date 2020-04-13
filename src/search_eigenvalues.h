#ifndef MARVEL_SEARCH_EIGENVALUES
#define MARVEL_SEARCH_EIGENVALUES

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <iostream>

namespace marvel{
int find_max_min_eigenvalues(const Eigen::SparseMatrix<double>& A, double& max_eigvalue, double& min_eigvalue);


double find_max_eigenvalue(const Eigen::SparseMatrix<double>& A, const size_t max_itrs = 1000);

double find_min_eigenvalue(const Eigen::SparseMatrix<double>& A, const size_t max_itrs = 1000);

double find_min_eigenvalue(const Eigen::SparseMatrix<double>& A, const double& max_eig, const size_t max_iters);
double find_condition_number(const Eigen::SparseMatrix<double>& A, const size_t max_itrs = 1000);
}
#endif
