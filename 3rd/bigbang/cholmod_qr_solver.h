#ifndef CHOLMOD_QR_SOLVER_H
#define CHOLMOD_QR_SOLVER_H

#include <Eigen/Sparse>

namespace bigbang {

int cholmod_qr_solve(const Eigen::SparseMatrix<double> &A, const Eigen::VectorXd &b, Eigen::VectorXd &x);

}

#endif
