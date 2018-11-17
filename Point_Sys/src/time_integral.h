#ifndef TIME_INTEGRAL_H
#define TIME_INTEGRAL_H

#include <vector>

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Eigen/IterativeLinearSolvers>

typedef Eigen::MatrixXd DeM_d;
typedef Eigen::SparseMatrix<double> SpM_d;
typedef std::vector<Eigen::Triplet<double>> Tris_d;

int implicit_PCG(const SpM_d &M, const Spm_d &K, DeM_d )
#endif
