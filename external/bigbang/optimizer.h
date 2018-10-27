#ifndef OPTIMIZER_H
#define OPTIMIZER_H

#include <memory>
#include <Eigen/Sparse>
#include <optimization.h>
#include <boost/property_tree/ptree.hpp>
#include <chrono>
#include <Eigen/CholmodSupport>

namespace bigbang {

template <typename T>
class Functional;

template <typename T>
class Constraint;

struct opt_args {
  size_t max_iter;
  double eps;
  bool lineseach;
};

using pfunc=std::shared_ptr<Functional<double>>;
using pcons=std::shared_ptr<Constraint<double>>;

class optimize_callback
{
public:
  virtual int operator() (const size_t iter, const double value, const double error) = 0;
  std::vector<double> bufferA_, bufferB_;
  std::chrono::time_point<std::chrono::system_clock> tic_, toc_;
};

class uncons_newton_solver
{
public:
  uncons_newton_solver(const std::shared_ptr<Functional<double>> &func,
                       boost::property_tree::ptree &opts);
  int optimize(double *x, optimize_callback *cb=NULL);
protected:
  const size_t dim_;
  const std::shared_ptr<Functional<double>> &func_;
  boost::property_tree::ptree &opts_;
  Eigen::SparseMatrix<double> H_;
  //  Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower|Eigen::Upper> solver_;
  Eigen::CholmodSimplicialLDLT<Eigen::SparseMatrix<double>> llt_solver_;
};

int newton_solve(double *x, const size_t dim, const pfunc &f, boost::property_tree::ptree &opts);

int newton_solve_with_constrained_dofs(double *x, const size_t dim, const pfunc &f, const std::vector<size_t> &g2l, const opt_args &args);

int lbfgs_solve(double *x, const size_t dim, const pfunc &f, const opt_args &args);

int constrained_newton_solve(double *x, const size_t dim, const pfunc &f, const pcons &c, const opt_args &args);

int gauss_newton_solve(double *x, const size_t dim, const pcons &f);

int apply_jacobi(const Eigen::SparseMatrix<double, Eigen::RowMajor> &A, const Eigen::VectorXd &rhs, Eigen::VectorXd &x);

int apply_gauss_seidel(const Eigen::SparseMatrix<double, Eigen::RowMajor> &A, const Eigen::VectorXd &rhs, Eigen::VectorXd &x, bool increase=true);

// apply gauss seidel to symmetric sparse matrix
void apply_gauss_seidel_sym(const Eigen::SparseMatrix<double> &A,
                            const Eigen::VectorXd &rhs,
                            Eigen::VectorXd &x,
                            const bool increase=true);

// apply damped jacobi to symmetric sparse matrix
void apply_damped_jacobi_sym(const Eigen::SparseMatrix<double> &A,
                             const Eigen::VectorXd &rhs,
                             Eigen::VectorXd &x,
                             Eigen::VectorXd &xtmp,
                             const double w=0.6666666666666666666);

typedef void (*alglib_rep_func)(const alglib::real_1d_array &x, double func, void *ptr);

struct alglib_solver_args
{
  double epsf_, epsx_, epsg_;
  size_t maxits_;
  bool precscale_;
  alglib_solver_args()
      : epsf_(0), epsx_(0), epsg_(0),
        maxits_(0), precscale_(false) {}
};

// Bleic solver with boundary constraints
int alglib_bleic_bc_solve(const std::shared_ptr<Functional<double>> &func,
                          const size_t dim, double *X,
                          const double *lbnd, const double *ubnd,
                          const alglib_solver_args &args,
                          alglib_rep_func rep_cb=NULL);

// Bleic solver with general LINEAR constraints
int alglib_bleic_lc_solve(const std::shared_ptr<Functional<double>> &func,
                          const std::shared_ptr<Constraint<double>> &cons,
                          const size_t dim, double *X,
                          const alglib_solver_args &args,
                          alglib_rep_func rep_cb=NULL);

int alglib_lbfgs_solve(const std::shared_ptr<Functional<double>> &func,
                       const size_t dim, double *X,
                       const alglib_solver_args &args,
                       alglib_rep_func rep_cb=NULL,
                       const double *prec=0,
                       const char option='D');

}
#endif
