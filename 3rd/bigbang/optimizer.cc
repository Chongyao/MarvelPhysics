#include "optimizer.h"

#include <iostream>
#include <iomanip>
#include <Eigen/UmfPackSupport>
#include <Eigen/CholmodSupport>
#include <lbfgs.h>
#include <zjucad/linear_solver/linear_solver.h>

#include "config.h"
#include "def.h"
#include "util.h"

using namespace std;
using namespace Eigen;
using namespace alglib;

namespace bigbang {

uncons_newton_solver::uncons_newton_solver(const shared_ptr<Functional<double>> &func,
                                           boost::property_tree::ptree &opts)
    : func_(func), opts_(opts), dim_(func->Nx()) {
  H_.resize(dim_, dim_);
}

int uncons_newton_solver::optimize(double *x, optimize_callback *cb) {
  Map<VectorXd> X(x, dim_);
  VectorXd xstar = X, dx(dim_), grad(dim_), xnew(dim_);

  const size_t maxits = opts_.get<size_t>("maxits.value", 1000);
  const bool need_linesearch = opts_.get<bool>("linesearch.value", true);
  const double epsg = opts_.get<double>("epsg.value", 1e-8);

  size_t iter = 0;
  
  for (iter = 0; iter < maxits; ++iter) {
    double value = 0; {
      func_->Val(&xstar[0], &value);
    }
    grad.setZero(); {
      func_->Gra(&xstar[0], &grad[0]);
    }
    H_.resize(dim_, dim_); {
      vector<Triplet<double>> trips;
      func_->Hes(&xstar[0], &trips);
      H_.reserve(trips.size());
      H_.setFromTriplets(trips.begin(), trips.end());
      H_.makeCompressed();
    }
    llt_solver_.compute(H_);
    ASSERT(llt_solver_.info() == Eigen::Success);
    dx = -llt_solver_.solve(grad);
    ASSERT(llt_solver_.info() == Eigen::Success);    

    double h = 1.0;
    // line search here
    if ( need_linesearch ) {
      h /= 0.5;
      double lhs = 0.0, rhs = 0.0;
      do {
        h *= 0.5;
        xnew = xstar+h*dx;
        lhs = 0;
        func_->Val(&xnew[0], &lhs);
        rhs = value+h*0.45*(grad.dot(dx));
      } while ( lhs >= rhs && h > 1e-12 );
    }
    xstar += h*dx;

    const double nlp_error = grad.lpNorm<Infinity>();
    if ( cb )
      if ( (*cb)(iter, value, nlp_error) )
        return __LINE__;
    
    if ( nlp_error < epsg ) {
      cout << "\t@CONVERGED" << endl;
      break;
    }
    if ( h <= 1e-12 ) {
      cout << "\t@line search stall." << endl;
      break;
    }
  }
  X = xstar;
  cout << "Number of Iterations....:  " << iter << endl;
  return 0;
}

int newton_solve(double *x, const size_t dim, const pfunc &f, boost::property_tree::ptree &opts) {
  if ( dim != f->Nx() ) {
    cerr << "[error] dim not match\n";
    return __LINE__;
  }
  
  Map<VectorXd> X(x, dim);
  VectorXd xstar = X, dx(dim), grad(dim), xnew(dim);
  SparseMatrix<double> H(dim, dim);

  const size_t maxits = opts.get<size_t>("maxits.value", 1000);
  const bool need_linesearch = opts.get<bool>("linesearch.value", true);
  const double epsg = opts.get<double>("epsg.value", 1e-8);

  for (size_t iter = 0; iter < maxits; ++iter) {
    double value = 0; {
      f->Val(&xstar[0], &value);
      if ( iter % 100 == 0 ) {
        cout << "\t@iter " << iter << endl;
        cout << "\t@energy value: " << value << endl;
      }
    }
    grad.setZero(); {
      f->Gra(&xstar[0], &grad[0]);
      if ( grad.norm() <= epsg ) {
        cout << "\t@GRADIENT CONVERGED\n\n";
        break;
      }
    }
    H.resize(dim, dim); {
      vector<Triplet<double>> trips;
      f->Hes(&xstar[0], &trips);
      H.reserve(trips.size());
      H.setFromTriplets(trips.begin(), trips.end());
      H.makeCompressed();
    }
    ConjugateGradient<SparseMatrix<double>, Lower|Upper> sol;
    //    CholmodSimplicialLLT<SparseMatrix<double>> sol;
    sol.compute(H);
    ASSERT(sol.info() == Eigen::Success);
    dx = -sol.solve(grad);
    ASSERT(sol.info() == Eigen::Success);    
    double xstar_norm = xstar.norm();
    double h = 1.0;
    // line search here
    if ( need_linesearch ) {
      h /= 0.5;
      double lhs = 0.0, rhs = 0.0;
      do {
        h *= 0.5;
        xnew = xstar+h*dx;
        lhs = 0;
        f->Val(&xnew[0], &lhs);
        rhs = value+h*0.45*(grad.dot(dx));
      } while ( lhs >= rhs && h > 1e-12 );
    }
    xstar += h*dx;
    if ( h <= 1e-12 ) {
      cout << "\t@line search stall." << endl;
      break;
    } else if ( h*dx.norm() <= epsg*xstar_norm ) {
      cout << "\t@CONVERGED after " << iter << " iterations\n\n";
      break;
    }
  }
  X = xstar;
  return 0;
}

int newton_solve_with_constrained_dofs(double *x, const size_t dim, const pfunc &f, const vector<size_t> &g2l, const opt_args &args) {
  if ( dim != f->Nx() ) {
    cerr << "[error] dim not match\n";
    return __LINE__;
  }
  // SimplicialCholesky<SparseMatrix<double>> sol;
  // sol.setMode(SimplicialCholeskyLDLT);
  CholmodSimplicialLLT<SparseMatrix<double>> sol;
  Map<VectorXd> X(x, dim);
  VectorXd xstar = X, dx, Dx(dim), xnew(dim);
  
  for (size_t iter = 0; iter < args.max_iter; ++iter) {
    double value = 0; {
      f->Val(&xstar[0], &value);
      if ( iter % 100 == 0 ) {
        cout << "\t@iter " << iter << endl;
        cout << "\t@energy value: " << value << endl;
      }
    }
    VectorXd grad = VectorXd::Zero(dim); {
      f->Gra(&xstar[0], &grad[0]);
      grad *= -1.0;
      if ( grad.norm() <= args.eps ) {
        cout << "\t@GRADIENT CONVERGED\n\n";
        break;
      }
    }
    SparseMatrix<double> H(dim, dim); {
      vector<Triplet<double>> trips;
      f->Hes(&xstar[0], &trips);
      H.reserve(trips.size());
      H.setFromTriplets(trips.begin(), trips.end());
    }

    rm_spmat_col_row(H, g2l);
    rm_vector_row(grad, g2l);
    sol.compute(H);
    ASSERT(sol.info() == Eigen::Success);
    dx = sol.solve(grad);
    ASSERT(sol.info() == Eigen::Success);
    Dx.setZero();
    rc_vector_row(dx, g2l, Dx);

    double xstar_norm = xstar.norm();
    double h = 1.0;
    // line search here
    if ( args.lineseach ) {
      h /= 0.5;
      double lhs = 0.0, rhs = 0.0;
      do {
        h *= 0.5;
        xnew = xstar+h*Dx;
        lhs = 0;
        f->Val(&xnew[0], &lhs);
        rhs = value+h*0.45*(-grad.dot(Dx));
      } while ( lhs >= rhs && h > 1e-12 );
    }
    xstar += h*Dx;
    if ( h <= 1e-12 ) {
      cout << "\t@line search stall." << endl;
      break;
    } else if ( h*Dx.norm() <= args.eps*xstar_norm ) {
      cout << "\t@CONVERGED after " << iter << " iterations\n\n";
      break;
    }
  }
  X = xstar;
  return 0;
}

// static shared_ptr<Functional<double>> energy;

static lbfgsfloatval_t evaluate(void *instance, const lbfgsfloatval_t *x,
                                lbfgsfloatval_t *g, const int n,
                                const lbfgsfloatval_t step) {
  Functional<double>* energy = static_cast<Functional<double>*>(instance);
  double f = 0;
  energy->Val(x, &f);
  std::fill(g, g+n, 0);
  energy->Gra(x, g);
  return f;
}

static int progress(void *instance,
                    const lbfgsfloatval_t *x,
                    const lbfgsfloatval_t *g,
                    const lbfgsfloatval_t fx,
                    const lbfgsfloatval_t xnorm,
                    const lbfgsfloatval_t gnorm,
                    const lbfgsfloatval_t step,
                    int n, int k, int ls) {
  if ( k % 100 == 0 )
    cout << "----- LBFGS: count (" << k << "), energy ("
         << fx << "), gnorm (" << gnorm << ")" << endl;
  return 0;
}

int lbfgs_solve(double *x, const size_t dim, const pfunc &f, const opt_args &args) {
  if ( dim != f->Nx() ) {
    cerr << "[error] dim not match\n";
    return __LINE__;
  }
  
  lbfgsfloatval_t fx;
  lbfgs_parameter_t param;
  lbfgs_parameter_init(&param);
  param.max_iterations = args.max_iter;

  int ret = lbfgs(dim, x, &fx, evaluate, progress, f.get(), &param);

  cout << "L-BFGS optimization terminated with status code = " << ret << endl;
  return 0;
}

int constrained_newton_solve(double *x, const size_t dim, pfunc &f, const pcons &c, const opt_args &args) {
  if ( dim != f->Nx() || dim != c->Nx() ) {
    cerr << "[error] dim not match\n";
    return __LINE__;
  }
  const size_t fdim = c->Nf();
  SimplicialCholesky<SparseMatrix<double>> solver;
  solver.setMode(SimplicialCholeskyLDLT);

  VectorXd xstar = VectorXd::Zero(dim+fdim);
  std::copy(x, x+dim, xstar.data());
  VectorXd xprev = xstar.head(dim);
  for (size_t iter = 0; iter < args.max_iter; ++iter) {
    VectorXd rhs(dim+fdim); {
      rhs.setZero();
      f->Gra(&xstar[0], &rhs[0]);
      c->Val(&xstar[0], &rhs[dim]);
    }
    SparseMatrix<double> LHS(dim+fdim, dim+fdim); {
      vector<Triplet<double>> trips;
      f->Hes(&xstar[0], &trips);
      const auto begin = trips.end();
      c->Jac(&xstar[0], dim, &trips);
      const auto end = trips.end();
      for (auto it = begin; it != end; ++it) {

      }
    }
    solver.compute(LHS);
    ASSERT(solver.info() == Success);
    VectorXd dx = solver.solve(rhs);
    ASSERT(solver.info() == Success);
    xstar += dx;
    if ( (xstar.head(dim)-xprev).norm() <= args.eps*xstar.head(dim).norm() ) {
      cout << "\t@CONVERGED";
      break;
    }
    xprev = xstar.head(dim);
  }
  std::copy(xstar.data(), xstar.data()+dim, x);
  return 0;
}

//int gauss_newton_solve(double *x, const size_t dim, const pcons &f) {
//  if ( dim != f->Nx() ) {
//    return __LINE__;
//  }
//  const size_t fdim = f->Nf();
//  for (size_t i = 0; i < MAX_ITER; ++i) {
//    VectorXd value(fdim); {
//      f->Val(x, &value[0]);
//    }
//    SparseMatrix<double> J(fdim, dim); {

//    }

//  }
//  return 0;
//}

int apply_jacobi(const SparseMatrix<double, RowMajor> &A, const VectorXd &rhs, VectorXd &x) {
  ASSERT(A.rows() == rhs.rows());
  ASSERT(A.cols() == x.rows());
  ASSERT(A.rows() == A.outerSize());
  VectorXd xtemp(x.rows());
#pragma omp parallel for
  for (size_t i = 0; i < A.outerSize(); ++i) {
    double temp = rhs[i];
    double diag = 1.0;
    for (SparseMatrix<double, RowMajor>::InnerIterator it(A, i); it; ++it) {
      if ( static_cast<size_t>(it.col()) == i )
        diag = it.value();
      else
        temp -= it.value()*x[it.col()];
    }
    xtemp[i] = temp/diag;
  }
  x = xtemp;
  return 0;
}

int apply_gauss_seidel(const SparseMatrix<double, RowMajor> &A, const VectorXd &rhs, VectorXd &x, bool increase) {
  ASSERT(A.rows() == rhs.rows());
  ASSERT(A.cols() == x.rows());
  ASSERT(A.rows() == A.outerSize());
  const int64_t rows = A.rows();
  const int64_t begin = increase ? 0 : rows-1;
  const int64_t end = increase ? rows : -1;
  const int delta = increase ? 1 : -1;
  for (int64_t i = begin; i != end; i += delta) {
    double temp = rhs[i];
    double diag = 1.0;
    for (SparseMatrix<double, RowMajor>::InnerIterator it(A, i); it; ++it) {
      if ( static_cast<size_t>(it.col()) == i )
        diag = it.value();
      else
        temp -= it.value()*x[it.col()];
    }
    x[i] = temp/diag;
  }
  return 0;
}

void apply_gauss_seidel_sym(const SparseMatrix<double> &A,
                            const VectorXd &rhs,
                            VectorXd &x,
                            const bool increase) {
  const ptrdiff_t rows = A.rows();
  const ptrdiff_t begin = increase ? 0 : rows-1;
  const ptrdiff_t end = increase ? rows : -1;
  const int delta = increase ? 1 : -1;

  for (ptrdiff_t j = begin; j != end; j += delta) {
    double tmp = rhs[j];
    double diag = 1.0;
    for (SparseMatrix<double>::InnerIterator it(A, j); it; ++it) {
      if ( it.row() == j )
        diag = it.value();
      else
        tmp -= it.value()*x[it.row()];
    }
    x[j] = tmp/diag;
  }
}

void apply_damped_jacobi_sym(const SparseMatrix<double> &A,
                             const VectorXd &rhs,
                             VectorXd &x,
                             VectorXd &xtmp,
                             const double w) {
  for (ptrdiff_t j = 0; j < A.outerSize(); ++j) {
    double temp = rhs[j];
    double diag = 1.0;
    for (SparseMatrix<double>::InnerIterator it(A, j); it; ++it) {
      if ( it.row() == j )
        diag = it.value();
      else
        temp -= it.value()*x[it.row()];
    }
    xtmp[j] = w*temp/diag+(1-w)*x[j];
  }
  x = xtmp;
}

//===============================================================================
static void function1_grad(const real_1d_array &x, double &func, real_1d_array &grad, void *ptr)
{
  Functional<double> *energy = static_cast<Functional<double>*>(ptr);
  
  const size_t dim = energy->Nx();
  const double *ptrx = x.getcontent();
  double *ptrg = grad.getcontent();
  
  func = 0;
  energy->Val(ptrx, &func);

  std::fill(ptrg, ptrg+dim, 0);
  energy->Gra(ptrx, ptrg);
}

int alglib_lbfgs_solve(const shared_ptr<Functional<double>> &func,
                       const size_t dim, double *X,
                       const alglib_solver_args &args,
                       alglib_rep_func rep_cb,
                       const double *prec,
                       const char option) {
  ASSERT(func.get() && dim == func->Nx());
  
  real_1d_array x;
  x.setcontent(dim, X);
  real_1d_array diag;
  real_2d_array upperL;

  minlbfgsstate state;
  minlbfgsreport rep;
    
  alglib::minlbfgscreate(dim, 3, x, state);
  alglib::minlbfgssetcond(state, args.epsg_, args.epsf_, args.epsx_, args.maxits_);
  alglib::minlbfgssetxrep(state, true);

  if ( args.precscale_ )
    alglib::minlbfgssetprecscale(state);

  if ( prec ) {
    if ( option == 'D' ) {
      //-> Diagonal preconditioner
      diag.setcontent(dim, prec);
      alglib::minlbfgssetprecdiag(state, diag);
    } else if ( option == 'L' ) {
      //-> cholesky preconditioner
      // input prec is the flattened lower triangular matrix in column major order
      // output upperL is row major upper triangular matrix
      upperL.setcontent(dim, dim, prec);
      const bool is_upper = true;
      alglib::minlbfgssetpreccholesky(state, upperL, is_upper);
    } else {
      cerr << "[!] invalid preconditioner type." << endl;
    }
  }

  alglib::minlbfgsoptimize(state, function1_grad, rep_cb, func.get());
  alglib::minlbfgsresults(state, x, rep);

  const double *ptrx = x.getcontent();
  std::copy(ptrx, ptrx+dim, X);

  cout << "[INFO] lbfgs status: " << int(rep.terminationtype) << endl; // EXPECTED: 4
  cout << "[INFO] lbfgs actual iter: " << int(rep.iterationscount) << endl;
  return 0;
}

int alglib_bleic_bc_solve(const shared_ptr<Functional<double>> &func,
                          const size_t dim, double *X,
                          const double *lbnd, const double *ubnd,
                          const alglib_solver_args &args,
                          alglib_rep_func rep_cb) {
  ASSERT(dim == func->Nx());
  ASSERT(func.get());
  
  real_1d_array x;
  x.setcontent(dim, X);
  real_1d_array bndl;
  bndl.setcontent(dim, lbnd);
  real_1d_array bndu;
  bndu.setcontent(dim, ubnd);

  minbleicstate state;
  minbleicreport rep;
  
  minbleiccreate(x, state);
  minbleicsetbc(state, bndl, bndu);
  minbleicsetcond(state, args.epsg_, args.epsf_, args.epsx_, args.maxits_);
  minbleicsetxrep(state, true);
  if ( args.precscale_ )
    alglib::minbleicsetprecscale(state);
  alglib::minbleicoptimize(state, function1_grad, rep_cb, func.get());
  minbleicresults(state, x, rep);

  const double *ptrx = x.getcontent();
  std::copy(ptrx, ptrx+dim, X);

  cout << "# Bleic status: " << int(rep.terminationtype) << endl; // EXPECTED: 4
  VectorXd glast = VectorXd::Zero(dim);
  func->Gra(X, glast.data());
  cout << "# Bleic gradient norm: " << glast.norm() << endl;
  
  return 0;
}

int alglib_bleic_lc_solve(const shared_ptr<Functional<double>> &func,
                          const shared_ptr<Constraint<double>> &cons,
                          const size_t dim, double *X,
                          const alglib_solver_args &args,
                          alglib_rep_func rep_cb) {
  ASSERT(dim == func->Nx() && dim == cons->Nx());
  ASSERT(func.get() && cons.get());
  
  real_1d_array x;
  x.setcontent(dim, X);

  const size_t K = cons->Nf();
  
  real_2d_array lc; {
    //-> get target constraint value
    VectorXd cv = VectorXd::Zero(K), tx = VectorXd::Zero(dim);
    cons->Val(tx.data(), cv.data());

    //-> get Jacobian
    SparseMatrix<double> J(K, dim);
    vector<Triplet<double>> trips;
    cons->Jac(X, 0, &trips);
    J.setFromTriplets(trips.begin(), trips.end());

    Matrix<double, -1, -1, RowMajor> C(K, dim+1);
    C.setZero();
    for (int j = 0; j < J.outerSize(); ++j) {
      for (SparseMatrix<double>::InnerIterator it(J, j); it; ++it) {
        C(it.row(), it.col()) = it.value();
      }
    }
    C.col(C.cols()-1) = -cv;
        
    lc.setcontent(K, dim+1, C.data()); 
  }

  integer_1d_array ct;
  ct.setlength(K);
  std::fill(ct.getcontent(), ct.getcontent()+ct.length(), 0);
  
  minbleicstate state;
  minbleicreport rep;
  
  minbleiccreate(x, state);
  minbleicsetlc(state, lc, ct);
  minbleicsetcond(state, args.epsg_, args.epsf_, args.epsx_, args.maxits_);
  minbleicsetxrep(state, true);
  if ( args.precscale_ )
    alglib::minbleicsetprecscale(state);
  alglib::minbleicoptimize(state, function1_grad, rep_cb, func.get());
  minbleicresults(state, x, rep);

  const double *ptrx = x.getcontent();
  std::copy(ptrx, ptrx+dim, X);

  cout << "# Bleic status: " << int(rep.terminationtype) << endl; // EXPECTED: 4
  VectorXd glast = VectorXd::Zero(dim);
  func->Gra(X, glast.data());
  cout << "# Bleic gradient norm: " << glast.norm() << endl;
  
  return 0;
}

}
