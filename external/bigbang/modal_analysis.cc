#include "modal_analysis.h"

#include <iostream>
#include <Eigen/UmfPackSupport>
#include <Eigen/CholmodSupport>

#include "arpaca.h"
#include "bigbang/config.h"
#include "util.h"
#include "def.h"
#include "energy.h"

using namespace std;
using namespace Eigen;
using namespace arpaca;

namespace bigbang {

//int gevp_solver::solve(const spmatd_t &A, const spmatd_t &B, const size_t num,
//                       const EigenvalueType mode, MatrixXd &U, VectorXd &lambda) {
//  const size_t dim = A.cols();
//  U.resize(dim, num);
//  lambda.resize(num);

//  SimplicialLLT<spmatd_t> llt_solver;
//  UmfPackLU<spmatd_t> lu_solver;
//  spmatd_t L, Linv, C, Id(dim, dim);
//  Id.setIdentity();

//  llt_solver.compute(B);
//  ASSERT(llt_solver.info() == Success);
//  L = llt_solver.matrixL();

//  lu_solver.compute(L);
//  ASSERT(lu_solver.info() == Success);
//  Linv = lu_solver.solve(Id);

//  C = Linv*A*Linv.transpose();
//  SymmetricEigenSolver<double> eig_solver;
//  eig_solver.SetEigenvalueType(static_cast<arpaca::EigenvalueType>(mode));
//  eig_solver.SetMaxIterations(10000);
//  eig_solver.SetNumLanczosVectors(0);
//  eig_solver.SetTolerance(-1);
//  eig_solver.Solve(C.cols(), num, MakeDefaultOperator(C));
//  printf("[INFO] arpack %d iter, %d converged, %s\n",
//         eig_solver.num_actual_iterations(), eig_solver.num_converged_eigenvalues(), eig_solver.GetInfo());

//  U = eig_solver.eigenvectors();
//  lambda = eig_solver.eigenvalues();
//  U = (Linv.transpose()*U).eval();

//  return 0;
//}

basis_builder::basis_builder(const spmatd_t &M, const spmatd_t &K, const unordered_set<size_t> &fixDOF)
  : dim_(M.cols()) {
  cout << "[info] init basis builder\n";
  SimplicialLLT<spmatd_t> llt_solver;
  UmfPackLU<spmatd_t> lu_solver;
  spmatd_t L, Id(dim_, dim_);
  Id.setIdentity();

  llt_solver.compute(M);
  ASSERT(llt_solver.info() == Success);
  L = llt_solver.matrixL();

  lu_solver.compute(L);
  ASSERT(lu_solver.info() == Success);
  Linv_ = lu_solver.solve(Id);

  C_ = Linv_*K*Linv_.transpose();

  g2l_.resize(dim_);
  size_t ptr = 0;
  for (size_t i = 0; i < dim_; ++i) {
    if ( fixDOF.find(i) != fixDOF.end() )
      g2l_[i] = -1;
    else
      g2l_[i] = ptr++;
  }
  rm_spmat_col_row(C_, g2l_);
}

int basis_builder::compute(const size_t num, MatrixXd &U, VectorXd &lambda) const {
  cout << "[info] compute basis" << endl;
  SymmetricEigenSolver<double> solver;
  solver.SetEigenvalueType(arpaca::ALGEBRAIC_SMALLEST);
  solver.SetMaxIterations(20000);
  solver.SetNumLanczosVectors(0);
  solver.SetTolerance(-1);
  solver.Solve(C_.cols(), num, MakeDefaultOperator(C_));
  printf("[INFO] arpack %d iter, %d converged, %s\n", solver.num_actual_iterations(), solver.num_converged_eigenvalues(), solver.GetInfo());

  MatrixXd Ut = solver.eigenvectors();
  lambda = solver.eigenvalues();

  U.resize(dim_, num);
#pragma omp parallel for
  for (size_t i = 0; i < U.rows(); ++i) {
    if ( g2l_[i] == -1 )
      U.row(i).setZero();
    else
      U.row(i) = Ut.row(g2l_[i]);
  }
  U = Linv_.transpose()*U.eval();
  return 0;
}

template <class Mat, class Vec>
static void contract(const vector<Mat> &T, const Vec &v, Mat &Tv) {
  ASSERT(T.size() == v.size());
  Tv.resize(T.front().rows(), T.front().cols());
  Tv.setZero();

  for (size_t i = 0; i < T.size(); ++i)
    Tv += T[i]*v[i];
}

template <class Mat, class Vec1, class Vec2, class Vec3>
static inline void project_with_metric(const Mat &M, const Vec1 &v, const Vec2 &u, Vec3 &vu) {
  vu = v.dot(M*u)*u/u.dot(M*u);
}

static void mass_gram_schmidt(const SparseMatrix<double> &M, const MatrixXd &U,
                              MatrixXd &Un) {
  Un = U;
  VectorXd tm;
  for (size_t i = 0; i < U.cols(); ++i) {
    for (size_t j = 0; j < i; ++j) {
      project_with_metric(M, U.col(i), Un.col(j), tm);
      Un.col(i) -= tm;
    }
    const double len = (Un.col(i).dot(M*Un.col(i)));
    Un.col(i) /= sqrt(len);
  }
}

// solve modal derivative of StVK materials
int basis_builder::extend(const shared_ptr<Functional<double>> &energy,
                          const double *x0, const SparseMatrix<double> &M, const MatrixXd &U,
                          MatrixXd &Ue) {
  ASSERT(dynamic_pointer_cast<elastic_potential>(energy)->type_ == elastic_potential::STVK);
  cout << "[INFO] solve modal derivatives..." << endl;
  const size_t dim = energy->Nx();
  ASSERT(dim == U.rows());
  Eigen::Map<const VectorXd> X0(x0, dim);

  const double dx = 0.1;
  vector<SparseMatrix<double>> Hijk(dim);

  #pragma omp parallel for
  for (size_t i = 0; i < dim; ++i) {
    VectorXd p = X0;

    p[i] += dx;
    SparseMatrix<double> K1(dim, dim); {
      vector<Triplet<double>> trips;
      energy->Hes(p.data(), &trips);
      K1.setFromTriplets(trips.begin(), trips.end());
    }

    p[i] -= 2*dx;
    SparseMatrix<double> K2(dim, dim); {
      vector<Triplet<double>> trips;
      energy->Hes(p.data(), &trips);
      K2.setFromTriplets(trips.begin(), trips.end());
    }

    Hijk[i] = (K1-K2)/(2*dx);
  }

  cout << "[INFO] solve psi_{ij}..." << endl;
  
  SparseMatrix<double> K(dim, dim); {
    vector<Triplet<double>> trips;
    energy->Hes(x0, &trips);
    K.setFromTriplets(trips.begin(), trips.end());
    cout << "# K(x0) norm: " << K.norm() << endl;
  }
  rm_spmat_col_row(K, g2l_);
  CholmodSimplicialLLT<SparseMatrix<double>> solver;
  solver.compute(K);
  ASSERT(solver.info() == Eigen::Success);

  const size_t orig_basis_num = U.cols();
  Ue = MatrixXd::Zero(dim, orig_basis_num*(orig_basis_num+1)/2);
  #pragma omp parallel for
  for (size_t i = 0; i < orig_basis_num; ++i) {
    SparseMatrix<double> HUi;
    VectorXd rhs, dx, Dx;
    for (size_t j = i; j < orig_basis_num; ++j) {
      const int id = (2*orig_basis_num+1-i)*i/2+j-i;
      contract(Hijk, U.col(i), HUi);
      rhs = -HUi*U.col(j);
      rm_vector_row(rhs, g2l_);
      dx = solver.solve(rhs);
      Dx = VectorXd::Zero(dim);
      rc_vector_row(dx, g2l_, Dx);
      Ue.col(id) += Dx;
    }
  }

  // mass-gram-schimdt normalization
  cout << "[INFO] mass gram schmidt normalization..." << endl;
  
  MatrixXd Un;
  mass_gram_schmidt(M, Ue, Un);
  Ue = Un;
  
  return 0;
}

}
