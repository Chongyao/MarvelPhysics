#include "util.h"

#include <iostream>
#include <SymGEigsSolver.h>
#include <MatOp/SparseSymMatProd.h>
#include <Eigen/Geometry>

#include "arpaca.h"

using namespace std;
using namespace Eigen;
using namespace Spectra;

namespace bigbang {

int solve_gen_eig_prob(const Eigen::SparseMatrix<double> &K,
                       const Eigen::SparseMatrix<double> &M,
                       const size_t eignum, const int options,
                       const string &solver,
                       Eigen::VectorXd *eigenvalues,
                       Eigen::MatrixXd *eigenvectors) {
  assert(options == 0 || options == 1);
  
  //-> M = LL^T
  SparseMatrix<double> invL(M.rows(), M.cols());
  for (size_t i = 0; i < M.rows(); ++i)
    invL.coeffRef(i, i) = 1.0/sqrt(M.coeff(i, i));

  const SparseMatrix<double> LMLT = invL*K*invL.transpose();

  if ( solver == "arpaca" ) {
    arpaca::SymmetricEigenSolver<double> solver;
    if ( options == 0 )
      solver.SetEigenvalueType(arpaca::ALGEBRAIC_SMALLEST);
    else if ( options == 1 )
      solver.SetEigenvalueType(arpaca::ALGEBRAIC_LARGEST);
    solver.SetMaxIterations(20000);
    solver.SetNumLanczosVectors(0);
    solver.SetTolerance(-1);
    solver.Solve(LMLT.cols(), eignum, arpaca::MakeDefaultOperator(LMLT));

    if ( eigenvalues )
      *eigenvalues = solver.eigenvalues();
    if ( eigenvectors )
      *eigenvectors = invL.transpose()*solver.eigenvectors();

    return 0;
  } 

  if ( solver == "spectra" ) {
    SparseSymMatProd<double> op(LMLT);

    const int maxits = 20000;
    const double tolerance = 1e-10;
    const int ncv = min((int)LMLT.cols(), 2*(int)eignum);
    
    if ( options == 0 ) {
      Spectra::SymEigsSolver<double, SMALLEST_ALGE, SparseSymMatProd<double>> solver(&op, eignum, ncv);
      solver.init();
      int nconv = solver.compute(maxits, tolerance, SMALLEST_ALGE);

      if ( eigenvalues )
        *eigenvalues = solver.eigenvalues();
      if ( eigenvectors )
        *eigenvectors = invL.transpose()*solver.eigenvectors();
    } else if ( options == 1 ) {
      Spectra::SymEigsSolver<double, LARGEST_ALGE, SparseSymMatProd<double>> solver(&op, eignum, ncv);
      solver.init();
      int nconv = solver.compute(maxits, tolerance, SMALLEST_ALGE);

      if ( eigenvalues )
        *eigenvalues = solver.eigenvalues();
      if ( eigenvectors )
        *eigenvectors = invL.transpose()*solver.eigenvectors();
    }

    return 0;
  }

  return __LINE__;
}

int solve_sym_eig_prob(const Eigen::SparseMatrix<double> &K,
                       const size_t eignum, const int options,
                       const string &solver,
                       Eigen::VectorXd *eigenvalues,
                       Eigen::MatrixXd *eigenvectors) {
  assert(options == 0 || options == 1);
  
  if ( solver == "arpaca" ) {
    arpaca::SymmetricEigenSolver<double> solver;
    if ( options == 0 )
      solver.SetEigenvalueType(arpaca::ALGEBRAIC_SMALLEST);
    else if ( options == 1 )
      solver.SetEigenvalueType(arpaca::ALGEBRAIC_LARGEST);
    solver.SetMaxIterations(20000);
    solver.SetNumLanczosVectors(0);
    solver.SetTolerance(-1);
    solver.Solve(K.cols(), eignum, arpaca::MakeDefaultOperator(K));

    if ( eigenvalues )
      *eigenvalues = solver.eigenvalues();
    if ( eigenvectors )
      *eigenvectors = solver.eigenvectors();

    return 0;
  } 

  if ( solver == "spectra" ) {
    SparseSymMatProd<double> op(K);

    const int maxits = 20000;
    const double tolerance = 1e-10;
    const int ncv = min((int)K.cols(), 2*(int)eignum);
    
    if ( options == 0 ) {
      Spectra::SymEigsSolver<double, SMALLEST_ALGE, SparseSymMatProd<double>> solver(&op, eignum, ncv);
      solver.init();
      int nconv = solver.compute(maxits, tolerance, SMALLEST_ALGE);

      if ( eigenvalues )
        *eigenvalues = solver.eigenvalues();
      if ( eigenvectors )
        *eigenvectors = solver.eigenvectors();
    } else if ( options == 1 ) {
      Spectra::SymEigsSolver<double, LARGEST_ALGE, SparseSymMatProd<double>> solver(&op, eignum, ncv);
      solver.init();
      int nconv = solver.compute(maxits, tolerance, SMALLEST_ALGE);

      if ( eigenvalues )
        *eigenvalues = solver.eigenvalues();
      if ( eigenvectors )
        *eigenvectors = solver.eigenvectors();
    }

    return 0;
  }

  return __LINE__;
}

int stable_polor_decomp(double *dF) {
  Matrix3d A(dF);
  Quaterniond q = Quaterniond(A);
  q.normalize();
  for(size_t iter=0;iter<5;iter++) {
    Matrix3d R=q.matrix();
    Vector3d omega=(R.col(0).cross(A.col(0))+R.col(1).cross(A.col(1))+R.col(2).cross(A.col(2)))*(1.0/fabs(R.col(0).dot(A.col(0))+R.col(1).dot(A.col(1))+R.col(2).dot(A.col(2)))+1e-9);
    double w=omega.norm();
    if ( w < 1e-9 ) {
      break;          
    }
    q=Quaterniond(AngleAxisd(w,(1.0/w)*omega))*q;
    q.normalize();
  }
  Map<Matrix3d> R(dF);
  R = q.matrix();
  return 0;  
}

}
