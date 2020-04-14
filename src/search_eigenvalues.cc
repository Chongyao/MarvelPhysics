#include "search_eigenvalues.h"
#include <Spectra/GenEigsSolver.h>
#include <Spectra/MatOp/SparseGenMatProd.h>
#include <Spectra/MatOp/SparseGenRealShiftSolve.h>
#include <Spectra/GenEigsRealShiftSolver.h>

using namespace std;
using namespace Eigen;
using namespace Spectra;

namespace marvel{

double find_max_eigenvalue(const Eigen::SparseMatrix<double>& A, const size_t max_itrs){
  VectorXd v = VectorXd::Random(A.cols());
  // v = v / v.norm();

  // double last_lambda = 1e40, new_lambda = 0;
  // for(size_t i = 0; i < max_itrs; ++i){
  //   v = A * v;
  //   new_lambda = v.norm();
  //   if(fabs((new_lambda - last_lambda)) > 1e-20)
  //     last_lambda = new_lambda;
  //   else{
  //     v = v / new_lambda;
  //     break;
  //   }
  //   v = v / new_lambda;
  // }
  double last_lambda = 1e40, new_lambda = 0;
  VectorXd w(v.size());
  for(size_t i = 0; i < max_itrs; ++i){
    w = A * v;
    new_lambda = v.dot(w);
    if(fabs(new_lambda - last_lambda) > 1e-20)
      last_lambda = new_lambda;
    else
      break;

    v.array() = w.array() -  w.sum() / w.size();
    v.array() /= v.norm();
  }

  return new_lambda;
}

double find_min_eigenvalue(const Eigen::SparseMatrix<double>& A, const size_t max_itrs){
  ConjugateGradient<SparseMatrix<double>, Lower|Upper> cg;
  cg.compute(A);

  VectorXd v = VectorXd::Ones(A.cols());
  v = v / v.norm();

  double last_v = 1e40, new_v = 0;
  for(size_t i = 0; i < max_itrs; ++i){
    v = cg.solve(v);
    new_v = v.norm();
    if(fabs(new_v - last_v) > 1e-20)
      last_v = new_v;
    else
      break;
    v = v / new_v;
  }
  return 1.0 / new_v;
}

double find_min_eigenvalue(const Eigen::SparseMatrix<double>& A, const double& max_eig, const size_t max_iters){
  SparseMatrix<double> B = A;
 #pragma omp parallel for
  for (int k=0; k<A.outerSize(); ++k)
    for (SparseMatrix<double>::InnerIterator it(B,k); it; ++it)
      if(it.index() == k){
        it.valueRef() -= max_eig;
        break;
      }
  double min_eig = find_max_eigenvalue(B);
  min_eig = fabs(max_eig) - fabs(min_eig);
  return min_eig;
}
double find_condition_number(const Eigen::SparseMatrix<double>& A, const size_t max_itrs){
  double
      max_eig = find_max_eigenvalue(A, max_itrs),
      min_eig = find_min_eigenvalue(A, max_eig, max_itrs);
  cout << "max eig value is " << max_eig << " min eig " << min_eig << endl;
  return max_eig / min_eig;
}


int find_max_min_eigenvalues(const Eigen::SparseMatrix<double>& A, double& max_eigvalue, double& min_eigvalue){
    SparseGenMatProd<double> op(A);

    // Construct eigen solver object, requesting the largest three eigenvalues
    GenEigsSolver< double, LARGEST_MAGN, SparseGenMatProd<double> > eigs(&op, 1, 3);

    // Initialize and compute
    eigs.init();
    int nconv = eigs.compute();

    // Retrieve results
    Eigen::VectorXcd evalues;
    if(eigs.info() == SUCCESSFUL)
        evalues = eigs.eigenvalues();



    max_eigvalue = evalues(0).real();

    SparseGenRealShiftSolve<double> op_min(A);
    // Construct eigen solver object with shift 0
    // This will find eigenvalues that are closest to 0
    GenEigsRealShiftSolver< double, LARGEST_MAGN,
                        SparseGenRealShiftSolve<double> > eigs_min(&op_min, 1, 3, 0.0);

    eigs_min.init();
    eigs_min.compute();

    if(eigs_min.info() == SUCCESSFUL)
    {
        Eigen::VectorXd evalues = eigs_min.eigenvalues().real();
        min_eigvalue = evalues(0);
    }

    cout << "max eigvalue " << max_eigvalue << " min eigvalue " << min_eigvalue << endl;
  return 0;
}
}
