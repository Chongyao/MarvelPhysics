#include "search_eigenvalues.h"
#include <Spectra/GenEigsSolver.h>
#include <Spectra/MatOp/SparseSymMatProd.h>
#include <Spectra/SymEigsSolver.h>

using namespace std;
using namespace Eigen;
using namespace Spectra;

namespace marvel{
int find_max_min_eigenvalues(const Eigen::SparseMatrix<double, Eigen::RowMajor>& A, double& max_eigvalue, double& min_eigvalue){
  // Construct matrix operation object using the wrapper class SparseGenMatProd
  SparseSymMatProd<double> op(A);

  // Construct eigen solver object, requesting the largest eigenvalues
  SymEigsSolver< double, LARGEST_ALGE, SparseSymMatProd<double> > eigs_max(&op, 1, 2);
  // Initialize and compute
  eigs_max.init();
  int nconv = eigs_max.compute();

  // Retrieve results
  Eigen::VectorXcd evalues;
  if(eigs_max.info() == SUCCESSFUL){
    evalues = eigs_max.eigenvalues();
    max_eigvalue = evalues(0).real();
  }

  SymEigsSolver< double, SMALLEST_ALGE, SparseSymMatProd<double> > eigs_min(&op, 1, 2);
  eigs_min.init();
  nconv = eigs_min.compute();
  if(eigs_min.info() == SUCCESSFUL){
    evalues = eigs_min.eigenvalues();
    min_eigvalue = evalues(0).real();
  }
  return 0;
}
}
