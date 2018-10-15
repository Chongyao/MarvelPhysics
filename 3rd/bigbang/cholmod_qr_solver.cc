#include "cholmod_qr_solver.h"

#include <iostream>
#include <SuiteSparseQR.hpp>

using namespace std;
using namespace Eigen;

namespace bigbang {

int cholmod_qr_solve(const SparseMatrix<double> &LHS, const VectorXd &rhs, VectorXd &xx) {
  cholmod_common common;
  cholmod_sparse *A;
  cholmod_dense *x, *b;

  cholmod_l_start(&common);

  //-> triplets to sparse matrix
  {
    const size_t rows = LHS.rows(), cols = LHS.cols(), nnz = LHS.nonZeros();
    cholmod_triplet *T = cholmod_l_allocate_triplet(rows, cols, nnz, 0, CHOLMOD_REAL, &common);
    if ( T == nullptr ) {
      cerr << "[ERR] cholmod_l_allocate_triplet" << endl;
      return EXIT_FAILURE;
    }

    size_t k = 0;
    long* i = (long*) T->i;
    long* j = (long*) T->j;
    double* v = (double*) T->x;
    for (size_t col = 0; col < LHS.outerSize(); ++col) {
      for (SparseMatrix<double>::InnerIterator it(LHS, col); it; ++it) {
        i[k] = it.row();
        j[k] = it.col();
        v[k] = it.value();
        ++k;
      }
    }
    T->nnz = k;
    
    A = cholmod_l_triplet_to_sparse(T, T->nnz, &common);
    if ( A == nullptr ) {
      cerr << "[ERR] cholmod_l_triplet_to_sparse" << endl;
      return EXIT_FAILURE;
    }
    
    cholmod_l_free_triplet(&T, &common);
  }

  //-> copy rhs to b
  b = cholmod_l_zeros(A->nrow, 1, A->xtype, &common);
  double *raw_b = (double *)b->x;
  std::copy(rhs.data(), rhs.data()+rhs.size(), raw_b);
  
  //-> solve
  x = SuiteSparseQR<double>(A, b, &common);

  //-> copy back solution
  double *raw_x = (double*)(x->x);
  xx = VectorXd::Zero(LHS.cols());
  std::copy(raw_x, raw_x+xx.size(), xx.data());

  cholmod_l_free_sparse(&A, &common);
  cholmod_l_free_dense(&x, &common);
  cholmod_l_free_dense(&b, &common);
  cholmod_l_finish(&common);

  return EXIT_SUCCESS;
}

}
