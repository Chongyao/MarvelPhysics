#ifndef MASS_MATRIX_H
#define MASS_MATRIX_H

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "eigen_ext.h"
#include <iostream>
namespace marvel {

// using mati_t=zjucad::matrix::matrix<size_t>;
// using matd_t=zjucad::matrix::matrix<double>;
// using spmat_t=Eigen::SparseMatrix<double>;

// int calc_mass_matrix(const mati_t &cell,
//                      const matd_t &nods,
//                      const double rho,
//                      const size_t dim,
//                      spmat_t *M,
//                      bool lumped);


// int calc_surf_mass_matrix(const mati_t &cell, const matd_t &nods,
//                           const double rho, spmat_t *M);

//TODO: integrate mass with baiss and quadrature
template<typename T, size_t dim_, size_t num_per_cell_>
int calc_mass_vector(const Matrix<T, dim_, -1> nods, const Matrix<size_t, num_per_cell_, -1> cells, const T& rho, Matrix<T, -1, 1>& mass_vector){
  std::cout << "use partial specilization" << std::endl;
}

template<typename T>
int calc_mass_vector(const Matrix<T, 3, -1> nods, const Matrix<int, 4, -1> cells, const T& rho, Matrix<T, -1, 1>& mass_vector){
  const size_t num_nods = nods.cols();
  const size_t dim  = nods.rows();
  const Matrix<int, 3, 1> all_rows = Matrix<int, 3, 1>::LinSpaced(dim, 0, dim - 1);

  mass_vector.resize(num_nods);
  mass_vector.setZero();
  #pragma omp parallel for
  for(size_t cell_id = 0; cell_id< cells.cols(); ++cell_id){
    
    Matrix<T, 3, 4> one_tet_ = indexing(nods, all_rows, cells.col(cell_id));
    Matrix<T, 3, 3> one_tet = one_tet_.block(0, 0, 3, 3) - nods.col(cells(3, cell_id)) * Matrix<T, 1, 3>::Ones();
    T volume = fabs(one_tet.determinant()) / 6.0;
    T coeff = rho * volume / 20.0 * 3;
    for(size_t p_id = 0; p_id < cells.rows(); ++p_id){
      #pragma omp atomic
      mass_vector(cells(p_id, cell_id)) += coeff;
    }
  }
}


}

#endif
