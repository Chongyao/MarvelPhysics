#ifndef MASS_MATRIX_H
#define MASS_MATRIX_H

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "eigen_ext.h"
#include <iostream>
#include <vector>
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
int calc_mass_vector(const Matrix<T, 3, -1> nods, const MatrixXi cells, const T& rho, Matrix<T, -1, 1>& mass_vector){
  const size_t num_nods = nods.cols();
  const size_t dim  = nods.rows();
  const Matrix<int, 3, 1> all_rows = Matrix<int, 3, 1>::LinSpaced(dim, 0, dim - 1);

  mass_vector.resize(num_nods);
  mass_vector.setZero();
  std::vector<Eigen::Triplet<T> > trips;

  // #pragma omp parallel for
  for(size_t cell_id = 0; cell_id< cells.cols(); ++cell_id){
    MatrixXd one_cell;
    T volume;
    T coeff;

    Matrix<T, 3, 4> one_tet_ = indexing(nods, all_rows, cells.col(cell_id));
    one_cell.resize(3, 3);
    one_cell = one_tet_.block(0, 0, 3, 3) - one_tet_.col(3) * Matrix<T, 1, 3>::Ones();
    volume = fabs(one_cell.determinant()) / 6.0;
    coeff = rho * volume / 20.0;


 // for(size_t p_id = 0; p_id < cells.rows(); ++p_id){
    //   #pragma omp atomic
    //   mass_vector(cells(p_id, cell_id)) += coeff;
    // }

    for (size_t p = 0; p < cells.rows(); ++p) {
      for (size_t q = p; q < cells.rows(); ++q) {
        trips.push_back(Triplet<T>(cells(p, cell_id), cells(q, cell_id), coeff));
        trips.push_back(Triplet<T>(cells(q, cell_id), cells(p, cell_id), coeff));
      }
    }
  }

#pragma omp parallel for
  for (size_t i = 0; i < trips.size(); ++i){
    trips[i] = Triplet<T>(trips[i].row(), trips[i].row(), trips[i].value());
  }

  SparseMatrix<T> mass(num_nods, num_nods);
  mass.setFromTriplets(trips.begin(), trips.end());
  #pragma omp parallel for
  for (size_t i = 0; i < num_nods; ++i){
    mass_vector[i] = mass.coeff(i, i);
  }
  return 0;
  
}


template<typename T, size_t dim_, size_t num_per_cell_, size_t bas_order_, size_t num_qdrt_, template<typename, size_t, size_t, size_t > class BASIS, //  basis
         template<typename, size_t, size_t, size_t> class QDRT> //
int mass_calculator(const Eigen::Matrix<T, dim_, -1>& nods, const Eigen::Matrix<int, num_per_cell_, -1>& cells, const T& rho, Matrix<T, -1, 1>& mass_vector){
  using basis = BASIS<T, dim_, bas_order_, num_per_cell_>;
  using qdrt = QDRT<T, dim_, num_qdrt_, num_per_cell_>;
  
  const size_t num_cells = cells.cols(), num_nods = nods.cols();
  const Matrix<int, dim_, 1> all_rows_ = Matrix<int, dim_, 1>::LinSpaced(dim_, 0, dim_ -1);
  
  const qdrt quadrature_ = qdrt();


  mass_vector.resize(num_nods);
  mass_vector.setZero();
  std::vector<Eigen::Triplet<T> > trips;

  #pragma omp parallel for
  for(size_t cell_id = 0; cell_id < num_cells; ++cell_id){
    const Matrix<T, dim_, num_per_cell_> X_cell = indexing(nods, all_rows_, cells.col(cell_id));
    T mass = 0.0;
    for(size_t qdrt_id = 0; qdrt_id < num_qdrt_; ++qdrt_id){
      Eigen::Matrix<T, num_per_cell_, dim_> Dphi_Dxi_tmp;
      Eigen::Matrix<T, dim_, dim_> Dm_inv_tmp;
      T jac_det_tmp;

      basis::calc_Dphi_Dxi(quadrature_.PNT_.col(qdrt_id), X_cell.data(), Dphi_Dxi_tmp);
      basis::calc_InvDm_Det(Dphi_Dxi_tmp, X_cell.data(), jac_det_tmp, Dm_inv_tmp);
      mass += quadrature_.WGT_[qdrt_id] * jac_det_tmp;

    }
    mass *= rho / num_per_cell_;
    for (size_t p = 0; p < cells.rows(); ++p) 
      for (size_t q = p; q < cells.rows(); ++q) {
        #pragma omp critical
        {
          trips.push_back(Triplet<T>(cells(p, cell_id), cells(q, cell_id), mass));
          trips.push_back(Triplet<T>(cells(q, cell_id), cells(p, cell_id), mass));
        }
      }
  }

  #pragma omp parallel for
  for (size_t i = 0; i < trips.size(); ++i){
    trips[i] = Triplet<T>(trips[i].row(), trips[i].row(), trips[i].value());
  }
  SparseMatrix<T> mass(num_nods, num_nods);
  mass.setFromTriplets(trips.begin(), trips.end());
 #pragma omp parallel for
  for (size_t i = 0; i < num_nods; ++i){
    mass_vector[i] = mass.coeff(i, i);
  }
  return 0;
}




}

#endif
