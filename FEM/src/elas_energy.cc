#include "elas_energy.h"
#include <Eigen/Dense>
#include "eigen_ext.h"
#include <iostream>
#include <set>
#include <fstream>

using namespace std;
using namespace Eigen;
namespace marvel{

ELAS_TEMP
ELAS_CLASS::BaseElas(const Matrix<T, dim_, -1>& nods, const Matrix<int, num_per_cell_, -1>& cells, const T& ym, const T&poi):
    finite_element<T, dim_, dim_, num_per_cell_, bas_order_, qdrt_axis_, CSTTT, BASIS, QDRT>(nods, cells){
  T mu, lambda;
  compute_lame_coeffs(ym, poi, mu, lambda);
  base_class::mtr_.resize(2, cells.cols());
  base_class::mtr_.row(0) = Matrix<T, 1, -1>::Ones(cells.cols()) * lambda;
  base_class::mtr_.row(1) = Matrix<T, 1, -1>::Ones(cells.cols()) * mu;
}


template class BaseElas<double, 3, 4, 1, 1, linear_csttt, basis_func, quadrature>;
template class BaseElas<double, 3, 8, 1, 2, linear_csttt, basis_func, quadrature>;

template class BaseElas<double, 3, 4, 1, 1, stvk, basis_func, quadrature>;
template class BaseElas<double, 3, 8, 1, 2, stvk, basis_func, quadrature>;


}
