#include "basis_func.h"
namespace marvel{
using namespace Eigen;
using namespace std;

template<typename T, size_t dim_, size_t order_, size_t num_per_cell_>
static void shape_func<T, dim_, order_, num_per_cell_>::calc_basis_value(const Eigen::Matrix<T, dim_, 1>& PNT, const T* X, Eigen::Matrix<T, 8, 1>& basis_value){
  assert(0);
  return;
}



}
