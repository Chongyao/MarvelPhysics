#ifndef MARVEL_POISSON_FEM
#define MARVEL_POISSON_FEM

#include "def.h"
#include "data_str_core.h"
#include "gaussian_quadrature.h"
#include "basis_func.h"

#include <Eigen/Dense>
#include "eigen_ext.h"
#include <iostream>
#include <set>
#include <fstream>

namespace marvel{
// Functional dimension is 1 since the solution is scalar function.
template<typename T, size_t dim_, size_t num_per_cell_, size_t bas_order_, size_t qdrt_axis_,
         template<typename, size_t, size_t, size_t > class BASIS, //  basis
         template<typename, size_t, size_t, size_t> class QDRT> //
class poisson : public Functional<T, 1>{
 public:
  poisson(const Eigen::Matrix<T, dim_, -1>& nods, const Eigen::Matrix<int, num_per_cell_, -1>& cells, const Eigen::Matrix<T, -1, 1>& k);
  size_t Nx() const;
  int Val(const T *x, std::shared_ptr<dat_str_core<T,dim_>>& data) const;
  int Gra(const T *x, std::shared_ptr<dat_str_core<T,dim_>>& data) const;
  int Hes(const T *x, std::shared_ptr<dat_str_core<T,dim_>>& data) const;
 private:
  void PreComputation();
  const size_t num_nods_, num_cells_, num_qdrt_;
  const Eigen::Matrix<T, dim_, -1> nods_; // vertices
  const Eigen::Matrix<int, num_per_cell_, -1> cells_; // elements
  const Eigen::Matrix<T, -1, 1>& k_;
  const Matrix<int, 1, 1> all_rows_;

 private:
  const qdrt quadrature_;
  std::vector<std::vector<Eigen::Matrix<T, dim_, dim_>>> Dm_inv_;
  std::vector<std::vector<Eigen::Matrix<T, num_per_cell_, dim_>>> Dphi_Dxi_;
  std::vector<std::vector<T>> Jac_det_;
};



#define POI_TEMP template<typename T, size_t dim_, size_t num_per_cell_, size_t bas_order_, size_t qdrt_axis_,template<typename, size_t, size_t, size_t > class BASIS,template<typename, size_t, size_t, size_t> class QDRT>
#define POI_CLASS BaseElas<T, dim_, num_per_cell_, bas_order_, qdrt_axis_, BASIS, QDRT>


}
#endif


