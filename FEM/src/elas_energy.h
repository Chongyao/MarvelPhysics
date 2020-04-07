#ifndef ELAS_ENERGY
#define ELAS_ENERGY
#include "def.h"
#include "data_str_core.h"
#include "constitutive.h"
#include "gaussian_quadrature.h"
#include "basis_func.h"

#include <Eigen/Dense>
#include "eigen_ext.h"
#include <iostream>
#include <set>
#include <fstream>
namespace marvel{

using namespace Eigen;
using namespace std;

template<typename T>
inline void compute_lame_coeffs(const T Ym, const T Pr,
                                T &mu, T &lambda) {
  mu = Ym/(2*(1+Pr));
  lambda = Ym*Pr/((1+Pr)*(1-2*Pr));
}
template<typename T, size_t dim_, size_t num_per_cell_, size_t bas_order_, size_t qdrt_axis_,
         template<typename, size_t> class CSTTT,  // constituitive function
         template<typename, size_t, size_t, size_t > class BASIS, //  basis
         template<typename, size_t, size_t, size_t> class QDRT> //
class BaseElas : public Functional<T, dim_>{
  using basis = BASIS<T, dim_, bas_order_, num_per_cell_>;
  using csttt = CSTTT<T, dim_>;
  using qdrt = QDRT<T, dim_, qdrt_axis_, num_per_cell_>;

 public:
  BaseElas(const Eigen::Matrix<T, dim_, -1>& nods, const Eigen::Matrix<int, num_per_cell_, -1>& cells, const T& ym, const T&poi);
  size_t Nx() const;
  
  int Val(const T *x, std::shared_ptr<dat_str_core<T,dim_>>& data) const;
  int Gra(const T *x, std::shared_ptr<dat_str_core<T,dim_>>& data) const;
  int Hes(const T *x, std::shared_ptr<dat_str_core<T,dim_>>& data) const;
  
  protected:
  void PreComputation();
  
 private:
  const size_t all_dim_, num_nods_, num_cells_, num_qdrt_;
  Eigen::Matrix<T, 2, -1> mtr_;
  const Eigen::Matrix<T, dim_, -1> nods_; // vertices
  const Eigen::Matrix<int, num_per_cell_, -1> cells_; // elements
  Matrix<int, dim_, 1> all_rows_;
  const qdrt quadrature_;

  
 private:  // precomputed values
  std::vector<std::vector<Eigen::Matrix<T, dim_, dim_>>> Dm_inv_;
  std::vector<std::vector<T>> Jac_det_;
  std::vector<std::vector<Eigen::Matrix<T, dim_ * dim_, dim_ * num_per_cell_>, Eigen::aligned_allocator<Eigen::Matrix<T, dim_ *dim_, dim_ * num_per_cell_>>>> Ddef_Dx_;
  std::vector<std::vector<Eigen::Matrix<T, num_per_cell_, dim_>>> Dphi_Dxi_;


};

#define ELAS_TEMP template<typename T, size_t dim_, size_t num_per_cell_, size_t bas_order_, size_t qdrt_axis_,template<typename, size_t> class CSTTT,template<typename, size_t, size_t, size_t > class BASIS,template<typename, size_t, size_t, size_t> class QDRT>
#define ELAS_CLASS BaseElas<T, dim_, num_per_cell_, bas_order_, qdrt_axis_, CSTTT, BASIS, QDRT>



}
#endif
