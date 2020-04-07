#ifndef ELAS_ENERGY
#define ELAS_ENERGY
#include "FEM.h"
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

template<typename T>
inline void compute_lame_coeffs(const T Ym, const T Pr,
                                T &mu, T &lambda) {
  mu = Ym/(2*(1+Pr));
  lambda = Ym*Pr/((1+Pr)*(1-2*Pr));
}

template<typename T, size_t dim_, size_t num_per_cell_, size_t bas_order_, size_t qdrt_axis_,
         template<typename, size_t, size_t> class CSTTT,  // constituitive function
         template<typename, size_t, size_t, size_t, size_t > class BASIS, //  basis
         template<typename, size_t, size_t, size_t> class QDRT> //
class BaseElas : public finite_element<T, dim_, dim_, num_per_cell_, bas_order_, qdrt_axis_, CSTTT, BASIS, QDRT>{
 public:
  using base_class = finite_element<T, dim_, dim_, num_per_cell_, bas_order_, qdrt_axis_, CSTTT, BASIS, QDRT>;
  BaseElas(const Eigen::Matrix<T, dim_, -1>& nods, const Eigen::Matrix<int, num_per_cell_, -1>& cells, const T& ym, const T&poi);

};

#define ELAS_CLASS BaseElas<T, dim_, num_per_cell_, bas_order_, qdrt_axis_, CSTTT, BASIS, QDRT>
#define ELAS_TEMP template<typename T, size_t dim_, size_t num_per_cell_, size_t bas_order_, size_t qdrt_axis_,template<typename, size_t, size_t> class CSTTT,template<typename, size_t, size_t, size_t, size_t > class BASIS, template<typename, size_t, size_t, size_t> class QDRT> //


}
#endif
