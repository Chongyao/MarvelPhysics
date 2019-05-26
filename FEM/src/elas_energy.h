#ifndef ELAS_ENERGY
#define ELAS_ENERGY
#include "def.h"
#include "data_str_core.h"

#include "constitutive.h"
#include "gaussian_quadrature.h"
#include "basis_func.h"

#include <Eigen/Dense>
namespace marvel{
template<typename T>
inline void compute_lame_coeffs(const T Ym, const T Pr,
                                T &mu, T &lambda) {
  mu = Ym/(2*(1+Pr));
  lambda = Ym*Pr/((1+Pr)*(1-2*Pr));
}


template<typename T, size_t dim_, size_t num_per_cell_, size_t bas_order_, size_t num_quad_,
         template<typename, size_t> class CSTTT,
         template<typename, size_t, size_t> class BASIS,
         template<typename, size_t, size_t> class GAUS>
class BaseElas : public Functional<T, dim_>{
 public:
  BaseElas(const Eigen::DenseBase<T>& nods, const Eigen::DenseBase<T>& cells,
           const double& ym, const double&poi):
      all_dim_(nods.size()), num_nods_(nods.cols()), nods_(nods), cells_(cells){
    
    static_assert(std::is_base_of<elas_csttt<T, dim_>, CSTTT<T, dim_>>::value, "CSTTT must derive from elas_csttt");
    static_assert(std::is_base_of<basis_func<T, dim_, bas_order_>, BASIS<T, dim_, bas_order_>>::value, "BASIS must derive from basis_func");
    static_assert(std::is_base_of<gaus_quad<T, dim_, num_quad_>, GAUS<T, dim_, num_quad_>>::value, "GAUS must derive from gaus_quad");
    
    T mu, lambda;
    compute_lame_coeffs(ym, poi, mu, lambda);
    mtr_.resize(cells_.cols(), 2);
    mtr_.col(0) = Eigen::Matrix<T, -1, 1>::Ones(cells_.cols()) * lambda;
    mtr_.col(1) = Eigen::Matrix<T, -1, 1>::Ones(cells_.cols()) * mu;
  }
  
  size_t Nx() const {return all_dim_;}
    
  int Val(const T *x, std::shared_ptr<dat_str_core<T,dim_>>& data) const {
    // Eigen::Map<const Eigen::Matrix<T, 
  }
  int Gra(const T *x, std::shared_ptr<dat_str_core<T,dim_>>& data) const ;
  int Hes(const T *x, std::shared_ptr<dat_str_core<T,dim_>>& data) const ;
  
 private:
  const size_t all_dim_;
  const size_t num_nods_;
  const Eigen::Matrix<T, -1, 2> mtr_;
  const Eigen::Matrix<T, dim_, -1> nods_;
  const Eigen::Matrix<size_t, num_per_cell_, -1> cells_;
  
};
}
#endif
