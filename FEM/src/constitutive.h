#ifndef FEM_CONSTITUTIVE
#define FEM_CONSTITUTIVE
#include <Eigen/Dense>
namespace marvel{
template<typename T, size_t dim_>
class elas_csttt{
 public:
  virtual  Eigen::Matrix<T, dim_, dim_>& gra(const Eigen::DenseBase<T>& F)const = 0;
  virtual  Eigen::Matrix<T, dim_ * dim_, dim_ * dim_>& hes(const Eigen::DenseBase<T>& F) const = 0;
};

}
#endif
