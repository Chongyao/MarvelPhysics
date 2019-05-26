#ifndef FEM_BASIS
#define FEM_BASIS
#include <Eigen/Dense>
namespace marvel{
template<typename T, size_t dim_, size_t order_>
class basis_func{
 public:
  virtual void get_def_gra(const Eigen::MatrixBase<T>&PNT, const T* const x, const T* const X, Eigen::DenseBase<T>& def_gra) const = 0;
  virtual void get_Ddef_Dx(const Eigen::MatrixBase<T>&PNT, const T* const x, const T* const X, const Eigen::DenseBase<T>& def_gra, Eigen::DenseBase<T>& Ddef_Dx) const = 0;
};
}
#endif
