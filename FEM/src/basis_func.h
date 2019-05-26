#ifndef FEM_BASIS
#define FEM_BASIS
#include <Eigen/Dense>
namespace marvel{
using namespace Eigen;

template<typename T, size_t dim_, size_t order_>
class basis_func{
 public:
  virtual void get_def_gra(const Eigen::MatrixBase<T>&PNT, const T* const x, const T* const X, Eigen::DenseBase<T>& def_gra) const = 0;
  virtual void get_Ddef_Dx(const Eigen::MatrixBase<T>&PNT, const T* const x, const T* const X, const Eigen::DenseBase<T>& def_gra, Eigen::DenseBase<T>& Ddef_Dx) const = 0;
};

template<typename T>
class tet_linear_basis : public basis_func<T, 3, 1>{
  void get_def_gra(const Eigen::MatrixBase<T>&PNT, const T* const x, const T* const X, Eigen::DenseBase<T>& def_gra) const{
    const Map<const Matrix<T, 3, 4>> deformed(x);
    const Map<const Matrix<T, 3, 4>> rest(X);
    const Matrix<T, 3, 3> Dx_D = deformed.block(0, 0, 3, 3) - deformed.col(3) * Matrix<T, 1, 3>::Ones();
    const Matrix<T, 3, 3> DX_D = rest.block(0, 0, 3, 3) - rest.col(3) * Matrix<T, 1, 3>::Ones();
    def_gra = Dx_D * DX_D.inverse();
    return;
  }

  //TODO : too many zeros in the matrix can acce
  void get_Ddef_Dx(const Eigen::MatrixBase<T>&PNT, const T* const x, const T* const X, const Eigen::DenseBase<T>& def_gra, Eigen::DenseBase<T>& Ddef_Dx) const{
    Ddef_Dx.setZero();
    const Map<const Matrix<T, 3, 4>> rest(X);
    const Matrix<T, 3, 3> Drest_D = rest.block(0, 0, 3, 3) - rest.col(3) * Matrix<T, 1, 3>::Ones();
    const Matrix<T, 3, 4> inv_Drest_D; {
      inv_Drest_D.block(0, 0, 3, 3) = Drest_D.inverse();
      inv_Drest_D.col(3) = - inv_Drest_D.col(0) - inv_Drest_D.col(1) - inv_Drest_D.col(2);
    }
    
    for(size_t d = 0; d < 3; ++d){
      Ddef_Dx.block(d * 3, d * 4, 3, 4) = inv_Drest_D;
    }
    return;
  }
  
};

  


}
#endif
