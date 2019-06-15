#ifndef FEM_BASIS
#define FEM_BASIS
#include <Eigen/Dense>
#include <unsupported/Eigen/KroneckerProduct>

namespace marvel{
using namespace Eigen;
using namespace std;



template<typename T, size_t dim_, size_t order_, size_t num_per_cell_>
class basis_func{
 public:

  static void calc_InvDm_Det(const Eigen::Matrix<T, dim_, 1>& PNT, const T* X, double& Jac_det, Matrix<T, dim_, dim_>& Dm_inv);
  static void get_def_gra(const Eigen::Matrix<T, dim_, 1>&PNT, const T* const x, const T* const X, Eigen::Matrix<T, dim_, dim_> & def_gra) ;
  static void get_Ddef_Dx(const Eigen::Matrix<T, dim_, 1>&PNT, const T* const x, const T* const X, const Eigen::Matrix<T, dim_, dim_>& def_gra, Eigen::Matrix<T, dim_ * dim_, dim_ * num_per_cell_>& Ddef_Dx);
  static void get_Ddef_Dx(const T* const X, Eigen::Matrix<T, 9, 12>& Ddef_Dx);
};

//TODO: integrate tri mesh
template<typename T>
class basis_func<T, 3, 1, 4>{
  //TODO:DX_D can be calculated once
 public:

  static void calc_InvDm_Det(const Eigen::Matrix<T, 3, 1>& PNT, const T* X, double& Jac_det, Matrix<T, 3, 3>& Dm_inv) {
    Dm_inv.setZero();

    const Map<const Matrix<T, 3, 4>> rest(X);
    Matrix<T, 3, 3> Dm;
    for(int i = 0; i < 3; i++)
      Dm.col(i) = rest.col(i) - rest.col(3);
    Dm_inv = Dm.inverse();
    Jac_det = Dm.determinant();
  }

  static void get_def_gra(const Eigen::Matrix<T, 3, 1>&PNT, const T* const x, const Matrix<T, 3, 3>& Dm_inv, Eigen::Matrix<T, 3, 3> & def_gra) {
    const Map<const Matrix<T, 3, 4>> deformed(x);
    const Matrix<T, 3, 3> Dx_D = deformed.block(0, 0, 3, 3) - deformed.col(3) * Matrix<T, 1, 3>::Ones();
    def_gra = Dx_D * Dm_inv;
    return;
  }

  //TODO : too many zeros in the matrix can acce
  static void get_Ddef_Dx(const Eigen::Matrix<T, 3, 1>&PNT, const T* const x, const T* const X, const Eigen::Matrix<T, 3, 3>& def_gra, Eigen::Matrix<T, 9, 12>& Ddef_Dx){
    Ddef_Dx.setZero();
    const Map<const Matrix<T, 3, 4>> rest(X);
    const Matrix<T, 3, 3> Drest_D = rest.block(0, 0, 3, 3) - rest.col(3) * Matrix<T, 1, 3>::Ones();
    Matrix<T, 3, 3> inv_Drest_D = Drest_D.inverse().transpose();
    Ddef_Dx.block(0, 0, 9, 9) = kroneckerProduct(inv_Drest_D, Matrix<T, 3, 3>::Identity());
    for(size_t i = 0; i < 3; ++i){
      for(size_t j = 0; j < 3; ++j){
        Ddef_Dx.col(9 + i) -= Ddef_Dx.col(j*3 + i);
      }
    }

    return;
  }

  static void get_Ddef_Dx(const T* const X, Eigen::Matrix<T, 9, 12>& Ddef_Dx){
    Ddef_Dx.setZero();
    const Map<const Matrix<T, 3, 4>> rest(X);
    const Matrix<T, 3, 3> Drest_D = rest.block(0, 0, 3, 3) - rest.col(3) * Matrix<T, 1, 3>::Ones();
    Matrix<T, 3, 3> inv_Drest_D = Drest_D.inverse().transpose();
    Ddef_Dx.block(0, 0, 9, 9) = kroneckerProduct(inv_Drest_D, Matrix<T, 3, 3>::Identity());
    for(size_t i = 0; i < 3; ++i){
      for(size_t j = 0; j < 3; ++j){
        Ddef_Dx.col(9 + i) -= Ddef_Dx.col(j*3 + i);
      }
    }

    return;
  }
  
};

  


}
#endif
