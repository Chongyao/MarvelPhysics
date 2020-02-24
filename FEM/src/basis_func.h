#ifndef FEM_BASIS
#define FEM_BASIS
#include <Eigen/Dense>
#include <unsupported/Eigen/KroneckerProduct>

namespace marvel{
using namespace Eigen;
using namespace std;

template<typename T, size_t dim_, size_t order_, size_t num_per_cell_>
struct impl_Dphi_Dxi{
  static void calc_Dhpi_Dxi(const Eigen::Matrix<T, dim_, 1>& PNT, const T* X,  Eigen::Matrix<T, num_per_cell_, dim_>& Dphi_Dxi){
    assert(0);
    return;
  }
};

template<typename T>
struct impl_Dphi_Dxi<T, 3, 1, 4>{
  static void calc_Dphi_Dxi(const Eigen::Matrix<T, 3, 1>& PNT, const T* X,  Eigen::Matrix<T, 4, 3>& Dphi_Dxi){
    Dphi_Dxi.setZero();
    Dphi_Dxi.template topRows<3>().setIdentity();
    Dphi_Dxi.row(3) = Eigen::Matrix<T, 1, 3>::Ones() * (-1);
    return;
  }
};

template<typename T>
struct impl_Dphi_Dxi<T, 3, 1, 8>{
  static void calc_Dphi_Dxi(const Eigen::Matrix<T, 3, 1>& PNT, const T* X,  Eigen::Matrix<T, 8, 3>& Dphi_Dxi){
    Dphi_Dxi.setZero();
    const double xi0 = PNT(0), xi1 = PNT(1), xi2 = PNT(2);

    Map<Matrix<T, 24, 1>> jac(Dphi_Dxi.data());
    const T
        tt1 = 1 - xi1, tt2 = 1 - xi2,
        tt3 = 1 + xi1, tt4 = 1 + xi2,
        tt5 = 1 - xi0, tt6 = 1 + xi0;

    jac[0] = -(tt1*tt2)/8.0;
    jac[1] = (tt1*tt2)/8.0;
    jac[2] = -(tt3*tt2)/8.0;
    jac[3] = (tt3*tt2)/8.0;
    jac[4] = -(tt1*tt4)/8.0;
    jac[5] = (tt1*tt4)/8.0;
    jac[6] = -(tt3*tt4)/8.0;
    jac[7] = (tt3*tt4)/8.0;
    jac[8] = -(tt5*tt2)/8.0;
    jac[9] = -(tt6*tt2)/8.0;
    jac[10] = (tt5*tt2)/8.0;
    jac[11] = (tt6*tt2)/8.0;
    jac[12] = -(tt5*tt4)/8.0;
    jac[13] = -(tt6*tt4)/8.0;
    jac[14] = (tt5*tt4)/8.0;
    jac[15] = (tt6*tt4)/8.0;
    jac[16] = -(tt5*tt1)/8.0;
    jac[17] = -(tt6*tt1)/8.0;
    jac[18] = -(tt5*tt3)/8.0;
    jac[19] = -(tt6*tt3)/8.0;
    jac[20] = (tt5*tt1)/8.0;
    jac[21] = (tt6*tt1)/8.0;
    jac[22] = (tt5*tt3)/8.0;
    jac[23] = (tt6*tt3)/8.0;

    Dphi_Dxi /= 8.0;
    return;

  }
};


  



template<typename T, size_t dim_, size_t order_, size_t num_per_cell_>
class basis_func{
 public:
  static void calc_Dphi_Dxi(const Eigen::Matrix<T, dim_, 1>& PNT, const T* X,  Eigen::Matrix<T, num_per_cell_, dim_>& Dphi_Dxi){
    return impl_Dphi_Dxi<T, dim_, order_, num_per_cell_>::calc_Dphi_Dxi(PNT, X, Dphi_Dxi);
  }
 
  static void calc_InvDm_Det(const Eigen::Matrix<T, num_per_cell_, dim_>& Dphi_Dxi, const T* X, T& Jac_det, Matrix<T, dim_, dim_>& Dm_inv){
    Dm_inv.setZero();
    const Eigen::Map<const Matrix<T, dim_, num_per_cell_>> rest(X);
    Matrix<T, dim_, dim_> Dm = rest * Dphi_Dxi;
    Dm_inv = Dm.inverse();
    Jac_det = fabs(Dm.determinant());
    return;
  }
  static void get_def_gra(const Eigen::Matrix<T, num_per_cell_, dim_>& Dphi_Dxi, const T* const x, const Matrix<T, dim_, dim_>& Dm_inv,  Eigen::Matrix<T, dim_, dim_> & def_gra) {
    const Eigen::Map<const Matrix<T, dim_, num_per_cell_>> deformed(x);
    def_gra = deformed * Dphi_Dxi * Dm_inv;
    return;
  }
  
  static void get_Ddef_Dx(const Eigen::Matrix<T, num_per_cell_, dim_>& Dphi_Dxi, const Eigen::Matrix<T, dim_, dim_>& Dm_inv, Eigen::Matrix<T, dim_ * dim_, dim_ * num_per_cell_>& Ddef_Dx){
    Ddef_Dx.setZero();
    const Eigen::Matrix<T, num_per_cell_, dim_> Ddef_Dx_compressed = Dphi_Dxi * Dm_inv;
    #pragma omp parallel for
    for(size_t i = 0; i < num_per_cell_; ++i)
      for(size_t j = 0; j < dim_; ++j)
        Ddef_Dx.block(j * dim_, i * dim_, dim_, dim_) = Eigen::Matrix<T, dim_, dim_>::Identity() * Ddef_Dx_compressed(i, j);
    return;
  }
  
  // static void get_Ddef_Dx(const T* const X, Eigen::Matrix<T, 9, 12>& Ddef_Dx);
};

//TODO: integrate tri mesh
// template<typename T>
// static void basis_func<T, 3, 1, 4>::calc_Dphi_Dxi(const Eigen::Matrix<T, 3, 1>& PNT, const T* X,  Eigen::Matrix<T, 4, 3>& Dphi_Dxi){
//     Dphi_Dxi.setZero();
//     Dphi_Dxi.template topRows<3>().setIdentity();
//     Dphi_Dxi.row(3) = Eigen::Matrix<T, 1, 3>::Ones() * (-1);
//     return;
// }

  //TODO : too many zeros in the matrix can acce
  // static void get_Ddef_Dx(const Eigen::Matrix<T, 3, 1>&PNT, const T* const x, const T* const X, const Eigen::Matrix<T, 3, 3>& def_gra, Eigen::Matrix<T, 9, 12>& Ddef_Dx){
  //   Ddef_Dx.setZero();
  //   const Eigen::Map<const Matrix<T, 3, 4>> rest(X);
  //   const Matrix<T, 3, 3> Drest_D = rest.block(0, 0, 3, 3) - rest.col(3) * Matrix<T, 1, 3>::Ones();
  //   Matrix<T, 3, 3> inv_Drest_D = Drest_D.inverse().transpose();
  //   Ddef_Dx.block(0, 0, 9, 9) = kroneckerProduct(inv_Drest_D, Matrix<T, 3, 3>::Identity());
  //   for(size_t i = 0; i < 3; ++i){
  //     for(size_t j = 0; j < 3; ++j){
  //       Ddef_Dx.col(9 + i) -= Ddef_Dx.col(j*3 + i);
  //     }
  //   }

  //   return;
  // }

  // static void get_Ddef_Dx(const T* const X, Eigen::Matrix<T, 9, 12>& Ddef_Dx){
  //   Ddef_Dx.setZero();
  //   const Eigen::Map<const Matrix<T, 3, 4>> rest(X);
  //   const Matrix<T, 3, 3> Drest_D = rest.block(0, 0, 3, 3) - rest.col(3) * Matrix<T, 1, 3>::Ones();
  //   Matrix<T, 3, 3> inv_Drest_D = Drest_D.inverse().transpose();
  //   Ddef_Dx.block(0, 0, 9, 9) = kroneckerProduct(inv_Drest_D, Matrix<T, 3, 3>::Identity());
  //   for(size_t i = 0; i < 3; ++i){
  //     for(size_t j = 0; j < 3; ++j){
  //       Ddef_Dx.col(9 + i) -= Ddef_Dx.col(j*3 + i);
  //     }
  //   }

  //   return;
  // }
  
// };

// template<typename T>
// class basis_func<T, 3, 1, 8>{
//   //The order of the  vertices of one hex is the same as that in vtk format.
//  public:
//     static void calc_Dphi_Dxi(const Eigen::Matrix<T, 3, 1>& PNT, const T* X,  Eigen::Matrix<T, 4, 3>& Dphi_Dxi){
//     Dphi_Dxi.setZero();
//     const double xi0 = PNT(0), xi1 = PNT(1), xi2 = PNT(2);
//     Dphi_Dxi <<
//         (xi1 - 1)*(1 - xi2), (xi0 - 1)*(1 - xi2), (xi0 - 1)*(1 - xi1),
//         (1 - xi1)*(1 - xi2), (xi0 - 1)*(1 - xi2), (xi0 - 1)*(1 - xi1),
//         (xi1 - 1)*(1 - xi2), (1 - xi0)*(1 - xi2), (xi0 - 1)*(1 - xi1),
//         (1 - xi1)*(1 - xi2), (1 - xi0)*(1 - xi2), (xi0 - 1)*(1 - xi1),
//         (xi1 - 1)*(1 - xi2), (xi0 - 1)*(1 - xi2), (1 - xi0)*(1 - xi1),
//         (1 - xi1)*(1 - xi2), (xi0 - 1)*(1 - xi2), (1 - xi0)*(1 - xi1),
//         (xi1 - 1)*(1 - xi2), (1 - xi0)*(1 - xi2), (1 - xi0)*(1 - xi1),
//         (1 - xi1)*(1 - xi2), (1 - xi0)*(1 - xi2), (1 - xi0)*(1 - xi1);
//     return;
//   }

//   //TODO : too many zeros in the matrix can acce
//   static void get_Ddef_Dx(const Eigen::Matrix<T, 3, 1>&PNT, const T* const x, const T* const X, const Eigen::Matrix<T, 3, 3>& def_gra, Eigen::Matrix<T, 9, 12>& Ddef_Dx){
//     Ddef_Dx.setZero();
//     const Eigen::Map<const Matrix<T, 3, 4>> rest(X);
//     const Matrix<T, 3, 3> Drest_D = rest.block(0, 0, 3, 3) - rest.col(3) * Matrix<T, 1, 3>::Ones();
//     Matrix<T, 3, 3> inv_Drest_D = Drest_D.inverse().transpose();
//     Ddef_Dx.block(0, 0, 9, 9) = kroneckerProduct(inv_Drest_D, Matrix<T, 3, 3>::Identity());iuuy8
//     for(size_t i = 0; i < 3; ++i){
//       for(size_t j = 0; j < 3; ++j){
//         Ddef_Dx.col(9 + i) -= Ddef_Dx.col(j*3 + i);
//       }
//     }

//     return;
//   }
  
  
// static void get_Ddef_Dx(const T* const X, Eigen::Matrix<T, 9, 12>& Ddef_Dx){
//     Ddef_Dx.setZero();
//     const Eigen::Map<const Matrix<T, 3, 4>> rest(X);
//     const Matrix<T, 3, 3> Drest_D = rest.block(0, 0, 3, 3) - rest.col(3) * Matrix<T, 1, 3>::Ones();
//     Matrix<T, 3, 3> inv_Drest_D = Drest_D.inverse().transpose();
//     Ddef_Dx.block(0, 0, 9, 9) = kroneckerProduct(inv_Drest_D, Matrix<T, 3, 3>::Identity());
//     for(size_t i = 0; i < 3; ++i){
//       for(size_t j = 0; j < 3; ++j){
//         Ddef_Dx.col(9 + i) -= Ddef_Dx.col(j*3 + i);
//       }
//     }

//     return;
//   }
  
// };

  


}
#endif
