#ifndef FEM_CONSTITUTIVE
#define FEM_CONSTITUTIVE
#include <Eigen/Dense>
#include <cmath>

namespace marvel{
using namespace Eigen;




template<typename T, size_t dim_>
class elas_csttt{
 public:
  // virtual ~elas_csttt(){}
  static T
  val(const Eigen::Matrix<T, dim_, dim_>& F, const double& lam, const double& mu) ;
  
  static  Eigen::Matrix<T, dim_ * dim_, 1>
  gra(const Eigen::Matrix<T, dim_, dim_>& F, const double& lam, const double& mu);
  
  static  Eigen::Matrix<T, dim_ * dim_, dim_ * dim_>
  hes(const Eigen::Matrix<T, dim_, dim_>& F, const double& lam, const double& mu) ;
};

template<typename T, size_t dim_>
class linear_csttt : public elas_csttt<T, dim_>{
 public:
  
  static T
  val(const Eigen::Matrix<T, dim_, dim_>& F, const double& lam, const double& mu) {
    Matrix<T, dim_, dim_> strain = 0.5 * (F + F.transpose()) - Matrix<T, dim_, dim_>::Identity();
    return mu * (strain.array() * strain.array()).sum() + 0.5 * lam * strain.trace() * strain.trace();
  }

  static Matrix<T, dim_ * dim_, 1>
  gra(const Eigen::Matrix<T, dim_, dim_>& F, const double& lam, const double& mu){
    const Matrix<T, dim_, dim_> Iden = Matrix<T, dim_, dim_>::Identity();
    Matrix<T, dim_, dim_> strain = 0.5 * (F + F.transpose()) - Iden;
    Matrix<T, dim_ , dim_> gra_mat = mu * (F + F.transpose() - 2 * Iden) + lam * (F - Iden).trace() * Iden;
    Map<Matrix<T, dim_ * dim_, 1>> gra_vec(gra_mat.data());
    return std::move(gra_vec);
  }

  static Eigen::Matrix<T, dim_ * dim_, dim_ * dim_>
  hes(const Eigen::Matrix<T, dim_, dim_>& F, const double& lam, const double& mu) {
    static Matrix<T, dim_ * dim_, dim_ * dim_> hes;
    static bool have_calc = false;
    if(!have_calc){
      const Matrix<T, dim_ * dim_, dim_ * dim_> Iden = Matrix<T, dim_ * dim_, dim_ * dim_>::Identity();
      Matrix<T, dim_, dim_> strain = 0.5 * (F + F.transpose()) - Matrix<T, dim_, dim_>::Identity();

      Matrix<T, dim_ * dim_, dim_ * dim_> DDtrace = Matrix<T, dim_ * dim_, dim_ * dim_>::Zero();
      for(size_t row = 0; row < dim_ * dim_; row += dim_ +1){
        for(size_t col = 0; col < dim_ * dim_; col += dim_ + 1){
          DDtrace(row, col) = 1;
        }
      }

      Matrix<T, dim_ * dim_, dim_ * dim_> Dsquare = Matrix<T, dim_ * dim_, dim_ * dim_>::Zero();{
        for(size_t row = 0; row < dim_; ++row){
          for(size_t col = 0; col < dim_; ++col){
            Dsquare(row * dim_ + col, col * dim_ + row) = 1;
          }
        }
        Dsquare += Iden;
      }
      hes = mu * Dsquare + lam * DDtrace;
      have_calc = true;
    }
    
    return std::move(hes);
  }
  
  
};

}
#endif
