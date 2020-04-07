#ifndef FEM_CONSTITUTIVE
#define FEM_CONSTITUTIVE
#include <Eigen/Dense>
#include <cmath>
#include <Eigen/Geometry>
#include <unsupported/Eigen/CXX11/Tensor>
#include "tensor.h"
#include <iostream>
namespace marvel{
using namespace Eigen;
using namespace std;



template<typename T, size_t dim_, size_t field_>
class constitutive{
  static T
  val(const Eigen::Matrix<T, field_, dim_>& F, const Eigen::Matrix<T, -1, 1>& mtr) ;
  static  Eigen::Matrix<T, field_ * dim_, 1>
  gra(const Eigen::Matrix<T, field_, dim_>& F, const Eigen::Matrix<T, -1, -1>& mtr);
  
  static  Eigen::Matrix<T, field_ * dim_, field_ * dim_>
  hes(const Eigen::Matrix<T, field_, dim_>& F,const Eigen::Matrix<T, -1, -1>& mtr);  
};



template<typename T, size_t dim_, size_t field_>
class elas_csttt : public constitutive<T, dim_, field_>{
};

template<typename T, size_t dim_, size_t field_>
class arap_csttt : public constitutive<T, dim_, field_>{

 public:
  
  static T
  val(const Eigen::Matrix<T, dim_, dim_>& F, const Eigen::Matrix<T, -1, -1>& mtr) {
    T lam = mtr(0), mu = mtr(1);
    Matrix<T, dim_, dim_> R;
    {
      Eigen::Matrix3d A; A=F.template cast<double>();
        Eigen::Quaterniond q=Eigen::Quaterniond(A);
        q.normalize();
        for(size_t iter=0;iter<5;iter++)
          {      
            Eigen::Matrix3d R=q.matrix();
            Eigen::Vector3d omega=(R.col(0).cross(A.col(0))+R.col(1).cross(A.col(1))+R.col(2).cross(A.col(2)))*(1.0/fabs(R.col(0).dot(A.col(0))+R.col(1).dot(A.col(1))+R.col(2).dot(A.col(2)))+1e-9);
            double w=omega.norm();
            if(w<1e-9)
        {
          break; 
        }
            q=Eigen::Quaterniond(AngleAxisd(w,(1.0/w)*omega))*q;
            q.normalize();
          }
        R = q.matrix().cast<T>();
    }
    return lam/2*(F-R).squaredNorm();
  }

  static Matrix<T, dim_ * dim_, 1>
  gra(const Eigen::Matrix<T, dim_, dim_>& F, const Eigen::Matrix<T, -1, 1>& mtr){
    T lam = mtr(0), mu = mtr(1);
    Matrix<T, dim_, dim_> R;
    {
      Matrix3d A; A=F.template cast<double>();
        Eigen::Quaterniond q=Eigen::Quaterniond(A);
        q.normalize();
        for(size_t iter=0;iter<5;iter++)
          {      
            Eigen::Matrix3d R=q.matrix();
            Eigen::Vector3d omega=(R.col(0).cross(A.col(0))+R.col(1).cross(A.col(1))+R.col(2).cross(A.col(2)))*(1.0/fabs(R.col(0).dot(A.col(0))+R.col(1).dot(A.col(1))+R.col(2).dot(A.col(2)))+1e-9);
            double w=omega.norm();
            if(w<1e-9)
        {
          break; 
        }
            q=Eigen::Quaterniond(AngleAxisd(w,(1.0/w)*omega))*q;
            q.normalize();
          }
        R = q.matrix().cast<T>();
    }
    Matrix<T, dim_*dim_, 1> gra_vec;
    Eigen::Map<Matrix<T, dim_, dim_>>(gra_vec.data()) = lam*(F-R);
    return std::move(gra_vec);
  }

  static Eigen::Matrix<T, dim_ * dim_, dim_ * dim_>
  hes(const Eigen::Matrix<T, dim_, dim_>& F, const Eigen::Matrix<T, -1, 1>& mtr) {
    T lam = mtr(0), mu = mtr(1);
    Matrix<T, dim_*dim_, dim_*dim_> hes;
    hes.setIdentity();
    hes *= lam;
    return std::move(hes);
  }
  
};


template<typename T, size_t dim_, size_t field_>
class linear_csttt : public elas_csttt<T, dim_, field_>{
 public:
  
  static T
  val(const Eigen::Matrix<T, dim_, dim_>& F, const Eigen::Matrix<T, -1, 1>& mtr) {
    const T lam = mtr(0), mu = mtr(1);
    Matrix<T, dim_, dim_> strain = 0.5 * (F + F.transpose()) - Matrix<T, dim_, dim_>::Identity();
    return mu * (strain.array() * strain.array()).sum() + 0.5 * lam * strain.trace() * strain.trace();
  }

  static Matrix<T, dim_ * dim_, 1>
  gra(const Eigen::Matrix<T, dim_, dim_>& F, const Eigen::Matrix<T, -1, 1>& mtr){
    const T lam = mtr(0), mu = mtr(1);
    const Matrix<T, dim_, dim_> Iden = Matrix<T, dim_, dim_>::Identity();
    Matrix<T, dim_, dim_> strain = 0.5 * (F + F.transpose()) - Iden;
    Matrix<T, dim_ , dim_> gra_mat = mu * (F + F.transpose() - 2 * Iden) + lam * (F - Iden).trace() * Iden;
    Eigen::Map<Matrix<T, dim_ * dim_, 1>> gra_vec(gra_mat.data());
    return std::move(gra_vec);
  }

  static Eigen::Matrix<T, dim_ * dim_, dim_ * dim_>
  hes(const Eigen::Matrix<T, dim_, dim_>& F, const Eigen::Matrix<T, -1, 1>& mtr) {
    const T lam = mtr(0), mu = mtr(1);
    static Matrix<T, dim_ * dim_, dim_ * dim_> hes;
    static bool have_calc = false;
    if(!have_calc){
      const Matrix<T, dim_ * dim_, dim_ * dim_> Iden = Matrix<T, dim_ * dim_, dim_ * dim_>::Identity();
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

template<typename T, size_t dim_, size_t field_>
class stvk : public elas_csttt<T, dim_, field_>{
 public:
  using tensor_type = fourth_tensor<T, dim_, dim_, dim_, dim_>;
  static T
  val(const Eigen::Matrix<T, dim_, dim_>& F, const Eigen::Matrix<T, -1, 1>& mtr) {
    const T lam = mtr(0), mu = mtr(1);
    Matrix<T, dim_, dim_> strain = 0.5 * (F.transpose() * F - Matrix<T, dim_, dim_>::Identity());
    return mu * (strain.array() * strain.array()).sum() + 0.5 * lam * strain.trace() * strain.trace();
  }

  static Matrix<T, dim_ * dim_, 1>
  gra(const Eigen::Matrix<T, dim_, dim_>& F, const Eigen::Matrix<T, -1, 1>& mtr){
    T lam = mtr(0), mu = mtr(1);
    const Matrix<T, dim_, dim_> Iden = Matrix<T, dim_, dim_>::Identity();
    Matrix<T, dim_, dim_> strain = 0.5 * (F.transpose() * F - Iden);
    Matrix<T, dim_ , dim_> gra_mat = F * (2 * mu * strain + lam * strain.trace() * Iden);
    Eigen::Map<Matrix<T, dim_ * dim_, 1>> gra_vec(gra_mat.data());
    return std::move(gra_vec);
  }

  static Eigen::Matrix<T, dim_ * dim_, dim_ * dim_>
  hes(const Eigen::Matrix<T, dim_, dim_>& F, const Eigen::Matrix<T, -1, 1>& mtr) {
    T lam = mtr(0), mu = mtr(1);
    Matrix<T, dim_, dim_> strain = 0.5 * (F.transpose() * F - Matrix<T, dim_, dim_>::Identity());
    const Matrix<T, dim_, dim_> Iden = Matrix<T, dim_, dim_>::Identity();
    static Matrix<T, dim_ * dim_, dim_ * dim_> hes;
    //TODO: fill this hes
     tensor_type dF_dF;{
      for(size_t row_out = 0; row_out < dim_; ++row_out){
        for(size_t col_out = 0; col_out < dim_; ++col_out){
          Matrix<T, dim_, dim_> zero = Matrix<T, dim_, dim_>::Zero();
          zero(row_out, col_out) = 1;
          dF_dF(row_out, col_out) = zero;
        }
      }
    }


    Matrix<T, dim_, dim_> rhs = 2 * mu * strain + lam * strain.trace() * Iden;


    tensor_type drhs_dF;{
      decltype(dF_dF) dstrain_dF;{
        for(size_t row_out = 0; row_out < dim_; ++row_out){
          for(size_t col_out = 0; col_out < dim_; ++col_out){
            Matrix<T, dim_, dim_> zero = Matrix<T, dim_, dim_>::Zero();
            zero.col(row_out) += F.col(col_out);
            zero.col(col_out) += F.col(row_out);
            dstrain_dF(row_out, col_out) = mu * zero;
            
          }
        }
      }

      tensor_type  dtrace_dF;{
        for(size_t row_out = 0; row_out < dim_; ++row_out){
          dtrace_dF(row_out, row_out) = lam * F;
        }
      }
      drhs_dF = dstrain_dF + dtrace_dF;
    }


    tensor_type hes_tensor = dF_dF * rhs + F * drhs_dF ;
    hes_tensor.Flatten(hes);

    return std::move(hes);
  }
  
  
};









}
#endif




