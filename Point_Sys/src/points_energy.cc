#include "points_energy.h"
#include "get_nn.h"
#include <cmath>
#include <Eigen/LU>
#define PI 3.14159265359

using namespace std;
using namespace Eigen;
namespace marvel{

point_sys::point_sys(const MatrixXd  &points, const double &rho, const double &vol_all, const size_t &nearest_num):
    points_(points), rho_(rho), vol_all_(vol_all), dim_(points_.cols()), nearest_num_(nearest_num),SH_(points, nearest_num){

  sup_radi_ = SH_.get_sup_radi();
  mass_i_.setZero(dim_);
  //init
  double mass_total = rho*vol_all;
  double mass_sigma = (sup_radi_/3).array().cube().sum();
  scal_fac_ = mass_total/mass_sigma;
  mass_i_ = (scal_fac_*rho_*sup_radi_).array().cube();
  // for(size_t i = 0; i < dim_; ++i){
  //   mass_i(i) = scal_fac_ * pow(sup_radi(i) , 3)*rho_;
  // }



}

size_t point_sys::Nx() const{
  return dim_;
}
int point_sys::calc_fri(){
  friends_.clear();
  weig_.clear();
  for(size_t i = 0; i < dim_; ++i){
    vector<size_t> one_fris;
    vector<double> weig_of_one_p;
    SH_.get_friends(i, sup_radi_(i), one_fris);
    friends_.push_back(one_fris);
    for(auto one_fri : one_fris){
      weig_of_one_p.push_back(kernel(i, one_fri));
    }
  }
}

int point_sys::calc_rhoi_vi(const double *x){
  //init
  Map<const MatrixXd, 0, OuterStride<> > points_curr(x, 3, dim_, OuterStride<>(dim_));
  rho_i_.setZero(dim_);
  vol_i_.setZero(dim_);

  calc_fri();
  assert(friends_.size() > 0 && weig_.size() > 0);
  for(size_t i = 0; i < dim_; ++i){
    for(size_t j = 0; j < friends_[i].size(); ++j){
      rho_i_(i) += mass_i_(i)*weig_[i][j];
    }
  }

  vol_i_ = mass_i_.array() / rho_i_.array();

  return 0;
}
double point_sys::kernel(const double &r, const double &h){
  if(r < h)
    return 315*pow((h*h - r*r), 3)/(64*PI*pow(h, 9));
  else
    return 0;
}
double point_sys::kernel(const Eigen::Vector3d &xj, const Eigen::Vector3d &xi, const double &h){
  double r = (xj - xi).norm();
  return kernel(r, h);
}

double point_sys::kernel(const size_t &i, const size_t &j){
  double r = (points_.col(j) - points_.col(i)).norm();
  return kernel(r, sup_radi_(i));
}


int point_sys::calc_defo_gra(const double *x, double *def_gra){
  Map<const MatrixXd, 0, OuterStride<> > points_curr(x, 3, dim_, OuterStride<>(dim_));
  Map<MatrixXd, 0, OuterStride<> > d_u(def_gra, 9, dim_, OuterStride<>(dim_));

  MatrixXd disp = points_curr - points_;
  
  for(size_t i = 0; i < dim_; ++i){
    VectorXd one_du(9);
    Matrix3d sys_mat;
    VectorXd b(9);
    sys_mat.setZero(3, 3);
    b.setZero(9);
    
    for(size_t j = 0; j < friends_[i].size(); ++j){
      Vector3d xij = points_curr.col(friends_[i][j]) - points_curr.col(i);
      sys_mat += weig_[i][j]*xij*xij.transpose();
      for(size_t k = 0; k < 3; ++k){
        b.segment(3*k, 3) += (disp(k, friends_[i][j]) - disp(k, i))*weig_[i][j]*xij;        
      }
    }

    for(size_t k = 0; k < 3; ++k){
      //TODO: use better solver considering the sysmertic
      one_du.segment(3*k, 3) = sys_mat.lu() .solve(b.segment(3*k, 3));      
    }
    d_u.col(i) = one_du;    
  }
  return 0;
}



}//namespace : marvel


