#include "points_energy.h"
#include "get_nn.h"
#include <cmath>
#include <Eigen/SVD>
#include <Eigen/LU>
#include <Eigen/Geometry>
#include <iostream>
#define PI 3.14159265359

using namespace std;
using namespace Eigen;
namespace marvel{



Matrix3d safe_inv(const MatrixXd& sys_mat){
  //TODO: use better solver considering the sysmertic
  JacobiSVD<MatrixXd> svd(sys_mat, ComputeThinU | ComputeThinV);

  auto sin_val = svd.singularValues();
  Matrix3d inv_sin_val;
  inv_sin_val.setZero(3, 3);
  
  for(size_t p = 0; p < 3; ++p){
    inv_sin_val(p, p) = sin_val(p)>0?1/sin_val(p):0;
  }
  
  return std::move(svd.matrixV() * inv_sin_val * svd.matrixU().transpose());
}

void cons_law(const Matrix3d &strain, Matrix3d &stress, const Matrix3d &def_gra, const double &You, const double &Poi){
  //TODO:add more consititutive law

  double G = You/(2 + 2*Poi);
  double lam = You*Poi/(1+Poi)/(1-2*Poi);
  //St.Venant-Kirchhof model
  // double trace = strain(0, 0) + strain(1, 1) + strain(2, 2);
  // stress = def_gra*(2*G*strain + lam*trace*MatrixXd::Identity(3, 3));

  //Linear model
  double trace = def_gra(0, 0) + def_gra(1, 1) + def_gra(2, 2) - 3;
  stress = G*(def_gra + def_gra.transpose() - 2*Matrix3d::Identity()) + lam*trace*Matrix3d::Identity();

}





point_sys::point_sys(const MatrixXd  &points, const double &rho, const double &Young, const double &Poission, const double &vol_all, const double &kv, const vector<vector<size_t>> &friends, const VectorXd &sup_radi):
    points_(points), rho_(rho), Young_(Young), Poission_(Poission), vol_all_(vol_all), dim_(points_.cols()), kv_(kv), friends_(friends), sup_radi_(sup_radi){

  // sup_radi_ = SH_.get_sup_radi();
  
  mass_i_.setZero(dim_);
  //init
  double mass_total = rho*vol_all;
  double mass_sigma = (sup_radi_/3).array().cube().sum();
  scal_fac_ = mass_total/mass_sigma;
  mass_i_ = (scal_fac_*rho_*sup_radi_).array().cube();
}

double point_sys::get_mass(const size_t &i) const{
  return mass_i_(i);
}

size_t point_sys::Nx() const{
  return dim_;
}
int point_sys::calc_weig() const{
  // friends_ = vector<vector<size_t>>(dim_);
  weig_ = vector<vector<double>>(dim_);
#pragma parallel omp for  
  for(size_t i = 0; i < dim_; ++i){
    vector<double> weig_of_one_p;
    for(auto one_fri : friends_[i]){
      weig_of_one_p.push_back(kernel(i, one_fri));
    }
    weig_[i] = weig_of_one_p;
  }
}

int point_sys::calc_rhoi_vi() const{
  //init
  rho_i_.setZero(dim_);
  vol_i_.setZero(dim_);

  calc_weig();
  assert(friends_.size() > 0 && weig_.size() > 0);
#pragma parallel omp for
  for(size_t i = 0; i < dim_; ++i){
    for(size_t j = 0; j < friends_[i].size(); ++j){
      rho_i_(i) += mass_i_(i)*weig_[i][j];
    }
  }

  vol_i_ = mass_i_.array() / rho_i_.array();
  return 0;
}
double point_sys::kernel(const double &r, const double &h) const {
  if(r < h)
    return 315*pow((h*h - r*r), 3)/(64*PI*pow(h, 9));
  else
    return 0;
}

double point_sys::kernel(const size_t &i, const size_t &j) const{
  double r = (points_.col(j) - points_.col(i)).norm();
  return kernel(r, sup_radi_(i));
}

int point_sys::calc_defo_gra(const double *disp, energy_dat &dat_str) const{
  Map<const Matrix<double, Dynamic, Dynamic> > _disp(disp, 3, dim_);
#pragma parallel omp for
  for(size_t i = 0; i < dim_; ++i){
    Matrix3d one_du;
    Matrix3d b = Matrix3d::Zero();
    for(size_t j = 0; j < friends_[i].size(); ++j){
      Vector3d xij = points_.col(friends_[i][j]) - points_.col(i);
      for(size_t k = 0; k < 3; ++k){
        b.col(k) += (_disp(k, friends_[i][j]) - _disp(k, i))*weig_[i][j]*xij;
        
      }
    }
    
    //clac inverse of sys_mat by SVD
    
    Map<Matrix3d> inv_A(dat_str.inv_A_all_.col(i).data());
    
    for(size_t k = 0; k < 3; ++k){
      one_du.row(k) = (inv_A * b.col(k)).transpose();
    }
    MatrixXd F = one_du + MatrixXd::Identity(3, 3);
    dat_str.save_ele_def_gra(i, F);

  }
  return 0;
}



int point_sys::pre_compute(energy_dat &dat_str) const {
  // Map<const Matrix<double, Dynamic, Dynamic> > points_curr(x, 3, dim_);
  //elasiticity do not updata
  // SH_.update_points(points_curr);
  //
  calc_rhoi_vi();
 #pragma parallel omp for
  for(size_t i = 0; i < dim_; ++i){
    Matrix3d sys_mat;
    sys_mat.setZero(3, 3);
    
    // assert(friends_[i].size() >= 3);
    for(size_t j = 0; j < friends_[i].size(); ++j){
      Vector3d xij = points_.col(friends_[i][j]) - points_.col(i);
      sys_mat += weig_[i][j]*xij*xij.transpose();
      
    }
    
    auto inv_A = safe_inv(sys_mat);
    
    dat_str.save_ele_inv_all(i, inv_A);
  }
  return 0;
}
int point_sys::Gra(const double *disp, energy_dat &dat_str) const{
#pragma parallel omp for
  for(size_t i = 0; i < dim_; ++i){
    
    
    Map<MatrixXd> def_gra(dat_str.def_gra_.col(i).data(), 3, 3);
    
    //calculate strain and stress
    Matrix3d strain = def_gra.transpose()*def_gra - Matrix3d::Identity();
    // Matrix3d strain = 0.5*(def_gra + def_gra.transpose()) - Matrix3d::Identity();
    
    dat_str.save_ele_strain(i, strain);
    Matrix3d stress;
    cons_law(strain, stress, def_gra, Young_, Poission_);
    
    dat_str.save_ele_stress(i, stress);
    //calculate Fv
    Matrix3d gra_def_gra;
    auto trans_def_gra = def_gra.transpose();
    for(int j = 0; j < 3; ++j){
      Vector3d cross1 = trans_def_gra.col((j+1)%3), cross2 = trans_def_gra.col((j+2)%3);
      gra_def_gra.col(j) = cross1.cross(cross2);
    }
    gra_def_gra.transposeInPlace();

    //assemble Fe and Fv
    Map<MatrixXd> inv_A(dat_str.inv_A_all_.col(i).data(), 3, 3);
    Matrix3d pre_F = -vol_i_(i)*(2*def_gra*stress + kv_*(def_gra.determinant() - 1)*gra_def_gra)*inv_A;
    dat_str.save_ele_pre_F(i, pre_F);
    
    //add to gra_
    Vector3d di;
    di.setZero(3);
    for(size_t j = 0; j < friends_[i].size(); ++j){
      double w = weig_[i][j];
      Vector3d xij = (points_.col(friends_[i][j]) - points_.col(i));
      di += -w*xij;
      
      dat_str.save_ele_gra(friends_[i][j], w*pre_F*xij);
    }
    

    dat_str.save_ele_gra(i, pre_F*di);
  }

  return 0;
}


int point_sys::gravity(const double *x, energy_dat &dat_str,  const double &gravity) const{

  MatrixXd g(3, dim_);
  g.setZero(3, dim_);
  g.row(2) = VectorXd::Constant(dim_, -gravity).cwiseProduct(mass_i_).transpose();
  Map<MatrixXd> Gra(dat_str.gra_.data(), 3, dim_);
  Gra += g;
  return 0;
}


}//namespace : marvel


