#include "data_stream.h"
#include <iostream>
using namespace Eigen;
using namespace std;
namespace marvel{

energy_dat::energy_dat(const size_t dim):dat_str<double>::dat_str(dim), dim_(dim),def_gra_(dim), hes_(3*dim_, 3*dim_){

  zero_mats = vector<Matrix3d>(dim, Matrix3d::Zero());
  inv_A_all_ = zero_mats;
  set_zero();
  sigma_w_points_.setZero(3, dim);


}

// energy_dat::energy_dat(const energy_dat &other){
  
// }

int energy_dat::save_val(const double val){
  dat_str<double>::val_ += val;
  return 0;
}

int energy_dat::save_ele_mat(const size_t &ele_id, const Matrix3d &ele_mat, vector<Matrix3d> &whole_mat, bool if_plus){
  if(if_plus)
    whole_mat[ele_id] = whole_mat[ele_id] + ele_mat;
  else
    whole_mat[ele_id] = ele_mat;
  return 0;
}

int energy_dat::save_ele_vol_cross(const size_t &ele_id, const MatrixXd &ele_mat){
  return save_ele_mat(ele_id, ele_mat, vol_cross_, false);
}

int energy_dat::save_ele_def_gra(const size_t &ele_id, const MatrixXd &ele_mat){
  return save_ele_mat(ele_id, ele_mat, def_gra_, false);
}

int energy_dat::save_ele_inv_all(const size_t &ele_id, const MatrixXd &ele_mat){
  return save_ele_mat(ele_id, ele_mat, inv_A_all_, false);
}


int energy_dat::save_ele_gra(const size_t &ele_id, const Vector3d &ele_mat){
  for(size_t i = 0; i < 3; ++i){
    dat_str<double>::gra_(ele_id * 3 + i) += ele_mat(i);
  }
  
  return 0;
}

int energy_dat::save_ele_strain(const size_t &ele_id, const MatrixXd &ele_mat){
  return save_ele_mat(ele_id, ele_mat, strain_, false);    
}

int energy_dat::save_ele_stress(const size_t &ele_id, const MatrixXd &ele_mat){
  return save_ele_mat(ele_id, ele_mat, stress_, false);    }


int energy_dat::set_zero(){

  vol_cross_ = zero_mats;
  def_gra_ = zero_mats;
  strain_ = zero_mats;
  stress_ = zero_mats;
  

  dat_str<double>::val_ = 0;
  dat_str<double>::gra_.setZero(3 * dim_);
  std::vector<Triplet<double>> to_swap(0);
  dat_str<double>::hes_trips.swap(to_swap);
  
  hes_.setZero();

  return 0;
}


//TODO:: look for a way to return a pointer to the 
// const MatrixXd& energy_dat::ele_mat(const size_t &ele_id, const size_t &rows, const size_t &cols, MatrixXd &mat){
//   return Map<MatrixXd>(mat.col(ele_id).data(), rows, cols);
// }

// MatrixXd& energy_dat::ele_def_gra(const size_t &ele_id){
//   return ele_mat(ele_id, 3, 3, def_gra_);
// }
// MatrixXd& energy_dat::ele_inv_all(const size_t &ele_id){
//   return ele_mat(ele_id, 3, 3, inv_A_all_);  
// }
// MatrixXd& energy_dat::ele_gra(const size_t &ele_id){
//   return ele_mat(ele_id, 3, 1, gra_);    
// }
// MatrixXd& energy_dat::ele_hes(const size_t &ele_id){
//   return ele_mat(ele_id, 3, 3, hes_);      
// }


}
