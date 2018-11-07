#include "data_stream.h"

using namespace Eigen;
namespace marvel{
energy_dat::energy_dat(const size_t &dim):dim_(dim),def_gra_(9, dim_), inv_A_all_(9, dim_), gra_(3, dim_), hes_(9, dim_), strain_(9, dim_), stress_(9, dim_), pre_F_(9, dim_){
  gra_.setZero(3, dim_);
  hes_.setZero(3, dim_);
}

// energy_dat::energy_dat(const energy_dat &other){
  
// }

int energy_dat::save_ele_mat(const size_t &ele_id, const size_t &rows, const MatrixXd &ele_mat, MatrixXd &whole_mat, bool if_plus){
  if(if_plus)
    whole_mat.col(ele_id) += Map<const VectorXd>(ele_mat.data(), rows);
  else
    whole_mat.col(ele_id) = Map<const VectorXd>(ele_mat.data(), rows);
  return 0;
}
int energy_dat::save_ele_def_gra(const size_t &ele_id, const MatrixXd &ele_mat){
  return save_ele_mat(ele_id, 9, ele_mat, def_gra_, false);
}
int energy_dat::save_ele_inv_all(const size_t &ele_id, const MatrixXd &ele_mat){
  return save_ele_mat(ele_id, 9, ele_mat, inv_A_all_, false);
  return 0;
}
int energy_dat::save_ele_gra(const size_t &ele_id, const MatrixXd &ele_mat){
  return save_ele_mat(ele_id, 3, ele_mat, gra_, true);    
}
int energy_dat::save_ele_hes(const size_t &ele_id, const MatrixXd &ele_mat){
  return save_ele_mat(ele_id, 3, ele_mat, hes_, true);    
}
int energy_dat::save_ele_strain(const size_t &ele_id, const MatrixXd &ele_mat){
  return save_ele_mat(ele_id, 9, ele_mat, strain_, false);    
}
int energy_dat::save_ele_stress(const size_t &ele_id, const MatrixXd &ele_mat){
  return save_ele_mat(ele_id, 9, ele_mat, stress_, false);    }
int energy_dat::save_ele_pre_F(const size_t &ele_id, const MatrixXd &ele_mat){
  return save_ele_mat(ele_id, 9, ele_mat, pre_F_, false);    }

int energy_dat::set_zero(){
  def_gra_.setZero(9, dim_);
  gra_.setZero(3, dim_);
  hes_.setZero(9, dim_);
  strain_.setZero(9, dim_);
  stress_.setZero(9, dim_);
  pre_F_.setZero(9, dim_);

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
