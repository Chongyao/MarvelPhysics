#ifndef DATA_STREAM_H
#define DATA_STREAM_H
#include <Eigen/Core>
#include <Eigen/SparseCore>
#include "data_str_core.h"
#include <memory>
namespace marvel{
using VM3 = std::vector<Eigen::Matrix3d>;
using SVM = std::shared_ptr<std::vector<Eigen::Matrix3d>>;


class energy_dat: public dat_str_core<double, 3>{
 public:

  energy_dat(const size_t dof);
  // ~energy_dat();



  
  // energy_dat(const energy_dzat &other);
  //save data by element
  int save_ele_vol_cross(const size_t &ele_id, const Eigen::MatrixXd &ele_mat);
  int save_ele_def_gra(const size_t &ele_id, const Eigen::MatrixXd &ele_mat);
  int save_ele_inv_all(const size_t &ele_id, const Eigen::MatrixXd &ele_mat);
  // int save_ele_gra(const size_t &ele_id, const Eigen::Vector3d &ele_mat);
  // int save_ele_hes(const size_t &ele_id, const Eigen::MatrixXd &ele_mat);
  int save_ele_strain(const size_t &ele_id, const Eigen::MatrixXd &ele_mat);
  int save_ele_stress(const size_t &ele_id, const Eigen::MatrixXd &ele_mat);
  int save_ele_pre_F(const size_t &ele_id, const Eigen::MatrixXd &ele_mat);
  // int save_val(const double val);
  int set_zero();
  //get data by element

 // private:
  
  VM3 vol_cross_;
  VM3 def_gra_;
  Eigen::Matrix<double, 3, Eigen::Dynamic> sigma_w_points_;
  VM3 inv_A_all_;
  VM3 strain_;
  VM3 stress_;


  Eigen::SparseMatrix<double> hes_;
  int save_ele_mat(const size_t &ele_id, const Eigen::Matrix3d &ele_mat, std::vector<Eigen::Matrix3d> &whole_mat, bool if_plus);
 private:
  VM3 zero_mats;

};

}
#endif
