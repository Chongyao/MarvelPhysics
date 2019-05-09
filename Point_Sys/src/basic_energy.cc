#include "basic_energy.h"
#include <Eigen/SparseCore>

using namespace std;
using namespace Eigen;

namespace marvel{
/******************************************momentum*******************************/
momentum::momentum(const size_t dim, const Eigen::SparseMatrix<double>& mass_sparse, const double& dt):mass_sparse_(mass_sparse), dim_(dim), dispk_(VectorXd::Zero(3 * dim)), vk_(VectorXd::Zero(3 *dim)), dt_(dt), d1dt_(1 / dt), d1dtdt_(1 / dt /dt){}

int momentum::Val(const double *disp, energy_dat &dat_str) const{
  Map<const VectorXd> _disp(disp, 3 * dim_);
  const VectorXd acce = (_disp  - dispk_) * d1dt_ - vk_;
  // dat_str.val_ += 0.5 * acce.dot(mass_sparse_ * acce);
  dat_str.save_val(0.5 * acce.dot(mass_sparse_ * acce));
  return 0;
}
int momentum::Gra(const double *disp, energy_dat &dat_str) const{
  Map<const VectorXd> _disp(disp, 3 * dim_);

  const VectorXd acce = (_disp  - dispk_) * d1dtdt_  - vk_ * d1dt_;
  // dat_str.gra_ += mass_sparse_ * acce;
  dat_str.save_gra(mass_sparse_ * acce);

  return 0;
}
int momentum::Hes(const double *disp, energy_dat &dat_str) const{
  for(size_t i = 0; i < 3 *dim_; ++i){
    dat_str.hes_trips.push_back(Triplet<double>(i, i, d1dtdt_ * mass_sparse_.coeff(i, i)));
    // dat_str->hes_trips.push_back(Triplet<double>(i, i, d1dtdt_ * mass_sparse_.coeff(i, i)));
  }
  return 0;
}
int momentum::update_location_and_velocity(const double *new_dispk_ptr) {
  Map<const VectorXd> new_dispk(new_dispk_ptr, 3 * dim_);
  vk_ = (new_dispk - dispk_) * d1dt_;
  dispk_ = new_dispk;
  return 0;
}



/******************************************momentum*******************************/
/******************************************position_constraint*******************************/
position_constraint::position_constraint(const size_t dim, const double &w, const vector<size_t> &cons):w_(w), cons_(cons), dim_(dim){}
int position_constraint::Val(const double *disp, energy_dat &dat_str) const{
  Map<const MatrixXd> _disp(disp, 3, dim_);
  for(auto iter_c = cons_.begin(); iter_c != cons_.end(); ++iter_c)
    dat_str.val_ += w_ * _disp.col(*iter_c).dot(_disp.col(*iter_c));
  return 0;
}
int position_constraint::Gra(const double *disp, energy_dat &dat_str)const {
  Map<const MatrixXd> _disp(disp, 3, dim_);
  for(auto iter_c = cons_.begin(); iter_c != cons_.end(); ++iter_c)
    dat_str.save_ele_gra(*iter_c, 2.0 * w_ * _disp.col(*iter_c));
  return 0;
}
int position_constraint::Hes(const double *disp, energy_dat &dat_str) const{
  Map<const MatrixXd> _disp(disp, 3, dim_);
  for(auto iter_c = cons_.begin(); iter_c != cons_.end(); ++iter_c){
    for(size_t j = 0; j < 3; ++j){
      dat_str.hes_trips.push_back(Triplet<double>(*iter_c*3 + j, *iter_c*3 + j, 2 * w_));
    }
  }
  return 0;
}
/******************************************position_constraint*******************************/

/******************************************gravity*******************************/
gravity_energy::gravity_energy(const size_t dim, const double &w_g, const double &gravity, const VectorXd &mass, const char &axis):
    w_g_(w_g), dim_(dim), gravity_(gravity), mass_(mass), axis_(axis){}
int gravity_energy::Val(const double *disp, energy_dat &dat_str) const{

  Map<const MatrixXd> _disp(disp, 3, dim_);
  size_t which_axis = size_t(axis_ - 'x');
  dat_str.val_ += (_disp.row(which_axis).transpose().array() * mass_.array()).sum() * w_g_ * gravity_;
  return 0;
}

int gravity_energy::Gra(const double *disp, energy_dat &dat_str) const{
  size_t which_axis = size_t(axis_ - 'x');

  MatrixXd g(1, dim_);
  //do not add mass!!!1
  // g.row(which_axis) = VectorXd::Constant(dim_, -gravity_ * w_g_).transpose();
  g = VectorXd::Constant(dim_, gravity_ * w_g_).cwiseProduct(mass_).transpose();
  Map<MatrixXd> Gra(dat_str.gra_.data(), 3, dim_);
  Gra.row(which_axis) += g;
  return 0;
}

int gravity_energy::Hes(const double *disp, energy_dat &dat_str) const{
  return 0;
}
/******************************************gravity*******************************/

/*************************************collision*********************************/
collision::collision(const size_t dim, const double &w_coll, const char &ground_axis, const double &ground_pos, const size_t &num_surf_point, const shared_ptr<MatrixXd>& init_points_ptr):ground_axis_(ground_axis), ground_pos_(ground_pos), w_coll_(w_coll), num_surf_point_(num_surf_point), dim_(dim), init_points_ptr_(init_points_ptr){}

int collision::Val(const double *disp, energy_dat &dat_str) const{
  const size_t which_axis = size_t(ground_axis_ - 'x');
  
  Map<const MatrixXd> _disp(disp, 3, dim_);
  for(size_t i = 0; i < dim_; ++i){
    const double position_now = _disp(which_axis, i) + (*init_points_ptr_)(which_axis, i);
    if (( position_now - ground_pos_) < 0){
      dat_str.val_ += w_coll_ * pow((ground_pos_ - position_now), 2);      
    }

  }
  return 0;
}

int collision::Gra(const double *disp, energy_dat &dat_str) const {
  const size_t which_axis = size_t(ground_axis_ - 'x');
  
  Map<const MatrixXd> _disp(disp, 3, dim_);
  for(size_t i = 0; i < dim_; ++i){
    const double position_now = _disp(which_axis, i) + (*init_points_ptr_)(which_axis, i);
    if (( position_now - ground_pos_) < 0){
      dat_str.gra_(which_axis, i) += 2 * w_coll_ * (position_now - ground_pos_);  
    }

  }
  return 0;

}
int collision::Hes(const double *disp, energy_dat &dat_str) const{
  const size_t which_axis = size_t(ground_axis_ - 'x');
  Map<const MatrixXd> _disp(disp, 3, dim_);
  #pragma omp parallel for
  for(size_t i = 0; i < dim_; ++i){
    const double position_now = _disp(which_axis, i) + (*init_points_ptr_)(which_axis, i);
    if (( position_now - ground_pos_) < 0){
      for(size_t j = 0; j < 3; ++j){
#pragma omp critical
        {
        dat_str.hes_trips.push_back(Triplet<double>(i * 3 + j, i * 3 + j, 2 * w_coll_));
        }
      }
    }
    
  }

  return 0;
}


/*************************************collision*********************************/


}//namespace marvel
