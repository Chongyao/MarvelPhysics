#include "basic_energy.h"
#include <Eigen/SparseCore>
#include <iostream>
using namespace std;
using namespace Eigen;

namespace marvel{
/******************************************position_constraint*******************************/
position_constraint::position_constraint(const double &w, const vector<size_t> &cons, const size_t dim):w_(w), cons_(cons), dim_(dim){}
int position_constraint::Gra(const double *disp, energy_dat &dat_str){
  Map<const MatrixXd> _disp(disp, 3, dim_);
  for(auto iter_c = cons_.begin(); iter_c != cons_.end(); ++iter_c)
    dat_str.save_ele_gra(*iter_c, -2.0 * w_ * _disp.col(*iter_c));
  return 0;
}

int position_constraint::Hes(const double *disp, energy_dat &dat_str){
  Map<const MatrixXd> _disp(disp, 3, dim_);
  for(auto iter_c = cons_.begin(); iter_c != cons_.end(); ++iter_c){
    for(size_t j = 0; j < 3; ++j){
      dat_str.hes_trips.push_back(Triplet<double>(*iter_c*3 + j, *iter_c*3 + j, -2 * w_));
    }
  }
  return 0;
}
/******************************************position_constraint*******************************/

/******************************************gravity*******************************/
gravity_energy::gravity_energy(const double &w_g, const double &gravity, const size_t &dim, const VectorXd &mass, const char &axis):
    w_g_(w_g), dim_(dim), gravity_(gravity), mass_(mass), axis_(axis){}

int gravity_energy::Val(const double *disp, energy_dat &dat_str){
  Map<const MatrixXd> _disp(disp, 3, dim_);
  size_t which_axis = size_t(axis_ - 'x');
  dat_str.Val_ += (_disp.row(which_axis).transpose().array() * mass_.array()).sum() * w_g_ * gravity_;
  return 0;
}
int gravity_energy::Gra(const double *disp, energy_dat &dat_str){
  size_t which_axis = size_t(axis_ - 'x');

  MatrixXd g(3, dim_);
  g.setZero(3, dim_);
  //do not add mass!!!1
  // g.row(which_axis) = VectorXd::Constant(dim_, -gravity_ * w_g_).transpose();
  g.row(which_axis) = VectorXd::Constant(dim_, -gravity_ * w_g_).cwiseProduct(mass_).transpose();
  Map<MatrixXd> Gra(dat_str.gra_.data(), 3, dim_);
  Gra += g;
  return 0;
}
int gravity_energy::Hes(const double *disp, energy_dat &dat_str){
  return 0;
}
/******************************************gravity*******************************/

/*************************************collision*********************************/
collision::collision(const double &w_coll, const char &ground_axis, const double &ground_pos, const size_t &num_surf_point, const size_t &dim):ground_axis_(ground_axis), ground_pos_(ground_pos), w_coll_(w_coll), num_surf_point_(num_surf_point), dim_(dim){}

int collision::Val(const double *init_points, const double *disp, energy_dat &dat_str){
  size_t which_axis = size_t(ground_axis_ - 'x');
  
  Map<const MatrixXd> _disp(disp, 3, dim_);
  Map<const MatrixXd> _init_points(init_points, 3, dim_);
  cout << num_surf_point_  <<" " << _disp.cols() << endl;
  for(size_t i = 0; i < dim_; ++i){
    double position_now = _disp(which_axis, i) + _init_points(which_axis, i);
    if (( position_now - ground_pos_) < 0){
      // cout <<" i is " << i << " coll val " << endl;
      // cout << w_coll_ * pow((ground_pos_ - position_now), 2) << endl;
      dat_str.Val_ += w_coll_ * pow((ground_pos_ - position_now), 2);      
    }

  }
  return 0;
}
int collision::Gra(const double *init_points, const double *disp, energy_dat &dat_str, const VectorXd& mass){
  size_t which_axis = size_t(ground_axis_ - 'x');
  
  Map<const MatrixXd> _disp(disp, 3, dim_);
  Map<const MatrixXd> _init_points(init_points, 3, dim_);

  for(size_t i = 0; i < dim_; ++i){
    double position_now = _disp(which_axis, i) + _init_points(which_axis, i);

    if (( position_now - ground_pos_) < 0){
      
      dat_str.gra_(which_axis, i) += 2 * w_coll_ * (ground_pos_ - position_now);  
    }

  }
  return 0;

}
int Hes(const double *disp, energy_dat &dat_str){
  return 0;
}


/*************************************collision*********************************/


}//namespace marvel
