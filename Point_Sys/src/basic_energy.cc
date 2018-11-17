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
  cout << dat_str.gra_ <<endl;
  for(auto iter_c = cons_.begin(); iter_c != cons_.end(); ++iter_c){
    dat_str.save_ele_gra(*iter_c, -2.0 * w_ * _disp.col(*iter_c));
    // cout << "teset norm" << endl;
    // cout << _disp.col(*iter_c).norm()<<endl;
    // dat_str.save_ele_gra(*iter_c,  -w_ * _disp.col(*iter_c) / _disp.col(*iter_c).norm());
  }
  cout << "[INFO]>>>>>>>>>>>>>>>>>>>after positon<<<<<<<<<<<<<<<<<<" << endl;
  cout << dat_str.gra_ << endl;
  return 0;
}

int position_constraint::Hes(const double *disp, energy_dat &dat_str){
  Map<const MatrixXd> _disp(disp, 3, dim_);
  for(auto iter_c = cons_.begin(); iter_c != cons_.end(); ++iter_c){
    for(size_t j = 0; j < 3; ++j){
      dat_str.hes_trips.push_back(Triplet<double>(*iter_c*3 + j, *iter_c*3 + j, -2 * w_));
      cout << "p = " << *iter_c*3 + j << " q = " << *iter_c*3 + j << 2*w_ << endl;;
    }
  }
  // for(auto iter_c = cons_.begin(); iter_c != cons_.end(); ++iter_c){
  //   for(size_t m = 0; m < 3; ++m){
  //     double norm = _disp.col(*iter_c).norm();
  //     double kmm = w_ * (norm*norm - _disp(m, *iter_c)*_disp(m, *iter_c)) /pow(norm, 3);
  //     dat_str.hes_trips.push_back(Triplet<double>(*iter_c*3 + m, *iter_c*3 + m, kmm));
      
  //     double km_m1 = -w_ *_disp(m, *iter_c) * _disp((m + 1)%3, *iter_c)/pow(norm, 3);
  //     dat_str.hes_trips.push_back(Triplet<double>(*iter_c*3 + m, *iter_c*3 + (m + 1)%3, km_m1));
      
  //     double km_m2 = -w_ *_disp(m, *iter_c) * _disp((m + 2)%3, *iter_c)/pow(norm, 3);
  //     dat_str.hes_trips.push_back(Triplet<double>(*iter_c*3 + m, *iter_c*3 + (m + 2)%3, km_m2));
  //   }
  // }
  return 0;
}
/******************************************position_constraint*******************************/

/******************************************gravity*******************************/
gravity_energy::gravity_energy(const double &w_g, const double &gravity, const size_t &dim, const VectorXd &mass, const char &axis):
    w_g_(w_g), dim_(dim), gravity_(gravity), mass_(mass), axis_(axis){}

int gravity_energy::Val(const double *disp, energy_dat &dat_str){
  Map<const MatrixXd> _disp(disp, 3, dim_);
  size_t which_axis = size_t(axis_ - 'x');
  cout << "gravity energy " << (_disp.row(which_axis).transpose().array() * mass_.array()).sum() * w_g_ * gravity_ << endl;
  dat_str.Val_ += (_disp.row(which_axis).transpose().array() * mass_.array()).sum() * w_g_ * gravity_;
  return 0;
}
int gravity_energy::Gra(const double *disp, energy_dat &dat_str){
  size_t which_axis = size_t(axis_ - 'x');

  MatrixXd g(3, dim_);
  g.setZero(3, dim_);
  g.row(which_axis) = VectorXd::Constant(dim_, -gravity_ * w_g_).cwiseProduct(mass_).transpose();
  Map<MatrixXd> Gra(dat_str.gra_.data(), 3, dim_);
  Gra += g;
  return 0;
}
int gravity_energy::Hes(const double *disp, energy_dat &dat_str){
  return 0;
}
}//namespace marvel
