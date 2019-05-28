#include "basic_energy.h"
#include <Eigen/SparseCore>

using namespace std;
using namespace Eigen;

namespace marvel{
/******************************************momentum*******************************/

template<size_t dim_>
size_t momentum<dim_>::Nx() const{
  return dim_ * dof_;
}

template<size_t dim_>
momentum<dim_>::momentum(const size_t dof, const VectorXd& mass_vec, const double& dt):dof_(dof), dispk_(VectorXd::Zero(dim_ * dof)), vk_(VectorXd::Zero(dim_ *dof)), dt_(dt), d1dt_(1 / dt), d1dtdt_(1 / dt /dt), mass_vec_(dim_ * dof){

  #pragma omp parallel for
  for(size_t i = 0; i < dof; ++i){
    for(size_t j = 0; j < dim_; ++j){
      mass_vec_(i * dim_ + j) = mass_vec(i);
    }
  }

}

template<size_t dim_>
int momentum<dim_>::Val(const double *x, data_ptr<dim_> &data) const{
  Map<const VectorXd> _x(x, dim_ * dof_);
  const VectorXd acce = (_x  - dispk_) * d1dt_ - vk_;
  // data.val_ += 0.5 * acce.dot(mass_sparse_ * acce);
  
  data->save_val(0.5 * acce.dot(mass_vec_.cwiseProduct(acce)));
  return 0;
}
template<size_t dim_>
int momentum<dim_>::Gra(const double *x, data_ptr<dim_> &data) const{
  Map<const VectorXd> _x(x, dim_ * dof_);

  const VectorXd acce = (_x  - dispk_) * d1dtdt_  - vk_ * d1dt_;
  // data.gra_ += mass_sparse_ * acce;
  data->save_gra(mass_vec_.cwiseProduct(acce));

  return 0;
}

template<size_t dim_>
int momentum<dim_>::Hes(const double *x, data_ptr<dim_> &data) const{
  #pragma omp parallel for
  for(size_t i = 0; i < dim_ *dof_; ++i){
    // data.hes_trips.push_back(Triplet<double>(i, i, d1dtdt_ * mass_sparse_.coeff(i, i)));
    data->save_hes(i, i, d1dtdt_ * mass_vec_(i));
  }
  return 0;
}

template<size_t dim_>
int momentum<dim_>::update_location_and_velocity(const double *new_dispk_ptr) {
  Map<const VectorXd> new_dispk(new_dispk_ptr, dim_ * dof_);
  vk_ = (new_dispk - dispk_) * d1dt_;
  dispk_ = new_dispk;
  return 0;
}



/******************************************momentum*******************************/
/******************************************position_constraint*******************************/
template<size_t dim_>
position_constraint<dim_>::position_constraint(const size_t dof, const double &w, const vector<size_t> &cons):w_(w), cons_(cons), dof_(dof){}

template<size_t dim_>
int position_constraint<dim_>::Val(const double *x, data_ptr<dim_> &data) const{
  Map<const MatrixXd> _x(x, dim_, dof_);
  for(auto iter_c = cons_.begin(); iter_c != cons_.end(); ++iter_c)
    data->save_val( w_ * _x.col(*iter_c).dot(_x.col(*iter_c)));
  return 0;
}

template<size_t dim_>
int position_constraint<dim_>::Gra(const double *x, data_ptr<dim_> &data)const {
  Map<const MatrixXd> _x(x, dim_, dof_);
  for(auto iter_c = cons_.begin(); iter_c != cons_.end(); ++iter_c)
    data->save_gra(*iter_c, 2.0 * w_ * _x.col(*iter_c));
  return 0;
}

template<size_t dim_>
int position_constraint<dim_>::Hes(const double *x, data_ptr<dim_> &data) const{
  Map<const MatrixXd> _x(x, dim_, dof_);
  for(auto iter_c = cons_.begin(); iter_c != cons_.end(); ++iter_c){
    for(size_t j = 0; j < dim_; ++j){
      data->save_hes(*iter_c*dim_ + j, *iter_c*dim_ + j, 2 * w_);
    }
  }
  return 0;
}

template<size_t dim_>
size_t position_constraint<dim_>::Nx() const{
  return dim_ * dof_;
}
/******************************************position_constraint*******************************/

/******************************************gravity*******************************/
//dof here is not about dim_
template<size_t dim_>
gravity_energy<dim_>::gravity_energy(const size_t dof, const double &w_g, const double &gravity, const VectorXd &mass, const char &axis):
    w_g_(w_g), dof_(dof), gravity_(gravity), mass_(mass), axis_(axis){}

template<size_t dim_>
int gravity_energy<dim_>::Val(const double *x, data_ptr<dim_> &data) const{

  Map<const MatrixXd> _x(x, dim_, dof_);
  size_t which_axis = size_t(axis_ - 'x');
  data->save_val((_x.row(which_axis).transpose().array() * mass_.array()).sum() * w_g_ * gravity_);
  return 0;
}

template<size_t dim_>
int gravity_energy<dim_>::Gra(const double *x, data_ptr<dim_> &data) const{
  size_t which_axis = size_t(axis_ - 'x');

  MatrixXd g(dim_, dof_);
  //do not add mass!!!1
  // g.row(which_axis) = VectorXd::Constant(dof_, -gravity_ * w_g_).transpose();
  g.row(which_axis) = VectorXd::Constant(dof_, gravity_ * w_g_).cwiseProduct(mass_).transpose();
  Map<const VectorXd> g_(g.data(), dim_ * dof_);
  data->save_gra(g_);
  // Gra.row(which_axis) += g;
  return 0;
}
template<size_t dim_>
int gravity_energy<dim_>::Hes(const double *x, data_ptr<dim_> &data) const{
  return 0;
}

template<size_t dim_>
size_t gravity_energy<dim_>::Nx() const{
  return dim_ * dof_;
}
/******************************************gravity*******************************/

/*************************************collision*********************************/
template<size_t dim_>
collision<dim_>::collision(const size_t dof_, const double &w_coll, const char &ground_axis, const double &ground_pos, const size_t &num_surf_point, const shared_ptr<MatrixXd>& init_points_ptr):ground_axis_(ground_axis), ground_pos_(ground_pos), w_coll_(w_coll), num_surf_point_(num_surf_point), dof_(dof_), init_points_ptr_(init_points_ptr){}

template<size_t dim_>
int collision<dim_>::Val(const double *x, data_ptr<dim_> &data) const{
  const size_t which_axis = size_t(ground_axis_ - 'x');
  
  Map<const MatrixXd> _x(x, dim_, dof_);
  for(size_t i = 0; i < dof_; ++i){
    const double position_now = _x(which_axis, i) + (*init_points_ptr_)(which_axis, i);
    if (( position_now - ground_pos_) < 0){
      data->save_val(w_coll_ * pow((ground_pos_ - position_now), 2));      
    }

  }
  return 0;
}
template<size_t dim_>
int collision<dim_>::Gra(const double *x, data_ptr<dim_> &data) const {
  const size_t which_axis = size_t(ground_axis_ - 'x');
  
  Map<const MatrixXd> _x(x, dim_, dof_);
  for(size_t i = 0; i < dof_; ++i){
    const double position_now = _x(which_axis, i) + (*init_points_ptr_)(which_axis, i);
    if (( position_now - ground_pos_) < 0){
      data->save_gra(i * dim_ + which_axis, 2 * w_coll_ * (position_now - ground_pos_));  
    }

  }
  return 0;

}
template<size_t dim_>
int collision<dim_>::Hes(const double *x, data_ptr<dim_> &data) const{
  const size_t which_axis = size_t(ground_axis_ - 'x');
  Map<const MatrixXd> _x(x, dim_, dof_);
  #pragma omp parallel for
  for(size_t i = 0; i < dof_; ++i){
    const double position_now = _x(which_axis, i) + (*init_points_ptr_)(which_axis, i);
    if (( position_now - ground_pos_) < 0){
      for(size_t j = 0; j < dim_; ++j){
// #pragma omp critical
//         {
//         data.hes_trips.push_back(Triplet<double>(i * 3 + j, i * 3 + j, 2 * w_coll_));
//         }

        data->save_hes(i * dim_ + j, i * dim_ + j, 2 * w_coll_);
      }
    }
    
  }

  return 0;
}
template<size_t dim_>
size_t collision<dim_>::Nx() const{
  return dim_ * dof_;
}

/*************************************collision*********************************/


}//namespace marvel
