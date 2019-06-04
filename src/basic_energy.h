#ifndef BASIC_ENERGY_H
#define BASIC_ENERGY_H
#include "def.h"
#include "data_str_core.h"
#include <vector>
#include <memory>

namespace marvel{

template<size_t dim_>
using data_ptr = std::shared_ptr<dat_str_core<double, dim_>>;

template<size_t dim_>
class position_constraint : public Functional<double, dim_>{
 public:
  //used for displace based
  position_constraint(const size_t dof, const double &w, const std::vector<size_t> &cons);

  //used for position based
  position_constraint(const double *rest, const size_t dof, const double &w, const std::vector<size_t> &cons);
  int Val(const double *x, data_ptr<dim_> &data) const;  
  int Gra(const double *x, data_ptr<dim_> &data) const ;
  int Hes(const double *x, data_ptr<dim_> &data) const ;
  size_t Nx() const;

 private:
  Eigen::MatrixXd rest_;
  const size_t dof_;
  const double w_;
  const std::vector<size_t> cons_;
};

template<size_t dim_>
class gravity_energy : public Functional<double, dim_>{
 public:
  gravity_energy(const size_t dof, const double &w_g, const double &gravity, const Eigen::VectorXd &mass, const char &axis);
  int Val(const double *disp, data_ptr<dim_> &data) const ;
  int Gra(const double *disp, data_ptr<dim_> &data) const;
  int Hes(const double *disp, data_ptr<dim_> &data) const;
  size_t Nx() const;
 private:
  const char axis_;
  const size_t dof_;
  const double w_g_;
  const double gravity_;
  const Eigen::VectorXd mass_;
};

//simple collision with ground
template<size_t dim_>
class collision : public Functional<double, dim_>{
 public:
  collision(const size_t dim, const double &w_coll, const char &ground_axis, const double &ground_pos, const size_t &num_surf_point , const std::shared_ptr<Eigen::MatrixXd>& init_points_ptr);
  int Val(const double *disp, data_ptr<dim_> &data) const;
  int Gra(const double *disp, data_ptr<dim_> &data) const;
  int Hes(const double *disp, data_ptr<dim_> &data) const;
  size_t Nx() const;
 private:
  const char ground_axis_;
  const double w_coll_;
  const double ground_pos_;
  const size_t num_surf_point_;
  const size_t dof_;
  const std::shared_ptr<Eigen::MatrixXd> init_points_ptr_;

};

template<size_t dim_>
class momentum : public  Functional<double, dim_>{
 public:
  //used for displacement based 
  momentum(const size_t dof,const Eigen::VectorXd& mass_vec, const double& dt);

  //used for position based
  momentum(const double* x, const size_t dof,const Eigen::VectorXd& mass_vec, const double& dt);
  int Val(const double *disp, data_ptr<dim_> &data) const ;
  int Gra(const double *disp, data_ptr<dim_> &data) const ;
  int Hes(const double *disp, data_ptr<dim_> &data) const ;
  int update_location_and_velocity(const double *new_dispk_ptr);
  size_t Nx() const;
 private:
  Eigen::VectorXd vk_, dispk_;
  const size_t dof_;
  Eigen::VectorXd mass_vec_;
  const double dt_;
  const double d1dt_;
  const double d1dtdt_;

};

template class position_constraint<3>;
template class momentum<3>;
template class gravity_energy<3>;
template class collision<3>;

}//namespace marvel
#endif
