#ifndef BASIC_ENERGY_H
#define BASIC_ENERGY_H
#include "def.h"
#include "data_stream.h"
#include <vector>
#include <memory>
namespace marvel{


class position_constraint {
 public:
  position_constraint(const size_t dim, const double &w, const std::vector<size_t> &cons);
  int Val(const double *disp, energy_dat &dat_str) const ;  
  int Gra(const double *disp, energy_dat &dat_str) const ;
  int Hes(const double *disp, energy_dat &dat_str) const;
 private:
  const size_t dim_;
  const double w_;
  const std::vector<size_t> cons_;
};


class gravity_energy{
 public:
  gravity_energy(const size_t dim, const double &w_g, const double &gravity, const Eigen::VectorXd &mass, const char &axis);
  int Val(const double *disp, energy_dat &dat_str) const ;
  int Gra(const double *disp, energy_dat &dat_str) const;
  int Hes(const double *disp, energy_dat &dat_str) const;
 private:
  const char axis_;
  const size_t dim_;
  const double w_g_;
  const double gravity_;
  const Eigen::VectorXd mass_;
};

//simple collision with ground

class collision{
 public:
  collision(const size_t dim, const double &w_coll, const char &ground_axis, const double &ground_pos, const size_t &num_surf_point , const std::shared_ptr<Eigen::MatrixXd>& init_points_ptr);
  int Val(const double *disp, energy_dat &dat_str) const;
  int Gra(const double *disp, energy_dat &dat_str) const;
  int Hes(const double *disp, energy_dat &dat_str) const;
 private:
  const char ground_axis_;
  const double w_coll_;
  const double ground_pos_;
  const size_t num_surf_point_;
  const size_t dim_;
  const std::shared_ptr<Eigen::MatrixXd> init_points_ptr_;

};


class momentum{
 public:
  momentum(const size_t dim,const Eigen::SparseMatrix<double>& mass_sparse, const double& dt);
  int Val(const double *disp, energy_dat &dat_str) const ;
  int Gra(const double *disp, energy_dat &dat_str) const ;
  int Hes(const double *disp, energy_dat &dat_str) const ;
  int update_location_and_velocity(const double *new_dispk_ptr);

 private:
  Eigen::VectorXd vk_, dispk_;
  const size_t dim_;
  const Eigen::SparseMatrix<double>& mass_sparse_;
  const double dt_;
  const double d1dt_;
  const double d1dtdt_;

  
  
};

}//namespace marvel
#endif
