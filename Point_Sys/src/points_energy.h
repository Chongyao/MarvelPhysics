#ifndef POINTS_ENERGY_H
#define POINTS_ENERGY_H

#include <src/def.h>
#include <Eigen/Core>
#include "get_nn.h"


namespace marvel{

// class point_sys: public bigbang::Functional<double>{
class point_sys{
 public:
  point_sys(const Eigen::MatrixXd &points, const double &rho, const double &vol_all, const size_t &nearest_num);
  size_t Nx() const ;
  // int Val(const double *x, double *val) const;
  // int Gra(const double *x, double *gra) const;
  // int Hes(const double *x,  std::vector<Eigen::Triplet<double>> *hes) const;
  int calc_defo_gra(const double *x, double *def_gra) const;
 private:
  const Eigen::MatrixXd points_;
  Eigen::MatrixXi NN_;
  
  const size_t dim_;
  const double rho_;
  double vol_all_;
  double scal_fac_;
  mutable spatial_hash SH_;
  size_t nearest_num_;

  int pre_compute(const double *x) const ;
  int calc_rhoi_vi(const double *x) const;
  int calc_fri() const;
  


  double kernel(const double &r, const double &h) const;
  double kernel(const Eigen::Vector3d &xj, const Eigen::Vector3d &xi, const double &h) const;
  double kernel(const size_t &i, const size_t &j) const;

  //point info
  Eigen::VectorXd mass_i_;
  Eigen::VectorXd sup_radi_;
  mutable Eigen::VectorXd rho_i_;
  mutable Eigen::VectorXd vol_i_;

  mutable std::vector<std::vector<size_t>> friends_;
  mutable std::vector<std::vector<double>> weig_;

};
}
#endif
