#ifndef POINTS_ENERGY_H
#define POINTS_ENERGY_H

#include <def.h>
#include <Eigen/Core>
#include "get_nn.h"


namespace marvel{

class point_sys: public bigbang::Functional<double>{
 public:
  point_sys(const Eigen::MatrixXd &points, const double &rho, const double &vol_all, const size_t &nearest_num);
  size_t Nx() const ;
  // int Val(const double *x, double *val) const;
  // int Gra(const double *x, double *gra) const;
  // int Hes(const double *x,  std::vector<Eigen::Triplet<double>> *hes) const;
 private:
  Eigen::MatrixXd points_;
  Eigen::MatrixXi NN_;

  
  const size_t dim_;
  const double rho_;

  size_t nearest_num_;
  double vol_all_;

  double scal_fac_;

  int calc_rhoi_vi(const double *x);
  int calc_defo_gra(const double *x, double *def_gra);

  spatial_hash SH_;

  double kernel(const double &r, const double &h);
  double kernel(const Eigen::Vector3d &xj, const Eigen::Vector3d &xi, const double &h);
  double kernel(const size_t &i, const size_t &j);

  //point info
  Eigen::VectorXd mass_i_;
  Eigen::VectorXd sup_radi_;
  
  Eigen::VectorXd rho_i_;
  Eigen::VectorXd vol_i_;

  
  
};
}
#endif
