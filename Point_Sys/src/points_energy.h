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
  Eigen::VectorXd sup_radi_;
  
  size_t dim_;
  double rho_;

  size_t nearest_num_;
  double vol_all_;

  double scal_fac_;

  int calc_rhoi_vi(const double *x, const Eigen::VectorXd &rho_i, const Eigen::VectorXd &vol_i);
  int calc_defo_gra(const double *x, double *def_gra);

  spatial_hash SH_;

  double kernel(const double &h, const double &r);
  double kernel(const Eigen::Vector3d &xj, const Eigen::Vector3d &xi, const double &r);
  
};
}
#endif
