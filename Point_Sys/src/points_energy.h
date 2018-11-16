#ifndef POINTS_ENERGY_H
#define POINTS_ENERGY_H
#include <utility>
#include <src/def.h>
#include <Eigen/Core>
#include <Eigen/SparseCore>
#include "get_nn.h"
#include "data_stream.h"
namespace marvel{







// class point_sys: public bigbang::Functional<double>{
class point_sys{
 public:
  point_sys(const Eigen::MatrixXd &points, const double &rho, const double &Young, const double &Poission, const double &vol_all, const double &kv, const std::vector<std::vector<size_t>> &friends, const Eigen::VectorXd &sup_radi );
  size_t Nx() const ;

  int pre_compute(energy_dat &dat_str) const ;//calc once
  int calc_defo_gra(const double *disp, energy_dat &dat_str) const;
  int Val(const double *disp, energy_dat &dat_str)const;
  int Gra(const double *disp, energy_dat &dat_str) const;
  int gravity(const double *disp, energy_dat &dat_str, const double &gravity) const;
  double get_mass(const size_t &i) const;
  int Hessian(const double*disp, energy_dat &dat_str);
  const Eigen::SparseMatrix<double>& get_Mass_Matrix();
  const Eigen::VectorXd& get_Mass_VectorXd();
  // int Gra(const double *x, double *gra) const;
  // int Val(const double *x, double *val) const;
  // int Gra(const double *x, double *gra) const;
  // int Hes(const double *x,  std::vector<Eigen::Triplet<double>> *hes) const;
 private:
  const Eigen::MatrixXd points_;
  Eigen::SparseMatrix<double> M_;
  const size_t dim_;
  const double rho_;
  const double Poission_;
  const double Young_;
  const double kv_;
  
  double vol_all_;
  double scal_fac_;
  // mutable spatial_hash SH_;


  int calc_rhoi_vi() const;
  int calc_weig() const;
  int calc_Mass_matrix();
  double kernel(const double &r, const double &h) const;
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
