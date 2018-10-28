#ifndef POINTS_ENERGY_H
#define POINTS_ENERGY_H

#include <def.h>
#include <Eigen/Core>



namespace marvel{


class point_sys: public bigbang::Functional<double>{
 public:
  point_sys(const Eigen::MatrixXd &points_, const double &rho_, const double &vol_all_);
  // size_t Nx() const ;
  // int Val(const double *x, double *val) const;
  // int Gra(const double *x, double *gra) const;
  // int Hes(const double *x,  std::vector<Eigen::Triplet<double>> *hes) const;
 private:
  Eigen::MatrixXd points;
  size_t dim;
  double rho;
  size_t nearest_num;
  double vol_all;

  
  
};
}
#endif
