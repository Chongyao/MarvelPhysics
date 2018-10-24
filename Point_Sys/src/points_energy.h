#ifndef POINTS_ENERGY_H
#define POINTS_ENERGY_H

#include <bigbang/def.h>

typedef zjucad::matrix::matrix<size_t> mati_t;
typedef zjucad::matrix::matrix<double> matd_t;

namespace marvel{


class point_sys: public bigbang::Functional<double>{
 public:
  point_sys(const matd_t &points_, const double &rho_);
  size_t Nx() const ;
  int Val(const double *x, double *val) const;
  int Gra(const double *x, double *gra) const;
  int Hes(const double *x,  std::vector<Eigen::Triplet<double>> *hes) const;
 private:
  matd_t points;
  size_t dim;
  double rho;
};
}
#endif
