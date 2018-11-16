#ifndef BASIC_ENERGY_H
#define BASIC_ENERGY_H
#include "data_stream.h"
#include <vector>
namespace marvel{
class position_constraint{
 public:
  position_constraint(const double &w, const std::vector<size_t> &cons, const size_t dim);
  int Gra(const double *disp, energy_dat &dat_str);
  int Hes(const double *disp, energy_dat &dat_str);
 private:
  const size_t dim_;
  const double w_;
  const std::vector<size_t> cons_;
};

class gravity_energy{
 public:
  gravity_energy(const double &w_g, const double &gravity, const size_t &dim, const Eigen::VectorXd &mass, const char &axis);
  int Val(const double *disp, energy_dat &dat_str);
  int Gra(const double *disp, energy_dat &dat_str);
  int Hes(const double *disp, energy_dat &dat_str);
 private:
  const char axis_;
  const size_t dim_;
  const double w_g_;
  const double gravity_;
  const Eigen::VectorXd mass_;
};
}//namespace marvel
#endif
