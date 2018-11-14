#ifndef BASIC_ENERGY_H
#define BASIC_ENERGY_H
#include "data_stream.h"

namespace marvel{
class position_constraint{
 public:
  position_constraint(const double &w, const vector<size_t> &cons, const size_t dim);
  int Gra(const double *disp, energy_dat &dat_str);
  int Hes(const double *disp, energy_dat &dat_str);
 private:
  const size_t dim_;
  const double w_;
  const std::vector<double> cons_;
}
}//namespace marvel
#endif
