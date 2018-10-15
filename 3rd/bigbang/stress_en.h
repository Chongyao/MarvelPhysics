#ifndef FEM_STENCIL_H
#define FEM_STENCIL_H

#include <zjucad/matrix/matrix.h>
#include <boost/property_tree/ptree.hpp>

#include "def.h"

namespace bigbang {

using mati_t=zjucad::matrix::matrix<size_t>;
using matd_t=zjucad::matrix::matrix<double>;
class fem_stencil_energy_d3 : public Functional<double>
{
public:
  enum Material {
    LINEAR,
    STVK,
    COROTATED,
    NEOHOOKEAN,
    STANEO,
    BOWERNEO
  };
  fem_stencil_energy_d3(const std::vector<std::shared_ptr<hexs_stencil>> &elem_H,
                        const mati_t &hexs_H, const matd_t &nods_H,
                        const Material type,
                        const mati_t &hexs_h, const matd_t &nods_h, const matd_t &face_lame_h,
                        const size_t qnum_per_cube,
                        const double w);
  size_t Nx() const;
  int Val(const double *x, double *val) const;
  int Gra(const double *x, double *gra) const;
  int Hes(const double *x, std::vector<Eigen::Triplet<double>> *hes) const;
public:
  const double w_;
  const size_t dim_;
  const std::vector<std::shared_ptr<hexs_stencil>> &elem_;
  const Material type_;
  mutable matd_t stencilE_;

  mutable matd_t stress_average;
};
}
