#ifndef PHX_UTIL_H
#define PHX_UTIL_H

#include <zjucad/matrix/matrix.h>
#include "def.h"
#include "energy.h"
#include <memory>

namespace bigbang {

using mati_t=zjucad::matrix::matrix<size_t>;
using matd_t=zjucad::matrix::matrix<double>;

inline void compute_lame_coeffs(const double Ym, const double Pr,
                                double &mu, double &lambda) {
  mu = Ym/(2*(1+Pr));
  lambda = Ym*Pr/((1+Pr)*(1-2*Pr));
}

inline void compute_young_poisson(const double mu, const double lambda,
                                  double &Ym, double &Pr) {
  Pr = 1.0/(2*(1+mu/lambda));
  Ym = 2*mu*(Pr+1);
}

template <class ELAS_POTEN>
std::shared_ptr<Functional<double>>
build_elastic_energy(const std::string &mtr, const mati_t &hexs, const matd_t &nods,
                     const matd_t &params) {
  if ( mtr == "coro" )
    return std::make_shared<ELAS_POTEN>(hexs, nods, ELAS_POTEN::COROTATED, params, 1.0);
  if ( mtr == "neohookean" )
    return std::make_shared<ELAS_POTEN>(hexs, nods, ELAS_POTEN::NEOHOOKEAN, params, 1.0);
  if ( mtr == "linear" )
    return std::make_shared<ELAS_POTEN>(hexs, nods, ELAS_POTEN::LINEAR, params, 1.0);
  if ( mtr == "stvk" )
    return std::make_shared<ELAS_POTEN>(hexs, nods, ELAS_POTEN::STVK, params, 1.0);
  if ( mtr == "stabneo" )
    return std::make_shared<ELAS_POTEN>(hexs, nods, ELAS_POTEN::STANEO, params, 1.0);
  if ( mtr == "bowerneo" )
    return std::make_shared<ELAS_POTEN>(hexs, nods, ELAS_POTEN::BOWERNEO, params, 1.0);
}

int solve_harmonic_disp_3d(const mati_t &cell, const matd_t &nods,
                           const std::shared_ptr<Functional<double>> &elas,
                           const double force_scale, matd_t &har_disp, const matd_t *deformed_nods = nullptr, const bool modify_stiff = false);

}

#endif
