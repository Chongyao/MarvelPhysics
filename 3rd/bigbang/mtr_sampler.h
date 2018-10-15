#ifndef MTR_SAMPLER_H
#define MTR_SAMPLER_H

#include <zjucad/matrix/matrix.h>

namespace bigbang {

using mati_t=zjucad::matrix::matrix<size_t>;
using matd_t=zjucad::matrix::matrix<double>;

// sampling material in [-1,1]x[-1,1] domain
// patch_H is the central element of a stencil
// extracted in advance.
//  1__3
//  |__|
//  0  2
class square_mtr_sampler
{
public:
  square_mtr_sampler(const mati_t &patch_H, const mati_t &patch_h,
                     const matd_t &patch_mtr_h, const matd_t &nods_h);
  void get_material(const matd_t &qr_pnt, double *qr_mtr);
private:
  matd_t bc_param_;
  const matd_t patch_mtr_h_;
};


// // sampling material in [-1,1]x[-1,1]x[-1,1]
// class cube_mtr_sampler
// {
// public:
//   cube_mtr_sampler();
//   void get_material(const matd_t &eps, matd_t &mtr);
// private:
//   matd_t bc_, bc_param_;
// };

int calc_inv_bilinear_interp_2d(const matd_t &verts, const matd_t &pos, matd_t &eps, const size_t maxits);

int calc_inv_trilinear_interp_3d(const matd_t &verts, const matd_t &pos, matd_t &eps, const size_t maxits);

}
#endif
