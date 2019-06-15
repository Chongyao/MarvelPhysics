#include "mtr_sampler.h"

#include <iostream>
#include <minpack.h>
#include <Eigen/Eigen>
#include <zjucad/matrix/itr_matrix.h>

#include "bigbang/config.h"

using namespace std;
using namespace zjucad::matrix;
using namespace Eigen;

namespace bigbang {

extern "C" void quad4_shape_func_val_(double *val, const double *eps);
extern "C" void quad4_shape_func_jac_(double *jac, const double *eps);
extern "C" void hex8_shape_func_(double *val, const double *eps);
extern "C" void hex8_shape_func_jac_(double *jac, const double *eps);

struct bilinear_interp_prb {
  int m, n;
  double *verts_and_pos;
};

struct trilinear_interp_prb {
  int m, n;
  double *verts_and_pos;
};

static void lmder_cb(int *m, int *n, double *x, double *fvec, double *fjac, int *ldfjac, int *iflag) {
  assert(*m == 2 && *n == 2);
  
  const double *VERT_AND_POS = *(double**)(n+1);
  Eigen::Map<const MatrixXd> VP(VERT_AND_POS, 2, 5);

  Vector4d ws;
  Matrix<double, 4, 2> ws_jac;
  
  if ( *iflag == 1 ) {
    quad4_shape_func_val_(ws.data(), x);
    Eigen::Map<VectorXd>(fvec, 2) = VP.topLeftCorner(2, 4)*ws-VP.col(4);
  } else if ( *iflag == 2 ) {
    quad4_shape_func_jac_(ws_jac.data(), x);
    Eigen::Map<MatrixXd>(fjac, 2, 2) = VP.topLeftCorner(2, 4)*ws_jac;
  }
}

static void lmder_cb2(int *m, int *n, double *x, double *fvec, double *fjac, int *ldfjac, int *iflag) {
  assert(*m == 3 && *n == 3);
  
  const double *VERT_AND_POS = *(double**)(n+1);
  Eigen::Map<const MatrixXd> VP(VERT_AND_POS, 3, 9);

  Matrix<double, 8, 1> ws;
  Matrix<double, 8, 3> ws_jac;

  if ( *iflag == 1 ) {
    hex8_shape_func_(ws.data(), x);
    Eigen::Map<VectorXd>(fvec, 3) = VP.topLeftCorner(3, 8)*ws-VP.col(8);
  } else if ( *iflag == 2 ) {
    hex8_shape_func_jac_(ws_jac.data(), x);
    Eigen::Map<MatrixXd>(fjac, 3, 3) = VP.topLeftCorner(3, 8)*ws_jac;
  }
}

int calc_inv_bilinear_interp_2d(const matd_t &verts, const matd_t &pos,
                                matd_t &eps, const size_t maxits) {
  assert(verts.size(1) == 2 && verts.size(2) == 4);
  assert(pos.size() == 2);
  
  matd_t verts_and_pos = zeros<double>(2, 5);
  verts_and_pos(colon(), colon(0, 3)) = verts;
  verts_and_pos(colon(), 4) = pos;
  
  int m = 2, n = 2;
  double fvec[2], fjac[4];
  int ipvt[2];
  int info;
  double tol = 1e-8;

  int lwa = 5*n+m+5;
  double wa[lwa];

  bilinear_interp_prb prb;
  prb.m = m;
  prb.n = n;
  prb.verts_and_pos = &verts_and_pos[0];

  for (size_t i = 0; i < maxits; ++i) {
    lmder1_(lmder_cb, &prb.m, &prb.n, &eps[0], fvec, fjac, &m, &tol, &info, &ipvt[0], &wa[0], &lwa);
    if ( info < 4 )
      break;
  }
  return info;  
}

int calc_inv_trilinear_interp_3d(const matd_t &verts, const matd_t &pos,
                                 matd_t &eps, const size_t maxits) {
  assert(verts.size(1) == 3 && verts.size(2) == 8);
  assert(pos.size() == 3);
  
  matd_t verts_and_pos = zeros<double>(3, 9);
  verts_and_pos(colon(), colon(0, 7)) = verts;
  verts_and_pos(colon(), 8) = pos;
  
  int m = 3, n = 3;
  double fvec[3], fjac[9];
  int ipvt[3];
  int info;
  double tol = 1e-8;

  int lwa = 5*n+m+5;
  double wa[lwa];

  trilinear_interp_prb prb;
  prb.m = m;
  prb.n = n;
  prb.verts_and_pos = &verts_and_pos[0];

  for (size_t i = 0; i < maxits; ++i) {
    lmder1_(lmder_cb2, &prb.m, &prb.n, &eps[0], fvec, fjac, &m, &tol, &info, &ipvt[0], &wa[0], &lwa);
    if ( info < 4 )
      break;
  }
  return info;  
}

// verts with same pos has same index in
// coarse patch and fine patch
square_mtr_sampler::square_mtr_sampler(const mati_t &patch_H, const mati_t &patch_h,
                                       const matd_t &patch_mtr_h, const matd_t &nods_h)
    : patch_mtr_h_(patch_mtr_h), bc_param_(zeros<double>(2, patch_h.size(2))) {
  ASSERT(patch_H.size() == 4 && patch_h.size(1) == 4);

  //-> corner vertices of one coarse element
  const matd_t vertH = nods_h(colon(), patch_H);

  matd_t bc = zeros<double>(2, patch_h.size(2));
  for (size_t i = 0; i < patch_h.size(2); ++i) {
    const matd_t verth = nods_h(colon(), patch_h(colon(), i));
    bc(colon(), i) = verth*ones<double>(4, 1)/4.0;
    matd_t tmp_bc = zeros<double>(2, 1);
    calc_inv_bilinear_interp_2d(vertH, bc(colon(), i), tmp_bc, 10);
    bc_param_(colon(), i) = tmp_bc;
  }
}

void square_mtr_sampler::get_material(const matd_t &qr_pnt, double *qr_mtr) {
  matd_t dist = zeros<double>(bc_param_.size(2), 1);
  for (size_t i = 0; i < bc_param_.size(2); ++i)
    dist[i] = norm(qr_pnt-bc_param_(colon(), i));

  const size_t idx = std::min_element(dist.begin(), dist.end())-dist.begin();
  itr_matrix<double*>(2, 1, qr_mtr) = patch_mtr_h_(colon(), idx);
}

}
