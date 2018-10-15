
#include "stencil.h"

#include <iostream>
#include <zjucad/matrix/itr_matrix.h>
#include <zjucad/matrix/io.h>
#include <hjlib/math/polar.h>
#include <chrono>

#include "util.h"
#include "config.h"
#include "mtr_sampler.h"
#include "gauss_quadrature.h"

using namespace std;
using namespace zjucad::matrix;
using namespace Eigen;

namespace bigbang {

extern "C" {
  void quad_neo_f_val_at_qr_(double *val, const double *F, const double *mu, const double *lam);
  void quad_neo_f_val_at_qr_jac_(double *jac, const double *F, const double *mu, const double *lam);
  void quad_neo_f_val_at_qr_hes_(double *hes, const double *F, const double *mu, const double *lam);

  void quad_stvk_f_val_at_qr_(double *val, const double *F, const double *mu, const double *lam);
  void quad_stvk_f_val_at_qr_jac_(double *jac, const double *F, const double *mu, const double *lam);
  void quad_stvk_f_val_at_qr_hes_(double *hes, const double *F, const double *mu, const double *lam);

  void quad_coro_f_val_at_qr_(double *val, const double *F, const double *R, const double *mu, const double *lam);
  void quad_coro_f_val_at_qr_jac_(double *jac, const double *F, const double *R, const double *mu, const double *lam);
  void quad_coro_f_val_at_qr_hes_(double *hes, const double *F, const double *R, const double *mu, const double *lam);

  void hex_coro_f_at_qr_(double *val, const double *F, const double *R, const double *mu, const double *lam);
  void hex_coro_f_at_qr_jac_(double *jac, const double *F, const double *R, const double *mu, const double *lam);
  void hex_coro_f_at_qr_hes_(double *hes, const double *F, const double *R, const double *mu, const double *lam);

  void hex_stvk_f_at_qr_(double *val, const double *F, const double *mu, const double *lam);
  void hex_stvk_f_at_qr_jac_(double *jac, const double *F, const double *mu, const double *lam);
  void hex_stvk_f_at_qr_hes_(double *hes, const double *F, const double *mu, const double *lam);

  void hex_neo_f_at_qr_(double *val, const double *F, const double *mu, const double *lam);
  void hex_neo_f_at_qr_jac_(double *jac, const double *F, const double *mu, const double *lam);
  void hex_neo_f_at_qr_hes_(double *hes, const double *F, const double *mu, const double *lam);

  void hex_sta_neo_f_at_qr_(double *val, const double *F, const double *mu, const double *lam);
  void hex_sta_neo_f_at_qr_jac_(double *jac, const double *F, const double *mu, const double *lam);
  void hex_sta_neo_f_at_qr_hes_(double *hes, const double *F, const double *mu, const double *lam);

  void hex_bower_neo_f_at_qr_(double *val, const double *F, const double *mu, const double *lam);
  void hex_bower_neo_f_at_qr_jac_(double *jac, const double *F, const double *mu, const double *lam);
  void hex_bower_neo_f_at_qr_hes_(double *hes, const double *F, const double *mu, const double *lam);
}

extern "C" {
  void quad4_shape_func_val_(double *val, const double *eps);
  void quad4_shape_func_jac_(double *jac, const double *eps);

  void quad9_shape_func_val_(double *val, const double *eps);
  void quad9_shape_func_jac_(double *jac, const double *eps);

  void hex8_shape_func_(double *val, const double *eps);
  void hex8_shape_func_jac_(double *jac, const double *eps);

  void hex27_shape_func_(double *val, const double *eps);
  void hex27_shape_func_jac_(double *jac, const double *eps);
}

void full_quad_dis_shape_func_val(double *val, const double *eps, const matd_t &basisXY) {
  const size_t basis_num = basisXY.size(1), sample_num = basisXY.size(2);
  const size_t sample_num_axis = sqrt(sample_num);
  ASSERT(sample_num_axis*sample_num_axis == sample_num);

  const double spacing = 2.0/(sample_num_axis-1);
  const matd_t sample_pos = -1*ones<double>(sample_num_axis, 1)+spacing*matd_t(colon(0, sample_num_axis-1));

  const size_t loc[2] = {(size_t)std::floor((eps[0]+1)/spacing), (size_t)std::floor((eps[1]+1)/spacing)};
  const double param_xy[2] = {(2*eps[0]-sample_pos[loc[0]]-sample_pos[loc[0]+1])/spacing,
                              (2*eps[1]-sample_pos[loc[1]]-sample_pos[loc[1]+1])/spacing};

  matd_t bi_coord = zeros<double>(4, 1);
  quad4_shape_func_val_(&bi_coord[0], param_xy);

  mati_t four_adj_node(4, 1);
  four_adj_node[0] = loc[0]*sample_num_axis+loc[1];
  four_adj_node[1] = four_adj_node[0]+1;
  four_adj_node[2] = (loc[0]+1)*sample_num_axis+loc[1];
  four_adj_node[3] = four_adj_node[2]+1;

  itr_matrix<double *>(basis_num, 1, val) = basisXY(colon(), four_adj_node)*bi_coord;
}

void full_quad_dis_shape_func_jac(double *jac, const double *eps, const matd_t &basisXY) {
  const size_t basis_num = basisXY.size(1), sample_num = basisXY.size(2);
  const size_t sample_num_axis = sqrt(sample_num);
  ASSERT(sample_num_axis*sample_num_axis == sample_num);

  const double spacing = 2.0/(sample_num_axis-1);
  const matd_t sample_pos = -1*ones<double>(sample_num_axis, 1)+spacing*matd_t(colon(0, sample_num_axis-1));

  const size_t loc[2] = {(size_t)std::floor((eps[0]+1)/spacing), (size_t)std::floor((eps[1]+1)/spacing)};
  const double param_xy[2] = {(2*eps[0]-sample_pos[loc[0]]-sample_pos[loc[0]+1])/spacing,
                              (2*eps[1]-sample_pos[loc[1]]-sample_pos[loc[1]+1])/spacing};

  matd_t bi_coord_jac = zeros<double>(4, 2);
  quad4_shape_func_jac_(&bi_coord_jac[0], param_xy);

  mati_t four_adj_node(4, 1);
  four_adj_node[0] = loc[0]*sample_num_axis+loc[1];
  four_adj_node[1] = four_adj_node[0]+1;
  four_adj_node[2] = (loc[0]+1)*sample_num_axis+loc[1];
  four_adj_node[3] = four_adj_node[2]+1;

  itr_matrix<double *>(basis_num, 2, jac) = basisXY(colon(), four_adj_node)*bi_coord_jac*2/spacing;
}

void full_hexs_dis_shape_func_val(double *val, const double *eps, const matd_t &basisXYZ) {
  const size_t basis_num = basisXYZ.size(1), sample_num = basisXYZ.size(2);
  const size_t sample_num_axis = cbrt(sample_num);
  ASSERT(sample_num_axis*sample_num_axis*sample_num_axis == sample_num);

  const double spacing = 2.0/(sample_num_axis-1);
  const matd_t sample_pos = -1*ones<double>(sample_num_axis, 1)+spacing*matd_t(colon(0, sample_num_axis-1));

  const size_t loc[3] = {(size_t)std::floor((eps[0]+1)/spacing),
                         (size_t)std::floor((eps[1]+1)/spacing),
                         (size_t)std::floor((eps[2]+1)/spacing)};
  const double param_xyz[3] = {(2*eps[0]-sample_pos[loc[0]]-sample_pos[loc[0]+1])/spacing,
                               (2*eps[1]-sample_pos[loc[1]]-sample_pos[loc[1]+1])/spacing,
                               (2*eps[2]-sample_pos[loc[2]]-sample_pos[loc[2]+1])/spacing};

  matd_t tri_coord = zeros<double>(8, 1);
  hex8_shape_func_(&tri_coord[0], param_xyz);

  mati_t eight_adj_node(8, 1); {
    eight_adj_node[0] = loc[0]*sample_num_axis*sample_num_axis+loc[1]*sample_num_axis+loc[2];
    eight_adj_node[1] = eight_adj_node[0]+1;
    eight_adj_node[2] = loc[0]*sample_num_axis*sample_num_axis+(loc[1]+1)*sample_num_axis+loc[2];
    eight_adj_node[3] = eight_adj_node[2]+1;
    eight_adj_node[4] = (loc[0]+1)*sample_num_axis*sample_num_axis+loc[1]*sample_num_axis+loc[2];
    eight_adj_node[5] = eight_adj_node[4]+1;
    eight_adj_node[6] = (loc[0]+1)*sample_num_axis*sample_num_axis+(loc[1]+1)*sample_num_axis+loc[2];
    eight_adj_node[7] = eight_adj_node[6]+1; 
  }

  itr_matrix<double *>(basis_num, 1, val) = basisXYZ(colon(), eight_adj_node)*tri_coord;
}

void full_hexs_dis_shape_func_jac(double *jac, const double *eps, const matd_t &basisXYZ) {
  const size_t basis_num = basisXYZ.size(1), sample_num = basisXYZ.size(2);
  const size_t sample_num_axis = cbrt(sample_num);
  ASSERT(sample_num_axis*sample_num_axis*sample_num_axis == sample_num);

  const double spacing = 2.0/(sample_num_axis-1);
  const matd_t sample_pos = -1*ones<double>(sample_num_axis, 1)+spacing*matd_t(colon(0, sample_num_axis-1));

  const size_t loc[3] = {(size_t)std::floor((eps[0]+1)/spacing),
                         (size_t)std::floor((eps[1]+1)/spacing),
                         (size_t)std::floor((eps[2]+1)/spacing)};
  const double param_xyz[3] = {(2*eps[0]-sample_pos[loc[0]]-sample_pos[loc[0]+1])/spacing,
                               (2*eps[1]-sample_pos[loc[1]]-sample_pos[loc[1]+1])/spacing,
                               (2*eps[2]-sample_pos[loc[2]]-sample_pos[loc[2]+1])/spacing};

  matd_t tri_coord_jac = zeros<double>(8, 3);
  hex8_shape_func_jac_(&tri_coord_jac[0], param_xyz);

  mati_t eight_adj_node(8, 1); {
    eight_adj_node[0] = loc[0]*sample_num_axis*sample_num_axis+loc[1]*sample_num_axis+loc[2];
    eight_adj_node[1] = eight_adj_node[0]+1;
    eight_adj_node[2] = loc[0]*sample_num_axis*sample_num_axis+(loc[1]+1)*sample_num_axis+loc[2];
    eight_adj_node[3] = eight_adj_node[2]+1;
    eight_adj_node[4] = (loc[0]+1)*sample_num_axis*sample_num_axis+loc[1]*sample_num_axis+loc[2];
    eight_adj_node[5] = eight_adj_node[4]+1;
    eight_adj_node[6] = (loc[0]+1)*sample_num_axis*sample_num_axis+(loc[1]+1)*sample_num_axis+loc[2];
    eight_adj_node[7] = eight_adj_node[6]+1; 
  }

  itr_matrix<double *>(basis_num, 3, jac) = basisXYZ(colon(), eight_adj_node)*tri_coord_jac*2/spacing;
}

//===============================================================================
quad_stencil::quad_stencil(const mati_t &adjc_elem, const vector<size_t> &faces)
    : adjc_elem_(adjc_elem) {
  adjc_vert_num_ = adjc_elem_.size();

  //-> record face index in the stencil
  faces_.resize(faces.size(), 1);
  std::copy(faces.begin(), faces.end(), faces_.begin());

  //-> set default shape functions for quad stencil
  if ( adjc_vert_num_ == 4 )
    fem_sf_jac_ = rest_sf_jac_ = quad4_shape_func_jac_;
  else if ( adjc_vert_num_ == 9 )
    fem_sf_jac_ = rest_sf_jac_ = quad9_shape_func_jac_;
  else
    ASSERT(0);

  this->dis_fem_sf_jac_  = nullptr;
  this->isotropic_       = true;
}

int quad_stencil::set_full_discrete_shape_func(const matd_t &basisXY, const bool isotropic,
                                               const bool augmented) {
  if ( basisXY.size(1) != adjc_vert_num_ ) {
    ASSERT(0);
    return __LINE__;
  }
  
  this->dis_fem_sf_jac_ = full_quad_dis_shape_func_jac;
  this->basisXY_   = basisXY;
  this->isotropic_ = isotropic;
  this->augmented_ = augmented;

  return 0;
}

int quad_stencil::set_def_grad_oper(const matd_t &adjc_nods) {
  //-> save rest shape
  adjc_rest_ = adjc_nods;
  
  if ( this->isotropic_ ) {
    //-> high order basis part
    adjc_H_invDmH_.resize(adjc_vert_num_, 2*qr_num_);
    dets_.resize(1, qr_num_);

    //-> one for bases for disp, the other for rest shape
    dNude_ = zeros<double>(adjc_vert_num_, 2*qr_num_);
    dNxde_ = zeros<double>(adjc_vert_num_, 2*qr_num_);
    DmHs_ = zeros<double>(2, 2*qr_num_);

    for (size_t i = 0; i < qr_num_; ++i) {
      matd_t H = zeros<double>(adjc_vert_num_, 2);
      if ( dis_fem_sf_jac_ )
        dis_fem_sf_jac_(&H[0], &qrs_(0, i), basisXY_);
      else 
        fem_sf_jac_(&H[0], &qrs_(0, i));
      dNude_(colon(), colon(2*i, 2*i+1)) = H;

      matd_t H0 = zeros<double>(adjc_vert_num_, 2);
      rest_sf_jac_(&H0[0], &qrs_(0, i));
      dNxde_(colon(), colon(2*i, 2*i+1)) = H0;      
      matd_t DmH = adjc_nods*H0;
      matd_t cpDmH = DmH;
      dets_[i] = fabs(det(cpDmH));

      if ( inv(DmH) ) std::cerr << "# inv fail" << std::endl;
      DmHs_(colon(), colon(2*i, 2*i+1)) = DmH;
      
      adjc_H_invDmH_(colon(), colon(2*i, 2*i+1)) = H*DmH;
    }

    adjc_df_op_ = zeros<double>(2*adjc_vert_num_, 4*qr_num_);
    for (size_t i = 0; i < qr_num_; ++i) {
      adjc_df_op_(colon(), colon(4*i, 4*i+3)) = kroneckerId<double, 2>(adjc_H_invDmH_(colon(), colon(2*i, 2*i+1)));
    }
  } else if ( !this->isotropic_ && !this->augmented_ ) {
    ASSERT(dis_fem_sf_jac_);
    
    //-> anistropic bases
    const size_t b_cols = basisXY_.size(2)/2;

    const matd_t basisX = basisXY_(colon(), colon(0*b_cols, 1*b_cols-1));
    const matd_t basisY = basisXY_(colon(), colon(1*b_cols, 2*b_cols-1));

    adjc_df_op_ = zeros<double>(2*adjc_vert_num_, 4*qr_num_);
    dets_.resize(1, qr_num_);

    //-> jacobian of bases for u and X
    dNude_ = zeros<double>(2*adjc_vert_num_, 2*qr_num_);
    dNxde_ = zeros<double>(adjc_vert_num_, 2*qr_num_);
    DmHs_ = zeros<double>(2, 2*qr_num_);
                         
    #pragma omp parallel for
    for (size_t i = 0; i < qr_num_; ++i) {
      matd_t Hx(adjc_vert_num_, 2), Hy(adjc_vert_num_, 2);
      {
        dis_fem_sf_jac_(&Hx[0], &qrs_(0, i), basisX);
        dis_fem_sf_jac_(&Hy[0], &qrs_(0, i), basisY);
      }

      matd_t H0 = zeros<double>(adjc_vert_num_, 2);
      rest_sf_jac_(&H0[0], &qrs_(0, i));
      dNxde_(colon(), colon(2*i, 2*i+1)) = H0;

      matd_t DmH = adjc_nods*H0;
      matd_t cpDmH = DmH;
      dets_[i] = fabs(det(cpDmH));

      if ( inv(DmH) ) std::cerr << "# inv fail" << std::endl;
      DmHs_(colon(), colon(2*i, 2*i+1)) = DmH;

      matd_t temp_df_op = zeros<double>(4, 2*adjc_vert_num_); {
        for (size_t p = 0; p < 2; ++p) {
          for (size_t q = 0; q < adjc_vert_num_; ++q) {
            temp_df_op(2*p+0, 2*q+0) = dot(Hx(q, colon()), DmH(colon(), p));
            temp_df_op(2*p+1, 2*q+1) = dot(Hy(q, colon()), DmH(colon(), p));
          }
        }
      }
      adjc_df_op_(colon(), colon(4*i, 4*i+3)) = trans(temp_df_op);

      //-> prepare for coro bases deformation gradient operator
      //-> only need to save Dm and dNdE
      for (size_t k = 0; k < adjc_vert_num_; ++k) {
        dNude_(2*k+0, colon(2*i, 2*i+1)) = Hx(k, colon());
        dNude_(2*k+1, colon(2*i, 2*i+1)) = Hy(k, colon());
      }
    }
  } else if ( !this->isotropic_ && this->augmented_ ) {
    ASSERT(dis_fem_sf_jac_);
    
    //-> anistropic bases
    const size_t b_cols = basisXY_.size(2)/(2*2);
    ASSERT(basisXY_.size(2)%4 == 0);

    const matd_t basis00 = basisXY_(colon(), colon(0*b_cols, 1*b_cols-1));
    const matd_t basis01 = basisXY_(colon(), colon(1*b_cols, 2*b_cols-1));
    const matd_t basis10 = basisXY_(colon(), colon(2*b_cols, 3*b_cols-1));
    const matd_t basis11 = basisXY_(colon(), colon(3*b_cols, 4*b_cols-1));

    adjc_df_op_ = zeros<double>(2*adjc_vert_num_, 4*qr_num_);
    dets_.resize(1, qr_num_);

    //-> jacobian of bases for u and X
    dNude_ = zeros<double>(4*adjc_vert_num_, 2*qr_num_);
    dNxde_ = zeros<double>(adjc_vert_num_, 2*qr_num_);
    DmHs_ = zeros<double>(2, 2*qr_num_);
                         
    #pragma omp parallel for
    for (size_t i = 0; i < qr_num_; ++i) {
      matd_t H00(adjc_vert_num_, 2), H01(adjc_vert_num_, 2);
      matd_t H10(adjc_vert_num_, 2), H11(adjc_vert_num_, 2);
      {
        dis_fem_sf_jac_(&H00[0], &qrs_(0, i), basis00);
        dis_fem_sf_jac_(&H01[0], &qrs_(0, i), basis01);
        dis_fem_sf_jac_(&H10[0], &qrs_(0, i), basis10);
        dis_fem_sf_jac_(&H11[0], &qrs_(0, i), basis11);
      }

      matd_t H0 = zeros<double>(adjc_vert_num_, 2);
      rest_sf_jac_(&H0[0], &qrs_(0, i));
      dNxde_(colon(), colon(2*i, 2*i+1)) = H0;
      matd_t DmH = adjc_nods*H0;
      matd_t cpDmH = DmH;
      dets_[i] = fabs(det(cpDmH));

      if ( inv(DmH) ) std::cerr << "# inv fail" << std::endl;
      DmHs_(colon(), colon(2*i, 2*i+1)) = DmH;

      //-> prepare for coro bases deformation gradient operator
      //-> only need to save Dm and dNdE
      for (size_t k = 0; k < adjc_vert_num_; ++k) {
        dNude_(4*k+0, colon(2*i, 2*i+1)) = H00(k, colon());
        dNude_(4*k+1, colon(2*i, 2*i+1)) = H10(k, colon());
        dNude_(4*k+2, colon(2*i, 2*i+1)) = H01(k, colon());
        dNude_(4*k+3, colon(2*i, 2*i+1)) = H11(k, colon());
      }
    }    
  }
  return 0;
}

int quad_stencil::set_mtrs_and_qrs(const size_t per_face_to_many,
                                   const mati_t &quad_h, const matd_t &nods_h,
                                   const matd_t &lame_h) {
  const size_t num = 2*sqrt(per_face_to_many*faces_.size());
  const double *qr_ptn, *qr_wgts;
  switch ( num ) {
    case 2: qr_ptn = G_QR_PNT_2; qr_wgts = G_QR_WGT_2; break;
    case 4: qr_ptn = G_QR_PNT_4; qr_wgts = G_QR_WGT_4; break;
    case 6: qr_ptn = G_QR_PNT_6; qr_wgts = G_QR_WGT_6; break;
    case 8: qr_ptn = G_QR_PNT_8; qr_wgts = G_QR_WGT_8; break;
    default: ASSERT(0);
  }
  qr_num_ = num*num;
  qrs_ = zeros<double>(2, qr_num_);
  qrw_ = zeros<double>(1, qr_num_);
  for (size_t i = 0; i < num; ++i) {
    for (size_t j = 0; j < num; ++j) {
      const size_t idx = i*num+j;
      qrs_(0, idx) = qr_ptn[i];
      qrs_(1, idx) = qr_ptn[j];
      qrw_[idx] = qr_wgts[i]*qr_wgts[j];
    }
  }
  
  mati_t fine_ids(faces_.size()*per_face_to_many, 1);
  for (size_t i = 0; i < faces_.size(); ++i) {
    fine_ids(colon(i*per_face_to_many, (i+1)*per_face_to_many-1))
        = colon(faces_[i]*per_face_to_many, (faces_[i]+1)*per_face_to_many-1);
  }
  const mati_t patch_h = quad_h(colon(), fine_ids);
  const matd_t patch_mtr_h = lame_h(colon(), fine_ids);

  const size_t v_stride = sqrt(adjc_elem_.size())-1;
  matd_t corner_vert = zeros<double>(nods_h.size(1), 4);
  for (size_t i = 0; i < 2; ++i) {
    for (size_t j = 0; j < 2; ++j) {
      const size_t idx = 2*i+j;
      const size_t off = i*v_stride*(v_stride+1)+j*v_stride;
      corner_vert(colon(), idx) = nods_h(colon(), adjc_elem_[off]);
    }
  }

  /// paramterize all BCs into [-1,1]x[-1,1]
  matd_t bc = zeros<double>(2, patch_h.size(2)), param_bc = bc;
  for (size_t i = 0; i < patch_h.size(2); ++i) {
    bc(colon(), i) = nods_h(colon(), patch_h(colon(), i))*ones<double>(4, 1)/4;
    matd_t tmp_res = zeros<double>(2, 1);
    calc_inv_bilinear_interp_2d(corner_vert, bc(colon(), i), tmp_res, 10);
    param_bc(colon(), i) = tmp_res;
  }
      
  lame_on_qrs_ = zeros<double>(lame_h.size(1), qr_num_);
  for (size_t j = 0; j < qr_num_; ++j) {
    matd_t dist = zeros<double>(param_bc.size(2), 1);
    for (size_t k = 0; k < param_bc.size(2); ++k)
      dist[k] = norm(qrs_(colon(), j)-param_bc(colon(), k));

    const size_t idx = std::min_element(dist.begin(), dist.end())-dist.begin();
    lame_on_qrs_(colon(), j) = patch_mtr_h(colon(), idx);
    // cout << qrs_(colon(), j) << endl;
    // cout << lame_on_qrs_(colon(), j) << endl << endl;
    // getchar();
  }

  // static matd_t buffer_lame = lame_on_qrs_;
  // static int cnt = 0;
  // cout << cnt++ << " " << norm(lame_on_qrs_-buffer_lame) << endl;
  // getchar();

  return 0;
}

extern "C" {
  /* x, x0 is for a single vert */
  void coro_ani_basis_df_2d_(double *val, const double *x, const double *x0, const double *R,
                             const double *Nu, const double *Nx, const double *dNude,
                             const double *dNxde, const double *Dm);
  void coro_ani_basis_df_2d_jac_(double *jac, const double *x, const double *x0, const double *R,
                                 const double *Nu, const double *Nx, const double *dNude,
                                 const double *dNxde, const double *Dm);
  void coro_iso_basis_df_2d_(double *val, const double *x, const double *x0, const double *R,
                             const double *Nu, const double *Nx, const double *dNude,
                             const double *dNxde, const double *Dm);
  void coro_iso_basis_df_2d_jac_(double *jac, const double *x, const double *x0, const double *R,
                                 const double *Nu, const double *Nx, const double *dNude,
                                 const double *dNxde, const double *Dm);
  void coro_aug_ani_df_2d_(double *val, const double *x, const double *x0, const double *R,
                           const double *Nu, const double *Nx, const double *dNude,
                           const double *dNxde, const double *Dm);
  void coro_aug_ani_df_2d_jac_(double *val, const double *x, const double *x0, const double *R,
                               const double *Nu, const double *Nx, const double *dNude,
                               const double *dNxde, const double *Dm);
}

int quad_stencil::query_local_frame(const matd_t &curr_nods, const size_t qrID,
                                    matd_t &R) {
  // ASSERT(qrID >= 0 && qrID < qr_num_);

  // matd_t H = zeros<double>(adjc_vert_num_, 2);

  // if ( this->isotropic_ ) {
  //   //-> isotropic case
  //   if ( dis_fem_sf_jac_ )
  //     dis_fem_sf_jac_(&H[0], &qrs_(0, qrID), basisXY_);
  //   else 
  //     fem_sf_jac_(&H[0], &qrs_(0, qrID));
  // } else {
  //   //-> anistropic case
  //   rest_sf_jac_(&H[0], &qrs_(0, qrID));
  // }

  // //-> polar the current def gradient
  // R = (curr_nods-adjc_rest_)*H*DmHs_(colon(), colon(2*qrID, 2*qrID+1))+eye<double>(2);
  // hj::polar2d rs_decomp;
  // rs_decomp(R);

  matd_t H = zeros<double>(adjc_vert_num_, 2);
  const double eps[2] = {0};
  if ( this->isotropic_ ) {
    //-> isotropic case
    if ( dis_fem_sf_jac_ )
      dis_fem_sf_jac_(&H[0], eps, basisXY_);
    else 
      fem_sf_jac_(&H[0], eps);
  } else {
    //-> anistropic case
    rest_sf_jac_(&H[0], eps);
  }

  matd_t Dm = adjc_rest_*H;
  inv(Dm);
  R = (curr_nods*H)*Dm;
  hj::polar2d rs_decomp;
  rs_decomp(R);

  return 0;
}

int quad_stencil::get_coro_bases_df(const matd_t &curr_nods, const matd_t &R,
                                    const size_t qrID, matd_t &dF, matd_t &G) {
  ASSERT(qrID >= 0 && qrID < qr_num_);
  
  dF = zeros<double>(2, 2);
  G = zeros<double>(4, 2*adjc_vert_num_);
  matd_t tmpF = zeros<double>(2, 2);
  matd_t tmpG = zeros<double>(4, 2);

  //-> summed over each vert
  if ( !this->isotropic_ && !this->augmented_ ) {
    matd_t dNude = zeros<double>(2, 2);
    matd_t dNxde = zeros<double>(1, 2);
    for (size_t i = 0; i < adjc_vert_num_; ++i) {
      dNude = dNude_(colon(2*i+0, 2*i+1), colon(2*qrID, 2*qrID+1));
      dNxde = dNxde_(i, colon(2*qrID, 2*qrID+1));

      coro_ani_basis_df_2d_(&tmpF[0], &curr_nods(0, i), &adjc_rest_(0, i), &R[0],
                            nullptr, nullptr, &dNude[0], &dNxde[0], &DmHs_(0, 2*qrID));      
      dF += tmpF;

      coro_ani_basis_df_2d_jac_(&tmpG[0], nullptr, nullptr, &R[0],
                                nullptr, nullptr, &dNude[0], &dNxde[0], &DmHs_(0, 2*qrID));
      G(colon(), colon(2*i+0, 2*i+1)) = tmpG;
    }
  } else if ( !this->isotropic_ && this->augmented_ ) {
    matd_t dNude = zeros<double>(4, 2);
    matd_t dNxde = zeros<double>(1, 2);
    for (size_t i = 0; i < adjc_vert_num_; ++i) {
      dNude = dNude_(colon(4*i+0, 4*i+3), colon(2*qrID, 2*qrID+1));
      dNxde = dNxde_(i, colon(2*qrID, 2*qrID+1));

      coro_aug_ani_df_2d_(&tmpF[0], &curr_nods(0, i), &adjc_rest_(0, i), &R[0],
                          nullptr, nullptr, &dNude[0], &dNxde[0], &DmHs_(0, 2*qrID));  
      dF += tmpF;

      coro_aug_ani_df_2d_jac_(&tmpG[0], nullptr, nullptr, &R[0],
                              nullptr, nullptr, &dNude[0], &dNxde[0], &DmHs_(0, 2*qrID));
      G(colon(), colon(2*i+0, 2*i+1)) = tmpG;
    }
  } else if ( this->isotropic_ ) {
    matd_t dNude = zeros<double>(1, 2);
    matd_t dNxde = zeros<double>(1, 2);
    for (size_t i = 0; i < adjc_vert_num_; ++i) {
      dNude = dNude_(i, colon(2*qrID, 2*qrID+1));
      dNxde = dNxde_(i, colon(2*qrID, 2*qrID+1));

      coro_iso_basis_df_2d_(&tmpF[0], &curr_nods(0, i), &adjc_rest_(0, i), &R[0],
                            nullptr, nullptr, &dNude[0], &dNxde[0], &DmHs_(0, 2*qrID));
      dF += tmpF;
      
      coro_iso_basis_df_2d_jac_(&tmpG[0], nullptr, &adjc_rest_(0, i), &R[0],
                                nullptr, nullptr, &dNude[0], &dNxde[0], &DmHs_(0, 2*qrID));
      G(colon(), colon(2*i+0, 2*i+1)) = tmpG;
    }
  }
  
  dF += eye<double>(2);
  return 0;
}

void quad_stencil::apply_local_LMA(const mati_t &quad_H, const matd_t &nods_H, const size_t eignum,
                                   const SparseMatrix<double> &K_h,
                                   const SparseMatrix<double> &M_h,
                                   const SparseMatrix<double> &P,
                                   const boost::property_tree::ptree &pt, 
                                   VectorXd &eigs_H,
                                   VectorXd &rqs_H) {
  //-> should be called after the energy is established
  ASSERT(adjc_rest_.size() > 0 && lame_on_qrs_.size() > 0);

  const string cons_law = pt.get<string>("elas.value");
  const size_t adjc_dim = 2*adjc_vert_num_;

  //-> for each quadrature points
  matd_t adjc_H = zeros<double>(adjc_dim, adjc_dim), HF = zeros<double>(4, 4);
  const matd_t R = eye<double>(2), dR = eye<double>(2);
  for (size_t j = 0; j < this->qr_num_; ++j) {
    matd_t G;
    matd_t dF = zeros<double>(2, 2);
    this->get_coro_bases_df(adjc_rest_, R, j, dF, G);
            
    HF = zeros<double>(4, 4);
    if ( cons_law == "coro" )
      quad_coro_f_val_at_qr_hes_(&HF[0], &dF[0], &dR[0], &lame_on_qrs_(0, j), &lame_on_qrs_(1, j));
    else if ( cons_law == "stvk" )
      quad_stvk_f_val_at_qr_hes_(&HF[0], &dF[0], &lame_on_qrs_(0, j), &lame_on_qrs_(1, j));
    else if ( cons_law == "neohookean" )
      quad_neo_f_val_at_qr_hes_(&HF[0], &dF[0], &lame_on_qrs_(0, j), &lame_on_qrs_(1, j));
      
    adjc_H += dets_[j]*qrw_[j]*trans(G)*HF*G;
  }

  //-> insert into a sparse K
  SparseMatrix<double> K_H = Map<const MatrixXd>(&adjc_H[0], adjc_H.size(1), adjc_H.size(2)).sparseView();

  //-> global to local mapping
  mati_t g2l = -1*ones<size_t>(nods_H.size(2), 1);
  g2l(adjc_elem_) = colon(0, adjc_elem_.size()-1);

  //-> get mass matrix
  SparseMatrix<double> M_H(adjc_H.size(1), adjc_H.size(2));
  M_H.setZero();
  for (size_t i = 0; i < faces_.size(); ++i) {
    const size_t fid = faces_[i];
    const double len = norm(nods_H(colon(), quad_H(0, fid))-nods_H(colon(), quad_H(1, fid)));
    const double area = len*len;
    for (size_t j = 0; j < quad_H.size(1); ++j) {
      const size_t pid = quad_H(j, fid);
      ASSERT(g2l[pid] != -1);
      M_H.coeffRef(2*g2l[pid]+0, 2*g2l[pid]+0) += area/4;
      M_H.coeffRef(2*g2l[pid]+1, 2*g2l[pid]+1) += area/4;
    }
  }

  //-> solve generalized eigen
  const int smallest = 0;
  MatrixXd eigv_H;
  solve_gen_eig_prob(K_H, M_H, eignum, smallest, "arpaca", &eigs_H, &eigv_H);

  //-> calc rayleigh quotient
  rqs_H = VectorXd::Zero(eigv_H.cols());
  for (size_t i = 0; i < eigv_H.cols(); ++i) {
    VectorXd eigv_h = P*eigv_H.col(i); // prolongate coarse eigenvector
    normalize_vec_metric(eigv_h, M_h);
    rqs_H[i] = eigv_h.dot(K_h*eigv_h);
  }
}

bool is_adjc_quad_stencil(const shared_ptr<quad_stencil> &x,
                          const shared_ptr<quad_stencil> &y,
                          mati_t &shared_verts) {
  mati_t elem_x = x->adjc_elem_;
  mati_t elem_y = y->adjc_elem_;

  find_intersections(elem_x, elem_y, shared_verts);

  if ( shared_verts.size() < 2 )
    return false;
  
  return true;
}

///=============== HEX STENCIL ===============///
hexs_stencil::hexs_stencil(const mati_t &adjc_elem, const vector<size_t> &faces)
    : adjc_elem_(adjc_elem) {
  adjc_vert_num_ = adjc_elem_.size();

  faces_.resize(faces.size(), 1);
  std::copy(faces.begin(), faces.end(), faces_.begin());

  if ( adjc_vert_num_ == 8 )
    fem_sf_jac_ = rest_sf_jac_ = hex8_shape_func_jac_;
  else if ( adjc_vert_num_ == 27 )
    fem_sf_jac_ = rest_sf_jac_ = hex27_shape_func_jac_;
  else
    ASSERT(0);

  this->dis_fem_sf_jac_ = nullptr;
  this->isotropic_      = true;
}

int hexs_stencil::set_full_discrete_shape_func(const matd_t &basisXYZ, const bool isotropic,
                                               const bool augmented) {
  if ( basisXYZ.size(1) != adjc_vert_num_ ) {
    ASSERT(0);
    return __LINE__;
  }

  this->dis_fem_sf_jac_ = full_hexs_dis_shape_func_jac;
  this->basisXYZ_  = basisXYZ;
  this->isotropic_ = isotropic;
  this->augmented_ = augmented;

  return 0;
}

int hexs_stencil::set_mtrs_and_qrs(const size_t per_face_to_many, const size_t qnum_per_cube,
                                   const mati_t &hexs_h, const matd_t &nods_h,
                                   const matd_t &lame_h) {
  const size_t num = qnum_per_cube*cbrt(faces_.size()*per_face_to_many); /*!!![ MODIFIED ] !!!*/
  const double *qr_ptn, *qr_wgts;
  switch ( num ) {
    case 1: qr_ptn = G_QR_PNT_1; qr_wgts = G_QR_WGT_1; break;
    case 2: qr_ptn = G_QR_PNT_2; qr_wgts = G_QR_WGT_2; break;
    case 4: qr_ptn = G_QR_PNT_4; qr_wgts = G_QR_WGT_4; break;
    case 6: qr_ptn = G_QR_PNT_6; qr_wgts = G_QR_WGT_6; break;
    case 8: qr_ptn = G_QR_PNT_8; qr_wgts = G_QR_WGT_8; break;
    default: ASSERT(0);
  }
  qr_num_ = num*num*num;
  qrs_ = zeros<double>(3, qr_num_);
  qrw_ = zeros<double>(1, qr_num_);
  for (size_t i = 0; i < num; ++i) {
    for (size_t j = 0; j < num; ++j) {
      for (size_t k = 0; k < num; ++k) {
        const size_t idx = i*num*num+j*num+k;
        qrs_(0, idx) = qr_ptn[i];
        qrs_(1, idx) = qr_ptn[j];
        qrs_(2, idx) = qr_ptn[k];
        qrw_[idx] = qr_wgts[i]*qr_wgts[j]*qr_wgts[k];
      }
    }
  }
  
  mati_t fine_ids(faces_.size()*per_face_to_many, 1);
  for (size_t i = 0; i < faces_.size(); ++i) {
    fine_ids(colon(i*per_face_to_many, (i+1)*per_face_to_many-1))
        = colon(faces_[i]*per_face_to_many, (faces_[i]+1)*per_face_to_many-1);
  }
  const mati_t patch_cell_h = hexs_h(colon(), fine_ids);
  const matd_t patch_matr_h = lame_h(colon(), fine_ids);

  const size_t v_stride = cbrt(adjc_elem_.size())-1;
  matd_t corner_vert = zeros<double>(nods_h.size(1), 8);
  for (size_t i = 0; i < 2; ++i) {
    for (size_t j = 0; j < 2; ++j) {
      for (size_t k = 0; k < 2; ++k) {
        const size_t idx = i*2*2+j*2+k;
        const size_t off = i*v_stride*(v_stride+1)*(v_stride+1)+j*v_stride*(v_stride+1)+k*v_stride;
        corner_vert(colon(), idx) = nods_h(colon(), adjc_elem_[off]);
      }
    }
  }

  ///-> paramterize all barycenters into [-1,+1]x[-1,+1]x[-1,+1]
  matd_t bc = zeros<double>(3, patch_cell_h.size(2)), param_bc = bc;
  for (size_t i = 0; i < patch_cell_h.size(2); ++i) {
    bc(colon(), i) = nods_h(colon(), patch_cell_h(colon(), i))*ones<double>(8, 1)/8;
    matd_t tmp_res = zeros<double>(3, 1);
    calc_inv_trilinear_interp_3d(corner_vert, bc(colon(), i), tmp_res, 10);
    param_bc(colon(), i) = tmp_res;
  }
      
  lame_on_qrs_ = zeros<double>(lame_h.size(1), qr_num_);
  for (size_t j = 0; j < qr_num_; ++j) {
    matd_t dist = zeros<double>(param_bc.size(2), 1);
    for (size_t k = 0; k < param_bc.size(2); ++k)
      dist[k] = norm(qrs_(colon(), j)-param_bc(colon(), k));
      
    const size_t idx = std::min_element(dist.begin(), dist.end())-dist.begin();
    lame_on_qrs_(colon(), j) = patch_matr_h(colon(), idx);
  }

  return 0;  
}

int hexs_stencil::set_def_grad_oper(const matd_t &adjc_nods) {
  adjc_rest_ = adjc_nods;
  {
    matd_t H = zeros<double>(adjc_vert_num_, 3);
    const double eps[3] = {0};
    if ( this->isotropic_ ) {
      //-> isotropic case
      if ( dis_fem_sf_jac_ )
        dis_fem_sf_jac_(&H[0], eps, basisXYZ_);
      else 
        fem_sf_jac_(&H[0], eps);
    } else {
      //-> anistropic case
      rest_sf_jac_(&H[0], eps);
    }
    matd_t Dm = adjc_rest_*H;
    inv(Dm);
    HinvDmH0_ = H*Dm;
  }
  
  if ( this->isotropic_ ) {
    adjc_H_invDmH_.resize(adjc_vert_num_, 3*qr_num_);
    dets_.resize(1, qr_num_);

    dNude_ = zeros<double>(adjc_vert_num_, 3*qr_num_);
    dNxde_ = zeros<double>(adjc_vert_num_, 3*qr_num_);
    DmHs_ = zeros<double>(3, 3*qr_num_);
    
    for (size_t i = 0; i < qr_num_; ++i) {
      matd_t H = zeros<double>(adjc_vert_num_, 3);
      if ( dis_fem_sf_jac_ )
        dis_fem_sf_jac_(&H[0], &qrs_(0, i), basisXYZ_);
      else
        fem_sf_jac_(&H[0], &qrs_(0, i));
      dNude_(colon(), colon(3*i, 3*i+2)) = H;

      matd_t H0 = zeros<double>(adjc_vert_num_, 3);
      rest_sf_jac_(&H0[0], &qrs_(0, i));
      matd_t DmH = adjc_nods*H0;
      matd_t cpDmH = DmH;
      dets_[i] = fabs(det(cpDmH));

      if ( inv(DmH) ) std::cerr << "# inv fail" << std::endl;
      DmHs_(colon(), colon(3*i, 3*i+2)) = DmH;
      
      adjc_H_invDmH_(colon(), colon(3*i, 3*i+2)) = H*DmH;
    }

    adjc_df_op_ = zeros<double>(3*adjc_vert_num_, 9*qr_num_);
    for (size_t i = 0; i < qr_num_; ++i) {
      adjc_df_op_(colon(), colon(9*i, 9*i+8)) = kroneckerId<double, 3>(adjc_H_invDmH_(colon(), colon(3*i, 3*i+2)));
    }
  } else if ( !this->isotropic_ && !this->augmented_ ) {
    ASSERT(dis_fem_sf_jac_);
    
    //-> anistropic bases
    const size_t b_cols = basisXYZ_.size(2)/3;

    const matd_t basisX = basisXYZ_(colon(), colon(0*b_cols, 1*b_cols-1));
    const matd_t basisY = basisXYZ_(colon(), colon(1*b_cols, 2*b_cols-1));
    const matd_t basisZ = basisXYZ_(colon(), colon(2*b_cols, 3*b_cols-1));

    adjc_df_op_ = zeros<double>(3*adjc_vert_num_, 9*qr_num_);
    dets_.resize(1, qr_num_);

    dNude_ = zeros<double>(3*adjc_vert_num_, 3*qr_num_);
    dNxde_ = zeros<double>(adjc_vert_num_, 3*qr_num_);
    DmHs_ = zeros<double>(3, 3*qr_num_);
        
    #pragma omp parallel for
    for (size_t i = 0; i < qr_num_; ++i) {
      matd_t Hx(adjc_vert_num_, 3), Hy(adjc_vert_num_, 3), Hz(adjc_vert_num_, 3);
      {
        dis_fem_sf_jac_(&Hx[0], &qrs_(0, i), basisX);
        dis_fem_sf_jac_(&Hy[0], &qrs_(0, i), basisY);
        dis_fem_sf_jac_(&Hz[0], &qrs_(0, i), basisZ);
      }
      
      matd_t H0 = zeros<double>(adjc_vert_num_, 3);
      rest_sf_jac_(&H0[0], &qrs_(0, i));
      dNxde_(colon(), colon(3*i, 3*i+2)) = H0;
      
      matd_t DmH = adjc_nods*H0;
      matd_t cpDmH = DmH;
      dets_[i] = fabs(det(cpDmH));

      if ( inv(DmH) ) std::cerr << "# inv fail" << std::endl;
      DmHs_(colon(), colon(3*i, 3*i+2)) = DmH;

      matd_t temp_df_op = zeros<double>(9, 3*adjc_vert_num_); {
        for (size_t p = 0; p < 3; ++p) {
          for (size_t q = 0; q < adjc_vert_num_; ++q) {
            temp_df_op(3*p+0, 3*q+0) = dot(Hx(q, colon()), DmH(colon(), p));
            temp_df_op(3*p+1, 3*q+1) = dot(Hy(q, colon()), DmH(colon(), p));
            temp_df_op(3*p+2, 3*q+2) = dot(Hz(q, colon()), DmH(colon(), p));
          }
        }
      }
      adjc_df_op_(colon(), colon(9*i, 9*i+8)) = trans(temp_df_op);

      for (size_t k = 0; k < adjc_vert_num_; ++k) {
        dNude_(3*k+0, colon(3*i, 3*i+2)) = Hx(k, colon());
        dNude_(3*k+1, colon(3*i, 3*i+2)) = Hy(k, colon());
        dNude_(3*k+2, colon(3*i, 3*i+2)) = Hz(k, colon());
      }
    }
  } else if ( !this->isotropic_ && this->augmented_ ) {
    ASSERT(dis_fem_sf_jac_);

    const size_t b_cols = basisXYZ_.size(2)/(3*3);
    const size_t b_num  = basisXYZ_.size(1);
    ASSERT(basisXYZ_.size(2)%9 == 0);

    //-> reshape the bases
    matd_t basis_entry = zeros<double>(9*b_num, b_cols);
    for (size_t p = 0; p < 3; ++p) {
      for (size_t q = 0; q < 3; ++q) {
        const size_t col_major_idx = 3*q+p;
        const size_t row_major_idx = 3*p+q;
        basis_entry(colon(col_major_idx*b_num, (col_major_idx+1)*b_num-1), colon())
            = basisXYZ_(colon(), colon(row_major_idx*b_cols, (row_major_idx+1)*b_cols-1));
      }
    }
    
    adjc_df_op_ = zeros<double>(3*adjc_vert_num_, 9*qr_num_);
    dets_.resize(1, qr_num_);

    //-> jacobian of bases for u and X
    dNude_ = zeros<double>(9*adjc_vert_num_, 3*qr_num_);
    dNxde_ = zeros<double>(adjc_vert_num_, 3*qr_num_);
    DmHs_ = zeros<double>(3, 3*qr_num_);

    for (size_t i = 0; i < qr_num_; ++i) {
      vector<matd_t> Hs(9, zeros<double>(adjc_vert_num_, 3));
      for (size_t k = 0; k < 9; ++k) {
        dis_fem_sf_jac_(&Hs[k](0, 0), &qrs_(0, i), basis_entry(colon(k*b_num, (k+1)*b_num-1), colon()));
      }
      
      matd_t H0 = zeros<double>(adjc_vert_num_, 3);
      rest_sf_jac_(&H0[0], &qrs_(0, i));
      dNxde_(colon(), colon(3*i, 3*i+2)) = H0;
      matd_t DmH = adjc_nods*H0;
      matd_t cpDmH = DmH;
      dets_[i] = fabs(det(cpDmH));

      if ( inv(DmH) ) std::cerr << "# inv fail" << std::endl;
      DmHs_(colon(), colon(3*i, 3*i+2)) = DmH;

      for (size_t k = 0; k < adjc_vert_num_; ++k) {
        for (size_t m = 0; m < 9; ++m) {
          dNude_(9*k+m, colon(3*i, 3*i+2)) = Hs[m](k, colon());
        }
      }
    }
  }
  return 0;  
}

int hexs_stencil::query_local_frame(const matd_t &curr_nods, const size_t qrID,
                                    matd_t &R) {
  matd_t H = zeros<double>(adjc_vert_num_, 3);
  const double eps[3] = {0};
  if ( this->isotropic_ ) {
    //-> isotropic case
    if ( dis_fem_sf_jac_ )
      dis_fem_sf_jac_(&H[0], eps, basisXYZ_);
    else 
      fem_sf_jac_(&H[0], eps);
  } else {
    //-> anistropic case
    rest_sf_jac_(&H[0], eps);
  }

  matd_t Dm = adjc_rest_*H;
  inv(Dm);
  R = (curr_nods*H)*Dm;
  hj::polar3d rs_decomp;
  rs_decomp(R);

  return 0;
}

bool is_adjc_quad_stencil(const shared_ptr<hexs_stencil> &x,
                          const shared_ptr<hexs_stencil> &y,
                          mati_t &shared_verts) {
  mati_t elem_x = x->adjc_elem_;
  mati_t elem_y = y->adjc_elem_;

  find_intersections(elem_x, elem_y, shared_verts);

  if ( shared_verts.size() < 4 )
    return false;
  
  return true;
}

extern "C" {
  /* x, x0 is for a single vert */
  void coro_ani_basis_df_3d_(double *val, const double *x, const double *x0, const double *R,
                             const double *Nu, const double *Nx, const double *dNude,
                             const double *dNxde, const double *Dm);
  void coro_ani_basis_df_3d_jac_(double *jac, const double *x, const double *x0, const double *R,
                                 const double *Nu, const double *Nx, const double *dNude,
                                 const double *dNxde, const double *Dm);
  void coro_iso_basis_df_3d_(double *val, const double *x, const double *x0, const double *R,
                             const double *Nu, const double *Nx, const double *dNude,
                             const double *dNxde, const double *Dm);
  void coro_iso_basis_df_3d_jac_(double *jac, const double *x, const double *x0, const double *R,
                                 const double *Nu, const double *Nx, const double *dNude,
                                 const double *dNxde, const double *Dm);
  void coro_aug_ani_df_3d_(double *val, const double *x, const double *x0, const double *R,
                           const double *Nu, const double *Nx, const double *dNude,
                           const double *dNxde, const double *Dm);
  void coro_aug_ani_df_3d_jac_(double *jac, const double *x, const double *x0, const double *R,
                               const double *Nu, const double *Nx, const double *dNude,
                               const double *dNxde, const double *Dm);
}

int hexs_stencil::get_coro_bases_df(const matd_t &curr_nods, const matd_t &R,
                                    const size_t qrID, matd_t &dF, matd_t *G) {
  assert(qrID >= 0 && qrID < qr_num_);
  
  dF = zeros<double>(3, 3);
  if ( G )
    *G = zeros<double>(9, 3*adjc_vert_num_);

  matd_t tmpF = zeros<double>(3, 3);
  //  matd_t tmpG = zeros<double>(9, 3);

  //-> summed over each vert
  if ( !this->isotropic_ && !this->augmented_ ) {
    matd_t dNude = zeros<double>(3, 3);
    matd_t dNxde = zeros<double>(1, 3);
    for (size_t i = 0; i < adjc_vert_num_; ++i) {
      dNude = dNude_(colon(3*i+0, 3*i+2), colon(3*qrID, 3*qrID+2));
      dNxde = dNxde_(i, colon(3*qrID, 3*qrID+2));

      coro_ani_basis_df_3d_(&tmpF[0], &curr_nods(0, i), &adjc_rest_(0, i), &R[0],
                            nullptr, nullptr, &dNude[0], &dNxde[0], &DmHs_(0, 3*qrID));      
      dF += tmpF;

      if ( G )
        coro_ani_basis_df_3d_jac_(&(*G)(0, 3*i), nullptr, nullptr, &R[0],
                                  nullptr, nullptr, &dNude[0], &dNxde[0], &DmHs_(0, 3*qrID));
      //      G(colon(), colon(3*i+0, 3*i+2)) = tmpG;
    }
  } else if ( !this->isotropic_ && this->augmented_ ) {
    double pNude[27], pNxde[3];
    itr_matrix<double*> dNude(9, 3, pNude);
    itr_matrix<double*> dNxde(1, 3, pNxde);

    for (size_t i = 0; i < adjc_vert_num_; ++i) {
      dNude = dNude_(colon(9*i+0, 9*i+8), colon(3*qrID, 3*qrID+2));
      dNxde = dNxde_(i, colon(3*qrID, 3*qrID+2));

      coro_aug_ani_df_3d_(&tmpF[0], &curr_nods(0, i), &adjc_rest_(0, i), &R[0],
                          nullptr, nullptr, &dNude[0], &dNxde[0], &DmHs_(0, 3*qrID));
      dF += tmpF;

      if ( G )
        coro_aug_ani_df_3d_jac_(&(*G)(0, 3*i), nullptr, nullptr, &R[0],
                                nullptr, nullptr, &dNude[0], &dNxde[0], &DmHs_(0, 3*qrID));
      //      G(colon(), colon(3*i+0, 3*i+2)) = tmpG;
    }
  } else if ( this->isotropic_ ) {
    matd_t dNude = zeros<double>(1, 3);
    matd_t dNxde = zeros<double>(1, 3);
    for (size_t i = 0; i < adjc_vert_num_; ++i) {
      dNude = dNude_(i, colon(3*qrID, 3*qrID+2));
      dNxde = dNxde_(i, colon(3*qrID, 3*qrID+2));

      coro_iso_basis_df_3d_(&tmpF[0], &curr_nods(0, i), &adjc_rest_(0, i), &R[0],
                            nullptr, nullptr, &dNude[0], &dNxde[0], &DmHs_(0, 3*qrID));
      dF += tmpF;

      if ( G )
        coro_iso_basis_df_3d_jac_(&(*G)(0, 3*i), nullptr, &adjc_rest_(0, i), &R[0],
                                  nullptr, nullptr, &dNude[0], &dNxde[0], &DmHs_(0, 3*qrID));
      //      G(colon(), colon(3*i+0, 3*i+2)) = tmpG;
    }
  }
  
  //-> dF=du+Id
  dF += eye<double>(3);
  return 0;
}

void hexs_stencil::apply_local_LMA(const mati_t &hexs_H, const matd_t &nods_H, const size_t eignum,
                                   const SparseMatrix<double> &K_h,
                                   const SparseMatrix<double> &M_h,
                                   const SparseMatrix<double> &P,
                                   const boost::property_tree::ptree &pt, 
                                   VectorXd &eigs_H,
                                   VectorXd &rqs_H) {
  //-> should be called after the energy is established
  ASSERT(adjc_rest_.size() > 0 && lame_on_qrs_.size() > 0);

  const string cons_law = pt.get<string>("elas.value");
  const size_t adjc_dim = 3*adjc_vert_num_;

  //-> for each quadrature points
  matd_t adjc_H = zeros<double>(adjc_dim, adjc_dim), HF = zeros<double>(9, 9);
  const matd_t R = eye<double>(3), dR = eye<double>(3);
  for (size_t j = 0; j < this->qr_num_; ++j) {
    matd_t G, dF = zeros<double>(3, 3);
    this->get_coro_bases_df(adjc_rest_, R, j, dF, &G);
            
    HF = zeros<double>(9, 9);
    if ( cons_law == "coro" )
      hex_coro_f_at_qr_hes_(&HF[0], &dF[0], &dR[0], &lame_on_qrs_(0, j), &lame_on_qrs_(1, j));
    else if ( cons_law == "stvk" )
      hex_stvk_f_at_qr_hes_(&HF[0], &dF[0], &lame_on_qrs_(0, j), &lame_on_qrs_(1, j));
    else if ( cons_law == "neohookean" )
      hex_neo_f_at_qr_hes_(&HF[0], &dF[0], &lame_on_qrs_(0, j), &lame_on_qrs_(1, j));
    else if ( cons_law == "stabneo" )
      hex_sta_neo_f_at_qr_hes_(&HF[0], &dF[0], &lame_on_qrs_(0, j), &lame_on_qrs_(1, j));
    else if ( cons_law == "bowerneo" )
      hex_bower_neo_f_at_qr_hes_(&HF[0], &dF[0], &lame_on_qrs_(0, j), &lame_on_qrs_(1, j));
      
    adjc_H += dets_[j]*qrw_[j]*trans(G)*HF*G;
  }

  //-> insert into a sparse K
  SparseMatrix<double> K_H = Map<const MatrixXd>(&adjc_H[0], adjc_H.size(1), adjc_H.size(2)).sparseView();

  //-> global to local mapping
  mati_t g2l = -1*ones<size_t>(nods_H.size(2), 1);
  g2l(adjc_elem_) = colon(0, adjc_elem_.size()-1);

  //-> get mass matrix
  SparseMatrix<double> M_H(adjc_H.size(1), adjc_H.size(2));
  M_H.setZero();
  for (size_t i = 0; i < faces_.size(); ++i) {
    const size_t fid = faces_[i];
    const double len = norm(nods_H(colon(), hexs_H(0, fid))-nods_H(colon(), hexs_H(1, fid)));
    const double area = len*len*len;
    for (size_t j = 0; j < hexs_H.size(1); ++j) {
      const size_t pid = hexs_H(j, fid);
      ASSERT(g2l[pid] != -1);
      M_H.coeffRef(3*g2l[pid]+0, 3*g2l[pid]+0) += area/8;
      M_H.coeffRef(3*g2l[pid]+1, 3*g2l[pid]+1) += area/8;
      M_H.coeffRef(3*g2l[pid]+2, 3*g2l[pid]+2) += area/8;
    }
  }

  //-> solve generalized eigen
  const int smallest = 0;
  MatrixXd eigv_H;
  solve_gen_eig_prob(K_H, M_H, eignum, smallest, "arpaca", &eigs_H, &eigv_H);

  //-> calc rayleigh quotient
  rqs_H = VectorXd::Zero(eigv_H.cols());
  for (size_t i = 0; i < eigv_H.cols(); ++i) {
    VectorXd eigv_h = P*eigv_H.col(i); // prolongate coarse eigenvector
    rqs_H[i] = eigv_h.dot(K_h*eigv_h)/eigv_h.dot(M_h*eigv_h);
  }
}
  
///===== stencil 2d energy =====///
fem_stencil_energy_d2::fem_stencil_energy_d2(const vector<shared_ptr<quad_stencil>> &elem_H,
                                             const mati_t &quad_H, const matd_t &nods_H,
                                             const Material type,
                                             const mati_t &quad_h, const matd_t &nods_h,
                                             const matd_t &face_lame_h, const double w)
    : elem_(elem_H), type_(type), w_(w), dim_(nods_H.size()) {
  const size_t one_to_many = quad_h.size(2)/quad_H.size(2);
  ASSERT(quad_h.size(2)%quad_H.size(2) == 0);

  //  #pragma omp parallel for
  for (size_t i = 0; i < elem_.size(); ++i) {
    elem_[i]->set_mtrs_and_qrs(one_to_many, quad_h, nods_h, face_lame_h);
    const matd_t adjc_nods = nods_H(colon(), elem_[i]->adjc_elem_);
    elem_[i]->set_def_grad_oper(adjc_nods);
  }
}

size_t fem_stencil_energy_d2::Nx() const {
  return dim_;
}

int fem_stencil_energy_d2::Val(const double *x, double *val) const {
  assert(x && val);
  
  RETURN_WITH_COND_TRUE(w_ == 0.0);
  itr_matrix<const double *> X(2, dim_/2, x);

  #pragma omp parallel for
  for (size_t i = 0; i < elem_.size(); ++i) {
    const auto &ptr = elem_[i];

    hj::polar2d rs_decomp;
    matd_t dR = zeros<double>(2, 2);
    
    const matd_t adjc_vert = X(colon(), ptr->adjc_elem_);
    const matd_t adjc_disp = adjc_vert-ptr->adjc_rest_;

    matd_t dF = zeros<double>(2, 2);  
    matd_t R = zeros<double>(2, 2);   // element wise rotation
    ptr->query_local_frame(adjc_vert, 0, R);

    double value = 0;
    for (size_t j = 0; j < ptr->qr_num_; ++j) {
      if ( ptr->isotropic_ ) {
        dF(colon()) = trans(ptr->adjc_df_op_(colon(), colon(4*j, 4*j+3)))*adjc_disp(colon())
            +eye<double>(2)(colon());      
      } else {
        matd_t G;
        ptr->get_coro_bases_df(adjc_vert, R, j, dF, G);
        ptr->adjc_df_op_(colon(), colon(4*j, 4*j+3)) = trans(G);
      }

      double vr = 0;
      switch ( type_ ) {
        case LINEAR:
          break;
        case COROTATED:
          dR = dF;
          rs_decomp(dR);          
          quad_coro_f_val_at_qr_(&vr, &dF[0], &dR[0], &ptr->lame_on_qrs_(0, j), &ptr->lame_on_qrs_(1, j));
          break;
        case STVK:
          quad_stvk_f_val_at_qr_(&vr, &dF[0], &ptr->lame_on_qrs_(0, j), &ptr->lame_on_qrs_(1, j));
          break;
        case NEOHOOKEAN:
          quad_neo_f_val_at_qr_(&vr, &dF[0], &ptr->lame_on_qrs_(0, j), &ptr->lame_on_qrs_(1, j));
          break;
        default:
          break;
      }

      value += vr*ptr->dets_[j]*ptr->qrw_[j];
    }

    #pragma omp critical
    {
      *val += w_*value;
    } 
  }
  return 0;
}

int fem_stencil_energy_d2::Gra(const double *x, double *gra) const {
  assert(x && gra);

  RETURN_WITH_COND_TRUE(w_ == 0.0);
  itr_matrix<const double *> X(2, dim_/2, x);
  itr_matrix<double *> G(2, dim_/2, gra);

  #pragma omp parallel for
  for (size_t i = 0; i < elem_.size(); ++i) {
    const auto &ptr = elem_[i];

    hj::polar2d rs_decomp;
    matd_t dR = zeros<double>(2, 2);
    
    const matd_t adjc_vert = X(colon(), ptr->adjc_elem_);
    const matd_t adjc_disp = adjc_vert-ptr->adjc_rest_;
      
    matd_t dF = zeros<double>(2, 2);
    matd_t R = zeros<double>(2, 2);
    ptr->query_local_frame(adjc_vert, 0, R);

    matd_t adjc_g = zeros<double>(ptr->adjc_vert_num_*2, 1), gF = zeros<double>(4, 1);
    for (size_t j = 0; j < ptr->qr_num_; ++j) {
      if ( ptr->isotropic_ ) {
        dF(colon()) = trans(ptr->adjc_df_op_(colon(), colon(4*j, 4*j+3)))*adjc_disp(colon())
            +eye<double>(2)(colon());
      } else {
        matd_t G;
        ptr->get_coro_bases_df(adjc_vert, R, j, dF, G);
        ptr->adjc_df_op_(colon(), colon(4*j, 4*j+3)) = trans(G);
      }
            
      gF = zeros<double>(4, 1);
      switch ( type_ ) {
        case LINEAR:
          break;
        case COROTATED:
          dR = dF;
          rs_decomp(dR);
          quad_coro_f_val_at_qr_jac_(&gF[0], &dF[0], &dR[0], &ptr->lame_on_qrs_(0, j), &ptr->lame_on_qrs_(1, j));
          break;
        case STVK:
          quad_stvk_f_val_at_qr_jac_(&gF[0], &dF[0], &ptr->lame_on_qrs_(0, j), &ptr->lame_on_qrs_(1, j));
          break;
        case NEOHOOKEAN:
          quad_neo_f_val_at_qr_jac_(&gF[0], &dF[0], &ptr->lame_on_qrs_(0, j), &ptr->lame_on_qrs_(1, j));
          break;
        default:
          break;
      }

      adjc_g += ptr->dets_[j]*ptr->qrw_[j]*ptr->adjc_df_op_(colon(), colon(4*j, 4*j+3))*gF;
    }

    #pragma omp critical
    {
      G(colon(), ptr->adjc_elem_) += w_*itr_matrix<const double *>(2, ptr->adjc_vert_num_, &adjc_g[0]);
    }
  }
  return 0;
}

int fem_stencil_energy_d2::Hes(const double *x, vector<Triplet<double>> *hes) const {
  assert(x && hes);

  RETURN_WITH_COND_TRUE(w_ == 0.0);
  itr_matrix<const double *> X(2, dim_/2, x);

  #pragma omp parallel for
  for (size_t i = 0; i < elem_.size(); ++i) {
    const auto &ptr = elem_[i];

    hj::polar2d rs_decomp;
    matd_t dR = zeros<double>(2, 2);
    const matd_t adjc_vert = X(colon(), ptr->adjc_elem_);
    const matd_t adjc_disp = adjc_vert-ptr->adjc_rest_;

    matd_t dF = zeros<double>(2, 2);
    matd_t R = zeros<double>(2, 2);
    ptr->query_local_frame(adjc_vert, 0, R);  

    const size_t adjc_dim = ptr->adjc_vert_num_*2;
    matd_t adjc_H = zeros<double>(adjc_dim, adjc_dim), HF = zeros<double>(4, 4);
    for (size_t j = 0; j < ptr->qr_num_; ++j) {
      if ( ptr->isotropic_ ) {
        dF(colon()) = trans(ptr->adjc_df_op_(colon(), colon(4*j, 4*j+3)))*adjc_disp(colon())
            +eye<double>(2)(colon());
      } else {
        matd_t G;
        ptr->get_coro_bases_df(adjc_vert, R, j, dF, G);
        ptr->adjc_df_op_(colon(), colon(4*j, 4*j+3)) = trans(G);
      }
            
      HF = zeros<double>(4, 4);
      switch ( type_ ) {
        case LINEAR:
          break;
        case COROTATED:
          dR = dF;
          rs_decomp(dR);
          quad_coro_f_val_at_qr_hes_(&HF[0], &dF[0], &dR[0], &ptr->lame_on_qrs_(0, j), &ptr->lame_on_qrs_(1, j));
          break;
        case STVK:
          quad_stvk_f_val_at_qr_hes_(&HF[0], &dF[0], &ptr->lame_on_qrs_(0, j), &ptr->lame_on_qrs_(1, j));
          break;
        case NEOHOOKEAN:          
          quad_neo_f_val_at_qr_hes_(&HF[0], &dF[0], &ptr->lame_on_qrs_(0, j), &ptr->lame_on_qrs_(1, j));
          break;
        default:
          break;
      }
      
      adjc_H += ptr->dets_[j]*ptr->qrw_[j]*ptr->adjc_df_op_(colon(), colon(4*j, 4*j+3))*HF*trans(ptr->adjc_df_op_(colon(), colon(4*j, 4*j+3)));      
    }

    #pragma omp critical
    {
      for (size_t p = 0; p < adjc_H.size(1); ++p) {
        for (size_t q = 0; q < adjc_H.size(2); ++q) {
          const size_t I = 2*ptr->adjc_elem_[p/2]+p%2;
          const size_t J = 2*ptr->adjc_elem_[q/2]+q%2;
          hes->push_back(Triplet<double>(I, J, w_*adjc_H(p, q)));
        }
      }
    }
  }
  return 0;
}

///===== stencil 3d energy =====///
fem_stencil_energy_d3::fem_stencil_energy_d3(const vector<shared_ptr<hexs_stencil>> &elem_H,
                                             const mati_t &hexs_H, const matd_t &nods_H,
                                             const Material type,
                                             const mati_t &hexs_h, const matd_t &nods_h,
                                             const matd_t &face_lame_h,
                                             const size_t qnum_per_cube,
                                             const double w)
    : w_(w), dim_(nods_H.size()), elem_(elem_H), type_(type) {
  const size_t one_to_many = hexs_h.size(2)/hexs_H.size(2);
  ASSERT(hexs_h.size(2)%hexs_H.size(2) == 0);

  #pragma omp parallel for
  for (size_t i = 0; i < elem_.size(); ++i) {
    elem_[i]->set_mtrs_and_qrs(one_to_many, qnum_per_cube, hexs_h, nods_h, face_lame_h);
    const matd_t adjc_nods = nods_H(colon(), elem_[i]->adjc_elem_);
    elem_[i]->set_def_grad_oper(adjc_nods);
  }
}

size_t fem_stencil_energy_d3::Nx() const {
  return dim_;
}

int fem_stencil_energy_d3::Val(const double *x, double *val) const {
  assert(x && val);
  
  RETURN_WITH_COND_TRUE(w_ == 0.0);
  itr_matrix<const double *> X(3, dim_/3, x);

  double &val_cp = *val;

  #pragma omp parallel for reduction(+:val_cp)
  for (size_t i = 0; i < elem_.size(); ++i) {
    const auto &ptr = elem_[i];
  
    const matd_t adjc_vert = X(colon(), ptr->adjc_elem_);
    const matd_t adjc_disp = adjc_vert-ptr->adjc_rest_;

    hj::polar3d rs_decomp;      
    matd_t dR = zeros<double>(3, 3);
    matd_t dF = zeros<double>(3, 3);
    matd_t R = zeros<double>(3, 3);

    if ( !ptr->isotropic_ ) {
      R = adjc_vert*ptr->HinvDmH0_;
      rs_decomp(R);
    }
    
    double value = 0;
    switch ( type_ ) {
      case LINEAR: {
        break;
      }
      case COROTATED: {
        for (size_t j = 0; j < ptr->qr_num_; ++j) {
          if ( ptr->isotropic_ ) {
            dF = adjc_disp*ptr->adjc_H_invDmH_(colon(), colon(3*j, 3*j+2))+eye<double>(3);
          } else {
            ptr->get_coro_bases_df(adjc_vert, R, j, dF, NULL);
          }
          double vr = 0;
          dR = dF;
          rs_decomp(dR);
          hex_coro_f_at_qr_(&vr, &dF[0], &dR[0], &ptr->lame_on_qrs_(0, j), &ptr->lame_on_qrs_(1, j));
          value += vr*ptr->dets_[j]*ptr->qrw_[j];
        }
        break;
      }
      case STVK: {
        for (size_t j = 0; j < ptr->qr_num_; ++j) {
          if ( ptr->isotropic_ ) {
            dF = adjc_disp*ptr->adjc_H_invDmH_(colon(), colon(3*j, 3*j+2))+eye<double>(3);
          } else {
            ptr->get_coro_bases_df(adjc_vert, R, j, dF, NULL);
          }
          double vr = 0;
          hex_stvk_f_at_qr_(&vr, &dF[0], &ptr->lame_on_qrs_(0, j), &ptr->lame_on_qrs_(1, j));
          value += vr*ptr->dets_[j]*ptr->qrw_[j];
        }
        break;
      }
      case NEOHOOKEAN: {
        for (size_t j = 0; j < ptr->qr_num_; ++j) {
          if ( ptr->isotropic_ ) {
            dF = adjc_disp*ptr->adjc_H_invDmH_(colon(), colon(3*j, 3*j+2))+eye<double>(3);
          } else {
            ptr->get_coro_bases_df(adjc_vert, R, j, dF, NULL);
          }
          double vr = 0;
          hex_neo_f_at_qr_(&vr, &dF[0], &ptr->lame_on_qrs_(0, j), &ptr->lame_on_qrs_(1, j));
          value += vr*ptr->dets_[j]*ptr->qrw_[j];
        }
        break;
      }
      case STANEO: {
        for (size_t j = 0; j < ptr->qr_num_; ++j) {
          if ( ptr->isotropic_ ) {
            dF = adjc_disp*ptr->adjc_H_invDmH_(colon(), colon(3*j, 3*j+2))+eye<double>(3);
          } else {
            ptr->get_coro_bases_df(adjc_vert, R, j, dF, NULL);
          }
          double vr = 0;
          hex_sta_neo_f_at_qr_(&vr, &dF[0], &ptr->lame_on_qrs_(0, j), &ptr->lame_on_qrs_(1, j));
          value += vr*ptr->dets_[j]*ptr->qrw_[j];
        }
        break;
      }
      case BOWERNEO: {
        for (size_t j = 0; j < ptr->qr_num_; ++j) {
          if ( ptr->isotropic_ ) {
            dF = adjc_disp*ptr->adjc_H_invDmH_(colon(), colon(3*j, 3*j+2))+eye<double>(3);
          } else {
            ptr->get_coro_bases_df(adjc_vert, R, j, dF, NULL);
          }
          double vr = 0;
          hex_bower_neo_f_at_qr_(&vr, &dF[0], &ptr->lame_on_qrs_(0, j), &ptr->lame_on_qrs_(1, j));
          value += vr*ptr->dets_[j]*ptr->qrw_[j];
        }
        break;
      }
      default:
        break;
    }

    val_cp = val_cp+w_*value;

  }

  return 0;
}

int fem_stencil_energy_d3::Gra(const double *x, double *gra) const {
  assert(x && gra);

  RETURN_WITH_COND_TRUE(w_ == 0.0);
  itr_matrix<const double *> X(3, dim_/3, x);
  itr_matrix<double *> G(3, dim_/3, gra);

  #pragma omp parallel for
  for (size_t i = 0; i < elem_.size(); ++i) {
    const auto &ptr = elem_[i];

    const matd_t adjc_vert = X(colon(), ptr->adjc_elem_);
    const matd_t adjc_disp = adjc_vert-ptr->adjc_rest_;

    hj::polar3d rs_decomp;
    matd_t dR = zeros<double>(3, 3);
    matd_t dF = zeros<double>(3, 3);
    matd_t R = zeros<double>(3, 3);

    matd_t GF;
    matd_t adjc_g, gF;

    if ( !ptr->isotropic_ ) {
      R = adjc_vert*ptr->HinvDmH0_;
      rs_decomp(R);
    }

    adjc_g = zeros<double>(ptr->adjc_vert_num_*3, 1);
    gF = zeros<double>(9, 1);
    for (size_t j = 0; j < ptr->qr_num_; ++j) {
      if ( ptr->isotropic_ ) {
        dF = adjc_disp*ptr->adjc_H_invDmH_(colon(), colon(3*j, 3*j+2))+eye<double>(3);
      } else {
        ptr->get_coro_bases_df(adjc_vert, R, j, dF, &GF);
        ptr->adjc_df_op_(colon(), colon(9*j, 9*j+8)) = trans(GF);
      }
            
      gF = zeros<double>(9, 1);
      switch ( type_ ) {
        case LINEAR:
          break;
        case COROTATED:
          dR = dF;
          rs_decomp(dR);
          hex_coro_f_at_qr_jac_(&gF[0], &dF[0], &dR[0], &ptr->lame_on_qrs_(0, j), &ptr->lame_on_qrs_(1, j));
          break;
        case STVK:
          hex_stvk_f_at_qr_jac_(&gF[0], &dF[0], &ptr->lame_on_qrs_(0, j), &ptr->lame_on_qrs_(1, j));
          break;
        case NEOHOOKEAN:
          hex_neo_f_at_qr_jac_(&gF[0], &dF[0], &ptr->lame_on_qrs_(0, j), &ptr->lame_on_qrs_(1, j));
          break;
        case STANEO:
          hex_sta_neo_f_at_qr_jac_(&gF[0], &dF[0], &ptr->lame_on_qrs_(0, j), &ptr->lame_on_qrs_(1, j));
          break;
        case BOWERNEO:
          hex_bower_neo_f_at_qr_jac_(&gF[0], &dF[0], &ptr->lame_on_qrs_(0, j), &ptr->lame_on_qrs_(1, j));
          break;
        default:
          break;
      }

      adjc_g += ptr->dets_[j]*ptr->qrw_[j]*zjucad::matrix::temp(ptr->adjc_df_op_(colon(), colon(9*j, 9*j+8))*gF); 
    }
    #pragma omp critical
    {
      G(colon(), ptr->adjc_elem_) += w_*itr_matrix<const double *>(3, ptr->adjc_vert_num_, &adjc_g[0]); 
    }
  }
  return 0;
}

int fem_stencil_energy_d3::Hes(const double *x, vector<Triplet<double>> *hes) const {
  assert(x && hes);

  RETURN_WITH_COND_TRUE(w_ == 0.0);
  itr_matrix<const double *> X(3, dim_/3, x);
  
  #if 0
  cout << "cube size: " << elem_.size() << endl;
  cout << "cueb quad num: " << elem_[0]->qr_num_ << endl;
  std::chrono::time_point<std::chrono::system_clock> tic, toc;
  tic = std::chrono::system_clock::now();
  #endif
  
  #pragma omp parallel for
  for (size_t i = 0; i < elem_.size(); ++i) {
    const auto &ptr = elem_[i];
    const matd_t adjc_vert = X(colon(), ptr->adjc_elem_);
    const matd_t adjc_disp = adjc_vert-ptr->adjc_rest_;

    hj::polar3d rs_decomp;
    matd_t dR = zeros<double>(3, 3);    
    matd_t dF = zeros<double>(3, 3);
    matd_t R = zeros<double>(3, 3);

    matd_t adjc_H, HF, GF;

    if ( !ptr->isotropic_ ) {
      R = adjc_vert*ptr->HinvDmH0_;
      rs_decomp(R);
    }

    const size_t adjc_dim = ptr->adjc_vert_num_*3;
    adjc_H = zeros<double>(adjc_dim, adjc_dim);
    HF = zeros<double>(9, 9);
    for (size_t j = 0; j < ptr->qr_num_; ++j) {
      if ( ptr->isotropic_ ) {
        dF = adjc_disp*ptr->adjc_H_invDmH_(colon(), colon(3*j, 3*j+2))+eye<double>(3);
      } else {
        ptr->get_coro_bases_df(adjc_vert, R, j, dF, &GF);
        ptr->adjc_df_op_(colon(), colon(9*j, 9*j+8)) = trans(GF);
      }        
            
      HF = zeros<double>(9, 9);
      switch ( type_ ) {
        case LINEAR:
          break;
        case COROTATED:
          dR = dF;
          rs_decomp(dR);
          hex_coro_f_at_qr_hes_(&HF[0], &dF[0], &dR[0], &ptr->lame_on_qrs_(0, j), &ptr->lame_on_qrs_(1, j));
          break;
        case STVK:
          hex_stvk_f_at_qr_hes_(&HF[0], &dF[0], &ptr->lame_on_qrs_(0, j), &ptr->lame_on_qrs_(1, j));
          break;
        case NEOHOOKEAN:          
          hex_neo_f_at_qr_hes_(&HF[0], &dF[0], &ptr->lame_on_qrs_(0, j), &ptr->lame_on_qrs_(1, j));
          break;
        case STANEO:
          hex_sta_neo_f_at_qr_hes_(&HF[0], &dF[0], &ptr->lame_on_qrs_(0, j), &ptr->lame_on_qrs_(1, j));
          break;
        case BOWERNEO:
          hex_bower_neo_f_at_qr_hes_(&HF[0], &dF[0], &ptr->lame_on_qrs_(0, j), &ptr->lame_on_qrs_(1, j));
          break;
        default:
          break;
      }

      adjc_H += ptr->dets_[j]*ptr->qrw_[j]*zjucad::matrix::temp(ptr->adjc_df_op_(colon(), colon(9*j, 9*j+8))*HF)*trans(ptr->adjc_df_op_(colon(), colon(9*j, 9*j+8)));
    }

    vector<Triplet<double>> local_hes;
    local_hes.reserve(adjc_H.size());
    for (size_t p = 0; p < adjc_H.size(1); ++p) {
      for (size_t q = 0; q < adjc_H.size(2); ++q) {
        const size_t I = 3*ptr->adjc_elem_[p/3]+p%3;
        const size_t J = 3*ptr->adjc_elem_[q/3]+q%3;
        local_hes.push_back(Triplet<double>(I, J, w_*adjc_H(p, q)));
      }
    }
          
    #pragma omp critical
    {
      hes->insert(hes->end(), local_hes.begin(), local_hes.end());
    }
  }

  #if 0
  toc = std::chrono::system_clock::now();
  std::chrono::duration<double, std::milli> time_diff = toc-tic;
  cout << "# time diff: " << time_diff.count() << endl;
  getchar();
  #endif
  return 0;
}

///===== build stencil energies =====///
shared_ptr<Functional<double>>
build_stencil_energy(const string &mtr, const vector<shared_ptr<quad_stencil>> &stencils,
                     const mati_t &quad_H, const matd_t &nods_H,
                     const mati_t &quad_h, const matd_t &nods_h, const matd_t &face_lame_h) {
  if ( mtr == "coro" )
    return make_shared<fem_stencil_energy_d2>(stencils, quad_H, nods_H, fem_stencil_energy_d2::COROTATED,
                                              quad_h, nods_h, face_lame_h, 1.0);
  if ( mtr == "neohookean" )
    return make_shared<fem_stencil_energy_d2>(stencils, quad_H, nods_H, fem_stencil_energy_d2::NEOHOOKEAN,
                                              quad_h, nods_h, face_lame_h, 1.0);
  if ( mtr == "linear" )
    return make_shared<fem_stencil_energy_d2>(stencils, quad_H, nods_H, fem_stencil_energy_d2::LINEAR,
                                              quad_h, nods_h, face_lame_h, 1.0);
  if ( mtr == "stvk" )
    return make_shared<fem_stencil_energy_d2>(stencils, quad_H, nods_H, fem_stencil_energy_d2::STVK,
                                              quad_h, nods_h, face_lame_h, 1.0);
}

shared_ptr<Functional<double>>
build_stencil_energy(const string &mtr, const vector<shared_ptr<hexs_stencil>> &stencils,
                     const size_t qnum_per_cube,
                     const mati_t &hexs_H, const matd_t &nods_H,
                     const mati_t &hexs_h, const matd_t &nods_h, const matd_t &face_lame_h) {
  if ( mtr == "coro" )
    return make_shared<fem_stencil_energy_d3>(stencils, hexs_H, nods_H, fem_stencil_energy_d3::COROTATED,
                                              hexs_h, nods_h, face_lame_h, qnum_per_cube, 1.0);
  if ( mtr == "neohookean" )
    return make_shared<fem_stencil_energy_d3>(stencils, hexs_H, nods_H, fem_stencil_energy_d3::NEOHOOKEAN,
                                              hexs_h, nods_h, face_lame_h, qnum_per_cube, 1.0);
  if ( mtr == "linear" )
    return make_shared<fem_stencil_energy_d3>(stencils, hexs_H, nods_H, fem_stencil_energy_d3::LINEAR,
                                              hexs_h, nods_h, face_lame_h, qnum_per_cube, 1.0);
  if ( mtr == "stvk" )
    return make_shared<fem_stencil_energy_d3>(stencils, hexs_H, nods_H, fem_stencil_energy_d3::STVK,
                                              hexs_h, nods_h, face_lame_h, qnum_per_cube, 1.0);
  if ( mtr == "stabneo" )
    return make_shared<fem_stencil_energy_d3>(stencils, hexs_H, nods_H, fem_stencil_energy_d3::STANEO,
                                              hexs_h, nods_h, face_lame_h, qnum_per_cube, 1.0);
  if ( mtr == "bowerneo" )
    return make_shared<fem_stencil_energy_d3>(stencils, hexs_H, nods_H, fem_stencil_energy_d3::BOWERNEO,
                                              hexs_h, nods_h, face_lame_h, qnum_per_cube, 1.0);
    
}


}



///===== deprecated =====///
// int fem_quad9_stencil::set_mtrs_and_qrs(const size_t per_face_to_many,
//                                         const mati_t &quad_h, const matd_t &nods_h,
//                                         const matd_t &lame_h) {
//   const size_t num = 2*sqrt(faces_.size()*per_face_to_many);
//   const double *qr_ptn, *qr_wgts;
//   switch ( num ) {
//     case 2: qr_ptn = G_QR_PNT_2; qr_wgts = G_QR_WGT_2; break;
//     case 4: qr_ptn = G_QR_PNT_4; qr_wgts = G_QR_WGT_4; break;
//     case 6: qr_ptn = G_QR_PNT_6; qr_wgts = G_QR_WGT_6; break;
//     case 8: qr_ptn = G_QR_PNT_8; qr_wgts = G_QR_WGT_8; break;
//     default: ASSERT(0);
//   }
//   qr_num_ = num*num;
//   qrs_ = zeros<double>(2, qr_num_);
//   qrw_ = zeros<double>(1, qr_num_);
//   for (size_t i = 0; i < num; ++i) {
//     for (size_t j = 0; j < num; ++j) {
//       const size_t idx = i*num+j;
//       qrs_(0, idx) = qr_ptn[i];
//       qrs_(1, idx) = qr_ptn[j];
//       qrw_[idx] = qr_wgts[i]*qr_wgts[j];
//     }
//   }
  
//   mati_t fine_ids(faces_.size()*per_face_to_many, 1);
//   for (size_t i = 0; i < faces_.size(); ++i) {
//     fine_ids(colon(i*per_face_to_many, (i+1)*per_face_to_many-1))
//         = colon(faces_[i]*per_face_to_many, (faces_[i]+1)*per_face_to_many-1);
//   }
//   const mati_t patch_h     = quad_h(colon(), fine_ids);
//   const matd_t patch_mtr_h = lame_h(colon(), fine_ids);

//   const size_t vert_strip_size = sqrt(adjc_elem_.size());  
//   matd_t corner_vert = zeros<double>(nods_h.size(1), 4);
//   corner_vert(colon(), 0) = nods_h(colon(), adjc_elem_[0]);
//   corner_vert(colon(), 1) = nods_h(colon(), adjc_elem_[0+vert_strip_size-1]);
//   corner_vert(colon(), 2) = nods_h(colon(), adjc_elem_[adjc_elem_.size()-vert_strip_size]);
//   corner_vert(colon(), 3) = nods_h(colon(), adjc_elem_[adjc_elem_.size()-1]);

//   /// paramterize all BCs into [-1,+1]x[-1,+1]
//   matd_t bc = zeros<double>(2, patch_h.size(2)), param_bc = bc;
//   for (size_t i = 0; i < patch_h.size(2); ++i) {
//     bc(colon(), i) = nods_h(colon(), patch_h(colon(), i))*ones<double>(4, 1)/4;
//     matd_t tmp_res = zeros<double>(2, 1);
//     calc_inv_bilinear_interp_2d(corner_vert, bc(colon(), i), tmp_res, 10);
//     param_bc(colon(), i) = tmp_res;
//   }
      
//   lame_on_qrs_ = zeros<double>(lame_h.size(1), qr_num_);
//   for (size_t j = 0; j < qr_num_; ++j) {
//     matd_t dist = zeros<double>(param_bc.size(2), 1);
//     for (size_t k = 0; k < param_bc.size(2); ++k)
//       dist[k] = norm(qrs_(colon(), j)-param_bc(colon(), k));
      
//     const size_t idx = std::min_element(dist.begin(), dist.end())-dist.begin();
//     lame_on_qrs_(colon(), j) = patch_mtr_h(colon(), idx);
//   }

//   return 0;
// }


// void quad_dis_shape_func_val(double *val, const double *eps, const matd_t &basisX, const matd_t &basisY) {
//   const size_t basis_num_x = basisX.size(1), sample_num_x = basisX.size(2);
//   const size_t basis_num_y = basisY.size(1), sample_num_y = basisY.size(2);

//   const double spacing_x = 2.0/(sample_num_x-1), spacing_y = 2.0/(sample_num_y-1);
//   const matd_t sample_pos_x = -1*ones<double>(sample_num_x, 1)+spacing_x*matd_t(colon(0, sample_num_x-1));
//   const matd_t sample_pos_y = -1*ones<double>(sample_num_y, 1)+spacing_y*matd_t(colon(0, sample_num_y-1));

//   const size_t loc_x = std::floor((eps[0]+1)/spacing_x);
//   const size_t loc_y = std::floor((eps[1]+1)/spacing_y);
  
//   for (size_t i = 0; i < basis_num_x; ++i) {
//     for (size_t j = 0; j < basis_num_y; ++j) {
//       const size_t idx = i*basis_num_y+j;

//       const double sf_val_i = two_point_lag_interp(eps[0],
//                                                    sample_pos_x[loc_x],   basisX(i, loc_x),
//                                                    sample_pos_x[loc_x+1], basisX(i, loc_x+1));      
//       const double sf_val_j = two_point_lag_interp(eps[1],
//                                                    sample_pos_y[loc_y],   basisY(j, loc_y),
//                                                    sample_pos_y[loc_y+1], basisY(j, loc_y+1));
      
//       val[idx] = sf_val_i*sf_val_j;
//     }
//   }
// }

// void quad_dis_shape_func_jac(double *jac, const double *eps, const matd_t &basisX, const matd_t &basisY) {
//   const size_t basis_num_x = basisX.size(1), sample_num_x = basisX.size(2);
//   const size_t basis_num_y = basisY.size(1), sample_num_y = basisY.size(2);

//   itr_matrix<double *> Jac(basis_num_x*basis_num_y, 2, jac);

//   const double spacing_x = 2.0/(sample_num_x-1), spacing_y = 2.0/(sample_num_y-1);
//   const matd_t sample_pos_x = -1*ones<double>(sample_num_x, 1)+spacing_x*matd_t(colon(0, sample_num_x-1));
//   const matd_t sample_pos_y = -1*ones<double>(sample_num_y, 1)+spacing_y*matd_t(colon(0, sample_num_y-1));

//   const size_t loc_x = std::floor((eps[0]+1)/spacing_x);
//   const size_t loc_y = std::floor((eps[1]+1)/spacing_y);
  
//   for (size_t i = 0; i < basis_num_x; ++i) {
//     for (size_t j = 0; j < basis_num_y; ++j) {
//       const size_t idx = i*basis_num_y+j;

//       const double sf_val_i = two_point_lag_interp(eps[0],
//                                                    sample_pos_x[loc_x],   basisX(i, loc_x),
//                                                    sample_pos_x[loc_x+1], basisX(i, loc_x+1));      
//       const double sf_val_j = two_point_lag_interp(eps[1],
//                                                    sample_pos_y[loc_y],   basisY(j, loc_y),
//                                                    sample_pos_y[loc_y+1], basisY(j, loc_y+1));

//       Jac(idx, 0) = (basisX(i, loc_x+1)-basisX(i, loc_x))/spacing_x*sf_val_j;
//       Jac(idx, 1) = sf_val_i*(basisY(j, loc_y+1)-basisY(j, loc_y))/spacing_y;
//     }
//   }
// }
  // void quad4_shape_function_jac_(double *jac, const double *eps, const double *px, const double *py);
  // void quad9_shape_function_jac_(double *jac, const double *eps, const double *px, const double *py);
  // void quad16_shape_function_jac_(double *jac, const double *eps, const double *px, const double *py);


//===============================================================================
// fem_quad9_stencil::fem_quad9_stencil(const mati_t &quad9, const vector<size_t> &faces)
//     : quad_stencil(0, quad9, quad9, faces, nullptr) {}

// int fem_quad9_stencil::set_def_grad_oper(const matd_t &adjc_nods) {
//   adjc_rest_ = adjc_nods;
  
//   //-> high order basis part
//   if ( this->isotropic_ ) {
//     adjc_H_invDmH_.resize(adjc_vert_num_, 2*qr_num_);
//     dets_.resize(1, qr_num_);
    
//     dNude_ = zeros<double>(adjc_vert_num_, 2*qr_num_);
//     dNxde_ = zeros<double>(adjc_vert_num_, 2*qr_num_);
//     DmHs_ = zeros<double>(2, 2*qr_num_);
    
//     for (size_t i = 0; i < qr_num_; ++i) {
//       matd_t H = zeros<double>(adjc_vert_num_, 2);
//       if ( full_dis_high_sf_jac_ )
//         full_dis_high_sf_jac_(&H[0], &qrs_(0, i), basisXY_);
//       // else if ( dis_high_order_sf_jac_ )
//       //   dis_high_order_sf_jac_(&H[0], &qrs_(0, i), basisX_, basisY_);
//       else 
//         high_order_sf_jac_(&H[0], &qrs_(0, i));

//       dNude_(colon(), colon(2*i, 2*i+1)) = H;

//       matd_t H0 = zeros<double>(adjc_vert_num_, 2);
//       rest_sf_jac_(&H0[0], &qrs_(0, i));
//       dNxde_(colon(), colon(2*i, 2*i+1)) = H0;
//       matd_t DmH = adjc_nods*H0;
//       matd_t cpDmH = DmH;
//       dets_[i] = fabs(det(cpDmH));

//       if ( inv(DmH) ) std::cerr << "# inv fail" << std::endl;
//       DmHs_(colon(), colon(2*i, 2*i+1)) = DmH;
      
//       adjc_H_invDmH_(colon(), colon(2*i, 2*i+1)) = H*DmH;
//     }

//     adjc_df_op_ = zeros<double>(2*adjc_vert_num_, 4*qr_num_);
//     for (size_t i = 0; i < qr_num_; ++i) {
//       adjc_df_op_(colon(), colon(4*i, 4*i+3)) = kroneckerId<double, 2>(adjc_H_invDmH_(colon(), colon(2*i, 2*i+1)));
//     }
//   } else {
//     //-> anistropic bases
//     const size_t b_cols = basisXY_.size(2)/2;

//     const matd_t basisX = basisXY_(colon(), colon(0*b_cols, 1*b_cols-1));
//     const matd_t basisY = basisXY_(colon(), colon(1*b_cols, 2*b_cols-1));

//     adjc_df_op_ = zeros<double>(2*adjc_vert_num_, 4*qr_num_);
//     dets_.resize(1, qr_num_);
    
//     dNude_ = zeros<double>(2*adjc_vert_num_, 2*qr_num_);
//     dNxde_ = zeros<double>(adjc_vert_num_, 2*qr_num_);
//     DmHs_ = zeros<double>(2, 2*qr_num_);
                         
//     #pragma omp parallel for
//     for (size_t i = 0; i < qr_num_; ++i) {
//       matd_t Hx(adjc_vert_num_, 2), Hy(adjc_vert_num_, 2);
//       {
//         full_dis_high_sf_jac_(&Hx[0], &qrs_(0, i), basisX);
//         full_dis_high_sf_jac_(&Hy[0], &qrs_(0, i), basisY);
//       }

//       matd_t H0 = zeros<double>(adjc_vert_num_, 2);
//       rest_sf_jac_(&H0[0], &qrs_(0, i));
//       dNxde_(colon(), colon(2*i, 2*i+1)) = H0;
//       matd_t DmH = adjc_nods*H0;
//       matd_t cpDmH = DmH;
//       dets_[i] = fabs(det(cpDmH));

//       if ( inv(DmH) ) std::cerr << "# inv fail" << std::endl;
//       DmHs_(colon(), colon(2*i, 2*i+1)) = DmH;

//       matd_t temp_df_op = zeros<double>(4, 2*adjc_vert_num_); {
//         for (size_t p = 0; p < 2; ++p) {
//           for (size_t q = 0; q < adjc_vert_num_; ++q) {
//             temp_df_op(2*p+0, 2*q+0) = dot(Hx(q, colon()), DmH(colon(), p));
//             temp_df_op(2*p+1, 2*q+1) = dot(Hy(q, colon()), DmH(colon(), p));
//           }
//         }
//       }
//       adjc_df_op_(colon(), colon(4*i, 4*i+3)) = trans(temp_df_op);

//       //-> prepare for coro bases deformation gradient operator
//       //-> only need to save Dm and dNdE
//       for (size_t k = 0; k < adjc_vert_num_; ++k) {
//         dNude_(2*k+0, colon(2*i, 2*i+1)) = Hx(k, colon());
//         dNude_(2*k+1, colon(2*i, 2*i+1)) = Hy(k, colon());
//       }
//     }
//   }
  
//   return 0;
// }

// static inline double two_point_lag_interp(const double x,
//                                           const double x1, const double y1,
//                                           const double x2, const double y2) {
//   return ((x-x1)*y2-(x-x2)*y1)/(x2-x1);
// }
// int quad_stencil::get_coro_bases_df_jac(const matd_t &R, const size_t qrID, matd_t &G) {
//   // G is [4 x 2#adjc_verts] matrix
//   ASSERT(qrID >= 0 && qrID < qr_num_);

//   G = zeros<double>(4, 2*adjc_vert_num_);
//   matd_t temp = zeros<double>(4, 2);

//   if ( !this->isotropic_ ) {
//     matd_t dNude = zeros<double>(2, 2);
//     matd_t dNxde = zeros<double>(1, 2);
//     for (size_t i = 0; i < adjc_vert_num_; ++i) {
//       dNude = dNude_(colon(2*i+0, 2*i+1), colon(2*qrID, 2*qrID+1));
//       dNxde = dNxde_(i, colon(2*qrID, 2*qrID+1));
//       coro_ani_basis_df_2d_jac_(&temp[0], nullptr, nullptr, &R[0],
//                                 nullptr, nullptr, &dNude[0], &dNxde[0], &DmHs_(0, 2*qrID));
//       G(colon(), colon(2*i+0, 2*i+1)) = temp;
//     }
//   } else {
//     matd_t dNude = zeros<double>(1, 2);
//     matd_t dNxde = zeros<double>(1, 2);
//     for (size_t i = 0; i < adjc_vert_num_; ++i) {
//       dNude = dNude_(i, colon(2*qrID, 2*qrID+1));
//       dNxde = dNxde_(i, colon(2*qrID, 2*qrID+1));
//       coro_iso_basis_df_2d_jac_(&temp[0], nullptr, &adjc_rest_(0, i), &R[0],
//                                 nullptr, nullptr, &dNude[0], &dNxde[0], &DmHs_(0, 2*qrID));
//       G(colon(), colon(2*i+0, 2*i+1)) = temp;
//     }
//   }

//   return 0;
// }
