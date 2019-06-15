#include "energy.h"

#include <iostream>
#include <hjlib/math/blas_lapack.h>
#include <zjucad/matrix/lapack.h>
#include <zjucad/matrix/itr_matrix.h>
#include <zjucad/matrix/io.h>
#include <unsupported/Eigen/KroneckerProduct>
#include <unsupported/Eigen/MatrixFunctions>
#include <jtflib/mesh/util.h>
#include <hjlib/math/polar.h>
#include <zjucad/matrix/io.h>
#include <chrono>

#include "mass_matrix.h"
#include "config.h"
#include "geom_util.h"
#include "util.h"
#include "mesh_partition.h"

using namespace std;
using namespace zjucad::matrix;
using namespace Eigen;

namespace bigbang {

extern "C" {

void tet_linear_(double *val, const double *x, const double *Dm, const double *vol, const double *lam, const double *miu);
void tet_linear_jac_(double *jac, const double *x, const double *Dm, const double *vol, const double *lam, const double *miu);
void tet_linear_hes_(double *hes, const double *x, const double *Dm, const double *vol, const double *lam, const double *miu);

void tet_coro_(const double *val, const double *x, const double *Dm, const double *R, const double *vol, const double *lam, const double *mu);
void tet_coro_jac_(const double *jac, const double *x, const double *Dm, const double *R, const double *vol, const double *lam, const double *mu);
void tet_coro_hes_(const double *hes, const double *x, const double *Dm, const double *R, const double *vol, const double *lam, const double *mu);
  
void tet_stvk_(double *val, const double *x, const double *Dm, const double *vol, const double *lam, const double *miu);
void tet_stvk_jac_(double *jac, const double *x, const double *Dm, const double *vol, const double *lam, const double *miu);
void tet_stvk_hes_(double *hes, const double *x, const double *Dm, const double *vol, const double *lam, const double *miu);

void tet_neohookean_(double *val, const double *x, const double *Dm, const double *vol, const double *lam, const double *miu);
void tet_neohookean_jac_(double *jac, const double *x, const double *Dm, const double *vol, const double *lam, const double *miu);
void tet_neohookean_hes_(double *hes, const double *x, const double *Dm, const double *vol, const double *lam, const double *miu);

void tet_arap_(double *val, const double *x, const double *Dm, const double *R, const double *vol);
void tet_arap_jac_(double *jac, const double *x, const double *Dm, const double *R, const double *vol);
void tet_arap_hes_(double *hes, const double *x, const double *Dm, const double *R, const double *vol);

void hex_linear_(double *val, const double *x, const double *h, const double *lam, const double *miu);
void hex_linear_jac_(double *jac, const double *x, const double *h, const double *lam, const double *miu);
void hex_linear_hes_(double *hes, const double *x, const double *h, const double *lam, const double *miu);

void hex_stvk_(double *val, const double *x, const double *h, const double *lam, const double *miu);
void hex_stvk_jac_(double *jac, const double *x, const double *h, const double *lam, const double *miu);
void hex_stvk_hes_(double *hes, const double *x, const double *h, const double *lam, const double *miu);

void calc_edge_length_(double *val, const double *x);
void calc_edge_length_jac_(double *jac, const double *x);
void calc_edge_length_hes_(double *hes, const double *x);

void poly_spring_(double *val, const double *x, const double *d);
void poly_spring_jac_(double *jac, const double *x, const double *d);
void poly_spring_hes_(double *hes, const double *x, const double *d);
  
void mass_spring_(double *val, const double *x, const double *d);
void mass_spring_jac_(double *jac, const double *x, const double *d);
void mass_spring_hes_(double *hes, const double *x, const double *d);

void calc_dih_angle_(double *val, const double *x);
void calc_dih_angle_jac_(double *jac, const double *x);
void calc_dih_angle_hes_(double *hes, const double *x);

void surf_bending_(double *val, const double *x, const double *d, const double *l, const double *area);
void surf_bending_jac_(double *jac, const double *x, const double *d, const double *l, const double *area);
void surf_bending_hes_(double *hes, const double *x, const double *d, const double *l, const double *area);

void bw98_stretch_(double *val, const double *x, const double *invUV, const double *area);
void bw98_stretch_jac_(double *jac, const double *x, const double *invUV, const double *area);
void bw98_stretch_hes_(double *hes, const double *x, const double *invUV, const double *area);

void bw98_shear_(double *val, const double *x, const double *invUV, const double *area);
void bw98_shear_jac_(double *jac, const double *x, const double *invUV, const double *area);
void bw98_shear_hes_(double *hes, const double *x, const double *invUV, const double *area);

void fem_stretch_(double *val, const double *x, const double *Dm, const double *area, const double *k);
void fem_stretch_jac_(double *jac, const double *x, const double *Dm, const double *area, const double *k);
void fem_stretch_hes_(double *hes, const double *x, const double *Dm, const double *area, const double *k);

void line_bending_(double *val, const double *x, const double *d1, const double *d2);
void line_bending_jac_(double *jac, const double *x, const double *d1, const double *d2);
void line_bending_hes_(double *hes, const double *x, const double *d1, const double *d2);

void rod_stretch_(double *val, const double *x, const double *d, const double *Es, const double *r);
void rod_stretch_jac_(double *jac, const double *x, const double *d, const double *Es, const double *r);
void rod_stretch_hes_(double *hes, const double *x, const double *d, const double *Es, const double *r);

void rod_bend_(double *val, const double *q, const double *u, const double *d, const double *E, const double *G, const double *r);
void rod_bend_jac_(double *jac, const double *q, const double *u, const double *d, const double *E, const double *G, const double *r);
void rod_bend_hes_(double *hes, const double *q, const double *u, const double *d, const double *E, const double *G, const double *r);

void rod_couple_(double *val, const double *xq, const double *d, const double *kappa);
void rod_couple_jac_(double *jac, const double *xq, const double *d, const double *kappa);
void rod_couple_hes_(double *hes, const double *xq, const double *d, const double *kappa);

void tri_linear_elas_(double *val, const double *x, const double *D, const double *area, const double *lam, const double *miu);
void tri_linear_elas_jac_(double *jac, const double *x, const double *D, const double *area, const double *lam, const double *miu);
void tri_linear_elas_hes_(double *hes, const double *x, const double *D, const double *area, const double *lam, const double *miu);

void tri_stvk_elas_(double *val, const double *x, const double *D, const double *area, const double *lam, const double *miu);
void tri_stvk_elas_jac_(double *jac, const double *x, const double *D, const double *area, const double *lam, const double *miu);
void tri_stvk_elas_hes_(double *hes, const double *x, const double *D, const double *area, const double *lam, const double *miu);

void polynomial_elas_(double *val, const double *x, const double *Dm, const double *coef);
void polynomial_elas_jac_(double *jac, const double *x, const double *Dm, const double *coef);
void polynomial_elas_hes_(double *hes, const double *x, const double *Dm, const double *coef);

void polynomial_elas_invar_(double *val, const double *x, const double *Dm, const double *coef);
void polynomial_elas_invar_jac_(double *jac, const double *x, const double *Dm, const double *coef);
void polynomial_elas_invar_hes_(double *hes, const double *x, const double *Dm, const double *coef);

void polynomial_elas_invar1_(double *val, const double *x, const double *Dm, const double *coef, const double *area);
void polynomial_elas_invar1_jac_(double *val, const double *x, const double *Dm, const double *coef, const double *area);
void polynomial_elas_invar1_hes_(double *val, const double *x, const double *Dm, const double *coef, const double *area);

  void nnet_tri_elas_(double *val, const double *x, const double *Dm, const double *netW);
  void nnet_tri_elas_jac_(double *jac, const double *x, const double *Dm, const double *netW);
  void nnet_tri_elas_hes_(double *hes, const double *x, const double *Dm, const double *netW);
}

typedef double scalarD;
void const_len_spring(double *out, const double *x, const double *l0) {
  //input
  scalarD a1 = x[0];
  scalarD a2 = x[1];
  scalarD a3 = x[2];
  scalarD b1 = x[3];
  scalarD b2 = x[4];
  scalarD b3 = x[5];
  scalarD r = *l0;

  //temp
  scalarD tt1;
  scalarD tt2;
  scalarD tt3;
  scalarD tt4;

  tt1=a1-b1;
  tt2=a2-b2;
  tt3=a3-b3;
  tt4=1/sqrt(pow(tt3,2)+pow(tt2,2)+pow(tt1,2));
  out[0]=tt1*tt4*r;
  out[1]=tt2*tt4*r;
  out[2]=tt4*tt3*r;
}
void const_len_spring_jac(double *out, const double *x, const double *l0) {
  //input
  scalarD a1 = x[0];
  scalarD a2 = x[1];
  scalarD a3 = x[2];
  scalarD b1 = x[3];
  scalarD b2 = x[4];
  scalarD b3 = x[5];
  scalarD r = *l0;

  //temp
  scalarD tt1;
  scalarD tt2;
  scalarD tt3;
  scalarD tt4;
  scalarD tt5;
  scalarD tt6;
  scalarD tt7;
  scalarD tt8;
  scalarD tt9;
  scalarD tt10;
  scalarD tt11;
  scalarD tt12;
  scalarD tt13;
  scalarD tt14;
  scalarD tt15;
  scalarD tt16;
  scalarD tt17;

  tt1=a1-b1;
  tt2=pow(tt1,2);
  tt3=a2-b2;
  tt4=pow(tt3,2);
  tt5=a3-b3;
  tt6=pow(tt5,2);
  tt7=sqrt(tt6+tt4+tt2);
  tt8=1/pow(tt7,3);
  tt9=1/tt7;
  tt10=tt9*r;
  tt11=-tt1*tt3*tt8*r;
  tt12=-tt1*tt8*tt5*r;
  tt13=-tt3*tt8*tt5*r;
  tt14=-tt9*r;
  tt15=tt1*tt3*tt8*r;
  tt16=tt1*tt8*tt5*r;
  tt17=tt3*tt8*tt5*r;
  out[0]=tt10-tt2*tt8*r;
  out[1]=tt11;
  out[2]=tt12;
  out[3]=tt11;
  out[4]=tt10-tt4*tt8*r;
  out[5]=tt13;
  out[6]=tt12;
  out[7]=tt13;
  out[8]=tt10-tt8*tt6*r;
  out[9]=tt14+tt2*tt8*r;
  out[10]=tt15;
  out[11]=tt16;
  out[12]=tt15;
  out[13]=tt14+tt4*tt8*r;
  out[14]=tt17;
  out[15]=tt16;
  out[16]=tt17;
  out[17]=tt8*tt6*r+tt14;
}
//==============================================================================
implicit_avf::implicit_avf(const mati_t &cell, const matd_t &nods,
                           const double rho, const double h,
                           const std::shared_ptr<Functional<double>> &Ep)
    : rho_(rho), h_(h), h2_(h_*h_), dim_(nods.size()), Ep_(Ep) {
  calc_mass_matrix(cell, nods, rho, nods.size(1), &M_, true);
  invM_.resize(dim_, dim_);
  invM_.setIdentity();
  for (size_t i = 0; i < invM_.rows(); ++i)
    invM_.coeffRef(i, i) = 1.0/M_.coeff(i, i);
  
  xn_ = Eigen::Map<const VectorXd>(&nods[0], dim_);
  vn_.setZero(dim_);
  fn_.setZero(dim_);
  Ep_->Gra(xn_.data(), fn_.data());
}

size_t implicit_avf::Nx() const {
  return dim_;
}

int implicit_avf::Val(const double *x, double *val) const {
  Eigen::Map<const VectorXd> X(x, dim_);

  const VectorXd dv = X-xn_-h_*vn_+h2_/12*invM_*fn_;
  double v1 = dv.dot(M_*dv);

  const VectorXd half_x = (X+xn_)/2;
  double v2 = 0;
  Ep_->Val(half_x.data(), &v2);

  *val += 6/h2_*v1+8*v2;

  return 0;
}

int implicit_avf::Gra(const double *x, double *gra) const {
  Eigen::Map<const VectorXd> X(x, dim_);
  Eigen::Map<VectorXd> G(gra, dim_);

  const VectorXd g1 = 2*M_*(X-xn_-h_*vn_+h2_/12*invM_*fn_);

  const VectorXd half_x = (X+xn_)/2;
  VectorXd g2 = VectorXd::Zero(dim_);
  Ep_->Gra(half_x.data(), g2.data());
  g2 *= 0.5;
    
  G += 6/h2_*g1+8*g2;

  return 0;
}

int implicit_avf::Hes(const double *x, vector<Triplet<double>> *hes) const {
  Eigen::Map<const VectorXd> X(x, dim_);

  const VectorXd half_x = (X+xn_)/2;
  vector<Triplet<double>> trips;
  Ep_->Hes(half_x.data(), &trips);

  hes->reserve(trips.size()+dim_);

  for (auto &it : trips)
    hes->push_back(Triplet<double>(it.row(), it.col(), 4*it.value()));

  for (size_t i = 0; i < M_.rows(); ++i)
    hes->push_back(Triplet<double>(i, i, 12/h2_*M_.coeff(i, i)));
  
  return 0;
}

void implicit_avf::Init(const double *x0, const double *v0) {
  if ( x0 != nullptr ) {
    xn_ = Eigen::Map<const VectorXd>(x0, dim_);
    fn_.setZero(dim_);
    Ep_->Gra(xn_.data(), fn_.data());
  }
  if ( v0 != nullptr )
    vn_ = Eigen::Map<const VectorXd>(v0, dim_);
}

void implicit_avf::Update(const double *x) {
   Eigen::Map<const VectorXd> X(x, dim_);
   vn_ = (2*(X-xn_)/h_-vn_).eval();
   xn_ = X;
   fn_.setZero(dim_);
   Ep_->Gra(xn_.data(), fn_.data());
}
//==============================================================================
momentum_potential_imp_euler::momentum_potential_imp_euler(const mati_t &cell, const matd_t &nods,
                                                           const double rho, const double h, const double w)
  : rho_(rho), h_(h), w_(w), dim_(nods.size()) {
  calc_mass_matrix(cell, nods, rho, nods.size(1), &M_, true);
  xn_ = Eigen::Map<const VectorXd>(&nods[0], dim_);
  vn_.setZero(dim_);
}

momentum_potential_imp_euler::momentum_potential_imp_euler(const matd_t &nods, const SparseMatrix<double> &M,
                                                           const double h, const double w)
    : h_(h), w_(w), dim_(nods.size()), M_(M) {
  xn_ = Eigen::Map<const VectorXd>(&nods[0], dim_);
  vn_.setZero(dim_);
}

size_t momentum_potential_imp_euler::Nx() const {
  return dim_;
}

int momentum_potential_imp_euler::Val(const double *x, double *val) const {
  RETURN_WITH_COND_TRUE(w_ == 0.0);
  Eigen::Map<const VectorXd> X(x, dim_);
  VectorXd dv = (X-xn_)/h_-vn_;
  *val += w_*0.5*dv.dot(M_*dv);
  return 0;
}

int momentum_potential_imp_euler::Gra(const double *x, double *gra) const {
  RETURN_WITH_COND_TRUE(w_ == 0.0);
  Eigen::Map<const VectorXd> X(x, dim_);
  Eigen::Map<VectorXd> g(gra, dim_);
  g += w_*M_*((X-xn_)/h_-vn_)/h_;
  return 0;
}

int momentum_potential_imp_euler::Hes(const double *x, vector<Triplet<double>> *hes) const {
  RETURN_WITH_COND_TRUE(w_ == 0.0);
  const double coeff = w_/(h_*h_);
  for (size_t j = 0; j < M_.outerSize(); ++j) {
    for (SparseMatrix<double>::InnerIterator it(M_, j); it; ++it)
      hes->push_back(Triplet<double>(it.row(), it.col(), coeff*it.value()));
  }
  return 0;
}

void momentum_potential_imp_euler::Init(const double *x0, const double *v0) {
  if ( x0 != nullptr )
    xn_ = Eigen::Map<const VectorXd>(x0, dim_);
  if ( v0 != nullptr )
    vn_ = Eigen::Map<const VectorXd>(v0, dim_);
}

void momentum_potential_imp_euler::Update(const double *x) {
  Eigen::Map<const VectorXd> X(x, dim_);
  vn_ = (X-xn_)/h_;
  xn_ = X;
}

double momentum_potential_imp_euler::QueryKineticEnergy() const {
  return 0.5*vn_.dot(M_*vn_);
}
//==============================================================================
momentum_potential_bdf2::momentum_potential_bdf2(const mati_t &cell, const matd_t &nods,
                                                 const double rho, const double h, const double w)
  : rho_(rho), h_(h), w_(w), dim_(nods.size()) {
  calc_mass_matrix(cell, nods, rho, nods.size(1), &M_, false);
  xn_ = Eigen::Map<const VectorXd>(&nods[0], dim_);
  xnn_ = xn_;
  vn_.setZero(dim_);
  vnn_.setZero(dim_);
}

momentum_potential_bdf2::momentum_potential_bdf2(const matd_t &nods, const SparseMatrix<double> &M,
                                                 const double h, const double w)
    : h_(h), w_(w), dim_(nods.size()), M_(M) {
  xn_ = Eigen::Map<const VectorXd>(&nods[0], dim_);
  xnn_ = xn_;
  vn_.setZero(dim_);
  vnn_.setZero(dim_);  
}

size_t momentum_potential_bdf2::Nx() const {
  return dim_;
}

int momentum_potential_bdf2::Val(const double *x, double *val) const {
  RETURN_WITH_COND_TRUE(w_ == 0.0);
  Eigen::Map<const VectorXd> X(x, dim_);
  VectorXd dv = 1.5/h_*X-2.0/h_*xn_+0.5/h_*xnn_-4.0/3*vn_+1.0/3*vnn_;
  *val += 0.5*w_*dv.dot(M_*dv);
  return 0;
}

int momentum_potential_bdf2::Gra(const double *x, double *gra) const {
  RETURN_WITH_COND_TRUE(w_ == 0.0);
  Eigen::Map<const VectorXd> X(x, dim_);
  Eigen::Map<VectorXd> g(gra, dim_);
  g += w_*1.5/h_*M_*(1.5/h_*X-2.0/h_*xn_+0.5/h_*xnn_-4.0/3*vn_+1.0/3*vnn_);
  return 0;
}

int momentum_potential_bdf2::Hes(const double *x, vector<Triplet<double>> *hes) const {
  RETURN_WITH_COND_TRUE(w_ == 0.0);
  const double coeff = w_*2.25/(h_*h_);
  for (size_t j = 0; j < M_.outerSize(); ++j) {
    for (SparseMatrix<double>::InnerIterator it(M_, j); it; ++it) {
      hes->push_back(Triplet<double>(it.row(), it.col(), coeff*it.value()));
    }
  }
  return 0;
}

void momentum_potential_bdf2::Init(const double *x0, const double *v0) {
  if ( v0 != nullptr ) {
    vn_ = Eigen::Map<const VectorXd>(v0, dim_);
    vnn_ = vn_;
  }
  if ( x0 != nullptr ) {
    xn_ = Eigen::Map<const VectorXd>(x0, dim_);
    xnn_ = xn_-vn_*h_;
  }
}

void momentum_potential_bdf2::Update(const double *x) {
  Eigen::Map<const VectorXd> X(x, dim_);
  vnn_ = vn_;
  vn_ = (3*X-4*xn_+xnn_)/(2*h_);
  xnn_ = xn_;
  xn_ = X;
}

double momentum_potential_bdf2::QueryKineticEnergy() const {
  return 0.5*vn_.dot(M_*vn_);
}
//==============================================================================
gravitational_potential::gravitational_potential(const mati_t &cell, const matd_t &nods,
                                                 const double rho, const double w, const char axis)
    : dim_(nods.size()), w_(w), rd_(nods.size(1)), direction_(axis-'X') {
  ASSERT(direction_ >= 0 && direction_ < rd_);
  calc_mass_matrix(cell, nods, rho, 1, &M_, true);
}

gravitational_potential::gravitational_potential(const matd_t &nods, const SparseMatrix<double> &M,
                                                 const double w, const char axis)
    : dim_(nods.size()), w_(w), rd_(nods.size(1)), direction_(axis-'X') {
  ASSERT(direction_ >= 0 && direction_ < rd_);
  ASSERT(nods.size() == M.cols());
  M_.resize(nods.size(2), nods.size(2));
  M_.setIdentity();
  for (size_t i = 0; i < M_.rows(); ++i)
    M_.coeffRef(i, i) = M.coeff(rd_*i, rd_*i);
}

size_t gravitational_potential::Nx() const {
  return dim_;
}

int gravitational_potential::Val(const double *x, double *val) const {
  RETURN_WITH_COND_TRUE(w_ == 0.0);
  Eigen::Map<const MatrixXd> X(x, rd_, dim_/rd_);

  for (size_t i = 0; i < X.cols(); ++i)
    *val += w_*9.8*M_.coeff(i, i)*X(direction_, i);

  return 0;
}

int gravitational_potential::Gra(const double *x, double *gra) const {
  RETURN_WITH_COND_TRUE(w_ == 0.0);
  Eigen::Map<MatrixXd> g(gra, rd_, dim_/rd_);

  for (size_t i = 0; i < g.cols(); ++i)
    g(direction_, i) += w_*9.8*M_.coeff(i, i);

  return 0;
}

//==============================================================================
torque_3d_potential::torque_3d_potential(const matd_t &nods, const SparseMatrix<double> &M,
                                         const Vector3d &omega, const double w)
    : dim_(nods.size()), w_(w), M_(M), rest_(nods) {
  ASSERT(nods.size(1) == 3);

  omega_ = zeros<double>(3, 1);
  omega_[0] = omega[0];
  omega_[1] = omega[1];
  omega_[2] = omega[2];
  omega_ /= norm(omega_);
  
  //-> calculate mass centers
  bc_ = rest_*ones<double>(rest_.size(2), 1)/rest_.size(2);
}

size_t torque_3d_potential::Nx() const {
  return dim_;
}

int torque_3d_potential::Val(const double *x, double *val) const {
  RETURN_WITH_COND_TRUE(w_ == 0.0);
  itr_matrix<const double *> X(3, Nx()/3, x);

  for (size_t i = 0; i < X.size(2); ++i) {
    const matd_t force = M_.coeff(3*i, 3*i)*cross(omega_, rest_(colon(), i)-bc_);
    *val += -0.5*w_*dot(force, X(colon(), i));
  }

  return 0;
}

int torque_3d_potential::Gra(const double *x, double *gra) const {
  RETURN_WITH_COND_TRUE(w_ == 0.0);
  itr_matrix<double *> G(3, Nx()/3, gra);

  for (size_t i = 0; i < G.size(2); ++i) {
    const matd_t force = M_.coeff(3*i, 3*i)*cross(omega_, rest_(colon(), i)-bc_);
    G(colon(), i) += -0.5*w_*force;
  }  
  return 0;
}

int torque_3d_potential::Hes(const double *x, vector<Triplet<double>> *hes) const {  
  return __LINE__;
}

//==============================================================================
poly_elastic_energy_2d::poly_elastic_energy_2d(const mati_t &tris, const matd_t &nods, const VectorXd &coef,
                                               const MTR mtr, const double mu, const double lam,
                                               const int reg_option,
                                               const double w)
    : coef_(coef), reg_option_(reg_option), fem_2d_energy(tris, nods, mtr, mu, lam, w) {}

size_t poly_elastic_energy_2d::Nx() const {
  return dim_;
}

int poly_elastic_energy_2d::Val(const double *x, double *val) const {
  itr_matrix<const double *> X(2, dim_/2, x);
  double value = 0; matd_t vert = zeros<double>(2, 3);
  for (size_t i = 0; i < tris_.size(2); ++i) {
    vert = X(colon(), tris_(colon(), i));

    switch ( reg_option_ ) {
      case 0: polynomial_elas_invar_(&value, &vert[0], &Dm_(0, i), coef_.data()); break;
      case 1: polynomial_elas_invar1_(&value, &vert[0], &Dm_(0, i), coef_.data(), &area_[i]); break;
      case 2: nnet_tri_elas_(&value, &vert[0], &Dm_(0, i), coef_.data()); break;
    }
    
    value *= w_;
    *val += value;
    
    *(const_cast<double*>(elem_energy_.data()+i)) = value;
  }
  return 0;
}

int poly_elastic_energy_2d::Gra(const double *x, double *gra) const {
  itr_matrix<const double *> X(2, dim_/2, x);
  itr_matrix<double *> G(2, dim_/2, gra);
  matd_t vert = zeros<double>(2, 3), g = zeros<double>(2, 3);
  for (size_t i = 0; i < tris_.size(2); ++i) {
    vert = X(colon(), tris_(colon(), i));

    switch ( reg_option_ ) {
      case 0: polynomial_elas_invar_jac_(&g[0], &vert[0], &Dm_(0, i), coef_.data()); break;
      case 1: polynomial_elas_invar1_jac_(&g[0], &vert[0], &Dm_(0, i), coef_.data(), &area_[i]); break;
      case 2: nnet_tri_elas_jac_(&g[0], &vert[0], &Dm_(0, i), coef_.data()); break;
    }

    g *= w_;
    G(colon(), tris_(colon(), i)) += g;
  }
  return 0;
}

int poly_elastic_energy_2d::Hes(const double *x, vector<Triplet<double>> *hes) const {
  itr_matrix<const double *> X(2, dim_/2, x);
  matd_t H = zeros<double>(6, 6), vert = zeros<double>(2, 3);
  for (size_t i = 0; i < tris_.size(2); ++i) {
    vert = X(colon(), tris_(colon(), i));

    switch ( reg_option_ ) {
      case 0: polynomial_elas_invar_hes_(&H[0], &vert[0], &Dm_(0, i), coef_.data()); break;
      case 1: polynomial_elas_invar1_hes_(&H[0], &vert[0], &Dm_(0, i), coef_.data(), &area_[i]); break;
      case 2: nnet_tri_elas_hes_(&H[0], &vert[0], &Dm_(0, i), coef_.data()); break;
    }

    H *= w_;
    for (size_t p = 0; p < 6; ++p) {
      for (size_t q = 0; q < 6; ++q) {
        const size_t I = 2*tris_(p/2, i)+p%2;
        const size_t J = 2*tris_(q/2, i)+q%2;
        hes->push_back(Triplet<double>(I, J, H(p, q)));
      }
    }
  }
  return 0;
}
//==============================================================================
fem_2d_energy::fem_2d_energy(const mati_t &tris, const matd_t &nods, const MTR mtr, const double mu, const double lam, const double w)
    : tris_(tris), dim_(nods.size()), mtr_(mtr), w_(w) {
  area_.resize(tris.size(2));
  Dm_.resize(4, tris.size(2));
  for (size_t i = 0; i < tris_.size(2); ++i) {
    const matd_t vert = nods(colon(), tris(colon(), i));
    matd_t dm = vert(colon(), colon(1, 2))-vert(colon(), 0)*ones<double>(1, 2);
    matd_t dm_cp = dm;
    inv(dm);
    Dm_(colon(), i) = dm(colon());
    area_[i] = 0.5*std::fabs(det(dm_cp));
  }
  lam_ = lam;
  miu_ = mu;
  elem_energy_ = VectorXd::Zero(tris_.size(2));
  elem_strain_ = VectorXd::Zero(3*tris_.size(2));
  elem_stress_ = VectorXd::Zero(4*tris_.size(2));
}

size_t fem_2d_energy::Nx() const {
  return dim_;
}

static void get_voigt_strain(const double *x, const double *Dm, double *strain) {
  itr_matrix<const double *> X(2, 3, x), DM(2, 2, Dm);
  itr_matrix<double *> E(3, 1, strain);
  matd_t DS = X(colon(), colon(1, 2))-X(colon(), 0)*ones<double>(1, 2);
  matd_t F = DS*DM;
  matd_t s = trans(F)*F; // Cauchy strain //(trans(F)*F-eye<double>(2))/2;
  E[0] = s(0, 0); E[1] = s(1, 1); E[2] = s(0, 1);
}

static void tri_linear_stress(double *P, const double *x, const double *Dm,
                              const double *lam, const double *mu) {
  Eigen::Map<const Matrix<double, 2, 3>> X(x);
  Eigen::Map<const Matrix2d> DM(Dm);
  Eigen::Map<Matrix2d> PS(P);
  Matrix2d Ds;
  Ds.col(0) = X.col(1)-X.col(0);
  Ds.col(1) = X.col(2)-X.col(0);
  Matrix2d F = Ds*DM;
  Matrix2d Id = Matrix2d::Identity();
  PS = (*mu)*(F+F.transpose()-2*Id)+(*lam)*(F-Id).trace()*Id;
}

static void tri_stvk_stress(double *P, const double *x, const double *Dm,
                            const double *lam, const double *mu) {
  Eigen::Map<const Matrix<double, 2, 3>> X(x);
  Eigen::Map<const Matrix2d> DM(Dm);
  Eigen::Map<Matrix2d> PS(P);
  Matrix2d Ds;
  Ds.col(0) = X.col(1)-X.col(0);
  Ds.col(1) = X.col(2)-X.col(0);
  Matrix2d F = Ds*DM;
  Matrix2d Id = Matrix2d::Identity();
  Matrix2d E = (F.transpose()*F-Id)/2;
  PS = F*(2*(*mu)*E  + (*lam)*E.trace()*Id);
}

int fem_2d_energy::Val(const double *x, double *val) const {
  itr_matrix<const double *> X(2, dim_/2, x);
  matd_t vert = zeros<double>(2, 3);  double value = 0;
  for (size_t i = 0; i < tris_.size(2); ++i) {
    vert = X(colon(), tris_(colon(), i));
    switch ( mtr_ ) {
      case LINEAR: tri_linear_elas_(&value, &vert[0], &Dm_(0, i), &area_[i], &lam_, &miu_); break;
      case STVK: tri_stvk_elas_(&value, &vert[0], &Dm_(0, i), &area_[i], &lam_, &miu_); break;
    }
    *(const_cast<double*>(elem_energy_.data()+i)) = w_*value;
    get_voigt_strain(&vert[0], &Dm_(0, i), const_cast<double*>(elem_strain_.data()+3*i));
    switch ( mtr_ ) {
      case LINEAR: tri_linear_stress(const_cast<double*>(elem_stress_.data()+4*i),
                                     &vert[0], &Dm_(0, i), &lam_, &miu_);
      case STVK:   tri_stvk_stress  (const_cast<double*>(elem_stress_.data()+4*i),
                                     &vert[0], &Dm_(0, i), &lam_, &miu_);
    }
    *val += w_*value;
  }
  return 0;
}

int fem_2d_energy::Gra(const double *x, double *gra) const {
  itr_matrix<const double *> X(2, dim_/2, x);
  itr_matrix<double *> G(2, dim_/2, gra);
  matd_t g = zeros<double>(2, 3), vert = zeros<double>(2, 3);
  for (size_t i = 0; i < tris_.size(2); ++i) {
    vert = X(colon(), tris_(colon(), i));
    switch ( mtr_ ) {
      case LINEAR: tri_linear_elas_jac_(&g[0], &vert[0], &Dm_(0, i), &area_[i], &lam_, &miu_); break;
      case STVK: tri_stvk_elas_jac_(&g[0], &vert[0], &Dm_(0, i), &area_[i], &lam_, &miu_); break;
    }
    G(colon(), tris_(colon(), i)) += w_*g;
  }
  return 0;
}

int fem_2d_energy::Hes(const double *x, vector<Triplet<double>> *hes) const {
  itr_matrix<const double *> X(2, dim_/2, x);
  matd_t H = zeros<double>(6, 6), vert = zeros<double>(2, 3);
  for (size_t i = 0; i < tris_.size(2); ++i) {
    vert = X(colon(), tris_(colon(), i));
    switch ( mtr_ ) {
      case LINEAR: tri_linear_elas_hes_(&H[0], nullptr, &Dm_(0, i), &area_[i], &lam_, &miu_); break;
      case STVK: tri_stvk_elas_hes_(&H[0], &vert[0], &Dm_(0, i), &area_[i], &lam_, &miu_); break;
    }
    for (size_t p = 0; p < 6; ++p) {
      for (size_t q = 0; q < 6; ++q) {
        const size_t I = 2*tris_(p/2, i)+p%2;
        const size_t J = 2*tris_(q/2, i)+q%2;
        hes->push_back(Triplet<double>(I, J, w_*H(p, q)));
      }
    }
  }
  return 0;
}
//==============================================================================
elastic_potential::elastic_potential(const mati_t &tets, const matd_t &nods, Material type,
                                     const double Ym, const double Pr, const double w)
  : tets_(tets), type_(type), dim_(nods.size()), w_(w), rest_(nods) {
  vol_.resize(1, tets_.size(2));
  Dm_.resize(9, tets_.size(2));
  #pragma omp parallel for
  for (size_t i = 0; i < tets_.size(2); ++i) {
    matd_t edge = nods(colon(), tets_(colon(1, 3), i))-nods(colon(), tets_(0, i))*ones<double>(1, 3);
    matd_t cp_edge = edge;
    vol_[i] = fabs(det(cp_edge))/6.0;
    if ( inv(edge) )
      cerr << "\tdegenerated tet " << i << endl;
    std::copy(edge.begin(), edge.end(), &Dm_(0, i));
  }
  // calculate \lambda and \miu according
  // to Young's modulus and Poisson ratio
  lam_ = Ym*Pr/((1.0+Pr)*(1.0-2.0*Pr));
  miu_ = Ym/(2.0*(1.0+Pr));
}

size_t elastic_potential::Nx() const {
  return dim_;
}

int elastic_potential::Val(const double *x, double *val) const {
  RETURN_WITH_COND_TRUE(w_ == 0.0);
  itr_matrix<const double*> X(3, dim_/3, x);
  matd_t df(3, 3);
  hj::polar3d decomp;
  for (size_t i = 0; i < tets_.size(2); ++i) {
    matd_t vert = X(colon(), tets_(colon(), i));
    double value = 0;
    switch ( type_ ) {
      case LINEAR:
        tet_linear_(&value, &vert[0], &Dm_(0, i), &vol_[i], &lam_, &miu_);
        break;
      case STVK:
        tet_stvk_(&value, &vert[0], &Dm_(0, i), &vol_[i], &lam_, &miu_);
        break;
      case COROTATIONAL:
        df = (vert(colon(), colon(1, 3))-vert(colon(), 0)*ones<double>(1, 3))*itr_matrix<const double *>(3, 3, &Dm_(0, i));
        decomp(df);
        tet_coro_(&value, &vert[0], &Dm_(0, i), &df[0], &vol_[i], &lam_, &miu_);
        break; 
      case NEOHOOKEAN:
        tet_neohookean_(&value, &vert[0], &Dm_(0, i), &vol_[i], &lam_, &miu_);
        break;
      default:
        break;
    }
    *val += w_*value;
  }
  return 0;
}

int elastic_potential::Gra(const double *x, double *gra) const {
  RETURN_WITH_COND_TRUE(w_ == 0.0);
  itr_matrix<const double *> X(3, dim_/3, x);
  itr_matrix<double *> G(3, dim_/3, gra);
  matd_t vert(3, 4), x0(3, 4), grad(3, 4);
  matd_t df(3, 3);
  hj::polar3d decomp;
  for (size_t i = 0; i < tets_.size(2); ++i) {
    vert = X(colon(), tets_(colon(), i));
    grad = zeros<double>(3, 4);
    switch ( type_ ) {
      case LINEAR:
        tet_linear_jac_(&grad[0], &vert[0], &Dm_(0, i), &vol_[i], &lam_, &miu_);
        break;
      case STVK:
        tet_stvk_jac_(&grad[0], &vert[0], &Dm_(0, i), &vol_[i], &lam_, &miu_);
        break;
      case COROTATIONAL:
        df = (vert(colon(), colon(1, 3))-vert(colon(), 0)*ones<double>(1, 3))*itr_matrix<const double *>(3, 3, &Dm_(0, i));
        decomp(df);
        tet_coro_jac_(&grad[0], &vert[0], &Dm_(0, i), &df[0], &vol_[i], &lam_, &miu_);
        break;
      case NEOHOOKEAN:
        tet_neohookean_jac_(&grad[0], &vert[0], &Dm_(0, i), &vol_[i], &lam_, &miu_);
      default:
        break;
    }
    G(colon(), tets_(colon(), i)) += w_*grad;
  }
  return 0;
}

int elastic_potential::Hes(const double *x, vector<Triplet<double>> *hes) const {
  RETURN_WITH_COND_TRUE(w_ == 0.0);
  itr_matrix<const double *> X(3, dim_/3, x);
  matd_t vert(3, 4), H(12, 12);
  matd_t df(3, 3);
  hj::polar3d decomp;
  for (size_t i = 0; i < tets_.size(2); ++i) {
    vert = X(colon(), tets_(colon(), i));
    H = zeros<double>(12, 12);
    switch ( type_ ) {
      case LINEAR:
        tet_linear_hes_(&H[0], nullptr, &Dm_(0, i), &vol_[i], &lam_, &miu_);
        break;
      case STVK:
        tet_stvk_hes_(&H[0], &vert[0], &Dm_(0, i), &vol_[i], &lam_, &miu_);
        break;
      case COROTATIONAL:
        df = (vert(colon(), colon(1, 3))-vert(colon(), 0)*ones<double>(1, 3))*itr_matrix<const double *>(3, 3, &Dm_(0, i));
        decomp(df);
        tet_coro_hes_(&H[0], &vert[0], &Dm_(0, i), &df[0], &vol_[i], &lam_, &miu_);
        break;
      case NEOHOOKEAN:
        tet_neohookean_hes_(&H[0], &vert[0], &Dm_(0, i), &vol_[i], &lam_, &miu_);
        break;
      default:
        break;
    }
    for (size_t p = 0; p < 12; ++p) {
      for (size_t q = 0; q < 12; ++q) {
        const size_t I = 3*tets_(p/3, i)+p%3;
        const size_t J = 3*tets_(q/3, i)+q%3;
        if ( H(p, q) != 0.0 )
          hes->push_back(Triplet<double>(I, J, w_*H(p, q)));
      }
    }
  }
  return 0;
}
//==============================================================================

void vox_SF(double *val, const double *epsilon);
void vox_SF_jac(double *jac, const double *epsilon);
void quad_SF(double *val, const double *epsilon);
void quad_SF_jac(double *jac, const double *epsilon);
void quad9_SF(double *val, const double *epsilon);
void quad9_SF_jac(double *jac, const double *epsilon);
// void quad9_BL_jac(double *jac, const double *epsilon);
// void quad9_BR_jac(double *jac, const double *epsilon);
// void quad9_TR_jac(double *jac, const double *epsilon);
// void quad9_TL_jac(double *jac, const double *epsilon);

extern "C" {
  void vox_stvk_at_quadr_(double *val, const double *x, const double *H_invDmH, const double *detDmH, const double *gw, const double *mu, const double *lam);
  void vox_stvk_at_quadr_jac_(double *jac, const double *x, const double *H_invDmH, const double *detDmH, const double *gw, const double *mu, const double *lam);
  void vox_stvk_at_quadr_hes_(double *hes, const double *x, const double *H_invDmH, const double *detDmH, const double *gw, const double *mu, const double *lam);

  void vox_neo_at_quadr_(double *val, const double *x, const double *H_invDmH, const double *detDmH, const double *gw, const double *mu, const double *lam);
  void vox_neo_at_quadr_jac_(double *jac, const double *x, const double *H_invDmH, const double *detDmH, const double *gw, const double *mu, const double *lam);
  void vox_neo_at_quadr_hes_(double *hes, const double *x, const double *H_invDmH, const double *detDmH, const double *gw, const double *mu, const double *lam);

  void vox_sta_neo_at_quadr_(double *val, const double *x, const double *H_invDmH, const double *detDmH, const double *gw, const double *mu, const double *lam);
  void vox_sta_neo_at_quadr_jac_(double *jac, const double *x, const double *H_invDmH, const double *detDmH, const double *gw, const double *mu, const double *lam);
  void vox_sta_neo_at_quadr_hes_(double *hes, const double *x, const double *H_invDmH, const double *detDmH, const double *gw, const double *mu, const double *lam);

  void vox_bower_neo_at_quadr_(double *val, const double *x, const double *H_invDmH, const double *detDmH, const double *gw, const double *mu, const double *lam);
  void vox_bower_neo_at_quadr_jac_(double *jac, const double *x, const double *H_invDmH, const double *detDmH, const double *gw, const double *mu, const double *lam);
  void vox_bower_neo_at_quadr_hes_(double *hes, const double *x, const double *H_invDmH, const double *detDmH, const double *gw, const double *mu, const double *lam);

  void vox_linear_at_quadr_(double *val, const double *x, const double *H_invDmH, const double *detDmH, const double *gw, const double *mu, const double *lam);
  void vox_linear_at_quadr_jac_(double *jac, const double *x, const double *H_invDmH, const double *detDmH, const double *gw, const double *mu, const double *lam);
  void vox_linear_at_quadr_hes_(double *hes, const double *x, const double *H_invDmH, const double *detDmH, const double *gw, const double *mu, const double *lam);

  void vox_corotated_at_quadr_(double *val, const double *x, const double *H_invDmH, const double *R, const double *detDmH, const double *gw, const double *mu, const double *lam);
  void vox_corotated_at_quadr_jac_(double *jac, const double *x, const double *H_invDmH, const double *R, const double *detDmH, const double *gw, const double *mu, const double *lam);
  void vox_corotated_at_quadr_hes_(double *hes, const double *x, const double *H_invDmH, const double *R, const double *detDmH, const double *gw, const double *mu, const double *lam);

  void quad_stvk_at_quadr_(double *val, const double *x, const double *H_invDmH, const double *detDmH, const double *gw, const double *mu, const double *lam);
  void quad_stvk_at_quadr_jac_(double *jac, const double *x, const double *H_invDmH, const double *detDmH, const double *gw, const double *mu, const double *lam);
  void quad_stvk_at_quadr_hes_(double *hes, const double *x, const double *H_invDmH, const double *detDmH, const double *gw, const double *mu, const double *lam);

  void quad_neo_at_quadr_(double *val, const double *x, const double *H_invDmH, const double *detDmH, const double *gw, const double *mu, const double *lam);
  void quad_neo_at_quadr_jac_(double *jac, const double *x, const double *H_invDmH, const double *detDmH, const double *gw, const double *mu, const double *lam);
  void quad_neo_at_quadr_hes_(double *hes, const double *x, const double *H_invDmH, const double *detDmH, const double *gw, const double *mu, const double *lam);

  void quad_linear_at_quadr_(double *val, const double *x, const double *H_invDmH, const double *detDmH, const double *gw, const double *mu, const double *lam);
  void quad_linear_at_quadr_jac_(double *jac, const double *x, const double *H_invDmH, const double *detDmH, const double *gw, const double *mu, const double *lam);
  void quad_linear_at_quadr_hes_(double *hes, const double *x, const double *H_invDmH, const double *detDmH, const double *gw, const double *mu, const double *lam);

  void quad_coro_at_quadr_(double *val, const double *x, const double *H_invDmH, const double *R, const double *detDmH, const double *gw, const double *mu, const double *lam);
  void quad_coro_at_quadr_jac_(double *jac, const double *x, const double *H_invDmH, const double *R, const double *detDmH, const double *gw, const double *mu, const double *lam);
  void quad_coro_at_quadr_hes_(double *hes, const double *x, const double *H_invDmH, const double *R, const double *detDmH, const double *gw, const double *mu, const double *lam);

  void quad9_stvk_at_quadr_(double *val, const double *x, const double *H_invDmH, const double *detDmH, const double *gw, const double *mu, const double *lam);
  void quad9_stvk_at_quadr_jac_(double *jac, const double *x, const double *H_invDmH, const double *detDmH, const double *gw, const double *mu, const double *lam);
  void quad9_stvk_at_quadr_hes_(double *hes, const double *x, const double *H_invDmH, const double *detDmH, const double *gw, const double *mu, const double *lam);

  void quad9_neo_at_quadr_(double *val, const double *x, const double *H_invDmH, const double *detDmH, const double *gw, const double *mu, const double *lam);
  void quad9_neo_at_quadr_jac_(double *jac, const double *x, const double *H_invDmH, const double *detDmH, const double *gw, const double *mu, const double *lam);
  void quad9_neo_at_quadr_hes_(double *hes, const double *x, const double *H_invDmH, const double *detDmH, const double *gw, const double *mu, const double *lam);

  void quad9_linear_at_quadr_(double *val, const double *x, const double *H_invDmH, const double *detDmH, const double *gw, const double *mu, const double *lam);
  void quad9_linear_at_quadr_jac_(double *jac, const double *x, const double *H_invDmH, const double *detDmH, const double *gw, const double *mu, const double *lam);
  void quad9_linear_at_quadr_hes_(double *hes, const double *x, const double *H_invDmH, const double *detDmH, const double *gw, const double *mu, const double *lam);

  void quad9_coro_at_quadr_(double *val, const double *x, const double *H_invDmH, const double *R, const double *detDmH, const double *gw, const double *mu, const double *lam);
  void quad9_coro_at_quadr_jac_(double *jac, const double *x, const double *H_invDmH, const double *R, const double *detDmH, const double *gw, const double *mu, const double *lam);
  void quad9_coro_at_quadr_hes_(double *hes, const double *x, const double *H_invDmH, const double *R, const double *detDmH, const double *gw, const double *mu, const double *lam);

}
//===============================================================================
extern "C" {
  void quad_stvk_psi_qr_(double *val, const double *F, const double *mu, const double *lam);
  void quad_stvk_psi_qr_jac_(double *jac, const double *F, const double *mu, const double *lam);
  void quad_stvk_psi_qr_hes_(double *hes, const double *F, const double *mu, const double *lam);

  void quad_neo_psi_qr_(double *val, const double *F, const double *mu, const double *lam);
  void quad_neo_psi_qr_jac_(double *jac, const double *F, const double *mu, const double *lam);
  void quad_neo_psi_qr_hes_(double *hes, const double *F, const double *mu, const double *lam);

  // void quad_lin_psi_qr_(double *val, const double *F, const double *mu, const double *lam);
  // void quad_lin_psi_qr_jac_(double *jac, const double *F, const double *mu, const double *lam);

  // void quad_coro_psi_qr_(double *val, const double *F, const double *mu, const double *lam);
  // void quad_coro_psi_qr_jac_(double *jac, const double *F, const double *mu, const double *lam);

  void quad_neo_q2l_psi_(double *val, const double *F, const double *mu, const double *lam);
  void quad_neo_q2l_psi_jac_(double *jac, const double *F, const double *mu, const double *lam);
  void quad_neo_q2l_psi_hes_(double *hes, const double *F, const double *mu, const double *lam);
  void quad_neo_f_val_at_qr_hes_(double *hes, const double *F, const double *mu, const double *lam);

  void quad_stvk_q2l_psi_(double *val, const double *F, const double *mu, const double *lam);
  void quad_stvk_q2l_psi_jac_(double *jac, const double *F, const double *mu, const double *lam);
  void quad_stvk_q2l_psi_hes_(double *hes, const double *F, const double *mu, const double *lam);
  void quad_stvk_f_val_at_qr_hes_(double *hes, const double *F, const double *mu, const double *lam);

  void quad_lin_q2l_psi_(double *val, const double *F, const double *mu, const double *lam);
  void quad_lin_q2l_psi_jac_(double *jac, const double *F, const double *mu, const double *lam);
  void quad_lin_q2l_psi_hes_(double *hes, const double *F, const double *mu, const double *lam) {
    std::fill(hes, hes+64, 0);
    hes[36] = hes[63] = 2*(*mu)+(*lam);
    hes[39] = hes[60] = (*lam);
    hes[46] = hes[45] = hes[54] = hes[53] = (*mu);
  }
}

//        _______
//       |   |   |
//       |_4_|_3_|
//       |   |   |
//       |_1_|_2_|
//   
quad49_hybrid_energy::quad49_hybrid_energy(const mati_t &quad4, const matd_t &nods, Material type,
                                           const matd_t &lame, const double w)
    : dim_(nods.size()), w_(w), type_(type), lame_(lame), quad4_(quad4) {
  ASSERT(quad4_.size(2) == lame_.size(2));

  //-> extract quad9
  get_quad9_elem(quad4_, quad9_);

  //-> select two quadratures and calc weights
  const double qrs[2] = {-1.0/sqrt(3.0), +1.0/sqrt(3.0)};
  const double qws[2] = {1.0, 1.0};
  QUAD_NUM_ = 2*2;
  quadr_ = zeros<double>(2, QUAD_NUM_);
  qw_ = ones<double>(QUAD_NUM_, 1);
  for (size_t i = 0; i < 2; ++i) {
    for (size_t j = 0; j < 2; ++j) {
      const size_t idx = 2*i+j;
      quadr_(0, idx) = qrs[i];
      quadr_(1, idx) = qrs[j];
      qw_[idx] = qws[i]*qws[j];
    }
  }

  dets_.resize(1, QUAD_NUM_*quad4.size(2));
  H_invDmH4_.resize(2*4, QUAD_NUM_*quad4_.size(2));
  H_invDmH9_.resize(2*9, QUAD_NUM_*quad4_.size(2));

  matd_t dets = zeros<double>(1, QUAD_NUM_*quad4.size(2));
  
  #pragma omp parallel for
  for (size_t i = 0; i < quad4_.size(2); ++i) {

    const matd_t X = nods(colon(), quad4_(colon(), i));
    for (size_t j = 0; j < QUAD_NUM_; ++j) { // for linear basis
      const size_t idx = QUAD_NUM_*i+j;
      matd_t H = zeros<double>(4, 2);
      quad_SF_jac(&H[0], &quadr_(0, j));
      H *= 2;

      matd_t DmH = X*H;
      dets[idx] = fabs(DmH(0,0)*DmH(1,1)-DmH(0,1)*DmH(1,0))/4;

      if ( inv(DmH) ) cerr << "# inv fail" << endl;
      const matd_t HDmH = H*DmH;
      H_invDmH4_(colon(), idx) = HDmH(colon());
    }

    const matd_t Y = nods(colon(), quad9_(colon(), i/4));
    for (size_t j = 0; j < QUAD_NUM_; ++j) { // for quadratic basis
      const size_t idx = QUAD_NUM_*i+j;

      double xi[2] = {0, 0};
      switch ( i%4 ) {
        case 0:
          xi[0] = (quadr_(0, j)-1)/2;
          xi[1] = (quadr_(1, j)-1)/2;
          break;
        case 1:
          xi[0] = (quadr_(0, j)+1)/2;
          xi[1] = (quadr_(1, j)-1)/2;
          break;
        case 2:
          xi[0] = (quadr_(0, j)+1)/2;
          xi[1] = (quadr_(1, j)+1)/2;
          break;
        case 3:
          xi[0] = (quadr_(0, j)-1)/2;
          xi[1] = (quadr_(1, j)+1)/2;
          break;
      }
      matd_t H = zeros<double>(9, 2);
      quad9_SF_jac(&H[0], xi);
      
      matd_t DmH = Y*H;
      dets_[idx] = fabs(DmH(0,0)*DmH(1,1)-DmH(0,1)*DmH(1,0))/4;

      if ( inv(DmH) ) cerr << "# inv fail" << endl;
      const matd_t HDmH = H*DmH;
      H_invDmH9_(colon(), idx) = HDmH(colon());
    }
  }
}

size_t quad49_hybrid_energy::Nx() const {
  return dim_;
}

int quad49_hybrid_energy::Val(const double *x, double *val) const {
  RETURN_WITH_COND_TRUE(w_ == 0.0);
  itr_matrix<const double *> X(2, dim_/2, x);

  const double ONE = 1.0;
  
  #pragma omp parallel for
  for (size_t i = 0; i < quad4_.size(2); ++i) {
    hj::polar2d decomp;
    matd_t DF(2,2), psiF(2,2), phiF(2,2), gradP(2,2), hessP(4,4), deltaF(2,2);
    const matd_t vert4 = X(colon(), quad4_(colon(), i));
    const matd_t vert9 = X(colon(), quad9_(colon(), i/4));

    double value = 0;
    switch ( type_ ) {
      case LINEAR: {
        break;
      }
      case COROTATED: {
        break;
      }
      case STVK: {
        for (size_t j = 0; j < QUAD_NUM_; ++j) {
          const size_t idx = QUAD_NUM_*i+j;

          double vr = 0;
          quad_stvk_at_quadr_(&vr, &vert4[0], &H_invDmH4_(0, idx), &ONE, &ONE, &lame_(0, i), &lame_(1, i));

          //-> calculate scaling
          psiF = vert4*itr_matrix<const double *>(4, 2, &H_invDmH4_(0, idx));
          phiF = vert9*itr_matrix<const double *>(9, 2, &H_invDmH9_(0, idx));
          deltaF = phiF-psiF;

          gradP = zeros<double>(2, 2);
          quad_stvk_psi_qr_jac_(&gradP[0], &psiF[0], &lame_(0, i), &lame_(1, i));

          //-> hessian
          hessP = zeros<double>(4, 4);
          quad_stvk_psi_qr_hes_(&hessP[0], &psiF[0], &lame_(0, i), &lame_(1, i));

          double scale = 0.0;
          if ( fabs(vr) >= 1e-8 )
            scale = 1.0+(dot(gradP(colon()), deltaF(colon()))+0.5*dot(deltaF(colon()), hessP*deltaF(colon())))/vr;
          
          value += vr*scale*dets_[idx]*qw_[j];
        }
        break;
      }
      case NEOHOOKEAN: {
        for (size_t j = 0; j < QUAD_NUM_; ++j) {
          const size_t idx = QUAD_NUM_*i+j;

          double vr = 0;
          quad_neo_at_quadr_(&vr, &vert4[0], &H_invDmH4_(0, idx), &ONE, &ONE, &lame_(0, i), &lame_(1, i));

          //-> calculate scaling
          psiF = vert4*itr_matrix<const double *>(4, 2, &H_invDmH4_(0, idx));
          phiF = vert9*itr_matrix<const double *>(9, 2, &H_invDmH9_(0, idx));
          deltaF = phiF-psiF;

          gradP = zeros<double>(2, 2);
          quad_neo_psi_qr_jac_(&gradP[0], &psiF[0], &lame_(0, i), &lame_(1, i));

          //-> hessian
          hessP = zeros<double>(4, 4);
          quad_neo_psi_qr_hes_(&hessP[0], &psiF[0], &lame_(0, i), &lame_(1, i));

          double scale = 0.0;
          if ( fabs(vr) >= 1e-8 )
            scale = 1.0+(dot(gradP(colon()), deltaF(colon()))+0.5*dot(deltaF(colon()), hessP*deltaF(colon())))/vr;
          
          value += vr*scale*dets_[idx]*qw_[j];
        }
        break;
      }
      default:
        break;
    }
    #pragma omp critical
    {
      *val += w_*value;
    }
  }
  return 0;
}

int quad49_hybrid_energy::Gra(const double *x, double *gra) const {
  RETURN_WITH_COND_TRUE(w_ == 0.0);
  itr_matrix<const double *> X(2, dim_/2, x);
  itr_matrix<double *> G(2, dim_/2, gra);
  
  const double ONE = 1.0;
  
  #pragma omp parallel for
  for (size_t i = 0; i < quad4_.size(2); ++i) {
    hj::polar2d decomp;
    matd_t DF(2,2), psiF(2,2), phiF(2,2), gradP(2,2), hessP(4,4), deltaF(2,2);
    const matd_t vert4 = X(colon(), quad4_(colon(), i));
    const matd_t vert9 = X(colon(), quad9_(colon(), i/4));

    matd_t g = zeros<double>(2, 4);
    switch ( type_ ) {
      case LINEAR: {
        break;
      }
      case COROTATED: {
        break;
      }
      case STVK: {
        for (size_t j = 0; j < QUAD_NUM_; ++j) {
          const size_t idx = QUAD_NUM_*i+j;

          double vr = 0;
          quad_stvk_at_quadr_(&vr, &vert4[0], &H_invDmH4_(0, idx), &ONE, &ONE, &lame_(0, i), &lame_(1, i));
          matd_t gr = zeros<double>(2, 4);
          quad_stvk_at_quadr_jac_(&gr[0], &vert4[0], &H_invDmH4_(0, idx), &ONE, &ONE, &lame_(0, i), &lame_(1, i));

          //-> calculate scaling
          psiF = vert4*itr_matrix<const double *>(4, 2, &H_invDmH4_(0, idx));
          phiF = vert9*itr_matrix<const double *>(9, 2, &H_invDmH9_(0, idx));
          deltaF = phiF-psiF;

          gradP = zeros<double>(2, 2);
          quad_stvk_psi_qr_jac_(&gradP[0], &psiF[0], &lame_(0, i), &lame_(1, i));

          //-> hessian
          hessP = zeros<double>(4, 4);
          quad_stvk_psi_qr_hes_(&hessP[0], &psiF[0], &lame_(0, i), &lame_(1, i));

          double scale = 0.0;
          if ( fabs(vr) >= 1e-8 )
            scale = 1.0+(dot(gradP(colon()), deltaF(colon()))+0.5*dot(deltaF(colon()), hessP*deltaF(colon())))/vr;
          
          g += gr*scale*dets_[idx]*qw_[j];
        }
        break;
      }
      case NEOHOOKEAN: {
        for (size_t j = 0; j < QUAD_NUM_; ++j) {
          const size_t idx = QUAD_NUM_*i+j;

          double vr = 0;
          quad_neo_at_quadr_(&vr, &vert4[0], &H_invDmH4_(0, idx), &ONE, &ONE, &lame_(0, i), &lame_(1, i));
          matd_t gr = zeros<double>(2, 4);
          quad_neo_at_quadr_jac_(&gr[0], &vert4[0], &H_invDmH4_(0, idx), &ONE, &ONE, &lame_(0, i), &lame_(1, i));

          //-> calculate scaling
          psiF = vert4*itr_matrix<const double *>(4, 2, &H_invDmH4_(0, idx));
          phiF = vert9*itr_matrix<const double *>(9, 2, &H_invDmH9_(0, idx));
          deltaF = phiF-psiF;

          gradP = zeros<double>(2, 2);
          quad_neo_psi_qr_jac_(&gradP[0], &psiF[0], &lame_(0, i), &lame_(1, i));

          //-> hessian
          hessP = zeros<double>(4, 4);
          quad_neo_psi_qr_hes_(&hessP[0], &psiF[0], &lame_(0, i), &lame_(1, i));

          double scale = 0.0;
          if ( fabs(vr) >= 1e-8 )
            scale = 1.0+(dot(gradP(colon()), deltaF(colon()))+0.5*dot(deltaF(colon()), hessP*deltaF(colon())))/vr;
          
          g += gr*scale*dets_[idx]*qw_[j];
        }
        break;
      }
      default:
        break;
    }
    #pragma omp critical
    {
      G(colon(), quad4_(colon(), i)) += w_*g;
    }
  }
  return 0;
}

int quad49_hybrid_energy::Hes(const double *x, vector<Triplet<double>> *hes) const {
  RETURN_WITH_COND_TRUE(w_ == 0.0);
  itr_matrix<const double *> X(2, dim_/2, x);

  const double ONE = 1.0;
  
  #pragma omp parallel for
  for (size_t i = 0; i < quad4_.size(2); ++i) {
    hj::polar2d decomp;
    matd_t DF(2,2), psiF(2,2), phiF(2,2), gradP(2,2), hessP(4,4), deltaF(2, 2);
    const matd_t vert4 = X(colon(), quad4_(colon(), i));
    const matd_t vert9 = X(colon(), quad9_(colon(), i/4));

    matd_t H = zeros<double>(8, 8);
    switch ( type_ ) {
      case LINEAR: {
        break;
      }
      case COROTATED: {
        break;
      }
      case STVK: {
        for (size_t j = 0; j < QUAD_NUM_; ++j) {
          const size_t idx = QUAD_NUM_*i+j;

          double vr = 0;
          quad_stvk_at_quadr_(&vr, &vert4[0], &H_invDmH4_(0, idx), &ONE, &ONE, &lame_(0, i), &lame_(1, i));
          matd_t Hr = zeros<double>(8, 8);
          quad_stvk_at_quadr_hes_(&Hr[0], &vert4[0], &H_invDmH4_(0, idx), &ONE, &ONE, &lame_(0, i), &lame_(1, i));
          
          //-> calculate scaling
          psiF = vert4*itr_matrix<const double *>(4, 2, &H_invDmH4_(0, idx));
          phiF = vert9*itr_matrix<const double *>(9, 2, &H_invDmH9_(0, idx));
          deltaF = phiF-psiF;

          gradP = zeros<double>(2, 2);
          quad_stvk_psi_qr_jac_(&gradP[0], &psiF[0], &lame_(0, i), &lame_(1, i));

          //-> hessian
          hessP = zeros<double>(4, 4);
          quad_stvk_psi_qr_hes_(&hessP[0], &psiF[0], &lame_(0, i), &lame_(1, i));

          double scale = 0.0;
          if ( fabs(vr) >= 1e-8 )
            scale = 1.0+(dot(gradP(colon()), deltaF(colon()))+0.5*dot(deltaF(colon()), hessP*deltaF(colon())))/vr;

          H += Hr*scale*dets_[idx]*qw_[j];
        }
        break;
      }
      case NEOHOOKEAN: {
        for (size_t j = 0; j < QUAD_NUM_; ++j) {
          const size_t idx = QUAD_NUM_*i+j;

          double vr = 0;
          quad_neo_at_quadr_(&vr, &vert4[0], &H_invDmH4_(0, idx), &ONE, &ONE, &lame_(0, i), &lame_(1, i));
          matd_t Hr = zeros<double>(8, 8);
          quad_neo_at_quadr_hes_(&Hr[0], &vert4[0], &H_invDmH4_(0, idx), &ONE, &ONE, &lame_(0, i), &lame_(1, i));

          //-> calculate scaling
          psiF = vert4*itr_matrix<const double *>(4, 2, &H_invDmH4_(0, idx));
          phiF = vert9*itr_matrix<const double *>(9, 2, &H_invDmH9_(0, idx));
          deltaF = phiF-psiF;

          gradP = zeros<double>(2, 2);
          quad_neo_psi_qr_jac_(&gradP[0], &psiF[0], &lame_(0, i), &lame_(1, i));

          //-> hessian
          hessP = zeros<double>(4, 4);
          quad_neo_psi_qr_hes_(&hessP[0], &psiF[0], &lame_(0, i), &lame_(1, i));

          double scale = 0.0;
          if ( fabs(vr) >= 1e-8 )
            scale = 1.0+(dot(gradP(colon()), deltaF(colon()))+0.5*dot(deltaF(colon()), hessP*deltaF(colon())))/vr;
          
          H += Hr*scale*dets_[idx]*qw_[j];
        }
        break;
      }
      default:
        break;
    }
    #pragma omp critical
    {
      for (size_t p = 0; p < 8; ++p) {
        for (size_t q = 0; q < 8; ++q) {
          const size_t I = 2*quad4_(p/2, i)+p%2;
          const size_t J = 2*quad4_(q/2, i)+q%2;
          hes->push_back(Triplet<double>(I, J, w_*H(p, q)));
        }
      }
    }
  }
  return 0;
}
//===============================================================================
quad4_taylor_energy::quad4_taylor_energy(const mati_t &quad4, const matd_t &nods, Material type,
                                         const matd_t &lame, const double w, const char hessian)
    : quad49_hybrid_energy(quad4, nods, type, lame, w), hessian_(hessian) {
  opG4_ = zeros<double>( 8, 4*H_invDmH4_.size(2));
  opG9_ = zeros<double>(18, 4*H_invDmH9_.size(2));
  #pragma omp parallel for
  for (size_t j = 0; j < H_invDmH4_.size(2); ++j) {
    matd_t tmp_mat = itr_matrix<const double *>(4, 2, &H_invDmH4_(0, j));
    opG4_(colon(), colon(4*j, 4*j+3)) = kroneckerId<double, 2>(tmp_mat);
    tmp_mat = itr_matrix<const double *>(9, 2, &H_invDmH9_(0, j));
    opG9_(colon(), colon(4*j, 4*j+3)) = kroneckerId<double, 2>(tmp_mat);
  }

  //-> determine whether to eig locally
  debug_flag_ = false;

  //-> record the eigenvalues of each element's hessian matrix
  eigval_ = zeros<double>(8, quad4_.size(2));

  //-> determine if the quad9 vertex is the one of the current quad4 vertex
  in_quad_ = zeros<int>(9, quad4_.size(2));
  in_quad9_ = zeros<int>(4, quad4_.size(2));
  #pragma omp parallel for
  for (size_t i = 0; i < quad4_.size(2); ++i) {
    const mati_t q4 = quad4_(colon(), i);
    const mati_t q9 = quad9_(colon(), i/4);
    for (size_t j = 0; j < q9.size(); ++j) {
      const auto it = std::find(q4.begin(), q4.end(), q9[j]);
      if ( it != q4.end() ) {
        in_quad_(j, i) = it-q4.begin()+1;
        in_quad9_(it-q4.begin(), i) = j;
      }
    }
    ASSERT(sum(in_quad_(colon(), i)) == 10); // 1+2+3+4=10
  }

  vector<Triplet<double>> trips;
  this->HesQuasiNewton(&nods[0], &trips);
  for (auto it : trips) {
    int p = it.row(), q = it.col();
    if ( p > q ) std::swap(p, q);
    sp_pattern_.insert(std::make_pair(p, q));
  }
}

size_t quad4_taylor_energy::Nx() const {
  return dim_;
}

int quad4_taylor_energy::Val(const double *x, double *val) const {
  RETURN_WITH_COND_TRUE(w_ == 0.0);
  itr_matrix<const double *> X(2, dim_/2, x);

  #pragma omp parallel for
  for (size_t i = 0; i < quad4_.size(2); ++i) {
    matd_t dF = zeros<double>(2, 4);
    const matd_t vert4 = X(colon(), quad4_(colon(), i));
    const matd_t vert9 = X(colon(), quad9_(colon(), i/4));

    double value = 0;
    
    for (size_t j = 0; j < QUAD_NUM_; ++j) { //-> for each quadrature
      const size_t idx = QUAD_NUM_*i+j;

      dF = zeros<double>(2, 4);
      dF(colon(), colon(0, 1)) = vert4*itr_matrix<const double *>(4, 2, &H_invDmH4_(0, idx));
      dF(colon(), colon(2, 3)) = vert9*itr_matrix<const double *>(9, 2, &H_invDmH9_(0, idx));

      double vr = 0;
      switch ( type_ ) {
        case LINEAR:
          quad_lin_q2l_psi_(&vr, &dF[0], &lame_(0, i), &lame_(1, i));
          break;
        case COROTATED:
          break;
        case STVK:
          quad_stvk_q2l_psi_(&vr, &dF[0], &lame_(0, i), &lame_(1, i));
          break;
        case NEOHOOKEAN:
          quad_neo_q2l_psi_(&vr, &dF[0], &lame_(0, i), &lame_(1, i));
          break;
        default:
          break;
      }
      
      value += vr*dets_[idx]*qw_[j];
    }
    
    #pragma omp critical
    {
      *val += w_*value;
    }
  }
  return 0;
}

int quad4_taylor_energy::Gra(const double *x, double *gra) const {
  RETURN_WITH_COND_TRUE(w_ == 0.0);
  itr_matrix<const double *> X(2, dim_/2, x);
  itr_matrix<double *> G(2, dim_/2, gra);
  
  #pragma omp parallel for
  for (size_t i = 0; i < quad4_.size(2); ++i) {
    matd_t dF = zeros<double>(2, 4);
    const matd_t vert4 = X(colon(), quad4_(colon(), i));
    const matd_t vert9 = X(colon(), quad9_(colon(), i/4));

    matd_t g4 = zeros<double>(8, 1), g9 = zeros<double>(18, 1);
    matd_t gx4 = g4, gx9 = g9;
    matd_t gf = zeros<double>(8, 1);
    
    for (size_t j = 0; j <QUAD_NUM_; ++j) {
      const size_t idx = QUAD_NUM_*i+j;

      dF = zeros<double>(2, 4);
      dF(colon(), colon(0, 1)) = vert4*itr_matrix<const double *>(4, 2, &H_invDmH4_(0, idx));
      dF(colon(), colon(2, 3)) = vert9*itr_matrix<const double *>(9, 2, &H_invDmH9_(0, idx));

      gf = zeros<double>(8, 1);
      switch ( type_ ) {
        case LINEAR:
          quad_lin_q2l_psi_jac_(&gf[0], &dF[0], &lame_(0, i), &lame_(1, i));
          break;
        case COROTATED:
          break;
        case STVK:
          quad_stvk_q2l_psi_jac_(&gf[0], &dF[0], &lame_(0, i), &lame_(1, i));
          break;
        case NEOHOOKEAN:
          quad_neo_q2l_psi_jac_(&gf[0], &dF[0], &lame_(0, i), &lame_(1, i));
          break;
        default:
          break;
      }

      gx4 = opG4_(colon(), colon(4*idx, 4*idx+3))*gf(colon(0, 3));
      gx9 = opG9_(colon(), colon(4*idx, 4*idx+3))*gf(colon(4, 7));
      g4 += gx4*dets_[idx]*qw_[j];
      g9 += gx9*dets_[idx]*qw_[j];
    }

    #pragma omp critical
    {
      G(colon(), quad4_(colon(), i))   += w_*itr_matrix<const double *>(2, 4, &g4[0]);
      G(colon(), quad9_(colon(), i/4)) += w_*itr_matrix<const double *>(2, 9, &g9[0]);
    }
  }
  return 0;
}

int quad4_taylor_energy::Hes(const double *x, vector<Triplet<double>> *hes) const {
  RETURN_WITH_COND_TRUE(w_ == 0.0);

  if ( hessian_ == 'N' )
    return HesExactNewton(x, hes);

  if ( hessian_ == 'Q' )
    return HesQuasiNewton(x, hes);

  if ( hessian_ == 'T' )
    return HesTruncNewton(x, hes);

  if ( hessian_ == 'S' )
    return HesStiffNewton(x, hes);

  if ( hessian_ == 'M' )
    return HesMixedNewton(x, hes);

  return __LINE__;
}

int quad4_taylor_energy::HesExactNewton(const double *x, vector<Triplet<double>> *hes) const {
  itr_matrix<const double *> X(2, dim_/2, x);

  #pragma omp parallel for
  for (size_t i = 0; i < quad4_.size(2); ++i) {
    matd_t dF = zeros<double>(2, 4);
    const matd_t vert4 = X(colon(), quad4_(colon(), i));
    const matd_t vert9 = X(colon(), quad9_(colon(), i/4));

    matd_t H44 = zeros<double>(8, 8), H99 = zeros<double>(18, 18);
    matd_t H94 = zeros<double>(18, 8), H49 = zeros<double>(8, 18);
    matd_t Hx44 = H44, Hx99 = H99, Hx94 = H94, Hx49 = H49;
    matd_t Hf = zeros<double>(8, 8);
    
    for (size_t j = 0; j < QUAD_NUM_; ++j) {
      const size_t idx = QUAD_NUM_*i+j;

      dF = zeros<double>(2, 4);
      dF(colon(), colon(0, 1)) = vert4*itr_matrix<const double *>(4, 2, &H_invDmH4_(0, idx));
      dF(colon(), colon(2, 3)) = vert9*itr_matrix<const double *>(9, 2, &H_invDmH9_(0, idx));

      Hf = zeros<double>(8, 8);
      switch ( type_ ) {
        case LINEAR:
          quad_lin_q2l_psi_hes_(&Hf[0], &dF[0], &lame_(0, i), &lame_(1, i));
          break;
        case COROTATED:
          break;
        case STVK:
          quad_stvk_q2l_psi_hes_(&Hf[0], &dF[0], &lame_(0, i), &lame_(1, i));
          break;
        case NEOHOOKEAN:
          quad_neo_q2l_psi_hes_(&Hf[0], &dF[0], &lame_(0, i), &lame_(1, i));
          break;
        default:
          break;
      }
      
      Hx44 = opG4_(colon(), colon(4*idx, 4*idx+3))*Hf(colon(0, 3), colon(0, 3))*trans(opG4_(colon(), colon(4*idx, 4*idx+3)));
      Hx99 = opG9_(colon(), colon(4*idx, 4*idx+3))*Hf(colon(4, 7), colon(4, 7))*trans(opG9_(colon(), colon(4*idx, 4*idx+3)));
      Hx94 = opG9_(colon(), colon(4*idx, 4*idx+3))*Hf(colon(0, 3), colon(4, 7))*trans(opG4_(colon(), colon(4*idx, 4*idx+3)));
      Hx49 = opG4_(colon(), colon(4*idx, 4*idx+3))*Hf(colon(4, 7), colon(0, 3))*trans(opG9_(colon(), colon(4*idx, 4*idx+3)));
      H44 += Hx44*dets_[idx]*qw_[j];
      H99 += Hx99*dets_[idx]*qw_[j];
      H94 += Hx94*dets_[idx]*qw_[j];
      H49 += Hx49*dets_[idx]*qw_[j];
    }
    #pragma omp critical
    {
      for (size_t p = 0; p < 8; ++p) {
        for (size_t q = 0; q < 8; ++q) {
          const size_t I = 2*quad4_(p/2, i)+p%2;
          const size_t J = 2*quad4_(q/2, i)+q%2;
          if ( H44(p, q) == 0.0 )
            continue;
          hes->push_back(Triplet<double>(I, J, w_*H44(p, q)));
        }
      }
      for (size_t p = 0; p < 18; ++p) {
        for (size_t q = 0; q < 18; ++q) {
          const size_t I = 2*quad9_(p/2, i/4)+p%2;
          const size_t J = 2*quad9_(q/2, i/4)+q%2;
          if ( H99(p, q) == 0.0 )
            continue;
          hes->push_back(Triplet<double>(I, J, w_*H99(p, q)));
        }
      }
      for (size_t p = 0; p < 18; ++p) {
        for (size_t q = 0; q < 8; ++q) {
          const size_t I = 2*quad9_(p/2, i/4)+p%2;
          const size_t J = 2*quad4_(q/2, i)+q%2;
          if ( H94(p, q) == 0.0 )
            continue;
          hes->push_back(Triplet<double>(I, J, w_*H94(p, q)));
        }
      }
      for (size_t p = 0; p < 8; ++p) {
        for (size_t q = 0; q < 18; ++q) {
          const size_t I = 2*quad4_(p/2, i)+p%2;
          const size_t J = 2*quad9_(q/2, i/4)+q%2;
          if ( H49(p, q) == 0.0 )
            continue;
          hes->push_back(Triplet<double>(I, J, w_*H49(p, q)));
        }
      }
    }
  }
  return 0;
}

int quad4_taylor_energy::HesQuasiNewton(const double *x, vector<Triplet<double>> *hes) const {
  itr_matrix<const double *> X(2, dim_/2, x);

  const double ONE = 1.0;
  
  #pragma omp parallel for
  for (size_t i = 0; i < quad4_.size(2); ++i) {
    const matd_t vert4 = X(colon(), quad4_(colon(), i));

    matd_t H  = zeros<double>(8, 8), Hr = zeros<double>(8, 8);
    for (size_t j = 0; j < QUAD_NUM_; ++j) {
      const size_t idx = QUAD_NUM_*i+j;
      Hr = zeros<double>(8, 8);
      switch ( type_ ) {
        case LINEAR:
          quad_linear_at_quadr_hes_(&Hr[0], &vert4[0], &H_invDmH4_(0, idx), &ONE, &ONE, &lame_(0, i), &lame_(1, i));
          break;
        case COROTATED:
          break;
        case STVK:
          quad_stvk_at_quadr_hes_(&Hr[0], &vert4[0], &H_invDmH4_(0, idx), &ONE, &ONE, &lame_(0, i), &lame_(1, i));
          break;
        case NEOHOOKEAN:
          quad_neo_at_quadr_hes_(&Hr[0], &vert4[0], &H_invDmH4_(0, idx), &ONE, &ONE, &lame_(0, i), &lame_(1, i));
          break;
        default:
          break;
      }      
      H += Hr*dets_[idx]*qw_[j];
#if 0
      const matd_t dF = vert4*itr_matrix<const double *>(4, 2, &H_invDmH4_(0, idx));
      matd_t Hf = zeros<double>(4, 4);
      quad_neo_f_val_at_qr_hes_(&Hf[0], &dF[0], &lame_(0, i), &lame_(1, i));
      const matd_t Hx = opG4_(colon(), colon(4*idx, 4*idx+3))*Hf*trans(opG4_(colon(), colon(4*idx, 4*idx+3)));
      cout << norm(Hx-Hr) << endl;
      getchar();
#endif
    }

    if ( debug_flag_ ) {
      matd_t ev = zeros<double>(H.size(2), 1);
      matd_t copyH = H;
      eig(copyH, ev);
      eigval_(colon(), i) = ev;
    }

    #pragma omp critical
    {
      for (size_t p = 0; p < 8; ++p) {
        for (size_t q = 0; q < 8; ++q) {
          const size_t I = 2*quad4_(p/2, i)+p%2;
          const size_t J = 2*quad4_(q/2, i)+q%2;
          hes->push_back(Triplet<double>(I, J, w_*H(p, q)));
        }
      }
    }
  }
  return 0;
}

int quad4_taylor_energy::HesTruncNewton(const double *x, vector<Triplet<double>> *hes) const {
  itr_matrix<const double *> X(2, dim_/2, x);
  
  #pragma omp parallel for
  for (size_t i = 0; i < quad4_.size(2); ++i) {
    matd_t dF = zeros<double>(2, 4);
    const matd_t vert4 = X(colon(), quad4_(colon(), i));
    const matd_t vert9 = X(colon(), quad9_(colon(), i/4));

    matd_t H44 = zeros<double>(8, 8), H99 = zeros<double>(18, 18);
    matd_t H94 = zeros<double>(18, 8), H49 = zeros<double>(8, 18);
    matd_t Hx44 = H44, Hx99 = H99, Hx94 = H94, Hx49 = H49;
    matd_t Hf = zeros<double>(8, 8);
    
    for (size_t j = 0; j < QUAD_NUM_; ++j) {
      const size_t idx = QUAD_NUM_*i+j;

      dF = zeros<double>(2, 4);
      dF(colon(), colon(0, 1)) = vert4*itr_matrix<const double *>(4, 2, &H_invDmH4_(0, idx));
      dF(colon(), colon(2, 3)) = vert9*itr_matrix<const double *>(9, 2, &H_invDmH9_(0, idx));

      Hf = zeros<double>(8, 8);
      switch ( type_ ) {
        case LINEAR:
          quad_lin_q2l_psi_hes_(&Hf[0], &dF[0], &lame_(0, i), &lame_(1, i));
          break;
        case COROTATED:
          break;
        case STVK:
          quad_stvk_q2l_psi_hes_(&Hf[0], &dF[0], &lame_(0, i), &lame_(1, i));
          break;
        case NEOHOOKEAN:
          quad_neo_q2l_psi_hes_(&Hf[0], &dF[0], &lame_(0, i), &lame_(1, i));
          break;
        default:
          break;
      }
      
      Hx44 = opG4_(colon(), colon(4*idx, 4*idx+3))*Hf(colon(0, 3), colon(0, 3))*trans(opG4_(colon(), colon(4*idx, 4*idx+3)));
      Hx99 = opG9_(colon(), colon(4*idx, 4*idx+3))*Hf(colon(4, 7), colon(4, 7))*trans(opG9_(colon(), colon(4*idx, 4*idx+3)));
      Hx94 = opG9_(colon(), colon(4*idx, 4*idx+3))*Hf(colon(0, 3), colon(4, 7))*trans(opG4_(colon(), colon(4*idx, 4*idx+3)));
      Hx49 = opG4_(colon(), colon(4*idx, 4*idx+3))*Hf(colon(4, 7), colon(0, 3))*trans(opG9_(colon(), colon(4*idx, 4*idx+3)));
      H44 += Hx44*dets_[idx]*qw_[j];
      H99 += Hx99*dets_[idx]*qw_[j];
      H94 += Hx94*dets_[idx]*qw_[j];
      H49 += Hx49*dets_[idx]*qw_[j];
    }

    matd_t H = H44;
    for (size_t p = 0; p < 18; ++p) {
      for (size_t q = 0; q < 18; ++q) {
        const int posI = in_quad_(p/2, i);
        const int posJ = in_quad_(q/2, i);
        if ( posI && posJ ) {
          H(2*(posI-1)+p%2, 2*(posJ-1)+q%2) += H99(p, q); 
        }
      }
    }
    for (size_t p = 0; p < 18; ++p) {
      for (size_t q = 0; q < 8; ++q) {
        const int posI = in_quad_(p/2, i);
        if ( posI )
          H(2*(posI-1)+p%2, q) += H94(p, q);
      }
    }
    for (size_t p = 0; p < 8; ++p) {
      for (size_t q = 0; q < 18; ++q) {
        const int posJ = in_quad_(q/2, i);
        if ( posJ )
          H(p, 2*(posJ-1)+q%2) += H49(p, q);
      }
    }

    if ( debug_flag_ ) {
      matd_t Ht  = zeros<double>(8, 8), Hr = zeros<double>(8, 8);
      const double ONE = 1.0;
      for (size_t j = 0; j < QUAD_NUM_; ++j) {
        const size_t idx = QUAD_NUM_*i+j;
        Hr = zeros<double>(8, 8);
        switch ( type_ ) {
          case LINEAR:
            quad_linear_at_quadr_hes_(&Hr[0], &vert4[0], &H_invDmH4_(0, idx), &ONE, &ONE, &lame_(0, i), &lame_(1, i));
            break;
          case COROTATED:
            break;
          case STVK:
            quad_stvk_at_quadr_hes_(&Hr[0], &vert4[0], &H_invDmH4_(0, idx), &ONE, &ONE, &lame_(0, i), &lame_(1, i));
            break;
          case NEOHOOKEAN:
            quad_neo_at_quadr_hes_(&Hr[0], &vert4[0], &H_invDmH4_(0, idx), &ONE, &ONE, &lame_(0, i), &lame_(1, i));
            break;
          default:
            break;
        }      
        Ht += Hr*dets_[idx]*qw_[j];
      }
      matd_t diffH = H-Ht;
      matd_t ev = zeros<double>(H.size(2), 1);
      eig(diffH, ev);
      eigval_(colon(), i) = ev;
    }
    
    #pragma omp critical
    {
      for (size_t p = 0; p < 8; ++p) {
        for (size_t q = 0; q < 8; ++q) {
          const size_t I = 2*quad4_(p/2, i)+p%2;
          const size_t J = 2*quad4_(q/2, i)+q%2;
          if ( H(p, q) == 0.0 )
            continue;
          hes->push_back(Triplet<double>(I, J, w_*H(p, q)));
        }
      }
    }
  }
  return 0;
}

int quad4_taylor_energy::HesStiffNewton(const double *x, vector<Triplet<double>> *hes) const {
  itr_matrix<const double *> X(2, dim_/2, x);

  trips_bak_.clear();
  
  #pragma omp parallel for
  for (size_t i = 0; i < quad4_.size(2); ++i) {
    matd_t dF = zeros<double>(2, 4);
    const matd_t vert4 = X(colon(), quad4_(colon(), i));
    const matd_t vert9 = X(colon(), quad9_(colon(), i/4));

    matd_t H44 = zeros<double>(8, 8), H99 = zeros<double>(18, 18);
    matd_t H94 = zeros<double>(18, 8), H49 = zeros<double>(8, 18);
    matd_t Hx44 = H44, Hx99 = H99, Hx94 = H94, Hx49 = H49;
    matd_t Hf = zeros<double>(8, 8);
    
    for (size_t j = 0; j < QUAD_NUM_; ++j) {
      const size_t idx = QUAD_NUM_*i+j;

      dF = zeros<double>(2, 4);
      dF(colon(), colon(0, 1)) = vert4*itr_matrix<const double *>(4, 2, &H_invDmH4_(0, idx));
      dF(colon(), colon(2, 3)) = vert9*itr_matrix<const double *>(9, 2, &H_invDmH9_(0, idx));

      Hf = zeros<double>(8, 8);
      switch ( type_ ) {
        case LINEAR:
          quad_lin_q2l_psi_hes_(&Hf[0], &dF[0], &lame_(0, i), &lame_(1, i));
          break;
        case COROTATED:
          break;
        case STVK:
          quad_stvk_q2l_psi_hes_(&Hf[0], &dF[0], &lame_(0, i), &lame_(1, i));
          break;
        case NEOHOOKEAN:
          quad_neo_q2l_psi_hes_(&Hf[0], &dF[0], &lame_(0, i), &lame_(1, i));
          break;
        default:
          break;
      }
      
      Hx44 = opG4_(colon(), colon(4*idx, 4*idx+3))*Hf(colon(0, 3), colon(0, 3))*trans(opG4_(colon(), colon(4*idx, 4*idx+3)));
      Hx99 = opG9_(colon(), colon(4*idx, 4*idx+3))*Hf(colon(4, 7), colon(4, 7))*trans(opG9_(colon(), colon(4*idx, 4*idx+3)));
      Hx94 = opG9_(colon(), colon(4*idx, 4*idx+3))*Hf(colon(0, 3), colon(4, 7))*trans(opG4_(colon(), colon(4*idx, 4*idx+3)));
      Hx49 = opG4_(colon(), colon(4*idx, 4*idx+3))*Hf(colon(4, 7), colon(0, 3))*trans(opG9_(colon(), colon(4*idx, 4*idx+3)));
      H44 += Hx44*dets_[idx]*qw_[j];
      H99 += Hx99*dets_[idx]*qw_[j];
      H94 += Hx94*dets_[idx]*qw_[j];
      H49 += Hx49*dets_[idx]*qw_[j];
    }

    matd_t H = H99;
    for (size_t p = 0; p < 8; ++p) {
      for (size_t q = 0; q < 8; ++q) {
        const int posI = 2*in_quad9_(p/2, i)+p%2;
        const int posJ = 2*in_quad9_(q/2, i)+q%2;
        H(posI, posJ) += H44(p, q);
      }
    }
    for (size_t p = 0; p < 18; ++p) {
      for (size_t q = 0; q < 8; ++q) {
        const int posJ = 2*in_quad9_(q/2, i)+q%2;
        H(p, posJ) += H94(p, q);
      }
    }
    for (size_t p = 0; p < 8; ++p) {
      for (size_t q = 0; q < 18; ++q) {
        const int posI = 2*in_quad9_(p/2, i)+p%2;
        H(posI, q) += H49(p, q);
      }
    }

    //-> remove related entries not in the sparse pattern
    matd_t cacheH = zeros<double>(18, 18);
    for (size_t p = 0; p < 18; ++p) {
      for (size_t q = p+1; q < 18; ++q) {
        int I = 2*quad9_(p/2, i/4)+p%2;
        int J = 2*quad9_(q/2, i/4)+q%2;
        if ( I > J ) std::swap(I, J);
        if ( sp_pattern_.find(make_pair(I, J)) == sp_pattern_.end() ) {
          cacheH(p, q) = cacheH(q, p) = H(p, q);
          H(p, q) = H(q, p) = 0;
        }
      }
    }
    
    #pragma omp critical
    {
      for (size_t p = 0; p < 18; ++p) {
        for (size_t q = 0; q < 18; ++q) {
          const size_t I = 2*quad9_(p/2, i/4)+p%2;
          const size_t J = 2*quad9_(q/2, i/4)+q%2;
          if ( H(p, q) == 0.0 )
            continue;
          hes->push_back(Triplet<double>(I, J, w_*H(p, q)));
        }
      }
      for (size_t p = 0; p < 18; ++p) {
        for (size_t q = 0; q < 18; ++q) {
          const size_t I = 2*quad9_(p/2, i/4)+p%2;
          const size_t J = 2*quad9_(q/2, i/4)+q%2;
          if ( cacheH(p, q) == 0.0 )
            continue;
          trips_bak_.push_back(Triplet<double>(I, J, w_*cacheH(p, q)));
        }
      }
    }
  }
  return 0;
}

int quad4_taylor_energy::HesMixedNewton(const double *x, vector<Triplet<double>> *hes) const {
  itr_matrix<const double *> X(2, dim_/2, x);

  trips_bak_M_.clear();

  const double ONE = 1.0;
  
  #pragma omp parallel for
  for (size_t i = 0; i < quad4_.size(2); ++i) {
    matd_t dF = zeros<double>(2, 4);
    const matd_t vert4 = X(colon(), quad4_(colon(), i));
    const matd_t vert9 = X(colon(), quad9_(colon(), i/4));

    //-> for H2
    matd_t H44 = zeros<double>(8, 8), H99 = zeros<double>(18, 18);
    matd_t H94 = zeros<double>(18, 8), H49 = zeros<double>(8, 18);
    matd_t Hx44 = H44, Hx99 = H99, Hx94 = H94, Hx49 = H49;
    matd_t Hf = zeros<double>(8, 8);

    //-> for H1
    matd_t Hq = zeros<double>(8, 8), Hr = zeros<double>(8, 8);
    
    for (size_t j = 0; j < QUAD_NUM_; ++j) {
      const size_t idx = QUAD_NUM_*i+j;

      dF = zeros<double>(2, 4);
      dF(colon(), colon(0, 1)) = vert4*itr_matrix<const double *>(4, 2, &H_invDmH4_(0, idx));
      dF(colon(), colon(2, 3)) = vert9*itr_matrix<const double *>(9, 2, &H_invDmH9_(0, idx));

      Hf = zeros<double>(8, 8);
      Hr = zeros<double>(8, 8);

      switch ( type_ ) {
        case LINEAR:
          quad_lin_q2l_psi_hes_(&Hf[0], &dF[0], &lame_(0, i), &lame_(1, i));
          quad_linear_at_quadr_hes_(&Hr[0], &vert4[0], &H_invDmH4_(0, idx), &ONE, &ONE, &lame_(0, i), &lame_(1, i));
          break;
        case COROTATED:
          break;
        case STVK:
          quad_stvk_q2l_psi_hes_(&Hf[0], &dF[0], &lame_(0, i), &lame_(1, i));
          quad_stvk_at_quadr_hes_(&Hr[0], &vert4[0], &H_invDmH4_(0, idx), &ONE, &ONE, &lame_(0, i), &lame_(1, i));
          break;
        case NEOHOOKEAN:
          quad_neo_q2l_psi_hes_(&Hf[0], &dF[0], &lame_(0, i), &lame_(1, i));
          quad_neo_at_quadr_hes_(&Hr[0], &vert4[0], &H_invDmH4_(0, idx), &ONE, &ONE, &lame_(0, i), &lame_(1, i));
          break;
        default:
          break;
      }
      
      Hq += Hr*dets_[idx]*qw_[j];
      
      Hx44 = opG4_(colon(), colon(4*idx, 4*idx+3))*Hf(colon(0, 3), colon(0, 3))*trans(opG4_(colon(), colon(4*idx, 4*idx+3)));
      Hx99 = opG9_(colon(), colon(4*idx, 4*idx+3))*Hf(colon(4, 7), colon(4, 7))*trans(opG9_(colon(), colon(4*idx, 4*idx+3)));
      Hx94 = opG9_(colon(), colon(4*idx, 4*idx+3))*Hf(colon(0, 3), colon(4, 7))*trans(opG4_(colon(), colon(4*idx, 4*idx+3)));
      Hx49 = opG4_(colon(), colon(4*idx, 4*idx+3))*Hf(colon(4, 7), colon(0, 3))*trans(opG9_(colon(), colon(4*idx, 4*idx+3)));
      H44 += Hx44*dets_[idx]*qw_[j];
      H99 += Hx99*dets_[idx]*qw_[j];
      H94 += Hx94*dets_[idx]*qw_[j];
      H49 += Hx49*dets_[idx]*qw_[j];
    }

    //-> exact H2
    matd_t H = H99;
    for (size_t p = 0; p < 8; ++p) {
      for (size_t q = 0; q < 8; ++q) {
        const int posI = 2*in_quad9_(p/2, i)+p%2;
        const int posJ = 2*in_quad9_(q/2, i)+q%2;
        H(posI, posJ) += H44(p, q);
      }
    }
    for (size_t p = 0; p < 18; ++p) {
      for (size_t q = 0; q < 8; ++q) {
        const int posJ = 2*in_quad9_(q/2, i)+q%2;
        H(p, posJ) += H94(p, q);
      }
    }
    for (size_t p = 0; p < 8; ++p) {
      for (size_t q = 0; q < 18; ++q) {
        const int posI = 2*in_quad9_(p/2, i)+p%2;
        H(posI, q) += H49(p, q);
      }
    }

    //-> quasi H1
    matd_t HQ = zeros<double>(18, 18);
    for (size_t p = 0; p < 8; ++p) {
      for (size_t q = 0; q < 8; ++q) {
        const int posI = 2*in_quad9_(p/2, i)+p%2;
        const int posJ = 2*in_quad9_(q/2, i)+q%2;
        HQ(posI, posJ) += Hq(p, q);
      }
    }

    //-> get the difference
    const matd_t diffH = H-HQ;
    
    #pragma omp critical
    {
      for (size_t p = 0; p < 8; ++p) {
        for (size_t q = 0; q < 8; ++q) {
          const size_t I = 2*quad4_(p/2, i)+p%2;
          const size_t J = 2*quad4_(q/2, i)+q%2;
          hes->push_back(Triplet<double>(I, J, w_*Hq(p, q)));
        }
      }
      for (size_t p = 0; p < 18; ++p) {
        for (size_t q = 0; q < 18; ++q) {
          const size_t I = 2*quad9_(p/2, i/4)+p%2;
          const size_t J = 2*quad9_(q/2, i/4)+q%2;
          if ( diffH(p, q) == 0.0 )
            continue;
          trips_bak_M_.push_back(Triplet<double>(I, J, w_*diffH(p, q)));
        }
      }
    }
  }
  return 0;
}

void quad4_taylor_energy::set_hess_flag(const char flag) {
  ASSERT(flag == 'N' || flag == 'Q' || flag == 'T' || flag == 'S' || flag == 'M' );
  hessian_ = flag;
}

const matd_t& quad4_taylor_energy::get_elem_eig_vals() const {
  return eigval_;
}
//===============================================================================
quad4_elastic_energy::quad4_elastic_energy(const mati_t &quad, const matd_t &nods, Material type,
                                           const matd_t &lame, const double w)
    : dim_(nods.size()), w_(w), type_(type), lame_(lame), quad_(quad) {
  ASSERT(quad_.size(2) == lame_.size(2));

  //-> select quadratures and calc weights
  const double qrs[2] = {-1.0/sqrt(3.0), +1.0/sqrt(3.0)};
  const double qws[2] = {1.0, 1.0};
  QUAD_NUM_ = 2*2;
  quadr_ = zeros<double>(2, QUAD_NUM_);
  qw_ = ones<double>(QUAD_NUM_, 1);
  for (size_t i = 0; i < 2; ++i) {
    for (size_t j = 0; j < 2; ++j) {
      const size_t idx = 2*i+j;
      quadr_(0, idx) = qrs[i];
      quadr_(1, idx) = qrs[j];
      qw_[idx] = qws[i]*qws[j];
    }
  }

  detDmH_.resize(1, QUAD_NUM_*quad_.size(2));
  H_invDmH_.resize(8, QUAD_NUM_*quad_.size(2));

  #pragma omp parallel for
  for (size_t i = 0; i < quad_.size(2); ++i) {
    const matd_t X = nods(colon(), quad_(colon(), i));

    for (size_t j = 0; j < QUAD_NUM_; ++j) { // for four quadratures
      const size_t idx = QUAD_NUM_*i+j;
      matd_t H = zeros<double>(4, 2);
      quad_SF_jac(&H[0], &quadr_(0, j));

      matd_t DmH = X*H;
      matd_t cpDmH = DmH;

      detDmH_[idx] = fabs(det(cpDmH));

      if ( inv(DmH) )
        cerr << "# inv fail" << endl;
      matd_t HDmH = H*DmH;
      H_invDmH_(colon(), idx) = HDmH(colon());
    }
  }
  cout << "# quad4 elastic energy assembled done" << endl;
}

size_t quad4_elastic_energy::Nx() const {
  return dim_;
}

int quad4_elastic_energy::Val(const double *x, double *val) const {
  RETURN_WITH_COND_TRUE(w_ == 0.0);
  itr_matrix<const double *> X(2, dim_/2, x);

  #pragma omp parallel for
  for (size_t i = 0; i < quad_.size(2); ++i) {
    hj::polar2d decomp;
    matd_t vert(2, 4), DF(2, 2);
    vert = X(colon(), quad_(colon(), i));
    double value = 0;
    switch ( type_ ) {
      case LINEAR: {
        for (size_t j = 0; j < QUAD_NUM_; ++j) {
          const size_t idx = QUAD_NUM_*i+j;
          double vr = 0;
          quad_linear_at_quadr_(&vr, &vert[0], &H_invDmH_(0, idx), &detDmH_[idx], &qw_[j], &lame_(0, i), &lame_(1, i));
          value += vr;
        }
        break;
      }
      case COROTATED: {
        for (size_t j = 0; j < QUAD_NUM_; ++j) {
          const size_t idx = QUAD_NUM_*i+j;
          double vr = 0;
          DF = vert*itr_matrix<const double *>(4, 2, &H_invDmH_(0, idx));
          decomp(DF);
          quad_coro_at_quadr_(&vr, &vert[0], &H_invDmH_(0, idx), &DF[0], &detDmH_[idx], &qw_[j], &lame_(0, i), &lame_(1, i));
          value += vr;
        }
        break;
      }
      case STVK: {
        for (size_t j = 0; j < QUAD_NUM_; ++j) {
          const size_t idx = QUAD_NUM_*i+j;
          double vr = 0;
          quad_stvk_at_quadr_(&vr, &vert[0], &H_invDmH_(0, idx), &detDmH_[idx], &qw_[j], &lame_(0, i), &lame_(1, i));
          value += vr;
        }
        break;
      }
      case NEOHOOKEAN: {
        for (size_t j = 0; j < QUAD_NUM_; ++j) {
          const size_t idx = QUAD_NUM_*i+j;
          double vr = 0;
          quad_neo_at_quadr_(&vr, &vert[0], &H_invDmH_(0, idx), &detDmH_[idx], &qw_[j], &lame_(0, i), &lame_(1, i));
          value += vr;
        }
        break;
      }
      default:
        break;
    }
    #pragma omp critical
    {
      *val += w_*value;
    }
  }
  return 0;
}

int quad4_elastic_energy::Gra(const double *x, double *gra) const {
  RETURN_WITH_COND_TRUE(w_ == 0.0);
  itr_matrix<const double *> X(2, dim_/2, x);
  itr_matrix<double *> G(2, dim_/2, gra);

  #pragma omp parallel for
  for (size_t i = 0; i < quad_.size(2); ++i) {
    hj::polar2d decomp;
    matd_t vert(2, 4), g(2, 4), gr(2, 4), DF(2, 2);
    vert = X(colon(), quad_(colon(), i));
    g = zeros<double>(2, 4);
    switch ( type_ ) {
      case LINEAR: {
        for (size_t j = 0; j < QUAD_NUM_; ++j) {
          const size_t idx = QUAD_NUM_*i+j;
          gr = zeros<double>(2, 4);
          quad_linear_at_quadr_jac_(&gr[0], &vert[0], &H_invDmH_(0, idx), &detDmH_[idx], &qw_[j], &lame_(0, i), &lame_(1, i));
          g += gr;
        }
        break;
      }
      case COROTATED: {
        for (size_t j = 0; j < QUAD_NUM_; ++j) {
          const size_t idx = QUAD_NUM_*i+j;
          gr = zeros<double>(2, 4);
          DF = vert*itr_matrix<const double *>(4, 2, &H_invDmH_(0, idx));
          decomp(DF);
          quad_coro_at_quadr_jac_(&gr[0], &vert[0], &H_invDmH_(0, idx), &DF[0], &detDmH_[idx], &qw_[j], &lame_(0, i), &lame_(1, i));
          g += gr;
        }
        break;
      }
      case STVK: {
        for (size_t j = 0; j < QUAD_NUM_; ++j) {
          const size_t idx = QUAD_NUM_*i+j;
          gr = zeros<double>(2, 4);
          quad_stvk_at_quadr_jac_(&gr[0], &vert[0], &H_invDmH_(0, idx), &detDmH_[idx], &qw_[j], &lame_(0, i), &lame_(1, i));
          g += gr;
        }
        break;
      }
      case NEOHOOKEAN: {
        for (size_t j = 0; j < QUAD_NUM_; ++j) {
          const size_t idx = QUAD_NUM_*i+j;
          gr = zeros<double>(2, 4);
          quad_neo_at_quadr_jac_(&gr[0], &vert[0], &H_invDmH_(0, idx), &detDmH_[idx], &qw_[j], &lame_(0, i), &lame_(1, i));
          g += gr;
        }
        break;
      }
      default:
        break;
    }
    #pragma omp critical
    {
      G(colon(), quad_(colon(), i)) += w_*g;
    }
  }
  return 0;
}

int quad4_elastic_energy::Hes(const double *x, vector<Triplet<double>> *hes) const {
  RETURN_WITH_COND_TRUE(w_ == 0.0);
  itr_matrix<const double *> X(2, dim_/2, x);

  #pragma omp parallel for
  for (size_t i = 0; i < quad_.size(2); ++i) {
    hj::polar2d decomp;
    matd_t vert(2, 4), H(8, 8), Hr(8, 8), DF(2, 2);
    vert = X(colon(), quad_(colon(), i));
    H = zeros<double>(8, 8);
    switch ( type_ ) {
      case LINEAR: {
        for (size_t j = 0; j < QUAD_NUM_; ++j) {
          const size_t idx = QUAD_NUM_*i+j;
          Hr = zeros<double>(8, 8);
          quad_linear_at_quadr_hes_(&Hr[0], &vert[0], &H_invDmH_(0, idx), &detDmH_[idx], &qw_[j], &lame_(0, i), &lame_(1, i));
          H += Hr;
        }
        break;
      }
      case COROTATED: {
        for (size_t j = 0; j < QUAD_NUM_; ++j) {
          const size_t idx = QUAD_NUM_*i+j;
          Hr = zeros<double>(8, 8);
          DF = vert*itr_matrix<const double *>(4, 2, &H_invDmH_(0, idx));
          decomp(DF);
          quad_coro_at_quadr_hes_(&Hr[0], &vert[0], &H_invDmH_(0, idx), &DF[0], &detDmH_[idx], &qw_[j], &lame_(0, i), &lame_(1, i));
          H += Hr;
        }
        break;
      }
      case STVK: {
        for (size_t j = 0; j < QUAD_NUM_; ++j) {
          const size_t idx = QUAD_NUM_*i+j;
          Hr = zeros<double>(8, 8);
          quad_stvk_at_quadr_hes_(&Hr[0], &vert[0], &H_invDmH_(0, idx), &detDmH_[idx], &qw_[j], &lame_(0, i), &lame_(1, i));
          H += Hr;
        }
        break;
      }
      case NEOHOOKEAN: {
        for (size_t j = 0; j < QUAD_NUM_; ++j) {
          const size_t idx = QUAD_NUM_*i+j;
          Hr = zeros<double>(8, 8);
          quad_neo_at_quadr_hes_(&Hr[0], &vert[0], &H_invDmH_(0, idx), &detDmH_[idx], &qw_[j], &lame_(0, i), &lame_(1, i));
          H += Hr;
        }
        break;
      }
      default:
        break;
    }
    #pragma omp critical
    {
      for (size_t p = 0; p < 8; ++p) {
        for (size_t q = 0; q < 8; ++q) {
          const size_t I = 2*quad_(p/2, i)+p%2;
          const size_t J = 2*quad_(q/2, i)+q%2;
          hes->push_back(Triplet<double>(I, J, w_*H(p, q)));
        }
      }
    }
  }
  return 0;
}

void quad4_elastic_energy::update_mtr_params(const matd_t &new_params) {
  ASSERT(new_params.size(1) == lame_.size(1) && new_params.size(2) == lame_.size(2));
  lame_ = new_params;
}
//===============================================================================
quad4_admm_energy::quad4_admm_energy(const mati_t &quad4, const matd_t &nods, Material type,
                                     const matd_t &lame, const double w, const char hessian)
    : dim_(nods.size()), w_(w), type_(type), lame_(lame), quad4_(quad4), hessian_(hessian) {
  ASSERT(quad4_.size(2) == lame_.size(2));

  //-> extract quad9
  get_quad9_elem(quad4_, quad9_);

  //-> select two quadratures and calc weights
  const double qrs[2] = {-1.0/sqrt(3.0), +1.0/sqrt(3.0)};
  const double qws[2] = {1.0, 1.0};
  QUAD_NUM_ = 2*2;
  quadr_ = zeros<double>(2, QUAD_NUM_);
  qw_ = ones<double>(QUAD_NUM_, 1);
  for (size_t i = 0; i < 2; ++i) {
    for (size_t j = 0; j < 2; ++j) {
      const size_t idx = 2*i+j;
      quadr_(0, idx) = qrs[i];
      quadr_(1, idx) = qrs[j];
      qw_[idx] = qws[i]*qws[j];
    }
  }

  dets_.resize(1, QUAD_NUM_*quad4.size(2));
  H_invDmH9_.resize(2*9, QUAD_NUM_*quad4_.size(2));
  
  #pragma omp parallel for
  for (size_t i = 0; i < quad4_.size(2); ++i) {
    const matd_t Y = nods(colon(), quad9_(colon(), i/4));
    for (size_t j = 0; j < QUAD_NUM_; ++j) { // for quadratic basis
      const size_t idx = QUAD_NUM_*i+j;

      double xi[2] = {0, 0};
      switch ( i%4 ) {
        case 0:
          xi[0] = (quadr_(0, j)-1)/2;
          xi[1] = (quadr_(1, j)-1)/2;
          break;
        case 1:
          xi[0] = (quadr_(0, j)+1)/2;
          xi[1] = (quadr_(1, j)-1)/2;
          break;
        case 2:
          xi[0] = (quadr_(0, j)+1)/2;
          xi[1] = (quadr_(1, j)+1)/2;
          break;
        case 3:
          xi[0] = (quadr_(0, j)-1)/2;
          xi[1] = (quadr_(1, j)+1)/2;
          break;
      }
      matd_t H = zeros<double>(9, 2);
      quad9_SF_jac(&H[0], xi);
      
      matd_t DmH = Y*H;
      dets_[idx] = fabs(DmH(0,0)*DmH(1,1)-DmH(0,1)*DmH(1,0))/4;

      if ( inv(DmH) ) cerr << "# inv fail" << endl;
      const matd_t HDmH = H*DmH;
      H_invDmH9_(colon(), idx) = HDmH(colon());
    }
  }

  //-> determine if the quad9 vertex is the one of the current quad4 vertex
  in_quad_ = zeros<int>(9, quad4_.size(2));
  #pragma omp parallel for
  for (size_t i = 0; i < quad4_.size(2); ++i) {
    const mati_t q4 = quad4_(colon(), i);
    const mati_t q9 = quad9_(colon(), i/4);
    for (size_t j = 0; j < q9.size(); ++j) {
      const auto it = std::find(q4.begin(), q4.end(), q9[j]);
      if ( it != q4.end() )
        in_quad_(j, i) = it-q4.begin()+1;
    }
    ASSERT(sum(in_quad_(colon(), i)) == 10);
  }

  exactH_ = truncH_ = zeros<double>(18*18, quad4_.size(2));
}

size_t quad4_admm_energy::Nx() const {
  return dim_;
}

int quad4_admm_energy::Val(const double *x, double *val) const {
  RETURN_WITH_COND_TRUE(w_ == 0.0);
  itr_matrix<const double *> X(2, dim_/2, x);

  #pragma omp parallel for
  for (size_t i = 0; i < quad4_.size(2); ++i) {
    hj::polar2d decomp;
    matd_t DF(2, 2);
    const matd_t vert4 = X(colon(), quad4_(colon(), i));
    const matd_t vert9 = X(colon(), quad9_(colon(), i/4));

    double value = 0;
    switch ( type_ ) {
      case LINEAR: {
        break;
      }
      case COROTATED: {
        break;
      }
      case STVK: {
        for (size_t j = 0; j < QUAD_NUM_; ++j) {
          const size_t idx = QUAD_NUM_*i+j;
          double vr = 0;
          quad9_stvk_at_quadr_(&vr, &vert9[0], &H_invDmH9_(0, idx), &dets_[idx], &qw_[j], &lame_(0, i), &lame_(1, i));
          value += vr;
        }
        break;
      }
      case NEOHOOKEAN: {
        for (size_t j = 0; j < QUAD_NUM_; ++j) {
          const size_t idx = QUAD_NUM_*i+j;
          double vr = 0;
          quad9_neo_at_quadr_(&vr, &vert9[0], &H_invDmH9_(0, idx), &dets_[idx], &qw_[j], &lame_(0, i), &lame_(1, i));
          value += vr;
        }
        break;
      }
      default:
        break;
    }
    #pragma omp critical
    {
      *val += w_*value;
    }
  }
  return 0;
}

int quad4_admm_energy::Gra(const double *x, double *gra) const {
  RETURN_WITH_COND_TRUE(w_ == 0.0);
  itr_matrix<const double *> X(2, dim_/2, x);
  itr_matrix<double *> G(2, dim_/2, gra);
  
  #pragma omp parallel for
  for (size_t i = 0; i < quad4_.size(2); ++i) {
    hj::polar2d decomp;
    matd_t DF(2, 2);
    const matd_t vert4 = X(colon(), quad4_(colon(), i));
    const matd_t vert9 = X(colon(), quad9_(colon(), i/4));

    matd_t g = zeros<double>(2, 9);
    switch ( type_ ) {
      case LINEAR: {
        break;
      }
      case COROTATED: {
        break;
      }
      case STVK: {
        for (size_t j = 0; j < QUAD_NUM_; ++j) {
          const size_t idx = QUAD_NUM_*i+j;
          matd_t gr = zeros<double>(2, 9);
          quad9_stvk_at_quadr_jac_(&gr[0], &vert9[0], &H_invDmH9_(0, idx), &dets_[idx], &qw_[j], &lame_(0, i), &lame_(1, i));
          g += gr;
        }
        break;
      }
      case NEOHOOKEAN: {
        for (size_t j = 0; j < QUAD_NUM_; ++j) {
          const size_t idx = QUAD_NUM_*i+j;
          matd_t gr = zeros<double>(2, 9);
          quad9_neo_at_quadr_jac_(&gr[0], &vert9[0], &H_invDmH9_(0, idx), &dets_[idx], &qw_[j], &lame_(0, i), &lame_(1, i));
          g += gr;
        }
        break;
      }
      default:
        break;
    }
    #pragma omp critical
    {
      G(colon(), quad9_(colon(), i/4)) += w_*g;
    }
  }
  return 0;
}

void quad4_admm_energy::set_hess_flag(const char flag) {
  ASSERT(flag == 'N' || flag == 'T');
  hessian_ = flag;
}

int quad4_admm_energy::Hes(const double *x, vector<Triplet<double>> *hes) const {
  RETURN_WITH_COND_TRUE(w_ == 0.0);

  if ( hessian_ == 'N' )
    return HesExactNewton(x, hes);

  if ( hessian_ == 'T' )
    return HesTruncNewton(x, hes);

  return __LINE__;
}

int quad4_admm_energy::HesExactNewton(const double *x, vector<Triplet<double>> *hes) const {
  itr_matrix<const double *> X(2, dim_/2, x);

  #pragma omp parallel for
  for (size_t i = 0; i < quad4_.size(2); ++i) {
    hj::polar2d decomp;
    matd_t DF(2,2);
    const matd_t vert4 = X(colon(), quad4_(colon(), i));
    const matd_t vert9 = X(colon(), quad9_(colon(), i/4));

    matd_t H = zeros<double>(18, 18);
    switch ( type_ ) {
      case LINEAR: {
        break;
      }
      case COROTATED: {
        break;
      }
      case STVK: {
        for (size_t j = 0; j < QUAD_NUM_; ++j) {
          const size_t idx = QUAD_NUM_*i+j;
          matd_t Hr = zeros<double>(18, 18);
          quad9_stvk_at_quadr_hes_(&Hr[0], &vert9[0], &H_invDmH9_(0, idx), &dets_[idx], &qw_[j], &lame_(0, i), &lame_(1, i));
          H += Hr;
        }
        break;
      }
      case NEOHOOKEAN: {
        for (size_t j = 0; j < QUAD_NUM_; ++j) {
          const size_t idx = QUAD_NUM_*i+j;
          matd_t Hr = zeros<double>(18, 18);
          quad9_neo_at_quadr_hes_(&Hr[0], &vert9[0], &H_invDmH9_(0, idx), &dets_[idx], &qw_[j], &lame_(0, i), &lame_(1, i));
          H += Hr;
        }
        break;
      }
      default:
        break;
    }
    #pragma omp critical
    {
      for (size_t p = 0; p < 18; ++p) {
        for (size_t q = 0; q < 18; ++q) {
          const size_t I = 2*quad9_(p/2, i/4)+p%2;
          const size_t J = 2*quad9_(q/2, i/4)+q%2;
          hes->push_back(Triplet<double>(I, J, w_*H(p, q)));
        }
      }
    }
  }
  return 0;
}

int quad4_admm_energy::HesTruncNewton(const double *x, vector<Triplet<double>> *hes) const {
  itr_matrix<const double *> X(2, dim_/2, x);

  #pragma omp parallel for
  for (size_t i = 0; i < quad4_.size(2); ++i) {
    hj::polar2d decomp;
    matd_t DF(2,2);
    const matd_t vert4 = X(colon(), quad4_(colon(), i));
    const matd_t vert9 = X(colon(), quad9_(colon(), i/4));

    matd_t H = zeros<double>(18, 18);
    switch ( type_ ) {
      case LINEAR: {
        break;
      }
      case COROTATED: {
        break;
      }
      case STVK: {
        for (size_t j = 0; j < QUAD_NUM_; ++j) {
          const size_t idx = QUAD_NUM_*i+j;
          matd_t Hr = zeros<double>(18, 18);
          quad9_stvk_at_quadr_hes_(&Hr[0], &vert9[0], &H_invDmH9_(0, idx), &dets_[idx], &qw_[j], &lame_(0, i), &lame_(1, i));
          H += Hr;
        }
        break;
      }
      case NEOHOOKEAN: {
        for (size_t j = 0; j < QUAD_NUM_; ++j) {
          const size_t idx = QUAD_NUM_*i+j;
          matd_t Hr = zeros<double>(18, 18);
          quad9_neo_at_quadr_hes_(&Hr[0], &vert9[0], &H_invDmH9_(0, idx), &dets_[idx], &qw_[j], &lame_(0, i), &lame_(1, i));
          H += Hr;
        }
        break;
      }
      default:
        break;
    }

    if ( debug_flag_ ) {
      matd_t Ht = H;
      for (size_t p = 0; p < 18; ++p) {
        for (size_t q = 0; q < 18; ++q) {
          if ( in_quad_(p/2, i) && in_quad_(q/2, i) )
            continue;
          Ht(p, q) = 0;
        }
      }
      exactH_(colon(), i) = H(colon());
      truncH_(colon(), i) = Ht(colon());
    }
    
    #pragma omp critical
    {
      for (size_t p = 0; p < 18; ++p) {
        for (size_t q = 0; q < 18; ++q) {
          const size_t I = 2*quad9_(p/2, i/4)+p%2;
          const size_t J = 2*quad9_(q/2, i/4)+q%2;
          if ( in_quad_(p/2, i) && in_quad_(q/2, i) )
            hes->push_back(Triplet<double>(I, J, w_*H(p, q)));
        }
      }
    }
  }
  return 0;
}
//==============================================================================
quad9_elastic_energy::quad9_elastic_energy(const mati_t &quad9, const matd_t &nods, Material type,
                                           const matd_t &lame, const double w)
    : dim_(nods.size()), w_(w), type_(type), lame_(lame), quad9_(quad9) {
  ASSERT(4*quad9.size(2) == lame.size(2));

  QUAD_NUM_ = 16;
  quadr_ = zeros<double>(2, QUAD_NUM_);
  qw_ = zeros<double>(QUAD_NUM_, 1);

  //-> generate quadrature points and its weights
  const double Q1 = sqrt(3.0/7-2.0/7*sqrt(6.0/5));
  const double Q2 = sqrt(3.0/7+2.0/7*sqrt(6.0/5));
  const double w1 = (18+sqrt(30.0))/36;
  const double w2 = (18-sqrt(30.0))/36;
  const double Qr[4] = {-Q2, -Q1, +Q1, +Q2};
  const double wr[4] = {w2, w1, w1, w2};
  for (size_t i = 0; i < 4; ++i) {
    for (size_t j = 0; j < 4; ++j) {
      quadr_(0, 4*i+j) = Qr[i];
      quadr_(1, 4*i+j) = Qr[j];
      qw_[4*i+j] = wr[i]*wr[j];
    }
  }

  detDmH_.resize(1, QUAD_NUM_*quad9_.size(2));
  H_invDmH_.resize(18, QUAD_NUM_*quad9_.size(2));

  #pragma omp parallel for
  for (size_t i = 0; i < quad9_.size(2); ++i) {
    const matd_t X = nods(colon(), quad9_(colon(), i));

    for (size_t j = 0; j < QUAD_NUM_; ++j) { // for four quadratures
      const size_t idx = QUAD_NUM_*i+j;
      matd_t H = zeros<double>(9, 2);
      quad9_SF_jac(&H[0], &quadr_(0, j));

      matd_t DmH = X*H;
      matd_t cpDmH = DmH;

      detDmH_[idx] = fabs(det(cpDmH));

      if ( inv(DmH) )
        cerr << "# inv fail" << endl;
      matd_t HDmH = H*DmH;
      H_invDmH_(colon(), idx) = HDmH(colon());
    }
  }
}

size_t quad9_elastic_energy::Nx() const {
  return dim_;
}

int quad9_elastic_energy::Val(const double *x, double *val) const {
  RETURN_WITH_COND_TRUE(w_ == 0.0);
  itr_matrix<const double *> X(2, dim_/2, x);

  #pragma omp parallel for
  for (size_t i = 0; i < quad9_.size(2); ++i) {
    hj::polar2d decomp;
    matd_t vert(2, 9), DF(2, 2);
    vert = X(colon(), quad9_(colon(), i));
    double value = 0;
    switch ( type_ ) {
      case LINEAR: {
        for (size_t j = 0; j < QUAD_NUM_; ++j) {
          const size_t idx = QUAD_NUM_*i+j;
          double vr = 0;
          quad9_linear_at_quadr_(&vr, &vert[0], &H_invDmH_(0, idx), &detDmH_[idx], &qw_[j], &lame_(0, i), &lame_(1, i));
          value += vr;
        }
        break;
      }
      case COROTATED: {
        for (size_t j = 0; j < QUAD_NUM_; ++j) {
          const size_t idx = QUAD_NUM_*i+j;
          double vr = 0;
          DF = vert*itr_matrix<const double *>(9, 2, &H_invDmH_(0, idx));
          decomp(DF);
          quad9_coro_at_quadr_(&vr, &vert[0], &H_invDmH_(0, idx), &DF[0], &detDmH_[idx], &qw_[j], &lame_(0, i), &lame_(1, i));
          value += vr;
        }
        break;
      }
      case STVK: {
        for (size_t j = 0; j < QUAD_NUM_; ++j) {
          const size_t idx = QUAD_NUM_*i+j;
          double vr = 0;
          quad9_stvk_at_quadr_(&vr, &vert[0], &H_invDmH_(0, idx), &detDmH_[idx], &qw_[j], &lame_(0, i), &lame_(1, i));
          value += vr;
        }
        break;
      }
      case NEOHOOKEAN: {
        for (size_t j = 0; j < QUAD_NUM_; ++j) {
          const size_t idx = QUAD_NUM_*i+j;
          double vr = 0;
          quad9_neo_at_quadr_(&vr, &vert[0], &H_invDmH_(0, idx), &detDmH_[idx], &qw_[j], &lame_(0, i), &lame_(1, i));
          value += vr;
        }
        break;
      }
      default:
        break;
    }
    #pragma omp critical
    {
      *val += w_*value;
    }
  }
  return 0;
}

int quad9_elastic_energy::Gra(const double *x, double *gra) const {
  RETURN_WITH_COND_TRUE(w_ == 0.0);
  itr_matrix<const double *> X(2, dim_/2, x);
  itr_matrix<double *> G(2, dim_/2, gra);

  #pragma omp parallel for
  for (size_t i = 0; i < quad9_.size(2); ++i) {
    hj::polar2d decomp;
    matd_t g(2, 9), gr(2, 9), DF(2, 2);
    const matd_t vert = X(colon(), quad9_(colon(), i));
    g = zeros<double>(2, 9);
    switch ( type_ ) {
      case LINEAR: {
        for (size_t j = 0; j < QUAD_NUM_; ++j) {
          const size_t idx = QUAD_NUM_*i+j;
          gr = zeros<double>(2, 9);
          quad9_linear_at_quadr_jac_(&gr[0], &vert[0], &H_invDmH_(0, idx), &detDmH_[idx], &qw_[j], &lame_(0, i), &lame_(1, i));
          g += gr;
        }
        break;
      }
      case COROTATED: {
        for (size_t j = 0; j < QUAD_NUM_; ++j) {
          const size_t idx = QUAD_NUM_*i+j;
          gr = zeros<double>(2, 9);
          DF = vert*itr_matrix<const double *>(9, 2, &H_invDmH_(0, idx));
          decomp(DF);
          quad9_coro_at_quadr_jac_(&gr[0], &vert[0], &H_invDmH_(0, idx), &DF[0], &detDmH_[idx], &qw_[j], &lame_(0, i), &lame_(1, i));
          g += gr;
        }
        break;
      }
      case STVK: {
        for (size_t j = 0; j < QUAD_NUM_; ++j) {
          const size_t idx = QUAD_NUM_*i+j;
          gr = zeros<double>(2, 9);
          quad9_stvk_at_quadr_jac_(&gr[0], &vert[0], &H_invDmH_(0, idx), &detDmH_[idx], &qw_[j], &lame_(0, i), &lame_(1, i));
          g += gr;
        }
        break;
      }
      case NEOHOOKEAN: {
        for (size_t j = 0; j < QUAD_NUM_; ++j) {
          const size_t idx = QUAD_NUM_*i+j;
          gr = zeros<double>(2, 9);
          quad9_neo_at_quadr_jac_(&gr[0], &vert[0], &H_invDmH_(0, idx), &detDmH_[idx], &qw_[j], &lame_(0, i), &lame_(1, i));
          g += gr;
        }
        break;
      }
      default:
        break;
    }
    #pragma omp critical
    {
      G(colon(), quad9_(colon(), i)) += w_*g;
    }
  }
  return 0;
}

int quad9_elastic_energy::Hes(const double *x, vector<Triplet<double>> *hes) const {
  RETURN_WITH_COND_TRUE(w_ == 0.0);
  itr_matrix<const double *> X(2, dim_/2, x);

  #pragma omp parallel for
  for (size_t i = 0; i < quad9_.size(2); ++i) {
    hj::polar2d decomp;
    matd_t H(18, 18), Hr(18, 18), DF(2, 2);
    const matd_t vert = X(colon(), quad9_(colon(), i));
    H = zeros<double>(18, 18);
    switch ( type_ ) {
      case LINEAR: {
        for (size_t j = 0; j < QUAD_NUM_; ++j) {
          const size_t idx = QUAD_NUM_*i+j;
          Hr = zeros<double>(18, 18);
          quad9_linear_at_quadr_hes_(&Hr[0], &vert[0], &H_invDmH_(0, idx), &detDmH_[idx], &qw_[j], &lame_(0, i), &lame_(1, i));
          H += Hr;
        }
        break;
      }
      case COROTATED: {
        for (size_t j = 0; j < QUAD_NUM_; ++j) {
          const size_t idx = QUAD_NUM_*i+j;
          Hr = zeros<double>(18, 18);
          DF = vert*itr_matrix<const double *>(9, 2, &H_invDmH_(0, idx));
          decomp(DF);
          quad9_coro_at_quadr_hes_(&Hr[0], &vert[0], &H_invDmH_(0, idx), &DF[0], &detDmH_[idx], &qw_[j], &lame_(0, i), &lame_(1, i));
          H += Hr;
        }
        break;
      }
      case STVK: {
        for (size_t j = 0; j < QUAD_NUM_; ++j) {
          const size_t idx = QUAD_NUM_*i+j;
          Hr = zeros<double>(18, 18);
          quad9_stvk_at_quadr_hes_(&Hr[0], &vert[0], &H_invDmH_(0, idx), &detDmH_[idx], &qw_[j], &lame_(0, i), &lame_(1, i));
          H += Hr;
        }
        break;
      }
      case NEOHOOKEAN: {
        for (size_t j = 0; j < QUAD_NUM_; ++j) {
          const size_t idx = QUAD_NUM_*i+j;
          Hr = zeros<double>(18, 18);
          quad9_neo_at_quadr_hes_(&Hr[0], &vert[0], &H_invDmH_(0, idx), &detDmH_[idx], &qw_[j], &lame_(0, i), &lame_(1, i));
          H += Hr;
        }
        break;
      }
      default:
        break;
    }
    #pragma omp critical
    {
      for (size_t p = 0; p < 18; ++p) {
        for (size_t q = 0; q < 18; ++q) {
          const size_t I = 2*quad9_(p/2, i)+p%2;
          const size_t J = 2*quad9_(q/2, i)+q%2;
          hes->push_back(Triplet<double>(I, J, w_*H(p, q)));
        }
      }
    }
  }
  return 0;
}

void quad9_elastic_energy::update_mtr_params(const matd_t &new_params) {
  
}
//===============================================================================
voxel_elastic_potential::voxel_elastic_potential(const mati_t &cube, const matd_t &nods, Material type, const matd_t &lame, const double w)
    : dim_(nods.size()), w_(w), type_(type), lame_(lame), cube_(cube) {
  ASSERT(cube_.size(2) == lame_.size(2));

  elemE_ = zeros<double>(cube_.size(2), 1);
  
  //-> select quadratures and calc weights
  const double qrs[2] = {-1.0/sqrt(3.0), +1.0/sqrt(3.0)};
  const double qws[2] = {1.0, 1.0};
  const size_t QUAD_NUM = 2*2*2;
  quadr_ = zeros<double>(3, QUAD_NUM);
  qw_ = ones<double>(QUAD_NUM, 1);
  for (size_t i = 0; i < 2; ++i) {
    for (size_t j = 0; j < 2; ++j) {
      for (size_t k = 0; k < 2; ++k) {
        const size_t idx = 4*i+2*j+k;
        quadr_(0, idx) = qrs[i];
        quadr_(1, idx) = qrs[j];
        quadr_(2, idx) = qrs[k];
        qw_[idx] = qws[i]*qws[j]*qws[k];
      }
    }
  }
  
  detDmH_.resize(1, QUAD_NUM*cube_.size(2));
  H_invDmH_.resize(24, QUAD_NUM*cube_.size(2));

  #pragma omp parallel for
  for (size_t i = 0; i < cube_.size(2); ++i) {
    const matd_t X = nods(colon(), cube_(colon(), i));

    for (size_t j = 0; j < QUAD_NUM; ++j) { // for eight quadratures
      const size_t idx = QUAD_NUM*i+j;
      matd_t H = zeros<double>(8, 3);
      vox_SF_jac(&H[0], &quadr_(0, j));

      matd_t DmH = X*H;
      matd_t cpDmH = DmH;

      detDmH_[idx] = fabs(det(cpDmH));

      if ( inv(DmH) )
        cerr << "# inv fail" << endl;
      matd_t HDmH = H*DmH;
      H_invDmH_(colon(), idx) = HDmH(colon());
    }
  }
}

size_t voxel_elastic_potential::Nx() const {
  return dim_;
}

int voxel_elastic_potential::Val(const double *x, double *val) const {
  RETURN_WITH_COND_TRUE(w_ == 0.0);
  itr_matrix<const double *> X(3, dim_/3, x);

  hj::polar3d decomp;
  matd_t vert(3, 8), DF(3, 3);
  
  for (size_t i = 0; i < cube_.size(2); ++i) {
    vert = X(colon(), cube_(colon(), i));
    double value = 0;
    
    switch ( type_ ) {
      case LINEAR: {
        for (size_t j = 0; j < 8; ++j) {
          const size_t idx = 8*i+j;
          double vr = 0;
          vox_linear_at_quadr_(&vr, &vert[0], &H_invDmH_(0, idx), &detDmH_[idx], &qw_[j], &lame_(0, i), &lame_(1, i));
          value += vr;
        }
        break;
      }
      case COROTATED: {
        for (size_t j = 0; j < 8; ++j) {
          const size_t idx = 8*i+j;
          double vr = 0;
          DF = vert*itr_matrix<const double *>(8, 3, &H_invDmH_(0, idx));
          decomp(DF);
          vox_corotated_at_quadr_(&vr, &vert[0], &H_invDmH_(0, idx), &DF[0], &detDmH_[idx], &qw_[j], &lame_(0, i), &lame_(1, i));
          value += vr;
        }
        break;
      }
      case STVK: {
        for (size_t j = 0; j < 8; ++j) {
          const size_t idx = 8*i+j;
          double vr = 0;
          vox_stvk_at_quadr_(&vr, &vert[0], &H_invDmH_(0, idx), &detDmH_[idx], &qw_[j], &lame_(0, i), &lame_(1, i));
          value += vr;
        }
        break;
      }
      case NEOHOOKEAN: {
        for (size_t j = 0; j < 8; ++j) {
          const size_t idx = 8*i+j;
          double vr = 0;
          vox_neo_at_quadr_(&vr, &vert[0], &H_invDmH_(0, idx), &detDmH_[idx], &qw_[j], &lame_(0, i), &lame_(1, i));
          value += vr;
        }
        break;
      }
      case STANEO: {
        for (size_t j = 0; j < 8; ++j) {
          const size_t idx = 8*i+j;
          double vr = 0;
          vox_sta_neo_at_quadr_(&vr, &vert[0], &H_invDmH_(0, idx), &detDmH_[idx], &qw_[j], &lame_(0, i), &lame_(1, i));
          value += vr;
        }
        break;
      }
      case BOWERNEO: {
        for (size_t j = 0; j < 8; ++j) {
          const size_t idx = 8*i+j;
          double vr = 0;
          vox_bower_neo_at_quadr_(&vr, &vert[0], &H_invDmH_(0, idx), &detDmH_[idx], &qw_[j], &lame_(0, i), &lame_(1, i));
          value += vr;
        }
        break;
      }
      default:
        break;
    }
    elemE_[i] = w_*value;
    {
      *val += w_*value;
    }
  }
    
  return 0;
}

int voxel_elastic_potential::Gra(const double *x, double *gra) const {
  RETURN_WITH_COND_TRUE(w_ == 0.0);
  itr_matrix<const double *> X(3, dim_/3, x);
  itr_matrix<double *> G(3, dim_/3, gra);

  matd_t vert(3, 8), g(3, 8), gr(3, 8), DF(3, 3);
  hj::polar3d decomp;
  
  for (size_t i = 0; i < cube_.size(2); ++i) {
    vert = X(colon(), cube_(colon(), i));
    g = zeros<double>(3, 8);
    switch ( type_ ) {
      case LINEAR: {
        for (size_t j = 0; j < 8; ++j) {
          const size_t idx = 8*i+j;
          gr = zeros<double>(3, 8);
          vox_linear_at_quadr_jac_(&gr[0], &vert[0], &H_invDmH_(0, idx), &detDmH_[idx], &qw_[j], &lame_(0, i), &lame_(1, i));
          g += gr;
        }
        break;
      }
      case COROTATED: {
        for (size_t j = 0; j < 8; ++j) {
          const size_t idx = 8*i+j;
          gr = zeros<double>(3, 8);
          DF = vert*itr_matrix<const double *>(8, 3, &H_invDmH_(0, idx));
          decomp(DF);
          vox_corotated_at_quadr_jac_(&gr[0], &vert[0], &H_invDmH_(0, idx), &DF[0], &detDmH_[idx], &qw_[j], &lame_(0, i), &lame_(1, i));
          g += gr;
        }
        break;
      }
      case STVK: {
        for (size_t j = 0; j < 8; ++j) {
          const size_t idx = 8*i+j;
          gr = zeros<double>(3, 8);
          vox_stvk_at_quadr_jac_(&gr[0], &vert[0], &H_invDmH_(0, idx), &detDmH_[idx], &qw_[j], &lame_(0, i), &lame_(1, i));
          g += gr;
        }
        break;
      }
      case NEOHOOKEAN: {
        for (size_t j = 0; j < 8; ++j) {
          const size_t idx = 8*i+j;
          gr = zeros<double>(3, 8);
          vox_neo_at_quadr_jac_(&gr[0], &vert[0], &H_invDmH_(0, idx), &detDmH_[idx], &qw_[j], &lame_(0, i), &lame_(1, i));
          g += gr;
        }
        break;
      }
      case STANEO: {
        for (size_t j = 0; j < 8; ++j) {
          const size_t idx = 8*i+j;
          gr = zeros<double>(3, 8);
          vox_sta_neo_at_quadr_jac_(&gr[0], &vert[0], &H_invDmH_(0, idx), &detDmH_[idx], &qw_[j], &lame_(0, i), &lame_(1, i));
          g += gr;
        }
        break;
      }
      case BOWERNEO: {
        for (size_t j = 0; j < 8; ++j) {
          const size_t idx = 8*i+j;
          gr = zeros<double>(3, 8);
          vox_bower_neo_at_quadr_jac_(&gr[0], &vert[0], &H_invDmH_(0, idx), &detDmH_[idx], &qw_[j], &lame_(0, i), &lame_(1, i));
          g += gr;
        }
        break;
      }
      default:
        break;
    }
    {
      G(colon(), cube_(colon(), i)) += w_*g;
    }
  }
  return 0;
}

int voxel_elastic_potential::Hes(const double *x, vector<Triplet<double>> *hes) const {
  RETURN_WITH_COND_TRUE(w_ == 0.0);
  itr_matrix<const double *> X(3, dim_/3, x);

  matd_t vert(3, 8), H(24, 24), Hr(24, 24), DF(3, 3);
  hj::polar3d decomp;
  
  for (size_t i = 0; i < cube_.size(2); ++i) {
    vert = X(colon(), cube_(colon(), i));
    H = zeros<double>(24, 24);
    switch ( type_ ) {
      case LINEAR: {
        for (size_t j = 0; j < 8; ++j) {
          const size_t idx = 8*i+j;
          Hr = zeros<double>(24, 24);
          vox_linear_at_quadr_hes_(&Hr[0], &vert[0], &H_invDmH_(0, idx), &detDmH_[idx], &qw_[j], &lame_(0, i), &lame_(1, i));
          H += Hr;
        }
        break;
      }
      case COROTATED: {
        for (size_t j = 0; j < 8; ++j) {
          const size_t idx = 8*i+j;
          Hr = zeros<double>(24, 24);
          DF = vert*itr_matrix<const double *>(8, 3, &H_invDmH_(0, idx));
          decomp(DF);
          vox_corotated_at_quadr_hes_(&Hr[0], &vert[0], &H_invDmH_(0, idx), &DF[0], &detDmH_[idx], &qw_[j], &lame_(0, i), &lame_(1, i));
          H += Hr;
        }
        break;
      }
      case STVK: {
        for (size_t j = 0; j < 8; ++j) {
          const size_t idx = 8*i+j;
          Hr = zeros<double>(24, 24);
          vox_stvk_at_quadr_hes_(&Hr[0], &vert[0], &H_invDmH_(0, idx), &detDmH_[idx], &qw_[j], &lame_(0, i), &lame_(1, i));
          H += Hr;
        }
        break;
      }
      case NEOHOOKEAN: {
        for (size_t j = 0; j < 8; ++j) {
          const size_t idx = 8*i+j;
          Hr = zeros<double>(24, 24);
          vox_neo_at_quadr_hes_(&Hr[0], &vert[0], &H_invDmH_(0, idx), &detDmH_[idx], &qw_[j], &lame_(0, i), &lame_(1, i));
          H += Hr;
        }
        break;
      }
      case STANEO: {
        for (size_t j = 0; j < 8; ++j) {
          const size_t idx = 8*i+j;
          Hr = zeros<double>(24, 24);
          vox_sta_neo_at_quadr_hes_(&Hr[0], &vert[0], &H_invDmH_(0, idx), &detDmH_[idx], &qw_[j], &lame_(0, i), &lame_(1, i));
          H += Hr;
        }
        break;
      }
      case BOWERNEO: {
        for (size_t j = 0; j < 8; ++j) {
          const size_t idx = 8*i+j;
          Hr = zeros<double>(24, 24);
          vox_bower_neo_at_quadr_hes_(&Hr[0], &vert[0], &H_invDmH_(0, idx), &detDmH_[idx], &qw_[j], &lame_(0, i), &lame_(1, i));
          H += Hr;
        }
        break;
      }
      default:
        break;
    }
    {
      for (size_t p = 0; p < 24; ++p) {
        for (size_t q = 0; q < 24; ++q) {
          const size_t I = 3*cube_(p/3, i)+p%3;
          const size_t J = 3*cube_(q/3, i)+q%3;
          if ( fabs(H(p, q)) != 0.0 )
            hes->push_back(Triplet<double>(I, J, w_*H(p, q)));
        }
      }
    }
  }

  return 0;
}

void voxel_elastic_potential::update_mtr_params(const matd_t &new_params) {
  ASSERT(new_params.size(1) == lame_.size(1) && new_params.size(2) == lame_.size(2));
  lame_ = new_params;
}
//==============================================================================
positional_potential::positional_potential(const matd_t &nods, const double w)
    : dim_(nods.size()), rd_(nods.size(1)), w_(w) {}

size_t positional_potential::Nx() const {
  return dim_;
}

int positional_potential::Val(const double *x, double *val) const {
  Eigen::Map<const MatrixXd> X(x, rd_, dim_/rd_);
  
  if ( rd_ == 3 ) {
    RETURN_WITH_COND_TRUE(w_ == 0.0 || fixed3d_.empty());
    for (auto &elem : fixed3d_) {
      const size_t id = elem.first;
      *val += 0.5*w_*(X.col(id)-elem.second).squaredNorm();
    }
    return 0;
  }

  if ( rd_ == 2 ) {
    RETURN_WITH_COND_TRUE(w_ == 0.0 || fixed2d_.empty());
    for (auto &elem : fixed2d_) {
      const size_t id = elem.first;
      *val += 0.5*w_*(X.col(id)-elem.second).squaredNorm();
    }
    return 0;
  }

  return __LINE__;
}

int positional_potential::Gra(const double *x, double *gra) const {
  Eigen::Map<const MatrixXd> X(x, rd_, dim_/rd_);
  Eigen::Map<MatrixXd> G(gra, rd_, dim_/rd_);
  
  if ( rd_ == 3 ) {
    RETURN_WITH_COND_TRUE(w_ == 0.0 || fixed3d_.empty());
    for (auto &elem : fixed3d_) {
      const size_t id = elem.first;
      G.col(id) += w_*(X.col(id)-elem.second);
    }
    return 0;
  }

  if ( rd_ == 2 ) {
    RETURN_WITH_COND_TRUE(w_ == 0.0 || fixed2d_.empty());
    for (auto &elem : fixed2d_) {
      const size_t id = elem.first;
      G.col(id) += w_*(X.col(id)-elem.second);
    }
    return 0;
  }
  
  return __LINE__;
}

int positional_potential::Hes(const double *x, vector<Triplet<double>> *hes) const {
  if ( rd_ == 3 ) {
    RETURN_WITH_COND_TRUE(w_ == 0.0 || fixed3d_.empty());
    for (auto &elem : fixed3d_) {
      const size_t id = elem.first;
      hes->push_back(Triplet<double>(3*id+0, 3*id+0, w_));
      hes->push_back(Triplet<double>(3*id+1, 3*id+1, w_));
      hes->push_back(Triplet<double>(3*id+2, 3*id+2, w_));
    }
    return 0;
  }

  if ( rd_ == 2 ) {
    RETURN_WITH_COND_TRUE(w_ == 0.0 || fixed2d_.empty());
    for (auto &elem : fixed2d_) {
      const size_t id = elem.first;
      hes->push_back(Triplet<double>(2*id+0, 2*id+0, w_));
      hes->push_back(Triplet<double>(2*id+1, 2*id+1, w_));
    }
    return 0;
  }

  return __LINE__;
}

void positional_potential::update_fixed_verts(const vector<size_t> &pids, const double *pos) {
  if ( rd_ == 3 ) {
    fixed3d_.clear();
    for (auto &id : pids) {
      ASSERT(id >= 0 && id < Nx()/rd_);
      fixed3d_.insert(make_pair(id, Vector3d(pos[3*id+0], pos[3*id+1], pos[3*id+2])));
    }
    return;
  }

  if ( rd_ == 2 ) {
    fixed2d_.clear();
    for (auto &id : pids) {
      ASSERT(id >= 0 && id < Nx()/rd_);
      fixed2d_.insert(make_pair(id, Vector2d(pos[2*id+0], pos[2*id+1])));
    }
    return;
  }
}

int positional_potential::Pin(const size_t id, const double *pos) {
  if ( id < 0 || id >= Nx()/rd_ )
    return __LINE__;

  if ( rd_ == 3) {
    fixed3d_[id] = Vector3d(pos);
    return 0;
  }

  if ( rd_ == 2 ) {
    fixed2d_[id] = Vector2d(pos);
    return 0;
  }

  return __LINE__;
}

int positional_potential::Release(const size_t id) {
  if ( id < 0 || id >= Nx()/rd_ ) {
    cerr << "[info] vertex index is out of range\n";
    return __LINE__;
  }

  if ( rd_ == 3 ) {
    auto it = fixed3d_.find(id);
    if ( it == fixed3d_.end() ) {
      cerr << "[info] vertex " << id << " is not fixed\n";
      return __LINE__;
    }
    fixed3d_.erase(it);
    return 0;
  }

  if ( rd_ == 2 ) {
    auto it = fixed2d_.find(id);
    if ( it == fixed2d_.end() ) {
      cerr << "[info] vertex " << id << " is not fixed\n";
      return __LINE__;
    }
    fixed2d_.erase(it);
    return 0;
  }

  return __LINE__;
}

//==============================================================================
spring_potential::spring_potential(const mati_t &edge, const matd_t &nods, const double w, char option)
  : dim_(nods.size()), edge_(edge), option_(option) {
  len_.resize(edge_.size(2), 1);
#pragma omp parallel for
  for (size_t i = 0; i < edge_.size(2); ++i) {
    matd_t vert = nods(colon(), edge_(colon(), i));
    calc_edge_length_(&len_[i], &vert[0]);
  }
  w_ = w*ones<double>(edge_.size(2), 1);
}

size_t spring_potential::Nx() const {
  return dim_;
}

int spring_potential::Val(const double *x, double *val) const {
  RETURN_WITH_COND_TRUE(max(w_) == 0.0);
  itr_matrix<const double *> X(3, dim_/3, x);
  for (size_t i = 0; i < edge_.size(2); ++i) {
    matd_t vert = X(colon(), edge_(colon(), i));
    double value = 0;
    mass_spring_(&value, &vert[0], &len_[i]);
    *val += w_[i]*value;
  }
  return 0;
}

int spring_potential::Gra(const double *x, double *gra) const {
  RETURN_WITH_COND_TRUE(max(w_) == 0.0);
  itr_matrix<const double *> X(3, dim_/3, x);
  itr_matrix<double *> G(3, dim_/3, gra);
  for (size_t i = 0; i < edge_.size(2); ++i) {
    matd_t vert = X(colon(), edge_(colon(), i));
    matd_t g = zeros<double>(3, 2);
    mass_spring_jac_(&g[0], &vert[0], &len_[i]);
    G(colon(), edge_(colon(), i)) += w_[i]*g;
  }
  return 0;
}

void spring_potential::Newton(const matd_t &vert, const double &rest_len, matd_t &H) const {
  mass_spring_hes_(&H[0], &vert[0], &rest_len);
}

void spring_potential::GaussNewton(const matd_t &vert, const double &rest_len, matd_t &H) const {
  matd_t g = zeros<double>(6, 1);
  double energy = 0;
  mass_spring_(&energy, &vert[0], &rest_len);
  if ( energy < 1e-16 ) {
    calc_edge_length_jac_(&g[0], &vert[0]);
    g /= sqrt(rest_len);
  } else {
    mass_spring_jac_(&g[0], &vert[0], &rest_len);
    g /= 2*sqrt(energy);
  }
  H = 2*g*trans(g);
}

void spring_potential::Mixed(const matd_t &vert, const double &rest_len, matd_t &H) const {
  double curr_length = 0.0;
  calc_edge_length_(&curr_length, &vert[0]);
  if ( curr_length < rest_len )
    GaussNewton(vert, rest_len, H);
  else
    Newton(vert, rest_len, H);
}

int spring_potential::Hes(const double *x, vector<Triplet<double>> *hes) const {
  RETURN_WITH_COND_TRUE(max(w_) == 0.0);
  itr_matrix<const double *> X(3, dim_/3, x);
  for (size_t i = 0; i < edge_.size(2); ++i) {
    matd_t vert = X(colon(), edge_(colon(), i));
    matd_t H = zeros<double>(6, 6);
    switch ( option_ ) {
      case 'N': Newton(vert, len_[i], H); break;
      case 'G': GaussNewton(vert, len_[i], H); break;
      case 'M': Mixed(vert, len_[i], H); break;
    }
    for (size_t p = 0; p < 6; ++p) {
      for (size_t q = 0; q < 6; ++q) {
        const size_t I = 3*edge_(p/3, i)+p%3;
        const size_t J = 3*edge_(q/3, i)+q%3;
        if ( H(p, q) != 0.0 )
          hes->push_back(Triplet<double>(I, J, w_[i]*H(p, q)));
      }
    }
  }
  return 0;
}

int spring_potential::ResetEdgeMaterial(const size_t p, const size_t q, const double scale_w, const double scale_len) {
  size_t i = 0;
  for ( ; i < edge_.size(2); ++i) {
    if ( (p == edge_(0, i) && q == edge_(1, i)) || (p == edge_(1, i) && q == edge_(0, i)) ) {
      len_[i] *= scale_len;
      w_[i] *= scale_w;
      break;
    }
  }
  if ( i == edge_.size(2) ) {
    cerr << "[Error] No such edge\n";
    return __LINE__;
  }
  return 0;
}
//==============================================================================
stvk_spring_potential::stvk_spring_potential(const mati_t &edge, const matd_t &nods, const double w)
    : dim_(nods.size()), edge_(edge), w_(w) {
  len_.resize(edge_.size(2), 1);
  #pragma omp parallel for
  for (size_t i = 0; i < edge_.size(2); ++i) {
    matd_t vert = nods(colon(), edge_(colon(), i));
    calc_edge_length_(&len_[i], &vert[0]);
  }
}

size_t stvk_spring_potential::Nx() const {
  return dim_;
}

int stvk_spring_potential::Val(const double *x, double *val) const {
  RETURN_WITH_COND_TRUE(w_ == 0.0);
  itr_matrix<const double *> X(3, dim_/3, x);
  for (size_t i = 0; i < edge_.size(2); ++i) {
    const matd_t vert = X(colon(), edge_(colon(), i));
    double value = 0;
    poly_spring_(&value, &vert[0], &len_[i]);
    *val += w_*value;
  }
  return 0;
}

int stvk_spring_potential::Gra(const double *x, double *gra) const {
  RETURN_WITH_COND_TRUE(w_ == 0.0);
  itr_matrix<const double *> X(3, dim_/3, x);
  itr_matrix<double *> G(3, dim_/3, gra);
  for (size_t i = 0; i < edge_.size(2); ++i) {
    const matd_t vert = X(colon(), edge_(colon(), i));
    matd_t g = zeros<double>(3, 2);
    poly_spring_jac_(&g[0], &vert[0], &len_[i]);
    G(colon(), edge_(colon(), i)) += w_*g;
  }
  return 0;
}

int stvk_spring_potential::Hes(const double *x, vector<Triplet<double>> *hes) const {
  RETURN_WITH_COND_TRUE(w_ == 0.0);
  itr_matrix<const double *> X(3, dim_/3, x);
  for (size_t i = 0; i < edge_.size(2); ++i) {
    const matd_t vert = X(colon(), edge_(colon(), i));
    matd_t H = zeros<double>(6, 6);
    poly_spring_hes_(&H[0], &vert[0], &len_[i]);
    for (size_t p = 0; p < 6; ++p) {
      for (size_t q = 0; q < 6; ++q) {
        const size_t I = 3*edge_(p/3, i)+p%3;
        const size_t J = 3*edge_(q/3, i)+q%3;
        hes->push_back(Triplet<double>(I, J, w_*H(p, q)));
      }
    }
  }
  return 0;
}
//==============================================================================
gauss_newton_spring::gauss_newton_spring(const mati_t &edge, const matd_t &nods, const double w)
  : dim_(nods.size()), edge_(edge) {
  len_.resize(edge_.size(2), 1);
#pragma omp parallel for
  for (size_t i = 0; i < edge_.size(2); ++i) {
    matd_t vert = nods(colon(), edge_(colon(), i));
    calc_edge_length_(&len_[i], &vert[0]);
  }
  w_ = w*ones<double>(edge_.size(2), 1);
}

size_t gauss_newton_spring::Nx() const {
  return dim_;
}

int gauss_newton_spring::Val(const double *x, double *val) const {
  RETURN_WITH_COND_TRUE(max(w_) == 0.0);
  itr_matrix<const double *> X(3, dim_/3, x);
  for (size_t i = 0; i < edge_.size(2); ++i) {
    matd_t vert = X(colon(), edge_(colon(), i));
    double value = 0;
    mass_spring_(&value, &vert[0], &len_[i]);
    *val += w_[i]*value;
  }
  return 0;
}

int gauss_newton_spring::Gra(const double *x, double *gra) const {
  RETURN_WITH_COND_TRUE(max(w_) == 0.0);
  itr_matrix<const double *> X(3, dim_/3, x);
  itr_matrix<double *> G(3, dim_/3, gra);
  for (size_t i = 0; i < edge_.size(2); ++i) {
    matd_t vert = X(colon(), edge_(colon(), i));
    matd_t g = zeros<double>(3, 2);
    mass_spring_jac_(&g[0], &vert[0], &len_[i]);
    G(colon(), edge_(colon(), i)) += w_[i]*g;
  }
  return 0;
}

int gauss_newton_spring::Hes(const double *x, vector<Triplet<double>> *hes) const {
  RETURN_WITH_COND_TRUE(max(w_) == 0.0);
  itr_matrix<const double *> X(3, dim_/3, x);
  for (size_t i = 0; i < edge_.size(2); ++i) {
    matd_t vert = X(colon(), edge_(colon(), i));
    double curr_len = norm(vert(colon(), 0)-vert(colon(), 1));
    matd_t g = zeros<double>(6, 1);
    g(colon(0, 2)) = (vert(colon(), 0)-vert(colon(), 1))/(sqrt(len_[i])*curr_len);
    g(colon(3, 5)) = -g(colon(0, 2));
    matd_t H = g*trans(g);
    for (size_t p = 0; p < 6; ++p) {
      for (size_t q = 0; q < 6; ++q) {
        const size_t I = 3*edge_(p/3, i)+p%3;
        const size_t J = 3*edge_(q/3, i)+q%3;
        if ( H(p, q) != 0.0 )
          hes->push_back(Triplet<double>(I, J, w_[i]*H(p, q)));
      }
    }
  }
  return 0;
}
//==============================================================================
fast_mass_spring::fast_mass_spring(const mati_t &edge, const matd_t &nods, const double w)
  : dim_(nods.size()), w_(w), edge_(edge) {
  len_.resize(edge_.size(2), 1);
  d_.resize(3, edge_.size(2));
#pragma omp parallel for
  for (size_t i = 0; i < edge_.size(2); ++i) {
    d_(colon(), i) = nods(colon(), edge_(0, i))-nods(colon(), edge_(1, i));
    len_[i] = norm(d_(colon(), i));
  }
  S_.resize(aux_dim(), Nx()); {
    vector<Triplet<double>> trips;
    for (size_t i = 0; i < edge_.size(2); ++i) {
      trips.push_back(Triplet<double>(3*i+0, 3*edge_(0, i)+0, 1.0));
      trips.push_back(Triplet<double>(3*i+0, 3*edge_(1, i)+0, -1.0));
      trips.push_back(Triplet<double>(3*i+1, 3*edge_(0, i)+1, 1.0));
      trips.push_back(Triplet<double>(3*i+1, 3*edge_(1, i)+1, -1.0));
      trips.push_back(Triplet<double>(3*i+2, 3*edge_(0, i)+2, 1.0));
      trips.push_back(Triplet<double>(3*i+2, 3*edge_(1, i)+2, -1.0));
    }
    S_.reserve(trips.size());
    S_.setFromTriplets(trips.begin(), trips.end());
  }
}

size_t fast_mass_spring::Nx() const {
  return dim_;
}

int fast_mass_spring::Val(const double *x, double *val) const {
  RETURN_WITH_COND_TRUE(w_ == 0.0);
  itr_matrix<const double *> X(3, dim_/3, x);
  matd_t dx(3, 1);
  for (size_t i = 0; i < edge_.size(2); ++i) {
    dx = X(colon(), edge_(0, i))-X(colon(), edge_(1, i))-d_(colon(), i);
    *val += 0.5*w_*dot(dx, dx);
  }
  return 0;
}

int fast_mass_spring::Gra(const double *x, double *gra) const {
  RETURN_WITH_COND_TRUE(w_ == 0.0);
  itr_matrix<const double *> X(3, dim_/3, x);
  itr_matrix<double *> G(3, dim_/3, gra);
  matd_t dx(3, 1);
  for (size_t i = 0; i < edge_.size(2); ++i) {
    dx = X(colon(), edge_(0, i))-X(colon(), edge_(1, i))-d_(colon(), i);
    G(colon(), edge_(0, i)) += w_*dx;
    G(colon(), edge_(1, i)) -= w_*dx;
  }
  return 0;
}

int fast_mass_spring::Hes(const double *x, vector<Triplet<double>> *hes) const {
  RETURN_WITH_COND_TRUE(w_ == 0.0);
  for (size_t i = 0; i < edge_.size(2); ++i) {
    add_diag_block<double, 3>(edge_(0, i), edge_(0, i), w_, hes);
    add_diag_block<double, 3>(edge_(0, i), edge_(1, i), -w_, hes);
    add_diag_block<double, 3>(edge_(1, i), edge_(0, i), -w_, hes);
    add_diag_block<double, 3>(edge_(1, i), edge_(1, i), w_, hes);
  }
  return 0;
}

void fast_mass_spring::LocalSolve(const double *x) {
  itr_matrix<const double *> X(3, dim_/3, x);
#pragma omp parallel for
  for (size_t i = 0; i < edge_.size(2); ++i) {
    d_(colon(), i) = X(colon(), edge_(0, i))-X(colon(), edge_(1, i));
    double dnorm = norm(d_(colon(), i));
    d_(colon(), i) *= len_[i]/dnorm;
  }
}

void fast_mass_spring::Project() {
#pragma omp parallel for
  for (size_t i = 0; i < d_.size(2); ++i) {
    double dnorm = norm(d_(colon(), i));
    d_(colon(), i) *= len_[i]/dnorm;
  }
}

void fast_mass_spring::build_jts_pattern() {
  vector<Triplet<double>> trips;
  for (size_t i = 0; i < edge_.size(2); ++i) {
    for (size_t p = 0; p < 6; ++p) {
      for (size_t q = 0; q < 6; ++q) {
        const size_t I = 3*edge_(p/3, i)+p%3;
        const size_t J = 3*edge_(q/3, i)+q%3;
        trips.push_back(Triplet<double>(I, J, 0.0));
      }
    }
  }
  JtS_.resize(dim_, dim_);
  JtS_.setFromTriplets(trips.begin(), trips.end());
  JtS_.makeCompressed();
  // J^TS is column major
  for (size_t i = 0; i < edge_.size(2); ++i) {
    for (size_t p = 0; p < 6; ++p) {
      for (size_t q = 0; q < 6; ++q) {
        const size_t I = 3*edge_(p/3, i)+p%3;
        const size_t J = 3*edge_(q/3, i)+q%3;
        for (size_t cnt = JtS_.outerIndexPtr()[J]; cnt < JtS_.outerIndexPtr()[J+1]; ++cnt) {
          if ( JtS_.innerIndexPtr()[cnt] == I )
            ijp_.insert(make_pair(make_pair(I, J), cnt));
        }
      }
    }
  }
}

SparseMatrix<double>& fast_mass_spring::get_jts(const double *x) {
  itr_matrix<const double *> X(3, dim_/3, x);
  std::fill(JtS_.valuePtr(), JtS_.valuePtr()+JtS_.nonZeros(), 0.0);
  for (size_t i = 0; i < edge_.size(2); ++i) {
    matd_t vert = X(colon(), edge_(colon(), i));
    matd_t J = zeros<double>(3, 6);
    const_len_spring_jac(&J[0], &vert[0], &len_[i]);
    matd_t jts(6, 6);
    jts(colon(), colon(0, 2)) = trans(J);
    jts(colon(), colon(3, 5)) = -trans(J);
    for (size_t p = 0; p < 6; ++p) {
      for (size_t q = 0; q < 6; ++q) {
        const size_t I = 3*edge_(p/3, i)+p%3;
        const size_t J = 3*edge_(q/3, i)+q%3;
        JtS_.valuePtr()[ijp_[make_pair(I, J)]] += w_*jts(p, q);
      }
    }
  }
  return JtS_;
}
//==============================================================================
surf_bending_potential::surf_bending_potential(const mati_t &diams, const matd_t &nods, const double w)
  : dim_(nods.size()), diams_(diams), w_(w) {
  len_.resize(diams_.size(2), 1);
  angle_.resize(diams_.size(2), 1);
  area_.resize(diams_.size(2), 1);
#pragma omp parallel for
  for (size_t i = 0; i < diams_.size(2); ++i) {
    matd_t vert = nods(colon(), diams_(colon(), i));
    calc_edge_length_(&len_[i], &vert(0, 1));
    calc_dih_angle_(&angle_[i], &vert[0]);
    area_[i] = calc_tri_area(vert(colon(), colon(0, 2)))+calc_tri_area(vert(colon(), colon(1, 3)));
  }
}

size_t surf_bending_potential::Nx() const {
  return dim_;
}

int surf_bending_potential::Val(const double *x, double *val) const {
  RETURN_WITH_COND_TRUE(w_ == 0.0);
  itr_matrix<const double *> X(3, Nx()/3, x);
  for (size_t i = 0; i < diams_.size(2); ++i) {
    matd_t vert = X(colon(), diams_(colon(), i));
    double value = 0;
    surf_bending_(&value, &vert[0], &angle_[i], &len_[i], &area_[i]);
    *val += w_*value;
  }
  return 0;
}

int surf_bending_potential::Gra(const double *x, double *gra) const {
  RETURN_WITH_COND_TRUE(w_ == 0.0);
  itr_matrix<const double *> X(3, Nx()/3, x);
  itr_matrix<double *> G(3, Nx()/3, gra);
  for (size_t i = 0; i < diams_.size(2); ++i) {
    matd_t vert = X(colon(), diams_(colon(), i));
    matd_t grad = zeros<double>(3, 4);
    surf_bending_jac_(&grad[0], &vert[0], &angle_[i], &len_[i], &area_[i]);
    G(colon(), diams_(colon(), i)) += w_*grad;
  }
  return 0;
}

int surf_bending_potential::Hes(const double *x, vector<Triplet<double>> *hes) const {
  RETURN_WITH_COND_TRUE(w_ == 0.0);
  itr_matrix<const double *> X(3, Nx()/3, x);
  for (size_t i = 0; i < diams_.size(2); ++i) {
    matd_t vert = X(colon(), diams_(colon(), i));
    Matrix<double, 12, 12> H = Matrix<double, 12, 12>::Zero();
//    {
//      double value = 0;
//      surf_bending_(&value, &vert[0], &angle_[i], &len_[i], &area_[i]);
//      Matrix<double, 12, 1> g = Matrix<double, 12, 1>::Zero();
//      surf_bending_jac_(g.data(), &vert[0], &angle_[i], &len_[i], &area_[i]);
//      g = 0.5*g/sqrt(value);
//      H = 2*g*g.transpose();
//    }
    surf_bending_hes_(H.data(), &vert[0], &angle_[i], &len_[i], &area_[i]);
    for (size_t p = 0; p < 12; ++p) {
      for (size_t q = 0; q < 12; ++q) {
        if ( H(p, q) != 0.0 ) {
          const size_t I = 3*diams_(p/3, i)+p%3;
          const size_t J = 3*diams_(q/3, i)+q%3;
          hes->push_back(Triplet<double>(I, J, w_*H(p, q)));
        }
      }
    }
  }
  return 0;
}
//==============================================================================
isometric_bending::isometric_bending(const mati_t &diams, const matd_t &nods, const double w)
  : diams_(diams), w_(w), dim_(nods.size()) {
  area_.resize(diams_.size(2), 1);
  cotv_.resize(4, diams_.size(2));
#pragma omp parallel for
  for (size_t i = 0; i < diams_.size(2); ++i) {
    matd_t vert = nods(colon(), diams_(colon(), i));
    area_[i] = calc_tri_area(vert(colon(), colon(0, 2)))+calc_tri_area(vert(colon(), colon(1, 3)));
    double v012 = cal_cot_val(&vert(0, 0), &vert(0, 1), &vert(0, 2));
    double v021 = cal_cot_val(&vert(0, 0), &vert(0, 2), &vert(0, 1));
    double v312 = cal_cot_val(&vert(0, 3), &vert(0, 1), &vert(0, 2));
    double v321 = cal_cot_val(&vert(0, 3), &vert(0, 2), &vert(0, 1));
    cotv_(0, i) = -v012-v021;
    cotv_(1, i) = v021+v321;
    cotv_(2, i) = v012+v312;
    cotv_(3, i) = -v312-v321;
  }
}

size_t isometric_bending::Nx() const {
  return dim_;
}

int isometric_bending::Val(const double *x, double *val) const {
  RETURN_WITH_COND_TRUE(w_ == 0.0);
  itr_matrix<const double *> X(3, dim_/3, x);
  for (size_t i = 0; i < diams_.size(2); ++i) {
    matd_t vert = X(colon(), diams_(colon(), i));
    matd_t lx = vert*cotv_(colon(), i);
    *val += 0.5*w_*3/area_[i]*dot(lx, lx);
  }
  return 0;
}

int isometric_bending::Gra(const double *x, double *gra) const {
  RETURN_WITH_COND_TRUE(w_ == 0.0);
  itr_matrix<const double *> X(3, dim_/3, x);
  itr_matrix<double *> G(3, dim_/3, gra);
  for (size_t i = 0; i < diams_.size(2); ++i) {
    matd_t vert = X(colon(), diams_(colon(), i));
    matd_t lx = vert*cotv_(colon(), i);
    G(colon(), diams_(colon(), i)) += w_*3/area_[i]*(lx*trans(cotv_(colon(), i)));
  }
  return 0;
}

int isometric_bending::Hes(const double *x, vector<Triplet<double>> *hes) const {
  RETURN_WITH_COND_TRUE(w_ == 0.0);
  for (size_t i = 0; i < diams_.size(2); ++i) {
    matd_t H = w_*3/area_[i]*cotv_(colon(), i)*temp(trans(cotv_(colon(), i)));
    for (size_t p = 0; p < 4; ++p) {
      for (size_t q = 0; q < 4; ++q) {
        add_diag_block<double, 3>(diams_(p, i), diams_(q, i), H(p, q), hes);
      }
    }
  }
  return 0;
}
//==============================================================================
ext_force_energy::ext_force_energy(const matd_t &nods, const double w)
    : dim_(nods.size()), rd_(nods.size(1)), w_(w) {
  force_ = zeros<double>(dim_, 1);
}

size_t ext_force_energy::Nx() const {
  return dim_;
}

int ext_force_energy::Val(const double *x, double *val) const {
  itr_matrix<const double *> X(Nx(), 1, x);
  *val += -w_*dot(force_, X);
  return 0;
}

int ext_force_energy::Gra(const double *x, double *gra) const {
  itr_matrix<const double *> X(Nx(), 1, x);
  itr_matrix<double *> G(Nx(), 1, gra);
  G += w_*-force_;
  return 0;
}

int ext_force_energy::ApplyForce(const size_t id, const double *f) {
  if ( id < 0 || id >= Nx()/rd_ )
    return __LINE__;
  std::copy(f, f+rd_, &force_[rd_*id]);
  return 0;
}

int ext_force_energy::RemoveForce(const size_t id) {
  if ( id < 0 || id >= Nx()/rd_ )
    return __LINE__;
  std::fill(&force_[rd_*id], &force_[rd_*id+rd_], 0);
  return 0;
}
//==============================================================================
line_bending_potential::line_bending_potential(const mati_t &edge, const matd_t &nods, const double w)
  : dim_(nods.size()), w_(w), edge_(edge) {
  len_.resize(edge_.size(2), 1);
#pragma omp parallel for
  for (size_t i = 0; i < edge_.size(2); ++i) {
    matd_t vert = nods(colon(), edge_(colon(), i));
    calc_edge_length_(&len_[i], &vert[0]);
  }
}

size_t line_bending_potential::Nx() const {
  return dim_;
}

int line_bending_potential::Val(const double *x, double *val) const {
  RETURN_WITH_COND_TRUE(w_ == 0.0);
  itr_matrix<const double *> X(3, dim_/3, x);
  for (size_t i = 0; i < edge_.size(2)-1; ++i) {
    matd_t vert(3, 4);
    vert(colon(), colon(0, 1)) = X(colon(), edge_(colon(), i));
    vert(colon(), colon(2, 3)) = X(colon(), edge_(colon(), i+1));
    double value = 0;
    line_bending_(&value, &vert[0], &len_[i], &len_[i+1]);
    *val += w_*value;
  }
  return 0;
}

int line_bending_potential::Gra(const double *x, double *gra) const {
  RETURN_WITH_COND_TRUE(w_ == 0.0);
  itr_matrix<const double *> X(3, dim_/3, x);
  itr_matrix<double *> G(3, dim_/3, gra);
  for (size_t i = 0; i < edge_.size(2)-1; ++i) {
    matd_t vert(3, 4);
    vert(colon(), colon(0, 1)) = X(colon(), edge_(colon(), i));
    vert(colon(), colon(2, 3)) = X(colon(), edge_(colon(), i+1));
    matd_t g(3, 4);
    line_bending_jac_(&g[0], &vert[0], &len_[i], &len_[i+1]);
    G(colon(), edge_(colon(), i)) += w_*g(colon(), colon(0, 1));
    G(colon(), edge_(colon(), i+1)) += w_*g(colon(), colon(2, 3));
  }
  return 0;
}

int line_bending_potential::Hes(const double *x, vector<Triplet<double>> *hes) const {
  RETURN_WITH_COND_TRUE(w_ == 0.0);
  itr_matrix<const double *> X(3, dim_/3, x);
  for (size_t i = 0; i < edge_.size(2)-1; ++i) {
    matd_t vert(3, 4);
    vert(colon(), colon(0, 1)) = X(colon(), edge_(colon(), i));
    vert(colon(), colon(2, 3)) = X(colon(), edge_(colon(), i+1));
    matd_t H(12, 12);
    line_bending_hes_(&H[0], &vert[0], &len_[i], &len_[i+1]);
    for (size_t p = 0; p < 12; ++p) {
      for (size_t q = 0; q < 12; ++q) {
        const size_t I = 3*edge_((p%6)/3, i+p/6)+p%3;
        const size_t J = 3*edge_((q%6)/3, i+q/6)+q%3;
        if ( H(p, q) != 0.0 )
          hes->push_back(Triplet<double>(I, J, w_*H(p, q)));
      }
    }
  }
  return 0;
}
//==============================================================================
tet_arap_energy::tet_arap_energy(const mati_t &tets, const matd_t &nods, const double w)
  : dim_(nods.size()), tets_(tets), w_(w) {
  vol_.resize(tets_.size(2), 1);
  R_.resize(9, tets_.size(2));
  D_.resize(9, tets_.size(2));
#pragma omp parallel for
  for (size_t i = 0; i < tets_.size(2); ++i) {
    matd_t basis = nods(colon(), tets_(colon(1, 3), i))-nods(colon(), tets_(0, i))*ones<double>(1, 3);
    Eigen::Map<Matrix3d> B(basis.begin());
    Eigen::Map<Matrix3d>(&D_(0, i)) = B.inverse();
    vol_[i] = std::fabs(B.determinant())/6.0;
    Eigen::Map<Matrix3d>(&R_(0, i)) = Matrix3d::Identity();
  }
}

size_t tet_arap_energy::Nx() const {
  return dim_;
}

int tet_arap_energy::Val(const double *x, double *val) const {
  RETURN_WITH_COND_TRUE(w_ == 0.0);
  itr_matrix<const double *> X(3, dim_/3, x);
  for (size_t i = 0; i < tets_.size(2); ++i) {
    matd_t vert = X(colon(), tets_(colon(), i));
    double value = 0;
    tet_arap_(&value, &vert[0], &D_(0, i), &R_(0, i), &vol_[i]);
    *val += w_*value;
  }
  return 0;
}

int tet_arap_energy::Gra(const double *x, double *gra) const {
  RETURN_WITH_COND_TRUE(w_ == 0.0);
  itr_matrix<const double *> X(3, dim_/3, x);
  itr_matrix<double *> G(3, dim_/3, gra);
  for (size_t i = 0; i < tets_.size(2); ++i) {
    matd_t vert = X(colon(), tets_(colon(), i));
    matd_t jac = zeros<double>(3, 4);
    tet_arap_jac_(&jac[0], &vert[0], &D_(0, i), &R_(0, i), &vol_[i]);
    G(colon(), tets_(colon(), i)) += w_*jac;
  }
  return 0;
}

int tet_arap_energy::Hes(const double *x, vector<Triplet<double>> *hes) const {
  RETURN_WITH_COND_TRUE(w_ == 0.0);
  for (size_t i = 0; i < tets_.size(2); ++i) {
    matd_t H = zeros<double>(12, 12);
    tet_arap_hes_(&H[0], nullptr, &D_(0, i), &R_(0, i), &vol_[i]);
    for (size_t p = 0; p < 12; ++p) {
      for (size_t q = 0; q < 12; ++q) {
        const size_t I = 3*tets_(p/3, i)+p%3;
        const size_t J = 3*tets_(q/3, i)+q%3;
        hes->push_back(Triplet<double>(I, J, w_*H(p, q)));
      }
    }
  }
  return 0;
}

void tet_arap_energy::LocalSolve(const double *x) {
  itr_matrix<const double *> X(3, dim_/3, x);
#pragma omp parallel for
  for (size_t i = 0; i < tets_.size(2); ++i) {
    matd_t Ds = X(colon(), tets_(colon(1, 3), i))-X(colon(), tets_(0, i))*ones<double>(1, 3);
    matd_t df = Ds*itr_matrix<const double *>(3, 3, &D_(0, i));
    Eigen::Map<Matrix3d> F(df.begin());
    JacobiSVD<Matrix3d> svd(F, ComputeFullU|ComputeFullV);
    Eigen::Map<Matrix3d>(&R_(0, i)) = svd.matrixU()*svd.matrixV().transpose();
  }
}

void tet_arap_energy::CalcLieAlgebraCoord(VectorXd &vec) const {
  vec.setZero(3*tets_.size(2));
#pragma omp parallel for
  for (size_t i = 0; i < R_.size(2); ++i) {
    Matrix3d logR = Matrix3d(&R_(0, i)).log();
    vec[3*i+0] = -logR(1, 2);
    vec[3*i+1] = logR(0, 2);
    vec[3*i+2] = -logR(0, 1);
  }
}

void tet_arap_energy::UpdateRotation(const VectorXd &vec) {
#pragma omp parallel for
  for (size_t i = 0; i < R_.size(2); ++i) {
    Matrix3d logR = Matrix3d::Zero();
    logR(0, 1) = -vec[3*i+2];
    logR(0, 2) = vec[3*i+1];
    logR(1, 2) = -vec[3*i+0];
    logR(1, 0) = -logR(0, 1);
    logR(2, 0) = -logR(0, 2);
    logR(2, 1) = -logR(1, 2);
    Matrix3d rot = logR.exp();
    std::copy(rot.data(), rot.data()+rot.size(), &R_(0, i));
  }
}
//==============================================================================
bw98_stretch_energy::bw98_stretch_energy(const mati_t &tris, const matd_t &nods, const double w)
  : dim_(nods.size()), w_(w), tris_(tris) {
  // build local frame
  matd_t O(3, tris_.size(2)), T(3, tris_.size(2)), B(3, tris_.size(2)), N(3, tris_.size(2));
  jtf::mesh::cal_face_normal(tris, nods, N, true);
#pragma omp parallel for
  for (size_t i = 0; i < tris_.size(2); ++i) {
    const matd_t vert = nods(colon(), tris_(colon(), i));
    O(colon(), i) = vert*ones<double>(3, 1)/3.0;
    T(colon(), i) = vert(colon(), 1)-vert(colon(), 0);
    T(colon(), i) /= norm(T(colon(), i));
    B(colon(), i) = cross(N(colon(), i), T(colon(), i));
    B(colon(), i) /= norm(B(colon(), i));
  }
  // calc inv(u, v) and area
  invUV_.resize(4, tris_.size(2));
  area_.resize(tris_.size(2), 1);
#pragma omp parallel for
  for (size_t i = 0; i < tris_.size(2); ++i) {
    const matd_t vert = nods(colon(), tris_(colon(), i));
    matd_t uv(2, 3);
    for (size_t j = 0; j < 3; ++j) {
      uv(0, j) = dot(vert(colon(), j)-O(colon(), i), T(colon(), i));
      uv(1, j) = dot(vert(colon(), j)-O(colon(), i), B(colon(), i));
    }
    matd_t base = uv(colon(), colon(1, 2))-uv(colon(), 0)*ones<double>(1, 2);
    if ( inv(base) ) {
      cerr << "\t@degenerated triangle " << i << endl;
      exit(EXIT_FAILURE);
    }
    std::copy(base.begin(), base.end(), &invUV_(0, i));
    area_[i] = calc_tri_area(vert);
  }
}

size_t bw98_stretch_energy::Nx() const {
  return dim_;
}

int bw98_stretch_energy::Val(const double *x, double *val) const {
  RETURN_WITH_COND_TRUE(w_ == 0.0);
  itr_matrix<const double *> X(3, dim_/3, x);
  for (size_t i = 0; i < tris_.size(2); ++i) {
    matd_t vert = X(colon(), tris_(colon(), i));
    double value = 0;
    bw98_stretch_(&value, &vert[0], &invUV_(0, i), &area_[i]);
    *val += w_*value;
  }
  return 0;
}

int bw98_stretch_energy::Gra(const double *x, double *gra) const {
  RETURN_WITH_COND_TRUE(w_ == 0.0);
  itr_matrix<const double *> X(3, dim_/3, x);
  itr_matrix<double *> G(3, dim_/3, gra);
  for (size_t i = 0; i < tris_.size(2); ++i) {
    matd_t vert = X(colon(), tris_(colon(), i));
    matd_t jac = zeros<double>(3, 3);
    bw98_stretch_jac_(&jac[0], &vert[0], &invUV_(0, i), &area_[i]);
    G(colon(), tris_(colon(), i)) += w_*jac;
  }
  return 0;
}

int bw98_stretch_energy::Hes(const double *x, vector<Triplet<double>> *hes) const {
  RETURN_WITH_COND_TRUE(w_ == 0.0);
  itr_matrix<const double *> X(3, dim_/3, x);
  for (size_t i = 0; i < tris_.size(2); ++i) {
    matd_t vert = X(colon(), tris_(colon(), i));
    matd_t H = zeros<double>(9, 9);
    bw98_stretch_hes_(&H[0], &vert[0], &invUV_(0, i), &area_[i]);
    for (size_t p = 0; p < 9; ++p) {
      for (size_t q = 0; q < 9; ++q) {
        const size_t I = 3*tris_(p/3, i)+p%3;
        const size_t J = 3*tris_(q/3, i)+q%3;
        hes->push_back(Triplet<double>(I, J, w_*H(p, q)));
      }
    }
  }
  return 0;
}
//==============================================================================
bw98_shear_energy::bw98_shear_energy(const mati_t &tris, const matd_t &nods, const double w)
  : dim_(nods.size()), w_(w), tris_(tris) {
  // build local frame
  matd_t O(3, tris_.size(2)), T(3, tris_.size(2)), B(3, tris_.size(2)), N(3, tris_.size(2));
  jtf::mesh::cal_face_normal(tris, nods, N, true);
#pragma omp parallel for
  for (size_t i = 0; i < tris_.size(2); ++i) {
    const matd_t vert = nods(colon(), tris_(colon(), i));
    O(colon(), i) = vert*ones<double>(3, 1)/3.0;
    T(colon(), i) = vert(colon(), 1)-vert(colon(), 0);
    T(colon(), i) /= norm(T(colon(), i));
    B(colon(), i) = cross(N(colon(), i), T(colon(), i));
    B(colon(), i) /= norm(B(colon(), i));
  }
  // calc inv(u, v) and area
  invUV_.resize(4, tris_.size(2));
  area_.resize(tris_.size(2), 1);
#pragma omp parallel for
  for (size_t i = 0; i < tris_.size(2); ++i) {
    const matd_t vert = nods(colon(), tris_(colon(), i));
    matd_t uv(2, 3);
    for (size_t j = 0; j < 3; ++j) {
      uv(0, j) = dot(vert(colon(), j)-O(colon(), i), T(colon(), i));
      uv(1, j) = dot(vert(colon(), j)-O(colon(), i), B(colon(), i));
    }
    matd_t base = uv(colon(), colon(1, 2))-uv(colon(), 0)*ones<double>(1, 2);
    inv(base);
    std::copy(base.begin(), base.end(), &invUV_(0, i));
    area_[i] = calc_tri_area(vert);
  }
}

size_t bw98_shear_energy::Nx() const {
  return dim_;
}

int bw98_shear_energy::Val(const double *x, double *val) const {
  RETURN_WITH_COND_TRUE(w_ == 0.0);
  itr_matrix<const double *> X(3, dim_/3, x);
  for (size_t i = 0; i < tris_.size(2); ++i) {
    matd_t vert = X(colon(), tris_(colon(), i));
    double value = 0;
    bw98_shear_(&value, &vert[0], &invUV_(0, i), &area_[i]);
    *val += w_*value;
  }
  return 0;
}

int bw98_shear_energy::Gra(const double *x, double *gra) const {
  RETURN_WITH_COND_TRUE(w_ == 0.0);
  itr_matrix<const double *> X(3, dim_/3, x);
  itr_matrix<double *> G(3, dim_/3, gra);
  for (size_t i = 0; i < tris_.size(2); ++i) {
    matd_t vert = X(colon(), tris_(colon(), i));
    matd_t jac = zeros<double>(3, 3);
    bw98_shear_jac_(&jac[0], &vert[0], &invUV_(0, i), &area_[i]);
    G(colon(), tris_(colon(), i)) += w_*jac;
  }
  return 0;
}

int bw98_shear_energy::Hes(const double *x, vector<Triplet<double>> *hes) const {
  RETURN_WITH_COND_TRUE(w_ == 0.0);
  itr_matrix<const double *> X(3, dim_/3, x);
  for (size_t i = 0; i < tris_.size(2); ++i) {
    matd_t vert = X(colon(), tris_(colon(), i));
    matd_t H = zeros<double>(9, 9);
    bw98_shear_hes_(&H[0], &vert[0], &invUV_(0, i), &area_[i]);
    for (size_t p = 0; p < 9; ++p) {
      for (size_t q = 0; q < 9; ++q) {
        const size_t I = 3*tris_(p/3, i)+p%3;
        const size_t J = 3*tris_(q/3, i)+q%3;
        hes->push_back(Triplet<double>(I, J, w_*H(p, q)));
      }
    }
  }
  return 0;
}
//==============================================================================
fem_stretch_energy::fem_stretch_energy(const mati_t &tris, const matd_t &nods, const double w)
  : dim_(nods.size()), w_(w), tris_(tris) {

}

size_t fem_stretch_energy::Nx() const {
  return dim_;
}

int fem_stretch_energy::Val(const double *x, double *val) const {
  RETURN_WITH_COND_TRUE(w_ == 0.0);
  itr_matrix<const double *> X(3, dim_/3, x);
  for (size_t i = 0; i < tris_.size(2); ++i) {
    matd_t vert = X(colon(), tris_(colon(), i));
    double value = 0;
    fem_stretch_(&value, &vert[0], &Dm_(0, i), &area_[i], &K_(0, i));
    *val += w_*value;
  }
  return 0;
}

int fem_stretch_energy::Gra(const double *x, double *gra) const {
  RETURN_WITH_COND_TRUE(w_ == 0.0);
  itr_matrix<const double *> X(3, dim_/3, x);
  itr_matrix<double *> G(3, dim_/3, gra);
  for (size_t i = 0; i < tris_.size(2); ++i) {
    matd_t vert = X(colon(), tris_(colon(), i));
    matd_t jac = zeros<double>(3, 3);
    fem_stretch_jac_(&jac[0], &vert[0], &Dm_(0, i), &area_[i], &K_(0, i));
    G(colon(), tris_(colon(), i)) += w_*jac;
  }
  return 0;
}

int fem_stretch_energy::Hes(const double *x, vector<Triplet<double>> *hes) const {
  RETURN_WITH_COND_TRUE(w_ == 0.0);
  itr_matrix<const double *> X(3, dim_/3, x);
  for (size_t i = 0; i < tris_.size(2); ++i) {
    // TODO
  }
  return 0;
}
//==============================================================================
low_pass_filter_energy::low_pass_filter_energy(const mati_t &tris, const matd_t &nods, const size_t patch_num, const double w)
  : dim_(nods.size()), w_(w), tris_(tris), patch_num_(patch_num) {
  mesh_partition mp(tris, nods);
  vector<ptn_to_patch> info;
  mp.init(info);
  mp.run(patch_num_, info);
  // allocate patches
  pat_.resize(mp.get_actual_patch_num());
  for (size_t i = 0; i < info.size(); ++i) {
    pat_[info[i].id_patch].push_back(make_pair(i, info[i].dist));
  }
  // calc coef
  double sigma = 10;//2*()*pow(1, 1.0/3);
#pragma omp parallel for
  for (size_t i = 0; i < pat_.size(); ++i) {
    double sum = 0;
    for (size_t j = 0; j < pat_[i].size(); ++j) {
      double d = pat_[i][j].second;
      if ( d < 2*sigma )
        pat_[i][j].second = exp(-0.5*d*d/(sigma*sigma));
      else
        pat_[i][j].second = 0;
      sum += pat_[i][j].second;
    }
    for (size_t j = 0; j < pat_[i].size(); ++j) {
      pat_[i][j].second /= sum;
    }
  }
}

size_t low_pass_filter_energy::Nx() const {
  return dim_;
}

int low_pass_filter_energy::Val(const double *x, double *val) const {
  RETURN_WITH_COND_TRUE(w_ == 0.0);
  Eigen::Map<const MatrixXd> X(x, 3, dim_/3);
  Eigen::Map<const MatrixXd> Q(ref_, 3, dim_/3);
  for (auto &arr : pat_) {
    Vector3d lhs = Vector3d::Zero(), rhs = Vector3d::Zero();
    for (auto &elem: arr) {
      const size_t pi = elem.first;
      lhs += elem.second*X.col(pi);
      rhs += elem.second*Q.col(pi);
    }
    *val += 0.5*w_*(lhs-rhs).squaredNorm();
  }
  return 0;
}

int low_pass_filter_energy::Gra(const double *x, double *gra) const {
  RETURN_WITH_COND_TRUE(w_ == 0.0);
  Eigen::Map<const MatrixXd> X(x, 3, dim_/3);
  Eigen::Map<const MatrixXd> Q(ref_, 3, dim_/3);
  Eigen::Map<MatrixXd> G(gra, 3, dim_/3);
  for (auto &arr : pat_) {
    Vector3d lhs = Vector3d::Zero(), rhs = Vector3d::Zero();
    for (auto &elem : arr) {
      const size_t pi = elem.first;
      lhs += elem.second*X.col(pi);
      rhs += elem.second*Q.col(pi);
    }
    for (auto &elem : arr) {
      G.col(elem.first) += w_*elem.second*(lhs-rhs);
    }
  }
  return 0;
}

int low_pass_filter_energy::Hes(const double *x, vector<Triplet<double>> *hes) const {
  RETURN_WITH_COND_TRUE(w_ == 0.0);
  for (auto &arr : pat_) {
    for (size_t i = 0; i < arr.size(); ++i) {
      for (size_t j = 0; j < arr.size(); ++j) {
        add_diag_block<double, 3>(arr[i].first, arr[j].first, w_*arr[i].second*arr[j].second, hes);
      }
    }
  }
  return 0;
}

void low_pass_filter_energy::Update(const double *ref) {
  ref_ = ref;
}
//==============================================================================
cosserat_stretch_energy::cosserat_stretch_energy(const mati_t &rod, const matd_t &nods,
                                                 const double Es, const double r, const double w)
  : rod_(rod), r_size_(nods.size()), q_size_(4*(rod.size()-1)), elem_num_(rod.size()-1), Es_(Es), r_(r), w_(w) {
  len_ = VectorXd::Zero(elem_num_);
  for (size_t i = 0; i < elem_num_; ++i) {
    len_(i) = norm(nods(colon(), rod_[i])-nods(colon(), rod_[i+1]));
  }
}

size_t cosserat_stretch_energy::Nx() const {
  return r_size_+q_size_;
}

int cosserat_stretch_energy::Val(const double *xq, double *val) const {
  Eigen::Map<const VectorXd> X(xq, r_size_+q_size_);
  for (size_t i = 0; i < elem_num_; ++i) {
    Matrix<double, 3, 2> rr;
    rr.col(0) = X.segment<3>(3*rod_[i]);
    rr.col(1) = X.segment<3>(3*rod_[i+1]);
    double value = 0;
    rod_stretch_(&value, rr.data(), &len_(i), &Es_, &r_);
    *val += w_*value;
  }
  return 0;
}

int cosserat_stretch_energy::Gra(const double *xq, double *gra) const {
  Eigen::Map<const VectorXd> X(xq, r_size_+q_size_);
  Eigen::Map<VectorXd> G(gra, r_size_+q_size_);
  for (size_t i = 0; i < elem_num_; ++i) {
    Matrix<double, 3, 2> rr;
    rr.col(0) = X.segment<3>(3*rod_[i]);
    rr.col(1) = X.segment<3>(3*rod_[i+1]);
    Matrix<double, 3, 2> g = Matrix<double, 3, 2>::Zero();
    rod_stretch_jac_(g.data(), rr.data(), &len_(i), &Es_, &r_);
    G.segment<3>(3*rod_[i]) += w_*g.col(0);
    G.segment<3>(3*rod_[i+1]) += w_*g.col(1);
  }
  return 0;
}

int cosserat_stretch_energy::Hes(const double *xq, vector<Triplet<double>> *hes) const {
  Eigen::Map<const VectorXd> X(xq, r_size_+q_size_);
  for (size_t i = 0; i < elem_num_; ++i) {
    Matrix<double, 3, 2> rr;
    rr.col(0) = X.segment<3>(3*rod_[i]);
    rr.col(1) = X.segment<3>(3*rod_[i+1]);
    Matrix<double, 6, 6> H = Matrix<double, 6, 6>::Zero();
    if ( (rr.col(0)-rr.col(1)).norm() >= len_(i) ) {
      rod_stretch_hes_(H.data(), rr.data(), &len_(i), &Es_, &r_);
    } else {
      double value = 0;
      rod_stretch_(&value, rr.data(), &len_(i), &Es_, &r_);
      Matrix<double, 6, 1> g = Matrix<double, 6, 1>::Zero();
      rod_stretch_jac_(g.data(), rr.data(), &len_(i), &Es_, &r_);
      g = 0.5*g/sqrt(value);
      H = 2*g*g.transpose();
    }
    for (size_t p = 0; p < 6; ++p) {
      for (size_t q = 0; q < 6; ++q) {
        const size_t I = 3*rod_[i+p/3]+p%3;
        const size_t J = 3*rod_[i+q/3]+q%3;
        hes->push_back(Triplet<double>(I, J, w_*H(p, q)));
      }
    }
  }
  return 0;
}
//==============================================================================
cosserat_bend_energy::cosserat_bend_energy(const mati_t &rod, const matd_t &nods,
                                           const double E, const double G, const double r, const double w)
  : rod_(rod), r_size_(nods.size()), q_size_(4*(rod.size()-1)), elem_num_(rod.size()-2), E_(E), G_(G), r_(r), w_(w) {
  len_ = VectorXd::Zero(elem_num_);
  for (size_t i = 0; i < elem_num_; ++i) {
    len_(i) = 0.5*(norm(nods(colon(), rod_[i])-nods(colon(), rod_[i+1]))
        +norm(nods(colon(), rod_[i+1])-nods(colon(), rod_[i+2])));
  }
}

size_t cosserat_bend_energy::Nx() const {
  return r_size_+q_size_;
}

int cosserat_bend_energy::Val(const double *xq, double *val) const {
  Eigen::Map<const VectorXd> X(xq, r_size_+q_size_);
  for (size_t i = 0; i < elem_num_; ++i) {
    Matrix<double, 4, 2> qq;
    qq.col(0) = X.segment<4>(r_size_+4*i);
    qq.col(1) = X.segment<4>(r_size_+4*(i+1));
    double value = 0;
    const double uk[3] = {0, 0, 0};
    rod_bend_(&value, qq.data(), &uk[0], &len_(i), &E_, &G_, &r_);
    *val += w_*value;
  }
  return 0;
}

int cosserat_bend_energy::Gra(const double *xq, double *gra) const {
  Eigen::Map<const VectorXd> X(xq, r_size_+q_size_);
  Eigen::Map<VectorXd> G(gra, r_size_+q_size_);
  for (size_t i = 0; i < elem_num_; ++i) {
    Matrix<double, 4, 2> qq;
    qq.col(0) = X.segment<4>(r_size_+4*i);
    qq.col(1) = X.segment<4>(r_size_+4*(i+1));
    Matrix<double, 4, 2> g = Matrix<double, 4, 2>::Zero();
    const double uk[3] = {0, 0, 0};
    rod_bend_jac_(g.data(), qq.data(), &uk[0], &len_(i), &E_, &G_, &r_);
    G.segment<4>(r_size_+4*i) += w_*g.col(0);
    G.segment<4>(r_size_+4*(i+1)) += w_*g.col(1);
  }
  return 0;
}

int cosserat_bend_energy::Hes(const double *xq, vector<Triplet<double>> *hes) const {
  Eigen::Map<const VectorXd> X(xq, r_size_+q_size_);
  for (size_t i = 0; i < elem_num_; ++i) {
    Matrix<double, 4, 2> qq;
    qq.col(0) = X.segment<4>(r_size_+4*i);
    qq.col(1) = X.segment<4>(r_size_+4*(i+1));
    Matrix<double, 8, 8> H = Matrix<double, 8, 8>::Zero();
    const double uk[3] = {0, 0, 0};
    rod_bend_hes_(H.data(), qq.data(), &uk[0], &len_(i), &E_, &G_, &r_);
    for (size_t p = 0; p < 8; ++p) {
      for (size_t q = 0; q < 8; ++q) {
        const size_t I = r_size_+4*(i+p/4)+p%4;
        const size_t J = r_size_+4*(i+q/4)+q%4;
        hes->push_back(Triplet<double>(I, J, w_*H(p, q)));
      }
    }
  }
  return 0;
}
//==============================================================================
cosserat_couple_energy::cosserat_couple_energy(const mati_t &rod, const matd_t &nods,
                                               const double kappa, const double w)
  : rod_(rod), r_size_(nods.size()), q_size_(4*(rod.size()-1)), elem_num_(rod.size()-1), kappa_(kappa), w_(w) {
  len_ = VectorXd::Zero(elem_num_);
  for (size_t i = 0; i < elem_num_; ++i) {
    len_(i) = norm(nods(colon(), rod[i])-nods(colon(), rod[i+1]));
  }
}

size_t cosserat_couple_energy::Nx() const {
  return r_size_+q_size_;
}

int cosserat_couple_energy::Val(const double *xq, double *val) const {
  Eigen::Map<const VectorXd> XQ(xq, r_size_+q_size_);
  for (size_t i = 0; i < elem_num_; ++i) {
    VectorXd rq = VectorXd::Zero(10);
    rq.segment<3>(0) = XQ.segment<3>(3*rod_[i]);
    rq.segment<3>(3) = XQ.segment<3>(3*rod_[i+1]);
    rq.segment<4>(6) = XQ.segment<4>(r_size_+4*i);
    double value = 0;
    rod_couple_(&value, rq.data(), &len_(i), &kappa_);
    *val += w_*value;
  }
  return 0;
}

int cosserat_couple_energy::Gra(const double *xq, double *gra) const {
  Eigen::Map<const VectorXd> XQ(xq, r_size_+q_size_);
  Eigen::Map<VectorXd> G(gra, r_size_+q_size_);
  for (size_t i = 0; i < elem_num_; ++i) {
    VectorXd rq = VectorXd::Zero(10);
    rq.segment<3>(0) = XQ.segment<3>(3*rod_[i]);
    rq.segment<3>(3) = XQ.segment<3>(3*rod_[i+1]);
    rq.segment<4>(6) = XQ.segment<4>(r_size_+4*i);
    VectorXd g = VectorXd::Zero(10);
    rod_couple_jac_(g.data(), rq.data(), &len_(i), &kappa_);
    G.segment<3>(3*rod_[i]) += w_*g.segment<3>(0);
    G.segment<3>(3*rod_[i+1]) += w_*g.segment<3>(3);
    G.segment<4>(r_size_+4*i) += w_*g.segment<4>(6);
  }
  return 0;
}

int cosserat_couple_energy::Hes(const double *xq, vector<Triplet<double>> *hes) const {
  Eigen::Map<const VectorXd> XQ(xq, r_size_+q_size_);
  for (size_t i = 0; i < elem_num_; ++i) {
    VectorXd rq = VectorXd::Zero(10);
    rq.segment<3>(0) = XQ.segment<3>(3*rod_[i]);
    rq.segment<3>(3) = XQ.segment<3>(3*rod_[i+1]);
    rq.segment<4>(6) = XQ.segment<4>(r_size_+4*i);
    Matrix<double, 10, 10> H = Matrix<double, 10, 10>::Zero();
    rod_couple_hes_(H.data(), rq.data(), &len_(i), &kappa_);
    for (size_t p = 0; p < 10; ++p) {
      for (size_t q = 0; q < 10; ++q) {
        const size_t I = p < 6 ? 3*rod_[i+p/3]+p%3 : r_size_+4*i+p-6;
        const size_t J = q < 6 ? 3*rod_[i+q/3]+q%3 : r_size_+4*i+q-6;
        hes->push_back(Triplet<double>(I, J, w_*H(p, q)));
      }
    }
  }
  return 0;
}
//==============================================================================
}
