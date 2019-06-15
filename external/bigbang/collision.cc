#include "collision.h"

#include <iostream>
#include <zjucad/matrix/itr_matrix.h>
#include "config.h"

using namespace std;
using namespace Eigen;
using namespace zjucad::matrix;

namespace bigbang {

//===============================================================================
lineSDF::lineSDF(const double *center, const double *n)
    : C_(Vector2d(center)), N_(Vector2d(n).normalized()) {}

bool lineSDF::inside(const double *x) const {
  Vector2d X(x);
  return (X-C_).dot(N_) < 0.0;
}

void lineSDF::Val(const double *x, double *val) const {
  Vector2d X(x);
  *val = pow((X-C_).dot(N_), 2);
}

void lineSDF::Gra(const double *x, double *gra) const {
  Vector2d X(x);
  Eigen::Map<Vector2d> G(gra);
  G = 2*(X-C_).dot(N_)*N_;
}

void lineSDF::Hes(const double *x, double *hes) const {
  Eigen::Map<Matrix2d> H(hes);
  H = 2*N_*N_.transpose();
}
//===============================================================================
circleSDF::circleSDF(const double *center, const double r)
    : C_(Vector2d(center)), R_(r) {}

bool circleSDF::inside(const double *x) const {
  Vector2d X(x);
  return (X-C_).norm()-R_ < 0;
}

void circleSDF::Val(const double *x, double *val) const {
  Vector2d X(x);
  *val = pow((X-C_).norm()-R_, 2);
}

void circleSDF::Gra(const double *x, double *gra) const {
  Vector2d X(x);
  Eigen::Map<Vector2d> G(gra);
  G = 2*(X-C_-R_*(X-C_)/(X-C_).norm());
}

void circleSDF::Hes(const double *x, double *hes) const {
  Vector2d X(x);
  Eigen::Map<Matrix2d> H(hes);
  Vector2d J = (X-C_)/(X-C_).norm();
  H = 2*J*J.transpose();
}
//===============================================================================
planeSDF::planeSDF(const double *center, const double *n)
    : C_(Vector3d(center)), N_(Vector3d(n).normalized()) {}

bool planeSDF::inside(const double *x) const {
  Vector3d X(x);
  return (X-C_).dot(N_) < 0.0;
}

void planeSDF::Val(const double *x, double *val) const {
  Vector3d X(x);
  *val = pow((X-C_).dot(N_), 2);
}

void planeSDF::Gra(const double *x, double *gra) const {
  Vector3d X(x);
  Eigen::Map<Vector3d> G(gra);
  G = 2*(X-C_).dot(N_)*N_;
}

void planeSDF::Hes(const double *x, double *hes) const {
  Eigen::Map<Matrix3d> H(hes);
  H = 2*N_*N_.transpose();
}
//===============================================================================
sphereSDF::sphereSDF(const double *center, const double r)
    : C_(Vector3d(center)), R_(r) {}

bool sphereSDF::inside(const double *x) const {
  Vector3d X(x);
  return (X-C_).norm()-R_ < 0;
}

void sphereSDF::Val(const double *x, double *val) const {
  Vector3d X(x);
  *val = pow((X-C_).norm()-R_, 2);
}

void sphereSDF::Gra(const double *x, double *gra) const {
  Vector3d X(x);
  Eigen::Map<Vector3d> G(gra);
  G = 2*(X-C_-R_*(X-C_)/(X-C_).norm());
}

void sphereSDF::Hes(const double *x, double *hes) const {
  Vector3d X(x);
  Eigen::Map<Matrix3d> H(hes);
  Vector3d J = (X-C_)/(X-C_).norm();
  H = 2*J*J.transpose();
}
//===============================================================================
extern "C" {
  void torus_sdf_(double *val, const double *X, const double *C, const double *N, const double *r, const double *R);
  void torus_sdf_jac_(double *jac, const double *X, const double *C, const double *N, const double *r, const double *R);  
}

torusSDF::torusSDF(const double *center, const double *n, const double r, const double R)
    : C_(Vector3d(center)), N_(Vector3d(n).normalized()), r_(r), R_(R) {
  ASSERT(0 < r_ && r_ < R_);
}

bool torusSDF::inside(const double *x) const {
  double d = 0;
  torus_sdf_(&d, x, C_.data(), N_.data(), &r_, &R_);
  return d < 0;
}

void torusSDF::Val(const double *x, double *val) const {
  double d = 0;
  torus_sdf_(&d, x, C_.data(), N_.data(), &r_, &R_);
  *val = d*d;
}

void torusSDF::Gra(const double *x, double *gra) const {
  Eigen::Map<Vector3d> G(gra);
  double d = 0;
  torus_sdf_(&d, x, C_.data(), N_.data(), &r_, &R_);
  Vector3d g;
  torus_sdf_jac_(g.data(), x, C_.data(), N_.data(), &r_, &R_);
  G = 2*d*g;
}

void torusSDF::Hes(const double *x, double *hes) const {
  Eigen::Map<Matrix3d> H(hes);
  Vector3d g;
  torus_sdf_jac_(g.data(), x, C_.data(), N_.data(), &r_, &R_);
  H = 2*g*g.transpose();
}
//===============================================================================
extern "C" {
  void cylinder_sdf_(double *val, const double *X, const double *C, const double *N, const double *r);
  void cylinder_sdf_jac_(double *jac, const double *X, const double *C, const double *N, const double *r);
}

cylinderSDF::cylinderSDF(const double *center, const double *n, const double r)
    : C_(Vector3d(center)), N_(Vector3d(n).normalized()), R_(r) {
  ASSERT(R_ > 0);
}

bool cylinderSDF::inside(const double *x) const {
  double d = 0;
  cylinder_sdf_(&d, x, C_.data(), N_.data(), &R_);
  return d < 0;
}

void cylinderSDF::Val(const double *x, double *val) const {
  double d = 0;
  cylinder_sdf_(&d, x, C_.data(), N_.data(), &R_);
  *val = d*d;
}

void cylinderSDF::Gra(const double *x, double *gra) const {
  Eigen::Map<Vector3d> G(gra);
  double d = 0;
  cylinder_sdf_(&d, x, C_.data(), N_.data(), &R_);
  Vector3d g;
  cylinder_sdf_jac_(g.data(), x, C_.data(), N_.data(), &R_);
  G = 2*d*g;  
}

void cylinderSDF::Hes(const double *x, double *hes) const {
  Eigen::Map<Matrix3d> H(hes);
  Vector3d g;
  cylinder_sdf_jac_(g.data(), x, C_.data(), N_.data(), &R_);
  H = 2*g*g.transpose();
}
//===============================================================================
geom_contact_energy::geom_contact_energy(const vector<shared_ptr<signed_dist_func>> &objs,
                                         const size_t dim, const double w, const size_t rd)
    : objs_(objs), dim_(dim), w_(w), rd_(rd) {}

size_t geom_contact_energy::Nx() const {
  return dim_;
}

int geom_contact_energy::Val(const double *x, double *val) const {
  RETURN_WITH_COND_TRUE(w_ == 0.0);
  for (size_t i = 0; i < dim_/rd_; ++i) {
    for (size_t j = 0; j < objs_.size(); ++j) {
      if ( !objs_[j].get() )
        continue;
      if ( objs_[j]->inside(x+rd_*i) ) {
        double value = 0;
        objs_[j]->Val(x+rd_*i, &value);
        *val += w_*value;
      }
    }
  }
  return 0;
}

int geom_contact_energy::Gra(const double *x, double *gra) const {
  RETURN_WITH_COND_TRUE(w_ == 0.0);
  Eigen::Map<MatrixXd> G(gra, rd_, dim_/rd_);
  for (size_t i = 0; i < dim_/rd_; ++i) {
    for (size_t j = 0; j < objs_.size(); ++j) {
      if ( !objs_[j].get() )
        continue;
      if ( objs_[j]->inside(x+rd_*i) ) {
        VectorXd g(rd_);
        objs_[j]->Gra(x+rd_*i, g.data());
        G.col(i) += w_*g;
      }
    }
  }
  return 0;
}

int geom_contact_energy::Hes(const double *x, vector<Triplet<double>> *hes) const {
  RETURN_WITH_COND_TRUE(w_ == 0.0);
  for (size_t i = 0; i < dim_/rd_; ++i) {
    for (size_t j = 0; j < objs_.size(); ++j) {
      if ( !objs_[j].get() )
        continue;
      if ( objs_[j]->inside(x+rd_*i) ) {
        MatrixXd H(rd_, rd_);
        objs_[j]->Hes(x+rd_*i, H.data());
        for (size_t p = 0; p < rd_; ++p)
          for (size_t q = 0; q < rd_; ++q)
            hes->push_back(Triplet<double>(rd_*i+p, rd_*i+q, w_*H(p, q)));
      }
    }
  }
  return 0;
}
//===============================================================================
}
