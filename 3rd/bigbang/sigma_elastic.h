#ifndef SIGMA_ELASTIC_H
#define SIGMA_ELASTIC_H

#include <zjucad/matrix/matrix.h>
#include <Eigen/Sparse>
#include <zjucad/matrix/itr_matrix.h>
#include <jtflib/mesh/mesh.h>
#include <Eigen/CholmodSupport>

#include "def.h"

using mati_t=zjucad::matrix::matrix<size_t>;
using matd_t=zjucad::matrix::matrix<double>;

namespace bigbang {

// Psi(F)=\mu*\|E\|^2+\lambda/2*tr(E)^2
// Psi(\sigma) = \mu*\sum_i(\sigma_i^2-1)^2/4+\lambda/2(\sum_i\sigma_i^2-3)^2/4

#define SQ(x) ((x)*(x))

class cons_law_func
{
public:
  virtual std::string name() const { return "empty"; }
  virtual size_t num_args() const = 0;
  virtual void Val(const double s1, const double s2, const double s3, const double *args, double *val) const = 0;
  virtual void Jac(const double s1, const double s2, const double s3, const double *args, double *jac) const = 0;
  virtual void Hes(const double s1, const double s2, const double s3, const double *args, double *hes) const = 0;
  virtual void dVdp(const double s1, const double s2, const double s3, const double *args, double *vp) const = 0;
  virtual void dJdp(const double s1, const double s2, const double s3, const double *args, double *jp) const = 0;
};

class coro_linear_psi : public cons_law_func
{
public:
  std::string name() const {
    return "coro";
  }
  size_t num_args() const {
    return 2;
  }
  void Val(const double s1, const double s2, const double s3, const double *args, double *val) const {
    const double mu = args[0], lambda = args[1];
    *val += mu*(SQ(s1-1)+SQ(s2-1)+SQ(s3-1))+0.5*lambda*SQ(s1+s2+s3-3);
  }
  void Jac(const double s1, const double s2, const double s3, const double *args, double *jac) const {
    const double mu = args[0], lambda = args[1];
    jac[0] += 2*mu*(s1-1)+lambda*(s1+s2+s3-3);
    jac[4] += 2*mu*(s2-1)+lambda*(s1+s2+s3-3);
    jac[8] += 2*mu*(s3-1)+lambda*(s1+s2+s3-3);
  }
  void Hes(const double s1, const double s2, const double s3, const double *args, double *hes) const {
    zjucad::matrix::itr_matrix<double *> H(3, 3, hes);
    const double mu = args[0], lambda = args[1];
    H(0, 0) += 2*mu+lambda; H(0, 1) += lambda;       H(0, 2) += lambda;
    H(1, 0) += lambda;      H(1, 1) += 2*mu+lambda;  H(1, 2) += lambda;
    H(2, 0) += lambda;      H(2, 1) += lambda;       H(2, 2) += 2*mu+lambda;
  }
  void dVdp(const double s1, const double s2, const double s3, const double *args, double *vp) const {
    vp[0] += SQ(s1-1)+SQ(s2-1)+SQ(s3-1);
    vp[1] += 0.5*SQ(s1+s2+s3-3);
  }
  void dJdp(const double s1, const double s2, const double s3, const double *args, double *jp) const {
    jp[0] += 2*(s1-1); jp[3] += (s1+s2+s3-3);
    jp[1] += 2*(s2-1); jp[4] += (s1+s2+s3-3);
    jp[2] += 2*(s3-1); jp[5] += (s1+s2+s3-3);    
  }
};

class stvk_psi : public cons_law_func
{
public:
  std::string name() const {
    return "stvk";
  }
  size_t num_args() const {
    return 2;
  }
  void Val(const double s1, const double s2, const double s3, const double *args, double *val) const {
    const double mu = args[0], lambda = args[1];
    *val += 0.25*mu*( SQ(SQ(s1)-1)+SQ(SQ(s2)-1)+SQ(SQ(s3)-1) )+0.125*lambda*( SQ( SQ(s1)+SQ(s2)+SQ(s3)-3 ) );
  }
  void Jac(const double s1, const double s2, const double s3, const double *args, double *jac) const {
    const double mu = args[0], lambda = args[1];
    jac[0] += mu*(SQ(s1)-1)*s1+0.5*lambda*(SQ(s1)+SQ(s2)+SQ(s3)-3)*s1;
    jac[4] += mu*(SQ(s2)-1)*s2+0.5*lambda*(SQ(s1)+SQ(s2)+SQ(s3)-3)*s2;
    jac[8] += mu*(SQ(s3)-1)*s3+0.5*lambda*(SQ(s1)+SQ(s2)+SQ(s3)-3)*s3;
  }
  void Hes(const double s1, const double s2, const double s3, const double *args, double *hes) const {
    zjucad::matrix::itr_matrix<double *> H(3, 3, hes);
    const double mu = args[0], lambda = args[1];
    H(0, 0) += mu*(3*SQ(s1)-1)+lambda*SQ(s1)+0.5*lambda*(SQ(s1)+SQ(s2)+SQ(s3)-3);
    H(0, 1) += lambda*s1*s2;
    H(0, 2) += lambda*s1*s3;
    H(1, 0) += lambda*s1*s2;
    H(1, 1) += mu*(3*SQ(s2)-1)+lambda*SQ(s2)+0.5*lambda*(SQ(s1)+SQ(s2)+SQ(s3)-3);
    H(1, 2) += lambda*s2*s3;
    H(2, 0) += lambda*s1*s3;
    H(2, 1) += lambda*s2*s3;
    H(2, 2) += mu*(3*SQ(s3)-1)+lambda*SQ(s3)+0.5*lambda*(SQ(s1)+SQ(s2)+SQ(s3)-3);
  }
  void dVdp(const double s1, const double s2, const double s3, const double *args, double *vp) const {
    vp[0] += 0.25*( SQ(SQ(s1)-1)+SQ(SQ(s2)-1)+SQ(SQ(s3)-1) );
    vp[1] += 0.125*( SQ( SQ(s1)+SQ(s2)+SQ(s3)-3 ) );
  }
  void dJdp(const double s1, const double s2, const double s3, const double *args, double *jp) const {
    jp[0] += (SQ(s1)-1)*s1; jp[3] += 0.5*(SQ(s1)+SQ(s2)+SQ(s3)-3)*s1;
    jp[1] += (SQ(s2)-1)*s2; jp[4] += 0.5*(SQ(s1)+SQ(s2)+SQ(s3)-3)*s2;
    jp[2] += (SQ(s3)-1)*s3; jp[5] += 0.5*(SQ(s1)+SQ(s2)+SQ(s3)-3)*s3;
  }
};

class neohookean_psi : public cons_law_func
{
public:
  std::string name() const {
    return "neohookean";
  }
  size_t num_args() const {
    return 2;
  }
  void Val(const double s1, const double s2, const double s3, const double *args, double *val) const {
    const double mu = args[0], lambda = args[1];
    *val += 0.5*mu*(SQ(s1)+SQ(s2)+SQ(s3)-2*log(s1*s2*s3)-3)+0.5*lambda*SQ(log(s1*s2*s3));
  }
  void Jac(const double s1, const double s2, const double s3, const double *args, double *jac) const {
    const double mu = args[0], lambda = args[1];
    jac[0] += mu*(s1-1.0/s1)+lambda*log(s1*s2*s3)/s1;
    jac[4] += mu*(s2-1.0/s2)+lambda*log(s1*s2*s3)/s2;
    jac[8] += mu*(s3-1.0/s3)+lambda*log(s1*s2*s3)/s3;
  }
  void Hes(const double s1, const double s2, const double s3, const double *args, double *hes) const {
    zjucad::matrix::itr_matrix<double *> H(3, 3, hes);
    const double mu = args[0], lambda = args[1];
    H(0, 0) += lambda/SQ(s1)*(1-log(s1*s2*s3))+mu*(1+1/SQ(s1));
    H(0, 1) += lambda/(s1*s2);
    H(0, 2) += lambda/(s1*s3);
    H(1, 0) += lambda/(s1*s2);
    H(1, 1) += lambda/SQ(s2)*(1-log(s1*s2*s3))+mu*(1+1/SQ(s2));
    H(1, 2) += lambda/(s2*s3);
    H(2, 0) += lambda/(s1*s3);
    H(2, 1) += lambda/(s2*s3);
    H(2, 2) += lambda/SQ(s3)*(1-log(s1*s2*s3))+mu*(1+1/SQ(s3));
  }
  void dVdp(const double s1, const double s2, const double s3, const double *args, double *vp) const {
    vp[0] += 0.5*(SQ(s1)+SQ(s2)+SQ(s3)-2*log(s1*s2*s3)-3);
    vp[1] += 0.5*SQ(log(s1*s2*s3));
  }
  void dJdp(const double s1, const double s2, const double s3, const double *args, double *jp) const {
    jp[0] += (s1-1.0/s1); jp[3] += log(s1*s2*s3)/s1;
    jp[1] += (s2-1.0/s2); jp[4] += log(s1*s2*s3)/s2;
    jp[2] += (s3-1.0/s3); jp[5] += log(s1*s2*s3)/s3;
  }
};

class blend_psi : public cons_law_func
{
public:
  std::string name() const {
    return "blend";
  }
  blend_psi(const std::shared_ptr<cons_law_func> &psiA,
            const std::shared_ptr<cons_law_func> &psiB)
      : psiA_(psiA), psiB_(psiB) {}
  size_t num_args() const {
    return psiA_->num_args()+psiB_->num_args();
  }
  void Val(const double s1, const double s2, const double s3, const double *args, double *val) const {
    double v1 = 0, v2 = 0;
    psiA_->Val(s1, s2, s3, args, &v1);
    psiB_->Val(s1, s2, s3, args+psiA_->num_args(), &v2);
    
    const double dist = SQ(s1-1)+SQ(s2-1)+SQ(s3-1);
    const double w = 1/(4096*pow(dist, 3)+1);

    *val += w*v1+(1-w)*v2;
  }
  void Jac(const double s1, const double s2, const double s3, const double *args, double *jac) const {
    double v1 = 0, v2 = 0;
    psiA_->Val(s1, s2, s3, args, &v1);
    psiB_->Val(s1, s2, s3, args+psiA_->num_args(), &v2);
    double g1[9] = {0}, g2[9] = {0};
    psiA_->Jac(s1, s2, s3, args, g1);
    psiB_->Jac(s1, s2, s3, args+psiA_->num_args(), g2);
    Eigen::Vector3d dv1 = Eigen::Vector3d(g1[0], g1[4], g1[8]);
    Eigen::Vector3d dv2 = Eigen::Vector3d(g2[0], g2[4], g2[8]);
    
    const double dist = SQ(s1-1)+SQ(s2-1)+SQ(s3-1);
    const double w = 1/(4096*pow(dist, 3)+1);
    Eigen::Vector3d dw = -24576*SQ(dist*w)*Eigen::Vector3d(s1-1, s2-1, s3-1);

    Eigen::Vector3d jacobian = dw*v1+w*dv1 + (-dw)*v2+(1-w)*dv2;
    jac[0] += jacobian[0];
    jac[4] += jacobian[1];
    jac[8] += jacobian[2];
  }
  void Hes(const double s1, const double s2, const double s3, const double *args, double *hes) const {}
  void dVdp(const double s1, const double s2, const double s3, const double *args, double *vp) const {
    Eigen::Map<Eigen::VectorXd> VP(vp, this->num_args());
    psiA_->dVdp(s1, s2, s3, args, VP.data());
    psiB_->dVdp(s1, s2, s3, args+psiA_->num_args(), VP.data()+psiA_->num_args());

    const double dist = SQ(s1-1)+SQ(s2-1)+SQ(s3-1);
    const double w = 1/(4096*pow(dist, 3)+1);

    VP.head(psiA_->num_args()) *= w;
    VP.tail(psiB_->num_args()) *= 1-w;    
  }
  void dJdp(const double s1, const double s2, const double s3, const double *args, double *jp) const {
    Eigen::Map<Eigen::MatrixXd> JP(jp, 3, this->num_args());

    const double dist = SQ(s1-1)+SQ(s2-1)+SQ(s3-1);
    const double w = 1/(4096*pow(dist, 3)+1);
    Eigen::Vector3d dw = -24576*SQ(dist*w)*Eigen::Vector3d(s1-1, s2-1, s3-1);

    Eigen::VectorXd Vpa = Eigen::VectorXd::Zero(psiA_->num_args());
    Eigen::MatrixXd Jpa = Eigen::MatrixXd::Zero(3, psiA_->num_args());
    psiA_->dVdp(s1, s2, s3, args, Vpa.data());
    psiA_->dJdp(s1, s2, s3, args, Jpa.data());
    JP.topLeftCorner(3, psiA_->num_args()) += dw*Vpa.transpose()+w*Jpa;

    Eigen::VectorXd Vpb = Eigen::VectorXd::Zero(psiB_->num_args());
    Eigen::MatrixXd Jpb = Eigen::MatrixXd::Zero(3, psiB_->num_args());
    psiB_->dVdp(s1, s2, s3, args+psiA_->num_args(), Vpb.data());
    psiB_->dJdp(s1, s2, s3, args+psiA_->num_args(), Jpb.data());
    JP.topRightCorner(3, psiB_->num_args()) += -dw*Vpb.transpose()+(1-w)*Jpb;
  }
private:
  std::shared_ptr<cons_law_func> psiA_, psiB_;
};

class hybrid_psi : public cons_law_func
{
public:
  std::string name() const {
    return "hybrid";
  }
  hybrid_psi() {
    law_vec_.resize(3);
    law_vec_[0] = std::make_shared<coro_linear_psi>();
    law_vec_[1] = std::make_shared<stvk_psi>();
    law_vec_[2] = std::make_shared<neohookean_psi>();
  }
  size_t num_args() const {
    size_t cnt = 0;
    for (auto &ptr : law_vec_)
      cnt += ptr->num_args();
    return cnt;
  }
  void Val(const double s1, const double s2, const double s3, const double *args, double *val) const {
    size_t curr = 0;
    for (auto &ptr : law_vec_) {
      ptr->Val(s1, s2, s3, &args[curr], val);
      curr += ptr->num_args();
    }
  }
  void Jac(const double s1, const double s2, const double s3, const double *args, double *jac) const {
    size_t curr = 0;
    for (auto &ptr : law_vec_) {
      ptr->Jac(s1, s2, s3, &args[curr], jac);
      curr += ptr->num_args();
    }
  }
  void Hes(const double s1, const double s2, const double s3, const double *args, double *hes) const {
    size_t curr = 0;
    for (auto &ptr : law_vec_) {
      ptr->Hes(s1, s2, s3, &args[curr], hes);
      curr += ptr->num_args();
    }
  }
  void dVdp(const double s1, const double s2, const double s3, const double *args, double *vp) const {
    size_t curr = 0;
    for (auto &ptr : law_vec_) {
      ptr->dVdp(s1, s2, s3, &args[curr], &vp[curr]);
      curr += ptr->num_args();
    }
  }
  void dJdp(const double s1, const double s2, const double s3, const double *args, double *jp) const {
    using zjucad::matrix::colon;
    size_t curr = 0;
    for (auto &ptr : law_vec_) {
      ptr->dJdp(s1, s2, s3, &args[curr], &jp[3*curr]);
      curr += ptr->num_args();
    }
  }
private:
  std::vector<std::shared_ptr<cons_law_func>> law_vec_;
};

#undef SQ

//-> used for regressing coarse material
class vox_force_matching_energy : public Functional<double>
{
public:
  vox_force_matching_energy(const mati_t &cube, const matd_t &nods,
                            const matd_t &Rx, const matd_t &Rf,
                            const matd_t &datw,
                            const std::shared_ptr<cons_law_func> &psi,
                            const std::vector<std::vector<size_t>> &fixv,
                            const std::string &metric);
  size_t Nx() const;
  int Val(const double *x, double *val) const;
  int Gra(const double *x, double *gra) const;
  int Hes(const double *x, std::vector<Eigen::Triplet<double>> *hes) const;
protected:
  const mati_t &cube_;
  const matd_t &nods_;
  const double w_;
  const matd_t &Rx_, &Rf_;  // restricted positions and forces
  const matd_t &datw_;      // weight for each data sample
  const std::vector<std::vector<size_t>> &fixv_;
  
  // computing forces
  matd_t quadrature_, qw_;
  matd_t detXH_;
  matd_t H_invXH_, GradOp_;

  // pattern selector
  const std::shared_ptr<cons_law_func> psi_;
  size_t args_num_;

  const double force_eps_;
  const std::string metric_;

  std::shared_ptr<Functional<double>> elas_;
  std::vector<std::shared_ptr<Eigen::CholmodSimplicialLDLT<Eigen::SparseMatrix<double>>>> solvers_;
  std::vector<Eigen::SparseMatrix<double>> K_;
};

//-> used for material regularization
class vox_args_smooth_energy : public Functional<double>
{
public:
  vox_args_smooth_energy(const mati_t &cube, const matd_t &nods,
                         const std::shared_ptr<cons_law_func> &psi, const double w);
  size_t Nx() const;
  int Val(const double *x, double *val) const;
  int Gra(const double *x, double *gra) const;
  int Hes(const double *x, std::vector<Eigen::Triplet<double>> *hes) const;
private:
  const mati_t &cube_;
  const double w_;
  double volume_;
  
  std::shared_ptr<jtf::mesh::face2hex_adjacent> f2v_;
  const std::shared_ptr<cons_law_func> psi_;
  size_t args_num_;
};

//-> used for resimulate
class vox_sigma_elastic_energy : public Functional<double>
{
public:
  vox_sigma_elastic_energy(const mati_t &cube, const matd_t &nods, const matd_t &args,
                           const std::shared_ptr<cons_law_func> &psi, const double w);
  size_t Nx() const;
  int Val(const double *x, double *val) const;
  int Gra(const double *x, double *gra) const;
  int Hes(const double *x, std::vector<Eigen::Triplet<double>> *hes) const;
private:
  const mati_t &cube_;
  const size_t dim_;
  const matd_t &args_;
  const double w_;

  matd_t quadrature_, qw_;
  matd_t detXH_;
  matd_t H_invXH_, GradOp_;

  const std::shared_ptr<cons_law_func> psi_;
};

}

#endif
