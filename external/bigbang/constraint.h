#ifndef CONSTRAINT_H
#define CONSTRAINT_H

#include <unordered_map>

#include "def.h"

namespace bigbang {

using mati_t=zjucad::matrix::matrix<size_t>;
using matd_t=zjucad::matrix::matrix<double>;

/// $(x-p)^2=0$, nonlinear constraint
class position_constraint : public constraint_piece<double>
{
public:
  position_constraint(const mati_t &point, const size_t dim, const double *pos);
  size_t dim() const { return dim_; }
  int eval_val(const double *x, double *val) const;
  int eval_jac(const double *x, double *jac) const;
private:
  const size_t dim_;
  matd_t pos_;
};

/// $(l_{curr}-l_{rest})^2/d-d=0$, nonlinear constraint
class inext_constraint : public constraint_piece<double>
{
public:
  inext_constraint(const mati_t &edge, const matd_t &nods);
  size_t dim() const { return dim_; }
  int eval_val(const double *x, double *val) const;
  int eval_jac(const double *x, double *jac) const;
private:
  const size_t dim_;
  double len_;
};

class asm_constraint : public Constraint<double>
{
public:
  asm_constraint(const std::vector<std::shared_ptr<constraint_piece<double>>> &buffer);
  size_t Nx() const { return dim_; }
  size_t Nf() const { return buffer_.size(); }
  int Val(const double *x, double *val) const;
  int Jac(const double *x, const size_t off, std::vector<Eigen::Triplet<double>> *jac) const;
private:
  size_t dim_;
  const std::vector<std::shared_ptr<constraint_piece<double>>> &buffer_;
};

class zero_moment_constraint : public Constraint<double>
{
public:
  typedef Eigen::Triplet<double> TPL;
  zero_moment_constraint(const matd_t &nods);
  size_t Nx() const;
  size_t Nf() const;
  int Val(const double *x, double *val) const;
  int Jac(const double *x, const size_t off, std::vector<TPL> *jac) const;
  int Hes(const double *x, const size_t off, std::vector<std::vector<TPL>> *hes) const;
private:
  const size_t dim_, rd_;
  const matd_t nods_;
};

class first_moment_constraint : public Constraint<double>
{
public:
  typedef Eigen::Triplet<double> TPL;
  first_moment_constraint(const matd_t &nods);
  size_t Nx() const;
  size_t Nf() const;
  int Val(const double *x, double *val) const;
  int Jac(const double *x, const size_t off, std::vector<TPL> *jac) const;
  int Hes(const double *x, const size_t off, std::vector<std::vector<TPL>> *hes) const;
private:
  const size_t dim_, rd_;
  const matd_t &nods_;
  matd_t xbar_;
};

class hard_position_constraint : public Constraint<double>
{
public:
  typedef Eigen::Triplet<double> TPL;
  hard_position_constraint(const matd_t &nods);
  size_t Nx() const;
  size_t Nf() const;
  int Val(const double *x, double *val) const;
  int Jac(const double *x, const size_t off, std::vector<TPL> *jac) const;
  int Hes(const double *x, const size_t off, std::vector<std::vector<TPL>> *hes) const;
  template <class Container>
  void update_fixed_verts(const Container &pids, const double *pos) {
    if ( rd_ == 3 ) {
      fixed3d_.clear();
      for (auto &id : pids) {
        assert(id >= 0 && id < Nx()/3);
        fixed3d_.insert(std::make_pair(id, Eigen::Vector3d(pos[3*id+0], pos[3*id+1], pos[3*id+2])));
      }
      return;
    }

    if ( rd_ == 2 ) {
      fixed2d_.clear();
      for (auto &id : pids) {
        assert(id >= 0 && id < Nx()/2);
        fixed2d_.insert(std::make_pair(id, Eigen::Vector2d(pos[2*id+0], pos[2*id+1])));
      }
      return;
    }    
  }
private:
  const size_t dim_, rd_;
  std::unordered_Eigen::Map<size_t, Eigen::Vector3d> fixed3d_;
  std::unordered_Eigen::Map<size_t, Eigen::Vector2d> fixed2d_;
};

class axis_restrict_constraint : public Constraint<double>
{
public:
  typedef Eigen::Triplet<double> TPL;
  axis_restrict_constraint(const matd_t &nods, const char axis='X');
  size_t Nx() const;
  size_t Nf() const;
  int Val(const double *x, double *val) const;
  int Jac(const double *x, const size_t off, std::vector<TPL> *jac) const;
  int Hes(const double *x, const size_t off, std::vector<std::vector<TPL>> *hes) const;
  void set_cons_axis(const char axis) { axis_ = axis-'X'; }
private:
  const size_t dim_, rd_;
  size_t axis_;
  const matd_t rest_;
};

class dofs_constraint : public Constraint<double>
{
public:
  typedef Eigen::Triplet<double> TPL;
  dofs_constraint(const matd_t &rest, const std::vector<size_t> &dofs)
      : rest_(rest), dofs_(dofs), dim_(rest.size()) {}
  size_t Nx() const {
    return dim_;
  }
  size_t Nf() const {
    return dofs_.size();
  }
  int Val(const double *x, double *val) const {
    using zjucad::matrix::colon;
    for (size_t i = 0; i < dofs_.size(); ++i) {
      val[i] += x[dofs_[i]]-rest_(colon())[dofs_[i]];
    }
    return 0;
  }
  int Jac(const double *x, const size_t off, std::vector<TPL> *jac) const {
    for (size_t i = 0; i < dofs_.size(); ++i)
      jac->push_back(TPL(off+i, dofs_[i], 1.0));
    return 0;
  }
  int Hes(const double *x, const size_t off, std::vector<std::vector<TPL>> *hes) const {
    return __LINE__;
  }
private:
  const size_t dim_;
  const matd_t &rest_;
  const std::vector<size_t> &dofs_;
};

class slip_constraint : public Constraint<double>
{
public:
  typedef Eigen::Triplet<double> TPL;
  slip_constraint(const matd_t &nods);
  size_t Nx() const;
  size_t Nf() const;
  int Val(const double *x, double *val) const;
  int Jac(const double *x, const size_t off, std::vector<TPL> *jac) const;
  int Hes(const double *x, const size_t off, std::vector<std::vector<TPL>> *hes) const;
  void set_slip_subspace(const matd_t &basis);
private:
  const size_t dim_, rd_;
  matd_t A_;
  const matd_t rest_;
};

}
#endif

