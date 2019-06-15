#include "constraint.h"

#include <iostream>
#include <zjucad/matrix/itr_matrix.h>

#include "config.h"

using namespace std;
using namespace Eigen;
using namespace zjucad::matrix;

namespace bigbang {

//==============================================================================
position_constraint::position_constraint(const mati_t &point, const size_t dim, const double *pos)
  : dim_(dim), constraint_piece(point, 1.0, constraint_piece::EQUAL) {
  pos_.resize(3, 1);
  std::copy(pos, pos+3, pos_.begin());
}

int position_constraint::eval_val(const double *x, double *val) const {
  itr_matrix<const double *> X(3, dim_/3, x);
  matd_t dx = X(colon(), pn_)-pos_;
  *val = dot(dx, dx);
  return 0;
}

int position_constraint::eval_jac(const double *x, double *jac) const {
  itr_matrix<const double *> X(3, dim_/3, x);
  matd_t dx = 2*(X(colon(), pn_)-pos_);
  std::copy(dx.begin(), dx.end(), jac);
  return 0;
}
//==============================================================================
inext_constraint::inext_constraint(const mati_t &edge, const matd_t &nods)
  : dim_(nods.size()), constraint_piece(edge, 1.0, constraint_piece::EQUAL) {
  matd_t vert = nods(colon(), pn_);
  len_ = norm(vert(colon(), 0)-vert(colon(), 1));
}

int inext_constraint::eval_val(const double *x, double *val) const {
  matd_t X = itr_matrix<const double *>(3, dim_/3, x)(colon(), pn_);
  matd_t ed = X(colon(), 0)-X(colon(), 1);
  *val = dot(ed, ed)/len_-len_;
  return 0;
}

int inext_constraint::eval_jac(const double *x, double *jac) const {
  matd_t X = itr_matrix<const double *>(3, dim_/3, x)(colon(), pn_);
  matd_t grad = 2*(X(colon(), 0)-X(colon(), 1))/len_;
  std::copy(grad.begin(), grad.end(), jac);
  grad *= -1.0;
  std::copy(grad.begin(), grad.end(), jac+3);
  return 0;
}
//===============================================================================
asm_constraint::asm_constraint(const vector<shared_ptr<constraint_piece<double>>> &buffer)
  : buffer_(buffer) {
  if ( buffer_.empty() ) {
    cerr << "[error] no input constraints\n";
    exit(EXIT_FAILURE);
  }
  dim_ = buffer_[0]->dim();
  for (auto &elem : buffer_) {
    if ( dim_ != elem->dim() ) {
      cerr << "[error] incompatible constraints\n";
      exit(EXIT_FAILURE);
    }
  }
}

int asm_constraint::Val(const double *x, double *val) const {
#pragma omp parallel for
  for (size_t i = 0; i < buffer_.size(); ++i) {
    buffer_[i]->eval_val(x, val+i);
  }
  return 0;
}

int asm_constraint::Jac(const double *x, const size_t off, vector<Triplet<double>> *jac) const {
  for (size_t i = 0; i < buffer_.size(); ++i) {
    matd_t grad(3*buffer_[i]->pn_.size(), 1);
    buffer_[i]->eval_jac(x, &grad[0]);
    for (size_t j = 0; j < grad.size(); ++j) {
      if ( grad[j]*grad[j] != 0.0 )
        jac->push_back(Triplet<double>(off+i, 3*buffer_[i]->pn_[j/3]+j%3, grad[j]));
    }
  }
  return 0;
}
//==============================================================================
zero_moment_constraint::zero_moment_constraint(const matd_t &nods)
    : nods_(nods), dim_(nods.size()), rd_(nods.size(1)) {}

size_t zero_moment_constraint::Nx() const {
  return dim_;
}

size_t zero_moment_constraint::Nf() const {
  return rd_;
}

int zero_moment_constraint::Val(const double *x, double *val) const {
  itr_matrix<const double *> X(rd_, dim_/rd_, x);
  itr_matrix<double *> V(Nf(), 1, val);
  for (size_t i = 0; i < nods_.size(2); ++i)
    V += X(colon(), i)-nods_(colon(), i);
  return 0;
}

int zero_moment_constraint::Jac(const double *x, const size_t off, vector<Triplet<double>> *jac) const {
  for (size_t i = 0; i < nods_.size(2); ++i) {
    for (size_t j = 0; j < rd_; ++j) {
      jac->push_back(Triplet<double>(off+j, rd_*i+j, 1));
    }
  }
  return 0;
}

int zero_moment_constraint::Hes(const double *x, const size_t off, vector<vector<Triplet<double>>> *hes) const {
  return __LINE__;
}
//===============================================================================
first_moment_constraint::first_moment_constraint(const matd_t &nods)
    : nods_(nods), dim_(nods.size()), rd_(nods.size(1)) {
  xbar_ = nods_*ones<double>(nods_.size(2), 1)/nods_.size(2);
}

size_t first_moment_constraint::Nx() const {
  return dim_;
}

size_t first_moment_constraint::Nf() const {
  if ( rd_ == 3 )
    return 3;
  if ( rd_ == 2 )
    return 1;
}

int first_moment_constraint::Val(const double *x, double *val) const {
  itr_matrix<const double *> X(rd_, Nx()/rd_, x);
  itr_matrix<double *> V(Nf(), 1, val);

  if ( rd_ == 3 ) {
    for (size_t i = 0; i < nods_.size(2); ++i) {
      const matd_t XXi = X(colon(), i)-nods_(colon(), i);
      const matd_t XXb = X(colon(), i)-xbar_;
      V += cross(XXi, XXb);
    }
  } else if ( rd_ == 2 ) {
    for (size_t i = 0; i < nods_.size(2); ++i) {
      const matd_t XXi = X(colon(), i)-nods_(colon(), i);
      const matd_t XXb = X(colon(), i)-xbar_;
      *val += XXi[0]*XXb[1]-XXi[1]*XXb[0];
    }
  }
  return 0;
}

int first_moment_constraint::Jac(const double *x, const size_t off, vector<Triplet<double>> *jac) const {
  if ( rd_ == 3 ) {
    for (size_t i = 0; i < nods_.size(2); ++i) {
      const double j01 = nods_(2, i)-xbar_[2];
      jac->push_back(Triplet<double>(off+0, 3*i+1, j01));
      jac->push_back(Triplet<double>(off+1, 3*i+0, -j01));

      const double j02 = xbar_[1]-nods_(1, i);
      jac->push_back(Triplet<double>(off+0, 3*i+2, j02));
      jac->push_back(Triplet<double>(off+2, 3*i+0, -j02));

      const double j12 = nods_(0, i)-xbar_[0];
      jac->push_back(Triplet<double>(off+1, 3*i+2, j12));
      jac->push_back(Triplet<double>(off+2, 3*i+1, -j12));
    }
  } else if ( rd_ == 2 ) {
    for (size_t i = 0; i < nods_.size(2); ++i) {
      jac->push_back(Triplet<double>(off+0, 2*i+0, nods_(1, i)-xbar_[1]));
      jac->push_back(Triplet<double>(off+0, 2*i+1, xbar_[0]-nods_(0, i)));
    }
  }
  return 0;
}

int first_moment_constraint::Hes(const double *x, const size_t off, vector<vector<Triplet<double>>> *hes) const {
  return __LINE__;
}
//===============================================================================
hard_position_constraint::hard_position_constraint(const matd_t &nods)
    : dim_(nods.size()), rd_(nods.size(1)) {}

size_t hard_position_constraint::Nx() const {
  return dim_;
}

size_t hard_position_constraint::Nf() const {
  if ( rd_ == 3 )
    return 3*fixed3d_.size();
  if ( rd_ == 2 )
    return 2*fixed2d_.size();
}

int hard_position_constraint::Val(const double *x, double *val) const {
  Eigen::Map<const MatrixXd> X(x, rd_, Nx()/rd_);
  Eigen::Map<MatrixXd> V(val, rd_, Nf()/rd_);

  if ( rd_ == 3 ) {
    size_t i = 0;
    for (auto &elem : fixed3d_) {
      const size_t pid = elem.first;
      V.col(i) += X.col(pid)-elem.second;
      ++i;
    }
    return 0;
  }

  if ( rd_ ==  2 ) {
    size_t i = 0;
    for (auto &elem : fixed2d_) {
      const size_t pid = elem.first;
      V.col(i) += X.col(pid)-elem.second;
      ++i;
    }
    return 0;
  }
  
  return __LINE__;
}

int hard_position_constraint::Jac(const double *x, const size_t off, vector<Triplet<double>> *jac) const {
  if ( rd_ == 3 ) {
    size_t i = 0;
    for (auto &elem : fixed3d_) {
      const size_t pid = elem.first;
      jac->push_back(Triplet<double>(off+3*i+0, 3*pid+0, 1));
      jac->push_back(Triplet<double>(off+3*i+1, 3*pid+1, 1));
      jac->push_back(Triplet<double>(off+3*i+2, 3*pid+2, 1));
      ++i;
    }
    return 0;
  }
  
  if ( rd_ == 2 ) {    
    size_t i = 0;
    for (auto &elem : fixed2d_) {
      const size_t pid = elem.first;
      jac->push_back(Triplet<double>(off+2*i+0, 2*pid+0, 1));
      jac->push_back(Triplet<double>(off+2*i+1, 2*pid+1, 1));
      ++i;
    }
    return 0;
  }
  
  return __LINE__;
}

int hard_position_constraint::Hes(const double *x, const size_t off, vector<vector<Triplet<double>>> *hes) const {
  return __LINE__;
}

//===============================================================================
axis_restrict_constraint::axis_restrict_constraint(const matd_t &nods, const char axis)
    : dim_(nods.size()), rd_(nods.size(1)), axis_(axis-'X'), rest_(nods) {
  ASSERT(axis_ >= 0 && axis_ < rd_);
}

size_t axis_restrict_constraint::Nx() const {
  return dim_;
}

size_t axis_restrict_constraint::Nf() const {
  return dim_/rd_;
}

int axis_restrict_constraint::Val(const double *x, double *val) const {
  itr_matrix<const double *> X(rd_, dim_/rd_, x);
  for (size_t i = 0; i < X.size(2); ++i) {
    val[i] += X(axis_, i)-rest_(axis_, i);
  }
  return 0;
}

int axis_restrict_constraint::Jac(const double *x, const size_t off, vector<Triplet<double>> *jac) const {
  for (size_t i = 0; i < dim_/rd_; ++i) {
    jac->push_back(Triplet<double>(off+i, rd_*i+axis_, 1.0));
  }
  return 0;
}

int axis_restrict_constraint::Hes(const double *x, const size_t off, vector<vector<Triplet<double>>> *hes) const {
  return __LINE__;
}

}
