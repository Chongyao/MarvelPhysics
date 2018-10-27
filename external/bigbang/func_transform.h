#ifndef FUNC_TRANSFORM_H
#define FUNC_TRANSFORM_H

#include <hjlib/math_func/math_func.h>
#include "def.h"

namespace bigbang {

using hj::math_func::coo2val_t;
using hj::math_func::func_ctx;
using hj::math_func::coo_set;
using hj::math_func::coo_l2g;

template <typename value_type, typename int_type>
class func_transformer : public hj::math_func::math_func_t<value_type, int_type>
{
public:
  func_transformer(const std::shared_ptr<Functional<value_type>> &func, const value_type *x)
      : func_(func), dim_(func->Nx()) {
    func->Hes(x, &trips_);
    H_.resize(dim_, dim_);
    H_.reserve(trips_.size());
    H_.setFromTriplets(trips_.begin(), trips_.end());
    H_.makeCompressed();
  }
  size_t nx() const {
    return dim_;
  }
  size_t nf() const {
    return 1;
  }
  int eval(size_t k, const value_type *x, const coo2val_t<value_type, int_type> &cv,
           func_ctx *ctx = 0) const {
    if ( k == 0 ) {
      value_type value = 0;
      func_->Val(x, &value);
      int_type c[] = {0};
      cv[c] += value;
      return 0;
    }
    if ( k == 1 ) {
      Eigen::Matrix<value_type, -1, 1> g;
      g.setZero(dim_);
      func_->Gra(x, g.data());
      for (int_type i = 0; i < g.size(); ++i) {
        int_type c[] = {0, i};
        cv[c] += g[i];
      }
      return 0;
    }
    if ( k == 2 ) {
      trips_.clear();
      func_->Hes(x, &trips_);
      for (auto &it : trips_) {
        int_type c[] = {0, it.row(), it.col()};
        cv[c] += it.value();
      }
      return 0;
    }
  }
  int patt(size_t k, coo_set<int_type> &cs, const coo_l2g &l2g,
           func_ctx *ctx = 0) const {
    if ( k == 2 ) {
      for (int_type j = 0; j < H_.outerSize(); ++j) {
        for (typename Eigen::SparseMatrix<value_type>::InnerIterator it(H_, j); it; ++it) {
          int_type c[] = {0, (int_type)it.row(), (int_type)it.col()};
          l2g.add(cs, c);
        }
      }
    }
    return 0;
  }
  size_t nnz(size_t k) const {
    if ( k == 0 )
      return -1;
    if ( k == 1 )
      return -1;
    if ( k == 2 )
      return H_.nonZeros();
  }
private:
  const std::shared_ptr<Functional<value_type>> func_;
  const int_type dim_;
  Eigen::SparseMatrix<value_type> H_;
  mutable std::vector<Eigen::Triplet<value_type>> trips_;
};

}
#endif
