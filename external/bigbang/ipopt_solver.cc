#include "ipopt_solver.h"

#include <iostream>

#include "config.h"
#include "def.h"

using namespace std;
using namespace Eigen;

namespace bigbang {

ipopt_opt_framework::ipopt_opt_framework(const std::shared_ptr<Functional<Number>> &obj,
                                         const std::shared_ptr<Constraint<Number>> &con,
                                         Number *x0)
    : obj_(obj), con_(con), x0_(x0), dim_(obj_->Nx())
{
  //=> obj can not be null, con can (unconstrained problem).
  K_.resize(dim_, dim_); {
    g_trips_.clear();
    obj_->Hes(x0, &g_trips_);
    K_.setFromTriplets(g_trips_.begin(), g_trips_.end());
  }

  SparseMatrix<Number> H_(dim_, dim_);
  H_.setZero();
  
  if ( con_.get() ) {
    J_.resize(con_->Nf(), con_->Nx());
    vector<Triplet<Number>> trips;
    con_->Jac(x0, 0, &trips);
    J_.setFromTriplets(trips.begin(), trips.end());

    vector<vector<Triplet<Number>>> tripsH;
    con_->Hes(x0, 0, &tripsH);
    for (auto &trisI : tripsH) {
      SparseMatrix<Number> h(dim_, dim_);
      h.setFromTriplets(trisI.begin(), trisI.end());
      H_ += h;
    }
  }

  lagH_ = K_+H_;
  nnz_lagH_ = 0;
  for (size_t j = 0; j < lagH_.outerSize(); ++j) {
    for (SparseMatrix<Number>::InnerIterator it(lagH_, j); it; ++it) {
      if ( it.col() <= it.row() )
        ++nnz_lagH_;
    }
  }
}

ipopt_opt_framework::~ipopt_opt_framework() {}

bool ipopt_opt_framework::get_nlp_info(Ipopt::Index& n, Ipopt::Index& m, Ipopt::Index& nnz_jac_g,
                                       Ipopt::Index& nnz_h_lag, IndexStyleEnum& index_style)
{
  n = obj_->Nx();
  m = con_.get() ? con_->Nf() : 0;
  nnz_jac_g = J_.nonZeros();
  nnz_h_lag = nnz_lagH_;
  index_style = TNLP::C_STYLE;

  return true;
}

bool ipopt_opt_framework::get_bounds_info(Ipopt::Index n, Number* x_l, Number* x_u,
                                          Ipopt::Index m, Number* g_l, Number* g_u)
{
  for (Ipopt::Index i = 0; i < n; ++i) {
    x_l[i] = -1e19;
    x_u[i] = 1e19;
  }

  if ( con_.get() ) {  // equality constraint
    for (Ipopt::Index i = 0; i < m; ++i)
      g_l[i] = g_u[i] = 0;
  }

  return true;
}

bool ipopt_opt_framework::get_starting_point(Ipopt::Index n, bool init_x, Number* x,
                                             bool init_z, Number* z_L, Number* z_U,
                                             Ipopt::Index m, bool init_lambda,
                                             Number* lambda)
{
  assert(init_x == true);
  assert(init_z == false);
  assert(init_lambda == false);

  std::copy(x0_, x0_+dim_, x);

  return true;
}

bool ipopt_opt_framework::eval_f(Ipopt::Index n, const Number* x, bool new_x, Number& obj_value)
{
  assert(n == obj_->Nx());

  obj_value = 0;
  obj_->Val(x, &obj_value);

  return true;
}
  
bool ipopt_opt_framework::eval_grad_f(Ipopt::Index n, const Number* x, bool new_x, Number* grad_f)
{
  assert(n == obj_->Nx());

  std::fill(grad_f, grad_f+n, 0);
  obj_->Gra(x, grad_f);

  return true;
}

bool ipopt_opt_framework::eval_g(Ipopt::Index n, const Number* x, bool new_x, Ipopt::Index m, Number* g)
{
  if ( !con_.get() )
    return false;

  std::fill(g, g+m, 0);
  con_->Val(x, g);

  return true;
}

bool ipopt_opt_framework::eval_jac_g(Ipopt::Index n, const Number* x, bool new_x,
                                     Ipopt::Index m, Ipopt::Index nele_jac, Ipopt::Index* iRow, Ipopt::Index *jCol,
                                     Number* values)
{
  if ( !con_.get() )
    return false;

  if ( values == NULL ) {
    size_t count = 0;
    for (size_t j = 0; j < J_.outerSize(); ++j) {
      for (SparseMatrix<Number>::InnerIterator it(J_, j); it; ++it) {
        iRow[count] = it.row();
        jCol[count] = it.col();
        ++count;
      }
    }
  } else {
    SparseMatrix<Number> J(m, n); {
      vector<Triplet<Number>> trips;
      con_->Jac(x, 0, &trips);
      J.setFromTriplets(trips.begin(), trips.end());
    }
    
    size_t count = 0;
    for (size_t j = 0; j < J.outerSize(); ++j) {
      for (SparseMatrix<Number>::InnerIterator it(J, j); it; ++it) {
        values[count++] = it.value();
      }
    }
  }
  
  return true;
}

bool ipopt_opt_framework::eval_h(Ipopt::Index n, const Number* x, bool new_x,
                                 Number obj_factor, Ipopt::Index m, const Number* lambda,
                                 bool new_lambda, Ipopt::Index nele_hess, Ipopt::Index* iRow,
                                 Ipopt::Index* jCol, Number* values)
{
  if (values == NULL) {
    size_t count = 0;
    for (size_t j = 0; j < lagH_.outerSize(); ++j) {
      for (SparseMatrix<Number>::InnerIterator it(lagH_, j); it; ++it) {
        if ( it.col() <= it.row() ) {
          iRow[count] = it.row();
          jCol[count] = it.col();
          ++count;
        }
      }
    }
  } else {
    if ( con_.get () ) {
      // query K
      SparseMatrix<Number> K(n, n); {
        vector<Triplet<Number>> trips;
        obj_->Hes(x, &trips);
        K.setFromTriplets(trips.begin(), trips.end());
        K *= obj_factor;
      }
      
      // query H
      vector<SparseMatrix<Number>> H(m); {
        vector<vector<Triplet<Number>>> trips(m);
        con_->Hes(x, 0, &trips);
        #pragma omp parallel for
        for (size_t i = 0; i < trips.size(); ++i) {
          H[i].resize(n, n);
          H[i].setFromTriplets(trips[i].begin(), trips[i].end());
          H[i] *= lambda[i];
        }
      }

      for (auto &mat : H)
        K += mat;

      size_t count = 0;
      for (size_t j = 0; j < K.outerSize(); ++j) {
        for (SparseMatrix<Number>::InnerIterator it(K, j); it; ++it) {
          if ( it.col() <= it.row() ) {
            values[count++] = it.value();
          }
        }
      }
    } else {
      //=> no constraints
      g_trips_.clear();
      obj_->Hes(x, &g_trips_);
      K_.resize(dim_, dim_);
      K_.setFromTriplets(g_trips_.begin(), g_trips_.end());
      K_ *= obj_factor;
      
      size_t count = 0;
      for (size_t j = 0; j < K_.outerSize(); ++j) {
        for (SparseMatrix<Number>::InnerIterator it(K_, j); it; ++it) {
          if ( it.col() <= it.row() ) {
            values[count++] = it.value();
          }
        }
      }
    }
  }
 
  return true;
}

void ipopt_opt_framework::finalize_solution(SolverReturn status,
                                            Ipopt::Index n, const Number* x, const Number* z_L, const Number* z_U,
                                            Ipopt::Index m, const Number* g, const Number* lambda,
                                            Number obj_value,
                                            const IpoptData* ip_data,
                                            IpoptCalculatedQuantities* ip_cq)
{
  std::copy(x, x+n, x0_);
}

}