#include "multigrid_base.h"

#include <iostream>
#include <Eigen/CholmodSupport>

#include "def.h"
#include "optimizer.h"
#include "config.h"

using namespace std;
using namespace Eigen;
using namespace zjucad::matrix;

namespace bigbang {

FASFEMLayer::FASFEMLayer() {
}

void FASFEMLayer::nlsmooth() {
  static const size_t jacobi_maxits = 0;//pt.get<size_t>("jacobi_maxits.value", 3);
  
  // linearization
  rhs_ = VectorXd::Zero(dim_); {
    energy_->Gra(x_.data(), rhs_.data());
    rhs_ *= -1;
  }
  LHS_.resize(dim_, dim_); {
    vector<Triplet<double>> trips;
    energy_->Hes(x_.data(), &trips);
    LHS_.setFromTriplets(trips.begin(), trips.end());
  }
  // damped Jacobi
  for (size_t iter = 0; iter < jacobi_maxits; ++iter)
    apply_damped_jacobi_sym(LHS_, rhs_, x_, xtmp_);
}

void FASFEMLayer::solve() {
  static const size_t direct_maxits = 0;//pt.get<size_t>("direct_maxits.value", 1000);
  static const double tolerance = 0;//pt.get<double>("tolerance.value", 1e-8);

  CholmodSimplicialLLT<SparseMatrix<double>> solver;
  VectorXd dx = VectorXd::Zero(dim_);

  for (size_t iter = 0; iter < direct_maxits; ++iter) {
    rhs_ = VectorXd::Zero(dim_); {
      energy_->Gra(x_.data(), rhs_.data());
      rhs_ *= -1;
    }
    LHS_.resize(dim_, dim_); {
      vector<Triplet<double>> trips;
      energy_->Hes(x_.data(), &trips);
      LHS_.setFromTriplets(trips.begin(), trips.end());
    }
    solver.compute(LHS_);
    ASSERT(solver.info() == Eigen::Success);
    dx = solver.solve(rhs_);
    ASSERT(solver.info() == Eigen::Success);
    x_ += dx;
  }
}

//=====================================================================

FASFEMSolver::FASFEMSolver(const ptree &pt)
    : pt_(pt) {
  
}

void FASFEMSolver::direct_solve() {
  levels_.front()->solve();
}

void FASFEMSolver::vcycle(const size_t curr) {
  const auto& curr_level = levels_[curr];
  const size_t relax_its = pt_.get<size_t>("relax_its.value", 3);
  
  if ( curr+1 == levels_.size() ) {
    curr_level->solve();
  } else {
    const auto& next_level = levels_[curr+1];
    
    // prev smooth
    for (size_t i = 0; i < relax_its; ++i)
      curr_level->nlsmooth();

    // restrict
    next_level->x_ = R_[curr]*curr_level->x_;
    VectorXd g = VectorXd::Zero(curr_level->dim_);
    curr_level->energy_->Gra(curr_level->x_.data(), g.data());
    next_level->rhs_ = R_[curr]*(curr_level->rhs_-g);

    // go to next level
    vcycle(curr+1);

    // prolongation: add interpolated error
    curr_level->x_ += P_[curr]*next_level->x_;

    // post smooth
    for (size_t i = 0; i < relax_its; ++i)
      curr_level->nlsmooth();
  }
}

}
