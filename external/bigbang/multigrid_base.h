#ifndef EMBEDDED_FEM_FAS_H
#define EMBEDDED_FEM_FAS_H

#include <zjucad/matrix/matrix.h>
#include <boost/property_tree/ptree.hpp>
#include <Eigen/Sparse>

namespace bigbang {

using mati_t=zjucad::matrix::matrix<size_t>;
using matd_t=zjucad::matrix::matrix<double>;
using boost::property_tree::ptree;

template <typename T>
class Functional;

class FASFEMLayer
{
public:
  enum ENERGY_TYPE {
    MOMENTUM,
    ELASTIC,
    GRAVITY,
    EXTFORCE,
    FIXVERT
  };

  FASFEMLayer();
  void nlsmooth();
  void solve();

  // geometry info
  mati_t voxs_;
  matd_t nods_;
  size_t dim_;

  // physics info

  // numeric info
  std::vector<std::shared_ptr<Functional<double>>> buff_;
  std::shared_ptr<Functional<double>> energy_;
  Eigen::VectorXd rhs_, x_;
protected:
  Eigen::SparseMatrix<double> LHS_;
  Eigen::VectorXd xtmp_;
};

class FASFEMSolver
{
public:
  FASFEMSolver(const ptree &pt);
  int construct_levels();
  void direct_solve();
  void vcycle(const size_t curr);
protected:
  const ptree &pt_;
  std::vector<std::shared_ptr<FASFEMLayer>> levels_;
  std::vector<Eigen::SparseMatrix<double>> R_, P_;
};

}

#endif
