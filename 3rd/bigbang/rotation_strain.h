#ifndef ROTATION_STRAIN_H
#define ROTATION_STRAIN_H

#include <zjucad/matrix/matrix.h>

#include "def.h"

using mati_t=zjucad::matrix::matrix<size_t>;
using matd_t=zjucad::matrix::matrix<double>;

namespace bigbang {

extern "C" {
  // rs: [9x1], u: [12x1]
  void u_to_rs(double *rs, const double *u, const double *Dm);
  // jac: [9x12]
  void u_to_rs_jac(double *jac, const double *u, const double *Dm);
  // df: [9x1], x: [12x1]
  void x_to_df(double *df, const double *x, const double *Dm);
  // jac: [9x12]
  void x_to_df_jac(double *jac, const double *x, const double *Dm);
}

typedef Eigen::Matrix<double, 12, 1> vec12;
typedef Eigen::Matrix<double, 9, 1> vec9;
typedef Eigen::Matrix<double, 6, 1> vec6;
typedef Eigen::Matrix<double, 3, 1> vec3;

class rs2euc_energy : public Functional<double>
{
public:
  rs2euc_energy(const mati_t &tets, const matd_t &nods, const matd_t &u);
  size_t Nx() const;
  int Val(const double *x, double *val) const;
  int Gra(const double *x, double *gra) const;
  int Hes(const double *x, std::vector<Eigen::Triplet<double>> *hes) const;
private:
  const mati_t &tets_;
  const matd_t &nods_;
  const matd_t &u_;
  const size_t dim_;
  matd_t Dm_;
  matd_t vol_;
  Eigen::Matrix3d Id_;
};

// W: [9*#tets, 3*#vert]
void get_euc2rs_map(const mati_t &tets, const matd_t &nods, Eigen::SparseMatrix<double> &W);

}
#endif
