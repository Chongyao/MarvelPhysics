#include "rotation_strain.h"

#include <algorithm>
#include <iostream>
#include <hjlib/math/blas_lapack.h>
#include <zjucad/matrix/lapack.h>
#include <unsupported/Eigen/MatrixFunctions>
#include <zjucad/matrix/itr_matrix.h>

#include "util.h"

using namespace std;
using namespace Eigen;
using namespace zjucad::matrix;

namespace bigbang {

extern "C" {
  void u_to_rs(double *rs, const double *u, const double *Dm) {
    const double
        u0 = u[0],
        u1 = u[1],
        u2 = u[2],
        u3 = u[3],
        u4 = u[4],
        u5 = u[5],
        u6 = u[6],
        u7 = u[7],
        u8 = u[8],
        u9 = u[9],
        u10 = u[10],
        u11 = u[11];
    const double
        dm0 = Dm[0],
        dm1 = Dm[1],
        dm2 = Dm[2],
        dm3 = Dm[3],
        dm4 = Dm[4],
        dm5 = Dm[5],
        dm6 = Dm[6],
        dm7 = Dm[7],
        dm8 = Dm[8];
    const double
        tt1=-u1,
        tt2=u10+tt1,
        tt3=-u0,
        tt4=u3+tt3,
        tt5=dm3*tt4,
        tt6=u4+tt1,
        tt7=u6+tt3,
        tt8=dm4*tt7,
        tt9=u7+tt1,
        tt10=u9+tt3,
        tt11=dm5*tt10,
        tt12=dm8*tt2,
        tt13=-u2,
        tt14=tt13+u11,
        tt15=dm6*tt6,
        tt16=u5+tt13,
        tt17=dm7*tt9,
        tt18=u8+tt13,
        tt19=dm6*tt4,
        tt20=dm7*tt7,
        tt21=dm8*tt10;
    const double rtn[] = {((-dm2*tt2)+tt5-dm0*tt6+tt8-dm1*tt9+tt11)/2.0,(tt12-dm5*tt14+tt15-dm3*tt16+tt17-dm4*tt18)/2.0,((-dm2*tt14)+tt19-dm0*tt16+tt20-dm1*tt18+tt21)/2.0,(2*dm0*tt4+2*dm1*tt7+2*dm2*tt10)/2.0,(2*dm5*tt2+2*dm3*tt6+2*dm4*tt9)/2.0,(2*dm8*tt14+2*dm6*tt16+2*dm7*tt18)/2.0,(tt12+dm5*tt14+tt15+dm3*tt16+tt17+dm4*tt18)/2.0,(dm2*tt14+tt19+dm0*tt16+tt20+dm1*tt18+tt21)/2.0,(dm2*tt2+tt5+dm0*tt6+tt8+dm1*tt9+tt11)/2.0};
    std::copy(rtn, rtn+9, rs);
  }
  void u_to_rs_jac(double *jac, const double *u, const double *Dm) {
    const double
        dm0 = Dm[0],
        dm1 = Dm[1],
        dm2 = Dm[2],
        dm3 = Dm[3],
        dm4 = Dm[4],
        dm5 = Dm[5],
        dm6 = Dm[6],
        dm7 = Dm[7],
        dm8 = Dm[8];
    const double
        tt1=((-dm3)-dm4-dm5)/2.0,
        tt2=((-dm6)-dm7-dm8)/2.0,
        tt3=(dm0+dm1+dm2)/2.0,
        tt4=((-dm0)-dm1-dm2)/2.0,
        tt5=dm3/2.0,
        tt6=dm6/2.0,
        tt7=-dm0/2.0,
        tt8=dm0/2.0,
        tt9=dm4/2.0,
        tt10=dm7/2.0,
        tt11=-dm1/2.0,
        tt12=dm1/2.0,
        tt13=dm5/2.0,
        tt14=dm8/2.0,
        tt15=-dm2/2.0,
        tt16=dm2/2.0;
    const double rtn[] = {tt1,0,tt2,((-2*dm0)-2*dm1-2*dm2)/2.0,0,0,0,tt2,tt1,tt3,tt2,0,0,((-2*dm3)-2*dm4-2*dm5)/2.0,0,tt2,0,tt4,0,(dm3+dm4+dm5)/2.0,tt3,0,0,((-2*dm6)-2*dm7-2*dm8)/2.0,tt1,tt4,0,tt5,0,tt6,dm0,0,0,0,tt6,tt5,tt7,tt6,0,0,dm3,0,tt6,0,tt8,0,-dm3/2.0,tt7,0,0,dm6,tt5,tt8,0,tt9,0,tt10,dm1,0,0,0,tt10,tt9,tt11,tt10,0,0,dm4,0,tt10,0,tt12,0,-dm4/2.0,tt11,0,0,dm7,tt9,tt12,0,tt13,0,tt14,dm2,0,0,0,tt14,tt13,tt15,tt14,0,0,dm5,0,tt14,0,tt16,0,-dm5/2.0,tt15,0,0,dm8,tt13,tt16,0};
    std::copy(rtn, rtn+108, jac);
  }
  void x_to_df(double *df, const double *x, const double *Dm) {
    const double
        u0 = x[0],
        u1 = x[1],
        u2 = x[2],
        u3 = x[3],
        u4 = x[4],
        u5 = x[5],
        u6 = x[6],
        u7 = x[7],
        u8 = x[8],
        u9 = x[9],
        u10 = x[10],
        u11 = x[11];
    const double
        dm0 = Dm[0],
        dm1 = Dm[1],
        dm2 = Dm[2],
        dm3 = Dm[3],
        dm4 = Dm[4],
        dm5 = Dm[5],
        dm6 = Dm[6],
        dm7 = Dm[7],
        dm8 = Dm[8];
    const double
        tt1=-u0,
        tt2=u3+tt1,
        tt3=u6+tt1,
        tt4=u9+tt1,
        tt5=-u1,
        tt6=u10+tt5,
        tt7=u4+tt5,
        tt8=u7+tt5,
        tt9=-u2,
        tt10=tt9+u11,
        tt11=u5+tt9,
        tt12=u8+tt9;
    const double rtn[]= {dm2*tt4+dm1*tt3+dm0*tt2,dm1*tt8+dm0*tt7+dm2*tt6,dm1*tt12+dm0*tt11+dm2*tt10,dm5*tt4+dm4*tt3+dm3*tt2,dm4*tt8+dm3*tt7+dm5*tt6,dm4*tt12+dm3*tt11+dm5*tt10,dm8*tt4+dm7*tt3+dm6*tt2,dm7*tt8+dm6*tt7+dm8*tt6,dm7*tt12+dm6*tt11+dm8*tt10};
    std::copy(rtn, rtn+9, df);
  }
  void x_to_df_jac(double *jac, const double *x, const double *Dm) {
    const double
        dm0 = Dm[0],
        dm1 = Dm[1],
        dm2 = Dm[2],
        dm3 = Dm[3],
        dm4 = Dm[4],
        dm5 = Dm[5],
        dm6 = Dm[6],
        dm7 = Dm[7],
        dm8 = Dm[8];
    const double
        tt1=(-dm2)-dm1-dm0,
        tt2=(-dm5)-dm4-dm3,
        tt3=(-dm8)-dm7-dm6;
    const double rtn[] = {tt1,0,0,tt2,0,0,tt3,0,0,0,tt1,0,0,tt2,0,0,tt3,0,0,0,tt1,0,0,tt2,0,0,tt3,dm0,0,0,dm3,0,0,dm6,0,0,0,dm0,0,0,dm3,0,0,dm6,0,0,0,dm0,0,0,dm3,0,0,dm6,dm1,0,0,dm4,0,0,dm7,0,0,0,dm1,0,0,dm4,0,0,dm7,0,0,0,dm1,0,0,dm4,0,0,dm7,dm2,0,0,dm5,0,0,dm8,0,0,0,dm2,0,0,dm5,0,0,dm8,0,0,0,dm2,0,0,dm5,0,0,dm8};
    std::copy(rtn, rtn+108, jac);
  }
}

template <class Vec3, class Mat3d>
static Mat3d vec3_to_skew_symm(const Vec3 &v) {
  Mat3d R;
  R(0, 0) = R(1, 1) = R(2, 2) = 0;
  R(0, 1) = v[0];
  R(1, 0) = -R(0, 1);
  R(1, 2) = v[1];
  R(2, 1) = -R(1, 2);
  R(0, 2) = v[2];
  R(2, 0) = -R(0, 2);
  return R;
}

template <class Vec6, class Mat3d>
static Mat3d vec6_to_symm(const Vec6 &v) {
  Mat3d S;
  S(0, 0) = v[0];
  S(1, 1) = v[1];
  S(2, 2) = v[2];
  S(1, 2) = v[3];
  S(0, 2) = v[4];
  S(0, 1) = v[5];
  S(2, 1) = S(1, 2);
  S(2, 0) = S(0, 2);
  S(1, 0) = S(0, 1);
  return S;
}

rs2euc_energy::rs2euc_energy(const mati_t &tets, const matd_t &nods, const matd_t &u)
    : tets_(tets), nods_(nods), u_(u), dim_(nods.size()), Id_(Matrix3d::Identity()) {
  Dm_.resize(9, tets_.size(2));
  vol_.resize(tets_.size(2));
#pragma omp parallel for
  for (size_t i = 0; i < tets.size(2); ++i) {
    matd_t dm = nods(colon(), tets(colon(1, 3), i))-nods(colon(), tets(0, i))*ones<double>(1, 3);
    matd_t tm_dm = dm;
    vol_[i] = fabs(det(tm_dm))/6;
    inv(dm);
    Dm_(colon(), i) = dm(colon());
  }
}

size_t rs2euc_energy::Nx() const {
  return dim_;
}

int rs2euc_energy::Val(const double *x, double *val) const {
  itr_matrix<const double *> X(3, Nx()/3, x);
  matd_t vert(3, 4), disp(3, 4);
  vec9 df, rs;
  for (size_t i = 0; i < tets_.size(2); ++i) {
    vert = X(colon(), tets_(colon(), i));
    disp = u_(colon(), tets_(colon(), i));
    x_to_df(df.data(), &vert[0], &Dm_(0, i));
    u_to_rs(rs.data(), &disp[0], &Dm_(0, i));
    *val += 0.5*vol_[i]*(Eigen::Map<const Matrix3d>(df.data())
                         -vec3_to_skew_symm<vec3, Matrix3d>(rs.head(3)).exp()
                         *(Id_+vec6_to_symm<vec6, Matrix3d>(rs.tail(6)))
                         ).squaredNorm();
    // *val += 0.5*vol_[i]*(Eigen::Map<const Matrix3d>(df.data())
    //                      -vec3_to_skew_symm<vec3, Matrix3d>(rs.head(3)).exp()
    //                       *vec6_to_symm<vec6, Matrix3d>(rs.tail(6)).exp()
    //                      ).squaredNorm();
  }
  return 0;
}

int rs2euc_energy::Gra(const double *x, double *gra) const {
  itr_matrix<const double *> X(3, Nx()/3, x);
  itr_matrix<double *> Grad(3, Nx()/3, gra);
  Matrix<double, 9, 12> G;
  matd_t vert(3, 4), disp(3, 4);
  vec9 df, rs;
  for (size_t i = 0; i < tets_.size(2); ++i) {
    vert = X(colon(), tets_(colon(), i));
    disp = u_(colon(), tets_(colon(), i));
    x_to_df(df.data(), &vert[0], &Dm_(0, i));
    x_to_df_jac(G.data(), nullptr, &Dm_(0, i));
    u_to_rs(rs.data(), &disp[0], &Dm_(0, i));
    // Matrix3d tm = Eigen::Map<const Matrix3d>(df.data())
    //     -vec3_to_skew_symm<vec3, Matrix3d>(rs.head(3)).exp()*vec6_to_symm<vec6, Matrix3d>(rs.tail(6)).exp();
    Matrix3d tm = Eigen::Map<const Matrix3d>(df.data())
        -vec3_to_skew_symm<vec3, Matrix3d>(rs.head(3)).exp()*(Id_+vec6_to_symm<vec6, Matrix3d>(rs.tail(6)));
    vec12 g = vol_[i]*G.transpose()*flatten<vec9, Matrix3d>(tm);
    Grad(colon(), tets_(colon(), i)) += itr_matrix<const double *>(3, 4, g.data());
  }
  return 0;
}

int rs2euc_energy:: Hes(const double *x, vector<Triplet<double>> *hes) const {
  matd_t G = zeros<double>(9, 12), H = zeros<double>(12, 12);
  for (size_t i = 0; i < tets_.size(2); ++i) {
    x_to_df_jac(&G[0], nullptr, &Dm_(0, i));
    H = vol_[i]*trans(G)*G;
    for (size_t p = 0; p < 12; ++p) {
      for (size_t q = 0; q < 12; ++q) {
        const size_t I = 3*tets_(p/3, i)+p%3;
        const size_t J = 3*tets_(q/3, i)+q%3;
        if ( H(p, q) != 0.0 )
          hes->push_back(Triplet<double>(I, J, H(p, q)));
      }
    }
  }
  return 0;
}

void get_euc2rs_map(const mati_t &tets, const matd_t &nods, SparseMatrix<double> &W) {
  matd_t Dm = zeros(9, tets.size(2));
  #pragma omp parallel for
  for (size_t i = 0; i < tets.size(2); ++i) {
    matd_t dm = nods(colon(), tets(colon(1, 3), i))-nods(colon(), tets(0, i))*ones<double>(1, 3);
    inv(dm);
    Dm(colon(), i) = dm(colon());
  }
  
  vector<Triplet<double>> trips;
  matd_t H = zeros<double>(9, 12);
  for (size_t i = 0; i < tets.size(2); ++i) {
    u_to_rs_jac(&H[0], nullptr, &Dm(0, i));
    for (size_t p = 0; p < 9; ++p) {
      for (size_t q = 0; q < 12; ++q) {
        if ( H(p, q) != 0.0 ) {
          const size_t I = 9*i+p;
          const size_t J = 3*tets(q/3, i)+q%3;
          trips.push_back(Triplet<double>(I, J, H(p, q)));
        }
      }
    }
  }
  W.resize(9*tets.size(2), 3*nods.size(2));
  W.reserve(trips.size());
  W.setFromTriplets(trips.begin(), trips.end());
}

}
