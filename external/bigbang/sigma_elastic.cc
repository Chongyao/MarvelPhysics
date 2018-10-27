#include "sigma_elastic.h"

#include <iostream>
#include <hjlib/math/blas_lapack.h>
#include <zjucad/matrix/lapack.h>
#include <zjucad/matrix/itr_matrix.h>
#include <zjucad/matrix/io.h>
#include <unsupported/Eigen/KroneckerProduct>
#include <jtflib/mesh/mesh.h>

#include "config.h"
#include "vox_subdivision.h"
#include "util.h"
#include "energy.h"
#include "phx_util.h"

using namespace std;
using namespace zjucad::matrix;
using namespace Eigen;

namespace bigbang {

void vox_SF(double *val, const double *epsilon);
void vox_SF_jac(double *jac, const double *epsilon);

extern "C" void vox_dfdpf_(double *jac, const double *PF, const double *U, const double *VT, const double *G);

static inline void QuadratureGradOperator(const double *HinvXH, double *Op) {
  Map<const Matrix<double, 8, 3>> Dm(HinvXH);
  Map<Matrix<double, 24, 9>> rtn(Op);
  rtn = kroneckerProduct(Dm, Matrix3d::Identity());
}

//===============================================================================
vox_sigma_elastic_energy::vox_sigma_elastic_energy(const mati_t &cube, const matd_t &nods, const matd_t &args,
                                                   const shared_ptr<cons_law_func> &psi, const double w)
    : cube_(cube), dim_(nods.size()), args_(args), w_(w), psi_(psi) {
  ASSERT(args.size(1) == psi_->num_args() && args.size(2) == cube_.size(2));
  const size_t QUAD_NUM = 8;
  
  // select eight quadrature points
  const double QU = 1.0/sqrt(3.0);
  quadrature_ = zeros<double>(3, QUAD_NUM);
  quadrature_(0, 0) = -QU; quadrature_(0, 1) = +QU; quadrature_(0, 2) = -QU; quadrature_(0, 3) = +QU;
  quadrature_(1, 0) = -QU; quadrature_(1, 1) = -QU; quadrature_(1, 2) = +QU; quadrature_(1, 3) = +QU;
  quadrature_(2, 0) = -QU; quadrature_(2, 1) = -QU; quadrature_(2, 2) = -QU; quadrature_(2, 3) = -QU;

  quadrature_(0, 4) = -QU; quadrature_(0, 5) = +QU; quadrature_(0, 6) = -QU; quadrature_(0, 7) = +QU;
  quadrature_(1, 4) = -QU; quadrature_(1, 5) = -QU; quadrature_(1, 6) = +QU; quadrature_(1, 7) = +QU;
  quadrature_(2, 4) = +QU; quadrature_(2, 5) = +QU; quadrature_(2, 6) = +QU; quadrature_(2, 7) = +QU;

  qw_ = ones<double>(QUAD_NUM, 1);

  detXH_.resize(1, QUAD_NUM*cube_.size(2));
  H_invXH_.resize(8*3, QUAD_NUM*cube_.size(2));
  GradOp_.resize(24*9, QUAD_NUM*cube_.size(2));

  #pragma omp parallel for
  for (size_t i = 0; i < cube_.size(2); ++i) {
    const matd_t X = nods(colon(), cube_(colon(), i));

    for (size_t j = 0; j < QUAD_NUM; ++j) { // for eight quadratures
      const size_t idx = QUAD_NUM*i+j;
      matd_t H = zeros<double>(8, 3);
      vox_SF_jac(&H[0], &quadrature_(0, j));

      matd_t XH = X*H;
      matd_t cpXH = XH;

      detXH_[idx] = det(cpXH);

      if ( inv(XH) )
        cerr << "# inv fail" << endl;
      matd_t HinvXH = H*XH;
      H_invXH_(colon(), idx) = HinvXH(colon());

      QuadratureGradOperator(&H_invXH_(0, idx), &GradOp_(0, idx));
    }
  }
}

size_t vox_sigma_elastic_energy::Nx() const {
  return dim_;
}

int vox_sigma_elastic_energy::Val(const double *x, double *val) const {
  RETURN_WITH_COND_TRUE(w_ == 0.0);
  itr_matrix<const double *> X(3, dim_/3, x);

  matd_t vert(3, 8), F(3, 3), S(3, 3), U(3, 3), VT(3, 3);

  for (size_t i = 0; i < cube_.size(2); ++i) {    
    vert = X(colon(), cube_(colon(), i));
    double value = 0;

    for (size_t j = 0; j < 8; ++j) {
      const size_t idx = 8*i+j;
      F = vert*itr_matrix<const double *>(8, 3, &H_invXH_(0, idx));
      svd(F, U, S, VT);
      double vr = 0;
      psi_->Val(S(0, 0), S(1, 1), S(2, 2), &args_(0, i), &vr);
      value += qw_[j]*detXH_[idx]*vr;
    }
    
    *val += w_*value;
  }
  return 0;
}

int vox_sigma_elastic_energy::Gra(const double *x, double *gra) const {
  RETURN_WITH_COND_TRUE(w_ == 0.0);
  itr_matrix<const double *> X(3, dim_/3, x);
  itr_matrix<double *> G(3, dim_/3, gra);
  
  matd_t vert(3, 8), F(3, 3), S(3, 3), U(3, 3), VT(3, 3), P(3, 3),
      jacp = zeros<double>(3, 3), g = zeros<double>(3, 8);

  for (size_t i = 0; i < cube_.size(2); ++i) {
    vert = X(colon(), cube_(colon(), i));
    g = zeros<double>(3, 8);

    for (size_t j = 0; j < 8; ++j) {
      const size_t idx = 8*i+j;
      F = vert*itr_matrix<const double *>(8, 3, &H_invXH_(0, idx));
      svd(F, U, S, VT);
      jacp = zeros<double>(3, 3);
      psi_->Jac(S(0, 0), S(1, 1), S(2, 2), &args_(0, i), &jacp[0]);
      P = U*jacp*VT;
      g(colon()) += qw_[j]*detXH_[idx]*(itr_matrix<const double *>(24, 9, &GradOp_(0, idx))*P(colon()));
    }

    G(colon(), cube_(colon(), i)) += w_*g;
  }
  return 0;
}

static inline void solve_2x2_linear_system(const double a00, const double a01,
                                           const double a10, const double a11,
                                           const double b0,  const double b1,
                                           double &x0,       double &x1) {
  const double de = a01*a10-a00*a11;
  x0 = -(a11*b0-a01*b1)/de;
  x1 = (a10*b0-a00*b1)/de;
  ASSERT(!std::isnan(x0) && !std::isnan(x1));
}

template <typename T>
zjucad::matrix::matrix<T> Diag(const zjucad::matrix::matrix<T> &A) {
  zjucad::matrix::matrix<T> rtn = zjucad::matrix::zeros<T>(A.size(1), A.size(2));
  for (size_t i = 0; i < A.size(1); ++i)
    rtn(i, i) = A(i, i);
  return rtn;
}

int vox_sigma_elastic_energy::Hes(const double *x, vector<Triplet<double>> *hes) const {
  RETURN_WITH_COND_TRUE(w_ == 0.0);
  itr_matrix<const double *> X(3, dim_/3, x);

  matd_t U(3, 3), S(3, 3), VT(3, 3);
  
  for (size_t i = 0; i < cube_.size(2); ++i) {
    matd_t vert = X(colon(), cube_(colon(), i));
    matd_t H = zeros<double>(24, 24);

    // for eight quadratures
    for (size_t j = 0; j < 8; ++j) {
      const size_t idx = 8*i+j;
      matd_t F = vert*itr_matrix<const double *>(8, 3, &H_invXH_(0, idx));

      svd(F, U, S, VT);
      matd_t jacp = zeros<double>(3, 3);
      psi_->Jac(S(0, 0), S(1, 1), S(2, 2), &args_(0, i), &jacp[0]);
      matd_t hesp = zeros<double>(3, 3);
      psi_->Hes(S(0, 0), S(1, 1), S(2, 2), &args_(0, i), &hesp[0]);

      matd_t dPdF = zeros<double>(9, 9);

      for (size_t di = 0; di < 3; ++di) {
        for (size_t dj = 0; dj < 3; ++dj) {
          const size_t col_id = 3*dj+di;

          // two anti-symmetric matrices
          matd_t omegaU = zeros<double>(3, 3), omegaVT = zeros<double>(3, 3);
          for (size_t k = 0; k < 3; ++k) {
            for (size_t l = k+1; l < 3; ++l) {
              solve_2x2_linear_system(S(l,l), S(k,k), S(k,k), S(l,l),
                                      U(di,k)*VT(l,dj), -U(di,l)*VT(k,dj),
                                      omegaU(k,l), omegaVT(k,l));
              omegaU (l, k) = -omegaU (k, l);
              omegaVT(l, k) = -omegaVT(k, l);
            }
          }
          matd_t dUdFij = U*omegaU, dVTdFij = omegaVT*VT;

          matd_t dFdFij = zeros<double>(3, 3);
          dFdFij(di, dj) = 1;
          
          matd_t dPFdFij = zeros<double>(3, 3);
          for (size_t d = 0; d < 3; ++d) {
            matd_t tmp = trans(U)*dFdFij*trans(VT);
            dPFdFij += hesp(d, d)*Diag(tmp);
          }
          
          matd_t tmp = dUdFij*jacp*VT + U*dPFdFij*VT + U*jacp*dVTdFij;
          dPdF(colon(), col_id) = tmp(colon());
        }
      }
      
      const matd_t Gg = itr_matrix<const double *>(24, 9, &GradOp_(0, idx));
      H += qw_[j]*detXH_[idx]*Gg*dPdF*trans(Gg);
    }

    for (size_t p = 0; p < 24; ++p) {
      for (size_t q = 0; q < 24; ++q) {
        const size_t I = 3*cube_(p/3, i)+p%3;
        const size_t J = 3*cube_(q/3, i)+q%3;
        hes->push_back(Triplet<double>(I, J, w_*H(p, q)));
      }
    }
  }
  return 0;
}
//===============================================================================

vox_force_matching_energy::vox_force_matching_energy(const mati_t &cube, const matd_t &nods,
                                                     const matd_t &Rx,   const matd_t &Rf, const matd_t &datw,
                                                     const shared_ptr<cons_law_func> &psi,
                                                     const vector<vector<size_t>> &fixv,
                                                     const string &metric)
    : cube_(cube), nods_(nods), Rx_(Rx), Rf_(Rf), datw_(datw), w_(1.0), psi_(psi), fixv_(fixv),
      force_eps_(1e-6), metric_(metric) {
  ASSERT(Rx.size(1) == Rf.size(1) && Rx.size(2) == Rf.size(2));
  ASSERT(fixv_.size() == Rx.size(2));
  
  args_num_ = psi_->num_args();

  const size_t QUAD_NUM = 8;
  
  // select eight quadrature points
  const double QU = 1.0/sqrt(3.0);
  quadrature_ = zeros<double>(3, QUAD_NUM);
  quadrature_(0, 0) = -QU; quadrature_(0, 1) = +QU; quadrature_(0, 2) = -QU; quadrature_(0, 3) = +QU;
  quadrature_(1, 0) = -QU; quadrature_(1, 1) = -QU; quadrature_(1, 2) = +QU; quadrature_(1, 3) = +QU;
  quadrature_(2, 0) = -QU; quadrature_(2, 1) = -QU; quadrature_(2, 2) = -QU; quadrature_(2, 3) = -QU;

  quadrature_(0, 4) = -QU; quadrature_(0, 5) = +QU; quadrature_(0, 6) = -QU; quadrature_(0, 7) = +QU;
  quadrature_(1, 4) = -QU; quadrature_(1, 5) = -QU; quadrature_(1, 6) = +QU; quadrature_(1, 7) = +QU;
  quadrature_(2, 4) = +QU; quadrature_(2, 5) = +QU; quadrature_(2, 6) = +QU; quadrature_(2, 7) = +QU;

  qw_ = ones<double>(QUAD_NUM, 1);

  detXH_.resize(1, QUAD_NUM*cube_.size(2));
  H_invXH_.resize(8*3, QUAD_NUM*cube_.size(2));
  GradOp_.resize(24*9, QUAD_NUM*cube_.size(2));

  #pragma omp parallel for
  for (size_t i = 0; i < cube_.size(2); ++i) {
    const matd_t X = nods(colon(), cube_(colon(), i));

    for (size_t j = 0; j < QUAD_NUM; ++j) { // for eight quadratures
      const size_t idx = QUAD_NUM*i+j;
      matd_t H = zeros<double>(8, 3);
      vox_SF_jac(&H[0], &quadrature_(0, j));

      matd_t XH = X*H;
      matd_t cpXH = XH;

      detXH_[idx] = det(cpXH);

      if ( inv(XH) )
        cerr << "# inv fail" << endl;
      matd_t HinvXH = H*XH;
      H_invXH_(colon(), idx) = HinvXH(colon());

      QuadratureGradOperator(&H_invXH_(0, idx), &GradOp_(0, idx));
    }
  }

  const string mtr_name = psi_->name();
  if ( mtr_name != "coro" && mtr_name != "stvk" && mtr_name != "neohookean" ) {
    cerr << "[ERR] unsupported material type for regression currently!" << endl;
    exit(1);    
  }
  matd_t tmp_params = ones<double>(psi_->num_args(), cube.size(2));
  elas_ = build_elastic_energy<voxel_elastic_potential>(mtr_name, cube, nods, tmp_params);

  solvers_.resize(Rx.size(2));
  for (size_t i = 0; i < solvers_.size(); ++i) {
    solvers_[i] = make_shared<CholmodSimplicialLDLT<SparseMatrix<double>>>();
  }
  K_.resize(Rx.size(2));
}

size_t vox_force_matching_energy::Nx() const {
  return args_num_*cube_.size(2)+1;
}

int vox_force_matching_energy::Val(const double *x, double *val) const {
  itr_matrix<const double *> X(args_num_, cube_.size(2), x);

  // update material
  dynamic_pointer_cast<voxel_elastic_potential>(elas_)->update_mtr_params(X);

  // for each data sample
  #pragma omp parallel for
  for (size_t i = 0; i < Rx_.size(2); ++i) {
    matd_t vert(3, 8), F(3, 3), S(3, 3), U(3, 3), VT(3, 3), P(3, 3),
        jacp = zeros<double>(3, 3), g = zeros<double>(3, 8);

    const matd_t nods = itr_matrix<const double *>(3, Rx_.size(1)/3, &Rx_(0, i));
    matd_t fint = zeros<double>(3, Rf_.size(1)/3);
    vector<Triplet<double>> trips;
    
    // for each voxel elastic force
    for (size_t j = 0; j < cube_.size(2); ++j) {
      vert = nods(colon(), cube_(colon(), j));
      g = zeros<double>(3, 8);

      // for each cubature points
      for (size_t k = 0; k < 8; ++k) {
        const size_t idx = 8*j+k;
        F = vert*itr_matrix<const double *>(8, 3, &H_invXH_(0, idx));
        svd(F, U, S, VT);
        jacp = zeros<double>(3, 3);
        psi_->Jac(S(0, 0), S(1, 1), S(2, 2), &X(0, j), &jacp[0]);
        P = U*jacp*VT;
        g(colon()) += qw_[k]*detXH_[idx]*(itr_matrix<const double *>(24, 9, &GradOp_(0, idx))*P(colon()));
      }

      fint(colon(), cube_(colon(), j)) += w_*g;
    }
    
    // for each fixed vertices
    for (size_t j = 0; j < fixv_[i].size(); ++j) {
      const size_t pid = fixv_[i][j], xid = Nx()-1;
      fint(colon(), pid) += x[xid]*(nods(colon(), pid)-nods_(colon(), pid));
      add_diag_block<double, 3>(pid, pid, x[xid], &trips);
    }

    const matd_t fdiff = fint(colon())-Rf_(colon(), i);
    VectorXd deltaX;

    if ( metric_ == "invK2" ) {
      elas_->Hes(&nods[0], &trips);
      SparseMatrix<double> K(elas_->Nx(), elas_->Nx());
      K.setFromTriplets(trips.begin(), trips.end());
      solvers_[i]->compute(K);
      ASSERT(!std::isnan(K.sum()));
      deltaX = solvers_[i]->solve(Map<const VectorXd>(&fdiff[0], fdiff.size()));
    } else if ( metric_ == "identity" ) {
      deltaX = Map<const VectorXd>(&fdiff[0], fdiff.size());
    } else {
      cerr << "[ERROR] not supported metric!" << endl;
      exit(EXIT_FAILURE);
    }
    
    #pragma omp critical
    {
      *val += 0.5*datw_[i]*deltaX.squaredNorm();
    }
  }
  return 0;
}

int vox_force_matching_energy::Gra(const double *x, double *gra) const {
  itr_matrix<const double *> X(args_num_, cube_.size(2), x);
  Map<VectorXd> G(gra, Nx());

  // for each data sample
  #pragma omp parallel for
  for (size_t i = 0; i < Rx_.size(2); ++i) {
    matd_t F(3, 3), S(3, 3), U(3, 3), VT(3, 3), P(3, 3);
      
    const matd_t nods = itr_matrix<const double *>(3, Rx_.size(1)/3, &Rx_(0, i));
    matd_t fint = zeros<double>(3, Rf_.size(1)/3);
    vector<Triplet<double>> trips;
    
    // for each voxel elastic force
    for (size_t j = 0; j < cube_.size(2); ++j) {
      const matd_t vert = nods(colon(), cube_(colon(), j));
      matd_t g = zeros<double>(3, 8);

      // for each cubature points
      for (size_t k = 0; k < 8; ++k) {
        const size_t idx = 8*j+k;
        F = vert*itr_matrix<const double *>(8, 3, &H_invXH_(0, idx));
        svd(F, U, S, VT);
        matd_t jacp = zeros<double>(3, 3);
        psi_->Jac(S(0, 0), S(1, 1), S(2, 2), &X(0, j), &jacp[0]);
        P = U*jacp*VT;
        g(colon()) += qw_[k]*detXH_[idx]*(itr_matrix<const double *>(24, 9, &GradOp_(0, idx))*P(colon()));

        matd_t dfdJ = zeros<double>(24, 3);
        const double PF[3] = {S(0, 0), S(1, 1), S(2, 2)};
        vox_dfdpf_(&dfdJ[0], &PF[0], &U[0], &VT[0], &GradOp_(0, idx));
        matd_t dJdp = zeros<double>(3, args_num_);
        psi_->dJdp(PF[0], PF[1], PF[2], &X(0, j), &dJdp[0]);
        matd_t dfdp = dfdJ*dJdp;
        for (size_t p = 0; p < dfdp.size(1); ++p) {
          for (size_t q = 0; q < dfdp.size(2); ++q) {
            const size_t I = 3*cube_(p/3, j)+p%3;
            const size_t J = args_num_*j+q;
            trips.push_back(Triplet<double>(I, J, qw_[k]*detXH_[idx]*dfdp(p, q)));
          }
        }
      }
      
      fint(colon(), cube_(colon(), j)) += w_*g;
    }

    // for each fixed vertices
    for (size_t j = 0; j < fixv_[i].size(); ++j) {
      const size_t pid = fixv_[i][j], xid = Nx()-1;
      fint(colon(), pid) += x[xid]*(nods(colon(), pid)-nods_(colon(), pid));
      trips.push_back(Triplet<double>(3*pid+0, xid, nods(0, pid)-nods_(0, pid)));
      trips.push_back(Triplet<double>(3*pid+1, xid, nods(1, pid)-nods_(1, pid)));
      trips.push_back(Triplet<double>(3*pid+2, xid, nods(2, pid)-nods_(2, pid)));
    }

    SparseMatrix<double> J(Rf_.size(1), Nx());
    J.reserve(trips.size());
    J.setFromTriplets(trips.begin(), trips.end());

    const matd_t fdiff = fint(colon())-Rf_(colon(), i);
    VectorXd deltaX;
    
    if ( metric_ == "invK2" ) {
      deltaX = solvers_[i]->solve(Map<const VectorXd>(&fdiff[0], fdiff.size()));
      deltaX = solvers_[i]->solve(deltaX.eval());
    } else if ( metric_ == "identity" ) {
      deltaX = Map<const VectorXd>(&fdiff[0], fdiff.size());
    } else {
      cerr << "[ERROR] not supported metric!" << endl;
      exit(EXIT_FAILURE);
    }
    
    #pragma omp critical
    {
      G += datw_[i]*J.transpose()*deltaX;
    }
  }  
  return 0;
}

int vox_force_matching_energy::Hes(const double *x, vector<Triplet<double>> *hes) const {
  return __LINE__;
}

// int vox_force_matching_energy::Hes(const double *x, vector<Triplet<double>> *hes) const {  
//   itr_matrix<const double *> X(args_num_, cube_.size(2), x);

//   matd_t F(3, 3), S(3, 3), U(3, 3), VT(3, 3), P(3, 3);

//   SparseMatrix<double> HES(Nx(), Nx());
//   HES.setZero();
  
//   // for each data sample
//   for (size_t i = 0; i < Rx_.size(2); ++i) {
//     const matd_t nods = itr_matrix<const double *>(3, Rx_.size(1)/3, &Rx_(0, i));
//     vector<Triplet<double>> trips;
    
//     // for each voxel elastic force
//     for (size_t j = 0; j < cube_.size(2); ++j) {
//       const matd_t vert = nods(colon(), cube_(colon(), j));

//       // for each cubature points
//       for (size_t k = 0; k < 8; ++k) {
//         const size_t idx = 8*j+k;
//         F = vert*itr_matrix<const double *>(8, 3, &H_invXH_(0, idx));
//         svd(F, U, S, VT);
//         matd_t jacp = zeros<double>(3, 3);
//         psi_->Jac(S(0, 0), S(1, 1), S(2, 2), &X(0, j), &jacp[0]);
//         P = U*jacp*VT;

//         matd_t dfdJ = zeros<double>(24, 3);
//         const double PF[3] = {S(0, 0), S(1, 1), S(2, 2)};
//         vox_dfdpf_(&dfdJ[0], &PF[0], &U[0], &VT[0], &GradOp_(0, idx));
//         matd_t dJdp = zeros<double>(3, args_num_);
//         psi_->dJdp(PF[0], PF[1], PF[2], &X(0, j), &dJdp[0]);
//         matd_t dfdp = dfdJ*dJdp;
//         for (size_t p = 0; p < dfdp.size(1); ++p) {
//           for (size_t q = 0; q < dfdp.size(2); ++q) {
//             const size_t I = 3*cube_(p/3, j)+p%3;
//             const size_t J = args_num_*j+q;
//             trips.push_back(Triplet<double>(I, J, qw_[k]*detXH_[idx]*dfdp(p, q)));
//           }
//         }
//       }
//     }

//     // for each fixed vertices
//     for (size_t j = 0; j < fixv_[i].size(); ++j) {
//       const size_t pid = fixv_[i][j], xid = Nx()-1;
//       trips.push_back(Triplet<double>(3*pid+0, xid, nods(0, pid)-nods_(0, pid)));
//       trips.push_back(Triplet<double>(3*pid+1, xid, nods(1, pid)-nods_(1, pid)));
//       trips.push_back(Triplet<double>(3*pid+2, xid, nods(2, pid)-nods_(2, pid)));
//     }

//     SparseMatrix<double> J(Rf_.size(1), Nx());
//     J.reserve(trips.size());
//     J.setFromTriplets(trips.begin(), trips.end());

//     if ( !ptn_normalized_ ) {
//       HES += datw_[i]*J.transpose()*J;
//     } else {
//       for (size_t r = 0; r < Rf_.size(1)/3; ++r) {
//         const double normf = norm(Rf_(colon(3*r, 3*r+2), i));
//         if ( normf < force_eps_ )
//           continue;
//         const double sqinv = 1.0/(normf*normf);
//         const auto blockJ = J.block(3*r, 0, 3, J.cols());
//         HES += datw_[i]*sqinv*blockJ.transpose()*blockJ;
//       }
//     }
//   }

//   for (size_t j = 0; j < HES.outerSize(); ++j) {
//     for (SparseMatrix<double>::InnerIterator it(HES, j); it; ++it) {
//       hes->push_back(Triplet<double>(it.row(), it.col(), it.value()));
//     }
//   }
  
//   return 0;
// }
//===============================================================================

vox_args_smooth_energy::vox_args_smooth_energy(const mati_t &cube, const matd_t &nods,
                                               const shared_ptr<cons_law_func> &psi,
                                               const double w)
    : cube_(cube), psi_(psi), w_(w) {
  args_num_ = psi->num_args();
  f2v_.reset(jtf::mesh::face2hex_adjacent::create(cube));

  const double length = norm(nods(colon(), cube(0, 0))-nods(colon(), cube(1, 0)));
  volume_ = pow(length, 3);
}

size_t vox_args_smooth_energy::Nx() const {
  return args_num_*cube_.size(2)+1;
}

int vox_args_smooth_energy::Val(const double *x, double *val) const {
  RETURN_WITH_COND_TRUE(w_ == 0.0);
  
  Map<const MatrixXd> X(x, args_num_, cube_.size(2));
  
  for (auto &it : f2v_->face2hex_) {
    if ( f2v_->is_outside_face(it) )
      continue;
    const size_t p = it.first, q = it.second;
    *val += 0.5*w_*volume_*(X.col(p)-X.col(q)).squaredNorm();
  }
  return 0;
}

int vox_args_smooth_energy::Gra(const double *x, double *gra) const {
  RETURN_WITH_COND_TRUE(w_ == 0.0);
  
  Map<const MatrixXd> X(x, args_num_, cube_.size(2));
  Map<MatrixXd> G(gra, args_num_, cube_.size(2));
  
  for (auto &it : f2v_->face2hex_) {
    if ( f2v_->is_outside_face(it) )
      continue;
    const size_t p = it.first, q = it.second;
    G.col(p) += w_*volume_*(X.col(p)-X.col(q));
    G.col(q) += w_*volume_*(X.col(q)-X.col(p));  
  }
  return 0;
}

int vox_args_smooth_energy::Hes(const double *x, vector<Triplet<double>> *hes) const {
  RETURN_WITH_COND_TRUE(w_ == 0.0);
  
  for (auto &it : f2v_->face2hex_) {
    if ( f2v_->is_outside_face(it) )
      continue;
    const size_t p = it.first, q = it.second;
    const double entry = w_*volume_;
    runtime_dim_add_diag_block(args_num_, p, p, entry, hes);
    runtime_dim_add_diag_block(args_num_, p, q, -entry, hes);
    runtime_dim_add_diag_block(args_num_, q, p, -entry, hes);
    runtime_dim_add_diag_block(args_num_, q, q, entry, hes);
  }
  return 0;
}
//===============================================================================
}
