#include "basis_opt_energy.h"

#include <iostream>
#include <set>
#include <unsupported/Eigen/KroneckerProduct>
#include <zjucad/matrix/io.h>
#include <zjucad/matrix/itr_matrix.h>
#include <hjlib/math/blas_lapack.h>
#include <zjucad/matrix/lapack.h>
#include <Eigen/SPQRSupport>

#include "config.h"

using namespace std;
using namespace zjucad::matrix;
using namespace Eigen;
using TPL=Eigen::Triplet<double>;
using ptree=boost::property_tree::ptree;

namespace bigbang {

///===== full kronecker delta constraint =====///
full_kronecker_delta_cons::full_kronecker_delta_cons(const fine_mesh_patch &patch) {
  verts_num_ = patch.nods_.size(2);
  basis_num_ = patch.basis_pnt_.size();

  std::set<size_t> is_basis_pt(patch.basis_pnt_.begin(), patch.basis_pnt_.end());
  basis_in_strip_.resize(basis_num_, 1);

  size_t cnt = 0;
  for (size_t i = 0; i < patch.vert_strip_.size(); ++i) {
    const size_t vid = patch.vert_strip_[i];
    if ( is_basis_pt.find(vid) != is_basis_pt.end() ) {
      basis_in_strip_[cnt++] = i;
    }
  }
}
size_t full_kronecker_delta_cons::Nx() const {
  return basis_num_*verts_num_;
}
size_t full_kronecker_delta_cons::Nf() const {
  return basis_num_*basis_num_;
}
int full_kronecker_delta_cons::Val(const double *x, double *val) const {
  itr_matrix<const double*> X(basis_num_, verts_num_, x);
  itr_matrix<double *> Jc(basis_num_, basis_num_, val);
  Jc = X(colon(), basis_in_strip_)-eye<double>(basis_num_);
  return 0;
}
int full_kronecker_delta_cons::Jac(const double *x, const size_t off, std::vector<TPL> *jac) const {
  for (size_t i = 0; i < basis_num_; ++i) {
    for (size_t j = 0; j < basis_num_; ++j) {
      const size_t idx_f = i*basis_num_+j;
      const size_t idx_x = basis_in_strip_[i]*basis_num_+j;
      jac->push_back(TPL(off+idx_f, idx_x, 1.0));
    }          
  }
  return 0;
}
int full_kronecker_delta_cons::Hes(const double *x, const size_t off, vector<vector<TPL>> *hes) const {
  return __LINE__;
}

///===== full partition to unity constraint =====///
full_partition_to_unity_cons::full_partition_to_unity_cons(const fine_mesh_patch &patch) {
  verts_num_ = patch.nods_.size(2);
  basis_num_ = patch.basis_pnt_.size();
}
size_t full_partition_to_unity_cons::Nx() const {
  return basis_num_*verts_num_;
}
size_t full_partition_to_unity_cons::Nf() const {
  return verts_num_;
}
int full_partition_to_unity_cons::Val(const double *x, double *val) const {
  Map<const MatrixXd> X(x, basis_num_, verts_num_);
  for (size_t j = 0; j < X.cols(); ++j)
    val[j] = X.col(j).sum()-1;
  return 0;
}
int full_partition_to_unity_cons::Jac(const double *x, const size_t off, std::vector<TPL> *jac) const {
  for (size_t i = 0; i < verts_num_; ++i) {
    for (size_t j = 0; j < basis_num_; ++j) {
      const size_t idx = i*basis_num_+j;
      jac->push_back(TPL(off+i, idx, 1.0));
    }
  }
  return 0;
}
int full_partition_to_unity_cons::Hes(const double *x, const size_t off, vector<vector<TPL>> *hes) const {
  return __LINE__;
}

///===== full patch interpolation constraint =====///
full_patch_interp_cons::full_patch_interp_cons(const fine_mesh_patch &patch,
                                               const double power,
                                               const boost::property_tree::ptree &pt,
                                               const MatrixXd *modes)
    : patch_(patch), eps_(1e-5) {
  ASSERT(power >= 0 && power <= 1.0);

  if ( modes ) {
    //-> specified modes
    modes_ = *modes;
    wgts_ = VectorXd::Ones(modes_.cols());
  } else {
    //-> used modal modes
    //    cout << "[INTERP CONS] eigenvalues: " << patch.freqs_.transpose() << endl;
    modes_ = patch.modes_;
    
    //  wgts_.resize(patch.freqs_.size());
    wgts_ = patch.freqs_;
    //  wgts_.normalize();
    for (size_t i = 0; i < wgts_.size(); ++i)
      wgts_[i] = 1.0/pow(wgts_[i], power);

    //  const double scale = pt.get<double>("rigid_scale.value");    
    // //-> large weight for rigid modes
    // if ( patch.nods_.size(1) == 2 )
    //   wgts_[0] = 1e8; //scale*wgts_[1];
    // else if ( patch.nods_.size(1) == 3 )
    //   wgts_[0] = wgts_[1] = wgts_[2] = 1e8; //scale*wgts_[3];
    // cout << "[INTERP CONS] weights for modes: " << wgts_.transpose() << endl;
  }
  
  //  mode_num_ = std::min(mode_num, (size_t)patch.modes_.cols());
  mode_num_ = modes_.cols();  
  cout << "[INTERP CONS] used mode num: " << mode_num_ << endl;

  verts_num_ = patch.nods_.size(2);
  basis_num_ = patch.basis_pnt_.size();

  nods_sz1_ = patch.nods_.size(1);
  nods_sz2_ = patch.nods_.size(2);
  nods_sz_  = patch.nods_.size();
}
size_t full_patch_interp_cons::Nx() const {
  return basis_num_*verts_num_;
}
size_t full_patch_interp_cons::Nf() const {
  return nods_sz_*mode_num_;
}
int full_patch_interp_cons::Val(const double *x, double *val) const {
  itr_matrix<const double *> basisXY(basis_num_, verts_num_, x);
    
  for (size_t i = 0; i < mode_num_; ++i) {
    itr_matrix<const double *> ref_disp(nods_sz1_, nods_sz2_, &modes_(0, i));
    itr_matrix<double *>       resd    (nods_sz1_, nods_sz2_, val+i*nods_sz_);

    const double mode_wgt = std::sqrt(wgts_[i]);

    for (size_t q = 0; q < patch_.vert_strip_.size(); ++q) {
      const size_t pid = patch_.vert_strip_[q];

      matd_t interp_disp = zeros<double>(nods_sz1_, 1);          
      for (size_t n = 0; n < patch_.basis_pnt_.size(); ++n) {
        const size_t bid = patch_.basis_pnt_[n];              
        interp_disp += basisXY(n, q)*ref_disp(colon(), bid);
      }

      resd(colon(), pid) = (interp_disp-ref_disp(colon(), pid))*mode_wgt;
    }
  }
  return 0;
}
int full_patch_interp_cons::Jac(const double *x, const size_t off, std::vector<TPL> *jac) const {
  for (size_t i = 0; i < mode_num_; ++i) {
    itr_matrix<const double *> ref_disp(nods_sz1_, nods_sz2_, &modes_(0, i));

    const double mode_wgt = std::sqrt(wgts_[i]);
      
    for (size_t q = 0; q < patch_.vert_strip_.size(); ++q) {
      const size_t pid = patch_.vert_strip_[q];
          
      for (size_t n = 0; n < patch_.basis_pnt_.size(); ++n) {
        const size_t bid = patch_.basis_pnt_[n];

        for (size_t k = 0; k < nods_sz1_; ++k)
          jac->push_back(TPL(off+i*nods_sz_+pid*nods_sz1_+k, q*basis_num_+n, mode_wgt*ref_disp(k, bid)));
      }
    }
  }
  return 0;
}
int full_patch_interp_cons::Hes(const double *x, const size_t off, std::vector<std::vector<TPL>> *hes) const {
  return __LINE__;
}
double full_patch_interp_cons::query_interp_err(const double *x, const double *mode) const {
  itr_matrix<const double *> basisXY(basis_num_, verts_num_, x);

  itr_matrix<const double *> ref_disp(nods_sz1_, nods_sz2_, mode);
  matd_t resd = zeros<double>(nods_sz1_, nods_sz2_);

  for (size_t q = 0; q < patch_.vert_strip_.size(); ++q) {
    const size_t pid = patch_.vert_strip_[q];

    matd_t interp_disp = zeros<double>(nods_sz1_, 1);          
    for (size_t n = 0; n < patch_.basis_pnt_.size(); ++n) {
      const size_t bid = patch_.basis_pnt_[n];              
      interp_disp += basisXY(n, q)*ref_disp(colon(), bid);
    }

    resd(colon(), pid) = interp_disp-ref_disp(colon(), pid);
  }
  
  return norm(resd);
}

///===== modified houman polyharmonic energy =====///
elemwise_polyharm_energy::elemwise_polyharm_energy(const fine_mesh_patch &patch, const double w)
    : patch_(patch), w_(w), rd_(patch.nods_.size(1)) {
  basis_num_ = patch_.basis_pnt_.size();
  verts_num_ = patch_.nods_.size(2);

  mati_t permu(verts_num_, 1);
  permu(patch_.vert_strip_(colon())) = colon(0, verts_num_-1);
  
  K_.resize(patch_.energy_->Nx(), patch_.energy_->Nx()); {
    vector<Triplet<double>> trips, perm_trips;
    patch_.energy_->Hes(&patch_.REST_[0], &trips);
    for (auto &it : trips) {
      const size_t pid_i = it.row()/rd_, off_i = it.row()%rd_;
      const size_t pid_j = it.col()/rd_, off_j = it.col()%rd_;
      const size_t perm_pid_i = permu[pid_i], perm_pid_j = permu[pid_j];
      perm_trips.push_back(Triplet<double>(rd_*perm_pid_i+off_i, rd_*perm_pid_j+off_j, it.value()));
    }
    K_.setFromTriplets(perm_trips.begin(), perm_trips.end());
  }  

  G_.resize(rd_);
  for (size_t i = 0; i < rd_; ++i) {
    G_[i] = make_shared<SparseMatrix<double>>(rd_*verts_num_, verts_num_);

    vector<Triplet<double>> trips;    
    for (size_t j = 0; j < verts_num_; ++j) {
      trips.push_back(Triplet<double>(rd_*j+i, j, 1));
    }

    G_[i]->setFromTriplets(trips.begin(), trips.end());
  }
}
size_t elemwise_polyharm_energy::Nx() const {
  return basis_num_*verts_num_;
}
int elemwise_polyharm_energy::Val(const double *x, double *val) const {
  Map<const MatrixXd> X(x, basis_num_, verts_num_);

  for (size_t i = 0; i < X.rows(); ++i) {
    for (const auto &pG : G_) {
      const VectorXd p = (*pG)*X.row(i).transpose();
      *val += 0.5*w_*p.dot(K_*p);
    }
  }
  
  return 0;
}
int elemwise_polyharm_energy::Gra(const double *x, double *gra) const {
  Map<const MatrixXd> X(x, basis_num_, verts_num_);
  Map<MatrixXd> G(gra, basis_num_, verts_num_);
    
  for (size_t i = 0; i < X.rows(); ++i) {
    for (const auto &pG : G_) {
      const VectorXd p = (*pG)*X.row(i).transpose();
      G.row(i) += w_*((*pG).transpose()*K_*p).transpose();
    }
  }
  
  return 0;
}
int elemwise_polyharm_energy::Hes(const double *x, vector<Triplet<double>> *hes) const {
  SparseMatrix<double> H(Nx(), Nx());
  
  for (size_t i = 0; i < basis_num_; ++i) {
    H.setZero();
    for (const auto &pG : G_) 
      H += w_*(*pG).transpose()*K_*(*pG);

    for (size_t j = 0; j < H.outerSize(); ++j) {
      for (SparseMatrix<double>::InnerIterator it(H, j); it; ++it) {
        const size_t I = it.row()*basis_num_+i;
        const size_t J = it.col()*basis_num_+i;
        hes->push_back(Triplet<double>(I, J, it.value()));
      }
    }
  }
  return 0;
}

///===== full anistropic interpolation constraint =====///
full_anis_interp_cons::full_anis_interp_cons(const fine_mesh_patch &patch, const ptree &pt,
                                             const MatrixXd *modes)
    : patch_(patch), eps_(1e-5) {
  const double power = pt.get<double>("eig_power.value");

  if ( modes ) {
    modes_ = *modes;
    wgts_ = VectorXd::Ones(modes_.cols());
  } else {
    //-> use modal modes
    //    cout << "[ANIS INTERP CONS] eigenvalues: " << patch_.freqs_.transpose() << endl;
    modes_ = patch.modes_;

    wgts_ = patch.freqs_.cwiseInverse();
    wgts_ /= wgts_[0];
    for (size_t i = 0; i < wgts_.size(); ++i)
      wgts_[i] = pow(wgts_[i], power);
  }
  
  //-> all modes are needed
  mode_num_ = modes_.cols();
  cout << "[ANIS INTERP CONS] used mode num: " << mode_num_ << endl;
  //  cout << "[ANIS INTERP CONS] weights: " << wgts_.transpose() << endl;

  verts_num_ = patch.nods_.size(2);
  basis_num_ = patch.basis_pnt_.size();

  nods_sz1_ = patch.nods_.size(1);
  nods_sz2_ = patch.nods_.size(2);
  nods_sz_  = patch.nods_.size();
}
size_t full_anis_interp_cons::Nx() const {
  return basis_num_*verts_num_*nods_sz1_;
}
size_t full_anis_interp_cons::Nf() const {
  return nods_sz_*mode_num_;
}
int full_anis_interp_cons::Val(const double *x, double *val) const {
  itr_matrix<const double *> basisXY(basis_num_, nods_sz1_*verts_num_, x);
  
  for (size_t i = 0; i < mode_num_; ++i) {
    itr_matrix<const double *> ref_disp(nods_sz1_, nods_sz2_, &modes_(0, i));
    itr_matrix<double *>       resd    (nods_sz1_, nods_sz2_, val+i*nods_sz_);

    const double mode_wgt = std::sqrt(wgts_[i]);

    for (size_t q = 0; q < patch_.vert_strip_.size(); ++q) {
      const size_t pid = patch_.vert_strip_[q];

      matd_t interp_disp = zeros<double>(nods_sz1_, 1);          
      for (size_t n = 0; n < patch_.basis_pnt_.size(); ++n) {
        const size_t bid = patch_.basis_pnt_[n];              

        for (size_t k = 0; k < nods_sz1_; ++k) 
          interp_disp[k] += basisXY(n, k*verts_num_+q)*ref_disp(k, bid);
      }

      resd(colon(), pid) = (interp_disp-ref_disp(colon(), pid))*mode_wgt;
    }
  }
  return 0;
}
int full_anis_interp_cons::Jac(const double *x, const size_t off, vector<TPL> *jac) const {
  for (size_t i = 0; i < mode_num_; ++i) {
    itr_matrix<const double *> ref_disp(nods_sz1_, nods_sz2_, &modes_(0, i));

    const double mode_wgt = std::sqrt(wgts_[i]);
      
    for (size_t q = 0; q < patch_.vert_strip_.size(); ++q) {
      const size_t pid = patch_.vert_strip_[q];
          
      for (size_t n = 0; n < patch_.basis_pnt_.size(); ++n) {
        const size_t bid = patch_.basis_pnt_[n];

        for (size_t k = 0; k < nods_sz1_; ++k) {
          const size_t I = off+i*nods_sz_+pid*nods_sz1_+k;
          const size_t J = (k*verts_num_+q)*basis_num_+n;
          jac->push_back(TPL(I, J, mode_wgt*ref_disp(k, bid)));
        }
      }
    }
  }
  return 0;
}
int full_anis_interp_cons::Hes(const double *x, const size_t off, vector<vector<TPL>> *hes) const {
  return __LINE__;
}
void full_anis_interp_cons::reset_mode_number(const size_t number) {
  ASSERT(number >= 0 && number <= patch_.modes_.cols());
  mode_num_ = number;
}

///===== full augmented ansitropic interp constraint =====///
full_aug_anis_interp_cons::full_aug_anis_interp_cons(const fine_mesh_patch &patch, const ptree &pt,
                                                     const MatrixXd *modes, const VectorXd *weights)
    : patch_(patch) {
  const double power = pt.get<double>("eig_power.value");

  if ( modes && weights ) {
    modes_ = *modes;
    wgts_  = *weights;
  } else {
    //-> use modal modes
    cout << "[ANIS INTERP CONS] first five eigenvalues: " << patch_.freqs_.transpose().head(5) << endl;
    cout << "[ANIS INTERP CONS] last five eigenvalues: " << patch_.freqs_.transpose().tail(5) << endl;
    modes_ = patch.modes_;

    wgts_ = patch.freqs_.cwiseInverse();
    wgts_ /= wgts_[0];

    cout << "[ANIS INTERP CONS] first five weights: " << wgts_.transpose().head(5) << endl;
    cout << "[ANIS INTERP CONS] last five weights: " << wgts_.transpose().tail(5) << endl;

    std::for_each(wgts_.data(), wgts_.data()+wgts_.size(), [&power](double &x){x = pow(x, power);});
    // for (size_t i = 0; i < wgts_.size(); ++i)
    //   wgts_[i] = pow(wgts_[i], power);
  }

  //-> all modes are needed
  mode_num_ = modes_.cols();
  cout << "[ANIS INTERP CONS] used mode num: " << mode_num_ << endl;

  verts_num_ = patch.nods_.size(2);
  basis_num_ = patch.basis_pnt_.size();

  nods_sz1_ = patch.nods_.size(1);
  nods_sz2_ = patch.nods_.size(2);
  nods_sz_  = patch.nods_.size();  
}
size_t full_aug_anis_interp_cons::Nx() const {
  return basis_num_*verts_num_*(nods_sz1_*nods_sz1_);
}
size_t full_aug_anis_interp_cons::Nf() const {
  return nods_sz_*mode_num_;
}
int full_aug_anis_interp_cons::Val(const double *x, double *val) const {
  itr_matrix<const double *> basisXY(basis_num_, nods_sz1_*nods_sz1_*verts_num_, x);
  
  for (size_t i = 0; i < mode_num_; ++i) {
    itr_matrix<const double *> ref_disp(nods_sz1_, nods_sz2_, &modes_(0, i));
    itr_matrix<double *>       resd    (nods_sz1_, nods_sz2_, val+i*nods_sz_);

    const double mode_wgt = std::sqrt(wgts_[i]);

    for (size_t q = 0; q < patch_.vert_strip_.size(); ++q) {
      const size_t pid = patch_.vert_strip_[q];

      matd_t interp_disp = zeros<double>(nods_sz1_, 1);          
      for (size_t n = 0; n < patch_.basis_pnt_.size(); ++n) {
        const size_t bid = patch_.basis_pnt_[n];              

        for (size_t k = 0; k < nods_sz1_; ++k) {
          for (size_t m = 0; m < nods_sz1_; ++m) {
            interp_disp[k] += basisXY(n, (m+k*nods_sz1_)*verts_num_+q)*ref_disp(m, bid);
          }
        }
      }

      resd(colon(), pid) = (interp_disp-ref_disp(colon(), pid))*mode_wgt;
    }
  }
  return 0;
}
int full_aug_anis_interp_cons::Jac(const double *x, const size_t off, vector<TPL> *jac) const {
  for (size_t i = 0; i < mode_num_; ++i) {
    itr_matrix<const double *> ref_disp(nods_sz1_, nods_sz2_, &modes_(0, i));

    const double mode_wgt = std::sqrt(wgts_[i]);
      
    for (size_t q = 0; q < patch_.vert_strip_.size(); ++q) {
      const size_t pid = patch_.vert_strip_[q];
          
      for (size_t n = 0; n < patch_.basis_pnt_.size(); ++n) {
        const size_t bid = patch_.basis_pnt_[n];

        for (size_t k = 0; k < nods_sz1_; ++k) {
          const size_t I = off+i*nods_sz_+pid*nods_sz1_+k;
          for (size_t m = 0; m < nods_sz1_; ++m) {
            const size_t J = ((m+k*nods_sz1_)*verts_num_+q)*basis_num_+n;
            jac->push_back(TPL(I, J, mode_wgt*ref_disp(m, bid)));
          }
        }
      }
    }
  }
  return 0;
}
int full_aug_anis_interp_cons::Hes(const double *x, const size_t off, vector<vector<TPL>> *hes) const {
  return __LINE__;
}

///===== full aug anis regularization energy =====///
full_aug_anis_reg_energy::full_aug_anis_reg_energy(const fine_mesh_patch &patch, const matd_t &refs,
                                                   const double w)
    : w_(w), refs_(refs) {
  verts_num_ = patch.nods_.size(2);
  basis_num_ = patch.basis_pnt_.size();

  nods_sz1_ = patch.nods_.size(1);
  nods_sz2_ = patch.nods_.size(2);
  nods_sz_  = patch.nods_.size();  
}
size_t full_aug_anis_reg_energy::Nx() const {
  return basis_num_*verts_num_*(nods_sz1_*nods_sz1_);
}
int full_aug_anis_reg_energy::Val(const double *x, double *val) const {
  RETURN_WITH_COND_TRUE(w_ == 0.0);
  itr_matrix<const double *> X(basis_num_, Nx()/basis_num_, x);
  for (size_t i = 0; i < nods_sz1_*nods_sz1_; ++i) {
    const matd_t curr_entry = X(colon(), colon(i*verts_num_, (i+1)*verts_num_-1));
     // *val += 0.5*w_*sqnorm(curr_entry-refs_);
  }
  return 0;
}
int full_aug_anis_reg_energy::Gra(const double *x, double *gra) const {
  RETURN_WITH_COND_TRUE(w_ == 0.0);
  itr_matrix<const double *> X(basis_num_, Nx()/basis_num_, x);
  itr_matrix<double *> G(basis_num_, Nx()/basis_num_, gra);
  for (size_t i = 0; i < nods_sz1_*nods_sz1_; ++i) {
    G(colon(), colon(i*verts_num_, (i+1)*verts_num_-1))
        += w_*(X(colon(), colon(i*verts_num_, (i+1)*verts_num_-1))-refs_);
  }
  return 0;
}
int full_aug_anis_reg_energy::Hes(const double *x, vector<TPL> *hes) const {
  RETURN_WITH_COND_TRUE(w_ == 0.0);
  for (size_t i = 0; i < Nx(); ++i) {
    hes->push_back(Triplet<double>(i, i, w_));
  }
  return 0;
}

///===== full aug anis smooth energy =====///
void quad_SF_jac(double *val, const double *eps);
void vox_SF_jac(double *val, const double *eps);
full_aug_anis_smooth_energy::full_aug_anis_smooth_energy(const fine_mesh_patch &patch, const double w)
    : w_(w), rd_(patch.nods_.size(1)) {
  verts_num_ = patch.nods_.size(2);
  basis_num_ = patch.basis_pnt_.size();

  nods_sz1_ = patch.nods_.size(1);
  nods_sz2_ = patch.nods_.size(2);
  nods_sz_  = patch.nods_.size();

  const mati_t &cell = patch.cell_;
  const matd_t &nods = patch.nods_;

  //-> calculate gradient operator using one quadrature will cause 
  //-> indefinite problem, thus we use two quadratures along each axis
  const double qrs[2] = {-1.0/sqrt(3.0), 1.0/sqrt(3.0)};
  const double qrw[2] = {1.0, 1.0};
  const size_t qr_num = pow(2, rd_);
  
  vector<Triplet<double>> trips;
  for (size_t i = 0; i < cell.size(2); ++i) {
    const matd_t vert = nods(colon(), cell(colon(), i));

    matd_t Gi = zeros<double>(rd_*qr_num, cell.size(1));
    for (size_t j = 0; j < qr_num; ++j) {    
      matd_t H;
      if ( rd_ == 2 ) {
        H = zeros<double>(4, 2);
        const size_t u = j/2;
        const size_t v = j%2;
        const double xi[2] = {qrs[u], qrs[v]};
        quad_SF_jac(&H[0], xi);
      } else if ( rd_ == 3 ) {
        H = zeros<double>(8, 3);
        const size_t u = j/4;
        const size_t v = (j%4)/2;
        const size_t w = (j%4)%2;
        const double xi[3] = {qrs[u], qrs[v], qrs[w]};
        vox_SF_jac(&H[0], xi);
      }

      matd_t DmH = vert*H;
      matd_t cpDmH = DmH;
      const double vol = sqrt(fabs(det(cpDmH)));
      if ( inv(DmH) ) cerr << "# degenerated element" << endl;

      Gi(colon(rd_*j, rd_*(j+1)-1), colon()) = trans(H*DmH*vol);
    }    

    // //-> svd
    // cout << "\tGi size: " << Gi.size(1) << ", " << Gi.size(2) << endl;
    // cout << "\t all one: " << Gi*ones<double>(Gi.size(2), 1) << endl;

    // matd_t GTG = trans(Gi)*Gi;
    // matd_t U, S, VT;
    // svd(GTG, U, S, VT);
    // cout << S << endl;
    // getchar();
    
    //-> collect
    for (size_t p = 0; p < Gi.size(1); ++p) {
      for (size_t q = 0; q < Gi.size(2); ++q) {
        trips.push_back(Triplet<double>(qr_num*rd_*i+p, cell(q, i), Gi(p, q)));
      }
    }
  }
  //-> permutate the grad operator
  for (auto &it : trips) {
    it = Triplet<double>(it.row(), patch.vert_strip_perm_[it.col()], it.value());
  }
  G_.resize(qr_num*rd_*cell.size(2), nods.size(2));
  G_.reserve(trips.size());
  G_.setFromTriplets(trips.begin(), trips.end());
}
size_t full_aug_anis_smooth_energy::Nx() const {
  return basis_num_*verts_num_*(nods_sz1_*nods_sz1_);
}
int full_aug_anis_smooth_energy::Val(const double *x, double *val) const {
  RETURN_WITH_COND_TRUE(w_ == 0.0);
  Map<const MatrixXd> X(x, basis_num_, Nx()/basis_num_);
  for (size_t i = 0; i < nods_sz1_*nods_sz1_; ++i) {
    for (size_t j = 0; j < basis_num_; ++j) {
      const VectorXd f = X.block(j, i*verts_num_, 1, verts_num_).transpose();
      *val += 0.5*w_*(G_*f).squaredNorm();
    }
  }
  return 0;
}
int full_aug_anis_smooth_energy::Gra(const double *x, double *gra) const {
  RETURN_WITH_COND_TRUE(w_ == 0.0);
  Map<const MatrixXd> X(x, basis_num_, Nx()/basis_num_);
  Map<MatrixXd> G(gra, basis_num_, Nx()/basis_num_);
  for (size_t i = 0; i < nods_sz1_*nods_sz1_; ++i) {
    for (size_t j = 0; j < basis_num_; ++j) {
      const VectorXd f = X.block(j, i*verts_num_, 1, verts_num_).transpose();
      const VectorXd gf = w_*G_.transpose()*G_*f;
      G.block(j, i*verts_num_, 1, verts_num_) += gf.transpose();
    }
  }
  return 0;
}
int full_aug_anis_smooth_energy::Hes(const double *x, vector<Triplet<double>> *hes) const {
  RETURN_WITH_COND_TRUE(w_ == 0.0);
  for (size_t i = 0; i < nods_sz1_*nods_sz1_; ++i) {
    for (size_t j = 0; j < basis_num_; ++j) {
      const MatrixXd Hf = w_*G_.transpose()*G_;
      for (size_t p = 0; p < Hf.rows(); ++p) {
        for (size_t q = 0; q < Hf.cols(); ++q) {
          const size_t I = basis_num_*(i*verts_num_+p)+j;
          const size_t J = basis_num_*(i*verts_num_+q)+j;
          hes->push_back(Triplet<double>(I, J, Hf(p, q)));
        }
      }
    }    
  }
  return 0;
}

///===== full anistropic PoU constraint =====///
full_anis_PoU_cons::full_anis_PoU_cons(const fine_mesh_patch &patch) {
  verts_num_ = patch.nods_.size(2);
  basis_num_ = patch.basis_pnt_.size();
  nods_sz1_  = patch.nods_.size(1);
}
size_t full_anis_PoU_cons::Nx() const {
  return basis_num_*verts_num_*nods_sz1_;
}
size_t full_anis_PoU_cons::Nf() const {
  return verts_num_*nods_sz1_;
}
int full_anis_PoU_cons::Val(const double *x, double *val) const {
  Map<const MatrixXd> X(x, basis_num_, verts_num_*nods_sz1_);
  for (size_t j = 0; j < X.cols(); ++j)
    val[j] = X.col(j).sum()-1;
  return 0;
}
int full_anis_PoU_cons::Jac(const double *x, const size_t off, std::vector<TPL> *jac) const {
  for (size_t i = 0; i < verts_num_*nods_sz1_; ++i) {
    for (size_t j = 0; j < basis_num_; ++j) {
      const size_t idx = i*basis_num_+j;
      jac->push_back(TPL(off+i, idx, 1.0));
    }
  }
  return 0;
}
int full_anis_PoU_cons::Hes(const double *x, const size_t off, vector<vector<TPL>> *hes) const {
  return __LINE__;
}

///===== full augmented ansitropic PoU constraint =====///
full_aug_anis_PoU_cons::full_aug_anis_PoU_cons(const fine_mesh_patch &patch) {
  verts_num_ = patch.nods_.size(2);
  basis_num_ = patch.basis_pnt_.size();
  nods_sz1_  = patch.nods_.size(1);
}
size_t full_aug_anis_PoU_cons::Nx() const {
  return basis_num_*verts_num_*nods_sz1_*nods_sz1_;
}
size_t full_aug_anis_PoU_cons::Nf() const {
  return verts_num_*nods_sz1_*nods_sz1_;
}
int full_aug_anis_PoU_cons::Val(const double *x, double *val) const {
  Map<const MatrixXd> X(x, basis_num_, Nx()/basis_num_);
  for (size_t m = 0; m < nods_sz1_; ++m) {
    for (size_t n = 0; n < nods_sz1_; ++n) {
      const size_t id = m*nods_sz1_+n;
      for (size_t k = 0; k < verts_num_; ++k) {
        if ( m == n )
          val[id*verts_num_+k] = X.col(id*verts_num_+k).sum()-1;
        else
          val[id*verts_num_+k] = X.col(id*verts_num_+k).sum();
      }
    }
  }
  return 0;
}
int full_aug_anis_PoU_cons::Jac(const double *x, const size_t off, vector<TPL> *jac) const {
  for (size_t m = 0; m < nods_sz1_; ++m) {
    for (size_t n = 0; n < nods_sz1_; ++n) {
      const size_t id = m*nods_sz1_+n;
      for (size_t k = 0; k < verts_num_; ++k) {
        for (size_t n = 0; n < basis_num_; ++n) {
          jac->push_back(TPL(off+id*verts_num_+k, basis_num_*(id*verts_num_+k)+n, 1.0));
        }
      }
    }
  }
  return 0;
}
int full_aug_anis_PoU_cons::Hes(const double *x, const size_t off, vector<vector<TPL>> *hes) const {
  return __LINE__;
}

///===== full aug ansi mini perturbation =====///
extern "C" {
  void basis_mini_trace_(double *val, const double *N, const double *K);
  void basis_mini_trace_jac_(double *jac, const double *N, const double *K);
  void basis_mini_trace_hes_(double *hes, const double *N, const double *K);
}
full_aug_anis_mini_perturbation::full_aug_anis_mini_perturbation(const fine_mesh_patch &patch, const double w)
    : rd_(patch.nods_.size(1)), w_(w), patch_(patch) {
  verts_num_ = patch.nods_.size(2);
  basis_num_ = patch.basis_pnt_.size();

  nods_sz1_ = patch.nods_.size(1);
  nods_sz2_ = patch.nods_.size(2);
  nods_sz_  = patch.nods_.size();
}
size_t full_aug_anis_mini_perturbation::Nx() const {
  return basis_num_*verts_num_*(nods_sz1_*nods_sz1_);
}
int full_aug_anis_mini_perturbation::Val(const double *x, double *val) const {
  RETURN_WITH_COND_TRUE(w_ == 0.0);  
  itr_matrix<const double *> X(basis_num_, Nx()/basis_num_, x);

  for (size_t i = 0; i < patch_.vert_strip_.size(); ++i) {
    for (size_t j = 0; j < patch_.vert_strip_.size(); ++j) {
      const size_t vid_i = patch_.vert_strip_[i];
      const size_t vid_j = patch_.vert_strip_[j];
      const Matrix3d Kij = patch_.K_.block(3*vid_i, 3*vid_j, 3, 3);
      if ( Kij.norm() < 1e-8 ) // no connection between i, j
        continue;
      
      for (size_t m = 0; m < patch_.basis_pnt_.size(); ++m) {
        const size_t n = m;

        Matrix3d Nni, Nmj;
        for (size_t p = 0; p < 3; ++p) {
          for (size_t q = 0; q < 3; ++q) {
            const size_t idx = 3*p+q;
            Nni(p, q) = X(n, i+idx*verts_num_);
            Nmj(p, q) = X(m, j+idx*verts_num_);
          }
        }

        *val += 0.5*(Nni.transpose()*Kij*Nmj).trace();
      }
    }
  }
  return 0;
}
int full_aug_anis_mini_perturbation::Gra(const double *x, double *gra) const {
  RETURN_WITH_COND_TRUE(w_ == 0.0);
  itr_matrix<const double *> X(basis_num_, Nx()/basis_num_, x);
  itr_matrix<double *> G(basis_num_, Nx()/basis_num_, gra);

  for (size_t i = 0; i < patch_.vert_strip_.size(); ++i) {
    for (size_t j = 0; j < patch_.vert_strip_.size(); ++j) {
      const size_t vid_i = patch_.vert_strip_[i];
      const size_t vid_j = patch_.vert_strip_[j];
      const Matrix3d Kij = patch_.K_.block(3*vid_i, 3*vid_j, 3, 3);
      if ( Kij.norm() < 1e-8 )
        continue;
      
      for (size_t m = 0; m < patch_.basis_pnt_.size(); ++m) {
        const size_t n = m;

        Matrix<double, 3, 6> N;
        for (size_t p = 0; p < 3; ++p) {
          for (size_t q = 0; q < 3; ++q) {
            const size_t idx = 3*p+q;
            N(p, q) = X(n, i+idx*verts_num_);
            N(p, q+3) = X(m, j+idx*verts_num_);
          }
        }

        Matrix<double, 3, 6> gN;
        basis_mini_trace_jac_(gN.data(), N.data(), Kij.data());
        gN *= 0.5;

        for (size_t p = 0; p < 3; ++p) {
          for (size_t q = 0; q < 3; ++q) {
            const size_t idx = 3*p+q;
            G(n, i+idx*verts_num_) += gN(p, q);
            G(m, j+idx*verts_num_) += gN(p, q+3);
          }          
        }
      }
    }
  }
  return 0;
}
int full_aug_anis_mini_perturbation::Hes(const double *x, vector<Triplet<double>> *hes) const {
  RETURN_WITH_COND_TRUE(w_ == 0.0);

  for (size_t i = 0; i < patch_.vert_strip_.size(); ++i) {
    for (size_t j = 0; j < patch_.vert_strip_.size(); ++j) {
      const size_t vid_i = patch_.vert_strip_[i];
      const size_t vid_j = patch_.vert_strip_[j];
      const Matrix3d Kij = patch_.K_.block(3*vid_i, 3*vid_j, 3, 3);
      if ( Kij.norm() < 1e-8 )
        continue;
      
      for (size_t m = 0; m < patch_.basis_pnt_.size(); ++m) {
        const size_t n = m;

        Matrix<double, 18, 18> HN;
        basis_mini_trace_hes_(HN.data(), nullptr, Kij.data());
        HN *= 0.5;

        for (size_t w = 0; w < 18; ++w) {
          for (size_t z = 0; z < 18; ++z) {
            if ( fabs(HN(w, z)) == 0.0 )
              continue;
            size_t p, q;

            p = (w%9)/3;
            q = (w%9)%3;              
            const size_t idx_1 = 3*q+p;

            p = (z%9)/3;
            q = (z%9)%3;
            const size_t idx_2 = 3*q+p;

            size_t I, J;
              
            if ( w < 9 )
              I = (i+idx_1*verts_num_)*basis_num_+n;
            else
              I = (j+idx_1*verts_num_)*basis_num_+m;

            if ( z < 9 )
              J = (i+idx_2*verts_num_)*basis_num_+n;
            else
              J = (j+idx_2*verts_num_)*basis_num_+m;
              
            hes->push_back(Triplet<double>(I, J, HN(w, z)));
          }
        }
      }
    }
  }  
  return 0;
}

///===== full aug anis diag cons =====///
full_aug_anis_diag_cons::full_aug_anis_diag_cons(const fine_mesh_patch &patch) {
  verts_num_ = patch.nods_.size(2);
  basis_num_ = patch.basis_pnt_.size();
  nods_sz1_  = patch.nods_.size(1);

  if ( nods_sz1_ == 2 ) {
    index_.push_back(make_pair(0, 1));
    index_.push_back(make_pair(1, 0));
  } else if ( nods_sz1_ == 3 ) {
    index_.push_back(make_pair(0, 1));
    index_.push_back(make_pair(0, 2));
    index_.push_back(make_pair(1, 0));
    index_.push_back(make_pair(1, 2));
    index_.push_back(make_pair(2, 0));
    index_.push_back(make_pair(2, 1));
  } else {
    ASSERT(0);
  }
}
size_t full_aug_anis_diag_cons::Nx() const {
  return basis_num_*verts_num_*nods_sz1_*nods_sz1_;
}
size_t full_aug_anis_diag_cons::Nf() const {
  return index_.size()*basis_num_*verts_num_;
}
int full_aug_anis_diag_cons::Val(const double *x, double *val) const {
  itr_matrix<const double *> X(basis_num_, Nx()/basis_num_, x);
  itr_matrix<double *> CV(basis_num_, Nf()/basis_num_, val);
  for (size_t i = 0; i < index_.size(); ++i) {
    const size_t p = index_[i].first, q = index_[i].second;
    const size_t id = nods_sz1_*p+q;
    CV(colon(), colon(i*verts_num_, (i+1)*verts_num_-1)) = X(colon(), colon(id*verts_num_, (id+1)*verts_num_-1));
  }
  return 0;
}
int full_aug_anis_diag_cons::Jac(const double *x, const size_t off, vector<TPL> *jac) const {
  for (size_t i = 0; i < index_.size(); ++i) {
    const size_t p = index_[i].first, q = index_[i].second;
    const size_t id = nods_sz1_*p+q;
    for (size_t j = 0; j < basis_num_; ++j) {
      for (size_t k = 0; k < verts_num_; ++k) {
        jac->push_back(Triplet<double>(off+(i*verts_num_+k)*basis_num_+j, (id*verts_num_+k)*basis_num_+j, 1.0));
      }
    }
  }
  return 0;
}
int full_aug_anis_diag_cons::Hes(const double *x, const size_t off, vector<vector<TPL>> *hes) const {
  return __LINE__;
}

///===== full aug anistropic RI constraint =====///
full_aug_anis_RI_cons::full_aug_anis_RI_cons(const fine_mesh_patch &patch)
    : patch_(patch), param_nods_(patch.param_nods_) {
  verts_num_ = patch.nods_.size(2);
  basis_num_ = patch.basis_pnt_.size();
  nods_sz1_  = patch.nods_.size(1);
  ASSERT(nods_sz1_ == 3);

  cross_mat_ = zeros<double>(nods_sz1_*nods_sz1_, verts_num_);
  for (size_t i = 0; i < param_nods_.size(2); ++i) {
    const matd_t coord = param_nods_(colon(), i);
    cross_mat_(0, i) = 0;
    cross_mat_(1, i) = coord[2];
    cross_mat_(2, i) = -coord[1];
    cross_mat_(3, i) = -coord[2];
    cross_mat_(4, i) = 0;
    cross_mat_(5, i) = coord[0];
    cross_mat_(6, i) = coord[1];
    cross_mat_(7, i) = -coord[0];
    cross_mat_(8, i) = 0;
  }
}
size_t full_aug_anis_RI_cons::Nx() const {
  return basis_num_*verts_num_*nods_sz1_*nods_sz1_;
}
size_t full_aug_anis_RI_cons::Nf() const {
  return verts_num_*nods_sz1_*nods_sz1_;
}
int full_aug_anis_RI_cons::Val(const double *x, double *val) const {
  itr_matrix<const double *> X(basis_num_, Nx()/basis_num_, x);
  itr_matrix<double *> Cv(nods_sz1_*nods_sz1_, verts_num_, val);

  for (size_t i = 0; i < patch_.vert_strip_.size(); ++i) {
    const size_t vert_id = patch_.vert_strip_[i];
    
    matd_t sum_cross_mat = zeros<double>(nods_sz1_, nods_sz1_);
    matd_t N = zeros<double>(nods_sz1_, nods_sz1_);
    for (size_t j = 0; j < basis_num_; ++j) {
      N(0, 0) = X(j, 0*verts_num_+i);
      N(0, 1) = X(j, 1*verts_num_+i);
      N(0, 2) = X(j, 2*verts_num_+i);
      N(1, 0) = X(j, 3*verts_num_+i);
      N(1, 1) = X(j, 4*verts_num_+i);
      N(1, 2) = X(j, 5*verts_num_+i);
      N(2, 0) = X(j, 6*verts_num_+i);
      N(2, 1) = X(j, 7*verts_num_+i);
      N(2, 2) = X(j, 8*verts_num_+i);

      const size_t basis_id = patch_.basis_pnt_[j];
      sum_cross_mat += N*itr_matrix<const double *>(3, 3, &cross_mat_(0, basis_id));
    }

    Cv(colon(), i) = sum_cross_mat(colon())-cross_mat_(colon(), vert_id);
  }
  return 0;
}
int full_aug_anis_RI_cons::Jac(const double *x, const size_t off, vector<TPL> *jac) const {
  for (size_t i = 0; i < patch_.vert_strip_.size(); ++i) {
    const size_t vert_id = patch_.vert_strip_[i];

    for (size_t j = 0; j < basis_num_; ++j) {
      size_t N_idx[9];
      for (size_t p = 0; p < 3; ++p) {
        for (size_t q = 0; q < 3; ++q) {
          const size_t row_major_idx = 3*p+q;
          const size_t col_major_idx = 3*q+p;
          N_idx[col_major_idx] = (row_major_idx*verts_num_+i)*basis_num_+j;
        }
      }

      const size_t basis_id = patch_.basis_pnt_[j];
      for (size_t p = 0; p < 3; ++p) {
        for (size_t q = 0; q < 3; ++q) {
          const size_t idx = 3*p+q;
          if ( cross_mat_(idx, basis_id) == 0.0 )
            continue;
          jac->push_back(Triplet<double>(off+9*i+3*p+0, N_idx[3*q+0], cross_mat_(idx, basis_id)));
          jac->push_back(Triplet<double>(off+9*i+3*p+1, N_idx[3*q+1], cross_mat_(idx, basis_id)));
          jac->push_back(Triplet<double>(off+9*i+3*p+2, N_idx[3*q+2], cross_mat_(idx, basis_id)));
        }
      }
    }
  }
  return 0;
}
int full_aug_anis_RI_cons::Hes(const double *x, const size_t off, vector<vector<TPL>> *hes) const {
  return __LINE__;
}

///===== full aug anis zero bnd cons =====///
full_aug_anis_zero_bnd_cons::full_aug_anis_zero_bnd_cons(const fine_mesh_patch &patch) {
  verts_num_ = patch.nods_.size(2);
  basis_num_ = patch.basis_pnt_.size();
  nods_sz1_  = patch.nods_.size(1);

  ind_ = patch.FEM_indicator_;
  cnt_zero_ = ind_.size()-sum(ind_);
}
size_t full_aug_anis_zero_bnd_cons::Nx() const {
  return basis_num_*verts_num_*nods_sz1_*nods_sz1_;
}
size_t full_aug_anis_zero_bnd_cons::Nf() const {
  return cnt_zero_*nods_sz1_*nods_sz1_;
}
int full_aug_anis_zero_bnd_cons::Val(const double *x, double *val) const {
  for (size_t i = 0; i < nods_sz1_*nods_sz1_; ++i) {
    itr_matrix<const double *> X(basis_num_, verts_num_, x+i*basis_num_*verts_num_);
    itr_matrix<double *> Cv(cnt_zero_, 1, val+i*cnt_zero_);

    size_t cnt = 0;
    for (size_t j = 0; j < ind_.size(); ++j) {
      if ( ind_[j] == 0 ) {
        Cv[cnt] = X[j];
        ++cnt;
      }
    }
  }
  return 0;
}
int full_aug_anis_zero_bnd_cons::Jac(const double *x, const size_t off, vector<TPL> *jac) const {
  for (size_t i = 0; i < nods_sz1_*nods_sz1_; ++i) {
    size_t cnt = 0;
    for (size_t j = 0; j < ind_.size(); ++j) {
      if ( ind_[j] == 0 ) {
        jac->push_back(TPL(off+i*cnt_zero_+cnt, i*basis_num_*verts_num_+j, 1.0));
        ++cnt;
      }
    }
  }
  return 0;
}
int full_aug_anis_zero_bnd_cons::Hes(const double *x, const size_t off, vector<vector<TPL>> *hes) const {
  return __LINE__;
}

///===== full anistropic kronecker constraint =====///
full_anis_kronecker_cons::full_anis_kronecker_cons(const fine_mesh_patch &patch) {
  verts_num_ = patch.nods_.size(2);
  basis_num_ = patch.basis_pnt_.size();
  nods_sz1_  = patch.nods_.size(1);

  std::set<size_t> is_basis_pt(patch.basis_pnt_.begin(), patch.basis_pnt_.end());
  basis_in_strip_.resize(basis_num_, 1);

  size_t cnt = 0;
  for (size_t i = 0; i < patch.vert_strip_.size(); ++i) {
    const size_t vid = patch.vert_strip_[i];
    if ( is_basis_pt.find(vid) != is_basis_pt.end() ) {
      basis_in_strip_[cnt++] = i;
    }
  }
}
size_t full_anis_kronecker_cons::Nx() const {
  return basis_num_*verts_num_*nods_sz1_;
}
size_t full_anis_kronecker_cons::Nf() const {
  return basis_num_*basis_num_*nods_sz1_;
}
int full_anis_kronecker_cons::Val(const double *x, double *val) const {
  // itr_matrix<const double*> X(basis_num_, verts_num_*nods_sz1_, x);
  // itr_matrix<double *> Jc(basis_num_, basis_num_*nods_sz1_, val);
  for (size_t i = 0; i < nods_sz1_; ++i) {
    itr_matrix<const double *> X(basis_num_, verts_num_, x+i*basis_num_*verts_num_);
    itr_matrix<double *> Jc(basis_num_, basis_num_, val+i*basis_num_*basis_num_);    
    Jc = X(colon(), basis_in_strip_)-eye<double>(basis_num_);
  }
  return 0;
}
int full_anis_kronecker_cons::Jac(const double *x, const size_t off, std::vector<TPL> *jac) const {
  for (size_t k = 0; k < nods_sz1_; ++k) {
    const size_t off_f = basis_num_*basis_num_*k;
    const size_t off_x = basis_num_*verts_num_*k;
    for (size_t i = 0; i < basis_num_; ++i) {
      for (size_t j = 0; j < basis_num_; ++j) {
        const size_t idx_f = off_f+i*basis_num_+j;
        const size_t idx_x = off_x+basis_in_strip_[i]*basis_num_+j;
        jac->push_back(TPL(off+idx_f, idx_x, 1.0));
      }          
    }
  }
  return 0;
}
int full_anis_kronecker_cons::Hes(const double *x, const size_t off, vector<vector<TPL>> *hes) const {
  return __LINE__;
}

///===== full augment anis kronecker delta =====///
full_aug_anis_kronecker_cons::full_aug_anis_kronecker_cons(const fine_mesh_patch &patch) {
  verts_num_ = patch.nods_.size(2);
  basis_num_ = patch.basis_pnt_.size();
  nods_sz1_  = patch.nods_.size(1);

  std::set<size_t> is_basis_pt(patch.basis_pnt_.begin(), patch.basis_pnt_.end());
  basis_in_strip_.resize(basis_num_, 1);

  size_t cnt = 0;
  for (size_t i = 0; i < patch.vert_strip_.size(); ++i) {
    const size_t vid = patch.vert_strip_[i];
    if ( is_basis_pt.find(vid) != is_basis_pt.end() ) {
      basis_in_strip_[cnt++] = i;
    }
  } 
}
size_t full_aug_anis_kronecker_cons::Nx() const {
  return basis_num_*verts_num_*nods_sz1_*nods_sz1_;
}
size_t full_aug_anis_kronecker_cons::Nf() const {
  return basis_num_*basis_num_*nods_sz1_*nods_sz1_;
}
int full_aug_anis_kronecker_cons::Val(const double *x, double *val) const {
  for (size_t i = 0; i < nods_sz1_*nods_sz1_; ++i) {
    const size_t m = i/nods_sz1_, n = i%nods_sz1_;
    itr_matrix<const double *> X(basis_num_, verts_num_, x+i*basis_num_*verts_num_);
    itr_matrix<double *> Jc(basis_num_, basis_num_, val+i*basis_num_*basis_num_);
    if ( m == n )
      Jc = X(colon(), basis_in_strip_)-eye<double>(basis_num_);
    else
      Jc = X(colon(), basis_in_strip_);
  }
  return 0;
}
int full_aug_anis_kronecker_cons::Jac(const double *x, const size_t off, vector<TPL> *jac) const {
  for (size_t k = 0; k < nods_sz1_*nods_sz1_; ++k) {
    const size_t off_f = basis_num_*basis_num_*k;
    const size_t off_x = basis_num_*verts_num_*k;
    for (size_t i = 0; i < basis_num_; ++i) {
      for (size_t j = 0; j < basis_num_; ++j) {
        const size_t idx_f = off_f+i*basis_num_+j;
        const size_t idx_x = off_x+basis_in_strip_[i]*basis_num_+j;
        jac->push_back(TPL(off+idx_f, idx_x, 1.0));
      }          
    }
  }
  return 0;
}
int full_aug_anis_kronecker_cons::Hes(const double *x, const size_t off, vector<vector<TPL>> *hes) const {
  return __LINE__;
}

///===== full anis to iso constraint =====///
full_anis_to_iso_cons::full_anis_to_iso_cons(const fine_mesh_patch &patch) {
  verts_num_ = patch.nods_.size(2);
  basis_num_ = patch.basis_pnt_.size();
  nods_sz1_  = patch.nods_.size(1);
}
size_t full_anis_to_iso_cons::Nx() const {
  return basis_num_*verts_num_*nods_sz1_;
}
size_t full_anis_to_iso_cons::Nf() const {
  return basis_num_*verts_num_*(nods_sz1_-1);
}
int full_anis_to_iso_cons::Val(const double *x, double *val) const {
  itr_matrix<const double *> X(basis_num_, verts_num_*nods_sz1_, x);
  itr_matrix<double *> Jc(basis_num_, verts_num_*(nods_sz1_-1), val);
  for (size_t i = 0; i < nods_sz1_-1; ++i) {
    Jc(colon(), colon(i*verts_num_, (i+1)*verts_num_-1))
        = X(colon(), colon(i*verts_num_, (i+1)*verts_num_-1))
        -X(colon(), colon((i+1)*verts_num_, (i+2)*verts_num_-1));
  }
  return 0;
}
int full_anis_to_iso_cons::Jac(const double *x, const size_t off, vector<TPL> *jac) const {
  const size_t unit_num = basis_num_*verts_num_;
  for (size_t i = 0; i < nods_sz1_-1; ++i) {
    for (size_t j = 0; j < unit_num; ++j) {
      jac->push_back(TPL(off+i*unit_num+j, i*unit_num+j, 1.0));
      jac->push_back(TPL(off+i*unit_num+j, (i+1)*unit_num+j, -1.0));
    }
  }
  return 0;
}
int full_anis_to_iso_cons::Hes(const double *x, const size_t off, vector<vector<TPL>> *hes) const {
  return __LINE__;
}

// ///===== kronecker delta energy =====///
// kronecker_delta_energy::kronecker_delta_energy(const fine_quad_patch &patch, const double w)
//     : w_(w) {
//   vert_num_1d_  = patch.vert_strip_.size(2);
//   basis_num_1d_ = sqrt(patch.basis_pnt_.size());
//   const size_t stride = sqrt(patch.ONE_TO_MANY_);
//   loc_ = stride*zjucad::matrix::colon(0, basis_num_1d_-1);
//   Id_ = Eigen::MatrixXd::Identity(basis_num_1d_, basis_num_1d_);
// }
// size_t kronecker_delta_energy::Nx() const {
//   return 2*basis_num_1d_*vert_num_1d_;
// }
// int kronecker_delta_energy::Val(const double *x, double *val) const {
//   Eigen::Map<const Eigen::MatrixXd> basisX(x,        basis_num_1d_, vert_num_1d_);
//   Eigen::Map<const Eigen::MatrixXd> basisY(x+Nx()/2, basis_num_1d_, vert_num_1d_);
//   for (size_t j = 0; j < loc_.size(); ++j) {
//     *val += 0.5*w_*(basisX.col(loc_[j])-Id_.col(j)).squaredNorm();
//     *val += 0.5*w_*(basisY.col(loc_[j])-Id_.col(j)).squaredNorm();
//   }
//   return 0;
// }
// int kronecker_delta_energy::Gra(const double *x, double *gra) const {
//   Eigen::Map<const Eigen::MatrixXd> basisX(x,        basis_num_1d_, vert_num_1d_);
//   Eigen::Map<const Eigen::MatrixXd> basisY(x+Nx()/2, basis_num_1d_, vert_num_1d_);
//   Eigen::Map<Eigen::MatrixXd> gradX(gra,        basis_num_1d_, vert_num_1d_);
//   Eigen::Map<Eigen::MatrixXd> gradY(gra+Nx()/2, basis_num_1d_, vert_num_1d_);
//   for (size_t j = 0; j < loc_.size(); ++j) {
//     gradX.col(loc_[j]) += w_*(basisX.col(loc_[j])-Id_.col(j));
//     gradY.col(loc_[j]) += w_*(basisY.col(loc_[j])-Id_.col(j));
//   }
//   return 0;
// }
// int kronecker_delta_energy::Hes(const double *x, std::vector<TPL> *hes) const {
//   for (size_t i = 0; i < loc_.size(); ++i) {
//     const size_t idx = basis_num_1d_*loc_[i];
//     for (size_t j = 0; j < basis_num_1d_; ++j) {
//       hes->push_back(TPL(idx+j, idx+j, w_));
//       hes->push_back(TPL(idx+j+Nx()/2, idx+j+Nx()/2, w_));
//     }
//   }
//   return 0;
// }

// ///===== partition to unity energy =====///
// partition_to_unity_energy::partition_to_unity_energy(const fine_quad_patch &patch, const double w)
//     : w_(w) {
//   vert_num_1d_  = patch.vert_strip_.size(2);
//   basis_num_1d_ = sqrt(patch.basis_pnt_.size());
// }
// size_t partition_to_unity_energy::Nx() const {
//   return 2*vert_num_1d_*basis_num_1d_;
// }
// int partition_to_unity_energy::Val(const double *x, double *val) const {
//   Eigen::Map<const Eigen::MatrixXd> X(x, basis_num_1d_, 2*vert_num_1d_);
//   for (size_t j = 0; j < X.cols(); ++j) {
//     *val += 0.5*w_*pow(X.col(j).sum()-1, 2);
//   }
//   return 0;
// }
// int partition_to_unity_energy::Gra(const double *x, double *gra) const {
//   Eigen::Map<const Eigen::MatrixXd> X(x, basis_num_1d_, 2*vert_num_1d_);
//   Eigen::Map<Eigen::MatrixXd> G(gra, basis_num_1d_, 2*vert_num_1d_);
//   for (size_t j = 0; j < X.cols(); ++j) {
//     G.col(j) += w_*(X.col(j).sum()-1)*Eigen::VectorXd::Ones(X.rows());
//   }
//   return 0;
// }
// int partition_to_unity_energy::Hes(const double *x, std::vector<TPL> *hes) const {
//   for (size_t j = 0; j < 2*vert_num_1d_; ++j) {
//     for (size_t p = 0; p < basis_num_1d_; ++p) {
//       for (size_t q = 0; q < basis_num_1d_; ++q) {
//         const size_t I = j*basis_num_1d_+p;
//         const size_t J = j*basis_num_1d_+q;
//         hes->push_back(TPL(I, J, w_));
//       }
//     }
//   }
//   return 0;
// }

// ///===== patch interpolation constraint (2D) =====///
// patch_interp_cons_2d::patch_interp_cons_2d(const fine_quad_patch &patch, const size_t mode_num,
//                                            const int wei_option)
//     : patch_(patch), eps_(1e-5) {
//   mode_num_ = std::min(mode_num, (size_t)patch.modes_.cols());

//   vert_num_1d_    = patch.vert_strip_.size(2);
//   basis_num_1d_   = sqrt(patch.basis_pnt_.size());
//   ASSERT(basis_num_1d_*basis_num_1d_ == patch.basis_pnt_.size());
//   basis_flat_dim_ = basis_num_1d_*vert_num_1d_;

//   basis_strip_    = itr_matrix<const size_t*>(basis_num_1d_, basis_num_1d_, &patch_.basis_pnt_[0]);
//   nods_sz1_       = patch.nods_.size(1);
//   nods_sz2_       = patch.nods_.size(2);
//   nods_sz_        = patch.nods_.size();

//   wgts_ = patch.freqs_;
//   wgts_.head(3).setZero();
//   wgts_ += Eigen::VectorXd::Constant(wgts_.size(), eps_);
//   if ( wei_option == -1 )
//     wgts_ = wgts_.cwiseInverse().eval();
//   else if ( wei_option == 0 )
//     wgts_ = Eigen::VectorXd::Ones(wgts_.size());
// }
// size_t patch_interp_cons_2d::Nx() const {
//   return 2*basis_flat_dim_;
// }
// size_t patch_interp_cons_2d::Nf() const {
//   return nods_sz_*mode_num_;
// }
// int patch_interp_cons_2d::Val(const double *x, double *val) const {
//   itr_matrix<const double *> sfx(basis_num_1d_, vert_num_1d_, x);
//   itr_matrix<const double *> sfy(basis_num_1d_, vert_num_1d_, x+basis_flat_dim_);

//   //-> for each modes
//   for (size_t i = 0; i < mode_num_; ++i) {
//     itr_matrix<const double *> ref_disp(nods_sz1_, nods_sz2_, &patch_.modes_(0, i));
//     itr_matrix<double *>       resd    (nods_sz1_, nods_sz2_, val+i*nods_sz_);

//     const double mode_wgt = std::sqrt(wgts_[i]);

//     for (size_t q = 0; q < vert_num_1d_; ++q) {
//       for (size_t p = 0; p < vert_num_1d_; ++p) {
//         //-> for every fine vertices
//         const size_t pid = patch_.vert_strip_(p, q);

//         matd_t interp_disp = zeros<double>(nods_sz1_, 1);
//         for (size_t n = 0; n < basis_num_1d_; ++n) {
//           for (size_t m = 0; m < basis_num_1d_; ++m) {
//             //-> for every basis points
//             const size_t bid = basis_strip_(m, n);
              
//             interp_disp += sfx(n, q)*sfy(m, p)*ref_disp(colon(), bid);
//           }
//         }

//         resd(colon(), pid) = (interp_disp-ref_disp(colon(), pid))*mode_wgt;
//       }
//     }
//   }
//   return 0;
// }
// int patch_interp_cons_2d::Jac(const double *x, const size_t off, std::vector<TPL> *jac) const {
//   itr_matrix<const double *> sfx(basis_num_1d_, vert_num_1d_, x);
//   itr_matrix<const double *> sfy(basis_num_1d_, vert_num_1d_, x+basis_flat_dim_);

//   for (size_t i = 0; i < mode_num_; ++i) {
//     itr_matrix<const double *> ref_disp(nods_sz1_, nods_sz2_, &patch_.modes_(0, i));

//     const double mode_wgt = std::sqrt(wgts_[i]);
      
//     for (size_t q = 0; q < vert_num_1d_; ++q) {
//       for (size_t p = 0; p < vert_num_1d_; ++p) {
//         //-> for every fine vertices
//         const size_t pid = patch_.vert_strip_(p, q);
          
//         for (size_t n = 0; n < basis_num_1d_; ++n) {
//           for (size_t m = 0; m < basis_num_1d_; ++m) {
//             //-> for every basis points
//             const size_t bid = basis_strip_(m, n);
              
//             jac->push_back(TPL(off+i*nods_sz_+pid*2+0, q*basis_num_1d_+n,
//                                mode_wgt*sfy(m, p)*ref_disp(0, bid)));
//             jac->push_back(TPL(off+i*nods_sz_+pid*2+0, p*basis_num_1d_+m+Nx()/2,
//                                mode_wgt*sfx(n, q)*ref_disp(0, bid)));
//             jac->push_back(TPL(off+i*nods_sz_+pid*2+1, q*basis_num_1d_+n,
//                                mode_wgt*sfy(m, p)*ref_disp(1, bid)));
//             jac->push_back(TPL(off+i*nods_sz_+pid*2+1, p*basis_num_1d_+m+Nx()/2,
//                                mode_wgt*sfx(n, q)*ref_disp(1, bid)));
//           }
//         }
//       }
//     }
//   }
//   return 0;
// }
// int patch_interp_cons_2d::Hes(const double *x, const size_t off, vector<vector<TPL>> *jes) const {
//   return __LINE__;
// }

}
