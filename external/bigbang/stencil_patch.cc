#include "stencil_patch.h"

#include <zjucad/matrix/io.h>
#include <IpIpoptApplication.hpp>
#include <IpSolveStatistics.hpp>
#include <IpIpoptCalculatedQuantities.hpp>
#include <Eigen/SparseQR>
#include <jtflib/optimizer/optimizer.h>
#include <hjlib/math/blas_lapack.h>
#include <zjucad/matrix/lapack.h>
#include <hjlib/math/polar.h>
#include <unsupported/Eigen/KroneckerProduct>
#include <Eigen/CholmodSupport>
#include <Eigen/SVD>
#include <Eigen/SPQRSupport>

#include "util.h"
#include "geom_util.h"
#include "stencil.h"
#include "energy.h"
#include "ipopt_solver.h"
#include "phx_util.h"
#include "constraint.h"
#include "io.h"
#include "optimizer.h"
#include "func_transform.h"
#include "vox_subdivision.h"
#include "arpaca.h"
#include "basis_opt_energy.h"
#include "grad_check.h"
#include "nodal_basis_optimizer.h"
#include "mass_matrix.h"

using namespace std;
using namespace zjucad::matrix;
using namespace Eigen;
using namespace Ipopt;
using boost::property_tree::ptree;

extern "C" {
  void quad4_shape_func_val_(double *val, const double *eps);
  void quad9_shape_func_val_(double *val, const double *eps);
  void hex8_shape_func_(double *val, const double *eps);
  void hex27_shape_func_(double *val, const double *eps);
}

namespace bigbang {

static shared_ptr<Functional<double>> g_func;

static void eval_val_gra(const int m, const int n, const double *x,
                         double *value, double *jacobian) {
  *value = 0;
  g_func->Val(x, value);
  std::fill(jacobian, jacobian+n, 0);
  g_func->Gra(x, jacobian);
}

static void eval_gra_hes(const int m, const int n, const double *x,
                         double *value, double *jacobian) {
  std::fill(value, value+n, 0);
  g_func->Gra(x, value);

  SparseMatrix<double> H(n, n); {
    vector<Triplet<double>> trips;
    g_func->Hes(x, &trips);
    H.setFromTriplets(trips.begin(), trips.end());
  }
  const MatrixXd Hd(H);
  std::copy(Hd.data(), Hd.data()+Hd.size(), jacobian);
}

///===== optimize anistropic basis =====///
void fine_mesh_patch::optimize_anis_basis(matd_t &basisXY) {
  //-> used to interpolate deformable modes
  shared_ptr<Constraint<double>> ip_cons;
  ip_cons = make_shared<full_anis_interp_cons>(*this, pt_);

  //-> extract the rigid modes
  MatrixXd rot_modes;
  if ( nods_.size(1) == 2 )
    rot_modes = bkp_modes_.col(2);
  else if ( nods_.size(1) == 3 )
    rot_modes = bkp_modes_.block(0, 3, bkp_modes_.rows(), 3);

  //-> used to interpolate rotation modes
  shared_ptr<Constraint<double>> re_cons;
  re_cons = make_shared<full_anis_interp_cons>(*this, pt_, &rot_modes);

  //-> two hard CONSTRAINTS
  vector<shared_ptr<Constraint<double>>> cbf(3); {
    cbf[0] = make_shared<full_anis_kronecker_cons>(*this);
    cbf[1] = make_shared<full_anis_PoU_cons>(*this);
  }
  shared_ptr<Constraint<double>> opt_con;
  try {
    opt_con = make_shared<constraint_t<double>>(cbf);
  } catch ( ... ) {
    ASSERT(0);
  }

#if 0
  do {
    VectorXd tmp_x = VectorXd::Random(opt_con->Nx());
    
    VectorXd cv = VectorXd::Zero(opt_con->Nf());
    opt_con->Val(tmp_x.data(), cv.data());

    SparseMatrix<double> Jc(opt_con->Nf(), opt_con->Nx()); {
      vector<Triplet<double>> trips;
      opt_con->Jac(nullptr, 0, &trips);
      Jc.setFromTriplets(trips.begin(), trips.end());
    }

    SparseQR<SparseMatrix<double>, COLAMDOrdering<int>> qr_solver;
    qr_solver.compute(Jc);
    const size_t rank_hes = qr_solver.rank();
    cout << "# [OPT BASIS] system size and rank: " << Jc.cols() << " " << rank_hes << endl;

    getchar();
  } while (0);
#endif

  //-> turn constraints into soft penalties
  vector<shared_ptr<Functional<double>>> ebf(3); {
    ebf[0] = make_shared<least_square_wrapper<double>>(re_cons, 1e8);
    ebf[1] = make_shared<least_square_wrapper<double>>(ip_cons, 1e0);
    ebf[2] = make_shared<least_square_wrapper<double>>(opt_con, 1e8);
  }
  shared_ptr<Functional<double>> opt_obj;
  try {
    opt_obj = make_shared<energy_t<double>>(ebf);
  } catch ( ... ) {
    ASSERT(0);
  }
  
  VectorXd basis = VectorXd::Zero(opt_obj->Nx());

#if 0
  do {
    g_func = opt_obj;
    int rtn = numeric_grad_check(eval_val_gra, 1, g_func->Nx(), &basis[0]);
    cout << "# [optimize basis] numerical check: " << rtn << endl;
    getchar();
  } while (0);
#endif

#if 0
  do {
    //-> check the rank of the hessian
    SparseMatrix<double> Hes(opt_obj->Nx(), opt_obj->Nx()); {
      vector<Triplet<double>> trips;
      opt_obj->Hes(&basis[0], &trips);
      Hes.setFromTriplets(trips.begin(), trips.end());
    }
    SparseQR<SparseMatrix<double>, COLAMDOrdering<int>> qr_solver;
    qr_solver.compute(Hes);
    const size_t rank_hes = qr_solver.rank();
    cout << "# [OPT BASIS] system size and rank: " << Hes.cols() << " " << rank_hes << endl;
  } while (0);
#endif

  double prev_value = 0; {
    opt_obj->Val(basis.data(), &prev_value);
    cout << "# prev basis opt energy: " << prev_value << endl;
  }
  VectorXd rhs = VectorXd::Zero(opt_obj->Nx()); {
    opt_obj->Gra(basis.data(), rhs.data());
    rhs *= -1;
  }
  SparseMatrix<double> LHS(opt_obj->Nx(), opt_obj->Nx()); {
    vector<Triplet<double>> trips;
    opt_obj->Hes(nullptr, &trips);
    LHS.setFromTriplets(trips.begin(), trips.end());
  }
  CholmodSimplicialLLT<SparseMatrix<double>> solver; {
    solver.compute(LHS);
    ASSERT(solver.info() == Eigen::Success);
    basis = solver.solve(rhs);
    ASSERT(solver.info() == Eigen::Success);
  }
  double post_value = 0; {
    opt_obj->Val(basis.data(), &post_value);
    cout << "# post basis opt energy: " << post_value << endl;
  }
  
  const size_t basis_num = this->basis_pnt_.size();
  basisXY = itr_matrix<const double *>(basis_num, opt_obj->Nx()/basis_num, &basis[0]);
}

void fine_mesh_patch::optimize_aug_anis_basis(matd_t &basisXY,
                                              const std::bitset<6> &flags,
                                              const matd_t *gl_har_disp,
                                              const MatrixXd *gl_eig_disp,
                                              const VectorXd *gl_eig_wgts) {
  cout << endl << "[INFO] optimize anistropic basis..." << endl;

  enum opt_hard_cons_type {
    TRANSLATION,
    ROTATION,
    KRONECKER,
    HARMONIC,
    DIAGONAL,
    VANISHING
  };
  const int NUM = VANISHING;
  vector<shared_ptr<Constraint<double>>> cbf(VANISHING+1);

  if ( flags[NUM-KRONECKER] )
    cbf[KRONECKER] = make_shared<full_aug_anis_kronecker_cons>(*this);

  if ( flags[NUM-TRANSLATION] )
    cbf[TRANSLATION] = make_shared<full_aug_anis_PoU_cons>(*this);

  if ( flags[NUM-VANISHING] )
    cbf[VANISHING] = make_shared<full_aug_anis_zero_bnd_cons>(*this);

  if ( flags[NUM-DIAGONAL] )
    cbf[DIAGONAL] = make_shared<full_aug_anis_diag_cons>(*this);

  //-> extract the rotation modes
  MatrixXd rot_modes; VectorXd rot_wgts;
  if ( nods_.size(1) == 2 ) {
    rot_modes = bkp_modes_.col(2);
    rot_wgts = VectorXd::Ones(rot_modes.cols());
    if ( flags[NUM-ROTATION] )
      cbf[ROTATION] = make_shared<full_aug_anis_interp_cons>(*this, pt_, &rot_modes, &rot_wgts);
  } else if ( nods_.size(1) == 3 ) {
    if ( flags[NUM-ROTATION] ) {
      cbf[ROTATION] = make_shared<full_aug_anis_RI_cons>(*this);
      // rot_modes =  bkp_modes_.block(0, 3, bkp_modes_.rows(), 3);
      // rot_wgts = VectorXd::Ones(rot_modes.cols());
      // cbf[ROTATION] = make_shared<full_aug_anis_interp_cons>(*this, pt_, &rot_modes, &rot_wgts);
    }
  }

  //-> used to interpolate global harmonic modes
  MatrixXd har_modes;
  if ( flags[NUM-HARMONIC] && gl_har_disp ) {
    har_modes = MatrixXd::Zero(nods_.size(), gl_har_disp->size(2));
    for (size_t i = 0; i < har_modes.cols(); ++i) {
      itr_matrix<const double *> full_mode(nods_.size(1), gl_har_disp->size(1)/nods_.size(1), &(*gl_har_disp)(0, i));
      itr_matrix<double *>(nods_.size(1), nods_.size(2), &har_modes(0, i)) = full_mode(colon(), new_to_orig_);
    }
    const VectorXd har_wgts = VectorXd::Ones(har_modes.cols());
    cbf[HARMONIC] = make_shared<full_aug_anis_interp_cons>(*this, pt_, &har_modes, &har_wgts);

    #if 0
    do {
      const matd_t curr = nods_+itr_matrix<const double *>(conf_modes.rows(), conf_modes.cols(), &conf_modes(0, 0));
      hex_mesh_write_to_vtk("./tmp_output/local_har.vtk", curr, cell_);
      getchar();
    } while (0);
    #endif
  }
  
  //-> collect all hard constraints
  shared_ptr<Constraint<double>> opt_con;
  try {
    opt_con = make_shared<constraint_t<double>>(cbf);
  } catch ( ... ) {
    ASSERT(0);
  }

  //-> find the kernel and special solution of constraints
  VectorXd us = VectorXd::Zero(opt_con->Nx());
  VectorXd cv = VectorXd::Zero(opt_con->Nf()); {
    opt_con->Val(us.data(), cv.data());
    cv *= -1;
  }
  SparseMatrix<double> Jc(opt_con->Nf(), opt_con->Nx()); {
    vector<Triplet<double>> trips;
    opt_con->Jac(nullptr, 0, &trips);
    Jc.reserve(trips.size());
    Jc.setFromTriplets(trips.begin(), trips.end());
    Jc.makeCompressed();
    cout << "\t@jacobian size: " << Jc.rows() << ", " << Jc.cols() << endl;
  }

  SPQR<decltype(Jc)> qr_solver;
  qr_solver.setPivotThreshold(1e-6);
  qr_solver.compute(Jc);
  const VectorXd spN = qr_solver.solve(cv);
  const size_t ker_dim = Jc.cols()-qr_solver.rank();
  cout << "\t@kernel dim: " << ker_dim << endl;

  if ( ker_dim == 0 ) {
    const size_t basis_num = this->basis_pnt_.size();
    basisXY = itr_matrix<const double *>(basis_num, spN.size()/basis_num, &spN[0]);  
    return;
  }

  SparseMatrix<double> KER; {
    qr_solver.compute(Jc.transpose());
    const MatrixXd Id = MatrixXd::Identity(Jc.cols(), Jc.cols());
    MatrixXd denseQ = qr_solver.matrixQ()*Id;
    denseQ = denseQ.topRightCorner(denseQ.rows(), ker_dim).eval();
    vector<Triplet<double>> trips;
    for (size_t i = 0; i < denseQ.rows(); ++i) {
      for (size_t j = 0; j < denseQ.cols(); ++j) {
        if ( fabs(denseQ(i, j)) >= 1e-8 )
          trips.push_back(Triplet<double>(i, j, denseQ(i, j)));
      }
    }
    KER.resize(denseQ.rows(), denseQ.cols());
    KER.reserve(trips.size());
    KER.setFromTriplets(trips.begin(), trips.end());
  }  
  const VectorXd tmp = VectorXd::Random(ker_dim);
  cout << "\t@test special solution: " << (Jc*spN-cv).norm() << endl;
  cout << "\t@test general solution: " << (Jc*(KER*tmp+spN)-cv).norm() << endl;

  shared_ptr<Constraint<double>> ip_cons;
  shared_ptr<Functional<double>> opt_obj;
  if ( pt_.get<string>("basis_opt_reg.value") == "spectrum" ) {
    ip_cons = make_shared<full_aug_anis_interp_cons>(*this, pt_, nullptr, nullptr);
    opt_obj = make_shared<least_square_wrapper<double>>(ip_cons, 1e0);
  } else if ( pt_.get<string>("basis_opt_reg.value") == "smooth" ) {
    opt_obj = make_shared<full_aug_anis_smooth_energy>(*this, 1e0);    
  } else if ( pt_.get<string>("basis_opt_reg.value") == "minenergy" ) {
    opt_obj = make_shared<full_aug_anis_mini_perturbation>(*this, 1e0);
  } else {
    cerr << "[ERR] inlegal regularization" << endl;
    ASSERT(0);
  }

  VectorXd q = VectorXd::Zero(ker_dim);
  VectorXd basis = KER*q+spN;

#if 0
  do {
    g_func = opt_obj;
    int rtn = numeric_grad_check(eval_val_gra, 1, g_func->Nx(), &basis[0]);
    cout << "# [optimize basis] numerical check: " << rtn << endl;
    getchar();
  } while (0);
#endif

  double prev_value = 0; {
    opt_obj->Val(basis.data(), &prev_value);
    cout << "\t@prev basis opt energy: " << prev_value << endl;
  }
  VectorXd rhs = VectorXd::Zero(ker_dim); {
    VectorXd grad = VectorXd::Zero(opt_obj->Nx());
    opt_obj->Gra(basis.data(), grad.data());
    cout << "\t@prev gradient norm: " << grad.norm() << endl;
    rhs = -KER.transpose()*grad;
  }
  SparseMatrix<double> LHS(ker_dim, ker_dim); {
    SparseMatrix<double> Hess(opt_obj->Nx(), opt_obj->Nx());
    vector<Triplet<double>> trips;
    opt_obj->Hes(nullptr, &trips);
    Hess.setFromTriplets(trips.begin(), trips.end());
    //    cout << "\t@Hessian norm: " << Hess.norm() << endl;
    LHS = KER.transpose()*Hess*KER;
    //    cout << "\t@LHS norm: " << LHS.norm() << endl;
  }
  #if 0
  {
    SPQR<decltype(LHS)> qr_solver;
    qr_solver.compute(LHS);
    cout << "# size remainigdof LHS: " << LHS.cols() << ", " << LHS.cols()-qr_solver.rank() << endl;
    
    VectorXd lhs_eigs; MatrixXd lhs_eigv;
    solve_sym_eig_prob(LHS, 0.5*LHS.cols(), 0, "arpaca", &lhs_eigs, &lhs_eigv);
    cout << "\tlhs_eigs: " << lhs_eigs.maxCoeff() << endl;

    hex_mesh_write_to_vtk("./tmp_output/A_rest_shape.vtk", nods_, cell_);

    for (size_t i = 0; i < 3; ++i) {
      const VectorXd m0 = lhs_eigv.col(0);
      const VectorXd b0 = KER*m0+spN;
      double v0 = 0;
      opt_obj->Val(b0.data(), &v0);
      cout << i << ", " << v0 << endl;
      {
        const size_t basis_num_ = basis_pnt_.size();
        const size_t rd_ = nods_.size(1);
        const size_t verts_num_ = nods_.size(2);
      
        itr_matrix<const double *> X(basis_num_, b0.size()/basis_num_, b0.data());

        matd_t uh = zeros<double>(nods_.size(1), nods_.size(2));
        matd_t uH = zeros<double>(nods_.size(1), basis_pnt_.size());
        uH(colon(), 0) = 0.01;

        for (size_t i = 0; i < this->vert_strip_.size(); ++i) {
          const size_t vid = this->vert_strip_[i];

          for (size_t j = 0; j < this->basis_pnt_.size(); ++j) {

            for (size_t k = 0; k < rd_; ++k) {
              for (size_t m = 0; m < rd_; ++m) {
                const size_t idx = k*rd_+m;
                uh(k, vid) += X(j, idx*verts_num_+i)*uH(m, j);
              }
            }
          }
        }

        matd_t curr = nods_+uh;
        char outf[256];
        sprintf(outf, "./tmp_output/A_nullN_interp_shape_%zu.vtk", i);
        hex_mesh_write_to_vtk(outf, curr, cell_);
      }
    }
    
    getchar();
    getchar();
  }
  #endif
  CholmodSimplicialLLT<SparseMatrix<double>> solver; {
    solver.compute(LHS);
    ASSERT(solver.info() == Eigen::Success);
    q = solver.solve(rhs);
    ASSERT(solver.info() == Eigen::Success);
  }
  basis = KER*q+spN;
  double post_value = 0; {
    opt_obj->Val(basis.data(), &post_value);
    cout << "\t@post basis opt energy: " << post_value << endl;
  }
  double cons_value = 0; {
    VectorXd cv = VectorXd::Zero(opt_con->Nf());
    opt_con->Val(basis.data(), cv.data());
    cons_value = cv.lpNorm<Infinity>();
    cout << "\t@hard constraints residual infnorm: " << cons_value << endl;
  }

  const size_t basis_num = this->basis_pnt_.size();
  basisXY = itr_matrix<const double *>(basis_num, opt_obj->Nx()/basis_num, &basis[0]);

#if 0
  do {
    //-> verify rotation invariant condition
    for (size_t m = 0; m < this->nods_.size(2); ++m) {
      cout << "basis on current vert " << m << endl;
      
      matd_t sumM = zeros<double>(3, 3);
      for (size_t i = 0; i < basisXY.size(1); ++i) {
        matd_t tmp(3, 3);
        tmp(0, 0) = basisXY(i, 0*this->nods_.size(2)+m);
        tmp(0, 1) = basisXY(i, 1*this->nods_.size(2)+m);
        tmp(0, 2) = basisXY(i, 2*this->nods_.size(2)+m);
        tmp(1, 0) = basisXY(i, 3*this->nods_.size(2)+m);
        tmp(1, 1) = basisXY(i, 4*this->nods_.size(2)+m);
        tmp(1, 2) = basisXY(i, 5*this->nods_.size(2)+m);
        tmp(2, 0) = basisXY(i, 6*this->nods_.size(2)+m);
        tmp(2, 1) = basisXY(i, 7*this->nods_.size(2)+m);
        tmp(2, 2) = basisXY(i, 8*this->nods_.size(2)+m);

        const matd_t pnods = param_nods_(colon(), basis_pnt_[i]);
        matd_t cross_mat = zeros<double>(3, 3);
        cross_mat(0, 1) = -pnods[2]; cross_mat(1, 0) = pnods[2];
        cross_mat(0, 2) = pnods[1]; cross_mat(2, 0) = -pnods[1];
        cross_mat(1, 2) = -pnods[0]; cross_mat(2, 1) = pnods[0];
        sumM += tmp*cross_mat;
      }

      cout << "# curr nods: " << param_nods_(colon(), vert_strip_[m]) << endl;
      cout << "# sum: " << sumM << endl;
      getchar();
    }
  } while (0);
#endif
}

///===== fine quad patch =====///
fine_quad_patch::fine_quad_patch(const quad_stencil &dom, const mati_t &quad_H,
                                 const mati_t &quad_h, const matd_t &nods_h, const matd_t &lame_h,
                                 const bool to_be_refine, const ptree &pt)
    : fine_mesh_patch(pt) {
  ONE_TO_MANY_ = quad_h.size(2)/quad_H.size(2);
  ASSERT(quad_h.size(2)%quad_H.size(2) == 0);

  //-> construct mesh patch for stencil $dom
  cell_.resize(quad_h.size(1), dom.faces_.size()*ONE_TO_MANY_);
  lame_.resize(lame_h.size(1), dom.faces_.size()*ONE_TO_MANY_);
  for (size_t i = 0; i < dom.faces_.size(); ++i) {
    const auto fid = dom.faces_[i];
    cell_(colon(), colon(i*ONE_TO_MANY_, (i+1)*ONE_TO_MANY_-1)) = quad_h(colon(), colon(fid*ONE_TO_MANY_, (fid+1)*ONE_TO_MANY_-1));
    lame_(colon(), colon(i*ONE_TO_MANY_, (i+1)*ONE_TO_MANY_-1)) = lame_h(colon(), colon(fid*ONE_TO_MANY_, (fid+1)*ONE_TO_MANY_-1));
  }
  nods_ = nods_h;
  remove_extra_verts(cell_, nods_, &new_to_orig_);

  //-> get the vertices also on the coarse mesh
  global_to_local_ = ones<size_t>(nods_h.size(2), 1)*-1;
  global_to_local_(new_to_orig_) = colon(0, new_to_orig_.size()-1);
  basis_pnt_ = global_to_local_(dom.adjc_elem_);

  //-> subdvide patch for smoother basis if neccessary
  if ( to_be_refine ) {
    // mati_t refine_cell;
    // matd_t refine_nods;
    // matd_t refine_lame;
    // subdivide_quad(cell_, nods_, refine_cell, refine_nods);
    // const size_t sbd_per_to_many = refine_cell.size(2)/cell_.size(2);
    // refine_lame.resize(lame_.size(1), refine_cell.size(2));
    // for (size_t i = 0; i < lame_.size(2); ++i)
    //   refine_lame(colon(), colon(i*sbd_per_to_many, (i+1)*sbd_per_to_many-1)) = lame_(colon(), i)*ones<double>(1, sbd_per_to_many);
    // cell_ = refine_cell;
    // nods_ = refine_nods;
    // lame_ = refine_lame;
    // ONE_TO_MANY_ *= sbd_per_to_many;
  }

  //-> record the rest shape
  REST_ = nods_;

  //-> init edge2cell
  e2c_.reset(jtf::mesh::edge2cell_adjacent::create(cell_, false));
  
  //-> extract vertex strip
  strip_ext_ = make_shared<straight_quad_strip_extractor>(cell_);
  this->extract_vert_strip();

  //-> calc vert strip permutation and basis points permutation
  vert_strip_perm_ = -1*ones<size_t>(nods_.size(2), 1);
  vert_strip_perm_(vert_strip_(colon())) = colon(0, vert_strip_.size()-1);
  basis_pnt_perm_ = -1*ones<size_t>(nods_.size(2), 1);
  basis_pnt_perm_(basis_pnt_(colon())) = colon(0, basis_pnt_.size()-1);

  //-> init elastic energy
  energy_ = build_elastic_energy<quad4_elastic_energy>(pt.get<string>("elas.value"), cell_, REST_, lame_); {
    vector<Triplet<double>> trips;
    energy_->Hes(&REST_[0], &trips);
    K_.resize(energy_->Nx(), energy_->Nx());
    K_.setFromTriplets(trips.begin(), trips.end());
  }

  #if 0
  quad_mesh_write_to_vtk("./tmp_output/debug_patch.vtk", nods_, cell_);
  #endif
}

void fine_quad_patch::extract_vert_strip() {
  const size_t START_FACE = 0;
  vector<size_t> start_faces;
  strip_ext_->extract(START_FACE, 'X', start_faces);

  vector<vector<size_t>> hori_face_strip(start_faces.size());
  for (size_t i = 0; i < start_faces.size(); ++i)
    strip_ext_->extract(start_faces[i], 'Y', hori_face_strip[i]);

  vert_strip_ = zeros<size_t>(start_faces.size()+1, start_faces.size()+1);
  for (size_t i = 0; i < hori_face_strip.size(); ++i) {
    for (size_t j = 0; j < hori_face_strip[i].size(); ++j) {
      vert_strip_(j, i) = cell_(0, hori_face_strip[i][j]);
      if ( j+1 == hori_face_strip[i].size() )
        vert_strip_(j+1, i) = cell_(3, hori_face_strip[i][j]);
      if ( i+1 == hori_face_strip.size() )
        vert_strip_(j, i+1) = cell_(1, hori_face_strip[i][j]);
      if ( j+1 == hori_face_strip[i].size() && i+1 == hori_face_strip.size() )
        vert_strip_(j+1, i+1) = cell_(2, hori_face_strip[i][j]);
    }
  }
}

void fine_quad_patch::modal_analysis() {
  ASSERT(energy_.get());
  SparseMatrix<double> K(energy_->Nx(), energy_->Nx()); {
    vector<Triplet<double>> trips;
    energy_->Hes(&REST_[0], &trips);
    K.setFromTriplets(trips.begin(), trips.end());
  }
  SparseMatrix<double> M;
  calc_surf_mass_matrix(cell_, REST_, 1.0, &M);
  
  const size_t eignum = K.cols()*0.8;
  
  const int smallest = 0;
  VectorXd eigenvalues; MatrixXd eigenvectors;
  int eig_rtn = solve_gen_eig_prob(K, M, eignum, smallest, "arpaca", &eigenvalues, &eigenvectors);
  ASSERT(eig_rtn == 0);
  cout << "[QUAD LMA] norm of head 3 lambda: " << eigenvalues.head(3).norm() << endl;

  bkp_freqs_ = eigenvalues;
  bkp_modes_ = eigenvectors;
  
  //-> separate translation and infinitesimal rotation
  VectorXd trans_x = VectorXd::Zero(bkp_modes_.rows()), trans_y = trans_x;
  for (size_t i = 0; i < trans_x.size()/2; ++i)
    trans_x[2*i+0] = trans_y[2*i+1] = 1;
  normalize_vec_metric(trans_x, M);
  normalize_vec_metric(trans_y, M);
  
  const MatrixXd kernel_K = bkp_modes_.topLeftCorner(bkp_modes_.rows(), 3);
  MatrixXd LHS = MatrixXd::Zero(2, 3); {
    LHS.row(0) = trans_x.transpose()*M*kernel_K;
    LHS.row(1) = trans_y.transpose()*M*kernel_K;
  }
  FullPivLU<MatrixXd> lu(LHS);
  ASSERT(lu.kernel().cols() == 1);
  VectorXd inf_rot_basis = kernel_K*lu.kernel();
  normalize_vec_metric(inf_rot_basis, M);

  bkp_modes_.col(0) = trans_x;
  bkp_modes_.col(1) = trans_y;
  bkp_modes_.col(2) = inf_rot_basis;

  //-> only preserve deformable bases
  const size_t non_rigid_mode_num = bkp_freqs_.size()-3;
  freqs_ = bkp_freqs_.tail(non_rigid_mode_num);
  modes_ = bkp_modes_.topRightCorner(bkp_modes_.rows(), non_rigid_mode_num);

#if 0
  do {
    const string outdir = pt_.get<string>("outdir.value");

    //-> write rotational basis
    static size_t count = 0;
    if ( count == 0 ) {
      char outfile[256];
      for (size_t i = 3; i < 6; ++i) {
        const matd_t disp = itr_matrix<const double *>(nods_.size(1), nods_.size(2), &bkp_modes_(0, i));
        const matd_t curr = REST_+0.001*disp;
        sprintf(outfile, "%s/ROT_patch%04zu_mode_%03zu.vtk", outdir.c_str(), count, i);
        quad_mesh_write_to_vtk(outfile, REST_, cell_, &disp, "POINT");
      }
    }
    ++count;
  } while (0);
#endif

#if 0
  cout << "# freqs: " << bkp_freqs_.transpose() << endl;
  do {
    //-> write modal basis
    static size_t count = 0;

    char outfile[256];
    for (size_t i = 0; i < 10; ++i) {
      const matd_t disp = itr_matrix<const double*>(nods_.size(1), nods_.size(2), &bkp_modes_(0, i));
      const matd_t def_nods = nods_+disp;

      sprintf(outfile, "./tmp_output/new_quad_patch%04zu_origin_mode%03zu.vtk", count, i);
      quad_mesh_write_to_vtk(outfile, def_nods, cell_, &lame_, "CELL");
    }
    ++count;
  } while (0);
#endif

#if 0
  do {
    static int count = 0;
    cout << "\n\n\n\nCOUNT: " << count++ << endl;
    const MatrixXd KK(K);
    cout << "K norm: " << K.norm() << endl;
    cout << "K sum: " << K.sum() << endl;
    cout << "K blue norm: " << K.blueNorm() << endl;
    cout << "K max and min: " << KK.minCoeff() << " " << KK.maxCoeff() << endl;
    cout << "# modes norm as identification: " << modes_.norm() << endl;
    cout << "# modes lp norm: " << modes_.lpNorm<Infinity>() << endl;
    cout << "# modes l1 norm: " << modes_.lpNorm<1>() << endl;
    cout << "freqs: " << freqs_.transpose() << endl;
    cout << "# UTKU-Lambda norm: " << (modes_.transpose()*K*modes_-MatrixXd(freqs_.asDiagonal())).norm() << endl;
    cout << "# UTMU-Id norm: " << (modes_.transpose()*M*modes_-MatrixXd::Identity(modes_.cols(), modes_.cols())).norm() << endl;
    //    cout << modes_.topLeftCorner(modes_.rows(), 3) << endl;
    getchar();
  } while (0);
#endif
}

// void fine_quad_patch::optimize_basis_houman(const SparseMatrix<double> &N, matd_t &basisXY) {
//   const size_t basis_num = this->basis_pnt_.size(), sample_num = this->nods_.size(2);
//   basisXY = zeros<double>(basis_num, sample_num);

//   for (size_t i = 0; i < basisXY.size(1); ++i) {
//     for (size_t j = 0; j < basisXY.size(2); ++j) {
//       const size_t basis_pt  = new_to_orig_[this->basis_pnt_[i]];
//       const size_t sample_pt = new_to_orig_[this->vert_strip_[j]];
//       basisXY(i, j) = N.coeff(basis_pt, sample_pt);
//     }
//   }
// }

void fine_quad_patch::optimize_basis(matd_t &basisXY) {
  ASSERT(0);
//   //-> remove first two translational bases
//   // modes_ = modes_.topRightCorner(modes_.rows(), modes_.cols()-2).eval();
//   // freqs_ = freqs_.tail(freqs_.size()-2).eval();

//   const size_t mode_num = basis_pnt_.size()+pt_.get<int>("eigmode_plus.value");
//   const double power = pt_.get<double>("eig_power.value");  

//   shared_ptr<Constraint<double>> ip_cons;
//   ip_cons = make_shared<full_patch_interp_cons>(*this, power, pt_);

//   //-> rest shape reconstruction constraint
//   const MatrixXd rest_shape = Eigen::Map<const MatrixXd>(&REST_[0], REST_.size(), 1);
//   shared_ptr<Constraint<double>> re_cons;
//   re_cons = make_shared<full_patch_interp_cons>(*this, power, pt_, &rest_shape);

//   //-> CONSTRAINTS
//   vector<shared_ptr<Constraint<double>>> cbf(3); {
//     cbf[0] = make_shared<full_partition_to_unity_cons>(*this);
//     cbf[1] = make_shared<full_kronecker_delta_cons>(*this);
//   }
//   shared_ptr<Constraint<double>> opt_con;
//   try {
//     opt_con = make_shared<constraint_t<double>>(cbf);
//   } catch ( ... ) {
//     ASSERT(0);
//   }

//   //-> OBJECTIVE
//   vector<shared_ptr<Functional<double>>> ebf(3); {
//     ebf[0] = make_shared<least_square_wrapper<double>>(ip_cons, 1e-1);
//     ebf[1] = make_shared<least_square_wrapper<double>>(re_cons, 1e8);
//     // ebf[1] = make_shared<least_square_wrapper<double>>(cbf[0], w_pu);
//     // ebf[2] = make_shared<least_square_wrapper<double>>(cbf[1], 1e4);
//   }
//   shared_ptr<Functional<double>> opt_obj;
//   try {
//     opt_obj = make_shared<energy_t<double>>(ebf);
//   } catch ( ... ) {
//     ASSERT(0);
//   }

//   ASSERT(opt_obj->Nx() == opt_con->Nx());
//   matd_t basis = ones<double>(opt_obj->Nx(), 1);

// #if 0
//   do {
//     g_func = opt_obj;
//     int rtn = numeric_grad_check(eval_val_gra, 1, g_func->Nx(), &basis[0]);
//     cout << "# [optimize basis] numerical check: " << rtn << endl;
//     getchar();
//   } while (0);
// #endif

// #if 0
//   do {
//     //-> check the rank of the hessian
//     SparseMatrix<double> Hes(opt_obj->Nx(), opt_obj->Nx()); {
//       vector<Triplet<double>> trips;
//       opt_obj->Hes(&basis[0], &trips);
//       Hes.setFromTriplets(trips.begin(), trips.end());
//     }
//     // SparseMatrix<double> Jac(opt_con->Nf(), opt_con->Nx()); {
//     //   vector<Triplet<double>> trips;
//     //   opt_con->Jac(&basis[0], 0, &trips);
//     //   Jac.setFromTriplets(trips.begin(), trips.end());
//     // }
    
//     SparseQR<SparseMatrix<double>, COLAMDOrdering<int>> qr_solver;
//     qr_solver.compute(Hes);
//     const size_t rank_hes = qr_solver.rank();
//     // qr_solver.compute(Jac);
//     // const size_t rank_jac = qr_solver.rank();    

//     cout << "# [OPT BASIS] system size and rank: " << Hes.cols() << " " << rank_hes << endl;
//     // cout << "# jac rank: " << rank_jac << endl;
//     // cout << "# hes rank: " << rank_hes << endl;
//     // getchar();
//   } while (0);
// #endif

//   //  shared_ptr<Constraint<double>> null_con;

//   SmartPtr<ipopt_opt_framework> opt_prb = new ipopt_opt_framework(opt_obj, opt_con, &basis[0]);
//   SmartPtr<IpoptApplication> app = IpoptApplicationFactory();
//   ApplicationReturnStatus status = app->Initialize();
//   ASSERT(status == Ipopt::Solve_Succeeded);
//   status = app->OptimizeTNLP(opt_prb);
//   ASSERT(status == Ipopt::Solve_Succeeded);

//   // shared_ptr<hj::math_func::math_func_t<double, int32_t>> objE;
//   // objE = make_shared<func_transformer<double, int32_t>>(opt_obj, &basis[0]);
//   // boost::property_tree::ptree opts; {
//   //   opts.put("package.value", "jtf");
//   //   opts.put("alg.value", "SQP");
//   //   opts.put("iter.value", 10);
//   //   opts.put("epsg.value", 1e-8);
//   //   opts.put("linear_solver/type.value", "PETsc");
//   // }
//   // jtf::optimize(*objE, basis, opts, nullptr, nullptr, NULL);

//   const size_t verts_num = this->nods_.size(2), basis_num = this->basis_pnt_.size();
//   basisXY = itr_matrix<const double *>(basis_num, verts_num, &basis[0]);
}

void fine_quad_patch::optimize_basis_modi_houman(matd_t &basisXY) {
  ASSERT(0);
//   //-> kronecker delta constraint
//   shared_ptr<Constraint<double>> lag;
//   lag = make_shared<full_kronecker_delta_cons>(*this);

//   //-> OBJECTIVE
//   vector<shared_ptr<Functional<double>>> ebf(2); {
//     ebf[0] = make_shared<elemwise_polyharm_energy>(*this, 1e0);
//     ebf[1] = make_shared<least_square_wrapper<double>>(lag, 1e4);
//   }
//   shared_ptr<Functional<double>> opt_obj;
//   try {
//     opt_obj = make_shared<energy_t<double>>(ebf);
//   } catch ( ... ) {
//     ASSERT(0);
//   }

//   //-> init basis
//   matd_t basis = ones<double>(opt_obj->Nx(), 1);

// #if 0
//   do {
//     g_func = opt_obj;
//     int rtn = numeric_grad_check(eval_val_gra, 1, g_func->Nx(), &basis[0]);
//     cout << "# [optimize basis] numerical check: " << rtn << endl;
//     getchar();
//   } while (0);
// #endif

//   shared_ptr<Constraint<double>> opt_con;
//   SmartPtr<ipopt_opt_framework> opt_prb = new ipopt_opt_framework(opt_obj, opt_con, &basis[0]);
//   SmartPtr<IpoptApplication> app = IpoptApplicationFactory();
//   ApplicationReturnStatus status = app->Initialize();
//   ASSERT(status == Ipopt::Solve_Succeeded);
//   status = app->OptimizeTNLP(opt_prb);
//   ASSERT(status == Ipopt::Solve_Succeeded);

//   // shared_ptr<hj::math_func::math_func_t<double, int32_t>> objE;
//   // objE = make_shared<func_transformer<double, int32_t>>(opt_obj, &basis[0]);
//   // boost::property_tree::ptree opts; {
//   //   opts.put("package.value", "jtf");
//   //   opts.put("alg.value", "SQP");
//   //   opts.put("iter.value", 10);
//   //   opts.put("epsg.value", 1e-8);
//   //   opts.put("linear_solver/type.value", "PETsc");
//   // }
//   // jtf::optimize(*objE, basis, opts, nullptr, nullptr, NULL);

//   const size_t verts_num = this->nods_.size(2), basis_num = this->basis_pnt_.size();
//   basisXY = itr_matrix<const double *>(basis_num, verts_num, &basis[0]);
}

int fine_quad_patch::write_optimized_basis(const char *filename, const matd_t &basisXY) {
  matd_t reorder_basis = zeros<double>(basisXY.size(1), basisXY.size(2));
  reorder_basis(colon(), vert_strip_(colon())) = basisXY;
  return quad_mesh_write_to_vtk(filename, nods_, cell_, &reorder_basis, "POINT");
}

void fine_quad_patch::sample_FEM_basis(matd_t &basisFEM) {
  typedef void (*FEM_shape_func)(double *val, const double *eps);

  FEM_shape_func sf;
  const size_t num_basis_pt_axis = sqrt(basis_pnt_.size());
  if ( num_basis_pt_axis == 2 )
    sf = quad4_shape_func_val_;
  else if ( num_basis_pt_axis == 3 )
    sf = quad9_shape_func_val_;

  const size_t num_sample_pt_axis = sqrt(nods_.size(2));
  matd_t sample_pt_axis = 2.0/(num_sample_pt_axis-1)*matd_t(colon(0, num_sample_pt_axis-1));
  sample_pt_axis += -1;

  basisFEM.resize(basis_pnt_.size(), nods_.size(2));
  for (size_t i = 0; i < num_sample_pt_axis; ++i) {
    for (size_t j = 0; j < num_sample_pt_axis; ++j) {
      const size_t idx = i*num_sample_pt_axis+j;
      const double pts[2] = {sample_pt_axis[i], sample_pt_axis[j]};
      sf(&basisFEM(0, idx), pts);
    }
  }
}

void fine_quad_patch::apply_local_LMA(const size_t eignum, const quad_stencil *ptr_sten,
                                      SparseMatrix<double> &K_h, SparseMatrix<double> &M_h,
                                      SparseMatrix<double> &P, VectorXd &eigs_h) {
  const size_t rd = nods_.size(1);
  const size_t basis_num = basis_pnt_.size();
  const size_t verts_num = vert_strip_.size();

  const matd_t &basis = ptr_sten->basisXY_;
  ASSERT(basis.size(1) == basis_num);
  const bool isotropic = ptr_sten->isotropic_;
  const bool augmented = ptr_sten->augmented_;

  //-> assign the prolongation matrix
  get_prolongation_matrix(ptr_sten, this, P);

  //-> calc local stiffness matrix
  ASSERT(energy_.get());
  K_h.resize(energy_->Nx(), energy_->Nx()); {
    vector<Triplet<double>> trips;
    energy_->Hes(&REST_[0], &trips);
    K_h.setFromTriplets(trips.begin(), trips.end());
  }
  
  //-> calc local mass matrix
  calc_surf_mass_matrix(cell_, REST_, 1.0, &M_h);

  //-> solve gev on local fine mesh
  const int smallest = 0;
  solve_gen_eig_prob(K_h, M_h, eignum, smallest, "arpaca", &eigs_h, NULL);
}

int stencil_prolongate_patch(const matd_t &nods_H, const quad_stencil* ptr_sten,
                             fine_mesh_patch *ptr_patch,
                             mati_t &cell_h, matd_t &nods_h) {
  const size_t rd = ptr_patch->nods_.size(1);
  const size_t basis_num = ptr_patch->basis_pnt_.size();
  const size_t verts_num = ptr_patch->vert_strip_.size();
  ASSERT(nods_H.size(2) == basis_num);

  const matd_t &basis = ptr_sten->basisXY_;
  const bool isotropic = ptr_sten->isotropic_;
  const bool augmented = ptr_sten->augmented_;

  if ( !ptr_patch->Pro_.get() ) {
    ptr_patch->Pro_ = make_shared<SparseMatrix<double>>();
    get_prolongation_matrix(ptr_sten, ptr_patch, *ptr_patch->Pro_);
  }

  //-> estimate rotation
  matd_t R; {
    matd_t H = zeros<double>(ptr_sten->adjc_vert_num_, 2);
    const double xi[2] = {0};
    ptr_sten->rest_sf_jac_(&H[0], xi);
    matd_t Ds = nods_H*H, Dm = ptr_sten->adjc_rest_*H;
    inv(Dm);
    R = Ds*Dm;
    hj::polar2d rs;
    rs(R);
  }

  const matd_t &adjc_rest = ptr_sten->adjc_rest_;
  matd_t u0 = trans(R)*nods_H-adjc_rest;
  matd_t Pu(ptr_patch->REST_.size(1), ptr_patch->REST_.size(2));
  Eigen::Map<VectorXd>(&Pu[0], Pu.size()) = (*ptr_patch->Pro_)*Eigen::Map<const VectorXd>(&u0[0], u0.size());
  nods_h = R*(ptr_patch->REST_+Pu);
  cell_h = ptr_patch->cell_;

  return 0;
}

int get_prolongation_matrix(const quad_stencil *ptr_sten, const fine_mesh_patch *ptr_patch,
                            SparseMatrix<double> &P) {
  ASSERT(ptr_sten && ptr_patch);
  
  const size_t rd = ptr_patch->nods_.size(1);
  const size_t basis_num = ptr_patch->basis_pnt_.size();
  const size_t verts_num = ptr_patch->vert_strip_.size();

  const matd_t &basis  = ptr_sten->basisXY_;
  const bool isotropic = ptr_sten->isotropic_;
  const bool augmented = ptr_sten->augmented_;

  //-> get prolongation matrix
  P.resize(rd*verts_num, rd*basis_num); {
    vector<Triplet<double>> trips;
    for (size_t i = 0; i < ptr_patch->vert_strip_.size(); ++i) {
      const size_t pid = ptr_patch->vert_strip_[i];    
      for (size_t j = 0; j < ptr_patch->basis_pnt_.size(); ++j) {
        for (size_t k = 0; k < rd; ++k) { // each dimension
          if ( !isotropic && !augmented ) {
            trips.push_back(Triplet<double>(rd*pid+k, rd*j+k, basis(j, i+k*verts_num)));
          } else if ( !isotropic && augmented ) {
            for (size_t m = 0; m < rd; ++m)
              trips.push_back(Triplet<double>(rd*pid+k, rd*j+m, basis(j, i+(rd*k+m)*verts_num)));
          } else if ( isotropic ) {
            trips.push_back(Triplet<double>(rd*pid+k, rd*j+k, basis(j, i)));
          }
        }
      }
    }
    P.reserve(trips.size());
    P.setFromTriplets(trips.begin(), trips.end());
  }
  return 0;
}

int find_shared_boundary_element(const fine_quad_patch *patch_x,
                                 const fine_quad_patch *patch_y,
                                 mati_t &shared_bnd_edges,
                                 mati_t &corres_faces) {
  ASSERT(patch_x && patch_y);
  
  mati_t bnd_edge_x, bnd_edge_y;
  jtf::mesh::get_boundary_edge(*(patch_x->e2c_), bnd_edge_x);
  jtf::mesh::get_boundary_edge(*(patch_y->e2c_), bnd_edge_y);

  std::set<pair<size_t, size_t>> set_x, set_y;
  for (size_t i = 0; i < bnd_edge_x.size(2); ++i) {
    const size_t max_id = patch_x->new_to_orig_[max(bnd_edge_x(colon(), i))];
    const size_t min_id = patch_x->new_to_orig_[min(bnd_edge_x(colon(), i))];
    set_x.insert(std::make_pair(min_id, max_id));
  }
  for (size_t i = 0; i < bnd_edge_y.size(2); ++i) {
    const size_t max_id = patch_y->new_to_orig_[max(bnd_edge_y(colon(), i))];
    const size_t min_id = patch_y->new_to_orig_[min(bnd_edge_y(colon(), i))];
    set_y.insert(std::make_pair(min_id, max_id));
  }

  vector<pair<size_t, size_t>> bnd_edges(set_x.size()+set_y.size());
  auto it = std::set_intersection(set_x.begin(), set_x.end(), set_y.begin(), set_y.end(), bnd_edges.begin());
  bnd_edges.resize(it-bnd_edges.begin());

  if ( bnd_edges.size() == 0 )
    return __LINE__;

  shared_bnd_edges.resize(2, bnd_edges.size());
  for (size_t i = 0; i < bnd_edges.size(); ++i) {
    shared_bnd_edges(0, i) = bnd_edges[i].first;
    shared_bnd_edges(1, i) = bnd_edges[i].second;
  }           

  corres_faces.resize(2, bnd_edges.size());
  size_t i = 0;
  for (const auto &e : bnd_edges) {
    size_t left, right;
    {
      const size_t p = patch_x->global_to_local_[e.first];
      const size_t q = patch_x->global_to_local_[e.second];
      const auto facet = patch_x->e2c_->query(p, q);
      left = ( facet.first == -1 ) ? facet.second : facet.first;
    }
    {
      const size_t p = patch_y->global_to_local_[e.first];
      const size_t q = patch_y->global_to_local_[e.second];
      const auto facet = patch_y->e2c_->query(p, q);
      right = ( facet.first == -1 ) ? facet.second : facet.first;
    }
    corres_faces(0, i) = left;
    corres_faces(1, i) = right;
    ++i;
  }
  
  return 0;
}
///===== nonconforming energy =====///
stencil_nonconforming_d2::stencil_nonconforming_d2(const mati_t &quad_H, const matd_t &rest_H,
                                                   const vector<shared_ptr<quad_stencil>> &stencils,
                                                   const vector<shared_ptr<fine_mesh_patch>> &patches,
                                                   const double w)
    : w_(w), dim_(rest_H.size()), rest_H_(rest_H), stencils_(stencils), patches_(patches) {
  const size_t rd = rest_H.size(1);

  for (size_t i = 0; i < patches.size(); ++i) {
    if ( !patches[i]->Pro_ ) {
      patches[i]->Pro_ = make_shared<SparseMatrix<double>>();
      get_prolongation_matrix(stencils[i].get(), patches[i].get(), *patches[i]->Pro_);
    }
  }

  const double len = norm(rest_H(colon(), quad_H(0, 0))-rest_H(colon(), quad_H(1, 0)));

  for (size_t i = 0; i < stencils_.size(); ++i) {
    for (size_t j = i+1; j < stencils_.size(); ++j) {
      const auto &stenI = stencils_[i];
      const auto &stenJ = stencils_[j];
      const auto &ptchI = patches_[i];
      const auto &ptchJ = patches_[j];

      //-> check if stenicl i, j are adjacent
      mati_t shared_bs_pnt;
      if ( !is_adjc_quad_stencil(stenI, stenJ, shared_bs_pnt) )
        continue;

      //-> line integral length
      elem_bnd_len_.push_back((shared_bs_pnt.size()-1)*len);

      //-> stencil i and j is adjacent
      adjc_sten_.push_back(make_pair(i, j));
      const mati_t ptch_node_i = ptchI->new_to_orig_(ptchI->vert_strip_);
      const mati_t ptch_node_j = ptchJ->new_to_orig_(ptchJ->vert_strip_);
      mati_t shared_nd_pnt;
      bigbang::find_intersections(ptch_node_i, ptch_node_j, shared_nd_pnt);

      cout << "shared nd pnt: " << trans(shared_nd_pnt) << endl;

      mati_t corres_faces, bnd_edges;
      find_shared_boundary_element(dynamic_cast<fine_quad_patch*>(ptchI.get()),
                                   dynamic_cast<fine_quad_patch*>(ptchJ.get()),
                                   bnd_edges, corres_faces);

      shared_ptr<VectorXd> sk = make_shared<VectorXd>();
      *sk = VectorXd::Zero(2*shared_nd_pnt.size());
      for (size_t f = 0; f < corres_faces.size(2); ++f) {
        const double mu1 = ptchI->lame_(0, corres_faces(0, f));
        const double mu2 = ptchJ->lame_(0, corres_faces(1, f));
        const double homo_mu = 1.0/(1.0/mu1+1.0/mu2);

        const size_t ei = bnd_edges(0, f);
        const size_t ej = bnd_edges(1, f);

        // cout << "edge num: " << f << endl;
        // cout << "ei ej: " << ei << ", " << ej << endl;
        // cout << "mu1 mu2: " << mu1 << ", " << mu2 << endl;

        int pos_ei = std::find(shared_nd_pnt.begin(), shared_nd_pnt.end(), ei)-shared_nd_pnt.begin();
        int pos_ej = std::find(shared_nd_pnt.begin(), shared_nd_pnt.end(), ej)-shared_nd_pnt.begin();

        // cout << "pos ei ej: " << pos_ei << ", " << pos_ej << endl << endl;
        // getchar();

        sk->segment<2>(2*pos_ei) += homo_mu*Vector2d::Ones();
        sk->segment<2>(2*pos_ej) += homo_mu*Vector2d::Ones();
      }
      sk->cwiseSqrt();
      K_.push_back(sk);
      
      #if 0
      cout << "share nodes of patch " << i << ", " << j << endl;
      cout << trans(shared_nd_pnt) << endl;
      #endif

      const mati_t fixedI = ptchI->global_to_local_(shared_nd_pnt);
      const mati_t fixedJ = ptchJ->global_to_local_(shared_nd_pnt);
      ASSERT(fixedI.size() == fixedJ.size());

      #if 0
      cout << "fixedI " << trans(ptchI->vert_strip_perm_(fixedI)) << endl;
      cout << "fixedJ " << trans(ptchJ->vert_strip_perm_(fixedJ)) << endl;
      getchar();
      #endif

      shared_ptr<SparseMatrix<double>> Si = make_shared<SparseMatrix<double>>();
      shared_ptr<SparseMatrix<double>> Sj = make_shared<SparseMatrix<double>>();
      ASSERT(build_selection_matrix(fixedI, ptchI->nods_.size(2), rd, *Si) == 0);
      ASSERT(build_selection_matrix(fixedJ, ptchJ->nods_.size(2), rd, *Sj) == 0);
      S_.push_back(make_pair(Si, Sj));
    }
  }
}
size_t stencil_nonconforming_d2::Nx() const {
  return dim_;
}
int stencil_nonconforming_d2::Val(const double *x, double *val) const {
  itr_matrix<const double *> X(2, Nx()/2, x);

  for (size_t p = 0; p < adjc_sten_.size(); ++p) {
    const size_t I = adjc_sten_[p].first;
    const size_t J = adjc_sten_[p].second;    
    const auto &SI = *S_[p].first;
    const auto &SJ = *S_[p].second;
    
    matd_t xi_h; {
      const auto &sten  = stencils_[I];
      const auto &patch = patches_[I];
      const auto &P = *patch->Pro_;
      
      const matd_t x_H = X(colon(), sten->adjc_elem_);
      matd_t R = zeros<double>(2, 2);
      sten->query_local_frame(x_H, 0, R);

      const matd_t u_H = trans(R)*x_H-rest_H_(colon(), sten->adjc_elem_);
      matd_t u_h(patch->REST_.size(1), patch->REST_.size(2));
      Eigen::Map<VectorXd>(&u_h[0], u_h.size()) = P*Eigen::Map<const VectorXd>(&u_H[0], u_H.size());
      xi_h = R*(patch->REST_+u_h);      
    }
    
    matd_t xj_h; {
      const auto &sten  = stencils_[J];
      const auto &patch = patches_[J];
      const auto &P = *patch->Pro_;
      
      const matd_t x_H = X(colon(), sten->adjc_elem_);
      matd_t R = zeros<double>(2, 2);
      sten->query_local_frame(x_H, 0, R);

      const matd_t u_H = trans(R)*x_H-rest_H_(colon(), sten->adjc_elem_);
      matd_t u_h(patch->REST_.size(1), patch->REST_.size(2));
      Eigen::Map<VectorXd>(&u_h[0], u_h.size()) = P*Eigen::Map<const VectorXd>(&u_H[0], u_H.size());
      xj_h = R*(patch->REST_+u_h);
    }

    const SparseMatrix<double> K(K_[p]->asDiagonal());
    *val += 0.5*w_*elem_bnd_len_[p]*(K*SI*Eigen::Map<const VectorXd>(&xi_h[0], xi_h.size())-
                                     K*SJ*Eigen::Map<const VectorXd>(&xj_h[0], xj_h.size())).squaredNorm();
  }
  
  return 0;
}
int stencil_nonconforming_d2::Gra(const double *x, double *gra) const {
  itr_matrix<const double *> X(2, Nx()/2, x);
  itr_matrix<double *> G(2, Nx()/2, gra);

  for (size_t p = 0; p < adjc_sten_.size(); ++p) {
    const size_t I = adjc_sten_[p].first;
    const size_t J = adjc_sten_[p].second;
    const auto &Si = *S_[p].first;
    const auto &Sj = *S_[p].second;
    
    matd_t xi_h, Ri; {
      const auto &sten = stencils_[I];
      const auto &patch = patches_[I];
      const auto &P = *patch->Pro_;
      
      const matd_t x_H = X(colon(), sten->adjc_elem_);
      Ri = zeros<double>(2, 2);
      sten->query_local_frame(x_H, 0, Ri);

      const matd_t u_H = trans(Ri)*x_H-rest_H_(colon(), sten->adjc_elem_);
      matd_t u_h(patch->REST_.size(1), patch->REST_.size(2));
      Eigen::Map<VectorXd>(&u_h[0], u_h.size()) = P*Eigen::Map<const VectorXd>(&u_H[0], u_H.size());
      xi_h = Ri*(patch->REST_+u_h);      
    }
    
    matd_t xj_h, Rj; {
      const auto &sten = stencils_[J];
      const auto &patch = patches_[J];
      const auto &P = *patch->Pro_;
      
      const matd_t x_H = X(colon(), sten->adjc_elem_);
      Rj = zeros<double>(2, 2);
      sten->query_local_frame(x_H, 0, Rj);

      const matd_t u_H = trans(Rj)*x_H-rest_H_(colon(), sten->adjc_elem_);
      matd_t u_h(patch->REST_.size(1), patch->REST_.size(2));
      Eigen::Map<VectorXd>(&u_h[0], u_h.size()) = P*Eigen::Map<const VectorXd>(&u_H[0], u_H.size());
      xj_h = Rj*(patch->REST_+u_h);
    }

    const SparseMatrix<double> K(K_[p]->asDiagonal());
    const VectorXd resd = K*(Si*Eigen::Map<const VectorXd>(&xi_h[0], xi_h.size())
                            -Sj*Eigen::Map<const VectorXd>(&xj_h[0], xj_h.size()));
    
    VectorXd gI, gJ;
    {
      const size_t sten_num = stencils_[I]->adjc_elem_.size();
      const size_t ptch_num = patches_[I]->nods_.size(2);
      SparseMatrix<double> RH, Rh;
      block_diagonalize_matrix(Ri, sten_num, RH);
      block_diagonalize_matrix(Ri, ptch_num, Rh);
      const auto jacI = K*Si*Rh*(*patches_[I]->Pro_)*RH.transpose();
      gI = w_*elem_bnd_len_[p]*jacI.transpose()*resd;
    }
    {
      const size_t sten_num = stencils_[J]->adjc_elem_.size();
      const size_t ptch_num = patches_[J]->nods_.size(2);
      SparseMatrix<double> RH, Rh;
      block_diagonalize_matrix(Rj, sten_num, RH);
      block_diagonalize_matrix(Rj, ptch_num, Rh);
      const auto jacJ = K*Sj*Rh*(*patches_[J]->Pro_)*RH.transpose();
      gJ = -w_*elem_bnd_len_[p]*jacJ.transpose()*resd;
    }
    
    G(colon(), stencils_[I]->adjc_elem_) += itr_matrix<const double *>(2, gI.size()/2, gI.data());
    G(colon(), stencils_[J]->adjc_elem_) += itr_matrix<const double *>(2, gJ.size()/2, gJ.data());
  }
  
  return 0;
}
int stencil_nonconforming_d2::Hes(const double *x, vector<Triplet<double>> *hes) const {
  itr_matrix<const double *> X(2, Nx()/2, x);

  for (size_t p = 0; p < adjc_sten_.size(); ++p) {
    const size_t I = adjc_sten_[p].first;
    const size_t J = adjc_sten_[p].second;    
    const auto &Si = *S_[p].first;
    const auto &Sj = *S_[p].second;
    
    matd_t xi_h, Ri; {
      const auto &sten = stencils_[I];
      const auto &patch = patches_[I];
      const auto &P = *patch->Pro_;
      
      const matd_t x_H = X(colon(), sten->adjc_elem_);
      Ri = zeros<double>(2, 2);
      sten->query_local_frame(x_H, 0, Ri);
    }
    
    matd_t xj_h, Rj; {
      const auto &sten = stencils_[J];
      const auto &patch = patches_[J];
      const auto &P = *patch->Pro_;
      
      const matd_t x_H = X(colon(), sten->adjc_elem_);
      Rj = zeros<double>(2, 2);
      sten->query_local_frame(x_H, 0, Rj);
    }

    const SparseMatrix<double> K(K_[p]->asDiagonal());    
    SparseMatrix<double> jacI, jacJ;
    {
      const size_t sten_num = stencils_[I]->adjc_elem_.size();
      const size_t ptch_num = patches_[I]->nods_.size(2);
      SparseMatrix<double> RH, Rh;
      block_diagonalize_matrix(Ri, sten_num, RH);
      block_diagonalize_matrix(Ri, ptch_num, Rh);
      jacI = K*Si*Rh*(*patches_[I]->Pro_)*RH.transpose();
    }
    {
      const size_t sten_num = stencils_[J]->adjc_elem_.size();
      const size_t ptch_num = patches_[J]->nods_.size(2);
      SparseMatrix<double> RH, Rh;
      block_diagonalize_matrix(Rj, sten_num, RH);
      block_diagonalize_matrix(Rj, ptch_num, Rh);
      jacJ = -K*Sj*Rh*(*patches_[J]->Pro_)*RH.transpose();
    }

    const size_t idx[2] = {I, J};
    const vector<SparseMatrix<double>*> Jac{&jacI, &jacJ};

    for (size_t m = 0; m < 2; ++m) {
      for (size_t n = 0; n < 2; ++n) {
        const SparseMatrix<double> H = w_*elem_bnd_len_[p]*Jac[m]->transpose()*(*Jac[n]);

        for (size_t k = 0; k < H.outerSize(); ++k) {
          for (SparseMatrix<double>::InnerIterator it(H, k); it; ++it) {
            const size_t row = it.row(), col = it.col();
            const size_t II = 2*stencils_[idx[m]]->adjc_elem_[row/2]+row%2;
            const size_t JJ = 2*stencils_[idx[n]]->adjc_elem_[col/2]+col%2;
            hes->push_back(Triplet<double>(II, JJ, it.value()));
          }
        }
      }
    }    
  }
  
  return 0;
}

///==================== fine hexs patch ====================///
fine_hexs_patch::fine_hexs_patch(const hexs_stencil &dom, const mati_t &hexs_H,
                                 const mati_t &hexs_h, const matd_t &nods_h, const matd_t &lame_h,
                                 const bool to_be_refine, const ptree &pt)
    : fine_mesh_patch(pt) {
  ONE_TO_MANY_ = hexs_h.size(2)/hexs_H.size(2);
  ASSERT(hexs_h.size(2)%hexs_H.size(2) == 0);

  //-> construct mesh patch for stencil $dom
  cell_.resize(hexs_h.size(1), dom.faces_.size()*ONE_TO_MANY_);
  lame_.resize(lame_h.size(1), dom.faces_.size()*ONE_TO_MANY_);
  for (size_t i = 0; i < dom.faces_.size(); ++i) {
    const auto fid = dom.faces_[i];
    cell_(colon(), colon(i*ONE_TO_MANY_, (i+1)*ONE_TO_MANY_-1)) = hexs_h(colon(), colon(fid*ONE_TO_MANY_, (fid+1)*ONE_TO_MANY_-1));
    lame_(colon(), colon(i*ONE_TO_MANY_, (i+1)*ONE_TO_MANY_-1)) = lame_h(colon(), colon(fid*ONE_TO_MANY_, (fid+1)*ONE_TO_MANY_-1));
  }
  nods_ = nods_h;
  remove_extra_verts(cell_, nods_, &new_to_orig_);

  //-> get the vertices also on the coarse mesh
  global_to_local_ = ones<size_t>(nods_h.size(2), 1)*-1;
  global_to_local_(new_to_orig_) = colon(0, new_to_orig_.size()-1);
  basis_pnt_ = global_to_local_(dom.adjc_elem_);

  //-> subdvide patch for smoother basis if neccessary
  if ( to_be_refine ) {
    ;
  }

  //-> record the rest shape
  REST_ = nods_;

  //-> init face2hex
  f2h_.reset(jtf::mesh::face2hex_adjacent::create(cell_));
  
  //-> extract verts strip according to basis order
  strip_ext_ = make_shared<straight_cube_strip_extractor>(cell_);
  this->extract_vert_strip();
  
  //-> calc vert strip permutation and basis points permutation
  vert_strip_perm_ = -1*ones<size_t>(nods_.size(2), 1);
  vert_strip_perm_(vert_strip_(colon())) = colon(0, vert_strip_.size()-1);
  basis_pnt_perm_ = -1*ones<size_t>(nods_.size(2), 1);
  basis_pnt_perm_(basis_pnt_(colon())) = colon(0, basis_pnt_.size()-1);

  //-> init elastic energy
  energy_ = build_elastic_energy<voxel_elastic_potential>(pt.get<string>("elas.value"), cell_, REST_, lame_); {
    vector<Triplet<double>> trips;
    energy_->Hes(&REST_[0], &trips);
    K_.resize(energy_->Nx(), energy_->Nx());
    K_.setFromTriplets(trips.begin(), trips.end());
  }
  
#if 0
  static size_t count = 0;
  char outfile[256];
  sprintf(outfile, "./tmp_output/patch_mesh_%zu.vtk", count);
  hex_mesh_write_to_vtk(outfile, nods_, cell_);
  cout << trans(vert_strip_) << endl;
  getchar();
  ++count;
#endif

  //-> init FEM basis indicator
  ASSERT(basis_pnt_.size() == 8);
  this->sample_FEM_basis(this->FEM_basis_);
  FEM_indicator_.resize(FEM_basis_.size(1), FEM_basis_.size(2));
  for (size_t i = 0; i < FEM_indicator_.size(); ++i)
    FEM_indicator_[i] = (fabs(FEM_basis_[i]) < 1e-8) ? 0 : 1;

#if 0
  const string outf("./tmp_output/indicator.vtk");
  const matd_t float_ind = FEM_indicator_;
  this->write_optimized_basis(outf.c_str(), float_ind);
  getchar();
#endif
}

void fine_hexs_patch::extract_vert_strip() {
  const size_t START_FACE = 0;

  const size_t vert_dim = cbrt(nods_.size(2));
  ASSERT(vert_dim*vert_dim*vert_dim == nods_.size(2));
  size_t vert_arr[vert_dim][vert_dim][vert_dim] = {0};

  vector<size_t> x_faces;
  strip_ext_->extract(START_FACE, 'X', x_faces);

  for (size_t i = 0; i < x_faces.size(); ++i) {
    vector<size_t> y_faces;
    strip_ext_->extract(x_faces[i], 'Y', y_faces);

    for (size_t j = 0; j < y_faces.size(); ++j) {
      vector<size_t> z_faces;
      strip_ext_->extract(y_faces[j], 'Z', z_faces);

      for (size_t k = 0; k < z_faces.size(); ++k) {
        vert_arr[i][j][k] = cell_(0, z_faces[k]);
        if ( k+1 == z_faces.size() )                   // touch z
          vert_arr[i][j][k+1] = cell_(4, z_faces[k]);
        if ( j+1 == y_faces.size() )                   // touch y
          vert_arr[i][j+1][k] = cell_(2, z_faces[k]);
        if ( i+1 == x_faces.size() )                   // touch x
          vert_arr[i+1][j][k] = cell_(1, z_faces[k]);
        if ( k+1 == z_faces.size() && j+1 == y_faces.size() ) // touch yz
          vert_arr[i][j+1][k+1] = cell_(6, z_faces[k]);
        if ( i+1 == x_faces.size() && j+1 == y_faces.size() ) // touch xy
          vert_arr[i+1][j+1][k] = cell_(3, z_faces[k]);
        if ( i+1 == x_faces.size() && k+1 == z_faces.size() ) // touch xz
          vert_arr[i+1][j][k+1] = cell_(5, z_faces[k]);
        if ( i+1 == x_faces.size() && j+1 == y_faces.size() && k+1 == z_faces.size() ) // touch xyz
          vert_arr[i+1][j+1][k+1] = cell_(7, z_faces[k]);
      }
    }
  }

  vert_strip_.resize(nods_.size(2), 1);
  std::copy(&vert_arr[0][0][0], &vert_arr[0][0][0]+nods_.size(2), vert_strip_.begin());

  //-> build parametric nodes
  matd_t coord_pts = 2.0/(vert_dim-1)*matd_t(colon(0, vert_dim-1));
  coord_pts += -1;
  param_nods_ = zeros<double>(3, nods_.size(2));
  for (size_t i = 0; i < vert_dim; ++i) {
    for (size_t j = 0; j < vert_dim; ++j) {
      for (size_t k = 0; k < vert_dim; ++k) {
        const size_t pid = vert_arr[i][j][k];
        param_nods_(0, pid) = coord_pts[i];
        param_nods_(1, pid) = coord_pts[j];
        param_nods_(2, pid) = coord_pts[k];
      }
    }
  }
}

void fine_hexs_patch::modal_analysis() {
  // //-> stiffness matrix
  // ASSERT(energy_.get());
  // SparseMatrix<double> K(energy_->Nx(), energy_->Nx()); {
  //   vector<Triplet<double>> trips;
  //   energy_->Hes(&REST_[0], &trips);
  //   K.setFromTriplets(trips.begin(), trips.end());
  // }

  //-> mass matrix
  M_.resize(energy_->Nx(), energy_->Nx());  
  const bool lumped = true;
  calc_mass_matrix(cell_, REST_, 1.0, 3, &M_, lumped);
  
  // const size_t eignum = 20; //K.cols()*0.8;
  // cout << "[LMA] eigen number to be solved: " << eignum << endl;

  // const int smallest = 0;
  // VectorXd eigenvalues; MatrixXd eigenvectors;
  // int eig_rtn = solve_gen_eig_prob(K, M_, eignum, smallest, "arpaca", &eigenvalues, &eigenvectors);
  // ASSERT(eig_rtn == 0);
  // cout << "[HEXS LMA] norm of head 6 lambda: " << eigenvalues.head(6).norm() << endl;

  // bkp_freqs_ = eigenvalues;
  // bkp_modes_ = eigenvectors;

  bkp_freqs_ = VectorXd::Zero(6);
  bkp_modes_ = MatrixXd::Zero(nods_.size(), 6);

#if 0
  do {
    cout << "# freqs: " << bkp_freqs_.transpose() << endl;
    //-> write modal basis
    static size_t count = 0;
    char outfile[256];
    for (size_t i = 0; i < 10; ++i) {
      const matd_t disp = itr_matrix<const double*>(nods_.size(1), nods_.size(2), &bkp_modes_(0, i));
      const matd_t def_nods = nods_+disp;

      sprintf(outfile, "./tmp_output/patch%zu_origin_mode%03zu.vtk", count, i);
      hex_mesh_write_to_vtk(outfile, def_nods, cell_, &lame_, "CELL");
      sprintf(outfile, "./tmp_output/patch%zu_dispfield%03zu.vtk", count, i);
      hex_mesh_write_to_vtk(outfile, nods_, cell_, &disp, "POINT");
    }
    ++count;
    getchar();
  } while (0);
#endif
  
#if 0
  do {
    cout << "## U^TKU-Lambda: "
         << (bkp_modes_.transpose()*K*bkp_modes_-MatrixXd(bkp_freqs_.asDiagonal())).norm() << endl;
    cout << "## U^TMU-Id: "
         << (bkp_modes_.transpose()*M*bkp_modes_-MatrixXd::Identity(eignum, eignum)).norm() << endl;
  } while (0);
#endif
  
  //-> separate the translational and rotational parts
  const size_t mode_rows = bkp_modes_.rows();

  MatrixXd tran_b = MatrixXd::Zero(mode_rows, 3);
  for (size_t i = 0; i < mode_rows/3; ++i)
    tran_b.block(3*i, 0, 3, 3) = Matrix3d::Identity();
  tran_b.col(0) /= sqrt(tran_b.col(0).dot(M_*tran_b.col(0)));
  tran_b.col(1) /= sqrt(tran_b.col(1).dot(M_*tran_b.col(1)));
  tran_b.col(2) /= sqrt(tran_b.col(2).dot(M_*tran_b.col(2)));

  //-> construct rotation modes
  const matd_t Id = eye<double>(3);
  MatrixXd rot_b = MatrixXd::Zero(mode_rows, 3);
  const matd_t bc = nods_*ones<double>(nods_.size(2), 1)/nods_.size(2);
  for (size_t i = 0; i < 3; ++i) {
    matd_t disp = zeros<double>(nods_.size(1), nods_.size(2));
    for (size_t j = 0; j < nods_.size(2); ++j) {
      disp(colon(), j) = cross(Id(colon(), i), nods_(colon(), j)-bc);
    }
    rot_b.col(i) = Eigen::Map<const VectorXd>(&disp[0], disp.size());
  }

  // const matd_t bbc = nods_*ones<double>(nods_.size(2), 1)/nods_.size(2);
  // cout << "real: " << endl;
  // matd_t tmp = nods_-bbc*ones<double>(1, nods_.size(2));
  // cout << tmp/fabs(max(tmp)) << endl;
  // cout << "param: " << endl;
  // cout << param_nods_ << endl;
  // getchar();

  // rot_b.col(0) /= sqrt(rot_b.col(0).dot(M_*rot_b.col(0)));
  // rot_b.col(1) /= sqrt(rot_b.col(1).dot(M_*rot_b.col(1)));
  // rot_b.col(2) /= sqrt(rot_b.col(2).dot(M_*rot_b.col(2)));
  
  // const MatrixXd kernel_K = bkp_modes_.topLeftCorner(mode_rows, 6);
  // MatrixXd LHS = MatrixXd::Zero(3, 6); {
  //   LHS.row(0) = tran_b.col(0).transpose()*M_*kernel_K;
  //   LHS.row(1) = tran_b.col(1).transpose()*M_*kernel_K;
  //   LHS.row(2) = tran_b.col(2).transpose()*M_*kernel_K;
  // }

  // FullPivLU<MatrixXd> lu(LHS);
  // MatrixXd rot_b = kernel_K*lu.kernel(); {
  //   ASSERT(lu.kernel().cols() == 3);
  //   //-> gram schimdt process with metric M
  //   rot_b.col(0) /= sqrt(rot_b.col(0).dot(M_*rot_b.col(0)));

  //   rot_b.col(1) -= (rot_b.col(1).dot(M_*rot_b.col(0))/rot_b.col(0).dot(M_*rot_b.col(0)))*rot_b.col(0);
  //   rot_b.col(1) /= sqrt(rot_b.col(1).dot(M_*rot_b.col(1)));

  //   rot_b.col(2) -=
  //       (rot_b.col(2).dot(M_*rot_b.col(0))/rot_b.col(0).dot(M_*rot_b.col(0)))*rot_b.col(0)+
  //       (rot_b.col(2).dot(M_*rot_b.col(1))/rot_b.col(1).dot(M_*rot_b.col(1)))*rot_b.col(1);
  //   rot_b.col(2) /= sqrt(rot_b.col(2).dot(M_*rot_b.col(2)));
  // }

  bkp_modes_.col(0) = tran_b.col(0);
  bkp_modes_.col(1) = tran_b.col(1);
  bkp_modes_.col(2) = tran_b.col(2);
  bkp_modes_.col(3) = rot_b.col(0);
  bkp_modes_.col(4) = rot_b.col(1);
  bkp_modes_.col(5) = rot_b.col(2);

  // //-> only preserve deformable bases
  // const size_t non_rigid_mode_num = bkp_freqs_.size()-6;
  // freqs_ = bkp_freqs_.tail(non_rigid_mode_num);
  // modes_ = bkp_modes_.topRightCorner(mode_rows, non_rigid_mode_num);

#if 0
  do {
    const MatrixXd rU = bkp_modes_.topLeftCorner(bkp_modes_.rows(), 6);
    const MatrixXd UTKU = rU.transpose()*K*rU;
    const MatrixXd UTMU = rU.transpose()*M*rU;
    cout << "U^TKU norm: " << UTKU.norm() << endl;
    cout << "U^TMU-Id norm: " << (UTMU-MatrixXd::Identity(6, 6)).norm() << endl;
  } while (0);
#endif
  
#if 1
  do {
    char outfile[256];
    static size_t count = 0;
    if ( count == 0 ) {
      //-> write rotational modes
      sprintf(outfile, "%s/REST_patch%04zu.vtk", pt_.get<string>("outdir.value").c_str(), count);
      hex_mesh_write_to_vtk(outfile, REST_, cell_, &lame_, "CELL");
      sprintf(outfile, "%s/REST_domain_patch%04zu.vtk", pt_.get<string>("outdir.value").c_str(), count);
      hex_mesh_write_to_vtk(outfile, param_nods_, cell_);
      for (size_t i = 3; i < 6; ++i) {
        const matd_t disp = itr_matrix<const double *>(nods_.size(1), nods_.size(2), &bkp_modes_(0, i));
        const matd_t def_nods = nods_+0.0001*disp;
        sprintf(outfile, "%s/ROT_patch%04zu_mode_%03zu.vtk", pt_.get<string>("outdir.value").c_str(), count, i-3);
        hex_mesh_write_to_vtk(outfile, def_nods, cell_, &lame_, "CELL");
      }
    }
    ++count;
  } while (0);
#endif
}

void fine_hexs_patch::optimize_basis(matd_t &basisXYZ) {
  ASSERT(0);
//   const size_t mode_num = basis_pnt_.size()+pt_.get<int>("eigmode_plus.value");

//   const double power = pt_.get<double>("eig_power.value");
//   shared_ptr<Constraint<double>> ip_cons;
//   ip_cons = make_shared<full_patch_interp_cons>(*this, power, pt_);

//   //-> CONSTRAINTS
//   vector<shared_ptr<Constraint<double>>> cbf(2); {
//     cbf[0] = make_shared<full_partition_to_unity_cons>(*this);
//     cbf[1] = make_shared<full_kronecker_delta_cons>(*this);
//   }
//   shared_ptr<Constraint<double>> opt_con;
//   try {
//     opt_con = make_shared<constraint_t<double>>(cbf);
//   } catch ( ... ) {
//     ASSERT(0);
//   }

//   //-> OBJECTIVE
//   vector<shared_ptr<Functional<double>>> ebf(3); {
//     ebf[0] = make_shared<least_square_wrapper<double>>(ip_cons, 1e0);
//     // ebf[1] = make_shared<least_square_wrapper<double>>(cbf[0], 1e4);
//     // ebf[2] = make_shared<least_square_wrapper<double>>(cbf[1], 1e4);
//   }
//   shared_ptr<Functional<double>> opt_obj;
//   try {
//     opt_obj = make_shared<energy_t<double>>(ebf);
//   } catch ( ... ) {
//     ASSERT(0);
//   }

//   ASSERT(opt_obj->Nx() == opt_con->Nx());
//   matd_t basis = ones<double>(opt_obj->Nx(), 1);

// #if 0
//   do {
//     g_func = opt_obj;
//     int rtn = numeric_grad_check(eval_val_gra, 1, g_func->Nx(), &basis[0]);
//     cout << "# [optimize basis] numerical check: " << rtn << endl;
//     getchar();
//   } while (0);
// #endif

// #if 0
//   do {
//     //-> check the rank of the hessian
//     SparseMatrix<double> Hes(opt_obj->Nx(), opt_obj->Nx()); {
//       vector<Triplet<double>> trips;
//       opt_obj->Hes(&basis[0], &trips);
//       Hes.setFromTriplets(trips.begin(), trips.end());
//     }
//     // SparseMatrix<double> Jac(opt_con->Nf(), opt_con->Nx()); {
//     //   vector<Triplet<double>> trips;
//     //   opt_con->Jac(&basis[0], 0, &trips);
//     //   Jac.setFromTriplets(trips.begin(), trips.end());
//     // }
    
//     SparseQR<SparseMatrix<double>, COLAMDOrdering<int>> qr_solver;
//     qr_solver.compute(Hes);
//     const size_t rank_hes = qr_solver.rank();
//     // qr_solver.compute(Jac);
//     // const size_t rank_jac = qr_solver.rank();    

//     cout << "[OPT BASIS] system size and rank: " << Hes.cols() << " " << rank_hes << endl;
//     //    cout << "# jac rank: " << rank_jac << endl;
//     //    cout << "# hes rank: " << rank_hes << endl;
//     //    getchar();
//   } while (0);
// #endif

//   // shared_ptr<hj::math_func::math_func_t<double, int32_t>> objE;
//   // objE = make_shared<func_transformer<double, int32_t>>(opt_obj, &basis[0]);
//   // boost::property_tree::ptree opts; {
//   //   opts.put("package.value", "jtf");
//   //   opts.put("alg.value", "SQP");
//   //   opts.put("iter.value", 10);
//   //   opts.put("epsg.value", 1e-8);
//   //   opts.put("linear_solver/type.value", "PETsc");
//   // }
//   // jtf::optimize(*objE, basis, opts, nullptr, nullptr, NULL);

//   shared_ptr<Constraint<double>> null_con;

//   SmartPtr<ipopt_opt_framework> opt_prb = new ipopt_opt_framework(opt_obj, opt_con, &basis[0]);
//   SmartPtr<IpoptApplication> app = IpoptApplicationFactory();
//   ApplicationReturnStatus status = app->Initialize();
//   ASSERT(status == Ipopt::Solve_Succeeded);
//   status = app->OptimizeTNLP(opt_prb);
//   ASSERT(status == Ipopt::Solve_Succeeded);

//   const size_t verts_num = this->nods_.size(2), basis_num = this->basis_pnt_.size();
//   basisXYZ = itr_matrix<const double *>(basis_num, verts_num, &basis[0]);  
}

void fine_hexs_patch::optimize_basis_modi_houman(matd_t &basisXYZ) {
  ASSERT(0);
}

int fine_hexs_patch::write_optimized_basis(const char *filename, const matd_t &basisXYZ) {
  matd_t reorder_basis = zeros<double>(basisXYZ.size(1),  basisXYZ.size(2));
  reorder_basis(colon(), vert_strip_(colon())) = basisXYZ;
  return hex_mesh_write_to_vtk(filename, nods_, cell_, &reorder_basis, "POINT");
}

void fine_hexs_patch::sample_FEM_basis(matd_t &basisFEM) {
  typedef void (*FEM_shape_func)(double *val, const double *eps);

  FEM_shape_func sf;
  const size_t num_basis_pt_axis = cbrt(basis_pnt_.size());
  if ( num_basis_pt_axis == 2 )
    sf = hex8_shape_func_;
  else if ( num_basis_pt_axis == 3 )
    sf = hex27_shape_func_;

  const size_t num_sample_pt_axis = cbrt(nods_.size(2));
  matd_t sample_pt_axis = 2.0/(num_sample_pt_axis-1)*matd_t(colon(0, num_sample_pt_axis-1));
  sample_pt_axis += -1;

  basisFEM.resize(basis_pnt_.size(), nods_.size(2));
  for (size_t i = 0; i < num_sample_pt_axis; ++i) {
    for (size_t j = 0; j < num_sample_pt_axis; ++j) {
      for (size_t k = 0; k < num_sample_pt_axis; ++k) {
        const size_t idx = i*num_sample_pt_axis*num_sample_pt_axis+j*num_sample_pt_axis+k;
        const double pts[3] = {sample_pt_axis[i], sample_pt_axis[j], sample_pt_axis[k]};
        sf(&basisFEM(0, idx), pts);
      }
    }
  }
}

void fine_hexs_patch::apply_local_LMA(const size_t eignum, const hexs_stencil *ptr_sten,
                                      SparseMatrix<double> &K_h, SparseMatrix<double> &M_h,
                                      SparseMatrix<double> &P, VectorXd &eigs_h) {
  const size_t rd = nods_.size(1);
  const size_t basis_num = basis_pnt_.size();
  const size_t verts_num = vert_strip_.size();

  const matd_t &basis = ptr_sten->basisXYZ_;
  ASSERT(basis.size(1) == basis_num);
  const bool isotropic = ptr_sten->isotropic_;
  const bool augmented = ptr_sten->augmented_;

  //-> assign the prolongation matrix
  get_prolongation_matrix(ptr_sten, this, P);

  //-> calc local stiffness matrix
  ASSERT(energy_.get());
  K_h.resize(energy_->Nx(), energy_->Nx()); {
    vector<Triplet<double>> trips;
    energy_->Hes(&REST_[0], &trips);
    K_h.setFromTriplets(trips.begin(), trips.end());
  }
  
  //-> calc local mass matrix
  const bool lumped = true;
  calc_mass_matrix(cell_, REST_, 1.0, REST_.size(1), &M_h, lumped);

  //-> solve gev on local fine mesh
  const int smallest = 0;
  solve_gen_eig_prob(K_h, M_h, eignum, smallest, "arpaca", &eigs_h, NULL);
}

int get_prolongation_matrix(const hexs_stencil *ptr_sten, const fine_mesh_patch *ptr_patch,
                            SparseMatrix<double> &P) {
  ASSERT(ptr_sten && ptr_patch);
  
  const size_t rd = ptr_patch->nods_.size(1);
  const size_t basis_num = ptr_patch->basis_pnt_.size();
  const size_t verts_num = ptr_patch->vert_strip_.size();

  const matd_t &basis  = ptr_sten->basisXYZ_;
  const bool isotropic = ptr_sten->isotropic_;
  const bool augmented = ptr_sten->augmented_;

  //-> get prolongation matrix
  P.resize(rd*verts_num, rd*basis_num); {
    vector<Triplet<double>> trips;
    for (size_t i = 0; i < ptr_patch->vert_strip_.size(); ++i) {
      const size_t pid = ptr_patch->vert_strip_[i];    
      for (size_t j = 0; j < ptr_patch->basis_pnt_.size(); ++j) {
        for (size_t k = 0; k < rd; ++k) { // each dimension
          if ( !isotropic && !augmented ) {
            trips.push_back(Triplet<double>(rd*pid+k, rd*j+k, basis(j, i+k*verts_num)));
          } else if ( !isotropic && augmented ) {
            for (size_t m = 0; m < rd; ++m)
              trips.push_back(Triplet<double>(rd*pid+k, rd*j+m, basis(j, i+(rd*k+m)*verts_num)));
          } else if ( isotropic ) {
            trips.push_back(Triplet<double>(rd*pid+k, rd*j+k, basis(j, i)));
          }
        }
      }
    }
    P.reserve(trips.size());
    P.setFromTriplets(trips.begin(), trips.end());
  }
  return 0;
}

int stencil_prolongate_patch(const matd_t &nods_H, const hexs_stencil* ptr_sten,
                             fine_mesh_patch *ptr_patch,
                             mati_t &cell_h, matd_t &nods_h) {
  const size_t rd = ptr_patch->nods_.size(1);
  const size_t basis_num = ptr_patch->basis_pnt_.size();
  const size_t verts_num = ptr_patch->vert_strip_.size();
  ASSERT(nods_H.size(2) == basis_num);

  const matd_t &basis = ptr_sten->basisXYZ_;
  const bool isotropic = ptr_sten->isotropic_;
  const bool augmented = ptr_sten->augmented_;

  if ( !ptr_patch->Pro_.get() ) {
    ptr_patch->Pro_ = make_shared<SparseMatrix<double>>();
    get_prolongation_matrix(ptr_sten, ptr_patch, *ptr_patch->Pro_);
  }

  //-> estimate rotation
  matd_t R; {
    matd_t H = zeros<double>(ptr_sten->adjc_vert_num_, 3);
    const double xi[3] = {0};
    ptr_sten->rest_sf_jac_(&H[0], xi);
    matd_t Ds = nods_H*H, Dm = ptr_sten->adjc_rest_*H;
    inv(Dm);
    R = Ds*Dm;
    hj::polar3d rs;
    rs(R);
  }

  const matd_t &adjc_rest = ptr_sten->adjc_rest_;
  matd_t u0 = trans(R)*nods_H-adjc_rest;
  matd_t Pu(ptr_patch->REST_.size(1), ptr_patch->REST_.size(2));
  Eigen::Map<VectorXd>(&Pu[0], Pu.size()) = (*ptr_patch->Pro_)*Eigen::Map<const VectorXd>(&u0[0], u0.size());
  nods_h = R*(ptr_patch->REST_+Pu);
  cell_h = ptr_patch->cell_;

  return 0;
}

int find_shared_boundary_element(const fine_hexs_patch *patch_x,
                                 const fine_hexs_patch *patch_y,
                                 mati_t &shared_bnd_faces,
                                 mati_t &corres_elems) {
  typedef std::tuple<size_t, size_t, size_t, size_t> tuple4_size_t;
  ASSERT(patch_x && patch_y);

  mati_t bnd_face_x, bnd_face_y;
  jtf::mesh::get_outside_face(*(patch_x->f2h_), bnd_face_x);
  jtf::mesh::get_outside_face(*(patch_y->f2h_), bnd_face_y);
  ASSERT(bnd_face_x.size(1) == 4 && bnd_face_y.size(1) == 4);
  bnd_face_x(colon()) = zjucad::matrix::temp(patch_x->new_to_orig_(bnd_face_x(colon())));
  bnd_face_y(colon()) = zjucad::matrix::temp(patch_y->new_to_orig_(bnd_face_y(colon())));

  std::set<tuple4_size_t> set_x, set_y;
  for (size_t i = 0; i < bnd_face_x.size(2); ++i) {
    std::sort(&bnd_face_x(0, i), &bnd_face_x(0, i)+4);
    set_x.insert(make_tuple(bnd_face_x(0, i), bnd_face_x(1, i), bnd_face_x(2, i), bnd_face_x(3, i)));
  }
  for (size_t i = 0; i < bnd_face_y.size(2); ++i) {
    std::sort(&bnd_face_y(0, i), &bnd_face_y(0, i)+4);
    set_y.insert(make_tuple(bnd_face_y(0, i), bnd_face_y(1, i), bnd_face_y(2, i), bnd_face_y(3, i)));
  }

  vector<tuple4_size_t> bnd_faces(set_x.size()+set_y.size());
  auto it = std::set_intersection(set_x.begin(), set_x.end(), set_y.begin(), set_y.end(), bnd_faces.begin());
  bnd_faces.resize(it-bnd_faces.begin());

  if ( bnd_faces.size() == 0 )
    return __LINE__;

  shared_bnd_faces.resize(4, bnd_faces.size());
  for (size_t i = 0; i < bnd_faces.size(); ++i) {
    shared_bnd_faces(0, i) = std::get<0>(bnd_faces[i]);
    shared_bnd_faces(1, i) = std::get<1>(bnd_faces[i]);
    shared_bnd_faces(2, i) = std::get<2>(bnd_faces[i]);
    shared_bnd_faces(3, i) = std::get<3>(bnd_faces[i]);
  }

  corres_elems.resize(2, bnd_faces.size());
  size_t i = 0;
  for (const auto &e : bnd_faces) {
    size_t left, right;
    {
      mati_t corner = patch_x->global_to_local_(shared_bnd_faces(colon(), i));
      const auto adjc_hexs = patch_x->f2h_->query(&corner[0]);
      left = ( adjc_hexs.first == -1 ) ? adjc_hexs.second : adjc_hexs.first;
    }
    {
      mati_t corner = patch_y->global_to_local_(shared_bnd_faces(colon(), i));
      const auto adjc_hexs = patch_y->f2h_->query(&corner[0]);
      right = ( adjc_hexs.first == -1 ) ? adjc_hexs.second : adjc_hexs.first;
    }
    corres_elems(0, i) = left;
    corres_elems(1, i) = right;
    ++i;
  }

  return 0;
}

stencil_nonconforming_d3::stencil_nonconforming_d3(const mati_t &hexs_H, const matd_t &rest_H,
                                                   const vector<shared_ptr<hexs_stencil>> &stencils,
                                                   const vector<shared_ptr<fine_mesh_patch>> &patches,
                                                   const double w)
    : w_(w), dim_(rest_H.size()), rest_H_(rest_H), stencils_(stencils), patches_(patches) {
  const size_t rd = rest_H.size(1);

  for (size_t i = 0; i < patches.size(); ++i) {
    if ( !patches[i]->Pro_ ) {
      patches[i]->Pro_ = make_shared<SparseMatrix<double>>();
      get_prolongation_matrix(stencils[i].get(), patches[i].get(), *patches[i]->Pro_);
    }
  }

  const double len = norm(rest_H(colon(), hexs_H(0, 0))-rest_H(colon(), hexs_H(1, 0)));
  const double area = len*len;

  for (size_t i = 0; i < stencils_.size(); ++i) {
    for (size_t j = i+1; j < stencils_.size(); ++j) {
      const auto &stenI = stencils_[i];
      const auto &stenJ = stencils_[j];
      const auto &ptchI = patches_[i];
      const auto &ptchJ = patches_[j];

      //-> check if stenicl i, j are adjacent
      mati_t shared_bs_pnt;
      if ( !is_adjc_quad_stencil(stenI, stenJ, shared_bs_pnt) )
        continue;
      
      //-> stencil i and j is adjacent
      adjc_sten_.push_back(make_pair(i, j));
      const mati_t ptch_node_i = ptchI->new_to_orig_(ptchI->vert_strip_);
      const mati_t ptch_node_j = ptchJ->new_to_orig_(ptchJ->vert_strip_);
      mati_t shared_nd_pnt;
      bigbang::find_intersections(ptch_node_i, ptch_node_j, shared_nd_pnt);
      shared_fine_nods_.push_back(make_shared<mati_t>(shared_nd_pnt));

      mati_t corres_hexs, bnd_faces;
      find_shared_boundary_element(dynamic_cast<fine_hexs_patch*>(ptchI.get()),
                                   dynamic_cast<fine_hexs_patch*>(ptchJ.get()),
                                   bnd_faces, corres_hexs);

      // //-> integral area
      // elem_bnd_len_.push_back(area*bnd_faces.size(2));

      unordered_Eigen::Map<size_t, size_t> bnd_to_local;
      for (size_t cnt = 0; cnt < shared_nd_pnt.size(); ++cnt)
        bnd_to_local.insert(make_pair(shared_nd_pnt[cnt], cnt));
      VectorXi node_count = VectorXi::Zero(3*shared_nd_pnt.size());
      for (const auto &fp : bnd_faces) {
        node_count.segment<3>(3*bnd_to_local[fp]) += Vector3i::Ones();
      }
      ASSERT(node_count.maxCoeff() == 4);
      VectorXd vono_area = VectorXd::Zero(3*shared_nd_pnt.size());
      for (size_t cnt = 0; cnt < vono_area.size(); ++cnt)
        vono_area[cnt] = sqrt(area/4.0*node_count[cnt]);
      A_.push_back(make_shared<VectorXd>(vono_area));
      
      VectorXd left_stiff  = VectorXd::Zero(shared_nd_pnt.size());
      VectorXd right_stiff = VectorXd::Zero(shared_nd_pnt.size());      
      for (size_t f = 0; f < corres_hexs.size(2); ++f) {
        const size_t left_elem = corres_hexs(0, f), right_elem = corres_hexs(1, f);      
        const double mu1 = ptchI->lame_(0, left_elem),  lambda1 = ptchI->lame_(1, left_elem);              
        const double mu2 = ptchJ->lame_(0, right_elem), lambda2 = ptchJ->lame_(1, right_elem);
        
        double ym1, pr1, ym2, pr2;
        compute_young_poisson(mu1, lambda1, ym1, pr1);
        compute_young_poisson(mu2, lambda2, ym2, pr2);
        double k1, k2;
        k1 = ym1*len;
        k2 = ym2*len;

        const size_t ei = bnd_faces(0, f);
        const size_t ej = bnd_faces(1, f);
        const size_t ep = bnd_faces(2, f);
        const size_t eq = bnd_faces(3, f);

        // cout << "edge num: " << f << endl;
        // cout << "ei ej: " << ei << ", " << ej << endl;
        // cout << "mu1 mu2: " << mu1 << ", " << mu2 << endl;

        int pos_ei = std::find(shared_nd_pnt.begin(), shared_nd_pnt.end(), ei)-shared_nd_pnt.begin();
        int pos_ej = std::find(shared_nd_pnt.begin(), shared_nd_pnt.end(), ej)-shared_nd_pnt.begin();
        int pos_ep = std::find(shared_nd_pnt.begin(), shared_nd_pnt.end(), ep)-shared_nd_pnt.begin();
        int pos_eq = std::find(shared_nd_pnt.begin(), shared_nd_pnt.end(), eq)-shared_nd_pnt.begin();

        // cout << "pos ei ej: " << pos_ei << ", " << pos_ej << endl << endl;
        // getchar();
        // sk->segment<3>(3*pos_ei) += homo_mu*Vector3d::Ones();
        // sk->segment<3>(3*pos_ej) += homo_mu*Vector3d::Ones();
        // sk->segment<3>(3*pos_ep) += homo_mu*Vector3d::Ones();
        // sk->segment<3>(3*pos_eq) += homo_mu*Vector3d::Ones();

        left_stiff[pos_ei] += k1; right_stiff[pos_ei] += k2;
        left_stiff[pos_ej] += k1; right_stiff[pos_ej] += k2;
        left_stiff[pos_ep] += k1; right_stiff[pos_ep] += k2;
        left_stiff[pos_eq] += k1; right_stiff[pos_eq] += k2;
      }

      shared_ptr<VectorXd> sk = make_shared<VectorXd>();
      *sk = VectorXd::Zero(3*shared_nd_pnt.size());
      for (size_t k = 0; k < shared_nd_pnt.size(); ++k)
        (*sk)[3*k+0] = (*sk)[3*k+1] = (*sk)[3*k+2] = left_stiff[k]*right_stiff[k]/(left_stiff[k]+right_stiff[k]);
      sk->cwiseSqrt();
      K_.push_back(sk);
      
      #if 0
      cout << "share nodes of patch " << i << ", " << j << endl;
      cout << trans(shared_nd_pnt) << endl;
      #endif

      const mati_t fixedI = ptchI->global_to_local_(shared_nd_pnt);
      const mati_t fixedJ = ptchJ->global_to_local_(shared_nd_pnt);
      ASSERT(fixedI.size() == fixedJ.size());

      #if 0
      cout << "fixedI " << trans(ptchI->vert_strip_perm_(fixedI)) << endl;
      cout << "fixedJ " << trans(ptchJ->vert_strip_perm_(fixedJ)) << endl;
      getchar();
      #endif

      shared_ptr<SparseMatrix<double>> Si = make_shared<SparseMatrix<double>>();
      shared_ptr<SparseMatrix<double>> Sj = make_shared<SparseMatrix<double>>();
      ASSERT(build_selection_matrix(fixedI, ptchI->nods_.size(2), rd, *Si) == 0);
      ASSERT(build_selection_matrix(fixedJ, ptchJ->nods_.size(2), rd, *Sj) == 0);
      S_.push_back(make_pair(Si, Sj));
    }
  }
}
size_t stencil_nonconforming_d3::Nx() const {
  return dim_;
}
int stencil_nonconforming_d3::Val(const double *x, double *val) const {
  RETURN_WITH_COND_TRUE(w_ == 0.0);  
  itr_matrix<const double *> X(3, Nx()/3, x);

  //-> clear the buffer to save the conforming energy
  adjc_sten_E_ = zeros<double>(adjc_sten_.size(), 1);
  
  for (size_t p = 0; p < adjc_sten_.size(); ++p) {
    const size_t I = adjc_sten_[p].first;
    const size_t J = adjc_sten_[p].second;    
    const auto &SI = *S_[p].first;
    const auto &SJ = *S_[p].second;
    
    matd_t xi_h; {
      const auto &sten  = stencils_[I];
      const auto &patch = patches_[I];
      const auto &P = *patch->Pro_;
      
      const matd_t x_H = X(colon(), sten->adjc_elem_);
      matd_t R = zeros<double>(3, 3);
      sten->query_local_frame(x_H, 0, R);

      const matd_t u_H = trans(R)*x_H-rest_H_(colon(), sten->adjc_elem_);
      matd_t u_h(patch->REST_.size(1), patch->REST_.size(2));
      Eigen::Map<VectorXd>(&u_h[0], u_h.size()) = P*Eigen::Map<const VectorXd>(&u_H[0], u_H.size());
      xi_h = R*(patch->REST_+u_h);      
    }
    
    matd_t xj_h; {
      const auto &sten  = stencils_[J];
      const auto &patch = patches_[J];
      const auto &P = *patch->Pro_;
      
      const matd_t x_H = X(colon(), sten->adjc_elem_);
      matd_t R = zeros<double>(3, 3);
      sten->query_local_frame(x_H, 0, R);

      const matd_t u_H = trans(R)*x_H-rest_H_(colon(), sten->adjc_elem_);
      matd_t u_h(patch->REST_.size(1), patch->REST_.size(2));
      Eigen::Map<VectorXd>(&u_h[0], u_h.size()) = P*Eigen::Map<const VectorXd>(&u_H[0], u_H.size());
      xj_h = R*(patch->REST_+u_h);
    }

    //-> dynamic update stiffness
    size_t sp = 0;
    for (const auto &sn : *shared_fine_nods_[p]) {
      const size_t idx_in_i = patches_[I]->global_to_local_[sn];
      const size_t idx_in_j = patches_[J]->global_to_local_[sn];
      const MatrixXd Ki = patches_[I]->K_.block(3*idx_in_i, 3*idx_in_i, 3, 3);
      const MatrixXd Kj = patches_[J]->K_.block(3*idx_in_j, 3*idx_in_j, 3, 3);

      Vector3d d;
      itr_matrix<double*>(3, 1, d.data()) = xi_h(colon(), idx_in_i)-xj_h(colon(), idx_in_j);
      if ( d.norm() < 1e-5 ) {
        K_[p]->segment<3>(3*sp).setZero();
        continue;
      }
      d.normalize();

      // cout << Ki << endl << endl;
      // cout << Kj << endl << endl;
      // VectorXd d2 = Vector3d::Random().normalized();
      // cout << Ki*d2 << endl << endl;
      // cout << Kj*(-d2) << endl << endl;
      // getchar(); getchar();
      
      const double stiff_i = 0.5*d.dot(Ki*d), stiff_j = 0.5*d.dot(Kj*d);
      // cout << stiff_i << ", " << stiff_j << endl;
      // getchar();getchar();
      
      const double homo_k = (stiff_i*stiff_j)/(stiff_i+stiff_j);
      K_[p]->segment<3>(3*sp) = sqrt(homo_k)*Vector3d::Ones();
      ++sp;
    }

    VectorXd metric = K_[p]->cwiseProduct(*A_[p]);
    const SparseMatrix<double> K(metric.asDiagonal());
    adjc_sten_E_[p] = 0.5*w_*(K*SI*Eigen::Map<const VectorXd>(&xi_h[0], xi_h.size())-
                                               K*SJ*Eigen::Map<const VectorXd>(&xj_h[0], xj_h.size())).squaredNorm();
    *val += adjc_sten_E_[p];
  
  }
  
  return 0;
}
int stencil_nonconforming_d3::Gra(const double *x, double *gra) const {
  RETURN_WITH_COND_TRUE(w_ == 0.0);
  itr_matrix<const double *> X(3, Nx()/3, x);
  itr_matrix<double *> G(3, Nx()/3, gra);

  for (size_t p = 0; p < adjc_sten_.size(); ++p) {
    const size_t I = adjc_sten_[p].first;
    const size_t J = adjc_sten_[p].second;
    const auto &Si = *S_[p].first;
    const auto &Sj = *S_[p].second;
    
    matd_t xi_h, Ri; {
      const auto &sten = stencils_[I];
      const auto &patch = patches_[I];
      const auto &P = *patch->Pro_;
      
      const matd_t x_H = X(colon(), sten->adjc_elem_);
      Ri = zeros<double>(3, 3);
      sten->query_local_frame(x_H, 0, Ri);

      const matd_t u_H = trans(Ri)*x_H-rest_H_(colon(), sten->adjc_elem_);
      matd_t u_h(patch->REST_.size(1), patch->REST_.size(2));
      Eigen::Map<VectorXd>(&u_h[0], u_h.size()) = P*Eigen::Map<const VectorXd>(&u_H[0], u_H.size());
      xi_h = Ri*(patch->REST_+u_h);      
    }
    
    matd_t xj_h, Rj; {
      const auto &sten = stencils_[J];
      const auto &patch = patches_[J];
      const auto &P = *patch->Pro_;
      
      const matd_t x_H = X(colon(), sten->adjc_elem_);
      Rj = zeros<double>(3, 3);
      sten->query_local_frame(x_H, 0, Rj);

      const matd_t u_H = trans(Rj)*x_H-rest_H_(colon(), sten->adjc_elem_);
      matd_t u_h(patch->REST_.size(1), patch->REST_.size(2));
      Eigen::Map<VectorXd>(&u_h[0], u_h.size()) = P*Eigen::Map<const VectorXd>(&u_H[0], u_H.size());
      xj_h = Rj*(patch->REST_+u_h);
    }

    //-> dynamic update stiffness
    size_t sp = 0;
    for (const auto &sn : *shared_fine_nods_[p]) {
      const size_t idx_in_i = patches_[I]->global_to_local_[sn];
      const size_t idx_in_j = patches_[J]->global_to_local_[sn];
      const MatrixXd Ki = patches_[I]->K_.block(3*idx_in_i, 3*idx_in_i, 3, 3);
      const MatrixXd Kj = patches_[J]->K_.block(3*idx_in_j, 3*idx_in_j, 3, 3);
      Vector3d d;
      itr_matrix<double*>(3, 1, d.data()) = xi_h(colon(), idx_in_i)-xj_h(colon(), idx_in_j);
      if ( d.norm() < 1e-5 ) {
        K_[p]->segment<3>(3*sp).setZero();
        continue;
      }
      d.normalize();
      const double stiff_i = 0.5*d.dot(Ki*d), stiff_j = 0.5*d.dot(Kj*d);
      const double homo_k = (stiff_i*stiff_j)/(stiff_i+stiff_j);
      K_[p]->segment<3>(3*sp) = sqrt(homo_k)*Vector3d::Ones();
      ++sp;
    }

    VectorXd metric = K_[p]->cwiseProduct(*A_[p]);
    const SparseMatrix<double> K(metric.asDiagonal());
    const VectorXd resd = K*(Si*Eigen::Map<const VectorXd>(&xi_h[0], xi_h.size())
                            -Sj*Eigen::Map<const VectorXd>(&xj_h[0], xj_h.size()));
    
    VectorXd gI, gJ;
    {
      const size_t sten_num = stencils_[I]->adjc_elem_.size();
      const size_t ptch_num = patches_[I]->nods_.size(2);
      SparseMatrix<double> RH, Rh;
      block_diagonalize_matrix(Ri, sten_num, RH);
      block_diagonalize_matrix(Ri, ptch_num, Rh);
      const auto jacI = K*Si*Rh*(*patches_[I]->Pro_)*RH.transpose();
      gI = w_*jacI.transpose()*resd;
    }
    {
      const size_t sten_num = stencils_[J]->adjc_elem_.size();
      const size_t ptch_num = patches_[J]->nods_.size(2);
      SparseMatrix<double> RH, Rh;
      block_diagonalize_matrix(Rj, sten_num, RH);
      block_diagonalize_matrix(Rj, ptch_num, Rh);
      const auto jacJ = K*Sj*Rh*(*patches_[J]->Pro_)*RH.transpose();
      gJ = -w_*jacJ.transpose()*resd;
    }
    
    G(colon(), stencils_[I]->adjc_elem_) += itr_matrix<const double *>(3, gI.size()/3, gI.data());
    G(colon(), stencils_[J]->adjc_elem_) += itr_matrix<const double *>(3, gJ.size()/3, gJ.data());
  }
  
  return 0;
}
int stencil_nonconforming_d3::Hes(const double *x, vector<Triplet<double>> *hes) const {
  RETURN_WITH_COND_TRUE(w_ == 0.0);
  itr_matrix<const double *> X(3, Nx()/3, x);

  for (size_t p = 0; p < adjc_sten_.size(); ++p) {
    const size_t I = adjc_sten_[p].first;
    const size_t J = adjc_sten_[p].second;    
    const auto &Si = *S_[p].first;
    const auto &Sj = *S_[p].second;
    
    matd_t xi_h, Ri; {
      const auto &sten = stencils_[I];
      const auto &patch = patches_[I];
      const auto &P = *patch->Pro_;
      
      const matd_t x_H = X(colon(), sten->adjc_elem_);
      Ri = zeros<double>(3, 3);
      sten->query_local_frame(x_H, 0, Ri);
    }
    
    matd_t xj_h, Rj; {
      const auto &sten = stencils_[J];
      const auto &patch = patches_[J];
      const auto &P = *patch->Pro_;
      
      const matd_t x_H = X(colon(), sten->adjc_elem_);
      Rj = zeros<double>(3, 3);
      sten->query_local_frame(x_H, 0, Rj);
    }

    //-> dynamic update stiffness
    size_t sp = 0;
    for (const auto &sn : *shared_fine_nods_[p]) {
      const size_t idx_in_i = patches_[I]->global_to_local_[sn];
      const size_t idx_in_j = patches_[J]->global_to_local_[sn];
      const MatrixXd Ki = patches_[I]->K_.block(3*idx_in_i, 3*idx_in_i, 3, 3);
      const MatrixXd Kj = patches_[J]->K_.block(3*idx_in_j, 3*idx_in_j, 3, 3);
      Vector3d d;
      itr_matrix<double*>(3, 1, d.data()) = xi_h(colon(), idx_in_i)-xj_h(colon(), idx_in_j);
      if ( d.norm() < 1e-5 ) {
        K_[p]->segment<3>(3*sp).setZero();
        continue;
      }
      d.normalize();
      const double stiff_i = 0.5*d.dot(Ki*d), stiff_j = 0.5*d.dot(Kj*d);
      const double homo_k = (stiff_i*stiff_j)/(stiff_i+stiff_j);
      K_[p]->segment<3>(3*sp) = sqrt(homo_k)*Vector3d::Ones();
      ++sp;
    }

    VectorXd metric = K_[p]->cwiseProduct(*A_[p]);
    const SparseMatrix<double> K(metric.asDiagonal());
    SparseMatrix<double> jacI, jacJ;
    {
      const size_t sten_num = stencils_[I]->adjc_elem_.size();
      const size_t ptch_num = patches_[I]->nods_.size(2);
      SparseMatrix<double> RH, Rh;
      block_diagonalize_matrix(Ri, sten_num, RH);
      block_diagonalize_matrix(Ri, ptch_num, Rh);
      jacI = K*Si*Rh*(*patches_[I]->Pro_)*RH.transpose();
    }
    {
      const size_t sten_num = stencils_[J]->adjc_elem_.size();
      const size_t ptch_num = patches_[J]->nods_.size(2);
      SparseMatrix<double> RH, Rh;
      block_diagonalize_matrix(Rj, sten_num, RH);
      block_diagonalize_matrix(Rj, ptch_num, Rh);
      jacJ = -K*Sj*Rh*(*patches_[J]->Pro_)*RH.transpose();
    }

    const size_t idx[2] = {I, J};
    const vector<SparseMatrix<double>*> Jac{&jacI, &jacJ};

    for (size_t m = 0; m < 2; ++m) {
      for (size_t n = 0; n < 2; ++n) {
        const SparseMatrix<double> H = w_*Jac[m]->transpose()*(*Jac[n]);

        for (size_t k = 0; k < H.outerSize(); ++k) {
          for (SparseMatrix<double>::InnerIterator it(H, k); it; ++it) {
            const size_t row = it.row(), col = it.col();
            const size_t II = 3*stencils_[idx[m]]->adjc_elem_[row/3]+row%3;
            const size_t JJ = 3*stencils_[idx[n]]->adjc_elem_[col/3]+col%3;
            hes->push_back(Triplet<double>(II, JJ, it.value()));
          }
        }
      }
    }    
  }
  return 0;
}

}

// #if 0
//   do {
//     VectorXd tmp_x = VectorXd::Random(opt_con->Nx());
    
//     VectorXd cv = VectorXd::Zero(opt_con->Nf());
//     opt_con->Val(tmp_x.data(), cv.data());

//     SparseMatrix<double> Jc(opt_con->Nf(), opt_con->Nx()); {
//       vector<Triplet<double>> trips;
//       opt_con->Jac(nullptr, 0, &trips);
//       Jc.setFromTriplets(trips.begin(), trips.end());
//     }

//     SparseQR<SparseMatrix<double>, COLAMDOrdering<int>> qr_solver;
//     qr_solver.compute(Jc);
//     const size_t rank_hes = qr_solver.rank();
//     cout << "# [OPT BASIS] system size and rank: " << Jc.cols() << " " << rank_hes << endl;

//     getchar();
//   } while (0);
// #endif

  //    cbf[2] = make_shared<full_aug_anis_zero_bnd_cons>(*this);

  //  try {
  //    opt_obj = make_shared<energy_t<double>>(ebf);
  //  } catch ( ... ) {
  //    ASSERT(0);
  //  }
  // static size_t count = 0;
  // if ( count >= 94 && count <= 98 ) {
  //   //    cout << "\tnull ratio: " << KER.nonZeros()*1.0/KER.rows()/KER.cols() << endl;
  //   cout << "\t@residual norm: " << (Jc*spN-cv).norm() << endl;
  //   cout << "\t@kernel size: " << KER.rows() << ", " << KER.cols() << endl;
  //   cout << "\t@test kernel: " << (Jc*KER).norm() << endl;
  //   getchar();
  // }
  // ++count;
  //-> turn constraints into soft penalties
  //  vector<shared_ptr<Functional<double>>> ebf(6); {
    //    ebf[0] = make_shared<least_square_wrapper<double>>(cbf[KRONECKER], 1e8);
    //    ebf[1] = make_shared<least_square_wrapper<double>>(cbf[TRANSLATION], 1e8);
    //    ebf[2] = make_shared<least_square_wrapper<double>>(cbf[ROTATION], 1e8);
    //    ebf[3] = make_shared<least_square_wrapper<double>>(cbf[HARMONIC], 1e8);
    //    ebf[4] = make_shared<least_square_wrapper<double>>(ip_cons, 1e-2);
    //    ebf[5] = make_shared<full_aug_anis_reg_energy>(*this, this->FEM_basis_, 1e2);
  //  }
    //-> used to interpolate deformable modes

  // //-> used to interpolate eigenvectors
  // MatrixXd eig_modes;
  // if ( gl_eig_disp ) {
  //   eig_modes = MatrixXd::Zero(nods_.size(), gl_eig_disp->cols());
  //   for (size_t i = 0; i < eig_modes.cols(); ++i) {
  //     itr_matrix<const double *> full_mode(nods_.size(1), gl_eig_disp->rows()/nods_.size(1), &(*gl_eig_disp)(0, i));
  //     itr_matrix<double *>(nods_.size(1), nods_.size(2), &eig_modes(0, i)) = full_mode(colon(), new_to_orig_);
  //   }
  //   const VectorXd eig_wgts = *gl_eig_wgts;
  //   cbf[MODAL] = make_shared<full_aug_anis_interp_cons>(*this, pt_, &eig_modes, &eig_wgts);
  // }

// //-> use zjucad optimizer
// shared_ptr<hj::math_func::math_func_t<double, int32_t>> objE;
// objE = make_shared<func_transformer<double, int32_t>>(opt_obj, &basis[0]);
// boost::property_tree::ptree opts; {
//   opts.put("package.value", "jtf");
//   opts.put("alg.value", "SQP");
//   opts.put("iter.value", 10);
//   opts.put("epsg.value", 1e-5);
//   opts.put("linear_solver/type.value", "PETsc");
// }
// jtf::optimize(*objE, basis, opts, nullptr, nullptr, NULL);

// shared_ptr<Constraint<double>> null_con;  
// SmartPtr<ipopt_opt_framework> opt_prb = new ipopt_opt_framework(opt_obj, opt_con, &basis[0]);
// SmartPtr<IpoptApplication> app = IpoptApplicationFactory();
// ApplicationReturnStatus status = app->Initialize();
// ASSERT(status == Ipopt::Solve_Succeeded);
// status = app->OptimizeTNLP(opt_prb);
//  ASSERT(status == Ipopt::Solve_Succeeded);
  
// // //-> solve GEV on galerkin projected matrix
// const SparseMatrix<double> PTKP = P.transpose()*K*P;
// SparseMatrix<double> M_H(rd*basis_num, rd*basis_num); {
//   vector<Triplet<double>> trips;
//   for (size_t i = 0; i < basis_pnt_mass_.size(); ++i) {
//     runtime_dim_add_diag_block(rd, i, i, basis_pnt_mass_[i], &trips);
//   }
//   M_H.reserve(trips.size());
//   M_H.setFromTriplets(trips.begin(), trips.end());
//   //    cout << "# total local mass coarse: " << M_H.sum() << endl;
// }
// cout << "M on patch: " <<  endl << M_H << endl;
// VectorXd eigs_H;
// solve_gen_eig_prob(PTKP, M_H, eignum, 0, "arpaca", &eigs_H, NULL);

// cout << eigs_H.transpose() << endl;
// getchar();

//  ASSERT(eigs_h.size() == eigs_H.size());

// //-> write out local spectrum
// static size_t count = 0;
// char outfile[256];
// sprintf(outfile, "%s/local_eigs_stencil_%04zu.csv", pt_.get<string>("outdir.value").c_str(), count++);

// ofstream ofs(outfile);
// ofs << "F, C" << endl;
// for (size_t i = 0; i < eigs_H.size(); ++i)
//   ofs << eigs_h[i] << ", " << eigs_H[i] << endl;
// ofs.close();
  
// matd_t disp_H = 0.5*ones<double>(rd, basis_pnt_.size());
// matd_t disp_h = zeros<double>(rd, verts_num);
// Eigen::Map<VectorXd>(&disp_h[0], disp_h.size()) = P*Eigen::Map<const VectorXd>(&disp_H[0], disp_H.size());

// static int count = 0;
// char outf[256];
// if ( cell_.size(1) == 4 ) {
//   sprintf(outf, "%s/mesh_patch_origin_%04d.vtk", pt_.get<string>("outdir.value").c_str(), count);
//   quad_mesh_write_to_vtk(outf, nods_, cell_);
//   sprintf(outf, "%s/mesh_patch_interp_%04d.vtk", pt_.get<string>("outdir.value").c_str(), count);
//   const matd_t nods_h = disp_h+nods_;
//   quad_mesh_write_to_vtk(outf, nods_h, cell_);
//   ++count;
// } else if ( cell_.size(1) == 8 ) {
//   sprintf(outf, "%s/mesh_patch_origin_%04d.vtk", pt_.get<string>("outdir.value").c_str(), count);
//   hex_mesh_write_to_vtk(outf, nods_, cell_);
//   sprintf(outf, "%s/mesh_patch_interp_%04d.vtk", pt_.get<string>("outdir.value").c_str(), count);
//   const matd_t nods_h = disp_h+nods_;
//   hex_mesh_write_to_vtk(outf, nods_h, cell_);
//   ++count;
// }

  
// //-> calc lumped area for each basis point
// basis_pnt_mass_ = zeros<double>(basis_pnt_.size(), 1);
// mati_t basis_to_local = -1*ones<size_t>(vert_strip_.size(), 1);
// basis_to_local(basis_pnt_) = colon(0, basis_pnt_.size()-1);
// for (size_t i = 0; i < dom.faces_.size(); ++i) {
//   const size_t fid = dom.faces_[i];

//   const double len = norm(nods_h(colon(), quad_H(0, fid))-nods_h(colon(), quad_H(1, fid)));
//   const double area = len*len;
    
//   const mati_t pid = global_to_local(quad_H(colon(), fid));
//   basis_pnt_mass_(basis_to_local(pid)) += area/4;
// }

// void fine_quad_patch::optimize_basis(matd_t &basisX, matd_t &basisY) {
//   size_t mode_num = pt_.get<size_t>("mode_num.value", this->modes_.cols());

//   //-> assemble energies for basis optimization
//   const int wei_option = pt_.get<int>("basis_weight.value");
//   shared_ptr<Constraint<double>> interp_cons = make_shared<patch_interp_cons_2d>(*this, mode_num, wei_option);
//   vector<shared_ptr<Functional<double>>> ebf(3); {
//     //    ebf[0] = make_shared<partition_to_unity_energy>(*this, 1e3);
//     ebf[1] = make_shared<kronecker_delta_energy>(*this, 1e3);
//     ebf[2] = make_shared<least_square_wrapper<double>>(interp_cons, 1e0);
//   }  
//   shared_ptr<Functional<double>> opt_obj;
//   try {
//     opt_obj = make_shared<energy_t<double>>(ebf);
//   } catch ( ... ) {
//     ASSERT(0);
//   }

//   matd_t basis = ones<double>(opt_obj->Nx(), 1);

// #if 0
//   do {
//     //-> check the rank of the hessian
//     SparseMatrix<double> Hes(opt_obj->Nx(), opt_obj->Nx()); {
//       vector<Triplet<double>> trips;
//       opt_obj->Hes(&basis[0], &trips);
//       Hes.setFromTriplets(trips.begin(), trips.end());
//     }
//     SparseQR<SparseMatrix<double>, COLAMDOrdering<int>> qr_solver;
//     qr_solver.compute(Hes);
//     cout << "# matrix size: " << Hes.cols() << endl;
//     cout << "# matrix rank: " << qr_solver.rank() << endl;
//     getchar();
//   } while (0);
// #endif

//   if ( pt_.get<string>("basis_optimizer.value") == "SQP" ) {
//     shared_ptr<hj::math_func::math_func_t<double, int32_t>> objE;
//     objE = make_shared<func_transformer<double, int32_t>>(opt_obj, &basis[0]);
//     boost::property_tree::ptree opts; {
//       opts.put("package.value", "jtf");
//       opts.put("alg.value", "SQP");
//       opts.put("iter.value", 20);
//       opts.put("epsg.value", 1e-8);    // the default epsx is 1e-20
//       opts.put("linear_solver/type.value", "PETsc");
//     }
//     jtf::optimize(*objE, basis, opts, nullptr, nullptr, NULL);
//   } else if ( pt_.get<string>("basis_optimizer.value") == "LBFGS" ) {
//     alglib_solver_args opt_args;
//     alglib_lbfgs_solve(opt_obj, basis.size(), &basis[0], opt_args);
//   }
  
//   const size_t vert_num_1d  = this->vert_strip_.size(2);
//   const size_t basis_num_1d = sqrt(this->basis_pnt_.size());  
//   basisX = itr_matrix<const double *>(basis_num_1d, vert_num_1d, &basis[0]);
//   basisY = itr_matrix<const double *>(basis_num_1d, vert_num_1d, &basis[0]+basis.size()/2);
// }

// void fine_quad_patch::analyse_dof() {
//   size_t mode_num = pt_.get<size_t>("mode_num.value", this->modes_.cols());
//   cout << "# mode_num: " << mode_num << endl;
//   cout << "# basis dim: " << basis_pnt_.size() << " " << nods_.size(2) << endl;
  
//   //-> CONSTRAINTS
//   vector<shared_ptr<Constraint<double>>> cbf(3); {
//     cbf[0] = make_shared<full_partition_to_unity_cons>(*this);
//     cbf[1] = make_shared<full_kronecker_delta_cons>(*this);
//     cbf[2] = make_shared<full_patch_interp_cons>(*this, mode_num, 0);
//   }
//   shared_ptr<Constraint<double>> opt_con;
//   try {
//     opt_con = make_shared<constraint_t<double>>(cbf);
//   } catch ( ... ) {
//     ASSERT(0);
//   }

//   matd_t basis = ones<double>(opt_con->Nx(), 1);
  
// #if 1
//   do {
//     //-> check the rank of the hessian
//     SparseMatrix<double> Jac(opt_con->Nf(), opt_con->Nx()); {
//       vector<Triplet<double>> trips;
//       opt_con->Jac(&basis[0], 0, &trips);
//       Jac.setFromTriplets(trips.begin(), trips.end());
//     }
//     SparseQR<SparseMatrix<double>, COLAMDOrdering<int>> qr_solver;
//     qr_solver.compute(Jac);
//     cout << "# matrix size: " << Jac.rows() << " " << Jac.cols() << endl;
//     cout << "# matrix rank: " << qr_solver.rank() << endl;
//     getchar();
//   } while (0);
// #endif

//   // //-> save out basis for a visualization
//   // const size_t verts_num = this->nods_.size(2), basis_num = this->basis_pnt_.size();
//   // basisXY = itr_matrix<const double *>(basis_num, verts_num, &basis[0]);
// }
  // //-> get prolongation matrix
  // P.resize(rd*verts_num, rd*basis_num); {
  //   vector<Triplet<double>> trips;
  //   for (size_t i = 0; i < vert_strip_.size(); ++i) {
  //     const size_t pid = vert_strip_[i];    
  //     for (size_t j = 0; j < basis_pnt_.size(); ++j) {
  //       for (size_t k = 0; k < rd; ++k) { // each dimension
  //         if ( !isotropic && !augmented ) {
  //           trips.push_back(Triplet<double>(rd*pid+k, rd*j+k, basis(j, i+k*verts_num)));
  //         } else if ( !isotropic && augmented ) {
  //           for (size_t m = 0; m < rd; ++m) {
  //             trips.push_back(Triplet<double>(rd*pid+k, rd*j+m, basis(j, i+(rd*k+m)*verts_num)));
  //           }
  //         } else if ( isotropic )
  //           trips.push_back(Triplet<double>(rd*pid+k, rd*j+k, basis(j, i)));
  //       }
  //     }
  //   }
  //   P.reserve(trips.size());
  //   P.setFromTriplets(trips.begin(), trips.end());
  // }
