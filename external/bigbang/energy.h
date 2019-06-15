#ifndef ENERGY_H
#define ENERGY_H

#include <zjucad/matrix/matrix.h>
#include <Eigen/Sparse>
#include <unordered_map>
#include <unordered_set>
#include <boost/functional/hash.hpp>

#include "def.h"

using mati_t=zjucad::matrix::matrix<size_t>;
using matd_t=zjucad::matrix::matrix<double>;

namespace bigbang {

class fem_2d_energy : public Functional<double>
{
public:
  enum MTR {LINEAR, STVK};
  fem_2d_energy(const mati_t &tris, const matd_t &nods) : tris_(tris), dim_(nods.size()) {}
  fem_2d_energy(const mati_t &tris, const matd_t &nods, const MTR mtr, const double mu, const double lam, const double w);
  virtual ~fem_2d_energy() {}
  size_t Nx() const;
  int Val(const double *x, double *val) const;
  int Gra(const double *x, double *gra) const;
  int Hes(const double *x, std::vector<Eigen::Triplet<double>> *hes) const;
public:
  Eigen::VectorXd elem_energy_, elem_strain_;
  Eigen::VectorXd elem_stress_;
protected:
  const mati_t &tris_;
  const size_t dim_;
  MTR mtr_;
  double lam_, miu_;
  double w_;
  matd_t area_;
  matd_t Dm_;
};

class poly_elastic_energy_2d : public fem_2d_energy
{
public:
  poly_elastic_energy_2d(const mati_t &tris, const matd_t &nods, const Eigen::VectorXd &coef,
                         const MTR mtr, const double mu, const double lam, const int reg_option,
                         const double w);
  size_t Nx() const;
  int Val(const double *x, double *val) const;
  int Gra(const double *x, double *gra) const;
  int Hes(const double *x, std::vector<Eigen::Triplet<double>> *hes) const;
private:
  const Eigen::VectorXd coef_;
  const int reg_option_;
};
  
class momentum_potential : public Functional<double>
{
public:
  ~momentum_potential() {}
  virtual size_t Nx() const = 0;
  virtual int Val(const double *x, double *val) const = 0;
  virtual int Gra(const double *x, double *gra) const = 0;
  virtual int Hes(const double *x, std::vector<Eigen::Triplet<double>> *hes) const = 0;
  virtual void Update(const double *x) = 0;
  virtual double QueryKineticEnergy() const = 0;
  virtual const Eigen::SparseMatrix<double>& MassMatrix() const = 0;
  virtual const Eigen::VectorXd& CurrVelocity() const = 0;
  virtual double timestep() const = 0;
  void ClearVelocity() { vn_.setZero(); }
  Eigen::VectorXd xn_, vn_;
};

class implicit_avf : public momentum_potential
{
public:
  implicit_avf(const mati_t &cell, const matd_t &nods, const double rho, const double h,
               const std::shared_ptr<Functional<double>> &Ep);
  size_t Nx() const;
  int Val(const double *x, double *val) const;
  int Gra(const double *x, double *gra) const;
  int Hes(const double *x, std::vector<Eigen::Triplet<double>> *hes) const;
  void Init(const double *x0, const double *v0);
  void Update(const double *x);
  double QueryKineticEnergy() const { return vn_.dot(M_*vn_)/2; }
  const Eigen::SparseMatrix<double>& MassMatrix() const { return M_; }
  const Eigen::VectorXd& CurrVelocity() const { return vn_; }
  double timestep() const { return h_; }
private:
  const double rho_, h_, h2_;
  const size_t dim_;
  const std::shared_ptr<Functional<double>> Ep_;
  Eigen::SparseMatrix<double> M_, invM_;
  Eigen::VectorXd fn_;
};

class momentum_potential_imp_euler : public momentum_potential
{
public:
  momentum_potential_imp_euler(const mati_t &cell, const matd_t &nods, const double rho, const double h, const double w=1.0);
  momentum_potential_imp_euler(const matd_t &nods, const Eigen::SparseMatrix<double> &M, const double h, const double w=1.0);
  size_t Nx() const;
  int Val(const double *x, double *val) const;
  int Gra(const double *x, double *gra) const;
  int Hes(const double *x, std::vector<Eigen::Triplet<double>> *hes) const;
  void Init(const double *x0, const double *v0);
  void Update(const double *x);
  double QueryKineticEnergy() const;
  const Eigen::SparseMatrix<double>& MassMatrix() const { return M_; }
  const Eigen::VectorXd& CurrVelocity() const { return vn_; }
  double timestep() const { return h_; }
private:
  double rho_, h_;
  const size_t dim_;
  double w_;
  Eigen::SparseMatrix<double> M_;
};

class momentum_potential_bdf2 : public momentum_potential
{
public:
  momentum_potential_bdf2(const mati_t &cell, const matd_t &nods, const double rho, const double h, const double w=1.0);
  momentum_potential_bdf2(const matd_t &nods, const Eigen::SparseMatrix<double> &M, const double h, const double w=1.0);
  size_t Nx() const;
  int Val(const double *x, double *val) const;
  int Gra(const double *x, double *gra) const;
  int Hes(const double *x, std::vector<Eigen::Triplet<double>> *hes) const;
  void Init(const double *x0, const double *v0);
  void Update(const double *x);
  double QueryKineticEnergy() const;
  const Eigen::SparseMatrix<double>& MassMatrix() const { return M_; }
  const Eigen::VectorXd& CurrVelocity() const { return vn_; }
  double timestep() const { return h_; }
private:
  double rho_, h_;
  const size_t dim_;
  double w_;
  Eigen::SparseMatrix<double> M_;
  Eigen::VectorXd xnn_, vnn_;
};

class elastic_potential : public Functional<double>
{
public:
  enum Material {
    LINEAR,
    STVK,
    COROTATIONAL,
    NEOHOOKEAN
  };
  elastic_potential(const mati_t &tets, const matd_t &nods, Material type, const double Ym, const double Pr, const double w);
  size_t Nx() const;
  int Val(const double *x, double *val) const;
  int Gra(const double *x, double *gra) const;
  int Hes(const double *x, std::vector<Eigen::Triplet<double>> *hes) const;
public:
  Material type_;
private:
  const size_t dim_;
  double w_;
  const mati_t &tets_;
  matd_t vol_;
  double lam_, miu_;
  matd_t Dm_;
  matd_t rest_;
};

class voxel_elastic_potential : public Functional<double>
{
public:
  enum Material {
    LINEAR,
    COROTATED,
    STVK,
    NEOHOOKEAN,
    STANEO,
    BOWERNEO
  };
  voxel_elastic_potential(const mati_t &cube, const matd_t &nods, Material type,
                          const matd_t &lame, const double w);
  size_t Nx() const;
  int Val(const double *x, double *val) const;
  int Gra(const double *x, double *gra) const;
  int Hes(const double *x, std::vector<Eigen::Triplet<double>> *hes) const;
  void update_mtr_params(const matd_t &new_params);
public:
  mutable matd_t elemE_;
private:
  const size_t dim_;
  double w_;
  Material type_;
  const mati_t &cube_;
  matd_t lame_;

  matd_t quadr_, qw_;
  matd_t detDmH_;
  matd_t H_invDmH_;
};

class quad49_hybrid_energy : public Functional<double>
{
public:
  enum Material {
    LINEAR,
    COROTATED,
    STVK,
    NEOHOOKEAN,
    STANEO,
    BOWERNEO
  };
  quad49_hybrid_energy(const mati_t &quad4, const matd_t &nods, Material type,
                       const matd_t &lame, const double w);
  size_t Nx() const;
  int Val(const double *x, double *val) const;
  int Gra(const double *x, double *gra) const;
  int Hes(const double *x, std::vector<Eigen::Triplet<double>> *hes) const;
protected:
  const size_t dim_;
  double w_;
  Material type_;
  const mati_t &quad4_;
  mati_t quad9_;
  matd_t lame_;

  matd_t quadr_, qw_;
  matd_t dets_;

  matd_t H_invDmH4_, H_invDmH9_;

  size_t QUAD_NUM_;
};

class quad4_taylor_energy : public quad49_hybrid_energy
{
public:
  quad4_taylor_energy(const mati_t &quad4, const matd_t &nods, Material type,
                      const matd_t &lame, const double w, const char hessian='N');
  size_t Nx() const;
  int Val(const double *x, double *val) const;
  int Gra(const double *x, double *gra) const;
  int Hes(const double *x, std::vector<Eigen::Triplet<double>> *hes) const;
  void set_hess_flag (const char flag);
  void set_debug_flag(const bool flag) { debug_flag_ = flag; }
  bool get_debug_flag() const { return debug_flag_; }
  const matd_t& get_elem_eig_vals() const;
public:
  mutable std::vector<Eigen::Triplet<double>> trips_bak_, trips_bak_M_;
private:
  int HesExactNewton(const double *x, std::vector<Eigen::Triplet<double>> *hes) const;
  int HesQuasiNewton(const double *x, std::vector<Eigen::Triplet<double>> *hes) const;
  int HesTruncNewton(const double *x, std::vector<Eigen::Triplet<double>> *hes) const;
  int HesStiffNewton(const double *x, std::vector<Eigen::Triplet<double>> *hes) const;
  int HesMixedNewton(const double *x, std::vector<Eigen::Triplet<double>> *hes) const;
private:
  char hessian_;
  matd_t opG4_, opG9_;
  mutable matd_t eigval_;
  bool debug_flag_;
  zjucad::matrix::matrix<int> in_quad_, in_quad9_;
  std::unordered_set<std::pair<int, int>, boost::hash<std::pair<int, int>>> sp_pattern_;
};

class quad4_admm_energy : public Functional<double>
{
public:
  enum Material {
    LINEAR,
    COROTATED,
    STVK,
    NEOHOOKEAN,
    STANEO,
    BOWERNEO
  };
  quad4_admm_energy(const mati_t &quad4, const matd_t &nods, Material type,
                    const matd_t &lame, const double w, const char hessian='N');
  size_t Nx() const;
  int Val(const double *x, double *val) const;
  int Gra(const double *x, double *gra) const;
  int Hes(const double *x, std::vector<Eigen::Triplet<double>> *hes) const;
  void set_hess_flag(const char flag);
  void set_debug_flag(const bool flag) { debug_flag_ = flag; }
  bool get_debug_flag() const { return debug_flag_; }
public:
  mutable matd_t exactH_, truncH_;
private:
  int HesExactNewton(const double *x, std::vector<Eigen::Triplet<double>> *hes) const;
  int HesTruncNewton(const double *x, std::vector<Eigen::Triplet<double>> *hes) const;
private:
  char hessian_;
  const size_t dim_;
  double w_;
  Material type_;
  const mati_t &quad4_;
  mati_t quad9_;
  matd_t lame_;

  matd_t quadr_, qw_;
  matd_t dets_;

  matd_t H_invDmH9_;
  size_t QUAD_NUM_;
  zjucad::matrix::matrix<int> in_quad_;

  bool debug_flag_=false;
};

class quad9_elastic_energy : public Functional<double>
{
public:
  enum Material {
    LINEAR,
    COROTATED,
    STVK,
    NEOHOOKEAN,
    STANEO,
    BOWERNEO
  };
  quad9_elastic_energy(const mati_t &quad9, const matd_t &nods, Material type,
                       const matd_t &lame, const double w);
  size_t Nx() const;
  int Val(const double *x, double *val) const;
  int Gra(const double *x, double *gra) const;
  int Hes(const double *x, std::vector<Eigen::Triplet<double>> *hes) const;
  void update_mtr_params(const matd_t &new_params);
private:
  const size_t dim_;
  double w_;
  Material type_;
  const mati_t &quad9_;
  matd_t lame_;

  matd_t quadr_, qw_;
  matd_t detDmH_;
  matd_t H_invDmH_;

  size_t QUAD_NUM_;
};

class quad4_elastic_energy : public Functional<double>
{
public:
  enum Material {
    LINEAR,
    COROTATED,
    STVK,
    NEOHOOKEAN,
    STANEO,
    BOWERNEO
  };
  quad4_elastic_energy(const mati_t &quad, const matd_t &nods, Material type,
                       const matd_t &lame, const double w);
  size_t Nx() const;
  int Val(const double *x, double *val) const;
  int Gra(const double *x, double *gra) const;
  int Hes(const double *x, std::vector<Eigen::Triplet<double>> *hes) const;
  void update_mtr_params(const matd_t &new_params);
private:
  const size_t dim_;
  double w_;
  Material type_;
  const mati_t &quad_;
  matd_t lame_;

  matd_t quadr_, qw_;
  matd_t detDmH_;
  matd_t H_invDmH_;

  size_t QUAD_NUM_;
};

class gravitational_potential : public Functional<double>
{
public:
  gravitational_potential(const mati_t &cell, const matd_t &nods, const double rho, const double w, const char axis='Y');
  gravitational_potential(const matd_t &nods, const Eigen::SparseMatrix<double> &M, const double w, const char axis='Y');
  size_t Nx() const;
  int Val(const double *x, double *val) const;
  int Gra(const double *x, double *gra) const;
  int Hes(const double *x, std::vector<Eigen::Triplet<double>> *hes) const { return __LINE__; }
  void update_weight(const double new_w) { w_ = new_w; }
  void update_direction(const char axis) { direction_ = axis-'X'; }
private:
  int direction_;
  const size_t dim_;
  const int rd_;
  double w_;
  Eigen::SparseMatrix<double> M_;
};

class torque_3d_potential : public Functional<double>
{
public:
  torque_3d_potential(const matd_t &nods, const Eigen::SparseMatrix<double> &M,
                      const Eigen::Vector3d &omega, const double w);
  size_t Nx() const;
  int Val(const double *x, double *val) const;
  int Gra(const double *x, double *gra) const;
  int Hes(const double *x, std::vector<Eigen::Triplet<double>> *hes) const;
  void update_weight(const double new_w) { w_ = new_w; }
  void update_omega(const double *new_omega) {
    omega_[0] = new_omega[0];
    omega_[1] = new_omega[1];
    omega_[2] = new_omega[2];
    omega_ /= norm(omega_);
  }
private:
  matd_t bc_, omega_;
  const matd_t &rest_;
  const size_t dim_;
  double w_;
  const Eigen::SparseMatrix<double> &M_;
};

class positional_potential : public Functional<double>
{
public:
  positional_potential(const matd_t &nods, const double w);
  size_t Nx() const;
  int Val(const double *x, double *val) const;
  int Gra(const double *x, double *gra) const;
  int Hes(const double *x, std::vector<Eigen::Triplet<double>> *hes) const;
  void update_fixed_verts(const std::vector<size_t> &pids, const double *pos);
  int Pin(const size_t id, const double *pos);
  int Release(const size_t id);
private:
  const size_t dim_, rd_;
  double w_;
  std::unordered_map<size_t, Eigen::Vector3d> fixed3d_;
  std::unordered_map<size_t, Eigen::Vector2d> fixed2d_;
};

class spring_potential : public Functional<double>
{
public:
  spring_potential(const mati_t &edge, const matd_t &nods, const double w, char option='M');
  size_t Nx() const;
  int Val(const double *x, double *val) const;
  int Gra(const double *x, double *gra) const;
  int Hes(const double *x, std::vector<Eigen::Triplet<double>> *hes) const;
  void ResetWeight(const double w) { w_ = w; }
  int ResetEdgeMaterial(const size_t p, const size_t q, const double scale_w, const double scale_len);
private:
  void Newton(const matd_t &vert, const double &rest_len, matd_t &H) const;
  void GaussNewton(const matd_t &vert, const double &rest_len, matd_t &H) const;
  void Mixed(const matd_t &vert, const double &rest_len, matd_t &H) const;
private:
  const size_t dim_;
  matd_t w_;
  const mati_t &edge_;
  matd_t len_;
  char option_;
};

class stvk_spring_potential : public Functional<double>
{
public:
  stvk_spring_potential(const mati_t &edge, const matd_t &nods, const double w);
  size_t Nx() const;
  int Val(const double *x, double *val) const;
  int Gra(const double *x, double *gra) const;
  int Hes(const double *x, std::vector<Eigen::Triplet<double>> *hes) const;
private:
  const size_t dim_;
  const double w_;
  const mati_t &edge_;
  matd_t len_;
};

class gauss_newton_spring : public Functional<double>
{
public:
  gauss_newton_spring(const mati_t &edge, const matd_t &nods, const double w);
  size_t Nx() const;
  int Val(const double *x, double *val) const;
  int Gra(const double *x, double *gra) const;
  int Hes(const double *x, std::vector<Eigen::Triplet<double>> *hes) const;
  void ResetWeight(const double w) { w_ = w; }
private:
  const size_t dim_;
  matd_t w_;
  const mati_t &edge_;
  matd_t len_;
};

template<typename T>
void hash_combine(size_t &seed, T const &v) {
  seed ^= std::hash<T>()(v)+0x9e3779b9+(seed<<6)+(seed>>2);
}

struct pair_hash {
public:
  template<typename T, typename U>
  size_t operator()(const std::pair<T, U> &rhs) const {
    size_t retval = std::hash<T>()(rhs.first);
    hash_combine(retval, rhs.second);
    return retval;
  }
};

class fast_mass_spring : public Functional<double>
{
public:
  fast_mass_spring(const mati_t &edge, const matd_t &nods, const double w);
  size_t Nx() const;
  int Val(const double *x, double *val) const;
  int Gra(const double *x, double *gra) const;
  int Hes(const double *x, std::vector<Eigen::Triplet<double>> *hes) const;
  void LocalSolve(const double *x);
  void Project();
  size_t aux_dim() const { return 3*edge_.size(2); }
  const double* get_aux_var() const { return d_.begin(); }
  const Eigen::SparseMatrix<double>& get_df_mat() const { return S_; }
  void build_jts_pattern();
  Eigen::SparseMatrix<double>& get_jts(const double *x);
public:
  const mati_t &edge_;
  matd_t len_;
  matd_t d_;
private:
  const size_t dim_;
  double w_;
  std::unordered_map<std::pair<size_t, size_t>, size_t, pair_hash> ijp_;
  Eigen::SparseMatrix<double> S_, JtS_;
};

class line_bending_potential : public Functional<double>
{
public:
  line_bending_potential(const mati_t &edge, const matd_t &nods, const double w);
  size_t Nx() const;
  int Val(const double *x, double *val) const;
  int Gra(const double *x, double *gra) const;
  int Hes(const double *x, std::vector<Eigen::Triplet<double>> *hes) const;
private:
  const size_t dim_;
  double w_;
  const mati_t &edge_;
  matd_t len_;
};

class surf_bending_potential : public Functional<double>
{
public:
  surf_bending_potential(const mati_t &diams, const matd_t &nods, const double w);
  size_t Nx() const;
  int Val(const double *x, double *val) const;
  int Gra(const double *x, double *gra) const;
  int Hes(const double *x, std::vector<Eigen::Triplet<double>> *hes) const;
  void ResetWeight(const double w) { w_ = w; }
private:
  const size_t dim_;
  double w_;
  const mati_t &diams_;
  matd_t len_, area_, angle_;
};

class ext_force_energy : public Functional<double>
{
public:
  ext_force_energy(const matd_t &nods, const double w);
  size_t Nx() const;
  int Val(const double *x, double *val) const;
  int Gra(const double *x, double *gra) const;
  int Hes(const double *x, std::vector<Eigen::Triplet<double>> *hes) const { return __LINE__; }
  int ApplyForce(const size_t id, const double *f);
  int RemoveForce(const size_t id);
private:
  const size_t dim_, rd_;
  double w_;
  matd_t force_;
};

class general_fext_energy : public Functional<double>
{
public:
  general_fext_energy(const size_t dim, const double w)
      : dim_(dim), w_(w) {
    fext_ = Eigen::VectorXd::Zero(dim);
  }
  size_t Nx() const {
    return dim_;
  }
  int Val(const double *x, double *val) const {
    Eigen::Map<const Eigen::VectorXd> X(x, dim_);
    *val += -w_*fext_.dot(X);
    return 0;
  }
  int Gra(const double *x, double *gra) const {
    Eigen::Map<Eigen::VectorXd> G(gra, dim_);
    G += -w_*fext_;
    return 0;
  }
  int Hes(const double *x, std::vector<Eigen::Triplet<double>> *hes) const {
    return __LINE__;
  }
  void update_ext_force(const Eigen::VectorXd &f) { fext_ = f; }
private:
  Eigen::VectorXd fext_;
  const size_t dim_;
  const double w_;
};

class isometric_bending : public Functional<double>
{
public:
  isometric_bending(const mati_t &diams, const matd_t &nods, const double w);
  size_t Nx() const;
  int Val(const double *x, double *val) const;
  int Gra(const double *x, double *gra) const;
  int Hes(const double *x, std::vector<Eigen::Triplet<double>> *hes) const;
private:
  const size_t dim_;
  const mati_t &diams_;
  double w_;
  matd_t cotv_, area_;
};

class tet_arap_energy : public Functional<double>
{
public:
  tet_arap_energy(const mati_t &tets, const matd_t &nods, const double w);
  size_t Nx() const;
  int Val(const double *x, double *val) const;
  int Gra(const double *x, double *gra) const;
  int Hes(const double *x, std::vector<Eigen::Triplet<double>> *hes) const;
  void LocalSolve(const double *x);
  void CalcLieAlgebraCoord(Eigen::VectorXd &vec) const;
  void UpdateRotation(const Eigen::VectorXd &vec);
  size_t aux_dim() const { return 9*tets_.size(2); }
  const double* get_aux_var() const { return R_.begin(); }
public:
  matd_t R_;
private:
  const size_t dim_;
  double w_;
  const mati_t &tets_;
  matd_t vol_;
  matd_t D_;
};

class bw98_stretch_energy : public Functional<double>
{
public:
  bw98_stretch_energy(const mati_t &tris, const matd_t &nods, const double w);
  size_t Nx() const;
  int Val(const double *x, double *val) const;
  int Gra(const double *x, double *gra) const;
  int Hes(const double *x, std::vector<Eigen::Triplet<double>> *hes) const;
private:
  const size_t dim_;
  double w_;
  const mati_t &tris_;
  matd_t area_;
  matd_t invUV_;
};

class bw98_shear_energy : public Functional<double>
{
public:
  bw98_shear_energy(const mati_t &tris, const matd_t &nods, const double w);
  size_t Nx() const;
  int Val(const double *x, double *val) const;
  int Gra(const double *x, double *gra) const;
  int Hes(const double *x, std::vector<Eigen::Triplet<double>> *hes) const;
private:
  const size_t dim_;
  double w_;
  const mati_t &tris_;
  matd_t area_;
  matd_t invUV_;
};

class fem_stretch_energy : public Functional<double>
{
public:
  fem_stretch_energy(const mati_t &tris, const matd_t &nods, const double w);
  size_t Nx() const;
  int Val(const double *x, double *val) const;
  int Gra(const double *x, double *gra) const;
  int Hes(const double *x, std::vector<Eigen::Triplet<double>> *hes) const;
private:
  const size_t dim_;
  double w_;
  const mati_t &tris_;
  matd_t area_, K_;
  matd_t Dm_;
};

class low_pass_filter_energy : public Functional<double>
{
public:
  low_pass_filter_energy(const mati_t &tris, const matd_t &nods, const size_t patch_num, const double w);
  size_t Nx() const;
  int Val(const double *x, double *val) const;
  int Gra(const double *x, double *gra) const;
  int Hes(const double *x, std::vector<Eigen::Triplet<double>> *hes) const;
  void Update(const double *ref);
private:
  const size_t dim_;
  double w_;
  const size_t patch_num_;
  const mati_t &tris_;
  const double *ref_;
  std::vector<std::vector<std::pair<size_t, double>>> pat_;
};

class cosserat_stretch_energy : public Functional<double>
{
public:
  cosserat_stretch_energy(const mati_t &rod, const matd_t &nods,
                          const double Es, const double r, const double w=1.0);
  size_t Nx() const;
  int Val(const double *xq, double *val) const;
  int Gra(const double *xq, double *gra) const;
  int Hes(const double *xq, std::vector<Eigen::Triplet<double>> *hes) const;
private:
  const mati_t &rod_;
  const double w_;
  const size_t r_size_, q_size_;
  const size_t elem_num_;
  const double Es_, r_;
  Eigen::VectorXd len_;
};

class cosserat_bend_energy : public Functional<double>
{
public:
  cosserat_bend_energy(const mati_t &rod, const matd_t &nods,
                       const double E, const double G, const double r, const double w=1.0);
  size_t Nx() const;
  int Val(const double *xq, double *val) const;
  int Gra(const double *xq, double *gra) const;
  int Hes(const double *xq, std::vector<Eigen::Triplet<double>> *hes) const;
private:
  const mati_t &rod_;
  const double w_;
  const size_t r_size_, q_size_;
  const size_t elem_num_;
  const double E_, G_, r_;
  Eigen::VectorXd len_;
};

class cosserat_couple_energy : public Functional<double>
{
public:
  cosserat_couple_energy(const mati_t &rod, const matd_t &nods,
                         const double kappa, const double w=1.0);
  size_t Nx() const;
  int Val(const double *xq, double *val) const;
  int Gra(const double *xq, double *gra) const;
  int Hes(const double *xq, std::vector<Eigen::Triplet<double>> *hes) const;
private:
  const mati_t &rod_;
  const double w_;
  const size_t r_size_, q_size_;
  const size_t elem_num_;
  const double kappa_;
  Eigen::VectorXd len_;
};

}

#endif
