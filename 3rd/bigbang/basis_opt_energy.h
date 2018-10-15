#ifndef BASIS_OPT_ENERGY_H
#define BASIS_OPT_ENERGY_H

#include "def.h"
#include "stencil_patch.h"

namespace bigbang {

class full_kronecker_delta_cons : public Constraint<double>
{
public:
  typedef Eigen::Triplet<double> TPL;
  full_kronecker_delta_cons(const fine_mesh_patch &patch);
  size_t Nx() const;
  size_t Nf() const;
  int Val(const double *x, double *val) const;
  int Jac(const double *x, const size_t off, std::vector<TPL> *jac) const;
  int Hes(const double *x, const size_t off, std::vector<std::vector<TPL>> *hes) const;
private:
  size_t verts_num_, basis_num_;
  mati_t basis_in_strip_;
};

class full_partition_to_unity_cons : public Constraint<double>
{
public:
  typedef Eigen::Triplet<double> TPL;
  full_partition_to_unity_cons(const fine_mesh_patch &patch);
  size_t Nx() const;
  size_t Nf() const;
  int Val(const double *x, double *val) const;
  int Jac(const double *x, const size_t off, std::vector<TPL> *jac) const;
  int Hes(const double *x, const size_t off, std::vector<std::vector<TPL>> *hes) const;
private:
  size_t verts_num_, basis_num_;
};

class full_patch_interp_cons : public Constraint<double>
{
public:
  typedef Eigen::Triplet<double> TPL;
  full_patch_interp_cons(const fine_mesh_patch &patch, const double power,
                         const boost::property_tree::ptree &pt,
                         const Eigen::MatrixXd *modes=nullptr);
  size_t Nx() const;
  size_t Nf() const;
  int Val(const double *x, double *val) const;
  int Jac(const double *x, const size_t off, std::vector<TPL> *jac) const;
  int Hes(const double *x, const size_t off, std::vector<std::vector<TPL>> *hes) const;
  double query_interp_err(const double *x, const double *mode) const;
private:
  size_t basis_num_, verts_num_;
  size_t mode_num_;
  size_t nods_sz1_, nods_sz2_, nods_sz_;
  const fine_mesh_patch &patch_;
  const double eps_;
  Eigen::VectorXd wgts_;
  Eigen::MatrixXd modes_;
};

class elemwise_polyharm_energy : public Functional<double>
{
public:
  elemwise_polyharm_energy(const fine_mesh_patch &patch, const double w);
  size_t Nx() const;
  int Val(const double *x, double *val) const;
  int Gra(const double *x, double *gra) const;
  int Hes(const double *x, std::vector<Eigen::Triplet<double>> *hes) const;
private:
  const double w_;
  const fine_mesh_patch &patch_;
  Eigen::SparseMatrix<double> K_;
  std::vector<std::shared_ptr<Eigen::SparseMatrix<double>>> G_;
  size_t basis_num_, verts_num_;
  const size_t rd_;
};

class full_anis_interp_cons : public Constraint<double>
{
public:
  typedef Eigen::Triplet<double> TPL;
  full_anis_interp_cons(const fine_mesh_patch &patch, const boost::property_tree::ptree &pt,
                        const Eigen::MatrixXd *modes=nullptr);
  size_t Nx() const;
  size_t Nf() const;
  int Val(const double *x, double *val) const;
  int Jac(const double *x, const size_t off, std::vector<TPL> *jac) const;
  int Hes(const double *x, const size_t off, std::vector<std::vector<TPL>> *hes) const;
  void reset_mode_number(const size_t number);
private:
  size_t basis_num_, verts_num_;
  size_t mode_num_;
  size_t nods_sz1_, nods_sz2_, nods_sz_;
  const fine_mesh_patch &patch_;
  const double eps_;
  Eigen::VectorXd wgts_;
  Eigen::MatrixXd modes_;
};

class full_aug_anis_interp_cons : public Constraint<double>
{
public:
  typedef Eigen::Triplet<double> TPL;
  full_aug_anis_interp_cons(const fine_mesh_patch &patch, const boost::property_tree::ptree &pt,
                            const Eigen::MatrixXd *modes, const Eigen::VectorXd *weights);
  size_t Nx() const;
  size_t Nf() const;
  int Val(const double *x, double *val) const;
  int Jac(const double *x, const size_t off, std::vector<TPL> *jac) const;
  int Hes(const double *x, const size_t off, std::vector<std::vector<TPL>> *hes) const;
private:
  size_t basis_num_, verts_num_;
  size_t mode_num_;
  size_t nods_sz1_, nods_sz2_, nods_sz_;
  const fine_mesh_patch &patch_;
  Eigen::VectorXd wgts_;
  Eigen::MatrixXd modes_;
};

class full_aug_anis_reg_energy : public Functional<double>
{
public:
  typedef Eigen::Triplet<double> TPL;
  full_aug_anis_reg_energy(const fine_mesh_patch &patch, const matd_t &refs, const double w);
  size_t Nx() const;
  int Val(const double *x, double *val) const;
  int Gra(const double *x, double *gra) const;
  int Hes(const double *x, std::vector<TPL> *hes) const;
private:
  size_t basis_num_, verts_num_;
  size_t nods_sz1_, nods_sz2_, nods_sz_;
  const double w_;
  const matd_t refs_;
};

class full_aug_anis_smooth_energy : public Functional<double>
{
public:
  typedef Eigen::Triplet<double> TPL;
  full_aug_anis_smooth_energy(const fine_mesh_patch &patch, const double w);
  size_t Nx() const;
  int Val(const double *x, double *val) const;
  int Gra(const double *x, double *gra) const;
  int Hes(const double *x, std::vector<Eigen::Triplet<double>> *hes) const;
private:
  size_t basis_num_, verts_num_;
  size_t nods_sz1_, nods_sz2_, nods_sz_;
  const double w_;
  Eigen::SparseMatrix<double> G_;
  const size_t rd_;
};

class full_aug_anis_mini_perturbation : public Functional<double>
{
public:
  typedef Eigen::Triplet<double> TPL;
  full_aug_anis_mini_perturbation(const fine_mesh_patch &patch, const double w);
  size_t Nx() const;
  int Val(const double *x, double *val) const;
  int Gra(const double *x, double *gra) const;
  int Hes(const double *x, std::vector<Eigen::Triplet<double>> *hes) const;
private:
  size_t basis_num_, verts_num_;
  size_t nods_sz1_, nods_sz2_, nods_sz_;
  const double w_;
  Eigen::SparseMatrix<double> P_;
  const size_t rd_;
  const fine_mesh_patch &patch_;
};

class full_anis_kronecker_cons : public Constraint<double>
{
public:
  typedef Eigen::Triplet<double> TPL;
  full_anis_kronecker_cons(const fine_mesh_patch &patch);
  size_t Nx() const;
  size_t Nf() const;
  int Val(const double *x, double *val) const;
  int Jac(const double *x, const size_t off, std::vector<TPL> *jac) const;
  int Hes(const double *x, const size_t off, std::vector<std::vector<TPL>> *hes) const;
private:
  size_t verts_num_, basis_num_, nods_sz1_;
  mati_t basis_in_strip_;
};

class full_aug_anis_kronecker_cons : public Constraint<double>
{
public:
  typedef Eigen::Triplet<double> TPL;
  full_aug_anis_kronecker_cons(const fine_mesh_patch &patch);
  size_t Nx() const;
  size_t Nf() const;
  int Val(const double *x, double *val) const;
  int Jac(const double *x, const size_t off, std::vector<TPL> *jac) const;
  int Hes(const double *x, const size_t off, std::vector<std::vector<TPL>> *hes) const;
private:
  size_t verts_num_, basis_num_, nods_sz1_;
  mati_t basis_in_strip_;
};

class full_aug_anis_zero_bnd_cons : public Constraint<double>
{
public:
  typedef Eigen::Triplet<double> TPL;
  full_aug_anis_zero_bnd_cons(const fine_mesh_patch &patch);
  size_t Nx() const;
  size_t Nf() const;
  int Val(const double *x, double *val) const;
  int Jac(const double *x, const size_t off, std::vector<TPL> *jac) const;
  int Hes(const double *x, const size_t off, std::vector<std::vector<TPL>> *hes) const;
private:
  size_t verts_num_, basis_num_, nods_sz1_;
  mati_t ind_;
  size_t cnt_zero_;
};

class full_anis_PoU_cons : public Constraint<double>
{
public:
  typedef Eigen::Triplet<double> TPL;
  full_anis_PoU_cons(const fine_mesh_patch &patch);
  size_t Nx() const;
  size_t Nf() const;
  int Val(const double *x, double *val) const;
  int Jac(const double *x, const size_t off, std::vector<TPL> *jac) const;
  int Hes(const double *x, const size_t off, std::vector<std::vector<TPL>> *hes) const;
private:
  size_t verts_num_, basis_num_, nods_sz1_;
};

class full_aug_anis_PoU_cons : public Constraint<double>
{
public:
  typedef Eigen::Triplet<double> TPL;
  full_aug_anis_PoU_cons(const fine_mesh_patch &patch);
  size_t Nx() const;
  size_t Nf() const;
  int Val(const double *x, double *val) const;
  int Jac(const double *x, const size_t off, std::vector<TPL> *jac) const;
  int Hes(const double *x, const size_t off, std::vector<std::vector<TPL>> *hes) const;
private:
  size_t verts_num_, basis_num_, nods_sz1_;
};

class full_aug_anis_diag_cons : public Constraint<double>
{
public:
  typedef Eigen::Triplet<double> TPL;
  full_aug_anis_diag_cons(const fine_mesh_patch &patch);
  size_t Nx() const;
  size_t Nf() const;
  int Val(const double *x, double *val) const;
  int Jac(const double *x, const size_t off, std::vector<TPL> *jac) const;
  int Hes(const double *x, const size_t off, std::vector<std::vector<TPL>> *hes) const;
private:
  size_t verts_num_, basis_num_, nods_sz1_;
  std::vector<std::pair<size_t, size_t>> index_;
};

class full_aug_anis_RI_cons : public Constraint<double>
{
public:
  typedef Eigen::Triplet<double> TPL;
  full_aug_anis_RI_cons(const fine_mesh_patch &patch);
  size_t Nx() const;
  size_t Nf() const;
  int Val(const double *x, double *val) const;
  int Jac(const double *x, const size_t off, std::vector<TPL> *jac) const;
  int Hes(const double *x, const size_t off, std::vector<std::vector<TPL>> *hes) const;
private:
  size_t verts_num_, basis_num_, nods_sz1_;
  const fine_mesh_patch &patch_;
  const matd_t &param_nods_;
  matd_t cross_mat_;
};

class full_anis_to_iso_cons : public Constraint<double>
{
public:
  /* To verify the equivalence to iso bases */
  typedef Eigen::Triplet<double> TPL;
  full_anis_to_iso_cons(const fine_mesh_patch &patch);
  size_t Nx() const;
  size_t Nf() const;
  int Val(const double *x, double *val) const;
  int Jac(const double *x, const size_t off, std::vector<TPL> *jac) const;
  int Hes(const double *x, const size_t off, std::vector<std::vector<TPL>> *hes) const;
private:
  size_t verts_num_, basis_num_, nods_sz1_;
};

// class kronecker_delta_energy : public Functional<double>
// {
// public:
//   typedef Eigen::Triplet<double> TPL;
//   kronecker_delta_energy(const fine_quad_patch &patch, const double w);
//   size_t Nx() const;
//   int Val(const double *x, double *val) const;
//   int Gra(const double *x, double *gra) const;
//   int Hes(const double *x, std::vector<TPL> *hes) const;
// private:
//   const double w_;
//   size_t vert_num_1d_, basis_num_1d_;
//   mati_t loc_;
//   Eigen::MatrixXd Id_;
// };

// class partition_to_unity_energy : public Functional<double>
// {
// public:
//   typedef Eigen::Triplet<double> TPL;
//   partition_to_unity_energy(const fine_quad_patch &patch, const double w);
//   size_t Nx() const;
//   int Val(const double *x, double *val) const;
//   int Gra(const double *x, double *gra) const;
//   int Hes(const double *x, std::vector<TPL> *hes) const;
// private:
//   const double w_;
//   size_t vert_num_1d_, basis_num_1d_;
// };

// class patch_interp_cons_2d : public Constraint<double>
// {
// public:
//   typedef Eigen::Triplet<double> TPL;
//   patch_interp_cons_2d(const fine_quad_patch &patch, const size_t mode_num, const int wei_option);
//   size_t Nx() const;
//   size_t Nf() const;
//   int Val(const double *x, double *val) const;
//   int Jac(const double *x, const size_t off, std::vector<TPL> *jac) const;
//   int Hes(const double *x, const size_t off, std::vector<std::vector<TPL>> *hes) const;
// private:
//   size_t mode_num_;     // # data set num
//   size_t vert_num_1d_, basis_num_1d_, basis_flat_dim_;
//   size_t nods_sz1_, nods_sz2_, nods_sz_;
//   const fine_quad_patch &patch_;
//   mati_t basis_strip_;
//   const double eps_;
//   Eigen::VectorXd wgts_;
// };

}

#endif
