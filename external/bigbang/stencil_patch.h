#ifndef STENCIL_PATCH_H
#define STENCIL_PATCH_H

#include <boost/property_tree/ptree.hpp>
#include <jtflib/mesh/mesh.h>
#include <bitset>
#include "stencil.h"

namespace bigbang {

struct fine_mesh_patch
{
  fine_mesh_patch(const boost::property_tree::ptree &pt) : pt_(pt) {}
  virtual ~fine_mesh_patch() {}
  virtual void extract_vert_strip() = 0;
  virtual void modal_analysis() = 0;
  virtual void optimize_basis(matd_t &basis) = 0;
  virtual void optimize_basis_modi_houman(matd_t &basis) = 0;
  virtual int write_optimized_basis(const char *filename, const matd_t &basis) = 0;
  virtual void sample_FEM_basis(matd_t &basisFEM) = 0;
  virtual void optimize_anis_basis(matd_t &basis);
  virtual void optimize_aug_anis_basis(matd_t &basis,
                                       const std::bitset<6> &flags,
                                       const matd_t *gl_har_disp=nullptr,
                                       const Eigen::MatrixXd *gl_eig_disp=nullptr,
                                       const Eigen::VectorXd *gl_eig_wgts=nullptr);
  
  mati_t cell_;
  matd_t nods_, REST_;
  matd_t lame_;

  //-> nodes of parametric domain
  matd_t param_nods_;
  
  mati_t new_to_orig_;
  mati_t global_to_local_;

  mati_t vert_strip_perm_;
  mati_t basis_pnt_perm_;

  size_t ONE_TO_MANY_;
  mati_t basis_pnt_;
  mati_t vert_strip_;

  Eigen::VectorXd freqs_;
  Eigen::MatrixXd modes_;
  
  Eigen::VectorXd bkp_freqs_;
  Eigen::MatrixXd bkp_modes_;

  const boost::property_tree::ptree &pt_;
  std::shared_ptr<Functional<double>> energy_;
  Eigen::SparseMatrix<double> K_, M_;

  matd_t FEM_basis_;
  mati_t FEM_indicator_;

  //-> prolongation matrix, responsible for
  //-> prolongate basis points to sample nodes
  std::shared_ptr<Eigen::SparseMatrix<double>> Pro_;
};

class straight_quad_strip_extractor;
class straight_cube_strip_extractor;

//===============================================================================

struct fine_quad_patch : fine_mesh_patch
{
  fine_quad_patch(const quad_stencil &dom, const mati_t &quad_H,
                  const mati_t &quad_h, const matd_t &nods_h, const matd_t &lame_h,
                  const bool to_be_refine, const boost::property_tree::ptree &pt);
  void extract_vert_strip();
  void modal_analysis();
  void optimize_basis(matd_t &basis);
  void optimize_basis_modi_houman(matd_t &basis);
  int write_optimized_basis(const char *filename, const matd_t &basisXY);
  void sample_FEM_basis(matd_t &basisFEM);
  void apply_local_LMA(const size_t eignum, const quad_stencil *ptr_sten,
                       Eigen::SparseMatrix<double> &K_h, Eigen::SparseMatrix<double> &M_h,
                       Eigen::SparseMatrix<double> &P, Eigen::VectorXd &eigs_h);
  std::shared_ptr<straight_quad_strip_extractor> strip_ext_;
  std::shared_ptr<jtf::mesh::edge2cell_adjacent> e2c_;
};

int get_prolongation_matrix(const quad_stencil *ptr_sten, const fine_mesh_patch *ptr_patch,
                            Eigen::SparseMatrix<double> &P);
int stencil_prolongate_patch(const matd_t &nods_H, const quad_stencil *ptr_sten,
                             fine_mesh_patch *ptr_patch,
                             mati_t &cell_h, matd_t &nods_h);
int find_shared_boundary_element(const fine_quad_patch *patch_x,
                                 const fine_quad_patch *patch_y,
                                 mati_t &shared_bnd_edges,      // in global sense
                                 mati_t &corres_faces);         // in local sense

class stencil_nonconforming_d2 : public Functional<double>
{
public:
  stencil_nonconforming_d2(const mati_t &quad_H, const matd_t &rest_H,
                           const std::vector<std::shared_ptr<quad_stencil>> &stencils,
                           const std::vector<std::shared_ptr<fine_mesh_patch>> &patches,
                           const double w);
  size_t Nx() const;
  int Val(const double *x, double *val) const;
  int Gra(const double *x, double *gra) const;
  int Hes(const double *x, std::vector<Eigen::Triplet<double>> *hes) const;
protected:
  const size_t dim_;
  const double w_;
  const matd_t &rest_H_;
  const std::vector<std::shared_ptr<quad_stencil>> &stencils_;
  const std::vector<std::shared_ptr<fine_mesh_patch>> &patches_;
  std::vector<std::pair<size_t, size_t>> adjc_sten_;
  std::vector<std::pair<std::shared_ptr<Eigen::SparseMatrix<double>>,
                        std::shared_ptr<Eigen::SparseMatrix<double>>>> S_;
  std::vector<std::shared_ptr<Eigen::VectorXd>> K_;
  std::vector<double> elem_bnd_len_;
};

//===============================================================================

struct fine_hexs_patch : fine_mesh_patch
{
  fine_hexs_patch(const hexs_stencil &dom, const mati_t &hexs_H,
                  const mati_t &hexs_h, const matd_t &nods_h, const matd_t &lame_h,
                  const bool to_be_refine, const boost::property_tree::ptree &pt);
  void extract_vert_strip();
  void modal_analysis();
  void optimize_basis(matd_t &basis);
  void optimize_basis_modi_houman(matd_t &basis);
  int write_optimized_basis(const char *filename, const matd_t &basisXYZ);
  void sample_FEM_basis(matd_t &basisFEM);
  void apply_local_LMA(const size_t eignum, const hexs_stencil *ptr_sten,
                       Eigen::SparseMatrix<double> &K_h, Eigen::SparseMatrix<double> &M_h,
                       Eigen::SparseMatrix<double> &P, Eigen::VectorXd &eigs_h);
  std::shared_ptr<straight_cube_strip_extractor> strip_ext_;
  std::shared_ptr<jtf::mesh::face2hex_adjacent> f2h_;
};

int get_prolongation_matrix(const hexs_stencil *ptr_sten, const fine_mesh_patch *ptr_patch,
                            Eigen::SparseMatrix<double> &P);
int stencil_prolongate_patch(const matd_t &nods_H, const hexs_stencil *ptr_sten,
                             fine_mesh_patch *ptr_patch,
                             mati_t &cell_h, matd_t &nods_h);
int find_shared_boundary_element(const fine_hexs_patch *patch_x,
                                 const fine_hexs_patch *patch_y,
                                 mati_t &shared_bnd_faces,
                                 mati_t &corres_elems);

class stencil_nonconforming_d3 : public Functional<double>
{
public:
  stencil_nonconforming_d3(const mati_t &hexs_H, const matd_t &rest_H,
                           const std::vector<std::shared_ptr<hexs_stencil>> &stencils,
                           const std::vector<std::shared_ptr<fine_mesh_patch>> &patches,
                           const double w);
  size_t Nx() const;
  int Val(const double *x, double *val) const;
  int Gra(const double *x, double *gra) const;
  int Hes(const double *x, std::vector<Eigen::Triplet<double>> *hes) const;
public:
  const size_t dim_;
  const double w_;
  const matd_t &rest_H_;
  const std::vector<std::shared_ptr<hexs_stencil>> &stencils_;
  const std::vector<std::shared_ptr<fine_mesh_patch>> &patches_;
  std::vector<std::pair<size_t, size_t>> adjc_sten_;
  mutable matd_t adjc_sten_E_;
  std::vector<std::pair<std::shared_ptr<Eigen::SparseMatrix<double>>,
                        std::shared_ptr<Eigen::SparseMatrix<double>>>> S_;
  std::vector<std::shared_ptr<mati_t>> shared_fine_nods_;
  std::vector<std::shared_ptr<Eigen::VectorXd>> K_;
  std::vector<std::shared_ptr<Eigen::VectorXd>> A_;
  //  std::vector<double> elem_bnd_len_;
};

}

#endif
