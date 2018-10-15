#ifndef FEM_STENCIL_H
#define FEM_STENCIL_H

#include <zjucad/matrix/matrix.h>
#include <boost/property_tree/ptree.hpp>

#include "def.h"

namespace bigbang {

using mati_t=zjucad::matrix::matrix<size_t>;
using matd_t=zjucad::matrix::matrix<double>;
using namespace std;
// All the basis functions follows the order
//
// N-1.. 2N-1 ..3N-1
//  : ... : ... :
//  : ... : ... :
//  2 .. N+2 .. :
//  1 .. N+1 .. :
//  0 ... N ... 2N
//
typedef void (*FEM_basis_func)(double*, const double*);
typedef void (*dis_FEM_basis_func)(double*, const double*, const matd_t &basisXYZ);

void full_quad_dis_shape_func_val(double *val, const double *eps, const matd_t &basisXY);
void full_quad_dis_shape_func_jac(double *jac, const double *eps, const matd_t &basisXY);
void full_hexs_dis_shape_func_val(double *val, const double *eps, const matd_t &basisXYZ);
void full_hexs_dis_shape_func_jac(double *jac, const double *eps, const matd_t &basisXYZ);

///=========== QUAD STENCIL ==========///
class quad_stencil
{
public:
  quad_stencil(const mati_t &adjc_elem, const std::vector<size_t> &faces);
  mati_t adjc_elem_;
  mati_t faces_;
  size_t adjc_vert_num_;
  FEM_basis_func fem_sf_jac_, rest_sf_jac_;

  virtual int set_full_discrete_shape_func(const matd_t &basisXY, const bool isotropic=true,
                                           const bool augmented=false);
  dis_FEM_basis_func dis_fem_sf_jac_;
  matd_t basisXY_;
  bool isotropic_, augmented_;
  
  virtual int set_mtrs_and_qrs(const size_t per_face_to_many,
                               const mati_t &quad_h, const matd_t &nods_h,
                               const matd_t &lame_h);
  size_t qr_num_;
  matd_t qrs_, qrw_;
  matd_t lame_on_qrs_;

  virtual int set_def_grad_oper(const matd_t &adjc_nods);
  matd_t adjc_H_invDmH_; // matrix: [adjc_vert_num_, 2*#quadratures]
  matd_t adjc_df_op_;    // matrix: [2*adjc_vert_num_, 4*#quadratures]
  matd_t dets_;          // matrix: [#quadratures]
  matd_t adjc_rest_;
  matd_t DmHs_;          // matrix: [2 x 2#quadratures]
  
  virtual int query_local_frame(const matd_t &curr_nods, const size_t qrID, matd_t &R);
  virtual int get_coro_bases_df(const matd_t &curr_nods, const matd_t &R, const size_t qrID,
                                matd_t &dF, matd_t &G);
  matd_t dNude_, dNxde_;

  virtual void apply_local_LMA(const mati_t &quad_H, const matd_t &nods_H, const size_t eignum,
                               const Eigen::SparseMatrix<double> &K_h,
                               const Eigen::SparseMatrix<double> &M_h,
                               const Eigen::SparseMatrix<double> &P,
                               const boost::property_tree::ptree &pt,
                               Eigen::VectorXd &eigs_H,
                               Eigen::VectorXd &rqs_H);

};

bool is_adjc_quad_stencil(const std::shared_ptr<quad_stencil> &x, const std::shared_ptr<quad_stencil> &y,
                          mati_t &shared_verts);

///========== HEX STENCIL ==========///
class hexs_stencil
{
public:
  hexs_stencil(const mati_t &adjc_elem, const std::vector<size_t> &faces);
  mati_t adjc_elem_;
  mati_t faces_;
  size_t adjc_vert_num_;
  FEM_basis_func fem_sf_jac_, rest_sf_jac_;

  virtual int set_full_discrete_shape_func(const matd_t &basisXYZ, const bool isotropic=true,
                                           const bool augmented=false);
  dis_FEM_basis_func dis_fem_sf_jac_;
  matd_t basisXYZ_;
  bool isotropic_, augmented_;

  /* qnum_per_cube: how many 1D quadrature number in EACH FINE cube */
  virtual int set_mtrs_and_qrs(const size_t per_face_to_many, const size_t qnum_per_cube,
                               const mati_t &hexs_h, const matd_t &nods_h,
                               const matd_t &lame_h);
  size_t qr_num_;
  matd_t qrs_, qrw_;
  matd_t lame_on_qrs_;

  virtual int set_def_grad_oper(const matd_t &adjc_nods);
  matd_t adjc_H_invDmH_; // matrix: [adjc_vert_num_, 3*#quadratures]
  matd_t adjc_df_op_;    // matrix: [3*adjc_vert_num_, 9*#quadratures]
  matd_t dets_;          // matrix: [#quadratures]
  matd_t adjc_rest_;
  matd_t DmHs_;
  matd_t HinvDmH0_;      // matrix: [adjc_vert_num, 3]

  virtual int query_local_frame(const matd_t &curr_nods, const size_t qrID, matd_t &R);
  virtual int get_coro_bases_df(const matd_t &curr_nods, const matd_t &R, const size_t qrID,
                                matd_t &dF, matd_t *G);
  matd_t dNude_, dNxde_;

  virtual void apply_local_LMA(const mati_t &hexs_H, const matd_t &nods_H, const size_t eignum,
                               const Eigen::SparseMatrix<double> &K_h,
                               const Eigen::SparseMatrix<double> &M_h,
                               const Eigen::SparseMatrix<double> &P,
                               const boost::property_tree::ptree &pt,
                               Eigen::VectorXd &eigs_H,
                               Eigen::VectorXd &rqs_H);
};

bool is_adjc_quad_stencil(const std::shared_ptr<hexs_stencil> &x, const std::shared_ptr<hexs_stencil> &y,
                          mati_t &shared_verts);

///========== STENCIL ENERGIES ==========///
class fem_stencil_energy_d2 : public Functional<double>
{
public:
  enum Material {
    LINEAR,
    STVK,
    COROTATED,
    NEOHOOKEAN
  };
  fem_stencil_energy_d2(const std::vector<std::shared_ptr<quad_stencil>> &elem_H,
                        const mati_t &quad_H, const matd_t &nods_H,
                        const Material type,
                        const mati_t &quad_h, const matd_t &nods_h, const matd_t &face_lame_h,
                        const double w);
  size_t Nx() const;
  int Val(const double *x, double *val) const;
  int Gra(const double *x, double *gra) const;
  int Hes(const double *x, std::vector<Eigen::Triplet<double>> *hes) const;
protected:
  const double w_;
  const size_t dim_;
  const std::vector<std::shared_ptr<quad_stencil>> &elem_;
  const Material type_;
};

class fem_stencil_energy_d3 : public Functional<double>
{
public:
  enum Material {
    LINEAR,
    STVK,
    COROTATED,
    NEOHOOKEAN,
    STANEO,
    BOWERNEO
  };
  fem_stencil_energy_d3(const std::vector<std::shared_ptr<hexs_stencil>> &elem_H,
                                             const mati_t &hexs_H, const matd_t &nods_H,
                                             const Material type,
                                             const mati_t &hexs_h, const matd_t &nods_h,
                                             const matd_t &face_lame_h,
                                             const size_t qnum_per_cube,
                                             const double w);
  size_t Nx() const;
  int Val(const double *x, double *val) const;
  int Gra(const double *x, double *gra) const;
  int Hes(const double *x, std::vector<Eigen::Triplet<double>> *hes) const;

    


public:
  const double w_;
  const size_t dim_;
  const std::vector<std::shared_ptr<hexs_stencil>> &elem_;
  const Material type_;
  mutable matd_t stencilE_;


};

std::shared_ptr<Functional<double>>
build_stencil_energy(const std::string &mtr, const std::vector<std::shared_ptr<quad_stencil>> &stencils,
                     const mati_t &quad_H, const matd_t &nods_H,
                     const mati_t &quad_h, const matd_t &nods_h, const matd_t &face_lame_h);

std::shared_ptr<Functional<double>>
build_stencil_energy(const std::string &mtr, const std::vector<std::shared_ptr<hexs_stencil>> &stencils,
                     const size_t qnum_per_cube,
                     const mati_t &hexs_H, const matd_t &nods_H,
                     const mati_t &hexs_h, const matd_t &nods_h, const matd_t &face_lame_h);

}

#endif
