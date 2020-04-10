#ifndef MARVEL_MULTIGRID
#define MARVEL_MULTIGRID

#include<vector>
#include<memory>
#include<Eigen/Dense>
#include<Eigen/Sparse>
#include<data_str_core.h>
#include"gauss_seidel.h"
namespace marvel{

template<typename T>

using VS = std::vector<std::shared_ptr<T>>;
using SPM = Eigen::SparseMatrix<double, Eigen::RowMajor>;

using MatrixXd = Eigen::MatrixXd;
using MatrixXi = Eigen::MatrixXi;

struct transfer {
  transfer(const SPM& I, const SPM&R);
  SPM I_;
  SPM R_;
};


transfer get_transfer(const MatrixXd& nods_H, const MatrixXi& cells_H, const MatrixXd& nods_h, const MatrixXi& cells_h); 


struct layer{
  layer(const SPM& A, const solver_type& type, const size_t itrs);
  int solve();
  Eigen::VectorXd get_residual()const;
  
  const SPM A_;
  Eigen::VectorXd rhs_;
  Eigen::VectorXd u_;
  // const bool if_direct_;
  const solver_type type_;
  const size_t itrs_;


  std::shared_ptr<Jacobi_like_solver> GS_solver_;

};

class multigrid_process{
 public:
  multigrid_process(const std::vector<int>& process, VS<layer>& layers, const VS<transfer>& transfers);
  int execute(double* solution);
  
  
 protected:
  int relax(const size_t layer_id);
  int restrict(const size_t layer_id);
  int correct(const size_t layer_id);
  
  std::vector<int> process_;
  VS<layer> layers_;
  VS<transfer> transfers_;
  
  const size_t dof_;
  const SPM A_h_;
  const Eigen::VectorXd b_h_;
  
  // size_t layer_id_{0};xo

  
};

}
#endif
