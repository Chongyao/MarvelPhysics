#ifndef MARVEL_PCG
#define MARVEL_PCG
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <functional>
#include "implicit_euler.h"
#include "config.h"
namespace marvel{
using precond_type = std::function<Eigen::VectorXd(const Eigen::VectorXd&)>;

class PCG{
 public:
  using SPM = Eigen::SparseMatrix<double, Eigen::RowMajor>;
  PCG(const SPM& A, const size_t& max_itrs, const double& tol, const precond_type& M = nullptr);
  Eigen::VectorXd solve(const Eigen::VectorXd& b);
  private:
  const size_t dim_;
  const SPM A_;
  const size_t max_itrs_;
  const double tol_;
  const precond_type M_;
};

template<typename T, size_t dim_>
class newton_iter_myPCG : public newton_iter<T, dim_>{
 public:
  using SMP_TYPE = typename dat_str_core<T, dim_>::SMP_TYPE;
  newton_iter_myPCG(std::shared_ptr<dat_str_core<T, dim_>>& dat_str,
                    std::shared_ptr<Functional<T, dim_>>& energy,
                    const T time_step = 0.01, const size_t max_iter = 20,
                    const T tol = 1e-4, const bool if_pre_compute_hes = false,
                    const bool if_line_search = true,
                    const bool if_hes_constant = false,
                    const precond_type& PD = nullptr) : newton_iter<T, dim_>(dat_str, energy, time_step, max_iter, tol, if_pre_compute_hes, if_line_search, if_hes_constant), M_(PD){}

  precond_type M_;
  int set_preconditioner(const precond_type M){
    M_ = M;
    return 0;
  }
  bool if_first_{true};
  int linear_solver(const SMP_TYPE* A, const Eigen::Matrix<T, -1, 1>& b, Eigen::Matrix<T, -1, 1>& solution){
    static std::shared_ptr<PCG> pcg;
    using base = newton_iter<T, dim_>;
    if(base::if_hes_constant_ && if_first_){
      pcg = std::make_shared<PCG>(*A, A->rows() * 2, base::tol_, M_);
      if_first_  = false;
    }
    else if(!base::if_hes_constant_)
      pcg = std::make_shared<PCG>(*A, A->rows() * 2, base::tol_, M_);
      
    __TIME_BEGIN__;
    solution = pcg->solve(b);
    __TIME_END__("[INFO] solve linear system cost");

    return 0;
  }
};


}
#endif
