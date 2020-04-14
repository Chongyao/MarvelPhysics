#ifndef MARVEL_IMPLICIT_EULER
#define MARVEL_IMPLICIT_EULER
#include "DEFINE_TYPE.h"
#include "def.h"
#include "data_str_core.h"
#include <memory>

namespace marvel{

template<typename T, size_t dim_>
class newton_iter{
 public:
  using SMP_TYPE = typename dat_str_core<T, dim_>::SMP_TYPE;
  newton_iter(std::shared_ptr<dat_str_core<T, dim_>>& dat_str,
              std::shared_ptr<Functional<T, dim_>>& energy,
              const T time_step = 0.01, const size_t max_iter = 20, const T tol = 1e-4, 
              const bool if_pre_compute_hes = false, const bool if_line_search = true, const bool if_hes_constant = false);
  virtual int solve(T* x);
  // public:

  virtual int linear_solver(const SMP_TYPE* A, const Eigen::Matrix<T, -1, 1>& b, Eigen::Matrix<T, -1, 1>& solution) ;
  
 protected:
  const T time_step_;
  const size_t max_iter_;
  const T tol_;
  const size_t dof_;

  const bool if_line_search_;
  const bool if_hes_constant_{false};
  bool has_hes_computed_{false};
  std::shared_ptr<Eigen::ConjugateGradient<Eigen::SparseMatrix<T>, Eigen::Lower|Eigen::Upper>> cg;
  SMP_TYPE constant_hes_;

  std::shared_ptr<dat_str_core<T, dim_>> dat_str_;
  std::shared_ptr<Functional<T, dim_>> energy_;
};



}

#endif
