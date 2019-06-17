#ifndef MARVEL_IMPLICIT_EULER
#define MARVEL_IMPLICIT_EULER
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
              const T time_step = 0.01, const size_t max_iter = 20, const T tol = 1e-4, const bool if_pre_compute_hes = false, const bool if_line_search = true);
  int solve(T* x);
 protected:
  template<typename Derived, typename OtherDerived>
  int linear_solver(const SMP_TYPE* A, const Eigen::MatrixBase<Derived>& b, Eigen::MatrixBase<OtherDerived>& solution);
  
 protected:
  const T time_step_;
  const size_t max_iter_;
  const T tol_;
  const size_t dof_;

  const bool if_line_search_;

  std::shared_ptr<dat_str_core<T, dim_>> dat_str_;
  std::shared_ptr<Functional<T, dim_>> energy_;
};

template class newton_iter<double, 3>;
template class newton_iter<float, 3>;
  
}

#include "implicit_euler.imp"

#endif
