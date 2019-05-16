#ifndef MARVEL_IMPLICIT_EULER
#define MARVEL_IMPLICIT_EULER
#include "def.h"
#include "data_str_core.h"
#include <memory>

namespace marvel{
template<typename T, size_t dim_>
class newton_iter{
 public:
  newton_iter(std::shared_ptr<dat_str_core<T, dim_>> data,
              std::shared_ptr<energy_t<T, dim_>> energy,
              const T time_step = 0.01, const size_t max_iter = 20, const T tol = 1e-4, const bool if_pre_compute_hes = false);
  int solve(T* x);
 private:
  const T time_step_;
  const size_t max_iter_;
  const T tol_;

  std::shared_ptr<dat_str_core<T, dim_>> data_;
  std::shared_ptr<energy_t<T, dim_>> energy_;
}
  
}
#endif
