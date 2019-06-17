#ifndef MARVEL_IMPLICIT_EULER_GPU
#define MARVEL_IMPLICIT_EULER_GPU
#include "implicit_euler.h"
#include ""
namespace marvel{
template<typename DOUBLE, size_t dim_>
class newton_iter_gpu : public newton_iter<DOUBLE, dim_>{
  newton_iter_gpu(std::shared_ptr<dat_str_core<T, dim_>>& dat_str,
                  std::shared_ptr<Functional<T, dim_>>& energy,
                  const T time_step = 0.01, const size_t max_iter = 20, const T tol = 1e-4, const bool if_pre_compute_hes = false, const bool if_line_search = true) : newton_iter<DOUBLE, dim_>(dat_str, energy, time_step, max_iter, tol, if_pre_compute_hes, if_line_search){}
  
  template<typename Derived, typename OtherDerived>
  int linear_solver(const Eigen::SparseMatrix<T>* A, const Eigen::MatrixBase<Derived>& b, Eigen::MatrixBase<OtherDerived>& solution) override{
    
  }



  
}
}
#endif
