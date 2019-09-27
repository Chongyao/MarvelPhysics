#ifndef MARVEL_IMPLICIT_EULER_GPU
#define MARVEL_IMPLICIT_EULER_GPU
#include "implicit_euler.h"
// #include "gpu_solve.h"
#include "pcg_gpu.h"
// #include "vec.h"
namespace marvel{
template<typename T, size_t dim_>
class newton_iter_gpu : public newton_iter<T, dim_>{
 public:
  using SMP_TYPE = typename dat_str_core<T, dim_>::SMP_TYPE;
  newton_iter_gpu(std::shared_ptr<dat_str_core<T, dim_>>& dat_str,
                  std::shared_ptr<Functional<T, dim_>>& energy,
                  const T time_step = 0.01, const size_t max_iter = 20, const T tol = 1e-4, const bool if_pre_compute_hes = false, const bool if_line_search = true, const bool if_hes_constant = false) : newton_iter<T, dim_>(dat_str, energy, time_step, max_iter, tol, if_pre_compute_hes, if_line_search, if_hes_constant){}

 // protected:
  bool if_first{true};
  
  virtual int linear_solver(const SMP_TYPE* A, const Eigen::Matrix<T, -1, 1>& b, Eigen::Matrix<T, -1, 1>& solution){
    cout << "in gpu pcg solver" << endl;
    __TIME_BEGIN__;
    if(if_first){
      readIAandJA(A->cols(), A->nonZeros(), A->outerIndexPtr(), A->innerIndexPtr());
      if_first = false;
    }

    PCG<T> solver;
    solver.function_pcg(A->cols(), A->nonZeros(), A->valuePtr(), b.data(), solution.data());
    __TIME_END__("[INFO] solve linear system cost");

    // for(size_t i = 0; i < num_rhs; ++i){
    //   #pragma omp parallel for
    //   for(size_t i = 0; i < num_rhs; ++i){
    //     solution(i) = solution_double[i];
    //   }      
    // }
    
    return 0;
  }
  


  
};


// template class newton_iter_gpu<float, 3>;
}
#endif
