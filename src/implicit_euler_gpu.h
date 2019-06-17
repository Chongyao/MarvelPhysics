#ifndef MARVEL_IMPLICIT_EULER_GPU
#define MARVEL_IMPLICIT_EULER_GPU
#include "implicit_euler.h"
#include "gpu_solve.h"
#include "vec.h"
namespace marvel{
template<typename T, size_t dim_>
class newton_iter_gpu : public newton_iter<T, dim_>{
 public:
  using SMP_TYPE = typename dat_str_core<T, dim_>::SMP_TYPE;
  newton_iter_gpu(std::shared_ptr<dat_str_core<T, dim_>>& dat_str,
                  std::shared_ptr<Functional<T, dim_>>& energy,
                  const T time_step = 0.01, const size_t max_iter = 20, const T tol = 1e-4, const bool if_pre_compute_hes = false, const bool if_line_search = true) : newton_iter<T, dim_>(dat_str, energy, time_step, max_iter, tol, if_pre_compute_hes, if_line_search){}

 // protected:

  virtual int linear_solver(const SMP_TYPE* A, const Eigen::Matrix<T, -1, 1>& b, Eigen::Matrix<T, -1, 1>& solution){
    cout << "in gpu " << endl;
    // cout << A->valuePtr() << endl;
    int num_rowptr = A->outerSize() + 1;
    vector<int> rowptr(num_rowptr);{
      #pragma omp parallel
      for(size_t i = 0; i < rowptr.size(); ++i){
        rowptr[i] = A->outerIndexPtr()[i] + 1;
      }
    }
    
     int num_colptr = A->nonZeros();
    vector<int> colptr(num_colptr);{
      for(size_t i = 0; i < colptr.size(); ++i){
        colptr[i] = A->innerIndexPtr()[i] + 1;
      }
    }
    int num_rhs = b.rows();

    gpucg_solve_(&rowptr[0], &num_rowptr, &colptr[0], &num_colptr, A->valuePtr(), &num_colptr, b.data(), &num_rhs, solution.data());
    return 0;
  }
  


  
};

template class newton_iter_gpu<FLOAT_TYPE, 3>;
// template class newton_iter_gpu<float, 3>;
}
#endif
