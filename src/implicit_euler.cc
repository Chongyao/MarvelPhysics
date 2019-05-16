#include "implicit_euler.h"
#include "line_search.h"

#include <Eigen/Core>
using namespace std;
using namespace std;
namespace marvel{

template<typename T, size_t dim_>
newton_iter<T,dim_>::newton_iter(shared_ptr<dat_str_core<T, dim_>>& data, shared_ptr<energy_t<T, dim_>>& energy,
                                 const T time_step, const size_t max_iter,
                                 const T tol, const bool if_pre_compute_hes)
    :time_step_(time_step), max_iter_(max_iter), tol_(tol), data_(data), energy_(energy){

  size_t dof = data->get_dof();
  Matrix<T, Dynamic, 1> random_x(dim_ * dof);{
  #pragma omp parallel for
    for(size_t i = 0; i < dim_ * dof; ++i){
      random_x(i) = i * 4.5 + i * i;
    }
  }
  dat_str->set_zero();
  energy_->Val(random_x.data(), dat_str_);
  energy_->Gra(random_x.data(), dat_str_);
  energy_->Hes(random_x.data(), dat_str_);

  cout<<"the number of nonzeros with comparison: \n"
      << (Map<VectorXd> (sm1.valuePtr(), sm1.nonZeros()).array() != 0).count()
      << endl;
  dat_str->set_zero_after_pre_compute();
  
}

template<typename T, size_t dim_>
int newton_iter<T, dim_>::solve(T* x){
  for(size_t newton_i = 0; newton_i < max_iter_; ++newton_i){
    cout << "[INFO]: newton iter is " << newton_i << endl;
    dat_str->set_zero();
    // dat_str->hes_reserve(nnzs);
    energy->Val(displace_plus.data(), dat_str);
    energy->Gra(displace_plus.data(), dat_str);
    auto start = system_clock::now();
    energy->Hes(displace_plus.data(), dat_str);
    auto end = system_clock::now();
    auto duration = duration_cast<microseconds>(end - start);
    cout <<  "assemble hessian cost" 
         << double(duration.count()) * microseconds::period::num / microseconds::period::den 
         << "seconds" << endl;
    dat_str->hes_compress();
    


    
  }

}


}

