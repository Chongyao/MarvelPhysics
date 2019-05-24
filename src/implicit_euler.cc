#include "implicit_euler.h"
#include "line_search.h"
#include "config.h"
#include <Eigen/Core>
using namespace std;
using namespace Eigen;
namespace marvel{

template<typename T, size_t dim_>
newton_iter<T,dim_>::newton_iter(shared_ptr<dat_str_core<T, dim_>>& dat_str, shared_ptr<Functional<T, dim_>>& energy,
                                 const T time_step, const size_t max_iter, const T tol,
                                 const bool if_pre_compute_hes, const bool if_line_search)
    :time_step_(time_step), max_iter_(max_iter), tol_(tol), dat_str_(dat_str), energy_(energy), if_line_search_(if_line_search), dof_(dat_str->get_dof()){


  Matrix<T, Dynamic, 1> random_x(dim_ * dof_);{
  #pragma omp parallel for
    for(size_t i = 0; i < dim_ * dof_; ++i){
      random_x(i) = i * 4.5 + i * i;
    }
  }
  dat_str_->set_zero();
  __TIME_BEGIN__
  energy_->Val(random_x.data(), dat_str_);
  energy_->Gra(random_x.data(), dat_str_);
  energy_->Hes(random_x.data(), dat_str_);

  const auto& sm1 = dat_str_->get_hes();
  cout<<"the number of nonzeros with comparison: \n"
      << (Map<const Matrix<T, -1, 1>> (sm1.valuePtr(), sm1.nonZeros()).array() != 0).count()
      << endl;
  dat_str_->set_zero_after_pre_compute();
  __TIME_END__("[INFO] Pre_compute_hes");
      
}

template<typename T, size_t dim_>
int newton_iter<T, dim_>::solve(T* x){
  
  Map<Matrix<T, Dynamic, 1>>X(x, dim_ * dof_);
  const Matrix<T, Dynamic, 1>& res = dat_str_->get_gra();
  SimplicialLLT<SparseMatrix<T>> llt;
  Matrix<T, -1, 1> solution;

  
  for(size_t newton_i = 0; newton_i < max_iter_; ++newton_i){
    cout << "[INFO]>>>>>>>>>>>>newton iter is " << newton_i << endl;
    dat_str_->set_zero();
    energy_->Val(x, dat_str_);
    energy_->Gra(x, dat_str_);
   
    __TIME_BEGIN__
    energy_->Hes(x, dat_str_);
    __TIME_END__("[INFO]Assemble hessian");
        
    dat_str_->hes_compress();
    
    const T res_value = res.array().square().sum();
    cout << "[INFO]: Newton res " << res_value << endl;
    cout << "[INFO]: ALL Energy: " << dat_str_->get_val() << endl;
    if(res_value < 1e-4){
      cout << endl;
      break;
    }

    //Solve linear system By Cholesky
    __TIME_BEGIN__
    llt.compute(dat_str_->get_hes());
    size_t time = 1;
    while(llt.info() != Eigen::Success){
      cout <<"lltinfo "<< llt.info() << endl;
      dat_str_->hes_add_diag(time);
      llt.compute(dat_str_->get_hes());
      time *= 2;
    }

    // const auto  solution = llt.solve(-dat_str_->get_gra());
    solution = llt.solve(-res).eval();
    __TIME_END__("[INFO]: Solve linear system by Cholesky");

    
    if(if_line_search_)
      X += line_search<T, dim_>(dat_str_->get_val(), static_cast<T>(res.dot(solution)), energy_, dat_str_, x, solution.data()) * solution;
    else
      X +=  solution;

  }
  return 0;
}
}//namespace
