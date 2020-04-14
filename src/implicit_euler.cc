#include "implicit_euler.h"
#include "line_search.h"
#include "config.h"
#include <fstream>
#include <Eigen/Core>
#include <typeinfo>
using namespace std;
using namespace Eigen;
namespace marvel{

template<typename T, size_t dim_>
newton_iter<T,dim_>::newton_iter(shared_ptr<dat_str_core<T, dim_>>& dat_str, shared_ptr<Functional<T, dim_>>& energy,
                                 const size_t max_iter, const T tol,
                                 const bool if_pre_compute_hes, const bool if_line_search,
                                  const bool if_hes_constant_)
    :max_iter_(max_iter), tol_(tol), dat_str_(dat_str), energy_(energy), if_line_search_(if_line_search), if_hes_constant_(if_hes_constant_), dof_(dat_str->get_dof()), cg(make_shared<ConjugateGradient<SparseMatrix<T>, Lower|Upper>>()){

  if(if_pre_compute_hes){
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
    dat_str_->setFromTriplets();
    const auto& sm1 = dat_str_->get_hes();
    cout<<"the number of nonzeros with comparison: \n"
        << (Eigen::Map<const Matrix<T, -1, 1>> (sm1.valuePtr(), sm1.nonZeros()).array() != 0).count()
        << endl;

    dat_str_->set_zero_after_pre_compute();
    __TIME_END__("[INFO] Pre_compute_hes");

  }
      
      
}

template<typename T, size_t dim_>
int newton_iter<T, dim_>::solve(T* x){
  Eigen::Map<Matrix<T, Dynamic, 1>>X(x, dim_ * dof_);
  Matrix<T, Dynamic, 1> res = Matrix<T, -1, 1>::Zero(dim_ * dof_);
  Matrix<T, -1, 1> solution(dim_ * dof_);
  T res_last_iter;
  
  Matrix<T, -1, 1> rhs;
  for(size_t newton_i = 0; newton_i < max_iter_; ++newton_i){
    cout << "[INFO]>>>>>>>>>>>>>>>newton iter is " << newton_i << endl;
    dat_str_->set_zero();
    energy_->Val(x, dat_str_);
    energy_->Gra(x, dat_str_);
    res = dat_str_->get_gra();
    // __TIME_BEGIN__
    if (!(if_hes_constant_ && has_hes_computed_))
      energy_->Hes(x, dat_str_);
    // __TIME_END__("[INFO]Assemble hessian");
        
    dat_str_->hes_compress();
    const T res_value = res.array().abs().sum();
    cout << "[INFO]: Newton res " << res_value << endl;
    cout << "[INFO]: ALL Energy: " << dat_str_->get_val() << endl;
    T eps;{
      if(typeid(res_value) == typeid(float))
        eps = 1e-2;
      else if(typeid(res_value) == typeid(double))
        eps = 1e-5;
      else
        throw runtime_error("run type");
    }

    
    if(res_value < eps){
      cout << endl;
      break;
    }

    if(newton_i > 1 && fabs(res_value / res_last_iter -1) < 1e-2 )
      break;

    res_last_iter = res_value;

    const  SMP_TYPE *A;
    if (!(if_hes_constant_ && has_hes_computed_)) {
      A = &dat_str_->get_hes();
      if (if_hes_constant_ && !has_hes_computed_)
        constant_hes_ = dat_str_->get_hes();
    }     
    else
      A = &constant_hes_;

    linear_solver(A, -res, solution);
    cout << "solution " << solution.array().abs().sum() << endl;    
    if(if_line_search_)
      X += line_search<T, dim_>(dat_str_->get_val(), static_cast<T>(res.dot(solution)), energy_, dat_str_, x, solution.data()) * solution;
    else
      X +=  solution;

  }

  has_hes_computed_ = true;
  return 0;
}

template<typename T, size_t dim_>
int newton_iter<T, dim_>::linear_solver(const SMP_TYPE* A, const Eigen::Matrix<T, -1, 1>& b, Eigen::Matrix<T, -1, 1>& solution){
  
  __TIME_BEGIN__
 if (!(if_hes_constant_ && has_hes_computed_)){
        cg->compute(*A);        
 }

  solution = cg->solve(b);
  __TIME_END__("[INFO]: Solve linear system by CG");
  return 0;
}
template class newton_iter<double, 3>;
template class newton_iter<float, 3>;
template class newton_iter<double, 1>;
template class newton_iter<float, 1>;

}//namespace
