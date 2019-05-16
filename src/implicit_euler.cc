#include "implicit_euler.h"
#include "line_search.h"
#include "config.h"
#include <Eigen/Core>
using namespace std;
using namespace std;
namespace marvel{

template<typename T, size_t dim_>
newton_iter<T,dim_>::newton_iter(shared_ptr<dat_str_core<T, dim_>>& data, shared_ptr<energy_t<T, dim_>>& energy,
                                 const T time_step, const size_t max_iter, const T tol,
                                 const bool if_pre_compute_hes, const bool if_line_search)
    :time_step_(time_step), max_iter_(max_iter), tol_(tol), data_(data), energy_(energy), if_line_search_(if_line_search), dof_(data->get_dof()){


  Matrix<T, Dynamic, 1> random_x(dim_ * dof_);{
  #pragma omp parallel for
    for(size_t i = 0; i < dim_ * dof; ++i){
      random_x(i) = i * 4.5 + i * i;
    }
  }
  data->set_zero();
  energy_->Val(random_x.data(), data_);
  energy_->Gra(random_x.data(), data_);
  energy_->Hes(random_x.data(), data_);

  cout<<"the number of nonzeros with comparison: \n"
      << (Map<VectorXd> (sm1.valuePtr(), sm1.nonZeros()).array() != 0).count()
      << endl;
  data->set_zero_after_pre_compute();
  
}

template<typename T, size_t dim_>
int newton_iter<T, dim_>::solve(T* x){
  
  Map<Matrix<T, Dynamic, 1>>X(x, dim_ * dof_);
  Map<const VectorXd>res(data->get_gra().data(), 3 * dim);
  SimplicialLLT<SparseMatrix<double>> llt;
  auto solution = Matrix<T, Dynamic, 1>::Zero(dim_ * dof_);
  
  for(size_t newton_i = 0; newton_i < max_iter_; ++newton_i){
    cout << "[INFO]>>>>>>>>>>>>newton iter is " << newton_i << endl;
    data->set_zero();
    // data->hes_reserve(nnzs);
    energy->Val(x, data_);
    energy->Gra(x, data);
    
    __TIME_BEGIN__
    energy->Hes(displace_plus.data(), data);
    __TIME_END__("[INFO]Assemble hessian")
        
    data->hes_compress();
    
    const T res_value = res.array().square().sum();
    cout << "[INFO]: Newton res " <<std::setprecision(9)<< res_value << endl;
    cout << "[INFO]: ALL Energy: " << data->get_val() << endl;
    if(res_value < 1e-4){
      cout << endl;
      break;
    }

    //Solve linear system By Cholesky
    __TIME_BEGIN__
    llt.compute(dat_str->get_hes());
    size_t time = 1;
    while(llt.info() != Eigen::Success){
      cout <<"lltinfo "<< llt.info() << endl;
      dat_str->hes_add_diag(time);
      llt.compute(dat_str->get_hes());
      time *= 2;
    }
    solution = llt.solve(-res);
    __TIME_END__("[INFO]: Solve linear system by Cholesky")

    if(if_line_search_)
      X += line_search<T, dim_>(data->get_val(), res.dot(solution), energy, data, X, solution) * solution;
    else
      X +=  solution;
  }
  return 0;
}

