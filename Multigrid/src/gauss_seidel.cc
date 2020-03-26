#include "config.h"

#include "gauss_seidel.h"
#include <iostream>
#include <Eigen/Eigenvalues>
#include <set>
#include "search_eigenvalues.h"
namespace marvel{
using namespace std;
using namespace Eigen;


int gauss_seidel_solver(const Eigen::SparseMatrix<double, Eigen::RowMajor>& A, const Eigen::VectorXd& b, Eigen::VectorXd& solution, const size_t itrs){
  assert(A.rows() == A.cols()&&A.rows() == b.size());
  const size_t dof = b.size();
  
  const auto& first_ids = A.outerIndexPtr();
  const auto& col_ids = A.innerIndexPtr();
  const auto& vals = A.valuePtr();
  
  std::vector<size_t> dig_ids(dof);{
    #pragma omp parallel for
    for(size_t i = 0; i < dof; ++i){
      const size_t start = first_ids[i], end = first_ids[i + 1];
      for(size_t j = start; j < end; ++j){
        if(col_ids[j] == i){
          dig_ids[i] =  j;
          break;
        }
      }
    }
  }
  auto solution_now_ptr = std::make_shared<Eigen::VectorXd>(solution);
  
  auto solution_next_ptr = std::make_shared<Eigen::VectorXd>(b);
  VectorXd dig_vals(dof);{
    #pragma omp parallel for
    for(size_t i = 0; i < dof; ++i)
      dig_vals[i] = vals[dig_ids[i]];
  }
 
  size_t itr_cnt = 0;
  do{
    auto& solution_next = *solution_next_ptr;
    auto& solution_now = *solution_now_ptr;
    solution_next = b.array() / dig_vals.array();
    //loop for upper triangle part of A
    #pragma omp parallel for
    for(size_t i = 0; i < dof - 1; ++i){
      const size_t start = dig_ids[i] + 1, end = first_ids[i + 1];

      double sum = 0;
      for(size_t j = start; j < end; ++j){
        sum += solution_now(col_ids[j]) * vals[j];
      }
      solution_next(i) -= sum / vals[dig_ids[i]];
    }

    //loop for lower triangle part of A
    for(size_t i = 1; i < dof; ++i){
      const size_t start = first_ids[i], end = dig_ids[i];

      double sum = 0;
      for(size_t j = start; j < end; ++j){
        sum += solution_next(col_ids[j]) * vals[j];
      }
      solution_next(i) -= sum / vals[dig_ids[i]];
    }

    swap(solution_next, solution_now);
    ++itr_cnt;
  }while(itr_cnt < itrs);
  solution = *solution_now_ptr;
  
  return 0;
}

//===========================Gauss seidel=============================//

Gauss_seidel::Gauss_seidel(const SPM& A, const size_t& max_itr):A_(A), max_itr_(max_itr), dof_(A.rows()){
  assert(A.rows() == A.cols());
  const auto& first_ids = A_.outerIndexPtr();
  const auto& col_ids = A_.innerIndexPtr();
  const auto& vals = A_.valuePtr();

  dig_ids_.resize(dof_);
  #pragma omp parallel for
  for(size_t i = 0; i < dof_; ++i){
    const size_t start = first_ids[i], end = first_ids[i + 1];
    for(size_t j = start; j < end; ++j){
      if(col_ids[j] == i){
        dig_ids_[i] =  j;
        break;
      }
    }
  }
  dig_vals_.resize(dof_);
  #pragma omp parallel for
  for(size_t i = 0; i < dof_; ++i)
    dig_vals_[i] = vals[dig_ids_[i]];
}


  
int Gauss_seidel::solve(const Eigen::VectorXd& b, Eigen::VectorXd& solution) const{
  const auto& first_ids = A_.outerIndexPtr();
  const auto& col_ids = A_.innerIndexPtr();
  const auto& vals = A_.valuePtr();

  auto solution_now_ptr = std::make_shared<Eigen::VectorXd>(solution);
  auto solution_next_ptr = std::make_shared<Eigen::VectorXd>(b);

  size_t itr_cnt = 0;
  double time_upper = 0, time_lower = 0;
  do{
    auto& solution_next = *solution_next_ptr;
    auto& solution_now = *solution_now_ptr;
    solution_next = b.array() / dig_vals_.array();
    //loop for upper triangle part of A
    __TIME_BEGIN__;
    #pragma omp parallel for
    for(size_t i = 0; i < dof_ - 1; ++i){
      const size_t start = dig_ids_[i] + 1, end = first_ids[i + 1];

      double sum = 0;
      for(size_t j = start; j < end; ++j)
        sum += solution_now(col_ids[j]) * (vals[j]);

      solution_next(i) -= sum / dig_vals_[i];
    }
    time_upper +=__TIME_END__("upper part ", false);
    __TIME_BEGIN__;
    //loop for lower triangle part of A
    for(size_t i = 1; i < dof_; ++i){
      const size_t start = first_ids[i], end = dig_ids_[i];

      double sum = 0;
      for(size_t j = start; j < end; ++j)
        sum += solution_next(col_ids[j]) * (vals[j]);
      solution_next(i) -= sum / dig_vals_[i];
    }
    time_lower += __TIME_END__("lower part ", false);
    swap(solution_next_ptr, solution_now_ptr);
    ++itr_cnt;
  }while(itr_cnt < max_itr_);
  cout << "time upper is " << time_upper << " seconds.\n"
       << "time lower is " << time_lower << " seconds.\n";
  solution = *solution_now_ptr;
  
  return 0;
}
#if 1
//=========================Weighted_Jacobi========================//
Weighted_Jacobi::Weighted_Jacobi(const SPM& A, const size_t& max_itr, const double& error, const double& w)
    :w_(w), error_(error), max_itr_(max_itr),dig_vals_(A.diagonal()),U_plus_L_(A.triangularView<StrictlyUpper>() + A.triangularView<StrictlyLower>()), A_(A){
  double max_eigval, min_eigval;
  find_max_min_eigenvalues(A, max_eigval, min_eigval);
  *(const_cast<double*>(&w_)) = 2.0 / (max_eigval + min_eigval);
}

int Weighted_Jacobi::solve(const VectorXd& b, VectorXd& solution) const{
   for(size_t i = 0; i < max_itr_; ++i)
     solution.array() = w_ * (b - U_plus_L_ * solution).array() / dig_vals_.array() + (1 - w_) * solution.array();


  return 0;
}
#endif



}
