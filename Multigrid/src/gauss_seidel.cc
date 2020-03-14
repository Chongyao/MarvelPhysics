#include "gauss_seidel.h"
#include <iostream>
#include <Eigen/Eigenvalues>
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


#if 0
//===========================Gauss seidel=============================//
Gauss_seidel::Gauss_seidel(const SPM& A, const double* b, const size_t& max_iter, const double& error, const double& SOR_weight) : max_iter_(max_iter), error_(error), w_(SOR_weight), A_(A), b_(Map<const VectorXd>(b, A.outerSize(), 1)), dim_(A.rows()){}

int Gauss_seidel::solve(const double* init_x, double* post_x) const {
  const Map<const VectorXd> init(init_x, dim_);
  Map<VectorXd> solution(post_x, dim_);

  VectorXd prev = init;
  VectorXd post(dim_), diag_ele(dim_);
  cout<<" Before gauss-seidel residual " << (A_ * init - b_).norm()<<endl;
  for(size_t iter_time = 0; iter_time < max_iter_; ++iter_time){
    const double res = (A_ * prev - b_).norm();
    if(fabs(res) < error_)
      break;

    post.setZero();    
    //loop for upper triangle part
    #pragma omp parallel for
    for(size_t row_id = 0; row_id < dim_; ++row_id){
      double upper_part = 0;
      for (SPM::InnerIterator it(A_,row_id); it; ++it){
        if(it.col() > row_id)
          upper_part += it.value() * prev(it.col());
        else if(it.col() == row_id)
          diag_ele(row_id) = it.value();
      }
      
      post(row_id) += (1 - w_) * prev(row_id) + w_ / diag_ele(row_id) * (b_(row_id) - upper_part);
    }
    
    for(size_t row_id = 0; row_id < dim_; ++row_id){
      double low_part = 0;
      for (SPM::InnerIterator it(A_, row_id); it; ++it){
        if(it.col() >= row_id)
          break;
        low_part += it.value() * post(it.col());
      }
      post(row_id) -= w_  / diag_ele(row_id) * low_part;
    }
    
    prev = post;
  }
  cout<<" After gauss-seidel residual " << (A_ * solution - b_).norm()<<endl;
  solution = post;
  return 0;
}

void Gauss_seidel::update_problem(const SPM& A, const double* b){
  A_ = A;
  b_ = Map<const VectorXd>(b, dim_);
  return;
}
#endif


//===========================Gauss seidel=============================//

}
