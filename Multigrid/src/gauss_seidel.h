#ifndef MARVEL_GAUSS_SEIDEL
#define MARVEL_GAUSS_SEIDEL
#include <Eigen/Sparse>
#include <memory>

namespace marvel{

int gauss_seidel(const Eigen::SparseMatrix<double, Eigen::RowMajor>& A, const Eigen::VectorXd& b, Eigen::VectorXd& solution, const size_t itrs){
  assert(A.rows() == A.cols() == b.size());
  const size_t dof = b.size();
  solution.resize(dof);
  solution.setZero();

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
  
  size_t itr_cnt = 0;
  do{
    auto& solution_next = *solution_next_ptr;
    auto& solution_now = *solution_now_ptr;
    solution_next = b;

    //loop for upper triangle part of A
    #pragma omp parallel for
    for(size_t i = 0; i < dof - 1; ++i){
      const size_t start = dig_ids[i] + 1, end = first_ids[i + 1];
      double sum = 0;
      for(size_t j = start; j < end; ++j){
        sum += solution_now(col_ids[j]) * vals[j];
      }
      solution_next(i) -= sum;
    }

    //loop for lower triangle part of A
    for(size_t i = 1; i < dof; ++i){
      const size_t start = first_ids[i], end = dig_ids[i];
      double sum = 0;
      for(size_t j = start; j < end; ++j){
        sum += solution_next(col_ids[j]) * vals[j];
      }
      solution_next(i) -= sum;
    }

    #pragma omp parallel for
    for(size_t i = 0; i < dof; ++i){
      solution_next(i) /= vals[dig_ids[i]];
    }
    swap(solution_next_ptr, solution_now_ptr);
    
    const double res = (A * solution_now - b).norm() / b.norm();
    if(res < 1e-6)
      break;
    
    ++itr_cnt;
  }while(itr_cnt < itrs);
  solution = *solution_now_ptr;
  
  return 0;
}
}
#endif
