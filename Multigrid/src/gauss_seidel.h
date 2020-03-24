#ifndef MARVEL_GAUSS_SEIDEL
#define MARVEL_GAUSS_SEIDEL
#include <Eigen/Sparse>
#include <memory>

namespace marvel{
int gauss_seidel_solver(const Eigen::SparseMatrix<double, Eigen::RowMajor>& A, const Eigen::VectorXd& b, Eigen::VectorXd& solution, const size_t itrs);

class Gauss_seidel{
 public:
  using SPM = Eigen::SparseMatrix<double, Eigen::RowMajor>;
  Gauss_seidel(const SPM& A, const size_t& max_itr);
  
  int solve(const Eigen::VectorXd& b, Eigen::VectorXd& solution) const;

 private:
  const SPM A_;
  const size_t max_itr_;
  const size_t dof_;

  std::vector<size_t> dig_ids_;
  Eigen::VectorXd dig_vals_;
  
  
  Eigen::VectorXd b_;
  
};




}
#endif
