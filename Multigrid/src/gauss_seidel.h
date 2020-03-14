#ifndef MARVEL_GAUSS_SEIDEL
#define MARVEL_GAUSS_SEIDEL
#include <Eigen/Sparse>
#include <memory>

namespace marvel{
int gauss_seidel_solver(const Eigen::SparseMatrix<double, Eigen::RowMajor>& A, const Eigen::VectorXd& b, Eigen::VectorXd& solution, const size_t itrs);
#if 0
class Gauss_seidel{
 public:
  using SPM = Eigen::SparseMatrix<double, Eigen::RowMajor>;
  Gauss_seidel(const SPM& A, const double* b, const size_t& max_iter, const double& error, const double& SOR_weight);
  int solve(const double* init_x, double* post_x) const;
  void update_problem(const SPM& A, const double* b);

 private:
  SPM A_;
  Eigen::VectorXd b_;
  const size_t max_iter_;
  const double error_;
  //SOR weight
  const double w_;
  const size_t dim_;
  
};
#endif



}
#endif
