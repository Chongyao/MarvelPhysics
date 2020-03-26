#ifndef MARVEL_GAUSS_SEIDEL
#define MARVEL_GAUSS_SEIDEL
#include <Eigen/Sparse>
#include <memory>

namespace marvel{
int gauss_seidel_solver(const Eigen::SparseMatrix<double, Eigen::RowMajor>& A, const Eigen::VectorXd& b, Eigen::VectorXd& solution, const size_t itrs);

class Jacobi_like_solver{
 public:
  ~Jacobi_like_solver(){}
  virtual int solve(const Eigen::VectorXd& b, Eigen::VectorXd& solution)const = 0;
};
enum class solver_type {GS, WJ, PCG};


class Gauss_seidel : public Jacobi_like_solver{
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
#if 1
class Weighted_Jacobi: public Jacobi_like_solver{
 public:
  using SPM = Eigen::SparseMatrix<double, Eigen::RowMajor>;
  Weighted_Jacobi(const SPM& A, const size_t& max_itr, const double& error, const double& w = 2.0/3.0);
  int solve(const Eigen::VectorXd& b, Eigen::VectorXd& solution) const;
  int set_optimal_w();
  
 private:
  const Eigen::VectorXd dig_vals_;
  const SPM U_plus_L_;
  const SPM A_;
  
  double w_;
  const double error_;
  const size_t max_itr_;
  
  
};
#endif



}
#endif
