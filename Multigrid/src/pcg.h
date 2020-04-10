#ifndef MARVEL_PCG
#define MARVEL_PCG
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <functional>

namespace marvel{
using precond_type = std::function<Eigen::VectorXd(Eigen::VectorXd&)>;

class PCG{
 public:
  using SPM = Eigen::SparseMatrix<double, Eigen::RowMajor>;
  PCG(const SPM& A, const size_t& max_itrs, const double& tol, const precond_type& M = nullptr);

  Eigen::VectorXd solve(const Eigen::VectorXd& b, const Eigen::VectorXd* x0 = nullptr);
  private:
  const size_t dim_;
  const SPM A_;
  const size_t max_itrs_;
  const double tol_;
  const precond_type M_;
};
}
#endif
