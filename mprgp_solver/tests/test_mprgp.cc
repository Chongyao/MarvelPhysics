#define USE_STOPWATCH_HELPER_FLAG
#define INFO_PRINT_FLAG
#define ERROR_PRINT_FLAG
#define WARN_PRINT_FLAG
#define ASSERT_EXT_FLAG
#define DEBUG_PRINT_FLAG
#include "mprgp_solver.h"
#include <iostream>

#include <Eigen/Dense>
#include <Eigen/Sparse>

using namespace std;
using namespace chaos::mprgp;
using namespace chaos::utils;
using namespace Eigen;



int main(int argc, char *argv[])
{
  // log_utils::lu_printer::get_instance().set_os(log_utils::DEBUG, "debug_info");
  const int n = 100;
  const MatrixXd M = MatrixXd::Random(n,n) + MatrixXd::Identity(n,n)*3.0f;
  const MatrixXd MtM = M.transpose()*M;
  const SparseMatrix<double> A = Dense2Sparse(MtM);

  const MatrixXd JM = MatrixXd::Identity(n,n);
  const SparseMatrix<double> J = Dense2Sparse(JM);

  const VectorXd B = VectorXd::Random(n);
  const VectorXd c = VectorXd::Random(n);
  VectorXd x = VectorXd::Random(n);

  {
    STW_START("general");
    cout << "general: " << endl;
    VectorXd y = x;
    const int rlst_code = MPRGPGeneralConSolver<double>::solve(MtM,B,J,c,y);
    cout << rlst_code << endl;
    // cout << "y: " << y << endl;
    STW_END("general", "general solver");
  }
  {
    STW_START("decoupled");
    cout << "decoupled: " << endl;
    VectorXd y = x;
    const int rlst_code = MPRGPDecoupleConSolver<double>::solve(A,B,J,c,y);
    cout << rlst_code << endl;
    // cout << "y: " << y << endl;
    STW_END("decoupled", "decoupled solver");
  }
  {
    STW_START("direct");
    cout << "direct: " << endl;
    VectorXd y = x;
    solveQP(MtM, B, J, c, y);
    // cout << "y: " << y << endl;
    STW_END("direct", "direct solver");
  }
  return 0;
}
