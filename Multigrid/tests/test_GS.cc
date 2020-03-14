#include "Multigrid/src/gauss_seidel.h"
#include <Eigen/Dense>
#include <iostream>
using namespace std;
using namespace Eigen;
using namespace marvel;

int main(int argc, char** argv){
  const size_t dim = stoi(argv[1]);
  MatrixXd A = MatrixXd::Random(dim, dim);
  A = (A.transpose() * A).eval();
  VectorXd b = VectorXd::Random(dim);

  cout <<"A is " <<endl<< A << endl << "b is" << endl << b << endl;

  VectorXd solution = A.lu().solve(b);
  cout << "residual of lu "<< (A * solution - b).norm() << endl;

  SparseMatrix<double> SP_A(dim, dim);{
    vector<Triplet<double>> trips;
    for(size_t i = 0; i < dim; ++i)
      for(size_t j = 0; j < dim; ++j){
        trips.push_back(Triplet<double>(i, j, A(i, j)));
      }
    SP_A.reserve(trips.size());
    SP_A.setFromTriplets(trips.begin(), trips.end());
  }
  
  VectorXd solu_gs = VectorXd::Zero(dim);
  gauss_seidel_solver(SP_A, b, solu_gs, stoi(argv[2]));
  cout << "residual of gs "<< (A * solu_gs - b).norm() << endl;
  
  return 0;
}


