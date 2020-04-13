#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>     // std::cout
#include <algorithm>    // std::random_shuffle
#include <vector>       // std::vector
#include <ctime>        // std::time
#include <cstdlib> 
using namespace Eigen;
using namespace std;

int main(int argc, char** argv){
  const size_t dof = 6;
  VectorXi per_vec = VectorXi::LinSpaced(dof, 0, dof - 1);
  std::random_shuffle(per_vec.data(), per_vec.data() + per_vec.size());
  
  cout<<"permutation vector " << endl<< per_vec << endl;
  PermutationMatrix<dof, dof, int> P = per_vec.asPermutation();

  MatrixXd A = MatrixXd::Random(dof, dof);
  A = A.transpose() * A;
  cout << "origin of A " << endl << A << endl;
  cout << "P * A " << endl << P * A << endl;

  SparseMatrix<double> SP_A(dof, dof);{
    vector<Triplet<double>> trips;
    for(size_t i = 0; i < dof; ++i){
      for(size_t j = 0; j < dof; ++j){
        trips.push_back(Triplet<double>(i, j, A(i, j)));
      }
    }
    SP_A.reserve(trips.size());
    SP_A.setFromTriplets(trips.begin(), trips.end());
  }
  cout <<"SP A " << endl<< MatrixXd(SP_A) << endl;
  SparseMatrix<double> SP_A_P1 = P.transpose() * SP_A * P;
  // SP_A.twistedBy(P).evalTo(SP_A_P1);b
  cout << "twisted by P" << endl <<  MatrixXd(SP_A_P1) <<endl;;
  // cout << "twisted by inverse P" << endl << SP_A.twistedBy(P.inverse()) <<endl;;

  return 0;
}

