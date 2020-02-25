#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/Eigenvalues>
#include <memory>
#include <iostream>

#include "Multigrid/src/Sparsification.h"
using namespace std;
using namespace Eigen;
using namespace marvel;
int main(){
  const size_t dim = 5;
  //construct A
  MatrixXd A = MatrixXd::Random(dim, dim);
  A = (A.transpose() * A).eval();
  for(size_t i = 0; i < dim; ++i){
    A(i, i) = 0;
    for(size_t j = i + 1; j < dim; ++j){
      if(A(i, j) > 0){
        A(i, j) *= -1;
        A(j, i) *= -1;
      }
    }
    for(size_t j = 0; j < dim; ++j){
      if(j != i)
        A(i, i) += A(i, j);
    }
    A(i, i) *= -1;
  }
  cout << A << endl;
  VectorXd b = VectorXd::Random(dim);
  MatrixXd E = VectorXd::Ones(dim).asDiagonal();
  MatrixXd L = A + E;
  const VectorXd x_gt = L.ldlt().solve(b);
  cout << "ground truth " <<endl <<x_gt << endl;
  cout << "L x - b " <<  (L * x_gt - b).norm()<< endl;

  //get eigen values
  VectorXd eig_vals = A.eigenvalues().real();
  sort(eig_vals.data(), eig_vals.data() + eig_vals.size());
  cout<<"eig values :\n" <<eig_vals << endl <<endl<<eig_vals.maxCoeff() / eig_vals(1) << endl<<endl;

  //=========================================test====================//
  Adjc_graph G(A);
  G.Sparsification();
  SparseMatrix<double> spm_A(A.rows(), A.cols());{
    vector<Triplet<double>> trips;
    G.build_mat_from_graph(trips);
    spm_A.setFromTriplets(trips.begin(), trips.end());
  }
  cout << MatrixXd(spm_A) << endl;;

  //get eigen values
  eig_vals = MatrixXd(spm_A).eigenvalues().real();
  sort(eig_vals.data(), eig_vals.data() + eig_vals.size());
  cout <<"eig values \n"<<eig_vals << endl <<endl<<eig_vals.maxCoeff() / eig_vals(1) << endl<<endl;
  getchar();

  //see solution
  L = MatrixXd(spm_A) + E;
  VectorXd x = L.ldlt().solve(b);
  cout << "L x- b" << (L * x - b).norm() << endl;

  cout << (x - x_gt) << endl;
  
  // // sparsification
  // cout << "=======================Sparsification=============================" <<endl;
  // auto sparsification = [&A](const size_t i, const size_t j)->double{
  //   const double w = A(i, j);
  //   A(i, i) += w;
  //   A(j, j) += w;
  //   A(i, j) = 0;
  //   A(j, i) = 0;
  //   return w;
  // };uyyyyyyyyyyyy
  // const double w = sparsification(0, 1);
  // cout << A << endl;
  // eig_vals = A.eigenvalues().real();
  // sort(eig_vals.data(), eig_vals.data() + eig_vals.size());
  // cout <<eig_vals << endl <<endl<<eig_vals.maxCoeff() / eig_vals(1) << endl<<endl;

  // {
  //   L = A + E;
  //   VectorXd x  = L.ldlt().solve(b);
  //     cout << " x - x_gt " << endl << x - x_gt <<endl;
  // }
  // //compensate
  // cout << "==============================Compensation========================="<<endl;
  // auto compensate = [&A](const double& w, const size_t i, const size_t j, const size_t k)->void{
  //   A(i, i) -= w;
  //   A(j, j) -= w;
  //   A(k, k) -= 2 * w;
  //   A(i, k) += w;
  //   A(k, i) += w;
  //   A(j, k) += w;
  //   A(k, j) += w;
  // };

  // compensate(w, 0, 1, 4);
  
  // cout << A << endl;
  // eig_vals = A.eigenvalues().real();
  // sort(eig_vals.data(), eig_vals.data() + eig_vals.size());
  // cout <<eig_vals << endl <<endl<<eig_vals.maxCoeff() / eig_vals(1) << endl<<endl;
  
  // {
  //   L = A + E;
  //   VectorXd x  = L.ldlt().solve(b);
  //   cout << " x - x_gt " << endl << x - x_gt <<endl;
  // }


  
  
  
  return 0;
}	
