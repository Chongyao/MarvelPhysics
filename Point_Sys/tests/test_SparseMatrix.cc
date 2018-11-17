#include <vector>
#include <iostream>
#include <Eigen/Sparse>

using namespace std;
using namespace Eigen;

int main(){
  SparseMatrix<double> Spmat(20,20);
  vector<Triplet<double>> TripletList;
  TripletList.push_back(Triplet<double>(2,3,4));
  TripletList.push_back(Triplet<double>(6,3,4));
  TripletList.push_back(Triplet<double>(6,3,4));
  Spmat.setFromTriplets(TripletList.begin(), TripletList.end());

  cout << Spmat;
  Spmat.setFromTriplets(TripletList.begin(), TripletList.end());
  cout << Spmat;
    
}

