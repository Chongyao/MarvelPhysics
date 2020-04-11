#include "Multigrid/src/pcg.h"
#include "Multigrid/src/hsc.h"
using namespace Eigen;
using namespace std;
using namespace marvel;
int main(int argc, char** argv){
  Matrix2d A;
  A << 4, 1, 1, 3;
  cout << "A "  << endl<< A << endl;

  SPM L = A.sparseView();
  Vector2d b;
  b << 1, 2;
  cout << "b " << endl << b << endl;
  PCG pcg(L, 2, 1e-14);
  VectorXd u = pcg.solve(b);
  cout << "u "<< endl<< u << endl;
  return 0;
}
