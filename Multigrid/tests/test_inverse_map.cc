#include "Multigrid/src/inverse_isoparametric_map.h"
using namespace std;
using namespace Eigen;
using namespace marvel;

int main(int argc, char** argv){
  MatrixXd nods(3, 8);
  nods <<
      -1, -1, -1,
      1, -1, -1,
      1, 1, -1,
      -1, 1, -1,
      -1, -1, 1,
      1, -1, 1,
      1, 1, 1,
      1, -1, 1;
  VectorXi cell(8) = VectorXi::LinSpaced(8, 0, 7);
  cout << cell << endl;

  Vector3d gl_coor;
  gl_coor<< -1,-1,-1;

  Vector3d iso_coor;
  inverse_isoparametric_hex(gl_coor, nods.data(), iso_coor);
  cout << iso_coor << endl;

  return 0;
}
