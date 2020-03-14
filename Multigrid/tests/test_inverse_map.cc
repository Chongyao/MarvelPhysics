#include "Multigrid/src/inverse_isoparametric_map.h"
#include <iostream>
using namespace std;
using namespace Eigen;
using namespace marvel;

int main(int argc, char** argv){
  MatrixXd nods(8, 3);
  nods <<
      -1, -1, -1,
      1, -1, -1,
      1, 1, -1,
      -1, 1, -1,
      -1, -1, 1,
      1, -1, 1,
      1, 1, 1,
      -1, 1, 1;
  nods.transposeInPlace();
  nods *= 2;
  VectorXi cell = VectorXi::LinSpaced(8, 0, 7);

  Vector3d gl_coor;
  gl_coor<< 0.5, -1.2 , 1;

  Vector3d iso_coor = inverse_isoparametric_hex(gl_coor, nods.data());
  
  cout << iso_coor << endl;

  return 0;
}
