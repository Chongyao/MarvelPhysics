#include "io.h"
#include "FEM/src/elas_energy.h"
#include <iostream>
using namespace std;
using namespace Eigen;
using namespace marvel;

int main(int argc, char** argv){
  const char* filename = argv[1];
  cout << filename << endl;
  
  MatrixXd nods(1, 1);
  MatrixXi tets(1, 1);
  tet_mesh_read_from_vtk(filename, nods, tets);
  cout <<"V"<< nods.rows() << " " << nods.cols() << endl << "T " << tets.rows() << " "<< tets.cols() << endl;

  cout << tets.col(tets.cols() -1) << endl;
  
  return 0;
}
