#include "DEFINE_TYPE.h"
#define EIGEN_USE_BLAS
// #include <omp.h>

#include "basic_energy.h"
#include "implicit_euler.h"
#include "io.h"
#include "extract_surface.imp"

#include "FEM/src/elas_energy.h"
#include "FEM/src/mass_matrix.h"
#include <iostream>

using namespace std;
using namespace Eigen;
using namespace marvel;

// using FLOAT_TYPE = double;

int main(int argc, char** argv){
  if (argc < 3)
  {
    cerr << "please run ./ext path_to_model type(tet, hex..)" << endl;
    return -1;
  } 
  const char* filename = argv[1];
  const string type = argv[2];

  Matrix<FLOAT_TYPE, -1, -1> nods(1, 1);
  MatrixXi cells(1, 1);

  if(type  == "tet")
    mesh_read_from_vtk<FLOAT_TYPE, 4>(filename, nods, cells);
  else if(type == "hex")
    mesh_read_from_vtk<FLOAT_TYPE, 8>(filename, nods, cells);
  
  MatrixXi surface;
  extract_surface(nods, cells, surface, type);
  MatrixXi s = Map<MatrixXi>(surface.data(), 3, surface.size()/3);
  cerr << "tri num:" << s.cols() << endl;
  tri_mesh_write_to_vtk("surface.vtk", nods, s);
  return 0;
}
