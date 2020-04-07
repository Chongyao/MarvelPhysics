#include "DEFINE_TYPE.h"
#define EIGEN_USE_BLAS
#include "basic_energy.h"
#include "implicit_euler.h"
#include "io.h"

#include "FEM/src/elas_energy.h"
#include "FEM/src/mass_matrix.h"
#include "FEM/src/poisson.h"
#include <iostream>

using namespace std;
using namespace Eigen;
using namespace marvel;

int main(int argc, char** argv){
  Eigen::initParallel();
  // std::cout.precision(17);
  const char* filename = argv[1];
  
  Matrix<FLOAT_TYPE, -1, -1> nods(1, 1);
  MatrixXi cells(1, 1);


  const string type = argv[4];
  mesh_read_from_vtk<FLOAT_TYPE, 8>(filename, nods, cells);
    
  const size_t num_nods = nods.cols();
  cout <<"V"<< nods.rows() << " " << nods.cols() << endl << "T " << cells.rows() << " "<< cells.cols() << endl;

  Matrix<FLOAT_TYPE, -1, 1> u(num_nods);

  const string outdir = argv[3];
  
  //set mtr
  constexpr  FLOAT_TYPE k = 10;
  const      FLOAT_TYPE w_pos = 1e6;

  Matrix<FLOAT_TYPE, 1, -1> mtr = Matrix<FLOAT_TYPE, 1, -1>::Ones(cells.cols()) * k;

  //read fixed points
  vector<size_t> cons(0);
  const char* cons_file_path = argv[2];
  read_fixed_verts_from_csv(cons_file_path, cons);
  cout << "constrint " << cons.size() << " points" << endl;
  Matrix<FLOAT_TYPE, -1, 1>  u_cons = Matrix<FLOAT_TYPE, -1, 1>::Random(num_nods);
  
  cout << "build energy" << endl;
  enum energy_type{POISSON, BD};
  vector<shared_ptr<Functional<FLOAT_TYPE, 1>>> ebf(BD + 1);{
    ebf[POISSON] = make_shared<HEX_POISSON>(nods, cells, mtr);
    ebf[BD] = make_shared<position_constraint<FLOAT_TYPE, 1>>(u_cons.data(), num_nods, w_pos, cons);

    }
  cout << "assemble energy" << endl;
  
  shared_ptr<Functional<FLOAT_TYPE, 1>> energy;
  try {
    energy = make_shared<energy_t<FLOAT_TYPE, 1>>(ebf);

  } catch ( std::exception &e ) {
    cerr << e.what() << endl;
    exit(EXIT_FAILURE);
  }

  //Sovle

  const string filename_tmp = outdir  + "/frame_origin.vtk";
  shared_ptr<dat_str_core<FLOAT_TYPE, 1>>  dat_str = make_shared<dat_str_core<FLOAT_TYPE, 1>>(num_nods);
  newton_iter<FLOAT_TYPE, 1> imp_euler(dat_str, energy, 10000, 20, 1e-4, true, false, true);

  for(size_t f_id = 0; f_id < 1; ++f_id){
    cout << "[frame " << f_id << "]" << endl;
    imp_euler.solve(u.data());
    
    const string filename = outdir  + "/frame_" + to_string(f_id) + ".vtk";
    mesh_write_to_vtk<FLOAT_TYPE, 8>(filename.c_str(), nods, cells);
  }
  return 0;
}



