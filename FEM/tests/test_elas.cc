
#include "DEFINE_TYPE.h"
#define EIGEN_USE_BLAS
// #define EIGEN_USE_LAPACKE
#include "basic_energy.h"
#include "implicit_euler.h"
#include "implicit_euler_gpu.h"
#include "io.h"

#include "FEM/src/elas_energy.h"
#include "FEM/src/mass_matrix.h"

#include <iostream>

using namespace std;
using namespace Eigen;
using namespace marvel;

// using FLOAT_TYPE = double;
using TET_ELAS = BaseElas<FLOAT_TYPE, 3, 4, 1, 1, linear_csttt, basis_func, quadrature>;
int main(int argc, char** argv){
  Eigen::initParallel();
  std::cout.precision(17);
  const char* filename = argv[1];
  
  Matrix<FLOAT_TYPE, -1, -1> nods(1, 1);
  MatrixXi tets(1, 1);
  tet_mesh_read_from_vtk<FLOAT_TYPE>(filename, nods, tets);
  const size_t num_nods = nods.cols();
  cout <<"V"<< nods.rows() << " " << nods.cols() << endl << "T " << tets.rows() << " "<< tets.cols() << endl;

  {//scale
    nods.row(0) += Matrix<FLOAT_TYPE, 1, -1>::Ones(nods.cols()) * 10;
    nods.row(1) += Matrix<FLOAT_TYPE, 1, -1>::Ones(nods.cols()) * 10;
    nods.row(2) += Matrix<FLOAT_TYPE, 1, -1>::Ones(nods.cols()) * 0.5;
    nods *= 0.05;
  }

  const string outdir = argv[3];
  
  //set mtr
  constexpr  FLOAT_TYPE rho = 100;
  constexpr  FLOAT_TYPE Young = 10000.0;
  constexpr  FLOAT_TYPE poi = 0.45;
  constexpr  FLOAT_TYPE gravity = 9.8;
  constexpr  FLOAT_TYPE dt = 0.01;
  const      FLOAT_TYPE w_pos = 1e4;
  const      size_t num_frame = 100;

  //read fixed points
  vector<size_t> cons(0);
  const char* cons_file_path = argv[2];
  // if ( boost::filesystem::exists(cons_file_path) )
  read_fixed_verts_from_csv(cons_file_path, cons);
  cout << "constrint " << cons.size() << " points" << endl;


  
  
  //calc mass vector
  Matrix<FLOAT_TYPE, -1, 1> mass_vec(nods.rows() * num_nods);

  calc_mass_vector<FLOAT_TYPE>(nods, tets, rho, mass_vec);
  cout <<"total mass is "<< mass_vec.sum() << endl;
  cout << "build energy" << endl;
  shared_ptr<Matrix<FLOAT_TYPE, -1, -1>> init_points_ptr  = make_shared<Matrix<FLOAT_TYPE, -1, -1>>(Matrix<FLOAT_TYPE, -1, -1>::Zero(nods.rows(), nods.cols()));
  enum energy_type{ELAS, GRAV, KIN, POS};
  vector<shared_ptr<Functional<FLOAT_TYPE, 3>>> ebf(POS + 1);{
    ebf[ELAS] = make_shared<TET_ELAS>(nods, tets, Young, poi);
    // ebf[ELAS] = nullptr;
    ebf[GRAV] = make_shared<gravity_energy<FLOAT_TYPE, 3>>(num_nods, 1, gravity, mass_vec, 'x');
    ebf[KIN] = make_shared<momentum<FLOAT_TYPE, 3>>(nods.data(), num_nods, mass_vec, dt);
    // ebf[POS] = make_shared<position_constraint<FLOAT_TYPE, 3>>(nods.data(), num_nods, w_pos, cons);
    // ebf[POS] = nullptr;
    ebf[POS] = make_shared<collision<FLOAT_TYPE, 3>>(nods.cols(), 1e5, 'x', 0.05, nods.cols(), init_points_ptr);
    
    }
  cout << "assemble energy" << endl;


  
  
  
  shared_ptr<Functional<FLOAT_TYPE, 3>> energy;
  try {
    energy = make_shared<energy_t<FLOAT_TYPE, 3>>(ebf);

  } catch ( std::exception &e ) {
    cerr << e.what() << endl;
    exit(EXIT_FAILURE);
  }

  //Sovle
  const string filename_tmp = outdir  + "/frame_origin.vtk";
  // tet_mesh_write_to_vtk<FLOAT_TYPE>(filename_tmp.c_str(), nods, tets);
  shared_ptr<dat_str_core<FLOAT_TYPE, 3>>  dat_str = make_shared<dat_str_core<FLOAT_TYPE, 3>>(num_nods);
  newton_iter<FLOAT_TYPE, 3> imp_euler(dat_str, energy, dt, 20, 1e-4, true, false);
  


  
  for(size_t f_id = 0; f_id < num_frame; ++f_id){
    cout << "[frame " << f_id << "]" << endl;
    imp_euler.solve(nods.data());

    dynamic_pointer_cast<momentum<FLOAT_TYPE, 3>>(ebf[KIN])->update_location_and_velocity(nods.data());

    const string filename = outdir  + "/frame_" + to_string(f_id) + ".vtk";
    tet_mesh_write_to_vtk<FLOAT_TYPE>(filename.c_str(), nods, tets);
  }

  


  



  
  
  return 0;
}
