#include "DEFINE_TYPE.h"
#define EIGEN_USE_BLAS
#include "basic_energy.h"
#include "implicit_euler.h"
#include "io.h"

#include "FEM/src/elas_energy.h"
#include "FEM/src/mass_matrix.h"
#include <iostream>

using namespace std;
using namespace Eigen;
using namespace marvel;


// using FLOAT_TYPE = double;

using TET_ELAS = BaseElas<FLOAT_TYPE, 3, 4, 1, 1, linear_csttt, basis_func, quadrature>;
using HEX_ELAS = BaseElas<FLOAT_TYPE, 3, 8, 1, 2, linear_csttt, basis_func, quadrature>;

int main(int argc, char** argv){
  Eigen::initParallel();
  std::cout.precision(17);
  const char* filename = argv[1];
  
  Matrix<FLOAT_TYPE, -1, -1> nods(1, 1);
  MatrixXi cells(1, 1);

  const string type = argv[4];
  if(type  == "tet")
    mesh_read_from_vtk<FLOAT_TYPE, 4>(filename, nods, cells);
  else if(type == "hex")
    mesh_read_from_vtk<FLOAT_TYPE, 8>(filename, nods, cells);
    
  

  const size_t num_nods = nods.cols();
  cout <<"V"<< nods.rows() << " " << nods.cols() << endl << "T " << cells.rows() << " "<< cells.cols() << endl;

  // {/scale
  //   nods.row(0) += Matrix<FLOAT_TYPE, 1, -1>::Ones(nods.cols()) * 10;
  //   nods.row(1) += Matrix<FLOAT_TYPE, 1, -1>::Ones(nods.cols()) * 10;
  //   nods.row(2) += Matrix<FLOAT_TYPE, 1, -1>::Ones(nods.cols()) * 0.5;
  //   nods *= 0.05;
  // }

  const string outdir = argv[3];
  
  //set mtr
  constexpr  FLOAT_TYPE rho = 0.5;
  constexpr  FLOAT_TYPE Young = 8000.0;
  constexpr  FLOAT_TYPE poi = 0.45;
  constexpr  FLOAT_TYPE gravity = 98;
  constexpr  FLOAT_TYPE dt = 0.1;
  const      FLOAT_TYPE w_pos = 1e4;
  const      size_t num_frame = 100;


  //read fixed points
  vector<size_t> cons(0);
  const char* cons_file_path = argv[2];
  read_fixed_verts_from_csv(cons_file_path, cons);
  cout << "constrint " << cons.size() << " points" << endl;
  

  //set collision
  // vector<shared_ptr<signed_dist_func<FLOAT_TYPE, 3>>>  objs(6);
  // Matrix<FLOAT_TYPE, 3, 1> plane_center;plane_center.setZero();
  // for(size_t i = 0; i < 2; ++i){
  //   plane_center += Matrix<FLOAT_TYPE, 3, 1>::Ones() * i;
  //   for(size_t j = 0; j < 3; ++j){
  //     Matrix<FLOAT_TYPE, 3, 1> plane_normal = Matrix<FLOAT_TYPE, 3, 1>::Zero();
  //     plane_normal(j) = pow(-1, i) * 1;
  //     objs[i * 3 + j] = make_shared<planeSDF<FLOAT_TYPE,3>>(plane_center.data(), plane_normal.data());
  //   }
  // }
  
  

  
  //calc mass vector
  Matrix<FLOAT_TYPE, -1, 1> mass_vec(nods.rows() * num_nods);
  calc_mass_vector<FLOAT_TYPE>(nods, cells, rho, mass_vec);

  cout << "build energy" << endl;
  shared_ptr<Matrix<FLOAT_TYPE, -1, -1>> init_points_ptr  = make_shared<Matrix<FLOAT_TYPE, -1, -1>>(Matrix<FLOAT_TYPE, -1, -1>::Zero(nods.rows(), nods.cols()));
  enum energy_type{ELAS, GRAV, KIN, POS};
  vector<shared_ptr<Functional<FLOAT_TYPE, 3>>> ebf(POS + 1);{
    if(type == "tet" )
      ebf[ELAS] = make_shared<TET_ELAS>(nods, cells, Young, poi);
    else if(type == "hex")
      ebf[ELAS] = make_shared<HEX_ELAS>(nods, cells, Young, poi);

    // ebf[ELAS] = nullptr;
    ebf[GRAV] = make_shared<gravity_energy<FLOAT_TYPE, 3>>(num_nods, 1, gravity, mass_vec, 'y');
    ebf[KIN] = make_shared<momentum<FLOAT_TYPE, 3>>(nods.data(), num_nods, mass_vec, dt);
    ebf[POS] = make_shared<position_constraint<FLOAT_TYPE, 3>>(nods.data(), num_nods, w_pos, cons);

    // ebf[POS] = make_shared<collision<FLOAT_TYPE, 3>>(nods.cols(), 1e5, 'x', 0.05, nods.cols(), init_points_ptr);
    // ebf[POS] = make_shared<geom_contact_energy<FLOAT_TYPE,3>>(objs, num_nods, 1e5);
    // ebf[POS] = nullptr;
    

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
  shared_ptr<dat_str_core<FLOAT_TYPE, 3>>  dat_str = make_shared<dat_str_core<FLOAT_TYPE, 3>>(num_nods);
  newton_iter<FLOAT_TYPE, 3> imp_euler(dat_str, energy, dt, 20, 1e-4, true, false, true);

  for(size_t f_id = 0; f_id < num_frame; ++f_id){
    cout << "[frame " << f_id << "]" << endl;
    imp_euler.solve(nods.data());
    dynamic_pointer_cast<momentum<FLOAT_TYPE, 3>>(ebf[KIN])->update_location_and_velocity(nods.data());

    const string filename = outdir  + "/frame_" + to_string(f_id) + ".vtk";
    tet_mesh_write_to_vtk<FLOAT_TYPE>(filename.c_str(), nods, cells);
  }
  return 0;
}
