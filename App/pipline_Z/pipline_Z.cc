#include "DEFINE_TYPE.h"
#define EIGEN_USE_BLAS
// #include <omp.h>

#include "basic_energy.h"
#include "implicit_euler.h"
#include "io.h"

#include "FEM/src/elas_energy.h"
#include "FEM/src/mass_matrix.h"
#include "mprgp_solver/src/collision_proxy.h"
#include "mprgp_solver/src/collision_plane.h"
#include <iostream>
#include "constrained_newton.h"

using namespace std;
using namespace Eigen;
using namespace marvel;
using namespace chaos::collision;
// using FLOAT_TYPE = double;

using TET_ELAS = BaseElas<FLOAT_TYPE, 3, 4, 1, 1, linear_csttt, basis_func, quadrature>;
using HEX_ELAS = BaseElas<FLOAT_TYPE, 3, 8, 1, 2, linear_csttt, basis_func, quadrature>;

int main(int argc, char** argv){
  Eigen::initParallel();
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
  
  const string outdir = argv[3];

  //set mtr
  constexpr  FLOAT_TYPE rho = 20;
  constexpr  FLOAT_TYPE Young = 2000.0;
  constexpr  FLOAT_TYPE poi = 0.3;
  constexpr  FLOAT_TYPE gravity = 98;
  constexpr  FLOAT_TYPE dt = 0.01;
  const      FLOAT_TYPE w_pos = 1e6;
  const      size_t num_frame = 100;

  //read fixed points
  vector<size_t> cons(0);
  const char *cons_file_path = argv[2];
  read_fixed_verts_from_csv(cons_file_path, cons);
  cout << "constrint " << cons.size() << " points" << endl;

  //set collision
  vector<vector<size_t>> cells_4_coll(1, vector<size_t>(cells.size()));{
    copy(cells.data(), cells.data() + cells.size(), cells_4_coll[0].begin());
  }
  vector<vector<FLOAT_TYPE>> nods_4_coll(1, vector<FLOAT_TYPE>(nods.size()));{
    copy(nods.data(), nods.data() + nods.size(), nods_4_coll[0].begin());
  }
  
  shared_ptr<collision_proxy<FLOAT_TYPE>> coll = make_shared<collision_zplane<FLOAT_TYPE>>() ;
  coll->init(cells_4_coll, nods_4_coll);


  
  //calc mass vector
  Matrix<FLOAT_TYPE, -1, 1> mass_vec(nods.rows() * num_nods);
  // calc_mass_vector<FLOAT_TYPE>(nods, cells, rho, mass_vec);
  if(type == "tet")
    mass_calculator<FLOAT_TYPE, 3, 4, 1, 1, basis_func, quadrature>(nods, cells, rho, mass_vec);
  else if (type == "hex")
    mass_calculator<FLOAT_TYPE, 3, 8, 1, 2, basis_func, quadrature>(nods, cells, rho, mass_vec);

  cout << "build energy" << endl;
  shared_ptr<Matrix<FLOAT_TYPE, -1, -1>> init_points_ptr  = make_shared<Matrix<FLOAT_TYPE, -1, -1>>(Matrix<FLOAT_TYPE, -1, -1>::Zero(nods.rows(), nods.cols()));
  enum energy_type{ELAS, GRAV, KIN, POS};
  vector<shared_ptr<Functional<FLOAT_TYPE, 3>>> ebf(POS + 1);{
    if(type == "tet" )
      ebf[ELAS] = make_shared<TET_ELAS>(nods, cells, Young, poi);
    else if(type == "hex")
      ebf[ELAS] = make_shared<HEX_ELAS>(nods, cells, Young, poi);
    ebf[GRAV] = make_shared<gravity_energy<FLOAT_TYPE, 3>>(num_nods, 1, gravity, mass_vec, 'y');
    ebf[KIN] = make_shared<momentum<FLOAT_TYPE, 3>>(nods.data(), num_nods, mass_vec, dt);
    ebf[POS] = make_shared<position_constraint<FLOAT_TYPE, 3>>(nods.data(), num_nods, w_pos, cons);
    }

  cout << "assemble energy" << endl;

  shared_ptr<Functional<FLOAT_TYPE, 3>> energy;
  try
  {
    energy = make_shared<energy_t<FLOAT_TYPE, 3>>(ebf);
  }
  catch (std::exception &e)
  {
    cerr << e.what() << endl;
    exit(EXIT_FAILURE);
  }

  //Sovle

  const string filename_tmp = outdir + "/frame_origin.vtk";
  shared_ptr<dat_str_core<FLOAT_TYPE, 3>> dat_str = make_shared<dat_str_core<FLOAT_TYPE, 3>>(num_nods);
  constrained_newton<FLOAT_TYPE, 3> imp_euler(coll, dat_str, energy, 20, 1e-4, true, true);

  for (size_t f_id = 0; f_id < num_frame; ++f_id){
    cout << "[frame " << f_id << "]" << endl;
    int res;
    do{
      imp_euler.solve(&nods_4_coll[0][0]);
      res = coll->detect(nods_4_coll);
    }while(!res);
    if(ebf[KIN]!= nullptr)
      dynamic_pointer_cast<momentum<FLOAT_TYPE, 3>>(ebf[KIN])->update_location_and_velocity(nods.data());

    const string filename = outdir  + "/frame_" + to_string(f_id) + ".vtk";
    if(type == "tet")
      mesh_write_to_vtk<FLOAT_TYPE, 4>(filename.c_str(), nods, cells);
    else if (type == "hex")
      mesh_write_to_vtk<FLOAT_TYPE, 8>(filename.c_str(), nods, cells);
  }
  return 0;
}
