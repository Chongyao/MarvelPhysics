#include "basic_energy.h"
#include "implicit_euler.h"
#include "io.h"

#include "FEM/src/elas_energy.h"
#include "FEM/src/mass_matrix.h"

#include <iostream>

using namespace std;
using namespace Eigen;
using namespace marvel;

using TET_ELAS = BaseElas<double, 3, 4, 1, 1, stvk, basis_func, quadrature>;
int main(int argc, char** argv){
  std::cout.precision(10);
  const char* filename = argv[1];
  
  MatrixXd nods(1, 1);
  MatrixXi tets(1, 1);
  tet_mesh_read_from_vtk(filename, nods, tets);
  const size_t num_nods = nods.cols();
  cout <<"V"<< nods.rows() << " " << nods.cols() << endl << "T " << tets.rows() << " "<< tets.cols() << endl;

  const string outdir = argv[3];
  
  //set mtr
  constexpr  double rho = 5;
  constexpr  double Young = 5000.0;
  constexpr  double poi = 0.45;
  constexpr  double gravity = 0;
  constexpr  double dt = 0.01;
  const      double w_pos = 1e4;
  const      size_t num_frame = 100;

  //read fixed points
  vector<size_t> cons(0);
  const char* cons_file_path = argv[2];
  // if ( boost::filesystem::exists(cons_file_path) )
  read_fixed_verts_from_csv(cons_file_path, cons);
  cout << "constrint " << cons.size() << " points" << endl;

  //read rotated points
  size_t center = 94;
  vector<size_t> rotate(0);
  const char* rotate_file_path = argv[4];
  cout << rotate_file_path << endl;
  read_fixed_verts_from_csv(rotate_file_path, rotate);
  cout << "rorate " << rotate.size() << " points" << endl;  
  
  
  //calc mass vector
  VectorXd mass_vec(nods.rows() * num_nods);

  calc_mass_vector<double>(nods, tets, rho, mass_vec);
  cout <<"total mass is "<< mass_vec.sum() << endl;
  cout << "build energy" << endl;
  enum energy_type{ELAS, GRAV, KIN, POS};
  vector<shared_ptr<Functional<double, 3>>> ebf(POS + 1);{
    ebf[ELAS] = make_shared<TET_ELAS>(nods, tets, Young, poi);
    ebf[GRAV] = make_shared<gravity_energy<3>>(num_nods, 1, gravity, mass_vec, 'z');
    ebf[KIN] = make_shared<momentum<3>>(nods.data(), num_nods, mass_vec, dt);
    ebf[POS] = make_shared<position_constraint<3>>(nods.data(), num_nods, w_pos, cons);
    }
  cout << "assemble energy" << endl;


  //give velocity
  Vector3d w = Vector3d::Zero(); w(0) = 10;
  for(size_t i = 0; i < rotate.size(); ++i){
    size_t rot_id = rotate[i];
    Vector3d r = nods.col(rot_id) - nods.col(center);
    dynamic_pointer_cast<momentum<3>>(ebf[KIN])->vk_.segment<3>(rot_id * 3) = r.cross(w);
  }
  
  
  
  shared_ptr<Functional<double, 3>> energy;
  try {
    energy = make_shared<energy_t<double, 3>>(ebf);

  } catch ( std::exception &e ) {
    cerr << e.what() << endl;
    exit(EXIT_FAILURE);
  }

  //Sovle
  const string filename_tmp = outdir  + "/frame_origin.vtk";
  tet_mesh_write_to_vtk(filename_tmp.c_str(), nods, tets);
  shared_ptr<dat_str_core<double, 3>>  dat_str = make_shared<dat_str_core<double, 3>>(num_nods);
  newton_iter<double, 3> imp_euler(dat_str, energy, dt, 20, 1e-4, true, false);
  
  Map<VectorXd> xk(nods.data(), nods.size() );
  VectorXd new_nods = xk;
  
  for(size_t f_id = 0; f_id < num_frame; ++f_id){
    cout << "[frame " << f_id << "]" << endl;
    imp_euler.solve(new_nods.data());
    xk = new_nods;
    dynamic_pointer_cast<momentum<3>>(ebf[KIN])->update_location_and_velocity(new_nods.data());
    //give velocity
    Vector3d w = Vector3d::Zero(); w(0) = 10;
    for(size_t i = 0; i < rotate.size(); ++i){
    size_t rot_id = rotate[i];
    Vector3d r = nods.col(rot_id) - nods.col(center);
    dynamic_pointer_cast<momentum<3>>(ebf[KIN])->vk_.segment<3>(rot_id * 3) = r.cross(w);

    }

    const string filename = outdir  + "/frame_" + to_string(f_id) + ".vtk";
    tet_mesh_write_to_vtk(filename.c_str(), nods, tets);
  }

  


  



  
  
  return 0;
}
