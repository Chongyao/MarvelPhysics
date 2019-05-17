#include <string>
#include <iostream>
#include <chrono>

#include <libigl/include/igl/readOBJ.h>
#include <libigl/include/igl/writeOBJ.h>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/SparseCore>
#include <Eigen/SparseCholesky>
#include<Eigen/SparseLU>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/filesystem.hpp>

#include "Point_Sys/src/geometry.h"
#include "Point_Sys/src/gen_points.h"
#include "Point_Sys/src/get_nn.h"
#include "Point_Sys/src/points_energy.h"
#include "Point_Sys/src/data_stream.h"
#include "Point_Sys/src/gen_surf.h"
#include "io.h"
#include "basic_energy.h"
#include "config.h"
#include "implicit_euler.h"

using namespace marvel;
using namespace std;
using namespace Eigen;
using namespace igl;
using namespace chrono;
using namespace boost;

int main(int argc, char** argv){


  
  __TIME_BEGIN__
  Eigen::initParallel();
  cout << "[INFO]>>>>>>>>>>>>>>>>>>>Eigen parallel<<<<<<<<<<<<<<<<<<" << endl;
  cout << "enable parallel in Eigen in " << nbThreads() << " threads" << endl;
  __TIME_END__("hey")
  
  cout << "[INFO]>>>>>>>>>>>>>>>>>>>READ JSON FILE<<<<<<<<<<<<<<<<<<" << endl;
  boost::property_tree::ptree pt;{
    const string jsonfile_path = argv[1];
    
    cout << jsonfile_path << endl;
    const size_t ext = jsonfile_path.rfind(".json");
    if (ext != std::string::npos){
      read_json(jsonfile_path, pt);
      cout << "read json successful" <<endl;
    }
    else{
      cout << "json file extension error" << endl;
      return 0;
    }
  }
  
  cout << "[INFO]>>>>>>>>>>>>>>>>>>>IMPORT MESH<<<<<<<<<<<<<<<<<<" << endl;
  const string mesh_name = pt.get<string>("surf");
  const string indir = pt.get<string>("indir");
  const string outdir = pt.get<string>("outdir") + mesh_name;
  //mkdir
  boost::filesystem::path outpath(outdir);
  if ( !boost::filesystem::exists(outdir) )
    boost::filesystem::create_directories(outdir);

  MatrixXi surf;
  MatrixXd nods;
  readOBJ((indir+mesh_name+".obj").c_str(), nods, surf);
  cout << "surf: " << surf.rows() << " " << surf.cols() << endl << "nods: " << nods.rows() << " " << nods.cols() << endl;
  
  surf.transposeInPlace();
  nods.transposeInPlace();


  cout << "[INFO]>>>>>>>>>>>>>>>>>>>Generate sampled points<<<<<<<<<<<<<<<<<<" << endl;
  MatrixXd points(3,3);
  MatrixXd test(3, 3);
  gen_points(nods, surf, pt.get<size_t>("num_in_axis"), points, true);
  cout << "[INFO]points num is" << points.cols() << endl;
  
  cout << "[INFO]>>>>>>>>>>>>>>>>>>>points<<<<<<<<<<<<<<<<<<" << endl;

  size_t dim = points.cols();
  cout <<"generate points done." << endl;


  cout << "[INFO] Assemble energies..." << endl;
  enum {POTS, CONS, GRAV, MOME};
  vector<std::shared_ptr<Functional<double, 3>>> ebf(MOME + 1);
  std::shared_ptr<dat_str_core<double, 3>>  dat_str = make_shared<energy_dat>(dim);
  

  
  
  cout << "[INFO]>>>>>>>>>>>>>>>>>>>Build spatial hash<<<<<<<<<<<<<<<<<<" << endl;
  spatial_hash SH(points, pt.get<size_t>("nn_num"));
  

  
  cout << "[INFO]>>>>>>>>>>>>>>>>>>>Build Point System<<<<<<<<<<<<<<<<<<" << endl;
  //calc volume 
  double volume = clo_surf_vol(nods, surf);
  //calc support radii
  VectorXd sup_radi = SH.get_sup_radi();
  
  cout << "[INFO]>>>>>>>>>>>>>>>>>>>sup_radi<<<<<<<<<<<<<<<<<<" << endl;

  
  //get friends of every point
  vector<vector<size_t>> friends_all(dim);
#pragma omp parallel for
  for(size_t i = 0; i < dim; ++i){
    SH.get_friends(points.col(i), sup_radi(i), friends_all[i]);
  }

  ebf[POTS] = make_shared<point_sys>(points, pt.get<double>("rho"), pt.get<double>("Young"), pt.get<double>("Poission"), volume, pt.get<double>("kv"), friends_all, sup_radi);
  dynamic_pointer_cast<point_sys>(ebf[POTS])->pre_compute(dat_str);

  
  cout << "[INFO]>>>>>>>>>>>>>>>>>>>Simple Constraint Points<<<<<<<<<<<<<<<<<<" << endl;
  vector<size_t> cons;
  auto cons_file_path = indir + mesh_name +".csv";
  if ( boost::filesystem::exists(cons_file_path) )
    read_fixed_verts_from_csv(cons_file_path.c_str(), cons);
  cout << "constrint " << cons.size() << " points" << endl;

  
  ebf[CONS] = std::make_shared<position_constraint<3>>(dim, pt.get<double>("position_weig"), cons);
    

  cout << "[INFO]>>>>>>>>>>>>>>>>>>>Gravity<<<<<<<<<<<<<<<<<<" << endl;
  double gravity = pt.get<double>("gravity");
  const auto mass_vector = dynamic_pointer_cast<point_sys>(ebf[POTS])->get_Mass_VectorXd();
  ebf[GRAV] = make_shared<gravity_energy<3>>(dim, pt.get<double>("w_g"), gravity,  mass_vector, 'y');
  
  
  cout << "[INFO]>>>>>>>>>>>>>>>>>>>MOMENTUM<<<<<<<<<<<<<<<<<<" << endl;
  double delt_t = pt.get<double>("time_step");
  // momentum MO(dim, PS.get_Mass_Matrix(), delt_t);
  ebf[MOME] = make_shared<momentum<3>>(dim, mass_vector, delt_t);

  
  // cout << "[INFO]>>>>>>>>>>>>>>>>>>>Collision<<<<<<<<<<<<<<<<<<" << endl;
  // collision COLL(pt.get<double>("w_coll"), 'y', pt.get<double>("coll_pos"), static_cast<size_t>(nods.cols()), dim);

  //energy all
  std::shared_ptr<Functional<double, 3>> energy;
  try {
    energy = make_shared<energy_t<double, 3>>(ebf);

  } catch ( std::exception &e ) {
    cerr << e.what() << endl;
    exit(EXIT_FAILURE);
  }

  
  cout << "[INFO]>>>>>>>>>>>>>>>>>>>SOlVE<<<<<<<<<<<<<<<<<<" << endl;
  //initilize variables in time integration


  MatrixXd displace;
  MatrixXd vet_displace;
  MatrixXd points_now;
  displace.setZero(3, dim);
  vet_displace.setZero(3, nods.cols());

  
  size_t iters_perframe = floor(1.0/delt_t/pt.get<size_t>("rate"));
  newton_iter<double, 3> imp_euler(dat_str, energy, delt_t);
  
  for(size_t i = 0; i < pt.get<size_t>("max_iter"); ++i){
    
    cout << "[INFO]>>>>>>>>>>>>>>>>>>>iter "<< i <<" <<<<<<<<<<<<<<<<<<" << endl;
    // cout << "displace is " << endl<< displace.block(0, 0, 3, 8) << endl;

    //newtown iter
    imp_euler.solve(displace.data());
    
    
    dynamic_pointer_cast<momentum<3>>(ebf[MOME])->update_location_and_velocity(displace.data());
    auto surf_filename = outdir  + "/" + mesh_name + "_" + to_string(i) + ".obj";
    auto point_filename = outdir + "/" + mesh_name + "_points_" + to_string(i) + ".vtk";

    points_now = points + displace;
    vet_displace = displace.block(0, 0, 3, nods.cols());
    point_write_to_vtk(point_filename.c_str(), points_now.data(), dim);
    writeOBJ(surf_filename.c_str(), (nods + vet_displace).transpose(), surf.transpose());
    


  }
  //done

  return 0;
}

  




