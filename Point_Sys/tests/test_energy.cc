#include <string>
#include <iostream>
#include <chrono>

#include <libigl/include/igl/readOBJ.h>
#include <libigl/include/igl/writeOBJ.h>
#include <Eigen/Core>
#include <Eigen/Geometry>
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



using namespace marvel;
using namespace std;
using namespace Eigen;
using namespace igl;
using namespace chrono;
using namespace boost;

int main(int argc, char** argv){

  
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
  MatrixXd points;
  gen_points(nods, surf, pt.get<size_t>("num_in_axis"), points);
  size_t dim = points.cols();

  
  cout << "[INFO]>>>>>>>>>>>>>>>>>>>Build spatial hash<<<<<<<<<<<<<<<<<<" << endl;
  spatial_hash SH(points);

  
  cout << "[INFO]>>>>>>>>>>>>>>>>>>>Build Point System<<<<<<<<<<<<<<<<<<" << endl;
  //calc volume 
  double volume = clo_surf_vol(nods, surf);
  //calc support radii
  VectorXd sup_radi = SH.get_sup_radi();
  //get friends of every point
  vector<vector<size_t>> friends_all(dim);
#pragma parallel omp for
  for(size_t i = 0; i < dim; ++i){
    SH.get_friends(points.col(i), sup_radi(i), friends_all[i]);
  }
  
  point_sys PS(points, pt.get<double>("rho"), pt.get<double>("Young"), pt.get<double>("Poission"), volume, pt.get<double>("kv"), friends_all, sup_radi);

  
  cout << "[INFO]>>>>>>>>>>>>>>>>>>>Construct Deform_Surf<<<<<<<<<<<<<<<<<<" << endl;
  //get friends of every vertex in mesh
  double kernel_cof = pt.get<double>("ker_cof");
  vector<vector<size_t>> vet_fris(nods.cols());
  for(size_t j = 0; j < nods.cols(); ++j){
    vector<size_t> fris;
    SH.get_friends(nods.col(j), kernel_cof, fris);
    vet_fris[j] = fris;
  }
  
  deform_surf_MLS<double> DS(surf, nods, points, vet_fris, kernel_cof);


  
  cout << "[INFO]>>>>>>>>>>>>>>>>>>>Simple Constraint Points<<<<<<<<<<<<<<<<<<" << endl;
  //add simple constraints
  //This should read from file. We loop for some points to restrain here.
  //Constraints vary from different models and situations.
  vector<size_t> cons;
  for(size_t i = 0; i < points.cols(); ++i){
    if(points(2, i) > 2){
      cons.push_back(i);
      cout << i << " ";
    }
  }
  cout << endl;

  
  
  cout << "[INFO]>>>>>>>>>>>>>>>>>>>SOlVE<<<<<<<<<<<<<<<<<<" << endl;
  //initilize variables in time integration
  energy_dat dat_str (dim);

  double delt_t = pt.get<double>("time_step");
  MatrixXd displace;
  MatrixXd velocity;
  MatrixXd acce;
  MatrixXd new_acce;
  MatrixXd gra;
  MatrixXd vet_displace;
  displace.setZero(3, dim);
  velocity.setZero(3, dim);
  acce.setZero(3, dim);
  gra.setZero(3, dim);
  new_acce.setZero(3, dim);
  vet_displace.setZero(3, nods.cols());



  for(size_t i = 0; i < pt.get<size_t>("max_iter"); ++i){
    cout << "iter is "<<endl<< i << endl;
    cout << "displace is " << endl<< displace.block(0, 0, 3, 8) << endl;
    cout << "velocity is "<<endl<< velocity.block(0, 0, 3, 8) << endl;

    PS.calc_defo_gra(displace.data(), dat_str);
    PS.Gra(displace.data(), dat_str);
    PS.gravity(displace.data(), dat_str, 9.8);
    
#pragma parallel omp for
    for(size_t j = 0; j < dim; ++j){
      assert(PS.get_mass(j) > 0);
      new_acce.col(j) = dat_str.gra_.col(j)/PS.get_mass(j);
    }
    cout << "new acce is "<<endl << new_acce.block(0, 0, 3, 8) << endl;
    velocity += 0.5*(new_acce + acce)*delt_t;
#pragma parallel omp for
    for(auto c : cons){
      velocity.col(c) = MatrixXd::Zero(3, 1);
    }
#pragma parallel omp for
    displace += velocity*delt_t + 0.5*acce*delt_t*delt_t;
    for(auto c : cons){
      displace.col(c) = MatrixXd::Zero(3, 1);
    }
    vet_displace = DS.update_surf(displace, dat_str.def_gra_);
    acce = new_acce;

    dat_str.set_zero();
    auto surf_filename = outdir  + "/" + mesh_name + "_" + to_string(i) + ".vtk";
    auto point_filename = outdir + "/" + mesh_name + "_points_" + to_string(i) + ".vtk";
    MatrixXd points_now = points + displace;
    point_write_to_vtk(point_filename.c_str(), points_now.data(), dim);
    tri_mesh_write_to_vtk(surf_filename.c_str(), nods + vet_displace, surf);
  }

  //done
  
}


