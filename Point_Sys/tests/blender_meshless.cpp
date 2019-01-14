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
#include "Point_Sys/src/basic_energy.cc"

#include "vtk2surf.h"



using namespace marvel;
using namespace std;
using namespace Eigen;
using namespace igl;
using namespace chrono;
using namespace boost;

int main(int argc, char** argv){
  
  Eigen::initParallel();
  cout << "[INFO]>>>>>>>>>>>>>>>>>>>Eigen parallel<<<<<<<<<<<<<<<<<<" << endl;
  cout << "enable parallel in Eigen in " << nbThreads() << " threads" << endl;
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
  
  auto common = pt.get_child("common");
  auto blender = pt.get_child("blender");
  auto physics_para = pt.get_child("physics_para");
  auto simulation_para = pt.get_child("simulation_para");
  
  cout << "[INFO]>>>>>>>>>>>>>>>>>>>IMPORT MESH<<<<<<<<<<<<<<<<<<" << endl;
  const string mesh_name = blender.get<string>("surf");
  const string indir = "../input";
  const string outdir = "../output/" + mesh_name;
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
  gen_points(nods, surf, simulation_para.get<size_t>("num_in_axis"), points, true);
  cout << points.rows() << " " << points.cols() << endl;
  // #if 1
 // points = nods;
  // #endif
  size_t dim = points.cols();
  cout <<"generate points done." << endl;
  
  cout << "[INFO]>>>>>>>>>>>>>>>>>>>Build spatial hash<<<<<<<<<<<<<<<<<<" << endl;
  spatial_hash SH(points, simulation_para.get<size_t>("nn_num"));

  
  cout << "[INFO]>>>>>>>>>>>>>>>>>>>Build Point System<<<<<<<<<<<<<<<<<<" << endl;
  //calc volume 
  double volume = clo_surf_vol(nods, surf);
  //calc support radii
  VectorXd sup_radi = SH.get_sup_radi();
  
  cout << "[INFO]>>>>>>>>>>>>>>>>>>>sup_radi<<<<<<<<<<<<<<<<<<" << endl;
  cout << sup_radi << endl;
  //get friends of every point
  vector<vector<size_t>> friends_all(dim);
// #pragma omp parallel for
  for(size_t i = 0; i < dim; ++i){
    SH.get_friends(points.col(i), sup_radi(i), friends_all[i]);
  }



  
  point_sys PS(points, common.get<double>("rho"), physics_para.get<double>("Young"), physics_para.get<double>("Poission"), volume, simulation_para.get<double>("kv"), friends_all, sup_radi);

  
  cout << "[INFO]>>>>>>>>>>>>>>>>>>>Construct Deform_Surf<<<<<<<<<<<<<<<<<<" << endl;
  // //get friends of every vertex in mesh
  // double kernel_cof = pt.get<double>("ker_cof");
  // vector<vector<size_t>> vet_fris(nods.cols());
  // for(size_t j = 0; j < nods.cols(); ++j){
  //   vector<size_t> fris;
  //   SH.get_friends(nods.col(j), kernel_cof, fris, false);
  //   vet_fris[j] = fris;
  // }
  
  // deform_surf_MLS<double> DS(surf, nods, points, vet_fris, kernel_cof);



  cout << "[INFO]>>>>>>>>>>>>>>>>>>>Simple Constraint Points<<<<<<<<<<<<<<<<<<" << endl;
  //add simple constraints
  //This should read from file. We loop for some points to restrain here.
  //Constraints vary from different models and situations.
  vector<size_t> cons(0);

  cout << endl;
  position_constraint pos_cons(simulation_para.get<double>("position_weig"), cons, dim);
  
  cout << "[INFO]>>>>>>>>>>>>>>>>>>>Gravity<<<<<<<<<<<<<<<<<<" << endl;
  double gravity = common.get<double>("gravity");
  gravity_energy GE(simulation_para.get<double>("w_g"), gravity, dim, PS.get_Mass_VectorXd(), 'y');

  
  cout << "[INFO]>>>>>>>>>>>>>>>>>>>COLLISION<<<<<<<<<<<<<<<<<<" << endl;
  collision COLL(simulation_para.get<double>("w_coll"),'y', simulation_para.get<double>("g_pos"), nods.cols(), dim);

  
  cout << "[INFO]>>>>>>>>>>>>>>>>>>>SOlVE<<<<<<<<<<<<<<<<<<" << endl;
  //initilize variables in time integration
  energy_dat dat_str (dim);


  double delt_t = common.get<double>("time_step");
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

  PS.pre_compute(dat_str);
  size_t iters_perframe = floor(1/delt_t/common.get<size_t>("frame_rate"));
  size_t max_iter  = static_cast<size_t>(ceil(common.get<size_t>("totao_time") / delt_t));
  double dump = simulation_para.get<double>("dump");
  double previous_step_Val = 0;

  auto start = system_clock::now();
  for(size_t i = 0; i < max_iter; ++i){
    cerr << "iter is "<<endl<< i << endl;
    cout << "displace is " << endl<< displace.block(0, 0, 3, 8) << endl;
    // cout << "velocity is "<<endl<< velocity.block(0, 0, 3, 8) << endl;

    PS.Val(displace.data(), dat_str);
    PS.Gra(displace.data(), dat_str);

    #if 0
    cout << "[INFO]>>>>>>>>>>>>>>>>>>>difference check<<<<<<<<<<<<<<<<<<" << endl;
    
    { 
      cout << "[INFO]>>>>>>>>>>>>>>>>>>>gra by formula<<<<<<<<<<<<<<<<<<" << endl;
      cout << dat_str.gra_ << endl;
      auto init_val = dat_str.Val_;
      cout << init_val << endl << endl << endl << endl;
      MatrixXd new_gra(3, dim);
      double delt_x = 1e-14;
    
      for(size_t i = 0; i < points.size(); ++i){
        dat_str.set_zero();
        displace(i) += delt_x;
        PS.Val(displace.data(), dat_str);
        new_gra(i) = dat_str.Val_ - init_val;
        displace(i) -= delt_x;
      }
      new_gra /= delt_x;
      cout << -new_gra << endl;
      PS.Val(displace.data(), dat_str);
      PS.Gra(displace.data(), dat_str);
      dat_str.gra_ = -new_gra;
    }
    #endif
    GE.Val(displace.data(), dat_str);
    GE.Gra(displace.data(), dat_str);
    cout << dat_str.gra_.array().square().sum() << endl;
    cout << "[INFO]>>>>>>>>>>>>>>>>>>>Gra norm<<<<<<<<<<<<<<<<<<" << endl;
    


    COLL.Val(nods.data(), displace.data(), dat_str);
    COLL.Gra(nods.data(), displace.data(), dat_str);
    // cout << new_acce - new_acce.col(0) * MatrixXd::Ones(1, dim);    
    // pos_cons.Gra(displace.data(), dat_str);
#pragma omp parallel for
    for(size_t j = 0; j < dim; ++j){
      assert(PS.get_mass(j) > 0);
      new_acce.col(j) = dat_str.gra_.col(j)/PS.get_mass(j) - velocity.col(j)*dump;
    }
    // cout << "new acce is "<<endl << new_acce.block(0, 0, 3, 8) << endl;


    
    // velocity += 0.5*(new_acce + acce)*delt_t;

    velocity += delt_t * new_acce;
// #pragma omp parallel for
//     for(size_t j = 0; j < cons.size(); ++j){
//       velocity.col(cons[j]) = MatrixXd::Zero(3, 1);    
//     }

    // displace += velocity*delt_t + 0.5*acce*delt_t*delt_t;
    displace += delt_t *velocity;
// #pragma omp parallel for
//     for(size_t j = 0; j < cons.size(); ++j){
//       displace.col(cons[j]) = MatrixXd::Zero(3, 1);    
//     }




    cout << "delt energy :"<< dat_str.Val_ - previous_step_Val << endl;
    if (i > 10 && fabs(dat_str.Val_ - previous_step_Val) < 1e-6)
      break;
    
    previous_step_Val = dat_str.Val_;
    cout << "[INFO]>>>>>>>>>>>>>>>>>>>Energy Val<<<<<<<<<<<<<<<<<<" << endl;
    cout << "totao energy: " << dat_str.Val_ << endl;
    // cout << "[INFO]>>>>>>>>>>>>>>>>>>>VOL conservation val<<<<<<<<<<<<<<<<<<" << endl;
    // cout << dat_str.vol_val_.transpose();
        

    acce = new_acce;

    if(i%iters_perframe == 0){ 
      auto surf_filename = outdir  + "/" + mesh_name + "_" + to_string(i) + ".vtk";
      auto point_filename = outdir + "/" + mesh_name + "_points_" + to_string(i) + ".vtk";
      MatrixXd points_now = points + displace;
      point_write_to_vtk(point_filename.c_str(), points_now.data(), dim);
      point_vector_append2vtk(false, point_filename.c_str(), velocity, dim, "velocity");
      point_vector_append2vtk(true, point_filename.c_str(), acce, dim, "accelarate");
      point_scalar_append2vtk(true, point_filename.c_str(), dat_str.ela_val_, dim, "strain_Energy");
      point_scalar_append2vtk(true, point_filename.c_str(), dat_str.vol_val_, dim, "vol_conservation_Energy");
      // vet_displace = DS.update_surf(displace, dat_str.def_gra_);
      vet_displace = displace.block(0, 0, 3, nods.cols());
      tri_mesh_write_to_vtk(surf_filename.c_str(), nods + vet_displace, surf);

    }

    dat_str.set_zero();
  }
  auto end = system_clock::now();
  auto duration = duration_cast<microseconds>(end - start);
  cout <<  "花费了" 
       << double(duration.count()) * microseconds::period::num / microseconds::period::den 
       << "秒" << endl;

  //done
}




