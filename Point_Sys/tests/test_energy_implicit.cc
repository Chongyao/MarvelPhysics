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
  gen_points(nods, surf, pt.get<size_t>("num_in_axis"), points);

  // #if 1
 // points = nods;
  // #endif
  size_t dim = points.cols();
  cout <<"generate points done." << endl;
  
  cout << "[INFO]>>>>>>>>>>>>>>>>>>>Build spatial hash<<<<<<<<<<<<<<<<<<" << endl;
  spatial_hash SH(points, pt.get<size_t>("nn_num"));

  
  cout << "[INFO]>>>>>>>>>>>>>>>>>>>Build Point System<<<<<<<<<<<<<<<<<<" << endl;
  //calc volume 
  double volume = clo_surf_vol(nods, surf);
  //calc support radii
  VectorXd sup_radi = SH.get_sup_radi();
  //get friends of every point
  vector<vector<size_t>> friends_all(dim);
#pragma omp parallel for
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
    SH.get_friends(nods.col(j), kernel_cof, fris, false);
    vet_fris[j] = fris;
  }
  
  deform_surf_MLS<double> DS(surf, nods, points, vet_fris, kernel_cof);



  cout << "[INFO]>>>>>>>>>>>>>>>>>>>Simple Constraint Points<<<<<<<<<<<<<<<<<<" << endl;

  //add simple constraints
  //This should read from file. We loop for some points to restrain here.
  //Constraints vary from different models and situations.
  vector<size_t> cons;
  for(size_t i = 0; i < points.cols(); ++i){
    if(points(2, i) > 0.7 ){
      cons.push_back(i);
      cout << i << " ";
    }
  }
  cout << endl;
  position_constraint pos_cons(pt.get<double>("position_weig"), cons, dim);


  
  
  cout << "[INFO]>>>>>>>>>>>>>>>>>>>SOlVE<<<<<<<<<<<<<<<<<<" << endl;
  //initilize variables in time integration
  energy_dat dat_str (dim);

  double delt_t = pt.get<double>("time_step");
  SparseVector displace(3 * dim);
  SparseVector velocity(3 * dim);
  MatrixXd acce;
  MatrixXd new_acce;
  MatrixXd gra;
  MatrixXd vet_displace;
  velocity.setZero(3, dim);
  acce.setZero(3, dim);
  gra.setZero(3, dim);
  new_acce.setZero(3, dim);
  vet_displace.setZero(3, nods.cols());

  PS.pre_compute(dat_str);
  size_t iters_perframe = floor(1/delt_t/20);

  SparseMatrix A_CG(dim * 3, dim * 3);
  SparseVector b_CG(dim * 3);
  SparseMatrix M = PS.get_Mass_Matrix();
  
  for(size_t i = 0; i < pt.get<size_t>("max_iter"); ++i){
    cout << "iter is "<<endl<< i << endl;
    cout << "displace is " << endl<< displace.block(0, 0, 3, 8) << endl;
    // cout << "velocity is "<<endl<< velocity.block(0, 0, 3, 8) << endl;

    PS.calc_defo_gra(displace.data(), dat_str);
    PS.Gra(displace.data(), dat_str);
    PS.gravity(displace.data(), dat_str, pt.get<double>("gravity"));
    pos_cons.Gra(displace.data(), dat_str);
    PS.Hessian(displace.data(), dat_str);
    pos_cons.Hes(displace.data(), dat_str);
    
    //implicit time integral
    { SparseVector<double> velocity(3 * dim)
      A.setZero();
      dat_str.hes_.setFromTriplets(dat_str.hes_trips.begin(), dat_str.hes_trips.end());
      A = M - delt_t*delt_t*dat_str.hes_;
      b = M * 
     
    }
    
    cout << "[INFO]>>>>>>>>>>>>>>>>>>>Elasticity Energy Val<<<<<<<<<<<<<<<<<<" << endl;
    cout << dat_str.ela_val_.transpose();
    cout << "[INFO]>>>>>>>>>>>>>>>>>>>VOL conservation val<<<<<<<<<<<<<<<<<<" << endl;
    cout << dat_str.vol_val_.transpose();
        
    vet_displace = DS.update_surf(displace, dat_str.def_gra_);
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
      
      tri_mesh_write_to_vtk(surf_filename.c_str(), nods + vet_displace, surf);
    }

    dat_str.set_zero();
  }
  //done
}




