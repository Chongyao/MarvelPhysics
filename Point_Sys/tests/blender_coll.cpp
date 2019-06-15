
#define EIGEN_USE_BLAS
// #define NDEBUG
// #define EIGEN_USE_LAPACKE
// #define EIGEN_USE_LAPACKE_STRICT


#include <string>
#include <iostream>
#include <chrono>

#include <memory>
#include <set>
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
#include "basic_energy.h"

#include "vtk2surf.h"

#include "coll_wrapper/coll_wrapper.h"
#include "implicit_euler.h"
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
  const string indir = "../input/" + mesh_name;
  const string outdir = "../output/" + mesh_name;
  //mkdir
  boost::filesystem::path outpath(outdir);
  if ( !boost::filesystem::exists(outdir) )
    boost::filesystem::create_directories(outdir);

  MatrixXi surf;
  MatrixXd nods;
  if(!readOBJ((indir + '/' +mesh_name+".obj").c_str(), nods, surf)){
    cout << "model OBJ file wrong" << endl;
    return 1;    
  }

  cout << "surf: " << surf.rows() << " " << surf.cols() << endl << "nods: " << nods.rows() << " " << nods.cols() << endl;
  
  surf.transposeInPlace();
  nods.transposeInPlace();
  

  
  cout << "[INFO]>>>>>>>>>>>>>>>>>>>Generate sampled points<<<<<<<<<<<<<<<<<<" << endl;
  MatrixXd points(3,3);
  MatrixXd test(3, 3);
  gen_points(nods, surf, simulation_para.get<size_t>("num_in_axis"), points, true);
  cout << points.rows() << " " << points.cols() << endl;

  const size_t dim = points.cols();
  cout <<"generate points done." << endl;
  std::shared_ptr<dat_str_core<double, 3>> dat_str = make_shared<energy_dat>(dim);  
  
  cout << "[INFO] Assemble energies..." << endl;
  enum {POTS, CONS,  GRAV,  MOME};
  vector<std::shared_ptr<Functional<double, 3>>> ebf(MOME + 1); 
  

  
  cout << "[INFO]>>>>>>>>>>>>>>>>>>>POINTS ENERGY<<<<<<<<<<<<<<<<<<" << endl;
  cout << "[INFO]Build spatial " << endl;
  spatial_hash SH(points, simulation_para.get<size_t>("nn_num"));
  cout << "[INFO]Build Point System" << endl;
  //calc volume 
  const double volume = clo_surf_vol(nods, surf);
  //calc support radii
  VectorXd sup_radi = SH.get_sup_radi();
  cout << "[INFO]sup_radi" << endl;
  //get friends of every point
  vector<vector<size_t>> friends_all(dim);
  #pragma omp parallel for
  for(size_t i = 0; i < dim; ++i){
    SH.get_friends(points.col(i), sup_radi(i), friends_all[i]);
  }

  ebf[POTS] = make_shared<point_sys>(points, common.get<double>("density"), physics_para.get<double>("Young"), physics_para.get<double>("Poission"), volume, simulation_para.get<double>("kv"), friends_all, sup_radi);

  dynamic_pointer_cast<point_sys>(ebf[POTS])->pre_compute(dat_str);



  cout << "[INFO]>>>>>>>>>>>>>>>>>>>Simple Constraint Points<<<<<<<<<<<<<<<<<<" << endl;

  vector<size_t> cons(0);
  auto cons_file_path = indir + "/" +  mesh_name +".csv";
  if ( boost::filesystem::exists(cons_file_path) )
    read_fixed_verts_from_csv(cons_file_path.c_str(), cons);
  cout << "constrint " << cons.size() << " points" << endl;
  ebf[CONS] = std::make_shared<position_constraint<3>>(dim, simulation_para.get<double>("position_weig"), cons);
  
  cout << "[INFO]>>>>>>>>>>>>>>>>>>>Gravity<<<<<<<<<<<<<<<<<<" << endl;
  const double gravity = common.get<double>("gravity");
  const auto mass_vector = dynamic_pointer_cast<point_sys>(ebf[POTS])->get_Mass_VectorXd();
  ebf[GRAV] = make_shared<gravity_energy<3>>(dim, common.get<double>("gravity"), gravity,  mass_vector, 'z');

  cout << "[INFO]>>>>>>>>>>>>>>>>>>>MOMENTUM<<<<<<<<<<<<<<<<<<" << endl;
  double delt_t = common.get<double>("time_step");
  // momentum MO(dim, PS.get_Mass_Matrix(), delt_t);
  ebf[MOME] = make_shared<momentum<3>>(dim, mass_vector, delt_t);

  
  cout << "[INFO]>>>>>>>>>>>>>>>>>>>COLLISION<<<<<<<<<<<<<<<<<<" << endl;
  vector<std::shared_ptr<MatrixXi> > obta_surfs;
  vector<std::shared_ptr<MatrixXd>> obta_nods;{
    filesystem::path obstacles_path(indir + "/" + "obstacles");
    filesystem::directory_iterator end_iter;
    for(filesystem::directory_iterator file_iter(obstacles_path); file_iter != end_iter; ++file_iter){
      cout << file_iter->path() << endl;
      if(filesystem::is_regular_file(file_iter->path())){
        string one_obstacle = file_iter->path().string();
        MatrixXd one_obta_nods;
        MatrixXi one_obta_surf;
        if(!readOBJ(one_obstacle.c_str(), one_obta_nods, one_obta_surf))
          continue;
        one_obta_surf.transposeInPlace();
        one_obta_nods.transposeInPlace();
        obta_surfs.push_back(std::move(make_shared<MatrixXi>(one_obta_surf)));
        obta_nods.push_back(std::move(make_shared<MatrixXd>(one_obta_nods)));
        cout << one_obstacle << endl;
      }
    }    
  }


  size_t num_fake_tris = dim%3 ? dim / 3 + 1: dim / 3 ;
  cout << " face surf tris num is " << num_fake_tris << endl;
  std::shared_ptr<MatrixXi> fake_surf_ptr = make_shared<MatrixXi>(3, num_fake_tris);{
    #pragma omp parallel for
    for(size_t i = 0; i < num_fake_tris; ++i){
      (*fake_surf_ptr)(0, i) = i * 3;
      (*fake_surf_ptr)(1, i) = i * 3 + 1 >= dim ? i * 3 - 1 : i * 3 + 1;
      (*fake_surf_ptr)(2, i) = i * 3 + 2 >= dim ? i * 3 - 2 : i * 3 + 2;
    }
  }
  
  // auto coll_ptr = make_shared<coll_wrapper>(obta_surfs, obta_nods, fake_surf_ptr, points);
  


  cout << "[INFO]>>>>>>>>>>>>>>>>>>>SOlVE<<<<<<<<<<<<<<<<<<" << endl;
  //initilize variables in time integration


  string solver = simulation_para.get<string>("solver");


  MatrixXd points_pos;
  MatrixXd displace;
  MatrixXd velocity;
  MatrixXd acce;

  MatrixXd vet_displace;
  points_pos.setZero(3, dim); 
  displace.setZero(3, dim); 
  velocity.setZero(3, dim); 
  acce.setZero(3, dim); 


  vet_displace.setZero(3, nods.cols());
  Eigen::Map<VectorXd> acce_vec(acce.data(), 3 * dim);
  Eigen::Map<VectorXd> velo_vec(velocity.data(), 3 * dim);
  Eigen::Map<MatrixXd> gra_vec(dat_str->get_gra().data(), 3, dim);

  size_t iters_perframe = static_cast<size_t>(round(1.0/delt_t/common.get<size_t>("frame_rate")));
  size_t max_iter  = static_cast<size_t>(ceil(common.get<double>("total_time") / delt_t));
  cout << "max iter is " << max_iter << endl;
  double dump = simulation_para.get<double>("dump");


  auto start = system_clock::now();
  size_t frame_id = 0;
  if(solver == "explicit"){
    ebf[MOME] = nullptr;
    // ebf[GRAV] = nullptr;
    ebf[CONS] = nullptr;
    std::shared_ptr<Functional<double, 3>> energy;
    try {
      energy = make_shared<energy_t<double, 3>>(ebf);

    } catch ( std::exception &e ) {
      cerr << e.what() << endl;
      exit(EXIT_FAILURE);
    }

    for(size_t i = 0; i < max_iter; ++i){
      dat_str->set_zero();
      energy->Val(displace.data(), dat_str);
      energy->Gra(displace.data(), dat_str);
         
     #pragma omp parallel for
      for(size_t j = 0; j < dim; ++j){
        
        acce.col(j) = -gra_vec.col(j)/mass_vector(j) - velocity.col(j)*dump;
      }
      for(const auto& c: cons){
        acce.col(c).setZero();
      }

      // acce_vec = (dat_str->get_gra()).array() / mass_vector.array() -  velo_vec.array() * dump;
      // cout << acce.row(0) << endl;



      velocity += delt_t * acce;
      displace += delt_t *velocity;
      // points_pos = points + displace;

      //>>>>>>>>>>>>>>>>>>COLLID<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<//
      
      // coll_ptr->Collide(obta_surfs, obta_nods, velocity, points_pos);
      // displace = points_pos - points;      
      //>>>>>>>>>>>>>>>>>>COLLID<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<//

    

      
      if(i%iters_perframe == 0){
        cout << "frame is " << frame_id << endl;
        auto surf_filename = outdir  + "/" + mesh_name + "_" + to_string(frame_id) + ".obj";
        auto point_filename = outdir + "/" + mesh_name + "_points_" + to_string(frame_id) + ".vtk";
        MatrixXd points_now = points + displace;
        point_write_to_vtk(point_filename.c_str(), points_now.data(), dim);
        
        vet_displace = displace.block(0, 0, 3, nods.cols());
        writeOBJ(surf_filename.c_str(), (nods + vet_displace).transpose(), surf.transpose());
        ++frame_id;

      }

    
    }

  }
  else{//TODO:need to be rewrite
    std::shared_ptr<Functional<double, 3>> energy;
    try {
      energy = make_shared<energy_t<double, 3>>(ebf);

    } catch ( std::exception &e ) {
      cerr << e.what() << endl;
      exit(EXIT_FAILURE);
    }

    newton_iter<double, 3> imp_euler(dat_str, energy, delt_t, 20, 1e-2, true, false);
    for(size_t i = 0; i < max_iter; ++i){
      cout << "iter is "<< i << endl;
      imp_euler.solve(displace.data());
      dynamic_pointer_cast<momentum<3>>(ebf[MOME])->update_location_and_velocity(displace.data());

      if(i%iters_perframe == 0){
        auto surf_filename = outdir  + "/" + mesh_name + "_" + to_string(frame_id) + ".obj";
        auto point_filename = outdir + "/" + mesh_name + "_points_" + to_string(frame_id) + ".vtk";
        MatrixXd points_now = points + displace;
        point_write_to_vtk(point_filename.c_str(), points_now.data(), dim);
        
        vet_displace = displace.block(0, 0, 3, nods.cols());
        writeOBJ(surf_filename.c_str(), (nods + vet_displace).transpose(), surf.transpose());
        ++frame_id;

      }
    }

    
    
  }
  auto end = system_clock::now();
  auto duration = duration_cast<microseconds>(end - start);
  cout <<  "花费了" 
       << double(duration.count()) * microseconds::period::num / microseconds::period::den 
       << "秒" << endl;

  return 0;
  //done
}





