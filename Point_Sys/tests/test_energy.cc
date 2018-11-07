#include <string>
#include <iostream>
#include <chrono>

#include <libigl/include/igl/readOBJ.h>
#include <libigl/include/igl/writeOBJ.h>
#include <Eigen/Core>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

#include "Point_Sys/src/geometry.h"
#include "Point_Sys/src/gen_points.h"
#include "Point_Sys/src/get_nn.h"
#include "Point_Sys/src/points_energy.h"
#include "Point_Sys/src/data_stream.h"
#include "Point_Sys/src/gen_surf.h"
#include "io.h"
#include <Eigen/Geometry>


using namespace marvel;
using namespace std;
using namespace Eigen;
using namespace igl;
using namespace chrono;
using namespace boost;

// const MatrixXd& ele_mat(const size_t &ele_id, const size_t &rows, const size_t &cols, MatrixXd &&mat){
//   return Map<MatrixXd>(mat.col(ele_id).data(), rows, cols);
// }

int main(int argc, char** argv){

  
  boost::property_tree::ptree pt;{
    const string jsonfile_path = argv[1];
    
    cout << "[INFO]>>>>>>>>>>>>>>>>>>>json file path<<<<<<<<<<<<<<<<<<" << endl;
    cout << jsonfile_path << endl;
    const size_t ext = jsonfile_path.rfind(".json");
    if (ext != std::string::npos){
      read_json(jsonfile_path, pt);
      cout << "read json successful" <<pt.get<string>("surf")<<endl;
    }
    else{
      cout << "json file extension error" << endl;
      return 0;
    }
  }
  

  MatrixXi surf;
  MatrixXd nods;
  // jtf::mesh::load_obj(pt.get<string>("surf.value").c_str(), surf, nods);
  readOBJ(pt.get<string>("surf").c_str(), nods, surf);
  cout << "surf: " << surf.rows() << " " << surf.cols() << endl << "nods: " << nods.rows() << " " << nods.cols() << endl;

  surf.transposeInPlace();
  nods.transposeInPlace();
  MatrixXd points;
  gen_points(nods, surf, pt.get<size_t>("num_in_axis"), points);
  size_t dim = points.cols();
  
  cout << "[INFO]>>>>>>>>>>>>>>>>>>>POINTS<<<<<<<<<<<<<<<<<<" << endl;
  cout << points.block(0, 0, 3, dim > 10?10:dim) << endl;
  auto points_curr = points;
  //calc volume 
  double volume = clo_surf_vol(nods, surf);

  spatial_hash SH(points, 4);
  VectorXd sup_radi = SH.get_sup_radi();
  vector<vector<size_t>> friends_all(dim);

#pragma parallel omp for
  for(size_t i = 0; i < dim; ++i){
    SH.get_friends(points.col(i), sup_radi(i), friends_all[i]);
  }
  
  point_sys PS(points, pt.get<double>("rho"), pt.get<double>("Young"), pt.get<double>("Poission"), volume, 1, friends_all, sup_radi);

  energy_dat dat_str (dim);
  
  PS.pre_compute(points_curr.data(), dat_str);

  //construct deform_surf_MLS
  double kernel_cof = pt.get<double>("ker_cof");
  vector<vector<size_t>> vet_fris(nods.cols());
  for(size_t j = 0; j < nods.cols(); ++j){
    vector<size_t> fris;
    SH.get_friends(nods.col(j), kernel_cof, fris);
    vet_fris[j] = fris;
  }
  deform_surf_MLS<double> DS(surf, nods, points, vet_fris, kernel_cof);

  //TODO: make below a function
  //test spatial_hash::find_NN
  #if 0
  MatrixXi near_nei = MatrixXi::Zero(4, nods.cols());
  for(size_t i = 0; i < nods.cols(); ++i){
    
    auto ver_nn_num = 4;
    Vector2i two_index;
    two_index << 2, 3;
    bool is_co_line = true, is_co_plane = true;
    do{
      auto NN_index = SH.get_NN(nods.col(i), ver_nn_num);
      Vector3d V1 = (points.col(NN_index(1)) - points.col(NN_index(0))).normalized();
      Vector3d V2 = (points.col(NN_index(two_index(0))) - points.col(NN_index(0))).normalized();
      Vector3d cross_ = V1.cross(V2);
      if(fabs(cross_.norm()) < 1e-5){
        two_index(0)++;
        two_index(1)++;
        ver_nn_num ++;
        is_co_line = true;
      }
      else{
        is_co_line = false;
        do{
          Vector3d V3 = (points.col(NN_index(two_index(1))) - points.col(NN_index(0))).normalized();

          if(fabs(cross_.cross(V3).norm()) < 1e-5){
            two_index(1)++;
            ver_nn_num ++;
            NN_index = SH.get_NN(nods.col(i), ver_nn_num);
            // assert(ver_nn_num < 10);
            is_co_plane = true;
          }
          else{
            is_co_plane = false;
            near_nei(0, i) = NN_index(0);
            near_nei(1, i) = NN_index(1);
            near_nei(2, i) = NN_index(two_index(0));
            near_nei(3, i) = NN_index(two_index(1));
            break;
          }
        }while(is_co_plane);
      }
    }while(is_co_line);

  }
  cout << near_nei.block(0, 0, 4, 10);
#endif

  
  
  double delt_t = 0.05;
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
  size_t max_iter = 50;

  //add simple constraints
  vector<size_t> cons;
  for(size_t i = 0; i < points.cols(); ++i){
    if(points(2, i) > 2){
      cons.push_back(i);
      cout << i << " ";
    }
  }
  cout << endl;

  for(size_t i = 0; i < max_iter; ++i){
    dat_str.gra_.setZero(3, dim);
    cout << "[INFO]>>>>>>>>>>>>>>>>>>>disp<<<<<<<<<<<<<<<<<<" << endl;
    cout << "iter is "<<endl<< i << endl;
    // cout << "displace is "<<endl<< displace.block(0, 0, 3, 5) << endl;
    cout << "displace is " << endl<< displace.transpose().block(0, 0, 7, 3) << endl;
    cout << "velocity is "<<endl<< velocity.block(0, 0, 3, 5) << endl;

    PS.calc_defo_gra(displace.data(), dat_str);
    cout << "def_gra " << dat_str.def_gra_.block(0, 0, 9, 5) << endl;
    auto comp_4debug = dat_str.gra_;
    PS.Gra(displace.data(), dat_str);
    cout << "elasticity acce is " <<endl<< (dat_str.gra_ - comp_4debug).block(0, 0, 3, 5) << endl;

    // cout << "elasitic force " << dat_str.gra_.block(0, 0, 3, 5) << endl;    
    PS.gravity(displace.data(), dat_str, 9.8);
#pragma parallel omp for
    for(size_t j = 0; j < dim; ++j){
      assert(PS.get_mass(j) > 0);
      new_acce.col(j) = dat_str.gra_.col(j)/PS.get_mass(j);
    }
    velocity += 0.5*(new_acce + acce)*delt_t;
    for(auto c : cons){
      velocity.col(c) = MatrixXd::Zero(3, 1);
    }
    displace += velocity*delt_t + 0.5*acce*delt_t*delt_t;
    for(auto c : cons){
      displace.col(c) = MatrixXd::Zero(3, 1);
    }      

    vet_displace = DS.update_surf(displace, dat_str.def_gra_);
    cout << "vet displace " <<endl<< vet_displace.block(0, 0, 3, 5) << endl;
    acce = new_acce;

    dat_str.set_zero();
    // auto filename = pt.get<string>("res") + "_" + to_string(i) + ".obj";
    auto filename = pt.get<string>("res") + "_" + to_string(i) + ".vtk";
    // writeOBJ(filename, (nods + vet_displace).transpose(), surf.transpose());
    
    tri_mesh_write_to_vtk(filename.c_str(), nods + vet_displace, surf);
  }

  //done
  points += displace;
  


  


  // auto tmp = MatrixXd::Random(9, 1);
  // cout <<"tmp mat is "<< tmp << endl;
  // auto tmp1 = ele_mat(0, 3, 3, tmp);
  // cout <<"after function " <<   tmp1;


  // Map<MatrixXd> tmp (def_gra.col(0).data(), 3, 3);
  // cout << "this is tmp" << tmp<<endl << endl;
  // tmp(2, 2) = 999;
  // cout << "this is tmp" << tmp<<endl << endl;
  //c out <<"this is def_gra.col(0) : " << def_gra.col(0) << endl;
  
}


