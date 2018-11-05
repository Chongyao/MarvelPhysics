#include <string>
#include <iostream>
#include <chrono>

#include <zjucad/ptree/ptree.h>
#include <libigl/include/igl/readOBJ.h>
#include <libigl/include/igl/writeOBJ.h>
#include <Eigen/Core>

#include "Point_Sys/src/geometry.h"
#include "Point_Sys/src/gen_points.h"
#include "Point_Sys/src/get_nn.h"
#include "Point_Sys/src/points_energy.h"
#include "Point_Sys/src/data_stream.h"
#include "Point_Sys/src/gen_surf.h"



using namespace marvel;
using namespace std;
using namespace Eigen;
using namespace igl;
using namespace chrono;

// const MatrixXd& ele_mat(const size_t &ele_id, const size_t &rows, const size_t &cols, MatrixXd &&mat){
//   return Map<MatrixXd>(mat.col(ele_id).data(), rows, cols);
// }

int main(int argc, char** argv){

  boost::property_tree::ptree pt;
  zjucad::read_cmdline(argc, argv, pt);

  MatrixXi surf;
  MatrixXd nods;
  // jtf::mesh::load_obj(pt.get<string>("surf.value").c_str(), surf, nods);
  readOBJ(pt.get<string>("surf.value").c_str(), nods, surf);
  cout << "surf: " << surf.rows() << " " << surf.cols() << endl << "nods: " << nods.rows() << " " << nods.cols() << endl;

  surf.transposeInPlace();
  nods.transposeInPlace();
  MatrixXd points;
  gen_points(nods, surf, pt.get<size_t>("num_in_axis.value"), points);
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
  
  point_sys PS(points, pt.get<double>("rho.value"), pt.get<double>("Young.value"), pt.get<double>("Poission.value"), volume, 1, friends_all, sup_radi);

  energy_dat dat_str (dim);
  PS.pre_compute(points_curr.data(), dat_str);

  //construct deform_surf_MLS
  double kernel_cof = pt.get<double>("ker_cof.value");
  vector<vector<size_t>> vet_fris(nods.cols());
  for(size_t j = 0; j < nods.cols(); ++j){
    vector<size_t> fris;
    SH.get_friends(nods.col(j), kernel_cof, fris);
    vet_fris[j] = fris;
  }
  deform_surf_MLS<double> DS(surf, nods, points, vet_fris, kernel_cof);
  
  
  double delt_t = 0.01;
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
  size_t max_iter = 5;
  

  for(size_t i = 0; i < max_iter; ++i){
    
    cout << "[INFO]>>>>>>>>>>>>>>>>>>>disp<<<<<<<<<<<<<<<<<<" << endl;
    cout << "iter is "<<endl<< i << endl;
    cout << "displace is "<<endl<< displace.block(0, 0, 3, 5) << endl;
    cout << "velocity is "<<endl<< velocity.block(0, 0, 3, 5) << endl;
    cout << "acce is " <<endl<< acce.block(0, 0, 3, 5) << endl;

    PS.calc_defo_gra(displace.data(), dat_str);
    // cout << "def_gra " << dat_str.def_gra_.block(0, 0, 9, 5) << endl;
    PS.Gra(displace.data(), dat_str);
    // cout << "elasitic force " << dat_str.gra_.block(0, 0, 3, 5) << endl;    
    PS.gravity(displace.data(), dat_str, 9.8);
    for(size_t j = 0; j < dim; ++j){
      assert(PS.get_mass(j) > 0);
      new_acce.col(j) = dat_str.gra_.col(j)/PS.get_mass(j);
    }
    
    displace += velocity*delt_t + 0.5*acce*delt_t*delt_t;
    vet_displace = DS.update_surf(displace, dat_str.def_gra_);
    velocity += 0.5*(new_acce + acce)*delt_t;
    acce = new_acce;

    dat_str.set_zero();

    writeOBJ(pt.get<string>("res.value").c_str() + to_string(i) + ".obj", nods + vet_displace, surf);    
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
  // cout <<"this is def_gra.col(0) : " << def_gra.col(0) << endl;
  
}


