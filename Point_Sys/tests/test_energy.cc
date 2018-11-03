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

  
  cout << "[INFO]>>>>>>>>>>>>>>>>>>>POINTS<<<<<<<<<<<<<<<<<<" << endl;
  cout << points.block(0, 0, 3, points.cols() > 10?10:points.cols()) << endl;
  //calc volume 
  double volume = clo_surf_vol(nods, surf);
  
  point_sys PS(points, pt.get<double>("rho.value"), pt.get<double>("Young.value"), pt.get<double>("Poission.value"), volume, 4, 1);
  Matrix3d change;
  change << 1, 0, 0,
      0, 2, 0,
      0, 0, 5;
  MatrixXd points_curr = change*points;
  // MatrixXd def_gra(9, points.cols());
  // MatrixXd inv_A_all(9, points.cols());

  energy_dat dat_str (points.cols());
  PS.pre_compute(points_curr.data(), dat_str);
  
  double delt_t = 0.01;
  MatrixXd displace;
  MatrixXd velocity;
  MatrixXd acce;
  MatrixXd new_acce;
  MatrixXd gra;
  displace.setZero(3, points.cols());
  velocity.setZero(3, points.cols());
  acce.setZero(3, points.cols());
  gra.setZero(3, points.cols());
  new_acce.setZero(3, points.cols());
  size_t max_iter = 5;
  for(size_t i = 0; i < max_iter; ++i){
    
    cout << "[INFO]>>>>>>>>>>>>>>>>>>>disp<<<<<<<<<<<<<<<<<<" << endl;
    cout << "iter is "<<endl<< i << endl;
    cout << "displace is "<<endl<< displace.block(0, 0, 3, 5) << endl;
    cout << "velocity is "<<endl<< velocity.block(0, 0, 3, 5) << endl;
    cout << "acce is " <<endl<< acce.block(0, 0, 3, 5) << endl;

    PS.calc_defo_gra(displace.data(), dat_str);
    cout << "def_gra " << dat_str.def_gra_.block(0, 0, 9, 5) << endl;
    PS.Gra(displace.data(), dat_str);
    cout << "elasitic force " << dat_str.gra_.block(0, 0, 3, 5) << endl;    
    PS.gravity(displace.data(), dat_str, 9.8);
    cout << "gravity " << dat_str.gra_.block(0, 0, 3, 5) << endl;    
    for(size_t j = 0; j < points.cols(); ++j){
      assert(PS.get_mass(j) > 0);
      new_acce.col(j) = dat_str.gra_.col(j)/PS.get_mass(j);
    }

    displace += velocity*delt_t + 0.5*acce*delt_t*delt_t;
    velocity += 0.5*(new_acce + acce)*delt_t;
    acce = new_acce;

    dat_str.set_zero();
  }

  //done
  points += displace;
  
  // writeOBJ(pt.get<string>("surf.value").c_str(), nods, surf);

  


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


