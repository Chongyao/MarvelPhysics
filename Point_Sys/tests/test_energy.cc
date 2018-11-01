#include <string>
#include <iostream>
#include <chrono>

#include <zjucad/ptree/ptree.h>
#include <libigl/include/igl/readOBJ.h>
#include <Eigen/Core>

#include "Point_Sys/src/geometry.h"
#include "Point_Sys/src/gen_points.h"
#include "Point_Sys/src/get_nn.h"
#include "Point_Sys/src/points_energy.h"




using namespace marvel;
using namespace std;
using namespace Eigen;
using namespace igl;
using namespace chrono;

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
  
  point_sys PS(points, pt.get<double>("rho.value"), volume, 4);
  MatrixXd points_curr = points.array()*0.5;
  MatrixXd def_gra(9, points.cols());
  
  PS.calc_defo_gra(points_curr.data(), def_gra.data());
  cout << def_gra.block(0, 0, 9, 10) << endl;
}


