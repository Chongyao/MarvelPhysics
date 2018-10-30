#include <string>
#include <iostream>
#include <zjucad/ptree/ptree.h>

#include <Eigen/Core>
#include "Point_Sys/src/gen_points.h"
#include <libigl/include/igl/readOBJ.h>
#include "Point_Sys/src/get_nn.h"

#include <chrono>

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

  



  cout << "spatial hash " <<endl;
  MatrixXi NN_;
  VectorXd sup_radii_;
  auto start = system_clock::now();
  spatial_hash SH(points, 10);
  auto end = system_clock::now();
  auto duration = duration_cast<microseconds>(end - start);
  cout <<  "花费了" 
       << double(duration.count()) * microseconds::period::num / microseconds::period::den 
       << "秒" << endl;
  

  cout << "all done " << endl;
                                                                                                                        

  
}

