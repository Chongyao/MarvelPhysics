#include <string>
#include <iostream>
#include <zjucad/ptree/ptree.h>

#include <Eigen/Core>
#include "Point_Sys/src/gen_points.h"
#include <libigl/include/igl/readOBJ.h>
using namespace marvel;
using namespace std;
using namespace Eigen;
using namespace igl;


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
  cout << pt.get<string>("points_out.value").c_str() << endl;
  // point_write_to_vtk(pt.get<string>("points_out.value").c_str(), points);
  cout << "all done " << endl;
                                                                                                                        

  
}

