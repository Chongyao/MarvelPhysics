#include <string>
#include <iostream>

// #include <zjucad/matrix/matrix.h>
// #include <zjucad/matrix/io.h>
// #include <zjucad/ptree/ptree.h>
// #include <jtflib/mesh/io.h>  

// #include "bigbang/io.h"
#include "Point_Sys/src/gen_points.h"
using namespace marvel;
using namespace std;
using namespace bigbang;
using namespace zjucad::matrix;
typedef zjucad::matrix::matrix<size_t> mati_t;
typedef zjucad::matrix::matrix<double> matd_t;

int main(int argc, char** argv){

  boost::property_tree::ptree pt;
  zjucad::read_cmdline(argc, argv, pt);

  mati_t surf;
  matd_t nods;
  jtf::mesh::load_obj(pt.get<string>("surf.value").c_str(), surf, nods);
  cout << "surf: " << surf.size(1) << " " << surf.size(2) << endl << "nods: " << nods.size(1) << " " << nods.size(2) << endl;

                                                                                                                        
  matd_t points;
  gen_points(nods, surf, pt.get<size_t>("num_in_axis.value"), points);
  cout << pt.get<string>("points_out.value").c_str() << endl;
  // cout << "points" << points<<endl;
  point_write_to_vtk(pt.get<string>("points_out.value").c_str(), points);
  cout << "all done " << endl;
                                                                                                                        

  
}

