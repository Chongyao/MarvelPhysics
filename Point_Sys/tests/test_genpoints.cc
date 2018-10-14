#include <string>

#include <zjucad/matrix/matrix.h>
#include <zjucad/ptree/ptree.h>
#include <jtflib/mesh/io.h>

#include "src/io.h"
#include "Point_Sys/src/init.h"
using namespace marval;
using namespace std;
typedef zjucad::matrix::matrix<size_t> mati_t;
typedef zjucad::matrix::matrix<double> matd_t;

int main(int argc, char** argv){

  boost::property_tree::ptree pt;
  zjucad::read_cmdline(argc, argv, pt);

  mati_t surf;
  matd_t nods;
  jtf::mesh::load_obj(pt.get<string>("surf.value").c_str(), surf, nods);

  matd_t points;
  get_inner_points(points, surf, nods);

  

  
  

  
}
