#include <WindingNumber/UT_Array.h>
#include <WindingNumber/UT_Array.cpp>
#include <WindingNumber/UT_SolidAngle.h>
#include <WindingNumber/UT_SolidAngle.cpp>

#include <zjucad/matrix/matrix.h>
#include <zjucad/matrix/io.h>
#include <zjucad/ptree/ptree.h>
#include <jtflib/mesh/io.h>
#include "bigbang/io.h"

using namespace std;
using namespace bigbang;
using namespace zjucad::matrix;
using namespace HDK_Sample;
typedef zjucad::matrix::matrix<size_t> mati_t;
typedef zjucad::matrix::matrix<double> matd_t;

int main(int argc, char** argv){
  boost::property_tree::ptree pt;
  zjucad::read_cmdline(argc, argv, pt);

  mati_t surf;
  matd_t nods;
  jtf::mesh::load_obj(pt.get<string>("surf.value").c_str(), surf, nods);
  cout << "surf: " << surf.size(1) << " " << surf.size(2) << endl << "nods: " << nods.size(1) << " " << nods.size(2) << endl;

  matrix<float> nods_flo = nods;
  matrix<int> surf_int = surf;
  using UT_Vector3T = UT_FixedVector<float,3>;
  UT_Vector3T UT_nods[nods.size(2)];
#pragma omp parallel for
  for(size_t i = 0; i < nods.size(2); ++i){
    UT_Vector3T vec_tmp;
    copy(nods_flo(colon(), i).begin(),nods_flo(colon(), i).end(), vec_tmp.vec);
    UT_nods[i] = vec_tmp;
  }
  UT_SolidAngle<float, float>  Comp_WN(int(surf.size(2)), &surf_int(0, 0), int(nods.size(2)), UT_nods);

  matd_t points(3, 420);
  
  for(size_t i = 0; i < points.size(2); ++i){
    points(1, i) = 1.1-0.01*i;
    UT_Vector3T query_point(&points(0, i));
    double sol_angle = Comp_WN.computeSolidAngle(query_point);
    cout << points(1, i) << " solid angle: " << sol_angle << endl;
  }

  


  


  
}
