#include "init.h"
#include <vector>
#include <WindingNumber/UT_Array.h>
#include <WindingNumber/UT_Array.cpp>
#include <WindingNumber/UT_SolidAngle.h>
#include <WindingNumber/UT_SolidAngle.cpp>
typedef zjucad::matrix::matrix<size_t> mati_t;
typedef zjucad::matrix::matrix<double> matd_t;
using namespace zjucad::matrix;
using namespace std;
using namespace HDK_Sample;

#define PI 3.14159265359
namespace marval{


int build_bdbox(const matd_t &nods, matd_t & bdbox){
  //simple bounding box
  bdbox = nods(colon(), 0)*ones<double>(1, 2);
  
  for(size_t i = 0; i < nods.size(2); ++i){
    for(size_t j = 0; j < nods.size(1); ++j){
      if(bdbox(j, 0) > nods(j, i))
        bdbox(j, 0) = nods(j ,i);
      if(bdbox(j, 1) < nods(j, i))
        bdbox(j, 1) = nods(j ,i);      
    }
  }

  return 0;
}
int get_inner_points(matd_t &points, const mati_t &surf, const matd_t &nods){
  // #include <WindingNumber/UT_SolidAngle.cpp>
  using UT_Vector3T = UT_FixedVector<float,3>;
  UT_Vector3T * UT_nods;

  matrix<float> nods_flo = nods;
  matrix<int> surf_int = surf;
#pragma omp parallel for 
  for(size_t i = 0; i < nods.size(2); ++i){
    UT_nods[i] = UT_Vector3T(&nods_flo(0, i));
  }
  
  UT_SolidAngle<float, float>  Comp_WN(int(surf.size(2)), &surf_int(0, 0), int(nods.size(2)), UT_nods);
  vector<int> inside_id;
#pragma omp parallel for
  for(size_t i = 0; i < points.size(2); ++i){
    UT_Vector3T query_point(&points(0, i));
    float sol_angle = Comp_WN.computeSolidAngle(query_point);

    if (sol_angle/(4*PI) - 1 < 1e-5)
      inside_id.push_back(i);
  }

  mati_t inside_mat(inside_id.size(), 1);
  points = points(colon(), inside_mat(colon(), 0));
  
  return 0;
}

int gen_points(const matd_t &nods, const mati_t &surf, const size_t &num_in_axis, matd_t points){
  assert(num_in_axis > 1);
  
  matd_t bdbox;
  int res = build_bdbox(nods, bdbox);
  vector<double> intervals(3);
#pragma omp parallel for 
  for(size_t i = 0; i < 3; ++i){
    intervals[i] = (bdbox(i, 1) - bdbox(i, 0))/(num_in_axis - 1);
  }

  points.resize(3, num_in_axis*num_in_axis*num_in_axis);
  size_t num_in_plane = num_in_axis*num_in_axis;
#pragma omp parallel for 
  for(size_t i = 0; i < num_in_axis; ++i){
    points(0, colon(num_in_plane*i, num_in_plane*(i + 1))) = ones<double>(1, num_in_plane)*(bdbox(0, 0)+ intervals[0]*i);
    for(size_t j = 0; j < num_in_axis; ++j){
      points(1, colon(num_in_plane*i + num_in_axis*j, num_in_plane + num_in_axis*(j + 1))) = ones<double>(1, num_in_axis)*(bdbox(1,0) + intervals[1]*j);
      for(size_t k = 0; k < num_in_axis; ++k){
        points(2, num_in_plane*i + num_in_axis*j + k) = bdbox(2, 0) + intervals[2]*k;
      }
    }
  }

  cout << "raw points "<< points.size(2) << endl;

  res = get_inner_points(points, surf, nods);
}

}
