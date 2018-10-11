#include "init.h"
#include <vector>

typedef zjucad::matrix::matrix<size_t> mati_t;
typedef zjucad::matrix::matrix<double> matd_t;
using namespace zjucad::matrix;
using namespace std;

namespace marval{


int build_bdbox(const matd_t &nods, matd_t & bdbox){
  //simple bounding box
  mdbox = nods(colon(), 0)*ones<double>(1, 2);
  
  for(size_t i = 0; i < nods.size(1); ++i){
    for(size_t j = 0; j < nods.size(0); ++j){
      if(bdbox(j, 0) > nods(j, i))
        bdbox(j, 0) = nods(j ,i);
      if(bdbox(j, 1) < nods(j, i))
        bdbox(j, 1) = nods(j ,i);      
    }
  }
}


int gen_points(const matd_t &nods, const mati_t &surf, const size_t &num_in_axis, madt_t points){
  assert(num_in_axis > 1);
  
  matd_t bdbox;
  build_bdbox(nods, bdbox);
  vector<double> intervals(3);
  for(size_t i = 0; i < 3; ++i){
    intervals[i] = (bdbox(i, 1) - bdbox(i, 0))/(num_in_axis - 1);
  }

  points.resize(3, num_in_axis*num_in_axis*num_in_axis);
  size_t num_in_plane = num_in_axis*num_in_axis;
  for(size_t i = 0; i < num_in_axis; ++i){
    points(0, colon(num_in_plane*i, num_in_plane*(i + 1))) = ones<double>(1, num_in_plane)*(bdbox(0, 0)+ intervals[0]*i);
    for(size_t j = 0; j < num_in_axis; ++j){
      points(1, colon(num_in_plane*i + num_in_axis*j, num_in_plane + num_in_axis*(j + 1))) = ones<double>(1, num_in_axis)*(bdbox(1,0) + intervals[1]*j);
      for(size_t k = 0; k < num_in_axis; ++k){
        points(2, num_in_plane*i + num_in_axis*j + k) = bdbox(2, 0) + intervals[2]*k;
      }
    }
  }
}

}

