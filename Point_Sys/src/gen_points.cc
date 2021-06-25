#include "gen_points.h"
#include "geometry.h"

#include <vector>
#include <iostream>
#include <cmath>
#include <WindingNumber/UT_Array.h>
#include <WindingNumber/UT_Array.cpp>
#include <WindingNumber/UT_SolidAngle.h>
#include <WindingNumber/UT_SolidAngle.cpp>
#include <Eigen/Core>
// typedef zjucad::matrix::matrix<size_t> mati_t;
// typedef zjucad::matrix::matrix<double> matd_t;
using namespace Eigen;
using namespace std;
using namespace HDK_Sample;



#define PI 3.14159265359
namespace marvel{



int get_inner_points(Eigen::MatrixXd &points, const Eigen::MatrixXi &surf, const Eigen::MatrixXd &nods, const bool save_surf_points){
  using UT_Vector3T = UT_FixedVector<float,3>;
  UT_Vector3T UT_nods[nods.cols()];

  MatrixXf nods_flo = nods.cast<float>();

// #pragma omp parallel for
  for(size_t i = 0; i < nods.cols(); ++i){
    UT_Vector3T vec_tmp;
    copy(nods_flo.col(i).data(), nods_flo.col(i).data() +  nods_flo.rows(), vec_tmp.vec);
    UT_nods[i] = vec_tmp;
  }

  UT_SolidAngle<float, float>  Comp_WN(int(surf.cols()), surf.data(), int(nods.cols()), UT_nods);
  vector<int> inside_id;
// #pragma omp parallel for
  for(size_t i = 0; i < points.cols(); ++i){
    UT_Vector3T query_point(&points(0, i));
    float sol_angle = Comp_WN.computeSolidAngle(query_point);

#pragma omp critical
    {
      if(!save_surf_points && fabs(sol_angle) > 1)
        inside_id.push_back(i);
      else if (save_surf_points && fabs(sol_angle - 4*PI) < 0.3)
        inside_id.push_back(i);

    }
  }
  
  MatrixXd points_tmp(points.rows(), inside_id.size());
// #pragma omp parallel for
  for(size_t i = 0; i < inside_id.size(); ++i){
    points_tmp.col(i) = points.col(inside_id[i]);
  }
  
  points = points_tmp;
  return 0;
}

int gen_points(const MatrixXd &nods, const MatrixXi &surf, const double interval, MatrixXd &points, const bool save_surf_points){
  //assert(num_in_axis > 1);
  MatrixXd bdbox;

  int res = build_bdbox(nods, bdbox);
  cout << "bdbox: " << bdbox << endl;
//   vector<double> intervals(3);

// // #pragma omp parallel for 
//   for(size_t i = 0; i < 3; ++i){
//     intervals[i] = (bdbox(i, 1) - bdbox(i, 0))/(num_in_axis - 1);
//   }

  vector<size_t> num_in_axis(3);
  for(size_t i = 0; i < 3; ++i){
    num_in_axis[i] = ceil((bdbox(i, 1) - bdbox(i, 0)) / interval);
  }

  points.setZero(3, num_in_axis[0]*num_in_axis[1]*num_in_axis[2]);
  const size_t num_in_plane = num_in_axis[1]*num_in_axis[2];
// #pragma omp parallel for 
  for(size_t i = 0; i < num_in_axis[0]; ++i){
    points.block(0, num_in_plane*i, 1, num_in_plane) = MatrixXd::Ones(1, num_in_plane)*(bdbox(0, 0)+ interval*i);
    for(size_t j = 0; j < num_in_axis[1]; ++j){
      points.block(1, num_in_plane*i + num_in_axis[2]*j, 1, num_in_axis[2]) = MatrixXd::Ones(1, num_in_axis[2])*(bdbox(1,0) + interval*j);
      for(size_t k = 0; k < num_in_axis[2]; ++k){
        points(2, num_in_plane*i + num_in_axis[2]*j + k) = bdbox(2, 0) + interval*k;
      }
    }
  }

  cout << "[INFO]genarate raw points done. size is " <<points.cols() << endl;
  res = get_inner_points(points, surf, nods, save_surf_points);
  cout << "after select: "  << points.cols() << endl;
  if (save_surf_points){
    MatrixXd points_inc_nods(3, nods.cols() + points.cols());
    points_inc_nods.block(0, 0, 3, nods.cols()) = nods;
    points_inc_nods.block(0, nods.cols(), 3, points.cols()) = points;
    points = points_inc_nods;
  }
  return 0;
}
int gen_points(const MatrixXd &nods, const MatrixXi &surf, const size_t &num_in_axis, MatrixXd &points, const bool save_surf_points){
  //assert(num_in_axis > 1);
  MatrixXd bdbox;

  int res = build_bdbox(nods, bdbox);
  cout << "bdbox: " << bdbox << endl;
  vector<double> intervals(3);

// #pragma omp parallel for 
  for(size_t i = 0; i < 3; ++i){
    intervals[i] = (bdbox(i, 1) - bdbox(i, 0))/(num_in_axis - 1);
  }

  points.setZero(3, num_in_axis*num_in_axis*num_in_axis);
  const size_t num_in_plane = num_in_axis*num_in_axis;
// #pragma omp parallel for 
  for(size_t i = 0; i < num_in_axis; ++i){
    points.block(0, num_in_plane*i, 1, num_in_plane) = MatrixXd::Ones(1, num_in_plane)*(bdbox(0, 0)+ intervals[0]*i);
    for(size_t j = 0; j < num_in_axis; ++j){
      points.block(1, num_in_plane*i + num_in_axis*j, 1, num_in_axis) = MatrixXd::Ones(1, num_in_axis)*(bdbox(1,0) + intervals[1]*j);
      for(size_t k = 0; k < num_in_axis; ++k){
        points(2, num_in_plane*i + num_in_axis*j + k) = bdbox(2, 0) + intervals[2]*k;
      }
    }
  }

  cout << "[INFO]genarate raw points done. size is " <<points.cols() << endl;
  res = get_inner_points(points, surf, nods, save_surf_points);
  cout << "after select: "  << points.cols() << endl;
  if (save_surf_points){
    MatrixXd points_inc_nods(3, nods.cols() + points.cols());
    points_inc_nods.block(0, 0, 3, nods.cols()) = nods;
    points_inc_nods.block(0, nods.cols(), 3, points.cols()) = points;
    points = points_inc_nods;
  }
  return 0;
}

  
}
