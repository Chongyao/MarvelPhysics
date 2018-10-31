#include "geometry.h"

#include <Eigen/LU>

using namespace Eigen;
namespace marvel{
double clo_surf_vol(const MatrixXd &nods, const MatrixXi &surf){
  //TODO:check if the surface is closed and manifold
  double volume = 0;
  for(size_t i = 0; i < surf.cols(); ++i){
    Matrix3d tet;
    for(size_t j = 0; j < 3; ++j){
      tet.row(j) = nods.col(surf(j, i));
    }
    //TODO:check
    volume += tet.determinant();
  }

  return volume;

  
}
}//namespace : marvel
