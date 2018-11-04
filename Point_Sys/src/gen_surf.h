#ifndef GEN_SURF_H
#define GEN_SURF_H

#include "get_nn.h"
#include <Eigen/Core>
namespace marvel{


class deform_surf{
 public:
  //TODO: do not use copy construct!!
  deform_surf(const Eigen::MatrixXi &surf, const Eigen::MatirxXd &nods, const Eigen::Matrix &sam_points, const spatial_hash &SH, const double );
 private:
  const Eigen::MatrixXi surf_;
  Eigen::MatrixXd nods_;
  const Eigen::MatrixXd sam_points_;
  spatial_hash SH_;
  
  
};//class deform surf



}//namespace: marvel
#endif
