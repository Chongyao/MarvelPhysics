#ifndef GEN_SURF_H
#define GEN_SURF_H

#include "get_nn.h"
#include <Eigen/Core>
namespace marvel{


class deform_surf{
 public:
  deform_surf(const Eigen::MatrixXi &surf, const Eigen::MatirxXd &nods, const Eigen::Matrix &sam_points);
 private:
  const Eigen::MatrixXi surf_;
  Eigen::MatrixXd nods_;
  const Eigen::MatrixXd sam_points_;
  spatial_hash SH_;
  
  
};//class deform surf



}//namespace: marvel
#endif
