#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <Eigen/Core>

namespace marvel{

double clo_surf_vol(const Eigen::MatrixXd &nods, const Eigen::MatrixXi &surf);

}//namespcae : marvel
#endif
