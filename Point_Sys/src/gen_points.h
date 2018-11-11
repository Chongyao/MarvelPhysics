#ifndef INIT_H
#define INIT_H
#include <iostream>

#include <Eigen/Core>


namespace marvel{
int build_bdbox(const Eigen::MatrixXd &nods, Eigen::MatrixXd & bdbox);
int get_inner_points(Eigen::MatrixXd &points, const Eigen::MatrixXi &surf, const Eigen::MatrixXd &nods);
int gen_points(const Eigen::MatrixXd &nods, const Eigen::MatrixXi &surf, const size_t &num_in_axis, Eigen::MatrixXd &points);
}
#endif

