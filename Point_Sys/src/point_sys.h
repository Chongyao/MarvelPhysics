#ifndef POINT_STS_H
#define POINT_STS_H
#include <eigen3/Eigen/Core>


namespace marvel{


class Point_system{
 public:
  Point_system();
 private:
  Eigen::MatrixXd points;
  size_t num;
  
};


}
#endif
