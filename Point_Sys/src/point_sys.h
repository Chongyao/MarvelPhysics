#ifndef POINT_STS_H
#define POINT_STS_H
#include <zjucad/martix/martix.h>

typedef zjucad::matrix::matrix<size_t> mati_t;
typedef zjucad::matrix::matrix<double> matd_t;

namespace marval{


class Point_system{
 public:
  Point_system();
 private:
  matd_t points;
  size_t num;
  
};


}
#endif
