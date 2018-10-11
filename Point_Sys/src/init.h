#ifndef INIT_H
#define INIT_H
#include <zjcuad/martix/martix.h>
#include <iostream>

typedef zjucad::matrix::matrix<size_t> mati_t;
typedef zjucad::matrix::matrix<double> matd_t;

namespace marval{



int build_bdbox(const matd_t &nods, matd_t & bdbox);
int gen_points(const matd_t &nods, const mati_t &surf);
#endif

}
