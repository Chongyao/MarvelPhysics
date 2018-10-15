#ifndef INIT_H
#define INIT_H
#include <zjucad/matrix/matrix.h>
#include <iostream>

typedef zjucad::matrix::matrix<size_t> mati_t;
typedef zjucad::matrix::matrix<double> matd_t;

namespace marval{
int build_bdbox(const matd_t &nods, matd_t & bdbox);
int get_inner_points(matd_t &points, const mati_t &surf, const matd_t &nods);
int gen_points(const matd_t &nods, const mati_t &surf, const size_t &num_in_axis, matd_t &points);

}
#endif
