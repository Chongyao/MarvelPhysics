#ifndef GEOMETRY_UTIL_H
#define GEOMETRY_UTIL_H

#include <Eigen/Sparse>
#include <zjucad/matrix/matrix.h>
#include <jtflib/mesh/mesh.h>
#include <unordered_set>

using mati_t=zjucad::matrix::matrix<size_t>;
using matd_t=zjucad::matrix::matrix<double>;
using csc_t=Eigen::SparseMatrix<double>;

namespace bigbang {

void get_edge_elem(const mati_t &tris, mati_t &edge);

void get_diam_elem(const mati_t &tris, mati_t &diam);

// Here we assume that the qaud mesh has subdivisional
// structure, and i-th quad9 element is composed of small
// quad elements indexed from 4*i+0 to 4*i+3.
void get_quad9_elem(const mati_t &quad, mati_t &quad9);

/// Here we assume that the hexs mesh has subdivisional
/// structure, and i-th hex27 element is composed of
/// small cubes indexed from 8*i+0 to 8*i+7
void get_hex27_elem(const mati_t &hexs, mati_t &hex27);

double calc_tri_area(const matd_t &vert);

void eval_tri_rot(const double *x0, const double *x1, double *R);

void remove_extra_verts(mati_t &cell, matd_t &nods, mati_t *new2orig_mapping=0);

void extract_tet_surface(const mati_t &tets, const matd_t &nods, mati_t &tris, matd_t &snods);

void shrink_surface(const mati_t &tris, matd_t &nods, const double dist);

void subdivide_surface(const mati_t &tris, const matd_t &nods, mati_t &out_tris, matd_t &out_nods);

int interp_pts_in_tets(const matd_t &v, const mati_t &tet, const matd_t &pts, csc_t &coef);

void normalize_point_cloud(matd_t &v);

struct tri_mesh_subdivider
{
  tri_mesh_subdivider(const mati_t &tris, const matd_t &nods);
  int subdivide(mati_t &ftris, matd_t &fnods, const int sbd_sign);
  const size_t dim_;
  const mati_t &tris_;
  const matd_t &nods_;
  Eigen::SparseMatrix<double> R_, P_;
  size_t one_to_many_;
};

//        3---------2
//        |         |
//        |         |
//        |         |
//        |         |
//        0---------1
//  y
//  |__x
//
class straight_quad_strip_extractor
{
public:
  straight_quad_strip_extractor(const mati_t &quad);
  int extract(const size_t start_face, const char direction, std::vector<size_t> &face_strip);
private:
  const mati_t &quad_;
  std::shared_ptr<jtf::mesh::edge2cell_adjacent> e2c_;
};

//
//          4--------6
//         /|       /|
//        5--------7 |
//        | |      | |
//        | |      | |
//        | 0------|-2
//        |/       |/
//   z    1--------3
//   |__y
//  /
// x
class straight_cube_strip_extractor
{
public:
  straight_cube_strip_extractor(const mati_t &hexs);
  int extract(const size_t start_cell, const char direction, std::vector<size_t> &cube_strip);
private:
  const mati_t &hexs_;
  std::shared_ptr<jtf::mesh::face2hex_adjacent> f2h_;
};

int get_bnd_edges(const mati_t &cell, mati_t &bnds);
int get_surf_mesh_bnd_nodes(const mati_t &cell, mati_t &bnd_nodes);
int get_surf_mesh_bnd_nodes(const mati_t &cell, std::unordered_set<size_t> &bnd_nodes);

}
#endif
