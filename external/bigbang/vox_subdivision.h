#ifndef VOX_SUBDIVISION_H
#define VOX_SUBDIVISION_H

#include <map>
#include <zjucad/matrix/matrix.h>
#include <Eigen/Sparse>

using mati_t=zjucad::matrix::matrix<size_t>;
using matd_t=zjucad::matrix::matrix<double>;

namespace bigbang {

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


//        3---------2
//        |         |
//        |         |
//        |         |
//        |         |
//        0---------1
//  y
//  |__x
//

void collect_one_quad_edges(const mati_t &quad,
                            std::vector<std::vector<size_t>> &edges);

void collect_one_vox_edges(const mati_t &vox,
                           std::vector<std::vector<size_t>> &edges);

void collect_one_vox_faces(const mati_t &vox,
                           std::vector<std::vector<size_t>> &faces);

class quad_edges
{
public:
  quad_edges(const mati_t &cell);
  size_t edges_num() const;
  size_t query_edge_idx(const size_t p, const size_t q) const;
  size_t query_edge_idx(const std::vector<size_t> &ab) const;
  std::Eigen::Map<std::vector<size_t>, size_t> edge2idx_;
};

class vox_edges
{
public:
  vox_edges(const mati_t &cell);
  size_t edges_num() const;
  size_t query_edge_idx(const size_t p, const size_t q) const;
  size_t query_edge_idx(const std::vector<size_t> &ab) const;
  std::Eigen::Map<std::vector<size_t>, size_t> edge2idx_;
};

class vox_faces
{
public:
  vox_faces(const mati_t &cell);
  size_t faces_num() const;
  size_t query_face_idx(const size_t p, const size_t q, const size_t m, const size_t n) const;
  size_t query_face_idx(const std::vector<size_t> &abcd) const;
  std::Eigen::Map<std::vector<size_t>, size_t> face2idx_;
  std::Eigen::Map<std::vector<size_t>, std::vector<size_t>> face2vox_;
};

/* 3D cube */
int subdivide_vox(const mati_t &voxs, const matd_t &nods,
                  mati_t &new_voxs,   matd_t &new_nods,
                  Eigen::SparseMatrix<double> *P=nullptr);

/* 2D quad */
int subdivide_quad(const mati_t &quad, const matd_t &nods,
                   mati_t &new_quad, matd_t &new_nods,
                   Eigen::SparseMatrix<double> *P=nullptr);

}

#endif
