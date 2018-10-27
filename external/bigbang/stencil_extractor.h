#ifndef STENCIL_EXTRACTOR_H
#define STENCIL_EXTRACTOR_H

#include <jtflib/mesh/mesh.h>

#include "stencil.h"

namespace bigbang {

typedef double (*large_to_small_mapping)(const int k, const double *eps, double *xi);
//        ___
//       |   |
//       |_0_|
//
inline double LS_2D_Q1_ID(const int k, const double *eps, double *xi) {
  if ( k == 0 ) return 1;
  if ( k == 1 ) {
    xi[0] = eps[0];
    xi[1] = eps[1];
  }
  return 0;
}
//        _______
//       |   |   |
//       |_3_|_2_|
//       |   |   |
//       |_0_|_1_|
//   
inline double LS_2D_Q4_E0(const int k, const double *eps, double *xi) {
  if ( k == 0 ) return 2;
  if ( k == 1 ) {
    xi[0] = (eps[0]-1)/2;
    xi[1] = (eps[1]-1)/2;
  }
  return 0;
}
inline double LS_2D_Q4_E1(const int k, const double *eps, double *xi) {
  if ( k == 0 ) return 2;
  if ( k == 1 ) {
    xi[0] = (eps[0]+1)/2;
    xi[1] = (eps[1]-1)/2;
  }
  return 0;
}
inline double LS_2D_Q4_E2(const int k, const double *eps, double *xi) {
  if ( k == 0 ) return 2;
  if ( k == 1 ) {
    xi[0] = (eps[0]+1)/2;
    xi[1] = (eps[1]+1)/2;
  }
  return 0;
}
inline double LS_2D_Q4_E3(const int k, const double *eps, double *xi) {
  if ( k == 0 ) return 2;
  if ( k == 1 ) {
    xi[0] = (eps[0]-1)/2;
    xi[1] = (eps[1]+1)/2;
  }
  return 0;
}
//        ___________
//       |   |   |   |
//       |_6_|_7_|_8_|
//       |   |   |   |
//       |_3_|_4_|_5_|
//       |   |   |   |
//       |_0_|_1_|_2_|
//
inline double LS_2D_Q9_E0(const int k, const double *eps, double *xi) {
  if ( k == 0 ) return 3;
  if ( k == 1 ) {
    xi[0] = (eps[0]-2)/3;
    xi[1] = (eps[1]-2)/3;
  }
  return 0;
}
inline double LS_2D_Q9_E1(const int k, const double *eps, double *xi) {
  if ( k == 0 ) return 3;
  if ( k == 1 ) {
    xi[0] = eps[0]/3;
    xi[1] = (eps[1]-2)/3;
  }
  return 0;
}
inline double LS_2D_Q9_E2(const int k, const double *eps, double *xi) {
  if ( k == 0 ) return 3;
  if ( k == 1 ) {
    xi[0] = (eps[0]+2)/3;
    xi[1] = (eps[1]-2)/3;
  }
  return 0;
}
inline double LS_2D_Q9_E3(const int k, const double *eps, double *xi) {
  if ( k == 0 ) return 3;
  if ( k == 1 ) {
    xi[0] = (eps[0]-2)/3;
    xi[1] = eps[1]/3;
  }
  return 0;
}
inline double LS_2D_Q9_E4(const int k, const double *eps, double *xi) {
  if ( k == 0 ) return 3;
  if ( k == 1 ) {
    xi[0] = eps[0]/3;
    xi[1] = eps[1]/3;
  }
  return 0;
}
inline double LS_2D_Q9_E5(const int k, const double *eps, double *xi) {
  if ( k == 0 ) return 3;
  if ( k == 1 ) {
    xi[0] = (eps[0]+2)/3;
    xi[1] = eps[1]/3;
  }
  return 0;
}
inline double LS_2D_Q9_E6(const int k, const double *eps, double *xi) {
  if ( k == 0 ) return 3;
  if ( k == 1 ) {
    xi[0] = (eps[0]-2)/3;
    xi[1] = (eps[1]+2)/3;
  }
  return 0;
}
inline double LS_2D_Q9_E7(const int k, const double *eps, double *xi) {
  if ( k == 0 ) return 3;
  if ( k == 1 ) {
    xi[0] = eps[0]/3;
    xi[1] = (eps[1]+2)/3;
  }
  return 0;
}
inline double LS_2D_Q9_E8(const int k, const double *eps, double *xi) {
  if ( k == 0 ) return 3;
  if ( k == 1 ) {
    xi[0] = (eps[0]+2)/3;
    xi[1] = (eps[1]+2)/3;
  }
  return 0;
}

class stencil_extractor
{
public:
  stencil_extractor(const mati_t &cell);
  int extract(std::vector<std::shared_ptr<quad_stencil>> &elem, const int max_radius) const;
  int extract_quad4_stencil(std::vector<std::shared_ptr<quad_stencil>> &elem) const;
  int extract_quad9_stencil(std::vector<std::shared_ptr<quad_stencil>> &elem) const;
  int extract_radius_2_quad9(std::vector<std::shared_ptr<quad_stencil>> &elem) const; //-> for debug use
private:
  const mati_t &cell_;
  mati_t ordered_curr_cell_;
};

class hex_stencil_extractor
{
public:
  hex_stencil_extractor(const mati_t &cell);
  int extract(std::vector<std::shared_ptr<hexs_stencil>> &elem, const int max_radius) const;
  int extract_hex8_stencil(std::vector<std::shared_ptr<hexs_stencil>> &elem) const;
  int extract_hex27_stencil(std::vector<std::shared_ptr<hexs_stencil>> &elem) const;
private:
  const mati_t &cell_;
};

class neighborhood_searcher
{
public:
  neighborhood_searcher(const mati_t &cell);
  int search(const size_t fid, const int max_radius, mati_t &neighbor,
             std::vector<size_t> &faces, large_to_small_mapping &l2s) const;
  bool search_r3(const size_t fid, mati_t &neigh, std::vector<size_t> &faces, large_to_small_mapping &l2s) const;
  bool search_r2(const size_t fid, mati_t &neigh, std::vector<size_t> &faces, large_to_small_mapping &l2s) const;
  bool search_r1(const size_t fid, mati_t &neigh, std::vector<size_t> &faces, large_to_small_mapping &l2s) const;
  bool is_center_of_quad9      (const size_t fid, mati_t &neigh, std::vector<size_t> &faces) const;
  bool is_top_left_of_quad4    (const size_t fid, mati_t &neigh, std::vector<size_t> &faces) const;
  bool is_top_right_of_quad4   (const size_t fid, mati_t &neigh, std::vector<size_t> &faces) const;
  bool is_bottom_left_of_quad4 (const size_t fid, mati_t &neigh, std::vector<size_t> &faces) const;
  bool is_bottom_right_of_quad4(const size_t fid, mati_t &neigh, std::vector<size_t> &faces) const;
protected:
  const mati_t &cell_;
  mati_t edges_;
  std::shared_ptr<jtf::mesh::edge2cell_adjacent> e2c_;
};

}

#endif
