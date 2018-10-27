#include "stencil_extractor.h"

#include <zjucad/matrix/io.h>

#include "config.h"
#include "util.h"
#include "geom_util.h"

using namespace std;
using namespace zjucad::matrix;

namespace bigbang {

//===============================================================================
/**
 *  Assumption: Each quad has a consistent order on the index permuatation:
 *  all start from the "bottom left" vertices with a counterclockwise order.
 *
 *  3__2
 *  |__|
 *  0  1
 *
 **/
//===============================================================================
neighborhood_searcher::neighborhood_searcher(const mati_t &cell)
    : cell_(cell) {
  //-> for quad mesh only, extract all edges
  edges_.resize(2, 4*cell_.size(2));
  for (size_t i = 0; i < cell_.size(2); ++i) {
    edges_(0, 4*i+0) = cell_(0, i);
    edges_(1, 4*i+0) = cell_(1, i);

    edges_(0, 4*i+1) = cell_(1, i);
    edges_(1, 4*i+1) = cell_(2, i);

    edges_(0, 4*i+2) = cell_(2, i);
    edges_(1, 4*i+2) = cell_(3, i);

    edges_(0, 4*i+3) = cell_(3, i);
    edges_(1, 4*i+3) = cell_(0, i);
  }
  e2c_.reset(jtf::mesh::edge2cell_adjacent::create(cell_));
}

int neighborhood_searcher::search(const size_t fid, const int max_radius, mati_t &neigh, vector<size_t> &faces, large_to_small_mapping &l2s) const {
  if ( max_radius == 3 ) {
    if ( search_r3(fid, neigh, faces, l2s) )
      return 3;

    if ( search_r2(fid, neigh, faces, l2s) )
      return 2;

    if ( search_r1(fid, neigh, faces, l2s) )
      return 1;
  }

  if ( max_radius == 2 ) {
    if ( search_r2(fid, neigh, faces, l2s) )
      return 2;

    if ( search_r1(fid, neigh, faces, l2s) )
      return 1;
  }

  if ( max_radius == 1 ) {
    if ( search_r1(fid, neigh, faces, l2s) )      
      return 1;    
  }

  return -1;
}

bool neighborhood_searcher::search_r1(const size_t fid, mati_t &neigh, vector<size_t> &faces, large_to_small_mapping &l2s) const {
  mati_t idx_perm(4, 1); {
    idx_perm[0] = 0;
    idx_perm[1] = 3;
    idx_perm[2] = 1;
    idx_perm[3] = 2;
  }
  neigh = cell_(colon(), fid)(idx_perm);
  l2s = LS_2D_Q1_ID;

  faces.resize(1);
  faces[0] = fid;

  return true;
}

bool neighborhood_searcher::search_r2(const size_t fid, mati_t &neigh, vector<size_t> &faces, large_to_small_mapping &l2s) const {
  // 1 category, 4 cases
  if ( is_top_left_of_quad4(fid, neigh, faces) ) {
    l2s = LS_2D_Q4_E3;
    return true;
  }

  if ( is_top_right_of_quad4(fid, neigh, faces) ) {
    l2s = LS_2D_Q4_E2;
    return true;
  }

  if ( is_bottom_left_of_quad4(fid, neigh, faces) ) {
    l2s = LS_2D_Q4_E0;
    return true;
  }

  if ( is_bottom_right_of_quad4(fid, neigh, faces) ) {
    l2s = LS_2D_Q4_E1;
    return true;
  }
  
  return false;
}

bool neighborhood_searcher::search_r3(const size_t fid, mati_t &neigh, vector<size_t> &faces, large_to_small_mapping &l2s) const {
  // 3 categories, 9 cases
  if ( is_center_of_quad9(fid, neigh, faces) ) {
    l2s = LS_2D_Q9_E4;
    return true;
  }

  for (size_t i = 0; i < 4; ++i) {
    const pair<size_t, size_t> adjf = e2c_->query(edges_(0, 4*fid+i), edges_(1, 4*fid+i));
    if ( e2c_->is_boundary_edge(adjf) )
      continue;
    const size_t next_face = adjf.first+adjf.second-fid;
    if ( is_center_of_quad9(next_face, neigh, faces) ) {
      if ( i == 0 )
        l2s = LS_2D_Q9_E7;
      else if ( i == 1 )
        l2s = LS_2D_Q9_E3;
      else if ( i == 2 )
        l2s = LS_2D_Q9_E1;
      else
        l2s = LS_2D_Q9_E5;
      return true;
    }
  }  
  
  return false;
}

bool neighborhood_searcher::is_center_of_quad9(const size_t fid, mati_t &neigh, vector<size_t> &adjc_face) const {
  //   P  U  N
  //   L  f  R
  //   Q  D  M
  pair<size_t, size_t> adjf;

  size_t face_DRUL[4];
  for (size_t i = 0; i < 4; ++i) {
    adjf = e2c_->query(edges_(0, 4*fid+i), edges_(1, 4*fid+i));
    if ( e2c_->is_boundary_edge(adjf) )
      return false;
    face_DRUL[i] = adjf.first+adjf.second-fid;
  }

  size_t face_MNPQ[4];
  for (size_t i = 0; i < 4; ++i) {
    const size_t fid = face_DRUL[i];
    adjf = e2c_->query(edges_(0, 4*fid+(i+1)%4), edges_(1, 4*fid+(i+1)%4));
    if ( e2c_->is_boundary_edge(adjf) )
      return false;
    face_MNPQ[i] = adjf.first+adjf.second-fid;
  }

  const size_t fD = face_DRUL[0], fR = face_DRUL[1], fU = face_DRUL[2], fL = face_DRUL[3];
  const size_t fM = face_MNPQ[0], fN = face_MNPQ[1], fP = face_MNPQ[2], fQ = face_MNPQ[3];

  adjc_face.resize(9);
  adjc_face[0] = fQ;
  adjc_face[1] = fL;
  adjc_face[2] = fP;
  adjc_face[3] = fD;
  adjc_face[4] = fid;
  adjc_face[5] = fU;
  adjc_face[6] = fM;
  adjc_face[7] = fR;
  adjc_face[8] = fN;

  neigh.resize(16, 1);
  neigh[0] = cell_(0, fQ);
  neigh[1] = cell_(0, fL);
  neigh[2] = cell_(0, fP);
  neigh[3] = cell_(3, fP);

  neigh[4] = cell_(0, fD);
  neigh[5] = cell_(0, fid);
  neigh[6] = cell_(0, fU);
  neigh[7] = cell_(3, fU);

  neigh[8] = cell_(0, fM);
  neigh[9] = cell_(0, fR);
  neigh[10] = cell_(0, fN);
  neigh[11] = cell_(3, fN);

  neigh[12] = cell_(1, fM);
  neigh[13] = cell_(1, fR);
  neigh[14] = cell_(1, fN);
  neigh[15] = cell_(2, fN);

  return true;
}

#define ASSIGN_QUAD4_INDEX(A, B, C, D) \
  do {                                          \
    neigh[0] = cell_(0, A);                     \
    neigh[1] = cell_(0, B);                     \
    neigh[2] = cell_(3, B);                     \
    neigh[3] = cell_(0, C);                     \
    neigh[4] = cell_(0, D);                     \
    neigh[5] = cell_(3, D);                     \
    neigh[6] = cell_(1, C);                     \
    neigh[7] = cell_(1, D);                     \
    neigh[8] = cell_(2, D);                     \
  } while (0);

bool neighborhood_searcher::is_top_left_of_quad4(const size_t fid, mati_t &neigh, vector<size_t> &adjc_face) const {
  //  f  1
  //  0  2
  pair<size_t, size_t> adjf;

  adjf = e2c_->query(&edges_(0, 4*fid+0));
  if ( e2c_->is_boundary_edge(adjf) )
    return false;
  const size_t f0 = adjf.first+adjf.second-fid;

  adjf = e2c_->query(&edges_(0, 4*fid+1));
  if ( e2c_->is_boundary_edge(adjf) )
    return false;
  const size_t f1 = adjf.first+adjf.second-fid;

  adjf = e2c_->query(&edges_(0, 4*f0+1));
  if ( e2c_->is_boundary_edge(adjf) )
    return false;
  const size_t f2 = adjf.first+adjf.second-f0;

  adjc_face.resize(4);
  adjc_face[0] = f0;
  adjc_face[1] = fid;
  adjc_face[2] = f2;
  adjc_face[3] = f1;
  
  neigh.resize(9, 1);
  ASSIGN_QUAD4_INDEX(f0, fid, f2, f1);

  return true;
}

bool neighborhood_searcher::is_bottom_left_of_quad4(const size_t fid, mati_t &neigh, vector<size_t> &adjc_face) const {
  //  1  2
  //  f  0
  pair<size_t, size_t> adjf;

  adjf = e2c_->query(&edges_(0, 4*fid+1));
  if ( e2c_->is_boundary_edge(adjf) )
    return false;
  const size_t f0 = adjf.first+adjf.second-fid;

  adjf = e2c_->query(&edges_(0, 4*fid+2));
  if ( e2c_->is_boundary_edge(adjf) )
    return false;
  const size_t f1 = adjf.first+adjf.second-fid;

  adjf = e2c_->query(&edges_(0, 4*f0+2));
  if ( e2c_->is_boundary_edge(adjf) )
    return false;
  const size_t f2 = adjf.first+adjf.second-f0;

  adjc_face.resize(4);
  adjc_face[0] = fid;
  adjc_face[1] = f1;
  adjc_face[2] = f0;
  adjc_face[3] = f2;
  
  neigh.resize(9, 1);
  ASSIGN_QUAD4_INDEX(fid, f1, f0, f2);

  return true;
}

bool neighborhood_searcher::is_bottom_right_of_quad4(const size_t fid, mati_t &neigh, vector<size_t> &adjc_face) const {
  //  2  0
  //  1  f
  pair<size_t, size_t> adjf;

  adjf = e2c_->query(&edges_(0, 4*fid+2));
  if ( e2c_->is_boundary_edge(adjf) )
    return false;
  const size_t f0 = adjf.first+adjf.second-fid;
  
  adjf = e2c_->query(&edges_(0, 4*fid+3));
  if ( e2c_->is_boundary_edge(adjf) )
    return false;
  const size_t f1 = adjf.first+adjf.second-fid;

  adjf = e2c_->query(&edges_(0, 4*f0+3));
  if ( e2c_->is_boundary_edge(adjf) )
    return false;
  const size_t f2 = adjf.first+adjf.second-f0;

  adjc_face.resize(4);
  adjc_face[0] = f1;
  adjc_face[1] = f2;
  adjc_face[2] = fid;
  adjc_face[3] = f0;

  neigh.resize(9, 1);
  ASSIGN_QUAD4_INDEX(f1, f2, fid, f0);

  return true;
}

bool neighborhood_searcher::is_top_right_of_quad4(const size_t fid, mati_t &neigh, vector<size_t> &adjc_face) const {
  //  0  f
  //  2  1
  pair<size_t, size_t> adjf;

  adjf = e2c_->query(&edges_(0, 4*fid+3));
  if ( e2c_->is_boundary_edge(adjf) )
    return false;
  const size_t f0 = adjf.first+adjf.second-fid;
  
  adjf = e2c_->query(&edges_(0, 4*fid+0));
  if ( e2c_->is_boundary_edge(adjf) )
    return false;
  const size_t f1 = adjf.first+adjf.second-fid;

  adjf = e2c_->query(&edges_(0, 4*f0+0));
  if ( e2c_->is_boundary_edge(adjf) )
    return false;
  const size_t f2 = adjf.first+adjf.second-f0;

  adjc_face.resize(4);
  adjc_face[0] = f2;
  adjc_face[1] = f0;
  adjc_face[2] = f1;
  adjc_face[3] = fid;

  neigh.resize(9, 1);
  ASSIGN_QUAD4_INDEX(f2, f0, f1, fid);

  return true;
}

//===============================================================================
stencil_extractor::stencil_extractor(const mati_t &cell)
    : cell_(cell) {
  // reorder the vertex index in order to map the correct shape functions
  ordered_curr_cell_.resize(cell.size(1), cell.size(2));
  for (size_t i = 0; i < cell.size(2); ++i) {    
    mati_t idx_perm(4, 1); {
      idx_perm[0] = 0;
      idx_perm[1] = 3;
      idx_perm[2] = 1;
      idx_perm[3] = 2;
    }
    ordered_curr_cell_(colon(), i) = cell(colon(), i)(idx_perm);
  }
}

int stencil_extractor::extract(vector<shared_ptr<quad_stencil>> &elem, const int max_radius) const {
  // ASSERT(max_radius >= 1 && max_radius <= 3);
  // elem.resize(cell_.size(2));

  // neighborhood_searcher handle(cell_);

  // mati_t nei_num_count = zeros<size_t>(3, 1);
  // for (size_t i = 0; i < elem.size(); ++i)
  // {
  //   mati_t ordered_adjc_cell;
  //   large_to_small_mapping l2s;
  //   vector<size_t> adjc_face;
  //   int nei_num = handle.search(i, max_radius, ordered_adjc_cell, adjc_face, l2s);
  //   ASSERT(nei_num != -1);
  //   ++nei_num_count[nei_num-1];

  //   elem[i] = make_shared<quad_stencil>(i, ordered_curr_cell_(colon(), i), ordered_adjc_cell, adjc_face, l2s);
  // }

  // cout << "# neighboring elements: " << trans(nei_num_count) << endl;
  ASSERT(max_radius == 1 || max_radius == 2);
  if ( max_radius == 1 )
    return extract_quad4_stencil(elem);
  if ( max_radius == 2 )
    return extract_quad9_stencil(elem);
  return __LINE__;
}

int stencil_extractor::extract_quad4_stencil(vector<shared_ptr<quad_stencil>> &elem) const {
  elem.resize(cell_.size(2));
  for (size_t i = 0; i < elem.size(); ++i) {
    mati_t r_cell = ordered_curr_cell_(colon(), i);

    vector<size_t> faces{i};
    elem[i] = make_shared<quad_stencil>(r_cell, faces);
  }
  
  return 0;
}

int stencil_extractor::extract_quad9_stencil(vector<shared_ptr<quad_stencil>> &elem) const {
  mati_t quad9;
  get_quad9_elem(cell_, quad9);

  /// reorder quad9 vert permutation
  /// to satisfy basis point order
  for (size_t i = 0; i < quad9.size(2); ++i) {
    mati_t idx_perm(9, 1); {
      idx_perm[0] = 0;
      idx_perm[1] = 7;
      idx_perm[2] = 3;
      idx_perm[3] = 4;
      idx_perm[4] = 8;
      idx_perm[5] = 6;
      idx_perm[6] = 1;
      idx_perm[7] = 5;
      idx_perm[8] = 2;
    }
    quad9(colon(), i) = zjucad::matrix::temp(quad9(colon(), i)(idx_perm));
  }

  elem.resize(quad9.size(2));
  for (size_t i = 0; i < quad9.size(2); ++i) {
    const mati_t curr_elem = quad9(colon(), i);

    vector<size_t> faces{4*i+0, 4*i+1, 4*i+2, 4*i+3};    
    elem[i] = make_shared<quad_stencil>(curr_elem, faces);
  }
  
  return 0;
}

int stencil_extractor::extract_radius_2_quad9(vector<shared_ptr<quad_stencil>> &elem) const {
  // assert(cell_.size(2) == ordered_curr_cell_.size(2));
  // elem.resize(cell_.size(2));

  // mati_t quad9;
  // get_quad9_elem(cell_, quad9);

  // for (size_t i = 0; i < quad9.size(2); ++i) {
  //   mati_t idx_perm(9, 1); {
  //     idx_perm[0] = 0;
  //     idx_perm[1] = 7;
  //     idx_perm[2] = 3;
  //     idx_perm[3] = 4;
  //     idx_perm[4] = 8;
  //     idx_perm[5] = 6;
  //     idx_perm[6] = 1;
  //     idx_perm[7] = 5;
  //     idx_perm[8] = 2;
  //   }
  //   quad9(colon(), i) = zjucad::matrix::temp(quad9(colon(), i)(idx_perm));
  // }
  
  // for (size_t i = 0; i < ordered_curr_cell_.size(2); ++i) {
  //   const mati_t curr_elem = ordered_curr_cell_(colon(), i);
  //   const mati_t adjc_elem = quad9(colon(), i/4);

  //   vector<size_t> adjc_face{i/4*4+0, i/4*4+1, i/4*4+2, i/4*4+3};
    
  //   if ( i%4 == 0 ) {
  //     elem[i] = make_shared<quad_stencil>(i, curr_elem, adjc_elem, adjc_face, LS_2D_Q4_E0);
  //     continue;
  //   }
  //   if ( i%4 == 1 ) {
  //     elem[i] = make_shared<quad_stencil>(i, curr_elem, adjc_elem, adjc_face, LS_2D_Q4_E1);
  //     continue;
  //   }
  //   if ( i%4 == 2 ) {
  //     elem[i] = make_shared<quad_stencil>(i, curr_elem, adjc_elem, adjc_face, LS_2D_Q4_E2);
  //     continue;
  //   }
  //   if ( i%4 == 3 ) {
  //     elem[i] = make_shared<quad_stencil>(i, curr_elem, adjc_elem, adjc_face, LS_2D_Q4_E3);
  //     continue;
  //   }
  // }  
  return 0;
}

///===== hex stencil extractor =====///
hex_stencil_extractor::hex_stencil_extractor(const mati_t &cell)
    : cell_(cell) {}

int hex_stencil_extractor::extract(vector<shared_ptr<hexs_stencil>> &elem, const int max_radius) const {
  ASSERT(max_radius == 1 || max_radius == 2);
  if ( max_radius == 1 )
    return extract_hex8_stencil(elem);
  if ( max_radius == 2 )
    return extract_hex27_stencil(elem);
  return __LINE__;
}

int hex_stencil_extractor::extract_hex8_stencil(vector<shared_ptr<hexs_stencil>> &elem) const {
  elem.resize(cell_.size(2));
  for (size_t i = 0; i < elem.size(); ++i) {
    mati_t r_cell(8, 1); {
      r_cell[0] = cell_(0, i);
      r_cell[1] = cell_(4, i);
      r_cell[2] = cell_(2, i);
      r_cell[3] = cell_(6, i);
      r_cell[4] = cell_(1, i);
      r_cell[5] = cell_(5, i);
      r_cell[6] = cell_(3, i);
      r_cell[7] = cell_(7, i);
    }
    vector<size_t> faces{i};
    elem[i] = make_shared<hexs_stencil>(r_cell, faces);
  }
  return 0;
}

int hex_stencil_extractor::extract_hex27_stencil(vector<shared_ptr<hexs_stencil>> &elem) const {
  mati_t hex27;
  get_hex27_elem(cell_, hex27);

  elem.resize(hex27.size(2));
  for (size_t i = 0; i < elem.size(); ++i) {
    const mati_t face_idx = colon(8*i, 8*i+7);
    vector<size_t> faces(face_idx.begin(), face_idx.end());

    const mati_t r_cell = hex27(colon(), i);
    elem[i] = make_shared<hexs_stencil>(r_cell, faces);
  }
  return 0;
}

}
