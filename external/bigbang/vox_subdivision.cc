#include "vox_subdivision.h"

#include <iostream>
#include "config.h"
#include "util.h"

using namespace std;
using namespace zjucad::matrix;

namespace bigbang {

static const int QUAD_EDGE_ORDER[4][2] = {{0, 1}, {1, 2}, {2, 3}, {3, 0}};

static const int VOX_EDGE_ORDER[12][2] = {{0, 1}, {1, 3}, {3, 2}, {2, 0},
                                          {4, 5}, {5, 7}, {7, 6}, {6, 4},
                                          {0, 4}, {1, 5}, {3, 7}, {2, 6}};

static const int VOX_FACE_ORDER[6][4] = {{0, 1, 3, 2}, {1, 0, 4, 5},
                                         {5, 4, 6, 7}, {7, 6, 2, 3},
                                         {2, 6, 4, 0}, {7, 3, 1, 5}};

void collect_one_vox_edges(const mati_t &vox,
                           vector<vector<size_t>> &edges) {
  edges.resize(12);
  for (size_t i = 0; i < 12; ++i) {
    for (size_t j = 0; j < 2; ++j)
      edges[i].push_back(vox[VOX_EDGE_ORDER[i][j]]);
    if ( edges[i][0] > edges[i][1] )
      std::swap(edges[i][0], edges[i][1]);
  }
}

void collect_one_vox_faces(const mati_t &vox,
                           vector<vector<size_t>> &faces) {
  faces.resize(6);
  for (size_t i = 0; i < 6; ++i) {
    for (size_t j = 0; j < 4; ++j)
      faces[i].push_back(vox[VOX_FACE_ORDER[i][j]]);
    std::sort(faces[i].begin(), faces[i].end());
  }
}

void collect_one_quad_edges(const mati_t &quad,
                            vector<vector<size_t>> &edges) {
  edges.resize(4);
  for (size_t i = 0; i < 4; ++i) {
    for (size_t j = 0; j < 2; ++j)
      edges[i].push_back(quad[QUAD_EDGE_ORDER[i][j]]);
    if ( edges[i][0] > edges[i][1] )
      std::swap(edges[i][0], edges[i][1]);
  }
}
//===============================================================================
quad_edges::quad_edges(const mati_t &cell) {
  size_t cnt = 0;
  for (size_t i = 0; i < cell.size(2); ++i) {
    vector<vector<size_t>> edges;
    collect_one_quad_edges(cell(colon(), i), edges);

    for (size_t j = 0; j < edges.size(); ++j) {
      if ( edge2idx_.find(edges[j]) == edge2idx_.end() )
        edge2idx_.insert(make_pair(edges[j], cnt++));
    }
  }
}

size_t quad_edges::edges_num() const {
  return edge2idx_.size();
}

size_t quad_edges::query_edge_idx(const size_t p, const size_t q) const {
  const vector<size_t> pq = {p, q};
  return this->query_edge_idx(pq);
}

size_t quad_edges::query_edge_idx(const vector<size_t> &ab) const {
  vector<size_t> key(ab);
  if ( key[0] > key[1] )
    std::swap(key[0], key[1]);
  const auto it = edge2idx_.find(key);
  return (it == edge2idx_.end() ? -1 : it->second);
}
//===============================================================================
vox_edges::vox_edges(const mati_t &cell) {
  size_t cnt = 0;
  for (size_t i = 0; i < cell.size(2); ++i) {
    vector<vector<size_t>> edges;
    collect_one_vox_edges(cell(colon(), i), edges);

    for (size_t j = 0; j < edges.size(); ++j) {
      if ( edge2idx_.find(edges[j]) == edge2idx_.end() )
        edge2idx_.insert(make_pair(edges[j], cnt++));
    }
  }
}

size_t vox_edges::edges_num() const {
  return edge2idx_.size();
}

size_t vox_edges::query_edge_idx(const size_t p, const size_t q) const {
  const vector<size_t> pq = {p, q};
  return this->query_edge_idx(pq);
}

size_t vox_edges::query_edge_idx(const vector<size_t> &ab) const {
  vector<size_t> key(ab);
  if ( key[0] > key[1] )
    std::swap(key[0], key[1]);
  const auto it = edge2idx_.find(key);
  return (it == edge2idx_.end() ? -1 : it->second);
}

//===============================================================================
vox_faces::vox_faces(const mati_t &cell) {
  size_t cnt = 0;
  for (size_t i = 0; i < cell.size(2); ++i) {
    vector<vector<size_t>> faces;
    collect_one_vox_faces(cell(colon(), i), faces);

    for (size_t j = 0; j < faces.size(); ++j) {
      if ( face2idx_.find(faces[j]) == face2idx_.end() ) {
        face2idx_.insert(make_pair(faces[j], cnt++));
        face2vox_.insert(make_pair(faces[j], vector<size_t>()));
      }
      face2vox_[faces[j]].push_back(i);
    }
  }
}

size_t vox_faces::faces_num() const {
  return face2idx_.size();
}

size_t vox_faces::query_face_idx(const size_t p, const size_t q, const size_t m,  const size_t n) const {
  const vector<size_t> pqmn = {p, q, m, n};
  return this->query_face_idx(pqmn);
}

size_t vox_faces::query_face_idx(const vector<size_t> &abcd) const {
  vector<size_t> key(abcd);
  std::sort(key.begin(), key.end());
  const auto it = face2idx_.find(key);
  return (it == face2idx_.end() ? -1 : it->second);
}

//===============================================================================
int subdivide_quad(const mati_t &quad, const matd_t &nods,
                   mati_t &new_quad, matd_t &new_nods,
                   Eigen::SparseMatrix<double> *P) {
  shared_ptr<quad_edges> ptre = make_shared<quad_edges>(quad);

  const size_t new_elem_num = 4*quad.size(2);
  const size_t new_nods_num = nods.size(2)+ptre->edges_num()+quad.size(2);

  new_quad.resize(4, new_elem_num);
  new_nods.resize(2, new_nods_num);

  vector<Eigen::Triplet<double>> trips;

  // construct nods
  size_t curr_nods_num = nods.size(2);
  new_nods(colon(), colon(0, curr_nods_num-1)) = nods;
  for (size_t i = 0; i < nods.size(2); ++i) {
    add_diag_block<double, 2>(i, i, 1.0, &trips);
  }
  for (auto &entry : ptre->edge2idx_) {
    new_nods(colon(), curr_nods_num+entry.second) = 0.5*(nods(colon(), entry.first[0])+nods(colon(), entry.first[1]));
    add_diag_block<double, 2>(curr_nods_num+entry.second, entry.first[0], 0.5, &trips);
    add_diag_block<double, 2>(curr_nods_num+entry.second, entry.first[1], 0.5, &trips);
  }
  curr_nods_num += ptre->edges_num();
  for (size_t i = 0; i < quad.size(2); ++i) {
    new_nods(colon(), curr_nods_num+i) = nods(colon(), quad(colon(), i))*ones<double>(4, 1)/4.0;
    for (size_t j = 0; j < 4; ++j)
      add_diag_block<double, 2>(curr_nods_num+i, quad(j, i), 0.25, &trips);
  }
  curr_nods_num += quad.size(2);
  ASSERT(curr_nods_num == new_nods_num);

  if ( P != nullptr ) {
    P->resize(2*new_nods_num, 2*nods.size(2));
    P->setFromTriplets(trips.begin(), trips.end());
  }

  const size_t edge_vert_offset = nods.size(2);
  const size_t elem_vert_offset = edge_vert_offset+ptre->edges_num();
  
  // construct topology
  #pragma omp parallel for
  for (size_t i = 0; i < quad.size(2); ++i) {
    const mati_t x = quad(colon(), i);
    
    new_quad(0, 4*i+0) = x[0];
    new_quad(1, 4*i+0) = ptre->query_edge_idx(x[0], x[1])+edge_vert_offset;
    new_quad(2, 4*i+0) = i+elem_vert_offset;
    new_quad(3, 4*i+0) = ptre->query_edge_idx(x[0], x[3])+edge_vert_offset;

    new_quad(0, 4*i+1) = ptre->query_edge_idx(x[0], x[1])+edge_vert_offset;
    new_quad(1, 4*i+1) = x[1];
    new_quad(2, 4*i+1) = ptre->query_edge_idx(x[1], x[2])+edge_vert_offset;
    new_quad(3, 4*i+1) = i+elem_vert_offset;

    new_quad(0, 4*i+2) = i+elem_vert_offset;
    new_quad(1, 4*i+2) = ptre->query_edge_idx(x[1], x[2])+edge_vert_offset;
    new_quad(2, 4*i+2) = x[2];
    new_quad(3, 4*i+2) = ptre->query_edge_idx(x[2], x[3])+edge_vert_offset;

    new_quad(0, 4*i+3) = ptre->query_edge_idx(x[0], x[3])+edge_vert_offset;
    new_quad(1, 4*i+3) = i+elem_vert_offset;
    new_quad(2, 4*i+3) = ptre->query_edge_idx(x[2], x[3])+edge_vert_offset;
    new_quad(3, 4*i+3) = x[3];
  }
  return 0;
}
//===============================================================================
int subdivide_vox(const mati_t &voxs, const matd_t &nods,
                  mati_t &new_voxs, matd_t &new_nods,
                  Eigen::SparseMatrix<double> *P) {
  shared_ptr<vox_edges> ptre = make_shared<vox_edges>(voxs);
  shared_ptr<vox_faces> ptrf = make_shared<vox_faces>(voxs);

  const size_t new_elem_num = 8*voxs.size(2);
  const size_t new_nods_num = nods.size(2)+ptre->edges_num()+ptrf->faces_num()+voxs.size(2);

  new_voxs.resize(8, new_elem_num);
  new_nods.resize(3, new_nods_num);

  vector<Eigen::Triplet<double>> trips;
  
  // construct nods
  size_t curr_nods_num = nods.size(2);
  new_nods(colon(), colon(0, curr_nods_num-1)) = nods;
  for (size_t i = 0; i < nods.size(2); ++i) {
    add_diag_block<double, 3>(i, i, 1.0, &trips);
  }
  for (auto &entry : ptre->edge2idx_) {
    new_nods(colon(), curr_nods_num+entry.second)
        = 0.5*(nods(colon(), entry.first[0])+nods(colon(), entry.first[1]));
    add_diag_block<double, 3>(curr_nods_num+entry.second, entry.first[0], 0.5, &trips);
    add_diag_block<double, 3>(curr_nods_num+entry.second, entry.first[1], 0.5, &trips);
  }
  curr_nods_num += ptre->edges_num();
  for (auto &entry : ptrf->face2idx_) {
    new_nods(colon(), curr_nods_num+entry.second)
        = 0.25*(nods(colon(), entry.first[0])+nods(colon(), entry.first[1])+nods(colon(), entry.first[2])+nods(colon(), entry.first[3]));
    add_diag_block<double, 3>(curr_nods_num+entry.second, entry.first[0], 0.25, &trips);
    add_diag_block<double, 3>(curr_nods_num+entry.second, entry.first[1], 0.25, &trips);
    add_diag_block<double, 3>(curr_nods_num+entry.second, entry.first[2], 0.25, &trips);
    add_diag_block<double, 3>(curr_nods_num+entry.second, entry.first[3], 0.25, &trips);
  }
  curr_nods_num += ptrf->faces_num();
  for (size_t i = 0; i < voxs.size(2); ++i) {
    new_nods(colon(), curr_nods_num+i) = nods(colon(), voxs(colon(), i))*ones<double>(8, 1)/8.0;
    for (size_t j = 0; j < 8; ++j)
      add_diag_block<double, 3>(curr_nods_num+i, voxs(j, i), 0.125, &trips);
  }
  curr_nods_num += voxs.size(2);
  ASSERT(curr_nods_num == new_nods_num);

  if ( P != nullptr ) {
    P->resize(3*new_nods_num, 3*nods.size(2));
    P->setFromTriplets(trips.begin(), trips.end());
  }

  const size_t edge_vert_offset = nods.size(2);
  const size_t face_vert_offset = edge_vert_offset+ptre->edges_num();
  const size_t elem_vert_offset = face_vert_offset+ptrf->faces_num();
  
  // construct topology
  #pragma omp parallel for
  for (size_t i = 0; i < voxs.size(2); ++i) {
    const mati_t x = voxs(colon(), i);
    
    new_voxs(0, 8*i+0) = x[0];
    new_voxs(1, 8*i+0) = ptre->query_edge_idx(x[0], x[1])+edge_vert_offset;
    new_voxs(2, 8*i+0) = ptre->query_edge_idx(x[0], x[2])+edge_vert_offset;
    new_voxs(3, 8*i+0) = ptrf->query_face_idx(x[0], x[1], x[2], x[3])+face_vert_offset;
    new_voxs(4, 8*i+0) = ptre->query_edge_idx(x[0], x[4])+edge_vert_offset;
    new_voxs(5, 8*i+0) = ptrf->query_face_idx(x[0], x[1], x[4], x[5])+face_vert_offset;
    new_voxs(6, 8*i+0) = ptrf->query_face_idx(x[0], x[2], x[6], x[4])+face_vert_offset;
    new_voxs(7, 8*i+0) = i+elem_vert_offset;

    new_voxs(0, 8*i+1) = ptre->query_edge_idx(x[0], x[1])+edge_vert_offset;
    new_voxs(1, 8*i+1) = x[1];
    new_voxs(2, 8*i+1) = ptrf->query_face_idx(x[0], x[1], x[2], x[3])+face_vert_offset;
    new_voxs(3, 8*i+1) = ptre->query_edge_idx(x[1], x[3])+edge_vert_offset;
    new_voxs(4, 8*i+1) = ptrf->query_face_idx(x[0], x[1], x[4], x[5])+face_vert_offset;
    new_voxs(5, 8*i+1) = ptre->query_edge_idx(x[1], x[5])+edge_vert_offset;
    new_voxs(6, 8*i+1) = i+elem_vert_offset;
    new_voxs(7, 8*i+1) = ptrf->query_face_idx(x[1], x[3], x[7], x[5])+face_vert_offset;

    new_voxs(0, 8*i+2) = ptre->query_edge_idx(x[0], x[2])+edge_vert_offset;
    new_voxs(1, 8*i+2) = ptrf->query_face_idx(x[0], x[1], x[2], x[3])+face_vert_offset;
    new_voxs(2, 8*i+2) = x[2];
    new_voxs(3, 8*i+2) = ptre->query_edge_idx(x[2], x[3])+edge_vert_offset;
    new_voxs(4, 8*i+2) = ptrf->query_face_idx(x[0], x[2], x[6], x[4])+face_vert_offset;
    new_voxs(5, 8*i+2) = i+elem_vert_offset;
    new_voxs(6, 8*i+2) = ptre->query_edge_idx(x[2], x[6])+edge_vert_offset;
    new_voxs(7, 8*i+2) = ptrf->query_face_idx(x[2], x[3], x[7], x[6])+face_vert_offset;

    new_voxs(0, 8*i+3) = ptrf->query_face_idx(x[0], x[1], x[2], x[3])+face_vert_offset;
    new_voxs(1, 8*i+3) = ptre->query_edge_idx(x[1], x[3])+edge_vert_offset;
    new_voxs(2, 8*i+3) = ptre->query_edge_idx(x[2], x[3])+edge_vert_offset;
    new_voxs(3, 8*i+3) = x[3];
    new_voxs(4, 8*i+3) = i+elem_vert_offset;
    new_voxs(5, 8*i+3) = ptrf->query_face_idx(x[1], x[3], x[7], x[5])+face_vert_offset;
    new_voxs(6, 8*i+3) = ptrf->query_face_idx(x[2], x[3], x[7], x[6])+face_vert_offset;
    new_voxs(7, 8*i+3) = ptre->query_edge_idx(x[3], x[7])+edge_vert_offset;

    new_voxs(0, 8*i+4) = ptre->query_edge_idx(x[0], x[4])+edge_vert_offset;
    new_voxs(1, 8*i+4) = ptrf->query_face_idx(x[0], x[1], x[5], x[4])+face_vert_offset;
    new_voxs(2, 8*i+4) = ptrf->query_face_idx(x[0], x[2], x[6], x[4])+face_vert_offset;
    new_voxs(3, 8*i+4) = i+elem_vert_offset;
    new_voxs(4, 8*i+4) = x[4];
    new_voxs(5, 8*i+4) = ptre->query_edge_idx(x[4], x[5])+edge_vert_offset;
    new_voxs(6, 8*i+4) = ptre->query_edge_idx(x[4], x[6])+edge_vert_offset;
    new_voxs(7, 8*i+4) = ptrf->query_face_idx(x[4], x[5], x[6], x[7])+face_vert_offset;

    new_voxs(0, 8*i+5) = ptrf->query_face_idx(x[0], x[1], x[4], x[5])+face_vert_offset;
    new_voxs(1, 8*i+5) = ptre->query_edge_idx(x[1], x[5])+edge_vert_offset;
    new_voxs(2, 8*i+5) = i+elem_vert_offset;
    new_voxs(3, 8*i+5) = ptrf->query_face_idx(x[1], x[3], x[7], x[5])+face_vert_offset;
    new_voxs(4, 8*i+5) = ptre->query_edge_idx(x[4], x[5])+edge_vert_offset;
    new_voxs(5, 8*i+5) = x[5];
    new_voxs(6, 8*i+5) = ptrf->query_face_idx(x[4], x[5], x[6], x[7])+face_vert_offset;
    new_voxs(7, 8*i+5) = ptre->query_edge_idx(x[5], x[7])+edge_vert_offset;

    new_voxs(0, 8*i+6) = ptrf->query_face_idx(x[0], x[2], x[6], x[4])+face_vert_offset;
    new_voxs(1, 8*i+6) = i+elem_vert_offset;
    new_voxs(2, 8*i+6) = ptre->query_edge_idx(x[2], x[6])+edge_vert_offset;
    new_voxs(3, 8*i+6) = ptrf->query_face_idx(x[2], x[3], x[7], x[6])+face_vert_offset;
    new_voxs(4, 8*i+6) = ptre->query_edge_idx(x[4], x[6])+edge_vert_offset;
    new_voxs(5, 8*i+6) = ptrf->query_face_idx(x[4], x[5], x[6], x[7])+face_vert_offset;
    new_voxs(6, 8*i+6) = x[6];
    new_voxs(7, 8*i+6) = ptre->query_edge_idx(x[6], x[7])+edge_vert_offset;

    new_voxs(0, 8*i+7) = i+elem_vert_offset;
    new_voxs(1, 8*i+7) = ptrf->query_face_idx(x[1], x[3], x[7], x[5])+face_vert_offset;
    new_voxs(2, 8*i+7) = ptrf->query_face_idx(x[2], x[3], x[7], x[6])+face_vert_offset;
    new_voxs(3, 8*i+7) = ptre->query_edge_idx(x[3], x[7])+edge_vert_offset;
    new_voxs(4, 8*i+7) = ptrf->query_face_idx(x[4], x[5], x[6], x[7])+face_vert_offset;
    new_voxs(5, 8*i+7) = ptre->query_edge_idx(x[5], x[7])+edge_vert_offset;
    new_voxs(6, 8*i+7) = ptre->query_edge_idx(x[6], x[7])+edge_vert_offset;
    new_voxs(7, 8*i+7) = x[7];
  }
  return 0;
}

}
