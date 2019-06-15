#include "geom_util.h"

#include <ANN/ANN.h>
#include <jtflib/mesh/mesh.h>
#include <jtflib/mesh/util.h>
#include <Eigen/Geometry>
#include <zjucad/matrix/io.h>
#include <hjlib/math/blas_lapack.h>
#include <zjucad/matrix/lapack.h>

#include "util.h"

using namespace std;
using namespace jtf::mesh;
using namespace zjucad::matrix;
using namespace Eigen;

namespace bigbang {

void get_edge_elem(const mati_t &tris, mati_t &edge) {
  edge2cell_adjacent *e2c = edge2cell_adjacent::create(tris, false);
  edge.resize(2, e2c->edges_.size());
#pragma omp parallel
  for (size_t i = 0; i < edge.size(2); ++i) {
    edge(0, i) = e2c->edges_[i].first;
    edge(1, i) = e2c->edges_[i].second;
  }
  delete e2c;
}

void get_diam_elem(const mati_t &tris, mati_t &diam) {
  edge2cell_adjacent *ea = edge2cell_adjacent::create(tris, false);
  mati_t bd_ed_id;
  get_boundary_edge_idx(*ea, bd_ed_id);
  diam.resize(4, ea->edges_.size()-bd_ed_id.size());
  for(size_t ei = 0, di = 0; ei < ea->edges_.size(); ++ei) {
    pair<size_t, size_t> nb_tr_id = ea->edge2cell_[ei];
    if( ea->is_boundary_edge(nb_tr_id) ) continue;
    diam(colon(1, 2), di) = ea->get_edge(ei);
    // orient
    bool need_swap = true;
    for(size_t k = 0; k < 3; ++k) {
      if( diam(1, di) == tris(k, nb_tr_id.first) ) {
        if( diam(2, di) != tris((k+1)%3, nb_tr_id.first) )
          need_swap = false;
      }
    }
    if( need_swap )
      swap(diam(1, di), diam(2, di));
    diam(0, di) = zjucad::matrix::sum(tris(colon(), nb_tr_id.first))
        - zjucad::matrix::sum(diam(colon(1, 2), di));
    diam(3, di) = zjucad::matrix::sum(tris(colon(), nb_tr_id.second))
        - zjucad::matrix::sum(diam(colon(1, 2), di));
    ++di;
  }
  delete ea;
}

// Here we assume that the qaud mesh has subdivisional
// structure, and i-th quad9 element is composed of small
// quad elements indexed from 4*i+0~4*i+3.
void get_quad9_elem(const mati_t &quad, mati_t &quad9) {
  ASSERT(quad.size(2)%4 == 0);
  const size_t quad9_num = quad.size(2)/4;
  quad9.resize(9, quad9_num);

  for (size_t i = 0; i < quad9_num; ++i) {
    quad9(0, i) = quad(0, 4*i+0);
    quad9(1, i) = quad(1, 4*i+1);
    quad9(2, i) = quad(2, 4*i+2);
    quad9(3, i) = quad(3, 4*i+3);
    quad9(4, i) = quad(1, 4*i+0);
    quad9(5, i) = quad(2, 4*i+1);
    quad9(6, i) = quad(3, 4*i+2);
    quad9(7, i) = quad(0, 4*i+3);
    quad9(8, i) = quad(2, 4*i+0);
  }
}

/// Here we assume that the hexs mesh has subdivisional
/// structure, and i-th hex27 element is composed of
/// small cubes indexed from 8*i+0 to 8*i+7
void get_hex27_elem(const mati_t &hexs, mati_t &hex27) {
  ASSERT(hexs.size(2)%8 == 0);
  const size_t hex27_num = hexs.size(2)/8;
  hex27.resize(27, hex27_num);

  for (size_t i = 0; i < hex27_num; ++i) {
    const size_t off = 8*i;
    hex27(0, i) = hexs(0, off+0);
    hex27(1, i) = hexs(4, off+0);
    hex27(2, i) = hexs(4, off+4);

    hex27(3, i) = hexs(2, off+0);
    hex27(4, i) = hexs(6, off+0);
    hex27(5, i) = hexs(6, off+4);

    hex27(6, i) = hexs(2, off+2);
    hex27(7, i) = hexs(6, off+2);
    hex27(8, i) = hexs(6, off+6);

    hex27(9, i) = hexs(1, off+0);
    hex27(10, i) = hexs(5, off+0);
    hex27(11, i) = hexs(5, off+4);

    hex27(12, i) = hexs(3, off+0);
    hex27(13, i) = hexs(7, off+0);
    hex27(14, i) = hexs(7, off+4);
    
    hex27(15, i) = hexs(3, off+2);
    hex27(16, i) = hexs(7, off+2);
    hex27(17, i) = hexs(7, off+6);

    hex27(18, i) = hexs(1, off+1);
    hex27(19, i) = hexs(5, off+1);
    hex27(20, i) = hexs(5, off+5);

    hex27(21, i) = hexs(3, off+1);
    hex27(22, i) = hexs(7, off+1);
    hex27(23, i) = hexs(7, off+5);

    hex27(24, i) = hexs(3, off+3);
    hex27(25, i) = hexs(7, off+3);
    hex27(26, i) = hexs(7, off+7);
  }
}

double calc_tri_area(const matd_t &vert) {
  const double a = norm(vert(colon(), 1)-vert(colon(), 0));
  const double b = norm(vert(colon(), 2)-vert(colon(), 1));
  const double c = norm(vert(colon(), 0)-vert(colon(), 2));
  const double s = (a+b+c)/2;
  return sqrt(s*(s-a)*(s-b)*(s-c));
}

void remove_extra_verts(mati_t &cell, matd_t &node, mati_t *new2orig_mapping) {
  set<size_t> used_node_idx(cell.begin(), cell.end());
  if(used_node_idx.size() == node.size(2)) {
    if ( new2orig_mapping != 0 )
      *new2orig_mapping = colon(0, node.size(2)-1);
    return;
  }
  matrixst used_node_mat(used_node_idx.size(),1);
  std::copy(used_node_idx.begin(), used_node_idx.end(), used_node_mat.begin());

  map<size_t,size_t> p2p;

  matrixd new_node(node.size(1), used_node_mat.size());
  for(size_t pi = 0; pi < used_node_mat.size(); ++pi){
    new_node(colon(),pi) = node(colon(), used_node_mat[pi]);
    p2p[used_node_mat[pi]] = pi;
  }
  for(size_t pi = 0; pi < cell.size(); ++pi)
    cell[pi] = p2p[cell[pi]];

  if(new2orig_mapping != 0){
    *new2orig_mapping = used_node_mat;
  }

  node = new_node;
}

void extract_tet_surface(const mati_t &tets, const matd_t &nods, mati_t &tris, matd_t &snods) {
  using jtf::mesh::face2tet_adjacent;
  shared_ptr<face2tet_adjacent> f2t(face2tet_adjacent::create(tets));
  jtf::mesh::get_outside_face(*f2t, tris, true);
  snods = nods;
  remove_extra_verts(tris, snods);
}

void shrink_surface(const mati_t &tris, matd_t &nods, const double dist) {
  matd_t normal;
  jtf::mesh::cal_point_normal(tris, nods, normal);
  nods -= dist*normal;
}

tri_mesh_subdivider::tri_mesh_subdivider(const mati_t &tris, const matd_t &nods)
    : tris_(tris), nods_(nods), dim_(nods.size(1)), one_to_many_(1) {
  R_.resize(nods_.size(), nods_.size());
  R_.setIdentity();
  P_.resize(nods_.size(), nods_.size());
  P_.setIdentity();
}

int tri_mesh_subdivider::subdivide(mati_t &ftris, matd_t &fnods, const int sbd_sign) {
  if ( sbd_sign == 0 ) {
    ftris = tris_;
    fnods = nods_;
    return 0;
  }

  int remain = sbd_sign;

  mati_t tm_tris = tris_;
  matd_t tm_nods = nods_;

  using jtf::mesh::edge2cell_adjacent;
  
  while ( remain-- ) {
    one_to_many_ *= 4;
    
    shared_ptr<edge2cell_adjacent> e2c(edge2cell_adjacent::create(tm_tris, false));
    ftris.resize(3, 4*tm_tris.size(2));
    fnods.resize(dim_, tm_nods.size(2)+e2c->edges_.size());

    vector<int64_t> vis(e2c->edges_.size());
    std::fill(vis.begin(), vis.end(), -1);
    fnods(colon(), colon(0, tm_nods.size(2)-1)) = tm_nods;

    vector<Triplet<double>> p_trips, r_trips;
    for (size_t i = 0; i < tm_nods.size(2); ++i) {
      runtime_dim_add_diag_block<double>(dim_, i, i, 1.0, &r_trips);
      runtime_dim_add_diag_block<double>(dim_, i, i, 1.0, &p_trips);
    }

    int64_t ptris = 0, pnods = tm_nods.size(2);
    for (size_t i = 0; i < tm_tris.size(2); ++i) {
      vector<size_t> mid(3);
      for (size_t j = 0; j < 3; ++j) {
        const size_t P = tm_tris(j, i);
        const size_t Q = tm_tris((j+1)%3, i);
        const size_t eid = e2c->get_edge_idx(P, Q);
        if ( vis[eid] == -1 ) {
          mid[j] = vis[eid] = pnods;
          fnods(colon(), pnods) = (tm_nods(colon(), P)+tm_nods(colon(), Q))/2;
          runtime_dim_add_diag_block<double>(dim_, pnods, P, 0.5, &p_trips);
          runtime_dim_add_diag_block<double>(dim_, pnods, Q, 0.5, &p_trips);
          ++pnods;
        } else {
          mid[j] = vis[eid];
        }
      }
      const size_t small_tris0[] = {mid[1], mid[2], mid[0]};
      ftris(colon(), ptris++) = itr_matrix<const size_t *>(3, 1, small_tris0);
      const size_t small_tris1[] = {tm_tris(0, i), mid[0], mid[2]};
      ftris(colon(), ptris++) = itr_matrix<const size_t *>(3, 1, small_tris1);
      const size_t small_tris2[] = {tm_tris(1, i), mid[1], mid[0]};
      ftris(colon(), ptris++) = itr_matrix<const size_t *>(3, 1, small_tris2);
      const size_t small_tris3[] = {tm_tris(2, i), mid[2], mid[1]};
      ftris(colon(), ptris++) = itr_matrix<const size_t *>(3, 1, small_tris3);
    }
    
    SparseMatrix<double> tm_R(dim_*tm_nods.size(2), dim_*fnods.size(2));
    tm_R.setFromTriplets(r_trips.begin(), r_trips.end());
    R_ = R_.eval()*tm_R;

    SparseMatrix<double> tm_P(dim_*fnods.size(2), dim_*tm_nods.size(2));
    tm_P.setFromTriplets(p_trips.begin(), p_trips.end());
    P_ = tm_P*P_.eval();

    tm_tris = ftris;
    tm_nods = fnods;
  }
  
  return 0;
}

void subdivide_surface(const mati_t &tris, const matd_t &nods, mati_t &ftris, matd_t &fnods) {
  using jtf::mesh::edge2cell_adjacent;
  unique_ptr<edge2cell_adjacent> e2c(edge2cell_adjacent::create(tris, false));
  ftris.resize(3, 4 * tris.size(2));
  fnods.resize(3, nods.size(2) + e2c->edges_.size());

  vector<int> vis(e2c->edges_.size());
  std::fill(vis.begin(), vis.end(), -1);
  fnods(colon(), colon(0, nods.size(2) - 1)) = nods;

  int ptris = 0, pnods = nods.size(2);
  for (size_t i = 0; i < tris.size(2); ++i) {
    vector<size_t> mid(3);
    for (size_t j = 0; j < 3; ++j) {
      size_t P = tris(j, i);
      size_t Q = tris((j + 1) % 3, i);
      size_t eid = e2c->get_edge_idx(P, Q);
      if ( vis[eid] == -1 ) {
        mid[j] = vis[eid] = pnods;
        fnods(colon(), pnods) = 0.5 * (nods(colon(), P) + nods(colon(), Q));
        ++pnods;
      } else {
        mid[j] = vis[eid];
      }
    }
    /// @ guarantee the order of surface wouldn't change during subdivision
    const size_t small_tris0[] = {mid[1], mid[2], mid[0]};
    ftris(colon(), ptris++) = itr_matrix<const size_t *>(3, 1, small_tris0);
    const size_t small_tris1[] = {tris(0, i), mid[0], mid[2]};
    ftris(colon(), ptris++) = itr_matrix<const size_t *>(3, 1, small_tris1);
    const size_t small_tris2[] = {tris(1, i), mid[1], mid[0]};
    ftris(colon(), ptris++) = itr_matrix<const size_t *>(3, 1, small_tris2);
    const size_t small_tris3[] = {tris(2, i), mid[2], mid[1]};
    ftris(colon(), ptris++) = itr_matrix<const size_t *>(3, 1, small_tris3);
  }
}

int interp_pts_in_tets(const matd_t &v, const mati_t &tet, const matd_t &pts, csc_t &coef)
{
  const size_t tn = tet.size(2), pn = pts.size(2);

  vector<double*> pv(tn);
  matrix<double> tet_center(3, tn); {
    for(int i = 0; i < tn; ++i) {
      tet_center(colon(), i) = v(colon(), tet(colon(), i))*ones<double>(4, 1)/4;
      pv[i] = &tet_center(0, i);
    }
  }

  auto_ptr<ANNkd_tree> kdt(new ANNkd_tree(&pv[0], tn, v.size(1), 32));
  matrix<matrix<double> > bary_op(tn); {
    for(int i = 0; i < tn; ++i) {
      matrix<double> v44 = ones<double>(4, 4);
      for(int j = 0; j < 4; ++j)
        v44(colon(0, 2), j) = v(colon(), tet(j, i));
      inv(v44);
      bary_op[i] = v44;
    }
    cout << "create bary-coor operators success." << endl;
  }

  vector<Triplet<double>> trips;

  matrix<double> pt, w;
  const int ave_k = 40, iter_n = 4;
  const int max_k = static_cast<int>(40*floor(pow(2.0, iter_n)+0.5));
  matrix<double> dist(max_k);
  matrix<int> idx(max_k);
  double min_good = 1;
  int outside_cnt = 0;

  for(int pi = 0; pi < pn; ++pi) {
    if((pi%1000) == 0)
      cerr << "process " << pi << endl;

    pt = pts(colon(), pi);
    pair<int, double> best_t(-1, -10);

    for(int ki = 0, k = ave_k; ki < iter_n && k < max_k; ++ki, k*=2) {
      if(k > max_k)
        k = max_k;
      const double r2 = 1e1;
      kdt->annkSearch(&pt[0], max_k, &idx[0], &dist[0], 1e-10);
      for(int ti = (k > 40)?k/2:0; ti < k; ++ti) {
        int t_idx = idx[ti];
        w = bary_op[t_idx](colon(0, 3), colon(0, 2))*pt + bary_op[t_idx](colon(), 3);
        double good = min(w);
        if(best_t.second < good) {
          best_t.second = good;
          best_t.first = t_idx;
        }
        if(best_t.second >= 0)
          break;
      }
      if(best_t.second >= 0)
        break;
    }

    if(best_t.second < 0)
      ++outside_cnt;
    if(best_t.second < min_good)
      min_good = best_t.second;
    if(best_t.first < 0) {
      cout << "Wow, very bad point!!" << endl;
      return __LINE__;
    }

    w = bary_op[best_t.first](colon(0, 3), colon(0, 2))*pt + bary_op[best_t.first](colon(), 3);

    if(fabs(sum(w)-1) > 1e-9) {

      cout << "strange weight." << trans(w);
      cout << "sum : " << sum(w) << endl;
    }
    trips.push_back(Triplet<double>(tet(0, best_t.first), pi, w[0]));
    trips.push_back(Triplet<double>(tet(1, best_t.first), pi, w[1]));
    trips.push_back(Triplet<double>(tet(2, best_t.first), pi, w[2]));
    trips.push_back(Triplet<double>(tet(3, best_t.first), pi, w[3]));
  }
  cout << "outside pt num is: " << outside_cnt << " min_good is: " << min_good << endl;
  coef.resize(v.size(2), pn);
  coef.reserve(trips.size());
  coef.setFromTriplets(trips.begin(), trips.end());
  return 0;
}

void normalize_point_cloud(matd_t &v) {
  const matd_t center = v*ones<double>(v.size(2), 1)/v.size(2);
  v = zjucad::matrix::temp(v-center*ones<double>(1, v.size(2)));
  Vector3d minc(min(v(0,colon())), min(v(1,colon())), min(v(2,colon())));
  Vector3d maxc(max(v(0,colon())), max(v(1,colon())), max(v(2,colon())));
  const double diag_len = (maxc-minc).norm();
  v /= diag_len;
}

int get_bnd_edges(const mati_t &cell, mati_t &bnds) {
  shared_ptr<jtf::mesh::edge2cell_adjacent> e2c(jtf::mesh::edge2cell_adjacent::create(cell));

  vector<size_t> bnd_buffer;  
  for (size_t i = 0; i < cell.size(2); ++i) {
    for (size_t j = 0; j < cell.size(1); ++j) {
      const size_t p = cell(j, i), q = cell((j+1)%cell.size(1), i);
      const auto adjf = e2c->query(p, q);
      if ( e2c->is_boundary_edge(adjf) ) {
        bnd_buffer.push_back(p);
        bnd_buffer.push_back(q);
      }
    }
  }
  bnds.resize(2, bnd_buffer.size()/2);
  std::copy(bnd_buffer.begin(), bnd_buffer.end(), bnds.begin());  
  return 0;
}

///===== quad strip extractor =====///
straight_quad_strip_extractor::straight_quad_strip_extractor(const mati_t &quad)
    : quad_(quad) {
  e2c_.reset(jtf::mesh::edge2cell_adjacent::create(quad_));
}

int straight_quad_strip_extractor::extract(const size_t start_face, const char direction,
                                           vector<size_t> &face_strip) {
  ASSERT(direction == 'X' || direction == 'Y');
  const size_t axis = direction-'X';
  
  matrix<int> vis = zeros<int>(quad_.size(2), 1);
  const size_t enum_order[2][2][2] = {{{1, 2}, {3, 0}}, {{2, 3}, {0, 1}}};
  
  deque<size_t> face_seq;
  stack<size_t> q;
  vis[start_face] = 1;
  q.push(start_face);
  face_seq.push_back(start_face);

  while ( !q.empty() ) {
    const size_t curr_face = q.top();
    q.pop();
    
    for (size_t i = 0; i < 2; ++i) {
      const size_t edge_p = quad_(enum_order[axis][i][0], curr_face);
      const size_t edge_q = quad_(enum_order[axis][i][1], curr_face);
      
      const auto adjf = e2c_->query(edge_p, edge_q);
      if ( e2c_->is_boundary_edge(adjf) )
        continue;

      const size_t next_face = adjf.first+adjf.second-curr_face;
      if ( vis[next_face] )
        continue;

      vis[next_face] = 1;
      q.push(next_face);
      if ( i == 0 )
        face_seq.push_back(next_face);
      else
        face_seq.push_front(next_face);
    }
  }
  face_strip = vector<size_t>(face_seq.begin(), face_seq.end());

  return 0;
}

///===== cube strip extractor =====///
straight_cube_strip_extractor::straight_cube_strip_extractor(const mati_t &hexs)
    : hexs_(hexs) {
  f2h_.reset(jtf::mesh::face2hex_adjacent::create(hexs_));
}

int straight_cube_strip_extractor::extract(const size_t start_cube, const char direction,
                                           vector<size_t> &cube_strip) {
  ASSERT(direction == 'X' || direction == 'Y' || direction == 'Z');
  const size_t axis = direction-'X';

  matrix<int> vis = zeros<int>(hexs_.size(2), 1);
  const size_t enum_order[3][2][4] = {{{1, 3, 7, 5}, {0, 2, 6, 4}},  // along x
                                      {{2, 3, 7, 6}, {0, 1, 5, 4}},  // along y
                                      {{4, 5, 7, 6}, {0, 1, 3, 2}}}; // along z

  deque<size_t> cube_seq;
  stack<size_t> q;
  vis[start_cube] = 1;
  q.push(start_cube);
  cube_seq.push_back(start_cube);

  while ( !q.empty() ) {
    const size_t curr_cube = q.top();
    q.pop();

    for (size_t i = 0; i < 2; ++i) {
      const size_t face_a = hexs_(enum_order[axis][i][0], curr_cube);
      const size_t face_b = hexs_(enum_order[axis][i][1], curr_cube);
      const size_t face_c = hexs_(enum_order[axis][i][2], curr_cube);
      const size_t face_d = hexs_(enum_order[axis][i][3], curr_cube);

      const auto adjf = f2h_->query(face_a, face_b, face_c, face_d);
      if ( adjf.first == -1 || adjf.second == -1 )
        continue;

      const size_t next_cube = adjf.first+adjf.second-curr_cube;
      if ( vis[next_cube] )
        continue;

      vis[next_cube] = 1;
      q.push(next_cube);
      if ( i == 0 )
        cube_seq.push_back(next_cube);
      else
        cube_seq.push_front(next_cube);
    }    
  }
  cube_strip = vector<size_t>(cube_seq.begin(), cube_seq.end());
  
  return 0;
}

///===== bouddary nodes extracting interface =====///
int get_surf_mesh_bnd_nodes(const mati_t &cell, mati_t &bnd_nodes) {
  shared_ptr<jtf::mesh::edge2cell_adjacent> e2c(jtf::mesh::edge2cell_adjacent::create(cell, false));

  mati_t bnd_edges;
  jtf::mesh::get_boundary_edge(*e2c, bnd_edges);

  std::unordered_set<size_t> bnd_node_set(bnd_edges.begin(), bnd_edges.end());
  bnd_nodes.resize(bnd_node_set.size(), 1);
  std::copy(bnd_node_set.begin(), bnd_node_set.end(), bnd_nodes.begin());
  
  return 0;
}

int get_surf_mesh_bnd_nodes(const mati_t &cell, unordered_set<size_t> &bnd_nodes) {
  shared_ptr<jtf::mesh::edge2cell_adjacent> e2c(jtf::mesh::edge2cell_adjacent::create(cell, false));

  mati_t bnd_edges;
  jtf::mesh::get_boundary_edge(*e2c, bnd_edges);

  bnd_nodes = std::unordered_set<size_t>(bnd_edges.begin(), bnd_edges.end());  

  return 0;
}

}
