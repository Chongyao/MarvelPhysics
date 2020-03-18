#include "Sparsification.h"
#include <iostream>
#include <list>
using namespace Eigen;
using namespace std;
namespace marvel{

//===================graph========================//
//construct from a dense Matrix
void Adjc_graph::init(){
  vertices_.resize(dof_, nullptr);
  labels_.resize(dof_, mark_state::unmarked);
  VectorXi all_v = VectorXi::LinSpaced(1, 0, dof_ - 1);
  unmarked_vertices_.insert(all_v.data(), all_v.data() + all_v.size());
  
  #pragma omp parallel for
  for(size_t i = 0; i < dof_; ++i){
    vertices_[i] = make_shared<unordered_set<size_t>>();
  }
  return;
}

Adjc_graph::Adjc_graph(const Eigen::MatrixXd& L)
    :dof_(L.rows()){
  init();
  dig_vals_ = L * VectorXd::Ones(dof_);

  size_t num_edges = 0;
  for(size_t i = 0; i < dof_; ++i){
    for(size_t j = i + 1; j < dof_; ++j){
      if(fabs(L(i, j)) < 1e-8)
        continue;
      vertices_[i]->insert(num_edges);
      vertices_[j]->insert(num_edges);
      edges_.push_back(make_shared<TPL>(i, j, -L(i, j)));
      ++num_edges;
    }
  }
}


Adjc_graph::Adjc_graph(const SparseMatrix<double>& L):dof_(L.rows()){
  init();

  dig_vals_ = L * VectorXd::Ones(dof_);
  size_t num_edges = 0;
  for(int k=0; k<L.outerSize(); ++k)
    for (SparseMatrix<double>::InnerIterator it(L,k); it; ++it){
      if(it.index() >= k)
        break;
      vertices_[it.row()]->insert(num_edges);
      vertices_[it.col()]->insert(num_edges);
      edges_.push_back(make_shared<TPL>(it.row(), it.col(), -it.value()));
      ++num_edges;
    }
}


int Adjc_graph::build_mat_from_graph(vector<TPL>& trips)const{
  #pragma omp parallel for
  for(size_t i = 0; i < edges_.size(); ++i){
    if(edges_[i] == nullptr)
      continue;
    
    vector<Triplet<double>> trips_e;
    const auto& trip = *edges_[i];
    trips_e.push_back(TPL(trip.row(), trip.col(), -trip.value()));
    trips_e.push_back(TPL(trip.col(), trip.row(), -trip.value()));
    trips_e.push_back(TPL(trip.row(), trip.row(), trip.value()));
    trips_e.push_back(TPL(trip.col(), trip.col(), trip.value()));

    #pragma omp critical
    {
      trips.insert(trips.end(), trips_e.begin(), trips_e.end());
    }
  }

  for(size_t i = 0; i < dof_; ++i){
    trips.push_back(TPL(i, i, dig_vals_[i]));
  }
  return 0;
}


void Adjc_graph::sparsify_one_edge(const size_t edge_id){
  const auto& trip = *edges_[edge_id];
  vertices_[trip.row()]->erase(edge_id);
  vertices_[trip.col()]->erase(edge_id);
  
  labels_[trip.row()] = mark_state::fine;
  labels_[trip.col()] = mark_state::fine;
  unmarked_vertices_.erase(trip.row());
  unmarked_vertices_.erase(trip.col());
  
  edges_[edge_id] = nullptr;
  return;
}

void Adjc_graph::compensate_one_edge(const size_t edge_id_ik, const size_t edge_id_jk, const double& w){
  const auto& trip_ik = edges_[edge_id_ik];
  edges_[edge_id_ik] = make_shared<TPL>(trip_ik->row(), trip_ik->col(), trip_ik->value() + w);
  const auto& trip_jk = edges_[edge_id_jk];
  edges_[edge_id_jk] = make_shared<TPL>(trip_jk->row(), trip_jk->col(), trip_jk->value() + w);
  return;
}

bool Adjc_graph::is_connect(const size_t i, const size_t j, size_t& edge_ij)const{
  const pair<size_t, size_t> host_target = vertices_[i]->size() < vertices_[j]->size() ? pair<size_t, size_t>({i, j}) : pair<size_t, size_t>({j, i});
  const size_t target = host_target.second;
  const size_t host = host_target.first;
  for(const auto& edge_id: *(vertices_[host])){
    const auto& trip = *edges_[edge_id];
    if(trip.row() == target || trip.col() == target){
      edge_ij = edge_id;
      return true;
    }

  }
  return false;
}

int Adjc_graph::sparsify_one_tri(const size_t edge_id_i, const size_t edge_id_j, const size_t edge_id_k, size_t& sparsified_edge_id){
  vector<pair<size_t, double>> ws;{
    ws.push_back({edge_id_i, fabs(edges_[edge_id_i]->value())});
    ws.push_back({edge_id_j, fabs(edges_[edge_id_j]->value())});
    ws.push_back({edge_id_k, fabs(edges_[edge_id_k]->value())});
  }
  sort(ws.begin(), ws.end(), [](const pair<size_t, double>&lhs, const pair<size_t,double>&rhs)->bool{
      return lhs.second < rhs.second;
    });
  sparsify_one_edge(ws[0].first);
  compensate_one_edge(ws[1].first, ws[2].first, ws[0].second);
  sparsified_edge_id = ws[0].first;
  return 0;
}

int Adjc_graph::Sparsification(){
  for(size_t i = 0; i < dof_; ++i){
    //loop for triangles of i
    const size_t num_adjc = vertices_[i]->size();
    // unordered_map<size_t, size_t> adjc_V_E;
    list<pair<size_t, size_t>>adjc_V_E;
    for(const auto& edge_id: *(vertices_[i])){
      adjc_V_E.push_back({edges_[edge_id]->row() == i ? edges_[edge_id]->col() : edges_[edge_id]->row(), edge_id});
    }
      

    for(auto it_a = adjc_V_E.begin();it_a!= adjc_V_E.end();++it_a)
      for(auto it_b = next(it_a) ; it_b != adjc_V_E.end(); ++ it_b){
        size_t edge_id = -1;
        if(is_connect(it_a->first, it_b->first, edge_id)){
          size_t sp_edge_id = -1;
          sparsify_one_tri(it_a->second, it_b->second, edge_id, sp_edge_id);
          if(sp_edge_id == it_a->second){
            adjc_V_E.erase(it_a);
            --it_a;
            break;
          }else if(sp_edge_id == it_b->second){
            adjc_V_E.erase(it_b);
            --it_b;
            continue;
          }
        }
      }

    //set unmarked neighbours as coarse
    for(const auto& it : adjc_V_E){
      if(labels_[it.first] == mark_state::unmarked){
        labels_[it.first] = mark_state::coarse;
        unmarked_vertices_.erase(it.first);
      }

    }
  }

  //set unmarked vertices that have fine neighbors to coarse, else fine
  for(const auto& v : unmarked_vertices_){
    bool have_fine_neighbors = false;
    for(const auto& edge_id : *(vertices_[v])){
      const size_t neighbor = edges_[edge_id]->row() == v ? edges_[edge_id]->col() : edges_[edge_id]->row();
      if(labels_[neighbor] == mark_state::fine){
        have_fine_neighbors = true;
        break;
      }
    }
    labels_[v] = have_fine_neighbors ? mark_state::coarse : mark_state::fine;
    //This may cause bug
    unmarked_vertices_.erase(v);
  }
  //for any fine-fine connection, set one of the endpoints to coarse
  for(const auto& edge_ptr : edges_)
    if(edge_ptr != nullptr
       && labels_[edge_ptr->row()] == mark_state::fine
       && labels_[edge_ptr->col()] == mark_state::fine)
      labels_[edge_ptr->row()] = mark_state::coarse;
  
  //Set any coarse variables connected only to coarse variable as fine
  for(size_t i = 0; i < dof_; ++i){
    if(labels_[i] == mark_state::fine)
      continue;
    size_t num_coarse_neighbors = 0;
    for(const auto& edge_id : *(vertices_[i])){
      const size_t neighbor = edges_[edge_id]->row() == i ? edges_[edge_id]->col() : edges_[edge_id]->row();
      if(labels_[neighbor] == mark_state::coarse)
        ++num_coarse_neighbors;
    }
    if(num_coarse_neighbors == 1)
      labels_[i] = mark_state::fine;
  }
  
  return 0;
}

int Adjc_graph::build_reordered_mat_from_graph(std::vector<TPL>& trips){
  assert(unmarked_vertices_.size() == 0);
  VectorXi perm_vec(dof_);
  VectorXi perm_vec_inverse(dof_);
  size_t coarse_id = 0, fine_id = dof_ - 1;
  for(size_t i = 0; i < dof_; ++i)
    if(labels_[i] == mark_state::coarse){
      perm_vec(coarse_id) = i;
      perm_vec_inverse(i) = coarse_id;
      ++coarse_id;
    }else{
      perm_vec(fine_id) = i;
      perm_vec_inverse(i) = fine_id;
      --fine_id;
    }

  num_coarse_ = coarse_id;


  #pragma omp parallel for
  for(size_t i = 0; i < edges_.size(); ++i){
    if(edges_[i] == nullptr)
      continue;
    
    vector<Triplet<double>> trips_e;
    const auto& trip = *edges_[i];
    const size_t
        row_id = perm_vec_inverse(trip.row()),
        col_id = perm_vec_inverse(trip.col());
    
    trips_e.push_back(TPL(row_id, col_id, -trip.value()));
    trips_e.push_back(TPL(col_id, row_id, -trip.value()));
    trips_e.push_back(TPL(row_id, row_id, trip.value()));
    trips_e.push_back(TPL(col_id, col_id, trip.value()));

    #pragma omp critical
    {
      trips.insert(trips.end(), trips_e.begin(), trips_e.end());
    }
  }

  for(size_t i = 0; i < dof_; ++i){
    trips.push_back(TPL(i, i, dig_vals_[i]));
  }
  return 0;
}
//===================graph========================//

int Schur_complement(const SparseMatrix<double>& L, const size_t& coarse_num,
                     SparseMatrix<double>& topleft,
                     SparseMatrix<double>& bottomright){
  const size_t fine_num = L.rows() - coarse_num;
  bottomright = L.bottomRightCorner(fine_num, fine_num);
  const SparseMatrix<double>
      L_cc = L.topLeftCorner(coarse_num, coarse_num),
      L_fc = L.bottomLeftCorner(fine_num, coarse_num);

  SimplicialLLT<SparseMatrix<double>> llt(bottomright);
  topleft = L_cc - L_fc.transpose() * llt.solve(L_fc);
  
  return 0;
}



}
