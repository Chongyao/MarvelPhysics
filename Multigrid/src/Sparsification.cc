#include "Sparsification.h"
#include <iostream>
#include <list>
using namespace Eigen;
using namespace std;
namespace marvel{

//===================graph========================//
//construct from a dense Matrix
Adjc_graph::Adjc_graph(const Eigen::MatrixXd& L)
    :dof_(L.rows()){
  vertices_.resize(dof_, nullptr);
  #pragma omp parallel for
  for(size_t i = 0; i < dof_; ++i){
    vertices_[i] = make_shared<unordered_set<size_t>>();
  }

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
  return 0;
}


void Adjc_graph::sparsify_one_edge(const size_t edge_id){
  const auto& trip = *edges_[edge_id];
  vertices_[trip.row()]->erase(edge_id);
  vertices_[trip.col()]->erase(edge_id);  
  edges_[edge_id] = nullptr;
  return;
}

void Adjc_graph::compensate_one_edge(const size_t edge_id_ik, const size_t edge_id_jk, const double& w){
  const auto& trip_ik = *edges_[edge_id_ik];
  edges_[edge_id_ik] = make_shared<TPL>(trip_ik.row(), trip_ik.col(), trip_ik.value() + w);
  const auto& trip_jk = *edges_[edge_id_jk];
  edges_[edge_id_jk] = make_shared<TPL>(trip_jk.row(), trip_jk.col(), trip_jk.value() + w);
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
    ws.push_back({edge_id_i, edges_[edge_id_i]->value()});
    ws.push_back({edge_id_j, edges_[edge_id_j]->value()});
    ws.push_back({edge_id_k, edges_[edge_id_k]->value()});
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
    for(const auto& edge_id: *(vertices_[i]))
      adjc_V_E.push_back({edges_[edge_id]->row() == i ? edges_[edge_id]->col() : i, edge_id});

    for(auto it_a = adjc_V_E.begin();it_a!= adjc_V_E.end();++it_a)
      for(auto it_b = next(it_a) ; it_b != adjc_V_E.end(); ++ it_b){
        size_t edge_id = -1;
        if(is_connect(it_a->second, it_b->second, edge_id)){
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
  }
}
//===================graph========================//



}
