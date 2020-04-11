#include "hsc.h"
#include <boost/property_tree/ptree.hpp>
namespace marvel{
using namespace std;
using namespace Eigen;
using namespace boost::property_tree;

HSC::HSC(VS<layer>& layers, const VS<transfer>& transfers,const bool pre_sm, const bool post_sm, const bool diag_PD, const size_t num_mu, const size_t num_V)
    :multigrid_process({0}, layers, transfers), num_layers_(layers.size()),pre_sm_(pre_sm), post_sm_(post_sm), diag_PD_(diag_PD), num_mu_(num_mu), num_V_(num_V){}


VectorXd HSC::solve(const VectorXd& b){
  layers_[0]->rhs_ = b;
  for(size_t i = 0; i < num_V_; ++i){
    unified_multigrid(0, pre_sm_, post_sm_, num_mu_, diag_PD_);
  }
  return layers_[0]->u_;
}

void HSC::unified_multigrid(
    const size_t& curr_layer_id, const bool& pre_sm,
    const bool& post_sm, const size_t& num_mu, const bool& diag_PD){
  auto& layer_ptr =  layers_[curr_layer_id];
  auto& e = layer_ptr->u_, r = layer_ptr->rhs_;
  auto& A = layer_ptr->A_; 
  // e.setZero();
  if(pre_sm){
    relax(curr_layer_id);
  }
  restrict(curr_layer_id + 1);
  if(curr_layer_id + 1 == num_layers_ - 1)
    relax(curr_layer_id + 1);
  else{
    layers_[curr_layer_id + 1]->u_.setZero();
    for(size_t i = 0; i < num_mu; ++i)
      unified_multigrid(curr_layer_id + 1, pre_sm, post_sm, num_V, diag_PD);
  }
  correct(curr_layer_id);
  
  if(diag_PD){
    VectorXd diag_vals = A.diagonal(), res = layer_ptr->get_residual();
    e.array() += res.array() / diag_vals.array();
  }
  if(post_sm){
    relax(curr_layer_id);
  }
  return;
}


HSC set_hierarchy(const SPM& L, const bool truncate, const ptree& pt){
  //get parameters
  const size_t gs_itrs = pt.get<size_t>("gs_itrs", 2);
  const string sol_type_str = pt.get<string>("solver_type", "WJ");
  const size_t coarest_num = pt.get<size_t>("coarest_num", 1024);
  solver_type sol_type;
  if(sol_type_str == "GS")
    sol_type = solver_type::GS;
  else if(sol_type_str == "WJ")
    sol_type = solver_type::WJ;

  
  VS<Adjc_graph> graphs;
  VS<Sparsify> sp_ops;

  VS<layer> layers;
  VS<transfer> transfers;
  size_t dim = L.rows();
  shared_ptr<SPM>
      L_now = make_shared<SPM>(L),
      L_next = make_shared<SPM>(L.rows(), L.rows());
  
  do{
    layers.push_back(make_shared<layer>(*L_now, sol_type, gs_itrs));
    graphs.push_back(make_shared<Adjc_graph>(L_now, 1));
    sp_ops.push_back(make_shared<Sparsify>(L_now.rows()));
    
    auto& sp_op = sp_ops[sp_ops.size() - 1];
    sp_op.sparsify_and_compensate(graph);
    sp_op.post_coloring(graph);
    sp_op.reorder_coarse_and_fine();
    VectorXd perm_inv = sp_op.get_perm_inv();

    {//get sparsified and reordered L
      vector<Triplet<double>> trips;
      graphs[graphs.size() - 1]->build_reordered_mat_from_graph(perm_inv, trips);
      add_dig_vals(*L_now, trips, perm_inv);
      L_now.setZero();
      L_now->reserve(trips.size());
      L_now->setFromTriplets(trips.begin(), trips.end());
    }
    {//get S and coraser L
      SPM S(L_now.rows(), sp_op->num_coarse_);
      schur_complement(sp_op->num_coarse_, *L_now, *L_next, S);
      PermutationMatrix<-1, -1> perm = perm_inv.asPermutation();
      S = perm.transpose() * S;
      transfers.push_back(make_shared<transfer>(S, S.transpose()));
      dim = sp_op->num_coarse_;
    }
    L_now = L_next;
  }while(dim > coarest_num);
 
  
}


}
