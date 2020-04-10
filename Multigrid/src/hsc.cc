#include "hsc.h"
namespace marvel{
using namespace std;
using namespace Eigen;

HSC::HSC(VS<layer>& layers, const VS<transfer>& transfers)
    :multigrid_process({0}, layers, transfers), num_layers_(layers.size()){
  process_.resize((layers.size() - 1) * 2, -1);
  for(size_t i = 0; i < layers.size() - 1; ++i){
    process_[i] = 1;
  }
}

VectorXd HSC::solve(const VectorXd& b){
  
}

void HSC::unified_multigrid(
    const size_t& curr_layer_id, const bool& pre_sm,
    const bool& post_sm, const size_t& num_V, const bool& diag_PD){
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
    for(size_t i = 0; i < num_V; ++i)
      unified_multigrid(curr_layer_id + 1, pre_sm, post_sm, num_V, diag_PD);
  }
  correct(curr_layer_id);
  
  
}



}
