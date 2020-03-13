#include "multigrid.h"
#include "gauss_seidel.h"

namespace marvel{
using namespace std;
using namespace Eigen;
//============================transfer===========================//
transfer::transfer(const SPM& I, const SPM& R): I_(I), R_(R){}
//============================transfer===========================//
//============================layer===========================//
layer::layer(const SPM& A, const bool if_direct, const size_t itrs):A_(A), u_(VectorXd::Zero(A.rows())), if_direct_(if_direct), itrs_(itrs){}

int layer::solve(){
  if(if_direct_){
    ConjugateGradient<SPM, Lower|Upper> cg;
    cg.compute(A_);
    u_ = cg.solve(rhs_);
  }else{
    gauss_seidel(A_, rhs_, u_, itrs_);
  }
  return 0;
}

VectorXd layer::get_residual()const{
  return rhs_ - A_ * u_;
}

//============================layer===========================//

//============================multigrid process===========================//

multigrid_process::multigrid_process(const vector<int>& process, VS<layer>& layers, const VS<transfer>& transfers):layers_(layers), transfers_(transfers), process_(process), dof_(layers[0]->A_.rows()){
  assert(layers_.size() == transfers_.size() + 1);
  int sum = 0;
  for(const auto& op : process){
    sum += op;
  }
  assert(sum == 0);
}

int multigrid_process::execute(double* solution){
  const size_t num_layers = layers_.size();
  for(const auto& op : process_){
    if(op == 1){
      relax();
      restrict();
    }else{
      correct();
      relax();
    }
  }
  return 0;
}

int multigrid_process::relax(){
  layers_[layer_id_]->solve();
  return 0;
}

int multigrid_process::restrict(){
  assert(layer_id_ + 1 < layers_.size());
  layers_[layer_id_ + 1]->rhs_
      = transfers_[layer_id_]->R_ * layers_[layer_id_]->get_residual();
  return 0;
}

int multigrid_process::correct(){
  assert(layer_id - 1 >= 0);
  layers_[layer_id_ - 1]->u_
      += transfers_[layer_id_ - 1]->I_ * layers_[layer_id_]->u_;
  return 0;
}


}
