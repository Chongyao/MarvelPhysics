#include "multigrid.h"
#include "gauss_seidel.h"
#include "inverse_isoparametric_map.h"
#include "eigen_ext.h"
#include <unsupported/Eigen/KroneckerProduct>
#include <Eigen/Eigenvalues>
#include <fstream>
#include <map>
namespace marvel{
using namespace std;
using namespace Eigen;
//============================transfer===========================//
transfer::transfer(const SPM& I, const SPM& R): I_(I), R_(R){}


transfer get_transfer(const MatrixXd& nods_H, const MatrixXi& cells_H, const MatrixXd& nods_h, const MatrixXi& cells_h){
  const size_t num_per_cell_ = 8;
  assert(cells_h.cols() % cells_H.cols() == 0);
  assert(cells_H.rows() == cells_h.rows() && cells_H.rows() == num_per_cell_);
  
  const size_t one_to_many = cells_h.cols() / cells_H.cols();
  
  const VectorXi all_rows = Matrix<int, 3, 1>::LinSpaced(3, 0, 2);
  
  vector<Triplet<double>> trips_I; //n_h X n_H
  
  vector<bool> if_calc (nods_h.cols(), false);
  for(size_t i = 0; i < cells_H.cols(); ++i){
    const Matrix<double, 3, num_per_cell_> x_H_cell = indexing(nods_H, all_rows, cells_H.col(i));
    vector<Triplet<double>> trips_patch;
    #pragma omp parallel for
    for(size_t j = 0; j < one_to_many; ++j){
      for(size_t k = 0; k < num_per_cell_; ++k){
        const size_t nod_h_id = cells_h(k, i * one_to_many + j);
        if(if_calc[nod_h_id])
          continue;
        else{
          if_calc[nod_h_id] = true;
          Matrix<double, num_per_cell_, 1> w;
          VectorXd xi = inverse_isoparametric_hex(nods_h.col(nod_h_id), x_H_cell.data(), w.data());
          for(size_t m = 0; m < num_per_cell_; ++m){
            trips_patch.push_back(Triplet<double>(nod_h_id, cells_H(m, i), w(m)));
          }
        }
      }
    }
    trips_I.insert(trips_I.end(), trips_patch.begin(), trips_patch.end());
  }
  SparseMatrix<double> I(nods_h.cols(), nods_H.cols());{
    I.reserve(trips_I.size());
    I.setFromTriplets(trips_I.begin(), trips_I.end());
    Matrix3d Iden = Matrix3d::Identity();
    I = kroneckerProduct(I, Iden).eval();
  }

  return transfer(I, I.transpose());
}
//============================transfer===========================//
//============================layer===========================//
layer::layer(const SPM& A, const bool if_direct, const size_t itrs):A_(A), u_(VectorXd::Zero(A.rows())), if_direct_(if_direct), itrs_(itrs){}

int layer::solve(){
  if(if_direct_){
    ConjugateGradient<SPM, Lower|Upper> cg;
    cg.compute(A_);
    u_ = cg.solve(rhs_);
  }else{
    // Gauss_seidel GS(A_, rhs_.data(), itrs_, 1e-4, 1.0);
    // VectorXd solution(rhs_.size());
    // GS.solve(u_.data(), solution.data());
    gauss_seidel_solver(A_, rhs_, u_, itrs_);
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
  size_t layer_id = 0;
  for(const auto& op : process_){
    layer_id += op;
    if(op == 1){
      relax(layer_id - 1);
      restrict(layer_id);
    }else{
      correct(layer_id);
      relax(layer_id);
    }
    if(layer_id == num_layers - 1){
      relax(layer_id);
    }

  }
  return 0;
}

int multigrid_process::relax(const size_t layer_id){
  cout << "relax layer " << layer_id << endl;
  cout << "residual befor GS " << layers_[layer_id]->get_residual().norm() << endl;
  // static size_t t = 0;
  // auto A = MatrixXd(layers_[layer_id]->A_);
  // EigenSolver<MatrixXd> s(A);
  // MatrixXd eigvecs = s.eigenvectors().real();
  // VectorXd eigvalues = s.eigenvalues().real();
  // map<double, size_t> order;{
  //   for(size_t i = 0; i < eigvalues.rows(); ++i){
  //     order.insert({eigvalues(i), i});
  //   }
  // }

  
  // {//check spectrum
  //   VectorXd r = layers_[layer_id]->get_residual();
    
  //   ofstream ofs(to_string(t)+"_before.txt");
  //   for(auto& v : order){
  //     const size_t i = v.second;
  //     VectorXd tmp = eigvecs.col(i) / eigvecs.col(i).norm();
  //     double w = r.dot(tmp);
  //     ofs << fabs(w) << endl;
  //     r -= tmp * w;
  //   }
  //   cout << "after proj " << r.norm() << endl;
  //   ofs.close();
  // }
  
  layers_[layer_id]->solve();
  cout << "residual after GS" << layers_[layer_id]->get_residual().norm() << endl;
  // {//check spectrum
  //   VectorXd r = layers_[layer_id]->get_residual();
  //   ofstream ofs(to_string(t)+"_after.txt");
  //   for(auto& v : order){
  //     const size_t i = v.second;
  //     VectorXd tmp = eigvecs.col(i) / eigvecs.col(i).norm();
  //     double w = r.dot(tmp);
  //     ofs << -fabs(w) << endl;
  //     r -= tmp * w;
  //   }
  //   cout << "after proj " << r.norm() << endl;
  //   ofs.close();
  //   ++t;
  // }
  return 0;
}

int multigrid_process::restrict(const size_t layer_id){
  assert(layer_id< layers_.size());
  assert(transfers_[layer_id - 1]->R_.cols() == layers_[layer_id - 1]->get_residual().rows());
  layers_[layer_id]->rhs_
      = transfers_[layer_id - 1]->R_ * layers_[layer_id - 1]->get_residual();
  return 0;
}

int multigrid_process::correct(const size_t layer_id){
  assert(layer_id - 1 >= 0);
  assert(transfers_[layer_id]->I_.cols() == layers_[layer_id + 1]->u_.rows());
  layers_[layer_id]->u_
      += transfers_[layer_id]->I_ * layers_[layer_id + 1]->u_;
  return 0;
}


}
