#include "pcg.h"
#include <iostream>
namespace marvel{
using namespace std;
using namespace Eigen;

PCG::PCG(const SPM& A, const size_t& max_itrs, const double& tol, const precond_type& M): A_(A), max_itrs_(max_itrs), tol_(tol), M_(M), dim_(A.rows()){}


Eigen::VectorXd PCG::solve(const VectorXd& b){
  assert(b.size() == dim_);

  
  VectorXd xk = VectorXd::Zero(dim_);
  VectorXd rk = b - A_ * xk;
  VectorXd dk = M_ == nullptr ? rk : M_(rk);
  VectorXd rk_p(dim_);
  VectorXd Adk(dim_);

  double dkAdk, alpha, beta, residual;

  bool convergence = false;
  size_t itrs = 0;
  
  if(M_ != nullptr){
    for(size_t i = 0; i < max_itrs_; ++i){
      Adk = A_ * dk;
      dkAdk = dk.dot(Adk);
      alpha = dk.dot(rk) / dkAdk;
      xk += alpha * dk;
      rk -= alpha * Adk;
      residual = fabs(rk.maxCoeff());
      if(residual < tol_){
        convergence = true;
        itrs = i;
        break;
      }

      rk_p = M_(rk);
      beta = rk_p.dot(Adk) / dkAdk;
      dk = rk_p - beta * dk;
    }
  }else{
    double rkrk = rk.dot(rk), rkrk_new;
    for(size_t i = 0; i < max_itrs_; ++i){
      Adk = A_ * dk;
      alpha = rkrk / dk.dot(Adk);
      xk += alpha * dk;
      rk -= alpha * Adk;
      residual = fabs(rk.maxCoeff());
      if(residual < tol_){
        convergence = true;
        itrs = i;
        break;
      }
      rkrk_new = rk.dot(rk);
      beta =rkrk_new / rkrk;
      dk = rk + beta * dk;
      rkrk = rkrk_new;
    }
  }
  if(convergence)
    cout << "PCG converge with " << itrs
         <<" steps and the residual is " << residual << endl;
  else
    cout << "PCG exceed max iterations and the rediual now is " << residual << endl;
  
  return xk;
}

    

}
