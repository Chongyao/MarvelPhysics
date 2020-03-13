#include "inverse_isoparametric_map.h"
#include "FEM/src/basis_func.h"
#include <iostream>
namespace marvel{
using namespace std;
using namespace Eigen;

Vector3d inverse_isoparametric_hex(const Eigen::VectorXd& Xp, const double* X){
  const Map<const Matrix<double, 3, 8>> nods(X);
  

  Vector3d xi_now = Vector3d::Zero();
  Matrix<double, 8 ,1> basis_value = shape_func<double, 3, 1, 8>::calc_basis_value(xi_now, X);

  Vector3d X_now = nods * basis_value;
  Matrix<double, 8, 3> Dphi_Dxi;
  
  auto diff = 999.0;
  do{
    shape_func<double, 3, 1, 8>::calc_Dphi_Dxi(xi_now, X, Dphi_Dxi);
    Matrix3d J = nods * Dphi_Dxi;
    
    xi_now += J.lu().solve(Xp - X_now);
    
    basis_value = shape_func<double, 3, 1, 8>::calc_basis_value(xi_now, X);
    X_now = nods * basis_value;
    diff = (Xp - X_now).maxCoeff();
  }while(fabs(diff) > 1e-3);
  
  return xi_now;
}

}

