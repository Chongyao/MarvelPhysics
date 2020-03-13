#include "inverse_isoparametric_map.h"
namespace marvel{
using namespace std;
using namespace Eigen;

int inverse_isoparametric_hex(const Eigen::VectorXd& Xp, const double* X, Vector3d& iso_coorxsa){
  const Map<const Matrix<double, 3, 8>> nods(X);
  
  Vector3d X_now = Vector3d::Zero();
  Vector3d xi_now = Vector3d::Zero();
  Matrix<double, 8 ,1> basis_value;
  Matrix<double, 8, 3> Dphi_Dxi;

  shape_func<double, 3, 1, 8>::calc_basis_value(xi_now, X, basis_value);
  X_now = nods * basis_value;
  
  auto diff = 999.0;
  do{
    shape_func<double, 3, 1, 8>::calc_Dphi_Dxi(xi_now, X, Dphi_Dxi);
    Matrix3d J = nods * Dphi_Dxi;
    
    xi_now += J.lu().solve(Xp - X_now);
    
    shape_func<double, 3, 1, 8>::calc_basis_value(xi_now, X, basis_value);
    X_now = nods * basis_value;
    diff = (Xp - X_now).maxCoeff();
  }while(fabs(diff) > 1e-3);
  
  return 0;
}

}

