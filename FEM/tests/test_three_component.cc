#include "FEM/src/elas_energy.h"
#include "FEM/src/mass_matrix.h"
#include <cmath>
#include <iostream>
using namespace std;
using namespace Eigen;
using namespace marvel;

using CSTTT = stvk<double, 3, 3>;
using BASIS = basis_func<double, 3,3, 1, 8>;
using QDRT = quadrature<double, 3, 2, 8>;

int main(){
  VectorXi one_cube(8);
  one_cube << 0, 1, 2, 3, 4, 5, 6, 7;

  MatrixXd nods(3, 8);
  nods <<
      -1, 1, -1, 1, -1, 1, -1, 1,
      -1, -1, 1, 1, -1, -1, 1, 1,
      -1, -1, -1, -1, 1, 1, 1, 1;

  MatrixXd deformed = nods * 2;

  const QDRT q = QDRT();
  const size_t qdrt_num = q.PNT_.cols();
  
  for(size_t qdrt_id = 0; qdrt_id < qdrt_num; ++qdrt_id){
    Matrix<double, 8, 3> Dphi_Dxi;
    BASIS::calc_Dphi_Dxi(q.PNT_.col(qdrt_id), nods.data(), Dphi_Dxi);
    
    double jac_det;
    Matrix3d Dm_inv;
    BASIS::calc_InvDm_Det(Dphi_Dxi, nods.data(), jac_det, Dm_inv);

    Matrix3d def_gra;
    BASIS::get_def_gra(Dphi_Dxi, deformed.data(), Dm_inv, def_gra);
    cout << def_gra << endl;
    getchar();
  }


  
  return 0;
}



          
