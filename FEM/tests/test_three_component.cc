#include "FEM/src/elas_energy.h"
#include "FEM/src/mass_matrix.h"
#include <cmath>
#include <iostream>
using namespace std;
using namespace Eigen;
using namespace marvel;

using CSTTT = stvk<double, 3>;
using BASIS = basis_func<double, 3, 1, 4>;
using QDRT = quadrature<double, 3, 1, 4>;

int main(){

  Matrix3d F;

  MatrixXd rest(3, 4), deformed(3, 4);
  rest << 2,0,0,-1,
      0,3,0,-1,
      0,0,5,-1;

  deformed << 1,2,3,4,5,6,7,12,11,10,9,8;

  cout <<"rest is " << endl << rest << endl << "deformed is "<<endl << deformed << endl;

  Matrix3d Dm, Dm_inv;

  for(int i = 0; i < 3; i++)
    Dm.col(i) = rest.col(i) - rest.col(3);
  Dm_inv = Dm.inverse();

  double jac_det = Dm.determinant();
  BASIS::get_def_gra(QDRT::PNT_, deformed.data(), Dm_inv, F);
  cout << "def gra is " << endl << F <<endl <<  "jac det is " << jac_det/6.0 << endl;
  Matrix<double, 9, 12> Ddef_Dx;
  BASIS::get_Ddef_Dx(QDRT::PNT_, deformed.data(), rest.data(), F, Ddef_Dx);
  cout << "Ddef_Dx is " <<endl << Ddef_Dx << endl;



  constexpr double Ym = 4000, Pr = 0.45;
  double mu, lam;
  compute_lame_coeffs(Ym, Pr, mu, lam);
    // lam= 0;

  double val = CSTTT::val(F, lam, mu);
  cout  << "val is " <<  val * jac_det / 6 << endl;

  auto gra = CSTTT::gra(F, lam, mu);
  cout << "f based gra is " << endl << gra  << endl;
  cout << "x based gra is " <<endl << Ddef_Dx.transpose() * gra * jac_det / 6 << endl; 

  auto hes = CSTTT::hes(F, lam, mu);
  cout << "f based hes is " << endl << hes *jac_det / 6<< endl;
  cout << "x based hes us "<< endl <<    Ddef_Dx.transpose() * hes * Ddef_Dx *  jac_det / 6.0 << endl;

  
 


  //test Gra
  // Matrix<double, 12, 1> gra_x_based = Ddef_Dx.transpose() * gra;
  // cout << "x based gra is " <<endl << gra_x_based << endl;
  
  // auto hes_x_based = Ddef_Dx.transpose() * hes * Ddef_Dx;
  // cout << "x based hes is "<< endl << hes_x_based << endl;

  // cout << "test "<<endl;
  // Matrix3d D;
  // for(size_t i = 0; i < 3; ++i){
  //   D.col(i) = rest.col(i + 1) - rest.col(0);
  // }
  // cout << D << jac_det / 6.0 << " " << lam << " " << mu<<endl;
  
  // cout << other_val(deformed, D, jac_det / 6.0, lam, mu);
  
  
}

// real::D(3,3)=reshape((/-2,-2,-3,3,0,-1,0,5,-1),(3,3));

          
