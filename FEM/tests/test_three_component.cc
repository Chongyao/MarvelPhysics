#include "FEM/src/elas_energy.h"
#include "FEM/src/mass_matrix.h"
#include <cmath>
#include <iostream>
using namespace std;
using namespace Eigen;
using namespace marvel;

using CSTTT = linear_csttt<double, 3>;
using BASIS = basis_func<double, 3, 1, 4>;
using QDRT = quadrature<double, 3, 1, 4>;
double other_val(const Ref<MatrixXd> X,  const Ref<Matrix3d> D, const double volume, const double lam, const double miu){
  const double  tt1 = -X(1,1), tt2 = -X(2,1) ,tt3 = -X(3,1);
  return volume *(0.5*lam*pow((0.5*(2*D(3,3)*(X(3,4)+tt3)+2*D(2,3)*(X(3,3)+tt3)+2*D(1,3)*(X(3,2)+tt3))+0.5*(2*(X(2,4)+tt2)*D(3,2)+2*D(2,2)*(X(2,3)+tt2)+2*D(1,2)*(X(2,2)+tt2))+0.5*(2*(X(1,4)+tt1)*D(3,1)+2*(X(1,3)+tt1)*D(2,1)+2*D(1,1)*(X(1,2)+tt1))-3),2)+miu*(pow((D(3,3)*X(3,4)+D(2,3)*X(3,3)-X(3,1)*D(3,3)+D(1,3)*X(3,2)-D(2,3)*X(3,1)-D(1,3)*X(3,1)-1),2)+0.5*pow((D(3,2)*X(3,4)+D(2,2)*X(3,3)+X(2,4)*D(3,3)-X(2,1)*D(3,3)+D(1,2)*X(3,2)-X(3,1)*D(3,2)-D(2,2)*X(3,1)-D(1,2)*X(3,1)+D(2,3)*X(2,3)-X(2,1)*D(2,3)+D(1,3)*X(2,2)-D(1,3)*X(2,1)),2)+0.5*pow((D(3,1)*X(3,4)+D(2,1)*X(3,3)+X(1,4)*D(3,3)-X(1,1)*D(3,3)+D(1,1)*X(3,2)-D(3,1)*X(3,1)-D(2,1)*X(3,1)-D(1,1)*X(3,1)+X(1,3)*D(2,3)-X(1,1)*D(2,3)+X(1,2)*D(1,3)-X(1,1)*D(1,3)),2)+pow((X(2,4)*D(3,2)-X(2,1)*D(3,2)+D(2,2)*X(2,3)+D(1,2)*X(2,2)-X(2,1)*D(2,2)-D(1,2)*X(2,1)-1),2)+0.5*pow((X(1,4)*D(3,2)-X(1,1)*D(3,2)+X(2,4)*D(3,1)-X(2,1)*D(3,1)+D(2,1)*X(2,3)+D(1,1)*X(2,2)+X(1,3)*D(2,2)-X(1,1)*D(2,2)-D(2,1)*X(2,1)-D(1,1)*X(2,1)+D(1,2)*X(1,2)-X(1,1)*D(1,2)),2)+pow((X(1,4)*D(3,1)-X(1,1)*D(3,1)+X(1,3)*D(2,1)-X(1,1)*D(2,1)+D(1,1)*X(1,2)-D(1,1)*X(1,1)-1),2)));
}
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

          
