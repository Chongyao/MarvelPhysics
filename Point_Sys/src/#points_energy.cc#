#include "points_energy.h"
#include "get_nn.h"
using namespace std;
using namespace Eigen;
namespace marvel{

point_sys::point_sys(const MatrixXd  &points_, const double &rho_, const double &vol_all_):
    points(points_), rho(rho_), vol_all(vol_all_), dim(points_.cols()){
  calc_NNN(points, NN, sup_radi, 10);

  //init
  double mass_total = rho*vol_all;
  double mass_sigma = (sup_radi/3).array().cube().sum();
  scal_fac = mass_total/mass_sigma;
  
}

size_t point_sys::Nx() const{
  return dim;
}

int point_sys::calc_rhoi_vi(const double *x, const VectorXd &rho_i, const VectorXd &vol_i){
  
}



}


