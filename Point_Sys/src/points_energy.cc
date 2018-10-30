#include "points_energy.h"
#include "get_nn.h"
using namespace std;
using namespace Eigen;
namespace marvel{

point_sys::point_sys(const MatrixXd  &points, const double &rho, const double &vol_all, const size_t &nearest_num):
    points_(points), rho_(rho), vol_all_(vol_all), dim_(points_.cols()), nearest_num_(nearest_num),SH_(points, nearest_num){

  sup_radi_ = SH_.get_sup_radi();
  //init
  double mass_total = rho*vol_all;
  double mass_sigma = (sup_radi_/3).array().cube().sum();
  scal_fac_ = mass_total/mass_sigma;


}

size_t point_sys::Nx() const{
  return dim_;
}

int point_sys::calc_rhoi_vi(const double *x, const VectorXd &rho_i, const VectorXd &vol_i){
  
}



}


