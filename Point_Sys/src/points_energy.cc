#include "points_energy.h"
#include "get_nn.h"
using namespace std;
using namespace Eigen;
namespace marvel{

point_sys::point_sys(const MatrixXd  &points_, const double &rho_, const double &vol_all_):
    points(points_), rho(rho_), vol_all(vol_all_){
  calc_NNN(points, NN, sup_radi);
  
}



}


