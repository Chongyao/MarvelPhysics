#include "basic_energy.h"
#include <Eigen/SparseCore>

using namespace std;
using namespace Eigen;

namespace marvel{
position_constraint::position_constraint(const double &w, const vector<size_t> &cons, const size_t dim):w_(w), cons_(cons), dim_(dim){}
int position_constraint::Gra(const double *disp, energy_dat &dat_str){
  Map<const MatrixXd> _disp(disp, 3, dim_);
  for(auto iter_c = cons_.begin(); iter_c != cons_.end(); ++iter_c)
    dat_str.save_ele_gra(*iter_c, -2.0 * w_ * _disp.col(*iter_c));
  return 0;
}

int position_constraint::Hes(const double *disp, energy_dat &dat_str){
  for(size_t i = 0; i < cons_.size(); ++i)
    for(size_t j = 0; j < 3; ++j)
      dat_str.hes_trips.push_back(Triplet<double>(cons_[i]*3 + j, cons_[i]*3 + j, 2 * w_));
  return 0;
}

}//namespace marvel
