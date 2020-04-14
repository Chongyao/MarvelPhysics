#ifndef MARVEL_CONSTRAINED_NEWTON
#define MARVEL_CONSTRAINED_NEWTON
#include "implicit_euler.h"
#include "mprgp_solver/src/collision_proxy.h"
namespace marvel{
template<typename T, size_t dim_>
class constrained_newton : public newton_iter<T, dim_>{
 public:
  using SMP_TYPE = typename dat_str_core<T, dim_>::SMP_TYPE;
  using coll_type = typename chaos::collision::collision_proxy<T>;
  constrained_newton(
      std::shared_ptr<coll_type>& coll,
      std::shared_ptr<dat_str_core<T, dim_>>& dat_str,
      std::shared_ptr<Functional<T, dim_>>& energy,
      const size_t max_iter = 20,const T tol = 1e-4,
      const bool if_pre_compute_hes = false,
      const bool if_hes_constant = false)
      :newton_iter<T, dim_>(dat_str, energy, max_iter, tol, if_pre_compute_hes, false, if_hes_constant), coll_(coll){}

  //TODO: need to adapt to multi body
  int linear_solver(const SMP_TYPE*A, const Eigen::Matrix<T, -1, 1>&b, Eigen::Matrix<T, -1, 1>& solution) override{
    if(coll_->get_last_coll_num()){
      vector<vector<T>> sol_coll;
      sol_coll = coll_->response(*A, b, false);
      copy(sol_coll[0].begin(), sol_coll[0].end(), solution.data());
    }else{
      newton_iter<T, dim_>::linear_solver(A, b, solution);
    }
    return 0;
  }
 private:
  std::shared_ptr<coll_type>& coll_;
  
};
}
#endif
