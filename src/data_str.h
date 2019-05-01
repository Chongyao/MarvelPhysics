#ifndef DATA_STR_H
#define DATA_STR_H
#include <Eigen/Sparse>
namespace marvel{
template <typename T>
class dat_str_core{
 public:
  dat_str_core(const size_t dim):val_(0), gra_(Eigen::VectorXd::Zero(3 * dim)), hes_trips(0){}
  
  T val_;
  Eigen::VectorXd gra_;
  std::vector<Eigen::Triplet<T>> hes_trips;
};
}
#endif
