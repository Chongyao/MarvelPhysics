#ifndef FEM_QUADRATURE
#define FEM_QUADRATURE
#include <array>
#include <Eigen/Core>

namespace marvel{
template<typename T, size_t dim_, size_t num_qdrt_>
class quadrature{
 public:
  static Eigen::Matrix<T, dim_, num_qdrt_> PNT_;
  static std::array<T, num_qdrt_> WGT_;
};
}
#endif
