#ifndef FEM_QUADRATURE
#define FEM_QUADRATURE
#include <array>
namespace marvel{
template<typename T, size_t dim_, size_t num_quad_>
class gaus_quad{
 public:
  static std::array<T, num_quad_> PNT_;
  static std::array<T, num_quad_> WGT_;
};
}
#endif
