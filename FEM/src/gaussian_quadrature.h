#ifndef FEM_QUADRATURE
#define FEM_QUADRATURE
#include <array>
#include <Eigen/Core>

namespace marvel{
using namespace Eigen;
template<typename T, size_t dim_, size_t num_qdrt_>
class quadrature{
 public:
  static Eigen::Matrix<T, dim_, num_qdrt_> PNT_;
  static std::array<T, num_qdrt_> WGT_;
};


template<typename T>
class tri_one_qdrt : public quadrature<T, 2, 1>{
 public:
  tri_one_qdrt(){
    this-> PNT_ = Matrix<T, 2, 1>::Ones() / 3;
    this-> WGT_[0] = 1.0;
  }
};

template<typename T>
class tet_one_qdrt : public  quadrature<T, 3, 1>{
 public:
  tet_one_qdrt(){
    this->PNT_ = Matrix<T, 3, 1>::Ones() /4;
    this->WGT_[0]  = 1.0;
  }
  
};


}
#endif
