#ifndef MARVEL_TENSOR
#define MARVEL_TENSOR
#include <Eigen/Dense>
#include <array>
namespace marvel{


template<typename T, size_t i, size_t j, size_t k, size_t l>
class fourth_tensor{
  using tensor_type = std::array<std::array<Eigen::Matrix<T, k, l> , j> , i>;
 public:
  fourth_tensor(){
    for(size_t i_id = 0; i_id < i; ++i_id){
      for(size_t j_id = 0; j_id < j; ++j_id){
        tensor_[i][j] = Eigen::Matrix<T, k, l>::Zero();
      }
    }
  }
  
  Eigen::Matrix<T, k, l>& operator()(const size_t row, const size_t col){
    return tensor_[row][col];
  }
  

  
  tensor_type operator * (const Eigen::Ref<Eigen::Matrix<T, i, j>>& mat){
    fourth_tensor res;
    for(size_t res_row = 0; res_row < i; ++res_row){
      for(size_t res_col = 0; res_col < j; ++res_col){
        for(size_t e_id = 0; e_id < i; ++e_id){
          res(res_row, res_col) += tensor_[res_row][e_id] * mat(e_id,res_col);
        }
      }
    }
    return res;
  }
  tensor_type mat_mult_this (const Eigen::Ref<Eigen::Matrix<T, i, j>>& mat){
    fourth_tensor res;
    for(size_t res_row = 0; res_row < i; ++res_row){
      for(size_t res_col = 0; res_col < j; ++res_col){
        for(size_t e_id = 0; e_id < i; ++e_id){
          res(res_row, res_col) += mat(res_row, e_id) * tensor_[e_id][res_col];
        }
      }
    }
    return res;
  }


  tensor_type operator + (const tensor_type& other){
    fourth_tensor res;
    #pragma omp parallel
    for(size_t res_row = 0; res_row < i; ++res_row){
      for(size_t res_col = 0; res_col < j; ++res_col){
        res(res_row, res_col) =  tensor_[res_row][res_col] + other(res_row, res_col);
      }
    }
  }
  void Flatten(Eigen::Matrix<T, i * j, k * l>& flatted){
    for(size_t row_out = 0; row_out < i; ++row_out){
      for(size_t col_out = 0; col_out < j; ++col_out){
        Eigen::Map<Eigen::Matrix<T, 1, k * l>> vec(tensor_[row_out][col_out].data());
        flatted.row(col_out * i + row_out) = vec;
      }
    }
  }

  tensor_type tensor_;
  
};
}
#endif
