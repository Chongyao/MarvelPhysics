#include "data_str_core.h"

using namespace std;
using namespace Eigen;
namespace marvel{

template <typename T, size_t dim_>
dat_str_core<T, dim_>::dat_str_core(const size_t& dof):dof_(dof), val_(0), gra_(dof), hes_(dof, dof), all_one_(Matrix<T, Dynamic, 1>::Ones(dof)){
  set_zero();
}

template <typename T, size_t dim_>
int dat_str_core<T, dim_>::set_zero(){
  val_ = 0;
  gra_.setZero();
  hes_.setZero();
  return 0;
}

template <typename T, size_t dim_>
int dat_str_core<T, dim_>::hes_reserve(const VectorXi& nnzs){
  hes_.reserve(nnzs);
  return 0;
}

template <typename T, size_t dim_>
int dat_str_core<T, dim_>::hes_compress(){
  hes_.makeCompressed();
  return 0;
}

template <typename T, size_t dim_>
int dat_str_core<T, dim_>::hes_add_diag(const size_t& time){
  hes_ += (time * all_one_).asDiagonal();
  return 0;
}


template <typename T, size_t dim_>
int dat_str_core<T, dim_>::save_val(const T& val){

  #pragma omp atomic
  val_ += val;
  return 0;
}

template <typename T, size_t dim_>
int dat_str_core<T, dim_>::save_gra(const Eigen::Matrix<T, Eigen::Dynamic, 1>& gra){
  gra_ += gra;
  return 0;
}
template <typename T, size_t dim_>
int dat_str_core<T, dim_>::save_gra(const size_t& pos, const Eigen::Matrix<T, dim_, 1>& one_gra){
  for(size_t d = 0; d < dim_; ++d){
    #pragma omp atomic
    gra_(dim_ * pos + d) += one_gra(d);      
  }
  return 0;
}


template <typename T, size_t dim_>
int dat_str_core<T, dim_>::save_hes(const size_t&m, const size_t& n, const Eigen::Matrix<T, dim_, dim_>& loc_hes){
  for(size_t row = 0; row < dim_; ++row){
    for(size_t col = 0; col < dim_; ++col){
      if(loc_hes(row, col)){
        #pragma omp atomic
        hes_.coeffRef(m * dim_ + row, n * dim_ + col) += loc_hes(row, col);
      }
                     
    }
  }
  return 0;
}

    
template <typename T, size_t dim_>
int dat_str_core<T, dim_>::save_hes(const size_t& row, const size_t& col, const T& value){
  #pragma omp atomic
  hes_.coeffRef(row, col) += value;
  return 0;
}

template <typename T, size_t dim_>
const T dat_str_core<T, dim_>::get_val() const{
  return val_;
}

template <typename T, size_t dim_>
const Eigen::Matrix<T, Eigen::Dynamic, 1>& dat_str_core<T, dim_>::get_gra() const{
  return gra_;
}

template <typename T, size_t dim_>
const Eigen::SparseMatrix<T>& dat_str_core<T, dim_>::get_hes()const{
  return hes_;
}
}//namespace

