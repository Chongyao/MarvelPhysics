#ifndef DATA_STR_H
#define DATA_STR_H
#include <Eigen/Dense>
#include <Eigen/Sparse>
namespace marvel{
template <typename T, size_t dim_>
class dat_str_core{
 public:
  virtual ~dat_str_core() {}
  dat_str_core(const size_t& dof);

  int set_zero();
  //!!!!!!!WARNING!!!!!!!!!:   reserve enough space 
  int hes_reserve(const Eigen::VectorXi& nnzs);
  int hes_compress();
  int hes_add_diag(const size_t& time);
  
  int save_val(const T& val);
  int save_gra(const Eigen::Matrix<T, Eigen::Dynamic, 1>& gra);
  int save_gra(const size_t& pos, const Eigen::Matrix<T, dim_, 1>& point_gra);
  int save_gra(const size_t& pos, const T& one_gra);
  

  int save_hes(const size_t&m, const size_t& n, const Eigen::Matrix<T, dim_, dim_>& loc_hes);
  int save_hes(const size_t& row, const size_t& col, const T& value);

  int set_zero_after_pre_compute();

  
  const T get_val() const;
  const Eigen::Matrix<T, Eigen::Dynamic, 1>& get_gra() const;
  const Eigen::SparseMatrix<T>& get_hes()const;


  //TODO:add Perfect Forwardincg
  

 private:
  const size_t dof_;
  T val_;
  Eigen::Matrix<T, Eigen::Dynamic, 1> gra_;
  Eigen::SparseMatrix<T> hes_;
  Eigen::Matrix<T, Eigen::Dynamic, 1> all_one_;
  bool if_pre_compute_hes_{false};

};

template class dat_str_core<double, 3>;
}
#endif
