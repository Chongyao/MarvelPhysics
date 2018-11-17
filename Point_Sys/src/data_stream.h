#ifndef DATA_STREAM_H
#define DATA_STREAM_H
#include <Eigen/Core>
#include <Eigen/SparseCore>
namespace marvel{


struct energy_dat{
 public:
  energy_dat(const size_t &dim);
  // ~energy_dat();

  //data stream
  const size_t dim_;
  
  // energy_dat(const energy_dat &other);
  //save data by element
  int save_ele_def_gra(const size_t &ele_id, const Eigen::MatrixXd &ele_mat);
  int save_ele_inv_all(const size_t &ele_id, const Eigen::MatrixXd &ele_mat);
  int save_ele_gra(const size_t &ele_id, const Eigen::MatrixXd &ele_mat);
  int save_ele_hes(const size_t &ele_id, const Eigen::MatrixXd &ele_mat);
  int save_ele_strain(const size_t &ele_id, const Eigen::MatrixXd &ele_mat);
  int save_ele_stress(const size_t &ele_id, const Eigen::MatrixXd &ele_mat);
  int save_ele_pre_F(const size_t &ele_id, const Eigen::MatrixXd &ele_mat);
  int save_val(const double val);
  int set_zero();
  //get data by element

  // Eigen::MatrixXd& ele_def_gra(const size_t &ele_id);
  // Eigen::MatrixXd& ele_inv_all(const size_t &ele_id);
  // Eigen::MatrixXd& ele_gra(const size_t &ele_id);
  // Eigen::MatrixXd& ele_hes(const size_t &ele_id);
 // private:
  Eigen::MatrixXd def_gra_;
  Eigen::MatrixXd sigma_w_points_;
  Eigen::MatrixXd inv_A_all_;
  Eigen::MatrixXd gra_;
  Eigen::MatrixXd strain_;
  Eigen::MatrixXd stress_;
  Eigen::MatrixXd pre_F_;
  Eigen::VectorXd ela_val_;
  Eigen::VectorXd vol_val_;
  double Val_;
  std::vector<Eigen::Triplet<double>> hes_trips;
  Eigen::SparseMatrix<double> hes_;
  // const Eigen::MatrixXd& ele_mat(const size_t &ele_id, const size_t &rows, const size_t &cols, Eigen::MatrixXd &mat);
  int save_ele_mat(const size_t &ele_id, const size_t &rows, const Eigen::MatrixXd &ele_mat, Eigen::MatrixXd &whole_mat, bool if_plus);
};

}
#endif
