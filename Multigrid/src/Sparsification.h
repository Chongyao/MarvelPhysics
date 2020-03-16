#ifndef MULTIGRID_SPARSIFICATION_ZCY
#define MULTIGRID_SPARSIFICATION_ZCY
#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <vector>
#include <memory>
#include <unordered_set>
#include <set>

namespace marvel{

//TODO: decouplt and adjc_graph and sparsification

class Adjc_graph{
  template<typename T>
  using VS = std::vector<std::shared_ptr<T>>;
  using TPL = Eigen::Triplet<double>;
 public:
  Adjc_graph(const Eigen::MatrixXd& L);
  Adjc_graph(const Eigen::SparseMatrix<double>& L);
  int build_mat_from_graph(Eigen::MatrixXd& L)const;
  int build_mat_from_graph(std::vector<TPL>& trips)const;

  int build_reordered_mat_from_graph(std::vector<TPL>& trips);
   
  size_t num_coarse_{0};
  int Sparsification();
 private:
  const size_t dof_;
  VS<std::unordered_set<size_t>> vertices_;
  VS<TPL> edges_;
  Eigen::VectorXd dig_vals_;

  enum class mark_state{fine, coarse, unmarked};
  std::vector<mark_state> labels_;
  std::set<size_t> unmarked_vertices_;


  
 private:
  void init();
  void sparsify_one_edge(const size_t edge_id);
  void compensate_one_edge(const size_t edge_id_ik, const size_t edge_id_jk, const double& w);
  bool is_connect(const size_t i, const size_t j, size_t& edge_ij)const;
  int sparsify_one_tri(const size_t edge_id_i, const size_t edge_id_j, const size_t edge_id_k, size_t& sparsified_edge_id); 

};

int Schur_complement(const Eigen::SparseMatrix<double>& L, const size_t& coarse_num,
                     Eigen::SparseMatrix<double>& topleft,
                     Eigen::SparseMatrix<double>& bottomright);

}
#endif
