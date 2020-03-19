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
template<typename T>
  using VS = std::vector<std::shared_ptr<T>>;
  using TPL = Eigen::Triplet<double>;

class Sparsify;

class Adjc_graph{
  
 public:
  Adjc_graph(const Eigen::MatrixXd& L);
  Adjc_graph(const Eigen::SparseMatrix<double>& L);
  int build_mat_from_graph(Eigen::MatrixXd& L)const;
  int build_mat_from_graph(std::vector<TPL>& trips)const;

  int build_reordered_mat_from_graph(const Eigen::VectorXi& perm_inv, std::vector<TPL>& trips);

 private:
  const size_t dof_;
  VS<std::unordered_set<size_t>> vertices_;
  VS<TPL> edges_;
  Eigen::VectorXd dig_vals_;

  friend class Sparsify;

  
 private:
  void init();
  bool is_connect(const size_t i, const size_t j, size_t& edge_ij)const;

};

class Sparsify{
 public:
  Sparsify(const size_t& dof);
  void sparsify_and_compensate(Adjc_graph& graph);
  void post_coloring(const Adjc_graph& graph);
  void reorder_coarse_and_fine();
  
  Eigen::VectorXi get_perm_inv()const;
  
  size_t num_coarse_{0};
 private:
  const size_t dof_;
  
  enum class mark_state{fine, coarse, unmarked};
  std::vector<mark_state> labels_;
  std::set<size_t> unmarked_vertices_;

  Eigen::VectorXi perm_vec_;
  Eigen::VectorXi perm_vec_inv_;
  bool has_reordered{false};

 private:
  void sparsify_one_edge(Adjc_graph& graph, const size_t edge_id);
  void compensate_one_edge(Adjc_graph& graph, const size_t edge_id_ik, const size_t edge_id_jk, const double& w);
  void sparsify_one_tri(Adjc_graph& graph, const size_t edge_id_i, const size_t edge_id_j, const size_t edge_id_k, size_t& sparsified_edge_id); 

  
  
};


int Schur_complement(const Eigen::SparseMatrix<double>& L, const size_t& coarse_num,
                     Eigen::SparseMatrix<double>& topleft,
                     Eigen::SparseMatrix<double>& bottomright);

}
#endif
