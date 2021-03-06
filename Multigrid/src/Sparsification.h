#ifndef MULTIGRID_SPARSIFICATION_ZCY
#define MULTIGRID_SPARSIFICATION_ZCY
#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <vector>
#include <memory>
#include <unordered_set>
#include <set>
#include <iostream>
namespace marvel{

//TODO: decouplt and adjc_graph and sparsification
template<typename T>
  using VS = std::vector<std::shared_ptr<T>>;
  using TPL = Eigen::Triplet<double>;

class Sparsify;

class Adjc_graph{
  
 public:
  //option 1 : positive w, -1 : negative w, 0 : all the w
  Adjc_graph(const Eigen::SparseMatrix<double>& L, const int option = 1);
  int build_mat_from_graph(Eigen::MatrixXd& L)const;
  int build_mat_from_graph(std::vector<TPL>& trips)const;

  int build_reordered_mat_from_graph(const Eigen::VectorXi& perm_inv, std::vector<TPL>& trips);

 private:
  const size_t dof_;
  VS<std::unordered_set<size_t>> vertices_;
  VS<TPL> edges_;

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
  Eigen::VectorXi get_perm() const;
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

int add_dig_vals(const Eigen::SparseMatrix<double>& L, std::vector<TPL>& trips, const Eigen::VectorXi* perm = nullptr);
int Schur_complement(const Eigen::SparseMatrix<double>& L, const size_t& coarse_num,
                     Eigen::SparseMatrix<double>& topleft,
                     Eigen::SparseMatrix<double>& bottomright);


template<int mode_L, int mode_C, int mode_S>
int Schur_complement(const size_t& coarse_num,
                     const Eigen::SparseMatrix<double, mode_L>& L, 
                     Eigen::SparseMatrix<double, mode_C>& topleft,
                     Eigen::SparseMatrix<double, mode_S>& S){

  const size_t fine_num = L.rows() - coarse_num;

  Eigen::SparseMatrix<double, Eigen::RowMajor>
      L_cc = L.topLeftCorner(coarse_num, coarse_num),
      L_fc = L.bottomLeftCorner(fine_num, coarse_num),
      L_cf = L.topRightCorner(coarse_num, fine_num),
      L_ff = L.bottomRightCorner(fine_num, fine_num);

  
  const Eigen::VectorXd L_ff_diag = L_ff.diagonal();

  {//Lff_inv * Lfc
    #pragma omp parallel for
    for (int k=0; k<L_fc.outerSize(); ++k)
      for (decltype(L_fc)::InnerIterator it(L_fc,k); it; ++it)
        it.valueRef() /= L_ff_diag(it.row()) ;
  }

  {//set S
    S.resize(L.rows(), coarse_num);
    S.setZero();
    S.topRows(coarse_num) = Eigen::MatrixXd::Identity(coarse_num, coarse_num).sparseView();
    S.bottomRows(fine_num) = -L_fc;
  }

  {//set LH
    decltype(L_cc) L_comp = L_cf * L_fc;
    topleft = L_cc - L_comp;
  }
  return 0;
}

}
#endif
