#ifndef NODAL_BASIS_OPTIMIZER_H
#define NODAL_BASIS_OPTIMIZER_H

#include <unordered_set>
#include <boost/property_tree/ptree.hpp>
#include <zjucad/matrix/matrix.h>
#include <Eigen/Sparse>

namespace jtf { namespace mesh {
class edge2cell_adjacent;
class N_ring_face_at_point;
}}

namespace bigbang {

using mati_t=zjucad::matrix::matrix<size_t>;
using matd_t=zjucad::matrix::matrix<double>;

class nodal_basis_optimizer
{
public:
  nodal_basis_optimizer(const mati_t &quad_H, const mati_t &quad_h, const matd_t &nods_h,
                        const matd_t &lame_h, const size_t max_rd);
  int optimize(const boost::property_tree::ptree &pt, Eigen::SparseMatrix<double> &N);
  int optimize_nodal_basis(const size_t pid, const boost::property_tree::ptree &pt,
                           std::vector<Eigen::Triplet<double>> &trips);
private:
  const size_t one_to_many_;
  const mati_t &quad_H_, &quad_h_;
  const matd_t &nods_h_, &lame_h_;
  std::unordered_set<size_t> bnd_node_h_;
  std::shared_ptr<jtf::mesh::edge2cell_adjacent> e2c_;
  std::shared_ptr<jtf::mesh::N_ring_face_at_point> p2f_;
};

}

#endif
