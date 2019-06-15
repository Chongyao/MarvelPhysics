#include "nodal_basis_optimizer.h"

#include <zjucad/matrix/io.h>
#include <jtflib/mesh/mesh.h>

#include "config.h"
#include "geom_util.h"
#include "io.h"
#include "energy.h"
#include "phx_util.h"
#include "util.h"
#include "grad_check.h"

using namespace std;
using namespace zjucad::matrix;
using namespace Eigen;
using namespace jtf::mesh;

#define DEBUG_FLAG 0

namespace bigbang {

class nodal_basis_opt_energy_2d : public Functional<double>
{
public:
  nodal_basis_opt_energy_2d(const SparseMatrix<double> &K, const double w=1.0)
      : K_(K), dim_(K.cols()/2), w_(w) {
    ASSERT(K_.cols()%2 == 0);

    const size_t verts_num_ = dim_;
    Gx_.resize(2*verts_num_, verts_num_); {
      vector<Triplet<double>> trips;
      for (size_t i = 0; i < verts_num_; ++i) {
        trips.push_back(Triplet<double>(2*i+0, i, 1));
      }
      Gx_.setFromTriplets(trips.begin(), trips.end());
    }

    Gy_.resize(2*verts_num_, verts_num_); {
      vector<Triplet<double>> trips;
      for (size_t i = 0; i < verts_num_; ++i) {
        trips.push_back(Triplet<double>(2*i+1, i, 1));
      }
      Gy_.setFromTriplets(trips.begin(), trips.end());
    }
  }
  size_t Nx() const {
    return dim_;
  }
  int Val(const double *x, double *val) const {
    Eigen::Map<const VectorXd> X(x, dim_);
    *val += 0.5*w_*(X.dot(Gx_.transpose()*K_*Gx_*X)+X.dot(Gy_.transpose()*K_*Gy_*X));
    return 0;
  }
  int Gra(const double *x, double *gra) const {
    Eigen::Map<const VectorXd> X(x, dim_);
    Eigen::Map<VectorXd> G(gra, dim_);
    G += w_*(Gx_.transpose()*K_*Gx_*X+Gy_.transpose()*K_*Gy_*X);
    return 0;
  }
  int Hes(const double *x, vector<Triplet<double>> *hes) const {
    const SparseMatrix<double> H = w_*(Gx_.transpose()*K_*Gx_+Gy_.transpose()*K_*Gy_);
    for (size_t j = 0; j < H.outerSize(); ++j) {
      for (SparseMatrix<double>::InnerIterator it(H, j); it; ++it) {
        hes->push_back(Triplet<double>(it.row(), it.col(), it.value()));
      }
    }
    return 0;
  }
private:
  const size_t dim_;
  const double w_;
  const SparseMatrix<double> &K_;
  SparseMatrix<double> Gx_, Gy_;
};

static shared_ptr<Functional<double>> g_func;

static void eval_val_gra(const int m, const int n, const double *x,
                         double *value, double *jacobian) {
  *value = 0;
  g_func->Val(x, value);
  std::fill(jacobian, jacobian+n, 0);
  g_func->Gra(x, jacobian);
}

static void eval_gra_hes(const int m, const int n, const double *x,
                         double *value, double *jacobian) {
  std::fill(value, value+n, 0);
  g_func->Gra(x, value);

  SparseMatrix<double> H(n, n); {
    vector<Triplet<double>> trips;
    g_func->Hes(x, &trips);
    H.setFromTriplets(trips.begin(), trips.end());
    MatrixXd Hd = MatrixXd(H);
    std::copy(Hd.data(), Hd.data()+Hd.size(), jacobian);
  }  
}

nodal_basis_optimizer::nodal_basis_optimizer(const mati_t &quad_H, const mati_t &quad_h, const matd_t &nods_h,
                                             const matd_t &lame_h, const size_t max_rd)
    : quad_H_(quad_H), quad_h_(quad_h), nods_h_(nods_h), lame_h_(lame_h),
      one_to_many_(quad_h.size(2)/quad_H.size(2)) {
  ASSERT(quad_h.size(2)%quad_H.size(2) == 0);
  const bool show_debug_info = false;
  e2c_.reset(edge2cell_adjacent::create(quad_H_, show_debug_info));

  p2f_.reset(new N_ring_face_at_point);
  p2f_->add_all_faces(quad_H_, *e2c_, max_rd);

  get_surf_mesh_bnd_nodes(quad_h, bnd_node_h_);
}

int nodal_basis_optimizer::optimize_nodal_basis(const size_t pid, const boost::property_tree::ptree &pt,
                                                vector<Triplet<double>> &trips) {
  //-> extract N ring face of pid
  mati_t N_faces, N_nodes; {
    unordered_set<size_t> N_ring_face_of_pid;
    const vector<vector<size_t>> &f_ring = p2f_->p2f_[pid];
    for (size_t i = 0; i < f_ring.size(); ++i) {
      for (size_t j = 0; j < f_ring[i].size(); ++j) {
        if ( f_ring[i][j] != -1 )
          N_ring_face_of_pid.insert(f_ring[i][j]);
      }
    }

    N_faces.resize(N_ring_face_of_pid.size());
    std::copy(N_ring_face_of_pid.begin(), N_ring_face_of_pid.end(), N_faces.begin());

    mati_t face_node = quad_H_(colon(), N_faces);
    std::unordered_set<size_t> node_set(face_node.begin(), face_node.end());
    N_nodes.resize(node_set.size(), 1);
    std::copy(node_set.begin(),  node_set.end(), N_nodes.begin());

    #if DEBUG_FLAG
    cout << "curr pid: " << pid << endl;
    cout << "N ring faces: " << trans(N_faces) << endl;
    cout << "N ring nodes: " << trans(N_nodes) << endl;
    #endif
  }

  //-> extract corresponding fine mesh patch
  mati_t cell, bnd_node, new_to_orig, orig_to_new;
  matd_t nods, lame; {
    cell.resize(quad_h_.size(1), N_faces.size()*one_to_many_);
    lame.resize(lame_h_.size(1), N_faces.size()*one_to_many_);
    for (size_t i = 0; i < N_faces.size(); ++i) {
      const size_t fid = N_faces[i];
      cell(colon(), colon(i*one_to_many_, (i+1)*one_to_many_-1)) = quad_h_(colon(), colon(fid*one_to_many_, (fid+1)*one_to_many_-1));
      lame(colon(), colon(i*one_to_many_, (i+1)*one_to_many_-1)) = lame_h_(colon(), colon(fid*one_to_many_, (fid+1)*one_to_many_-1));
    }
    nods = nods_h_;
    remove_extra_verts(cell, nods, &new_to_orig);

    //-> build global to local mapping
    orig_to_new = -1*ones<size_t>(nods_h_.size(2), 1);
    orig_to_new(new_to_orig) = colon(0, new_to_orig.size()-1);

    //-> find boundary of the patch
    get_surf_mesh_bnd_nodes(cell, bnd_node);

    #if DEBUG_FLAG
    quad_mesh_write_to_vtk("./nodal_patch.vtk", nods, cell);
    cout << "boundary node: " << trans(bnd_node) << endl;
    #endif
  }

  //-> query hessian and dofs need to be constrained
  shared_ptr<Functional<double>> elas_energy;
  elas_energy = build_elastic_energy<quad4_elastic_energy>(pt.get<string>("elas.value"), cell, nods, lame);
  SparseMatrix<double> K(elas_energy->Nx(), elas_energy->Nx()); {
    vector<Triplet<double>> trips;
    elas_energy->Hes(&nods[0], &trips);
    K.setFromTriplets(trips.begin(), trips.end());
  }

  //-> get all constrained verts and generate global to local mapping
  std::unordered_set<size_t> fix_dofs;
  for (auto &p : bnd_node) {
    const size_t orig_id = new_to_orig[p];
    if ( bnd_node_h_.find(orig_id) != bnd_node_h_.end() ) // if it is on the global boundary
      continue;
    fix_dofs.insert(p);
  }
  for (auto &p : N_nodes)  // if it is coarse mesh verts
    fix_dofs.insert(orig_to_new[p]);

  vector<size_t> g2l(nods.size(2)); {
    size_t cnt = 0;
    for (size_t i = 0; i < g2l.size(); ++i) {
      if ( fix_dofs.find(i) != fix_dofs.end() )
        g2l[i] = -1;
      else
        g2l[i] = cnt++;
    }
  }

  #if DEBUG_FLAG
  cout << "# constrained dofs: " << endl;
  std::for_each(fix_dofs.begin(), fix_dofs.end(), [](size_t it){ cout << it << " "; });
  cout << "# top node: " << orig_to_new[pid] << endl;
  #endif

  #if DEBUG_FLAG
  do {
    g_func = make_shared<nodal_basis_opt_energy_2d>(K, 1.0);
    VectorXd x_tmp = VectorXd::Random(nods.size(2));
    int rtn = numeric_grad_check(eval_val_gra, 1, g_func->Nx(), x_tmp.data());
    cout << "# val gra check: " << rtn << endl;
    rtn = numeric_grad_check(eval_gra_hes, g_func->Nx(), g_func->Nx(), x_tmp.data());
    cout << "# gra hes check: " << rtn << endl;
  } while (0);
  #endif

  //-> optimzie
  VectorXd basis_pid = VectorXd::Zero(nods.size(2), 1);
  basis_pid[orig_to_new[pid]] = 1; {
    shared_ptr<Functional<double>> opt_obj;
    opt_obj = make_shared<nodal_basis_opt_energy_2d>(K, 1.0);
  
    VectorXd g = VectorXd::Zero(opt_obj->Nx()); {
      opt_obj->Gra(basis_pid.data(), g.data());
      g *= -1;
    }
    SparseMatrix<double> H(opt_obj->Nx(), opt_obj->Nx()); {
      vector<Triplet<double>> trips;
      opt_obj->Hes(nullptr, &trips);
      H.setFromTriplets(trips.begin(), trips.end());
    }

    rm_spmat_col_row(H, g2l);
    rm_vector_row(g, g2l);
  
    SimplicialLLT<SparseMatrix<double>> solver;
    solver.compute(H);
    ASSERT(solver.info() == Eigen::Success);
    VectorXd dx = solver.solve(g);
    ASSERT(solver.info() == Eigen::Success);

    VectorXd Dx = VectorXd::Zero(opt_obj->Nx(), 1);
    rc_vector_row(dx, g2l, Dx);
    basis_pid += Dx;
  }

  #if DEBUG_FLAG
  matd_t basis_tmp = itr_matrix<const double*>(1, basis_pid.size(), basis_pid.data());
  quad_mesh_write_to_vtk("./nodal_basis_for_patch.vtk", nods, cell, &basis_tmp, "POINT");
  #endif
  
  //-> save optimized basis into interp matrix N
  for (size_t i = 0; i < basis_pid.size(); ++i)
    trips.push_back(Triplet<double>(pid, new_to_orig[i], basis_pid[i]));
  
  return 0;
}

int nodal_basis_optimizer::optimize(const boost::property_tree::ptree &pt,
                                    SparseMatrix<double> &N) {
  const size_t basis_num  = max(quad_H_)+1;
  const size_t sample_num = nods_h_.size(2);

  vector<Triplet<double>> all_trips;
  #pragma omp parallel for
  for (size_t i = 0; i < basis_num; ++i) {
    vector<Triplet<double>> trips;
    this->optimize_nodal_basis(i, pt, trips);

    #pragma omp critical
    {
      all_trips.insert(all_trips.end(), trips.begin(), trips.end());
    }
  }

  N.resize(basis_num, sample_num);
  N.setFromTriplets(all_trips.begin(), all_trips.end());
  
  return 0;
}

}
