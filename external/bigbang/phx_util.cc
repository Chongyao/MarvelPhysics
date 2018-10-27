#include "phx_util.h"

#include <zjucad/matrix/matrix.h>
#include <zjucad/matrix/itr_matrix.h>
#include <jtflib/mesh/mesh.h>
#include <jtflib/mesh/util.h>
#include <IpIpoptApplication.hpp>
#include <IpSolveStatistics.hpp>
#include <IpIpoptCalculatedQuantities.hpp>
#include <zjucad/matrix/io.h>
#include "constraint.h"
#include "io.h"
#include "func_transform.h"
#include "ipopt_solver.h"
#include "mass_matrix.h"
#include "util.h"
using namespace std;
using namespace zjucad::matrix;
using namespace Eigen;
using namespace Ipopt;

namespace bigbang {

static void gen_surface_traction(const mati_t &surf, const matd_t &normals,
                                 const int i, const int j, VectorXd &force) {
  assert(i >= 0 && i < 3 && j >= 0 && j < 3);
  itr_matrix<double *> extf(3, force.size()/3, force.data());
  extf = zeros<double>(3, force.size()/3);

  const matd_t Id = eye<double>(3);
  const matd_t P = (Id(colon(), i)*trans(Id(colon(), j))+
                    Id(colon(), j)*trans(Id(colon(), i)))/2;

  for (size_t fi = 0; fi < surf.size(2); ++fi) {
    const matd_t f = P*normals(colon(), fi);
    extf(colon(), surf(colon(), fi)) += f*ones<double>(1, surf.size(1))/surf.size(1);
  }
}

int solve_harmonic_disp_3d(const mati_t &cell, const matd_t &nods,
                           const shared_ptr<Functional<double>> &elas,
                           const double wf, matd_t &hdisp, const matd_t *deformed_nods, const bool modify_stiff) {

  cout << "[INFO] solve harmonic displacement field..." << endl;

  mati_t surf;
  if ( cell.size(1) == 4 ) {
    // tet mesh
    using jtf::mesh::face2tet_adjacent;
    shared_ptr<face2tet_adjacent> f2t(face2tet_adjacent::create(cell));
    jtf::mesh::get_outside_face(*f2t, surf, true, &nods);
  } else if ( cell.size(1) == 8 ) {
    // hex mesh
    using jtf::mesh::face2hex_adjacent;
    shared_ptr<face2hex_adjacent> f2h(face2hex_adjacent::create(cell));
    jtf::mesh::get_outside_face(*f2h, surf);
  }

  matd_t normal;
  jtf::mesh::cal_face_normal(surf, nods, normal, false);

  //-> visualize the normal
  #if 0
  write_surf_vec_field("./tmp_output/A_normals.vtk", surf, nods, normal);
  #endif

  //-> config constraints
  vector<shared_ptr<Constraint<double>>> cbf(2); {
    cbf[0] = make_shared<zero_moment_constraint>(*deformed_nods);
    cbf[1] = make_shared<first_moment_constraint>(*deformed_nods);
  }
  shared_ptr<Constraint<double>> opt_con = make_shared<constraint_t<double>>(cbf);

  //-> config energies
  vector<shared_ptr<Functional<double>>> ebf(2);
  if(!modify_stiff)
    ebf[0] = make_shared<linearization_wrapper<double>>(elas, &nods[0]);
  else{
    ebf[0] = make_shared<linearization_wrapper<double>>(elas, &((*deformed_nods)[0]));
    // ebf[0] = make_shared<linearization_wrapper<double>>(elas, &nods[0]);
  }
  
  
    // vector<Triplet<double>> tmp_trips;
    // ebf[0]->Hes(&((*deformed_nods)[0]), &tmp_trips);
    // SparseMatrix<double> tmp_K(ebf[0]->Nx(), ebf[0]->Nx());
    // tmp_K.setFromTriplets(tmp_trips.begin(), tmp_trips.end());
    // cerr << "hes dim: " << tmp_K.rows() << " " << tmp_K.cols() << endl;

    // MatrixXd dense_K(tmp_K);
    // matd_t dense_KK = itr_matrix<const double *>(dense_K.rows(), dense_K.cols(), dense_K.data());
    // matd_t S(dense_KK.size(1), 1);

    // eig(dense_KK, S);
    // // cerr  << S(colon(0, 10)) << endl;
    // cout << "S " << S<<endl;

    // cerr << "here" << endl;
  

  SparseMatrix<double> K_h, M_h;
  VectorXd ev_h; {
    vector<Triplet<double>> trips;
    ebf[0]->Hes(&((*deformed_nods)[0]), &trips);
    K_h.resize(ebf[0]->Nx(), ebf[0]->Nx());
    K_h.setFromTriplets(trips.begin(), trips.end());
    calc_mass_matrix(cell, nods, 1.0, 3, &M_h, true);
    solve_gen_eig_prob(K_h, M_h, 100, 0, "arpaca", &ev_h, nullptr);
  }

  for(size_t i = 0; i < 10; ++i){
    cout << ev_h[i]<< endl;
  }
 
  // int PD = 1;
  // for(size_t i = 0; i < S.size(); ++i){

  //   if(S[i]<0)
  //     PD *= -1;
  // }
  // cout << "PD is " << PD;

  
  ebf[1] = make_shared<general_fext_energy>(nods.size(), wf);
  shared_ptr<Functional<double>> opt_obj = make_shared<energy_t<double>>(ebf);
  
  //-> allocate harmonic displacement field
  hdisp = zeros<double>(nods.size(), 6);
  matd_t curr;
  VectorXd ext_force;

  size_t cnt = 0;
  for (int i = 0; i < 3; ++i) {
    for (int j = i; j < 3; ++j) {      
      ext_force = VectorXd::Zero(nods.size(), 1);
      gen_surface_traction(surf, normal, i, j, ext_force);
      dynamic_pointer_cast<general_fext_energy>(ebf[1])->update_ext_force(ext_force);

      curr = *deformed_nods;
      SmartPtr<ipopt_opt_framework> optprb = new ipopt_opt_framework(opt_obj, opt_con, &curr[0]);      
      SmartPtr<IpoptApplication> app = IpoptApplicationFactory();
      app->Options()->SetIntegerValue("max_iter", 10);
           app->Options()->SetStringValue("linear_solver", "ma86");
      app->Initialize();
      app->OptimizeTNLP(optprb);
      
      hdisp(colon(), cnt) = curr(colon())-nods(colon());
      ++cnt;
    }
  }

  cout << "[INFO] ...solve done!" << endl;
  return 0;
}

}
