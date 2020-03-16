#include "DEFINE_TYPE.h"
#include "basic_energy.h"
#include "implicit_euler.h"
#include "io.h"

#include "FEM/src/elas_energy.h"
#include "FEM/src/mass_matrix.h"
#include "Multigrid/src/Sparsification.h"
#include <iostream>


using namespace std;
using namespace Eigen;
using namespace marvel;


// using FLOAT_TYPE = double;

using TET_ELAS = BaseElas<FLOAT_TYPE, 3, 4, 1, 1, linear_csttt, basis_func, quadrature>;
using HEX_ELAS = BaseElas<FLOAT_TYPE, 3, 8, 1, 2, linear_csttt, basis_func, quadrature>;

int main(int argc, char** argv){
  Eigen::initParallel();
  // std::cout.precision(17);
  const char* filename = argv[1];
  
  Matrix<FLOAT_TYPE, -1, -1> nods(1, 1);
  MatrixXi cells(1, 1);

  const string type = argv[4];
  if(type  == "tet")
    mesh_read_from_vtk<FLOAT_TYPE, 4>(filename, nods, cells);
  else if(type == "hex")
    mesh_read_from_vtk<FLOAT_TYPE, 8>(filename, nods, cells);
    
  const size_t num_nods = nods.cols();
  cout <<"V"<< nods.rows() << " " << nods.cols() << endl << "T " << cells.rows() << " "<< cells.cols() << endl;

  const string outdir = argv[3];
  
  //set mtr
  constexpr  FLOAT_TYPE rho = 20;
  constexpr  FLOAT_TYPE Young = 80.0;
  constexpr  FLOAT_TYPE poi = 0.3;
  constexpr  FLOAT_TYPE gravity = 98;
  constexpr  FLOAT_TYPE dt = 0.01;
  const      FLOAT_TYPE w_pos = 1e6;
  const      size_t num_frame = 100;


  //read fixed points
  vector<size_t> cons(0);
  const char* cons_file_path = argv[2];
  read_fixed_verts_from_csv(cons_file_path, cons);
  cout << "constrint " << cons.size() << " points" << endl;
  
  //calc mass vector
  Matrix<FLOAT_TYPE, -1, 1> mass_vec(nods.rows() * num_nods);
  // calc_mass_vector<FLOAT_TYPE>(nods, cells, rho, mass_vec);
  if(type == "tet")
    mass_calculator<FLOAT_TYPE, 3, 4, 1, 1, basis_func, quadrature>(nods, cells, rho, mass_vec);
  else if (type == "hex")
    mass_calculator<FLOAT_TYPE, 3, 8, 1, 2, basis_func, quadrature>(nods, cells, rho, mass_vec);

  cout << "build energy" << endl;
  shared_ptr<Matrix<FLOAT_TYPE, -1, -1>> init_points_ptr  = make_shared<Matrix<FLOAT_TYPE, -1, -1>>(Matrix<FLOAT_TYPE, -1, -1>::Zero(nods.rows(), nods.cols()));
  enum energy_type{ELAS, GRAV, KIN, POS};
  vector<shared_ptr<Functional<FLOAT_TYPE, 3>>> ebf(POS + 1);{
    if(type == "tet" )
      ebf[ELAS] = make_shared<TET_ELAS>(nods, cells, Young, poi);
    else if(type == "hex")
      ebf[ELAS] = make_shared<HEX_ELAS>(nods, cells, Young, poi);

    // ebf[GRAV] = make_shared<gravity_energy<FLOAT_TYPE, 3>>(num_nods, 1, gravity, mass_vec, 'y');
    // ebf[KIN] = make_shared<momentum<FLOAT_TYPE, 3>>(nods.data(), num_nods, mass_vec, dt);
    ebf[KIN] = nullptr;
    ebf[POS] = make_shared<position_constraint<FLOAT_TYPE, 3>>(nods.data(), num_nods, w_pos, cons);

  }
  cout << "assemble energy" << endl;
  
  shared_ptr<Functional<FLOAT_TYPE, 3>> energy;
  try {
    energy = make_shared<energy_t<FLOAT_TYPE, 3>>(ebf);

  } catch ( std::exception &e ) {
    cerr << e.what() << endl;
    exit(EXIT_FAILURE);
  }

  //Sovle

  const string filename_tmp = outdir  + "/frame_origin.vtk";
  shared_ptr<dat_str_core<FLOAT_TYPE, 3>>  dat_str = make_shared<dat_str_core<FLOAT_TYPE, 3>>(num_nods);
  newton_iter<FLOAT_TYPE, 3> imp_euler(dat_str, energy, dt, 20, 1e-4, true, false, true);

  
  //sparsification
  energy->Hes(nods.data(), dat_str);
  SparseMatrix<double> L = dat_str->get_hes();
  cout << L.topLeftCorner(10, 10) << endl;
  getchar();
  
  
  Adjc_graph graph(L);
  graph.Sparsification();
  
  VectorXd eigvals = MatrixXd(L).eigenvalues().real();
  sort(eigvals.data(), eigvals.data() + eigvals.size());
  cout << "condition number " << eigvals.array().abs().maxCoeff() / eigvals.array().abs().minCoeff() << endl;;
  
  
  
  vector<Triplet<double>> trips;
  graph.build_reordered_mat_from_graph(trips);
  SparseMatrix<double> L_new(L.rows(), L.cols());{
    L_new.reserve(trips.size());
    L_new.setFromTriplets(trips.begin(), trips.end());
  }


  eigvals = MatrixXd(L_new).eigenvalues().real();
  sort(eigvals.data(), eigvals.data() + eigvals.size());
  // cout << eigvals << endl; 
  cout << "condition number " << eigvals.array().abs().maxCoeff() / eigvals.array().abs().minCoeff() << endl;;

  {
    SparseLU<SparseMatrix<double>> lu(L_new);
    MatrixXd test_CN = MatrixXd(L) * lu.solve(MatrixXd::Identity(L.rows(), L.rows()));
    cout << "test CN" << endl;
    eigvals = test_CN.eigenvalues().real();
    sort(eigvals.data(), eigvals.data() + eigvals.size());
    cout << eigvals << endl;
    cout << "condition number " << eigvals.array().abs().maxCoeff() / eigvals.array().abs().minCoeff() << endl;;

  }
  cout <<  "Schur_complement"  << endl;
  SparseMatrix<double> L_H, L_ff;
  Schur_complement(L_new, graph.num_coarse_, L_H, L_ff);
  eigvals = MatrixXd(L_H).eigenvalues().real();
  sort(eigvals.data(), eigvals.data() + eigvals.size());
  cout << eigvals << endl;
  cout << "condition number " << eigvals.array().abs().maxCoeff() / eigvals.array().abs().minCoeff() << endl;;

  return 0;
}

