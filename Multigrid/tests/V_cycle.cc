#include "Multigrid/src/multigrid.h"

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/filesystem.hpp>
#include "config.h"
#include "DEFINE_TYPE.h"
#include "basic_energy.h"
#include "implicit_euler.h"
#include "io.h"
#include "FEM/src/elas_energy.h"
#include "FEM/src/mass_matrix.h"
#include <iostream>
#include <string>
#include<Eigen/SparseCholesky>
#include <math.h>

using namespace std;
using namespace Eigen;
using namespace marvel;

using HEX_ELAS = BaseElas<FLOAT_TYPE, 3, 8, 1, 2, linear_csttt, basis_func, quadrature>;

int main(int argc, char** argv){
  Eigen::initParallel();
  boost::property_tree::ptree pt;{
    const string jsonfile_path = argv[1];
    cout << jsonfile_path << endl;
    const size_t ext = jsonfile_path.rfind(".json");
    if (ext != std::string::npos){
      read_json(jsonfile_path, pt);
      cout << "read json successful" <<endl;
    }
    else{
      cout << "json file extension error" << endl;
      return 0;
    }
  }

  const size_t num_layers = pt.get<size_t>("num_layers");
  const size_t finest_layer_id = pt.get<size_t>("finest_layer");
  const string mesh_name = pt.get<string>("mesh_name");
  const string path = pt.get<string>("path");
  
  //read the finest mesh
  Matrix<FLOAT_TYPE, -1, -1> nods(1, 1);
  MatrixXi cells(1, 1);
  const string finest_mesh_fle = path+mesh_name + ".sub"+to_string(finest_layer_id)+".vtk";
  mesh_read_from_vtk<FLOAT_TYPE, 8>(finest_mesh_fle.c_str(), nods, cells);
  
  const size_t num_nods = nods.cols();
  cout <<"V"<< nods.rows() << " " << nods.cols() << endl << "T " << cells.rows() << " "<< cells.cols() << endl;

  
  //read the parameters
  const auto rho = pt.get<FLOAT_TYPE>("rho", 20);
  const auto Young = pt.get<FLOAT_TYPE>("Young", 4000);
  const auto poi = pt.get<FLOAT_TYPE>("poi", 0.3);
  const auto gravity = pt.get<FLOAT_TYPE>("gravity", 9.8);
  const auto w_pos = pt.get<FLOAT_TYPE>("w_pos", 1e6);

  //read fixed points
  vector<size_t> cons(0);
  const string cons_file_path = path+mesh_name + ".sub"+to_string(finest_layer_id)+"-fixed.csv";
  read_fixed_verts_from_csv(cons_file_path.c_str(), cons);
  cout << "constrint " << cons.size() << " points" << endl;

  //calc mass vector
  Matrix<FLOAT_TYPE, -1, 1> mass_vec(nods.rows() * num_nods);
  mass_calculator<FLOAT_TYPE, 3, 8, 1, 2, basis_func, quadrature>(nods, cells, rho, mass_vec);


  //build energy
  enum energy_type{ELAS, GRAV, POS};
  vector<shared_ptr<Functional<FLOAT_TYPE, 3>>> ebf(POS + 1);{
    ebf[ELAS] = make_shared<HEX_ELAS>(nods, cells, Young, poi);
    ebf[GRAV] = make_shared<gravity_energy<FLOAT_TYPE, 3>>(num_nods, 1, gravity, mass_vec, 'y');
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

  shared_ptr<dat_str_core<FLOAT_TYPE, 3>>  dat_str = make_shared<dat_str_core<FLOAT_TYPE, 3>>(num_nods);
  newton_iter<FLOAT_TYPE, 3> imp_euler(dat_str, energy, 0.01, 20, 1e-4, true, false, true);

  //set hierarchy
  VS<transfer> transfers(num_layers - 1);{
    shared_ptr<MatrixXd> nods_h_ptr = make_shared<MatrixXd>(nods);
    shared_ptr<MatrixXd> nods_H_ptr = make_shared<MatrixXd>(MatrixXd::Zero(1, 1));
    shared_ptr<MatrixXi> cells_h_ptr = make_shared<MatrixXi>(cells);
    shared_ptr<MatrixXi> cells_H_ptr = make_shared<MatrixXi>(MatrixXi::Zero(1, 1));

    
    for(size_t i = 0; i < num_layers - 1; ++i){
      const string coarse_mesh_file = path + mesh_name + ".sub" + to_string(finest_layer_id - i - 1) + ".vtk";
      cout << coarse_mesh_file << endl;
      mesh_read_from_vtk<FLOAT_TYPE, 8>(coarse_mesh_file.c_str(), *nods_H_ptr, *cells_H_ptr);
      transfers[i] = make_shared<transfer>(get_transfer(*nods_H_ptr, *cells_H_ptr, *nods_h_ptr, *cells_h_ptr));
      {
        Map<const VectorXd> fine_nods(nods_h_ptr->data(), nods_h_ptr->size());
        Map<const VectorXd> coar_nods(nods_H_ptr->data(), nods_H_ptr->size());
        assert((transfers[i]->I_ * coar_nods  - fine_nods).norm() < 1e-6);
      }
      swap(nods_h_ptr, nods_H_ptr);
      swap(cells_h_ptr, cells_H_ptr);

    }
  }

  cout << "================set transfers done================" << endl;
  const string sol_type_str = pt.get<string>("solver_type", "GS");
  solver_type sol_type;
  if(sol_type_str == "GS")
    sol_type = solver_type::GS;
  else if(sol_type_str == "WJ")
    sol_type = solver_type::WJ;
          
  const size_t gs_itrs = pt.get<size_t>("gs_itrs", 40);
  VS<layer> layers(num_layers);{
    dat_str->set_zero();
    energy->Gra(nods.data(), dat_str);
    energy->Hes(nods.data(), dat_str);
    layers[0] = make_shared<layer>(dat_str->get_hes(), sol_type, gs_itrs);
    layers[0]->rhs_ = - dat_str->get_gra();
    for(size_t i = 1; i < num_layers; ++i){
      layers[i] = make_shared<layer>(transfers[i - 1]->R_ * layers[i - 1]->A_ * transfers[i - 1]->I_, sol_type, gs_itrs);
    }
    
  }

  cout << "=================compare to CG=================="<<endl;
  ConjugateGradient<SparseMatrix<FLOAT_TYPE>, Lower|Upper> cg;
  const SparseMatrix<double> K = dat_str->get_hes();
  const VectorXd rhs = -dat_str->get_gra();
  cg.setMaxIterations(pt.get<size_t>("cg_itrs", 2 * rhs.size()));
  __TIME_BEGIN__;
  cg.compute(K);
  VectorXd solu_cg = cg.solve(rhs);
  __TIME_END__("cg ");
  cout << "residual is " << (rhs - K * solu_cg).norm() / rhs.norm()<< endl;

  cout << "================set layers done================" << endl;
  vector<int> one_V(2 * (num_layers - 1), -1);
  for(size_t i = 0; i < num_layers - 1; ++i){
    one_V[i] = 1;
  }

  const size_t num_V = pt.get<size_t>("num_V", 1);
  // vector<int> process(num_V * one_V.size());
  // #pragma omp parallel for
  // for(size_t i = 0; i < num_V; ++i){
  //   copy(one_V.begin(), one_V.end(), process.begin() + i * one_V.size());
  // }
  
  // cout << "================set process done================" << endl;
  // multigrid_process MP(process, layers, transfers);
  // VectorXd solution = VectorXd::Zero(nods.cols());
  // __TIME_BEGIN__;
  // MP.execute(solution.data());
  // __TIME_END__("V cycle ");

  multigrid_process MP(one_V, layers, transfers);
  VectorXd solution = VectorXd::Zero(nods.size());
  VectorXd last_solution = solution;
  vector<double> diff;
  vector<double> error;
  for(size_t i = 0; i < num_V; ++i){
    __TIME_BEGIN__;
    MP.execute(solution.data());
    __TIME_END__(to_string(i) + " V ");
    diff.push_back((solution - last_solution).norm());
    error.push_back((solution - solu_cg).norm());
    last_solution = solution;
  }

  cout << "diff " << endl;
  for(auto& v : diff)
    cout << v << endl;

  // calculate rate of convergence
  vector<double> order;
  for(size_t i = 0; i < diff.size() - 2; ++i){
    order.push_back(log(diff[i + 2] / diff[i + 1]) / log(diff[i + 1] / diff[i]));
  }
  cout << "order " << endl;
  for(auto& v : order)
    cout << v << endl;


  //try linear
  cout << "try linear convergence." << endl;
  vector<double> linear_rate;
  for(size_t i = 0; i < error.size() - 1; ++i){
    linear_rate.push_back(error[i + 1] / error[i]);
  }
  for(auto& v : linear_rate)
    cout << v << endl;

  {//write error
    const string prefix = "../results/V_cycle/layers_" + to_string(num_layers) + "_gs_itrs_" + to_string(gs_itrs);
    boost::filesystem::path outpath(prefix);
    if ( !boost::filesystem::exists(outpath) )
      boost::filesystem::create_directories(outpath);
    
    auto writefile =
        [&](const vector<double>& vec, const string& name){
          ofstream ofs(prefix + "/" + name +".txt");
          for(auto& v : vec){
            ofs << v << "\n";
          }
          ofs.close();
        };
    writefile(error, "error");
    writefile(diff, "diff");
    writefile(order, "order");
    writefile(linear_rate, "linear_rate");

  }



  // cout << "=========compare to llt===============" << endl;
  // SimplicialLDLT<SparseMatrix<double>> llt;
  // __TIME_BEGIN__;
  // llt.compute(K);
  // VectorXd solu_llt = llt.solve(rhs);
  // __TIME_END__("LLT ");
  // cout << "residual is " << (rhs - K * solu_llt).norm() << endl;
  
  return 0;
}


