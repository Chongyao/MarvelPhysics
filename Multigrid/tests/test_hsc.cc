#include "DEFINE_TYPE.h"
#define EIGEN_USE_BLAS
#include "basic_energy.h"
#include "implicit_euler.h"
#include "io.h"

#include "FEM/src/elas_energy.h"
#include "FEM/src/poisson.h"
#include <iostream>
#include <functional>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

#include "Multigrid/src/hsc.h"
#include "Multigrid/src/pcg.h"

using namespace std;
using namespace std::placeholders;
using namespace Eigen;
using namespace marvel;
// using boost::property_tree::ptree;

int main(int argc, char** argv){
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
  
  Eigen::initParallel();

  auto common = pt.get_child("common");
  auto phy_paras = pt.get_child("physics");
  
  const char* filename = common.get<string>("mesh").c_str();
  Matrix<FLOAT_TYPE, -1, -1> nods(1, 1);
  MatrixXi cells(1, 1);

  const string type = common.get<string>("type");
  if(type == "hex")
    mesh_read_from_vtk<FLOAT_TYPE, 8>(filename, nods, cells);
    
  const size_t num_nods = nods.cols();
  cout <<"V"<< nods.rows() << " " << nods.cols() << endl << "T " << cells.rows() << " "<< cells.cols() << endl;

  Matrix<FLOAT_TYPE, 1, -1> u = Matrix<FLOAT_TYPE, 1, -1>::Zero(1, num_nods);
  const string outdir = common.get<string>("outdir");
  
  //set mtr
  const FLOAT_TYPE k = phy_paras.get<FLOAT_TYPE>("k", 10);
  const FLOAT_TYPE w_pos = phy_paras.get<FLOAT_TYPE>("w_pos", 1e6);

  Matrix<FLOAT_TYPE, 1, -1> mtr = Matrix<FLOAT_TYPE, 1, -1>::Ones(cells.cols()) * k;

  //read fixed points
  vector<size_t> cons(0);
  Matrix<FLOAT_TYPE, -1, 1>  u_cons = Matrix<FLOAT_TYPE, 1, -1>::Zero(num_nods);{
    vector<size_t> cons_top(0), cons_bottom(0);
    const char* cons_file_path = common.get<string>("top_fixed").c_str();
    read_fixed_verts_from_csv(cons_file_path, cons_top);
    cons_file_path = common.get<string>("bottom_fixed").c_str();
    read_fixed_verts_from_csv(cons_file_path, cons_bottom);
    cons.resize(cons_bottom.size() + cons_top.size());
    copy(cons_top.begin(), cons_top.end(), cons.begin());
    copy(cons_bottom.begin(), cons_bottom.end(), cons.begin() + cons_top.size());
    
    #pragma omp parallel for
    for(size_t i = 0; i < cons_top.size(); ++i)
      u_cons(cons_top[i]) = 10;
    #pragma omp parallel for
    for(size_t i = 0; i < cons_bottom.size(); ++i)
      u_cons(cons_bottom[i]) =  5 * sin(nods(1, cons_bottom[i])) + 5 * cos(nods(2, cons_bottom[i]));
    cout << "constraint " << cons.size() << " points" << endl;
  }
  
  
  cout << "build energy" << endl;
  enum energy_type{POISSON, BD};
  vector<shared_ptr<Functional<FLOAT_TYPE, 1>>> ebf(BD + 1);{
    ebf[POISSON] = make_shared<HEX_POISSON>(nods, cells, mtr);
    ebf[BD] = make_shared<position_constraint<FLOAT_TYPE, 1>>(u_cons.data(), num_nods, w_pos, cons);

  }
  cout << "assemble energy" << endl;
  
  shared_ptr<Functional<FLOAT_TYPE, 1>> energy;
  try {
    energy = make_shared<energy_t<FLOAT_TYPE, 1>>(ebf);

  } catch ( std::exception &e ) {
    cerr << e.what() << endl;
    exit(EXIT_FAILURE);
  }


  
  //Sovle
  const string filename_tmp = outdir  + "/frame_origin.vtk";
  shared_ptr<dat_str_core<FLOAT_TYPE, 1>>  dat_str = make_shared<dat_str_core<FLOAT_TYPE, 1>>(num_nods);
  newton_iter_myPCG<FLOAT_TYPE, 1> imp_euler(dat_str, energy, 10000, 20, 1e-10, true, false, true);

  energy->Hes(u.data(), dat_str);
  #if 0
  {//test M matrix
    const SPM M = dat_str->get_hes();
    M.pruned(1e-15);
    #pragma omp parallel for
    for (int k=0; k<M.outerSize(); ++k)
      for (decltype(M)::InnerIterator it(M,k); it; ++it){
        bool test =(it.index() == k && it.value() > 0) || ( it.index() != k && it.value() < 0);
        if(!test)
          cout << it.row() << " " << it.col() << " " << it.value() << endl;
      }

  }
  #endif

  if(pt.get<bool>("PD", true)){
    HSC hsc = set_hierarchy(dat_str->get_hes(), pt);
    imp_euler.set_preconditioner(std::bind(&HSC::solve, hsc, std::__1::placeholders::_1));

  }

  for(size_t f_id = 0; f_id < 1; ++f_id){
    cout << "[frame " << f_id << "]" << endl;
    imp_euler.solve(u.data());
    
    const string filename = outdir  + "/frame_" + to_string(f_id) + ".vtk";
    mesh_write_to_vtk<FLOAT_TYPE, 8>(filename.c_str(), nods, cells);
    point_scalar_append2vtk(false, filename.c_str(), u.transpose(), u.size(), "u");
  }
  return 0;
}




