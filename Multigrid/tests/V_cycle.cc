#include "multigrid.h"

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/filesystem.hpp>

#include "DEFINE_TYPE.h"
#include "basic_energy.h"
#include "implicit_euler.h"
#include "io.h"
#include "FEM/src/elas_energy.h"
#include "FEM/src/mass_matrix.h"
#include <iostream>
#include <string>

using namespace std;
using namespace Eigen;
using namespace marvel;

using HEX_ELAS = BaseElas<FLOAT_TYPE, 3, 8, 1, 2, linear_csttt, basis_func, quadrature>;

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

  const size_t num_layers = pt.get<string>("num_layers");
  const string mesh_name = pt.get<string>("mesh_name");
  
  //read the finest mesh
  Matrix<FLOAT_TYPE, -1, -1> nods(1, 1);
  MatrixXi cells(1, 1);
  mesh_read_from_vtk<FLOAT_TYPE, 8>(mesh_name + ".vtk", nods, cells);

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
  const string cons_file_path = mesh_name + ".csv";
  read_fixed_verts_from_csv(cons_file_path.c_str(), cons);
  cout << "constrint " << cons.size() << " points" << endl;

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
  newton_iter<FLOAT_TYPE, 3> imp_euler(dat_str, energy, dt, 20, 1e-4, true, false, true);


  
  //set hierarchy
  VS<transfer> transfers(num_layers - 1);{
    
  }
  VS<layer> layers(num_layers);{
    layers[0] = make_shared<layer>(dat_str->get_hes(), false, 20);
    #pragma omp parallel for
    for(size_t i = 1; i < num_layers; ++i)
      layers[i] = make_shared<layer>(transfers[i - 1]->R_ * layers[i - 1]->A_ * transfers[i - 1]->I_, false, 20);
  }
  
  
  return 0;
}


