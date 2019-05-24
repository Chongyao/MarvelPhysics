#include <string>
#include <iostream>
#include <chrono>

#include <libigl/include/igl/readOBJ.h>
#include <libigl/include/igl/writeOBJ.h>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/SparseCore>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/filesystem.hpp>

#include "Point_Sys/src/geometry.h"
#include "Point_Sys/src/gen_points.h"
#include "Point_Sys/src/get_nn.h"
#include "Point_Sys/src/points_energy.h"
#include "Point_Sys/src/data_stream.h"
#include "Point_Sys/src/gen_surf.h"
#include "io.h"



using namespace marvel;
using namespace std;
using namespace Eigen;
using namespace igl;
using namespace chrono;
using namespace boost;

int main(int argc, char** argv){

  //cout.precision(20);  
  cout << "[INFO]>>>>>>>>>>>>>>>>>>>READ JSON FILE<<<<<<<<<<<<<<<<<<" << endl;
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
  
  cout << "[INFO]>>>>>>>>>>>>>>>>>>>IMPORT MESH<<<<<<<<<<<<<<<<<<" << endl;
  const string mesh_name = pt.get<string>("surf");
  const string indir = pt.get<string>("indir");
  const string outdir = pt.get<string>("outdir") + mesh_name;
  //mkdir
  boost::filesystem::path outpath(outdir);
  if ( !boost::filesystem::exists(outdir) )
    boost::filesystem::create_directories(outdir);
  MatrixXi surf;
  MatrixXd nods;
  readOBJ((indir+mesh_name+".obj").c_str(), nods, surf);
  cout << "surf: " << surf.rows() << " " << surf.cols() << endl << "nods: " << nods.rows() << " " << nods.cols() << endl;
  
  surf.transposeInPlace();
  nods.transposeInPlace();

    cout << "[INFO]>>>>>>>>>>>>>>>>>>>Generate sampled points<<<<<<<<<<<<<<<<<<" << endl;
  MatrixXd points(3,3);
  MatrixXd test(3, 3);
  gen_points(nods, surf, pt.get<size_t>("num_in_axis"), points);

  // #if 1
  points = nods;
  // #endif
  size_t dim = points.cols();
  cout <<"generate points done." << endl;
  
  cout << "[INFO]>>>>>>>>>>>>>>>>>>>Build spatial hash<<<<<<<<<<<<<<<<<<" << endl;
  spatial_hash SH(points, pt.get<size_t>("nn_num"));

  
  cout << "[INFO]>>>>>>>>>>>>>>>>>>>Build Point System<<<<<<<<<<<<<<<<<<" << endl;
  //calc volume 
  double volume = clo_surf_vol(nods, surf);
  //calc support radii
  VectorXd sup_radi = SH.get_sup_radi();
  
  cout << "[INFO]>>>>>>>>>>>>>>>>>>>sup_radii<<<<<<<<<<<<<<<<<<" << endl;
  cout << sup_radi << endl;
  // assert(0);
  //get friends of every point
  vector<vector<size_t>> friends_all(dim);
#pragma omp parallel for
  for(size_t i = 0; i < dim; ++i){
    SH.get_friends(points.col(i), sup_radi(i), friends_all[i]);
  }
  
  // cout << "[INFO]>>>>>>>>>>>>>>>>>>>friends<<<<<<<<<<<<<<<<<<" << endl;
  // for(size_t i = 0; i < dim; ++i){
  //   cout << "i = " << i <<endl;
  //   for(size_t j = 0; j < friends_all[i].size(); ++j){
  //     cout << friends_all[i][j] << " ";
  //   }
  //   cout << endl;
  // }

  point_sys PS(points, pt.get<double>("rho"), pt.get<double>("Young"), pt.get<double>("Poission"), volume, pt.get<double>("kv"), friends_all, sup_radi);


  cout << "[INFO]>>>>>>>>>>>>>>>>>>>TEST<<<<<<<<<<<<<<<<<<" << endl;
  
  cout << "[INFO]>>>>>>>>>>>>>>>>>>>initial position<<<<<<<<<<<<<<<<<<" << endl;
  cout << points << endl;
  
  Matrix3d change;
  change << 2, 0, 0,
      0, 2, 0,
      0, 0, 2;
  cout << change << endl;

  MatrixXd displace = change * points - points;
  
  cout << "[INFO]>>>>>>>>>>>>>>>>>>>displace<<<<<<<<<<<<<<<<<<" << endl;
  cout << displace << endl;
  std::shared_ptr<dat_str_core<double, 3>>  dat_str = make_shared<energy_dat>(dim);


  PS.pre_compute(dat_str);

  
  PS.Val(displace.data(), dat_str);
  PS.Gra(displace.data(), dat_str);
  cout << "[INFO]>>>>>>>>>>>>>>>>>>>gra<<<<<<<<<<<<<<<<<<" << endl;
  cout << dat_str->get_gra() << endl << endl;

  //check gra by difference

  
  cout << "[INFO]>>>>>>>>>>>>>>>>>>>difference check<<<<<<<<<<<<<<<<<<" << endl;
  { auto init_val = dat_str->get_val();
    auto init_gra = dat_str->get_gra();
    cout << init_val << endl << endl << endl << endl;
    MatrixXd new_gra(3, dim);
    double delt_x = 1e-10;
    for(size_t i = 0; i < points.size(); ++i){
      dat_str->set_zero();
      displace(i) += delt_x;
      PS.Val(displace.data(), dat_str);
      new_gra(i) = dat_str->get_val() - init_val;
      displace(i) -= delt_x;
    }
    new_gra /= delt_x;
    cout << -new_gra << endl <<endl << endl;
    auto delt_gra = init_gra - (- new_gra);
    cout<< delt_gra.array() / init_gra.array() << endl << endl;
    
  }
  
  PS.Val(displace.data(), dat_str);
  PS.Gra(displace.data(), dat_str);
  auto gra_now = dat_str->get_gra();  
  PS.Hes(displace.data(), dat_str);
  
  // dat_str->get_hes().setFromTriplets(dat_str.hes_trips.begin(), dat_str.hes_trips.end());
  
  cout << MatrixXd(dat_str->get_hes()) << endl;
  MatrixXd init_hes = MatrixXd(dat_str->get_hes());

  cout << "[INFO]>>>>>>>>>>>>>>>>>>>numeric difference hes<<<<<<<<<<<<<<<<<<" << endl;
  //check Hessian by difference
  {
    
    double delt_x = 1e-8;
    
    MatrixXd hes_dif(3 * dim, 3 * dim);
    for(size_t i = 0; i < points.size(); ++i){
      dat_str->set_zero();
      displace(i) += delt_x;
      PS.Val(displace.data(), dat_str);
      PS.Gra(displace.data(), dat_str);
      MatrixXd one_hes = (dat_str->get_gra() - gra_now) / delt_x;
      hes_dif.col(i) = Map<VectorXd>(one_hes.data(), 3*dim);
      displace(i) -= delt_x;
    }
    cout << -hes_dif <<endl << endl;
    MatrixXd delt_hes = init_hes- (-hes_dif);
    cout << delt_hes <<endl << endl;
    cout << delt_hes.array() / init_hes.array()<< endl;

  }

  




  cout << "[INFO]>>>>>>>>>>>>>>>>>>>Hessian<<<<<<<<<<<<<<<<<<" << endl;
  cout << MatrixXd(dat_str->get_hes());
  
  
  // displace += dat_str->get_gra();
  

  // dat_str.set_zero();
  // PS.calc_defo_gra(displace.data(), dat_str);
  // PS.Gra(displace.data(), dat_str);
  // cout << dat_str.ela_val_ << endl;
  //numeric difference
  // assert(0);
  
      
      

}
