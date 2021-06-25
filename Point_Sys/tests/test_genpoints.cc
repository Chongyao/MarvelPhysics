#include <string>
#include <iostream>

#include <Eigen/Core>
#include "Point_Sys/src/gen_points.h"
#include <libigl/include/igl/readOBJ.h>
#include <libigl/include/igl/writeOBJ.h>
#include "Point_Sys/src/get_nn.h"

#include "io.h"
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/filesystem.hpp>
#include <chrono>

using namespace marvel;
using namespace std;
using namespace Eigen;
using namespace igl;
using namespace chrono;

int main(int argc, char** argv){
  boost::property_tree::ptree pt;{
    // const string jsonfile_path = argv[1];
    
    // cout << jsonfile_path << endl;
    // const size_t ext = jsonfile_path.rfind(".json");
    // if (ext != std::string::npos){
    //   read_json(jsonfile_path, pt);
    //   cout << "read json successful" <<endl;
    // }
    // else{
    //   cout << "json file extension error" << endl;
    //   return 0;
    // }
  }
  
  cout << "[INFO]>>>>>>>>>>>>>>>>>>>IMPORT MESH<<<<<<<<<<<<<<<<<<" << endl;
  const string mesh_name = pt.get<string>("surf","dragon");
  const string indir = pt.get<string>("indir", "../../Point_Sys/data/");
  const string outdir = pt.get<string>("outdir", "../../Point_Sys/data/") + mesh_name;
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
  cout << surf.rows() << " " << surf.cols() << endl;
  cout << nods.rows() << " " << nods.cols() << endl;                                               
  
  cout << "[INFO]>>>>>>>>>>>>>>>>>>>Generate sampled points<<<<<<<<<<<<<<<<<<" << endl;
  MatrixXd points;
  // const size_t num_in_axis = pt.get<size_t>("num_in_axis", 25);
  // cout << "[INFO]: num in axis: " << num_in_axis << endl;
  // gen_points(nods, surf, num_in_axis, points, false);

  const double interval = 0.01;
  gen_points(nods, surf, interval, points, false);
  auto point_filename = outdir + "/" + mesh_name + "_points_" +to_string(points.cols());
  cout << "[INFO]: output file :" << point_filename << endl;
  point_write_to_vtk((point_filename + ".vtk").c_str(), points.data(), points.cols());
  
  MatrixXi surf_points = MatrixXi::Zero(0,0);
  writeOBJ((point_filename + ".obj").c_str(), points.transpose(), surf_points);

  cout << points.rows() << " " << points.cols() << endl;
                                                                                  



  cout << "spatial hash " <<endl;
  MatrixXi NN_;
  VectorXd sup_radii_;
  auto start = system_clock::now();
  spatial_hash SH(points, 10);
  auto end = system_clock::now();
  auto duration = duration_cast<microseconds>(end - start);
  cout <<  "花费了" 
       << double(duration.count()) * microseconds::period::num / microseconds::period::den 
       << "秒" << endl;

  cout << "all done " << endl;
                                                                                                                        

  
}

