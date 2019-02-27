#include <string>
#include <iostream>

#include <Eigen/Core>
#include <libigl/include/igl/readOBJ.h>

#include "io.h"
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/filesystem.hpp>
#include <chrono>
#include <limits>
#include <Collision/CollisionDetect-rigid/src/Collision_eigen.h>
using namespace marvel;
using namespace std;
using namespace Eigen;
using namespace igl;
using namespace chrono;

// const double DOUBLE_MAX = std::numeric_limits<double>::max();
const double DOUBLE_MAX = 100;
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
  cout << surf.rows() << " " << surf.cols() << endl;
  cout << nods.rows() << " " << nods.cols() << endl;

  const auto num_surf = static_cast<size_t>(surf.cols());
  const auto num_nods = static_cast<size_t>(nods.cols());  
  cout << "[INFO]>>>>>>>>>>>>>>>>>>>COLL<<<<<<<<<<<<<<<<<<" << endl;
  
  cout << surf << endl;

  //set plane
  MatrixXi plane_surf(3, 1);
  plane_surf << 0 , 1 , 2;
  cout << plane_surf << endl;
  MatrixXd plane_nods(3, 3);
  plane_nods << 0, -DOUBLE_MAX, DOUBLE_MAX,
      DOUBLE_MAX, -DOUBLE_MAX, -DOUBLE_MAX,
      -3, -3, -3;
  cout << plane_nods << endl;
  

  auto COLL_ptr = Collision_zcy::getInstance();
  // COLL_ptr->Transform_Pair(0, 1);

  COLL_ptr->Transform_Mesh(3, 1, plane_surf.data(), plane_nods.data(), plane_nods.data(), 0);
  COLL_ptr->Transform_Mesh(num_nods, num_surf,
                           surf.data(), nods.data(), nods.data(), 1);

  auto nods_pre = nods;
  for(size_t i = 0; i < pt.get<size_t>("times"); ++i){
    
    nods.row(2) -= MatrixXd::Ones(1, nods.cols()) * 0.333;
    
    COLL_ptr->Collid();
    auto pairs = COLL_ptr->getContactPairs();
    auto times = COLL_ptr->getContactTimes();
    cout <<" times size is " <<  times.size() << endl;
    for(size_t i = 0; i < pairs.size(); ++i){
      cout << "ContactTime is "  << times[i] << endl;
      cout << pairs[i].size() << endl;
      cout << "pair " << i << " :" << endl;
      uint mesh_id, face_id;
      pairs[i][0].get(mesh_id, face_id);
      cout << "mesh id is " << mesh_id << " face id is " << face_id << endl;
      pairs[i][1].get(mesh_id, face_id);
      cout << "mesh id is " << mesh_id << " face id is " << face_id << endl;      
    }
    COLL_ptr->Transform_Mesh(num_nods, num_surf,
                             surf.data(), nods.data(), nods_pre.data(), 1, false);

    nods_pre = nods;
    cout << nods << endl << endl;
  }

  
  

  
  
  auto start = system_clock::now();
  auto end = system_clock::now();
  auto duration = duration_cast<microseconds>(end - start);
  cout <<  "花费了" 
       << double(duration.count()) * microseconds::period::num / microseconds::period::den 
       << "秒" << endl;

  cout << "all done " << endl;

  
}

