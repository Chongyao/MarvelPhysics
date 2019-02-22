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
#include <Collision/CollisionDetect-cloth/src/Collision.h>
using namespace marvel;
using namespace std;
using namespace Eigen;
using namespace igl;
using namespace chrono;

const double DOUBLE_MAX = std::numeric_limits<double>::max();
// const double DOUBLE_MAX = 9999;
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
  
  vector<unsigned int> vec_plane_surf(plane_surf.size());
  vector<double> vec_plane_nods(plane_nods.size());
  // vector<double> vec_plane_nods_next(plane_nods.size());
  copy(plane_nods.data(), plane_nods.data() + plane_nods.size(), &vec_plane_nods[0]);
  copy(plane_surf.data(), plane_surf.data() + plane_surf.size(), &vec_plane_surf[0]);

  vector<unsigned int> vec_obj_surf(surf.size());
  vector<double> vec_obj_nods(nods.size());
  vector<double> vec_obj_nods_next(nods.size());
  copy(nods.data(), nods.data() + nods.size(), &vec_obj_nods[0]);
  copy(surf.data(), surf.data() + surf.size(), &vec_obj_surf[0]);
  
  //coll
   // Collision COLL;
  auto COLL_ptr = Collision::getInstance();
  // COLL_ptr->Transform_Pair(0, 1);
  COLL_ptr->Transform_Mesh(3, 1, vec_plane_surf, vec_plane_nods, vec_plane_nods, 0);
  COLL_ptr->Transform_Mesh(num_nods, num_surf,
                      vec_obj_surf, vec_obj_nods, vec_obj_nods, 1);

  cout << nods << endl << endl;
  for(size_t i = 0; i < pt.get<size_t>("times"); ++i){
    
    nods.row(2) -= MatrixXd::Ones(1, nods.cols());
    copy(nods.data(), nods.data() + nods.size(), &vec_obj_nods_next[0]);
    
    COLL_ptr->Collid();
    COLL_ptr->Transform_Mesh(num_nods, num_surf,
                        vec_obj_surf, vec_obj_nods_next, vec_obj_nods, 1, false);
    
    vec_obj_nods = vec_obj_nods_next;
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

