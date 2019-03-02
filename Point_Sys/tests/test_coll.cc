#include <string>
#include <iostream>

#include <Eigen/Core>
#include <libigl/include/igl/readOBJ.h>
#include <libigl/include/igl/writeOBJ.h>

#include "io.h"
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/filesystem.hpp>
#include <chrono>
#include <limits>
#include <Collision/CollisionDetect-rigid/src/Collision_eigen.h>
#include "coll_response.h"
using namespace marvel;
using namespace std;
using namespace Eigen;
using namespace igl;
using namespace chrono;

// const double DOUBLE_MAX = std::numeric_limits<double>::max();
Matrix3d get_tri_pos(const MatrixXi& tris, const MatrixXd& verts, const size_t& face_id){
  Matrix3d tri;
  for(size_t i = 0; i < 3; ++i){
    size_t vert_id = tris(i, face_id);
    tri.col(i) = verts.col(vert_id);
  }
  return std::move(tri);
}
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
  cout << surf << endl;
  const auto num_surf = static_cast<size_t>(surf.cols());
  const auto num_nods = static_cast<size_t>(nods.cols());  
  cout << "[INFO]>>>>>>>>>>>>>>>>>>>COLL<<<<<<<<<<<<<<<<<<" << endl;
  

  //set plane
  MatrixXi plane_surf(3, 1);
  plane_surf << 0 , 1 , 2;
  cout << plane_surf << endl;
  const double plane_z = 0.5;
  MatrixXd plane_nods(3, 3);
  plane_nods << 0, -DOUBLE_MAX, DOUBLE_MAX,
      DOUBLE_MAX, -DOUBLE_MAX, -DOUBLE_MAX,
      plane_z, plane_z, plane_z;
  cout << plane_nods << endl;
  

  auto COLL_ptr = Collision_zcy::getInstance();
  COLL_ptr->Transform_Pair(0, 1);
  

  COLL_ptr->Transform_Mesh(3, 1, plane_surf.data(), plane_nods.data(), plane_nods.data(), 1);
  COLL_ptr->Transform_Mesh(num_nods, num_surf,
                             surf.data(), nods.data(), nods.data(), 0);


  COLL_ptr->Collid();
  double gravity = pt.get<double>("gravity");
  MatrixXd velo = MatrixXd::Zero(3, num_nods),
      new_velo = velo;
  MatrixXd new_nods = nods;
  MatrixXd acce = MatrixXd::Zero(3, num_nods);
  acce.row(2) = MatrixXd::Ones(1, num_nods) * (-gravity);
  double delt_t = pt.get<double>("time_step");
  size_t max_iter = pt.get<size_t>("times");
  cout << "max iter is " << max_iter <<" " << num_surf << " " << num_nods << endl;
  for(size_t i = 0; i < max_iter; ++i){
    cout << "iter is " << i << endl;
    // cout << "position" << endl;
    cout << new_nods << endl;
    // cout << "velocity" << endl;
    // cout << new_velo << endl;
    new_velo += acce * delt_t;
    new_nods += new_velo * delt_t;

    cout << "position" << endl;
    cout << new_nods << endl;
    cout << "velocity" << endl;
    cout << new_velo << endl;

    COLL_ptr->Transform_Mesh(num_nods, num_surf,
                             surf.data(), new_nods.data(), nods.data(), 0);
    
    
    COLL_ptr->Collid();
    auto pairs = COLL_ptr->getContactPairs();
    auto times = COLL_ptr->getContactTimes();
    cout <<" times size is " <<  times.size() << endl;
    // assert(pairs.size() == 0);
    vector<size_t> if_response(num_nods, false);
    map<size_t , pair<size_t, size_t>> candidates;
    map<size_t, double> get_time;
    
    for(size_t j = 0; j < pairs.size(); ++j){
      unsigned int mesh_id1, face_id1, mesh_id2, face_id2;{
        cout <<endl<<endl<<endl<< "j is " << j <<endl;
        pairs[j][0].get(mesh_id1, face_id1);
        pairs[j][1].get(mesh_id2, face_id2);
        cout << mesh_id1 << " " << mesh_id2 << " " << face_id1 << " " << face_id2 <<endl;
        if(mesh_id1 == mesh_id2)
          continue;
        if(mesh_id2 == 0){
          auto exchange = mesh_id2;
          mesh_id2 = mesh_id1;
          mesh_id1 = exchange;

          exchange = face_id2;
          face_id2 = face_id1;
          face_id1 = exchange;
        }
      }//mesh_id,face_id...

      for(size_t tri_dim = 0; tri_dim < 3; ++tri_dim){
        candidates.insert({surf(tri_dim, face_id1), {mesh_id2, face_id2}});
        get_time.insert({surf(tri_dim), times[j]});
      }
    }//loop for pairs

    for(auto iter = candidates.begin(); iter != candidates.end(); ++iter){
      size_t vert_id =  iter->first,
          obta_id = iter->second.first, coll_plane_id = iter->second.second;
      cout << "vert id is " <<  vert_id << endl;
      point_response(plane_nods.data(), get_time[vert_id],
                     nods.col(vert_id).data(), new_nods.col(vert_id).data(),
                     velo.col(vert_id).data(), new_velo.col(vert_id).data(),0, 0.5);
      
    }
    for(size_t i = 0; i < num_nods; ++i){
      if(new_nods(2, i) < plane_z)
        cout << "i is " << i << endl << new_nods.col(i) << endl;
      assert(new_nods(2, i) >= plane_z);
    }

    
    if(i%50 == 0){
      auto surf_filename = outdir  + "/" + mesh_name + "_" + to_string(i) + ".obj";
      writeOBJ(surf_filename.c_str(), new_nods.transpose(), surf.transpose());      
    }
    nods = new_nods;
    velo = new_velo;
  }

  


  
}

