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

int main(int argc, char** argv){

  cout << "[INFO]>>>>>>>>>>>>>>>>>>>obstacle<<<<<<<<<<<<<<<<<<" << endl;

  MatrixXi plane_surf;
  MatrixXd plane_nods;
  readOBJ("Plane.obj", plane_nods, plane_surf);
  plane_surf.transposeInPlace();
  plane_nods.transposeInPlace();
  cout << plane_surf.size() << endl;
  cout << "[INFO]>>>>>>>>>>>>>>>>>>>COLL<<<<<<<<<<<<<<<<<<" << endl;
  

  auto COLL_ptr = Collision_zcy::getInstance();
  COLL_ptr->Transform_Pair(0, 1);


  

  for(size_t file_id = 0; file_id < 798; ++file_id){
    cout << "file id is " << file_id << endl;  

  cout << "[INFO]>>>>>>>>>>>>>>>>>>>IMPORT MESH<<<<<<<<<<<<<<<<<<" << endl;
  MatrixXi surf;
  MatrixXd nods;
  readOBJ("mesh_test_old_" + to_string(file_id) + ".obj", nods, surf);
  cout << "surf: " << surf.rows() << " " << surf.cols() << endl << "nods: " << nods.rows() << " " << nods.cols() << endl;
  surf.transposeInPlace();
  nods.transposeInPlace();
  const auto num_surf = static_cast<size_t>(surf.cols());
  const auto num_nods = static_cast<size_t>(nods.cols());

  
  cout << "[INFO]>>>>>>>>>>>>>>>>>>>new mesh<<<<<<<<<<<<<<<<<<" << endl;

  MatrixXi surf_new;
  MatrixXd nods_new;
  readOBJ("mesh_test_new_" + to_string(file_id) + ".obj", nods_new, surf_new);
  surf_new.transposeInPlace();
  nods_new.transposeInPlace();
  // cout << nods_new.size() << endl;

  
  COLL_ptr->Transform_Mesh(num_nods, num_surf,
                           (unsigned*)surf.data(), nods_new.data(), nods.data(), 0);
  COLL_ptr->Transform_Mesh(4, 2, (unsigned*)plane_surf.data(), plane_nods.data(), plane_nods.data(), 1);
  
  COLL_ptr->Collid();
  auto pairs = COLL_ptr->getContactPairs();
  auto times = COLL_ptr->getContactTimes();
    
  for(size_t j = 0; j < pairs.size(); ++j){
    unsigned int mesh_id1, face_id1, mesh_id2, face_id2;{
      // cout << "j is " << j <<endl;
      pairs[j][0].get(mesh_id1, face_id1);
      pairs[j][1].get(mesh_id2, face_id2);
      if(mesh_id1 != mesh_id2)
        cout << mesh_id1 << " " << mesh_id2 << " " << face_id1 << " " << face_id2 <<endl;
    }//mesh_id,face_id...
    
  }//loop for pairs

  nods = nods_new;

  
}

}
