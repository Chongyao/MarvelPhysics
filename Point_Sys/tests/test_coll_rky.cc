#include <iostream>
#include <limits>
#include <Collision/CollisionDetect-cloth/src/Collision_zcy.h>
using namespace std;

// const double DOUBLE_MAX = std::numeric_limits<double>::max(); 
const double DOUBLE_MAX = 100;
int main(int argc, char** argv){


  
  cout << "[INFO]>>>>>>>>>>>>>>>>>>>SET CUBE<<<<<<<<<<<<<<<<<<" << endl;
  vector<unsigned int> cube_surf =
      // {1,3,0,
      //  7,5,4,
      //  4,1,0,
      //  5,2,1,
      //  2,7,3,
      //  0,7,4,
      //  1,2,3,
      //  7,6,5,
      //  4,5,1,
      //  5,6,2,
      //  2,6,7,
      //  0,3,7};
      {
        0,1,2
      };
  
  vector<double> cube_nods =
      {
        0.510261,0.441394,0.52167,
        0.0849153, 0.752542,0.3,
        -0.87304,0.441254,0.3
      };
      // {1,-1,-1,
      //  1,-1,1,
      //  -1,-1,1,
      //  -1,-1,-1,
      //  1,1,-1,
      //  1,1,1,
      //  -1,1,1,
      //  -1,1,-1};
  vector<double> cube_nods_new = {
    0.510297,0.441397,0.520889,
    0.0849433,0.752562,0.300636,
    -0.873071,0.441254,0.29936
  };
  vector<double> cube_nods_pre = cube_nods;
  // const auto num_surf = 12;
  // const auto num_nods = 8;
  const auto num_surf = 1;
  const auto num_nods = 3 ; 
  cout << "[INFO]>>>>>>>>>>>>>>>>>>>SET PLANE<<<<<<<<<<<<<<<<<<" << endl;
  vector<unsigned int> plane_surf =
      {1, 2, 0,
       1, 3, 2};
      // {0, 1, 2};
  const double plane_z = -1.5;
  vector<double> plane_nods ={
    -40.8575,-40.956,0.3,
    50.1425,-40.956,0.3,
    -40.8575,50.044,0.3,
    50.1425,50.044,0.3
  };
      // {0, DOUBLE_MAX, plane_z,
      //  -DOUBLE_MAX, -DOUBLE_MAX, plane_z,
      //  DOUBLE_MAX,-DOUBLE_MAX, plane_z};
  cout << "[INFO]>>>>>>>>>>>>>>>>>>>COLLISION<<<<<<<<<<<<<<<<<<" << endl;
  auto COLL_ptr = Collision_zcy::getInstance();
  COLL_ptr->Transform_Pair(0, 1);
  COLL_ptr->Transform_Mesh(num_nods, num_surf,
                             cube_surf.data(), cube_nods.data(), cube_nods_pre.data(), 0);
  // COLL_ptr->Transform_Mesh(3, 1, plane_surf.data(), plane_nods.data(), plane_nods.data(), 1);
  COLL_ptr->Transform_Mesh(4, 2, plane_surf.data(), plane_nods.data(), plane_nods.data(), 1);  

  COLL_ptr->Collid();
  COLL_ptr->Transform_Mesh(num_nods, num_surf,
                             cube_surf.data(), cube_nods_new.data(), cube_nods_pre.data(), 0);
  COLL_ptr->Collid();

#if 1
    auto pairs = COLL_ptr->getContactPairs();
    auto times = COLL_ptr->getContactTimes();
    for(size_t i = 0; i < pairs.size(); ++i){
      // cout << "ContactTime is "  << times[i] << endl;
      cout << pairs[i].size() << endl;
      cout << "pair " << i << " :" << endl;
      uint mesh_id, face_id;
      pairs[i][0].get(mesh_id, face_id);
      cout << "mesh id is " << mesh_id << " face id is " << face_id << endl;
      pairs[i][1].get(mesh_id, face_id);
      cout << "mesh id is " << mesh_id << " face id is " << face_id << endl;      
    }    
  
#endif

#if 0

  for(size_t i = 0; i < 999999; ++i){
    cout <<endl<<endl<< "time is " << i << endl << " before: " << endl;
    for(auto& p : cube_nods){cout << p << " ";}
    auto iter = i;
    cout << iter << endl;
    for(size_t j = 0; j < num_nods; ++j){
      cube_nods[j*3 + 2] -= 0.0003;
    }
    cout <<endl<< "after : " << endl;
    for(auto& p : cube_nods){cout << p << " ";}
    cout << endl;

    
    COLL_ptr->Collid();

    auto pairs = COLL_ptr->getContactPairs();
    auto times = COLL_ptr->getContactTimes();
    cout <<" times size is " <<  times.size() << endl;
    assert(pairs.size() == 0);
    for(size_t i = 0; i < pairs.size(); ++i){
      // cout << "ContactTime is "  << times[i] << endl;
      cout << pairs[i].size() << endl;
      cout << "pair " << i << " :" << endl;
      uint mesh_id, face_id;
      pairs[i][0].get(mesh_id, face_id);
      cout << "mesh id is " << mesh_id << " face id is " << face_id << endl;
      pairs[i][1].get(mesh_id, face_id);
      cout << "mesh id is " << mesh_id << " face id is " << face_id << endl;      
    }


    COLL_ptr->Transform_Mesh(num_nods, num_surf,
                             cube_surf.data(), cube_nods.data(), cube_nods_pre.data(), 0, false);
    
    cube_nods_pre = cube_nods;
  }

  #endif
  return 0;
}

