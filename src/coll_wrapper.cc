#include "coll_wrapper.h"
#include <map>


using namespace std;
using namespace Eigen;
namespace marvel{
Matrix3d get_tri_pos(const MatrixXi& tris, const MatrixXd& verts, const size_t& face_id){
  Matrix3d tri;
  for(size_t i = 0; i < 3; ++i){
    size_t vert_id = tris(i, face_id);
    tri.col(i) = verts.col(vert_id);
  }
  return std::move(tri);
}
coll_wrapper::coll_wrapper(const std::vector<std::shared_ptr<Eigen::MatrixXi>>&obta_surfs,
               const std::vector<std::shared_ptr<Eigen::MatrixXd>>& obta_nods,
               const std::shared_ptr<Eigen::MatrixXi>& core_surf_ptr,
                           const double* core_nods, const size_t& num_nods):
    core_num_nods(num_nods), core_num_tris(core_surf_ptr->cols()),
    core_tris_ptr(core_surf_ptr), COLL_ptr(Collision_zcy::getInstance()){
    Map<const MatrixXd> core_nods_map(core_nods, 3, num_nods);
    pre_core_nods = core_nods_map;
    pre_core_velo = MatrixXd::Zero(3, core_num_nods);
    
    size_t obta_num = obta_surfs.size();

    //make mesh pair
    for(size_t i = 0; i < obta_num; ++i){
      COLL_ptr->Transform_Pair(0, i + 1);
    }

    //pass mesh info
    COLL_ptr->Transform_Mesh(core_num_nods, core_num_tris, (unsigned*)core_tris_ptr->data(), core_nods, core_nods, 0, false);
    for(size_t i = 0; i < obta_num; ++i){
      COLL_ptr->Transform_Mesh(obta_nods[i]->cols(), obta_surfs[i]->cols(), (unsigned*)obta_surfs[i]->data(), obta_nods[i]->data(), obta_nods[i]->data(), i + 1, false);    
    }
  
    COLL_ptr->Collid();    
  }
coll_wrapper::coll_wrapper(const std::vector<std::shared_ptr<Eigen::MatrixXi>>& obta_surfs,
                           const std::vector<std::shared_ptr<Eigen::MatrixXd>>& obta_nods,
                           const std::shared_ptr<Eigen::MatrixXi>& core_surf_ptr,
                           const Eigen::MatrixXd& core_nods):
    core_num_nods(core_nods.cols()), core_num_tris(core_surf_ptr->cols()),
    core_tris_ptr(core_surf_ptr), COLL_ptr(Collision_zcy::getInstance()){

  
  pre_core_nods = core_nods;
  pre_core_velo = MatrixXd::Zero(3, core_num_nods);

  

  size_t obta_num = obta_surfs.size();

  //make mesh pair
  for(size_t i = 0; i < obta_num; ++i){
    COLL_ptr->Transform_Pair(0, i + 1);
  }

  //pass mesh info
  COLL_ptr->Transform_Mesh(core_num_nods, core_num_tris, (unsigned*)core_tris_ptr->data(), core_nods.data(), core_nods.data(), 0, false);
  for(size_t i = 0; i < obta_num; ++i){
    COLL_ptr->Transform_Mesh(obta_nods[i]->cols(), obta_surfs[i]->cols(), (unsigned*)obta_surfs[i]->data(), obta_nods[i]->data(), obta_nods[i]->data(), i + 1, false);    
  }
  
  COLL_ptr->Collid();
}

#if 0
int coll_wrapper::Collide(const std::vector<std::shared_ptr<Eigen::MatrixXi>>& obta_surfs,
                          const std::vector<std::shared_ptr<Eigen::MatrixXd>>& obta_nods,
                          double* new_core_velo,  double* new_core_nods){
  Map<MatrixXd> new_core_velo_map(new_core_velo, 3, core_num_nods);
  Map<MatrixXd> new_core_nods_map(new_core_nods, 3, core_num_nods);
  
  COLL_ptr->Transform_Mesh(core_num_nods, core_num_tris, (unsigned *)core_tris_ptr->data(), new_core_nods, pre_core_nods.data(), 0, false);
  COLL_ptr->Collid();
  
  auto pairs = COLL_ptr->getContactPairs();
  auto times = COLL_ptr->getContactTimes();


  if(pairs.size() != 0){
    // map<size_t , pair<size_t, size_t>> candidates;
    // map<size_t, double> get_time;
    auto coll_comp = [](const coll_info& one, const coll_info& other)->bool{
      if(one.point_id_ != other.point_id_){
        return one.point_id_ < other.point_id_;                        
      }
      else if(one.mesh_id_ != other.mesh_id_){
        return one.mesh_id_ < other.mesh_id_;
      }
      else {
        return(one.face_id_ < other.face_id_);             
      }


    };
    auto candidates = set<coll_info, decltype(coll_comp)>(coll_comp) ;
    map<size_t, bool> get_if_response;
    
    for(size_t j = 0; j < pairs.size(); ++j){
      unsigned int mesh_id1, face_id1, mesh_id2, face_id2;{

        pairs[j][0].get(mesh_id1, face_id1);
        pairs[j][1].get(mesh_id2, face_id2);
        // cout << ">>>>>>>>>>>>>>>>collid<<<<<<<<<<<<<<<"<<endl;

        if(mesh_id1 == mesh_id2)
          continue;
        else
          cout << mesh_id1 << " " << mesh_id2 << " " << face_id1 << " " << face_id2 << endl;
        if(mesh_id2 == 0){
          mesh_id2 = mesh_id1;
          mesh_id1 = 0;

          auto exchange = face_id2;
          face_id2 = face_id1;
          face_id1 = exchange;
        }
      }//mesh_id,face_id...


      for(size_t tri_dim = 0; tri_dim < 3; ++tri_dim){
        size_t point_id  = (*core_tris_ptr)(tri_dim, face_id1);
        coll_info one_info(point_id, mesh_id2, face_id2, times[j]);
        
        candidates.insert(one_info);
        get_if_response.insert({point_id, false});
      }
    }//loop for pairs

    
    for(auto iter = candidates.begin(); iter != candidates.end(); ++iter){
      size_t vert_id = iter->point_id_, obta_id = iter->mesh_id_, coll_plane_id = iter->face_id_;
      cout << "vert_id is " << vert_id << " " << obta_id << " " << coll_plane_id << endl;
      auto plane_nods = get_tri_pos(*(obta_surfs[obta_id - 1]), *(obta_nods[obta_id - 1]), coll_plane_id);

      if(!get_if_response[vert_id] && point_response(plane_nods.data(), iter->time_,
                        pre_core_nods.col(vert_id).data(), new_core_nods_map.col(vert_id).data(),
                        pre_core_velo.col(vert_id).data(), new_core_velo_map.col(vert_id).data())){

        get_if_response[vert_id] = true;
      }
      
      
    }
        
  }//IF COLLIDE
  pre_core_nods = new_core_nods_map;
  pre_core_velo = new_core_velo_map;

  return 0;

}
#endif

int coll_wrapper::Collide(const std::vector<std::shared_ptr<Eigen::MatrixXi>>& obta_surfs,
                      const std::vector<std::shared_ptr<Eigen::MatrixXd>>& obta_nods,
                       MatrixXd& new_core_velo,  MatrixXd& new_core_nods){
  COLL_ptr->Transform_Mesh(core_num_nods, core_num_tris, (unsigned *)core_tris_ptr->data(), new_core_nods.data(), pre_core_nods.data(), 0, false);
  COLL_ptr->Collid();
  
  auto pairs = COLL_ptr->getContactPairs();
  auto times = COLL_ptr->getContactTimes();


  if(pairs.size() != 0){
    // map<size_t , pair<size_t, size_t>> candidates;
    // map<size_t, double> get_time;
    auto coll_comp = [](const coll_info& one, const coll_info& other)->bool{
      if(one.point_id_ != other.point_id_){
        return one.point_id_ < other.point_id_;                        
      }
      else if(one.mesh_id_ != other.mesh_id_){
        return one.mesh_id_ < other.mesh_id_;
      }
      else {
        return(one.face_id_ < other.face_id_);             
      }


    };
    auto candidates = set<coll_info, decltype(coll_comp)>(coll_comp) ;
    map<size_t, bool> get_if_response;
    
    for(size_t j = 0; j < pairs.size(); ++j){
      unsigned int mesh_id1, face_id1, mesh_id2, face_id2;{

        pairs[j][0].get(mesh_id1, face_id1);
        pairs[j][1].get(mesh_id2, face_id2);
        // cout << ">>>>>>>>>>>>>>>>collid<<<<<<<<<<<<<<<"<<endl;

        if(mesh_id1 == mesh_id2)
          continue;
        else
          cout << mesh_id1 << " " << mesh_id2 << " " << face_id1 << " " << face_id2 << endl;
        if(mesh_id2 == 0){
          mesh_id2 = mesh_id1;
          mesh_id1 = 0;

          auto exchange = face_id2;
          face_id2 = face_id1;
          face_id1 = exchange;
        }
      }//mesh_id,face_id...


      for(size_t tri_dim = 0; tri_dim < 3; ++tri_dim){
        size_t point_id  = (*core_tris_ptr)(tri_dim, face_id1);
        coll_info one_info(point_id, mesh_id2, face_id2, times[j]);
        
        candidates.insert(one_info);
        get_if_response.insert({point_id, false});
      }
    }//loop for pairs

    
    for(auto iter = candidates.begin(); iter != candidates.end(); ++iter){
      size_t vert_id = iter->point_id_, obta_id = iter->mesh_id_, coll_plane_id = iter->face_id_;
      cout << "vert_id is " << vert_id << " " << obta_id << " " << coll_plane_id << endl;
      auto plane_nods = get_tri_pos(*(obta_surfs[obta_id - 1]), *(obta_nods[obta_id - 1]), coll_plane_id);

      if(!get_if_response[vert_id] && point_response(plane_nods.data(), iter->time_,
                        pre_core_nods.col(vert_id).data(), new_core_nods.col(vert_id).data(),
                        pre_core_velo.col(vert_id).data(), new_core_velo.col(vert_id).data())){

        get_if_response[vert_id] = true;
      }
      
      
    }
        
  }//IF COLLIDE



  pre_core_nods = new_core_nods;
  pre_core_velo = new_core_velo;
  
  return 0;
  
}

int coll_wrapper::Collide_imp(const std::vector<std::shared_ptr<Eigen::MatrixXi>>& obta_surfs,
                              const std::vector<std::shared_ptr<Eigen::MatrixXd>>& obta_nods,
                              Eigen::MatrixXd& new_core_nods, Eigen::MatrixXd& new_core_velo,
                              energy_dat& dat_str,
                              const VectorXd& mass){
  COLL_ptr->Transform_Mesh(core_num_nods, core_num_tris, (unsigned *)core_tris_ptr->data(), new_core_nods.data(), pre_core_nods.data(), 0, false);
  COLL_ptr->Collid();
  
  auto pairs = COLL_ptr->getContactPairs();
  auto times = COLL_ptr->getContactTimes();


  if(pairs.size() != 0){
    // map<size_t , pair<size_t, size_t>> candidates;
    // map<size_t, double> get_time;
    auto coll_comp = [](const coll_info& one, const coll_info& other)->bool{
      if(one.point_id_ != other.point_id_){
        return one.point_id_ < other.point_id_;                        
      }
      else if(one.mesh_id_ != other.mesh_id_){
        return one.mesh_id_ < other.mesh_id_;
      }
      else {
        return(one.face_id_ < other.face_id_);             
      }


    };
    auto candidates = set<coll_info, decltype(coll_comp)>(coll_comp) ;
    map<size_t, bool> get_if_response;
    
    for(size_t j = 0; j < pairs.size(); ++j){
      unsigned int mesh_id1, face_id1, mesh_id2, face_id2;{

        pairs[j][0].get(mesh_id1, face_id1);
        pairs[j][1].get(mesh_id2, face_id2);
        // cout << ">>>>>>>>>>>>>>>>collid<<<<<<<<<<<<<<<"<<endl;

        if(mesh_id1 == mesh_id2)
          continue;
        else
          cout << mesh_id1 << " " << mesh_id2 << " " << face_id1 << " " << face_id2 << endl;
        if(mesh_id2 == 0){
          mesh_id2 = mesh_id1;
          mesh_id1 = 0;

          auto exchange = face_id2;
          face_id2 = face_id1;
          face_id1 = exchange;
        }
      }//mesh_id,face_id...


      for(size_t tri_dim = 0; tri_dim < 3; ++tri_dim){
        size_t point_id  = (*core_tris_ptr)(tri_dim, face_id1);
        coll_info one_info(point_id, mesh_id2, face_id2, times[j]);
        
        candidates.insert(one_info);
        get_if_response.insert({point_id, false});
      }
    }//loop for pairs

    
    for(auto iter = candidates.begin(); iter != candidates.end(); ++iter){
      size_t vert_id = iter->point_id_, obta_id = iter->mesh_id_, coll_plane_id = iter->face_id_;
      cout << "vert_id is " << vert_id << " " << obta_id << " " << coll_plane_id << endl;
      auto plane_nods = get_tri_pos(*(obta_surfs[obta_id - 1]), *(obta_nods[obta_id - 1]), coll_plane_id);

      Vector3d force;
      if(!get_if_response[vert_id] &&
         imp_response2(plane_nods.data(), iter->time_,
                       pre_core_nods.col(vert_id).data(), new_core_nods.col(vert_id).data(),
                       pre_core_velo.col(vert_id).data(), new_core_velo.col(vert_id).data(),
                      force.data())){
        force *= mass(vert_id);
        get_if_response[vert_id] = true;
        // dat_str.save_ele_gra(vert_id, force);
      }
      
      
    }
        
  }//IF COLLIDE


  pre_core_nods = new_core_nods;
  pre_core_velo = new_core_velo;
  
  return 0;
  

}




}//namespace marvel
