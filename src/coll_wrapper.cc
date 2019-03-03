#include "coll_wrapper.h"


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
  COLL_ptr->Transform_Mesh(core_num_nods, core_num_tris, core_tris_ptr->data(), core_nods.data(), core_nods.data(), 0, false);
  for(size_t i = 0; i < obta_num; ++i){
    COLL_ptr->Transform_Mesh(obta_nods[i]->cols(), obta_surfs[i]->cols(), obta_surfs[i]->data(), obta_nods[i]->data(), obta_nods[i]->data(), i + 1, false);    
  }
  
  COLL_ptr->Collid();
}


int coll_wrapper::Collide(const std::vector<std::shared_ptr<Eigen::MatrixXi>>& obta_surfs,
                      const std::vector<std::shared_ptr<Eigen::MatrixXd>>& obta_nods,
                      const MatrixXd& new_core_velo, const MatrixXd& new_core_nods){
  COLL_ptr->Transform_Mesh(core_num_nods, core_num_tris, core_tris_ptr->data(), new_core_nods.data(), pre_core_nods.data(), 0, false);
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
        
    for(size_t j = 0; j < pairs.size(); ++j){
      unsigned int mesh_id1, face_id1, mesh_id2, face_id2;{

        pairs[j][0].get(mesh_id1, face_id1);
        pairs[j][1].get(mesh_id2, face_id2);
        // cout << ">>>>>>>>>>>>>>>>collid<<<<<<<<<<<<<<<"<<endl;

        if(mesh_id1 == mesh_id2)
          continue;
        if(mesh_id2 == 0){
          mesh_id2 = mesh_id1;
          mesh_id1 = 0;

          auto exchange = face_id2;
          face_id2 = face_id1;
          face_id1 = exchange;
        }
      }//mesh_id,face_id...
      cout << mesh_id1 << " " << mesh_id2 << " " << face_id1 << " " << face_id2 << endl;

      for(size_t tri_dim = 0; tri_dim < 3; ++tri_dim){
        coll_info one_info((*core_tris_ptr)(tri_dim, face_id1), mesh_id2, face_id2, times[j]);
        candidates.insert(one_info);
      }
    }//loop for pairs

    for(auto iter = candidates.begin(); iter != candidates.end(); ++iter){
      size_t vert_id = iter->point_id_, obta_id = iter->mesh_id_, coll_plane_id = iter->face_id_;
      cout << "vert_id is " << vert_id << " " << obta_id << " " << coll_plane_id << endl;
      auto plane_nods = get_tri_pos(*(obta_surfs[obta_id - 1]), *(obta_nods[obta_id - 1]), coll_plane_id);

      point_response(plane_nods.data(), iter->time_,
                     pre_core_nods.col(vert_id).data(), new_core_nods.col(vert_id).data(),
                     pre_core_velo.col(vert_id).data(), new_core_velo.col(vert_id).data());
      
    }
#pragma omp parallel for
    for(size_t j = 0; j < core_num_nods; ++j){
      assert(new_core_nods(2, j) >= 0.3);
    }
        

        
  }//IF COLLIDE
  pre_core_nods = new_core_nods;
  pre_core_velo = new_core_velo;

  return 0;
  
}

}//namespace marvel
