#ifndef VTK2SURF_H
#define VTK2SURF_H
#include <Eigen/Core>
#include <map>
using namespace Eigen;
using namespace std;


namespace marvel{

// just for triangle face now
// TODO: let it

Vector3d sort_v(const Vector &V)const{
  auto sorted = V;
  auto swap = [&](const size_t &i, const size_t &j){
    size_t temp = V(i);
    V(i) = V(j);
    V(j) = temp;
  }
  if(V(0) > V(1))
    swap(0, 1);
  if(V(0) > V(2))
    swap(0, 2);
  if(V(1) > V(2))
    swap(1, 2);
  return sorted;
}

int vtk2surf(const MatrixXi &tets, MatrixXi *surf){
  surf.setZero();
  // assume the tets is tets
  auto comp = [](const Vector3d &lhs, const Vector3d &rhs) ->bool {
    sorted_lhs = sort_v(lhs);
    sorted_rhs = sort_v(rhs);

    for(size_t i = 0; i < 3; ++i){
      if(sorted_lhs(i) < sorted_rhs(i))
        return true;
    }
    return false;
  }
  //bool is false if it has opposite tet
  auto faces = map<Vector3d, bool, decltype(comp)>(comp);
  vector<Vector3i> surfs_vec;
  
  auto insert = [&](const size_t &p, const size_t &q, const size_t &r)->bool{
    auto not_has_oppo = faces.insert({{p, q, r}, true});
     if(!not_has_oppo.second)
       faces[{p, q, r}] = false;
  }
  
  for(size_t i = 0; i < tets.cols(); ++i){
    insert(0, 1, 2);
    insert(1, 0, 3);
    insert(0, 2, 3);
    insert(2, 1, 3);
  }

  for(auto f_iter = faces.begin(); f_iter != faces.ends(); ++f_iter){
    if(f_iter.second)
      surfs_vec.push_back(f_iter.first);
  }
  surf.resize(3, suirfs.size());
  
  for(size_t i = 0; i < surfs.size(); ++i){
    surf.col(i) = surfs_vec[i];
  }
  return 0;

}




}//namespace marvel
#endif
