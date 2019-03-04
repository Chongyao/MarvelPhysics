#ifndef COLL_RESPONSE
#define COLL_RESPONSE
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <iostream>
using namespace Eigen;
using namespace std;
//only for triangles

enum axis{x, y, z};

inline bool project_dim(const Eigen::Vector3d& normal, axis& dim1, axis& dim2){
  dim1 = x;
  dim2 = y;
  
  size_t count = 0;
  for(size_t i = 0; i < 3; ++i){
    if(fabs(normal(i)) < 1e-6){
      if(count == 0){
        dim1 = static_cast<axis>(i);
        dim2 = static_cast<axis>((i + 1) % 3);
        ++count;
      }
      else if(count == 1){
        dim2 = static_cast<axis>(i);
        ++count;
      }
    }
  }
  if(count > 2)
    return false;
  else
    return true;

}

bool check_inside(const axis& dim1, const axis& dim2, const Eigen::Vector3d& point, const Eigen::MatrixXd& obstacle, const double& area){
  
  const double s = ( (point(dim1) - obstacle(dim1, 2)) * (obstacle(dim2, 0) - obstacle(dim2, 2))
                 - (point(dim2) - obstacle(dim2, 2)) * (obstacle(dim1, 0) - obstacle(dim1, 2)) ) /
             ( (obstacle(dim1, 1) - obstacle(dim1, 2)) * (obstacle(dim2, 0) - obstacle(dim2, 2))
           -(obstacle(dim2, 1) - obstacle(dim2, 2)) * (obstacle(dim1, 0) - obstacle(dim1, 2)) );
  const double t = ( (point(dim1) - obstacle(dim1, 2)) * (obstacle(dim2, 1) - obstacle(dim2, 2))
             - (point(dim2) - obstacle(dim2, 2)) * (obstacle(dim1, 1) - obstacle(dim1, 2)) ) /
      ((obstacle(dim1, 0) - obstacle(dim1, 2)) * (obstacle(dim2, 1) - obstacle(dim2, 2))
       - (obstacle(dim2, 0) - obstacle(dim2, 2)) * (obstacle(dim1, 1) - obstacle(dim1, 2)));
  
  cout<< " check through : " << s << " " << t << " " << 1 - s - t << endl;
  if(s != s){
    cout << "nan" << endl << dim1 << dim2 << endl<< point << endl << endl << obstacle;
    cout << " jerer" << endl;
    assert(s == s);
  }

  // cout << obstacle << endl;
  return ( s > 0 && s < 1 && t > 0 && t < 1 && (1 - s - t) > 0 && (1 - s - t) < 1 );
        
  
}
    
Vector3d get_cross_point(const Eigen::Vector3d& before, const Eigen::Vector3d& after,
                         const Eigen::Vector3d& norm, const double& d, bool& is_through, double offset_ = 0){
  const Vector3d line = after - before;
  
  //may divede zero but OK;
  const double offset = ( -d - norm(0) * before(0) - norm(1) * before(1) - norm(2) * before(2) )
      / ( norm(0) * line(0) + norm(1) * line(1) + norm(2) * line(2) );

  // cout << offset << endl;
  if(offset <= 1 && offset >= 0){
    is_through = true;
    offset_ = offset * 0.999;
    return (before + offset * line);
  }
  else{
    is_through = false;
    // cout << "cross point is " << before + offset * line << endl;
    return after;
  }
      
}
#if 0
int response(const double* const obstacle,
             const double time, 
             const double* const pre_pos, const double* const next_pos,
             const double* const pre_velo, const double* const next_velo,
             double* const res_pos, double* const res_velo,
             const double& friction = 0, const double& res = 0.5){
  Map<const MatrixXd> obstacle_(obstacle, 3, 3);
  Map<const MatrixXd> pre_pos_(pre_pos, 3, 3);
  Map<const MatrixXd> next_pos_(next_pos, 3, 3);
  Map<const MatrixXd> next_velo_(next_velo, 3, 3);
  Map<const MatrixXd> pre_velo_(pre_velo,3, 3);
  Map<MatrixXd> res_pos_(res_pos, 3, 3);
  res_pos_ = next_pos_;
  Map<MatrixXd> res_velo_(res_velo, 3, 3);
  res_velo_ = next_velo_;



  // cout << "before response" << endl << "pre pos : " <<endl << pre_pos_ << endl << "after_pos :" <<endl<< next_pos_ <<endl<< "pre velo :" << endl << pre_velo_ << endl << "after velo : " << endl << next_velo_ << endl;
  // cout << "beore  response" <<endl<< res_pos_ << endl << endl << res_velo_ << endl;
  
  double area;
  Vector3d plane_normal;{
    Vector3d one_edge = obstacle_.col(1) - obstacle_.col(0);
    Vector3d other_edge = obstacle_.col(2) - obstacle_.col(0);
    plane_normal = one_edge.cross(other_edge);
    double norm  = plane_normal.norm();
    area = norm / 2;
    if(norm > 1e-6)
      plane_normal = plane_normal / norm;
    else
      plane_normal = Vector3d::Zero();
  }
  
  const double para_d = -plane_normal.dot(obstacle_.col(0)) ;

  for(size_t i = 0; i < 3; ++i){
    bool is_through= false;
    auto cross_point = get_cross_point(pre_pos_.col(i), next_pos_.col(i), plane_normal, para_d, is_through);
    if(is_through && check_inside(cross_point, obstacle_, area)){
      res_pos_.col(i) = cross_point;
      res_velo_.col(i) = pre_velo_.col(i) * (1 - time) + time * next_velo_.col(i);
      Vector3d projection = res_velo_.col(i).dot(plane_normal) * plane_normal;
      res_velo_.col(i) = (res_velo_.col(i) - projection) * (1 - friction) - projection * res;
    }
  }
  // cout << "after  response" <<endl<< res_pos_ << endl << endl << res_velo_ << endl;
  return 0;
}
#endif
bool point_response(const double* const obstacle,
             const double time, 
             const double* const pre_pos, double* const next_pos,
             const double* const pre_velo, double* const next_velo,
             // double* const res_pos, double* const res_velo,
             const double& friction = 0, const double& res = 1){
  
  Map<const MatrixXd> obstacle_(obstacle, 3, 3);
  Map<const Vector3d> pre_pos_(pre_pos);
  Map<Vector3d> next_pos_(next_pos);
  Map<Vector3d> next_velo_(next_velo);
  Map<const Vector3d> pre_velo_(pre_velo);
  
  Vector3d res_pos_ = next_pos_;
  Vector3d res_velo_ = next_velo_;







  
  double area;
  //TODO:can pre_store to speed up
  Vector3d plane_normal;{
    Vector3d one_edge = obstacle_.col(1) - obstacle_.col(0);
    Vector3d other_edge = obstacle_.col(2) - obstacle_.col(0);
    plane_normal = one_edge.cross(other_edge);
    double norm  = plane_normal.norm();
    area = norm / 2;
    if(norm > 1e-6)
      plane_normal = plane_normal / norm;
    else
      plane_normal = Vector3d::Zero();
  }
  
  const double para_d = -plane_normal.dot(obstacle_.col(0)) ;
  

  bool is_through= false;
  double offset_ = 0;
  auto cross_point = get_cross_point(pre_pos_, next_pos_, plane_normal, para_d, is_through, offset_);
  if(!is_through)
    return false;

  axis dim1, dim2;
  if(project_dim(plane_normal, dim1, dim2) && check_inside(dim1, dim2, cross_point, obstacle_, area)){
    cout << "before response" << endl << "pre pos : " <<endl << pre_pos_ << endl << "after_pos :" <<endl<< next_pos_ <<endl;
    cout << "pre velo :" << endl << pre_velo_ << endl << "after velo : " << endl << next_velo_ << endl;
    // res_pos_ = cross_point ;
    res_pos_ = pre_pos_ + offset_ * (next_pos_ - pre_pos_);
    // res_velo_ = pre_velo_ * (1 - time) + time * next_velo_;
    res_velo_ = next_velo_;
    Vector3d projection = res_velo_.dot(plane_normal) * plane_normal;
    res_velo_ = (res_velo_ - projection) * (1 - friction) - projection * res;

    next_pos_ = res_pos_;
    next_velo_ = res_velo_;
    cout << "after  response" <<endl<< res_pos_ << endl << res_velo_ << endl; ;
    return true;        
  }
  else
    return false;



}


#endif
