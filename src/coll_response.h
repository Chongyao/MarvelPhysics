#ifndef COLL_RESPONSE
#define COLL_RESPONSE
#include <Eigen/Core>
#include <Eigen/Geometry>

using namespace Eigen;
//only for triangles
bool check_inside(const Eigen::Vector3d& point, const Eigen::MatrixXd& obstacle, const double& area){
  double s = (obstacle(1, 0)*obstacle(0, 2) - obstacle(0, 0)*obstacle(1, 2) + (obstacle(1, 2) - obstacle(1, 0))*point(0) + (obstacle(0, 0) - obstacle(0, 2))*point(1)) / (2*area);
  double t = (obstacle(0, 0)*obstacle(1, 1) - obstacle(1, 0)*obstacle(0, 1) + (obstacle(1, 0) - obstacle(1, 1))*point(0) + (obstacle(0, 1) - obstacle(0, 0))*point(1)) / (2*area);
  
  // cout<< " check through : " << s << " " << t << " " << 1 - s - t << endl;
  return ( s > 0 && s < 1 && t > 0 && t < 1 && (1 - s - t) > 0 && (1 - s - t) < 1 );
        
  
}
    
Vector3d get_cross_point(const Eigen::Vector3d& before, const Eigen::Vector3d& after,
                         const Eigen::Vector3d& norm, const double& d, bool& is_through){
  const Vector3d line = after - before;
  const double offset = ( -d - norm(0) * before(0) - norm(1) * before(1) - norm(2) * before(2) )
      / ( norm(0) * line(0) + norm(1) * line(1) + norm(2) * line(2) );

  // cout << offset << endl;
  if(offset > 1 || offset < 0){
    is_through = false;
    return after;
  }
  else{
    is_through = true;
    return (before + offset * line);
  }
      
}
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
int point_response(const double* const obstacle,
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




  // cout << "before response" << endl << "pre pos : " <<endl << pre_pos_ << endl << "after_pos :" <<endl<< next_pos_ <<endl<< "pre velo :" << endl << pre_velo_ << endl << "after velo : " << endl << next_velo_ << endl;
  // cout << "beore  response" <<endl<< res_pos_ << endl << endl << res_velo_ << endl;
  
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
  auto cross_point = get_cross_point(pre_pos_, next_pos_, plane_normal, para_d, is_through);
  if(is_through && check_inside(cross_point, obstacle_, area)){
    res_pos_ = cross_point;
    res_velo_ = pre_velo_ * (1 - time) + time * next_velo_;
    Vector3d projection = res_velo_.dot(plane_normal) * plane_normal;
    res_velo_ = (res_velo_ - projection) * (1 - friction) - projection * res;

    next_pos_ = res_pos_;
    next_velo_ = res_velo_;
        
  }

  // cout << "after  response" <<endl<< res_pos_ << endl << endl << res_velo_ << endl;

  return 0;
}


#endif
