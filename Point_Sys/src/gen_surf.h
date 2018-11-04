#ifndef GEN_SURF_H
#define GEN_SURF_H

#include "get_nn.h"
#include <Eigen/Core>
namespace marvel{


template<typename T>
class deform_surf{
 public:
  deform_surf(const Eigen::MatrixXi &surf, const Eigen::Matrix<T,-1,-1> &nods, const Eigen::Matrix<T, -1, -1> &sam_points);
  virtual ~deform_surf(){}
  // Eigen::Matrix<T, -1, -1>& update_surf(const Eigen::Matrix<T, -1, -1>& moved_points) const = 0;

 protected:
  const Eigen::MatrixXi surf_;
  const Eigen::Matrix<T, -1, -1> nods_;
  const Eigen::Matrix<T, -1, -1> sam_points_;
};


//TODO::complete this update
/********************************CLASS deform_surf_MLS******************************/
template<typename T>
class deform_surf_MLS : public deform_surf<T>{
 public:
  deform_surf_MLS(const Eigen::MatrixXi &surf, const Eigen::Matrix<T,-1,-1> &nods, const Eigen::Matrix<T, -1, -1> &sam_points, const std::vector<std::vector<size_t>> &friends, const T &kernel_cof);
  Eigen::Matrix<T, -1, -1> update_surf(const Eigen::Matrix<T, -1 ,-1>& moved_points_dis, const Eigen::Matrix<T, -1, -1> &def_gra_all) const;
 private:
  std::vector<std::vector<size_t>> friends_;
  T kernel_cof_;
  T kernel(const T &r, const T &h) const;
  T kernel(const size_t &vertex_id, const size_t &point_id, const T &h) const;
};


/********************************CLASS deform_surf_MLS******************************/





/********************************CLASS deform_surf_LI******************************/

//linear interpolate vertex by four non-co-plane sampled points
template<typename T>
class deform_surf_LI : public deform_surf<T>{
 public:
  //TODO: do not use copy construct!!
  deform_surf_LI(const Eigen::MatrixXi &surf, const Eigen::Matrix<T,-1,-1> &nods, const Eigen::Matrix<T, -1, -1> &sam_points, const Eigen::MatrixXi &four_NN);
  Eigen::Matrix<T, -1, -1>& update_surf(const Eigen::Matrix<T, -1, -1> &moved_points) const;
 private:
  const Eigen::MatrixXi surf_;
  const Eigen::Matrix<T, -1, -1> nods_;
  const Eigen::Matrix<T, -1, -1> sam_points_;
  const Eigen::MatrixXi four_NN_;

  Eigen::Matrix<T, -1, -1> basis_weig;
  
  
  
};
//class deform_surf_LI
/********************************CLASS deform_surf_LI******************************/


}//namespace: marvel
#include "gen_surf.imp"
#endif
