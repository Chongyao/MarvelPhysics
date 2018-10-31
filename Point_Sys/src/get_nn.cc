#include "get_nn.h"
#include <cmath>
#include <algorithm>
#include <iostream>
#include <unordered_map>
using namespace std;
using namespace Eigen;

namespace marvel{



spatial_hash::spatial_hash(const MatrixXd &points_, const size_t &nn_num_):points(points_), points_num(points_.cols()), nn_num(nn_num_){
  hash_NNN();
}


int spatial_hash::get_shell(const Eigen::Vector3i &query, const int &radi, std::vector<Vector3i> &shell) const {
  assert(radi > -1);
  //init
  shell.clear();

  auto loop = [&](const int &face, const int &face_axis, const int &radi_1, const int &radi_2){
    if(face <= max_id(face_axis) && face >= min_id(face_axis) ){
      int axis_1 = (face_axis + 1)%3, axis_2 = (face_axis + 2)%3;
      for(int j = query(axis_1) - radi_1; j < query(axis_1) + radi_1 + 1; ++j){
        for(int k = query(axis_2) - radi_2; k < query(axis_2) + radi_2 + 1; ++k){
          Vector3i one_grid;
          one_grid[face_axis] = face;
          one_grid[axis_1] = j;
          one_grid[axis_2] = k;
          shell.push_back(one_grid);
        }
      }    

    }
  };
  if(radi > 0){
    for(size_t i = 0; i < 6; ++i){
      loop(i%2 == 0?query(i/2) - radi:query(i/2) + radi, i/2, i<4?radi : radi - 1, i<2?radi : radi - 1); 
    }
  }
  else{
    shell.push_back(query);
  }

  return 0;
}

int spatial_hash::find_NN(const size_t &point_id, vector<pair_dis> &NN_cand){

  size_t cand_num = 0;
  int grid_delt = 0;
  // bool once_more = false;
  int once_more = 0;
  //count for grid_delt;
  do{
    vector<Vector3i> shell;
    get_shell(points_dis.col(point_id), grid_delt, shell);
    for(auto &grid : shell){
      auto range = points_hash.equal_range(grid);
      if( range.first != range.second){
        for_each(range.first, range.second, [&](unordered_multimap<Vector3i,size_t>::value_type  &one_point){
            NN_cand.push_back({point_id, one_point.second, (points.col(point_id) - points.col(one_point.second)).norm()});
          });
      }
    }
    ++grid_delt;
    if (NN_cand.size() > nn_num + 2)
      // once_more = !once_more;
      once_more++;
  }while(NN_cand.size() < nn_num + 2 || once_more<2);
  return 0;
}


int spatial_hash::hash_NNN(){
  //init
  NN.setZero(nn_num, points_num);
  points_hash.clear();
  //set hash parameter
  double cell_size = pow(points.cols()/nn_num, 1/3);
  size_t table_size = size_t(floor(pow(points.cols(), 0.5)));
  
  //generate discretized 3D position
  points_dis = floor(points.array()/cell_size).cast<int>();
  max_id = {points_dis.row(0).maxCoeff(), points_dis.row(1).maxCoeff(), points_dis.row(2).maxCoeff()};
  min_id = {points_dis.row(0).minCoeff(), points_dis.row(1).minCoeff(), points_dis.row(2).minCoeff()};
  
  

  //build hash_map
  points_hash = unordered_multimap<Vector3i,size_t>(table_size);

  //insert elements
#pragma parallel omp for
  for(size_t i = 0; i < points_num; ++i){
    points_hash.insert({points_dis.col(i), i});
  }
  //calc NN

  return 0;
}

const MatrixXi& spatial_hash::get_NN() const {
  return NN;
}
const VectorXd& spatial_hash::get_sup_radi() {
  assert(points.cols() > nn_num);

  //init data

  sup_radi.setZero(points_num);
  
  #pragma parallel omp for
  for(size_t i = 0; i < points_num; ++i){
    vector<pair_dis> NN_cand;
    find_NN(i, NN_cand);
    sort(NN_cand.begin(), NN_cand.end(), [](const pair_dis &a, const pair_dis &b){return a.dis < b.dis;});
    
    for(size_t j = 0; j < nn_num; ++j){
      NN(j, i) = NN_cand[j + 1].n;
      sup_radi[i] += NN_cand[j + 1].dis;
    }
  }
  sup_radi*= 3.0/nn_num;
  cout << NN.block(0,0,nn_num, 10) << endl;
  return sup_radi;
}
int spatial_hash::get_friends(const size_t &point_id, const double &sup_radi, vector<size_t> &friends) const{
  friends.clear();
  int grid_delt = 0;
  bool has_friends = false;
  bool has_points = false;
  int once_more = 0;
  do{
    has_friends = false;
    has_points = false;
    vector<Vector3i> shell;
    get_shell(points_dis.col(point_id), grid_delt, shell);


    for(auto &one_grid : shell){
      auto range = points_hash.equal_range(one_grid);

      if( range.first != range.second){
        has_points = true;
        for_each(range.first, range.second, [&](const unordered_multimap<Vector3i,size_t>::value_type  &one_point){
            double dis = (points.col(one_point.second) - points.col(point_id)).norm();
            if(dis < sup_radi && dis != 0){
              friends.push_back(one_point.second);
              has_friends = true;
            }
          });
      }
    }
    ++grid_delt;

    if(has_points && !has_friends)
      ++once_more;
  }while(has_friends || once_more < 2);
  return 0;
}
int spatial_hash::update_points(const Eigen::MatrixXd &points_){
  points = points_;
  hash_NNN();
  return 0;
}

}//namespace:marvel
