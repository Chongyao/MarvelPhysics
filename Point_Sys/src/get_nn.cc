#include "get_nn.h"
#include <cmath>
#include <algorithm>
#include <iostream>
#include <unordered_map>
using namespace std;
using namespace Eigen;

namespace marvel{



//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<NEW VERSION>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
spatial_hash::spatial_hash(const MatrixXd &points_, const size_t &nn_num_):points(points_), nn_num(nn_num_){
  hash_NNN();
}

//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<SIMPLE VVERSION>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
int calc_NNN(const MatrixXd &points, MatrixXi &NN, VectorXd &sup_radi, const size_t &nn_num){
  assert(points.cols() > nn_num);
  
  NN.setZero(nn_num, points.cols());
  sup_radi.setZero(points.cols());
  
  size_t points_num = points.cols();
  vector<pair_dis> pairs(points_num*(points_num - 1) / 2);
  auto id_calc = [=](size_t i, size_t j){return ((2*points_num - i - 1)*i/2) + j - i - 1; };
#pragma omp parallel for
  for(size_t i = 0; i < points_num - 1; ++i){
    for(size_t j = i + 1; j < points_num; ++j){      
      pairs[id_calc(i ,j)] = {i, j, (points.col(i) - points.col(j)).norm()};
    }
  }

  sort(pairs.begin(), pairs.end(), [](const pair_dis &a, const pair_dis &b){return a.dis < b.dis;});

  vector<int> count(points_num);
  size_t total = 0;
  for(size_t i = 0; i < pairs.size(); ++i){
    if(total == nn_num*points_num) break;
    if(count[pairs[i].m] != -1){
      NN(count[pairs[i].m], pairs[i].m) = pairs[i].n;
      sup_radi(pairs[i].m) += pairs[i].dis;
      count[pairs[i].m]++;
      total++;
      if(count[pairs[i].m] == nn_num)
        count[pairs[i].m] = -1;
    }
    if(count[pairs[i].n] != -1){
      NN(count[pairs[i].n], pairs[i].n) = pairs[i].m;
      sup_radi(pairs[i].n) += pairs[i].dis;
      count[pairs[i].n]++;
      total++;
      if(count[pairs[i].n] == nn_num)
        count[pairs[i].n] = -1;
    }
  }

  sup_radi*= 3.0/nn_num;
  cout << "done " << endl;
  cout << NN.block(0,0,nn_num, 10) << endl;
  cout << sup_radi.head(10) << endl;


  return 0;
}
//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<SIMPLE VVERSION>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<NEW VVERSION>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

int spatial_hash::get_shell(const Eigen::Vector3i &query, const int &radi, std::vector<Vector3i> &shell){
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
  assert(points.cols() > nn_num);
  //init data
  const size_t points_num = points.cols();
  NN.setZero(nn_num, points_num);
  sup_radi.setZero(points_num);

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
  for(size_t i = 0; i < points_num; ++i){
    points_hash.insert({points_dis.col(i), i});
  }
  //calc NN
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
  cout << "A hehre" << endl;
  sup_radi*= 3.0/nn_num;
  cout << "done " << endl;
  cout << NN.block(0,0,nn_num, 10) << endl;
  cout << sup_radi.head(10) << endl;

  return 0;
}

const MatrixXi& spatial_hash::get_NN() const {
  return NN;
}
const VectorXd& spatial_hash::get_sup_radi() const{
  return sup_radi;
}
}//namespace:marvel
