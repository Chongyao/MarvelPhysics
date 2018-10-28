#include "get_nn.h"
#include <cmath>
#include <algorithm>
using namespace std;
using namespace Eigen;

namespace marvel{

spatial_hash::spatial_hash(const size_t &table_size_):table_size(table_size_){
  prime_num = vector<size_t> {7373856093, 19349663, 83492791};
  hash_table = make_shared<table_type>(table_size_);
}

bool spatial_hash::find_ele(const list<value_type> &bucket, const value_type &val){
  for(auto &v : bucket){
    if (v == val )
      return true;
  }
  return false;  
}

size_t spatial_hash::hash_func(const key_type &key){
  return  ( (key(0)*prime_num[0]) ^ (key(1)*prime_num[1]) ^ (key(2)*prime_num[2]) ) % table_size;  
}

int spatial_hash::insert(const key_type &key, const value_type &value){
  size_t table_id = hash_func(key);
  if( !find_ele((*hash_table)[table_id], value) )
    (*hash_table)[table_id].push_front(value);
  return 0;
}


int spatial_hash::get_val(const key_type &key, std::list<value_type> &vals){
  size_t table_id = hash_func(key);
  vals = (*hash_table)[table_id];
}
//TODO: use hash method or exist library
// int calc_NNN(const MatrixXd &points, MatrixXi &NN, VectorXd &sup_radi, const size_t &nn_num){
//   //init
//   NN = MatrixXi(10, points.cols());
//   sup_radi = VectorXd(points.cols());

//   size_t cell_size = size_t(floor(pow(points.cols()/nn_num, 1/3)));
//   size_t table_size = size_t(floor(pow(points.cols(), 0.5)));

//   MatrixXi points_dis(points.rows(), points.cols()); 
//   points_dis = points.cast<int>();
  

//   //make hash table
//   spatial_hash points_hash(table_size);
//   for(size_t i = 0; i < points.cols(); ++i){
//     points_hash.insert(points_dis.col(i), i);
//   }

//   //
//   for(size_t i = 0; i < points.cols(); ++i){
//   }
// }

struct pair_dis{
  size_t m;
  size_t n;
  double dis;
};

int calc_NNN(const MatrixXd &points, MatrixXi &NN, VectorXd &sup_radi){
  NN = MatrixXi(10, points.cols());
  sup_radi = VectorXd(points.cols());
  
  size_t points_num = points.cols();
  vector<pair_dis> pairs(points_num*(points_num - 1) / 2);
  
  auto id_calc = [=](size_t i, size_t j){return ((2*points_num - i - 1)*i/2) + j - i; };
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
    if(total == 10*points_num) break;
    if(count[pairs[i].m] != -1){
      NN(count[pairs[i].m], pairs[i].m) = pairs[i].n;
      sup_radi(pairs[i].m) += pairs[i].dis; 
      count[pairs[i].m]++;
      total++;
      if(count[pairs[i].m] == 10)
        count[pairs[i].m] = -1;
    }
    if(count[pairs[i].n] != -1){
      NN(count[pairs[i].n], pairs[i].n) = pairs[i].m;
      sup_radi(pairs[i].n) += pairs[i].dis;
      count[pairs[i].n]++;
      total++;
      if(count[pairs[i].n] == 10)
        count[pairs[i].n] = -1;
    }
  }

  sup_radi*= 0.3;
  return 0;
}



}//namespace:marvel
