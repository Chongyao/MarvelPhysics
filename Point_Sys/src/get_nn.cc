#include "get_nn.h"
#include <cmath>
using namespace std;
using namespace Eigen;

namespace marvel{

spatial_hash::spatial_hash(const size_t &table_size_):table_size(table_size_){
  prime_num = vector<size_t> {7373856093, 19349663, 83492791};
  hash_table = make_shared<table_type>(table_size_);
}

bool spatial_hash::find_ele(const forward_list<value_type> &bucket, const value_type &val){
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


int spatial_hash::get_val(const key_type &key, std::forward_list<value_type> &vals){
  size_t table_id = hash_func(key);
  vals = (*hash_table)[table_id];
}

int calc_NNN(const MatrixXd &points, MatrixXi &NN, VectorXd &sup_radi, const size_t &nn_num){
  //init
  NN = MatrixXi(10, points.cols());
  sup_radi = VectorXd(points.cols());

  size_t cell_size = size_t(floor(pow(points.cols()/nn_num, 1/3)));
  size_t table_size = size_t(floor(pow(points.cols(), 0.5)));

  MatrixXi points_dis(points.rows(), points.cols()); 
  points_dis = points.cast<int>();
  

  //make hash table
  spatial_hash points_hash(table_size);
  for(size_t i = 0; i < points.cols(); ++i){
    points_hash.insert(points_dis.col(i), i);
  }

  //
  

}





}
