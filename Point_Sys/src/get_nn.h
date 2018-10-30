#ifndef GET_NN_H
#define GET_NN_H

#include<vector>
#include<memory>
#include<list>
#include<Eigen/Core>

#include <unordered_map>

namespace std{
template <>
struct hash<Eigen::Vector3i>{
  typedef size_t result_type;
  typedef Eigen::Vector3i argument_type;
  size_t operator()(const Eigen::Vector3i &key) const{
    // return  ( (key(0)*73856093) ^ (key(1)*19349663) ^ (key(2)*83492791) ) % 5999;
    return  ( (key(0)*73856093) ^ (key(1)*19349663) ^ (key(2)*83492791) );          
  }
};
  
}//namespace: std

namespace marvel{
struct pair_dis{
  size_t m;
  size_t n;
  double dis;
};
//hash table only for nearest neighbours
//TODO:: make it a template

// struct hlist_node {
//     struct hlist_node *next, **pprev;
  
// };
//>>>>>>>>>>>>>>>>>>>>>>>>>>OLD VERISON <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
// struct hlist_head{
//   struct hlist_node *first;
// };
// typedef size_t value_type;
// typedef Eigen::Vector3i key_type;
// typedef std::vector<std::list<value_type>> table_type ;

// class spatial_hash{
//  public:
//   spatial_hash(const size_t &table_size_);
//   int insert(const key_type &key, const value_type &value);
//   int get_val(const key_type &key, std::list<value_type> &vals);
//  private:
//   size_t hash_func(const key_type &key);
//   bool find_ele(const std::list<value_type> &bucket, const value_type &val);
//   size_t table_size;
//   size_t cell_size;

//   std::vector<size_t> prime_num;
//   std::shared_ptr<table_type> hash_table;
// };

//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<OLD VERSION>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

class spatial_hash{
 public:
  spatial_hash(const Eigen::MatrixXd &points_, const size_t &nn_num_);
  const Eigen::MatrixXi& get_NN() const;
  const Eigen::VectorXd& get_sup_radi() const;  
 private:

  Eigen::MatrixXd points;
  size_t nn_num;

  std::unordered_multimap<Eigen::Vector3i, size_t> points_hash;
  Eigen::MatrixXi points_dis;
  Eigen::MatrixXi NN;
  Eigen::VectorXd sup_radi;
  // std::vector<pair_dis> NN_cand;
  int find_NN(const size_t &point_id, std::vector<pair_dis> &NN_cand);
  int hash_NNN();


};
//calculate nearest neighbours of N points
int calc_NNN(const Eigen::MatrixXd &points, Eigen::MatrixXi &NN, Eigen::VectorXd &sup_radi, const size_t &nn_num);

int hash_NNN(const Eigen::MatrixXd &points, Eigen::MatrixXi &NN, Eigen::VectorXd &sup_radi, const size_t &nn_num);


}
#endif
