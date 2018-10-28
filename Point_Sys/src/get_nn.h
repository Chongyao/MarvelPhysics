#ifndef GET_NN_H
#define GET_NN_H

#include<vector>
#include<memory>
#include<forward_list>
#include<Eigen/Core>
namespace marvel{

//hash table only for nearest neighbours
//TODO:: make it a template

// struct hlist_node {
//     struct hlist_node *next, **pprev;
  
// };

// struct hlist_head{
//   struct hlist_node *first;
// };
typedef size_t value_type;
typedef Eigen::Vector3i key_type;
typedef std::vector<std::forward_list<value_type>> table_type ;



class spatial_hash{
 public:
  spatial_hash(const size_t &table_size_);
  int insert(const key_type &key, const value_type &value);
  int get_val(const key_type &key, std::forward_list<value_type> &vals);
 private:
  size_t hash_func(const key_type &key);
  bool find_ele(const std::forward_list<value_type> &bucket, const value_type &val);
  size_t table_size;
  size_t cell_size;

  std::vector<size_t> prime_num;
  std::shared_ptr<table_type> hash_table;
 
};


//calculate nearest neighbours of N points
int calc_NNN(const Eigen::MatrixXd &points, Eigen::MatrixXi &NN, Eigen::VectorXd &sup_radi);
  
  
  

}

#endif
