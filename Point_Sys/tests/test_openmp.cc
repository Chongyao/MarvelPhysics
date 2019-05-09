#include<Eigen/Core>
#include <vector>
#include <iostream>
#include <random>
#include <chrono>
using namespace std;
using namespace Eigen;
using namespace chrono;

class mysum{
 public:
  mysum():sum(0){}
  size_t sum;

  int operator()(const size_t& value){
    #pragma omp atomic
    sum += value;
  }
  
};



int main(int argc, char** argv){
  const size_t vec_size = 10000000;
  vector<size_t> test_vec(vec_size, 0);
  size_t  min = 0, max = 200, sum_cri = 0, sum_atomic = 0, sum_true = 0;



  //set value
#pragma omp parallel for
  for(size_t i = 0; i < test_vec.size(); ++i){
    test_vec[i] = min + (rand() % static_cast<int>(max - min + 1));
  }

  //get groundtruth
  for(size_t i = 0; i < test_vec.size(); ++i){
    sum_true += test_vec[i];
  }
  cout << "sum true is " << sum_true << endl;


  //test critical
  auto start = system_clock::now();
#pragma omp parallel for
  for(size_t i = 0; i < test_vec.size(); ++i){
#pragma omp critical
    {
      sum_cri += test_vec[i];
    }
  }

  auto end = system_clock::now();
  auto duration = duration_cast<microseconds>(end - start);
  cout <<  "critical sum 花费了" 
       << double(duration.count()) * microseconds::period::num / microseconds::period::den 
       << "秒" << endl;
  cout << "sum critical is " << sum_cri << endl;


  //test atomic
  start = system_clock::now();
#pragma omp parallel for
  for(size_t i = 0; i < test_vec.size(); ++i){
#pragma omp atomic
    sum_atomic += test_vec[i];
  }

  end = system_clock::now();
  duration = duration_cast<microseconds>(end - start);
  cout <<  "atomic sum 花费了" 
       << double(duration.count()) * microseconds::period::num / microseconds::period::den 
       << "秒" << endl;
  cout << "sum atomic is " << sum_atomic << endl;

  
  //test func
  mysum sum_func;
  start = system_clock::now();
#pragma omp parallel for
  for(size_t i = 0; i < test_vec.size(); ++i){
    sum_func(test_vec[i]);
  }

  end = system_clock::now();
  duration = duration_cast<microseconds>(end - start);
  cout <<  "func sum 花费了" 
       << double(duration.count()) * microseconds::period::num / microseconds::period::den 
       << "秒" << endl;
  cout << "sum func is " << sum_func.sum << endl;  
  

  //test without block
  size_t sum_without_block;
  start = system_clock::now();
#pragma omp parallel for
  for(size_t i = 0; i < test_vec.size(); ++i){
    sum_without_block += test_vec[i];
  }

  end = system_clock::now();
  duration = duration_cast<microseconds>(end - start);
  cout <<  "sum without block 花费了" 
       << double(duration.count()) * microseconds::period::num / microseconds::period::den 
       << "秒" << endl;
  cout << "sum without block is " << sum_without_block << endl;  

  // for(const auto& v : test_vec){
  //   cout << v << endl;
  // }
  
}
