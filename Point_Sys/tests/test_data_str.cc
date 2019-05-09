#define EIGEN_USE_BLAS
#include "data_str.h"

#include <chrono>
#include <random>
#include <memory>
#include <iostream>
#include <time.h>
using namespace std;
using namespace marvel;
using namespace Eigen;
using namespace chrono;



int main(int argc, char** argv){
  Eigen::initParallel();
  srand (time(NULL));
  const size_t dof = 200000 * 3;
  
  auto mydata = make_shared<dat_str_core<double, 3>>(dof);
  
  //>>>>>>>>>>>test val
  double sum_val = 0;
  size_t min = 0, max = 4;
  vector<double> vals(dof);
  
  #pragma omp parallel for
  for(size_t i = 0; i < dof; ++i){
    double val_tmp = min + (rand() % static_cast<int>(max - min + 1));
    vals[i] = val_tmp;
  }
  
  auto start = system_clock::now();
  #pragma omp parallel for
  for(size_t i = 0; i < dof; ++i){
    #pragma omp atomic
    sum_val += vals[i];
  }
  auto end = system_clock::now();
  auto duration = duration_cast<microseconds>(end - start);
  cout <<  "dir sum 花费了" 
       << double(duration.count()) * microseconds::period::num / microseconds::period::den 
       << "秒" << endl;
  cout << "sum dir is " << sum_val << endl;

  
  start = system_clock::now();
  #pragma omp parallel for
  for(size_t i = 0; i < dof; ++i){
    mydata->save_val(vals[i]);
  }
  end = system_clock::now();
  duration = duration_cast<microseconds>(end - start);
  cout <<  "dat sum 花费了" 
       << double(duration.count()) * microseconds::period::num / microseconds::period::den 
       << "秒" << endl;
  cout << "sum dat is " << mydata->get_val() << endl;


  //test>>>>>>>val




  cout << "[INFO]>>>>>>>>>>>>>>>>>>>test gra<<<<<<<<<<<<<<<<<<" << endl;
  mydata->set_zero();
  VectorXd gra_add = VectorXd::Random(dof);
  cout << "ground truth " << gra_add.sum() << endl;
  start = system_clock::now();
  mydata->save_gra(gra_add);
  end = system_clock::now();
  duration = duration_cast<microseconds>(end - start);
  cout <<  "save gra all  花费了" 
       << double(duration.count()) * microseconds::period::num / microseconds::period::den 
       << "秒" << endl;
  cout << "gra sum " << (mydata->get_gra()).sum() << endl;


  mydata->set_zero();
  start = system_clock::now();
  #pragma omp parallel for
  for(size_t i = 0; i < dof / 3; ++i){
    mydata->save_gra(i, gra_add.segment(i * 3, 3));
  }

  end = system_clock::now();
  duration = duration_cast<microseconds>(end - start);
  cout <<  "save gra parallel  花费了" 
       << double(duration.count()) * microseconds::period::num / microseconds::period::den 
       << "秒" << endl;
  cout << "gra sum " << (mydata->get_gra()).sum() << endl;
  
  
  cout << "[INFO]>>>>>>>>>>>>>>>>>>>hes<<<<<<<<<<<<<<<<<<" << endl;
  mydata->set_zero();
  // mydata->hes_reserve(VectorXi::Constant(dof, 4));
  MatrixXi row_ids(6, dof);
  MatrixXd nz_vals = MatrixXd::Random(6, dof);
  cout << "ground truth is " << nz_vals.sum()<< endl;
  #pragma omp parallel for
  for(size_t i = 0; i < dof; ++i){
    for(size_t j = 0; j < 6; ++j){
      row_ids(j, i) = min + (rand() % static_cast<int>(max - min + 1));
    }
  }
  //use class dat_str_core
  start = system_clock::now();
  // #pragma omp parallel for
  for(size_t i = 0; i < dof; ++i){
    for(size_t j = 0; j < 6; ++j){
      mydata->save_hes(row_ids(j, i), i, nz_vals(j, i));
    }
  }

  end = system_clock::now();
  duration = duration_cast<microseconds>(end - start);
  cout <<  "save hes by class花费了" 
       << double(duration.count()) * microseconds::period::num / microseconds::period::den 
       << "秒" << endl;
  cout << "hes sum in class is " << (mydata->get_hes()).sum() << endl;

  //use triplets
  vector<Triplet<double>> trips;
  SparseMatrix<double> HES(dof, dof);
  trips.reserve(6 * dof);
  start = system_clock::now();
  for(size_t i = 0; i < dof; ++i){
    for(size_t j = 0; j < 6; ++j){
      trips.push_back(Triplet<double>(row_ids(j,i), i, nz_vals(j, i)));
    }
  }
  HES.setFromTriplets(trips.begin(), trips.end());
  end = system_clock::now();
  duration = duration_cast<microseconds>(end - start);
  cout <<  "save hes by class花费了" 
       << double(duration.count()) * microseconds::period::num / microseconds::period::den 
       << "秒" << endl;
  cout << "hes sum in class is " << HES.sum() << endl;
  
  return 0;
}
