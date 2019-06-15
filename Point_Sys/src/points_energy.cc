#include <math.h>
#include <iostream>
#include <memory>
#include <chrono>


#include "points_energy.h"
#include "get_nn.h"

#include <Eigen/SVD>
#include <Eigen/LU>
#include <Eigen/Geometry>


using namespace std;
using namespace Eigen;
using namespace chrono;



#define PI 3.14159265359


namespace marvel{

Matrix3d safe_inv(const MatrixXd& sys_mat){
  //TODO: use better solver considering the sysmertic
  JacobiSVD<MatrixXd> svd(sys_mat, ComputeThinU | ComputeThinV);

  auto sin_val = svd.singularValues();
  Matrix3d inv_sin_val;
  inv_sin_val.setZero(3, 3);
  #pragma omp parallel for
  for(size_t p = 0; p < 3; ++p)
    inv_sin_val(p, p) = sin_val(p)>0?1/sin_val(p):0;

  return std::move(svd.matrixV() * inv_sin_val * svd.matrixU().transpose());
}

Matrix3d cons_law(const Matrix3d &strain, const double &You, const double &Poi){
  assert(Poi > 0 && Poi < 0.5);
  Matrix3d stress = Matrix3d::Zero();
  
  for(size_t i = 0; i < 3; ++i){
    for(size_t j = 0; j < 3; ++j){
      stress(i, i) += j==i ? (1-Poi) * strain(j, j) : Poi*strain(j, j);
    }
    stress((i + 1)%3, i) += (1 - 2*Poi)*strain((i + 1)%3, i);
    stress((i + 2)%3, i) += (1 - 2*Poi)*strain((i + 2)%3, i);
  }
  stress *= You/(1+Poi)/(1-2*Poi);
  return stress;
}

void cons_law(const Matrix3d &strain, Matrix3d &stress, const Matrix3d &def_gra, const double &You, const double &Poi){
  //TODO:add more consititutive law

  double G = You/(2 + 2*Poi);
  double lam = You*Poi/(1+Poi)/(1-2*Poi);
  //St.Venant-Kirchhof model
  // double trace = strain(0, 0) + strain(1, 1) + strain(2, 2);
  // stress = def_gra*(2*G*strain + lam*trace*MatrixXd::Identity(3, 3));
  //Linear model
  double trace = def_gra(0, 0) + def_gra(1, 1) + def_gra(2, 2) - 3; 
  // stress = G*(def_gra + def_gra.transpose() - 2*Matrix3d::Identity()) + lam*trace*Matrix3d::Identity();
  stress.noalias() = def_gra*(2*G*strain + lam*trace*Matrix3d::Identity());

}




//set dim_ to 1 to debug

point_sys::point_sys(const MatrixXd  &points, const double &rho, const double &Young, const double &Poission, const double &vol_all, const double &kv, const vector<vector<size_t>> &friends, const VectorXd &sup_radi):
    dim_(points.cols()), points_(points), rho_(rho), Young_(Young), Poission_(Poission), vol_all_(vol_all), kv_(kv), friends_(friends), sup_radi_(sup_radi){
  
  mass_i_.setZero(dim_);
  //init
  double mass_total = rho*vol_all;
  double mass_sigma = (sup_radi_/3).array().cube().sum();

  mass_i_ =  (sup_radi_/3).array().cube();
  
  calc_weig();
  calc_rhoi_vi();
  scal_fac_ = dim_/rho_i_.sum();

  mass_i_ *= scal_fac_*rho_;
  rho_i_ *= scal_fac_*rho_;

}

double point_sys::get_mass(const size_t &i) const{
  return mass_i_(i);
}

size_t point_sys::Nx() const{
  return dim_ * 3;
}

int point_sys::calc_weig() const{
  // friends_ = vector<vector<size_t>>(dim_);
  weig_ = vector<vector<double>>(dim_);
 #pragma omp parallel for  
  for(size_t i = 0; i < dim_; ++i){
    vector<double> weig_of_one_p(friends_[i].size());
    for(size_t j = 0; j < friends_[i].size(); ++j){
      weig_of_one_p[j] = kernel(i, friends_[i][j]);
    }
    weig_[i] = weig_of_one_p;
  }
  
  return 0;
}

int point_sys::calc_rhoi_vi() const{
  //init
  rho_i_.setZero(dim_);
  vol_i_.setZero(dim_);

  
  assert(friends_.size() > 0 && weig_.size() > 0);
#pragma omp parallel for
  for(size_t i = 0; i < dim_; ++i){
    for(size_t j = 0; j < friends_[i].size(); ++j){
      rho_i_(i) += mass_i_(friends_[i][j])*weig_[i][j];
    }
  }  
  vol_i_ = mass_i_.array() / rho_i_.array();
  return 0;
}

double point_sys::kernel(const double &r, const double &h) const {
  if(r < h){
    double old_res =  315*pow((h*h - r*r), 3)/(64*PI*pow(h, 9));
    return old_res;
  }

  else
    return 0;
}

double point_sys::kernel(const size_t &i, const size_t &j) const{
  double r = (points_.col(j) - points_.col(i)).norm();
  return kernel(r, sup_radi_(i));
}


int point_sys::calc_defo_gra(const double *x, dat_ptr &data) const{
  Eigen::Map<const Matrix<double, Dynamic, Dynamic> > _x(x, 3, dim_);
#pragma omp parallel for
  for(size_t i = 0; i < dim_; ++i){
    Matrix3d one_du = Matrix3d::Ones();
    Matrix3d b = Matrix3d::Zero();
    for(size_t iter_j = 0; iter_j < friends_[i].size(); ++iter_j){
      Vector3d xij = points_.col(friends_[i][iter_j]) - points_.col(i);
      for(size_t k = 0; k < 3; ++k){
        b.col(k) += (_x(k, friends_[i][iter_j]) - _x(k, i))*weig_[i][iter_j]*xij;
      }
    }


    //clac inverse of sys_mat by SVD
    const Matrix3d& inv_A = dynamic_pointer_cast<energy_dat>(data)->inv_A_all_[i];
    for(size_t k = 0; k < 3; ++k){
      one_du.row(k) = inv_A * b.col(k);
    }
    one_du.transposeInPlace();
    // if(one_du.sum() != 0){
    //   cout << one_du << endl;
    // }
      
    Matrix3d F = one_du + Matrix3d::Identity();

    dynamic_pointer_cast<energy_dat>(data)->save_ele_def_gra(i, F);

  }
  return 0;
}


//calculate inv_A

int point_sys::pre_compute(dat_ptr &data) const {
 #pragma omp parallel for
  for(size_t i = 0; i < dim_; ++i){
    Matrix3d sys_mat;
    sys_mat.setZero(3, 3);
    
    
    // assert(friends_[i].size() >= 3);
    for(size_t iter_j = 0; iter_j < friends_[i].size(); ++iter_j){
      Vector3d xij = points_.col(friends_[i][iter_j]) - points_.col(i);
      dynamic_pointer_cast<energy_dat>(data)->sigma_w_points_.col(i) += -weig_[i][iter_j] * xij;
      sys_mat += weig_[i][iter_j]*xij*(xij.transpose());
    }
    auto inv_A = safe_inv(sys_mat);
    dynamic_pointer_cast<energy_dat>(data)->save_ele_inv_all(i, inv_A);
  }

  return 0;
}


int point_sys::Val(const double *x, dat_ptr &data)const {


  calc_defo_gra(x, data);

 #pragma omp parallel for 
  for(size_t i = 0; i < dim_; ++i){
    const Matrix3d& def_gra = dynamic_pointer_cast<energy_dat>(data)->def_gra_[i];

    //calculate strain and stress
    Matrix3d strain = def_gra.transpose()*def_gra - Matrix3d::Identity();
    Matrix3d stress = cons_law(strain, Young_, Poission_);

    //save energy
    double energy = 0.5*vol_i_(i)*(stress.array()*strain.array()).sum();
    energy += 0.5 * kv_ * pow(def_gra.determinant() - 1, 2);

    dynamic_pointer_cast<energy_dat>(data)->save_ele_strain(i, strain);
    dynamic_pointer_cast<energy_dat>(data)->save_ele_stress(i, stress);
    dynamic_pointer_cast<energy_dat>(data)->save_val(energy);    

  }

  return 0;
}


int point_sys::Gra(const double *x, dat_ptr &data) const{


#pragma omp parallel for
  for(size_t i = 0; i < dim_; ++i){

    const Matrix3d& def_gra = dynamic_pointer_cast<energy_dat>(data)->def_gra_[i];
    const Matrix3d& stress = dynamic_pointer_cast<energy_dat>(data)->stress_[i];
    const Matrix3d& strain = dynamic_pointer_cast<energy_dat>(data)->strain_[i];


    
    //calculate Fv
    Matrix3d gra_def_gra;
    auto trans_def_gra = def_gra.transpose();
    for(int j = 0; j < 3; ++j){
      Vector3d cross1 = trans_def_gra.col((j+1)%3), cross2 = trans_def_gra.col((j+2)%3);
      gra_def_gra.col(j) = cross1.cross(cross2);
    }
    gra_def_gra.transposeInPlace();

    dynamic_pointer_cast<energy_dat>(data)->save_ele_vol_cross(i, gra_def_gra);
    
    //assemble Fe and Fv
    const Matrix3d& inv_A = dynamic_pointer_cast<energy_dat>(data)->inv_A_all_[i];
    
    Matrix3d pre_F = vol_i_(i)*(2*def_gra*stress + kv_*(def_gra.determinant() - 1)*gra_def_gra)*inv_A;

    
    //add to gra_
    Vector3d di;
    
    di.setZero(3);
    Vector3d xij, force_j,force_i;
    for(size_t iter_j = 0; iter_j < friends_[i].size(); ++iter_j){
      double w = weig_[i][iter_j];
      xij = (points_.col(friends_[i][iter_j]) - points_.col(i));
      di += -w*xij;
      force_j = w*pre_F*xij;
      data->save_gra(friends_[i][iter_j], force_j);
    }
    
    force_i = pre_F*di;
    data->save_gra(i, force_i);

  }

  return 0;
}

//This is origin version
#if 1

int point_sys::Hes(const double*x, dat_ptr &data) const{
  
  auto start = system_clock::now();
  // omp_set_num_threads(4);  
#pragma omp parallel for
  for(size_t i = 0; i < dim_; ++i){
    // cout << "dim is " << i << endl;
    Matrix3d Kpq = Matrix3d::Zero();
    Matrix3d Kpq_vol = Matrix3d::Zero();
    Matrix3d K_local = Matrix3d::Zero();
    Matrix3d one_line;
    Matrix3d def_mult_dsig;

    Matrix3d dsigma_duk;
    Vector3d dp, dq, xij, xik;
    
    Matrix3d left_deri_det_matrix;
    Matrix3d right_deri;
    Vector3d next_line;
    Vector3d next_next_line;
    


    const auto& stress = dynamic_pointer_cast<energy_dat>(data)->stress_[i];
    const auto& def_gra = dynamic_pointer_cast<energy_dat>(data)->def_gra_[i];
    const auto& inv_A = dynamic_pointer_cast<energy_dat>(data)->inv_A_all_[i];
    const auto& vol_cross = dynamic_pointer_cast<energy_dat>(data)->vol_cross_[i];
    

    
    double def_gra_det = def_gra.determinant();

    vector<Matrix3d> ela_col(3);
    vector<Matrix3d> vol_col(3);

    for(size_t iter_k = 0; iter_k < friends_[i].size(); ++iter_k){
      xik = points_.col(friends_[i][iter_k]) - points_.col(i);
      dq = ( friends_[i][iter_k]== i) ? Vector3d(inv_A * dynamic_pointer_cast<energy_dat>(data)->sigma_w_points_.col(i)): inv_A * xik * weig_[i][iter_k];

      //elastic hessian
      for(size_t l = 0; l < 3; ++l){
        one_line.setZero(3, 3);
        one_line.row(l) = dq.transpose();

        dsigma_duk =
            cons_law(def_gra.row(l).transpose()*(dq.transpose()) + dq*def_gra.row(l), Young_, Poission_);
        ela_col[l] = 2*vol_i_(i)*(one_line*stress + def_gra*dsigma_duk);


      }

      //volume conserving forcing
      for(size_t l = 0; l < 3; ++l){
        left_deri_det_matrix = def_gra;
        left_deri_det_matrix.row(l) = dq.transpose();
        double left_deri = left_deri_det_matrix.determinant();


        right_deri = Matrix3d::Zero();{
          next_line = def_gra.row( (1 + 1) % 3 ).transpose();
          next_next_line = def_gra.row( (1 + 2) % 3 ).transpose();
          right_deri.row( (l + 1) % 3 ) = (next_next_line.cross(dq)).transpose();
          right_deri.row( (l + 2) % 3 ) = (dq.cross(next_line)).transpose();

        }
        vol_col[l] = kv_ * vol_i_(i) * (left_deri * vol_cross + (def_gra_det - 1) * right_deri);
      }

      for(size_t iter_j = 0; iter_j <= iter_k; ++iter_j){
 
        xij = points_.col(friends_[i][iter_j]) - points_.col(i);
        dp = (friends_[i][iter_j] == i) ? Vector3d(inv_A * dynamic_pointer_cast<energy_dat>(data)->sigma_w_points_.col(i)) : inv_A * xij * weig_[i][iter_j];
       
        for(size_t l = 0; l < 3; ++l){
          Kpq.col(l) = ela_col[l] * dp;
          Kpq_vol.col(l) = vol_col[l] * dp;
        
        }
        
        K_local = Kpq + Kpq_vol;

        
        data->save_hes(friends_[i][iter_j], friends_[i][iter_k], K_local);
        if(iter_k != iter_j)
          data->save_hes(friends_[i][iter_k], friends_[i][iter_j], K_local.transpose());

        
      }//for loop: p friends_[i].size()
    }//for_loop: q friends_[i].size()
    
  }//for_loop: i dim

  auto end = system_clock::now();
  auto duration = duration_cast<microseconds>(end - start);
  cout <<  "hessian花费了" 
       << double(duration.count()) * microseconds::period::num / microseconds::period::den 
       << "秒" << endl;


  return 0;

}//point_sys::Hessian
#endif


int point_sys::calc_Mass_matrix(){
  vector<Triplet<double>> mass_triplets(3 * dim_);
// #pragma omp parallel for
  for(size_t i = 0; i < dim_; ++i){
    for(size_t j = 0; j < 3; ++j){
      mass_triplets[i*3 + j] = Triplet<double>(i * 3 + j, i * 3 + j, mass_i_(i));
    }
  }
  M_ = SparseMatrix<double>(3*dim_, 3*dim_);
  M_.setFromTriplets(mass_triplets.begin(), mass_triplets.end());
  return 0;
 }

const SparseMatrix<double>& point_sys::get_Mass_Matrix(){
  calc_Mass_matrix();
  return M_;
}


const VectorXd point_sys::get_Mass_VectorXd(){
  return mass_i_;
}
// #endif


}//namespace : marvel


