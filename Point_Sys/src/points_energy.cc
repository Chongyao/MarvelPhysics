#include <math.h>
#include <iostream>

#include "points_energy.h"
#include "get_nn.h"

#include <Eigen/SVD>
#include <Eigen/LU>
#include <Eigen/Geometry>

using namespace std;
using namespace Eigen;

#define AUTODIFF_ENABLE_EIGEN_SUPPORT
#include <autodiff/autodiff/autodiff.hpp>
using namespace autodiff;

#define PI 3.14159265359


namespace marvel{
#if 0
var STRAIN_ENERGY(const VectorXv& u_all, const Matrix3d &inv_A, const MatrixXd &x){
  size_t num_friends = u_all.size()/3;
  
  Matrix3v def_u;
  for(size_t i = 0; i < num_friends; ++i){
    Vector3d xij = x.col(i) - x.col(0);
    def_u.col(i) = 
  }

}
#endif


Matrix3d safe_inv(const MatrixXd& sys_mat){
  //TODO: use better solver considering the sysmertic
  JacobiSVD<MatrixXd> svd(sys_mat, ComputeThinU | ComputeThinV);

  auto sin_val = svd.singularValues();
  Matrix3d inv_sin_val;
  inv_sin_val.setZero(3, 3);
  
  for(size_t p = 0; p < 3; ++p){
    inv_sin_val(p, p) = sin_val(p)>0?1/sin_val(p):0;
  }
  
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
  stress = def_gra*(2*G*strain + lam*trace*Matrix3d::Identity());

}




//set dim_ to 1 to debug
point_sys::point_sys(const MatrixXd  &points, const double &rho, const double &Young, const double &Poission, const double &vol_all, const double &kv, const vector<vector<size_t>> &friends, const VectorXd &sup_radi):
    points_(points), rho_(rho), Young_(Young), Poission_(Poission), vol_all_(vol_all), dim_(points.cols()), kv_(kv), friends_(friends), sup_radi_(sup_radi){
  
  mass_i_.setZero(dim_);
  //init
  double mass_total = rho*vol_all;
  double mass_sigma = (sup_radi_/3).array().cube().sum();

  // scal_fac_ = mass_total/mass_sigma;
  // scal_fac_ /= 19.2649;
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
  return dim_;
}
int point_sys::calc_weig() const{
  // friends_ = vector<vector<size_t>>(dim_);
  weig_ = vector<vector<double>>(dim_);
//#pragma omp parallel for  
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
//#pragma omp parallel for
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

int point_sys::calc_defo_gra(const double *disp, energy_dat &dat_str) const{
  Map<const Matrix<double, Dynamic, Dynamic> > _disp(disp, 3, dim_);
//#pragma omp parallel for
  for(size_t i = 0; i < dim_; ++i){
    Matrix3d one_du;
    Matrix3d b = Matrix3d::Zero();
    for(size_t iter_j = 0; iter_j < friends_[i].size(); ++iter_j){
      Vector3d xij = points_.col(friends_[i][iter_j]) - points_.col(i);
      for(size_t k = 0; k < 3; ++k){
        b.col(k) += (_disp(k, friends_[i][iter_j]) - _disp(k, i))*weig_[i][iter_j]*xij;
        
      }
    }
    
    //clac inverse of sys_mat by SVD
    
    Map<Matrix3d> inv_A(dat_str.inv_A_all_.col(i).data());
    
    
    for(size_t k = 0; k < 3; ++k){
      one_du.row(k) = (inv_A * b.col(k)).transpose();
    }
    MatrixXd F = one_du + MatrixXd::Identity(3, 3);
    dat_str.save_ele_def_gra(i, F);

  }
  return 0;
}


//calculate inv_A
int point_sys::pre_compute(energy_dat &dat_str) const {
 #pragma omp parallel for
  for(size_t i = 0; i < dim_; ++i){
    Matrix3d sys_mat;
    sys_mat.setZero(3, 3);
    
    
    // assert(friends_[i].size() >= 3);
    for(size_t iter_j = 0; iter_j < friends_[i].size(); ++iter_j){
      Vector3d xij = points_.col(friends_[i][iter_j]) - points_.col(i);
      dat_str.sigma_w_points_.col(i) += -weig_[i][iter_j] * xij;
      sys_mat += weig_[i][iter_j]*xij*(xij.transpose());
    }
    auto inv_A = safe_inv(sys_mat);
    
    dat_str.save_ele_inv_all(i, inv_A);
  }
  return 0;
}
int point_sys::Val(const double *disp, energy_dat &dat_str)const {
  calc_defo_gra(disp, dat_str);
  // dat_str.Val_ = 0;
#pragma omp parallel for 
  for(size_t i = 0; i < dim_; ++i){
    Map<MatrixXd> def_gra(dat_str.def_gra_.col(i).data(), 3, 3);
    
    //calculate strain and stress
    Matrix3d strain = def_gra.transpose()*def_gra - Matrix3d::Identity();
    // Matrix3d strain = 0.5*(def_gra + def_gra.transpose()) - Matrix3d::Identity();
    
    dat_str.save_ele_strain(i, strain);
    // Matrix3d stress;
    // cons_law(strain, stress, def_gra, Young_, Poission_);

    Matrix3d stress = cons_law(strain, Young_, Poission_);
    // assert(stress_test == stress);
    dat_str.save_ele_stress(i, stress);

    //save energy
    double energy = 0.5*vol_i_(i)*(stress.array()*strain.array()).sum();
#pragma omp critical
    {
      dat_str.save_val(energy);      
    }


    // dat_str.Val_ += 0.5*kv_*pow((def_gra.determinant() - 1), 2);
    // dat_str.ela_val_(i) += 0.5*vol_i_(i)*(stress.array()*strain.array());
    // dat_str.vol_val_(i) += 0.5*kv_*pow((def_gra.determinant() - 1), 2);
  }
  cout << "elasiticy energy" << dat_str.Val_ << endl;
 
}
int point_sys::Gra(const double *disp, energy_dat &dat_str) const{
#pragma omp parallel for
  for(size_t i = 0; i < dim_; ++i){
    Map<const MatrixXd> def_gra(dat_str.def_gra_.col(i).data(), 3, 3);
    Map<const MatrixXd> stress(dat_str.stress_.col(i).data(), 3, 3);
    Map<const MatrixXd> strain(dat_str.strain_.col(i).data(), 3, 3);
    
    //calculate Fv
    Matrix3d gra_def_gra;
    auto trans_def_gra = def_gra.transpose();
    for(int j = 0; j < 3; ++j){
      Vector3d cross1 = trans_def_gra.col((j+1)%3), cross2 = trans_def_gra.col((j+2)%3);
      gra_def_gra.col(j) = cross1.cross(cross2);
    }
    gra_def_gra.transposeInPlace();

    //assemble Fe and Fv
    Map<MatrixXd> inv_A(dat_str.inv_A_all_.col(i).data(), 3, 3);
    Matrix3d pre_F = -vol_i_(i)*(2*def_gra*stress + kv_*(def_gra.determinant() - 1)*gra_def_gra)*inv_A;
#pragma omp critical
    {
    dat_str.save_ele_pre_F(i, pre_F);      
    }

    
    //add to gra_
    Vector3d di;
    di.setZero(3);
    Vector3d xij, force_j,force_i;
    for(size_t iter_j = 0; iter_j < friends_[i].size(); ++iter_j){
      double w = weig_[i][iter_j];
      xij = (points_.col(friends_[i][iter_j]) - points_.col(i));
      di += -w*xij;
      force_j = w*pre_F*xij;
#pragma omp critical
      {
      dat_str.save_ele_gra(friends_[i][iter_j], force_j);      
      }

    }
    
    force_i = pre_F*di;
#pragma omp critical
    {
    dat_str.save_ele_gra(i, force_i);      
    }

  }

  return 0;
}
int point_sys::gravity(const double *x, energy_dat &dat_str,  const double &gravity) const{

  MatrixXd g(3, dim_);
  g.setZero(3, dim_);
  g.row(2) = VectorXd::Constant(dim_, -gravity).cwiseProduct(mass_i_).transpose();

                                                       
  Map<MatrixXd> Gra(dat_str.gra_.data(), 3, dim_);
  Gra += g;
  return 0;
}

// #ipf 0
int point_sys::Hessian(const double*disp, energy_dat &dat_str){
  
  vector<Triplet<double>> TripletList;

  //TODO: estimate entries
  // TripletList.reserve()

  //TODO:consider sysmetric
  
  Matrix3d Kpq = Matrix3d::Zero();
  Matrix3d one_line;
#pragma omp parallel for
  for(size_t i = 0; i < dim_; ++i){
    Map<const MatrixXd> stress(dat_str.stress_.col(i).data(), 3, 3);
    Map<const MatrixXd> def_gra(dat_str.def_gra_.col(i).data(), 3, 3);
    Map<const MatrixXd> inv_A(dat_str.inv_A_all_.col(i).data(), 3, 3);
    
    // vector<Triplet<double>> test_Trip;
    
    for(size_t iter_j = 0; iter_j < friends_[i].size(); ++iter_j){
      auto xij = points_.col(friends_[i][iter_j]) - points_.col(i);
      auto dp = (friends_[i][iter_j] == i) ? Vector3d(inv_A * dat_str.sigma_w_points_.col(i)) : inv_A * xij * weig_[i][iter_j];
      
      for(size_t iter_k = 0; iter_k < friends_[i].size(); ++iter_k){
        auto xik = points_.col(friends_[i][iter_k]) - points_.col(i);
        auto dq = ( friends_[i][iter_k]== i) ? Vector3d(inv_A * dat_str.sigma_w_points_.col(i)): inv_A * xik * weig_[i][iter_k];
        
        for(size_t l = 0; l < 3; ++l){
          one_line.setZero(3, 3);
          one_line.row(l) = dq.transpose();
          Matrix3d dsigma_duk =
              cons_law(def_gra.row(l).transpose()*dq.transpose() + dq*def_gra.row(l), Young_, Poission_);
          

          Kpq.col(l) = 2*vol_i_(i)*(one_line*stress + def_gra*dsigma_duk)*dp;
        }



#pragma omp critical
        {
        for(size_t m = 0; m < 3; ++m){
          for(size_t n = 0; n < 3; ++n){
            if (Kpq(m, n) != 0){
              dat_str.hes_trips.push_back(Triplet<double>(friends_[i][iter_j]*3 + m, friends_[i][iter_k]*3 + n, Kpq(m,n)));
              // test_Trip.push_back(Triplet<double>(friends_[i][iter_j]*3 + m, friends_[i][iter_k]*3 + n, Kpq(m,n)));
               
            }
          }
        }          
        }
        

      }
    }

  //   SparseMatrix<double> EL_Hessian_test(dim_*3, dim_*3);
  // EL_Hessian_test.setFromTriplets(test_Trip.begin(), test_Trip.end());
  // MatrixXd dense_test = MatrixXd(EL_Hessian_test);
  // cout << "[INFO]>>>>>>>>>>>>>>>>>>>test elas hessian<<<<<<<<<<<<<<<<<<" << endl;
  // cout << dense_test << endl << endl<<dense_test - dense_test.transpose() << endl;
  }

}

int point_sys::calc_Mass_matrix(){
  vector<Triplet<double>> mass_triplets(3 * dim_);
//#pragma omp parallel for
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


const VectorXd& point_sys::get_Mass_VectorXd(){
  return mass_i_;
}
// #endif


}//namespace : marvel


