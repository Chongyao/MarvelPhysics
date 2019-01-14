#include <string>
#include <iostream>
#include <chrono>

#include <libigl/include/igl/readOBJ.h>
#include <libigl/include/igl/writeOBJ.h>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/SparseCore>
#include <Eigen/SparseCholesky>
#include<Eigen/SparseLU>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/filesystem.hpp>

#include "Point_Sys/src/geometry.h"
#include "Point_Sys/src/gen_points.h"
#include "Point_Sys/src/get_nn.h"
#include "Point_Sys/src/points_energy.h"
#include "Point_Sys/src/data_stream.h"
#include "Point_Sys/src/gen_surf.h"
#include "io.h"
#include "Point_Sys/src/basic_energy.h"



using namespace marvel;
using namespace std;
using namespace Eigen;
using namespace igl;
using namespace chrono;
using namespace boost;

int main(int argc, char** argv){

  Eigen::initParallel();
  cout << "[INFO]>>>>>>>>>>>>>>>>>>>Eigen parallel<<<<<<<<<<<<<<<<<<" << endl;
  cout << "enable parallel in Eigen in " << nbThreads() << " threads" << endl;
  
  cout << "[INFO]>>>>>>>>>>>>>>>>>>>READ JSON FILE<<<<<<<<<<<<<<<<<<" << endl;
  boost::property_tree::ptree pt;{
    const string jsonfile_path = argv[1];
    
    cout << jsonfile_path << endl;
    const size_t ext = jsonfile_path.rfind(".json");
    if (ext != std::string::npos){
      read_json(jsonfile_path, pt);
      cout << "read json successful" <<endl;
    }
    else{
      cout << "json file extension error" << endl;
      return 0;
    }
  }
  
  cout << "[INFO]>>>>>>>>>>>>>>>>>>>IMPORT MESH<<<<<<<<<<<<<<<<<<" << endl;
  const string mesh_name = pt.get<string>("surf");
  const string indir = pt.get<string>("indir");
  const string outdir = pt.get<string>("outdir") + mesh_name;
  //mkdir
  boost::filesystem::path outpath(outdir);
  if ( !boost::filesystem::exists(outdir) )
    boost::filesystem::create_directories(outdir);

  MatrixXi surf;
  MatrixXd nods;
  readOBJ((indir+mesh_name+".obj").c_str(), nods, surf);
  cout << "surf: " << surf.rows() << " " << surf.cols() << endl << "nods: " << nods.rows() << " " << nods.cols() << endl;
  
  surf.transposeInPlace();
  nods.transposeInPlace();


  cout << "[INFO]>>>>>>>>>>>>>>>>>>>Generate sampled points<<<<<<<<<<<<<<<<<<" << endl;
  MatrixXd points(3,3);
  MatrixXd test(3, 3);
  gen_points(nods, surf, pt.get<size_t>("num_in_axis"), points, false);

  
  cout << "[INFO]>>>>>>>>>>>>>>>>>>>points<<<<<<<<<<<<<<<<<<" << endl;
  cout << points << endl;
  size_t dim = points.cols();
  cout <<"generate points done." << endl;
  
  cout << "[INFO]>>>>>>>>>>>>>>>>>>>Build spatial hash<<<<<<<<<<<<<<<<<<" << endl;
  spatial_hash SH(points, pt.get<size_t>("nn_num"));

  
  cout << "[INFO]>>>>>>>>>>>>>>>>>>>Build Point System<<<<<<<<<<<<<<<<<<" << endl;
  //calc volume 
  double volume = clo_surf_vol(nods, surf);
  //calc support radii
  VectorXd sup_radi = SH.get_sup_radi();
  
  cout << "[INFO]>>>>>>>>>>>>>>>>>>>sup_radi<<<<<<<<<<<<<<<<<<" << endl;
  cout << sup_radi << endl;
  
  //get friends of every point
  vector<vector<size_t>> friends_all(dim);
#pragma omp parallel for
  for(size_t i = 0; i < dim; ++i){
    SH.get_friends(points.col(i), sup_radi(i), friends_all[i]);
  }

  point_sys PS(points, pt.get<double>("rho"), pt.get<double>("Young"), pt.get<double>("Poission"), volume, pt.get<double>("kv"), friends_all, sup_radi);

  
  cout << "[INFO]>>>>>>>>>>>>>>>>>>>Construct Deform_Surf<<<<<<<<<<<<<<<<<<" << endl;
  //get friends of every vertex in mesh
  double kernel_cof = pt.get<double>("ker_cof");
  vector<vector<size_t>> vet_fris(nods.cols());
  for(size_t j = 0; j < nods.cols(); ++j){
    vector<size_t> fris;
    SH.get_friends(nods.col(j), kernel_cof, fris, false);
    vet_fris[j] = fris;
  }
  
  deform_surf_MLS<double> DS(surf, nods, points, vet_fris, kernel_cof);



  cout << "[INFO]>>>>>>>>>>>>>>>>>>>Simple Constraint Points<<<<<<<<<<<<<<<<<<" << endl;

  //add simple constraints
  //This should read from file. We loop for some points to restrain here.
  //Constraints vary from different models and situations.
  vector<size_t> cons;
  for(size_t i = 0; i < points.cols(); ++i){
    if(points(2, i) > 0.7 ){
      cons.push_back(i);
      cout << i << " ";
    }
  }
  cout << endl;
  position_constraint pos_cons(pt.get<double>("position_weig"), cons, dim);

  cout << "[INFO]>>>>>>>>>>>>>>>>>>>Gravity<<<<<<<<<<<<<<<<<<" << endl;
  double gravity = pt.get<double>("gravity");
  gravity_energy GE(pt.get<double>("w_g"), gravity, dim, PS.get_Mass_VectorXd(), 'z');
  
  
  cout << "[INFO]>>>>>>>>>>>>>>>>>>>SOlVE<<<<<<<<<<<<<<<<<<" << endl;
  //initilize variables in time integration
  energy_dat dat_str (dim);

  double delt_t = pt.get<double>("time_step");
  MatrixXd displace;
  MatrixXd velocity;
  MatrixXd acce;
  MatrixXd new_acce;
  MatrixXd gra;
  MatrixXd vet_displace;
  displace.setZero(3, dim);
  velocity.setZero(3, dim);
  acce.setZero(3, dim);
  gra.setZero(3, dim);
  new_acce.setZero(3, dim);
  vet_displace.setZero(3, nods.cols());

  PS.pre_compute(dat_str);
  size_t iters_perframe = floor(1.0/delt_t/pt.get<size_t>("rate"));

  SparseMatrix<double> A_CG(dim * 3, dim * 3);
  VectorXd b_CG(dim * 3);
  SparseMatrix<double> M = PS.get_Mass_Matrix();

  MatrixXd displace_dyna(3, dim - cons.size());
  
  for(size_t i = 0; i < pt.get<size_t>("max_iter"); ++i){
    cout << "iter is "<<endl<< i << endl;


    auto point_cons_ptr = cons.begin();
    size_t count = 0;
    for(size_t i = 0; i < dim - cons.size(); ++i){
      if(*point_cons_ptr != count){
        displace_dyna.col(i) = displace.col(count);
        --i;
      }
      else
        ++point_cons_ptr;
      ++count;
    }
    cout << "displace is " << endl<< displace_dyna << endl;
    // cout << "velocity is "<<endl<< velocity.block(0, 0, 3, 8) << endl;

    PS.Val(displace.data(), dat_str);
    
    MatrixXd temp_gra = dat_str.gra_;
    PS.Gra(displace.data(), dat_str);
    temp_gra = dat_str.gra_ - temp_gra;
    cout << "elasticity gra " << temp_gra.array().square().sum()<<endl;
    PS.Hessian(displace.data(), dat_str);
    
    GE.Val(displace.data(), dat_str);
    temp_gra = dat_str.gra_;
    GE.Gra(displace.data(), dat_str);
    temp_gra = dat_str.gra_ - temp_gra;
    cout << "gravity gra " << temp_gra.array().square().sum()<<endl;

    temp_gra = dat_str.gra_;
    pos_cons.Gra(displace.data(), dat_str);
    temp_gra = dat_str.gra_ - temp_gra;
    cout << "position gra " << temp_gra.array().square().sum()<<endl;
    
    pos_cons.Hes(displace.data(),dat_str);
    cout << "[INFO]>>>>>>>>>>>>>>>>>>>gra<<<<<<<<<<<<<<<<<<" << endl;
    // cout << dat_str.gra_.block(0, 0, 3, 8) << endl;
    // cout << dat_str.gra_.array().square().sum() << endl;
    

    //implicit time integral
    #if 0
    { 
      A_CG.setZero();
      Map<VectorXd> _velo(velocity.data(), 3*dim);
      Map<VectorXd> _F(dat_str.gra_.data(), 3*dim);
      dat_str.hes_.setFromTriplets(dat_str.hes_trips.begin(), dat_str.hes_trips.end());
      A_CG = M + delt_t*delt_t*dat_str.hes_;
      b_CG = M * _velo + delt_t * _F;

      // cout << "[INFO]>>>>>>>>>>>>>>>>>>>A_CG<<<<<<<<<<<<<<<<<<" << endl;      
      // ConjugateGradient<SparseMatrix<double>, Lower|Upper> cg;
      // // cg.setMaxIterations(50);
      // // cg.setTolerance(1e-40);
      // cg.compute(A_CG);
      // _velo = cg.solve(b_CG);
      // cout << "#iterations:     " << cg.iterations() << endl;
      // cout << "estimated error: " << cg.error()      << endl;


      
      // cout << "[INFO]>>>>>>>>>>>>>>>>>>>LLT<<<<<<<<<<<<<<<<<<" << endl;
      // SimplicialLLT<SparseMatrix<double>,Lower> llt;
      // llt.compute(A_CG);
      // if(llt.info() != Success){
      //   cout << "Not SPD !!!" << endl;
      // }
      // _velo = llt.solve(b_CG);
      // if(llt.info()!=Success){
      //   cout << "Solve fail" << endl;
      // }

      
      cout << "[INFO]>>>>>>>>>>>>>>>>>>>LU<<<<<<<<<<<<<<<<<<" << endl;
      SparseLU<SparseMatrix<double>> lu;
      lu.compute(A_CG);
      _velo = lu.solve(b_CG);
      
      cout << velocity.block(0, 0, 3, 10) << endl;
      displace += delt_t * velocity;
    }
    #endif

    {
      A_CG.setZero();
      Map<VectorXd> _disp(displace.data(), 3*dim);
      Map<VectorXd> _velo(velocity.data(), 3*dim);
      Map<VectorXd> _F(dat_str.gra_.data(), 3*dim);
      dat_str.hes_.setFromTriplets(dat_str.hes_trips.begin(), dat_str.hes_trips.end());
      A_CG = M + delt_t*delt_t*dat_str.hes_;
      b_CG = delt_t * M * _velo + delt_t * delt_t *  _F + (M + delt_t * delt_t * dat_str.hes_) * _disp;

      cout << "[INFO]>>>>>>>>>>>>>>>>>>>LU<<<<<<<<<<<<<<<<<<" << endl;
      SparseLU<SparseMatrix<double>> lu;
      lu.compute(A_CG);
      _disp = lu.solve(b_CG);

      
    }

    
    
    cout << "[INFO]>>>>>>>>>>>>>>>>>>>Elasticity Energy Val<<<<<<<<<<<<<<<<<<" << endl;
    cout << dat_str.Val_ << endl;
    
    cout << "[INFO]>>>>>>>>>>>>>>>>>>>GRA<<<<<<<<<<<<<<<<<<" << endl;
    cout << dat_str.gra_.array().square().sum() << endl;
    cout << endl << endl << endl;
    // cout << "[INFO]>>>>>>>>>>>>>>>>>>>VOL conservation val<<<<<<<<<<<<<<<<<<" << endl;
    // cout << dat_str.vol_val_.transpose();
        
    // vet_displace = DS.update_surf(displace, dat_str.def_gra_);
    // if(i%iters_perframe == 0){
      auto surf_filename = outdir  + "/" + mesh_name + "_" + to_string(i) + ".vtk";
      auto point_filename = outdir + "/" + mesh_name + "_points_" + to_string(i) + ".vtk";

      MatrixXd points_now = points + displace;

      point_write_to_vtk(point_filename.c_str(), points_now.data(), dim);
      // point_vector_append2vtk(false, point_filename.c_str(), velocity, dim, "velocity");
      point_scalar_append2vtk(true, point_filename.c_str(), dat_str.ela_val_, dim, "strain_Energy");
      point_scalar_append2vtk(true, point_filename.c_str(), dat_str.vol_val_, dim, "vol_conservation_Energy");

      // vet_displace = displace.block(0, 0, 3, nods.cols());
      // tri_mesh_write_to_vtk(surf_filename.c_str(), nods + vet_displace, surf);
    // }

    dat_str.set_zero();
  }
  //done
}



