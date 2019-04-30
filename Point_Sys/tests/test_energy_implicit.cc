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
  gen_points(nods, surf, pt.get<size_t>("num_in_axis"), points, true);

  
  cout << "[INFO]>>>>>>>>>>>>>>>>>>>points<<<<<<<<<<<<<<<<<<" << endl;

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

  
  //get friends of every point
  vector<vector<size_t>> friends_all(dim);
#pragma omp parallel for
  for(size_t i = 0; i < dim; ++i){
    SH.get_friends(points.col(i), sup_radi(i), friends_all[i]);
  }

  point_sys PS(points, pt.get<double>("rho"), pt.get<double>("Young"), pt.get<double>("Poission"), volume, pt.get<double>("kv"), friends_all, sup_radi);

  
  cout << "[INFO]>>>>>>>>>>>>>>>>>>>Simple Constraint Points<<<<<<<<<<<<<<<<<<" << endl;

  //add simple constraints
  //This should read from file. We loop for some points to restrain here.
  //Constraints vary from different models and situations.
  vector<size_t> cons;
  auto cons_file_path = indir + mesh_name +".csv";
  if ( boost::filesystem::exists(cons_file_path) )
    read_fixed_verts_from_csv(cons_file_path.c_str(), cons);
  cout << "constrint " << cons.size() << " points" << endl;
  position_constraint pos_cons(pt.get<double>("position_weig"), cons, dim);

  cout << "[INFO]>>>>>>>>>>>>>>>>>>>Gravity<<<<<<<<<<<<<<<<<<" << endl;
  double gravity = pt.get<double>("gravity");
  gravity_energy GE(pt.get<double>("w_g"), gravity, dim, PS.get_Mass_VectorXd(), 'y');
  
  
  cout << "[INFO]>>>>>>>>>>>>>>>>>>>MOMENTUM<<<<<<<<<<<<<<<<<<" << endl;
  double delt_t = pt.get<double>("time_step");
  momentum MO(dim, PS.get_Mass_Matrix(), delt_t);
  



  
  cout << "[INFO]>>>>>>>>>>>>>>>>>>>SOlVE<<<<<<<<<<<<<<<<<<" << endl;
  //initilize variables in time integration
  energy_dat dat_str (dim);


  MatrixXd displace;
  MatrixXd vet_displace;
  displace.setZero(3, dim);
  vet_displace.setZero(3, nods.cols());


  PS.pre_compute(dat_str);
  size_t iters_perframe = floor(1.0/delt_t/pt.get<size_t>("rate"));
  VectorXd solution = VectorXd::Zero(3 * dim);
  VectorXd displace_search = VectorXd::Zero(3 * dim);

  double d1dtdt = 1 / delt_t /delt_t, d1dt = 1 / delt_t;

  Map<const VectorXd>res(dat_str.gra_.data(), 3 * dim);
  MatrixXd points_now;
  
  for(size_t i = 0; i < pt.get<size_t>("max_iter"); ++i){
    cout << "iter is "<<endl<< i << endl;
    // cout << "displace is " << endl<< displace.block(0, 0, 3, 8) << endl;

    //newtown iter
    Map<VectorXd> displace_plus(displace.data(), 3*dim);

    for(size_t newton_i = 0; newton_i < 20; ++newton_i){
      dat_str.set_zero();
      cout << "newton iter is "<< newton_i << endl;
      
      MO.Val(displace_plus.data(), dat_str);
      MO.Gra(displace_plus.data(), dat_str);
      MO.Hes(displace_plus.data(), dat_str);
      
      PS.Val(displace_plus.data(), dat_str);
      PS.Gra(displace_plus.data(), dat_str);
      PS.Hessian(displace_plus.data(), dat_str);


      GE.Val(displace_plus.data(), dat_str);
      GE.Gra(displace_plus.data(), dat_str);

      pos_cons.Gra(displace_plus.data(), dat_str);
      pos_cons.Hes(displace_plus.data(),dat_str);



      double res_value = res.array().square().sum();
      if(res_value < 1e-10){
        cout << "[INFO]Newton res " <<endl << res_value << endl;;
        cout << "[INFO]>>>>>>>>>>>>>>>>>>>ALL Energy Val<<<<<<<<<<<<<<<<<<" << endl;
        cout << dat_str.Val_ << endl;
        dat_str.set_zero();
        break;
      }
      
      
      //implicit time integral
      dat_str.hes_.setFromTriplets(dat_str.hes_trips.begin(), dat_str.hes_trips.end());
             
      // cout << "[INFO]>>>>>>>>>>>>>>>>>>>LLT<<<<<<<<<<<<<<<<<<" << endl;
      SimplicialLLT<SparseMatrix<double>> llt;
      llt.compute(dat_str.hes_);
      VectorXd all_one = VectorXd::Ones(3 * dim);
      while(llt.info() != Eigen::Success){
        cout <<"lltinfo "<< llt.info() << endl;
        dat_str.hes_ += all_one.asDiagonal();
        llt.compute(dat_str.hes_);
        all_one *= 2;
      }
      solution += llt.solve(-res);


      // cout << "[INFO]>>>>>>>>>>>>>>>>>>>A_CG<<<<<<<<<<<<<<<<<<" << endl;
      // ConjugateGradient<SparseMatrix<double>, Lower|Upper> cg;
      // cg.setMaxIterations(3*dim);
      // cg.setTolerance(1e-8);
      // cg.compute(A_CG);
      // disp_t_plus += cg.solve(b_CG);

      
      // {//Line search
      //   const double c = 1e-2;
      //   double alpha = 1 / 0.5;

      //   double Val_init = dat_str.Val_, down = res.dot(solution), Val_upbound, Val_func;
      //   cout << "down is "<< down << endl;

      //   do{
      //     alpha *= 0.5;
      //     Val_upbound = Val_init - c * alpha * down;
      //     displace_search = displace_plus + alpha * solution;
      //     dat_str.Val_ = 0;
          
      //     MO.Val(displace_search.data(), dat_str);
      //     PS.Val(displace_search.data(), dat_str);
      //     GE.Val(displace_search.data(), dat_str);
      //     pos_cons.Val(displace_search.data(), dat_str);

      //     Val_func = dat_str.Val_;
      //     cout << "Val_func is " << Val_func << " UP bound is " << Val_upbound << endl;

      //   }while(Val_func > Val_upbound);
      //   displace_plus = displace_search;
          
      // }

      displace_plus += solution;
      

    }
    
    MO.update_location_and_velocity(displace_plus.data());


    auto surf_filename = outdir  + "/" + mesh_name + "_" + to_string(i) + ".vtk";
    auto point_filename = outdir + "/" + mesh_name + "_points_" + to_string(i) + ".vtk";

    MatrixXd points_now = points + displace;

    point_write_to_vtk(point_filename.c_str(), points_now.data(), dim);




  }
  //done

  return 0;
}

  




