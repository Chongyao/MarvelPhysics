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
#include "basic_energy.h"
#include "config.h"


using namespace marvel;
using namespace std;
using namespace Eigen;
using namespace igl;
using namespace chrono;
using namespace boost;

int main(int argc, char** argv){


  
  __TIME_BEGIN__
  Eigen::initParallel();
  cout << "[INFO]>>>>>>>>>>>>>>>>>>>Eigen parallel<<<<<<<<<<<<<<<<<<" << endl;
  cout << "enable parallel in Eigen in " << nbThreads() << " threads" << endl;
  __TIME_END__("hey")
  
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
  cout << "[INFO]points num is" << points.cols() << endl;
  
  cout << "[INFO]>>>>>>>>>>>>>>>>>>>points<<<<<<<<<<<<<<<<<<" << endl;

  size_t dim = points.cols();
  cout <<"generate points done." << endl;


  cout << "[INFO] Assemble energies..." << endl;
  enum {POTS, CONS, GRAV, MOME};
  vector<std::shared_ptr<Functional<double, 3>>> ebf(MOME + 1); 
  

  
  
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

  ebf[POTS] = make_shared<point_sys>(points, pt.get<double>("rho"), pt.get<double>("Young"), pt.get<double>("Poission"), volume, pt.get<double>("kv"), friends_all, sup_radi);


  
  cout << "[INFO]>>>>>>>>>>>>>>>>>>>Simple Constraint Points<<<<<<<<<<<<<<<<<<" << endl;
  vector<size_t> cons;
  auto cons_file_path = indir + mesh_name +".csv";
  if ( boost::filesystem::exists(cons_file_path) )
    read_fixed_verts_from_csv(cons_file_path.c_str(), cons);
  cout << "constrint " << cons.size() << " points" << endl;

  
  ebf[CONS] = std::make_shared<position_constraint<3>>(dim, pt.get<double>("position_weig"), cons);
    

  cout << "[INFO]>>>>>>>>>>>>>>>>>>>Gravity<<<<<<<<<<<<<<<<<<" << endl;
  double gravity = pt.get<double>("gravity");
  const auto mass_vector = dynamic_pointer_cast<point_sys>(ebf[POTS])->get_Mass_VectorXd();
  ebf[GRAV] = make_shared<gravity_energy<3>>(dim, pt.get<double>("w_g"), gravity,  mass_vector, 'y');
  
  
  cout << "[INFO]>>>>>>>>>>>>>>>>>>>MOMENTUM<<<<<<<<<<<<<<<<<<" << endl;
  double delt_t = pt.get<double>("time_step");
  // momentum MO(dim, PS.get_Mass_Matrix(), delt_t);
  ebf[MOME] = make_shared<momentum<3>>(dim, mass_vector, delt_t);

  
  // cout << "[INFO]>>>>>>>>>>>>>>>>>>>Collision<<<<<<<<<<<<<<<<<<" << endl;
  // collision COLL(pt.get<double>("w_coll"), 'y', pt.get<double>("coll_pos"), static_cast<size_t>(nods.cols()), dim);

  //energy all
  std::shared_ptr<Functional<double, 3>> energy;
  try {
    energy = make_shared<energy_t<double, 3>>(ebf);

  } catch ( std::exception &e ) {
    cerr << e.what() << endl;
    exit(EXIT_FAILURE);
  }






  
  cout << "[INFO]>>>>>>>>>>>>>>>>>>>SOlVE<<<<<<<<<<<<<<<<<<" << endl;
  //initilize variables in time integration
  std::shared_ptr<dat_str_core<double, 3>>  dat_str = make_shared<energy_dat>(dim);
  dynamic_pointer_cast<point_sys>(ebf[POTS])->pre_compute(dat_str);{
    VectorXd random_x(3 * dim);{
      #pragma omp parallel for
      for(size_t i = 0; i < 3 * dim; ++i){
        random_x(i) = i * 4.5 + i * i ;
      }
    }
    dat_str->set_zero();
    const auto& sm1 = dat_str->get_hes();
    cout<<"the number of nonzeros with comparison: \n"
        << (Eigen::Map<Eigen::VectorXd> (sm1.valuePtr(), sm1.nonZeros()).array() != 0).count()
        << endl;

    energy->Val(random_x.data(), dat_str);
    energy->Gra(random_x.data(), dat_str);
    energy->Hes(random_x.data(), dat_str);
    cout<<"the number of nonzeros with comparison: \n"
        << (Eigen::Map<Eigen::VectorXd> (sm1.valuePtr(), sm1.nonZeros()).array() != 0).count()
        << endl;
    dat_str->set_zero_after_pre_compute();
    dat_str->set_zero();

  }


  MatrixXd displace;
  MatrixXd vet_displace;
  displace.setZero(3, dim);
  vet_displace.setZero(3, nods.cols());

  
  size_t iters_perframe = floor(1.0/delt_t/pt.get<size_t>("rate"));
  VectorXd solution = VectorXd::Zero(3 * dim);
  VectorXd displace_search = VectorXd::Zero(3 * dim);
  double d1dtdt = 1 / delt_t /delt_t, d1dt = 1 / delt_t;
  Map<const VectorXd>res(dat_str->get_gra().data(), 3 * dim);


  MatrixXd points_now;
  
  
  for(size_t i = 0; i < pt.get<size_t>("max_iter"); ++i){
    
    cout << "[INFO]>>>>>>>>>>>>>>>>>>>iter "<< i <<" <<<<<<<<<<<<<<<<<<" << endl;
    // cout << "displace is " << endl<< displace.block(0, 0, 3, 8) << endl;

    //newtown iter
    Map<VectorXd> displace_plus(displace.data(), 3*dim);
    for(size_t newton_i = 0; newton_i < 20; ++newton_i){
      cout << "newton iter is " << newton_i << endl;

      dat_str->set_zero();
      // dat_str->hes_reserve(nnzs);
      energy->Val(displace_plus.data(), dat_str);
      energy->Gra(displace_plus.data(), dat_str);
      energy->Hes(displace_plus.data(), dat_str);
      dat_str->hes_compress();

      const double res_value = res.array().square().sum();
      cout << "[INFO]Newton res " <<std::setprecision(9)<< res_value << endl;
      cout << "[INFO] ALL Energy: " << dat_str->get_val() << endl;

      if(res_value < 1e-4){
        cout << endl;
        break;
      }
      
      
      //implicit time integral

             
      // cout << "[INFO]>>>>>>>>>>>>>>>>>>>LLT<<<<<<<<<<<<<<<<<<" << endl;
      auto start = system_clock::now();
      SimplicialLLT<SparseMatrix<double>> llt;
      llt.compute(dat_str->get_hes());
      size_t time = 1;
      while(llt.info() != Eigen::Success){
        cout <<"lltinfo "<< llt.info() << endl;
        dat_str->hes_add_diag(time);
        llt.compute(dat_str->get_hes());
        time *= 2;
      }

      solution = llt.solve(-res);
      auto end = system_clock::now();
      auto duration = duration_cast<microseconds>(end - start);
      cout <<  "solve linear system花费了" 
           << double(duration.count()) * microseconds::period::num / microseconds::period::den 
           << "秒" << endl;



      #if 0
      {//Line search

        const double c = 1e-4, c2 = 0.9;

        double Val_init = dat_str.val_, down = res.dot(solution), Val_upbound, Val_func, down_new =0;
        // cout <<endl << endl<<"Val init is " << Val_init << "down is "<< down << endl;

        auto cal_val_gra = [&](const double alp){
          displace_search = displace_plus + alp * solution;
          dat_str.set_zero();
          MO.Val(displace_search.data(), dat_str);
          PS.Val(displace_search.data(), dat_str);
          GE.Val(displace_search.data(), dat_str);
          pos_cons.Val(displace_search.data(), dat_str);
          
          MO.Gra(displace_search.data(), dat_str);
          PS.Gra(displace_search.data(), dat_str);
          GE.Gra(displace_search.data(), dat_str);
          pos_cons.Gra(displace_search.data(), dat_str);
        };
        auto cal_val = [&](const double alp)->double{
          displace_search = displace_plus + alp * solution;
          dat_str.set_zero();
          MO.Val(displace_search.data(), dat_str);
          PS.Val(displace_search.data(), dat_str);
          GE.Val(displace_search.data(), dat_str);
          pos_cons.Val(displace_search.data(), dat_str);
          return dat_str.val_;
        };
        auto cal_gra = [&](){//use little
          MO.Gra(displace_search.data(), dat_str);
          PS.Gra(displace_search.data(), dat_str);
          GE.Gra(displace_search.data(), dat_str);
          pos_cons.Gra(displace_search.data(), dat_str);
        };

        auto zoom = [&](double alpha_low, double alpha_high, double val_low)->double{
          double alpha_star = alpha_high;
          size_t count_j = 1;
          do{
            double alpha_j = 0.5 * (alpha_low + alpha_high);
            cout << "alpha j is "<< alpha_j << endl;
            double val_j = cal_val(alpha_j);
            if(val_j > Val_init + c * alpha_j * down || val_j > val_low){
              alpha_high = alpha_j;
              cout << " here " << endl;
            }
              
            else{
              cal_gra();
              double deri = solution.dot(dat_str.gra_);
              if(fabs(deri) <= -c2 * down){
                alpha_star = alpha_j;
                break;
              }
              if(deri * (alpha_high - alpha_low) >= 0)
                alpha_high = alpha_low;
              alpha_low = alpha_j;
              val_low = val_j;
            }
            ++count_j;
          }while(count_j < 10);
          return alpha_star;
        };
        double val_now, val_before = Val_init, alpha_now = 1, alpha_before = 0, alpha_fin, alpha_max = 2;
        size_t count = 1;
        do{
          cout << "alpha now is "<< alpha_now << endl;
          val_now = cal_val(alpha_now);
          if(val_now > Val_init + c * alpha_now * down || (val_now >= val_before && count > 1)){
            alpha_fin = zoom(alpha_before, alpha_now, val_before);
            break;
          }
          cal_gra();
          double deri = solution.dot(dat_str.gra_);          
          if(fabs(deri) <= -c2 * down){
            alpha_fin = alpha_now;
            break;
          }
          if(deri >= 0){
            alpha_fin = zoom(alpha_now, alpha_before, val_now);
            break;
          }
          alpha_before = alpha_now;
          alpha_now = 0.5 * (alpha_now + alpha_max);

          val_before = val_now;
          ++count;
        }while(1);
        cout<< "line search alpha is "<<  alpha_fin<<endl;        
        displace_plus += alpha_fin * solution;

      }

      #else
      displace_plus += solution;
      #endif

      cout << endl;
    }
    
    dynamic_pointer_cast<momentum<3>>(ebf[MOME])->update_location_and_velocity(displace_plus.data());


    auto surf_filename = outdir  + "/" + mesh_name + "_" + to_string(i) + ".obj";
    auto point_filename = outdir + "/" + mesh_name + "_points_" + to_string(i) + ".vtk";

    points_now = points + displace;
    vet_displace = displace.block(0, 0, 3, nods.cols());
    point_write_to_vtk(point_filename.c_str(), points_now.data(), dim);
    writeOBJ(surf_filename.c_str(), (nods + vet_displace).transpose(), surf.transpose());
    


  }
  //done

  return 0;
}

  




