#include <string>
#include <iostream>
#include <chrono>

#include <memory>

#include <libigl/include/igl/readOBJ.h>
#include <libigl/include/igl/writeOBJ.h>
#include <Eigen/Core>
#include <Eigen/Geometry>
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
#include "Point_Sys/src/basic_energy.cc"

#include "vtk2surf.h"


#include "coll_response.h"
#include <Collision/CollisionDetect-rigid/src/Collision_eigen.h>


using namespace marvel;
using namespace std;
using namespace Eigen;
using namespace igl;
using namespace chrono;
using namespace boost;


Matrix3d get_tri_pos(const MatrixXi& tris, const MatrixXd& verts, const size_t& face_id){
  Matrix3d tri;
  for(size_t i = 0; i < 3; ++i){
    size_t vert_id = tris(i, face_id);
    tri.col(i) = verts.col(vert_id);
  }
  return std::move(tri);
}

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
  
  auto common = pt.get_child("common");
  auto blender = pt.get_child("blender");
  auto physics_para = pt.get_child("physics_para");
  auto simulation_para = pt.get_child("simulation_para");
  
  cout << "[INFO]>>>>>>>>>>>>>>>>>>>IMPORT MESH<<<<<<<<<<<<<<<<<<" << endl;
  const string mesh_name = blender.get<string>("surf");
  const string indir = "../input";
  const string outdir = "../output/" + mesh_name;
  //mkdir
  boost::filesystem::path outpath(outdir);
  if ( !boost::filesystem::exists(outdir) )
    boost::filesystem::create_directories(outdir);

  MatrixXi surf;
  MatrixXd nods;
  readOBJ((indir + '/' +mesh_name+".obj").c_str(), nods, surf);
  cout << "surf: " << surf.rows() << " " << surf.cols() << endl << "nods: " << nods.rows() << " " << nods.cols() << endl;
  
  surf.transposeInPlace();
  nods.transposeInPlace();
  

  cout << "[INFO]>>>>>>>>>>>>>>>>>>>Generate sampled points<<<<<<<<<<<<<<<<<<" << endl;
  MatrixXd points(3,3);
  MatrixXd test(3, 3);
  gen_points(nods, surf, simulation_para.get<size_t>("num_in_axis"), points, true);
  cout << points.rows() << " " << points.cols() << endl;
  // #if 1
  // points = nods;
  // #endif
  size_t dim = points.cols();
  cout <<"generate points done." << endl;
  
  cout << "[INFO]>>>>>>>>>>>>>>>>>>>Build spatial hash<<<<<<<<<<<<<<<<<<" << endl;
  spatial_hash SH(points, simulation_para.get<size_t>("nn_num"));

  
  cout << "[INFO]>>>>>>>>>>>>>>>>>>>Build Point System<<<<<<<<<<<<<<<<<<" << endl;
  //calc volume 
  const double volume = clo_surf_vol(nods, surf);
  //calc support radii
  VectorXd sup_radi = SH.get_sup_radi();
  
  cout << "[INFO]>>>>>>>>>>>>>>>>>>>sup_radi<<<<<<<<<<<<<<<<<<" << endl;

  //get friends of every point
  vector<vector<size_t>> friends_all(dim);
  // #pragma omp parallel for
  for(size_t i = 0; i < dim; ++i){
    SH.get_friends(points.col(i), sup_radi(i), friends_all[i]);
  }



  
  point_sys PS(points, common.get<double>("density"), physics_para.get<double>("Young"), physics_para.get<double>("Poission"), volume, simulation_para.get<double>("kv"), friends_all, sup_radi);


  cout << "[INFO]>>>>>>>>>>>>>>>>>>>Simple Constraint Points<<<<<<<<<<<<<<<<<<" << endl;
  //add simple constraints
  //This should read from file. We loop for some points to restrain here.
  //Constraints vary from different models and situations.
  vector<size_t> cons(0);
  auto cons_file_path = indir + "/" +  mesh_name +".csv";

  if ( boost::filesystem::exists(cons_file_path) )
    read_fixed_verts_from_csv(cons_file_path.c_str(), cons);
  cout << "constrint " << cons.size() << " points" << endl;
  #if 1
  for(auto con : cons){
    cout << con << " ";
  }
  cout << endl;
  #endif

  position_constraint pos_cons(simulation_para.get<double>("position_weig"), cons, dim);

  cout << "[INFO]>>>>>>>>>>>>>>>>>>>Gravity<<<<<<<<<<<<<<<<<<" << endl;
  const double gravity = common.get<double>("gravity");
  gravity_energy GE(simulation_para.get<double>("w_g"), gravity, dim, PS.get_Mass_VectorXd(), 'y');

  
  cout << "[INFO]>>>>>>>>>>>>>>>>>>>COLLISION<<<<<<<<<<<<<<<<<<" << endl;
  // collision COLL(simulation_para.get<double>("w_coll"),'y', simulation_para.get<double>("g_pos"), nods.cols(), dim);
  vector<std::shared_ptr<MatrixXi> > obta_surfs;
  vector<std::shared_ptr<MatrixXd>> obta_nods;{
    filesystem::path obstacles_path(indir + "/" + "obstacles");
    filesystem::directory_iterator end_iter;
    for(filesystem::directory_iterator file_iter(obstacles_path); file_iter != end_iter; ++file_iter){
      if(filesystem::is_regular_file(file_iter->path())){
        string one_obstacle = file_iter->path().string();
        MatrixXd one_obta_nods;
        MatrixXi one_obta_surf;
        readOBJ(one_obstacle.c_str(), one_obta_nods, one_obta_surf);
        one_obta_surf.transposeInPlace();
        one_obta_nods.transposeInPlace();
        obta_surfs.push_back(std::move(make_shared<MatrixXi>(one_obta_surf)));
        obta_nods.push_back(std::move(make_shared<MatrixXd>(one_obta_nods)));
        cout << one_obstacle << endl;
      }
    }    
  }

  size_t obta_num = obta_surfs.size();
  // vector<vector<double>> obta_areas(obta_num);{
  //   for(auto& obta : obta_areas){
  //     obta.resize(obta_surfs[i]->cols());
      
  //   }

  // }
  
  size_t num_fake_tris = dim%3 ? dim / 3 + 1: dim / 3 ;
  MatrixXi fake_surf(3, num_fake_tris);{
    #pragma omp parallel for
    for(size_t i = 0; i < num_fake_tris; ++i){
      fake_surf(0, i) = i * 3;
      fake_surf(1, i) = i * 3 + 1 > dim ? i * 3 - 1 : i * 3 + 1;
      fake_surf(2, i) = i * 3 + 2 > dim ? i * 3 - 2 : i * 3 + 2;
    }
  }
  
  
  auto COLL_ptr = Collision_zcy::getInstance();
  COLL_ptr->Transform_Mesh(dim, num_fake_tris, fake_surf.data(), points.data(), points.data(), 0, false);
  for(size_t i = 0; i < obta_num; ++i){
    COLL_ptr->Transform_Mesh(obta_nods[i]->cols(), obta_surfs[i]->cols(), obta_surfs[i]->data(), obta_nods[i]->data(), obta_nods[i]->data(), i + 1, false);
    COLL_ptr->Transform_Pair(0, i + 1);
  }

  COLL_ptr->Collid();

  

             


  
  cout << "[INFO]>>>>>>>>>>>>>>>>>>>SOlVE<<<<<<<<<<<<<<<<<<" << endl;
  //initilize variables in time integration
  energy_dat dat_str (dim);

  string solver = simulation_para.get<string>("solver");

  double delt_t = common.get<double>("time_step");
  MatrixXd points_pos, new_pos;
  MatrixXd displace, new_displace;
  MatrixXd velocity, new_velocity;
  MatrixXd acce;
  MatrixXd new_acce;
  MatrixXd gra;
  MatrixXd vet_displace;
  points_pos.setZero(3, dim); new_pos.setZero(3, dim);
  displace.setZero(3, dim); new_displace.setZero(3, dim);
  velocity.setZero(3, dim); new_velocity.setZero(3, dim);
  acce.setZero(3, dim); new_acce.setZero(3, dim);
  gra.setZero(3, dim);

  vet_displace.setZero(3, nods.cols());

  PS.pre_compute(dat_str);
  size_t iters_perframe = static_cast<size_t>(round(1.0/delt_t/common.get<size_t>("frame_rate")));
  size_t max_iter  = static_cast<size_t>(ceil(common.get<double>("total_time") / delt_t));
  cout << "max iter is " << max_iter << endl;
  double dump = simulation_para.get<double>("dump");
  double previous_step_Val = 0;

  auto start = system_clock::now();
  size_t frame_id = 0;
  if(solver == "explicit"){
    for(size_t i = 0; i < max_iter; ++i){
      cerr << "iter is "<<endl<< i << endl;
      cout << "displace is " << endl<< displace.block(0, 0, 3, 7) << endl;
      cout << "velocity is "<<endl<< velocity.block(0, 0, 3, 7) << endl;
      // cout << "acce is " << endl << acce.block(0, 0, 3, 8) << endl;

      points_pos = points + displace;
      GE.Val(displace.data(), dat_str);
      GE.Gra(displace.data(), dat_str);

      PS.Val(displace.data(), dat_str);
      PS.Gra(displace.data(), dat_str);

      // COLL.Val(points.data(), displace.data(), dat_str);
      // COLL.Gra(points.data(), displace.data(), dat_str, PS.get_Mass_VectorXd());
      pos_cons.Gra(displace.data(), dat_str);
      pos_cons.Hes(displace.data(),dat_str);
         

      


      for(size_t j = 0; j < dim; ++j){
        assert(PS.get_mass(j) > 0);
        new_acce.col(j) = dat_str.gra_.col(j)/PS.get_mass(j) - velocity.col(j)*dump;
      }
    
      new_velocity += delt_t * new_acce;
      new_displace += delt_t *velocity;
      new_pos = points + new_displace;



      //>>>>>>>>>>>>>>>>>>COLLID<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<//
      COLL_ptr->Transform_Mesh(dim, num_fake_tris, fake_surf.data(), new_pos.data(), points_pos.data(), 0, false);
      COLL_ptr->Collid();

      auto pairs = COLL_ptr->getContactPairs();
      auto times = COLL_ptr->getContactTimes();
      // assert(pairs.size() == 0);
      for(size_t j = 0; j < pairs.size(); ++j){

        uint mesh_id1, face_id1, mesh_id2, face_id2;{
          cout << "j is " << j << endl;
          pairs[j][0].get(mesh_id1, face_id1);
          pairs[j][1].get(mesh_id2, face_id2);
          cout << mesh_id1 << " " << mesh_id2 << " " << face_id1 << " " << face_id2;
          if(mesh_id2 == 0){
            mesh_id2 = mesh_id1;
            mesh_id1 = 0;

            auto exchange = face_id2;
            face_id2 = face_id1;
            face_id1 = exchange;
          }
        }
        if(mesh_id2 == 0)
          continue;

        //TODO: can be faster
        auto coll_plane = get_tri_pos(*(obta_surfs[mesh_id2 - 1]), *(obta_nods[mesh_id2 - 1]), face_id2);
        auto pre_pos = get_tri_pos(fake_surf, points_pos, face_id1);
        auto pre_velo = get_tri_pos(fake_surf, velocity, face_id1);
        auto next_pos = get_tri_pos(fake_surf, new_pos, face_id1);
        auto next_velo = get_tri_pos(fake_surf, new_velocity, face_id1);

        cout << "before response" << endl << "pre pos : " <<endl << pre_pos << endl << "after_pos :" <<endl<< next_pos << "pre velo :" << endl << pre_velo << endl << "after pos : " << endl << next_velo << endl;
        Matrix3d res_pos = Matrix3d::Zero();
        Matrix3d res_velo = Matrix3d::Zero();

        response(coll_plane.data(), times[j], pre_pos.data(), next_pos.data(),
                 pre_velo.data(), next_velo.data(),
                 res_pos.data(), res_velo.data());

        cout << "after response" << res_pos << endl << endl << res_velo << endl;
        for(size_t k = 0; k < 3; ++k){
          size_t vert_id = fake_surf(k, face_id1);
          new_velocity.col(vert_id) = res_velo.col(k);
          new_displace.col(vert_id) = res_pos.col(k) - points.col(k);
          new_pos.col(vert_id) = res_pos.col(k);
        }
      }//TODO:make it a new class
      //>>>>>>>>>>>>>>>>>>COLLID<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<//
      
      if (pairs.size() != 0)
        return 0;
      
      
      if (i > 10 && fabs(dat_str.Val_ - previous_step_Val) < 1e-6)
        break;
    
      previous_step_Val = dat_str.Val_;
      
      acce = new_acce;
      velocity = new_velocity;
      displace = new_displace;
      
      if(i%iters_perframe == 0){
        auto surf_filename = outdir  + "/" + mesh_name + "_" + to_string(frame_id) + ".obj";
        auto point_filename = outdir + "/" + mesh_name + "_points_" + to_string(frame_id) + ".vtk";
        MatrixXd points_now = points + displace;
        point_write_to_vtk(point_filename.c_str(), points_now.data(), dim);
        point_vector_append2vtk(false, point_filename.c_str(), velocity, dim, "velocity");
        point_vector_append2vtk(true, point_filename.c_str(), acce, dim, "accelarate");
        point_scalar_append2vtk(true, point_filename.c_str(), dat_str.ela_val_, dim, "strain_Energy");
        point_scalar_append2vtk(true, point_filename.c_str(), dat_str.vol_val_, dim, "vol_conservation_Energy");

        vet_displace = displace.block(0, 0, 3, nods.cols());
        writeOBJ(surf_filename.c_str(), (nods + vet_displace).transpose(), surf.transpose());
        ++frame_id;

      }

      dat_str.set_zero();
    }
  }
  else{
    SparseMatrix<double> A_CG(dim * 3, dim * 3);
    VectorXd b_CG(dim * 3);
    SparseMatrix<double> M = PS.get_Mass_Matrix();

    for(size_t i = 0; i < max_iter; ++i){
      cout << "iter is "<<endl<< i << endl;
      cout << "displace is " << endl<< displace.block(0, 0, 3, 8) << endl;
      // cout << "velocity is "<<endl<< velocity.block(0, 0, 3, 8) << endl;

      //newtown iter
      auto displace_plus = displace;
      Map<VectorXd> disp_t_plus(displace_plus.data(), 3*dim);
      Map<VectorXd> disp_t(displace.data(), 3*dim);
      Map<VectorXd> velo_t(velocity.data(), 3*dim);
      Map<VectorXd> _F(dat_str.gra_.data(), 3*dim);


      for(size_t newton_i = 0; newton_i < 999; ++newton_i){
        cout << "newton iter " << newton_i << endl;
      

    
        PS.Val(displace_plus.data(), dat_str);
        PS.Gra(displace_plus.data(), dat_str);
        PS.Hessian(displace_plus.data(), dat_str);

        GE.Val(displace_plus.data(), dat_str);
        GE.Gra(displace_plus.data(), dat_str);

        pos_cons.Gra(displace_plus.data(), dat_str);
        pos_cons.Hes(displace_plus.data(),dat_str);

        // COLL.Val(points.data(), displace_plus.data(), dat_str);
        // COLL.Gra(points.data(), displace_plus.data(), dat_str, PS.get_Mass_VectorXd());
        // COLL.Hes(displace_plus.data(), dat_str);
              

        //test  convergence
        auto res = M * ((disp_t_plus - disp_t) / delt_t - velo_t) - delt_t * _F;

        double res_value = res.array().square().sum();
        if(res_value < 1e-10){
          cout << "[INFO]Newton res " <<endl << res_value << endl;;
          cout << "[INFO]>>>>>>>>>>>>>>>>>>>Elasticity Energy Val<<<<<<<<<<<<<<<<<<" << endl;
          cout << dat_str.Val_ << endl;
          cout << "[INFO]>>>>>>>>>>>>>>>>>>>GRA<<<<<<<<<<<<<<<<<<" << endl;
          cout << dat_str.gra_.array().square().sum() << endl;
          cout << endl<<endl;
          break;
        }
      
      

    
        //implicit time integral
    
        A_CG.setZero();
        dat_str.hes_.setFromTriplets(dat_str.hes_trips.begin(), dat_str.hes_trips.end());
    
    
        A_CG = M + delt_t*delt_t*dat_str.hes_;
        b_CG = M * (delt_t * velo_t + disp_t - disp_t_plus) + delt_t * delt_t * _F;  
    
      
        cout << "[INFO]>>>>>>>>>>>>>>>>>>>A_CG<<<<<<<<<<<<<<<<<<" << endl;          ConjugateGradient<SparseMatrix<double>, Lower|Upper> cg;
        cg.setMaxIterations(3*dim);
        cg.setTolerance(1e-8);
        cg.compute(A_CG);
        disp_t_plus += cg.solve(b_CG);
        dat_str.set_zero();
        cout << "#iterations:     " << cg.iterations() << endl;
        cout << "estimated error: " << cg.error()      << endl;

      }
    
    
      velocity = (displace_plus - displace)/delt_t;
      displace = displace_plus;
      if(i%iters_perframe == 0){
        auto surf_filename = outdir  + "/" + mesh_name + "_" + to_string(frame_id) + ".obj";
        auto point_filename = outdir + "/" + mesh_name + "_points_" + to_string(frame_id) + ".vtk";
        MatrixXd points_now = points + displace;
        point_write_to_vtk(point_filename.c_str(), points_now.data(), dim);
        point_vector_append2vtk(false, point_filename.c_str(), velocity, dim, "velocity");
        point_vector_append2vtk(true, point_filename.c_str(), acce, dim, "accelarate");
        point_scalar_append2vtk(true, point_filename.c_str(), dat_str.ela_val_, dim, "strain_Energy");
        point_scalar_append2vtk(true, point_filename.c_str(), dat_str.vol_val_, dim, "vol_conservation_Energy");

        vet_displace = displace.block(0, 0, 3, nods.cols());
        writeOBJ(surf_filename.c_str(), (nods + vet_displace).transpose(), surf.transpose());
        ++frame_id;

      }
    }
  }
  auto end = system_clock::now();
  auto duration = duration_cast<microseconds>(end - start);
  cout <<  "花费了" 
       << double(duration.count()) * microseconds::period::num / microseconds::period::den 
       << "秒" << endl;

  //done
}





