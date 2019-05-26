#ifndef IO_H
#define IO_H

#include <fstream>
#include <iomanip>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>

// using mati_t=zjucad::matrix::matrix<size_t>;
// using matd_t=zjucad::matrix::matrix<double>;

namespace marvel {



int read_fixed_verts_from_csv(const char *filename, std::vector<size_t> &fixed, Eigen::MatrixXd *pos=nullptr);



int tri_mesh_write_to_vtk(const char *path, const Eigen::MatrixXd &nods, const Eigen::MatrixXi &tris, const Eigen::MatrixXd *mtr=nullptr);
// int quad_mesh_write_to_vtk(const char *path, const matd_t &nods, const mati_t &quad,
//                            const matd_t *mtr=nullptr, const char *type="CELL");
int point_write_to_vtk(const char *path, const double *nods, const size_t num_points);
int point_vector_append2vtk(const bool is_append, const char* path, const Eigen::MatrixXd &vectors, const size_t num_vecs, const char* vector_name);
int point_scalar_append2vtk(const bool is_append, const char* path, const Eigen::VectorXd &scalars, const size_t num_sca, const char* scalar_name);


int tet_mesh_read_from_vtk(const char* filename,  Eigen::MatrixXd& nods, Eigen::MatrixXi & tets,  double*   mtr= nullptr){
  std::ifstream ifs(filename);
  if(ifs.fail()) {
    std::cerr << "[info] " << "can not open file" << filename << std::endl;
    return __LINE__;
  }




  std::string str;
  int point_num = 0,cell_num = 0;
  
  while(!ifs.eof()){
    ifs >> str;
    if(str == "POINTS"){
      ifs >> point_num >> str;
      nods.derived().resize(3, point_num);
      double item;
      for(size_t i = 0;i < point_num; ++i){
        for(size_t j = 0;j < 3; ++j){
          ifs >> nods(j, i);
        }
          
      }
      continue;
    }
    if(str == "CELLS"){
      ifs >> cell_num >> str;
      int point_number_of_cell = 0;
      Eigen::Matrix<int, -1, -1> tet_temp(4, cell_num);
      size_t true_cell_num = 0;
      for(size_t ci = 0; ci < cell_num; ++ci){
        ifs >> point_number_of_cell;
        if(point_number_of_cell != 4){
          for(size_t i = 0; i < point_number_of_cell; ++i)
            ifs >> str;
        }else{
          int p;
          for(size_t i = 0; i < point_number_of_cell; ++i){
            ifs >> p;
            tet_temp(i, true_cell_num) = p;
          }
            
            
          ++true_cell_num;
        }
      }
      tets.resize(4, true_cell_num);
      tets = tet_temp.block(0, 0, 4, true_cell_num);
      break;
    }
  }





  std::vector<double> tmp;
  double mtrval;
  if ( mtr != nullptr ) {
    
    while ( !ifs.eof() ) {
      ifs >> str;
      if ( str == "LOOKUP_TABLE" ) {
        ifs >> str;
        for (size_t i = 0; i < tets.cols(); ++i) {
          ifs >> mtrval;
          tmp.push_back(mtrval);
        }
      }
    }

    if ( tmp.size() > 0 ) {
      assert(tmp.size() % tets.cols() == 0);
      Eigen::Map<Eigen::MatrixXd>mtr_mat(mtr, tets.cols(), tmp.size()/tets.cols());
      std::copy(tmp.begin(), tmp.end(), mtr_mat.data());
      mtr_mat.transposeInPlace();
    }
  }
  
  ifs.close();

  return 0;
}
// int tet_mesh_write_to_vtk(const char *path, const matd_t &nods, const mati_t &tets, const matd_t *mtr=nullptr);

}

#endif
