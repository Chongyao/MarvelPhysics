#ifndef IO_H
#define IO_H

#include <fstream>
#include <iomanip>
#include <Eigen/Dense>
#include <Eigen/Sparse>

// using mati_t=zjucad::matrix::matrix<size_t>;
// using matd_t=zjucad::matrix::matrix<double>;

namespace marvel {

// int read_fixed_verts(const char *filename, std::vector<size_t> &fixed);

// int read_fixed_verts_from_csv(const char *filename, std::vector<size_t> &fixed, matd_t *pos=nullptr);
// int hex_mesh_read_from_vtk(const char *path, matd_t *node=nullptr, mati_t *hex=nullptr, matd_t *mtr=nullptr);

// int hex_mesh_write_to_vtk(const char *path, const matd_t &nods, const mati_t &hexs,
//                           const matd_t *mtr=nullptr, const char *type="CELL");

// int tet_mesh_write_to_vtk(const char *path, const matd_t &nods, const mati_t &tets, const matd_t *mtr=nullptr);

int tri_mesh_write_to_vtk(const char *path, const Eigen::MatrixXd &nods, const Eigen::MatrixXi &tris, const Eigen::MatrixXd *mtr=nullptr);
// int quad_mesh_write_to_vtk(const char *path, const matd_t &nods, const mati_t &quad,
//                            const matd_t *mtr=nullptr, const char *type="CELL");
int point_write_to_vtk(const char *path, const double *nods, const size_t num_points);
int vector_append_to_vtk(const bool is_append, const char* path, const Eigen::MatrixXd &vectors, const size_t num_vecs, const char* vector_name);
/*
int point_write_to_vtk(const char *path, const matd_t &nods);

int write_surf_vec_field(const char *path, const mati_t &surf, const matd_t &nods, const matd_t &vf);

template <typename T>
int read_dense_matrix(const char *path, zjucad::matrix::matrix<T> &dat) {
  std::ifstream ifs(path, std::ios::binary);
  if ( ifs.fail() )
    return __LINE__;
  int64_t mat_size[2];
  ifs.read((char *)&mat_size[0], 2*sizeof(int64_t));
  dat.resize(mat_size[0], mat_size[1]);
  ifs.read((char *)&dat[0], dat.size()*sizeof(T));
  ifs.close();
  return 0;
}

template <typename T>
int write_dense_matrix(const char *path, const zjucad::matrix::matrix<T> &dat ) {
  std::ofstream ofs(path, std::ios::binary);
  if ( ofs.fail() )
    return __LINE__;
  const int64_t mat_size[2] = {dat.size(1), dat.size(2)};
  ofs.write((const char *)&mat_size[0], 2*sizeof(int64_t));
  ofs.write((const char *)&dat[0], dat.size()*sizeof(T));
  ofs.close();
  return 0;
}

template <class Mat>
struct io_mat_trait {
  typedef typename Mat::value_type value_type;
};

template <typename T>
struct io_mat_trait<Eigen::Matrix<T, -1, -1>> {
  typedef T value_type;
};
  
template <typename T>
struct io_mat_trait<Eigen::Matrix<T, -1, 1>> {
  typedef T value_type;
};

template <class Mat>
int read_dense_matrix(const char *path, Mat &dat) {
  typedef typename io_mat_trait<Mat>::value_type T;
  std::ifstream ifs(path, std::ios::binary);
  if ( ifs.fail() )
    return __LINE__;
  int64_t mat_size[2];
  ifs.read((char *)&mat_size[0], 2*sizeof(int64_t));
  dat.resize(mat_size[0], mat_size[1]);
  ifs.read((char *)dat.data(), dat.size()*sizeof(T));
  ifs.close();
  return 0;
}

template <class Mat>
int write_dense_matrix(const char *path, const Mat &dat) {
  typedef typename io_mat_trait<Mat>::value_type T;
  std::ofstream ofs(path, std::ios::binary);
  if ( ofs.fail() )
    return __LINE__;
  const int64_t mat_size[2] = {dat.rows(), dat.cols()};
  ofs.write((const char *)&mat_size[0], 2*sizeof(int64_t));
  ofs.write((const char *)dat.data(), dat.size()*sizeof(T));
  ofs.close();
  return 0;
}

template <typename T>
int write_sparse_matrix(const char *path, const Eigen::SparseMatrix<T> &mat) {
  //-> only for col major sparse matrix
  Eigen::SparseMatrix<T> mat_tm = mat;
  if ( !mat_tm.isCompressed() )
    mat_tm.makeCompressed();

  std::ofstream ofs(path, std::ios::binary);
  if ( ofs.fail() )
    return __LINE__;
  
  const Eigen::Index rows = mat_tm.rows(), cols = mat_tm.cols(), nnz = mat_tm.nonZeros();
  ofs.write((const char *)&rows, sizeof(Eigen::Index));
  ofs.write((const char *)&cols, sizeof(Eigen::Index));
  ofs.write((const char *)&nnz,  sizeof(Eigen::Index));
  //-> the default SparseMatrix::StorageIndex is int
  ofs.write((const char *)mat_tm.outerIndexPtr(), (cols+1)*sizeof(int));
  ofs.write((const char *)mat_tm.innerIndexPtr(), nnz*sizeof(int));
  ofs.write((const char *)mat_tm.valuePtr(),      nnz*sizeof(T));
  ofs.close();

  return 0;
}

template <typename T>
int read_sparse_matrix(const char *path, Eigen::SparseMatrix<T> &mat) {
  std::ifstream ifs(path, std::ios::binary);
  if ( ifs.fail() )
    return __LINE__;

  Eigen::Index rows, cols, nnz;
  ifs.read((char *)&rows, sizeof(Eigen::Index));
  ifs.read((char *)&cols, sizeof(Eigen::Index));
  ifs.read((char *)&nnz,  sizeof(Eigen::Index));

  std::vector<int> ptr(cols+1), idx(nnz);
  std::vector<T> val(nnz);
  ifs.read((char *)&ptr[0], ptr.size()*sizeof(int));
  ifs.read((char *)&idx[0], idx.size()*sizeof(int));
  ifs.read((char *)&val[0], val.size()*sizeof(T));
  ifs.close();
  
  mat = Eigen::Map<Eigen::SparseMatrix<T>>(rows, cols, nnz, &ptr[0], &idx[0], &val[0]);

  return 0;
}

template <typename T>
int read_vector_vector(const char *path, std::vector<std::vector<T>> &object) {
  std::ifstream ifs(path);
  if ( ifs.fail() )
    return __LINE__;
  size_t object_num;  ifs >> object_num;
  object.resize(object_num);
  for (size_t i = 0; i < object_num; ++i) {
    size_t curr_num;  ifs >> curr_num;
    object[i].resize(curr_num);
    for (size_t j = 0; j < curr_num; ++j)
      ifs >> object[i][j];
  }
  return 0;
}

template <typename T>
int write_vector_vector(const char *path, const std::vector<std::vector<T>> &object) {
  std::ofstream ofs(path);
  if ( ofs.fail() )
    return __LINE__;
  ofs << object.size() << std::endl;
  for (size_t i = 0; i < object.size(); ++i) {
    ofs << object[i].size() << std::endl;
    for (size_t j = 0; j < object[i].size(); ++j)
      ofs << object[i][j] << " ";
    ofs << std::endl;
  }
  return 0;
}

template <typename T>
int write_idx_val(const char *path, const std::vector<T> &val,
                  const int64_t start, const int64_t end) {
  if ( 0 > start || start > end || end >= val.size() )
    return __LINE__;
  std::ofstream ofs(path);
  if ( ofs.fail() )
    return __LINE__;
  ofs << std::setprecision(12);
  ofs << "idx, value" << std::endl;
  for (size_t i = start; i <= end; ++i) {
    ofs << i << ", " << val[i] << std::endl;
  }
  ofs.close();
  return 0;
}

template <typename T1, typename T2>
int write_val_val(const char *path,
                  const std::vector<T1> &val1, const std::vector<T2> &val2,
                  const int64_t start, const int64_t end) {
  if ( 0 > start || start > end || end >= val1.size() )
    return __LINE__;
  if ( 0 > start || start > end || end >= val2.size() )
    return __LINE__;
  std::ofstream ofs(path);
  if ( ofs.fail() )
    return __LINE__;
  ofs << "value1, value2" << std::endl;
  ofs << std::setprecision(12);
  for (size_t i = start; i <= end; ++i) {
    ofs << val1[i] << ", " << val2[i] << std::endl;
  }
  ofs.close();
  return 0;
}

int read_json_file(const char *file, Json::Value &json);

class signed_dist_func;

int read_geom_objs_from_json(const char *file, std::vector<std::shared_ptr<signed_dist_func>> &objs);
*/
}

#endif
