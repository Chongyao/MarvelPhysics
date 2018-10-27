#ifndef PARAM_UTIL_H
#define PARAM_UTIL_H

#include <zjucad/matrix/matrix.h>
#include <zjucad/matrix/itr_matrix.h>
#include <Eigen/Eigen>
#include <Eigen/Sparse>
#include <unordered_set>
#include <hjlib/math/blas_lapack.h>
#include <zjucad/matrix/lapack.h>
#include <set>
#include <numeric>
#include "config.h"

namespace bigbang {

template <typename T, size_t dim = 3>
void add_diag_block(const size_t row, const size_t col, const T val, std::vector<Eigen::Triplet<T>> *mat) {
  const size_t row_start = row * dim;
  const size_t col_start = col * dim;
  for (size_t offset = 0; offset < dim; ++offset)
    mat->push_back(Eigen::Triplet<T>(row_start+offset, col_start+offset, val));
}

template <typename T>
void runtime_dim_add_diag_block(const size_t dim, const size_t row, const size_t col, const T val,
                                std::vector<Eigen::Triplet<T>> *trip) {
  const size_t row_start = row * dim;
  const size_t col_start = col * dim;
  for (size_t off = 0; off < dim; ++off)
    trip->push_back(Eigen::Triplet<T>(row_start+off, col_start+off, val));
}

template <typename T>
void insert_block(const size_t start_x, const size_t start_y, const T *data,
                  const size_t rows, const size_t cols,
                  std::vector<Eigen::Triplet<T>> *trips) {
  // inserted block is column major
  for (size_t j = 0; j < cols; ++j) {
    for (size_t i = 0; i < rows; ++i) {
      size_t idx = i+j*rows;
      trips->push_back(Eigen::Triplet<T>(start_x+i, start_y+j, data[idx]));
    }
  }
}

template <typename T>
void rm_spmat_col_row(Eigen::SparseMatrix<T> &A,
                      const std::unordered_set<size_t> &idx) {
  std::vector<size_t> g2l(A.cols());
  size_t ptr = 0;
  for (size_t i = 0; i < g2l.size(); ++i) {
    if ( idx.find(i) != idx.end() )
      g2l[i] = -1;
    else
      g2l[i] = ptr++;
  }
  rm_spmat_col_row<T>(A, g2l);
}

template <typename T, class Con>
void rm_spmat_col_row(Eigen::SparseMatrix<T> &A,
                      const Con &g2l) {
  size_t new_size = 0;
  for (size_t i = 0; i < g2l.size(); ++i) {
    if ( g2l[i] != -1)
      ++new_size;
  }
  std::vector<Eigen::Triplet<T>> trips;
  for (size_t j = 0; j < A.outerSize(); ++j) {
    for (typename Eigen::SparseMatrix<T>::InnerIterator it(A, j); it; ++it) {
      if ( g2l[it.row()] != -1 && g2l[it.col()] != -1 )
        trips.push_back(Eigen::Triplet<T>(g2l[it.row()], g2l[it.col()], it.value()));
    }
  }
  A.resize(new_size, new_size);
  A.reserve(trips.size());
  A.setFromTriplets(trips.begin(), trips.end());
}

template <typename T, class Con>
void rm_vector_row(Eigen::Matrix<T, -1, 1> &b,
                   const Con &g2l) {
  size_t new_size = 0;
  for (size_t i = 0; i < g2l.size(); ++i) {
    if ( g2l[i] != -1 )
      ++new_size;
  }
  Eigen::Matrix<T, -1, 1> sub;
  sub.resize(new_size);
#pragma omp parallel for
  for (size_t i = 0; i < g2l.size(); ++i)
    if ( g2l[i] != -1 )
      sub[g2l[i]] = b[i];
  b = sub;
}

template <typename T, class Con>
void rc_vector_row(const Eigen::Matrix<T, -1, 1> &l, const Con &g2l, Eigen::Matrix<T, -1, 1> &g) {
#pragma omp parallel for
  for (size_t i = 0; i < g2l.size(); ++i) {
    if ( g2l[i] != -1 )
      g[i] = l[g2l[i]];
  }
}

template <class Vec, class Con>
void rc_vector_row(const Vec &l, const Con &g2l, Vec &g) {
  #pragma omp parallel for
  for (size_t i = 0; i < g2l.size(); ++i) {
    if ( g2l[i] != -1 )
      g[i] = l[g2l[i]];
  }
}

template <typename T>
Eigen::SparseMatrix<T> sparse_diag_matrix(const Eigen::DiagonalMatrix<T, -1> &diag) {
  std::vector<Eigen::Triplet<T>> trips;
  for (size_t i = 0; i < diag.cols(); ++i)
    trips.push_back(Eigen::Triplet<T>(i, i, diag.diagonal()[i]));
  Eigen::SparseMatrix<T> rtn;
  rtn.resize(diag.rows(), diag.cols());
  rtn.reserve(trips.size());
  rtn.setFromTriplets(trips.begin(), trips.end());
  return rtn;
}

template <typename T>
bool is_symm(const Eigen::SparseMatrix<T> &A) {
  Eigen::SparseMatrix<T> AT = A.transpose();
  if ( (AT - A).squaredNorm() < 1e-16 )
    return true;
  return false;
}

template <typename T>
inline T cal_cot_val(const T* a, const T* b, const T* c) {
  Eigen::Matrix<T, 3, 1> ab(b[0]-a[0], b[1]-a[1], b[2]-a[2]);
  Eigen::Matrix<T, 3, 1> bc(c[0]-b[0], c[1]-b[1], c[2]-b[2]);
  Eigen::Matrix<T, 3, 1> ca(a[0]-c[0], a[1]-c[1], a[2]-c[2]);
  return 0.5 * (ab.dot(ab) + bc.dot(bc) - ca.dot(ca)) / ab.cross(bc).norm();
}

template <typename T>
inline T cal_cot_val(const T *point_set) {
  return cal_cot_val<T>(&point_set[0], &point_set[3], &point_set[6]);
}

template <typename T>
inline T safe_acos(const T cosval) {
  return std::acos(std::min(1.0, std::max(-1.0, cosval)));
}

template <typename T>
inline T safe_sqrt(const T val) {
  return std::sqrt(std::max(static_cast<T>(0), val));
}

template <typename INT=size_t>
INT build_global_local_mapping(const INT dim, const std::unordered_set<INT> &fixDOF, std::vector<INT> &g2l) {
  g2l.resize(dim);
  INT ptr = static_cast<INT>(0);
  for (INT i = 0; i < dim; ++i) {
    if ( fixDOF.find(i) != fixDOF.end() )
      g2l[i] = static_cast<INT>(-1);
    else
      g2l[i] = ptr++;
  }
  return ptr;
}

template <typename T, int Option>
int extract_triplets_from_spmat(const Eigen::SparseMatrix<T, Option> &A, std::vector<Eigen::Triplet<T>> *trips) {
  if ( !trips )
    return __LINE__;
  for (size_t j = 0; j < A.outerSize(); ++j) {
    for (typename Eigen::SparseMatrix<T, Option>::InnerIterator it(A, j); it; ++it) {
      trips->push_back(Eigen::Triplet<T>(it.row(), it.col(), it.value()));
    }
  }
  return 0;
}

template <typename T>
void extract_rotation(const T *df, T *R) {
  zjucad::matrix::matrix<T> F = zjucad::matrix::itr_matrix<const T *>(3, 3, df);
  zjucad::matrix::matrix<T> U(3, 3), S(3, 3), VT(3, 3);
  svd(F, U, S, VT);
  zjucad::matrix::itr_matrix<T*>(3, 3, R) = U*VT;
  // Eigen::Matrix<T, 3, 3> J(df);
  // Eigen::JacobiSVD<Eigen::Matrix<T, 3, 3>> svd(J, Eigen::ComputeFullU|Eigen::ComputeFullV);
  // Eigen::Map<Eigen::Matrix<T, 3, 3>>(R, 3, 3) = svd.matrixU()*svd.matrixV().transpose();
}

template <class Vec, class Mat>
Vec flatten(const Mat &A) {
  Vec v;
  v.resize(A.size());
  std::copy(&A(0, 0), &A(0, 0)+A.size(), &v[0]);
  return v;
}

template <typename T, int dim>
zjucad::matrix::matrix<T>
kroneckerId(const zjucad::matrix::matrix<T> &A) {
  using zjucad::matrix::colon;
  using zjucad::matrix::eye;
  using zjucad::matrix::zeros;
  zjucad::matrix::matrix<T> tmp = zeros<T>(dim*A.size(1), dim*A.size(2));
  for (size_t i = 0; i < A.size(1); ++i) {
    for (size_t j = 0; j < A.size(2); ++j) {
      tmp(colon(dim*i, dim*i+dim-1), colon(dim*j, dim*j+dim-1)) = A(i, j)*eye<T>(dim);
    }
  }
  return tmp;
}

template <class Con1, class Con2>
int sort_with_indices(const Con1 &vals, Con2 &idx) {
  typedef typename Con1::value_type val_type;
  typedef typename Con2::value_type int_type;
  idx.resize(vals.size());
  std::iota(idx.begin(), idx.end(), static_cast<int_type>(0));
  std::sort(idx.begin(), idx.end(), [&vals](int_type i, int_type j) { return vals[i] < vals[j]; });
  return 0;
}

// A: symmetric matrix
// B: diagonal matrix
// option: 0 for smallest, 1 for largest
int solve_gen_eig_prob(const Eigen::SparseMatrix<double> &K,
                       const Eigen::SparseMatrix<double> &M,
                       const size_t eignum, const int options,
                       const std::string &solver,
                       Eigen::VectorXd *eigenvalues,
                       Eigen::MatrixXd *eigenvectors);

// A: symmetric matrix
// option: 0 for smallest, 1 for largest
int solve_sym_eig_prob(const Eigen::SparseMatrix<double> &K,
                       const size_t eignum, const int options,
                       const std::string &solver,
                       Eigen::VectorXd *eigenvalues,
                       Eigen::MatrixXd *eigenvectors);
  
template <class Mat>
void normalize_vec_metric(Eigen::VectorXd &x, const Mat &M) {
  assert(M.rows() == M.cols() && M.cols() == x.size());
  const double normM = sqrt(x.dot(M*x));
  x /= normM;
}

template <class Con1, class Con2, class Con3>
void find_intersections(const Con1 &a, const Con2 &b, Con3 &intersection) {
  typedef typename Con1::value_type value_type_1;
  typedef typename Con2::value_type value_type_2;
  typedef typename Con3::value_type value_type_3;
  
  std::set<value_type_1> set_a(a.begin(), a.end());
  std::set<value_type_2> set_b(b.begin(), b.end());

  std::vector<value_type_3> res(a.size()+b.size());
  const auto it = std::set_intersection(set_a.begin(), set_a.end(), set_b.begin(), set_b.end(), res.begin());

  intersection.resize(it-res.begin());
  std::copy(res.begin(), it, intersection.begin());
}

template <class Vec, class Mat>
struct selection_matrix_builder
{
  static int build(const Vec &idx, const size_t dim, const size_t rd, Mat &S) {
    for (auto &i : idx)
      if ( i >= dim )
        return __LINE__;

    S.resize(rd*idx.size(), rd*dim);
    std::fill(&S(0, 0), &S(0, 0)+S.size(), 0);
    for (size_t i = 0; i < idx.size(); ++i) {
      for (size_t j = 0; j < rd; ++j) {
        S(rd*i+j, rd*idx[i]+j) = 1;
      }
    }
    return 0;
  }
};

template <class Vec>
struct selection_matrix_builder<Vec, Eigen::SparseMatrix<double>>
{
  static int build(const Vec &idx, const size_t dim, const size_t rd, Eigen::SparseMatrix<double> &S) {
    for (auto &i : idx)
      if ( i >= dim )
        return __LINE__;

    S.resize(rd*idx.size(), rd*dim);
    S.setZero();
    for (size_t i = 0; i < idx.size(); ++i) {
      for (size_t j = 0; j < rd; ++j) {
        S.insert(rd*i+j, rd*idx[i]+j) = 1;
      }
    }
    return 0;
  }
};

template <class Vec, class Mat>
int build_selection_matrix(const Vec &idx, const size_t dim, const size_t rd, Mat &S) {
  return selection_matrix_builder<Vec, Mat>::build(idx, dim, rd, S);
}

template <typename T>
void block_diagonalize_matrix(const zjucad::matrix::matrix<T> &mat, const size_t rd,
                              Eigen::SparseMatrix<T> &diagM) {
  const size_t rows = rd*mat.size(1), cols = rd*mat.size(2);
  diagM.resize(rows, cols);
  
  for (size_t k = 0; k < rd; ++k) {
    for (size_t i = 0; i < mat.size(1); ++i) {
      for (size_t j = 0; j < mat.size(2); ++j) {
        diagM.insert(k*mat.size(1)+i, k*mat.size(2)+j) = mat(i, j);
      }
    }
  }
}

}
#endif
