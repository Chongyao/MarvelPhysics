/** -*- mode: c++ -*-
 * @file mprgp_utils.h
 * @author LamKamhang (Cool_Lam@outlook.com)
 * @brief A Documented file.
 * @version 1.0
 * @date Mon Nov 18 10:56:07 CST 2019
 *
 * Detailed description
 *
 *
 * @copyright Copyright (c) 2019
 */
#pragma once

#include "utils/logger/log_utils.h"
#include "utils/timer/timer_utils.h"
#include "utils/logger/assert_utils.h"

#include <Eigen/Sparse>
#include <Eigen/Dense>

#include <vector>

namespace chaos
{
  namespace mprgp
  {
    template<typename T> inline constexpr T scalar_max() {return std::numeric_limits<T>::max() / 10;}
    template<typename T> inline constexpr T scalar_eps() {return 1;}

    template<> inline constexpr double scalar_eps() {return 1e-9f;}
    template<> inline constexpr float scalar_eps() {return 1e-5f;}

    template<typename Vec>
    inline void mask_face(const Vec &in, Vec &out, const std::vector<char> &face)
    {
      assert_ext((size_t) in.size () == (size_t) face.size (),
                 "in mask face: check in.size == face.size");
      out.resize (in.size ());
      assert_ext ((size_t) in.size () == (size_t) out.size (),
                  "in mask face: check in.size == out.size");
      for (size_t i = 0; i < face.size (); ++i)
        face [i] == 0 ? out [i] = in [i] : out [i] = 0.0f;
    }

    template<typename Vec>
    inline int count_constraints (const Vec &face)
    {
      int cnt = 0;

      for (size_t i = 0; i < face.size (); ++i)
        {
          assert_ext ((int)face [i] >= 0,
                      "in count constraints: check face[i] >= 0");
          cnt += (int)face [i];
        }

      return cnt;
    }

    template<class T>
    Eigen::SparseMatrix<T> Dense2Sparse(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &Dense,
                                        const T tol = scalar_eps<T>())
    {
      using Tri = Eigen::Triplet<T>;

      std::vector<Tri> striplet;
      for (int i = 0; i < Dense.rows(); ++i)
        for (int j = 0; j < Dense.cols(); ++j)
          if (std::abs(Dense(i, j)) >= tol)
            striplet.push_back(Tri(i, j, Dense(i, j)));

      Eigen::SparseMatrix<T> Sparse(Dense.rows(), Dense.cols());
      Sparse.setFromTriplets(striplet.begin(), striplet.end());
      return Sparse;
    }

    template<typename T, typename Vec>
    inline void get_diagonal (const Eigen::SparseMatrix<T> &A, Vec &diag_A)
    {
      diag_A.setConstant(std::min(A.rows(), A.cols()), 0);
      assert_ext ((size_t) diag_A.size () == (size_t) A.rows () ||
                  (size_t) diag_A.size () == (size_t) A.cols (),
                  "in get diagonal: check diag_A size");

      for (int i = 0; i < A.outerSize (); ++i)
        for (typename Eigen::SparseMatrix<T>::InnerIterator it (A, i); it; ++it)
          {
            if (it.col () == it.row ())
              {
                assert_ext (it.row () == i && it.col () == i,
                            "in get_diagonal: check item col equal to row");
                diag_A [i] = it.value ();
                break;
              }
          }
    }

    template<typename T, typename Vec>
    inline void get_diagonal (const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &A, Vec &diag_A)
    {
      return get_diagonal(Dense2Sparse(A), diag_A);
    }
  }  // mprgp
}  // chaos
