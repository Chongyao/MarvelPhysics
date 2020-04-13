/** -*- mode: c++ -*-
 * @file mprgp_precondition.h
 * @author LamKamhang (Cool_Lam@outlook.com)
 * @brief A Documented file.
 * @version 1.0
 * @date Sat Nov 16 16:39:53 CST 2019
 *
 * Detailed description
 *
 *
 * @copyright Copyright (c) 2019
 */
#pragma once

#include "utils/logger/log_utils.h"
#include "utils/timer/timer_utils.h"
#include "mprgp_utils.h"

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <memory>
#include <vector>

namespace chaos
{
  namespace mprgp
  {
    template <typename T, typename Mat>
    class MPRGPPreconditionBase
    {
      using Vec = Eigen::Matrix<T, Eigen::Dynamic, 1>;

    public:
      virtual ~MPRGPPreconditionBase() = default;

    public:
      virtual int solve(const Vec &g, Vec &z, const Vec &phi) const = 0;
    };

    // no precondition
    template <typename T, typename Mat>
    class InFaceNoPrecondition : public MPRGPPreconditionBase<T, Mat>
    {
      using Vec = Eigen::Matrix<T, Eigen::Dynamic, 1>;

    public:
      InFaceNoPrecondition(const std::vector<char> &face);

    public:
      int solve(const Vec &g, Vec &z, const Vec &phi) const;

    private:
      const std::vector<char> &face;
    };

    // Jacobian precondition
    template <typename T, typename Mat, bool Precond = true>
    class DiagonalInFacePrecondition : public MPRGPPreconditionBase<T, Mat>
    {
      using Vec = Eigen::Matrix<T, Eigen::Dynamic, 1>;

    public:
      DiagonalInFacePrecondition(const Mat &M, const std::vector<char> &face);

    public:
      int solve(const Vec &g, Vec &z, const Vec &phi) const;

    protected:
      Vec invert_diag;
      const std::vector<char> &face;
    };

    // decoupled constraints
    template <typename T, typename Mat, bool Precond = true>
    class DiagonalDecouplePrecondition : public DiagonalInFacePrecondition<T, Mat, Precond>
    {
      using Vec = Eigen::Matrix<T, Eigen::Dynamic, 1>;

    public:
      DiagonalDecouplePrecondition(const Mat &M,
                                   const std::vector<char> &face,
                                   const Eigen::SparseMatrix<T> &J);
      int solve(const Vec &g, Vec &z, const Vec &phi) const;

    private:
      Eigen::SparseMatrix<T> JMinv; // J*M^{-1}
      Vec invert_JMinvJt;
    };

    /////////////////////////////////////////////////////////////////////////////
    //                         template implementation                         //
    /////////////////////////////////////////////////////////////////////////////
    template <typename T, typename Mat>
    inline InFaceNoPrecondition<T, Mat>::InFaceNoPrecondition(const std::vector<char> &face)
      : face(face)
    {}

    template <typename T, typename Mat>
    inline int InFaceNoPrecondition<T, Mat>
    ::solve(const Vec &g, Vec &z, const Vec __attribute__((__unused__)) &phi) const
    {
      mask_face(g, z, face);
      return 0;
    }

    template <typename T, typename Mat, bool Precond>
    inline DiagonalInFacePrecondition<T, Mat, Precond>
    ::DiagonalInFacePrecondition(const Mat &M,
                                 const std::vector<char> &face)
      : face(face)
    {
      invert_diag.resize(M.rows());
      assert_ext(invert_diag.size() == M.rows(),
                 "DiagonalInFacePrecondition: check invert_diag.size == M.rows");
      if (Precond)
        {
          get_diagonal(M, invert_diag);
          for (int i = 0; i < invert_diag.size(); ++i)
            {
              T item = invert_diag[i];
              // T item = M.diag(i);
              assert_ext(item >= scalar_eps<T>(),
                         "DiagonalInFacePrecondition: diag(i) >= eps");
              invert_diag[i] = 1.0 / item;
            }
        }
    }

    template <typename T, typename Mat, bool Precond>
    inline int DiagonalInFacePrecondition<T, Mat, Precond>
    ::solve(__attribute__((__unused__)) const Vec &g, Vec &z, const Vec &phi) const
    {
      if (Precond)
        {
          z.resize(phi.size());
          for (int i = 0; i < z.size(); ++i)
            {
              z[i] = phi[i] * invert_diag[i];
            }
        }
      else
        {
          z = phi;
        }
      return 0;
    }

    template <typename T, typename Mat, bool Precond>
    inline DiagonalDecouplePrecondition<T, Mat, Precond>
    ::DiagonalDecouplePrecondition(const Mat &M,
                                   const std::vector<char> &face,
                                   const Eigen::SparseMatrix<T> &J)
      : DiagonalInFacePrecondition<T, Mat, Precond>(M, face)
    {
      if (Precond && J.rows() > 0)
        {
          JMinv = J * this->invert_diag.asDiagonal();
          const Eigen::SparseMatrix<T> JMinvJt = JMinv * J.transpose();
          get_diagonal(JMinvJt, invert_JMinvJt);
          for (int i = 0; i < invert_JMinvJt.size(); ++i)
            {
              assert_ext(invert_JMinvJt[i] >= scalar_eps<T>(),
                         "DiagonalDecouplePrecondition: invert_JMinvJt[i] >= eps");
              invert_JMinvJt[i] = 1.0 / invert_JMinvJt[i];
            }
        }
    }

    template <typename T, typename Mat, bool Precond>
    inline int DiagonalDecouplePrecondition<T, Mat, Precond>::solve(const Vec &g, Vec &z, const Vec &phi) const
    {
      if (Precond && JMinv.rows() > 0)
        {
          Vec lambda = invert_JMinvJt.asDiagonal() * (JMinv * g);
          assert_ext((size_t)lambda.size() == this->face.size(),
                     "DiagonalDecouplePrecondition::solve: check lambda.size == face.size");
          for (size_t i = 0; i < this->face.size(); ++i)
            {
              if (this->face[i] == 0)
                lambda[i] = 0;
            }
          z = this->invert_diag.asDiagonal() * g - JMinv.transpose() * lambda;
        }
      else
        {
          z = phi;
        }
      return 0;
    }
  } // namespace mprgp
}  // chaos
