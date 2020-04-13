/** -*- mode: c++ -*-
 * @file mprgp_projector.h
 * @author LamKamhang (Cool_Lam@outlook.com)
 * @brief A Documented file.
 * @version 1.0
 * @date Mon Nov 18 20:16:42 CST 2019
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

#include "SimpleQP/SimpleQP.h"

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <vector>

namespace chaos
{
  namespace mprgp
  {
    // mprgp Projector methods
    template <typename T>
    class MPRGPProjectorBase
    {
      using Vec = Eigen::Matrix<T, Eigen::Dynamic, 1>;

    public:
      virtual ~MPRGPProjectorBase() = default;
      virtual const std::vector<char> &get_face() const;

      // largest step in the -d direction
      virtual T step_limit(const Vec &x, const Vec &d,
                           const T alpha_cg = scalar_max<T>()) const = 0;

      // the three major operations.
      // project the point to the feasible domain
      virtual void project(const Vec &x, Vec &y) const = 0;
      // free gradient
      virtual void get_phi(const Vec &g, Vec &phi) const = 0;
      // chopped gradient
      virtual void get_beta(const Vec &g, const Vec &phi, Vec &beta) const = 0;

      virtual void decide_face(const Vec &x) = 0;

      // calculate phi^T * phi under the Projector method.
      virtual T phiTphi(const Vec &x, const T &alpha_bar, const Vec &phi) const = 0;

      // whether the point in the feasible domain.
      virtual bool is_feasible(const Vec &x) const = 0;

    protected:
      std::vector<char> face;
    };

    template <typename T>
    class LowerBoundProjector : public MPRGPProjectorBase<T>
    {
      using Vec = Eigen::Matrix<T, Eigen::Dynamic, 1>;
      using MPRGPProjectorBase<T>::face;

    public:
      LowerBoundProjector(const Vec &L);
      virtual ~LowerBoundProjector() = default;

      // largest step in the -d direction
      virtual T step_limit(const Vec &x, const Vec &d,
                           const T alpha_cg = scalar_max<T>()) const;

      // the three major operations.
      // project the point to the feasible domain
      virtual void project(const Vec &x, Vec &y) const;
      // free gradient
      virtual void get_phi(const Vec &g, Vec &phi) const;
      // chopped gradient
      virtual void get_beta(const Vec &g, const Vec &phi, Vec &beta) const;

      virtual void decide_face(const Vec &x);

      // calculate phi^T * phi under the Projector method.
      virtual T phiTphi(const Vec &x, const T &alpha_bar, const Vec &phi) const;

      // whether the point in the feasible domain.
      virtual bool is_feasible(const Vec &x) const;

    private:
      const Vec &L;
    };

    template <typename T>
    class BoxBoundProjector : public MPRGPProjectorBase<T>
    {
      using Vec = Eigen::Matrix<T, Eigen::Dynamic, 1>;
      using MPRGPProjectorBase<T>::face;

    public:
      BoxBoundProjector(const Vec &L, const Vec &H);
      virtual ~BoxBoundProjector() = default;

      // largest step in the -d direction
      virtual T step_limit(const Vec &x, const Vec &d,
                           const T alpha_cg = scalar_max<T>()) const;

      // the three major operations.
      // project the point to the feasible domain
      virtual void project(const Vec &x, Vec &y) const;
      // free gradient
      virtual void get_phi(const Vec &g, Vec &phi) const;
      // chopped gradient
      virtual void get_beta(const Vec &g, const Vec &phi, Vec &beta) const;

      virtual void decide_face(const Vec &x);

      // calculate phi^T * phi under the Projector method.
      virtual T phiTphi(const Vec &x, const T &alpha_bar, const Vec &phi) const;

      // whether the point in the feasible domain.
      virtual bool is_feasible(const Vec &x) const;

    private:
      const Vec &L;
      const Vec &H;
    };

    template <typename T>
    class MPRGPGeneralConProjectorBase : public MPRGPProjectorBase<T>
    {
      using Vec = Eigen::Matrix<T, Eigen::Dynamic, 1>;

    protected:
      using MPRGPProjectorBase<T>::face;

    public:
      MPRGPGeneralConProjectorBase(const Eigen::SparseMatrix<T> &J, const Vec &c);
      virtual ~MPRGPGeneralConProjectorBase() = default;

      const Eigen::SparseMatrix<T> &get_J() const;

      // largest step in the -d direction
      virtual T step_limit(const Vec &x, const Vec &d,
                           const T alpha_cg = scalar_max<T>()) const;

      virtual void decide_face(const Vec &x);

      // calculate phi^T * phi under the Projector method.
      virtual T phiTphi(const Vec &x, const T &alpha_bar, const Vec &phi) const;

      // whether the point in the feasible domain.
      virtual bool is_feasible(const Vec &x) const;

    public:
      // the three major operations.
      // project the point to the feasible domain
      virtual void project(const Vec &x, Vec &y) const = 0;
      // free gradient
      virtual void get_phi(const Vec &g, Vec &phi) const = 0;
      // chopped gradient
      virtual void get_beta(const Vec &g, const Vec &phi, Vec &beta) const = 0;

    protected:
      const Eigen::SparseMatrix<T> &J;
      const Vec &c;
      size_t num_active_cons;
    };

    template <typename T>
    class GeneralConProjector : public MPRGPGeneralConProjectorBase<T>
    {
      using Vec = Eigen::Matrix<T, Eigen::Dynamic, 1>;
      using Mat = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
      using Tri = Eigen::Triplet<T>;

      using MPRGPGeneralConProjectorBase<T>::face;
      using MPRGPGeneralConProjectorBase<T>::J;
      using MPRGPGeneralConProjectorBase<T>::c;
      using MPRGPGeneralConProjectorBase<T>::num_active_cons;

    public:
      GeneralConProjector(const Eigen::SparseMatrix<T> &J, const Vec &c);
      virtual ~GeneralConProjector() = default;

      // the three major operations.
      // project the point to the feasible domain
      virtual void project(const Vec &x, Vec &y) const;
      // free gradient
      virtual void get_phi(const Vec &g, Vec &phi) const;
      // chopped gradient
      virtual void get_beta(const Vec &g, const Vec &phi, Vec &beta) const;

      virtual void decide_face(const Vec &x);

    public:
      const Eigen::SparseMatrix<T> &get_active_J() const;

    private:
      Eigen::SparseMatrix<T> active_J, active_JJt, JJt;
    };

    template <typename T>
    class DecoupleConProjector : public MPRGPGeneralConProjectorBase<T>
    {
      using Vec = Eigen::Matrix<T, Eigen::Dynamic, 1>;
      using Mat = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
      using Tri = Eigen::Triplet<T>;

      using MPRGPGeneralConProjectorBase<T>::face;
      using MPRGPGeneralConProjectorBase<T>::J;
      using MPRGPGeneralConProjectorBase<T>::c;
      using MPRGPGeneralConProjectorBase<T>::num_active_cons;

    public:
      DecoupleConProjector(const Eigen::SparseMatrix<T> &J,
                           const Vec &JJt,
                           const Vec &c);
      virtual ~DecoupleConProjector() = default;

      // the three major operations.
      // project the point to the feasible domain
      virtual void project(const Vec &x, Vec &y) const;
      // free gradient
      virtual void get_phi(const Vec &g, Vec &phi) const;
      // chopped gradient
      virtual void get_beta(const Vec &g, const Vec &phi, Vec &beta) const;

    private:
      // diagonal elements of J*J^T
      const Vec &JJt;
    };

    /////////////////////////////////////////////////////////////////////////////
    //                         template implementation                         //
    /////////////////////////////////////////////////////////////////////////////
    template <typename T>
    inline const std::vector<char> &MPRGPProjectorBase<T>::get_face() const { return face; }

    template <typename T>
    inline LowerBoundProjector<T>::LowerBoundProjector(const Vec &L)
      : L(L)
    {
      face.assign(L.size(), 0);
    }

    template <typename T>
    inline T LowerBoundProjector<T>::step_limit(const Vec &x,
                                                const Vec &d,
                                                const T __attribute__((__unused__)) alpha_cg)
      const
    {
      assert_ext(x.size() == d.size() &&
                 x.size() == L.size(),
                 "LBP::step_limit: check x.size() == d.size() == L.size()");

      T res = scalar_max<T>();
      for (size_t i = 0; i < x.size(); ++i)
        // avoid dividing a small float number
        if (d[i] > scalar_eps<T>() && x[i] > L[i])
          res = std::min<T>(res, (x[i] - L[i]) / d[i]);
      return res;
    }

    template <typename T>
    inline void LowerBoundProjector<T>::project(const Vec &x, Vec &y) const
    {
      assert_ext(x.size() == L.size(),
                 "LBP::project: check x.size() == L.size()");
      y.resize(x.size());
      assert_ext(x.size() == y.size(),
                 "LBP::project: check x.size() == y.size()");
      for (size_t i = 0; i < x.size(); ++i)
        y[i] = std::max<T>(x[i], L[i]);
    }

    template <typename T>
    inline void LowerBoundProjector<T>::get_phi(const Vec &g, Vec &phi) const
    {
      mask_face(g, phi, face);
    }

    template <typename T>
    inline void LowerBoundProjector<T>::get_beta(const Vec &g,
                                                 const Vec __attribute__((__unused__)) & phi,
                                                 Vec &beta) const
    {
      assert_ext(g.size() == face.size(),
                 "LBP::get_beta: check g.size() == face.size()");
      beta.resize(g.size());
      assert_ext(g.size() == beta.size(),
                 "LBP::get_beta: check g.size() == beta.size()");
      for (size_t i = 0; i < face.size(); ++i)
        face[i] == 0 ? beta[i] = 0.0f : beta[i] = std::min<T>(g[i], 0.0f);
    }

    template <typename T>
    inline void LowerBoundProjector<T>::decide_face(const Vec &x)
    {
      assert_ext(x.size() == L.size() &&
                 x.size() == face.size(),
                 "LBP::decide_face: check x.size() == L.size() == face.size()");
      face.assign(x.size(), 0);
      for (size_t i = 0; i < face.size(); ++i)
        if (std::abs(x[i] - L[i]) < scalar_eps<T>())
          face[i] = 2;
    }

    template <typename T>
    inline T LowerBoundProjector<T>::phiTphi(const Vec &x,
                                             const T &alpha_bar,
                                             const Vec &phi) const
    {
      assert_ext(x.size() == L.size() &&
                 x.size() == phi.size(),
                 "LBP::phiTphi: check x.size() == L.size() and x.size() == phi.size()");
      T res = 0;
      for (size_t i = 0; i < x.size(); ++i)
        {
          T phi_tilde = 0;
          if (phi[i] > scalar_eps<T>() && x[i] > L[i])
            phi_tilde = std::min<T>((x[i] - L[i]) / alpha_bar, phi[i]);
          assert_ext(phi_tilde * phi[i] >= 0,
                     "LBP::phiTphi: check phi_tilde * phi[i] >= 0");
          res += phi_tilde * phi[i];
        }
      return res;
    }

    template <typename T>
    inline bool LowerBoundProjector<T>::is_feasible(const Vec &x) const
    {
      for (size_t i = 0; i < x.size(); ++i)
        {
          if (x[i] < L[i] - scalar_eps<T>())
            return false;
        }
      return true;
    }

    template <typename T>
    inline BoxBoundProjector<T>::BoxBoundProjector(const Vec &L, const Vec &H)
      : L(L), H(H)
    {
      assert_ext(L.size() == H.size(),
                 "BBP: check L.size() == H.size()");
      assert_ext(([L, H]() {
            for (size_t i = 0; i < L.size(); ++i)
              if (L[i] > H[i])
                return false;
            return true;
          }()),
        "BBP: check L <= H");
      face.assign(L.size(), 0);
    }

    template <typename T>
    inline T BoxBoundProjector<T>::step_limit(const Vec &x,
                                              const Vec &d,
                                              const T __attribute__((__unused__)) alpha_cg)
      const
    {
      assert_ext(x.size() == d.size() &&
                 x.size() == L.size() &&
                 x.size() == H.size(),
                 "BBP::step_limit: check d == x == L == H(in size range)");
      T res = scalar_max<T>();
      for (size_t i = 0; i < x.size(); ++i)
        {
          if (d[i] > scalar_eps<T>() && x[i] > L[i])
            res = std::min<T>(res, (x[i] - L[i]) / d[i]);
          else if (d[i] < -scalar_eps<T>() && x[i] < H[i])
            res = std::min<T>(res, (x[i] - H[i]) / d[i]);
        }
      return res;
    }

    template <typename T>
    inline void BoxBoundProjector<T>::project(const Vec &x, Vec &y) const
    {
      assert_ext(x.size() == L.size() &&
                 x.size() == H.size(),
                 "BBP::project: check x.size() == L.size() == H.size()");
      y.resize(x.size());
      assert_ext(x.size() == y.size(),
                 "BBP::project: check x.size() == y.size()");
      for (size_t i = 0; i < x.size(); ++i)
        y[i] = x[i] > H[i] ? H[i] : x[i] < L[i] ? L[i] : x[i];
    }

    template <typename T>
    inline void BoxBoundProjector<T>::get_phi(const Vec &g, Vec &phi) const
    {
      mask_face(g, phi, face);
    }

    template <typename T>
    inline void BoxBoundProjector<T>::get_beta(const Vec &g,
                                               const Vec __attribute__((__unused__)) & phi,
                                               Vec &beta) const
    {
      assert_ext(g.size() == face.size(),
                 "BBP::get_beta: check g.size() == face.size()");
      beta.resize(g.size());
      assert_ext(g.size() == beta.size(),
                 "BBP::get_beta: check g.size() == beta.size()");
      for (size_t i = 0; i < face.size(); ++i)
        beta[i] = face[i] == 0 ? 0 : face[i] == 1 ? std::max<T>(g[i], 0) : std::min<T>(g[i], 0);
    }

    template <typename T>
    inline void BoxBoundProjector<T>::decide_face(const Vec &x)
    {
      assert_ext(x.size() == L.size() &&
                 x.size() == H.size() &&
                 x.size() == face.size(),
                 "BBP::decide_face: check x.size() == L.size() == H.size() == face.size()");
      face.assign(x.size(), 0);
      for (size_t i = 0; i < face.size(); ++i)
        if (abs(x[i] - L[i]) < scalar_eps<T>())
          face[i] = 2;
        else if (abs(x[i] - H[i]) < scalar_eps<T>())
          face[i] = 1;
    }

    template <typename T>
    inline T BoxBoundProjector<T>::phiTphi(const Vec &x,
                                           const T &alpha_bar,
                                           const Vec &phi) const
    {
      assert_ext(x.size() == L.size() &&
                 x.size() == H.size() &&
                 x.size() == phi.size(),
                 "BBP::phiTphi: check x == L == H == phi (in size range)");
      T res = 0;
      for (size_t i = 0; i < x.size(); ++i)
        {
          T phi_tilde = 0;
          if (phi[i] > scalar_eps<T>() && x[i] > L[i])
            phi_tilde = std::min<T>((x[i] - L[i]) / alpha_bar, phi[i]);
          else if (phi[i] < -scalar_eps<T>() && x[i] < H[i])
            phi_tilde = std::max<T>((x[i] - H[i]) / alpha_bar, phi[i]);
          assert_ext(phi_tilde * phi[i] >= 0,
                     "BBP::phiTphi: phi_tilde * phi[i] >= 0");
          res += phi_tilde * phi[i];
        }
      return res;
    }

    template <typename T>
    inline bool BoxBoundProjector<T>::is_feasible(const Vec &x) const
    {
      for (size_t i = 0; i < x.size(); ++i)
        if (x[i] < L[i] - scalar_eps<T>() ||
            x[i] > H[i] + scalar_eps<T>())
          return false;
      return true;
    }

    template <typename T>
    inline MPRGPGeneralConProjectorBase<T>::MPRGPGeneralConProjectorBase(const Eigen::SparseMatrix<T> &J,
                                                                         const Vec &c)
      : J(J), c(c), num_active_cons(0)
    {
      assert_ext(c.size() == J.rows(),
                 "MPRGPGeneralConProjectorbase: check c.size() == J.rows()");
      face.assign(c.size(), 0);
    }

    template <typename T>
    inline const Eigen::SparseMatrix<T> &
    MPRGPGeneralConProjectorBase<T>::get_J() const
    {
      return J;
    }

    template <typename T>
    inline T MPRGPGeneralConProjectorBase<T>::step_limit(const Vec &x,
                                                         const Vec &d,
                                                         const T alpha_cg)
      const
    {
      T alpha = alpha_cg == scalar_max<T>() ? alpha_cg : (alpha_cg + scalar_eps<T>());
      assert_ext(alpha > 0,
                 "MPRGPGeneralConProjectorBase::step_limit: checek alpha > 0");
      assert_ext(x.size() == d.size() &&
                 x.size() == J.cols(),
                 "MPRGPGeneralConProjectorBase::step_limit: check x.size() == d.size() == J.cols()");
      const Vec Jx = J * x;
      const Vec Jd = J * d;

      for (int i = 0; i < Jd.size(); ++i)
        if (Jd[i] > scalar_eps<T>() && Jx[i] > c[i])
          alpha = std::min<T>(alpha, (Jx[i] - c[i]) / Jd[i]);
      assert_ext(alpha >= 0,
                 "MPRGPGeneralConProjectorBase::step_limit: checek alpha >= 0");
      return alpha;
    }

    template <typename T>
    inline void MPRGPGeneralConProjectorBase<T>::decide_face(const Vec &x)
    {
      assert_ext(x.size() == J.cols(),
                 "MPRGPGeneralConProjectorBase::decide_face: x.size == J.cols");
      face.assign(c.size(), 0);
      const Vec Jx = J * x;
      num_active_cons = 0;
      for (size_t i = 0; i < face.size(); ++i)
        {
          if (abs(Jx[i] - c[i]) < scalar_eps<T>())
            {
              face[i] = 1;
              num_active_cons++;
            }
        }
    }

    template <typename T>
    inline T MPRGPGeneralConProjectorBase<T>::phiTphi(const Vec &x,
                                                      const T &alpha_bar,
                                                      const Vec &phi) const
    {
      assert_ext(alpha_bar > scalar_eps<T>(),
                 "MPRGPGeneralConProjectorBase::phiTphi: check alpha_bar > eps");
      assert_ext(x.size() == phi.size(),
                 "MPRGPGeneralConProjectorBase::phiTphi: x.size() == phi.size()");
      const Vec x_alpha_phi = x - alpha_bar * phi;
      Vec px;
      project(x_alpha_phi, px);
      const T res = ((x - px).dot(phi)) * (1.0 / alpha_bar);
      return std::max<T>(0, res);
    }

    template <typename T>
    inline bool MPRGPGeneralConProjectorBase<T>::is_feasible(const Vec &x) const
    {
      assert_ext(x.size() == J.cols(),
                 "MPRGPGeneralConProjectorBase::is_feasible: check x.size() == J.cols()");
      const Vec Jx = J * x;
      for (int i = 0; i < c.size(); ++i)
        if (Jx[i] < c[i] - scalar_eps<T>())
          return false;
      return true;
    }

    template <typename T>
    inline GeneralConProjector<T>::GeneralConProjector(const Eigen::SparseMatrix<T> &J,
                                                       const Vec &c)
      : MPRGPGeneralConProjectorBase<T>(J, c)
    {
      JJt = J * J.transpose();
    }

    template <typename T>
    inline void GeneralConProjector<T>::project(const Vec &x, Vec &y) const
    {
      if (J.rows() == 0)
        {
          y = x;
        }
      else
        {
          Vec lambda(c.size());
          const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> A = JJt;
          assert_ext(c.size() == J.rows() &&
                     x.size() == J.cols(),
                     "GCP::project: c.size == J.rows and x.size == J.cols");
          const Vec b = c - J * x;
          const Mat I = Mat::Identity(lambda.size(), lambda.size());
          Vec c0(lambda.size());
          c0.setZero();

          // TODO
          // solve a simple QP problem.
          // result is stored in lambda.
          // chaos::utils::time_point_id_t id;
          // chaos::utils::TIMER_BEGIN(id);
          solveQP(A, b, I, c0, lambda);
          // chaos::utils::TIMER_END(id, __PRETTY_FUNCTION__);

          y = x + J.transpose() * lambda;
          assert_ext(this->is_feasible(y), "GeneralConProjector::project: test y is feasible.");
        }
    }

    template <typename T>
    inline void GeneralConProjector<T>::get_phi(const Vec &g, Vec &phi) const
    {
      if (active_J.rows() == 0)
        {
          phi = g;
        }
      else
        {
          Eigen::SimplicialLLT<Eigen::SparseMatrix<T>> solver;
          solver.compute(active_JJt);
          error_msg_cond(solver.info() != Eigen::Success,
                         "LLT decomposition failed for compute phi");
          const Vec active_Jg = active_J * g;
          const Vec active_lambda = solver.solve(-active_Jg);
          error_msg_cond(solver.info() != Eigen::Success,
                         "LLT solving failed for compute phi");
          phi = g + active_J.transpose() * active_lambda;
        }
    }

    template <typename T>
    inline void GeneralConProjector<T>::get_beta(const Vec &g, const Vec &phi, Vec &beta) const
    {
      if (active_J.rows() == 0)
        {
          beta = g - phi;
        }
      else
        {
          const Vec gphi = g - phi;
          Vec lambda(active_J.rows());
          const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> A = active_JJt;
          const Vec b = active_J * gphi;
          const Mat I = Mat::Identity(lambda.size(), lambda.size());
          Vec c0(lambda.size());
          c0.setZero();

          // TODO
          // solve a simple QP problem.
          // result is stored in lambda.
          // chaos::utils::time_point_id_t id;
          // chaos::utils::TIMER_BEGIN(id);
          solveQP(A, b, I, c0, lambda);
          // chaos::utils::TIMER_END(id, __PRETTY_FUNCTION__);

          beta = gphi - active_J.transpose() * lambda;
        }
    }

    template <typename T>
    inline const Eigen::SparseMatrix<T> &
    GeneralConProjector<T>::get_active_J() const
    {
      return active_J;
    }

    template <typename T>
    inline void GeneralConProjector<T>::decide_face(const Vec &x)
    {
      assert_ext(x.size() == J.cols(),
                 "GeneralConProjector::decide_face: x.size == J.cols");
      face.assign(c.size(), 0);
      const Vec Jx = J * x;
      num_active_cons = 0;
      std::vector<Tri> triplets;
      for (size_t i = 0; i < face.size(); ++i)
        {
          if (abs(Jx[i] - c[i]) < scalar_eps<T>())
            {
              face[i] = 1;
              triplets.push_back(Tri(num_active_cons, i, 1));
              num_active_cons++;
            }
        }
      Eigen::SparseMatrix<T> P(num_active_cons, face.size());
      P.setFromTriplets(triplets.begin(), triplets.end());
      active_J = P * J;
      active_JJt = active_J * active_J.transpose();
    }

    template <typename T>
    inline DecoupleConProjector<T>::DecoupleConProjector(const Eigen::SparseMatrix<T> &J,
                                                         const Vec &JJt,
                                                         const Vec &c)
      : MPRGPGeneralConProjectorBase<T>(J, c), JJt(JJt)
    {
      assert_ext(J.rows() == JJt.size(), "DecoupleConProjector:: J.rows == JJt.size()");
    }

    template <typename T>
    inline void DecoupleConProjector<T>::project(const Vec &x, Vec &y) const
    {
      if (J.rows() == 0)
        {
          y = x;
        }
      else
        {
          assert_ext(x.size() == J.cols(), "DCP::project: x.size() == J.cols()");
          Vec lambda = c - J * x;
          for (int i = 0; i < lambda.size(); ++i)
            {
              assert_ext(JJt[i] >= scalar_eps<T>(), "DCP::project: check JJt[i] >= eps");
              lambda[i] = std::max<T>(lambda[i] / JJt[i], 0);
            }
          y = x + J.transpose() * lambda;
          assert_ext(this->is_feasible(y), "DCP::project: check whether y in feasible domain");
        }
    }

    template <typename T>
    inline void DecoupleConProjector<T>::get_phi(const Vec &g, Vec &phi) const
    {
      assert_ext(num_active_cons >= 0 &&
                 num_active_cons <= face.size(),
                 "DCP::get_phi: check num_active_cons in [0, face.size()]");
      if (num_active_cons == 0)
        {
          phi = g;
        }
      else
        {
          // lambda = -Jg / JJt
          Vec lambda = -J * g;
          assert_ext(J.cols() == g.size() &&
                     (size_t)lambda.size() == face.size(),
                     "DCP::get_phi: check J.cols() == g.size() and lambda.size() == face.size()");
          for (size_t i = 0; i < face.size(); ++i)
            {
              if (face[i] == 0)
                lambda[i] = 0;
              else
                lambda[i] /= JJt[i];
            }
          phi = g + J.transpose() * lambda;
        }
    }

    template <typename T>
    inline void DecoupleConProjector<T>::get_beta(const Vec &g, const Vec &phi, Vec &beta) const
    {
      assert_ext(num_active_cons >= 0 &&
                 num_active_cons <= face.size(),
                 "DCP::get_beta: check num_active_cons in [0, face.size()]");
      if (num_active_cons == 0)
        {
          beta = g - phi;
        }
      else
        {
          // need to optimize
          beta = g - phi;
          Vec lambda = J * beta;
          for (size_t i = 0; i < face.size(); ++i)
            {
              if (face[i] != 0)
                lambda[i] = std::max<T>(lambda[i] / JJt[i], 0);
            }
          beta -= J.transpose() * lambda;
        }
    }
  } // namespace mprgp
}  // chaos
