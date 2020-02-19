/** -*- mode: c++ -*-
 * @file mprgp_solver.h
 * @author LamKamhang (Cool_Lam@outlook.com)
 * @brief A Documented file.
 * @version 1.0
 * @date Mon Nov 11 16:34:34 CST 2019
 *
 * Detailed description
 *
 *
 * @copyright Copyright (c) 2019
 */

#pragma once

#include "mprgp_precondition.h"
#include "mprgp_projector.h"
#include "mprgp_utils.h"

#include "utils/logger/log_utils.h"

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <fstream>
#include <sstream>
#include <string>
#include <vector>

namespace chaos
{
  namespace mprgp
  {
    /**
     * MPRGP framework
     * has three important operations
     *  - free gradient
     *  - chopped gradient
     *  - projection
     * Assume to solve the QP problem
     *   \forall x ---> min 1/2 x^TAx - b^Tx
     *
     * T : type of the matrix and vector
     * Mat: matrix type
     */
    template <typename T, typename Mat>
    class MPRGPBase
    {
      using Vec = Eigen::Matrix<T, Eigen::Dynamic, 1>;

    public:
      MPRGPBase(const Mat &A, const Vec &b,
                MPRGPPreconditionBase<T, Mat> &preconditioner,
                MPRGPProjectorBase<T> &projector,
                const T tol = 1e-3,
                const size_t max_it = 1000,
                Vec *ev = NULL);
      virtual ~MPRGPBase() = default;

    public:
      virtual int solve(Vec &result) = 0;

    public:
      void set_tol(const T tol);
      void set_max_it(const size_t max_it);

      size_t get_it () const;
      T get_residual() const;

      bool write_vtk(const Vec &x, const std::string &filename) const;
      bool print_face() const;
      void print_solve_info() const;

      static T specRad(const Mat &A, Vec *ev = NULL, const T &eps = 1e-3);

    protected:
      void initialize(const Vec &result);
      T computeGradients(const Vec &g, const bool comp_phi = true);
      bool proportional(const Vec &result) const;

      void CGStep(const Vec &Ap, T alpha_cg, Vec &result);
      void ExpStep(const Vec &Ap, T alpha_f, Vec &result);
      void PropStep(Vec &result);

      T projectFuncValue(const Vec &x) const;
      T getFuncValue(const Vec &x) const;
      void ExpMonotonicStep(Vec &result);
      T CGMonotonicStep(Vec &result, bool &result_is_prop);

    protected:
      // the QP target function parameters.
      Mat A;
      Vec b;

      // MPRGP precondion and MPRGP projector
      MPRGPPreconditionBase<T, Mat> &preconditioner;
      MPRGPProjectorBase<T> &projector;

      // mprgp internal parameters
      T tol;
      size_t max_it;
      T gamma, alpha_bar;

      // temporary parameters
      Vec g, p, z, beta, phi, gp, y;
      size_t iteration;
      T residual;

      // solving info
      int result_code;
      int num_cg, num_exp, num_prop;
      T func_val;
    };

    // traditional MPRGP method framework
    template <typename T, typename Mat>
    class MPRGPTraditional : public MPRGPBase<T, Mat>
    {
      using Vec = Eigen::Matrix<T, Eigen::Dynamic, 1>;
      using super = MPRGPBase<T, Mat>;

    public:
      MPRGPTraditional(const Mat &A,
                       const Vec &b,
                       MPRGPPreconditionBase<T, Mat> &preconditioner,
                       MPRGPProjectorBase<T> &projector,
                       const T tol = 1e-3,
                       const size_t max_it = 1000,
                       Vec *ev = NULL);
      virtual ~MPRGPTraditional() = default;

    public:
      virtual int solve(Vec &result);
    };

    template <typename T, typename Mat>
    class MPRGPMonotonic : public MPRGPBase<T, Mat>
    {
      using Vec = Eigen::Matrix<T, Eigen::Dynamic, 1>;
      using super = MPRGPBase<T, Mat>;

    public:
      MPRGPMonotonic(const Mat &A,
                     const Vec &b,
                     MPRGPPreconditionBase<T, Mat> &preconditioner,
                     MPRGPProjectorBase<T> &projector,
                     const T tol = 1e-3,
                     const size_t max_it = 1000,
                     Vec *ev = NULL);
      virtual ~MPRGPMonotonic() = default;

    public:
      virtual int solve(Vec &result);
    };

    template <typename T = double>
    class MPRGPLowerBoundSolver
    {
      using Vec = Eigen::Matrix<T, Eigen::Dynamic, 1>;

    public:
      template <typename Mat>
      static int solve(const Mat &A,
                       const Vec &b,
                       const Vec &L,
                       Vec &x,
                       const T tol = 1e-3,
                       const size_t max_it = 1000,
                       Vec *ev = NULL);

      template <typename Mat>
      static int solve(const Mat &A,
                       const Vec &b,
                       LowerBoundProjector<T> &projector,
                       Vec &x,
                       const T tol = 1e-3,
                       const size_t max_it = 1000,
                       Vec *ev = NULL);
    };

    template <typename T = double>
    class MPRGPBoxBoundSolver
    {
      using Vec = Eigen::Matrix<T, Eigen::Dynamic, 1>;

    public:
      template <typename Mat>
      static int solve(const Mat &A,
                       const Vec &b,
                       const Vec &L,
                       const Vec &H,
                       Vec &x,
                       const T tol = 1e-3,
                       const size_t max_it = 1000,
                       Vec *ev = NULL);

      template <typename Mat>
      static int solve(const Mat &A,
                       const Vec &b,
                       BoxBoundProjector<T> &projector,
                       Vec &x,
                       const T tol = 1e-3,
                       const size_t max_it = 1000,
                       Vec *ev = NULL);
    };

    template <typename T = double>
    class MPRGPGeneralConSolver
    {
      using Vec = Eigen::Matrix<T, Eigen::Dynamic, 1>;

    public:
      template <typename Mat, bool precond = true>
      static int solve(const Mat &A,
                       const Vec &b,
                       const Eigen::SparseMatrix<T> &J,
                       const Vec &c,
                       Vec &x,
                       const T tol = 1e-3,
                       const size_t max_it = 1000,
                       Vec *ev = NULL);

      template <typename Mat, bool precond = true>
      static int solve(const Mat &A,
                       const Vec &b,
                       GeneralConProjector<T> &projector,
                       Vec &x,
                       const T tol = 1e-3,
                       const size_t max_it = 1000,
                       Vec *ev = NULL);
    };

    template <typename T = double>
    class MPRGPDecoupleConSolver
    {
      using Vec = Eigen::Matrix<T, Eigen::Dynamic, 1>;

    public:
      template <typename Mat, bool precond = true>
      static int solve(const Mat &A,
                       const Vec &b,
                       const Eigen::SparseMatrix<T> &J,
                       const Vec &c,
                       Vec &x,
                       const T tol = 1e-3,
                       const size_t max_it = 1000,
                       Vec *ev = NULL);

      template <typename Mat, bool precond = true>
      static int solve(const Mat &A,
                       const Vec &b,
                       DecoupleConProjector<T> &projector,
                       Vec &x,
                       const T tol = 1e-3,
                       const size_t max_it = 1000,
                       Vec *ev = NULL);
    };

    /////////////////////////////////////////////////////////////////////////////
    //                         template implementation                         //
    /////////////////////////////////////////////////////////////////////////////
    template <typename T, typename Mat>
    MPRGPBase<T, Mat>::MPRGPBase(const Mat &A,
                                 const Vec &b,
                                 MPRGPPreconditionBase<T, Mat> &preconditioner,
                                 MPRGPProjectorBase<T> &projector,
                                 const T tol,
                                 const size_t max_it,
                                 Vec *ev)
      : A(A)
      , b(b)
      , preconditioner(preconditioner)
      , projector(projector)
      , tol(std::max<T>(tol, 1e-30))
      , max_it(max_it)
      , gamma(1.0f)
      , alpha_bar(1.0 / specRad(A, ev, tol))
      , iteration(0)
      , residual(0.0f)
    {
      // TODO print alpha_bar
      debug_msg("alpha_bar: %lf", alpha_bar);
    }

    template <typename T, typename Mat>
    inline void MPRGPBase<T, Mat>::set_tol(const T tol)
    {
      this->tol = std::max<T>(tol, 1e-30);
    }

    template <typename T, typename Mat>
    inline void MPRGPBase<T, Mat>::set_max_it(const size_t max_it)
    {
      this->max_it = max_it;
    }

    template <typename T, typename Mat>
    inline size_t MPRGPBase<T, Mat>::get_it() const
    {
      return iteration;
    }

    template <typename T, typename Mat>
    inline T MPRGPBase<T, Mat>::get_residual() const
    {
      return residual;
    }

    template <typename T, typename Mat>
    inline bool MPRGPBase<T, Mat>::write_vtk(const Vec &x, const std::string &filename) const
    {
      Vec points(x.size() * 4);
      points.head(x.size()) = x;
      points.segment(x.size(), x.size()) = phi + x;
      points.segment(x.size() * 2, x.size()) = beta + x;
      points.segment(x.size() * 3, x.size()) = g + x;

      std::ofstream fout(filename);
      if (!fout.is_open())
        {
          error_msg_ext("write_vtk: cannot open file: %s", filename.c_str ());
          return false;
        }

      // write head
      fout << "# vtk DataFile Version 3.1 \n";
      fout << "write phi(x), beta(x) and g(x) \nASCII\nDATASET UNSTRUCTURED_GRID\n";

      // points
      fout << "POINTS " << points.size() / 3 << " FLOAT\n";
      for (size_t i = 0; i < points.size(); i += 3)
        {
          fout << points[i + 0] << " " << points[i + 1] << " " << points[i + 2] << "\n";
        }

      // lines
      const int xp = x.size() / 3;
      const int num_lines = 3 * xp;
      fout << "CELLS " << num_lines << " " << 3 * num_lines << "\n";
      for (int i = 0; i < xp; i++)
        {
          fout << 2 << " " << i << " " << i + xp << "\n";
          fout << 2 << " " << i << " " << i + 2 * xp << "\n";
          fout << 2 << " " << i << " " << i + 3 * xp << "\n";
        }

      fout << "CELL_TYPES " << num_lines << "\n";
      for (int i = 0; i < num_lines; ++i)
        {
          fout << 3;
          if (i == num_lines - 1)
            fout << "\n";
          else
            fout << " ";
        }

      const bool succ = fout.good();
      fout.close();
      return succ;
    }

    template <typename T, typename Mat>
    inline bool MPRGPBase<T, Mat>::print_face() const
    {
      const std::vector<char> &f = projector.getFace();
      std::ostringstream oss;
      oss << "face:";
      for (const auto &item : f)
        oss << " " << item;
      info_msg("%s", oss.str().c_str());
      return true;
    }

    template <typename T, typename Mat>
    inline void MPRGPBase<T, Mat>::print_solve_info() const
    {
      // TODO print
      info_msg ("MPRGP cg steps: %d", num_cg);
      info_msg ("      exp steps: %d", num_exp);
      info_msg ("      prop steps: %d", num_prop);
      info_msg ("      total_iters: %lu", iteration);
      info_msg ("      residual: %lf", residual);
      info_msg ("      constraints: %d", count_constraints (projector.get_face ()));
    }

    template <typename T, typename Mat>
    T MPRGPBase<T, Mat>::specRad(const Mat &A, Vec *ev, const T &eps)
    {
      T delta;
      Vec tmp, tmpOut;

      if (ev != NULL && ev->size() == A.rows())
        {
          tmp = (*ev);
        }
      else
        {
          tmp.resize(A.rows());
          tmp.setRandom(); ///@todo warm start ?
          tmp.normalize();
        }

      T normTmpOut;
      //power method
      for (size_t iter = 0; iter <= 1000; iter++)
        {
          tmpOut = A * tmp;

          normTmpOut = tmpOut.norm();
          if (normTmpOut < scalar_eps<T>())
            {
              if (ev)
                *ev = tmp;
              normTmpOut = scalar_eps<T>();
              break;
            }
          // normalize
          tmpOut /= normTmpOut;
          delta = (tmpOut - tmp).norm();

          if (delta <= eps)
            {
              if (ev)
                *ev = tmp;
              break;
            }
          tmp = tmpOut;
        }
      return normTmpOut;
    }

    template <typename T, typename Mat>
    inline void MPRGPBase<T, Mat>::initialize(const Vec &result)
    {
      g = A * result - b;
      projector.decide_face(result);
      projector.get_phi(g, phi);
      preconditioner.solve(g, z, phi);

      assert_ext(g.dot(z) >= 0,
                 "initialize: g.dot(z) >= 0");
      p = z;
      result_code = -1;
      num_cg = num_exp = num_prop = iteration = 0;

      // TODO PRINT FUNC VALUE
    }

    template <typename T, typename Mat>
    inline T MPRGPBase<T, Mat>::computeGradients(const Vec &g, const bool comp_phi)
    {
      if (comp_phi)
        projector.get_phi(g, phi);
      projector.get_beta(g, phi, beta);
      gp = phi + beta;
      return gp.norm();
    }

    template <typename T, typename Mat>
    inline bool MPRGPBase<T, Mat>::proportional(const Vec &result) const
    {
      //test proportional x: beta*beta <= gamma*gamma*phi*phiTilde
      const T left = beta.dot(beta);
      const T phitphi = projector.phiTphi(result, alpha_bar, phi);
      const T right = gamma * gamma * phitphi;
      return left <= right;
    }

    template <typename T, typename Mat>
    inline void MPRGPBase<T, Mat>::CGStep(const Vec &Ap, T alpha_cg, Vec &result)
    {
      num_cg++;

      result -= alpha_cg * p;
      g -= alpha_cg * Ap;

      projector.get_phi(g, phi);
      preconditioner.solve(g, z, phi);

      const T beta = (z.dot(Ap)) / (p.dot(Ap));

      p = z - beta * p;

      // TODO PRINT FUNC VALUE
    }

    template <typename T, typename Mat>
    inline void MPRGPBase<T, Mat>::ExpStep(const Vec &Ap, T alpha_f, Vec &result)
    {
      num_exp++;

      Vec &xTmp = beta;
      xTmp = result - alpha_f * p;
      g -= alpha_f * Ap;
      projector.decide_face(xTmp);
      projector.get_phi(g, phi);

      xTmp -= alpha_bar * phi;
      projector.project(xTmp, result);

      g = A * result - b;

      projector.decide_face(result);
      projector.get_phi(g, phi);
      preconditioner.solve(g, z, phi);

      p = z;

      // TODO PRINT FUNC VALUE
    }

    template <typename T, typename Mat>
    inline void MPRGPBase<T, Mat>::PropStep(Vec &result)
    {
      num_prop++;

      const Vec &D = beta;
      Vec &AD = gp;

      AD = A * D;

      const T ddad = D.dot(AD);
      // assert_ne_ext(ddad,0,"|beta| = " << beta.norm()
      //               << ", |phi| = " << phi.norm()
      //               << ", |face| = " << COUNT_CONSTRAINTS(projector.getFace()));
      const T alpha_cg = (g.dot(D)) / ddad; // assert_eq(alpha_cg, alpha_cg);
      result -= alpha_cg * D;

      Vec &xTmp = beta;
      xTmp = result;
      projector.project(xTmp, result);
      g -= alpha_cg * AD;

      projector.decide_face(result);
      projector.get_phi(g, phi);

      preconditioner.solve(g, z, phi); // assert_eq(z, z);
      p = z;

      // TODO PRINT FUNC VALUE
    }

    // 1/2 x^T*A*x - b^T*x
    template <typename T, typename Mat>
    inline T MPRGPBase<T, Mat>::getFuncValue(const Vec &x) const
    {
      return 0.5f * x.dot(A * x) - x.dot(b);
    }

    template <typename T, typename Mat>
    inline T MPRGPBase<T, Mat>::projectFuncValue(const Vec &x) const
    {
      static Vec tempx;
      projector.project(x, tempx);
      return getFuncValue(tempx);
    }

    template <typename T, typename Mat>
    inline void MPRGPBase<T, Mat>::ExpMonotonicStep(Vec &result)
    {
      num_exp++;

      projector.project(result, y);
      g = A * y - b;
      projector.decide_face(y);
      projector.get_phi(g, phi);

      y -= alpha_bar * phi;
      projector.project(y, result);

      g = A * result - b;
      projector.decide_face(result);
      projector.get_phi(g, phi);

      preconditioner.solve(g, z, phi); // assert_eq(z, z);
      p = z;

      // TODO PRINT FUNC VALUE
    }

    template <typename T, typename Mat>
    T MPRGPBase<T, Mat>::CGMonotonicStep(Vec &result, bool &result_is_prop)
    {
      Vec &Ap = gp;
      Ap = A * p;

      const T pd = p.dot(Ap); // assert_ge(pd, 0); // pd = p^t*A*p > 0
      // assert_ge_ext(g.dot(z),0,"\ng.dot(phi) = "<<g.dot(phi)<<"\n|phi|="<<phi.norm()<<"\n|beta|="<<beta.norm()<<"\n|z|="<<z.norm());
      T alpha_cg = (z.dot(g)) / pd; // assert_ge(alpha_cg, 0);

      y = result - alpha_cg * p;
      T fy = projectFuncValue(y);
      T fx = getFuncValue(result);
      // DEBUG_LOG("fx-fy: " << fx-fy);

      result_is_prop = true;
      while (result_is_prop && residual > tol && fy <= fx && iteration < max_it)
        {
          result = y;
          g -= alpha_cg * Ap;
          projector.get_phi(g, phi);       // assert_ge(g.dot(phi),0);
          preconditioner.solve(g, z, phi); // assert_eq(z, z);
          // assert_ge_ext(g.dot(z),0,"\ng.dot(phi) = "<<g.dot(phi)<<"\n|phi|="<<phi.norm()<<"\n|beta|="<<beta.norm());
          const T beta = (z.dot(Ap)) / (p.dot(Ap)); // assert_eq(beta, beta);
          p = z - beta * p;
          residual = computeGradients(g, false);
          // DEBUG_LOG(name << " residual = " << residual);
          result_is_prop = proportional(result);

          Ap = A * p;
          const T pd = p.dot(Ap);     // assert_ge(pd, 0); // pd = p^t*A*p > 0
          alpha_cg = (z.dot(g)) / pd; // assert_ge(alpha_cg, 0);
          y = result - alpha_cg * p;
          fx = fy;
          fy = projectFuncValue(y);

          num_cg++;
          iteration++;
        }

      return residual;
    }

    template <typename T, typename Mat>
    MPRGPTraditional<T, Mat>::MPRGPTraditional(const Mat &A,
                                               const Vec &b,
                                               MPRGPPreconditionBase<T, Mat> &preconditioner,
                                               MPRGPProjectorBase<T> &projector,
                                               const T tol,
                                               const size_t max_it,
                                               Vec *ev)
      : super(A, b, preconditioner, projector, tol, max_it, ev)
    {
    }

    template <typename T, typename Mat>
    int MPRGPTraditional<T, Mat>::solve(Vec &result)
    {
      this->initialize(result);
      for (this->iteration = 0; this->iteration < this->max_it; this->iteration++)
        {
          // DEBUG_LOG(MB::name << " step "<<MB::iteration);
          this->residual = this->computeGradients(this->g);

          if (this->residual <= this->tol)
            {
              this->result_code = 0;
              break;
            }

          if (this->proportional(result))
            {

              Vec &Ap = this->gp;
              Ap = this->A * this->p;

              const T pd = this->p.dot(Ap);                   // assert_ge(pd, 0); // pd = p^t*A*p > 0
              const T alpha_cg = (this->z.dot(this->g)) / pd; // assert_ge(alpha_cg, 0);

              const T alpha_f = this->projector.step_limit(result, this->p, alpha_cg);
              // assert_ge(alpha_f, 0);

              if (alpha_cg <= alpha_f)
                {
                  this->CGStep(Ap, alpha_cg, result);
                }
              else
                {
                  this->ExpStep(Ap, alpha_f, result);
                }
            }
          else
            {
              this->PropStep(result);
            }
        }

      this->result_code = this->residual <= this->tol ? 0 : -1;
      // ERROR_LOG_COND(this->name << " is not convergent with "<< this->iteration <<
      //                " iterations."<<endl,this->result_code >= 0);

      this->print_solve_info();
      return this->result_code;
    }

    template <typename T, typename Mat>
    MPRGPMonotonic<T, Mat>::MPRGPMonotonic(const Mat &A,
                                           const Vec &b,
                                           MPRGPPreconditionBase<T, Mat> &preconditioner,
                                           MPRGPProjectorBase<T> &projector,
                                           const T tol,
                                           const size_t max_it,
                                           Vec *ev)
      : super(A, b, preconditioner, projector, tol, max_it, ev)
    {
    }

    template <typename T, typename Mat>
    int MPRGPMonotonic<T, Mat>::solve(Vec &result)
    {
      this->initialize(result);

      for (this->iteration = 0; this->iteration < this->max_it;)
        {
          this->residual = this->computeGradients(this->g);
          // DEBUG_LOG(this->name << " residual = " << this->residual);
          if (this->residual <= this->tol)
            break;

          bool result_is_prop = this->proportional(result);
          if (result_is_prop)
            {
              this->residual = this->CGMonotonicStep(result, result_is_prop);
            }

          if (this->residual <= this->tol)
            break;

          if (result_is_prop || this->y.size() <= 0 || !(this->projector.is_feasible(this->y)))
            {
              this->ExpMonotonicStep(result);
              this->iteration++;
            }
          else if (!result_is_prop)
            {
              this->g = this->A * result - this->b;

              this->PropStep(result);
              this->iteration++;
            }
        }

      this->z = result;
      this->projector.project(this->z, result);

      this->result_code = this->residual <= this->tol ? 0 : -1;

      // ERROR_LOG_COND(this->name << " is not convergent with "<< this->iteration <<
      //                " iterations."<<endl,this->result_code >= 0);

      // debug_fun(this->projector.DECIDE_FACE(result));
      this->print_solve_info();
      return this->result_code;
    }

    template <typename T>
    template <typename Mat>
    inline int MPRGPLowerBoundSolver<T>::solve(const Mat &A,
                                               const Vec &b,
                                               const Vec &L,
                                               Vec &x,
                                               const T tol,
                                               const size_t max_it,
                                               Vec *ev)
    {
      LowerBoundProjector<T> projector(L);
      const Vec init_x = x;
      projector.project(init_x, x);
      return solve(A, b, projector, x, tol, max_it, ev);
    }

    template <typename T>
    template <typename Mat>
    inline int MPRGPLowerBoundSolver<T>::solve(const Mat &A,
                                               const Vec &b,
                                               LowerBoundProjector<T> &projector,
                                               Vec &x,
                                               const T tol,
                                               const size_t max_it,
                                               Vec *ev)
    {
      DiagonalInFacePrecondition<T, Mat> preconditioner(A, projector.get_face());
      MPRGPMonotonic<T, Mat> solver(A, b, preconditioner, projector, tol, max_it, ev);
      return solver.solve(x);
    }

    template <typename T>
    template <typename Mat>
    inline int MPRGPBoxBoundSolver<T>::solve(const Mat &A,
                                             const Vec &b,
                                             const Vec &L,
                                             const Vec &H,
                                             Vec &x,
                                             const T tol,
                                             const size_t max_it,
                                             Vec *ev)
    {
      BoxBoundProjector<T> projector(L, H);
      const Vec init_x = x;
      projector.project(init_x, x);
      return solve(A, b, projector, x, tol, max_it, ev);
    }

    template <typename T>
    template <typename Mat>
    inline int MPRGPBoxBoundSolver<T>::solve(const Mat &A,
                                             const Vec &b,
                                             BoxBoundProjector<T> &projector,
                                             Vec &x,
                                             const T tol,
                                             const size_t max_it,
                                             Vec *ev)
    {
      DiagonalInFacePrecondition<T, Mat> preconditioner(A, projector.get_face());
      MPRGPTraditional<T, Mat> solver(A, b, preconditioner, projector, tol, max_it, ev);
      return solver.solve(x);
    }

    template <typename T>
    template <typename Mat, bool precond>
    inline int MPRGPGeneralConSolver<T>::solve(const Mat &A,
                                               const Vec &b,
                                               const Eigen::SparseMatrix<T> &J,
                                               const Vec &c,
                                               Vec &x,
                                               const T tol,
                                               const size_t max_it,
                                               Vec *ev)
    {
      GeneralConProjector<T> projector(J, c);
      const Vec init_x = x;
      projector.project(init_x, x);
      return solve(A, b, projector, x, tol, max_it, ev);
    }

    template <typename T>
    template <typename Mat, bool precond>
    inline int MPRGPGeneralConSolver<T>::solve(const Mat &A,
                                               const Vec &b,
                                               GeneralConProjector<T> &projector,
                                               Vec &x,
                                               const T tol,
                                               const size_t max_it,
                                               Vec *ev)
    {
      DiagonalInFacePrecondition<T, Mat, precond> precondition(A, projector.get_face());
      MPRGPMonotonic<T, Mat> solver(A, b, precondition, projector, tol, max_it, ev);
      return solver.solve(x);
    }

    template <typename T>
    template <typename Mat, bool precond>
    inline int MPRGPDecoupleConSolver<T>::solve(const Mat &A,
                                                const Vec &b,
                                                const Eigen::SparseMatrix<T> &J,
                                                const Vec &c,
                                                Vec &x,
                                                const T tol,
                                                const size_t max_it,
                                                Vec *ev)
    {
      const Eigen::SparseMatrix<T> JJt_mat = J * J.transpose();
      // assert_eq_ext(JJt_mat.nonZeros(), J.rows(), "Matrix J is not decoupled.\n" << J);
      Vec JJt;
      get_diagonal(JJt_mat, JJt);
      DecoupleConProjector<T> projector(J, JJt, c);
      const Vec init_x = x;
      projector.project(init_x, x);
      return solve(A, b, projector, x, tol, max_it, ev);
    }

    template <typename T>
    template <typename Mat, bool precond>
    inline int MPRGPDecoupleConSolver<T>::solve(const Mat &A,
                                                const Vec &b,
                                                DecoupleConProjector<T> &projector,
                                                Vec &x,
                                                const T tol,
                                                const size_t max_it,
                                                Vec *ev)
    {
      DiagonalDecouplePrecondition<T, Mat, precond>
        preconditioner(A, projector.get_face(), projector.get_J());
      MPRGPMonotonic<T, Mat> solver(A, b, preconditioner, projector, tol, max_it, ev);
      return solver.solve(x);
    }
  } // namespace mprgp
}  // chaos
