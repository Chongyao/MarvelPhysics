#ifndef WRAPPED_IPOPT_SOLVER_H
#define WRAPPED_IPOPT_SOLVER_H

#include <IpTNLP.hpp>
#include <memory>
#include <Eigen/Sparse>
#include <chrono>

namespace bigbang {

using Ipopt::Number;
using Ipopt::IpoptData;
using Ipopt::SolverReturn;
using Ipopt::IpoptCalculatedQuantities;

template <typename T>
class Functional;

template <typename T>
class Constraint;

class ipopt_opt_framework : public Ipopt::TNLP
{
public:
  ipopt_opt_framework(const std::shared_ptr<Functional<Number>> &obj,
                      const std::shared_ptr<Constraint<Number>> &con,
                      Number *x0);
  virtual ~ipopt_opt_framework();
  virtual bool get_nlp_info(Ipopt::Index& n, Ipopt::Index& m, Ipopt::Index& nnz_jac_g,
                            Ipopt::Index& nnz_h_lag, IndexStyleEnum& index_style);
  virtual bool get_bounds_info(Ipopt::Index n, Number* x_l, Number* x_u,
                               Ipopt::Index m, Number* g_l, Number* g_u);
  virtual bool get_starting_point(Ipopt::Index n, bool init_x, Number* x,
                                  bool init_z, Number* z_L, Number* z_U,
                                  Ipopt::Index m, bool init_lambda,
                                  Number* lambda);
  virtual bool eval_f(Ipopt::Index n, const Number* x, bool new_x, Number& obj_value);
  virtual bool eval_grad_f(Ipopt::Index n, const Number* x, bool new_x, Number* grad_f);
  virtual bool eval_g(Ipopt::Index n, const Number* x, bool new_x, Ipopt::Index m, Number* g);
  virtual bool eval_jac_g(Ipopt::Index n, const Number* x, bool new_x,
                          Ipopt::Index m, Ipopt::Index nele_jac, Ipopt::Index* iRow, Ipopt::Index *jCol,
                          Number* values);
  virtual bool eval_h(Ipopt::Index n, const Number* x, bool new_x,
                      Number obj_factor, Ipopt::Index m, const Number* lambda,
                      bool new_lambda, Ipopt::Index nele_hess, Ipopt::Index* iRow,
                      Ipopt::Index* jCol, Number* values);
  virtual void finalize_solution(SolverReturn status,
                                 Ipopt::Index n, const Number* x, const Number* z_L, const Number* z_U,
                                 Ipopt::Index m, const Number* g, const Number* lambda,
                                 Number obj_value,
				 const IpoptData* ip_data,
				 IpoptCalculatedQuantities* ip_cq);
public:
  std::vector<double> bufferA_, bufferB_;
  std::chrono::time_point<std::chrono::system_clock> tic_, toc_;
protected:
  const std::shared_ptr<Functional<Number>> &obj_;
  const std::shared_ptr<Constraint<Number>> &con_;
  const size_t dim_;
  Number *x0_;
  Eigen::SparseMatrix<Number> K_, J_, lagH_;
  size_t nnz_lagH_;
  std::vector<Eigen::Triplet<Number>> g_trips_;
};

}

#endif
