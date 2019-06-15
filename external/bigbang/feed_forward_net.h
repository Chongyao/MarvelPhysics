#ifndef FEED_FORWARD_NET_H
#define FEED_FORWARD_NET_H

#include <iostream>
#include <Eigen/Eigen>

#include "src/def.h"
#include "src/optimizer.h"

namespace bigbang {

typedef double (*ffnet_act_func)(const double);

class ffnet_l2_reg : public Functional<double>
{
public:
  typedef Eigen::MatrixXd Mat;
  typedef Eigen::VectorXd Vec;
  typedef Eigen::Triplet<double> TPL;

  ffnet_l2_reg(const size_t dim, const double w) : dim_(dim), w_(w) {}
  size_t Nx() const {
    return dim_;
  }
  int Val(const double *x, double *val) const {
    Map<const Vec> X(x, dim_);
    *val += 0.5*w_*X.squaredNorm();
    return 0;
  }
  int Gra(const double *x, double *gra) const {
    Map<const Vec> X(x, dim_);
    Map<Vec> G(gra, dim_);
    G += w_*X;
    return 0;
  }
  int Hes(const double *x, std::vector<TPL> *hes) const {
    return __LINE__;
  }
private:
  const size_t dim_;
  const double w_;
};

class ffnet_mse : public Functional<double>
{
public:
  typedef Eigen::MatrixXd Mat;
  typedef Eigen::VectorXd Vec;
  typedef Eigen::Triplet<double> TPL;
  /**
   * in:    feature vectors
   * out:   function vectors
   * dimH:  dim of hidder layer
   * act:   activation function
   * actdf: differential of activation function
   */
  ffnet_mse(const Mat &in, const Mat &out, const size_t dimH,
            ffnet_act_func act, ffnet_act_func actdf)
      :dimI_(in.rows()), dimH_(dimH), dimO_(out.rows()), act_(act), actdf_(actdf) {
    assert(in.cols() == out.cols());
    
    in_ = Mat::Ones(dimI_+1, in.cols());
    in_.topLeftCorner(in.rows(), in.cols()) = in;
    out_ = out;
    
    szA_ = dimH_*(dimI_+1);
    szB_ = dimO_*(dimH_+1);
  }
  size_t Nx() const {
    return szA_+szB_;
  }
  int Val(const double *w, double *val) const {
    Map<const Mat> A(w, dimH_, dimI_+1), B(w+szA_, dimO_, dimH_+1);
    Vec layer_hid = Vec::Ones(dimH_+1), layer_out;
    for (size_t i = 0; i < in_.cols(); ++i) {
      layer_hid.head(dimH_) = A*in_.col(i);
      layer_hid.head(dimH_) = layer_hid.head(dimH_).unaryExpr(act_).eval();
      layer_out = B*layer_hid;
      *val += 0.5*(out_.col(i)-layer_out).squaredNorm();
    }
    return 0;
  }
  int Gra(const double *w, double *gra) const {
    Map<const Mat> A(w, dimH_, dimI_+1), B(w+szA_, dimO_, dimH_+1);
    Map<Mat> gA(gra, dimH_, dimI_+1), gB(gra+szA_, dimO_, dimH_+1);
    Vec layer_hid = Vec::Ones(dimH_+1), layer_hid_df = Vec::Zero(dimH_+1), layer_out;
    for (size_t i = 0; i < in_.cols(); ++i) {
      layer_hid.head(dimH_) = layer_hid_df.head(dimH_) = A*in_.col(i);      
      layer_hid.head(dimH_) = layer_hid.head(dimH_).unaryExpr(act_).eval();
      layer_hid_df.head(dimH_) = layer_hid_df.head(dimH_).unaryExpr(actdf_).eval();
      
      layer_out = B*layer_hid;
      
      for (size_t j = 0; j < dimO_; ++j) {
        gB.row(j) += -(out_(j, i)-layer_out(j))*layer_hid.transpose(); 
        for (size_t k = 0; k < dimH_; ++k) {
          gA.row(k) += -(out_(j, i)-layer_out(j))*B(j, k)*layer_hid_df(k)*in_.col(i).transpose();
        }
      }
    }
    return 0;
  }
  int Hes(const double *w, std::vector<TPL> *hes) const {
    return __LINE__;
  }
private:
  const size_t dimI_, dimO_, dimH_;
  size_t szA_, szB_;
  ffnet_act_func act_, actdf_;
  Mat in_, out_;
};

inline double SIGMOID(const double x) {
  return 1/(1+exp(-x));
}

inline double SIGMOID_DIFF(const double x) {
  const double exp_x = exp(-x);
  return exp_x/pow(1+exp_x, 2);
}

class feed_forward_net
{
public:
  typedef Eigen::MatrixXd Mat;
  typedef Eigen::VectorXd Vec;
  feed_forward_net(const size_t dimI, const size_t dimO, const size_t dimH)
      : dimI_(dimI), dimO_(dimO), dimH_(dimH) {
    dim_ = dimH_*(dimI_+1)+dimO_*(dimH_+1);
  }
  void load(const Mat &in, const Mat &out) {
    assert(in_.rows() == dimI_ && out_.rows() == dimO_);
    assert(in_.cols() == out_.cols());
    in_ = in;
    out_ = out;
  }
  int train() {
    std::vector<std::shared_ptr<Functional<double>>> buffer(2);
    buffer[0] = std::make_shared<ffnet_l2_reg>(dim_, 0.0);
    buffer[1] = std::make_shared<ffnet_mse>(in_, out_, dimH_, &SIGMOID, &SIGMOID_DIFF);
    std::shared_ptr<Functional<double>> obj;
    obj = std::make_shared<energy_t<double>>(buffer);

    W_ = Vec::Random(obj->Nx());
    opt_args option = {500000, 1e-8, true};
    lbfgs_solve(W_.data(), obj->Nx(), obj, option);
    
    return 0;
  }
  void predict(const Vec &x, Vec &y) {
    const size_t szA = dimH_*(dimI_+1);
    Map<const Mat> A(W_.data(), dimH_, dimI_+1), B(W_.data()+szA, dimO_, dimH_+1);
    Vec x1 = Vec::Ones(dimI_+1), layer_hid = Vec::Ones(dimH_+1);
    x1.head(x.size()) = x;
    layer_hid.head(dimH_) = A*x1;
    layer_hid.head(dimH_) = layer_hid.head(dimH_).unaryExpr(&SIGMOID).eval();
    y = B*layer_hid;
  }
private:
  const size_t dimI_, dimO_, dimH_;
  size_t dim_;
  Mat in_, out_;
  Vec W_;
};

}
#endif
