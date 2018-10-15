#ifndef MINI_NNET_H
#define MINI_NNET_H

#include <Eigen/Dense>

#include "config.h"
#include "def.h"
#include "optimizer.h"
#include <math.h>
namespace nnet {

using Mat = Eigen::MatrixXd;
using Vec = Eigen::VectorXd;
using bigbang::Functional;

typedef double (*nnet_act_func)(const double x, const short order);

inline double NNET_CONSTANT(const double x, const short order) {
  switch ( order ) {
    case 0: return 1;
    case 1: return 0;
  }
}

inline double NNET_IDENTITY(const double x, const short order) {
  switch ( order ) {
    case 0: return x;
    case 1: return 1;
  }
}

inline double NNET_SIGMOID (const double x, const short order) {
  switch ( order ) {
    case 0: return 1/(1+exp(-x));
    case 1: return [](const double z)->double { return z/pow(1+z, 2); }(exp(-x));
  }
}

inline double NNET_TANH    (const double x, const short order) {
  switch ( order ) {
    case 0: return tanh(x);
    case 1: return 1.0-pow(tanh(x), 2);
  }
}

inline double NNET_RELU    (const double x, const short order) {
  switch ( order ) {
    case 0: return std::max(x, 0.0);
    case 1: return x < 0 ? 0 : 1;
  }
}

typedef double (*nnet_loss)(const double fx, const double y, const short order);

inline double NNET_LOSS_MSE(const double fx, const double y, const short order) {
  switch ( order ) {
    case 0: return pow(fx-y, 2)/2.0;
    case 1: return fx-y;
  }
}

inline double NNET_LOSS_CE (const double fx, const double y, const short order) {
  switch ( order ) {
    case 0: return -y*log(fx);
    case 1: return -y/fx;
  }
}

struct nnet_layer
{
  nnet_layer(const size_t dim) : dim_(dim) {
    x_ = y_ = dEdx_ = dEdy_ = Vec::Zero(dim_);
    dydx_ = Mat::Identity(dim_, dim_);
  }
  virtual void forward() = 0;
  virtual void backward() = 0;
  const size_t dim_;
  Vec x_, y_, dEdx_, dEdy_;
  Mat dydx_;
  nnet_act_func activate_;
};

struct input_layer : nnet_layer
{
  input_layer(const size_t dim) : nnet_layer(dim) {}
  void forward() {
    y_ = x_;
    y_(dim_-1) = 1;
  }
  void backward() {
    dydx_.setIdentity();
    dydx_(dim_-1, dim_-1) = 0;
  }
};

struct hidden_layer : nnet_layer
{
  hidden_layer(const size_t dim, nnet_act_func sigma) : nnet_layer(dim) {
    activate_ = sigma;
  }
  void forward() {
    for (size_t i = 0; i < dim_-1; ++i)
      y_[i] = activate_(x_[i], 0);
    y_[dim_-1] = 1;
  }
  void backward() {
    dydx_.setZero();
    for (size_t i = 0; i < dim_-1; ++i)
      dydx_(i, i) = activate_(x_[i], 1);
  }
};

struct output_layer_identity : nnet_layer
{
  output_layer_identity(const size_t dim) : nnet_layer(dim) {}
  void forward() {
    y_ = x_;
  }
  void backward() {
    dydx_.setIdentity();
  }
};

struct output_layer_softmax : nnet_layer
{
  output_layer_softmax(const size_t dim) : nnet_layer(dim) {}
  void forward() {
    const double expsum = x_.unaryExpr( [](double r)->double {return exp(r);} ).sum();
    y_ = x_.unaryExpr( [](double r)->double {return exp(r);})/expsum;
  }
  void backward() {
    const double expsum = x_.unaryExpr( [](double r)->double {return exp(r);} ).sum();
    for (size_t i = 0; i < dim_; ++i) {
      for (size_t j = 0; j < dim_; ++j) {
        const double expx = exp(x_[i]);
        if ( j == i )
          dydx_(i, i) = (expsum-expx)*expx/(expsum*expsum);
        else
          dydx_(i, j) = -exp(x_[j])*expx/(expsum*expsum);
      }
    }
  }
};

class neural_network
{
public:
  neural_network() : dimW_(0) {}
  virtual ~neural_network() {}
  // default configuration for regression
  virtual size_t configure_net(const size_t dimI, const size_t dimO) {
    // design your own network structure here
    loss_ = NNET_LOSS_MSE;

    // set layers of the nnet
    layers_.resize(3);
    size_t curr_layer;
    
    // input layer
    curr_layer = 0;
    layers_[curr_layer] = std::make_shared<input_layer>(dimI+1);

    // hidden layer
    curr_layer = 1;
    layers_[curr_layer] = std::make_shared<hidden_layer>(10+1, &NNET_SIGMOID);
     
    // output layer
    curr_layer = 2;
    layers_[curr_layer] = std::make_shared<output_layer_identity>(dimO);
    
    return dimW_ = allocate_memory_for_weights();
  }
  size_t get_weight_num() const {
    return dimW_;
  }
  void set_weights(const double *w, const size_t dim) {
    ASSERT(dim == dimW_);

    size_t cnt = 0;
    for (size_t i = 0; i < W_.size(); ++i) {
      std::copy(w+cnt, w+cnt+W_[i]->size(), W_[i]->data());
      cnt += W_[i]->size();
    }
  }
  int predict(const Vec &input, Vec &out) {
    ASSERT(input.size()+1 == layers_.front()->dim_);

    layers_.front()->x_.setZero();
    layers_.front()->x_.head(input.size()) = input;
    layers_.front()->forward();
    
    for (size_t i = 1; i < layers_.size(); ++i) {
      layers_[i]->x_ = (*W_[i-1])*layers_[i-1]->y_;
      layers_[i]->forward();
    }

    out = layers_.back()->y_;
  }
  void forward(const Vec &input, const Vec &target, double *val) {
    ASSERT(target.size() == layers_.back()->dim_);
    
    Vec output;
    this->predict(input, output);

    // dump loss function value
    *val = 0;
    for (size_t i = 0; i < target.size(); ++i)
      *val += loss_(output[i], target[i], 0);
  }
  void backward(const Vec &target, double *jac, const size_t dim) {
    ASSERT(target.size() == layers_.back()->dim_);
    ASSERT(dim == dimW_);

    for (size_t i = 0; i < layers_.back()->dEdy_.size(); ++i)
      layers_.back()->dEdy_[i] = loss_(layers_.back()->y_[i], target[i], 1);
    layers_.back()->backward();
    layers_.back()->dEdx_ = layers_.back()->dydx_.transpose()*layers_.back()->dEdy_;

    for (int i = dEdW_.size()-1; i >= 0; --i) {
      *dEdW_[i] = layers_[i+1]->dEdx_*layers_[i]->y_.transpose();      
      // compute previous layer dE/dy
      layers_[i]->dEdy_ = W_[i]->transpose()*layers_[i+1]->dEdx_;
      // compute previous layer dy/dx
      layers_[i]->backward();
      // compute previous layer dE/dx_j = dE/dy_i*dy_i/dx_j
      layers_[i]->dEdx_ = layers_[i]->dydx_.transpose()*layers_[i]->dEdy_;
    }

    // dump loss function gradient
    size_t cnt = 0;
    for (size_t i = 0; i < dEdW_.size(); ++i) {
      std::copy(dEdW_[i]->data(), dEdW_[i]->data()+dEdW_[i]->size(), jac+cnt);
      cnt += dEdW_[i]->size();
    }
  }
protected:
  size_t dimW_;
  std::vector<std::shared_ptr<nnet_layer>> layers_;
  nnet_loss loss_;
  size_t allocate_memory_for_weights() {
    W_.resize(layers_.size()-1);
    dEdW_.resize(layers_.size()-1);
    size_t cnt = 0;
    for (size_t i = 0; i < W_.size(); ++i) {
      W_[i]    = std::make_shared<Mat>(layers_[i+1]->dim_, layers_[i]->dim_); W_[i]->setZero();
      dEdW_[i] = std::make_shared<Mat>(layers_[i+1]->dim_, layers_[i]->dim_); dEdW_[i]->setZero();
      cnt += W_[i]->size();
    }
    return cnt;
  }
private:
  std::vector<std::shared_ptr<Mat>> W_, dEdW_;
};

class nnet_loss_func : public bigbang::Functional<double>
{
public:
  typedef Eigen::Triplet<double> TPL;
  nnet_loss_func(const std::shared_ptr<neural_network> &net, const Mat &in, const Mat &out)
      : net_(net), in_(in), out_(out), dim_(net->get_weight_num()) {
    ASSERT(in_.cols() == out_.cols());
    g_ = Vec::Zero(dim_);
  }
  size_t Nx() const {
    return dim_;
  }
  int Val(const double *w, double *val) const {
    net_->set_weights(w, dim_);

    double loss = 0; Vec tmp_grad = Vec::Zero(dim_);
    const_cast<Vec &>(g_).setZero();
    for (size_t i = 0; i < in_.cols(); ++i) {
      net_->forward(in_.col(i), out_.col(i), &loss);
      *val += loss;
      net_->backward(out_.col(i), tmp_grad.data(), tmp_grad.size());
      const_cast<Vec &>(g_) += tmp_grad;
    }
    
    return 0;
  }
  int Gra(const double *w, double *gra) const {
    Eigen::Map<Vec> G(gra, dim_);
    G += g_;
    return 0;
  }
  int Hes(const double *w, std::vector<TPL> *hes) const {
    return __LINE__;
  }
private:
  const size_t dim_;
  const std::shared_ptr<neural_network> &net_;
  const Mat &in_, &out_;
  Vec g_;
};

int nnet_train(const std::shared_ptr<neural_network> &net, const Mat &in, const Mat &out, Vec &w) {
  std::shared_ptr<Functional<double>> func
      = std::make_shared<nnet_loss_func>(net, in, out);

  double prev_err = 0;
  func->Val(w.data(), &prev_err);
  std::cout << "# prev error: " << prev_err << std::endl;
  
  bigbang::opt_args args = {500000, 1e-8, true};
  bigbang::lbfgs_solve(w.data(), w.size(), func, args);

  double post_err = 0;
  func->Val(w.data(), &post_err);
  std::cout << "# post error: " << post_err << std::endl;
  
  return 0;
}

int nnet_predict(const std::shared_ptr<neural_network> &net, const Vec &w, const Vec &in, Vec &out) {
  ASSERT(w.size() == net->get_weight_num());
  net->set_weights(w.data(), w.size());
  return net->predict(in, out);
}

}

#endif
