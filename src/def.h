#ifndef NUMERIC_DEF_H
#define NUMERIC_DEF_H

#include <Eigen/Sparse>
#include <memory>

namespace bigbang {

template <typename T>
class Functional
{
public:
  virtual ~Functional() {}
  virtual size_t Nx() const = 0;
  virtual int Val(const T *x, T *val) const = 0;
  virtual int Gra(const T *x, T *gra) const = 0;
  virtual int Hes(const T *x, std::vector<Eigen::Triplet<T>> *hes) const = 0;
  virtual void ResetWeight(const double w) {}
  virtual int operator ()(const T *x, T *val, T *gra, const T step, bool graON) { // for LBFGS
    this->Val(x, val);
    if ( graON )
      this->Gra(x, gra);
    return 0;
  }
};

template <typename T>
class Constraint
{
public:
  virtual ~Constraint() {}
  virtual size_t Nx() const = 0;
  virtual size_t Nf() const = 0;
  virtual int Val(const T *x, T *val) const = 0;
  virtual int Jac(const T *x, const size_t off, std::vector<Eigen::Triplet<T>> *jac) const = 0;
  virtual int Hes(const T *x, const size_t off, std::vector<std::vector<Eigen::Triplet<T>>> *hes) const {
    return __LINE__;
  }
};

template <typename T>
class least_square_wrapper : public Functional<T>
{
public:
  least_square_wrapper(const std::shared_ptr<Constraint<T>> &cons, const T w)
      : cons_(cons), w_(w) {}
  size_t Nx() const {
    return cons_->Nx();
  }
  int Val(const T *x, T *val) const {
    Eigen::Matrix<T, -1, 1> Cv;
    Cv.setZero(cons_->Nf());
    cons_->Val(x, Cv.data());

    *val += 0.5*w_*Cv.squaredNorm();
    return 0;
  }
  int Gra(const T *x, T *gra) const {
    Eigen::Matrix<T, -1, 1> Cv;
    Cv.setZero(cons_->Nf());
    cons_->Val(x, Cv.data());

    Eigen::SparseMatrix<T> Jc(cons_->Nf(), cons_->Nx()); {
      std::vector<Eigen::Triplet<T>> trips;
      cons_->Jac(x, 0, &trips);
      Jc.setFromTriplets(trips.begin(), trips.end());
    }

    Eigen::Map<Eigen::Matrix<T, -1, 1>>(gra, Nx()) += w_*Jc.transpose()*Cv;
    return 0;
  }
  int Hes(const T *x, std::vector<Eigen::Triplet<T>> *hes) const {
    Eigen::SparseMatrix<T> Jc(cons_->Nf(), cons_->Nx()); {
      std::vector<Eigen::Triplet<T>> trips;
      cons_->Jac(x, 0, &trips);
      Jc.setFromTriplets(trips.begin(), trips.end());
    }

    Eigen::SparseMatrix<T> JcTJc = Jc.transpose()*Jc;
    for (size_t j = 0; j < JcTJc.outerSize(); ++j) {
      for (typename Eigen::SparseMatrix<T>::InnerIterator it(JcTJc, j); it; ++it) {
        hes->push_back(Eigen::Triplet<T>(it.row(), it.col(), w_*it.value()));
      }
    }

    return 0;
  }
protected:
  const std::shared_ptr<Constraint<T>> &cons_;
  const T w_;
};

template <typename T>
class linearization_wrapper : public Functional<T>
{
public:
  linearization_wrapper(const std::shared_ptr<Functional<T>> &f, const double *x)
      : dim_(f->Nx()) {
    x0_ = Eigen::Map<const Eigen::Matrix<T, -1, 1>>(x, dim_);

    val_x0_ = 0;
    f->Val(x, &val_x0_);

    gra_x0_.setZero(dim_);
    f->Gra(x, gra_x0_.data());

    hes_x0_.resize(dim_, dim_);
    std::vector<Eigen::Triplet<T>> trips;
    f->Hes(x, &trips);
    hes_x0_.setFromTriplets(trips.begin(), trips.end());    
  }
  size_t Nx() const {
    return dim_;
  }
  int Val(const double *x, double *val) const {
    const Eigen::Matrix<T, -1, 1> dx = Eigen::Map<const Eigen::Matrix<T, -1, 1>>(x, dim_)-x0_;
    *val += val_x0_+gra_x0_.dot(dx)+0.5*dx.dot(hes_x0_*dx);
    return 0;
  }
  int Gra(const double *x, double *gra) const {
    const Eigen::Matrix<T, -1, 1> dx = Eigen::Map<const Eigen::Matrix<T, -1, 1>>(x, dim_)-x0_;
    Eigen::Map<Eigen::Matrix<T, -1, 1>>(gra, dim_) += hes_x0_*dx+gra_x0_;
    return 0;
  }
  int Hes(const double *x, std::vector<Eigen::Triplet<T>> *hes) const {
    for (size_t j = 0; j < hes_x0_.outerSize(); ++j) {
      for (typename Eigen::SparseMatrix<T>::InnerIterator it(hes_x0_, j); it; ++it) {
        hes->push_back(Eigen::Triplet<T>(it.row(), it.col(), it.value()));
      }
    }

    // //-->PD
    // SparseMatrix<double> K_h, M_h;
    // VectorXd ev_h; {
    //   K_h.resize(this->Nx(), this->Nx());
    //   K_h.setFromTriplets(hes->begin(), hes->end());
    //   calc_mass_matrix(hexs_h, nods_h, 1, rd, &M_h, lumped);
    //   solve_gen_eig_prob(K_h, M_h, dim_*3, 0, "arpaca", &ev_h, nullptr);
    // }


    
    return 0;
  }
private:
  const size_t dim_;
  Eigen::Matrix<T, -1, 1> x0_;

  T val_x0_;
  Eigen::Matrix<T, -1, 1> gra_x0_;
  Eigen::SparseMatrix<T> hes_x0_;
};

template <typename T>
class energy_t : public Functional<T>
{
public:
  class null_input_exception : public std::exception {
  public :
    const char* what() const throw() {
      return "null input exception";
    }
  };
  class compatibility_exception : public std::exception {
  public :
    const char* what() const throw() {
      return "compatibility exception";
    }
  };
  energy_t(const std::vector<std::shared_ptr<Functional<T>>> &buffer)
    : buffer_(buffer), dim_(-1) {
    for (auto &e : buffer_) {
      if ( e.get() ) {
        dim_ = e->Nx();
        break;
      }
    }
    if ( dim_ == -1 ) {
      throw null_input_exception();
    }
    for (auto &e : buffer_) {
      if ( e.get() && e->Nx() != dim_ ) {
        throw compatibility_exception();
      }
    }
  }
  size_t Nx() const {
    return dim_;
  }
  int Val(const T *x, T *val) const {
    for (auto &e : buffer_) {
      if ( e.get() ) {
        e->Val(x, val);
      }
    }
    return 0;
  }
  int Gra(const T *x, T *gra) const {
    for (auto &e : buffer_) {
      if ( e.get() ) {
        e->Gra(x, gra);
      }
    }
    return 0;
  }
  int Hes(const T *x, std::vector<Eigen::Triplet<T>> *hes) const {
    for (auto &e : buffer_) {
      if ( e.get() ) {
        e->Hes(x, hes);
      }
    }
    return 0;
  }
protected:
  const std::vector<std::shared_ptr<Functional<T>>> &buffer_;
  size_t dim_;
};

template <typename T>
class constraint_t : public Constraint<T>
{
public:
  class null_input_exception : std::exception {
  public :
    const char* what() const throw() {
      return "null input exception";
    }
  };
  class compatibility_exception : std::exception {
  public :
    const char* what() const throw() {
      return "compatibility exception";
    }
  };
  constraint_t(const std::vector<std::shared_ptr<Constraint<T>>> &buffer)
    : buffer_(buffer), xdim_(-1) {
    for (auto &e : buffer_) {
      if ( e.get() ) {
        xdim_ = e->Nx();
        break;
      }
    }
    if ( xdim_ == -1 )
      throw null_input_exception();
    bool compatible = true;
    for (auto &c : buffer_) {
      if ( c.get() ) {
        if ( c->Nx() != xdim_ )
          compatible = false;
      }
    }
    if ( !compatible )
      throw compatibility_exception();
  }
  size_t Nx() const {
    return xdim_;
  }
  size_t Nf() const {
    size_t fdim = 0;
    for (auto &c : buffer_) {
      if ( c.get() )
        fdim += c->Nf();
    }
    return fdim;
  }
  int Val(const T *x, T *val) const {
    Eigen::Map<Eigen::Matrix<T, -1, 1>> v(val, Nf());
    size_t offset = 0;
    for (auto &c : buffer_) {
      if ( c.get() ) {
        const size_t nf = c->Nf();
        Eigen::Matrix<T, -1, 1> value(nf);
        value.setZero();
        c->Val(x, value.data());
        v.segment(offset, nf) += value;
        offset += nf;
      }
    }
    return 0;
  }
  int Jac(const T *x, const size_t off, std::vector<Eigen::Triplet<T>> *jac) const {
    size_t offset = off;
    for (auto &c : buffer_) {
      if ( c.get() ) {
        c->Jac(x, offset, jac);
        offset += c->Nf();
      }
    }
    return 0;
  }
  int Hes(const T *x, const size_t off, std::vector<std::vector<Eigen::Triplet<T>>> *hes) const {
    const size_t fdim = Nf();
    if ( hes->size() != fdim )
      hes->resize(fdim);
    size_t offset = 0;
    for (auto &c : buffer_) {
      if ( c.get() ) {
        c->Hes(x, offset, hes);
        offset += c->Nf();
      }
    }
    return 0;
  }
protected:
  const std::vector<std::shared_ptr<Constraint<T>>> &buffer_;
  size_t xdim_;
};



}

#endif // NUMERIC_DEF_H
