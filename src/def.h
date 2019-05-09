#ifndef NUMERIC_DEF_H
#define NUMERIC_DEF_H

#include <Eigen/Sparse>
#include <memory>
#include "data_str.h"
namespace marvel {

template <typename T>
class Functional
{
public:
  virtual ~Functional() {}
 // virtual size_t Nx() const = 0;
  virtual int Val(const T *x, std::shared_ptr<dat_str_core<T>>& data) const = 0;
  virtual int Gra(const T *x, std::shared_ptr<dat_str_core<T>>& data) const = 0;
  virtual int Hes(const T *x, std::shared_ptr<dat_str_core<T>>& data) const = 0;
};
#if 0
template <typename T, size_t dim>
class energy_t : public Functional<T, dim>
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

#endif

}

#endif // NUMERIC_DEF_H
