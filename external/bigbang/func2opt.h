#ifndef HJ_FUNC_TO_OPT_H_
#define HJ_FUNC_TO_OPT_H_

#include <jtflib/function/function.h>
#include <hjlib/math_func/math_func.h>

enum SMART_PTR{STD_SHARED, RAW};

template <typename T, SMART_PTR sp>
struct math_func_ptr{};

template <typename T>
struct math_func_ptr<T, STD_SHARED>
{typedef std::shared_ptr<T> math_func_ptr_type;};

template <typename T>
struct math_func_ptr<T, RAW>
{typedef T* math_func_ptr_type;};

template <typename VAL_TYPE, typename INT_TYPE, SMART_PTR SP>
class func2opt : public jtf::function::functionN1_t<VAL_TYPE, INT_TYPE>
{
public:
  typedef hj::math_func::math_func math_func_type;
  typedef typename math_func_ptr<math_func_type, SP>::math_func_ptr_type MFPT;

  func2opt(MFPT fptr):f_(fptr){
    assert(f_->nf() == 1);
    for(size_t k = 0; k < 3; ++k)
      cp_[k].reset(hj::math_func::patt<int32_t>(*f_, k));
    assert(cp_[0]->is_dense());
    assert(cp_[1]->is_dense());
    assert(!cp_[2]->is_dense());
  }

  virtual size_t dim(void) const{
    return f_->nx();
  }
  virtual int val(const double *x, double &v){
    return f_->eval(0,x,hj::math_func::coo2val(*cp_[0], &v));
  }

  virtual int gra(const double *x, double *g){
    return f_->eval(1, x, hj::math_func::coo2val(*cp_[1], g));
  }
  virtual int hes(const double *x, size_t &nnz, size_t &format, double *h,
                  int32_t *ptr, int32_t *idx, double alpha = 1){
    if(fabs(alpha- 1.0)) {
        std::cout << "no support for alpha." << std::endl;
        return -1;
      }
    format = 1;
    if(h == 0 && ptr == 0 && idx == 0) {// query nnz
        nnz = cp_[2]->nnz();
        return 0;
      }
    if(h == 0 && ptr != 0 && idx != 0) {// query patten
        int32_t leading[] = {0};
        hj::math_func::coo2csc(*cp_[2], ptr, idx, leading);
        return 0;
      }
    if(h != 0 && ptr != 0 && idx != 0) {// accumulate
        f_->eval(2, x, hj::math_func::coo2val(*cp_[2], h));
        return 0;
      }
  }

  virtual int gra(const double *x, size_t &nnz, double *g, int32_t *idx)
  {
    if(g == 0 && idx == 0){
        nnz = cp_[1]->nnz();
        return 0;
      }
    INT_TYPE ptr[2];
    hj::math_func::coo2csc(*cp_[1], ptr, idx);
    return f_->eval(1, x, hj::math_func::coo2val(*cp_[1],g));
  }
  virtual int hes_block(const double *x, double *h, double alpha=1)
  { assert(0); return -1; }
private:
  MFPT f_;
  std::shared_ptr<hj::math_func::coo_pat<int32_t> > cp_[3];
};

#endif
