#ifndef FEM_COLLISION_H
#define FEM_COLLISION_H
#include "def.h"
#include "data_str_core.h"

namespace marvel {
template<typename T, size_t dim_>
using data_ptr = std::shared_ptr<dat_str_core<T, dim_>>;

template<typename T, size_t dim_>
class signed_dist_func
{
public:
  virtual bool inside(const T *x) const = 0;
  virtual void Val(const T *x, T * val) const = 0;  // return d*d
  virtual void Gra(const T *x, T * gra) const = 0;
  virtual void Hes(const T *x, T * hes) const = 0;
  virtual void update_center_position(const T* delt_x) = 0;
};
template<typename T, size_t dim_>
class planeSDF : public signed_dist_func<T, dim_>
{
public:
  planeSDF(const T *center, const T *n);
  bool inside(const T *x) const;
  void Val(const T *x, T *val) const;
  void Gra(const T *x, T *gra) const;
  void Hes(const T *x, T *hes) const;
  void update_center_position(const T* delt_x);
  
private:
   Eigen::Matrix<T, dim_, 1> C_, N_;
};

template<typename T, size_t dim_>
class sphereSDF : public signed_dist_func<T, dim_>
{
public:
  sphereSDF(const T *center, const T r);
  bool inside(const T *x) const;
  void Val(const T *x, T * val) const ;  // return d*d
  void Gra(const T *x, T * gra) const ;
  void Hes(const T *x, T * hes) const ;
  void update_center_position(const T* delt_x);

private:
  Eigen::Matrix<T, dim_, 1> C_;
  const T R_;
};

template<typename T, size_t dim_>
class geom_contact_energy : public Functional<T, dim_>
{
public:
  geom_contact_energy(const std::vector<std::shared_ptr<signed_dist_func<T, dim_>>> &objs,
                      const size_t dof, const T w);
  size_t Nx() const;
  int Val(const T *x, data_ptr<T, dim_> &data) const;
  int Gra(const T *x, data_ptr<T, dim_> &data) const;
  int Hes(const T *x, data_ptr<T, dim_> &data) const;
  int update_center_position(const size_t obj_id, const T* delt_x);
private:
  const size_t dof_;
  const T w_;
  const std::vector<std::shared_ptr<signed_dist_func<T, dim_>>> &objs_;
};

}
#include "SDF.imp"

#endif
