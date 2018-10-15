#ifndef FEM_COLLISION_H
#define FEM_COLLISION_H

#include "def.h"

namespace bigbang {

class signed_dist_func
{
public:
  virtual bool inside(const double *x) const = 0;
  virtual void Val(const double *x, double *val) const = 0;  // return d*d
  virtual void Gra(const double *x, double *gra) const = 0;
  virtual void Hes(const double *x, double *hes) const = 0;
};

class lineSDF : public signed_dist_func
{
public:
  lineSDF(const double *center, const double *n);
  bool inside(const double *x) const;
  void Val(const double *x, double *val) const;
  void Gra(const double *x, double *gra) const;
  void Hes(const double *x, double *hes) const;
private:
  const Eigen::Vector2d C_, N_;
};

class circleSDF : public signed_dist_func
{
public:
  circleSDF(const double *center, const double r);
  bool inside(const double *x) const;
  void Val(const double *x, double *val) const;
  void Gra(const double *x, double *gra) const;
  void Hes(const double *x, double *hes) const;
private:
  const Eigen::Vector2d C_;
  const double R_;
};

class planeSDF : public signed_dist_func
{
public:
  planeSDF(const double *center, const double *n);
  bool inside(const double *x) const;
  void Val(const double *x, double *val) const;
  void Gra(const double *x, double *gra) const;
  void Hes(const double *x, double *hes) const;
private:
  const Eigen::Vector3d C_, N_;
};

class sphereSDF : public signed_dist_func
{
public:
  sphereSDF(const double *center, const double r);
  bool inside(const double *x) const;
  void Val(const double *x, double *val) const;
  void Gra(const double *x, double *gra) const;
  void Hes(const double *x, double *hes) const;
private:
  const Eigen::Vector3d C_;
  const double R_;
};

class torusSDF : public signed_dist_func
{
public:
  torusSDF(const double *center, const double *n, const double r, const double R);
  bool inside(const double *x) const;
  void Val(const double *x, double *val) const;
  void Gra(const double *x, double *gra) const;
  void Hes(const double *x, double *hes) const;
private:
  const Eigen::Vector3d C_, N_;
  const double r_, R_;
};

class cylinderSDF : public signed_dist_func
{
public:
  cylinderSDF(const double *center, const double *n, const double r);
  bool inside(const double *x) const;
  void Val(const double *x, double *val) const;
  void Gra(const double *x, double *gra) const;
  void Hes(const double *x, double *hes) const;
private:
  const Eigen::Vector3d C_, N_;
  const double R_;
};

class geom_contact_energy : public Functional<double>
{
public:
  geom_contact_energy(const std::vector<std::shared_ptr<signed_dist_func>> &objs,
                      const size_t dim, const double w, const size_t rd=3);
  size_t Nx() const;
  int Val(const double *x, double *val) const;
  int Gra(const double *x, double *gra) const;
  int Hes(const double *x, std::vector<Eigen::Triplet<double>> *hes) const;
private:
  const size_t dim_, rd_;
  const double w_;
  const std::vector<std::shared_ptr<signed_dist_func>> &objs_;
};

}

#endif
