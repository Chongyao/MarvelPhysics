#include "mass_matrix.h"

#include <iostream>
#include <hjlib/math/blas_lapack.h>
#include <zjucad/matrix/lapack.h>
#include <unsupported/Eigen/KroneckerProduct>

#include "config.h"
#include "geom_util.h"

using namespace std;
using namespace zjucad::matrix;
using namespace Eigen;

namespace bigbang {

int calc_mass_matrix(const mati_t &cell,
                     const matd_t &nods,
                     const double rho,
                     const size_t dim,
                     spmat_t *M,
                     bool lumped) {
  const size_t SIMPLEX = cell.size(1);
  vector<Triplet<double>> trips;
  for (size_t i = 0; i < cell.size(2); ++i) {
    double coeff = 0;
    switch ( SIMPLEX ) {
      case 2: {
        matd_t edge = nods(colon(), cell(0, i))-nods(colon(), cell(1, i));
        double length = norm(edge);
        coeff = rho*length/6.0;
        break;
      }
      case 3: {
        matd_t edge = nods(colon(), cell(colon(1, 2), i))-nods(colon(), cell(0, i))*ones<double>(1, 2);
        double area = 0.5*norm(cross(edge(colon(), 0), edge(colon(), 1)));
        coeff = rho*area/12.0;
        break;
      }
      case 4: {
        matd_t edge = nods(colon(), cell(colon(1, 3), i))-nods(colon(), cell(0, i))*ones<double>(1, 3);
        double volume = fabs(det(edge))/6.0;
        coeff = rho*volume/20.0;
        break;
      }
      case 8: {
        double h = norm(nods(colon(), cell(0, i))-nods(colon(), cell(1, i)));
        double volume = h*h*h;
        coeff = rho*volume/72.0;
        break;
      }
      default: {
        cerr << "[info] mesh type not supported\n";
        exit(EXIT_FAILURE);
      }
    }
    for (size_t p = 0; p < cell.size(1); ++p) {
      for (size_t q = p; q < cell.size(1); ++q) {
        trips.push_back(Triplet<double>(cell(p, i), cell(q, i), coeff));
        trips.push_back(Triplet<double>(cell(q, i), cell(p, i), coeff));
      }
    }
  }
  if ( lumped ) {
#pragma omp parallel for
    for (size_t i = 0; i < trips.size(); ++i)
      trips[i] = Triplet<double>(trips[i].row(), trips[i].row(), trips[i].value());
  }
  M->resize(nods.size(2), nods.size(2));
  M->reserve(trips.size());
  M->setFromTriplets(trips.begin(), trips.end());
  if ( dim > 1 ) {
    SparseMatrix<double> I(dim, dim);
    I.setIdentity();
    SparseMatrix<double> Mtemp = kroneckerProduct(*M, I);
    *M = Mtemp;
  }
  return 0;
}

// mass matrix for 2D (3D) triangle or quad mesh
int calc_surf_mass_matrix(const mati_t &cell, const matd_t &nods,
                          const double rho, spmat_t *M) {
  ASSERT(cell.size(1) == 3 || cell.size(1) == 4);
  ASSERT(nods.size(1) == 2 || nods.size(1) == 3);

  matd_t vert_mass = zeros<double>(nods.size(2), 1);
  
  const size_t edge_num = cell.size(1);

  if ( edge_num == 3 ) { // for triangle
    for (size_t i = 0; i < cell.size(2); ++i) {
      const double area = calc_tri_area(nods(colon(), cell(colon(), i)));
      vert_mass(cell(colon(), i)) += rho*area/3;
    }    
  }

  if ( edge_num == 4 ) { // for quad
    for (size_t i = 0; i < cell.size(2); ++i) {
      const double area = calc_tri_area(nods(colon(), cell(colon(0, 2), i)))
          +calc_tri_area(nods(colon(), cell(colon(1, 3), i)));
      vert_mass(cell(colon(), i)) += rho*area/4;
    }
  }

  M->resize(nods.size(), nods.size());
  M->setIdentity();

  if ( nods.size(1) == 2 ) { // for two dimension
    for (size_t i = 0; i < vert_mass.size(); ++i)
      M->coeffRef(2*i+0, 2*i+0) = M->coeffRef(2*i+1, 2*i+1) = vert_mass[i];
  }

  if ( nods.size(1) == 3 ) { // for three dimension
    for (size_t i = 0; i < vert_mass.size(); ++i)
      M->coeffRef(3*i+0, 3*i+0) = M->coeffRef(3*i+1, 3*i+1) = M->coeffRef(3*i+2, 3*i+2) = vert_mass[i];
  }
  
  return 0;
}

}
