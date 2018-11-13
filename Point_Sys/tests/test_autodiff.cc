// C++ includes
#include <iostream>
using namespace std;

// Eigen includes
#include <eigen3/Eigen/Core>
using namespace Eigen;

// autodiff include
#define AUTODIFF_ENABLE_EIGEN_SUPPORT
#include <autodiff/autodiff/autodiff.hpp>
using namespace autodiff;

// The scalar function for which the gradient is needed
var f(const VectorXv& x, const MatrixXv &M)
{
    // return sqrt(x.cwiseProduct(x).sum()); // sqrt(sum([x(i) * x(i) for i = 1:5]))
  var res = 0;

  VectorXv test = M*x;
  cout << M << test << endl;
  res = test.sum();
  return res;
}

int main()
{
    VectorXv x(3);                         // the input vector x with 5 variables
    x << 1, 2, 3;                    // x = [1, 2, 3, 4, 5]
    MatrixXv M  (3, 3);
    M << 1, 2, 3, 4, 5, 6,7,8,9;
    var y = f(x, M);                          // the output variable y

    auto dydx = grad(y, x);            // evaluate the gradient vector dy/dx
    cout << "y = " << y << endl;           // print the evaluated output y
    cout << "dy/dx = \n" << dydx << endl;  // print the evaluated gradient vector dy/dxAP

}
