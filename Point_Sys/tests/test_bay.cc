#include <Eigen/Core>
#include <iostream>

using namespace std;
using namespace Eigen;

enum axis{x, y, z};
bool check_inside(const axis& dim1, const axis& dim2, const Eigen::Vector3d& point, const Eigen::MatrixXd& obstacle, const double& area){
  
  const double s = ( (point(dim1) - obstacle(dim1, 2)) * (obstacle(dim2, 0) - obstacle(dim2, 2))
                 - (point(dim2) - obstacle(dim2, 2)) * (obstacle(dim1, 0) - obstacle(dim1, 2)) ) /
             ( (obstacle(dim1, 1) - obstacle(dim1, 2)) * (obstacle(dim2, 0) - obstacle(dim2, 2))
           -(obstacle(dim2, 1) - obstacle(dim2, 2)) * (obstacle(dim1, 0) - obstacle(dim1, 2)) );
  const double t = ( (point(dim1) - obstacle(dim1, 2)) * (obstacle(dim2, 1) - obstacle(dim2, 2))
             - (point(dim2) - obstacle(dim2, 2)) * (obstacle(dim1, 1) - obstacle(dim1, 2)) ) /
      ((obstacle(dim1, 0) - obstacle(dim1, 2)) * (obstacle(dim2, 1) - obstacle(dim2, 2))
       - (obstacle(dim2, 0) - obstacle(dim2, 2)) * (obstacle(dim1, 1) - obstacle(dim1, 2)));
  
  cout<< " check through : " << s << " " << t << " " << 1 - s - t << endl;
  if(s != s){
    cout << "nan" << endl << dim1 << dim2 << endl << point << endl << endl << obstacle;
    cout << " jerer" << endl;
    assert(s == s);
  }

  // cout << obstacle << endl;
  return ( s > 0 && s < 1 && t > 0 && t < 1 && (1 - s - t) > 0 && (1 - s - t) < 1 );
        
  
}
bool check_inside(const Eigen::Vector3d& point, const Eigen::MatrixXd& obstacle, const double& area){
  // double s = (obstacle(1, 0)*obstacle(0, 2) - obstacle(0, 0)*obstacle(1, 2) + (obstacle(1, 2) - obstacle(1, 0))*point(0) + (obstacle(0, 0) - obstacle(0, 2))*point(1)) / (2*area);
  // double t = (obstacle(0, 0)*obstacle(1, 1) - obstacle(1, 0)*obstacle(0, 1) + (obstacle(1, 0) - obstacle(1, 1))*point(0) + (obstacle(0, 1) - obstacle(0, 0))*point(1)) / (2*area);
  //TODO:  divide zero
  const double s = ( (point(0) - obstacle(0, 2)) * (obstacle(1, 0) - obstacle(1, 2))
                     - (point(1) - obstacle(1, 2)) * (obstacle(0, 0) - obstacle(0, 2)) ) /
                    ( (obstacle(0, 1) - obstacle(0, 2)) * (obstacle(1, 0) - obstacle(1, 2)) -
                     (obstacle(1, 1) - obstacle(1, 2)) * (obstacle(0, 0) - obstacle(0, 2)) );
  const double t = ( (point(0) - obstacle(0, 2)) * (obstacle(1, 1) - obstacle(1, 2))
                     - (point(1) - obstacle(1, 2)) * (obstacle(0, 1) - obstacle(0, 2)) ) /
      ((obstacle(0, 0) - obstacle(0, 2)) * (obstacle(1, 1) - obstacle(1, 2))
       - (obstacle(1, 0) - obstacle(1, 2)) * (obstacle(0, 1) - obstacle(0, 2)));
  
  cout<< " check through : " << s << " " << t << " " << 1 - s - t << endl;


  return ( s > 0 && s < 1 && t > 0 && t < 1 && (1 - s - t) > 0 && (1 - s - t) < 1 );
        
  
}

int main(int argc, char** argv){
  Matrix3d tri;
  tri <<  50.1425 ,-40.8575, -40.8575,
      -40.956,   50.044,  -40.956,
      0.3,      0.3,      0.3;
  // tri <<  50.1425 ,50.1425,-40.8575
  //     -40.956 ,50.044, 50.044
  //     0.3,  0.3,  0.3;

  cout << tri << endl;
  double k = (50.044 + 40.956) / -(40.8575 + 50.1425);
  double b =50.044 - k * (-40.8575);
  
  Vector3d point;
  
  point <<-0.343676,
      0.441254,
      0.3;

  cout << point << endl;

  cout << "kx + b" << k*point(0) + b << endl;
  
  cout << check_inside(point, tri, 0);
  
}
