#include <Eigen/Core>
#include <iostream>

using namespace std;
using namespace Eigen;

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
  tri << 5.14254, 5.14254, -4.85746,
      -4.95604, 5.04396, 5.04396,
      0.3, 0.3, 0.3;
  cout << tri << endl;

  Vector3d point;
  point << -0.343676, 0.441254,0.3;
  cout << check_inside(point, tri, 0);
  
}
