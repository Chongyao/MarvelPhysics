#include "gen_surf.h"
#include <Eigen/Core>


deform_surf::deform_surf(const Eigen::MatrixXi &surf, const Eigen::MatirxXd &nods, const Eigen::Matrix &sam_points):surf_(surf), nods_(nods), sam_points_(sam_points), SH_(SH){}

// deform
