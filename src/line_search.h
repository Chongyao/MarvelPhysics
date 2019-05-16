#ifndef MARVEL_LINE_SEARCH
#define MARVEL_LINE_SEARCH
#include "def.h"
#include "data_str_core.h"
#include "Eigen/Core"
#include "memory"

namespace marvel{
template<typename T, size_t dim_>
double line_search(const double& val_init, const double& down,
                   const std::shared_ptr<energy_t<T, dim_>>& energy,
                   std::shared_ptr<dat_str_core<T, dim_>>& data,
                   const Eigen::VectorXd& xk, const Eigen::VectorXd& pk);
}//namespace marvel
#endif
