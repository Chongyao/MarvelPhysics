#ifndef MARVEL_LINE_SEARCH
#define MARVEL_LINE_SEARCH
#include "def.h"
#include "data_str_core.h"
#include <Eigen/Core>
#include <memory>

namespace marvel{
template<typename T, size_t dim_>
T line_search(const T& val_init, const T& down,
              const std::shared_ptr<Functional<T, dim_>>& energy,
              std::shared_ptr<dat_str_core<T, dim_>>& data,
              const Eigen::Matrix<T, -1, 1>& xk, const Eigen::Matrix<T, -1, 1>& pk);

// extern template  double line_search(const double& val_init, const double& down,
//               const std::shared_ptr<Functional<double, 3>>& energy,
//               std::shared_ptr<dat_str_core<double, 3>>& data,
//               const Eigen::Matrix<double, -1, 1>& xk, const Eigen::Matrix<double, -1, 1>& pk);




}//namespace marvel
#include "line_search_imp.h"

#endif
