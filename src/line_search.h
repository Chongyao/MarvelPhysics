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
              const T* const xk, const T* const  pk);

}//namespace marvel

#endif
