#ifndef MARVEL_HSC
#define MARVEL_HSC
#include "Sparsification.h"
#include "multigrid.h"

namespace marvel{

class HSC : public multigrid_process{
 public:
  HSC(VS<layer>& layers, const VS<transfer>& transfers);
  Eigen::VectorXd solve(const Eigen::VectorXd& b);
private:
  Eigen::VectorXd unified_multigrid(
      const size_t& curr_layer, const bool& pre_sm,
      const bool& post_sm, const size_t& num_V, const bool& diag_PD);
  
  const size_t num_layers_;
};




}
#endif
