#ifndef MARVEL_HSC
#define MARVEL_HSC
#include "Sparsification.h"
#include "multigrid.h"
#include <boost/property_tree/ptree.hpp>

namespace marvel{

class HSC : public multigrid_process{
 public:
  HSC(VS<layer>& layers, const VS<transfer>& transfers,const bool pre_sm, const bool post_sm, const bool diag_PD, const size_t num_mu, const size_t num_V);
  Eigen::VectorXd solve(const Eigen::VectorXd& b);
  
 private:
  void unified_multigrid(
    const size_t& curr_layer_id, const bool& pre_sm,
    const bool& post_sm, const size_t& num_mu, const bool& diag_PD);
  
  const size_t num_layers_;

  const bool pre_sm_, post_sm_, diag_PD_;
  const size_t num_mu_, num_V_;
  
};


HSC set_hierarchy(const SPM& L, const boost::property_tree::ptree& pt);




}
#endif
