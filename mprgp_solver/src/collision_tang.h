/** -*- mode: c++ -*-
 * @file collision_tang.h
 * @author LamKamhang (Cool_Lam@outlook.com)
 * @brief A Documented file.
 * @version 1.0
 * @date Sun Apr 12 17:55:25 CST 2020
 *
 * Detailed description
 *
 *
 * @copyright Copyright (c) 2020
 */
#ifndef COLLISION_TANG_H
#define COLLISION_TANG_H

#include <memory>
#include <vector>

#include "collision_proxy.h"
// Tang's Collision Detect.
#include "collision_detect/src/CollisionAPI.h"
#include "utils/logger/assert_utils.h"

namespace chaos
{
  namespace collision
  {
    template<typename T>
    class collision_tang: public collision_proxy<T>
    {
      using sMat = Eigen::SparseMatrix<T>;
      using meshType = std::vector<size_t>;
      using vtxType = std::vector<T>;

      using collision_proxy<T>::meshes;
      using collision_proxy<T>::vtxs;
      using collision_proxy<T>::J;
      using collision_proxy<T>::c;
      using collision_proxy<T>::starts;
      using collision_proxy<T>::total_v;
    public:
      collision_tang() = default;
      collision_tang(
        const std::vector<meshType> &meshes,
        const std::vector<vtxType> &vtxs)
        : collision_proxy<T>(meshes, vtxs)
        {}

    public:
      virtual size_t add_model(
        const meshType &mesh,
        const vtxType &vtxs);

      virtual int detect(
        const std::vector<vtxType> &nvtxs);

    private:
      virtual void _setup(
        const std::vector<meshType> &meshes,
        const std::vector<vtxType> &vtxs);

    private:
      std::shared_ptr<Collision> coll_ptr;
    };
  } // namespace collision
} // namespace chaos


////////////////////////////////////////////////////////////////////////
//                    template implementation                         //
////////////////////////////////////////////////////////////////////////
namespace chaos
{
  namespace collision
  {
    template<typename T>
    size_t collision_tang<T>::add_model(
      const meshType &mesh,
      const vtxType &vtxs)
    {
      size_t id = this->meshes.size();
      this->meshes.push_back(mesh);
      this->vtxs.push_back(vtxs);
      this->starts.push_back(this->total_v);
      this->total_v += vtxs.size();

      for (size_t i = 0; i < id; ++i) {
        coll_ptr->Transform_Pair(i, id);
      }

      coll_ptr->Transform_Mesh(
        vtxs.size()/3, mesh.size()/3,
        mesh, vtxs, vtxs, id, true);
      return id;
    }

    template<typename T>
    int collision_tang<T>::detect(const std::vector<vtxType> &nvtxs)
    {
      assert(nvtxs.size() == vtxs.size());
      for (size_t i = 0; i < meshes.size(); ++i) {
        coll_ptr->Transform_Mesh(
          vtxs[i].size()/3, meshes[i].size()/3,
          meshes[i], nvtxs[i], vtxs[i], i, true);
      }
      coll_ptr->Collid();
      int res = coll_ptr->getCCD_res();
      if (res == 0) {
        return 0;
      } else {
        // TODO::setup constraint J and c.
        // get VF and EE from coll_ptr.
        std::vector<IDs<4>> vf;
        std::vector<IDs<6>> ee;

        this->_setup_constraint(vf, ee);
        return 1;
      }
    }

    template<typename T>
    void collision_tang<T>::_setup(
      const std::vector<meshType> &meshes,
      const std::vector<vtxType> &vtxs)
    {
      assert_ext(meshes.size() == vtxs.size(),
                 "the num of model mesh should be "
                 "equal the num of model vtx.");

      this->meshes = meshes;
      this->vtxs = vtxs;
      this->total_v = 0;
      coll_ptr = Collision::getInstance();

      // transform pairs and initial meshes.
      for (size_t i = 0; i < meshes.size(); ++i) {
        for (size_t j = 0; j < i; ++j) {
          coll_ptr->Transform_Pair(j, i);
        }
        // TODO::assume that the vtxType is just a simple 1-dim array.
        // TODO::assume that the meshType is also a simple 1-dim array.
        coll_ptr->Transform_Mesh(
          vtxs[i].size()/3, meshes[i].size()/3,
          meshes[i], vtxs[i], vtxs[i], i, true);
        this->starts.push_back(this->total_v);
        this->total_v += vtxs[i].size();
      }
    }
  } // namespace collision
} // namespace chaos


#endif /* COLLISION_TANG_H */
