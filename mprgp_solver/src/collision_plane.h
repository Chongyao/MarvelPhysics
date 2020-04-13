/** -*- mode: c++ -*-
 * @file collision_plane.h
 * @author LamKamhang (Cool_Lam@outlook.com)
 * @brief A Documented file.
 * @version 1.0
 * @date Mon Apr 13 16:05:32 CST 2020
 *
 * Detailed description
 *
 *
 * @copyright Copyright (c) 2020
 */
#ifndef COLLISION_PLANE_H
#define COLLISION_PLANE_H

#include <memory>
#include <vector>

#include "collision_proxy.h"

namespace chaos
{
  namespace collision
  {
    template<typename T>
    class collision_zplane: public collision_proxy<T>
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
      collision_zplane() = default;
      collision_zplane(const std::vector<meshType> &meshes,
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
    size_t collision_zplane<T>::add_model(
      const meshType &mesh,
      const vtxType &vtxs)
    {
      size_t id = this->meshes.size();
      this->meshes.push_back(mesh);
      this->vtxs.push_back(vtxs);
      this->starts.push_back(this->total_v);
      this->total_v += vtxs.size();

      return id;
    }

    template<typename T>
    int collision_zplane<T>::detect(const std::vector<vtxType> &nvtxs)
    {
      using Tri = Eigen::Triplet<T>;
      std::vector<Tri> triplets;
      int cnt = 0;
      auto set_dim3 = [](std::vector<Tri> &triplets,
                         size_t row, size_t col,
                         T val1, T val2, T val3){
        triplets.emplace_back(row, col, val1);
        triplets.emplace_back(row, col+1, val2);
        triplets.emplace_back(row, col+2, val3);
      };

      for (size_t i = 0; i < nvtxs.size() - 1; ++i) {
        const vtxType& nvtxM = nvtxs[i];
        size_t num_v = nvtxM.size() / 3;
        for (size_t j = 0; j < num_v; ++j) {
          if (nvtxM[3*j+2] < 0) {
            // setup the constraint J.
            set_dim3(triplets, cnt, starts[i]+3*j, 0, 0, 1);
            set_dim3(triplets, cnt, total_v-9, 0, 0, 1.0/3);
            set_dim3(triplets, cnt, total_v-6, 0, 0, 1.0/3);
            set_dim3(triplets, cnt, total_v-3, 0, 0, 1.0/3);
            cnt++;
          }
        }
      }
      return cnt != 0;
    }

    template<typename T>
    void collision_zplane<T>::_setup(
      const std::vector<meshType> &meshes,
      const std::vector<vtxType> &vtxs)
    {
      assert_ext(meshes.size() == vtxs.size(),
                 "the num of model mesh should be "
                 "equal the num of model vtx.");

      this->meshes = meshes;
      this->vtxs = vtxs;
      this->total_v = 0;
      for (size_t i = 0; i < meshes.size(); ++i) {
        this->starts.push_back(this->total_v);
        this->total_v += vtxs[i].size();
      }

      // add a z-plane.
      // set a plane
      vtxType plane_verts = {
        1, 0, 0,
        0, 1, 0,
        1, 1, 0,
      };
      meshType plane_mesh = {
        0, 1, 2
      };
      this->meshes.emplace_back(plane_mesh);
      this->vtxs.emplace_back(plane_verts);
      this->starts.emplace_back(this->total_v);
      this->total_v += 9;

      debug_msg("total_v: %lu", total_v);
      debug_msg("meshes.size: %lu", this->meshes.size());
    }
  } // namespace collision
} // namespace chaos


#endif /* COLLISION_PLANE_H */
