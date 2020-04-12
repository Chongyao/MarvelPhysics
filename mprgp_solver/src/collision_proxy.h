/** -*- mode: c++ -*-
 * @file collision_proxy.h
 * @author LamKamhang (Cool_Lam@outlook.com)
 * @brief A Documented file.
 * @version 1.0
 * @date Sun Apr 12 13:52:57 CST 2020
 *
 * Detailed description
 *
 *
 * @copyright Copyright (c) 2020
 */
#ifndef COLLISION_PROXY_H
#define COLLISION_PROXY_H

#include <vector>
#include <memory>

#include "collision_detect/src/CollisionAPI.h"
#include "mprgp_solver.h"

namespace chaos
{
  namespace collision
  {
    template<typename T>
    class collision_proxy
    {
      using sMat = Eigen::SparseMatrix<T>;
      using meshType = std::vector<size_t>;
      using vtxType = std::vector<T>;
    public:
      collision_proxy() = default;
      collision_proxy(
        const std::vector<meshType> &meshes,
        const std::vector<vtxType> &vtxs)
        {_setup(meshes, vtxs);}

      void init(
        const std::vector<meshType> &meshes,
        const std::vector<vtxType> &vtxs)
        {_setup(meshes, vtxs);}

      void update(const std::vector<vtxType> &nvtxs)
        {this->vtxs = nvtxs;}

      // before calling this function,
      // J and c should be setuped.
      template<typename Mat, typename Vec>
      std::vector<vtxType>
      response(const std::vector<Mat> A,
               const std::vector<Vec> b,
               bool do_update = true);

      // before calling this function,
      // J and c should be setuped.
      template<typename Mat, typename Vec>
      std::vector<vtxType>
      response(const Mat &eA,
               const Vec &eb,
               bool do_update = true);
    protected:
      template<size_t num>
      struct IDs{
        size_t id[num];
      };

      // this function used to setup J and c.
      void _setup_constraint(
        // the ID4 of vf is:
        // ID4 => {mid_of_V, vid_of_V, mid_of_F, fid_of_F}
        const std::vector<IDs<4>> &vf,
        // the ID6 of ee is:
        // Assume Vij is the j's vertex of the i's edge.
        // ID6 => {mid_of_E1, vid_of_V11, vid_of_V12,
        //         mid_of_E2, vid_of_V21, vid_of_V22}
        const std::vector<IDs<6>> &ee);

    public:
      virtual size_t add_model(
        const meshType &mesh,
        const vtxType &vtxs) = 0;

      virtual int detect(
        const std::vector<vtxType> &nvtxs) = 0;

    private:
      virtual void _setup(
        const std::vector<meshType> &meshes,
        const std::vector<vtxType> &vtxs) = 0;

    protected:
      std::vector<meshType> meshes;
      std::vector<vtxType> vtxs;
      sMat J;
      Eigen::Matrix<T, Eigen::Dynamic, 1> c;
      size_t total_v;
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
    template<typename Mat, typename Vec>
    // std::vector<vtxType>
    std::vector<std::vector<T> >
    collision_proxy<T>::response(
      const std::vector<Mat> A,
      const std::vector<Vec> b,
      bool do_update)
    {
      // Assert that all sparse matrices are column major.
      // build the environment stiffness Matrix eA.
      sMat eA(total_v, total_v);
      size_t head_dim = 0, tail_dim = total_v;
      // TODO::optimize it.
      // assemble As to environmentA
      // A1, A2, ... , An => A = diag{A1, A2, ... , An};
      for (size_t i = 0; i < A.size(); ++i) {
        size_t cur_dim = A[i].cols();
        tail_dim -= cur_dim;
        Eigen::SparseMatrix<double, Eigen::RowMajor>
          col_block(total_v, cur_dim);
        col_block.middleRows(head_dim, cur_dim) = A[i];
        eA.middleCols(head_dim, cur_dim) = col_block;
        head_dim += cur_dim;
      }

      // build the environment b.
      Vec eb(total_v);
      size_t pos = 0;
      for (size_t i = 0; i < b.size(); ++i) {
        size_t num = b[i].size();
        eb.segment<num>(pos) = b[i];
        pos += num;
      }

      return response(eA, eb, do_update);
    }

    template<typename T>
    template<typename Mat, typename Vec>
    // std::vector<vtxType>
    std::vector<std::vector<T> >
    collision_proxy<T>::response(
      const Mat &eA,
      const Vec &eb,
      bool do_update)
    {
      Vec x;
      // 0 means ok.
      int res = chaos::mprgp::MPRGPDecoupleConSolver<T>::solve(
        eA, eb, J, c, x);

      if (res == 0) {
        // success and return the new vtxs.
        // TODO::optimize it.
        size_t pos = 0;
        if (do_update) {
          for (size_t i = 0; i < vtxs.size(); ++i) {
            for (size_t j = 0; j < vtxs[i].size(); ++j) {
              vtxs[i][j] = x[pos++];
            }
          }
          return vtxs;
        } else {
          std::vector<vtxType> nvtxs;
          for (size_t i = 0; i < vtxs.size(); ++i) {
            vtxType v;
            for (size_t j = 0; j < vtxs[i].size(); ++j) {
              v.push_back(x[pos++]);
            }
            nvtxs.push_back(v);
          }
          return nvtxs;
        }
      } else {
        // failed and return the last vtxs.
        return this->vtxs;
      }
    }

    template<typename T>
    void collision_proxy<T>::_setup_constraint(
      // the ID4 of vf is:
      // ID4 => {mid_of_V, vid_of_V, mid_of_F, fid_of_F}
      const std::vector<IDs<4>> &vf,
      // the ID6 of ee is:
      // Assume Vij is the j's vertex of the i's edge.
      // ID6 => {mid_of_E1, vid_of_V11, vid_of_V12,
      //         mid_of_E2, vid_of_V21, vid_of_V22}
      const std::vector<IDs<6>> &ee)
    {
      // clear the data of J and c.
      // VF to Mat
      // EE to Mat
    }
  } // namespace collision
} // namespace chaos

#include "collision_tang.h"

#endif /* COLLISION_PROXY_H */
