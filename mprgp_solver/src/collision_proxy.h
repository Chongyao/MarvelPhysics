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

#include "mprgp_solver.h"

namespace chaos
{
  namespace collision
  {
    template<size_t num>
    struct IDs{
      size_t id[num];
    };

    // the ID4 of vf is:
    // ID4 => {mid_of_V, vid_of_V, mid_of_F, fid_of_F}
    enum ID4{
      MID_OF_V = 0, VID_OF_V, MID_OF_F, FID_OF_F
    };
    // the ID6 of ee is:
    // Assume Vij is the j's vertex of the i's edge.
    // ID6 => {mid_of_E1, vid_of_V11, vid_of_V12,
    //         mid_of_E2, vid_of_V21, vid_of_V22}
    enum ID6{
      MID_OF_E1 = 0, VID_OF_V11, VID_OF_V12,
      MID_OF_E2    , VID_OF_V21, VID_OF_V22
    };


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

      size_t get_last_coll_num()
      {return last_coll_num;}
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
      // each model's start point of the vtxs.
      std::vector<size_t> starts;
      sMat J;
      Eigen::Matrix<T, Eigen::Dynamic, 1> c;
      // Attention::total_v is not only total number of vtxs,
      //            total_v is dim * (num of vtxs),
      //            which means 3 times of num of vtxs.
      size_t total_v;
      size_t last_coll_num;
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
      size_t vf_num = vf.size(), ee_num = ee.size();
      J = sMat(vf_num + ee_num, total_v);
      c = Eigen::Matrix<T, Eigen::Dynamic, 1>(vf_num+ee_num);

      auto get_point = [](const std::vector<T> &vtx,
                          size_t start)->Eigen::Matrix<T, 3, 1>{
        Eigen::Matrix<T, 3, 1> V;
        V << vtx[start], vtx[start+1], vtx[start+2];
        return V;
      };

      auto set_dim3 = [](sMat &J, size_t row, size_t col,
                         T val1, T val2, T val3){
        J.insert(row, col) = val1;
        J.insert(row, col+1) = val2;
        J.insert(row, col+2) = val3;
      };

      // VF to Mat
      for (size_t i = 0; i < vf_num; ++i) {
        // calculate the face normal.
        // TODO:: If the collision is severe, which means there
        //        are many same faces in the constraint vf,
        //        we can use a map to speed up this calculation.
        size_t midF = vf[i].id[MID_OF_F];
        size_t fid = vf[i].id[FID_OF_F];
        size_t vpos = starts[vf[i].id[MID_OF_V]] + vf[i].id[VID_OF_V]*3;
        size_t fstart = starts[midF];
        const auto &mesh = meshes[midF];
        const auto &vtx  = vtxs[midF];

        size_t fv1id = mesh[fid*3];
        size_t fv2id = mesh[fid*3+1];
        size_t fv3id = mesh[fid*3+2];
        Eigen::Matrix<T, 3, 1> FV1 = get_point(vtx, fv1id*3);
        Eigen::Matrix<T, 3, 1> FV2 = get_point(vtx, fv2id*3);
        Eigen::Matrix<T, 3, 1> FV3 = get_point(vtx, fv3id*3);
        Eigen::Matrix<T, 3, 1> V = get_point(vtxs[vf[i].id[MID_OF_V]], vf[i].id[VID_OF_V]*3);

        Eigen::Matrix<T, 3, 1> e1 = FV2 - FV1, e2 = FV3 - FV1;

        Eigen::Matrix<T, 3, 1> fn = e1.cross(e2);
        // suppose that last frame does not penetrate.
        Eigen::Matrix<T, 3, 1> e = V - FV1;
        if (e.dot(fn) < 0) {
          fn = -fn;
        }
        fn.normalize();

        // get the project point of V in the F plane.
        set_dim3(J, i, vpos, fn[0], fn[1], fn[2]);
        set_dim3(J, i, fstart+fv1id*3, -fn[0]/3, -fn[1]/3, -fn[2]/3);
        set_dim3(J, i, fstart+fv2id*3, -fn[0]/3, -fn[1]/3, -fn[2]/3);
        set_dim3(J, i, fstart+fv3id*3, -fn[0]/3, -fn[1]/3, -fn[2]/3);
      }

      // EE to Mat
      // TODO
      for (size_t i = 0; i < ee_num; ++i) {

      }
    }
  } // namespace collision
} // namespace chaos

#endif /* COLLISION_PROXY_H */
