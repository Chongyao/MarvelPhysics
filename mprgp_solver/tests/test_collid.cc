#define INFO_PRINT_FLAG
#include <iostream>
#include <vector>

#include "collision_detect/src/CollisionAPI.h"
#include "logger/log_utils.h"

using namespace std;

// auto res = DEFAULT_PRINTER.set_os(chaos::utils::INFO, "infoout");
void advance(double x, vector<double> & verts, int offset = 2) {
  for (size_t i = 0; i < verts.size() / 3; ++i) {
      verts[3*i+offset] += x;
    }
}

void print(const vector<double> &verts) {
  for (size_t i = 0; i < verts.size() / 3; ++i) {
    cout << verts[3*i] << "\t" << verts[3*i+1] << "\t" << verts[3*i+2] << endl;
  }
}

int main(int argc, char *argv[])
{
  // assemble a cube
  vector<double> cube_verts = {
    0.0, 0.0, 0.0,
    0.1, 0.0, 0.0,
    0.1, 0.1, 0.0,
    0.0, 0.1, 0.0,
    0.0, 0.0, 0.1,
    0.1, 0.0, 0.1,
    0.1, 0.1, 0.1,
    0.0, 0.1, 0.1,
  };

  vector<double> pre_cube_verts;

  vector<unsigned int> cube_mesh = {
    0, 3, 2, 2, 1, 0,
    5, 4, 0, 0, 1, 5,
    2, 6, 5, 5, 1, 2,
    4, 5, 6, 6, 7, 4,
    4, 7, 3, 3, 0, 4,
    7, 6, 2, 2, 3, 7,
  };

  // set a plane
  vector<double> plane_verts = {
    -1, -1, -0.15,
     1, -1, -0.15,
     1,  1, -0.15,
    -1,  1, -0.15,
    // 0.0, 0.0, -0.15,
    // 0.1, 0.0, -0.15,
    // 0.1, 0.1, -0.15,
    // 0.0, 0.1, -0.15,
  };
  vector<unsigned int> plane_mesh = {
    // 0, 1, 2, 2, 3, 0,
    3, 0, 1, 1, 2, 3,
  };

  auto coll_ptr = Collision::getInstance();
  coll_ptr->Transform_Pair(0, 1);

  for (int cnt = 0; cnt < 10; cnt++) {
    info_msg("test#%d", cnt);
    pre_cube_verts = cube_verts;
    // cube1
    advance(-0.04, cube_verts);
    coll_ptr->Transform_Mesh(cube_verts.size() / 3, cube_mesh.size() / 3, cube_mesh, cube_verts, pre_cube_verts, 0, false);

    // plane
    info_msg("add plane");
    coll_ptr->Transform_Mesh(plane_verts.size() / 3, plane_mesh.size()/3, plane_mesh, plane_verts, plane_verts, 1, false);

    info_msg("verts");
    print(cube_verts);
    info_msg("pre_verts");
    print(pre_cube_verts);
    cout << endl;

    coll_ptr->Collid();
    auto pairs = coll_ptr->getContactPairs();
    int res = coll_ptr->getCCD_res();
    info_msg("ccd result: %d", res);
    info_msg("pairs: %d", pairs.size());

    if (res) {
      auto infos = coll_ptr->getImpactInfo();
      // for (auto & info : infos) {
      for (int i = 0; i < pairs.size(); ++i) {
        auto pair = pairs[i];
        unsigned int mesh_id1, face_id1, mesh_id2, face_id2;
        pair[0].get(mesh_id1, face_id1);
        pair[1].get(mesh_id2, face_id2);
        auto info = infos[i];
        info_msg ("\nmesh id: %u face id: %u\n"
                  "mesh id: %u face id: %u\n"
                  "vf_ee %d\n"
                  "face_id: %d %d\n"
                  "verts_id: %d %d %d %d\n"
                  "dist: %lf\n"
                  "time: %lf\n"
                  "ccdres: %d",
                  mesh_id1, face_id1, mesh_id2, face_id2,
                  info.IsVF_OR_EE,
                  info.f_id[0], info.f_id[1],
                  info.vertex_id[0], info.vertex_id[1],
                  info.vertex_id[2], info.vertex_id[3],
                  info.dist,
                  info.time,
                  info.CCDres);
        cout << "------------------" << endl;
      }
      break;
    }
    // for (auto &pair : pairs) {
    //   info_msg("pair_size: %d", pair.size());
    //   unsigned int mesh_id1, face_id1, mesh_id2, face_id2;
    //   pair[0].get(mesh_id1, face_id1);
    //   pair[1].get(mesh_id2, face_id2);
    //   info_msg ("mesh id: %u face id: %u | "
    //           "mesh id: %u face id: %u\n",
    //           mesh_id1, face_id1, mesh_id2, face_id2);
    // }
  }
  return 0;
}
