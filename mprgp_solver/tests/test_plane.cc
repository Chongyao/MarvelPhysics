#define DEBUG_PRINT_FLAG
#define INFO_PRINT_FLAG
#include "collision_plane.h"
// #include "collision_tang.h"
using namespace std;
using namespace chaos::collision;

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
  vector<size_t> cube_mesh = {
    0, 3, 2, 2, 1, 0,
    5, 4, 0, 0, 1, 5,
    2, 6, 5, 5, 1, 2,
    4, 5, 6, 6, 7, 4,
    4, 7, 3, 3, 0, 4,
    7, 6, 2, 2, 3, 7,
  };



  advance(0.2, cube_verts);

  collision_zplane<double> coll;
  coll.init({cube_mesh}, {cube_verts});

  for (int cnt = 0; cnt < 10; cnt++) {
    info_msg("test#%d", cnt);
    // cube1
    advance(-0.04, cube_verts);
    int res = coll.detect({cube_verts});
    print(cube_verts);
    cout << res << endl;
    if (res) {
      cout << "collision!" << endl;
      break;
    }
  }

  return 0;
}
