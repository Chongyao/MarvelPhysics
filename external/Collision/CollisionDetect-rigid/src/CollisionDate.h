#ifndef SELFCOLLISIONDATE
#define SELFCOLLISIONDATE
#include "cmesh.h"

struct CollisionDate{
	CollisionDate(mesh* m, bool flag) {
		ms = m;
		enable_selfcollision = flag;
	}
	mesh* ms;
	bool enable_selfcollision;
};

#endif