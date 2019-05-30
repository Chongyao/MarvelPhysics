#include "cmesh.h"
#include <set>
#include <iostream>
#include <stdio.h>

using namespace std;
#include "mat3f.h"
#include "box.h"
#include "tmbvh.hpp"

//rky
#include "collid.h"
//rky

bool findd;

#include <omp.h>
//rky
double tmp_timing_start;
# define	TIMING_BEGIN \
	tmp_timing_start = omp_get_wtime();\

# define	TIMING_END(message) \
	{double tmp_timing_finish = omp_get_wtime();\
	double  tmp_timing_duration = tmp_timing_finish - tmp_timing_start;\
	printf("%s: %2.5f seconds\n", (message), tmp_timing_duration);}

//#define POVRAY_EXPORT
#define OBJ_DIR "e:\\temp\\output-objs\\"

//#define VEC_CLOTH

#pragma warning(disable: 4996)


inline double fmax(double a, double b, double c)
{
  double t = a;
  if (b > t) t = b;
  if (c > t) t = c;
  return t;
}

inline double fmin(double a, double b, double c)
{
  double t = a;
  if (b < t) t = b;
  if (c < t) t = c;
  return t;
}

inline int project3(const vec3f &ax, 
	const vec3f &p1, const vec3f &p2, const vec3f &p3)
{
  double P1 = ax.dot(p1);
  double P2 = ax.dot(p2);
  double P3 = ax.dot(p3);
  
  double mx1 = fmax(P1, P2, P3);
  double mn1 = fmin(P1, P2, P3);

  if (mn1 > 0) return 0;
  if (0 > mx1) return 0;
  return 1;
}

inline int project6(vec3f &ax, 
	 vec3f &p1, vec3f &p2, vec3f &p3, 
	 vec3f &q1, vec3f &q2, vec3f &q3)
{
  double P1 = ax.dot(p1);
  double P2 = ax.dot(p2);
  double P3 = ax.dot(p3);
  double Q1 = ax.dot(q1);
  double Q2 = ax.dot(q2);
  double Q3 = ax.dot(q3);
  
  double mx1 = fmax(P1, P2, P3);
  double mn1 = fmin(P1, P2, P3);
  double mx2 = fmax(Q1, Q2, Q3);
  double mn2 = fmin(Q1, Q2, Q3);

  if (mn1 > mx2) return 0;
  if (mn2 > mx1) return 0;
  return 1;
}


// very robust triangle intersection test
// uses no divisions
// works on coplanar triangles

bool 
tri_contact (vec3f &P1, vec3f &P2, vec3f &P3, vec3f &Q1, vec3f &Q2, vec3f &Q3) 
{
  vec3f p1;
  vec3f p2 = P2-P1;
  vec3f p3 = P3-P1;
  vec3f q1 = Q1-P1;
  vec3f q2 = Q2-P1;
  vec3f q3 = Q3-P1;
  
  vec3f e1 = p2-p1;
  vec3f e2 = p3-p2;
  vec3f e3 = p1-p3;

  vec3f f1 = q2-q1;
  vec3f f2 = q3-q2;
  vec3f f3 = q1-q3;

  vec3f n1 = e1.cross(e2);
  vec3f m1 = f1.cross(f2);

  vec3f g1 = e1.cross(n1);
  vec3f g2 = e2.cross(n1);
  vec3f g3 = e3.cross(n1);

  vec3f  h1 = f1.cross(m1);
  vec3f h2 = f2.cross(m1);
  vec3f h3 = f3.cross(m1);

  vec3f ef11 = e1.cross(f1);
  vec3f ef12 = e1.cross(f2);
  vec3f ef13 = e1.cross(f3);
  vec3f ef21 = e2.cross(f1);
  vec3f ef22 = e2.cross(f2);
  vec3f ef23 = e2.cross(f3);
  vec3f ef31 = e3.cross(f1);
  vec3f ef32 = e3.cross(f2);
  vec3f ef33 = e3.cross(f3);

  // now begin the series of tests
  if (!project3(n1, q1, q2, q3)) return false;
  if (!project3(m1, -q1, p2-q1, p3-q1)) return false;

  if (!project6(ef11, p1, p2, p3, q1, q2, q3)) return false;
  if (!project6(ef12, p1, p2, p3, q1, q2, q3)) return false;
  if (!project6(ef13, p1, p2, p3, q1, q2, q3)) return false;
  if (!project6(ef21, p1, p2, p3, q1, q2, q3)) return false;
  if (!project6(ef22, p1, p2, p3, q1, q2, q3)) return false;
  if (!project6(ef23, p1, p2, p3, q1, q2, q3)) return false;
  if (!project6(ef31, p1, p2, p3, q1, q2, q3)) return false;
  if (!project6(ef32, p1, p2, p3, q1, q2, q3)) return false;
  if (!project6(ef33, p1, p2, p3, q1, q2, q3)) return false;
  if (!project6(g1, p1, p2, p3, q1, q2, q3)) return false;
  if (!project6(g2, p1, p2, p3, q1, q2, q3)) return false;
  if (!project6(g3, p1, p2, p3, q1, q2, q3)) return false;
  if (!project6(h1, p1, p2, p3, q1, q2, q3)) return false;
  if (!project6(h2, p1, p2, p3, q1, q2, q3)) return false;
  if (!project6(h3, p1, p2, p3, q1, q2, q3)) return false;

  return true;
}

extern void mesh_id(int id, vector<mesh *> &m, int &mid, int &fid);

extern void pushMesh2GPU(int  numFace, int numVert, void *faces, void *nodes);
extern void updateMesh2GPU(void *nodes,void *prenode);
static tri3f *s_faces;
static vec3f *s_nodes;
static int s_numFace = 0, s_numVert = 0;

void updateMesh2GPU(vector <mesh *> &ms)
{
	vec3f *curVert = s_nodes;

	//rky	
	vec3f *preVert = new vec3f[s_numVert];
	vector<vec3f> tem;
	vec3f *oldcurVert = preVert;
	for (int i = 0; i < ms.size(); i++) {
		mesh *m = ms[i];
		cout << m->_num_vtx << endl;
		memcpy(oldcurVert, m->_ovtxs, sizeof(vec3f)*m->_num_vtx);
		oldcurVert += m->_num_vtx;
	}
    // for (int j = 0; j < s_numVert; ++j) {
    //     cout << preVert[j].x << " "<< preVert[j].y << " " << preVert[j].z << endl;
    // }

	for (int i = 0; i < ms.size(); i++) {
		mesh *m = ms[i];
		memcpy(curVert, m->_vtxs, sizeof(vec3f)*m->_num_vtx);
		curVert += m->_num_vtx;
	}

	updateMesh2GPU(s_nodes, preVert);
	//updateMesh2GPU(s_nodes);
	//delete[] preVert;
}

void pushMesh2GPU(vector<mesh *> &ms)
{
	for (int i = 0; i < ms.size(); i++) {
		s_numFace += ms[i]->_num_tri;
		s_numVert += ms[i]->_num_vtx;
	}

	s_faces = new tri3f[s_numFace];
	s_nodes = new vec3f[s_numVert];

	int curFace = 0;
	int vertCount = 0;
	vec3f *curVert = s_nodes;
	for (int i = 0; i < ms.size(); i++) {
		mesh *m = ms[i];
		for (int j = 0; j < m->_num_tri; j++) {
			tri3f &t = m->_tris[j];
			s_faces[curFace++] = tri3f(t.id0() + vertCount, t.id1() + vertCount, t.id2() + vertCount);
		}
		vertCount += m->_num_vtx;

		memcpy(curVert, m->_vtxs, sizeof(vec3f)*m->_num_vtx);
		curVert += m->_num_vtx;
	}

	pushMesh2GPU(s_numFace, s_numVert, s_faces, s_nodes);
}

extern int getCollisionsGPU(int *,double *);
extern void initGPU();

//rky
extern void init(int num);
extern void front_clear();

bool cmp(vector<tri_pair> a, vector<tri_pair> b){
	unsigned int ta[4],tb[4];
	for (int i = 0; i < 2; i++){
		a[i].get(ta[i * 2],ta[i * 2 + 1]);
		b[i].get(tb[i * 2], tb[i * 2 + 1]);
	}
	if (ta[0] != tb[0])
		return ta[0] < tb[0];
	else if (ta[2] != tb[2])
		return ta[2] < tb[2];
	else if (ta[1] != tb[1])
		return ta[1] < tb[1];
	else
		return ta[3] < tb[3];
}

//static vector<bvh*> bvhC;
//static std::vector<mesh *> meshes;


void body_collide_gpu(vector<mesh_pair> mpair, vector<CollisionDate> bodys, vector<vector<tri_pair>> &contacts,vector<double>& contact_time){
	static vector<bvh*> bvhC;
	//front_list fIntra;
	static std::vector<mesh *> meshes;

	//FILE* f = fopen("bvh.txt", "a");

	vector<tri_pair> fret;
	static vector<int> _tri_offset;
	if (bvhC.size() == 0) {

		double tt= omp_get_wtime();

		for (int i = 0; i < bodys.size(); i++) {
			meshes.push_back(bodys[i].ms);
			_tri_offset.push_back(i == 0 ? bodys[0].ms->_num_tri : (_tri_offset[i - 1] + bodys[i].ms->_num_tri));
			bvh* tem = new bvh(meshes, i);
			bvhC.push_back(tem);
			//for (int st = 0; st < tem->get_nums() * 2 - 1; st++)
			//{
				//fprintf(f, "%d\n", (tem->root() + st)->triID());
			//}
		}
		double ttt = omp_get_wtime();
		printf(" BVH shijian   %f\n", ttt - tt);

		initGPU();
		//rky
		init(bvhC.size());
		pushMesh2GPU(meshes);
		
		

		for (int i = 0; i < bvhC.size(); i++) {
			bvhC[i]->push2GPU(i, i == 0 ? 0 : _tri_offset[i - 1]);
		}
		
		


		int mesh_1, mesh_2;
		for (int i = 0; i < mpair.size(); i++)
		{
			mesh_1 = mpair[i].first;
			mesh_2 = mpair[i].second;
			front_list tem_fl;

			if (mesh_2 != -1)
			{
				//fret.clear();
				bvhC[mesh_1]->collide(bvhC[mesh_2], tem_fl);
				tem_fl.push2GPU(bvhC[mesh_1]->root(), bvhC[mesh_2]->root(), mesh_1, mesh_2);
			}

		}
	}//if

	//fclose(f);

	updateMesh2GPU(meshes);

#ifdef FOR_DEBUG
	vec3f *pts = meshes[0]->getVtxs() + 3126;
	printf("XXXXXXXXXXXX3126: %lf, %lf, %lf\n", pts->x, pts->y, pts->z);
#endif

	int *buffer = new int[10240 * 2];
	double *time_buffer = new double[10240 * 2];

	int count = getCollisionsGPU(buffer, time_buffer);
	printf("Collision num %d\n",count);

	tri_pair *pairs = (tri_pair *)buffer;
	vector<tri_pair> ret(pairs, pairs + count);

	//Find mesh id and face id
	for (int i = 0; i < count; i++){
		vector<tri_pair> tem;
		int mid1, mid2;
		unsigned int fid1, fid2;
		ret[i].get(fid1, fid2);

		for (int j = 0; j < _tri_offset.size(); j++){
			if (fid1 < _tri_offset[j]){
			    //if (fid1 <= _tri_offset[j]){
				mid1 = j == 0 ? 0 : j;
				break;
			}			
		}
		tem.push_back(tri_pair(mid1, fid1 - (mid1 == 0 ? 0 : _tri_offset[mid1 - 1])));
        //tem.push_back(tri_pair(mid1, fid1 - 1-(mid1 == 0 ? 0 : _tri_offset[mid1 - 1])));
		for (int j = 0; j < _tri_offset.size(); j++){
			if (fid2 < _tri_offset[j]){
                //if (fid2 <= _tri_offset[j]){
				mid2 = j == 0 ? 0 : j;
				break;
			}
		}

		tem.push_back(tri_pair(mid2, fid2 - (mid2 == 0 ? 0 : _tri_offset[mid2 - 1])));
        //tem.push_back(tri_pair(mid2, fid2 -1- (mid2 == 0 ? 0 : _tri_offset[mid2 - 1])));
		contacts.push_back(tem);
		contact_time.push_back(time_buffer[i]);
	}

	//std::sort(contacts.begin(), contacts.end(), cmp);
	
	delete[] buffer;
	delete[] time_buffer;
}
