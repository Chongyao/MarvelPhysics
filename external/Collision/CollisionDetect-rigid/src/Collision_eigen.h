#ifndef SELFCOLLISION_EIGEN
#define SELFCOLLISION_EIGEN

#include "cmesh.h"
#include "tmbvh.hpp"
#include "CollisionDate.h"

#include "iostream"
#include<vector>

typedef pair<int, int> mesh_pair;

class Collision_zcy{
 public:

  ~Collision_zcy(){
    for (int i = 0; i < dl_mesh.size(); i++) {
      delete dl_mesh[i];
    }
  }

  //µ÷ÓÃº¯Êý½Ó¿Ú£¬µ÷ÓÃÄÚ²¿Åö×²¼ì²âËã·¨
  void Collid();

  //Åö×²¶ÔÊäÈë½Ó¿Ú
  void Transform_Pair(unsigned int a,unsigned int b);
	
  //Ä£ÐÍÍø¸ñÊäÈë½Ó¿Ú£¬ÊäÈëÄ£ÐÍÍø¸ñµÄµã¼¯ºÍÃæ¼¯
  void Transform_Mesh(unsigned int numVtx, unsigned int numTri, const unsigned int* tris, 
                      const double* vtxs,
                      const double* pre_vtxs,
                      int m_id,bool able_selfcollision=false
                      );

  //Êä³ö½Ó¿Ú£¬·µ»Ø·¢ÉúÅö×²µÄÄ£ÐÍÍø¸ñºÍÈý½ÇÐÎÃæÆ¬µÄ¼¯ºÏ
  vector<vector<tri_pair>> getContactPairs(){ return contact_pairs; }

  //Êä³ö½Ó¿Ú£¬·µ»Ø·¢ÉúÅö×²µÄÅö×²¶ÔÊýÁ¿
  int getNumContacts(){ return contact_pairs.size(); }
	
  //Êä³ö½Ó¿Ú£¬·µ»ØÅö×²¶Ô·¢ÉúÅö×²µÄÊ±¼ä
  vector<double> getContactTimes() { return contact_time; }

  static Collision_zcy* getInstance()			
  {
    if (instance == NULL) {
      instance = new Collision_zcy();
      return instance;
    }
    else
      return instance;
  }


  static Collision_zcy* instance;


 private:
  Collision_zcy():_is_first(true){}
  vector<CollisionDate> bodys;     
  vector<mesh_pair> mesh_pairs;
  vector<vector<tri_pair>> contact_pairs;
  vector<mesh*> dl_mesh;//delete mesh points
  vector<double> contact_time;



  bool _is_first;//ÊÇ·ñµÚÒ»´Î´«ÈëÍø¸ñ

};


#endif
