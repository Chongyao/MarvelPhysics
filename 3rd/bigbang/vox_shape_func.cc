#include <algorithm>

namespace bigbang {

void vox_SF(double *val, const double *epsilon) {
  const double 
      tt1 = 1-epsilon[0],
      tt2 = 1-epsilon[1],
      tt3 = 1-epsilon[2],
      tt4 = epsilon[0]+1,
      tt5 = epsilon[1]+1,
      tt6 = epsilon[2]+1;
  
  val[0] = (tt1*tt2*tt3)/8.0;
  val[1] = (tt4*tt2*tt3)/8.0;
  val[2] = (tt1*tt5*tt3)/8.0;
  val[3] = (tt4*tt5*tt3)/8.0;
  val[4] = (tt1*tt2*tt6)/8.0;
  val[5] = (tt4*tt2*tt6)/8.0;
  val[6] = (tt1*tt5*tt6)/8.0;
  val[7] = (tt4*tt5*tt6)/8.0;
}


void vox_SF_jac(double *jac, const double *epsilon) {
  const double
      tt1 = 1-epsilon[1],
      tt2 = 1-epsilon[2],
      tt3 = epsilon[1]+1,
      tt4 = epsilon[2]+1,
      tt5 = 1-epsilon[0],
      tt6 = epsilon[0]+1;

  jac[0] = -(tt1*tt2)/8.0;
  jac[1] = (tt1*tt2)/8.0;
  jac[2] = -(tt3*tt2)/8.0;
  jac[3] = (tt3*tt2)/8.0;
  jac[4] = -(tt1*tt4)/8.0;
  jac[5] = (tt1*tt4)/8.0;
  jac[6] = -(tt3*tt4)/8.0;
  jac[7] = (tt3*tt4)/8.0;
  jac[8] = -(tt5*tt2)/8.0;
  jac[9] = -(tt6*tt2)/8.0;
  jac[10] = (tt5*tt2)/8.0;
  jac[11] = (tt6*tt2)/8.0;
  jac[12] = -(tt5*tt4)/8.0;
  jac[13] = -(tt6*tt4)/8.0;
  jac[14] = (tt5*tt4)/8.0;
  jac[15] = (tt6*tt4)/8.0;
  jac[16] = -(tt5*tt1)/8.0;
  jac[17] = -(tt6*tt1)/8.0;
  jac[18] = -(tt5*tt3)/8.0;
  jac[19] = -(tt6*tt3)/8.0;
  jac[20] = (tt5*tt1)/8.0;
  jac[21] = (tt6*tt1)/8.0;
  jac[22] = (tt5*tt3)/8.0;
  jac[23] = (tt6*tt3)/8.0;
}

//        quad4
//
//        3---------2
//        |         |
//        |         |
//        |         |
//        |         |
//        0---------1
//  y
//  |__x
//

void quad_SF(double *val, const double *epsilon) {
  const double
      tt1 = 1-epsilon[0],
      tt2 = 1-epsilon[1],
      tt3 = epsilon[0]+1,
      tt4 = epsilon[1]+1;

  val[0] = (tt1*tt2)/4.0;
  val[1] = (tt3*tt2)/4.0;
  val[2] = (tt3*tt4)/4.0;
  val[3] = (tt1*tt4)/4.0;
}

void quad_SF_jac(double *jac, const double *epsilon) {
  const double
      tt1 = 1-epsilon[1],
      tt2 = epsilon[1]+1,
      tt3 = 1-epsilon[0],
      tt4 = epsilon[0]+1;

  jac[0] = -tt1/4.0;
  jac[1] = tt1/4.0;
  jac[2] = tt2/4.0;
  jac[3] = -tt2/4.0;
  jac[4] = -tt3/4.0;
  jac[5] = -tt4/4.0;
  jac[6] = tt4/4.0;
  jac[7] = tt3/4.0;
}

//        quad9
//
//        3----6----2
//        |    |    |
//        |    |    |
//        7----8----5
//        |    |    |
//        |    |    |
//        0----4----1
//  y
//  |__x
//

void quad9_SF(double *val, const double *epsilon) {
  const double xi = epsilon[0], eta = epsilon[1];
  const double
      tt1 = eta-1,
      tt2 = xi-1,
      tt3 = xi+1,
      tt4 = eta+1,
      tt5 = 1-xi*xi,
      tt6 = 1-eta*eta;

  val[0] = (tt1*eta*tt2*xi)/4.0;
  val[1] = (tt1*eta*xi*tt3)/4.0;
  val[2] = (eta*tt4*xi*tt3)/4.0;
  val[3] = (eta*tt4*tt2*xi)/4.0;
  val[4] = (tt1*eta*tt5)/2.0;
  val[5] = (tt6*xi*tt3)/2.0;
  val[6] = (eta*tt4*tt5)/2.0;
  val[7] = (tt6*tt2*xi)/2.0;
  val[8] = tt6*tt5;
}

void quad9_SF_jac(double *jac, const double *epsilon) {
  const double xi = epsilon[0], eta = epsilon[1];
  const double
      tt1 = eta-1,
      tt2 = xi-1,
      tt3 = (tt1*eta*xi)/4.0,
      tt4 = xi+1,
      tt5=eta+1,
      tt6=(eta*tt5*xi)/4.0,
      tt7=1-eta*eta,
      tt8=(tt7*xi)/2.0,
      tt9=(eta*tt2*xi)/4.0,
      tt10=(eta*xi*tt4)/4.0,
      tt11=1-xi*xi,
      tt12=(eta*tt11)/2.0;

  const double rtn[] = {tt3+(tt1*eta*tt2)/4.0,(tt1*eta*tt4)/4.0+tt3,(eta*tt5*tt4)/4.0+tt6,tt6+(eta*tt5*tt2)/4.0,-tt1*eta*xi,(tt7*tt4)/2.0+tt8,-eta*tt5*xi,tt8+(tt7*tt2)/2.0,-2*tt7*xi,tt9+(tt1*tt2*xi)/4.0,tt10+(tt1*xi*tt4)/4.0,(tt5*xi*tt4)/4.0+tt10,(tt5*tt2*xi)/4.0+tt9,tt12+(tt1*tt11)/2.0,-eta*xi*tt4,(tt5*tt11)/2.0+tt12,-eta*tt2*xi,-2*eta*tt11};
  std::copy(rtn, rtn+18, jac);
}

void quad9_BL_jac(double *jac, const double *epsilon) {
  const double u = epsilon[0], v = epsilon[1];
  const double
      tt1=u-1,
      tt2=tt1/2.0,
      tt3=tt2-1,
      tt4=v-1,
      tt5=tt4/2.0,
      tt6=tt5-1,
      tt7=(tt1*tt6*tt4)/16.0,
      tt8=tt2+1,
      tt9=tt5+1,
      tt10=(tt1*tt9*tt4)/16.0,
      tt11=1-tt4*tt4/4.0,
      tt12=(tt1*tt11)/4.0,
      tt13=(tt3*tt1*tt4)/16.0,
      tt14=(tt8*tt1*tt4)/16.0,
      tt15=1-tt1*tt1/4.0,
      tt16=(tt15*tt4)/4.0;
  const double rtn[] = {tt7+(tt3*tt6*tt4)/8.0,tt7+(tt8*tt6*tt4)/8.0,tt10+(tt8*tt9*tt4)/8.0,tt10+(tt3*tt9*tt4)/8.0,-(tt1*tt6*tt4)/4.0,tt12+(tt8*tt11)/2.0,-(tt1*tt9*tt4)/4.0,tt12+(tt3*tt11)/2.0,-tt1*tt11,tt13+(tt3*tt1*tt6)/8.0,tt14+(tt8*tt1*tt6)/8.0,tt14+(tt8*tt1*tt9)/8.0,tt13+(tt3*tt1*tt9)/8.0,tt16+(tt15*tt6)/2.0,-(tt8*tt1*tt4)/4.0,tt16+(tt15*tt9)/2.0,-(tt3*tt1*tt4)/4.0,-tt15*tt4};
  std::copy(rtn, rtn+18, jac);
}

void quad9_BR_jac(double *jac, const double *epsilon) {
  const double u = epsilon[0], v = epsilon[1];
  const double
      tt1=u+1,
      tt2=v-1,
      tt3=tt2/2.0,
      tt4=tt3-1,
      tt5=(tt1*tt4*tt2)/16.0,
      tt6=tt1/2.0,
      tt7=tt6-1,
      tt8=tt6+1,
      tt9=tt3+1,
      tt10=(tt1*tt9*tt2)/16.0,
      tt11=1-tt2*tt2/4.0,
      tt12=(tt1*tt11)/4.0,
      tt13=(tt1*tt7*tt2)/16.0,
      tt14=(tt1*tt8*tt2)/16.0,
      tt15=1-tt1*tt1/4.0,
      tt16=(tt15*tt2)/4.0;
  const double rtn[] = {(tt7*tt4*tt2)/8.0+tt5,(tt8*tt4*tt2)/8.0+tt5,(tt8*tt9*tt2)/8.0+tt10,(tt7*tt9*tt2)/8.0+tt10,-(tt1*tt4*tt2)/4.0,(tt8*tt11)/2.0+tt12,-(tt1*tt9*tt2)/4.0,(tt7*tt11)/2.0+tt12,-tt1*tt11,tt13+(tt1*tt7*tt4)/8.0,tt14+(tt1*tt8*tt4)/8.0,tt14+(tt1*tt8*tt9)/8.0,tt13+(tt1*tt7*tt9)/8.0,tt16+(tt15*tt4)/2.0,-(tt1*tt8*tt2)/4.0,tt16+(tt15*tt9)/2.0,-(tt1*tt7*tt2)/4.0,-tt15*tt2};
  std::copy(rtn, rtn+18, jac);
}

void quad9_TR_jac(double *jac, const double *epsilon) {
  const double u = epsilon[0], v = epsilon[1];
  const double
      tt1=u+1,
      tt2=v+1,
      tt3=tt2/2.0,
      tt4=tt3-1,
      tt5=(tt1*tt2*tt4)/16.0,
      tt6=tt1/2.0,
      tt7=tt6-1,
      tt8=tt6+1,
      tt9=tt3+1,
      tt10=(tt1*tt2*tt9)/16.0,
      tt11=1-tt2*tt2/4.0,
      tt12=(tt1*tt11)/4.0,
      tt13=(tt1*tt7*tt2)/16.0,
      tt14=(tt1*tt8*tt2)/16.0,
      tt15=1-tt1*tt1/4.0,
      tt16=(tt15*tt2)/4.0;
  const double rtn[] = {(tt7*tt2*tt4)/8.0+tt5,(tt8*tt2*tt4)/8.0+tt5,(tt8*tt2*tt9)/8.0+tt10,(tt7*tt2*tt9)/8.0+tt10,-(tt1*tt2*tt4)/4.0,(tt8*tt11)/2.0+tt12,-(tt1*tt2*tt9)/4.0,(tt7*tt11)/2.0+tt12,-tt1*tt11,(tt1*tt7*tt4)/8.0+tt13,(tt1*tt8*tt4)/8.0+tt14,(tt1*tt8*tt9)/8.0+tt14,(tt1*tt7*tt9)/8.0+tt13,(tt15*tt4)/2.0+tt16,-(tt1*tt8*tt2)/4.0,(tt15*tt9)/2.0+tt16,-(tt1*tt7*tt2)/4.0,-tt15*tt2};
  std::copy(rtn, rtn+18, jac);
}

void quad9_TL_jac(double *jac, const double *epsilon) {
  const double u = epsilon[0], v = epsilon[1];
  const double
      tt1=u-1,
      tt2=tt1/2.0,
      tt3=tt2-1,
      tt4=v+1,
      tt5=tt4/2.0,
      tt6=tt5-1,
      tt7=(tt1*tt4*tt6)/16.0,
      tt8=tt2+1,
      tt9=tt5+1,
      tt10=(tt1*tt4*tt9)/16.0,
      tt11=1-tt4*tt4/4.0,
      tt12=(tt1*tt11)/4.0,
      tt13=(tt3*tt1*tt4)/16.0,
      tt14=(tt8*tt1*tt4)/16.0,
      tt15=1-tt1*tt1/4.0,
      tt16=(tt15*tt4)/4.0;
  const double rtn[] = {tt7+(tt3*tt4*tt6)/8.0,tt7+(tt8*tt4*tt6)/8.0,tt10+(tt8*tt4*tt9)/8.0,tt10+(tt3*tt4*tt9)/8.0,-(tt1*tt4*tt6)/4.0,tt12+(tt8*tt11)/2.0,-(tt1*tt4*tt9)/4.0,tt12+(tt3*tt11)/2.0,-tt1*tt11,(tt3*tt1*tt6)/8.0+tt13,(tt8*tt1*tt6)/8.0+tt14,(tt8*tt1*tt9)/8.0+tt14,(tt3*tt1*tt9)/8.0+tt13,(tt15*tt6)/2.0+tt16,-(tt8*tt1*tt4)/4.0,(tt15*tt9)/2.0+tt16,-(tt3*tt1*tt4)/4.0,-tt15*tt4};
  std::copy(rtn, rtn+18, jac);
}

}
