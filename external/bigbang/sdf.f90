SUBROUTINE torus_sdf(val, X, C, N, r1, r2) 
IMPLICIT NONE 
REAL(KIND=8) val(1, 1) 
REAL(KIND=8) X(3, 1) 
REAL(KIND=8) C(3, 1) 
REAL(KIND=8) N(3, 1) 
REAL(KIND=8) r1(1, 1) 
REAL(KIND=8) r2(1, 1) 
REAL(KIND=8)  tt1 
REAL(KIND=8)  tt2 
REAL(KIND=8)  tt3 
REAL(KIND=8)  tt4 
REAL(KIND=8)  tt5 
REAL(KIND=8)  tt6 
REAL(KIND=8)  tt7 
REAL(KIND=8)  tt8 
REAL(KIND=8)  tt9 
REAL(KIND=8)  tt10 
tt1 = -r1(1,1)
tt2 = -C(1,1)
tt3 = tt1-r2(1,1)
tt4 = (-N(3,1)*X(3,1))+C(3,1)*N(3,1)-N(2,1)*X(2,1)+C(2,1)*N(2,1)-&
&N(1,1)*X(1,1)+C(1,1)*N(1,1)
tt5 = N(1,1)*tt4+X(1,1)+tt2
tt6 = -C(2,1)
tt7 = N(2,1)*tt4+X(2,1)+tt6
tt8 = -C(3,1)
tt9 = N(3,1)*tt4+X(3,1)+tt8
tt10 = 1/sqrt(tt9**2+tt7**2+tt5**2)
val(1,1) = sqrt(((tt3*tt9*tt10)/2.0E+0+X(3,1)+tt8)**2+((tt3*tt7*t&
&t10)/2.0E+0+X(2,1)+tt6)**2+((tt3*tt5*tt10)/2.0E+0+X(1,1)+tt2)**2)&
&-(r2(1,1)+tt1)/2.0E+0
END 
SUBROUTINE torus_sdf_jac(jac, X, C, N, r1, r2) 
IMPLICIT NONE 
REAL(KIND=8) jac(1, 3) 
REAL(KIND=8) X(3, 1) 
REAL(KIND=8) C(3, 1) 
REAL(KIND=8) N(3, 1) 
REAL(KIND=8) r1(1, 1) 
REAL(KIND=8) r2(1, 1) 
REAL(KIND=8)  tt1 
REAL(KIND=8)  tt2 
REAL(KIND=8)  tt3 
REAL(KIND=8)  tt4 
REAL(KIND=8)  tt5 
REAL(KIND=8)  tt6 
REAL(KIND=8)  tt7 
REAL(KIND=8)  tt8 
REAL(KIND=8)  tt9 
REAL(KIND=8)  tt10 
REAL(KIND=8)  tt11 
REAL(KIND=8)  tt12 
REAL(KIND=8)  tt13 
REAL(KIND=8)  tt14 
REAL(KIND=8)  tt15 
REAL(KIND=8)  tt16 
REAL(KIND=8)  tt17 
REAL(KIND=8)  tt18 
REAL(KIND=8)  tt19 
REAL(KIND=8)  tt20 
REAL(KIND=8)  tt21 
REAL(KIND=8)  tt22 
REAL(KIND=8)  tt23 
REAL(KIND=8)  tt24 
tt1 = (-r2(1,1))-r1(1,1)
tt2 = -C(1,1)
tt3 = (-N(3,1)*X(3,1))+C(3,1)*N(3,1)-N(2,1)*X(2,1)+C(2,1)*N(2,1)-&
&N(1,1)*X(1,1)+C(1,1)*N(1,1)
tt4 = N(1,1)*tt3+X(1,1)+tt2
tt5 = 1-N(1,1)**2
tt6 = -C(2,1)
tt7 = N(2,1)*tt3+X(2,1)+tt6
tt8 = -C(3,1)
tt9 = N(3,1)*tt3+X(3,1)+tt8
tt10 = (-2*N(1,1)*N(3,1)*tt9)-2*N(1,1)*N(2,1)*tt7+2*tt5*tt4
tt11 = sqrt(tt9**2+tt7**2+tt4**2)
tt12 = 1/tt11**3
tt13 = 1/tt11
tt14 = (tt1*tt4*tt13)/2.0E+0+X(1,1)+tt2
tt15 = -(N(1,1)*tt1*N(2,1)*tt13)/2.0E+0
tt16 = (tt1*tt7*tt13)/2.0E+0+X(2,1)+tt6
tt17 = -(N(1,1)*tt1*N(3,1)*tt13)/2.0E+0
tt18 = (tt1*tt9*tt13)/2.0E+0+X(3,1)+tt8
tt19 = 1/sqrt(tt18**2+tt16**2+tt14**2)
tt20 = 1-N(2,1)**2
tt21 = (-2*N(2,1)*N(3,1)*tt9)+2*tt20*tt7-2*N(1,1)*N(2,1)*tt4
tt22 = -(tt1*N(2,1)*N(3,1)*tt13)/2.0E+0
tt23 = 1-N(3,1)**2
tt24 = 2*tt23*tt9-2*N(2,1)*N(3,1)*tt7-2*N(1,1)*N(3,1)*tt4
jac(1,1) = ((2*(tt17-(tt1*tt9*tt10*tt12)/4.0E+0)*tt18+2*(tt15-(tt&
&1*tt7*tt10*tt12)/4.0E+0)*tt16+2*((tt5*tt1*tt13)/2.0E+0-(tt1*tt4*t&
&t10*tt12)/4.0E+0+1)*tt14)*tt19)/2.0E+0
jac(1,2) = ((2*(tt22-(tt1*tt9*tt21*tt12)/4.0E+0)*tt18+2*((tt1*tt2&
&0*tt13)/2.0E+0-(tt1*tt7*tt21*tt12)/4.0E+0+1)*tt16+2*(tt15-(tt1*tt&
&4*tt21*tt12)/4.0E+0)*tt14)*tt19)/2.0E+0
jac(1,3) = ((2*((tt1*tt23*tt13)/2.0E+0-(tt1*tt9*tt24*tt12)/4.0E+0&
&+1)*tt18+2*(tt22-(tt1*tt7*tt24*tt12)/4.0E+0)*tt16+2*(tt17-(tt1*tt&
&4*tt24*tt12)/4.0E+0)*tt14)*tt19)/2.0E+0
END 
SUBROUTINE torus_sdf_hes(hes, X, C, N, r1, r2) 
IMPLICIT NONE 
REAL(KIND=8) hes(3, 3) 
REAL(KIND=8) X(3, 1) 
REAL(KIND=8) C(3, 1) 
REAL(KIND=8) N(3, 1) 
REAL(KIND=8) r1(1, 1) 
REAL(KIND=8) r2(1, 1) 
REAL(KIND=8)  tt1 
REAL(KIND=8)  tt2 
REAL(KIND=8)  tt3 
REAL(KIND=8)  tt4 
REAL(KIND=8)  tt5 
REAL(KIND=8)  tt6 
REAL(KIND=8)  tt7 
REAL(KIND=8)  tt8 
REAL(KIND=8)  tt9 
REAL(KIND=8)  tt10 
REAL(KIND=8)  tt11 
REAL(KIND=8)  tt12 
REAL(KIND=8)  tt13 
REAL(KIND=8)  tt14 
REAL(KIND=8)  tt15 
REAL(KIND=8)  tt16 
REAL(KIND=8)  tt17 
REAL(KIND=8)  tt18 
REAL(KIND=8)  tt19 
REAL(KIND=8)  tt20 
REAL(KIND=8)  tt21 
REAL(KIND=8)  tt22 
REAL(KIND=8)  tt23 
REAL(KIND=8)  tt24 
REAL(KIND=8)  tt25 
REAL(KIND=8)  tt26 
REAL(KIND=8)  tt27 
REAL(KIND=8)  tt28 
REAL(KIND=8)  tt29 
REAL(KIND=8)  tt30 
REAL(KIND=8)  tt31 
REAL(KIND=8)  tt32 
REAL(KIND=8)  tt33 
REAL(KIND=8)  tt34 
REAL(KIND=8)  tt35 
REAL(KIND=8)  tt36 
REAL(KIND=8)  tt37 
REAL(KIND=8)  tt38 
REAL(KIND=8)  tt39 
REAL(KIND=8)  tt40 
REAL(KIND=8)  tt41 
REAL(KIND=8)  tt42 
REAL(KIND=8)  tt43 
REAL(KIND=8)  tt44 
REAL(KIND=8)  tt45 
REAL(KIND=8)  tt46 
REAL(KIND=8)  tt47 
REAL(KIND=8)  tt48 
REAL(KIND=8)  tt49 
REAL(KIND=8)  tt50 
REAL(KIND=8)  tt51 
REAL(KIND=8)  tt52 
REAL(KIND=8)  tt53 
REAL(KIND=8)  tt54 
REAL(KIND=8)  tt55 
REAL(KIND=8)  tt56 
REAL(KIND=8)  tt57 
REAL(KIND=8)  tt58 
REAL(KIND=8)  tt59 
REAL(KIND=8)  tt60 
tt1 = (-r2(1,1))-r1(1,1)
tt2 = -C(1,1)
tt3 = (-N(3,1)*X(3,1))+C(3,1)*N(3,1)-N(2,1)*X(2,1)+C(2,1)*N(2,1)-&
&N(1,1)*X(1,1)+C(1,1)*N(1,1)
tt4 = N(1,1)*tt3+X(1,1)+tt2
tt5 = N(1,1)**2
tt6 = 1-tt5
tt7 = -C(2,1)
tt8 = N(2,1)*tt3+X(2,1)+tt7
tt9 = -C(3,1)
tt10 = N(3,1)*tt3+X(3,1)+tt9
tt11 = (-2*N(1,1)*N(3,1)*tt10)-2*N(1,1)*N(2,1)*tt8+2*tt6*tt4
tt12 = sqrt(tt10**2+tt8**2+tt4**2)
tt13 = 1/tt12**3
tt14 = 1/tt12
tt15 = (tt6*tt1*tt14)/2.0E+0-(tt1*tt4*tt11*tt13)/4.0E+0+1
tt16 = (tt1*tt4*tt14)/2.0E+0+X(1,1)+tt2
tt17 = -(N(1,1)*tt1*N(2,1)*tt14)/2.0E+0
tt18 = tt17-(tt1*tt8*tt11*tt13)/4.0E+0
tt19 = (tt1*tt8*tt14)/2.0E+0+X(2,1)+tt7
tt20 = -(N(1,1)*tt1*N(3,1)*tt14)/2.0E+0
tt21 = tt20-(tt1*tt10*tt11*tt13)/4.0E+0
tt22 = (tt1*tt10*tt14)/2.0E+0+X(3,1)+tt9
tt23 = 2*tt21*tt22+2*tt18*tt19+2*tt15*tt16
tt24 = sqrt(tt22**2+tt19**2+tt16**2)
tt25 = 1/tt24**3
tt26 = tt11**2
tt27 = 1/tt12**5
tt28 = N(2,1)**2
tt29 = 2*tt5*tt28
tt30 = N(3,1)**2
tt31 = 2*tt5*tt30
tt32 = tt31+tt29+2*tt6**2
tt33 = 1/tt24
tt34 = 1-tt28
tt35 = (-2*N(2,1)*N(3,1)*tt10)+2*tt34*tt8-2*N(1,1)*N(2,1)*tt4
tt36 = tt17-(tt1*tt4*tt35*tt13)/4.0E+0
tt37 = (tt1*tt34*tt14)/2.0E+0-(tt1*tt8*tt35*tt13)/4.0E+0+1
tt38 = -(tt1*N(2,1)*N(3,1)*tt14)/2.0E+0
tt39 = tt38-(tt1*tt10*tt35*tt13)/4.0E+0
tt40 = 2*tt39*tt22+2*tt37*tt19+2*tt36*tt16
tt41 = 2*N(1,1)*N(2,1)*tt30-2*N(1,1)*N(2,1)*tt34-2*N(1,1)*tt6*N(2&
&,1)
tt42 = (tt1*N(2,1)*N(3,1)*tt11*tt13)/4.0E+0
tt43 = (N(1,1)*tt1*N(3,1)*tt35*tt13)/4.0E+0
tt44 = ((2*(tt43+tt42-(tt1*tt41*tt10*tt13)/4.0E+0+(3.0E+0*tt1*tt1&
&0*tt11*tt35*tt27)/8.0E+0)*tt22+2*((N(1,1)*tt1*N(2,1)*tt35*tt13)/4&
&.0E+0-(tt1*tt34*tt11*tt13)/4.0E+0-(tt1*tt41*tt8*tt13)/4.0E+0+(3.0&
&E+0*tt1*tt8*tt11*tt35*tt27)/8.0E+0)*tt19+2*((-(tt6*tt1*tt35*tt13)&
&/4.0E+0)+(N(1,1)*tt1*N(2,1)*tt11*tt13)/4.0E+0-(tt1*tt41*tt4*tt13)&
&/4.0E+0+(3.0E+0*tt1*tt4*tt11*tt35*tt27)/8.0E+0)*tt16+2*tt21*tt39+&
&2*tt18*tt37+2*tt15*tt36)*tt33)/2.0E+0-(tt23*tt40*tt25)/4.0E+0
tt45 = 1-tt30
tt46 = 2*tt45*tt10-2*N(2,1)*N(3,1)*tt8-2*N(1,1)*N(3,1)*tt4
tt47 = tt20-(tt1*tt4*tt46*tt13)/4.0E+0
tt48 = tt38-(tt1*tt8*tt46*tt13)/4.0E+0
tt49 = (tt1*tt45*tt14)/2.0E+0-(tt1*tt10*tt46*tt13)/4.0E+0+1
tt50 = 2*tt49*tt22+2*tt48*tt19+2*tt47*tt16
tt51 = (-2*N(1,1)*N(3,1)*tt45)+2*N(1,1)*tt28*N(3,1)-2*N(1,1)*tt6*&
&N(3,1)
tt52 = (N(1,1)*tt1*N(2,1)*tt46*tt13)/4.0E+0
tt53 = ((2*((N(1,1)*tt1*N(3,1)*tt46*tt13)/4.0E+0-(tt1*tt45*tt11*t&
&t13)/4.0E+0-(tt1*tt51*tt10*tt13)/4.0E+0+(3.0E+0*tt1*tt10*tt11*tt4&
&6*tt27)/8.0E+0)*tt22+2*(tt52+tt42-(tt1*tt51*tt8*tt13)/4.0E+0+(3.0&
&E+0*tt1*tt8*tt11*tt46*tt27)/8.0E+0)*tt19+2*((-(tt6*tt1*tt46*tt13)&
&/4.0E+0)+(N(1,1)*tt1*N(3,1)*tt11*tt13)/4.0E+0-(tt1*tt51*tt4*tt13)&
&/4.0E+0+(3.0E+0*tt1*tt4*tt11*tt46*tt27)/8.0E+0)*tt16+2*tt21*tt49+&
&2*tt18*tt48+2*tt15*tt47)*tt33)/2.0E+0-(tt23*tt50*tt25)/4.0E+0
tt54 = tt35**2
tt55 = 2*tt28*tt30
tt56 = tt55+2*tt34**2+tt29
tt57 = (-2*N(2,1)*N(3,1)*tt45)-2*N(2,1)*tt34*N(3,1)+2*tt5*N(2,1)*&
&N(3,1)
tt58 = ((2*((tt1*N(2,1)*N(3,1)*tt46*tt13)/4.0E+0-(tt1*tt45*tt35*t&
&t13)/4.0E+0-(tt1*tt57*tt10*tt13)/4.0E+0+(3.0E+0*tt1*tt10*tt35*tt4&
&6*tt27)/8.0E+0)*tt22+2*((-(tt1*tt34*tt46*tt13)/4.0E+0)+(tt1*N(2,1&
&)*N(3,1)*tt35*tt13)/4.0E+0-(tt1*tt57*tt8*tt13)/4.0E+0+(3.0E+0*tt1&
&*tt8*tt35*tt46*tt27)/8.0E+0)*tt19+2*(tt52+tt43-(tt1*tt57*tt4*tt13&
&)/4.0E+0+(3.0E+0*tt1*tt4*tt35*tt46*tt27)/8.0E+0)*tt16+2*tt39*tt49&
&+2*tt37*tt48+2*tt36*tt47)*tt33)/2.0E+0-(tt40*tt50*tt25)/4.0E+0
tt59 = tt46**2
tt60 = 2*tt45**2+tt55+tt31
hes(1,1) = ((2*tt21**2+2*tt18**2+2*tt15**2+2*((N(1,1)*tt1*N(3,1)*&
&tt11*tt13)/2.0E+0-(tt1*tt32*tt10*tt13)/4.0E+0+(3.0E+0*tt1*tt10*tt&
&26*tt27)/8.0E+0)*tt22+2*((N(1,1)*tt1*N(2,1)*tt11*tt13)/2.0E+0-(tt&
&1*tt32*tt8*tt13)/4.0E+0+(3.0E+0*tt1*tt8*tt26*tt27)/8.0E+0)*tt19+2&
&*((-(tt6*tt1*tt11*tt13)/2.0E+0)-(tt1*tt32*tt4*tt13)/4.0E+0+(3.0E+&
&0*tt1*tt4*tt26*tt27)/8.0E+0)*tt16)*tt33)/2.0E+0-(tt23**2*tt25)/4.&
&0E+0
hes(1,2) = tt44
hes(1,3) = tt53
hes(2,1) = tt44
hes(2,2) = ((2*tt39**2+2*tt37**2+2*tt36**2+2*((tt1*N(2,1)*N(3,1)*&
&tt35*tt13)/2.0E+0-(tt1*tt56*tt10*tt13)/4.0E+0+(3.0E+0*tt1*tt10*tt&
&54*tt27)/8.0E+0)*tt22+2*((-(tt1*tt34*tt35*tt13)/2.0E+0)-(tt1*tt56&
&*tt8*tt13)/4.0E+0+(3.0E+0*tt1*tt8*tt54*tt27)/8.0E+0)*tt19+2*((N(1&
&,1)*tt1*N(2,1)*tt35*tt13)/2.0E+0-(tt1*tt56*tt4*tt13)/4.0E+0+(3.0E&
&+0*tt1*tt4*tt54*tt27)/8.0E+0)*tt16)*tt33)/2.0E+0-(tt40**2*tt25)/4&
&.0E+0
hes(2,3) = tt58
hes(3,1) = tt53
hes(3,2) = tt58
hes(3,3) = ((2*tt49**2+2*tt48**2+2*tt47**2+2*((-(tt1*tt45*tt46*tt&
&13)/2.0E+0)-(tt1*tt60*tt10*tt13)/4.0E+0+(3.0E+0*tt1*tt10*tt59*tt2&
&7)/8.0E+0)*tt22+2*((tt1*N(2,1)*N(3,1)*tt46*tt13)/2.0E+0-(tt1*tt60&
&*tt8*tt13)/4.0E+0+(3.0E+0*tt1*tt8*tt59*tt27)/8.0E+0)*tt19+2*((N(1&
&,1)*tt1*N(3,1)*tt46*tt13)/2.0E+0-(tt1*tt60*tt4*tt13)/4.0E+0+(3.0E&
&+0*tt1*tt4*tt59*tt27)/8.0E+0)*tt16)*tt33)/2.0E+0-(tt50**2*tt25)/4&
&.0E+0
END 
SUBROUTINE cylinder_sdf(val, X, C, N, r1) 
IMPLICIT NONE 
REAL(KIND=8) val(1, 1) 
REAL(KIND=8) X(3, 1) 
REAL(KIND=8) C(3, 1) 
REAL(KIND=8) N(3, 1) 
REAL(KIND=8) r1(1, 1) 
REAL(KIND=8)  tt1 
tt1 = (-N(3,1)*X(3,1))+C(3,1)*N(3,1)-N(2,1)*X(2,1)+C(2,1)*N(2,1)-&
&N(1,1)*X(1,1)+C(1,1)*N(1,1)
val(1,1) = sqrt((N(3,1)*tt1+X(3,1)-C(3,1))**2+(N(2,1)*tt1+X(2,1)-&
&C(2,1))**2+(N(1,1)*tt1+X(1,1)-C(1,1))**2)-r1(1,1)
END 
SUBROUTINE cylinder_sdf_jac(jac, X, C, N, r1) 
IMPLICIT NONE 
REAL(KIND=8) jac(1, 3) 
REAL(KIND=8) X(3, 1) 
REAL(KIND=8) C(3, 1) 
REAL(KIND=8) N(3, 1) 
REAL(KIND=8) r1(1, 1) 
REAL(KIND=8)  tt1 
REAL(KIND=8)  tt2 
REAL(KIND=8)  tt3 
REAL(KIND=8)  tt4 
REAL(KIND=8)  tt5 
tt1 = (-N(3,1)*X(3,1))+C(3,1)*N(3,1)-N(2,1)*X(2,1)+C(2,1)*N(2,1)-&
&N(1,1)*X(1,1)+C(1,1)*N(1,1)
tt2 = N(1,1)*tt1+X(1,1)-C(1,1)
tt3 = N(2,1)*tt1+X(2,1)-C(2,1)
tt4 = N(3,1)*tt1+X(3,1)-C(3,1)
tt5 = 1/sqrt(tt4**2+tt3**2+tt2**2)
jac(1,1) = (((-2*N(1,1)*N(3,1)*tt4)-2*N(1,1)*N(2,1)*tt3+2*(1-N(1,&
&1)**2)*tt2)*tt5)/2.0E+0
jac(1,2) = (((-2*N(2,1)*N(3,1)*tt4)+2*(1-N(2,1)**2)*tt3-2*N(1,1)*&
&N(2,1)*tt2)*tt5)/2.0E+0
jac(1,3) = ((2*(1-N(3,1)**2)*tt4-2*N(2,1)*N(3,1)*tt3-2*N(1,1)*N(3&
&,1)*tt2)*tt5)/2.0E+0
END 
SUBROUTINE cylinder_sdf_hes(hes, X, C, N, r1) 
IMPLICIT NONE 
REAL(KIND=8) hes(3, 3) 
REAL(KIND=8) X(3, 1) 
REAL(KIND=8) C(3, 1) 
REAL(KIND=8) N(3, 1) 
REAL(KIND=8) r1(1, 1) 
REAL(KIND=8)  tt1 
REAL(KIND=8)  tt2 
REAL(KIND=8)  tt3 
REAL(KIND=8)  tt4 
REAL(KIND=8)  tt5 
REAL(KIND=8)  tt6 
REAL(KIND=8)  tt7 
REAL(KIND=8)  tt8 
REAL(KIND=8)  tt9 
REAL(KIND=8)  tt10 
REAL(KIND=8)  tt11 
REAL(KIND=8)  tt12 
REAL(KIND=8)  tt13 
REAL(KIND=8)  tt14 
REAL(KIND=8)  tt15 
REAL(KIND=8)  tt16 
REAL(KIND=8)  tt17 
REAL(KIND=8)  tt18 
REAL(KIND=8)  tt19 
REAL(KIND=8)  tt20 
REAL(KIND=8)  tt21 
REAL(KIND=8)  tt22 
tt1 = N(1,1)**2
tt2 = 1-tt1
tt3 = (-N(3,1)*X(3,1))+C(3,1)*N(3,1)-N(2,1)*X(2,1)+C(2,1)*N(2,1)-&
&N(1,1)*X(1,1)+C(1,1)*N(1,1)
tt4 = N(1,1)*tt3+X(1,1)-C(1,1)
tt5 = N(2,1)*tt3+X(2,1)-C(2,1)
tt6 = N(3,1)*tt3+X(3,1)-C(3,1)
tt7 = (-2*N(1,1)*N(3,1)*tt6)-2*N(1,1)*N(2,1)*tt5+2*tt2*tt4
tt8 = sqrt(tt6**2+tt5**2+tt4**2)
tt9 = 1/tt8**3
tt10 = N(2,1)**2
tt11 = 2*tt1*tt10
tt12 = N(3,1)**2
tt13 = 2*tt1*tt12
tt14 = 1/tt8
tt15 = 1-tt10
tt16 = (-2*N(2,1)*N(3,1)*tt6)+2*tt15*tt5-2*N(1,1)*N(2,1)*tt4
tt17 = ((2*N(1,1)*N(2,1)*tt12-2*N(1,1)*N(2,1)*tt15-2*N(1,1)*tt2*N&
&(2,1))*tt14)/2.0E+0-(tt7*tt16*tt9)/4.0E+0
tt18 = 1-tt12
tt19 = 2*tt18*tt6-2*N(2,1)*N(3,1)*tt5-2*N(1,1)*N(3,1)*tt4
tt20 = (((-2*N(1,1)*N(3,1)*tt18)+2*N(1,1)*tt10*N(3,1)-2*N(1,1)*tt&
&2*N(3,1))*tt14)/2.0E+0-(tt7*tt19*tt9)/4.0E+0
tt21 = 2*tt10*tt12
tt22 = (((-2*N(2,1)*N(3,1)*tt18)-2*N(2,1)*tt15*N(3,1)+2*tt1*N(2,1&
&)*N(3,1))*tt14)/2.0E+0-(tt16*tt19*tt9)/4.0E+0
hes(1,1) = ((tt13+tt11+2*tt2**2)*tt14)/2.0E+0-(tt7**2*tt9)/4.0E+0&
&
hes(1,2) = tt17
hes(1,3) = tt20
hes(2,1) = tt17
hes(2,2) = ((tt21+2*tt15**2+tt11)*tt14)/2.0E+0-(tt16**2*tt9)/4.0E&
&+0
hes(2,3) = tt22
hes(3,1) = tt20
hes(3,2) = tt22
hes(3,3) = ((2*tt18**2+tt21+tt13)*tt14)/2.0E+0-(tt19**2*tt9)/4.0E&
&+0
END 
