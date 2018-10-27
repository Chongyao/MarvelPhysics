SUBROUTINE vox_stvk_at_quadr_mtr(val, mtr, X, H_invDmH, detDmH, gw) 
IMPLICIT NONE 
REAL(KIND=8) val(1, 1) 
REAL(KIND=8) mtr(2, 1) 
REAL(KIND=8) X(3, 8) 
REAL(KIND=8) H_invDmH(8, 3) 
REAL(KIND=8) detDmH(1, 1) 
REAL(KIND=8) gw(1, 1) 
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
tt1 = (X(3,8)*H_invDmH(8,1)+X(3,7)*H_invDmH(7,1)+X(3,6)*H_invDmH(&
&6,1)+X(3,5)*H_invDmH(5,1)+X(3,4)*H_invDmH(4,1)+H_invDmH(3,1)*X(3,&
&3)+H_invDmH(2,1)*X(3,2)+H_invDmH(1,1)*X(3,1))**2+(X(2,8)*H_invDmH&
&(8,1)+X(2,7)*H_invDmH(7,1)+X(2,6)*H_invDmH(6,1)+X(2,5)*H_invDmH(5&
&,1)+X(2,4)*H_invDmH(4,1)+X(2,3)*H_invDmH(3,1)+H_invDmH(2,1)*X(2,2&
&)+H_invDmH(1,1)*X(2,1))**2+(X(1,8)*H_invDmH(8,1)+X(1,7)*H_invDmH(&
&7,1)+X(1,6)*H_invDmH(6,1)+X(1,5)*H_invDmH(5,1)+X(1,4)*H_invDmH(4,&
&1)+X(1,3)*H_invDmH(3,1)+X(1,2)*H_invDmH(2,1)+H_invDmH(1,1)*X(1,1)&
&)**2-1
tt2 = X(1,1)**2
tt3 = X(2,1)**2
tt4 = X(1,2)**2
tt5 = X(2,2)**2
tt6 = X(3,1)**2
tt7 = X(1,3)**2
tt8 = X(2,3)**2
tt9 = X(3,2)**2
tt10 = X(3,3)**2
tt11 = X(1,4)**2
tt12 = X(2,4)**2
tt13 = X(3,4)**2
tt14 = X(1,5)**2
tt15 = X(2,5)**2
tt16 = X(3,5)**2
tt17 = X(1,6)**2
tt18 = X(2,6)**2
tt19 = X(3,6)**2
tt20 = X(1,7)**2
tt21 = X(2,7)**2
tt22 = X(3,7)**2
tt23 = X(1,8)**2
tt24 = X(2,8)**2
tt25 = X(3,8)**2
tt26 = (X(3,8)*H_invDmH(8,2)+X(3,7)*H_invDmH(7,2)+X(3,6)*H_invDmH&
&(6,2)+X(3,5)*H_invDmH(5,2)+X(3,4)*H_invDmH(4,2)+H_invDmH(3,2)*X(3&
&,3)+H_invDmH(2,2)*X(3,2)+H_invDmH(1,2)*X(3,1))**2+(X(2,8)*H_invDm&
&H(8,2)+X(2,7)*H_invDmH(7,2)+X(2,6)*H_invDmH(6,2)+X(2,5)*H_invDmH(&
&5,2)+X(2,4)*H_invDmH(4,2)+X(2,3)*H_invDmH(3,2)+H_invDmH(2,2)*X(2,&
&2)+H_invDmH(1,2)*X(2,1))**2+(X(1,8)*H_invDmH(8,2)+X(1,7)*H_invDmH&
&(7,2)+X(1,6)*H_invDmH(6,2)+X(1,5)*H_invDmH(5,2)+X(1,4)*H_invDmH(4&
&,2)+X(1,3)*H_invDmH(3,2)+X(1,2)*H_invDmH(2,2)+X(1,1)*H_invDmH(1,2&
&))**2-1
tt27 = (X(3,8)*H_invDmH(8,3)+X(3,7)*H_invDmH(7,3)+X(3,6)*H_invDmH&
&(6,3)+X(3,5)*H_invDmH(5,3)+X(3,4)*H_invDmH(4,3)+H_invDmH(3,3)*X(3&
&,3)+H_invDmH(2,3)*X(3,2)+H_invDmH(1,3)*X(3,1))**2+(X(2,8)*H_invDm&
&H(8,3)+X(2,7)*H_invDmH(7,3)+X(2,6)*H_invDmH(6,3)+X(2,5)*H_invDmH(&
&5,3)+X(2,4)*H_invDmH(4,3)+X(2,3)*H_invDmH(3,3)+X(2,2)*H_invDmH(2,&
&3)+H_invDmH(1,3)*X(2,1))**2+(X(1,8)*H_invDmH(8,3)+X(1,7)*H_invDmH&
&(7,3)+X(1,6)*H_invDmH(6,3)+X(1,5)*H_invDmH(5,3)+X(1,4)*H_invDmH(4&
&,3)+X(1,3)*H_invDmH(3,3)+X(1,2)*H_invDmH(2,3)+X(1,1)*H_invDmH(1,3&
&))**2-1
val(1,1) = detDmH(1,1)*gw(1,1)*((mtr(2,1)*(tt27/2.0E+0+tt26/2.0E+&
&0+tt1/2.0E+0)**2)/2.0E+0+mtr(1,1)*(tt27**2/4.0E+0+(tt25*H_invDmH(&
&8,2)*H_invDmH(8,3)+tt24*H_invDmH(8,2)*H_invDmH(8,3)+tt23*H_invDmH&
&(8,2)*H_invDmH(8,3)+X(3,7)*X(3,8)*H_invDmH(7,2)*H_invDmH(8,3)+X(2&
&,7)*X(2,8)*H_invDmH(7,2)*H_invDmH(8,3)+X(1,7)*X(1,8)*H_invDmH(7,2&
&)*H_invDmH(8,3)+X(3,6)*X(3,8)*H_invDmH(6,2)*H_invDmH(8,3)+X(2,6)*&
&X(2,8)*H_invDmH(6,2)*H_invDmH(8,3)+X(1,6)*X(1,8)*H_invDmH(6,2)*H_&
&invDmH(8,3)+X(3,5)*X(3,8)*H_invDmH(5,2)*H_invDmH(8,3)+X(2,5)*X(2,&
&8)*H_invDmH(5,2)*H_invDmH(8,3)+X(1,5)*X(1,8)*H_invDmH(5,2)*H_invD&
&mH(8,3)+X(3,4)*X(3,8)*H_invDmH(4,2)*H_invDmH(8,3)+X(2,4)*X(2,8)*H&
&_invDmH(4,2)*H_invDmH(8,3)+X(1,4)*X(1,8)*H_invDmH(4,2)*H_invDmH(8&
&,3)+H_invDmH(3,2)*X(3,3)*X(3,8)*H_invDmH(8,3)+H_invDmH(2,2)*X(3,2&
&)*X(3,8)*H_invDmH(8,3)+H_invDmH(1,2)*X(3,1)*X(3,8)*H_invDmH(8,3)+&
&X(2,3)*X(2,8)*H_invDmH(3,2)*H_invDmH(8,3)+X(1,3)*X(1,8)*H_invDmH(&
&3,2)*H_invDmH(8,3)+H_invDmH(2,2)*X(2,2)*X(2,8)*H_invDmH(8,3)+H_in&
&vDmH(1,2)*X(2,1)*X(2,8)*H_invDmH(8,3)+X(1,2)*X(1,8)*H_invDmH(2,2)&
&*H_invDmH(8,3)+X(1,1)*H_invDmH(1,2)*X(1,8)*H_invDmH(8,3)+X(3,7)*X&
&(3,8)*H_invDmH(7,3)*H_invDmH(8,2)+X(2,7)*X(2,8)*H_invDmH(7,3)*H_i&
&nvDmH(8,2)+X(1,7)*X(1,8)*H_invDmH(7,3)*H_invDmH(8,2)+X(3,6)*X(3,8&
&)*H_invDmH(6,3)*H_invDmH(8,2)+X(2,6)*X(2,8)*H_invDmH(6,3)*H_invDm&
&H(8,2)+X(1,6)*X(1,8)*H_invDmH(6,3)*H_invDmH(8,2)+X(3,5)*X(3,8)*H_&
&invDmH(5,3)*H_invDmH(8,2)+X(2,5)*X(2,8)*H_invDmH(5,3)*H_invDmH(8,&
&2)+X(1,5)*X(1,8)*H_invDmH(5,3)*H_invDmH(8,2)+X(3,4)*X(3,8)*H_invD&
&mH(4,3)*H_invDmH(8,2)+X(2,4)*X(2,8)*H_invDmH(4,3)*H_invDmH(8,2)+X&
&(1,4)*X(1,8)*H_invDmH(4,3)*H_invDmH(8,2)+H_invDmH(3,3)*X(3,3)*X(3&
&,8)*H_invDmH(8,2)+H_invDmH(2,3)*X(3,2)*X(3,8)*H_invDmH(8,2)+H_inv&
&DmH(1,3)*X(3,1)*X(3,8)*H_invDmH(8,2)+X(2,3)*X(2,8)*H_invDmH(3,3)*&
&H_invDmH(8,2)+X(1,3)*X(1,8)*H_invDmH(3,3)*H_invDmH(8,2)+X(2,2)*H_&
&invDmH(2,3)*X(2,8)*H_invDmH(8,2)+H_invDmH(1,3)*X(2,1)*X(2,8)*H_in&
&vDmH(8,2)+X(1,2)*X(1,8)*H_invDmH(2,3)*H_invDmH(8,2)+X(1,1)*H_invD&
&mH(1,3)*X(1,8)*H_invDmH(8,2)+tt22*H_invDmH(7,2)*H_invDmH(7,3)+tt2&
&1*H_invDmH(7,2)*H_invDmH(7,3)+tt20*H_invDmH(7,2)*H_invDmH(7,3)+X(&
&3,6)*X(3,7)*H_invDmH(6,2)*H_invDmH(7,3)+X(2,6)*X(2,7)*H_invDmH(6,&
&2)*H_invDmH(7,3)+X(1,6)*X(1,7)*H_invDmH(6,2)*H_invDmH(7,3)+X(3,5)&
&*X(3,7)*H_invDmH(5,2)*H_invDmH(7,3)+X(2,5)*X(2,7)*H_invDmH(5,2)*H&
&_invDmH(7,3)+X(1,5)*X(1,7)*H_invDmH(5,2)*H_invDmH(7,3)+X(3,4)*X(3&
&,7)*H_invDmH(4,2)*H_invDmH(7,3)+X(2,4)*X(2,7)*H_invDmH(4,2)*H_inv&
&DmH(7,3)+X(1,4)*X(1,7)*H_invDmH(4,2)*H_invDmH(7,3)+H_invDmH(3,2)*&
&X(3,3)*X(3,7)*H_invDmH(7,3)+H_invDmH(2,2)*X(3,2)*X(3,7)*H_invDmH(&
&7,3)+H_invDmH(1,2)*X(3,1)*X(3,7)*H_invDmH(7,3)+X(2,3)*X(2,7)*H_in&
&vDmH(3,2)*H_invDmH(7,3)+X(1,3)*X(1,7)*H_invDmH(3,2)*H_invDmH(7,3)&
&+H_invDmH(2,2)*X(2,2)*X(2,7)*H_invDmH(7,3)+H_invDmH(1,2)*X(2,1)*X&
&(2,7)*H_invDmH(7,3)+X(1,2)*X(1,7)*H_invDmH(2,2)*H_invDmH(7,3)+X(1&
&,1)*H_invDmH(1,2)*X(1,7)*H_invDmH(7,3)+X(3,6)*X(3,7)*H_invDmH(6,3&
&)*H_invDmH(7,2)+X(2,6)*X(2,7)*H_invDmH(6,3)*H_invDmH(7,2)+X(1,6)*&
&X(1,7)*H_invDmH(6,3)*H_invDmH(7,2)+X(3,5)*X(3,7)*H_invDmH(5,3)*H_&
&invDmH(7,2)+X(2,5)*X(2,7)*H_invDmH(5,3)*H_invDmH(7,2)+X(1,5)*X(1,&
&7)*H_invDmH(5,3)*H_invDmH(7,2)+X(3,4)*X(3,7)*H_invDmH(4,3)*H_invD&
&mH(7,2)+X(2,4)*X(2,7)*H_invDmH(4,3)*H_invDmH(7,2)+X(1,4)*X(1,7)*H&
&_invDmH(4,3)*H_invDmH(7,2)+H_invDmH(3,3)*X(3,3)*X(3,7)*H_invDmH(7&
&,2)+H_invDmH(2,3)*X(3,2)*X(3,7)*H_invDmH(7,2)+H_invDmH(1,3)*X(3,1&
&)*X(3,7)*H_invDmH(7,2)+X(2,3)*X(2,7)*H_invDmH(3,3)*H_invDmH(7,2)+&
&X(1,3)*X(1,7)*H_invDmH(3,3)*H_invDmH(7,2)+X(2,2)*H_invDmH(2,3)*X(&
&2,7)*H_invDmH(7,2)+H_invDmH(1,3)*X(2,1)*X(2,7)*H_invDmH(7,2)+X(1,&
&2)*X(1,7)*H_invDmH(2,3)*H_invDmH(7,2)+X(1,1)*H_invDmH(1,3)*X(1,7)&
&*H_invDmH(7,2)+tt19*H_invDmH(6,2)*H_invDmH(6,3)+tt18*H_invDmH(6,2&
&)*H_invDmH(6,3)+tt17*H_invDmH(6,2)*H_invDmH(6,3)+X(3,5)*X(3,6)*H_&
&invDmH(5,2)*H_invDmH(6,3)+X(2,5)*X(2,6)*H_invDmH(5,2)*H_invDmH(6,&
&3)+X(1,5)*X(1,6)*H_invDmH(5,2)*H_invDmH(6,3)+X(3,4)*X(3,6)*H_invD&
&mH(4,2)*H_invDmH(6,3)+X(2,4)*X(2,6)*H_invDmH(4,2)*H_invDmH(6,3)+X&
&(1,4)*X(1,6)*H_invDmH(4,2)*H_invDmH(6,3)+H_invDmH(3,2)*X(3,3)*X(3&
&,6)*H_invDmH(6,3)+H_invDmH(2,2)*X(3,2)*X(3,6)*H_invDmH(6,3)+H_inv&
&DmH(1,2)*X(3,1)*X(3,6)*H_invDmH(6,3)+X(2,3)*X(2,6)*H_invDmH(3,2)*&
&H_invDmH(6,3)+X(1,3)*X(1,6)*H_invDmH(3,2)*H_invDmH(6,3)+H_invDmH(&
&2,2)*X(2,2)*X(2,6)*H_invDmH(6,3)+H_invDmH(1,2)*X(2,1)*X(2,6)*H_in&
&vDmH(6,3)+X(1,2)*X(1,6)*H_invDmH(2,2)*H_invDmH(6,3)+X(1,1)*H_invD&
&mH(1,2)*X(1,6)*H_invDmH(6,3)+X(3,5)*X(3,6)*H_invDmH(5,3)*H_invDmH&
&(6,2)+X(2,5)*X(2,6)*H_invDmH(5,3)*H_invDmH(6,2)+X(1,5)*X(1,6)*H_i&
&nvDmH(5,3)*H_invDmH(6,2)+X(3,4)*X(3,6)*H_invDmH(4,3)*H_invDmH(6,2&
&)+X(2,4)*X(2,6)*H_invDmH(4,3)*H_invDmH(6,2)+X(1,4)*X(1,6)*H_invDm&
&H(4,3)*H_invDmH(6,2)+H_invDmH(3,3)*X(3,3)*X(3,6)*H_invDmH(6,2)+H_&
&invDmH(2,3)*X(3,2)*X(3,6)*H_invDmH(6,2)+H_invDmH(1,3)*X(3,1)*X(3,&
&6)*H_invDmH(6,2)+X(2,3)*X(2,6)*H_invDmH(3,3)*H_invDmH(6,2)+X(1,3)&
&*X(1,6)*H_invDmH(3,3)*H_invDmH(6,2)+X(2,2)*H_invDmH(2,3)*X(2,6)*H&
&_invDmH(6,2)+H_invDmH(1,3)*X(2,1)*X(2,6)*H_invDmH(6,2)+X(1,2)*X(1&
&,6)*H_invDmH(2,3)*H_invDmH(6,2)+X(1,1)*H_invDmH(1,3)*X(1,6)*H_inv&
&DmH(6,2)+tt16*H_invDmH(5,2)*H_invDmH(5,3)+tt15*H_invDmH(5,2)*H_in&
&vDmH(5,3)+tt14*H_invDmH(5,2)*H_invDmH(5,3)+X(3,4)*X(3,5)*H_invDmH&
&(4,2)*H_invDmH(5,3)+X(2,4)*X(2,5)*H_invDmH(4,2)*H_invDmH(5,3)+X(1&
&,4)*X(1,5)*H_invDmH(4,2)*H_invDmH(5,3)+H_invDmH(3,2)*X(3,3)*X(3,5&
&)*H_invDmH(5,3)+H_invDmH(2,2)*X(3,2)*X(3,5)*H_invDmH(5,3)+H_invDm&
&H(1,2)*X(3,1)*X(3,5)*H_invDmH(5,3)+X(2,3)*X(2,5)*H_invDmH(3,2)*H_&
&invDmH(5,3)+X(1,3)*X(1,5)*H_invDmH(3,2)*H_invDmH(5,3)+H_invDmH(2,&
&2)*X(2,2)*X(2,5)*H_invDmH(5,3)+H_invDmH(1,2)*X(2,1)*X(2,5)*H_invD&
&mH(5,3)+X(1,2)*X(1,5)*H_invDmH(2,2)*H_invDmH(5,3)+X(1,1)*H_invDmH&
&(1,2)*X(1,5)*H_invDmH(5,3)+X(3,4)*X(3,5)*H_invDmH(4,3)*H_invDmH(5&
&,2)+X(2,4)*X(2,5)*H_invDmH(4,3)*H_invDmH(5,2)+X(1,4)*X(1,5)*H_inv&
&DmH(4,3)*H_invDmH(5,2)+H_invDmH(3,3)*X(3,3)*X(3,5)*H_invDmH(5,2)+&
&H_invDmH(2,3)*X(3,2)*X(3,5)*H_invDmH(5,2)+H_invDmH(1,3)*X(3,1)*X(&
&3,5)*H_invDmH(5,2)+X(2,3)*X(2,5)*H_invDmH(3,3)*H_invDmH(5,2)+X(1,&
&3)*X(1,5)*H_invDmH(3,3)*H_invDmH(5,2)+X(2,2)*H_invDmH(2,3)*X(2,5)&
&*H_invDmH(5,2)+H_invDmH(1,3)*X(2,1)*X(2,5)*H_invDmH(5,2)+X(1,2)*X&
&(1,5)*H_invDmH(2,3)*H_invDmH(5,2)+X(1,1)*H_invDmH(1,3)*X(1,5)*H_i&
&nvDmH(5,2)+tt13*H_invDmH(4,2)*H_invDmH(4,3)+tt12*H_invDmH(4,2)*H_&
&invDmH(4,3)+tt11*H_invDmH(4,2)*H_invDmH(4,3)+H_invDmH(3,2)*X(3,3)&
&*X(3,4)*H_invDmH(4,3)+H_invDmH(2,2)*X(3,2)*X(3,4)*H_invDmH(4,3)+H&
&_invDmH(1,2)*X(3,1)*X(3,4)*H_invDmH(4,3)+X(2,3)*X(2,4)*H_invDmH(3&
&,2)*H_invDmH(4,3)+X(1,3)*X(1,4)*H_invDmH(3,2)*H_invDmH(4,3)+H_inv&
&DmH(2,2)*X(2,2)*X(2,4)*H_invDmH(4,3)+H_invDmH(1,2)*X(2,1)*X(2,4)*&
&H_invDmH(4,3)+X(1,2)*X(1,4)*H_invDmH(2,2)*H_invDmH(4,3)+X(1,1)*H_&
&invDmH(1,2)*X(1,4)*H_invDmH(4,3)+H_invDmH(3,3)*X(3,3)*X(3,4)*H_in&
&vDmH(4,2)+H_invDmH(2,3)*X(3,2)*X(3,4)*H_invDmH(4,2)+H_invDmH(1,3)&
&*X(3,1)*X(3,4)*H_invDmH(4,2)+X(2,3)*X(2,4)*H_invDmH(3,3)*H_invDmH&
&(4,2)+X(1,3)*X(1,4)*H_invDmH(3,3)*H_invDmH(4,2)+X(2,2)*H_invDmH(2&
&,3)*X(2,4)*H_invDmH(4,2)+H_invDmH(1,3)*X(2,1)*X(2,4)*H_invDmH(4,2&
&)+X(1,2)*X(1,4)*H_invDmH(2,3)*H_invDmH(4,2)+X(1,1)*H_invDmH(1,3)*&
&X(1,4)*H_invDmH(4,2)+H_invDmH(3,2)*H_invDmH(3,3)*tt10+H_invDmH(2,&
&2)*X(3,2)*H_invDmH(3,3)*X(3,3)+H_invDmH(1,2)*X(3,1)*H_invDmH(3,3)&
&*X(3,3)+H_invDmH(2,3)*H_invDmH(3,2)*X(3,2)*X(3,3)+H_invDmH(1,3)*X&
&(3,1)*H_invDmH(3,2)*X(3,3)+tt8*H_invDmH(3,2)*H_invDmH(3,3)+tt7*H_&
&invDmH(3,2)*H_invDmH(3,3)+H_invDmH(2,2)*X(2,2)*X(2,3)*H_invDmH(3,&
&3)+H_invDmH(1,2)*X(2,1)*X(2,3)*H_invDmH(3,3)+X(1,2)*X(1,3)*H_invD&
&mH(2,2)*H_invDmH(3,3)+X(1,1)*H_invDmH(1,2)*X(1,3)*H_invDmH(3,3)+H&
&_invDmH(2,2)*H_invDmH(2,3)*tt9+H_invDmH(1,2)*H_invDmH(2,3)*X(3,1)&
&*X(3,2)+H_invDmH(1,3)*H_invDmH(2,2)*X(3,1)*X(3,2)+X(2,2)*H_invDmH&
&(2,3)*X(2,3)*H_invDmH(3,2)+H_invDmH(1,3)*X(2,1)*X(2,3)*H_invDmH(3&
&,2)+X(1,2)*X(1,3)*H_invDmH(2,3)*H_invDmH(3,2)+X(1,1)*H_invDmH(1,3&
&)*X(1,3)*H_invDmH(3,2)+H_invDmH(1,2)*H_invDmH(1,3)*tt6+H_invDmH(2&
&,2)*tt5*H_invDmH(2,3)+H_invDmH(1,2)*X(2,1)*X(2,2)*H_invDmH(2,3)+t&
&t4*H_invDmH(2,2)*H_invDmH(2,3)+X(1,1)*H_invDmH(1,2)*X(1,2)*H_invD&
&mH(2,3)+H_invDmH(1,3)*X(2,1)*H_invDmH(2,2)*X(2,2)+X(1,1)*X(1,2)*H&
&_invDmH(1,3)*H_invDmH(2,2)+H_invDmH(1,2)*H_invDmH(1,3)*tt3+tt2*H_&
&invDmH(1,2)*H_invDmH(1,3))**2/2.0E+0+(tt25*H_invDmH(8,1)*H_invDmH&
&(8,3)+tt24*H_invDmH(8,1)*H_invDmH(8,3)+tt23*H_invDmH(8,1)*H_invDm&
&H(8,3)+X(3,7)*X(3,8)*H_invDmH(7,1)*H_invDmH(8,3)+X(2,7)*X(2,8)*H_&
&invDmH(7,1)*H_invDmH(8,3)+X(1,7)*X(1,8)*H_invDmH(7,1)*H_invDmH(8,&
&3)+X(3,6)*X(3,8)*H_invDmH(6,1)*H_invDmH(8,3)+X(2,6)*X(2,8)*H_invD&
&mH(6,1)*H_invDmH(8,3)+X(1,6)*X(1,8)*H_invDmH(6,1)*H_invDmH(8,3)+X&
&(3,5)*X(3,8)*H_invDmH(5,1)*H_invDmH(8,3)+X(2,5)*X(2,8)*H_invDmH(5&
&,1)*H_invDmH(8,3)+X(1,5)*X(1,8)*H_invDmH(5,1)*H_invDmH(8,3)+X(3,4&
&)*X(3,8)*H_invDmH(4,1)*H_invDmH(8,3)+X(2,4)*X(2,8)*H_invDmH(4,1)*&
&H_invDmH(8,3)+X(1,4)*X(1,8)*H_invDmH(4,1)*H_invDmH(8,3)+H_invDmH(&
&3,1)*X(3,3)*X(3,8)*H_invDmH(8,3)+H_invDmH(2,1)*X(3,2)*X(3,8)*H_in&
&vDmH(8,3)+H_invDmH(1,1)*X(3,1)*X(3,8)*H_invDmH(8,3)+X(2,3)*X(2,8)&
&*H_invDmH(3,1)*H_invDmH(8,3)+X(1,3)*X(1,8)*H_invDmH(3,1)*H_invDmH&
&(8,3)+H_invDmH(2,1)*X(2,2)*X(2,8)*H_invDmH(8,3)+H_invDmH(1,1)*X(2&
&,1)*X(2,8)*H_invDmH(8,3)+X(1,2)*X(1,8)*H_invDmH(2,1)*H_invDmH(8,3&
&)+H_invDmH(1,1)*X(1,1)*X(1,8)*H_invDmH(8,3)+X(3,7)*X(3,8)*H_invDm&
&H(7,3)*H_invDmH(8,1)+X(2,7)*X(2,8)*H_invDmH(7,3)*H_invDmH(8,1)+X(&
&1,7)*X(1,8)*H_invDmH(7,3)*H_invDmH(8,1)+X(3,6)*X(3,8)*H_invDmH(6,&
&3)*H_invDmH(8,1)+X(2,6)*X(2,8)*H_invDmH(6,3)*H_invDmH(8,1)+X(1,6)&
&*X(1,8)*H_invDmH(6,3)*H_invDmH(8,1)+X(3,5)*X(3,8)*H_invDmH(5,3)*H&
&_invDmH(8,1)+X(2,5)*X(2,8)*H_invDmH(5,3)*H_invDmH(8,1)+X(1,5)*X(1&
&,8)*H_invDmH(5,3)*H_invDmH(8,1)+X(3,4)*X(3,8)*H_invDmH(4,3)*H_inv&
&DmH(8,1)+X(2,4)*X(2,8)*H_invDmH(4,3)*H_invDmH(8,1)+X(1,4)*X(1,8)*&
&H_invDmH(4,3)*H_invDmH(8,1)+H_invDmH(3,3)*X(3,3)*X(3,8)*H_invDmH(&
&8,1)+H_invDmH(2,3)*X(3,2)*X(3,8)*H_invDmH(8,1)+H_invDmH(1,3)*X(3,&
&1)*X(3,8)*H_invDmH(8,1)+X(2,3)*X(2,8)*H_invDmH(3,3)*H_invDmH(8,1)&
&+X(1,3)*X(1,8)*H_invDmH(3,3)*H_invDmH(8,1)+X(2,2)*H_invDmH(2,3)*X&
&(2,8)*H_invDmH(8,1)+H_invDmH(1,3)*X(2,1)*X(2,8)*H_invDmH(8,1)+X(1&
&,2)*X(1,8)*H_invDmH(2,3)*H_invDmH(8,1)+X(1,1)*H_invDmH(1,3)*X(1,8&
&)*H_invDmH(8,1)+tt22*H_invDmH(7,1)*H_invDmH(7,3)+tt21*H_invDmH(7,&
&1)*H_invDmH(7,3)+tt20*H_invDmH(7,1)*H_invDmH(7,3)+X(3,6)*X(3,7)*H&
&_invDmH(6,1)*H_invDmH(7,3)+X(2,6)*X(2,7)*H_invDmH(6,1)*H_invDmH(7&
&,3)+X(1,6)*X(1,7)*H_invDmH(6,1)*H_invDmH(7,3)+X(3,5)*X(3,7)*H_inv&
&DmH(5,1)*H_invDmH(7,3)+X(2,5)*X(2,7)*H_invDmH(5,1)*H_invDmH(7,3)+&
&X(1,5)*X(1,7)*H_invDmH(5,1)*H_invDmH(7,3)+X(3,4)*X(3,7)*H_invDmH(&
&4,1)*H_invDmH(7,3)+X(2,4)*X(2,7)*H_invDmH(4,1)*H_invDmH(7,3)+X(1,&
&4)*X(1,7)*H_invDmH(4,1)*H_invDmH(7,3)+H_invDmH(3,1)*X(3,3)*X(3,7)&
&*H_invDmH(7,3)+H_invDmH(2,1)*X(3,2)*X(3,7)*H_invDmH(7,3)+H_invDmH&
&(1,1)*X(3,1)*X(3,7)*H_invDmH(7,3)+X(2,3)*X(2,7)*H_invDmH(3,1)*H_i&
&nvDmH(7,3)+X(1,3)*X(1,7)*H_invDmH(3,1)*H_invDmH(7,3)+H_invDmH(2,1&
&)*X(2,2)*X(2,7)*H_invDmH(7,3)+H_invDmH(1,1)*X(2,1)*X(2,7)*H_invDm&
&H(7,3)+X(1,2)*X(1,7)*H_invDmH(2,1)*H_invDmH(7,3)+H_invDmH(1,1)*X(&
&1,1)*X(1,7)*H_invDmH(7,3)+X(3,6)*X(3,7)*H_invDmH(6,3)*H_invDmH(7,&
&1)+X(2,6)*X(2,7)*H_invDmH(6,3)*H_invDmH(7,1)+X(1,6)*X(1,7)*H_invD&
&mH(6,3)*H_invDmH(7,1)+X(3,5)*X(3,7)*H_invDmH(5,3)*H_invDmH(7,1)+X&
&(2,5)*X(2,7)*H_invDmH(5,3)*H_invDmH(7,1)+X(1,5)*X(1,7)*H_invDmH(5&
&,3)*H_invDmH(7,1)+X(3,4)*X(3,7)*H_invDmH(4,3)*H_invDmH(7,1)+X(2,4&
&)*X(2,7)*H_invDmH(4,3)*H_invDmH(7,1)+X(1,4)*X(1,7)*H_invDmH(4,3)*&
&H_invDmH(7,1)+H_invDmH(3,3)*X(3,3)*X(3,7)*H_invDmH(7,1)+H_invDmH(&
&2,3)*X(3,2)*X(3,7)*H_invDmH(7,1)+H_invDmH(1,3)*X(3,1)*X(3,7)*H_in&
&vDmH(7,1)+X(2,3)*X(2,7)*H_invDmH(3,3)*H_invDmH(7,1)+X(1,3)*X(1,7)&
&*H_invDmH(3,3)*H_invDmH(7,1)+X(2,2)*H_invDmH(2,3)*X(2,7)*H_invDmH&
&(7,1)+H_invDmH(1,3)*X(2,1)*X(2,7)*H_invDmH(7,1)+X(1,2)*X(1,7)*H_i&
&nvDmH(2,3)*H_invDmH(7,1)+X(1,1)*H_invDmH(1,3)*X(1,7)*H_invDmH(7,1&
&)+tt19*H_invDmH(6,1)*H_invDmH(6,3)+tt18*H_invDmH(6,1)*H_invDmH(6,&
&3)+tt17*H_invDmH(6,1)*H_invDmH(6,3)+X(3,5)*X(3,6)*H_invDmH(5,1)*H&
&_invDmH(6,3)+X(2,5)*X(2,6)*H_invDmH(5,1)*H_invDmH(6,3)+X(1,5)*X(1&
&,6)*H_invDmH(5,1)*H_invDmH(6,3)+X(3,4)*X(3,6)*H_invDmH(4,1)*H_inv&
&DmH(6,3)+X(2,4)*X(2,6)*H_invDmH(4,1)*H_invDmH(6,3)+X(1,4)*X(1,6)*&
&H_invDmH(4,1)*H_invDmH(6,3)+H_invDmH(3,1)*X(3,3)*X(3,6)*H_invDmH(&
&6,3)+H_invDmH(2,1)*X(3,2)*X(3,6)*H_invDmH(6,3)+H_invDmH(1,1)*X(3,&
&1)*X(3,6)*H_invDmH(6,3)+X(2,3)*X(2,6)*H_invDmH(3,1)*H_invDmH(6,3)&
&+X(1,3)*X(1,6)*H_invDmH(3,1)*H_invDmH(6,3)+H_invDmH(2,1)*X(2,2)*X&
&(2,6)*H_invDmH(6,3)+H_invDmH(1,1)*X(2,1)*X(2,6)*H_invDmH(6,3)+X(1&
&,2)*X(1,6)*H_invDmH(2,1)*H_invDmH(6,3)+H_invDmH(1,1)*X(1,1)*X(1,6&
&)*H_invDmH(6,3)+X(3,5)*X(3,6)*H_invDmH(5,3)*H_invDmH(6,1)+X(2,5)*&
&X(2,6)*H_invDmH(5,3)*H_invDmH(6,1)+X(1,5)*X(1,6)*H_invDmH(5,3)*H_&
&invDmH(6,1)+X(3,4)*X(3,6)*H_invDmH(4,3)*H_invDmH(6,1)+X(2,4)*X(2,&
&6)*H_invDmH(4,3)*H_invDmH(6,1)+X(1,4)*X(1,6)*H_invDmH(4,3)*H_invD&
&mH(6,1)+H_invDmH(3,3)*X(3,3)*X(3,6)*H_invDmH(6,1)+H_invDmH(2,3)*X&
&(3,2)*X(3,6)*H_invDmH(6,1)+H_invDmH(1,3)*X(3,1)*X(3,6)*H_invDmH(6&
&,1)+X(2,3)*X(2,6)*H_invDmH(3,3)*H_invDmH(6,1)+X(1,3)*X(1,6)*H_inv&
&DmH(3,3)*H_invDmH(6,1)+X(2,2)*H_invDmH(2,3)*X(2,6)*H_invDmH(6,1)+&
&H_invDmH(1,3)*X(2,1)*X(2,6)*H_invDmH(6,1)+X(1,2)*X(1,6)*H_invDmH(&
&2,3)*H_invDmH(6,1)+X(1,1)*H_invDmH(1,3)*X(1,6)*H_invDmH(6,1)+tt16&
&*H_invDmH(5,1)*H_invDmH(5,3)+tt15*H_invDmH(5,1)*H_invDmH(5,3)+tt1&
&4*H_invDmH(5,1)*H_invDmH(5,3)+X(3,4)*X(3,5)*H_invDmH(4,1)*H_invDm&
&H(5,3)+X(2,4)*X(2,5)*H_invDmH(4,1)*H_invDmH(5,3)+X(1,4)*X(1,5)*H_&
&invDmH(4,1)*H_invDmH(5,3)+H_invDmH(3,1)*X(3,3)*X(3,5)*H_invDmH(5,&
&3)+H_invDmH(2,1)*X(3,2)*X(3,5)*H_invDmH(5,3)+H_invDmH(1,1)*X(3,1)&
&*X(3,5)*H_invDmH(5,3)+X(2,3)*X(2,5)*H_invDmH(3,1)*H_invDmH(5,3)+X&
&(1,3)*X(1,5)*H_invDmH(3,1)*H_invDmH(5,3)+H_invDmH(2,1)*X(2,2)*X(2&
&,5)*H_invDmH(5,3)+H_invDmH(1,1)*X(2,1)*X(2,5)*H_invDmH(5,3)+X(1,2&
&)*X(1,5)*H_invDmH(2,1)*H_invDmH(5,3)+H_invDmH(1,1)*X(1,1)*X(1,5)*&
&H_invDmH(5,3)+X(3,4)*X(3,5)*H_invDmH(4,3)*H_invDmH(5,1)+X(2,4)*X(&
&2,5)*H_invDmH(4,3)*H_invDmH(5,1)+X(1,4)*X(1,5)*H_invDmH(4,3)*H_in&
&vDmH(5,1)+H_invDmH(3,3)*X(3,3)*X(3,5)*H_invDmH(5,1)+H_invDmH(2,3)&
&*X(3,2)*X(3,5)*H_invDmH(5,1)+H_invDmH(1,3)*X(3,1)*X(3,5)*H_invDmH&
&(5,1)+X(2,3)*X(2,5)*H_invDmH(3,3)*H_invDmH(5,1)+X(1,3)*X(1,5)*H_i&
&nvDmH(3,3)*H_invDmH(5,1)+X(2,2)*H_invDmH(2,3)*X(2,5)*H_invDmH(5,1&
&)+H_invDmH(1,3)*X(2,1)*X(2,5)*H_invDmH(5,1)+X(1,2)*X(1,5)*H_invDm&
&H(2,3)*H_invDmH(5,1)+X(1,1)*H_invDmH(1,3)*X(1,5)*H_invDmH(5,1)+tt&
&13*H_invDmH(4,1)*H_invDmH(4,3)+tt12*H_invDmH(4,1)*H_invDmH(4,3)+t&
&t11*H_invDmH(4,1)*H_invDmH(4,3)+H_invDmH(3,1)*X(3,3)*X(3,4)*H_inv&
&DmH(4,3)+H_invDmH(2,1)*X(3,2)*X(3,4)*H_invDmH(4,3)+H_invDmH(1,1)*&
&X(3,1)*X(3,4)*H_invDmH(4,3)+X(2,3)*X(2,4)*H_invDmH(3,1)*H_invDmH(&
&4,3)+X(1,3)*X(1,4)*H_invDmH(3,1)*H_invDmH(4,3)+H_invDmH(2,1)*X(2,&
&2)*X(2,4)*H_invDmH(4,3)+H_invDmH(1,1)*X(2,1)*X(2,4)*H_invDmH(4,3)&
&+X(1,2)*X(1,4)*H_invDmH(2,1)*H_invDmH(4,3)+H_invDmH(1,1)*X(1,1)*X&
&(1,4)*H_invDmH(4,3)+H_invDmH(3,3)*X(3,3)*X(3,4)*H_invDmH(4,1)+H_i&
&nvDmH(2,3)*X(3,2)*X(3,4)*H_invDmH(4,1)+H_invDmH(1,3)*X(3,1)*X(3,4&
&)*H_invDmH(4,1)+X(2,3)*X(2,4)*H_invDmH(3,3)*H_invDmH(4,1)+X(1,3)*&
&X(1,4)*H_invDmH(3,3)*H_invDmH(4,1)+X(2,2)*H_invDmH(2,3)*X(2,4)*H_&
&invDmH(4,1)+H_invDmH(1,3)*X(2,1)*X(2,4)*H_invDmH(4,1)+X(1,2)*X(1,&
&4)*H_invDmH(2,3)*H_invDmH(4,1)+X(1,1)*H_invDmH(1,3)*X(1,4)*H_invD&
&mH(4,1)+H_invDmH(3,1)*H_invDmH(3,3)*tt10+H_invDmH(2,1)*X(3,2)*H_i&
&nvDmH(3,3)*X(3,3)+H_invDmH(1,1)*X(3,1)*H_invDmH(3,3)*X(3,3)+H_inv&
&DmH(2,3)*H_invDmH(3,1)*X(3,2)*X(3,3)+H_invDmH(1,3)*H_invDmH(3,1)*&
&X(3,1)*X(3,3)+tt8*H_invDmH(3,1)*H_invDmH(3,3)+tt7*H_invDmH(3,1)*H&
&_invDmH(3,3)+H_invDmH(2,1)*X(2,2)*X(2,3)*H_invDmH(3,3)+H_invDmH(1&
&,1)*X(2,1)*X(2,3)*H_invDmH(3,3)+X(1,2)*X(1,3)*H_invDmH(2,1)*H_inv&
&DmH(3,3)+H_invDmH(1,1)*X(1,1)*X(1,3)*H_invDmH(3,3)+H_invDmH(2,1)*&
&H_invDmH(2,3)*tt9+H_invDmH(1,1)*H_invDmH(2,3)*X(3,1)*X(3,2)+H_inv&
&DmH(1,3)*H_invDmH(2,1)*X(3,1)*X(3,2)+H_invDmH(1,1)*H_invDmH(1,3)*&
&tt6+X(2,2)*H_invDmH(2,3)*X(2,3)*H_invDmH(3,1)+H_invDmH(1,3)*X(2,1&
&)*X(2,3)*H_invDmH(3,1)+X(1,2)*X(1,3)*H_invDmH(2,3)*H_invDmH(3,1)+&
&X(1,1)*H_invDmH(1,3)*X(1,3)*H_invDmH(3,1)+H_invDmH(2,1)*tt5*H_inv&
&DmH(2,3)+H_invDmH(1,1)*X(2,1)*X(2,2)*H_invDmH(2,3)+tt4*H_invDmH(2&
&,1)*H_invDmH(2,3)+H_invDmH(1,1)*X(1,1)*X(1,2)*H_invDmH(2,3)+H_inv&
&DmH(1,3)*H_invDmH(2,1)*X(2,1)*X(2,2)+H_invDmH(1,1)*H_invDmH(1,3)*&
&tt3+X(1,1)*X(1,2)*H_invDmH(1,3)*H_invDmH(2,1)+H_invDmH(1,1)*tt2*H&
&_invDmH(1,3))**2/2.0E+0+tt26**2/4.0E+0+(tt25*H_invDmH(8,1)*H_invD&
&mH(8,2)+tt24*H_invDmH(8,1)*H_invDmH(8,2)+tt23*H_invDmH(8,1)*H_inv&
&DmH(8,2)+X(3,7)*X(3,8)*H_invDmH(7,1)*H_invDmH(8,2)+X(2,7)*X(2,8)*&
&H_invDmH(7,1)*H_invDmH(8,2)+X(1,7)*X(1,8)*H_invDmH(7,1)*H_invDmH(&
&8,2)+X(3,6)*X(3,8)*H_invDmH(6,1)*H_invDmH(8,2)+X(2,6)*X(2,8)*H_in&
&vDmH(6,1)*H_invDmH(8,2)+X(1,6)*X(1,8)*H_invDmH(6,1)*H_invDmH(8,2)&
&+X(3,5)*X(3,8)*H_invDmH(5,1)*H_invDmH(8,2)+X(2,5)*X(2,8)*H_invDmH&
&(5,1)*H_invDmH(8,2)+X(1,5)*X(1,8)*H_invDmH(5,1)*H_invDmH(8,2)+X(3&
&,4)*X(3,8)*H_invDmH(4,1)*H_invDmH(8,2)+X(2,4)*X(2,8)*H_invDmH(4,1&
&)*H_invDmH(8,2)+X(1,4)*X(1,8)*H_invDmH(4,1)*H_invDmH(8,2)+H_invDm&
&H(3,1)*X(3,3)*X(3,8)*H_invDmH(8,2)+H_invDmH(2,1)*X(3,2)*X(3,8)*H_&
&invDmH(8,2)+H_invDmH(1,1)*X(3,1)*X(3,8)*H_invDmH(8,2)+X(2,3)*X(2,&
&8)*H_invDmH(3,1)*H_invDmH(8,2)+X(1,3)*X(1,8)*H_invDmH(3,1)*H_invD&
&mH(8,2)+H_invDmH(2,1)*X(2,2)*X(2,8)*H_invDmH(8,2)+H_invDmH(1,1)*X&
&(2,1)*X(2,8)*H_invDmH(8,2)+X(1,2)*X(1,8)*H_invDmH(2,1)*H_invDmH(8&
&,2)+H_invDmH(1,1)*X(1,1)*X(1,8)*H_invDmH(8,2)+X(3,7)*X(3,8)*H_inv&
&DmH(7,2)*H_invDmH(8,1)+X(2,7)*X(2,8)*H_invDmH(7,2)*H_invDmH(8,1)+&
&X(1,7)*X(1,8)*H_invDmH(7,2)*H_invDmH(8,1)+X(3,6)*X(3,8)*H_invDmH(&
&6,2)*H_invDmH(8,1)+X(2,6)*X(2,8)*H_invDmH(6,2)*H_invDmH(8,1)+X(1,&
&6)*X(1,8)*H_invDmH(6,2)*H_invDmH(8,1)+X(3,5)*X(3,8)*H_invDmH(5,2)&
&*H_invDmH(8,1)+X(2,5)*X(2,8)*H_invDmH(5,2)*H_invDmH(8,1)+X(1,5)*X&
&(1,8)*H_invDmH(5,2)*H_invDmH(8,1)+X(3,4)*X(3,8)*H_invDmH(4,2)*H_i&
&nvDmH(8,1)+X(2,4)*X(2,8)*H_invDmH(4,2)*H_invDmH(8,1)+X(1,4)*X(1,8&
&)*H_invDmH(4,2)*H_invDmH(8,1)+H_invDmH(3,2)*X(3,3)*X(3,8)*H_invDm&
&H(8,1)+H_invDmH(2,2)*X(3,2)*X(3,8)*H_invDmH(8,1)+H_invDmH(1,2)*X(&
&3,1)*X(3,8)*H_invDmH(8,1)+X(2,3)*X(2,8)*H_invDmH(3,2)*H_invDmH(8,&
&1)+X(1,3)*X(1,8)*H_invDmH(3,2)*H_invDmH(8,1)+H_invDmH(2,2)*X(2,2)&
&*X(2,8)*H_invDmH(8,1)+H_invDmH(1,2)*X(2,1)*X(2,8)*H_invDmH(8,1)+X&
&(1,2)*X(1,8)*H_invDmH(2,2)*H_invDmH(8,1)+X(1,1)*H_invDmH(1,2)*X(1&
&,8)*H_invDmH(8,1)+tt22*H_invDmH(7,1)*H_invDmH(7,2)+tt21*H_invDmH(&
&7,1)*H_invDmH(7,2)+tt20*H_invDmH(7,1)*H_invDmH(7,2)+X(3,6)*X(3,7)&
&*H_invDmH(6,1)*H_invDmH(7,2)+X(2,6)*X(2,7)*H_invDmH(6,1)*H_invDmH&
&(7,2)+X(1,6)*X(1,7)*H_invDmH(6,1)*H_invDmH(7,2)+X(3,5)*X(3,7)*H_i&
&nvDmH(5,1)*H_invDmH(7,2)+X(2,5)*X(2,7)*H_invDmH(5,1)*H_invDmH(7,2&
&)+X(1,5)*X(1,7)*H_invDmH(5,1)*H_invDmH(7,2)+X(3,4)*X(3,7)*H_invDm&
&H(4,1)*H_invDmH(7,2)+X(2,4)*X(2,7)*H_invDmH(4,1)*H_invDmH(7,2)+X(&
&1,4)*X(1,7)*H_invDmH(4,1)*H_invDmH(7,2)+H_invDmH(3,1)*X(3,3)*X(3,&
&7)*H_invDmH(7,2)+H_invDmH(2,1)*X(3,2)*X(3,7)*H_invDmH(7,2)+H_invD&
&mH(1,1)*X(3,1)*X(3,7)*H_invDmH(7,2)+X(2,3)*X(2,7)*H_invDmH(3,1)*H&
&_invDmH(7,2)+X(1,3)*X(1,7)*H_invDmH(3,1)*H_invDmH(7,2)+H_invDmH(2&
&,1)*X(2,2)*X(2,7)*H_invDmH(7,2)+H_invDmH(1,1)*X(2,1)*X(2,7)*H_inv&
&DmH(7,2)+X(1,2)*X(1,7)*H_invDmH(2,1)*H_invDmH(7,2)+H_invDmH(1,1)*&
&X(1,1)*X(1,7)*H_invDmH(7,2)+X(3,6)*X(3,7)*H_invDmH(6,2)*H_invDmH(&
&7,1)+X(2,6)*X(2,7)*H_invDmH(6,2)*H_invDmH(7,1)+X(1,6)*X(1,7)*H_in&
&vDmH(6,2)*H_invDmH(7,1)+X(3,5)*X(3,7)*H_invDmH(5,2)*H_invDmH(7,1)&
&+X(2,5)*X(2,7)*H_invDmH(5,2)*H_invDmH(7,1)+X(1,5)*X(1,7)*H_invDmH&
&(5,2)*H_invDmH(7,1)+X(3,4)*X(3,7)*H_invDmH(4,2)*H_invDmH(7,1)+X(2&
&,4)*X(2,7)*H_invDmH(4,2)*H_invDmH(7,1)+X(1,4)*X(1,7)*H_invDmH(4,2&
&)*H_invDmH(7,1)+H_invDmH(3,2)*X(3,3)*X(3,7)*H_invDmH(7,1)+H_invDm&
&H(2,2)*X(3,2)*X(3,7)*H_invDmH(7,1)+H_invDmH(1,2)*X(3,1)*X(3,7)*H_&
&invDmH(7,1)+X(2,3)*X(2,7)*H_invDmH(3,2)*H_invDmH(7,1)+X(1,3)*X(1,&
&7)*H_invDmH(3,2)*H_invDmH(7,1)+H_invDmH(2,2)*X(2,2)*X(2,7)*H_invD&
&mH(7,1)+H_invDmH(1,2)*X(2,1)*X(2,7)*H_invDmH(7,1)+X(1,2)*X(1,7)*H&
&_invDmH(2,2)*H_invDmH(7,1)+X(1,1)*H_invDmH(1,2)*X(1,7)*H_invDmH(7&
&,1)+tt19*H_invDmH(6,1)*H_invDmH(6,2)+tt18*H_invDmH(6,1)*H_invDmH(&
&6,2)+tt17*H_invDmH(6,1)*H_invDmH(6,2)+X(3,5)*X(3,6)*H_invDmH(5,1)&
&*H_invDmH(6,2)+X(2,5)*X(2,6)*H_invDmH(5,1)*H_invDmH(6,2)+X(1,5)*X&
&(1,6)*H_invDmH(5,1)*H_invDmH(6,2)+X(3,4)*X(3,6)*H_invDmH(4,1)*H_i&
&nvDmH(6,2)+X(2,4)*X(2,6)*H_invDmH(4,1)*H_invDmH(6,2)+X(1,4)*X(1,6&
&)*H_invDmH(4,1)*H_invDmH(6,2)+H_invDmH(3,1)*X(3,3)*X(3,6)*H_invDm&
&H(6,2)+H_invDmH(2,1)*X(3,2)*X(3,6)*H_invDmH(6,2)+H_invDmH(1,1)*X(&
&3,1)*X(3,6)*H_invDmH(6,2)+X(2,3)*X(2,6)*H_invDmH(3,1)*H_invDmH(6,&
&2)+X(1,3)*X(1,6)*H_invDmH(3,1)*H_invDmH(6,2)+H_invDmH(2,1)*X(2,2)&
&*X(2,6)*H_invDmH(6,2)+H_invDmH(1,1)*X(2,1)*X(2,6)*H_invDmH(6,2)+X&
&(1,2)*X(1,6)*H_invDmH(2,1)*H_invDmH(6,2)+H_invDmH(1,1)*X(1,1)*X(1&
&,6)*H_invDmH(6,2)+X(3,5)*X(3,6)*H_invDmH(5,2)*H_invDmH(6,1)+X(2,5&
&)*X(2,6)*H_invDmH(5,2)*H_invDmH(6,1)+X(1,5)*X(1,6)*H_invDmH(5,2)*&
&H_invDmH(6,1)+X(3,4)*X(3,6)*H_invDmH(4,2)*H_invDmH(6,1)+X(2,4)*X(&
&2,6)*H_invDmH(4,2)*H_invDmH(6,1)+X(1,4)*X(1,6)*H_invDmH(4,2)*H_in&
&vDmH(6,1)+H_invDmH(3,2)*X(3,3)*X(3,6)*H_invDmH(6,1)+H_invDmH(2,2)&
&*X(3,2)*X(3,6)*H_invDmH(6,1)+H_invDmH(1,2)*X(3,1)*X(3,6)*H_invDmH&
&(6,1)+X(2,3)*X(2,6)*H_invDmH(3,2)*H_invDmH(6,1)+X(1,3)*X(1,6)*H_i&
&nvDmH(3,2)*H_invDmH(6,1)+H_invDmH(2,2)*X(2,2)*X(2,6)*H_invDmH(6,1&
&)+H_invDmH(1,2)*X(2,1)*X(2,6)*H_invDmH(6,1)+X(1,2)*X(1,6)*H_invDm&
&H(2,2)*H_invDmH(6,1)+X(1,1)*H_invDmH(1,2)*X(1,6)*H_invDmH(6,1)+tt&
&16*H_invDmH(5,1)*H_invDmH(5,2)+tt15*H_invDmH(5,1)*H_invDmH(5,2)+t&
&t14*H_invDmH(5,1)*H_invDmH(5,2)+X(3,4)*X(3,5)*H_invDmH(4,1)*H_inv&
&DmH(5,2)+X(2,4)*X(2,5)*H_invDmH(4,1)*H_invDmH(5,2)+X(1,4)*X(1,5)*&
&H_invDmH(4,1)*H_invDmH(5,2)+H_invDmH(3,1)*X(3,3)*X(3,5)*H_invDmH(&
&5,2)+H_invDmH(2,1)*X(3,2)*X(3,5)*H_invDmH(5,2)+H_invDmH(1,1)*X(3,&
&1)*X(3,5)*H_invDmH(5,2)+X(2,3)*X(2,5)*H_invDmH(3,1)*H_invDmH(5,2)&
&+X(1,3)*X(1,5)*H_invDmH(3,1)*H_invDmH(5,2)+H_invDmH(2,1)*X(2,2)*X&
&(2,5)*H_invDmH(5,2)+H_invDmH(1,1)*X(2,1)*X(2,5)*H_invDmH(5,2)+X(1&
&,2)*X(1,5)*H_invDmH(2,1)*H_invDmH(5,2)+H_invDmH(1,1)*X(1,1)*X(1,5&
&)*H_invDmH(5,2)+X(3,4)*X(3,5)*H_invDmH(4,2)*H_invDmH(5,1)+X(2,4)*&
&X(2,5)*H_invDmH(4,2)*H_invDmH(5,1)+X(1,4)*X(1,5)*H_invDmH(4,2)*H_&
&invDmH(5,1)+H_invDmH(3,2)*X(3,3)*X(3,5)*H_invDmH(5,1)+H_invDmH(2,&
&2)*X(3,2)*X(3,5)*H_invDmH(5,1)+H_invDmH(1,2)*X(3,1)*X(3,5)*H_invD&
&mH(5,1)+X(2,3)*X(2,5)*H_invDmH(3,2)*H_invDmH(5,1)+X(1,3)*X(1,5)*H&
&_invDmH(3,2)*H_invDmH(5,1)+H_invDmH(2,2)*X(2,2)*X(2,5)*H_invDmH(5&
&,1)+H_invDmH(1,2)*X(2,1)*X(2,5)*H_invDmH(5,1)+X(1,2)*X(1,5)*H_inv&
&DmH(2,2)*H_invDmH(5,1)+X(1,1)*H_invDmH(1,2)*X(1,5)*H_invDmH(5,1)+&
&tt13*H_invDmH(4,1)*H_invDmH(4,2)+tt12*H_invDmH(4,1)*H_invDmH(4,2)&
&+tt11*H_invDmH(4,1)*H_invDmH(4,2)+H_invDmH(3,1)*X(3,3)*X(3,4)*H_i&
&nvDmH(4,2)+H_invDmH(2,1)*X(3,2)*X(3,4)*H_invDmH(4,2)+H_invDmH(1,1&
&)*X(3,1)*X(3,4)*H_invDmH(4,2)+X(2,3)*X(2,4)*H_invDmH(3,1)*H_invDm&
&H(4,2)+X(1,3)*X(1,4)*H_invDmH(3,1)*H_invDmH(4,2)+H_invDmH(2,1)*X(&
&2,2)*X(2,4)*H_invDmH(4,2)+H_invDmH(1,1)*X(2,1)*X(2,4)*H_invDmH(4,&
&2)+X(1,2)*X(1,4)*H_invDmH(2,1)*H_invDmH(4,2)+H_invDmH(1,1)*X(1,1)&
&*X(1,4)*H_invDmH(4,2)+H_invDmH(3,2)*X(3,3)*X(3,4)*H_invDmH(4,1)+H&
&_invDmH(2,2)*X(3,2)*X(3,4)*H_invDmH(4,1)+H_invDmH(1,2)*X(3,1)*X(3&
&,4)*H_invDmH(4,1)+X(2,3)*X(2,4)*H_invDmH(3,2)*H_invDmH(4,1)+X(1,3&
&)*X(1,4)*H_invDmH(3,2)*H_invDmH(4,1)+H_invDmH(2,2)*X(2,2)*X(2,4)*&
&H_invDmH(4,1)+H_invDmH(1,2)*X(2,1)*X(2,4)*H_invDmH(4,1)+X(1,2)*X(&
&1,4)*H_invDmH(2,2)*H_invDmH(4,1)+X(1,1)*H_invDmH(1,2)*X(1,4)*H_in&
&vDmH(4,1)+H_invDmH(3,1)*H_invDmH(3,2)*tt10+H_invDmH(2,1)*H_invDmH&
&(3,2)*X(3,2)*X(3,3)+H_invDmH(2,2)*H_invDmH(3,1)*X(3,2)*X(3,3)+H_i&
&nvDmH(1,1)*X(3,1)*H_invDmH(3,2)*X(3,3)+H_invDmH(1,2)*H_invDmH(3,1&
&)*X(3,1)*X(3,3)+H_invDmH(2,1)*H_invDmH(2,2)*tt9+H_invDmH(1,1)*H_i&
&nvDmH(2,2)*X(3,1)*X(3,2)+H_invDmH(1,2)*H_invDmH(2,1)*X(3,1)*X(3,2&
&)+tt8*H_invDmH(3,1)*H_invDmH(3,2)+tt7*H_invDmH(3,1)*H_invDmH(3,2)&
&+H_invDmH(2,1)*X(2,2)*X(2,3)*H_invDmH(3,2)+H_invDmH(1,1)*X(2,1)*X&
&(2,3)*H_invDmH(3,2)+X(1,2)*X(1,3)*H_invDmH(2,1)*H_invDmH(3,2)+H_i&
&nvDmH(1,1)*X(1,1)*X(1,3)*H_invDmH(3,2)+H_invDmH(1,1)*H_invDmH(1,2&
&)*tt6+H_invDmH(2,2)*X(2,2)*X(2,3)*H_invDmH(3,1)+H_invDmH(1,2)*X(2&
&,1)*X(2,3)*H_invDmH(3,1)+X(1,2)*X(1,3)*H_invDmH(2,2)*H_invDmH(3,1&
&)+X(1,1)*H_invDmH(1,2)*X(1,3)*H_invDmH(3,1)+H_invDmH(2,1)*H_invDm&
&H(2,2)*tt5+H_invDmH(1,1)*X(2,1)*H_invDmH(2,2)*X(2,2)+H_invDmH(1,2&
&)*H_invDmH(2,1)*X(2,1)*X(2,2)+tt4*H_invDmH(2,1)*H_invDmH(2,2)+H_i&
&nvDmH(1,1)*X(1,1)*X(1,2)*H_invDmH(2,2)+H_invDmH(1,1)*H_invDmH(1,2&
&)*tt3+X(1,1)*H_invDmH(1,2)*X(1,2)*H_invDmH(2,1)+H_invDmH(1,1)*tt2&
&*H_invDmH(1,2))**2/2.0E+0+tt1**2/4.0E+0))
END 
SUBROUTINE vox_stvk_at_quadr_mtr_jac(jac, mtr, X, H_invDmH, detDmH, gw) 
IMPLICIT NONE 
REAL(KIND=8) jac(1, 2) 
REAL(KIND=8) mtr(2, 1) 
REAL(KIND=8) X(3, 8) 
REAL(KIND=8) H_invDmH(8, 3) 
REAL(KIND=8) detDmH(1, 1) 
REAL(KIND=8) gw(1, 1) 
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
tt1 = (X(3,8)*H_invDmH(8,1)+X(3,7)*H_invDmH(7,1)+X(3,6)*H_invDmH(&
&6,1)+X(3,5)*H_invDmH(5,1)+X(3,4)*H_invDmH(4,1)+H_invDmH(3,1)*X(3,&
&3)+H_invDmH(2,1)*X(3,2)+H_invDmH(1,1)*X(3,1))**2+(X(2,8)*H_invDmH&
&(8,1)+X(2,7)*H_invDmH(7,1)+X(2,6)*H_invDmH(6,1)+X(2,5)*H_invDmH(5&
&,1)+X(2,4)*H_invDmH(4,1)+X(2,3)*H_invDmH(3,1)+H_invDmH(2,1)*X(2,2&
&)+H_invDmH(1,1)*X(2,1))**2+(X(1,8)*H_invDmH(8,1)+X(1,7)*H_invDmH(&
&7,1)+X(1,6)*H_invDmH(6,1)+X(1,5)*H_invDmH(5,1)+X(1,4)*H_invDmH(4,&
&1)+X(1,3)*H_invDmH(3,1)+X(1,2)*H_invDmH(2,1)+H_invDmH(1,1)*X(1,1)&
&)**2-1
tt2 = X(1,1)**2
tt3 = X(2,1)**2
tt4 = X(1,2)**2
tt5 = X(2,2)**2
tt6 = X(3,1)**2
tt7 = X(1,3)**2
tt8 = X(2,3)**2
tt9 = X(3,2)**2
tt10 = X(3,3)**2
tt11 = X(1,4)**2
tt12 = X(2,4)**2
tt13 = X(3,4)**2
tt14 = X(1,5)**2
tt15 = X(2,5)**2
tt16 = X(3,5)**2
tt17 = X(1,6)**2
tt18 = X(2,6)**2
tt19 = X(3,6)**2
tt20 = X(1,7)**2
tt21 = X(2,7)**2
tt22 = X(3,7)**2
tt23 = X(1,8)**2
tt24 = X(2,8)**2
tt25 = X(3,8)**2
tt26 = (X(3,8)*H_invDmH(8,2)+X(3,7)*H_invDmH(7,2)+X(3,6)*H_invDmH&
&(6,2)+X(3,5)*H_invDmH(5,2)+X(3,4)*H_invDmH(4,2)+H_invDmH(3,2)*X(3&
&,3)+H_invDmH(2,2)*X(3,2)+H_invDmH(1,2)*X(3,1))**2+(X(2,8)*H_invDm&
&H(8,2)+X(2,7)*H_invDmH(7,2)+X(2,6)*H_invDmH(6,2)+X(2,5)*H_invDmH(&
&5,2)+X(2,4)*H_invDmH(4,2)+X(2,3)*H_invDmH(3,2)+H_invDmH(2,2)*X(2,&
&2)+H_invDmH(1,2)*X(2,1))**2+(X(1,8)*H_invDmH(8,2)+X(1,7)*H_invDmH&
&(7,2)+X(1,6)*H_invDmH(6,2)+X(1,5)*H_invDmH(5,2)+X(1,4)*H_invDmH(4&
&,2)+X(1,3)*H_invDmH(3,2)+X(1,2)*H_invDmH(2,2)+X(1,1)*H_invDmH(1,2&
&))**2-1
tt27 = (X(3,8)*H_invDmH(8,3)+X(3,7)*H_invDmH(7,3)+X(3,6)*H_invDmH&
&(6,3)+X(3,5)*H_invDmH(5,3)+X(3,4)*H_invDmH(4,3)+H_invDmH(3,3)*X(3&
&,3)+H_invDmH(2,3)*X(3,2)+H_invDmH(1,3)*X(3,1))**2+(X(2,8)*H_invDm&
&H(8,3)+X(2,7)*H_invDmH(7,3)+X(2,6)*H_invDmH(6,3)+X(2,5)*H_invDmH(&
&5,3)+X(2,4)*H_invDmH(4,3)+X(2,3)*H_invDmH(3,3)+X(2,2)*H_invDmH(2,&
&3)+H_invDmH(1,3)*X(2,1))**2+(X(1,8)*H_invDmH(8,3)+X(1,7)*H_invDmH&
&(7,3)+X(1,6)*H_invDmH(6,3)+X(1,5)*H_invDmH(5,3)+X(1,4)*H_invDmH(4&
&,3)+X(1,3)*H_invDmH(3,3)+X(1,2)*H_invDmH(2,3)+X(1,1)*H_invDmH(1,3&
&))**2-1
jac(1,1) = detDmH(1,1)*gw(1,1)*(tt27**2/4.0E+0+(tt25*H_invDmH(8,2&
&)*H_invDmH(8,3)+tt24*H_invDmH(8,2)*H_invDmH(8,3)+tt23*H_invDmH(8,&
&2)*H_invDmH(8,3)+X(3,7)*X(3,8)*H_invDmH(7,2)*H_invDmH(8,3)+X(2,7)&
&*X(2,8)*H_invDmH(7,2)*H_invDmH(8,3)+X(1,7)*X(1,8)*H_invDmH(7,2)*H&
&_invDmH(8,3)+X(3,6)*X(3,8)*H_invDmH(6,2)*H_invDmH(8,3)+X(2,6)*X(2&
&,8)*H_invDmH(6,2)*H_invDmH(8,3)+X(1,6)*X(1,8)*H_invDmH(6,2)*H_inv&
&DmH(8,3)+X(3,5)*X(3,8)*H_invDmH(5,2)*H_invDmH(8,3)+X(2,5)*X(2,8)*&
&H_invDmH(5,2)*H_invDmH(8,3)+X(1,5)*X(1,8)*H_invDmH(5,2)*H_invDmH(&
&8,3)+X(3,4)*X(3,8)*H_invDmH(4,2)*H_invDmH(8,3)+X(2,4)*X(2,8)*H_in&
&vDmH(4,2)*H_invDmH(8,3)+X(1,4)*X(1,8)*H_invDmH(4,2)*H_invDmH(8,3)&
&+H_invDmH(3,2)*X(3,3)*X(3,8)*H_invDmH(8,3)+H_invDmH(2,2)*X(3,2)*X&
&(3,8)*H_invDmH(8,3)+H_invDmH(1,2)*X(3,1)*X(3,8)*H_invDmH(8,3)+X(2&
&,3)*X(2,8)*H_invDmH(3,2)*H_invDmH(8,3)+X(1,3)*X(1,8)*H_invDmH(3,2&
&)*H_invDmH(8,3)+H_invDmH(2,2)*X(2,2)*X(2,8)*H_invDmH(8,3)+H_invDm&
&H(1,2)*X(2,1)*X(2,8)*H_invDmH(8,3)+X(1,2)*X(1,8)*H_invDmH(2,2)*H_&
&invDmH(8,3)+X(1,1)*H_invDmH(1,2)*X(1,8)*H_invDmH(8,3)+X(3,7)*X(3,&
&8)*H_invDmH(7,3)*H_invDmH(8,2)+X(2,7)*X(2,8)*H_invDmH(7,3)*H_invD&
&mH(8,2)+X(1,7)*X(1,8)*H_invDmH(7,3)*H_invDmH(8,2)+X(3,6)*X(3,8)*H&
&_invDmH(6,3)*H_invDmH(8,2)+X(2,6)*X(2,8)*H_invDmH(6,3)*H_invDmH(8&
&,2)+X(1,6)*X(1,8)*H_invDmH(6,3)*H_invDmH(8,2)+X(3,5)*X(3,8)*H_inv&
&DmH(5,3)*H_invDmH(8,2)+X(2,5)*X(2,8)*H_invDmH(5,3)*H_invDmH(8,2)+&
&X(1,5)*X(1,8)*H_invDmH(5,3)*H_invDmH(8,2)+X(3,4)*X(3,8)*H_invDmH(&
&4,3)*H_invDmH(8,2)+X(2,4)*X(2,8)*H_invDmH(4,3)*H_invDmH(8,2)+X(1,&
&4)*X(1,8)*H_invDmH(4,3)*H_invDmH(8,2)+H_invDmH(3,3)*X(3,3)*X(3,8)&
&*H_invDmH(8,2)+H_invDmH(2,3)*X(3,2)*X(3,8)*H_invDmH(8,2)+H_invDmH&
&(1,3)*X(3,1)*X(3,8)*H_invDmH(8,2)+X(2,3)*X(2,8)*H_invDmH(3,3)*H_i&
&nvDmH(8,2)+X(1,3)*X(1,8)*H_invDmH(3,3)*H_invDmH(8,2)+X(2,2)*H_inv&
&DmH(2,3)*X(2,8)*H_invDmH(8,2)+H_invDmH(1,3)*X(2,1)*X(2,8)*H_invDm&
&H(8,2)+X(1,2)*X(1,8)*H_invDmH(2,3)*H_invDmH(8,2)+X(1,1)*H_invDmH(&
&1,3)*X(1,8)*H_invDmH(8,2)+tt22*H_invDmH(7,2)*H_invDmH(7,3)+tt21*H&
&_invDmH(7,2)*H_invDmH(7,3)+tt20*H_invDmH(7,2)*H_invDmH(7,3)+X(3,6&
&)*X(3,7)*H_invDmH(6,2)*H_invDmH(7,3)+X(2,6)*X(2,7)*H_invDmH(6,2)*&
&H_invDmH(7,3)+X(1,6)*X(1,7)*H_invDmH(6,2)*H_invDmH(7,3)+X(3,5)*X(&
&3,7)*H_invDmH(5,2)*H_invDmH(7,3)+X(2,5)*X(2,7)*H_invDmH(5,2)*H_in&
&vDmH(7,3)+X(1,5)*X(1,7)*H_invDmH(5,2)*H_invDmH(7,3)+X(3,4)*X(3,7)&
&*H_invDmH(4,2)*H_invDmH(7,3)+X(2,4)*X(2,7)*H_invDmH(4,2)*H_invDmH&
&(7,3)+X(1,4)*X(1,7)*H_invDmH(4,2)*H_invDmH(7,3)+H_invDmH(3,2)*X(3&
&,3)*X(3,7)*H_invDmH(7,3)+H_invDmH(2,2)*X(3,2)*X(3,7)*H_invDmH(7,3&
&)+H_invDmH(1,2)*X(3,1)*X(3,7)*H_invDmH(7,3)+X(2,3)*X(2,7)*H_invDm&
&H(3,2)*H_invDmH(7,3)+X(1,3)*X(1,7)*H_invDmH(3,2)*H_invDmH(7,3)+H_&
&invDmH(2,2)*X(2,2)*X(2,7)*H_invDmH(7,3)+H_invDmH(1,2)*X(2,1)*X(2,&
&7)*H_invDmH(7,3)+X(1,2)*X(1,7)*H_invDmH(2,2)*H_invDmH(7,3)+X(1,1)&
&*H_invDmH(1,2)*X(1,7)*H_invDmH(7,3)+X(3,6)*X(3,7)*H_invDmH(6,3)*H&
&_invDmH(7,2)+X(2,6)*X(2,7)*H_invDmH(6,3)*H_invDmH(7,2)+X(1,6)*X(1&
&,7)*H_invDmH(6,3)*H_invDmH(7,2)+X(3,5)*X(3,7)*H_invDmH(5,3)*H_inv&
&DmH(7,2)+X(2,5)*X(2,7)*H_invDmH(5,3)*H_invDmH(7,2)+X(1,5)*X(1,7)*&
&H_invDmH(5,3)*H_invDmH(7,2)+X(3,4)*X(3,7)*H_invDmH(4,3)*H_invDmH(&
&7,2)+X(2,4)*X(2,7)*H_invDmH(4,3)*H_invDmH(7,2)+X(1,4)*X(1,7)*H_in&
&vDmH(4,3)*H_invDmH(7,2)+H_invDmH(3,3)*X(3,3)*X(3,7)*H_invDmH(7,2)&
&+H_invDmH(2,3)*X(3,2)*X(3,7)*H_invDmH(7,2)+H_invDmH(1,3)*X(3,1)*X&
&(3,7)*H_invDmH(7,2)+X(2,3)*X(2,7)*H_invDmH(3,3)*H_invDmH(7,2)+X(1&
&,3)*X(1,7)*H_invDmH(3,3)*H_invDmH(7,2)+X(2,2)*H_invDmH(2,3)*X(2,7&
&)*H_invDmH(7,2)+H_invDmH(1,3)*X(2,1)*X(2,7)*H_invDmH(7,2)+X(1,2)*&
&X(1,7)*H_invDmH(2,3)*H_invDmH(7,2)+X(1,1)*H_invDmH(1,3)*X(1,7)*H_&
&invDmH(7,2)+tt19*H_invDmH(6,2)*H_invDmH(6,3)+tt18*H_invDmH(6,2)*H&
&_invDmH(6,3)+tt17*H_invDmH(6,2)*H_invDmH(6,3)+X(3,5)*X(3,6)*H_inv&
&DmH(5,2)*H_invDmH(6,3)+X(2,5)*X(2,6)*H_invDmH(5,2)*H_invDmH(6,3)+&
&X(1,5)*X(1,6)*H_invDmH(5,2)*H_invDmH(6,3)+X(3,4)*X(3,6)*H_invDmH(&
&4,2)*H_invDmH(6,3)+X(2,4)*X(2,6)*H_invDmH(4,2)*H_invDmH(6,3)+X(1,&
&4)*X(1,6)*H_invDmH(4,2)*H_invDmH(6,3)+H_invDmH(3,2)*X(3,3)*X(3,6)&
&*H_invDmH(6,3)+H_invDmH(2,2)*X(3,2)*X(3,6)*H_invDmH(6,3)+H_invDmH&
&(1,2)*X(3,1)*X(3,6)*H_invDmH(6,3)+X(2,3)*X(2,6)*H_invDmH(3,2)*H_i&
&nvDmH(6,3)+X(1,3)*X(1,6)*H_invDmH(3,2)*H_invDmH(6,3)+H_invDmH(2,2&
&)*X(2,2)*X(2,6)*H_invDmH(6,3)+H_invDmH(1,2)*X(2,1)*X(2,6)*H_invDm&
&H(6,3)+X(1,2)*X(1,6)*H_invDmH(2,2)*H_invDmH(6,3)+X(1,1)*H_invDmH(&
&1,2)*X(1,6)*H_invDmH(6,3)+X(3,5)*X(3,6)*H_invDmH(5,3)*H_invDmH(6,&
&2)+X(2,5)*X(2,6)*H_invDmH(5,3)*H_invDmH(6,2)+X(1,5)*X(1,6)*H_invD&
&mH(5,3)*H_invDmH(6,2)+X(3,4)*X(3,6)*H_invDmH(4,3)*H_invDmH(6,2)+X&
&(2,4)*X(2,6)*H_invDmH(4,3)*H_invDmH(6,2)+X(1,4)*X(1,6)*H_invDmH(4&
&,3)*H_invDmH(6,2)+H_invDmH(3,3)*X(3,3)*X(3,6)*H_invDmH(6,2)+H_inv&
&DmH(2,3)*X(3,2)*X(3,6)*H_invDmH(6,2)+H_invDmH(1,3)*X(3,1)*X(3,6)*&
&H_invDmH(6,2)+X(2,3)*X(2,6)*H_invDmH(3,3)*H_invDmH(6,2)+X(1,3)*X(&
&1,6)*H_invDmH(3,3)*H_invDmH(6,2)+X(2,2)*H_invDmH(2,3)*X(2,6)*H_in&
&vDmH(6,2)+H_invDmH(1,3)*X(2,1)*X(2,6)*H_invDmH(6,2)+X(1,2)*X(1,6)&
&*H_invDmH(2,3)*H_invDmH(6,2)+X(1,1)*H_invDmH(1,3)*X(1,6)*H_invDmH&
&(6,2)+tt16*H_invDmH(5,2)*H_invDmH(5,3)+tt15*H_invDmH(5,2)*H_invDm&
&H(5,3)+tt14*H_invDmH(5,2)*H_invDmH(5,3)+X(3,4)*X(3,5)*H_invDmH(4,&
&2)*H_invDmH(5,3)+X(2,4)*X(2,5)*H_invDmH(4,2)*H_invDmH(5,3)+X(1,4)&
&*X(1,5)*H_invDmH(4,2)*H_invDmH(5,3)+H_invDmH(3,2)*X(3,3)*X(3,5)*H&
&_invDmH(5,3)+H_invDmH(2,2)*X(3,2)*X(3,5)*H_invDmH(5,3)+H_invDmH(1&
&,2)*X(3,1)*X(3,5)*H_invDmH(5,3)+X(2,3)*X(2,5)*H_invDmH(3,2)*H_inv&
&DmH(5,3)+X(1,3)*X(1,5)*H_invDmH(3,2)*H_invDmH(5,3)+H_invDmH(2,2)*&
&X(2,2)*X(2,5)*H_invDmH(5,3)+H_invDmH(1,2)*X(2,1)*X(2,5)*H_invDmH(&
&5,3)+X(1,2)*X(1,5)*H_invDmH(2,2)*H_invDmH(5,3)+X(1,1)*H_invDmH(1,&
&2)*X(1,5)*H_invDmH(5,3)+X(3,4)*X(3,5)*H_invDmH(4,3)*H_invDmH(5,2)&
&+X(2,4)*X(2,5)*H_invDmH(4,3)*H_invDmH(5,2)+X(1,4)*X(1,5)*H_invDmH&
&(4,3)*H_invDmH(5,2)+H_invDmH(3,3)*X(3,3)*X(3,5)*H_invDmH(5,2)+H_i&
&nvDmH(2,3)*X(3,2)*X(3,5)*H_invDmH(5,2)+H_invDmH(1,3)*X(3,1)*X(3,5&
&)*H_invDmH(5,2)+X(2,3)*X(2,5)*H_invDmH(3,3)*H_invDmH(5,2)+X(1,3)*&
&X(1,5)*H_invDmH(3,3)*H_invDmH(5,2)+X(2,2)*H_invDmH(2,3)*X(2,5)*H_&
&invDmH(5,2)+H_invDmH(1,3)*X(2,1)*X(2,5)*H_invDmH(5,2)+X(1,2)*X(1,&
&5)*H_invDmH(2,3)*H_invDmH(5,2)+X(1,1)*H_invDmH(1,3)*X(1,5)*H_invD&
&mH(5,2)+tt13*H_invDmH(4,2)*H_invDmH(4,3)+tt12*H_invDmH(4,2)*H_inv&
&DmH(4,3)+tt11*H_invDmH(4,2)*H_invDmH(4,3)+H_invDmH(3,2)*X(3,3)*X(&
&3,4)*H_invDmH(4,3)+H_invDmH(2,2)*X(3,2)*X(3,4)*H_invDmH(4,3)+H_in&
&vDmH(1,2)*X(3,1)*X(3,4)*H_invDmH(4,3)+X(2,3)*X(2,4)*H_invDmH(3,2)&
&*H_invDmH(4,3)+X(1,3)*X(1,4)*H_invDmH(3,2)*H_invDmH(4,3)+H_invDmH&
&(2,2)*X(2,2)*X(2,4)*H_invDmH(4,3)+H_invDmH(1,2)*X(2,1)*X(2,4)*H_i&
&nvDmH(4,3)+X(1,2)*X(1,4)*H_invDmH(2,2)*H_invDmH(4,3)+X(1,1)*H_inv&
&DmH(1,2)*X(1,4)*H_invDmH(4,3)+H_invDmH(3,3)*X(3,3)*X(3,4)*H_invDm&
&H(4,2)+H_invDmH(2,3)*X(3,2)*X(3,4)*H_invDmH(4,2)+H_invDmH(1,3)*X(&
&3,1)*X(3,4)*H_invDmH(4,2)+X(2,3)*X(2,4)*H_invDmH(3,3)*H_invDmH(4,&
&2)+X(1,3)*X(1,4)*H_invDmH(3,3)*H_invDmH(4,2)+X(2,2)*H_invDmH(2,3)&
&*X(2,4)*H_invDmH(4,2)+H_invDmH(1,3)*X(2,1)*X(2,4)*H_invDmH(4,2)+X&
&(1,2)*X(1,4)*H_invDmH(2,3)*H_invDmH(4,2)+X(1,1)*H_invDmH(1,3)*X(1&
&,4)*H_invDmH(4,2)+H_invDmH(3,2)*H_invDmH(3,3)*tt10+H_invDmH(2,2)*&
&X(3,2)*H_invDmH(3,3)*X(3,3)+H_invDmH(1,2)*X(3,1)*H_invDmH(3,3)*X(&
&3,3)+H_invDmH(2,3)*H_invDmH(3,2)*X(3,2)*X(3,3)+H_invDmH(1,3)*X(3,&
&1)*H_invDmH(3,2)*X(3,3)+tt8*H_invDmH(3,2)*H_invDmH(3,3)+tt7*H_inv&
&DmH(3,2)*H_invDmH(3,3)+H_invDmH(2,2)*X(2,2)*X(2,3)*H_invDmH(3,3)+&
&H_invDmH(1,2)*X(2,1)*X(2,3)*H_invDmH(3,3)+X(1,2)*X(1,3)*H_invDmH(&
&2,2)*H_invDmH(3,3)+X(1,1)*H_invDmH(1,2)*X(1,3)*H_invDmH(3,3)+H_in&
&vDmH(2,2)*H_invDmH(2,3)*tt9+H_invDmH(1,2)*H_invDmH(2,3)*X(3,1)*X(&
&3,2)+H_invDmH(1,3)*H_invDmH(2,2)*X(3,1)*X(3,2)+X(2,2)*H_invDmH(2,&
&3)*X(2,3)*H_invDmH(3,2)+H_invDmH(1,3)*X(2,1)*X(2,3)*H_invDmH(3,2)&
&+X(1,2)*X(1,3)*H_invDmH(2,3)*H_invDmH(3,2)+X(1,1)*H_invDmH(1,3)*X&
&(1,3)*H_invDmH(3,2)+H_invDmH(1,2)*H_invDmH(1,3)*tt6+H_invDmH(2,2)&
&*tt5*H_invDmH(2,3)+H_invDmH(1,2)*X(2,1)*X(2,2)*H_invDmH(2,3)+tt4*&
&H_invDmH(2,2)*H_invDmH(2,3)+X(1,1)*H_invDmH(1,2)*X(1,2)*H_invDmH(&
&2,3)+H_invDmH(1,3)*X(2,1)*H_invDmH(2,2)*X(2,2)+X(1,1)*X(1,2)*H_in&
&vDmH(1,3)*H_invDmH(2,2)+H_invDmH(1,2)*H_invDmH(1,3)*tt3+tt2*H_inv&
&DmH(1,2)*H_invDmH(1,3))**2/2.0E+0+(tt25*H_invDmH(8,1)*H_invDmH(8,&
&3)+tt24*H_invDmH(8,1)*H_invDmH(8,3)+tt23*H_invDmH(8,1)*H_invDmH(8&
&,3)+X(3,7)*X(3,8)*H_invDmH(7,1)*H_invDmH(8,3)+X(2,7)*X(2,8)*H_inv&
&DmH(7,1)*H_invDmH(8,3)+X(1,7)*X(1,8)*H_invDmH(7,1)*H_invDmH(8,3)+&
&X(3,6)*X(3,8)*H_invDmH(6,1)*H_invDmH(8,3)+X(2,6)*X(2,8)*H_invDmH(&
&6,1)*H_invDmH(8,3)+X(1,6)*X(1,8)*H_invDmH(6,1)*H_invDmH(8,3)+X(3,&
&5)*X(3,8)*H_invDmH(5,1)*H_invDmH(8,3)+X(2,5)*X(2,8)*H_invDmH(5,1)&
&*H_invDmH(8,3)+X(1,5)*X(1,8)*H_invDmH(5,1)*H_invDmH(8,3)+X(3,4)*X&
&(3,8)*H_invDmH(4,1)*H_invDmH(8,3)+X(2,4)*X(2,8)*H_invDmH(4,1)*H_i&
&nvDmH(8,3)+X(1,4)*X(1,8)*H_invDmH(4,1)*H_invDmH(8,3)+H_invDmH(3,1&
&)*X(3,3)*X(3,8)*H_invDmH(8,3)+H_invDmH(2,1)*X(3,2)*X(3,8)*H_invDm&
&H(8,3)+H_invDmH(1,1)*X(3,1)*X(3,8)*H_invDmH(8,3)+X(2,3)*X(2,8)*H_&
&invDmH(3,1)*H_invDmH(8,3)+X(1,3)*X(1,8)*H_invDmH(3,1)*H_invDmH(8,&
&3)+H_invDmH(2,1)*X(2,2)*X(2,8)*H_invDmH(8,3)+H_invDmH(1,1)*X(2,1)&
&*X(2,8)*H_invDmH(8,3)+X(1,2)*X(1,8)*H_invDmH(2,1)*H_invDmH(8,3)+H&
&_invDmH(1,1)*X(1,1)*X(1,8)*H_invDmH(8,3)+X(3,7)*X(3,8)*H_invDmH(7&
&,3)*H_invDmH(8,1)+X(2,7)*X(2,8)*H_invDmH(7,3)*H_invDmH(8,1)+X(1,7&
&)*X(1,8)*H_invDmH(7,3)*H_invDmH(8,1)+X(3,6)*X(3,8)*H_invDmH(6,3)*&
&H_invDmH(8,1)+X(2,6)*X(2,8)*H_invDmH(6,3)*H_invDmH(8,1)+X(1,6)*X(&
&1,8)*H_invDmH(6,3)*H_invDmH(8,1)+X(3,5)*X(3,8)*H_invDmH(5,3)*H_in&
&vDmH(8,1)+X(2,5)*X(2,8)*H_invDmH(5,3)*H_invDmH(8,1)+X(1,5)*X(1,8)&
&*H_invDmH(5,3)*H_invDmH(8,1)+X(3,4)*X(3,8)*H_invDmH(4,3)*H_invDmH&
&(8,1)+X(2,4)*X(2,8)*H_invDmH(4,3)*H_invDmH(8,1)+X(1,4)*X(1,8)*H_i&
&nvDmH(4,3)*H_invDmH(8,1)+H_invDmH(3,3)*X(3,3)*X(3,8)*H_invDmH(8,1&
&)+H_invDmH(2,3)*X(3,2)*X(3,8)*H_invDmH(8,1)+H_invDmH(1,3)*X(3,1)*&
&X(3,8)*H_invDmH(8,1)+X(2,3)*X(2,8)*H_invDmH(3,3)*H_invDmH(8,1)+X(&
&1,3)*X(1,8)*H_invDmH(3,3)*H_invDmH(8,1)+X(2,2)*H_invDmH(2,3)*X(2,&
&8)*H_invDmH(8,1)+H_invDmH(1,3)*X(2,1)*X(2,8)*H_invDmH(8,1)+X(1,2)&
&*X(1,8)*H_invDmH(2,3)*H_invDmH(8,1)+X(1,1)*H_invDmH(1,3)*X(1,8)*H&
&_invDmH(8,1)+tt22*H_invDmH(7,1)*H_invDmH(7,3)+tt21*H_invDmH(7,1)*&
&H_invDmH(7,3)+tt20*H_invDmH(7,1)*H_invDmH(7,3)+X(3,6)*X(3,7)*H_in&
&vDmH(6,1)*H_invDmH(7,3)+X(2,6)*X(2,7)*H_invDmH(6,1)*H_invDmH(7,3)&
&+X(1,6)*X(1,7)*H_invDmH(6,1)*H_invDmH(7,3)+X(3,5)*X(3,7)*H_invDmH&
&(5,1)*H_invDmH(7,3)+X(2,5)*X(2,7)*H_invDmH(5,1)*H_invDmH(7,3)+X(1&
&,5)*X(1,7)*H_invDmH(5,1)*H_invDmH(7,3)+X(3,4)*X(3,7)*H_invDmH(4,1&
&)*H_invDmH(7,3)+X(2,4)*X(2,7)*H_invDmH(4,1)*H_invDmH(7,3)+X(1,4)*&
&X(1,7)*H_invDmH(4,1)*H_invDmH(7,3)+H_invDmH(3,1)*X(3,3)*X(3,7)*H_&
&invDmH(7,3)+H_invDmH(2,1)*X(3,2)*X(3,7)*H_invDmH(7,3)+H_invDmH(1,&
&1)*X(3,1)*X(3,7)*H_invDmH(7,3)+X(2,3)*X(2,7)*H_invDmH(3,1)*H_invD&
&mH(7,3)+X(1,3)*X(1,7)*H_invDmH(3,1)*H_invDmH(7,3)+H_invDmH(2,1)*X&
&(2,2)*X(2,7)*H_invDmH(7,3)+H_invDmH(1,1)*X(2,1)*X(2,7)*H_invDmH(7&
&,3)+X(1,2)*X(1,7)*H_invDmH(2,1)*H_invDmH(7,3)+H_invDmH(1,1)*X(1,1&
&)*X(1,7)*H_invDmH(7,3)+X(3,6)*X(3,7)*H_invDmH(6,3)*H_invDmH(7,1)+&
&X(2,6)*X(2,7)*H_invDmH(6,3)*H_invDmH(7,1)+X(1,6)*X(1,7)*H_invDmH(&
&6,3)*H_invDmH(7,1)+X(3,5)*X(3,7)*H_invDmH(5,3)*H_invDmH(7,1)+X(2,&
&5)*X(2,7)*H_invDmH(5,3)*H_invDmH(7,1)+X(1,5)*X(1,7)*H_invDmH(5,3)&
&*H_invDmH(7,1)+X(3,4)*X(3,7)*H_invDmH(4,3)*H_invDmH(7,1)+X(2,4)*X&
&(2,7)*H_invDmH(4,3)*H_invDmH(7,1)+X(1,4)*X(1,7)*H_invDmH(4,3)*H_i&
&nvDmH(7,1)+H_invDmH(3,3)*X(3,3)*X(3,7)*H_invDmH(7,1)+H_invDmH(2,3&
&)*X(3,2)*X(3,7)*H_invDmH(7,1)+H_invDmH(1,3)*X(3,1)*X(3,7)*H_invDm&
&H(7,1)+X(2,3)*X(2,7)*H_invDmH(3,3)*H_invDmH(7,1)+X(1,3)*X(1,7)*H_&
&invDmH(3,3)*H_invDmH(7,1)+X(2,2)*H_invDmH(2,3)*X(2,7)*H_invDmH(7,&
&1)+H_invDmH(1,3)*X(2,1)*X(2,7)*H_invDmH(7,1)+X(1,2)*X(1,7)*H_invD&
&mH(2,3)*H_invDmH(7,1)+X(1,1)*H_invDmH(1,3)*X(1,7)*H_invDmH(7,1)+t&
&t19*H_invDmH(6,1)*H_invDmH(6,3)+tt18*H_invDmH(6,1)*H_invDmH(6,3)+&
&tt17*H_invDmH(6,1)*H_invDmH(6,3)+X(3,5)*X(3,6)*H_invDmH(5,1)*H_in&
&vDmH(6,3)+X(2,5)*X(2,6)*H_invDmH(5,1)*H_invDmH(6,3)+X(1,5)*X(1,6)&
&*H_invDmH(5,1)*H_invDmH(6,3)+X(3,4)*X(3,6)*H_invDmH(4,1)*H_invDmH&
&(6,3)+X(2,4)*X(2,6)*H_invDmH(4,1)*H_invDmH(6,3)+X(1,4)*X(1,6)*H_i&
&nvDmH(4,1)*H_invDmH(6,3)+H_invDmH(3,1)*X(3,3)*X(3,6)*H_invDmH(6,3&
&)+H_invDmH(2,1)*X(3,2)*X(3,6)*H_invDmH(6,3)+H_invDmH(1,1)*X(3,1)*&
&X(3,6)*H_invDmH(6,3)+X(2,3)*X(2,6)*H_invDmH(3,1)*H_invDmH(6,3)+X(&
&1,3)*X(1,6)*H_invDmH(3,1)*H_invDmH(6,3)+H_invDmH(2,1)*X(2,2)*X(2,&
&6)*H_invDmH(6,3)+H_invDmH(1,1)*X(2,1)*X(2,6)*H_invDmH(6,3)+X(1,2)&
&*X(1,6)*H_invDmH(2,1)*H_invDmH(6,3)+H_invDmH(1,1)*X(1,1)*X(1,6)*H&
&_invDmH(6,3)+X(3,5)*X(3,6)*H_invDmH(5,3)*H_invDmH(6,1)+X(2,5)*X(2&
&,6)*H_invDmH(5,3)*H_invDmH(6,1)+X(1,5)*X(1,6)*H_invDmH(5,3)*H_inv&
&DmH(6,1)+X(3,4)*X(3,6)*H_invDmH(4,3)*H_invDmH(6,1)+X(2,4)*X(2,6)*&
&H_invDmH(4,3)*H_invDmH(6,1)+X(1,4)*X(1,6)*H_invDmH(4,3)*H_invDmH(&
&6,1)+H_invDmH(3,3)*X(3,3)*X(3,6)*H_invDmH(6,1)+H_invDmH(2,3)*X(3,&
&2)*X(3,6)*H_invDmH(6,1)+H_invDmH(1,3)*X(3,1)*X(3,6)*H_invDmH(6,1)&
&+X(2,3)*X(2,6)*H_invDmH(3,3)*H_invDmH(6,1)+X(1,3)*X(1,6)*H_invDmH&
&(3,3)*H_invDmH(6,1)+X(2,2)*H_invDmH(2,3)*X(2,6)*H_invDmH(6,1)+H_i&
&nvDmH(1,3)*X(2,1)*X(2,6)*H_invDmH(6,1)+X(1,2)*X(1,6)*H_invDmH(2,3&
&)*H_invDmH(6,1)+X(1,1)*H_invDmH(1,3)*X(1,6)*H_invDmH(6,1)+tt16*H_&
&invDmH(5,1)*H_invDmH(5,3)+tt15*H_invDmH(5,1)*H_invDmH(5,3)+tt14*H&
&_invDmH(5,1)*H_invDmH(5,3)+X(3,4)*X(3,5)*H_invDmH(4,1)*H_invDmH(5&
&,3)+X(2,4)*X(2,5)*H_invDmH(4,1)*H_invDmH(5,3)+X(1,4)*X(1,5)*H_inv&
&DmH(4,1)*H_invDmH(5,3)+H_invDmH(3,1)*X(3,3)*X(3,5)*H_invDmH(5,3)+&
&H_invDmH(2,1)*X(3,2)*X(3,5)*H_invDmH(5,3)+H_invDmH(1,1)*X(3,1)*X(&
&3,5)*H_invDmH(5,3)+X(2,3)*X(2,5)*H_invDmH(3,1)*H_invDmH(5,3)+X(1,&
&3)*X(1,5)*H_invDmH(3,1)*H_invDmH(5,3)+H_invDmH(2,1)*X(2,2)*X(2,5)&
&*H_invDmH(5,3)+H_invDmH(1,1)*X(2,1)*X(2,5)*H_invDmH(5,3)+X(1,2)*X&
&(1,5)*H_invDmH(2,1)*H_invDmH(5,3)+H_invDmH(1,1)*X(1,1)*X(1,5)*H_i&
&nvDmH(5,3)+X(3,4)*X(3,5)*H_invDmH(4,3)*H_invDmH(5,1)+X(2,4)*X(2,5&
&)*H_invDmH(4,3)*H_invDmH(5,1)+X(1,4)*X(1,5)*H_invDmH(4,3)*H_invDm&
&H(5,1)+H_invDmH(3,3)*X(3,3)*X(3,5)*H_invDmH(5,1)+H_invDmH(2,3)*X(&
&3,2)*X(3,5)*H_invDmH(5,1)+H_invDmH(1,3)*X(3,1)*X(3,5)*H_invDmH(5,&
&1)+X(2,3)*X(2,5)*H_invDmH(3,3)*H_invDmH(5,1)+X(1,3)*X(1,5)*H_invD&
&mH(3,3)*H_invDmH(5,1)+X(2,2)*H_invDmH(2,3)*X(2,5)*H_invDmH(5,1)+H&
&_invDmH(1,3)*X(2,1)*X(2,5)*H_invDmH(5,1)+X(1,2)*X(1,5)*H_invDmH(2&
&,3)*H_invDmH(5,1)+X(1,1)*H_invDmH(1,3)*X(1,5)*H_invDmH(5,1)+tt13*&
&H_invDmH(4,1)*H_invDmH(4,3)+tt12*H_invDmH(4,1)*H_invDmH(4,3)+tt11&
&*H_invDmH(4,1)*H_invDmH(4,3)+H_invDmH(3,1)*X(3,3)*X(3,4)*H_invDmH&
&(4,3)+H_invDmH(2,1)*X(3,2)*X(3,4)*H_invDmH(4,3)+H_invDmH(1,1)*X(3&
&,1)*X(3,4)*H_invDmH(4,3)+X(2,3)*X(2,4)*H_invDmH(3,1)*H_invDmH(4,3&
&)+X(1,3)*X(1,4)*H_invDmH(3,1)*H_invDmH(4,3)+H_invDmH(2,1)*X(2,2)*&
&X(2,4)*H_invDmH(4,3)+H_invDmH(1,1)*X(2,1)*X(2,4)*H_invDmH(4,3)+X(&
&1,2)*X(1,4)*H_invDmH(2,1)*H_invDmH(4,3)+H_invDmH(1,1)*X(1,1)*X(1,&
&4)*H_invDmH(4,3)+H_invDmH(3,3)*X(3,3)*X(3,4)*H_invDmH(4,1)+H_invD&
&mH(2,3)*X(3,2)*X(3,4)*H_invDmH(4,1)+H_invDmH(1,3)*X(3,1)*X(3,4)*H&
&_invDmH(4,1)+X(2,3)*X(2,4)*H_invDmH(3,3)*H_invDmH(4,1)+X(1,3)*X(1&
&,4)*H_invDmH(3,3)*H_invDmH(4,1)+X(2,2)*H_invDmH(2,3)*X(2,4)*H_inv&
&DmH(4,1)+H_invDmH(1,3)*X(2,1)*X(2,4)*H_invDmH(4,1)+X(1,2)*X(1,4)*&
&H_invDmH(2,3)*H_invDmH(4,1)+X(1,1)*H_invDmH(1,3)*X(1,4)*H_invDmH(&
&4,1)+H_invDmH(3,1)*H_invDmH(3,3)*tt10+H_invDmH(2,1)*X(3,2)*H_invD&
&mH(3,3)*X(3,3)+H_invDmH(1,1)*X(3,1)*H_invDmH(3,3)*X(3,3)+H_invDmH&
&(2,3)*H_invDmH(3,1)*X(3,2)*X(3,3)+H_invDmH(1,3)*H_invDmH(3,1)*X(3&
&,1)*X(3,3)+tt8*H_invDmH(3,1)*H_invDmH(3,3)+tt7*H_invDmH(3,1)*H_in&
&vDmH(3,3)+H_invDmH(2,1)*X(2,2)*X(2,3)*H_invDmH(3,3)+H_invDmH(1,1)&
&*X(2,1)*X(2,3)*H_invDmH(3,3)+X(1,2)*X(1,3)*H_invDmH(2,1)*H_invDmH&
&(3,3)+H_invDmH(1,1)*X(1,1)*X(1,3)*H_invDmH(3,3)+H_invDmH(2,1)*H_i&
&nvDmH(2,3)*tt9+H_invDmH(1,1)*H_invDmH(2,3)*X(3,1)*X(3,2)+H_invDmH&
&(1,3)*H_invDmH(2,1)*X(3,1)*X(3,2)+H_invDmH(1,1)*H_invDmH(1,3)*tt6&
&+X(2,2)*H_invDmH(2,3)*X(2,3)*H_invDmH(3,1)+H_invDmH(1,3)*X(2,1)*X&
&(2,3)*H_invDmH(3,1)+X(1,2)*X(1,3)*H_invDmH(2,3)*H_invDmH(3,1)+X(1&
&,1)*H_invDmH(1,3)*X(1,3)*H_invDmH(3,1)+H_invDmH(2,1)*tt5*H_invDmH&
&(2,3)+H_invDmH(1,1)*X(2,1)*X(2,2)*H_invDmH(2,3)+tt4*H_invDmH(2,1)&
&*H_invDmH(2,3)+H_invDmH(1,1)*X(1,1)*X(1,2)*H_invDmH(2,3)+H_invDmH&
&(1,3)*H_invDmH(2,1)*X(2,1)*X(2,2)+H_invDmH(1,1)*H_invDmH(1,3)*tt3&
&+X(1,1)*X(1,2)*H_invDmH(1,3)*H_invDmH(2,1)+H_invDmH(1,1)*tt2*H_in&
&vDmH(1,3))**2/2.0E+0+tt26**2/4.0E+0+(tt25*H_invDmH(8,1)*H_invDmH(&
&8,2)+tt24*H_invDmH(8,1)*H_invDmH(8,2)+tt23*H_invDmH(8,1)*H_invDmH&
&(8,2)+X(3,7)*X(3,8)*H_invDmH(7,1)*H_invDmH(8,2)+X(2,7)*X(2,8)*H_i&
&nvDmH(7,1)*H_invDmH(8,2)+X(1,7)*X(1,8)*H_invDmH(7,1)*H_invDmH(8,2&
&)+X(3,6)*X(3,8)*H_invDmH(6,1)*H_invDmH(8,2)+X(2,6)*X(2,8)*H_invDm&
&H(6,1)*H_invDmH(8,2)+X(1,6)*X(1,8)*H_invDmH(6,1)*H_invDmH(8,2)+X(&
&3,5)*X(3,8)*H_invDmH(5,1)*H_invDmH(8,2)+X(2,5)*X(2,8)*H_invDmH(5,&
&1)*H_invDmH(8,2)+X(1,5)*X(1,8)*H_invDmH(5,1)*H_invDmH(8,2)+X(3,4)&
&*X(3,8)*H_invDmH(4,1)*H_invDmH(8,2)+X(2,4)*X(2,8)*H_invDmH(4,1)*H&
&_invDmH(8,2)+X(1,4)*X(1,8)*H_invDmH(4,1)*H_invDmH(8,2)+H_invDmH(3&
&,1)*X(3,3)*X(3,8)*H_invDmH(8,2)+H_invDmH(2,1)*X(3,2)*X(3,8)*H_inv&
&DmH(8,2)+H_invDmH(1,1)*X(3,1)*X(3,8)*H_invDmH(8,2)+X(2,3)*X(2,8)*&
&H_invDmH(3,1)*H_invDmH(8,2)+X(1,3)*X(1,8)*H_invDmH(3,1)*H_invDmH(&
&8,2)+H_invDmH(2,1)*X(2,2)*X(2,8)*H_invDmH(8,2)+H_invDmH(1,1)*X(2,&
&1)*X(2,8)*H_invDmH(8,2)+X(1,2)*X(1,8)*H_invDmH(2,1)*H_invDmH(8,2)&
&+H_invDmH(1,1)*X(1,1)*X(1,8)*H_invDmH(8,2)+X(3,7)*X(3,8)*H_invDmH&
&(7,2)*H_invDmH(8,1)+X(2,7)*X(2,8)*H_invDmH(7,2)*H_invDmH(8,1)+X(1&
&,7)*X(1,8)*H_invDmH(7,2)*H_invDmH(8,1)+X(3,6)*X(3,8)*H_invDmH(6,2&
&)*H_invDmH(8,1)+X(2,6)*X(2,8)*H_invDmH(6,2)*H_invDmH(8,1)+X(1,6)*&
&X(1,8)*H_invDmH(6,2)*H_invDmH(8,1)+X(3,5)*X(3,8)*H_invDmH(5,2)*H_&
&invDmH(8,1)+X(2,5)*X(2,8)*H_invDmH(5,2)*H_invDmH(8,1)+X(1,5)*X(1,&
&8)*H_invDmH(5,2)*H_invDmH(8,1)+X(3,4)*X(3,8)*H_invDmH(4,2)*H_invD&
&mH(8,1)+X(2,4)*X(2,8)*H_invDmH(4,2)*H_invDmH(8,1)+X(1,4)*X(1,8)*H&
&_invDmH(4,2)*H_invDmH(8,1)+H_invDmH(3,2)*X(3,3)*X(3,8)*H_invDmH(8&
&,1)+H_invDmH(2,2)*X(3,2)*X(3,8)*H_invDmH(8,1)+H_invDmH(1,2)*X(3,1&
&)*X(3,8)*H_invDmH(8,1)+X(2,3)*X(2,8)*H_invDmH(3,2)*H_invDmH(8,1)+&
&X(1,3)*X(1,8)*H_invDmH(3,2)*H_invDmH(8,1)+H_invDmH(2,2)*X(2,2)*X(&
&2,8)*H_invDmH(8,1)+H_invDmH(1,2)*X(2,1)*X(2,8)*H_invDmH(8,1)+X(1,&
&2)*X(1,8)*H_invDmH(2,2)*H_invDmH(8,1)+X(1,1)*H_invDmH(1,2)*X(1,8)&
&*H_invDmH(8,1)+tt22*H_invDmH(7,1)*H_invDmH(7,2)+tt21*H_invDmH(7,1&
&)*H_invDmH(7,2)+tt20*H_invDmH(7,1)*H_invDmH(7,2)+X(3,6)*X(3,7)*H_&
&invDmH(6,1)*H_invDmH(7,2)+X(2,6)*X(2,7)*H_invDmH(6,1)*H_invDmH(7,&
&2)+X(1,6)*X(1,7)*H_invDmH(6,1)*H_invDmH(7,2)+X(3,5)*X(3,7)*H_invD&
&mH(5,1)*H_invDmH(7,2)+X(2,5)*X(2,7)*H_invDmH(5,1)*H_invDmH(7,2)+X&
&(1,5)*X(1,7)*H_invDmH(5,1)*H_invDmH(7,2)+X(3,4)*X(3,7)*H_invDmH(4&
&,1)*H_invDmH(7,2)+X(2,4)*X(2,7)*H_invDmH(4,1)*H_invDmH(7,2)+X(1,4&
&)*X(1,7)*H_invDmH(4,1)*H_invDmH(7,2)+H_invDmH(3,1)*X(3,3)*X(3,7)*&
&H_invDmH(7,2)+H_invDmH(2,1)*X(3,2)*X(3,7)*H_invDmH(7,2)+H_invDmH(&
&1,1)*X(3,1)*X(3,7)*H_invDmH(7,2)+X(2,3)*X(2,7)*H_invDmH(3,1)*H_in&
&vDmH(7,2)+X(1,3)*X(1,7)*H_invDmH(3,1)*H_invDmH(7,2)+H_invDmH(2,1)&
&*X(2,2)*X(2,7)*H_invDmH(7,2)+H_invDmH(1,1)*X(2,1)*X(2,7)*H_invDmH&
&(7,2)+X(1,2)*X(1,7)*H_invDmH(2,1)*H_invDmH(7,2)+H_invDmH(1,1)*X(1&
&,1)*X(1,7)*H_invDmH(7,2)+X(3,6)*X(3,7)*H_invDmH(6,2)*H_invDmH(7,1&
&)+X(2,6)*X(2,7)*H_invDmH(6,2)*H_invDmH(7,1)+X(1,6)*X(1,7)*H_invDm&
&H(6,2)*H_invDmH(7,1)+X(3,5)*X(3,7)*H_invDmH(5,2)*H_invDmH(7,1)+X(&
&2,5)*X(2,7)*H_invDmH(5,2)*H_invDmH(7,1)+X(1,5)*X(1,7)*H_invDmH(5,&
&2)*H_invDmH(7,1)+X(3,4)*X(3,7)*H_invDmH(4,2)*H_invDmH(7,1)+X(2,4)&
&*X(2,7)*H_invDmH(4,2)*H_invDmH(7,1)+X(1,4)*X(1,7)*H_invDmH(4,2)*H&
&_invDmH(7,1)+H_invDmH(3,2)*X(3,3)*X(3,7)*H_invDmH(7,1)+H_invDmH(2&
&,2)*X(3,2)*X(3,7)*H_invDmH(7,1)+H_invDmH(1,2)*X(3,1)*X(3,7)*H_inv&
&DmH(7,1)+X(2,3)*X(2,7)*H_invDmH(3,2)*H_invDmH(7,1)+X(1,3)*X(1,7)*&
&H_invDmH(3,2)*H_invDmH(7,1)+H_invDmH(2,2)*X(2,2)*X(2,7)*H_invDmH(&
&7,1)+H_invDmH(1,2)*X(2,1)*X(2,7)*H_invDmH(7,1)+X(1,2)*X(1,7)*H_in&
&vDmH(2,2)*H_invDmH(7,1)+X(1,1)*H_invDmH(1,2)*X(1,7)*H_invDmH(7,1)&
&+tt19*H_invDmH(6,1)*H_invDmH(6,2)+tt18*H_invDmH(6,1)*H_invDmH(6,2&
&)+tt17*H_invDmH(6,1)*H_invDmH(6,2)+X(3,5)*X(3,6)*H_invDmH(5,1)*H_&
&invDmH(6,2)+X(2,5)*X(2,6)*H_invDmH(5,1)*H_invDmH(6,2)+X(1,5)*X(1,&
&6)*H_invDmH(5,1)*H_invDmH(6,2)+X(3,4)*X(3,6)*H_invDmH(4,1)*H_invD&
&mH(6,2)+X(2,4)*X(2,6)*H_invDmH(4,1)*H_invDmH(6,2)+X(1,4)*X(1,6)*H&
&_invDmH(4,1)*H_invDmH(6,2)+H_invDmH(3,1)*X(3,3)*X(3,6)*H_invDmH(6&
&,2)+H_invDmH(2,1)*X(3,2)*X(3,6)*H_invDmH(6,2)+H_invDmH(1,1)*X(3,1&
&)*X(3,6)*H_invDmH(6,2)+X(2,3)*X(2,6)*H_invDmH(3,1)*H_invDmH(6,2)+&
&X(1,3)*X(1,6)*H_invDmH(3,1)*H_invDmH(6,2)+H_invDmH(2,1)*X(2,2)*X(&
&2,6)*H_invDmH(6,2)+H_invDmH(1,1)*X(2,1)*X(2,6)*H_invDmH(6,2)+X(1,&
&2)*X(1,6)*H_invDmH(2,1)*H_invDmH(6,2)+H_invDmH(1,1)*X(1,1)*X(1,6)&
&*H_invDmH(6,2)+X(3,5)*X(3,6)*H_invDmH(5,2)*H_invDmH(6,1)+X(2,5)*X&
&(2,6)*H_invDmH(5,2)*H_invDmH(6,1)+X(1,5)*X(1,6)*H_invDmH(5,2)*H_i&
&nvDmH(6,1)+X(3,4)*X(3,6)*H_invDmH(4,2)*H_invDmH(6,1)+X(2,4)*X(2,6&
&)*H_invDmH(4,2)*H_invDmH(6,1)+X(1,4)*X(1,6)*H_invDmH(4,2)*H_invDm&
&H(6,1)+H_invDmH(3,2)*X(3,3)*X(3,6)*H_invDmH(6,1)+H_invDmH(2,2)*X(&
&3,2)*X(3,6)*H_invDmH(6,1)+H_invDmH(1,2)*X(3,1)*X(3,6)*H_invDmH(6,&
&1)+X(2,3)*X(2,6)*H_invDmH(3,2)*H_invDmH(6,1)+X(1,3)*X(1,6)*H_invD&
&mH(3,2)*H_invDmH(6,1)+H_invDmH(2,2)*X(2,2)*X(2,6)*H_invDmH(6,1)+H&
&_invDmH(1,2)*X(2,1)*X(2,6)*H_invDmH(6,1)+X(1,2)*X(1,6)*H_invDmH(2&
&,2)*H_invDmH(6,1)+X(1,1)*H_invDmH(1,2)*X(1,6)*H_invDmH(6,1)+tt16*&
&H_invDmH(5,1)*H_invDmH(5,2)+tt15*H_invDmH(5,1)*H_invDmH(5,2)+tt14&
&*H_invDmH(5,1)*H_invDmH(5,2)+X(3,4)*X(3,5)*H_invDmH(4,1)*H_invDmH&
&(5,2)+X(2,4)*X(2,5)*H_invDmH(4,1)*H_invDmH(5,2)+X(1,4)*X(1,5)*H_i&
&nvDmH(4,1)*H_invDmH(5,2)+H_invDmH(3,1)*X(3,3)*X(3,5)*H_invDmH(5,2&
&)+H_invDmH(2,1)*X(3,2)*X(3,5)*H_invDmH(5,2)+H_invDmH(1,1)*X(3,1)*&
&X(3,5)*H_invDmH(5,2)+X(2,3)*X(2,5)*H_invDmH(3,1)*H_invDmH(5,2)+X(&
&1,3)*X(1,5)*H_invDmH(3,1)*H_invDmH(5,2)+H_invDmH(2,1)*X(2,2)*X(2,&
&5)*H_invDmH(5,2)+H_invDmH(1,1)*X(2,1)*X(2,5)*H_invDmH(5,2)+X(1,2)&
&*X(1,5)*H_invDmH(2,1)*H_invDmH(5,2)+H_invDmH(1,1)*X(1,1)*X(1,5)*H&
&_invDmH(5,2)+X(3,4)*X(3,5)*H_invDmH(4,2)*H_invDmH(5,1)+X(2,4)*X(2&
&,5)*H_invDmH(4,2)*H_invDmH(5,1)+X(1,4)*X(1,5)*H_invDmH(4,2)*H_inv&
&DmH(5,1)+H_invDmH(3,2)*X(3,3)*X(3,5)*H_invDmH(5,1)+H_invDmH(2,2)*&
&X(3,2)*X(3,5)*H_invDmH(5,1)+H_invDmH(1,2)*X(3,1)*X(3,5)*H_invDmH(&
&5,1)+X(2,3)*X(2,5)*H_invDmH(3,2)*H_invDmH(5,1)+X(1,3)*X(1,5)*H_in&
&vDmH(3,2)*H_invDmH(5,1)+H_invDmH(2,2)*X(2,2)*X(2,5)*H_invDmH(5,1)&
&+H_invDmH(1,2)*X(2,1)*X(2,5)*H_invDmH(5,1)+X(1,2)*X(1,5)*H_invDmH&
&(2,2)*H_invDmH(5,1)+X(1,1)*H_invDmH(1,2)*X(1,5)*H_invDmH(5,1)+tt1&
&3*H_invDmH(4,1)*H_invDmH(4,2)+tt12*H_invDmH(4,1)*H_invDmH(4,2)+tt&
&11*H_invDmH(4,1)*H_invDmH(4,2)+H_invDmH(3,1)*X(3,3)*X(3,4)*H_invD&
&mH(4,2)+H_invDmH(2,1)*X(3,2)*X(3,4)*H_invDmH(4,2)+H_invDmH(1,1)*X&
&(3,1)*X(3,4)*H_invDmH(4,2)+X(2,3)*X(2,4)*H_invDmH(3,1)*H_invDmH(4&
&,2)+X(1,3)*X(1,4)*H_invDmH(3,1)*H_invDmH(4,2)+H_invDmH(2,1)*X(2,2&
&)*X(2,4)*H_invDmH(4,2)+H_invDmH(1,1)*X(2,1)*X(2,4)*H_invDmH(4,2)+&
&X(1,2)*X(1,4)*H_invDmH(2,1)*H_invDmH(4,2)+H_invDmH(1,1)*X(1,1)*X(&
&1,4)*H_invDmH(4,2)+H_invDmH(3,2)*X(3,3)*X(3,4)*H_invDmH(4,1)+H_in&
&vDmH(2,2)*X(3,2)*X(3,4)*H_invDmH(4,1)+H_invDmH(1,2)*X(3,1)*X(3,4)&
&*H_invDmH(4,1)+X(2,3)*X(2,4)*H_invDmH(3,2)*H_invDmH(4,1)+X(1,3)*X&
&(1,4)*H_invDmH(3,2)*H_invDmH(4,1)+H_invDmH(2,2)*X(2,2)*X(2,4)*H_i&
&nvDmH(4,1)+H_invDmH(1,2)*X(2,1)*X(2,4)*H_invDmH(4,1)+X(1,2)*X(1,4&
&)*H_invDmH(2,2)*H_invDmH(4,1)+X(1,1)*H_invDmH(1,2)*X(1,4)*H_invDm&
&H(4,1)+H_invDmH(3,1)*H_invDmH(3,2)*tt10+H_invDmH(2,1)*H_invDmH(3,&
&2)*X(3,2)*X(3,3)+H_invDmH(2,2)*H_invDmH(3,1)*X(3,2)*X(3,3)+H_invD&
&mH(1,1)*X(3,1)*H_invDmH(3,2)*X(3,3)+H_invDmH(1,2)*H_invDmH(3,1)*X&
&(3,1)*X(3,3)+H_invDmH(2,1)*H_invDmH(2,2)*tt9+H_invDmH(1,1)*H_invD&
&mH(2,2)*X(3,1)*X(3,2)+H_invDmH(1,2)*H_invDmH(2,1)*X(3,1)*X(3,2)+t&
&t8*H_invDmH(3,1)*H_invDmH(3,2)+tt7*H_invDmH(3,1)*H_invDmH(3,2)+H_&
&invDmH(2,1)*X(2,2)*X(2,3)*H_invDmH(3,2)+H_invDmH(1,1)*X(2,1)*X(2,&
&3)*H_invDmH(3,2)+X(1,2)*X(1,3)*H_invDmH(2,1)*H_invDmH(3,2)+H_invD&
&mH(1,1)*X(1,1)*X(1,3)*H_invDmH(3,2)+H_invDmH(1,1)*H_invDmH(1,2)*t&
&t6+H_invDmH(2,2)*X(2,2)*X(2,3)*H_invDmH(3,1)+H_invDmH(1,2)*X(2,1)&
&*X(2,3)*H_invDmH(3,1)+X(1,2)*X(1,3)*H_invDmH(2,2)*H_invDmH(3,1)+X&
&(1,1)*H_invDmH(1,2)*X(1,3)*H_invDmH(3,1)+H_invDmH(2,1)*H_invDmH(2&
&,2)*tt5+H_invDmH(1,1)*X(2,1)*H_invDmH(2,2)*X(2,2)+H_invDmH(1,2)*H&
&_invDmH(2,1)*X(2,1)*X(2,2)+tt4*H_invDmH(2,1)*H_invDmH(2,2)+H_invD&
&mH(1,1)*X(1,1)*X(1,2)*H_invDmH(2,2)+H_invDmH(1,1)*H_invDmH(1,2)*t&
&t3+X(1,1)*H_invDmH(1,2)*X(1,2)*H_invDmH(2,1)+H_invDmH(1,1)*tt2*H_&
&invDmH(1,2))**2/2.0E+0+tt1**2/4.0E+0)
jac(1,2) = (detDmH(1,1)*gw(1,1)*(tt27/2.0E+0+tt26/2.0E+0+tt1/2.0E&
&+0)**2)/2.0E+0
END 
SUBROUTINE vox_neo_at_quadr_mtr(val, mtr, X, H_invDmH, detDmH, gw) 
IMPLICIT NONE 
REAL(KIND=8) val(1, 1) 
REAL(KIND=8) mtr(2, 1) 
REAL(KIND=8) X(3, 8) 
REAL(KIND=8) H_invDmH(8, 3) 
REAL(KIND=8) detDmH(1, 1) 
REAL(KIND=8) gw(1, 1) 
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
tt1 = X(1,8)*H_invDmH(8,1)+X(1,7)*H_invDmH(7,1)+X(1,6)*H_invDmH(6&
&,1)+X(1,5)*H_invDmH(5,1)+X(1,4)*H_invDmH(4,1)+X(1,3)*H_invDmH(3,1&
&)+X(1,2)*H_invDmH(2,1)+H_invDmH(1,1)*X(1,1)
tt2 = X(2,8)*H_invDmH(8,1)+X(2,7)*H_invDmH(7,1)+X(2,6)*H_invDmH(6&
&,1)+X(2,5)*H_invDmH(5,1)+X(2,4)*H_invDmH(4,1)+X(2,3)*H_invDmH(3,1&
&)+H_invDmH(2,1)*X(2,2)+H_invDmH(1,1)*X(2,1)
tt3 = X(3,8)*H_invDmH(8,1)+X(3,7)*H_invDmH(7,1)+X(3,6)*H_invDmH(6&
&,1)+X(3,5)*H_invDmH(5,1)+X(3,4)*H_invDmH(4,1)+H_invDmH(3,1)*X(3,3&
&)+H_invDmH(2,1)*X(3,2)+H_invDmH(1,1)*X(3,1)
tt4 = X(1,8)*H_invDmH(8,2)+X(1,7)*H_invDmH(7,2)+X(1,6)*H_invDmH(6&
&,2)+X(1,5)*H_invDmH(5,2)+X(1,4)*H_invDmH(4,2)+X(1,3)*H_invDmH(3,2&
&)+X(1,2)*H_invDmH(2,2)+X(1,1)*H_invDmH(1,2)
tt5 = X(2,8)*H_invDmH(8,2)+X(2,7)*H_invDmH(7,2)+X(2,6)*H_invDmH(6&
&,2)+X(2,5)*H_invDmH(5,2)+X(2,4)*H_invDmH(4,2)+X(2,3)*H_invDmH(3,2&
&)+H_invDmH(2,2)*X(2,2)+H_invDmH(1,2)*X(2,1)
tt6 = X(3,8)*H_invDmH(8,2)+X(3,7)*H_invDmH(7,2)+X(3,6)*H_invDmH(6&
&,2)+X(3,5)*H_invDmH(5,2)+X(3,4)*H_invDmH(4,2)+H_invDmH(3,2)*X(3,3&
&)+H_invDmH(2,2)*X(3,2)+H_invDmH(1,2)*X(3,1)
tt7 = X(1,8)*H_invDmH(8,3)+X(1,7)*H_invDmH(7,3)+X(1,6)*H_invDmH(6&
&,3)+X(1,5)*H_invDmH(5,3)+X(1,4)*H_invDmH(4,3)+X(1,3)*H_invDmH(3,3&
&)+X(1,2)*H_invDmH(2,3)+X(1,1)*H_invDmH(1,3)
tt8 = X(2,8)*H_invDmH(8,3)+X(2,7)*H_invDmH(7,3)+X(2,6)*H_invDmH(6&
&,3)+X(2,5)*H_invDmH(5,3)+X(2,4)*H_invDmH(4,3)+X(2,3)*H_invDmH(3,3&
&)+X(2,2)*H_invDmH(2,3)+H_invDmH(1,3)*X(2,1)
tt9 = X(3,8)*H_invDmH(8,3)+X(3,7)*H_invDmH(7,3)+X(3,6)*H_invDmH(6&
&,3)+X(3,5)*H_invDmH(5,3)+X(3,4)*H_invDmH(4,3)+H_invDmH(3,3)*X(3,3&
&)+H_invDmH(2,3)*X(3,2)+H_invDmH(1,3)*X(3,1)
tt10 = log(tt1*(tt5*tt9-tt6*tt8)-tt4*(tt2*tt9-tt3*tt8)+(tt2*tt6-t&
&t3*tt5)*tt7)
val(1,1) = 5.0E-1*detDmH(1,1)*gw(1,1)*(mtr(2,1)*tt10**2+mtr(1,1)*&
&((-2*tt10)+tt9**2+tt8**2+tt7**2+tt6**2+tt5**2+tt4**2+tt3**2+tt2**&
&2+tt1**2-3))
END 
SUBROUTINE vox_neo_at_quadr_mtr_jac(jac, mtr, X, H_invDmH, detDmH, gw) 
IMPLICIT NONE 
REAL(KIND=8) jac(1, 2) 
REAL(KIND=8) mtr(2, 1) 
REAL(KIND=8) X(3, 8) 
REAL(KIND=8) H_invDmH(8, 3) 
REAL(KIND=8) detDmH(1, 1) 
REAL(KIND=8) gw(1, 1) 
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
tt1 = X(1,8)*H_invDmH(8,1)+X(1,7)*H_invDmH(7,1)+X(1,6)*H_invDmH(6&
&,1)+X(1,5)*H_invDmH(5,1)+X(1,4)*H_invDmH(4,1)+X(1,3)*H_invDmH(3,1&
&)+X(1,2)*H_invDmH(2,1)+H_invDmH(1,1)*X(1,1)
tt2 = X(2,8)*H_invDmH(8,1)+X(2,7)*H_invDmH(7,1)+X(2,6)*H_invDmH(6&
&,1)+X(2,5)*H_invDmH(5,1)+X(2,4)*H_invDmH(4,1)+X(2,3)*H_invDmH(3,1&
&)+H_invDmH(2,1)*X(2,2)+H_invDmH(1,1)*X(2,1)
tt3 = X(3,8)*H_invDmH(8,1)+X(3,7)*H_invDmH(7,1)+X(3,6)*H_invDmH(6&
&,1)+X(3,5)*H_invDmH(5,1)+X(3,4)*H_invDmH(4,1)+H_invDmH(3,1)*X(3,3&
&)+H_invDmH(2,1)*X(3,2)+H_invDmH(1,1)*X(3,1)
tt4 = X(1,8)*H_invDmH(8,2)+X(1,7)*H_invDmH(7,2)+X(1,6)*H_invDmH(6&
&,2)+X(1,5)*H_invDmH(5,2)+X(1,4)*H_invDmH(4,2)+X(1,3)*H_invDmH(3,2&
&)+X(1,2)*H_invDmH(2,2)+X(1,1)*H_invDmH(1,2)
tt5 = X(2,8)*H_invDmH(8,2)+X(2,7)*H_invDmH(7,2)+X(2,6)*H_invDmH(6&
&,2)+X(2,5)*H_invDmH(5,2)+X(2,4)*H_invDmH(4,2)+X(2,3)*H_invDmH(3,2&
&)+H_invDmH(2,2)*X(2,2)+H_invDmH(1,2)*X(2,1)
tt6 = X(3,8)*H_invDmH(8,2)+X(3,7)*H_invDmH(7,2)+X(3,6)*H_invDmH(6&
&,2)+X(3,5)*H_invDmH(5,2)+X(3,4)*H_invDmH(4,2)+H_invDmH(3,2)*X(3,3&
&)+H_invDmH(2,2)*X(3,2)+H_invDmH(1,2)*X(3,1)
tt7 = X(1,8)*H_invDmH(8,3)+X(1,7)*H_invDmH(7,3)+X(1,6)*H_invDmH(6&
&,3)+X(1,5)*H_invDmH(5,3)+X(1,4)*H_invDmH(4,3)+X(1,3)*H_invDmH(3,3&
&)+X(1,2)*H_invDmH(2,3)+X(1,1)*H_invDmH(1,3)
tt8 = X(2,8)*H_invDmH(8,3)+X(2,7)*H_invDmH(7,3)+X(2,6)*H_invDmH(6&
&,3)+X(2,5)*H_invDmH(5,3)+X(2,4)*H_invDmH(4,3)+X(2,3)*H_invDmH(3,3&
&)+X(2,2)*H_invDmH(2,3)+H_invDmH(1,3)*X(2,1)
tt9 = X(3,8)*H_invDmH(8,3)+X(3,7)*H_invDmH(7,3)+X(3,6)*H_invDmH(6&
&,3)+X(3,5)*H_invDmH(5,3)+X(3,4)*H_invDmH(4,3)+H_invDmH(3,3)*X(3,3&
&)+H_invDmH(2,3)*X(3,2)+H_invDmH(1,3)*X(3,1)
tt10 = log(tt1*(tt5*tt9-tt6*tt8)-tt4*(tt2*tt9-tt3*tt8)+(tt2*tt6-t&
&t3*tt5)*tt7)
jac(1,1) = 5.0E-1*detDmH(1,1)*gw(1,1)*((-2*tt10)+tt9**2+tt8**2+tt&
&7**2+tt6**2+tt5**2+tt4**2+tt3**2+tt2**2+tt1**2-3)
jac(1,2) = 5.0E-1*detDmH(1,1)*gw(1,1)*tt10**2
END 
SUBROUTINE vox_linear_at_quadr_mtr(val, mtr, X, H_invDmH, detDmH, gw) 
IMPLICIT NONE 
REAL(KIND=8) val(1, 1) 
REAL(KIND=8) mtr(2, 1) 
REAL(KIND=8) X(3, 8) 
REAL(KIND=8) H_invDmH(8, 3) 
REAL(KIND=8) detDmH(1, 1) 
REAL(KIND=8) gw(1, 1) 
val(1,1) = detDmH(1,1)*gw(1,1)*((mtr(2,1)*((2*X(3,8)*H_invDmH(8,3&
&)+2*X(3,7)*H_invDmH(7,3)+2*X(3,6)*H_invDmH(6,3)+2*X(3,5)*H_invDmH&
&(5,3)+2*X(3,4)*H_invDmH(4,3)+2*H_invDmH(3,3)*X(3,3)+2*H_invDmH(2,&
&3)*X(3,2)+2*H_invDmH(1,3)*X(3,1))/2.0E+0+(2*X(2,8)*H_invDmH(8,2)+&
&2*X(2,7)*H_invDmH(7,2)+2*X(2,6)*H_invDmH(6,2)+2*X(2,5)*H_invDmH(5&
&,2)+2*X(2,4)*H_invDmH(4,2)+2*X(2,3)*H_invDmH(3,2)+2*H_invDmH(2,2)&
&*X(2,2)+2*H_invDmH(1,2)*X(2,1))/2.0E+0+(2*X(1,8)*H_invDmH(8,1)+2*&
&X(1,7)*H_invDmH(7,1)+2*X(1,6)*H_invDmH(6,1)+2*X(1,5)*H_invDmH(5,1&
&)+2*X(1,4)*H_invDmH(4,1)+2*X(1,3)*H_invDmH(3,1)+2*X(1,2)*H_invDmH&
&(2,1)+2*H_invDmH(1,1)*X(1,1))/2.0E+0-3)**2)/2.0E+0+mtr(1,1)*((X(3&
&,8)*H_invDmH(8,3)+X(3,7)*H_invDmH(7,3)+X(3,6)*H_invDmH(6,3)+X(3,5&
&)*H_invDmH(5,3)+X(3,4)*H_invDmH(4,3)+H_invDmH(3,3)*X(3,3)+H_invDm&
&H(2,3)*X(3,2)+H_invDmH(1,3)*X(3,1)-1)**2+(X(2,8)*H_invDmH(8,3)+X(&
&3,8)*H_invDmH(8,2)+X(2,7)*H_invDmH(7,3)+X(3,7)*H_invDmH(7,2)+X(2,&
&6)*H_invDmH(6,3)+X(3,6)*H_invDmH(6,2)+X(2,5)*H_invDmH(5,3)+X(3,5)&
&*H_invDmH(5,2)+X(2,4)*H_invDmH(4,3)+X(3,4)*H_invDmH(4,2)+H_invDmH&
&(3,2)*X(3,3)+X(2,3)*H_invDmH(3,3)+H_invDmH(2,2)*X(3,2)+H_invDmH(1&
&,2)*X(3,1)+X(2,2)*H_invDmH(2,3)+H_invDmH(1,3)*X(2,1))**2/2.0E+0+(&
&X(1,8)*H_invDmH(8,3)+X(3,8)*H_invDmH(8,1)+X(1,7)*H_invDmH(7,3)+X(&
&3,7)*H_invDmH(7,1)+X(1,6)*H_invDmH(6,3)+X(3,6)*H_invDmH(6,1)+X(1,&
&5)*H_invDmH(5,3)+X(3,5)*H_invDmH(5,1)+X(1,4)*H_invDmH(4,3)+X(3,4)&
&*H_invDmH(4,1)+H_invDmH(3,1)*X(3,3)+X(1,3)*H_invDmH(3,3)+H_invDmH&
&(2,1)*X(3,2)+H_invDmH(1,1)*X(3,1)+X(1,2)*H_invDmH(2,3)+X(1,1)*H_i&
&nvDmH(1,3))**2/2.0E+0+(X(2,8)*H_invDmH(8,2)+X(2,7)*H_invDmH(7,2)+&
&X(2,6)*H_invDmH(6,2)+X(2,5)*H_invDmH(5,2)+X(2,4)*H_invDmH(4,2)+X(&
&2,3)*H_invDmH(3,2)+H_invDmH(2,2)*X(2,2)+H_invDmH(1,2)*X(2,1)-1)**&
&2+(X(1,8)*H_invDmH(8,2)+X(2,8)*H_invDmH(8,1)+X(1,7)*H_invDmH(7,2)&
&+X(2,7)*H_invDmH(7,1)+X(1,6)*H_invDmH(6,2)+X(2,6)*H_invDmH(6,1)+X&
&(1,5)*H_invDmH(5,2)+X(2,5)*H_invDmH(5,1)+X(1,4)*H_invDmH(4,2)+X(2&
&,4)*H_invDmH(4,1)+X(1,3)*H_invDmH(3,2)+X(2,3)*H_invDmH(3,1)+H_inv&
&DmH(2,1)*X(2,2)+X(1,2)*H_invDmH(2,2)+H_invDmH(1,1)*X(2,1)+X(1,1)*&
&H_invDmH(1,2))**2/2.0E+0+(X(1,8)*H_invDmH(8,1)+X(1,7)*H_invDmH(7,&
&1)+X(1,6)*H_invDmH(6,1)+X(1,5)*H_invDmH(5,1)+X(1,4)*H_invDmH(4,1)&
&+X(1,3)*H_invDmH(3,1)+X(1,2)*H_invDmH(2,1)+H_invDmH(1,1)*X(1,1)-1&
&)**2))
END 
SUBROUTINE vox_linear_at_quadr_mtr_jac(jac, mtr, X, H_invDmH, detDmH, gw) 
IMPLICIT NONE 
REAL(KIND=8) jac(1, 2) 
REAL(KIND=8) mtr(2, 1) 
REAL(KIND=8) X(3, 8) 
REAL(KIND=8) H_invDmH(8, 3) 
REAL(KIND=8) detDmH(1, 1) 
REAL(KIND=8) gw(1, 1) 
jac(1,1) = detDmH(1,1)*gw(1,1)*((X(3,8)*H_invDmH(8,3)+X(3,7)*H_in&
&vDmH(7,3)+X(3,6)*H_invDmH(6,3)+X(3,5)*H_invDmH(5,3)+X(3,4)*H_invD&
&mH(4,3)+H_invDmH(3,3)*X(3,3)+H_invDmH(2,3)*X(3,2)+H_invDmH(1,3)*X&
&(3,1)-1)**2+(X(2,8)*H_invDmH(8,3)+X(3,8)*H_invDmH(8,2)+X(2,7)*H_i&
&nvDmH(7,3)+X(3,7)*H_invDmH(7,2)+X(2,6)*H_invDmH(6,3)+X(3,6)*H_inv&
&DmH(6,2)+X(2,5)*H_invDmH(5,3)+X(3,5)*H_invDmH(5,2)+X(2,4)*H_invDm&
&H(4,3)+X(3,4)*H_invDmH(4,2)+H_invDmH(3,2)*X(3,3)+X(2,3)*H_invDmH(&
&3,3)+H_invDmH(2,2)*X(3,2)+H_invDmH(1,2)*X(3,1)+X(2,2)*H_invDmH(2,&
&3)+H_invDmH(1,3)*X(2,1))**2/2.0E+0+(X(1,8)*H_invDmH(8,3)+X(3,8)*H&
&_invDmH(8,1)+X(1,7)*H_invDmH(7,3)+X(3,7)*H_invDmH(7,1)+X(1,6)*H_i&
&nvDmH(6,3)+X(3,6)*H_invDmH(6,1)+X(1,5)*H_invDmH(5,3)+X(3,5)*H_inv&
&DmH(5,1)+X(1,4)*H_invDmH(4,3)+X(3,4)*H_invDmH(4,1)+H_invDmH(3,1)*&
&X(3,3)+X(1,3)*H_invDmH(3,3)+H_invDmH(2,1)*X(3,2)+H_invDmH(1,1)*X(&
&3,1)+X(1,2)*H_invDmH(2,3)+X(1,1)*H_invDmH(1,3))**2/2.0E+0+(X(2,8)&
&*H_invDmH(8,2)+X(2,7)*H_invDmH(7,2)+X(2,6)*H_invDmH(6,2)+X(2,5)*H&
&_invDmH(5,2)+X(2,4)*H_invDmH(4,2)+X(2,3)*H_invDmH(3,2)+H_invDmH(2&
&,2)*X(2,2)+H_invDmH(1,2)*X(2,1)-1)**2+(X(1,8)*H_invDmH(8,2)+X(2,8&
&)*H_invDmH(8,1)+X(1,7)*H_invDmH(7,2)+X(2,7)*H_invDmH(7,1)+X(1,6)*&
&H_invDmH(6,2)+X(2,6)*H_invDmH(6,1)+X(1,5)*H_invDmH(5,2)+X(2,5)*H_&
&invDmH(5,1)+X(1,4)*H_invDmH(4,2)+X(2,4)*H_invDmH(4,1)+X(1,3)*H_in&
&vDmH(3,2)+X(2,3)*H_invDmH(3,1)+H_invDmH(2,1)*X(2,2)+X(1,2)*H_invD&
&mH(2,2)+H_invDmH(1,1)*X(2,1)+X(1,1)*H_invDmH(1,2))**2/2.0E+0+(X(1&
&,8)*H_invDmH(8,1)+X(1,7)*H_invDmH(7,1)+X(1,6)*H_invDmH(6,1)+X(1,5&
&)*H_invDmH(5,1)+X(1,4)*H_invDmH(4,1)+X(1,3)*H_invDmH(3,1)+X(1,2)*&
&H_invDmH(2,1)+H_invDmH(1,1)*X(1,1)-1)**2)
jac(1,2) = (detDmH(1,1)*gw(1,1)*((2*X(3,8)*H_invDmH(8,3)+2*X(3,7)&
&*H_invDmH(7,3)+2*X(3,6)*H_invDmH(6,3)+2*X(3,5)*H_invDmH(5,3)+2*X(&
&3,4)*H_invDmH(4,3)+2*H_invDmH(3,3)*X(3,3)+2*H_invDmH(2,3)*X(3,2)+&
&2*H_invDmH(1,3)*X(3,1))/2.0E+0+(2*X(2,8)*H_invDmH(8,2)+2*X(2,7)*H&
&_invDmH(7,2)+2*X(2,6)*H_invDmH(6,2)+2*X(2,5)*H_invDmH(5,2)+2*X(2,&
&4)*H_invDmH(4,2)+2*X(2,3)*H_invDmH(3,2)+2*H_invDmH(2,2)*X(2,2)+2*&
&H_invDmH(1,2)*X(2,1))/2.0E+0+(2*X(1,8)*H_invDmH(8,1)+2*X(1,7)*H_i&
&nvDmH(7,1)+2*X(1,6)*H_invDmH(6,1)+2*X(1,5)*H_invDmH(5,1)+2*X(1,4)&
&*H_invDmH(4,1)+2*X(1,3)*H_invDmH(3,1)+2*X(1,2)*H_invDmH(2,1)+2*H_&
&invDmH(1,1)*X(1,1))/2.0E+0-3)**2)/2.0E+0
END 
SUBROUTINE vox_corotated_at_quadr_mtr(val, mtr, X, H_invDmH, R, detDmH, gw) 
IMPLICIT NONE 
REAL(KIND=8) val(1, 1) 
REAL(KIND=8) mtr(2, 1) 
REAL(KIND=8) X(3, 8) 
REAL(KIND=8) H_invDmH(8, 3) 
REAL(KIND=8) R(3, 3) 
REAL(KIND=8) detDmH(1, 1) 
REAL(KIND=8) gw(1, 1) 
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
REAL(KIND=8)  tt61 
REAL(KIND=8)  tt62 
REAL(KIND=8)  tt63 
REAL(KIND=8)  tt64 
REAL(KIND=8)  tt65 
REAL(KIND=8)  tt66 
REAL(KIND=8)  tt67 
REAL(KIND=8)  tt68 
REAL(KIND=8)  tt69 
REAL(KIND=8)  tt70 
REAL(KIND=8)  tt71 
REAL(KIND=8)  tt72 
tt1 = H_invDmH(1,1)*X(1,1)
tt2 = X(1,2)*H_invDmH(2,1)
tt3 = X(1,3)*H_invDmH(3,1)
tt4 = X(1,4)*H_invDmH(4,1)
tt5 = X(1,5)*H_invDmH(5,1)
tt6 = X(1,6)*H_invDmH(6,1)
tt7 = X(1,7)*H_invDmH(7,1)
tt8 = X(1,8)*H_invDmH(8,1)
tt9 = H_invDmH(1,1)*X(2,1)
tt10 = H_invDmH(2,1)*X(2,2)
tt11 = X(2,3)*H_invDmH(3,1)
tt12 = X(2,4)*H_invDmH(4,1)
tt13 = X(2,5)*H_invDmH(5,1)
tt14 = X(2,6)*H_invDmH(6,1)
tt15 = X(2,7)*H_invDmH(7,1)
tt16 = X(2,8)*H_invDmH(8,1)
tt17 = H_invDmH(1,1)*X(3,1)
tt18 = H_invDmH(2,1)*X(3,2)
tt19 = H_invDmH(3,1)*X(3,3)
tt20 = X(3,4)*H_invDmH(4,1)
tt21 = X(3,5)*H_invDmH(5,1)
tt22 = X(3,6)*H_invDmH(6,1)
tt23 = X(3,7)*H_invDmH(7,1)
tt24 = X(3,8)*H_invDmH(8,1)
tt25 = X(1,1)*H_invDmH(1,2)
tt26 = X(1,2)*H_invDmH(2,2)
tt27 = X(1,3)*H_invDmH(3,2)
tt28 = X(1,4)*H_invDmH(4,2)
tt29 = X(1,5)*H_invDmH(5,2)
tt30 = X(1,6)*H_invDmH(6,2)
tt31 = X(1,7)*H_invDmH(7,2)
tt32 = X(1,8)*H_invDmH(8,2)
tt33 = H_invDmH(1,2)*X(2,1)
tt34 = H_invDmH(2,2)*X(2,2)
tt35 = X(2,3)*H_invDmH(3,2)
tt36 = X(2,4)*H_invDmH(4,2)
tt37 = X(2,5)*H_invDmH(5,2)
tt38 = X(2,6)*H_invDmH(6,2)
tt39 = X(2,7)*H_invDmH(7,2)
tt40 = X(2,8)*H_invDmH(8,2)
tt41 = H_invDmH(1,2)*X(3,1)
tt42 = H_invDmH(2,2)*X(3,2)
tt43 = H_invDmH(3,2)*X(3,3)
tt44 = X(3,4)*H_invDmH(4,2)
tt45 = X(3,5)*H_invDmH(5,2)
tt46 = X(3,6)*H_invDmH(6,2)
tt47 = X(3,7)*H_invDmH(7,2)
tt48 = X(3,8)*H_invDmH(8,2)
tt49 = X(1,1)*H_invDmH(1,3)
tt50 = X(1,2)*H_invDmH(2,3)
tt51 = X(1,3)*H_invDmH(3,3)
tt52 = X(1,4)*H_invDmH(4,3)
tt53 = X(1,5)*H_invDmH(5,3)
tt54 = X(1,6)*H_invDmH(6,3)
tt55 = X(1,7)*H_invDmH(7,3)
tt56 = X(1,8)*H_invDmH(8,3)
tt57 = H_invDmH(1,3)*X(2,1)
tt58 = X(2,2)*H_invDmH(2,3)
tt59 = X(2,3)*H_invDmH(3,3)
tt60 = X(2,4)*H_invDmH(4,3)
tt61 = X(2,5)*H_invDmH(5,3)
tt62 = X(2,6)*H_invDmH(6,3)
tt63 = X(2,7)*H_invDmH(7,3)
tt64 = X(2,8)*H_invDmH(8,3)
tt65 = H_invDmH(1,3)*X(3,1)
tt66 = H_invDmH(2,3)*X(3,2)
tt67 = H_invDmH(3,3)*X(3,3)
tt68 = X(3,4)*H_invDmH(4,3)
tt69 = X(3,5)*H_invDmH(5,3)
tt70 = X(3,6)*H_invDmH(6,3)
tt71 = X(3,7)*H_invDmH(7,3)
tt72 = X(3,8)*H_invDmH(8,3)
val(1,1) = detDmH(1,1)*gw(1,1)*((mtr(2,1)*(R(3,3)*(tt72+tt71+tt70&
&+tt69+tt68+tt67+tt66+tt65)+R(2,3)*(tt64+tt63+tt62+tt61+tt60+tt59+&
&tt58+tt57)+R(1,3)*(tt56+tt55+tt54+tt53+tt52+tt51+tt50+tt49)+R(3,2&
&)*(tt48+tt47+tt46+tt45+tt44+tt43+tt42+tt41)+R(2,2)*(tt40+tt39+tt3&
&8+tt37+tt36+tt35+tt34+tt33)+R(1,2)*(tt32+tt31+tt30+tt29+tt28+tt27&
&+tt26+tt25)+R(3,1)*(tt24+tt23+tt22+tt21+tt20+tt19+tt18+tt17)+R(2,&
&1)*(tt16+tt15+tt14+tt13+tt12+tt11+tt10+tt9)+R(1,1)*(tt8+tt7+tt6+t&
&t5+tt4+tt3+tt2+tt1)-3)**2)/2.0E+0+mtr(1,1)*((tt72+tt71+tt70+tt69+&
&tt68+tt67-R(3,3)+tt66+tt65)**2+(tt64+tt63+tt62+tt61+tt60+tt59-R(2&
&,3)+tt58+tt57)**2+(tt56+tt55+tt54+tt53+tt52+tt51+tt50-R(1,3)+tt49&
&)**2+(tt48+tt47+tt46+tt45+tt44+tt43+tt42-R(3,2)+tt41)**2+(tt40+tt&
&39+tt38+tt37+tt36+tt35+tt34-R(2,2)+tt33)**2+(tt32+tt31+tt30+tt29+&
&tt28+tt27+tt26-R(1,2)+tt25)**2+(tt24+tt23+tt22+tt21+tt20+tt19+tt1&
&8+tt17-R(3,1))**2+(tt16+tt15+tt14+tt13+tt12+tt11+tt10+tt9-R(2,1))&
&**2+(tt8+tt7+tt6+tt5+tt4+tt3+tt2+tt1-R(1,1))**2))
END 
SUBROUTINE vox_corotated_at_quadr_mtr_jac(jac, mtr, X, H_invDmH, R, detDmH, gw)
                                                                               
IMPLICIT NONE 
REAL(KIND=8) jac(1, 2) 
REAL(KIND=8) mtr(2, 1) 
REAL(KIND=8) X(3, 8) 
REAL(KIND=8) H_invDmH(8, 3) 
REAL(KIND=8) R(3, 3) 
REAL(KIND=8) detDmH(1, 1) 
REAL(KIND=8) gw(1, 1) 
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
REAL(KIND=8)  tt61 
REAL(KIND=8)  tt62 
REAL(KIND=8)  tt63 
REAL(KIND=8)  tt64 
REAL(KIND=8)  tt65 
REAL(KIND=8)  tt66 
REAL(KIND=8)  tt67 
REAL(KIND=8)  tt68 
REAL(KIND=8)  tt69 
REAL(KIND=8)  tt70 
REAL(KIND=8)  tt71 
REAL(KIND=8)  tt72 
tt1 = H_invDmH(1,1)*X(1,1)
tt2 = X(1,2)*H_invDmH(2,1)
tt3 = X(1,3)*H_invDmH(3,1)
tt4 = X(1,4)*H_invDmH(4,1)
tt5 = X(1,5)*H_invDmH(5,1)
tt6 = X(1,6)*H_invDmH(6,1)
tt7 = X(1,7)*H_invDmH(7,1)
tt8 = X(1,8)*H_invDmH(8,1)
tt9 = H_invDmH(1,1)*X(2,1)
tt10 = H_invDmH(2,1)*X(2,2)
tt11 = X(2,3)*H_invDmH(3,1)
tt12 = X(2,4)*H_invDmH(4,1)
tt13 = X(2,5)*H_invDmH(5,1)
tt14 = X(2,6)*H_invDmH(6,1)
tt15 = X(2,7)*H_invDmH(7,1)
tt16 = X(2,8)*H_invDmH(8,1)
tt17 = H_invDmH(1,1)*X(3,1)
tt18 = H_invDmH(2,1)*X(3,2)
tt19 = H_invDmH(3,1)*X(3,3)
tt20 = X(3,4)*H_invDmH(4,1)
tt21 = X(3,5)*H_invDmH(5,1)
tt22 = X(3,6)*H_invDmH(6,1)
tt23 = X(3,7)*H_invDmH(7,1)
tt24 = X(3,8)*H_invDmH(8,1)
tt25 = X(1,1)*H_invDmH(1,2)
tt26 = X(1,2)*H_invDmH(2,2)
tt27 = X(1,3)*H_invDmH(3,2)
tt28 = X(1,4)*H_invDmH(4,2)
tt29 = X(1,5)*H_invDmH(5,2)
tt30 = X(1,6)*H_invDmH(6,2)
tt31 = X(1,7)*H_invDmH(7,2)
tt32 = X(1,8)*H_invDmH(8,2)
tt33 = H_invDmH(1,2)*X(2,1)
tt34 = H_invDmH(2,2)*X(2,2)
tt35 = X(2,3)*H_invDmH(3,2)
tt36 = X(2,4)*H_invDmH(4,2)
tt37 = X(2,5)*H_invDmH(5,2)
tt38 = X(2,6)*H_invDmH(6,2)
tt39 = X(2,7)*H_invDmH(7,2)
tt40 = X(2,8)*H_invDmH(8,2)
tt41 = H_invDmH(1,2)*X(3,1)
tt42 = H_invDmH(2,2)*X(3,2)
tt43 = H_invDmH(3,2)*X(3,3)
tt44 = X(3,4)*H_invDmH(4,2)
tt45 = X(3,5)*H_invDmH(5,2)
tt46 = X(3,6)*H_invDmH(6,2)
tt47 = X(3,7)*H_invDmH(7,2)
tt48 = X(3,8)*H_invDmH(8,2)
tt49 = X(1,1)*H_invDmH(1,3)
tt50 = X(1,2)*H_invDmH(2,3)
tt51 = X(1,3)*H_invDmH(3,3)
tt52 = X(1,4)*H_invDmH(4,3)
tt53 = X(1,5)*H_invDmH(5,3)
tt54 = X(1,6)*H_invDmH(6,3)
tt55 = X(1,7)*H_invDmH(7,3)
tt56 = X(1,8)*H_invDmH(8,3)
tt57 = H_invDmH(1,3)*X(2,1)
tt58 = X(2,2)*H_invDmH(2,3)
tt59 = X(2,3)*H_invDmH(3,3)
tt60 = X(2,4)*H_invDmH(4,3)
tt61 = X(2,5)*H_invDmH(5,3)
tt62 = X(2,6)*H_invDmH(6,3)
tt63 = X(2,7)*H_invDmH(7,3)
tt64 = X(2,8)*H_invDmH(8,3)
tt65 = H_invDmH(1,3)*X(3,1)
tt66 = H_invDmH(2,3)*X(3,2)
tt67 = H_invDmH(3,3)*X(3,3)
tt68 = X(3,4)*H_invDmH(4,3)
tt69 = X(3,5)*H_invDmH(5,3)
tt70 = X(3,6)*H_invDmH(6,3)
tt71 = X(3,7)*H_invDmH(7,3)
tt72 = X(3,8)*H_invDmH(8,3)
jac(1,1) = detDmH(1,1)*gw(1,1)*((tt72+tt71+tt70+tt69+tt68+tt67-R(&
&3,3)+tt66+tt65)**2+(tt64+tt63+tt62+tt61+tt60+tt59-R(2,3)+tt58+tt5&
&7)**2+(tt56+tt55+tt54+tt53+tt52+tt51+tt50-R(1,3)+tt49)**2+(tt48+t&
&t47+tt46+tt45+tt44+tt43+tt42-R(3,2)+tt41)**2+(tt40+tt39+tt38+tt37&
&+tt36+tt35+tt34-R(2,2)+tt33)**2+(tt32+tt31+tt30+tt29+tt28+tt27+tt&
&26-R(1,2)+tt25)**2+(tt24+tt23+tt22+tt21+tt20+tt19+tt18+tt17-R(3,1&
&))**2+(tt16+tt15+tt14+tt13+tt12+tt11+tt10+tt9-R(2,1))**2+(tt8+tt7&
&+tt6+tt5+tt4+tt3+tt2+tt1-R(1,1))**2)
jac(1,2) = (detDmH(1,1)*gw(1,1)*(R(3,3)*(tt72+tt71+tt70+tt69+tt68&
&+tt67+tt66+tt65)+R(2,3)*(tt64+tt63+tt62+tt61+tt60+tt59+tt58+tt57)&
&+R(1,3)*(tt56+tt55+tt54+tt53+tt52+tt51+tt50+tt49)+R(3,2)*(tt48+tt&
&47+tt46+tt45+tt44+tt43+tt42+tt41)+R(2,2)*(tt40+tt39+tt38+tt37+tt3&
&6+tt35+tt34+tt33)+R(1,2)*(tt32+tt31+tt30+tt29+tt28+tt27+tt26+tt25&
&)+R(3,1)*(tt24+tt23+tt22+tt21+tt20+tt19+tt18+tt17)+R(2,1)*(tt16+t&
&t15+tt14+tt13+tt12+tt11+tt10+tt9)+R(1,1)*(tt8+tt7+tt6+tt5+tt4+tt3&
&+tt2+tt1)-3)**2)/2.0E+0
END 
