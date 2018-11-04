SUBROUTINE vox_stvk_at_quadr(val, X, H_invDmH, detDmH, gw, mu, lam) 
IMPLICIT NONE 
REAL(KIND=8) val(1, 1) 
REAL(KIND=8) X(3, 8) 
REAL(KIND=8) H_invDmH(8, 3) 
REAL(KIND=8) detDmH(1, 1) 
REAL(KIND=8) gw(1, 1) 
REAL(KIND=8) mu(1, 1) 
REAL(KIND=8) lam(1, 1) 
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
val(1,1) = detDmH(1,1)*gw(1,1)*((lam(1,1)*(tt27/2.0E+0+tt26/2.0E+&
&0+tt1/2.0E+0)**2)/2.0E+0+mu(1,1)*(tt27**2/4.0E+0+(tt25*H_invDmH(8&
&,2)*H_invDmH(8,3)+tt24*H_invDmH(8,2)*H_invDmH(8,3)+tt23*H_invDmH(&
&8,2)*H_invDmH(8,3)+X(3,7)*X(3,8)*H_invDmH(7,2)*H_invDmH(8,3)+X(2,&
&7)*X(2,8)*H_invDmH(7,2)*H_invDmH(8,3)+X(1,7)*X(1,8)*H_invDmH(7,2)&
&*H_invDmH(8,3)+X(3,6)*X(3,8)*H_invDmH(6,2)*H_invDmH(8,3)+X(2,6)*X&
&(2,8)*H_invDmH(6,2)*H_invDmH(8,3)+X(1,6)*X(1,8)*H_invDmH(6,2)*H_i&
&nvDmH(8,3)+X(3,5)*X(3,8)*H_invDmH(5,2)*H_invDmH(8,3)+X(2,5)*X(2,8&
&)*H_invDmH(5,2)*H_invDmH(8,3)+X(1,5)*X(1,8)*H_invDmH(5,2)*H_invDm&
&H(8,3)+X(3,4)*X(3,8)*H_invDmH(4,2)*H_invDmH(8,3)+X(2,4)*X(2,8)*H_&
&invDmH(4,2)*H_invDmH(8,3)+X(1,4)*X(1,8)*H_invDmH(4,2)*H_invDmH(8,&
&3)+H_invDmH(3,2)*X(3,3)*X(3,8)*H_invDmH(8,3)+H_invDmH(2,2)*X(3,2)&
&*X(3,8)*H_invDmH(8,3)+H_invDmH(1,2)*X(3,1)*X(3,8)*H_invDmH(8,3)+X&
&(2,3)*X(2,8)*H_invDmH(3,2)*H_invDmH(8,3)+X(1,3)*X(1,8)*H_invDmH(3&
&,2)*H_invDmH(8,3)+H_invDmH(2,2)*X(2,2)*X(2,8)*H_invDmH(8,3)+H_inv&
&DmH(1,2)*X(2,1)*X(2,8)*H_invDmH(8,3)+X(1,2)*X(1,8)*H_invDmH(2,2)*&
&H_invDmH(8,3)+X(1,1)*H_invDmH(1,2)*X(1,8)*H_invDmH(8,3)+X(3,7)*X(&
&3,8)*H_invDmH(7,3)*H_invDmH(8,2)+X(2,7)*X(2,8)*H_invDmH(7,3)*H_in&
&vDmH(8,2)+X(1,7)*X(1,8)*H_invDmH(7,3)*H_invDmH(8,2)+X(3,6)*X(3,8)&
&*H_invDmH(6,3)*H_invDmH(8,2)+X(2,6)*X(2,8)*H_invDmH(6,3)*H_invDmH&
&(8,2)+X(1,6)*X(1,8)*H_invDmH(6,3)*H_invDmH(8,2)+X(3,5)*X(3,8)*H_i&
&nvDmH(5,3)*H_invDmH(8,2)+X(2,5)*X(2,8)*H_invDmH(5,3)*H_invDmH(8,2&
&)+X(1,5)*X(1,8)*H_invDmH(5,3)*H_invDmH(8,2)+X(3,4)*X(3,8)*H_invDm&
&H(4,3)*H_invDmH(8,2)+X(2,4)*X(2,8)*H_invDmH(4,3)*H_invDmH(8,2)+X(&
&1,4)*X(1,8)*H_invDmH(4,3)*H_invDmH(8,2)+H_invDmH(3,3)*X(3,3)*X(3,&
&8)*H_invDmH(8,2)+H_invDmH(2,3)*X(3,2)*X(3,8)*H_invDmH(8,2)+H_invD&
&mH(1,3)*X(3,1)*X(3,8)*H_invDmH(8,2)+X(2,3)*X(2,8)*H_invDmH(3,3)*H&
&_invDmH(8,2)+X(1,3)*X(1,8)*H_invDmH(3,3)*H_invDmH(8,2)+X(2,2)*H_i&
&nvDmH(2,3)*X(2,8)*H_invDmH(8,2)+H_invDmH(1,3)*X(2,1)*X(2,8)*H_inv&
&DmH(8,2)+X(1,2)*X(1,8)*H_invDmH(2,3)*H_invDmH(8,2)+X(1,1)*H_invDm&
&H(1,3)*X(1,8)*H_invDmH(8,2)+tt22*H_invDmH(7,2)*H_invDmH(7,3)+tt21&
&*H_invDmH(7,2)*H_invDmH(7,3)+tt20*H_invDmH(7,2)*H_invDmH(7,3)+X(3&
&,6)*X(3,7)*H_invDmH(6,2)*H_invDmH(7,3)+X(2,6)*X(2,7)*H_invDmH(6,2&
&)*H_invDmH(7,3)+X(1,6)*X(1,7)*H_invDmH(6,2)*H_invDmH(7,3)+X(3,5)*&
&X(3,7)*H_invDmH(5,2)*H_invDmH(7,3)+X(2,5)*X(2,7)*H_invDmH(5,2)*H_&
&invDmH(7,3)+X(1,5)*X(1,7)*H_invDmH(5,2)*H_invDmH(7,3)+X(3,4)*X(3,&
&7)*H_invDmH(4,2)*H_invDmH(7,3)+X(2,4)*X(2,7)*H_invDmH(4,2)*H_invD&
&mH(7,3)+X(1,4)*X(1,7)*H_invDmH(4,2)*H_invDmH(7,3)+H_invDmH(3,2)*X&
&(3,3)*X(3,7)*H_invDmH(7,3)+H_invDmH(2,2)*X(3,2)*X(3,7)*H_invDmH(7&
&,3)+H_invDmH(1,2)*X(3,1)*X(3,7)*H_invDmH(7,3)+X(2,3)*X(2,7)*H_inv&
&DmH(3,2)*H_invDmH(7,3)+X(1,3)*X(1,7)*H_invDmH(3,2)*H_invDmH(7,3)+&
&H_invDmH(2,2)*X(2,2)*X(2,7)*H_invDmH(7,3)+H_invDmH(1,2)*X(2,1)*X(&
&2,7)*H_invDmH(7,3)+X(1,2)*X(1,7)*H_invDmH(2,2)*H_invDmH(7,3)+X(1,&
&1)*H_invDmH(1,2)*X(1,7)*H_invDmH(7,3)+X(3,6)*X(3,7)*H_invDmH(6,3)&
&*H_invDmH(7,2)+X(2,6)*X(2,7)*H_invDmH(6,3)*H_invDmH(7,2)+X(1,6)*X&
&(1,7)*H_invDmH(6,3)*H_invDmH(7,2)+X(3,5)*X(3,7)*H_invDmH(5,3)*H_i&
&nvDmH(7,2)+X(2,5)*X(2,7)*H_invDmH(5,3)*H_invDmH(7,2)+X(1,5)*X(1,7&
&)*H_invDmH(5,3)*H_invDmH(7,2)+X(3,4)*X(3,7)*H_invDmH(4,3)*H_invDm&
&H(7,2)+X(2,4)*X(2,7)*H_invDmH(4,3)*H_invDmH(7,2)+X(1,4)*X(1,7)*H_&
&invDmH(4,3)*H_invDmH(7,2)+H_invDmH(3,3)*X(3,3)*X(3,7)*H_invDmH(7,&
&2)+H_invDmH(2,3)*X(3,2)*X(3,7)*H_invDmH(7,2)+H_invDmH(1,3)*X(3,1)&
&*X(3,7)*H_invDmH(7,2)+X(2,3)*X(2,7)*H_invDmH(3,3)*H_invDmH(7,2)+X&
&(1,3)*X(1,7)*H_invDmH(3,3)*H_invDmH(7,2)+X(2,2)*H_invDmH(2,3)*X(2&
&,7)*H_invDmH(7,2)+H_invDmH(1,3)*X(2,1)*X(2,7)*H_invDmH(7,2)+X(1,2&
&)*X(1,7)*H_invDmH(2,3)*H_invDmH(7,2)+X(1,1)*H_invDmH(1,3)*X(1,7)*&
&H_invDmH(7,2)+tt19*H_invDmH(6,2)*H_invDmH(6,3)+tt18*H_invDmH(6,2)&
&*H_invDmH(6,3)+tt17*H_invDmH(6,2)*H_invDmH(6,3)+X(3,5)*X(3,6)*H_i&
&nvDmH(5,2)*H_invDmH(6,3)+X(2,5)*X(2,6)*H_invDmH(5,2)*H_invDmH(6,3&
&)+X(1,5)*X(1,6)*H_invDmH(5,2)*H_invDmH(6,3)+X(3,4)*X(3,6)*H_invDm&
&H(4,2)*H_invDmH(6,3)+X(2,4)*X(2,6)*H_invDmH(4,2)*H_invDmH(6,3)+X(&
&1,4)*X(1,6)*H_invDmH(4,2)*H_invDmH(6,3)+H_invDmH(3,2)*X(3,3)*X(3,&
&6)*H_invDmH(6,3)+H_invDmH(2,2)*X(3,2)*X(3,6)*H_invDmH(6,3)+H_invD&
&mH(1,2)*X(3,1)*X(3,6)*H_invDmH(6,3)+X(2,3)*X(2,6)*H_invDmH(3,2)*H&
&_invDmH(6,3)+X(1,3)*X(1,6)*H_invDmH(3,2)*H_invDmH(6,3)+H_invDmH(2&
&,2)*X(2,2)*X(2,6)*H_invDmH(6,3)+H_invDmH(1,2)*X(2,1)*X(2,6)*H_inv&
&DmH(6,3)+X(1,2)*X(1,6)*H_invDmH(2,2)*H_invDmH(6,3)+X(1,1)*H_invDm&
&H(1,2)*X(1,6)*H_invDmH(6,3)+X(3,5)*X(3,6)*H_invDmH(5,3)*H_invDmH(&
&6,2)+X(2,5)*X(2,6)*H_invDmH(5,3)*H_invDmH(6,2)+X(1,5)*X(1,6)*H_in&
&vDmH(5,3)*H_invDmH(6,2)+X(3,4)*X(3,6)*H_invDmH(4,3)*H_invDmH(6,2)&
&+X(2,4)*X(2,6)*H_invDmH(4,3)*H_invDmH(6,2)+X(1,4)*X(1,6)*H_invDmH&
&(4,3)*H_invDmH(6,2)+H_invDmH(3,3)*X(3,3)*X(3,6)*H_invDmH(6,2)+H_i&
&nvDmH(2,3)*X(3,2)*X(3,6)*H_invDmH(6,2)+H_invDmH(1,3)*X(3,1)*X(3,6&
&)*H_invDmH(6,2)+X(2,3)*X(2,6)*H_invDmH(3,3)*H_invDmH(6,2)+X(1,3)*&
&X(1,6)*H_invDmH(3,3)*H_invDmH(6,2)+X(2,2)*H_invDmH(2,3)*X(2,6)*H_&
&invDmH(6,2)+H_invDmH(1,3)*X(2,1)*X(2,6)*H_invDmH(6,2)+X(1,2)*X(1,&
&6)*H_invDmH(2,3)*H_invDmH(6,2)+X(1,1)*H_invDmH(1,3)*X(1,6)*H_invD&
&mH(6,2)+tt16*H_invDmH(5,2)*H_invDmH(5,3)+tt15*H_invDmH(5,2)*H_inv&
&DmH(5,3)+tt14*H_invDmH(5,2)*H_invDmH(5,3)+X(3,4)*X(3,5)*H_invDmH(&
&4,2)*H_invDmH(5,3)+X(2,4)*X(2,5)*H_invDmH(4,2)*H_invDmH(5,3)+X(1,&
&4)*X(1,5)*H_invDmH(4,2)*H_invDmH(5,3)+H_invDmH(3,2)*X(3,3)*X(3,5)&
&*H_invDmH(5,3)+H_invDmH(2,2)*X(3,2)*X(3,5)*H_invDmH(5,3)+H_invDmH&
&(1,2)*X(3,1)*X(3,5)*H_invDmH(5,3)+X(2,3)*X(2,5)*H_invDmH(3,2)*H_i&
&nvDmH(5,3)+X(1,3)*X(1,5)*H_invDmH(3,2)*H_invDmH(5,3)+H_invDmH(2,2&
&)*X(2,2)*X(2,5)*H_invDmH(5,3)+H_invDmH(1,2)*X(2,1)*X(2,5)*H_invDm&
&H(5,3)+X(1,2)*X(1,5)*H_invDmH(2,2)*H_invDmH(5,3)+X(1,1)*H_invDmH(&
&1,2)*X(1,5)*H_invDmH(5,3)+X(3,4)*X(3,5)*H_invDmH(4,3)*H_invDmH(5,&
&2)+X(2,4)*X(2,5)*H_invDmH(4,3)*H_invDmH(5,2)+X(1,4)*X(1,5)*H_invD&
&mH(4,3)*H_invDmH(5,2)+H_invDmH(3,3)*X(3,3)*X(3,5)*H_invDmH(5,2)+H&
&_invDmH(2,3)*X(3,2)*X(3,5)*H_invDmH(5,2)+H_invDmH(1,3)*X(3,1)*X(3&
&,5)*H_invDmH(5,2)+X(2,3)*X(2,5)*H_invDmH(3,3)*H_invDmH(5,2)+X(1,3&
&)*X(1,5)*H_invDmH(3,3)*H_invDmH(5,2)+X(2,2)*H_invDmH(2,3)*X(2,5)*&
&H_invDmH(5,2)+H_invDmH(1,3)*X(2,1)*X(2,5)*H_invDmH(5,2)+X(1,2)*X(&
&1,5)*H_invDmH(2,3)*H_invDmH(5,2)+X(1,1)*H_invDmH(1,3)*X(1,5)*H_in&
&vDmH(5,2)+tt13*H_invDmH(4,2)*H_invDmH(4,3)+tt12*H_invDmH(4,2)*H_i&
&nvDmH(4,3)+tt11*H_invDmH(4,2)*H_invDmH(4,3)+H_invDmH(3,2)*X(3,3)*&
&X(3,4)*H_invDmH(4,3)+H_invDmH(2,2)*X(3,2)*X(3,4)*H_invDmH(4,3)+H_&
&invDmH(1,2)*X(3,1)*X(3,4)*H_invDmH(4,3)+X(2,3)*X(2,4)*H_invDmH(3,&
&2)*H_invDmH(4,3)+X(1,3)*X(1,4)*H_invDmH(3,2)*H_invDmH(4,3)+H_invD&
&mH(2,2)*X(2,2)*X(2,4)*H_invDmH(4,3)+H_invDmH(1,2)*X(2,1)*X(2,4)*H&
&_invDmH(4,3)+X(1,2)*X(1,4)*H_invDmH(2,2)*H_invDmH(4,3)+X(1,1)*H_i&
&nvDmH(1,2)*X(1,4)*H_invDmH(4,3)+H_invDmH(3,3)*X(3,3)*X(3,4)*H_inv&
&DmH(4,2)+H_invDmH(2,3)*X(3,2)*X(3,4)*H_invDmH(4,2)+H_invDmH(1,3)*&
&X(3,1)*X(3,4)*H_invDmH(4,2)+X(2,3)*X(2,4)*H_invDmH(3,3)*H_invDmH(&
&4,2)+X(1,3)*X(1,4)*H_invDmH(3,3)*H_invDmH(4,2)+X(2,2)*H_invDmH(2,&
&3)*X(2,4)*H_invDmH(4,2)+H_invDmH(1,3)*X(2,1)*X(2,4)*H_invDmH(4,2)&
&+X(1,2)*X(1,4)*H_invDmH(2,3)*H_invDmH(4,2)+X(1,1)*H_invDmH(1,3)*X&
&(1,4)*H_invDmH(4,2)+H_invDmH(3,2)*H_invDmH(3,3)*tt10+H_invDmH(2,2&
&)*X(3,2)*H_invDmH(3,3)*X(3,3)+H_invDmH(1,2)*X(3,1)*H_invDmH(3,3)*&
&X(3,3)+H_invDmH(2,3)*H_invDmH(3,2)*X(3,2)*X(3,3)+H_invDmH(1,3)*X(&
&3,1)*H_invDmH(3,2)*X(3,3)+tt8*H_invDmH(3,2)*H_invDmH(3,3)+tt7*H_i&
&nvDmH(3,2)*H_invDmH(3,3)+H_invDmH(2,2)*X(2,2)*X(2,3)*H_invDmH(3,3&
&)+H_invDmH(1,2)*X(2,1)*X(2,3)*H_invDmH(3,3)+X(1,2)*X(1,3)*H_invDm&
&H(2,2)*H_invDmH(3,3)+X(1,1)*H_invDmH(1,2)*X(1,3)*H_invDmH(3,3)+H_&
&invDmH(2,2)*H_invDmH(2,3)*tt9+H_invDmH(1,2)*H_invDmH(2,3)*X(3,1)*&
&X(3,2)+H_invDmH(1,3)*H_invDmH(2,2)*X(3,1)*X(3,2)+X(2,2)*H_invDmH(&
&2,3)*X(2,3)*H_invDmH(3,2)+H_invDmH(1,3)*X(2,1)*X(2,3)*H_invDmH(3,&
&2)+X(1,2)*X(1,3)*H_invDmH(2,3)*H_invDmH(3,2)+X(1,1)*H_invDmH(1,3)&
&*X(1,3)*H_invDmH(3,2)+H_invDmH(1,2)*H_invDmH(1,3)*tt6+H_invDmH(2,&
&2)*tt5*H_invDmH(2,3)+H_invDmH(1,2)*X(2,1)*X(2,2)*H_invDmH(2,3)+tt&
&4*H_invDmH(2,2)*H_invDmH(2,3)+X(1,1)*H_invDmH(1,2)*X(1,2)*H_invDm&
&H(2,3)+H_invDmH(1,3)*X(2,1)*H_invDmH(2,2)*X(2,2)+X(1,1)*X(1,2)*H_&
&invDmH(1,3)*H_invDmH(2,2)+H_invDmH(1,2)*H_invDmH(1,3)*tt3+tt2*H_i&
&nvDmH(1,2)*H_invDmH(1,3))**2/2.0E+0+(tt25*H_invDmH(8,1)*H_invDmH(&
&8,3)+tt24*H_invDmH(8,1)*H_invDmH(8,3)+tt23*H_invDmH(8,1)*H_invDmH&
&(8,3)+X(3,7)*X(3,8)*H_invDmH(7,1)*H_invDmH(8,3)+X(2,7)*X(2,8)*H_i&
&nvDmH(7,1)*H_invDmH(8,3)+X(1,7)*X(1,8)*H_invDmH(7,1)*H_invDmH(8,3&
&)+X(3,6)*X(3,8)*H_invDmH(6,1)*H_invDmH(8,3)+X(2,6)*X(2,8)*H_invDm&
&H(6,1)*H_invDmH(8,3)+X(1,6)*X(1,8)*H_invDmH(6,1)*H_invDmH(8,3)+X(&
&3,5)*X(3,8)*H_invDmH(5,1)*H_invDmH(8,3)+X(2,5)*X(2,8)*H_invDmH(5,&
&1)*H_invDmH(8,3)+X(1,5)*X(1,8)*H_invDmH(5,1)*H_invDmH(8,3)+X(3,4)&
&*X(3,8)*H_invDmH(4,1)*H_invDmH(8,3)+X(2,4)*X(2,8)*H_invDmH(4,1)*H&
&_invDmH(8,3)+X(1,4)*X(1,8)*H_invDmH(4,1)*H_invDmH(8,3)+H_invDmH(3&
&,1)*X(3,3)*X(3,8)*H_invDmH(8,3)+H_invDmH(2,1)*X(3,2)*X(3,8)*H_inv&
&DmH(8,3)+H_invDmH(1,1)*X(3,1)*X(3,8)*H_invDmH(8,3)+X(2,3)*X(2,8)*&
&H_invDmH(3,1)*H_invDmH(8,3)+X(1,3)*X(1,8)*H_invDmH(3,1)*H_invDmH(&
&8,3)+H_invDmH(2,1)*X(2,2)*X(2,8)*H_invDmH(8,3)+H_invDmH(1,1)*X(2,&
&1)*X(2,8)*H_invDmH(8,3)+X(1,2)*X(1,8)*H_invDmH(2,1)*H_invDmH(8,3)&
&+H_invDmH(1,1)*X(1,1)*X(1,8)*H_invDmH(8,3)+X(3,7)*X(3,8)*H_invDmH&
&(7,3)*H_invDmH(8,1)+X(2,7)*X(2,8)*H_invDmH(7,3)*H_invDmH(8,1)+X(1&
&,7)*X(1,8)*H_invDmH(7,3)*H_invDmH(8,1)+X(3,6)*X(3,8)*H_invDmH(6,3&
&)*H_invDmH(8,1)+X(2,6)*X(2,8)*H_invDmH(6,3)*H_invDmH(8,1)+X(1,6)*&
&X(1,8)*H_invDmH(6,3)*H_invDmH(8,1)+X(3,5)*X(3,8)*H_invDmH(5,3)*H_&
&invDmH(8,1)+X(2,5)*X(2,8)*H_invDmH(5,3)*H_invDmH(8,1)+X(1,5)*X(1,&
&8)*H_invDmH(5,3)*H_invDmH(8,1)+X(3,4)*X(3,8)*H_invDmH(4,3)*H_invD&
&mH(8,1)+X(2,4)*X(2,8)*H_invDmH(4,3)*H_invDmH(8,1)+X(1,4)*X(1,8)*H&
&_invDmH(4,3)*H_invDmH(8,1)+H_invDmH(3,3)*X(3,3)*X(3,8)*H_invDmH(8&
&,1)+H_invDmH(2,3)*X(3,2)*X(3,8)*H_invDmH(8,1)+H_invDmH(1,3)*X(3,1&
&)*X(3,8)*H_invDmH(8,1)+X(2,3)*X(2,8)*H_invDmH(3,3)*H_invDmH(8,1)+&
&X(1,3)*X(1,8)*H_invDmH(3,3)*H_invDmH(8,1)+X(2,2)*H_invDmH(2,3)*X(&
&2,8)*H_invDmH(8,1)+H_invDmH(1,3)*X(2,1)*X(2,8)*H_invDmH(8,1)+X(1,&
&2)*X(1,8)*H_invDmH(2,3)*H_invDmH(8,1)+X(1,1)*H_invDmH(1,3)*X(1,8)&
&*H_invDmH(8,1)+tt22*H_invDmH(7,1)*H_invDmH(7,3)+tt21*H_invDmH(7,1&
&)*H_invDmH(7,3)+tt20*H_invDmH(7,1)*H_invDmH(7,3)+X(3,6)*X(3,7)*H_&
&invDmH(6,1)*H_invDmH(7,3)+X(2,6)*X(2,7)*H_invDmH(6,1)*H_invDmH(7,&
&3)+X(1,6)*X(1,7)*H_invDmH(6,1)*H_invDmH(7,3)+X(3,5)*X(3,7)*H_invD&
&mH(5,1)*H_invDmH(7,3)+X(2,5)*X(2,7)*H_invDmH(5,1)*H_invDmH(7,3)+X&
&(1,5)*X(1,7)*H_invDmH(5,1)*H_invDmH(7,3)+X(3,4)*X(3,7)*H_invDmH(4&
&,1)*H_invDmH(7,3)+X(2,4)*X(2,7)*H_invDmH(4,1)*H_invDmH(7,3)+X(1,4&
&)*X(1,7)*H_invDmH(4,1)*H_invDmH(7,3)+H_invDmH(3,1)*X(3,3)*X(3,7)*&
&H_invDmH(7,3)+H_invDmH(2,1)*X(3,2)*X(3,7)*H_invDmH(7,3)+H_invDmH(&
&1,1)*X(3,1)*X(3,7)*H_invDmH(7,3)+X(2,3)*X(2,7)*H_invDmH(3,1)*H_in&
&vDmH(7,3)+X(1,3)*X(1,7)*H_invDmH(3,1)*H_invDmH(7,3)+H_invDmH(2,1)&
&*X(2,2)*X(2,7)*H_invDmH(7,3)+H_invDmH(1,1)*X(2,1)*X(2,7)*H_invDmH&
&(7,3)+X(1,2)*X(1,7)*H_invDmH(2,1)*H_invDmH(7,3)+H_invDmH(1,1)*X(1&
&,1)*X(1,7)*H_invDmH(7,3)+X(3,6)*X(3,7)*H_invDmH(6,3)*H_invDmH(7,1&
&)+X(2,6)*X(2,7)*H_invDmH(6,3)*H_invDmH(7,1)+X(1,6)*X(1,7)*H_invDm&
&H(6,3)*H_invDmH(7,1)+X(3,5)*X(3,7)*H_invDmH(5,3)*H_invDmH(7,1)+X(&
&2,5)*X(2,7)*H_invDmH(5,3)*H_invDmH(7,1)+X(1,5)*X(1,7)*H_invDmH(5,&
&3)*H_invDmH(7,1)+X(3,4)*X(3,7)*H_invDmH(4,3)*H_invDmH(7,1)+X(2,4)&
&*X(2,7)*H_invDmH(4,3)*H_invDmH(7,1)+X(1,4)*X(1,7)*H_invDmH(4,3)*H&
&_invDmH(7,1)+H_invDmH(3,3)*X(3,3)*X(3,7)*H_invDmH(7,1)+H_invDmH(2&
&,3)*X(3,2)*X(3,7)*H_invDmH(7,1)+H_invDmH(1,3)*X(3,1)*X(3,7)*H_inv&
&DmH(7,1)+X(2,3)*X(2,7)*H_invDmH(3,3)*H_invDmH(7,1)+X(1,3)*X(1,7)*&
&H_invDmH(3,3)*H_invDmH(7,1)+X(2,2)*H_invDmH(2,3)*X(2,7)*H_invDmH(&
&7,1)+H_invDmH(1,3)*X(2,1)*X(2,7)*H_invDmH(7,1)+X(1,2)*X(1,7)*H_in&
&vDmH(2,3)*H_invDmH(7,1)+X(1,1)*H_invDmH(1,3)*X(1,7)*H_invDmH(7,1)&
&+tt19*H_invDmH(6,1)*H_invDmH(6,3)+tt18*H_invDmH(6,1)*H_invDmH(6,3&
&)+tt17*H_invDmH(6,1)*H_invDmH(6,3)+X(3,5)*X(3,6)*H_invDmH(5,1)*H_&
&invDmH(6,3)+X(2,5)*X(2,6)*H_invDmH(5,1)*H_invDmH(6,3)+X(1,5)*X(1,&
&6)*H_invDmH(5,1)*H_invDmH(6,3)+X(3,4)*X(3,6)*H_invDmH(4,1)*H_invD&
&mH(6,3)+X(2,4)*X(2,6)*H_invDmH(4,1)*H_invDmH(6,3)+X(1,4)*X(1,6)*H&
&_invDmH(4,1)*H_invDmH(6,3)+H_invDmH(3,1)*X(3,3)*X(3,6)*H_invDmH(6&
&,3)+H_invDmH(2,1)*X(3,2)*X(3,6)*H_invDmH(6,3)+H_invDmH(1,1)*X(3,1&
&)*X(3,6)*H_invDmH(6,3)+X(2,3)*X(2,6)*H_invDmH(3,1)*H_invDmH(6,3)+&
&X(1,3)*X(1,6)*H_invDmH(3,1)*H_invDmH(6,3)+H_invDmH(2,1)*X(2,2)*X(&
&2,6)*H_invDmH(6,3)+H_invDmH(1,1)*X(2,1)*X(2,6)*H_invDmH(6,3)+X(1,&
&2)*X(1,6)*H_invDmH(2,1)*H_invDmH(6,3)+H_invDmH(1,1)*X(1,1)*X(1,6)&
&*H_invDmH(6,3)+X(3,5)*X(3,6)*H_invDmH(5,3)*H_invDmH(6,1)+X(2,5)*X&
&(2,6)*H_invDmH(5,3)*H_invDmH(6,1)+X(1,5)*X(1,6)*H_invDmH(5,3)*H_i&
&nvDmH(6,1)+X(3,4)*X(3,6)*H_invDmH(4,3)*H_invDmH(6,1)+X(2,4)*X(2,6&
&)*H_invDmH(4,3)*H_invDmH(6,1)+X(1,4)*X(1,6)*H_invDmH(4,3)*H_invDm&
&H(6,1)+H_invDmH(3,3)*X(3,3)*X(3,6)*H_invDmH(6,1)+H_invDmH(2,3)*X(&
&3,2)*X(3,6)*H_invDmH(6,1)+H_invDmH(1,3)*X(3,1)*X(3,6)*H_invDmH(6,&
&1)+X(2,3)*X(2,6)*H_invDmH(3,3)*H_invDmH(6,1)+X(1,3)*X(1,6)*H_invD&
&mH(3,3)*H_invDmH(6,1)+X(2,2)*H_invDmH(2,3)*X(2,6)*H_invDmH(6,1)+H&
&_invDmH(1,3)*X(2,1)*X(2,6)*H_invDmH(6,1)+X(1,2)*X(1,6)*H_invDmH(2&
&,3)*H_invDmH(6,1)+X(1,1)*H_invDmH(1,3)*X(1,6)*H_invDmH(6,1)+tt16*&
&H_invDmH(5,1)*H_invDmH(5,3)+tt15*H_invDmH(5,1)*H_invDmH(5,3)+tt14&
&*H_invDmH(5,1)*H_invDmH(5,3)+X(3,4)*X(3,5)*H_invDmH(4,1)*H_invDmH&
&(5,3)+X(2,4)*X(2,5)*H_invDmH(4,1)*H_invDmH(5,3)+X(1,4)*X(1,5)*H_i&
&nvDmH(4,1)*H_invDmH(5,3)+H_invDmH(3,1)*X(3,3)*X(3,5)*H_invDmH(5,3&
&)+H_invDmH(2,1)*X(3,2)*X(3,5)*H_invDmH(5,3)+H_invDmH(1,1)*X(3,1)*&
&X(3,5)*H_invDmH(5,3)+X(2,3)*X(2,5)*H_invDmH(3,1)*H_invDmH(5,3)+X(&
&1,3)*X(1,5)*H_invDmH(3,1)*H_invDmH(5,3)+H_invDmH(2,1)*X(2,2)*X(2,&
&5)*H_invDmH(5,3)+H_invDmH(1,1)*X(2,1)*X(2,5)*H_invDmH(5,3)+X(1,2)&
&*X(1,5)*H_invDmH(2,1)*H_invDmH(5,3)+H_invDmH(1,1)*X(1,1)*X(1,5)*H&
&_invDmH(5,3)+X(3,4)*X(3,5)*H_invDmH(4,3)*H_invDmH(5,1)+X(2,4)*X(2&
&,5)*H_invDmH(4,3)*H_invDmH(5,1)+X(1,4)*X(1,5)*H_invDmH(4,3)*H_inv&
&DmH(5,1)+H_invDmH(3,3)*X(3,3)*X(3,5)*H_invDmH(5,1)+H_invDmH(2,3)*&
&X(3,2)*X(3,5)*H_invDmH(5,1)+H_invDmH(1,3)*X(3,1)*X(3,5)*H_invDmH(&
&5,1)+X(2,3)*X(2,5)*H_invDmH(3,3)*H_invDmH(5,1)+X(1,3)*X(1,5)*H_in&
&vDmH(3,3)*H_invDmH(5,1)+X(2,2)*H_invDmH(2,3)*X(2,5)*H_invDmH(5,1)&
&+H_invDmH(1,3)*X(2,1)*X(2,5)*H_invDmH(5,1)+X(1,2)*X(1,5)*H_invDmH&
&(2,3)*H_invDmH(5,1)+X(1,1)*H_invDmH(1,3)*X(1,5)*H_invDmH(5,1)+tt1&
&3*H_invDmH(4,1)*H_invDmH(4,3)+tt12*H_invDmH(4,1)*H_invDmH(4,3)+tt&
&11*H_invDmH(4,1)*H_invDmH(4,3)+H_invDmH(3,1)*X(3,3)*X(3,4)*H_invD&
&mH(4,3)+H_invDmH(2,1)*X(3,2)*X(3,4)*H_invDmH(4,3)+H_invDmH(1,1)*X&
&(3,1)*X(3,4)*H_invDmH(4,3)+X(2,3)*X(2,4)*H_invDmH(3,1)*H_invDmH(4&
&,3)+X(1,3)*X(1,4)*H_invDmH(3,1)*H_invDmH(4,3)+H_invDmH(2,1)*X(2,2&
&)*X(2,4)*H_invDmH(4,3)+H_invDmH(1,1)*X(2,1)*X(2,4)*H_invDmH(4,3)+&
&X(1,2)*X(1,4)*H_invDmH(2,1)*H_invDmH(4,3)+H_invDmH(1,1)*X(1,1)*X(&
&1,4)*H_invDmH(4,3)+H_invDmH(3,3)*X(3,3)*X(3,4)*H_invDmH(4,1)+H_in&
&vDmH(2,3)*X(3,2)*X(3,4)*H_invDmH(4,1)+H_invDmH(1,3)*X(3,1)*X(3,4)&
&*H_invDmH(4,1)+X(2,3)*X(2,4)*H_invDmH(3,3)*H_invDmH(4,1)+X(1,3)*X&
&(1,4)*H_invDmH(3,3)*H_invDmH(4,1)+X(2,2)*H_invDmH(2,3)*X(2,4)*H_i&
&nvDmH(4,1)+H_invDmH(1,3)*X(2,1)*X(2,4)*H_invDmH(4,1)+X(1,2)*X(1,4&
&)*H_invDmH(2,3)*H_invDmH(4,1)+X(1,1)*H_invDmH(1,3)*X(1,4)*H_invDm&
&H(4,1)+H_invDmH(3,1)*H_invDmH(3,3)*tt10+H_invDmH(2,1)*X(3,2)*H_in&
&vDmH(3,3)*X(3,3)+H_invDmH(1,1)*X(3,1)*H_invDmH(3,3)*X(3,3)+H_invD&
&mH(2,3)*H_invDmH(3,1)*X(3,2)*X(3,3)+H_invDmH(1,3)*H_invDmH(3,1)*X&
&(3,1)*X(3,3)+tt8*H_invDmH(3,1)*H_invDmH(3,3)+tt7*H_invDmH(3,1)*H_&
&invDmH(3,3)+H_invDmH(2,1)*X(2,2)*X(2,3)*H_invDmH(3,3)+H_invDmH(1,&
&1)*X(2,1)*X(2,3)*H_invDmH(3,3)+X(1,2)*X(1,3)*H_invDmH(2,1)*H_invD&
&mH(3,3)+H_invDmH(1,1)*X(1,1)*X(1,3)*H_invDmH(3,3)+H_invDmH(2,1)*H&
&_invDmH(2,3)*tt9+H_invDmH(1,1)*H_invDmH(2,3)*X(3,1)*X(3,2)+H_invD&
&mH(1,3)*H_invDmH(2,1)*X(3,1)*X(3,2)+H_invDmH(1,1)*H_invDmH(1,3)*t&
&t6+X(2,2)*H_invDmH(2,3)*X(2,3)*H_invDmH(3,1)+H_invDmH(1,3)*X(2,1)&
&*X(2,3)*H_invDmH(3,1)+X(1,2)*X(1,3)*H_invDmH(2,3)*H_invDmH(3,1)+X&
&(1,1)*H_invDmH(1,3)*X(1,3)*H_invDmH(3,1)+H_invDmH(2,1)*tt5*H_invD&
&mH(2,3)+H_invDmH(1,1)*X(2,1)*X(2,2)*H_invDmH(2,3)+tt4*H_invDmH(2,&
&1)*H_invDmH(2,3)+H_invDmH(1,1)*X(1,1)*X(1,2)*H_invDmH(2,3)+H_invD&
&mH(1,3)*H_invDmH(2,1)*X(2,1)*X(2,2)+H_invDmH(1,1)*H_invDmH(1,3)*t&
&t3+X(1,1)*X(1,2)*H_invDmH(1,3)*H_invDmH(2,1)+H_invDmH(1,1)*tt2*H_&
&invDmH(1,3))**2/2.0E+0+tt26**2/4.0E+0+(tt25*H_invDmH(8,1)*H_invDm&
&H(8,2)+tt24*H_invDmH(8,1)*H_invDmH(8,2)+tt23*H_invDmH(8,1)*H_invD&
&mH(8,2)+X(3,7)*X(3,8)*H_invDmH(7,1)*H_invDmH(8,2)+X(2,7)*X(2,8)*H&
&_invDmH(7,1)*H_invDmH(8,2)+X(1,7)*X(1,8)*H_invDmH(7,1)*H_invDmH(8&
&,2)+X(3,6)*X(3,8)*H_invDmH(6,1)*H_invDmH(8,2)+X(2,6)*X(2,8)*H_inv&
&DmH(6,1)*H_invDmH(8,2)+X(1,6)*X(1,8)*H_invDmH(6,1)*H_invDmH(8,2)+&
&X(3,5)*X(3,8)*H_invDmH(5,1)*H_invDmH(8,2)+X(2,5)*X(2,8)*H_invDmH(&
&5,1)*H_invDmH(8,2)+X(1,5)*X(1,8)*H_invDmH(5,1)*H_invDmH(8,2)+X(3,&
&4)*X(3,8)*H_invDmH(4,1)*H_invDmH(8,2)+X(2,4)*X(2,8)*H_invDmH(4,1)&
&*H_invDmH(8,2)+X(1,4)*X(1,8)*H_invDmH(4,1)*H_invDmH(8,2)+H_invDmH&
&(3,1)*X(3,3)*X(3,8)*H_invDmH(8,2)+H_invDmH(2,1)*X(3,2)*X(3,8)*H_i&
&nvDmH(8,2)+H_invDmH(1,1)*X(3,1)*X(3,8)*H_invDmH(8,2)+X(2,3)*X(2,8&
&)*H_invDmH(3,1)*H_invDmH(8,2)+X(1,3)*X(1,8)*H_invDmH(3,1)*H_invDm&
&H(8,2)+H_invDmH(2,1)*X(2,2)*X(2,8)*H_invDmH(8,2)+H_invDmH(1,1)*X(&
&2,1)*X(2,8)*H_invDmH(8,2)+X(1,2)*X(1,8)*H_invDmH(2,1)*H_invDmH(8,&
&2)+H_invDmH(1,1)*X(1,1)*X(1,8)*H_invDmH(8,2)+X(3,7)*X(3,8)*H_invD&
&mH(7,2)*H_invDmH(8,1)+X(2,7)*X(2,8)*H_invDmH(7,2)*H_invDmH(8,1)+X&
&(1,7)*X(1,8)*H_invDmH(7,2)*H_invDmH(8,1)+X(3,6)*X(3,8)*H_invDmH(6&
&,2)*H_invDmH(8,1)+X(2,6)*X(2,8)*H_invDmH(6,2)*H_invDmH(8,1)+X(1,6&
&)*X(1,8)*H_invDmH(6,2)*H_invDmH(8,1)+X(3,5)*X(3,8)*H_invDmH(5,2)*&
&H_invDmH(8,1)+X(2,5)*X(2,8)*H_invDmH(5,2)*H_invDmH(8,1)+X(1,5)*X(&
&1,8)*H_invDmH(5,2)*H_invDmH(8,1)+X(3,4)*X(3,8)*H_invDmH(4,2)*H_in&
&vDmH(8,1)+X(2,4)*X(2,8)*H_invDmH(4,2)*H_invDmH(8,1)+X(1,4)*X(1,8)&
&*H_invDmH(4,2)*H_invDmH(8,1)+H_invDmH(3,2)*X(3,3)*X(3,8)*H_invDmH&
&(8,1)+H_invDmH(2,2)*X(3,2)*X(3,8)*H_invDmH(8,1)+H_invDmH(1,2)*X(3&
&,1)*X(3,8)*H_invDmH(8,1)+X(2,3)*X(2,8)*H_invDmH(3,2)*H_invDmH(8,1&
&)+X(1,3)*X(1,8)*H_invDmH(3,2)*H_invDmH(8,1)+H_invDmH(2,2)*X(2,2)*&
&X(2,8)*H_invDmH(8,1)+H_invDmH(1,2)*X(2,1)*X(2,8)*H_invDmH(8,1)+X(&
&1,2)*X(1,8)*H_invDmH(2,2)*H_invDmH(8,1)+X(1,1)*H_invDmH(1,2)*X(1,&
&8)*H_invDmH(8,1)+tt22*H_invDmH(7,1)*H_invDmH(7,2)+tt21*H_invDmH(7&
&,1)*H_invDmH(7,2)+tt20*H_invDmH(7,1)*H_invDmH(7,2)+X(3,6)*X(3,7)*&
&H_invDmH(6,1)*H_invDmH(7,2)+X(2,6)*X(2,7)*H_invDmH(6,1)*H_invDmH(&
&7,2)+X(1,6)*X(1,7)*H_invDmH(6,1)*H_invDmH(7,2)+X(3,5)*X(3,7)*H_in&
&vDmH(5,1)*H_invDmH(7,2)+X(2,5)*X(2,7)*H_invDmH(5,1)*H_invDmH(7,2)&
&+X(1,5)*X(1,7)*H_invDmH(5,1)*H_invDmH(7,2)+X(3,4)*X(3,7)*H_invDmH&
&(4,1)*H_invDmH(7,2)+X(2,4)*X(2,7)*H_invDmH(4,1)*H_invDmH(7,2)+X(1&
&,4)*X(1,7)*H_invDmH(4,1)*H_invDmH(7,2)+H_invDmH(3,1)*X(3,3)*X(3,7&
&)*H_invDmH(7,2)+H_invDmH(2,1)*X(3,2)*X(3,7)*H_invDmH(7,2)+H_invDm&
&H(1,1)*X(3,1)*X(3,7)*H_invDmH(7,2)+X(2,3)*X(2,7)*H_invDmH(3,1)*H_&
&invDmH(7,2)+X(1,3)*X(1,7)*H_invDmH(3,1)*H_invDmH(7,2)+H_invDmH(2,&
&1)*X(2,2)*X(2,7)*H_invDmH(7,2)+H_invDmH(1,1)*X(2,1)*X(2,7)*H_invD&
&mH(7,2)+X(1,2)*X(1,7)*H_invDmH(2,1)*H_invDmH(7,2)+H_invDmH(1,1)*X&
&(1,1)*X(1,7)*H_invDmH(7,2)+X(3,6)*X(3,7)*H_invDmH(6,2)*H_invDmH(7&
&,1)+X(2,6)*X(2,7)*H_invDmH(6,2)*H_invDmH(7,1)+X(1,6)*X(1,7)*H_inv&
&DmH(6,2)*H_invDmH(7,1)+X(3,5)*X(3,7)*H_invDmH(5,2)*H_invDmH(7,1)+&
&X(2,5)*X(2,7)*H_invDmH(5,2)*H_invDmH(7,1)+X(1,5)*X(1,7)*H_invDmH(&
&5,2)*H_invDmH(7,1)+X(3,4)*X(3,7)*H_invDmH(4,2)*H_invDmH(7,1)+X(2,&
&4)*X(2,7)*H_invDmH(4,2)*H_invDmH(7,1)+X(1,4)*X(1,7)*H_invDmH(4,2)&
&*H_invDmH(7,1)+H_invDmH(3,2)*X(3,3)*X(3,7)*H_invDmH(7,1)+H_invDmH&
&(2,2)*X(3,2)*X(3,7)*H_invDmH(7,1)+H_invDmH(1,2)*X(3,1)*X(3,7)*H_i&
&nvDmH(7,1)+X(2,3)*X(2,7)*H_invDmH(3,2)*H_invDmH(7,1)+X(1,3)*X(1,7&
&)*H_invDmH(3,2)*H_invDmH(7,1)+H_invDmH(2,2)*X(2,2)*X(2,7)*H_invDm&
&H(7,1)+H_invDmH(1,2)*X(2,1)*X(2,7)*H_invDmH(7,1)+X(1,2)*X(1,7)*H_&
&invDmH(2,2)*H_invDmH(7,1)+X(1,1)*H_invDmH(1,2)*X(1,7)*H_invDmH(7,&
&1)+tt19*H_invDmH(6,1)*H_invDmH(6,2)+tt18*H_invDmH(6,1)*H_invDmH(6&
&,2)+tt17*H_invDmH(6,1)*H_invDmH(6,2)+X(3,5)*X(3,6)*H_invDmH(5,1)*&
&H_invDmH(6,2)+X(2,5)*X(2,6)*H_invDmH(5,1)*H_invDmH(6,2)+X(1,5)*X(&
&1,6)*H_invDmH(5,1)*H_invDmH(6,2)+X(3,4)*X(3,6)*H_invDmH(4,1)*H_in&
&vDmH(6,2)+X(2,4)*X(2,6)*H_invDmH(4,1)*H_invDmH(6,2)+X(1,4)*X(1,6)&
&*H_invDmH(4,1)*H_invDmH(6,2)+H_invDmH(3,1)*X(3,3)*X(3,6)*H_invDmH&
&(6,2)+H_invDmH(2,1)*X(3,2)*X(3,6)*H_invDmH(6,2)+H_invDmH(1,1)*X(3&
&,1)*X(3,6)*H_invDmH(6,2)+X(2,3)*X(2,6)*H_invDmH(3,1)*H_invDmH(6,2&
&)+X(1,3)*X(1,6)*H_invDmH(3,1)*H_invDmH(6,2)+H_invDmH(2,1)*X(2,2)*&
&X(2,6)*H_invDmH(6,2)+H_invDmH(1,1)*X(2,1)*X(2,6)*H_invDmH(6,2)+X(&
&1,2)*X(1,6)*H_invDmH(2,1)*H_invDmH(6,2)+H_invDmH(1,1)*X(1,1)*X(1,&
&6)*H_invDmH(6,2)+X(3,5)*X(3,6)*H_invDmH(5,2)*H_invDmH(6,1)+X(2,5)&
&*X(2,6)*H_invDmH(5,2)*H_invDmH(6,1)+X(1,5)*X(1,6)*H_invDmH(5,2)*H&
&_invDmH(6,1)+X(3,4)*X(3,6)*H_invDmH(4,2)*H_invDmH(6,1)+X(2,4)*X(2&
&,6)*H_invDmH(4,2)*H_invDmH(6,1)+X(1,4)*X(1,6)*H_invDmH(4,2)*H_inv&
&DmH(6,1)+H_invDmH(3,2)*X(3,3)*X(3,6)*H_invDmH(6,1)+H_invDmH(2,2)*&
&X(3,2)*X(3,6)*H_invDmH(6,1)+H_invDmH(1,2)*X(3,1)*X(3,6)*H_invDmH(&
&6,1)+X(2,3)*X(2,6)*H_invDmH(3,2)*H_invDmH(6,1)+X(1,3)*X(1,6)*H_in&
&vDmH(3,2)*H_invDmH(6,1)+H_invDmH(2,2)*X(2,2)*X(2,6)*H_invDmH(6,1)&
&+H_invDmH(1,2)*X(2,1)*X(2,6)*H_invDmH(6,1)+X(1,2)*X(1,6)*H_invDmH&
&(2,2)*H_invDmH(6,1)+X(1,1)*H_invDmH(1,2)*X(1,6)*H_invDmH(6,1)+tt1&
&6*H_invDmH(5,1)*H_invDmH(5,2)+tt15*H_invDmH(5,1)*H_invDmH(5,2)+tt&
&14*H_invDmH(5,1)*H_invDmH(5,2)+X(3,4)*X(3,5)*H_invDmH(4,1)*H_invD&
&mH(5,2)+X(2,4)*X(2,5)*H_invDmH(4,1)*H_invDmH(5,2)+X(1,4)*X(1,5)*H&
&_invDmH(4,1)*H_invDmH(5,2)+H_invDmH(3,1)*X(3,3)*X(3,5)*H_invDmH(5&
&,2)+H_invDmH(2,1)*X(3,2)*X(3,5)*H_invDmH(5,2)+H_invDmH(1,1)*X(3,1&
&)*X(3,5)*H_invDmH(5,2)+X(2,3)*X(2,5)*H_invDmH(3,1)*H_invDmH(5,2)+&
&X(1,3)*X(1,5)*H_invDmH(3,1)*H_invDmH(5,2)+H_invDmH(2,1)*X(2,2)*X(&
&2,5)*H_invDmH(5,2)+H_invDmH(1,1)*X(2,1)*X(2,5)*H_invDmH(5,2)+X(1,&
&2)*X(1,5)*H_invDmH(2,1)*H_invDmH(5,2)+H_invDmH(1,1)*X(1,1)*X(1,5)&
&*H_invDmH(5,2)+X(3,4)*X(3,5)*H_invDmH(4,2)*H_invDmH(5,1)+X(2,4)*X&
&(2,5)*H_invDmH(4,2)*H_invDmH(5,1)+X(1,4)*X(1,5)*H_invDmH(4,2)*H_i&
&nvDmH(5,1)+H_invDmH(3,2)*X(3,3)*X(3,5)*H_invDmH(5,1)+H_invDmH(2,2&
&)*X(3,2)*X(3,5)*H_invDmH(5,1)+H_invDmH(1,2)*X(3,1)*X(3,5)*H_invDm&
&H(5,1)+X(2,3)*X(2,5)*H_invDmH(3,2)*H_invDmH(5,1)+X(1,3)*X(1,5)*H_&
&invDmH(3,2)*H_invDmH(5,1)+H_invDmH(2,2)*X(2,2)*X(2,5)*H_invDmH(5,&
&1)+H_invDmH(1,2)*X(2,1)*X(2,5)*H_invDmH(5,1)+X(1,2)*X(1,5)*H_invD&
&mH(2,2)*H_invDmH(5,1)+X(1,1)*H_invDmH(1,2)*X(1,5)*H_invDmH(5,1)+t&
&t13*H_invDmH(4,1)*H_invDmH(4,2)+tt12*H_invDmH(4,1)*H_invDmH(4,2)+&
&tt11*H_invDmH(4,1)*H_invDmH(4,2)+H_invDmH(3,1)*X(3,3)*X(3,4)*H_in&
&vDmH(4,2)+H_invDmH(2,1)*X(3,2)*X(3,4)*H_invDmH(4,2)+H_invDmH(1,1)&
&*X(3,1)*X(3,4)*H_invDmH(4,2)+X(2,3)*X(2,4)*H_invDmH(3,1)*H_invDmH&
&(4,2)+X(1,3)*X(1,4)*H_invDmH(3,1)*H_invDmH(4,2)+H_invDmH(2,1)*X(2&
&,2)*X(2,4)*H_invDmH(4,2)+H_invDmH(1,1)*X(2,1)*X(2,4)*H_invDmH(4,2&
&)+X(1,2)*X(1,4)*H_invDmH(2,1)*H_invDmH(4,2)+H_invDmH(1,1)*X(1,1)*&
&X(1,4)*H_invDmH(4,2)+H_invDmH(3,2)*X(3,3)*X(3,4)*H_invDmH(4,1)+H_&
&invDmH(2,2)*X(3,2)*X(3,4)*H_invDmH(4,1)+H_invDmH(1,2)*X(3,1)*X(3,&
&4)*H_invDmH(4,1)+X(2,3)*X(2,4)*H_invDmH(3,2)*H_invDmH(4,1)+X(1,3)&
&*X(1,4)*H_invDmH(3,2)*H_invDmH(4,1)+H_invDmH(2,2)*X(2,2)*X(2,4)*H&
&_invDmH(4,1)+H_invDmH(1,2)*X(2,1)*X(2,4)*H_invDmH(4,1)+X(1,2)*X(1&
&,4)*H_invDmH(2,2)*H_invDmH(4,1)+X(1,1)*H_invDmH(1,2)*X(1,4)*H_inv&
&DmH(4,1)+H_invDmH(3,1)*H_invDmH(3,2)*tt10+H_invDmH(2,1)*H_invDmH(&
&3,2)*X(3,2)*X(3,3)+H_invDmH(2,2)*H_invDmH(3,1)*X(3,2)*X(3,3)+H_in&
&vDmH(1,1)*X(3,1)*H_invDmH(3,2)*X(3,3)+H_invDmH(1,2)*H_invDmH(3,1)&
&*X(3,1)*X(3,3)+H_invDmH(2,1)*H_invDmH(2,2)*tt9+H_invDmH(1,1)*H_in&
&vDmH(2,2)*X(3,1)*X(3,2)+H_invDmH(1,2)*H_invDmH(2,1)*X(3,1)*X(3,2)&
&+tt8*H_invDmH(3,1)*H_invDmH(3,2)+tt7*H_invDmH(3,1)*H_invDmH(3,2)+&
&H_invDmH(2,1)*X(2,2)*X(2,3)*H_invDmH(3,2)+H_invDmH(1,1)*X(2,1)*X(&
&2,3)*H_invDmH(3,2)+X(1,2)*X(1,3)*H_invDmH(2,1)*H_invDmH(3,2)+H_in&
&vDmH(1,1)*X(1,1)*X(1,3)*H_invDmH(3,2)+H_invDmH(1,1)*H_invDmH(1,2)&
&*tt6+H_invDmH(2,2)*X(2,2)*X(2,3)*H_invDmH(3,1)+H_invDmH(1,2)*X(2,&
&1)*X(2,3)*H_invDmH(3,1)+X(1,2)*X(1,3)*H_invDmH(2,2)*H_invDmH(3,1)&
&+X(1,1)*H_invDmH(1,2)*X(1,3)*H_invDmH(3,1)+H_invDmH(2,1)*H_invDmH&
&(2,2)*tt5+H_invDmH(1,1)*X(2,1)*H_invDmH(2,2)*X(2,2)+H_invDmH(1,2)&
&*H_invDmH(2,1)*X(2,1)*X(2,2)+tt4*H_invDmH(2,1)*H_invDmH(2,2)+H_in&
&vDmH(1,1)*X(1,1)*X(1,2)*H_invDmH(2,2)+H_invDmH(1,1)*H_invDmH(1,2)&
&*tt3+X(1,1)*H_invDmH(1,2)*X(1,2)*H_invDmH(2,1)+H_invDmH(1,1)*tt2*&
&H_invDmH(1,2))**2/2.0E+0+tt1**2/4.0E+0))
END 
SUBROUTINE vox_stvk_at_quadr_jac(jac, X, H_invDmH, detDmH, gw, mu, lam) 
IMPLICIT NONE 
REAL(KIND=8) jac(1, 24) 
REAL(KIND=8) X(3, 8) 
REAL(KIND=8) H_invDmH(8, 3) 
REAL(KIND=8) detDmH(1, 1) 
REAL(KIND=8) gw(1, 1) 
REAL(KIND=8) mu(1, 1) 
REAL(KIND=8) lam(1, 1) 
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
tt1 = X(1,8)*H_invDmH(8,1)+X(1,7)*H_invDmH(7,1)+X(1,6)*H_invDmH(6&
&,1)+X(1,5)*H_invDmH(5,1)+X(1,4)*H_invDmH(4,1)+X(1,3)*H_invDmH(3,1&
&)+X(1,2)*H_invDmH(2,1)+H_invDmH(1,1)*X(1,1)
tt2 = X(1,8)*H_invDmH(8,2)+X(1,7)*H_invDmH(7,2)+X(1,6)*H_invDmH(6&
&,2)+X(1,5)*H_invDmH(5,2)+X(1,4)*H_invDmH(4,2)+X(1,3)*H_invDmH(3,2&
&)+X(1,2)*H_invDmH(2,2)+X(1,1)*H_invDmH(1,2)
tt3 = X(1,8)*H_invDmH(8,3)+X(1,7)*H_invDmH(7,3)+X(1,6)*H_invDmH(6&
&,3)+X(1,5)*H_invDmH(5,3)+X(1,4)*H_invDmH(4,3)+X(1,3)*H_invDmH(3,3&
&)+X(1,2)*H_invDmH(2,3)+X(1,1)*H_invDmH(1,3)
tt4 = X(2,8)*H_invDmH(8,1)+X(2,7)*H_invDmH(7,1)+X(2,6)*H_invDmH(6&
&,1)+X(2,5)*H_invDmH(5,1)+X(2,4)*H_invDmH(4,1)+X(2,3)*H_invDmH(3,1&
&)+H_invDmH(2,1)*X(2,2)+H_invDmH(1,1)*X(2,1)
tt5 = X(3,8)*H_invDmH(8,1)+X(3,7)*H_invDmH(7,1)+X(3,6)*H_invDmH(6&
&,1)+X(3,5)*H_invDmH(5,1)+X(3,4)*H_invDmH(4,1)+H_invDmH(3,1)*X(3,3&
&)+H_invDmH(2,1)*X(3,2)+H_invDmH(1,1)*X(3,1)
tt6 = tt5**2+tt4**2+tt1**2-1
tt7 = X(2,8)*H_invDmH(8,2)+X(2,7)*H_invDmH(7,2)+X(2,6)*H_invDmH(6&
&,2)+X(2,5)*H_invDmH(5,2)+X(2,4)*H_invDmH(4,2)+X(2,3)*H_invDmH(3,2&
&)+H_invDmH(2,2)*X(2,2)+H_invDmH(1,2)*X(2,1)
tt8 = X(3,8)*H_invDmH(8,2)+X(3,7)*H_invDmH(7,2)+X(3,6)*H_invDmH(6&
&,2)+X(3,5)*H_invDmH(5,2)+X(3,4)*H_invDmH(4,2)+H_invDmH(3,2)*X(3,3&
&)+H_invDmH(2,2)*X(3,2)+H_invDmH(1,2)*X(3,1)
tt9 = tt8**2+tt7**2+tt2**2-1
tt10 = X(2,8)*H_invDmH(8,3)+X(2,7)*H_invDmH(7,3)+X(2,6)*H_invDmH(&
&6,3)+X(2,5)*H_invDmH(5,3)+X(2,4)*H_invDmH(4,3)+X(2,3)*H_invDmH(3,&
&3)+X(2,2)*H_invDmH(2,3)+H_invDmH(1,3)*X(2,1)
tt11 = X(3,8)*H_invDmH(8,3)+X(3,7)*H_invDmH(7,3)+X(3,6)*H_invDmH(&
&6,3)+X(3,5)*H_invDmH(5,3)+X(3,4)*H_invDmH(4,3)+H_invDmH(3,3)*X(3,&
&3)+H_invDmH(2,3)*X(3,2)+H_invDmH(1,3)*X(3,1)
tt12 = tt11**2+tt10**2+tt3**2-1
tt13 = tt12/2.0E+0+tt9/2.0E+0+tt6/2.0E+0
tt14 = X(1,1)**2
tt15 = X(2,1)**2
tt16 = X(1,2)**2
tt17 = X(2,2)**2
tt18 = X(3,1)**2
tt19 = X(1,3)**2
tt20 = X(2,3)**2
tt21 = X(3,2)**2
tt22 = X(3,3)**2
tt23 = X(1,4)**2
tt24 = X(2,4)**2
tt25 = X(3,4)**2
tt26 = X(1,5)**2
tt27 = X(2,5)**2
tt28 = X(3,5)**2
tt29 = X(1,6)**2
tt30 = X(2,6)**2
tt31 = X(3,6)**2
tt32 = X(1,7)**2
tt33 = X(2,7)**2
tt34 = X(3,7)**2
tt35 = X(1,8)**2
tt36 = X(2,8)**2
tt37 = X(3,8)**2
tt38 = tt37*H_invDmH(8,1)*H_invDmH(8,2)+tt36*H_invDmH(8,1)*H_invD&
&mH(8,2)+tt35*H_invDmH(8,1)*H_invDmH(8,2)+X(3,7)*X(3,8)*H_invDmH(7&
&,1)*H_invDmH(8,2)+X(2,7)*X(2,8)*H_invDmH(7,1)*H_invDmH(8,2)+X(1,7&
&)*X(1,8)*H_invDmH(7,1)*H_invDmH(8,2)+X(3,6)*X(3,8)*H_invDmH(6,1)*&
&H_invDmH(8,2)+X(2,6)*X(2,8)*H_invDmH(6,1)*H_invDmH(8,2)+X(1,6)*X(&
&1,8)*H_invDmH(6,1)*H_invDmH(8,2)+X(3,5)*X(3,8)*H_invDmH(5,1)*H_in&
&vDmH(8,2)+X(2,5)*X(2,8)*H_invDmH(5,1)*H_invDmH(8,2)+X(1,5)*X(1,8)&
&*H_invDmH(5,1)*H_invDmH(8,2)+X(3,4)*X(3,8)*H_invDmH(4,1)*H_invDmH&
&(8,2)+X(2,4)*X(2,8)*H_invDmH(4,1)*H_invDmH(8,2)+X(1,4)*X(1,8)*H_i&
&nvDmH(4,1)*H_invDmH(8,2)+H_invDmH(3,1)*X(3,3)*X(3,8)*H_invDmH(8,2&
&)+H_invDmH(2,1)*X(3,2)*X(3,8)*H_invDmH(8,2)+H_invDmH(1,1)*X(3,1)*&
&X(3,8)*H_invDmH(8,2)+X(2,3)*X(2,8)*H_invDmH(3,1)*H_invDmH(8,2)+X(&
&1,3)*X(1,8)*H_invDmH(3,1)*H_invDmH(8,2)+H_invDmH(2,1)*X(2,2)*X(2,&
&8)*H_invDmH(8,2)+H_invDmH(1,1)*X(2,1)*X(2,8)*H_invDmH(8,2)+X(1,2)&
&*X(1,8)*H_invDmH(2,1)*H_invDmH(8,2)+H_invDmH(1,1)*X(1,1)*X(1,8)*H&
&_invDmH(8,2)+X(3,7)*X(3,8)*H_invDmH(7,2)*H_invDmH(8,1)+X(2,7)*X(2&
&,8)*H_invDmH(7,2)*H_invDmH(8,1)+X(1,7)*X(1,8)*H_invDmH(7,2)*H_inv&
&DmH(8,1)+X(3,6)*X(3,8)*H_invDmH(6,2)*H_invDmH(8,1)+X(2,6)*X(2,8)*&
&H_invDmH(6,2)*H_invDmH(8,1)+X(1,6)*X(1,8)*H_invDmH(6,2)*H_invDmH(&
&8,1)+X(3,5)*X(3,8)*H_invDmH(5,2)*H_invDmH(8,1)+X(2,5)*X(2,8)*H_in&
&vDmH(5,2)*H_invDmH(8,1)+X(1,5)*X(1,8)*H_invDmH(5,2)*H_invDmH(8,1)&
&+X(3,4)*X(3,8)*H_invDmH(4,2)*H_invDmH(8,1)+X(2,4)*X(2,8)*H_invDmH&
&(4,2)*H_invDmH(8,1)+X(1,4)*X(1,8)*H_invDmH(4,2)*H_invDmH(8,1)+H_i&
&nvDmH(3,2)*X(3,3)*X(3,8)*H_invDmH(8,1)+H_invDmH(2,2)*X(3,2)*X(3,8&
&)*H_invDmH(8,1)+H_invDmH(1,2)*X(3,1)*X(3,8)*H_invDmH(8,1)+X(2,3)*&
&X(2,8)*H_invDmH(3,2)*H_invDmH(8,1)+X(1,3)*X(1,8)*H_invDmH(3,2)*H_&
&invDmH(8,1)+H_invDmH(2,2)*X(2,2)*X(2,8)*H_invDmH(8,1)+H_invDmH(1,&
&2)*X(2,1)*X(2,8)*H_invDmH(8,1)+X(1,2)*X(1,8)*H_invDmH(2,2)*H_invD&
&mH(8,1)+X(1,1)*H_invDmH(1,2)*X(1,8)*H_invDmH(8,1)+tt34*H_invDmH(7&
&,1)*H_invDmH(7,2)+tt33*H_invDmH(7,1)*H_invDmH(7,2)+tt32*H_invDmH(&
&7,1)*H_invDmH(7,2)+X(3,6)*X(3,7)*H_invDmH(6,1)*H_invDmH(7,2)+X(2,&
&6)*X(2,7)*H_invDmH(6,1)*H_invDmH(7,2)+X(1,6)*X(1,7)*H_invDmH(6,1)&
&*H_invDmH(7,2)+X(3,5)*X(3,7)*H_invDmH(5,1)*H_invDmH(7,2)+X(2,5)*X&
&(2,7)*H_invDmH(5,1)*H_invDmH(7,2)+X(1,5)*X(1,7)*H_invDmH(5,1)*H_i&
&nvDmH(7,2)+X(3,4)*X(3,7)*H_invDmH(4,1)*H_invDmH(7,2)+X(2,4)*X(2,7&
&)*H_invDmH(4,1)*H_invDmH(7,2)+X(1,4)*X(1,7)*H_invDmH(4,1)*H_invDm&
&H(7,2)+H_invDmH(3,1)*X(3,3)*X(3,7)*H_invDmH(7,2)+H_invDmH(2,1)*X(&
&3,2)*X(3,7)*H_invDmH(7,2)+H_invDmH(1,1)*X(3,1)*X(3,7)*H_invDmH(7,&
&2)+X(2,3)*X(2,7)*H_invDmH(3,1)*H_invDmH(7,2)+X(1,3)*X(1,7)*H_invD&
&mH(3,1)*H_invDmH(7,2)+H_invDmH(2,1)*X(2,2)*X(2,7)*H_invDmH(7,2)+H&
&_invDmH(1,1)*X(2,1)*X(2,7)*H_invDmH(7,2)+X(1,2)*X(1,7)*H_invDmH(2&
&,1)*H_invDmH(7,2)+H_invDmH(1,1)*X(1,1)*X(1,7)*H_invDmH(7,2)+X(3,6&
&)*X(3,7)*H_invDmH(6,2)*H_invDmH(7,1)+X(2,6)*X(2,7)*H_invDmH(6,2)*&
&H_invDmH(7,1)+X(1,6)*X(1,7)*H_invDmH(6,2)*H_invDmH(7,1)+X(3,5)*X(&
&3,7)*H_invDmH(5,2)*H_invDmH(7,1)+X(2,5)*X(2,7)*H_invDmH(5,2)*H_in&
&vDmH(7,1)+X(1,5)*X(1,7)*H_invDmH(5,2)*H_invDmH(7,1)+X(3,4)*X(3,7)&
&*H_invDmH(4,2)*H_invDmH(7,1)+X(2,4)*X(2,7)*H_invDmH(4,2)*H_invDmH&
&(7,1)+X(1,4)*X(1,7)*H_invDmH(4,2)*H_invDmH(7,1)+H_invDmH(3,2)*X(3&
&,3)*X(3,7)*H_invDmH(7,1)+H_invDmH(2,2)*X(3,2)*X(3,7)*H_invDmH(7,1&
&)+H_invDmH(1,2)*X(3,1)*X(3,7)*H_invDmH(7,1)+X(2,3)*X(2,7)*H_invDm&
&H(3,2)*H_invDmH(7,1)+X(1,3)*X(1,7)*H_invDmH(3,2)*H_invDmH(7,1)+H_&
&invDmH(2,2)*X(2,2)*X(2,7)*H_invDmH(7,1)+H_invDmH(1,2)*X(2,1)*X(2,&
&7)*H_invDmH(7,1)+X(1,2)*X(1,7)*H_invDmH(2,2)*H_invDmH(7,1)+X(1,1)&
&*H_invDmH(1,2)*X(1,7)*H_invDmH(7,1)+tt31*H_invDmH(6,1)*H_invDmH(6&
&,2)+tt30*H_invDmH(6,1)*H_invDmH(6,2)+tt29*H_invDmH(6,1)*H_invDmH(&
&6,2)+X(3,5)*X(3,6)*H_invDmH(5,1)*H_invDmH(6,2)+X(2,5)*X(2,6)*H_in&
&vDmH(5,1)*H_invDmH(6,2)+X(1,5)*X(1,6)*H_invDmH(5,1)*H_invDmH(6,2)&
&+X(3,4)*X(3,6)*H_invDmH(4,1)*H_invDmH(6,2)+X(2,4)*X(2,6)*H_invDmH&
&(4,1)*H_invDmH(6,2)+X(1,4)*X(1,6)*H_invDmH(4,1)*H_invDmH(6,2)+H_i&
&nvDmH(3,1)*X(3,3)*X(3,6)*H_invDmH(6,2)+H_invDmH(2,1)*X(3,2)*X(3,6&
&)*H_invDmH(6,2)+H_invDmH(1,1)*X(3,1)*X(3,6)*H_invDmH(6,2)+X(2,3)*&
&X(2,6)*H_invDmH(3,1)*H_invDmH(6,2)+X(1,3)*X(1,6)*H_invDmH(3,1)*H_&
&invDmH(6,2)+H_invDmH(2,1)*X(2,2)*X(2,6)*H_invDmH(6,2)+H_invDmH(1,&
&1)*X(2,1)*X(2,6)*H_invDmH(6,2)+X(1,2)*X(1,6)*H_invDmH(2,1)*H_invD&
&mH(6,2)+H_invDmH(1,1)*X(1,1)*X(1,6)*H_invDmH(6,2)+X(3,5)*X(3,6)*H&
&_invDmH(5,2)*H_invDmH(6,1)+X(2,5)*X(2,6)*H_invDmH(5,2)*H_invDmH(6&
&,1)+X(1,5)*X(1,6)*H_invDmH(5,2)*H_invDmH(6,1)+X(3,4)*X(3,6)*H_inv&
&DmH(4,2)*H_invDmH(6,1)+X(2,4)*X(2,6)*H_invDmH(4,2)*H_invDmH(6,1)+&
&X(1,4)*X(1,6)*H_invDmH(4,2)*H_invDmH(6,1)+H_invDmH(3,2)*X(3,3)*X(&
&3,6)*H_invDmH(6,1)+H_invDmH(2,2)*X(3,2)*X(3,6)*H_invDmH(6,1)+H_in&
&vDmH(1,2)*X(3,1)*X(3,6)*H_invDmH(6,1)+X(2,3)*X(2,6)*H_invDmH(3,2)&
&*H_invDmH(6,1)+X(1,3)*X(1,6)*H_invDmH(3,2)*H_invDmH(6,1)+H_invDmH&
&(2,2)*X(2,2)*X(2,6)*H_invDmH(6,1)+H_invDmH(1,2)*X(2,1)*X(2,6)*H_i&
&nvDmH(6,1)+X(1,2)*X(1,6)*H_invDmH(2,2)*H_invDmH(6,1)+X(1,1)*H_inv&
&DmH(1,2)*X(1,6)*H_invDmH(6,1)+tt28*H_invDmH(5,1)*H_invDmH(5,2)+tt&
&27*H_invDmH(5,1)*H_invDmH(5,2)+tt26*H_invDmH(5,1)*H_invDmH(5,2)+X&
&(3,4)*X(3,5)*H_invDmH(4,1)*H_invDmH(5,2)+X(2,4)*X(2,5)*H_invDmH(4&
&,1)*H_invDmH(5,2)+X(1,4)*X(1,5)*H_invDmH(4,1)*H_invDmH(5,2)+H_inv&
&DmH(3,1)*X(3,3)*X(3,5)*H_invDmH(5,2)+H_invDmH(2,1)*X(3,2)*X(3,5)*&
&H_invDmH(5,2)+H_invDmH(1,1)*X(3,1)*X(3,5)*H_invDmH(5,2)+X(2,3)*X(&
&2,5)*H_invDmH(3,1)*H_invDmH(5,2)+X(1,3)*X(1,5)*H_invDmH(3,1)*H_in&
&vDmH(5,2)+H_invDmH(2,1)*X(2,2)*X(2,5)*H_invDmH(5,2)+H_invDmH(1,1)&
&*X(2,1)*X(2,5)*H_invDmH(5,2)+X(1,2)*X(1,5)*H_invDmH(2,1)*H_invDmH&
&(5,2)+H_invDmH(1,1)*X(1,1)*X(1,5)*H_invDmH(5,2)+X(3,4)*X(3,5)*H_i&
&nvDmH(4,2)*H_invDmH(5,1)+X(2,4)*X(2,5)*H_invDmH(4,2)*H_invDmH(5,1&
&)+X(1,4)*X(1,5)*H_invDmH(4,2)*H_invDmH(5,1)+H_invDmH(3,2)*X(3,3)*&
&X(3,5)*H_invDmH(5,1)+H_invDmH(2,2)*X(3,2)*X(3,5)*H_invDmH(5,1)+H_&
&invDmH(1,2)*X(3,1)*X(3,5)*H_invDmH(5,1)+X(2,3)*X(2,5)*H_invDmH(3,&
&2)*H_invDmH(5,1)+X(1,3)*X(1,5)*H_invDmH(3,2)*H_invDmH(5,1)+H_invD&
&mH(2,2)*X(2,2)*X(2,5)*H_invDmH(5,1)+H_invDmH(1,2)*X(2,1)*X(2,5)*H&
&_invDmH(5,1)+X(1,2)*X(1,5)*H_invDmH(2,2)*H_invDmH(5,1)+X(1,1)*H_i&
&nvDmH(1,2)*X(1,5)*H_invDmH(5,1)+tt25*H_invDmH(4,1)*H_invDmH(4,2)+&
&tt24*H_invDmH(4,1)*H_invDmH(4,2)+tt23*H_invDmH(4,1)*H_invDmH(4,2)&
&+H_invDmH(3,1)*X(3,3)*X(3,4)*H_invDmH(4,2)+H_invDmH(2,1)*X(3,2)*X&
&(3,4)*H_invDmH(4,2)+H_invDmH(1,1)*X(3,1)*X(3,4)*H_invDmH(4,2)+X(2&
&,3)*X(2,4)*H_invDmH(3,1)*H_invDmH(4,2)+X(1,3)*X(1,4)*H_invDmH(3,1&
&)*H_invDmH(4,2)+H_invDmH(2,1)*X(2,2)*X(2,4)*H_invDmH(4,2)+H_invDm&
&H(1,1)*X(2,1)*X(2,4)*H_invDmH(4,2)+X(1,2)*X(1,4)*H_invDmH(2,1)*H_&
&invDmH(4,2)+H_invDmH(1,1)*X(1,1)*X(1,4)*H_invDmH(4,2)+H_invDmH(3,&
&2)*X(3,3)*X(3,4)*H_invDmH(4,1)+H_invDmH(2,2)*X(3,2)*X(3,4)*H_invD&
&mH(4,1)+H_invDmH(1,2)*X(3,1)*X(3,4)*H_invDmH(4,1)+X(2,3)*X(2,4)*H&
&_invDmH(3,2)*H_invDmH(4,1)+X(1,3)*X(1,4)*H_invDmH(3,2)*H_invDmH(4&
&,1)+H_invDmH(2,2)*X(2,2)*X(2,4)*H_invDmH(4,1)+H_invDmH(1,2)*X(2,1&
&)*X(2,4)*H_invDmH(4,1)+X(1,2)*X(1,4)*H_invDmH(2,2)*H_invDmH(4,1)+&
&X(1,1)*H_invDmH(1,2)*X(1,4)*H_invDmH(4,1)+H_invDmH(3,1)*H_invDmH(&
&3,2)*tt22+H_invDmH(2,1)*H_invDmH(3,2)*X(3,2)*X(3,3)+H_invDmH(2,2)&
&*H_invDmH(3,1)*X(3,2)*X(3,3)+H_invDmH(1,1)*X(3,1)*H_invDmH(3,2)*X&
&(3,3)+H_invDmH(1,2)*H_invDmH(3,1)*X(3,1)*X(3,3)+H_invDmH(2,1)*H_i&
&nvDmH(2,2)*tt21+H_invDmH(1,1)*H_invDmH(2,2)*X(3,1)*X(3,2)+H_invDm&
&H(1,2)*H_invDmH(2,1)*X(3,1)*X(3,2)+tt20*H_invDmH(3,1)*H_invDmH(3,&
&2)+tt19*H_invDmH(3,1)*H_invDmH(3,2)+H_invDmH(2,1)*X(2,2)*X(2,3)*H&
&_invDmH(3,2)+H_invDmH(1,1)*X(2,1)*X(2,3)*H_invDmH(3,2)+X(1,2)*X(1&
&,3)*H_invDmH(2,1)*H_invDmH(3,2)+H_invDmH(1,1)*X(1,1)*X(1,3)*H_inv&
&DmH(3,2)+H_invDmH(1,1)*H_invDmH(1,2)*tt18+H_invDmH(2,2)*X(2,2)*X(&
&2,3)*H_invDmH(3,1)+H_invDmH(1,2)*X(2,1)*X(2,3)*H_invDmH(3,1)+X(1,&
&2)*X(1,3)*H_invDmH(2,2)*H_invDmH(3,1)+X(1,1)*H_invDmH(1,2)*X(1,3)&
&*H_invDmH(3,1)+H_invDmH(2,1)*H_invDmH(2,2)*tt17+H_invDmH(1,1)*X(2&
&,1)*H_invDmH(2,2)*X(2,2)+H_invDmH(1,2)*H_invDmH(2,1)*X(2,1)*X(2,2&
&)+tt16*H_invDmH(2,1)*H_invDmH(2,2)+H_invDmH(1,1)*X(1,1)*X(1,2)*H_&
&invDmH(2,2)+H_invDmH(1,1)*H_invDmH(1,2)*tt15+X(1,1)*H_invDmH(1,2)&
&*X(1,2)*H_invDmH(2,1)+H_invDmH(1,1)*tt14*H_invDmH(1,2)
tt39 = tt37*H_invDmH(8,1)*H_invDmH(8,3)+tt36*H_invDmH(8,1)*H_invD&
&mH(8,3)+tt35*H_invDmH(8,1)*H_invDmH(8,3)+X(3,7)*X(3,8)*H_invDmH(7&
&,1)*H_invDmH(8,3)+X(2,7)*X(2,8)*H_invDmH(7,1)*H_invDmH(8,3)+X(1,7&
&)*X(1,8)*H_invDmH(7,1)*H_invDmH(8,3)+X(3,6)*X(3,8)*H_invDmH(6,1)*&
&H_invDmH(8,3)+X(2,6)*X(2,8)*H_invDmH(6,1)*H_invDmH(8,3)+X(1,6)*X(&
&1,8)*H_invDmH(6,1)*H_invDmH(8,3)+X(3,5)*X(3,8)*H_invDmH(5,1)*H_in&
&vDmH(8,3)+X(2,5)*X(2,8)*H_invDmH(5,1)*H_invDmH(8,3)+X(1,5)*X(1,8)&
&*H_invDmH(5,1)*H_invDmH(8,3)+X(3,4)*X(3,8)*H_invDmH(4,1)*H_invDmH&
&(8,3)+X(2,4)*X(2,8)*H_invDmH(4,1)*H_invDmH(8,3)+X(1,4)*X(1,8)*H_i&
&nvDmH(4,1)*H_invDmH(8,3)+H_invDmH(3,1)*X(3,3)*X(3,8)*H_invDmH(8,3&
&)+H_invDmH(2,1)*X(3,2)*X(3,8)*H_invDmH(8,3)+H_invDmH(1,1)*X(3,1)*&
&X(3,8)*H_invDmH(8,3)+X(2,3)*X(2,8)*H_invDmH(3,1)*H_invDmH(8,3)+X(&
&1,3)*X(1,8)*H_invDmH(3,1)*H_invDmH(8,3)+H_invDmH(2,1)*X(2,2)*X(2,&
&8)*H_invDmH(8,3)+H_invDmH(1,1)*X(2,1)*X(2,8)*H_invDmH(8,3)+X(1,2)&
&*X(1,8)*H_invDmH(2,1)*H_invDmH(8,3)+H_invDmH(1,1)*X(1,1)*X(1,8)*H&
&_invDmH(8,3)+X(3,7)*X(3,8)*H_invDmH(7,3)*H_invDmH(8,1)+X(2,7)*X(2&
&,8)*H_invDmH(7,3)*H_invDmH(8,1)+X(1,7)*X(1,8)*H_invDmH(7,3)*H_inv&
&DmH(8,1)+X(3,6)*X(3,8)*H_invDmH(6,3)*H_invDmH(8,1)+X(2,6)*X(2,8)*&
&H_invDmH(6,3)*H_invDmH(8,1)+X(1,6)*X(1,8)*H_invDmH(6,3)*H_invDmH(&
&8,1)+X(3,5)*X(3,8)*H_invDmH(5,3)*H_invDmH(8,1)+X(2,5)*X(2,8)*H_in&
&vDmH(5,3)*H_invDmH(8,1)+X(1,5)*X(1,8)*H_invDmH(5,3)*H_invDmH(8,1)&
&+X(3,4)*X(3,8)*H_invDmH(4,3)*H_invDmH(8,1)+X(2,4)*X(2,8)*H_invDmH&
&(4,3)*H_invDmH(8,1)+X(1,4)*X(1,8)*H_invDmH(4,3)*H_invDmH(8,1)+H_i&
&nvDmH(3,3)*X(3,3)*X(3,8)*H_invDmH(8,1)+H_invDmH(2,3)*X(3,2)*X(3,8&
&)*H_invDmH(8,1)+H_invDmH(1,3)*X(3,1)*X(3,8)*H_invDmH(8,1)+X(2,3)*&
&X(2,8)*H_invDmH(3,3)*H_invDmH(8,1)+X(1,3)*X(1,8)*H_invDmH(3,3)*H_&
&invDmH(8,1)+X(2,2)*H_invDmH(2,3)*X(2,8)*H_invDmH(8,1)+H_invDmH(1,&
&3)*X(2,1)*X(2,8)*H_invDmH(8,1)+X(1,2)*X(1,8)*H_invDmH(2,3)*H_invD&
&mH(8,1)+X(1,1)*H_invDmH(1,3)*X(1,8)*H_invDmH(8,1)+tt34*H_invDmH(7&
&,1)*H_invDmH(7,3)+tt33*H_invDmH(7,1)*H_invDmH(7,3)+tt32*H_invDmH(&
&7,1)*H_invDmH(7,3)+X(3,6)*X(3,7)*H_invDmH(6,1)*H_invDmH(7,3)+X(2,&
&6)*X(2,7)*H_invDmH(6,1)*H_invDmH(7,3)+X(1,6)*X(1,7)*H_invDmH(6,1)&
&*H_invDmH(7,3)+X(3,5)*X(3,7)*H_invDmH(5,1)*H_invDmH(7,3)+X(2,5)*X&
&(2,7)*H_invDmH(5,1)*H_invDmH(7,3)+X(1,5)*X(1,7)*H_invDmH(5,1)*H_i&
&nvDmH(7,3)+X(3,4)*X(3,7)*H_invDmH(4,1)*H_invDmH(7,3)+X(2,4)*X(2,7&
&)*H_invDmH(4,1)*H_invDmH(7,3)+X(1,4)*X(1,7)*H_invDmH(4,1)*H_invDm&
&H(7,3)+H_invDmH(3,1)*X(3,3)*X(3,7)*H_invDmH(7,3)+H_invDmH(2,1)*X(&
&3,2)*X(3,7)*H_invDmH(7,3)+H_invDmH(1,1)*X(3,1)*X(3,7)*H_invDmH(7,&
&3)+X(2,3)*X(2,7)*H_invDmH(3,1)*H_invDmH(7,3)+X(1,3)*X(1,7)*H_invD&
&mH(3,1)*H_invDmH(7,3)+H_invDmH(2,1)*X(2,2)*X(2,7)*H_invDmH(7,3)+H&
&_invDmH(1,1)*X(2,1)*X(2,7)*H_invDmH(7,3)+X(1,2)*X(1,7)*H_invDmH(2&
&,1)*H_invDmH(7,3)+H_invDmH(1,1)*X(1,1)*X(1,7)*H_invDmH(7,3)+X(3,6&
&)*X(3,7)*H_invDmH(6,3)*H_invDmH(7,1)+X(2,6)*X(2,7)*H_invDmH(6,3)*&
&H_invDmH(7,1)+X(1,6)*X(1,7)*H_invDmH(6,3)*H_invDmH(7,1)+X(3,5)*X(&
&3,7)*H_invDmH(5,3)*H_invDmH(7,1)+X(2,5)*X(2,7)*H_invDmH(5,3)*H_in&
&vDmH(7,1)+X(1,5)*X(1,7)*H_invDmH(5,3)*H_invDmH(7,1)+X(3,4)*X(3,7)&
&*H_invDmH(4,3)*H_invDmH(7,1)+X(2,4)*X(2,7)*H_invDmH(4,3)*H_invDmH&
&(7,1)+X(1,4)*X(1,7)*H_invDmH(4,3)*H_invDmH(7,1)+H_invDmH(3,3)*X(3&
&,3)*X(3,7)*H_invDmH(7,1)+H_invDmH(2,3)*X(3,2)*X(3,7)*H_invDmH(7,1&
&)+H_invDmH(1,3)*X(3,1)*X(3,7)*H_invDmH(7,1)+X(2,3)*X(2,7)*H_invDm&
&H(3,3)*H_invDmH(7,1)+X(1,3)*X(1,7)*H_invDmH(3,3)*H_invDmH(7,1)+X(&
&2,2)*H_invDmH(2,3)*X(2,7)*H_invDmH(7,1)+H_invDmH(1,3)*X(2,1)*X(2,&
&7)*H_invDmH(7,1)+X(1,2)*X(1,7)*H_invDmH(2,3)*H_invDmH(7,1)+X(1,1)&
&*H_invDmH(1,3)*X(1,7)*H_invDmH(7,1)+tt31*H_invDmH(6,1)*H_invDmH(6&
&,3)+tt30*H_invDmH(6,1)*H_invDmH(6,3)+tt29*H_invDmH(6,1)*H_invDmH(&
&6,3)+X(3,5)*X(3,6)*H_invDmH(5,1)*H_invDmH(6,3)+X(2,5)*X(2,6)*H_in&
&vDmH(5,1)*H_invDmH(6,3)+X(1,5)*X(1,6)*H_invDmH(5,1)*H_invDmH(6,3)&
&+X(3,4)*X(3,6)*H_invDmH(4,1)*H_invDmH(6,3)+X(2,4)*X(2,6)*H_invDmH&
&(4,1)*H_invDmH(6,3)+X(1,4)*X(1,6)*H_invDmH(4,1)*H_invDmH(6,3)+H_i&
&nvDmH(3,1)*X(3,3)*X(3,6)*H_invDmH(6,3)+H_invDmH(2,1)*X(3,2)*X(3,6&
&)*H_invDmH(6,3)+H_invDmH(1,1)*X(3,1)*X(3,6)*H_invDmH(6,3)+X(2,3)*&
&X(2,6)*H_invDmH(3,1)*H_invDmH(6,3)+X(1,3)*X(1,6)*H_invDmH(3,1)*H_&
&invDmH(6,3)+H_invDmH(2,1)*X(2,2)*X(2,6)*H_invDmH(6,3)+H_invDmH(1,&
&1)*X(2,1)*X(2,6)*H_invDmH(6,3)+X(1,2)*X(1,6)*H_invDmH(2,1)*H_invD&
&mH(6,3)+H_invDmH(1,1)*X(1,1)*X(1,6)*H_invDmH(6,3)+X(3,5)*X(3,6)*H&
&_invDmH(5,3)*H_invDmH(6,1)+X(2,5)*X(2,6)*H_invDmH(5,3)*H_invDmH(6&
&,1)+X(1,5)*X(1,6)*H_invDmH(5,3)*H_invDmH(6,1)+X(3,4)*X(3,6)*H_inv&
&DmH(4,3)*H_invDmH(6,1)+X(2,4)*X(2,6)*H_invDmH(4,3)*H_invDmH(6,1)+&
&X(1,4)*X(1,6)*H_invDmH(4,3)*H_invDmH(6,1)+H_invDmH(3,3)*X(3,3)*X(&
&3,6)*H_invDmH(6,1)+H_invDmH(2,3)*X(3,2)*X(3,6)*H_invDmH(6,1)+H_in&
&vDmH(1,3)*X(3,1)*X(3,6)*H_invDmH(6,1)+X(2,3)*X(2,6)*H_invDmH(3,3)&
&*H_invDmH(6,1)+X(1,3)*X(1,6)*H_invDmH(3,3)*H_invDmH(6,1)+X(2,2)*H&
&_invDmH(2,3)*X(2,6)*H_invDmH(6,1)+H_invDmH(1,3)*X(2,1)*X(2,6)*H_i&
&nvDmH(6,1)+X(1,2)*X(1,6)*H_invDmH(2,3)*H_invDmH(6,1)+X(1,1)*H_inv&
&DmH(1,3)*X(1,6)*H_invDmH(6,1)+tt28*H_invDmH(5,1)*H_invDmH(5,3)+tt&
&27*H_invDmH(5,1)*H_invDmH(5,3)+tt26*H_invDmH(5,1)*H_invDmH(5,3)+X&
&(3,4)*X(3,5)*H_invDmH(4,1)*H_invDmH(5,3)+X(2,4)*X(2,5)*H_invDmH(4&
&,1)*H_invDmH(5,3)+X(1,4)*X(1,5)*H_invDmH(4,1)*H_invDmH(5,3)+H_inv&
&DmH(3,1)*X(3,3)*X(3,5)*H_invDmH(5,3)+H_invDmH(2,1)*X(3,2)*X(3,5)*&
&H_invDmH(5,3)+H_invDmH(1,1)*X(3,1)*X(3,5)*H_invDmH(5,3)+X(2,3)*X(&
&2,5)*H_invDmH(3,1)*H_invDmH(5,3)+X(1,3)*X(1,5)*H_invDmH(3,1)*H_in&
&vDmH(5,3)+H_invDmH(2,1)*X(2,2)*X(2,5)*H_invDmH(5,3)+H_invDmH(1,1)&
&*X(2,1)*X(2,5)*H_invDmH(5,3)+X(1,2)*X(1,5)*H_invDmH(2,1)*H_invDmH&
&(5,3)+H_invDmH(1,1)*X(1,1)*X(1,5)*H_invDmH(5,3)+X(3,4)*X(3,5)*H_i&
&nvDmH(4,3)*H_invDmH(5,1)+X(2,4)*X(2,5)*H_invDmH(4,3)*H_invDmH(5,1&
&)+X(1,4)*X(1,5)*H_invDmH(4,3)*H_invDmH(5,1)+H_invDmH(3,3)*X(3,3)*&
&X(3,5)*H_invDmH(5,1)+H_invDmH(2,3)*X(3,2)*X(3,5)*H_invDmH(5,1)+H_&
&invDmH(1,3)*X(3,1)*X(3,5)*H_invDmH(5,1)+X(2,3)*X(2,5)*H_invDmH(3,&
&3)*H_invDmH(5,1)+X(1,3)*X(1,5)*H_invDmH(3,3)*H_invDmH(5,1)+X(2,2)&
&*H_invDmH(2,3)*X(2,5)*H_invDmH(5,1)+H_invDmH(1,3)*X(2,1)*X(2,5)*H&
&_invDmH(5,1)+X(1,2)*X(1,5)*H_invDmH(2,3)*H_invDmH(5,1)+X(1,1)*H_i&
&nvDmH(1,3)*X(1,5)*H_invDmH(5,1)+tt25*H_invDmH(4,1)*H_invDmH(4,3)+&
&tt24*H_invDmH(4,1)*H_invDmH(4,3)+tt23*H_invDmH(4,1)*H_invDmH(4,3)&
&+H_invDmH(3,1)*X(3,3)*X(3,4)*H_invDmH(4,3)+H_invDmH(2,1)*X(3,2)*X&
&(3,4)*H_invDmH(4,3)+H_invDmH(1,1)*X(3,1)*X(3,4)*H_invDmH(4,3)+X(2&
&,3)*X(2,4)*H_invDmH(3,1)*H_invDmH(4,3)+X(1,3)*X(1,4)*H_invDmH(3,1&
&)*H_invDmH(4,3)+H_invDmH(2,1)*X(2,2)*X(2,4)*H_invDmH(4,3)+H_invDm&
&H(1,1)*X(2,1)*X(2,4)*H_invDmH(4,3)+X(1,2)*X(1,4)*H_invDmH(2,1)*H_&
&invDmH(4,3)+H_invDmH(1,1)*X(1,1)*X(1,4)*H_invDmH(4,3)+H_invDmH(3,&
&3)*X(3,3)*X(3,4)*H_invDmH(4,1)+H_invDmH(2,3)*X(3,2)*X(3,4)*H_invD&
&mH(4,1)+H_invDmH(1,3)*X(3,1)*X(3,4)*H_invDmH(4,1)+X(2,3)*X(2,4)*H&
&_invDmH(3,3)*H_invDmH(4,1)+X(1,3)*X(1,4)*H_invDmH(3,3)*H_invDmH(4&
&,1)+X(2,2)*H_invDmH(2,3)*X(2,4)*H_invDmH(4,1)+H_invDmH(1,3)*X(2,1&
&)*X(2,4)*H_invDmH(4,1)+X(1,2)*X(1,4)*H_invDmH(2,3)*H_invDmH(4,1)+&
&X(1,1)*H_invDmH(1,3)*X(1,4)*H_invDmH(4,1)+H_invDmH(3,1)*H_invDmH(&
&3,3)*tt22+H_invDmH(2,1)*X(3,2)*H_invDmH(3,3)*X(3,3)+H_invDmH(1,1)&
&*X(3,1)*H_invDmH(3,3)*X(3,3)+H_invDmH(2,3)*H_invDmH(3,1)*X(3,2)*X&
&(3,3)+H_invDmH(1,3)*H_invDmH(3,1)*X(3,1)*X(3,3)+tt20*H_invDmH(3,1&
&)*H_invDmH(3,3)+tt19*H_invDmH(3,1)*H_invDmH(3,3)+H_invDmH(2,1)*X(&
&2,2)*X(2,3)*H_invDmH(3,3)+H_invDmH(1,1)*X(2,1)*X(2,3)*H_invDmH(3,&
&3)+X(1,2)*X(1,3)*H_invDmH(2,1)*H_invDmH(3,3)+H_invDmH(1,1)*X(1,1)&
&*X(1,3)*H_invDmH(3,3)+H_invDmH(2,1)*H_invDmH(2,3)*tt21+H_invDmH(1&
&,1)*H_invDmH(2,3)*X(3,1)*X(3,2)+H_invDmH(1,3)*H_invDmH(2,1)*X(3,1&
&)*X(3,2)+H_invDmH(1,1)*H_invDmH(1,3)*tt18+X(2,2)*H_invDmH(2,3)*X(&
&2,3)*H_invDmH(3,1)+H_invDmH(1,3)*X(2,1)*X(2,3)*H_invDmH(3,1)+X(1,&
&2)*X(1,3)*H_invDmH(2,3)*H_invDmH(3,1)+X(1,1)*H_invDmH(1,3)*X(1,3)&
&*H_invDmH(3,1)+H_invDmH(2,1)*tt17*H_invDmH(2,3)+H_invDmH(1,1)*X(2&
&,1)*X(2,2)*H_invDmH(2,3)+tt16*H_invDmH(2,1)*H_invDmH(2,3)+H_invDm&
&H(1,1)*X(1,1)*X(1,2)*H_invDmH(2,3)+H_invDmH(1,3)*H_invDmH(2,1)*X(&
&2,1)*X(2,2)+H_invDmH(1,1)*H_invDmH(1,3)*tt15+X(1,1)*X(1,2)*H_invD&
&mH(1,3)*H_invDmH(2,1)+H_invDmH(1,1)*tt14*H_invDmH(1,3)
tt40 = tt37*H_invDmH(8,2)*H_invDmH(8,3)+tt36*H_invDmH(8,2)*H_invD&
&mH(8,3)+tt35*H_invDmH(8,2)*H_invDmH(8,3)+X(3,7)*X(3,8)*H_invDmH(7&
&,2)*H_invDmH(8,3)+X(2,7)*X(2,8)*H_invDmH(7,2)*H_invDmH(8,3)+X(1,7&
&)*X(1,8)*H_invDmH(7,2)*H_invDmH(8,3)+X(3,6)*X(3,8)*H_invDmH(6,2)*&
&H_invDmH(8,3)+X(2,6)*X(2,8)*H_invDmH(6,2)*H_invDmH(8,3)+X(1,6)*X(&
&1,8)*H_invDmH(6,2)*H_invDmH(8,3)+X(3,5)*X(3,8)*H_invDmH(5,2)*H_in&
&vDmH(8,3)+X(2,5)*X(2,8)*H_invDmH(5,2)*H_invDmH(8,3)+X(1,5)*X(1,8)&
&*H_invDmH(5,2)*H_invDmH(8,3)+X(3,4)*X(3,8)*H_invDmH(4,2)*H_invDmH&
&(8,3)+X(2,4)*X(2,8)*H_invDmH(4,2)*H_invDmH(8,3)+X(1,4)*X(1,8)*H_i&
&nvDmH(4,2)*H_invDmH(8,3)+H_invDmH(3,2)*X(3,3)*X(3,8)*H_invDmH(8,3&
&)+H_invDmH(2,2)*X(3,2)*X(3,8)*H_invDmH(8,3)+H_invDmH(1,2)*X(3,1)*&
&X(3,8)*H_invDmH(8,3)+X(2,3)*X(2,8)*H_invDmH(3,2)*H_invDmH(8,3)+X(&
&1,3)*X(1,8)*H_invDmH(3,2)*H_invDmH(8,3)+H_invDmH(2,2)*X(2,2)*X(2,&
&8)*H_invDmH(8,3)+H_invDmH(1,2)*X(2,1)*X(2,8)*H_invDmH(8,3)+X(1,2)&
&*X(1,8)*H_invDmH(2,2)*H_invDmH(8,3)+X(1,1)*H_invDmH(1,2)*X(1,8)*H&
&_invDmH(8,3)+X(3,7)*X(3,8)*H_invDmH(7,3)*H_invDmH(8,2)+X(2,7)*X(2&
&,8)*H_invDmH(7,3)*H_invDmH(8,2)+X(1,7)*X(1,8)*H_invDmH(7,3)*H_inv&
&DmH(8,2)+X(3,6)*X(3,8)*H_invDmH(6,3)*H_invDmH(8,2)+X(2,6)*X(2,8)*&
&H_invDmH(6,3)*H_invDmH(8,2)+X(1,6)*X(1,8)*H_invDmH(6,3)*H_invDmH(&
&8,2)+X(3,5)*X(3,8)*H_invDmH(5,3)*H_invDmH(8,2)+X(2,5)*X(2,8)*H_in&
&vDmH(5,3)*H_invDmH(8,2)+X(1,5)*X(1,8)*H_invDmH(5,3)*H_invDmH(8,2)&
&+X(3,4)*X(3,8)*H_invDmH(4,3)*H_invDmH(8,2)+X(2,4)*X(2,8)*H_invDmH&
&(4,3)*H_invDmH(8,2)+X(1,4)*X(1,8)*H_invDmH(4,3)*H_invDmH(8,2)+H_i&
&nvDmH(3,3)*X(3,3)*X(3,8)*H_invDmH(8,2)+H_invDmH(2,3)*X(3,2)*X(3,8&
&)*H_invDmH(8,2)+H_invDmH(1,3)*X(3,1)*X(3,8)*H_invDmH(8,2)+X(2,3)*&
&X(2,8)*H_invDmH(3,3)*H_invDmH(8,2)+X(1,3)*X(1,8)*H_invDmH(3,3)*H_&
&invDmH(8,2)+X(2,2)*H_invDmH(2,3)*X(2,8)*H_invDmH(8,2)+H_invDmH(1,&
&3)*X(2,1)*X(2,8)*H_invDmH(8,2)+X(1,2)*X(1,8)*H_invDmH(2,3)*H_invD&
&mH(8,2)+X(1,1)*H_invDmH(1,3)*X(1,8)*H_invDmH(8,2)+tt34*H_invDmH(7&
&,2)*H_invDmH(7,3)+tt33*H_invDmH(7,2)*H_invDmH(7,3)+tt32*H_invDmH(&
&7,2)*H_invDmH(7,3)+X(3,6)*X(3,7)*H_invDmH(6,2)*H_invDmH(7,3)+X(2,&
&6)*X(2,7)*H_invDmH(6,2)*H_invDmH(7,3)+X(1,6)*X(1,7)*H_invDmH(6,2)&
&*H_invDmH(7,3)+X(3,5)*X(3,7)*H_invDmH(5,2)*H_invDmH(7,3)+X(2,5)*X&
&(2,7)*H_invDmH(5,2)*H_invDmH(7,3)+X(1,5)*X(1,7)*H_invDmH(5,2)*H_i&
&nvDmH(7,3)+X(3,4)*X(3,7)*H_invDmH(4,2)*H_invDmH(7,3)+X(2,4)*X(2,7&
&)*H_invDmH(4,2)*H_invDmH(7,3)+X(1,4)*X(1,7)*H_invDmH(4,2)*H_invDm&
&H(7,3)+H_invDmH(3,2)*X(3,3)*X(3,7)*H_invDmH(7,3)+H_invDmH(2,2)*X(&
&3,2)*X(3,7)*H_invDmH(7,3)+H_invDmH(1,2)*X(3,1)*X(3,7)*H_invDmH(7,&
&3)+X(2,3)*X(2,7)*H_invDmH(3,2)*H_invDmH(7,3)+X(1,3)*X(1,7)*H_invD&
&mH(3,2)*H_invDmH(7,3)+H_invDmH(2,2)*X(2,2)*X(2,7)*H_invDmH(7,3)+H&
&_invDmH(1,2)*X(2,1)*X(2,7)*H_invDmH(7,3)+X(1,2)*X(1,7)*H_invDmH(2&
&,2)*H_invDmH(7,3)+X(1,1)*H_invDmH(1,2)*X(1,7)*H_invDmH(7,3)+X(3,6&
&)*X(3,7)*H_invDmH(6,3)*H_invDmH(7,2)+X(2,6)*X(2,7)*H_invDmH(6,3)*&
&H_invDmH(7,2)+X(1,6)*X(1,7)*H_invDmH(6,3)*H_invDmH(7,2)+X(3,5)*X(&
&3,7)*H_invDmH(5,3)*H_invDmH(7,2)+X(2,5)*X(2,7)*H_invDmH(5,3)*H_in&
&vDmH(7,2)+X(1,5)*X(1,7)*H_invDmH(5,3)*H_invDmH(7,2)+X(3,4)*X(3,7)&
&*H_invDmH(4,3)*H_invDmH(7,2)+X(2,4)*X(2,7)*H_invDmH(4,3)*H_invDmH&
&(7,2)+X(1,4)*X(1,7)*H_invDmH(4,3)*H_invDmH(7,2)+H_invDmH(3,3)*X(3&
&,3)*X(3,7)*H_invDmH(7,2)+H_invDmH(2,3)*X(3,2)*X(3,7)*H_invDmH(7,2&
&)+H_invDmH(1,3)*X(3,1)*X(3,7)*H_invDmH(7,2)+X(2,3)*X(2,7)*H_invDm&
&H(3,3)*H_invDmH(7,2)+X(1,3)*X(1,7)*H_invDmH(3,3)*H_invDmH(7,2)+X(&
&2,2)*H_invDmH(2,3)*X(2,7)*H_invDmH(7,2)+H_invDmH(1,3)*X(2,1)*X(2,&
&7)*H_invDmH(7,2)+X(1,2)*X(1,7)*H_invDmH(2,3)*H_invDmH(7,2)+X(1,1)&
&*H_invDmH(1,3)*X(1,7)*H_invDmH(7,2)+tt31*H_invDmH(6,2)*H_invDmH(6&
&,3)+tt30*H_invDmH(6,2)*H_invDmH(6,3)+tt29*H_invDmH(6,2)*H_invDmH(&
&6,3)+X(3,5)*X(3,6)*H_invDmH(5,2)*H_invDmH(6,3)+X(2,5)*X(2,6)*H_in&
&vDmH(5,2)*H_invDmH(6,3)+X(1,5)*X(1,6)*H_invDmH(5,2)*H_invDmH(6,3)&
&+X(3,4)*X(3,6)*H_invDmH(4,2)*H_invDmH(6,3)+X(2,4)*X(2,6)*H_invDmH&
&(4,2)*H_invDmH(6,3)+X(1,4)*X(1,6)*H_invDmH(4,2)*H_invDmH(6,3)+H_i&
&nvDmH(3,2)*X(3,3)*X(3,6)*H_invDmH(6,3)+H_invDmH(2,2)*X(3,2)*X(3,6&
&)*H_invDmH(6,3)+H_invDmH(1,2)*X(3,1)*X(3,6)*H_invDmH(6,3)+X(2,3)*&
&X(2,6)*H_invDmH(3,2)*H_invDmH(6,3)+X(1,3)*X(1,6)*H_invDmH(3,2)*H_&
&invDmH(6,3)+H_invDmH(2,2)*X(2,2)*X(2,6)*H_invDmH(6,3)+H_invDmH(1,&
&2)*X(2,1)*X(2,6)*H_invDmH(6,3)+X(1,2)*X(1,6)*H_invDmH(2,2)*H_invD&
&mH(6,3)+X(1,1)*H_invDmH(1,2)*X(1,6)*H_invDmH(6,3)+X(3,5)*X(3,6)*H&
&_invDmH(5,3)*H_invDmH(6,2)+X(2,5)*X(2,6)*H_invDmH(5,3)*H_invDmH(6&
&,2)+X(1,5)*X(1,6)*H_invDmH(5,3)*H_invDmH(6,2)+X(3,4)*X(3,6)*H_inv&
&DmH(4,3)*H_invDmH(6,2)+X(2,4)*X(2,6)*H_invDmH(4,3)*H_invDmH(6,2)+&
&X(1,4)*X(1,6)*H_invDmH(4,3)*H_invDmH(6,2)+H_invDmH(3,3)*X(3,3)*X(&
&3,6)*H_invDmH(6,2)+H_invDmH(2,3)*X(3,2)*X(3,6)*H_invDmH(6,2)+H_in&
&vDmH(1,3)*X(3,1)*X(3,6)*H_invDmH(6,2)+X(2,3)*X(2,6)*H_invDmH(3,3)&
&*H_invDmH(6,2)+X(1,3)*X(1,6)*H_invDmH(3,3)*H_invDmH(6,2)+X(2,2)*H&
&_invDmH(2,3)*X(2,6)*H_invDmH(6,2)+H_invDmH(1,3)*X(2,1)*X(2,6)*H_i&
&nvDmH(6,2)+X(1,2)*X(1,6)*H_invDmH(2,3)*H_invDmH(6,2)+X(1,1)*H_inv&
&DmH(1,3)*X(1,6)*H_invDmH(6,2)+tt28*H_invDmH(5,2)*H_invDmH(5,3)+tt&
&27*H_invDmH(5,2)*H_invDmH(5,3)+tt26*H_invDmH(5,2)*H_invDmH(5,3)+X&
&(3,4)*X(3,5)*H_invDmH(4,2)*H_invDmH(5,3)+X(2,4)*X(2,5)*H_invDmH(4&
&,2)*H_invDmH(5,3)+X(1,4)*X(1,5)*H_invDmH(4,2)*H_invDmH(5,3)+H_inv&
&DmH(3,2)*X(3,3)*X(3,5)*H_invDmH(5,3)+H_invDmH(2,2)*X(3,2)*X(3,5)*&
&H_invDmH(5,3)+H_invDmH(1,2)*X(3,1)*X(3,5)*H_invDmH(5,3)+X(2,3)*X(&
&2,5)*H_invDmH(3,2)*H_invDmH(5,3)+X(1,3)*X(1,5)*H_invDmH(3,2)*H_in&
&vDmH(5,3)+H_invDmH(2,2)*X(2,2)*X(2,5)*H_invDmH(5,3)+H_invDmH(1,2)&
&*X(2,1)*X(2,5)*H_invDmH(5,3)+X(1,2)*X(1,5)*H_invDmH(2,2)*H_invDmH&
&(5,3)+X(1,1)*H_invDmH(1,2)*X(1,5)*H_invDmH(5,3)+X(3,4)*X(3,5)*H_i&
&nvDmH(4,3)*H_invDmH(5,2)+X(2,4)*X(2,5)*H_invDmH(4,3)*H_invDmH(5,2&
&)+X(1,4)*X(1,5)*H_invDmH(4,3)*H_invDmH(5,2)+H_invDmH(3,3)*X(3,3)*&
&X(3,5)*H_invDmH(5,2)+H_invDmH(2,3)*X(3,2)*X(3,5)*H_invDmH(5,2)+H_&
&invDmH(1,3)*X(3,1)*X(3,5)*H_invDmH(5,2)+X(2,3)*X(2,5)*H_invDmH(3,&
&3)*H_invDmH(5,2)+X(1,3)*X(1,5)*H_invDmH(3,3)*H_invDmH(5,2)+X(2,2)&
&*H_invDmH(2,3)*X(2,5)*H_invDmH(5,2)+H_invDmH(1,3)*X(2,1)*X(2,5)*H&
&_invDmH(5,2)+X(1,2)*X(1,5)*H_invDmH(2,3)*H_invDmH(5,2)+X(1,1)*H_i&
&nvDmH(1,3)*X(1,5)*H_invDmH(5,2)+tt25*H_invDmH(4,2)*H_invDmH(4,3)+&
&tt24*H_invDmH(4,2)*H_invDmH(4,3)+tt23*H_invDmH(4,2)*H_invDmH(4,3)&
&+H_invDmH(3,2)*X(3,3)*X(3,4)*H_invDmH(4,3)+H_invDmH(2,2)*X(3,2)*X&
&(3,4)*H_invDmH(4,3)+H_invDmH(1,2)*X(3,1)*X(3,4)*H_invDmH(4,3)+X(2&
&,3)*X(2,4)*H_invDmH(3,2)*H_invDmH(4,3)+X(1,3)*X(1,4)*H_invDmH(3,2&
&)*H_invDmH(4,3)+H_invDmH(2,2)*X(2,2)*X(2,4)*H_invDmH(4,3)+H_invDm&
&H(1,2)*X(2,1)*X(2,4)*H_invDmH(4,3)+X(1,2)*X(1,4)*H_invDmH(2,2)*H_&
&invDmH(4,3)+X(1,1)*H_invDmH(1,2)*X(1,4)*H_invDmH(4,3)+H_invDmH(3,&
&3)*X(3,3)*X(3,4)*H_invDmH(4,2)+H_invDmH(2,3)*X(3,2)*X(3,4)*H_invD&
&mH(4,2)+H_invDmH(1,3)*X(3,1)*X(3,4)*H_invDmH(4,2)+X(2,3)*X(2,4)*H&
&_invDmH(3,3)*H_invDmH(4,2)+X(1,3)*X(1,4)*H_invDmH(3,3)*H_invDmH(4&
&,2)+X(2,2)*H_invDmH(2,3)*X(2,4)*H_invDmH(4,2)+H_invDmH(1,3)*X(2,1&
&)*X(2,4)*H_invDmH(4,2)+X(1,2)*X(1,4)*H_invDmH(2,3)*H_invDmH(4,2)+&
&X(1,1)*H_invDmH(1,3)*X(1,4)*H_invDmH(4,2)+H_invDmH(3,2)*H_invDmH(&
&3,3)*tt22+H_invDmH(2,2)*X(3,2)*H_invDmH(3,3)*X(3,3)+H_invDmH(1,2)&
&*X(3,1)*H_invDmH(3,3)*X(3,3)+H_invDmH(2,3)*H_invDmH(3,2)*X(3,2)*X&
&(3,3)+H_invDmH(1,3)*X(3,1)*H_invDmH(3,2)*X(3,3)+tt20*H_invDmH(3,2&
&)*H_invDmH(3,3)+tt19*H_invDmH(3,2)*H_invDmH(3,3)+H_invDmH(2,2)*X(&
&2,2)*X(2,3)*H_invDmH(3,3)+H_invDmH(1,2)*X(2,1)*X(2,3)*H_invDmH(3,&
&3)+X(1,2)*X(1,3)*H_invDmH(2,2)*H_invDmH(3,3)+X(1,1)*H_invDmH(1,2)&
&*X(1,3)*H_invDmH(3,3)+H_invDmH(2,2)*H_invDmH(2,3)*tt21+H_invDmH(1&
&,2)*H_invDmH(2,3)*X(3,1)*X(3,2)+H_invDmH(1,3)*H_invDmH(2,2)*X(3,1&
&)*X(3,2)+X(2,2)*H_invDmH(2,3)*X(2,3)*H_invDmH(3,2)+H_invDmH(1,3)*&
&X(2,1)*X(2,3)*H_invDmH(3,2)+X(1,2)*X(1,3)*H_invDmH(2,3)*H_invDmH(&
&3,2)+X(1,1)*H_invDmH(1,3)*X(1,3)*H_invDmH(3,2)+H_invDmH(1,2)*H_in&
&vDmH(1,3)*tt18+H_invDmH(2,2)*tt17*H_invDmH(2,3)+H_invDmH(1,2)*X(2&
&,1)*X(2,2)*H_invDmH(2,3)+tt16*H_invDmH(2,2)*H_invDmH(2,3)+X(1,1)*&
&H_invDmH(1,2)*X(1,2)*H_invDmH(2,3)+H_invDmH(1,3)*X(2,1)*H_invDmH(&
&2,2)*X(2,2)+X(1,1)*X(1,2)*H_invDmH(1,3)*H_invDmH(2,2)+H_invDmH(1,&
&2)*H_invDmH(1,3)*tt15+tt14*H_invDmH(1,2)*H_invDmH(1,3)
jac(1,1) = detDmH(1,1)*gw(1,1)*(mu(1,1)*(H_invDmH(1,3)*tt3*tt12+(&
&H_invDmH(1,2)*X(1,8)*H_invDmH(8,3)+H_invDmH(1,3)*X(1,8)*H_invDmH(&
&8,2)+H_invDmH(1,2)*X(1,7)*H_invDmH(7,3)+H_invDmH(1,3)*X(1,7)*H_in&
&vDmH(7,2)+H_invDmH(1,2)*X(1,6)*H_invDmH(6,3)+H_invDmH(1,3)*X(1,6)&
&*H_invDmH(6,2)+H_invDmH(1,2)*X(1,5)*H_invDmH(5,3)+H_invDmH(1,3)*X&
&(1,5)*H_invDmH(5,2)+H_invDmH(1,2)*X(1,4)*H_invDmH(4,3)+H_invDmH(1&
&,3)*X(1,4)*H_invDmH(4,2)+H_invDmH(1,2)*X(1,3)*H_invDmH(3,3)+H_inv&
&DmH(1,3)*X(1,3)*H_invDmH(3,2)+H_invDmH(1,2)*X(1,2)*H_invDmH(2,3)+&
&X(1,2)*H_invDmH(1,3)*H_invDmH(2,2)+2*X(1,1)*H_invDmH(1,2)*H_invDm&
&H(1,3))*tt40+(H_invDmH(1,1)*X(1,8)*H_invDmH(8,3)+H_invDmH(1,3)*X(&
&1,8)*H_invDmH(8,1)+H_invDmH(1,1)*X(1,7)*H_invDmH(7,3)+H_invDmH(1,&
&3)*X(1,7)*H_invDmH(7,1)+H_invDmH(1,1)*X(1,6)*H_invDmH(6,3)+H_invD&
&mH(1,3)*X(1,6)*H_invDmH(6,1)+H_invDmH(1,1)*X(1,5)*H_invDmH(5,3)+H&
&_invDmH(1,3)*X(1,5)*H_invDmH(5,1)+H_invDmH(1,1)*X(1,4)*H_invDmH(4&
&,3)+H_invDmH(1,3)*X(1,4)*H_invDmH(4,1)+H_invDmH(1,1)*X(1,3)*H_inv&
&DmH(3,3)+H_invDmH(1,3)*X(1,3)*H_invDmH(3,1)+H_invDmH(1,1)*X(1,2)*&
&H_invDmH(2,3)+X(1,2)*H_invDmH(1,3)*H_invDmH(2,1)+2*H_invDmH(1,1)*&
&X(1,1)*H_invDmH(1,3))*tt39+H_invDmH(1,2)*tt2*tt9+(H_invDmH(1,1)*X&
&(1,8)*H_invDmH(8,2)+H_invDmH(1,2)*X(1,8)*H_invDmH(8,1)+H_invDmH(1&
&,1)*X(1,7)*H_invDmH(7,2)+H_invDmH(1,2)*X(1,7)*H_invDmH(7,1)+H_inv&
&DmH(1,1)*X(1,6)*H_invDmH(6,2)+H_invDmH(1,2)*X(1,6)*H_invDmH(6,1)+&
&H_invDmH(1,1)*X(1,5)*H_invDmH(5,2)+H_invDmH(1,2)*X(1,5)*H_invDmH(&
&5,1)+H_invDmH(1,1)*X(1,4)*H_invDmH(4,2)+H_invDmH(1,2)*X(1,4)*H_in&
&vDmH(4,1)+H_invDmH(1,1)*X(1,3)*H_invDmH(3,2)+H_invDmH(1,2)*X(1,3)&
&*H_invDmH(3,1)+H_invDmH(1,1)*X(1,2)*H_invDmH(2,2)+H_invDmH(1,2)*X&
&(1,2)*H_invDmH(2,1)+2*H_invDmH(1,1)*X(1,1)*H_invDmH(1,2))*tt38+H_&
&invDmH(1,1)*tt1*tt6)+lam(1,1)*(H_invDmH(1,3)*tt3+H_invDmH(1,2)*tt&
&2+H_invDmH(1,1)*tt1)*tt13)
jac(1,2) = detDmH(1,1)*gw(1,1)*(mu(1,1)*(H_invDmH(1,3)*tt10*tt12+&
&(H_invDmH(1,2)*X(2,8)*H_invDmH(8,3)+H_invDmH(1,3)*X(2,8)*H_invDmH&
&(8,2)+H_invDmH(1,2)*X(2,7)*H_invDmH(7,3)+H_invDmH(1,3)*X(2,7)*H_i&
&nvDmH(7,2)+H_invDmH(1,2)*X(2,6)*H_invDmH(6,3)+H_invDmH(1,3)*X(2,6&
&)*H_invDmH(6,2)+H_invDmH(1,2)*X(2,5)*H_invDmH(5,3)+H_invDmH(1,3)*&
&X(2,5)*H_invDmH(5,2)+H_invDmH(1,2)*X(2,4)*H_invDmH(4,3)+H_invDmH(&
&1,3)*X(2,4)*H_invDmH(4,2)+H_invDmH(1,2)*X(2,3)*H_invDmH(3,3)+H_in&
&vDmH(1,3)*X(2,3)*H_invDmH(3,2)+H_invDmH(1,2)*X(2,2)*H_invDmH(2,3)&
&+H_invDmH(1,3)*H_invDmH(2,2)*X(2,2)+2*H_invDmH(1,2)*H_invDmH(1,3)&
&*X(2,1))*tt40+(H_invDmH(1,1)*X(2,8)*H_invDmH(8,3)+H_invDmH(1,3)*X&
&(2,8)*H_invDmH(8,1)+H_invDmH(1,1)*X(2,7)*H_invDmH(7,3)+H_invDmH(1&
&,3)*X(2,7)*H_invDmH(7,1)+H_invDmH(1,1)*X(2,6)*H_invDmH(6,3)+H_inv&
&DmH(1,3)*X(2,6)*H_invDmH(6,1)+H_invDmH(1,1)*X(2,5)*H_invDmH(5,3)+&
&H_invDmH(1,3)*X(2,5)*H_invDmH(5,1)+H_invDmH(1,1)*X(2,4)*H_invDmH(&
&4,3)+H_invDmH(1,3)*X(2,4)*H_invDmH(4,1)+H_invDmH(1,1)*X(2,3)*H_in&
&vDmH(3,3)+H_invDmH(1,3)*X(2,3)*H_invDmH(3,1)+H_invDmH(1,1)*X(2,2)&
&*H_invDmH(2,3)+H_invDmH(1,3)*H_invDmH(2,1)*X(2,2)+2*H_invDmH(1,1)&
&*H_invDmH(1,3)*X(2,1))*tt39+H_invDmH(1,2)*tt7*tt9+(H_invDmH(1,1)*&
&X(2,8)*H_invDmH(8,2)+H_invDmH(1,2)*X(2,8)*H_invDmH(8,1)+H_invDmH(&
&1,1)*X(2,7)*H_invDmH(7,2)+H_invDmH(1,2)*X(2,7)*H_invDmH(7,1)+H_in&
&vDmH(1,1)*X(2,6)*H_invDmH(6,2)+H_invDmH(1,2)*X(2,6)*H_invDmH(6,1)&
&+H_invDmH(1,1)*X(2,5)*H_invDmH(5,2)+H_invDmH(1,2)*X(2,5)*H_invDmH&
&(5,1)+H_invDmH(1,1)*X(2,4)*H_invDmH(4,2)+H_invDmH(1,2)*X(2,4)*H_i&
&nvDmH(4,1)+H_invDmH(1,1)*X(2,3)*H_invDmH(3,2)+H_invDmH(1,2)*X(2,3&
&)*H_invDmH(3,1)+H_invDmH(1,1)*H_invDmH(2,2)*X(2,2)+H_invDmH(1,2)*&
&H_invDmH(2,1)*X(2,2)+2*H_invDmH(1,1)*H_invDmH(1,2)*X(2,1))*tt38+H&
&_invDmH(1,1)*tt4*tt6)+lam(1,1)*(H_invDmH(1,3)*tt10+H_invDmH(1,2)*&
&tt7+H_invDmH(1,1)*tt4)*tt13)
jac(1,3) = detDmH(1,1)*gw(1,1)*(mu(1,1)*(H_invDmH(1,3)*tt11*tt12+&
&(H_invDmH(1,2)*X(3,8)*H_invDmH(8,3)+H_invDmH(1,3)*X(3,8)*H_invDmH&
&(8,2)+H_invDmH(1,2)*X(3,7)*H_invDmH(7,3)+H_invDmH(1,3)*X(3,7)*H_i&
&nvDmH(7,2)+H_invDmH(1,2)*X(3,6)*H_invDmH(6,3)+H_invDmH(1,3)*X(3,6&
&)*H_invDmH(6,2)+H_invDmH(1,2)*X(3,5)*H_invDmH(5,3)+H_invDmH(1,3)*&
&X(3,5)*H_invDmH(5,2)+H_invDmH(1,2)*X(3,4)*H_invDmH(4,3)+H_invDmH(&
&1,3)*X(3,4)*H_invDmH(4,2)+H_invDmH(1,2)*H_invDmH(3,3)*X(3,3)+H_in&
&vDmH(1,3)*H_invDmH(3,2)*X(3,3)+H_invDmH(1,2)*H_invDmH(2,3)*X(3,2)&
&+H_invDmH(1,3)*H_invDmH(2,2)*X(3,2)+2*H_invDmH(1,2)*H_invDmH(1,3)&
&*X(3,1))*tt40+(H_invDmH(1,1)*X(3,8)*H_invDmH(8,3)+H_invDmH(1,3)*X&
&(3,8)*H_invDmH(8,1)+H_invDmH(1,1)*X(3,7)*H_invDmH(7,3)+H_invDmH(1&
&,3)*X(3,7)*H_invDmH(7,1)+H_invDmH(1,1)*X(3,6)*H_invDmH(6,3)+H_inv&
&DmH(1,3)*X(3,6)*H_invDmH(6,1)+H_invDmH(1,1)*X(3,5)*H_invDmH(5,3)+&
&H_invDmH(1,3)*X(3,5)*H_invDmH(5,1)+H_invDmH(1,1)*X(3,4)*H_invDmH(&
&4,3)+H_invDmH(1,3)*X(3,4)*H_invDmH(4,1)+H_invDmH(1,1)*H_invDmH(3,&
&3)*X(3,3)+H_invDmH(1,3)*H_invDmH(3,1)*X(3,3)+H_invDmH(1,1)*H_invD&
&mH(2,3)*X(3,2)+H_invDmH(1,3)*H_invDmH(2,1)*X(3,2)+2*H_invDmH(1,1)&
&*H_invDmH(1,3)*X(3,1))*tt39+H_invDmH(1,2)*tt8*tt9+(H_invDmH(1,1)*&
&X(3,8)*H_invDmH(8,2)+H_invDmH(1,2)*X(3,8)*H_invDmH(8,1)+H_invDmH(&
&1,1)*X(3,7)*H_invDmH(7,2)+H_invDmH(1,2)*X(3,7)*H_invDmH(7,1)+H_in&
&vDmH(1,1)*X(3,6)*H_invDmH(6,2)+H_invDmH(1,2)*X(3,6)*H_invDmH(6,1)&
&+H_invDmH(1,1)*X(3,5)*H_invDmH(5,2)+H_invDmH(1,2)*X(3,5)*H_invDmH&
&(5,1)+H_invDmH(1,1)*X(3,4)*H_invDmH(4,2)+H_invDmH(1,2)*X(3,4)*H_i&
&nvDmH(4,1)+H_invDmH(1,1)*H_invDmH(3,2)*X(3,3)+H_invDmH(1,2)*H_inv&
&DmH(3,1)*X(3,3)+H_invDmH(1,1)*H_invDmH(2,2)*X(3,2)+H_invDmH(1,2)*&
&H_invDmH(2,1)*X(3,2)+2*H_invDmH(1,1)*H_invDmH(1,2)*X(3,1))*tt38+H&
&_invDmH(1,1)*tt5*tt6)+lam(1,1)*(H_invDmH(1,3)*tt11+H_invDmH(1,2)*&
&tt8+H_invDmH(1,1)*tt5)*tt13)
jac(1,4) = detDmH(1,1)*gw(1,1)*(mu(1,1)*(H_invDmH(2,3)*tt3*tt12+(&
&X(1,8)*H_invDmH(2,2)*H_invDmH(8,3)+X(1,8)*H_invDmH(2,3)*H_invDmH(&
&8,2)+X(1,7)*H_invDmH(2,2)*H_invDmH(7,3)+X(1,7)*H_invDmH(2,3)*H_in&
&vDmH(7,2)+X(1,6)*H_invDmH(2,2)*H_invDmH(6,3)+X(1,6)*H_invDmH(2,3)&
&*H_invDmH(6,2)+X(1,5)*H_invDmH(2,2)*H_invDmH(5,3)+X(1,5)*H_invDmH&
&(2,3)*H_invDmH(5,2)+X(1,4)*H_invDmH(2,2)*H_invDmH(4,3)+X(1,4)*H_i&
&nvDmH(2,3)*H_invDmH(4,2)+X(1,3)*H_invDmH(2,2)*H_invDmH(3,3)+X(1,3&
&)*H_invDmH(2,3)*H_invDmH(3,2)+2*X(1,2)*H_invDmH(2,2)*H_invDmH(2,3&
&)+X(1,1)*H_invDmH(1,2)*H_invDmH(2,3)+X(1,1)*H_invDmH(1,3)*H_invDm&
&H(2,2))*tt40+(X(1,8)*H_invDmH(2,1)*H_invDmH(8,3)+X(1,8)*H_invDmH(&
&2,3)*H_invDmH(8,1)+X(1,7)*H_invDmH(2,1)*H_invDmH(7,3)+X(1,7)*H_in&
&vDmH(2,3)*H_invDmH(7,1)+X(1,6)*H_invDmH(2,1)*H_invDmH(6,3)+X(1,6)&
&*H_invDmH(2,3)*H_invDmH(6,1)+X(1,5)*H_invDmH(2,1)*H_invDmH(5,3)+X&
&(1,5)*H_invDmH(2,3)*H_invDmH(5,1)+X(1,4)*H_invDmH(2,1)*H_invDmH(4&
&,3)+X(1,4)*H_invDmH(2,3)*H_invDmH(4,1)+X(1,3)*H_invDmH(2,1)*H_inv&
&DmH(3,3)+X(1,3)*H_invDmH(2,3)*H_invDmH(3,1)+2*X(1,2)*H_invDmH(2,1&
&)*H_invDmH(2,3)+H_invDmH(1,1)*X(1,1)*H_invDmH(2,3)+X(1,1)*H_invDm&
&H(1,3)*H_invDmH(2,1))*tt39+H_invDmH(2,2)*tt2*tt9+(X(1,8)*H_invDmH&
&(2,1)*H_invDmH(8,2)+X(1,8)*H_invDmH(2,2)*H_invDmH(8,1)+X(1,7)*H_i&
&nvDmH(2,1)*H_invDmH(7,2)+X(1,7)*H_invDmH(2,2)*H_invDmH(7,1)+X(1,6&
&)*H_invDmH(2,1)*H_invDmH(6,2)+X(1,6)*H_invDmH(2,2)*H_invDmH(6,1)+&
&X(1,5)*H_invDmH(2,1)*H_invDmH(5,2)+X(1,5)*H_invDmH(2,2)*H_invDmH(&
&5,1)+X(1,4)*H_invDmH(2,1)*H_invDmH(4,2)+X(1,4)*H_invDmH(2,2)*H_in&
&vDmH(4,1)+X(1,3)*H_invDmH(2,1)*H_invDmH(3,2)+X(1,3)*H_invDmH(2,2)&
&*H_invDmH(3,1)+2*X(1,2)*H_invDmH(2,1)*H_invDmH(2,2)+H_invDmH(1,1)&
&*X(1,1)*H_invDmH(2,2)+X(1,1)*H_invDmH(1,2)*H_invDmH(2,1))*tt38+H_&
&invDmH(2,1)*tt1*tt6)+lam(1,1)*(H_invDmH(2,3)*tt3+H_invDmH(2,2)*tt&
&2+H_invDmH(2,1)*tt1)*tt13)
jac(1,5) = detDmH(1,1)*gw(1,1)*(mu(1,1)*(H_invDmH(2,3)*tt10*tt12+&
&(H_invDmH(2,2)*X(2,8)*H_invDmH(8,3)+H_invDmH(2,3)*X(2,8)*H_invDmH&
&(8,2)+H_invDmH(2,2)*X(2,7)*H_invDmH(7,3)+H_invDmH(2,3)*X(2,7)*H_i&
&nvDmH(7,2)+H_invDmH(2,2)*X(2,6)*H_invDmH(6,3)+H_invDmH(2,3)*X(2,6&
&)*H_invDmH(6,2)+H_invDmH(2,2)*X(2,5)*H_invDmH(5,3)+H_invDmH(2,3)*&
&X(2,5)*H_invDmH(5,2)+H_invDmH(2,2)*X(2,4)*H_invDmH(4,3)+H_invDmH(&
&2,3)*X(2,4)*H_invDmH(4,2)+H_invDmH(2,2)*X(2,3)*H_invDmH(3,3)+H_in&
&vDmH(2,3)*X(2,3)*H_invDmH(3,2)+2*H_invDmH(2,2)*X(2,2)*H_invDmH(2,&
&3)+H_invDmH(1,2)*X(2,1)*H_invDmH(2,3)+H_invDmH(1,3)*X(2,1)*H_invD&
&mH(2,2))*tt40+(H_invDmH(2,1)*X(2,8)*H_invDmH(8,3)+H_invDmH(2,3)*X&
&(2,8)*H_invDmH(8,1)+H_invDmH(2,1)*X(2,7)*H_invDmH(7,3)+H_invDmH(2&
&,3)*X(2,7)*H_invDmH(7,1)+H_invDmH(2,1)*X(2,6)*H_invDmH(6,3)+H_inv&
&DmH(2,3)*X(2,6)*H_invDmH(6,1)+H_invDmH(2,1)*X(2,5)*H_invDmH(5,3)+&
&H_invDmH(2,3)*X(2,5)*H_invDmH(5,1)+H_invDmH(2,1)*X(2,4)*H_invDmH(&
&4,3)+H_invDmH(2,3)*X(2,4)*H_invDmH(4,1)+H_invDmH(2,1)*X(2,3)*H_in&
&vDmH(3,3)+H_invDmH(2,3)*X(2,3)*H_invDmH(3,1)+2*H_invDmH(2,1)*X(2,&
&2)*H_invDmH(2,3)+H_invDmH(1,1)*X(2,1)*H_invDmH(2,3)+H_invDmH(1,3)&
&*H_invDmH(2,1)*X(2,1))*tt39+H_invDmH(2,2)*tt7*tt9+(H_invDmH(2,1)*&
&X(2,8)*H_invDmH(8,2)+H_invDmH(2,2)*X(2,8)*H_invDmH(8,1)+H_invDmH(&
&2,1)*X(2,7)*H_invDmH(7,2)+H_invDmH(2,2)*X(2,7)*H_invDmH(7,1)+H_in&
&vDmH(2,1)*X(2,6)*H_invDmH(6,2)+H_invDmH(2,2)*X(2,6)*H_invDmH(6,1)&
&+H_invDmH(2,1)*X(2,5)*H_invDmH(5,2)+H_invDmH(2,2)*X(2,5)*H_invDmH&
&(5,1)+H_invDmH(2,1)*X(2,4)*H_invDmH(4,2)+H_invDmH(2,2)*X(2,4)*H_i&
&nvDmH(4,1)+H_invDmH(2,1)*X(2,3)*H_invDmH(3,2)+H_invDmH(2,2)*X(2,3&
&)*H_invDmH(3,1)+2*H_invDmH(2,1)*H_invDmH(2,2)*X(2,2)+H_invDmH(1,1&
&)*X(2,1)*H_invDmH(2,2)+H_invDmH(1,2)*H_invDmH(2,1)*X(2,1))*tt38+H&
&_invDmH(2,1)*tt4*tt6)+lam(1,1)*(H_invDmH(2,3)*tt10+H_invDmH(2,2)*&
&tt7+H_invDmH(2,1)*tt4)*tt13)
jac(1,6) = detDmH(1,1)*gw(1,1)*(mu(1,1)*(H_invDmH(2,3)*tt11*tt12+&
&(H_invDmH(2,2)*X(3,8)*H_invDmH(8,3)+H_invDmH(2,3)*X(3,8)*H_invDmH&
&(8,2)+H_invDmH(2,2)*X(3,7)*H_invDmH(7,3)+H_invDmH(2,3)*X(3,7)*H_i&
&nvDmH(7,2)+H_invDmH(2,2)*X(3,6)*H_invDmH(6,3)+H_invDmH(2,3)*X(3,6&
&)*H_invDmH(6,2)+H_invDmH(2,2)*X(3,5)*H_invDmH(5,3)+H_invDmH(2,3)*&
&X(3,5)*H_invDmH(5,2)+H_invDmH(2,2)*X(3,4)*H_invDmH(4,3)+H_invDmH(&
&2,3)*X(3,4)*H_invDmH(4,2)+H_invDmH(2,2)*H_invDmH(3,3)*X(3,3)+H_in&
&vDmH(2,3)*H_invDmH(3,2)*X(3,3)+2*H_invDmH(2,2)*H_invDmH(2,3)*X(3,&
&2)+H_invDmH(1,2)*H_invDmH(2,3)*X(3,1)+H_invDmH(1,3)*H_invDmH(2,2)&
&*X(3,1))*tt40+(H_invDmH(2,1)*X(3,8)*H_invDmH(8,3)+H_invDmH(2,3)*X&
&(3,8)*H_invDmH(8,1)+H_invDmH(2,1)*X(3,7)*H_invDmH(7,3)+H_invDmH(2&
&,3)*X(3,7)*H_invDmH(7,1)+H_invDmH(2,1)*X(3,6)*H_invDmH(6,3)+H_inv&
&DmH(2,3)*X(3,6)*H_invDmH(6,1)+H_invDmH(2,1)*X(3,5)*H_invDmH(5,3)+&
&H_invDmH(2,3)*X(3,5)*H_invDmH(5,1)+H_invDmH(2,1)*X(3,4)*H_invDmH(&
&4,3)+H_invDmH(2,3)*X(3,4)*H_invDmH(4,1)+H_invDmH(2,1)*H_invDmH(3,&
&3)*X(3,3)+H_invDmH(2,3)*H_invDmH(3,1)*X(3,3)+2*H_invDmH(2,1)*H_in&
&vDmH(2,3)*X(3,2)+H_invDmH(1,1)*H_invDmH(2,3)*X(3,1)+H_invDmH(1,3)&
&*H_invDmH(2,1)*X(3,1))*tt39+H_invDmH(2,2)*tt8*tt9+(H_invDmH(2,1)*&
&X(3,8)*H_invDmH(8,2)+H_invDmH(2,2)*X(3,8)*H_invDmH(8,1)+H_invDmH(&
&2,1)*X(3,7)*H_invDmH(7,2)+H_invDmH(2,2)*X(3,7)*H_invDmH(7,1)+H_in&
&vDmH(2,1)*X(3,6)*H_invDmH(6,2)+H_invDmH(2,2)*X(3,6)*H_invDmH(6,1)&
&+H_invDmH(2,1)*X(3,5)*H_invDmH(5,2)+H_invDmH(2,2)*X(3,5)*H_invDmH&
&(5,1)+H_invDmH(2,1)*X(3,4)*H_invDmH(4,2)+H_invDmH(2,2)*X(3,4)*H_i&
&nvDmH(4,1)+H_invDmH(2,1)*H_invDmH(3,2)*X(3,3)+H_invDmH(2,2)*H_inv&
&DmH(3,1)*X(3,3)+2*H_invDmH(2,1)*H_invDmH(2,2)*X(3,2)+H_invDmH(1,1&
&)*H_invDmH(2,2)*X(3,1)+H_invDmH(1,2)*H_invDmH(2,1)*X(3,1))*tt38+H&
&_invDmH(2,1)*tt5*tt6)+lam(1,1)*(H_invDmH(2,3)*tt11+H_invDmH(2,2)*&
&tt8+H_invDmH(2,1)*tt5)*tt13)
jac(1,7) = detDmH(1,1)*gw(1,1)*(mu(1,1)*(H_invDmH(3,3)*tt3*tt12+(&
&X(1,8)*H_invDmH(3,2)*H_invDmH(8,3)+X(1,8)*H_invDmH(3,3)*H_invDmH(&
&8,2)+X(1,7)*H_invDmH(3,2)*H_invDmH(7,3)+X(1,7)*H_invDmH(3,3)*H_in&
&vDmH(7,2)+X(1,6)*H_invDmH(3,2)*H_invDmH(6,3)+X(1,6)*H_invDmH(3,3)&
&*H_invDmH(6,2)+X(1,5)*H_invDmH(3,2)*H_invDmH(5,3)+X(1,5)*H_invDmH&
&(3,3)*H_invDmH(5,2)+X(1,4)*H_invDmH(3,2)*H_invDmH(4,3)+X(1,4)*H_i&
&nvDmH(3,3)*H_invDmH(4,2)+2*X(1,3)*H_invDmH(3,2)*H_invDmH(3,3)+X(1&
&,2)*H_invDmH(2,2)*H_invDmH(3,3)+X(1,1)*H_invDmH(1,2)*H_invDmH(3,3&
&)+X(1,2)*H_invDmH(2,3)*H_invDmH(3,2)+X(1,1)*H_invDmH(1,3)*H_invDm&
&H(3,2))*tt40+(X(1,8)*H_invDmH(3,1)*H_invDmH(8,3)+X(1,8)*H_invDmH(&
&3,3)*H_invDmH(8,1)+X(1,7)*H_invDmH(3,1)*H_invDmH(7,3)+X(1,7)*H_in&
&vDmH(3,3)*H_invDmH(7,1)+X(1,6)*H_invDmH(3,1)*H_invDmH(6,3)+X(1,6)&
&*H_invDmH(3,3)*H_invDmH(6,1)+X(1,5)*H_invDmH(3,1)*H_invDmH(5,3)+X&
&(1,5)*H_invDmH(3,3)*H_invDmH(5,1)+X(1,4)*H_invDmH(3,1)*H_invDmH(4&
&,3)+X(1,4)*H_invDmH(3,3)*H_invDmH(4,1)+2*X(1,3)*H_invDmH(3,1)*H_i&
&nvDmH(3,3)+X(1,2)*H_invDmH(2,1)*H_invDmH(3,3)+H_invDmH(1,1)*X(1,1&
&)*H_invDmH(3,3)+X(1,2)*H_invDmH(2,3)*H_invDmH(3,1)+X(1,1)*H_invDm&
&H(1,3)*H_invDmH(3,1))*tt39+H_invDmH(3,2)*tt2*tt9+(X(1,8)*H_invDmH&
&(3,1)*H_invDmH(8,2)+X(1,8)*H_invDmH(3,2)*H_invDmH(8,1)+X(1,7)*H_i&
&nvDmH(3,1)*H_invDmH(7,2)+X(1,7)*H_invDmH(3,2)*H_invDmH(7,1)+X(1,6&
&)*H_invDmH(3,1)*H_invDmH(6,2)+X(1,6)*H_invDmH(3,2)*H_invDmH(6,1)+&
&X(1,5)*H_invDmH(3,1)*H_invDmH(5,2)+X(1,5)*H_invDmH(3,2)*H_invDmH(&
&5,1)+X(1,4)*H_invDmH(3,1)*H_invDmH(4,2)+X(1,4)*H_invDmH(3,2)*H_in&
&vDmH(4,1)+2*X(1,3)*H_invDmH(3,1)*H_invDmH(3,2)+X(1,2)*H_invDmH(2,&
&1)*H_invDmH(3,2)+H_invDmH(1,1)*X(1,1)*H_invDmH(3,2)+X(1,2)*H_invD&
&mH(2,2)*H_invDmH(3,1)+X(1,1)*H_invDmH(1,2)*H_invDmH(3,1))*tt38+H_&
&invDmH(3,1)*tt1*tt6)+lam(1,1)*(H_invDmH(3,3)*tt3+H_invDmH(3,2)*tt&
&2+H_invDmH(3,1)*tt1)*tt13)
jac(1,8) = detDmH(1,1)*gw(1,1)*(mu(1,1)*(H_invDmH(3,3)*tt10*tt12+&
&(X(2,8)*H_invDmH(3,2)*H_invDmH(8,3)+X(2,8)*H_invDmH(3,3)*H_invDmH&
&(8,2)+X(2,7)*H_invDmH(3,2)*H_invDmH(7,3)+X(2,7)*H_invDmH(3,3)*H_i&
&nvDmH(7,2)+X(2,6)*H_invDmH(3,2)*H_invDmH(6,3)+X(2,6)*H_invDmH(3,3&
&)*H_invDmH(6,2)+X(2,5)*H_invDmH(3,2)*H_invDmH(5,3)+X(2,5)*H_invDm&
&H(3,3)*H_invDmH(5,2)+X(2,4)*H_invDmH(3,2)*H_invDmH(4,3)+X(2,4)*H_&
&invDmH(3,3)*H_invDmH(4,2)+2*X(2,3)*H_invDmH(3,2)*H_invDmH(3,3)+H_&
&invDmH(2,2)*X(2,2)*H_invDmH(3,3)+H_invDmH(1,2)*X(2,1)*H_invDmH(3,&
&3)+X(2,2)*H_invDmH(2,3)*H_invDmH(3,2)+H_invDmH(1,3)*X(2,1)*H_invD&
&mH(3,2))*tt40+(X(2,8)*H_invDmH(3,1)*H_invDmH(8,3)+X(2,8)*H_invDmH&
&(3,3)*H_invDmH(8,1)+X(2,7)*H_invDmH(3,1)*H_invDmH(7,3)+X(2,7)*H_i&
&nvDmH(3,3)*H_invDmH(7,1)+X(2,6)*H_invDmH(3,1)*H_invDmH(6,3)+X(2,6&
&)*H_invDmH(3,3)*H_invDmH(6,1)+X(2,5)*H_invDmH(3,1)*H_invDmH(5,3)+&
&X(2,5)*H_invDmH(3,3)*H_invDmH(5,1)+X(2,4)*H_invDmH(3,1)*H_invDmH(&
&4,3)+X(2,4)*H_invDmH(3,3)*H_invDmH(4,1)+2*X(2,3)*H_invDmH(3,1)*H_&
&invDmH(3,3)+H_invDmH(2,1)*X(2,2)*H_invDmH(3,3)+H_invDmH(1,1)*X(2,&
&1)*H_invDmH(3,3)+X(2,2)*H_invDmH(2,3)*H_invDmH(3,1)+H_invDmH(1,3)&
&*X(2,1)*H_invDmH(3,1))*tt39+H_invDmH(3,2)*tt7*tt9+(X(2,8)*H_invDm&
&H(3,1)*H_invDmH(8,2)+X(2,8)*H_invDmH(3,2)*H_invDmH(8,1)+X(2,7)*H_&
&invDmH(3,1)*H_invDmH(7,2)+X(2,7)*H_invDmH(3,2)*H_invDmH(7,1)+X(2,&
&6)*H_invDmH(3,1)*H_invDmH(6,2)+X(2,6)*H_invDmH(3,2)*H_invDmH(6,1)&
&+X(2,5)*H_invDmH(3,1)*H_invDmH(5,2)+X(2,5)*H_invDmH(3,2)*H_invDmH&
&(5,1)+X(2,4)*H_invDmH(3,1)*H_invDmH(4,2)+X(2,4)*H_invDmH(3,2)*H_i&
&nvDmH(4,1)+2*X(2,3)*H_invDmH(3,1)*H_invDmH(3,2)+H_invDmH(2,1)*X(2&
&,2)*H_invDmH(3,2)+H_invDmH(1,1)*X(2,1)*H_invDmH(3,2)+H_invDmH(2,2&
&)*X(2,2)*H_invDmH(3,1)+H_invDmH(1,2)*X(2,1)*H_invDmH(3,1))*tt38+H&
&_invDmH(3,1)*tt4*tt6)+lam(1,1)*(H_invDmH(3,3)*tt10+H_invDmH(3,2)*&
&tt7+H_invDmH(3,1)*tt4)*tt13)
jac(1,9) = detDmH(1,1)*gw(1,1)*(mu(1,1)*(H_invDmH(3,3)*tt11*tt12+&
&(H_invDmH(3,2)*X(3,8)*H_invDmH(8,3)+H_invDmH(3,3)*X(3,8)*H_invDmH&
&(8,2)+H_invDmH(3,2)*X(3,7)*H_invDmH(7,3)+H_invDmH(3,3)*X(3,7)*H_i&
&nvDmH(7,2)+H_invDmH(3,2)*X(3,6)*H_invDmH(6,3)+H_invDmH(3,3)*X(3,6&
&)*H_invDmH(6,2)+H_invDmH(3,2)*X(3,5)*H_invDmH(5,3)+H_invDmH(3,3)*&
&X(3,5)*H_invDmH(5,2)+H_invDmH(3,2)*X(3,4)*H_invDmH(4,3)+H_invDmH(&
&3,3)*X(3,4)*H_invDmH(4,2)+2*H_invDmH(3,2)*H_invDmH(3,3)*X(3,3)+H_&
&invDmH(2,2)*X(3,2)*H_invDmH(3,3)+H_invDmH(1,2)*X(3,1)*H_invDmH(3,&
&3)+H_invDmH(2,3)*H_invDmH(3,2)*X(3,2)+H_invDmH(1,3)*X(3,1)*H_invD&
&mH(3,2))*tt40+(H_invDmH(3,1)*X(3,8)*H_invDmH(8,3)+H_invDmH(3,3)*X&
&(3,8)*H_invDmH(8,1)+H_invDmH(3,1)*X(3,7)*H_invDmH(7,3)+H_invDmH(3&
&,3)*X(3,7)*H_invDmH(7,1)+H_invDmH(3,1)*X(3,6)*H_invDmH(6,3)+H_inv&
&DmH(3,3)*X(3,6)*H_invDmH(6,1)+H_invDmH(3,1)*X(3,5)*H_invDmH(5,3)+&
&H_invDmH(3,3)*X(3,5)*H_invDmH(5,1)+H_invDmH(3,1)*X(3,4)*H_invDmH(&
&4,3)+H_invDmH(3,3)*X(3,4)*H_invDmH(4,1)+2*H_invDmH(3,1)*H_invDmH(&
&3,3)*X(3,3)+H_invDmH(2,1)*X(3,2)*H_invDmH(3,3)+H_invDmH(1,1)*X(3,&
&1)*H_invDmH(3,3)+H_invDmH(2,3)*H_invDmH(3,1)*X(3,2)+H_invDmH(1,3)&
&*H_invDmH(3,1)*X(3,1))*tt39+H_invDmH(3,2)*tt8*tt9+(H_invDmH(3,1)*&
&X(3,8)*H_invDmH(8,2)+H_invDmH(3,2)*X(3,8)*H_invDmH(8,1)+H_invDmH(&
&3,1)*X(3,7)*H_invDmH(7,2)+H_invDmH(3,2)*X(3,7)*H_invDmH(7,1)+H_in&
&vDmH(3,1)*X(3,6)*H_invDmH(6,2)+H_invDmH(3,2)*X(3,6)*H_invDmH(6,1)&
&+H_invDmH(3,1)*X(3,5)*H_invDmH(5,2)+H_invDmH(3,2)*X(3,5)*H_invDmH&
&(5,1)+H_invDmH(3,1)*X(3,4)*H_invDmH(4,2)+H_invDmH(3,2)*X(3,4)*H_i&
&nvDmH(4,1)+2*H_invDmH(3,1)*H_invDmH(3,2)*X(3,3)+H_invDmH(2,1)*H_i&
&nvDmH(3,2)*X(3,2)+H_invDmH(2,2)*H_invDmH(3,1)*X(3,2)+H_invDmH(1,1&
&)*X(3,1)*H_invDmH(3,2)+H_invDmH(1,2)*H_invDmH(3,1)*X(3,1))*tt38+H&
&_invDmH(3,1)*tt5*tt6)+lam(1,1)*(H_invDmH(3,3)*tt11+H_invDmH(3,2)*&
&tt8+H_invDmH(3,1)*tt5)*tt13)
jac(1,10) = detDmH(1,1)*gw(1,1)*(mu(1,1)*(H_invDmH(4,3)*tt3*tt12+&
&(X(1,8)*H_invDmH(4,2)*H_invDmH(8,3)+X(1,8)*H_invDmH(4,3)*H_invDmH&
&(8,2)+X(1,7)*H_invDmH(4,2)*H_invDmH(7,3)+X(1,7)*H_invDmH(4,3)*H_i&
&nvDmH(7,2)+X(1,6)*H_invDmH(4,2)*H_invDmH(6,3)+X(1,6)*H_invDmH(4,3&
&)*H_invDmH(6,2)+X(1,5)*H_invDmH(4,2)*H_invDmH(5,3)+X(1,5)*H_invDm&
&H(4,3)*H_invDmH(5,2)+2*X(1,4)*H_invDmH(4,2)*H_invDmH(4,3)+X(1,3)*&
&H_invDmH(3,2)*H_invDmH(4,3)+X(1,2)*H_invDmH(2,2)*H_invDmH(4,3)+X(&
&1,1)*H_invDmH(1,2)*H_invDmH(4,3)+X(1,3)*H_invDmH(3,3)*H_invDmH(4,&
&2)+X(1,2)*H_invDmH(2,3)*H_invDmH(4,2)+X(1,1)*H_invDmH(1,3)*H_invD&
&mH(4,2))*tt40+(X(1,8)*H_invDmH(4,1)*H_invDmH(8,3)+X(1,8)*H_invDmH&
&(4,3)*H_invDmH(8,1)+X(1,7)*H_invDmH(4,1)*H_invDmH(7,3)+X(1,7)*H_i&
&nvDmH(4,3)*H_invDmH(7,1)+X(1,6)*H_invDmH(4,1)*H_invDmH(6,3)+X(1,6&
&)*H_invDmH(4,3)*H_invDmH(6,1)+X(1,5)*H_invDmH(4,1)*H_invDmH(5,3)+&
&X(1,5)*H_invDmH(4,3)*H_invDmH(5,1)+2*X(1,4)*H_invDmH(4,1)*H_invDm&
&H(4,3)+X(1,3)*H_invDmH(3,1)*H_invDmH(4,3)+X(1,2)*H_invDmH(2,1)*H_&
&invDmH(4,3)+H_invDmH(1,1)*X(1,1)*H_invDmH(4,3)+X(1,3)*H_invDmH(3,&
&3)*H_invDmH(4,1)+X(1,2)*H_invDmH(2,3)*H_invDmH(4,1)+X(1,1)*H_invD&
&mH(1,3)*H_invDmH(4,1))*tt39+H_invDmH(4,2)*tt2*tt9+(X(1,8)*H_invDm&
&H(4,1)*H_invDmH(8,2)+X(1,8)*H_invDmH(4,2)*H_invDmH(8,1)+X(1,7)*H_&
&invDmH(4,1)*H_invDmH(7,2)+X(1,7)*H_invDmH(4,2)*H_invDmH(7,1)+X(1,&
&6)*H_invDmH(4,1)*H_invDmH(6,2)+X(1,6)*H_invDmH(4,2)*H_invDmH(6,1)&
&+X(1,5)*H_invDmH(4,1)*H_invDmH(5,2)+X(1,5)*H_invDmH(4,2)*H_invDmH&
&(5,1)+2*X(1,4)*H_invDmH(4,1)*H_invDmH(4,2)+X(1,3)*H_invDmH(3,1)*H&
&_invDmH(4,2)+X(1,2)*H_invDmH(2,1)*H_invDmH(4,2)+H_invDmH(1,1)*X(1&
&,1)*H_invDmH(4,2)+X(1,3)*H_invDmH(3,2)*H_invDmH(4,1)+X(1,2)*H_inv&
&DmH(2,2)*H_invDmH(4,1)+X(1,1)*H_invDmH(1,2)*H_invDmH(4,1))*tt38+H&
&_invDmH(4,1)*tt1*tt6)+lam(1,1)*(H_invDmH(4,3)*tt3+H_invDmH(4,2)*t&
&t2+H_invDmH(4,1)*tt1)*tt13)
jac(1,11) = detDmH(1,1)*gw(1,1)*(mu(1,1)*(H_invDmH(4,3)*tt10*tt12&
&+(X(2,8)*H_invDmH(4,2)*H_invDmH(8,3)+X(2,8)*H_invDmH(4,3)*H_invDm&
&H(8,2)+X(2,7)*H_invDmH(4,2)*H_invDmH(7,3)+X(2,7)*H_invDmH(4,3)*H_&
&invDmH(7,2)+X(2,6)*H_invDmH(4,2)*H_invDmH(6,3)+X(2,6)*H_invDmH(4,&
&3)*H_invDmH(6,2)+X(2,5)*H_invDmH(4,2)*H_invDmH(5,3)+X(2,5)*H_invD&
&mH(4,3)*H_invDmH(5,2)+2*X(2,4)*H_invDmH(4,2)*H_invDmH(4,3)+X(2,3)&
&*H_invDmH(3,2)*H_invDmH(4,3)+H_invDmH(2,2)*X(2,2)*H_invDmH(4,3)+H&
&_invDmH(1,2)*X(2,1)*H_invDmH(4,3)+X(2,3)*H_invDmH(3,3)*H_invDmH(4&
&,2)+X(2,2)*H_invDmH(2,3)*H_invDmH(4,2)+H_invDmH(1,3)*X(2,1)*H_inv&
&DmH(4,2))*tt40+(X(2,8)*H_invDmH(4,1)*H_invDmH(8,3)+X(2,8)*H_invDm&
&H(4,3)*H_invDmH(8,1)+X(2,7)*H_invDmH(4,1)*H_invDmH(7,3)+X(2,7)*H_&
&invDmH(4,3)*H_invDmH(7,1)+X(2,6)*H_invDmH(4,1)*H_invDmH(6,3)+X(2,&
&6)*H_invDmH(4,3)*H_invDmH(6,1)+X(2,5)*H_invDmH(4,1)*H_invDmH(5,3)&
&+X(2,5)*H_invDmH(4,3)*H_invDmH(5,1)+2*X(2,4)*H_invDmH(4,1)*H_invD&
&mH(4,3)+X(2,3)*H_invDmH(3,1)*H_invDmH(4,3)+H_invDmH(2,1)*X(2,2)*H&
&_invDmH(4,3)+H_invDmH(1,1)*X(2,1)*H_invDmH(4,3)+X(2,3)*H_invDmH(3&
&,3)*H_invDmH(4,1)+X(2,2)*H_invDmH(2,3)*H_invDmH(4,1)+H_invDmH(1,3&
&)*X(2,1)*H_invDmH(4,1))*tt39+H_invDmH(4,2)*tt7*tt9+(X(2,8)*H_invD&
&mH(4,1)*H_invDmH(8,2)+X(2,8)*H_invDmH(4,2)*H_invDmH(8,1)+X(2,7)*H&
&_invDmH(4,1)*H_invDmH(7,2)+X(2,7)*H_invDmH(4,2)*H_invDmH(7,1)+X(2&
&,6)*H_invDmH(4,1)*H_invDmH(6,2)+X(2,6)*H_invDmH(4,2)*H_invDmH(6,1&
&)+X(2,5)*H_invDmH(4,1)*H_invDmH(5,2)+X(2,5)*H_invDmH(4,2)*H_invDm&
&H(5,1)+2*X(2,4)*H_invDmH(4,1)*H_invDmH(4,2)+X(2,3)*H_invDmH(3,1)*&
&H_invDmH(4,2)+H_invDmH(2,1)*X(2,2)*H_invDmH(4,2)+H_invDmH(1,1)*X(&
&2,1)*H_invDmH(4,2)+X(2,3)*H_invDmH(3,2)*H_invDmH(4,1)+H_invDmH(2,&
&2)*X(2,2)*H_invDmH(4,1)+H_invDmH(1,2)*X(2,1)*H_invDmH(4,1))*tt38+&
&H_invDmH(4,1)*tt4*tt6)+lam(1,1)*(H_invDmH(4,3)*tt10+H_invDmH(4,2)&
&*tt7+H_invDmH(4,1)*tt4)*tt13)
jac(1,12) = detDmH(1,1)*gw(1,1)*(mu(1,1)*(H_invDmH(4,3)*tt11*tt12&
&+(X(3,8)*H_invDmH(4,2)*H_invDmH(8,3)+X(3,8)*H_invDmH(4,3)*H_invDm&
&H(8,2)+X(3,7)*H_invDmH(4,2)*H_invDmH(7,3)+X(3,7)*H_invDmH(4,3)*H_&
&invDmH(7,2)+X(3,6)*H_invDmH(4,2)*H_invDmH(6,3)+X(3,6)*H_invDmH(4,&
&3)*H_invDmH(6,2)+X(3,5)*H_invDmH(4,2)*H_invDmH(5,3)+X(3,5)*H_invD&
&mH(4,3)*H_invDmH(5,2)+2*X(3,4)*H_invDmH(4,2)*H_invDmH(4,3)+H_invD&
&mH(3,2)*X(3,3)*H_invDmH(4,3)+H_invDmH(2,2)*X(3,2)*H_invDmH(4,3)+H&
&_invDmH(1,2)*X(3,1)*H_invDmH(4,3)+H_invDmH(3,3)*X(3,3)*H_invDmH(4&
&,2)+H_invDmH(2,3)*X(3,2)*H_invDmH(4,2)+H_invDmH(1,3)*X(3,1)*H_inv&
&DmH(4,2))*tt40+(X(3,8)*H_invDmH(4,1)*H_invDmH(8,3)+X(3,8)*H_invDm&
&H(4,3)*H_invDmH(8,1)+X(3,7)*H_invDmH(4,1)*H_invDmH(7,3)+X(3,7)*H_&
&invDmH(4,3)*H_invDmH(7,1)+X(3,6)*H_invDmH(4,1)*H_invDmH(6,3)+X(3,&
&6)*H_invDmH(4,3)*H_invDmH(6,1)+X(3,5)*H_invDmH(4,1)*H_invDmH(5,3)&
&+X(3,5)*H_invDmH(4,3)*H_invDmH(5,1)+2*X(3,4)*H_invDmH(4,1)*H_invD&
&mH(4,3)+H_invDmH(3,1)*X(3,3)*H_invDmH(4,3)+H_invDmH(2,1)*X(3,2)*H&
&_invDmH(4,3)+H_invDmH(1,1)*X(3,1)*H_invDmH(4,3)+H_invDmH(3,3)*X(3&
&,3)*H_invDmH(4,1)+H_invDmH(2,3)*X(3,2)*H_invDmH(4,1)+H_invDmH(1,3&
&)*X(3,1)*H_invDmH(4,1))*tt39+H_invDmH(4,2)*tt8*tt9+(X(3,8)*H_invD&
&mH(4,1)*H_invDmH(8,2)+X(3,8)*H_invDmH(4,2)*H_invDmH(8,1)+X(3,7)*H&
&_invDmH(4,1)*H_invDmH(7,2)+X(3,7)*H_invDmH(4,2)*H_invDmH(7,1)+X(3&
&,6)*H_invDmH(4,1)*H_invDmH(6,2)+X(3,6)*H_invDmH(4,2)*H_invDmH(6,1&
&)+X(3,5)*H_invDmH(4,1)*H_invDmH(5,2)+X(3,5)*H_invDmH(4,2)*H_invDm&
&H(5,1)+2*X(3,4)*H_invDmH(4,1)*H_invDmH(4,2)+H_invDmH(3,1)*X(3,3)*&
&H_invDmH(4,2)+H_invDmH(2,1)*X(3,2)*H_invDmH(4,2)+H_invDmH(1,1)*X(&
&3,1)*H_invDmH(4,2)+H_invDmH(3,2)*X(3,3)*H_invDmH(4,1)+H_invDmH(2,&
&2)*X(3,2)*H_invDmH(4,1)+H_invDmH(1,2)*X(3,1)*H_invDmH(4,1))*tt38+&
&H_invDmH(4,1)*tt5*tt6)+lam(1,1)*(H_invDmH(4,3)*tt11+H_invDmH(4,2)&
&*tt8+H_invDmH(4,1)*tt5)*tt13)
jac(1,13) = detDmH(1,1)*gw(1,1)*(mu(1,1)*(H_invDmH(5,3)*tt3*tt12+&
&(X(1,8)*H_invDmH(5,2)*H_invDmH(8,3)+X(1,8)*H_invDmH(5,3)*H_invDmH&
&(8,2)+X(1,7)*H_invDmH(5,2)*H_invDmH(7,3)+X(1,7)*H_invDmH(5,3)*H_i&
&nvDmH(7,2)+X(1,6)*H_invDmH(5,2)*H_invDmH(6,3)+X(1,6)*H_invDmH(5,3&
&)*H_invDmH(6,2)+2*X(1,5)*H_invDmH(5,2)*H_invDmH(5,3)+X(1,4)*H_inv&
&DmH(4,2)*H_invDmH(5,3)+X(1,3)*H_invDmH(3,2)*H_invDmH(5,3)+X(1,2)*&
&H_invDmH(2,2)*H_invDmH(5,3)+X(1,1)*H_invDmH(1,2)*H_invDmH(5,3)+X(&
&1,4)*H_invDmH(4,3)*H_invDmH(5,2)+X(1,3)*H_invDmH(3,3)*H_invDmH(5,&
&2)+X(1,2)*H_invDmH(2,3)*H_invDmH(5,2)+X(1,1)*H_invDmH(1,3)*H_invD&
&mH(5,2))*tt40+(X(1,8)*H_invDmH(5,1)*H_invDmH(8,3)+X(1,8)*H_invDmH&
&(5,3)*H_invDmH(8,1)+X(1,7)*H_invDmH(5,1)*H_invDmH(7,3)+X(1,7)*H_i&
&nvDmH(5,3)*H_invDmH(7,1)+X(1,6)*H_invDmH(5,1)*H_invDmH(6,3)+X(1,6&
&)*H_invDmH(5,3)*H_invDmH(6,1)+2*X(1,5)*H_invDmH(5,1)*H_invDmH(5,3&
&)+X(1,4)*H_invDmH(4,1)*H_invDmH(5,3)+X(1,3)*H_invDmH(3,1)*H_invDm&
&H(5,3)+X(1,2)*H_invDmH(2,1)*H_invDmH(5,3)+H_invDmH(1,1)*X(1,1)*H_&
&invDmH(5,3)+X(1,4)*H_invDmH(4,3)*H_invDmH(5,1)+X(1,3)*H_invDmH(3,&
&3)*H_invDmH(5,1)+X(1,2)*H_invDmH(2,3)*H_invDmH(5,1)+X(1,1)*H_invD&
&mH(1,3)*H_invDmH(5,1))*tt39+H_invDmH(5,2)*tt2*tt9+(X(1,8)*H_invDm&
&H(5,1)*H_invDmH(8,2)+X(1,8)*H_invDmH(5,2)*H_invDmH(8,1)+X(1,7)*H_&
&invDmH(5,1)*H_invDmH(7,2)+X(1,7)*H_invDmH(5,2)*H_invDmH(7,1)+X(1,&
&6)*H_invDmH(5,1)*H_invDmH(6,2)+X(1,6)*H_invDmH(5,2)*H_invDmH(6,1)&
&+2*X(1,5)*H_invDmH(5,1)*H_invDmH(5,2)+X(1,4)*H_invDmH(4,1)*H_invD&
&mH(5,2)+X(1,3)*H_invDmH(3,1)*H_invDmH(5,2)+X(1,2)*H_invDmH(2,1)*H&
&_invDmH(5,2)+H_invDmH(1,1)*X(1,1)*H_invDmH(5,2)+X(1,4)*H_invDmH(4&
&,2)*H_invDmH(5,1)+X(1,3)*H_invDmH(3,2)*H_invDmH(5,1)+X(1,2)*H_inv&
&DmH(2,2)*H_invDmH(5,1)+X(1,1)*H_invDmH(1,2)*H_invDmH(5,1))*tt38+H&
&_invDmH(5,1)*tt1*tt6)+lam(1,1)*(H_invDmH(5,3)*tt3+H_invDmH(5,2)*t&
&t2+H_invDmH(5,1)*tt1)*tt13)
jac(1,14) = detDmH(1,1)*gw(1,1)*(mu(1,1)*(H_invDmH(5,3)*tt10*tt12&
&+(X(2,8)*H_invDmH(5,2)*H_invDmH(8,3)+X(2,8)*H_invDmH(5,3)*H_invDm&
&H(8,2)+X(2,7)*H_invDmH(5,2)*H_invDmH(7,3)+X(2,7)*H_invDmH(5,3)*H_&
&invDmH(7,2)+X(2,6)*H_invDmH(5,2)*H_invDmH(6,3)+X(2,6)*H_invDmH(5,&
&3)*H_invDmH(6,2)+2*X(2,5)*H_invDmH(5,2)*H_invDmH(5,3)+X(2,4)*H_in&
&vDmH(4,2)*H_invDmH(5,3)+X(2,3)*H_invDmH(3,2)*H_invDmH(5,3)+H_invD&
&mH(2,2)*X(2,2)*H_invDmH(5,3)+H_invDmH(1,2)*X(2,1)*H_invDmH(5,3)+X&
&(2,4)*H_invDmH(4,3)*H_invDmH(5,2)+X(2,3)*H_invDmH(3,3)*H_invDmH(5&
&,2)+X(2,2)*H_invDmH(2,3)*H_invDmH(5,2)+H_invDmH(1,3)*X(2,1)*H_inv&
&DmH(5,2))*tt40+(X(2,8)*H_invDmH(5,1)*H_invDmH(8,3)+X(2,8)*H_invDm&
&H(5,3)*H_invDmH(8,1)+X(2,7)*H_invDmH(5,1)*H_invDmH(7,3)+X(2,7)*H_&
&invDmH(5,3)*H_invDmH(7,1)+X(2,6)*H_invDmH(5,1)*H_invDmH(6,3)+X(2,&
&6)*H_invDmH(5,3)*H_invDmH(6,1)+2*X(2,5)*H_invDmH(5,1)*H_invDmH(5,&
&3)+X(2,4)*H_invDmH(4,1)*H_invDmH(5,3)+X(2,3)*H_invDmH(3,1)*H_invD&
&mH(5,3)+H_invDmH(2,1)*X(2,2)*H_invDmH(5,3)+H_invDmH(1,1)*X(2,1)*H&
&_invDmH(5,3)+X(2,4)*H_invDmH(4,3)*H_invDmH(5,1)+X(2,3)*H_invDmH(3&
&,3)*H_invDmH(5,1)+X(2,2)*H_invDmH(2,3)*H_invDmH(5,1)+H_invDmH(1,3&
&)*X(2,1)*H_invDmH(5,1))*tt39+H_invDmH(5,2)*tt7*tt9+(X(2,8)*H_invD&
&mH(5,1)*H_invDmH(8,2)+X(2,8)*H_invDmH(5,2)*H_invDmH(8,1)+X(2,7)*H&
&_invDmH(5,1)*H_invDmH(7,2)+X(2,7)*H_invDmH(5,2)*H_invDmH(7,1)+X(2&
&,6)*H_invDmH(5,1)*H_invDmH(6,2)+X(2,6)*H_invDmH(5,2)*H_invDmH(6,1&
&)+2*X(2,5)*H_invDmH(5,1)*H_invDmH(5,2)+X(2,4)*H_invDmH(4,1)*H_inv&
&DmH(5,2)+X(2,3)*H_invDmH(3,1)*H_invDmH(5,2)+H_invDmH(2,1)*X(2,2)*&
&H_invDmH(5,2)+H_invDmH(1,1)*X(2,1)*H_invDmH(5,2)+X(2,4)*H_invDmH(&
&4,2)*H_invDmH(5,1)+X(2,3)*H_invDmH(3,2)*H_invDmH(5,1)+H_invDmH(2,&
&2)*X(2,2)*H_invDmH(5,1)+H_invDmH(1,2)*X(2,1)*H_invDmH(5,1))*tt38+&
&H_invDmH(5,1)*tt4*tt6)+lam(1,1)*(H_invDmH(5,3)*tt10+H_invDmH(5,2)&
&*tt7+H_invDmH(5,1)*tt4)*tt13)
jac(1,15) = detDmH(1,1)*gw(1,1)*(mu(1,1)*(H_invDmH(5,3)*tt11*tt12&
&+(X(3,8)*H_invDmH(5,2)*H_invDmH(8,3)+X(3,8)*H_invDmH(5,3)*H_invDm&
&H(8,2)+X(3,7)*H_invDmH(5,2)*H_invDmH(7,3)+X(3,7)*H_invDmH(5,3)*H_&
&invDmH(7,2)+X(3,6)*H_invDmH(5,2)*H_invDmH(6,3)+X(3,6)*H_invDmH(5,&
&3)*H_invDmH(6,2)+2*X(3,5)*H_invDmH(5,2)*H_invDmH(5,3)+X(3,4)*H_in&
&vDmH(4,2)*H_invDmH(5,3)+H_invDmH(3,2)*X(3,3)*H_invDmH(5,3)+H_invD&
&mH(2,2)*X(3,2)*H_invDmH(5,3)+H_invDmH(1,2)*X(3,1)*H_invDmH(5,3)+X&
&(3,4)*H_invDmH(4,3)*H_invDmH(5,2)+H_invDmH(3,3)*X(3,3)*H_invDmH(5&
&,2)+H_invDmH(2,3)*X(3,2)*H_invDmH(5,2)+H_invDmH(1,3)*X(3,1)*H_inv&
&DmH(5,2))*tt40+(X(3,8)*H_invDmH(5,1)*H_invDmH(8,3)+X(3,8)*H_invDm&
&H(5,3)*H_invDmH(8,1)+X(3,7)*H_invDmH(5,1)*H_invDmH(7,3)+X(3,7)*H_&
&invDmH(5,3)*H_invDmH(7,1)+X(3,6)*H_invDmH(5,1)*H_invDmH(6,3)+X(3,&
&6)*H_invDmH(5,3)*H_invDmH(6,1)+2*X(3,5)*H_invDmH(5,1)*H_invDmH(5,&
&3)+X(3,4)*H_invDmH(4,1)*H_invDmH(5,3)+H_invDmH(3,1)*X(3,3)*H_invD&
&mH(5,3)+H_invDmH(2,1)*X(3,2)*H_invDmH(5,3)+H_invDmH(1,1)*X(3,1)*H&
&_invDmH(5,3)+X(3,4)*H_invDmH(4,3)*H_invDmH(5,1)+H_invDmH(3,3)*X(3&
&,3)*H_invDmH(5,1)+H_invDmH(2,3)*X(3,2)*H_invDmH(5,1)+H_invDmH(1,3&
&)*X(3,1)*H_invDmH(5,1))*tt39+H_invDmH(5,2)*tt8*tt9+(X(3,8)*H_invD&
&mH(5,1)*H_invDmH(8,2)+X(3,8)*H_invDmH(5,2)*H_invDmH(8,1)+X(3,7)*H&
&_invDmH(5,1)*H_invDmH(7,2)+X(3,7)*H_invDmH(5,2)*H_invDmH(7,1)+X(3&
&,6)*H_invDmH(5,1)*H_invDmH(6,2)+X(3,6)*H_invDmH(5,2)*H_invDmH(6,1&
&)+2*X(3,5)*H_invDmH(5,1)*H_invDmH(5,2)+X(3,4)*H_invDmH(4,1)*H_inv&
&DmH(5,2)+H_invDmH(3,1)*X(3,3)*H_invDmH(5,2)+H_invDmH(2,1)*X(3,2)*&
&H_invDmH(5,2)+H_invDmH(1,1)*X(3,1)*H_invDmH(5,2)+X(3,4)*H_invDmH(&
&4,2)*H_invDmH(5,1)+H_invDmH(3,2)*X(3,3)*H_invDmH(5,1)+H_invDmH(2,&
&2)*X(3,2)*H_invDmH(5,1)+H_invDmH(1,2)*X(3,1)*H_invDmH(5,1))*tt38+&
&H_invDmH(5,1)*tt5*tt6)+lam(1,1)*(H_invDmH(5,3)*tt11+H_invDmH(5,2)&
&*tt8+H_invDmH(5,1)*tt5)*tt13)
jac(1,16) = detDmH(1,1)*gw(1,1)*(mu(1,1)*(H_invDmH(6,3)*tt3*tt12+&
&(X(1,8)*H_invDmH(6,2)*H_invDmH(8,3)+X(1,8)*H_invDmH(6,3)*H_invDmH&
&(8,2)+X(1,7)*H_invDmH(6,2)*H_invDmH(7,3)+X(1,7)*H_invDmH(6,3)*H_i&
&nvDmH(7,2)+2*X(1,6)*H_invDmH(6,2)*H_invDmH(6,3)+X(1,5)*H_invDmH(5&
&,2)*H_invDmH(6,3)+X(1,4)*H_invDmH(4,2)*H_invDmH(6,3)+X(1,3)*H_inv&
&DmH(3,2)*H_invDmH(6,3)+X(1,2)*H_invDmH(2,2)*H_invDmH(6,3)+X(1,1)*&
&H_invDmH(1,2)*H_invDmH(6,3)+X(1,5)*H_invDmH(5,3)*H_invDmH(6,2)+X(&
&1,4)*H_invDmH(4,3)*H_invDmH(6,2)+X(1,3)*H_invDmH(3,3)*H_invDmH(6,&
&2)+X(1,2)*H_invDmH(2,3)*H_invDmH(6,2)+X(1,1)*H_invDmH(1,3)*H_invD&
&mH(6,2))*tt40+(X(1,8)*H_invDmH(6,1)*H_invDmH(8,3)+X(1,8)*H_invDmH&
&(6,3)*H_invDmH(8,1)+X(1,7)*H_invDmH(6,1)*H_invDmH(7,3)+X(1,7)*H_i&
&nvDmH(6,3)*H_invDmH(7,1)+2*X(1,6)*H_invDmH(6,1)*H_invDmH(6,3)+X(1&
&,5)*H_invDmH(5,1)*H_invDmH(6,3)+X(1,4)*H_invDmH(4,1)*H_invDmH(6,3&
&)+X(1,3)*H_invDmH(3,1)*H_invDmH(6,3)+X(1,2)*H_invDmH(2,1)*H_invDm&
&H(6,3)+H_invDmH(1,1)*X(1,1)*H_invDmH(6,3)+X(1,5)*H_invDmH(5,3)*H_&
&invDmH(6,1)+X(1,4)*H_invDmH(4,3)*H_invDmH(6,1)+X(1,3)*H_invDmH(3,&
&3)*H_invDmH(6,1)+X(1,2)*H_invDmH(2,3)*H_invDmH(6,1)+X(1,1)*H_invD&
&mH(1,3)*H_invDmH(6,1))*tt39+H_invDmH(6,2)*tt2*tt9+(X(1,8)*H_invDm&
&H(6,1)*H_invDmH(8,2)+X(1,8)*H_invDmH(6,2)*H_invDmH(8,1)+X(1,7)*H_&
&invDmH(6,1)*H_invDmH(7,2)+X(1,7)*H_invDmH(6,2)*H_invDmH(7,1)+2*X(&
&1,6)*H_invDmH(6,1)*H_invDmH(6,2)+X(1,5)*H_invDmH(5,1)*H_invDmH(6,&
&2)+X(1,4)*H_invDmH(4,1)*H_invDmH(6,2)+X(1,3)*H_invDmH(3,1)*H_invD&
&mH(6,2)+X(1,2)*H_invDmH(2,1)*H_invDmH(6,2)+H_invDmH(1,1)*X(1,1)*H&
&_invDmH(6,2)+X(1,5)*H_invDmH(5,2)*H_invDmH(6,1)+X(1,4)*H_invDmH(4&
&,2)*H_invDmH(6,1)+X(1,3)*H_invDmH(3,2)*H_invDmH(6,1)+X(1,2)*H_inv&
&DmH(2,2)*H_invDmH(6,1)+X(1,1)*H_invDmH(1,2)*H_invDmH(6,1))*tt38+H&
&_invDmH(6,1)*tt1*tt6)+lam(1,1)*(H_invDmH(6,3)*tt3+H_invDmH(6,2)*t&
&t2+H_invDmH(6,1)*tt1)*tt13)
jac(1,17) = detDmH(1,1)*gw(1,1)*(mu(1,1)*(H_invDmH(6,3)*tt10*tt12&
&+(X(2,8)*H_invDmH(6,2)*H_invDmH(8,3)+X(2,8)*H_invDmH(6,3)*H_invDm&
&H(8,2)+X(2,7)*H_invDmH(6,2)*H_invDmH(7,3)+X(2,7)*H_invDmH(6,3)*H_&
&invDmH(7,2)+2*X(2,6)*H_invDmH(6,2)*H_invDmH(6,3)+X(2,5)*H_invDmH(&
&5,2)*H_invDmH(6,3)+X(2,4)*H_invDmH(4,2)*H_invDmH(6,3)+X(2,3)*H_in&
&vDmH(3,2)*H_invDmH(6,3)+H_invDmH(2,2)*X(2,2)*H_invDmH(6,3)+H_invD&
&mH(1,2)*X(2,1)*H_invDmH(6,3)+X(2,5)*H_invDmH(5,3)*H_invDmH(6,2)+X&
&(2,4)*H_invDmH(4,3)*H_invDmH(6,2)+X(2,3)*H_invDmH(3,3)*H_invDmH(6&
&,2)+X(2,2)*H_invDmH(2,3)*H_invDmH(6,2)+H_invDmH(1,3)*X(2,1)*H_inv&
&DmH(6,2))*tt40+(X(2,8)*H_invDmH(6,1)*H_invDmH(8,3)+X(2,8)*H_invDm&
&H(6,3)*H_invDmH(8,1)+X(2,7)*H_invDmH(6,1)*H_invDmH(7,3)+X(2,7)*H_&
&invDmH(6,3)*H_invDmH(7,1)+2*X(2,6)*H_invDmH(6,1)*H_invDmH(6,3)+X(&
&2,5)*H_invDmH(5,1)*H_invDmH(6,3)+X(2,4)*H_invDmH(4,1)*H_invDmH(6,&
&3)+X(2,3)*H_invDmH(3,1)*H_invDmH(6,3)+H_invDmH(2,1)*X(2,2)*H_invD&
&mH(6,3)+H_invDmH(1,1)*X(2,1)*H_invDmH(6,3)+X(2,5)*H_invDmH(5,3)*H&
&_invDmH(6,1)+X(2,4)*H_invDmH(4,3)*H_invDmH(6,1)+X(2,3)*H_invDmH(3&
&,3)*H_invDmH(6,1)+X(2,2)*H_invDmH(2,3)*H_invDmH(6,1)+H_invDmH(1,3&
&)*X(2,1)*H_invDmH(6,1))*tt39+H_invDmH(6,2)*tt7*tt9+(X(2,8)*H_invD&
&mH(6,1)*H_invDmH(8,2)+X(2,8)*H_invDmH(6,2)*H_invDmH(8,1)+X(2,7)*H&
&_invDmH(6,1)*H_invDmH(7,2)+X(2,7)*H_invDmH(6,2)*H_invDmH(7,1)+2*X&
&(2,6)*H_invDmH(6,1)*H_invDmH(6,2)+X(2,5)*H_invDmH(5,1)*H_invDmH(6&
&,2)+X(2,4)*H_invDmH(4,1)*H_invDmH(6,2)+X(2,3)*H_invDmH(3,1)*H_inv&
&DmH(6,2)+H_invDmH(2,1)*X(2,2)*H_invDmH(6,2)+H_invDmH(1,1)*X(2,1)*&
&H_invDmH(6,2)+X(2,5)*H_invDmH(5,2)*H_invDmH(6,1)+X(2,4)*H_invDmH(&
&4,2)*H_invDmH(6,1)+X(2,3)*H_invDmH(3,2)*H_invDmH(6,1)+H_invDmH(2,&
&2)*X(2,2)*H_invDmH(6,1)+H_invDmH(1,2)*X(2,1)*H_invDmH(6,1))*tt38+&
&H_invDmH(6,1)*tt4*tt6)+lam(1,1)*(H_invDmH(6,3)*tt10+H_invDmH(6,2)&
&*tt7+H_invDmH(6,1)*tt4)*tt13)
jac(1,18) = detDmH(1,1)*gw(1,1)*(mu(1,1)*(H_invDmH(6,3)*tt11*tt12&
&+(X(3,8)*H_invDmH(6,2)*H_invDmH(8,3)+X(3,8)*H_invDmH(6,3)*H_invDm&
&H(8,2)+X(3,7)*H_invDmH(6,2)*H_invDmH(7,3)+X(3,7)*H_invDmH(6,3)*H_&
&invDmH(7,2)+2*X(3,6)*H_invDmH(6,2)*H_invDmH(6,3)+X(3,5)*H_invDmH(&
&5,2)*H_invDmH(6,3)+X(3,4)*H_invDmH(4,2)*H_invDmH(6,3)+H_invDmH(3,&
&2)*X(3,3)*H_invDmH(6,3)+H_invDmH(2,2)*X(3,2)*H_invDmH(6,3)+H_invD&
&mH(1,2)*X(3,1)*H_invDmH(6,3)+X(3,5)*H_invDmH(5,3)*H_invDmH(6,2)+X&
&(3,4)*H_invDmH(4,3)*H_invDmH(6,2)+H_invDmH(3,3)*X(3,3)*H_invDmH(6&
&,2)+H_invDmH(2,3)*X(3,2)*H_invDmH(6,2)+H_invDmH(1,3)*X(3,1)*H_inv&
&DmH(6,2))*tt40+(X(3,8)*H_invDmH(6,1)*H_invDmH(8,3)+X(3,8)*H_invDm&
&H(6,3)*H_invDmH(8,1)+X(3,7)*H_invDmH(6,1)*H_invDmH(7,3)+X(3,7)*H_&
&invDmH(6,3)*H_invDmH(7,1)+2*X(3,6)*H_invDmH(6,1)*H_invDmH(6,3)+X(&
&3,5)*H_invDmH(5,1)*H_invDmH(6,3)+X(3,4)*H_invDmH(4,1)*H_invDmH(6,&
&3)+H_invDmH(3,1)*X(3,3)*H_invDmH(6,3)+H_invDmH(2,1)*X(3,2)*H_invD&
&mH(6,3)+H_invDmH(1,1)*X(3,1)*H_invDmH(6,3)+X(3,5)*H_invDmH(5,3)*H&
&_invDmH(6,1)+X(3,4)*H_invDmH(4,3)*H_invDmH(6,1)+H_invDmH(3,3)*X(3&
&,3)*H_invDmH(6,1)+H_invDmH(2,3)*X(3,2)*H_invDmH(6,1)+H_invDmH(1,3&
&)*X(3,1)*H_invDmH(6,1))*tt39+H_invDmH(6,2)*tt8*tt9+(X(3,8)*H_invD&
&mH(6,1)*H_invDmH(8,2)+X(3,8)*H_invDmH(6,2)*H_invDmH(8,1)+X(3,7)*H&
&_invDmH(6,1)*H_invDmH(7,2)+X(3,7)*H_invDmH(6,2)*H_invDmH(7,1)+2*X&
&(3,6)*H_invDmH(6,1)*H_invDmH(6,2)+X(3,5)*H_invDmH(5,1)*H_invDmH(6&
&,2)+X(3,4)*H_invDmH(4,1)*H_invDmH(6,2)+H_invDmH(3,1)*X(3,3)*H_inv&
&DmH(6,2)+H_invDmH(2,1)*X(3,2)*H_invDmH(6,2)+H_invDmH(1,1)*X(3,1)*&
&H_invDmH(6,2)+X(3,5)*H_invDmH(5,2)*H_invDmH(6,1)+X(3,4)*H_invDmH(&
&4,2)*H_invDmH(6,1)+H_invDmH(3,2)*X(3,3)*H_invDmH(6,1)+H_invDmH(2,&
&2)*X(3,2)*H_invDmH(6,1)+H_invDmH(1,2)*X(3,1)*H_invDmH(6,1))*tt38+&
&H_invDmH(6,1)*tt5*tt6)+lam(1,1)*(H_invDmH(6,3)*tt11+H_invDmH(6,2)&
&*tt8+H_invDmH(6,1)*tt5)*tt13)
jac(1,19) = detDmH(1,1)*gw(1,1)*(mu(1,1)*(H_invDmH(7,3)*tt3*tt12+&
&(X(1,8)*H_invDmH(7,2)*H_invDmH(8,3)+X(1,8)*H_invDmH(7,3)*H_invDmH&
&(8,2)+2*X(1,7)*H_invDmH(7,2)*H_invDmH(7,3)+X(1,6)*H_invDmH(6,2)*H&
&_invDmH(7,3)+X(1,5)*H_invDmH(5,2)*H_invDmH(7,3)+X(1,4)*H_invDmH(4&
&,2)*H_invDmH(7,3)+X(1,3)*H_invDmH(3,2)*H_invDmH(7,3)+X(1,2)*H_inv&
&DmH(2,2)*H_invDmH(7,3)+X(1,1)*H_invDmH(1,2)*H_invDmH(7,3)+X(1,6)*&
&H_invDmH(6,3)*H_invDmH(7,2)+X(1,5)*H_invDmH(5,3)*H_invDmH(7,2)+X(&
&1,4)*H_invDmH(4,3)*H_invDmH(7,2)+X(1,3)*H_invDmH(3,3)*H_invDmH(7,&
&2)+X(1,2)*H_invDmH(2,3)*H_invDmH(7,2)+X(1,1)*H_invDmH(1,3)*H_invD&
&mH(7,2))*tt40+(X(1,8)*H_invDmH(7,1)*H_invDmH(8,3)+X(1,8)*H_invDmH&
&(7,3)*H_invDmH(8,1)+2*X(1,7)*H_invDmH(7,1)*H_invDmH(7,3)+X(1,6)*H&
&_invDmH(6,1)*H_invDmH(7,3)+X(1,5)*H_invDmH(5,1)*H_invDmH(7,3)+X(1&
&,4)*H_invDmH(4,1)*H_invDmH(7,3)+X(1,3)*H_invDmH(3,1)*H_invDmH(7,3&
&)+X(1,2)*H_invDmH(2,1)*H_invDmH(7,3)+H_invDmH(1,1)*X(1,1)*H_invDm&
&H(7,3)+X(1,6)*H_invDmH(6,3)*H_invDmH(7,1)+X(1,5)*H_invDmH(5,3)*H_&
&invDmH(7,1)+X(1,4)*H_invDmH(4,3)*H_invDmH(7,1)+X(1,3)*H_invDmH(3,&
&3)*H_invDmH(7,1)+X(1,2)*H_invDmH(2,3)*H_invDmH(7,1)+X(1,1)*H_invD&
&mH(1,3)*H_invDmH(7,1))*tt39+H_invDmH(7,2)*tt2*tt9+(X(1,8)*H_invDm&
&H(7,1)*H_invDmH(8,2)+X(1,8)*H_invDmH(7,2)*H_invDmH(8,1)+2*X(1,7)*&
&H_invDmH(7,1)*H_invDmH(7,2)+X(1,6)*H_invDmH(6,1)*H_invDmH(7,2)+X(&
&1,5)*H_invDmH(5,1)*H_invDmH(7,2)+X(1,4)*H_invDmH(4,1)*H_invDmH(7,&
&2)+X(1,3)*H_invDmH(3,1)*H_invDmH(7,2)+X(1,2)*H_invDmH(2,1)*H_invD&
&mH(7,2)+H_invDmH(1,1)*X(1,1)*H_invDmH(7,2)+X(1,6)*H_invDmH(6,2)*H&
&_invDmH(7,1)+X(1,5)*H_invDmH(5,2)*H_invDmH(7,1)+X(1,4)*H_invDmH(4&
&,2)*H_invDmH(7,1)+X(1,3)*H_invDmH(3,2)*H_invDmH(7,1)+X(1,2)*H_inv&
&DmH(2,2)*H_invDmH(7,1)+X(1,1)*H_invDmH(1,2)*H_invDmH(7,1))*tt38+H&
&_invDmH(7,1)*tt1*tt6)+lam(1,1)*(H_invDmH(7,3)*tt3+H_invDmH(7,2)*t&
&t2+H_invDmH(7,1)*tt1)*tt13)
jac(1,20) = detDmH(1,1)*gw(1,1)*(mu(1,1)*(H_invDmH(7,3)*tt10*tt12&
&+(X(2,8)*H_invDmH(7,2)*H_invDmH(8,3)+X(2,8)*H_invDmH(7,3)*H_invDm&
&H(8,2)+2*X(2,7)*H_invDmH(7,2)*H_invDmH(7,3)+X(2,6)*H_invDmH(6,2)*&
&H_invDmH(7,3)+X(2,5)*H_invDmH(5,2)*H_invDmH(7,3)+X(2,4)*H_invDmH(&
&4,2)*H_invDmH(7,3)+X(2,3)*H_invDmH(3,2)*H_invDmH(7,3)+H_invDmH(2,&
&2)*X(2,2)*H_invDmH(7,3)+H_invDmH(1,2)*X(2,1)*H_invDmH(7,3)+X(2,6)&
&*H_invDmH(6,3)*H_invDmH(7,2)+X(2,5)*H_invDmH(5,3)*H_invDmH(7,2)+X&
&(2,4)*H_invDmH(4,3)*H_invDmH(7,2)+X(2,3)*H_invDmH(3,3)*H_invDmH(7&
&,2)+X(2,2)*H_invDmH(2,3)*H_invDmH(7,2)+H_invDmH(1,3)*X(2,1)*H_inv&
&DmH(7,2))*tt40+(X(2,8)*H_invDmH(7,1)*H_invDmH(8,3)+X(2,8)*H_invDm&
&H(7,3)*H_invDmH(8,1)+2*X(2,7)*H_invDmH(7,1)*H_invDmH(7,3)+X(2,6)*&
&H_invDmH(6,1)*H_invDmH(7,3)+X(2,5)*H_invDmH(5,1)*H_invDmH(7,3)+X(&
&2,4)*H_invDmH(4,1)*H_invDmH(7,3)+X(2,3)*H_invDmH(3,1)*H_invDmH(7,&
&3)+H_invDmH(2,1)*X(2,2)*H_invDmH(7,3)+H_invDmH(1,1)*X(2,1)*H_invD&
&mH(7,3)+X(2,6)*H_invDmH(6,3)*H_invDmH(7,1)+X(2,5)*H_invDmH(5,3)*H&
&_invDmH(7,1)+X(2,4)*H_invDmH(4,3)*H_invDmH(7,1)+X(2,3)*H_invDmH(3&
&,3)*H_invDmH(7,1)+X(2,2)*H_invDmH(2,3)*H_invDmH(7,1)+H_invDmH(1,3&
&)*X(2,1)*H_invDmH(7,1))*tt39+H_invDmH(7,2)*tt7*tt9+(X(2,8)*H_invD&
&mH(7,1)*H_invDmH(8,2)+X(2,8)*H_invDmH(7,2)*H_invDmH(8,1)+2*X(2,7)&
&*H_invDmH(7,1)*H_invDmH(7,2)+X(2,6)*H_invDmH(6,1)*H_invDmH(7,2)+X&
&(2,5)*H_invDmH(5,1)*H_invDmH(7,2)+X(2,4)*H_invDmH(4,1)*H_invDmH(7&
&,2)+X(2,3)*H_invDmH(3,1)*H_invDmH(7,2)+H_invDmH(2,1)*X(2,2)*H_inv&
&DmH(7,2)+H_invDmH(1,1)*X(2,1)*H_invDmH(7,2)+X(2,6)*H_invDmH(6,2)*&
&H_invDmH(7,1)+X(2,5)*H_invDmH(5,2)*H_invDmH(7,1)+X(2,4)*H_invDmH(&
&4,2)*H_invDmH(7,1)+X(2,3)*H_invDmH(3,2)*H_invDmH(7,1)+H_invDmH(2,&
&2)*X(2,2)*H_invDmH(7,1)+H_invDmH(1,2)*X(2,1)*H_invDmH(7,1))*tt38+&
&H_invDmH(7,1)*tt4*tt6)+lam(1,1)*(H_invDmH(7,3)*tt10+H_invDmH(7,2)&
&*tt7+H_invDmH(7,1)*tt4)*tt13)
jac(1,21) = detDmH(1,1)*gw(1,1)*(mu(1,1)*(H_invDmH(7,3)*tt11*tt12&
&+(X(3,8)*H_invDmH(7,2)*H_invDmH(8,3)+X(3,8)*H_invDmH(7,3)*H_invDm&
&H(8,2)+2*X(3,7)*H_invDmH(7,2)*H_invDmH(7,3)+X(3,6)*H_invDmH(6,2)*&
&H_invDmH(7,3)+X(3,5)*H_invDmH(5,2)*H_invDmH(7,3)+X(3,4)*H_invDmH(&
&4,2)*H_invDmH(7,3)+H_invDmH(3,2)*X(3,3)*H_invDmH(7,3)+H_invDmH(2,&
&2)*X(3,2)*H_invDmH(7,3)+H_invDmH(1,2)*X(3,1)*H_invDmH(7,3)+X(3,6)&
&*H_invDmH(6,3)*H_invDmH(7,2)+X(3,5)*H_invDmH(5,3)*H_invDmH(7,2)+X&
&(3,4)*H_invDmH(4,3)*H_invDmH(7,2)+H_invDmH(3,3)*X(3,3)*H_invDmH(7&
&,2)+H_invDmH(2,3)*X(3,2)*H_invDmH(7,2)+H_invDmH(1,3)*X(3,1)*H_inv&
&DmH(7,2))*tt40+(X(3,8)*H_invDmH(7,1)*H_invDmH(8,3)+X(3,8)*H_invDm&
&H(7,3)*H_invDmH(8,1)+2*X(3,7)*H_invDmH(7,1)*H_invDmH(7,3)+X(3,6)*&
&H_invDmH(6,1)*H_invDmH(7,3)+X(3,5)*H_invDmH(5,1)*H_invDmH(7,3)+X(&
&3,4)*H_invDmH(4,1)*H_invDmH(7,3)+H_invDmH(3,1)*X(3,3)*H_invDmH(7,&
&3)+H_invDmH(2,1)*X(3,2)*H_invDmH(7,3)+H_invDmH(1,1)*X(3,1)*H_invD&
&mH(7,3)+X(3,6)*H_invDmH(6,3)*H_invDmH(7,1)+X(3,5)*H_invDmH(5,3)*H&
&_invDmH(7,1)+X(3,4)*H_invDmH(4,3)*H_invDmH(7,1)+H_invDmH(3,3)*X(3&
&,3)*H_invDmH(7,1)+H_invDmH(2,3)*X(3,2)*H_invDmH(7,1)+H_invDmH(1,3&
&)*X(3,1)*H_invDmH(7,1))*tt39+H_invDmH(7,2)*tt8*tt9+(X(3,8)*H_invD&
&mH(7,1)*H_invDmH(8,2)+X(3,8)*H_invDmH(7,2)*H_invDmH(8,1)+2*X(3,7)&
&*H_invDmH(7,1)*H_invDmH(7,2)+X(3,6)*H_invDmH(6,1)*H_invDmH(7,2)+X&
&(3,5)*H_invDmH(5,1)*H_invDmH(7,2)+X(3,4)*H_invDmH(4,1)*H_invDmH(7&
&,2)+H_invDmH(3,1)*X(3,3)*H_invDmH(7,2)+H_invDmH(2,1)*X(3,2)*H_inv&
&DmH(7,2)+H_invDmH(1,1)*X(3,1)*H_invDmH(7,2)+X(3,6)*H_invDmH(6,2)*&
&H_invDmH(7,1)+X(3,5)*H_invDmH(5,2)*H_invDmH(7,1)+X(3,4)*H_invDmH(&
&4,2)*H_invDmH(7,1)+H_invDmH(3,2)*X(3,3)*H_invDmH(7,1)+H_invDmH(2,&
&2)*X(3,2)*H_invDmH(7,1)+H_invDmH(1,2)*X(3,1)*H_invDmH(7,1))*tt38+&
&H_invDmH(7,1)*tt5*tt6)+lam(1,1)*(H_invDmH(7,3)*tt11+H_invDmH(7,2)&
&*tt8+H_invDmH(7,1)*tt5)*tt13)
jac(1,22) = detDmH(1,1)*gw(1,1)*(mu(1,1)*(H_invDmH(8,3)*tt3*tt12+&
&(2*X(1,8)*H_invDmH(8,2)*H_invDmH(8,3)+X(1,7)*H_invDmH(7,2)*H_invD&
&mH(8,3)+X(1,6)*H_invDmH(6,2)*H_invDmH(8,3)+X(1,5)*H_invDmH(5,2)*H&
&_invDmH(8,3)+X(1,4)*H_invDmH(4,2)*H_invDmH(8,3)+X(1,3)*H_invDmH(3&
&,2)*H_invDmH(8,3)+X(1,2)*H_invDmH(2,2)*H_invDmH(8,3)+X(1,1)*H_inv&
&DmH(1,2)*H_invDmH(8,3)+X(1,7)*H_invDmH(7,3)*H_invDmH(8,2)+X(1,6)*&
&H_invDmH(6,3)*H_invDmH(8,2)+X(1,5)*H_invDmH(5,3)*H_invDmH(8,2)+X(&
&1,4)*H_invDmH(4,3)*H_invDmH(8,2)+X(1,3)*H_invDmH(3,3)*H_invDmH(8,&
&2)+X(1,2)*H_invDmH(2,3)*H_invDmH(8,2)+X(1,1)*H_invDmH(1,3)*H_invD&
&mH(8,2))*tt40+(2*X(1,8)*H_invDmH(8,1)*H_invDmH(8,3)+X(1,7)*H_invD&
&mH(7,1)*H_invDmH(8,3)+X(1,6)*H_invDmH(6,1)*H_invDmH(8,3)+X(1,5)*H&
&_invDmH(5,1)*H_invDmH(8,3)+X(1,4)*H_invDmH(4,1)*H_invDmH(8,3)+X(1&
&,3)*H_invDmH(3,1)*H_invDmH(8,3)+X(1,2)*H_invDmH(2,1)*H_invDmH(8,3&
&)+H_invDmH(1,1)*X(1,1)*H_invDmH(8,3)+X(1,7)*H_invDmH(7,3)*H_invDm&
&H(8,1)+X(1,6)*H_invDmH(6,3)*H_invDmH(8,1)+X(1,5)*H_invDmH(5,3)*H_&
&invDmH(8,1)+X(1,4)*H_invDmH(4,3)*H_invDmH(8,1)+X(1,3)*H_invDmH(3,&
&3)*H_invDmH(8,1)+X(1,2)*H_invDmH(2,3)*H_invDmH(8,1)+X(1,1)*H_invD&
&mH(1,3)*H_invDmH(8,1))*tt39+H_invDmH(8,2)*tt2*tt9+(2*X(1,8)*H_inv&
&DmH(8,1)*H_invDmH(8,2)+X(1,7)*H_invDmH(7,1)*H_invDmH(8,2)+X(1,6)*&
&H_invDmH(6,1)*H_invDmH(8,2)+X(1,5)*H_invDmH(5,1)*H_invDmH(8,2)+X(&
&1,4)*H_invDmH(4,1)*H_invDmH(8,2)+X(1,3)*H_invDmH(3,1)*H_invDmH(8,&
&2)+X(1,2)*H_invDmH(2,1)*H_invDmH(8,2)+H_invDmH(1,1)*X(1,1)*H_invD&
&mH(8,2)+X(1,7)*H_invDmH(7,2)*H_invDmH(8,1)+X(1,6)*H_invDmH(6,2)*H&
&_invDmH(8,1)+X(1,5)*H_invDmH(5,2)*H_invDmH(8,1)+X(1,4)*H_invDmH(4&
&,2)*H_invDmH(8,1)+X(1,3)*H_invDmH(3,2)*H_invDmH(8,1)+X(1,2)*H_inv&
&DmH(2,2)*H_invDmH(8,1)+X(1,1)*H_invDmH(1,2)*H_invDmH(8,1))*tt38+H&
&_invDmH(8,1)*tt1*tt6)+lam(1,1)*(H_invDmH(8,3)*tt3+H_invDmH(8,2)*t&
&t2+H_invDmH(8,1)*tt1)*tt13)
jac(1,23) = detDmH(1,1)*gw(1,1)*(mu(1,1)*(H_invDmH(8,3)*tt10*tt12&
&+(2*X(2,8)*H_invDmH(8,2)*H_invDmH(8,3)+X(2,7)*H_invDmH(7,2)*H_inv&
&DmH(8,3)+X(2,6)*H_invDmH(6,2)*H_invDmH(8,3)+X(2,5)*H_invDmH(5,2)*&
&H_invDmH(8,3)+X(2,4)*H_invDmH(4,2)*H_invDmH(8,3)+X(2,3)*H_invDmH(&
&3,2)*H_invDmH(8,3)+H_invDmH(2,2)*X(2,2)*H_invDmH(8,3)+H_invDmH(1,&
&2)*X(2,1)*H_invDmH(8,3)+X(2,7)*H_invDmH(7,3)*H_invDmH(8,2)+X(2,6)&
&*H_invDmH(6,3)*H_invDmH(8,2)+X(2,5)*H_invDmH(5,3)*H_invDmH(8,2)+X&
&(2,4)*H_invDmH(4,3)*H_invDmH(8,2)+X(2,3)*H_invDmH(3,3)*H_invDmH(8&
&,2)+X(2,2)*H_invDmH(2,3)*H_invDmH(8,2)+H_invDmH(1,3)*X(2,1)*H_inv&
&DmH(8,2))*tt40+(2*X(2,8)*H_invDmH(8,1)*H_invDmH(8,3)+X(2,7)*H_inv&
&DmH(7,1)*H_invDmH(8,3)+X(2,6)*H_invDmH(6,1)*H_invDmH(8,3)+X(2,5)*&
&H_invDmH(5,1)*H_invDmH(8,3)+X(2,4)*H_invDmH(4,1)*H_invDmH(8,3)+X(&
&2,3)*H_invDmH(3,1)*H_invDmH(8,3)+H_invDmH(2,1)*X(2,2)*H_invDmH(8,&
&3)+H_invDmH(1,1)*X(2,1)*H_invDmH(8,3)+X(2,7)*H_invDmH(7,3)*H_invD&
&mH(8,1)+X(2,6)*H_invDmH(6,3)*H_invDmH(8,1)+X(2,5)*H_invDmH(5,3)*H&
&_invDmH(8,1)+X(2,4)*H_invDmH(4,3)*H_invDmH(8,1)+X(2,3)*H_invDmH(3&
&,3)*H_invDmH(8,1)+X(2,2)*H_invDmH(2,3)*H_invDmH(8,1)+H_invDmH(1,3&
&)*X(2,1)*H_invDmH(8,1))*tt39+H_invDmH(8,2)*tt7*tt9+(2*X(2,8)*H_in&
&vDmH(8,1)*H_invDmH(8,2)+X(2,7)*H_invDmH(7,1)*H_invDmH(8,2)+X(2,6)&
&*H_invDmH(6,1)*H_invDmH(8,2)+X(2,5)*H_invDmH(5,1)*H_invDmH(8,2)+X&
&(2,4)*H_invDmH(4,1)*H_invDmH(8,2)+X(2,3)*H_invDmH(3,1)*H_invDmH(8&
&,2)+H_invDmH(2,1)*X(2,2)*H_invDmH(8,2)+H_invDmH(1,1)*X(2,1)*H_inv&
&DmH(8,2)+X(2,7)*H_invDmH(7,2)*H_invDmH(8,1)+X(2,6)*H_invDmH(6,2)*&
&H_invDmH(8,1)+X(2,5)*H_invDmH(5,2)*H_invDmH(8,1)+X(2,4)*H_invDmH(&
&4,2)*H_invDmH(8,1)+X(2,3)*H_invDmH(3,2)*H_invDmH(8,1)+H_invDmH(2,&
&2)*X(2,2)*H_invDmH(8,1)+H_invDmH(1,2)*X(2,1)*H_invDmH(8,1))*tt38+&
&H_invDmH(8,1)*tt4*tt6)+lam(1,1)*(H_invDmH(8,3)*tt10+H_invDmH(8,2)&
&*tt7+H_invDmH(8,1)*tt4)*tt13)
jac(1,24) = detDmH(1,1)*gw(1,1)*(mu(1,1)*(H_invDmH(8,3)*tt11*tt12&
&+(2*X(3,8)*H_invDmH(8,2)*H_invDmH(8,3)+X(3,7)*H_invDmH(7,2)*H_inv&
&DmH(8,3)+X(3,6)*H_invDmH(6,2)*H_invDmH(8,3)+X(3,5)*H_invDmH(5,2)*&
&H_invDmH(8,3)+X(3,4)*H_invDmH(4,2)*H_invDmH(8,3)+H_invDmH(3,2)*X(&
&3,3)*H_invDmH(8,3)+H_invDmH(2,2)*X(3,2)*H_invDmH(8,3)+H_invDmH(1,&
&2)*X(3,1)*H_invDmH(8,3)+X(3,7)*H_invDmH(7,3)*H_invDmH(8,2)+X(3,6)&
&*H_invDmH(6,3)*H_invDmH(8,2)+X(3,5)*H_invDmH(5,3)*H_invDmH(8,2)+X&
&(3,4)*H_invDmH(4,3)*H_invDmH(8,2)+H_invDmH(3,3)*X(3,3)*H_invDmH(8&
&,2)+H_invDmH(2,3)*X(3,2)*H_invDmH(8,2)+H_invDmH(1,3)*X(3,1)*H_inv&
&DmH(8,2))*tt40+(2*X(3,8)*H_invDmH(8,1)*H_invDmH(8,3)+X(3,7)*H_inv&
&DmH(7,1)*H_invDmH(8,3)+X(3,6)*H_invDmH(6,1)*H_invDmH(8,3)+X(3,5)*&
&H_invDmH(5,1)*H_invDmH(8,3)+X(3,4)*H_invDmH(4,1)*H_invDmH(8,3)+H_&
&invDmH(3,1)*X(3,3)*H_invDmH(8,3)+H_invDmH(2,1)*X(3,2)*H_invDmH(8,&
&3)+H_invDmH(1,1)*X(3,1)*H_invDmH(8,3)+X(3,7)*H_invDmH(7,3)*H_invD&
&mH(8,1)+X(3,6)*H_invDmH(6,3)*H_invDmH(8,1)+X(3,5)*H_invDmH(5,3)*H&
&_invDmH(8,1)+X(3,4)*H_invDmH(4,3)*H_invDmH(8,1)+H_invDmH(3,3)*X(3&
&,3)*H_invDmH(8,1)+H_invDmH(2,3)*X(3,2)*H_invDmH(8,1)+H_invDmH(1,3&
&)*X(3,1)*H_invDmH(8,1))*tt39+H_invDmH(8,2)*tt8*tt9+(2*X(3,8)*H_in&
&vDmH(8,1)*H_invDmH(8,2)+X(3,7)*H_invDmH(7,1)*H_invDmH(8,2)+X(3,6)&
&*H_invDmH(6,1)*H_invDmH(8,2)+X(3,5)*H_invDmH(5,1)*H_invDmH(8,2)+X&
&(3,4)*H_invDmH(4,1)*H_invDmH(8,2)+H_invDmH(3,1)*X(3,3)*H_invDmH(8&
&,2)+H_invDmH(2,1)*X(3,2)*H_invDmH(8,2)+H_invDmH(1,1)*X(3,1)*H_inv&
&DmH(8,2)+X(3,7)*H_invDmH(7,2)*H_invDmH(8,1)+X(3,6)*H_invDmH(6,2)*&
&H_invDmH(8,1)+X(3,5)*H_invDmH(5,2)*H_invDmH(8,1)+X(3,4)*H_invDmH(&
&4,2)*H_invDmH(8,1)+H_invDmH(3,2)*X(3,3)*H_invDmH(8,1)+H_invDmH(2,&
&2)*X(3,2)*H_invDmH(8,1)+H_invDmH(1,2)*X(3,1)*H_invDmH(8,1))*tt38+&
&H_invDmH(8,1)*tt5*tt6)+lam(1,1)*(H_invDmH(8,3)*tt11+H_invDmH(8,2)&
&*tt8+H_invDmH(8,1)*tt5)*tt13)
END 
SUBROUTINE vox_stvk_at_quadr_hes(hes, X, H_invDmH, detDmH, gw, mu, lam) 
IMPLICIT NONE 
REAL(KIND=8) hes(24, 24) 
REAL(KIND=8) X(3, 8) 
REAL(KIND=8) H_invDmH(8, 3) 
REAL(KIND=8) detDmH(1, 1) 
REAL(KIND=8) gw(1, 1) 
REAL(KIND=8) mu(1, 1) 
REAL(KIND=8) lam(1, 1) 
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
REAL(KIND=8)  tt73 
REAL(KIND=8)  tt74 
REAL(KIND=8)  tt75 
REAL(KIND=8)  tt76 
REAL(KIND=8)  tt77 
REAL(KIND=8)  tt78 
REAL(KIND=8)  tt79 
REAL(KIND=8)  tt80 
REAL(KIND=8)  tt81 
REAL(KIND=8)  tt82 
REAL(KIND=8)  tt83 
REAL(KIND=8)  tt84 
REAL(KIND=8)  tt85 
REAL(KIND=8)  tt86 
REAL(KIND=8)  tt87 
REAL(KIND=8)  tt88 
REAL(KIND=8)  tt89 
REAL(KIND=8)  tt90 
REAL(KIND=8)  tt91 
REAL(KIND=8)  tt92 
REAL(KIND=8)  tt93 
REAL(KIND=8)  tt94 
REAL(KIND=8)  tt95 
REAL(KIND=8)  tt96 
REAL(KIND=8)  tt97 
REAL(KIND=8)  tt98 
REAL(KIND=8)  tt99 
REAL(KIND=8)  tt100 
REAL(KIND=8)  tt101 
REAL(KIND=8)  tt102 
REAL(KIND=8)  tt103 
REAL(KIND=8)  tt104 
REAL(KIND=8)  tt105 
REAL(KIND=8)  tt106 
REAL(KIND=8)  tt107 
REAL(KIND=8)  tt108 
REAL(KIND=8)  tt109 
REAL(KIND=8)  tt110 
REAL(KIND=8)  tt111 
REAL(KIND=8)  tt112 
REAL(KIND=8)  tt113 
REAL(KIND=8)  tt114 
REAL(KIND=8)  tt115 
REAL(KIND=8)  tt116 
REAL(KIND=8)  tt117 
REAL(KIND=8)  tt118 
REAL(KIND=8)  tt119 
REAL(KIND=8)  tt120 
REAL(KIND=8)  tt121 
REAL(KIND=8)  tt122 
REAL(KIND=8)  tt123 
REAL(KIND=8)  tt124 
REAL(KIND=8)  tt125 
REAL(KIND=8)  tt126 
REAL(KIND=8)  tt127 
REAL(KIND=8)  tt128 
REAL(KIND=8)  tt129 
REAL(KIND=8)  tt130 
REAL(KIND=8)  tt131 
REAL(KIND=8)  tt132 
REAL(KIND=8)  tt133 
REAL(KIND=8)  tt134 
REAL(KIND=8)  tt135 
REAL(KIND=8)  tt136 
REAL(KIND=8)  tt137 
REAL(KIND=8)  tt138 
REAL(KIND=8)  tt139 
REAL(KIND=8)  tt140 
REAL(KIND=8)  tt141 
REAL(KIND=8)  tt142 
REAL(KIND=8)  tt143 
REAL(KIND=8)  tt144 
REAL(KIND=8)  tt145 
REAL(KIND=8)  tt146 
REAL(KIND=8)  tt147 
REAL(KIND=8)  tt148 
REAL(KIND=8)  tt149 
REAL(KIND=8)  tt150 
REAL(KIND=8)  tt151 
REAL(KIND=8)  tt152 
REAL(KIND=8)  tt153 
REAL(KIND=8)  tt154 
REAL(KIND=8)  tt155 
REAL(KIND=8)  tt156 
REAL(KIND=8)  tt157 
REAL(KIND=8)  tt158 
REAL(KIND=8)  tt159 
REAL(KIND=8)  tt160 
REAL(KIND=8)  tt161 
REAL(KIND=8)  tt162 
REAL(KIND=8)  tt163 
REAL(KIND=8)  tt164 
REAL(KIND=8)  tt165 
REAL(KIND=8)  tt166 
REAL(KIND=8)  tt167 
REAL(KIND=8)  tt168 
REAL(KIND=8)  tt169 
REAL(KIND=8)  tt170 
REAL(KIND=8)  tt171 
REAL(KIND=8)  tt172 
REAL(KIND=8)  tt173 
REAL(KIND=8)  tt174 
REAL(KIND=8)  tt175 
REAL(KIND=8)  tt176 
REAL(KIND=8)  tt177 
REAL(KIND=8)  tt178 
REAL(KIND=8)  tt179 
REAL(KIND=8)  tt180 
REAL(KIND=8)  tt181 
REAL(KIND=8)  tt182 
REAL(KIND=8)  tt183 
REAL(KIND=8)  tt184 
REAL(KIND=8)  tt185 
REAL(KIND=8)  tt186 
REAL(KIND=8)  tt187 
REAL(KIND=8)  tt188 
REAL(KIND=8)  tt189 
REAL(KIND=8)  tt190 
REAL(KIND=8)  tt191 
REAL(KIND=8)  tt192 
REAL(KIND=8)  tt193 
REAL(KIND=8)  tt194 
REAL(KIND=8)  tt195 
REAL(KIND=8)  tt196 
REAL(KIND=8)  tt197 
REAL(KIND=8)  tt198 
REAL(KIND=8)  tt199 
REAL(KIND=8)  tt200 
REAL(KIND=8)  tt201 
REAL(KIND=8)  tt202 
REAL(KIND=8)  tt203 
REAL(KIND=8)  tt204 
REAL(KIND=8)  tt205 
REAL(KIND=8)  tt206 
REAL(KIND=8)  tt207 
REAL(KIND=8)  tt208 
REAL(KIND=8)  tt209 
REAL(KIND=8)  tt210 
REAL(KIND=8)  tt211 
REAL(KIND=8)  tt212 
REAL(KIND=8)  tt213 
REAL(KIND=8)  tt214 
REAL(KIND=8)  tt215 
REAL(KIND=8)  tt216 
REAL(KIND=8)  tt217 
REAL(KIND=8)  tt218 
REAL(KIND=8)  tt219 
REAL(KIND=8)  tt220 
REAL(KIND=8)  tt221 
REAL(KIND=8)  tt222 
REAL(KIND=8)  tt223 
REAL(KIND=8)  tt224 
REAL(KIND=8)  tt225 
REAL(KIND=8)  tt226 
REAL(KIND=8)  tt227 
REAL(KIND=8)  tt228 
REAL(KIND=8)  tt229 
REAL(KIND=8)  tt230 
REAL(KIND=8)  tt231 
REAL(KIND=8)  tt232 
REAL(KIND=8)  tt233 
REAL(KIND=8)  tt234 
REAL(KIND=8)  tt235 
REAL(KIND=8)  tt236 
REAL(KIND=8)  tt237 
REAL(KIND=8)  tt238 
REAL(KIND=8)  tt239 
REAL(KIND=8)  tt240 
REAL(KIND=8)  tt241 
REAL(KIND=8)  tt242 
REAL(KIND=8)  tt243 
REAL(KIND=8)  tt244 
REAL(KIND=8)  tt245 
REAL(KIND=8)  tt246 
REAL(KIND=8)  tt247 
REAL(KIND=8)  tt248 
REAL(KIND=8)  tt249 
REAL(KIND=8)  tt250 
REAL(KIND=8)  tt251 
REAL(KIND=8)  tt252 
REAL(KIND=8)  tt253 
REAL(KIND=8)  tt254 
REAL(KIND=8)  tt255 
REAL(KIND=8)  tt256 
REAL(KIND=8)  tt257 
REAL(KIND=8)  tt258 
REAL(KIND=8)  tt259 
REAL(KIND=8)  tt260 
REAL(KIND=8)  tt261 
REAL(KIND=8)  tt262 
REAL(KIND=8)  tt263 
REAL(KIND=8)  tt264 
REAL(KIND=8)  tt265 
REAL(KIND=8)  tt266 
REAL(KIND=8)  tt267 
REAL(KIND=8)  tt268 
REAL(KIND=8)  tt269 
REAL(KIND=8)  tt270 
REAL(KIND=8)  tt271 
REAL(KIND=8)  tt272 
REAL(KIND=8)  tt273 
REAL(KIND=8)  tt274 
REAL(KIND=8)  tt275 
REAL(KIND=8)  tt276 
REAL(KIND=8)  tt277 
REAL(KIND=8)  tt278 
REAL(KIND=8)  tt279 
REAL(KIND=8)  tt280 
REAL(KIND=8)  tt281 
REAL(KIND=8)  tt282 
REAL(KIND=8)  tt283 
REAL(KIND=8)  tt284 
REAL(KIND=8)  tt285 
REAL(KIND=8)  tt286 
REAL(KIND=8)  tt287 
REAL(KIND=8)  tt288 
REAL(KIND=8)  tt289 
REAL(KIND=8)  tt290 
REAL(KIND=8)  tt291 
REAL(KIND=8)  tt292 
REAL(KIND=8)  tt293 
REAL(KIND=8)  tt294 
REAL(KIND=8)  tt295 
REAL(KIND=8)  tt296 
REAL(KIND=8)  tt297 
REAL(KIND=8)  tt298 
REAL(KIND=8)  tt299 
REAL(KIND=8)  tt300 
REAL(KIND=8)  tt301 
REAL(KIND=8)  tt302 
REAL(KIND=8)  tt303 
REAL(KIND=8)  tt304 
REAL(KIND=8)  tt305 
REAL(KIND=8)  tt306 
REAL(KIND=8)  tt307 
REAL(KIND=8)  tt308 
REAL(KIND=8)  tt309 
REAL(KIND=8)  tt310 
REAL(KIND=8)  tt311 
REAL(KIND=8)  tt312 
REAL(KIND=8)  tt313 
REAL(KIND=8)  tt314 
REAL(KIND=8)  tt315 
REAL(KIND=8)  tt316 
REAL(KIND=8)  tt317 
REAL(KIND=8)  tt318 
REAL(KIND=8)  tt319 
REAL(KIND=8)  tt320 
REAL(KIND=8)  tt321 
REAL(KIND=8)  tt322 
REAL(KIND=8)  tt323 
REAL(KIND=8)  tt324 
REAL(KIND=8)  tt325 
REAL(KIND=8)  tt326 
REAL(KIND=8)  tt327 
REAL(KIND=8)  tt328 
REAL(KIND=8)  tt329 
REAL(KIND=8)  tt330 
REAL(KIND=8)  tt331 
REAL(KIND=8)  tt332 
REAL(KIND=8)  tt333 
REAL(KIND=8)  tt334 
REAL(KIND=8)  tt335 
REAL(KIND=8)  tt336 
REAL(KIND=8)  tt337 
REAL(KIND=8)  tt338 
REAL(KIND=8)  tt339 
REAL(KIND=8)  tt340 
REAL(KIND=8)  tt341 
REAL(KIND=8)  tt342 
REAL(KIND=8)  tt343 
REAL(KIND=8)  tt344 
REAL(KIND=8)  tt345 
REAL(KIND=8)  tt346 
REAL(KIND=8)  tt347 
REAL(KIND=8)  tt348 
REAL(KIND=8)  tt349 
REAL(KIND=8)  tt350 
REAL(KIND=8)  tt351 
REAL(KIND=8)  tt352 
REAL(KIND=8)  tt353 
REAL(KIND=8)  tt354 
REAL(KIND=8)  tt355 
REAL(KIND=8)  tt356 
REAL(KIND=8)  tt357 
REAL(KIND=8)  tt358 
REAL(KIND=8)  tt359 
REAL(KIND=8)  tt360 
REAL(KIND=8)  tt361 
REAL(KIND=8)  tt362 
REAL(KIND=8)  tt363 
REAL(KIND=8)  tt364 
REAL(KIND=8)  tt365 
REAL(KIND=8)  tt366 
REAL(KIND=8)  tt367 
REAL(KIND=8)  tt368 
REAL(KIND=8)  tt369 
REAL(KIND=8)  tt370 
REAL(KIND=8)  tt371 
REAL(KIND=8)  tt372 
REAL(KIND=8)  tt373 
REAL(KIND=8)  tt374 
REAL(KIND=8)  tt375 
REAL(KIND=8)  tt376 
REAL(KIND=8)  tt377 
REAL(KIND=8)  tt378 
REAL(KIND=8)  tt379 
REAL(KIND=8)  tt380 
REAL(KIND=8)  tt381 
REAL(KIND=8)  tt382 
REAL(KIND=8)  tt383 
REAL(KIND=8)  tt384 
REAL(KIND=8)  tt385 
REAL(KIND=8)  tt386 
REAL(KIND=8)  tt387 
REAL(KIND=8)  tt388 
REAL(KIND=8)  tt389 
REAL(KIND=8)  tt390 
REAL(KIND=8)  tt391 
REAL(KIND=8)  tt392 
REAL(KIND=8)  tt393 
REAL(KIND=8)  tt394 
REAL(KIND=8)  tt395 
REAL(KIND=8)  tt396 
REAL(KIND=8)  tt397 
REAL(KIND=8)  tt398 
REAL(KIND=8)  tt399 
REAL(KIND=8)  tt400 
REAL(KIND=8)  tt401 
REAL(KIND=8)  tt402 
REAL(KIND=8)  tt403 
REAL(KIND=8)  tt404 
REAL(KIND=8)  tt405 
REAL(KIND=8)  tt406 
REAL(KIND=8)  tt407 
REAL(KIND=8)  tt408 
REAL(KIND=8)  tt409 
REAL(KIND=8)  tt410 
REAL(KIND=8)  tt411 
REAL(KIND=8)  tt412 
REAL(KIND=8)  tt413 
REAL(KIND=8)  tt414 
REAL(KIND=8)  tt415 
REAL(KIND=8)  tt416 
REAL(KIND=8)  tt417 
REAL(KIND=8)  tt418 
REAL(KIND=8)  tt419 
REAL(KIND=8)  tt420 
REAL(KIND=8)  tt421 
REAL(KIND=8)  tt422 
REAL(KIND=8)  tt423 
REAL(KIND=8)  tt424 
REAL(KIND=8)  tt425 
REAL(KIND=8)  tt426 
REAL(KIND=8)  tt427 
REAL(KIND=8)  tt428 
REAL(KIND=8)  tt429 
REAL(KIND=8)  tt430 
REAL(KIND=8)  tt431 
REAL(KIND=8)  tt432 
REAL(KIND=8)  tt433 
REAL(KIND=8)  tt434 
REAL(KIND=8)  tt435 
REAL(KIND=8)  tt436 
REAL(KIND=8)  tt437 
REAL(KIND=8)  tt438 
REAL(KIND=8)  tt439 
REAL(KIND=8)  tt440 
REAL(KIND=8)  tt441 
REAL(KIND=8)  tt442 
REAL(KIND=8)  tt443 
REAL(KIND=8)  tt444 
REAL(KIND=8)  tt445 
REAL(KIND=8)  tt446 
REAL(KIND=8)  tt447 
REAL(KIND=8)  tt448 
REAL(KIND=8)  tt449 
REAL(KIND=8)  tt450 
REAL(KIND=8)  tt451 
REAL(KIND=8)  tt452 
REAL(KIND=8)  tt453 
REAL(KIND=8)  tt454 
REAL(KIND=8)  tt455 
REAL(KIND=8)  tt456 
REAL(KIND=8)  tt457 
REAL(KIND=8)  tt458 
REAL(KIND=8)  tt459 
REAL(KIND=8)  tt460 
REAL(KIND=8)  tt461 
REAL(KIND=8)  tt462 
REAL(KIND=8)  tt463 
REAL(KIND=8)  tt464 
REAL(KIND=8)  tt465 
REAL(KIND=8)  tt466 
REAL(KIND=8)  tt467 
REAL(KIND=8)  tt468 
REAL(KIND=8)  tt469 
REAL(KIND=8)  tt470 
REAL(KIND=8)  tt471 
REAL(KIND=8)  tt472 
REAL(KIND=8)  tt473 
REAL(KIND=8)  tt474 
REAL(KIND=8)  tt475 
REAL(KIND=8)  tt476 
REAL(KIND=8)  tt477 
REAL(KIND=8)  tt478 
REAL(KIND=8)  tt479 
REAL(KIND=8)  tt480 
REAL(KIND=8)  tt481 
REAL(KIND=8)  tt482 
REAL(KIND=8)  tt483 
REAL(KIND=8)  tt484 
REAL(KIND=8)  tt485 
REAL(KIND=8)  tt486 
REAL(KIND=8)  tt487 
REAL(KIND=8)  tt488 
REAL(KIND=8)  tt489 
REAL(KIND=8)  tt490 
REAL(KIND=8)  tt491 
REAL(KIND=8)  tt492 
REAL(KIND=8)  tt493 
REAL(KIND=8)  tt494 
REAL(KIND=8)  tt495 
REAL(KIND=8)  tt496 
REAL(KIND=8)  tt497 
REAL(KIND=8)  tt498 
REAL(KIND=8)  tt499 
REAL(KIND=8)  tt500 
REAL(KIND=8)  tt501 
REAL(KIND=8)  tt502 
REAL(KIND=8)  tt503 
REAL(KIND=8)  tt504 
REAL(KIND=8)  tt505 
REAL(KIND=8)  tt506 
REAL(KIND=8)  tt507 
REAL(KIND=8)  tt508 
REAL(KIND=8)  tt509 
REAL(KIND=8)  tt510 
REAL(KIND=8)  tt511 
REAL(KIND=8)  tt512 
REAL(KIND=8)  tt513 
REAL(KIND=8)  tt514 
REAL(KIND=8)  tt515 
REAL(KIND=8)  tt516 
REAL(KIND=8)  tt517 
REAL(KIND=8)  tt518 
REAL(KIND=8)  tt519 
REAL(KIND=8)  tt520 
REAL(KIND=8)  tt521 
REAL(KIND=8)  tt522 
REAL(KIND=8)  tt523 
REAL(KIND=8)  tt524 
REAL(KIND=8)  tt525 
REAL(KIND=8)  tt526 
REAL(KIND=8)  tt527 
REAL(KIND=8)  tt528 
REAL(KIND=8)  tt529 
REAL(KIND=8)  tt530 
REAL(KIND=8)  tt531 
REAL(KIND=8)  tt532 
REAL(KIND=8)  tt533 
REAL(KIND=8)  tt534 
REAL(KIND=8)  tt535 
REAL(KIND=8)  tt536 
REAL(KIND=8)  tt537 
REAL(KIND=8)  tt538 
REAL(KIND=8)  tt539 
REAL(KIND=8)  tt540 
REAL(KIND=8)  tt541 
REAL(KIND=8)  tt542 
REAL(KIND=8)  tt543 
REAL(KIND=8)  tt544 
REAL(KIND=8)  tt545 
REAL(KIND=8)  tt546 
REAL(KIND=8)  tt547 
REAL(KIND=8)  tt548 
REAL(KIND=8)  tt549 
REAL(KIND=8)  tt550 
REAL(KIND=8)  tt551 
REAL(KIND=8)  tt552 
REAL(KIND=8)  tt553 
REAL(KIND=8)  tt554 
REAL(KIND=8)  tt555 
REAL(KIND=8)  tt556 
REAL(KIND=8)  tt557 
REAL(KIND=8)  tt558 
REAL(KIND=8)  tt559 
REAL(KIND=8)  tt560 
REAL(KIND=8)  tt561 
REAL(KIND=8)  tt562 
REAL(KIND=8)  tt563 
REAL(KIND=8)  tt564 
REAL(KIND=8)  tt565 
REAL(KIND=8)  tt566 
REAL(KIND=8)  tt567 
REAL(KIND=8)  tt568 
REAL(KIND=8)  tt569 
REAL(KIND=8)  tt570 
REAL(KIND=8)  tt571 
REAL(KIND=8)  tt572 
REAL(KIND=8)  tt573 
REAL(KIND=8)  tt574 
REAL(KIND=8)  tt575 
REAL(KIND=8)  tt576 
REAL(KIND=8)  tt577 
REAL(KIND=8)  tt578 
REAL(KIND=8)  tt579 
REAL(KIND=8)  tt580 
REAL(KIND=8)  tt581 
REAL(KIND=8)  tt582 
REAL(KIND=8)  tt583 
REAL(KIND=8)  tt584 
REAL(KIND=8)  tt585 
REAL(KIND=8)  tt586 
REAL(KIND=8)  tt587 
REAL(KIND=8)  tt588 
REAL(KIND=8)  tt589 
REAL(KIND=8)  tt590 
REAL(KIND=8)  tt591 
REAL(KIND=8)  tt592 
REAL(KIND=8)  tt593 
REAL(KIND=8)  tt594 
REAL(KIND=8)  tt595 
REAL(KIND=8)  tt596 
REAL(KIND=8)  tt597 
REAL(KIND=8)  tt598 
REAL(KIND=8)  tt599 
REAL(KIND=8)  tt600 
REAL(KIND=8)  tt601 
REAL(KIND=8)  tt602 
REAL(KIND=8)  tt603 
REAL(KIND=8)  tt604 
REAL(KIND=8)  tt605 
REAL(KIND=8)  tt606 
REAL(KIND=8)  tt607 
REAL(KIND=8)  tt608 
REAL(KIND=8)  tt609 
REAL(KIND=8)  tt610 
REAL(KIND=8)  tt611 
REAL(KIND=8)  tt612 
REAL(KIND=8)  tt613 
REAL(KIND=8)  tt614 
REAL(KIND=8)  tt615 
REAL(KIND=8)  tt616 
REAL(KIND=8)  tt617 
REAL(KIND=8)  tt618 
REAL(KIND=8)  tt619 
REAL(KIND=8)  tt620 
REAL(KIND=8)  tt621 
REAL(KIND=8)  tt622 
REAL(KIND=8)  tt623 
REAL(KIND=8)  tt624 
REAL(KIND=8)  tt625 
REAL(KIND=8)  tt626 
REAL(KIND=8)  tt627 
REAL(KIND=8)  tt628 
REAL(KIND=8)  tt629 
REAL(KIND=8)  tt630 
REAL(KIND=8)  tt631 
REAL(KIND=8)  tt632 
REAL(KIND=8)  tt633 
REAL(KIND=8)  tt634 
REAL(KIND=8)  tt635 
REAL(KIND=8)  tt636 
REAL(KIND=8)  tt637 
REAL(KIND=8)  tt638 
REAL(KIND=8)  tt639 
REAL(KIND=8)  tt640 
REAL(KIND=8)  tt641 
REAL(KIND=8)  tt642 
REAL(KIND=8)  tt643 
REAL(KIND=8)  tt644 
REAL(KIND=8)  tt645 
REAL(KIND=8)  tt646 
REAL(KIND=8)  tt647 
REAL(KIND=8)  tt648 
REAL(KIND=8)  tt649 
REAL(KIND=8)  tt650 
REAL(KIND=8)  tt651 
REAL(KIND=8)  tt652 
REAL(KIND=8)  tt653 
REAL(KIND=8)  tt654 
REAL(KIND=8)  tt655 
REAL(KIND=8)  tt656 
REAL(KIND=8)  tt657 
REAL(KIND=8)  tt658 
REAL(KIND=8)  tt659 
REAL(KIND=8)  tt660 
REAL(KIND=8)  tt661 
REAL(KIND=8)  tt662 
REAL(KIND=8)  tt663 
REAL(KIND=8)  tt664 
REAL(KIND=8)  tt665 
REAL(KIND=8)  tt666 
REAL(KIND=8)  tt667 
REAL(KIND=8)  tt668 
REAL(KIND=8)  tt669 
REAL(KIND=8)  tt670 
REAL(KIND=8)  tt671 
REAL(KIND=8)  tt672 
REAL(KIND=8)  tt673 
REAL(KIND=8)  tt674 
REAL(KIND=8)  tt675 
REAL(KIND=8)  tt676 
REAL(KIND=8)  tt677 
REAL(KIND=8)  tt678 
REAL(KIND=8)  tt679 
REAL(KIND=8)  tt680 
REAL(KIND=8)  tt681 
REAL(KIND=8)  tt682 
REAL(KIND=8)  tt683 
REAL(KIND=8)  tt684 
REAL(KIND=8)  tt685 
REAL(KIND=8)  tt686 
REAL(KIND=8)  tt687 
REAL(KIND=8)  tt688 
REAL(KIND=8)  tt689 
REAL(KIND=8)  tt690 
REAL(KIND=8)  tt691 
REAL(KIND=8)  tt692 
REAL(KIND=8)  tt693 
REAL(KIND=8)  tt694 
REAL(KIND=8)  tt695 
REAL(KIND=8)  tt696 
REAL(KIND=8)  tt697 
REAL(KIND=8)  tt698 
REAL(KIND=8)  tt699 
REAL(KIND=8)  tt700 
REAL(KIND=8)  tt701 
REAL(KIND=8)  tt702 
REAL(KIND=8)  tt703 
REAL(KIND=8)  tt704 
REAL(KIND=8)  tt705 
REAL(KIND=8)  tt706 
REAL(KIND=8)  tt707 
REAL(KIND=8)  tt708 
REAL(KIND=8)  tt709 
REAL(KIND=8)  tt710 
REAL(KIND=8)  tt711 
REAL(KIND=8)  tt712 
REAL(KIND=8)  tt713 
REAL(KIND=8)  tt714 
REAL(KIND=8)  tt715 
REAL(KIND=8)  tt716 
REAL(KIND=8)  tt717 
REAL(KIND=8)  tt718 
REAL(KIND=8)  tt719 
REAL(KIND=8)  tt720 
REAL(KIND=8)  tt721 
REAL(KIND=8)  tt722 
REAL(KIND=8)  tt723 
REAL(KIND=8)  tt724 
REAL(KIND=8)  tt725 
REAL(KIND=8)  tt726 
REAL(KIND=8)  tt727 
REAL(KIND=8)  tt728 
REAL(KIND=8)  tt729 
REAL(KIND=8)  tt730 
REAL(KIND=8)  tt731 
REAL(KIND=8)  tt732 
REAL(KIND=8)  tt733 
REAL(KIND=8)  tt734 
REAL(KIND=8)  tt735 
REAL(KIND=8)  tt736 
REAL(KIND=8)  tt737 
REAL(KIND=8)  tt738 
REAL(KIND=8)  tt739 
REAL(KIND=8)  tt740 
REAL(KIND=8)  tt741 
REAL(KIND=8)  tt742 
REAL(KIND=8)  tt743 
REAL(KIND=8)  tt744 
REAL(KIND=8)  tt745 
REAL(KIND=8)  tt746 
REAL(KIND=8)  tt747 
REAL(KIND=8)  tt748 
REAL(KIND=8)  tt749 
REAL(KIND=8)  tt750 
REAL(KIND=8)  tt751 
REAL(KIND=8)  tt752 
REAL(KIND=8)  tt753 
REAL(KIND=8)  tt754 
REAL(KIND=8)  tt755 
REAL(KIND=8)  tt756 
REAL(KIND=8)  tt757 
REAL(KIND=8)  tt758 
REAL(KIND=8)  tt759 
REAL(KIND=8)  tt760 
REAL(KIND=8)  tt761 
REAL(KIND=8)  tt762 
REAL(KIND=8)  tt763 
REAL(KIND=8)  tt764 
REAL(KIND=8)  tt765 
REAL(KIND=8)  tt766 
REAL(KIND=8)  tt767 
REAL(KIND=8)  tt768 
REAL(KIND=8)  tt769 
REAL(KIND=8)  tt770 
REAL(KIND=8)  tt771 
REAL(KIND=8)  tt772 
REAL(KIND=8)  tt773 
REAL(KIND=8)  tt774 
REAL(KIND=8)  tt775 
REAL(KIND=8)  tt776 
REAL(KIND=8)  tt777 
REAL(KIND=8)  tt778 
REAL(KIND=8)  tt779 
REAL(KIND=8)  tt780 
REAL(KIND=8)  tt781 
REAL(KIND=8)  tt782 
REAL(KIND=8)  tt783 
REAL(KIND=8)  tt784 
REAL(KIND=8)  tt785 
REAL(KIND=8)  tt786 
REAL(KIND=8)  tt787 
REAL(KIND=8)  tt788 
REAL(KIND=8)  tt789 
REAL(KIND=8)  tt790 
REAL(KIND=8)  tt791 
REAL(KIND=8)  tt792 
REAL(KIND=8)  tt793 
REAL(KIND=8)  tt794 
REAL(KIND=8)  tt795 
REAL(KIND=8)  tt796 
REAL(KIND=8)  tt797 
REAL(KIND=8)  tt798 
REAL(KIND=8)  tt799 
REAL(KIND=8)  tt800 
REAL(KIND=8)  tt801 
REAL(KIND=8)  tt802 
REAL(KIND=8)  tt803 
REAL(KIND=8)  tt804 
REAL(KIND=8)  tt805 
REAL(KIND=8)  tt806 
REAL(KIND=8)  tt807 
REAL(KIND=8)  tt808 
REAL(KIND=8)  tt809 
REAL(KIND=8)  tt810 
REAL(KIND=8)  tt811 
REAL(KIND=8)  tt812 
REAL(KIND=8)  tt813 
REAL(KIND=8)  tt814 
REAL(KIND=8)  tt815 
REAL(KIND=8)  tt816 
REAL(KIND=8)  tt817 
REAL(KIND=8)  tt818 
REAL(KIND=8)  tt819 
REAL(KIND=8)  tt820 
REAL(KIND=8)  tt821 
REAL(KIND=8)  tt822 
REAL(KIND=8)  tt823 
REAL(KIND=8)  tt824 
REAL(KIND=8)  tt825 
REAL(KIND=8)  tt826 
REAL(KIND=8)  tt827 
REAL(KIND=8)  tt828 
REAL(KIND=8)  tt829 
REAL(KIND=8)  tt830 
REAL(KIND=8)  tt831 
REAL(KIND=8)  tt832 
REAL(KIND=8)  tt833 
REAL(KIND=8)  tt834 
REAL(KIND=8)  tt835 
REAL(KIND=8)  tt836 
REAL(KIND=8)  tt837 
REAL(KIND=8)  tt838 
REAL(KIND=8)  tt839 
REAL(KIND=8)  tt840 
REAL(KIND=8)  tt841 
REAL(KIND=8)  tt842 
REAL(KIND=8)  tt843 
REAL(KIND=8)  tt844 
REAL(KIND=8)  tt845 
REAL(KIND=8)  tt846 
REAL(KIND=8)  tt847 
REAL(KIND=8)  tt848 
REAL(KIND=8)  tt849 
REAL(KIND=8)  tt850 
REAL(KIND=8)  tt851 
REAL(KIND=8)  tt852 
REAL(KIND=8)  tt853 
REAL(KIND=8)  tt854 
REAL(KIND=8)  tt855 
REAL(KIND=8)  tt856 
REAL(KIND=8)  tt857 
REAL(KIND=8)  tt858 
REAL(KIND=8)  tt859 
REAL(KIND=8)  tt860 
REAL(KIND=8)  tt861 
REAL(KIND=8)  tt862 
REAL(KIND=8)  tt863 
REAL(KIND=8)  tt864 
REAL(KIND=8)  tt865 
REAL(KIND=8)  tt866 
REAL(KIND=8)  tt867 
REAL(KIND=8)  tt868 
REAL(KIND=8)  tt869 
REAL(KIND=8)  tt870 
REAL(KIND=8)  tt871 
REAL(KIND=8)  tt872 
REAL(KIND=8)  tt873 
REAL(KIND=8)  tt874 
REAL(KIND=8)  tt875 
REAL(KIND=8)  tt876 
REAL(KIND=8)  tt877 
REAL(KIND=8)  tt878 
REAL(KIND=8)  tt879 
REAL(KIND=8)  tt880 
REAL(KIND=8)  tt881 
REAL(KIND=8)  tt882 
REAL(KIND=8)  tt883 
REAL(KIND=8)  tt884 
REAL(KIND=8)  tt885 
REAL(KIND=8)  tt886 
REAL(KIND=8)  tt887 
REAL(KIND=8)  tt888 
REAL(KIND=8)  tt889 
REAL(KIND=8)  tt890 
REAL(KIND=8)  tt891 
REAL(KIND=8)  tt892 
REAL(KIND=8)  tt893 
REAL(KIND=8)  tt894 
REAL(KIND=8)  tt895 
REAL(KIND=8)  tt896 
REAL(KIND=8)  tt897 
REAL(KIND=8)  tt898 
REAL(KIND=8)  tt899 
REAL(KIND=8)  tt900 
REAL(KIND=8)  tt901 
REAL(KIND=8)  tt902 
REAL(KIND=8)  tt903 
REAL(KIND=8)  tt904 
REAL(KIND=8)  tt905 
REAL(KIND=8)  tt906 
REAL(KIND=8)  tt907 
REAL(KIND=8)  tt908 
REAL(KIND=8)  tt909 
REAL(KIND=8)  tt910 
REAL(KIND=8)  tt911 
REAL(KIND=8)  tt912 
REAL(KIND=8)  tt913 
REAL(KIND=8)  tt914 
REAL(KIND=8)  tt915 
REAL(KIND=8)  tt916 
REAL(KIND=8)  tt917 
REAL(KIND=8)  tt918 
REAL(KIND=8)  tt919 
REAL(KIND=8)  tt920 
REAL(KIND=8)  tt921 
REAL(KIND=8)  tt922 
REAL(KIND=8)  tt923 
REAL(KIND=8)  tt924 
REAL(KIND=8)  tt925 
REAL(KIND=8)  tt926 
REAL(KIND=8)  tt927 
REAL(KIND=8)  tt928 
REAL(KIND=8)  tt929 
REAL(KIND=8)  tt930 
REAL(KIND=8)  tt931 
REAL(KIND=8)  tt932 
REAL(KIND=8)  tt933 
REAL(KIND=8)  tt934 
REAL(KIND=8)  tt935 
REAL(KIND=8)  tt936 
REAL(KIND=8)  tt937 
REAL(KIND=8)  tt938 
REAL(KIND=8)  tt939 
REAL(KIND=8)  tt940 
REAL(KIND=8)  tt941 
REAL(KIND=8)  tt942 
REAL(KIND=8)  tt943 
REAL(KIND=8)  tt944 
REAL(KIND=8)  tt945 
REAL(KIND=8)  tt946 
REAL(KIND=8)  tt947 
REAL(KIND=8)  tt948 
REAL(KIND=8)  tt949 
tt1 = H_invDmH(1,1)**2
tt2 = H_invDmH(1,2)**2
tt3 = H_invDmH(1,3)**2
tt4 = X(1,8)*H_invDmH(8,1)+X(1,7)*H_invDmH(7,1)+X(1,6)*H_invDmH(6&
&,1)+X(1,5)*H_invDmH(5,1)+X(1,4)*H_invDmH(4,1)+X(1,3)*H_invDmH(3,1&
&)+X(1,2)*H_invDmH(2,1)+H_invDmH(1,1)*X(1,1)
tt5 = tt4**2
tt6 = X(2,8)*H_invDmH(8,1)+X(2,7)*H_invDmH(7,1)+X(2,6)*H_invDmH(6&
&,1)+X(2,5)*H_invDmH(5,1)+X(2,4)*H_invDmH(4,1)+X(2,3)*H_invDmH(3,1&
&)+H_invDmH(2,1)*X(2,2)+H_invDmH(1,1)*X(2,1)
tt7 = tt6**2
tt8 = X(3,8)*H_invDmH(8,1)+X(3,7)*H_invDmH(7,1)+X(3,6)*H_invDmH(6&
&,1)+X(3,5)*H_invDmH(5,1)+X(3,4)*H_invDmH(4,1)+H_invDmH(3,1)*X(3,3&
&)+H_invDmH(2,1)*X(3,2)+H_invDmH(1,1)*X(3,1)
tt9 = tt8**2
tt10 = tt9+tt7+tt5-1
tt11 = X(1,8)*H_invDmH(8,2)+X(1,7)*H_invDmH(7,2)+X(1,6)*H_invDmH(&
&6,2)+X(1,5)*H_invDmH(5,2)+X(1,4)*H_invDmH(4,2)+X(1,3)*H_invDmH(3,&
&2)+X(1,2)*H_invDmH(2,2)+X(1,1)*H_invDmH(1,2)
tt12 = tt11**2
tt13 = X(2,8)*H_invDmH(8,2)+X(2,7)*H_invDmH(7,2)+X(2,6)*H_invDmH(&
&6,2)+X(2,5)*H_invDmH(5,2)+X(2,4)*H_invDmH(4,2)+X(2,3)*H_invDmH(3,&
&2)+H_invDmH(2,2)*X(2,2)+H_invDmH(1,2)*X(2,1)
tt14 = tt13**2
tt15 = X(3,8)*H_invDmH(8,2)+X(3,7)*H_invDmH(7,2)+X(3,6)*H_invDmH(&
&6,2)+X(3,5)*H_invDmH(5,2)+X(3,4)*H_invDmH(4,2)+H_invDmH(3,2)*X(3,&
&3)+H_invDmH(2,2)*X(3,2)+H_invDmH(1,2)*X(3,1)
tt16 = tt15**2
tt17 = tt16+tt14+tt12-1
tt18 = X(1,8)*H_invDmH(8,3)+X(1,7)*H_invDmH(7,3)+X(1,6)*H_invDmH(&
&6,3)+X(1,5)*H_invDmH(5,3)+X(1,4)*H_invDmH(4,3)+X(1,3)*H_invDmH(3,&
&3)+X(1,2)*H_invDmH(2,3)+X(1,1)*H_invDmH(1,3)
tt19 = tt18**2
tt20 = X(2,8)*H_invDmH(8,3)+X(2,7)*H_invDmH(7,3)+X(2,6)*H_invDmH(&
&6,3)+X(2,5)*H_invDmH(5,3)+X(2,4)*H_invDmH(4,3)+X(2,3)*H_invDmH(3,&
&3)+X(2,2)*H_invDmH(2,3)+H_invDmH(1,3)*X(2,1)
tt21 = tt20**2
tt22 = X(3,8)*H_invDmH(8,3)+X(3,7)*H_invDmH(7,3)+X(3,6)*H_invDmH(&
&6,3)+X(3,5)*H_invDmH(5,3)+X(3,4)*H_invDmH(4,3)+H_invDmH(3,3)*X(3,&
&3)+H_invDmH(2,3)*X(3,2)+H_invDmH(1,3)*X(3,1)
tt23 = tt22**2
tt24 = tt23+tt21+tt19-1
tt25 = tt24/2.0E+0+tt17/2.0E+0+tt10/2.0E+0
tt26 = lam(1,1)*(tt3+tt2+tt1)*tt25
tt27 = tt1*tt10
tt28 = X(1,1)**2
tt29 = X(2,1)**2
tt30 = X(1,2)**2
tt31 = X(2,2)**2
tt32 = X(3,1)**2
tt33 = X(1,3)**2
tt34 = X(2,3)**2
tt35 = X(3,2)**2
tt36 = X(3,3)**2
tt37 = X(1,4)**2
tt38 = X(2,4)**2
tt39 = X(3,4)**2
tt40 = X(1,5)**2
tt41 = X(2,5)**2
tt42 = X(3,5)**2
tt43 = X(1,6)**2
tt44 = X(2,6)**2
tt45 = X(3,6)**2
tt46 = X(1,7)**2
tt47 = X(2,7)**2
tt48 = X(3,7)**2
tt49 = X(1,8)**2
tt50 = X(2,8)**2
tt51 = X(3,8)**2
tt52 = tt51*H_invDmH(8,1)*H_invDmH(8,2)+tt50*H_invDmH(8,1)*H_invD&
&mH(8,2)+tt49*H_invDmH(8,1)*H_invDmH(8,2)+X(3,7)*X(3,8)*H_invDmH(7&
&,1)*H_invDmH(8,2)+X(2,7)*X(2,8)*H_invDmH(7,1)*H_invDmH(8,2)+X(1,7&
&)*X(1,8)*H_invDmH(7,1)*H_invDmH(8,2)+X(3,6)*X(3,8)*H_invDmH(6,1)*&
&H_invDmH(8,2)+X(2,6)*X(2,8)*H_invDmH(6,1)*H_invDmH(8,2)+X(1,6)*X(&
&1,8)*H_invDmH(6,1)*H_invDmH(8,2)+X(3,5)*X(3,8)*H_invDmH(5,1)*H_in&
&vDmH(8,2)+X(2,5)*X(2,8)*H_invDmH(5,1)*H_invDmH(8,2)+X(1,5)*X(1,8)&
&*H_invDmH(5,1)*H_invDmH(8,2)+X(3,4)*X(3,8)*H_invDmH(4,1)*H_invDmH&
&(8,2)+X(2,4)*X(2,8)*H_invDmH(4,1)*H_invDmH(8,2)+X(1,4)*X(1,8)*H_i&
&nvDmH(4,1)*H_invDmH(8,2)+H_invDmH(3,1)*X(3,3)*X(3,8)*H_invDmH(8,2&
&)+H_invDmH(2,1)*X(3,2)*X(3,8)*H_invDmH(8,2)+H_invDmH(1,1)*X(3,1)*&
&X(3,8)*H_invDmH(8,2)+X(2,3)*X(2,8)*H_invDmH(3,1)*H_invDmH(8,2)+X(&
&1,3)*X(1,8)*H_invDmH(3,1)*H_invDmH(8,2)+H_invDmH(2,1)*X(2,2)*X(2,&
&8)*H_invDmH(8,2)+H_invDmH(1,1)*X(2,1)*X(2,8)*H_invDmH(8,2)+X(1,2)&
&*X(1,8)*H_invDmH(2,1)*H_invDmH(8,2)+H_invDmH(1,1)*X(1,1)*X(1,8)*H&
&_invDmH(8,2)+X(3,7)*X(3,8)*H_invDmH(7,2)*H_invDmH(8,1)+X(2,7)*X(2&
&,8)*H_invDmH(7,2)*H_invDmH(8,1)+X(1,7)*X(1,8)*H_invDmH(7,2)*H_inv&
&DmH(8,1)+X(3,6)*X(3,8)*H_invDmH(6,2)*H_invDmH(8,1)+X(2,6)*X(2,8)*&
&H_invDmH(6,2)*H_invDmH(8,1)+X(1,6)*X(1,8)*H_invDmH(6,2)*H_invDmH(&
&8,1)+X(3,5)*X(3,8)*H_invDmH(5,2)*H_invDmH(8,1)+X(2,5)*X(2,8)*H_in&
&vDmH(5,2)*H_invDmH(8,1)+X(1,5)*X(1,8)*H_invDmH(5,2)*H_invDmH(8,1)&
&+X(3,4)*X(3,8)*H_invDmH(4,2)*H_invDmH(8,1)+X(2,4)*X(2,8)*H_invDmH&
&(4,2)*H_invDmH(8,1)+X(1,4)*X(1,8)*H_invDmH(4,2)*H_invDmH(8,1)+H_i&
&nvDmH(3,2)*X(3,3)*X(3,8)*H_invDmH(8,1)+H_invDmH(2,2)*X(3,2)*X(3,8&
&)*H_invDmH(8,1)+H_invDmH(1,2)*X(3,1)*X(3,8)*H_invDmH(8,1)+X(2,3)*&
&X(2,8)*H_invDmH(3,2)*H_invDmH(8,1)+X(1,3)*X(1,8)*H_invDmH(3,2)*H_&
&invDmH(8,1)+H_invDmH(2,2)*X(2,2)*X(2,8)*H_invDmH(8,1)+H_invDmH(1,&
&2)*X(2,1)*X(2,8)*H_invDmH(8,1)+X(1,2)*X(1,8)*H_invDmH(2,2)*H_invD&
&mH(8,1)+X(1,1)*H_invDmH(1,2)*X(1,8)*H_invDmH(8,1)+tt48*H_invDmH(7&
&,1)*H_invDmH(7,2)+tt47*H_invDmH(7,1)*H_invDmH(7,2)+tt46*H_invDmH(&
&7,1)*H_invDmH(7,2)+X(3,6)*X(3,7)*H_invDmH(6,1)*H_invDmH(7,2)+X(2,&
&6)*X(2,7)*H_invDmH(6,1)*H_invDmH(7,2)+X(1,6)*X(1,7)*H_invDmH(6,1)&
&*H_invDmH(7,2)+X(3,5)*X(3,7)*H_invDmH(5,1)*H_invDmH(7,2)+X(2,5)*X&
&(2,7)*H_invDmH(5,1)*H_invDmH(7,2)+X(1,5)*X(1,7)*H_invDmH(5,1)*H_i&
&nvDmH(7,2)+X(3,4)*X(3,7)*H_invDmH(4,1)*H_invDmH(7,2)+X(2,4)*X(2,7&
&)*H_invDmH(4,1)*H_invDmH(7,2)+X(1,4)*X(1,7)*H_invDmH(4,1)*H_invDm&
&H(7,2)+H_invDmH(3,1)*X(3,3)*X(3,7)*H_invDmH(7,2)+H_invDmH(2,1)*X(&
&3,2)*X(3,7)*H_invDmH(7,2)+H_invDmH(1,1)*X(3,1)*X(3,7)*H_invDmH(7,&
&2)+X(2,3)*X(2,7)*H_invDmH(3,1)*H_invDmH(7,2)+X(1,3)*X(1,7)*H_invD&
&mH(3,1)*H_invDmH(7,2)+H_invDmH(2,1)*X(2,2)*X(2,7)*H_invDmH(7,2)+H&
&_invDmH(1,1)*X(2,1)*X(2,7)*H_invDmH(7,2)+X(1,2)*X(1,7)*H_invDmH(2&
&,1)*H_invDmH(7,2)+H_invDmH(1,1)*X(1,1)*X(1,7)*H_invDmH(7,2)+X(3,6&
&)*X(3,7)*H_invDmH(6,2)*H_invDmH(7,1)+X(2,6)*X(2,7)*H_invDmH(6,2)*&
&H_invDmH(7,1)+X(1,6)*X(1,7)*H_invDmH(6,2)*H_invDmH(7,1)+X(3,5)*X(&
&3,7)*H_invDmH(5,2)*H_invDmH(7,1)+X(2,5)*X(2,7)*H_invDmH(5,2)*H_in&
&vDmH(7,1)+X(1,5)*X(1,7)*H_invDmH(5,2)*H_invDmH(7,1)+X(3,4)*X(3,7)&
&*H_invDmH(4,2)*H_invDmH(7,1)+X(2,4)*X(2,7)*H_invDmH(4,2)*H_invDmH&
&(7,1)+X(1,4)*X(1,7)*H_invDmH(4,2)*H_invDmH(7,1)+H_invDmH(3,2)*X(3&
&,3)*X(3,7)*H_invDmH(7,1)+H_invDmH(2,2)*X(3,2)*X(3,7)*H_invDmH(7,1&
&)+H_invDmH(1,2)*X(3,1)*X(3,7)*H_invDmH(7,1)+X(2,3)*X(2,7)*H_invDm&
&H(3,2)*H_invDmH(7,1)+X(1,3)*X(1,7)*H_invDmH(3,2)*H_invDmH(7,1)+H_&
&invDmH(2,2)*X(2,2)*X(2,7)*H_invDmH(7,1)+H_invDmH(1,2)*X(2,1)*X(2,&
&7)*H_invDmH(7,1)+X(1,2)*X(1,7)*H_invDmH(2,2)*H_invDmH(7,1)+X(1,1)&
&*H_invDmH(1,2)*X(1,7)*H_invDmH(7,1)+tt45*H_invDmH(6,1)*H_invDmH(6&
&,2)+tt44*H_invDmH(6,1)*H_invDmH(6,2)+tt43*H_invDmH(6,1)*H_invDmH(&
&6,2)+X(3,5)*X(3,6)*H_invDmH(5,1)*H_invDmH(6,2)+X(2,5)*X(2,6)*H_in&
&vDmH(5,1)*H_invDmH(6,2)+X(1,5)*X(1,6)*H_invDmH(5,1)*H_invDmH(6,2)&
&+X(3,4)*X(3,6)*H_invDmH(4,1)*H_invDmH(6,2)+X(2,4)*X(2,6)*H_invDmH&
&(4,1)*H_invDmH(6,2)+X(1,4)*X(1,6)*H_invDmH(4,1)*H_invDmH(6,2)+H_i&
&nvDmH(3,1)*X(3,3)*X(3,6)*H_invDmH(6,2)+H_invDmH(2,1)*X(3,2)*X(3,6&
&)*H_invDmH(6,2)+H_invDmH(1,1)*X(3,1)*X(3,6)*H_invDmH(6,2)+X(2,3)*&
&X(2,6)*H_invDmH(3,1)*H_invDmH(6,2)+X(1,3)*X(1,6)*H_invDmH(3,1)*H_&
&invDmH(6,2)+H_invDmH(2,1)*X(2,2)*X(2,6)*H_invDmH(6,2)+H_invDmH(1,&
&1)*X(2,1)*X(2,6)*H_invDmH(6,2)+X(1,2)*X(1,6)*H_invDmH(2,1)*H_invD&
&mH(6,2)+H_invDmH(1,1)*X(1,1)*X(1,6)*H_invDmH(6,2)+X(3,5)*X(3,6)*H&
&_invDmH(5,2)*H_invDmH(6,1)+X(2,5)*X(2,6)*H_invDmH(5,2)*H_invDmH(6&
&,1)+X(1,5)*X(1,6)*H_invDmH(5,2)*H_invDmH(6,1)+X(3,4)*X(3,6)*H_inv&
&DmH(4,2)*H_invDmH(6,1)+X(2,4)*X(2,6)*H_invDmH(4,2)*H_invDmH(6,1)+&
&X(1,4)*X(1,6)*H_invDmH(4,2)*H_invDmH(6,1)+H_invDmH(3,2)*X(3,3)*X(&
&3,6)*H_invDmH(6,1)+H_invDmH(2,2)*X(3,2)*X(3,6)*H_invDmH(6,1)+H_in&
&vDmH(1,2)*X(3,1)*X(3,6)*H_invDmH(6,1)+X(2,3)*X(2,6)*H_invDmH(3,2)&
&*H_invDmH(6,1)+X(1,3)*X(1,6)*H_invDmH(3,2)*H_invDmH(6,1)+H_invDmH&
&(2,2)*X(2,2)*X(2,6)*H_invDmH(6,1)+H_invDmH(1,2)*X(2,1)*X(2,6)*H_i&
&nvDmH(6,1)+X(1,2)*X(1,6)*H_invDmH(2,2)*H_invDmH(6,1)+X(1,1)*H_inv&
&DmH(1,2)*X(1,6)*H_invDmH(6,1)+tt42*H_invDmH(5,1)*H_invDmH(5,2)+tt&
&41*H_invDmH(5,1)*H_invDmH(5,2)+tt40*H_invDmH(5,1)*H_invDmH(5,2)+X&
&(3,4)*X(3,5)*H_invDmH(4,1)*H_invDmH(5,2)+X(2,4)*X(2,5)*H_invDmH(4&
&,1)*H_invDmH(5,2)+X(1,4)*X(1,5)*H_invDmH(4,1)*H_invDmH(5,2)+H_inv&
&DmH(3,1)*X(3,3)*X(3,5)*H_invDmH(5,2)+H_invDmH(2,1)*X(3,2)*X(3,5)*&
&H_invDmH(5,2)+H_invDmH(1,1)*X(3,1)*X(3,5)*H_invDmH(5,2)+X(2,3)*X(&
&2,5)*H_invDmH(3,1)*H_invDmH(5,2)+X(1,3)*X(1,5)*H_invDmH(3,1)*H_in&
&vDmH(5,2)+H_invDmH(2,1)*X(2,2)*X(2,5)*H_invDmH(5,2)+H_invDmH(1,1)&
&*X(2,1)*X(2,5)*H_invDmH(5,2)+X(1,2)*X(1,5)*H_invDmH(2,1)*H_invDmH&
&(5,2)+H_invDmH(1,1)*X(1,1)*X(1,5)*H_invDmH(5,2)+X(3,4)*X(3,5)*H_i&
&nvDmH(4,2)*H_invDmH(5,1)+X(2,4)*X(2,5)*H_invDmH(4,2)*H_invDmH(5,1&
&)+X(1,4)*X(1,5)*H_invDmH(4,2)*H_invDmH(5,1)+H_invDmH(3,2)*X(3,3)*&
&X(3,5)*H_invDmH(5,1)+H_invDmH(2,2)*X(3,2)*X(3,5)*H_invDmH(5,1)+H_&
&invDmH(1,2)*X(3,1)*X(3,5)*H_invDmH(5,1)+X(2,3)*X(2,5)*H_invDmH(3,&
&2)*H_invDmH(5,1)+X(1,3)*X(1,5)*H_invDmH(3,2)*H_invDmH(5,1)+H_invD&
&mH(2,2)*X(2,2)*X(2,5)*H_invDmH(5,1)+H_invDmH(1,2)*X(2,1)*X(2,5)*H&
&_invDmH(5,1)+X(1,2)*X(1,5)*H_invDmH(2,2)*H_invDmH(5,1)+X(1,1)*H_i&
&nvDmH(1,2)*X(1,5)*H_invDmH(5,1)+tt39*H_invDmH(4,1)*H_invDmH(4,2)+&
&tt38*H_invDmH(4,1)*H_invDmH(4,2)+tt37*H_invDmH(4,1)*H_invDmH(4,2)&
&+H_invDmH(3,1)*X(3,3)*X(3,4)*H_invDmH(4,2)+H_invDmH(2,1)*X(3,2)*X&
&(3,4)*H_invDmH(4,2)+H_invDmH(1,1)*X(3,1)*X(3,4)*H_invDmH(4,2)+X(2&
&,3)*X(2,4)*H_invDmH(3,1)*H_invDmH(4,2)+X(1,3)*X(1,4)*H_invDmH(3,1&
&)*H_invDmH(4,2)+H_invDmH(2,1)*X(2,2)*X(2,4)*H_invDmH(4,2)+H_invDm&
&H(1,1)*X(2,1)*X(2,4)*H_invDmH(4,2)+X(1,2)*X(1,4)*H_invDmH(2,1)*H_&
&invDmH(4,2)+H_invDmH(1,1)*X(1,1)*X(1,4)*H_invDmH(4,2)+H_invDmH(3,&
&2)*X(3,3)*X(3,4)*H_invDmH(4,1)+H_invDmH(2,2)*X(3,2)*X(3,4)*H_invD&
&mH(4,1)+H_invDmH(1,2)*X(3,1)*X(3,4)*H_invDmH(4,1)+X(2,3)*X(2,4)*H&
&_invDmH(3,2)*H_invDmH(4,1)+X(1,3)*X(1,4)*H_invDmH(3,2)*H_invDmH(4&
&,1)+H_invDmH(2,2)*X(2,2)*X(2,4)*H_invDmH(4,1)+H_invDmH(1,2)*X(2,1&
&)*X(2,4)*H_invDmH(4,1)+X(1,2)*X(1,4)*H_invDmH(2,2)*H_invDmH(4,1)+&
&X(1,1)*H_invDmH(1,2)*X(1,4)*H_invDmH(4,1)+H_invDmH(3,1)*H_invDmH(&
&3,2)*tt36+H_invDmH(2,1)*H_invDmH(3,2)*X(3,2)*X(3,3)+H_invDmH(2,2)&
&*H_invDmH(3,1)*X(3,2)*X(3,3)+H_invDmH(1,1)*X(3,1)*H_invDmH(3,2)*X&
&(3,3)+H_invDmH(1,2)*H_invDmH(3,1)*X(3,1)*X(3,3)+H_invDmH(2,1)*H_i&
&nvDmH(2,2)*tt35+H_invDmH(1,1)*H_invDmH(2,2)*X(3,1)*X(3,2)+H_invDm&
&H(1,2)*H_invDmH(2,1)*X(3,1)*X(3,2)+tt34*H_invDmH(3,1)*H_invDmH(3,&
&2)+tt33*H_invDmH(3,1)*H_invDmH(3,2)+H_invDmH(2,1)*X(2,2)*X(2,3)*H&
&_invDmH(3,2)+H_invDmH(1,1)*X(2,1)*X(2,3)*H_invDmH(3,2)+X(1,2)*X(1&
&,3)*H_invDmH(2,1)*H_invDmH(3,2)+H_invDmH(1,1)*X(1,1)*X(1,3)*H_inv&
&DmH(3,2)+H_invDmH(1,1)*H_invDmH(1,2)*tt32+H_invDmH(2,2)*X(2,2)*X(&
&2,3)*H_invDmH(3,1)+H_invDmH(1,2)*X(2,1)*X(2,3)*H_invDmH(3,1)+X(1,&
&2)*X(1,3)*H_invDmH(2,2)*H_invDmH(3,1)+X(1,1)*H_invDmH(1,2)*X(1,3)&
&*H_invDmH(3,1)+H_invDmH(2,1)*H_invDmH(2,2)*tt31+H_invDmH(1,1)*X(2&
&,1)*H_invDmH(2,2)*X(2,2)+H_invDmH(1,2)*H_invDmH(2,1)*X(2,1)*X(2,2&
&)+tt30*H_invDmH(2,1)*H_invDmH(2,2)+H_invDmH(1,1)*X(1,1)*X(1,2)*H_&
&invDmH(2,2)+H_invDmH(1,1)*H_invDmH(1,2)*tt29+X(1,1)*H_invDmH(1,2)&
&*X(1,2)*H_invDmH(2,1)+H_invDmH(1,1)*tt28*H_invDmH(1,2)
tt53 = 2*H_invDmH(1,1)*H_invDmH(1,2)*tt52
tt54 = H_invDmH(1,1)*X(1,8)*H_invDmH(8,2)+H_invDmH(1,2)*X(1,8)*H_&
&invDmH(8,1)+H_invDmH(1,1)*X(1,7)*H_invDmH(7,2)+H_invDmH(1,2)*X(1,&
&7)*H_invDmH(7,1)+H_invDmH(1,1)*X(1,6)*H_invDmH(6,2)+H_invDmH(1,2)&
&*X(1,6)*H_invDmH(6,1)+H_invDmH(1,1)*X(1,5)*H_invDmH(5,2)+H_invDmH&
&(1,2)*X(1,5)*H_invDmH(5,1)+H_invDmH(1,1)*X(1,4)*H_invDmH(4,2)+H_i&
&nvDmH(1,2)*X(1,4)*H_invDmH(4,1)+H_invDmH(1,1)*X(1,3)*H_invDmH(3,2&
&)+H_invDmH(1,2)*X(1,3)*H_invDmH(3,1)+H_invDmH(1,1)*X(1,2)*H_invDm&
&H(2,2)+H_invDmH(1,2)*X(1,2)*H_invDmH(2,1)+2*H_invDmH(1,1)*X(1,1)*&
&H_invDmH(1,2)
tt55 = tt2*tt17
tt56 = tt51*H_invDmH(8,1)*H_invDmH(8,3)+tt50*H_invDmH(8,1)*H_invD&
&mH(8,3)+tt49*H_invDmH(8,1)*H_invDmH(8,3)+X(3,7)*X(3,8)*H_invDmH(7&
&,1)*H_invDmH(8,3)+X(2,7)*X(2,8)*H_invDmH(7,1)*H_invDmH(8,3)+X(1,7&
&)*X(1,8)*H_invDmH(7,1)*H_invDmH(8,3)+X(3,6)*X(3,8)*H_invDmH(6,1)*&
&H_invDmH(8,3)+X(2,6)*X(2,8)*H_invDmH(6,1)*H_invDmH(8,3)+X(1,6)*X(&
&1,8)*H_invDmH(6,1)*H_invDmH(8,3)+X(3,5)*X(3,8)*H_invDmH(5,1)*H_in&
&vDmH(8,3)+X(2,5)*X(2,8)*H_invDmH(5,1)*H_invDmH(8,3)+X(1,5)*X(1,8)&
&*H_invDmH(5,1)*H_invDmH(8,3)+X(3,4)*X(3,8)*H_invDmH(4,1)*H_invDmH&
&(8,3)+X(2,4)*X(2,8)*H_invDmH(4,1)*H_invDmH(8,3)+X(1,4)*X(1,8)*H_i&
&nvDmH(4,1)*H_invDmH(8,3)+H_invDmH(3,1)*X(3,3)*X(3,8)*H_invDmH(8,3&
&)+H_invDmH(2,1)*X(3,2)*X(3,8)*H_invDmH(8,3)+H_invDmH(1,1)*X(3,1)*&
&X(3,8)*H_invDmH(8,3)+X(2,3)*X(2,8)*H_invDmH(3,1)*H_invDmH(8,3)+X(&
&1,3)*X(1,8)*H_invDmH(3,1)*H_invDmH(8,3)+H_invDmH(2,1)*X(2,2)*X(2,&
&8)*H_invDmH(8,3)+H_invDmH(1,1)*X(2,1)*X(2,8)*H_invDmH(8,3)+X(1,2)&
&*X(1,8)*H_invDmH(2,1)*H_invDmH(8,3)+H_invDmH(1,1)*X(1,1)*X(1,8)*H&
&_invDmH(8,3)+X(3,7)*X(3,8)*H_invDmH(7,3)*H_invDmH(8,1)+X(2,7)*X(2&
&,8)*H_invDmH(7,3)*H_invDmH(8,1)+X(1,7)*X(1,8)*H_invDmH(7,3)*H_inv&
&DmH(8,1)+X(3,6)*X(3,8)*H_invDmH(6,3)*H_invDmH(8,1)+X(2,6)*X(2,8)*&
&H_invDmH(6,3)*H_invDmH(8,1)+X(1,6)*X(1,8)*H_invDmH(6,3)*H_invDmH(&
&8,1)+X(3,5)*X(3,8)*H_invDmH(5,3)*H_invDmH(8,1)+X(2,5)*X(2,8)*H_in&
&vDmH(5,3)*H_invDmH(8,1)+X(1,5)*X(1,8)*H_invDmH(5,3)*H_invDmH(8,1)&
&+X(3,4)*X(3,8)*H_invDmH(4,3)*H_invDmH(8,1)+X(2,4)*X(2,8)*H_invDmH&
&(4,3)*H_invDmH(8,1)+X(1,4)*X(1,8)*H_invDmH(4,3)*H_invDmH(8,1)+H_i&
&nvDmH(3,3)*X(3,3)*X(3,8)*H_invDmH(8,1)+H_invDmH(2,3)*X(3,2)*X(3,8&
&)*H_invDmH(8,1)+H_invDmH(1,3)*X(3,1)*X(3,8)*H_invDmH(8,1)+X(2,3)*&
&X(2,8)*H_invDmH(3,3)*H_invDmH(8,1)+X(1,3)*X(1,8)*H_invDmH(3,3)*H_&
&invDmH(8,1)+X(2,2)*H_invDmH(2,3)*X(2,8)*H_invDmH(8,1)+H_invDmH(1,&
&3)*X(2,1)*X(2,8)*H_invDmH(8,1)+X(1,2)*X(1,8)*H_invDmH(2,3)*H_invD&
&mH(8,1)+X(1,1)*H_invDmH(1,3)*X(1,8)*H_invDmH(8,1)+tt48*H_invDmH(7&
&,1)*H_invDmH(7,3)+tt47*H_invDmH(7,1)*H_invDmH(7,3)+tt46*H_invDmH(&
&7,1)*H_invDmH(7,3)+X(3,6)*X(3,7)*H_invDmH(6,1)*H_invDmH(7,3)+X(2,&
&6)*X(2,7)*H_invDmH(6,1)*H_invDmH(7,3)+X(1,6)*X(1,7)*H_invDmH(6,1)&
&*H_invDmH(7,3)+X(3,5)*X(3,7)*H_invDmH(5,1)*H_invDmH(7,3)+X(2,5)*X&
&(2,7)*H_invDmH(5,1)*H_invDmH(7,3)+X(1,5)*X(1,7)*H_invDmH(5,1)*H_i&
&nvDmH(7,3)+X(3,4)*X(3,7)*H_invDmH(4,1)*H_invDmH(7,3)+X(2,4)*X(2,7&
&)*H_invDmH(4,1)*H_invDmH(7,3)+X(1,4)*X(1,7)*H_invDmH(4,1)*H_invDm&
&H(7,3)+H_invDmH(3,1)*X(3,3)*X(3,7)*H_invDmH(7,3)+H_invDmH(2,1)*X(&
&3,2)*X(3,7)*H_invDmH(7,3)+H_invDmH(1,1)*X(3,1)*X(3,7)*H_invDmH(7,&
&3)+X(2,3)*X(2,7)*H_invDmH(3,1)*H_invDmH(7,3)+X(1,3)*X(1,7)*H_invD&
&mH(3,1)*H_invDmH(7,3)+H_invDmH(2,1)*X(2,2)*X(2,7)*H_invDmH(7,3)+H&
&_invDmH(1,1)*X(2,1)*X(2,7)*H_invDmH(7,3)+X(1,2)*X(1,7)*H_invDmH(2&
&,1)*H_invDmH(7,3)+H_invDmH(1,1)*X(1,1)*X(1,7)*H_invDmH(7,3)+X(3,6&
&)*X(3,7)*H_invDmH(6,3)*H_invDmH(7,1)+X(2,6)*X(2,7)*H_invDmH(6,3)*&
&H_invDmH(7,1)+X(1,6)*X(1,7)*H_invDmH(6,3)*H_invDmH(7,1)+X(3,5)*X(&
&3,7)*H_invDmH(5,3)*H_invDmH(7,1)+X(2,5)*X(2,7)*H_invDmH(5,3)*H_in&
&vDmH(7,1)+X(1,5)*X(1,7)*H_invDmH(5,3)*H_invDmH(7,1)+X(3,4)*X(3,7)&
&*H_invDmH(4,3)*H_invDmH(7,1)+X(2,4)*X(2,7)*H_invDmH(4,3)*H_invDmH&
&(7,1)+X(1,4)*X(1,7)*H_invDmH(4,3)*H_invDmH(7,1)+H_invDmH(3,3)*X(3&
&,3)*X(3,7)*H_invDmH(7,1)+H_invDmH(2,3)*X(3,2)*X(3,7)*H_invDmH(7,1&
&)+H_invDmH(1,3)*X(3,1)*X(3,7)*H_invDmH(7,1)+X(2,3)*X(2,7)*H_invDm&
&H(3,3)*H_invDmH(7,1)+X(1,3)*X(1,7)*H_invDmH(3,3)*H_invDmH(7,1)+X(&
&2,2)*H_invDmH(2,3)*X(2,7)*H_invDmH(7,1)+H_invDmH(1,3)*X(2,1)*X(2,&
&7)*H_invDmH(7,1)+X(1,2)*X(1,7)*H_invDmH(2,3)*H_invDmH(7,1)+X(1,1)&
&*H_invDmH(1,3)*X(1,7)*H_invDmH(7,1)+tt45*H_invDmH(6,1)*H_invDmH(6&
&,3)+tt44*H_invDmH(6,1)*H_invDmH(6,3)+tt43*H_invDmH(6,1)*H_invDmH(&
&6,3)+X(3,5)*X(3,6)*H_invDmH(5,1)*H_invDmH(6,3)+X(2,5)*X(2,6)*H_in&
&vDmH(5,1)*H_invDmH(6,3)+X(1,5)*X(1,6)*H_invDmH(5,1)*H_invDmH(6,3)&
&+X(3,4)*X(3,6)*H_invDmH(4,1)*H_invDmH(6,3)+X(2,4)*X(2,6)*H_invDmH&
&(4,1)*H_invDmH(6,3)+X(1,4)*X(1,6)*H_invDmH(4,1)*H_invDmH(6,3)+H_i&
&nvDmH(3,1)*X(3,3)*X(3,6)*H_invDmH(6,3)+H_invDmH(2,1)*X(3,2)*X(3,6&
&)*H_invDmH(6,3)+H_invDmH(1,1)*X(3,1)*X(3,6)*H_invDmH(6,3)+X(2,3)*&
&X(2,6)*H_invDmH(3,1)*H_invDmH(6,3)+X(1,3)*X(1,6)*H_invDmH(3,1)*H_&
&invDmH(6,3)+H_invDmH(2,1)*X(2,2)*X(2,6)*H_invDmH(6,3)+H_invDmH(1,&
&1)*X(2,1)*X(2,6)*H_invDmH(6,3)+X(1,2)*X(1,6)*H_invDmH(2,1)*H_invD&
&mH(6,3)+H_invDmH(1,1)*X(1,1)*X(1,6)*H_invDmH(6,3)+X(3,5)*X(3,6)*H&
&_invDmH(5,3)*H_invDmH(6,1)+X(2,5)*X(2,6)*H_invDmH(5,3)*H_invDmH(6&
&,1)+X(1,5)*X(1,6)*H_invDmH(5,3)*H_invDmH(6,1)+X(3,4)*X(3,6)*H_inv&
&DmH(4,3)*H_invDmH(6,1)+X(2,4)*X(2,6)*H_invDmH(4,3)*H_invDmH(6,1)+&
&X(1,4)*X(1,6)*H_invDmH(4,3)*H_invDmH(6,1)+H_invDmH(3,3)*X(3,3)*X(&
&3,6)*H_invDmH(6,1)+H_invDmH(2,3)*X(3,2)*X(3,6)*H_invDmH(6,1)+H_in&
&vDmH(1,3)*X(3,1)*X(3,6)*H_invDmH(6,1)+X(2,3)*X(2,6)*H_invDmH(3,3)&
&*H_invDmH(6,1)+X(1,3)*X(1,6)*H_invDmH(3,3)*H_invDmH(6,1)+X(2,2)*H&
&_invDmH(2,3)*X(2,6)*H_invDmH(6,1)+H_invDmH(1,3)*X(2,1)*X(2,6)*H_i&
&nvDmH(6,1)+X(1,2)*X(1,6)*H_invDmH(2,3)*H_invDmH(6,1)+X(1,1)*H_inv&
&DmH(1,3)*X(1,6)*H_invDmH(6,1)+tt42*H_invDmH(5,1)*H_invDmH(5,3)+tt&
&41*H_invDmH(5,1)*H_invDmH(5,3)+tt40*H_invDmH(5,1)*H_invDmH(5,3)+X&
&(3,4)*X(3,5)*H_invDmH(4,1)*H_invDmH(5,3)+X(2,4)*X(2,5)*H_invDmH(4&
&,1)*H_invDmH(5,3)+X(1,4)*X(1,5)*H_invDmH(4,1)*H_invDmH(5,3)+H_inv&
&DmH(3,1)*X(3,3)*X(3,5)*H_invDmH(5,3)+H_invDmH(2,1)*X(3,2)*X(3,5)*&
&H_invDmH(5,3)+H_invDmH(1,1)*X(3,1)*X(3,5)*H_invDmH(5,3)+X(2,3)*X(&
&2,5)*H_invDmH(3,1)*H_invDmH(5,3)+X(1,3)*X(1,5)*H_invDmH(3,1)*H_in&
&vDmH(5,3)+H_invDmH(2,1)*X(2,2)*X(2,5)*H_invDmH(5,3)+H_invDmH(1,1)&
&*X(2,1)*X(2,5)*H_invDmH(5,3)+X(1,2)*X(1,5)*H_invDmH(2,1)*H_invDmH&
&(5,3)+H_invDmH(1,1)*X(1,1)*X(1,5)*H_invDmH(5,3)+X(3,4)*X(3,5)*H_i&
&nvDmH(4,3)*H_invDmH(5,1)+X(2,4)*X(2,5)*H_invDmH(4,3)*H_invDmH(5,1&
&)+X(1,4)*X(1,5)*H_invDmH(4,3)*H_invDmH(5,1)+H_invDmH(3,3)*X(3,3)*&
&X(3,5)*H_invDmH(5,1)+H_invDmH(2,3)*X(3,2)*X(3,5)*H_invDmH(5,1)+H_&
&invDmH(1,3)*X(3,1)*X(3,5)*H_invDmH(5,1)+X(2,3)*X(2,5)*H_invDmH(3,&
&3)*H_invDmH(5,1)+X(1,3)*X(1,5)*H_invDmH(3,3)*H_invDmH(5,1)+X(2,2)&
&*H_invDmH(2,3)*X(2,5)*H_invDmH(5,1)+H_invDmH(1,3)*X(2,1)*X(2,5)*H&
&_invDmH(5,1)+X(1,2)*X(1,5)*H_invDmH(2,3)*H_invDmH(5,1)+X(1,1)*H_i&
&nvDmH(1,3)*X(1,5)*H_invDmH(5,1)+tt39*H_invDmH(4,1)*H_invDmH(4,3)+&
&tt38*H_invDmH(4,1)*H_invDmH(4,3)+tt37*H_invDmH(4,1)*H_invDmH(4,3)&
&+H_invDmH(3,1)*X(3,3)*X(3,4)*H_invDmH(4,3)+H_invDmH(2,1)*X(3,2)*X&
&(3,4)*H_invDmH(4,3)+H_invDmH(1,1)*X(3,1)*X(3,4)*H_invDmH(4,3)+X(2&
&,3)*X(2,4)*H_invDmH(3,1)*H_invDmH(4,3)+X(1,3)*X(1,4)*H_invDmH(3,1&
&)*H_invDmH(4,3)+H_invDmH(2,1)*X(2,2)*X(2,4)*H_invDmH(4,3)+H_invDm&
&H(1,1)*X(2,1)*X(2,4)*H_invDmH(4,3)+X(1,2)*X(1,4)*H_invDmH(2,1)*H_&
&invDmH(4,3)+H_invDmH(1,1)*X(1,1)*X(1,4)*H_invDmH(4,3)+H_invDmH(3,&
&3)*X(3,3)*X(3,4)*H_invDmH(4,1)+H_invDmH(2,3)*X(3,2)*X(3,4)*H_invD&
&mH(4,1)+H_invDmH(1,3)*X(3,1)*X(3,4)*H_invDmH(4,1)+X(2,3)*X(2,4)*H&
&_invDmH(3,3)*H_invDmH(4,1)+X(1,3)*X(1,4)*H_invDmH(3,3)*H_invDmH(4&
&,1)+X(2,2)*H_invDmH(2,3)*X(2,4)*H_invDmH(4,1)+H_invDmH(1,3)*X(2,1&
&)*X(2,4)*H_invDmH(4,1)+X(1,2)*X(1,4)*H_invDmH(2,3)*H_invDmH(4,1)+&
&X(1,1)*H_invDmH(1,3)*X(1,4)*H_invDmH(4,1)+H_invDmH(3,1)*H_invDmH(&
&3,3)*tt36+H_invDmH(2,1)*X(3,2)*H_invDmH(3,3)*X(3,3)+H_invDmH(1,1)&
&*X(3,1)*H_invDmH(3,3)*X(3,3)+H_invDmH(2,3)*H_invDmH(3,1)*X(3,2)*X&
&(3,3)+H_invDmH(1,3)*H_invDmH(3,1)*X(3,1)*X(3,3)+tt34*H_invDmH(3,1&
&)*H_invDmH(3,3)+tt33*H_invDmH(3,1)*H_invDmH(3,3)+H_invDmH(2,1)*X(&
&2,2)*X(2,3)*H_invDmH(3,3)+H_invDmH(1,1)*X(2,1)*X(2,3)*H_invDmH(3,&
&3)+X(1,2)*X(1,3)*H_invDmH(2,1)*H_invDmH(3,3)+H_invDmH(1,1)*X(1,1)&
&*X(1,3)*H_invDmH(3,3)+H_invDmH(2,1)*H_invDmH(2,3)*tt35+H_invDmH(1&
&,1)*H_invDmH(2,3)*X(3,1)*X(3,2)+H_invDmH(1,3)*H_invDmH(2,1)*X(3,1&
&)*X(3,2)+H_invDmH(1,1)*H_invDmH(1,3)*tt32+X(2,2)*H_invDmH(2,3)*X(&
&2,3)*H_invDmH(3,1)+H_invDmH(1,3)*X(2,1)*X(2,3)*H_invDmH(3,1)+X(1,&
&2)*X(1,3)*H_invDmH(2,3)*H_invDmH(3,1)+X(1,1)*H_invDmH(1,3)*X(1,3)&
&*H_invDmH(3,1)+H_invDmH(2,1)*tt31*H_invDmH(2,3)+H_invDmH(1,1)*X(2&
&,1)*X(2,2)*H_invDmH(2,3)+tt30*H_invDmH(2,1)*H_invDmH(2,3)+H_invDm&
&H(1,1)*X(1,1)*X(1,2)*H_invDmH(2,3)+H_invDmH(1,3)*H_invDmH(2,1)*X(&
&2,1)*X(2,2)+H_invDmH(1,1)*H_invDmH(1,3)*tt29+X(1,1)*X(1,2)*H_invD&
&mH(1,3)*H_invDmH(2,1)+H_invDmH(1,1)*tt28*H_invDmH(1,3)
tt57 = 2*H_invDmH(1,1)*H_invDmH(1,3)*tt56
tt58 = tt51*H_invDmH(8,2)*H_invDmH(8,3)+tt50*H_invDmH(8,2)*H_invD&
&mH(8,3)+tt49*H_invDmH(8,2)*H_invDmH(8,3)+X(3,7)*X(3,8)*H_invDmH(7&
&,2)*H_invDmH(8,3)+X(2,7)*X(2,8)*H_invDmH(7,2)*H_invDmH(8,3)+X(1,7&
&)*X(1,8)*H_invDmH(7,2)*H_invDmH(8,3)+X(3,6)*X(3,8)*H_invDmH(6,2)*&
&H_invDmH(8,3)+X(2,6)*X(2,8)*H_invDmH(6,2)*H_invDmH(8,3)+X(1,6)*X(&
&1,8)*H_invDmH(6,2)*H_invDmH(8,3)+X(3,5)*X(3,8)*H_invDmH(5,2)*H_in&
&vDmH(8,3)+X(2,5)*X(2,8)*H_invDmH(5,2)*H_invDmH(8,3)+X(1,5)*X(1,8)&
&*H_invDmH(5,2)*H_invDmH(8,3)+X(3,4)*X(3,8)*H_invDmH(4,2)*H_invDmH&
&(8,3)+X(2,4)*X(2,8)*H_invDmH(4,2)*H_invDmH(8,3)+X(1,4)*X(1,8)*H_i&
&nvDmH(4,2)*H_invDmH(8,3)+H_invDmH(3,2)*X(3,3)*X(3,8)*H_invDmH(8,3&
&)+H_invDmH(2,2)*X(3,2)*X(3,8)*H_invDmH(8,3)+H_invDmH(1,2)*X(3,1)*&
&X(3,8)*H_invDmH(8,3)+X(2,3)*X(2,8)*H_invDmH(3,2)*H_invDmH(8,3)+X(&
&1,3)*X(1,8)*H_invDmH(3,2)*H_invDmH(8,3)+H_invDmH(2,2)*X(2,2)*X(2,&
&8)*H_invDmH(8,3)+H_invDmH(1,2)*X(2,1)*X(2,8)*H_invDmH(8,3)+X(1,2)&
&*X(1,8)*H_invDmH(2,2)*H_invDmH(8,3)+X(1,1)*H_invDmH(1,2)*X(1,8)*H&
&_invDmH(8,3)+X(3,7)*X(3,8)*H_invDmH(7,3)*H_invDmH(8,2)+X(2,7)*X(2&
&,8)*H_invDmH(7,3)*H_invDmH(8,2)+X(1,7)*X(1,8)*H_invDmH(7,3)*H_inv&
&DmH(8,2)+X(3,6)*X(3,8)*H_invDmH(6,3)*H_invDmH(8,2)+X(2,6)*X(2,8)*&
&H_invDmH(6,3)*H_invDmH(8,2)+X(1,6)*X(1,8)*H_invDmH(6,3)*H_invDmH(&
&8,2)+X(3,5)*X(3,8)*H_invDmH(5,3)*H_invDmH(8,2)+X(2,5)*X(2,8)*H_in&
&vDmH(5,3)*H_invDmH(8,2)+X(1,5)*X(1,8)*H_invDmH(5,3)*H_invDmH(8,2)&
&+X(3,4)*X(3,8)*H_invDmH(4,3)*H_invDmH(8,2)+X(2,4)*X(2,8)*H_invDmH&
&(4,3)*H_invDmH(8,2)+X(1,4)*X(1,8)*H_invDmH(4,3)*H_invDmH(8,2)+H_i&
&nvDmH(3,3)*X(3,3)*X(3,8)*H_invDmH(8,2)+H_invDmH(2,3)*X(3,2)*X(3,8&
&)*H_invDmH(8,2)+H_invDmH(1,3)*X(3,1)*X(3,8)*H_invDmH(8,2)+X(2,3)*&
&X(2,8)*H_invDmH(3,3)*H_invDmH(8,2)+X(1,3)*X(1,8)*H_invDmH(3,3)*H_&
&invDmH(8,2)+X(2,2)*H_invDmH(2,3)*X(2,8)*H_invDmH(8,2)+H_invDmH(1,&
&3)*X(2,1)*X(2,8)*H_invDmH(8,2)+X(1,2)*X(1,8)*H_invDmH(2,3)*H_invD&
&mH(8,2)+X(1,1)*H_invDmH(1,3)*X(1,8)*H_invDmH(8,2)+tt48*H_invDmH(7&
&,2)*H_invDmH(7,3)+tt47*H_invDmH(7,2)*H_invDmH(7,3)+tt46*H_invDmH(&
&7,2)*H_invDmH(7,3)+X(3,6)*X(3,7)*H_invDmH(6,2)*H_invDmH(7,3)+X(2,&
&6)*X(2,7)*H_invDmH(6,2)*H_invDmH(7,3)+X(1,6)*X(1,7)*H_invDmH(6,2)&
&*H_invDmH(7,3)+X(3,5)*X(3,7)*H_invDmH(5,2)*H_invDmH(7,3)+X(2,5)*X&
&(2,7)*H_invDmH(5,2)*H_invDmH(7,3)+X(1,5)*X(1,7)*H_invDmH(5,2)*H_i&
&nvDmH(7,3)+X(3,4)*X(3,7)*H_invDmH(4,2)*H_invDmH(7,3)+X(2,4)*X(2,7&
&)*H_invDmH(4,2)*H_invDmH(7,3)+X(1,4)*X(1,7)*H_invDmH(4,2)*H_invDm&
&H(7,3)+H_invDmH(3,2)*X(3,3)*X(3,7)*H_invDmH(7,3)+H_invDmH(2,2)*X(&
&3,2)*X(3,7)*H_invDmH(7,3)+H_invDmH(1,2)*X(3,1)*X(3,7)*H_invDmH(7,&
&3)+X(2,3)*X(2,7)*H_invDmH(3,2)*H_invDmH(7,3)+X(1,3)*X(1,7)*H_invD&
&mH(3,2)*H_invDmH(7,3)+H_invDmH(2,2)*X(2,2)*X(2,7)*H_invDmH(7,3)+H&
&_invDmH(1,2)*X(2,1)*X(2,7)*H_invDmH(7,3)+X(1,2)*X(1,7)*H_invDmH(2&
&,2)*H_invDmH(7,3)+X(1,1)*H_invDmH(1,2)*X(1,7)*H_invDmH(7,3)+X(3,6&
&)*X(3,7)*H_invDmH(6,3)*H_invDmH(7,2)+X(2,6)*X(2,7)*H_invDmH(6,3)*&
&H_invDmH(7,2)+X(1,6)*X(1,7)*H_invDmH(6,3)*H_invDmH(7,2)+X(3,5)*X(&
&3,7)*H_invDmH(5,3)*H_invDmH(7,2)+X(2,5)*X(2,7)*H_invDmH(5,3)*H_in&
&vDmH(7,2)+X(1,5)*X(1,7)*H_invDmH(5,3)*H_invDmH(7,2)+X(3,4)*X(3,7)&
&*H_invDmH(4,3)*H_invDmH(7,2)+X(2,4)*X(2,7)*H_invDmH(4,3)*H_invDmH&
&(7,2)+X(1,4)*X(1,7)*H_invDmH(4,3)*H_invDmH(7,2)+H_invDmH(3,3)*X(3&
&,3)*X(3,7)*H_invDmH(7,2)+H_invDmH(2,3)*X(3,2)*X(3,7)*H_invDmH(7,2&
&)+H_invDmH(1,3)*X(3,1)*X(3,7)*H_invDmH(7,2)+X(2,3)*X(2,7)*H_invDm&
&H(3,3)*H_invDmH(7,2)+X(1,3)*X(1,7)*H_invDmH(3,3)*H_invDmH(7,2)+X(&
&2,2)*H_invDmH(2,3)*X(2,7)*H_invDmH(7,2)+H_invDmH(1,3)*X(2,1)*X(2,&
&7)*H_invDmH(7,2)+X(1,2)*X(1,7)*H_invDmH(2,3)*H_invDmH(7,2)+X(1,1)&
&*H_invDmH(1,3)*X(1,7)*H_invDmH(7,2)+tt45*H_invDmH(6,2)*H_invDmH(6&
&,3)+tt44*H_invDmH(6,2)*H_invDmH(6,3)+tt43*H_invDmH(6,2)*H_invDmH(&
&6,3)+X(3,5)*X(3,6)*H_invDmH(5,2)*H_invDmH(6,3)+X(2,5)*X(2,6)*H_in&
&vDmH(5,2)*H_invDmH(6,3)+X(1,5)*X(1,6)*H_invDmH(5,2)*H_invDmH(6,3)&
&+X(3,4)*X(3,6)*H_invDmH(4,2)*H_invDmH(6,3)+X(2,4)*X(2,6)*H_invDmH&
&(4,2)*H_invDmH(6,3)+X(1,4)*X(1,6)*H_invDmH(4,2)*H_invDmH(6,3)+H_i&
&nvDmH(3,2)*X(3,3)*X(3,6)*H_invDmH(6,3)+H_invDmH(2,2)*X(3,2)*X(3,6&
&)*H_invDmH(6,3)+H_invDmH(1,2)*X(3,1)*X(3,6)*H_invDmH(6,3)+X(2,3)*&
&X(2,6)*H_invDmH(3,2)*H_invDmH(6,3)+X(1,3)*X(1,6)*H_invDmH(3,2)*H_&
&invDmH(6,3)+H_invDmH(2,2)*X(2,2)*X(2,6)*H_invDmH(6,3)+H_invDmH(1,&
&2)*X(2,1)*X(2,6)*H_invDmH(6,3)+X(1,2)*X(1,6)*H_invDmH(2,2)*H_invD&
&mH(6,3)+X(1,1)*H_invDmH(1,2)*X(1,6)*H_invDmH(6,3)+X(3,5)*X(3,6)*H&
&_invDmH(5,3)*H_invDmH(6,2)+X(2,5)*X(2,6)*H_invDmH(5,3)*H_invDmH(6&
&,2)+X(1,5)*X(1,6)*H_invDmH(5,3)*H_invDmH(6,2)+X(3,4)*X(3,6)*H_inv&
&DmH(4,3)*H_invDmH(6,2)+X(2,4)*X(2,6)*H_invDmH(4,3)*H_invDmH(6,2)+&
&X(1,4)*X(1,6)*H_invDmH(4,3)*H_invDmH(6,2)+H_invDmH(3,3)*X(3,3)*X(&
&3,6)*H_invDmH(6,2)+H_invDmH(2,3)*X(3,2)*X(3,6)*H_invDmH(6,2)+H_in&
&vDmH(1,3)*X(3,1)*X(3,6)*H_invDmH(6,2)+X(2,3)*X(2,6)*H_invDmH(3,3)&
&*H_invDmH(6,2)+X(1,3)*X(1,6)*H_invDmH(3,3)*H_invDmH(6,2)+X(2,2)*H&
&_invDmH(2,3)*X(2,6)*H_invDmH(6,2)+H_invDmH(1,3)*X(2,1)*X(2,6)*H_i&
&nvDmH(6,2)+X(1,2)*X(1,6)*H_invDmH(2,3)*H_invDmH(6,2)+X(1,1)*H_inv&
&DmH(1,3)*X(1,6)*H_invDmH(6,2)+tt42*H_invDmH(5,2)*H_invDmH(5,3)+tt&
&41*H_invDmH(5,2)*H_invDmH(5,3)+tt40*H_invDmH(5,2)*H_invDmH(5,3)+X&
&(3,4)*X(3,5)*H_invDmH(4,2)*H_invDmH(5,3)+X(2,4)*X(2,5)*H_invDmH(4&
&,2)*H_invDmH(5,3)+X(1,4)*X(1,5)*H_invDmH(4,2)*H_invDmH(5,3)+H_inv&
&DmH(3,2)*X(3,3)*X(3,5)*H_invDmH(5,3)+H_invDmH(2,2)*X(3,2)*X(3,5)*&
&H_invDmH(5,3)+H_invDmH(1,2)*X(3,1)*X(3,5)*H_invDmH(5,3)+X(2,3)*X(&
&2,5)*H_invDmH(3,2)*H_invDmH(5,3)+X(1,3)*X(1,5)*H_invDmH(3,2)*H_in&
&vDmH(5,3)+H_invDmH(2,2)*X(2,2)*X(2,5)*H_invDmH(5,3)+H_invDmH(1,2)&
&*X(2,1)*X(2,5)*H_invDmH(5,3)+X(1,2)*X(1,5)*H_invDmH(2,2)*H_invDmH&
&(5,3)+X(1,1)*H_invDmH(1,2)*X(1,5)*H_invDmH(5,3)+X(3,4)*X(3,5)*H_i&
&nvDmH(4,3)*H_invDmH(5,2)+X(2,4)*X(2,5)*H_invDmH(4,3)*H_invDmH(5,2&
&)+X(1,4)*X(1,5)*H_invDmH(4,3)*H_invDmH(5,2)+H_invDmH(3,3)*X(3,3)*&
&X(3,5)*H_invDmH(5,2)+H_invDmH(2,3)*X(3,2)*X(3,5)*H_invDmH(5,2)+H_&
&invDmH(1,3)*X(3,1)*X(3,5)*H_invDmH(5,2)+X(2,3)*X(2,5)*H_invDmH(3,&
&3)*H_invDmH(5,2)+X(1,3)*X(1,5)*H_invDmH(3,3)*H_invDmH(5,2)+X(2,2)&
&*H_invDmH(2,3)*X(2,5)*H_invDmH(5,2)+H_invDmH(1,3)*X(2,1)*X(2,5)*H&
&_invDmH(5,2)+X(1,2)*X(1,5)*H_invDmH(2,3)*H_invDmH(5,2)+X(1,1)*H_i&
&nvDmH(1,3)*X(1,5)*H_invDmH(5,2)+tt39*H_invDmH(4,2)*H_invDmH(4,3)+&
&tt38*H_invDmH(4,2)*H_invDmH(4,3)+tt37*H_invDmH(4,2)*H_invDmH(4,3)&
&+H_invDmH(3,2)*X(3,3)*X(3,4)*H_invDmH(4,3)+H_invDmH(2,2)*X(3,2)*X&
&(3,4)*H_invDmH(4,3)+H_invDmH(1,2)*X(3,1)*X(3,4)*H_invDmH(4,3)+X(2&
&,3)*X(2,4)*H_invDmH(3,2)*H_invDmH(4,3)+X(1,3)*X(1,4)*H_invDmH(3,2&
&)*H_invDmH(4,3)+H_invDmH(2,2)*X(2,2)*X(2,4)*H_invDmH(4,3)+H_invDm&
&H(1,2)*X(2,1)*X(2,4)*H_invDmH(4,3)+X(1,2)*X(1,4)*H_invDmH(2,2)*H_&
&invDmH(4,3)+X(1,1)*H_invDmH(1,2)*X(1,4)*H_invDmH(4,3)+H_invDmH(3,&
&3)*X(3,3)*X(3,4)*H_invDmH(4,2)+H_invDmH(2,3)*X(3,2)*X(3,4)*H_invD&
&mH(4,2)+H_invDmH(1,3)*X(3,1)*X(3,4)*H_invDmH(4,2)+X(2,3)*X(2,4)*H&
&_invDmH(3,3)*H_invDmH(4,2)+X(1,3)*X(1,4)*H_invDmH(3,3)*H_invDmH(4&
&,2)+X(2,2)*H_invDmH(2,3)*X(2,4)*H_invDmH(4,2)+H_invDmH(1,3)*X(2,1&
&)*X(2,4)*H_invDmH(4,2)+X(1,2)*X(1,4)*H_invDmH(2,3)*H_invDmH(4,2)+&
&X(1,1)*H_invDmH(1,3)*X(1,4)*H_invDmH(4,2)+H_invDmH(3,2)*H_invDmH(&
&3,3)*tt36+H_invDmH(2,2)*X(3,2)*H_invDmH(3,3)*X(3,3)+H_invDmH(1,2)&
&*X(3,1)*H_invDmH(3,3)*X(3,3)+H_invDmH(2,3)*H_invDmH(3,2)*X(3,2)*X&
&(3,3)+H_invDmH(1,3)*X(3,1)*H_invDmH(3,2)*X(3,3)+tt34*H_invDmH(3,2&
&)*H_invDmH(3,3)+tt33*H_invDmH(3,2)*H_invDmH(3,3)+H_invDmH(2,2)*X(&
&2,2)*X(2,3)*H_invDmH(3,3)+H_invDmH(1,2)*X(2,1)*X(2,3)*H_invDmH(3,&
&3)+X(1,2)*X(1,3)*H_invDmH(2,2)*H_invDmH(3,3)+X(1,1)*H_invDmH(1,2)&
&*X(1,3)*H_invDmH(3,3)+H_invDmH(2,2)*H_invDmH(2,3)*tt35+H_invDmH(1&
&,2)*H_invDmH(2,3)*X(3,1)*X(3,2)+H_invDmH(1,3)*H_invDmH(2,2)*X(3,1&
&)*X(3,2)+X(2,2)*H_invDmH(2,3)*X(2,3)*H_invDmH(3,2)+H_invDmH(1,3)*&
&X(2,1)*X(2,3)*H_invDmH(3,2)+X(1,2)*X(1,3)*H_invDmH(2,3)*H_invDmH(&
&3,2)+X(1,1)*H_invDmH(1,3)*X(1,3)*H_invDmH(3,2)+H_invDmH(1,2)*H_in&
&vDmH(1,3)*tt32+H_invDmH(2,2)*tt31*H_invDmH(2,3)+H_invDmH(1,2)*X(2&
&,1)*X(2,2)*H_invDmH(2,3)+tt30*H_invDmH(2,2)*H_invDmH(2,3)+X(1,1)*&
&H_invDmH(1,2)*X(1,2)*H_invDmH(2,3)+H_invDmH(1,3)*X(2,1)*H_invDmH(&
&2,2)*X(2,2)+X(1,1)*X(1,2)*H_invDmH(1,3)*H_invDmH(2,2)+H_invDmH(1,&
&2)*H_invDmH(1,3)*tt29+tt28*H_invDmH(1,2)*H_invDmH(1,3)
tt59 = 2*H_invDmH(1,2)*H_invDmH(1,3)*tt58
tt60 = H_invDmH(1,1)*X(1,8)*H_invDmH(8,3)+H_invDmH(1,3)*X(1,8)*H_&
&invDmH(8,1)+H_invDmH(1,1)*X(1,7)*H_invDmH(7,3)+H_invDmH(1,3)*X(1,&
&7)*H_invDmH(7,1)+H_invDmH(1,1)*X(1,6)*H_invDmH(6,3)+H_invDmH(1,3)&
&*X(1,6)*H_invDmH(6,1)+H_invDmH(1,1)*X(1,5)*H_invDmH(5,3)+H_invDmH&
&(1,3)*X(1,5)*H_invDmH(5,1)+H_invDmH(1,1)*X(1,4)*H_invDmH(4,3)+H_i&
&nvDmH(1,3)*X(1,4)*H_invDmH(4,1)+H_invDmH(1,1)*X(1,3)*H_invDmH(3,3&
&)+H_invDmH(1,3)*X(1,3)*H_invDmH(3,1)+H_invDmH(1,1)*X(1,2)*H_invDm&
&H(2,3)+X(1,2)*H_invDmH(1,3)*H_invDmH(2,1)+2*H_invDmH(1,1)*X(1,1)*&
&H_invDmH(1,3)
tt61 = H_invDmH(1,2)*X(1,8)*H_invDmH(8,3)+H_invDmH(1,3)*X(1,8)*H_&
&invDmH(8,2)+H_invDmH(1,2)*X(1,7)*H_invDmH(7,3)+H_invDmH(1,3)*X(1,&
&7)*H_invDmH(7,2)+H_invDmH(1,2)*X(1,6)*H_invDmH(6,3)+H_invDmH(1,3)&
&*X(1,6)*H_invDmH(6,2)+H_invDmH(1,2)*X(1,5)*H_invDmH(5,3)+H_invDmH&
&(1,3)*X(1,5)*H_invDmH(5,2)+H_invDmH(1,2)*X(1,4)*H_invDmH(4,3)+H_i&
&nvDmH(1,3)*X(1,4)*H_invDmH(4,2)+H_invDmH(1,2)*X(1,3)*H_invDmH(3,3&
&)+H_invDmH(1,3)*X(1,3)*H_invDmH(3,2)+H_invDmH(1,2)*X(1,2)*H_invDm&
&H(2,3)+X(1,2)*H_invDmH(1,3)*H_invDmH(2,2)+2*X(1,1)*H_invDmH(1,2)*&
&H_invDmH(1,3)
tt62 = tt3*tt24
tt63 = H_invDmH(1,3)*tt18+H_invDmH(1,2)*tt11+H_invDmH(1,1)*tt4
tt64 = H_invDmH(1,3)*tt20+H_invDmH(1,2)*tt13+H_invDmH(1,1)*tt6
tt65 = H_invDmH(1,1)*X(2,8)*H_invDmH(8,2)+H_invDmH(1,2)*X(2,8)*H_&
&invDmH(8,1)+H_invDmH(1,1)*X(2,7)*H_invDmH(7,2)+H_invDmH(1,2)*X(2,&
&7)*H_invDmH(7,1)+H_invDmH(1,1)*X(2,6)*H_invDmH(6,2)+H_invDmH(1,2)&
&*X(2,6)*H_invDmH(6,1)+H_invDmH(1,1)*X(2,5)*H_invDmH(5,2)+H_invDmH&
&(1,2)*X(2,5)*H_invDmH(5,1)+H_invDmH(1,1)*X(2,4)*H_invDmH(4,2)+H_i&
&nvDmH(1,2)*X(2,4)*H_invDmH(4,1)+H_invDmH(1,1)*X(2,3)*H_invDmH(3,2&
&)+H_invDmH(1,2)*X(2,3)*H_invDmH(3,1)+H_invDmH(1,1)*H_invDmH(2,2)*&
&X(2,2)+H_invDmH(1,2)*H_invDmH(2,1)*X(2,2)+2*H_invDmH(1,1)*H_invDm&
&H(1,2)*X(2,1)
tt66 = H_invDmH(1,1)*X(2,8)*H_invDmH(8,3)+H_invDmH(1,3)*X(2,8)*H_&
&invDmH(8,1)+H_invDmH(1,1)*X(2,7)*H_invDmH(7,3)+H_invDmH(1,3)*X(2,&
&7)*H_invDmH(7,1)+H_invDmH(1,1)*X(2,6)*H_invDmH(6,3)+H_invDmH(1,3)&
&*X(2,6)*H_invDmH(6,1)+H_invDmH(1,1)*X(2,5)*H_invDmH(5,3)+H_invDmH&
&(1,3)*X(2,5)*H_invDmH(5,1)+H_invDmH(1,1)*X(2,4)*H_invDmH(4,3)+H_i&
&nvDmH(1,3)*X(2,4)*H_invDmH(4,1)+H_invDmH(1,1)*X(2,3)*H_invDmH(3,3&
&)+H_invDmH(1,3)*X(2,3)*H_invDmH(3,1)+H_invDmH(1,1)*X(2,2)*H_invDm&
&H(2,3)+H_invDmH(1,3)*H_invDmH(2,1)*X(2,2)+2*H_invDmH(1,1)*H_invDm&
&H(1,3)*X(2,1)
tt67 = H_invDmH(1,2)*X(2,8)*H_invDmH(8,3)+H_invDmH(1,3)*X(2,8)*H_&
&invDmH(8,2)+H_invDmH(1,2)*X(2,7)*H_invDmH(7,3)+H_invDmH(1,3)*X(2,&
&7)*H_invDmH(7,2)+H_invDmH(1,2)*X(2,6)*H_invDmH(6,3)+H_invDmH(1,3)&
&*X(2,6)*H_invDmH(6,2)+H_invDmH(1,2)*X(2,5)*H_invDmH(5,3)+H_invDmH&
&(1,3)*X(2,5)*H_invDmH(5,2)+H_invDmH(1,2)*X(2,4)*H_invDmH(4,3)+H_i&
&nvDmH(1,3)*X(2,4)*H_invDmH(4,2)+H_invDmH(1,2)*X(2,3)*H_invDmH(3,3&
&)+H_invDmH(1,3)*X(2,3)*H_invDmH(3,2)+H_invDmH(1,2)*X(2,2)*H_invDm&
&H(2,3)+H_invDmH(1,3)*H_invDmH(2,2)*X(2,2)+2*H_invDmH(1,2)*H_invDm&
&H(1,3)*X(2,1)
tt68 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt61*tt67+tt60*tt66+2*tt3*tt&
&18*tt20+tt54*tt65+2*tt2*tt11*tt13+2*tt1*tt4*tt6)+lam(1,1)*tt63*tt&
&64)
tt69 = H_invDmH(1,3)*tt22+H_invDmH(1,2)*tt15+H_invDmH(1,1)*tt8
tt70 = H_invDmH(1,1)*X(3,8)*H_invDmH(8,2)+H_invDmH(1,2)*X(3,8)*H_&
&invDmH(8,1)+H_invDmH(1,1)*X(3,7)*H_invDmH(7,2)+H_invDmH(1,2)*X(3,&
&7)*H_invDmH(7,1)+H_invDmH(1,1)*X(3,6)*H_invDmH(6,2)+H_invDmH(1,2)&
&*X(3,6)*H_invDmH(6,1)+H_invDmH(1,1)*X(3,5)*H_invDmH(5,2)+H_invDmH&
&(1,2)*X(3,5)*H_invDmH(5,1)+H_invDmH(1,1)*X(3,4)*H_invDmH(4,2)+H_i&
&nvDmH(1,2)*X(3,4)*H_invDmH(4,1)+H_invDmH(1,1)*H_invDmH(3,2)*X(3,3&
&)+H_invDmH(1,2)*H_invDmH(3,1)*X(3,3)+H_invDmH(1,1)*H_invDmH(2,2)*&
&X(3,2)+H_invDmH(1,2)*H_invDmH(2,1)*X(3,2)+2*H_invDmH(1,1)*H_invDm&
&H(1,2)*X(3,1)
tt71 = H_invDmH(1,1)*X(3,8)*H_invDmH(8,3)+H_invDmH(1,3)*X(3,8)*H_&
&invDmH(8,1)+H_invDmH(1,1)*X(3,7)*H_invDmH(7,3)+H_invDmH(1,3)*X(3,&
&7)*H_invDmH(7,1)+H_invDmH(1,1)*X(3,6)*H_invDmH(6,3)+H_invDmH(1,3)&
&*X(3,6)*H_invDmH(6,1)+H_invDmH(1,1)*X(3,5)*H_invDmH(5,3)+H_invDmH&
&(1,3)*X(3,5)*H_invDmH(5,1)+H_invDmH(1,1)*X(3,4)*H_invDmH(4,3)+H_i&
&nvDmH(1,3)*X(3,4)*H_invDmH(4,1)+H_invDmH(1,1)*H_invDmH(3,3)*X(3,3&
&)+H_invDmH(1,3)*H_invDmH(3,1)*X(3,3)+H_invDmH(1,1)*H_invDmH(2,3)*&
&X(3,2)+H_invDmH(1,3)*H_invDmH(2,1)*X(3,2)+2*H_invDmH(1,1)*H_invDm&
&H(1,3)*X(3,1)
tt72 = H_invDmH(1,2)*X(3,8)*H_invDmH(8,3)+H_invDmH(1,3)*X(3,8)*H_&
&invDmH(8,2)+H_invDmH(1,2)*X(3,7)*H_invDmH(7,3)+H_invDmH(1,3)*X(3,&
&7)*H_invDmH(7,2)+H_invDmH(1,2)*X(3,6)*H_invDmH(6,3)+H_invDmH(1,3)&
&*X(3,6)*H_invDmH(6,2)+H_invDmH(1,2)*X(3,5)*H_invDmH(5,3)+H_invDmH&
&(1,3)*X(3,5)*H_invDmH(5,2)+H_invDmH(1,2)*X(3,4)*H_invDmH(4,3)+H_i&
&nvDmH(1,3)*X(3,4)*H_invDmH(4,2)+H_invDmH(1,2)*H_invDmH(3,3)*X(3,3&
&)+H_invDmH(1,3)*H_invDmH(3,2)*X(3,3)+H_invDmH(1,2)*H_invDmH(2,3)*&
&X(3,2)+H_invDmH(1,3)*H_invDmH(2,2)*X(3,2)+2*H_invDmH(1,2)*H_invDm&
&H(1,3)*X(3,1)
tt73 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt61*tt72+tt60*tt71+2*tt3*tt&
&18*tt22+tt54*tt70+2*tt2*tt11*tt15+2*tt1*tt4*tt8)+lam(1,1)*tt63*tt&
&69)
tt74 = H_invDmH(2,3)*tt18+H_invDmH(2,2)*tt11+H_invDmH(2,1)*tt4
tt75 = lam(1,1)*(H_invDmH(1,3)*H_invDmH(2,3)+H_invDmH(1,2)*H_invD&
&mH(2,2)+H_invDmH(1,1)*H_invDmH(2,1))*tt25
tt76 = H_invDmH(1,1)*H_invDmH(2,1)*tt10
tt77 = X(1,8)*H_invDmH(2,1)*H_invDmH(8,2)+X(1,8)*H_invDmH(2,2)*H_&
&invDmH(8,1)+X(1,7)*H_invDmH(2,1)*H_invDmH(7,2)+X(1,7)*H_invDmH(2,&
&2)*H_invDmH(7,1)+X(1,6)*H_invDmH(2,1)*H_invDmH(6,2)+X(1,6)*H_invD&
&mH(2,2)*H_invDmH(6,1)+X(1,5)*H_invDmH(2,1)*H_invDmH(5,2)+X(1,5)*H&
&_invDmH(2,2)*H_invDmH(5,1)+X(1,4)*H_invDmH(2,1)*H_invDmH(4,2)+X(1&
&,4)*H_invDmH(2,2)*H_invDmH(4,1)+X(1,3)*H_invDmH(2,1)*H_invDmH(3,2&
&)+X(1,3)*H_invDmH(2,2)*H_invDmH(3,1)+2*X(1,2)*H_invDmH(2,1)*H_inv&
&DmH(2,2)+H_invDmH(1,1)*X(1,1)*H_invDmH(2,2)+X(1,1)*H_invDmH(1,2)*&
&H_invDmH(2,1)
tt78 = (H_invDmH(1,1)*H_invDmH(2,2)+H_invDmH(1,2)*H_invDmH(2,1))*&
&tt52
tt79 = H_invDmH(1,2)*H_invDmH(2,2)*tt17
tt80 = X(1,8)*H_invDmH(2,1)*H_invDmH(8,3)+X(1,8)*H_invDmH(2,3)*H_&
&invDmH(8,1)+X(1,7)*H_invDmH(2,1)*H_invDmH(7,3)+X(1,7)*H_invDmH(2,&
&3)*H_invDmH(7,1)+X(1,6)*H_invDmH(2,1)*H_invDmH(6,3)+X(1,6)*H_invD&
&mH(2,3)*H_invDmH(6,1)+X(1,5)*H_invDmH(2,1)*H_invDmH(5,3)+X(1,5)*H&
&_invDmH(2,3)*H_invDmH(5,1)+X(1,4)*H_invDmH(2,1)*H_invDmH(4,3)+X(1&
&,4)*H_invDmH(2,3)*H_invDmH(4,1)+X(1,3)*H_invDmH(2,1)*H_invDmH(3,3&
&)+X(1,3)*H_invDmH(2,3)*H_invDmH(3,1)+2*X(1,2)*H_invDmH(2,1)*H_inv&
&DmH(2,3)+H_invDmH(1,1)*X(1,1)*H_invDmH(2,3)+X(1,1)*H_invDmH(1,3)*&
&H_invDmH(2,1)
tt81 = X(1,8)*H_invDmH(2,2)*H_invDmH(8,3)+X(1,8)*H_invDmH(2,3)*H_&
&invDmH(8,2)+X(1,7)*H_invDmH(2,2)*H_invDmH(7,3)+X(1,7)*H_invDmH(2,&
&3)*H_invDmH(7,2)+X(1,6)*H_invDmH(2,2)*H_invDmH(6,3)+X(1,6)*H_invD&
&mH(2,3)*H_invDmH(6,2)+X(1,5)*H_invDmH(2,2)*H_invDmH(5,3)+X(1,5)*H&
&_invDmH(2,3)*H_invDmH(5,2)+X(1,4)*H_invDmH(2,2)*H_invDmH(4,3)+X(1&
&,4)*H_invDmH(2,3)*H_invDmH(4,2)+X(1,3)*H_invDmH(2,2)*H_invDmH(3,3&
&)+X(1,3)*H_invDmH(2,3)*H_invDmH(3,2)+2*X(1,2)*H_invDmH(2,2)*H_inv&
&DmH(2,3)+X(1,1)*H_invDmH(1,2)*H_invDmH(2,3)+X(1,1)*H_invDmH(1,3)*&
&H_invDmH(2,2)
tt82 = (H_invDmH(1,1)*H_invDmH(2,3)+H_invDmH(1,3)*H_invDmH(2,1))*&
&tt56
tt83 = (H_invDmH(1,2)*H_invDmH(2,3)+H_invDmH(1,3)*H_invDmH(2,2))*&
&tt58
tt84 = H_invDmH(1,3)*H_invDmH(2,3)*tt24
tt85 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt84+2*H_invDmH(1,3)*H_invDm&
&H(2,3)*tt19+tt83+tt82+tt61*tt81+tt60*tt80+tt79+2*H_invDmH(1,2)*H_&
&invDmH(2,2)*tt12+tt78+tt54*tt77+tt76+2*H_invDmH(1,1)*H_invDmH(2,1&
&)*tt5)+tt75+lam(1,1)*tt63*tt74)
tt86 = H_invDmH(2,3)*tt20+H_invDmH(2,2)*tt13+H_invDmH(2,1)*tt6
tt87 = 2*H_invDmH(1,1)*H_invDmH(2,1)*tt4*tt6
tt88 = 2*H_invDmH(1,2)*H_invDmH(2,2)*tt11*tt13
tt89 = H_invDmH(2,1)*X(2,8)*H_invDmH(8,2)+H_invDmH(2,2)*X(2,8)*H_&
&invDmH(8,1)+H_invDmH(2,1)*X(2,7)*H_invDmH(7,2)+H_invDmH(2,2)*X(2,&
&7)*H_invDmH(7,1)+H_invDmH(2,1)*X(2,6)*H_invDmH(6,2)+H_invDmH(2,2)&
&*X(2,6)*H_invDmH(6,1)+H_invDmH(2,1)*X(2,5)*H_invDmH(5,2)+H_invDmH&
&(2,2)*X(2,5)*H_invDmH(5,1)+H_invDmH(2,1)*X(2,4)*H_invDmH(4,2)+H_i&
&nvDmH(2,2)*X(2,4)*H_invDmH(4,1)+H_invDmH(2,1)*X(2,3)*H_invDmH(3,2&
&)+H_invDmH(2,2)*X(2,3)*H_invDmH(3,1)+2*H_invDmH(2,1)*H_invDmH(2,2&
&)*X(2,2)+H_invDmH(1,1)*X(2,1)*H_invDmH(2,2)+H_invDmH(1,2)*H_invDm&
&H(2,1)*X(2,1)
tt90 = 2*H_invDmH(1,3)*H_invDmH(2,3)*tt18*tt20
tt91 = H_invDmH(2,1)*X(2,8)*H_invDmH(8,3)+H_invDmH(2,3)*X(2,8)*H_&
&invDmH(8,1)+H_invDmH(2,1)*X(2,7)*H_invDmH(7,3)+H_invDmH(2,3)*X(2,&
&7)*H_invDmH(7,1)+H_invDmH(2,1)*X(2,6)*H_invDmH(6,3)+H_invDmH(2,3)&
&*X(2,6)*H_invDmH(6,1)+H_invDmH(2,1)*X(2,5)*H_invDmH(5,3)+H_invDmH&
&(2,3)*X(2,5)*H_invDmH(5,1)+H_invDmH(2,1)*X(2,4)*H_invDmH(4,3)+H_i&
&nvDmH(2,3)*X(2,4)*H_invDmH(4,1)+H_invDmH(2,1)*X(2,3)*H_invDmH(3,3&
&)+H_invDmH(2,3)*X(2,3)*H_invDmH(3,1)+2*H_invDmH(2,1)*X(2,2)*H_inv&
&DmH(2,3)+H_invDmH(1,1)*X(2,1)*H_invDmH(2,3)+H_invDmH(1,3)*H_invDm&
&H(2,1)*X(2,1)
tt92 = H_invDmH(2,2)*X(2,8)*H_invDmH(8,3)+H_invDmH(2,3)*X(2,8)*H_&
&invDmH(8,2)+H_invDmH(2,2)*X(2,7)*H_invDmH(7,3)+H_invDmH(2,3)*X(2,&
&7)*H_invDmH(7,2)+H_invDmH(2,2)*X(2,6)*H_invDmH(6,3)+H_invDmH(2,3)&
&*X(2,6)*H_invDmH(6,2)+H_invDmH(2,2)*X(2,5)*H_invDmH(5,3)+H_invDmH&
&(2,3)*X(2,5)*H_invDmH(5,2)+H_invDmH(2,2)*X(2,4)*H_invDmH(4,3)+H_i&
&nvDmH(2,3)*X(2,4)*H_invDmH(4,2)+H_invDmH(2,2)*X(2,3)*H_invDmH(3,3&
&)+H_invDmH(2,3)*X(2,3)*H_invDmH(3,2)+2*H_invDmH(2,2)*X(2,2)*H_inv&
&DmH(2,3)+H_invDmH(1,2)*X(2,1)*H_invDmH(2,3)+H_invDmH(1,3)*X(2,1)*&
&H_invDmH(2,2)
tt93 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt61*tt92+tt60*tt91+tt90+tt5&
&4*tt89+tt88+tt87)+lam(1,1)*tt63*tt86)
tt94 = H_invDmH(2,3)*tt22+H_invDmH(2,2)*tt15+H_invDmH(2,1)*tt8
tt95 = 2*H_invDmH(1,1)*H_invDmH(2,1)*tt4*tt8
tt96 = 2*H_invDmH(1,2)*H_invDmH(2,2)*tt11*tt15
tt97 = H_invDmH(2,1)*X(3,8)*H_invDmH(8,2)+H_invDmH(2,2)*X(3,8)*H_&
&invDmH(8,1)+H_invDmH(2,1)*X(3,7)*H_invDmH(7,2)+H_invDmH(2,2)*X(3,&
&7)*H_invDmH(7,1)+H_invDmH(2,1)*X(3,6)*H_invDmH(6,2)+H_invDmH(2,2)&
&*X(3,6)*H_invDmH(6,1)+H_invDmH(2,1)*X(3,5)*H_invDmH(5,2)+H_invDmH&
&(2,2)*X(3,5)*H_invDmH(5,1)+H_invDmH(2,1)*X(3,4)*H_invDmH(4,2)+H_i&
&nvDmH(2,2)*X(3,4)*H_invDmH(4,1)+H_invDmH(2,1)*H_invDmH(3,2)*X(3,3&
&)+H_invDmH(2,2)*H_invDmH(3,1)*X(3,3)+2*H_invDmH(2,1)*H_invDmH(2,2&
&)*X(3,2)+H_invDmH(1,1)*H_invDmH(2,2)*X(3,1)+H_invDmH(1,2)*H_invDm&
&H(2,1)*X(3,1)
tt98 = 2*H_invDmH(1,3)*H_invDmH(2,3)*tt18*tt22
tt99 = H_invDmH(2,1)*X(3,8)*H_invDmH(8,3)+H_invDmH(2,3)*X(3,8)*H_&
&invDmH(8,1)+H_invDmH(2,1)*X(3,7)*H_invDmH(7,3)+H_invDmH(2,3)*X(3,&
&7)*H_invDmH(7,1)+H_invDmH(2,1)*X(3,6)*H_invDmH(6,3)+H_invDmH(2,3)&
&*X(3,6)*H_invDmH(6,1)+H_invDmH(2,1)*X(3,5)*H_invDmH(5,3)+H_invDmH&
&(2,3)*X(3,5)*H_invDmH(5,1)+H_invDmH(2,1)*X(3,4)*H_invDmH(4,3)+H_i&
&nvDmH(2,3)*X(3,4)*H_invDmH(4,1)+H_invDmH(2,1)*H_invDmH(3,3)*X(3,3&
&)+H_invDmH(2,3)*H_invDmH(3,1)*X(3,3)+2*H_invDmH(2,1)*H_invDmH(2,3&
&)*X(3,2)+H_invDmH(1,1)*H_invDmH(2,3)*X(3,1)+H_invDmH(1,3)*H_invDm&
&H(2,1)*X(3,1)
tt100 = H_invDmH(2,2)*X(3,8)*H_invDmH(8,3)+H_invDmH(2,3)*X(3,8)*H&
&_invDmH(8,2)+H_invDmH(2,2)*X(3,7)*H_invDmH(7,3)+H_invDmH(2,3)*X(3&
&,7)*H_invDmH(7,2)+H_invDmH(2,2)*X(3,6)*H_invDmH(6,3)+H_invDmH(2,3&
&)*X(3,6)*H_invDmH(6,2)+H_invDmH(2,2)*X(3,5)*H_invDmH(5,3)+H_invDm&
&H(2,3)*X(3,5)*H_invDmH(5,2)+H_invDmH(2,2)*X(3,4)*H_invDmH(4,3)+H_&
&invDmH(2,3)*X(3,4)*H_invDmH(4,2)+H_invDmH(2,2)*H_invDmH(3,3)*X(3,&
&3)+H_invDmH(2,3)*H_invDmH(3,2)*X(3,3)+2*H_invDmH(2,2)*H_invDmH(2,&
&3)*X(3,2)+H_invDmH(1,2)*H_invDmH(2,3)*X(3,1)+H_invDmH(1,3)*H_invD&
&mH(2,2)*X(3,1)
tt101 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt61*tt100+tt60*tt99+tt98+t&
&t54*tt97+tt96+tt95)+lam(1,1)*tt63*tt94)
tt102 = H_invDmH(3,3)*tt18+H_invDmH(3,2)*tt11+H_invDmH(3,1)*tt4
tt103 = lam(1,1)*(H_invDmH(1,3)*H_invDmH(3,3)+H_invDmH(1,2)*H_inv&
&DmH(3,2)+H_invDmH(1,1)*H_invDmH(3,1))*tt25
tt104 = H_invDmH(1,1)*H_invDmH(3,1)*tt10
tt105 = X(1,8)*H_invDmH(3,1)*H_invDmH(8,2)+X(1,8)*H_invDmH(3,2)*H&
&_invDmH(8,1)+X(1,7)*H_invDmH(3,1)*H_invDmH(7,2)+X(1,7)*H_invDmH(3&
&,2)*H_invDmH(7,1)+X(1,6)*H_invDmH(3,1)*H_invDmH(6,2)+X(1,6)*H_inv&
&DmH(3,2)*H_invDmH(6,1)+X(1,5)*H_invDmH(3,1)*H_invDmH(5,2)+X(1,5)*&
&H_invDmH(3,2)*H_invDmH(5,1)+X(1,4)*H_invDmH(3,1)*H_invDmH(4,2)+X(&
&1,4)*H_invDmH(3,2)*H_invDmH(4,1)+2*X(1,3)*H_invDmH(3,1)*H_invDmH(&
&3,2)+X(1,2)*H_invDmH(2,1)*H_invDmH(3,2)+H_invDmH(1,1)*X(1,1)*H_in&
&vDmH(3,2)+X(1,2)*H_invDmH(2,2)*H_invDmH(3,1)+X(1,1)*H_invDmH(1,2)&
&*H_invDmH(3,1)
tt106 = (H_invDmH(1,1)*H_invDmH(3,2)+H_invDmH(1,2)*H_invDmH(3,1))&
&*tt52
tt107 = H_invDmH(1,2)*H_invDmH(3,2)*tt17
tt108 = X(1,8)*H_invDmH(3,1)*H_invDmH(8,3)+X(1,8)*H_invDmH(3,3)*H&
&_invDmH(8,1)+X(1,7)*H_invDmH(3,1)*H_invDmH(7,3)+X(1,7)*H_invDmH(3&
&,3)*H_invDmH(7,1)+X(1,6)*H_invDmH(3,1)*H_invDmH(6,3)+X(1,6)*H_inv&
&DmH(3,3)*H_invDmH(6,1)+X(1,5)*H_invDmH(3,1)*H_invDmH(5,3)+X(1,5)*&
&H_invDmH(3,3)*H_invDmH(5,1)+X(1,4)*H_invDmH(3,1)*H_invDmH(4,3)+X(&
&1,4)*H_invDmH(3,3)*H_invDmH(4,1)+2*X(1,3)*H_invDmH(3,1)*H_invDmH(&
&3,3)+X(1,2)*H_invDmH(2,1)*H_invDmH(3,3)+H_invDmH(1,1)*X(1,1)*H_in&
&vDmH(3,3)+X(1,2)*H_invDmH(2,3)*H_invDmH(3,1)+X(1,1)*H_invDmH(1,3)&
&*H_invDmH(3,1)
tt109 = X(1,8)*H_invDmH(3,2)*H_invDmH(8,3)+X(1,8)*H_invDmH(3,3)*H&
&_invDmH(8,2)+X(1,7)*H_invDmH(3,2)*H_invDmH(7,3)+X(1,7)*H_invDmH(3&
&,3)*H_invDmH(7,2)+X(1,6)*H_invDmH(3,2)*H_invDmH(6,3)+X(1,6)*H_inv&
&DmH(3,3)*H_invDmH(6,2)+X(1,5)*H_invDmH(3,2)*H_invDmH(5,3)+X(1,5)*&
&H_invDmH(3,3)*H_invDmH(5,2)+X(1,4)*H_invDmH(3,2)*H_invDmH(4,3)+X(&
&1,4)*H_invDmH(3,3)*H_invDmH(4,2)+2*X(1,3)*H_invDmH(3,2)*H_invDmH(&
&3,3)+X(1,2)*H_invDmH(2,2)*H_invDmH(3,3)+X(1,1)*H_invDmH(1,2)*H_in&
&vDmH(3,3)+X(1,2)*H_invDmH(2,3)*H_invDmH(3,2)+X(1,1)*H_invDmH(1,3)&
&*H_invDmH(3,2)
tt110 = (H_invDmH(1,1)*H_invDmH(3,3)+H_invDmH(1,3)*H_invDmH(3,1))&
&*tt56
tt111 = (H_invDmH(1,2)*H_invDmH(3,3)+H_invDmH(1,3)*H_invDmH(3,2))&
&*tt58
tt112 = H_invDmH(1,3)*H_invDmH(3,3)*tt24
tt113 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt112+2*H_invDmH(1,3)*H_inv&
&DmH(3,3)*tt19+tt111+tt110+tt61*tt109+tt60*tt108+tt107+2*H_invDmH(&
&1,2)*H_invDmH(3,2)*tt12+tt106+tt54*tt105+tt104+2*H_invDmH(1,1)*H_&
&invDmH(3,1)*tt5)+tt103+lam(1,1)*tt63*tt102)
tt114 = H_invDmH(3,3)*tt20+H_invDmH(3,2)*tt13+H_invDmH(3,1)*tt6
tt115 = 2*H_invDmH(1,1)*H_invDmH(3,1)*tt4*tt6
tt116 = 2*H_invDmH(1,2)*H_invDmH(3,2)*tt11*tt13
tt117 = X(2,8)*H_invDmH(3,1)*H_invDmH(8,2)+X(2,8)*H_invDmH(3,2)*H&
&_invDmH(8,1)+X(2,7)*H_invDmH(3,1)*H_invDmH(7,2)+X(2,7)*H_invDmH(3&
&,2)*H_invDmH(7,1)+X(2,6)*H_invDmH(3,1)*H_invDmH(6,2)+X(2,6)*H_inv&
&DmH(3,2)*H_invDmH(6,1)+X(2,5)*H_invDmH(3,1)*H_invDmH(5,2)+X(2,5)*&
&H_invDmH(3,2)*H_invDmH(5,1)+X(2,4)*H_invDmH(3,1)*H_invDmH(4,2)+X(&
&2,4)*H_invDmH(3,2)*H_invDmH(4,1)+2*X(2,3)*H_invDmH(3,1)*H_invDmH(&
&3,2)+H_invDmH(2,1)*X(2,2)*H_invDmH(3,2)+H_invDmH(1,1)*X(2,1)*H_in&
&vDmH(3,2)+H_invDmH(2,2)*X(2,2)*H_invDmH(3,1)+H_invDmH(1,2)*X(2,1)&
&*H_invDmH(3,1)
tt118 = 2*H_invDmH(1,3)*H_invDmH(3,3)*tt18*tt20
tt119 = X(2,8)*H_invDmH(3,1)*H_invDmH(8,3)+X(2,8)*H_invDmH(3,3)*H&
&_invDmH(8,1)+X(2,7)*H_invDmH(3,1)*H_invDmH(7,3)+X(2,7)*H_invDmH(3&
&,3)*H_invDmH(7,1)+X(2,6)*H_invDmH(3,1)*H_invDmH(6,3)+X(2,6)*H_inv&
&DmH(3,3)*H_invDmH(6,1)+X(2,5)*H_invDmH(3,1)*H_invDmH(5,3)+X(2,5)*&
&H_invDmH(3,3)*H_invDmH(5,1)+X(2,4)*H_invDmH(3,1)*H_invDmH(4,3)+X(&
&2,4)*H_invDmH(3,3)*H_invDmH(4,1)+2*X(2,3)*H_invDmH(3,1)*H_invDmH(&
&3,3)+H_invDmH(2,1)*X(2,2)*H_invDmH(3,3)+H_invDmH(1,1)*X(2,1)*H_in&
&vDmH(3,3)+X(2,2)*H_invDmH(2,3)*H_invDmH(3,1)+H_invDmH(1,3)*X(2,1)&
&*H_invDmH(3,1)
tt120 = X(2,8)*H_invDmH(3,2)*H_invDmH(8,3)+X(2,8)*H_invDmH(3,3)*H&
&_invDmH(8,2)+X(2,7)*H_invDmH(3,2)*H_invDmH(7,3)+X(2,7)*H_invDmH(3&
&,3)*H_invDmH(7,2)+X(2,6)*H_invDmH(3,2)*H_invDmH(6,3)+X(2,6)*H_inv&
&DmH(3,3)*H_invDmH(6,2)+X(2,5)*H_invDmH(3,2)*H_invDmH(5,3)+X(2,5)*&
&H_invDmH(3,3)*H_invDmH(5,2)+X(2,4)*H_invDmH(3,2)*H_invDmH(4,3)+X(&
&2,4)*H_invDmH(3,3)*H_invDmH(4,2)+2*X(2,3)*H_invDmH(3,2)*H_invDmH(&
&3,3)+H_invDmH(2,2)*X(2,2)*H_invDmH(3,3)+H_invDmH(1,2)*X(2,1)*H_in&
&vDmH(3,3)+X(2,2)*H_invDmH(2,3)*H_invDmH(3,2)+H_invDmH(1,3)*X(2,1)&
&*H_invDmH(3,2)
tt121 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt61*tt120+tt60*tt119+tt118&
&+tt54*tt117+tt116+tt115)+lam(1,1)*tt63*tt114)
tt122 = H_invDmH(3,3)*tt22+H_invDmH(3,2)*tt15+H_invDmH(3,1)*tt8
tt123 = 2*H_invDmH(1,1)*H_invDmH(3,1)*tt4*tt8
tt124 = 2*H_invDmH(1,2)*H_invDmH(3,2)*tt11*tt15
tt125 = H_invDmH(3,1)*X(3,8)*H_invDmH(8,2)+H_invDmH(3,2)*X(3,8)*H&
&_invDmH(8,1)+H_invDmH(3,1)*X(3,7)*H_invDmH(7,2)+H_invDmH(3,2)*X(3&
&,7)*H_invDmH(7,1)+H_invDmH(3,1)*X(3,6)*H_invDmH(6,2)+H_invDmH(3,2&
&)*X(3,6)*H_invDmH(6,1)+H_invDmH(3,1)*X(3,5)*H_invDmH(5,2)+H_invDm&
&H(3,2)*X(3,5)*H_invDmH(5,1)+H_invDmH(3,1)*X(3,4)*H_invDmH(4,2)+H_&
&invDmH(3,2)*X(3,4)*H_invDmH(4,1)+2*H_invDmH(3,1)*H_invDmH(3,2)*X(&
&3,3)+H_invDmH(2,1)*H_invDmH(3,2)*X(3,2)+H_invDmH(2,2)*H_invDmH(3,&
&1)*X(3,2)+H_invDmH(1,1)*X(3,1)*H_invDmH(3,2)+H_invDmH(1,2)*H_invD&
&mH(3,1)*X(3,1)
tt126 = 2*H_invDmH(1,3)*H_invDmH(3,3)*tt18*tt22
tt127 = H_invDmH(3,1)*X(3,8)*H_invDmH(8,3)+H_invDmH(3,3)*X(3,8)*H&
&_invDmH(8,1)+H_invDmH(3,1)*X(3,7)*H_invDmH(7,3)+H_invDmH(3,3)*X(3&
&,7)*H_invDmH(7,1)+H_invDmH(3,1)*X(3,6)*H_invDmH(6,3)+H_invDmH(3,3&
&)*X(3,6)*H_invDmH(6,1)+H_invDmH(3,1)*X(3,5)*H_invDmH(5,3)+H_invDm&
&H(3,3)*X(3,5)*H_invDmH(5,1)+H_invDmH(3,1)*X(3,4)*H_invDmH(4,3)+H_&
&invDmH(3,3)*X(3,4)*H_invDmH(4,1)+2*H_invDmH(3,1)*H_invDmH(3,3)*X(&
&3,3)+H_invDmH(2,1)*X(3,2)*H_invDmH(3,3)+H_invDmH(1,1)*X(3,1)*H_in&
&vDmH(3,3)+H_invDmH(2,3)*H_invDmH(3,1)*X(3,2)+H_invDmH(1,3)*H_invD&
&mH(3,1)*X(3,1)
tt128 = H_invDmH(3,2)*X(3,8)*H_invDmH(8,3)+H_invDmH(3,3)*X(3,8)*H&
&_invDmH(8,2)+H_invDmH(3,2)*X(3,7)*H_invDmH(7,3)+H_invDmH(3,3)*X(3&
&,7)*H_invDmH(7,2)+H_invDmH(3,2)*X(3,6)*H_invDmH(6,3)+H_invDmH(3,3&
&)*X(3,6)*H_invDmH(6,2)+H_invDmH(3,2)*X(3,5)*H_invDmH(5,3)+H_invDm&
&H(3,3)*X(3,5)*H_invDmH(5,2)+H_invDmH(3,2)*X(3,4)*H_invDmH(4,3)+H_&
&invDmH(3,3)*X(3,4)*H_invDmH(4,2)+2*H_invDmH(3,2)*H_invDmH(3,3)*X(&
&3,3)+H_invDmH(2,2)*X(3,2)*H_invDmH(3,3)+H_invDmH(1,2)*X(3,1)*H_in&
&vDmH(3,3)+H_invDmH(2,3)*H_invDmH(3,2)*X(3,2)+H_invDmH(1,3)*X(3,1)&
&*H_invDmH(3,2)
tt129 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt61*tt128+tt60*tt127+tt126&
&+tt54*tt125+tt124+tt123)+lam(1,1)*tt63*tt122)
tt130 = H_invDmH(4,3)*tt18+H_invDmH(4,2)*tt11+H_invDmH(4,1)*tt4
tt131 = lam(1,1)*(H_invDmH(1,3)*H_invDmH(4,3)+H_invDmH(1,2)*H_inv&
&DmH(4,2)+H_invDmH(1,1)*H_invDmH(4,1))*tt25
tt132 = H_invDmH(1,1)*H_invDmH(4,1)*tt10
tt133 = X(1,8)*H_invDmH(4,1)*H_invDmH(8,2)+X(1,8)*H_invDmH(4,2)*H&
&_invDmH(8,1)+X(1,7)*H_invDmH(4,1)*H_invDmH(7,2)+X(1,7)*H_invDmH(4&
&,2)*H_invDmH(7,1)+X(1,6)*H_invDmH(4,1)*H_invDmH(6,2)+X(1,6)*H_inv&
&DmH(4,2)*H_invDmH(6,1)+X(1,5)*H_invDmH(4,1)*H_invDmH(5,2)+X(1,5)*&
&H_invDmH(4,2)*H_invDmH(5,1)+2*X(1,4)*H_invDmH(4,1)*H_invDmH(4,2)+&
&X(1,3)*H_invDmH(3,1)*H_invDmH(4,2)+X(1,2)*H_invDmH(2,1)*H_invDmH(&
&4,2)+H_invDmH(1,1)*X(1,1)*H_invDmH(4,2)+X(1,3)*H_invDmH(3,2)*H_in&
&vDmH(4,1)+X(1,2)*H_invDmH(2,2)*H_invDmH(4,1)+X(1,1)*H_invDmH(1,2)&
&*H_invDmH(4,1)
tt134 = (H_invDmH(1,1)*H_invDmH(4,2)+H_invDmH(1,2)*H_invDmH(4,1))&
&*tt52
tt135 = H_invDmH(1,2)*H_invDmH(4,2)*tt17
tt136 = X(1,8)*H_invDmH(4,1)*H_invDmH(8,3)+X(1,8)*H_invDmH(4,3)*H&
&_invDmH(8,1)+X(1,7)*H_invDmH(4,1)*H_invDmH(7,3)+X(1,7)*H_invDmH(4&
&,3)*H_invDmH(7,1)+X(1,6)*H_invDmH(4,1)*H_invDmH(6,3)+X(1,6)*H_inv&
&DmH(4,3)*H_invDmH(6,1)+X(1,5)*H_invDmH(4,1)*H_invDmH(5,3)+X(1,5)*&
&H_invDmH(4,3)*H_invDmH(5,1)+2*X(1,4)*H_invDmH(4,1)*H_invDmH(4,3)+&
&X(1,3)*H_invDmH(3,1)*H_invDmH(4,3)+X(1,2)*H_invDmH(2,1)*H_invDmH(&
&4,3)+H_invDmH(1,1)*X(1,1)*H_invDmH(4,3)+X(1,3)*H_invDmH(3,3)*H_in&
&vDmH(4,1)+X(1,2)*H_invDmH(2,3)*H_invDmH(4,1)+X(1,1)*H_invDmH(1,3)&
&*H_invDmH(4,1)
tt137 = X(1,8)*H_invDmH(4,2)*H_invDmH(8,3)+X(1,8)*H_invDmH(4,3)*H&
&_invDmH(8,2)+X(1,7)*H_invDmH(4,2)*H_invDmH(7,3)+X(1,7)*H_invDmH(4&
&,3)*H_invDmH(7,2)+X(1,6)*H_invDmH(4,2)*H_invDmH(6,3)+X(1,6)*H_inv&
&DmH(4,3)*H_invDmH(6,2)+X(1,5)*H_invDmH(4,2)*H_invDmH(5,3)+X(1,5)*&
&H_invDmH(4,3)*H_invDmH(5,2)+2*X(1,4)*H_invDmH(4,2)*H_invDmH(4,3)+&
&X(1,3)*H_invDmH(3,2)*H_invDmH(4,3)+X(1,2)*H_invDmH(2,2)*H_invDmH(&
&4,3)+X(1,1)*H_invDmH(1,2)*H_invDmH(4,3)+X(1,3)*H_invDmH(3,3)*H_in&
&vDmH(4,2)+X(1,2)*H_invDmH(2,3)*H_invDmH(4,2)+X(1,1)*H_invDmH(1,3)&
&*H_invDmH(4,2)
tt138 = (H_invDmH(1,1)*H_invDmH(4,3)+H_invDmH(1,3)*H_invDmH(4,1))&
&*tt56
tt139 = (H_invDmH(1,2)*H_invDmH(4,3)+H_invDmH(1,3)*H_invDmH(4,2))&
&*tt58
tt140 = H_invDmH(1,3)*H_invDmH(4,3)*tt24
tt141 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt140+2*H_invDmH(1,3)*H_inv&
&DmH(4,3)*tt19+tt139+tt138+tt61*tt137+tt60*tt136+tt135+2*H_invDmH(&
&1,2)*H_invDmH(4,2)*tt12+tt134+tt54*tt133+tt132+2*H_invDmH(1,1)*H_&
&invDmH(4,1)*tt5)+tt131+lam(1,1)*tt63*tt130)
tt142 = H_invDmH(4,3)*tt20+H_invDmH(4,2)*tt13+H_invDmH(4,1)*tt6
tt143 = 2*H_invDmH(1,1)*H_invDmH(4,1)*tt4*tt6
tt144 = 2*H_invDmH(1,2)*H_invDmH(4,2)*tt11*tt13
tt145 = X(2,8)*H_invDmH(4,1)*H_invDmH(8,2)+X(2,8)*H_invDmH(4,2)*H&
&_invDmH(8,1)+X(2,7)*H_invDmH(4,1)*H_invDmH(7,2)+X(2,7)*H_invDmH(4&
&,2)*H_invDmH(7,1)+X(2,6)*H_invDmH(4,1)*H_invDmH(6,2)+X(2,6)*H_inv&
&DmH(4,2)*H_invDmH(6,1)+X(2,5)*H_invDmH(4,1)*H_invDmH(5,2)+X(2,5)*&
&H_invDmH(4,2)*H_invDmH(5,1)+2*X(2,4)*H_invDmH(4,1)*H_invDmH(4,2)+&
&X(2,3)*H_invDmH(3,1)*H_invDmH(4,2)+H_invDmH(2,1)*X(2,2)*H_invDmH(&
&4,2)+H_invDmH(1,1)*X(2,1)*H_invDmH(4,2)+X(2,3)*H_invDmH(3,2)*H_in&
&vDmH(4,1)+H_invDmH(2,2)*X(2,2)*H_invDmH(4,1)+H_invDmH(1,2)*X(2,1)&
&*H_invDmH(4,1)
tt146 = 2*H_invDmH(1,3)*H_invDmH(4,3)*tt18*tt20
tt147 = X(2,8)*H_invDmH(4,1)*H_invDmH(8,3)+X(2,8)*H_invDmH(4,3)*H&
&_invDmH(8,1)+X(2,7)*H_invDmH(4,1)*H_invDmH(7,3)+X(2,7)*H_invDmH(4&
&,3)*H_invDmH(7,1)+X(2,6)*H_invDmH(4,1)*H_invDmH(6,3)+X(2,6)*H_inv&
&DmH(4,3)*H_invDmH(6,1)+X(2,5)*H_invDmH(4,1)*H_invDmH(5,3)+X(2,5)*&
&H_invDmH(4,3)*H_invDmH(5,1)+2*X(2,4)*H_invDmH(4,1)*H_invDmH(4,3)+&
&X(2,3)*H_invDmH(3,1)*H_invDmH(4,3)+H_invDmH(2,1)*X(2,2)*H_invDmH(&
&4,3)+H_invDmH(1,1)*X(2,1)*H_invDmH(4,3)+X(2,3)*H_invDmH(3,3)*H_in&
&vDmH(4,1)+X(2,2)*H_invDmH(2,3)*H_invDmH(4,1)+H_invDmH(1,3)*X(2,1)&
&*H_invDmH(4,1)
tt148 = X(2,8)*H_invDmH(4,2)*H_invDmH(8,3)+X(2,8)*H_invDmH(4,3)*H&
&_invDmH(8,2)+X(2,7)*H_invDmH(4,2)*H_invDmH(7,3)+X(2,7)*H_invDmH(4&
&,3)*H_invDmH(7,2)+X(2,6)*H_invDmH(4,2)*H_invDmH(6,3)+X(2,6)*H_inv&
&DmH(4,3)*H_invDmH(6,2)+X(2,5)*H_invDmH(4,2)*H_invDmH(5,3)+X(2,5)*&
&H_invDmH(4,3)*H_invDmH(5,2)+2*X(2,4)*H_invDmH(4,2)*H_invDmH(4,3)+&
&X(2,3)*H_invDmH(3,2)*H_invDmH(4,3)+H_invDmH(2,2)*X(2,2)*H_invDmH(&
&4,3)+H_invDmH(1,2)*X(2,1)*H_invDmH(4,3)+X(2,3)*H_invDmH(3,3)*H_in&
&vDmH(4,2)+X(2,2)*H_invDmH(2,3)*H_invDmH(4,2)+H_invDmH(1,3)*X(2,1)&
&*H_invDmH(4,2)
tt149 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt61*tt148+tt60*tt147+tt146&
&+tt54*tt145+tt144+tt143)+lam(1,1)*tt63*tt142)
tt150 = H_invDmH(4,3)*tt22+H_invDmH(4,2)*tt15+H_invDmH(4,1)*tt8
tt151 = 2*H_invDmH(1,1)*H_invDmH(4,1)*tt4*tt8
tt152 = 2*H_invDmH(1,2)*H_invDmH(4,2)*tt11*tt15
tt153 = X(3,8)*H_invDmH(4,1)*H_invDmH(8,2)+X(3,8)*H_invDmH(4,2)*H&
&_invDmH(8,1)+X(3,7)*H_invDmH(4,1)*H_invDmH(7,2)+X(3,7)*H_invDmH(4&
&,2)*H_invDmH(7,1)+X(3,6)*H_invDmH(4,1)*H_invDmH(6,2)+X(3,6)*H_inv&
&DmH(4,2)*H_invDmH(6,1)+X(3,5)*H_invDmH(4,1)*H_invDmH(5,2)+X(3,5)*&
&H_invDmH(4,2)*H_invDmH(5,1)+2*X(3,4)*H_invDmH(4,1)*H_invDmH(4,2)+&
&H_invDmH(3,1)*X(3,3)*H_invDmH(4,2)+H_invDmH(2,1)*X(3,2)*H_invDmH(&
&4,2)+H_invDmH(1,1)*X(3,1)*H_invDmH(4,2)+H_invDmH(3,2)*X(3,3)*H_in&
&vDmH(4,1)+H_invDmH(2,2)*X(3,2)*H_invDmH(4,1)+H_invDmH(1,2)*X(3,1)&
&*H_invDmH(4,1)
tt154 = 2*H_invDmH(1,3)*H_invDmH(4,3)*tt18*tt22
tt155 = X(3,8)*H_invDmH(4,1)*H_invDmH(8,3)+X(3,8)*H_invDmH(4,3)*H&
&_invDmH(8,1)+X(3,7)*H_invDmH(4,1)*H_invDmH(7,3)+X(3,7)*H_invDmH(4&
&,3)*H_invDmH(7,1)+X(3,6)*H_invDmH(4,1)*H_invDmH(6,3)+X(3,6)*H_inv&
&DmH(4,3)*H_invDmH(6,1)+X(3,5)*H_invDmH(4,1)*H_invDmH(5,3)+X(3,5)*&
&H_invDmH(4,3)*H_invDmH(5,1)+2*X(3,4)*H_invDmH(4,1)*H_invDmH(4,3)+&
&H_invDmH(3,1)*X(3,3)*H_invDmH(4,3)+H_invDmH(2,1)*X(3,2)*H_invDmH(&
&4,3)+H_invDmH(1,1)*X(3,1)*H_invDmH(4,3)+H_invDmH(3,3)*X(3,3)*H_in&
&vDmH(4,1)+H_invDmH(2,3)*X(3,2)*H_invDmH(4,1)+H_invDmH(1,3)*X(3,1)&
&*H_invDmH(4,1)
tt156 = X(3,8)*H_invDmH(4,2)*H_invDmH(8,3)+X(3,8)*H_invDmH(4,3)*H&
&_invDmH(8,2)+X(3,7)*H_invDmH(4,2)*H_invDmH(7,3)+X(3,7)*H_invDmH(4&
&,3)*H_invDmH(7,2)+X(3,6)*H_invDmH(4,2)*H_invDmH(6,3)+X(3,6)*H_inv&
&DmH(4,3)*H_invDmH(6,2)+X(3,5)*H_invDmH(4,2)*H_invDmH(5,3)+X(3,5)*&
&H_invDmH(4,3)*H_invDmH(5,2)+2*X(3,4)*H_invDmH(4,2)*H_invDmH(4,3)+&
&H_invDmH(3,2)*X(3,3)*H_invDmH(4,3)+H_invDmH(2,2)*X(3,2)*H_invDmH(&
&4,3)+H_invDmH(1,2)*X(3,1)*H_invDmH(4,3)+H_invDmH(3,3)*X(3,3)*H_in&
&vDmH(4,2)+H_invDmH(2,3)*X(3,2)*H_invDmH(4,2)+H_invDmH(1,3)*X(3,1)&
&*H_invDmH(4,2)
tt157 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt61*tt156+tt60*tt155+tt154&
&+tt54*tt153+tt152+tt151)+lam(1,1)*tt63*tt150)
tt158 = H_invDmH(5,3)*tt18+H_invDmH(5,2)*tt11+H_invDmH(5,1)*tt4
tt159 = lam(1,1)*(H_invDmH(1,3)*H_invDmH(5,3)+H_invDmH(1,2)*H_inv&
&DmH(5,2)+H_invDmH(1,1)*H_invDmH(5,1))*tt25
tt160 = H_invDmH(1,1)*H_invDmH(5,1)*tt10
tt161 = X(1,8)*H_invDmH(5,1)*H_invDmH(8,2)+X(1,8)*H_invDmH(5,2)*H&
&_invDmH(8,1)+X(1,7)*H_invDmH(5,1)*H_invDmH(7,2)+X(1,7)*H_invDmH(5&
&,2)*H_invDmH(7,1)+X(1,6)*H_invDmH(5,1)*H_invDmH(6,2)+X(1,6)*H_inv&
&DmH(5,2)*H_invDmH(6,1)+2*X(1,5)*H_invDmH(5,1)*H_invDmH(5,2)+X(1,4&
&)*H_invDmH(4,1)*H_invDmH(5,2)+X(1,3)*H_invDmH(3,1)*H_invDmH(5,2)+&
&X(1,2)*H_invDmH(2,1)*H_invDmH(5,2)+H_invDmH(1,1)*X(1,1)*H_invDmH(&
&5,2)+X(1,4)*H_invDmH(4,2)*H_invDmH(5,1)+X(1,3)*H_invDmH(3,2)*H_in&
&vDmH(5,1)+X(1,2)*H_invDmH(2,2)*H_invDmH(5,1)+X(1,1)*H_invDmH(1,2)&
&*H_invDmH(5,1)
tt162 = (H_invDmH(1,1)*H_invDmH(5,2)+H_invDmH(1,2)*H_invDmH(5,1))&
&*tt52
tt163 = H_invDmH(1,2)*H_invDmH(5,2)*tt17
tt164 = X(1,8)*H_invDmH(5,1)*H_invDmH(8,3)+X(1,8)*H_invDmH(5,3)*H&
&_invDmH(8,1)+X(1,7)*H_invDmH(5,1)*H_invDmH(7,3)+X(1,7)*H_invDmH(5&
&,3)*H_invDmH(7,1)+X(1,6)*H_invDmH(5,1)*H_invDmH(6,3)+X(1,6)*H_inv&
&DmH(5,3)*H_invDmH(6,1)+2*X(1,5)*H_invDmH(5,1)*H_invDmH(5,3)+X(1,4&
&)*H_invDmH(4,1)*H_invDmH(5,3)+X(1,3)*H_invDmH(3,1)*H_invDmH(5,3)+&
&X(1,2)*H_invDmH(2,1)*H_invDmH(5,3)+H_invDmH(1,1)*X(1,1)*H_invDmH(&
&5,3)+X(1,4)*H_invDmH(4,3)*H_invDmH(5,1)+X(1,3)*H_invDmH(3,3)*H_in&
&vDmH(5,1)+X(1,2)*H_invDmH(2,3)*H_invDmH(5,1)+X(1,1)*H_invDmH(1,3)&
&*H_invDmH(5,1)
tt165 = X(1,8)*H_invDmH(5,2)*H_invDmH(8,3)+X(1,8)*H_invDmH(5,3)*H&
&_invDmH(8,2)+X(1,7)*H_invDmH(5,2)*H_invDmH(7,3)+X(1,7)*H_invDmH(5&
&,3)*H_invDmH(7,2)+X(1,6)*H_invDmH(5,2)*H_invDmH(6,3)+X(1,6)*H_inv&
&DmH(5,3)*H_invDmH(6,2)+2*X(1,5)*H_invDmH(5,2)*H_invDmH(5,3)+X(1,4&
&)*H_invDmH(4,2)*H_invDmH(5,3)+X(1,3)*H_invDmH(3,2)*H_invDmH(5,3)+&
&X(1,2)*H_invDmH(2,2)*H_invDmH(5,3)+X(1,1)*H_invDmH(1,2)*H_invDmH(&
&5,3)+X(1,4)*H_invDmH(4,3)*H_invDmH(5,2)+X(1,3)*H_invDmH(3,3)*H_in&
&vDmH(5,2)+X(1,2)*H_invDmH(2,3)*H_invDmH(5,2)+X(1,1)*H_invDmH(1,3)&
&*H_invDmH(5,2)
tt166 = (H_invDmH(1,1)*H_invDmH(5,3)+H_invDmH(1,3)*H_invDmH(5,1))&
&*tt56
tt167 = (H_invDmH(1,2)*H_invDmH(5,3)+H_invDmH(1,3)*H_invDmH(5,2))&
&*tt58
tt168 = H_invDmH(1,3)*H_invDmH(5,3)*tt24
tt169 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt168+2*H_invDmH(1,3)*H_inv&
&DmH(5,3)*tt19+tt167+tt166+tt61*tt165+tt60*tt164+tt163+2*H_invDmH(&
&1,2)*H_invDmH(5,2)*tt12+tt162+tt54*tt161+tt160+2*H_invDmH(1,1)*H_&
&invDmH(5,1)*tt5)+tt159+lam(1,1)*tt63*tt158)
tt170 = H_invDmH(5,3)*tt20+H_invDmH(5,2)*tt13+H_invDmH(5,1)*tt6
tt171 = 2*H_invDmH(1,1)*H_invDmH(5,1)*tt4*tt6
tt172 = 2*H_invDmH(1,2)*H_invDmH(5,2)*tt11*tt13
tt173 = X(2,8)*H_invDmH(5,1)*H_invDmH(8,2)+X(2,8)*H_invDmH(5,2)*H&
&_invDmH(8,1)+X(2,7)*H_invDmH(5,1)*H_invDmH(7,2)+X(2,7)*H_invDmH(5&
&,2)*H_invDmH(7,1)+X(2,6)*H_invDmH(5,1)*H_invDmH(6,2)+X(2,6)*H_inv&
&DmH(5,2)*H_invDmH(6,1)+2*X(2,5)*H_invDmH(5,1)*H_invDmH(5,2)+X(2,4&
&)*H_invDmH(4,1)*H_invDmH(5,2)+X(2,3)*H_invDmH(3,1)*H_invDmH(5,2)+&
&H_invDmH(2,1)*X(2,2)*H_invDmH(5,2)+H_invDmH(1,1)*X(2,1)*H_invDmH(&
&5,2)+X(2,4)*H_invDmH(4,2)*H_invDmH(5,1)+X(2,3)*H_invDmH(3,2)*H_in&
&vDmH(5,1)+H_invDmH(2,2)*X(2,2)*H_invDmH(5,1)+H_invDmH(1,2)*X(2,1)&
&*H_invDmH(5,1)
tt174 = 2*H_invDmH(1,3)*H_invDmH(5,3)*tt18*tt20
tt175 = X(2,8)*H_invDmH(5,1)*H_invDmH(8,3)+X(2,8)*H_invDmH(5,3)*H&
&_invDmH(8,1)+X(2,7)*H_invDmH(5,1)*H_invDmH(7,3)+X(2,7)*H_invDmH(5&
&,3)*H_invDmH(7,1)+X(2,6)*H_invDmH(5,1)*H_invDmH(6,3)+X(2,6)*H_inv&
&DmH(5,3)*H_invDmH(6,1)+2*X(2,5)*H_invDmH(5,1)*H_invDmH(5,3)+X(2,4&
&)*H_invDmH(4,1)*H_invDmH(5,3)+X(2,3)*H_invDmH(3,1)*H_invDmH(5,3)+&
&H_invDmH(2,1)*X(2,2)*H_invDmH(5,3)+H_invDmH(1,1)*X(2,1)*H_invDmH(&
&5,3)+X(2,4)*H_invDmH(4,3)*H_invDmH(5,1)+X(2,3)*H_invDmH(3,3)*H_in&
&vDmH(5,1)+X(2,2)*H_invDmH(2,3)*H_invDmH(5,1)+H_invDmH(1,3)*X(2,1)&
&*H_invDmH(5,1)
tt176 = X(2,8)*H_invDmH(5,2)*H_invDmH(8,3)+X(2,8)*H_invDmH(5,3)*H&
&_invDmH(8,2)+X(2,7)*H_invDmH(5,2)*H_invDmH(7,3)+X(2,7)*H_invDmH(5&
&,3)*H_invDmH(7,2)+X(2,6)*H_invDmH(5,2)*H_invDmH(6,3)+X(2,6)*H_inv&
&DmH(5,3)*H_invDmH(6,2)+2*X(2,5)*H_invDmH(5,2)*H_invDmH(5,3)+X(2,4&
&)*H_invDmH(4,2)*H_invDmH(5,3)+X(2,3)*H_invDmH(3,2)*H_invDmH(5,3)+&
&H_invDmH(2,2)*X(2,2)*H_invDmH(5,3)+H_invDmH(1,2)*X(2,1)*H_invDmH(&
&5,3)+X(2,4)*H_invDmH(4,3)*H_invDmH(5,2)+X(2,3)*H_invDmH(3,3)*H_in&
&vDmH(5,2)+X(2,2)*H_invDmH(2,3)*H_invDmH(5,2)+H_invDmH(1,3)*X(2,1)&
&*H_invDmH(5,2)
tt177 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt61*tt176+tt60*tt175+tt174&
&+tt54*tt173+tt172+tt171)+lam(1,1)*tt63*tt170)
tt178 = H_invDmH(5,3)*tt22+H_invDmH(5,2)*tt15+H_invDmH(5,1)*tt8
tt179 = 2*H_invDmH(1,1)*H_invDmH(5,1)*tt4*tt8
tt180 = 2*H_invDmH(1,2)*H_invDmH(5,2)*tt11*tt15
tt181 = X(3,8)*H_invDmH(5,1)*H_invDmH(8,2)+X(3,8)*H_invDmH(5,2)*H&
&_invDmH(8,1)+X(3,7)*H_invDmH(5,1)*H_invDmH(7,2)+X(3,7)*H_invDmH(5&
&,2)*H_invDmH(7,1)+X(3,6)*H_invDmH(5,1)*H_invDmH(6,2)+X(3,6)*H_inv&
&DmH(5,2)*H_invDmH(6,1)+2*X(3,5)*H_invDmH(5,1)*H_invDmH(5,2)+X(3,4&
&)*H_invDmH(4,1)*H_invDmH(5,2)+H_invDmH(3,1)*X(3,3)*H_invDmH(5,2)+&
&H_invDmH(2,1)*X(3,2)*H_invDmH(5,2)+H_invDmH(1,1)*X(3,1)*H_invDmH(&
&5,2)+X(3,4)*H_invDmH(4,2)*H_invDmH(5,1)+H_invDmH(3,2)*X(3,3)*H_in&
&vDmH(5,1)+H_invDmH(2,2)*X(3,2)*H_invDmH(5,1)+H_invDmH(1,2)*X(3,1)&
&*H_invDmH(5,1)
tt182 = 2*H_invDmH(1,3)*H_invDmH(5,3)*tt18*tt22
tt183 = X(3,8)*H_invDmH(5,1)*H_invDmH(8,3)+X(3,8)*H_invDmH(5,3)*H&
&_invDmH(8,1)+X(3,7)*H_invDmH(5,1)*H_invDmH(7,3)+X(3,7)*H_invDmH(5&
&,3)*H_invDmH(7,1)+X(3,6)*H_invDmH(5,1)*H_invDmH(6,3)+X(3,6)*H_inv&
&DmH(5,3)*H_invDmH(6,1)+2*X(3,5)*H_invDmH(5,1)*H_invDmH(5,3)+X(3,4&
&)*H_invDmH(4,1)*H_invDmH(5,3)+H_invDmH(3,1)*X(3,3)*H_invDmH(5,3)+&
&H_invDmH(2,1)*X(3,2)*H_invDmH(5,3)+H_invDmH(1,1)*X(3,1)*H_invDmH(&
&5,3)+X(3,4)*H_invDmH(4,3)*H_invDmH(5,1)+H_invDmH(3,3)*X(3,3)*H_in&
&vDmH(5,1)+H_invDmH(2,3)*X(3,2)*H_invDmH(5,1)+H_invDmH(1,3)*X(3,1)&
&*H_invDmH(5,1)
tt184 = X(3,8)*H_invDmH(5,2)*H_invDmH(8,3)+X(3,8)*H_invDmH(5,3)*H&
&_invDmH(8,2)+X(3,7)*H_invDmH(5,2)*H_invDmH(7,3)+X(3,7)*H_invDmH(5&
&,3)*H_invDmH(7,2)+X(3,6)*H_invDmH(5,2)*H_invDmH(6,3)+X(3,6)*H_inv&
&DmH(5,3)*H_invDmH(6,2)+2*X(3,5)*H_invDmH(5,2)*H_invDmH(5,3)+X(3,4&
&)*H_invDmH(4,2)*H_invDmH(5,3)+H_invDmH(3,2)*X(3,3)*H_invDmH(5,3)+&
&H_invDmH(2,2)*X(3,2)*H_invDmH(5,3)+H_invDmH(1,2)*X(3,1)*H_invDmH(&
&5,3)+X(3,4)*H_invDmH(4,3)*H_invDmH(5,2)+H_invDmH(3,3)*X(3,3)*H_in&
&vDmH(5,2)+H_invDmH(2,3)*X(3,2)*H_invDmH(5,2)+H_invDmH(1,3)*X(3,1)&
&*H_invDmH(5,2)
tt185 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt61*tt184+tt60*tt183+tt182&
&+tt54*tt181+tt180+tt179)+lam(1,1)*tt63*tt178)
tt186 = H_invDmH(6,3)*tt18+H_invDmH(6,2)*tt11+H_invDmH(6,1)*tt4
tt187 = lam(1,1)*(H_invDmH(1,3)*H_invDmH(6,3)+H_invDmH(1,2)*H_inv&
&DmH(6,2)+H_invDmH(1,1)*H_invDmH(6,1))*tt25
tt188 = H_invDmH(1,1)*H_invDmH(6,1)*tt10
tt189 = X(1,8)*H_invDmH(6,1)*H_invDmH(8,2)+X(1,8)*H_invDmH(6,2)*H&
&_invDmH(8,1)+X(1,7)*H_invDmH(6,1)*H_invDmH(7,2)+X(1,7)*H_invDmH(6&
&,2)*H_invDmH(7,1)+2*X(1,6)*H_invDmH(6,1)*H_invDmH(6,2)+X(1,5)*H_i&
&nvDmH(5,1)*H_invDmH(6,2)+X(1,4)*H_invDmH(4,1)*H_invDmH(6,2)+X(1,3&
&)*H_invDmH(3,1)*H_invDmH(6,2)+X(1,2)*H_invDmH(2,1)*H_invDmH(6,2)+&
&H_invDmH(1,1)*X(1,1)*H_invDmH(6,2)+X(1,5)*H_invDmH(5,2)*H_invDmH(&
&6,1)+X(1,4)*H_invDmH(4,2)*H_invDmH(6,1)+X(1,3)*H_invDmH(3,2)*H_in&
&vDmH(6,1)+X(1,2)*H_invDmH(2,2)*H_invDmH(6,1)+X(1,1)*H_invDmH(1,2)&
&*H_invDmH(6,1)
tt190 = (H_invDmH(1,1)*H_invDmH(6,2)+H_invDmH(1,2)*H_invDmH(6,1))&
&*tt52
tt191 = H_invDmH(1,2)*H_invDmH(6,2)*tt17
tt192 = X(1,8)*H_invDmH(6,1)*H_invDmH(8,3)+X(1,8)*H_invDmH(6,3)*H&
&_invDmH(8,1)+X(1,7)*H_invDmH(6,1)*H_invDmH(7,3)+X(1,7)*H_invDmH(6&
&,3)*H_invDmH(7,1)+2*X(1,6)*H_invDmH(6,1)*H_invDmH(6,3)+X(1,5)*H_i&
&nvDmH(5,1)*H_invDmH(6,3)+X(1,4)*H_invDmH(4,1)*H_invDmH(6,3)+X(1,3&
&)*H_invDmH(3,1)*H_invDmH(6,3)+X(1,2)*H_invDmH(2,1)*H_invDmH(6,3)+&
&H_invDmH(1,1)*X(1,1)*H_invDmH(6,3)+X(1,5)*H_invDmH(5,3)*H_invDmH(&
&6,1)+X(1,4)*H_invDmH(4,3)*H_invDmH(6,1)+X(1,3)*H_invDmH(3,3)*H_in&
&vDmH(6,1)+X(1,2)*H_invDmH(2,3)*H_invDmH(6,1)+X(1,1)*H_invDmH(1,3)&
&*H_invDmH(6,1)
tt193 = X(1,8)*H_invDmH(6,2)*H_invDmH(8,3)+X(1,8)*H_invDmH(6,3)*H&
&_invDmH(8,2)+X(1,7)*H_invDmH(6,2)*H_invDmH(7,3)+X(1,7)*H_invDmH(6&
&,3)*H_invDmH(7,2)+2*X(1,6)*H_invDmH(6,2)*H_invDmH(6,3)+X(1,5)*H_i&
&nvDmH(5,2)*H_invDmH(6,3)+X(1,4)*H_invDmH(4,2)*H_invDmH(6,3)+X(1,3&
&)*H_invDmH(3,2)*H_invDmH(6,3)+X(1,2)*H_invDmH(2,2)*H_invDmH(6,3)+&
&X(1,1)*H_invDmH(1,2)*H_invDmH(6,3)+X(1,5)*H_invDmH(5,3)*H_invDmH(&
&6,2)+X(1,4)*H_invDmH(4,3)*H_invDmH(6,2)+X(1,3)*H_invDmH(3,3)*H_in&
&vDmH(6,2)+X(1,2)*H_invDmH(2,3)*H_invDmH(6,2)+X(1,1)*H_invDmH(1,3)&
&*H_invDmH(6,2)
tt194 = (H_invDmH(1,1)*H_invDmH(6,3)+H_invDmH(1,3)*H_invDmH(6,1))&
&*tt56
tt195 = (H_invDmH(1,2)*H_invDmH(6,3)+H_invDmH(1,3)*H_invDmH(6,2))&
&*tt58
tt196 = H_invDmH(1,3)*H_invDmH(6,3)*tt24
tt197 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt196+2*H_invDmH(1,3)*H_inv&
&DmH(6,3)*tt19+tt195+tt194+tt61*tt193+tt60*tt192+tt191+2*H_invDmH(&
&1,2)*H_invDmH(6,2)*tt12+tt190+tt54*tt189+tt188+2*H_invDmH(1,1)*H_&
&invDmH(6,1)*tt5)+tt187+lam(1,1)*tt63*tt186)
tt198 = H_invDmH(6,3)*tt20+H_invDmH(6,2)*tt13+H_invDmH(6,1)*tt6
tt199 = 2*H_invDmH(1,1)*H_invDmH(6,1)*tt4*tt6
tt200 = 2*H_invDmH(1,2)*H_invDmH(6,2)*tt11*tt13
tt201 = X(2,8)*H_invDmH(6,1)*H_invDmH(8,2)+X(2,8)*H_invDmH(6,2)*H&
&_invDmH(8,1)+X(2,7)*H_invDmH(6,1)*H_invDmH(7,2)+X(2,7)*H_invDmH(6&
&,2)*H_invDmH(7,1)+2*X(2,6)*H_invDmH(6,1)*H_invDmH(6,2)+X(2,5)*H_i&
&nvDmH(5,1)*H_invDmH(6,2)+X(2,4)*H_invDmH(4,1)*H_invDmH(6,2)+X(2,3&
&)*H_invDmH(3,1)*H_invDmH(6,2)+H_invDmH(2,1)*X(2,2)*H_invDmH(6,2)+&
&H_invDmH(1,1)*X(2,1)*H_invDmH(6,2)+X(2,5)*H_invDmH(5,2)*H_invDmH(&
&6,1)+X(2,4)*H_invDmH(4,2)*H_invDmH(6,1)+X(2,3)*H_invDmH(3,2)*H_in&
&vDmH(6,1)+H_invDmH(2,2)*X(2,2)*H_invDmH(6,1)+H_invDmH(1,2)*X(2,1)&
&*H_invDmH(6,1)
tt202 = 2*H_invDmH(1,3)*H_invDmH(6,3)*tt18*tt20
tt203 = X(2,8)*H_invDmH(6,1)*H_invDmH(8,3)+X(2,8)*H_invDmH(6,3)*H&
&_invDmH(8,1)+X(2,7)*H_invDmH(6,1)*H_invDmH(7,3)+X(2,7)*H_invDmH(6&
&,3)*H_invDmH(7,1)+2*X(2,6)*H_invDmH(6,1)*H_invDmH(6,3)+X(2,5)*H_i&
&nvDmH(5,1)*H_invDmH(6,3)+X(2,4)*H_invDmH(4,1)*H_invDmH(6,3)+X(2,3&
&)*H_invDmH(3,1)*H_invDmH(6,3)+H_invDmH(2,1)*X(2,2)*H_invDmH(6,3)+&
&H_invDmH(1,1)*X(2,1)*H_invDmH(6,3)+X(2,5)*H_invDmH(5,3)*H_invDmH(&
&6,1)+X(2,4)*H_invDmH(4,3)*H_invDmH(6,1)+X(2,3)*H_invDmH(3,3)*H_in&
&vDmH(6,1)+X(2,2)*H_invDmH(2,3)*H_invDmH(6,1)+H_invDmH(1,3)*X(2,1)&
&*H_invDmH(6,1)
tt204 = X(2,8)*H_invDmH(6,2)*H_invDmH(8,3)+X(2,8)*H_invDmH(6,3)*H&
&_invDmH(8,2)+X(2,7)*H_invDmH(6,2)*H_invDmH(7,3)+X(2,7)*H_invDmH(6&
&,3)*H_invDmH(7,2)+2*X(2,6)*H_invDmH(6,2)*H_invDmH(6,3)+X(2,5)*H_i&
&nvDmH(5,2)*H_invDmH(6,3)+X(2,4)*H_invDmH(4,2)*H_invDmH(6,3)+X(2,3&
&)*H_invDmH(3,2)*H_invDmH(6,3)+H_invDmH(2,2)*X(2,2)*H_invDmH(6,3)+&
&H_invDmH(1,2)*X(2,1)*H_invDmH(6,3)+X(2,5)*H_invDmH(5,3)*H_invDmH(&
&6,2)+X(2,4)*H_invDmH(4,3)*H_invDmH(6,2)+X(2,3)*H_invDmH(3,3)*H_in&
&vDmH(6,2)+X(2,2)*H_invDmH(2,3)*H_invDmH(6,2)+H_invDmH(1,3)*X(2,1)&
&*H_invDmH(6,2)
tt205 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt61*tt204+tt60*tt203+tt202&
&+tt54*tt201+tt200+tt199)+lam(1,1)*tt63*tt198)
tt206 = H_invDmH(6,3)*tt22+H_invDmH(6,2)*tt15+H_invDmH(6,1)*tt8
tt207 = 2*H_invDmH(1,1)*H_invDmH(6,1)*tt4*tt8
tt208 = 2*H_invDmH(1,2)*H_invDmH(6,2)*tt11*tt15
tt209 = X(3,8)*H_invDmH(6,1)*H_invDmH(8,2)+X(3,8)*H_invDmH(6,2)*H&
&_invDmH(8,1)+X(3,7)*H_invDmH(6,1)*H_invDmH(7,2)+X(3,7)*H_invDmH(6&
&,2)*H_invDmH(7,1)+2*X(3,6)*H_invDmH(6,1)*H_invDmH(6,2)+X(3,5)*H_i&
&nvDmH(5,1)*H_invDmH(6,2)+X(3,4)*H_invDmH(4,1)*H_invDmH(6,2)+H_inv&
&DmH(3,1)*X(3,3)*H_invDmH(6,2)+H_invDmH(2,1)*X(3,2)*H_invDmH(6,2)+&
&H_invDmH(1,1)*X(3,1)*H_invDmH(6,2)+X(3,5)*H_invDmH(5,2)*H_invDmH(&
&6,1)+X(3,4)*H_invDmH(4,2)*H_invDmH(6,1)+H_invDmH(3,2)*X(3,3)*H_in&
&vDmH(6,1)+H_invDmH(2,2)*X(3,2)*H_invDmH(6,1)+H_invDmH(1,2)*X(3,1)&
&*H_invDmH(6,1)
tt210 = 2*H_invDmH(1,3)*H_invDmH(6,3)*tt18*tt22
tt211 = X(3,8)*H_invDmH(6,1)*H_invDmH(8,3)+X(3,8)*H_invDmH(6,3)*H&
&_invDmH(8,1)+X(3,7)*H_invDmH(6,1)*H_invDmH(7,3)+X(3,7)*H_invDmH(6&
&,3)*H_invDmH(7,1)+2*X(3,6)*H_invDmH(6,1)*H_invDmH(6,3)+X(3,5)*H_i&
&nvDmH(5,1)*H_invDmH(6,3)+X(3,4)*H_invDmH(4,1)*H_invDmH(6,3)+H_inv&
&DmH(3,1)*X(3,3)*H_invDmH(6,3)+H_invDmH(2,1)*X(3,2)*H_invDmH(6,3)+&
&H_invDmH(1,1)*X(3,1)*H_invDmH(6,3)+X(3,5)*H_invDmH(5,3)*H_invDmH(&
&6,1)+X(3,4)*H_invDmH(4,3)*H_invDmH(6,1)+H_invDmH(3,3)*X(3,3)*H_in&
&vDmH(6,1)+H_invDmH(2,3)*X(3,2)*H_invDmH(6,1)+H_invDmH(1,3)*X(3,1)&
&*H_invDmH(6,1)
tt212 = X(3,8)*H_invDmH(6,2)*H_invDmH(8,3)+X(3,8)*H_invDmH(6,3)*H&
&_invDmH(8,2)+X(3,7)*H_invDmH(6,2)*H_invDmH(7,3)+X(3,7)*H_invDmH(6&
&,3)*H_invDmH(7,2)+2*X(3,6)*H_invDmH(6,2)*H_invDmH(6,3)+X(3,5)*H_i&
&nvDmH(5,2)*H_invDmH(6,3)+X(3,4)*H_invDmH(4,2)*H_invDmH(6,3)+H_inv&
&DmH(3,2)*X(3,3)*H_invDmH(6,3)+H_invDmH(2,2)*X(3,2)*H_invDmH(6,3)+&
&H_invDmH(1,2)*X(3,1)*H_invDmH(6,3)+X(3,5)*H_invDmH(5,3)*H_invDmH(&
&6,2)+X(3,4)*H_invDmH(4,3)*H_invDmH(6,2)+H_invDmH(3,3)*X(3,3)*H_in&
&vDmH(6,2)+H_invDmH(2,3)*X(3,2)*H_invDmH(6,2)+H_invDmH(1,3)*X(3,1)&
&*H_invDmH(6,2)
tt213 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt61*tt212+tt60*tt211+tt210&
&+tt54*tt209+tt208+tt207)+lam(1,1)*tt63*tt206)
tt214 = H_invDmH(7,3)*tt18+H_invDmH(7,2)*tt11+H_invDmH(7,1)*tt4
tt215 = lam(1,1)*(H_invDmH(1,3)*H_invDmH(7,3)+H_invDmH(1,2)*H_inv&
&DmH(7,2)+H_invDmH(1,1)*H_invDmH(7,1))*tt25
tt216 = H_invDmH(1,1)*H_invDmH(7,1)*tt10
tt217 = X(1,8)*H_invDmH(7,1)*H_invDmH(8,2)+X(1,8)*H_invDmH(7,2)*H&
&_invDmH(8,1)+2*X(1,7)*H_invDmH(7,1)*H_invDmH(7,2)+X(1,6)*H_invDmH&
&(6,1)*H_invDmH(7,2)+X(1,5)*H_invDmH(5,1)*H_invDmH(7,2)+X(1,4)*H_i&
&nvDmH(4,1)*H_invDmH(7,2)+X(1,3)*H_invDmH(3,1)*H_invDmH(7,2)+X(1,2&
&)*H_invDmH(2,1)*H_invDmH(7,2)+H_invDmH(1,1)*X(1,1)*H_invDmH(7,2)+&
&X(1,6)*H_invDmH(6,2)*H_invDmH(7,1)+X(1,5)*H_invDmH(5,2)*H_invDmH(&
&7,1)+X(1,4)*H_invDmH(4,2)*H_invDmH(7,1)+X(1,3)*H_invDmH(3,2)*H_in&
&vDmH(7,1)+X(1,2)*H_invDmH(2,2)*H_invDmH(7,1)+X(1,1)*H_invDmH(1,2)&
&*H_invDmH(7,1)
tt218 = (H_invDmH(1,1)*H_invDmH(7,2)+H_invDmH(1,2)*H_invDmH(7,1))&
&*tt52
tt219 = H_invDmH(1,2)*H_invDmH(7,2)*tt17
tt220 = X(1,8)*H_invDmH(7,1)*H_invDmH(8,3)+X(1,8)*H_invDmH(7,3)*H&
&_invDmH(8,1)+2*X(1,7)*H_invDmH(7,1)*H_invDmH(7,3)+X(1,6)*H_invDmH&
&(6,1)*H_invDmH(7,3)+X(1,5)*H_invDmH(5,1)*H_invDmH(7,3)+X(1,4)*H_i&
&nvDmH(4,1)*H_invDmH(7,3)+X(1,3)*H_invDmH(3,1)*H_invDmH(7,3)+X(1,2&
&)*H_invDmH(2,1)*H_invDmH(7,3)+H_invDmH(1,1)*X(1,1)*H_invDmH(7,3)+&
&X(1,6)*H_invDmH(6,3)*H_invDmH(7,1)+X(1,5)*H_invDmH(5,3)*H_invDmH(&
&7,1)+X(1,4)*H_invDmH(4,3)*H_invDmH(7,1)+X(1,3)*H_invDmH(3,3)*H_in&
&vDmH(7,1)+X(1,2)*H_invDmH(2,3)*H_invDmH(7,1)+X(1,1)*H_invDmH(1,3)&
&*H_invDmH(7,1)
tt221 = X(1,8)*H_invDmH(7,2)*H_invDmH(8,3)+X(1,8)*H_invDmH(7,3)*H&
&_invDmH(8,2)+2*X(1,7)*H_invDmH(7,2)*H_invDmH(7,3)+X(1,6)*H_invDmH&
&(6,2)*H_invDmH(7,3)+X(1,5)*H_invDmH(5,2)*H_invDmH(7,3)+X(1,4)*H_i&
&nvDmH(4,2)*H_invDmH(7,3)+X(1,3)*H_invDmH(3,2)*H_invDmH(7,3)+X(1,2&
&)*H_invDmH(2,2)*H_invDmH(7,3)+X(1,1)*H_invDmH(1,2)*H_invDmH(7,3)+&
&X(1,6)*H_invDmH(6,3)*H_invDmH(7,2)+X(1,5)*H_invDmH(5,3)*H_invDmH(&
&7,2)+X(1,4)*H_invDmH(4,3)*H_invDmH(7,2)+X(1,3)*H_invDmH(3,3)*H_in&
&vDmH(7,2)+X(1,2)*H_invDmH(2,3)*H_invDmH(7,2)+X(1,1)*H_invDmH(1,3)&
&*H_invDmH(7,2)
tt222 = (H_invDmH(1,1)*H_invDmH(7,3)+H_invDmH(1,3)*H_invDmH(7,1))&
&*tt56
tt223 = (H_invDmH(1,2)*H_invDmH(7,3)+H_invDmH(1,3)*H_invDmH(7,2))&
&*tt58
tt224 = H_invDmH(1,3)*H_invDmH(7,3)*tt24
tt225 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt224+2*H_invDmH(1,3)*H_inv&
&DmH(7,3)*tt19+tt223+tt222+tt61*tt221+tt60*tt220+tt219+2*H_invDmH(&
&1,2)*H_invDmH(7,2)*tt12+tt218+tt54*tt217+tt216+2*H_invDmH(1,1)*H_&
&invDmH(7,1)*tt5)+tt215+lam(1,1)*tt63*tt214)
tt226 = H_invDmH(7,3)*tt20+H_invDmH(7,2)*tt13+H_invDmH(7,1)*tt6
tt227 = 2*H_invDmH(1,1)*H_invDmH(7,1)*tt4*tt6
tt228 = 2*H_invDmH(1,2)*H_invDmH(7,2)*tt11*tt13
tt229 = X(2,8)*H_invDmH(7,1)*H_invDmH(8,2)+X(2,8)*H_invDmH(7,2)*H&
&_invDmH(8,1)+2*X(2,7)*H_invDmH(7,1)*H_invDmH(7,2)+X(2,6)*H_invDmH&
&(6,1)*H_invDmH(7,2)+X(2,5)*H_invDmH(5,1)*H_invDmH(7,2)+X(2,4)*H_i&
&nvDmH(4,1)*H_invDmH(7,2)+X(2,3)*H_invDmH(3,1)*H_invDmH(7,2)+H_inv&
&DmH(2,1)*X(2,2)*H_invDmH(7,2)+H_invDmH(1,1)*X(2,1)*H_invDmH(7,2)+&
&X(2,6)*H_invDmH(6,2)*H_invDmH(7,1)+X(2,5)*H_invDmH(5,2)*H_invDmH(&
&7,1)+X(2,4)*H_invDmH(4,2)*H_invDmH(7,1)+X(2,3)*H_invDmH(3,2)*H_in&
&vDmH(7,1)+H_invDmH(2,2)*X(2,2)*H_invDmH(7,1)+H_invDmH(1,2)*X(2,1)&
&*H_invDmH(7,1)
tt230 = 2*H_invDmH(1,3)*H_invDmH(7,3)*tt18*tt20
tt231 = X(2,8)*H_invDmH(7,1)*H_invDmH(8,3)+X(2,8)*H_invDmH(7,3)*H&
&_invDmH(8,1)+2*X(2,7)*H_invDmH(7,1)*H_invDmH(7,3)+X(2,6)*H_invDmH&
&(6,1)*H_invDmH(7,3)+X(2,5)*H_invDmH(5,1)*H_invDmH(7,3)+X(2,4)*H_i&
&nvDmH(4,1)*H_invDmH(7,3)+X(2,3)*H_invDmH(3,1)*H_invDmH(7,3)+H_inv&
&DmH(2,1)*X(2,2)*H_invDmH(7,3)+H_invDmH(1,1)*X(2,1)*H_invDmH(7,3)+&
&X(2,6)*H_invDmH(6,3)*H_invDmH(7,1)+X(2,5)*H_invDmH(5,3)*H_invDmH(&
&7,1)+X(2,4)*H_invDmH(4,3)*H_invDmH(7,1)+X(2,3)*H_invDmH(3,3)*H_in&
&vDmH(7,1)+X(2,2)*H_invDmH(2,3)*H_invDmH(7,1)+H_invDmH(1,3)*X(2,1)&
&*H_invDmH(7,1)
tt232 = X(2,8)*H_invDmH(7,2)*H_invDmH(8,3)+X(2,8)*H_invDmH(7,3)*H&
&_invDmH(8,2)+2*X(2,7)*H_invDmH(7,2)*H_invDmH(7,3)+X(2,6)*H_invDmH&
&(6,2)*H_invDmH(7,3)+X(2,5)*H_invDmH(5,2)*H_invDmH(7,3)+X(2,4)*H_i&
&nvDmH(4,2)*H_invDmH(7,3)+X(2,3)*H_invDmH(3,2)*H_invDmH(7,3)+H_inv&
&DmH(2,2)*X(2,2)*H_invDmH(7,3)+H_invDmH(1,2)*X(2,1)*H_invDmH(7,3)+&
&X(2,6)*H_invDmH(6,3)*H_invDmH(7,2)+X(2,5)*H_invDmH(5,3)*H_invDmH(&
&7,2)+X(2,4)*H_invDmH(4,3)*H_invDmH(7,2)+X(2,3)*H_invDmH(3,3)*H_in&
&vDmH(7,2)+X(2,2)*H_invDmH(2,3)*H_invDmH(7,2)+H_invDmH(1,3)*X(2,1)&
&*H_invDmH(7,2)
tt233 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt61*tt232+tt60*tt231+tt230&
&+tt54*tt229+tt228+tt227)+lam(1,1)*tt63*tt226)
tt234 = H_invDmH(7,3)*tt22+H_invDmH(7,2)*tt15+H_invDmH(7,1)*tt8
tt235 = 2*H_invDmH(1,1)*H_invDmH(7,1)*tt4*tt8
tt236 = 2*H_invDmH(1,2)*H_invDmH(7,2)*tt11*tt15
tt237 = X(3,8)*H_invDmH(7,1)*H_invDmH(8,2)+X(3,8)*H_invDmH(7,2)*H&
&_invDmH(8,1)+2*X(3,7)*H_invDmH(7,1)*H_invDmH(7,2)+X(3,6)*H_invDmH&
&(6,1)*H_invDmH(7,2)+X(3,5)*H_invDmH(5,1)*H_invDmH(7,2)+X(3,4)*H_i&
&nvDmH(4,1)*H_invDmH(7,2)+H_invDmH(3,1)*X(3,3)*H_invDmH(7,2)+H_inv&
&DmH(2,1)*X(3,2)*H_invDmH(7,2)+H_invDmH(1,1)*X(3,1)*H_invDmH(7,2)+&
&X(3,6)*H_invDmH(6,2)*H_invDmH(7,1)+X(3,5)*H_invDmH(5,2)*H_invDmH(&
&7,1)+X(3,4)*H_invDmH(4,2)*H_invDmH(7,1)+H_invDmH(3,2)*X(3,3)*H_in&
&vDmH(7,1)+H_invDmH(2,2)*X(3,2)*H_invDmH(7,1)+H_invDmH(1,2)*X(3,1)&
&*H_invDmH(7,1)
tt238 = 2*H_invDmH(1,3)*H_invDmH(7,3)*tt18*tt22
tt239 = X(3,8)*H_invDmH(7,1)*H_invDmH(8,3)+X(3,8)*H_invDmH(7,3)*H&
&_invDmH(8,1)+2*X(3,7)*H_invDmH(7,1)*H_invDmH(7,3)+X(3,6)*H_invDmH&
&(6,1)*H_invDmH(7,3)+X(3,5)*H_invDmH(5,1)*H_invDmH(7,3)+X(3,4)*H_i&
&nvDmH(4,1)*H_invDmH(7,3)+H_invDmH(3,1)*X(3,3)*H_invDmH(7,3)+H_inv&
&DmH(2,1)*X(3,2)*H_invDmH(7,3)+H_invDmH(1,1)*X(3,1)*H_invDmH(7,3)+&
&X(3,6)*H_invDmH(6,3)*H_invDmH(7,1)+X(3,5)*H_invDmH(5,3)*H_invDmH(&
&7,1)+X(3,4)*H_invDmH(4,3)*H_invDmH(7,1)+H_invDmH(3,3)*X(3,3)*H_in&
&vDmH(7,1)+H_invDmH(2,3)*X(3,2)*H_invDmH(7,1)+H_invDmH(1,3)*X(3,1)&
&*H_invDmH(7,1)
tt240 = X(3,8)*H_invDmH(7,2)*H_invDmH(8,3)+X(3,8)*H_invDmH(7,3)*H&
&_invDmH(8,2)+2*X(3,7)*H_invDmH(7,2)*H_invDmH(7,3)+X(3,6)*H_invDmH&
&(6,2)*H_invDmH(7,3)+X(3,5)*H_invDmH(5,2)*H_invDmH(7,3)+X(3,4)*H_i&
&nvDmH(4,2)*H_invDmH(7,3)+H_invDmH(3,2)*X(3,3)*H_invDmH(7,3)+H_inv&
&DmH(2,2)*X(3,2)*H_invDmH(7,3)+H_invDmH(1,2)*X(3,1)*H_invDmH(7,3)+&
&X(3,6)*H_invDmH(6,3)*H_invDmH(7,2)+X(3,5)*H_invDmH(5,3)*H_invDmH(&
&7,2)+X(3,4)*H_invDmH(4,3)*H_invDmH(7,2)+H_invDmH(3,3)*X(3,3)*H_in&
&vDmH(7,2)+H_invDmH(2,3)*X(3,2)*H_invDmH(7,2)+H_invDmH(1,3)*X(3,1)&
&*H_invDmH(7,2)
tt241 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt61*tt240+tt60*tt239+tt238&
&+tt54*tt237+tt236+tt235)+lam(1,1)*tt63*tt234)
tt242 = H_invDmH(8,3)*tt18+H_invDmH(8,2)*tt11+H_invDmH(8,1)*tt4
tt243 = lam(1,1)*(H_invDmH(1,3)*H_invDmH(8,3)+H_invDmH(1,2)*H_inv&
&DmH(8,2)+H_invDmH(1,1)*H_invDmH(8,1))*tt25
tt244 = H_invDmH(1,1)*H_invDmH(8,1)*tt10
tt245 = 2*X(1,8)*H_invDmH(8,1)*H_invDmH(8,2)+X(1,7)*H_invDmH(7,1)&
&*H_invDmH(8,2)+X(1,6)*H_invDmH(6,1)*H_invDmH(8,2)+X(1,5)*H_invDmH&
&(5,1)*H_invDmH(8,2)+X(1,4)*H_invDmH(4,1)*H_invDmH(8,2)+X(1,3)*H_i&
&nvDmH(3,1)*H_invDmH(8,2)+X(1,2)*H_invDmH(2,1)*H_invDmH(8,2)+H_inv&
&DmH(1,1)*X(1,1)*H_invDmH(8,2)+X(1,7)*H_invDmH(7,2)*H_invDmH(8,1)+&
&X(1,6)*H_invDmH(6,2)*H_invDmH(8,1)+X(1,5)*H_invDmH(5,2)*H_invDmH(&
&8,1)+X(1,4)*H_invDmH(4,2)*H_invDmH(8,1)+X(1,3)*H_invDmH(3,2)*H_in&
&vDmH(8,1)+X(1,2)*H_invDmH(2,2)*H_invDmH(8,1)+X(1,1)*H_invDmH(1,2)&
&*H_invDmH(8,1)
tt246 = (H_invDmH(1,1)*H_invDmH(8,2)+H_invDmH(1,2)*H_invDmH(8,1))&
&*tt52
tt247 = H_invDmH(1,2)*H_invDmH(8,2)*tt17
tt248 = 2*X(1,8)*H_invDmH(8,1)*H_invDmH(8,3)+X(1,7)*H_invDmH(7,1)&
&*H_invDmH(8,3)+X(1,6)*H_invDmH(6,1)*H_invDmH(8,3)+X(1,5)*H_invDmH&
&(5,1)*H_invDmH(8,3)+X(1,4)*H_invDmH(4,1)*H_invDmH(8,3)+X(1,3)*H_i&
&nvDmH(3,1)*H_invDmH(8,3)+X(1,2)*H_invDmH(2,1)*H_invDmH(8,3)+H_inv&
&DmH(1,1)*X(1,1)*H_invDmH(8,3)+X(1,7)*H_invDmH(7,3)*H_invDmH(8,1)+&
&X(1,6)*H_invDmH(6,3)*H_invDmH(8,1)+X(1,5)*H_invDmH(5,3)*H_invDmH(&
&8,1)+X(1,4)*H_invDmH(4,3)*H_invDmH(8,1)+X(1,3)*H_invDmH(3,3)*H_in&
&vDmH(8,1)+X(1,2)*H_invDmH(2,3)*H_invDmH(8,1)+X(1,1)*H_invDmH(1,3)&
&*H_invDmH(8,1)
tt249 = (H_invDmH(1,1)*H_invDmH(8,3)+H_invDmH(1,3)*H_invDmH(8,1))&
&*tt56
tt250 = 2*X(1,8)*H_invDmH(8,2)*H_invDmH(8,3)+X(1,7)*H_invDmH(7,2)&
&*H_invDmH(8,3)+X(1,6)*H_invDmH(6,2)*H_invDmH(8,3)+X(1,5)*H_invDmH&
&(5,2)*H_invDmH(8,3)+X(1,4)*H_invDmH(4,2)*H_invDmH(8,3)+X(1,3)*H_i&
&nvDmH(3,2)*H_invDmH(8,3)+X(1,2)*H_invDmH(2,2)*H_invDmH(8,3)+X(1,1&
&)*H_invDmH(1,2)*H_invDmH(8,3)+X(1,7)*H_invDmH(7,3)*H_invDmH(8,2)+&
&X(1,6)*H_invDmH(6,3)*H_invDmH(8,2)+X(1,5)*H_invDmH(5,3)*H_invDmH(&
&8,2)+X(1,4)*H_invDmH(4,3)*H_invDmH(8,2)+X(1,3)*H_invDmH(3,3)*H_in&
&vDmH(8,2)+X(1,2)*H_invDmH(2,3)*H_invDmH(8,2)+X(1,1)*H_invDmH(1,3)&
&*H_invDmH(8,2)
tt251 = (H_invDmH(1,2)*H_invDmH(8,3)+H_invDmH(1,3)*H_invDmH(8,2))&
&*tt58
tt252 = H_invDmH(1,3)*H_invDmH(8,3)*tt24
tt253 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt252+2*H_invDmH(1,3)*H_inv&
&DmH(8,3)*tt19+tt251+tt61*tt250+tt249+tt60*tt248+tt247+2*H_invDmH(&
&1,2)*H_invDmH(8,2)*tt12+tt246+tt54*tt245+tt244+2*H_invDmH(1,1)*H_&
&invDmH(8,1)*tt5)+tt243+lam(1,1)*tt63*tt242)
tt254 = H_invDmH(8,3)*tt20+H_invDmH(8,2)*tt13+H_invDmH(8,1)*tt6
tt255 = 2*H_invDmH(1,1)*H_invDmH(8,1)*tt4*tt6
tt256 = 2*H_invDmH(1,2)*H_invDmH(8,2)*tt11*tt13
tt257 = 2*X(2,8)*H_invDmH(8,1)*H_invDmH(8,2)+X(2,7)*H_invDmH(7,1)&
&*H_invDmH(8,2)+X(2,6)*H_invDmH(6,1)*H_invDmH(8,2)+X(2,5)*H_invDmH&
&(5,1)*H_invDmH(8,2)+X(2,4)*H_invDmH(4,1)*H_invDmH(8,2)+X(2,3)*H_i&
&nvDmH(3,1)*H_invDmH(8,2)+H_invDmH(2,1)*X(2,2)*H_invDmH(8,2)+H_inv&
&DmH(1,1)*X(2,1)*H_invDmH(8,2)+X(2,7)*H_invDmH(7,2)*H_invDmH(8,1)+&
&X(2,6)*H_invDmH(6,2)*H_invDmH(8,1)+X(2,5)*H_invDmH(5,2)*H_invDmH(&
&8,1)+X(2,4)*H_invDmH(4,2)*H_invDmH(8,1)+X(2,3)*H_invDmH(3,2)*H_in&
&vDmH(8,1)+H_invDmH(2,2)*X(2,2)*H_invDmH(8,1)+H_invDmH(1,2)*X(2,1)&
&*H_invDmH(8,1)
tt258 = 2*H_invDmH(1,3)*H_invDmH(8,3)*tt18*tt20
tt259 = 2*X(2,8)*H_invDmH(8,1)*H_invDmH(8,3)+X(2,7)*H_invDmH(7,1)&
&*H_invDmH(8,3)+X(2,6)*H_invDmH(6,1)*H_invDmH(8,3)+X(2,5)*H_invDmH&
&(5,1)*H_invDmH(8,3)+X(2,4)*H_invDmH(4,1)*H_invDmH(8,3)+X(2,3)*H_i&
&nvDmH(3,1)*H_invDmH(8,3)+H_invDmH(2,1)*X(2,2)*H_invDmH(8,3)+H_inv&
&DmH(1,1)*X(2,1)*H_invDmH(8,3)+X(2,7)*H_invDmH(7,3)*H_invDmH(8,1)+&
&X(2,6)*H_invDmH(6,3)*H_invDmH(8,1)+X(2,5)*H_invDmH(5,3)*H_invDmH(&
&8,1)+X(2,4)*H_invDmH(4,3)*H_invDmH(8,1)+X(2,3)*H_invDmH(3,3)*H_in&
&vDmH(8,1)+X(2,2)*H_invDmH(2,3)*H_invDmH(8,1)+H_invDmH(1,3)*X(2,1)&
&*H_invDmH(8,1)
tt260 = 2*X(2,8)*H_invDmH(8,2)*H_invDmH(8,3)+X(2,7)*H_invDmH(7,2)&
&*H_invDmH(8,3)+X(2,6)*H_invDmH(6,2)*H_invDmH(8,3)+X(2,5)*H_invDmH&
&(5,2)*H_invDmH(8,3)+X(2,4)*H_invDmH(4,2)*H_invDmH(8,3)+X(2,3)*H_i&
&nvDmH(3,2)*H_invDmH(8,3)+H_invDmH(2,2)*X(2,2)*H_invDmH(8,3)+H_inv&
&DmH(1,2)*X(2,1)*H_invDmH(8,3)+X(2,7)*H_invDmH(7,3)*H_invDmH(8,2)+&
&X(2,6)*H_invDmH(6,3)*H_invDmH(8,2)+X(2,5)*H_invDmH(5,3)*H_invDmH(&
&8,2)+X(2,4)*H_invDmH(4,3)*H_invDmH(8,2)+X(2,3)*H_invDmH(3,3)*H_in&
&vDmH(8,2)+X(2,2)*H_invDmH(2,3)*H_invDmH(8,2)+H_invDmH(1,3)*X(2,1)&
&*H_invDmH(8,2)
tt261 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt61*tt260+tt60*tt259+tt258&
&+tt54*tt257+tt256+tt255)+lam(1,1)*tt63*tt254)
tt262 = H_invDmH(8,3)*tt22+H_invDmH(8,2)*tt15+H_invDmH(8,1)*tt8
tt263 = 2*H_invDmH(1,1)*H_invDmH(8,1)*tt4*tt8
tt264 = 2*H_invDmH(1,2)*H_invDmH(8,2)*tt11*tt15
tt265 = 2*X(3,8)*H_invDmH(8,1)*H_invDmH(8,2)+X(3,7)*H_invDmH(7,1)&
&*H_invDmH(8,2)+X(3,6)*H_invDmH(6,1)*H_invDmH(8,2)+X(3,5)*H_invDmH&
&(5,1)*H_invDmH(8,2)+X(3,4)*H_invDmH(4,1)*H_invDmH(8,2)+H_invDmH(3&
&,1)*X(3,3)*H_invDmH(8,2)+H_invDmH(2,1)*X(3,2)*H_invDmH(8,2)+H_inv&
&DmH(1,1)*X(3,1)*H_invDmH(8,2)+X(3,7)*H_invDmH(7,2)*H_invDmH(8,1)+&
&X(3,6)*H_invDmH(6,2)*H_invDmH(8,1)+X(3,5)*H_invDmH(5,2)*H_invDmH(&
&8,1)+X(3,4)*H_invDmH(4,2)*H_invDmH(8,1)+H_invDmH(3,2)*X(3,3)*H_in&
&vDmH(8,1)+H_invDmH(2,2)*X(3,2)*H_invDmH(8,1)+H_invDmH(1,2)*X(3,1)&
&*H_invDmH(8,1)
tt266 = 2*H_invDmH(1,3)*H_invDmH(8,3)*tt18*tt22
tt267 = 2*X(3,8)*H_invDmH(8,1)*H_invDmH(8,3)+X(3,7)*H_invDmH(7,1)&
&*H_invDmH(8,3)+X(3,6)*H_invDmH(6,1)*H_invDmH(8,3)+X(3,5)*H_invDmH&
&(5,1)*H_invDmH(8,3)+X(3,4)*H_invDmH(4,1)*H_invDmH(8,3)+H_invDmH(3&
&,1)*X(3,3)*H_invDmH(8,3)+H_invDmH(2,1)*X(3,2)*H_invDmH(8,3)+H_inv&
&DmH(1,1)*X(3,1)*H_invDmH(8,3)+X(3,7)*H_invDmH(7,3)*H_invDmH(8,1)+&
&X(3,6)*H_invDmH(6,3)*H_invDmH(8,1)+X(3,5)*H_invDmH(5,3)*H_invDmH(&
&8,1)+X(3,4)*H_invDmH(4,3)*H_invDmH(8,1)+H_invDmH(3,3)*X(3,3)*H_in&
&vDmH(8,1)+H_invDmH(2,3)*X(3,2)*H_invDmH(8,1)+H_invDmH(1,3)*X(3,1)&
&*H_invDmH(8,1)
tt268 = 2*X(3,8)*H_invDmH(8,2)*H_invDmH(8,3)+X(3,7)*H_invDmH(7,2)&
&*H_invDmH(8,3)+X(3,6)*H_invDmH(6,2)*H_invDmH(8,3)+X(3,5)*H_invDmH&
&(5,2)*H_invDmH(8,3)+X(3,4)*H_invDmH(4,2)*H_invDmH(8,3)+H_invDmH(3&
&,2)*X(3,3)*H_invDmH(8,3)+H_invDmH(2,2)*X(3,2)*H_invDmH(8,3)+H_inv&
&DmH(1,2)*X(3,1)*H_invDmH(8,3)+X(3,7)*H_invDmH(7,3)*H_invDmH(8,2)+&
&X(3,6)*H_invDmH(6,3)*H_invDmH(8,2)+X(3,5)*H_invDmH(5,3)*H_invDmH(&
&8,2)+X(3,4)*H_invDmH(4,3)*H_invDmH(8,2)+H_invDmH(3,3)*X(3,3)*H_in&
&vDmH(8,2)+H_invDmH(2,3)*X(3,2)*H_invDmH(8,2)+H_invDmH(1,3)*X(3,1)&
&*H_invDmH(8,2)
tt269 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt61*tt268+tt60*tt267+tt266&
&+tt54*tt265+tt264+tt263)+lam(1,1)*tt63*tt262)
tt270 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt67*tt72+tt66*tt71+2*tt3*t&
&t20*tt22+tt65*tt70+2*tt2*tt13*tt15+2*tt1*tt6*tt8)+lam(1,1)*tt64*t&
&t69)
tt271 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt81*tt67+tt80*tt66+tt90+tt&
&77*tt65+tt88+tt87)+lam(1,1)*tt74*tt64)
tt272 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt84+2*H_invDmH(1,3)*H_invD&
&mH(2,3)*tt21+tt83+tt82+tt67*tt92+tt66*tt91+tt79+2*H_invDmH(1,2)*H&
&_invDmH(2,2)*tt14+tt78+tt65*tt89+tt76+2*H_invDmH(1,1)*H_invDmH(2,&
&1)*tt7)+tt75+lam(1,1)*tt64*tt86)
tt273 = 2*H_invDmH(1,1)*H_invDmH(2,1)*tt6*tt8
tt274 = 2*H_invDmH(1,2)*H_invDmH(2,2)*tt13*tt15
tt275 = 2*H_invDmH(1,3)*H_invDmH(2,3)*tt20*tt22
tt276 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt67*tt100+tt66*tt99+tt275+&
&tt65*tt97+tt274+tt273)+lam(1,1)*tt64*tt94)
tt277 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt67*tt109+tt66*tt108+tt118&
&+tt65*tt105+tt116+tt115)+lam(1,1)*tt102*tt64)
tt278 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt112+2*H_invDmH(1,3)*H_inv&
&DmH(3,3)*tt21+tt111+tt110+tt67*tt120+tt66*tt119+tt107+2*H_invDmH(&
&1,2)*H_invDmH(3,2)*tt14+tt106+tt65*tt117+tt104+2*H_invDmH(1,1)*H_&
&invDmH(3,1)*tt7)+tt103+lam(1,1)*tt64*tt114)
tt279 = 2*H_invDmH(1,1)*H_invDmH(3,1)*tt6*tt8
tt280 = 2*H_invDmH(1,2)*H_invDmH(3,2)*tt13*tt15
tt281 = 2*H_invDmH(1,3)*H_invDmH(3,3)*tt20*tt22
tt282 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt67*tt128+tt66*tt127+tt281&
&+tt65*tt125+tt280+tt279)+lam(1,1)*tt64*tt122)
tt283 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt67*tt137+tt66*tt136+tt146&
&+tt65*tt133+tt144+tt143)+lam(1,1)*tt130*tt64)
tt284 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt140+2*H_invDmH(1,3)*H_inv&
&DmH(4,3)*tt21+tt139+tt138+tt67*tt148+tt66*tt147+tt135+2*H_invDmH(&
&1,2)*H_invDmH(4,2)*tt14+tt134+tt65*tt145+tt132+2*H_invDmH(1,1)*H_&
&invDmH(4,1)*tt7)+tt131+lam(1,1)*tt64*tt142)
tt285 = 2*H_invDmH(1,1)*H_invDmH(4,1)*tt6*tt8
tt286 = 2*H_invDmH(1,2)*H_invDmH(4,2)*tt13*tt15
tt287 = 2*H_invDmH(1,3)*H_invDmH(4,3)*tt20*tt22
tt288 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt67*tt156+tt66*tt155+tt287&
&+tt65*tt153+tt286+tt285)+lam(1,1)*tt64*tt150)
tt289 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt67*tt165+tt66*tt164+tt174&
&+tt65*tt161+tt172+tt171)+lam(1,1)*tt158*tt64)
tt290 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt168+2*H_invDmH(1,3)*H_inv&
&DmH(5,3)*tt21+tt167+tt166+tt67*tt176+tt66*tt175+tt163+2*H_invDmH(&
&1,2)*H_invDmH(5,2)*tt14+tt162+tt65*tt173+tt160+2*H_invDmH(1,1)*H_&
&invDmH(5,1)*tt7)+tt159+lam(1,1)*tt64*tt170)
tt291 = 2*H_invDmH(1,1)*H_invDmH(5,1)*tt6*tt8
tt292 = 2*H_invDmH(1,2)*H_invDmH(5,2)*tt13*tt15
tt293 = 2*H_invDmH(1,3)*H_invDmH(5,3)*tt20*tt22
tt294 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt67*tt184+tt66*tt183+tt293&
&+tt65*tt181+tt292+tt291)+lam(1,1)*tt64*tt178)
tt295 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt67*tt193+tt66*tt192+tt202&
&+tt65*tt189+tt200+tt199)+lam(1,1)*tt186*tt64)
tt296 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt196+2*H_invDmH(1,3)*H_inv&
&DmH(6,3)*tt21+tt195+tt194+tt67*tt204+tt66*tt203+tt191+2*H_invDmH(&
&1,2)*H_invDmH(6,2)*tt14+tt190+tt65*tt201+tt188+2*H_invDmH(1,1)*H_&
&invDmH(6,1)*tt7)+tt187+lam(1,1)*tt64*tt198)
tt297 = 2*H_invDmH(1,1)*H_invDmH(6,1)*tt6*tt8
tt298 = 2*H_invDmH(1,2)*H_invDmH(6,2)*tt13*tt15
tt299 = 2*H_invDmH(1,3)*H_invDmH(6,3)*tt20*tt22
tt300 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt67*tt212+tt66*tt211+tt299&
&+tt65*tt209+tt298+tt297)+lam(1,1)*tt64*tt206)
tt301 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt67*tt221+tt66*tt220+tt230&
&+tt65*tt217+tt228+tt227)+lam(1,1)*tt214*tt64)
tt302 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt224+2*H_invDmH(1,3)*H_inv&
&DmH(7,3)*tt21+tt223+tt222+tt67*tt232+tt66*tt231+tt219+2*H_invDmH(&
&1,2)*H_invDmH(7,2)*tt14+tt218+tt65*tt229+tt216+2*H_invDmH(1,1)*H_&
&invDmH(7,1)*tt7)+tt215+lam(1,1)*tt64*tt226)
tt303 = 2*H_invDmH(1,1)*H_invDmH(7,1)*tt6*tt8
tt304 = 2*H_invDmH(1,2)*H_invDmH(7,2)*tt13*tt15
tt305 = 2*H_invDmH(1,3)*H_invDmH(7,3)*tt20*tt22
tt306 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt67*tt240+tt66*tt239+tt305&
&+tt65*tt237+tt304+tt303)+lam(1,1)*tt64*tt234)
tt307 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt67*tt250+tt66*tt248+tt258&
&+tt65*tt245+tt256+tt255)+lam(1,1)*tt242*tt64)
tt308 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt252+2*H_invDmH(1,3)*H_inv&
&DmH(8,3)*tt21+tt251+tt67*tt260+tt249+tt66*tt259+tt247+2*H_invDmH(&
&1,2)*H_invDmH(8,2)*tt14+tt246+tt65*tt257+tt244+2*H_invDmH(1,1)*H_&
&invDmH(8,1)*tt7)+tt243+lam(1,1)*tt64*tt254)
tt309 = 2*H_invDmH(1,1)*H_invDmH(8,1)*tt6*tt8
tt310 = 2*H_invDmH(1,2)*H_invDmH(8,2)*tt13*tt15
tt311 = 2*H_invDmH(1,3)*H_invDmH(8,3)*tt20*tt22
tt312 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt67*tt268+tt66*tt267+tt311&
&+tt65*tt265+tt310+tt309)+lam(1,1)*tt64*tt262)
tt313 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt81*tt72+tt80*tt71+tt98+tt&
&77*tt70+tt96+tt95)+lam(1,1)*tt74*tt69)
tt314 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt92*tt72+tt91*tt71+tt275+t&
&t89*tt70+tt274+tt273)+lam(1,1)*tt86*tt69)
tt315 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt84+2*H_invDmH(1,3)*H_invD&
&mH(2,3)*tt23+tt83+tt82+tt72*tt100+tt71*tt99+tt79+2*H_invDmH(1,2)*&
&H_invDmH(2,2)*tt16+tt78+tt70*tt97+tt76+2*H_invDmH(1,1)*H_invDmH(2&
&,1)*tt9)+tt75+lam(1,1)*tt69*tt94)
tt316 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt109*tt72+tt108*tt71+tt126&
&+tt105*tt70+tt124+tt123)+lam(1,1)*tt102*tt69)
tt317 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt120*tt72+tt119*tt71+tt281&
&+tt117*tt70+tt280+tt279)+lam(1,1)*tt114*tt69)
tt318 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt112+2*H_invDmH(1,3)*H_inv&
&DmH(3,3)*tt23+tt111+tt110+tt72*tt128+tt71*tt127+tt107+2*H_invDmH(&
&1,2)*H_invDmH(3,2)*tt16+tt106+tt70*tt125+tt104+2*H_invDmH(1,1)*H_&
&invDmH(3,1)*tt9)+tt103+lam(1,1)*tt69*tt122)
tt319 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt72*tt137+tt71*tt136+tt154&
&+tt70*tt133+tt152+tt151)+lam(1,1)*tt130*tt69)
tt320 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt72*tt148+tt71*tt147+tt287&
&+tt70*tt145+tt286+tt285)+lam(1,1)*tt142*tt69)
tt321 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt140+2*H_invDmH(1,3)*H_inv&
&DmH(4,3)*tt23+tt139+tt138+tt72*tt156+tt71*tt155+tt135+2*H_invDmH(&
&1,2)*H_invDmH(4,2)*tt16+tt134+tt70*tt153+tt132+2*H_invDmH(1,1)*H_&
&invDmH(4,1)*tt9)+tt131+lam(1,1)*tt69*tt150)
tt322 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt72*tt165+tt71*tt164+tt182&
&+tt70*tt161+tt180+tt179)+lam(1,1)*tt158*tt69)
tt323 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt72*tt176+tt71*tt175+tt293&
&+tt70*tt173+tt292+tt291)+lam(1,1)*tt170*tt69)
tt324 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt168+2*H_invDmH(1,3)*H_inv&
&DmH(5,3)*tt23+tt167+tt166+tt72*tt184+tt71*tt183+tt163+2*H_invDmH(&
&1,2)*H_invDmH(5,2)*tt16+tt162+tt70*tt181+tt160+2*H_invDmH(1,1)*H_&
&invDmH(5,1)*tt9)+tt159+lam(1,1)*tt69*tt178)
tt325 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt72*tt193+tt71*tt192+tt210&
&+tt70*tt189+tt208+tt207)+lam(1,1)*tt186*tt69)
tt326 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt72*tt204+tt71*tt203+tt299&
&+tt70*tt201+tt298+tt297)+lam(1,1)*tt198*tt69)
tt327 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt196+2*H_invDmH(1,3)*H_inv&
&DmH(6,3)*tt23+tt195+tt194+tt72*tt212+tt71*tt211+tt191+2*H_invDmH(&
&1,2)*H_invDmH(6,2)*tt16+tt190+tt70*tt209+tt188+2*H_invDmH(1,1)*H_&
&invDmH(6,1)*tt9)+tt187+lam(1,1)*tt69*tt206)
tt328 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt72*tt221+tt71*tt220+tt238&
&+tt70*tt217+tt236+tt235)+lam(1,1)*tt214*tt69)
tt329 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt72*tt232+tt71*tt231+tt305&
&+tt70*tt229+tt304+tt303)+lam(1,1)*tt226*tt69)
tt330 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt224+2*H_invDmH(1,3)*H_inv&
&DmH(7,3)*tt23+tt223+tt222+tt72*tt240+tt71*tt239+tt219+2*H_invDmH(&
&1,2)*H_invDmH(7,2)*tt16+tt218+tt70*tt237+tt216+2*H_invDmH(1,1)*H_&
&invDmH(7,1)*tt9)+tt215+lam(1,1)*tt69*tt234)
tt331 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt72*tt250+tt71*tt248+tt266&
&+tt70*tt245+tt264+tt263)+lam(1,1)*tt242*tt69)
tt332 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt72*tt260+tt71*tt259+tt311&
&+tt70*tt257+tt310+tt309)+lam(1,1)*tt254*tt69)
tt333 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt252+2*H_invDmH(1,3)*H_inv&
&DmH(8,3)*tt23+tt251+tt72*tt268+tt249+tt71*tt267+tt247+2*H_invDmH(&
&1,2)*H_invDmH(8,2)*tt16+tt246+tt70*tt265+tt244+2*H_invDmH(1,1)*H_&
&invDmH(8,1)*tt9)+tt243+lam(1,1)*tt69*tt262)
tt334 = H_invDmH(2,1)**2
tt335 = H_invDmH(2,2)**2
tt336 = H_invDmH(2,3)**2
tt337 = lam(1,1)*(tt336+tt335+tt334)*tt25
tt338 = tt334*tt10
tt339 = 2*H_invDmH(2,1)*H_invDmH(2,2)*tt52
tt340 = tt335*tt17
tt341 = 2*H_invDmH(2,1)*H_invDmH(2,3)*tt56
tt342 = 2*H_invDmH(2,2)*H_invDmH(2,3)*tt58
tt343 = tt336*tt24
tt344 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt81*tt92+tt80*tt91+2*tt336&
&*tt18*tt20+tt77*tt89+2*tt335*tt11*tt13+2*tt334*tt4*tt6)+lam(1,1)*&
&tt74*tt86)
tt345 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt81*tt100+tt80*tt99+2*tt33&
&6*tt18*tt22+tt77*tt97+2*tt335*tt11*tt15+2*tt334*tt4*tt8)+lam(1,1)&
&*tt74*tt94)
tt346 = lam(1,1)*(H_invDmH(2,3)*H_invDmH(3,3)+H_invDmH(2,2)*H_inv&
&DmH(3,2)+H_invDmH(2,1)*H_invDmH(3,1))*tt25
tt347 = H_invDmH(2,1)*H_invDmH(3,1)*tt10
tt348 = (H_invDmH(2,1)*H_invDmH(3,2)+H_invDmH(2,2)*H_invDmH(3,1))&
&*tt52
tt349 = H_invDmH(2,2)*H_invDmH(3,2)*tt17
tt350 = (H_invDmH(2,1)*H_invDmH(3,3)+H_invDmH(2,3)*H_invDmH(3,1))&
&*tt56
tt351 = (H_invDmH(2,2)*H_invDmH(3,3)+H_invDmH(2,3)*H_invDmH(3,2))&
&*tt58
tt352 = H_invDmH(2,3)*H_invDmH(3,3)*tt24
tt353 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt352+2*H_invDmH(2,3)*H_inv&
&DmH(3,3)*tt19+tt351+tt350+tt81*tt109+tt80*tt108+tt349+2*H_invDmH(&
&2,2)*H_invDmH(3,2)*tt12+tt348+tt77*tt105+tt347+2*H_invDmH(2,1)*H_&
&invDmH(3,1)*tt5)+tt346+lam(1,1)*tt74*tt102)
tt354 = 2*H_invDmH(2,1)*H_invDmH(3,1)*tt4*tt6
tt355 = 2*H_invDmH(2,2)*H_invDmH(3,2)*tt11*tt13
tt356 = 2*H_invDmH(2,3)*H_invDmH(3,3)*tt18*tt20
tt357 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt81*tt120+tt80*tt119+tt356&
&+tt77*tt117+tt355+tt354)+lam(1,1)*tt74*tt114)
tt358 = 2*H_invDmH(2,1)*H_invDmH(3,1)*tt4*tt8
tt359 = 2*H_invDmH(2,2)*H_invDmH(3,2)*tt11*tt15
tt360 = 2*H_invDmH(2,3)*H_invDmH(3,3)*tt18*tt22
tt361 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt81*tt128+tt80*tt127+tt360&
&+tt77*tt125+tt359+tt358)+lam(1,1)*tt74*tt122)
tt362 = lam(1,1)*(H_invDmH(2,3)*H_invDmH(4,3)+H_invDmH(2,2)*H_inv&
&DmH(4,2)+H_invDmH(2,1)*H_invDmH(4,1))*tt25
tt363 = H_invDmH(2,1)*H_invDmH(4,1)*tt10
tt364 = (H_invDmH(2,1)*H_invDmH(4,2)+H_invDmH(2,2)*H_invDmH(4,1))&
&*tt52
tt365 = H_invDmH(2,2)*H_invDmH(4,2)*tt17
tt366 = (H_invDmH(2,1)*H_invDmH(4,3)+H_invDmH(2,3)*H_invDmH(4,1))&
&*tt56
tt367 = (H_invDmH(2,2)*H_invDmH(4,3)+H_invDmH(2,3)*H_invDmH(4,2))&
&*tt58
tt368 = H_invDmH(2,3)*H_invDmH(4,3)*tt24
tt369 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt368+2*H_invDmH(2,3)*H_inv&
&DmH(4,3)*tt19+tt367+tt366+tt81*tt137+tt80*tt136+tt365+2*H_invDmH(&
&2,2)*H_invDmH(4,2)*tt12+tt364+tt77*tt133+tt363+2*H_invDmH(2,1)*H_&
&invDmH(4,1)*tt5)+tt362+lam(1,1)*tt74*tt130)
tt370 = 2*H_invDmH(2,1)*H_invDmH(4,1)*tt4*tt6
tt371 = 2*H_invDmH(2,2)*H_invDmH(4,2)*tt11*tt13
tt372 = 2*H_invDmH(2,3)*H_invDmH(4,3)*tt18*tt20
tt373 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt81*tt148+tt80*tt147+tt372&
&+tt77*tt145+tt371+tt370)+lam(1,1)*tt74*tt142)
tt374 = 2*H_invDmH(2,1)*H_invDmH(4,1)*tt4*tt8
tt375 = 2*H_invDmH(2,2)*H_invDmH(4,2)*tt11*tt15
tt376 = 2*H_invDmH(2,3)*H_invDmH(4,3)*tt18*tt22
tt377 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt81*tt156+tt80*tt155+tt376&
&+tt77*tt153+tt375+tt374)+lam(1,1)*tt74*tt150)
tt378 = lam(1,1)*(H_invDmH(2,3)*H_invDmH(5,3)+H_invDmH(2,2)*H_inv&
&DmH(5,2)+H_invDmH(2,1)*H_invDmH(5,1))*tt25
tt379 = H_invDmH(2,1)*H_invDmH(5,1)*tt10
tt380 = (H_invDmH(2,1)*H_invDmH(5,2)+H_invDmH(2,2)*H_invDmH(5,1))&
&*tt52
tt381 = H_invDmH(2,2)*H_invDmH(5,2)*tt17
tt382 = (H_invDmH(2,1)*H_invDmH(5,3)+H_invDmH(2,3)*H_invDmH(5,1))&
&*tt56
tt383 = (H_invDmH(2,2)*H_invDmH(5,3)+H_invDmH(2,3)*H_invDmH(5,2))&
&*tt58
tt384 = H_invDmH(2,3)*H_invDmH(5,3)*tt24
tt385 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt384+2*H_invDmH(2,3)*H_inv&
&DmH(5,3)*tt19+tt383+tt382+tt81*tt165+tt80*tt164+tt381+2*H_invDmH(&
&2,2)*H_invDmH(5,2)*tt12+tt380+tt77*tt161+tt379+2*H_invDmH(2,1)*H_&
&invDmH(5,1)*tt5)+tt378+lam(1,1)*tt74*tt158)
tt386 = 2*H_invDmH(2,1)*H_invDmH(5,1)*tt4*tt6
tt387 = 2*H_invDmH(2,2)*H_invDmH(5,2)*tt11*tt13
tt388 = 2*H_invDmH(2,3)*H_invDmH(5,3)*tt18*tt20
tt389 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt81*tt176+tt80*tt175+tt388&
&+tt77*tt173+tt387+tt386)+lam(1,1)*tt74*tt170)
tt390 = 2*H_invDmH(2,1)*H_invDmH(5,1)*tt4*tt8
tt391 = 2*H_invDmH(2,2)*H_invDmH(5,2)*tt11*tt15
tt392 = 2*H_invDmH(2,3)*H_invDmH(5,3)*tt18*tt22
tt393 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt81*tt184+tt80*tt183+tt392&
&+tt77*tt181+tt391+tt390)+lam(1,1)*tt74*tt178)
tt394 = lam(1,1)*(H_invDmH(2,3)*H_invDmH(6,3)+H_invDmH(2,2)*H_inv&
&DmH(6,2)+H_invDmH(2,1)*H_invDmH(6,1))*tt25
tt395 = H_invDmH(2,1)*H_invDmH(6,1)*tt10
tt396 = (H_invDmH(2,1)*H_invDmH(6,2)+H_invDmH(2,2)*H_invDmH(6,1))&
&*tt52
tt397 = H_invDmH(2,2)*H_invDmH(6,2)*tt17
tt398 = (H_invDmH(2,1)*H_invDmH(6,3)+H_invDmH(2,3)*H_invDmH(6,1))&
&*tt56
tt399 = (H_invDmH(2,2)*H_invDmH(6,3)+H_invDmH(2,3)*H_invDmH(6,2))&
&*tt58
tt400 = H_invDmH(2,3)*H_invDmH(6,3)*tt24
tt401 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt400+2*H_invDmH(2,3)*H_inv&
&DmH(6,3)*tt19+tt399+tt398+tt81*tt193+tt80*tt192+tt397+2*H_invDmH(&
&2,2)*H_invDmH(6,2)*tt12+tt396+tt77*tt189+tt395+2*H_invDmH(2,1)*H_&
&invDmH(6,1)*tt5)+tt394+lam(1,1)*tt74*tt186)
tt402 = 2*H_invDmH(2,1)*H_invDmH(6,1)*tt4*tt6
tt403 = 2*H_invDmH(2,2)*H_invDmH(6,2)*tt11*tt13
tt404 = 2*H_invDmH(2,3)*H_invDmH(6,3)*tt18*tt20
tt405 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt81*tt204+tt80*tt203+tt404&
&+tt77*tt201+tt403+tt402)+lam(1,1)*tt74*tt198)
tt406 = 2*H_invDmH(2,1)*H_invDmH(6,1)*tt4*tt8
tt407 = 2*H_invDmH(2,2)*H_invDmH(6,2)*tt11*tt15
tt408 = 2*H_invDmH(2,3)*H_invDmH(6,3)*tt18*tt22
tt409 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt81*tt212+tt80*tt211+tt408&
&+tt77*tt209+tt407+tt406)+lam(1,1)*tt74*tt206)
tt410 = lam(1,1)*(H_invDmH(2,3)*H_invDmH(7,3)+H_invDmH(2,2)*H_inv&
&DmH(7,2)+H_invDmH(2,1)*H_invDmH(7,1))*tt25
tt411 = H_invDmH(2,1)*H_invDmH(7,1)*tt10
tt412 = (H_invDmH(2,1)*H_invDmH(7,2)+H_invDmH(2,2)*H_invDmH(7,1))&
&*tt52
tt413 = H_invDmH(2,2)*H_invDmH(7,2)*tt17
tt414 = (H_invDmH(2,1)*H_invDmH(7,3)+H_invDmH(2,3)*H_invDmH(7,1))&
&*tt56
tt415 = (H_invDmH(2,2)*H_invDmH(7,3)+H_invDmH(2,3)*H_invDmH(7,2))&
&*tt58
tt416 = H_invDmH(2,3)*H_invDmH(7,3)*tt24
tt417 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt416+2*H_invDmH(2,3)*H_inv&
&DmH(7,3)*tt19+tt415+tt414+tt81*tt221+tt80*tt220+tt413+2*H_invDmH(&
&2,2)*H_invDmH(7,2)*tt12+tt412+tt77*tt217+tt411+2*H_invDmH(2,1)*H_&
&invDmH(7,1)*tt5)+tt410+lam(1,1)*tt74*tt214)
tt418 = 2*H_invDmH(2,1)*H_invDmH(7,1)*tt4*tt6
tt419 = 2*H_invDmH(2,2)*H_invDmH(7,2)*tt11*tt13
tt420 = 2*H_invDmH(2,3)*H_invDmH(7,3)*tt18*tt20
tt421 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt81*tt232+tt80*tt231+tt420&
&+tt77*tt229+tt419+tt418)+lam(1,1)*tt74*tt226)
tt422 = 2*H_invDmH(2,1)*H_invDmH(7,1)*tt4*tt8
tt423 = 2*H_invDmH(2,2)*H_invDmH(7,2)*tt11*tt15
tt424 = 2*H_invDmH(2,3)*H_invDmH(7,3)*tt18*tt22
tt425 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt81*tt240+tt80*tt239+tt424&
&+tt77*tt237+tt423+tt422)+lam(1,1)*tt74*tt234)
tt426 = lam(1,1)*(H_invDmH(2,3)*H_invDmH(8,3)+H_invDmH(2,2)*H_inv&
&DmH(8,2)+H_invDmH(2,1)*H_invDmH(8,1))*tt25
tt427 = H_invDmH(2,1)*H_invDmH(8,1)*tt10
tt428 = (H_invDmH(2,1)*H_invDmH(8,2)+H_invDmH(2,2)*H_invDmH(8,1))&
&*tt52
tt429 = H_invDmH(2,2)*H_invDmH(8,2)*tt17
tt430 = (H_invDmH(2,1)*H_invDmH(8,3)+H_invDmH(2,3)*H_invDmH(8,1))&
&*tt56
tt431 = (H_invDmH(2,2)*H_invDmH(8,3)+H_invDmH(2,3)*H_invDmH(8,2))&
&*tt58
tt432 = H_invDmH(2,3)*H_invDmH(8,3)*tt24
tt433 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt432+2*H_invDmH(2,3)*H_inv&
&DmH(8,3)*tt19+tt431+tt81*tt250+tt430+tt80*tt248+tt429+2*H_invDmH(&
&2,2)*H_invDmH(8,2)*tt12+tt428+tt77*tt245+tt427+2*H_invDmH(2,1)*H_&
&invDmH(8,1)*tt5)+tt426+lam(1,1)*tt74*tt242)
tt434 = 2*H_invDmH(2,1)*H_invDmH(8,1)*tt4*tt6
tt435 = 2*H_invDmH(2,2)*H_invDmH(8,2)*tt11*tt13
tt436 = 2*H_invDmH(2,3)*H_invDmH(8,3)*tt18*tt20
tt437 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt81*tt260+tt80*tt259+tt436&
&+tt77*tt257+tt435+tt434)+lam(1,1)*tt74*tt254)
tt438 = 2*H_invDmH(2,1)*H_invDmH(8,1)*tt4*tt8
tt439 = 2*H_invDmH(2,2)*H_invDmH(8,2)*tt11*tt15
tt440 = 2*H_invDmH(2,3)*H_invDmH(8,3)*tt18*tt22
tt441 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt81*tt268+tt80*tt267+tt440&
&+tt77*tt265+tt439+tt438)+lam(1,1)*tt74*tt262)
tt442 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt92*tt100+tt91*tt99+2*tt33&
&6*tt20*tt22+tt89*tt97+2*tt335*tt13*tt15+2*tt334*tt6*tt8)+lam(1,1)&
&*tt86*tt94)
tt443 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt92*tt109+tt91*tt108+tt356&
&+tt89*tt105+tt355+tt354)+lam(1,1)*tt102*tt86)
tt444 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt352+2*H_invDmH(2,3)*H_inv&
&DmH(3,3)*tt21+tt351+tt350+tt92*tt120+tt91*tt119+tt349+2*H_invDmH(&
&2,2)*H_invDmH(3,2)*tt14+tt348+tt89*tt117+tt347+2*H_invDmH(2,1)*H_&
&invDmH(3,1)*tt7)+tt346+lam(1,1)*tt86*tt114)
tt445 = 2*H_invDmH(2,1)*H_invDmH(3,1)*tt6*tt8
tt446 = 2*H_invDmH(2,2)*H_invDmH(3,2)*tt13*tt15
tt447 = 2*H_invDmH(2,3)*H_invDmH(3,3)*tt20*tt22
tt448 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt92*tt128+tt91*tt127+tt447&
&+tt89*tt125+tt446+tt445)+lam(1,1)*tt86*tt122)
tt449 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt92*tt137+tt91*tt136+tt372&
&+tt89*tt133+tt371+tt370)+lam(1,1)*tt130*tt86)
tt450 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt368+2*H_invDmH(2,3)*H_inv&
&DmH(4,3)*tt21+tt367+tt366+tt92*tt148+tt91*tt147+tt365+2*H_invDmH(&
&2,2)*H_invDmH(4,2)*tt14+tt364+tt89*tt145+tt363+2*H_invDmH(2,1)*H_&
&invDmH(4,1)*tt7)+tt362+lam(1,1)*tt86*tt142)
tt451 = 2*H_invDmH(2,1)*H_invDmH(4,1)*tt6*tt8
tt452 = 2*H_invDmH(2,2)*H_invDmH(4,2)*tt13*tt15
tt453 = 2*H_invDmH(2,3)*H_invDmH(4,3)*tt20*tt22
tt454 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt92*tt156+tt91*tt155+tt453&
&+tt89*tt153+tt452+tt451)+lam(1,1)*tt86*tt150)
tt455 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt92*tt165+tt91*tt164+tt388&
&+tt89*tt161+tt387+tt386)+lam(1,1)*tt158*tt86)
tt456 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt384+2*H_invDmH(2,3)*H_inv&
&DmH(5,3)*tt21+tt383+tt382+tt92*tt176+tt91*tt175+tt381+2*H_invDmH(&
&2,2)*H_invDmH(5,2)*tt14+tt380+tt89*tt173+tt379+2*H_invDmH(2,1)*H_&
&invDmH(5,1)*tt7)+tt378+lam(1,1)*tt86*tt170)
tt457 = 2*H_invDmH(2,1)*H_invDmH(5,1)*tt6*tt8
tt458 = 2*H_invDmH(2,2)*H_invDmH(5,2)*tt13*tt15
tt459 = 2*H_invDmH(2,3)*H_invDmH(5,3)*tt20*tt22
tt460 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt92*tt184+tt91*tt183+tt459&
&+tt89*tt181+tt458+tt457)+lam(1,1)*tt86*tt178)
tt461 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt92*tt193+tt91*tt192+tt404&
&+tt89*tt189+tt403+tt402)+lam(1,1)*tt186*tt86)
tt462 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt400+2*H_invDmH(2,3)*H_inv&
&DmH(6,3)*tt21+tt399+tt398+tt92*tt204+tt91*tt203+tt397+2*H_invDmH(&
&2,2)*H_invDmH(6,2)*tt14+tt396+tt89*tt201+tt395+2*H_invDmH(2,1)*H_&
&invDmH(6,1)*tt7)+tt394+lam(1,1)*tt86*tt198)
tt463 = 2*H_invDmH(2,1)*H_invDmH(6,1)*tt6*tt8
tt464 = 2*H_invDmH(2,2)*H_invDmH(6,2)*tt13*tt15
tt465 = 2*H_invDmH(2,3)*H_invDmH(6,3)*tt20*tt22
tt466 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt92*tt212+tt91*tt211+tt465&
&+tt89*tt209+tt464+tt463)+lam(1,1)*tt86*tt206)
tt467 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt92*tt221+tt91*tt220+tt420&
&+tt89*tt217+tt419+tt418)+lam(1,1)*tt214*tt86)
tt468 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt416+2*H_invDmH(2,3)*H_inv&
&DmH(7,3)*tt21+tt415+tt414+tt92*tt232+tt91*tt231+tt413+2*H_invDmH(&
&2,2)*H_invDmH(7,2)*tt14+tt412+tt89*tt229+tt411+2*H_invDmH(2,1)*H_&
&invDmH(7,1)*tt7)+tt410+lam(1,1)*tt86*tt226)
tt469 = 2*H_invDmH(2,1)*H_invDmH(7,1)*tt6*tt8
tt470 = 2*H_invDmH(2,2)*H_invDmH(7,2)*tt13*tt15
tt471 = 2*H_invDmH(2,3)*H_invDmH(7,3)*tt20*tt22
tt472 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt92*tt240+tt91*tt239+tt471&
&+tt89*tt237+tt470+tt469)+lam(1,1)*tt86*tt234)
tt473 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt92*tt250+tt91*tt248+tt436&
&+tt89*tt245+tt435+tt434)+lam(1,1)*tt242*tt86)
tt474 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt432+2*H_invDmH(2,3)*H_inv&
&DmH(8,3)*tt21+tt431+tt92*tt260+tt430+tt91*tt259+tt429+2*H_invDmH(&
&2,2)*H_invDmH(8,2)*tt14+tt428+tt89*tt257+tt427+2*H_invDmH(2,1)*H_&
&invDmH(8,1)*tt7)+tt426+lam(1,1)*tt86*tt254)
tt475 = 2*H_invDmH(2,1)*H_invDmH(8,1)*tt6*tt8
tt476 = 2*H_invDmH(2,2)*H_invDmH(8,2)*tt13*tt15
tt477 = 2*H_invDmH(2,3)*H_invDmH(8,3)*tt20*tt22
tt478 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt92*tt268+tt91*tt267+tt477&
&+tt89*tt265+tt476+tt475)+lam(1,1)*tt86*tt262)
tt479 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt109*tt100+tt108*tt99+tt36&
&0+tt105*tt97+tt359+tt358)+lam(1,1)*tt102*tt94)
tt480 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt120*tt100+tt119*tt99+tt44&
&7+tt117*tt97+tt446+tt445)+lam(1,1)*tt114*tt94)
tt481 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt352+2*H_invDmH(2,3)*H_inv&
&DmH(3,3)*tt23+tt351+tt350+tt100*tt128+tt99*tt127+tt349+2*H_invDmH&
&(2,2)*H_invDmH(3,2)*tt16+tt348+tt97*tt125+tt347+2*H_invDmH(2,1)*H&
&_invDmH(3,1)*tt9)+tt346+lam(1,1)*tt94*tt122)
tt482 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt100*tt137+tt99*tt136+tt37&
&6+tt97*tt133+tt375+tt374)+lam(1,1)*tt130*tt94)
tt483 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt100*tt148+tt99*tt147+tt45&
&3+tt97*tt145+tt452+tt451)+lam(1,1)*tt142*tt94)
tt484 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt368+2*H_invDmH(2,3)*H_inv&
&DmH(4,3)*tt23+tt367+tt366+tt100*tt156+tt99*tt155+tt365+2*H_invDmH&
&(2,2)*H_invDmH(4,2)*tt16+tt364+tt97*tt153+tt363+2*H_invDmH(2,1)*H&
&_invDmH(4,1)*tt9)+tt362+lam(1,1)*tt94*tt150)
tt485 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt100*tt165+tt99*tt164+tt39&
&2+tt97*tt161+tt391+tt390)+lam(1,1)*tt158*tt94)
tt486 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt100*tt176+tt99*tt175+tt45&
&9+tt97*tt173+tt458+tt457)+lam(1,1)*tt170*tt94)
tt487 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt384+2*H_invDmH(2,3)*H_inv&
&DmH(5,3)*tt23+tt383+tt382+tt100*tt184+tt99*tt183+tt381+2*H_invDmH&
&(2,2)*H_invDmH(5,2)*tt16+tt380+tt97*tt181+tt379+2*H_invDmH(2,1)*H&
&_invDmH(5,1)*tt9)+tt378+lam(1,1)*tt94*tt178)
tt488 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt100*tt193+tt99*tt192+tt40&
&8+tt97*tt189+tt407+tt406)+lam(1,1)*tt186*tt94)
tt489 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt100*tt204+tt99*tt203+tt46&
&5+tt97*tt201+tt464+tt463)+lam(1,1)*tt198*tt94)
tt490 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt400+2*H_invDmH(2,3)*H_inv&
&DmH(6,3)*tt23+tt399+tt398+tt100*tt212+tt99*tt211+tt397+2*H_invDmH&
&(2,2)*H_invDmH(6,2)*tt16+tt396+tt97*tt209+tt395+2*H_invDmH(2,1)*H&
&_invDmH(6,1)*tt9)+tt394+lam(1,1)*tt94*tt206)
tt491 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt100*tt221+tt99*tt220+tt42&
&4+tt97*tt217+tt423+tt422)+lam(1,1)*tt214*tt94)
tt492 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt100*tt232+tt99*tt231+tt47&
&1+tt97*tt229+tt470+tt469)+lam(1,1)*tt226*tt94)
tt493 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt416+2*H_invDmH(2,3)*H_inv&
&DmH(7,3)*tt23+tt415+tt414+tt100*tt240+tt99*tt239+tt413+2*H_invDmH&
&(2,2)*H_invDmH(7,2)*tt16+tt412+tt97*tt237+tt411+2*H_invDmH(2,1)*H&
&_invDmH(7,1)*tt9)+tt410+lam(1,1)*tt94*tt234)
tt494 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt100*tt250+tt99*tt248+tt44&
&0+tt97*tt245+tt439+tt438)+lam(1,1)*tt242*tt94)
tt495 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt100*tt260+tt99*tt259+tt47&
&7+tt97*tt257+tt476+tt475)+lam(1,1)*tt254*tt94)
tt496 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt432+2*H_invDmH(2,3)*H_inv&
&DmH(8,3)*tt23+tt431+tt100*tt268+tt430+tt99*tt267+tt429+2*H_invDmH&
&(2,2)*H_invDmH(8,2)*tt16+tt428+tt97*tt265+tt427+2*H_invDmH(2,1)*H&
&_invDmH(8,1)*tt9)+tt426+lam(1,1)*tt94*tt262)
tt497 = H_invDmH(3,1)**2
tt498 = H_invDmH(3,2)**2
tt499 = H_invDmH(3,3)**2
tt500 = lam(1,1)*(tt499+tt498+tt497)*tt25
tt501 = tt497*tt10
tt502 = 2*H_invDmH(3,1)*H_invDmH(3,2)*tt52
tt503 = tt498*tt17
tt504 = 2*H_invDmH(3,1)*H_invDmH(3,3)*tt56
tt505 = 2*H_invDmH(3,2)*H_invDmH(3,3)*tt58
tt506 = tt499*tt24
tt507 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt109*tt120+tt108*tt119+2*t&
&t499*tt18*tt20+tt105*tt117+2*tt498*tt11*tt13+2*tt497*tt4*tt6)+lam&
&(1,1)*tt102*tt114)
tt508 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt109*tt128+tt108*tt127+2*t&
&t499*tt18*tt22+tt105*tt125+2*tt498*tt11*tt15+2*tt497*tt4*tt8)+lam&
&(1,1)*tt102*tt122)
tt509 = lam(1,1)*(H_invDmH(3,3)*H_invDmH(4,3)+H_invDmH(3,2)*H_inv&
&DmH(4,2)+H_invDmH(3,1)*H_invDmH(4,1))*tt25
tt510 = H_invDmH(3,1)*H_invDmH(4,1)*tt10
tt511 = (H_invDmH(3,1)*H_invDmH(4,2)+H_invDmH(3,2)*H_invDmH(4,1))&
&*tt52
tt512 = H_invDmH(3,2)*H_invDmH(4,2)*tt17
tt513 = (H_invDmH(3,1)*H_invDmH(4,3)+H_invDmH(3,3)*H_invDmH(4,1))&
&*tt56
tt514 = (H_invDmH(3,2)*H_invDmH(4,3)+H_invDmH(3,3)*H_invDmH(4,2))&
&*tt58
tt515 = H_invDmH(3,3)*H_invDmH(4,3)*tt24
tt516 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt515+2*H_invDmH(3,3)*H_inv&
&DmH(4,3)*tt19+tt514+tt513+tt109*tt137+tt108*tt136+tt512+2*H_invDm&
&H(3,2)*H_invDmH(4,2)*tt12+tt511+tt105*tt133+tt510+2*H_invDmH(3,1)&
&*H_invDmH(4,1)*tt5)+tt509+lam(1,1)*tt102*tt130)
tt517 = 2*H_invDmH(3,1)*H_invDmH(4,1)*tt4*tt6
tt518 = 2*H_invDmH(3,2)*H_invDmH(4,2)*tt11*tt13
tt519 = 2*H_invDmH(3,3)*H_invDmH(4,3)*tt18*tt20
tt520 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt109*tt148+tt108*tt147+tt5&
&19+tt105*tt145+tt518+tt517)+lam(1,1)*tt102*tt142)
tt521 = 2*H_invDmH(3,1)*H_invDmH(4,1)*tt4*tt8
tt522 = 2*H_invDmH(3,2)*H_invDmH(4,2)*tt11*tt15
tt523 = 2*H_invDmH(3,3)*H_invDmH(4,3)*tt18*tt22
tt524 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt109*tt156+tt108*tt155+tt5&
&23+tt105*tt153+tt522+tt521)+lam(1,1)*tt102*tt150)
tt525 = lam(1,1)*(H_invDmH(3,3)*H_invDmH(5,3)+H_invDmH(3,2)*H_inv&
&DmH(5,2)+H_invDmH(3,1)*H_invDmH(5,1))*tt25
tt526 = H_invDmH(3,1)*H_invDmH(5,1)*tt10
tt527 = (H_invDmH(3,1)*H_invDmH(5,2)+H_invDmH(3,2)*H_invDmH(5,1))&
&*tt52
tt528 = H_invDmH(3,2)*H_invDmH(5,2)*tt17
tt529 = (H_invDmH(3,1)*H_invDmH(5,3)+H_invDmH(3,3)*H_invDmH(5,1))&
&*tt56
tt530 = (H_invDmH(3,2)*H_invDmH(5,3)+H_invDmH(3,3)*H_invDmH(5,2))&
&*tt58
tt531 = H_invDmH(3,3)*H_invDmH(5,3)*tt24
tt532 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt531+2*H_invDmH(3,3)*H_inv&
&DmH(5,3)*tt19+tt530+tt529+tt109*tt165+tt108*tt164+tt528+2*H_invDm&
&H(3,2)*H_invDmH(5,2)*tt12+tt527+tt105*tt161+tt526+2*H_invDmH(3,1)&
&*H_invDmH(5,1)*tt5)+tt525+lam(1,1)*tt102*tt158)
tt533 = 2*H_invDmH(3,1)*H_invDmH(5,1)*tt4*tt6
tt534 = 2*H_invDmH(3,2)*H_invDmH(5,2)*tt11*tt13
tt535 = 2*H_invDmH(3,3)*H_invDmH(5,3)*tt18*tt20
tt536 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt109*tt176+tt108*tt175+tt5&
&35+tt105*tt173+tt534+tt533)+lam(1,1)*tt102*tt170)
tt537 = 2*H_invDmH(3,1)*H_invDmH(5,1)*tt4*tt8
tt538 = 2*H_invDmH(3,2)*H_invDmH(5,2)*tt11*tt15
tt539 = 2*H_invDmH(3,3)*H_invDmH(5,3)*tt18*tt22
tt540 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt109*tt184+tt108*tt183+tt5&
&39+tt105*tt181+tt538+tt537)+lam(1,1)*tt102*tt178)
tt541 = lam(1,1)*(H_invDmH(3,3)*H_invDmH(6,3)+H_invDmH(3,2)*H_inv&
&DmH(6,2)+H_invDmH(3,1)*H_invDmH(6,1))*tt25
tt542 = H_invDmH(3,1)*H_invDmH(6,1)*tt10
tt543 = (H_invDmH(3,1)*H_invDmH(6,2)+H_invDmH(3,2)*H_invDmH(6,1))&
&*tt52
tt544 = H_invDmH(3,2)*H_invDmH(6,2)*tt17
tt545 = (H_invDmH(3,1)*H_invDmH(6,3)+H_invDmH(3,3)*H_invDmH(6,1))&
&*tt56
tt546 = (H_invDmH(3,2)*H_invDmH(6,3)+H_invDmH(3,3)*H_invDmH(6,2))&
&*tt58
tt547 = H_invDmH(3,3)*H_invDmH(6,3)*tt24
tt548 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt547+2*H_invDmH(3,3)*H_inv&
&DmH(6,3)*tt19+tt546+tt545+tt109*tt193+tt108*tt192+tt544+2*H_invDm&
&H(3,2)*H_invDmH(6,2)*tt12+tt543+tt105*tt189+tt542+2*H_invDmH(3,1)&
&*H_invDmH(6,1)*tt5)+tt541+lam(1,1)*tt102*tt186)
tt549 = 2*H_invDmH(3,1)*H_invDmH(6,1)*tt4*tt6
tt550 = 2*H_invDmH(3,2)*H_invDmH(6,2)*tt11*tt13
tt551 = 2*H_invDmH(3,3)*H_invDmH(6,3)*tt18*tt20
tt552 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt109*tt204+tt108*tt203+tt5&
&51+tt105*tt201+tt550+tt549)+lam(1,1)*tt102*tt198)
tt553 = 2*H_invDmH(3,1)*H_invDmH(6,1)*tt4*tt8
tt554 = 2*H_invDmH(3,2)*H_invDmH(6,2)*tt11*tt15
tt555 = 2*H_invDmH(3,3)*H_invDmH(6,3)*tt18*tt22
tt556 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt109*tt212+tt108*tt211+tt5&
&55+tt105*tt209+tt554+tt553)+lam(1,1)*tt102*tt206)
tt557 = lam(1,1)*(H_invDmH(3,3)*H_invDmH(7,3)+H_invDmH(3,2)*H_inv&
&DmH(7,2)+H_invDmH(3,1)*H_invDmH(7,1))*tt25
tt558 = H_invDmH(3,1)*H_invDmH(7,1)*tt10
tt559 = (H_invDmH(3,1)*H_invDmH(7,2)+H_invDmH(3,2)*H_invDmH(7,1))&
&*tt52
tt560 = H_invDmH(3,2)*H_invDmH(7,2)*tt17
tt561 = (H_invDmH(3,1)*H_invDmH(7,3)+H_invDmH(3,3)*H_invDmH(7,1))&
&*tt56
tt562 = (H_invDmH(3,2)*H_invDmH(7,3)+H_invDmH(3,3)*H_invDmH(7,2))&
&*tt58
tt563 = H_invDmH(3,3)*H_invDmH(7,3)*tt24
tt564 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt563+2*H_invDmH(3,3)*H_inv&
&DmH(7,3)*tt19+tt562+tt561+tt109*tt221+tt108*tt220+tt560+2*H_invDm&
&H(3,2)*H_invDmH(7,2)*tt12+tt559+tt105*tt217+tt558+2*H_invDmH(3,1)&
&*H_invDmH(7,1)*tt5)+tt557+lam(1,1)*tt102*tt214)
tt565 = 2*H_invDmH(3,1)*H_invDmH(7,1)*tt4*tt6
tt566 = 2*H_invDmH(3,2)*H_invDmH(7,2)*tt11*tt13
tt567 = 2*H_invDmH(3,3)*H_invDmH(7,3)*tt18*tt20
tt568 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt109*tt232+tt108*tt231+tt5&
&67+tt105*tt229+tt566+tt565)+lam(1,1)*tt102*tt226)
tt569 = 2*H_invDmH(3,1)*H_invDmH(7,1)*tt4*tt8
tt570 = 2*H_invDmH(3,2)*H_invDmH(7,2)*tt11*tt15
tt571 = 2*H_invDmH(3,3)*H_invDmH(7,3)*tt18*tt22
tt572 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt109*tt240+tt108*tt239+tt5&
&71+tt105*tt237+tt570+tt569)+lam(1,1)*tt102*tt234)
tt573 = lam(1,1)*(H_invDmH(3,3)*H_invDmH(8,3)+H_invDmH(3,2)*H_inv&
&DmH(8,2)+H_invDmH(3,1)*H_invDmH(8,1))*tt25
tt574 = H_invDmH(3,1)*H_invDmH(8,1)*tt10
tt575 = (H_invDmH(3,1)*H_invDmH(8,2)+H_invDmH(3,2)*H_invDmH(8,1))&
&*tt52
tt576 = H_invDmH(3,2)*H_invDmH(8,2)*tt17
tt577 = (H_invDmH(3,1)*H_invDmH(8,3)+H_invDmH(3,3)*H_invDmH(8,1))&
&*tt56
tt578 = (H_invDmH(3,2)*H_invDmH(8,3)+H_invDmH(3,3)*H_invDmH(8,2))&
&*tt58
tt579 = H_invDmH(3,3)*H_invDmH(8,3)*tt24
tt580 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt579+2*H_invDmH(3,3)*H_inv&
&DmH(8,3)*tt19+tt578+tt109*tt250+tt577+tt108*tt248+tt576+2*H_invDm&
&H(3,2)*H_invDmH(8,2)*tt12+tt575+tt105*tt245+tt574+2*H_invDmH(3,1)&
&*H_invDmH(8,1)*tt5)+tt573+lam(1,1)*tt102*tt242)
tt581 = 2*H_invDmH(3,1)*H_invDmH(8,1)*tt4*tt6
tt582 = 2*H_invDmH(3,2)*H_invDmH(8,2)*tt11*tt13
tt583 = 2*H_invDmH(3,3)*H_invDmH(8,3)*tt18*tt20
tt584 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt109*tt260+tt108*tt259+tt5&
&83+tt105*tt257+tt582+tt581)+lam(1,1)*tt102*tt254)
tt585 = 2*H_invDmH(3,1)*H_invDmH(8,1)*tt4*tt8
tt586 = 2*H_invDmH(3,2)*H_invDmH(8,2)*tt11*tt15
tt587 = 2*H_invDmH(3,3)*H_invDmH(8,3)*tt18*tt22
tt588 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt109*tt268+tt108*tt267+tt5&
&87+tt105*tt265+tt586+tt585)+lam(1,1)*tt102*tt262)
tt589 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt120*tt128+tt119*tt127+2*t&
&t499*tt20*tt22+tt117*tt125+2*tt498*tt13*tt15+2*tt497*tt6*tt8)+lam&
&(1,1)*tt114*tt122)
tt590 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt120*tt137+tt119*tt136+tt5&
&19+tt117*tt133+tt518+tt517)+lam(1,1)*tt130*tt114)
tt591 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt515+2*H_invDmH(3,3)*H_inv&
&DmH(4,3)*tt21+tt514+tt513+tt120*tt148+tt119*tt147+tt512+2*H_invDm&
&H(3,2)*H_invDmH(4,2)*tt14+tt511+tt117*tt145+tt510+2*H_invDmH(3,1)&
&*H_invDmH(4,1)*tt7)+tt509+lam(1,1)*tt114*tt142)
tt592 = 2*H_invDmH(3,1)*H_invDmH(4,1)*tt6*tt8
tt593 = 2*H_invDmH(3,2)*H_invDmH(4,2)*tt13*tt15
tt594 = 2*H_invDmH(3,3)*H_invDmH(4,3)*tt20*tt22
tt595 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt120*tt156+tt119*tt155+tt5&
&94+tt117*tt153+tt593+tt592)+lam(1,1)*tt114*tt150)
tt596 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt120*tt165+tt119*tt164+tt5&
&35+tt117*tt161+tt534+tt533)+lam(1,1)*tt158*tt114)
tt597 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt531+2*H_invDmH(3,3)*H_inv&
&DmH(5,3)*tt21+tt530+tt529+tt120*tt176+tt119*tt175+tt528+2*H_invDm&
&H(3,2)*H_invDmH(5,2)*tt14+tt527+tt117*tt173+tt526+2*H_invDmH(3,1)&
&*H_invDmH(5,1)*tt7)+tt525+lam(1,1)*tt114*tt170)
tt598 = 2*H_invDmH(3,1)*H_invDmH(5,1)*tt6*tt8
tt599 = 2*H_invDmH(3,2)*H_invDmH(5,2)*tt13*tt15
tt600 = 2*H_invDmH(3,3)*H_invDmH(5,3)*tt20*tt22
tt601 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt120*tt184+tt119*tt183+tt6&
&00+tt117*tt181+tt599+tt598)+lam(1,1)*tt114*tt178)
tt602 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt120*tt193+tt119*tt192+tt5&
&51+tt117*tt189+tt550+tt549)+lam(1,1)*tt186*tt114)
tt603 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt547+2*H_invDmH(3,3)*H_inv&
&DmH(6,3)*tt21+tt546+tt545+tt120*tt204+tt119*tt203+tt544+2*H_invDm&
&H(3,2)*H_invDmH(6,2)*tt14+tt543+tt117*tt201+tt542+2*H_invDmH(3,1)&
&*H_invDmH(6,1)*tt7)+tt541+lam(1,1)*tt114*tt198)
tt604 = 2*H_invDmH(3,1)*H_invDmH(6,1)*tt6*tt8
tt605 = 2*H_invDmH(3,2)*H_invDmH(6,2)*tt13*tt15
tt606 = 2*H_invDmH(3,3)*H_invDmH(6,3)*tt20*tt22
tt607 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt120*tt212+tt119*tt211+tt6&
&06+tt117*tt209+tt605+tt604)+lam(1,1)*tt114*tt206)
tt608 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt120*tt221+tt119*tt220+tt5&
&67+tt117*tt217+tt566+tt565)+lam(1,1)*tt214*tt114)
tt609 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt563+2*H_invDmH(3,3)*H_inv&
&DmH(7,3)*tt21+tt562+tt561+tt120*tt232+tt119*tt231+tt560+2*H_invDm&
&H(3,2)*H_invDmH(7,2)*tt14+tt559+tt117*tt229+tt558+2*H_invDmH(3,1)&
&*H_invDmH(7,1)*tt7)+tt557+lam(1,1)*tt114*tt226)
tt610 = 2*H_invDmH(3,1)*H_invDmH(7,1)*tt6*tt8
tt611 = 2*H_invDmH(3,2)*H_invDmH(7,2)*tt13*tt15
tt612 = 2*H_invDmH(3,3)*H_invDmH(7,3)*tt20*tt22
tt613 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt120*tt240+tt119*tt239+tt6&
&12+tt117*tt237+tt611+tt610)+lam(1,1)*tt114*tt234)
tt614 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt120*tt250+tt119*tt248+tt5&
&83+tt117*tt245+tt582+tt581)+lam(1,1)*tt242*tt114)
tt615 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt579+2*H_invDmH(3,3)*H_inv&
&DmH(8,3)*tt21+tt578+tt120*tt260+tt577+tt119*tt259+tt576+2*H_invDm&
&H(3,2)*H_invDmH(8,2)*tt14+tt575+tt117*tt257+tt574+2*H_invDmH(3,1)&
&*H_invDmH(8,1)*tt7)+tt573+lam(1,1)*tt114*tt254)
tt616 = 2*H_invDmH(3,1)*H_invDmH(8,1)*tt6*tt8
tt617 = 2*H_invDmH(3,2)*H_invDmH(8,2)*tt13*tt15
tt618 = 2*H_invDmH(3,3)*H_invDmH(8,3)*tt20*tt22
tt619 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt120*tt268+tt119*tt267+tt6&
&18+tt117*tt265+tt617+tt616)+lam(1,1)*tt114*tt262)
tt620 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt128*tt137+tt127*tt136+tt5&
&23+tt125*tt133+tt522+tt521)+lam(1,1)*tt130*tt122)
tt621 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt128*tt148+tt127*tt147+tt5&
&94+tt125*tt145+tt593+tt592)+lam(1,1)*tt142*tt122)
tt622 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt515+2*H_invDmH(3,3)*H_inv&
&DmH(4,3)*tt23+tt514+tt513+tt128*tt156+tt127*tt155+tt512+2*H_invDm&
&H(3,2)*H_invDmH(4,2)*tt16+tt511+tt125*tt153+tt510+2*H_invDmH(3,1)&
&*H_invDmH(4,1)*tt9)+tt509+lam(1,1)*tt122*tt150)
tt623 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt128*tt165+tt127*tt164+tt5&
&39+tt125*tt161+tt538+tt537)+lam(1,1)*tt158*tt122)
tt624 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt128*tt176+tt127*tt175+tt6&
&00+tt125*tt173+tt599+tt598)+lam(1,1)*tt170*tt122)
tt625 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt531+2*H_invDmH(3,3)*H_inv&
&DmH(5,3)*tt23+tt530+tt529+tt128*tt184+tt127*tt183+tt528+2*H_invDm&
&H(3,2)*H_invDmH(5,2)*tt16+tt527+tt125*tt181+tt526+2*H_invDmH(3,1)&
&*H_invDmH(5,1)*tt9)+tt525+lam(1,1)*tt122*tt178)
tt626 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt128*tt193+tt127*tt192+tt5&
&55+tt125*tt189+tt554+tt553)+lam(1,1)*tt186*tt122)
tt627 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt128*tt204+tt127*tt203+tt6&
&06+tt125*tt201+tt605+tt604)+lam(1,1)*tt198*tt122)
tt628 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt547+2*H_invDmH(3,3)*H_inv&
&DmH(6,3)*tt23+tt546+tt545+tt128*tt212+tt127*tt211+tt544+2*H_invDm&
&H(3,2)*H_invDmH(6,2)*tt16+tt543+tt125*tt209+tt542+2*H_invDmH(3,1)&
&*H_invDmH(6,1)*tt9)+tt541+lam(1,1)*tt122*tt206)
tt629 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt128*tt221+tt127*tt220+tt5&
&71+tt125*tt217+tt570+tt569)+lam(1,1)*tt214*tt122)
tt630 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt128*tt232+tt127*tt231+tt6&
&12+tt125*tt229+tt611+tt610)+lam(1,1)*tt226*tt122)
tt631 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt563+2*H_invDmH(3,3)*H_inv&
&DmH(7,3)*tt23+tt562+tt561+tt128*tt240+tt127*tt239+tt560+2*H_invDm&
&H(3,2)*H_invDmH(7,2)*tt16+tt559+tt125*tt237+tt558+2*H_invDmH(3,1)&
&*H_invDmH(7,1)*tt9)+tt557+lam(1,1)*tt122*tt234)
tt632 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt128*tt250+tt127*tt248+tt5&
&87+tt125*tt245+tt586+tt585)+lam(1,1)*tt242*tt122)
tt633 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt128*tt260+tt127*tt259+tt6&
&18+tt125*tt257+tt617+tt616)+lam(1,1)*tt254*tt122)
tt634 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt579+2*H_invDmH(3,3)*H_inv&
&DmH(8,3)*tt23+tt578+tt128*tt268+tt577+tt127*tt267+tt576+2*H_invDm&
&H(3,2)*H_invDmH(8,2)*tt16+tt575+tt125*tt265+tt574+2*H_invDmH(3,1)&
&*H_invDmH(8,1)*tt9)+tt573+lam(1,1)*tt122*tt262)
tt635 = H_invDmH(4,1)**2
tt636 = H_invDmH(4,2)**2
tt637 = H_invDmH(4,3)**2
tt638 = lam(1,1)*(tt637+tt636+tt635)*tt25
tt639 = tt635*tt10
tt640 = 2*H_invDmH(4,1)*H_invDmH(4,2)*tt52
tt641 = tt636*tt17
tt642 = 2*H_invDmH(4,1)*H_invDmH(4,3)*tt56
tt643 = 2*H_invDmH(4,2)*H_invDmH(4,3)*tt58
tt644 = tt637*tt24
tt645 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt137*tt148+tt136*tt147+2*t&
&t637*tt18*tt20+tt133*tt145+2*tt636*tt11*tt13+2*tt635*tt4*tt6)+lam&
&(1,1)*tt130*tt142)
tt646 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt137*tt156+tt136*tt155+2*t&
&t637*tt18*tt22+tt133*tt153+2*tt636*tt11*tt15+2*tt635*tt4*tt8)+lam&
&(1,1)*tt130*tt150)
tt647 = lam(1,1)*(H_invDmH(4,3)*H_invDmH(5,3)+H_invDmH(4,2)*H_inv&
&DmH(5,2)+H_invDmH(4,1)*H_invDmH(5,1))*tt25
tt648 = H_invDmH(4,1)*H_invDmH(5,1)*tt10
tt649 = (H_invDmH(4,1)*H_invDmH(5,2)+H_invDmH(4,2)*H_invDmH(5,1))&
&*tt52
tt650 = H_invDmH(4,2)*H_invDmH(5,2)*tt17
tt651 = (H_invDmH(4,1)*H_invDmH(5,3)+H_invDmH(4,3)*H_invDmH(5,1))&
&*tt56
tt652 = (H_invDmH(4,2)*H_invDmH(5,3)+H_invDmH(4,3)*H_invDmH(5,2))&
&*tt58
tt653 = H_invDmH(4,3)*H_invDmH(5,3)*tt24
tt654 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt653+2*H_invDmH(4,3)*H_inv&
&DmH(5,3)*tt19+tt652+tt651+tt137*tt165+tt136*tt164+tt650+2*H_invDm&
&H(4,2)*H_invDmH(5,2)*tt12+tt649+tt133*tt161+tt648+2*H_invDmH(4,1)&
&*H_invDmH(5,1)*tt5)+tt647+lam(1,1)*tt130*tt158)
tt655 = 2*H_invDmH(4,1)*H_invDmH(5,1)*tt4*tt6
tt656 = 2*H_invDmH(4,2)*H_invDmH(5,2)*tt11*tt13
tt657 = 2*H_invDmH(4,3)*H_invDmH(5,3)*tt18*tt20
tt658 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt137*tt176+tt136*tt175+tt6&
&57+tt133*tt173+tt656+tt655)+lam(1,1)*tt130*tt170)
tt659 = 2*H_invDmH(4,1)*H_invDmH(5,1)*tt4*tt8
tt660 = 2*H_invDmH(4,2)*H_invDmH(5,2)*tt11*tt15
tt661 = 2*H_invDmH(4,3)*H_invDmH(5,3)*tt18*tt22
tt662 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt137*tt184+tt136*tt183+tt6&
&61+tt133*tt181+tt660+tt659)+lam(1,1)*tt130*tt178)
tt663 = lam(1,1)*(H_invDmH(4,3)*H_invDmH(6,3)+H_invDmH(4,2)*H_inv&
&DmH(6,2)+H_invDmH(4,1)*H_invDmH(6,1))*tt25
tt664 = H_invDmH(4,1)*H_invDmH(6,1)*tt10
tt665 = (H_invDmH(4,1)*H_invDmH(6,2)+H_invDmH(4,2)*H_invDmH(6,1))&
&*tt52
tt666 = H_invDmH(4,2)*H_invDmH(6,2)*tt17
tt667 = (H_invDmH(4,1)*H_invDmH(6,3)+H_invDmH(4,3)*H_invDmH(6,1))&
&*tt56
tt668 = (H_invDmH(4,2)*H_invDmH(6,3)+H_invDmH(4,3)*H_invDmH(6,2))&
&*tt58
tt669 = H_invDmH(4,3)*H_invDmH(6,3)*tt24
tt670 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt669+2*H_invDmH(4,3)*H_inv&
&DmH(6,3)*tt19+tt668+tt667+tt137*tt193+tt136*tt192+tt666+2*H_invDm&
&H(4,2)*H_invDmH(6,2)*tt12+tt665+tt133*tt189+tt664+2*H_invDmH(4,1)&
&*H_invDmH(6,1)*tt5)+tt663+lam(1,1)*tt130*tt186)
tt671 = 2*H_invDmH(4,1)*H_invDmH(6,1)*tt4*tt6
tt672 = 2*H_invDmH(4,2)*H_invDmH(6,2)*tt11*tt13
tt673 = 2*H_invDmH(4,3)*H_invDmH(6,3)*tt18*tt20
tt674 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt137*tt204+tt136*tt203+tt6&
&73+tt133*tt201+tt672+tt671)+lam(1,1)*tt130*tt198)
tt675 = 2*H_invDmH(4,1)*H_invDmH(6,1)*tt4*tt8
tt676 = 2*H_invDmH(4,2)*H_invDmH(6,2)*tt11*tt15
tt677 = 2*H_invDmH(4,3)*H_invDmH(6,3)*tt18*tt22
tt678 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt137*tt212+tt136*tt211+tt6&
&77+tt133*tt209+tt676+tt675)+lam(1,1)*tt130*tt206)
tt679 = lam(1,1)*(H_invDmH(4,3)*H_invDmH(7,3)+H_invDmH(4,2)*H_inv&
&DmH(7,2)+H_invDmH(4,1)*H_invDmH(7,1))*tt25
tt680 = H_invDmH(4,1)*H_invDmH(7,1)*tt10
tt681 = (H_invDmH(4,1)*H_invDmH(7,2)+H_invDmH(4,2)*H_invDmH(7,1))&
&*tt52
tt682 = H_invDmH(4,2)*H_invDmH(7,2)*tt17
tt683 = (H_invDmH(4,1)*H_invDmH(7,3)+H_invDmH(4,3)*H_invDmH(7,1))&
&*tt56
tt684 = (H_invDmH(4,2)*H_invDmH(7,3)+H_invDmH(4,3)*H_invDmH(7,2))&
&*tt58
tt685 = H_invDmH(4,3)*H_invDmH(7,3)*tt24
tt686 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt685+2*H_invDmH(4,3)*H_inv&
&DmH(7,3)*tt19+tt684+tt683+tt137*tt221+tt136*tt220+tt682+2*H_invDm&
&H(4,2)*H_invDmH(7,2)*tt12+tt681+tt133*tt217+tt680+2*H_invDmH(4,1)&
&*H_invDmH(7,1)*tt5)+tt679+lam(1,1)*tt130*tt214)
tt687 = 2*H_invDmH(4,1)*H_invDmH(7,1)*tt4*tt6
tt688 = 2*H_invDmH(4,2)*H_invDmH(7,2)*tt11*tt13
tt689 = 2*H_invDmH(4,3)*H_invDmH(7,3)*tt18*tt20
tt690 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt137*tt232+tt136*tt231+tt6&
&89+tt133*tt229+tt688+tt687)+lam(1,1)*tt130*tt226)
tt691 = 2*H_invDmH(4,1)*H_invDmH(7,1)*tt4*tt8
tt692 = 2*H_invDmH(4,2)*H_invDmH(7,2)*tt11*tt15
tt693 = 2*H_invDmH(4,3)*H_invDmH(7,3)*tt18*tt22
tt694 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt137*tt240+tt136*tt239+tt6&
&93+tt133*tt237+tt692+tt691)+lam(1,1)*tt130*tt234)
tt695 = lam(1,1)*(H_invDmH(4,3)*H_invDmH(8,3)+H_invDmH(4,2)*H_inv&
&DmH(8,2)+H_invDmH(4,1)*H_invDmH(8,1))*tt25
tt696 = H_invDmH(4,1)*H_invDmH(8,1)*tt10
tt697 = (H_invDmH(4,1)*H_invDmH(8,2)+H_invDmH(4,2)*H_invDmH(8,1))&
&*tt52
tt698 = H_invDmH(4,2)*H_invDmH(8,2)*tt17
tt699 = (H_invDmH(4,1)*H_invDmH(8,3)+H_invDmH(4,3)*H_invDmH(8,1))&
&*tt56
tt700 = (H_invDmH(4,2)*H_invDmH(8,3)+H_invDmH(4,3)*H_invDmH(8,2))&
&*tt58
tt701 = H_invDmH(4,3)*H_invDmH(8,3)*tt24
tt702 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt701+2*H_invDmH(4,3)*H_inv&
&DmH(8,3)*tt19+tt700+tt137*tt250+tt699+tt136*tt248+tt698+2*H_invDm&
&H(4,2)*H_invDmH(8,2)*tt12+tt697+tt133*tt245+tt696+2*H_invDmH(4,1)&
&*H_invDmH(8,1)*tt5)+tt695+lam(1,1)*tt130*tt242)
tt703 = 2*H_invDmH(4,1)*H_invDmH(8,1)*tt4*tt6
tt704 = 2*H_invDmH(4,2)*H_invDmH(8,2)*tt11*tt13
tt705 = 2*H_invDmH(4,3)*H_invDmH(8,3)*tt18*tt20
tt706 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt137*tt260+tt136*tt259+tt7&
&05+tt133*tt257+tt704+tt703)+lam(1,1)*tt130*tt254)
tt707 = 2*H_invDmH(4,1)*H_invDmH(8,1)*tt4*tt8
tt708 = 2*H_invDmH(4,2)*H_invDmH(8,2)*tt11*tt15
tt709 = 2*H_invDmH(4,3)*H_invDmH(8,3)*tt18*tt22
tt710 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt137*tt268+tt136*tt267+tt7&
&09+tt133*tt265+tt708+tt707)+lam(1,1)*tt130*tt262)
tt711 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt148*tt156+tt147*tt155+2*t&
&t637*tt20*tt22+tt145*tt153+2*tt636*tt13*tt15+2*tt635*tt6*tt8)+lam&
&(1,1)*tt142*tt150)
tt712 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt148*tt165+tt147*tt164+tt6&
&57+tt145*tt161+tt656+tt655)+lam(1,1)*tt158*tt142)
tt713 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt653+2*H_invDmH(4,3)*H_inv&
&DmH(5,3)*tt21+tt652+tt651+tt148*tt176+tt147*tt175+tt650+2*H_invDm&
&H(4,2)*H_invDmH(5,2)*tt14+tt649+tt145*tt173+tt648+2*H_invDmH(4,1)&
&*H_invDmH(5,1)*tt7)+tt647+lam(1,1)*tt142*tt170)
tt714 = 2*H_invDmH(4,1)*H_invDmH(5,1)*tt6*tt8
tt715 = 2*H_invDmH(4,2)*H_invDmH(5,2)*tt13*tt15
tt716 = 2*H_invDmH(4,3)*H_invDmH(5,3)*tt20*tt22
tt717 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt148*tt184+tt147*tt183+tt7&
&16+tt145*tt181+tt715+tt714)+lam(1,1)*tt142*tt178)
tt718 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt148*tt193+tt147*tt192+tt6&
&73+tt145*tt189+tt672+tt671)+lam(1,1)*tt186*tt142)
tt719 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt669+2*H_invDmH(4,3)*H_inv&
&DmH(6,3)*tt21+tt668+tt667+tt148*tt204+tt147*tt203+tt666+2*H_invDm&
&H(4,2)*H_invDmH(6,2)*tt14+tt665+tt145*tt201+tt664+2*H_invDmH(4,1)&
&*H_invDmH(6,1)*tt7)+tt663+lam(1,1)*tt142*tt198)
tt720 = 2*H_invDmH(4,1)*H_invDmH(6,1)*tt6*tt8
tt721 = 2*H_invDmH(4,2)*H_invDmH(6,2)*tt13*tt15
tt722 = 2*H_invDmH(4,3)*H_invDmH(6,3)*tt20*tt22
tt723 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt148*tt212+tt147*tt211+tt7&
&22+tt145*tt209+tt721+tt720)+lam(1,1)*tt142*tt206)
tt724 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt148*tt221+tt147*tt220+tt6&
&89+tt145*tt217+tt688+tt687)+lam(1,1)*tt214*tt142)
tt725 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt685+2*H_invDmH(4,3)*H_inv&
&DmH(7,3)*tt21+tt684+tt683+tt148*tt232+tt147*tt231+tt682+2*H_invDm&
&H(4,2)*H_invDmH(7,2)*tt14+tt681+tt145*tt229+tt680+2*H_invDmH(4,1)&
&*H_invDmH(7,1)*tt7)+tt679+lam(1,1)*tt142*tt226)
tt726 = 2*H_invDmH(4,1)*H_invDmH(7,1)*tt6*tt8
tt727 = 2*H_invDmH(4,2)*H_invDmH(7,2)*tt13*tt15
tt728 = 2*H_invDmH(4,3)*H_invDmH(7,3)*tt20*tt22
tt729 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt148*tt240+tt147*tt239+tt7&
&28+tt145*tt237+tt727+tt726)+lam(1,1)*tt142*tt234)
tt730 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt148*tt250+tt147*tt248+tt7&
&05+tt145*tt245+tt704+tt703)+lam(1,1)*tt242*tt142)
tt731 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt701+2*H_invDmH(4,3)*H_inv&
&DmH(8,3)*tt21+tt700+tt148*tt260+tt699+tt147*tt259+tt698+2*H_invDm&
&H(4,2)*H_invDmH(8,2)*tt14+tt697+tt145*tt257+tt696+2*H_invDmH(4,1)&
&*H_invDmH(8,1)*tt7)+tt695+lam(1,1)*tt142*tt254)
tt732 = 2*H_invDmH(4,1)*H_invDmH(8,1)*tt6*tt8
tt733 = 2*H_invDmH(4,2)*H_invDmH(8,2)*tt13*tt15
tt734 = 2*H_invDmH(4,3)*H_invDmH(8,3)*tt20*tt22
tt735 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt148*tt268+tt147*tt267+tt7&
&34+tt145*tt265+tt733+tt732)+lam(1,1)*tt142*tt262)
tt736 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt156*tt165+tt155*tt164+tt6&
&61+tt153*tt161+tt660+tt659)+lam(1,1)*tt158*tt150)
tt737 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt156*tt176+tt155*tt175+tt7&
&16+tt153*tt173+tt715+tt714)+lam(1,1)*tt170*tt150)
tt738 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt653+2*H_invDmH(4,3)*H_inv&
&DmH(5,3)*tt23+tt652+tt651+tt156*tt184+tt155*tt183+tt650+2*H_invDm&
&H(4,2)*H_invDmH(5,2)*tt16+tt649+tt153*tt181+tt648+2*H_invDmH(4,1)&
&*H_invDmH(5,1)*tt9)+tt647+lam(1,1)*tt150*tt178)
tt739 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt156*tt193+tt155*tt192+tt6&
&77+tt153*tt189+tt676+tt675)+lam(1,1)*tt186*tt150)
tt740 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt156*tt204+tt155*tt203+tt7&
&22+tt153*tt201+tt721+tt720)+lam(1,1)*tt198*tt150)
tt741 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt669+2*H_invDmH(4,3)*H_inv&
&DmH(6,3)*tt23+tt668+tt667+tt156*tt212+tt155*tt211+tt666+2*H_invDm&
&H(4,2)*H_invDmH(6,2)*tt16+tt665+tt153*tt209+tt664+2*H_invDmH(4,1)&
&*H_invDmH(6,1)*tt9)+tt663+lam(1,1)*tt150*tt206)
tt742 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt156*tt221+tt155*tt220+tt6&
&93+tt153*tt217+tt692+tt691)+lam(1,1)*tt214*tt150)
tt743 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt156*tt232+tt155*tt231+tt7&
&28+tt153*tt229+tt727+tt726)+lam(1,1)*tt226*tt150)
tt744 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt685+2*H_invDmH(4,3)*H_inv&
&DmH(7,3)*tt23+tt684+tt683+tt156*tt240+tt155*tt239+tt682+2*H_invDm&
&H(4,2)*H_invDmH(7,2)*tt16+tt681+tt153*tt237+tt680+2*H_invDmH(4,1)&
&*H_invDmH(7,1)*tt9)+tt679+lam(1,1)*tt150*tt234)
tt745 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt156*tt250+tt155*tt248+tt7&
&09+tt153*tt245+tt708+tt707)+lam(1,1)*tt242*tt150)
tt746 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt156*tt260+tt155*tt259+tt7&
&34+tt153*tt257+tt733+tt732)+lam(1,1)*tt254*tt150)
tt747 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt701+2*H_invDmH(4,3)*H_inv&
&DmH(8,3)*tt23+tt700+tt156*tt268+tt699+tt155*tt267+tt698+2*H_invDm&
&H(4,2)*H_invDmH(8,2)*tt16+tt697+tt153*tt265+tt696+2*H_invDmH(4,1)&
&*H_invDmH(8,1)*tt9)+tt695+lam(1,1)*tt150*tt262)
tt748 = H_invDmH(5,1)**2
tt749 = H_invDmH(5,2)**2
tt750 = H_invDmH(5,3)**2
tt751 = lam(1,1)*(tt750+tt749+tt748)*tt25
tt752 = tt748*tt10
tt753 = 2*H_invDmH(5,1)*H_invDmH(5,2)*tt52
tt754 = tt749*tt17
tt755 = 2*H_invDmH(5,1)*H_invDmH(5,3)*tt56
tt756 = 2*H_invDmH(5,2)*H_invDmH(5,3)*tt58
tt757 = tt750*tt24
tt758 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt165*tt176+tt164*tt175+2*t&
&t750*tt18*tt20+tt161*tt173+2*tt749*tt11*tt13+2*tt748*tt4*tt6)+lam&
&(1,1)*tt158*tt170)
tt759 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt165*tt184+tt164*tt183+2*t&
&t750*tt18*tt22+tt161*tt181+2*tt749*tt11*tt15+2*tt748*tt4*tt8)+lam&
&(1,1)*tt158*tt178)
tt760 = lam(1,1)*(H_invDmH(5,3)*H_invDmH(6,3)+H_invDmH(5,2)*H_inv&
&DmH(6,2)+H_invDmH(5,1)*H_invDmH(6,1))*tt25
tt761 = H_invDmH(5,1)*H_invDmH(6,1)*tt10
tt762 = (H_invDmH(5,1)*H_invDmH(6,2)+H_invDmH(5,2)*H_invDmH(6,1))&
&*tt52
tt763 = H_invDmH(5,2)*H_invDmH(6,2)*tt17
tt764 = (H_invDmH(5,1)*H_invDmH(6,3)+H_invDmH(5,3)*H_invDmH(6,1))&
&*tt56
tt765 = (H_invDmH(5,2)*H_invDmH(6,3)+H_invDmH(5,3)*H_invDmH(6,2))&
&*tt58
tt766 = H_invDmH(5,3)*H_invDmH(6,3)*tt24
tt767 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt766+2*H_invDmH(5,3)*H_inv&
&DmH(6,3)*tt19+tt765+tt764+tt165*tt193+tt164*tt192+tt763+2*H_invDm&
&H(5,2)*H_invDmH(6,2)*tt12+tt762+tt161*tt189+tt761+2*H_invDmH(5,1)&
&*H_invDmH(6,1)*tt5)+tt760+lam(1,1)*tt158*tt186)
tt768 = 2*H_invDmH(5,1)*H_invDmH(6,1)*tt4*tt6
tt769 = 2*H_invDmH(5,2)*H_invDmH(6,2)*tt11*tt13
tt770 = 2*H_invDmH(5,3)*H_invDmH(6,3)*tt18*tt20
tt771 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt165*tt204+tt164*tt203+tt7&
&70+tt161*tt201+tt769+tt768)+lam(1,1)*tt158*tt198)
tt772 = 2*H_invDmH(5,1)*H_invDmH(6,1)*tt4*tt8
tt773 = 2*H_invDmH(5,2)*H_invDmH(6,2)*tt11*tt15
tt774 = 2*H_invDmH(5,3)*H_invDmH(6,3)*tt18*tt22
tt775 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt165*tt212+tt164*tt211+tt7&
&74+tt161*tt209+tt773+tt772)+lam(1,1)*tt158*tt206)
tt776 = lam(1,1)*(H_invDmH(5,3)*H_invDmH(7,3)+H_invDmH(5,2)*H_inv&
&DmH(7,2)+H_invDmH(5,1)*H_invDmH(7,1))*tt25
tt777 = H_invDmH(5,1)*H_invDmH(7,1)*tt10
tt778 = (H_invDmH(5,1)*H_invDmH(7,2)+H_invDmH(5,2)*H_invDmH(7,1))&
&*tt52
tt779 = H_invDmH(5,2)*H_invDmH(7,2)*tt17
tt780 = (H_invDmH(5,1)*H_invDmH(7,3)+H_invDmH(5,3)*H_invDmH(7,1))&
&*tt56
tt781 = (H_invDmH(5,2)*H_invDmH(7,3)+H_invDmH(5,3)*H_invDmH(7,2))&
&*tt58
tt782 = H_invDmH(5,3)*H_invDmH(7,3)*tt24
tt783 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt782+2*H_invDmH(5,3)*H_inv&
&DmH(7,3)*tt19+tt781+tt780+tt165*tt221+tt164*tt220+tt779+2*H_invDm&
&H(5,2)*H_invDmH(7,2)*tt12+tt778+tt161*tt217+tt777+2*H_invDmH(5,1)&
&*H_invDmH(7,1)*tt5)+tt776+lam(1,1)*tt158*tt214)
tt784 = 2*H_invDmH(5,1)*H_invDmH(7,1)*tt4*tt6
tt785 = 2*H_invDmH(5,2)*H_invDmH(7,2)*tt11*tt13
tt786 = 2*H_invDmH(5,3)*H_invDmH(7,3)*tt18*tt20
tt787 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt165*tt232+tt164*tt231+tt7&
&86+tt161*tt229+tt785+tt784)+lam(1,1)*tt158*tt226)
tt788 = 2*H_invDmH(5,1)*H_invDmH(7,1)*tt4*tt8
tt789 = 2*H_invDmH(5,2)*H_invDmH(7,2)*tt11*tt15
tt790 = 2*H_invDmH(5,3)*H_invDmH(7,3)*tt18*tt22
tt791 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt165*tt240+tt164*tt239+tt7&
&90+tt161*tt237+tt789+tt788)+lam(1,1)*tt158*tt234)
tt792 = lam(1,1)*(H_invDmH(5,3)*H_invDmH(8,3)+H_invDmH(5,2)*H_inv&
&DmH(8,2)+H_invDmH(5,1)*H_invDmH(8,1))*tt25
tt793 = H_invDmH(5,1)*H_invDmH(8,1)*tt10
tt794 = (H_invDmH(5,1)*H_invDmH(8,2)+H_invDmH(5,2)*H_invDmH(8,1))&
&*tt52
tt795 = H_invDmH(5,2)*H_invDmH(8,2)*tt17
tt796 = (H_invDmH(5,1)*H_invDmH(8,3)+H_invDmH(5,3)*H_invDmH(8,1))&
&*tt56
tt797 = (H_invDmH(5,2)*H_invDmH(8,3)+H_invDmH(5,3)*H_invDmH(8,2))&
&*tt58
tt798 = H_invDmH(5,3)*H_invDmH(8,3)*tt24
tt799 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt798+2*H_invDmH(5,3)*H_inv&
&DmH(8,3)*tt19+tt797+tt165*tt250+tt796+tt164*tt248+tt795+2*H_invDm&
&H(5,2)*H_invDmH(8,2)*tt12+tt794+tt161*tt245+tt793+2*H_invDmH(5,1)&
&*H_invDmH(8,1)*tt5)+tt792+lam(1,1)*tt158*tt242)
tt800 = 2*H_invDmH(5,1)*H_invDmH(8,1)*tt4*tt6
tt801 = 2*H_invDmH(5,2)*H_invDmH(8,2)*tt11*tt13
tt802 = 2*H_invDmH(5,3)*H_invDmH(8,3)*tt18*tt20
tt803 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt165*tt260+tt164*tt259+tt8&
&02+tt161*tt257+tt801+tt800)+lam(1,1)*tt158*tt254)
tt804 = 2*H_invDmH(5,1)*H_invDmH(8,1)*tt4*tt8
tt805 = 2*H_invDmH(5,2)*H_invDmH(8,2)*tt11*tt15
tt806 = 2*H_invDmH(5,3)*H_invDmH(8,3)*tt18*tt22
tt807 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt165*tt268+tt164*tt267+tt8&
&06+tt161*tt265+tt805+tt804)+lam(1,1)*tt158*tt262)
tt808 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt176*tt184+tt175*tt183+2*t&
&t750*tt20*tt22+tt173*tt181+2*tt749*tt13*tt15+2*tt748*tt6*tt8)+lam&
&(1,1)*tt170*tt178)
tt809 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt176*tt193+tt175*tt192+tt7&
&70+tt173*tt189+tt769+tt768)+lam(1,1)*tt186*tt170)
tt810 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt766+2*H_invDmH(5,3)*H_inv&
&DmH(6,3)*tt21+tt765+tt764+tt176*tt204+tt175*tt203+tt763+2*H_invDm&
&H(5,2)*H_invDmH(6,2)*tt14+tt762+tt173*tt201+tt761+2*H_invDmH(5,1)&
&*H_invDmH(6,1)*tt7)+tt760+lam(1,1)*tt170*tt198)
tt811 = 2*H_invDmH(5,1)*H_invDmH(6,1)*tt6*tt8
tt812 = 2*H_invDmH(5,2)*H_invDmH(6,2)*tt13*tt15
tt813 = 2*H_invDmH(5,3)*H_invDmH(6,3)*tt20*tt22
tt814 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt176*tt212+tt175*tt211+tt8&
&13+tt173*tt209+tt812+tt811)+lam(1,1)*tt170*tt206)
tt815 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt176*tt221+tt175*tt220+tt7&
&86+tt173*tt217+tt785+tt784)+lam(1,1)*tt214*tt170)
tt816 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt782+2*H_invDmH(5,3)*H_inv&
&DmH(7,3)*tt21+tt781+tt780+tt176*tt232+tt175*tt231+tt779+2*H_invDm&
&H(5,2)*H_invDmH(7,2)*tt14+tt778+tt173*tt229+tt777+2*H_invDmH(5,1)&
&*H_invDmH(7,1)*tt7)+tt776+lam(1,1)*tt170*tt226)
tt817 = 2*H_invDmH(5,1)*H_invDmH(7,1)*tt6*tt8
tt818 = 2*H_invDmH(5,2)*H_invDmH(7,2)*tt13*tt15
tt819 = 2*H_invDmH(5,3)*H_invDmH(7,3)*tt20*tt22
tt820 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt176*tt240+tt175*tt239+tt8&
&19+tt173*tt237+tt818+tt817)+lam(1,1)*tt170*tt234)
tt821 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt176*tt250+tt175*tt248+tt8&
&02+tt173*tt245+tt801+tt800)+lam(1,1)*tt242*tt170)
tt822 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt798+2*H_invDmH(5,3)*H_inv&
&DmH(8,3)*tt21+tt797+tt176*tt260+tt796+tt175*tt259+tt795+2*H_invDm&
&H(5,2)*H_invDmH(8,2)*tt14+tt794+tt173*tt257+tt793+2*H_invDmH(5,1)&
&*H_invDmH(8,1)*tt7)+tt792+lam(1,1)*tt170*tt254)
tt823 = 2*H_invDmH(5,1)*H_invDmH(8,1)*tt6*tt8
tt824 = 2*H_invDmH(5,2)*H_invDmH(8,2)*tt13*tt15
tt825 = 2*H_invDmH(5,3)*H_invDmH(8,3)*tt20*tt22
tt826 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt176*tt268+tt175*tt267+tt8&
&25+tt173*tt265+tt824+tt823)+lam(1,1)*tt170*tt262)
tt827 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt184*tt193+tt183*tt192+tt7&
&74+tt181*tt189+tt773+tt772)+lam(1,1)*tt186*tt178)
tt828 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt184*tt204+tt183*tt203+tt8&
&13+tt181*tt201+tt812+tt811)+lam(1,1)*tt198*tt178)
tt829 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt766+2*H_invDmH(5,3)*H_inv&
&DmH(6,3)*tt23+tt765+tt764+tt184*tt212+tt183*tt211+tt763+2*H_invDm&
&H(5,2)*H_invDmH(6,2)*tt16+tt762+tt181*tt209+tt761+2*H_invDmH(5,1)&
&*H_invDmH(6,1)*tt9)+tt760+lam(1,1)*tt178*tt206)
tt830 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt184*tt221+tt183*tt220+tt7&
&90+tt181*tt217+tt789+tt788)+lam(1,1)*tt214*tt178)
tt831 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt184*tt232+tt183*tt231+tt8&
&19+tt181*tt229+tt818+tt817)+lam(1,1)*tt226*tt178)
tt832 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt782+2*H_invDmH(5,3)*H_inv&
&DmH(7,3)*tt23+tt781+tt780+tt184*tt240+tt183*tt239+tt779+2*H_invDm&
&H(5,2)*H_invDmH(7,2)*tt16+tt778+tt181*tt237+tt777+2*H_invDmH(5,1)&
&*H_invDmH(7,1)*tt9)+tt776+lam(1,1)*tt178*tt234)
tt833 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt184*tt250+tt183*tt248+tt8&
&06+tt181*tt245+tt805+tt804)+lam(1,1)*tt242*tt178)
tt834 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt184*tt260+tt183*tt259+tt8&
&25+tt181*tt257+tt824+tt823)+lam(1,1)*tt254*tt178)
tt835 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt798+2*H_invDmH(5,3)*H_inv&
&DmH(8,3)*tt23+tt797+tt184*tt268+tt796+tt183*tt267+tt795+2*H_invDm&
&H(5,2)*H_invDmH(8,2)*tt16+tt794+tt181*tt265+tt793+2*H_invDmH(5,1)&
&*H_invDmH(8,1)*tt9)+tt792+lam(1,1)*tt178*tt262)
tt836 = H_invDmH(6,1)**2
tt837 = H_invDmH(6,2)**2
tt838 = H_invDmH(6,3)**2
tt839 = lam(1,1)*(tt838+tt837+tt836)*tt25
tt840 = tt836*tt10
tt841 = 2*H_invDmH(6,1)*H_invDmH(6,2)*tt52
tt842 = tt837*tt17
tt843 = 2*H_invDmH(6,1)*H_invDmH(6,3)*tt56
tt844 = 2*H_invDmH(6,2)*H_invDmH(6,3)*tt58
tt845 = tt838*tt24
tt846 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt193*tt204+tt192*tt203+2*t&
&t838*tt18*tt20+tt189*tt201+2*tt837*tt11*tt13+2*tt836*tt4*tt6)+lam&
&(1,1)*tt186*tt198)
tt847 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt193*tt212+tt192*tt211+2*t&
&t838*tt18*tt22+tt189*tt209+2*tt837*tt11*tt15+2*tt836*tt4*tt8)+lam&
&(1,1)*tt186*tt206)
tt848 = lam(1,1)*(H_invDmH(6,3)*H_invDmH(7,3)+H_invDmH(6,2)*H_inv&
&DmH(7,2)+H_invDmH(6,1)*H_invDmH(7,1))*tt25
tt849 = H_invDmH(6,1)*H_invDmH(7,1)*tt10
tt850 = (H_invDmH(6,1)*H_invDmH(7,2)+H_invDmH(6,2)*H_invDmH(7,1))&
&*tt52
tt851 = H_invDmH(6,2)*H_invDmH(7,2)*tt17
tt852 = (H_invDmH(6,1)*H_invDmH(7,3)+H_invDmH(6,3)*H_invDmH(7,1))&
&*tt56
tt853 = (H_invDmH(6,2)*H_invDmH(7,3)+H_invDmH(6,3)*H_invDmH(7,2))&
&*tt58
tt854 = H_invDmH(6,3)*H_invDmH(7,3)*tt24
tt855 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt854+2*H_invDmH(6,3)*H_inv&
&DmH(7,3)*tt19+tt853+tt852+tt193*tt221+tt192*tt220+tt851+2*H_invDm&
&H(6,2)*H_invDmH(7,2)*tt12+tt850+tt189*tt217+tt849+2*H_invDmH(6,1)&
&*H_invDmH(7,1)*tt5)+tt848+lam(1,1)*tt186*tt214)
tt856 = 2*H_invDmH(6,1)*H_invDmH(7,1)*tt4*tt6
tt857 = 2*H_invDmH(6,2)*H_invDmH(7,2)*tt11*tt13
tt858 = 2*H_invDmH(6,3)*H_invDmH(7,3)*tt18*tt20
tt859 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt193*tt232+tt192*tt231+tt8&
&58+tt189*tt229+tt857+tt856)+lam(1,1)*tt186*tt226)
tt860 = 2*H_invDmH(6,1)*H_invDmH(7,1)*tt4*tt8
tt861 = 2*H_invDmH(6,2)*H_invDmH(7,2)*tt11*tt15
tt862 = 2*H_invDmH(6,3)*H_invDmH(7,3)*tt18*tt22
tt863 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt193*tt240+tt192*tt239+tt8&
&62+tt189*tt237+tt861+tt860)+lam(1,1)*tt186*tt234)
tt864 = lam(1,1)*(H_invDmH(6,3)*H_invDmH(8,3)+H_invDmH(6,2)*H_inv&
&DmH(8,2)+H_invDmH(6,1)*H_invDmH(8,1))*tt25
tt865 = H_invDmH(6,1)*H_invDmH(8,1)*tt10
tt866 = (H_invDmH(6,1)*H_invDmH(8,2)+H_invDmH(6,2)*H_invDmH(8,1))&
&*tt52
tt867 = H_invDmH(6,2)*H_invDmH(8,2)*tt17
tt868 = (H_invDmH(6,1)*H_invDmH(8,3)+H_invDmH(6,3)*H_invDmH(8,1))&
&*tt56
tt869 = (H_invDmH(6,2)*H_invDmH(8,3)+H_invDmH(6,3)*H_invDmH(8,2))&
&*tt58
tt870 = H_invDmH(6,3)*H_invDmH(8,3)*tt24
tt871 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt870+2*H_invDmH(6,3)*H_inv&
&DmH(8,3)*tt19+tt869+tt193*tt250+tt868+tt192*tt248+tt867+2*H_invDm&
&H(6,2)*H_invDmH(8,2)*tt12+tt866+tt189*tt245+tt865+2*H_invDmH(6,1)&
&*H_invDmH(8,1)*tt5)+tt864+lam(1,1)*tt186*tt242)
tt872 = 2*H_invDmH(6,1)*H_invDmH(8,1)*tt4*tt6
tt873 = 2*H_invDmH(6,2)*H_invDmH(8,2)*tt11*tt13
tt874 = 2*H_invDmH(6,3)*H_invDmH(8,3)*tt18*tt20
tt875 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt193*tt260+tt192*tt259+tt8&
&74+tt189*tt257+tt873+tt872)+lam(1,1)*tt186*tt254)
tt876 = 2*H_invDmH(6,1)*H_invDmH(8,1)*tt4*tt8
tt877 = 2*H_invDmH(6,2)*H_invDmH(8,2)*tt11*tt15
tt878 = 2*H_invDmH(6,3)*H_invDmH(8,3)*tt18*tt22
tt879 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt193*tt268+tt192*tt267+tt8&
&78+tt189*tt265+tt877+tt876)+lam(1,1)*tt186*tt262)
tt880 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt204*tt212+tt203*tt211+2*t&
&t838*tt20*tt22+tt201*tt209+2*tt837*tt13*tt15+2*tt836*tt6*tt8)+lam&
&(1,1)*tt198*tt206)
tt881 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt204*tt221+tt203*tt220+tt8&
&58+tt201*tt217+tt857+tt856)+lam(1,1)*tt214*tt198)
tt882 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt854+2*H_invDmH(6,3)*H_inv&
&DmH(7,3)*tt21+tt853+tt852+tt204*tt232+tt203*tt231+tt851+2*H_invDm&
&H(6,2)*H_invDmH(7,2)*tt14+tt850+tt201*tt229+tt849+2*H_invDmH(6,1)&
&*H_invDmH(7,1)*tt7)+tt848+lam(1,1)*tt198*tt226)
tt883 = 2*H_invDmH(6,1)*H_invDmH(7,1)*tt6*tt8
tt884 = 2*H_invDmH(6,2)*H_invDmH(7,2)*tt13*tt15
tt885 = 2*H_invDmH(6,3)*H_invDmH(7,3)*tt20*tt22
tt886 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt204*tt240+tt203*tt239+tt8&
&85+tt201*tt237+tt884+tt883)+lam(1,1)*tt198*tt234)
tt887 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt204*tt250+tt203*tt248+tt8&
&74+tt201*tt245+tt873+tt872)+lam(1,1)*tt242*tt198)
tt888 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt870+2*H_invDmH(6,3)*H_inv&
&DmH(8,3)*tt21+tt869+tt204*tt260+tt868+tt203*tt259+tt867+2*H_invDm&
&H(6,2)*H_invDmH(8,2)*tt14+tt866+tt201*tt257+tt865+2*H_invDmH(6,1)&
&*H_invDmH(8,1)*tt7)+tt864+lam(1,1)*tt198*tt254)
tt889 = 2*H_invDmH(6,1)*H_invDmH(8,1)*tt6*tt8
tt890 = 2*H_invDmH(6,2)*H_invDmH(8,2)*tt13*tt15
tt891 = 2*H_invDmH(6,3)*H_invDmH(8,3)*tt20*tt22
tt892 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt204*tt268+tt203*tt267+tt8&
&91+tt201*tt265+tt890+tt889)+lam(1,1)*tt198*tt262)
tt893 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt212*tt221+tt211*tt220+tt8&
&62+tt209*tt217+tt861+tt860)+lam(1,1)*tt214*tt206)
tt894 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt212*tt232+tt211*tt231+tt8&
&85+tt209*tt229+tt884+tt883)+lam(1,1)*tt226*tt206)
tt895 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt854+2*H_invDmH(6,3)*H_inv&
&DmH(7,3)*tt23+tt853+tt852+tt212*tt240+tt211*tt239+tt851+2*H_invDm&
&H(6,2)*H_invDmH(7,2)*tt16+tt850+tt209*tt237+tt849+2*H_invDmH(6,1)&
&*H_invDmH(7,1)*tt9)+tt848+lam(1,1)*tt206*tt234)
tt896 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt212*tt250+tt211*tt248+tt8&
&78+tt209*tt245+tt877+tt876)+lam(1,1)*tt242*tt206)
tt897 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt212*tt260+tt211*tt259+tt8&
&91+tt209*tt257+tt890+tt889)+lam(1,1)*tt254*tt206)
tt898 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt870+2*H_invDmH(6,3)*H_inv&
&DmH(8,3)*tt23+tt869+tt212*tt268+tt868+tt211*tt267+tt867+2*H_invDm&
&H(6,2)*H_invDmH(8,2)*tt16+tt866+tt209*tt265+tt865+2*H_invDmH(6,1)&
&*H_invDmH(8,1)*tt9)+tt864+lam(1,1)*tt206*tt262)
tt899 = H_invDmH(7,1)**2
tt900 = H_invDmH(7,2)**2
tt901 = H_invDmH(7,3)**2
tt902 = lam(1,1)*(tt901+tt900+tt899)*tt25
tt903 = tt899*tt10
tt904 = 2*H_invDmH(7,1)*H_invDmH(7,2)*tt52
tt905 = tt900*tt17
tt906 = 2*H_invDmH(7,1)*H_invDmH(7,3)*tt56
tt907 = 2*H_invDmH(7,2)*H_invDmH(7,3)*tt58
tt908 = tt901*tt24
tt909 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt221*tt232+tt220*tt231+2*t&
&t901*tt18*tt20+tt217*tt229+2*tt900*tt11*tt13+2*tt899*tt4*tt6)+lam&
&(1,1)*tt214*tt226)
tt910 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt221*tt240+tt220*tt239+2*t&
&t901*tt18*tt22+tt217*tt237+2*tt900*tt11*tt15+2*tt899*tt4*tt8)+lam&
&(1,1)*tt214*tt234)
tt911 = lam(1,1)*(H_invDmH(7,3)*H_invDmH(8,3)+H_invDmH(7,2)*H_inv&
&DmH(8,2)+H_invDmH(7,1)*H_invDmH(8,1))*tt25
tt912 = H_invDmH(7,1)*H_invDmH(8,1)*tt10
tt913 = (H_invDmH(7,1)*H_invDmH(8,2)+H_invDmH(7,2)*H_invDmH(8,1))&
&*tt52
tt914 = H_invDmH(7,2)*H_invDmH(8,2)*tt17
tt915 = (H_invDmH(7,1)*H_invDmH(8,3)+H_invDmH(7,3)*H_invDmH(8,1))&
&*tt56
tt916 = (H_invDmH(7,2)*H_invDmH(8,3)+H_invDmH(7,3)*H_invDmH(8,2))&
&*tt58
tt917 = H_invDmH(7,3)*H_invDmH(8,3)*tt24
tt918 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt917+2*H_invDmH(7,3)*H_inv&
&DmH(8,3)*tt19+tt916+tt221*tt250+tt915+tt220*tt248+tt914+2*H_invDm&
&H(7,2)*H_invDmH(8,2)*tt12+tt913+tt217*tt245+tt912+2*H_invDmH(7,1)&
&*H_invDmH(8,1)*tt5)+tt911+lam(1,1)*tt214*tt242)
tt919 = 2*H_invDmH(7,1)*H_invDmH(8,1)*tt4*tt6
tt920 = 2*H_invDmH(7,2)*H_invDmH(8,2)*tt11*tt13
tt921 = 2*H_invDmH(7,3)*H_invDmH(8,3)*tt18*tt20
tt922 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt221*tt260+tt220*tt259+tt9&
&21+tt217*tt257+tt920+tt919)+lam(1,1)*tt214*tt254)
tt923 = 2*H_invDmH(7,1)*H_invDmH(8,1)*tt4*tt8
tt924 = 2*H_invDmH(7,2)*H_invDmH(8,2)*tt11*tt15
tt925 = 2*H_invDmH(7,3)*H_invDmH(8,3)*tt18*tt22
tt926 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt221*tt268+tt220*tt267+tt9&
&25+tt217*tt265+tt924+tt923)+lam(1,1)*tt214*tt262)
tt927 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt232*tt240+tt231*tt239+2*t&
&t901*tt20*tt22+tt229*tt237+2*tt900*tt13*tt15+2*tt899*tt6*tt8)+lam&
&(1,1)*tt226*tt234)
tt928 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt232*tt250+tt231*tt248+tt9&
&21+tt229*tt245+tt920+tt919)+lam(1,1)*tt242*tt226)
tt929 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt917+2*H_invDmH(7,3)*H_inv&
&DmH(8,3)*tt21+tt916+tt232*tt260+tt915+tt231*tt259+tt914+2*H_invDm&
&H(7,2)*H_invDmH(8,2)*tt14+tt913+tt229*tt257+tt912+2*H_invDmH(7,1)&
&*H_invDmH(8,1)*tt7)+tt911+lam(1,1)*tt226*tt254)
tt930 = 2*H_invDmH(7,1)*H_invDmH(8,1)*tt6*tt8
tt931 = 2*H_invDmH(7,2)*H_invDmH(8,2)*tt13*tt15
tt932 = 2*H_invDmH(7,3)*H_invDmH(8,3)*tt20*tt22
tt933 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt232*tt268+tt231*tt267+tt9&
&32+tt229*tt265+tt931+tt930)+lam(1,1)*tt226*tt262)
tt934 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt240*tt250+tt239*tt248+tt9&
&25+tt237*tt245+tt924+tt923)+lam(1,1)*tt242*tt234)
tt935 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt240*tt260+tt239*tt259+tt9&
&32+tt237*tt257+tt931+tt930)+lam(1,1)*tt254*tt234)
tt936 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt917+2*H_invDmH(7,3)*H_inv&
&DmH(8,3)*tt23+tt916+tt240*tt268+tt915+tt239*tt267+tt914+2*H_invDm&
&H(7,2)*H_invDmH(8,2)*tt16+tt913+tt237*tt265+tt912+2*H_invDmH(7,1)&
&*H_invDmH(8,1)*tt9)+tt911+lam(1,1)*tt234*tt262)
tt937 = H_invDmH(8,1)**2
tt938 = H_invDmH(8,2)**2
tt939 = H_invDmH(8,3)**2
tt940 = lam(1,1)*(tt939+tt938+tt937)*tt25
tt941 = tt937*tt10
tt942 = 2*H_invDmH(8,1)*H_invDmH(8,2)*tt52
tt943 = tt938*tt17
tt944 = 2*H_invDmH(8,1)*H_invDmH(8,3)*tt56
tt945 = 2*H_invDmH(8,2)*H_invDmH(8,3)*tt58
tt946 = tt939*tt24
tt947 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt250*tt260+tt248*tt259+2*t&
&t939*tt18*tt20+tt245*tt257+2*tt938*tt11*tt13+2*tt937*tt4*tt6)+lam&
&(1,1)*tt242*tt254)
tt948 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt250*tt268+tt248*tt267+2*t&
&t939*tt18*tt22+tt245*tt265+2*tt938*tt11*tt15+2*tt937*tt4*tt8)+lam&
&(1,1)*tt242*tt262)
tt949 = detDmH(1,1)*gw(1,1)*(mu(1,1)*(tt260*tt268+tt259*tt267+2*t&
&t939*tt20*tt22+tt257*tt265+2*tt938*tt13*tt15+2*tt937*tt6*tt8)+lam&
&(1,1)*tt254*tt262)
hes(1,1) = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt63**2+mu(1,1)*(tt62+tt&
&61**2+tt60**2+2*tt3*tt19+tt59+tt57+tt55+tt54**2+2*tt2*tt12+tt53+t&
&t27+2*tt1*tt5)+tt26)
hes(1,2) = tt68
hes(1,3) = tt73
hes(1,4) = tt85
hes(1,5) = tt93
hes(1,6) = tt101
hes(1,7) = tt113
hes(1,8) = tt121
hes(1,9) = tt129
hes(1,10) = tt141
hes(1,11) = tt149
hes(1,12) = tt157
hes(1,13) = tt169
hes(1,14) = tt177
hes(1,15) = tt185
hes(1,16) = tt197
hes(1,17) = tt205
hes(1,18) = tt213
hes(1,19) = tt225
hes(1,20) = tt233
hes(1,21) = tt241
hes(1,22) = tt253
hes(1,23) = tt261
hes(1,24) = tt269
hes(2,1) = tt68
hes(2,2) = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt64**2+mu(1,1)*(tt62+tt&
&67**2+tt66**2+2*tt3*tt21+tt59+tt57+tt55+tt65**2+2*tt2*tt14+tt53+t&
&t27+2*tt1*tt7)+tt26)
hes(2,3) = tt270
hes(2,4) = tt271
hes(2,5) = tt272
hes(2,6) = tt276
hes(2,7) = tt277
hes(2,8) = tt278
hes(2,9) = tt282
hes(2,10) = tt283
hes(2,11) = tt284
hes(2,12) = tt288
hes(2,13) = tt289
hes(2,14) = tt290
hes(2,15) = tt294
hes(2,16) = tt295
hes(2,17) = tt296
hes(2,18) = tt300
hes(2,19) = tt301
hes(2,20) = tt302
hes(2,21) = tt306
hes(2,22) = tt307
hes(2,23) = tt308
hes(2,24) = tt312
hes(3,1) = tt73
hes(3,2) = tt270
hes(3,3) = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt69**2+mu(1,1)*(tt72**2&
&+tt71**2+tt62+2*tt3*tt23+tt59+tt57+tt70**2+tt55+2*tt2*tt16+tt53+t&
&t27+2*tt1*tt9)+tt26)
hes(3,4) = tt313
hes(3,5) = tt314
hes(3,6) = tt315
hes(3,7) = tt316
hes(3,8) = tt317
hes(3,9) = tt318
hes(3,10) = tt319
hes(3,11) = tt320
hes(3,12) = tt321
hes(3,13) = tt322
hes(3,14) = tt323
hes(3,15) = tt324
hes(3,16) = tt325
hes(3,17) = tt326
hes(3,18) = tt327
hes(3,19) = tt328
hes(3,20) = tt329
hes(3,21) = tt330
hes(3,22) = tt331
hes(3,23) = tt332
hes(3,24) = tt333
hes(4,1) = tt85
hes(4,2) = tt271
hes(4,3) = tt313
hes(4,4) = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt74**2+mu(1,1)*(tt343+t&
&t81**2+tt80**2+2*tt336*tt19+tt342+tt341+tt340+tt77**2+2*tt335*tt1&
&2+tt339+tt338+2*tt334*tt5)+tt337)
hes(4,5) = tt344
hes(4,6) = tt345
hes(4,7) = tt353
hes(4,8) = tt357
hes(4,9) = tt361
hes(4,10) = tt369
hes(4,11) = tt373
hes(4,12) = tt377
hes(4,13) = tt385
hes(4,14) = tt389
hes(4,15) = tt393
hes(4,16) = tt401
hes(4,17) = tt405
hes(4,18) = tt409
hes(4,19) = tt417
hes(4,20) = tt421
hes(4,21) = tt425
hes(4,22) = tt433
hes(4,23) = tt437
hes(4,24) = tt441
hes(5,1) = tt93
hes(5,2) = tt272
hes(5,3) = tt314
hes(5,4) = tt344
hes(5,5) = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt86**2+mu(1,1)*(tt343+t&
&t92**2+tt91**2+2*tt336*tt21+tt342+tt341+tt340+tt89**2+2*tt335*tt1&
&4+tt339+tt338+2*tt334*tt7)+tt337)
hes(5,6) = tt442
hes(5,7) = tt443
hes(5,8) = tt444
hes(5,9) = tt448
hes(5,10) = tt449
hes(5,11) = tt450
hes(5,12) = tt454
hes(5,13) = tt455
hes(5,14) = tt456
hes(5,15) = tt460
hes(5,16) = tt461
hes(5,17) = tt462
hes(5,18) = tt466
hes(5,19) = tt467
hes(5,20) = tt468
hes(5,21) = tt472
hes(5,22) = tt473
hes(5,23) = tt474
hes(5,24) = tt478
hes(6,1) = tt101
hes(6,2) = tt276
hes(6,3) = tt315
hes(6,4) = tt345
hes(6,5) = tt442
hes(6,6) = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt94**2+mu(1,1)*(tt100**&
&2+tt99**2+tt343+2*tt336*tt23+tt342+tt341+tt97**2+tt340+2*tt335*tt&
&16+tt339+tt338+2*tt334*tt9)+tt337)
hes(6,7) = tt479
hes(6,8) = tt480
hes(6,9) = tt481
hes(6,10) = tt482
hes(6,11) = tt483
hes(6,12) = tt484
hes(6,13) = tt485
hes(6,14) = tt486
hes(6,15) = tt487
hes(6,16) = tt488
hes(6,17) = tt489
hes(6,18) = tt490
hes(6,19) = tt491
hes(6,20) = tt492
hes(6,21) = tt493
hes(6,22) = tt494
hes(6,23) = tt495
hes(6,24) = tt496
hes(7,1) = tt113
hes(7,2) = tt277
hes(7,3) = tt316
hes(7,4) = tt353
hes(7,5) = tt443
hes(7,6) = tt479
hes(7,7) = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt102**2+mu(1,1)*(tt506+&
&tt109**2+tt108**2+2*tt499*tt19+tt505+tt504+tt503+tt105**2+2*tt498&
&*tt12+tt502+tt501+2*tt497*tt5)+tt500)
hes(7,8) = tt507
hes(7,9) = tt508
hes(7,10) = tt516
hes(7,11) = tt520
hes(7,12) = tt524
hes(7,13) = tt532
hes(7,14) = tt536
hes(7,15) = tt540
hes(7,16) = tt548
hes(7,17) = tt552
hes(7,18) = tt556
hes(7,19) = tt564
hes(7,20) = tt568
hes(7,21) = tt572
hes(7,22) = tt580
hes(7,23) = tt584
hes(7,24) = tt588
hes(8,1) = tt121
hes(8,2) = tt278
hes(8,3) = tt317
hes(8,4) = tt357
hes(8,5) = tt444
hes(8,6) = tt480
hes(8,7) = tt507
hes(8,8) = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt114**2+mu(1,1)*(tt506+&
&tt120**2+tt119**2+2*tt499*tt21+tt505+tt504+tt503+tt117**2+2*tt498&
&*tt14+tt502+tt501+2*tt497*tt7)+tt500)
hes(8,9) = tt589
hes(8,10) = tt590
hes(8,11) = tt591
hes(8,12) = tt595
hes(8,13) = tt596
hes(8,14) = tt597
hes(8,15) = tt601
hes(8,16) = tt602
hes(8,17) = tt603
hes(8,18) = tt607
hes(8,19) = tt608
hes(8,20) = tt609
hes(8,21) = tt613
hes(8,22) = tt614
hes(8,23) = tt615
hes(8,24) = tt619
hes(9,1) = tt129
hes(9,2) = tt282
hes(9,3) = tt318
hes(9,4) = tt361
hes(9,5) = tt448
hes(9,6) = tt481
hes(9,7) = tt508
hes(9,8) = tt589
hes(9,9) = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt122**2+mu(1,1)*(tt128*&
&*2+tt127**2+tt506+2*tt499*tt23+tt505+tt504+tt125**2+tt503+2*tt498&
&*tt16+tt502+tt501+2*tt497*tt9)+tt500)
hes(9,10) = tt620
hes(9,11) = tt621
hes(9,12) = tt622
hes(9,13) = tt623
hes(9,14) = tt624
hes(9,15) = tt625
hes(9,16) = tt626
hes(9,17) = tt627
hes(9,18) = tt628
hes(9,19) = tt629
hes(9,20) = tt630
hes(9,21) = tt631
hes(9,22) = tt632
hes(9,23) = tt633
hes(9,24) = tt634
hes(10,1) = tt141
hes(10,2) = tt283
hes(10,3) = tt319
hes(10,4) = tt369
hes(10,5) = tt449
hes(10,6) = tt482
hes(10,7) = tt516
hes(10,8) = tt590
hes(10,9) = tt620
hes(10,10) = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt130**2+mu(1,1)*(tt13&
&7**2+tt136**2+tt644+2*tt637*tt19+tt643+tt642+tt133**2+tt641+2*tt6&
&36*tt12+tt640+tt639+2*tt635*tt5)+tt638)
hes(10,11) = tt645
hes(10,12) = tt646
hes(10,13) = tt654
hes(10,14) = tt658
hes(10,15) = tt662
hes(10,16) = tt670
hes(10,17) = tt674
hes(10,18) = tt678
hes(10,19) = tt686
hes(10,20) = tt690
hes(10,21) = tt694
hes(10,22) = tt702
hes(10,23) = tt706
hes(10,24) = tt710
hes(11,1) = tt149
hes(11,2) = tt284
hes(11,3) = tt320
hes(11,4) = tt373
hes(11,5) = tt450
hes(11,6) = tt483
hes(11,7) = tt520
hes(11,8) = tt591
hes(11,9) = tt621
hes(11,10) = tt645
hes(11,11) = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt142**2+mu(1,1)*(tt14&
&8**2+tt147**2+tt644+2*tt637*tt21+tt643+tt642+tt145**2+tt641+2*tt6&
&36*tt14+tt640+tt639+2*tt635*tt7)+tt638)
hes(11,12) = tt711
hes(11,13) = tt712
hes(11,14) = tt713
hes(11,15) = tt717
hes(11,16) = tt718
hes(11,17) = tt719
hes(11,18) = tt723
hes(11,19) = tt724
hes(11,20) = tt725
hes(11,21) = tt729
hes(11,22) = tt730
hes(11,23) = tt731
hes(11,24) = tt735
hes(12,1) = tt157
hes(12,2) = tt288
hes(12,3) = tt321
hes(12,4) = tt377
hes(12,5) = tt454
hes(12,6) = tt484
hes(12,7) = tt524
hes(12,8) = tt595
hes(12,9) = tt622
hes(12,10) = tt646
hes(12,11) = tt711
hes(12,12) = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt150**2+mu(1,1)*(tt15&
&6**2+tt155**2+tt644+2*tt637*tt23+tt643+tt642+tt153**2+tt641+2*tt6&
&36*tt16+tt640+tt639+2*tt635*tt9)+tt638)
hes(12,13) = tt736
hes(12,14) = tt737
hes(12,15) = tt738
hes(12,16) = tt739
hes(12,17) = tt740
hes(12,18) = tt741
hes(12,19) = tt742
hes(12,20) = tt743
hes(12,21) = tt744
hes(12,22) = tt745
hes(12,23) = tt746
hes(12,24) = tt747
hes(13,1) = tt169
hes(13,2) = tt289
hes(13,3) = tt322
hes(13,4) = tt385
hes(13,5) = tt455
hes(13,6) = tt485
hes(13,7) = tt532
hes(13,8) = tt596
hes(13,9) = tt623
hes(13,10) = tt654
hes(13,11) = tt712
hes(13,12) = tt736
hes(13,13) = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt158**2+mu(1,1)*(tt16&
&5**2+tt164**2+tt757+2*tt750*tt19+tt756+tt755+tt161**2+tt754+2*tt7&
&49*tt12+tt753+tt752+2*tt748*tt5)+tt751)
hes(13,14) = tt758
hes(13,15) = tt759
hes(13,16) = tt767
hes(13,17) = tt771
hes(13,18) = tt775
hes(13,19) = tt783
hes(13,20) = tt787
hes(13,21) = tt791
hes(13,22) = tt799
hes(13,23) = tt803
hes(13,24) = tt807
hes(14,1) = tt177
hes(14,2) = tt290
hes(14,3) = tt323
hes(14,4) = tt389
hes(14,5) = tt456
hes(14,6) = tt486
hes(14,7) = tt536
hes(14,8) = tt597
hes(14,9) = tt624
hes(14,10) = tt658
hes(14,11) = tt713
hes(14,12) = tt737
hes(14,13) = tt758
hes(14,14) = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt170**2+mu(1,1)*(tt17&
&6**2+tt175**2+tt757+2*tt750*tt21+tt756+tt755+tt173**2+tt754+2*tt7&
&49*tt14+tt753+tt752+2*tt748*tt7)+tt751)
hes(14,15) = tt808
hes(14,16) = tt809
hes(14,17) = tt810
hes(14,18) = tt814
hes(14,19) = tt815
hes(14,20) = tt816
hes(14,21) = tt820
hes(14,22) = tt821
hes(14,23) = tt822
hes(14,24) = tt826
hes(15,1) = tt185
hes(15,2) = tt294
hes(15,3) = tt324
hes(15,4) = tt393
hes(15,5) = tt460
hes(15,6) = tt487
hes(15,7) = tt540
hes(15,8) = tt601
hes(15,9) = tt625
hes(15,10) = tt662
hes(15,11) = tt717
hes(15,12) = tt738
hes(15,13) = tt759
hes(15,14) = tt808
hes(15,15) = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt178**2+mu(1,1)*(tt18&
&4**2+tt183**2+tt757+2*tt750*tt23+tt756+tt755+tt181**2+tt754+2*tt7&
&49*tt16+tt753+tt752+2*tt748*tt9)+tt751)
hes(15,16) = tt827
hes(15,17) = tt828
hes(15,18) = tt829
hes(15,19) = tt830
hes(15,20) = tt831
hes(15,21) = tt832
hes(15,22) = tt833
hes(15,23) = tt834
hes(15,24) = tt835
hes(16,1) = tt197
hes(16,2) = tt295
hes(16,3) = tt325
hes(16,4) = tt401
hes(16,5) = tt461
hes(16,6) = tt488
hes(16,7) = tt548
hes(16,8) = tt602
hes(16,9) = tt626
hes(16,10) = tt670
hes(16,11) = tt718
hes(16,12) = tt739
hes(16,13) = tt767
hes(16,14) = tt809
hes(16,15) = tt827
hes(16,16) = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt186**2+mu(1,1)*(tt19&
&3**2+tt192**2+tt845+2*tt838*tt19+tt844+tt843+tt189**2+tt842+2*tt8&
&37*tt12+tt841+tt840+2*tt836*tt5)+tt839)
hes(16,17) = tt846
hes(16,18) = tt847
hes(16,19) = tt855
hes(16,20) = tt859
hes(16,21) = tt863
hes(16,22) = tt871
hes(16,23) = tt875
hes(16,24) = tt879
hes(17,1) = tt205
hes(17,2) = tt296
hes(17,3) = tt326
hes(17,4) = tt405
hes(17,5) = tt462
hes(17,6) = tt489
hes(17,7) = tt552
hes(17,8) = tt603
hes(17,9) = tt627
hes(17,10) = tt674
hes(17,11) = tt719
hes(17,12) = tt740
hes(17,13) = tt771
hes(17,14) = tt810
hes(17,15) = tt828
hes(17,16) = tt846
hes(17,17) = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt198**2+mu(1,1)*(tt20&
&4**2+tt203**2+tt845+2*tt838*tt21+tt844+tt843+tt201**2+tt842+2*tt8&
&37*tt14+tt841+tt840+2*tt836*tt7)+tt839)
hes(17,18) = tt880
hes(17,19) = tt881
hes(17,20) = tt882
hes(17,21) = tt886
hes(17,22) = tt887
hes(17,23) = tt888
hes(17,24) = tt892
hes(18,1) = tt213
hes(18,2) = tt300
hes(18,3) = tt327
hes(18,4) = tt409
hes(18,5) = tt466
hes(18,6) = tt490
hes(18,7) = tt556
hes(18,8) = tt607
hes(18,9) = tt628
hes(18,10) = tt678
hes(18,11) = tt723
hes(18,12) = tt741
hes(18,13) = tt775
hes(18,14) = tt814
hes(18,15) = tt829
hes(18,16) = tt847
hes(18,17) = tt880
hes(18,18) = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt206**2+mu(1,1)*(tt21&
&2**2+tt211**2+tt845+2*tt838*tt23+tt844+tt843+tt209**2+tt842+2*tt8&
&37*tt16+tt841+tt840+2*tt836*tt9)+tt839)
hes(18,19) = tt893
hes(18,20) = tt894
hes(18,21) = tt895
hes(18,22) = tt896
hes(18,23) = tt897
hes(18,24) = tt898
hes(19,1) = tt225
hes(19,2) = tt301
hes(19,3) = tt328
hes(19,4) = tt417
hes(19,5) = tt467
hes(19,6) = tt491
hes(19,7) = tt564
hes(19,8) = tt608
hes(19,9) = tt629
hes(19,10) = tt686
hes(19,11) = tt724
hes(19,12) = tt742
hes(19,13) = tt783
hes(19,14) = tt815
hes(19,15) = tt830
hes(19,16) = tt855
hes(19,17) = tt881
hes(19,18) = tt893
hes(19,19) = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt214**2+mu(1,1)*(tt22&
&1**2+tt220**2+tt908+2*tt901*tt19+tt907+tt906+tt217**2+tt905+2*tt9&
&00*tt12+tt904+tt903+2*tt899*tt5)+tt902)
hes(19,20) = tt909
hes(19,21) = tt910
hes(19,22) = tt918
hes(19,23) = tt922
hes(19,24) = tt926
hes(20,1) = tt233
hes(20,2) = tt302
hes(20,3) = tt329
hes(20,4) = tt421
hes(20,5) = tt468
hes(20,6) = tt492
hes(20,7) = tt568
hes(20,8) = tt609
hes(20,9) = tt630
hes(20,10) = tt690
hes(20,11) = tt725
hes(20,12) = tt743
hes(20,13) = tt787
hes(20,14) = tt816
hes(20,15) = tt831
hes(20,16) = tt859
hes(20,17) = tt882
hes(20,18) = tt894
hes(20,19) = tt909
hes(20,20) = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt226**2+mu(1,1)*(tt23&
&2**2+tt231**2+tt908+2*tt901*tt21+tt907+tt906+tt229**2+tt905+2*tt9&
&00*tt14+tt904+tt903+2*tt899*tt7)+tt902)
hes(20,21) = tt927
hes(20,22) = tt928
hes(20,23) = tt929
hes(20,24) = tt933
hes(21,1) = tt241
hes(21,2) = tt306
hes(21,3) = tt330
hes(21,4) = tt425
hes(21,5) = tt472
hes(21,6) = tt493
hes(21,7) = tt572
hes(21,8) = tt613
hes(21,9) = tt631
hes(21,10) = tt694
hes(21,11) = tt729
hes(21,12) = tt744
hes(21,13) = tt791
hes(21,14) = tt820
hes(21,15) = tt832
hes(21,16) = tt863
hes(21,17) = tt886
hes(21,18) = tt895
hes(21,19) = tt910
hes(21,20) = tt927
hes(21,21) = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt234**2+mu(1,1)*(tt24&
&0**2+tt239**2+tt908+2*tt901*tt23+tt907+tt906+tt237**2+tt905+2*tt9&
&00*tt16+tt904+tt903+2*tt899*tt9)+tt902)
hes(21,22) = tt934
hes(21,23) = tt935
hes(21,24) = tt936
hes(22,1) = tt253
hes(22,2) = tt307
hes(22,3) = tt331
hes(22,4) = tt433
hes(22,5) = tt473
hes(22,6) = tt494
hes(22,7) = tt580
hes(22,8) = tt614
hes(22,9) = tt632
hes(22,10) = tt702
hes(22,11) = tt730
hes(22,12) = tt745
hes(22,13) = tt799
hes(22,14) = tt821
hes(22,15) = tt833
hes(22,16) = tt871
hes(22,17) = tt887
hes(22,18) = tt896
hes(22,19) = tt918
hes(22,20) = tt928
hes(22,21) = tt934
hes(22,22) = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt242**2+mu(1,1)*(tt25&
&0**2+tt248**2+tt946+2*tt939*tt19+tt945+tt944+tt245**2+tt943+2*tt9&
&38*tt12+tt942+tt941+2*tt937*tt5)+tt940)
hes(22,23) = tt947
hes(22,24) = tt948
hes(23,1) = tt261
hes(23,2) = tt308
hes(23,3) = tt332
hes(23,4) = tt437
hes(23,5) = tt474
hes(23,6) = tt495
hes(23,7) = tt584
hes(23,8) = tt615
hes(23,9) = tt633
hes(23,10) = tt706
hes(23,11) = tt731
hes(23,12) = tt746
hes(23,13) = tt803
hes(23,14) = tt822
hes(23,15) = tt834
hes(23,16) = tt875
hes(23,17) = tt888
hes(23,18) = tt897
hes(23,19) = tt922
hes(23,20) = tt929
hes(23,21) = tt935
hes(23,22) = tt947
hes(23,23) = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt254**2+mu(1,1)*(tt26&
&0**2+tt259**2+tt946+2*tt939*tt21+tt945+tt944+tt257**2+tt943+2*tt9&
&38*tt14+tt942+tt941+2*tt937*tt7)+tt940)
hes(23,24) = tt949
hes(24,1) = tt269
hes(24,2) = tt312
hes(24,3) = tt333
hes(24,4) = tt441
hes(24,5) = tt478
hes(24,6) = tt496
hes(24,7) = tt588
hes(24,8) = tt619
hes(24,9) = tt634
hes(24,10) = tt710
hes(24,11) = tt735
hes(24,12) = tt747
hes(24,13) = tt807
hes(24,14) = tt826
hes(24,15) = tt835
hes(24,16) = tt879
hes(24,17) = tt892
hes(24,18) = tt898
hes(24,19) = tt926
hes(24,20) = tt933
hes(24,21) = tt936
hes(24,22) = tt948
hes(24,23) = tt949
hes(24,24) = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt262**2+mu(1,1)*(tt26&
&8**2+tt267**2+tt946+2*tt939*tt23+tt945+tt944+tt265**2+tt943+2*tt9&
&38*tt16+tt942+tt941+2*tt937*tt9)+tt940)
END 
SUBROUTINE vox_neo_at_quadr(val, X, H_invDmH, detDmH, gw, mu, lam) 
IMPLICIT NONE 
REAL(KIND=8) val(1, 1) 
REAL(KIND=8) X(3, 8) 
REAL(KIND=8) H_invDmH(8, 3) 
REAL(KIND=8) detDmH(1, 1) 
REAL(KIND=8) gw(1, 1) 
REAL(KIND=8) mu(1, 1) 
REAL(KIND=8) lam(1, 1) 
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
val(1,1) = detDmH(1,1)*gw(1,1)*((lam(1,1)*tt10**2)/2.0E+0+(mu(1,1&
&)*((-2*tt10)+tt9**2+tt8**2+tt7**2+tt6**2+tt5**2+tt4**2+tt3**2+tt2&
&**2+tt1**2-3))/2.0E+0)
END 
SUBROUTINE vox_neo_at_quadr_jac(jac, X, H_invDmH, detDmH, gw, mu, lam) 
IMPLICIT NONE 
REAL(KIND=8) jac(1, 24) 
REAL(KIND=8) X(3, 8) 
REAL(KIND=8) H_invDmH(8, 3) 
REAL(KIND=8) detDmH(1, 1) 
REAL(KIND=8) gw(1, 1) 
REAL(KIND=8) mu(1, 1) 
REAL(KIND=8) lam(1, 1) 
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
tt1 = X(1,8)*H_invDmH(8,1)+X(1,7)*H_invDmH(7,1)+X(1,6)*H_invDmH(6&
&,1)+X(1,5)*H_invDmH(5,1)+X(1,4)*H_invDmH(4,1)+X(1,3)*H_invDmH(3,1&
&)+X(1,2)*H_invDmH(2,1)+H_invDmH(1,1)*X(1,1)
tt2 = X(1,8)*H_invDmH(8,2)+X(1,7)*H_invDmH(7,2)+X(1,6)*H_invDmH(6&
&,2)+X(1,5)*H_invDmH(5,2)+X(1,4)*H_invDmH(4,2)+X(1,3)*H_invDmH(3,2&
&)+X(1,2)*H_invDmH(2,2)+X(1,1)*H_invDmH(1,2)
tt3 = X(1,8)*H_invDmH(8,3)+X(1,7)*H_invDmH(7,3)+X(1,6)*H_invDmH(6&
&,3)+X(1,5)*H_invDmH(5,3)+X(1,4)*H_invDmH(4,3)+X(1,3)*H_invDmH(3,3&
&)+X(1,2)*H_invDmH(2,3)+X(1,1)*H_invDmH(1,3)
tt4 = X(3,8)*H_invDmH(8,1)+X(3,7)*H_invDmH(7,1)+X(3,6)*H_invDmH(6&
&,1)+X(3,5)*H_invDmH(5,1)+X(3,4)*H_invDmH(4,1)+H_invDmH(3,1)*X(3,3&
&)+H_invDmH(2,1)*X(3,2)+H_invDmH(1,1)*X(3,1)
tt5 = X(2,8)*H_invDmH(8,2)+X(2,7)*H_invDmH(7,2)+X(2,6)*H_invDmH(6&
&,2)+X(2,5)*H_invDmH(5,2)+X(2,4)*H_invDmH(4,2)+X(2,3)*H_invDmH(3,2&
&)+H_invDmH(2,2)*X(2,2)+H_invDmH(1,2)*X(2,1)
tt6 = X(2,8)*H_invDmH(8,1)+X(2,7)*H_invDmH(7,1)+X(2,6)*H_invDmH(6&
&,1)+X(2,5)*H_invDmH(5,1)+X(2,4)*H_invDmH(4,1)+X(2,3)*H_invDmH(3,1&
&)+H_invDmH(2,1)*X(2,2)+H_invDmH(1,1)*X(2,1)
tt7 = X(3,8)*H_invDmH(8,2)+X(3,7)*H_invDmH(7,2)+X(3,6)*H_invDmH(6&
&,2)+X(3,5)*H_invDmH(5,2)+X(3,4)*H_invDmH(4,2)+H_invDmH(3,2)*X(3,3&
&)+H_invDmH(2,2)*X(3,2)+H_invDmH(1,2)*X(3,1)
tt8 = tt6*tt7-tt4*tt5
tt9 = X(2,8)*H_invDmH(8,3)+X(2,7)*H_invDmH(7,3)+X(2,6)*H_invDmH(6&
&,3)+X(2,5)*H_invDmH(5,3)+X(2,4)*H_invDmH(4,3)+X(2,3)*H_invDmH(3,3&
&)+X(2,2)*H_invDmH(2,3)+H_invDmH(1,3)*X(2,1)
tt10 = X(3,8)*H_invDmH(8,3)+X(3,7)*H_invDmH(7,3)+X(3,6)*H_invDmH(&
&6,3)+X(3,5)*H_invDmH(5,3)+X(3,4)*H_invDmH(4,3)+H_invDmH(3,3)*X(3,&
&3)+H_invDmH(2,3)*X(3,2)+H_invDmH(1,3)*X(3,1)
tt11 = tt6*tt10-tt4*tt9
tt12 = tt5*tt10-tt7*tt9
tt13 = H_invDmH(1,1)*tt12-H_invDmH(1,2)*tt11+H_invDmH(1,3)*tt8
tt14 = tt1*tt12-tt2*tt11+tt8*tt3
tt15 = 1/tt14
tt16 = log(tt14)
tt17 = tt1*(H_invDmH(1,2)*tt10-H_invDmH(1,3)*tt7)-tt2*(H_invDmH(1&
&,1)*tt10-H_invDmH(1,3)*tt4)+(H_invDmH(1,1)*tt7-H_invDmH(1,2)*tt4)&
&*tt3
tt18 = tt1*(H_invDmH(1,3)*tt5-H_invDmH(1,2)*tt9)-tt2*(H_invDmH(1,&
&3)*tt6-H_invDmH(1,1)*tt9)+(H_invDmH(1,2)*tt6-H_invDmH(1,1)*tt5)*t&
&t3
tt19 = H_invDmH(2,1)*tt12-H_invDmH(2,2)*tt11+H_invDmH(2,3)*tt8
tt20 = tt1*(H_invDmH(2,2)*tt10-H_invDmH(2,3)*tt7)-tt2*(H_invDmH(2&
&,1)*tt10-H_invDmH(2,3)*tt4)+(H_invDmH(2,1)*tt7-H_invDmH(2,2)*tt4)&
&*tt3
tt21 = tt1*(H_invDmH(2,3)*tt5-H_invDmH(2,2)*tt9)-tt2*(H_invDmH(2,&
&3)*tt6-H_invDmH(2,1)*tt9)+(H_invDmH(2,2)*tt6-H_invDmH(2,1)*tt5)*t&
&t3
tt22 = H_invDmH(3,1)*tt12-H_invDmH(3,2)*tt11+H_invDmH(3,3)*tt8
tt23 = tt1*(H_invDmH(3,2)*tt10-H_invDmH(3,3)*tt7)-tt2*(H_invDmH(3&
&,1)*tt10-H_invDmH(3,3)*tt4)+(H_invDmH(3,1)*tt7-H_invDmH(3,2)*tt4)&
&*tt3
tt24 = tt1*(H_invDmH(3,3)*tt5-H_invDmH(3,2)*tt9)-tt2*(H_invDmH(3,&
&3)*tt6-H_invDmH(3,1)*tt9)+(H_invDmH(3,2)*tt6-H_invDmH(3,1)*tt5)*t&
&t3
tt25 = H_invDmH(4,1)*tt12-H_invDmH(4,2)*tt11+H_invDmH(4,3)*tt8
tt26 = tt1*(H_invDmH(4,2)*tt10-H_invDmH(4,3)*tt7)-tt2*(H_invDmH(4&
&,1)*tt10-H_invDmH(4,3)*tt4)+(H_invDmH(4,1)*tt7-H_invDmH(4,2)*tt4)&
&*tt3
tt27 = tt1*(H_invDmH(4,3)*tt5-H_invDmH(4,2)*tt9)-tt2*(H_invDmH(4,&
&3)*tt6-H_invDmH(4,1)*tt9)+(H_invDmH(4,2)*tt6-H_invDmH(4,1)*tt5)*t&
&t3
tt28 = H_invDmH(5,1)*tt12-H_invDmH(5,2)*tt11+H_invDmH(5,3)*tt8
tt29 = tt1*(H_invDmH(5,2)*tt10-H_invDmH(5,3)*tt7)-tt2*(H_invDmH(5&
&,1)*tt10-H_invDmH(5,3)*tt4)+(H_invDmH(5,1)*tt7-H_invDmH(5,2)*tt4)&
&*tt3
tt30 = tt1*(H_invDmH(5,3)*tt5-H_invDmH(5,2)*tt9)-tt2*(H_invDmH(5,&
&3)*tt6-H_invDmH(5,1)*tt9)+(H_invDmH(5,2)*tt6-H_invDmH(5,1)*tt5)*t&
&t3
tt31 = H_invDmH(6,1)*tt12-H_invDmH(6,2)*tt11+H_invDmH(6,3)*tt8
tt32 = tt1*(H_invDmH(6,2)*tt10-H_invDmH(6,3)*tt7)-tt2*(H_invDmH(6&
&,1)*tt10-H_invDmH(6,3)*tt4)+(H_invDmH(6,1)*tt7-H_invDmH(6,2)*tt4)&
&*tt3
tt33 = tt1*(H_invDmH(6,3)*tt5-H_invDmH(6,2)*tt9)-tt2*(H_invDmH(6,&
&3)*tt6-H_invDmH(6,1)*tt9)+(H_invDmH(6,2)*tt6-H_invDmH(6,1)*tt5)*t&
&t3
tt34 = H_invDmH(7,1)*tt12-H_invDmH(7,2)*tt11+H_invDmH(7,3)*tt8
tt35 = tt1*(H_invDmH(7,2)*tt10-H_invDmH(7,3)*tt7)-tt2*(H_invDmH(7&
&,1)*tt10-H_invDmH(7,3)*tt4)+(H_invDmH(7,1)*tt7-H_invDmH(7,2)*tt4)&
&*tt3
tt36 = tt1*(H_invDmH(7,3)*tt5-H_invDmH(7,2)*tt9)-tt2*(H_invDmH(7,&
&3)*tt6-H_invDmH(7,1)*tt9)+(H_invDmH(7,2)*tt6-H_invDmH(7,1)*tt5)*t&
&t3
tt37 = H_invDmH(8,1)*tt12-H_invDmH(8,2)*tt11+tt8*H_invDmH(8,3)
tt38 = tt1*(H_invDmH(8,2)*tt10-tt7*H_invDmH(8,3))-tt2*(H_invDmH(8&
&,1)*tt10-tt4*H_invDmH(8,3))+(H_invDmH(8,1)*tt7-tt4*H_invDmH(8,2))&
&*tt3
tt39 = tt1*(tt5*H_invDmH(8,3)-H_invDmH(8,2)*tt9)-tt2*(tt6*H_invDm&
&H(8,3)-H_invDmH(8,1)*tt9)+(tt6*H_invDmH(8,2)-H_invDmH(8,1)*tt5)*t&
&t3
jac(1,1) = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt13*tt15*tt16+(mu(1,1)*&
&((-2*tt13*tt15)+2*H_invDmH(1,3)*tt3+2*H_invDmH(1,2)*tt2+2*H_invDm&
&H(1,1)*tt1))/2.0E+0)
jac(1,2) = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt17*tt15*tt16+(mu(1,1)*&
&((-2*tt17*tt15)+2*H_invDmH(1,3)*tt9+2*H_invDmH(1,2)*tt5+2*H_invDm&
&H(1,1)*tt6))/2.0E+0)
jac(1,3) = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt18*tt15*tt16+(mu(1,1)*&
&((-2*tt18*tt15)+2*H_invDmH(1,3)*tt10+2*H_invDmH(1,2)*tt7+2*H_invD&
&mH(1,1)*tt4))/2.0E+0)
jac(1,4) = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt19*tt15*tt16+(mu(1,1)*&
&((-2*tt19*tt15)+2*H_invDmH(2,3)*tt3+2*H_invDmH(2,2)*tt2+2*H_invDm&
&H(2,1)*tt1))/2.0E+0)
jac(1,5) = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt20*tt15*tt16+(mu(1,1)*&
&((-2*tt20*tt15)+2*H_invDmH(2,3)*tt9+2*H_invDmH(2,2)*tt5+2*H_invDm&
&H(2,1)*tt6))/2.0E+0)
jac(1,6) = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt21*tt15*tt16+(mu(1,1)*&
&((-2*tt21*tt15)+2*H_invDmH(2,3)*tt10+2*H_invDmH(2,2)*tt7+2*H_invD&
&mH(2,1)*tt4))/2.0E+0)
jac(1,7) = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt22*tt15*tt16+(mu(1,1)*&
&((-2*tt22*tt15)+2*H_invDmH(3,3)*tt3+2*H_invDmH(3,2)*tt2+2*H_invDm&
&H(3,1)*tt1))/2.0E+0)
jac(1,8) = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt23*tt15*tt16+(mu(1,1)*&
&((-2*tt23*tt15)+2*H_invDmH(3,3)*tt9+2*H_invDmH(3,2)*tt5+2*H_invDm&
&H(3,1)*tt6))/2.0E+0)
jac(1,9) = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt24*tt15*tt16+(mu(1,1)*&
&((-2*tt24*tt15)+2*H_invDmH(3,3)*tt10+2*H_invDmH(3,2)*tt7+2*H_invD&
&mH(3,1)*tt4))/2.0E+0)
jac(1,10) = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt25*tt15*tt16+(mu(1,1)&
&*((-2*tt25*tt15)+2*H_invDmH(4,3)*tt3+2*H_invDmH(4,2)*tt2+2*H_invD&
&mH(4,1)*tt1))/2.0E+0)
jac(1,11) = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt26*tt15*tt16+(mu(1,1)&
&*((-2*tt26*tt15)+2*H_invDmH(4,3)*tt9+2*H_invDmH(4,2)*tt5+2*H_invD&
&mH(4,1)*tt6))/2.0E+0)
jac(1,12) = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt27*tt15*tt16+(mu(1,1)&
&*((-2*tt27*tt15)+2*H_invDmH(4,3)*tt10+2*H_invDmH(4,2)*tt7+2*H_inv&
&DmH(4,1)*tt4))/2.0E+0)
jac(1,13) = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt28*tt15*tt16+(mu(1,1)&
&*((-2*tt28*tt15)+2*H_invDmH(5,3)*tt3+2*H_invDmH(5,2)*tt2+2*H_invD&
&mH(5,1)*tt1))/2.0E+0)
jac(1,14) = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt29*tt15*tt16+(mu(1,1)&
&*((-2*tt29*tt15)+2*H_invDmH(5,3)*tt9+2*H_invDmH(5,2)*tt5+2*H_invD&
&mH(5,1)*tt6))/2.0E+0)
jac(1,15) = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt30*tt15*tt16+(mu(1,1)&
&*((-2*tt30*tt15)+2*H_invDmH(5,3)*tt10+2*H_invDmH(5,2)*tt7+2*H_inv&
&DmH(5,1)*tt4))/2.0E+0)
jac(1,16) = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt31*tt15*tt16+(mu(1,1)&
&*((-2*tt31*tt15)+2*H_invDmH(6,3)*tt3+2*H_invDmH(6,2)*tt2+2*H_invD&
&mH(6,1)*tt1))/2.0E+0)
jac(1,17) = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt32*tt15*tt16+(mu(1,1)&
&*((-2*tt32*tt15)+2*H_invDmH(6,3)*tt9+2*H_invDmH(6,2)*tt5+2*H_invD&
&mH(6,1)*tt6))/2.0E+0)
jac(1,18) = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt33*tt15*tt16+(mu(1,1)&
&*((-2*tt33*tt15)+2*H_invDmH(6,3)*tt10+2*H_invDmH(6,2)*tt7+2*H_inv&
&DmH(6,1)*tt4))/2.0E+0)
jac(1,19) = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt34*tt15*tt16+(mu(1,1)&
&*((-2*tt34*tt15)+2*H_invDmH(7,3)*tt3+2*H_invDmH(7,2)*tt2+2*H_invD&
&mH(7,1)*tt1))/2.0E+0)
jac(1,20) = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt35*tt15*tt16+(mu(1,1)&
&*((-2*tt35*tt15)+2*H_invDmH(7,3)*tt9+2*H_invDmH(7,2)*tt5+2*H_invD&
&mH(7,1)*tt6))/2.0E+0)
jac(1,21) = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt36*tt15*tt16+(mu(1,1)&
&*((-2*tt36*tt15)+2*H_invDmH(7,3)*tt10+2*H_invDmH(7,2)*tt7+2*H_inv&
&DmH(7,1)*tt4))/2.0E+0)
jac(1,22) = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt37*tt15*tt16+(mu(1,1)&
&*((-2*tt37*tt15)+2*H_invDmH(8,3)*tt3+2*H_invDmH(8,2)*tt2+2*H_invD&
&mH(8,1)*tt1))/2.0E+0)
jac(1,23) = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt38*tt15*tt16+(mu(1,1)&
&*((-2*tt38*tt15)+2*H_invDmH(8,3)*tt9+2*H_invDmH(8,2)*tt5+2*H_invD&
&mH(8,1)*tt6))/2.0E+0)
jac(1,24) = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt39*tt15*tt16+(mu(1,1)&
&*((-2*tt39*tt15)+2*H_invDmH(8,3)*tt10+2*H_invDmH(8,2)*tt7+2*H_inv&
&DmH(8,1)*tt4))/2.0E+0)
END 
SUBROUTINE vox_neo_at_quadr_hes(hes, X, H_invDmH, detDmH, gw, mu, lam) 
IMPLICIT NONE 
REAL(KIND=8) hes(24, 24) 
REAL(KIND=8) X(3, 8) 
REAL(KIND=8) H_invDmH(8, 3) 
REAL(KIND=8) detDmH(1, 1) 
REAL(KIND=8) gw(1, 1) 
REAL(KIND=8) mu(1, 1) 
REAL(KIND=8) lam(1, 1) 
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
REAL(KIND=8)  tt73 
REAL(KIND=8)  tt74 
REAL(KIND=8)  tt75 
REAL(KIND=8)  tt76 
REAL(KIND=8)  tt77 
REAL(KIND=8)  tt78 
REAL(KIND=8)  tt79 
REAL(KIND=8)  tt80 
REAL(KIND=8)  tt81 
REAL(KIND=8)  tt82 
REAL(KIND=8)  tt83 
REAL(KIND=8)  tt84 
REAL(KIND=8)  tt85 
REAL(KIND=8)  tt86 
REAL(KIND=8)  tt87 
REAL(KIND=8)  tt88 
REAL(KIND=8)  tt89 
REAL(KIND=8)  tt90 
REAL(KIND=8)  tt91 
REAL(KIND=8)  tt92 
REAL(KIND=8)  tt93 
REAL(KIND=8)  tt94 
REAL(KIND=8)  tt95 
REAL(KIND=8)  tt96 
REAL(KIND=8)  tt97 
REAL(KIND=8)  tt98 
REAL(KIND=8)  tt99 
REAL(KIND=8)  tt100 
REAL(KIND=8)  tt101 
REAL(KIND=8)  tt102 
REAL(KIND=8)  tt103 
REAL(KIND=8)  tt104 
REAL(KIND=8)  tt105 
REAL(KIND=8)  tt106 
REAL(KIND=8)  tt107 
REAL(KIND=8)  tt108 
REAL(KIND=8)  tt109 
REAL(KIND=8)  tt110 
REAL(KIND=8)  tt111 
REAL(KIND=8)  tt112 
REAL(KIND=8)  tt113 
REAL(KIND=8)  tt114 
REAL(KIND=8)  tt115 
REAL(KIND=8)  tt116 
REAL(KIND=8)  tt117 
REAL(KIND=8)  tt118 
REAL(KIND=8)  tt119 
REAL(KIND=8)  tt120 
REAL(KIND=8)  tt121 
REAL(KIND=8)  tt122 
REAL(KIND=8)  tt123 
REAL(KIND=8)  tt124 
REAL(KIND=8)  tt125 
REAL(KIND=8)  tt126 
REAL(KIND=8)  tt127 
REAL(KIND=8)  tt128 
REAL(KIND=8)  tt129 
REAL(KIND=8)  tt130 
REAL(KIND=8)  tt131 
REAL(KIND=8)  tt132 
REAL(KIND=8)  tt133 
REAL(KIND=8)  tt134 
REAL(KIND=8)  tt135 
REAL(KIND=8)  tt136 
REAL(KIND=8)  tt137 
REAL(KIND=8)  tt138 
REAL(KIND=8)  tt139 
REAL(KIND=8)  tt140 
REAL(KIND=8)  tt141 
REAL(KIND=8)  tt142 
REAL(KIND=8)  tt143 
REAL(KIND=8)  tt144 
REAL(KIND=8)  tt145 
REAL(KIND=8)  tt146 
REAL(KIND=8)  tt147 
REAL(KIND=8)  tt148 
REAL(KIND=8)  tt149 
REAL(KIND=8)  tt150 
REAL(KIND=8)  tt151 
REAL(KIND=8)  tt152 
REAL(KIND=8)  tt153 
REAL(KIND=8)  tt154 
REAL(KIND=8)  tt155 
REAL(KIND=8)  tt156 
REAL(KIND=8)  tt157 
REAL(KIND=8)  tt158 
REAL(KIND=8)  tt159 
REAL(KIND=8)  tt160 
REAL(KIND=8)  tt161 
REAL(KIND=8)  tt162 
REAL(KIND=8)  tt163 
REAL(KIND=8)  tt164 
REAL(KIND=8)  tt165 
REAL(KIND=8)  tt166 
REAL(KIND=8)  tt167 
REAL(KIND=8)  tt168 
REAL(KIND=8)  tt169 
REAL(KIND=8)  tt170 
REAL(KIND=8)  tt171 
REAL(KIND=8)  tt172 
REAL(KIND=8)  tt173 
REAL(KIND=8)  tt174 
REAL(KIND=8)  tt175 
REAL(KIND=8)  tt176 
REAL(KIND=8)  tt177 
REAL(KIND=8)  tt178 
REAL(KIND=8)  tt179 
REAL(KIND=8)  tt180 
REAL(KIND=8)  tt181 
REAL(KIND=8)  tt182 
REAL(KIND=8)  tt183 
REAL(KIND=8)  tt184 
REAL(KIND=8)  tt185 
REAL(KIND=8)  tt186 
REAL(KIND=8)  tt187 
REAL(KIND=8)  tt188 
REAL(KIND=8)  tt189 
REAL(KIND=8)  tt190 
REAL(KIND=8)  tt191 
REAL(KIND=8)  tt192 
REAL(KIND=8)  tt193 
REAL(KIND=8)  tt194 
REAL(KIND=8)  tt195 
REAL(KIND=8)  tt196 
REAL(KIND=8)  tt197 
REAL(KIND=8)  tt198 
REAL(KIND=8)  tt199 
REAL(KIND=8)  tt200 
REAL(KIND=8)  tt201 
REAL(KIND=8)  tt202 
REAL(KIND=8)  tt203 
REAL(KIND=8)  tt204 
REAL(KIND=8)  tt205 
REAL(KIND=8)  tt206 
REAL(KIND=8)  tt207 
REAL(KIND=8)  tt208 
REAL(KIND=8)  tt209 
REAL(KIND=8)  tt210 
REAL(KIND=8)  tt211 
REAL(KIND=8)  tt212 
REAL(KIND=8)  tt213 
REAL(KIND=8)  tt214 
REAL(KIND=8)  tt215 
REAL(KIND=8)  tt216 
REAL(KIND=8)  tt217 
REAL(KIND=8)  tt218 
REAL(KIND=8)  tt219 
REAL(KIND=8)  tt220 
REAL(KIND=8)  tt221 
REAL(KIND=8)  tt222 
REAL(KIND=8)  tt223 
REAL(KIND=8)  tt224 
REAL(KIND=8)  tt225 
REAL(KIND=8)  tt226 
REAL(KIND=8)  tt227 
REAL(KIND=8)  tt228 
REAL(KIND=8)  tt229 
REAL(KIND=8)  tt230 
REAL(KIND=8)  tt231 
REAL(KIND=8)  tt232 
REAL(KIND=8)  tt233 
REAL(KIND=8)  tt234 
REAL(KIND=8)  tt235 
REAL(KIND=8)  tt236 
REAL(KIND=8)  tt237 
REAL(KIND=8)  tt238 
REAL(KIND=8)  tt239 
REAL(KIND=8)  tt240 
REAL(KIND=8)  tt241 
REAL(KIND=8)  tt242 
REAL(KIND=8)  tt243 
REAL(KIND=8)  tt244 
REAL(KIND=8)  tt245 
REAL(KIND=8)  tt246 
REAL(KIND=8)  tt247 
REAL(KIND=8)  tt248 
REAL(KIND=8)  tt249 
REAL(KIND=8)  tt250 
REAL(KIND=8)  tt251 
REAL(KIND=8)  tt252 
REAL(KIND=8)  tt253 
REAL(KIND=8)  tt254 
REAL(KIND=8)  tt255 
REAL(KIND=8)  tt256 
REAL(KIND=8)  tt257 
REAL(KIND=8)  tt258 
REAL(KIND=8)  tt259 
REAL(KIND=8)  tt260 
REAL(KIND=8)  tt261 
REAL(KIND=8)  tt262 
REAL(KIND=8)  tt263 
REAL(KIND=8)  tt264 
REAL(KIND=8)  tt265 
REAL(KIND=8)  tt266 
REAL(KIND=8)  tt267 
REAL(KIND=8)  tt268 
REAL(KIND=8)  tt269 
REAL(KIND=8)  tt270 
REAL(KIND=8)  tt271 
REAL(KIND=8)  tt272 
REAL(KIND=8)  tt273 
REAL(KIND=8)  tt274 
REAL(KIND=8)  tt275 
REAL(KIND=8)  tt276 
REAL(KIND=8)  tt277 
REAL(KIND=8)  tt278 
REAL(KIND=8)  tt279 
REAL(KIND=8)  tt280 
REAL(KIND=8)  tt281 
REAL(KIND=8)  tt282 
REAL(KIND=8)  tt283 
REAL(KIND=8)  tt284 
REAL(KIND=8)  tt285 
REAL(KIND=8)  tt286 
REAL(KIND=8)  tt287 
REAL(KIND=8)  tt288 
REAL(KIND=8)  tt289 
REAL(KIND=8)  tt290 
REAL(KIND=8)  tt291 
REAL(KIND=8)  tt292 
REAL(KIND=8)  tt293 
REAL(KIND=8)  tt294 
REAL(KIND=8)  tt295 
REAL(KIND=8)  tt296 
REAL(KIND=8)  tt297 
REAL(KIND=8)  tt298 
REAL(KIND=8)  tt299 
REAL(KIND=8)  tt300 
REAL(KIND=8)  tt301 
REAL(KIND=8)  tt302 
REAL(KIND=8)  tt303 
REAL(KIND=8)  tt304 
REAL(KIND=8)  tt305 
REAL(KIND=8)  tt306 
REAL(KIND=8)  tt307 
REAL(KIND=8)  tt308 
REAL(KIND=8)  tt309 
REAL(KIND=8)  tt310 
REAL(KIND=8)  tt311 
REAL(KIND=8)  tt312 
REAL(KIND=8)  tt313 
REAL(KIND=8)  tt314 
REAL(KIND=8)  tt315 
REAL(KIND=8)  tt316 
REAL(KIND=8)  tt317 
REAL(KIND=8)  tt318 
REAL(KIND=8)  tt319 
REAL(KIND=8)  tt320 
REAL(KIND=8)  tt321 
REAL(KIND=8)  tt322 
REAL(KIND=8)  tt323 
REAL(KIND=8)  tt324 
REAL(KIND=8)  tt325 
REAL(KIND=8)  tt326 
REAL(KIND=8)  tt327 
REAL(KIND=8)  tt328 
REAL(KIND=8)  tt329 
REAL(KIND=8)  tt330 
REAL(KIND=8)  tt331 
REAL(KIND=8)  tt332 
REAL(KIND=8)  tt333 
REAL(KIND=8)  tt334 
REAL(KIND=8)  tt335 
REAL(KIND=8)  tt336 
REAL(KIND=8)  tt337 
REAL(KIND=8)  tt338 
REAL(KIND=8)  tt339 
REAL(KIND=8)  tt340 
REAL(KIND=8)  tt341 
REAL(KIND=8)  tt342 
REAL(KIND=8)  tt343 
REAL(KIND=8)  tt344 
REAL(KIND=8)  tt345 
REAL(KIND=8)  tt346 
REAL(KIND=8)  tt347 
REAL(KIND=8)  tt348 
REAL(KIND=8)  tt349 
REAL(KIND=8)  tt350 
REAL(KIND=8)  tt351 
REAL(KIND=8)  tt352 
REAL(KIND=8)  tt353 
REAL(KIND=8)  tt354 
REAL(KIND=8)  tt355 
REAL(KIND=8)  tt356 
REAL(KIND=8)  tt357 
REAL(KIND=8)  tt358 
REAL(KIND=8)  tt359 
REAL(KIND=8)  tt360 
REAL(KIND=8)  tt361 
REAL(KIND=8)  tt362 
REAL(KIND=8)  tt363 
REAL(KIND=8)  tt364 
REAL(KIND=8)  tt365 
REAL(KIND=8)  tt366 
REAL(KIND=8)  tt367 
REAL(KIND=8)  tt368 
REAL(KIND=8)  tt369 
REAL(KIND=8)  tt370 
REAL(KIND=8)  tt371 
REAL(KIND=8)  tt372 
REAL(KIND=8)  tt373 
REAL(KIND=8)  tt374 
REAL(KIND=8)  tt375 
REAL(KIND=8)  tt376 
REAL(KIND=8)  tt377 
REAL(KIND=8)  tt378 
REAL(KIND=8)  tt379 
REAL(KIND=8)  tt380 
REAL(KIND=8)  tt381 
REAL(KIND=8)  tt382 
REAL(KIND=8)  tt383 
REAL(KIND=8)  tt384 
REAL(KIND=8)  tt385 
REAL(KIND=8)  tt386 
REAL(KIND=8)  tt387 
REAL(KIND=8)  tt388 
REAL(KIND=8)  tt389 
REAL(KIND=8)  tt390 
REAL(KIND=8)  tt391 
REAL(KIND=8)  tt392 
REAL(KIND=8)  tt393 
REAL(KIND=8)  tt394 
REAL(KIND=8)  tt395 
REAL(KIND=8)  tt396 
REAL(KIND=8)  tt397 
REAL(KIND=8)  tt398 
REAL(KIND=8)  tt399 
REAL(KIND=8)  tt400 
REAL(KIND=8)  tt401 
REAL(KIND=8)  tt402 
REAL(KIND=8)  tt403 
REAL(KIND=8)  tt404 
REAL(KIND=8)  tt405 
REAL(KIND=8)  tt406 
REAL(KIND=8)  tt407 
REAL(KIND=8)  tt408 
REAL(KIND=8)  tt409 
REAL(KIND=8)  tt410 
REAL(KIND=8)  tt411 
REAL(KIND=8)  tt412 
REAL(KIND=8)  tt413 
REAL(KIND=8)  tt414 
REAL(KIND=8)  tt415 
REAL(KIND=8)  tt416 
REAL(KIND=8)  tt417 
REAL(KIND=8)  tt418 
REAL(KIND=8)  tt419 
REAL(KIND=8)  tt420 
REAL(KIND=8)  tt421 
REAL(KIND=8)  tt422 
REAL(KIND=8)  tt423 
REAL(KIND=8)  tt424 
REAL(KIND=8)  tt425 
REAL(KIND=8)  tt426 
REAL(KIND=8)  tt427 
REAL(KIND=8)  tt428 
REAL(KIND=8)  tt429 
REAL(KIND=8)  tt430 
REAL(KIND=8)  tt431 
REAL(KIND=8)  tt432 
REAL(KIND=8)  tt433 
REAL(KIND=8)  tt434 
REAL(KIND=8)  tt435 
REAL(KIND=8)  tt436 
REAL(KIND=8)  tt437 
REAL(KIND=8)  tt438 
REAL(KIND=8)  tt439 
REAL(KIND=8)  tt440 
REAL(KIND=8)  tt441 
REAL(KIND=8)  tt442 
REAL(KIND=8)  tt443 
REAL(KIND=8)  tt444 
REAL(KIND=8)  tt445 
REAL(KIND=8)  tt446 
REAL(KIND=8)  tt447 
REAL(KIND=8)  tt448 
REAL(KIND=8)  tt449 
REAL(KIND=8)  tt450 
REAL(KIND=8)  tt451 
REAL(KIND=8)  tt452 
REAL(KIND=8)  tt453 
REAL(KIND=8)  tt454 
REAL(KIND=8)  tt455 
REAL(KIND=8)  tt456 
REAL(KIND=8)  tt457 
REAL(KIND=8)  tt458 
REAL(KIND=8)  tt459 
REAL(KIND=8)  tt460 
REAL(KIND=8)  tt461 
REAL(KIND=8)  tt462 
REAL(KIND=8)  tt463 
REAL(KIND=8)  tt464 
REAL(KIND=8)  tt465 
REAL(KIND=8)  tt466 
REAL(KIND=8)  tt467 
REAL(KIND=8)  tt468 
REAL(KIND=8)  tt469 
REAL(KIND=8)  tt470 
REAL(KIND=8)  tt471 
REAL(KIND=8)  tt472 
REAL(KIND=8)  tt473 
REAL(KIND=8)  tt474 
REAL(KIND=8)  tt475 
REAL(KIND=8)  tt476 
REAL(KIND=8)  tt477 
REAL(KIND=8)  tt478 
REAL(KIND=8)  tt479 
REAL(KIND=8)  tt480 
REAL(KIND=8)  tt481 
REAL(KIND=8)  tt482 
REAL(KIND=8)  tt483 
REAL(KIND=8)  tt484 
REAL(KIND=8)  tt485 
REAL(KIND=8)  tt486 
REAL(KIND=8)  tt487 
REAL(KIND=8)  tt488 
REAL(KIND=8)  tt489 
REAL(KIND=8)  tt490 
REAL(KIND=8)  tt491 
REAL(KIND=8)  tt492 
REAL(KIND=8)  tt493 
REAL(KIND=8)  tt494 
REAL(KIND=8)  tt495 
REAL(KIND=8)  tt496 
REAL(KIND=8)  tt497 
REAL(KIND=8)  tt498 
REAL(KIND=8)  tt499 
REAL(KIND=8)  tt500 
REAL(KIND=8)  tt501 
REAL(KIND=8)  tt502 
REAL(KIND=8)  tt503 
REAL(KIND=8)  tt504 
REAL(KIND=8)  tt505 
REAL(KIND=8)  tt506 
REAL(KIND=8)  tt507 
REAL(KIND=8)  tt508 
REAL(KIND=8)  tt509 
REAL(KIND=8)  tt510 
REAL(KIND=8)  tt511 
REAL(KIND=8)  tt512 
REAL(KIND=8)  tt513 
REAL(KIND=8)  tt514 
REAL(KIND=8)  tt515 
REAL(KIND=8)  tt516 
REAL(KIND=8)  tt517 
REAL(KIND=8)  tt518 
REAL(KIND=8)  tt519 
REAL(KIND=8)  tt520 
REAL(KIND=8)  tt521 
REAL(KIND=8)  tt522 
REAL(KIND=8)  tt523 
REAL(KIND=8)  tt524 
REAL(KIND=8)  tt525 
REAL(KIND=8)  tt526 
REAL(KIND=8)  tt527 
REAL(KIND=8)  tt528 
REAL(KIND=8)  tt529 
REAL(KIND=8)  tt530 
REAL(KIND=8)  tt531 
REAL(KIND=8)  tt532 
REAL(KIND=8)  tt533 
REAL(KIND=8)  tt534 
REAL(KIND=8)  tt535 
REAL(KIND=8)  tt536 
REAL(KIND=8)  tt537 
REAL(KIND=8)  tt538 
REAL(KIND=8)  tt539 
REAL(KIND=8)  tt540 
REAL(KIND=8)  tt541 
REAL(KIND=8)  tt542 
REAL(KIND=8)  tt543 
REAL(KIND=8)  tt544 
REAL(KIND=8)  tt545 
REAL(KIND=8)  tt546 
REAL(KIND=8)  tt547 
REAL(KIND=8)  tt548 
REAL(KIND=8)  tt549 
REAL(KIND=8)  tt550 
REAL(KIND=8)  tt551 
REAL(KIND=8)  tt552 
REAL(KIND=8)  tt553 
REAL(KIND=8)  tt554 
REAL(KIND=8)  tt555 
REAL(KIND=8)  tt556 
REAL(KIND=8)  tt557 
REAL(KIND=8)  tt558 
REAL(KIND=8)  tt559 
REAL(KIND=8)  tt560 
REAL(KIND=8)  tt561 
REAL(KIND=8)  tt562 
REAL(KIND=8)  tt563 
REAL(KIND=8)  tt564 
REAL(KIND=8)  tt565 
REAL(KIND=8)  tt566 
REAL(KIND=8)  tt567 
REAL(KIND=8)  tt568 
REAL(KIND=8)  tt569 
REAL(KIND=8)  tt570 
REAL(KIND=8)  tt571 
REAL(KIND=8)  tt572 
REAL(KIND=8)  tt573 
REAL(KIND=8)  tt574 
REAL(KIND=8)  tt575 
REAL(KIND=8)  tt576 
REAL(KIND=8)  tt577 
REAL(KIND=8)  tt578 
REAL(KIND=8)  tt579 
REAL(KIND=8)  tt580 
REAL(KIND=8)  tt581 
REAL(KIND=8)  tt582 
REAL(KIND=8)  tt583 
REAL(KIND=8)  tt584 
REAL(KIND=8)  tt585 
REAL(KIND=8)  tt586 
REAL(KIND=8)  tt587 
REAL(KIND=8)  tt588 
REAL(KIND=8)  tt589 
REAL(KIND=8)  tt590 
REAL(KIND=8)  tt591 
REAL(KIND=8)  tt592 
REAL(KIND=8)  tt593 
REAL(KIND=8)  tt594 
REAL(KIND=8)  tt595 
REAL(KIND=8)  tt596 
REAL(KIND=8)  tt597 
REAL(KIND=8)  tt598 
REAL(KIND=8)  tt599 
REAL(KIND=8)  tt600 
REAL(KIND=8)  tt601 
REAL(KIND=8)  tt602 
REAL(KIND=8)  tt603 
REAL(KIND=8)  tt604 
REAL(KIND=8)  tt605 
REAL(KIND=8)  tt606 
REAL(KIND=8)  tt607 
REAL(KIND=8)  tt608 
REAL(KIND=8)  tt609 
REAL(KIND=8)  tt610 
REAL(KIND=8)  tt611 
REAL(KIND=8)  tt612 
REAL(KIND=8)  tt613 
REAL(KIND=8)  tt614 
REAL(KIND=8)  tt615 
REAL(KIND=8)  tt616 
REAL(KIND=8)  tt617 
REAL(KIND=8)  tt618 
REAL(KIND=8)  tt619 
REAL(KIND=8)  tt620 
REAL(KIND=8)  tt621 
REAL(KIND=8)  tt622 
REAL(KIND=8)  tt623 
REAL(KIND=8)  tt624 
REAL(KIND=8)  tt625 
REAL(KIND=8)  tt626 
REAL(KIND=8)  tt627 
REAL(KIND=8)  tt628 
REAL(KIND=8)  tt629 
REAL(KIND=8)  tt630 
REAL(KIND=8)  tt631 
REAL(KIND=8)  tt632 
REAL(KIND=8)  tt633 
REAL(KIND=8)  tt634 
REAL(KIND=8)  tt635 
REAL(KIND=8)  tt636 
REAL(KIND=8)  tt637 
REAL(KIND=8)  tt638 
REAL(KIND=8)  tt639 
REAL(KIND=8)  tt640 
REAL(KIND=8)  tt641 
REAL(KIND=8)  tt642 
REAL(KIND=8)  tt643 
REAL(KIND=8)  tt644 
REAL(KIND=8)  tt645 
REAL(KIND=8)  tt646 
REAL(KIND=8)  tt647 
REAL(KIND=8)  tt648 
REAL(KIND=8)  tt649 
REAL(KIND=8)  tt650 
REAL(KIND=8)  tt651 
REAL(KIND=8)  tt652 
REAL(KIND=8)  tt653 
REAL(KIND=8)  tt654 
REAL(KIND=8)  tt655 
REAL(KIND=8)  tt656 
REAL(KIND=8)  tt657 
REAL(KIND=8)  tt658 
REAL(KIND=8)  tt659 
REAL(KIND=8)  tt660 
REAL(KIND=8)  tt661 
REAL(KIND=8)  tt662 
REAL(KIND=8)  tt663 
REAL(KIND=8)  tt664 
REAL(KIND=8)  tt665 
REAL(KIND=8)  tt666 
REAL(KIND=8)  tt667 
REAL(KIND=8)  tt668 
REAL(KIND=8)  tt669 
REAL(KIND=8)  tt670 
REAL(KIND=8)  tt671 
REAL(KIND=8)  tt672 
REAL(KIND=8)  tt673 
REAL(KIND=8)  tt674 
REAL(KIND=8)  tt675 
REAL(KIND=8)  tt676 
REAL(KIND=8)  tt677 
REAL(KIND=8)  tt678 
REAL(KIND=8)  tt679 
REAL(KIND=8)  tt680 
tt1 = X(3,8)*H_invDmH(8,1)+X(3,7)*H_invDmH(7,1)+X(3,6)*H_invDmH(6&
&,1)+X(3,5)*H_invDmH(5,1)+X(3,4)*H_invDmH(4,1)+H_invDmH(3,1)*X(3,3&
&)+H_invDmH(2,1)*X(3,2)+H_invDmH(1,1)*X(3,1)
tt2 = X(2,8)*H_invDmH(8,2)+X(2,7)*H_invDmH(7,2)+X(2,6)*H_invDmH(6&
&,2)+X(2,5)*H_invDmH(5,2)+X(2,4)*H_invDmH(4,2)+X(2,3)*H_invDmH(3,2&
&)+H_invDmH(2,2)*X(2,2)+H_invDmH(1,2)*X(2,1)
tt3 = X(2,8)*H_invDmH(8,1)+X(2,7)*H_invDmH(7,1)+X(2,6)*H_invDmH(6&
&,1)+X(2,5)*H_invDmH(5,1)+X(2,4)*H_invDmH(4,1)+X(2,3)*H_invDmH(3,1&
&)+H_invDmH(2,1)*X(2,2)+H_invDmH(1,1)*X(2,1)
tt4 = X(3,8)*H_invDmH(8,2)+X(3,7)*H_invDmH(7,2)+X(3,6)*H_invDmH(6&
&,2)+X(3,5)*H_invDmH(5,2)+X(3,4)*H_invDmH(4,2)+H_invDmH(3,2)*X(3,3&
&)+H_invDmH(2,2)*X(3,2)+H_invDmH(1,2)*X(3,1)
tt5 = tt3*tt4-tt1*tt2
tt6 = X(2,8)*H_invDmH(8,3)+X(2,7)*H_invDmH(7,3)+X(2,6)*H_invDmH(6&
&,3)+X(2,5)*H_invDmH(5,3)+X(2,4)*H_invDmH(4,3)+X(2,3)*H_invDmH(3,3&
&)+X(2,2)*H_invDmH(2,3)+H_invDmH(1,3)*X(2,1)
tt7 = X(3,8)*H_invDmH(8,3)+X(3,7)*H_invDmH(7,3)+X(3,6)*H_invDmH(6&
&,3)+X(3,5)*H_invDmH(5,3)+X(3,4)*H_invDmH(4,3)+H_invDmH(3,3)*X(3,3&
&)+H_invDmH(2,3)*X(3,2)+H_invDmH(1,3)*X(3,1)
tt8 = tt3*tt7-tt1*tt6
tt9 = tt2*tt7-tt4*tt6
tt10 = H_invDmH(1,1)*tt9-H_invDmH(1,2)*tt8+H_invDmH(1,3)*tt5
tt11 = tt10**2
tt12 = X(1,8)*H_invDmH(8,3)+X(1,7)*H_invDmH(7,3)+X(1,6)*H_invDmH(&
&6,3)+X(1,5)*H_invDmH(5,3)+X(1,4)*H_invDmH(4,3)+X(1,3)*H_invDmH(3,&
&3)+X(1,2)*H_invDmH(2,3)+X(1,1)*H_invDmH(1,3)
tt13 = X(1,8)*H_invDmH(8,2)+X(1,7)*H_invDmH(7,2)+X(1,6)*H_invDmH(&
&6,2)+X(1,5)*H_invDmH(5,2)+X(1,4)*H_invDmH(4,2)+X(1,3)*H_invDmH(3,&
&2)+X(1,2)*H_invDmH(2,2)+X(1,1)*H_invDmH(1,2)
tt14 = X(1,8)*H_invDmH(8,1)+X(1,7)*H_invDmH(7,1)+X(1,6)*H_invDmH(&
&6,1)+X(1,5)*H_invDmH(5,1)+X(1,4)*H_invDmH(4,1)+X(1,3)*H_invDmH(3,&
&1)+X(1,2)*H_invDmH(2,1)+H_invDmH(1,1)*X(1,1)
tt15 = tt14*tt9-tt13*tt8+tt5*tt12
tt16 = 1/tt15**2
tt17 = 2*H_invDmH(1,1)**2
tt18 = 2*H_invDmH(1,2)**2
tt19 = 2*H_invDmH(1,3)**2
tt20 = log(tt15)
tt21 = H_invDmH(1,1)*tt4-H_invDmH(1,2)*tt1
tt22 = H_invDmH(1,1)*tt7-H_invDmH(1,3)*tt1
tt23 = H_invDmH(1,2)*tt7-H_invDmH(1,3)*tt4
tt24 = tt14*tt23-tt13*tt22+tt21*tt12
tt25 = H_invDmH(1,1)*tt23-H_invDmH(1,2)*tt22+H_invDmH(1,3)*tt21
tt26 = 1/tt15
tt27 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt25*tt26*tt20-lam(1,1)*tt24&
&*tt10*tt16*tt20+(mu(1,1)*(2*tt24*tt10*tt16-2*tt25*tt26))/2.0E+0+l&
&am(1,1)*tt24*tt10*tt16)
tt28 = H_invDmH(1,2)*tt3-H_invDmH(1,1)*tt2
tt29 = H_invDmH(1,3)*tt3-H_invDmH(1,1)*tt6
tt30 = H_invDmH(1,3)*tt2-H_invDmH(1,2)*tt6
tt31 = tt14*tt30-tt13*tt29+tt28*tt12
tt32 = H_invDmH(1,1)*tt30-H_invDmH(1,2)*tt29+H_invDmH(1,3)*tt28
tt33 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt32*tt26*tt20-lam(1,1)*tt31&
&*tt10*tt16*tt20+(mu(1,1)*(2*tt31*tt10*tt16-2*tt32*tt26))/2.0E+0+l&
&am(1,1)*tt31*tt10*tt16)
tt34 = H_invDmH(2,1)*tt9-H_invDmH(2,2)*tt8+H_invDmH(2,3)*tt5
tt35 = 2*H_invDmH(1,1)*H_invDmH(2,1)
tt36 = 2*H_invDmH(1,2)*H_invDmH(2,2)
tt37 = 2*H_invDmH(1,3)*H_invDmH(2,3)
tt38 = detDmH(1,1)*gw(1,1)*((-lam(1,1)*tt10*tt34*tt16*tt20)+(mu(1&
&,1)*(2*tt10*tt34*tt16+tt37+tt36+tt35))/2.0E+0+lam(1,1)*tt10*tt34*&
&tt16)
tt39 = H_invDmH(2,1)*tt4-H_invDmH(2,2)*tt1
tt40 = H_invDmH(2,1)*tt7-H_invDmH(2,3)*tt1
tt41 = H_invDmH(2,2)*tt7-H_invDmH(2,3)*tt4
tt42 = tt14*tt41-tt13*tt40+tt39*tt12
tt43 = H_invDmH(1,1)*tt41-H_invDmH(1,2)*tt40+H_invDmH(1,3)*tt39
tt44 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt43*tt26*tt20-lam(1,1)*tt42&
&*tt10*tt16*tt20+(mu(1,1)*(2*tt42*tt10*tt16-2*tt43*tt26))/2.0E+0+l&
&am(1,1)*tt42*tt10*tt16)
tt45 = H_invDmH(2,2)*tt3-H_invDmH(2,1)*tt2
tt46 = H_invDmH(2,3)*tt3-H_invDmH(2,1)*tt6
tt47 = H_invDmH(2,3)*tt2-H_invDmH(2,2)*tt6
tt48 = tt14*tt47-tt13*tt46+tt45*tt12
tt49 = H_invDmH(1,1)*tt47-H_invDmH(1,2)*tt46+H_invDmH(1,3)*tt45
tt50 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt49*tt26*tt20-lam(1,1)*tt48&
&*tt10*tt16*tt20+(mu(1,1)*(2*tt48*tt10*tt16-2*tt49*tt26))/2.0E+0+l&
&am(1,1)*tt48*tt10*tt16)
tt51 = H_invDmH(3,1)*tt9-H_invDmH(3,2)*tt8+H_invDmH(3,3)*tt5
tt52 = 2*H_invDmH(1,1)*H_invDmH(3,1)
tt53 = 2*H_invDmH(1,2)*H_invDmH(3,2)
tt54 = 2*H_invDmH(1,3)*H_invDmH(3,3)
tt55 = detDmH(1,1)*gw(1,1)*((-lam(1,1)*tt10*tt51*tt16*tt20)+(mu(1&
&,1)*(2*tt10*tt51*tt16+tt54+tt53+tt52))/2.0E+0+lam(1,1)*tt10*tt51*&
&tt16)
tt56 = H_invDmH(3,1)*tt4-H_invDmH(3,2)*tt1
tt57 = H_invDmH(3,1)*tt7-H_invDmH(3,3)*tt1
tt58 = H_invDmH(3,2)*tt7-H_invDmH(3,3)*tt4
tt59 = tt14*tt58-tt13*tt57+tt56*tt12
tt60 = H_invDmH(1,1)*tt58-H_invDmH(1,2)*tt57+H_invDmH(1,3)*tt56
tt61 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt60*tt26*tt20-lam(1,1)*tt59&
&*tt10*tt16*tt20+(mu(1,1)*(2*tt59*tt10*tt16-2*tt60*tt26))/2.0E+0+l&
&am(1,1)*tt59*tt10*tt16)
tt62 = H_invDmH(3,2)*tt3-H_invDmH(3,1)*tt2
tt63 = H_invDmH(3,3)*tt3-H_invDmH(3,1)*tt6
tt64 = H_invDmH(3,3)*tt2-H_invDmH(3,2)*tt6
tt65 = tt14*tt64-tt13*tt63+tt62*tt12
tt66 = H_invDmH(1,1)*tt64-H_invDmH(1,2)*tt63+H_invDmH(1,3)*tt62
tt67 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt66*tt26*tt20-lam(1,1)*tt65&
&*tt10*tt16*tt20+(mu(1,1)*(2*tt65*tt10*tt16-2*tt66*tt26))/2.0E+0+l&
&am(1,1)*tt65*tt10*tt16)
tt68 = H_invDmH(4,1)*tt9-H_invDmH(4,2)*tt8+H_invDmH(4,3)*tt5
tt69 = 2*H_invDmH(1,1)*H_invDmH(4,1)
tt70 = 2*H_invDmH(1,2)*H_invDmH(4,2)
tt71 = 2*H_invDmH(1,3)*H_invDmH(4,3)
tt72 = detDmH(1,1)*gw(1,1)*((-lam(1,1)*tt10*tt68*tt16*tt20)+(mu(1&
&,1)*(2*tt10*tt68*tt16+tt71+tt70+tt69))/2.0E+0+lam(1,1)*tt10*tt68*&
&tt16)
tt73 = H_invDmH(4,1)*tt4-H_invDmH(4,2)*tt1
tt74 = H_invDmH(4,1)*tt7-H_invDmH(4,3)*tt1
tt75 = H_invDmH(4,2)*tt7-H_invDmH(4,3)*tt4
tt76 = tt14*tt75-tt13*tt74+tt73*tt12
tt77 = H_invDmH(1,1)*tt75-H_invDmH(1,2)*tt74+H_invDmH(1,3)*tt73
tt78 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt77*tt26*tt20-lam(1,1)*tt76&
&*tt10*tt16*tt20+(mu(1,1)*(2*tt76*tt10*tt16-2*tt77*tt26))/2.0E+0+l&
&am(1,1)*tt76*tt10*tt16)
tt79 = H_invDmH(4,2)*tt3-H_invDmH(4,1)*tt2
tt80 = H_invDmH(4,3)*tt3-H_invDmH(4,1)*tt6
tt81 = H_invDmH(4,3)*tt2-H_invDmH(4,2)*tt6
tt82 = tt14*tt81-tt13*tt80+tt79*tt12
tt83 = H_invDmH(1,1)*tt81-H_invDmH(1,2)*tt80+H_invDmH(1,3)*tt79
tt84 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt83*tt26*tt20-lam(1,1)*tt82&
&*tt10*tt16*tt20+(mu(1,1)*(2*tt82*tt10*tt16-2*tt83*tt26))/2.0E+0+l&
&am(1,1)*tt82*tt10*tt16)
tt85 = H_invDmH(5,1)*tt9-H_invDmH(5,2)*tt8+H_invDmH(5,3)*tt5
tt86 = 2*H_invDmH(1,1)*H_invDmH(5,1)
tt87 = 2*H_invDmH(1,2)*H_invDmH(5,2)
tt88 = 2*H_invDmH(1,3)*H_invDmH(5,3)
tt89 = detDmH(1,1)*gw(1,1)*((-lam(1,1)*tt10*tt85*tt16*tt20)+(mu(1&
&,1)*(2*tt10*tt85*tt16+tt88+tt87+tt86))/2.0E+0+lam(1,1)*tt10*tt85*&
&tt16)
tt90 = H_invDmH(5,1)*tt4-H_invDmH(5,2)*tt1
tt91 = H_invDmH(5,1)*tt7-H_invDmH(5,3)*tt1
tt92 = H_invDmH(5,2)*tt7-H_invDmH(5,3)*tt4
tt93 = tt14*tt92-tt13*tt91+tt90*tt12
tt94 = H_invDmH(1,1)*tt92-H_invDmH(1,2)*tt91+H_invDmH(1,3)*tt90
tt95 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt94*tt26*tt20-lam(1,1)*tt93&
&*tt10*tt16*tt20+(mu(1,1)*(2*tt93*tt10*tt16-2*tt94*tt26))/2.0E+0+l&
&am(1,1)*tt93*tt10*tt16)
tt96 = H_invDmH(5,2)*tt3-H_invDmH(5,1)*tt2
tt97 = H_invDmH(5,3)*tt3-H_invDmH(5,1)*tt6
tt98 = H_invDmH(5,3)*tt2-H_invDmH(5,2)*tt6
tt99 = tt14*tt98-tt13*tt97+tt96*tt12
tt100 = H_invDmH(1,1)*tt98-H_invDmH(1,2)*tt97+H_invDmH(1,3)*tt96
tt101 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt100*tt26*tt20-lam(1,1)*tt&
&99*tt10*tt16*tt20+(mu(1,1)*(2*tt99*tt10*tt16-2*tt100*tt26))/2.0E+&
&0+lam(1,1)*tt99*tt10*tt16)
tt102 = H_invDmH(6,1)*tt9-H_invDmH(6,2)*tt8+H_invDmH(6,3)*tt5
tt103 = 2*H_invDmH(1,1)*H_invDmH(6,1)
tt104 = 2*H_invDmH(1,2)*H_invDmH(6,2)
tt105 = 2*H_invDmH(1,3)*H_invDmH(6,3)
tt106 = detDmH(1,1)*gw(1,1)*((-lam(1,1)*tt10*tt102*tt16*tt20)+(mu&
&(1,1)*(2*tt10*tt102*tt16+tt105+tt104+tt103))/2.0E+0+lam(1,1)*tt10&
&*tt102*tt16)
tt107 = H_invDmH(6,1)*tt4-H_invDmH(6,2)*tt1
tt108 = H_invDmH(6,1)*tt7-H_invDmH(6,3)*tt1
tt109 = H_invDmH(6,2)*tt7-H_invDmH(6,3)*tt4
tt110 = tt14*tt109-tt13*tt108+tt107*tt12
tt111 = H_invDmH(1,1)*tt109-H_invDmH(1,2)*tt108+H_invDmH(1,3)*tt1&
&07
tt112 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt111*tt26*tt20-lam(1,1)*tt&
&110*tt10*tt16*tt20+(mu(1,1)*(2*tt110*tt10*tt16-2*tt111*tt26))/2.0&
&E+0+lam(1,1)*tt110*tt10*tt16)
tt113 = H_invDmH(6,2)*tt3-H_invDmH(6,1)*tt2
tt114 = H_invDmH(6,3)*tt3-H_invDmH(6,1)*tt6
tt115 = H_invDmH(6,3)*tt2-H_invDmH(6,2)*tt6
tt116 = tt14*tt115-tt13*tt114+tt113*tt12
tt117 = H_invDmH(1,1)*tt115-H_invDmH(1,2)*tt114+H_invDmH(1,3)*tt1&
&13
tt118 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt117*tt26*tt20-lam(1,1)*tt&
&116*tt10*tt16*tt20+(mu(1,1)*(2*tt116*tt10*tt16-2*tt117*tt26))/2.0&
&E+0+lam(1,1)*tt116*tt10*tt16)
tt119 = H_invDmH(7,1)*tt9-H_invDmH(7,2)*tt8+H_invDmH(7,3)*tt5
tt120 = 2*H_invDmH(1,1)*H_invDmH(7,1)
tt121 = 2*H_invDmH(1,2)*H_invDmH(7,2)
tt122 = 2*H_invDmH(1,3)*H_invDmH(7,3)
tt123 = detDmH(1,1)*gw(1,1)*((-lam(1,1)*tt10*tt119*tt16*tt20)+(mu&
&(1,1)*(2*tt10*tt119*tt16+tt122+tt121+tt120))/2.0E+0+lam(1,1)*tt10&
&*tt119*tt16)
tt124 = H_invDmH(7,1)*tt4-H_invDmH(7,2)*tt1
tt125 = H_invDmH(7,1)*tt7-H_invDmH(7,3)*tt1
tt126 = H_invDmH(7,2)*tt7-H_invDmH(7,3)*tt4
tt127 = tt14*tt126-tt13*tt125+tt124*tt12
tt128 = H_invDmH(1,1)*tt126-H_invDmH(1,2)*tt125+H_invDmH(1,3)*tt1&
&24
tt129 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt128*tt26*tt20-lam(1,1)*tt&
&127*tt10*tt16*tt20+(mu(1,1)*(2*tt127*tt10*tt16-2*tt128*tt26))/2.0&
&E+0+lam(1,1)*tt127*tt10*tt16)
tt130 = H_invDmH(7,2)*tt3-H_invDmH(7,1)*tt2
tt131 = H_invDmH(7,3)*tt3-H_invDmH(7,1)*tt6
tt132 = H_invDmH(7,3)*tt2-H_invDmH(7,2)*tt6
tt133 = tt14*tt132-tt13*tt131+tt130*tt12
tt134 = H_invDmH(1,1)*tt132-H_invDmH(1,2)*tt131+H_invDmH(1,3)*tt1&
&30
tt135 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt134*tt26*tt20-lam(1,1)*tt&
&133*tt10*tt16*tt20+(mu(1,1)*(2*tt133*tt10*tt16-2*tt134*tt26))/2.0&
&E+0+lam(1,1)*tt133*tt10*tt16)
tt136 = H_invDmH(8,1)*tt9-H_invDmH(8,2)*tt8+tt5*H_invDmH(8,3)
tt137 = 2*H_invDmH(1,1)*H_invDmH(8,1)
tt138 = 2*H_invDmH(1,2)*H_invDmH(8,2)
tt139 = 2*H_invDmH(1,3)*H_invDmH(8,3)
tt140 = detDmH(1,1)*gw(1,1)*((-lam(1,1)*tt10*tt136*tt16*tt20)+(mu&
&(1,1)*(2*tt10*tt136*tt16+tt139+tt138+tt137))/2.0E+0+lam(1,1)*tt10&
&*tt136*tt16)
tt141 = H_invDmH(8,1)*tt4-tt1*H_invDmH(8,2)
tt142 = H_invDmH(8,1)*tt7-tt1*H_invDmH(8,3)
tt143 = H_invDmH(8,2)*tt7-tt4*H_invDmH(8,3)
tt144 = tt14*tt143-tt13*tt142+tt141*tt12
tt145 = H_invDmH(1,1)*tt143-H_invDmH(1,2)*tt142+H_invDmH(1,3)*tt1&
&41
tt146 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt145*tt26*tt20-lam(1,1)*tt&
&144*tt10*tt16*tt20+(mu(1,1)*(2*tt144*tt10*tt16-2*tt145*tt26))/2.0&
&E+0+lam(1,1)*tt144*tt10*tt16)
tt147 = tt3*H_invDmH(8,2)-H_invDmH(8,1)*tt2
tt148 = tt3*H_invDmH(8,3)-H_invDmH(8,1)*tt6
tt149 = tt2*H_invDmH(8,3)-H_invDmH(8,2)*tt6
tt150 = tt14*tt149-tt13*tt148+tt147*tt12
tt151 = H_invDmH(1,1)*tt149-H_invDmH(1,2)*tt148+H_invDmH(1,3)*tt1&
&47
tt152 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt151*tt26*tt20-lam(1,1)*tt&
&150*tt10*tt16*tt20+(mu(1,1)*(2*tt150*tt10*tt16-2*tt151*tt26))/2.0&
&E+0+lam(1,1)*tt150*tt10*tt16)
tt153 = tt24**2
tt154 = detDmH(1,1)*gw(1,1)*((-lam(1,1)*tt31*tt24*tt16*tt20)+mu(1&
&,1)*tt31*tt24*tt16+lam(1,1)*tt31*tt24*tt16)
tt155 = H_invDmH(2,1)*tt23-H_invDmH(2,2)*tt22+H_invDmH(2,3)*tt21
tt156 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt155*tt26*tt20-lam(1,1)*tt&
&24*tt34*tt16*tt20+(mu(1,1)*(2*tt24*tt34*tt16-2*tt155*tt26))/2.0E+&
&0+lam(1,1)*tt24*tt34*tt16)
tt157 = detDmH(1,1)*gw(1,1)*((-lam(1,1)*tt24*tt42*tt16*tt20)+(mu(&
&1,1)*(2*tt24*tt42*tt16+tt37+tt36+tt35))/2.0E+0+lam(1,1)*tt24*tt42&
&*tt16)
tt158 = (H_invDmH(1,1)*H_invDmH(2,2)-H_invDmH(1,2)*H_invDmH(2,1))&
&*tt12-(H_invDmH(1,1)*H_invDmH(2,3)-H_invDmH(1,3)*H_invDmH(2,1))*t&
&t13+(H_invDmH(1,2)*H_invDmH(2,3)-H_invDmH(1,3)*H_invDmH(2,2))*tt1&
&4
tt159 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt158*tt26*tt20-lam(1,1)*tt&
&48*tt24*tt16*tt20+(mu(1,1)*(2*tt48*tt24*tt16-2*tt158*tt26))/2.0E+&
&0+lam(1,1)*tt48*tt24*tt16)
tt160 = H_invDmH(3,1)*tt23-H_invDmH(3,2)*tt22+H_invDmH(3,3)*tt21
tt161 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt160*tt26*tt20-lam(1,1)*tt&
&24*tt51*tt16*tt20+(mu(1,1)*(2*tt24*tt51*tt16-2*tt160*tt26))/2.0E+&
&0+lam(1,1)*tt24*tt51*tt16)
tt162 = detDmH(1,1)*gw(1,1)*((-lam(1,1)*tt24*tt59*tt16*tt20)+(mu(&
&1,1)*(2*tt24*tt59*tt16+tt54+tt53+tt52))/2.0E+0+lam(1,1)*tt24*tt59&
&*tt16)
tt163 = (H_invDmH(1,1)*H_invDmH(3,2)-H_invDmH(1,2)*H_invDmH(3,1))&
&*tt12-(H_invDmH(1,1)*H_invDmH(3,3)-H_invDmH(1,3)*H_invDmH(3,1))*t&
&t13+(H_invDmH(1,2)*H_invDmH(3,3)-H_invDmH(1,3)*H_invDmH(3,2))*tt1&
&4
tt164 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt163*tt26*tt20-lam(1,1)*tt&
&65*tt24*tt16*tt20+(mu(1,1)*(2*tt65*tt24*tt16-2*tt163*tt26))/2.0E+&
&0+lam(1,1)*tt65*tt24*tt16)
tt165 = H_invDmH(4,1)*tt23-H_invDmH(4,2)*tt22+H_invDmH(4,3)*tt21
tt166 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt165*tt26*tt20-lam(1,1)*tt&
&24*tt68*tt16*tt20+(mu(1,1)*(2*tt24*tt68*tt16-2*tt165*tt26))/2.0E+&
&0+lam(1,1)*tt24*tt68*tt16)
tt167 = detDmH(1,1)*gw(1,1)*((-lam(1,1)*tt24*tt76*tt16*tt20)+(mu(&
&1,1)*(2*tt24*tt76*tt16+tt71+tt70+tt69))/2.0E+0+lam(1,1)*tt24*tt76&
&*tt16)
tt168 = (H_invDmH(1,1)*H_invDmH(4,2)-H_invDmH(1,2)*H_invDmH(4,1))&
&*tt12-(H_invDmH(1,1)*H_invDmH(4,3)-H_invDmH(1,3)*H_invDmH(4,1))*t&
&t13+(H_invDmH(1,2)*H_invDmH(4,3)-H_invDmH(1,3)*H_invDmH(4,2))*tt1&
&4
tt169 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt168*tt26*tt20-lam(1,1)*tt&
&82*tt24*tt16*tt20+(mu(1,1)*(2*tt82*tt24*tt16-2*tt168*tt26))/2.0E+&
&0+lam(1,1)*tt82*tt24*tt16)
tt170 = H_invDmH(5,1)*tt23-H_invDmH(5,2)*tt22+H_invDmH(5,3)*tt21
tt171 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt170*tt26*tt20-lam(1,1)*tt&
&24*tt85*tt16*tt20+(mu(1,1)*(2*tt24*tt85*tt16-2*tt170*tt26))/2.0E+&
&0+lam(1,1)*tt24*tt85*tt16)
tt172 = detDmH(1,1)*gw(1,1)*((-lam(1,1)*tt24*tt93*tt16*tt20)+(mu(&
&1,1)*(2*tt24*tt93*tt16+tt88+tt87+tt86))/2.0E+0+lam(1,1)*tt24*tt93&
&*tt16)
tt173 = (H_invDmH(1,1)*H_invDmH(5,2)-H_invDmH(1,2)*H_invDmH(5,1))&
&*tt12-(H_invDmH(1,1)*H_invDmH(5,3)-H_invDmH(1,3)*H_invDmH(5,1))*t&
&t13+(H_invDmH(1,2)*H_invDmH(5,3)-H_invDmH(1,3)*H_invDmH(5,2))*tt1&
&4
tt174 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt173*tt26*tt20-lam(1,1)*tt&
&99*tt24*tt16*tt20+(mu(1,1)*(2*tt99*tt24*tt16-2*tt173*tt26))/2.0E+&
&0+lam(1,1)*tt99*tt24*tt16)
tt175 = H_invDmH(6,1)*tt23-H_invDmH(6,2)*tt22+H_invDmH(6,3)*tt21
tt176 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt175*tt26*tt20-lam(1,1)*tt&
&24*tt102*tt16*tt20+(mu(1,1)*(2*tt24*tt102*tt16-2*tt175*tt26))/2.0&
&E+0+lam(1,1)*tt24*tt102*tt16)
tt177 = detDmH(1,1)*gw(1,1)*((-lam(1,1)*tt24*tt110*tt16*tt20)+(mu&
&(1,1)*(2*tt24*tt110*tt16+tt105+tt104+tt103))/2.0E+0+lam(1,1)*tt24&
&*tt110*tt16)
tt178 = (H_invDmH(1,1)*H_invDmH(6,2)-H_invDmH(1,2)*H_invDmH(6,1))&
&*tt12-(H_invDmH(1,1)*H_invDmH(6,3)-H_invDmH(1,3)*H_invDmH(6,1))*t&
&t13+(H_invDmH(1,2)*H_invDmH(6,3)-H_invDmH(1,3)*H_invDmH(6,2))*tt1&
&4
tt179 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt178*tt26*tt20-lam(1,1)*tt&
&116*tt24*tt16*tt20+(mu(1,1)*(2*tt116*tt24*tt16-2*tt178*tt26))/2.0&
&E+0+lam(1,1)*tt116*tt24*tt16)
tt180 = H_invDmH(7,1)*tt23-H_invDmH(7,2)*tt22+H_invDmH(7,3)*tt21
tt181 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt180*tt26*tt20-lam(1,1)*tt&
&24*tt119*tt16*tt20+(mu(1,1)*(2*tt24*tt119*tt16-2*tt180*tt26))/2.0&
&E+0+lam(1,1)*tt24*tt119*tt16)
tt182 = detDmH(1,1)*gw(1,1)*((-lam(1,1)*tt24*tt127*tt16*tt20)+(mu&
&(1,1)*(2*tt24*tt127*tt16+tt122+tt121+tt120))/2.0E+0+lam(1,1)*tt24&
&*tt127*tt16)
tt183 = (H_invDmH(1,1)*H_invDmH(7,2)-H_invDmH(1,2)*H_invDmH(7,1))&
&*tt12-(H_invDmH(1,1)*H_invDmH(7,3)-H_invDmH(1,3)*H_invDmH(7,1))*t&
&t13+(H_invDmH(1,2)*H_invDmH(7,3)-H_invDmH(1,3)*H_invDmH(7,2))*tt1&
&4
tt184 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt183*tt26*tt20-lam(1,1)*tt&
&133*tt24*tt16*tt20+(mu(1,1)*(2*tt133*tt24*tt16-2*tt183*tt26))/2.0&
&E+0+lam(1,1)*tt133*tt24*tt16)
tt185 = H_invDmH(8,1)*tt23-H_invDmH(8,2)*tt22+tt21*H_invDmH(8,3)
tt186 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt185*tt26*tt20-lam(1,1)*tt&
&24*tt136*tt16*tt20+(mu(1,1)*(2*tt24*tt136*tt16-2*tt185*tt26))/2.0&
&E+0+lam(1,1)*tt24*tt136*tt16)
tt187 = detDmH(1,1)*gw(1,1)*((-lam(1,1)*tt24*tt144*tt16*tt20)+(mu&
&(1,1)*(2*tt24*tt144*tt16+tt139+tt138+tt137))/2.0E+0+lam(1,1)*tt24&
&*tt144*tt16)
tt188 = (H_invDmH(1,1)*H_invDmH(8,2)-H_invDmH(1,2)*H_invDmH(8,1))&
&*tt12+tt14*(H_invDmH(1,2)*H_invDmH(8,3)-H_invDmH(1,3)*H_invDmH(8,&
&2))-tt13*(H_invDmH(1,1)*H_invDmH(8,3)-H_invDmH(1,3)*H_invDmH(8,1)&
&)
tt189 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt188*tt26*tt20-lam(1,1)*tt&
&150*tt24*tt16*tt20+(mu(1,1)*(2*tt150*tt24*tt16-2*tt188*tt26))/2.0&
&E+0+lam(1,1)*tt150*tt24*tt16)
tt190 = tt31**2
tt191 = H_invDmH(2,1)*tt30-H_invDmH(2,2)*tt29+H_invDmH(2,3)*tt28
tt192 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt191*tt26*tt20-lam(1,1)*tt&
&31*tt34*tt16*tt20+(mu(1,1)*(2*tt31*tt34*tt16-2*tt191*tt26))/2.0E+&
&0+lam(1,1)*tt31*tt34*tt16)
tt193 = (H_invDmH(1,2)*H_invDmH(2,1)-H_invDmH(1,1)*H_invDmH(2,2))&
&*tt12-(H_invDmH(1,3)*H_invDmH(2,1)-H_invDmH(1,1)*H_invDmH(2,3))*t&
&t13+(H_invDmH(1,3)*H_invDmH(2,2)-H_invDmH(1,2)*H_invDmH(2,3))*tt1&
&4
tt194 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt193*tt26*tt20-lam(1,1)*tt&
&31*tt42*tt16*tt20+(mu(1,1)*(2*tt31*tt42*tt16-2*tt193*tt26))/2.0E+&
&0+lam(1,1)*tt31*tt42*tt16)
tt195 = detDmH(1,1)*gw(1,1)*((-lam(1,1)*tt31*tt48*tt16*tt20)+(mu(&
&1,1)*(2*tt31*tt48*tt16+tt37+tt36+tt35))/2.0E+0+lam(1,1)*tt31*tt48&
&*tt16)
tt196 = H_invDmH(3,1)*tt30-H_invDmH(3,2)*tt29+H_invDmH(3,3)*tt28
tt197 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt196*tt26*tt20-lam(1,1)*tt&
&31*tt51*tt16*tt20+(mu(1,1)*(2*tt31*tt51*tt16-2*tt196*tt26))/2.0E+&
&0+lam(1,1)*tt31*tt51*tt16)
tt198 = (H_invDmH(1,2)*H_invDmH(3,1)-H_invDmH(1,1)*H_invDmH(3,2))&
&*tt12-(H_invDmH(1,3)*H_invDmH(3,1)-H_invDmH(1,1)*H_invDmH(3,3))*t&
&t13+(H_invDmH(1,3)*H_invDmH(3,2)-H_invDmH(1,2)*H_invDmH(3,3))*tt1&
&4
tt199 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt198*tt26*tt20-lam(1,1)*tt&
&31*tt59*tt16*tt20+(mu(1,1)*(2*tt31*tt59*tt16-2*tt198*tt26))/2.0E+&
&0+lam(1,1)*tt31*tt59*tt16)
tt200 = detDmH(1,1)*gw(1,1)*((-lam(1,1)*tt31*tt65*tt16*tt20)+(mu(&
&1,1)*(2*tt31*tt65*tt16+tt54+tt53+tt52))/2.0E+0+lam(1,1)*tt31*tt65&
&*tt16)
tt201 = H_invDmH(4,1)*tt30-H_invDmH(4,2)*tt29+H_invDmH(4,3)*tt28
tt202 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt201*tt26*tt20-lam(1,1)*tt&
&31*tt68*tt16*tt20+(mu(1,1)*(2*tt31*tt68*tt16-2*tt201*tt26))/2.0E+&
&0+lam(1,1)*tt31*tt68*tt16)
tt203 = (H_invDmH(1,2)*H_invDmH(4,1)-H_invDmH(1,1)*H_invDmH(4,2))&
&*tt12-(H_invDmH(1,3)*H_invDmH(4,1)-H_invDmH(1,1)*H_invDmH(4,3))*t&
&t13+(H_invDmH(1,3)*H_invDmH(4,2)-H_invDmH(1,2)*H_invDmH(4,3))*tt1&
&4
tt204 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt203*tt26*tt20-lam(1,1)*tt&
&31*tt76*tt16*tt20+(mu(1,1)*(2*tt31*tt76*tt16-2*tt203*tt26))/2.0E+&
&0+lam(1,1)*tt31*tt76*tt16)
tt205 = detDmH(1,1)*gw(1,1)*((-lam(1,1)*tt31*tt82*tt16*tt20)+(mu(&
&1,1)*(2*tt31*tt82*tt16+tt71+tt70+tt69))/2.0E+0+lam(1,1)*tt31*tt82&
&*tt16)
tt206 = H_invDmH(5,1)*tt30-H_invDmH(5,2)*tt29+H_invDmH(5,3)*tt28
tt207 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt206*tt26*tt20-lam(1,1)*tt&
&31*tt85*tt16*tt20+(mu(1,1)*(2*tt31*tt85*tt16-2*tt206*tt26))/2.0E+&
&0+lam(1,1)*tt31*tt85*tt16)
tt208 = (H_invDmH(1,2)*H_invDmH(5,1)-H_invDmH(1,1)*H_invDmH(5,2))&
&*tt12-(H_invDmH(1,3)*H_invDmH(5,1)-H_invDmH(1,1)*H_invDmH(5,3))*t&
&t13+(H_invDmH(1,3)*H_invDmH(5,2)-H_invDmH(1,2)*H_invDmH(5,3))*tt1&
&4
tt209 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt208*tt26*tt20-lam(1,1)*tt&
&31*tt93*tt16*tt20+(mu(1,1)*(2*tt31*tt93*tt16-2*tt208*tt26))/2.0E+&
&0+lam(1,1)*tt31*tt93*tt16)
tt210 = detDmH(1,1)*gw(1,1)*((-lam(1,1)*tt31*tt99*tt16*tt20)+(mu(&
&1,1)*(2*tt31*tt99*tt16+tt88+tt87+tt86))/2.0E+0+lam(1,1)*tt31*tt99&
&*tt16)
tt211 = H_invDmH(6,1)*tt30-H_invDmH(6,2)*tt29+H_invDmH(6,3)*tt28
tt212 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt211*tt26*tt20-lam(1,1)*tt&
&31*tt102*tt16*tt20+(mu(1,1)*(2*tt31*tt102*tt16-2*tt211*tt26))/2.0&
&E+0+lam(1,1)*tt31*tt102*tt16)
tt213 = (H_invDmH(1,2)*H_invDmH(6,1)-H_invDmH(1,1)*H_invDmH(6,2))&
&*tt12-(H_invDmH(1,3)*H_invDmH(6,1)-H_invDmH(1,1)*H_invDmH(6,3))*t&
&t13+(H_invDmH(1,3)*H_invDmH(6,2)-H_invDmH(1,2)*H_invDmH(6,3))*tt1&
&4
tt214 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt213*tt26*tt20-lam(1,1)*tt&
&31*tt110*tt16*tt20+(mu(1,1)*(2*tt31*tt110*tt16-2*tt213*tt26))/2.0&
&E+0+lam(1,1)*tt31*tt110*tt16)
tt215 = detDmH(1,1)*gw(1,1)*((-lam(1,1)*tt31*tt116*tt16*tt20)+(mu&
&(1,1)*(2*tt31*tt116*tt16+tt105+tt104+tt103))/2.0E+0+lam(1,1)*tt31&
&*tt116*tt16)
tt216 = H_invDmH(7,1)*tt30-H_invDmH(7,2)*tt29+H_invDmH(7,3)*tt28
tt217 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt216*tt26*tt20-lam(1,1)*tt&
&31*tt119*tt16*tt20+(mu(1,1)*(2*tt31*tt119*tt16-2*tt216*tt26))/2.0&
&E+0+lam(1,1)*tt31*tt119*tt16)
tt218 = (H_invDmH(1,2)*H_invDmH(7,1)-H_invDmH(1,1)*H_invDmH(7,2))&
&*tt12-(H_invDmH(1,3)*H_invDmH(7,1)-H_invDmH(1,1)*H_invDmH(7,3))*t&
&t13+(H_invDmH(1,3)*H_invDmH(7,2)-H_invDmH(1,2)*H_invDmH(7,3))*tt1&
&4
tt219 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt218*tt26*tt20-lam(1,1)*tt&
&31*tt127*tt16*tt20+(mu(1,1)*(2*tt31*tt127*tt16-2*tt218*tt26))/2.0&
&E+0+lam(1,1)*tt31*tt127*tt16)
tt220 = detDmH(1,1)*gw(1,1)*((-lam(1,1)*tt31*tt133*tt16*tt20)+(mu&
&(1,1)*(2*tt31*tt133*tt16+tt122+tt121+tt120))/2.0E+0+lam(1,1)*tt31&
&*tt133*tt16)
tt221 = H_invDmH(8,1)*tt30-H_invDmH(8,2)*tt29+tt28*H_invDmH(8,3)
tt222 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt221*tt26*tt20-lam(1,1)*tt&
&31*tt136*tt16*tt20+(mu(1,1)*(2*tt31*tt136*tt16-2*tt221*tt26))/2.0&
&E+0+lam(1,1)*tt31*tt136*tt16)
tt223 = (H_invDmH(1,2)*H_invDmH(8,1)-H_invDmH(1,1)*H_invDmH(8,2))&
&*tt12+tt14*(H_invDmH(1,3)*H_invDmH(8,2)-H_invDmH(1,2)*H_invDmH(8,&
&3))-tt13*(H_invDmH(1,3)*H_invDmH(8,1)-H_invDmH(1,1)*H_invDmH(8,3)&
&)
tt224 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt223*tt26*tt20-lam(1,1)*tt&
&31*tt144*tt16*tt20+(mu(1,1)*(2*tt31*tt144*tt16-2*tt223*tt26))/2.0&
&E+0+lam(1,1)*tt31*tt144*tt16)
tt225 = detDmH(1,1)*gw(1,1)*((-lam(1,1)*tt31*tt150*tt16*tt20)+(mu&
&(1,1)*(2*tt31*tt150*tt16+tt139+tt138+tt137))/2.0E+0+lam(1,1)*tt31&
&*tt150*tt16)
tt226 = tt34**2
tt227 = 2*H_invDmH(2,1)**2
tt228 = 2*H_invDmH(2,2)**2
tt229 = 2*H_invDmH(2,3)**2
tt230 = H_invDmH(2,1)*tt41-H_invDmH(2,2)*tt40+H_invDmH(2,3)*tt39
tt231 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt230*tt26*tt20-lam(1,1)*tt&
&42*tt34*tt16*tt20+(mu(1,1)*(2*tt42*tt34*tt16-2*tt230*tt26))/2.0E+&
&0+lam(1,1)*tt42*tt34*tt16)
tt232 = H_invDmH(2,1)*tt47-H_invDmH(2,2)*tt46+H_invDmH(2,3)*tt45
tt233 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt232*tt26*tt20-lam(1,1)*tt&
&48*tt34*tt16*tt20+(mu(1,1)*(2*tt48*tt34*tt16-2*tt232*tt26))/2.0E+&
&0+lam(1,1)*tt48*tt34*tt16)
tt234 = 2*H_invDmH(2,1)*H_invDmH(3,1)
tt235 = 2*H_invDmH(2,2)*H_invDmH(3,2)
tt236 = 2*H_invDmH(2,3)*H_invDmH(3,3)
tt237 = detDmH(1,1)*gw(1,1)*((-lam(1,1)*tt34*tt51*tt16*tt20)+(mu(&
&1,1)*(2*tt34*tt51*tt16+tt236+tt235+tt234))/2.0E+0+lam(1,1)*tt34*t&
&t51*tt16)
tt238 = H_invDmH(2,1)*tt58-H_invDmH(2,2)*tt57+H_invDmH(2,3)*tt56
tt239 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt238*tt26*tt20-lam(1,1)*tt&
&59*tt34*tt16*tt20+(mu(1,1)*(2*tt59*tt34*tt16-2*tt238*tt26))/2.0E+&
&0+lam(1,1)*tt59*tt34*tt16)
tt240 = H_invDmH(2,1)*tt64-H_invDmH(2,2)*tt63+H_invDmH(2,3)*tt62
tt241 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt240*tt26*tt20-lam(1,1)*tt&
&65*tt34*tt16*tt20+(mu(1,1)*(2*tt65*tt34*tt16-2*tt240*tt26))/2.0E+&
&0+lam(1,1)*tt65*tt34*tt16)
tt242 = 2*H_invDmH(2,1)*H_invDmH(4,1)
tt243 = 2*H_invDmH(2,2)*H_invDmH(4,2)
tt244 = 2*H_invDmH(2,3)*H_invDmH(4,3)
tt245 = detDmH(1,1)*gw(1,1)*((-lam(1,1)*tt34*tt68*tt16*tt20)+(mu(&
&1,1)*(2*tt34*tt68*tt16+tt244+tt243+tt242))/2.0E+0+lam(1,1)*tt34*t&
&t68*tt16)
tt246 = H_invDmH(2,1)*tt75-H_invDmH(2,2)*tt74+H_invDmH(2,3)*tt73
tt247 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt246*tt26*tt20-lam(1,1)*tt&
&76*tt34*tt16*tt20+(mu(1,1)*(2*tt76*tt34*tt16-2*tt246*tt26))/2.0E+&
&0+lam(1,1)*tt76*tt34*tt16)
tt248 = H_invDmH(2,1)*tt81-H_invDmH(2,2)*tt80+H_invDmH(2,3)*tt79
tt249 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt248*tt26*tt20-lam(1,1)*tt&
&82*tt34*tt16*tt20+(mu(1,1)*(2*tt82*tt34*tt16-2*tt248*tt26))/2.0E+&
&0+lam(1,1)*tt82*tt34*tt16)
tt250 = 2*H_invDmH(2,1)*H_invDmH(5,1)
tt251 = 2*H_invDmH(2,2)*H_invDmH(5,2)
tt252 = 2*H_invDmH(2,3)*H_invDmH(5,3)
tt253 = detDmH(1,1)*gw(1,1)*((-lam(1,1)*tt34*tt85*tt16*tt20)+(mu(&
&1,1)*(2*tt34*tt85*tt16+tt252+tt251+tt250))/2.0E+0+lam(1,1)*tt34*t&
&t85*tt16)
tt254 = H_invDmH(2,1)*tt92-H_invDmH(2,2)*tt91+H_invDmH(2,3)*tt90
tt255 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt254*tt26*tt20-lam(1,1)*tt&
&93*tt34*tt16*tt20+(mu(1,1)*(2*tt93*tt34*tt16-2*tt254*tt26))/2.0E+&
&0+lam(1,1)*tt93*tt34*tt16)
tt256 = H_invDmH(2,1)*tt98-H_invDmH(2,2)*tt97+H_invDmH(2,3)*tt96
tt257 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt256*tt26*tt20-lam(1,1)*tt&
&99*tt34*tt16*tt20+(mu(1,1)*(2*tt99*tt34*tt16-2*tt256*tt26))/2.0E+&
&0+lam(1,1)*tt99*tt34*tt16)
tt258 = 2*H_invDmH(2,1)*H_invDmH(6,1)
tt259 = 2*H_invDmH(2,2)*H_invDmH(6,2)
tt260 = 2*H_invDmH(2,3)*H_invDmH(6,3)
tt261 = detDmH(1,1)*gw(1,1)*((-lam(1,1)*tt34*tt102*tt16*tt20)+(mu&
&(1,1)*(2*tt34*tt102*tt16+tt260+tt259+tt258))/2.0E+0+lam(1,1)*tt34&
&*tt102*tt16)
tt262 = H_invDmH(2,1)*tt109-H_invDmH(2,2)*tt108+H_invDmH(2,3)*tt1&
&07
tt263 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt262*tt26*tt20-lam(1,1)*tt&
&110*tt34*tt16*tt20+(mu(1,1)*(2*tt110*tt34*tt16-2*tt262*tt26))/2.0&
&E+0+lam(1,1)*tt110*tt34*tt16)
tt264 = H_invDmH(2,1)*tt115-H_invDmH(2,2)*tt114+H_invDmH(2,3)*tt1&
&13
tt265 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt264*tt26*tt20-lam(1,1)*tt&
&116*tt34*tt16*tt20+(mu(1,1)*(2*tt116*tt34*tt16-2*tt264*tt26))/2.0&
&E+0+lam(1,1)*tt116*tt34*tt16)
tt266 = 2*H_invDmH(2,1)*H_invDmH(7,1)
tt267 = 2*H_invDmH(2,2)*H_invDmH(7,2)
tt268 = 2*H_invDmH(2,3)*H_invDmH(7,3)
tt269 = detDmH(1,1)*gw(1,1)*((-lam(1,1)*tt34*tt119*tt16*tt20)+(mu&
&(1,1)*(2*tt34*tt119*tt16+tt268+tt267+tt266))/2.0E+0+lam(1,1)*tt34&
&*tt119*tt16)
tt270 = H_invDmH(2,1)*tt126-H_invDmH(2,2)*tt125+H_invDmH(2,3)*tt1&
&24
tt271 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt270*tt26*tt20-lam(1,1)*tt&
&127*tt34*tt16*tt20+(mu(1,1)*(2*tt127*tt34*tt16-2*tt270*tt26))/2.0&
&E+0+lam(1,1)*tt127*tt34*tt16)
tt272 = H_invDmH(2,1)*tt132-H_invDmH(2,2)*tt131+H_invDmH(2,3)*tt1&
&30
tt273 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt272*tt26*tt20-lam(1,1)*tt&
&133*tt34*tt16*tt20+(mu(1,1)*(2*tt133*tt34*tt16-2*tt272*tt26))/2.0&
&E+0+lam(1,1)*tt133*tt34*tt16)
tt274 = 2*H_invDmH(2,1)*H_invDmH(8,1)
tt275 = 2*H_invDmH(2,2)*H_invDmH(8,2)
tt276 = 2*H_invDmH(2,3)*H_invDmH(8,3)
tt277 = detDmH(1,1)*gw(1,1)*((-lam(1,1)*tt34*tt136*tt16*tt20)+(mu&
&(1,1)*(2*tt34*tt136*tt16+tt276+tt275+tt274))/2.0E+0+lam(1,1)*tt34&
&*tt136*tt16)
tt278 = H_invDmH(2,1)*tt143-H_invDmH(2,2)*tt142+H_invDmH(2,3)*tt1&
&41
tt279 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt278*tt26*tt20-lam(1,1)*tt&
&144*tt34*tt16*tt20+(mu(1,1)*(2*tt144*tt34*tt16-2*tt278*tt26))/2.0&
&E+0+lam(1,1)*tt144*tt34*tt16)
tt280 = H_invDmH(2,1)*tt149-H_invDmH(2,2)*tt148+H_invDmH(2,3)*tt1&
&47
tt281 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt280*tt26*tt20-lam(1,1)*tt&
&150*tt34*tt16*tt20+(mu(1,1)*(2*tt150*tt34*tt16-2*tt280*tt26))/2.0&
&E+0+lam(1,1)*tt150*tt34*tt16)
tt282 = tt42**2
tt283 = detDmH(1,1)*gw(1,1)*((-lam(1,1)*tt48*tt42*tt16*tt20)+mu(1&
&,1)*tt48*tt42*tt16+lam(1,1)*tt48*tt42*tt16)
tt284 = H_invDmH(3,1)*tt41-H_invDmH(3,2)*tt40+H_invDmH(3,3)*tt39
tt285 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt284*tt26*tt20-lam(1,1)*tt&
&42*tt51*tt16*tt20+(mu(1,1)*(2*tt42*tt51*tt16-2*tt284*tt26))/2.0E+&
&0+lam(1,1)*tt42*tt51*tt16)
tt286 = detDmH(1,1)*gw(1,1)*((-lam(1,1)*tt42*tt59*tt16*tt20)+(mu(&
&1,1)*(2*tt42*tt59*tt16+tt236+tt235+tt234))/2.0E+0+lam(1,1)*tt42*t&
&t59*tt16)
tt287 = (H_invDmH(2,1)*H_invDmH(3,2)-H_invDmH(2,2)*H_invDmH(3,1))&
&*tt12-(H_invDmH(2,1)*H_invDmH(3,3)-H_invDmH(2,3)*H_invDmH(3,1))*t&
&t13+(H_invDmH(2,2)*H_invDmH(3,3)-H_invDmH(2,3)*H_invDmH(3,2))*tt1&
&4
tt288 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt287*tt26*tt20-lam(1,1)*tt&
&65*tt42*tt16*tt20+(mu(1,1)*(2*tt65*tt42*tt16-2*tt287*tt26))/2.0E+&
&0+lam(1,1)*tt65*tt42*tt16)
tt289 = H_invDmH(4,1)*tt41-H_invDmH(4,2)*tt40+H_invDmH(4,3)*tt39
tt290 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt289*tt26*tt20-lam(1,1)*tt&
&42*tt68*tt16*tt20+(mu(1,1)*(2*tt42*tt68*tt16-2*tt289*tt26))/2.0E+&
&0+lam(1,1)*tt42*tt68*tt16)
tt291 = detDmH(1,1)*gw(1,1)*((-lam(1,1)*tt42*tt76*tt16*tt20)+(mu(&
&1,1)*(2*tt42*tt76*tt16+tt244+tt243+tt242))/2.0E+0+lam(1,1)*tt42*t&
&t76*tt16)
tt292 = (H_invDmH(2,1)*H_invDmH(4,2)-H_invDmH(2,2)*H_invDmH(4,1))&
&*tt12-(H_invDmH(2,1)*H_invDmH(4,3)-H_invDmH(2,3)*H_invDmH(4,1))*t&
&t13+(H_invDmH(2,2)*H_invDmH(4,3)-H_invDmH(2,3)*H_invDmH(4,2))*tt1&
&4
tt293 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt292*tt26*tt20-lam(1,1)*tt&
&82*tt42*tt16*tt20+(mu(1,1)*(2*tt82*tt42*tt16-2*tt292*tt26))/2.0E+&
&0+lam(1,1)*tt82*tt42*tt16)
tt294 = H_invDmH(5,1)*tt41-H_invDmH(5,2)*tt40+H_invDmH(5,3)*tt39
tt295 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt294*tt26*tt20-lam(1,1)*tt&
&42*tt85*tt16*tt20+(mu(1,1)*(2*tt42*tt85*tt16-2*tt294*tt26))/2.0E+&
&0+lam(1,1)*tt42*tt85*tt16)
tt296 = detDmH(1,1)*gw(1,1)*((-lam(1,1)*tt42*tt93*tt16*tt20)+(mu(&
&1,1)*(2*tt42*tt93*tt16+tt252+tt251+tt250))/2.0E+0+lam(1,1)*tt42*t&
&t93*tt16)
tt297 = (H_invDmH(2,1)*H_invDmH(5,2)-H_invDmH(2,2)*H_invDmH(5,1))&
&*tt12-(H_invDmH(2,1)*H_invDmH(5,3)-H_invDmH(2,3)*H_invDmH(5,1))*t&
&t13+(H_invDmH(2,2)*H_invDmH(5,3)-H_invDmH(2,3)*H_invDmH(5,2))*tt1&
&4
tt298 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt297*tt26*tt20-lam(1,1)*tt&
&99*tt42*tt16*tt20+(mu(1,1)*(2*tt99*tt42*tt16-2*tt297*tt26))/2.0E+&
&0+lam(1,1)*tt99*tt42*tt16)
tt299 = H_invDmH(6,1)*tt41-H_invDmH(6,2)*tt40+H_invDmH(6,3)*tt39
tt300 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt299*tt26*tt20-lam(1,1)*tt&
&42*tt102*tt16*tt20+(mu(1,1)*(2*tt42*tt102*tt16-2*tt299*tt26))/2.0&
&E+0+lam(1,1)*tt42*tt102*tt16)
tt301 = detDmH(1,1)*gw(1,1)*((-lam(1,1)*tt42*tt110*tt16*tt20)+(mu&
&(1,1)*(2*tt42*tt110*tt16+tt260+tt259+tt258))/2.0E+0+lam(1,1)*tt42&
&*tt110*tt16)
tt302 = (H_invDmH(2,1)*H_invDmH(6,2)-H_invDmH(2,2)*H_invDmH(6,1))&
&*tt12-(H_invDmH(2,1)*H_invDmH(6,3)-H_invDmH(2,3)*H_invDmH(6,1))*t&
&t13+(H_invDmH(2,2)*H_invDmH(6,3)-H_invDmH(2,3)*H_invDmH(6,2))*tt1&
&4
tt303 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt302*tt26*tt20-lam(1,1)*tt&
&116*tt42*tt16*tt20+(mu(1,1)*(2*tt116*tt42*tt16-2*tt302*tt26))/2.0&
&E+0+lam(1,1)*tt116*tt42*tt16)
tt304 = H_invDmH(7,1)*tt41-H_invDmH(7,2)*tt40+H_invDmH(7,3)*tt39
tt305 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt304*tt26*tt20-lam(1,1)*tt&
&42*tt119*tt16*tt20+(mu(1,1)*(2*tt42*tt119*tt16-2*tt304*tt26))/2.0&
&E+0+lam(1,1)*tt42*tt119*tt16)
tt306 = detDmH(1,1)*gw(1,1)*((-lam(1,1)*tt42*tt127*tt16*tt20)+(mu&
&(1,1)*(2*tt42*tt127*tt16+tt268+tt267+tt266))/2.0E+0+lam(1,1)*tt42&
&*tt127*tt16)
tt307 = (H_invDmH(2,1)*H_invDmH(7,2)-H_invDmH(2,2)*H_invDmH(7,1))&
&*tt12-(H_invDmH(2,1)*H_invDmH(7,3)-H_invDmH(2,3)*H_invDmH(7,1))*t&
&t13+(H_invDmH(2,2)*H_invDmH(7,3)-H_invDmH(2,3)*H_invDmH(7,2))*tt1&
&4
tt308 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt307*tt26*tt20-lam(1,1)*tt&
&133*tt42*tt16*tt20+(mu(1,1)*(2*tt133*tt42*tt16-2*tt307*tt26))/2.0&
&E+0+lam(1,1)*tt133*tt42*tt16)
tt309 = H_invDmH(8,1)*tt41-H_invDmH(8,2)*tt40+tt39*H_invDmH(8,3)
tt310 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt309*tt26*tt20-lam(1,1)*tt&
&42*tt136*tt16*tt20+(mu(1,1)*(2*tt42*tt136*tt16-2*tt309*tt26))/2.0&
&E+0+lam(1,1)*tt42*tt136*tt16)
tt311 = detDmH(1,1)*gw(1,1)*((-lam(1,1)*tt42*tt144*tt16*tt20)+(mu&
&(1,1)*(2*tt42*tt144*tt16+tt276+tt275+tt274))/2.0E+0+lam(1,1)*tt42&
&*tt144*tt16)
tt312 = tt14*(H_invDmH(2,2)*H_invDmH(8,3)-H_invDmH(2,3)*H_invDmH(&
&8,2))-tt13*(H_invDmH(2,1)*H_invDmH(8,3)-H_invDmH(2,3)*H_invDmH(8,&
&1))+(H_invDmH(2,1)*H_invDmH(8,2)-H_invDmH(2,2)*H_invDmH(8,1))*tt1&
&2
tt313 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt312*tt26*tt20-lam(1,1)*tt&
&150*tt42*tt16*tt20+(mu(1,1)*(2*tt150*tt42*tt16-2*tt312*tt26))/2.0&
&E+0+lam(1,1)*tt150*tt42*tt16)
tt314 = tt48**2
tt315 = H_invDmH(3,1)*tt47-H_invDmH(3,2)*tt46+H_invDmH(3,3)*tt45
tt316 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt315*tt26*tt20-lam(1,1)*tt&
&48*tt51*tt16*tt20+(mu(1,1)*(2*tt48*tt51*tt16-2*tt315*tt26))/2.0E+&
&0+lam(1,1)*tt48*tt51*tt16)
tt317 = (H_invDmH(2,2)*H_invDmH(3,1)-H_invDmH(2,1)*H_invDmH(3,2))&
&*tt12-(H_invDmH(2,3)*H_invDmH(3,1)-H_invDmH(2,1)*H_invDmH(3,3))*t&
&t13+(H_invDmH(2,3)*H_invDmH(3,2)-H_invDmH(2,2)*H_invDmH(3,3))*tt1&
&4
tt318 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt317*tt26*tt20-lam(1,1)*tt&
&48*tt59*tt16*tt20+(mu(1,1)*(2*tt48*tt59*tt16-2*tt317*tt26))/2.0E+&
&0+lam(1,1)*tt48*tt59*tt16)
tt319 = detDmH(1,1)*gw(1,1)*((-lam(1,1)*tt48*tt65*tt16*tt20)+(mu(&
&1,1)*(2*tt48*tt65*tt16+tt236+tt235+tt234))/2.0E+0+lam(1,1)*tt48*t&
&t65*tt16)
tt320 = H_invDmH(4,1)*tt47-H_invDmH(4,2)*tt46+H_invDmH(4,3)*tt45
tt321 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt320*tt26*tt20-lam(1,1)*tt&
&48*tt68*tt16*tt20+(mu(1,1)*(2*tt48*tt68*tt16-2*tt320*tt26))/2.0E+&
&0+lam(1,1)*tt48*tt68*tt16)
tt322 = (H_invDmH(2,2)*H_invDmH(4,1)-H_invDmH(2,1)*H_invDmH(4,2))&
&*tt12-(H_invDmH(2,3)*H_invDmH(4,1)-H_invDmH(2,1)*H_invDmH(4,3))*t&
&t13+(H_invDmH(2,3)*H_invDmH(4,2)-H_invDmH(2,2)*H_invDmH(4,3))*tt1&
&4
tt323 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt322*tt26*tt20-lam(1,1)*tt&
&48*tt76*tt16*tt20+(mu(1,1)*(2*tt48*tt76*tt16-2*tt322*tt26))/2.0E+&
&0+lam(1,1)*tt48*tt76*tt16)
tt324 = detDmH(1,1)*gw(1,1)*((-lam(1,1)*tt48*tt82*tt16*tt20)+(mu(&
&1,1)*(2*tt48*tt82*tt16+tt244+tt243+tt242))/2.0E+0+lam(1,1)*tt48*t&
&t82*tt16)
tt325 = H_invDmH(5,1)*tt47-H_invDmH(5,2)*tt46+H_invDmH(5,3)*tt45
tt326 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt325*tt26*tt20-lam(1,1)*tt&
&48*tt85*tt16*tt20+(mu(1,1)*(2*tt48*tt85*tt16-2*tt325*tt26))/2.0E+&
&0+lam(1,1)*tt48*tt85*tt16)
tt327 = (H_invDmH(2,2)*H_invDmH(5,1)-H_invDmH(2,1)*H_invDmH(5,2))&
&*tt12-(H_invDmH(2,3)*H_invDmH(5,1)-H_invDmH(2,1)*H_invDmH(5,3))*t&
&t13+(H_invDmH(2,3)*H_invDmH(5,2)-H_invDmH(2,2)*H_invDmH(5,3))*tt1&
&4
tt328 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt327*tt26*tt20-lam(1,1)*tt&
&48*tt93*tt16*tt20+(mu(1,1)*(2*tt48*tt93*tt16-2*tt327*tt26))/2.0E+&
&0+lam(1,1)*tt48*tt93*tt16)
tt329 = detDmH(1,1)*gw(1,1)*((-lam(1,1)*tt48*tt99*tt16*tt20)+(mu(&
&1,1)*(2*tt48*tt99*tt16+tt252+tt251+tt250))/2.0E+0+lam(1,1)*tt48*t&
&t99*tt16)
tt330 = H_invDmH(6,1)*tt47-H_invDmH(6,2)*tt46+H_invDmH(6,3)*tt45
tt331 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt330*tt26*tt20-lam(1,1)*tt&
&48*tt102*tt16*tt20+(mu(1,1)*(2*tt48*tt102*tt16-2*tt330*tt26))/2.0&
&E+0+lam(1,1)*tt48*tt102*tt16)
tt332 = (H_invDmH(2,2)*H_invDmH(6,1)-H_invDmH(2,1)*H_invDmH(6,2))&
&*tt12-(H_invDmH(2,3)*H_invDmH(6,1)-H_invDmH(2,1)*H_invDmH(6,3))*t&
&t13+(H_invDmH(2,3)*H_invDmH(6,2)-H_invDmH(2,2)*H_invDmH(6,3))*tt1&
&4
tt333 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt332*tt26*tt20-lam(1,1)*tt&
&48*tt110*tt16*tt20+(mu(1,1)*(2*tt48*tt110*tt16-2*tt332*tt26))/2.0&
&E+0+lam(1,1)*tt48*tt110*tt16)
tt334 = detDmH(1,1)*gw(1,1)*((-lam(1,1)*tt48*tt116*tt16*tt20)+(mu&
&(1,1)*(2*tt48*tt116*tt16+tt260+tt259+tt258))/2.0E+0+lam(1,1)*tt48&
&*tt116*tt16)
tt335 = H_invDmH(7,1)*tt47-H_invDmH(7,2)*tt46+H_invDmH(7,3)*tt45
tt336 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt335*tt26*tt20-lam(1,1)*tt&
&48*tt119*tt16*tt20+(mu(1,1)*(2*tt48*tt119*tt16-2*tt335*tt26))/2.0&
&E+0+lam(1,1)*tt48*tt119*tt16)
tt337 = (H_invDmH(2,2)*H_invDmH(7,1)-H_invDmH(2,1)*H_invDmH(7,2))&
&*tt12-(H_invDmH(2,3)*H_invDmH(7,1)-H_invDmH(2,1)*H_invDmH(7,3))*t&
&t13+(H_invDmH(2,3)*H_invDmH(7,2)-H_invDmH(2,2)*H_invDmH(7,3))*tt1&
&4
tt338 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt337*tt26*tt20-lam(1,1)*tt&
&48*tt127*tt16*tt20+(mu(1,1)*(2*tt48*tt127*tt16-2*tt337*tt26))/2.0&
&E+0+lam(1,1)*tt48*tt127*tt16)
tt339 = detDmH(1,1)*gw(1,1)*((-lam(1,1)*tt48*tt133*tt16*tt20)+(mu&
&(1,1)*(2*tt48*tt133*tt16+tt268+tt267+tt266))/2.0E+0+lam(1,1)*tt48&
&*tt133*tt16)
tt340 = H_invDmH(8,1)*tt47-H_invDmH(8,2)*tt46+tt45*H_invDmH(8,3)
tt341 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt340*tt26*tt20-lam(1,1)*tt&
&48*tt136*tt16*tt20+(mu(1,1)*(2*tt48*tt136*tt16-2*tt340*tt26))/2.0&
&E+0+lam(1,1)*tt48*tt136*tt16)
tt342 = tt14*(H_invDmH(2,3)*H_invDmH(8,2)-H_invDmH(2,2)*H_invDmH(&
&8,3))-tt13*(H_invDmH(2,3)*H_invDmH(8,1)-H_invDmH(2,1)*H_invDmH(8,&
&3))+(H_invDmH(2,2)*H_invDmH(8,1)-H_invDmH(2,1)*H_invDmH(8,2))*tt1&
&2
tt343 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt342*tt26*tt20-lam(1,1)*tt&
&48*tt144*tt16*tt20+(mu(1,1)*(2*tt48*tt144*tt16-2*tt342*tt26))/2.0&
&E+0+lam(1,1)*tt48*tt144*tt16)
tt344 = detDmH(1,1)*gw(1,1)*((-lam(1,1)*tt48*tt150*tt16*tt20)+(mu&
&(1,1)*(2*tt48*tt150*tt16+tt276+tt275+tt274))/2.0E+0+lam(1,1)*tt48&
&*tt150*tt16)
tt345 = tt51**2
tt346 = 2*H_invDmH(3,1)**2
tt347 = 2*H_invDmH(3,2)**2
tt348 = 2*H_invDmH(3,3)**2
tt349 = H_invDmH(3,1)*tt58-H_invDmH(3,2)*tt57+H_invDmH(3,3)*tt56
tt350 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt349*tt26*tt20-lam(1,1)*tt&
&59*tt51*tt16*tt20+(mu(1,1)*(2*tt59*tt51*tt16-2*tt349*tt26))/2.0E+&
&0+lam(1,1)*tt59*tt51*tt16)
tt351 = H_invDmH(3,1)*tt64-H_invDmH(3,2)*tt63+H_invDmH(3,3)*tt62
tt352 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt351*tt26*tt20-lam(1,1)*tt&
&65*tt51*tt16*tt20+(mu(1,1)*(2*tt65*tt51*tt16-2*tt351*tt26))/2.0E+&
&0+lam(1,1)*tt65*tt51*tt16)
tt353 = 2*H_invDmH(3,1)*H_invDmH(4,1)
tt354 = 2*H_invDmH(3,2)*H_invDmH(4,2)
tt355 = 2*H_invDmH(3,3)*H_invDmH(4,3)
tt356 = detDmH(1,1)*gw(1,1)*((-lam(1,1)*tt51*tt68*tt16*tt20)+(mu(&
&1,1)*(2*tt51*tt68*tt16+tt355+tt354+tt353))/2.0E+0+lam(1,1)*tt51*t&
&t68*tt16)
tt357 = H_invDmH(3,1)*tt75-H_invDmH(3,2)*tt74+H_invDmH(3,3)*tt73
tt358 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt357*tt26*tt20-lam(1,1)*tt&
&76*tt51*tt16*tt20+(mu(1,1)*(2*tt76*tt51*tt16-2*tt357*tt26))/2.0E+&
&0+lam(1,1)*tt76*tt51*tt16)
tt359 = H_invDmH(3,1)*tt81-H_invDmH(3,2)*tt80+H_invDmH(3,3)*tt79
tt360 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt359*tt26*tt20-lam(1,1)*tt&
&82*tt51*tt16*tt20+(mu(1,1)*(2*tt82*tt51*tt16-2*tt359*tt26))/2.0E+&
&0+lam(1,1)*tt82*tt51*tt16)
tt361 = 2*H_invDmH(3,1)*H_invDmH(5,1)
tt362 = 2*H_invDmH(3,2)*H_invDmH(5,2)
tt363 = 2*H_invDmH(3,3)*H_invDmH(5,3)
tt364 = detDmH(1,1)*gw(1,1)*((-lam(1,1)*tt51*tt85*tt16*tt20)+(mu(&
&1,1)*(2*tt51*tt85*tt16+tt363+tt362+tt361))/2.0E+0+lam(1,1)*tt51*t&
&t85*tt16)
tt365 = H_invDmH(3,1)*tt92-H_invDmH(3,2)*tt91+H_invDmH(3,3)*tt90
tt366 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt365*tt26*tt20-lam(1,1)*tt&
&93*tt51*tt16*tt20+(mu(1,1)*(2*tt93*tt51*tt16-2*tt365*tt26))/2.0E+&
&0+lam(1,1)*tt93*tt51*tt16)
tt367 = H_invDmH(3,1)*tt98-H_invDmH(3,2)*tt97+H_invDmH(3,3)*tt96
tt368 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt367*tt26*tt20-lam(1,1)*tt&
&99*tt51*tt16*tt20+(mu(1,1)*(2*tt99*tt51*tt16-2*tt367*tt26))/2.0E+&
&0+lam(1,1)*tt99*tt51*tt16)
tt369 = 2*H_invDmH(3,1)*H_invDmH(6,1)
tt370 = 2*H_invDmH(3,2)*H_invDmH(6,2)
tt371 = 2*H_invDmH(3,3)*H_invDmH(6,3)
tt372 = detDmH(1,1)*gw(1,1)*((-lam(1,1)*tt51*tt102*tt16*tt20)+(mu&
&(1,1)*(2*tt51*tt102*tt16+tt371+tt370+tt369))/2.0E+0+lam(1,1)*tt51&
&*tt102*tt16)
tt373 = H_invDmH(3,1)*tt109-H_invDmH(3,2)*tt108+H_invDmH(3,3)*tt1&
&07
tt374 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt373*tt26*tt20-lam(1,1)*tt&
&110*tt51*tt16*tt20+(mu(1,1)*(2*tt110*tt51*tt16-2*tt373*tt26))/2.0&
&E+0+lam(1,1)*tt110*tt51*tt16)
tt375 = H_invDmH(3,1)*tt115-H_invDmH(3,2)*tt114+H_invDmH(3,3)*tt1&
&13
tt376 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt375*tt26*tt20-lam(1,1)*tt&
&116*tt51*tt16*tt20+(mu(1,1)*(2*tt116*tt51*tt16-2*tt375*tt26))/2.0&
&E+0+lam(1,1)*tt116*tt51*tt16)
tt377 = 2*H_invDmH(3,1)*H_invDmH(7,1)
tt378 = 2*H_invDmH(3,2)*H_invDmH(7,2)
tt379 = 2*H_invDmH(3,3)*H_invDmH(7,3)
tt380 = detDmH(1,1)*gw(1,1)*((-lam(1,1)*tt51*tt119*tt16*tt20)+(mu&
&(1,1)*(2*tt51*tt119*tt16+tt379+tt378+tt377))/2.0E+0+lam(1,1)*tt51&
&*tt119*tt16)
tt381 = H_invDmH(3,1)*tt126-H_invDmH(3,2)*tt125+H_invDmH(3,3)*tt1&
&24
tt382 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt381*tt26*tt20-lam(1,1)*tt&
&127*tt51*tt16*tt20+(mu(1,1)*(2*tt127*tt51*tt16-2*tt381*tt26))/2.0&
&E+0+lam(1,1)*tt127*tt51*tt16)
tt383 = H_invDmH(3,1)*tt132-H_invDmH(3,2)*tt131+H_invDmH(3,3)*tt1&
&30
tt384 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt383*tt26*tt20-lam(1,1)*tt&
&133*tt51*tt16*tt20+(mu(1,1)*(2*tt133*tt51*tt16-2*tt383*tt26))/2.0&
&E+0+lam(1,1)*tt133*tt51*tt16)
tt385 = 2*H_invDmH(3,1)*H_invDmH(8,1)
tt386 = 2*H_invDmH(3,2)*H_invDmH(8,2)
tt387 = 2*H_invDmH(3,3)*H_invDmH(8,3)
tt388 = detDmH(1,1)*gw(1,1)*((-lam(1,1)*tt51*tt136*tt16*tt20)+(mu&
&(1,1)*(2*tt51*tt136*tt16+tt387+tt386+tt385))/2.0E+0+lam(1,1)*tt51&
&*tt136*tt16)
tt389 = H_invDmH(3,1)*tt143-H_invDmH(3,2)*tt142+H_invDmH(3,3)*tt1&
&41
tt390 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt389*tt26*tt20-lam(1,1)*tt&
&144*tt51*tt16*tt20+(mu(1,1)*(2*tt144*tt51*tt16-2*tt389*tt26))/2.0&
&E+0+lam(1,1)*tt144*tt51*tt16)
tt391 = H_invDmH(3,1)*tt149-H_invDmH(3,2)*tt148+H_invDmH(3,3)*tt1&
&47
tt392 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt391*tt26*tt20-lam(1,1)*tt&
&150*tt51*tt16*tt20+(mu(1,1)*(2*tt150*tt51*tt16-2*tt391*tt26))/2.0&
&E+0+lam(1,1)*tt150*tt51*tt16)
tt393 = tt59**2
tt394 = detDmH(1,1)*gw(1,1)*((-lam(1,1)*tt65*tt59*tt16*tt20)+mu(1&
&,1)*tt65*tt59*tt16+lam(1,1)*tt65*tt59*tt16)
tt395 = H_invDmH(4,1)*tt58-H_invDmH(4,2)*tt57+H_invDmH(4,3)*tt56
tt396 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt395*tt26*tt20-lam(1,1)*tt&
&59*tt68*tt16*tt20+(mu(1,1)*(2*tt59*tt68*tt16-2*tt395*tt26))/2.0E+&
&0+lam(1,1)*tt59*tt68*tt16)
tt397 = detDmH(1,1)*gw(1,1)*((-lam(1,1)*tt59*tt76*tt16*tt20)+(mu(&
&1,1)*(2*tt59*tt76*tt16+tt355+tt354+tt353))/2.0E+0+lam(1,1)*tt59*t&
&t76*tt16)
tt398 = (H_invDmH(3,1)*H_invDmH(4,2)-H_invDmH(3,2)*H_invDmH(4,1))&
&*tt12-(H_invDmH(3,1)*H_invDmH(4,3)-H_invDmH(3,3)*H_invDmH(4,1))*t&
&t13+(H_invDmH(3,2)*H_invDmH(4,3)-H_invDmH(3,3)*H_invDmH(4,2))*tt1&
&4
tt399 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt398*tt26*tt20-lam(1,1)*tt&
&82*tt59*tt16*tt20+(mu(1,1)*(2*tt82*tt59*tt16-2*tt398*tt26))/2.0E+&
&0+lam(1,1)*tt82*tt59*tt16)
tt400 = H_invDmH(5,1)*tt58-H_invDmH(5,2)*tt57+H_invDmH(5,3)*tt56
tt401 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt400*tt26*tt20-lam(1,1)*tt&
&59*tt85*tt16*tt20+(mu(1,1)*(2*tt59*tt85*tt16-2*tt400*tt26))/2.0E+&
&0+lam(1,1)*tt59*tt85*tt16)
tt402 = detDmH(1,1)*gw(1,1)*((-lam(1,1)*tt59*tt93*tt16*tt20)+(mu(&
&1,1)*(2*tt59*tt93*tt16+tt363+tt362+tt361))/2.0E+0+lam(1,1)*tt59*t&
&t93*tt16)
tt403 = (H_invDmH(3,1)*H_invDmH(5,2)-H_invDmH(3,2)*H_invDmH(5,1))&
&*tt12-(H_invDmH(3,1)*H_invDmH(5,3)-H_invDmH(3,3)*H_invDmH(5,1))*t&
&t13+(H_invDmH(3,2)*H_invDmH(5,3)-H_invDmH(3,3)*H_invDmH(5,2))*tt1&
&4
tt404 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt403*tt26*tt20-lam(1,1)*tt&
&99*tt59*tt16*tt20+(mu(1,1)*(2*tt99*tt59*tt16-2*tt403*tt26))/2.0E+&
&0+lam(1,1)*tt99*tt59*tt16)
tt405 = H_invDmH(6,1)*tt58-H_invDmH(6,2)*tt57+H_invDmH(6,3)*tt56
tt406 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt405*tt26*tt20-lam(1,1)*tt&
&59*tt102*tt16*tt20+(mu(1,1)*(2*tt59*tt102*tt16-2*tt405*tt26))/2.0&
&E+0+lam(1,1)*tt59*tt102*tt16)
tt407 = detDmH(1,1)*gw(1,1)*((-lam(1,1)*tt59*tt110*tt16*tt20)+(mu&
&(1,1)*(2*tt59*tt110*tt16+tt371+tt370+tt369))/2.0E+0+lam(1,1)*tt59&
&*tt110*tt16)
tt408 = (H_invDmH(3,1)*H_invDmH(6,2)-H_invDmH(3,2)*H_invDmH(6,1))&
&*tt12-(H_invDmH(3,1)*H_invDmH(6,3)-H_invDmH(3,3)*H_invDmH(6,1))*t&
&t13+(H_invDmH(3,2)*H_invDmH(6,3)-H_invDmH(3,3)*H_invDmH(6,2))*tt1&
&4
tt409 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt408*tt26*tt20-lam(1,1)*tt&
&116*tt59*tt16*tt20+(mu(1,1)*(2*tt116*tt59*tt16-2*tt408*tt26))/2.0&
&E+0+lam(1,1)*tt116*tt59*tt16)
tt410 = H_invDmH(7,1)*tt58-H_invDmH(7,2)*tt57+H_invDmH(7,3)*tt56
tt411 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt410*tt26*tt20-lam(1,1)*tt&
&59*tt119*tt16*tt20+(mu(1,1)*(2*tt59*tt119*tt16-2*tt410*tt26))/2.0&
&E+0+lam(1,1)*tt59*tt119*tt16)
tt412 = detDmH(1,1)*gw(1,1)*((-lam(1,1)*tt59*tt127*tt16*tt20)+(mu&
&(1,1)*(2*tt59*tt127*tt16+tt379+tt378+tt377))/2.0E+0+lam(1,1)*tt59&
&*tt127*tt16)
tt413 = (H_invDmH(3,1)*H_invDmH(7,2)-H_invDmH(3,2)*H_invDmH(7,1))&
&*tt12-(H_invDmH(3,1)*H_invDmH(7,3)-H_invDmH(3,3)*H_invDmH(7,1))*t&
&t13+(H_invDmH(3,2)*H_invDmH(7,3)-H_invDmH(3,3)*H_invDmH(7,2))*tt1&
&4
tt414 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt413*tt26*tt20-lam(1,1)*tt&
&133*tt59*tt16*tt20+(mu(1,1)*(2*tt133*tt59*tt16-2*tt413*tt26))/2.0&
&E+0+lam(1,1)*tt133*tt59*tt16)
tt415 = H_invDmH(8,1)*tt58-H_invDmH(8,2)*tt57+tt56*H_invDmH(8,3)
tt416 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt415*tt26*tt20-lam(1,1)*tt&
&59*tt136*tt16*tt20+(mu(1,1)*(2*tt59*tt136*tt16-2*tt415*tt26))/2.0&
&E+0+lam(1,1)*tt59*tt136*tt16)
tt417 = detDmH(1,1)*gw(1,1)*((-lam(1,1)*tt59*tt144*tt16*tt20)+(mu&
&(1,1)*(2*tt59*tt144*tt16+tt387+tt386+tt385))/2.0E+0+lam(1,1)*tt59&
&*tt144*tt16)
tt418 = tt14*(H_invDmH(3,2)*H_invDmH(8,3)-H_invDmH(3,3)*H_invDmH(&
&8,2))-tt13*(H_invDmH(3,1)*H_invDmH(8,3)-H_invDmH(3,3)*H_invDmH(8,&
&1))+(H_invDmH(3,1)*H_invDmH(8,2)-H_invDmH(3,2)*H_invDmH(8,1))*tt1&
&2
tt419 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt418*tt26*tt20-lam(1,1)*tt&
&150*tt59*tt16*tt20+(mu(1,1)*(2*tt150*tt59*tt16-2*tt418*tt26))/2.0&
&E+0+lam(1,1)*tt150*tt59*tt16)
tt420 = tt65**2
tt421 = H_invDmH(4,1)*tt64-H_invDmH(4,2)*tt63+H_invDmH(4,3)*tt62
tt422 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt421*tt26*tt20-lam(1,1)*tt&
&65*tt68*tt16*tt20+(mu(1,1)*(2*tt65*tt68*tt16-2*tt421*tt26))/2.0E+&
&0+lam(1,1)*tt65*tt68*tt16)
tt423 = (H_invDmH(3,2)*H_invDmH(4,1)-H_invDmH(3,1)*H_invDmH(4,2))&
&*tt12-(H_invDmH(3,3)*H_invDmH(4,1)-H_invDmH(3,1)*H_invDmH(4,3))*t&
&t13+(H_invDmH(3,3)*H_invDmH(4,2)-H_invDmH(3,2)*H_invDmH(4,3))*tt1&
&4
tt424 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt423*tt26*tt20-lam(1,1)*tt&
&65*tt76*tt16*tt20+(mu(1,1)*(2*tt65*tt76*tt16-2*tt423*tt26))/2.0E+&
&0+lam(1,1)*tt65*tt76*tt16)
tt425 = detDmH(1,1)*gw(1,1)*((-lam(1,1)*tt65*tt82*tt16*tt20)+(mu(&
&1,1)*(2*tt65*tt82*tt16+tt355+tt354+tt353))/2.0E+0+lam(1,1)*tt65*t&
&t82*tt16)
tt426 = H_invDmH(5,1)*tt64-H_invDmH(5,2)*tt63+H_invDmH(5,3)*tt62
tt427 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt426*tt26*tt20-lam(1,1)*tt&
&65*tt85*tt16*tt20+(mu(1,1)*(2*tt65*tt85*tt16-2*tt426*tt26))/2.0E+&
&0+lam(1,1)*tt65*tt85*tt16)
tt428 = (H_invDmH(3,2)*H_invDmH(5,1)-H_invDmH(3,1)*H_invDmH(5,2))&
&*tt12-(H_invDmH(3,3)*H_invDmH(5,1)-H_invDmH(3,1)*H_invDmH(5,3))*t&
&t13+(H_invDmH(3,3)*H_invDmH(5,2)-H_invDmH(3,2)*H_invDmH(5,3))*tt1&
&4
tt429 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt428*tt26*tt20-lam(1,1)*tt&
&65*tt93*tt16*tt20+(mu(1,1)*(2*tt65*tt93*tt16-2*tt428*tt26))/2.0E+&
&0+lam(1,1)*tt65*tt93*tt16)
tt430 = detDmH(1,1)*gw(1,1)*((-lam(1,1)*tt65*tt99*tt16*tt20)+(mu(&
&1,1)*(2*tt65*tt99*tt16+tt363+tt362+tt361))/2.0E+0+lam(1,1)*tt65*t&
&t99*tt16)
tt431 = H_invDmH(6,1)*tt64-H_invDmH(6,2)*tt63+H_invDmH(6,3)*tt62
tt432 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt431*tt26*tt20-lam(1,1)*tt&
&65*tt102*tt16*tt20+(mu(1,1)*(2*tt65*tt102*tt16-2*tt431*tt26))/2.0&
&E+0+lam(1,1)*tt65*tt102*tt16)
tt433 = (H_invDmH(3,2)*H_invDmH(6,1)-H_invDmH(3,1)*H_invDmH(6,2))&
&*tt12-(H_invDmH(3,3)*H_invDmH(6,1)-H_invDmH(3,1)*H_invDmH(6,3))*t&
&t13+(H_invDmH(3,3)*H_invDmH(6,2)-H_invDmH(3,2)*H_invDmH(6,3))*tt1&
&4
tt434 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt433*tt26*tt20-lam(1,1)*tt&
&65*tt110*tt16*tt20+(mu(1,1)*(2*tt65*tt110*tt16-2*tt433*tt26))/2.0&
&E+0+lam(1,1)*tt65*tt110*tt16)
tt435 = detDmH(1,1)*gw(1,1)*((-lam(1,1)*tt65*tt116*tt16*tt20)+(mu&
&(1,1)*(2*tt65*tt116*tt16+tt371+tt370+tt369))/2.0E+0+lam(1,1)*tt65&
&*tt116*tt16)
tt436 = H_invDmH(7,1)*tt64-H_invDmH(7,2)*tt63+H_invDmH(7,3)*tt62
tt437 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt436*tt26*tt20-lam(1,1)*tt&
&65*tt119*tt16*tt20+(mu(1,1)*(2*tt65*tt119*tt16-2*tt436*tt26))/2.0&
&E+0+lam(1,1)*tt65*tt119*tt16)
tt438 = (H_invDmH(3,2)*H_invDmH(7,1)-H_invDmH(3,1)*H_invDmH(7,2))&
&*tt12-(H_invDmH(3,3)*H_invDmH(7,1)-H_invDmH(3,1)*H_invDmH(7,3))*t&
&t13+(H_invDmH(3,3)*H_invDmH(7,2)-H_invDmH(3,2)*H_invDmH(7,3))*tt1&
&4
tt439 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt438*tt26*tt20-lam(1,1)*tt&
&65*tt127*tt16*tt20+(mu(1,1)*(2*tt65*tt127*tt16-2*tt438*tt26))/2.0&
&E+0+lam(1,1)*tt65*tt127*tt16)
tt440 = detDmH(1,1)*gw(1,1)*((-lam(1,1)*tt65*tt133*tt16*tt20)+(mu&
&(1,1)*(2*tt65*tt133*tt16+tt379+tt378+tt377))/2.0E+0+lam(1,1)*tt65&
&*tt133*tt16)
tt441 = H_invDmH(8,1)*tt64-H_invDmH(8,2)*tt63+tt62*H_invDmH(8,3)
tt442 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt441*tt26*tt20-lam(1,1)*tt&
&65*tt136*tt16*tt20+(mu(1,1)*(2*tt65*tt136*tt16-2*tt441*tt26))/2.0&
&E+0+lam(1,1)*tt65*tt136*tt16)
tt443 = tt14*(H_invDmH(3,3)*H_invDmH(8,2)-H_invDmH(3,2)*H_invDmH(&
&8,3))-tt13*(H_invDmH(3,3)*H_invDmH(8,1)-H_invDmH(3,1)*H_invDmH(8,&
&3))+(H_invDmH(3,2)*H_invDmH(8,1)-H_invDmH(3,1)*H_invDmH(8,2))*tt1&
&2
tt444 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt443*tt26*tt20-lam(1,1)*tt&
&65*tt144*tt16*tt20+(mu(1,1)*(2*tt65*tt144*tt16-2*tt443*tt26))/2.0&
&E+0+lam(1,1)*tt65*tt144*tt16)
tt445 = detDmH(1,1)*gw(1,1)*((-lam(1,1)*tt65*tt150*tt16*tt20)+(mu&
&(1,1)*(2*tt65*tt150*tt16+tt387+tt386+tt385))/2.0E+0+lam(1,1)*tt65&
&*tt150*tt16)
tt446 = tt68**2
tt447 = 2*H_invDmH(4,1)**2
tt448 = 2*H_invDmH(4,2)**2
tt449 = 2*H_invDmH(4,3)**2
tt450 = H_invDmH(4,1)*tt75-H_invDmH(4,2)*tt74+H_invDmH(4,3)*tt73
tt451 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt450*tt26*tt20-lam(1,1)*tt&
&76*tt68*tt16*tt20+(mu(1,1)*(2*tt76*tt68*tt16-2*tt450*tt26))/2.0E+&
&0+lam(1,1)*tt76*tt68*tt16)
tt452 = H_invDmH(4,1)*tt81-H_invDmH(4,2)*tt80+H_invDmH(4,3)*tt79
tt453 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt452*tt26*tt20-lam(1,1)*tt&
&82*tt68*tt16*tt20+(mu(1,1)*(2*tt82*tt68*tt16-2*tt452*tt26))/2.0E+&
&0+lam(1,1)*tt82*tt68*tt16)
tt454 = 2*H_invDmH(4,1)*H_invDmH(5,1)
tt455 = 2*H_invDmH(4,2)*H_invDmH(5,2)
tt456 = 2*H_invDmH(4,3)*H_invDmH(5,3)
tt457 = detDmH(1,1)*gw(1,1)*((-lam(1,1)*tt68*tt85*tt16*tt20)+(mu(&
&1,1)*(2*tt68*tt85*tt16+tt456+tt455+tt454))/2.0E+0+lam(1,1)*tt68*t&
&t85*tt16)
tt458 = H_invDmH(4,1)*tt92-H_invDmH(4,2)*tt91+H_invDmH(4,3)*tt90
tt459 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt458*tt26*tt20-lam(1,1)*tt&
&93*tt68*tt16*tt20+(mu(1,1)*(2*tt93*tt68*tt16-2*tt458*tt26))/2.0E+&
&0+lam(1,1)*tt93*tt68*tt16)
tt460 = H_invDmH(4,1)*tt98-H_invDmH(4,2)*tt97+H_invDmH(4,3)*tt96
tt461 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt460*tt26*tt20-lam(1,1)*tt&
&99*tt68*tt16*tt20+(mu(1,1)*(2*tt99*tt68*tt16-2*tt460*tt26))/2.0E+&
&0+lam(1,1)*tt99*tt68*tt16)
tt462 = 2*H_invDmH(4,1)*H_invDmH(6,1)
tt463 = 2*H_invDmH(4,2)*H_invDmH(6,2)
tt464 = 2*H_invDmH(4,3)*H_invDmH(6,3)
tt465 = detDmH(1,1)*gw(1,1)*((-lam(1,1)*tt68*tt102*tt16*tt20)+(mu&
&(1,1)*(2*tt68*tt102*tt16+tt464+tt463+tt462))/2.0E+0+lam(1,1)*tt68&
&*tt102*tt16)
tt466 = H_invDmH(4,1)*tt109-H_invDmH(4,2)*tt108+H_invDmH(4,3)*tt1&
&07
tt467 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt466*tt26*tt20-lam(1,1)*tt&
&110*tt68*tt16*tt20+(mu(1,1)*(2*tt110*tt68*tt16-2*tt466*tt26))/2.0&
&E+0+lam(1,1)*tt110*tt68*tt16)
tt468 = H_invDmH(4,1)*tt115-H_invDmH(4,2)*tt114+H_invDmH(4,3)*tt1&
&13
tt469 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt468*tt26*tt20-lam(1,1)*tt&
&116*tt68*tt16*tt20+(mu(1,1)*(2*tt116*tt68*tt16-2*tt468*tt26))/2.0&
&E+0+lam(1,1)*tt116*tt68*tt16)
tt470 = 2*H_invDmH(4,1)*H_invDmH(7,1)
tt471 = 2*H_invDmH(4,2)*H_invDmH(7,2)
tt472 = 2*H_invDmH(4,3)*H_invDmH(7,3)
tt473 = detDmH(1,1)*gw(1,1)*((-lam(1,1)*tt68*tt119*tt16*tt20)+(mu&
&(1,1)*(2*tt68*tt119*tt16+tt472+tt471+tt470))/2.0E+0+lam(1,1)*tt68&
&*tt119*tt16)
tt474 = H_invDmH(4,1)*tt126-H_invDmH(4,2)*tt125+H_invDmH(4,3)*tt1&
&24
tt475 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt474*tt26*tt20-lam(1,1)*tt&
&127*tt68*tt16*tt20+(mu(1,1)*(2*tt127*tt68*tt16-2*tt474*tt26))/2.0&
&E+0+lam(1,1)*tt127*tt68*tt16)
tt476 = H_invDmH(4,1)*tt132-H_invDmH(4,2)*tt131+H_invDmH(4,3)*tt1&
&30
tt477 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt476*tt26*tt20-lam(1,1)*tt&
&133*tt68*tt16*tt20+(mu(1,1)*(2*tt133*tt68*tt16-2*tt476*tt26))/2.0&
&E+0+lam(1,1)*tt133*tt68*tt16)
tt478 = 2*H_invDmH(4,1)*H_invDmH(8,1)
tt479 = 2*H_invDmH(4,2)*H_invDmH(8,2)
tt480 = 2*H_invDmH(4,3)*H_invDmH(8,3)
tt481 = detDmH(1,1)*gw(1,1)*((-lam(1,1)*tt68*tt136*tt16*tt20)+(mu&
&(1,1)*(2*tt68*tt136*tt16+tt480+tt479+tt478))/2.0E+0+lam(1,1)*tt68&
&*tt136*tt16)
tt482 = H_invDmH(4,1)*tt143-H_invDmH(4,2)*tt142+H_invDmH(4,3)*tt1&
&41
tt483 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt482*tt26*tt20-lam(1,1)*tt&
&144*tt68*tt16*tt20+(mu(1,1)*(2*tt144*tt68*tt16-2*tt482*tt26))/2.0&
&E+0+lam(1,1)*tt144*tt68*tt16)
tt484 = H_invDmH(4,1)*tt149-H_invDmH(4,2)*tt148+H_invDmH(4,3)*tt1&
&47
tt485 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt484*tt26*tt20-lam(1,1)*tt&
&150*tt68*tt16*tt20+(mu(1,1)*(2*tt150*tt68*tt16-2*tt484*tt26))/2.0&
&E+0+lam(1,1)*tt150*tt68*tt16)
tt486 = tt76**2
tt487 = detDmH(1,1)*gw(1,1)*((-lam(1,1)*tt82*tt76*tt16*tt20)+mu(1&
&,1)*tt82*tt76*tt16+lam(1,1)*tt82*tt76*tt16)
tt488 = H_invDmH(5,1)*tt75-H_invDmH(5,2)*tt74+H_invDmH(5,3)*tt73
tt489 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt488*tt26*tt20-lam(1,1)*tt&
&76*tt85*tt16*tt20+(mu(1,1)*(2*tt76*tt85*tt16-2*tt488*tt26))/2.0E+&
&0+lam(1,1)*tt76*tt85*tt16)
tt490 = detDmH(1,1)*gw(1,1)*((-lam(1,1)*tt76*tt93*tt16*tt20)+(mu(&
&1,1)*(2*tt76*tt93*tt16+tt456+tt455+tt454))/2.0E+0+lam(1,1)*tt76*t&
&t93*tt16)
tt491 = (H_invDmH(4,1)*H_invDmH(5,2)-H_invDmH(4,2)*H_invDmH(5,1))&
&*tt12-(H_invDmH(4,1)*H_invDmH(5,3)-H_invDmH(4,3)*H_invDmH(5,1))*t&
&t13+(H_invDmH(4,2)*H_invDmH(5,3)-H_invDmH(4,3)*H_invDmH(5,2))*tt1&
&4
tt492 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt491*tt26*tt20-lam(1,1)*tt&
&99*tt76*tt16*tt20+(mu(1,1)*(2*tt99*tt76*tt16-2*tt491*tt26))/2.0E+&
&0+lam(1,1)*tt99*tt76*tt16)
tt493 = H_invDmH(6,1)*tt75-H_invDmH(6,2)*tt74+H_invDmH(6,3)*tt73
tt494 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt493*tt26*tt20-lam(1,1)*tt&
&76*tt102*tt16*tt20+(mu(1,1)*(2*tt76*tt102*tt16-2*tt493*tt26))/2.0&
&E+0+lam(1,1)*tt76*tt102*tt16)
tt495 = detDmH(1,1)*gw(1,1)*((-lam(1,1)*tt76*tt110*tt16*tt20)+(mu&
&(1,1)*(2*tt76*tt110*tt16+tt464+tt463+tt462))/2.0E+0+lam(1,1)*tt76&
&*tt110*tt16)
tt496 = (H_invDmH(4,1)*H_invDmH(6,2)-H_invDmH(4,2)*H_invDmH(6,1))&
&*tt12-(H_invDmH(4,1)*H_invDmH(6,3)-H_invDmH(4,3)*H_invDmH(6,1))*t&
&t13+(H_invDmH(4,2)*H_invDmH(6,3)-H_invDmH(4,3)*H_invDmH(6,2))*tt1&
&4
tt497 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt496*tt26*tt20-lam(1,1)*tt&
&116*tt76*tt16*tt20+(mu(1,1)*(2*tt116*tt76*tt16-2*tt496*tt26))/2.0&
&E+0+lam(1,1)*tt116*tt76*tt16)
tt498 = H_invDmH(7,1)*tt75-H_invDmH(7,2)*tt74+H_invDmH(7,3)*tt73
tt499 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt498*tt26*tt20-lam(1,1)*tt&
&76*tt119*tt16*tt20+(mu(1,1)*(2*tt76*tt119*tt16-2*tt498*tt26))/2.0&
&E+0+lam(1,1)*tt76*tt119*tt16)
tt500 = detDmH(1,1)*gw(1,1)*((-lam(1,1)*tt76*tt127*tt16*tt20)+(mu&
&(1,1)*(2*tt76*tt127*tt16+tt472+tt471+tt470))/2.0E+0+lam(1,1)*tt76&
&*tt127*tt16)
tt501 = (H_invDmH(4,1)*H_invDmH(7,2)-H_invDmH(4,2)*H_invDmH(7,1))&
&*tt12-(H_invDmH(4,1)*H_invDmH(7,3)-H_invDmH(4,3)*H_invDmH(7,1))*t&
&t13+(H_invDmH(4,2)*H_invDmH(7,3)-H_invDmH(4,3)*H_invDmH(7,2))*tt1&
&4
tt502 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt501*tt26*tt20-lam(1,1)*tt&
&133*tt76*tt16*tt20+(mu(1,1)*(2*tt133*tt76*tt16-2*tt501*tt26))/2.0&
&E+0+lam(1,1)*tt133*tt76*tt16)
tt503 = H_invDmH(8,1)*tt75-H_invDmH(8,2)*tt74+tt73*H_invDmH(8,3)
tt504 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt503*tt26*tt20-lam(1,1)*tt&
&76*tt136*tt16*tt20+(mu(1,1)*(2*tt76*tt136*tt16-2*tt503*tt26))/2.0&
&E+0+lam(1,1)*tt76*tt136*tt16)
tt505 = detDmH(1,1)*gw(1,1)*((-lam(1,1)*tt76*tt144*tt16*tt20)+(mu&
&(1,1)*(2*tt76*tt144*tt16+tt480+tt479+tt478))/2.0E+0+lam(1,1)*tt76&
&*tt144*tt16)
tt506 = tt14*(H_invDmH(4,2)*H_invDmH(8,3)-H_invDmH(4,3)*H_invDmH(&
&8,2))-tt13*(H_invDmH(4,1)*H_invDmH(8,3)-H_invDmH(4,3)*H_invDmH(8,&
&1))+(H_invDmH(4,1)*H_invDmH(8,2)-H_invDmH(4,2)*H_invDmH(8,1))*tt1&
&2
tt507 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt506*tt26*tt20-lam(1,1)*tt&
&150*tt76*tt16*tt20+(mu(1,1)*(2*tt150*tt76*tt16-2*tt506*tt26))/2.0&
&E+0+lam(1,1)*tt150*tt76*tt16)
tt508 = tt82**2
tt509 = H_invDmH(5,1)*tt81-H_invDmH(5,2)*tt80+H_invDmH(5,3)*tt79
tt510 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt509*tt26*tt20-lam(1,1)*tt&
&82*tt85*tt16*tt20+(mu(1,1)*(2*tt82*tt85*tt16-2*tt509*tt26))/2.0E+&
&0+lam(1,1)*tt82*tt85*tt16)
tt511 = (H_invDmH(4,2)*H_invDmH(5,1)-H_invDmH(4,1)*H_invDmH(5,2))&
&*tt12-(H_invDmH(4,3)*H_invDmH(5,1)-H_invDmH(4,1)*H_invDmH(5,3))*t&
&t13+(H_invDmH(4,3)*H_invDmH(5,2)-H_invDmH(4,2)*H_invDmH(5,3))*tt1&
&4
tt512 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt511*tt26*tt20-lam(1,1)*tt&
&82*tt93*tt16*tt20+(mu(1,1)*(2*tt82*tt93*tt16-2*tt511*tt26))/2.0E+&
&0+lam(1,1)*tt82*tt93*tt16)
tt513 = detDmH(1,1)*gw(1,1)*((-lam(1,1)*tt82*tt99*tt16*tt20)+(mu(&
&1,1)*(2*tt82*tt99*tt16+tt456+tt455+tt454))/2.0E+0+lam(1,1)*tt82*t&
&t99*tt16)
tt514 = H_invDmH(6,1)*tt81-H_invDmH(6,2)*tt80+H_invDmH(6,3)*tt79
tt515 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt514*tt26*tt20-lam(1,1)*tt&
&82*tt102*tt16*tt20+(mu(1,1)*(2*tt82*tt102*tt16-2*tt514*tt26))/2.0&
&E+0+lam(1,1)*tt82*tt102*tt16)
tt516 = (H_invDmH(4,2)*H_invDmH(6,1)-H_invDmH(4,1)*H_invDmH(6,2))&
&*tt12-(H_invDmH(4,3)*H_invDmH(6,1)-H_invDmH(4,1)*H_invDmH(6,3))*t&
&t13+(H_invDmH(4,3)*H_invDmH(6,2)-H_invDmH(4,2)*H_invDmH(6,3))*tt1&
&4
tt517 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt516*tt26*tt20-lam(1,1)*tt&
&82*tt110*tt16*tt20+(mu(1,1)*(2*tt82*tt110*tt16-2*tt516*tt26))/2.0&
&E+0+lam(1,1)*tt82*tt110*tt16)
tt518 = detDmH(1,1)*gw(1,1)*((-lam(1,1)*tt82*tt116*tt16*tt20)+(mu&
&(1,1)*(2*tt82*tt116*tt16+tt464+tt463+tt462))/2.0E+0+lam(1,1)*tt82&
&*tt116*tt16)
tt519 = H_invDmH(7,1)*tt81-H_invDmH(7,2)*tt80+H_invDmH(7,3)*tt79
tt520 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt519*tt26*tt20-lam(1,1)*tt&
&82*tt119*tt16*tt20+(mu(1,1)*(2*tt82*tt119*tt16-2*tt519*tt26))/2.0&
&E+0+lam(1,1)*tt82*tt119*tt16)
tt521 = (H_invDmH(4,2)*H_invDmH(7,1)-H_invDmH(4,1)*H_invDmH(7,2))&
&*tt12-(H_invDmH(4,3)*H_invDmH(7,1)-H_invDmH(4,1)*H_invDmH(7,3))*t&
&t13+(H_invDmH(4,3)*H_invDmH(7,2)-H_invDmH(4,2)*H_invDmH(7,3))*tt1&
&4
tt522 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt521*tt26*tt20-lam(1,1)*tt&
&82*tt127*tt16*tt20+(mu(1,1)*(2*tt82*tt127*tt16-2*tt521*tt26))/2.0&
&E+0+lam(1,1)*tt82*tt127*tt16)
tt523 = detDmH(1,1)*gw(1,1)*((-lam(1,1)*tt82*tt133*tt16*tt20)+(mu&
&(1,1)*(2*tt82*tt133*tt16+tt472+tt471+tt470))/2.0E+0+lam(1,1)*tt82&
&*tt133*tt16)
tt524 = H_invDmH(8,1)*tt81-H_invDmH(8,2)*tt80+tt79*H_invDmH(8,3)
tt525 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt524*tt26*tt20-lam(1,1)*tt&
&82*tt136*tt16*tt20+(mu(1,1)*(2*tt82*tt136*tt16-2*tt524*tt26))/2.0&
&E+0+lam(1,1)*tt82*tt136*tt16)
tt526 = tt14*(H_invDmH(4,3)*H_invDmH(8,2)-H_invDmH(4,2)*H_invDmH(&
&8,3))-tt13*(H_invDmH(4,3)*H_invDmH(8,1)-H_invDmH(4,1)*H_invDmH(8,&
&3))+(H_invDmH(4,2)*H_invDmH(8,1)-H_invDmH(4,1)*H_invDmH(8,2))*tt1&
&2
tt527 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt526*tt26*tt20-lam(1,1)*tt&
&82*tt144*tt16*tt20+(mu(1,1)*(2*tt82*tt144*tt16-2*tt526*tt26))/2.0&
&E+0+lam(1,1)*tt82*tt144*tt16)
tt528 = detDmH(1,1)*gw(1,1)*((-lam(1,1)*tt82*tt150*tt16*tt20)+(mu&
&(1,1)*(2*tt82*tt150*tt16+tt480+tt479+tt478))/2.0E+0+lam(1,1)*tt82&
&*tt150*tt16)
tt529 = tt85**2
tt530 = 2*H_invDmH(5,1)**2
tt531 = 2*H_invDmH(5,2)**2
tt532 = 2*H_invDmH(5,3)**2
tt533 = H_invDmH(5,1)*tt92-H_invDmH(5,2)*tt91+H_invDmH(5,3)*tt90
tt534 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt533*tt26*tt20-lam(1,1)*tt&
&93*tt85*tt16*tt20+(mu(1,1)*(2*tt93*tt85*tt16-2*tt533*tt26))/2.0E+&
&0+lam(1,1)*tt93*tt85*tt16)
tt535 = H_invDmH(5,1)*tt98-H_invDmH(5,2)*tt97+H_invDmH(5,3)*tt96
tt536 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt535*tt26*tt20-lam(1,1)*tt&
&99*tt85*tt16*tt20+(mu(1,1)*(2*tt99*tt85*tt16-2*tt535*tt26))/2.0E+&
&0+lam(1,1)*tt99*tt85*tt16)
tt537 = 2*H_invDmH(5,1)*H_invDmH(6,1)
tt538 = 2*H_invDmH(5,2)*H_invDmH(6,2)
tt539 = 2*H_invDmH(5,3)*H_invDmH(6,3)
tt540 = detDmH(1,1)*gw(1,1)*((-lam(1,1)*tt85*tt102*tt16*tt20)+(mu&
&(1,1)*(2*tt85*tt102*tt16+tt539+tt538+tt537))/2.0E+0+lam(1,1)*tt85&
&*tt102*tt16)
tt541 = H_invDmH(5,1)*tt109-H_invDmH(5,2)*tt108+H_invDmH(5,3)*tt1&
&07
tt542 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt541*tt26*tt20-lam(1,1)*tt&
&110*tt85*tt16*tt20+(mu(1,1)*(2*tt110*tt85*tt16-2*tt541*tt26))/2.0&
&E+0+lam(1,1)*tt110*tt85*tt16)
tt543 = H_invDmH(5,1)*tt115-H_invDmH(5,2)*tt114+H_invDmH(5,3)*tt1&
&13
tt544 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt543*tt26*tt20-lam(1,1)*tt&
&116*tt85*tt16*tt20+(mu(1,1)*(2*tt116*tt85*tt16-2*tt543*tt26))/2.0&
&E+0+lam(1,1)*tt116*tt85*tt16)
tt545 = 2*H_invDmH(5,1)*H_invDmH(7,1)
tt546 = 2*H_invDmH(5,2)*H_invDmH(7,2)
tt547 = 2*H_invDmH(5,3)*H_invDmH(7,3)
tt548 = detDmH(1,1)*gw(1,1)*((-lam(1,1)*tt85*tt119*tt16*tt20)+(mu&
&(1,1)*(2*tt85*tt119*tt16+tt547+tt546+tt545))/2.0E+0+lam(1,1)*tt85&
&*tt119*tt16)
tt549 = H_invDmH(5,1)*tt126-H_invDmH(5,2)*tt125+H_invDmH(5,3)*tt1&
&24
tt550 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt549*tt26*tt20-lam(1,1)*tt&
&127*tt85*tt16*tt20+(mu(1,1)*(2*tt127*tt85*tt16-2*tt549*tt26))/2.0&
&E+0+lam(1,1)*tt127*tt85*tt16)
tt551 = H_invDmH(5,1)*tt132-H_invDmH(5,2)*tt131+H_invDmH(5,3)*tt1&
&30
tt552 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt551*tt26*tt20-lam(1,1)*tt&
&133*tt85*tt16*tt20+(mu(1,1)*(2*tt133*tt85*tt16-2*tt551*tt26))/2.0&
&E+0+lam(1,1)*tt133*tt85*tt16)
tt553 = 2*H_invDmH(5,1)*H_invDmH(8,1)
tt554 = 2*H_invDmH(5,2)*H_invDmH(8,2)
tt555 = 2*H_invDmH(5,3)*H_invDmH(8,3)
tt556 = detDmH(1,1)*gw(1,1)*((-lam(1,1)*tt85*tt136*tt16*tt20)+(mu&
&(1,1)*(2*tt85*tt136*tt16+tt555+tt554+tt553))/2.0E+0+lam(1,1)*tt85&
&*tt136*tt16)
tt557 = H_invDmH(5,1)*tt143-H_invDmH(5,2)*tt142+H_invDmH(5,3)*tt1&
&41
tt558 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt557*tt26*tt20-lam(1,1)*tt&
&144*tt85*tt16*tt20+(mu(1,1)*(2*tt144*tt85*tt16-2*tt557*tt26))/2.0&
&E+0+lam(1,1)*tt144*tt85*tt16)
tt559 = H_invDmH(5,1)*tt149-H_invDmH(5,2)*tt148+H_invDmH(5,3)*tt1&
&47
tt560 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt559*tt26*tt20-lam(1,1)*tt&
&150*tt85*tt16*tt20+(mu(1,1)*(2*tt150*tt85*tt16-2*tt559*tt26))/2.0&
&E+0+lam(1,1)*tt150*tt85*tt16)
tt561 = tt93**2
tt562 = detDmH(1,1)*gw(1,1)*((-lam(1,1)*tt99*tt93*tt16*tt20)+mu(1&
&,1)*tt99*tt93*tt16+lam(1,1)*tt99*tt93*tt16)
tt563 = H_invDmH(6,1)*tt92-H_invDmH(6,2)*tt91+H_invDmH(6,3)*tt90
tt564 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt563*tt26*tt20-lam(1,1)*tt&
&93*tt102*tt16*tt20+(mu(1,1)*(2*tt93*tt102*tt16-2*tt563*tt26))/2.0&
&E+0+lam(1,1)*tt93*tt102*tt16)
tt565 = detDmH(1,1)*gw(1,1)*((-lam(1,1)*tt93*tt110*tt16*tt20)+(mu&
&(1,1)*(2*tt93*tt110*tt16+tt539+tt538+tt537))/2.0E+0+lam(1,1)*tt93&
&*tt110*tt16)
tt566 = (H_invDmH(5,1)*H_invDmH(6,2)-H_invDmH(5,2)*H_invDmH(6,1))&
&*tt12-(H_invDmH(5,1)*H_invDmH(6,3)-H_invDmH(5,3)*H_invDmH(6,1))*t&
&t13+(H_invDmH(5,2)*H_invDmH(6,3)-H_invDmH(5,3)*H_invDmH(6,2))*tt1&
&4
tt567 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt566*tt26*tt20-lam(1,1)*tt&
&116*tt93*tt16*tt20+(mu(1,1)*(2*tt116*tt93*tt16-2*tt566*tt26))/2.0&
&E+0+lam(1,1)*tt116*tt93*tt16)
tt568 = H_invDmH(7,1)*tt92-H_invDmH(7,2)*tt91+H_invDmH(7,3)*tt90
tt569 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt568*tt26*tt20-lam(1,1)*tt&
&93*tt119*tt16*tt20+(mu(1,1)*(2*tt93*tt119*tt16-2*tt568*tt26))/2.0&
&E+0+lam(1,1)*tt93*tt119*tt16)
tt570 = detDmH(1,1)*gw(1,1)*((-lam(1,1)*tt93*tt127*tt16*tt20)+(mu&
&(1,1)*(2*tt93*tt127*tt16+tt547+tt546+tt545))/2.0E+0+lam(1,1)*tt93&
&*tt127*tt16)
tt571 = (H_invDmH(5,1)*H_invDmH(7,2)-H_invDmH(5,2)*H_invDmH(7,1))&
&*tt12-(H_invDmH(5,1)*H_invDmH(7,3)-H_invDmH(5,3)*H_invDmH(7,1))*t&
&t13+(H_invDmH(5,2)*H_invDmH(7,3)-H_invDmH(5,3)*H_invDmH(7,2))*tt1&
&4
tt572 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt571*tt26*tt20-lam(1,1)*tt&
&133*tt93*tt16*tt20+(mu(1,1)*(2*tt133*tt93*tt16-2*tt571*tt26))/2.0&
&E+0+lam(1,1)*tt133*tt93*tt16)
tt573 = H_invDmH(8,1)*tt92-H_invDmH(8,2)*tt91+tt90*H_invDmH(8,3)
tt574 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt573*tt26*tt20-lam(1,1)*tt&
&93*tt136*tt16*tt20+(mu(1,1)*(2*tt93*tt136*tt16-2*tt573*tt26))/2.0&
&E+0+lam(1,1)*tt93*tt136*tt16)
tt575 = detDmH(1,1)*gw(1,1)*((-lam(1,1)*tt93*tt144*tt16*tt20)+(mu&
&(1,1)*(2*tt93*tt144*tt16+tt555+tt554+tt553))/2.0E+0+lam(1,1)*tt93&
&*tt144*tt16)
tt576 = tt14*(H_invDmH(5,2)*H_invDmH(8,3)-H_invDmH(5,3)*H_invDmH(&
&8,2))-tt13*(H_invDmH(5,1)*H_invDmH(8,3)-H_invDmH(5,3)*H_invDmH(8,&
&1))+(H_invDmH(5,1)*H_invDmH(8,2)-H_invDmH(5,2)*H_invDmH(8,1))*tt1&
&2
tt577 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt576*tt26*tt20-lam(1,1)*tt&
&150*tt93*tt16*tt20+(mu(1,1)*(2*tt150*tt93*tt16-2*tt576*tt26))/2.0&
&E+0+lam(1,1)*tt150*tt93*tt16)
tt578 = tt99**2
tt579 = H_invDmH(6,1)*tt98-H_invDmH(6,2)*tt97+H_invDmH(6,3)*tt96
tt580 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt579*tt26*tt20-lam(1,1)*tt&
&99*tt102*tt16*tt20+(mu(1,1)*(2*tt99*tt102*tt16-2*tt579*tt26))/2.0&
&E+0+lam(1,1)*tt99*tt102*tt16)
tt581 = (H_invDmH(5,2)*H_invDmH(6,1)-H_invDmH(5,1)*H_invDmH(6,2))&
&*tt12-(H_invDmH(5,3)*H_invDmH(6,1)-H_invDmH(5,1)*H_invDmH(6,3))*t&
&t13+(H_invDmH(5,3)*H_invDmH(6,2)-H_invDmH(5,2)*H_invDmH(6,3))*tt1&
&4
tt582 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt581*tt26*tt20-lam(1,1)*tt&
&99*tt110*tt16*tt20+(mu(1,1)*(2*tt99*tt110*tt16-2*tt581*tt26))/2.0&
&E+0+lam(1,1)*tt99*tt110*tt16)
tt583 = detDmH(1,1)*gw(1,1)*((-lam(1,1)*tt99*tt116*tt16*tt20)+(mu&
&(1,1)*(2*tt99*tt116*tt16+tt539+tt538+tt537))/2.0E+0+lam(1,1)*tt99&
&*tt116*tt16)
tt584 = H_invDmH(7,1)*tt98-H_invDmH(7,2)*tt97+H_invDmH(7,3)*tt96
tt585 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt584*tt26*tt20-lam(1,1)*tt&
&99*tt119*tt16*tt20+(mu(1,1)*(2*tt99*tt119*tt16-2*tt584*tt26))/2.0&
&E+0+lam(1,1)*tt99*tt119*tt16)
tt586 = (H_invDmH(5,2)*H_invDmH(7,1)-H_invDmH(5,1)*H_invDmH(7,2))&
&*tt12-(H_invDmH(5,3)*H_invDmH(7,1)-H_invDmH(5,1)*H_invDmH(7,3))*t&
&t13+(H_invDmH(5,3)*H_invDmH(7,2)-H_invDmH(5,2)*H_invDmH(7,3))*tt1&
&4
tt587 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt586*tt26*tt20-lam(1,1)*tt&
&99*tt127*tt16*tt20+(mu(1,1)*(2*tt99*tt127*tt16-2*tt586*tt26))/2.0&
&E+0+lam(1,1)*tt99*tt127*tt16)
tt588 = detDmH(1,1)*gw(1,1)*((-lam(1,1)*tt99*tt133*tt16*tt20)+(mu&
&(1,1)*(2*tt99*tt133*tt16+tt547+tt546+tt545))/2.0E+0+lam(1,1)*tt99&
&*tt133*tt16)
tt589 = H_invDmH(8,1)*tt98-H_invDmH(8,2)*tt97+tt96*H_invDmH(8,3)
tt590 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt589*tt26*tt20-lam(1,1)*tt&
&99*tt136*tt16*tt20+(mu(1,1)*(2*tt99*tt136*tt16-2*tt589*tt26))/2.0&
&E+0+lam(1,1)*tt99*tt136*tt16)
tt591 = tt14*(H_invDmH(5,3)*H_invDmH(8,2)-H_invDmH(5,2)*H_invDmH(&
&8,3))-tt13*(H_invDmH(5,3)*H_invDmH(8,1)-H_invDmH(5,1)*H_invDmH(8,&
&3))+(H_invDmH(5,2)*H_invDmH(8,1)-H_invDmH(5,1)*H_invDmH(8,2))*tt1&
&2
tt592 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt591*tt26*tt20-lam(1,1)*tt&
&99*tt144*tt16*tt20+(mu(1,1)*(2*tt99*tt144*tt16-2*tt591*tt26))/2.0&
&E+0+lam(1,1)*tt99*tt144*tt16)
tt593 = detDmH(1,1)*gw(1,1)*((-lam(1,1)*tt99*tt150*tt16*tt20)+(mu&
&(1,1)*(2*tt99*tt150*tt16+tt555+tt554+tt553))/2.0E+0+lam(1,1)*tt99&
&*tt150*tt16)
tt594 = tt102**2
tt595 = 2*H_invDmH(6,1)**2
tt596 = 2*H_invDmH(6,2)**2
tt597 = 2*H_invDmH(6,3)**2
tt598 = H_invDmH(6,1)*tt109-H_invDmH(6,2)*tt108+H_invDmH(6,3)*tt1&
&07
tt599 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt598*tt26*tt20-lam(1,1)*tt&
&110*tt102*tt16*tt20+(mu(1,1)*(2*tt110*tt102*tt16-2*tt598*tt26))/2&
&.0E+0+lam(1,1)*tt110*tt102*tt16)
tt600 = H_invDmH(6,1)*tt115-H_invDmH(6,2)*tt114+H_invDmH(6,3)*tt1&
&13
tt601 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt600*tt26*tt20-lam(1,1)*tt&
&116*tt102*tt16*tt20+(mu(1,1)*(2*tt116*tt102*tt16-2*tt600*tt26))/2&
&.0E+0+lam(1,1)*tt116*tt102*tt16)
tt602 = 2*H_invDmH(6,1)*H_invDmH(7,1)
tt603 = 2*H_invDmH(6,2)*H_invDmH(7,2)
tt604 = 2*H_invDmH(6,3)*H_invDmH(7,3)
tt605 = detDmH(1,1)*gw(1,1)*((-lam(1,1)*tt102*tt119*tt16*tt20)+(m&
&u(1,1)*(2*tt102*tt119*tt16+tt604+tt603+tt602))/2.0E+0+lam(1,1)*tt&
&102*tt119*tt16)
tt606 = H_invDmH(6,1)*tt126-H_invDmH(6,2)*tt125+H_invDmH(6,3)*tt1&
&24
tt607 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt606*tt26*tt20-lam(1,1)*tt&
&127*tt102*tt16*tt20+(mu(1,1)*(2*tt127*tt102*tt16-2*tt606*tt26))/2&
&.0E+0+lam(1,1)*tt127*tt102*tt16)
tt608 = H_invDmH(6,1)*tt132-H_invDmH(6,2)*tt131+H_invDmH(6,3)*tt1&
&30
tt609 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt608*tt26*tt20-lam(1,1)*tt&
&133*tt102*tt16*tt20+(mu(1,1)*(2*tt133*tt102*tt16-2*tt608*tt26))/2&
&.0E+0+lam(1,1)*tt133*tt102*tt16)
tt610 = 2*H_invDmH(6,1)*H_invDmH(8,1)
tt611 = 2*H_invDmH(6,2)*H_invDmH(8,2)
tt612 = 2*H_invDmH(6,3)*H_invDmH(8,3)
tt613 = detDmH(1,1)*gw(1,1)*((-lam(1,1)*tt102*tt136*tt16*tt20)+(m&
&u(1,1)*(2*tt102*tt136*tt16+tt612+tt611+tt610))/2.0E+0+lam(1,1)*tt&
&102*tt136*tt16)
tt614 = H_invDmH(6,1)*tt143-H_invDmH(6,2)*tt142+H_invDmH(6,3)*tt1&
&41
tt615 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt614*tt26*tt20-lam(1,1)*tt&
&144*tt102*tt16*tt20+(mu(1,1)*(2*tt144*tt102*tt16-2*tt614*tt26))/2&
&.0E+0+lam(1,1)*tt144*tt102*tt16)
tt616 = H_invDmH(6,1)*tt149-H_invDmH(6,2)*tt148+H_invDmH(6,3)*tt1&
&47
tt617 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt616*tt26*tt20-lam(1,1)*tt&
&150*tt102*tt16*tt20+(mu(1,1)*(2*tt150*tt102*tt16-2*tt616*tt26))/2&
&.0E+0+lam(1,1)*tt150*tt102*tt16)
tt618 = tt110**2
tt619 = detDmH(1,1)*gw(1,1)*((-lam(1,1)*tt116*tt110*tt16*tt20)+mu&
&(1,1)*tt116*tt110*tt16+lam(1,1)*tt116*tt110*tt16)
tt620 = H_invDmH(7,1)*tt109-H_invDmH(7,2)*tt108+H_invDmH(7,3)*tt1&
&07
tt621 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt620*tt26*tt20-lam(1,1)*tt&
&110*tt119*tt16*tt20+(mu(1,1)*(2*tt110*tt119*tt16-2*tt620*tt26))/2&
&.0E+0+lam(1,1)*tt110*tt119*tt16)
tt622 = detDmH(1,1)*gw(1,1)*((-lam(1,1)*tt110*tt127*tt16*tt20)+(m&
&u(1,1)*(2*tt110*tt127*tt16+tt604+tt603+tt602))/2.0E+0+lam(1,1)*tt&
&110*tt127*tt16)
tt623 = (H_invDmH(6,1)*H_invDmH(7,2)-H_invDmH(6,2)*H_invDmH(7,1))&
&*tt12-(H_invDmH(6,1)*H_invDmH(7,3)-H_invDmH(6,3)*H_invDmH(7,1))*t&
&t13+(H_invDmH(6,2)*H_invDmH(7,3)-H_invDmH(6,3)*H_invDmH(7,2))*tt1&
&4
tt624 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt623*tt26*tt20-lam(1,1)*tt&
&133*tt110*tt16*tt20+(mu(1,1)*(2*tt133*tt110*tt16-2*tt623*tt26))/2&
&.0E+0+lam(1,1)*tt133*tt110*tt16)
tt625 = H_invDmH(8,1)*tt109-H_invDmH(8,2)*tt108+tt107*H_invDmH(8,&
&3)
tt626 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt625*tt26*tt20-lam(1,1)*tt&
&110*tt136*tt16*tt20+(mu(1,1)*(2*tt110*tt136*tt16-2*tt625*tt26))/2&
&.0E+0+lam(1,1)*tt110*tt136*tt16)
tt627 = detDmH(1,1)*gw(1,1)*((-lam(1,1)*tt110*tt144*tt16*tt20)+(m&
&u(1,1)*(2*tt110*tt144*tt16+tt612+tt611+tt610))/2.0E+0+lam(1,1)*tt&
&110*tt144*tt16)
tt628 = tt14*(H_invDmH(6,2)*H_invDmH(8,3)-H_invDmH(6,3)*H_invDmH(&
&8,2))-tt13*(H_invDmH(6,1)*H_invDmH(8,3)-H_invDmH(6,3)*H_invDmH(8,&
&1))+(H_invDmH(6,1)*H_invDmH(8,2)-H_invDmH(6,2)*H_invDmH(8,1))*tt1&
&2
tt629 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt628*tt26*tt20-lam(1,1)*tt&
&150*tt110*tt16*tt20+(mu(1,1)*(2*tt150*tt110*tt16-2*tt628*tt26))/2&
&.0E+0+lam(1,1)*tt150*tt110*tt16)
tt630 = tt116**2
tt631 = H_invDmH(7,1)*tt115-H_invDmH(7,2)*tt114+H_invDmH(7,3)*tt1&
&13
tt632 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt631*tt26*tt20-lam(1,1)*tt&
&116*tt119*tt16*tt20+(mu(1,1)*(2*tt116*tt119*tt16-2*tt631*tt26))/2&
&.0E+0+lam(1,1)*tt116*tt119*tt16)
tt633 = (H_invDmH(6,2)*H_invDmH(7,1)-H_invDmH(6,1)*H_invDmH(7,2))&
&*tt12-(H_invDmH(6,3)*H_invDmH(7,1)-H_invDmH(6,1)*H_invDmH(7,3))*t&
&t13+(H_invDmH(6,3)*H_invDmH(7,2)-H_invDmH(6,2)*H_invDmH(7,3))*tt1&
&4
tt634 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt633*tt26*tt20-lam(1,1)*tt&
&116*tt127*tt16*tt20+(mu(1,1)*(2*tt116*tt127*tt16-2*tt633*tt26))/2&
&.0E+0+lam(1,1)*tt116*tt127*tt16)
tt635 = detDmH(1,1)*gw(1,1)*((-lam(1,1)*tt116*tt133*tt16*tt20)+(m&
&u(1,1)*(2*tt116*tt133*tt16+tt604+tt603+tt602))/2.0E+0+lam(1,1)*tt&
&116*tt133*tt16)
tt636 = H_invDmH(8,1)*tt115-H_invDmH(8,2)*tt114+tt113*H_invDmH(8,&
&3)
tt637 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt636*tt26*tt20-lam(1,1)*tt&
&116*tt136*tt16*tt20+(mu(1,1)*(2*tt116*tt136*tt16-2*tt636*tt26))/2&
&.0E+0+lam(1,1)*tt116*tt136*tt16)
tt638 = tt14*(H_invDmH(6,3)*H_invDmH(8,2)-H_invDmH(6,2)*H_invDmH(&
&8,3))-tt13*(H_invDmH(6,3)*H_invDmH(8,1)-H_invDmH(6,1)*H_invDmH(8,&
&3))+(H_invDmH(6,2)*H_invDmH(8,1)-H_invDmH(6,1)*H_invDmH(8,2))*tt1&
&2
tt639 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt638*tt26*tt20-lam(1,1)*tt&
&116*tt144*tt16*tt20+(mu(1,1)*(2*tt116*tt144*tt16-2*tt638*tt26))/2&
&.0E+0+lam(1,1)*tt116*tt144*tt16)
tt640 = detDmH(1,1)*gw(1,1)*((-lam(1,1)*tt116*tt150*tt16*tt20)+(m&
&u(1,1)*(2*tt116*tt150*tt16+tt612+tt611+tt610))/2.0E+0+lam(1,1)*tt&
&116*tt150*tt16)
tt641 = tt119**2
tt642 = 2*H_invDmH(7,1)**2
tt643 = 2*H_invDmH(7,2)**2
tt644 = 2*H_invDmH(7,3)**2
tt645 = H_invDmH(7,1)*tt126-H_invDmH(7,2)*tt125+H_invDmH(7,3)*tt1&
&24
tt646 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt645*tt26*tt20-lam(1,1)*tt&
&127*tt119*tt16*tt20+(mu(1,1)*(2*tt127*tt119*tt16-2*tt645*tt26))/2&
&.0E+0+lam(1,1)*tt127*tt119*tt16)
tt647 = H_invDmH(7,1)*tt132-H_invDmH(7,2)*tt131+H_invDmH(7,3)*tt1&
&30
tt648 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt647*tt26*tt20-lam(1,1)*tt&
&133*tt119*tt16*tt20+(mu(1,1)*(2*tt133*tt119*tt16-2*tt647*tt26))/2&
&.0E+0+lam(1,1)*tt133*tt119*tt16)
tt649 = 2*H_invDmH(7,1)*H_invDmH(8,1)
tt650 = 2*H_invDmH(7,2)*H_invDmH(8,2)
tt651 = 2*H_invDmH(7,3)*H_invDmH(8,3)
tt652 = detDmH(1,1)*gw(1,1)*((-lam(1,1)*tt119*tt136*tt16*tt20)+(m&
&u(1,1)*(2*tt119*tt136*tt16+tt651+tt650+tt649))/2.0E+0+lam(1,1)*tt&
&119*tt136*tt16)
tt653 = H_invDmH(7,1)*tt143-H_invDmH(7,2)*tt142+H_invDmH(7,3)*tt1&
&41
tt654 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt653*tt26*tt20-lam(1,1)*tt&
&144*tt119*tt16*tt20+(mu(1,1)*(2*tt144*tt119*tt16-2*tt653*tt26))/2&
&.0E+0+lam(1,1)*tt144*tt119*tt16)
tt655 = H_invDmH(7,1)*tt149-H_invDmH(7,2)*tt148+H_invDmH(7,3)*tt1&
&47
tt656 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt655*tt26*tt20-lam(1,1)*tt&
&150*tt119*tt16*tt20+(mu(1,1)*(2*tt150*tt119*tt16-2*tt655*tt26))/2&
&.0E+0+lam(1,1)*tt150*tt119*tt16)
tt657 = tt127**2
tt658 = detDmH(1,1)*gw(1,1)*((-lam(1,1)*tt133*tt127*tt16*tt20)+mu&
&(1,1)*tt133*tt127*tt16+lam(1,1)*tt133*tt127*tt16)
tt659 = H_invDmH(8,1)*tt126-H_invDmH(8,2)*tt125+tt124*H_invDmH(8,&
&3)
tt660 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt659*tt26*tt20-lam(1,1)*tt&
&127*tt136*tt16*tt20+(mu(1,1)*(2*tt127*tt136*tt16-2*tt659*tt26))/2&
&.0E+0+lam(1,1)*tt127*tt136*tt16)
tt661 = detDmH(1,1)*gw(1,1)*((-lam(1,1)*tt127*tt144*tt16*tt20)+(m&
&u(1,1)*(2*tt127*tt144*tt16+tt651+tt650+tt649))/2.0E+0+lam(1,1)*tt&
&127*tt144*tt16)
tt662 = tt14*(H_invDmH(7,2)*H_invDmH(8,3)-H_invDmH(7,3)*H_invDmH(&
&8,2))-tt13*(H_invDmH(7,1)*H_invDmH(8,3)-H_invDmH(7,3)*H_invDmH(8,&
&1))+(H_invDmH(7,1)*H_invDmH(8,2)-H_invDmH(7,2)*H_invDmH(8,1))*tt1&
&2
tt663 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt662*tt26*tt20-lam(1,1)*tt&
&150*tt127*tt16*tt20+(mu(1,1)*(2*tt150*tt127*tt16-2*tt662*tt26))/2&
&.0E+0+lam(1,1)*tt150*tt127*tt16)
tt664 = tt133**2
tt665 = H_invDmH(8,1)*tt132-H_invDmH(8,2)*tt131+tt130*H_invDmH(8,&
&3)
tt666 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt665*tt26*tt20-lam(1,1)*tt&
&133*tt136*tt16*tt20+(mu(1,1)*(2*tt133*tt136*tt16-2*tt665*tt26))/2&
&.0E+0+lam(1,1)*tt133*tt136*tt16)
tt667 = tt14*(H_invDmH(7,3)*H_invDmH(8,2)-H_invDmH(7,2)*H_invDmH(&
&8,3))-tt13*(H_invDmH(7,3)*H_invDmH(8,1)-H_invDmH(7,1)*H_invDmH(8,&
&3))+(H_invDmH(7,2)*H_invDmH(8,1)-H_invDmH(7,1)*H_invDmH(8,2))*tt1&
&2
tt668 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt667*tt26*tt20-lam(1,1)*tt&
&133*tt144*tt16*tt20+(mu(1,1)*(2*tt133*tt144*tt16-2*tt667*tt26))/2&
&.0E+0+lam(1,1)*tt133*tt144*tt16)
tt669 = detDmH(1,1)*gw(1,1)*((-lam(1,1)*tt133*tt150*tt16*tt20)+(m&
&u(1,1)*(2*tt133*tt150*tt16+tt651+tt650+tt649))/2.0E+0+lam(1,1)*tt&
&133*tt150*tt16)
tt670 = tt136**2
tt671 = 2*H_invDmH(8,1)**2
tt672 = 2*H_invDmH(8,2)**2
tt673 = 2*H_invDmH(8,3)**2
tt674 = H_invDmH(8,1)*tt143-H_invDmH(8,2)*tt142+tt141*H_invDmH(8,&
&3)
tt675 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt674*tt26*tt20-lam(1,1)*tt&
&144*tt136*tt16*tt20+(mu(1,1)*(2*tt144*tt136*tt16-2*tt674*tt26))/2&
&.0E+0+lam(1,1)*tt144*tt136*tt16)
tt676 = H_invDmH(8,1)*tt149-H_invDmH(8,2)*tt148+tt147*H_invDmH(8,&
&3)
tt677 = detDmH(1,1)*gw(1,1)*(lam(1,1)*tt676*tt26*tt20-lam(1,1)*tt&
&150*tt136*tt16*tt20+(mu(1,1)*(2*tt150*tt136*tt16-2*tt676*tt26))/2&
&.0E+0+lam(1,1)*tt150*tt136*tt16)
tt678 = tt144**2
tt679 = detDmH(1,1)*gw(1,1)*((-lam(1,1)*tt150*tt144*tt16*tt20)+mu&
&(1,1)*tt150*tt144*tt16+lam(1,1)*tt150*tt144*tt16)
tt680 = tt150**2
hes(1,1) = detDmH(1,1)*gw(1,1)*((-lam(1,1)*tt11*tt16*tt20)+(mu(1,&
&1)*(2*tt11*tt16+tt19+tt18+tt17))/2.0E+0+lam(1,1)*tt11*tt16)
hes(1,2) = tt27
hes(1,3) = tt33
hes(1,4) = tt38
hes(1,5) = tt44
hes(1,6) = tt50
hes(1,7) = tt55
hes(1,8) = tt61
hes(1,9) = tt67
hes(1,10) = tt72
hes(1,11) = tt78
hes(1,12) = tt84
hes(1,13) = tt89
hes(1,14) = tt95
hes(1,15) = tt101
hes(1,16) = tt106
hes(1,17) = tt112
hes(1,18) = tt118
hes(1,19) = tt123
hes(1,20) = tt129
hes(1,21) = tt135
hes(1,22) = tt140
hes(1,23) = tt146
hes(1,24) = tt152
hes(2,1) = tt27
hes(2,2) = detDmH(1,1)*gw(1,1)*((-lam(1,1)*tt153*tt16*tt20)+(mu(1&
&,1)*(2*tt153*tt16+tt19+tt18+tt17))/2.0E+0+lam(1,1)*tt153*tt16)
hes(2,3) = tt154
hes(2,4) = tt156
hes(2,5) = tt157
hes(2,6) = tt159
hes(2,7) = tt161
hes(2,8) = tt162
hes(2,9) = tt164
hes(2,10) = tt166
hes(2,11) = tt167
hes(2,12) = tt169
hes(2,13) = tt171
hes(2,14) = tt172
hes(2,15) = tt174
hes(2,16) = tt176
hes(2,17) = tt177
hes(2,18) = tt179
hes(2,19) = tt181
hes(2,20) = tt182
hes(2,21) = tt184
hes(2,22) = tt186
hes(2,23) = tt187
hes(2,24) = tt189
hes(3,1) = tt33
hes(3,2) = tt154
hes(3,3) = detDmH(1,1)*gw(1,1)*((-lam(1,1)*tt190*tt16*tt20)+(mu(1&
&,1)*(2*tt190*tt16+tt19+tt18+tt17))/2.0E+0+lam(1,1)*tt190*tt16)
hes(3,4) = tt192
hes(3,5) = tt194
hes(3,6) = tt195
hes(3,7) = tt197
hes(3,8) = tt199
hes(3,9) = tt200
hes(3,10) = tt202
hes(3,11) = tt204
hes(3,12) = tt205
hes(3,13) = tt207
hes(3,14) = tt209
hes(3,15) = tt210
hes(3,16) = tt212
hes(3,17) = tt214
hes(3,18) = tt215
hes(3,19) = tt217
hes(3,20) = tt219
hes(3,21) = tt220
hes(3,22) = tt222
hes(3,23) = tt224
hes(3,24) = tt225
hes(4,1) = tt38
hes(4,2) = tt156
hes(4,3) = tt192
hes(4,4) = detDmH(1,1)*gw(1,1)*((-lam(1,1)*tt226*tt16*tt20)+(mu(1&
&,1)*(2*tt226*tt16+tt229+tt228+tt227))/2.0E+0+lam(1,1)*tt226*tt16)
hes(4,5) = tt231
hes(4,6) = tt233
hes(4,7) = tt237
hes(4,8) = tt239
hes(4,9) = tt241
hes(4,10) = tt245
hes(4,11) = tt247
hes(4,12) = tt249
hes(4,13) = tt253
hes(4,14) = tt255
hes(4,15) = tt257
hes(4,16) = tt261
hes(4,17) = tt263
hes(4,18) = tt265
hes(4,19) = tt269
hes(4,20) = tt271
hes(4,21) = tt273
hes(4,22) = tt277
hes(4,23) = tt279
hes(4,24) = tt281
hes(5,1) = tt44
hes(5,2) = tt157
hes(5,3) = tt194
hes(5,4) = tt231
hes(5,5) = detDmH(1,1)*gw(1,1)*((-lam(1,1)*tt282*tt16*tt20)+(mu(1&
&,1)*(2*tt282*tt16+tt229+tt228+tt227))/2.0E+0+lam(1,1)*tt282*tt16)
hes(5,6) = tt283
hes(5,7) = tt285
hes(5,8) = tt286
hes(5,9) = tt288
hes(5,10) = tt290
hes(5,11) = tt291
hes(5,12) = tt293
hes(5,13) = tt295
hes(5,14) = tt296
hes(5,15) = tt298
hes(5,16) = tt300
hes(5,17) = tt301
hes(5,18) = tt303
hes(5,19) = tt305
hes(5,20) = tt306
hes(5,21) = tt308
hes(5,22) = tt310
hes(5,23) = tt311
hes(5,24) = tt313
hes(6,1) = tt50
hes(6,2) = tt159
hes(6,3) = tt195
hes(6,4) = tt233
hes(6,5) = tt283
hes(6,6) = detDmH(1,1)*gw(1,1)*((-lam(1,1)*tt314*tt16*tt20)+(mu(1&
&,1)*(2*tt314*tt16+tt229+tt228+tt227))/2.0E+0+lam(1,1)*tt314*tt16)
hes(6,7) = tt316
hes(6,8) = tt318
hes(6,9) = tt319
hes(6,10) = tt321
hes(6,11) = tt323
hes(6,12) = tt324
hes(6,13) = tt326
hes(6,14) = tt328
hes(6,15) = tt329
hes(6,16) = tt331
hes(6,17) = tt333
hes(6,18) = tt334
hes(6,19) = tt336
hes(6,20) = tt338
hes(6,21) = tt339
hes(6,22) = tt341
hes(6,23) = tt343
hes(6,24) = tt344
hes(7,1) = tt55
hes(7,2) = tt161
hes(7,3) = tt197
hes(7,4) = tt237
hes(7,5) = tt285
hes(7,6) = tt316
hes(7,7) = detDmH(1,1)*gw(1,1)*((-lam(1,1)*tt345*tt16*tt20)+(mu(1&
&,1)*(2*tt345*tt16+tt348+tt347+tt346))/2.0E+0+lam(1,1)*tt345*tt16)
hes(7,8) = tt350
hes(7,9) = tt352
hes(7,10) = tt356
hes(7,11) = tt358
hes(7,12) = tt360
hes(7,13) = tt364
hes(7,14) = tt366
hes(7,15) = tt368
hes(7,16) = tt372
hes(7,17) = tt374
hes(7,18) = tt376
hes(7,19) = tt380
hes(7,20) = tt382
hes(7,21) = tt384
hes(7,22) = tt388
hes(7,23) = tt390
hes(7,24) = tt392
hes(8,1) = tt61
hes(8,2) = tt162
hes(8,3) = tt199
hes(8,4) = tt239
hes(8,5) = tt286
hes(8,6) = tt318
hes(8,7) = tt350
hes(8,8) = detDmH(1,1)*gw(1,1)*((-lam(1,1)*tt393*tt16*tt20)+(mu(1&
&,1)*(2*tt393*tt16+tt348+tt347+tt346))/2.0E+0+lam(1,1)*tt393*tt16)
hes(8,9) = tt394
hes(8,10) = tt396
hes(8,11) = tt397
hes(8,12) = tt399
hes(8,13) = tt401
hes(8,14) = tt402
hes(8,15) = tt404
hes(8,16) = tt406
hes(8,17) = tt407
hes(8,18) = tt409
hes(8,19) = tt411
hes(8,20) = tt412
hes(8,21) = tt414
hes(8,22) = tt416
hes(8,23) = tt417
hes(8,24) = tt419
hes(9,1) = tt67
hes(9,2) = tt164
hes(9,3) = tt200
hes(9,4) = tt241
hes(9,5) = tt288
hes(9,6) = tt319
hes(9,7) = tt352
hes(9,8) = tt394
hes(9,9) = detDmH(1,1)*gw(1,1)*((-lam(1,1)*tt420*tt16*tt20)+(mu(1&
&,1)*(2*tt420*tt16+tt348+tt347+tt346))/2.0E+0+lam(1,1)*tt420*tt16)
hes(9,10) = tt422
hes(9,11) = tt424
hes(9,12) = tt425
hes(9,13) = tt427
hes(9,14) = tt429
hes(9,15) = tt430
hes(9,16) = tt432
hes(9,17) = tt434
hes(9,18) = tt435
hes(9,19) = tt437
hes(9,20) = tt439
hes(9,21) = tt440
hes(9,22) = tt442
hes(9,23) = tt444
hes(9,24) = tt445
hes(10,1) = tt72
hes(10,2) = tt166
hes(10,3) = tt202
hes(10,4) = tt245
hes(10,5) = tt290
hes(10,6) = tt321
hes(10,7) = tt356
hes(10,8) = tt396
hes(10,9) = tt422
hes(10,10) = detDmH(1,1)*gw(1,1)*((-lam(1,1)*tt446*tt16*tt20)+(mu&
&(1,1)*(2*tt446*tt16+tt449+tt448+tt447))/2.0E+0+lam(1,1)*tt446*tt1&
&6)
hes(10,11) = tt451
hes(10,12) = tt453
hes(10,13) = tt457
hes(10,14) = tt459
hes(10,15) = tt461
hes(10,16) = tt465
hes(10,17) = tt467
hes(10,18) = tt469
hes(10,19) = tt473
hes(10,20) = tt475
hes(10,21) = tt477
hes(10,22) = tt481
hes(10,23) = tt483
hes(10,24) = tt485
hes(11,1) = tt78
hes(11,2) = tt167
hes(11,3) = tt204
hes(11,4) = tt247
hes(11,5) = tt291
hes(11,6) = tt323
hes(11,7) = tt358
hes(11,8) = tt397
hes(11,9) = tt424
hes(11,10) = tt451
hes(11,11) = detDmH(1,1)*gw(1,1)*((-lam(1,1)*tt486*tt16*tt20)+(mu&
&(1,1)*(2*tt486*tt16+tt449+tt448+tt447))/2.0E+0+lam(1,1)*tt486*tt1&
&6)
hes(11,12) = tt487
hes(11,13) = tt489
hes(11,14) = tt490
hes(11,15) = tt492
hes(11,16) = tt494
hes(11,17) = tt495
hes(11,18) = tt497
hes(11,19) = tt499
hes(11,20) = tt500
hes(11,21) = tt502
hes(11,22) = tt504
hes(11,23) = tt505
hes(11,24) = tt507
hes(12,1) = tt84
hes(12,2) = tt169
hes(12,3) = tt205
hes(12,4) = tt249
hes(12,5) = tt293
hes(12,6) = tt324
hes(12,7) = tt360
hes(12,8) = tt399
hes(12,9) = tt425
hes(12,10) = tt453
hes(12,11) = tt487
hes(12,12) = detDmH(1,1)*gw(1,1)*((-lam(1,1)*tt508*tt16*tt20)+(mu&
&(1,1)*(2*tt508*tt16+tt449+tt448+tt447))/2.0E+0+lam(1,1)*tt508*tt1&
&6)
hes(12,13) = tt510
hes(12,14) = tt512
hes(12,15) = tt513
hes(12,16) = tt515
hes(12,17) = tt517
hes(12,18) = tt518
hes(12,19) = tt520
hes(12,20) = tt522
hes(12,21) = tt523
hes(12,22) = tt525
hes(12,23) = tt527
hes(12,24) = tt528
hes(13,1) = tt89
hes(13,2) = tt171
hes(13,3) = tt207
hes(13,4) = tt253
hes(13,5) = tt295
hes(13,6) = tt326
hes(13,7) = tt364
hes(13,8) = tt401
hes(13,9) = tt427
hes(13,10) = tt457
hes(13,11) = tt489
hes(13,12) = tt510
hes(13,13) = detDmH(1,1)*gw(1,1)*((-lam(1,1)*tt529*tt16*tt20)+(mu&
&(1,1)*(2*tt529*tt16+tt532+tt531+tt530))/2.0E+0+lam(1,1)*tt529*tt1&
&6)
hes(13,14) = tt534
hes(13,15) = tt536
hes(13,16) = tt540
hes(13,17) = tt542
hes(13,18) = tt544
hes(13,19) = tt548
hes(13,20) = tt550
hes(13,21) = tt552
hes(13,22) = tt556
hes(13,23) = tt558
hes(13,24) = tt560
hes(14,1) = tt95
hes(14,2) = tt172
hes(14,3) = tt209
hes(14,4) = tt255
hes(14,5) = tt296
hes(14,6) = tt328
hes(14,7) = tt366
hes(14,8) = tt402
hes(14,9) = tt429
hes(14,10) = tt459
hes(14,11) = tt490
hes(14,12) = tt512
hes(14,13) = tt534
hes(14,14) = detDmH(1,1)*gw(1,1)*((-lam(1,1)*tt561*tt16*tt20)+(mu&
&(1,1)*(2*tt561*tt16+tt532+tt531+tt530))/2.0E+0+lam(1,1)*tt561*tt1&
&6)
hes(14,15) = tt562
hes(14,16) = tt564
hes(14,17) = tt565
hes(14,18) = tt567
hes(14,19) = tt569
hes(14,20) = tt570
hes(14,21) = tt572
hes(14,22) = tt574
hes(14,23) = tt575
hes(14,24) = tt577
hes(15,1) = tt101
hes(15,2) = tt174
hes(15,3) = tt210
hes(15,4) = tt257
hes(15,5) = tt298
hes(15,6) = tt329
hes(15,7) = tt368
hes(15,8) = tt404
hes(15,9) = tt430
hes(15,10) = tt461
hes(15,11) = tt492
hes(15,12) = tt513
hes(15,13) = tt536
hes(15,14) = tt562
hes(15,15) = detDmH(1,1)*gw(1,1)*((-lam(1,1)*tt578*tt16*tt20)+(mu&
&(1,1)*(2*tt578*tt16+tt532+tt531+tt530))/2.0E+0+lam(1,1)*tt578*tt1&
&6)
hes(15,16) = tt580
hes(15,17) = tt582
hes(15,18) = tt583
hes(15,19) = tt585
hes(15,20) = tt587
hes(15,21) = tt588
hes(15,22) = tt590
hes(15,23) = tt592
hes(15,24) = tt593
hes(16,1) = tt106
hes(16,2) = tt176
hes(16,3) = tt212
hes(16,4) = tt261
hes(16,5) = tt300
hes(16,6) = tt331
hes(16,7) = tt372
hes(16,8) = tt406
hes(16,9) = tt432
hes(16,10) = tt465
hes(16,11) = tt494
hes(16,12) = tt515
hes(16,13) = tt540
hes(16,14) = tt564
hes(16,15) = tt580
hes(16,16) = detDmH(1,1)*gw(1,1)*((-lam(1,1)*tt594*tt16*tt20)+(mu&
&(1,1)*(2*tt594*tt16+tt597+tt596+tt595))/2.0E+0+lam(1,1)*tt594*tt1&
&6)
hes(16,17) = tt599
hes(16,18) = tt601
hes(16,19) = tt605
hes(16,20) = tt607
hes(16,21) = tt609
hes(16,22) = tt613
hes(16,23) = tt615
hes(16,24) = tt617
hes(17,1) = tt112
hes(17,2) = tt177
hes(17,3) = tt214
hes(17,4) = tt263
hes(17,5) = tt301
hes(17,6) = tt333
hes(17,7) = tt374
hes(17,8) = tt407
hes(17,9) = tt434
hes(17,10) = tt467
hes(17,11) = tt495
hes(17,12) = tt517
hes(17,13) = tt542
hes(17,14) = tt565
hes(17,15) = tt582
hes(17,16) = tt599
hes(17,17) = detDmH(1,1)*gw(1,1)*((-lam(1,1)*tt618*tt16*tt20)+(mu&
&(1,1)*(2*tt618*tt16+tt597+tt596+tt595))/2.0E+0+lam(1,1)*tt618*tt1&
&6)
hes(17,18) = tt619
hes(17,19) = tt621
hes(17,20) = tt622
hes(17,21) = tt624
hes(17,22) = tt626
hes(17,23) = tt627
hes(17,24) = tt629
hes(18,1) = tt118
hes(18,2) = tt179
hes(18,3) = tt215
hes(18,4) = tt265
hes(18,5) = tt303
hes(18,6) = tt334
hes(18,7) = tt376
hes(18,8) = tt409
hes(18,9) = tt435
hes(18,10) = tt469
hes(18,11) = tt497
hes(18,12) = tt518
hes(18,13) = tt544
hes(18,14) = tt567
hes(18,15) = tt583
hes(18,16) = tt601
hes(18,17) = tt619
hes(18,18) = detDmH(1,1)*gw(1,1)*((-lam(1,1)*tt630*tt16*tt20)+(mu&
&(1,1)*(2*tt630*tt16+tt597+tt596+tt595))/2.0E+0+lam(1,1)*tt630*tt1&
&6)
hes(18,19) = tt632
hes(18,20) = tt634
hes(18,21) = tt635
hes(18,22) = tt637
hes(18,23) = tt639
hes(18,24) = tt640
hes(19,1) = tt123
hes(19,2) = tt181
hes(19,3) = tt217
hes(19,4) = tt269
hes(19,5) = tt305
hes(19,6) = tt336
hes(19,7) = tt380
hes(19,8) = tt411
hes(19,9) = tt437
hes(19,10) = tt473
hes(19,11) = tt499
hes(19,12) = tt520
hes(19,13) = tt548
hes(19,14) = tt569
hes(19,15) = tt585
hes(19,16) = tt605
hes(19,17) = tt621
hes(19,18) = tt632
hes(19,19) = detDmH(1,1)*gw(1,1)*((-lam(1,1)*tt641*tt16*tt20)+(mu&
&(1,1)*(2*tt641*tt16+tt644+tt643+tt642))/2.0E+0+lam(1,1)*tt641*tt1&
&6)
hes(19,20) = tt646
hes(19,21) = tt648
hes(19,22) = tt652
hes(19,23) = tt654
hes(19,24) = tt656
hes(20,1) = tt129
hes(20,2) = tt182
hes(20,3) = tt219
hes(20,4) = tt271
hes(20,5) = tt306
hes(20,6) = tt338
hes(20,7) = tt382
hes(20,8) = tt412
hes(20,9) = tt439
hes(20,10) = tt475
hes(20,11) = tt500
hes(20,12) = tt522
hes(20,13) = tt550
hes(20,14) = tt570
hes(20,15) = tt587
hes(20,16) = tt607
hes(20,17) = tt622
hes(20,18) = tt634
hes(20,19) = tt646
hes(20,20) = detDmH(1,1)*gw(1,1)*((-lam(1,1)*tt657*tt16*tt20)+(mu&
&(1,1)*(2*tt657*tt16+tt644+tt643+tt642))/2.0E+0+lam(1,1)*tt657*tt1&
&6)
hes(20,21) = tt658
hes(20,22) = tt660
hes(20,23) = tt661
hes(20,24) = tt663
hes(21,1) = tt135
hes(21,2) = tt184
hes(21,3) = tt220
hes(21,4) = tt273
hes(21,5) = tt308
hes(21,6) = tt339
hes(21,7) = tt384
hes(21,8) = tt414
hes(21,9) = tt440
hes(21,10) = tt477
hes(21,11) = tt502
hes(21,12) = tt523
hes(21,13) = tt552
hes(21,14) = tt572
hes(21,15) = tt588
hes(21,16) = tt609
hes(21,17) = tt624
hes(21,18) = tt635
hes(21,19) = tt648
hes(21,20) = tt658
hes(21,21) = detDmH(1,1)*gw(1,1)*((-lam(1,1)*tt664*tt16*tt20)+(mu&
&(1,1)*(2*tt664*tt16+tt644+tt643+tt642))/2.0E+0+lam(1,1)*tt664*tt1&
&6)
hes(21,22) = tt666
hes(21,23) = tt668
hes(21,24) = tt669
hes(22,1) = tt140
hes(22,2) = tt186
hes(22,3) = tt222
hes(22,4) = tt277
hes(22,5) = tt310
hes(22,6) = tt341
hes(22,7) = tt388
hes(22,8) = tt416
hes(22,9) = tt442
hes(22,10) = tt481
hes(22,11) = tt504
hes(22,12) = tt525
hes(22,13) = tt556
hes(22,14) = tt574
hes(22,15) = tt590
hes(22,16) = tt613
hes(22,17) = tt626
hes(22,18) = tt637
hes(22,19) = tt652
hes(22,20) = tt660
hes(22,21) = tt666
hes(22,22) = detDmH(1,1)*gw(1,1)*((-lam(1,1)*tt670*tt16*tt20)+(mu&
&(1,1)*(2*tt670*tt16+tt673+tt672+tt671))/2.0E+0+lam(1,1)*tt670*tt1&
&6)
hes(22,23) = tt675
hes(22,24) = tt677
hes(23,1) = tt146
hes(23,2) = tt187
hes(23,3) = tt224
hes(23,4) = tt279
hes(23,5) = tt311
hes(23,6) = tt343
hes(23,7) = tt390
hes(23,8) = tt417
hes(23,9) = tt444
hes(23,10) = tt483
hes(23,11) = tt505
hes(23,12) = tt527
hes(23,13) = tt558
hes(23,14) = tt575
hes(23,15) = tt592
hes(23,16) = tt615
hes(23,17) = tt627
hes(23,18) = tt639
hes(23,19) = tt654
hes(23,20) = tt661
hes(23,21) = tt668
hes(23,22) = tt675
hes(23,23) = detDmH(1,1)*gw(1,1)*((-lam(1,1)*tt678*tt16*tt20)+(mu&
&(1,1)*(2*tt678*tt16+tt673+tt672+tt671))/2.0E+0+lam(1,1)*tt678*tt1&
&6)
hes(23,24) = tt679
hes(24,1) = tt152
hes(24,2) = tt189
hes(24,3) = tt225
hes(24,4) = tt281
hes(24,5) = tt313
hes(24,6) = tt344
hes(24,7) = tt392
hes(24,8) = tt419
hes(24,9) = tt445
hes(24,10) = tt485
hes(24,11) = tt507
hes(24,12) = tt528
hes(24,13) = tt560
hes(24,14) = tt577
hes(24,15) = tt593
hes(24,16) = tt617
hes(24,17) = tt629
hes(24,18) = tt640
hes(24,19) = tt656
hes(24,20) = tt663
hes(24,21) = tt669
hes(24,22) = tt677
hes(24,23) = tt679
hes(24,24) = detDmH(1,1)*gw(1,1)*((-lam(1,1)*tt680*tt16*tt20)+(mu&
&(1,1)*(2*tt680*tt16+tt673+tt672+tt671))/2.0E+0+lam(1,1)*tt680*tt1&
&6)
END 
SUBROUTINE vox_sta_neo_at_quadr(val, X, H_invDmH, detDmH, gw, mu, lam) 
IMPLICIT NONE 
REAL(KIND=8) val(1, 1) 
REAL(KIND=8) X(3, 8) 
REAL(KIND=8) H_invDmH(8, 3) 
REAL(KIND=8) detDmH(1, 1) 
REAL(KIND=8) gw(1, 1) 
REAL(KIND=8) mu(1, 1) 
REAL(KIND=8) lam(1, 1) 
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
tt1 = X(3,8)*H_invDmH(8,1)+X(3,7)*H_invDmH(7,1)+X(3,6)*H_invDmH(6&
&,1)+X(3,5)*H_invDmH(5,1)+X(3,4)*H_invDmH(4,1)+H_invDmH(3,1)*X(3,3&
&)+H_invDmH(2,1)*X(3,2)+H_invDmH(1,1)*X(3,1)
tt2 = X(2,8)*H_invDmH(8,2)+X(2,7)*H_invDmH(7,2)+X(2,6)*H_invDmH(6&
&,2)+X(2,5)*H_invDmH(5,2)+X(2,4)*H_invDmH(4,2)+X(2,3)*H_invDmH(3,2&
&)+H_invDmH(2,2)*X(2,2)+H_invDmH(1,2)*X(2,1)
tt3 = X(2,8)*H_invDmH(8,1)+X(2,7)*H_invDmH(7,1)+X(2,6)*H_invDmH(6&
&,1)+X(2,5)*H_invDmH(5,1)+X(2,4)*H_invDmH(4,1)+X(2,3)*H_invDmH(3,1&
&)+H_invDmH(2,1)*X(2,2)+H_invDmH(1,1)*X(2,1)
tt4 = X(3,8)*H_invDmH(8,2)+X(3,7)*H_invDmH(7,2)+X(3,6)*H_invDmH(6&
&,2)+X(3,5)*H_invDmH(5,2)+X(3,4)*H_invDmH(4,2)+H_invDmH(3,2)*X(3,3&
&)+H_invDmH(2,2)*X(3,2)+H_invDmH(1,2)*X(3,1)
tt5 = X(1,8)*H_invDmH(8,3)+X(1,7)*H_invDmH(7,3)+X(1,6)*H_invDmH(6&
&,3)+X(1,5)*H_invDmH(5,3)+X(1,4)*H_invDmH(4,3)+X(1,3)*H_invDmH(3,3&
&)+X(1,2)*H_invDmH(2,3)+X(1,1)*H_invDmH(1,3)
tt6 = X(1,8)*H_invDmH(8,2)+X(1,7)*H_invDmH(7,2)+X(1,6)*H_invDmH(6&
&,2)+X(1,5)*H_invDmH(5,2)+X(1,4)*H_invDmH(4,2)+X(1,3)*H_invDmH(3,2&
&)+X(1,2)*H_invDmH(2,2)+X(1,1)*H_invDmH(1,2)
tt7 = X(2,8)*H_invDmH(8,3)+X(2,7)*H_invDmH(7,3)+X(2,6)*H_invDmH(6&
&,3)+X(2,5)*H_invDmH(5,3)+X(2,4)*H_invDmH(4,3)+X(2,3)*H_invDmH(3,3&
&)+X(2,2)*H_invDmH(2,3)+H_invDmH(1,3)*X(2,1)
tt8 = X(3,8)*H_invDmH(8,3)+X(3,7)*H_invDmH(7,3)+X(3,6)*H_invDmH(6&
&,3)+X(3,5)*H_invDmH(5,3)+X(3,4)*H_invDmH(4,3)+H_invDmH(3,3)*X(3,3&
&)+H_invDmH(2,3)*X(3,2)+H_invDmH(1,3)*X(3,1)
tt9 = X(1,8)*H_invDmH(8,1)+X(1,7)*H_invDmH(7,1)+X(1,6)*H_invDmH(6&
&,1)+X(1,5)*H_invDmH(5,1)+X(1,4)*H_invDmH(4,1)+X(1,3)*H_invDmH(3,1&
&)+X(1,2)*H_invDmH(2,1)+H_invDmH(1,1)*X(1,1)
tt10 = tt9*(tt2*tt8-tt4*tt7)-tt6*(tt3*tt8-tt1*tt7)+(tt3*tt4-tt1*t&
&t2)*tt5-1
val(1,1) = detDmH(1,1)*gw(1,1)*(((mu(1,1)+lam(1,1))*tt10**2)/2.0E&
&+0+(mu(1,1)*(tt8**2+tt7**2+tt5**2+tt4**2+tt2**2+tt6**2+tt1**2+tt3&
&**2+tt9**2-3))/2.0E+0-mu(1,1)*tt10)
END 
SUBROUTINE vox_sta_neo_at_quadr_jac(jac, X, H_invDmH, detDmH, gw, mu, lam) 
IMPLICIT NONE 
REAL(KIND=8) jac(1, 24) 
REAL(KIND=8) X(3, 8) 
REAL(KIND=8) H_invDmH(8, 3) 
REAL(KIND=8) detDmH(1, 1) 
REAL(KIND=8) gw(1, 1) 
REAL(KIND=8) mu(1, 1) 
REAL(KIND=8) lam(1, 1) 
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
tt1 = X(1,8)*H_invDmH(8,1)+X(1,7)*H_invDmH(7,1)+X(1,6)*H_invDmH(6&
&,1)+X(1,5)*H_invDmH(5,1)+X(1,4)*H_invDmH(4,1)+X(1,3)*H_invDmH(3,1&
&)+X(1,2)*H_invDmH(2,1)+H_invDmH(1,1)*X(1,1)
tt2 = X(1,8)*H_invDmH(8,2)+X(1,7)*H_invDmH(7,2)+X(1,6)*H_invDmH(6&
&,2)+X(1,5)*H_invDmH(5,2)+X(1,4)*H_invDmH(4,2)+X(1,3)*H_invDmH(3,2&
&)+X(1,2)*H_invDmH(2,2)+X(1,1)*H_invDmH(1,2)
tt3 = X(1,8)*H_invDmH(8,3)+X(1,7)*H_invDmH(7,3)+X(1,6)*H_invDmH(6&
&,3)+X(1,5)*H_invDmH(5,3)+X(1,4)*H_invDmH(4,3)+X(1,3)*H_invDmH(3,3&
&)+X(1,2)*H_invDmH(2,3)+X(1,1)*H_invDmH(1,3)
tt4 = X(3,8)*H_invDmH(8,1)+X(3,7)*H_invDmH(7,1)+X(3,6)*H_invDmH(6&
&,1)+X(3,5)*H_invDmH(5,1)+X(3,4)*H_invDmH(4,1)+H_invDmH(3,1)*X(3,3&
&)+H_invDmH(2,1)*X(3,2)+H_invDmH(1,1)*X(3,1)
tt5 = X(2,8)*H_invDmH(8,2)+X(2,7)*H_invDmH(7,2)+X(2,6)*H_invDmH(6&
&,2)+X(2,5)*H_invDmH(5,2)+X(2,4)*H_invDmH(4,2)+X(2,3)*H_invDmH(3,2&
&)+H_invDmH(2,2)*X(2,2)+H_invDmH(1,2)*X(2,1)
tt6 = X(2,8)*H_invDmH(8,1)+X(2,7)*H_invDmH(7,1)+X(2,6)*H_invDmH(6&
&,1)+X(2,5)*H_invDmH(5,1)+X(2,4)*H_invDmH(4,1)+X(2,3)*H_invDmH(3,1&
&)+H_invDmH(2,1)*X(2,2)+H_invDmH(1,1)*X(2,1)
tt7 = X(3,8)*H_invDmH(8,2)+X(3,7)*H_invDmH(7,2)+X(3,6)*H_invDmH(6&
&,2)+X(3,5)*H_invDmH(5,2)+X(3,4)*H_invDmH(4,2)+H_invDmH(3,2)*X(3,3&
&)+H_invDmH(2,2)*X(3,2)+H_invDmH(1,2)*X(3,1)
tt8 = tt6*tt7-tt4*tt5
tt9 = X(2,8)*H_invDmH(8,3)+X(2,7)*H_invDmH(7,3)+X(2,6)*H_invDmH(6&
&,3)+X(2,5)*H_invDmH(5,3)+X(2,4)*H_invDmH(4,3)+X(2,3)*H_invDmH(3,3&
&)+X(2,2)*H_invDmH(2,3)+H_invDmH(1,3)*X(2,1)
tt10 = X(3,8)*H_invDmH(8,3)+X(3,7)*H_invDmH(7,3)+X(3,6)*H_invDmH(&
&6,3)+X(3,5)*H_invDmH(5,3)+X(3,4)*H_invDmH(4,3)+H_invDmH(3,3)*X(3,&
&3)+H_invDmH(2,3)*X(3,2)+H_invDmH(1,3)*X(3,1)
tt11 = tt6*tt10-tt4*tt9
tt12 = tt5*tt10-tt7*tt9
tt13 = H_invDmH(1,1)*tt12-H_invDmH(1,2)*tt11+H_invDmH(1,3)*tt8
tt14 = mu(1,1)+lam(1,1)
tt15 = tt1*tt12-tt2*tt11+tt8*tt3-1
tt16 = tt1*(H_invDmH(1,2)*tt10-H_invDmH(1,3)*tt7)-tt2*(H_invDmH(1&
&,1)*tt10-H_invDmH(1,3)*tt4)+(H_invDmH(1,1)*tt7-H_invDmH(1,2)*tt4)&
&*tt3
tt17 = tt1*(H_invDmH(1,3)*tt5-H_invDmH(1,2)*tt9)-tt2*(H_invDmH(1,&
&3)*tt6-H_invDmH(1,1)*tt9)+(H_invDmH(1,2)*tt6-H_invDmH(1,1)*tt5)*t&
&t3
tt18 = H_invDmH(2,1)*tt12-H_invDmH(2,2)*tt11+H_invDmH(2,3)*tt8
tt19 = tt1*(H_invDmH(2,2)*tt10-H_invDmH(2,3)*tt7)-tt2*(H_invDmH(2&
&,1)*tt10-H_invDmH(2,3)*tt4)+(H_invDmH(2,1)*tt7-H_invDmH(2,2)*tt4)&
&*tt3
tt20 = tt1*(H_invDmH(2,3)*tt5-H_invDmH(2,2)*tt9)-tt2*(H_invDmH(2,&
&3)*tt6-H_invDmH(2,1)*tt9)+(H_invDmH(2,2)*tt6-H_invDmH(2,1)*tt5)*t&
&t3
tt21 = H_invDmH(3,1)*tt12-H_invDmH(3,2)*tt11+H_invDmH(3,3)*tt8
tt22 = tt1*(H_invDmH(3,2)*tt10-H_invDmH(3,3)*tt7)-tt2*(H_invDmH(3&
&,1)*tt10-H_invDmH(3,3)*tt4)+(H_invDmH(3,1)*tt7-H_invDmH(3,2)*tt4)&
&*tt3
tt23 = tt1*(H_invDmH(3,3)*tt5-H_invDmH(3,2)*tt9)-tt2*(H_invDmH(3,&
&3)*tt6-H_invDmH(3,1)*tt9)+(H_invDmH(3,2)*tt6-H_invDmH(3,1)*tt5)*t&
&t3
tt24 = H_invDmH(4,1)*tt12-H_invDmH(4,2)*tt11+H_invDmH(4,3)*tt8
tt25 = tt1*(H_invDmH(4,2)*tt10-H_invDmH(4,3)*tt7)-tt2*(H_invDmH(4&
&,1)*tt10-H_invDmH(4,3)*tt4)+(H_invDmH(4,1)*tt7-H_invDmH(4,2)*tt4)&
&*tt3
tt26 = tt1*(H_invDmH(4,3)*tt5-H_invDmH(4,2)*tt9)-tt2*(H_invDmH(4,&
&3)*tt6-H_invDmH(4,1)*tt9)+(H_invDmH(4,2)*tt6-H_invDmH(4,1)*tt5)*t&
&t3
tt27 = H_invDmH(5,1)*tt12-H_invDmH(5,2)*tt11+H_invDmH(5,3)*tt8
tt28 = tt1*(H_invDmH(5,2)*tt10-H_invDmH(5,3)*tt7)-tt2*(H_invDmH(5&
&,1)*tt10-H_invDmH(5,3)*tt4)+(H_invDmH(5,1)*tt7-H_invDmH(5,2)*tt4)&
&*tt3
tt29 = tt1*(H_invDmH(5,3)*tt5-H_invDmH(5,2)*tt9)-tt2*(H_invDmH(5,&
&3)*tt6-H_invDmH(5,1)*tt9)+(H_invDmH(5,2)*tt6-H_invDmH(5,1)*tt5)*t&
&t3
tt30 = H_invDmH(6,1)*tt12-H_invDmH(6,2)*tt11+H_invDmH(6,3)*tt8
tt31 = tt1*(H_invDmH(6,2)*tt10-H_invDmH(6,3)*tt7)-tt2*(H_invDmH(6&
&,1)*tt10-H_invDmH(6,3)*tt4)+(H_invDmH(6,1)*tt7-H_invDmH(6,2)*tt4)&
&*tt3
tt32 = tt1*(H_invDmH(6,3)*tt5-H_invDmH(6,2)*tt9)-tt2*(H_invDmH(6,&
&3)*tt6-H_invDmH(6,1)*tt9)+(H_invDmH(6,2)*tt6-H_invDmH(6,1)*tt5)*t&
&t3
tt33 = H_invDmH(7,1)*tt12-H_invDmH(7,2)*tt11+H_invDmH(7,3)*tt8
tt34 = tt1*(H_invDmH(7,2)*tt10-H_invDmH(7,3)*tt7)-tt2*(H_invDmH(7&
&,1)*tt10-H_invDmH(7,3)*tt4)+(H_invDmH(7,1)*tt7-H_invDmH(7,2)*tt4)&
&*tt3
tt35 = tt1*(H_invDmH(7,3)*tt5-H_invDmH(7,2)*tt9)-tt2*(H_invDmH(7,&
&3)*tt6-H_invDmH(7,1)*tt9)+(H_invDmH(7,2)*tt6-H_invDmH(7,1)*tt5)*t&
&t3
tt36 = H_invDmH(8,1)*tt12-H_invDmH(8,2)*tt11+tt8*H_invDmH(8,3)
tt37 = tt1*(H_invDmH(8,2)*tt10-tt7*H_invDmH(8,3))-tt2*(H_invDmH(8&
&,1)*tt10-tt4*H_invDmH(8,3))+(H_invDmH(8,1)*tt7-tt4*H_invDmH(8,2))&
&*tt3
tt38 = tt1*(tt5*H_invDmH(8,3)-H_invDmH(8,2)*tt9)-tt2*(tt6*H_invDm&
&H(8,3)-H_invDmH(8,1)*tt9)+(tt6*H_invDmH(8,2)-H_invDmH(8,1)*tt5)*t&
&t3
jac(1,1) = detDmH(1,1)*gw(1,1)*(tt14*tt13*tt15-mu(1,1)*tt13+(mu(1&
&,1)*(2*H_invDmH(1,3)*tt3+2*H_invDmH(1,2)*tt2+2*H_invDmH(1,1)*tt1)&
&)/2.0E+0)
jac(1,2) = detDmH(1,1)*gw(1,1)*(tt14*tt16*tt15-mu(1,1)*tt16+(mu(1&
&,1)*(2*H_invDmH(1,3)*tt9+2*H_invDmH(1,2)*tt5+2*H_invDmH(1,1)*tt6)&
&)/2.0E+0)
jac(1,3) = detDmH(1,1)*gw(1,1)*(tt14*tt17*tt15-mu(1,1)*tt17+(mu(1&
&,1)*(2*H_invDmH(1,3)*tt10+2*H_invDmH(1,2)*tt7+2*H_invDmH(1,1)*tt4&
&))/2.0E+0)
jac(1,4) = detDmH(1,1)*gw(1,1)*(tt14*tt18*tt15-mu(1,1)*tt18+(mu(1&
&,1)*(2*H_invDmH(2,3)*tt3+2*H_invDmH(2,2)*tt2+2*H_invDmH(2,1)*tt1)&
&)/2.0E+0)
jac(1,5) = detDmH(1,1)*gw(1,1)*(tt14*tt19*tt15-mu(1,1)*tt19+(mu(1&
&,1)*(2*H_invDmH(2,3)*tt9+2*H_invDmH(2,2)*tt5+2*H_invDmH(2,1)*tt6)&
&)/2.0E+0)
jac(1,6) = detDmH(1,1)*gw(1,1)*(tt14*tt20*tt15-mu(1,1)*tt20+(mu(1&
&,1)*(2*H_invDmH(2,3)*tt10+2*H_invDmH(2,2)*tt7+2*H_invDmH(2,1)*tt4&
&))/2.0E+0)
jac(1,7) = detDmH(1,1)*gw(1,1)*(tt14*tt21*tt15-mu(1,1)*tt21+(mu(1&
&,1)*(2*H_invDmH(3,3)*tt3+2*H_invDmH(3,2)*tt2+2*H_invDmH(3,1)*tt1)&
&)/2.0E+0)
jac(1,8) = detDmH(1,1)*gw(1,1)*(tt14*tt22*tt15-mu(1,1)*tt22+(mu(1&
&,1)*(2*H_invDmH(3,3)*tt9+2*H_invDmH(3,2)*tt5+2*H_invDmH(3,1)*tt6)&
&)/2.0E+0)
jac(1,9) = detDmH(1,1)*gw(1,1)*(tt14*tt23*tt15-mu(1,1)*tt23+(mu(1&
&,1)*(2*H_invDmH(3,3)*tt10+2*H_invDmH(3,2)*tt7+2*H_invDmH(3,1)*tt4&
&))/2.0E+0)
jac(1,10) = detDmH(1,1)*gw(1,1)*(tt14*tt24*tt15-mu(1,1)*tt24+(mu(&
&1,1)*(2*H_invDmH(4,3)*tt3+2*H_invDmH(4,2)*tt2+2*H_invDmH(4,1)*tt1&
&))/2.0E+0)
jac(1,11) = detDmH(1,1)*gw(1,1)*(tt14*tt25*tt15-mu(1,1)*tt25+(mu(&
&1,1)*(2*H_invDmH(4,3)*tt9+2*H_invDmH(4,2)*tt5+2*H_invDmH(4,1)*tt6&
&))/2.0E+0)
jac(1,12) = detDmH(1,1)*gw(1,1)*(tt14*tt26*tt15-mu(1,1)*tt26+(mu(&
&1,1)*(2*H_invDmH(4,3)*tt10+2*H_invDmH(4,2)*tt7+2*H_invDmH(4,1)*tt&
&4))/2.0E+0)
jac(1,13) = detDmH(1,1)*gw(1,1)*(tt14*tt27*tt15-mu(1,1)*tt27+(mu(&
&1,1)*(2*H_invDmH(5,3)*tt3+2*H_invDmH(5,2)*tt2+2*H_invDmH(5,1)*tt1&
&))/2.0E+0)
jac(1,14) = detDmH(1,1)*gw(1,1)*(tt14*tt28*tt15-mu(1,1)*tt28+(mu(&
&1,1)*(2*H_invDmH(5,3)*tt9+2*H_invDmH(5,2)*tt5+2*H_invDmH(5,1)*tt6&
&))/2.0E+0)
jac(1,15) = detDmH(1,1)*gw(1,1)*(tt14*tt29*tt15-mu(1,1)*tt29+(mu(&
&1,1)*(2*H_invDmH(5,3)*tt10+2*H_invDmH(5,2)*tt7+2*H_invDmH(5,1)*tt&
&4))/2.0E+0)
jac(1,16) = detDmH(1,1)*gw(1,1)*(tt14*tt30*tt15-mu(1,1)*tt30+(mu(&
&1,1)*(2*H_invDmH(6,3)*tt3+2*H_invDmH(6,2)*tt2+2*H_invDmH(6,1)*tt1&
&))/2.0E+0)
jac(1,17) = detDmH(1,1)*gw(1,1)*(tt14*tt31*tt15-mu(1,1)*tt31+(mu(&
&1,1)*(2*H_invDmH(6,3)*tt9+2*H_invDmH(6,2)*tt5+2*H_invDmH(6,1)*tt6&
&))/2.0E+0)
jac(1,18) = detDmH(1,1)*gw(1,1)*(tt14*tt32*tt15-mu(1,1)*tt32+(mu(&
&1,1)*(2*H_invDmH(6,3)*tt10+2*H_invDmH(6,2)*tt7+2*H_invDmH(6,1)*tt&
&4))/2.0E+0)
jac(1,19) = detDmH(1,1)*gw(1,1)*(tt14*tt33*tt15-mu(1,1)*tt33+(mu(&
&1,1)*(2*H_invDmH(7,3)*tt3+2*H_invDmH(7,2)*tt2+2*H_invDmH(7,1)*tt1&
&))/2.0E+0)
jac(1,20) = detDmH(1,1)*gw(1,1)*(tt14*tt34*tt15-mu(1,1)*tt34+(mu(&
&1,1)*(2*H_invDmH(7,3)*tt9+2*H_invDmH(7,2)*tt5+2*H_invDmH(7,1)*tt6&
&))/2.0E+0)
jac(1,21) = detDmH(1,1)*gw(1,1)*(tt14*tt35*tt15-mu(1,1)*tt35+(mu(&
&1,1)*(2*H_invDmH(7,3)*tt10+2*H_invDmH(7,2)*tt7+2*H_invDmH(7,1)*tt&
&4))/2.0E+0)
jac(1,22) = detDmH(1,1)*gw(1,1)*(tt14*tt36*tt15-mu(1,1)*tt36+(mu(&
&1,1)*(2*H_invDmH(8,3)*tt3+2*H_invDmH(8,2)*tt2+2*H_invDmH(8,1)*tt1&
&))/2.0E+0)
jac(1,23) = detDmH(1,1)*gw(1,1)*(tt14*tt37*tt15-mu(1,1)*tt37+(mu(&
&1,1)*(2*H_invDmH(8,3)*tt9+2*H_invDmH(8,2)*tt5+2*H_invDmH(8,1)*tt6&
&))/2.0E+0)
jac(1,24) = detDmH(1,1)*gw(1,1)*(tt14*tt38*tt15-mu(1,1)*tt38+(mu(&
&1,1)*(2*H_invDmH(8,3)*tt10+2*H_invDmH(8,2)*tt7+2*H_invDmH(8,1)*tt&
&4))/2.0E+0)
END 
SUBROUTINE vox_sta_neo_at_quadr_hes(hes, X, H_invDmH, detDmH, gw, mu, lam) 
IMPLICIT NONE 
REAL(KIND=8) hes(24, 24) 
REAL(KIND=8) X(3, 8) 
REAL(KIND=8) H_invDmH(8, 3) 
REAL(KIND=8) detDmH(1, 1) 
REAL(KIND=8) gw(1, 1) 
REAL(KIND=8) mu(1, 1) 
REAL(KIND=8) lam(1, 1) 
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
REAL(KIND=8)  tt73 
REAL(KIND=8)  tt74 
REAL(KIND=8)  tt75 
REAL(KIND=8)  tt76 
REAL(KIND=8)  tt77 
REAL(KIND=8)  tt78 
REAL(KIND=8)  tt79 
REAL(KIND=8)  tt80 
REAL(KIND=8)  tt81 
REAL(KIND=8)  tt82 
REAL(KIND=8)  tt83 
REAL(KIND=8)  tt84 
REAL(KIND=8)  tt85 
REAL(KIND=8)  tt86 
REAL(KIND=8)  tt87 
REAL(KIND=8)  tt88 
REAL(KIND=8)  tt89 
REAL(KIND=8)  tt90 
REAL(KIND=8)  tt91 
REAL(KIND=8)  tt92 
REAL(KIND=8)  tt93 
REAL(KIND=8)  tt94 
REAL(KIND=8)  tt95 
REAL(KIND=8)  tt96 
REAL(KIND=8)  tt97 
REAL(KIND=8)  tt98 
REAL(KIND=8)  tt99 
REAL(KIND=8)  tt100 
REAL(KIND=8)  tt101 
REAL(KIND=8)  tt102 
REAL(KIND=8)  tt103 
REAL(KIND=8)  tt104 
REAL(KIND=8)  tt105 
REAL(KIND=8)  tt106 
REAL(KIND=8)  tt107 
REAL(KIND=8)  tt108 
REAL(KIND=8)  tt109 
REAL(KIND=8)  tt110 
REAL(KIND=8)  tt111 
REAL(KIND=8)  tt112 
REAL(KIND=8)  tt113 
REAL(KIND=8)  tt114 
REAL(KIND=8)  tt115 
REAL(KIND=8)  tt116 
REAL(KIND=8)  tt117 
REAL(KIND=8)  tt118 
REAL(KIND=8)  tt119 
REAL(KIND=8)  tt120 
REAL(KIND=8)  tt121 
REAL(KIND=8)  tt122 
REAL(KIND=8)  tt123 
REAL(KIND=8)  tt124 
REAL(KIND=8)  tt125 
REAL(KIND=8)  tt126 
REAL(KIND=8)  tt127 
REAL(KIND=8)  tt128 
REAL(KIND=8)  tt129 
REAL(KIND=8)  tt130 
REAL(KIND=8)  tt131 
REAL(KIND=8)  tt132 
REAL(KIND=8)  tt133 
REAL(KIND=8)  tt134 
REAL(KIND=8)  tt135 
REAL(KIND=8)  tt136 
REAL(KIND=8)  tt137 
REAL(KIND=8)  tt138 
REAL(KIND=8)  tt139 
REAL(KIND=8)  tt140 
REAL(KIND=8)  tt141 
REAL(KIND=8)  tt142 
REAL(KIND=8)  tt143 
REAL(KIND=8)  tt144 
REAL(KIND=8)  tt145 
REAL(KIND=8)  tt146 
REAL(KIND=8)  tt147 
REAL(KIND=8)  tt148 
REAL(KIND=8)  tt149 
REAL(KIND=8)  tt150 
REAL(KIND=8)  tt151 
REAL(KIND=8)  tt152 
REAL(KIND=8)  tt153 
REAL(KIND=8)  tt154 
REAL(KIND=8)  tt155 
REAL(KIND=8)  tt156 
REAL(KIND=8)  tt157 
REAL(KIND=8)  tt158 
REAL(KIND=8)  tt159 
REAL(KIND=8)  tt160 
REAL(KIND=8)  tt161 
REAL(KIND=8)  tt162 
REAL(KIND=8)  tt163 
REAL(KIND=8)  tt164 
REAL(KIND=8)  tt165 
REAL(KIND=8)  tt166 
REAL(KIND=8)  tt167 
REAL(KIND=8)  tt168 
REAL(KIND=8)  tt169 
REAL(KIND=8)  tt170 
REAL(KIND=8)  tt171 
REAL(KIND=8)  tt172 
REAL(KIND=8)  tt173 
REAL(KIND=8)  tt174 
REAL(KIND=8)  tt175 
REAL(KIND=8)  tt176 
REAL(KIND=8)  tt177 
REAL(KIND=8)  tt178 
REAL(KIND=8)  tt179 
REAL(KIND=8)  tt180 
REAL(KIND=8)  tt181 
REAL(KIND=8)  tt182 
REAL(KIND=8)  tt183 
REAL(KIND=8)  tt184 
REAL(KIND=8)  tt185 
REAL(KIND=8)  tt186 
REAL(KIND=8)  tt187 
REAL(KIND=8)  tt188 
REAL(KIND=8)  tt189 
REAL(KIND=8)  tt190 
REAL(KIND=8)  tt191 
REAL(KIND=8)  tt192 
REAL(KIND=8)  tt193 
REAL(KIND=8)  tt194 
REAL(KIND=8)  tt195 
REAL(KIND=8)  tt196 
REAL(KIND=8)  tt197 
REAL(KIND=8)  tt198 
REAL(KIND=8)  tt199 
REAL(KIND=8)  tt200 
REAL(KIND=8)  tt201 
REAL(KIND=8)  tt202 
REAL(KIND=8)  tt203 
REAL(KIND=8)  tt204 
REAL(KIND=8)  tt205 
REAL(KIND=8)  tt206 
REAL(KIND=8)  tt207 
REAL(KIND=8)  tt208 
REAL(KIND=8)  tt209 
REAL(KIND=8)  tt210 
REAL(KIND=8)  tt211 
REAL(KIND=8)  tt212 
REAL(KIND=8)  tt213 
REAL(KIND=8)  tt214 
REAL(KIND=8)  tt215 
REAL(KIND=8)  tt216 
REAL(KIND=8)  tt217 
REAL(KIND=8)  tt218 
REAL(KIND=8)  tt219 
REAL(KIND=8)  tt220 
REAL(KIND=8)  tt221 
REAL(KIND=8)  tt222 
REAL(KIND=8)  tt223 
REAL(KIND=8)  tt224 
REAL(KIND=8)  tt225 
REAL(KIND=8)  tt226 
REAL(KIND=8)  tt227 
REAL(KIND=8)  tt228 
REAL(KIND=8)  tt229 
REAL(KIND=8)  tt230 
REAL(KIND=8)  tt231 
REAL(KIND=8)  tt232 
REAL(KIND=8)  tt233 
REAL(KIND=8)  tt234 
REAL(KIND=8)  tt235 
REAL(KIND=8)  tt236 
REAL(KIND=8)  tt237 
REAL(KIND=8)  tt238 
REAL(KIND=8)  tt239 
REAL(KIND=8)  tt240 
REAL(KIND=8)  tt241 
REAL(KIND=8)  tt242 
REAL(KIND=8)  tt243 
REAL(KIND=8)  tt244 
REAL(KIND=8)  tt245 
REAL(KIND=8)  tt246 
REAL(KIND=8)  tt247 
REAL(KIND=8)  tt248 
REAL(KIND=8)  tt249 
REAL(KIND=8)  tt250 
REAL(KIND=8)  tt251 
REAL(KIND=8)  tt252 
REAL(KIND=8)  tt253 
REAL(KIND=8)  tt254 
REAL(KIND=8)  tt255 
REAL(KIND=8)  tt256 
REAL(KIND=8)  tt257 
REAL(KIND=8)  tt258 
REAL(KIND=8)  tt259 
REAL(KIND=8)  tt260 
REAL(KIND=8)  tt261 
REAL(KIND=8)  tt262 
REAL(KIND=8)  tt263 
REAL(KIND=8)  tt264 
REAL(KIND=8)  tt265 
REAL(KIND=8)  tt266 
REAL(KIND=8)  tt267 
REAL(KIND=8)  tt268 
REAL(KIND=8)  tt269 
REAL(KIND=8)  tt270 
REAL(KIND=8)  tt271 
REAL(KIND=8)  tt272 
REAL(KIND=8)  tt273 
REAL(KIND=8)  tt274 
REAL(KIND=8)  tt275 
REAL(KIND=8)  tt276 
REAL(KIND=8)  tt277 
REAL(KIND=8)  tt278 
REAL(KIND=8)  tt279 
REAL(KIND=8)  tt280 
REAL(KIND=8)  tt281 
REAL(KIND=8)  tt282 
REAL(KIND=8)  tt283 
REAL(KIND=8)  tt284 
REAL(KIND=8)  tt285 
REAL(KIND=8)  tt286 
REAL(KIND=8)  tt287 
REAL(KIND=8)  tt288 
REAL(KIND=8)  tt289 
REAL(KIND=8)  tt290 
REAL(KIND=8)  tt291 
REAL(KIND=8)  tt292 
REAL(KIND=8)  tt293 
REAL(KIND=8)  tt294 
REAL(KIND=8)  tt295 
REAL(KIND=8)  tt296 
REAL(KIND=8)  tt297 
REAL(KIND=8)  tt298 
REAL(KIND=8)  tt299 
REAL(KIND=8)  tt300 
REAL(KIND=8)  tt301 
REAL(KIND=8)  tt302 
REAL(KIND=8)  tt303 
REAL(KIND=8)  tt304 
REAL(KIND=8)  tt305 
REAL(KIND=8)  tt306 
REAL(KIND=8)  tt307 
REAL(KIND=8)  tt308 
REAL(KIND=8)  tt309 
REAL(KIND=8)  tt310 
REAL(KIND=8)  tt311 
REAL(KIND=8)  tt312 
REAL(KIND=8)  tt313 
REAL(KIND=8)  tt314 
REAL(KIND=8)  tt315 
REAL(KIND=8)  tt316 
REAL(KIND=8)  tt317 
REAL(KIND=8)  tt318 
REAL(KIND=8)  tt319 
REAL(KIND=8)  tt320 
REAL(KIND=8)  tt321 
REAL(KIND=8)  tt322 
REAL(KIND=8)  tt323 
REAL(KIND=8)  tt324 
REAL(KIND=8)  tt325 
REAL(KIND=8)  tt326 
REAL(KIND=8)  tt327 
REAL(KIND=8)  tt328 
REAL(KIND=8)  tt329 
REAL(KIND=8)  tt330 
REAL(KIND=8)  tt331 
REAL(KIND=8)  tt332 
REAL(KIND=8)  tt333 
REAL(KIND=8)  tt334 
REAL(KIND=8)  tt335 
REAL(KIND=8)  tt336 
REAL(KIND=8)  tt337 
REAL(KIND=8)  tt338 
REAL(KIND=8)  tt339 
REAL(KIND=8)  tt340 
REAL(KIND=8)  tt341 
REAL(KIND=8)  tt342 
REAL(KIND=8)  tt343 
REAL(KIND=8)  tt344 
REAL(KIND=8)  tt345 
REAL(KIND=8)  tt346 
REAL(KIND=8)  tt347 
REAL(KIND=8)  tt348 
REAL(KIND=8)  tt349 
REAL(KIND=8)  tt350 
REAL(KIND=8)  tt351 
REAL(KIND=8)  tt352 
REAL(KIND=8)  tt353 
REAL(KIND=8)  tt354 
REAL(KIND=8)  tt355 
REAL(KIND=8)  tt356 
REAL(KIND=8)  tt357 
REAL(KIND=8)  tt358 
REAL(KIND=8)  tt359 
REAL(KIND=8)  tt360 
REAL(KIND=8)  tt361 
REAL(KIND=8)  tt362 
REAL(KIND=8)  tt363 
REAL(KIND=8)  tt364 
REAL(KIND=8)  tt365 
REAL(KIND=8)  tt366 
REAL(KIND=8)  tt367 
REAL(KIND=8)  tt368 
REAL(KIND=8)  tt369 
REAL(KIND=8)  tt370 
REAL(KIND=8)  tt371 
REAL(KIND=8)  tt372 
REAL(KIND=8)  tt373 
REAL(KIND=8)  tt374 
REAL(KIND=8)  tt375 
REAL(KIND=8)  tt376 
REAL(KIND=8)  tt377 
REAL(KIND=8)  tt378 
REAL(KIND=8)  tt379 
REAL(KIND=8)  tt380 
REAL(KIND=8)  tt381 
REAL(KIND=8)  tt382 
REAL(KIND=8)  tt383 
REAL(KIND=8)  tt384 
REAL(KIND=8)  tt385 
REAL(KIND=8)  tt386 
REAL(KIND=8)  tt387 
REAL(KIND=8)  tt388 
REAL(KIND=8)  tt389 
REAL(KIND=8)  tt390 
REAL(KIND=8)  tt391 
REAL(KIND=8)  tt392 
REAL(KIND=8)  tt393 
REAL(KIND=8)  tt394 
REAL(KIND=8)  tt395 
REAL(KIND=8)  tt396 
REAL(KIND=8)  tt397 
REAL(KIND=8)  tt398 
REAL(KIND=8)  tt399 
REAL(KIND=8)  tt400 
REAL(KIND=8)  tt401 
REAL(KIND=8)  tt402 
REAL(KIND=8)  tt403 
REAL(KIND=8)  tt404 
REAL(KIND=8)  tt405 
REAL(KIND=8)  tt406 
REAL(KIND=8)  tt407 
REAL(KIND=8)  tt408 
REAL(KIND=8)  tt409 
REAL(KIND=8)  tt410 
REAL(KIND=8)  tt411 
REAL(KIND=8)  tt412 
REAL(KIND=8)  tt413 
REAL(KIND=8)  tt414 
REAL(KIND=8)  tt415 
REAL(KIND=8)  tt416 
REAL(KIND=8)  tt417 
REAL(KIND=8)  tt418 
REAL(KIND=8)  tt419 
REAL(KIND=8)  tt420 
REAL(KIND=8)  tt421 
REAL(KIND=8)  tt422 
REAL(KIND=8)  tt423 
REAL(KIND=8)  tt424 
REAL(KIND=8)  tt425 
REAL(KIND=8)  tt426 
REAL(KIND=8)  tt427 
REAL(KIND=8)  tt428 
REAL(KIND=8)  tt429 
REAL(KIND=8)  tt430 
REAL(KIND=8)  tt431 
REAL(KIND=8)  tt432 
REAL(KIND=8)  tt433 
REAL(KIND=8)  tt434 
REAL(KIND=8)  tt435 
REAL(KIND=8)  tt436 
REAL(KIND=8)  tt437 
REAL(KIND=8)  tt438 
REAL(KIND=8)  tt439 
REAL(KIND=8)  tt440 
REAL(KIND=8)  tt441 
REAL(KIND=8)  tt442 
REAL(KIND=8)  tt443 
REAL(KIND=8)  tt444 
REAL(KIND=8)  tt445 
REAL(KIND=8)  tt446 
REAL(KIND=8)  tt447 
REAL(KIND=8)  tt448 
REAL(KIND=8)  tt449 
REAL(KIND=8)  tt450 
REAL(KIND=8)  tt451 
REAL(KIND=8)  tt452 
REAL(KIND=8)  tt453 
REAL(KIND=8)  tt454 
REAL(KIND=8)  tt455 
REAL(KIND=8)  tt456 
REAL(KIND=8)  tt457 
REAL(KIND=8)  tt458 
REAL(KIND=8)  tt459 
REAL(KIND=8)  tt460 
REAL(KIND=8)  tt461 
REAL(KIND=8)  tt462 
REAL(KIND=8)  tt463 
REAL(KIND=8)  tt464 
REAL(KIND=8)  tt465 
REAL(KIND=8)  tt466 
REAL(KIND=8)  tt467 
REAL(KIND=8)  tt468 
REAL(KIND=8)  tt469 
REAL(KIND=8)  tt470 
REAL(KIND=8)  tt471 
REAL(KIND=8)  tt472 
REAL(KIND=8)  tt473 
REAL(KIND=8)  tt474 
REAL(KIND=8)  tt475 
REAL(KIND=8)  tt476 
REAL(KIND=8)  tt477 
REAL(KIND=8)  tt478 
REAL(KIND=8)  tt479 
REAL(KIND=8)  tt480 
REAL(KIND=8)  tt481 
REAL(KIND=8)  tt482 
REAL(KIND=8)  tt483 
REAL(KIND=8)  tt484 
REAL(KIND=8)  tt485 
REAL(KIND=8)  tt486 
REAL(KIND=8)  tt487 
REAL(KIND=8)  tt488 
REAL(KIND=8)  tt489 
REAL(KIND=8)  tt490 
REAL(KIND=8)  tt491 
REAL(KIND=8)  tt492 
REAL(KIND=8)  tt493 
REAL(KIND=8)  tt494 
REAL(KIND=8)  tt495 
REAL(KIND=8)  tt496 
REAL(KIND=8)  tt497 
REAL(KIND=8)  tt498 
REAL(KIND=8)  tt499 
REAL(KIND=8)  tt500 
REAL(KIND=8)  tt501 
REAL(KIND=8)  tt502 
REAL(KIND=8)  tt503 
REAL(KIND=8)  tt504 
REAL(KIND=8)  tt505 
REAL(KIND=8)  tt506 
REAL(KIND=8)  tt507 
REAL(KIND=8)  tt508 
REAL(KIND=8)  tt509 
REAL(KIND=8)  tt510 
REAL(KIND=8)  tt511 
REAL(KIND=8)  tt512 
REAL(KIND=8)  tt513 
REAL(KIND=8)  tt514 
REAL(KIND=8)  tt515 
REAL(KIND=8)  tt516 
REAL(KIND=8)  tt517 
REAL(KIND=8)  tt518 
REAL(KIND=8)  tt519 
REAL(KIND=8)  tt520 
REAL(KIND=8)  tt521 
REAL(KIND=8)  tt522 
REAL(KIND=8)  tt523 
REAL(KIND=8)  tt524 
REAL(KIND=8)  tt525 
REAL(KIND=8)  tt526 
REAL(KIND=8)  tt527 
REAL(KIND=8)  tt528 
REAL(KIND=8)  tt529 
REAL(KIND=8)  tt530 
REAL(KIND=8)  tt531 
REAL(KIND=8)  tt532 
REAL(KIND=8)  tt533 
REAL(KIND=8)  tt534 
REAL(KIND=8)  tt535 
REAL(KIND=8)  tt536 
REAL(KIND=8)  tt537 
REAL(KIND=8)  tt538 
REAL(KIND=8)  tt539 
REAL(KIND=8)  tt540 
REAL(KIND=8)  tt541 
REAL(KIND=8)  tt542 
REAL(KIND=8)  tt543 
REAL(KIND=8)  tt544 
REAL(KIND=8)  tt545 
REAL(KIND=8)  tt546 
REAL(KIND=8)  tt547 
REAL(KIND=8)  tt548 
REAL(KIND=8)  tt549 
REAL(KIND=8)  tt550 
REAL(KIND=8)  tt551 
REAL(KIND=8)  tt552 
REAL(KIND=8)  tt553 
REAL(KIND=8)  tt554 
REAL(KIND=8)  tt555 
REAL(KIND=8)  tt556 
REAL(KIND=8)  tt557 
REAL(KIND=8)  tt558 
REAL(KIND=8)  tt559 
REAL(KIND=8)  tt560 
REAL(KIND=8)  tt561 
REAL(KIND=8)  tt562 
REAL(KIND=8)  tt563 
REAL(KIND=8)  tt564 
REAL(KIND=8)  tt565 
REAL(KIND=8)  tt566 
REAL(KIND=8)  tt567 
REAL(KIND=8)  tt568 
REAL(KIND=8)  tt569 
REAL(KIND=8)  tt570 
REAL(KIND=8)  tt571 
REAL(KIND=8)  tt572 
REAL(KIND=8)  tt573 
REAL(KIND=8)  tt574 
REAL(KIND=8)  tt575 
REAL(KIND=8)  tt576 
REAL(KIND=8)  tt577 
REAL(KIND=8)  tt578 
REAL(KIND=8)  tt579 
REAL(KIND=8)  tt580 
REAL(KIND=8)  tt581 
REAL(KIND=8)  tt582 
tt1 = (mu(1,1)*(2*H_invDmH(1,3)**2+2*H_invDmH(1,2)**2+2*H_invDmH(&
&1,1)**2))/2.0E+0
tt2 = mu(1,1)+lam(1,1)
tt3 = X(3,8)*H_invDmH(8,1)+X(3,7)*H_invDmH(7,1)+X(3,6)*H_invDmH(6&
&,1)+X(3,5)*H_invDmH(5,1)+X(3,4)*H_invDmH(4,1)+H_invDmH(3,1)*X(3,3&
&)+H_invDmH(2,1)*X(3,2)+H_invDmH(1,1)*X(3,1)
tt4 = X(2,8)*H_invDmH(8,2)+X(2,7)*H_invDmH(7,2)+X(2,6)*H_invDmH(6&
&,2)+X(2,5)*H_invDmH(5,2)+X(2,4)*H_invDmH(4,2)+X(2,3)*H_invDmH(3,2&
&)+H_invDmH(2,2)*X(2,2)+H_invDmH(1,2)*X(2,1)
tt5 = X(2,8)*H_invDmH(8,1)+X(2,7)*H_invDmH(7,1)+X(2,6)*H_invDmH(6&
&,1)+X(2,5)*H_invDmH(5,1)+X(2,4)*H_invDmH(4,1)+X(2,3)*H_invDmH(3,1&
&)+H_invDmH(2,1)*X(2,2)+H_invDmH(1,1)*X(2,1)
tt6 = X(3,8)*H_invDmH(8,2)+X(3,7)*H_invDmH(7,2)+X(3,6)*H_invDmH(6&
&,2)+X(3,5)*H_invDmH(5,2)+X(3,4)*H_invDmH(4,2)+H_invDmH(3,2)*X(3,3&
&)+H_invDmH(2,2)*X(3,2)+H_invDmH(1,2)*X(3,1)
tt7 = tt5*tt6-tt3*tt4
tt8 = X(2,8)*H_invDmH(8,3)+X(2,7)*H_invDmH(7,3)+X(2,6)*H_invDmH(6&
&,3)+X(2,5)*H_invDmH(5,3)+X(2,4)*H_invDmH(4,3)+X(2,3)*H_invDmH(3,3&
&)+X(2,2)*H_invDmH(2,3)+H_invDmH(1,3)*X(2,1)
tt9 = X(3,8)*H_invDmH(8,3)+X(3,7)*H_invDmH(7,3)+X(3,6)*H_invDmH(6&
&,3)+X(3,5)*H_invDmH(5,3)+X(3,4)*H_invDmH(4,3)+H_invDmH(3,3)*X(3,3&
&)+H_invDmH(2,3)*X(3,2)+H_invDmH(1,3)*X(3,1)
tt10 = tt5*tt9-tt3*tt8
tt11 = tt4*tt9-tt6*tt8
tt12 = H_invDmH(1,1)*tt11-H_invDmH(1,2)*tt10+H_invDmH(1,3)*tt7
tt13 = H_invDmH(1,1)*tt6-H_invDmH(1,2)*tt3
tt14 = H_invDmH(1,1)*tt9-H_invDmH(1,3)*tt3
tt15 = H_invDmH(1,2)*tt9-H_invDmH(1,3)*tt6
tt16 = H_invDmH(1,1)*tt15-H_invDmH(1,2)*tt14+H_invDmH(1,3)*tt13
tt17 = X(1,8)*H_invDmH(8,3)+X(1,7)*H_invDmH(7,3)+X(1,6)*H_invDmH(&
&6,3)+X(1,5)*H_invDmH(5,3)+X(1,4)*H_invDmH(4,3)+X(1,3)*H_invDmH(3,&
&3)+X(1,2)*H_invDmH(2,3)+X(1,1)*H_invDmH(1,3)
tt18 = X(1,8)*H_invDmH(8,2)+X(1,7)*H_invDmH(7,2)+X(1,6)*H_invDmH(&
&6,2)+X(1,5)*H_invDmH(5,2)+X(1,4)*H_invDmH(4,2)+X(1,3)*H_invDmH(3,&
&2)+X(1,2)*H_invDmH(2,2)+X(1,1)*H_invDmH(1,2)
tt19 = X(1,8)*H_invDmH(8,1)+X(1,7)*H_invDmH(7,1)+X(1,6)*H_invDmH(&
&6,1)+X(1,5)*H_invDmH(5,1)+X(1,4)*H_invDmH(4,1)+X(1,3)*H_invDmH(3,&
&1)+X(1,2)*H_invDmH(2,1)+H_invDmH(1,1)*X(1,1)
tt20 = tt19*tt15-tt18*tt14+tt13*tt17
tt21 = tt19*tt11-tt18*tt10+tt7*tt17-1
tt22 = detDmH(1,1)*gw(1,1)*(tt2*tt16*tt21+tt2*tt20*tt12-mu(1,1)*t&
&t16)
tt23 = H_invDmH(1,2)*tt5-H_invDmH(1,1)*tt4
tt24 = H_invDmH(1,3)*tt5-H_invDmH(1,1)*tt8
tt25 = H_invDmH(1,3)*tt4-H_invDmH(1,2)*tt8
tt26 = H_invDmH(1,1)*tt25-H_invDmH(1,2)*tt24+H_invDmH(1,3)*tt23
tt27 = tt19*tt25-tt18*tt24+tt23*tt17
tt28 = detDmH(1,1)*gw(1,1)*(tt2*tt26*tt21+tt2*tt27*tt12-mu(1,1)*t&
&t26)
tt29 = (mu(1,1)*(2*H_invDmH(1,3)*H_invDmH(2,3)+2*H_invDmH(1,2)*H_&
&invDmH(2,2)+2*H_invDmH(1,1)*H_invDmH(2,1)))/2.0E+0
tt30 = H_invDmH(2,1)*tt11-H_invDmH(2,2)*tt10+H_invDmH(2,3)*tt7
tt31 = detDmH(1,1)*gw(1,1)*(tt2*tt12*tt30+tt29)
tt32 = H_invDmH(2,1)*tt6-H_invDmH(2,2)*tt3
tt33 = H_invDmH(2,1)*tt9-H_invDmH(2,3)*tt3
tt34 = H_invDmH(2,2)*tt9-H_invDmH(2,3)*tt6
tt35 = H_invDmH(1,1)*tt34-H_invDmH(1,2)*tt33+H_invDmH(1,3)*tt32
tt36 = tt19*tt34-tt18*tt33+tt32*tt17
tt37 = detDmH(1,1)*gw(1,1)*(tt2*tt35*tt21+tt2*tt36*tt12-mu(1,1)*t&
&t35)
tt38 = H_invDmH(2,2)*tt5-H_invDmH(2,1)*tt4
tt39 = H_invDmH(2,3)*tt5-H_invDmH(2,1)*tt8
tt40 = H_invDmH(2,3)*tt4-H_invDmH(2,2)*tt8
tt41 = H_invDmH(1,1)*tt40-H_invDmH(1,2)*tt39+H_invDmH(1,3)*tt38
tt42 = tt19*tt40-tt18*tt39+tt38*tt17
tt43 = detDmH(1,1)*gw(1,1)*(tt2*tt41*tt21+tt2*tt42*tt12-mu(1,1)*t&
&t41)
tt44 = (mu(1,1)*(2*H_invDmH(1,3)*H_invDmH(3,3)+2*H_invDmH(1,2)*H_&
&invDmH(3,2)+2*H_invDmH(1,1)*H_invDmH(3,1)))/2.0E+0
tt45 = H_invDmH(3,1)*tt11-H_invDmH(3,2)*tt10+H_invDmH(3,3)*tt7
tt46 = detDmH(1,1)*gw(1,1)*(tt2*tt12*tt45+tt44)
tt47 = H_invDmH(3,1)*tt6-H_invDmH(3,2)*tt3
tt48 = H_invDmH(3,1)*tt9-H_invDmH(3,3)*tt3
tt49 = H_invDmH(3,2)*tt9-H_invDmH(3,3)*tt6
tt50 = H_invDmH(1,1)*tt49-H_invDmH(1,2)*tt48+H_invDmH(1,3)*tt47
tt51 = tt19*tt49-tt18*tt48+tt47*tt17
tt52 = detDmH(1,1)*gw(1,1)*(tt2*tt50*tt21+tt2*tt51*tt12-mu(1,1)*t&
&t50)
tt53 = H_invDmH(3,2)*tt5-H_invDmH(3,1)*tt4
tt54 = H_invDmH(3,3)*tt5-H_invDmH(3,1)*tt8
tt55 = H_invDmH(3,3)*tt4-H_invDmH(3,2)*tt8
tt56 = H_invDmH(1,1)*tt55-H_invDmH(1,2)*tt54+H_invDmH(1,3)*tt53
tt57 = tt19*tt55-tt18*tt54+tt53*tt17
tt58 = detDmH(1,1)*gw(1,1)*(tt2*tt56*tt21+tt2*tt57*tt12-mu(1,1)*t&
&t56)
tt59 = (mu(1,1)*(2*H_invDmH(1,3)*H_invDmH(4,3)+2*H_invDmH(1,2)*H_&
&invDmH(4,2)+2*H_invDmH(1,1)*H_invDmH(4,1)))/2.0E+0
tt60 = H_invDmH(4,1)*tt11-H_invDmH(4,2)*tt10+H_invDmH(4,3)*tt7
tt61 = detDmH(1,1)*gw(1,1)*(tt2*tt12*tt60+tt59)
tt62 = H_invDmH(4,1)*tt6-H_invDmH(4,2)*tt3
tt63 = H_invDmH(4,1)*tt9-H_invDmH(4,3)*tt3
tt64 = H_invDmH(4,2)*tt9-H_invDmH(4,3)*tt6
tt65 = H_invDmH(1,1)*tt64-H_invDmH(1,2)*tt63+H_invDmH(1,3)*tt62
tt66 = tt19*tt64-tt18*tt63+tt62*tt17
tt67 = detDmH(1,1)*gw(1,1)*(tt2*tt65*tt21+tt2*tt66*tt12-mu(1,1)*t&
&t65)
tt68 = H_invDmH(4,2)*tt5-H_invDmH(4,1)*tt4
tt69 = H_invDmH(4,3)*tt5-H_invDmH(4,1)*tt8
tt70 = H_invDmH(4,3)*tt4-H_invDmH(4,2)*tt8
tt71 = H_invDmH(1,1)*tt70-H_invDmH(1,2)*tt69+H_invDmH(1,3)*tt68
tt72 = tt19*tt70-tt18*tt69+tt68*tt17
tt73 = detDmH(1,1)*gw(1,1)*(tt2*tt71*tt21+tt2*tt72*tt12-mu(1,1)*t&
&t71)
tt74 = (mu(1,1)*(2*H_invDmH(1,3)*H_invDmH(5,3)+2*H_invDmH(1,2)*H_&
&invDmH(5,2)+2*H_invDmH(1,1)*H_invDmH(5,1)))/2.0E+0
tt75 = H_invDmH(5,1)*tt11-H_invDmH(5,2)*tt10+H_invDmH(5,3)*tt7
tt76 = detDmH(1,1)*gw(1,1)*(tt2*tt12*tt75+tt74)
tt77 = H_invDmH(5,1)*tt6-H_invDmH(5,2)*tt3
tt78 = H_invDmH(5,1)*tt9-H_invDmH(5,3)*tt3
tt79 = H_invDmH(5,2)*tt9-H_invDmH(5,3)*tt6
tt80 = H_invDmH(1,1)*tt79-H_invDmH(1,2)*tt78+H_invDmH(1,3)*tt77
tt81 = tt19*tt79-tt18*tt78+tt77*tt17
tt82 = detDmH(1,1)*gw(1,1)*(tt2*tt80*tt21+tt2*tt81*tt12-mu(1,1)*t&
&t80)
tt83 = H_invDmH(5,2)*tt5-H_invDmH(5,1)*tt4
tt84 = H_invDmH(5,3)*tt5-H_invDmH(5,1)*tt8
tt85 = H_invDmH(5,3)*tt4-H_invDmH(5,2)*tt8
tt86 = H_invDmH(1,1)*tt85-H_invDmH(1,2)*tt84+H_invDmH(1,3)*tt83
tt87 = tt19*tt85-tt18*tt84+tt83*tt17
tt88 = detDmH(1,1)*gw(1,1)*(tt2*tt86*tt21+tt2*tt87*tt12-mu(1,1)*t&
&t86)
tt89 = (mu(1,1)*(2*H_invDmH(1,3)*H_invDmH(6,3)+2*H_invDmH(1,2)*H_&
&invDmH(6,2)+2*H_invDmH(1,1)*H_invDmH(6,1)))/2.0E+0
tt90 = H_invDmH(6,1)*tt11-H_invDmH(6,2)*tt10+H_invDmH(6,3)*tt7
tt91 = detDmH(1,1)*gw(1,1)*(tt2*tt12*tt90+tt89)
tt92 = H_invDmH(6,1)*tt6-H_invDmH(6,2)*tt3
tt93 = H_invDmH(6,1)*tt9-H_invDmH(6,3)*tt3
tt94 = H_invDmH(6,2)*tt9-H_invDmH(6,3)*tt6
tt95 = H_invDmH(1,1)*tt94-H_invDmH(1,2)*tt93+H_invDmH(1,3)*tt92
tt96 = tt19*tt94-tt18*tt93+tt92*tt17
tt97 = detDmH(1,1)*gw(1,1)*(tt2*tt95*tt21+tt2*tt96*tt12-mu(1,1)*t&
&t95)
tt98 = H_invDmH(6,2)*tt5-H_invDmH(6,1)*tt4
tt99 = H_invDmH(6,3)*tt5-H_invDmH(6,1)*tt8
tt100 = H_invDmH(6,3)*tt4-H_invDmH(6,2)*tt8
tt101 = H_invDmH(1,1)*tt100-H_invDmH(1,2)*tt99+H_invDmH(1,3)*tt98&
&
tt102 = tt19*tt100-tt18*tt99+tt98*tt17
tt103 = detDmH(1,1)*gw(1,1)*(tt2*tt101*tt21+tt2*tt102*tt12-mu(1,1&
&)*tt101)
tt104 = (mu(1,1)*(2*H_invDmH(1,3)*H_invDmH(7,3)+2*H_invDmH(1,2)*H&
&_invDmH(7,2)+2*H_invDmH(1,1)*H_invDmH(7,1)))/2.0E+0
tt105 = H_invDmH(7,1)*tt11-H_invDmH(7,2)*tt10+H_invDmH(7,3)*tt7
tt106 = detDmH(1,1)*gw(1,1)*(tt2*tt12*tt105+tt104)
tt107 = H_invDmH(7,1)*tt6-H_invDmH(7,2)*tt3
tt108 = H_invDmH(7,1)*tt9-H_invDmH(7,3)*tt3
tt109 = H_invDmH(7,2)*tt9-H_invDmH(7,3)*tt6
tt110 = H_invDmH(1,1)*tt109-H_invDmH(1,2)*tt108+H_invDmH(1,3)*tt1&
&07
tt111 = tt19*tt109-tt18*tt108+tt107*tt17
tt112 = detDmH(1,1)*gw(1,1)*(tt2*tt110*tt21+tt2*tt111*tt12-mu(1,1&
&)*tt110)
tt113 = H_invDmH(7,2)*tt5-H_invDmH(7,1)*tt4
tt114 = H_invDmH(7,3)*tt5-H_invDmH(7,1)*tt8
tt115 = H_invDmH(7,3)*tt4-H_invDmH(7,2)*tt8
tt116 = H_invDmH(1,1)*tt115-H_invDmH(1,2)*tt114+H_invDmH(1,3)*tt1&
&13
tt117 = tt19*tt115-tt18*tt114+tt113*tt17
tt118 = detDmH(1,1)*gw(1,1)*(tt2*tt116*tt21+tt2*tt117*tt12-mu(1,1&
&)*tt116)
tt119 = (mu(1,1)*(2*H_invDmH(1,3)*H_invDmH(8,3)+2*H_invDmH(1,2)*H&
&_invDmH(8,2)+2*H_invDmH(1,1)*H_invDmH(8,1)))/2.0E+0
tt120 = H_invDmH(8,1)*tt11-H_invDmH(8,2)*tt10+tt7*H_invDmH(8,3)
tt121 = detDmH(1,1)*gw(1,1)*(tt2*tt12*tt120+tt119)
tt122 = H_invDmH(8,1)*tt6-tt3*H_invDmH(8,2)
tt123 = H_invDmH(8,1)*tt9-tt3*H_invDmH(8,3)
tt124 = H_invDmH(8,2)*tt9-tt6*H_invDmH(8,3)
tt125 = H_invDmH(1,1)*tt124-H_invDmH(1,2)*tt123+H_invDmH(1,3)*tt1&
&22
tt126 = tt19*tt124-tt18*tt123+tt122*tt17
tt127 = detDmH(1,1)*gw(1,1)*(tt2*tt125*tt21+tt2*tt126*tt12-mu(1,1&
&)*tt125)
tt128 = tt5*H_invDmH(8,2)-H_invDmH(8,1)*tt4
tt129 = tt5*H_invDmH(8,3)-H_invDmH(8,1)*tt8
tt130 = tt4*H_invDmH(8,3)-H_invDmH(8,2)*tt8
tt131 = H_invDmH(1,1)*tt130-H_invDmH(1,2)*tt129+H_invDmH(1,3)*tt1&
&28
tt132 = tt19*tt130-tt18*tt129+tt128*tt17
tt133 = detDmH(1,1)*gw(1,1)*(tt2*tt131*tt21+tt2*tt132*tt12-mu(1,1&
&)*tt131)
tt134 = detDmH(1,1)*gw(1,1)*tt2*tt27*tt20
tt135 = H_invDmH(2,1)*tt15-H_invDmH(2,2)*tt14+H_invDmH(2,3)*tt13
tt136 = detDmH(1,1)*gw(1,1)*(tt2*tt135*tt21+tt2*tt20*tt30-mu(1,1)&
&*tt135)
tt137 = detDmH(1,1)*gw(1,1)*(tt2*tt20*tt36+tt29)
tt138 = (H_invDmH(1,1)*H_invDmH(2,2)-H_invDmH(1,2)*H_invDmH(2,1))&
&*tt17-(H_invDmH(1,1)*H_invDmH(2,3)-H_invDmH(1,3)*H_invDmH(2,1))*t&
&t18+(H_invDmH(1,2)*H_invDmH(2,3)-H_invDmH(1,3)*H_invDmH(2,2))*tt1&
&9
tt139 = detDmH(1,1)*gw(1,1)*(tt2*tt138*tt21+tt2*tt42*tt20-mu(1,1)&
&*tt138)
tt140 = H_invDmH(3,1)*tt15-H_invDmH(3,2)*tt14+H_invDmH(3,3)*tt13
tt141 = detDmH(1,1)*gw(1,1)*(tt2*tt140*tt21+tt2*tt20*tt45-mu(1,1)&
&*tt140)
tt142 = detDmH(1,1)*gw(1,1)*(tt2*tt20*tt51+tt44)
tt143 = (H_invDmH(1,1)*H_invDmH(3,2)-H_invDmH(1,2)*H_invDmH(3,1))&
&*tt17-(H_invDmH(1,1)*H_invDmH(3,3)-H_invDmH(1,3)*H_invDmH(3,1))*t&
&t18+(H_invDmH(1,2)*H_invDmH(3,3)-H_invDmH(1,3)*H_invDmH(3,2))*tt1&
&9
tt144 = detDmH(1,1)*gw(1,1)*(tt2*tt143*tt21+tt2*tt57*tt20-mu(1,1)&
&*tt143)
tt145 = H_invDmH(4,1)*tt15-H_invDmH(4,2)*tt14+H_invDmH(4,3)*tt13
tt146 = detDmH(1,1)*gw(1,1)*(tt2*tt145*tt21+tt2*tt20*tt60-mu(1,1)&
&*tt145)
tt147 = detDmH(1,1)*gw(1,1)*(tt2*tt20*tt66+tt59)
tt148 = (H_invDmH(1,1)*H_invDmH(4,2)-H_invDmH(1,2)*H_invDmH(4,1))&
&*tt17-(H_invDmH(1,1)*H_invDmH(4,3)-H_invDmH(1,3)*H_invDmH(4,1))*t&
&t18+(H_invDmH(1,2)*H_invDmH(4,3)-H_invDmH(1,3)*H_invDmH(4,2))*tt1&
&9
tt149 = detDmH(1,1)*gw(1,1)*(tt2*tt148*tt21+tt2*tt72*tt20-mu(1,1)&
&*tt148)
tt150 = H_invDmH(5,1)*tt15-H_invDmH(5,2)*tt14+H_invDmH(5,3)*tt13
tt151 = detDmH(1,1)*gw(1,1)*(tt2*tt150*tt21+tt2*tt20*tt75-mu(1,1)&
&*tt150)
tt152 = detDmH(1,1)*gw(1,1)*(tt2*tt20*tt81+tt74)
tt153 = (H_invDmH(1,1)*H_invDmH(5,2)-H_invDmH(1,2)*H_invDmH(5,1))&
&*tt17-(H_invDmH(1,1)*H_invDmH(5,3)-H_invDmH(1,3)*H_invDmH(5,1))*t&
&t18+(H_invDmH(1,2)*H_invDmH(5,3)-H_invDmH(1,3)*H_invDmH(5,2))*tt1&
&9
tt154 = detDmH(1,1)*gw(1,1)*(tt2*tt153*tt21+tt2*tt87*tt20-mu(1,1)&
&*tt153)
tt155 = H_invDmH(6,1)*tt15-H_invDmH(6,2)*tt14+H_invDmH(6,3)*tt13
tt156 = detDmH(1,1)*gw(1,1)*(tt2*tt155*tt21+tt2*tt20*tt90-mu(1,1)&
&*tt155)
tt157 = detDmH(1,1)*gw(1,1)*(tt2*tt20*tt96+tt89)
tt158 = (H_invDmH(1,1)*H_invDmH(6,2)-H_invDmH(1,2)*H_invDmH(6,1))&
&*tt17-(H_invDmH(1,1)*H_invDmH(6,3)-H_invDmH(1,3)*H_invDmH(6,1))*t&
&t18+(H_invDmH(1,2)*H_invDmH(6,3)-H_invDmH(1,3)*H_invDmH(6,2))*tt1&
&9
tt159 = detDmH(1,1)*gw(1,1)*(tt2*tt158*tt21+tt2*tt102*tt20-mu(1,1&
&)*tt158)
tt160 = H_invDmH(7,1)*tt15-H_invDmH(7,2)*tt14+H_invDmH(7,3)*tt13
tt161 = detDmH(1,1)*gw(1,1)*(tt2*tt160*tt21+tt2*tt20*tt105-mu(1,1&
&)*tt160)
tt162 = detDmH(1,1)*gw(1,1)*(tt2*tt20*tt111+tt104)
tt163 = (H_invDmH(1,1)*H_invDmH(7,2)-H_invDmH(1,2)*H_invDmH(7,1))&
&*tt17-(H_invDmH(1,1)*H_invDmH(7,3)-H_invDmH(1,3)*H_invDmH(7,1))*t&
&t18+(H_invDmH(1,2)*H_invDmH(7,3)-H_invDmH(1,3)*H_invDmH(7,2))*tt1&
&9
tt164 = detDmH(1,1)*gw(1,1)*(tt2*tt163*tt21+tt2*tt117*tt20-mu(1,1&
&)*tt163)
tt165 = H_invDmH(8,1)*tt15-H_invDmH(8,2)*tt14+tt13*H_invDmH(8,3)
tt166 = detDmH(1,1)*gw(1,1)*(tt2*tt165*tt21+tt2*tt20*tt120-mu(1,1&
&)*tt165)
tt167 = detDmH(1,1)*gw(1,1)*(tt2*tt20*tt126+tt119)
tt168 = (H_invDmH(1,1)*H_invDmH(8,2)-H_invDmH(1,2)*H_invDmH(8,1))&
&*tt17+tt19*(H_invDmH(1,2)*H_invDmH(8,3)-H_invDmH(1,3)*H_invDmH(8,&
&2))-tt18*(H_invDmH(1,1)*H_invDmH(8,3)-H_invDmH(1,3)*H_invDmH(8,1)&
&)
tt169 = detDmH(1,1)*gw(1,1)*(tt2*tt168*tt21+tt2*tt132*tt20-mu(1,1&
&)*tt168)
tt170 = H_invDmH(2,1)*tt25-H_invDmH(2,2)*tt24+H_invDmH(2,3)*tt23
tt171 = detDmH(1,1)*gw(1,1)*(tt2*tt170*tt21+tt2*tt27*tt30-mu(1,1)&
&*tt170)
tt172 = (H_invDmH(1,2)*H_invDmH(2,1)-H_invDmH(1,1)*H_invDmH(2,2))&
&*tt17-(H_invDmH(1,3)*H_invDmH(2,1)-H_invDmH(1,1)*H_invDmH(2,3))*t&
&t18+(H_invDmH(1,3)*H_invDmH(2,2)-H_invDmH(1,2)*H_invDmH(2,3))*tt1&
&9
tt173 = detDmH(1,1)*gw(1,1)*(tt2*tt172*tt21+tt2*tt27*tt36-mu(1,1)&
&*tt172)
tt174 = detDmH(1,1)*gw(1,1)*(tt2*tt27*tt42+tt29)
tt175 = H_invDmH(3,1)*tt25-H_invDmH(3,2)*tt24+H_invDmH(3,3)*tt23
tt176 = detDmH(1,1)*gw(1,1)*(tt2*tt175*tt21+tt2*tt27*tt45-mu(1,1)&
&*tt175)
tt177 = (H_invDmH(1,2)*H_invDmH(3,1)-H_invDmH(1,1)*H_invDmH(3,2))&
&*tt17-(H_invDmH(1,3)*H_invDmH(3,1)-H_invDmH(1,1)*H_invDmH(3,3))*t&
&t18+(H_invDmH(1,3)*H_invDmH(3,2)-H_invDmH(1,2)*H_invDmH(3,3))*tt1&
&9
tt178 = detDmH(1,1)*gw(1,1)*(tt2*tt177*tt21+tt2*tt27*tt51-mu(1,1)&
&*tt177)
tt179 = detDmH(1,1)*gw(1,1)*(tt2*tt27*tt57+tt44)
tt180 = H_invDmH(4,1)*tt25-H_invDmH(4,2)*tt24+H_invDmH(4,3)*tt23
tt181 = detDmH(1,1)*gw(1,1)*(tt2*tt180*tt21+tt2*tt27*tt60-mu(1,1)&
&*tt180)
tt182 = (H_invDmH(1,2)*H_invDmH(4,1)-H_invDmH(1,1)*H_invDmH(4,2))&
&*tt17-(H_invDmH(1,3)*H_invDmH(4,1)-H_invDmH(1,1)*H_invDmH(4,3))*t&
&t18+(H_invDmH(1,3)*H_invDmH(4,2)-H_invDmH(1,2)*H_invDmH(4,3))*tt1&
&9
tt183 = detDmH(1,1)*gw(1,1)*(tt2*tt182*tt21+tt2*tt27*tt66-mu(1,1)&
&*tt182)
tt184 = detDmH(1,1)*gw(1,1)*(tt2*tt27*tt72+tt59)
tt185 = H_invDmH(5,1)*tt25-H_invDmH(5,2)*tt24+H_invDmH(5,3)*tt23
tt186 = detDmH(1,1)*gw(1,1)*(tt2*tt185*tt21+tt2*tt27*tt75-mu(1,1)&
&*tt185)
tt187 = (H_invDmH(1,2)*H_invDmH(5,1)-H_invDmH(1,1)*H_invDmH(5,2))&
&*tt17-(H_invDmH(1,3)*H_invDmH(5,1)-H_invDmH(1,1)*H_invDmH(5,3))*t&
&t18+(H_invDmH(1,3)*H_invDmH(5,2)-H_invDmH(1,2)*H_invDmH(5,3))*tt1&
&9
tt188 = detDmH(1,1)*gw(1,1)*(tt2*tt187*tt21+tt2*tt27*tt81-mu(1,1)&
&*tt187)
tt189 = detDmH(1,1)*gw(1,1)*(tt2*tt27*tt87+tt74)
tt190 = H_invDmH(6,1)*tt25-H_invDmH(6,2)*tt24+H_invDmH(6,3)*tt23
tt191 = detDmH(1,1)*gw(1,1)*(tt2*tt190*tt21+tt2*tt27*tt90-mu(1,1)&
&*tt190)
tt192 = (H_invDmH(1,2)*H_invDmH(6,1)-H_invDmH(1,1)*H_invDmH(6,2))&
&*tt17-(H_invDmH(1,3)*H_invDmH(6,1)-H_invDmH(1,1)*H_invDmH(6,3))*t&
&t18+(H_invDmH(1,3)*H_invDmH(6,2)-H_invDmH(1,2)*H_invDmH(6,3))*tt1&
&9
tt193 = detDmH(1,1)*gw(1,1)*(tt2*tt192*tt21+tt2*tt27*tt96-mu(1,1)&
&*tt192)
tt194 = detDmH(1,1)*gw(1,1)*(tt2*tt27*tt102+tt89)
tt195 = H_invDmH(7,1)*tt25-H_invDmH(7,2)*tt24+H_invDmH(7,3)*tt23
tt196 = detDmH(1,1)*gw(1,1)*(tt2*tt195*tt21+tt2*tt27*tt105-mu(1,1&
&)*tt195)
tt197 = (H_invDmH(1,2)*H_invDmH(7,1)-H_invDmH(1,1)*H_invDmH(7,2))&
&*tt17-(H_invDmH(1,3)*H_invDmH(7,1)-H_invDmH(1,1)*H_invDmH(7,3))*t&
&t18+(H_invDmH(1,3)*H_invDmH(7,2)-H_invDmH(1,2)*H_invDmH(7,3))*tt1&
&9
tt198 = detDmH(1,1)*gw(1,1)*(tt2*tt197*tt21+tt2*tt27*tt111-mu(1,1&
&)*tt197)
tt199 = detDmH(1,1)*gw(1,1)*(tt2*tt27*tt117+tt104)
tt200 = H_invDmH(8,1)*tt25-H_invDmH(8,2)*tt24+tt23*H_invDmH(8,3)
tt201 = detDmH(1,1)*gw(1,1)*(tt2*tt200*tt21+tt2*tt27*tt120-mu(1,1&
&)*tt200)
tt202 = (H_invDmH(1,2)*H_invDmH(8,1)-H_invDmH(1,1)*H_invDmH(8,2))&
&*tt17+tt19*(H_invDmH(1,3)*H_invDmH(8,2)-H_invDmH(1,2)*H_invDmH(8,&
&3))-tt18*(H_invDmH(1,3)*H_invDmH(8,1)-H_invDmH(1,1)*H_invDmH(8,3)&
&)
tt203 = detDmH(1,1)*gw(1,1)*(tt2*tt202*tt21+tt2*tt27*tt126-mu(1,1&
&)*tt202)
tt204 = detDmH(1,1)*gw(1,1)*(tt2*tt27*tt132+tt119)
tt205 = (mu(1,1)*(2*H_invDmH(2,3)**2+2*H_invDmH(2,2)**2+2*H_invDm&
&H(2,1)**2))/2.0E+0
tt206 = H_invDmH(2,1)*tt34-H_invDmH(2,2)*tt33+H_invDmH(2,3)*tt32
tt207 = detDmH(1,1)*gw(1,1)*(tt2*tt206*tt21+tt2*tt36*tt30-mu(1,1)&
&*tt206)
tt208 = H_invDmH(2,1)*tt40-H_invDmH(2,2)*tt39+H_invDmH(2,3)*tt38
tt209 = detDmH(1,1)*gw(1,1)*(tt2*tt208*tt21+tt2*tt42*tt30-mu(1,1)&
&*tt208)
tt210 = (mu(1,1)*(2*H_invDmH(2,3)*H_invDmH(3,3)+2*H_invDmH(2,2)*H&
&_invDmH(3,2)+2*H_invDmH(2,1)*H_invDmH(3,1)))/2.0E+0
tt211 = detDmH(1,1)*gw(1,1)*(tt2*tt30*tt45+tt210)
tt212 = H_invDmH(2,1)*tt49-H_invDmH(2,2)*tt48+H_invDmH(2,3)*tt47
tt213 = detDmH(1,1)*gw(1,1)*(tt2*tt212*tt21+tt2*tt51*tt30-mu(1,1)&
&*tt212)
tt214 = H_invDmH(2,1)*tt55-H_invDmH(2,2)*tt54+H_invDmH(2,3)*tt53
tt215 = detDmH(1,1)*gw(1,1)*(tt2*tt214*tt21+tt2*tt57*tt30-mu(1,1)&
&*tt214)
tt216 = (mu(1,1)*(2*H_invDmH(2,3)*H_invDmH(4,3)+2*H_invDmH(2,2)*H&
&_invDmH(4,2)+2*H_invDmH(2,1)*H_invDmH(4,1)))/2.0E+0
tt217 = detDmH(1,1)*gw(1,1)*(tt2*tt30*tt60+tt216)
tt218 = H_invDmH(2,1)*tt64-H_invDmH(2,2)*tt63+H_invDmH(2,3)*tt62
tt219 = detDmH(1,1)*gw(1,1)*(tt2*tt218*tt21+tt2*tt66*tt30-mu(1,1)&
&*tt218)
tt220 = H_invDmH(2,1)*tt70-H_invDmH(2,2)*tt69+H_invDmH(2,3)*tt68
tt221 = detDmH(1,1)*gw(1,1)*(tt2*tt220*tt21+tt2*tt72*tt30-mu(1,1)&
&*tt220)
tt222 = (mu(1,1)*(2*H_invDmH(2,3)*H_invDmH(5,3)+2*H_invDmH(2,2)*H&
&_invDmH(5,2)+2*H_invDmH(2,1)*H_invDmH(5,1)))/2.0E+0
tt223 = detDmH(1,1)*gw(1,1)*(tt2*tt30*tt75+tt222)
tt224 = H_invDmH(2,1)*tt79-H_invDmH(2,2)*tt78+H_invDmH(2,3)*tt77
tt225 = detDmH(1,1)*gw(1,1)*(tt2*tt224*tt21+tt2*tt81*tt30-mu(1,1)&
&*tt224)
tt226 = H_invDmH(2,1)*tt85-H_invDmH(2,2)*tt84+H_invDmH(2,3)*tt83
tt227 = detDmH(1,1)*gw(1,1)*(tt2*tt226*tt21+tt2*tt87*tt30-mu(1,1)&
&*tt226)
tt228 = (mu(1,1)*(2*H_invDmH(2,3)*H_invDmH(6,3)+2*H_invDmH(2,2)*H&
&_invDmH(6,2)+2*H_invDmH(2,1)*H_invDmH(6,1)))/2.0E+0
tt229 = detDmH(1,1)*gw(1,1)*(tt2*tt30*tt90+tt228)
tt230 = H_invDmH(2,1)*tt94-H_invDmH(2,2)*tt93+H_invDmH(2,3)*tt92
tt231 = detDmH(1,1)*gw(1,1)*(tt2*tt230*tt21+tt2*tt96*tt30-mu(1,1)&
&*tt230)
tt232 = H_invDmH(2,1)*tt100-H_invDmH(2,2)*tt99+H_invDmH(2,3)*tt98&
&
tt233 = detDmH(1,1)*gw(1,1)*(tt2*tt232*tt21+tt2*tt102*tt30-mu(1,1&
&)*tt232)
tt234 = (mu(1,1)*(2*H_invDmH(2,3)*H_invDmH(7,3)+2*H_invDmH(2,2)*H&
&_invDmH(7,2)+2*H_invDmH(2,1)*H_invDmH(7,1)))/2.0E+0
tt235 = detDmH(1,1)*gw(1,1)*(tt2*tt30*tt105+tt234)
tt236 = H_invDmH(2,1)*tt109-H_invDmH(2,2)*tt108+H_invDmH(2,3)*tt1&
&07
tt237 = detDmH(1,1)*gw(1,1)*(tt2*tt236*tt21+tt2*tt111*tt30-mu(1,1&
&)*tt236)
tt238 = H_invDmH(2,1)*tt115-H_invDmH(2,2)*tt114+H_invDmH(2,3)*tt1&
&13
tt239 = detDmH(1,1)*gw(1,1)*(tt2*tt238*tt21+tt2*tt117*tt30-mu(1,1&
&)*tt238)
tt240 = (mu(1,1)*(2*H_invDmH(2,3)*H_invDmH(8,3)+2*H_invDmH(2,2)*H&
&_invDmH(8,2)+2*H_invDmH(2,1)*H_invDmH(8,1)))/2.0E+0
tt241 = detDmH(1,1)*gw(1,1)*(tt2*tt30*tt120+tt240)
tt242 = H_invDmH(2,1)*tt124-H_invDmH(2,2)*tt123+H_invDmH(2,3)*tt1&
&22
tt243 = detDmH(1,1)*gw(1,1)*(tt2*tt242*tt21+tt2*tt126*tt30-mu(1,1&
&)*tt242)
tt244 = H_invDmH(2,1)*tt130-H_invDmH(2,2)*tt129+H_invDmH(2,3)*tt1&
&28
tt245 = detDmH(1,1)*gw(1,1)*(tt2*tt244*tt21+tt2*tt132*tt30-mu(1,1&
&)*tt244)
tt246 = detDmH(1,1)*gw(1,1)*tt2*tt42*tt36
tt247 = H_invDmH(3,1)*tt34-H_invDmH(3,2)*tt33+H_invDmH(3,3)*tt32
tt248 = detDmH(1,1)*gw(1,1)*(tt2*tt247*tt21+tt2*tt36*tt45-mu(1,1)&
&*tt247)
tt249 = detDmH(1,1)*gw(1,1)*(tt2*tt36*tt51+tt210)
tt250 = (H_invDmH(2,1)*H_invDmH(3,2)-H_invDmH(2,2)*H_invDmH(3,1))&
&*tt17-(H_invDmH(2,1)*H_invDmH(3,3)-H_invDmH(2,3)*H_invDmH(3,1))*t&
&t18+(H_invDmH(2,2)*H_invDmH(3,3)-H_invDmH(2,3)*H_invDmH(3,2))*tt1&
&9
tt251 = detDmH(1,1)*gw(1,1)*(tt2*tt250*tt21+tt2*tt57*tt36-mu(1,1)&
&*tt250)
tt252 = H_invDmH(4,1)*tt34-H_invDmH(4,2)*tt33+H_invDmH(4,3)*tt32
tt253 = detDmH(1,1)*gw(1,1)*(tt2*tt252*tt21+tt2*tt36*tt60-mu(1,1)&
&*tt252)
tt254 = detDmH(1,1)*gw(1,1)*(tt2*tt36*tt66+tt216)
tt255 = (H_invDmH(2,1)*H_invDmH(4,2)-H_invDmH(2,2)*H_invDmH(4,1))&
&*tt17-(H_invDmH(2,1)*H_invDmH(4,3)-H_invDmH(2,3)*H_invDmH(4,1))*t&
&t18+(H_invDmH(2,2)*H_invDmH(4,3)-H_invDmH(2,3)*H_invDmH(4,2))*tt1&
&9
tt256 = detDmH(1,1)*gw(1,1)*(tt2*tt255*tt21+tt2*tt72*tt36-mu(1,1)&
&*tt255)
tt257 = H_invDmH(5,1)*tt34-H_invDmH(5,2)*tt33+H_invDmH(5,3)*tt32
tt258 = detDmH(1,1)*gw(1,1)*(tt2*tt257*tt21+tt2*tt36*tt75-mu(1,1)&
&*tt257)
tt259 = detDmH(1,1)*gw(1,1)*(tt2*tt36*tt81+tt222)
tt260 = (H_invDmH(2,1)*H_invDmH(5,2)-H_invDmH(2,2)*H_invDmH(5,1))&
&*tt17-(H_invDmH(2,1)*H_invDmH(5,3)-H_invDmH(2,3)*H_invDmH(5,1))*t&
&t18+(H_invDmH(2,2)*H_invDmH(5,3)-H_invDmH(2,3)*H_invDmH(5,2))*tt1&
&9
tt261 = detDmH(1,1)*gw(1,1)*(tt2*tt260*tt21+tt2*tt87*tt36-mu(1,1)&
&*tt260)
tt262 = H_invDmH(6,1)*tt34-H_invDmH(6,2)*tt33+H_invDmH(6,3)*tt32
tt263 = detDmH(1,1)*gw(1,1)*(tt2*tt262*tt21+tt2*tt36*tt90-mu(1,1)&
&*tt262)
tt264 = detDmH(1,1)*gw(1,1)*(tt2*tt36*tt96+tt228)
tt265 = (H_invDmH(2,1)*H_invDmH(6,2)-H_invDmH(2,2)*H_invDmH(6,1))&
&*tt17-(H_invDmH(2,1)*H_invDmH(6,3)-H_invDmH(2,3)*H_invDmH(6,1))*t&
&t18+(H_invDmH(2,2)*H_invDmH(6,3)-H_invDmH(2,3)*H_invDmH(6,2))*tt1&
&9
tt266 = detDmH(1,1)*gw(1,1)*(tt2*tt265*tt21+tt2*tt102*tt36-mu(1,1&
&)*tt265)
tt267 = H_invDmH(7,1)*tt34-H_invDmH(7,2)*tt33+H_invDmH(7,3)*tt32
tt268 = detDmH(1,1)*gw(1,1)*(tt2*tt267*tt21+tt2*tt36*tt105-mu(1,1&
&)*tt267)
tt269 = detDmH(1,1)*gw(1,1)*(tt2*tt36*tt111+tt234)
tt270 = (H_invDmH(2,1)*H_invDmH(7,2)-H_invDmH(2,2)*H_invDmH(7,1))&
&*tt17-(H_invDmH(2,1)*H_invDmH(7,3)-H_invDmH(2,3)*H_invDmH(7,1))*t&
&t18+(H_invDmH(2,2)*H_invDmH(7,3)-H_invDmH(2,3)*H_invDmH(7,2))*tt1&
&9
tt271 = detDmH(1,1)*gw(1,1)*(tt2*tt270*tt21+tt2*tt117*tt36-mu(1,1&
&)*tt270)
tt272 = H_invDmH(8,1)*tt34-H_invDmH(8,2)*tt33+tt32*H_invDmH(8,3)
tt273 = detDmH(1,1)*gw(1,1)*(tt2*tt272*tt21+tt2*tt36*tt120-mu(1,1&
&)*tt272)
tt274 = detDmH(1,1)*gw(1,1)*(tt2*tt36*tt126+tt240)
tt275 = tt19*(H_invDmH(2,2)*H_invDmH(8,3)-H_invDmH(2,3)*H_invDmH(&
&8,2))-tt18*(H_invDmH(2,1)*H_invDmH(8,3)-H_invDmH(2,3)*H_invDmH(8,&
&1))+(H_invDmH(2,1)*H_invDmH(8,2)-H_invDmH(2,2)*H_invDmH(8,1))*tt1&
&7
tt276 = detDmH(1,1)*gw(1,1)*(tt2*tt275*tt21+tt2*tt132*tt36-mu(1,1&
&)*tt275)
tt277 = H_invDmH(3,1)*tt40-H_invDmH(3,2)*tt39+H_invDmH(3,3)*tt38
tt278 = detDmH(1,1)*gw(1,1)*(tt2*tt277*tt21+tt2*tt42*tt45-mu(1,1)&
&*tt277)
tt279 = (H_invDmH(2,2)*H_invDmH(3,1)-H_invDmH(2,1)*H_invDmH(3,2))&
&*tt17-(H_invDmH(2,3)*H_invDmH(3,1)-H_invDmH(2,1)*H_invDmH(3,3))*t&
&t18+(H_invDmH(2,3)*H_invDmH(3,2)-H_invDmH(2,2)*H_invDmH(3,3))*tt1&
&9
tt280 = detDmH(1,1)*gw(1,1)*(tt2*tt279*tt21+tt2*tt42*tt51-mu(1,1)&
&*tt279)
tt281 = detDmH(1,1)*gw(1,1)*(tt2*tt42*tt57+tt210)
tt282 = H_invDmH(4,1)*tt40-H_invDmH(4,2)*tt39+H_invDmH(4,3)*tt38
tt283 = detDmH(1,1)*gw(1,1)*(tt2*tt282*tt21+tt2*tt42*tt60-mu(1,1)&
&*tt282)
tt284 = (H_invDmH(2,2)*H_invDmH(4,1)-H_invDmH(2,1)*H_invDmH(4,2))&
&*tt17-(H_invDmH(2,3)*H_invDmH(4,1)-H_invDmH(2,1)*H_invDmH(4,3))*t&
&t18+(H_invDmH(2,3)*H_invDmH(4,2)-H_invDmH(2,2)*H_invDmH(4,3))*tt1&
&9
tt285 = detDmH(1,1)*gw(1,1)*(tt2*tt284*tt21+tt2*tt42*tt66-mu(1,1)&
&*tt284)
tt286 = detDmH(1,1)*gw(1,1)*(tt2*tt42*tt72+tt216)
tt287 = H_invDmH(5,1)*tt40-H_invDmH(5,2)*tt39+H_invDmH(5,3)*tt38
tt288 = detDmH(1,1)*gw(1,1)*(tt2*tt287*tt21+tt2*tt42*tt75-mu(1,1)&
&*tt287)
tt289 = (H_invDmH(2,2)*H_invDmH(5,1)-H_invDmH(2,1)*H_invDmH(5,2))&
&*tt17-(H_invDmH(2,3)*H_invDmH(5,1)-H_invDmH(2,1)*H_invDmH(5,3))*t&
&t18+(H_invDmH(2,3)*H_invDmH(5,2)-H_invDmH(2,2)*H_invDmH(5,3))*tt1&
&9
tt290 = detDmH(1,1)*gw(1,1)*(tt2*tt289*tt21+tt2*tt42*tt81-mu(1,1)&
&*tt289)
tt291 = detDmH(1,1)*gw(1,1)*(tt2*tt42*tt87+tt222)
tt292 = H_invDmH(6,1)*tt40-H_invDmH(6,2)*tt39+H_invDmH(6,3)*tt38
tt293 = detDmH(1,1)*gw(1,1)*(tt2*tt292*tt21+tt2*tt42*tt90-mu(1,1)&
&*tt292)
tt294 = (H_invDmH(2,2)*H_invDmH(6,1)-H_invDmH(2,1)*H_invDmH(6,2))&
&*tt17-(H_invDmH(2,3)*H_invDmH(6,1)-H_invDmH(2,1)*H_invDmH(6,3))*t&
&t18+(H_invDmH(2,3)*H_invDmH(6,2)-H_invDmH(2,2)*H_invDmH(6,3))*tt1&
&9
tt295 = detDmH(1,1)*gw(1,1)*(tt2*tt294*tt21+tt2*tt42*tt96-mu(1,1)&
&*tt294)
tt296 = detDmH(1,1)*gw(1,1)*(tt2*tt42*tt102+tt228)
tt297 = H_invDmH(7,1)*tt40-H_invDmH(7,2)*tt39+H_invDmH(7,3)*tt38
tt298 = detDmH(1,1)*gw(1,1)*(tt2*tt297*tt21+tt2*tt42*tt105-mu(1,1&
&)*tt297)
tt299 = (H_invDmH(2,2)*H_invDmH(7,1)-H_invDmH(2,1)*H_invDmH(7,2))&
&*tt17-(H_invDmH(2,3)*H_invDmH(7,1)-H_invDmH(2,1)*H_invDmH(7,3))*t&
&t18+(H_invDmH(2,3)*H_invDmH(7,2)-H_invDmH(2,2)*H_invDmH(7,3))*tt1&
&9
tt300 = detDmH(1,1)*gw(1,1)*(tt2*tt299*tt21+tt2*tt42*tt111-mu(1,1&
&)*tt299)
tt301 = detDmH(1,1)*gw(1,1)*(tt2*tt42*tt117+tt234)
tt302 = H_invDmH(8,1)*tt40-H_invDmH(8,2)*tt39+tt38*H_invDmH(8,3)
tt303 = detDmH(1,1)*gw(1,1)*(tt2*tt302*tt21+tt2*tt42*tt120-mu(1,1&
&)*tt302)
tt304 = tt19*(H_invDmH(2,3)*H_invDmH(8,2)-H_invDmH(2,2)*H_invDmH(&
&8,3))-tt18*(H_invDmH(2,3)*H_invDmH(8,1)-H_invDmH(2,1)*H_invDmH(8,&
&3))+(H_invDmH(2,2)*H_invDmH(8,1)-H_invDmH(2,1)*H_invDmH(8,2))*tt1&
&7
tt305 = detDmH(1,1)*gw(1,1)*(tt2*tt304*tt21+tt2*tt42*tt126-mu(1,1&
&)*tt304)
tt306 = detDmH(1,1)*gw(1,1)*(tt2*tt42*tt132+tt240)
tt307 = (mu(1,1)*(2*H_invDmH(3,3)**2+2*H_invDmH(3,2)**2+2*H_invDm&
&H(3,1)**2))/2.0E+0
tt308 = H_invDmH(3,1)*tt49-H_invDmH(3,2)*tt48+H_invDmH(3,3)*tt47
tt309 = detDmH(1,1)*gw(1,1)*(tt2*tt308*tt21+tt2*tt51*tt45-mu(1,1)&
&*tt308)
tt310 = H_invDmH(3,1)*tt55-H_invDmH(3,2)*tt54+H_invDmH(3,3)*tt53
tt311 = detDmH(1,1)*gw(1,1)*(tt2*tt310*tt21+tt2*tt57*tt45-mu(1,1)&
&*tt310)
tt312 = (mu(1,1)*(2*H_invDmH(3,3)*H_invDmH(4,3)+2*H_invDmH(3,2)*H&
&_invDmH(4,2)+2*H_invDmH(3,1)*H_invDmH(4,1)))/2.0E+0
tt313 = detDmH(1,1)*gw(1,1)*(tt2*tt45*tt60+tt312)
tt314 = H_invDmH(3,1)*tt64-H_invDmH(3,2)*tt63+H_invDmH(3,3)*tt62
tt315 = detDmH(1,1)*gw(1,1)*(tt2*tt314*tt21+tt2*tt66*tt45-mu(1,1)&
&*tt314)
tt316 = H_invDmH(3,1)*tt70-H_invDmH(3,2)*tt69+H_invDmH(3,3)*tt68
tt317 = detDmH(1,1)*gw(1,1)*(tt2*tt316*tt21+tt2*tt72*tt45-mu(1,1)&
&*tt316)
tt318 = (mu(1,1)*(2*H_invDmH(3,3)*H_invDmH(5,3)+2*H_invDmH(3,2)*H&
&_invDmH(5,2)+2*H_invDmH(3,1)*H_invDmH(5,1)))/2.0E+0
tt319 = detDmH(1,1)*gw(1,1)*(tt2*tt45*tt75+tt318)
tt320 = H_invDmH(3,1)*tt79-H_invDmH(3,2)*tt78+H_invDmH(3,3)*tt77
tt321 = detDmH(1,1)*gw(1,1)*(tt2*tt320*tt21+tt2*tt81*tt45-mu(1,1)&
&*tt320)
tt322 = H_invDmH(3,1)*tt85-H_invDmH(3,2)*tt84+H_invDmH(3,3)*tt83
tt323 = detDmH(1,1)*gw(1,1)*(tt2*tt322*tt21+tt2*tt87*tt45-mu(1,1)&
&*tt322)
tt324 = (mu(1,1)*(2*H_invDmH(3,3)*H_invDmH(6,3)+2*H_invDmH(3,2)*H&
&_invDmH(6,2)+2*H_invDmH(3,1)*H_invDmH(6,1)))/2.0E+0
tt325 = detDmH(1,1)*gw(1,1)*(tt2*tt45*tt90+tt324)
tt326 = H_invDmH(3,1)*tt94-H_invDmH(3,2)*tt93+H_invDmH(3,3)*tt92
tt327 = detDmH(1,1)*gw(1,1)*(tt2*tt326*tt21+tt2*tt96*tt45-mu(1,1)&
&*tt326)
tt328 = H_invDmH(3,1)*tt100-H_invDmH(3,2)*tt99+H_invDmH(3,3)*tt98&
&
tt329 = detDmH(1,1)*gw(1,1)*(tt2*tt328*tt21+tt2*tt102*tt45-mu(1,1&
&)*tt328)
tt330 = (mu(1,1)*(2*H_invDmH(3,3)*H_invDmH(7,3)+2*H_invDmH(3,2)*H&
&_invDmH(7,2)+2*H_invDmH(3,1)*H_invDmH(7,1)))/2.0E+0
tt331 = detDmH(1,1)*gw(1,1)*(tt2*tt45*tt105+tt330)
tt332 = H_invDmH(3,1)*tt109-H_invDmH(3,2)*tt108+H_invDmH(3,3)*tt1&
&07
tt333 = detDmH(1,1)*gw(1,1)*(tt2*tt332*tt21+tt2*tt111*tt45-mu(1,1&
&)*tt332)
tt334 = H_invDmH(3,1)*tt115-H_invDmH(3,2)*tt114+H_invDmH(3,3)*tt1&
&13
tt335 = detDmH(1,1)*gw(1,1)*(tt2*tt334*tt21+tt2*tt117*tt45-mu(1,1&
&)*tt334)
tt336 = (mu(1,1)*(2*H_invDmH(3,3)*H_invDmH(8,3)+2*H_invDmH(3,2)*H&
&_invDmH(8,2)+2*H_invDmH(3,1)*H_invDmH(8,1)))/2.0E+0
tt337 = detDmH(1,1)*gw(1,1)*(tt2*tt45*tt120+tt336)
tt338 = H_invDmH(3,1)*tt124-H_invDmH(3,2)*tt123+H_invDmH(3,3)*tt1&
&22
tt339 = detDmH(1,1)*gw(1,1)*(tt2*tt338*tt21+tt2*tt126*tt45-mu(1,1&
&)*tt338)
tt340 = H_invDmH(3,1)*tt130-H_invDmH(3,2)*tt129+H_invDmH(3,3)*tt1&
&28
tt341 = detDmH(1,1)*gw(1,1)*(tt2*tt340*tt21+tt2*tt132*tt45-mu(1,1&
&)*tt340)
tt342 = detDmH(1,1)*gw(1,1)*tt2*tt57*tt51
tt343 = H_invDmH(4,1)*tt49-H_invDmH(4,2)*tt48+H_invDmH(4,3)*tt47
tt344 = detDmH(1,1)*gw(1,1)*(tt2*tt343*tt21+tt2*tt51*tt60-mu(1,1)&
&*tt343)
tt345 = detDmH(1,1)*gw(1,1)*(tt2*tt51*tt66+tt312)
tt346 = (H_invDmH(3,1)*H_invDmH(4,2)-H_invDmH(3,2)*H_invDmH(4,1))&
&*tt17-(H_invDmH(3,1)*H_invDmH(4,3)-H_invDmH(3,3)*H_invDmH(4,1))*t&
&t18+(H_invDmH(3,2)*H_invDmH(4,3)-H_invDmH(3,3)*H_invDmH(4,2))*tt1&
&9
tt347 = detDmH(1,1)*gw(1,1)*(tt2*tt346*tt21+tt2*tt72*tt51-mu(1,1)&
&*tt346)
tt348 = H_invDmH(5,1)*tt49-H_invDmH(5,2)*tt48+H_invDmH(5,3)*tt47
tt349 = detDmH(1,1)*gw(1,1)*(tt2*tt348*tt21+tt2*tt51*tt75-mu(1,1)&
&*tt348)
tt350 = detDmH(1,1)*gw(1,1)*(tt2*tt51*tt81+tt318)
tt351 = (H_invDmH(3,1)*H_invDmH(5,2)-H_invDmH(3,2)*H_invDmH(5,1))&
&*tt17-(H_invDmH(3,1)*H_invDmH(5,3)-H_invDmH(3,3)*H_invDmH(5,1))*t&
&t18+(H_invDmH(3,2)*H_invDmH(5,3)-H_invDmH(3,3)*H_invDmH(5,2))*tt1&
&9
tt352 = detDmH(1,1)*gw(1,1)*(tt2*tt351*tt21+tt2*tt87*tt51-mu(1,1)&
&*tt351)
tt353 = H_invDmH(6,1)*tt49-H_invDmH(6,2)*tt48+H_invDmH(6,3)*tt47
tt354 = detDmH(1,1)*gw(1,1)*(tt2*tt353*tt21+tt2*tt51*tt90-mu(1,1)&
&*tt353)
tt355 = detDmH(1,1)*gw(1,1)*(tt2*tt51*tt96+tt324)
tt356 = (H_invDmH(3,1)*H_invDmH(6,2)-H_invDmH(3,2)*H_invDmH(6,1))&
&*tt17-(H_invDmH(3,1)*H_invDmH(6,3)-H_invDmH(3,3)*H_invDmH(6,1))*t&
&t18+(H_invDmH(3,2)*H_invDmH(6,3)-H_invDmH(3,3)*H_invDmH(6,2))*tt1&
&9
tt357 = detDmH(1,1)*gw(1,1)*(tt2*tt356*tt21+tt2*tt102*tt51-mu(1,1&
&)*tt356)
tt358 = H_invDmH(7,1)*tt49-H_invDmH(7,2)*tt48+H_invDmH(7,3)*tt47
tt359 = detDmH(1,1)*gw(1,1)*(tt2*tt358*tt21+tt2*tt51*tt105-mu(1,1&
&)*tt358)
tt360 = detDmH(1,1)*gw(1,1)*(tt2*tt51*tt111+tt330)
tt361 = (H_invDmH(3,1)*H_invDmH(7,2)-H_invDmH(3,2)*H_invDmH(7,1))&
&*tt17-(H_invDmH(3,1)*H_invDmH(7,3)-H_invDmH(3,3)*H_invDmH(7,1))*t&
&t18+(H_invDmH(3,2)*H_invDmH(7,3)-H_invDmH(3,3)*H_invDmH(7,2))*tt1&
&9
tt362 = detDmH(1,1)*gw(1,1)*(tt2*tt361*tt21+tt2*tt117*tt51-mu(1,1&
&)*tt361)
tt363 = H_invDmH(8,1)*tt49-H_invDmH(8,2)*tt48+tt47*H_invDmH(8,3)
tt364 = detDmH(1,1)*gw(1,1)*(tt2*tt363*tt21+tt2*tt51*tt120-mu(1,1&
&)*tt363)
tt365 = detDmH(1,1)*gw(1,1)*(tt2*tt51*tt126+tt336)
tt366 = tt19*(H_invDmH(3,2)*H_invDmH(8,3)-H_invDmH(3,3)*H_invDmH(&
&8,2))-tt18*(H_invDmH(3,1)*H_invDmH(8,3)-H_invDmH(3,3)*H_invDmH(8,&
&1))+(H_invDmH(3,1)*H_invDmH(8,2)-H_invDmH(3,2)*H_invDmH(8,1))*tt1&
&7
tt367 = detDmH(1,1)*gw(1,1)*(tt2*tt366*tt21+tt2*tt132*tt51-mu(1,1&
&)*tt366)
tt368 = H_invDmH(4,1)*tt55-H_invDmH(4,2)*tt54+H_invDmH(4,3)*tt53
tt369 = detDmH(1,1)*gw(1,1)*(tt2*tt368*tt21+tt2*tt57*tt60-mu(1,1)&
&*tt368)
tt370 = (H_invDmH(3,2)*H_invDmH(4,1)-H_invDmH(3,1)*H_invDmH(4,2))&
&*tt17-(H_invDmH(3,3)*H_invDmH(4,1)-H_invDmH(3,1)*H_invDmH(4,3))*t&
&t18+(H_invDmH(3,3)*H_invDmH(4,2)-H_invDmH(3,2)*H_invDmH(4,3))*tt1&
&9
tt371 = detDmH(1,1)*gw(1,1)*(tt2*tt370*tt21+tt2*tt57*tt66-mu(1,1)&
&*tt370)
tt372 = detDmH(1,1)*gw(1,1)*(tt2*tt57*tt72+tt312)
tt373 = H_invDmH(5,1)*tt55-H_invDmH(5,2)*tt54+H_invDmH(5,3)*tt53
tt374 = detDmH(1,1)*gw(1,1)*(tt2*tt373*tt21+tt2*tt57*tt75-mu(1,1)&
&*tt373)
tt375 = (H_invDmH(3,2)*H_invDmH(5,1)-H_invDmH(3,1)*H_invDmH(5,2))&
&*tt17-(H_invDmH(3,3)*H_invDmH(5,1)-H_invDmH(3,1)*H_invDmH(5,3))*t&
&t18+(H_invDmH(3,3)*H_invDmH(5,2)-H_invDmH(3,2)*H_invDmH(5,3))*tt1&
&9
tt376 = detDmH(1,1)*gw(1,1)*(tt2*tt375*tt21+tt2*tt57*tt81-mu(1,1)&
&*tt375)
tt377 = detDmH(1,1)*gw(1,1)*(tt2*tt57*tt87+tt318)
tt378 = H_invDmH(6,1)*tt55-H_invDmH(6,2)*tt54+H_invDmH(6,3)*tt53
tt379 = detDmH(1,1)*gw(1,1)*(tt2*tt378*tt21+tt2*tt57*tt90-mu(1,1)&
&*tt378)
tt380 = (H_invDmH(3,2)*H_invDmH(6,1)-H_invDmH(3,1)*H_invDmH(6,2))&
&*tt17-(H_invDmH(3,3)*H_invDmH(6,1)-H_invDmH(3,1)*H_invDmH(6,3))*t&
&t18+(H_invDmH(3,3)*H_invDmH(6,2)-H_invDmH(3,2)*H_invDmH(6,3))*tt1&
&9
tt381 = detDmH(1,1)*gw(1,1)*(tt2*tt380*tt21+tt2*tt57*tt96-mu(1,1)&
&*tt380)
tt382 = detDmH(1,1)*gw(1,1)*(tt2*tt57*tt102+tt324)
tt383 = H_invDmH(7,1)*tt55-H_invDmH(7,2)*tt54+H_invDmH(7,3)*tt53
tt384 = detDmH(1,1)*gw(1,1)*(tt2*tt383*tt21+tt2*tt57*tt105-mu(1,1&
&)*tt383)
tt385 = (H_invDmH(3,2)*H_invDmH(7,1)-H_invDmH(3,1)*H_invDmH(7,2))&
&*tt17-(H_invDmH(3,3)*H_invDmH(7,1)-H_invDmH(3,1)*H_invDmH(7,3))*t&
&t18+(H_invDmH(3,3)*H_invDmH(7,2)-H_invDmH(3,2)*H_invDmH(7,3))*tt1&
&9
tt386 = detDmH(1,1)*gw(1,1)*(tt2*tt385*tt21+tt2*tt57*tt111-mu(1,1&
&)*tt385)
tt387 = detDmH(1,1)*gw(1,1)*(tt2*tt57*tt117+tt330)
tt388 = H_invDmH(8,1)*tt55-H_invDmH(8,2)*tt54+tt53*H_invDmH(8,3)
tt389 = detDmH(1,1)*gw(1,1)*(tt2*tt388*tt21+tt2*tt57*tt120-mu(1,1&
&)*tt388)
tt390 = tt19*(H_invDmH(3,3)*H_invDmH(8,2)-H_invDmH(3,2)*H_invDmH(&
&8,3))-tt18*(H_invDmH(3,3)*H_invDmH(8,1)-H_invDmH(3,1)*H_invDmH(8,&
&3))+(H_invDmH(3,2)*H_invDmH(8,1)-H_invDmH(3,1)*H_invDmH(8,2))*tt1&
&7
tt391 = detDmH(1,1)*gw(1,1)*(tt2*tt390*tt21+tt2*tt57*tt126-mu(1,1&
&)*tt390)
tt392 = detDmH(1,1)*gw(1,1)*(tt2*tt57*tt132+tt336)
tt393 = (mu(1,1)*(2*H_invDmH(4,3)**2+2*H_invDmH(4,2)**2+2*H_invDm&
&H(4,1)**2))/2.0E+0
tt394 = H_invDmH(4,1)*tt64-H_invDmH(4,2)*tt63+H_invDmH(4,3)*tt62
tt395 = detDmH(1,1)*gw(1,1)*(tt2*tt394*tt21+tt2*tt66*tt60-mu(1,1)&
&*tt394)
tt396 = H_invDmH(4,1)*tt70-H_invDmH(4,2)*tt69+H_invDmH(4,3)*tt68
tt397 = detDmH(1,1)*gw(1,1)*(tt2*tt396*tt21+tt2*tt72*tt60-mu(1,1)&
&*tt396)
tt398 = (mu(1,1)*(2*H_invDmH(4,3)*H_invDmH(5,3)+2*H_invDmH(4,2)*H&
&_invDmH(5,2)+2*H_invDmH(4,1)*H_invDmH(5,1)))/2.0E+0
tt399 = detDmH(1,1)*gw(1,1)*(tt2*tt60*tt75+tt398)
tt400 = H_invDmH(4,1)*tt79-H_invDmH(4,2)*tt78+H_invDmH(4,3)*tt77
tt401 = detDmH(1,1)*gw(1,1)*(tt2*tt400*tt21+tt2*tt81*tt60-mu(1,1)&
&*tt400)
tt402 = H_invDmH(4,1)*tt85-H_invDmH(4,2)*tt84+H_invDmH(4,3)*tt83
tt403 = detDmH(1,1)*gw(1,1)*(tt2*tt402*tt21+tt2*tt87*tt60-mu(1,1)&
&*tt402)
tt404 = (mu(1,1)*(2*H_invDmH(4,3)*H_invDmH(6,3)+2*H_invDmH(4,2)*H&
&_invDmH(6,2)+2*H_invDmH(4,1)*H_invDmH(6,1)))/2.0E+0
tt405 = detDmH(1,1)*gw(1,1)*(tt2*tt60*tt90+tt404)
tt406 = H_invDmH(4,1)*tt94-H_invDmH(4,2)*tt93+H_invDmH(4,3)*tt92
tt407 = detDmH(1,1)*gw(1,1)*(tt2*tt406*tt21+tt2*tt96*tt60-mu(1,1)&
&*tt406)
tt408 = H_invDmH(4,1)*tt100-H_invDmH(4,2)*tt99+H_invDmH(4,3)*tt98&
&
tt409 = detDmH(1,1)*gw(1,1)*(tt2*tt408*tt21+tt2*tt102*tt60-mu(1,1&
&)*tt408)
tt410 = (mu(1,1)*(2*H_invDmH(4,3)*H_invDmH(7,3)+2*H_invDmH(4,2)*H&
&_invDmH(7,2)+2*H_invDmH(4,1)*H_invDmH(7,1)))/2.0E+0
tt411 = detDmH(1,1)*gw(1,1)*(tt2*tt60*tt105+tt410)
tt412 = H_invDmH(4,1)*tt109-H_invDmH(4,2)*tt108+H_invDmH(4,3)*tt1&
&07
tt413 = detDmH(1,1)*gw(1,1)*(tt2*tt412*tt21+tt2*tt111*tt60-mu(1,1&
&)*tt412)
tt414 = H_invDmH(4,1)*tt115-H_invDmH(4,2)*tt114+H_invDmH(4,3)*tt1&
&13
tt415 = detDmH(1,1)*gw(1,1)*(tt2*tt414*tt21+tt2*tt117*tt60-mu(1,1&
&)*tt414)
tt416 = (mu(1,1)*(2*H_invDmH(4,3)*H_invDmH(8,3)+2*H_invDmH(4,2)*H&
&_invDmH(8,2)+2*H_invDmH(4,1)*H_invDmH(8,1)))/2.0E+0
tt417 = detDmH(1,1)*gw(1,1)*(tt2*tt60*tt120+tt416)
tt418 = H_invDmH(4,1)*tt124-H_invDmH(4,2)*tt123+H_invDmH(4,3)*tt1&
&22
tt419 = detDmH(1,1)*gw(1,1)*(tt2*tt418*tt21+tt2*tt126*tt60-mu(1,1&
&)*tt418)
tt420 = H_invDmH(4,1)*tt130-H_invDmH(4,2)*tt129+H_invDmH(4,3)*tt1&
&28
tt421 = detDmH(1,1)*gw(1,1)*(tt2*tt420*tt21+tt2*tt132*tt60-mu(1,1&
&)*tt420)
tt422 = detDmH(1,1)*gw(1,1)*tt2*tt72*tt66
tt423 = H_invDmH(5,1)*tt64-H_invDmH(5,2)*tt63+H_invDmH(5,3)*tt62
tt424 = detDmH(1,1)*gw(1,1)*(tt2*tt423*tt21+tt2*tt66*tt75-mu(1,1)&
&*tt423)
tt425 = detDmH(1,1)*gw(1,1)*(tt2*tt66*tt81+tt398)
tt426 = (H_invDmH(4,1)*H_invDmH(5,2)-H_invDmH(4,2)*H_invDmH(5,1))&
&*tt17-(H_invDmH(4,1)*H_invDmH(5,3)-H_invDmH(4,3)*H_invDmH(5,1))*t&
&t18+(H_invDmH(4,2)*H_invDmH(5,3)-H_invDmH(4,3)*H_invDmH(5,2))*tt1&
&9
tt427 = detDmH(1,1)*gw(1,1)*(tt2*tt426*tt21+tt2*tt87*tt66-mu(1,1)&
&*tt426)
tt428 = H_invDmH(6,1)*tt64-H_invDmH(6,2)*tt63+H_invDmH(6,3)*tt62
tt429 = detDmH(1,1)*gw(1,1)*(tt2*tt428*tt21+tt2*tt66*tt90-mu(1,1)&
&*tt428)
tt430 = detDmH(1,1)*gw(1,1)*(tt2*tt66*tt96+tt404)
tt431 = (H_invDmH(4,1)*H_invDmH(6,2)-H_invDmH(4,2)*H_invDmH(6,1))&
&*tt17-(H_invDmH(4,1)*H_invDmH(6,3)-H_invDmH(4,3)*H_invDmH(6,1))*t&
&t18+(H_invDmH(4,2)*H_invDmH(6,3)-H_invDmH(4,3)*H_invDmH(6,2))*tt1&
&9
tt432 = detDmH(1,1)*gw(1,1)*(tt2*tt431*tt21+tt2*tt102*tt66-mu(1,1&
&)*tt431)
tt433 = H_invDmH(7,1)*tt64-H_invDmH(7,2)*tt63+H_invDmH(7,3)*tt62
tt434 = detDmH(1,1)*gw(1,1)*(tt2*tt433*tt21+tt2*tt66*tt105-mu(1,1&
&)*tt433)
tt435 = detDmH(1,1)*gw(1,1)*(tt2*tt66*tt111+tt410)
tt436 = (H_invDmH(4,1)*H_invDmH(7,2)-H_invDmH(4,2)*H_invDmH(7,1))&
&*tt17-(H_invDmH(4,1)*H_invDmH(7,3)-H_invDmH(4,3)*H_invDmH(7,1))*t&
&t18+(H_invDmH(4,2)*H_invDmH(7,3)-H_invDmH(4,3)*H_invDmH(7,2))*tt1&
&9
tt437 = detDmH(1,1)*gw(1,1)*(tt2*tt436*tt21+tt2*tt117*tt66-mu(1,1&
&)*tt436)
tt438 = H_invDmH(8,1)*tt64-H_invDmH(8,2)*tt63+tt62*H_invDmH(8,3)
tt439 = detDmH(1,1)*gw(1,1)*(tt2*tt438*tt21+tt2*tt66*tt120-mu(1,1&
&)*tt438)
tt440 = detDmH(1,1)*gw(1,1)*(tt2*tt66*tt126+tt416)
tt441 = tt19*(H_invDmH(4,2)*H_invDmH(8,3)-H_invDmH(4,3)*H_invDmH(&
&8,2))-tt18*(H_invDmH(4,1)*H_invDmH(8,3)-H_invDmH(4,3)*H_invDmH(8,&
&1))+(H_invDmH(4,1)*H_invDmH(8,2)-H_invDmH(4,2)*H_invDmH(8,1))*tt1&
&7
tt442 = detDmH(1,1)*gw(1,1)*(tt2*tt441*tt21+tt2*tt132*tt66-mu(1,1&
&)*tt441)
tt443 = H_invDmH(5,1)*tt70-H_invDmH(5,2)*tt69+H_invDmH(5,3)*tt68
tt444 = detDmH(1,1)*gw(1,1)*(tt2*tt443*tt21+tt2*tt72*tt75-mu(1,1)&
&*tt443)
tt445 = (H_invDmH(4,2)*H_invDmH(5,1)-H_invDmH(4,1)*H_invDmH(5,2))&
&*tt17-(H_invDmH(4,3)*H_invDmH(5,1)-H_invDmH(4,1)*H_invDmH(5,3))*t&
&t18+(H_invDmH(4,3)*H_invDmH(5,2)-H_invDmH(4,2)*H_invDmH(5,3))*tt1&
&9
tt446 = detDmH(1,1)*gw(1,1)*(tt2*tt445*tt21+tt2*tt72*tt81-mu(1,1)&
&*tt445)
tt447 = detDmH(1,1)*gw(1,1)*(tt2*tt72*tt87+tt398)
tt448 = H_invDmH(6,1)*tt70-H_invDmH(6,2)*tt69+H_invDmH(6,3)*tt68
tt449 = detDmH(1,1)*gw(1,1)*(tt2*tt448*tt21+tt2*tt72*tt90-mu(1,1)&
&*tt448)
tt450 = (H_invDmH(4,2)*H_invDmH(6,1)-H_invDmH(4,1)*H_invDmH(6,2))&
&*tt17-(H_invDmH(4,3)*H_invDmH(6,1)-H_invDmH(4,1)*H_invDmH(6,3))*t&
&t18+(H_invDmH(4,3)*H_invDmH(6,2)-H_invDmH(4,2)*H_invDmH(6,3))*tt1&
&9
tt451 = detDmH(1,1)*gw(1,1)*(tt2*tt450*tt21+tt2*tt72*tt96-mu(1,1)&
&*tt450)
tt452 = detDmH(1,1)*gw(1,1)*(tt2*tt72*tt102+tt404)
tt453 = H_invDmH(7,1)*tt70-H_invDmH(7,2)*tt69+H_invDmH(7,3)*tt68
tt454 = detDmH(1,1)*gw(1,1)*(tt2*tt453*tt21+tt2*tt72*tt105-mu(1,1&
&)*tt453)
tt455 = (H_invDmH(4,2)*H_invDmH(7,1)-H_invDmH(4,1)*H_invDmH(7,2))&
&*tt17-(H_invDmH(4,3)*H_invDmH(7,1)-H_invDmH(4,1)*H_invDmH(7,3))*t&
&t18+(H_invDmH(4,3)*H_invDmH(7,2)-H_invDmH(4,2)*H_invDmH(7,3))*tt1&
&9
tt456 = detDmH(1,1)*gw(1,1)*(tt2*tt455*tt21+tt2*tt72*tt111-mu(1,1&
&)*tt455)
tt457 = detDmH(1,1)*gw(1,1)*(tt2*tt72*tt117+tt410)
tt458 = H_invDmH(8,1)*tt70-H_invDmH(8,2)*tt69+tt68*H_invDmH(8,3)
tt459 = detDmH(1,1)*gw(1,1)*(tt2*tt458*tt21+tt2*tt72*tt120-mu(1,1&
&)*tt458)
tt460 = tt19*(H_invDmH(4,3)*H_invDmH(8,2)-H_invDmH(4,2)*H_invDmH(&
&8,3))-tt18*(H_invDmH(4,3)*H_invDmH(8,1)-H_invDmH(4,1)*H_invDmH(8,&
&3))+(H_invDmH(4,2)*H_invDmH(8,1)-H_invDmH(4,1)*H_invDmH(8,2))*tt1&
&7
tt461 = detDmH(1,1)*gw(1,1)*(tt2*tt460*tt21+tt2*tt72*tt126-mu(1,1&
&)*tt460)
tt462 = detDmH(1,1)*gw(1,1)*(tt2*tt72*tt132+tt416)
tt463 = (mu(1,1)*(2*H_invDmH(5,3)**2+2*H_invDmH(5,2)**2+2*H_invDm&
&H(5,1)**2))/2.0E+0
tt464 = H_invDmH(5,1)*tt79-H_invDmH(5,2)*tt78+H_invDmH(5,3)*tt77
tt465 = detDmH(1,1)*gw(1,1)*(tt2*tt464*tt21+tt2*tt81*tt75-mu(1,1)&
&*tt464)
tt466 = H_invDmH(5,1)*tt85-H_invDmH(5,2)*tt84+H_invDmH(5,3)*tt83
tt467 = detDmH(1,1)*gw(1,1)*(tt2*tt466*tt21+tt2*tt87*tt75-mu(1,1)&
&*tt466)
tt468 = (mu(1,1)*(2*H_invDmH(5,3)*H_invDmH(6,3)+2*H_invDmH(5,2)*H&
&_invDmH(6,2)+2*H_invDmH(5,1)*H_invDmH(6,1)))/2.0E+0
tt469 = detDmH(1,1)*gw(1,1)*(tt2*tt75*tt90+tt468)
tt470 = H_invDmH(5,1)*tt94-H_invDmH(5,2)*tt93+H_invDmH(5,3)*tt92
tt471 = detDmH(1,1)*gw(1,1)*(tt2*tt470*tt21+tt2*tt96*tt75-mu(1,1)&
&*tt470)
tt472 = H_invDmH(5,1)*tt100-H_invDmH(5,2)*tt99+H_invDmH(5,3)*tt98&
&
tt473 = detDmH(1,1)*gw(1,1)*(tt2*tt472*tt21+tt2*tt102*tt75-mu(1,1&
&)*tt472)
tt474 = (mu(1,1)*(2*H_invDmH(5,3)*H_invDmH(7,3)+2*H_invDmH(5,2)*H&
&_invDmH(7,2)+2*H_invDmH(5,1)*H_invDmH(7,1)))/2.0E+0
tt475 = detDmH(1,1)*gw(1,1)*(tt2*tt75*tt105+tt474)
tt476 = H_invDmH(5,1)*tt109-H_invDmH(5,2)*tt108+H_invDmH(5,3)*tt1&
&07
tt477 = detDmH(1,1)*gw(1,1)*(tt2*tt476*tt21+tt2*tt111*tt75-mu(1,1&
&)*tt476)
tt478 = H_invDmH(5,1)*tt115-H_invDmH(5,2)*tt114+H_invDmH(5,3)*tt1&
&13
tt479 = detDmH(1,1)*gw(1,1)*(tt2*tt478*tt21+tt2*tt117*tt75-mu(1,1&
&)*tt478)
tt480 = (mu(1,1)*(2*H_invDmH(5,3)*H_invDmH(8,3)+2*H_invDmH(5,2)*H&
&_invDmH(8,2)+2*H_invDmH(5,1)*H_invDmH(8,1)))/2.0E+0
tt481 = detDmH(1,1)*gw(1,1)*(tt2*tt75*tt120+tt480)
tt482 = H_invDmH(5,1)*tt124-H_invDmH(5,2)*tt123+H_invDmH(5,3)*tt1&
&22
tt483 = detDmH(1,1)*gw(1,1)*(tt2*tt482*tt21+tt2*tt126*tt75-mu(1,1&
&)*tt482)
tt484 = H_invDmH(5,1)*tt130-H_invDmH(5,2)*tt129+H_invDmH(5,3)*tt1&
&28
tt485 = detDmH(1,1)*gw(1,1)*(tt2*tt484*tt21+tt2*tt132*tt75-mu(1,1&
&)*tt484)
tt486 = detDmH(1,1)*gw(1,1)*tt2*tt87*tt81
tt487 = H_invDmH(6,1)*tt79-H_invDmH(6,2)*tt78+H_invDmH(6,3)*tt77
tt488 = detDmH(1,1)*gw(1,1)*(tt2*tt487*tt21+tt2*tt81*tt90-mu(1,1)&
&*tt487)
tt489 = detDmH(1,1)*gw(1,1)*(tt2*tt81*tt96+tt468)
tt490 = (H_invDmH(5,1)*H_invDmH(6,2)-H_invDmH(5,2)*H_invDmH(6,1))&
&*tt17-(H_invDmH(5,1)*H_invDmH(6,3)-H_invDmH(5,3)*H_invDmH(6,1))*t&
&t18+(H_invDmH(5,2)*H_invDmH(6,3)-H_invDmH(5,3)*H_invDmH(6,2))*tt1&
&9
tt491 = detDmH(1,1)*gw(1,1)*(tt2*tt490*tt21+tt2*tt102*tt81-mu(1,1&
&)*tt490)
tt492 = H_invDmH(7,1)*tt79-H_invDmH(7,2)*tt78+H_invDmH(7,3)*tt77
tt493 = detDmH(1,1)*gw(1,1)*(tt2*tt492*tt21+tt2*tt81*tt105-mu(1,1&
&)*tt492)
tt494 = detDmH(1,1)*gw(1,1)*(tt2*tt81*tt111+tt474)
tt495 = (H_invDmH(5,1)*H_invDmH(7,2)-H_invDmH(5,2)*H_invDmH(7,1))&
&*tt17-(H_invDmH(5,1)*H_invDmH(7,3)-H_invDmH(5,3)*H_invDmH(7,1))*t&
&t18+(H_invDmH(5,2)*H_invDmH(7,3)-H_invDmH(5,3)*H_invDmH(7,2))*tt1&
&9
tt496 = detDmH(1,1)*gw(1,1)*(tt2*tt495*tt21+tt2*tt117*tt81-mu(1,1&
&)*tt495)
tt497 = H_invDmH(8,1)*tt79-H_invDmH(8,2)*tt78+tt77*H_invDmH(8,3)
tt498 = detDmH(1,1)*gw(1,1)*(tt2*tt497*tt21+tt2*tt81*tt120-mu(1,1&
&)*tt497)
tt499 = detDmH(1,1)*gw(1,1)*(tt2*tt81*tt126+tt480)
tt500 = tt19*(H_invDmH(5,2)*H_invDmH(8,3)-H_invDmH(5,3)*H_invDmH(&
&8,2))-tt18*(H_invDmH(5,1)*H_invDmH(8,3)-H_invDmH(5,3)*H_invDmH(8,&
&1))+(H_invDmH(5,1)*H_invDmH(8,2)-H_invDmH(5,2)*H_invDmH(8,1))*tt1&
&7
tt501 = detDmH(1,1)*gw(1,1)*(tt2*tt500*tt21+tt2*tt132*tt81-mu(1,1&
&)*tt500)
tt502 = H_invDmH(6,1)*tt85-H_invDmH(6,2)*tt84+H_invDmH(6,3)*tt83
tt503 = detDmH(1,1)*gw(1,1)*(tt2*tt502*tt21+tt2*tt87*tt90-mu(1,1)&
&*tt502)
tt504 = (H_invDmH(5,2)*H_invDmH(6,1)-H_invDmH(5,1)*H_invDmH(6,2))&
&*tt17-(H_invDmH(5,3)*H_invDmH(6,1)-H_invDmH(5,1)*H_invDmH(6,3))*t&
&t18+(H_invDmH(5,3)*H_invDmH(6,2)-H_invDmH(5,2)*H_invDmH(6,3))*tt1&
&9
tt505 = detDmH(1,1)*gw(1,1)*(tt2*tt504*tt21+tt2*tt87*tt96-mu(1,1)&
&*tt504)
tt506 = detDmH(1,1)*gw(1,1)*(tt2*tt87*tt102+tt468)
tt507 = H_invDmH(7,1)*tt85-H_invDmH(7,2)*tt84+H_invDmH(7,3)*tt83
tt508 = detDmH(1,1)*gw(1,1)*(tt2*tt507*tt21+tt2*tt87*tt105-mu(1,1&
&)*tt507)
tt509 = (H_invDmH(5,2)*H_invDmH(7,1)-H_invDmH(5,1)*H_invDmH(7,2))&
&*tt17-(H_invDmH(5,3)*H_invDmH(7,1)-H_invDmH(5,1)*H_invDmH(7,3))*t&
&t18+(H_invDmH(5,3)*H_invDmH(7,2)-H_invDmH(5,2)*H_invDmH(7,3))*tt1&
&9
tt510 = detDmH(1,1)*gw(1,1)*(tt2*tt509*tt21+tt2*tt87*tt111-mu(1,1&
&)*tt509)
tt511 = detDmH(1,1)*gw(1,1)*(tt2*tt87*tt117+tt474)
tt512 = H_invDmH(8,1)*tt85-H_invDmH(8,2)*tt84+tt83*H_invDmH(8,3)
tt513 = detDmH(1,1)*gw(1,1)*(tt2*tt512*tt21+tt2*tt87*tt120-mu(1,1&
&)*tt512)
tt514 = tt19*(H_invDmH(5,3)*H_invDmH(8,2)-H_invDmH(5,2)*H_invDmH(&
&8,3))-tt18*(H_invDmH(5,3)*H_invDmH(8,1)-H_invDmH(5,1)*H_invDmH(8,&
&3))+(H_invDmH(5,2)*H_invDmH(8,1)-H_invDmH(5,1)*H_invDmH(8,2))*tt1&
&7
tt515 = detDmH(1,1)*gw(1,1)*(tt2*tt514*tt21+tt2*tt87*tt126-mu(1,1&
&)*tt514)
tt516 = detDmH(1,1)*gw(1,1)*(tt2*tt87*tt132+tt480)
tt517 = (mu(1,1)*(2*H_invDmH(6,3)**2+2*H_invDmH(6,2)**2+2*H_invDm&
&H(6,1)**2))/2.0E+0
tt518 = H_invDmH(6,1)*tt94-H_invDmH(6,2)*tt93+H_invDmH(6,3)*tt92
tt519 = detDmH(1,1)*gw(1,1)*(tt2*tt518*tt21+tt2*tt96*tt90-mu(1,1)&
&*tt518)
tt520 = H_invDmH(6,1)*tt100-H_invDmH(6,2)*tt99+H_invDmH(6,3)*tt98&
&
tt521 = detDmH(1,1)*gw(1,1)*(tt2*tt520*tt21+tt2*tt102*tt90-mu(1,1&
&)*tt520)
tt522 = (mu(1,1)*(2*H_invDmH(6,3)*H_invDmH(7,3)+2*H_invDmH(6,2)*H&
&_invDmH(7,2)+2*H_invDmH(6,1)*H_invDmH(7,1)))/2.0E+0
tt523 = detDmH(1,1)*gw(1,1)*(tt2*tt90*tt105+tt522)
tt524 = H_invDmH(6,1)*tt109-H_invDmH(6,2)*tt108+H_invDmH(6,3)*tt1&
&07
tt525 = detDmH(1,1)*gw(1,1)*(tt2*tt524*tt21+tt2*tt111*tt90-mu(1,1&
&)*tt524)
tt526 = H_invDmH(6,1)*tt115-H_invDmH(6,2)*tt114+H_invDmH(6,3)*tt1&
&13
tt527 = detDmH(1,1)*gw(1,1)*(tt2*tt526*tt21+tt2*tt117*tt90-mu(1,1&
&)*tt526)
tt528 = (mu(1,1)*(2*H_invDmH(6,3)*H_invDmH(8,3)+2*H_invDmH(6,2)*H&
&_invDmH(8,2)+2*H_invDmH(6,1)*H_invDmH(8,1)))/2.0E+0
tt529 = detDmH(1,1)*gw(1,1)*(tt2*tt90*tt120+tt528)
tt530 = H_invDmH(6,1)*tt124-H_invDmH(6,2)*tt123+H_invDmH(6,3)*tt1&
&22
tt531 = detDmH(1,1)*gw(1,1)*(tt2*tt530*tt21+tt2*tt126*tt90-mu(1,1&
&)*tt530)
tt532 = H_invDmH(6,1)*tt130-H_invDmH(6,2)*tt129+H_invDmH(6,3)*tt1&
&28
tt533 = detDmH(1,1)*gw(1,1)*(tt2*tt532*tt21+tt2*tt132*tt90-mu(1,1&
&)*tt532)
tt534 = detDmH(1,1)*gw(1,1)*tt2*tt102*tt96
tt535 = H_invDmH(7,1)*tt94-H_invDmH(7,2)*tt93+H_invDmH(7,3)*tt92
tt536 = detDmH(1,1)*gw(1,1)*(tt2*tt535*tt21+tt2*tt96*tt105-mu(1,1&
&)*tt535)
tt537 = detDmH(1,1)*gw(1,1)*(tt2*tt96*tt111+tt522)
tt538 = (H_invDmH(6,1)*H_invDmH(7,2)-H_invDmH(6,2)*H_invDmH(7,1))&
&*tt17-(H_invDmH(6,1)*H_invDmH(7,3)-H_invDmH(6,3)*H_invDmH(7,1))*t&
&t18+(H_invDmH(6,2)*H_invDmH(7,3)-H_invDmH(6,3)*H_invDmH(7,2))*tt1&
&9
tt539 = detDmH(1,1)*gw(1,1)*(tt2*tt538*tt21+tt2*tt117*tt96-mu(1,1&
&)*tt538)
tt540 = H_invDmH(8,1)*tt94-H_invDmH(8,2)*tt93+tt92*H_invDmH(8,3)
tt541 = detDmH(1,1)*gw(1,1)*(tt2*tt540*tt21+tt2*tt96*tt120-mu(1,1&
&)*tt540)
tt542 = detDmH(1,1)*gw(1,1)*(tt2*tt96*tt126+tt528)
tt543 = tt19*(H_invDmH(6,2)*H_invDmH(8,3)-H_invDmH(6,3)*H_invDmH(&
&8,2))-tt18*(H_invDmH(6,1)*H_invDmH(8,3)-H_invDmH(6,3)*H_invDmH(8,&
&1))+(H_invDmH(6,1)*H_invDmH(8,2)-H_invDmH(6,2)*H_invDmH(8,1))*tt1&
&7
tt544 = detDmH(1,1)*gw(1,1)*(tt2*tt543*tt21+tt2*tt132*tt96-mu(1,1&
&)*tt543)
tt545 = H_invDmH(7,1)*tt100-H_invDmH(7,2)*tt99+H_invDmH(7,3)*tt98&
&
tt546 = detDmH(1,1)*gw(1,1)*(tt2*tt545*tt21+tt2*tt102*tt105-mu(1,&
&1)*tt545)
tt547 = (H_invDmH(6,2)*H_invDmH(7,1)-H_invDmH(6,1)*H_invDmH(7,2))&
&*tt17-(H_invDmH(6,3)*H_invDmH(7,1)-H_invDmH(6,1)*H_invDmH(7,3))*t&
&t18+(H_invDmH(6,3)*H_invDmH(7,2)-H_invDmH(6,2)*H_invDmH(7,3))*tt1&
&9
tt548 = detDmH(1,1)*gw(1,1)*(tt2*tt547*tt21+tt2*tt102*tt111-mu(1,&
&1)*tt547)
tt549 = detDmH(1,1)*gw(1,1)*(tt2*tt102*tt117+tt522)
tt550 = H_invDmH(8,1)*tt100-H_invDmH(8,2)*tt99+tt98*H_invDmH(8,3)&
&
tt551 = detDmH(1,1)*gw(1,1)*(tt2*tt550*tt21+tt2*tt102*tt120-mu(1,&
&1)*tt550)
tt552 = tt19*(H_invDmH(6,3)*H_invDmH(8,2)-H_invDmH(6,2)*H_invDmH(&
&8,3))-tt18*(H_invDmH(6,3)*H_invDmH(8,1)-H_invDmH(6,1)*H_invDmH(8,&
&3))+(H_invDmH(6,2)*H_invDmH(8,1)-H_invDmH(6,1)*H_invDmH(8,2))*tt1&
&7
tt553 = detDmH(1,1)*gw(1,1)*(tt2*tt552*tt21+tt2*tt102*tt126-mu(1,&
&1)*tt552)
tt554 = detDmH(1,1)*gw(1,1)*(tt2*tt102*tt132+tt528)
tt555 = (mu(1,1)*(2*H_invDmH(7,3)**2+2*H_invDmH(7,2)**2+2*H_invDm&
&H(7,1)**2))/2.0E+0
tt556 = H_invDmH(7,1)*tt109-H_invDmH(7,2)*tt108+H_invDmH(7,3)*tt1&
&07
tt557 = detDmH(1,1)*gw(1,1)*(tt2*tt556*tt21+tt2*tt111*tt105-mu(1,&
&1)*tt556)
tt558 = H_invDmH(7,1)*tt115-H_invDmH(7,2)*tt114+H_invDmH(7,3)*tt1&
&13
tt559 = detDmH(1,1)*gw(1,1)*(tt2*tt558*tt21+tt2*tt117*tt105-mu(1,&
&1)*tt558)
tt560 = (mu(1,1)*(2*H_invDmH(7,3)*H_invDmH(8,3)+2*H_invDmH(7,2)*H&
&_invDmH(8,2)+2*H_invDmH(7,1)*H_invDmH(8,1)))/2.0E+0
tt561 = detDmH(1,1)*gw(1,1)*(tt2*tt105*tt120+tt560)
tt562 = H_invDmH(7,1)*tt124-H_invDmH(7,2)*tt123+H_invDmH(7,3)*tt1&
&22
tt563 = detDmH(1,1)*gw(1,1)*(tt2*tt562*tt21+tt2*tt126*tt105-mu(1,&
&1)*tt562)
tt564 = H_invDmH(7,1)*tt130-H_invDmH(7,2)*tt129+H_invDmH(7,3)*tt1&
&28
tt565 = detDmH(1,1)*gw(1,1)*(tt2*tt564*tt21+tt2*tt132*tt105-mu(1,&
&1)*tt564)
tt566 = detDmH(1,1)*gw(1,1)*tt2*tt117*tt111
tt567 = H_invDmH(8,1)*tt109-H_invDmH(8,2)*tt108+tt107*H_invDmH(8,&
&3)
tt568 = detDmH(1,1)*gw(1,1)*(tt2*tt567*tt21+tt2*tt111*tt120-mu(1,&
&1)*tt567)
tt569 = detDmH(1,1)*gw(1,1)*(tt2*tt111*tt126+tt560)
tt570 = tt19*(H_invDmH(7,2)*H_invDmH(8,3)-H_invDmH(7,3)*H_invDmH(&
&8,2))-tt18*(H_invDmH(7,1)*H_invDmH(8,3)-H_invDmH(7,3)*H_invDmH(8,&
&1))+(H_invDmH(7,1)*H_invDmH(8,2)-H_invDmH(7,2)*H_invDmH(8,1))*tt1&
&7
tt571 = detDmH(1,1)*gw(1,1)*(tt2*tt570*tt21+tt2*tt132*tt111-mu(1,&
&1)*tt570)
tt572 = H_invDmH(8,1)*tt115-H_invDmH(8,2)*tt114+tt113*H_invDmH(8,&
&3)
tt573 = detDmH(1,1)*gw(1,1)*(tt2*tt572*tt21+tt2*tt117*tt120-mu(1,&
&1)*tt572)
tt574 = tt19*(H_invDmH(7,3)*H_invDmH(8,2)-H_invDmH(7,2)*H_invDmH(&
&8,3))-tt18*(H_invDmH(7,3)*H_invDmH(8,1)-H_invDmH(7,1)*H_invDmH(8,&
&3))+(H_invDmH(7,2)*H_invDmH(8,1)-H_invDmH(7,1)*H_invDmH(8,2))*tt1&
&7
tt575 = detDmH(1,1)*gw(1,1)*(tt2*tt574*tt21+tt2*tt117*tt126-mu(1,&
&1)*tt574)
tt576 = detDmH(1,1)*gw(1,1)*(tt2*tt117*tt132+tt560)
tt577 = (mu(1,1)*(2*H_invDmH(8,3)**2+2*H_invDmH(8,2)**2+2*H_invDm&
&H(8,1)**2))/2.0E+0
tt578 = H_invDmH(8,1)*tt124-H_invDmH(8,2)*tt123+tt122*H_invDmH(8,&
&3)
tt579 = detDmH(1,1)*gw(1,1)*(tt2*tt578*tt21+tt2*tt126*tt120-mu(1,&
&1)*tt578)
tt580 = H_invDmH(8,1)*tt130-H_invDmH(8,2)*tt129+tt128*H_invDmH(8,&
&3)
tt581 = detDmH(1,1)*gw(1,1)*(tt2*tt580*tt21+tt2*tt132*tt120-mu(1,&
&1)*tt580)
tt582 = detDmH(1,1)*gw(1,1)*tt2*tt132*tt126
hes(1,1) = detDmH(1,1)*gw(1,1)*(tt2*tt12**2+tt1)
hes(1,2) = tt22
hes(1,3) = tt28
hes(1,4) = tt31
hes(1,5) = tt37
hes(1,6) = tt43
hes(1,7) = tt46
hes(1,8) = tt52
hes(1,9) = tt58
hes(1,10) = tt61
hes(1,11) = tt67
hes(1,12) = tt73
hes(1,13) = tt76
hes(1,14) = tt82
hes(1,15) = tt88
hes(1,16) = tt91
hes(1,17) = tt97
hes(1,18) = tt103
hes(1,19) = tt106
hes(1,20) = tt112
hes(1,21) = tt118
hes(1,22) = tt121
hes(1,23) = tt127
hes(1,24) = tt133
hes(2,1) = tt22
hes(2,2) = detDmH(1,1)*gw(1,1)*(tt2*tt20**2+tt1)
hes(2,3) = tt134
hes(2,4) = tt136
hes(2,5) = tt137
hes(2,6) = tt139
hes(2,7) = tt141
hes(2,8) = tt142
hes(2,9) = tt144
hes(2,10) = tt146
hes(2,11) = tt147
hes(2,12) = tt149
hes(2,13) = tt151
hes(2,14) = tt152
hes(2,15) = tt154
hes(2,16) = tt156
hes(2,17) = tt157
hes(2,18) = tt159
hes(2,19) = tt161
hes(2,20) = tt162
hes(2,21) = tt164
hes(2,22) = tt166
hes(2,23) = tt167
hes(2,24) = tt169
hes(3,1) = tt28
hes(3,2) = tt134
hes(3,3) = detDmH(1,1)*gw(1,1)*(tt2*tt27**2+tt1)
hes(3,4) = tt171
hes(3,5) = tt173
hes(3,6) = tt174
hes(3,7) = tt176
hes(3,8) = tt178
hes(3,9) = tt179
hes(3,10) = tt181
hes(3,11) = tt183
hes(3,12) = tt184
hes(3,13) = tt186
hes(3,14) = tt188
hes(3,15) = tt189
hes(3,16) = tt191
hes(3,17) = tt193
hes(3,18) = tt194
hes(3,19) = tt196
hes(3,20) = tt198
hes(3,21) = tt199
hes(3,22) = tt201
hes(3,23) = tt203
hes(3,24) = tt204
hes(4,1) = tt31
hes(4,2) = tt136
hes(4,3) = tt171
hes(4,4) = detDmH(1,1)*gw(1,1)*(tt2*tt30**2+tt205)
hes(4,5) = tt207
hes(4,6) = tt209
hes(4,7) = tt211
hes(4,8) = tt213
hes(4,9) = tt215
hes(4,10) = tt217
hes(4,11) = tt219
hes(4,12) = tt221
hes(4,13) = tt223
hes(4,14) = tt225
hes(4,15) = tt227
hes(4,16) = tt229
hes(4,17) = tt231
hes(4,18) = tt233
hes(4,19) = tt235
hes(4,20) = tt237
hes(4,21) = tt239
hes(4,22) = tt241
hes(4,23) = tt243
hes(4,24) = tt245
hes(5,1) = tt37
hes(5,2) = tt137
hes(5,3) = tt173
hes(5,4) = tt207
hes(5,5) = detDmH(1,1)*gw(1,1)*(tt2*tt36**2+tt205)
hes(5,6) = tt246
hes(5,7) = tt248
hes(5,8) = tt249
hes(5,9) = tt251
hes(5,10) = tt253
hes(5,11) = tt254
hes(5,12) = tt256
hes(5,13) = tt258
hes(5,14) = tt259
hes(5,15) = tt261
hes(5,16) = tt263
hes(5,17) = tt264
hes(5,18) = tt266
hes(5,19) = tt268
hes(5,20) = tt269
hes(5,21) = tt271
hes(5,22) = tt273
hes(5,23) = tt274
hes(5,24) = tt276
hes(6,1) = tt43
hes(6,2) = tt139
hes(6,3) = tt174
hes(6,4) = tt209
hes(6,5) = tt246
hes(6,6) = detDmH(1,1)*gw(1,1)*(tt2*tt42**2+tt205)
hes(6,7) = tt278
hes(6,8) = tt280
hes(6,9) = tt281
hes(6,10) = tt283
hes(6,11) = tt285
hes(6,12) = tt286
hes(6,13) = tt288
hes(6,14) = tt290
hes(6,15) = tt291
hes(6,16) = tt293
hes(6,17) = tt295
hes(6,18) = tt296
hes(6,19) = tt298
hes(6,20) = tt300
hes(6,21) = tt301
hes(6,22) = tt303
hes(6,23) = tt305
hes(6,24) = tt306
hes(7,1) = tt46
hes(7,2) = tt141
hes(7,3) = tt176
hes(7,4) = tt211
hes(7,5) = tt248
hes(7,6) = tt278
hes(7,7) = detDmH(1,1)*gw(1,1)*(tt2*tt45**2+tt307)
hes(7,8) = tt309
hes(7,9) = tt311
hes(7,10) = tt313
hes(7,11) = tt315
hes(7,12) = tt317
hes(7,13) = tt319
hes(7,14) = tt321
hes(7,15) = tt323
hes(7,16) = tt325
hes(7,17) = tt327
hes(7,18) = tt329
hes(7,19) = tt331
hes(7,20) = tt333
hes(7,21) = tt335
hes(7,22) = tt337
hes(7,23) = tt339
hes(7,24) = tt341
hes(8,1) = tt52
hes(8,2) = tt142
hes(8,3) = tt178
hes(8,4) = tt213
hes(8,5) = tt249
hes(8,6) = tt280
hes(8,7) = tt309
hes(8,8) = detDmH(1,1)*gw(1,1)*(tt2*tt51**2+tt307)
hes(8,9) = tt342
hes(8,10) = tt344
hes(8,11) = tt345
hes(8,12) = tt347
hes(8,13) = tt349
hes(8,14) = tt350
hes(8,15) = tt352
hes(8,16) = tt354
hes(8,17) = tt355
hes(8,18) = tt357
hes(8,19) = tt359
hes(8,20) = tt360
hes(8,21) = tt362
hes(8,22) = tt364
hes(8,23) = tt365
hes(8,24) = tt367
hes(9,1) = tt58
hes(9,2) = tt144
hes(9,3) = tt179
hes(9,4) = tt215
hes(9,5) = tt251
hes(9,6) = tt281
hes(9,7) = tt311
hes(9,8) = tt342
hes(9,9) = detDmH(1,1)*gw(1,1)*(tt2*tt57**2+tt307)
hes(9,10) = tt369
hes(9,11) = tt371
hes(9,12) = tt372
hes(9,13) = tt374
hes(9,14) = tt376
hes(9,15) = tt377
hes(9,16) = tt379
hes(9,17) = tt381
hes(9,18) = tt382
hes(9,19) = tt384
hes(9,20) = tt386
hes(9,21) = tt387
hes(9,22) = tt389
hes(9,23) = tt391
hes(9,24) = tt392
hes(10,1) = tt61
hes(10,2) = tt146
hes(10,3) = tt181
hes(10,4) = tt217
hes(10,5) = tt253
hes(10,6) = tt283
hes(10,7) = tt313
hes(10,8) = tt344
hes(10,9) = tt369
hes(10,10) = detDmH(1,1)*gw(1,1)*(tt2*tt60**2+tt393)
hes(10,11) = tt395
hes(10,12) = tt397
hes(10,13) = tt399
hes(10,14) = tt401
hes(10,15) = tt403
hes(10,16) = tt405
hes(10,17) = tt407
hes(10,18) = tt409
hes(10,19) = tt411
hes(10,20) = tt413
hes(10,21) = tt415
hes(10,22) = tt417
hes(10,23) = tt419
hes(10,24) = tt421
hes(11,1) = tt67
hes(11,2) = tt147
hes(11,3) = tt183
hes(11,4) = tt219
hes(11,5) = tt254
hes(11,6) = tt285
hes(11,7) = tt315
hes(11,8) = tt345
hes(11,9) = tt371
hes(11,10) = tt395
hes(11,11) = detDmH(1,1)*gw(1,1)*(tt2*tt66**2+tt393)
hes(11,12) = tt422
hes(11,13) = tt424
hes(11,14) = tt425
hes(11,15) = tt427
hes(11,16) = tt429
hes(11,17) = tt430
hes(11,18) = tt432
hes(11,19) = tt434
hes(11,20) = tt435
hes(11,21) = tt437
hes(11,22) = tt439
hes(11,23) = tt440
hes(11,24) = tt442
hes(12,1) = tt73
hes(12,2) = tt149
hes(12,3) = tt184
hes(12,4) = tt221
hes(12,5) = tt256
hes(12,6) = tt286
hes(12,7) = tt317
hes(12,8) = tt347
hes(12,9) = tt372
hes(12,10) = tt397
hes(12,11) = tt422
hes(12,12) = detDmH(1,1)*gw(1,1)*(tt2*tt72**2+tt393)
hes(12,13) = tt444
hes(12,14) = tt446
hes(12,15) = tt447
hes(12,16) = tt449
hes(12,17) = tt451
hes(12,18) = tt452
hes(12,19) = tt454
hes(12,20) = tt456
hes(12,21) = tt457
hes(12,22) = tt459
hes(12,23) = tt461
hes(12,24) = tt462
hes(13,1) = tt76
hes(13,2) = tt151
hes(13,3) = tt186
hes(13,4) = tt223
hes(13,5) = tt258
hes(13,6) = tt288
hes(13,7) = tt319
hes(13,8) = tt349
hes(13,9) = tt374
hes(13,10) = tt399
hes(13,11) = tt424
hes(13,12) = tt444
hes(13,13) = detDmH(1,1)*gw(1,1)*(tt2*tt75**2+tt463)
hes(13,14) = tt465
hes(13,15) = tt467
hes(13,16) = tt469
hes(13,17) = tt471
hes(13,18) = tt473
hes(13,19) = tt475
hes(13,20) = tt477
hes(13,21) = tt479
hes(13,22) = tt481
hes(13,23) = tt483
hes(13,24) = tt485
hes(14,1) = tt82
hes(14,2) = tt152
hes(14,3) = tt188
hes(14,4) = tt225
hes(14,5) = tt259
hes(14,6) = tt290
hes(14,7) = tt321
hes(14,8) = tt350
hes(14,9) = tt376
hes(14,10) = tt401
hes(14,11) = tt425
hes(14,12) = tt446
hes(14,13) = tt465
hes(14,14) = detDmH(1,1)*gw(1,1)*(tt2*tt81**2+tt463)
hes(14,15) = tt486
hes(14,16) = tt488
hes(14,17) = tt489
hes(14,18) = tt491
hes(14,19) = tt493
hes(14,20) = tt494
hes(14,21) = tt496
hes(14,22) = tt498
hes(14,23) = tt499
hes(14,24) = tt501
hes(15,1) = tt88
hes(15,2) = tt154
hes(15,3) = tt189
hes(15,4) = tt227
hes(15,5) = tt261
hes(15,6) = tt291
hes(15,7) = tt323
hes(15,8) = tt352
hes(15,9) = tt377
hes(15,10) = tt403
hes(15,11) = tt427
hes(15,12) = tt447
hes(15,13) = tt467
hes(15,14) = tt486
hes(15,15) = detDmH(1,1)*gw(1,1)*(tt2*tt87**2+tt463)
hes(15,16) = tt503
hes(15,17) = tt505
hes(15,18) = tt506
hes(15,19) = tt508
hes(15,20) = tt510
hes(15,21) = tt511
hes(15,22) = tt513
hes(15,23) = tt515
hes(15,24) = tt516
hes(16,1) = tt91
hes(16,2) = tt156
hes(16,3) = tt191
hes(16,4) = tt229
hes(16,5) = tt263
hes(16,6) = tt293
hes(16,7) = tt325
hes(16,8) = tt354
hes(16,9) = tt379
hes(16,10) = tt405
hes(16,11) = tt429
hes(16,12) = tt449
hes(16,13) = tt469
hes(16,14) = tt488
hes(16,15) = tt503
hes(16,16) = detDmH(1,1)*gw(1,1)*(tt2*tt90**2+tt517)
hes(16,17) = tt519
hes(16,18) = tt521
hes(16,19) = tt523
hes(16,20) = tt525
hes(16,21) = tt527
hes(16,22) = tt529
hes(16,23) = tt531
hes(16,24) = tt533
hes(17,1) = tt97
hes(17,2) = tt157
hes(17,3) = tt193
hes(17,4) = tt231
hes(17,5) = tt264
hes(17,6) = tt295
hes(17,7) = tt327
hes(17,8) = tt355
hes(17,9) = tt381
hes(17,10) = tt407
hes(17,11) = tt430
hes(17,12) = tt451
hes(17,13) = tt471
hes(17,14) = tt489
hes(17,15) = tt505
hes(17,16) = tt519
hes(17,17) = detDmH(1,1)*gw(1,1)*(tt2*tt96**2+tt517)
hes(17,18) = tt534
hes(17,19) = tt536
hes(17,20) = tt537
hes(17,21) = tt539
hes(17,22) = tt541
hes(17,23) = tt542
hes(17,24) = tt544
hes(18,1) = tt103
hes(18,2) = tt159
hes(18,3) = tt194
hes(18,4) = tt233
hes(18,5) = tt266
hes(18,6) = tt296
hes(18,7) = tt329
hes(18,8) = tt357
hes(18,9) = tt382
hes(18,10) = tt409
hes(18,11) = tt432
hes(18,12) = tt452
hes(18,13) = tt473
hes(18,14) = tt491
hes(18,15) = tt506
hes(18,16) = tt521
hes(18,17) = tt534
hes(18,18) = detDmH(1,1)*gw(1,1)*(tt2*tt102**2+tt517)
hes(18,19) = tt546
hes(18,20) = tt548
hes(18,21) = tt549
hes(18,22) = tt551
hes(18,23) = tt553
hes(18,24) = tt554
hes(19,1) = tt106
hes(19,2) = tt161
hes(19,3) = tt196
hes(19,4) = tt235
hes(19,5) = tt268
hes(19,6) = tt298
hes(19,7) = tt331
hes(19,8) = tt359
hes(19,9) = tt384
hes(19,10) = tt411
hes(19,11) = tt434
hes(19,12) = tt454
hes(19,13) = tt475
hes(19,14) = tt493
hes(19,15) = tt508
hes(19,16) = tt523
hes(19,17) = tt536
hes(19,18) = tt546
hes(19,19) = detDmH(1,1)*gw(1,1)*(tt2*tt105**2+tt555)
hes(19,20) = tt557
hes(19,21) = tt559
hes(19,22) = tt561
hes(19,23) = tt563
hes(19,24) = tt565
hes(20,1) = tt112
hes(20,2) = tt162
hes(20,3) = tt198
hes(20,4) = tt237
hes(20,5) = tt269
hes(20,6) = tt300
hes(20,7) = tt333
hes(20,8) = tt360
hes(20,9) = tt386
hes(20,10) = tt413
hes(20,11) = tt435
hes(20,12) = tt456
hes(20,13) = tt477
hes(20,14) = tt494
hes(20,15) = tt510
hes(20,16) = tt525
hes(20,17) = tt537
hes(20,18) = tt548
hes(20,19) = tt557
hes(20,20) = detDmH(1,1)*gw(1,1)*(tt2*tt111**2+tt555)
hes(20,21) = tt566
hes(20,22) = tt568
hes(20,23) = tt569
hes(20,24) = tt571
hes(21,1) = tt118
hes(21,2) = tt164
hes(21,3) = tt199
hes(21,4) = tt239
hes(21,5) = tt271
hes(21,6) = tt301
hes(21,7) = tt335
hes(21,8) = tt362
hes(21,9) = tt387
hes(21,10) = tt415
hes(21,11) = tt437
hes(21,12) = tt457
hes(21,13) = tt479
hes(21,14) = tt496
hes(21,15) = tt511
hes(21,16) = tt527
hes(21,17) = tt539
hes(21,18) = tt549
hes(21,19) = tt559
hes(21,20) = tt566
hes(21,21) = detDmH(1,1)*gw(1,1)*(tt2*tt117**2+tt555)
hes(21,22) = tt573
hes(21,23) = tt575
hes(21,24) = tt576
hes(22,1) = tt121
hes(22,2) = tt166
hes(22,3) = tt201
hes(22,4) = tt241
hes(22,5) = tt273
hes(22,6) = tt303
hes(22,7) = tt337
hes(22,8) = tt364
hes(22,9) = tt389
hes(22,10) = tt417
hes(22,11) = tt439
hes(22,12) = tt459
hes(22,13) = tt481
hes(22,14) = tt498
hes(22,15) = tt513
hes(22,16) = tt529
hes(22,17) = tt541
hes(22,18) = tt551
hes(22,19) = tt561
hes(22,20) = tt568
hes(22,21) = tt573
hes(22,22) = detDmH(1,1)*gw(1,1)*(tt2*tt120**2+tt577)
hes(22,23) = tt579
hes(22,24) = tt581
hes(23,1) = tt127
hes(23,2) = tt167
hes(23,3) = tt203
hes(23,4) = tt243
hes(23,5) = tt274
hes(23,6) = tt305
hes(23,7) = tt339
hes(23,8) = tt365
hes(23,9) = tt391
hes(23,10) = tt419
hes(23,11) = tt440
hes(23,12) = tt461
hes(23,13) = tt483
hes(23,14) = tt499
hes(23,15) = tt515
hes(23,16) = tt531
hes(23,17) = tt542
hes(23,18) = tt553
hes(23,19) = tt563
hes(23,20) = tt569
hes(23,21) = tt575
hes(23,22) = tt579
hes(23,23) = detDmH(1,1)*gw(1,1)*(tt2*tt126**2+tt577)
hes(23,24) = tt582
hes(24,1) = tt133
hes(24,2) = tt169
hes(24,3) = tt204
hes(24,4) = tt245
hes(24,5) = tt276
hes(24,6) = tt306
hes(24,7) = tt341
hes(24,8) = tt367
hes(24,9) = tt392
hes(24,10) = tt421
hes(24,11) = tt442
hes(24,12) = tt462
hes(24,13) = tt485
hes(24,14) = tt501
hes(24,15) = tt516
hes(24,16) = tt533
hes(24,17) = tt544
hes(24,18) = tt554
hes(24,19) = tt565
hes(24,20) = tt571
hes(24,21) = tt576
hes(24,22) = tt581
hes(24,23) = tt582
hes(24,24) = detDmH(1,1)*gw(1,1)*(tt2*tt132**2+tt577)
END 
SUBROUTINE vox_bower_neo_at_quadr(val, X, H_invDmH, detDmH, gw, mu, lam) 
IMPLICIT NONE 
REAL(KIND=8) val(1, 1) 
REAL(KIND=8) X(3, 8) 
REAL(KIND=8) H_invDmH(8, 3) 
REAL(KIND=8) detDmH(1, 1) 
REAL(KIND=8) gw(1, 1) 
REAL(KIND=8) mu(1, 1) 
REAL(KIND=8) lam(1, 1) 
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
tt1 = X(3,8)*H_invDmH(8,1)+X(3,7)*H_invDmH(7,1)+X(3,6)*H_invDmH(6&
&,1)+X(3,5)*H_invDmH(5,1)+X(3,4)*H_invDmH(4,1)+H_invDmH(3,1)*X(3,3&
&)+H_invDmH(2,1)*X(3,2)+H_invDmH(1,1)*X(3,1)
tt2 = X(2,8)*H_invDmH(8,2)+X(2,7)*H_invDmH(7,2)+X(2,6)*H_invDmH(6&
&,2)+X(2,5)*H_invDmH(5,2)+X(2,4)*H_invDmH(4,2)+X(2,3)*H_invDmH(3,2&
&)+H_invDmH(2,2)*X(2,2)+H_invDmH(1,2)*X(2,1)
tt3 = X(2,8)*H_invDmH(8,1)+X(2,7)*H_invDmH(7,1)+X(2,6)*H_invDmH(6&
&,1)+X(2,5)*H_invDmH(5,1)+X(2,4)*H_invDmH(4,1)+X(2,3)*H_invDmH(3,1&
&)+H_invDmH(2,1)*X(2,2)+H_invDmH(1,1)*X(2,1)
tt4 = X(3,8)*H_invDmH(8,2)+X(3,7)*H_invDmH(7,2)+X(3,6)*H_invDmH(6&
&,2)+X(3,5)*H_invDmH(5,2)+X(3,4)*H_invDmH(4,2)+H_invDmH(3,2)*X(3,3&
&)+H_invDmH(2,2)*X(3,2)+H_invDmH(1,2)*X(3,1)
tt5 = X(1,8)*H_invDmH(8,3)+X(1,7)*H_invDmH(7,3)+X(1,6)*H_invDmH(6&
&,3)+X(1,5)*H_invDmH(5,3)+X(1,4)*H_invDmH(4,3)+X(1,3)*H_invDmH(3,3&
&)+X(1,2)*H_invDmH(2,3)+X(1,1)*H_invDmH(1,3)
tt6 = (tt3*tt4-tt1*tt2)*tt5
tt7 = X(1,8)*H_invDmH(8,2)+X(1,7)*H_invDmH(7,2)+X(1,6)*H_invDmH(6&
&,2)+X(1,5)*H_invDmH(5,2)+X(1,4)*H_invDmH(4,2)+X(1,3)*H_invDmH(3,2&
&)+X(1,2)*H_invDmH(2,2)+X(1,1)*H_invDmH(1,2)
tt8 = X(2,8)*H_invDmH(8,3)+X(2,7)*H_invDmH(7,3)+X(2,6)*H_invDmH(6&
&,3)+X(2,5)*H_invDmH(5,3)+X(2,4)*H_invDmH(4,3)+X(2,3)*H_invDmH(3,3&
&)+X(2,2)*H_invDmH(2,3)+H_invDmH(1,3)*X(2,1)
tt9 = X(3,8)*H_invDmH(8,3)+X(3,7)*H_invDmH(7,3)+X(3,6)*H_invDmH(6&
&,3)+X(3,5)*H_invDmH(5,3)+X(3,4)*H_invDmH(4,3)+H_invDmH(3,3)*X(3,3&
&)+H_invDmH(2,3)*X(3,2)+H_invDmH(1,3)*X(3,1)
tt10 = -tt7*(tt3*tt9-tt1*tt8)
tt11 = X(1,8)*H_invDmH(8,1)+X(1,7)*H_invDmH(7,1)+X(1,6)*H_invDmH(&
&6,1)+X(1,5)*H_invDmH(5,1)+X(1,4)*H_invDmH(4,1)+X(1,3)*H_invDmH(3,&
&1)+X(1,2)*H_invDmH(2,1)+H_invDmH(1,1)*X(1,1)
tt12 = tt11*(tt2*tt9-tt4*tt8)
val(1,1) = detDmH(1,1)*gw(1,1)*((lam(1,1)*(tt12+tt10+tt6-1)**2)/2&
&.0E+0+(mu(1,1)*((tt9**2+tt8**2+tt5**2+tt4**2+tt2**2+tt7**2+tt1**2&
&+tt3**2+tt11**2)/(tt12+tt10+tt6)**(2.0E+0/3.0E+0)-3))/2.0E+0)
END 
SUBROUTINE vox_bower_neo_at_quadr_jac(jac, X, H_invDmH, detDmH, gw, mu, lam) 
IMPLICIT NONE 
REAL(KIND=8) jac(1, 24) 
REAL(KIND=8) X(3, 8) 
REAL(KIND=8) H_invDmH(8, 3) 
REAL(KIND=8) detDmH(1, 1) 
REAL(KIND=8) gw(1, 1) 
REAL(KIND=8) mu(1, 1) 
REAL(KIND=8) lam(1, 1) 
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
tt1 = X(3,8)*H_invDmH(8,1)+X(3,7)*H_invDmH(7,1)+X(3,6)*H_invDmH(6&
&,1)+X(3,5)*H_invDmH(5,1)+X(3,4)*H_invDmH(4,1)+H_invDmH(3,1)*X(3,3&
&)+H_invDmH(2,1)*X(3,2)+H_invDmH(1,1)*X(3,1)
tt2 = X(2,8)*H_invDmH(8,2)+X(2,7)*H_invDmH(7,2)+X(2,6)*H_invDmH(6&
&,2)+X(2,5)*H_invDmH(5,2)+X(2,4)*H_invDmH(4,2)+X(2,3)*H_invDmH(3,2&
&)+H_invDmH(2,2)*X(2,2)+H_invDmH(1,2)*X(2,1)
tt3 = X(2,8)*H_invDmH(8,1)+X(2,7)*H_invDmH(7,1)+X(2,6)*H_invDmH(6&
&,1)+X(2,5)*H_invDmH(5,1)+X(2,4)*H_invDmH(4,1)+X(2,3)*H_invDmH(3,1&
&)+H_invDmH(2,1)*X(2,2)+H_invDmH(1,1)*X(2,1)
tt4 = X(3,8)*H_invDmH(8,2)+X(3,7)*H_invDmH(7,2)+X(3,6)*H_invDmH(6&
&,2)+X(3,5)*H_invDmH(5,2)+X(3,4)*H_invDmH(4,2)+H_invDmH(3,2)*X(3,3&
&)+H_invDmH(2,2)*X(3,2)+H_invDmH(1,2)*X(3,1)
tt5 = tt3*tt4-tt1*tt2
tt6 = X(2,8)*H_invDmH(8,3)+X(2,7)*H_invDmH(7,3)+X(2,6)*H_invDmH(6&
&,3)+X(2,5)*H_invDmH(5,3)+X(2,4)*H_invDmH(4,3)+X(2,3)*H_invDmH(3,3&
&)+X(2,2)*H_invDmH(2,3)+H_invDmH(1,3)*X(2,1)
tt7 = X(3,8)*H_invDmH(8,3)+X(3,7)*H_invDmH(7,3)+X(3,6)*H_invDmH(6&
&,3)+X(3,5)*H_invDmH(5,3)+X(3,4)*H_invDmH(4,3)+H_invDmH(3,3)*X(3,3&
&)+H_invDmH(2,3)*X(3,2)+H_invDmH(1,3)*X(3,1)
tt8 = tt3*tt7-tt1*tt6
tt9 = tt2*tt7-tt4*tt6
tt10 = H_invDmH(1,1)*tt9-H_invDmH(1,2)*tt8+H_invDmH(1,3)*tt5
tt11 = X(1,8)*H_invDmH(8,3)+X(1,7)*H_invDmH(7,3)+X(1,6)*H_invDmH(&
&6,3)+X(1,5)*H_invDmH(5,3)+X(1,4)*H_invDmH(4,3)+X(1,3)*H_invDmH(3,&
&3)+X(1,2)*H_invDmH(2,3)+X(1,1)*H_invDmH(1,3)
tt12 = tt5*tt11
tt13 = X(1,8)*H_invDmH(8,2)+X(1,7)*H_invDmH(7,2)+X(1,6)*H_invDmH(&
&6,2)+X(1,5)*H_invDmH(5,2)+X(1,4)*H_invDmH(4,2)+X(1,3)*H_invDmH(3,&
&2)+X(1,2)*H_invDmH(2,2)+X(1,1)*H_invDmH(1,2)
tt14 = -tt13*tt8
tt15 = X(1,8)*H_invDmH(8,1)+X(1,7)*H_invDmH(7,1)+X(1,6)*H_invDmH(&
&6,1)+X(1,5)*H_invDmH(5,1)+X(1,4)*H_invDmH(4,1)+X(1,3)*H_invDmH(3,&
&1)+X(1,2)*H_invDmH(2,1)+H_invDmH(1,1)*X(1,1)
tt16 = tt15*tt9
tt17 = tt16+tt14+tt12-1
tt18 = tt16+tt14+tt12
tt19 = tt18**((-5.0E+0)/3.0E+0)
tt20 = tt7**2+tt6**2+tt11**2+tt4**2+tt2**2+tt13**2+tt1**2+tt3**2+&
&tt15**2
tt21 = tt18**((-2.0E+0)/3.0E+0)
tt22 = tt15*(H_invDmH(1,2)*tt7-H_invDmH(1,3)*tt4)-tt13*(H_invDmH(&
&1,1)*tt7-H_invDmH(1,3)*tt1)+(H_invDmH(1,1)*tt4-H_invDmH(1,2)*tt1)&
&*tt11
tt23 = tt15*(H_invDmH(1,3)*tt2-H_invDmH(1,2)*tt6)-tt13*(H_invDmH(&
&1,3)*tt3-H_invDmH(1,1)*tt6)+(H_invDmH(1,2)*tt3-H_invDmH(1,1)*tt2)&
&*tt11
tt24 = H_invDmH(2,1)*tt9-H_invDmH(2,2)*tt8+H_invDmH(2,3)*tt5
tt25 = tt15*(H_invDmH(2,2)*tt7-H_invDmH(2,3)*tt4)-tt13*(H_invDmH(&
&2,1)*tt7-H_invDmH(2,3)*tt1)+(H_invDmH(2,1)*tt4-H_invDmH(2,2)*tt1)&
&*tt11
tt26 = tt15*(H_invDmH(2,3)*tt2-H_invDmH(2,2)*tt6)-tt13*(H_invDmH(&
&2,3)*tt3-H_invDmH(2,1)*tt6)+(H_invDmH(2,2)*tt3-H_invDmH(2,1)*tt2)&
&*tt11
tt27 = H_invDmH(3,1)*tt9-H_invDmH(3,2)*tt8+H_invDmH(3,3)*tt5
tt28 = tt15*(H_invDmH(3,2)*tt7-H_invDmH(3,3)*tt4)-tt13*(H_invDmH(&
&3,1)*tt7-H_invDmH(3,3)*tt1)+(H_invDmH(3,1)*tt4-H_invDmH(3,2)*tt1)&
&*tt11
tt29 = tt15*(H_invDmH(3,3)*tt2-H_invDmH(3,2)*tt6)-tt13*(H_invDmH(&
&3,3)*tt3-H_invDmH(3,1)*tt6)+(H_invDmH(3,2)*tt3-H_invDmH(3,1)*tt2)&
&*tt11
tt30 = H_invDmH(4,1)*tt9-H_invDmH(4,2)*tt8+H_invDmH(4,3)*tt5
tt31 = tt15*(H_invDmH(4,2)*tt7-H_invDmH(4,3)*tt4)-tt13*(H_invDmH(&
&4,1)*tt7-H_invDmH(4,3)*tt1)+(H_invDmH(4,1)*tt4-H_invDmH(4,2)*tt1)&
&*tt11
tt32 = tt15*(H_invDmH(4,3)*tt2-H_invDmH(4,2)*tt6)-tt13*(H_invDmH(&
&4,3)*tt3-H_invDmH(4,1)*tt6)+(H_invDmH(4,2)*tt3-H_invDmH(4,1)*tt2)&
&*tt11
tt33 = H_invDmH(5,1)*tt9-H_invDmH(5,2)*tt8+H_invDmH(5,3)*tt5
tt34 = tt15*(H_invDmH(5,2)*tt7-H_invDmH(5,3)*tt4)-tt13*(H_invDmH(&
&5,1)*tt7-H_invDmH(5,3)*tt1)+(H_invDmH(5,1)*tt4-H_invDmH(5,2)*tt1)&
&*tt11
tt35 = tt15*(H_invDmH(5,3)*tt2-H_invDmH(5,2)*tt6)-tt13*(H_invDmH(&
&5,3)*tt3-H_invDmH(5,1)*tt6)+(H_invDmH(5,2)*tt3-H_invDmH(5,1)*tt2)&
&*tt11
tt36 = H_invDmH(6,1)*tt9-H_invDmH(6,2)*tt8+H_invDmH(6,3)*tt5
tt37 = tt15*(H_invDmH(6,2)*tt7-H_invDmH(6,3)*tt4)-tt13*(H_invDmH(&
&6,1)*tt7-H_invDmH(6,3)*tt1)+(H_invDmH(6,1)*tt4-H_invDmH(6,2)*tt1)&
&*tt11
tt38 = tt15*(H_invDmH(6,3)*tt2-H_invDmH(6,2)*tt6)-tt13*(H_invDmH(&
&6,3)*tt3-H_invDmH(6,1)*tt6)+(H_invDmH(6,2)*tt3-H_invDmH(6,1)*tt2)&
&*tt11
tt39 = H_invDmH(7,1)*tt9-H_invDmH(7,2)*tt8+H_invDmH(7,3)*tt5
tt40 = tt15*(H_invDmH(7,2)*tt7-H_invDmH(7,3)*tt4)-tt13*(H_invDmH(&
&7,1)*tt7-H_invDmH(7,3)*tt1)+(H_invDmH(7,1)*tt4-H_invDmH(7,2)*tt1)&
&*tt11
tt41 = tt15*(H_invDmH(7,3)*tt2-H_invDmH(7,2)*tt6)-tt13*(H_invDmH(&
&7,3)*tt3-H_invDmH(7,1)*tt6)+(H_invDmH(7,2)*tt3-H_invDmH(7,1)*tt2)&
&*tt11
tt42 = H_invDmH(8,1)*tt9-H_invDmH(8,2)*tt8+tt5*H_invDmH(8,3)
tt43 = tt15*(H_invDmH(8,2)*tt7-tt4*H_invDmH(8,3))-tt13*(H_invDmH(&
&8,1)*tt7-tt1*H_invDmH(8,3))+(H_invDmH(8,1)*tt4-tt1*H_invDmH(8,2))&
&*tt11
tt44 = tt15*(tt2*H_invDmH(8,3)-H_invDmH(8,2)*tt6)-tt13*(tt3*H_inv&
&DmH(8,3)-H_invDmH(8,1)*tt6)+(tt3*H_invDmH(8,2)-H_invDmH(8,1)*tt2)&
&*tt11
jac(1,1) = detDmH(1,1)*gw(1,1)*((mu(1,1)*((2*H_invDmH(1,3)*tt11+2&
&*H_invDmH(1,2)*tt13+2*H_invDmH(1,1)*tt15)*tt21+((-2.0E+0)*tt10*tt&
&19*tt20)/3.0E+0))/2.0E+0+lam(1,1)*tt10*tt17)
jac(1,2) = detDmH(1,1)*gw(1,1)*((mu(1,1)*((2*H_invDmH(1,3)*tt6+2*&
&H_invDmH(1,2)*tt2+2*H_invDmH(1,1)*tt3)*tt21+((-2.0E+0)*tt22*tt19*&
&tt20)/3.0E+0))/2.0E+0+lam(1,1)*tt22*tt17)
jac(1,3) = detDmH(1,1)*gw(1,1)*((mu(1,1)*((2*H_invDmH(1,3)*tt7+2*&
&H_invDmH(1,2)*tt4+2*H_invDmH(1,1)*tt1)*tt21+((-2.0E+0)*tt23*tt19*&
&tt20)/3.0E+0))/2.0E+0+lam(1,1)*tt23*tt17)
jac(1,4) = detDmH(1,1)*gw(1,1)*((mu(1,1)*((2*H_invDmH(2,3)*tt11+2&
&*H_invDmH(2,2)*tt13+2*H_invDmH(2,1)*tt15)*tt21+((-2.0E+0)*tt24*tt&
&19*tt20)/3.0E+0))/2.0E+0+lam(1,1)*tt24*tt17)
jac(1,5) = detDmH(1,1)*gw(1,1)*((mu(1,1)*((2*H_invDmH(2,3)*tt6+2*&
&H_invDmH(2,2)*tt2+2*H_invDmH(2,1)*tt3)*tt21+((-2.0E+0)*tt25*tt19*&
&tt20)/3.0E+0))/2.0E+0+lam(1,1)*tt25*tt17)
jac(1,6) = detDmH(1,1)*gw(1,1)*((mu(1,1)*((2*H_invDmH(2,3)*tt7+2*&
&H_invDmH(2,2)*tt4+2*H_invDmH(2,1)*tt1)*tt21+((-2.0E+0)*tt26*tt19*&
&tt20)/3.0E+0))/2.0E+0+lam(1,1)*tt26*tt17)
jac(1,7) = detDmH(1,1)*gw(1,1)*((mu(1,1)*((2*H_invDmH(3,3)*tt11+2&
&*H_invDmH(3,2)*tt13+2*H_invDmH(3,1)*tt15)*tt21+((-2.0E+0)*tt27*tt&
&19*tt20)/3.0E+0))/2.0E+0+lam(1,1)*tt27*tt17)
jac(1,8) = detDmH(1,1)*gw(1,1)*((mu(1,1)*((2*H_invDmH(3,3)*tt6+2*&
&H_invDmH(3,2)*tt2+2*H_invDmH(3,1)*tt3)*tt21+((-2.0E+0)*tt28*tt19*&
&tt20)/3.0E+0))/2.0E+0+lam(1,1)*tt28*tt17)
jac(1,9) = detDmH(1,1)*gw(1,1)*((mu(1,1)*((2*H_invDmH(3,3)*tt7+2*&
&H_invDmH(3,2)*tt4+2*H_invDmH(3,1)*tt1)*tt21+((-2.0E+0)*tt29*tt19*&
&tt20)/3.0E+0))/2.0E+0+lam(1,1)*tt29*tt17)
jac(1,10) = detDmH(1,1)*gw(1,1)*((mu(1,1)*((2*H_invDmH(4,3)*tt11+&
&2*H_invDmH(4,2)*tt13+2*H_invDmH(4,1)*tt15)*tt21+((-2.0E+0)*tt30*t&
&t19*tt20)/3.0E+0))/2.0E+0+lam(1,1)*tt30*tt17)
jac(1,11) = detDmH(1,1)*gw(1,1)*((mu(1,1)*((2*H_invDmH(4,3)*tt6+2&
&*H_invDmH(4,2)*tt2+2*H_invDmH(4,1)*tt3)*tt21+((-2.0E+0)*tt31*tt19&
&*tt20)/3.0E+0))/2.0E+0+lam(1,1)*tt31*tt17)
jac(1,12) = detDmH(1,1)*gw(1,1)*((mu(1,1)*((2*H_invDmH(4,3)*tt7+2&
&*H_invDmH(4,2)*tt4+2*H_invDmH(4,1)*tt1)*tt21+((-2.0E+0)*tt32*tt19&
&*tt20)/3.0E+0))/2.0E+0+lam(1,1)*tt32*tt17)
jac(1,13) = detDmH(1,1)*gw(1,1)*((mu(1,1)*((2*H_invDmH(5,3)*tt11+&
&2*H_invDmH(5,2)*tt13+2*H_invDmH(5,1)*tt15)*tt21+((-2.0E+0)*tt33*t&
&t19*tt20)/3.0E+0))/2.0E+0+lam(1,1)*tt33*tt17)
jac(1,14) = detDmH(1,1)*gw(1,1)*((mu(1,1)*((2*H_invDmH(5,3)*tt6+2&
&*H_invDmH(5,2)*tt2+2*H_invDmH(5,1)*tt3)*tt21+((-2.0E+0)*tt34*tt19&
&*tt20)/3.0E+0))/2.0E+0+lam(1,1)*tt34*tt17)
jac(1,15) = detDmH(1,1)*gw(1,1)*((mu(1,1)*((2*H_invDmH(5,3)*tt7+2&
&*H_invDmH(5,2)*tt4+2*H_invDmH(5,1)*tt1)*tt21+((-2.0E+0)*tt35*tt19&
&*tt20)/3.0E+0))/2.0E+0+lam(1,1)*tt35*tt17)
jac(1,16) = detDmH(1,1)*gw(1,1)*((mu(1,1)*((2*H_invDmH(6,3)*tt11+&
&2*H_invDmH(6,2)*tt13+2*H_invDmH(6,1)*tt15)*tt21+((-2.0E+0)*tt36*t&
&t19*tt20)/3.0E+0))/2.0E+0+lam(1,1)*tt36*tt17)
jac(1,17) = detDmH(1,1)*gw(1,1)*((mu(1,1)*((2*H_invDmH(6,3)*tt6+2&
&*H_invDmH(6,2)*tt2+2*H_invDmH(6,1)*tt3)*tt21+((-2.0E+0)*tt37*tt19&
&*tt20)/3.0E+0))/2.0E+0+lam(1,1)*tt37*tt17)
jac(1,18) = detDmH(1,1)*gw(1,1)*((mu(1,1)*((2*H_invDmH(6,3)*tt7+2&
&*H_invDmH(6,2)*tt4+2*H_invDmH(6,1)*tt1)*tt21+((-2.0E+0)*tt38*tt19&
&*tt20)/3.0E+0))/2.0E+0+lam(1,1)*tt38*tt17)
jac(1,19) = detDmH(1,1)*gw(1,1)*((mu(1,1)*((2*H_invDmH(7,3)*tt11+&
&2*H_invDmH(7,2)*tt13+2*H_invDmH(7,1)*tt15)*tt21+((-2.0E+0)*tt39*t&
&t19*tt20)/3.0E+0))/2.0E+0+lam(1,1)*tt39*tt17)
jac(1,20) = detDmH(1,1)*gw(1,1)*((mu(1,1)*((2*H_invDmH(7,3)*tt6+2&
&*H_invDmH(7,2)*tt2+2*H_invDmH(7,1)*tt3)*tt21+((-2.0E+0)*tt40*tt19&
&*tt20)/3.0E+0))/2.0E+0+lam(1,1)*tt40*tt17)
jac(1,21) = detDmH(1,1)*gw(1,1)*((mu(1,1)*((2*H_invDmH(7,3)*tt7+2&
&*H_invDmH(7,2)*tt4+2*H_invDmH(7,1)*tt1)*tt21+((-2.0E+0)*tt41*tt19&
&*tt20)/3.0E+0))/2.0E+0+lam(1,1)*tt41*tt17)
jac(1,22) = detDmH(1,1)*gw(1,1)*((mu(1,1)*((2*H_invDmH(8,3)*tt11+&
&2*H_invDmH(8,2)*tt13+2*H_invDmH(8,1)*tt15)*tt21+((-2.0E+0)*tt42*t&
&t19*tt20)/3.0E+0))/2.0E+0+lam(1,1)*tt42*tt17)
jac(1,23) = detDmH(1,1)*gw(1,1)*((mu(1,1)*((2*H_invDmH(8,3)*tt6+2&
&*H_invDmH(8,2)*tt2+2*H_invDmH(8,1)*tt3)*tt21+((-2.0E+0)*tt43*tt19&
&*tt20)/3.0E+0))/2.0E+0+lam(1,1)*tt43*tt17)
jac(1,24) = detDmH(1,1)*gw(1,1)*((mu(1,1)*((2*H_invDmH(8,3)*tt7+2&
&*H_invDmH(8,2)*tt4+2*H_invDmH(8,1)*tt1)*tt21+((-2.0E+0)*tt44*tt19&
&*tt20)/3.0E+0))/2.0E+0+lam(1,1)*tt44*tt17)
END 
SUBROUTINE vox_bower_neo_at_quadr_hes(hes, X, H_invDmH, detDmH, gw, mu, lam) 
IMPLICIT NONE 
