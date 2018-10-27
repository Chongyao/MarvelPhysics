SUBROUTINE quad_lin_q2l_psi(val, F, mu, lam) 
IMPLICIT NONE 
REAL(KIND=8) val(1, 1) 
REAL(KIND=8) F(2, 4) 
REAL(KIND=8) mu(1, 1) 
REAL(KIND=8) lam(1, 1) 
REAL(KIND=8)  tt1 
REAL(KIND=8)  tt2 
REAL(KIND=8)  tt3 
tt1 = -4*lam(1,1)
tt2 = -4*mu(1,1)
tt3 = 2*mu(1,1)+lam(1,1)
val(1,1) = (tt3*F(2,4)**2+(2*lam(1,1)*F(1,3)+tt2+tt1)*F(2,4)+mu(1&
&,1)*F(2,3)**2+2*mu(1,1)*F(1,4)*F(2,3)+mu(1,1)*F(1,4)**2+tt3*F(1,3&
&)**2+(tt2+tt1)*F(1,3)+4*mu(1,1)+4*lam(1,1))/2.0E+0
END 
SUBROUTINE quad_lin_q2l_psi_jac(jac, F, mu, lam) 
IMPLICIT NONE 
REAL(KIND=8) jac(1, 8) 
REAL(KIND=8) F(2, 4) 
REAL(KIND=8) mu(1, 1) 
REAL(KIND=8) lam(1, 1) 
REAL(KIND=8)  tt1 
REAL(KIND=8)  tt2 
REAL(KIND=8)  tt3 
REAL(KIND=8)  tt4 
tt1 = -4*lam(1,1)
tt2 = -4*mu(1,1)
tt3 = 2*mu(1,1)+lam(1,1)
tt4 = (2*mu(1,1)*F(2,3)+2*mu(1,1)*F(1,4))/2.0E+0
jac(1,1) = 0
jac(1,2) = 0
jac(1,3) = 0
jac(1,4) = 0
jac(1,5) = (2*lam(1,1)*F(2,4)+2*tt3*F(1,3)+tt2+tt1)/2.0E+0
jac(1,6) = tt4
jac(1,7) = tt4
jac(1,8) = (2*tt3*F(2,4)+2*lam(1,1)*F(1,3)+tt2+tt1)/2.0E+0
END 
SUBROUTINE quad_coro_F_val_at_qr(val, F, R, mu, lam) 
IMPLICIT NONE 
REAL(KIND=8) val(1, 1) 
REAL(KIND=8) F(2, 2) 
REAL(KIND=8) R(2, 2) 
REAL(KIND=8) mu(1, 1) 
REAL(KIND=8) lam(1, 1) 
val(1,1) = (lam(1,1)*((2*F(2,2)*R(2,2)+2*F(1,2)*R(1,2))/2.0E+0+(2&
&*F(2,1)*R(2,1)+2*F(1,1)*R(1,1))/2.0E+0-2)**2)/2.0E+0+mu(1,1)*((F(&
&2,2)*R(2,2)+F(1,2)*R(1,2)-1)**2+(F(2,1)*R(2,2)+R(2,1)*F(2,2)+F(1,&
&1)*R(1,2)+R(1,1)*F(1,2))**2/2.0E+0+(F(2,1)*R(2,1)+F(1,1)*R(1,1)-1&
&)**2)
END 
SUBROUTINE quad_coro_F_val_at_qr_jac(jac, F, R, mu, lam) 
IMPLICIT NONE 
REAL(KIND=8) jac(1, 4) 
REAL(KIND=8) F(2, 2) 
REAL(KIND=8) R(2, 2) 
REAL(KIND=8) mu(1, 1) 
REAL(KIND=8) lam(1, 1) 
REAL(KIND=8)  tt1 
REAL(KIND=8)  tt2 
REAL(KIND=8)  tt3 
REAL(KIND=8)  tt4 
tt1 = F(2,1)*R(2,1)+F(1,1)*R(1,1)-1
tt2 = F(2,1)*R(2,2)+R(2,1)*F(2,2)+F(1,1)*R(1,2)+R(1,1)*F(1,2)
tt3 = (2*F(2,2)*R(2,2)+2*F(1,2)*R(1,2))/2.0E+0+(2*F(2,1)*R(2,1)+2&
&*F(1,1)*R(1,1))/2.0E+0-2
tt4 = F(2,2)*R(2,2)+F(1,2)*R(1,2)-1
jac(1,1) = R(1,1)*lam(1,1)*tt3+mu(1,1)*(R(1,2)*tt2+2*R(1,1)*tt1)
jac(1,2) = lam(1,1)*R(2,1)*tt3+mu(1,1)*(R(2,2)*tt2+2*R(2,1)*tt1)
jac(1,3) = lam(1,1)*R(1,2)*tt3+mu(1,1)*(2*R(1,2)*tt4+R(1,1)*tt2)
jac(1,4) = lam(1,1)*R(2,2)*tt3+mu(1,1)*(2*R(2,2)*tt4+R(2,1)*tt2)
END 
SUBROUTINE quad_coro_F_val_at_qr_hes(hes, F, R, mu, lam) 
IMPLICIT NONE 
REAL(KIND=8) hes(4, 4) 
REAL(KIND=8) F(2, 2) 
REAL(KIND=8) R(2, 2) 
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
tt1 = R(1,1)**2
tt2 = R(1,2)**2
tt3 = mu(1,1)*(R(1,2)*R(2,2)+2*R(1,1)*R(2,1))+R(1,1)*lam(1,1)*R(2&
&,1)
tt4 = R(1,1)*mu(1,1)*R(1,2)+R(1,1)*lam(1,1)*R(1,2)
tt5 = R(1,1)*lam(1,1)*R(2,2)+mu(1,1)*R(1,2)*R(2,1)
tt6 = R(2,1)**2
tt7 = R(2,2)**2
tt8 = R(1,1)*mu(1,1)*R(2,2)+lam(1,1)*R(1,2)*R(2,1)
tt9 = mu(1,1)*R(2,1)*R(2,2)+lam(1,1)*R(2,1)*R(2,2)
tt10 = mu(1,1)*(2*R(1,2)*R(2,2)+R(1,1)*R(2,1))+lam(1,1)*R(1,2)*R(&
&2,2)
hes(1,1) = mu(1,1)*(tt2+2*tt1)+tt1*lam(1,1)
hes(1,2) = tt3
hes(1,3) = tt4
hes(1,4) = tt5
hes(2,1) = tt3
hes(2,2) = mu(1,1)*(tt7+2*tt6)+lam(1,1)*tt6
hes(2,3) = tt8
hes(2,4) = tt9
hes(3,1) = tt4
hes(3,2) = tt8
hes(3,3) = mu(1,1)*(2*tt2+tt1)+lam(1,1)*tt2
hes(3,4) = tt10
hes(4,1) = tt5
hes(4,2) = tt9
hes(4,3) = tt10
hes(4,4) = mu(1,1)*(2*tt7+tt6)+lam(1,1)*tt7
END 
SUBROUTINE quad_neo_F_val_at_qr(val, F2, mu, lam) 
IMPLICIT NONE 
REAL(KIND=8) val(1, 1) 
REAL(KIND=8) F2(2, 2) 
REAL(KIND=8) mu(1, 1) 
REAL(KIND=8) lam(1, 1) 
REAL(KIND=8)  tt1 
REAL(KIND=8)  tt2 
tt1 = -F2(1,2)*F2(2,1)
tt2 = F2(1,1)*F2(2,2)
val(1,1) = (lam(1,1)*(tt2+tt1-1)**2+mu(1,1)*((F2(2,2)**2+F2(2,1)*&
&*2+F2(1,2)**2+F2(1,1)**2)/(tt2+tt1)-2))/2.0E+0
END 
SUBROUTINE quad_neo_F_val_at_qr_jac(jac, F2, mu, lam) 
IMPLICIT NONE 
REAL(KIND=8) jac(1, 4) 
REAL(KIND=8) F2(2, 2) 
REAL(KIND=8) mu(1, 1) 
REAL(KIND=8) lam(1, 1) 
REAL(KIND=8)  tt1 
REAL(KIND=8)  tt2 
REAL(KIND=8)  tt3 
REAL(KIND=8)  tt4 
REAL(KIND=8)  tt5 
REAL(KIND=8)  tt6 
REAL(KIND=8)  tt7 
tt1 = -F2(1,2)*F2(2,1)
tt2 = F2(1,1)*F2(2,2)
tt3 = tt2+tt1-1
tt4 = tt2+tt1
tt5 = 1/tt4**2
tt6 = F2(2,2)**2+F2(2,1)**2+F2(1,2)**2+F2(1,1)**2
tt7 = 1/tt4
jac(1,1) = (mu(1,1)*(2*F2(1,1)*tt7-F2(2,2)*tt5*tt6)+2*lam(1,1)*F2&
&(2,2)*tt3)/2.0E+0
jac(1,2) = (mu(1,1)*(2*F2(2,1)*tt7+F2(1,2)*tt5*tt6)-2*lam(1,1)*F2&
&(1,2)*tt3)/2.0E+0
jac(1,3) = (mu(1,1)*(2*F2(1,2)*tt7+F2(2,1)*tt5*tt6)-2*lam(1,1)*F2&
&(2,1)*tt3)/2.0E+0
jac(1,4) = (mu(1,1)*(2*F2(2,2)*tt7-F2(1,1)*tt5*tt6)+2*F2(1,1)*lam&
&(1,1)*tt3)/2.0E+0
END 
SUBROUTINE quad_neo_F_val_at_qr_hes(hes, F2, mu, lam) 
IMPLICIT NONE 
REAL(KIND=8) hes(4, 4) 
REAL(KIND=8) F2(2, 2) 
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
tt1 = F2(2,2)**2
tt2 = -F2(1,2)*F2(2,1)
tt3 = F2(1,1)*F2(2,2)
tt4 = tt3+tt2
tt5 = 1/tt4**3
tt6 = F2(1,1)**2
tt7 = F2(1,2)**2
tt8 = F2(2,1)**2
tt9 = tt1+tt8+tt7+tt6
tt10 = 1/tt4**2
tt11 = -4*F2(1,1)*F2(2,2)*tt10
tt12 = 2/tt4
tt13 = (mu(1,1)*((-2*F2(2,1)*F2(2,2)*tt10)+2*F2(1,1)*F2(1,2)*tt10&
&-2*F2(1,2)*F2(2,2)*tt5*tt9)-2*lam(1,1)*F2(1,2)*F2(2,2))/2.0E+0
tt14 = (mu(1,1)*((-2*F2(1,2)*F2(2,2)*tt10)+2*F2(1,1)*F2(2,1)*tt10&
&-2*F2(2,1)*F2(2,2)*tt5*tt9)-2*lam(1,1)*F2(2,1)*F2(2,2))/2.0E+0
tt15 = tt3+tt2-1
tt16 = (mu(1,1)*((-2*tt1*tt10)-2*tt6*tt10-tt10*tt9+2*F2(1,1)*F2(2&
&,2)*tt5*tt9)+2*lam(1,1)*tt15+2*F2(1,1)*lam(1,1)*F2(2,2))/2.0E+0
tt17 = 4*F2(1,2)*F2(2,1)*tt10
tt18 = (mu(1,1)*(2*tt8*tt10+2*tt7*tt10+tt10*tt9+2*F2(1,2)*F2(2,1)&
&*tt5*tt9)-2*lam(1,1)*tt15+2*lam(1,1)*F2(1,2)*F2(2,1))/2.0E+0
tt19 = (mu(1,1)*(2*F2(1,2)*F2(2,2)*tt10-2*F2(1,1)*F2(2,1)*tt10-2*&
&F2(1,1)*F2(1,2)*tt5*tt9)-2*F2(1,1)*lam(1,1)*F2(1,2))/2.0E+0
tt20 = (mu(1,1)*(2*F2(2,1)*F2(2,2)*tt10-2*F2(1,1)*F2(1,2)*tt10-2*&
&F2(1,1)*F2(2,1)*tt5*tt9)-2*F2(1,1)*lam(1,1)*F2(2,1))/2.0E+0
hes(1,1) = (mu(1,1)*(tt12+tt11+2*tt1*tt5*tt9)+2*lam(1,1)*tt1)/2.0&
&E+0
hes(1,2) = tt13
hes(1,3) = tt14
hes(1,4) = tt16
hes(2,1) = tt13
hes(2,2) = (mu(1,1)*(tt12+tt17+2*tt7*tt5*tt9)+2*lam(1,1)*tt7)/2.0&
&E+0
hes(2,3) = tt18
hes(2,4) = tt19
hes(3,1) = tt14
hes(3,2) = tt18
hes(3,3) = (mu(1,1)*(tt12+tt17+2*tt8*tt5*tt9)+2*lam(1,1)*tt8)/2.0&
&E+0
hes(3,4) = tt20
hes(4,1) = tt16
hes(4,2) = tt19
hes(4,3) = tt20
hes(4,4) = (mu(1,1)*(tt12+tt11+2*tt6*tt5*tt9)+2*tt6*lam(1,1))/2.0&
&E+0
END 
SUBROUTINE quad_neo_q2l_psi(val, F, mu, lam) 
IMPLICIT NONE 
REAL(KIND=8) val(1, 1) 
REAL(KIND=8) F(2, 4) 
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
tt1 = F(1,2)**3
tt2 = F(2,1)**3
tt3 = F(1,2)**2
tt4 = F(2,1)**2
tt5 = F(1,1)**2
tt6 = F(2,2)**2
tt7 = F(1,1)**3
tt8 = F(2,2)**3
tt9 = tt5*mu(1,1)*tt3
tt10 = F(1,2)**4
tt11 = mu(1,1)*tt10
tt12 = F(1,3)**2
tt13 = mu(1,1)*tt3*tt12
tt14 = 2*mu(1,1)*tt1
tt15 = -2*F(1,1)*mu(1,1)*F(1,2)*F(1,3)
tt16 = F(1,4)**2
tt17 = tt5*mu(1,1)*tt16
tt18 = mu(1,1)*tt3
tt19 = -2*mu(1,1)*F(1,2)*F(1,4)
tt20 = mu(1,1)*tt16
tt21 = F(2,1)**4
tt22 = F(1,2)**5
tt23 = -2*tt7*mu(1,1)*F(1,2)
tt24 = -2*F(1,1)*mu(1,1)*tt1
tt25 = 2*tt7*mu(1,1)
tt26 = (2*F(1,1)*mu(1,1)-2*mu(1,1)*F(1,3))*F(1,4)+2*mu(1,1)*F(1,2&
&)*F(1,3)-2*F(1,1)*mu(1,1)*F(1,2)
tt27 = F(1,1)**4
tt28 = tt27*mu(1,1)
tt29 = tt5*mu(1,1)
tt30 = -2*F(1,1)*mu(1,1)*F(1,3)
tt31 = mu(1,1)*tt12
tt32 = F(2,2)**4
tt33 = F(1,1)**5
tt34 = -2*tt5*mu(1,1)*tt3
tt35 = 2*tt7*mu(1,1)*F(1,2)
tt36 = 2*F(1,1)*mu(1,1)*tt1
tt37 = -6*F(1,1)*mu(1,1)*F(1,2)
tt38 = 2*lam(1,1)*tt10*F(1,3)
tt39 = -2*mu(1,1)*tt3
tt40 = 2*F(1,1)*mu(1,1)*F(1,2)
tt41 = F(1,1)*mu(1,1)
tt42 = -2*tt27*lam(1,1)*F(1,4)
tt43 = -2*tt5*mu(1,1)
val(1,1) = ((tt33*lam(1,1)*tt8-3*tt27*lam(1,1)*F(1,2)*F(2,1)*tt6+&
&3*tt7*lam(1,1)*tt3*tt4*F(2,2)-tt5*lam(1,1)*tt1*tt2+(tt18+tt29)*tt&
&4+tt9+tt28)*F(2,4)**2+(((-2*tt27*lam(1,1)*F(1,2)*tt8)+6*tt7*lam(1&
&,1)*tt3*F(2,1)*tt6+((tt39+tt43)*F(2,1)-6*tt5*lam(1,1)*tt1*tt4)*F(&
&2,2)+2*F(1,1)*lam(1,1)*tt10*tt2+tt24+tt23)*F(2,3)+(4*tt27*lam(1,1&
&)*F(1,3)-4*tt33*lam(1,1))*tt32+((tt42-14*tt7*lam(1,1)*F(1,2)*F(1,&
&3)+16*tt27*lam(1,1)*F(1,2))*F(2,1)+((-F(1,1)*mu(1,1))-2*tt7*lam(1&
&,1))*F(1,3)+2*tt5*mu(1,1))*tt8+((6*tt7*lam(1,1)*F(1,2)*F(1,4)+18*&
&tt5*lam(1,1)*tt3*F(1,3)-24*tt7*lam(1,1)*tt3)*tt4+((3*mu(1,1)+6*tt&
&5*lam(1,1))*F(1,2)*F(1,3)+tt37)*F(2,1))*tt6+(((-6*tt5*lam(1,1)*tt&
&3*F(1,4))-10*F(1,1)*lam(1,1)*tt1*F(1,3)+16*tt5*lam(1,1)*tt1)*tt2+&
&(tt19+(tt41-6*F(1,1)*lam(1,1)*tt3)*F(1,3)+4*mu(1,1)*tt3+tt43)*tt4&
&-2*tt5*mu(1,1)*F(1,2)*F(1,4)+(F(1,1)*mu(1,1)*tt3-tt7*mu(1,1))*F(1&
&,3)+tt34-2*tt27*mu(1,1))*F(2,2)+(2*F(1,1)*lam(1,1)*tt1*F(1,4)+tt3&
&8-4*F(1,1)*lam(1,1)*tt10)*tt21+((-2*F(1,1)*mu(1,1)*F(1,4))+(2*lam&
&(1,1)*tt1+mu(1,1)*F(1,2))*F(1,3)+tt40)*tt2+((-2*tt7*mu(1,1)*F(1,4&
&))+(mu(1,1)*tt1+3*tt5*mu(1,1)*F(1,2))*F(1,3)+tt36+tt35)*F(2,1))*F&
&(2,4)+(tt7*lam(1,1)*tt3*tt8+((-3*tt5*lam(1,1)*tt1*F(2,1))+tt18+tt&
&29)*tt6+3*F(1,1)*lam(1,1)*tt10*tt4*F(2,2)-lam(1,1)*tt22*tt2+tt11+&
&tt9)*F(2,3)**2+((tt42-2*tt7*lam(1,1)*F(1,2)*F(1,3)+4*tt27*lam(1,1&
&)*F(1,2))*tt32+((10*tt7*lam(1,1)*F(1,2)*F(1,4)+6*tt5*lam(1,1)*tt3&
&*F(1,3)-16*tt7*lam(1,1)*tt3)*F(2,1)+(tt41+2*tt7*lam(1,1))*F(1,4)-&
&2*mu(1,1)*F(1,2)*F(1,3)+tt40)*tt8+(((-18*tt5*lam(1,1)*tt3*F(1,4))&
&-6*F(1,1)*lam(1,1)*tt1*F(1,3)+24*tt5*lam(1,1)*tt1)*tt4+((mu(1,1)-&
&6*tt5*lam(1,1))*F(1,2)*F(1,4)+tt30+tt39+4*tt5*mu(1,1))*F(2,1))*tt&
&6+((14*F(1,1)*lam(1,1)*tt1*F(1,4)+tt38-16*F(1,1)*lam(1,1)*tt10)*t&
&t2+((6*F(1,1)*lam(1,1)*tt3+3*F(1,1)*mu(1,1))*F(1,4)+tt37)*tt4+(3*&
&F(1,1)*mu(1,1)*tt3+tt7*mu(1,1))*F(1,4)-2*mu(1,1)*tt1*F(1,3)+tt36+&
&tt35)*F(2,2)+(4*lam(1,1)*tt22-4*lam(1,1)*tt10*F(1,4))*tt21+(((-2*&
&lam(1,1)*tt1)-mu(1,1)*F(1,2))*F(1,4)+2*mu(1,1)*tt3)*tt2+((tt5*mu(&
&1,1)*F(1,2)-mu(1,1)*tt1)*F(1,4)-2*F(1,1)*mu(1,1)*tt3*F(1,3)-2*mu(&
&1,1)*tt10+tt34)*F(2,1))*F(2,3)+(tt7*lam(1,1)*tt12-4*tt27*lam(1,1)&
&*F(1,3)+3*tt33*lam(1,1))*F(2,2)**5+(((4*tt27*lam(1,1)-2*tt7*lam(1&
&,1)*F(1,3))*F(1,4)-3*tt5*lam(1,1)*F(1,2)*tt12+16*tt7*lam(1,1)*F(1&
&,2)*F(1,3)-15*tt27*lam(1,1)*F(1,2))*F(2,1)+tt31+tt30+tt29)*tt32+(&
&(tt7*lam(1,1)*tt16+(6*tt5*lam(1,1)*F(1,2)*F(1,3)-16*tt7*lam(1,1)*&
&F(1,2))*F(1,4)+3*F(1,1)*lam(1,1)*tt3*tt12-24*tt5*lam(1,1)*tt3*F(1&
&,3)+30*tt7*lam(1,1)*tt3)*tt4+tt26*F(2,1)-2*tt7*mu(1,1)+tt7*lam(1,&
&1))*tt8+(((-3*tt5*lam(1,1)*F(1,2)*tt16)+(24*tt5*lam(1,1)*tt3-6*F(&
&1,1)*lam(1,1)*tt3*F(1,3))*F(1,4)-lam(1,1)*tt1*tt12+16*F(1,1)*lam(&
&1,1)*tt1*F(1,3)-30*tt5*lam(1,1)*tt1)*tt2+(tt20+tt19+tt31+tt30+tt1&
&8+tt29)*tt4+(6*tt5*mu(1,1)-3*tt5*lam(1,1))*F(1,2)*F(2,1)+tt17+(tt&
&15+4*tt5*mu(1,1)*F(1,2))*F(1,4)+tt13+(tt25-2*F(1,1)*mu(1,1)*tt3)*&
&F(1,3)+tt9+tt28)*tt6+((3*F(1,1)*lam(1,1)*tt3*tt16+(2*lam(1,1)*tt1&
&*F(1,3)-16*F(1,1)*lam(1,1)*tt1)*F(1,4)-4*lam(1,1)*tt10*F(1,3)+15*&
&F(1,1)*lam(1,1)*tt10)*tt21+tt26*tt2+(3*F(1,1)*lam(1,1)-6*F(1,1)*m&
&u(1,1))*tt3*tt4+((tt25-6*F(1,1)*mu(1,1)*tt3)*F(1,4)+(tt14-6*tt5*m&
&u(1,1)*F(1,2))*F(1,3)+tt24+tt23)*F(2,1))*F(2,2)+((-lam(1,1)*tt1*t&
&t16)+4*lam(1,1)*tt10*F(1,4)-3*lam(1,1)*tt22)*F(2,1)**5+(tt20+tt19&
&+tt18)*tt21+(2*mu(1,1)-lam(1,1))*tt1*tt2+(tt17+(tt15+tt14-2*tt5*m&
&u(1,1)*F(1,2))*F(1,4)+tt13+4*F(1,1)*mu(1,1)*tt3*F(1,3)+tt11+tt9)*&
&tt4)/(2*tt7*tt8-6*tt5*F(1,2)*F(2,1)*tt6+6*F(1,1)*tt3*tt4*F(2,2)-2&
&*tt1*tt2)
END 
SUBROUTINE quad_neo_q2l_psi_jac(jac, F, mu, lam) 
IMPLICIT NONE 
REAL(KIND=8) jac(1, 8) 
REAL(KIND=8) F(2, 4) 
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
tt1 = F(1,2)**3
tt2 = F(2,1)**3
tt3 = F(1,2)**2
tt4 = F(2,1)**2
tt5 = F(1,1)**2
tt6 = F(2,2)**2
tt7 = F(1,1)**3
tt8 = F(2,2)**3
tt9 = 2*tt7*tt8-6*tt5*F(1,2)*F(2,1)*tt6+6*F(1,1)*tt3*tt4*F(2,2)-2&
&*tt1*tt2
tt10 = 1/tt9
tt11 = 2*F(1,1)*mu(1,1)*tt3
tt12 = -2*mu(1,1)*F(1,2)*F(1,3)
tt13 = F(1,4)**2
tt14 = 2*F(1,1)*mu(1,1)*tt13
tt15 = -6*tt5*mu(1,1)*F(1,2)
tt16 = -2*mu(1,1)*tt1
tt17 = 6*tt5*mu(1,1)
tt18 = -2*mu(1,1)*F(1,2)
tt19 = 2*mu(1,1)*F(1,4)+tt18
tt20 = F(1,2)**4
tt21 = -16*lam(1,1)*tt1*F(1,4)
tt22 = F(2,1)**4
tt23 = 4*tt7*mu(1,1)
tt24 = -2*mu(1,1)*tt3
tt25 = 2*F(1,1)*mu(1,1)
tt26 = -2*mu(1,1)*F(1,3)
tt27 = tt26+tt25
tt28 = -6*tt5*mu(1,1)
tt29 = F(1,3)**2
tt30 = 12*F(1,1)*lam(1,1)*F(1,2)*F(1,3)
tt31 = F(2,2)**4
tt32 = F(1,1)**4
tt33 = F(2,2)**5
tt34 = -4*F(1,1)*mu(1,1)*tt3
tt35 = 6*tt5*mu(1,1)*F(1,2)
tt36 = 2*mu(1,1)*tt1
tt37 = 3*mu(1,1)*tt3+3*tt5*mu(1,1)
tt38 = -6*mu(1,1)*F(1,2)
tt39 = 3*mu(1,1)
tt40 = 6*lam(1,1)*tt3
tt41 = 2*mu(1,1)*F(1,2)
tt42 = 6*tt5*lam(1,1)
tt43 = 16*tt7*lam(1,1)*F(1,2)
tt44 = -6*tt5*lam(1,1)*F(1,2)*F(1,3)
tt45 = -8*tt7*lam(1,1)*F(1,4)
tt46 = F(2,3)**2
tt47 = -2*mu(1,1)*F(1,4)
tt48 = tt47+tt41
tt49 = tt48*tt2
tt50 = (2*lam(1,1)*tt1*F(1,4)-4*lam(1,1)*tt20)*tt22
tt51 = mu(1,1)*tt3
tt52 = -6*lam(1,1)*tt3
tt53 = -6*tt5*lam(1,1)
tt54 = -mu(1,1)
tt55 = 16*tt7*lam(1,1)*F(1,3)
tt56 = 2*lam(1,1)*tt20*tt2
tt57 = F(2,4)**2
tt58 = 1/tt9**2
tt59 = tt5*mu(1,1)*tt3
tt60 = mu(1,1)*tt20
tt61 = mu(1,1)*tt3*tt29
tt62 = -2*tt5*mu(1,1)*F(1,2)
tt63 = -2*F(1,1)*mu(1,1)*F(1,2)*F(1,3)
tt64 = tt5*mu(1,1)*tt13
tt65 = tt64+(tt63+tt36+tt62)*F(1,4)+tt61+4*F(1,1)*mu(1,1)*tt3*F(1&
&,3)+tt60+tt59
tt66 = 2*mu(1,1)-lam(1,1)
tt67 = -2*mu(1,1)*F(1,2)*F(1,4)
tt68 = mu(1,1)*tt13
tt69 = tt68+tt67+tt51
tt70 = F(1,2)**5
tt71 = (-lam(1,1)*tt1*tt13)+4*lam(1,1)*tt20*F(1,4)-3*lam(1,1)*tt7&
&0
tt72 = F(2,1)**5
tt73 = -2*tt7*mu(1,1)*F(1,2)
tt74 = -2*F(1,1)*mu(1,1)*tt1
tt75 = tt36+tt15
tt76 = tt75*F(1,3)
tt77 = 2*tt7*mu(1,1)
tt78 = -6*F(1,1)*mu(1,1)*tt3
tt79 = tt78+tt77
tt80 = tt79*F(1,4)
tt81 = (tt80+tt76+tt74+tt73)*F(2,1)
tt82 = -6*F(1,1)*mu(1,1)
tt83 = tt82+3*F(1,1)*lam(1,1)
tt84 = tt83*tt3*tt4
tt85 = -2*F(1,1)*mu(1,1)*F(1,2)
tt86 = 2*mu(1,1)*F(1,2)*F(1,3)
tt87 = tt27*F(1,4)
tt88 = tt87+tt86+tt85
tt89 = tt88*tt2
tt90 = -16*F(1,1)*lam(1,1)*tt1
tt91 = 2*lam(1,1)*tt1*F(1,3)
tt92 = 3*F(1,1)*lam(1,1)*tt3*tt13+(tt91+tt90)*F(1,4)-4*lam(1,1)*t&
&t20*F(1,3)+15*F(1,1)*lam(1,1)*tt20
tt93 = tt92*tt22
tt94 = tt32*mu(1,1)
tt95 = -2*F(1,1)*mu(1,1)*tt3
tt96 = 4*tt5*mu(1,1)*F(1,2)
tt97 = tt17-3*tt5*lam(1,1)
tt98 = tt5*mu(1,1)
tt99 = -2*F(1,1)*mu(1,1)*F(1,3)
tt100 = mu(1,1)*tt29
tt101 = tt68+tt67+tt100+tt99+tt51+tt98
tt102 = 24*tt5*lam(1,1)*tt3
tt103 = -6*F(1,1)*lam(1,1)*tt3*F(1,3)
tt104 = (-3*tt5*lam(1,1)*F(1,2)*tt13)+(tt103+tt102)*F(1,4)-lam(1,&
&1)*tt1*tt29+16*F(1,1)*lam(1,1)*tt1*F(1,3)-30*tt5*lam(1,1)*tt1
tt105 = tt104*tt2+tt101*tt4+tt97*F(1,2)*F(2,1)+tt64+(tt63+tt96)*F&
&(1,4)+tt61+(tt95+tt77)*F(1,3)+tt59+tt94
tt106 = -2*tt7*mu(1,1)
tt107 = -16*tt7*lam(1,1)*F(1,2)
tt108 = 6*tt5*lam(1,1)*F(1,2)*F(1,3)
tt109 = tt7*lam(1,1)*tt13+(tt108+tt107)*F(1,4)+3*F(1,1)*lam(1,1)*&
&tt3*tt29-24*tt5*lam(1,1)*tt3*F(1,3)+30*tt7*lam(1,1)*tt3
tt110 = tt109*tt4+tt88*F(2,1)+tt106+tt7*lam(1,1)
tt111 = 4*tt32*lam(1,1)-2*tt7*lam(1,1)*F(1,3)
tt112 = tt111*F(1,4)-3*tt5*lam(1,1)*F(1,2)*tt29+16*tt7*lam(1,1)*F&
&(1,2)*F(1,3)-15*tt32*lam(1,1)*F(1,2)
tt113 = tt112*F(2,1)+tt100+tt99+tt98
tt114 = F(1,1)**5
tt115 = tt7*lam(1,1)*tt29-4*tt32*lam(1,1)*F(1,3)+3*tt114*lam(1,1)&
&
tt116 = -2*tt5*mu(1,1)*tt3
tt117 = -2*mu(1,1)*tt20
tt118 = -2*F(1,1)*mu(1,1)*tt3*F(1,3)
tt119 = tt5*mu(1,1)*F(1,2)-mu(1,1)*tt1
tt120 = tt119*F(1,4)
tt121 = (tt120+tt118+tt117+tt116)*F(2,1)
tt122 = (-2*lam(1,1)*tt1)-mu(1,1)*F(1,2)
tt123 = tt122*F(1,4)+2*mu(1,1)*tt3
tt124 = tt123*tt2
tt125 = 4*lam(1,1)*tt70-4*lam(1,1)*tt20*F(1,4)
tt126 = tt125*tt22
tt127 = 2*tt7*mu(1,1)*F(1,2)
tt128 = 2*F(1,1)*mu(1,1)*tt1
tt129 = -2*mu(1,1)*tt1*F(1,3)
tt130 = tt7*mu(1,1)
tt131 = 3*F(1,1)*mu(1,1)*tt3
tt132 = (tt131+tt130)*F(1,4)
tt133 = -6*F(1,1)*mu(1,1)*F(1,2)
tt134 = 6*F(1,1)*lam(1,1)*tt3+3*F(1,1)*mu(1,1)
tt135 = tt134*F(1,4)+tt133
tt136 = tt135*tt4
tt137 = 2*lam(1,1)*tt20*F(1,3)
tt138 = 14*F(1,1)*lam(1,1)*tt1*F(1,4)+tt137-16*F(1,1)*lam(1,1)*tt&
&20
tt139 = tt138*tt2
tt140 = (tt139+tt136+tt132+tt129+tt128+tt127)*F(2,2)
tt141 = 4*tt5*mu(1,1)
tt142 = mu(1,1)+tt53
tt143 = tt142*F(1,2)*F(1,4)
tt144 = (-18*tt5*lam(1,1)*tt3*F(1,4))-6*F(1,1)*lam(1,1)*tt1*F(1,3&
&)+24*tt5*lam(1,1)*tt1
tt145 = tt144*tt4+(tt143+tt99+tt24+tt141)*F(2,1)
tt146 = tt145*tt6
tt147 = 2*F(1,1)*mu(1,1)*F(1,2)
tt148 = 2*tt7*lam(1,1)
tt149 = F(1,1)*mu(1,1)
tt150 = 10*tt7*lam(1,1)*F(1,2)*F(1,4)+6*tt5*lam(1,1)*tt3*F(1,3)-1&
&6*tt7*lam(1,1)*tt3
tt151 = tt150*F(2,1)+(tt149+tt148)*F(1,4)+tt12+tt147
tt152 = tt151*tt8
tt153 = -2*tt32*lam(1,1)*F(1,4)
tt154 = tt153-2*tt7*lam(1,1)*F(1,2)*F(1,3)+4*tt32*lam(1,1)*F(1,2)&
&
tt155 = tt154*tt31
tt156 = (-3*tt5*lam(1,1)*tt1*F(2,1))+tt51+tt98
tt157 = tt7*lam(1,1)*tt3*tt8+tt156*tt6+3*F(1,1)*lam(1,1)*tt20*tt4&
&*F(2,2)-lam(1,1)*tt70*tt2+tt60+tt59
tt158 = mu(1,1)*tt1+3*tt5*mu(1,1)*F(1,2)
tt159 = tt158*F(1,3)
tt160 = -2*tt7*mu(1,1)*F(1,4)
tt161 = (tt160+tt159+tt128+tt127)*F(2,1)
tt162 = 2*lam(1,1)*tt1+mu(1,1)*F(1,2)
tt163 = (-2*F(1,1)*mu(1,1)*F(1,4))+tt162*F(1,3)+tt147
tt164 = tt163*tt2
tt165 = 2*F(1,1)*lam(1,1)*tt1*F(1,4)+tt137-4*F(1,1)*lam(1,1)*tt20&
&
tt166 = tt165*tt22
tt167 = -2*tt32*mu(1,1)
tt168 = -tt7*mu(1,1)
tt169 = F(1,1)*mu(1,1)*tt3
tt170 = (tt169+tt168)*F(1,3)
tt171 = -2*tt5*mu(1,1)*F(1,2)*F(1,4)
tt172 = -2*tt5*mu(1,1)
tt173 = tt149-6*F(1,1)*lam(1,1)*tt3
tt174 = tt67+tt173*F(1,3)+4*mu(1,1)*tt3+tt172
tt175 = tt174*tt4
tt176 = (-6*tt5*lam(1,1)*tt3*F(1,4))-10*F(1,1)*lam(1,1)*tt1*F(1,3&
&)+16*tt5*lam(1,1)*tt1
tt177 = tt176*tt2
tt178 = (tt177+tt175+tt171+tt170+tt116+tt167)*F(2,2)
tt179 = tt39+tt42
tt180 = tt179*F(1,2)*F(1,3)
tt181 = 6*tt7*lam(1,1)*F(1,2)*F(1,4)+18*tt5*lam(1,1)*tt3*F(1,3)-2&
&4*tt7*lam(1,1)*tt3
tt182 = tt181*tt4+(tt180+tt133)*F(2,1)
tt183 = tt182*tt6
tt184 = -2*tt7*lam(1,1)
tt185 = -F(1,1)*mu(1,1)
tt186 = tt153-14*tt7*lam(1,1)*F(1,2)*F(1,3)+16*tt32*lam(1,1)*F(1,&
&2)
tt187 = tt186*F(2,1)+(tt185+tt184)*F(1,3)+2*tt5*mu(1,1)
tt188 = tt187*tt8
tt189 = 4*tt32*lam(1,1)*F(1,3)-4*tt114*lam(1,1)
tt190 = tt189*tt31
tt191 = (tt24+tt172)*F(2,1)
tt192 = -6*tt5*lam(1,1)*tt1*tt4
tt193 = (-2*tt32*lam(1,1)*F(1,2)*tt8)+6*tt7*lam(1,1)*tt3*F(2,1)*t&
&t6+(tt192+tt191)*F(2,2)+2*F(1,1)*lam(1,1)*tt20*tt2+tt74+tt73
tt194 = tt193*F(2,3)
tt195 = tt51+tt98
tt196 = tt114*lam(1,1)*tt8-3*tt32*lam(1,1)*F(1,2)*F(2,1)*tt6+3*tt&
&7*lam(1,1)*tt3*tt4*F(2,2)-tt5*lam(1,1)*tt1*tt2+tt195*tt4+tt59+tt9&
&4
tt197 = tt196*tt57+(tt194+tt190+tt188+tt183+tt178+tt166+tt164+tt1&
&61)*F(2,4)+tt157*tt46+(tt155+tt152+tt146+tt140+tt126+tt124+tt121)&
&*F(2,3)+tt115*tt33+tt113*tt31+tt110*tt8+tt105*tt6+(tt93+tt89+tt84&
&+tt81)*F(2,2)+tt71*tt72+tt69*tt22+tt66*tt1*tt2+tt65*tt4
tt198 = 2*tt5*mu(1,1)*F(1,2)
tt199 = 4*mu(1,1)*tt1
tt200 = 2*mu(1,1)*F(1,2)*tt29
tt201 = 6*mu(1,1)*tt3
tt202 = -2*F(1,1)*mu(1,1)
tt203 = 2*mu(1,1)*F(1,3)
tt204 = tt203+tt202
tt205 = -4*F(1,1)*mu(1,1)*F(1,2)*F(1,3)
tt206 = -12*F(1,1)*lam(1,1)*F(1,2)*F(1,3)
tt207 = -4*tt5*mu(1,1)*F(1,2)
tt208 = 6*F(1,1)*mu(1,1)*tt3
tt209 = 8*lam(1,1)*tt1*F(1,3)
tt210 = 6*F(1,1)*lam(1,1)*tt3*F(1,4)
tt211 = -18*tt5*lam(1,1)*tt3*tt4
tt212 = 2*mu(1,1)*tt3*F(1,3)
tt213 = -2*F(1,1)*mu(1,1)*F(1,2)*F(1,4)
tt214 = 2*tt5*mu(1,1)*F(1,4)
jac(1,1) = tt10*((5*tt32*lam(1,1)*tt8-12*tt7*lam(1,1)*F(1,2)*F(2,&
&1)*tt6+9*tt5*lam(1,1)*tt3*tt4*F(2,2)-2*F(1,1)*lam(1,1)*tt1*tt2+2*&
&F(1,1)*mu(1,1)*tt4+tt11+tt23)*tt57+(((-8*tt7*lam(1,1)*F(1,2)*tt8)&
&+18*tt5*lam(1,1)*tt3*F(2,1)*tt6+((-12*F(1,1)*lam(1,1)*tt1*tt4)-4*&
&F(1,1)*mu(1,1)*F(2,1))*F(2,2)+tt56+tt16+tt15)*F(2,3)+(tt55-20*tt3&
&2*lam(1,1))*tt31+((tt45-42*tt5*lam(1,1)*F(1,2)*F(1,3)+64*tt7*lam(&
&1,1)*F(1,2))*F(2,1)+(tt54+tt53)*F(1,3)+4*F(1,1)*mu(1,1))*tt8+((18&
&*tt5*lam(1,1)*F(1,2)*F(1,4)+36*F(1,1)*lam(1,1)*tt3*F(1,3)-72*tt5*&
&lam(1,1)*tt3)*tt4+(tt30+tt38)*F(2,1))*tt6+(((-12*F(1,1)*lam(1,1)*&
&tt3*F(1,4))-10*lam(1,1)*tt1*F(1,3)+32*F(1,1)*lam(1,1)*tt1)*tt2+((&
&tt52+mu(1,1))*F(1,3)-4*F(1,1)*mu(1,1))*tt4-4*F(1,1)*mu(1,1)*F(1,2&
&)*F(1,4)+(tt51-3*tt5*mu(1,1))*F(1,3)+tt34-8*tt7*mu(1,1))*F(2,2)+t&
&t50+tt49+((-6*tt5*mu(1,1)*F(1,4))+6*F(1,1)*mu(1,1)*F(1,2)*F(1,3)+&
&tt36+tt35)*F(2,1))*F(2,4)+(3*tt5*lam(1,1)*tt3*tt8+(tt25-6*F(1,1)*&
&lam(1,1)*tt1*F(2,1))*tt6+3*lam(1,1)*tt20*tt4*F(2,2)+tt11)*tt46+((&
&tt45+tt44+tt43)*tt31+((30*tt5*lam(1,1)*F(1,2)*F(1,4)+12*F(1,1)*la&
&m(1,1)*tt3*F(1,3)-48*tt5*lam(1,1)*tt3)*F(2,1)+(mu(1,1)+tt42)*F(1,&
&4)+tt41)*tt8+(((-36*F(1,1)*lam(1,1)*tt3*F(1,4))-6*lam(1,1)*tt1*F(&
&1,3)+48*F(1,1)*lam(1,1)*tt1)*tt4+((-12*F(1,1)*lam(1,1)*F(1,2)*F(1&
&,4))+tt26+8*F(1,1)*mu(1,1))*F(2,1))*tt6+((14*lam(1,1)*tt1*F(1,4)-&
&16*lam(1,1)*tt20)*tt2+((tt40+tt39)*F(1,4)+tt38)*tt4+tt37*F(1,4)+t&
&t36+tt35)*F(2,2)+(2*F(1,1)*mu(1,1)*F(1,2)*F(1,4)-2*mu(1,1)*tt3*F(&
&1,3)+tt34)*F(2,1))*F(2,3)+(3*tt5*lam(1,1)*tt29-16*tt7*lam(1,1)*F(&
&1,3)+15*tt32*lam(1,1))*tt33+(((16*tt7*lam(1,1)-6*tt5*lam(1,1)*F(1&
&,3))*F(1,4)-6*F(1,1)*lam(1,1)*F(1,2)*tt29+48*tt5*lam(1,1)*F(1,2)*&
&F(1,3)-60*tt7*lam(1,1)*F(1,2))*F(2,1)+tt26+tt25)*tt31+((3*tt5*lam&
&(1,1)*tt13+(tt30-48*tt5*lam(1,1)*F(1,2))*F(1,4)+3*lam(1,1)*tt3*tt&
&29-48*F(1,1)*lam(1,1)*tt3*F(1,3)+90*tt5*lam(1,1)*tt3)*tt4+tt19*F(&
&2,1)+tt28+3*tt5*lam(1,1))*tt8+(((-6*F(1,1)*lam(1,1)*F(1,2)*tt13)+&
&(48*F(1,1)*lam(1,1)*tt3-6*lam(1,1)*tt3*F(1,3))*F(1,4)+16*lam(1,1)&
&*tt1*F(1,3)-60*F(1,1)*lam(1,1)*tt1)*tt2+tt27*tt4+(12*F(1,1)*mu(1,&
&1)-6*F(1,1)*lam(1,1))*F(1,2)*F(2,1)+tt14+(tt12+8*F(1,1)*mu(1,1)*F&
&(1,2))*F(1,4)+(tt24+tt17)*F(1,3)+tt11+tt23)*tt6+((3*lam(1,1)*tt3*&
&tt13+tt21+15*lam(1,1)*tt20)*tt22+tt19*tt2+(3*lam(1,1)-6*mu(1,1))*&
&tt3*tt4+((tt17-6*mu(1,1)*tt3)*F(1,4)-12*F(1,1)*mu(1,1)*F(1,2)*F(1&
&,3)+tt16+tt15)*F(2,1))*F(2,2)+(tt14+(tt12-4*F(1,1)*mu(1,1)*F(1,2)&
&)*F(1,4)+4*mu(1,1)*tt3*F(1,3)+tt11)*tt4)-(6*tt5*tt8-12*F(1,1)*F(1&
&,2)*F(2,1)*tt6+6*tt3*tt4*F(2,2))*tt58*tt197
jac(1,2) = tt10*(((-3*tt32*lam(1,1)*F(1,2)*tt6)+6*tt7*lam(1,1)*tt&
&3*F(2,1)*F(2,2)-3*tt5*lam(1,1)*tt1*tt4+2*tt195*F(2,1))*tt57+((6*t&
&t7*lam(1,1)*tt3*tt6+((-12*tt5*lam(1,1)*tt1*F(2,1))+tt24+tt172)*F(&
&2,2)+6*F(1,1)*lam(1,1)*tt20*tt4)*F(2,3)+tt186*tt8+(2*tt181*F(2,1)&
&+tt180+tt133)*tt6+(3*tt176*tt4+2*tt174*F(2,1))*F(2,2)+4*tt165*tt2&
&+3*tt163*tt4+tt160+tt159+tt128+tt127)*F(2,4)+((-3*tt5*lam(1,1)*tt&
&1*tt6)+6*F(1,1)*lam(1,1)*tt20*F(2,1)*F(2,2)-3*lam(1,1)*tt70*tt4)*&
&tt46+(tt150*tt8+(2*tt144*F(2,1)+tt143+tt99+tt24+tt141)*tt6+(3*tt1&
&38*tt4+2*tt135*F(2,1))*F(2,2)+4*tt125*tt2+3*tt123*tt4+tt120+tt118&
&+tt117+tt116)*F(2,3)+tt112*tt31+(2*tt109*F(2,1)+tt87+tt86+tt85)*t&
&t8+(3*tt104*tt4+2*tt101*F(2,1)+tt97*F(1,2))*tt6+(4*tt92*tt2+3*tt8&
&8*tt4+2*tt83*tt3*F(2,1)+tt80+tt76+tt74+tt73)*F(2,2)+5*tt71*tt22+4&
&*tt69*tt2+3*tt66*tt1*tt4+2*tt65*F(2,1))-((-6*tt5*F(1,2)*tt6)+12*F&
&(1,1)*tt3*F(2,1)*F(2,2)-6*tt1*tt4)*tt58*tt197
jac(1,3) = tt10*(((-3*tt32*lam(1,1)*F(2,1)*tt6)+6*tt7*lam(1,1)*F(&
&1,2)*tt4*F(2,2)-3*tt5*lam(1,1)*tt3*tt2+2*mu(1,1)*F(1,2)*tt4+tt198&
&)*tt57+(((-2*tt32*lam(1,1)*tt8)+12*tt7*lam(1,1)*F(1,2)*F(2,1)*tt6&
&+(tt211-4*mu(1,1)*F(1,2)*F(2,1))*F(2,2)+8*F(1,1)*lam(1,1)*tt1*tt2&
&+tt78+tt106)*F(2,3)+(16*tt32*lam(1,1)-14*tt7*lam(1,1)*F(1,3))*F(2&
&,1)*tt8+((6*tt7*lam(1,1)*F(1,4)+36*tt5*lam(1,1)*F(1,2)*F(1,3)-48*&
&tt7*lam(1,1)*F(1,2))*tt4+(tt179*F(1,3)+tt82)*F(2,1))*tt6+(((-12*t&
&t5*lam(1,1)*F(1,2)*F(1,4))-30*F(1,1)*lam(1,1)*tt3*F(1,3)+48*tt5*l&
&am(1,1)*tt3)*tt2+(tt47+tt206+8*mu(1,1)*F(1,2))*tt4-2*tt5*mu(1,1)*&
&F(1,4)+2*F(1,1)*mu(1,1)*F(1,2)*F(1,3)+tt207)*F(2,2)+(tt210+tt209+&
&tt90)*tt22+((tt40+mu(1,1))*F(1,3)+tt25)*tt2+(tt37*F(1,3)+tt208+tt&
&77)*F(2,1))*F(2,4)+(2*tt7*lam(1,1)*F(1,2)*tt8+(tt41-9*tt5*lam(1,1&
&)*tt3*F(2,1))*tt6+12*F(1,1)*lam(1,1)*tt1*tt4*F(2,2)-5*lam(1,1)*tt&
&20*tt2+tt199+tt198)*tt46+(tt111*tt31+((10*tt7*lam(1,1)*F(1,4)+12*&
&tt5*lam(1,1)*F(1,2)*F(1,3)-32*tt7*lam(1,1)*F(1,2))*F(2,1)+tt26+tt&
&25)*tt8+(((-36*tt5*lam(1,1)*F(1,2)*F(1,4))-18*F(1,1)*lam(1,1)*tt3&
&*F(1,3)+72*tt5*lam(1,1)*tt3)*tt4+(tt142*F(1,4)-4*mu(1,1)*F(1,2))*&
&F(2,1))*tt6+((42*F(1,1)*lam(1,1)*tt3*F(1,4)+tt209-64*F(1,1)*lam(1&
&,1)*tt1)*tt2+(12*F(1,1)*lam(1,1)*F(1,2)*F(1,4)+tt82)*tt4+6*F(1,1)&
&*mu(1,1)*F(1,2)*F(1,4)-6*mu(1,1)*tt3*F(1,3)+tt208+tt77)*F(2,2)+(t&
&t21+20*lam(1,1)*tt20)*tt22+((tt52+tt54)*F(1,4)+4*mu(1,1)*F(1,2))*&
&tt2+((tt98-3*mu(1,1)*tt3)*F(1,4)+tt205-8*mu(1,1)*tt1+tt207)*F(2,1&
&))*F(2,3)+((-3*tt5*lam(1,1)*tt29)+tt55-15*tt32*lam(1,1))*F(2,1)*t&
&t31+(((6*tt5*lam(1,1)*F(1,3)-16*tt7*lam(1,1))*F(1,4)+6*F(1,1)*lam&
&(1,1)*F(1,2)*tt29-48*tt5*lam(1,1)*F(1,2)*F(1,3)+60*tt7*lam(1,1)*F&
&(1,2))*tt4+tt204*F(2,1))*tt8+(((-3*tt5*lam(1,1)*tt13)+(tt206+48*t&
&t5*lam(1,1)*F(1,2))*F(1,4)-3*lam(1,1)*tt3*tt29+48*F(1,1)*lam(1,1)&
&*tt3*F(1,3)-90*tt5*lam(1,1)*tt3)*tt2+tt48*tt4+tt97*F(2,1)+(tt99+t&
&t141)*F(1,4)+tt200+tt205+tt198)*tt6+((6*F(1,1)*lam(1,1)*F(1,2)*tt&
&13+(6*lam(1,1)*tt3*F(1,3)-48*F(1,1)*lam(1,1)*tt3)*F(1,4)-16*lam(1&
&,1)*tt1*F(1,3)+60*F(1,1)*lam(1,1)*tt1)*tt22+tt204*tt2+2*tt83*F(1,&
&2)*tt4+((-12*F(1,1)*mu(1,1)*F(1,2)*F(1,4))+(tt201+tt28)*F(1,3)+tt&
&78+tt106)*F(2,1))*F(2,2)+((-3*lam(1,1)*tt3*tt13)+16*lam(1,1)*tt1*&
&F(1,4)-15*lam(1,1)*tt20)*tt72+tt48*tt22+3*tt66*tt3*tt2+((tt99+tt2&
&01+tt172)*F(1,4)+tt200+8*F(1,1)*mu(1,1)*F(1,2)*F(1,3)+tt199+tt198&
&)*tt4)-((-6*tt5*F(2,1)*tt6)+12*F(1,1)*F(1,2)*tt4*F(2,2)-6*tt3*tt2&
&)*tt58*tt197
jac(1,4) = tt10*((3*tt114*lam(1,1)*tt6-6*tt32*lam(1,1)*F(1,2)*F(2&
&,1)*F(2,2)+3*tt7*lam(1,1)*tt3*tt4)*tt57+(((-6*tt32*lam(1,1)*F(1,2&
&)*tt6)+12*tt7*lam(1,1)*tt3*F(2,1)*F(2,2)+tt192+tt191)*F(2,3)+4*tt&
&189*tt8+3*tt187*tt6+2*tt182*F(2,2)+tt177+tt175+tt171+tt170+tt116+&
&tt167)*F(2,4)+(3*tt7*lam(1,1)*tt3*tt6+2*tt156*F(2,2)+3*F(1,1)*lam&
&(1,1)*tt20*tt4)*tt46+(4*tt154*tt8+3*tt151*tt6+2*tt145*F(2,2)+tt13&
&9+tt136+tt132+tt129+tt128+tt127)*F(2,3)+5*tt115*tt31+4*tt113*tt8+&
&3*tt110*tt6+2*tt105*F(2,2)+tt93+tt89+tt84+tt81)-(6*tt7*tt6-12*tt5&
&*F(1,2)*F(2,1)*F(2,2)+6*F(1,1)*tt3*tt4)*tt58*tt197
jac(1,5) = tt10*((4*tt32*lam(1,1)*tt31+((-14*tt7*lam(1,1)*F(1,2)*&
&F(2,1))+tt185+tt184)*tt8+(18*tt5*lam(1,1)*tt3*tt4+tt179*F(1,2)*F(&
&2,1))*tt6+((-10*F(1,1)*lam(1,1)*tt1*tt2)+tt173*tt4+tt169+tt168)*F&
&(2,2)+2*lam(1,1)*tt20*tt22+tt162*tt2+tt158*F(2,1))*F(2,4)+((-2*tt&
&7*lam(1,1)*F(1,2)*tt31)+(6*tt5*lam(1,1)*tt3*F(2,1)+tt18)*tt8+((-6&
&*F(1,1)*lam(1,1)*tt1*tt4)-2*F(1,1)*mu(1,1)*F(2,1))*tt6+(tt56+tt16&
&)*F(2,2)-2*F(1,1)*mu(1,1)*tt3*F(2,1))*F(2,3)+(2*tt7*lam(1,1)*F(1,&
&3)-4*tt32*lam(1,1))*tt33+(((-2*tt7*lam(1,1)*F(1,4))+tt44+tt43)*F(&
&2,1)+tt203+tt202)*tt31+((6*tt5*lam(1,1)*F(1,2)*F(1,4)+6*F(1,1)*la&
&m(1,1)*tt3*F(1,3)-24*tt5*lam(1,1)*tt3)*tt4+tt48*F(2,1))*tt8+(((-6&
&*F(1,1)*lam(1,1)*tt3*F(1,4))-2*lam(1,1)*tt1*F(1,3)+16*F(1,1)*lam(&
&1,1)*tt1)*tt2+tt204*tt4+tt213+tt212+tt95+tt77)*tt6+(tt50+tt49+tt7&
&5*F(2,1))*F(2,2)+(tt213+tt212+4*F(1,1)*mu(1,1)*tt3)*tt4)
jac(1,6) = tt10*(tt193*F(2,4)+2*tt157*F(2,3)+tt155+tt152+tt146+tt&
&140+tt126+tt124+tt121)
jac(1,7) = tt10*(((-2*tt32*lam(1,1)*F(2,1)*tt8)+6*tt7*lam(1,1)*F(&
&1,2)*tt4*tt6+((-6*tt5*lam(1,1)*tt3*tt2)-2*mu(1,1)*F(1,2)*tt4+tt62&
&)*F(2,2)+2*F(1,1)*lam(1,1)*tt1*tt22-2*F(1,1)*mu(1,1)*tt2-2*tt7*mu&
&(1,1)*F(2,1))*F(2,4)+((-2*tt32*lam(1,1)*tt31)+(10*tt7*lam(1,1)*F(&
&1,2)*F(2,1)+tt149+tt148)*tt8+(tt211+tt142*F(1,2)*F(2,1))*tt6+(14*&
&F(1,1)*lam(1,1)*tt1*tt2+tt134*tt4+tt131+tt130)*F(2,2)-4*lam(1,1)*&
&tt20*tt22+tt122*tt2+tt119*F(2,1))*F(2,3)+tt111*F(2,1)*tt31+((2*tt&
&7*lam(1,1)*F(1,4)+tt108+tt107)*tt4+tt27*F(2,1))*tt8+(((-6*tt5*lam&
&(1,1)*F(1,2)*F(1,4))+tt103+tt102)*tt2+tt19*tt4+tt214+tt63+tt96)*t&
&t6+((tt210+tt91+tt90)*tt22+tt27*tt2+tt79*F(2,1))*F(2,2)+(4*lam(1,&
&1)*tt20-2*lam(1,1)*tt1*F(1,4))*tt72+tt19*tt22+(tt214+tt63+tt36+tt&
&62)*tt4)
jac(1,8) = tt10*(2*tt196*F(2,4)+tt194+tt190+tt188+tt183+tt178+tt1&
&66+tt164+tt161)
END 
SUBROUTINE quad_neo_q2l_psi_hes(hes, F, mu, lam) 
IMPLICIT NONE 
REAL(KIND=8) hes(8, 8) 
REAL(KIND=8) F(2, 4) 
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
tt1 = F(1,2)**3
tt2 = F(2,1)**3
tt3 = F(1,2)**2
tt4 = F(2,1)**2
tt5 = F(1,1)**2
tt6 = F(2,2)**2
tt7 = F(1,1)**3
tt8 = F(2,2)**3
tt9 = 2*tt7*tt8-6*tt5*F(1,2)*F(2,1)*tt6+6*F(1,1)*tt3*tt4*F(2,2)-2&
&*tt1*tt2
tt10 = 1/tt9
tt11 = 2*mu(1,1)*tt3
tt12 = -4*mu(1,1)*F(1,2)*F(1,4)
tt13 = F(1,4)**2
tt14 = 2*mu(1,1)*tt13
tt15 = -12*F(1,1)*mu(1,1)*F(1,2)
tt16 = -12*mu(1,1)*F(1,2)*F(1,3)
tt17 = 12*tt5*mu(1,1)
tt18 = 2*mu(1,1)*tt4
tt19 = (-6*lam(1,1)*F(1,2)*tt13)+48*lam(1,1)*tt3*F(1,4)-60*lam(1,&
&1)*tt1
tt20 = 6*F(1,1)*lam(1,1)*tt13+(12*lam(1,1)*F(1,2)*F(1,3)-96*F(1,1&
&)*lam(1,1)*F(1,2))*F(1,4)-48*lam(1,1)*tt3*F(1,3)+180*F(1,1)*lam(1&
&,1)*tt3
tt21 = 2*mu(1,1)
tt22 = F(1,3)**2
tt23 = -12*F(1,1)*lam(1,1)*F(1,3)
tt24 = (tt23+48*tt5*lam(1,1))*F(1,4)-6*lam(1,1)*F(1,2)*tt22+96*F(&
&1,1)*lam(1,1)*F(1,2)*F(1,3)-180*tt5*lam(1,1)*F(1,2)
tt25 = F(2,2)**4
tt26 = 6*F(1,1)*lam(1,1)*tt22-48*tt5*lam(1,1)*F(1,3)+60*tt7*lam(1&
&,1)
tt27 = F(2,2)**5
tt28 = -4*mu(1,1)*tt3
tt29 = 12*F(1,1)*mu(1,1)*F(1,2)
tt30 = 6*F(1,1)*mu(1,1)*F(1,4)
tt31 = 8*mu(1,1)
tt32 = -12*lam(1,1)*F(1,2)*F(1,4)
tt33 = 48*tt5*lam(1,1)*F(1,2)
tt34 = -12*F(1,1)*lam(1,1)*F(1,2)*F(1,3)
tt35 = -24*tt5*lam(1,1)*F(1,4)
tt36 = F(2,3)**2
tt37 = 6*mu(1,1)*F(1,2)*F(1,3)
tt38 = -12*F(1,1)*mu(1,1)*F(1,4)
tt39 = 4*mu(1,1)
tt40 = 48*tt5*lam(1,1)*F(1,3)
tt41 = -4*mu(1,1)*F(2,1)
tt42 = -2*lam(1,1)*tt1*tt2
tt43 = F(2,4)**2
tt44 = 6*tt5*tt8-12*F(1,1)*F(1,2)*F(2,1)*tt6+6*tt3*tt4*F(2,2)
tt45 = 1/tt9**2
tt46 = 2*F(1,1)*mu(1,1)*tt3
tt47 = -4*F(1,1)*mu(1,1)*F(1,2)
tt48 = -2*mu(1,1)*F(1,2)*F(1,3)
tt49 = 2*F(1,1)*mu(1,1)*tt13
tt50 = tt49+(tt48+tt47)*F(1,4)+4*mu(1,1)*tt3*F(1,3)+tt46
tt51 = -6*tt5*mu(1,1)*F(1,2)
tt52 = -2*mu(1,1)*tt1
tt53 = -12*F(1,1)*mu(1,1)*F(1,2)*F(1,3)
tt54 = 6*tt5*mu(1,1)
tt55 = -6*mu(1,1)*tt3
tt56 = tt55+tt54
tt57 = tt56*F(1,4)
tt58 = (tt57+tt53+tt52+tt51)*F(2,1)
tt59 = -6*mu(1,1)
tt60 = tt59+3*lam(1,1)
tt61 = tt60*tt3*tt4
tt62 = -2*mu(1,1)*F(1,2)
tt63 = 2*mu(1,1)*F(1,4)
tt64 = tt63+tt62
tt65 = tt64*tt2
tt66 = F(1,2)**4
tt67 = -16*lam(1,1)*tt1*F(1,4)
tt68 = 3*lam(1,1)*tt3*tt13+tt67+15*lam(1,1)*tt66
tt69 = F(2,1)**4
tt70 = tt68*tt69
tt71 = 4*tt7*mu(1,1)
tt72 = -2*mu(1,1)*tt3
tt73 = 8*F(1,1)*mu(1,1)*F(1,2)
tt74 = 12*F(1,1)*mu(1,1)-6*F(1,1)*lam(1,1)
tt75 = 2*F(1,1)*mu(1,1)
tt76 = -2*mu(1,1)*F(1,3)
tt77 = tt76+tt75
tt78 = 48*F(1,1)*lam(1,1)*tt3
tt79 = -6*lam(1,1)*tt3*F(1,3)
tt80 = (-6*F(1,1)*lam(1,1)*F(1,2)*tt13)+(tt79+tt78)*F(1,4)+16*lam&
&(1,1)*tt1*F(1,3)-60*F(1,1)*lam(1,1)*tt1
tt81 = tt80*tt2+tt77*tt4+tt74*F(1,2)*F(2,1)+tt49+(tt48+tt73)*F(1,&
&4)+(tt72+tt54)*F(1,3)+tt46+tt71
tt82 = -6*tt5*mu(1,1)
tt83 = -48*tt5*lam(1,1)*F(1,2)
tt84 = 12*F(1,1)*lam(1,1)*F(1,2)*F(1,3)
tt85 = 3*tt5*lam(1,1)*tt13+(tt84+tt83)*F(1,4)+3*lam(1,1)*tt3*tt22&
&-48*F(1,1)*lam(1,1)*tt3*F(1,3)+90*tt5*lam(1,1)*tt3
tt86 = tt85*tt4+tt64*F(2,1)+tt82+3*tt5*lam(1,1)
tt87 = 16*tt7*lam(1,1)-6*tt5*lam(1,1)*F(1,3)
tt88 = tt87*F(1,4)-6*F(1,1)*lam(1,1)*F(1,2)*tt22+48*tt5*lam(1,1)*&
&F(1,2)*F(1,3)-60*tt7*lam(1,1)*F(1,2)
tt89 = tt88*F(2,1)+tt76+tt75
tt90 = F(1,1)**4
tt91 = 3*tt5*lam(1,1)*tt22-16*tt7*lam(1,1)*F(1,3)+15*tt90*lam(1,1&
&)
tt92 = -4*F(1,1)*mu(1,1)*tt3
tt93 = -2*mu(1,1)*tt3*F(1,3)
tt94 = 2*F(1,1)*mu(1,1)*F(1,2)*F(1,4)
tt95 = (tt94+tt93+tt92)*F(2,1)
tt96 = 6*tt5*mu(1,1)*F(1,2)
tt97 = 2*mu(1,1)*tt1
tt98 = 3*tt5*mu(1,1)
tt99 = 3*mu(1,1)*tt3
tt100 = tt99+tt98
tt101 = tt100*F(1,4)
tt102 = -6*mu(1,1)*F(1,2)
tt103 = 3*mu(1,1)
tt104 = 6*lam(1,1)*tt3
tt105 = tt104+tt103
tt106 = tt105*F(1,4)+tt102
tt107 = tt106*tt4
tt108 = 14*lam(1,1)*tt1*F(1,4)-16*lam(1,1)*tt66
tt109 = tt108*tt2
tt110 = (tt109+tt107+tt101+tt97+tt96)*F(2,2)
tt111 = 8*F(1,1)*mu(1,1)
tt112 = -12*F(1,1)*lam(1,1)*F(1,2)*F(1,4)
tt113 = (-36*F(1,1)*lam(1,1)*tt3*F(1,4))-6*lam(1,1)*tt1*F(1,3)+48&
&*F(1,1)*lam(1,1)*tt1
tt114 = tt113*tt4+(tt112+tt76+tt111)*F(2,1)
tt115 = tt114*tt6
tt116 = 2*mu(1,1)*F(1,2)
tt117 = 6*tt5*lam(1,1)
tt118 = 30*tt5*lam(1,1)*F(1,2)*F(1,4)+12*F(1,1)*lam(1,1)*tt3*F(1,&
&3)-48*tt5*lam(1,1)*tt3
tt119 = tt118*F(2,1)+(mu(1,1)+tt117)*F(1,4)+tt116
tt120 = tt119*tt8
tt121 = 16*tt7*lam(1,1)*F(1,2)
tt122 = -6*tt5*lam(1,1)*F(1,2)*F(1,3)
tt123 = -8*tt7*lam(1,1)*F(1,4)
tt124 = tt123+tt122+tt121
tt125 = tt124*tt25
tt126 = tt75-6*F(1,1)*lam(1,1)*tt1*F(2,1)
tt127 = 3*tt5*lam(1,1)*tt3*tt8+tt126*tt6+3*lam(1,1)*tt66*tt4*F(2,&
&2)+tt46
tt128 = 6*F(1,1)*mu(1,1)*F(1,2)*F(1,3)
tt129 = -6*tt5*mu(1,1)*F(1,4)
tt130 = (tt129+tt128+tt97+tt96)*F(2,1)
tt131 = -2*mu(1,1)*F(1,4)
tt132 = tt131+tt116
tt133 = tt132*tt2
tt134 = 2*lam(1,1)*tt1*F(1,4)-4*lam(1,1)*tt66
tt135 = tt134*tt69
tt136 = -8*tt7*mu(1,1)
tt137 = -3*tt5*mu(1,1)
tt138 = mu(1,1)*tt3
tt139 = (tt138+tt137)*F(1,3)
tt140 = -4*F(1,1)*mu(1,1)*F(1,2)*F(1,4)
tt141 = -4*F(1,1)*mu(1,1)
tt142 = -6*lam(1,1)*tt3
tt143 = tt142+mu(1,1)
tt144 = tt143*F(1,3)+tt141
tt145 = tt144*tt4
tt146 = (-12*F(1,1)*lam(1,1)*tt3*F(1,4))-10*lam(1,1)*tt1*F(1,3)+3&
&2*F(1,1)*lam(1,1)*tt1
tt147 = tt146*tt2
tt148 = (tt147+tt145+tt140+tt139+tt92+tt136)*F(2,2)
tt149 = 18*tt5*lam(1,1)*F(1,2)*F(1,4)+36*F(1,1)*lam(1,1)*tt3*F(1,&
&3)-72*tt5*lam(1,1)*tt3
tt150 = tt149*tt4+(tt84+tt102)*F(2,1)
tt151 = tt150*tt6
tt152 = -6*tt5*lam(1,1)
tt153 = -mu(1,1)
tt154 = tt123-42*tt5*lam(1,1)*F(1,2)*F(1,3)+64*tt7*lam(1,1)*F(1,2&
&)
tt155 = tt154*F(2,1)+(tt153+tt152)*F(1,3)+4*F(1,1)*mu(1,1)
tt156 = tt155*tt8
tt157 = 16*tt7*lam(1,1)*F(1,3)
tt158 = tt157-20*tt90*lam(1,1)
tt159 = tt158*tt25
tt160 = 2*lam(1,1)*tt66*tt2
tt161 = -4*F(1,1)*mu(1,1)*F(2,1)
tt162 = -12*F(1,1)*lam(1,1)*tt1*tt4
tt163 = -8*tt7*lam(1,1)*F(1,2)*tt8
tt164 = tt163+18*tt5*lam(1,1)*tt3*F(2,1)*tt6+(tt162+tt161)*F(2,2)&
&+tt160+tt52+tt51
tt165 = tt164*F(2,3)
tt166 = 5*tt90*lam(1,1)*tt8-12*tt7*lam(1,1)*F(1,2)*F(2,1)*tt6+9*t&
&t5*lam(1,1)*tt3*tt4*F(2,2)-2*F(1,1)*lam(1,1)*tt1*tt2+2*F(1,1)*mu(&
&1,1)*tt4+tt46+tt71
tt167 = tt166*tt43+(tt165+tt159+tt156+tt151+tt148+tt135+tt133+tt1&
&30)*F(2,4)+tt127*tt36+(tt125+tt120+tt115+tt110+tt95)*F(2,3)+tt91*&
&tt27+tt89*tt25+tt86*tt8+tt81*tt6+(tt70+tt65+tt61+tt58)*F(2,2)+tt5&
&0*tt4
tt168 = 1/tt9**3
tt169 = tt5*mu(1,1)*tt3
tt170 = mu(1,1)*tt66
tt171 = mu(1,1)*tt3*tt22
tt172 = -2*tt5*mu(1,1)*F(1,2)
tt173 = -2*F(1,1)*mu(1,1)*F(1,2)*F(1,3)
tt174 = tt5*mu(1,1)*tt13
tt175 = tt174+(tt173+tt97+tt172)*F(1,4)+tt171+4*F(1,1)*mu(1,1)*tt&
&3*F(1,3)+tt170+tt169
tt176 = tt21-lam(1,1)
tt177 = -2*mu(1,1)*F(1,2)*F(1,4)
tt178 = mu(1,1)*tt13
tt179 = tt178+tt177+tt138
tt180 = F(1,2)**5
tt181 = (-lam(1,1)*tt1*tt13)+4*lam(1,1)*tt66*F(1,4)-3*lam(1,1)*tt&
&180
tt182 = F(2,1)**5
tt183 = -2*tt7*mu(1,1)*F(1,2)
tt184 = -2*F(1,1)*mu(1,1)*tt1
tt185 = tt97+tt51
tt186 = tt185*F(1,3)
tt187 = 2*tt7*mu(1,1)
tt188 = -6*F(1,1)*mu(1,1)*tt3
tt189 = tt188+tt187
tt190 = tt189*F(1,4)
tt191 = (tt190+tt186+tt184+tt183)*F(2,1)
tt192 = -6*F(1,1)*mu(1,1)
tt193 = tt192+3*F(1,1)*lam(1,1)
tt194 = tt193*tt3*tt4
tt195 = -2*F(1,1)*mu(1,1)*F(1,2)
tt196 = 2*mu(1,1)*F(1,2)*F(1,3)
tt197 = tt77*F(1,4)
tt198 = tt197+tt196+tt195
tt199 = tt198*tt2
tt200 = -16*F(1,1)*lam(1,1)*tt1
tt201 = 2*lam(1,1)*tt1*F(1,3)
tt202 = 3*F(1,1)*lam(1,1)*tt3*tt13+(tt201+tt200)*F(1,4)-4*lam(1,1&
&)*tt66*F(1,3)+15*F(1,1)*lam(1,1)*tt66
tt203 = tt202*tt69
tt204 = tt90*mu(1,1)
tt205 = -2*F(1,1)*mu(1,1)*tt3
tt206 = 4*tt5*mu(1,1)*F(1,2)
tt207 = -3*tt5*lam(1,1)
tt208 = tt54+tt207
tt209 = tt5*mu(1,1)
tt210 = -2*F(1,1)*mu(1,1)*F(1,3)
tt211 = mu(1,1)*tt22
tt212 = tt178+tt177+tt211+tt210+tt138+tt209
tt213 = 24*tt5*lam(1,1)*tt3
tt214 = -6*F(1,1)*lam(1,1)*tt3*F(1,3)
tt215 = (-3*tt5*lam(1,1)*F(1,2)*tt13)+(tt214+tt213)*F(1,4)-lam(1,&
&1)*tt1*tt22+16*F(1,1)*lam(1,1)*tt1*F(1,3)-30*tt5*lam(1,1)*tt1
tt216 = tt215*tt2+tt212*tt4+tt208*F(1,2)*F(2,1)+tt174+(tt173+tt20&
&6)*F(1,4)+tt171+(tt205+tt187)*F(1,3)+tt169+tt204
tt217 = -2*tt7*mu(1,1)
tt218 = -16*tt7*lam(1,1)*F(1,2)
tt219 = 6*tt5*lam(1,1)*F(1,2)*F(1,3)
tt220 = tt7*lam(1,1)*tt13+(tt219+tt218)*F(1,4)+3*F(1,1)*lam(1,1)*&
&tt3*tt22-24*tt5*lam(1,1)*tt3*F(1,3)+30*tt7*lam(1,1)*tt3
tt221 = tt220*tt4+tt198*F(2,1)+tt217+tt7*lam(1,1)
tt222 = 4*tt90*lam(1,1)-2*tt7*lam(1,1)*F(1,3)
tt223 = tt222*F(1,4)-3*tt5*lam(1,1)*F(1,2)*tt22+16*tt7*lam(1,1)*F&
&(1,2)*F(1,3)-15*tt90*lam(1,1)*F(1,2)
tt224 = tt223*F(2,1)+tt211+tt210+tt209
tt225 = F(1,1)**5
tt226 = tt7*lam(1,1)*tt22-4*tt90*lam(1,1)*F(1,3)+3*tt225*lam(1,1)&
&
tt227 = -2*tt5*mu(1,1)*tt3
tt228 = -2*mu(1,1)*tt66
tt229 = -2*F(1,1)*mu(1,1)*tt3*F(1,3)
tt230 = tt5*mu(1,1)*F(1,2)
tt231 = -mu(1,1)*tt1
tt232 = tt231+tt230
tt233 = tt232*F(1,4)
tt234 = (tt233+tt229+tt228+tt227)*F(2,1)
tt235 = (-2*lam(1,1)*tt1)-mu(1,1)*F(1,2)
tt236 = tt235*F(1,4)+tt11
tt237 = tt236*tt2
tt238 = 4*lam(1,1)*tt180-4*lam(1,1)*tt66*F(1,4)
tt239 = tt238*tt69
tt240 = 2*tt7*mu(1,1)*F(1,2)
tt241 = 2*F(1,1)*mu(1,1)*tt1
tt242 = -2*mu(1,1)*tt1*F(1,3)
tt243 = tt7*mu(1,1)
tt244 = 3*F(1,1)*mu(1,1)*tt3
tt245 = (tt244+tt243)*F(1,4)
tt246 = -6*F(1,1)*mu(1,1)*F(1,2)
tt247 = 6*F(1,1)*lam(1,1)*tt3+3*F(1,1)*mu(1,1)
tt248 = tt247*F(1,4)+tt246
tt249 = tt248*tt4
tt250 = 2*lam(1,1)*tt66*F(1,3)
tt251 = 14*F(1,1)*lam(1,1)*tt1*F(1,4)+tt250-16*F(1,1)*lam(1,1)*tt&
&66
tt252 = tt251*tt2
tt253 = (tt252+tt249+tt245+tt242+tt241+tt240)*F(2,2)
tt254 = 4*tt5*mu(1,1)
tt255 = mu(1,1)+tt152
tt256 = tt255*F(1,2)*F(1,4)
tt257 = (-18*tt5*lam(1,1)*tt3*F(1,4))-6*F(1,1)*lam(1,1)*tt1*F(1,3&
&)+24*tt5*lam(1,1)*tt1
tt258 = tt257*tt4+(tt256+tt210+tt72+tt254)*F(2,1)
tt259 = tt258*tt6
tt260 = 2*F(1,1)*mu(1,1)*F(1,2)
tt261 = 2*tt7*lam(1,1)
tt262 = F(1,1)*mu(1,1)
tt263 = 10*tt7*lam(1,1)*F(1,2)*F(1,4)+6*tt5*lam(1,1)*tt3*F(1,3)-1&
&6*tt7*lam(1,1)*tt3
tt264 = tt263*F(2,1)+(tt262+tt261)*F(1,4)+tt48+tt260
tt265 = tt264*tt8
tt266 = -2*tt90*lam(1,1)*F(1,4)
tt267 = tt266-2*tt7*lam(1,1)*F(1,2)*F(1,3)+4*tt90*lam(1,1)*F(1,2)&
&
tt268 = tt267*tt25
tt269 = (-3*tt5*lam(1,1)*tt1*F(2,1))+tt138+tt209
tt270 = tt7*lam(1,1)*tt3*tt8+tt269*tt6+3*F(1,1)*lam(1,1)*tt66*tt4&
&*F(2,2)-lam(1,1)*tt180*tt2+tt170+tt169
tt271 = 3*tt5*mu(1,1)*F(1,2)
tt272 = mu(1,1)*tt1
tt273 = tt272+tt271
tt274 = tt273*F(1,3)
tt275 = -2*tt7*mu(1,1)*F(1,4)
tt276 = (tt275+tt274+tt241+tt240)*F(2,1)
tt277 = 2*lam(1,1)*tt1+mu(1,1)*F(1,2)
tt278 = -2*F(1,1)*mu(1,1)*F(1,4)
tt279 = tt278+tt277*F(1,3)+tt260
tt280 = tt279*tt2
tt281 = 2*F(1,1)*lam(1,1)*tt1*F(1,4)+tt250-4*F(1,1)*lam(1,1)*tt66&
&
tt282 = tt281*tt69
tt283 = -2*tt90*mu(1,1)
tt284 = -tt7*mu(1,1)
tt285 = F(1,1)*mu(1,1)*tt3
tt286 = (tt285+tt284)*F(1,3)
tt287 = -2*tt5*mu(1,1)*F(1,2)*F(1,4)
tt288 = -2*tt5*mu(1,1)
tt289 = 4*mu(1,1)*tt3
tt290 = tt262-6*F(1,1)*lam(1,1)*tt3
tt291 = tt177+tt290*F(1,3)+tt289+tt288
tt292 = tt291*tt4
tt293 = (-6*tt5*lam(1,1)*tt3*F(1,4))-10*F(1,1)*lam(1,1)*tt1*F(1,3&
&)+16*tt5*lam(1,1)*tt1
tt294 = tt293*tt2
tt295 = (tt294+tt292+tt287+tt286+tt227+tt283)*F(2,2)
tt296 = tt103+tt117
tt297 = tt296*F(1,2)*F(1,3)
tt298 = 6*tt7*lam(1,1)*F(1,2)*F(1,4)+18*tt5*lam(1,1)*tt3*F(1,3)-2&
&4*tt7*lam(1,1)*tt3
tt299 = tt298*tt4+(tt297+tt246)*F(2,1)
tt300 = tt299*tt6
tt301 = 2*tt5*mu(1,1)
tt302 = -2*tt7*lam(1,1)
tt303 = -F(1,1)*mu(1,1)
tt304 = tt266-14*tt7*lam(1,1)*F(1,2)*F(1,3)+16*tt90*lam(1,1)*F(1,&
&2)
tt305 = tt304*F(2,1)+(tt303+tt302)*F(1,3)+tt301
tt306 = tt305*tt8
tt307 = 4*tt90*lam(1,1)*F(1,3)-4*tt225*lam(1,1)
tt308 = tt307*tt25
tt309 = (tt72+tt288)*F(2,1)
tt310 = -6*tt5*lam(1,1)*tt1*tt4
tt311 = (-2*tt90*lam(1,1)*F(1,2)*tt8)+6*tt7*lam(1,1)*tt3*F(2,1)*t&
&t6+(tt310+tt309)*F(2,2)+2*F(1,1)*lam(1,1)*tt66*tt2+tt184+tt183
tt312 = tt311*F(2,3)
tt313 = tt138+tt209
tt314 = tt225*lam(1,1)*tt8-3*tt90*lam(1,1)*F(1,2)*F(2,1)*tt6+3*tt&
&7*lam(1,1)*tt3*tt4*F(2,2)-tt5*lam(1,1)*tt1*tt2+tt313*tt4+tt169+tt&
&204
tt315 = tt314*tt43+(tt312+tt308+tt306+tt300+tt295+tt282+tt280+tt2&
&76)*F(2,4)+tt270*tt36+(tt268+tt265+tt259+tt253+tt239+tt237+tt234)&
&*F(2,3)+tt226*tt27+tt224*tt25+tt221*tt8+tt216*tt6+(tt203+tt199+tt&
&194+tt191)*F(2,2)+tt181*tt182+tt179*tt69+tt176*tt1*tt2+tt175*tt4
tt316 = 3*tt132*tt4
tt317 = 4*tt134*tt2
tt318 = -6*F(1,1)*lam(1,1)*tt1*tt4
tt319 = 2*tt193*tt3*F(2,1)
tt320 = 3*tt198*tt4
tt321 = 4*tt202*tt2
tt322 = 3*tt215*tt4+2*tt212*F(2,1)+tt208*F(1,2)
tt323 = 2*tt220*F(2,1)+tt197+tt196+tt195
tt324 = 3*tt236*tt4
tt325 = 4*tt238*tt2
tt326 = 2*tt248*F(2,1)
tt327 = 3*tt251*tt4
tt328 = (tt327+tt326)*F(2,2)
tt329 = 2*tt257*F(2,1)+tt256+tt210+tt72+tt254
tt330 = tt329*tt6
tt331 = tt263*tt8
tt332 = (-3*tt5*lam(1,1)*tt1*tt6)+6*F(1,1)*lam(1,1)*tt66*F(2,1)*F&
&(2,2)-3*lam(1,1)*tt180*tt4
tt333 = 3*tt279*tt4
tt334 = 4*tt281*tt2
tt335 = 2*tt291*F(2,1)
tt336 = 3*tt293*tt4
tt337 = (tt336+tt335)*F(2,2)
tt338 = 2*tt298*F(2,1)+tt297+tt246
tt339 = tt338*tt6
tt340 = tt304*tt8
tt341 = -12*tt5*lam(1,1)*tt1*F(2,1)
tt342 = 6*tt7*lam(1,1)*tt3*tt6+(tt341+tt72+tt288)*F(2,2)+6*F(1,1)&
&*lam(1,1)*tt66*tt4
tt343 = tt342*F(2,3)
tt344 = (-3*tt90*lam(1,1)*F(1,2)*tt6)+6*tt7*lam(1,1)*tt3*F(2,1)*F&
&(2,2)-3*tt5*lam(1,1)*tt1*tt4+2*tt313*F(2,1)
tt345 = tt344*tt43+(tt343+tt340+tt339+tt337+tt334+tt333+tt275+tt2&
&74+tt241+tt240)*F(2,4)+tt332*tt36+(tt331+tt330+tt328+tt325+tt324+&
&tt233+tt229+tt228+tt227)*F(2,3)+tt223*tt25+tt323*tt8+tt322*tt6+(t&
&t321+tt320+tt319+tt190+tt186+tt184+tt183)*F(2,2)+5*tt181*tt69+4*t&
&t179*tt2+3*tt176*tt1*tt4+2*tt175*F(2,1)
tt346 = (-6*tt5*F(1,2)*tt6)+12*F(1,1)*tt3*F(2,1)*F(2,2)-6*tt1*tt4&
&
tt347 = (-(12*tt3*F(2,1)*F(2,2)-12*F(1,1)*F(1,2)*tt6)*tt45*tt315)&
&+2*tt346*tt44*tt168*tt315-tt346*tt45*tt167-tt44*tt45*tt345+tt10*(&
&((-12*tt7*lam(1,1)*F(1,2)*tt6)+18*tt5*lam(1,1)*tt3*F(2,1)*F(2,2)+&
&tt318+4*F(1,1)*mu(1,1)*F(2,1))*tt43+((18*tt5*lam(1,1)*tt3*tt6+(tt&
&141-24*F(1,1)*lam(1,1)*tt1*F(2,1))*F(2,2)+6*lam(1,1)*tt66*tt4)*F(&
&2,3)+tt154*tt8+(2*tt149*F(2,1)+tt84+tt102)*tt6+(3*tt146*tt4+2*tt1&
&44*F(2,1))*F(2,2)+tt317+tt316+tt129+tt128+tt97+tt96)*F(2,4)+(6*la&
&m(1,1)*tt66*F(2,1)*F(2,2)-6*F(1,1)*lam(1,1)*tt1*tt6)*tt36+(tt118*&
&tt8+(2*tt113*F(2,1)+tt112+tt76+tt111)*tt6+(3*tt108*tt4+2*tt106*F(&
&2,1))*F(2,2)+tt94+tt93+tt92)*F(2,3)+tt88*tt25+(2*tt85*F(2,1)+tt63&
&+tt62)*tt8+(3*tt80*tt4+2*tt77*F(2,1)+tt74*F(1,2))*tt6+(4*tt68*tt2&
&+3*tt64*tt4+2*tt60*tt3*F(2,1)+tt57+tt53+tt52+tt51)*F(2,2)+2*tt50*&
&F(2,1))
tt348 = 4*F(1,1)*mu(1,1)*F(1,2)
tt349 = -2*mu(1,1)*tt2
tt350 = -48*lam(1,1)*tt3*F(1,4)
tt351 = -4*mu(1,1)*F(1,2)*F(1,3)
tt352 = -2*mu(1,1)*F(2,1)
tt353 = 12*F(1,1)*lam(1,1)*F(1,3)
tt354 = -8*F(1,1)*mu(1,1)*F(1,2)
tt355 = 6*mu(1,1)*tt3
tt356 = 2*mu(1,1)*tt2
tt357 = (6*lam(1,1)*tt3*F(1,4)-16*lam(1,1)*tt1)*tt69
tt358 = 8*lam(1,1)*tt1*tt2
tt359 = -6*F(1,1)*lam(1,1)*tt3*tt2
tt360 = 2*tt5*mu(1,1)*F(1,2)
tt361 = 4*mu(1,1)*tt1
tt362 = 2*mu(1,1)*F(1,2)*tt22
tt363 = tt210+tt355+tt288
tt364 = tt363*F(1,4)+tt362+8*F(1,1)*mu(1,1)*F(1,2)*F(1,3)+tt361+t&
&t360
tt365 = (-3*lam(1,1)*tt3*tt13)+16*lam(1,1)*tt1*F(1,4)-15*lam(1,1)&
&*tt66
tt366 = tt355+tt82
tt367 = tt366*F(1,3)
tt368 = -12*F(1,1)*mu(1,1)*F(1,2)*F(1,4)
tt369 = (tt368+tt367+tt188+tt217)*F(2,1)
tt370 = 2*tt193*F(1,2)*tt4
tt371 = -2*F(1,1)*mu(1,1)
tt372 = 2*mu(1,1)*F(1,3)
tt373 = tt372+tt371
tt374 = tt373*tt2
tt375 = -48*F(1,1)*lam(1,1)*tt3
tt376 = 6*lam(1,1)*tt3*F(1,3)
tt377 = 6*F(1,1)*lam(1,1)*F(1,2)*tt13+(tt376+tt375)*F(1,4)-16*lam&
&(1,1)*tt1*F(1,3)+60*F(1,1)*lam(1,1)*tt1
tt378 = tt377*tt69
tt379 = -4*F(1,1)*mu(1,1)*F(1,2)*F(1,3)
tt380 = (-3*tt5*lam(1,1)*tt13)+(tt34+tt33)*F(1,4)-3*lam(1,1)*tt3*&
&tt22+48*F(1,1)*lam(1,1)*tt3*F(1,3)-90*tt5*lam(1,1)*tt3
tt381 = tt380*tt2+tt132*tt4+tt208*F(2,1)+(tt210+tt254)*F(1,4)+tt3&
&62+tt379+tt360
tt382 = 6*tt5*lam(1,1)*F(1,3)-16*tt7*lam(1,1)
tt383 = tt382*F(1,4)+6*F(1,1)*lam(1,1)*F(1,2)*tt22-48*tt5*lam(1,1&
&)*F(1,2)*F(1,3)+60*tt7*lam(1,1)*F(1,2)
tt384 = tt383*tt4+tt373*F(2,1)
tt385 = (-3*tt5*lam(1,1)*tt22)+tt157-15*tt90*lam(1,1)
tt386 = -4*tt5*mu(1,1)*F(1,2)
tt387 = -8*mu(1,1)*tt1
tt388 = tt209-3*mu(1,1)*tt3
tt389 = tt388*F(1,4)
tt390 = (tt389+tt379+tt387+tt386)*F(2,1)
tt391 = tt142+tt153
tt392 = tt391*F(1,4)+4*mu(1,1)*F(1,2)
tt393 = tt392*tt2
tt394 = tt67+20*lam(1,1)*tt66
tt395 = tt394*tt69
tt396 = 6*F(1,1)*mu(1,1)*tt3
tt397 = -6*mu(1,1)*tt3*F(1,3)
tt398 = 6*F(1,1)*mu(1,1)*F(1,2)*F(1,4)
tt399 = 12*F(1,1)*lam(1,1)*F(1,2)*F(1,4)
tt400 = tt399+tt192
tt401 = tt400*tt4
tt402 = 8*lam(1,1)*tt1*F(1,3)
tt403 = 42*F(1,1)*lam(1,1)*tt3*F(1,4)+tt402-64*F(1,1)*lam(1,1)*tt&
&1
tt404 = tt403*tt2
tt405 = (tt404+tt401+tt398+tt397+tt396+tt187)*F(2,2)
tt406 = -4*mu(1,1)*F(1,2)
tt407 = tt255*F(1,4)
tt408 = (-36*tt5*lam(1,1)*F(1,2)*F(1,4))-18*F(1,1)*lam(1,1)*tt3*F&
&(1,3)+72*tt5*lam(1,1)*tt3
tt409 = tt408*tt4+(tt407+tt406)*F(2,1)
tt410 = tt409*tt6
tt411 = 10*tt7*lam(1,1)*F(1,4)+12*tt5*lam(1,1)*F(1,2)*F(1,3)-32*t&
&t7*lam(1,1)*F(1,2)
tt412 = tt411*F(2,1)+tt76+tt75
tt413 = tt412*tt8
tt414 = tt222*tt25
tt415 = tt116-9*tt5*lam(1,1)*tt3*F(2,1)
tt416 = 2*tt7*lam(1,1)*F(1,2)*tt8+tt415*tt6+12*F(1,1)*lam(1,1)*tt&
&1*tt4*F(2,2)-5*lam(1,1)*tt66*tt2+tt361+tt360
tt417 = tt100*F(1,3)
tt418 = (tt417+tt396+tt187)*F(2,1)
tt419 = tt104+mu(1,1)
tt420 = tt419*F(1,3)+tt75
tt421 = tt420*tt2
tt422 = 6*F(1,1)*lam(1,1)*tt3*F(1,4)
tt423 = tt422+tt402+tt200
tt424 = tt423*tt69
tt425 = 2*F(1,1)*mu(1,1)*F(1,2)*F(1,3)
tt426 = -2*tt5*mu(1,1)*F(1,4)
tt427 = tt131+tt34+8*mu(1,1)*F(1,2)
tt428 = tt427*tt4
tt429 = (-12*tt5*lam(1,1)*F(1,2)*F(1,4))-30*F(1,1)*lam(1,1)*tt3*F&
&(1,3)+48*tt5*lam(1,1)*tt3
tt430 = tt429*tt2
tt431 = (tt430+tt428+tt426+tt425+tt386)*F(2,2)
tt432 = tt296*F(1,3)
tt433 = 6*tt7*lam(1,1)*F(1,4)+36*tt5*lam(1,1)*F(1,2)*F(1,3)-48*tt&
&7*lam(1,1)*F(1,2)
tt434 = tt433*tt4+(tt432+tt192)*F(2,1)
tt435 = tt434*tt6
tt436 = 16*tt90*lam(1,1)-14*tt7*lam(1,1)*F(1,3)
tt437 = tt436*F(2,1)*tt8
tt438 = 8*F(1,1)*lam(1,1)*tt1*tt2
tt439 = -4*mu(1,1)*F(1,2)*F(2,1)
tt440 = -18*tt5*lam(1,1)*tt3*tt4
tt441 = (tt440+tt439)*F(2,2)
tt442 = 12*tt7*lam(1,1)*F(1,2)*F(2,1)*tt6
tt443 = -2*tt90*lam(1,1)*tt8
tt444 = tt443+tt442+tt441+tt438+tt188+tt217
tt445 = tt444*F(2,3)
tt446 = (-3*tt90*lam(1,1)*F(2,1)*tt6)+6*tt7*lam(1,1)*F(1,2)*tt4*F&
&(2,2)-3*tt5*lam(1,1)*tt3*tt2+2*mu(1,1)*F(1,2)*tt4+tt360
tt447 = tt446*tt43+(tt445+tt437+tt435+tt431+tt424+tt421+tt418)*F(&
&2,4)+tt416*tt36+(tt414+tt413+tt410+tt405+tt395+tt393+tt390)*F(2,3&
&)+tt385*F(2,1)*tt25+tt384*tt8+tt381*tt6+(tt378+tt374+tt370+tt369)&
&*F(2,2)+tt365*tt182+tt132*tt69+3*tt176*tt3*tt2+tt364*tt4
tt448 = (-6*tt5*F(2,1)*tt6)+12*F(1,1)*F(1,2)*tt4*F(2,2)-6*tt3*tt2&
&
tt449 = (-(12*F(1,2)*tt4*F(2,2)-12*F(1,1)*F(2,1)*tt6)*tt45*tt315)&
&+2*tt448*tt44*tt168*tt315-tt448*tt45*tt167-tt44*tt45*tt447+tt10*(&
&((-12*tt7*lam(1,1)*F(2,1)*tt6)+18*tt5*lam(1,1)*F(1,2)*tt4*F(2,2)+&
&tt359+tt348)*tt43+(((-8*tt7*lam(1,1)*tt8)+36*tt5*lam(1,1)*F(1,2)*&
&F(2,1)*tt6-36*F(1,1)*lam(1,1)*tt3*tt4*F(2,2)+tt358+tt55+tt82)*F(2&
&,3)+(64*tt7*lam(1,1)-42*tt5*lam(1,1)*F(1,3))*F(2,1)*tt8+((18*tt5*&
&lam(1,1)*F(1,4)+72*F(1,1)*lam(1,1)*F(1,2)*F(1,3)-144*tt5*lam(1,1)&
&*F(1,2))*tt4+(tt353+tt59)*F(2,1))*tt6+(((-24*F(1,1)*lam(1,1)*F(1,&
&2)*F(1,4))-30*lam(1,1)*tt3*F(1,3)+96*F(1,1)*lam(1,1)*tt3)*tt2-12*&
&lam(1,1)*F(1,2)*F(1,3)*tt4-4*F(1,1)*mu(1,1)*F(1,4)+tt196+tt354)*F&
&(2,2)+tt357+tt356+(6*F(1,1)*mu(1,1)*F(1,3)+tt355+tt54)*F(2,1))*F(&
&2,4)+(6*tt5*lam(1,1)*F(1,2)*tt8-18*F(1,1)*lam(1,1)*tt3*F(2,1)*tt6&
&+12*lam(1,1)*tt1*tt4*F(2,2)+tt348)*tt36+(tt87*tt25+((30*tt5*lam(1&
&,1)*F(1,4)+24*F(1,1)*lam(1,1)*F(1,2)*F(1,3)-96*tt5*lam(1,1)*F(1,2&
&))*F(2,1)+tt21)*tt8+(((-72*F(1,1)*lam(1,1)*F(1,2)*F(1,4))-18*lam(&
&1,1)*tt3*F(1,3)+144*F(1,1)*lam(1,1)*tt3)*tt4-12*F(1,1)*lam(1,1)*F&
&(1,4)*F(2,1))*tt6+((42*lam(1,1)*tt3*F(1,4)-64*lam(1,1)*tt1)*tt2+(&
&12*lam(1,1)*F(1,2)*F(1,4)+tt59)*tt4+6*mu(1,1)*F(1,2)*F(1,4)+tt355&
&+tt54)*F(2,2)+(2*F(1,1)*mu(1,1)*F(1,4)+tt351+tt354)*F(2,1))*F(2,3&
&)+((-6*F(1,1)*lam(1,1)*tt22)+tt40-60*tt7*lam(1,1))*F(2,1)*tt25+((&
&(tt353-48*tt5*lam(1,1))*F(1,4)+6*lam(1,1)*F(1,2)*tt22-96*F(1,1)*l&
&am(1,1)*F(1,2)*F(1,3)+180*tt5*lam(1,1)*F(1,2))*tt4+tt352)*tt8+(((&
&-6*F(1,1)*lam(1,1)*tt13)+(96*F(1,1)*lam(1,1)*F(1,2)-12*lam(1,1)*F&
&(1,2)*F(1,3))*F(1,4)+48*lam(1,1)*tt3*F(1,3)-180*F(1,1)*lam(1,1)*t&
&t3)*tt2+tt74*F(2,1)+(tt76+tt111)*F(1,4)+tt351+tt348)*tt6+((6*lam(&
&1,1)*F(1,2)*tt13+tt350+60*lam(1,1)*tt1)*tt69+tt349+2*tt60*F(1,2)*&
&tt4+((-12*mu(1,1)*F(1,2)*F(1,4))-12*F(1,1)*mu(1,1)*F(1,3)+tt55+tt&
&82)*F(2,1))*F(2,2)+((tt76+tt141)*F(1,4)+8*mu(1,1)*F(1,2)*F(1,3)+t&
&t348)*tt4)
tt450 = 2*tt258*F(2,2)
tt451 = 3*tt264*tt6
tt452 = 4*tt267*tt8
tt453 = 3*tt7*lam(1,1)*tt3*tt6+2*tt269*F(2,2)+3*F(1,1)*lam(1,1)*t&
&t66*tt4
tt454 = 2*tt299*F(2,2)
tt455 = 3*tt305*tt6
tt456 = 4*tt307*tt8
tt457 = (-6*tt90*lam(1,1)*F(1,2)*tt6)+12*tt7*lam(1,1)*tt3*F(2,1)*&
&F(2,2)+tt310+tt309
tt458 = tt457*F(2,3)
tt459 = 3*tt225*lam(1,1)*tt6-6*tt90*lam(1,1)*F(1,2)*F(2,1)*F(2,2)&
&+3*tt7*lam(1,1)*tt3*tt4
tt460 = tt459*tt43+(tt458+tt456+tt455+tt454+tt294+tt292+tt287+tt2&
&86+tt227+tt283)*F(2,4)+tt453*tt36+(tt452+tt451+tt450+tt252+tt249+&
&tt245+tt242+tt241+tt240)*F(2,3)+5*tt226*tt25+4*tt224*tt8+3*tt221*&
&tt6+2*tt216*F(2,2)+tt203+tt199+tt194+tt191
tt461 = 6*tt7*tt6-12*tt5*F(1,2)*F(2,1)*F(2,2)+6*F(1,1)*tt3*tt4
tt462 = (-(18*tt5*tt6-24*F(1,1)*F(1,2)*F(2,1)*F(2,2)+6*tt3*tt4)*t&
&t45*tt315)+2*tt461*tt44*tt168*tt315-tt461*tt45*tt167-tt44*tt45*tt&
&460+tt10*((15*tt90*lam(1,1)*tt6-24*tt7*lam(1,1)*F(1,2)*F(2,1)*F(2&
&,2)+9*tt5*lam(1,1)*tt3*tt4)*tt43+(((-24*tt7*lam(1,1)*F(1,2)*tt6)+&
&36*tt5*lam(1,1)*tt3*F(2,1)*F(2,2)+tt162+tt161)*F(2,3)+4*tt158*tt8&
&+3*tt155*tt6+2*tt150*F(2,2)+tt147+tt145+tt140+tt139+tt92+tt136)*F&
&(2,4)+(9*tt5*lam(1,1)*tt3*tt6+2*tt126*F(2,2)+3*lam(1,1)*tt66*tt4)&
&*tt36+(4*tt124*tt8+3*tt119*tt6+2*tt114*F(2,2)+tt109+tt107+tt101+t&
&t97+tt96)*F(2,3)+5*tt91*tt25+4*tt89*tt8+3*tt86*tt6+2*tt81*F(2,2)+&
&tt70+tt65+tt61+tt58)
tt463 = -2*mu(1,1)*tt4
tt464 = 16*lam(1,1)*tt1-6*lam(1,1)*tt3*F(1,4)
tt465 = tt399+tt376+tt375
tt466 = -2*mu(1,1)
tt467 = (-6*tt5*lam(1,1)*F(1,4))+tt34+tt33
tt468 = 2*mu(1,1)*tt3*F(1,3)
tt469 = -2*F(1,1)*mu(1,1)*F(1,2)*F(1,4)
tt470 = tt469+tt468+4*F(1,1)*mu(1,1)*tt3
tt471 = tt185*F(2,1)
tt472 = (-6*F(1,1)*lam(1,1)*tt3*F(1,4))-2*lam(1,1)*tt1*F(1,3)+16*&
&F(1,1)*lam(1,1)*tt1
tt473 = tt472*tt2+tt373*tt4+tt469+tt468+tt205+tt187
tt474 = 6*tt5*lam(1,1)*F(1,2)*F(1,4)+6*F(1,1)*lam(1,1)*tt3*F(1,3)&
&-24*tt5*lam(1,1)*tt3
tt475 = tt474*tt4+tt132*F(2,1)
tt476 = (-2*tt7*lam(1,1)*F(1,4))+tt122+tt121
tt477 = tt476*F(2,1)+tt372+tt371
tt478 = 2*tt7*lam(1,1)*F(1,3)-4*tt90*lam(1,1)
tt479 = tt318-2*F(1,1)*mu(1,1)*F(2,1)
tt480 = 6*tt5*lam(1,1)*tt3*F(2,1)+tt62
tt481 = (-2*tt7*lam(1,1)*F(1,2)*tt25)+tt480*tt8+tt479*tt6+(tt160+&
&tt52)*F(2,2)-2*F(1,1)*mu(1,1)*tt3*F(2,1)
tt482 = tt290*tt4
tt483 = -10*F(1,1)*lam(1,1)*tt1*tt2
tt484 = 18*tt5*lam(1,1)*tt3*tt4+tt296*F(1,2)*F(2,1)
tt485 = (-14*tt7*lam(1,1)*F(1,2)*F(2,1))+tt303+tt302
tt486 = 4*tt90*lam(1,1)*tt25+tt485*tt8+tt484*tt6+(tt483+tt482+tt2&
&85+tt284)*F(2,2)+2*lam(1,1)*tt66*tt69+tt277*tt2+tt273*F(2,1)
tt487 = tt486*F(2,4)+tt481*F(2,3)+tt478*tt27+tt477*tt25+tt475*tt8&
&+tt473*tt6+(tt135+tt133+tt471)*F(2,2)+tt470*tt4
tt488 = tt10*((16*tt7*lam(1,1)*tt25+((-42*tt5*lam(1,1)*F(1,2)*F(2&
&,1))+tt153+tt152)*tt8+(36*F(1,1)*lam(1,1)*tt3*tt4+12*F(1,1)*lam(1&
&,1)*F(1,2)*F(2,1))*tt6+((-10*lam(1,1)*tt1*tt2)+tt143*tt4+tt138+tt&
&137)*F(2,2)+6*F(1,1)*mu(1,1)*F(1,2)*F(2,1))*F(2,4)+((-6*tt5*lam(1&
&,1)*F(1,2)*tt25)+12*F(1,1)*lam(1,1)*tt3*F(2,1)*tt8+(tt352-6*lam(1&
&,1)*tt1*tt4)*tt6-2*mu(1,1)*tt3*F(2,1))*F(2,3)+tt382*tt27+(tt467*F&
&(2,1)+tt466)*tt25+tt465*tt4*tt8+(tt464*tt2+tt463+tt177+tt72+tt54)&
&*tt6-12*F(1,1)*mu(1,1)*F(1,2)*F(2,1)*F(2,2)+(tt177+tt289)*tt4)-tt&
&44*tt45*tt487
tt489 = tt311*F(2,4)+2*tt270*F(2,3)+tt268+tt265+tt259+tt253+tt239&
&+tt237+tt234
tt490 = tt10*(tt164*F(2,4)+2*tt127*F(2,3)+tt125+tt120+tt115+tt110&
&+tt95)-tt44*tt45*tt489
tt491 = 4*F(1,1)*mu(1,1)*F(1,4)
tt492 = (tt112+tt79+tt78)*tt2
tt493 = ((6*tt5*lam(1,1)*F(1,4)+tt84+tt83)*tt4+2*mu(1,1)*F(2,1))*&
&tt8
tt494 = tt87*F(2,1)*tt25
tt495 = 2*lam(1,1)*tt1*tt69
tt496 = 2*tt5*mu(1,1)*F(1,4)
tt497 = tt496+tt173+tt97+tt172
tt498 = 4*lam(1,1)*tt66-2*lam(1,1)*tt1*F(1,4)
tt499 = tt189*F(2,1)
tt500 = tt77*tt2
tt501 = tt422+tt201+tt200
tt502 = tt501*tt69
tt503 = (-6*tt5*lam(1,1)*F(1,2)*F(1,4))+tt214+tt213
tt504 = tt503*tt2+tt64*tt4+tt496+tt173+tt206
tt505 = 2*tt7*lam(1,1)*F(1,4)+tt219+tt218
tt506 = tt505*tt4+tt77*F(2,1)
tt507 = tt247*tt4
tt508 = 14*F(1,1)*lam(1,1)*tt1*tt2
tt509 = tt440+tt255*F(1,2)*F(2,1)
tt510 = 10*tt7*lam(1,1)*F(1,2)*F(2,1)+tt262+tt261
tt511 = (-2*tt90*lam(1,1)*tt25)+tt510*tt8+tt509*tt6+(tt508+tt507+&
&tt244+tt243)*F(2,2)-4*lam(1,1)*tt66*tt69+tt235*tt2+tt232*F(2,1)
tt512 = -2*mu(1,1)*F(1,2)*tt4
tt513 = -6*tt5*lam(1,1)*tt3*tt2
tt514 = (-2*tt90*lam(1,1)*F(2,1)*tt8)+6*tt7*lam(1,1)*F(1,2)*tt4*t&
&t6+(tt513+tt512+tt172)*F(2,2)+2*F(1,1)*lam(1,1)*tt1*tt69-2*F(1,1)&
&*mu(1,1)*tt2-2*tt7*mu(1,1)*F(2,1)
tt515 = tt514*F(2,4)+tt511*F(2,3)+tt222*F(2,1)*tt25+tt506*tt8+tt5&
&04*tt6+(tt502+tt500+tt499)*F(2,2)+tt498*tt182+tt64*tt69+tt497*tt4
tt516 = tt10*(((-8*tt7*lam(1,1)*F(2,1)*tt8)+18*tt5*lam(1,1)*F(1,2&
&)*tt4*tt6+(tt47-12*F(1,1)*lam(1,1)*tt3*tt2)*F(2,2)+tt495+tt349-6*&
&tt5*mu(1,1)*F(2,1))*F(2,4)+((-8*tt7*lam(1,1)*tt25)+(30*tt5*lam(1,&
&1)*F(1,2)*F(2,1)+mu(1,1)+tt117)*tt8+((-36*F(1,1)*lam(1,1)*tt3*tt4&
&)-12*F(1,1)*lam(1,1)*F(1,2)*F(2,1))*tt6+(14*lam(1,1)*tt1*tt2+tt10&
&5*tt4+tt99+tt98)*F(2,2)+2*F(1,1)*mu(1,1)*F(1,2)*F(2,1))*F(2,3)+tt&
&494+tt493+(tt492+tt491+tt48+tt73)*tt6+(tt357+tt356+tt56*F(2,1))*F&
&(2,2)+(tt491+tt48+tt47)*tt4)-tt44*tt45*tt515
tt517 = 2*tt314*F(2,4)+tt312+tt308+tt306+tt300+tt295+tt282+tt280+&
&tt276
tt518 = tt10*(2*tt166*F(2,4)+tt165+tt159+tt156+tt151+tt148+tt135+&
&tt133+tt130)-tt44*tt45*tt517
tt519 = 6*tt7*lam(1,1)*tt3*F(2,2)
tt520 = -36*tt5*lam(1,1)*tt3*F(2,1)
tt521 = (-((-6*tt5*tt6)+24*F(1,1)*F(1,2)*F(2,1)*F(2,2)-18*tt3*tt4&
&)*tt45*tt315)+2*tt346*tt448*tt168*tt315-tt346*tt45*tt447-tt448*tt&
&45*tt345+tt10*(((-3*tt90*lam(1,1)*tt6)+12*tt7*lam(1,1)*F(1,2)*F(2&
&,1)*F(2,2)-9*tt5*lam(1,1)*tt3*tt4+4*mu(1,1)*F(1,2)*F(2,1))*tt43+(&
&(12*tt7*lam(1,1)*F(1,2)*tt6+(tt520+tt406)*F(2,2)+24*F(1,1)*lam(1,&
&1)*tt1*tt4)*F(2,3)+tt436*tt8+(2*tt433*F(2,1)+tt432+tt192)*tt6+(3*&
&tt429*tt4+2*tt427*F(2,1))*F(2,2)+4*tt423*tt2+3*tt420*tt4+tt417+tt&
&396+tt187)*F(2,4)+((-9*tt5*lam(1,1)*tt3*tt6)+24*F(1,1)*lam(1,1)*t&
&t1*F(2,1)*F(2,2)-15*lam(1,1)*tt66*tt4)*tt36+(tt411*tt8+(2*tt408*F&
&(2,1)+tt407+tt406)*tt6+(3*tt403*tt4+2*tt400*F(2,1))*F(2,2)+4*tt39&
&4*tt2+3*tt392*tt4+tt389+tt379+tt387+tt386)*F(2,3)+tt385*tt25+(2*t&
&t383*F(2,1)+tt372+tt371)*tt8+(3*tt380*tt4+2*tt132*F(2,1)+tt54+tt2&
&07)*tt6+(4*tt377*tt2+3*tt373*tt4+4*tt193*F(1,2)*F(2,1)+tt368+tt36&
&7+tt188+tt217)*F(2,2)+5*tt365*tt69+4*tt132*tt2+9*tt176*tt3*tt4+2*&
&tt364*F(2,1))
tt522 = (-(12*F(1,1)*tt3*F(2,1)-12*tt5*F(1,2)*F(2,2))*tt45*tt315)&
&+2*tt461*tt346*tt168*tt315-tt461*tt45*tt345-tt346*tt45*tt460+tt10&
&*((6*tt7*lam(1,1)*tt3*F(2,1)-6*tt90*lam(1,1)*F(1,2)*F(2,2))*tt43+&
&((12*tt7*lam(1,1)*tt3*F(2,2)+tt341+tt72+tt288)*F(2,3)+3*tt304*tt6&
&+2*tt338*F(2,2)+tt336+tt335)*F(2,4)+(6*F(1,1)*lam(1,1)*tt66*F(2,1&
&)-6*tt5*lam(1,1)*tt1*F(2,2))*tt36+(3*tt263*tt6+2*tt329*F(2,2)+tt3&
&27+tt326)*F(2,3)+4*tt223*tt8+3*tt323*tt6+2*tt322*F(2,2)+tt321+tt3&
&20+tt319+tt190+tt186+tt184+tt183)
tt523 = tt10*(((-14*tt7*lam(1,1)*F(1,2)*tt8)+(36*tt5*lam(1,1)*tt3&
&*F(2,1)+tt296*F(1,2))*tt6+(2*tt290*F(2,1)-30*F(1,1)*lam(1,1)*tt1*&
&tt4)*F(2,2)+8*lam(1,1)*tt66*tt2+3*tt277*tt4+tt272+tt271)*F(2,4)+(&
&6*tt5*lam(1,1)*tt3*tt8+(tt371-12*F(1,1)*lam(1,1)*tt1*F(2,1))*tt6+&
&6*lam(1,1)*tt66*tt4*F(2,2)+tt205)*F(2,3)+tt476*tt25+(2*tt474*F(2,&
&1)+tt131+tt116)*tt8+(3*tt472*tt4+2*tt373*F(2,1))*tt6+(tt317+tt316&
&+tt97+tt51)*F(2,2)+2*tt470*F(2,1))-tt346*tt45*tt487
tt524 = tt10*(tt342*F(2,4)+2*tt332*F(2,3)+tt331+tt330+tt328+tt325&
&+tt324+tt233+tt229+tt228+tt227)-tt346*tt45*tt489
tt525 = tt10*((tt443+tt442+tt441+tt438-6*F(1,1)*mu(1,1)*tt4+tt217&
&)*F(2,4)+(10*tt7*lam(1,1)*F(1,2)*tt8+(tt520+tt255*F(1,2))*tt6+(42&
&*F(1,1)*lam(1,1)*tt1*tt4+2*tt247*F(2,1))*F(2,2)-16*lam(1,1)*tt66*&
&tt2+3*tt235*tt4+tt231+tt230)*F(2,3)+tt414+(2*tt505*F(2,1)+tt76+tt&
&75)*tt8+(3*tt503*tt4+2*tt64*F(2,1))*tt6+(4*tt501*tt2+3*tt77*tt4+t&
&t188+tt187)*F(2,2)+5*tt498*tt69+4*tt64*tt2+2*tt497*F(2,1))-tt346*&
&tt45*tt515
tt526 = tt10*(2*tt344*F(2,4)+tt343+tt340+tt339+tt337+tt334+tt333+&
&tt275+tt274+tt241+tt240)-tt346*tt45*tt517
tt527 = 12*mu(1,1)*tt3
tt528 = 2*mu(1,1)*tt22
tt529 = 2*mu(1,1)*tt69
tt530 = -4*F(1,1)*mu(1,1)*F(1,3)
tt531 = -4*tt5*mu(1,1)
tt532 = 24*lam(1,1)*tt3*F(1,3)
tt533 = -36*tt5*lam(1,1)*F(1,2)*tt4
tt534 = -6*tt5*lam(1,1)*F(1,2)*tt2
tt535 = (-(12*F(1,1)*F(1,2)*tt4-12*tt5*F(2,1)*F(2,2))*tt45*tt315)&
&+2*tt461*tt448*tt168*tt315-tt461*tt45*tt447-tt448*tt45*tt460+tt10&
&*((6*tt7*lam(1,1)*F(1,2)*tt4-6*tt90*lam(1,1)*F(2,1)*F(2,2))*tt43+&
&(((-6*tt90*lam(1,1)*tt6)+24*tt7*lam(1,1)*F(1,2)*F(2,1)*F(2,2)+tt4&
&40+tt439)*F(2,3)+3*tt436*F(2,1)*tt6+2*tt434*F(2,2)+tt430+tt428+tt&
&426+tt425+tt386)*F(2,4)+(6*tt7*lam(1,1)*F(1,2)*tt6+2*tt415*F(2,2)&
&+12*F(1,1)*lam(1,1)*tt1*tt4)*tt36+(4*tt222*tt8+3*tt412*tt6+2*tt40&
&9*F(2,2)+tt404+tt401+tt398+tt397+tt396+tt187)*F(2,3)+4*tt385*F(2,&
&1)*tt8+3*tt384*tt6+2*tt381*F(2,2)+tt378+tt374+tt370+tt369)
tt536 = 4*mu(1,1)*F(1,2)*F(1,3)
tt537 = tt10*(((-14*tt7*lam(1,1)*F(2,1)*tt8)+(36*tt5*lam(1,1)*F(1&
&,2)*tt4+tt296*F(2,1))*tt6+((-30*F(1,1)*lam(1,1)*tt3*tt2)-12*F(1,1&
&)*lam(1,1)*F(1,2)*tt4+tt260)*F(2,2)+8*lam(1,1)*tt1*tt69+tt419*tt2&
&+tt100*F(2,1))*F(2,4)+((-2*tt7*lam(1,1)*tt25)+(12*tt5*lam(1,1)*F(&
&1,2)*F(2,1)+tt466)*tt8-18*F(1,1)*lam(1,1)*tt3*tt4*tt6+(tt358+tt55&
&)*F(2,2)-4*F(1,1)*mu(1,1)*F(1,2)*F(2,1))*F(2,3)+tt494+tt493+(tt49&
&2+tt278+tt536+tt47)*tt6+(tt357+tt356+tt366*F(2,1))*F(2,2)+(tt278+&
&tt536+tt73)*tt4)-tt448*tt45*tt487
tt538 = tt10*(tt444*F(2,4)+2*tt416*F(2,3)+tt414+tt413+tt410+tt405&
&+tt395+tt393+tt390)-tt448*tt45*tt489
tt539 = tt10*((6*tt7*lam(1,1)*tt4*tt6+((-12*tt5*lam(1,1)*F(1,2)*t&
&t2)+tt463+tt288)*F(2,2)+6*F(1,1)*lam(1,1)*tt3*tt69)*F(2,4)+(10*tt&
&7*lam(1,1)*F(2,1)*tt8+(tt533+tt255*F(2,1))*tt6+(42*F(1,1)*lam(1,1&
&)*tt3*tt2+12*F(1,1)*lam(1,1)*F(1,2)*tt4+6*F(1,1)*mu(1,1)*F(1,2))*&
&F(2,2)-16*lam(1,1)*tt1*tt69+tt391*tt2+tt388*F(2,1))*F(2,3)+tt382*&
&tt4*tt8+(tt467*tt2+tt463+tt210+tt254)*tt6+(tt465*tt69-12*F(1,1)*m&
&u(1,1)*F(1,2)*F(2,1))*F(2,2)+tt464*tt182-2*mu(1,1)*tt69+tt363*tt4&
&)-tt448*tt45*tt515
tt540 = tt10*(2*tt446*F(2,4)+tt445+tt437+tt435+tt431+tt424+tt421+&
&tt418)-tt448*tt45*tt517
tt541 = tt10*((16*tt90*lam(1,1)*tt8+3*tt485*tt6+2*tt484*F(2,2)+tt&
&483+tt482+tt285+tt284)*F(2,4)+(tt163+3*tt480*tt6+2*tt479*F(2,2)+t&
&t160+tt52)*F(2,3)+5*tt478*tt25+4*tt477*tt8+3*tt475*tt6+2*tt473*F(&
&2,2)+tt135+tt133+tt471)-tt461*tt45*tt487
tt542 = tt10*(tt457*F(2,4)+2*tt453*F(2,3)+tt452+tt451+tt450+tt252&
&+tt249+tt245+tt242+tt241+tt240)-tt461*tt45*tt489
tt543 = tt10*(((-6*tt90*lam(1,1)*F(2,1)*tt6)+12*tt7*lam(1,1)*F(1,&
&2)*tt4*F(2,2)+tt513+tt512+tt172)*F(2,4)+((-8*tt90*lam(1,1)*tt8)+3&
&*tt510*tt6+2*tt509*F(2,2)+tt508+tt507+tt244+tt243)*F(2,3)+4*tt222&
&*F(2,1)*tt8+3*tt506*tt6+2*tt504*F(2,2)+tt502+tt500+tt499)-tt461*t&
&t45*tt515
tt544 = tt10*(2*tt459*F(2,4)+tt458+tt456+tt455+tt454+tt294+tt292+&
&tt287+tt286+tt227+tt283)-tt461*tt45*tt517
tt545 = tt10*tt481
tt546 = tt10*((-2*tt7*lam(1,1)*F(2,1)*tt25)+(6*tt5*lam(1,1)*F(1,2&
&)*tt4+tt352)*tt8+(tt359+tt195)*tt6+(tt495+tt349)*F(2,2)-2*F(1,1)*&
&mu(1,1)*F(1,2)*tt4)
tt547 = tt10*tt486
tt548 = tt10*tt511
tt549 = tt10*tt311
tt550 = tt10*tt514
hes(1,1) = (-(12*F(1,1)*tt8-12*F(1,2)*F(2,1)*tt6)*tt45*tt315)+2*t&
&t44**2*tt168*tt315-2*tt44*tt45*tt167+tt10*((20*tt7*lam(1,1)*tt8-3&
&6*tt5*lam(1,1)*F(1,2)*F(2,1)*tt6+18*F(1,1)*lam(1,1)*tt3*tt4*F(2,2&
&)+tt42+tt18+tt11+tt17)*tt43+(((-24*tt5*lam(1,1)*F(1,2)*tt8)+36*F(&
&1,1)*lam(1,1)*tt3*F(2,1)*tt6+(tt41-12*lam(1,1)*tt1*tt4)*F(2,2)+tt&
&15)*F(2,3)+(tt40-80*tt7*lam(1,1))*tt25+((tt35-84*F(1,1)*lam(1,1)*&
&F(1,2)*F(1,3)+192*tt5*lam(1,1)*F(1,2))*F(2,1)+tt23+tt39)*tt8+((36&
&*F(1,1)*lam(1,1)*F(1,2)*F(1,4)+36*lam(1,1)*tt3*F(1,3)-144*F(1,1)*&
&lam(1,1)*tt3)*tt4+12*lam(1,1)*F(1,2)*F(1,3)*F(2,1))*tt6+((32*lam(&
&1,1)*tt1-12*lam(1,1)*tt3*F(1,4))*tt2-4*mu(1,1)*tt4+tt12-6*F(1,1)*&
&mu(1,1)*F(1,3)+tt28-24*tt5*mu(1,1))*F(2,2)+(tt38+tt37+tt29)*F(2,1&
&))*F(2,4)+(6*F(1,1)*lam(1,1)*tt3*tt8+(tt21-6*lam(1,1)*tt1*F(2,1))&
&*tt6+tt11)*tt36+((tt35+tt34+tt33)*tt25+((60*F(1,1)*lam(1,1)*F(1,2&
&)*F(1,4)+12*lam(1,1)*tt3*F(1,3)-96*F(1,1)*lam(1,1)*tt3)*F(2,1)+12&
&*F(1,1)*lam(1,1)*F(1,4))*tt8+((48*lam(1,1)*tt1-36*lam(1,1)*tt3*F(&
&1,4))*tt4+(tt32+tt31)*F(2,1))*tt6+(tt30+tt29)*F(2,2)+(2*mu(1,1)*F&
&(1,2)*F(1,4)+tt28)*F(2,1))*F(2,3)+tt26*tt27+(tt24*F(2,1)+tt21)*tt&
&25+(tt20*tt4-12*F(1,1)*mu(1,1)+6*F(1,1)*lam(1,1))*tt8+(tt19*tt2+t&
&t18+(12*mu(1,1)-6*lam(1,1))*F(1,2)*F(2,1)+tt14+8*mu(1,1)*F(1,2)*F&
&(1,4)+12*F(1,1)*mu(1,1)*F(1,3)+tt11+tt17)*tt6+(12*F(1,1)*mu(1,1)*&
&F(1,4)+tt16+tt15)*F(2,1)*F(2,2)+(tt14+tt12+tt11)*tt4)
hes(1,2) = tt347
hes(1,3) = tt449
hes(1,4) = tt462
hes(1,5) = tt488
hes(1,6) = tt490
hes(1,7) = tt516
hes(1,8) = tt518
hes(2,1) = tt347
hes(2,2) = (-(12*F(1,1)*tt3*F(2,2)-12*tt1*F(2,1))*tt45*tt315)+2*t&
&t346**2*tt168*tt315-2*tt346*tt45*tt345+tt10*((tt519-6*tt5*lam(1,1&
&)*tt1*F(2,1)+2*tt313)*tt43+((12*F(1,1)*lam(1,1)*tt66*F(2,1)-12*tt&
&5*lam(1,1)*tt1*F(2,2))*F(2,3)+2*tt298*tt6+(6*tt293*F(2,1)+2*tt291&
&)*F(2,2)+12*tt281*tt4+6*tt279*F(2,1))*F(2,4)+(6*F(1,1)*lam(1,1)*t&
&t66*F(2,2)-6*lam(1,1)*tt180*F(2,1))*tt36+(2*tt257*tt6+(6*tt251*F(&
&2,1)+2*tt248)*F(2,2)+12*tt238*tt4+6*tt236*F(2,1))*F(2,3)+2*tt220*&
&tt8+(6*tt215*F(2,1)+2*tt212)*tt6+(12*tt202*tt4+6*tt198*F(2,1)+2*t&
&t193*tt3)*F(2,2)+20*tt181*tt2+12*tt179*tt4+6*tt176*tt1*F(2,1)+2*t&
&t175)
hes(2,3) = tt521
hes(2,4) = tt522
hes(2,5) = tt523
hes(2,6) = tt524
hes(2,7) = tt525
hes(2,8) = tt526
hes(3,1) = tt449
hes(3,2) = tt521
hes(3,3) = (-(12*F(1,1)*tt4*F(2,2)-12*F(1,2)*tt2)*tt45*tt315)+2*t&
&t448**2*tt168*tt315-2*tt448*tt45*tt447+tt10*((6*tt7*lam(1,1)*tt4*&
&F(2,2)+tt534+tt18+tt301)*tt43+((12*tt7*lam(1,1)*F(2,1)*tt6+(tt533&
&+tt41)*F(2,2)+24*F(1,1)*lam(1,1)*tt3*tt2+tt15)*F(2,3)+(36*tt5*lam&
&(1,1)*F(1,3)-48*tt7*lam(1,1))*tt4*tt6+(((-12*tt5*lam(1,1)*F(1,4))&
&-60*F(1,1)*lam(1,1)*F(1,2)*F(1,3)+96*tt5*lam(1,1)*F(1,2))*tt2+(tt&
&23+tt31)*tt4+2*F(1,1)*mu(1,1)*F(1,3)+tt531)*F(2,2)+(tt399+tt532+t&
&t375)*tt69+12*lam(1,1)*F(1,2)*F(1,3)*tt2+(tt37+tt29)*F(2,1))*F(2,&
&4)+(2*tt7*lam(1,1)*tt8+(tt21-18*tt5*lam(1,1)*F(1,2)*F(2,1))*tt6+3&
&6*F(1,1)*lam(1,1)*tt3*tt4*F(2,2)-20*lam(1,1)*tt1*tt2+tt527+tt301)&
&*tt36+((12*tt5*lam(1,1)*F(1,3)-32*tt7*lam(1,1))*F(2,1)*tt8+(((-36&
&*tt5*lam(1,1)*F(1,4))-36*F(1,1)*lam(1,1)*F(1,2)*F(1,3)+144*tt5*la&
&m(1,1)*F(1,2))*tt4+tt41)*tt6+((84*F(1,1)*lam(1,1)*F(1,2)*F(1,4)+t&
&t532-192*F(1,1)*lam(1,1)*tt3)*tt2+12*F(1,1)*lam(1,1)*F(1,4)*tt4+t&
&t30+tt16+tt29)*F(2,2)+(tt350+80*lam(1,1)*tt1)*tt69+(tt32+tt39)*tt&
&2+((-6*mu(1,1)*F(1,2)*F(1,4))+tt530-24*mu(1,1)*tt3+tt531)*F(2,1))&
&*F(2,3)+tt26*tt4*tt8+(tt24*tt2+tt18+tt528+tt530+tt301)*tt6+(tt20*&
&tt69+2*tt193*tt4+(tt38+12*mu(1,1)*F(1,2)*F(1,3)+tt15)*F(2,1))*F(2&
&,2)+tt19*tt182+tt529+6*tt176*F(1,2)*tt2+(12*mu(1,1)*F(1,2)*F(1,4)&
&+tt528+8*F(1,1)*mu(1,1)*F(1,3)+tt527+tt301)*tt4)
hes(3,4) = tt535
hes(3,5) = tt537
hes(3,6) = tt538
hes(3,7) = tt539
hes(3,8) = tt540
hes(4,1) = tt462
hes(4,2) = tt522
hes(4,3) = tt535
hes(4,4) = (-(12*tt7*F(2,2)-12*tt5*F(1,2)*F(2,1))*tt45*tt315)+2*t&
&t461**2*tt168*tt315-2*tt461*tt45*tt460+tt10*((6*tt225*lam(1,1)*F(&
&2,2)-6*tt90*lam(1,1)*F(1,2)*F(2,1))*tt43+((12*tt7*lam(1,1)*tt3*F(&
&2,1)-12*tt90*lam(1,1)*F(1,2)*F(2,2))*F(2,3)+12*tt307*tt6+6*tt305*&
&F(2,2)+2*tt299)*F(2,4)+(tt519+2*tt269)*tt36+(12*tt267*tt6+6*tt264&
&*F(2,2)+2*tt258)*F(2,3)+20*tt226*tt8+12*tt224*tt6+6*tt221*F(2,2)+&
&2*tt216)
hes(4,5) = tt541
hes(4,6) = tt542
hes(4,7) = tt543
hes(4,8) = tt544
hes(5,1) = tt488
hes(5,2) = tt523
hes(5,3) = tt537
hes(5,4) = tt541
hes(5,5) = tt10*(2*tt7*lam(1,1)*tt27+(tt21-6*tt5*lam(1,1)*F(1,2)*&
&F(2,1))*tt25+6*F(1,1)*lam(1,1)*tt3*tt4*tt8+(tt42+tt18+tt11)*tt6+2&
&*mu(1,1)*tt3*tt4)
hes(5,6) = tt545
hes(5,7) = tt546
hes(5,8) = tt547
hes(6,1) = tt490
hes(6,2) = tt524
hes(6,3) = tt538
hes(6,4) = tt542
hes(6,5) = tt545
hes(6,6) = 2*tt10*tt270
hes(6,7) = tt548
hes(6,8) = tt549
hes(7,1) = tt516
hes(7,2) = tt525
hes(7,3) = tt539
hes(7,4) = tt543
hes(7,5) = tt546
hes(7,6) = tt548
hes(7,7) = tt10*(2*tt7*lam(1,1)*tt4*tt8+(tt534+tt18+tt301)*tt6+6*&
&F(1,1)*lam(1,1)*tt3*tt69*F(2,2)-2*lam(1,1)*tt1*tt182+tt529+2*tt5*&
&mu(1,1)*tt4)
hes(7,8) = tt550
hes(8,1) = tt518
hes(8,2) = tt526
hes(8,3) = tt540
hes(8,4) = tt544
hes(8,5) = tt547
hes(8,6) = tt549
hes(8,7) = tt550
hes(8,8) = 2*tt10*tt314
END 
SUBROUTINE quad_stvk_F_val_at_qr(val, F2, mu, lam) 
IMPLICIT NONE 
REAL(KIND=8) val(1, 1) 
REAL(KIND=8) F2(2, 2) 
REAL(KIND=8) mu(1, 1) 
REAL(KIND=8) lam(1, 1) 
REAL(KIND=8)  tt1 
REAL(KIND=8)  tt2 
tt1 = F2(2,1)**2+F2(1,1)**2-1
tt2 = F2(2,2)**2+F2(1,2)**2-1
val(1,1) = (lam(1,1)*(tt2/2.0E+0+tt1/2.0E+0)**2)/2.0E+0+mu(1,1)*(&
&tt2**2/4.0E+0+(F2(2,1)*F2(2,2)+F2(1,1)*F2(1,2))**2/2.0E+0+tt1**2/&
&4.0E+0)
END 
SUBROUTINE quad_stvk_F_val_at_qr_jac(jac, F2, mu, lam) 
IMPLICIT NONE 
REAL(KIND=8) jac(1, 4) 
REAL(KIND=8) F2(2, 2) 
REAL(KIND=8) mu(1, 1) 
REAL(KIND=8) lam(1, 1) 
REAL(KIND=8)  tt1 
REAL(KIND=8)  tt2 
REAL(KIND=8)  tt3 
REAL(KIND=8)  tt4 
tt1 = F2(2,1)**2+F2(1,1)**2-1
tt2 = F2(2,1)*F2(2,2)+F2(1,1)*F2(1,2)
tt3 = F2(2,2)**2+F2(1,2)**2-1
tt4 = tt3/2.0E+0+tt1/2.0E+0
jac(1,1) = F2(1,1)*lam(1,1)*tt4+mu(1,1)*(F2(1,2)*tt2+F2(1,1)*tt1)&
&
jac(1,2) = lam(1,1)*F2(2,1)*tt4+mu(1,1)*(F2(2,2)*tt2+F2(2,1)*tt1)&
&
jac(1,3) = mu(1,1)*(F2(1,2)*tt3+F2(1,1)*tt2)+lam(1,1)*F2(1,2)*tt4&
&
jac(1,4) = mu(1,1)*(F2(2,2)*tt3+F2(2,1)*tt2)+lam(1,1)*F2(2,2)*tt4&
&
END 
SUBROUTINE quad_stvk_F_val_at_qr_hes(hes, F2, mu, lam) 
IMPLICIT NONE 
REAL(KIND=8) hes(4, 4) 
REAL(KIND=8) F2(2, 2) 
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
tt1 = F2(1,1)**2
tt2 = F2(1,2)**2
tt3 = F2(2,1)**2
tt4 = F2(2,2)**2
tt5 = lam(1,1)*((tt4+tt2-1)/2.0E+0+(tt3+tt1-1)/2.0E+0)
tt6 = mu(1,1)*(F2(1,2)*F2(2,2)+2*F2(1,1)*F2(2,1))+F2(1,1)*lam(1,1&
&)*F2(2,1)
tt7 = mu(1,1)*(F2(2,1)*F2(2,2)+2*F2(1,1)*F2(1,2))+F2(1,1)*lam(1,1&
&)*F2(1,2)
tt8 = F2(1,1)*lam(1,1)*F2(2,2)+mu(1,1)*F2(1,2)*F2(2,1)
tt9 = F2(1,1)*mu(1,1)*F2(2,2)+lam(1,1)*F2(1,2)*F2(2,1)
tt10 = mu(1,1)*(2*F2(2,1)*F2(2,2)+F2(1,1)*F2(1,2))+lam(1,1)*F2(2,&
&1)*F2(2,2)
tt11 = mu(1,1)*(2*F2(1,2)*F2(2,2)+F2(1,1)*F2(2,1))+lam(1,1)*F2(1,&
&2)*F2(2,2)
hes(1,1) = tt5+mu(1,1)*(tt3+tt2+3*tt1-1)+tt1*lam(1,1)
hes(1,2) = tt6
hes(1,3) = tt7
hes(1,4) = tt8
hes(2,1) = tt6
hes(2,2) = tt5+mu(1,1)*(tt4+3*tt3+tt1-1)+lam(1,1)*tt3
hes(2,3) = tt9
hes(2,4) = tt10
hes(3,1) = tt7
hes(3,2) = tt9
hes(3,3) = tt5+mu(1,1)*(tt4+3*tt2+tt1-1)+lam(1,1)*tt2
hes(3,4) = tt11
hes(4,1) = tt8
hes(4,2) = tt10
hes(4,3) = tt11
hes(4,4) = tt5+mu(1,1)*(3*tt4+tt3+tt2-1)+lam(1,1)*tt4
END 
SUBROUTINE quad_stvk_q2l_psi(val, F, mu, lam) 
IMPLICIT NONE 
REAL(KIND=8) val(1, 1) 
REAL(KIND=8) F(2, 4) 
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
tt1 = F(1,1)**4
tt2 = F(1,1)**2
tt3 = 6*tt2*lam(1,1)
tt4 = 12*tt2*mu(1,1)
tt5 = F(1,2)**2
tt6 = 6*mu(1,1)+3*lam(1,1)
tt7 = F(1,1)**3
tt8 = -16*F(1,1)*mu(1,1)
tt9 = tt8-8*F(1,1)*lam(1,1)
tt10 = 4*mu(1,1)+2*lam(1,1)
tt11 = tt10*tt5
tt12 = F(1,3)**2
tt13 = -8*tt2*lam(1,1)
tt14 = -16*tt2*mu(1,1)
tt15 = (-16*mu(1,1))-8*lam(1,1)
tt16 = 16*F(1,1)*mu(1,1)+8*F(1,1)*lam(1,1)
tt17 = (2*tt2-4)*lam(1,1)
tt18 = (4*tt2-4)*mu(1,1)
tt19 = 12*mu(1,1)+6*lam(1,1)
tt20 = tt19*tt5
tt21 = F(1,4)**2
tt22 = F(2,1)**2
tt23 = tt19*tt22
tt24 = F(2,2)**2
tt25 = 8*F(1,1)*mu(1,1)*F(1,4)+8*mu(1,1)*F(1,2)*F(1,3)-16*F(1,1)*&
&mu(1,1)*F(1,2)
tt26 = 16*mu(1,1)+8*lam(1,1)
val(1,1) = ((tt19*tt24+tt10*tt22+tt11-4*mu(1,1)+tt17)*F(2,4)**2+(&
&(tt26*F(2,1)*F(2,2)+8*F(1,1)*mu(1,1)*F(1,2))*F(2,3)+tt15*F(2,2)**&
&3+(tt15*tt22+tt26*F(1,2)*F(1,4)+8*F(1,1)*lam(1,1)*F(1,3)+tt15*tt5&
&+tt13)*F(2,2)+tt25*F(2,1))*F(2,4)+(tt10*tt24+tt23+2*lam(1,1)*tt5+&
&tt18+tt17)*F(2,3)**2+(tt15*F(2,1)*tt24+tt25*F(2,2)+tt15*F(2,1)**3&
&+(8*lam(1,1)*F(1,2)*F(1,4)+tt16*F(1,3)-8*lam(1,1)*tt5+tt14+tt13)*&
&F(2,1))*F(2,3)+tt6*F(2,2)**4+(tt23+tt10*tt21+tt15*F(1,2)*F(1,4)+2&
&*lam(1,1)*tt12-8*F(1,1)*lam(1,1)*F(1,3)+tt20+tt3)*tt24+((8*mu(1,1&
&)*F(1,3)+tt8)*F(1,4)-16*mu(1,1)*F(1,2)*F(1,3)+24*F(1,1)*mu(1,1)*F&
&(1,2))*F(2,1)*F(2,2)+tt6*F(2,1)**4+(2*lam(1,1)*tt21-8*lam(1,1)*F(&
&1,2)*F(1,4)+tt10*tt12+tt9*F(1,3)+6*lam(1,1)*tt5+tt4+tt3)*tt22+(tt&
&20+tt18+tt17)*tt21+(tt16*F(1,2)*F(1,3)+tt15*F(1,2)**3+(tt14+tt13)&
&*F(1,2))*F(1,4)+(tt11+(12*tt2-4)*mu(1,1)+(6*tt2-4)*lam(1,1))*tt12&
&+(tt9*tt5-16*tt7*mu(1,1)-8*tt7*lam(1,1))*F(1,3)+tt6*F(1,2)**4+(tt&
&4+tt3)*tt5+(6*tt1+4)*mu(1,1)+(3*tt1+4)*lam(1,1))/8.0E+0
END 
SUBROUTINE quad_stvk_q2l_psi_jac(jac, F, mu, lam) 
IMPLICIT NONE 
REAL(KIND=8) jac(1, 8) 
REAL(KIND=8) F(2, 4) 
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
tt1 = F(1,1)**3
tt2 = 12*F(1,1)*lam(1,1)
tt3 = 24*F(1,1)*mu(1,1)
tt4 = tt3+tt2
tt5 = F(1,2)**2
tt6 = F(1,1)**2
tt7 = (-16*mu(1,1))-8*lam(1,1)
tt8 = tt7*tt5
tt9 = F(1,3)**2
tt10 = -16*F(1,1)*lam(1,1)
tt11 = -32*F(1,1)*mu(1,1)
tt12 = 16*mu(1,1)+8*lam(1,1)
tt13 = 8*F(1,1)*mu(1,1)+4*F(1,1)*lam(1,1)
tt14 = F(1,4)**2
tt15 = F(2,1)**2
tt16 = F(2,2)**2
tt17 = 8*mu(1,1)*F(1,4)-16*mu(1,1)*F(1,2)
tt18 = F(2,3)**2
tt19 = F(2,4)**2
tt20 = 6*tt6*lam(1,1)
tt21 = 12*tt6*mu(1,1)
tt22 = -8*F(1,1)*lam(1,1)
tt23 = -16*F(1,1)*mu(1,1)
tt24 = tt23+tt22
tt25 = 4*mu(1,1)+2*lam(1,1)
tt26 = 6*mu(1,1)+3*lam(1,1)
tt27 = F(2,1)**3
tt28 = 8*mu(1,1)*F(1,3)+tt23
tt29 = tt28*F(1,4)-16*mu(1,1)*F(1,2)*F(1,3)+24*F(1,1)*mu(1,1)*F(1&
&,2)
tt30 = 12*mu(1,1)+6*lam(1,1)
tt31 = -8*tt6*lam(1,1)
tt32 = -16*tt6*mu(1,1)
tt33 = -8*lam(1,1)*tt5
tt34 = 16*F(1,1)*mu(1,1)+8*F(1,1)*lam(1,1)
tt35 = tt34*F(1,3)
tt36 = 8*lam(1,1)*F(1,2)*F(1,4)
tt37 = -16*F(1,1)*mu(1,1)*F(1,2)
tt38 = 8*mu(1,1)*F(1,2)*F(1,3)
tt39 = 8*F(1,1)*mu(1,1)*F(1,4)
tt40 = 2*tt7*F(2,1)*F(2,2)
tt41 = F(1,2)**3
tt42 = tt30*tt5
tt43 = tt30*tt15
tt44 = F(2,2)**3
tt45 = 8*F(1,1)*lam(1,1)*F(1,3)
tt46 = tt12*F(1,2)*F(1,4)
tt47 = tt7*tt15
tt48 = tt25*tt5
tt49 = tt39+tt38+tt37
tt50 = (2*tt6-4)*lam(1,1)
tt51 = (4*tt6-4)*mu(1,1)
tt52 = tt12*F(2,1)*F(2,2)+8*F(1,1)*mu(1,1)*F(1,2)
jac(1,1) = (4*F(1,1)*lam(1,1)*tt19+(8*mu(1,1)*F(1,2)*F(2,3)+(8*la&
&m(1,1)*F(1,3)+tt10)*F(2,2)+tt17*F(2,1))*F(2,4)+tt13*tt18+(tt17*F(&
&2,2)+(tt12*F(1,3)+tt11+tt10)*F(2,1))*F(2,3)+(tt2-8*lam(1,1)*F(1,3&
&))*tt16+(24*mu(1,1)*F(1,2)-16*mu(1,1)*F(1,4))*F(2,1)*F(2,2)+(tt7*&
&F(1,3)+tt3+tt2)*tt15+tt13*tt14+(tt12*F(1,2)*F(1,3)+(tt11+tt10)*F(&
&1,2))*F(1,4)+tt4*tt9+(tt8-48*tt6*mu(1,1)-24*tt6*lam(1,1))*F(1,3)+&
&tt4*tt5+24*tt1*mu(1,1)+12*tt1*lam(1,1))/8.0E+0
jac(1,2) = (2*tt25*F(2,1)*tt19+(tt12*F(2,2)*F(2,3)+tt40+tt39+tt38&
&+tt37)*F(2,4)+2*tt30*F(2,1)*tt18+(tt7*tt16+3*tt7*tt15+tt36+tt35+t&
&t33+tt32+tt31)*F(2,3)+2*tt30*F(2,1)*tt16+tt29*F(2,2)+4*tt26*tt27+&
&2*(2*lam(1,1)*tt14-8*lam(1,1)*F(1,2)*F(1,4)+tt25*tt9+tt24*F(1,3)+&
&6*lam(1,1)*tt5+tt21+tt20)*F(2,1))/8.0E+0
jac(1,3) = (2*tt25*F(1,2)*tt19+(8*F(1,1)*mu(1,1)*F(2,3)+(tt12*F(1&
&,4)+2*tt7*F(1,2))*F(2,2)+tt28*F(2,1))*F(2,4)+4*lam(1,1)*F(1,2)*tt&
&18+(tt28*F(2,2)+(8*lam(1,1)*F(1,4)-16*lam(1,1)*F(1,2))*F(2,1))*F(&
&2,3)+(tt7*F(1,4)+2*tt30*F(1,2))*tt16+(tt3-16*mu(1,1)*F(1,3))*F(2,&
&1)*F(2,2)+(12*lam(1,1)*F(1,2)-8*lam(1,1)*F(1,4))*tt15+2*tt30*F(1,&
&2)*tt14+(tt35+3*tt7*tt5+tt32+tt31)*F(1,4)+2*tt25*F(1,2)*tt9+2*tt2&
&4*F(1,2)*F(1,3)+4*tt26*tt41+2*(tt21+tt20)*F(1,2))/8.0E+0
jac(1,4) = (2*tt30*F(2,2)*tt19+(tt12*F(2,1)*F(2,3)+3*tt7*tt16+tt4&
&7+tt46+tt45+tt8+tt31)*F(2,4)+2*tt25*F(2,2)*tt18+(tt40+tt39+tt38+t&
&t37)*F(2,3)+4*tt26*tt44+2*(tt43+tt25*tt14+tt7*F(1,2)*F(1,4)+2*lam&
&(1,1)*tt9-8*F(1,1)*lam(1,1)*F(1,3)+tt42+tt20)*F(2,2)+tt29*F(2,1))&
&/8.0E+0
jac(1,5) = ((8*F(1,1)*lam(1,1)*F(2,2)+8*mu(1,1)*F(1,2)*F(2,1))*F(&
&2,4)+(8*mu(1,1)*F(1,2)*F(2,2)+tt34*F(2,1))*F(2,3)+(4*lam(1,1)*F(1&
&,3)+tt22)*tt16+tt17*F(2,1)*F(2,2)+(2*tt25*F(1,3)+tt23+tt22)*tt15+&
&tt34*F(1,2)*F(1,4)+2*(tt48+(12*tt6-4)*mu(1,1)+(6*tt6-4)*lam(1,1))&
&*F(1,3)+tt24*tt5-16*tt1*mu(1,1)-8*tt1*lam(1,1))/8.0E+0
jac(1,6) = (tt52*F(2,4)+2*(tt25*tt16+tt43+2*lam(1,1)*tt5+tt51+tt5&
&0)*F(2,3)+tt7*F(2,1)*tt16+tt49*F(2,2)+tt7*tt27+(tt36+tt35+tt33+tt&
&32+tt31)*F(2,1))/8.0E+0
jac(1,7) = ((tt12*F(1,2)*F(2,2)+8*F(1,1)*mu(1,1)*F(2,1))*F(2,4)+(&
&8*F(1,1)*mu(1,1)*F(2,2)+8*lam(1,1)*F(1,2)*F(2,1))*F(2,3)+(2*tt25*&
&F(1,4)+tt7*F(1,2))*tt16+tt28*F(2,1)*F(2,2)+(4*lam(1,1)*F(1,4)-8*l&
&am(1,1)*F(1,2))*tt15+2*(tt42+tt51+tt50)*F(1,4)+tt34*F(1,2)*F(1,3)&
&+tt7*tt41+(tt32+tt31)*F(1,2))/8.0E+0
jac(1,8) = (2*(tt30*tt16+tt25*tt15+tt48-4*mu(1,1)+tt50)*F(2,4)+tt&
&52*F(2,3)+tt7*tt44+(tt47+tt46+tt45+tt8+tt31)*F(2,2)+tt49*F(2,1))/&
&8.0E+0
END 
SUBROUTINE quad_stvk_q2l_psi_hes(hes, F, mu, lam) 
IMPLICIT NONE 
REAL(KIND=8) hes(8, 8) 
REAL(KIND=8) F(2, 4) 
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
tt1 = F(1,1)**2
tt2 = 24*mu(1,1)+12*lam(1,1)
tt3 = F(1,2)**2
tt4 = F(1,3)**2
tt5 = (-32*mu(1,1))-16*lam(1,1)
tt6 = 8*mu(1,1)+4*lam(1,1)
tt7 = F(1,4)**2
tt8 = F(2,1)**2
tt9 = F(2,2)**2
tt10 = F(2,3)**2
tt11 = F(2,4)**2
tt12 = 12*F(1,1)*lam(1,1)
tt13 = 24*F(1,1)*mu(1,1)
tt14 = (-16*mu(1,1))-8*lam(1,1)
tt15 = 24*mu(1,1)*F(1,2)-16*mu(1,1)*F(1,4)
tt16 = -16*F(1,1)*lam(1,1)
tt17 = -32*F(1,1)*mu(1,1)
tt18 = 16*mu(1,1)+8*lam(1,1)
tt19 = tt18*F(1,3)+tt17+tt16
tt20 = 8*mu(1,1)*F(1,4)-16*mu(1,1)*F(1,2)
tt21 = (tt20*F(2,4)+tt19*F(2,3)+tt15*F(2,2)+2*(tt14*F(1,3)+tt13+t&
&t12)*F(2,1))/8.0E+0
tt22 = tt13+tt12
tt23 = ((8*mu(1,1)*F(2,3)-16*mu(1,1)*F(2,1))*F(2,4)-16*mu(1,1)*F(&
&2,2)*F(2,3)+24*mu(1,1)*F(2,1)*F(2,2)+tt19*F(1,4)+2*tt14*F(1,2)*F(&
&1,3)+2*tt22*F(1,2))/8.0E+0
tt24 = 8*lam(1,1)*F(1,3)+tt16
tt25 = (tt24*F(2,4)+tt20*F(2,3)+2*(tt12-8*lam(1,1)*F(1,3))*F(2,2)&
&+tt15*F(2,1))/8.0E+0
tt26 = tt14*tt3
tt27 = tt18*F(1,2)*F(1,4)
tt28 = tt14*tt8
tt29 = tt18*F(2,1)*F(2,3)
tt30 = (8*lam(1,1)*F(2,2)*F(2,4)+tt29-8*lam(1,1)*tt9+tt28+tt27+2*&
&tt22*F(1,3)+tt26-48*tt1*mu(1,1)-24*tt1*lam(1,1))/8.0E+0
tt31 = tt20*F(2,2)
tt32 = 8*F(1,1)*mu(1,1)+4*F(1,1)*lam(1,1)
tt33 = 8*mu(1,1)*F(1,2)*F(2,4)
tt34 = (tt33+2*tt32*F(2,3)+tt31+tt19*F(2,1))/8.0E+0
tt35 = -16*mu(1,1)*F(2,1)*F(2,2)
tt36 = 8*mu(1,1)*F(2,2)*F(2,3)
tt37 = 8*mu(1,1)*F(2,1)*F(2,4)
tt38 = (tt37+tt36+tt35+2*tt32*F(1,4)+tt18*F(1,2)*F(1,3)+(tt17+tt1&
&6)*F(1,2))/8.0E+0
tt39 = tt20*F(2,1)
tt40 = 8*mu(1,1)*F(1,2)*F(2,3)
tt41 = 8*F(1,1)*lam(1,1)*F(2,4)
tt42 = (tt41+tt40+tt24*F(2,2)+tt39)/8.0E+0
tt43 = 6*tt1*lam(1,1)
tt44 = 12*tt1*mu(1,1)
tt45 = -8*F(1,1)*lam(1,1)
tt46 = -16*F(1,1)*mu(1,1)
tt47 = tt46+tt45
tt48 = 4*mu(1,1)+2*lam(1,1)
tt49 = 6*mu(1,1)+3*lam(1,1)
tt50 = 12*mu(1,1)+6*lam(1,1)
tt51 = 2*tt50*tt9
tt52 = 2*tt14*F(2,2)*F(2,4)
tt53 = 2*tt48*tt11
tt54 = tt13-16*mu(1,1)*F(1,3)
tt55 = 8*lam(1,1)*F(1,4)-16*lam(1,1)*F(1,2)
tt56 = 8*mu(1,1)*F(1,3)+tt46
tt57 = (tt56*F(2,4)+tt55*F(2,3)+tt54*F(2,2)+2*(12*lam(1,1)*F(1,2)&
&-8*lam(1,1)*F(1,4))*F(2,1))/8.0E+0
tt58 = ((tt18*F(2,3)+2*tt14*F(2,1))*F(2,4)+2*tt14*F(2,2)*F(2,3)+4&
&*tt50*F(2,1)*F(2,2)+tt56*F(1,4)-16*mu(1,1)*F(1,2)*F(1,3)+24*F(1,1&
&)*mu(1,1)*F(1,2))/8.0E+0
tt59 = 16*F(1,1)*mu(1,1)+8*F(1,1)*lam(1,1)
tt60 = (tt33+tt59*F(2,3)+tt31+2*(2*tt48*F(1,3)+tt46+tt45)*F(2,1))&
&/8.0E+0
tt61 = -8*tt1*lam(1,1)
tt62 = -16*tt1*mu(1,1)
tt63 = tt59*F(1,3)
tt64 = tt14*tt9
tt65 = tt18*F(2,2)*F(2,4)
tt66 = (tt65+4*tt50*F(2,1)*F(2,3)+tt64+3*tt14*tt8+8*lam(1,1)*F(1,&
&2)*F(1,4)+tt63-8*lam(1,1)*tt3+tt62+tt61)/8.0E+0
tt67 = tt56*F(2,2)
tt68 = 8*lam(1,1)*F(1,2)*F(2,3)
tt69 = 8*F(1,1)*mu(1,1)*F(2,4)
tt70 = (tt69+tt68+tt67+2*(4*lam(1,1)*F(1,4)-8*lam(1,1)*F(1,2))*F(&
&2,1))/8.0E+0
tt71 = -16*F(1,1)*mu(1,1)*F(1,2)
tt72 = 8*mu(1,1)*F(1,2)*F(1,3)
tt73 = 8*F(1,1)*mu(1,1)*F(1,4)
tt74 = 2*tt14*F(2,1)*F(2,2)
tt75 = (4*tt48*F(2,1)*F(2,4)+tt18*F(2,2)*F(2,3)+tt74+tt73+tt72+tt&
&71)/8.0E+0
tt76 = tt18*F(1,4)+2*tt14*F(1,2)
tt77 = (tt76*F(2,4)+tt56*F(2,3)+2*(tt14*F(1,4)+2*tt50*F(1,2))*F(2&
&,2)+tt54*F(2,1))/8.0E+0
tt78 = (tt37+tt36+tt35+tt59*F(1,4)+4*tt48*F(1,2)*F(1,3)+2*tt47*F(&
&1,2))/8.0E+0
tt79 = (tt69+tt68+tt67+tt55*F(2,1))/8.0E+0
tt80 = (tt65+8*lam(1,1)*F(2,1)*F(2,3)+tt64-8*lam(1,1)*tt8+4*tt50*&
&F(1,2)*F(1,4)+tt63+3*tt14*tt3+tt62+tt61)/8.0E+0
tt81 = tt56*F(2,1)
tt82 = 8*F(1,1)*mu(1,1)*F(2,3)
tt83 = (4*tt48*F(1,2)*F(2,4)+tt82+tt76*F(2,2)+tt81)/8.0E+0
tt84 = tt50*tt3
tt85 = tt50*tt8
tt86 = (tt41+tt40+2*(4*lam(1,1)*F(1,3)+tt45)*F(2,2)+tt39)/8.0E+0
tt87 = (tt18*F(2,1)*F(2,4)+4*tt48*F(2,2)*F(2,3)+tt74+tt73+tt72+tt&
&71)/8.0E+0
tt88 = (tt18*F(1,2)*F(2,4)+tt82+2*(2*tt48*F(1,4)+tt14*F(1,2))*F(2&
&,2)+tt81)/8.0E+0
tt89 = (4*tt50*F(2,2)*F(2,4)+tt29+3*tt14*tt9+tt28+tt27+8*F(1,1)*l&
&am(1,1)*F(1,3)+tt26+tt61)/8.0E+0
tt90 = tt48*tt3
tt91 = (8*mu(1,1)*F(1,2)*F(2,2)+tt59*F(2,1))/8.0E+0
tt92 = (8*mu(1,1)*F(2,1)*F(2,2)+tt59*F(1,2))/8.0E+0
tt93 = (8*F(1,1)*lam(1,1)*F(2,2)+8*mu(1,1)*F(1,2)*F(2,1))/8.0E+0
tt94 = (2*tt1-4)*lam(1,1)
tt95 = (4*tt1-4)*mu(1,1)
tt96 = (8*F(1,1)*mu(1,1)*F(2,2)+8*lam(1,1)*F(1,2)*F(2,1))/8.0E+0
tt97 = (tt18*F(2,1)*F(2,2)+8*F(1,1)*mu(1,1)*F(1,2))/8.0E+0
tt98 = (tt18*F(1,2)*F(2,2)+8*F(1,1)*mu(1,1)*F(2,1))/8.0E+0
hes(1,1) = (4*lam(1,1)*tt11-16*lam(1,1)*F(2,2)*F(2,4)+tt6*tt10+tt&
&5*F(2,1)*F(2,3)+12*lam(1,1)*tt9+tt2*tt8+tt6*tt7+tt5*F(1,2)*F(1,4)&
&+tt2*tt4+((-96*F(1,1)*mu(1,1))-48*F(1,1)*lam(1,1))*F(1,3)+tt2*tt3&
&+72*tt1*mu(1,1)+36*tt1*lam(1,1))/8.0E+0
hes(1,2) = tt21
hes(1,3) = tt23
hes(1,4) = tt25
hes(1,5) = tt30
hes(1,6) = tt34
hes(1,7) = tt38
hes(1,8) = tt42
hes(2,1) = tt21
hes(2,2) = (tt53+tt52+2*tt50*tt10+6*tt14*F(2,1)*F(2,3)+tt51+12*tt&
&49*tt8+2*(2*lam(1,1)*tt7-8*lam(1,1)*F(1,2)*F(1,4)+tt48*tt4+tt47*F&
&(1,3)+6*lam(1,1)*tt3+tt44+tt43))/8.0E+0
hes(2,3) = tt57
hes(2,4) = tt58
hes(2,5) = tt60
hes(2,6) = tt66
hes(2,7) = tt70
hes(2,8) = tt75
hes(3,1) = tt23
hes(3,2) = tt57
hes(3,3) = (tt53+tt52+4*lam(1,1)*tt10-16*lam(1,1)*F(2,1)*F(2,3)+t&
&t51+12*lam(1,1)*tt8+2*tt50*tt7+6*tt14*F(1,2)*F(1,4)+2*tt48*tt4+2*&
&tt47*F(1,3)+12*tt49*tt3+2*(tt44+tt43))/8.0E+0
hes(3,4) = tt77
hes(3,5) = tt78
hes(3,6) = tt79
hes(3,7) = tt80
hes(3,8) = tt83
hes(4,1) = tt25
hes(4,2) = tt58
hes(4,3) = tt77
hes(4,4) = (2*tt50*tt11+6*tt14*F(2,2)*F(2,4)+2*tt48*tt10+2*tt14*F&
&(2,1)*F(2,3)+12*tt49*tt9+2*(tt85+tt48*tt7+tt14*F(1,2)*F(1,4)+2*la&
&m(1,1)*tt4-8*F(1,1)*lam(1,1)*F(1,3)+tt84+tt43))/8.0E+0
hes(4,5) = tt86
hes(4,6) = tt87
hes(4,7) = tt88
hes(4,8) = tt89
hes(5,1) = tt30
hes(5,2) = tt60
hes(5,3) = tt78
hes(5,4) = tt86
hes(5,5) = (4*lam(1,1)*tt9+2*tt48*tt8+2*(tt90+(12*tt1-4)*mu(1,1)+&
&(6*tt1-4)*lam(1,1)))/8.0E+0
hes(5,6) = tt91
hes(5,7) = tt92
hes(5,8) = tt93
hes(6,1) = tt34
hes(6,2) = tt66
hes(6,3) = tt79
hes(6,4) = tt87
hes(6,5) = tt91
hes(6,6) = (tt48*tt9+tt85+2*lam(1,1)*tt3+tt95+tt94)/4.0E+0
hes(6,7) = tt96
hes(6,8) = tt97
hes(7,1) = tt38
hes(7,2) = tt70
hes(7,3) = tt80
hes(7,4) = tt88
hes(7,5) = tt92
hes(7,6) = tt96
hes(7,7) = (2*tt48*tt9+4*lam(1,1)*tt8+2*(tt84+tt95+tt94))/8.0E+0
hes(7,8) = tt98
hes(8,1) = tt42
hes(8,2) = tt75
hes(8,3) = tt83
hes(8,4) = tt89
hes(8,5) = tt93
hes(8,6) = tt97
hes(8,7) = tt98
hes(8,8) = (tt50*tt9+tt48*tt8+tt90-4*mu(1,1)+tt94)/4.0E+0
END 
