SUBROUTINE hex_lin_F_at_qr(val, F, mu, lam) 
IMPLICIT NONE 
REAL(KIND=8) val(1, 1) 
REAL(KIND=8) F(3, 3) 
REAL(KIND=8) mu(1, 1) 
REAL(KIND=8) lam(1, 1) 
val(1,1) = (lam(1,1)*(F(3,3)+F(2,2)+F(1,1)-3)**2)/2.0E+0+mu(1,1)*&
&((F(3,3)-1)**2+(F(3,2)+F(2,3))**2/2.0E+0+(F(3,1)+F(1,3))**2/2.0E+&
&0+(F(2,2)-1)**2+(F(2,1)+F(1,2))**2/2.0E+0+(F(1,1)-1)**2)
END 
SUBROUTINE hex_lin_F_at_qr_jac(jac, F, mu, lam) 
IMPLICIT NONE 
REAL(KIND=8) jac(1, 9) 
REAL(KIND=8) F(3, 3) 
REAL(KIND=8) mu(1, 1) 
REAL(KIND=8) lam(1, 1) 
REAL(KIND=8)  tt1 
REAL(KIND=8)  tt2 
REAL(KIND=8)  tt3 
REAL(KIND=8)  tt4 
tt1 = lam(1,1)*(F(3,3)+F(2,2)+F(1,1)-3)
tt2 = mu(1,1)*(F(2,1)+F(1,2))
tt3 = mu(1,1)*(F(3,1)+F(1,3))
tt4 = mu(1,1)*(F(3,2)+F(2,3))
jac(1,1) = tt1+2*(F(1,1)-1)*mu(1,1)
jac(1,2) = tt2
jac(1,3) = tt3
jac(1,4) = tt2
jac(1,5) = tt1+2*mu(1,1)*(F(2,2)-1)
jac(1,6) = tt4
jac(1,7) = tt3
jac(1,8) = tt4
jac(1,9) = tt1+2*mu(1,1)*(F(3,3)-1)
END 
SUBROUTINE hex_coro_F_at_qr(val, F, R, mu, lam) 
IMPLICIT NONE 
REAL(KIND=8) val(1, 1) 
REAL(KIND=8) F(3, 3) 
REAL(KIND=8) R(3, 3) 
REAL(KIND=8) mu(1, 1) 
REAL(KIND=8) lam(1, 1) 
val(1,1) = (lam(1,1)*((2*F(3,3)*R(3,3)+2*F(2,3)*R(2,3)+2*F(1,3)*R&
&(1,3))/2.0E+0+(2*F(3,2)*R(3,2)+2*F(2,2)*R(2,2)+2*F(1,2)*R(1,2))/2&
&.0E+0+(2*F(3,1)*R(3,1)+2*F(2,1)*R(2,1)+2*F(1,1)*R(1,1))/2.0E+0-3)&
&**2)/2.0E+0+mu(1,1)*((F(3,3)*R(3,3)+F(2,3)*R(2,3)+F(1,3)*R(1,3)-1&
&)**2+(F(3,2)*R(3,3)+R(3,2)*F(3,3)+F(2,2)*R(2,3)+R(2,2)*F(2,3)+F(1&
&,2)*R(1,3)+R(1,2)*F(1,3))**2/2.0E+0+(F(3,1)*R(3,3)+R(3,1)*F(3,3)+&
&F(2,1)*R(2,3)+R(2,1)*F(2,3)+F(1,1)*R(1,3)+R(1,1)*F(1,3))**2/2.0E+&
&0+(F(3,2)*R(3,2)+F(2,2)*R(2,2)+F(1,2)*R(1,2)-1)**2+(F(3,1)*R(3,2)&
&+R(3,1)*F(3,2)+F(2,1)*R(2,2)+R(2,1)*F(2,2)+F(1,1)*R(1,2)+R(1,1)*F&
&(1,2))**2/2.0E+0+(F(3,1)*R(3,1)+F(2,1)*R(2,1)+F(1,1)*R(1,1)-1)**2&
&)
END 
SUBROUTINE hex_coro_F_at_qr_jac(jac, F, R, mu, lam) 
IMPLICIT NONE 
REAL(KIND=8) jac(1, 9) 
REAL(KIND=8) F(3, 3) 
REAL(KIND=8) R(3, 3) 
REAL(KIND=8) mu(1, 1) 
REAL(KIND=8) lam(1, 1) 
REAL(KIND=8)  tt1 
REAL(KIND=8)  tt2 
REAL(KIND=8)  tt3 
REAL(KIND=8)  tt4 
REAL(KIND=8)  tt5 
REAL(KIND=8)  tt6 
REAL(KIND=8)  tt7 
tt1 = F(3,1)*R(3,1)+F(2,1)*R(2,1)+F(1,1)*R(1,1)-1
tt2 = F(3,1)*R(3,2)+R(3,1)*F(3,2)+F(2,1)*R(2,2)+R(2,1)*F(2,2)+F(1&
&,1)*R(1,2)+R(1,1)*F(1,2)
tt3 = F(3,1)*R(3,3)+R(3,1)*F(3,3)+F(2,1)*R(2,3)+R(2,1)*F(2,3)+F(1&
&,1)*R(1,3)+R(1,1)*F(1,3)
tt4 = (2*F(3,3)*R(3,3)+2*F(2,3)*R(2,3)+2*F(1,3)*R(1,3))/2.0E+0+(2&
&*F(3,2)*R(3,2)+2*F(2,2)*R(2,2)+2*F(1,2)*R(1,2))/2.0E+0+(2*F(3,1)*&
&R(3,1)+2*F(2,1)*R(2,1)+2*F(1,1)*R(1,1))/2.0E+0-3
tt5 = F(3,2)*R(3,2)+F(2,2)*R(2,2)+F(1,2)*R(1,2)-1
tt6 = F(3,2)*R(3,3)+R(3,2)*F(3,3)+F(2,2)*R(2,3)+R(2,2)*F(2,3)+F(1&
&,2)*R(1,3)+R(1,2)*F(1,3)
tt7 = F(3,3)*R(3,3)+F(2,3)*R(2,3)+F(1,3)*R(1,3)-1
jac(1,1) = R(1,1)*lam(1,1)*tt4+mu(1,1)*(R(1,3)*tt3+R(1,2)*tt2+2*R&
&(1,1)*tt1)
jac(1,2) = lam(1,1)*R(2,1)*tt4+mu(1,1)*(R(2,3)*tt3+R(2,2)*tt2+2*R&
&(2,1)*tt1)
jac(1,3) = lam(1,1)*R(3,1)*tt4+mu(1,1)*(R(3,3)*tt3+R(3,2)*tt2+2*R&
&(3,1)*tt1)
jac(1,4) = lam(1,1)*R(1,2)*tt4+mu(1,1)*(R(1,3)*tt6+2*R(1,2)*tt5+R&
&(1,1)*tt2)
jac(1,5) = lam(1,1)*R(2,2)*tt4+mu(1,1)*(R(2,3)*tt6+2*R(2,2)*tt5+R&
&(2,1)*tt2)
jac(1,6) = lam(1,1)*R(3,2)*tt4+mu(1,1)*(R(3,3)*tt6+2*R(3,2)*tt5+R&
&(3,1)*tt2)
jac(1,7) = lam(1,1)*R(1,3)*tt4+mu(1,1)*(2*R(1,3)*tt7+R(1,2)*tt6+R&
&(1,1)*tt3)
jac(1,8) = lam(1,1)*R(2,3)*tt4+mu(1,1)*(2*R(2,3)*tt7+R(2,2)*tt6+R&
&(2,1)*tt3)
jac(1,9) = lam(1,1)*R(3,3)*tt4+mu(1,1)*(2*R(3,3)*tt7+R(3,2)*tt6+R&
&(3,1)*tt3)
END 
SUBROUTINE hex_coro_F_at_qr_hes(hes, F, R, mu, lam) 
IMPLICIT NONE 
REAL(KIND=8) hes(9, 9) 
REAL(KIND=8) F(3, 3) 
REAL(KIND=8) R(3, 3) 
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
tt1 = R(1,1)**2
tt2 = R(1,2)**2
tt3 = R(1,3)**2
tt4 = R(1,2)*R(2,2)
tt5 = R(1,3)*R(2,3)
tt6 = mu(1,1)*(tt5+tt4+2*R(1,1)*R(2,1))+R(1,1)*lam(1,1)*R(2,1)
tt7 = R(1,2)*R(3,2)
tt8 = R(1,3)*R(3,3)
tt9 = mu(1,1)*(tt8+tt7+2*R(1,1)*R(3,1))+R(1,1)*lam(1,1)*R(3,1)
tt10 = R(1,1)*mu(1,1)*R(1,2)+R(1,1)*lam(1,1)*R(1,2)
tt11 = R(1,1)*lam(1,1)*R(2,2)+mu(1,1)*R(1,2)*R(2,1)
tt12 = R(1,1)*lam(1,1)*R(3,2)+mu(1,1)*R(1,2)*R(3,1)
tt13 = R(1,1)*mu(1,1)*R(1,3)+R(1,1)*lam(1,1)*R(1,3)
tt14 = R(1,1)*lam(1,1)*R(2,3)+mu(1,1)*R(1,3)*R(2,1)
tt15 = R(1,1)*lam(1,1)*R(3,3)+mu(1,1)*R(1,3)*R(3,1)
tt16 = R(2,1)**2
tt17 = R(2,2)**2
tt18 = R(2,3)**2
tt19 = R(2,2)*R(3,2)
tt20 = R(2,3)*R(3,3)
tt21 = mu(1,1)*(tt20+tt19+2*R(2,1)*R(3,1))+lam(1,1)*R(2,1)*R(3,1)&
&
tt22 = R(1,1)*mu(1,1)*R(2,2)+lam(1,1)*R(1,2)*R(2,1)
tt23 = mu(1,1)*R(2,1)*R(2,2)+lam(1,1)*R(2,1)*R(2,2)
tt24 = lam(1,1)*R(2,1)*R(3,2)+mu(1,1)*R(2,2)*R(3,1)
tt25 = R(1,1)*mu(1,1)*R(2,3)+lam(1,1)*R(1,3)*R(2,1)
tt26 = mu(1,1)*R(2,1)*R(2,3)+lam(1,1)*R(2,1)*R(2,3)
tt27 = lam(1,1)*R(2,1)*R(3,3)+mu(1,1)*R(2,3)*R(3,1)
tt28 = R(3,1)**2
tt29 = R(3,2)**2
tt30 = R(3,3)**2
tt31 = R(1,1)*mu(1,1)*R(3,2)+lam(1,1)*R(1,2)*R(3,1)
tt32 = mu(1,1)*R(2,1)*R(3,2)+lam(1,1)*R(2,2)*R(3,1)
tt33 = mu(1,1)*R(3,1)*R(3,2)+lam(1,1)*R(3,1)*R(3,2)
tt34 = R(1,1)*mu(1,1)*R(3,3)+lam(1,1)*R(1,3)*R(3,1)
tt35 = mu(1,1)*R(2,1)*R(3,3)+lam(1,1)*R(2,3)*R(3,1)
tt36 = mu(1,1)*R(3,1)*R(3,3)+lam(1,1)*R(3,1)*R(3,3)
tt37 = R(1,1)*R(2,1)
tt38 = mu(1,1)*(tt5+2*R(1,2)*R(2,2)+tt37)+lam(1,1)*R(1,2)*R(2,2)
tt39 = R(1,1)*R(3,1)
tt40 = mu(1,1)*(tt8+2*R(1,2)*R(3,2)+tt39)+lam(1,1)*R(1,2)*R(3,2)
tt41 = mu(1,1)*R(1,2)*R(1,3)+lam(1,1)*R(1,2)*R(1,3)
tt42 = lam(1,1)*R(1,2)*R(2,3)+mu(1,1)*R(1,3)*R(2,2)
tt43 = lam(1,1)*R(1,2)*R(3,3)+mu(1,1)*R(1,3)*R(3,2)
tt44 = R(2,1)*R(3,1)
tt45 = mu(1,1)*(tt20+2*R(2,2)*R(3,2)+tt44)+lam(1,1)*R(2,2)*R(3,2)&
&
tt46 = mu(1,1)*R(1,2)*R(2,3)+lam(1,1)*R(1,3)*R(2,2)
tt47 = mu(1,1)*R(2,2)*R(2,3)+lam(1,1)*R(2,2)*R(2,3)
tt48 = lam(1,1)*R(2,2)*R(3,3)+mu(1,1)*R(2,3)*R(3,2)
tt49 = mu(1,1)*R(1,2)*R(3,3)+lam(1,1)*R(1,3)*R(3,2)
tt50 = mu(1,1)*R(2,2)*R(3,3)+lam(1,1)*R(2,3)*R(3,2)
tt51 = mu(1,1)*R(3,2)*R(3,3)+lam(1,1)*R(3,2)*R(3,3)
tt52 = mu(1,1)*(2*R(1,3)*R(2,3)+tt4+tt37)+lam(1,1)*R(1,3)*R(2,3)
tt53 = mu(1,1)*(2*R(1,3)*R(3,3)+tt7+tt39)+lam(1,1)*R(1,3)*R(3,3)
tt54 = mu(1,1)*(2*R(2,3)*R(3,3)+tt19+tt44)+lam(1,1)*R(2,3)*R(3,3)&
&
hes(1,1) = mu(1,1)*(tt3+tt2+2*tt1)+tt1*lam(1,1)
hes(1,2) = tt6
hes(1,3) = tt9
hes(1,4) = tt10
hes(1,5) = tt11
hes(1,6) = tt12
hes(1,7) = tt13
hes(1,8) = tt14
hes(1,9) = tt15
hes(2,1) = tt6
hes(2,2) = mu(1,1)*(tt18+tt17+2*tt16)+lam(1,1)*tt16
hes(2,3) = tt21
hes(2,4) = tt22
hes(2,5) = tt23
hes(2,6) = tt24
hes(2,7) = tt25
hes(2,8) = tt26
hes(2,9) = tt27
hes(3,1) = tt9
hes(3,2) = tt21
hes(3,3) = mu(1,1)*(tt30+tt29+2*tt28)+lam(1,1)*tt28
hes(3,4) = tt31
hes(3,5) = tt32
hes(3,6) = tt33
hes(3,7) = tt34
hes(3,8) = tt35
hes(3,9) = tt36
hes(4,1) = tt10
hes(4,2) = tt22
hes(4,3) = tt31
hes(4,4) = mu(1,1)*(tt3+2*tt2+tt1)+lam(1,1)*tt2
hes(4,5) = tt38
hes(4,6) = tt40
hes(4,7) = tt41
hes(4,8) = tt42
hes(4,9) = tt43
hes(5,1) = tt11
hes(5,2) = tt23
hes(5,3) = tt32
hes(5,4) = tt38
hes(5,5) = mu(1,1)*(tt18+2*tt17+tt16)+lam(1,1)*tt17
hes(5,6) = tt45
hes(5,7) = tt46
hes(5,8) = tt47
hes(5,9) = tt48
hes(6,1) = tt12
hes(6,2) = tt24
hes(6,3) = tt33
hes(6,4) = tt40
hes(6,5) = tt45
hes(6,6) = mu(1,1)*(tt30+2*tt29+tt28)+lam(1,1)*tt29
hes(6,7) = tt49
hes(6,8) = tt50
hes(6,9) = tt51
hes(7,1) = tt13
hes(7,2) = tt25
hes(7,3) = tt34
hes(7,4) = tt41
hes(7,5) = tt46
hes(7,6) = tt49
hes(7,7) = mu(1,1)*(2*tt3+tt2+tt1)+lam(1,1)*tt3
hes(7,8) = tt52
hes(7,9) = tt53
hes(8,1) = tt14
hes(8,2) = tt26
hes(8,3) = tt35
hes(8,4) = tt42
hes(8,5) = tt47
hes(8,6) = tt50
hes(8,7) = tt52
hes(8,8) = mu(1,1)*(2*tt18+tt17+tt16)+lam(1,1)*tt18
hes(8,9) = tt54
hes(9,1) = tt15
hes(9,2) = tt27
hes(9,3) = tt36
hes(9,4) = tt43
hes(9,5) = tt48
hes(9,6) = tt51
hes(9,7) = tt53
hes(9,8) = tt54
hes(9,9) = mu(1,1)*(2*tt30+tt29+tt28)+lam(1,1)*tt30
END 
SUBROUTINE hex_neo_F_at_qr(val, F, mu, lam) 
IMPLICIT NONE 
REAL(KIND=8) val(1, 1) 
REAL(KIND=8) F(3, 3) 
REAL(KIND=8) mu(1, 1) 
REAL(KIND=8) lam(1, 1) 
REAL(KIND=8)  tt1 
tt1 = log(F(1,1)*(F(2,2)*F(3,3)-F(2,3)*F(3,2))-F(1,2)*(F(2,1)*F(3&
&,3)-F(2,3)*F(3,1))+F(1,3)*(F(2,1)*F(3,2)-F(2,2)*F(3,1)))
val(1,1) = (lam(1,1)*tt1**2)/2.0E+0+(mu(1,1)*((-2*tt1)+F(3,3)**2+&
&F(3,2)**2+F(3,1)**2+F(2,3)**2+F(2,2)**2+F(2,1)**2+F(1,3)**2+F(1,2&
&)**2+F(1,1)**2-3))/2.0E+0
END 
SUBROUTINE hex_neo_F_at_qr_jac(jac, F, mu, lam) 
IMPLICIT NONE 
REAL(KIND=8) jac(1, 9) 
REAL(KIND=8) F(3, 3) 
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
tt1 = F(2,2)*F(3,3)-F(2,3)*F(3,2)
tt2 = F(2,1)*F(3,2)-F(2,2)*F(3,1)
tt3 = F(1,1)*tt1-F(1,2)*(F(2,1)*F(3,3)-F(2,3)*F(3,1))+F(1,3)*tt2
tt4 = 1/tt3
tt5 = log(tt3)
tt6 = F(1,3)*F(3,2)-F(1,2)*F(3,3)
tt7 = F(1,2)*F(2,3)-F(1,3)*F(2,2)
tt8 = F(2,3)*F(3,1)-F(2,1)*F(3,3)
tt9 = F(1,1)*F(3,3)-F(1,3)*F(3,1)
tt10 = F(1,3)*F(2,1)-F(1,1)*F(2,3)
tt11 = F(1,2)*F(3,1)-F(1,1)*F(3,2)
tt12 = F(1,1)*F(2,2)-F(1,2)*F(2,1)
jac(1,1) = lam(1,1)*tt1*tt4*tt5+(mu(1,1)*(2*F(1,1)-2*tt1*tt4))/2.&
&0E+0
jac(1,2) = lam(1,1)*tt6*tt4*tt5+(mu(1,1)*(2*F(2,1)-2*tt6*tt4))/2.&
&0E+0
jac(1,3) = lam(1,1)*tt7*tt4*tt5+(mu(1,1)*(2*F(3,1)-2*tt7*tt4))/2.&
&0E+0
jac(1,4) = lam(1,1)*tt8*tt4*tt5+(mu(1,1)*(2*F(1,2)-2*tt8*tt4))/2.&
&0E+0
jac(1,5) = lam(1,1)*tt9*tt4*tt5+(mu(1,1)*(2*F(2,2)-2*tt9*tt4))/2.&
&0E+0
jac(1,6) = lam(1,1)*tt10*tt4*tt5+(mu(1,1)*(2*F(3,2)-2*tt10*tt4))/&
&2.0E+0
jac(1,7) = lam(1,1)*tt2*tt4*tt5+(mu(1,1)*(2*F(1,3)-2*tt2*tt4))/2.&
&0E+0
jac(1,8) = lam(1,1)*tt11*tt4*tt5+(mu(1,1)*(2*F(2,3)-2*tt11*tt4))/&
&2.0E+0
jac(1,9) = lam(1,1)*tt12*tt4*tt5+(mu(1,1)*(2*F(3,3)-2*tt12*tt4))/&
&2.0E+0
END 
SUBROUTINE hex_neo_F_at_qr_hes(hes, F, mu, lam) 
IMPLICIT NONE 
REAL(KIND=8) hes(9, 9) 
REAL(KIND=8) F(3, 3) 
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
tt1 = F(2,2)*F(3,3)-F(2,3)*F(3,2)
tt2 = tt1**2
tt3 = F(2,1)*F(3,2)-F(2,2)*F(3,1)
tt4 = F(1,1)*tt1-F(1,2)*(F(2,1)*F(3,3)-F(2,3)*F(3,1))+F(1,3)*tt3
tt5 = 1/tt4**2
tt6 = log(tt4)
tt7 = F(1,3)*F(3,2)-F(1,2)*F(3,3)
tt8 = (-lam(1,1)*tt7*tt1*tt5*tt6)+mu(1,1)*tt7*tt1*tt5+lam(1,1)*tt&
&7*tt1*tt5
tt9 = F(1,2)*F(2,3)-F(1,3)*F(2,2)
tt10 = (-lam(1,1)*tt9*tt1*tt5*tt6)+mu(1,1)*tt9*tt1*tt5+lam(1,1)*t&
&t9*tt1*tt5
tt11 = F(2,3)*F(3,1)-F(2,1)*F(3,3)
tt12 = (-lam(1,1)*tt11*tt1*tt5*tt6)+mu(1,1)*tt11*tt1*tt5+lam(1,1)&
&*tt11*tt1*tt5
tt13 = F(1,1)*F(3,3)-F(1,3)*F(3,1)
tt14 = 1/tt4
tt15 = lam(1,1)*F(3,3)*tt14*tt6-lam(1,1)*tt13*tt1*tt5*tt6+(mu(1,1&
&)*(2*tt13*tt1*tt5-2*F(3,3)*tt14))/2.0E+0+lam(1,1)*tt13*tt1*tt5
tt16 = F(1,3)*F(2,1)-F(1,1)*F(2,3)
tt17 = (-lam(1,1)*F(2,3)*tt14*tt6)-lam(1,1)*tt16*tt1*tt5*tt6+(mu(&
&1,1)*(2*F(2,3)*tt14+2*tt16*tt1*tt5))/2.0E+0+lam(1,1)*tt16*tt1*tt5
tt18 = (-lam(1,1)*tt3*tt1*tt5*tt6)+mu(1,1)*tt3*tt1*tt5+lam(1,1)*t&
&t3*tt1*tt5
tt19 = F(1,2)*F(3,1)-F(1,1)*F(3,2)
tt20 = (-lam(1,1)*F(3,2)*tt14*tt6)-lam(1,1)*tt19*tt1*tt5*tt6+(mu(&
&1,1)*(2*F(3,2)*tt14+2*tt19*tt1*tt5))/2.0E+0+lam(1,1)*tt19*tt1*tt5
tt21 = F(1,1)*F(2,2)-F(1,2)*F(2,1)
tt22 = lam(1,1)*F(2,2)*tt14*tt6-lam(1,1)*tt21*tt1*tt5*tt6+(mu(1,1&
&)*(2*tt21*tt1*tt5-2*F(2,2)*tt14))/2.0E+0+lam(1,1)*tt21*tt1*tt5
tt23 = tt7**2
tt24 = (-lam(1,1)*tt9*tt7*tt5*tt6)+mu(1,1)*tt9*tt7*tt5+lam(1,1)*t&
&t9*tt7*tt5
tt25 = (-lam(1,1)*F(3,3)*tt14*tt6)-lam(1,1)*tt7*tt11*tt5*tt6+(mu(&
&1,1)*(2*F(3,3)*tt14+2*tt7*tt11*tt5))/2.0E+0+lam(1,1)*tt7*tt11*tt5
tt26 = (-lam(1,1)*tt13*tt7*tt5*tt6)+mu(1,1)*tt13*tt7*tt5+lam(1,1)&
&*tt13*tt7*tt5
tt27 = lam(1,1)*F(1,3)*tt14*tt6-lam(1,1)*tt16*tt7*tt5*tt6+(mu(1,1&
&)*(2*tt16*tt7*tt5-2*F(1,3)*tt14))/2.0E+0+lam(1,1)*tt16*tt7*tt5
tt28 = lam(1,1)*F(3,2)*tt14*tt6-lam(1,1)*tt3*tt7*tt5*tt6+(mu(1,1)&
&*(2*tt3*tt7*tt5-2*F(3,2)*tt14))/2.0E+0+lam(1,1)*tt3*tt7*tt5
tt29 = (-lam(1,1)*tt19*tt7*tt5*tt6)+mu(1,1)*tt19*tt7*tt5+lam(1,1)&
&*tt19*tt7*tt5
tt30 = (-lam(1,1)*F(1,2)*tt14*tt6)-lam(1,1)*tt21*tt7*tt5*tt6+(mu(&
&1,1)*(2*F(1,2)*tt14+2*tt21*tt7*tt5))/2.0E+0+lam(1,1)*tt21*tt7*tt5
tt31 = tt9**2
tt32 = lam(1,1)*F(2,3)*tt14*tt6-lam(1,1)*tt9*tt11*tt5*tt6+(mu(1,1&
&)*(2*tt9*tt11*tt5-2*F(2,3)*tt14))/2.0E+0+lam(1,1)*tt9*tt11*tt5
tt33 = (-lam(1,1)*F(1,3)*tt14*tt6)-lam(1,1)*tt9*tt13*tt5*tt6+(mu(&
&1,1)*(2*F(1,3)*tt14+2*tt9*tt13*tt5))/2.0E+0+lam(1,1)*tt9*tt13*tt5
tt34 = (-lam(1,1)*tt16*tt9*tt5*tt6)+mu(1,1)*tt16*tt9*tt5+lam(1,1)&
&*tt16*tt9*tt5
tt35 = (-lam(1,1)*F(2,2)*tt14*tt6)-lam(1,1)*tt9*tt3*tt5*tt6+(mu(1&
&,1)*(2*F(2,2)*tt14+2*tt9*tt3*tt5))/2.0E+0+lam(1,1)*tt9*tt3*tt5
tt36 = lam(1,1)*F(1,2)*tt14*tt6-lam(1,1)*tt9*tt19*tt5*tt6+(mu(1,1&
&)*(2*tt9*tt19*tt5-2*F(1,2)*tt14))/2.0E+0+lam(1,1)*tt9*tt19*tt5
tt37 = (-lam(1,1)*tt21*tt9*tt5*tt6)+mu(1,1)*tt21*tt9*tt5+lam(1,1)&
&*tt21*tt9*tt5
tt38 = tt11**2
tt39 = (-lam(1,1)*tt13*tt11*tt5*tt6)+mu(1,1)*tt13*tt11*tt5+lam(1,&
&1)*tt13*tt11*tt5
tt40 = (-lam(1,1)*tt16*tt11*tt5*tt6)+mu(1,1)*tt16*tt11*tt5+lam(1,&
&1)*tt16*tt11*tt5
tt41 = (-lam(1,1)*tt3*tt11*tt5*tt6)+mu(1,1)*tt3*tt11*tt5+lam(1,1)&
&*tt3*tt11*tt5
tt42 = lam(1,1)*F(3,1)*tt14*tt6-lam(1,1)*tt19*tt11*tt5*tt6+(mu(1,&
&1)*(2*tt19*tt11*tt5-2*F(3,1)*tt14))/2.0E+0+lam(1,1)*tt19*tt11*tt5
tt43 = (-lam(1,1)*F(2,1)*tt14*tt6)-lam(1,1)*tt21*tt11*tt5*tt6+(mu&
&(1,1)*(2*F(2,1)*tt14+2*tt21*tt11*tt5))/2.0E+0+lam(1,1)*tt21*tt11*&
&tt5
tt44 = tt13**2
tt45 = (-lam(1,1)*tt16*tt13*tt5*tt6)+mu(1,1)*tt16*tt13*tt5+lam(1,&
&1)*tt16*tt13*tt5
tt46 = (-lam(1,1)*F(3,1)*tt14*tt6)-lam(1,1)*tt3*tt13*tt5*tt6+(mu(&
&1,1)*(2*F(3,1)*tt14+2*tt3*tt13*tt5))/2.0E+0+lam(1,1)*tt3*tt13*tt5
tt47 = (-lam(1,1)*tt19*tt13*tt5*tt6)+mu(1,1)*tt19*tt13*tt5+lam(1,&
&1)*tt19*tt13*tt5
tt48 = F(1,1)*lam(1,1)*tt14*tt6-lam(1,1)*tt21*tt13*tt5*tt6+(mu(1,&
&1)*(2*tt21*tt13*tt5-2*F(1,1)*tt14))/2.0E+0+lam(1,1)*tt21*tt13*tt5
tt49 = tt16**2
tt50 = lam(1,1)*F(2,1)*tt14*tt6-lam(1,1)*tt16*tt3*tt5*tt6+(mu(1,1&
&)*(2*tt16*tt3*tt5-2*F(2,1)*tt14))/2.0E+0+lam(1,1)*tt16*tt3*tt5
tt51 = (-F(1,1)*lam(1,1)*tt14*tt6)-lam(1,1)*tt16*tt19*tt5*tt6+(mu&
&(1,1)*(2*F(1,1)*tt14+2*tt16*tt19*tt5))/2.0E+0+lam(1,1)*tt16*tt19*&
&tt5
tt52 = (-lam(1,1)*tt21*tt16*tt5*tt6)+mu(1,1)*tt21*tt16*tt5+lam(1,&
&1)*tt21*tt16*tt5
tt53 = tt3**2
tt54 = (-lam(1,1)*tt19*tt3*tt5*tt6)+mu(1,1)*tt19*tt3*tt5+lam(1,1)&
&*tt19*tt3*tt5
tt55 = (-lam(1,1)*tt21*tt3*tt5*tt6)+mu(1,1)*tt21*tt3*tt5+lam(1,1)&
&*tt21*tt3*tt5
tt56 = tt19**2
tt57 = (-lam(1,1)*tt21*tt19*tt5*tt6)+mu(1,1)*tt21*tt19*tt5+lam(1,&
&1)*tt21*tt19*tt5
tt58 = tt21**2
hes(1,1) = (-lam(1,1)*tt2*tt5*tt6)+(mu(1,1)*(2*tt2*tt5+2))/2.0E+0&
&+lam(1,1)*tt2*tt5
hes(1,2) = tt8
hes(1,3) = tt10
hes(1,4) = tt12
hes(1,5) = tt15
hes(1,6) = tt17
hes(1,7) = tt18
hes(1,8) = tt20
hes(1,9) = tt22
hes(2,1) = tt8
hes(2,2) = (-lam(1,1)*tt23*tt5*tt6)+(mu(1,1)*(2*tt23*tt5+2))/2.0E&
&+0+lam(1,1)*tt23*tt5
hes(2,3) = tt24
hes(2,4) = tt25
hes(2,5) = tt26
hes(2,6) = tt27
hes(2,7) = tt28
hes(2,8) = tt29
hes(2,9) = tt30
hes(3,1) = tt10
hes(3,2) = tt24
hes(3,3) = (-lam(1,1)*tt31*tt5*tt6)+(mu(1,1)*(2*tt31*tt5+2))/2.0E&
&+0+lam(1,1)*tt31*tt5
hes(3,4) = tt32
hes(3,5) = tt33
hes(3,6) = tt34
hes(3,7) = tt35
hes(3,8) = tt36
hes(3,9) = tt37
hes(4,1) = tt12
hes(4,2) = tt25
hes(4,3) = tt32
hes(4,4) = (-lam(1,1)*tt38*tt5*tt6)+(mu(1,1)*(2*tt38*tt5+2))/2.0E&
&+0+lam(1,1)*tt38*tt5
hes(4,5) = tt39
hes(4,6) = tt40
hes(4,7) = tt41
hes(4,8) = tt42
hes(4,9) = tt43
hes(5,1) = tt15
hes(5,2) = tt26
hes(5,3) = tt33
hes(5,4) = tt39
hes(5,5) = (-lam(1,1)*tt44*tt5*tt6)+(mu(1,1)*(2*tt44*tt5+2))/2.0E&
&+0+lam(1,1)*tt44*tt5
hes(5,6) = tt45
hes(5,7) = tt46
hes(5,8) = tt47
hes(5,9) = tt48
hes(6,1) = tt17
hes(6,2) = tt27
hes(6,3) = tt34
hes(6,4) = tt40
hes(6,5) = tt45
hes(6,6) = (-lam(1,1)*tt49*tt5*tt6)+(mu(1,1)*(2*tt49*tt5+2))/2.0E&
&+0+lam(1,1)*tt49*tt5
hes(6,7) = tt50
hes(6,8) = tt51
hes(6,9) = tt52
hes(7,1) = tt18
hes(7,2) = tt28
hes(7,3) = tt35
hes(7,4) = tt41
hes(7,5) = tt46
hes(7,6) = tt50
hes(7,7) = (-lam(1,1)*tt53*tt5*tt6)+(mu(1,1)*(2*tt53*tt5+2))/2.0E&
&+0+lam(1,1)*tt53*tt5
hes(7,8) = tt54
hes(7,9) = tt55
hes(8,1) = tt20
hes(8,2) = tt29
hes(8,3) = tt36
hes(8,4) = tt42
hes(8,5) = tt47
hes(8,6) = tt51
hes(8,7) = tt54
hes(8,8) = (-lam(1,1)*tt56*tt5*tt6)+(mu(1,1)*(2*tt56*tt5+2))/2.0E&
&+0+lam(1,1)*tt56*tt5
hes(8,9) = tt57
hes(9,1) = tt22
hes(9,2) = tt30
hes(9,3) = tt37
hes(9,4) = tt43
hes(9,5) = tt48
hes(9,6) = tt52
hes(9,7) = tt55
hes(9,8) = tt57
hes(9,9) = (-lam(1,1)*tt58*tt5*tt6)+(mu(1,1)*(2*tt58*tt5+2))/2.0E&
&+0+lam(1,1)*tt58*tt5
END 
SUBROUTINE hex_stvk_F_at_qr(val, F, mu, lam) 
IMPLICIT NONE 
REAL(KIND=8) val(1, 1) 
REAL(KIND=8) F(3, 3) 
REAL(KIND=8) mu(1, 1) 
REAL(KIND=8) lam(1, 1) 
REAL(KIND=8)  tt1 
REAL(KIND=8)  tt2 
REAL(KIND=8)  tt3 
tt1 = F(3,1)**2+F(2,1)**2+F(1,1)**2-1
tt2 = F(3,2)**2+F(2,2)**2+F(1,2)**2-1
tt3 = F(3,3)**2+F(2,3)**2+F(1,3)**2-1
val(1,1) = (lam(1,1)*(tt3/2.0E+0+tt2/2.0E+0+tt1/2.0E+0)**2)/2.0E+&
&0+mu(1,1)*(tt3**2/4.0E+0+(F(3,2)*F(3,3)+F(2,2)*F(2,3)+F(1,2)*F(1,&
&3))**2/2.0E+0+(F(3,1)*F(3,3)+F(2,1)*F(2,3)+F(1,1)*F(1,3))**2/2.0E&
&+0+tt2**2/4.0E+0+(F(3,1)*F(3,2)+F(2,1)*F(2,2)+F(1,1)*F(1,2))**2/2&
&.0E+0+tt1**2/4.0E+0)
END 
SUBROUTINE hex_stvk_F_at_qr_jac(jac, F, mu, lam) 
IMPLICIT NONE 
REAL(KIND=8) jac(1, 9) 
REAL(KIND=8) F(3, 3) 
REAL(KIND=8) mu(1, 1) 
REAL(KIND=8) lam(1, 1) 
REAL(KIND=8)  tt1 
REAL(KIND=8)  tt2 
REAL(KIND=8)  tt3 
REAL(KIND=8)  tt4 
REAL(KIND=8)  tt5 
REAL(KIND=8)  tt6 
REAL(KIND=8)  tt7 
tt1 = F(3,1)**2+F(2,1)**2+F(1,1)**2-1
tt2 = F(3,1)*F(3,2)+F(2,1)*F(2,2)+F(1,1)*F(1,2)
tt3 = F(3,1)*F(3,3)+F(2,1)*F(2,3)+F(1,1)*F(1,3)
tt4 = F(3,2)**2+F(2,2)**2+F(1,2)**2-1
tt5 = F(3,3)**2+F(2,3)**2+F(1,3)**2-1
tt6 = tt5/2.0E+0+tt4/2.0E+0+tt1/2.0E+0
tt7 = F(3,2)*F(3,3)+F(2,2)*F(2,3)+F(1,2)*F(1,3)
jac(1,1) = F(1,1)*lam(1,1)*tt6+mu(1,1)*(F(1,3)*tt3+F(1,2)*tt2+F(1&
&,1)*tt1)
jac(1,2) = lam(1,1)*F(2,1)*tt6+mu(1,1)*(F(2,3)*tt3+F(2,2)*tt2+F(2&
&,1)*tt1)
jac(1,3) = lam(1,1)*F(3,1)*tt6+mu(1,1)*(F(3,3)*tt3+F(3,2)*tt2+F(3&
&,1)*tt1)
jac(1,4) = lam(1,1)*F(1,2)*tt6+mu(1,1)*(F(1,3)*tt7+F(1,2)*tt4+F(1&
&,1)*tt2)
jac(1,5) = lam(1,1)*F(2,2)*tt6+mu(1,1)*(F(2,3)*tt7+F(2,2)*tt4+F(2&
&,1)*tt2)
jac(1,6) = lam(1,1)*F(3,2)*tt6+mu(1,1)*(F(3,3)*tt7+F(3,2)*tt4+F(3&
&,1)*tt2)
jac(1,7) = mu(1,1)*(F(1,3)*tt5+F(1,2)*tt7+F(1,1)*tt3)+lam(1,1)*F(&
&1,3)*tt6
jac(1,8) = mu(1,1)*(F(2,3)*tt5+F(2,2)*tt7+F(2,1)*tt3)+lam(1,1)*F(&
&2,3)*tt6
jac(1,9) = mu(1,1)*(F(3,3)*tt5+F(3,2)*tt7+F(3,1)*tt3)+lam(1,1)*F(&
&3,3)*tt6
END 
SUBROUTINE hex_stvk_F_at_qr_hes(hes, F, mu, lam) 
IMPLICIT NONE 
REAL(KIND=8) hes(9, 9) 
REAL(KIND=8) F(3, 3) 
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
tt1 = F(1,1)**2
tt2 = F(1,2)**2
tt3 = F(1,3)**2
tt4 = F(2,1)**2
tt5 = F(3,1)**2
tt6 = F(2,2)**2
tt7 = F(3,2)**2
tt8 = F(2,3)**2
tt9 = F(3,3)**2
tt10 = lam(1,1)*((tt9+tt8+tt3-1)/2.0E+0+(tt7+tt6+tt2-1)/2.0E+0+(t&
&t5+tt4+tt1-1)/2.0E+0)
tt11 = F(1,2)*F(2,2)
tt12 = F(1,3)*F(2,3)
tt13 = mu(1,1)*(tt12+tt11+2*F(1,1)*F(2,1))+F(1,1)*lam(1,1)*F(2,1)&
&
tt14 = F(1,2)*F(3,2)
tt15 = F(1,3)*F(3,3)
tt16 = mu(1,1)*(tt15+tt14+2*F(1,1)*F(3,1))+F(1,1)*lam(1,1)*F(3,1)&
&
tt17 = F(2,1)*F(2,2)
tt18 = F(3,1)*F(3,2)
tt19 = mu(1,1)*(tt18+tt17+2*F(1,1)*F(1,2))+F(1,1)*lam(1,1)*F(1,2)&
&
tt20 = F(1,1)*lam(1,1)*F(2,2)+mu(1,1)*F(1,2)*F(2,1)
tt21 = F(1,1)*lam(1,1)*F(3,2)+mu(1,1)*F(1,2)*F(3,1)
tt22 = F(2,1)*F(2,3)
tt23 = F(3,1)*F(3,3)
tt24 = mu(1,1)*(tt23+tt22+2*F(1,1)*F(1,3))+F(1,1)*lam(1,1)*F(1,3)&
&
tt25 = F(1,1)*lam(1,1)*F(2,3)+mu(1,1)*F(1,3)*F(2,1)
tt26 = F(1,1)*lam(1,1)*F(3,3)+mu(1,1)*F(1,3)*F(3,1)
tt27 = F(2,2)*F(3,2)
tt28 = F(2,3)*F(3,3)
tt29 = mu(1,1)*(tt28+tt27+2*F(2,1)*F(3,1))+lam(1,1)*F(2,1)*F(3,1)&
&
tt30 = F(1,1)*mu(1,1)*F(2,2)+lam(1,1)*F(1,2)*F(2,1)
tt31 = F(1,1)*F(1,2)
tt32 = mu(1,1)*(tt18+2*F(2,1)*F(2,2)+tt31)+lam(1,1)*F(2,1)*F(2,2)&
&
tt33 = lam(1,1)*F(2,1)*F(3,2)+mu(1,1)*F(2,2)*F(3,1)
tt34 = F(1,1)*mu(1,1)*F(2,3)+lam(1,1)*F(1,3)*F(2,1)
tt35 = F(1,1)*F(1,3)
tt36 = mu(1,1)*(tt23+2*F(2,1)*F(2,3)+tt35)+lam(1,1)*F(2,1)*F(2,3)&
&
tt37 = lam(1,1)*F(2,1)*F(3,3)+mu(1,1)*F(2,3)*F(3,1)
tt38 = F(1,1)*mu(1,1)*F(3,2)+lam(1,1)*F(1,2)*F(3,1)
tt39 = mu(1,1)*F(2,1)*F(3,2)+lam(1,1)*F(2,2)*F(3,1)
tt40 = mu(1,1)*(2*F(3,1)*F(3,2)+tt17+tt31)+lam(1,1)*F(3,1)*F(3,2)&
&
tt41 = F(1,1)*mu(1,1)*F(3,3)+lam(1,1)*F(1,3)*F(3,1)
tt42 = mu(1,1)*F(2,1)*F(3,3)+lam(1,1)*F(2,3)*F(3,1)
tt43 = mu(1,1)*(2*F(3,1)*F(3,3)+tt22+tt35)+lam(1,1)*F(3,1)*F(3,3)&
&
tt44 = F(1,1)*F(2,1)
tt45 = mu(1,1)*(tt12+2*F(1,2)*F(2,2)+tt44)+lam(1,1)*F(1,2)*F(2,2)&
&
tt46 = F(1,1)*F(3,1)
tt47 = mu(1,1)*(tt15+2*F(1,2)*F(3,2)+tt46)+lam(1,1)*F(1,2)*F(3,2)&
&
tt48 = F(2,2)*F(2,3)
tt49 = F(3,2)*F(3,3)
tt50 = mu(1,1)*(tt49+tt48+2*F(1,2)*F(1,3))+lam(1,1)*F(1,2)*F(1,3)&
&
tt51 = lam(1,1)*F(1,2)*F(2,3)+mu(1,1)*F(1,3)*F(2,2)
tt52 = lam(1,1)*F(1,2)*F(3,3)+mu(1,1)*F(1,3)*F(3,2)
tt53 = F(2,1)*F(3,1)
tt54 = mu(1,1)*(tt28+2*F(2,2)*F(3,2)+tt53)+lam(1,1)*F(2,2)*F(3,2)&
&
tt55 = mu(1,1)*F(1,2)*F(2,3)+lam(1,1)*F(1,3)*F(2,2)
tt56 = F(1,2)*F(1,3)
tt57 = mu(1,1)*(tt49+2*F(2,2)*F(2,3)+tt56)+lam(1,1)*F(2,2)*F(2,3)&
&
tt58 = lam(1,1)*F(2,2)*F(3,3)+mu(1,1)*F(2,3)*F(3,2)
tt59 = mu(1,1)*F(1,2)*F(3,3)+lam(1,1)*F(1,3)*F(3,2)
tt60 = mu(1,1)*F(2,2)*F(3,3)+lam(1,1)*F(2,3)*F(3,2)
tt61 = mu(1,1)*(2*F(3,2)*F(3,3)+tt48+tt56)+lam(1,1)*F(3,2)*F(3,3)&
&
tt62 = mu(1,1)*(2*F(1,3)*F(2,3)+tt11+tt44)+lam(1,1)*F(1,3)*F(2,3)&
&
tt63 = mu(1,1)*(2*F(1,3)*F(3,3)+tt14+tt46)+lam(1,1)*F(1,3)*F(3,3)&
&
tt64 = mu(1,1)*(2*F(2,3)*F(3,3)+tt27+tt53)+lam(1,1)*F(2,3)*F(3,3)&
&
hes(1,1) = tt10+mu(1,1)*(tt5+tt4+tt3+tt2+3*tt1-1)+tt1*lam(1,1)
hes(1,2) = tt13
hes(1,3) = tt16
hes(1,4) = tt19
hes(1,5) = tt20
hes(1,6) = tt21
hes(1,7) = tt24
hes(1,8) = tt25
hes(1,9) = tt26
hes(2,1) = tt13
hes(2,2) = tt10+mu(1,1)*(tt5+tt8+tt6+3*tt4+tt1-1)+lam(1,1)*tt4
hes(2,3) = tt29
hes(2,4) = tt30
hes(2,5) = tt32
hes(2,6) = tt33
hes(2,7) = tt34
hes(2,8) = tt36
hes(2,9) = tt37
hes(3,1) = tt16
hes(3,2) = tt29
hes(3,3) = tt10+mu(1,1)*(tt9+tt7+3*tt5+tt4+tt1-1)+lam(1,1)*tt5
hes(3,4) = tt38
hes(3,5) = tt39
hes(3,6) = tt40
hes(3,7) = tt41
hes(3,8) = tt42
hes(3,9) = tt43
hes(4,1) = tt19
hes(4,2) = tt30
hes(4,3) = tt38
hes(4,4) = tt10+mu(1,1)*(tt7+tt6+tt3+3*tt2+tt1-1)+lam(1,1)*tt2
hes(4,5) = tt45
hes(4,6) = tt47
hes(4,7) = tt50
hes(4,8) = tt51
hes(4,9) = tt52
hes(5,1) = tt20
hes(5,2) = tt32
hes(5,3) = tt39
hes(5,4) = tt45
hes(5,5) = tt10+mu(1,1)*(tt7+tt8+3*tt6+tt4+tt2-1)+lam(1,1)*tt6
hes(5,6) = tt54
hes(5,7) = tt55
hes(5,8) = tt57
hes(5,9) = tt58
hes(6,1) = tt21
hes(6,2) = tt33
hes(6,3) = tt40
hes(6,4) = tt47
hes(6,5) = tt54
hes(6,6) = tt10+mu(1,1)*(tt9+3*tt7+tt5+tt6+tt2-1)+lam(1,1)*tt7
hes(6,7) = tt59
hes(6,8) = tt60
hes(6,9) = tt61
hes(7,1) = tt24
hes(7,2) = tt34
hes(7,3) = tt41
hes(7,4) = tt50
hes(7,5) = tt55
hes(7,6) = tt59
hes(7,7) = tt10+mu(1,1)*(tt9+tt8+3*tt3+tt2+tt1-1)+lam(1,1)*tt3
hes(7,8) = tt62
hes(7,9) = tt63
hes(8,1) = tt25
hes(8,2) = tt36
hes(8,3) = tt42
hes(8,4) = tt51
hes(8,5) = tt57
hes(8,6) = tt60
hes(8,7) = tt62
hes(8,8) = tt10+mu(1,1)*(tt9+3*tt8+tt6+tt4+tt3-1)+lam(1,1)*tt8
hes(8,9) = tt64
hes(9,1) = tt26
hes(9,2) = tt37
hes(9,3) = tt43
hes(9,4) = tt52
hes(9,5) = tt58
hes(9,6) = tt61
hes(9,7) = tt63
hes(9,8) = tt64
hes(9,9) = tt10+mu(1,1)*(3*tt9+tt7+tt5+tt8+tt3-1)+lam(1,1)*tt9
END 
SUBROUTINE hex_sta_neo_F_at_qr(val, F, mu, lam) 
IMPLICIT NONE 
REAL(KIND=8) val(1, 1) 
REAL(KIND=8) F(3, 3) 
REAL(KIND=8) mu(1, 1) 
REAL(KIND=8) lam(1, 1) 
REAL(KIND=8)  tt1 
tt1 = F(1,1)*(F(2,2)*F(3,3)-F(2,3)*F(3,2))-F(1,2)*(F(2,1)*F(3,3)-&
&F(2,3)*F(3,1))+F(1,3)*(F(2,1)*F(3,2)-F(2,2)*F(3,1))-1
val(1,1) = ((mu(1,1)+lam(1,1))*tt1**2)/2.0E+0+(mu(1,1)*(F(3,3)**2&
&+F(3,2)**2+F(3,1)**2+F(2,3)**2+F(2,2)**2+F(2,1)**2+F(1,3)**2+F(1,&
&2)**2+F(1,1)**2-3))/2.0E+0-mu(1,1)*tt1
END 
SUBROUTINE hex_sta_neo_F_at_qr_jac(jac, F, mu, lam) 
IMPLICIT NONE 
REAL(KIND=8) jac(1, 9) 
REAL(KIND=8) F(3, 3) 
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
tt1 = F(2,2)*F(3,3)-F(2,3)*F(3,2)
tt2 = mu(1,1)+lam(1,1)
tt3 = F(2,1)*F(3,2)-F(2,2)*F(3,1)
tt4 = F(1,1)*tt1-F(1,2)*(F(2,1)*F(3,3)-F(2,3)*F(3,1))+F(1,3)*tt3-&
&1
tt5 = F(1,3)*F(3,2)-F(1,2)*F(3,3)
tt6 = F(1,2)*F(2,3)-F(1,3)*F(2,2)
tt7 = F(2,3)*F(3,1)-F(2,1)*F(3,3)
tt8 = F(1,1)*F(3,3)-F(1,3)*F(3,1)
tt9 = F(1,3)*F(2,1)-F(1,1)*F(2,3)
tt10 = F(1,2)*F(3,1)-F(1,1)*F(3,2)
tt11 = F(1,1)*F(2,2)-F(1,2)*F(2,1)
jac(1,1) = tt2*tt1*tt4-mu(1,1)*tt1+F(1,1)*mu(1,1)
jac(1,2) = tt2*tt5*tt4-mu(1,1)*tt5+mu(1,1)*F(2,1)
jac(1,3) = tt2*tt6*tt4+mu(1,1)*F(3,1)-mu(1,1)*tt6
jac(1,4) = tt2*tt7*tt4-mu(1,1)*tt7+mu(1,1)*F(1,2)
jac(1,5) = tt2*tt8*tt4-mu(1,1)*tt8+mu(1,1)*F(2,2)
jac(1,6) = tt2*tt9*tt4+mu(1,1)*F(3,2)-mu(1,1)*tt9
jac(1,7) = tt2*tt3*tt4-mu(1,1)*tt3+mu(1,1)*F(1,3)
jac(1,8) = tt2*tt10*tt4-mu(1,1)*tt10+mu(1,1)*F(2,3)
jac(1,9) = tt2*tt11*tt4+mu(1,1)*F(3,3)-mu(1,1)*tt11
END 
SUBROUTINE hex_sta_neo_F_at_qr_hes(hes, F, mu, lam) 
IMPLICIT NONE 
REAL(KIND=8) hes(9, 9) 
REAL(KIND=8) F(3, 3) 
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
tt1 = mu(1,1)+lam(1,1)
tt2 = F(2,2)*F(3,3)-F(2,3)*F(3,2)
tt3 = F(1,3)*F(3,2)-F(1,2)*F(3,3)
tt4 = tt1*tt3*tt2
tt5 = F(1,2)*F(2,3)-F(1,3)*F(2,2)
tt6 = tt1*tt5*tt2
tt7 = F(2,3)*F(3,1)-F(2,1)*F(3,3)
tt8 = tt1*tt7*tt2
tt9 = F(1,1)*F(3,3)-F(1,3)*F(3,1)
tt10 = F(2,1)*F(3,2)-F(2,2)*F(3,1)
tt11 = F(1,1)*tt2-F(1,2)*(F(2,1)*F(3,3)-F(2,3)*F(3,1))+F(1,3)*tt1&
&0-1
tt12 = tt1*F(3,3)*tt11+tt1*tt9*tt2-mu(1,1)*F(3,3)
tt13 = F(1,3)*F(2,1)-F(1,1)*F(2,3)
tt14 = (-tt1*F(2,3)*tt11)+tt1*tt13*tt2+mu(1,1)*F(2,3)
tt15 = tt1*tt10*tt2
tt16 = F(1,2)*F(3,1)-F(1,1)*F(3,2)
tt17 = (-tt1*F(3,2)*tt11)+tt1*tt16*tt2+mu(1,1)*F(3,2)
tt18 = F(1,1)*F(2,2)-F(1,2)*F(2,1)
tt19 = tt1*F(2,2)*tt11+tt1*tt18*tt2-mu(1,1)*F(2,2)
tt20 = tt1*tt5*tt3
tt21 = (-tt1*F(3,3)*tt11)+tt1*tt3*tt7+mu(1,1)*F(3,3)
tt22 = tt1*tt9*tt3
tt23 = tt1*F(1,3)*tt11+tt1*tt13*tt3-mu(1,1)*F(1,3)
tt24 = tt1*F(3,2)*tt11+tt1*tt10*tt3-mu(1,1)*F(3,2)
tt25 = tt1*tt16*tt3
tt26 = (-tt1*F(1,2)*tt11)+tt1*tt18*tt3+mu(1,1)*F(1,2)
tt27 = tt1*F(2,3)*tt11+tt1*tt5*tt7-mu(1,1)*F(2,3)
tt28 = (-tt1*F(1,3)*tt11)+tt1*tt5*tt9+mu(1,1)*F(1,3)
tt29 = tt1*tt13*tt5
tt30 = (-tt1*F(2,2)*tt11)+tt1*tt5*tt10+mu(1,1)*F(2,2)
tt31 = tt1*F(1,2)*tt11+tt1*tt5*tt16-mu(1,1)*F(1,2)
tt32 = tt1*tt18*tt5
tt33 = tt1*tt9*tt7
tt34 = tt1*tt13*tt7
tt35 = tt1*tt10*tt7
tt36 = tt1*F(3,1)*tt11+tt1*tt16*tt7-mu(1,1)*F(3,1)
tt37 = (-tt1*F(2,1)*tt11)+tt1*tt18*tt7+mu(1,1)*F(2,1)
tt38 = tt1*tt13*tt9
tt39 = (-tt1*F(3,1)*tt11)+tt1*tt10*tt9+mu(1,1)*F(3,1)
tt40 = tt1*tt16*tt9
tt41 = F(1,1)*tt1*tt11+tt1*tt18*tt9-F(1,1)*mu(1,1)
tt42 = tt1*F(2,1)*tt11+tt1*tt13*tt10-mu(1,1)*F(2,1)
tt43 = (-F(1,1)*tt1*tt11)+tt1*tt13*tt16+F(1,1)*mu(1,1)
tt44 = tt1*tt18*tt13
tt45 = tt1*tt16*tt10
tt46 = tt1*tt18*tt10
tt47 = tt1*tt18*tt16
hes(1,1) = tt1*tt2**2+mu(1,1)
hes(1,2) = tt4
hes(1,3) = tt6
hes(1,4) = tt8
hes(1,5) = tt12
hes(1,6) = tt14
hes(1,7) = tt15
hes(1,8) = tt17
hes(1,9) = tt19
hes(2,1) = tt4
hes(2,2) = tt1*tt3**2+mu(1,1)
hes(2,3) = tt20
hes(2,4) = tt21
hes(2,5) = tt22
hes(2,6) = tt23
hes(2,7) = tt24
hes(2,8) = tt25
hes(2,9) = tt26
hes(3,1) = tt6
hes(3,2) = tt20
hes(3,3) = tt1*tt5**2+mu(1,1)
hes(3,4) = tt27
hes(3,5) = tt28
hes(3,6) = tt29
hes(3,7) = tt30
hes(3,8) = tt31
hes(3,9) = tt32
hes(4,1) = tt8
hes(4,2) = tt21
hes(4,3) = tt27
hes(4,4) = tt1*tt7**2+mu(1,1)
hes(4,5) = tt33
hes(4,6) = tt34
hes(4,7) = tt35
hes(4,8) = tt36
hes(4,9) = tt37
hes(5,1) = tt12
hes(5,2) = tt22
hes(5,3) = tt28
hes(5,4) = tt33
hes(5,5) = tt1*tt9**2+mu(1,1)
hes(5,6) = tt38
hes(5,7) = tt39
hes(5,8) = tt40
hes(5,9) = tt41
hes(6,1) = tt14
hes(6,2) = tt23
hes(6,3) = tt29
hes(6,4) = tt34
hes(6,5) = tt38
hes(6,6) = tt1*tt13**2+mu(1,1)
hes(6,7) = tt42
hes(6,8) = tt43
hes(6,9) = tt44
hes(7,1) = tt15
hes(7,2) = tt24
hes(7,3) = tt30
hes(7,4) = tt35
hes(7,5) = tt39
hes(7,6) = tt42
hes(7,7) = tt1*tt10**2+mu(1,1)
hes(7,8) = tt45
hes(7,9) = tt46
hes(8,1) = tt17
hes(8,2) = tt25
hes(8,3) = tt31
hes(8,4) = tt36
hes(8,5) = tt40
hes(8,6) = tt43
hes(8,7) = tt45
hes(8,8) = tt1*tt16**2+mu(1,1)
hes(8,9) = tt47
hes(9,1) = tt19
hes(9,2) = tt26
hes(9,3) = tt32
hes(9,4) = tt37
hes(9,5) = tt41
hes(9,6) = tt44
hes(9,7) = tt46
hes(9,8) = tt47
hes(9,9) = tt1*tt18**2+mu(1,1)
END 
SUBROUTINE hex_bower_neo_F_at_qr(val, F, mu, lam) 
IMPLICIT NONE 
REAL(KIND=8) val(1, 1) 
REAL(KIND=8) F(3, 3) 
REAL(KIND=8) mu(1, 1) 
REAL(KIND=8) lam(1, 1) 
REAL(KIND=8)  tt1 
REAL(KIND=8)  tt2 
REAL(KIND=8)  tt3 
tt1 = F(1,3)*(F(2,1)*F(3,2)-F(2,2)*F(3,1))
tt2 = -F(1,2)*(F(2,1)*F(3,3)-F(2,3)*F(3,1))
tt3 = F(1,1)*(F(2,2)*F(3,3)-F(2,3)*F(3,2))
val(1,1) = (lam(1,1)*(tt3+tt2+tt1-1)**2)/2.0E+0+(mu(1,1)*((F(3,3)&
&**2+F(3,2)**2+F(3,1)**2+F(2,3)**2+F(2,2)**2+F(2,1)**2+F(1,3)**2+F&
&(1,2)**2+F(1,1)**2)/(tt3+tt2+tt1)**(2.0E+0/3.0E+0)-3))/2.0E+0
END 
SUBROUTINE hex_bower_neo_F_at_qr_jac(jac, F, mu, lam) 
IMPLICIT NONE 
REAL(KIND=8) jac(1, 9) 
REAL(KIND=8) F(3, 3) 
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
tt1 = F(2,2)*F(3,3)-F(2,3)*F(3,2)
tt2 = F(2,1)*F(3,2)-F(2,2)*F(3,1)
tt3 = F(1,3)*tt2
tt4 = -F(1,2)*(F(2,1)*F(3,3)-F(2,3)*F(3,1))
tt5 = F(1,1)*tt1
tt6 = tt5+tt4+tt3-1
tt7 = tt5+tt4+tt3
tt8 = tt7**((-5.0E+0)/3.0E+0)
tt9 = F(3,3)**2+F(3,2)**2+F(3,1)**2+F(2,3)**2+F(2,2)**2+F(2,1)**2&
&+F(1,3)**2+F(1,2)**2+F(1,1)**2
tt10 = tt7**((-2.0E+0)/3.0E+0)
tt11 = F(1,3)*F(3,2)-F(1,2)*F(3,3)
tt12 = F(1,2)*F(2,3)-F(1,3)*F(2,2)
tt13 = F(2,3)*F(3,1)-F(2,1)*F(3,3)
tt14 = F(1,1)*F(3,3)-F(1,3)*F(3,1)
tt15 = F(1,3)*F(2,1)-F(1,1)*F(2,3)
tt16 = F(1,2)*F(3,1)-F(1,1)*F(3,2)
tt17 = F(1,1)*F(2,2)-F(1,2)*F(2,1)
jac(1,1) = (mu(1,1)*(2*F(1,1)*tt10+((-2.0E+0)*tt1*tt8*tt9)/3.0E+0&
&))/2.0E+0+lam(1,1)*tt1*tt6
jac(1,2) = (mu(1,1)*(2*F(2,1)*tt10+((-2.0E+0)*tt11*tt8*tt9)/3.0E+&
&0))/2.0E+0+lam(1,1)*tt11*tt6
jac(1,3) = (mu(1,1)*(2*F(3,1)*tt10+((-2.0E+0)*tt12*tt8*tt9)/3.0E+&
&0))/2.0E+0+lam(1,1)*tt12*tt6
jac(1,4) = (mu(1,1)*(2*F(1,2)*tt10+((-2.0E+0)*tt13*tt8*tt9)/3.0E+&
&0))/2.0E+0+lam(1,1)*tt13*tt6
jac(1,5) = (mu(1,1)*(2*F(2,2)*tt10+((-2.0E+0)*tt14*tt8*tt9)/3.0E+&
&0))/2.0E+0+lam(1,1)*tt14*tt6
jac(1,6) = (mu(1,1)*(2*F(3,2)*tt10+((-2.0E+0)*tt15*tt8*tt9)/3.0E+&
&0))/2.0E+0+lam(1,1)*tt15*tt6
jac(1,7) = (mu(1,1)*(2*F(1,3)*tt10+((-2.0E+0)*tt2*tt8*tt9)/3.0E+0&
&))/2.0E+0+lam(1,1)*tt2*tt6
jac(1,8) = (mu(1,1)*(2*F(2,3)*tt10+((-2.0E+0)*tt16*tt8*tt9)/3.0E+&
&0))/2.0E+0+lam(1,1)*tt16*tt6
jac(1,9) = (mu(1,1)*(2*F(3,3)*tt10+((-2.0E+0)*tt17*tt8*tt9)/3.0E+&
&0))/2.0E+0+lam(1,1)*tt17*tt6
END 
SUBROUTINE hex_bower_neo_F_at_qr_hes(hes, F, mu, lam) 
IMPLICIT NONE 
REAL(KIND=8) hes(9, 9) 
REAL(KIND=8) F(3, 3) 
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
tt1 = F(2,2)*F(3,3)-F(2,3)*F(3,2)
tt2 = tt1**2
tt3 = F(2,1)*F(3,2)-F(2,2)*F(3,1)
tt4 = F(1,3)*tt3
tt5 = -F(1,2)*(F(2,1)*F(3,3)-F(2,3)*F(3,1))
tt6 = F(1,1)*tt1
tt7 = tt6+tt5+tt4
tt8 = tt7**((-8.0E+0)/3.0E+0)
tt9 = F(3,3)**2+F(3,2)**2+F(3,1)**2+F(2,3)**2+F(2,2)**2+F(2,1)**2&
&+F(1,3)**2+F(1,2)**2+F(1,1)**2
tt10 = tt7**((-5.0E+0)/3.0E+0)
tt11 = 2/tt7**(2.0E+0/3.0E+0)
tt12 = F(1,3)*F(3,2)-F(1,2)*F(3,3)
tt13 = (mu(1,1)*(((-4.0E+0)*F(2,1)*tt1*tt10)/3.0E+0+((-4.0E+0)*F(&
&1,1)*tt12*tt10)/3.0E+0+(1.0E+1*tt12*tt1*tt8*tt9)/9.0E+0))/2.0E+0+&
&lam(1,1)*tt12*tt1
tt14 = F(1,2)*F(2,3)-F(1,3)*F(2,2)
tt15 = (mu(1,1)*(((-4.0E+0)*F(3,1)*tt1*tt10)/3.0E+0+((-4.0E+0)*F(&
&1,1)*tt14*tt10)/3.0E+0+(1.0E+1*tt14*tt1*tt8*tt9)/9.0E+0))/2.0E+0+&
&lam(1,1)*tt14*tt1
tt16 = F(2,3)*F(3,1)-F(2,1)*F(3,3)
tt17 = (mu(1,1)*(((-4.0E+0)*F(1,2)*tt1*tt10)/3.0E+0+((-4.0E+0)*F(&
&1,1)*tt16*tt10)/3.0E+0+(1.0E+1*tt16*tt1*tt8*tt9)/9.0E+0))/2.0E+0+&
&lam(1,1)*tt16*tt1
tt18 = F(1,1)*F(3,3)-F(1,3)*F(3,1)
tt19 = tt6+tt5+tt4-1
tt20 = (mu(1,1)*(((-4.0E+0)*F(2,2)*tt1*tt10)/3.0E+0+((-4.0E+0)*F(&
&1,1)*tt18*tt10)/3.0E+0+((-2.0E+0)*F(3,3)*tt10*tt9)/3.0E+0+(1.0E+1&
&*tt18*tt1*tt8*tt9)/9.0E+0))/2.0E+0+lam(1,1)*F(3,3)*tt19+lam(1,1)*&
&tt18*tt1
tt21 = F(1,3)*F(2,1)-F(1,1)*F(2,3)
tt22 = (mu(1,1)*(((-4.0E+0)*F(3,2)*tt1*tt10)/3.0E+0+((-4.0E+0)*F(&
&1,1)*tt21*tt10)/3.0E+0+(2.0E+0*F(2,3)*tt10*tt9)/3.0E+0+(1.0E+1*tt&
&21*tt1*tt8*tt9)/9.0E+0))/2.0E+0-lam(1,1)*F(2,3)*tt19+lam(1,1)*tt2&
&1*tt1
tt23 = (mu(1,1)*(((-4.0E+0)*F(1,3)*tt1*tt10)/3.0E+0+((-4.0E+0)*F(&
&1,1)*tt3*tt10)/3.0E+0+(1.0E+1*tt3*tt1*tt8*tt9)/9.0E+0))/2.0E+0+la&
&m(1,1)*tt3*tt1
tt24 = F(1,2)*F(3,1)-F(1,1)*F(3,2)
tt25 = (mu(1,1)*(((-4.0E+0)*F(2,3)*tt1*tt10)/3.0E+0+((-4.0E+0)*F(&
&1,1)*tt24*tt10)/3.0E+0+(2.0E+0*F(3,2)*tt10*tt9)/3.0E+0+(1.0E+1*tt&
&24*tt1*tt8*tt9)/9.0E+0))/2.0E+0-lam(1,1)*F(3,2)*tt19+lam(1,1)*tt2&
&4*tt1
tt26 = F(1,1)*F(2,2)-F(1,2)*F(2,1)
tt27 = (mu(1,1)*(((-4.0E+0)*F(3,3)*tt1*tt10)/3.0E+0+((-4.0E+0)*F(&
&1,1)*tt26*tt10)/3.0E+0+((-2.0E+0)*F(2,2)*tt10*tt9)/3.0E+0+(1.0E+1&
&*tt26*tt1*tt8*tt9)/9.0E+0))/2.0E+0+lam(1,1)*F(2,2)*tt19+lam(1,1)*&
&tt26*tt1
tt28 = tt12**2
tt29 = (mu(1,1)*(((-4.0E+0)*F(3,1)*tt12*tt10)/3.0E+0+((-4.0E+0)*F&
&(2,1)*tt14*tt10)/3.0E+0+(1.0E+1*tt14*tt12*tt8*tt9)/9.0E+0))/2.0E+&
&0+lam(1,1)*tt14*tt12
tt30 = (mu(1,1)*(((-4.0E+0)*F(2,1)*tt16*tt10)/3.0E+0+((-4.0E+0)*F&
&(1,2)*tt12*tt10)/3.0E+0+(2.0E+0*F(3,3)*tt10*tt9)/3.0E+0+(1.0E+1*t&
&t12*tt16*tt8*tt9)/9.0E+0))/2.0E+0-lam(1,1)*F(3,3)*tt19+lam(1,1)*t&
&t12*tt16
tt31 = (mu(1,1)*(((-4.0E+0)*F(2,2)*tt12*tt10)/3.0E+0+((-4.0E+0)*F&
&(2,1)*tt18*tt10)/3.0E+0+(1.0E+1*tt18*tt12*tt8*tt9)/9.0E+0))/2.0E+&
&0+lam(1,1)*tt18*tt12
tt32 = (mu(1,1)*(((-4.0E+0)*F(3,2)*tt12*tt10)/3.0E+0+((-4.0E+0)*F&
&(2,1)*tt21*tt10)/3.0E+0+((-2.0E+0)*F(1,3)*tt10*tt9)/3.0E+0+(1.0E+&
&1*tt21*tt12*tt8*tt9)/9.0E+0))/2.0E+0+lam(1,1)*F(1,3)*tt19+lam(1,1&
&)*tt21*tt12
tt33 = (mu(1,1)*(((-4.0E+0)*F(1,3)*tt12*tt10)/3.0E+0+((-4.0E+0)*F&
&(2,1)*tt3*tt10)/3.0E+0+((-2.0E+0)*F(3,2)*tt10*tt9)/3.0E+0+(1.0E+1&
&*tt3*tt12*tt8*tt9)/9.0E+0))/2.0E+0+lam(1,1)*F(3,2)*tt19+lam(1,1)*&
&tt3*tt12
tt34 = (mu(1,1)*(((-4.0E+0)*F(2,3)*tt12*tt10)/3.0E+0+((-4.0E+0)*F&
&(2,1)*tt24*tt10)/3.0E+0+(1.0E+1*tt24*tt12*tt8*tt9)/9.0E+0))/2.0E+&
&0+lam(1,1)*tt24*tt12
tt35 = (mu(1,1)*(((-4.0E+0)*F(3,3)*tt12*tt10)/3.0E+0+((-4.0E+0)*F&
&(2,1)*tt26*tt10)/3.0E+0+(2.0E+0*F(1,2)*tt10*tt9)/3.0E+0+(1.0E+1*t&
&t26*tt12*tt8*tt9)/9.0E+0))/2.0E+0-lam(1,1)*F(1,2)*tt19+lam(1,1)*t&
&t26*tt12
tt36 = tt14**2
tt37 = (mu(1,1)*(((-4.0E+0)*F(3,1)*tt16*tt10)/3.0E+0+((-4.0E+0)*F&
&(1,2)*tt14*tt10)/3.0E+0+((-2.0E+0)*F(2,3)*tt10*tt9)/3.0E+0+(1.0E+&
&1*tt14*tt16*tt8*tt9)/9.0E+0))/2.0E+0+lam(1,1)*F(2,3)*tt19+lam(1,1&
&)*tt14*tt16
tt38 = (mu(1,1)*(((-4.0E+0)*F(3,1)*tt18*tt10)/3.0E+0+((-4.0E+0)*F&
&(2,2)*tt14*tt10)/3.0E+0+(2.0E+0*F(1,3)*tt10*tt9)/3.0E+0+(1.0E+1*t&
&t14*tt18*tt8*tt9)/9.0E+0))/2.0E+0-lam(1,1)*F(1,3)*tt19+lam(1,1)*t&
&t14*tt18
tt39 = (mu(1,1)*(((-4.0E+0)*tt14*F(3,2)*tt10)/3.0E+0+((-4.0E+0)*t&
&t21*F(3,1)*tt10)/3.0E+0+(1.0E+1*tt21*tt14*tt8*tt9)/9.0E+0))/2.0E+&
&0+lam(1,1)*tt21*tt14
tt40 = (mu(1,1)*(((-4.0E+0)*F(3,1)*tt3*tt10)/3.0E+0+((-4.0E+0)*F(&
&1,3)*tt14*tt10)/3.0E+0+(2.0E+0*F(2,2)*tt10*tt9)/3.0E+0+(1.0E+1*tt&
&14*tt3*tt8*tt9)/9.0E+0))/2.0E+0-lam(1,1)*F(2,2)*tt19+lam(1,1)*tt1&
&4*tt3
tt41 = (mu(1,1)*(((-4.0E+0)*F(3,1)*tt24*tt10)/3.0E+0+((-4.0E+0)*F&
&(2,3)*tt14*tt10)/3.0E+0+((-2.0E+0)*F(1,2)*tt10*tt9)/3.0E+0+(1.0E+&
&1*tt14*tt24*tt8*tt9)/9.0E+0))/2.0E+0+lam(1,1)*F(1,2)*tt19+lam(1,1&
&)*tt14*tt24
tt42 = (mu(1,1)*(((-4.0E+0)*tt14*F(3,3)*tt10)/3.0E+0+((-4.0E+0)*t&
&t26*F(3,1)*tt10)/3.0E+0+(1.0E+1*tt26*tt14*tt8*tt9)/9.0E+0))/2.0E+&
&0+lam(1,1)*tt26*tt14
tt43 = tt16**2
tt44 = (mu(1,1)*(((-4.0E+0)*F(2,2)*tt16*tt10)/3.0E+0+((-4.0E+0)*F&
&(1,2)*tt18*tt10)/3.0E+0+(1.0E+1*tt18*tt16*tt8*tt9)/9.0E+0))/2.0E+&
&0+lam(1,1)*tt18*tt16
tt45 = (mu(1,1)*(((-4.0E+0)*F(3,2)*tt16*tt10)/3.0E+0+((-4.0E+0)*F&
&(1,2)*tt21*tt10)/3.0E+0+(1.0E+1*tt21*tt16*tt8*tt9)/9.0E+0))/2.0E+&
&0+lam(1,1)*tt21*tt16
tt46 = (mu(1,1)*(((-4.0E+0)*F(1,3)*tt16*tt10)/3.0E+0+((-4.0E+0)*F&
&(1,2)*tt3*tt10)/3.0E+0+(1.0E+1*tt3*tt16*tt8*tt9)/9.0E+0))/2.0E+0+&
&lam(1,1)*tt3*tt16
tt47 = (mu(1,1)*(((-4.0E+0)*F(2,3)*tt16*tt10)/3.0E+0+((-4.0E+0)*F&
&(1,2)*tt24*tt10)/3.0E+0+((-2.0E+0)*F(3,1)*tt10*tt9)/3.0E+0+(1.0E+&
&1*tt24*tt16*tt8*tt9)/9.0E+0))/2.0E+0+lam(1,1)*F(3,1)*tt19+lam(1,1&
&)*tt24*tt16
tt48 = (mu(1,1)*(((-4.0E+0)*F(3,3)*tt16*tt10)/3.0E+0+((-4.0E+0)*F&
&(1,2)*tt26*tt10)/3.0E+0+(2.0E+0*F(2,1)*tt10*tt9)/3.0E+0+(1.0E+1*t&
&t26*tt16*tt8*tt9)/9.0E+0))/2.0E+0-lam(1,1)*F(2,1)*tt19+lam(1,1)*t&
&t26*tt16
tt49 = tt18**2
tt50 = (mu(1,1)*(((-4.0E+0)*F(3,2)*tt18*tt10)/3.0E+0+((-4.0E+0)*F&
&(2,2)*tt21*tt10)/3.0E+0+(1.0E+1*tt21*tt18*tt8*tt9)/9.0E+0))/2.0E+&
&0+lam(1,1)*tt21*tt18
tt51 = (mu(1,1)*(((-4.0E+0)*F(1,3)*tt18*tt10)/3.0E+0+((-4.0E+0)*F&
&(2,2)*tt3*tt10)/3.0E+0+(2.0E+0*F(3,1)*tt10*tt9)/3.0E+0+(1.0E+1*tt&
&3*tt18*tt8*tt9)/9.0E+0))/2.0E+0-lam(1,1)*F(3,1)*tt19+lam(1,1)*tt3&
&*tt18
tt52 = (mu(1,1)*(((-4.0E+0)*F(2,3)*tt18*tt10)/3.0E+0+((-4.0E+0)*F&
&(2,2)*tt24*tt10)/3.0E+0+(1.0E+1*tt24*tt18*tt8*tt9)/9.0E+0))/2.0E+&
&0+lam(1,1)*tt24*tt18
tt53 = (mu(1,1)*(((-4.0E+0)*F(3,3)*tt18*tt10)/3.0E+0+((-4.0E+0)*F&
&(2,2)*tt26*tt10)/3.0E+0+((-2.0E+0)*F(1,1)*tt10*tt9)/3.0E+0+(1.0E+&
&1*tt26*tt18*tt8*tt9)/9.0E+0))/2.0E+0+F(1,1)*lam(1,1)*tt19+lam(1,1&
&)*tt26*tt18
tt54 = tt21**2
tt55 = (mu(1,1)*(((-4.0E+0)*F(3,2)*tt3*tt10)/3.0E+0+((-4.0E+0)*F(&
&1,3)*tt21*tt10)/3.0E+0+((-2.0E+0)*F(2,1)*tt10*tt9)/3.0E+0+(1.0E+1&
&*tt21*tt3*tt8*tt9)/9.0E+0))/2.0E+0+lam(1,1)*F(2,1)*tt19+lam(1,1)*&
&tt21*tt3
tt56 = (mu(1,1)*(((-4.0E+0)*F(3,2)*tt24*tt10)/3.0E+0+((-4.0E+0)*F&
&(2,3)*tt21*tt10)/3.0E+0+(2.0E+0*F(1,1)*tt10*tt9)/3.0E+0+(1.0E+1*t&
&t21*tt24*tt8*tt9)/9.0E+0))/2.0E+0-F(1,1)*lam(1,1)*tt19+lam(1,1)*t&
&t21*tt24
tt57 = (mu(1,1)*(((-4.0E+0)*tt21*F(3,3)*tt10)/3.0E+0+((-4.0E+0)*t&
&t26*F(3,2)*tt10)/3.0E+0+(1.0E+1*tt26*tt21*tt8*tt9)/9.0E+0))/2.0E+&
&0+lam(1,1)*tt26*tt21
tt58 = tt3**2
tt59 = (mu(1,1)*(((-4.0E+0)*F(2,3)*tt3*tt10)/3.0E+0+((-4.0E+0)*F(&
&1,3)*tt24*tt10)/3.0E+0+(1.0E+1*tt24*tt3*tt8*tt9)/9.0E+0))/2.0E+0+&
&lam(1,1)*tt24*tt3
tt60 = (mu(1,1)*(((-4.0E+0)*tt3*F(3,3)*tt10)/3.0E+0+((-4.0E+0)*F(&
&1,3)*tt26*tt10)/3.0E+0+(1.0E+1*tt26*tt3*tt8*tt9)/9.0E+0))/2.0E+0+&
&lam(1,1)*tt26*tt3
tt61 = tt24**2
tt62 = (mu(1,1)*(((-4.0E+0)*tt24*F(3,3)*tt10)/3.0E+0+((-4.0E+0)*t&
&t26*F(2,3)*tt10)/3.0E+0+(1.0E+1*tt26*tt24*tt8*tt9)/9.0E+0))/2.0E+&
&0+lam(1,1)*tt26*tt24
tt63 = tt26**2
hes(1,1) = (mu(1,1)*(tt11+((-8.0E+0)*F(1,1)*tt1*tt10)/3.0E+0+(1.0&
&E+1*tt2*tt8*tt9)/9.0E+0))/2.0E+0+lam(1,1)*tt2
hes(1,2) = tt13
hes(1,3) = tt15
hes(1,4) = tt17
hes(1,5) = tt20
hes(1,6) = tt22
hes(1,7) = tt23
hes(1,8) = tt25
hes(1,9) = tt27
hes(2,1) = tt13
hes(2,2) = (mu(1,1)*(tt11+((-8.0E+0)*F(2,1)*tt12*tt10)/3.0E+0+(1.&
&0E+1*tt28*tt8*tt9)/9.0E+0))/2.0E+0+lam(1,1)*tt28
hes(2,3) = tt29
hes(2,4) = tt30
hes(2,5) = tt31
hes(2,6) = tt32
hes(2,7) = tt33
hes(2,8) = tt34
hes(2,9) = tt35
hes(3,1) = tt15
hes(3,2) = tt29
hes(3,3) = (mu(1,1)*(tt11+((-8.0E+0)*tt14*F(3,1)*tt10)/3.0E+0+(1.&
&0E+1*tt36*tt8*tt9)/9.0E+0))/2.0E+0+lam(1,1)*tt36
hes(3,4) = tt37
hes(3,5) = tt38
hes(3,6) = tt39
hes(3,7) = tt40
hes(3,8) = tt41
hes(3,9) = tt42
hes(4,1) = tt17
hes(4,2) = tt30
hes(4,3) = tt37
hes(4,4) = (mu(1,1)*(tt11+((-8.0E+0)*F(1,2)*tt16*tt10)/3.0E+0+(1.&
&0E+1*tt43*tt8*tt9)/9.0E+0))/2.0E+0+lam(1,1)*tt43
hes(4,5) = tt44
hes(4,6) = tt45
hes(4,7) = tt46
hes(4,8) = tt47
hes(4,9) = tt48
hes(5,1) = tt20
hes(5,2) = tt31
hes(5,3) = tt38
hes(5,4) = tt44
hes(5,5) = (mu(1,1)*(tt11+((-8.0E+0)*F(2,2)*tt18*tt10)/3.0E+0+(1.&
&0E+1*tt49*tt8*tt9)/9.0E+0))/2.0E+0+lam(1,1)*tt49
hes(5,6) = tt50
hes(5,7) = tt51
hes(5,8) = tt52
hes(5,9) = tt53
hes(6,1) = tt22
hes(6,2) = tt32
hes(6,3) = tt39
hes(6,4) = tt45
hes(6,5) = tt50
hes(6,6) = (mu(1,1)*(tt11+((-8.0E+0)*tt21*F(3,2)*tt10)/3.0E+0+(1.&
&0E+1*tt54*tt8*tt9)/9.0E+0))/2.0E+0+lam(1,1)*tt54
hes(6,7) = tt55
hes(6,8) = tt56
hes(6,9) = tt57
hes(7,1) = tt23
hes(7,2) = tt33
hes(7,3) = tt40
hes(7,4) = tt46
hes(7,5) = tt51
hes(7,6) = tt55
hes(7,7) = (mu(1,1)*(tt11+((-8.0E+0)*F(1,3)*tt3*tt10)/3.0E+0+(1.0&
&E+1*tt58*tt8*tt9)/9.0E+0))/2.0E+0+lam(1,1)*tt58
hes(7,8) = tt59
hes(7,9) = tt60
hes(8,1) = tt25
hes(8,2) = tt34
hes(8,3) = tt41
hes(8,4) = tt47
hes(8,5) = tt52
hes(8,6) = tt56
hes(8,7) = tt59
hes(8,8) = (mu(1,1)*(tt11+((-8.0E+0)*F(2,3)*tt24*tt10)/3.0E+0+(1.&
&0E+1*tt61*tt8*tt9)/9.0E+0))/2.0E+0+lam(1,1)*tt61
hes(8,9) = tt62
hes(9,1) = tt27
hes(9,2) = tt35
hes(9,3) = tt42
hes(9,4) = tt48
hes(9,5) = tt53
hes(9,6) = tt57
hes(9,7) = tt60
hes(9,8) = tt62
hes(9,9) = (mu(1,1)*(tt11+((-8.0E+0)*tt26*F(3,3)*tt10)/3.0E+0+(1.&
&0E+1*tt63*tt8*tt9)/9.0E+0))/2.0E+0+lam(1,1)*tt63
END 
SUBROUTINE basis_mini_trace(val, N, K) 
IMPLICIT NONE 
REAL(KIND=8) val(1, 1) 
REAL(KIND=8) N(3, 6) 
REAL(KIND=8) K(3, 3) 
val(1,1) = N(3,3)*(K(3,3)*N(3,6)+N(2,6)*K(3,2)+N(1,6)*K(3,1))+N(2&
&,3)*(K(2,3)*N(3,6)+K(2,2)*N(2,6)+N(1,6)*K(2,1))+N(1,3)*(K(1,3)*N(&
&3,6)+K(1,2)*N(2,6)+K(1,1)*N(1,6))+N(3,2)*(K(3,3)*N(3,5)+N(2,5)*K(&
&3,2)+N(1,5)*K(3,1))+N(2,2)*(K(2,3)*N(3,5)+K(2,2)*N(2,5)+N(1,5)*K(&
&2,1))+N(1,2)*(K(1,3)*N(3,5)+K(1,2)*N(2,5)+K(1,1)*N(1,5))+N(3,1)*(&
&K(3,3)*N(3,4)+N(2,4)*K(3,2)+N(1,4)*K(3,1))+N(2,1)*(K(2,3)*N(3,4)+&
&K(2,2)*N(2,4)+N(1,4)*K(2,1))+N(1,1)*(K(1,3)*N(3,4)+K(1,2)*N(2,4)+&
&K(1,1)*N(1,4))
END 
SUBROUTINE basis_mini_trace_jac(jac, N, K) 
IMPLICIT NONE 
REAL(KIND=8) jac(1, 18) 
REAL(KIND=8) N(3, 6) 
REAL(KIND=8) K(3, 3) 
jac(1,1) = K(1,3)*N(3,4)+K(1,2)*N(2,4)+K(1,1)*N(1,4)
jac(1,2) = K(2,3)*N(3,4)+K(2,2)*N(2,4)+N(1,4)*K(2,1)
jac(1,3) = K(3,3)*N(3,4)+N(2,4)*K(3,2)+N(1,4)*K(3,1)
jac(1,4) = K(1,3)*N(3,5)+K(1,2)*N(2,5)+K(1,1)*N(1,5)
jac(1,5) = K(2,3)*N(3,5)+K(2,2)*N(2,5)+N(1,5)*K(2,1)
jac(1,6) = K(3,3)*N(3,5)+N(2,5)*K(3,2)+N(1,5)*K(3,1)
jac(1,7) = K(1,3)*N(3,6)+K(1,2)*N(2,6)+K(1,1)*N(1,6)
jac(1,8) = K(2,3)*N(3,6)+K(2,2)*N(2,6)+N(1,6)*K(2,1)
jac(1,9) = K(3,3)*N(3,6)+N(2,6)*K(3,2)+N(1,6)*K(3,1)
jac(1,10) = K(3,1)*N(3,1)+K(2,1)*N(2,1)+K(1,1)*N(1,1)
jac(1,11) = N(3,1)*K(3,2)+N(2,1)*K(2,2)+N(1,1)*K(1,2)
jac(1,12) = N(3,1)*K(3,3)+N(2,1)*K(2,3)+N(1,1)*K(1,3)
jac(1,13) = K(3,1)*N(3,2)+K(2,1)*N(2,2)+K(1,1)*N(1,2)
jac(1,14) = K(3,2)*N(3,2)+K(2,2)*N(2,2)+K(1,2)*N(1,2)
jac(1,15) = N(3,2)*K(3,3)+N(2,2)*K(2,3)+N(1,2)*K(1,3)
jac(1,16) = K(3,1)*N(3,3)+K(2,1)*N(2,3)+K(1,1)*N(1,3)
jac(1,17) = K(3,2)*N(3,3)+K(2,2)*N(2,3)+K(1,2)*N(1,3)
jac(1,18) = K(3,3)*N(3,3)+K(2,3)*N(2,3)+K(1,3)*N(1,3)
END 
SUBROUTINE basis_mini_trace_hes(hes, N, K) 
IMPLICIT NONE 
REAL(KIND=8) hes(18, 18) 
REAL(KIND=8) N(3, 6) 
REAL(KIND=8) K(3, 3) 
hes(1,1) = 0
hes(1,2) = 0
hes(1,3) = 0
hes(1,4) = 0
hes(1,5) = 0
hes(1,6) = 0
hes(1,7) = 0
hes(1,8) = 0
hes(1,9) = 0
hes(1,10) = K(1,1)
hes(1,11) = K(1,2)
hes(1,12) = K(1,3)
hes(1,13) = 0
hes(1,14) = 0
hes(1,15) = 0
hes(1,16) = 0
hes(1,17) = 0
hes(1,18) = 0
hes(2,1) = 0
hes(2,2) = 0
hes(2,3) = 0
hes(2,4) = 0
hes(2,5) = 0
hes(2,6) = 0
hes(2,7) = 0
hes(2,8) = 0
hes(2,9) = 0
hes(2,10) = K(2,1)
hes(2,11) = K(2,2)
hes(2,12) = K(2,3)
hes(2,13) = 0
hes(2,14) = 0
hes(2,15) = 0
hes(2,16) = 0
hes(2,17) = 0
hes(2,18) = 0
hes(3,1) = 0
hes(3,2) = 0
hes(3,3) = 0
hes(3,4) = 0
hes(3,5) = 0
hes(3,6) = 0
hes(3,7) = 0
hes(3,8) = 0
hes(3,9) = 0
hes(3,10) = K(3,1)
hes(3,11) = K(3,2)
hes(3,12) = K(3,3)
hes(3,13) = 0
hes(3,14) = 0
hes(3,15) = 0
hes(3,16) = 0
hes(3,17) = 0
hes(3,18) = 0
hes(4,1) = 0
hes(4,2) = 0
hes(4,3) = 0
hes(4,4) = 0
hes(4,5) = 0
hes(4,6) = 0
hes(4,7) = 0
hes(4,8) = 0
hes(4,9) = 0
hes(4,10) = 0
hes(4,11) = 0
hes(4,12) = 0
hes(4,13) = K(1,1)
hes(4,14) = K(1,2)
hes(4,15) = K(1,3)
hes(4,16) = 0
hes(4,17) = 0
hes(4,18) = 0
hes(5,1) = 0
hes(5,2) = 0
hes(5,3) = 0
hes(5,4) = 0
hes(5,5) = 0
hes(5,6) = 0
hes(5,7) = 0
hes(5,8) = 0
hes(5,9) = 0
hes(5,10) = 0
hes(5,11) = 0
hes(5,12) = 0
hes(5,13) = K(2,1)
hes(5,14) = K(2,2)
hes(5,15) = K(2,3)
hes(5,16) = 0
hes(5,17) = 0
hes(5,18) = 0
hes(6,1) = 0
hes(6,2) = 0
hes(6,3) = 0
hes(6,4) = 0
hes(6,5) = 0
hes(6,6) = 0
hes(6,7) = 0
hes(6,8) = 0
hes(6,9) = 0
hes(6,10) = 0
hes(6,11) = 0
hes(6,12) = 0
hes(6,13) = K(3,1)
hes(6,14) = K(3,2)
hes(6,15) = K(3,3)
hes(6,16) = 0
hes(6,17) = 0
hes(6,18) = 0
hes(7,1) = 0
hes(7,2) = 0
hes(7,3) = 0
hes(7,4) = 0
hes(7,5) = 0
hes(7,6) = 0
hes(7,7) = 0
hes(7,8) = 0
hes(7,9) = 0
hes(7,10) = 0
hes(7,11) = 0
hes(7,12) = 0
hes(7,13) = 0
hes(7,14) = 0
hes(7,15) = 0
hes(7,16) = K(1,1)
hes(7,17) = K(1,2)
hes(7,18) = K(1,3)
hes(8,1) = 0
hes(8,2) = 0
hes(8,3) = 0
hes(8,4) = 0
hes(8,5) = 0
hes(8,6) = 0
hes(8,7) = 0
hes(8,8) = 0
hes(8,9) = 0
hes(8,10) = 0
hes(8,11) = 0
hes(8,12) = 0
hes(8,13) = 0
hes(8,14) = 0
hes(8,15) = 0
hes(8,16) = K(2,1)
hes(8,17) = K(2,2)
hes(8,18) = K(2,3)
hes(9,1) = 0
hes(9,2) = 0
hes(9,3) = 0
hes(9,4) = 0
hes(9,5) = 0
hes(9,6) = 0
hes(9,7) = 0
hes(9,8) = 0
hes(9,9) = 0
hes(9,10) = 0
hes(9,11) = 0
hes(9,12) = 0
hes(9,13) = 0
hes(9,14) = 0
hes(9,15) = 0
hes(9,16) = K(3,1)
hes(9,17) = K(3,2)
hes(9,18) = K(3,3)
hes(10,1) = K(1,1)
hes(10,2) = K(2,1)
hes(10,3) = K(3,1)
hes(10,4) = 0
hes(10,5) = 0
hes(10,6) = 0
hes(10,7) = 0
hes(10,8) = 0
hes(10,9) = 0
hes(10,10) = 0
hes(10,11) = 0
hes(10,12) = 0
hes(10,13) = 0
hes(10,14) = 0
hes(10,15) = 0
hes(10,16) = 0
hes(10,17) = 0
hes(10,18) = 0
hes(11,1) = K(1,2)
hes(11,2) = K(2,2)
hes(11,3) = K(3,2)
hes(11,4) = 0
hes(11,5) = 0
hes(11,6) = 0
hes(11,7) = 0
hes(11,8) = 0
hes(11,9) = 0
hes(11,10) = 0
hes(11,11) = 0
hes(11,12) = 0
hes(11,13) = 0
hes(11,14) = 0
hes(11,15) = 0
hes(11,16) = 0
hes(11,17) = 0
hes(11,18) = 0
hes(12,1) = K(1,3)
hes(12,2) = K(2,3)
hes(12,3) = K(3,3)
hes(12,4) = 0
hes(12,5) = 0
hes(12,6) = 0
hes(12,7) = 0
hes(12,8) = 0
hes(12,9) = 0
hes(12,10) = 0
hes(12,11) = 0
hes(12,12) = 0
hes(12,13) = 0
hes(12,14) = 0
hes(12,15) = 0
hes(12,16) = 0
hes(12,17) = 0
hes(12,18) = 0
hes(13,1) = 0
hes(13,2) = 0
hes(13,3) = 0
hes(13,4) = K(1,1)
hes(13,5) = K(2,1)
hes(13,6) = K(3,1)
hes(13,7) = 0
hes(13,8) = 0
hes(13,9) = 0
hes(13,10) = 0
hes(13,11) = 0
hes(13,12) = 0
hes(13,13) = 0
hes(13,14) = 0
hes(13,15) = 0
hes(13,16) = 0
hes(13,17) = 0
hes(13,18) = 0
hes(14,1) = 0
hes(14,2) = 0
hes(14,3) = 0
hes(14,4) = K(1,2)
hes(14,5) = K(2,2)
hes(14,6) = K(3,2)
hes(14,7) = 0
hes(14,8) = 0
hes(14,9) = 0
hes(14,10) = 0
hes(14,11) = 0
hes(14,12) = 0
hes(14,13) = 0
hes(14,14) = 0
hes(14,15) = 0
hes(14,16) = 0
hes(14,17) = 0
hes(14,18) = 0
hes(15,1) = 0
hes(15,2) = 0
hes(15,3) = 0
hes(15,4) = K(1,3)
hes(15,5) = K(2,3)
hes(15,6) = K(3,3)
hes(15,7) = 0
hes(15,8) = 0
hes(15,9) = 0
hes(15,10) = 0
hes(15,11) = 0
hes(15,12) = 0
hes(15,13) = 0
hes(15,14) = 0
hes(15,15) = 0
hes(15,16) = 0
hes(15,17) = 0
hes(15,18) = 0
hes(16,1) = 0
hes(16,2) = 0
hes(16,3) = 0
hes(16,4) = 0
hes(16,5) = 0
hes(16,6) = 0
hes(16,7) = K(1,1)
hes(16,8) = K(2,1)
hes(16,9) = K(3,1)
hes(16,10) = 0
hes(16,11) = 0
hes(16,12) = 0
hes(16,13) = 0
hes(16,14) = 0
hes(16,15) = 0
hes(16,16) = 0
hes(16,17) = 0
hes(16,18) = 0
hes(17,1) = 0
hes(17,2) = 0
hes(17,3) = 0
hes(17,4) = 0
hes(17,5) = 0
hes(17,6) = 0
hes(17,7) = K(1,2)
hes(17,8) = K(2,2)
hes(17,9) = K(3,2)
hes(17,10) = 0
hes(17,11) = 0
hes(17,12) = 0
hes(17,13) = 0
hes(17,14) = 0
hes(17,15) = 0
hes(17,16) = 0
hes(17,17) = 0
hes(17,18) = 0
hes(18,1) = 0
hes(18,2) = 0
hes(18,3) = 0
hes(18,4) = 0
hes(18,5) = 0
hes(18,6) = 0
hes(18,7) = K(1,3)
hes(18,8) = K(2,3)
hes(18,9) = K(3,3)
hes(18,10) = 0
hes(18,11) = 0
hes(18,12) = 0
hes(18,13) = 0
hes(18,14) = 0
hes(18,15) = 0
hes(18,16) = 0
hes(18,17) = 0
hes(18,18) = 0
END 
