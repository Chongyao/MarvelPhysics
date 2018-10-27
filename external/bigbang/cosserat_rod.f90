SUBROUTINE rod_stretch(val, X, d, Es, r) 
IMPLICIT NONE 
REAL(KIND=8) val(1, 1) 
REAL(KIND=8) X(3, 2) 
REAL(KIND=8) d(1, 1) 
REAL(KIND=8) Es(1, 1) 
REAL(KIND=8) r(1, 1) 
val(1,1) = 1.5707963267948966E+0*Es(1,1)*d(1,1)*r(1,1)**2*(sqrt((&
&X(3,2)-X(3,1))**2+(X(2,2)-X(2,1))**2+(X(1,2)-X(1,1))**2)/d(1,1)-1&
&)**2
END 
SUBROUTINE rod_stretch_jac(jac, X, d, Es, r) 
IMPLICIT NONE 
REAL(KIND=8) jac(1, 6) 
REAL(KIND=8) X(3, 2) 
REAL(KIND=8) d(1, 1) 
REAL(KIND=8) Es(1, 1) 
REAL(KIND=8) r(1, 1) 
REAL(KIND=8)  tt1 
REAL(KIND=8)  tt2 
REAL(KIND=8)  tt3 
REAL(KIND=8)  tt4 
REAL(KIND=8)  tt5 
REAL(KIND=8)  tt6 
REAL(KIND=8)  tt7 
tt1 = r(1,1)**2
tt2 = X(1,2)-X(1,1)
tt3 = X(2,2)-X(2,1)
tt4 = X(3,2)-X(3,1)
tt5 = sqrt(tt4**2+tt3**2+tt2**2)
tt6 = 1/tt5
tt7 = tt5/d(1,1)-1
jac(1,1) = -3.141592653589793E+0*Es(1,1)*tt1*tt2*tt6*tt7
jac(1,2) = -3.141592653589793E+0*Es(1,1)*tt1*tt3*tt6*tt7
jac(1,3) = -3.141592653589793E+0*Es(1,1)*tt1*tt4*tt6*tt7
jac(1,4) = 3.141592653589793E+0*Es(1,1)*tt1*tt2*tt6*tt7
jac(1,5) = 3.141592653589793E+0*Es(1,1)*tt1*tt3*tt6*tt7
jac(1,6) = 3.141592653589793E+0*Es(1,1)*tt1*tt4*tt6*tt7
END 
SUBROUTINE rod_stretch_hes(hes, X, d, Es, r) 
IMPLICIT NONE 
REAL(KIND=8) hes(6, 6) 
REAL(KIND=8) X(3, 2) 
REAL(KIND=8) d(1, 1) 
REAL(KIND=8) Es(1, 1) 
REAL(KIND=8) r(1, 1) 
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
tt1 = 1/d(1,1)
tt2 = r(1,1)**2
tt3 = X(1,2)-X(1,1)
tt4 = tt3**2
tt5 = X(2,2)-X(2,1)
tt6 = tt5**2
tt7 = X(3,2)-X(3,1)
tt8 = tt7**2
tt9 = tt8+tt6+tt4
tt10 = 1/tt9
tt11 = sqrt(tt9)
tt12 = 1/tt11**3
tt13 = tt1*tt11-1
tt14 = 1/tt11
tt15 = 3.141592653589793E+0*Es(1,1)*tt2*tt14*tt13
tt16 = tt15-3.141592653589793E+0*Es(1,1)*tt2*tt4*tt12*tt13+3.1415&
&92653589793E+0*Es(1,1)*tt1*tt2*tt4*tt10
tt17 = 3.141592653589793E+0*Es(1,1)*tt1*tt2*tt3*tt5*tt10-3.141592&
&653589793E+0*Es(1,1)*tt2*tt3*tt5*tt12*tt13
tt18 = 3.141592653589793E+0*Es(1,1)*tt1*tt2*tt3*tt7*tt10-3.141592&
&653589793E+0*Es(1,1)*tt2*tt3*tt7*tt12*tt13
tt19 = -3.141592653589793E+0*Es(1,1)*tt2*tt14*tt13
tt20 = tt19+3.141592653589793E+0*Es(1,1)*tt2*tt4*tt12*tt13-3.1415&
&92653589793E+0*Es(1,1)*tt1*tt2*tt4*tt10
tt21 = 3.141592653589793E+0*Es(1,1)*tt2*tt3*tt5*tt12*tt13-3.14159&
&2653589793E+0*Es(1,1)*tt1*tt2*tt3*tt5*tt10
tt22 = 3.141592653589793E+0*Es(1,1)*tt2*tt3*tt7*tt12*tt13-3.14159&
&2653589793E+0*Es(1,1)*tt1*tt2*tt3*tt7*tt10
tt23 = tt15-3.141592653589793E+0*Es(1,1)*tt2*tt6*tt12*tt13+3.1415&
&92653589793E+0*Es(1,1)*tt1*tt2*tt6*tt10
tt24 = 3.141592653589793E+0*Es(1,1)*tt1*tt2*tt5*tt7*tt10-3.141592&
&653589793E+0*Es(1,1)*tt2*tt5*tt7*tt12*tt13
tt25 = tt19+3.141592653589793E+0*Es(1,1)*tt2*tt6*tt12*tt13-3.1415&
&92653589793E+0*Es(1,1)*tt1*tt2*tt6*tt10
tt26 = 3.141592653589793E+0*Es(1,1)*tt2*tt5*tt7*tt12*tt13-3.14159&
&2653589793E+0*Es(1,1)*tt1*tt2*tt5*tt7*tt10
tt27 = tt15-3.141592653589793E+0*Es(1,1)*tt2*tt8*tt12*tt13+3.1415&
&92653589793E+0*Es(1,1)*tt1*tt2*tt8*tt10
tt28 = tt19+3.141592653589793E+0*Es(1,1)*tt2*tt8*tt12*tt13-3.1415&
&92653589793E+0*Es(1,1)*tt1*tt2*tt8*tt10
hes(1,1) = tt16
hes(1,2) = tt17
hes(1,3) = tt18
hes(1,4) = tt20
hes(1,5) = tt21
hes(1,6) = tt22
hes(2,1) = tt17
hes(2,2) = tt23
hes(2,3) = tt24
hes(2,4) = tt21
hes(2,5) = tt25
hes(2,6) = tt26
hes(3,1) = tt18
hes(3,2) = tt24
hes(3,3) = tt27
hes(3,4) = tt22
hes(3,5) = tt26
hes(3,6) = tt28
hes(4,1) = tt20
hes(4,2) = tt21
hes(4,3) = tt22
hes(4,4) = tt16
hes(4,5) = tt17
hes(4,6) = tt18
hes(5,1) = tt21
hes(5,2) = tt25
hes(5,3) = tt26
hes(5,4) = tt17
hes(5,5) = tt23
hes(5,6) = tt24
hes(6,1) = tt22
hes(6,2) = tt26
hes(6,3) = tt28
hes(6,4) = tt18
hes(6,5) = tt24
hes(6,6) = tt27
END 
SUBROUTINE rod_bend(val, Q, u, d, E, G, r) 
IMPLICIT NONE 
REAL(KIND=8) val(1, 1) 
REAL(KIND=8) Q(4, 2) 
REAL(KIND=8) u(3, 1) 
REAL(KIND=8) d(1, 1) 
REAL(KIND=8) E(1, 1) 
REAL(KIND=8) G(1, 1) 
REAL(KIND=8) r(1, 1) 
REAL(KIND=8)  tt1 
REAL(KIND=8)  tt2 
tt1 = r(1,1)**2
tt2 = 1/d(1,1)
val(1,1) = 7.853981633974483E-1*E(1,1)*d(1,1)*tt1*((-tt2*(2*Q(3,1&
&)*Q(4,2)-2*Q(3,2)*Q(4,1)+2*Q(1,1)*Q(2,2)-2*Q(1,2)*Q(2,1)))-u(3,1)&
&)**2+3.9269908169872416E-1*E(1,1)*d(1,1)*tt1*((-tt2*(2*Q(2,1)*Q(4&
&,2)-2*Q(2,2)*Q(4,1)-2*Q(1,1)*Q(3,2)+2*Q(1,2)*Q(3,1)))-u(2,1))**2+&
&3.9269908169872416E-1*E(1,1)*d(1,1)*tt1*((-tt2*(2*Q(1,1)*Q(4,2)-2&
&*Q(1,2)*Q(4,1)+2*Q(2,1)*Q(3,2)-2*Q(2,2)*Q(3,1)))-u(1,1))**2
END 
SUBROUTINE rod_bend_jac(jac, Q, u, d, E, G, r) 
IMPLICIT NONE 
REAL(KIND=8) jac(1, 8) 
REAL(KIND=8) Q(4, 2) 
REAL(KIND=8) u(3, 1) 
REAL(KIND=8) d(1, 1) 
REAL(KIND=8) E(1, 1) 
REAL(KIND=8) G(1, 1) 
REAL(KIND=8) r(1, 1) 
REAL(KIND=8)  tt1 
REAL(KIND=8)  tt2 
REAL(KIND=8)  tt3 
REAL(KIND=8)  tt4 
REAL(KIND=8)  tt5 
tt1 = r(1,1)**2
tt2 = 1/d(1,1)
tt3 = (-tt2*(2*Q(1,1)*Q(4,2)-2*Q(1,2)*Q(4,1)+2*Q(2,1)*Q(3,2)-2*Q(&
&2,2)*Q(3,1)))-u(1,1)
tt4 = (-tt2*(2*Q(2,1)*Q(4,2)-2*Q(2,2)*Q(4,1)-2*Q(1,1)*Q(3,2)+2*Q(&
&1,2)*Q(3,1)))-u(2,1)
tt5 = (-tt2*(2*Q(3,1)*Q(4,2)-2*Q(3,2)*Q(4,1)+2*Q(1,1)*Q(2,2)-2*Q(&
&1,2)*Q(2,1)))-u(3,1)
jac(1,1) = (-3.141592653589793E+0*E(1,1)*tt1*Q(2,2)*tt5)+1.570796&
&3267948966E+0*E(1,1)*tt1*Q(3,2)*tt4-1.5707963267948966E+0*E(1,1)*&
&tt1*Q(4,2)*tt3
jac(1,2) = 3.141592653589793E+0*E(1,1)*tt1*Q(1,2)*tt5-1.570796326&
&7948966E+0*E(1,1)*tt1*Q(4,2)*tt4-1.5707963267948966E+0*E(1,1)*tt1&
&*Q(3,2)*tt3
jac(1,3) = (-3.141592653589793E+0*E(1,1)*tt1*Q(4,2)*tt5)-1.570796&
&3267948966E+0*E(1,1)*tt1*Q(1,2)*tt4+1.5707963267948966E+0*E(1,1)*&
&tt1*Q(2,2)*tt3
jac(1,4) = 3.141592653589793E+0*E(1,1)*tt1*Q(3,2)*tt5+1.570796326&
&7948966E+0*E(1,1)*tt1*Q(2,2)*tt4+1.5707963267948966E+0*E(1,1)*tt1&
&*Q(1,2)*tt3
jac(1,5) = 3.141592653589793E+0*E(1,1)*tt1*Q(2,1)*tt5-1.570796326&
&7948966E+0*E(1,1)*tt1*Q(3,1)*tt4+1.5707963267948966E+0*E(1,1)*tt1&
&*Q(4,1)*tt3
jac(1,6) = (-3.141592653589793E+0*E(1,1)*Q(1,1)*tt1*tt5)+1.570796&
&3267948966E+0*E(1,1)*tt1*Q(4,1)*tt4+1.5707963267948966E+0*E(1,1)*&
&tt1*Q(3,1)*tt3
jac(1,7) = 3.141592653589793E+0*E(1,1)*tt1*Q(4,1)*tt5+1.570796326&
&7948966E+0*E(1,1)*Q(1,1)*tt1*tt4-1.5707963267948966E+0*E(1,1)*tt1&
&*Q(2,1)*tt3
jac(1,8) = (-3.141592653589793E+0*E(1,1)*tt1*Q(3,1)*tt5)-1.570796&
&3267948966E+0*E(1,1)*tt1*Q(2,1)*tt4-1.5707963267948966E+0*E(1,1)*&
&Q(1,1)*tt1*tt3
END 
SUBROUTINE rod_bend_hes(hes, Q, u, d, E, G, r) 
IMPLICIT NONE 
REAL(KIND=8) hes(8, 8) 
REAL(KIND=8) Q(4, 2) 
REAL(KIND=8) u(3, 1) 
REAL(KIND=8) d(1, 1) 
REAL(KIND=8) E(1, 1) 
REAL(KIND=8) G(1, 1) 
REAL(KIND=8) r(1, 1) 
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
tt1 = 1/d(1,1)
tt2 = r(1,1)**2
tt3 = Q(2,2)**2
tt4 = Q(3,2)**2
tt5 = 3.141592653589793E+0*E(1,1)*tt1*tt2*tt4
tt6 = Q(4,2)**2
tt7 = 3.141592653589793E+0*E(1,1)*tt1*tt2*tt6
tt8 = -6.283185307179586E+0*E(1,1)*tt1*tt2*Q(1,2)*Q(2,2)
tt9 = 3.141592653589793E+0*E(1,1)*tt1*tt2*Q(2,2)*Q(4,2)-3.1415926&
&53589793E+0*E(1,1)*tt1*tt2*Q(1,2)*Q(3,2)
tt10 = (-3.141592653589793E+0*E(1,1)*tt1*tt2*Q(1,2)*Q(4,2))-3.141&
&592653589793E+0*E(1,1)*tt1*tt2*Q(2,2)*Q(3,2)
tt11 = -3.141592653589793E+0*E(1,1)*tt1*tt2*Q(3,1)*Q(3,2)
tt12 = -3.141592653589793E+0*E(1,1)*tt1*tt2*Q(4,1)*Q(4,2)
tt13 = tt12+tt11-6.283185307179586E+0*E(1,1)*tt1*tt2*Q(2,1)*Q(2,2&
&)
tt14 = (-tt1*(2*Q(3,1)*Q(4,2)-2*Q(3,2)*Q(4,1)+2*Q(1,1)*Q(2,2)-2*Q&
&(1,2)*Q(2,1)))-u(3,1)
tt15 = -3.141592653589793E+0*E(1,1)*tt2*tt14
tt16 = tt15-3.141592653589793E+0*E(1,1)*tt1*tt2*Q(3,1)*Q(4,2)+3.1&
&41592653589793E+0*E(1,1)*tt1*tt2*Q(3,2)*Q(4,1)+6.283185307179586E&
&+0*E(1,1)*Q(1,1)*tt1*tt2*Q(2,2)
tt17 = 3.141592653589793E+0*E(1,1)*Q(1,1)*tt1*tt2*Q(3,2)
tt18 = 3.141592653589793E+0*E(1,1)*tt1*tt2*Q(2,1)*Q(4,2)
tt19 = (-tt1*(2*Q(2,1)*Q(4,2)-2*Q(2,2)*Q(4,1)-2*Q(1,1)*Q(3,2)+2*Q&
&(1,2)*Q(3,1)))-u(2,1)
tt20 = 1.5707963267948966E+0*E(1,1)*tt2*tt19
tt21 = tt20+tt18-6.283185307179586E+0*E(1,1)*tt1*tt2*Q(2,2)*Q(4,1&
&)+tt17
tt22 = (-tt1*(2*Q(1,1)*Q(4,2)-2*Q(1,2)*Q(4,1)+2*Q(2,1)*Q(3,2)-2*Q&
&(2,2)*Q(3,1)))-u(1,1)
tt23 = -1.5707963267948966E+0*E(1,1)*tt2*tt22
tt24 = tt23+3.141592653589793E+0*E(1,1)*Q(1,1)*tt1*tt2*Q(4,2)-3.1&
&41592653589793E+0*E(1,1)*tt1*tt2*Q(2,1)*Q(3,2)+6.283185307179586E&
&+0*E(1,1)*tt1*tt2*Q(2,2)*Q(3,1)
tt25 = Q(1,2)**2
tt26 = 3.141592653589793E+0*E(1,1)*tt1*tt2*Q(1,2)*Q(3,2)-3.141592&
&653589793E+0*E(1,1)*tt1*tt2*Q(2,2)*Q(4,2)
tt27 = 3.141592653589793E+0*E(1,1)*tt2*tt14
tt28 = tt27+3.141592653589793E+0*E(1,1)*tt1*tt2*Q(3,1)*Q(4,2)-3.1&
&41592653589793E+0*E(1,1)*tt1*tt2*Q(3,2)*Q(4,1)+6.283185307179586E&
&+0*E(1,1)*tt1*tt2*Q(1,2)*Q(2,1)
tt29 = tt12+tt11-6.283185307179586E+0*E(1,1)*Q(1,1)*tt1*tt2*Q(1,2&
&)
tt30 = tt23-3.141592653589793E+0*E(1,1)*Q(1,1)*tt1*tt2*Q(4,2)+6.2&
&83185307179586E+0*E(1,1)*tt1*tt2*Q(1,2)*Q(4,1)+3.141592653589793E&
&+0*E(1,1)*tt1*tt2*Q(2,1)*Q(3,2)
tt31 = -1.5707963267948966E+0*E(1,1)*tt2*tt19
tt32 = tt31+tt18+tt17-6.283185307179586E+0*E(1,1)*tt1*tt2*Q(1,2)*&
&Q(3,1)
tt33 = 3.141592653589793E+0*E(1,1)*tt1*tt2*tt25
tt34 = 3.141592653589793E+0*E(1,1)*tt1*tt2*tt3
tt35 = -6.283185307179586E+0*E(1,1)*tt1*tt2*Q(3,2)*Q(4,2)
tt36 = 3.141592653589793E+0*E(1,1)*tt1*tt2*Q(1,2)*Q(3,1)
tt37 = 3.141592653589793E+0*E(1,1)*tt1*tt2*Q(2,2)*Q(4,1)
tt38 = tt31-6.283185307179586E+0*E(1,1)*tt1*tt2*Q(2,1)*Q(4,2)+tt3&
&7+tt36
tt39 = 1.5707963267948966E+0*E(1,1)*tt2*tt22
tt40 = tt39+6.283185307179586E+0*E(1,1)*Q(1,1)*tt1*tt2*Q(4,2)-3.1&
&41592653589793E+0*E(1,1)*tt1*tt2*Q(1,2)*Q(4,1)+3.141592653589793E&
&+0*E(1,1)*tt1*tt2*Q(2,2)*Q(3,1)
tt41 = -3.141592653589793E+0*E(1,1)*Q(1,1)*tt1*tt2*Q(1,2)
tt42 = -3.141592653589793E+0*E(1,1)*tt1*tt2*Q(2,1)*Q(2,2)
tt43 = (-6.283185307179586E+0*E(1,1)*tt1*tt2*Q(4,1)*Q(4,2))+tt42+&
&tt41
tt44 = tt15+6.283185307179586E+0*E(1,1)*tt1*tt2*Q(3,1)*Q(4,2)-3.1&
&41592653589793E+0*E(1,1)*Q(1,1)*tt1*tt2*Q(2,2)+3.141592653589793E&
&+0*E(1,1)*tt1*tt2*Q(1,2)*Q(2,1)
tt45 = tt39+3.141592653589793E+0*E(1,1)*tt1*tt2*Q(1,2)*Q(4,1)+6.2&
&83185307179586E+0*E(1,1)*tt1*tt2*Q(2,1)*Q(3,2)-3.141592653589793E&
&+0*E(1,1)*tt1*tt2*Q(2,2)*Q(3,1)
tt46 = tt20+tt37-6.283185307179586E+0*E(1,1)*Q(1,1)*tt1*tt2*Q(3,2&
&)+tt36
tt47 = tt27+6.283185307179586E+0*E(1,1)*tt1*tt2*Q(3,2)*Q(4,1)+3.1&
&41592653589793E+0*E(1,1)*Q(1,1)*tt1*tt2*Q(2,2)-3.141592653589793E&
&+0*E(1,1)*tt1*tt2*Q(1,2)*Q(2,1)
tt48 = (-6.283185307179586E+0*E(1,1)*tt1*tt2*Q(3,1)*Q(3,2))+tt42+&
&tt41
tt49 = Q(2,1)**2
tt50 = Q(3,1)**2
tt51 = 3.141592653589793E+0*E(1,1)*tt1*tt2*tt50
tt52 = Q(4,1)**2
tt53 = 3.141592653589793E+0*E(1,1)*tt1*tt2*tt52
tt54 = -6.283185307179586E+0*E(1,1)*Q(1,1)*tt1*tt2*Q(2,1)
tt55 = 3.141592653589793E+0*E(1,1)*tt1*tt2*Q(2,1)*Q(4,1)-3.141592&
&653589793E+0*E(1,1)*Q(1,1)*tt1*tt2*Q(3,1)
tt56 = (-3.141592653589793E+0*E(1,1)*Q(1,1)*tt1*tt2*Q(4,1))-3.141&
&592653589793E+0*E(1,1)*tt1*tt2*Q(2,1)*Q(3,1)
tt57 = Q(1,1)**2
tt58 = 3.141592653589793E+0*E(1,1)*Q(1,1)*tt1*tt2*Q(3,1)-3.141592&
&653589793E+0*E(1,1)*tt1*tt2*Q(2,1)*Q(4,1)
tt59 = 3.141592653589793E+0*E(1,1)*tt57*tt1*tt2
tt60 = 3.141592653589793E+0*E(1,1)*tt1*tt2*tt49
tt61 = -6.283185307179586E+0*E(1,1)*tt1*tt2*Q(3,1)*Q(4,1)
hes(1,1) = tt7+tt5+6.283185307179586E+0*E(1,1)*tt1*tt2*tt3
hes(1,2) = tt8
hes(1,3) = tt9
hes(1,4) = tt10
hes(1,5) = tt13
hes(1,6) = tt16
hes(1,7) = tt21
hes(1,8) = tt24
hes(2,1) = tt8
hes(2,2) = tt7+tt5+6.283185307179586E+0*E(1,1)*tt1*tt2*tt25
hes(2,3) = tt10
hes(2,4) = tt26
hes(2,5) = tt28
hes(2,6) = tt29
hes(2,7) = tt30
hes(2,8) = tt32
hes(3,1) = tt9
hes(3,2) = tt10
hes(3,3) = 6.283185307179586E+0*E(1,1)*tt1*tt2*tt6+tt34+tt33
hes(3,4) = tt35
hes(3,5) = tt38
hes(3,6) = tt40
hes(3,7) = tt43
hes(3,8) = tt44
hes(4,1) = tt10
hes(4,2) = tt26
hes(4,3) = tt35
hes(4,4) = 6.283185307179586E+0*E(1,1)*tt1*tt2*tt4+tt34+tt33
hes(4,5) = tt45
hes(4,6) = tt46
hes(4,7) = tt47
hes(4,8) = tt48
hes(5,1) = tt13
hes(5,2) = tt28
hes(5,3) = tt38
hes(5,4) = tt45
hes(5,5) = tt53+tt51+6.283185307179586E+0*E(1,1)*tt1*tt2*tt49
hes(5,6) = tt54
hes(5,7) = tt55
hes(5,8) = tt56
hes(6,1) = tt16
hes(6,2) = tt29
hes(6,3) = tt40
hes(6,4) = tt46
hes(6,5) = tt54
hes(6,6) = tt53+tt51+6.283185307179586E+0*E(1,1)*tt57*tt1*tt2
hes(6,7) = tt56
hes(6,8) = tt58
hes(7,1) = tt21
hes(7,2) = tt30
hes(7,3) = tt43
hes(7,4) = tt47
hes(7,5) = tt55
hes(7,6) = tt56
hes(7,7) = 6.283185307179586E+0*E(1,1)*tt1*tt2*tt52+tt60+tt59
hes(7,8) = tt61
hes(8,1) = tt24
hes(8,2) = tt32
hes(8,3) = tt44
hes(8,4) = tt48
hes(8,5) = tt56
hes(8,6) = tt58
hes(8,7) = tt61
hes(8,8) = 6.283185307179586E+0*E(1,1)*tt1*tt2*tt50+tt60+tt59
END 
SUBROUTINE rod_couple(val, XQ, d, kappa) 
IMPLICIT NONE 
REAL(KIND=8) val(1, 1) 
REAL(KIND=8) XQ(10, 1) 
REAL(KIND=8) d(1, 1) 
REAL(KIND=8) kappa(1, 1) 
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
tt1 = XQ(1,1)**2
tt2 = XQ(2,1)**2
tt3 = XQ(3,1)**2
tt4 = XQ(4,1)**2
tt5 = XQ(5,1)**2
tt6 = XQ(6,1)**2
tt7 = tt6-2*XQ(3,1)*XQ(6,1)+tt5-2*XQ(2,1)*XQ(5,1)+tt4-2*XQ(1,1)*X&
&Q(4,1)+tt3+tt2+tt1
tt8 = -tt1
tt9 = -tt2
tt10 = -tt3
tt11 = 2*XQ(1,1)*XQ(4,1)
tt12 = -tt4
tt13 = 2*XQ(2,1)*XQ(5,1)
tt14 = -tt5
tt15 = 2*XQ(3,1)*XQ(6,1)
tt16 = -tt6
tt17 = tt16+tt15+tt14+tt13+tt12+tt11+tt10+tt9+tt8
tt18 = (-2*tt6)+4*XQ(3,1)*XQ(6,1)-2*tt5+4*XQ(2,1)*XQ(5,1)-2*tt4+4&
&*XQ(1,1)*XQ(4,1)-2*tt3-2*tt2-2*tt1
tt19 = XQ(7,1)**2
tt20 = XQ(8,1)**2
tt21 = tt18*tt19
tt22 = tt18*tt20
tt23 = XQ(9,1)**2
tt24 = XQ(10,1)**2
tt25 = 2*XQ(3,1)-2*XQ(6,1)
tt26 = 4*XQ(4,1)-4*XQ(1,1)
tt27 = 2*XQ(6,1)-2*XQ(3,1)
val(1,1) = -(5.0E-1*d(1,1)*kappa(1,1)*(tt17*XQ(10,1)**4+sqrt(tt7)&
&*(tt27*tt24+(tt26*XQ(8,1)+(4*XQ(2,1)-4*XQ(5,1))*XQ(7,1))*XQ(10,1)&
&+tt27*tt23+((4*XQ(5,1)-4*XQ(2,1))*XQ(8,1)+tt26*XQ(7,1))*XQ(9,1)+t&
&t25*tt20+tt25*tt19)+(tt18*tt23+tt22+tt21)*tt24+tt17*XQ(9,1)**4+(t&
&t22+tt21)*tt23+tt17*XQ(8,1)**4+tt18*tt19*tt20+tt17*XQ(7,1)**4+tt1&
&6+tt15+tt14+tt13+tt12+tt11+tt10+tt9+tt8))/tt7
END 
SUBROUTINE rod_couple_jac(jac, XQ, d, kappa) 
IMPLICIT NONE 
REAL(KIND=8) jac(1, 10) 
REAL(KIND=8) XQ(10, 1) 
REAL(KIND=8) d(1, 1) 
REAL(KIND=8) kappa(1, 1) 
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
tt1 = XQ(1,1)**2
tt2 = XQ(2,1)**2
tt3 = XQ(3,1)**2
tt4 = XQ(4,1)**2
tt5 = XQ(5,1)**2
tt6 = XQ(6,1)**2
tt7 = tt6-2*XQ(3,1)*XQ(6,1)+tt5-2*XQ(2,1)*XQ(5,1)+tt4-2*XQ(1,1)*X&
&Q(4,1)+tt3+tt2+tt1
tt8 = 1/tt7
tt9 = -2*XQ(1,1)
tt10 = 2*XQ(4,1)
tt11 = tt10+tt9
tt12 = XQ(7,1)**4
tt13 = 4*XQ(4,1)-4*XQ(1,1)
tt14 = XQ(7,1)**2
tt15 = XQ(8,1)**2
tt16 = XQ(8,1)**4
tt17 = tt13*tt14
tt18 = tt13*tt15
tt19 = XQ(9,1)**2
tt20 = XQ(9,1)**4
tt21 = sqrt(tt7)
tt22 = XQ(10,1)**2
tt23 = 2*XQ(1,1)
tt24 = -2*XQ(4,1)
tt25 = tt24+tt23
tt26 = 1/tt21
tt27 = 2*XQ(3,1)
tt28 = -2*XQ(6,1)
tt29 = tt28+tt27
tt30 = tt13*XQ(7,1)
tt31 = 4*XQ(5,1)-4*XQ(2,1)
tt32 = tt31*XQ(8,1)
tt33 = -2*XQ(3,1)
tt34 = 2*XQ(6,1)
tt35 = tt34+tt33
tt36 = 4*XQ(2,1)-4*XQ(5,1)
tt37 = tt36*XQ(7,1)
tt38 = tt13*XQ(8,1)
tt39 = tt35*tt22+(tt38+tt37)*XQ(10,1)+tt35*tt19+(tt32+tt30)*XQ(9,&
&1)+tt29*tt15+tt29*tt14
tt40 = XQ(10,1)**4
tt41 = 1/tt7**2
tt42 = -tt1
tt43 = -tt2
tt44 = -tt3
tt45 = 2*XQ(1,1)*XQ(4,1)
tt46 = -tt4
tt47 = 2*XQ(2,1)*XQ(5,1)
tt48 = -tt5
tt49 = 2*XQ(3,1)*XQ(6,1)
tt50 = -tt6
tt51 = tt50+tt49+tt48+tt47+tt46+tt45+tt44+tt43+tt42
tt52 = (-2*tt6)+4*XQ(3,1)*XQ(6,1)-2*tt5+4*XQ(2,1)*XQ(5,1)-2*tt4+4&
&*XQ(1,1)*XQ(4,1)-2*tt3-2*tt2-2*tt1
tt53 = tt52*tt14
tt54 = tt52*tt15
tt55 = tt54+tt53
tt56 = tt52*tt19+tt54+tt53
tt57 = tt51*tt40+tt21*tt39+tt56*tt22+tt51*tt20+tt55*tt19+tt51*tt1&
&6+tt52*tt14*tt15+tt51*tt12+tt50+tt49+tt48+tt47+tt46+tt45+tt44+tt4&
&3+tt42
tt58 = -2*XQ(2,1)
tt59 = 2*XQ(5,1)
tt60 = tt59+tt58
tt61 = tt31*tt14
tt62 = tt31*tt15
tt63 = 2*XQ(2,1)
tt64 = -2*XQ(5,1)
tt65 = tt64+tt63
tt66 = 4*XQ(6,1)-4*XQ(3,1)
tt67 = tt66*tt14
tt68 = tt66*tt15
tt69 = 4*XQ(1,1)-4*XQ(4,1)
tt70 = tt69*tt14
tt71 = tt69*tt15
tt72 = tt36*tt14
tt73 = tt36*tt15
tt74 = 4*XQ(3,1)-4*XQ(6,1)
tt75 = tt74*tt14
tt76 = tt74*tt15
jac(1,1) = 5.0E-1*d(1,1)*kappa(1,1)*tt25*tt41*tt57-5.0E-1*d(1,1)*&
&kappa(1,1)*tt8*(tt11*tt40+(tt25*tt26*tt39)/2.0E+0+(tt13*tt19+tt18&
&+tt17)*tt22+tt21*((-4*XQ(8,1)*XQ(10,1))-4*XQ(7,1)*XQ(9,1))+tt11*t&
&t20+(tt18+tt17)*tt19+tt11*tt16+tt13*tt14*tt15+tt11*tt12+tt10+tt9)
jac(1,2) = 5.0E-1*d(1,1)*kappa(1,1)*tt65*tt41*tt57-5.0E-1*d(1,1)*&
&kappa(1,1)*tt8*(tt60*tt40+(tt65*tt26*tt39)/2.0E+0+(tt31*tt19+tt62&
&+tt61)*tt22+tt21*(4*XQ(7,1)*XQ(10,1)-4*XQ(8,1)*XQ(9,1))+tt60*tt20&
&+(tt62+tt61)*tt19+tt60*tt16+tt31*tt14*tt15+tt60*tt12+tt59+tt58)
jac(1,3) = 5.0E-1*d(1,1)*kappa(1,1)*tt29*tt41*tt57-5.0E-1*d(1,1)*&
&kappa(1,1)*tt8*(tt35*tt40+(tt29*tt26*tt39)/2.0E+0+(tt66*tt19+tt68&
&+tt67)*tt22+tt21*((-2*tt22)-2*tt19+2*tt15+2*tt14)+tt35*tt20+(tt68&
&+tt67)*tt19+tt35*tt16+tt66*tt14*tt15+tt35*tt12+tt34+tt33)
jac(1,4) = 5.0E-1*d(1,1)*kappa(1,1)*tt11*tt41*tt57-5.0E-1*d(1,1)*&
&kappa(1,1)*tt8*(tt25*tt40+(tt11*tt26*tt39)/2.0E+0+(tt69*tt19+tt71&
&+tt70)*tt22+tt21*(4*XQ(8,1)*XQ(10,1)+4*XQ(7,1)*XQ(9,1))+tt25*tt20&
&+(tt71+tt70)*tt19+tt25*tt16+tt69*tt14*tt15+tt25*tt12+tt24+tt23)
jac(1,5) = 5.0E-1*d(1,1)*kappa(1,1)*tt60*tt41*tt57-5.0E-1*d(1,1)*&
&kappa(1,1)*tt8*(tt65*tt40+(tt60*tt26*tt39)/2.0E+0+(tt36*tt19+tt73&
&+tt72)*tt22+tt21*(4*XQ(8,1)*XQ(9,1)-4*XQ(7,1)*XQ(10,1))+tt65*tt20&
&+(tt73+tt72)*tt19+tt65*tt16+tt36*tt14*tt15+tt65*tt12+tt64+tt63)
jac(1,6) = 5.0E-1*d(1,1)*kappa(1,1)*tt35*tt41*tt57-5.0E-1*d(1,1)*&
&kappa(1,1)*tt8*(tt29*tt40+(tt35*tt26*tt39)/2.0E+0+tt21*(2*tt22+2*&
&tt19-2*tt15-2*tt14)+(tt74*tt19+tt76+tt75)*tt22+tt29*tt20+(tt76+tt&
&75)*tt19+tt29*tt16+tt74*tt14*tt15+tt29*tt12+tt28+tt27)
jac(1,7) = -5.0E-1*d(1,1)*kappa(1,1)*tt8*(2*tt52*XQ(7,1)*tt22+tt2&
&1*(tt36*XQ(10,1)+tt13*XQ(9,1)+2*tt29*XQ(7,1))+2*tt52*XQ(7,1)*tt19&
&+2*tt52*XQ(7,1)*tt15+4*tt51*XQ(7,1)**3)
jac(1,8) = -5.0E-1*d(1,1)*kappa(1,1)*tt8*(2*tt52*XQ(8,1)*tt22+tt2&
&1*(tt13*XQ(10,1)+tt31*XQ(9,1)+2*tt29*XQ(8,1))+2*tt52*XQ(8,1)*tt19&
&+4*tt51*XQ(8,1)**3+2*tt52*tt14*XQ(8,1))
jac(1,9) = -5.0E-1*d(1,1)*kappa(1,1)*tt8*(2*tt52*XQ(9,1)*tt22+4*t&
&t51*XQ(9,1)**3+tt21*(2*tt35*XQ(9,1)+tt32+tt30)+2*tt55*XQ(9,1))
jac(1,10) = -5.0E-1*d(1,1)*kappa(1,1)*tt8*(4*tt51*XQ(10,1)**3+tt2&
&1*(2*tt35*XQ(10,1)+tt38+tt37)+2*tt56*XQ(10,1))
END 
SUBROUTINE rod_couple_hes(hes, XQ, d, kappa) 
IMPLICIT NONE 
REAL(KIND=8) hes(10, 10) 
REAL(KIND=8) XQ(10, 1) 
REAL(KIND=8) d(1, 1) 
REAL(KIND=8) kappa(1, 1) 
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
tt1 = XQ(1,1)**2
tt2 = XQ(2,1)**2
tt3 = XQ(3,1)**2
tt4 = XQ(4,1)**2
tt5 = XQ(5,1)**2
tt6 = XQ(6,1)**2
tt7 = tt6-2*XQ(3,1)*XQ(6,1)+tt5-2*XQ(2,1)*XQ(5,1)+tt4-2*XQ(1,1)*X&
&Q(4,1)+tt3+tt2+tt1
tt8 = 1/tt7
tt9 = XQ(7,1)**4
tt10 = -2*tt9
tt11 = XQ(7,1)**2
tt12 = XQ(8,1)**2
tt13 = -4*tt11*tt12
tt14 = XQ(8,1)**4
tt15 = -2*tt14
tt16 = -4*tt11
tt17 = -4*tt12
tt18 = XQ(9,1)**2
tt19 = (tt17+tt16)*tt18
tt20 = XQ(9,1)**4
tt21 = -2*tt20
tt22 = 2*XQ(1,1)
tt23 = -2*XQ(4,1)
tt24 = tt23+tt22
tt25 = sqrt(tt7)
tt26 = 1/tt25
tt27 = (-4*XQ(8,1)*XQ(10,1))-4*XQ(7,1)*XQ(9,1)
tt28 = XQ(10,1)**2
tt29 = ((-4*tt18)+tt17+tt16)*tt28
tt30 = tt24**2
tt31 = 1/tt25**3
tt32 = 2*XQ(3,1)
tt33 = -2*XQ(6,1)
tt34 = tt33+tt32
tt35 = 4*XQ(4,1)-4*XQ(1,1)
tt36 = tt35*XQ(7,1)
tt37 = 4*XQ(5,1)-4*XQ(2,1)
tt38 = tt37*XQ(8,1)
tt39 = -2*XQ(3,1)
tt40 = 2*XQ(6,1)
tt41 = tt40+tt39
tt42 = 4*XQ(2,1)-4*XQ(5,1)
tt43 = tt42*XQ(7,1)
tt44 = tt35*XQ(8,1)
tt45 = tt41*tt28+(tt44+tt43)*XQ(10,1)+tt41*tt18+(tt38+tt36)*XQ(9,&
&1)+tt34*tt12+tt34*tt11
tt46 = tt26*tt45
tt47 = XQ(10,1)**4
tt48 = -2*tt47
tt49 = 1/tt7**2
tt50 = -2*XQ(1,1)
tt51 = 2*XQ(4,1)
tt52 = tt51+tt50
tt53 = tt35*tt11
tt54 = tt35*tt12
tt55 = tt54+tt53
tt56 = tt35*tt18+tt54+tt53
tt57 = tt52*tt47+(tt24*tt26*tt45)/2.0E+0+tt56*tt28+tt25*tt27+tt52&
&*tt20+tt55*tt18+tt52*tt14+tt35*tt11*tt12+tt52*tt9+tt51+tt50
tt58 = 1/tt7**3
tt59 = -tt1
tt60 = -tt2
tt61 = -tt3
tt62 = 2*XQ(1,1)*XQ(4,1)
tt63 = -tt4
tt64 = 2*XQ(2,1)*XQ(5,1)
tt65 = -tt5
tt66 = 2*XQ(3,1)*XQ(6,1)
tt67 = -tt6
tt68 = tt67+tt66+tt65+tt64+tt63+tt62+tt61+tt60+tt59
tt69 = (-2*tt6)+4*XQ(3,1)*XQ(6,1)-2*tt5+4*XQ(2,1)*XQ(5,1)-2*tt4+4&
&*XQ(1,1)*XQ(4,1)-2*tt3-2*tt2-2*tt1
tt70 = tt69*tt11
tt71 = tt69*tt12
tt72 = tt71+tt70
tt73 = tt69*tt18+tt71+tt70
tt74 = tt68*tt47+tt25*tt45+tt73*tt28+tt68*tt20+tt72*tt18+tt68*tt1&
&4+tt69*tt11*tt12+tt68*tt9+tt67+tt66+tt65+tt64+tt63+tt62+tt61+tt60&
&+tt59
tt75 = 1.0E+0*d(1,1)*kappa(1,1)*tt49*tt74
tt76 = 4*XQ(7,1)*XQ(10,1)-4*XQ(8,1)*XQ(9,1)
tt77 = 2*XQ(2,1)
tt78 = -2*XQ(5,1)
tt79 = tt78+tt77
tt80 = -2*XQ(2,1)
tt81 = 2*XQ(5,1)
tt82 = tt81+tt80
tt83 = tt37*tt11
tt84 = tt37*tt12
tt85 = tt84+tt83
tt86 = tt37*tt18+tt84+tt83
tt87 = tt82*tt47+(tt79*tt26*tt45)/2.0E+0+tt86*tt28+tt25*tt76+tt82&
&*tt20+tt85*tt18+tt82*tt14+tt37*tt11*tt12+tt82*tt9+tt81+tt80
tt88 = (-1.0E+0*d(1,1)*kappa(1,1)*tt24*tt79*tt58*tt74)+5.0E-1*d(1&
&,1)*kappa(1,1)*tt24*tt49*tt87+5.0E-1*d(1,1)*kappa(1,1)*tt79*tt49*&
&tt57-5.0E-1*d(1,1)*kappa(1,1)*tt8*((-(tt24*tt79*tt31*tt45)/4.0E+0&
&)+(tt79*tt26*tt27)/2.0E+0+(tt24*tt26*tt76)/2.0E+0)
tt89 = (-2*tt28)-2*tt18+2*tt12+2*tt11
tt90 = 4*XQ(6,1)-4*XQ(3,1)
tt91 = tt90*tt11
tt92 = tt90*tt12
tt93 = tt92+tt91
tt94 = tt90*tt18+tt92+tt91
tt95 = tt41*tt47+(tt34*tt26*tt45)/2.0E+0+tt94*tt28+tt25*tt89+tt41&
&*tt20+tt93*tt18+tt41*tt14+tt90*tt11*tt12+tt41*tt9+tt40+tt39
tt96 = (-1.0E+0*d(1,1)*kappa(1,1)*tt24*tt34*tt58*tt74)+5.0E-1*d(1&
&,1)*kappa(1,1)*tt24*tt49*tt95+5.0E-1*d(1,1)*kappa(1,1)*tt34*tt49*&
&tt57-5.0E-1*d(1,1)*kappa(1,1)*tt8*((-(tt24*tt34*tt31*tt45)/4.0E+0&
&)+(tt24*tt26*tt89)/2.0E+0+(tt34*tt26*tt27)/2.0E+0)
tt97 = 2*tt9
tt98 = 4*tt11*tt12
tt99 = 2*tt14
tt100 = 4*tt11
tt101 = 4*tt12
tt102 = (tt101+tt100)*tt18
tt103 = 2*tt20
tt104 = 4*XQ(8,1)*XQ(10,1)+4*XQ(7,1)*XQ(9,1)
tt105 = (4*tt18+tt101+tt100)*tt28
tt106 = -tt26*tt45
tt107 = 2*tt47
tt108 = 4*XQ(1,1)-4*XQ(4,1)
tt109 = tt108*tt11
tt110 = tt108*tt12
tt111 = tt110+tt109
tt112 = tt108*tt18+tt110+tt109
tt113 = tt24*tt47+(tt52*tt26*tt45)/2.0E+0+tt112*tt28+tt25*tt104+t&
&t24*tt20+tt111*tt18+tt24*tt14+tt108*tt11*tt12+tt24*tt9+tt23+tt22
tt114 = -1.0E+0*d(1,1)*kappa(1,1)*tt49*tt74
tt115 = tt114-1.0E+0*d(1,1)*kappa(1,1)*tt24*tt52*tt58*tt74+5.0E-1&
&*d(1,1)*kappa(1,1)*tt52*tt49*tt57+5.0E-1*d(1,1)*kappa(1,1)*tt24*t&
&t49*tt113-5.0E-1*d(1,1)*kappa(1,1)*tt8*(tt107+tt106-(tt24*tt52*tt&
&31*tt45)/4.0E+0+tt105+(tt24*tt26*tt104)/2.0E+0+(tt52*tt26*tt27)/2&
&.0E+0+tt103+tt102+tt99+tt98+tt97+2)
tt116 = 4*XQ(8,1)*XQ(9,1)-4*XQ(7,1)*XQ(10,1)
tt117 = tt42*tt11
tt118 = tt42*tt12
tt119 = tt118+tt117
tt120 = tt42*tt18+tt118+tt117
tt121 = tt79*tt47+(tt82*tt26*tt45)/2.0E+0+tt120*tt28+tt25*tt116+t&
&t79*tt20+tt119*tt18+tt79*tt14+tt42*tt11*tt12+tt79*tt9+tt78+tt77
tt122 = (-1.0E+0*d(1,1)*kappa(1,1)*tt24*tt82*tt58*tt74)+5.0E-1*d(&
&1,1)*kappa(1,1)*tt24*tt49*tt121+5.0E-1*d(1,1)*kappa(1,1)*tt82*tt4&
&9*tt57-5.0E-1*d(1,1)*kappa(1,1)*tt8*((-(tt24*tt82*tt31*tt45)/4.0E&
&+0)+(tt82*tt26*tt27)/2.0E+0+(tt24*tt26*tt116)/2.0E+0)
tt123 = 2*tt28+2*tt18-2*tt12-2*tt11
tt124 = 4*XQ(3,1)-4*XQ(6,1)
tt125 = tt124*tt11
tt126 = tt124*tt12
tt127 = tt126+tt125
tt128 = tt124*tt18+tt126+tt125
tt129 = tt34*tt47+(tt41*tt26*tt45)/2.0E+0+tt25*tt123+tt128*tt28+t&
&t34*tt20+tt127*tt18+tt34*tt14+tt124*tt11*tt12+tt34*tt9+tt33+tt32
tt130 = (-1.0E+0*d(1,1)*kappa(1,1)*tt24*tt41*tt58*tt74)+5.0E-1*d(&
&1,1)*kappa(1,1)*tt24*tt49*tt129+5.0E-1*d(1,1)*kappa(1,1)*tt41*tt4&
&9*tt57-5.0E-1*d(1,1)*kappa(1,1)*tt8*((-(tt24*tt41*tt31*tt45)/4.0E&
&+0)+(tt24*tt26*tt123)/2.0E+0+(tt41*tt26*tt27)/2.0E+0)
tt131 = XQ(7,1)**3
tt132 = -4*tt25*XQ(9,1)
tt133 = tt42*XQ(10,1)+tt35*XQ(9,1)+2*tt34*XQ(7,1)
tt134 = 2*tt69*XQ(7,1)*tt28+tt25*tt133+2*tt69*XQ(7,1)*tt18+2*tt69&
&*XQ(7,1)*tt12+4*tt68*tt131
tt135 = 5.0E-1*d(1,1)*kappa(1,1)*tt24*tt49*tt134-5.0E-1*d(1,1)*ka&
&ppa(1,1)*tt8*(2*tt35*XQ(7,1)*tt28+(tt24*tt26*tt133)/2.0E+0+2*tt35&
&*XQ(7,1)*tt18+tt132+2*tt35*XQ(7,1)*tt12+4*tt52*tt131)
tt136 = XQ(8,1)**3
tt137 = -4*tt25*XQ(10,1)
tt138 = tt35*XQ(10,1)+tt37*XQ(9,1)+2*tt34*XQ(8,1)
tt139 = 2*tt69*XQ(8,1)*tt28+tt25*tt138+2*tt69*XQ(8,1)*tt18+4*tt68&
&*tt136+2*tt69*tt11*XQ(8,1)
tt140 = 5.0E-1*d(1,1)*kappa(1,1)*tt24*tt49*tt139-5.0E-1*d(1,1)*ka&
&ppa(1,1)*tt8*(2*tt35*XQ(8,1)*tt28+(tt24*tt26*tt138)/2.0E+0+tt137+&
&2*tt35*XQ(8,1)*tt18+4*tt52*tt136+2*tt35*tt11*XQ(8,1))
tt141 = -4*tt25*XQ(7,1)
tt142 = 2*tt41*XQ(9,1)+tt38+tt36
tt143 = XQ(9,1)**3
tt144 = 2*tt69*XQ(9,1)*tt28+4*tt68*tt143+tt25*tt142+2*tt72*XQ(9,1&
&)
tt145 = 5.0E-1*d(1,1)*kappa(1,1)*tt24*tt49*tt144-5.0E-1*d(1,1)*ka&
&ppa(1,1)*tt8*(2*tt35*XQ(9,1)*tt28+4*tt52*tt143+(tt24*tt26*tt142)/&
&2.0E+0+2*tt55*XQ(9,1)+tt141)
tt146 = -4*tt25*XQ(8,1)
tt147 = 2*tt41*XQ(10,1)+tt44+tt43
tt148 = XQ(10,1)**3
tt149 = 4*tt68*tt148+tt25*tt147+2*tt73*XQ(10,1)
tt150 = 5.0E-1*d(1,1)*kappa(1,1)*tt24*tt49*tt149-5.0E-1*d(1,1)*ka&
&ppa(1,1)*tt8*(4*tt52*tt148+(tt24*tt26*tt147)/2.0E+0+2*tt56*XQ(10,&
&1)+tt146)
tt151 = tt79**2
tt152 = (-1.0E+0*d(1,1)*kappa(1,1)*tt79*tt34*tt58*tt74)+5.0E-1*d(&
&1,1)*kappa(1,1)*tt79*tt49*tt95+5.0E-1*d(1,1)*kappa(1,1)*tt34*tt49&
&*tt87-5.0E-1*d(1,1)*kappa(1,1)*tt8*((-(tt79*tt34*tt31*tt45)/4.0E+&
&0)+(tt79*tt26*tt89)/2.0E+0+(tt34*tt26*tt76)/2.0E+0)
tt153 = (-1.0E+0*d(1,1)*kappa(1,1)*tt52*tt79*tt58*tt74)+5.0E-1*d(&
&1,1)*kappa(1,1)*tt52*tt49*tt87+5.0E-1*d(1,1)*kappa(1,1)*tt79*tt49&
&*tt113-5.0E-1*d(1,1)*kappa(1,1)*tt8*((-(tt52*tt79*tt31*tt45)/4.0E&
&+0)+(tt79*tt26*tt104)/2.0E+0+(tt52*tt26*tt76)/2.0E+0)
tt154 = tt114-1.0E+0*d(1,1)*kappa(1,1)*tt79*tt82*tt58*tt74+5.0E-1&
&*d(1,1)*kappa(1,1)*tt82*tt49*tt87+5.0E-1*d(1,1)*kappa(1,1)*tt79*t&
&t49*tt121-5.0E-1*d(1,1)*kappa(1,1)*tt8*(tt107+tt106-(tt79*tt82*tt&
&31*tt45)/4.0E+0+tt105+(tt82*tt26*tt76)/2.0E+0+(tt79*tt26*tt116)/2&
&.0E+0+tt103+tt102+tt99+tt98+tt97+2)
tt155 = (-1.0E+0*d(1,1)*kappa(1,1)*tt79*tt41*tt58*tt74)+5.0E-1*d(&
&1,1)*kappa(1,1)*tt79*tt49*tt129+5.0E-1*d(1,1)*kappa(1,1)*tt41*tt4&
&9*tt87-5.0E-1*d(1,1)*kappa(1,1)*tt8*((-(tt79*tt41*tt31*tt45)/4.0E&
&+0)+(tt79*tt26*tt123)/2.0E+0+(tt41*tt26*tt76)/2.0E+0)
tt156 = 4*tt25*XQ(10,1)
tt157 = 5.0E-1*d(1,1)*kappa(1,1)*tt79*tt49*tt134-5.0E-1*d(1,1)*ka&
&ppa(1,1)*tt8*(2*tt37*XQ(7,1)*tt28+(tt79*tt26*tt133)/2.0E+0+tt156+&
&2*tt37*XQ(7,1)*tt18+2*tt37*XQ(7,1)*tt12+4*tt82*tt131)
tt158 = 5.0E-1*d(1,1)*kappa(1,1)*tt79*tt49*tt139-5.0E-1*d(1,1)*ka&
&ppa(1,1)*tt8*(2*tt37*XQ(8,1)*tt28+(tt79*tt26*tt138)/2.0E+0+2*tt37&
&*XQ(8,1)*tt18+tt132+4*tt82*tt136+2*tt37*tt11*XQ(8,1))
tt159 = 5.0E-1*d(1,1)*kappa(1,1)*tt79*tt49*tt144-5.0E-1*d(1,1)*ka&
&ppa(1,1)*tt8*(2*tt37*XQ(9,1)*tt28+4*tt82*tt143+(tt79*tt26*tt142)/&
&2.0E+0+2*tt85*XQ(9,1)+tt146)
tt160 = 4*tt25*XQ(7,1)
tt161 = 5.0E-1*d(1,1)*kappa(1,1)*tt79*tt49*tt149-5.0E-1*d(1,1)*ka&
&ppa(1,1)*tt8*(4*tt82*tt148+(tt79*tt26*tt147)/2.0E+0+2*tt86*XQ(10,&
&1)+tt160)
tt162 = tt34**2
tt163 = (-1.0E+0*d(1,1)*kappa(1,1)*tt52*tt34*tt58*tt74)+5.0E-1*d(&
&1,1)*kappa(1,1)*tt52*tt49*tt95+5.0E-1*d(1,1)*kappa(1,1)*tt34*tt49&
&*tt113-5.0E-1*d(1,1)*kappa(1,1)*tt8*((-(tt52*tt34*tt31*tt45)/4.0E&
&+0)+(tt52*tt26*tt89)/2.0E+0+(tt34*tt26*tt104)/2.0E+0)
tt164 = (-1.0E+0*d(1,1)*kappa(1,1)*tt82*tt34*tt58*tt74)+5.0E-1*d(&
&1,1)*kappa(1,1)*tt82*tt49*tt95+5.0E-1*d(1,1)*kappa(1,1)*tt34*tt49&
&*tt121-5.0E-1*d(1,1)*kappa(1,1)*tt8*((-(tt82*tt34*tt31*tt45)/4.0E&
&+0)+(tt82*tt26*tt89)/2.0E+0+(tt34*tt26*tt116)/2.0E+0)
tt165 = tt114-1.0E+0*d(1,1)*kappa(1,1)*tt34*tt41*tt58*tt74+5.0E-1&
&*d(1,1)*kappa(1,1)*tt41*tt49*tt95+5.0E-1*d(1,1)*kappa(1,1)*tt34*t&
&t49*tt129-5.0E-1*d(1,1)*kappa(1,1)*tt8*(tt107+tt106-(tt34*tt41*tt&
&31*tt45)/4.0E+0+(tt34*tt26*tt123)/2.0E+0+tt105+(tt41*tt26*tt89)/2&
&.0E+0+tt103+tt102+tt99+tt98+tt97+2)
tt166 = 5.0E-1*d(1,1)*kappa(1,1)*tt34*tt49*tt134-5.0E-1*d(1,1)*ka&
&ppa(1,1)*tt8*(2*tt90*XQ(7,1)*tt28+(tt34*tt26*tt133)/2.0E+0+2*tt90&
&*XQ(7,1)*tt18+2*tt90*XQ(7,1)*tt12+4*tt41*tt131+tt160)
tt167 = 4*tt25*XQ(8,1)
tt168 = 5.0E-1*d(1,1)*kappa(1,1)*tt34*tt49*tt139-5.0E-1*d(1,1)*ka&
&ppa(1,1)*tt8*(2*tt90*XQ(8,1)*tt28+(tt34*tt26*tt138)/2.0E+0+2*tt90&
&*XQ(8,1)*tt18+4*tt41*tt136+2*tt90*tt11*XQ(8,1)+tt167)
tt169 = 5.0E-1*d(1,1)*kappa(1,1)*tt34*tt49*tt144-5.0E-1*d(1,1)*ka&
&ppa(1,1)*tt8*(2*tt90*XQ(9,1)*tt28+4*tt41*tt143+(tt34*tt26*tt142)/&
&2.0E+0+2*tt93*XQ(9,1)+tt132)
tt170 = 5.0E-1*d(1,1)*kappa(1,1)*tt34*tt49*tt149-5.0E-1*d(1,1)*ka&
&ppa(1,1)*tt8*(4*tt41*tt148+(tt34*tt26*tt147)/2.0E+0+2*tt94*XQ(10,&
&1)+tt137)
tt171 = tt52**2
tt172 = (-1.0E+0*d(1,1)*kappa(1,1)*tt52*tt82*tt58*tt74)+5.0E-1*d(&
&1,1)*kappa(1,1)*tt52*tt49*tt121+5.0E-1*d(1,1)*kappa(1,1)*tt82*tt4&
&9*tt113-5.0E-1*d(1,1)*kappa(1,1)*tt8*((-(tt52*tt82*tt31*tt45)/4.0&
&E+0)+(tt82*tt26*tt104)/2.0E+0+(tt52*tt26*tt116)/2.0E+0)
tt173 = (-1.0E+0*d(1,1)*kappa(1,1)*tt52*tt41*tt58*tt74)+5.0E-1*d(&
&1,1)*kappa(1,1)*tt52*tt49*tt129+5.0E-1*d(1,1)*kappa(1,1)*tt41*tt4&
&9*tt113-5.0E-1*d(1,1)*kappa(1,1)*tt8*((-(tt52*tt41*tt31*tt45)/4.0&
&E+0)+(tt52*tt26*tt123)/2.0E+0+(tt41*tt26*tt104)/2.0E+0)
tt174 = 4*tt25*XQ(9,1)
tt175 = 5.0E-1*d(1,1)*kappa(1,1)*tt52*tt49*tt134-5.0E-1*d(1,1)*ka&
&ppa(1,1)*tt8*(2*tt108*XQ(7,1)*tt28+(tt52*tt26*tt133)/2.0E+0+2*tt1&
&08*XQ(7,1)*tt18+tt174+2*tt108*XQ(7,1)*tt12+4*tt24*tt131)
tt176 = 5.0E-1*d(1,1)*kappa(1,1)*tt52*tt49*tt139-5.0E-1*d(1,1)*ka&
&ppa(1,1)*tt8*(2*tt108*XQ(8,1)*tt28+(tt52*tt26*tt138)/2.0E+0+tt156&
&+2*tt108*XQ(8,1)*tt18+4*tt24*tt136+2*tt108*tt11*XQ(8,1))
tt177 = 5.0E-1*d(1,1)*kappa(1,1)*tt52*tt49*tt144-5.0E-1*d(1,1)*ka&
&ppa(1,1)*tt8*(2*tt108*XQ(9,1)*tt28+4*tt24*tt143+(tt52*tt26*tt142)&
&/2.0E+0+2*tt111*XQ(9,1)+tt160)
tt178 = 5.0E-1*d(1,1)*kappa(1,1)*tt52*tt49*tt149-5.0E-1*d(1,1)*ka&
&ppa(1,1)*tt8*(4*tt24*tt148+(tt52*tt26*tt147)/2.0E+0+2*tt112*XQ(10&
&,1)+tt167)
tt179 = tt82**2
tt180 = (-1.0E+0*d(1,1)*kappa(1,1)*tt82*tt41*tt58*tt74)+5.0E-1*d(&
&1,1)*kappa(1,1)*tt82*tt49*tt129+5.0E-1*d(1,1)*kappa(1,1)*tt41*tt4&
&9*tt121-5.0E-1*d(1,1)*kappa(1,1)*tt8*((-(tt82*tt41*tt31*tt45)/4.0&
&E+0)+(tt82*tt26*tt123)/2.0E+0+(tt41*tt26*tt116)/2.0E+0)
tt181 = 5.0E-1*d(1,1)*kappa(1,1)*tt82*tt49*tt134-5.0E-1*d(1,1)*ka&
&ppa(1,1)*tt8*(2*tt42*XQ(7,1)*tt28+(tt82*tt26*tt133)/2.0E+0+tt137+&
&2*tt42*XQ(7,1)*tt18+2*tt42*XQ(7,1)*tt12+4*tt79*tt131)
tt182 = 5.0E-1*d(1,1)*kappa(1,1)*tt82*tt49*tt139-5.0E-1*d(1,1)*ka&
&ppa(1,1)*tt8*(2*tt42*XQ(8,1)*tt28+(tt82*tt26*tt138)/2.0E+0+2*tt42&
&*XQ(8,1)*tt18+tt174+4*tt79*tt136+2*tt42*tt11*XQ(8,1))
tt183 = 5.0E-1*d(1,1)*kappa(1,1)*tt82*tt49*tt144-5.0E-1*d(1,1)*ka&
&ppa(1,1)*tt8*(2*tt42*XQ(9,1)*tt28+4*tt79*tt143+(tt82*tt26*tt142)/&
&2.0E+0+2*tt119*XQ(9,1)+tt167)
tt184 = 5.0E-1*d(1,1)*kappa(1,1)*tt82*tt49*tt149-5.0E-1*d(1,1)*ka&
&ppa(1,1)*tt8*(4*tt79*tt148+(tt82*tt26*tt147)/2.0E+0+2*tt120*XQ(10&
&,1)+tt141)
tt185 = tt41**2
tt186 = 5.0E-1*d(1,1)*kappa(1,1)*tt41*tt49*tt134-5.0E-1*d(1,1)*ka&
&ppa(1,1)*tt8*(2*tt124*XQ(7,1)*tt28+(tt41*tt26*tt133)/2.0E+0+2*tt1&
&24*XQ(7,1)*tt18+2*tt124*XQ(7,1)*tt12+4*tt34*tt131+tt141)
tt187 = 5.0E-1*d(1,1)*kappa(1,1)*tt41*tt49*tt139-5.0E-1*d(1,1)*ka&
&ppa(1,1)*tt8*(2*tt124*XQ(8,1)*tt28+(tt41*tt26*tt138)/2.0E+0+2*tt1&
&24*XQ(8,1)*tt18+4*tt34*tt136+2*tt124*tt11*XQ(8,1)+tt146)
tt188 = 5.0E-1*d(1,1)*kappa(1,1)*tt41*tt49*tt144-5.0E-1*d(1,1)*ka&
&ppa(1,1)*tt8*(2*tt124*XQ(9,1)*tt28+4*tt34*tt143+(tt41*tt26*tt142)&
&/2.0E+0+2*tt127*XQ(9,1)+tt174)
tt189 = 5.0E-1*d(1,1)*kappa(1,1)*tt41*tt49*tt149-5.0E-1*d(1,1)*ka&
&ppa(1,1)*tt8*(4*tt34*tt148+(tt41*tt26*tt147)/2.0E+0+2*tt128*XQ(10&
&,1)+tt156)
tt190 = 2*tt34*tt25
tt191 = 2*tt69*tt18
tt192 = 2*tt69*tt28
tt193 = -2.0E+0*d(1,1)*kappa(1,1)*tt69*tt8*XQ(7,1)*XQ(8,1)
tt194 = tt35*tt25
tt195 = -5.0E-1*d(1,1)*kappa(1,1)*tt8*(4*tt69*XQ(7,1)*XQ(9,1)+tt1&
&94)
tt196 = -5.0E-1*d(1,1)*kappa(1,1)*tt8*(4*tt69*XQ(7,1)*XQ(10,1)+tt&
&42*tt25)
tt197 = -5.0E-1*d(1,1)*kappa(1,1)*tt8*(4*tt69*XQ(8,1)*XQ(9,1)+tt3&
&7*tt25)
tt198 = -5.0E-1*d(1,1)*kappa(1,1)*tt8*(4*tt69*XQ(8,1)*XQ(10,1)+tt&
&194)
tt199 = 2*tt41*tt25
tt200 = -2.0E+0*d(1,1)*kappa(1,1)*tt69*tt8*XQ(9,1)*XQ(10,1)
hes(1,1) = tt75-1.0E+0*d(1,1)*kappa(1,1)*tt30*tt58*tt74+1.0E+0*d(&
&1,1)*kappa(1,1)*tt24*tt49*tt57-5.0E-1*d(1,1)*kappa(1,1)*tt8*(tt48&
&+tt46-(tt30*tt31*tt45)/4.0E+0+tt29+tt24*tt26*tt27+tt21+tt19+tt15+&
&tt13+tt10-2)
hes(1,2) = tt88
hes(1,3) = tt96
hes(1,4) = tt115
hes(1,5) = tt122
hes(1,6) = tt130
hes(1,7) = tt135
hes(1,8) = tt140
hes(1,9) = tt145
hes(1,10) = tt150
hes(2,1) = tt88
hes(2,2) = tt75-1.0E+0*d(1,1)*kappa(1,1)*tt151*tt58*tt74+1.0E+0*d&
&(1,1)*kappa(1,1)*tt79*tt49*tt87-5.0E-1*d(1,1)*kappa(1,1)*tt8*(tt4&
&8+tt46-(tt151*tt31*tt45)/4.0E+0+tt29+tt79*tt26*tt76+tt21+tt19+tt1&
&5+tt13+tt10-2)
hes(2,3) = tt152
hes(2,4) = tt153
hes(2,5) = tt154
hes(2,6) = tt155
hes(2,7) = tt157
hes(2,8) = tt158
hes(2,9) = tt159
hes(2,10) = tt161
hes(3,1) = tt96
hes(3,2) = tt152
hes(3,3) = tt75-1.0E+0*d(1,1)*kappa(1,1)*tt162*tt58*tt74+1.0E+0*d&
&(1,1)*kappa(1,1)*tt34*tt49*tt95-5.0E-1*d(1,1)*kappa(1,1)*tt8*(tt4&
&8+tt46-(tt162*tt31*tt45)/4.0E+0+tt29+tt34*tt26*tt89+tt21+tt19+tt1&
&5+tt13+tt10-2)
hes(3,4) = tt163
hes(3,5) = tt164
hes(3,6) = tt165
hes(3,7) = tt166
hes(3,8) = tt168
hes(3,9) = tt169
hes(3,10) = tt170
hes(4,1) = tt115
hes(4,2) = tt153
hes(4,3) = tt163
hes(4,4) = tt75-1.0E+0*d(1,1)*kappa(1,1)*tt171*tt58*tt74+1.0E+0*d&
&(1,1)*kappa(1,1)*tt52*tt49*tt113-5.0E-1*d(1,1)*kappa(1,1)*tt8*(tt&
&48+tt46-(tt171*tt31*tt45)/4.0E+0+tt29+tt52*tt26*tt104+tt21+tt19+t&
&t15+tt13+tt10-2)
hes(4,5) = tt172
hes(4,6) = tt173
hes(4,7) = tt175
hes(4,8) = tt176
hes(4,9) = tt177
hes(4,10) = tt178
hes(5,1) = tt122
hes(5,2) = tt154
hes(5,3) = tt164
hes(5,4) = tt172
hes(5,5) = tt75-1.0E+0*d(1,1)*kappa(1,1)*tt179*tt58*tt74+1.0E+0*d&
&(1,1)*kappa(1,1)*tt82*tt49*tt121-5.0E-1*d(1,1)*kappa(1,1)*tt8*(tt&
&48+tt46-(tt179*tt31*tt45)/4.0E+0+tt29+tt82*tt26*tt116+tt21+tt19+t&
&t15+tt13+tt10-2)
hes(5,6) = tt180
hes(5,7) = tt181
hes(5,8) = tt182
hes(5,9) = tt183
hes(5,10) = tt184
hes(6,1) = tt130
hes(6,2) = tt155
hes(6,3) = tt165
hes(6,4) = tt173
hes(6,5) = tt180
hes(6,6) = tt75-1.0E+0*d(1,1)*kappa(1,1)*tt185*tt58*tt74+1.0E+0*d&
&(1,1)*kappa(1,1)*tt41*tt49*tt129-5.0E-1*d(1,1)*kappa(1,1)*tt8*(tt&
&48+tt46-(tt185*tt31*tt45)/4.0E+0+tt41*tt26*tt123+tt29+tt21+tt19+t&
&t15+tt13+tt10-2)
hes(6,7) = tt186
hes(6,8) = tt187
hes(6,9) = tt188
hes(6,10) = tt189
hes(7,1) = tt135
hes(7,2) = tt157
hes(7,3) = tt166
hes(7,4) = tt175
hes(7,5) = tt181
hes(7,6) = tt186
hes(7,7) = -5.0E-1*d(1,1)*kappa(1,1)*tt8*(tt192+tt191+2*tt69*tt12&
&+12*tt68*tt11+tt190)
hes(7,8) = tt193
hes(7,9) = tt195
hes(7,10) = tt196
hes(8,1) = tt140
hes(8,2) = tt158
hes(8,3) = tt168
hes(8,4) = tt176
hes(8,5) = tt182
hes(8,6) = tt187
hes(8,7) = tt193
hes(8,8) = -5.0E-1*d(1,1)*kappa(1,1)*tt8*(tt192+tt191+12*tt68*tt1&
&2+2*tt69*tt11+tt190)
hes(8,9) = tt197
hes(8,10) = tt198
hes(9,1) = tt145
hes(9,2) = tt159
hes(9,3) = tt169
hes(9,4) = tt177
hes(9,5) = tt183
hes(9,6) = tt188
hes(9,7) = tt195
hes(9,8) = tt197
hes(9,9) = -5.0E-1*d(1,1)*kappa(1,1)*tt8*(tt192+12*tt68*tt18+2*tt&
&72+tt199)
hes(9,10) = tt200
hes(10,1) = tt150
hes(10,2) = tt161
hes(10,3) = tt170
hes(10,4) = tt178
hes(10,5) = tt184
hes(10,6) = tt189
hes(10,7) = tt196
hes(10,8) = tt198
hes(10,9) = tt200
hes(10,10) = -5.0E-1*d(1,1)*kappa(1,1)*tt8*(12*tt68*tt28+2*tt73+t&
&t199)
END 
