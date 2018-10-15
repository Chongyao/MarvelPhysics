SUBROUTINE bw98_stretch(val, X, invUV, area) 
IMPLICIT NONE 
REAL(KIND=8) val(1, 1) 
REAL(KIND=8) X(3, 3) 
REAL(KIND=8) invUV(2, 2) 
REAL(KIND=8) area(1, 1) 
REAL(KIND=8)  tt1 
REAL(KIND=8)  tt2 
REAL(KIND=8)  tt3 
REAL(KIND=8)  tt4 
REAL(KIND=8)  tt5 
REAL(KIND=8)  tt6 
REAL(KIND=8)  tt7 
REAL(KIND=8)  tt8 
REAL(KIND=8)  tt9 
tt1 = -X(1,1)
tt2 = X(1,2)+tt1
tt3 = X(1,3)+tt1
tt4 = -X(2,1)
tt5 = X(2,2)+tt4
tt6 = X(2,3)+tt4
tt7 = -X(3,1)
tt8 = X(3,2)+tt7
tt9 = X(3,3)+tt7
val(1,1) = 5.0E-1*area(1,1)*((sqrt((invUV(2,2)*tt9+invUV(1,2)*tt8&
&)**2+(invUV(2,2)*tt6+invUV(1,2)*tt5)**2+(tt3*invUV(2,2)+tt2*invUV&
&(1,2))**2)-1.0E+0)**2+(sqrt((invUV(2,1)*tt9+invUV(1,1)*tt8)**2+(i&
&nvUV(2,1)*tt6+invUV(1,1)*tt5)**2+(tt3*invUV(2,1)+invUV(1,1)*tt2)*&
&*2)-1.0E+0)**2)
END 
SUBROUTINE bw98_stretch_jac(jac, X, invUV, area) 
IMPLICIT NONE 
REAL(KIND=8) jac(1, 9) 
REAL(KIND=8) X(3, 3) 
REAL(KIND=8) invUV(2, 2) 
REAL(KIND=8) area(1, 1) 
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
tt1 = (-invUV(2,1))-invUV(1,1)
tt2 = -X(1,1)
tt3 = X(1,2)+tt2
tt4 = X(1,3)+tt2
tt5 = tt4*invUV(2,1)+invUV(1,1)*tt3
tt6 = -X(2,1)
tt7 = X(2,2)+tt6
tt8 = X(2,3)+tt6
tt9 = invUV(2,1)*tt8+invUV(1,1)*tt7
tt10 = -X(3,1)
tt11 = X(3,2)+tt10
tt12 = X(3,3)+tt10
tt13 = invUV(2,1)*tt12+invUV(1,1)*tt11
tt14 = sqrt(tt13**2+tt9**2+tt5**2)
tt15 = 1/tt14
tt16 = tt14-1.0E+0
tt17 = (-invUV(2,2))-invUV(1,2)
tt18 = tt4*invUV(2,2)+tt3*invUV(1,2)
tt19 = invUV(2,2)*tt8+invUV(1,2)*tt7
tt20 = invUV(2,2)*tt12+invUV(1,2)*tt11
tt21 = sqrt(tt20**2+tt19**2+tt18**2)
tt22 = 1/tt21
tt23 = tt21-1.0E+0
jac(1,1) = 5.0E-1*area(1,1)*(2*tt17*tt18*tt22*tt23+2*tt1*tt5*tt15&
&*tt16)
jac(1,2) = 5.0E-1*area(1,1)*(2*tt17*tt19*tt22*tt23+2*tt1*tt9*tt15&
&*tt16)
jac(1,3) = 5.0E-1*area(1,1)*(2*tt17*tt20*tt22*tt23+2*tt1*tt13*tt1&
&5*tt16)
jac(1,4) = 5.0E-1*area(1,1)*(2*invUV(1,2)*tt18*tt22*tt23+2*invUV(&
&1,1)*tt5*tt15*tt16)
jac(1,5) = 5.0E-1*area(1,1)*(2*invUV(1,2)*tt19*tt22*tt23+2*invUV(&
&1,1)*tt9*tt15*tt16)
jac(1,6) = 5.0E-1*area(1,1)*(2*invUV(1,2)*tt20*tt22*tt23+2*invUV(&
&1,1)*tt13*tt15*tt16)
jac(1,7) = 5.0E-1*area(1,1)*(2*invUV(2,2)*tt18*tt22*tt23+2*invUV(&
&2,1)*tt5*tt15*tt16)
jac(1,8) = 5.0E-1*area(1,1)*(2*invUV(2,2)*tt19*tt22*tt23+2*invUV(&
&2,1)*tt9*tt15*tt16)
jac(1,9) = 5.0E-1*area(1,1)*(2*invUV(2,2)*tt20*tt22*tt23+2*invUV(&
&2,1)*tt13*tt15*tt16)
END 
SUBROUTINE bw98_stretch_hes(hes, X, invUV, area) 
IMPLICIT NONE 
REAL(KIND=8) hes(9, 9) 
REAL(KIND=8) X(3, 3) 
REAL(KIND=8) invUV(2, 2) 
REAL(KIND=8) area(1, 1) 
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
tt1 = (-invUV(2,1))-invUV(1,1)
tt2 = tt1**2
tt3 = -X(1,1)
tt4 = X(1,2)+tt3
tt5 = X(1,3)+tt3
tt6 = tt5*invUV(2,1)+invUV(1,1)*tt4
tt7 = tt6**2
tt8 = -X(2,1)
tt9 = X(2,2)+tt8
tt10 = X(2,3)+tt8
tt11 = invUV(2,1)*tt10+invUV(1,1)*tt9
tt12 = tt11**2
tt13 = -X(3,1)
tt14 = X(3,2)+tt13
tt15 = X(3,3)+tt13
tt16 = invUV(2,1)*tt15+invUV(1,1)*tt14
tt17 = tt16**2
tt18 = tt17+tt12+tt7
tt19 = 1/tt18
tt20 = sqrt(tt18)
tt21 = 1/tt20**3
tt22 = tt20-1.0E+0
tt23 = 1/tt20
tt24 = 2*tt2*tt23*tt22
tt25 = (-invUV(2,2))-invUV(1,2)
tt26 = tt25**2
tt27 = tt5*invUV(2,2)+tt4*invUV(1,2)
tt28 = tt27**2
tt29 = invUV(2,2)*tt10+invUV(1,2)*tt9
tt30 = tt29**2
tt31 = invUV(2,2)*tt15+invUV(1,2)*tt14
tt32 = tt31**2
tt33 = tt32+tt30+tt28
tt34 = 1/tt33
tt35 = sqrt(tt33)
tt36 = 1/tt35**3
tt37 = tt35-1.0E+0
tt38 = 1/tt35
tt39 = 2*tt26*tt38*tt37
tt40 = 5.0E-1*area(1,1)*((-2*tt26*tt27*tt29*tt36*tt37)+2*tt26*tt2&
&7*tt29*tt34-2*tt2*tt6*tt11*tt21*tt22+2*tt2*tt6*tt11*tt19)
tt41 = 5.0E-1*area(1,1)*((-2*tt26*tt27*tt31*tt36*tt37)+2*tt26*tt2&
&7*tt31*tt34-2*tt2*tt6*tt16*tt21*tt22+2*tt2*tt6*tt16*tt19)
tt42 = 2*invUV(1,1)*tt1*tt23*tt22
tt43 = 2*invUV(1,2)*tt25*tt38*tt37
tt44 = 5.0E-1*area(1,1)*(tt43-2*invUV(1,2)*tt25*tt28*tt36*tt37+2*&
&invUV(1,2)*tt25*tt28*tt34+tt42-2*invUV(1,1)*tt1*tt7*tt21*tt22+2*i&
&nvUV(1,1)*tt1*tt7*tt19)
tt45 = 5.0E-1*area(1,1)*((-2*invUV(1,2)*tt25*tt27*tt29*tt36*tt37)&
&+2*invUV(1,2)*tt25*tt27*tt29*tt34-2*invUV(1,1)*tt1*tt6*tt11*tt21*&
&tt22+2*invUV(1,1)*tt1*tt6*tt11*tt19)
tt46 = 5.0E-1*area(1,1)*((-2*invUV(1,2)*tt25*tt27*tt31*tt36*tt37)&
&+2*invUV(1,2)*tt25*tt27*tt31*tt34-2*invUV(1,1)*tt1*tt6*tt16*tt21*&
&tt22+2*invUV(1,1)*tt1*tt6*tt16*tt19)
tt47 = 2*tt1*invUV(2,1)*tt23*tt22
tt48 = 2*tt25*invUV(2,2)*tt38*tt37
tt49 = 5.0E-1*area(1,1)*(tt48-2*tt25*invUV(2,2)*tt28*tt36*tt37+2*&
&tt25*invUV(2,2)*tt28*tt34+tt47-2*tt1*invUV(2,1)*tt7*tt21*tt22+2*t&
&t1*invUV(2,1)*tt7*tt19)
tt50 = 5.0E-1*area(1,1)*((-2*tt25*invUV(2,2)*tt27*tt29*tt36*tt37)&
&+2*tt25*invUV(2,2)*tt27*tt29*tt34-2*tt1*invUV(2,1)*tt6*tt11*tt21*&
&tt22+2*tt1*invUV(2,1)*tt6*tt11*tt19)
tt51 = 5.0E-1*area(1,1)*((-2*tt25*invUV(2,2)*tt27*tt31*tt36*tt37)&
&+2*tt25*invUV(2,2)*tt27*tt31*tt34-2*tt1*invUV(2,1)*tt6*tt16*tt21*&
&tt22+2*tt1*invUV(2,1)*tt6*tt16*tt19)
tt52 = 5.0E-1*area(1,1)*((-2*tt26*tt29*tt31*tt36*tt37)+2*tt26*tt2&
&9*tt31*tt34-2*tt2*tt11*tt16*tt21*tt22+2*tt2*tt11*tt16*tt19)
tt53 = 5.0E-1*area(1,1)*(tt43-2*invUV(1,2)*tt25*tt30*tt36*tt37+2*&
&invUV(1,2)*tt25*tt30*tt34+tt42-2*invUV(1,1)*tt1*tt12*tt21*tt22+2*&
&invUV(1,1)*tt1*tt12*tt19)
tt54 = 5.0E-1*area(1,1)*((-2*invUV(1,2)*tt25*tt29*tt31*tt36*tt37)&
&+2*invUV(1,2)*tt25*tt29*tt31*tt34-2*invUV(1,1)*tt1*tt11*tt16*tt21&
&*tt22+2*invUV(1,1)*tt1*tt11*tt16*tt19)
tt55 = 5.0E-1*area(1,1)*(tt48-2*tt25*invUV(2,2)*tt30*tt36*tt37+2*&
&tt25*invUV(2,2)*tt30*tt34+tt47-2*tt1*invUV(2,1)*tt12*tt21*tt22+2*&
&tt1*invUV(2,1)*tt12*tt19)
tt56 = 5.0E-1*area(1,1)*((-2*tt25*invUV(2,2)*tt29*tt31*tt36*tt37)&
&+2*tt25*invUV(2,2)*tt29*tt31*tt34-2*tt1*invUV(2,1)*tt11*tt16*tt21&
&*tt22+2*tt1*invUV(2,1)*tt11*tt16*tt19)
tt57 = 5.0E-1*area(1,1)*(tt43-2*invUV(1,2)*tt25*tt32*tt36*tt37+2*&
&invUV(1,2)*tt25*tt32*tt34+tt42-2*invUV(1,1)*tt1*tt17*tt21*tt22+2*&
&invUV(1,1)*tt1*tt17*tt19)
tt58 = 5.0E-1*area(1,1)*(tt48-2*tt25*invUV(2,2)*tt32*tt36*tt37+2*&
&tt25*invUV(2,2)*tt32*tt34+tt47-2*tt1*invUV(2,1)*tt17*tt21*tt22+2*&
&tt1*invUV(2,1)*tt17*tt19)
tt59 = invUV(1,1)**2
tt60 = 2*tt59*tt23*tt22
tt61 = invUV(1,2)**2
tt62 = 2*tt61*tt38*tt37
tt63 = 5.0E-1*area(1,1)*((-2*tt61*tt27*tt29*tt36*tt37)+2*tt61*tt2&
&7*tt29*tt34-2*tt59*tt6*tt11*tt21*tt22+2*tt59*tt6*tt11*tt19)
tt64 = 5.0E-1*area(1,1)*((-2*tt61*tt27*tt31*tt36*tt37)+2*tt61*tt2&
&7*tt31*tt34-2*tt59*tt6*tt16*tt21*tt22+2*tt59*tt6*tt16*tt19)
tt65 = 2*invUV(1,1)*invUV(2,1)*tt23*tt22
tt66 = 2*invUV(1,2)*invUV(2,2)*tt38*tt37
tt67 = 5.0E-1*area(1,1)*(tt66-2*invUV(1,2)*invUV(2,2)*tt28*tt36*t&
&t37+2*invUV(1,2)*invUV(2,2)*tt28*tt34+tt65-2*invUV(1,1)*invUV(2,1&
&)*tt7*tt21*tt22+2*invUV(1,1)*invUV(2,1)*tt7*tt19)
tt68 = 5.0E-1*area(1,1)*((-2*invUV(1,2)*invUV(2,2)*tt27*tt29*tt36&
&*tt37)+2*invUV(1,2)*invUV(2,2)*tt27*tt29*tt34-2*invUV(1,1)*invUV(&
&2,1)*tt6*tt11*tt21*tt22+2*invUV(1,1)*invUV(2,1)*tt6*tt11*tt19)
tt69 = 5.0E-1*area(1,1)*((-2*invUV(1,2)*invUV(2,2)*tt27*tt31*tt36&
&*tt37)+2*invUV(1,2)*invUV(2,2)*tt27*tt31*tt34-2*invUV(1,1)*invUV(&
&2,1)*tt6*tt16*tt21*tt22+2*invUV(1,1)*invUV(2,1)*tt6*tt16*tt19)
tt70 = 5.0E-1*area(1,1)*((-2*tt61*tt29*tt31*tt36*tt37)+2*tt61*tt2&
&9*tt31*tt34-2*tt59*tt11*tt16*tt21*tt22+2*tt59*tt11*tt16*tt19)
tt71 = 5.0E-1*area(1,1)*(tt66-2*invUV(1,2)*invUV(2,2)*tt30*tt36*t&
&t37+2*invUV(1,2)*invUV(2,2)*tt30*tt34+tt65-2*invUV(1,1)*invUV(2,1&
&)*tt12*tt21*tt22+2*invUV(1,1)*invUV(2,1)*tt12*tt19)
tt72 = 5.0E-1*area(1,1)*((-2*invUV(1,2)*invUV(2,2)*tt29*tt31*tt36&
&*tt37)+2*invUV(1,2)*invUV(2,2)*tt29*tt31*tt34-2*invUV(1,1)*invUV(&
&2,1)*tt11*tt16*tt21*tt22+2*invUV(1,1)*invUV(2,1)*tt11*tt16*tt19)
tt73 = 5.0E-1*area(1,1)*(tt66-2*invUV(1,2)*invUV(2,2)*tt32*tt36*t&
&t37+2*invUV(1,2)*invUV(2,2)*tt32*tt34+tt65-2*invUV(1,1)*invUV(2,1&
&)*tt17*tt21*tt22+2*invUV(1,1)*invUV(2,1)*tt17*tt19)
tt74 = invUV(2,1)**2
tt75 = 2*tt74*tt23*tt22
tt76 = invUV(2,2)**2
tt77 = 2*tt76*tt38*tt37
tt78 = 5.0E-1*area(1,1)*((-2*tt76*tt27*tt29*tt36*tt37)+2*tt76*tt2&
&7*tt29*tt34-2*tt74*tt6*tt11*tt21*tt22+2*tt74*tt6*tt11*tt19)
tt79 = 5.0E-1*area(1,1)*((-2*tt76*tt27*tt31*tt36*tt37)+2*tt76*tt2&
&7*tt31*tt34-2*tt74*tt6*tt16*tt21*tt22+2*tt74*tt6*tt16*tt19)
tt80 = 5.0E-1*area(1,1)*((-2*tt76*tt29*tt31*tt36*tt37)+2*tt76*tt2&
&9*tt31*tt34-2*tt74*tt11*tt16*tt21*tt22+2*tt74*tt11*tt16*tt19)
hes(1,1) = 5.0E-1*area(1,1)*(tt39-2*tt26*tt28*tt36*tt37+2*tt26*tt&
&28*tt34+tt24-2*tt2*tt7*tt21*tt22+2*tt2*tt7*tt19)
hes(1,2) = tt40
hes(1,3) = tt41
hes(1,4) = tt44
hes(1,5) = tt45
hes(1,6) = tt46
hes(1,7) = tt49
hes(1,8) = tt50
hes(1,9) = tt51
hes(2,1) = tt40
hes(2,2) = 5.0E-1*area(1,1)*(tt39-2*tt26*tt30*tt36*tt37+2*tt26*tt&
&30*tt34+tt24-2*tt2*tt12*tt21*tt22+2*tt2*tt12*tt19)
hes(2,3) = tt52
hes(2,4) = tt45
hes(2,5) = tt53
hes(2,6) = tt54
hes(2,7) = tt50
hes(2,8) = tt55
hes(2,9) = tt56
hes(3,1) = tt41
hes(3,2) = tt52
hes(3,3) = 5.0E-1*area(1,1)*(tt39-2*tt26*tt32*tt36*tt37+2*tt26*tt&
&32*tt34+tt24-2*tt2*tt17*tt21*tt22+2*tt2*tt17*tt19)
hes(3,4) = tt46
hes(3,5) = tt54
hes(3,6) = tt57
hes(3,7) = tt51
hes(3,8) = tt56
hes(3,9) = tt58
hes(4,1) = tt44
hes(4,2) = tt45
hes(4,3) = tt46
hes(4,4) = 5.0E-1*area(1,1)*(tt62-2*tt61*tt28*tt36*tt37+2*tt61*tt&
&28*tt34+tt60-2*tt59*tt7*tt21*tt22+2*tt59*tt7*tt19)
hes(4,5) = tt63
hes(4,6) = tt64
hes(4,7) = tt67
hes(4,8) = tt68
hes(4,9) = tt69
hes(5,1) = tt45
hes(5,2) = tt53
hes(5,3) = tt54
hes(5,4) = tt63
hes(5,5) = 5.0E-1*area(1,1)*(tt62-2*tt61*tt30*tt36*tt37+2*tt61*tt&
&30*tt34+tt60-2*tt59*tt12*tt21*tt22+2*tt59*tt12*tt19)
hes(5,6) = tt70
hes(5,7) = tt68
hes(5,8) = tt71
hes(5,9) = tt72
hes(6,1) = tt46
hes(6,2) = tt54
hes(6,3) = tt57
hes(6,4) = tt64
hes(6,5) = tt70
hes(6,6) = 5.0E-1*area(1,1)*(tt62-2*tt61*tt32*tt36*tt37+2*tt61*tt&
&32*tt34+tt60-2*tt59*tt17*tt21*tt22+2*tt59*tt17*tt19)
hes(6,7) = tt69
hes(6,8) = tt72
hes(6,9) = tt73
hes(7,1) = tt49
hes(7,2) = tt50
hes(7,3) = tt51
hes(7,4) = tt67
hes(7,5) = tt68
hes(7,6) = tt69
hes(7,7) = 5.0E-1*area(1,1)*(tt77-2*tt76*tt28*tt36*tt37+2*tt76*tt&
&28*tt34+tt75-2*tt74*tt7*tt21*tt22+2*tt74*tt7*tt19)
hes(7,8) = tt78
hes(7,9) = tt79
hes(8,1) = tt50
hes(8,2) = tt55
hes(8,3) = tt56
hes(8,4) = tt68
hes(8,5) = tt71
hes(8,6) = tt72
hes(8,7) = tt78
hes(8,8) = 5.0E-1*area(1,1)*(tt77-2*tt76*tt30*tt36*tt37+2*tt76*tt&
&30*tt34+tt75-2*tt74*tt12*tt21*tt22+2*tt74*tt12*tt19)
hes(8,9) = tt80
hes(9,1) = tt51
hes(9,2) = tt56
hes(9,3) = tt58
hes(9,4) = tt69
hes(9,5) = tt72
hes(9,6) = tt73
hes(9,7) = tt79
hes(9,8) = tt80
hes(9,9) = 5.0E-1*area(1,1)*(tt77-2*tt76*tt32*tt36*tt37+2*tt76*tt&
&32*tt34+tt75-2*tt74*tt17*tt21*tt22+2*tt74*tt17*tt19)
END 
SUBROUTINE bw98_shear(val, X, invUV, area) 
IMPLICIT NONE 
REAL(KIND=8) val(1, 1) 
REAL(KIND=8) X(3, 3) 
REAL(KIND=8) invUV(2, 2) 
REAL(KIND=8) area(1, 1) 
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
tt1 = -X(1,1)
tt2 = X(1,2)+tt1
tt3 = invUV(1,1)*tt2
tt4 = X(1,3)+tt1
tt5 = tt4*invUV(2,1)
tt6 = -X(2,1)
tt7 = X(2,2)+tt6
tt8 = invUV(1,1)*tt7
tt9 = X(2,3)+tt6
tt10 = invUV(2,1)*tt9
tt11 = -X(3,1)
tt12 = X(3,2)+tt11
tt13 = invUV(1,1)*tt12
tt14 = X(3,3)+tt11
tt15 = invUV(2,1)*tt14
val(1,1) = 5.0E-1*area(1,1)*((sqrt((invUV(2,2)*tt14+tt15+invUV(1,&
&2)*tt12+tt13)**2+(invUV(2,2)*tt9+tt10+invUV(1,2)*tt7+tt8)**2+(tt4&
&*invUV(2,2)+tt5+tt2*invUV(1,2)+tt3)**2)-1.4142135623730952E+0)**2&
&+(sqrt(((-invUV(2,2)*tt14)+tt15-invUV(1,2)*tt12+tt13)**2+((-invUV&
&(2,2)*tt9)+tt10-invUV(1,2)*tt7+tt8)**2+((-tt4*invUV(2,2))+tt5-tt2&
&*invUV(1,2)+tt3)**2)-1.4142135623730952E+0)**2)
END 
SUBROUTINE bw98_shear_jac(jac, X, invUV, area) 
IMPLICIT NONE 
REAL(KIND=8) jac(1, 9) 
REAL(KIND=8) X(3, 3) 
REAL(KIND=8) invUV(2, 2) 
REAL(KIND=8) area(1, 1) 
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
tt1 = -invUV(1,1)
tt2 = -invUV(2,1)
tt3 = invUV(2,2)+tt2+invUV(1,2)+tt1
tt4 = -X(1,1)
tt5 = X(1,2)+tt4
tt6 = invUV(1,1)*tt5
tt7 = X(1,3)+tt4
tt8 = tt7*invUV(2,1)
tt9 = (-tt7*invUV(2,2))+tt8-tt5*invUV(1,2)+tt6
tt10 = -X(2,1)
tt11 = X(2,2)+tt10
tt12 = invUV(1,1)*tt11
tt13 = X(2,3)+tt10
tt14 = invUV(2,1)*tt13
tt15 = (-invUV(2,2)*tt13)+tt14-invUV(1,2)*tt11+tt12
tt16 = -X(3,1)
tt17 = X(3,2)+tt16
tt18 = invUV(1,1)*tt17
tt19 = X(3,3)+tt16
tt20 = invUV(2,1)*tt19
tt21 = (-invUV(2,2)*tt19)+tt20-invUV(1,2)*tt17+tt18
tt22 = sqrt(tt21**2+tt15**2+tt9**2)
tt23 = 1/tt22
tt24 = tt22-1.4142135623730952E+0
tt25 = -invUV(1,2)
tt26 = -invUV(2,2)
tt27 = tt26+tt2+tt25+tt1
tt28 = tt7*invUV(2,2)+tt8+tt5*invUV(1,2)+tt6
tt29 = invUV(2,2)*tt13+tt14+invUV(1,2)*tt11+tt12
tt30 = invUV(2,2)*tt19+tt20+invUV(1,2)*tt17+tt18
tt31 = sqrt(tt30**2+tt29**2+tt28**2)
tt32 = 1/tt31
tt33 = tt31-1.4142135623730952E+0
tt34 = tt25+invUV(1,1)
tt35 = invUV(1,2)+invUV(1,1)
tt36 = tt26+invUV(2,1)
tt37 = invUV(2,2)+invUV(2,1)
jac(1,1) = 5.0E-1*area(1,1)*(2*tt27*tt28*tt32*tt33+2*tt3*tt9*tt23&
&*tt24)
jac(1,2) = 5.0E-1*area(1,1)*(2*tt27*tt29*tt32*tt33+2*tt3*tt15*tt2&
&3*tt24)
jac(1,3) = 5.0E-1*area(1,1)*(2*tt27*tt30*tt32*tt33+2*tt3*tt21*tt2&
&3*tt24)
jac(1,4) = 5.0E-1*area(1,1)*(2*tt35*tt28*tt32*tt33+2*tt34*tt9*tt2&
&3*tt24)
jac(1,5) = 5.0E-1*area(1,1)*(2*tt35*tt29*tt32*tt33+2*tt34*tt15*tt&
&23*tt24)
jac(1,6) = 5.0E-1*area(1,1)*(2*tt35*tt30*tt32*tt33+2*tt34*tt21*tt&
&23*tt24)
jac(1,7) = 5.0E-1*area(1,1)*(2*tt37*tt28*tt32*tt33+2*tt36*tt9*tt2&
&3*tt24)
jac(1,8) = 5.0E-1*area(1,1)*(2*tt37*tt29*tt32*tt33+2*tt36*tt15*tt&
&23*tt24)
jac(1,9) = 5.0E-1*area(1,1)*(2*tt37*tt30*tt32*tt33+2*tt36*tt21*tt&
&23*tt24)
END 
SUBROUTINE bw98_shear_hes(hes, X, invUV, area) 
IMPLICIT NONE 
REAL(KIND=8) hes(9, 9) 
REAL(KIND=8) X(3, 3) 
REAL(KIND=8) invUV(2, 2) 
REAL(KIND=8) area(1, 1) 
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
tt1 = -invUV(1,1)
tt2 = -invUV(2,1)
tt3 = invUV(2,2)+tt2+invUV(1,2)+tt1
tt4 = tt3**2
tt5 = -X(1,1)
tt6 = X(1,2)+tt5
tt7 = invUV(1,1)*tt6
tt8 = X(1,3)+tt5
tt9 = tt8*invUV(2,1)
tt10 = (-tt8*invUV(2,2))+tt9-tt6*invUV(1,2)+tt7
tt11 = tt10**2
tt12 = -X(2,1)
tt13 = X(2,2)+tt12
tt14 = invUV(1,1)*tt13
tt15 = X(2,3)+tt12
tt16 = invUV(2,1)*tt15
tt17 = (-invUV(2,2)*tt15)+tt16-invUV(1,2)*tt13+tt14
tt18 = tt17**2
tt19 = -X(3,1)
tt20 = X(3,2)+tt19
tt21 = invUV(1,1)*tt20
tt22 = X(3,3)+tt19
tt23 = invUV(2,1)*tt22
tt24 = (-invUV(2,2)*tt22)+tt23-invUV(1,2)*tt20+tt21
tt25 = tt24**2
tt26 = tt25+tt18+tt11
tt27 = 1/tt26
tt28 = sqrt(tt26)
tt29 = 1/tt28**3
tt30 = tt28-1.4142135623730952E+0
tt31 = 1/tt28
tt32 = 2*tt4*tt31*tt30
tt33 = -invUV(1,2)
tt34 = -invUV(2,2)
tt35 = tt34+tt2+tt33+tt1
tt36 = tt35**2
tt37 = tt8*invUV(2,2)+tt9+tt6*invUV(1,2)+tt7
tt38 = tt37**2
tt39 = invUV(2,2)*tt15+tt16+invUV(1,2)*tt13+tt14
tt40 = tt39**2
tt41 = invUV(2,2)*tt22+tt23+invUV(1,2)*tt20+tt21
tt42 = tt41**2
tt43 = tt42+tt40+tt38
tt44 = 1/tt43
tt45 = sqrt(tt43)
tt46 = 1/tt45**3
tt47 = tt45-1.4142135623730952E+0
tt48 = 1/tt45
tt49 = 2*tt36*tt48*tt47
tt50 = 5.0E-1*area(1,1)*((-2*tt36*tt37*tt39*tt46*tt47)+2*tt36*tt3&
&7*tt39*tt44-2*tt4*tt10*tt17*tt29*tt30+2*tt4*tt10*tt17*tt27)
tt51 = 5.0E-1*area(1,1)*((-2*tt36*tt37*tt41*tt46*tt47)+2*tt36*tt3&
&7*tt41*tt44-2*tt4*tt10*tt24*tt29*tt30+2*tt4*tt10*tt24*tt27)
tt52 = tt33+invUV(1,1)
tt53 = 2*tt52*tt3*tt31*tt30
tt54 = invUV(1,2)+invUV(1,1)
tt55 = 2*tt54*tt35*tt48*tt47
tt56 = 5.0E-1*area(1,1)*(tt55-2*tt54*tt35*tt38*tt46*tt47+2*tt54*t&
&t35*tt38*tt44+tt53-2*tt52*tt3*tt11*tt29*tt30+2*tt52*tt3*tt11*tt27&
&)
tt57 = 5.0E-1*area(1,1)*((-2*tt54*tt35*tt37*tt39*tt46*tt47)+2*tt5&
&4*tt35*tt37*tt39*tt44-2*tt52*tt3*tt10*tt17*tt29*tt30+2*tt52*tt3*t&
&t10*tt17*tt27)
tt58 = 5.0E-1*area(1,1)*((-2*tt54*tt35*tt37*tt41*tt46*tt47)+2*tt5&
&4*tt35*tt37*tt41*tt44-2*tt52*tt3*tt10*tt24*tt29*tt30+2*tt52*tt3*t&
&t10*tt24*tt27)
tt59 = tt34+invUV(2,1)
tt60 = 2*tt59*tt3*tt31*tt30
tt61 = invUV(2,2)+invUV(2,1)
tt62 = 2*tt35*tt61*tt48*tt47
tt63 = 5.0E-1*area(1,1)*(tt62-2*tt35*tt61*tt38*tt46*tt47+2*tt35*t&
&t61*tt38*tt44+tt60-2*tt59*tt3*tt11*tt29*tt30+2*tt59*tt3*tt11*tt27&
&)
tt64 = 5.0E-1*area(1,1)*((-2*tt35*tt61*tt37*tt39*tt46*tt47)+2*tt3&
&5*tt61*tt37*tt39*tt44-2*tt59*tt3*tt10*tt17*tt29*tt30+2*tt59*tt3*t&
&t10*tt17*tt27)
tt65 = 5.0E-1*area(1,1)*((-2*tt35*tt61*tt37*tt41*tt46*tt47)+2*tt3&
&5*tt61*tt37*tt41*tt44-2*tt59*tt3*tt10*tt24*tt29*tt30+2*tt59*tt3*t&
&t10*tt24*tt27)
tt66 = 5.0E-1*area(1,1)*((-2*tt36*tt39*tt41*tt46*tt47)+2*tt36*tt3&
&9*tt41*tt44-2*tt4*tt17*tt24*tt29*tt30+2*tt4*tt17*tt24*tt27)
tt67 = 5.0E-1*area(1,1)*(tt55-2*tt54*tt35*tt40*tt46*tt47+2*tt54*t&
&t35*tt40*tt44+tt53-2*tt52*tt3*tt18*tt29*tt30+2*tt52*tt3*tt18*tt27&
&)
tt68 = 5.0E-1*area(1,1)*((-2*tt54*tt35*tt39*tt41*tt46*tt47)+2*tt5&
&4*tt35*tt39*tt41*tt44-2*tt52*tt3*tt17*tt24*tt29*tt30+2*tt52*tt3*t&
&t17*tt24*tt27)
tt69 = 5.0E-1*area(1,1)*(tt62-2*tt35*tt61*tt40*tt46*tt47+2*tt35*t&
&t61*tt40*tt44+tt60-2*tt59*tt3*tt18*tt29*tt30+2*tt59*tt3*tt18*tt27&
&)
tt70 = 5.0E-1*area(1,1)*((-2*tt35*tt61*tt39*tt41*tt46*tt47)+2*tt3&
&5*tt61*tt39*tt41*tt44-2*tt59*tt3*tt17*tt24*tt29*tt30+2*tt59*tt3*t&
&t17*tt24*tt27)
tt71 = 5.0E-1*area(1,1)*(tt55-2*tt54*tt35*tt42*tt46*tt47+2*tt54*t&
&t35*tt42*tt44+tt53-2*tt52*tt3*tt25*tt29*tt30+2*tt52*tt3*tt25*tt27&
&)
tt72 = 5.0E-1*area(1,1)*(tt62-2*tt35*tt61*tt42*tt46*tt47+2*tt35*t&
&t61*tt42*tt44+tt60-2*tt59*tt3*tt25*tt29*tt30+2*tt59*tt3*tt25*tt27&
&)
tt73 = tt52**2
tt74 = 2*tt73*tt31*tt30
tt75 = tt54**2
tt76 = 2*tt75*tt48*tt47
tt77 = 5.0E-1*area(1,1)*((-2*tt75*tt37*tt39*tt46*tt47)+2*tt75*tt3&
&7*tt39*tt44-2*tt73*tt10*tt17*tt29*tt30+2*tt73*tt10*tt17*tt27)
tt78 = 5.0E-1*area(1,1)*((-2*tt75*tt37*tt41*tt46*tt47)+2*tt75*tt3&
&7*tt41*tt44-2*tt73*tt10*tt24*tt29*tt30+2*tt73*tt10*tt24*tt27)
tt79 = 2*tt52*tt59*tt31*tt30
tt80 = 2*tt54*tt61*tt48*tt47
tt81 = 5.0E-1*area(1,1)*(tt80-2*tt54*tt61*tt38*tt46*tt47+2*tt54*t&
&t61*tt38*tt44+tt79-2*tt52*tt59*tt11*tt29*tt30+2*tt52*tt59*tt11*tt&
&27)
tt82 = 5.0E-1*area(1,1)*((-2*tt54*tt61*tt37*tt39*tt46*tt47)+2*tt5&
&4*tt61*tt37*tt39*tt44-2*tt52*tt59*tt10*tt17*tt29*tt30+2*tt52*tt59&
&*tt10*tt17*tt27)
tt83 = 5.0E-1*area(1,1)*((-2*tt54*tt61*tt37*tt41*tt46*tt47)+2*tt5&
&4*tt61*tt37*tt41*tt44-2*tt52*tt59*tt10*tt24*tt29*tt30+2*tt52*tt59&
&*tt10*tt24*tt27)
tt84 = 5.0E-1*area(1,1)*((-2*tt75*tt39*tt41*tt46*tt47)+2*tt75*tt3&
&9*tt41*tt44-2*tt73*tt17*tt24*tt29*tt30+2*tt73*tt17*tt24*tt27)
tt85 = 5.0E-1*area(1,1)*(tt80-2*tt54*tt61*tt40*tt46*tt47+2*tt54*t&
&t61*tt40*tt44+tt79-2*tt52*tt59*tt18*tt29*tt30+2*tt52*tt59*tt18*tt&
&27)
tt86 = 5.0E-1*area(1,1)*((-2*tt54*tt61*tt39*tt41*tt46*tt47)+2*tt5&
&4*tt61*tt39*tt41*tt44-2*tt52*tt59*tt17*tt24*tt29*tt30+2*tt52*tt59&
&*tt17*tt24*tt27)
tt87 = 5.0E-1*area(1,1)*(tt80-2*tt54*tt61*tt42*tt46*tt47+2*tt54*t&
&t61*tt42*tt44+tt79-2*tt52*tt59*tt25*tt29*tt30+2*tt52*tt59*tt25*tt&
&27)
tt88 = tt59**2
tt89 = 2*tt88*tt31*tt30
tt90 = tt61**2
tt91 = 2*tt90*tt48*tt47
tt92 = 5.0E-1*area(1,1)*((-2*tt90*tt37*tt39*tt46*tt47)+2*tt90*tt3&
&7*tt39*tt44-2*tt88*tt10*tt17*tt29*tt30+2*tt88*tt10*tt17*tt27)
tt93 = 5.0E-1*area(1,1)*((-2*tt90*tt37*tt41*tt46*tt47)+2*tt90*tt3&
&7*tt41*tt44-2*tt88*tt10*tt24*tt29*tt30+2*tt88*tt10*tt24*tt27)
tt94 = 5.0E-1*area(1,1)*((-2*tt90*tt39*tt41*tt46*tt47)+2*tt90*tt3&
&9*tt41*tt44-2*tt88*tt17*tt24*tt29*tt30+2*tt88*tt17*tt24*tt27)
hes(1,1) = 5.0E-1*area(1,1)*(tt49-2*tt36*tt38*tt46*tt47+2*tt36*tt&
&38*tt44+tt32-2*tt4*tt11*tt29*tt30+2*tt4*tt11*tt27)
hes(1,2) = tt50
hes(1,3) = tt51
hes(1,4) = tt56
hes(1,5) = tt57
hes(1,6) = tt58
hes(1,7) = tt63
hes(1,8) = tt64
hes(1,9) = tt65
hes(2,1) = tt50
hes(2,2) = 5.0E-1*area(1,1)*(tt49-2*tt36*tt40*tt46*tt47+2*tt36*tt&
&40*tt44+tt32-2*tt4*tt18*tt29*tt30+2*tt4*tt18*tt27)
hes(2,3) = tt66
hes(2,4) = tt57
hes(2,5) = tt67
hes(2,6) = tt68
hes(2,7) = tt64
hes(2,8) = tt69
hes(2,9) = tt70
hes(3,1) = tt51
hes(3,2) = tt66
hes(3,3) = 5.0E-1*area(1,1)*(tt49-2*tt36*tt42*tt46*tt47+2*tt36*tt&
&42*tt44+tt32-2*tt4*tt25*tt29*tt30+2*tt4*tt25*tt27)
hes(3,4) = tt58
hes(3,5) = tt68
hes(3,6) = tt71
hes(3,7) = tt65
hes(3,8) = tt70
hes(3,9) = tt72
hes(4,1) = tt56
hes(4,2) = tt57
hes(4,3) = tt58
hes(4,4) = 5.0E-1*area(1,1)*(tt76-2*tt75*tt38*tt46*tt47+2*tt75*tt&
&38*tt44+tt74-2*tt73*tt11*tt29*tt30+2*tt73*tt11*tt27)
hes(4,5) = tt77
hes(4,6) = tt78
hes(4,7) = tt81
hes(4,8) = tt82
hes(4,9) = tt83
hes(5,1) = tt57
hes(5,2) = tt67
hes(5,3) = tt68
hes(5,4) = tt77
hes(5,5) = 5.0E-1*area(1,1)*(tt76-2*tt75*tt40*tt46*tt47+2*tt75*tt&
&40*tt44+tt74-2*tt73*tt18*tt29*tt30+2*tt73*tt18*tt27)
hes(5,6) = tt84
hes(5,7) = tt82
hes(5,8) = tt85
hes(5,9) = tt86
hes(6,1) = tt58
hes(6,2) = tt68
hes(6,3) = tt71
hes(6,4) = tt78
hes(6,5) = tt84
hes(6,6) = 5.0E-1*area(1,1)*(tt76-2*tt75*tt42*tt46*tt47+2*tt75*tt&
&42*tt44+tt74-2*tt73*tt25*tt29*tt30+2*tt73*tt25*tt27)
hes(6,7) = tt83
hes(6,8) = tt86
hes(6,9) = tt87
hes(7,1) = tt63
hes(7,2) = tt64
hes(7,3) = tt65
hes(7,4) = tt81
hes(7,5) = tt82
hes(7,6) = tt83
hes(7,7) = 5.0E-1*area(1,1)*(tt91-2*tt90*tt38*tt46*tt47+2*tt90*tt&
&38*tt44+tt89-2*tt88*tt11*tt29*tt30+2*tt88*tt11*tt27)
hes(7,8) = tt92
hes(7,9) = tt93
hes(8,1) = tt64
hes(8,2) = tt69
hes(8,3) = tt70
hes(8,4) = tt82
hes(8,5) = tt85
hes(8,6) = tt86
hes(8,7) = tt92
hes(8,8) = 5.0E-1*area(1,1)*(tt91-2*tt90*tt40*tt46*tt47+2*tt90*tt&
&40*tt44+tt89-2*tt88*tt18*tt29*tt30+2*tt88*tt18*tt27)
hes(8,9) = tt94
hes(9,1) = tt65
hes(9,2) = tt70
hes(9,3) = tt72
hes(9,4) = tt83
hes(9,5) = tt86
hes(9,6) = tt87
hes(9,7) = tt93
hes(9,8) = tt94
hes(9,9) = 5.0E-1*area(1,1)*(tt91-2*tt90*tt42*tt46*tt47+2*tt90*tt&
&42*tt44+tt89-2*tt88*tt25*tt29*tt30+2*tt88*tt25*tt27)
END 
SUBROUTINE fem_stretch(val, X, Dm, area, k) 
IMPLICIT NONE 
REAL(KIND=8) val(1, 1) 
REAL(KIND=8) X(3, 3) 
REAL(KIND=8) Dm(2, 2) 
REAL(KIND=8) area(1, 1) 
REAL(KIND=8) k(4, 1) 
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
tt1 = -X(1,1)
tt2 = X(1,2)+tt1
tt3 = X(1,3)+tt1
tt4 = tt3*Dm(2,1)+Dm(1,1)*tt2
tt5 = -X(2,1)
tt6 = X(2,2)+tt5
tt7 = X(2,3)+tt5
tt8 = Dm(2,1)*tt7+Dm(1,1)*tt6
tt9 = -X(3,1)
tt10 = X(3,2)+tt9
tt11 = X(3,3)+tt9
tt12 = Dm(2,1)*tt11+Dm(1,1)*tt10
tt13 = tt12**2+tt8**2+tt4**2-1
tt14 = tt3*Dm(2,2)+Dm(1,2)*tt2
tt15 = Dm(2,2)*tt7+Dm(1,2)*tt6
tt16 = Dm(2,2)*tt11+Dm(1,2)*tt10
tt17 = tt16**2+tt15**2+tt14**2-1
val(1,1) = 5.0E-1*area(1,1)*(2.5E-1*(tt12*tt16+tt8*tt15+tt4*tt14)&
&**2*k(4,1)+2.5E-1*k(3,1)*tt17**2+2.5E-1*k(1,1)*tt13**2+5.0E-1*k(2&
&,1)*tt13*tt17)
END 
SUBROUTINE fem_stretch_jac(jac, X, Dm, area, k) 
IMPLICIT NONE 
REAL(KIND=8) jac(1, 9) 
REAL(KIND=8) X(3, 3) 
REAL(KIND=8) Dm(2, 2) 
REAL(KIND=8) area(1, 1) 
REAL(KIND=8) k(4, 1) 
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
tt1 = (-Dm(2,1))-Dm(1,1)
tt2 = -X(1,1)
tt3 = X(1,2)+tt2
tt4 = X(1,3)+tt2
tt5 = tt4*Dm(2,1)+Dm(1,1)*tt3
tt6 = -X(2,1)
tt7 = X(2,2)+tt6
tt8 = X(2,3)+tt6
tt9 = Dm(2,1)*tt8+Dm(1,1)*tt7
tt10 = -X(3,1)
tt11 = X(3,2)+tt10
tt12 = X(3,3)+tt10
tt13 = Dm(2,1)*tt12+Dm(1,1)*tt11
tt14 = tt13**2+tt9**2+tt5**2-1
tt15 = (-Dm(2,2))-Dm(1,2)
tt16 = tt4*Dm(2,2)+Dm(1,2)*tt3
tt17 = Dm(2,2)*tt8+Dm(1,2)*tt7
tt18 = Dm(2,2)*tt12+Dm(1,2)*tt11
tt19 = tt18**2+tt17**2+tt16**2-1
tt20 = tt13*tt18+tt9*tt17+tt5*tt16
jac(1,1) = 5.0E-1*area(1,1)*(5.0E-1*(tt1*tt16+tt5*tt15)*tt20*k(4,&
&1)+1.0E+0*tt15*tt16*k(3,1)*tt19+1.0E+0*tt1*tt5*k(2,1)*tt19+1.0E+0&
&*k(2,1)*tt15*tt16*tt14+1.0E+0*k(1,1)*tt1*tt5*tt14)
jac(1,2) = 5.0E-1*area(1,1)*(5.0E-1*(tt1*tt17+tt15*tt9)*tt20*k(4,&
&1)+1.0E+0*tt15*tt17*k(3,1)*tt19+1.0E+0*tt1*k(2,1)*tt9*tt19+1.0E+0&
&*k(2,1)*tt15*tt17*tt14+1.0E+0*k(1,1)*tt1*tt9*tt14)
jac(1,3) = 5.0E-1*area(1,1)*(5.0E-1*(tt1*tt18+tt15*tt13)*tt20*k(4&
&,1)+1.0E+0*tt15*k(3,1)*tt18*tt19+1.0E+0*tt1*k(2,1)*tt13*tt19+1.0E&
&+0*k(2,1)*tt15*tt18*tt14+1.0E+0*k(1,1)*tt1*tt13*tt14)
jac(1,4) = 5.0E-1*area(1,1)*(5.0E-1*(Dm(1,1)*tt16+Dm(1,2)*tt5)*tt&
&20*k(4,1)+1.0E+0*Dm(1,2)*tt16*k(3,1)*tt19+1.0E+0*Dm(1,1)*tt5*k(2,&
&1)*tt19+1.0E+0*Dm(1,2)*k(2,1)*tt16*tt14+1.0E+0*Dm(1,1)*k(1,1)*tt5&
&*tt14)
jac(1,5) = 5.0E-1*area(1,1)*(5.0E-1*(Dm(1,1)*tt17+Dm(1,2)*tt9)*tt&
&20*k(4,1)+1.0E+0*Dm(1,2)*tt17*k(3,1)*tt19+1.0E+0*Dm(1,1)*k(2,1)*t&
&t9*tt19+1.0E+0*Dm(1,2)*k(2,1)*tt17*tt14+1.0E+0*Dm(1,1)*k(1,1)*tt9&
&*tt14)
jac(1,6) = 5.0E-1*area(1,1)*(5.0E-1*(Dm(1,1)*tt18+Dm(1,2)*tt13)*t&
&t20*k(4,1)+1.0E+0*Dm(1,2)*k(3,1)*tt18*tt19+1.0E+0*Dm(1,1)*k(2,1)*&
&tt13*tt19+1.0E+0*Dm(1,2)*k(2,1)*tt18*tt14+1.0E+0*Dm(1,1)*k(1,1)*t&
&t13*tt14)
jac(1,7) = 5.0E-1*area(1,1)*(5.0E-1*(Dm(2,1)*tt16+tt5*Dm(2,2))*tt&
&20*k(4,1)+1.0E+0*Dm(2,2)*tt16*k(3,1)*tt19+1.0E+0*Dm(2,1)*tt5*k(2,&
&1)*tt19+1.0E+0*k(2,1)*Dm(2,2)*tt16*tt14+1.0E+0*k(1,1)*Dm(2,1)*tt5&
&*tt14)
jac(1,8) = 5.0E-1*area(1,1)*(5.0E-1*(Dm(2,1)*tt17+Dm(2,2)*tt9)*tt&
&20*k(4,1)+1.0E+0*Dm(2,2)*tt17*k(3,1)*tt19+1.0E+0*Dm(2,1)*k(2,1)*t&
&t9*tt19+1.0E+0*k(2,1)*Dm(2,2)*tt17*tt14+1.0E+0*k(1,1)*Dm(2,1)*tt9&
&*tt14)
jac(1,9) = 5.0E-1*area(1,1)*(5.0E-1*(Dm(2,1)*tt18+Dm(2,2)*tt13)*t&
&t20*k(4,1)+1.0E+0*Dm(2,2)*k(3,1)*tt18*tt19+1.0E+0*Dm(2,1)*k(2,1)*&
&tt13*tt19+1.0E+0*k(2,1)*Dm(2,2)*tt18*tt14+1.0E+0*k(1,1)*Dm(2,1)*t&
&t13*tt14)
END 
SUBROUTINE fem_stretch_hes(hes, X, Dm, area, k) 
IMPLICIT NONE 
REAL(KIND=8) hes(9, 9) 
REAL(KIND=8) X(3, 3) 
REAL(KIND=8) Dm(2, 2) 
REAL(KIND=8) area(1, 1) 
REAL(KIND=8) k(4, 1) 
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
tt1 = (-Dm(2,1))-Dm(1,1)
tt2 = tt1**2
tt3 = -X(1,1)
tt4 = X(1,2)+tt3
tt5 = X(1,3)+tt3
tt6 = tt5*Dm(2,1)+Dm(1,1)*tt4
tt7 = tt6**2
tt8 = (-Dm(2,2))-Dm(1,2)
tt9 = tt5*Dm(2,2)+Dm(1,2)*tt4
tt10 = tt8**2
tt11 = tt9**2
tt12 = -X(2,1)
tt13 = X(2,2)+tt12
tt14 = X(2,3)+tt12
tt15 = Dm(2,1)*tt14+Dm(1,1)*tt13
tt16 = tt15**2
tt17 = -X(3,1)
tt18 = X(3,2)+tt17
tt19 = X(3,3)+tt17
tt20 = Dm(2,1)*tt19+Dm(1,1)*tt18
tt21 = tt20**2
tt22 = tt21+tt16+tt7-1
tt23 = 1.0E+0*k(1,1)*tt2*tt22
tt24 = 1.0E+0*k(2,1)*tt10*tt22
tt25 = Dm(2,2)*tt14+Dm(1,2)*tt13
tt26 = tt25**2
tt27 = Dm(2,2)*tt19+Dm(1,2)*tt18
tt28 = tt27**2
tt29 = tt28+tt26+tt11-1
tt30 = 1.0E+0*tt2*k(2,1)*tt29
tt31 = 1.0E+0*tt10*k(3,1)*tt29
tt32 = tt1*tt9+tt6*tt8
tt33 = tt20*tt27+tt15*tt25+tt6*tt9
tt34 = 1.0E+0*tt1*tt8*tt33*k(4,1)
tt35 = tt1*tt25+tt8*tt15
tt36 = 5.0E-1*area(1,1)*(5.0E-1*tt32*tt35*k(4,1)+2.0E+0*tt10*tt9*&
&tt25*k(3,1)+2.0E+0*tt1*tt6*k(2,1)*tt8*tt25+2.0E+0*tt1*k(2,1)*tt8*&
&tt9*tt15+2.0E+0*k(1,1)*tt2*tt6*tt15)
tt37 = tt1*tt27+tt8*tt20
tt38 = 5.0E-1*area(1,1)*(5.0E-1*tt32*tt37*k(4,1)+2.0E+0*tt10*tt9*&
&k(3,1)*tt27+2.0E+0*tt1*tt6*k(2,1)*tt8*tt27+2.0E+0*tt1*k(2,1)*tt8*&
&tt9*tt20+2.0E+0*k(1,1)*tt2*tt6*tt20)
tt39 = 1.0E+0*Dm(1,1)*k(1,1)*tt1*tt22
tt40 = 1.0E+0*Dm(1,2)*k(2,1)*tt8*tt22
tt41 = 1.0E+0*Dm(1,1)*tt1*k(2,1)*tt29
tt42 = 1.0E+0*Dm(1,2)*tt8*k(3,1)*tt29
tt43 = Dm(1,1)*tt9+Dm(1,2)*tt6
tt44 = 5.0E-1*(Dm(1,1)*tt8+Dm(1,2)*tt1)*tt33*k(4,1)
tt45 = 5.0E-1*area(1,1)*(tt44+5.0E-1*tt43*tt32*k(4,1)+tt42+tt41+t&
&t40+tt39+2.0E+0*Dm(1,2)*tt8*tt11*k(3,1)+2.0E+0*Dm(1,1)*tt6*k(2,1)&
&*tt8*tt9+2.0E+0*Dm(1,2)*tt1*tt6*k(2,1)*tt9+2.0E+0*Dm(1,1)*k(1,1)*&
&tt1*tt7)
tt46 = 2.0E+0*Dm(1,1)*k(1,1)*tt1*tt6*tt15
tt47 = 2.0E+0*Dm(1,2)*tt8*tt9*tt25*k(3,1)
tt48 = Dm(1,1)*tt25+Dm(1,2)*tt15
tt49 = 5.0E-1*area(1,1)*(5.0E-1*tt32*tt48*k(4,1)+tt47+2.0E+0*Dm(1&
&,2)*tt1*tt6*k(2,1)*tt25+2.0E+0*Dm(1,1)*k(2,1)*tt8*tt9*tt15+tt46)
tt50 = 2.0E+0*Dm(1,1)*k(1,1)*tt1*tt6*tt20
tt51 = 2.0E+0*Dm(1,2)*tt8*tt9*k(3,1)*tt27
tt52 = Dm(1,1)*tt27+Dm(1,2)*tt20
tt53 = 5.0E-1*area(1,1)*(5.0E-1*tt32*tt52*k(4,1)+tt51+2.0E+0*Dm(1&
&,2)*tt1*tt6*k(2,1)*tt27+2.0E+0*Dm(1,1)*k(2,1)*tt8*tt9*tt20+tt50)
tt54 = 1.0E+0*k(1,1)*tt1*Dm(2,1)*tt22
tt55 = 1.0E+0*k(2,1)*tt8*Dm(2,2)*tt22
tt56 = 1.0E+0*tt1*Dm(2,1)*k(2,1)*tt29
tt57 = 1.0E+0*tt8*Dm(2,2)*k(3,1)*tt29
tt58 = Dm(2,1)*tt9+tt6*Dm(2,2)
tt59 = 5.0E-1*(tt1*Dm(2,2)+Dm(2,1)*tt8)*tt33*k(4,1)
tt60 = 5.0E-1*area(1,1)*(tt59+5.0E-1*tt32*tt58*k(4,1)+tt57+tt56+t&
&t55+tt54+2.0E+0*tt8*Dm(2,2)*tt11*k(3,1)+2.0E+0*tt1*tt6*k(2,1)*Dm(&
&2,2)*tt9+2.0E+0*Dm(2,1)*tt6*k(2,1)*tt8*tt9+2.0E+0*k(1,1)*tt1*Dm(2&
&,1)*tt7)
tt61 = 2.0E+0*k(1,1)*tt1*Dm(2,1)*tt6*tt15
tt62 = 2.0E+0*tt8*Dm(2,2)*tt9*tt25*k(3,1)
tt63 = Dm(2,1)*tt25+Dm(2,2)*tt15
tt64 = 5.0E-1*area(1,1)*(5.0E-1*tt32*tt63*k(4,1)+tt62+2.0E+0*tt1*&
&tt6*k(2,1)*Dm(2,2)*tt25+2.0E+0*Dm(2,1)*k(2,1)*tt8*tt9*tt15+tt61)
tt65 = 2.0E+0*k(1,1)*tt1*Dm(2,1)*tt6*tt20
tt66 = 2.0E+0*tt8*Dm(2,2)*tt9*k(3,1)*tt27
tt67 = Dm(2,1)*tt27+Dm(2,2)*tt20
tt68 = 5.0E-1*area(1,1)*(5.0E-1*tt32*tt67*k(4,1)+tt66+2.0E+0*tt1*&
&tt6*k(2,1)*Dm(2,2)*tt27+2.0E+0*Dm(2,1)*k(2,1)*tt8*tt9*tt20+tt65)
tt69 = 5.0E-1*area(1,1)*(5.0E-1*tt35*tt37*k(4,1)+2.0E+0*tt10*tt25&
&*k(3,1)*tt27+2.0E+0*tt1*k(2,1)*tt8*tt15*tt27+2.0E+0*tt1*k(2,1)*tt&
&8*tt25*tt20+2.0E+0*k(1,1)*tt2*tt15*tt20)
tt70 = 5.0E-1*area(1,1)*(5.0E-1*tt43*tt35*k(4,1)+tt47+2.0E+0*Dm(1&
&,1)*tt6*k(2,1)*tt8*tt25+2.0E+0*Dm(1,2)*tt1*k(2,1)*tt9*tt15+tt46)
tt71 = 5.0E-1*area(1,1)*(tt44+5.0E-1*tt48*tt35*k(4,1)+tt42+tt41+t&
&t40+tt39+2.0E+0*Dm(1,2)*tt8*tt26*k(3,1)+2.0E+0*Dm(1,1)*k(1,1)*tt1&
&*tt16+2.0E+0*Dm(1,1)*k(2,1)*tt8*tt15*tt25+2.0E+0*Dm(1,2)*tt1*k(2,&
&1)*tt15*tt25)
tt72 = 2.0E+0*Dm(1,1)*k(1,1)*tt1*tt15*tt20
tt73 = 2.0E+0*Dm(1,2)*tt8*tt25*k(3,1)*tt27
tt74 = 5.0E-1*area(1,1)*(5.0E-1*tt35*tt52*k(4,1)+tt73+2.0E+0*Dm(1&
&,2)*tt1*k(2,1)*tt15*tt27+2.0E+0*Dm(1,1)*k(2,1)*tt8*tt25*tt20+tt72&
&)
tt75 = 5.0E-1*area(1,1)*(5.0E-1*tt58*tt35*k(4,1)+tt62+2.0E+0*Dm(2&
&,1)*tt6*k(2,1)*tt8*tt25+2.0E+0*tt1*k(2,1)*Dm(2,2)*tt9*tt15+tt61)
tt76 = 5.0E-1*area(1,1)*(tt59+5.0E-1*tt35*tt63*k(4,1)+tt57+tt56+t&
&t55+tt54+2.0E+0*tt8*Dm(2,2)*tt26*k(3,1)+2.0E+0*k(1,1)*tt1*Dm(2,1)&
&*tt16+2.0E+0*tt1*k(2,1)*Dm(2,2)*tt15*tt25+2.0E+0*Dm(2,1)*k(2,1)*t&
&t8*tt15*tt25)
tt77 = 2.0E+0*k(1,1)*tt1*Dm(2,1)*tt15*tt20
tt78 = 2.0E+0*tt8*Dm(2,2)*tt25*k(3,1)*tt27
tt79 = 5.0E-1*area(1,1)*(5.0E-1*tt35*tt67*k(4,1)+tt78+2.0E+0*tt1*&
&k(2,1)*Dm(2,2)*tt15*tt27+2.0E+0*Dm(2,1)*k(2,1)*tt8*tt25*tt20+tt77&
&)
tt80 = 5.0E-1*area(1,1)*(5.0E-1*tt43*tt37*k(4,1)+tt51+2.0E+0*Dm(1&
&,1)*tt6*k(2,1)*tt8*tt27+2.0E+0*Dm(1,2)*tt1*k(2,1)*tt9*tt20+tt50)
tt81 = 5.0E-1*area(1,1)*(5.0E-1*tt48*tt37*k(4,1)+tt73+2.0E+0*Dm(1&
&,1)*k(2,1)*tt8*tt15*tt27+2.0E+0*Dm(1,2)*tt1*k(2,1)*tt25*tt20+tt72&
&)
tt82 = 5.0E-1*area(1,1)*(tt44+5.0E-1*tt52*tt37*k(4,1)+tt42+tt41+2&
&.0E+0*Dm(1,2)*tt8*k(3,1)*tt28+tt40+tt39+2.0E+0*Dm(1,1)*k(1,1)*tt1&
&*tt21+2.0E+0*Dm(1,1)*k(2,1)*tt8*tt20*tt27+2.0E+0*Dm(1,2)*tt1*k(2,&
&1)*tt20*tt27)
tt83 = 5.0E-1*area(1,1)*(5.0E-1*tt58*tt37*k(4,1)+tt66+2.0E+0*Dm(2&
&,1)*tt6*k(2,1)*tt8*tt27+2.0E+0*tt1*k(2,1)*Dm(2,2)*tt9*tt20+tt65)
tt84 = 5.0E-1*area(1,1)*(5.0E-1*tt63*tt37*k(4,1)+tt78+2.0E+0*Dm(2&
&,1)*k(2,1)*tt8*tt15*tt27+2.0E+0*tt1*k(2,1)*Dm(2,2)*tt25*tt20+tt77&
&)
tt85 = 5.0E-1*area(1,1)*(tt59+5.0E-1*tt37*tt67*k(4,1)+tt57+tt56+2&
&.0E+0*tt8*Dm(2,2)*k(3,1)*tt28+tt55+tt54+2.0E+0*k(1,1)*tt1*Dm(2,1)&
&*tt21+2.0E+0*tt1*k(2,1)*Dm(2,2)*tt20*tt27+2.0E+0*Dm(2,1)*k(2,1)*t&
&t8*tt20*tt27)
tt86 = Dm(1,1)**2
tt87 = Dm(1,2)**2
tt88 = 1.0E+0*tt86*k(1,1)*tt22
tt89 = 1.0E+0*tt87*k(2,1)*tt22
tt90 = 1.0E+0*tt86*k(2,1)*tt29
tt91 = 1.0E+0*tt87*k(3,1)*tt29
tt92 = 1.0E+0*Dm(1,1)*Dm(1,2)*tt33*k(4,1)
tt93 = 5.0E-1*area(1,1)*(5.0E-1*tt43*tt48*k(4,1)+2.0E+0*tt87*tt9*&
&tt25*k(3,1)+2.0E+0*Dm(1,1)*Dm(1,2)*tt6*k(2,1)*tt25+2.0E+0*Dm(1,1)&
&*Dm(1,2)*k(2,1)*tt9*tt15+2.0E+0*tt86*k(1,1)*tt6*tt15)
tt94 = 5.0E-1*area(1,1)*(5.0E-1*tt43*tt52*k(4,1)+2.0E+0*tt87*tt9*&
&k(3,1)*tt27+2.0E+0*Dm(1,1)*Dm(1,2)*tt6*k(2,1)*tt27+2.0E+0*Dm(1,1)&
&*Dm(1,2)*k(2,1)*tt9*tt20+2.0E+0*tt86*k(1,1)*tt6*tt20)
tt95 = 1.0E+0*Dm(1,1)*k(1,1)*Dm(2,1)*tt22
tt96 = 1.0E+0*Dm(1,2)*k(2,1)*Dm(2,2)*tt22
tt97 = 1.0E+0*Dm(1,1)*Dm(2,1)*k(2,1)*tt29
tt98 = 1.0E+0*Dm(1,2)*Dm(2,2)*k(3,1)*tt29
tt99 = 5.0E-1*(Dm(1,1)*Dm(2,2)+Dm(1,2)*Dm(2,1))*tt33*k(4,1)
tt100 = 5.0E-1*area(1,1)*(tt99+5.0E-1*tt43*tt58*k(4,1)+tt98+tt97+&
&tt96+tt95+2.0E+0*Dm(1,2)*Dm(2,2)*tt11*k(3,1)+2.0E+0*Dm(1,1)*tt6*k&
&(2,1)*Dm(2,2)*tt9+2.0E+0*Dm(1,2)*Dm(2,1)*tt6*k(2,1)*tt9+2.0E+0*Dm&
&(1,1)*k(1,1)*Dm(2,1)*tt7)
tt101 = 2.0E+0*Dm(1,1)*k(1,1)*Dm(2,1)*tt6*tt15
tt102 = 2.0E+0*Dm(1,2)*Dm(2,2)*tt9*tt25*k(3,1)
tt103 = 5.0E-1*area(1,1)*(5.0E-1*tt43*tt63*k(4,1)+tt102+2.0E+0*Dm&
&(1,1)*tt6*k(2,1)*Dm(2,2)*tt25+2.0E+0*Dm(1,2)*Dm(2,1)*k(2,1)*tt9*t&
&t15+tt101)
tt104 = 2.0E+0*Dm(1,1)*k(1,1)*Dm(2,1)*tt6*tt20
tt105 = 2.0E+0*Dm(1,2)*Dm(2,2)*tt9*k(3,1)*tt27
tt106 = 5.0E-1*area(1,1)*(5.0E-1*tt43*tt67*k(4,1)+tt105+2.0E+0*Dm&
&(1,1)*tt6*k(2,1)*Dm(2,2)*tt27+2.0E+0*Dm(1,2)*Dm(2,1)*k(2,1)*tt9*t&
&t20+tt104)
tt107 = 5.0E-1*area(1,1)*(5.0E-1*tt48*tt52*k(4,1)+2.0E+0*tt87*tt2&
&5*k(3,1)*tt27+2.0E+0*Dm(1,1)*Dm(1,2)*k(2,1)*tt15*tt27+2.0E+0*Dm(1&
&,1)*Dm(1,2)*k(2,1)*tt25*tt20+2.0E+0*tt86*k(1,1)*tt15*tt20)
tt108 = 5.0E-1*area(1,1)*(5.0E-1*tt58*tt48*k(4,1)+tt102+2.0E+0*Dm&
&(1,2)*Dm(2,1)*tt6*k(2,1)*tt25+2.0E+0*Dm(1,1)*k(2,1)*Dm(2,2)*tt9*t&
&t15+tt101)
tt109 = 5.0E-1*area(1,1)*(tt99+5.0E-1*tt48*tt63*k(4,1)+tt98+tt97+&
&tt96+tt95+2.0E+0*Dm(1,2)*Dm(2,2)*tt26*k(3,1)+2.0E+0*Dm(1,1)*k(1,1&
&)*Dm(2,1)*tt16+2.0E+0*Dm(1,1)*k(2,1)*Dm(2,2)*tt15*tt25+2.0E+0*Dm(&
&1,2)*Dm(2,1)*k(2,1)*tt15*tt25)
tt110 = 2.0E+0*Dm(1,1)*k(1,1)*Dm(2,1)*tt15*tt20
tt111 = 2.0E+0*Dm(1,2)*Dm(2,2)*tt25*k(3,1)*tt27
tt112 = 5.0E-1*area(1,1)*(5.0E-1*tt48*tt67*k(4,1)+tt111+2.0E+0*Dm&
&(1,1)*k(2,1)*Dm(2,2)*tt15*tt27+2.0E+0*Dm(1,2)*Dm(2,1)*k(2,1)*tt25&
&*tt20+tt110)
tt113 = 5.0E-1*area(1,1)*(5.0E-1*tt58*tt52*k(4,1)+tt105+2.0E+0*Dm&
&(1,2)*Dm(2,1)*tt6*k(2,1)*tt27+2.0E+0*Dm(1,1)*k(2,1)*Dm(2,2)*tt9*t&
&t20+tt104)
tt114 = 5.0E-1*area(1,1)*(5.0E-1*tt63*tt52*k(4,1)+tt111+2.0E+0*Dm&
&(1,2)*Dm(2,1)*k(2,1)*tt15*tt27+2.0E+0*Dm(1,1)*k(2,1)*Dm(2,2)*tt25&
&*tt20+tt110)
tt115 = 5.0E-1*area(1,1)*(tt99+5.0E-1*tt52*tt67*k(4,1)+tt98+tt97+&
&2.0E+0*Dm(1,2)*Dm(2,2)*k(3,1)*tt28+tt96+tt95+2.0E+0*Dm(1,1)*k(1,1&
&)*Dm(2,1)*tt21+2.0E+0*Dm(1,1)*k(2,1)*Dm(2,2)*tt20*tt27+2.0E+0*Dm(&
&1,2)*Dm(2,1)*k(2,1)*tt20*tt27)
tt116 = Dm(2,1)**2
tt117 = Dm(2,2)**2
tt118 = 1.0E+0*k(1,1)*tt116*tt22
tt119 = 1.0E+0*k(2,1)*tt117*tt22
tt120 = 1.0E+0*tt116*k(2,1)*tt29
tt121 = 1.0E+0*tt117*k(3,1)*tt29
tt122 = 1.0E+0*Dm(2,1)*Dm(2,2)*tt33*k(4,1)
tt123 = 5.0E-1*area(1,1)*(5.0E-1*tt58*tt63*k(4,1)+2.0E+0*tt117*tt&
&9*tt25*k(3,1)+2.0E+0*Dm(2,1)*tt6*k(2,1)*Dm(2,2)*tt25+2.0E+0*Dm(2,&
&1)*k(2,1)*Dm(2,2)*tt9*tt15+2.0E+0*k(1,1)*tt116*tt6*tt15)
tt124 = 5.0E-1*area(1,1)*(5.0E-1*tt58*tt67*k(4,1)+2.0E+0*tt117*tt&
&9*k(3,1)*tt27+2.0E+0*Dm(2,1)*tt6*k(2,1)*Dm(2,2)*tt27+2.0E+0*Dm(2,&
&1)*k(2,1)*Dm(2,2)*tt9*tt20+2.0E+0*k(1,1)*tt116*tt6*tt20)
tt125 = 5.0E-1*area(1,1)*(5.0E-1*tt63*tt67*k(4,1)+2.0E+0*tt117*tt&
&25*k(3,1)*tt27+2.0E+0*Dm(2,1)*k(2,1)*Dm(2,2)*tt15*tt27+2.0E+0*Dm(&
&2,1)*k(2,1)*Dm(2,2)*tt25*tt20+2.0E+0*k(1,1)*tt116*tt15*tt20)
hes(1,1) = 5.0E-1*area(1,1)*(tt34+5.0E-1*tt32**2*k(4,1)+tt31+tt30&
&+tt24+tt23+2.0E+0*tt10*tt11*k(3,1)+4.0E+0*tt1*tt6*k(2,1)*tt8*tt9+&
&2.0E+0*k(1,1)*tt2*tt7)
hes(1,2) = tt36
hes(1,3) = tt38
hes(1,4) = tt45
hes(1,5) = tt49
hes(1,6) = tt53
hes(1,7) = tt60
hes(1,8) = tt64
hes(1,9) = tt68
hes(2,1) = tt36
hes(2,2) = 5.0E-1*area(1,1)*(tt34+5.0E-1*tt35**2*k(4,1)+tt31+tt30&
&+tt24+tt23+2.0E+0*tt10*tt26*k(3,1)+2.0E+0*k(1,1)*tt2*tt16+4.0E+0*&
&tt1*k(2,1)*tt8*tt15*tt25)
hes(2,3) = tt69
hes(2,4) = tt70
hes(2,5) = tt71
hes(2,6) = tt74
hes(2,7) = tt75
hes(2,8) = tt76
hes(2,9) = tt79
hes(3,1) = tt38
hes(3,2) = tt69
hes(3,3) = 5.0E-1*area(1,1)*(5.0E-1*tt37**2*k(4,1)+tt34+tt31+tt30&
&+2.0E+0*tt10*k(3,1)*tt28+tt24+tt23+2.0E+0*k(1,1)*tt2*tt21+4.0E+0*&
&tt1*k(2,1)*tt8*tt20*tt27)
hes(3,4) = tt80
hes(3,5) = tt81
hes(3,6) = tt82
hes(3,7) = tt83
hes(3,8) = tt84
hes(3,9) = tt85
hes(4,1) = tt45
hes(4,2) = tt70
hes(4,3) = tt80
hes(4,4) = 5.0E-1*area(1,1)*(tt92+5.0E-1*tt43**2*k(4,1)+tt91+tt90&
&+tt89+tt88+2.0E+0*tt87*tt11*k(3,1)+4.0E+0*Dm(1,1)*Dm(1,2)*tt6*k(2&
&,1)*tt9+2.0E+0*tt86*k(1,1)*tt7)
hes(4,5) = tt93
hes(4,6) = tt94
hes(4,7) = tt100
hes(4,8) = tt103
hes(4,9) = tt106
hes(5,1) = tt49
hes(5,2) = tt71
hes(5,3) = tt81
hes(5,4) = tt93
hes(5,5) = 5.0E-1*area(1,1)*(tt92+5.0E-1*tt48**2*k(4,1)+tt91+tt90&
&+tt89+tt88+2.0E+0*tt87*tt26*k(3,1)+2.0E+0*tt86*k(1,1)*tt16+4.0E+0&
&*Dm(1,1)*Dm(1,2)*k(2,1)*tt15*tt25)
hes(5,6) = tt107
hes(5,7) = tt108
hes(5,8) = tt109
hes(5,9) = tt112
hes(6,1) = tt53
hes(6,2) = tt74
hes(6,3) = tt82
hes(6,4) = tt94
hes(6,5) = tt107
hes(6,6) = 5.0E-1*area(1,1)*(5.0E-1*tt52**2*k(4,1)+tt92+tt91+tt90&
&+2.0E+0*tt87*k(3,1)*tt28+tt89+tt88+2.0E+0*tt86*k(1,1)*tt21+4.0E+0&
&*Dm(1,1)*Dm(1,2)*k(2,1)*tt20*tt27)
hes(6,7) = tt113
hes(6,8) = tt114
hes(6,9) = tt115
hes(7,1) = tt60
hes(7,2) = tt75
hes(7,3) = tt83
hes(7,4) = tt100
hes(7,5) = tt108
hes(7,6) = tt113
hes(7,7) = 5.0E-1*area(1,1)*(tt122+5.0E-1*tt58**2*k(4,1)+tt121+tt&
&120+tt119+tt118+2.0E+0*tt117*tt11*k(3,1)+4.0E+0*Dm(2,1)*tt6*k(2,1&
&)*Dm(2,2)*tt9+2.0E+0*k(1,1)*tt116*tt7)
hes(7,8) = tt123
hes(7,9) = tt124
hes(8,1) = tt64
hes(8,2) = tt76
hes(8,3) = tt84
hes(8,4) = tt103
hes(8,5) = tt109
hes(8,6) = tt114
hes(8,7) = tt123
hes(8,8) = 5.0E-1*area(1,1)*(tt122+5.0E-1*tt63**2*k(4,1)+tt121+tt&
&120+tt119+tt118+2.0E+0*tt117*tt26*k(3,1)+2.0E+0*k(1,1)*tt116*tt16&
&+4.0E+0*Dm(2,1)*k(2,1)*Dm(2,2)*tt15*tt25)
hes(8,9) = tt125
hes(9,1) = tt68
hes(9,2) = tt79
hes(9,3) = tt85
hes(9,4) = tt106
hes(9,5) = tt112
hes(9,6) = tt115
hes(9,7) = tt124
hes(9,8) = tt125
hes(9,9) = 5.0E-1*area(1,1)*(5.0E-1*tt67**2*k(4,1)+tt122+tt121+tt&
&120+2.0E+0*tt117*k(3,1)*tt28+tt119+tt118+2.0E+0*k(1,1)*tt116*tt21&
&+4.0E+0*Dm(2,1)*k(2,1)*Dm(2,2)*tt20*tt27)
END 
SUBROUTINE calc_edge_length(val, X) 
IMPLICIT NONE 
REAL(KIND=8) val(1, 1) 
REAL(KIND=8) X(3, 2) 
val(1,1) = sqrt((X(3,1)-X(3,2))**2+(X(2,1)-X(2,2))**2+(X(1,1)-X(1&
&,2))**2)
END 
SUBROUTINE calc_edge_length_jac(jac, X) 
IMPLICIT NONE 
REAL(KIND=8) jac(1, 6) 
REAL(KIND=8) X(3, 2) 
REAL(KIND=8)  tt1 
REAL(KIND=8)  tt2 
REAL(KIND=8)  tt3 
REAL(KIND=8)  tt4 
tt1 = X(1,1)-X(1,2)
tt2 = X(2,1)-X(2,2)
tt3 = X(3,1)-X(3,2)
tt4 = 1/sqrt(tt3**2+tt2**2+tt1**2)
jac(1,1) = tt1*tt4
jac(1,2) = tt2*tt4
jac(1,3) = tt4*tt3
jac(1,4) = -tt1*tt4
jac(1,5) = -tt2*tt4
jac(1,6) = -tt4*tt3
END 
SUBROUTINE calc_edge_length_hes(hes, X) 
IMPLICIT NONE 
REAL(KIND=8) hes(6, 6) 
REAL(KIND=8) X(3, 2) 
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
tt1 = X(1,1)-X(1,2)
tt2 = tt1**2
tt3 = X(2,1)-X(2,2)
tt4 = tt3**2
tt5 = X(3,1)-X(3,2)
tt6 = tt5**2
tt7 = sqrt(tt6+tt4+tt2)
tt8 = 1/tt7**3
tt9 = 1/tt7
tt10 = tt9-tt2*tt8
tt11 = -tt1*tt3*tt8
tt12 = -tt1*tt8*tt5
tt13 = -tt9
tt14 = tt13+tt2*tt8
tt15 = tt1*tt3*tt8
tt16 = tt1*tt8*tt5
tt17 = tt9-tt4*tt8
tt18 = -tt3*tt8*tt5
tt19 = tt13+tt4*tt8
tt20 = tt3*tt8*tt5
tt21 = tt9-tt8*tt6
tt22 = tt8*tt6+tt13
hes(1,1) = tt10
hes(1,2) = tt11
hes(1,3) = tt12
hes(1,4) = tt14
hes(1,5) = tt15
hes(1,6) = tt16
hes(2,1) = tt11
hes(2,2) = tt17
hes(2,3) = tt18
hes(2,4) = tt15
hes(2,5) = tt19
hes(2,6) = tt20
hes(3,1) = tt12
hes(3,2) = tt18
hes(3,3) = tt21
hes(3,4) = tt16
hes(3,5) = tt20
hes(3,6) = tt22
hes(4,1) = tt14
hes(4,2) = tt15
hes(4,3) = tt16
hes(4,4) = tt10
hes(4,5) = tt11
hes(4,6) = tt12
hes(5,1) = tt15
hes(5,2) = tt19
hes(5,3) = tt20
hes(5,4) = tt11
hes(5,5) = tt17
hes(5,6) = tt18
hes(6,1) = tt16
hes(6,2) = tt20
hes(6,3) = tt22
hes(6,4) = tt12
hes(6,5) = tt18
hes(6,6) = tt21
END 
SUBROUTINE mass_spring(val, X, d) 
IMPLICIT NONE 
REAL(KIND=8) val(1, 1) 
REAL(KIND=8) X(3, 2) 
REAL(KIND=8) d(1, 1) 
val(1,1) = d(1,1)*(sqrt((X(3,1)-X(3,2))**2+(X(2,1)-X(2,2))**2+(X(&
&1,1)-X(1,2))**2)/d(1,1)-1)**2
END 
SUBROUTINE mass_spring_jac(jac, X, d) 
IMPLICIT NONE 
REAL(KIND=8) jac(1, 6) 
REAL(KIND=8) X(3, 2) 
REAL(KIND=8) d(1, 1) 
REAL(KIND=8)  tt1 
REAL(KIND=8)  tt2 
REAL(KIND=8)  tt3 
REAL(KIND=8)  tt4 
REAL(KIND=8)  tt5 
REAL(KIND=8)  tt6 
tt1 = X(1,1)-X(1,2)
tt2 = X(2,1)-X(2,2)
tt3 = X(3,1)-X(3,2)
tt4 = sqrt(tt3**2+tt2**2+tt1**2)
tt5 = tt4/d(1,1)-1
tt6 = 1/tt4
jac(1,1) = 2*tt1*tt5*tt6
jac(1,2) = 2*tt2*tt5*tt6
jac(1,3) = 2*tt5*tt6*tt3
jac(1,4) = -2*tt1*tt5*tt6
jac(1,5) = -2*tt2*tt5*tt6
jac(1,6) = -2*tt5*tt6*tt3
END 
SUBROUTINE mass_spring_hes(hes, X, d) 
IMPLICIT NONE 
REAL(KIND=8) hes(6, 6) 
REAL(KIND=8) X(3, 2) 
REAL(KIND=8) d(1, 1) 
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
tt1 = X(1,1)-X(1,2)
tt2 = tt1**2
tt3 = 1/d(1,1)
tt4 = X(2,1)-X(2,2)
tt5 = tt4**2
tt6 = X(3,1)-X(3,2)
tt7 = tt6**2
tt8 = tt7+tt5+tt2
tt9 = sqrt(tt8)
tt10 = tt3*tt9-1
tt11 = 1/tt9**3
tt12 = 1/tt8
tt13 = 1/tt9
tt14 = 2*tt10*tt13
tt15 = tt14+2*tt3*tt2*tt12-2*tt2*tt10*tt11
tt16 = 2*tt3*tt1*tt4*tt12-2*tt1*tt4*tt10*tt11
tt17 = 2*tt3*tt1*tt12*tt6-2*tt1*tt10*tt11*tt6
tt18 = -2*tt10*tt13
tt19 = tt18-2*tt3*tt2*tt12+2*tt2*tt10*tt11
tt20 = 2*tt1*tt4*tt10*tt11-2*tt3*tt1*tt4*tt12
tt21 = 2*tt1*tt10*tt11*tt6-2*tt3*tt1*tt12*tt6
tt22 = tt14+2*tt3*tt5*tt12-2*tt5*tt10*tt11
tt23 = 2*tt3*tt4*tt12*tt6-2*tt4*tt10*tt11*tt6
tt24 = tt18-2*tt3*tt5*tt12+2*tt5*tt10*tt11
tt25 = 2*tt4*tt10*tt11*tt6-2*tt3*tt4*tt12*tt6
tt26 = 2*tt3*tt12*tt7-2*tt10*tt11*tt7+tt14
tt27 = (-2*tt3*tt12*tt7)+2*tt10*tt11*tt7+tt18
hes(1,1) = tt15
hes(1,2) = tt16
hes(1,3) = tt17
hes(1,4) = tt19
hes(1,5) = tt20
hes(1,6) = tt21
hes(2,1) = tt16
hes(2,2) = tt22
hes(2,3) = tt23
hes(2,4) = tt20
hes(2,5) = tt24
hes(2,6) = tt25
hes(3,1) = tt17
hes(3,2) = tt23
hes(3,3) = tt26
hes(3,4) = tt21
hes(3,5) = tt25
hes(3,6) = tt27
hes(4,1) = tt19
hes(4,2) = tt20
hes(4,3) = tt21
hes(4,4) = tt15
hes(4,5) = tt16
hes(4,6) = tt17
hes(5,1) = tt20
hes(5,2) = tt24
hes(5,3) = tt25
hes(5,4) = tt16
hes(5,5) = tt22
hes(5,6) = tt23
hes(6,1) = tt21
hes(6,2) = tt25
hes(6,3) = tt27
hes(6,4) = tt17
hes(6,5) = tt23
hes(6,6) = tt26
END 
SUBROUTINE poly_spring(val, X, d) 
IMPLICIT NONE 
REAL(KIND=8) val(1, 1) 
REAL(KIND=8) X(3, 2) 
REAL(KIND=8) d(1, 1) 
val(1,1) = ((X(3,1)-X(3,2))**2+(X(2,1)-X(2,2))**2+(X(1,1)-X(1,2))&
&**2-d(1,1)**2)**2
END 
SUBROUTINE poly_spring_jac(jac, X, d) 
IMPLICIT NONE 
REAL(KIND=8) jac(1, 6) 
REAL(KIND=8) X(3, 2) 
REAL(KIND=8) d(1, 1) 
REAL(KIND=8)  tt1 
REAL(KIND=8)  tt2 
REAL(KIND=8)  tt3 
REAL(KIND=8)  tt4 
tt1 = X(1,1)-X(1,2)
tt2 = X(2,1)-X(2,2)
tt3 = X(3,1)-X(3,2)
tt4 = tt3**2+tt2**2+tt1**2-d(1,1)**2
jac(1,1) = 4*tt1*tt4
jac(1,2) = 4*tt2*tt4
jac(1,3) = 4*tt4*tt3
jac(1,4) = -4*tt1*tt4
jac(1,5) = -4*tt2*tt4
jac(1,6) = -4*tt4*tt3
END 
SUBROUTINE poly_spring_hes(hes, X, d) 
IMPLICIT NONE 
REAL(KIND=8) hes(6, 6) 
REAL(KIND=8) X(3, 2) 
REAL(KIND=8) d(1, 1) 
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
tt1 = X(1,1)-X(1,2)
tt2 = tt1**2
tt3 = X(2,1)-X(2,2)
tt4 = tt3**2
tt5 = X(3,1)-X(3,2)
tt6 = tt5**2
tt7 = tt6+tt4+tt2-d(1,1)**2
tt8 = 4*tt7
tt9 = tt8+8*tt2
tt10 = 8*tt1*tt3
tt11 = 8*tt1*tt5
tt12 = -4*tt7
tt13 = tt12-8*tt2
tt14 = -8*tt1*tt3
tt15 = -8*tt1*tt5
tt16 = tt8+8*tt4
tt17 = 8*tt3*tt5
tt18 = tt12-8*tt4
tt19 = -8*tt3*tt5
tt20 = tt8+8*tt6
tt21 = tt12-8*tt6
hes(1,1) = tt9
hes(1,2) = tt10
hes(1,3) = tt11
hes(1,4) = tt13
hes(1,5) = tt14
hes(1,6) = tt15
hes(2,1) = tt10
hes(2,2) = tt16
hes(2,3) = tt17
hes(2,4) = tt14
hes(2,5) = tt18
hes(2,6) = tt19
hes(3,1) = tt11
hes(3,2) = tt17
hes(3,3) = tt20
hes(3,4) = tt15
hes(3,5) = tt19
hes(3,6) = tt21
hes(4,1) = tt13
hes(4,2) = tt14
hes(4,3) = tt15
hes(4,4) = tt9
hes(4,5) = tt10
hes(4,6) = tt11
hes(5,1) = tt14
hes(5,2) = tt18
hes(5,3) = tt19
hes(5,4) = tt10
hes(5,5) = tt16
hes(5,6) = tt17
hes(6,1) = tt15
hes(6,2) = tt19
hes(6,3) = tt21
hes(6,4) = tt11
hes(6,5) = tt17
hes(6,6) = tt20
END 
SUBROUTINE calc_dih_angle(val, X) 
IMPLICIT NONE 
REAL(KIND=8) val(1, 1) 
REAL(KIND=8) X(3, 4) 
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
tt1 = X(1,2)**2
tt2 = -2*X(1,2)*X(1,3)
tt3 = X(1,3)**2
tt4 = X(2,2)**2
tt5 = -2*X(2,2)*X(2,3)
tt6 = X(2,3)**2
tt7 = X(3,2)**2
tt8 = X(3,3)**2
tt9 = X(1,3)**3
tt10 = -tt1
tt11 = -tt3
tt12 = (tt11+2*X(1,2)*X(1,3)+tt10)*X(1,4)
tt13 = -X(1,4)
tt14 = tt13+X(1,3)
tt15 = X(2,2)**3
tt16 = X(1,2)**3
tt17 = -tt16
tt18 = tt3+tt2+tt1
tt19 = tt18*X(1,4)
tt20 = -X(1,2)
tt21 = -2*X(1,3)
tt22 = 3*X(1,4)
tt23 = 2*X(1,2)
tt24 = -3*X(1,4)
tt25 = X(1,4)+tt20
tt26 = X(2,3)**3
tt27 = -tt9
tt28 = -X(1,3)
tt29 = tt28+X(1,2)
tt30 = -2*X(1,2)
tt31 = 2*X(1,3)
tt32 = tt31+tt30
tt33 = X(1,4)+tt28
tt34 = X(1,1)*tt1
tt35 = -2*X(1,4)
tt36 = tt13+X(1,1)
tt37 = tt33*X(2,1)
tt38 = -2*X(1,1)
tt39 = 2*X(1,4)
tt40 = tt39+tt38
tt41 = -X(1,1)*tt1
tt42 = 2*X(1,1)*X(1,2)
tt43 = -X(1,1)
tt44 = X(1,3)+tt43
tt45 = 2*X(1,1)
tt46 = tt29*X(2,4)+tt25*X(2,3)+tt14*X(2,2)
tt47 = tt36*X(2,3)
tt48 = tt13+X(1,2)
tt49 = X(1,4)+tt43
tt50 = tt35+tt45
tt51 = tt48*X(2,1)
tt52 = tt49*X(2,2)
tt53 = -2*X(1,1)*X(1,2)
tt54 = tt20+X(1,1)
tt55 = 3*X(1,1)
tt56 = -3*X(1,1)
tt57 = X(1,3)+tt20
tt58 = tt28+X(1,1)
tt59 = X(1,2)+tt43
tt60 = tt21+tt23
tt61 = tt57*X(2,1)
tt62 = tt59*X(2,3)+tt58*X(2,2)+tt61
tt63 = X(1,2)*X(1,3)
tt64 = tt57*X(1,4)
tt65 = -X(1,1)*X(1,3)
tt66 = tt58*X(1,4)
tt67 = tt29*X(1,4)
tt68 = X(1,1)*X(1,2)
tt69 = (tt30+X(1,1))*X(1,3)
tt70 = (X(1,3)+X(1,2)+tt38)*X(1,4)
tt71 = -X(1,1)*X(1,2)
tt72 = tt54*X(1,4)
tt73 = (X(1,2)+X(1,1))*X(1,3)
tt74 = tt59*X(1,3)
tt75 = X(2,2)*X(2,3)
tt76 = -tt6
tt77 = -X(2,2)
tt78 = -X(2,3)
tt79 = -tt4
tt80 = X(2,1)*X(2,2)
tt81 = -X(2,1)*X(2,2)
val(1,1) = atan2(((tt62*tt8+((tt30+tt45)*X(2,3)+(tt31+tt38)*X(2,2&
&)+tt60*X(2,1))*X(3,2)*X(3,3)+tt62*tt7+tt59*tt26+((tt28+tt30+tt55)&
&*X(2,2)+tt61)*tt6+((tt31+X(1,2)+tt56)*tt4+tt60*X(2,1)*X(2,2)+tt59&
&*tt3+(tt42-2*tt1)*X(1,3)+tt16+tt41)*X(2,3)+tt58*tt15+tt57*X(2,1)*&
&tt4+(tt27+(tt23+X(1,1))*tt3+(tt10+tt53)*X(1,3)+tt34)*X(2,2)+(tt9-&
&3*X(1,2)*tt3+3*tt1*X(1,3)+tt17)*X(2,1))*X(3,4)+(tt54*X(2,4)+tt52+&
&tt51)*X(3,3)**3+(((X(1,3)+tt23+tt56)*X(2,4)+tt47+tt50*X(2,2)+(tt2&
&2+tt28+tt30)*X(2,1))*X(3,2)+tt46*X(3,1))*tt8+(((tt21+tt20+tt55)*X&
&(2,4)+tt40*X(2,3)+tt52+(tt24+tt31+X(1,2))*X(2,1))*tt7+(tt32*X(2,4&
&)+(tt35+tt23)*X(2,3)+(tt39+tt21)*X(2,2))*X(3,1)*X(3,2)+(tt54*tt6+&
&(tt23+tt38)*X(2,2)*X(2,3)+tt54*tt4+tt54*tt3+(2*tt1+tt53)*X(1,3)+t&
&t17+tt34)*X(2,4)+(tt52+tt51)*tt6+(tt50*tt4+(tt39+tt30)*X(2,1)*X(2&
&,2))*X(2,3)+tt49*tt15+tt48*X(2,1)*tt4+(tt19-X(1,1)*tt3+2*X(1,1)*X&
&(1,2)*X(1,3)+tt41)*X(2,2)+(tt12+X(1,2)*tt3-2*tt1*X(1,3)+tt16)*X(2&
&,1))*X(3,3)+(tt44*X(2,4)+tt47+tt37)*X(3,2)**3+tt46*X(3,1)*tt7+((t&
&t44*tt6+(tt21+tt45)*X(2,2)*X(2,3)+tt44*tt4+tt9+(tt30+tt43)*tt3+(t&
&t1+tt42)*X(1,3)+tt41)*X(2,4)+tt36*tt26+(tt40*X(2,2)+tt37)*tt6+(tt&
&36*tt4+(tt35+tt31)*X(2,1)*X(2,2)+tt12+X(1,1)*tt3-2*X(1,1)*X(1,2)*&
&X(1,3)+tt34)*X(2,3)+tt33*X(2,1)*tt4+(tt19+tt27+2*X(1,2)*tt3-tt1*X&
&(1,3))*X(2,1))*X(3,2)+((tt29*tt6+tt32*X(2,2)*X(2,3)+tt29*tt4+tt27&
&+3*X(1,2)*tt3-3*tt1*X(1,3)+tt16)*X(2,4)+tt25*tt26+(tt24+X(1,3)+tt&
&23)*X(2,2)*tt6+((tt22+tt21+tt20)*tt4+tt19-X(1,2)*tt3+2*tt1*X(1,3)&
&+tt17)*X(2,3)+tt14*tt15+(tt12+tt9-2*X(1,2)*tt3+tt1*X(1,3))*X(2,2)&
&)*X(3,1))/sqrt(tt8-2*X(3,2)*X(3,3)+tt7+tt6+tt5+tt4+tt3+tt2+tt1),(&
&-(((X(2,2)-X(2,1))*X(2,3)+tt79+tt80+tt74+tt10+tt68)*X(3,3)+(tt76+&
&(X(2,2)+X(2,1))*X(2,3)+tt81+tt11+tt73+tt71)*X(3,2)+(tt6+tt5+tt4+t&
&t3+tt2+tt1)*X(3,1))*X(3,4))-((tt77+X(2,1))*X(2,4)+tt4+tt81+tt72+t&
&t1+tt71)*tt8-(((X(2,3)+X(2,2)-2*X(2,1))*X(2,4)+(X(2,1)-2*X(2,2))*&
&X(2,3)+tt80+tt70+tt69+tt68)*X(3,2)+((tt78+X(2,2))*X(2,4)+tt75+tt7&
&9+tt67+tt63+tt10)*X(3,1))*X(3,3)-((tt78+X(2,1))*X(2,4)+tt6-X(2,1)&
&*X(2,3)+tt66+tt3+tt65)*tt7-((X(2,3)+tt77)*X(2,4)+tt76+tt75+tt64+t&
&t11+tt63)*X(3,1)*X(3,2)-((tt74+tt10+tt68)*X(2,3)+(tt11+tt73+tt71)&
&*X(2,2)+tt18*X(2,1))*X(2,4)-(tt72+tt1+tt71)*tt6-((tt70+tt69+tt68)&
&*X(2,2)+(tt67+tt63+tt10)*X(2,1))*X(2,3)-(tt66+tt3+tt65)*tt4-(tt64&
&+tt11+tt63)*X(2,1)*X(2,2))
END 
SUBROUTINE calc_dih_angle_jac(jac, X) 
IMPLICIT NONE 
REAL(KIND=8) jac(1, 12) 
REAL(KIND=8) X(3, 4) 
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
tt1 = X(1,2)**2
tt2 = -2*X(1,2)*X(1,3)
tt3 = X(1,3)**2
tt4 = X(2,2)**2
tt5 = -2*X(2,2)*X(2,3)
tt6 = X(2,3)**2
tt7 = X(3,2)**2
tt8 = X(3,3)**2
tt9 = tt8-2*X(3,2)*X(3,3)+tt7+tt6+tt5+tt4+tt3+tt2+tt1
tt10 = sqrt(tt9)
tt11 = 1/tt10
tt12 = X(1,2)*X(1,3)
tt13 = -tt3
tt14 = -X(1,2)
tt15 = X(1,3)+tt14
tt16 = tt15*X(1,4)
tt17 = tt16+tt13+tt12
tt18 = -X(1,1)*X(1,3)
tt19 = -X(1,3)
tt20 = tt19+X(1,1)
tt21 = tt20*X(1,4)
tt22 = tt21+tt3+tt18
tt23 = -tt1
tt24 = tt19+X(1,2)
tt25 = tt24*X(1,4)
tt26 = tt25+tt12+tt23
tt27 = X(1,1)*X(1,2)
tt28 = -2*X(1,2)
tt29 = (tt28+X(1,1))*X(1,3)
tt30 = -2*X(1,1)
tt31 = X(1,3)+X(1,2)+tt30
tt32 = tt31*X(1,4)
tt33 = tt32+tt29+tt27
tt34 = -X(1,1)*X(1,2)
tt35 = tt14+X(1,1)
tt36 = tt35*X(1,4)
tt37 = tt36+tt1+tt34
tt38 = tt3+tt2+tt1
tt39 = tt38*X(2,1)
tt40 = (X(1,2)+X(1,1))*X(1,3)
tt41 = tt13+tt40+tt34
tt42 = -X(1,1)
tt43 = X(1,2)+tt42
tt44 = tt43*X(1,3)
tt45 = tt44+tt23+tt27
tt46 = X(2,2)*X(2,3)
tt47 = -tt6
tt48 = -X(2,2)
tt49 = X(2,3)+tt48
tt50 = tt49*X(2,4)+tt47+tt46+tt16+tt13+tt12
tt51 = -X(2,3)
tt52 = tt51+X(2,1)
tt53 = tt52*X(2,4)+tt6-X(2,1)*X(2,3)+tt21+tt3+tt18
tt54 = -tt4
tt55 = tt51+X(2,2)
tt56 = tt55*X(2,4)+tt46+tt54+tt25+tt12+tt23
tt57 = X(2,1)*X(2,2)
tt58 = -2*X(2,2)
tt59 = -2*X(2,1)
tt60 = X(2,3)+X(2,2)+tt59
tt61 = tt60*X(2,4)+(tt58+X(2,1))*X(2,3)+tt57+tt32+tt29+tt27
tt62 = -X(2,1)*X(2,2)
tt63 = tt48+X(2,1)
tt64 = tt63*X(2,4)+tt4+tt62+tt36+tt1+tt34
tt65 = tt6+tt5+tt4+tt3+tt2+tt1
tt66 = tt47+(X(2,2)+X(2,1))*X(2,3)+tt62+tt13+tt40+tt34
tt67 = -X(2,1)
tt68 = X(2,2)+tt67
tt69 = tt68*X(2,3)+tt54+tt57+tt44+tt23+tt27
tt70 = (-(tt69*X(3,3)+tt66*X(3,2)+tt65*X(3,1))*X(3,4))-tt64*tt8-(&
&tt61*X(3,2)+tt56*X(3,1))*X(3,3)-tt53*tt7-tt50*X(3,1)*X(3,2)-(tt45&
&*X(2,3)+tt41*X(2,2)+tt39)*X(2,4)-tt37*tt6-(tt33*X(2,2)+tt26*X(2,1&
&))*X(2,3)-tt22*tt4-tt17*X(2,1)*X(2,2)
tt71 = X(2,3)**3
tt72 = 2*X(1,2)*X(1,3)
tt73 = 2*X(2,2)*X(2,3)
tt74 = -X(2,4)
tt75 = tt74+X(2,3)
tt76 = X(3,2)**3
tt77 = tt13+tt72+tt23
tt78 = tt77*X(2,2)
tt79 = X(2,2)**3
tt80 = -tt79
tt81 = -2*X(2,3)
tt82 = 2*X(2,2)
tt83 = X(2,4)+tt48
tt84 = X(3,3)**3
tt85 = tt38*X(2,2)
tt86 = -tt71
tt87 = 2*X(2,3)
tt88 = tt87+tt58
tt89 = tt1*X(1,3)
tt90 = -2*X(1,2)*tt3
tt91 = X(1,3)**3
tt92 = tt77*X(1,4)
tt93 = (tt92+tt91+tt90+tt89)*X(2,2)
tt94 = -X(1,4)
tt95 = tt94+X(1,3)
tt96 = tt95*tt79
tt97 = X(1,2)**3
tt98 = -tt97
tt99 = 2*tt1*X(1,3)
tt100 = -X(1,2)*tt3
tt101 = tt38*X(1,4)
tt102 = -2*X(1,3)
tt103 = 3*X(1,4)
tt104 = tt103+tt102+tt14
tt105 = tt104*tt4
tt106 = (tt105+tt101+tt100+tt99+tt98)*X(2,3)
tt107 = 2*X(1,2)
tt108 = -3*X(1,4)
tt109 = tt108+X(1,3)+tt107
tt110 = tt109*X(2,2)*tt6
tt111 = X(1,4)+tt14
tt112 = tt111*tt71
tt113 = -tt91
tt114 = 2*X(1,3)
tt115 = tt114+tt28
tt116 = tt24*tt6+tt115*X(2,2)*X(2,3)+tt24*tt4+tt113+3*X(1,2)*tt3-&
&3*tt1*X(1,3)+tt97
tt117 = tt116*X(2,4)
tt118 = -tt1*X(1,3)
tt119 = 2*X(1,2)*tt3
tt120 = (tt101+tt113+tt119+tt118)*X(2,1)
tt121 = X(1,4)+tt19
tt122 = tt121*X(2,1)*tt4
tt123 = X(1,1)*tt1
tt124 = -2*X(1,1)*X(1,2)*X(1,3)
tt125 = X(1,1)*tt3
tt126 = -2*X(1,4)
tt127 = tt126+tt114
tt128 = tt127*X(2,1)*X(2,2)
tt129 = tt94+X(1,1)
tt130 = tt129*tt4
tt131 = (tt130+tt128+tt92+tt125+tt124+tt123)*X(2,3)
tt132 = tt121*X(2,1)
tt133 = 2*X(1,4)
tt134 = tt133+tt30
tt135 = tt134*X(2,2)+tt132
tt136 = tt135*tt6
tt137 = tt129*tt71
tt138 = -X(1,1)*tt1
tt139 = 2*X(1,1)*X(1,2)
tt140 = tt28+tt42
tt141 = X(1,3)+tt42
tt142 = 2*X(1,1)
tt143 = tt102+tt142
tt144 = tt141*tt6+tt143*X(2,2)*X(2,3)+tt141*tt4+tt91+tt140*tt3+(t&
&t1+tt139)*X(1,3)+tt138
tt145 = tt144*X(2,4)
tt146 = tt24*X(2,4)+tt111*X(2,3)+tt95*X(2,2)
tt147 = tt129*X(2,3)
tt148 = tt141*X(2,4)+tt147+tt132
tt149 = -2*tt1*X(1,3)
tt150 = X(1,2)*tt3
tt151 = (tt92+tt150+tt149+tt97)*X(2,1)
tt152 = 2*X(1,1)*X(1,2)*X(1,3)
tt153 = -X(1,1)*tt3
tt154 = (tt101+tt153+tt152+tt138)*X(2,2)
tt155 = tt94+X(1,2)
tt156 = tt155*X(2,1)*tt4
tt157 = X(1,4)+tt42
tt158 = tt157*tt79
tt159 = tt133+tt28
tt160 = tt159*X(2,1)*X(2,2)
tt161 = tt126+tt142
tt162 = tt161*tt4
tt163 = (tt162+tt160)*X(2,3)
tt164 = tt155*X(2,1)
tt165 = tt157*X(2,2)
tt166 = tt165+tt164
tt167 = tt166*tt6
tt168 = -2*X(1,1)*X(1,2)
tt169 = 2*tt1
tt170 = (tt169+tt168)*X(1,3)
tt171 = tt35*tt3
tt172 = tt35*tt4
tt173 = tt107+tt30
tt174 = tt173*X(2,2)*X(2,3)
tt175 = tt35*tt6
tt176 = (tt175+tt174+tt172+tt171+tt170+tt98+tt123)*X(2,4)
tt177 = tt133+tt102
tt178 = tt126+tt107
tt179 = tt115*X(2,4)+tt178*X(2,3)+tt177*X(2,2)
tt180 = tt179*X(3,1)*X(3,2)
tt181 = tt108+tt114+X(1,2)
tt182 = 3*X(1,1)
tt183 = tt102+tt14+tt182
tt184 = tt183*X(2,4)+tt134*X(2,3)+tt165+tt181*X(2,1)
tt185 = tt184*tt7
tt186 = tt103+tt19+tt28
tt187 = -3*X(1,1)
tt188 = X(1,3)+tt107+tt187
tt189 = tt188*X(2,4)+tt147+tt161*X(2,2)+tt186*X(2,1)
tt190 = tt189*X(3,2)+tt146*X(3,1)
tt191 = tt35*X(2,4)+tt165+tt164
tt192 = 3*tt1*X(1,3)
tt193 = -3*X(1,2)*tt3
tt194 = (tt23+tt168)*X(1,3)
tt195 = tt107+X(1,1)
tt196 = tt195*tt3
tt197 = -2*tt1
tt198 = (tt197+tt139)*X(1,3)
tt199 = tt43*tt3
tt200 = tt102+tt107
tt201 = tt200*X(2,1)*X(2,2)
tt202 = tt114+X(1,2)+tt187
tt203 = tt202*tt4
tt204 = tt15*X(2,1)
tt205 = tt19+tt28+tt182
tt206 = tt205*X(2,2)+tt204
tt207 = tt43*X(2,3)
tt208 = tt207+tt20*X(2,2)+tt204
tt209 = tt200*X(2,1)
tt210 = tt114+tt30
tt211 = tt28+tt142
tt212 = tt211*X(2,3)+tt210*X(2,2)+tt209
tt213 = tt208*tt8+tt212*X(3,2)*X(3,3)+tt208*tt7+tt43*tt71+tt206*t&
&t6+(tt203+tt201+tt199+tt198+tt97+tt138)*X(2,3)+tt20*tt79+tt15*X(2&
&,1)*tt4+(tt113+tt196+tt194+tt123)*X(2,2)+(tt91+tt193+tt192+tt98)*&
&X(2,1)
tt214 = tt213*X(3,4)+tt191*tt84+tt190*tt8+(tt185+tt180+tt176+tt16&
&7+tt163+tt158+tt156+tt154+tt151)*X(3,3)+tt148*tt76+tt146*X(3,1)*t&
&t7+(tt145+tt137+tt136+tt131+tt122+tt120)*X(3,2)+(tt117+tt112+tt11&
&0+tt106+tt96+tt93)*X(3,1)
tt215 = 1/(tt214**2/tt9+tt70**2)
tt216 = tt126+X(1,3)+X(1,2)
tt217 = X(2,4)+tt51
tt218 = -2*X(2,4)
tt219 = X(1,4)+X(1,3)+tt28
tt220 = X(1,4)+tt102+X(1,1)
tt221 = tt94+tt107+tt42
tt222 = X(1,3)+tt28+X(1,1)
tt223 = -2*tt3
tt224 = tt115*X(1,4)
tt225 = -3*tt1
tt226 = 4*X(1,2)*X(1,3)
tt227 = tt200*X(1,4)
tt228 = 3*tt1
tt229 = -6*X(1,2)*X(1,3)
tt230 = 3*tt3
tt231 = 2*tt3
tt232 = tt224-2*X(1,1)*X(1,3)+tt139
tt233 = -4*X(1,2)*X(1,3)
tt234 = 2*X(1,1)*X(1,3)
tt235 = X(2,1)*tt4
tt236 = tt74+X(2,1)
tt237 = 2*X(2,4)
tt238 = 6*X(1,2)*X(1,3)
tt239 = -3*tt3
tt240 = -X(2,1)*tt4
tt241 = 2*X(2,1)*X(2,2)
tt242 = X(2,3)+tt67
tt243 = 2*X(2,1)
tt244 = 1/tt10**3
tt245 = tt81+tt82
tt246 = X(1,4)+tt102+X(1,2)
tt247 = tt94+tt114+tt42
tt248 = X(1,4)+tt28+X(1,1)
tt249 = tt102+X(1,2)+X(1,1)
tt250 = -2*tt4
tt251 = tt74+X(2,2)
tt252 = X(2,4)+tt67
tt253 = tt251*X(3,1)
tt254 = -2*X(2,1)*X(2,2)
tt255 = tt24*X(3,1)
jac(1,1) = tt11*tt70*((tt55*tt8+tt88*X(3,2)*X(3,3)+tt55*tt7+tt86+&
&3*X(2,2)*tt6+((-3*tt4)+tt13+tt72+tt23)*X(2,3)+tt79+tt85)*X(3,4)+t&
&t83*tt84+((-3*X(2,4))+X(2,3)+tt82)*X(3,2)*tt8+((3*X(2,4)+tt81+tt4&
&8)*tt7+tt65*X(2,4)-X(2,2)*tt6+2*tt4*X(2,3)+tt80+tt78)*X(3,3)+tt75&
&*tt76+((tt47+tt73+tt54+tt13+tt72+tt23)*X(2,4)+tt71-2*X(2,2)*tt6+(&
&tt4+tt3+tt2+tt1)*X(2,3))*X(3,2))*tt215-tt11*((-(tt24*X(3,3)+tt15*&
&X(3,2))*X(3,4))-tt111*tt8-tt216*X(3,2)*X(3,3)-tt121*tt7-(tt24*X(2&
&,3)+tt15*X(2,2))*X(2,4)-tt111*tt6-tt216*X(2,2)*X(2,3)-tt121*tt4)*&
&tt214*tt215
jac(1,2) = tt11*tt70*((tt15*tt8+tt200*X(3,2)*X(3,3)+tt15*tt7+tt15&
&*tt6+tt200*X(2,2)*X(2,3)+tt15*tt4+tt91+tt193+tt192+tt98)*X(3,4)+t&
&t155*tt84+tt186*X(3,2)*tt8+(tt181*tt7+tt155*tt6+tt159*X(2,2)*X(2,&
&3)+tt155*tt4+tt92+tt150+tt149+tt97)*X(3,3)+tt121*tt76+(tt121*tt6+&
&tt127*X(2,2)*X(2,3)+tt121*tt4+tt101+tt113+tt119+tt118)*X(3,2))*tt&
&215-tt11*((-(tt55*X(3,3)+tt49*X(3,2))*X(3,4))-tt83*tt8-(tt218+X(2&
&,3)+X(2,2))*X(3,2)*X(3,3)-tt217*tt7-tt38*X(2,4)-tt26*X(2,3)-tt17*&
&X(2,2))*tt214*tt215
jac(1,3) = tt11*(tt146*tt8+tt179*X(3,2)*X(3,3)+tt146*tt7+tt117+tt&
&112+tt110+tt106+tt96+tt93)*tt70*tt215-tt11*((-tt65*X(3,4))-tt56*X&
&(3,3)-tt50*X(3,2))*tt214*tt215
jac(1,4) = tt70*(tt11*((tt242*tt8+(tt81+tt243)*X(3,2)*X(3,3)+tt24&
&2*tt7+tt71+(tt58+tt67)*tt6+(tt4+tt241+tt3+(tt142-4*X(1,2))*X(1,3)&
&+tt228+tt168)*X(2,3)+tt240+(tt231+(tt28+tt30)*X(1,3)+tt139)*X(2,2&
&)+(tt239+tt238+tt225)*X(2,1))*X(3,4)+tt236*tt84+((tt237+tt59)*X(3&
&,2)+tt217*X(3,1))*tt8+(tt236*tt7+(tt218+tt87)*X(3,1)*X(3,2)+(tt47&
&+tt73+tt54+tt13+(4*X(1,2)+tt30)*X(1,3)+tt225+tt139)*X(2,4)+X(2,1)&
&*tt6-2*X(2,1)*X(2,2)*X(2,3)+tt235+(tt227+tt234+tt168)*X(2,2)+(tt2&
&24+tt3+tt233+tt228)*X(2,1))*X(3,3)+tt217*X(3,1)*tt7+((tt223+(tt10&
&7+tt142)*X(1,3)+tt168)*X(2,4)+tt232*X(2,3)+(tt227+tt231+tt2)*X(2,&
&1))*X(3,2)+((tt6+tt5+tt4+tt230+tt229+tt228)*X(2,4)+tt86+2*X(2,2)*&
&tt6+(tt54+tt227+tt13+tt226+tt225)*X(2,3)+(tt224+tt223+tt72)*X(2,2&
&))*X(3,1))-(tt200*tt244*tt214)/2.0E+0)*tt215-tt11*((-(tt222*X(3,3&
&)+tt141*X(3,2)+tt200*X(3,1))*X(3,4))-tt221*tt8-(tt220*X(3,2)+tt21&
&9*X(3,1))*X(3,3)-tt95*X(3,1)*X(3,2)-(tt222*X(2,3)+tt141*X(2,2)+tt&
&209)*X(2,4)-tt221*tt6-(tt220*X(2,2)+tt219*X(2,1))*X(2,3)-tt95*X(2&
&,1)*X(2,2))*tt214*tt215
jac(1,5) = tt70*(tt11*((tt20*tt8+tt210*X(3,2)*X(3,3)+tt20*tt7+tt2&
&05*tt6+(2*tt202*X(2,2)+tt209)*X(2,3)+3*tt20*tt4+2*tt15*X(2,1)*X(2&
&,2)+tt113+tt196+tt194+tt123)*X(3,4)+tt157*tt84+(tt161*X(3,2)+tt95&
&*X(3,1))*tt8+(tt157*tt7+tt177*X(3,1)*X(3,2)+(tt173*X(2,3)+2*tt35*&
&X(2,2))*X(2,4)+tt157*tt6+(2*tt161*X(2,2)+tt159*X(2,1))*X(2,3)+3*t&
&t157*tt4+2*tt155*X(2,1)*X(2,2)+tt101+tt153+tt152+tt138)*X(3,3)+tt&
&95*X(3,1)*tt7+((tt143*X(2,3)+2*tt141*X(2,2))*X(2,4)+tt134*tt6+(2*&
&tt129*X(2,2)+tt127*X(2,1))*X(2,3)+2*tt121*X(2,1)*X(2,2))*X(3,2)+(&
&(tt115*X(2,3)+2*tt24*X(2,2))*X(2,4)+tt109*tt6+2*tt104*X(2,2)*X(2,&
&3)+3*tt95*tt4+tt92+tt91+tt90+tt89)*X(3,1))-(tt245*tt244*tt214)/2.&
&0E+0)*tt215-tt11*((-((X(2,3)+tt58+X(2,1))*X(3,3)+tt242*X(3,2)+tt2&
&45*X(3,1))*X(3,4))-(tt74+tt82+tt67)*tt8-((X(2,4)+tt81+X(2,1))*X(3&
&,2)+(X(2,4)+X(2,3)+tt58)*X(3,1))*X(3,3)-tt75*X(3,1)*X(3,2)-tt41*X&
&(2,4)-tt33*X(2,3)-2*tt22*X(2,2)-tt17*X(2,1))*tt214*tt215
jac(1,6) = tt70*(tt11*((tt212*X(3,3)+2*tt208*X(3,2))*X(3,4)+tt189&
&*tt8+(2*tt184*X(3,2)+tt179*X(3,1))*X(3,3)+3*tt148*tt7+2*tt146*X(3&
&,1)*X(3,2)+tt145+tt137+tt136+tt131+tt122+tt120)-((2*X(3,2)-2*X(3,&
&3))*tt244*tt214)/2.0E+0)*tt215-tt11*((-tt66*X(3,4))-tt61*X(3,3)-2&
&*tt53*X(3,2)-tt50*X(3,1))*tt214*tt215
jac(1,7) = tt70*(tt11*((tt63*tt8+(tt82+tt59)*X(3,2)*X(3,3)+tt63*t&
&t7+tt63*tt6+(2*tt4+tt254+2*tt43*X(1,3)+tt197+tt139)*X(2,3)+tt80+t&
&t235+(tt239+2*tt195*X(1,3)+tt23+tt168)*X(2,2)+(tt230+tt229+tt228)&
&*X(2,1))*X(3,4)+(tt252*X(3,2)+tt253)*tt8+((tt218+tt243)*tt7+(tt23&
&7+tt58)*X(3,1)*X(3,2)+(2*tt35*X(1,3)+tt169+tt168)*X(2,4)+tt232*X(&
&2,2)+(tt227+tt72+tt197)*X(2,1))*X(3,3)+tt252*tt76+tt251*X(3,1)*tt&
&7+((tt6+tt5+tt4+tt230+2*tt140*X(1,3)+tt1+tt139)*X(2,4)-X(2,1)*tt6&
&+(tt241+tt227+tt234+tt168)*X(2,3)+tt240+(tt224+tt239+tt226+tt23)*&
&X(2,1))*X(3,2)+((tt47+tt73+tt54+tt239+tt238+tt225)*X(2,4)+X(2,2)*&
&tt6+(tt250+tt224+tt2+tt169)*X(2,3)+tt79+(tt227+tt230+tt233+tt1)*X&
&(2,2))*X(3,1))-(tt115*tt244*tt214)/2.0E+0)*tt215-tt11*((-(tt43*X(&
&3,3)+tt249*X(3,2)+tt115*X(3,1))*X(3,4))-(tt248*X(3,2)+tt155*X(3,1&
&))*X(3,3)-tt247*tt7-tt246*X(3,1)*X(3,2)-(tt207+tt249*X(2,2)+tt115&
&*X(2,1))*X(2,4)-(tt248*X(2,2)+tt164)*X(2,3)-tt247*tt4-tt246*X(2,1&
&)*X(2,2))*tt214*tt215
jac(1,8) = tt70*(tt11*((tt43*tt8+tt211*X(3,2)*X(3,3)+tt43*tt7+3*t&
&t43*tt6+2*tt206*X(2,3)+tt203+tt201+tt199+tt198+tt97+tt138)*X(3,4)&
&+(tt129*X(3,2)+tt111*X(3,1))*tt8+(tt134*tt7+tt178*X(3,1)*X(3,2)+(&
&2*tt35*X(2,3)+tt173*X(2,2))*X(2,4)+2*tt166*X(2,3)+tt162+tt160)*X(&
&3,3)+tt129*tt76+tt111*X(3,1)*tt7+((2*tt141*X(2,3)+tt143*X(2,2))*X&
&(2,4)+3*tt129*tt6+2*tt135*X(2,3)+tt130+tt128+tt92+tt125+tt124+tt1&
&23)*X(3,2)+((2*tt24*X(2,3)+tt115*X(2,2))*X(2,4)+3*tt111*tt6+2*tt1&
&09*X(2,2)*X(2,3)+tt105+tt101+tt100+tt99+tt98)*X(3,1))-(tt88*tt244&
&*tt214)/2.0E+0)*tt215-tt11*((-(tt68*X(3,3)+(tt81+X(2,2)+X(2,1))*X&
&(3,2)+tt88*X(3,1))*X(3,4))-((X(2,4)+tt58+X(2,1))*X(3,2)+tt253)*X(&
&3,3)-(tt74+tt87+tt67)*tt7-(X(2,4)+tt81+X(2,2))*X(3,1)*X(3,2)-tt45&
&*X(2,4)-2*tt37*X(2,3)-tt33*X(2,2)-tt26*X(2,1))*tt214*tt215
jac(1,9) = tt70*(tt11*((2*tt208*X(3,3)+tt212*X(3,2))*X(3,4)+3*tt1&
&91*tt8+2*tt190*X(3,3)+tt185+tt180+tt176+tt167+tt163+tt158+tt156+t&
&t154+tt151)-((2*X(3,3)-2*X(3,2))*tt244*tt214)/2.0E+0)*tt215-tt11*&
&((-tt69*X(3,4))-2*tt64*X(3,3)-tt61*X(3,2)-tt56*X(3,1))*tt214*tt21&
&5
jac(1,10) = tt11*(tt68*tt84+((tt51+tt58+3*X(2,1))*X(3,2)+tt49*X(3&
&,1))*tt8+((tt87+X(2,2)-3*X(2,1))*tt7+tt245*X(3,1)*X(3,2)+tt68*tt6&
&+(tt250+tt241)*X(2,3)+tt79+tt240+tt85+tt77*X(2,1))*X(3,3)+tt52*tt&
&76+tt49*X(3,1)*tt7+(tt86+(tt82+X(2,1))*tt6+(tt54+tt254+tt13+tt72+&
&tt23)*X(2,3)+tt235+tt39)*X(3,2)+(tt71-3*X(2,2)*tt6+(3*tt4+tt3+tt2&
&+tt1)*X(2,3)+tt80+tt78)*X(3,1))*tt70*tt215-tt11*((-tt35*tt8)-(tt3&
&1*X(3,2)+tt255)*X(3,3)-tt20*tt7-tt15*X(3,1)*X(3,2)-tt35*tt6-(tt31&
&*X(2,2)+tt24*X(2,1))*X(2,3)-tt20*tt4-tt15*X(2,1)*X(2,2))*tt214*tt&
&215
jac(1,11) = tt11*(tt35*tt84+(tt188*X(3,2)+tt255)*tt8+(tt183*tt7+t&
&t115*X(3,1)*X(3,2)+tt175+tt174+tt172+tt171+tt170+tt98+tt123)*X(3,&
&3)+tt141*tt76+tt24*X(3,1)*tt7+tt144*X(3,2)+tt116*X(3,1))*tt70*tt2&
&15-tt11*((-tt63*tt8)-(tt60*X(3,2)+tt55*X(3,1))*X(3,3)-tt52*tt7-tt&
&49*X(3,1)*X(3,2)-tt45*X(2,3)-tt41*X(2,2)-tt38*X(2,1))*tt214*tt215
jac(1,12) = tt11*tt213*tt70*tt215-((-tt69*X(3,3))-tt66*X(3,2)-tt6&
&5*X(3,1))*tt11*tt214*tt215
END 
SUBROUTINE calc_dih_angle_hes(hes, X) 
IMPLICIT NONE 
REAL(KIND=8) hes(12, 12) 
REAL(KIND=8) X(3, 4) 
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
tt1 = X(1,2)**2
tt2 = -2*X(1,2)*X(1,3)
tt3 = X(1,3)**2
tt4 = X(2,2)**2
tt5 = -2*X(2,2)*X(2,3)
tt6 = X(2,3)**2
tt7 = X(3,2)**2
tt8 = -2*X(3,2)*X(3,3)
tt9 = X(3,3)**2
tt10 = tt9+tt8+tt7+tt6+tt5+tt4+tt3+tt2+tt1
tt11 = sqrt(tt10)
tt12 = 1/tt11
tt13 = X(1,2)*X(1,3)
tt14 = -tt3
tt15 = -X(1,2)
tt16 = X(1,3)+tt15
tt17 = tt16*X(1,4)
tt18 = tt17+tt14+tt13
tt19 = -X(1,1)*X(1,3)
tt20 = -X(1,3)
tt21 = tt20+X(1,1)
tt22 = tt21*X(1,4)
tt23 = tt22+tt3+tt19
tt24 = -tt1
tt25 = tt20+X(1,2)
tt26 = tt25*X(1,4)
tt27 = tt26+tt13+tt24
tt28 = X(1,1)*X(1,2)
tt29 = -2*X(1,2)
tt30 = tt29+X(1,1)
tt31 = tt30*X(1,3)
tt32 = -2*X(1,1)
tt33 = X(1,3)+X(1,2)+tt32
tt34 = tt33*X(1,4)
tt35 = tt34+tt31+tt28
tt36 = -X(1,1)*X(1,2)
tt37 = tt15+X(1,1)
tt38 = tt37*X(1,4)
tt39 = tt38+tt1+tt36
tt40 = tt3+tt2+tt1
tt41 = tt40*X(2,1)
tt42 = X(1,2)+X(1,1)
tt43 = tt42*X(1,3)
tt44 = tt14+tt43+tt36
tt45 = -X(1,1)
tt46 = X(1,2)+tt45
tt47 = tt46*X(1,3)
tt48 = tt47+tt24+tt28
tt49 = X(2,2)*X(2,3)
tt50 = -tt6
tt51 = -X(2,2)
tt52 = X(2,3)+tt51
tt53 = tt52*X(2,4)+tt50+tt49+tt17+tt14+tt13
tt54 = -X(2,3)
tt55 = tt54+X(2,1)
tt56 = tt55*X(2,4)+tt6-X(2,1)*X(2,3)+tt22+tt3+tt19
tt57 = -tt4
tt58 = tt54+X(2,2)
tt59 = tt58*X(2,4)+tt49+tt57+tt26+tt13+tt24
tt60 = X(2,1)*X(2,2)
tt61 = -2*X(2,2)
tt62 = tt61+X(2,1)
tt63 = -2*X(2,1)
tt64 = X(2,3)+X(2,2)+tt63
tt65 = tt64*X(2,4)+tt62*X(2,3)+tt60+tt34+tt31+tt28
tt66 = -X(2,1)*X(2,2)
tt67 = tt51+X(2,1)
tt68 = tt67*X(2,4)+tt4+tt66+tt38+tt1+tt36
tt69 = tt6+tt5+tt4+tt3+tt2+tt1
tt70 = X(2,2)+X(2,1)
tt71 = tt50+tt70*X(2,3)+tt66+tt14+tt43+tt36
tt72 = -X(2,1)
tt73 = X(2,2)+tt72
tt74 = tt73*X(2,3)+tt57+tt60+tt47+tt24+tt28
tt75 = (-(tt74*X(3,3)+tt71*X(3,2)+tt69*X(3,1))*X(3,4))-tt68*tt9-(&
&tt65*X(3,2)+tt59*X(3,1))*X(3,3)-tt56*tt7-tt53*X(3,1)*X(3,2)-(tt48&
&*X(2,3)+tt44*X(2,2)+tt41)*X(2,4)-tt39*tt6-(tt35*X(2,2)+tt27*X(2,1&
&))*X(2,3)-tt23*tt4-tt18*X(2,1)*X(2,2)
tt76 = (tt4+tt3+tt2+tt1)*X(2,3)
tt77 = -2*X(2,2)*tt6
tt78 = X(2,3)**3
tt79 = 2*X(1,2)*X(1,3)
tt80 = 2*X(2,2)*X(2,3)
tt81 = tt50+tt80+tt57+tt14+tt79+tt24
tt82 = tt81*X(2,4)
tt83 = -X(2,4)
tt84 = tt83+X(2,3)
tt85 = X(3,2)**3
tt86 = tt14+tt79+tt24
tt87 = tt86*X(2,2)
tt88 = X(2,2)**3
tt89 = -tt88
tt90 = 2*tt4*X(2,3)
tt91 = -X(2,2)*tt6
tt92 = tt69*X(2,4)
tt93 = -2*X(2,3)
tt94 = 3*X(2,4)+tt93+tt51
tt95 = tt94*tt7
tt96 = 2*X(2,2)
tt97 = (-3*X(2,4))+X(2,3)+tt96
tt98 = X(2,4)+tt51
tt99 = X(3,3)**3
tt100 = tt40*X(2,2)
tt101 = -3*tt4
tt102 = -tt78
tt103 = 2*X(2,3)
tt104 = tt103+tt61
tt105 = tt58*tt9+tt104*X(3,2)*X(3,3)+tt58*tt7+tt102+3*X(2,2)*tt6+&
&(tt101+tt14+tt79+tt24)*X(2,3)+tt88+tt100
tt106 = tt105*X(3,4)+tt98*tt99+tt97*X(3,2)*tt9+(tt95+tt92+tt91+tt&
&90+tt89+tt87)*X(3,3)+tt84*tt85+(tt82+tt78+tt77+tt76)*X(3,2)
tt107 = X(1,4)+tt20
tt108 = -2*X(1,4)
tt109 = tt108+X(1,3)+X(1,2)
tt110 = X(1,4)+tt15
tt111 = (-(tt25*X(3,3)+tt16*X(3,2))*X(3,4))-tt110*tt9-tt109*X(3,2&
&)*X(3,3)-tt107*tt7-(tt25*X(2,3)+tt16*X(2,2))*X(2,4)-tt110*tt6-tt1&
&09*X(2,2)*X(2,3)-tt107*tt4
tt112 = 1/tt10
tt113 = tt1*X(1,3)
tt114 = -2*X(1,2)*tt3
tt115 = X(1,3)**3
tt116 = tt86*X(1,4)
tt117 = (tt116+tt115+tt114+tt113)*X(2,2)
tt118 = -X(1,4)
tt119 = tt118+X(1,3)
tt120 = tt119*tt88
tt121 = X(1,2)**3
tt122 = -tt121
tt123 = 2*tt1*X(1,3)
tt124 = -X(1,2)*tt3
tt125 = tt40*X(1,4)
tt126 = -2*X(1,3)
tt127 = 3*X(1,4)
tt128 = tt127+tt126+tt15
tt129 = tt128*tt4
tt130 = (tt129+tt125+tt124+tt123+tt122)*X(2,3)
tt131 = 2*X(1,2)
tt132 = -3*X(1,4)
tt133 = tt132+X(1,3)+tt131
tt134 = tt133*X(2,2)*tt6
tt135 = tt110*tt78
tt136 = -3*tt1*X(1,3)
tt137 = 3*X(1,2)*tt3
tt138 = -tt115
tt139 = tt25*tt4
tt140 = 2*X(1,3)
tt141 = tt140+tt29
tt142 = tt141*X(2,2)*X(2,3)
tt143 = tt25*tt6
tt144 = tt143+tt142+tt139+tt138+tt137+tt136+tt121
tt145 = tt144*X(2,4)
tt146 = -tt1*X(1,3)
tt147 = 2*X(1,2)*tt3
tt148 = (tt125+tt138+tt147+tt146)*X(2,1)
tt149 = tt107*X(2,1)*tt4
tt150 = X(1,1)*tt1
tt151 = -2*X(1,1)*X(1,2)*X(1,3)
tt152 = X(1,1)*tt3
tt153 = tt108+tt140
tt154 = tt153*X(2,1)*X(2,2)
tt155 = tt118+X(1,1)
tt156 = tt155*tt4
tt157 = (tt156+tt154+tt116+tt152+tt151+tt150)*X(2,3)
tt158 = tt107*X(2,1)
tt159 = 2*X(1,4)
tt160 = tt159+tt32
tt161 = tt160*X(2,2)
tt162 = tt161+tt158
tt163 = tt162*tt6
tt164 = tt155*tt78
tt165 = -X(1,1)*tt1
tt166 = 2*X(1,1)*X(1,2)
tt167 = (tt1+tt166)*X(1,3)
tt168 = tt29+tt45
tt169 = tt168*tt3
tt170 = X(1,3)+tt45
tt171 = tt170*tt4
tt172 = 2*X(1,1)
tt173 = tt126+tt172
tt174 = tt173*X(2,2)*X(2,3)
tt175 = tt170*tt6
tt176 = tt175+tt174+tt171+tt115+tt169+tt167+tt165
tt177 = tt176*X(2,4)
tt178 = tt25*X(2,4)+tt110*X(2,3)+tt119*X(2,2)
tt179 = tt155*X(2,3)
tt180 = tt170*X(2,4)+tt179+tt158
tt181 = -2*tt1*X(1,3)
tt182 = X(1,2)*tt3
tt183 = (tt116+tt182+tt181+tt121)*X(2,1)
tt184 = 2*X(1,1)*X(1,2)*X(1,3)
tt185 = -X(1,1)*tt3
tt186 = (tt125+tt185+tt184+tt165)*X(2,2)
tt187 = tt118+X(1,2)
tt188 = tt187*X(2,1)*tt4
tt189 = X(1,4)+tt45
tt190 = tt189*tt88
tt191 = tt159+tt29
tt192 = tt191*X(2,1)*X(2,2)
tt193 = tt108+tt172
tt194 = tt193*tt4
tt195 = (tt194+tt192)*X(2,3)
tt196 = tt187*X(2,1)
tt197 = tt189*X(2,2)
tt198 = tt197+tt196
tt199 = tt198*tt6
tt200 = -2*X(1,1)*X(1,2)
tt201 = 2*tt1
tt202 = (tt201+tt200)*X(1,3)
tt203 = tt37*tt3
tt204 = tt37*tt4
tt205 = tt131+tt32
tt206 = tt205*X(2,2)*X(2,3)
tt207 = tt37*tt6
tt208 = (tt207+tt206+tt204+tt203+tt202+tt122+tt150)*X(2,4)
tt209 = tt159+tt126
tt210 = tt108+tt131
tt211 = tt141*X(2,4)
tt212 = tt211+tt210*X(2,3)+tt209*X(2,2)
tt213 = tt212*X(3,1)*X(3,2)
tt214 = tt132+tt140+X(1,2)
tt215 = tt160*X(2,3)
tt216 = 3*X(1,1)
tt217 = tt126+tt15+tt216
tt218 = tt217*X(2,4)+tt215+tt197+tt214*X(2,1)
tt219 = tt218*tt7
tt220 = tt127+tt20+tt29
tt221 = tt193*X(2,2)
tt222 = -3*X(1,1)
tt223 = X(1,3)+tt131+tt222
tt224 = tt223*X(2,4)+tt179+tt221+tt220*X(2,1)
tt225 = tt224*X(3,2)+tt178*X(3,1)
tt226 = tt37*X(2,4)+tt197+tt196
tt227 = 3*tt1*X(1,3)
tt228 = -3*X(1,2)*tt3
tt229 = (tt24+tt200)*X(1,3)
tt230 = tt131+X(1,1)
tt231 = tt230*tt3
tt232 = -2*tt1
tt233 = (tt232+tt166)*X(1,3)
tt234 = tt46*tt3
tt235 = tt126+tt131
tt236 = tt235*X(2,1)*X(2,2)
tt237 = tt140+X(1,2)+tt222
tt238 = tt237*tt4
tt239 = tt16*X(2,1)
tt240 = tt20+tt29+tt216
tt241 = tt240*X(2,2)+tt239
tt242 = tt46*X(2,3)
tt243 = tt242+tt21*X(2,2)+tt239
tt244 = tt235*X(2,1)
tt245 = tt140+tt32
tt246 = tt29+tt172
tt247 = tt246*X(2,3)+tt245*X(2,2)+tt244
tt248 = tt243*tt9+tt247*X(3,2)*X(3,3)+tt243*tt7+tt46*tt78+tt241*t&
&t6+(tt238+tt236+tt234+tt233+tt121+tt165)*X(2,3)+tt21*tt88+tt16*X(&
&2,1)*tt4+(tt138+tt231+tt229+tt150)*X(2,2)+(tt115+tt228+tt227+tt12&
&2)*X(2,1)
tt249 = tt248*X(3,4)+tt226*tt99+tt225*tt9+(tt219+tt213+tt208+tt19&
&9+tt195+tt190+tt188+tt186+tt183)*X(3,3)+tt180*tt85+tt178*X(3,1)*t&
&t7+(tt177+tt164+tt163+tt157+tt149+tt148)*X(3,2)+(tt145+tt135+tt13&
&4+tt130+tt120+tt117)*X(3,1)
tt250 = 2*tt112*tt106*tt249+2*tt111*tt75
tt251 = tt249**2
tt252 = tt112*tt251+tt75**2
tt253 = 1/tt252**2
tt254 = X(2,4)+tt54
tt255 = -2*X(2,4)
tt256 = tt255+X(2,3)+X(2,2)
tt257 = (-(tt58*X(3,3)+tt52*X(3,2))*X(3,4))-tt98*tt9-tt256*X(3,2)&
&*X(3,3)-tt254*tt7-tt40*X(2,4)-tt27*X(2,3)-tt18*X(2,2)
tt258 = tt107*tt4
tt259 = tt153*X(2,2)*X(2,3)
tt260 = tt107*tt6
tt261 = tt187*tt4
tt262 = tt191*X(2,2)*X(2,3)
tt263 = tt187*tt6
tt264 = tt214*tt7
tt265 = tt16*tt9+tt235*X(3,2)*X(3,3)+tt16*tt7+tt16*tt6+tt235*X(2,&
&2)*X(2,3)+tt16*tt4+tt115+tt228+tt227+tt122
tt266 = tt265*X(3,4)+tt187*tt99+tt220*X(3,2)*tt9+(tt264+tt263+tt2&
&62+tt261+tt116+tt182+tt181+tt121)*X(3,3)+tt107*tt85+(tt260+tt259+&
&tt258+tt125+tt138+tt147+tt146)*X(3,2)
tt267 = 2*tt112*tt266*tt249+2*tt257*tt75
tt268 = 1/tt252
tt269 = (-tt69*X(3,4))-tt59*X(3,3)-tt53*X(3,2)
tt270 = tt178*tt9+tt212*X(3,2)*X(3,3)+tt178*tt7+tt145+tt135+tt134&
&+tt130+tt120+tt117
tt271 = 2*tt112*tt270*tt249+2*tt269*tt75
tt272 = X(1,4)+X(1,3)+tt29
tt273 = X(1,4)+tt126+X(1,1)
tt274 = tt118+tt131+tt45
tt275 = X(1,3)+tt29+X(1,1)
tt276 = (-(tt275*X(3,3)+tt170*X(3,2)+tt235*X(3,1))*X(3,4))-tt274*&
&tt9-(tt273*X(3,2)+tt272*X(3,1))*X(3,3)-tt119*X(3,1)*X(3,2)-(tt275&
&*X(2,3)+tt170*X(2,2)+tt244)*X(2,4)-tt274*tt6-(tt273*X(2,2)+tt272*&
&X(2,1))*X(2,3)-tt119*X(2,1)*X(2,2)
tt277 = -2*tt3
tt278 = tt141*X(1,4)
tt279 = (tt278+tt277+tt79)*X(2,2)
tt280 = -3*tt1
tt281 = 4*X(1,2)*X(1,3)
tt282 = tt235*X(1,4)
tt283 = (tt57+tt282+tt14+tt281+tt280)*X(2,3)
tt284 = 2*X(2,2)*tt6
tt285 = 3*tt1
tt286 = -6*X(1,2)*X(1,3)
tt287 = 3*tt3
tt288 = tt6+tt5+tt4+tt287+tt286+tt285
tt289 = tt288*X(2,4)
tt290 = 2*tt3
tt291 = tt282+tt290+tt2
tt292 = tt291*X(2,1)
tt293 = -2*X(1,1)*X(1,3)
tt294 = tt278+tt293+tt166
tt295 = tt294*X(2,3)
tt296 = tt277+(tt131+tt172)*X(1,3)+tt200
tt297 = tt296*X(2,4)
tt298 = -4*X(1,2)*X(1,3)
tt299 = (tt278+tt3+tt298+tt285)*X(2,1)
tt300 = 2*X(1,1)*X(1,3)
tt301 = (tt282+tt300+tt200)*X(2,2)
tt302 = X(2,1)*tt4
tt303 = -2*X(2,1)*X(2,2)*X(2,3)
tt304 = X(2,1)*tt6
tt305 = 4*X(1,2)
tt306 = (tt305+tt32)*X(1,3)
tt307 = (tt50+tt80+tt57+tt14+tt306+tt280+tt166)*X(2,4)
tt308 = tt255+tt103
tt309 = tt308*X(3,1)*X(3,2)
tt310 = tt83+X(2,1)
tt311 = tt310*tt7
tt312 = 2*X(2,4)
tt313 = tt312+tt63
tt314 = tt313*X(3,2)+tt254*X(3,1)
tt315 = 6*X(1,2)*X(1,3)
tt316 = -3*tt3
tt317 = (tt29+tt32)*X(1,3)
tt318 = -X(2,1)*tt4
tt319 = -4*X(1,2)
tt320 = (tt319+tt172)*X(1,3)
tt321 = 2*X(2,1)*X(2,2)
tt322 = tt61+tt72
tt323 = X(2,3)+tt72
tt324 = 2*X(2,1)
tt325 = tt93+tt324
tt326 = tt323*tt9+tt325*X(3,2)*X(3,3)+tt323*tt7+tt78+tt322*tt6+(t&
&t4+tt321+tt3+tt320+tt285+tt200)*X(2,3)+tt318+(tt290+tt317+tt166)*&
&X(2,2)+(tt316+tt315+tt280)*X(2,1)
tt327 = tt326*X(3,4)+tt310*tt99+tt314*tt9+(tt311+tt309+tt307+tt30&
&4+tt303+tt302+tt301+tt299)*X(3,3)+tt254*X(3,1)*tt7+(tt297+tt295+t&
&t292)*X(3,2)+(tt289+tt102+tt284+tt283+tt279)*X(3,1)
tt328 = 1/tt10**2
tt329 = (-tt235*tt328*tt251)+2*tt112*tt327*tt249+2*tt276*tt75
tt330 = tt235*X(2,3)
tt331 = tt141*X(2,2)
tt332 = tt235*X(2,4)
tt333 = tt235*X(2,2)
tt334 = tt141*X(2,3)
tt335 = tt334+tt333
tt336 = tt335*X(3,4)+(tt332+tt331)*X(3,3)+(tt211+tt330)*X(3,2)
tt337 = 1/tt11**3
tt338 = -X(2,2)*X(2,3)
tt339 = -tt52*X(2,4)
tt340 = -X(3,2)*X(3,3)
tt341 = -X(3,2)
tt342 = -(X(3,3)+tt341)*X(3,4)
tt343 = -tt12*(tt342+tt9+tt340+tt339+tt6+tt338)*tt249*tt268
tt344 = X(2,4)+X(2,3)+tt61
tt345 = X(2,4)+tt93+X(2,1)
tt346 = tt83+tt96+tt72
tt347 = tt93+tt96
tt348 = tt347*X(3,1)
tt349 = X(2,3)+tt61+X(2,1)
tt350 = (-(tt349*X(3,3)+tt323*X(3,2)+tt348)*X(3,4))-tt346*tt9-(tt&
&345*X(3,2)+tt344*X(3,1))*X(3,3)-tt84*X(3,1)*X(3,2)-tt44*X(2,4)-tt&
&35*X(2,3)-2*tt23*X(2,2)-tt18*X(2,1)
tt351 = 3*tt119*tt4
tt352 = 2*tt128*X(2,2)*X(2,3)
tt353 = tt133*tt6
tt354 = tt334+2*tt25*X(2,2)
tt355 = tt354*X(2,4)
tt356 = 2*tt107*X(2,1)*X(2,2)
tt357 = tt153*X(2,1)
tt358 = 2*tt155*X(2,2)
tt359 = (tt358+tt357)*X(2,3)
tt360 = tt160*tt6
tt361 = tt173*X(2,3)+2*tt170*X(2,2)
tt362 = tt361*X(2,4)
tt363 = 2*tt187*X(2,1)*X(2,2)
tt364 = 3*tt189*tt4
tt365 = tt191*X(2,1)
tt366 = 2*tt193*X(2,2)
tt367 = (tt366+tt365)*X(2,3)
tt368 = tt189*tt6
tt369 = tt205*X(2,3)+2*tt37*X(2,2)
tt370 = tt369*X(2,4)
tt371 = tt209*X(3,1)*X(3,2)
tt372 = tt189*tt7
tt373 = tt193*X(3,2)+tt119*X(3,1)
tt374 = 2*tt237*X(2,2)
tt375 = tt21*tt9+tt245*X(3,2)*X(3,3)+tt21*tt7+tt240*tt6+(tt374+tt&
&244)*X(2,3)+3*tt21*tt4+2*tt16*X(2,1)*X(2,2)+tt138+tt231+tt229+tt1&
&50
tt376 = tt375*X(3,4)+tt189*tt99+tt373*tt9+(tt372+tt371+tt370+tt36&
&8+tt367+tt364+tt363+tt125+tt185+tt184+tt165)*X(3,3)+tt119*X(3,1)*&
&tt7+(tt362+tt360+tt359+tt356)*X(3,2)+(tt355+tt353+tt352+tt351+tt1&
&16+tt115+tt114+tt113)*X(3,1)
tt377 = (-tt347*tt328*tt251)+2*tt112*tt376*tt249+2*tt350*tt75
tt378 = -2*tt6
tt379 = tt104*X(2,4)
tt380 = 4*X(2,2)*X(2,3)
tt381 = tt347*X(2,4)
tt382 = -tt7
tt383 = -tt99
tt384 = 3*tt4
tt385 = -6*X(2,2)*X(2,3)
tt386 = 3*tt6
tt387 = (tt9+tt8+tt7+tt386+tt385+tt384+tt3+tt2+tt1)*X(3,4)+tt383+&
&2*X(3,2)*tt9+(tt382+tt381+tt50+tt380+tt101+tt14+tt79+tt24)*X(3,3)&
&+(tt379+tt378+tt80)*X(3,2)
tt388 = -((-tt16*X(2,4))-tt109*X(2,3)-2*tt107*X(2,2))*tt12*tt249*&
&tt268
tt389 = (-tt71*X(3,4))-tt65*X(3,3)-2*tt56*X(3,2)-tt53*X(3,1)
tt390 = tt212*X(3,1)
tt391 = 2*tt218*X(3,2)
tt392 = tt247*X(3,3)+2*tt243*X(3,2)
tt393 = tt392*X(3,4)+tt224*tt9+(tt391+tt390)*X(3,3)+3*tt180*tt7+2&
&*tt178*X(3,1)*X(3,2)+tt177+tt164+tt163+tt157+tt149+tt148
tt394 = 2*X(3,2)
tt395 = -2*X(3,3)
tt396 = tt395+tt394
tt397 = (-tt396*tt328*tt251)+2*tt112*tt393*tt249+2*tt389*tt75
tt398 = (tt104*X(3,3)+2*tt58*X(3,2))*X(3,4)+tt97*tt9+2*tt94*X(3,2&
&)*X(3,3)+3*tt84*tt7+tt82+tt78+tt77+tt76
tt399 = -tt12*((-tt16*X(3,4))-tt109*X(3,3)-2*tt107*X(3,2))*tt249*&
&tt268
tt400 = X(1,4)+tt126+X(1,2)
tt401 = tt118+tt140+tt45
tt402 = X(1,4)+tt29+X(1,1)
tt403 = tt141*X(2,1)
tt404 = tt126+X(1,2)+X(1,1)
tt405 = tt141*X(3,1)
tt406 = (-(tt46*X(3,3)+tt404*X(3,2)+tt405)*X(3,4))-(tt402*X(3,2)+&
&tt187*X(3,1))*X(3,3)-tt401*tt7-tt400*X(3,1)*X(3,2)-(tt242+tt404*X&
&(2,2)+tt403)*X(2,4)-(tt402*X(2,2)+tt196)*X(2,3)-tt401*tt4-tt400*X&
&(2,1)*X(2,2)
tt407 = (tt282+tt287+tt298+tt1)*X(2,2)
tt408 = -2*tt4
tt409 = (tt408+tt278+tt2+tt201)*X(2,3)
tt410 = X(2,2)*tt6
tt411 = tt50+tt80+tt57+tt316+tt315+tt280
tt412 = tt411*X(2,4)
tt413 = (tt278+tt316+tt281+tt24)*X(2,1)
tt414 = (tt321+tt282+tt300+tt200)*X(2,3)
tt415 = -X(2,1)*tt6
tt416 = tt6+tt5+tt4+tt287+2*tt168*X(1,3)+tt1+tt166
tt417 = tt416*X(2,4)
tt418 = tt83+X(2,2)
tt419 = X(2,4)+tt72
tt420 = (tt282+tt79+tt232)*X(2,1)
tt421 = tt294*X(2,2)
tt422 = 2*tt37*X(1,3)
tt423 = (tt422+tt201+tt200)*X(2,4)
tt424 = tt312+tt61
tt425 = tt424*X(3,1)*X(3,2)
tt426 = tt255+tt324
tt427 = tt426*tt7
tt428 = tt418*X(3,1)
tt429 = tt419*X(3,2)+tt428
tt430 = 2*tt230*X(1,3)
tt431 = 2*tt46*X(1,3)
tt432 = -2*X(2,1)*X(2,2)
tt433 = 2*tt4
tt434 = tt96+tt63
tt435 = tt67*tt9+tt434*X(3,2)*X(3,3)+tt67*tt7+tt67*tt6+(tt433+tt4&
&32+tt431+tt232+tt166)*X(2,3)+tt89+tt302+(tt316+tt430+tt24+tt200)*&
&X(2,2)+(tt287+tt286+tt285)*X(2,1)
tt436 = tt435*X(3,4)+tt429*tt9+(tt427+tt425+tt423+tt421+tt420)*X(&
&3,3)+tt419*tt85+tt418*X(3,1)*tt7+(tt417+tt415+tt414+tt318+tt413)*&
&X(3,2)+(tt412+tt410+tt409+tt88+tt407)*X(3,1)
tt437 = (-tt141*tt328*tt251)+2*tt112*tt436*tt249+2*tt406*tt75
tt438 = tt330+tt331
tt439 = tt438*X(3,4)+(tt211+tt333)*X(3,3)+(tt332+tt334)*X(3,2)
tt440 = -tt58*X(2,4)
tt441 = -(X(3,2)-X(3,3))*X(3,4)
tt442 = -tt12*(tt441+tt340+tt7+tt440+tt338+tt4)*tt249*tt268
tt443 = X(2,4)+tt93+X(2,2)
tt444 = tt83+tt103+tt72
tt445 = X(2,4)+tt61+X(2,1)
tt446 = tt93+X(2,2)+X(2,1)
tt447 = (-(tt73*X(3,3)+tt446*X(3,2)+tt104*X(3,1))*X(3,4))-(tt445*&
&X(3,2)+tt428)*X(3,3)-tt444*tt7-tt443*X(3,1)*X(3,2)-tt48*X(2,4)-2*&
&tt39*X(2,3)-tt35*X(2,2)-tt27*X(2,1)
tt448 = 2*tt133*X(2,2)*X(2,3)
tt449 = 3*tt110*tt6
tt450 = 2*tt25*X(2,3)+tt331
tt451 = tt450*X(2,4)
tt452 = 2*tt162*X(2,3)
tt453 = 3*tt155*tt6
tt454 = tt173*X(2,2)
tt455 = 2*tt170*X(2,3)+tt454
tt456 = tt455*X(2,4)
tt457 = 2*tt198*X(2,3)
tt458 = 2*tt37*X(2,3)+tt205*X(2,2)
tt459 = tt458*X(2,4)
tt460 = tt210*X(3,1)*X(3,2)
tt461 = tt160*tt7
tt462 = tt155*X(3,2)+tt110*X(3,1)
tt463 = tt46*tt9+tt246*X(3,2)*X(3,3)+tt46*tt7+3*tt46*tt6+2*tt241*&
&X(2,3)+tt238+tt236+tt234+tt233+tt121+tt165
tt464 = tt463*X(3,4)+tt462*tt9+(tt461+tt460+tt459+tt457+tt194+tt1&
&92)*X(3,3)+tt155*tt85+tt110*X(3,1)*tt7+(tt456+tt453+tt452+tt156+t&
&t154+tt116+tt152+tt151+tt150)*X(3,2)+(tt451+tt449+tt448+tt129+tt1&
&25+tt124+tt123+tt122)*X(3,1)
tt465 = (-tt104*tt328*tt251)+2*tt112*tt464*tt249+2*tt447*tt75
tt466 = -4*X(2,2)*X(2,3)
tt467 = -2*tt7
tt468 = 6*X(2,2)*X(2,3)
tt469 = -3*tt6
tt470 = 2*X(3,2)*X(3,3)
tt471 = -tt9
tt472 = (tt471+tt470+tt382+tt469+tt468+tt101+tt14+tt79+tt24)*X(3,&
&4)+X(3,2)*tt9+(tt467+tt379+tt5+tt433)*X(3,3)+tt85+(tt381+tt386+tt&
&466+tt4+tt3+tt2+tt1)*X(3,2)
tt473 = -((-tt25*X(2,4))-2*tt110*X(2,3)-tt109*X(2,2))*tt12*tt249*&
&tt268
tt474 = (-tt74*X(3,4))-2*tt68*X(3,3)-tt65*X(3,2)-tt59*X(3,1)
tt475 = 2*tt243*X(3,3)+tt247*X(3,2)
tt476 = tt475*X(3,4)+3*tt226*tt9+2*tt225*X(3,3)+tt219+tt213+tt208&
&+tt199+tt195+tt190+tt188+tt186+tt183
tt477 = -2*X(3,2)
tt478 = 2*X(3,3)+tt477
tt479 = (-tt478*tt328*tt251)+2*tt112*tt476*tt249+2*tt474*tt75
tt480 = (2*tt58*X(3,3)+tt104*X(3,2))*X(3,4)+3*tt98*tt9+2*tt97*X(3&
&,2)*X(3,3)+tt95+tt92+tt91+tt90+tt89+tt87
tt481 = -tt12*((-tt25*X(3,4))-2*tt110*X(3,3)-tt109*X(3,2))*tt249*&
&tt268
tt482 = tt25*X(3,1)
tt483 = (-tt37*tt9)-(tt33*X(3,2)+tt482)*X(3,3)-tt21*tt7-tt16*X(3,&
&1)*X(3,2)-tt37*tt6-(tt33*X(2,2)+tt25*X(2,1))*X(2,3)-tt21*tt4-tt16&
&*X(2,1)*X(2,2)
tt484 = (tt384+tt3+tt2+tt1)*X(2,3)
tt485 = -3*X(2,2)*tt6
tt486 = (tt57+tt432+tt14+tt79+tt24)*X(2,3)
tt487 = tt96+X(2,1)
tt488 = tt487*tt6
tt489 = tt86*X(2,1)
tt490 = (tt408+tt321)*X(2,3)
tt491 = tt73*tt6
tt492 = tt347*X(3,1)*X(3,2)
tt493 = tt103+X(2,2)-3*X(2,1)
tt494 = tt493*tt7
tt495 = tt54+tt61+3*X(2,1)
tt496 = tt495*X(3,2)+tt52*X(3,1)
tt497 = tt73*tt99+tt496*tt9+(tt494+tt492+tt491+tt490+tt88+tt318+t&
&t100+tt489)*X(3,3)+tt55*tt85+tt52*X(3,1)*tt7+(tt102+tt488+tt486+t&
&t302+tt41)*X(3,2)+(tt78+tt485+tt484+tt89+tt87)*X(3,1)
tt498 = 2*tt112*tt497*tt249+2*tt483*tt75
tt499 = -(tt471+tt470+tt382+tt50+tt80+tt57)*tt12*tt249*tt268
tt500 = (-tt67*tt9)-(tt64*X(3,2)+tt58*X(3,1))*X(3,3)-tt55*tt7-tt5&
&2*X(3,1)*X(3,2)-tt48*X(2,3)-tt44*X(2,2)-tt40*X(2,1)
tt501 = tt141*X(3,1)*X(3,2)
tt502 = tt217*tt7
tt503 = tt223*X(3,2)+tt482
tt504 = tt37*tt99+tt503*tt9+(tt502+tt501+tt207+tt206+tt204+tt203+&
&tt202+tt122+tt150)*X(3,3)+tt170*tt85+tt25*X(3,1)*tt7+tt176*X(3,2)&
&+tt144*X(3,1)
tt505 = 2*tt112*tt504*tt249+2*tt500*tt75
tt506 = -tt85
tt507 = tt12*(tt99-3*X(3,2)*tt9+(3*tt7+tt6+tt5+tt4+tt3+tt2+tt1)*X&
&(3,3)+tt506+tt81*X(3,2))*tt75*tt268
tt508 = -((-tt25*X(2,3))-tt16*X(2,2))*tt12*tt249*tt268
tt509 = (-tt74*X(3,3))-tt71*X(3,2)-tt69*X(3,1)
tt510 = 2*tt112*tt248*tt249+2*tt509*tt75
tt511 = tt12*tt105*tt75*tt268
tt512 = -((-tt25*X(3,3))-tt16*X(3,2))*tt12*tt249*tt268
tt513 = (tt471+tt470+tt382+tt50+tt80+tt57+tt316+tt315+tt280)*X(3,&
&4)+tt99-2*X(3,2)*tt9+(tt7+tt6+tt5+tt4+tt278+tt3+tt298+tt285)*X(3,&
&3)+tt291*X(3,2)
tt514 = -((-tt235*X(2,4))-tt272*X(2,3)-tt119*X(2,2))*tt12*tt249*t&
&t268
tt515 = tt191*X(2,3)
tt516 = (tt330+2*tt16*X(2,2))*X(3,4)+(tt515+2*tt187*X(2,2))*X(3,3&
&)+(tt153*X(2,3)+2*tt107*X(2,2))*X(3,2)
tt517 = -X(1,2)*X(1,3)
tt518 = -tt16*X(1,4)
tt519 = -tt12*(tt342+tt9+tt340+tt518+tt3+tt517)*tt249*tt268
tt520 = (tt235*X(3,3)+2*tt16*X(3,2))*X(3,4)+tt220*tt9+2*tt214*X(3&
&,2)*X(3,3)+3*tt107*tt7+tt260+tt259+tt258+tt125+tt138+tt147+tt146
tt521 = -tt12*((-tt52*X(3,4))-tt256*X(3,3)-2*tt254*X(3,2))*tt249*&
&tt268
tt522 = 2*tt7
tt523 = (tt9+tt8+tt7+tt6+tt5+tt4+tt287+tt286+tt285)*X(3,4)-X(3,2)&
&*tt9+(tt522+tt282+tt79+tt232)*X(3,3)+tt506+(tt50+tt80+tt57+tt278+&
&tt316+tt281+tt24)*X(3,2)
tt524 = -((-tt141*X(2,4))-tt187*X(2,3)-tt400*X(2,2))*tt12*tt249*t&
&t268
tt525 = tt153*X(2,2)
tt526 = (2*tt16*X(2,3)+tt333)*X(3,4)+(2*tt187*X(2,3)+tt191*X(2,2)&
&)*X(3,3)+(2*tt107*X(2,3)+tt525)*X(3,2)
tt527 = -tt25*X(1,4)
tt528 = -tt12*(tt441+tt340+tt7+tt527+tt517+tt1)*tt249*tt268
tt529 = (2*tt16*X(3,3)+tt235*X(3,2))*X(3,4)+3*tt187*tt9+2*tt220*X&
&(3,2)*X(3,3)+tt264+tt263+tt262+tt261+tt116+tt182+tt181+tt121
tt530 = -tt12*((-tt58*X(3,4))-2*tt98*X(3,3)-tt256*X(3,2))*tt249*t&
&t268
tt531 = tt12*(tt383+3*X(3,2)*tt9+((-3*tt7)+tt50+tt80+tt57+tt14+tt&
&79+tt24)*X(3,3)+tt85+tt69*X(3,2))*tt75*tt268
tt532 = -(tt471+tt470+tt382+tt14+tt79+tt24)*tt12*tt249*tt268
tt533 = tt12*tt265*tt75*tt268
tt534 = -((-tt58*X(3,3))-tt52*X(3,2))*tt12*tt249*tt268
tt535 = tt254*tt9+tt308*X(3,2)*X(3,3)+tt254*tt7+tt289+tt102+tt284&
&+tt283+tt279
tt536 = -tt12*((-tt235*X(3,4))-tt272*X(3,3)-tt119*X(3,2))*tt249*t&
&t268
tt537 = tt119*tt9+tt209*X(3,2)*X(3,3)+tt119*tt7+tt355+tt353+tt352&
&+tt351+tt116+tt115+tt114+tt113
tt538 = -tt12*((-tt347*X(3,4))-tt344*X(3,3)-tt84*X(3,2))*tt249*tt&
&268
tt539 = tt212*X(3,3)+2*tt178*X(3,2)
tt540 = -(tt339+tt6+tt338+tt518+tt3+tt517)*tt12*tt249*tt268
tt541 = tt418*tt9+tt424*X(3,2)*X(3,3)+tt418*tt7+tt412+tt410+tt409&
&+tt88+tt407
tt542 = -tt12*((-tt141*X(3,4))-tt187*X(3,3)-tt400*X(3,2))*tt249*t&
&t268
tt543 = tt110*tt9+tt210*X(3,2)*X(3,3)+tt110*tt7+tt451+tt449+tt448&
&+tt129+tt125+tt124+tt123+tt122
tt544 = -tt12*((-tt104*X(3,4))-tt418*X(3,3)-tt443*X(3,2))*tt249*t&
&t268
tt545 = 2*tt178*X(3,3)+tt212*X(3,2)
tt546 = -(tt440+tt338+tt4+tt527+tt517+tt1)*tt12*tt249*tt268
tt547 = tt12*(tt52*tt9+tt347*X(3,2)*X(3,3)+tt52*tt7+tt78+tt485+tt&
&484+tt89+tt87)*tt75*tt268
tt548 = tt12*(tt25*tt9+tt141*X(3,2)*X(3,3)+tt25*tt7+tt143+tt142+t&
&t139+tt138+tt137+tt136+tt121)*tt75*tt268
tt549 = -tt81*tt12*tt249*tt268
tt550 = tt12*tt327-(tt235*tt337*tt249)/2.0E+0
tt551 = 2*X(2,1)*X(2,3)
tt552 = 2*X(3,1)*X(3,3)
tt553 = -2*tt9
tt554 = 2*X(3,1)
tt555 = -(tt395+tt554)*X(3,4)
tt556 = -6*X(1,2)
tt557 = 4*X(1,3)
tt558 = 6*X(1,2)
tt559 = -6*X(1,3)
tt560 = tt559+tt558
tt561 = tt560*X(2,4)
tt562 = tt193*X(2,3)
tt563 = -4*X(1,3)
tt564 = 6*X(1,3)
tt565 = tt564+tt556
tt566 = tt565*X(2,1)
tt567 = 1/tt11**5
tt568 = -tt337*tt249
tt569 = -((-tt170*X(2,4))-tt273*X(2,3)-tt119*X(2,1))*tt12*tt249*t&
&t268
tt570 = 2*tt6
tt571 = tt379-2*X(2,1)*X(2,3)+tt321+tt282+tt300+tt200
tt572 = tt75*((3.0E+0*tt235*tt347*tt567*tt249)/4.0E+0-(tt347*tt33&
&7*tt327)/2.0E+0-(tt235*tt337*tt376)/2.0E+0+tt12*((tt378+(tt96+tt3&
&24)*X(2,3)+tt432+tt290+tt317+tt166)*X(3,4)+tt571*X(3,3)+(tt381+tt&
&570+tt5+tt278+tt277+tt79)*X(3,1)))*tt268
tt573 = -tt12*((-tt170*X(3,4))-tt273*X(3,3)-tt119*X(3,1))*tt249*t&
&t268
tt574 = tt75*((3.0E+0*tt235*tt396*tt567*tt249)/4.0E+0-(tt396*tt33&
&7*tt327)/2.0E+0-(tt235*tt337*tt393)/2.0E+0+tt12*((tt325*X(3,3)+2*&
&tt323*X(3,2))*X(3,4)+tt313*tt9+(2*tt310*X(3,2)+tt308*X(3,1))*X(3,&
&3)+2*tt254*X(3,1)*X(3,2)+tt297+tt295+tt292))*tt268
tt575 = -tt62*X(2,3)
tt576 = -tt64*X(2,4)
tt577 = -X(3,1)*X(3,2)
tt578 = -(tt477+X(3,1))*X(3,3)
tt579 = -(X(3,3)+X(3,2)-2*X(3,1))*X(3,4)
tt580 = -tt12*(tt579+tt578+tt577+tt576+tt575+tt66)*tt249*tt268
tt581 = tt337*tt249
tt582 = tt75*(tt581+(3.0E+0*tt235*tt141*tt567*tt249)/4.0E+0-(tt14&
&1*tt337*tt327)/2.0E+0-(tt235*tt337*tt436)/2.0E+0+tt12*(((tt140+tt&
&319+tt172)*X(2,3)+(tt557+tt29+tt32)*X(2,2)+tt560*X(2,1))*X(3,4)+(&
&(tt126+tt305+tt32)*X(2,4)+tt221+(tt159+tt140+tt319)*X(2,1))*X(3,3&
&)+((tt563+tt131+tt172)*X(2,4)+tt215+(tt108+tt557+tt29)*X(2,1))*X(&
&3,2)+(tt565*X(2,4)+(tt108+tt126+tt305)*X(2,3)+(tt159+tt563+tt131)&
&*X(2,2))*X(3,1)))*tt268
tt583 = -((-tt275*X(2,4))-2*tt274*X(2,3)-tt273*X(2,2)-tt272*X(2,1&
&))*tt12*tt249*tt268
tt584 = -X(3,1)*tt7
tt585 = 2*X(3,1)*X(3,2)
tt586 = tt75*((3.0E+0*tt235*tt104*tt567*tt249)/4.0E+0-(tt104*tt33&
&7*tt327)/2.0E+0-(tt235*tt337*tt464)/2.0E+0+tt12*((tt9+tt8+tt7+tt3&
&86+2*tt322*X(2,3)+tt4+tt321+tt3+tt320+tt285+tt200)*X(3,4)-X(3,1)*&
&tt9+(tt585+tt381+tt551+tt432)*X(3,3)+tt584+tt294*X(3,2)+(tt379+tt&
&469+tt380+tt57+tt282+tt14+tt281+tt280)*X(3,1)))*tt268
tt587 = -tt12*((-tt275*X(3,4))-2*tt274*X(3,3)-tt273*X(3,2)-tt272*&
&X(3,1))*tt249*tt268
tt588 = tt75*((3.0E+0*tt235*tt478*tt567*tt249)/4.0E+0-(tt478*tt33&
&7*tt327)/2.0E+0-(tt235*tt337*tt476)/2.0E+0+tt12*((2*tt323*X(3,3)+&
&tt325*X(3,2))*X(3,4)+3*tt310*tt9+2*tt314*X(3,3)+tt311+tt309+tt307&
&+tt304+tt303+tt302+tt301+tt299))*tt268
tt589 = (tt333+tt403)*X(3,3)+(tt334+tt244)*X(3,2)+tt438*X(3,1)
tt590 = -tt70*X(2,3)
tt591 = X(3,1)*X(3,2)
tt592 = -(X(3,2)+X(3,1))*X(3,3)
tt593 = -tt12*(tt9+tt592+tt591+tt6+tt590+tt60)*tt249*tt268
tt594 = X(3,1)*tt7
tt595 = -2*X(3,1)*X(3,2)
tt596 = tt383+(tt394+X(3,1))*tt9+(tt382+tt595+tt50+tt80+tt57+tt14&
&+tt306+tt280+tt166)*X(3,3)+tt594+tt296*X(3,2)+tt288*X(3,1)
tt597 = -((-tt275*X(2,3))-tt170*X(2,2)-tt235*X(2,1))*tt12*tt249*t&
&t268
tt598 = -((-tt275*X(3,3))-tt170*X(3,2)-tt235*X(3,1))*tt12*tt249*t&
&t268
tt599 = tt12*tt376-(tt347*tt337*tt249)/2.0E+0
tt600 = 2*tt25*X(2,4)
tt601 = 2*tt170*X(2,4)
tt602 = 2*tt37*X(2,4)
tt603 = -tt12*((-tt323*X(3,4))-tt345*X(3,3)-tt84*X(3,1))*tt249*tt&
&268
tt604 = tt75*((3.0E+0*tt347*tt396*tt567*tt249)/4.0E+0-(tt396*tt33&
&7*tt376)/2.0E+0-(tt347*tt337*tt393)/2.0E+0+tt12*((tt245*X(3,3)+2*&
&tt21*X(3,2))*X(3,4)+tt193*tt9+(2*tt189*X(3,2)+tt209*X(3,1))*X(3,3&
&)+2*tt119*X(3,1)*X(3,2)+tt362+tt360+tt359+tt356))*tt268
tt605 = -((-tt404*X(2,4))-tt402*X(2,3)-2*tt401*X(2,2)-tt400*X(2,1&
&))*tt12*tt249*tt268
tt606 = tt75*((3.0E+0*tt141*tt347*tt567*tt249)/4.0E+0-(tt347*tt33&
&7*tt436)/2.0E+0-(tt141*tt337*tt376)/2.0E+0+tt12*((tt471+tt470+tt3&
&82+tt50+(4*X(2,2)+tt63)*X(2,3)+tt101+tt321+tt316+tt430+tt24+tt200&
&)*X(3,4)+X(3,1)*tt9+(tt595+tt278+tt293+tt166)*X(3,3)+tt594+(tt381&
&+tt551+tt432)*X(3,2)+(tt379+tt6+tt466+tt384+tt282+tt287+tt298+tt1&
&)*X(3,1)))*tt268
tt607 = -tt30*X(1,3)
tt608 = -tt33*X(1,4)
tt609 = -tt12*(tt579+tt578+tt577+tt608+tt607+tt36)*tt249*tt268
tt610 = tt75*(tt581+(3.0E+0*tt347*tt104*tt567*tt249)/4.0E+0-(tt10&
&4*tt337*tt376)/2.0E+0-(tt347*tt337*tt464)/2.0E+0+tt12*((2*tt240*X&
&(2,3)+tt374+tt244)*X(3,4)+(tt205*X(2,4)+2*tt189*X(2,3)+tt366+tt36&
&5)*X(3,3)+(tt173*X(2,4)+2*tt160*X(2,3)+tt358+tt357)*X(3,2)+(tt211&
&+2*tt133*X(2,3)+2*tt128*X(2,2))*X(3,1)))*tt268
tt611 = -tt12*((-tt349*X(3,4))-2*tt346*X(3,3)-tt345*X(3,2)-tt344*&
&X(3,1))*tt249*tt268
tt612 = tt75*((3.0E+0*tt347*tt478*tt567*tt249)/4.0E+0-(tt478*tt33&
&7*tt376)/2.0E+0-(tt347*tt337*tt476)/2.0E+0+tt12*((2*tt21*X(3,3)+t&
&t245*X(3,2))*X(3,4)+3*tt189*tt9+2*tt373*X(3,3)+tt372+tt371+tt370+&
&tt368+tt367+tt364+tt363+tt125+tt185+tt184+tt165))*tt268
tt613 = -X(3,1)
tt614 = tt99+(tt477+tt613)*tt9+(tt7+tt585+tt6+(tt324-4*X(2,2))*X(&
&2,3)+tt384+tt432+tt3+tt2+tt1)*X(3,3)+tt584+(tt570+(tt61+tt63)*X(2&
&,3)+tt321)*X(3,2)+(tt469+tt468+tt101+tt14+tt79+tt24)*X(3,1)
tt615 = -((-tt33*X(2,3))-2*tt21*X(2,2)-tt16*X(2,1))*tt12*tt249*tt&
&268
tt616 = tt369*X(3,3)+tt361*X(3,2)+tt354*X(3,1)
tt617 = -tt42*X(1,3)
tt618 = -tt12*(tt9+tt592+tt591+tt3+tt617+tt28)*tt249*tt268
tt619 = -((-tt349*X(3,3))-tt323*X(3,2)-tt347*X(3,1))*tt12*tt249*t&
&t268
tt620 = tt12*tt393-(tt396*tt337*tt249)/2.0E+0
tt621 = 2*tt243*X(3,4)
tt622 = -tt12*((-tt404*X(3,4))-tt402*X(3,3)-2*tt401*X(3,2)-tt400*&
&X(3,1))*tt249*tt268
tt623 = tt75*((3.0E+0*tt141*tt396*tt567*tt249)/4.0E+0-(tt396*tt33&
&7*tt436)/2.0E+0-(tt141*tt337*tt393)/2.0E+0+tt12*((tt434*X(3,3)+2*&
&tt67*X(3,2))*X(3,4)+tt419*tt9+(2*tt426*X(3,2)+tt424*X(3,1))*X(3,3&
&)+3*tt419*tt7+2*tt418*X(3,1)*X(3,2)+tt417+tt415+tt414+tt318+tt413&
&))*tt268
tt624 = -tt12*((-tt446*X(3,4))-tt445*X(3,3)-2*tt444*X(3,2)-tt443*&
&X(3,1))*tt249*tt268
tt625 = tt75*((3.0E+0*tt104*tt396*tt567*tt249)/4.0E+0-(tt396*tt33&
&7*tt464)/2.0E+0-(tt104*tt337*tt393)/2.0E+0+tt12*((tt246*X(3,3)+2*&
&tt46*X(3,2))*X(3,4)+tt155*tt9+(2*tt160*X(3,2)+tt210*X(3,1))*X(3,3&
&)+3*tt155*tt7+2*tt110*X(3,1)*X(3,2)+tt456+tt453+tt452+tt156+tt154&
&+tt116+tt152+tt151+tt150))*tt268
tt626 = -(tt576+tt575+tt66+tt608+tt607+tt36)*tt12*tt249*tt268
tt627 = tt75*(tt581+(3.0E+0*tt396*tt478*tt567*tt249)/4.0E+0-(tt39&
&6*tt337*tt476)/2.0E+0-(tt478*tt337*tt393)/2.0E+0+tt12*(tt247*X(3,&
&4)+2*tt224*X(3,3)+tt391+tt390))*tt268
tt628 = tt495*tt9+(2*tt493*X(3,2)+tt348)*X(3,3)+3*tt55*tt7+2*tt52&
&*X(3,1)*X(3,2)+tt102+tt488+tt486+tt302+tt41
tt629 = -((-tt33*X(3,3))-2*tt21*X(3,2)-tt16*X(3,1))*tt12*tt249*tt&
&268
tt630 = tt223*tt9+(2*tt217*X(3,2)+tt405)*X(3,3)+3*tt170*tt7+2*tt2&
&5*X(3,1)*X(3,2)+tt175+tt174+tt171+tt115+tt169+tt167+tt165
tt631 = -((-tt64*X(3,3))-2*tt55*X(3,2)-tt52*X(3,1))*tt12*tt249*tt&
&268
tt632 = -(tt6+tt590+tt60+tt3+tt617+tt28)*tt12*tt249*tt268
tt633 = tt12*tt436-(tt141*tt337*tt249)/2.0E+0
tt634 = -(tt477+tt554)*X(3,4)
tt635 = -((-tt46*X(2,4))-tt402*X(2,2)-tt187*X(2,1))*tt12*tt249*tt&
&268
tt636 = tt75*((3.0E+0*tt141*tt104*tt567*tt249)/4.0E+0-(tt104*tt33&
&7*tt436)/2.0E+0-(tt141*tt337*tt464)/2.0E+0+tt12*((2*tt67*X(2,3)+t&
&t433+tt432+tt431+tt232+tt166)*X(3,4)+tt571*X(3,2)+(tt381+tt80+tt4&
&08+tt278+tt2+tt201)*X(3,1)))*tt268
tt637 = -tt12*((-tt46*X(3,4))-tt402*X(3,2)-tt187*X(3,1))*tt249*tt&
&268
tt638 = tt75*((3.0E+0*tt141*tt478*tt567*tt249)/4.0E+0-(tt478*tt33&
&7*tt436)/2.0E+0-(tt141*tt337*tt476)/2.0E+0+tt12*((2*tt67*X(3,3)+t&
&t434*X(3,2))*X(3,4)+2*tt429*X(3,3)+tt427+tt425+tt423+tt421+tt420)&
&)*tt268
tt639 = (tt331+tt244)*X(3,3)+(tt330+tt403)*X(3,2)+tt335*X(3,1)
tt640 = -tt73*X(2,3)
tt641 = X(3,2)+tt613
tt642 = -tt641*X(3,3)
tt643 = -(tt642+tt7+tt577+tt640+tt4+tt66)*tt12*tt249*tt268
tt644 = tt641*tt9+(tt467+tt585+tt422+tt201+tt200)*X(3,3)+tt85+tt5&
&84+tt416*X(3,2)+tt411*X(3,1)
tt645 = -((-tt46*X(2,3))-tt404*X(2,2)-tt141*X(2,1))*tt12*tt249*tt&
&268
tt646 = -((-tt46*X(3,3))-tt404*X(3,2)-tt141*X(3,1))*tt12*tt249*tt&
&268
tt647 = tt12*tt464-(tt104*tt337*tt249)/2.0E+0
tt648 = -tt12*((-tt73*X(3,4))-tt445*X(3,2)-tt418*X(3,1))*tt249*tt&
&268
tt649 = tt75*((3.0E+0*tt104*tt478*tt567*tt249)/4.0E+0-(tt478*tt33&
&7*tt464)/2.0E+0-(tt104*tt337*tt476)/2.0E+0+tt12*((2*tt46*X(3,3)+t&
&t246*X(3,2))*X(3,4)+2*tt462*X(3,3)+tt461+tt460+tt459+tt457+tt194+&
&tt192))*tt268
tt650 = (tt341+X(3,1))*tt9+(tt522+tt595+2*tt73*X(2,3)+tt408+tt321&
&)*X(3,3)+tt506+tt594+(tt469+2*tt487*X(2,3)+tt57+tt432+tt14+tt79+t&
&t24)*X(3,2)+(tt386+tt385+tt384+tt3+tt2+tt1)*X(3,1)
tt651 = -((-2*tt37*X(2,3))-tt33*X(2,2)-tt25*X(2,1))*tt12*tt249*tt&
&268
tt652 = tt458*X(3,3)+tt455*X(3,2)+tt450*X(3,1)
tt653 = -tt46*X(1,3)
tt654 = -(tt642+tt7+tt577+tt653+tt1+tt36)*tt12*tt249*tt268
tt655 = -((-tt73*X(3,3))-tt446*X(3,2)-tt104*X(3,1))*tt12*tt249*tt&
&268
tt656 = tt12*tt476-(tt478*tt337*tt249)/2.0E+0
tt657 = 3*tt73*tt9+2*tt496*X(3,3)+tt494+tt492+tt491+tt490+tt88+tt&
&318+tt100+tt489
tt658 = -((-2*tt37*X(3,3))-tt33*X(3,2)-tt25*X(3,1))*tt12*tt249*tt&
&268
tt659 = 3*tt37*tt9+2*tt503*X(3,3)+tt502+tt501+tt207+tt206+tt204+t&
&t203+tt202+tt122+tt150
tt660 = -((-2*tt67*X(3,3))-tt64*X(3,2)-tt58*X(3,1))*tt12*tt249*tt&
&268
tt661 = -(tt640+tt4+tt66+tt653+tt1+tt36)*tt12*tt249*tt268
hes(1,1) = tt12*tt111*tt249*tt250*tt253-tt12*tt75*tt106*tt250*tt2&
&53
hes(1,2) = tt12*tt257*tt106*tt268-tt12*tt111*tt266*tt268+tt12*tt1&
&11*tt249*tt267*tt253-tt12*tt75*tt106*tt267*tt253
hes(1,3) = tt12*tt269*tt106*tt268-tt12*tt270*tt111*tt268+tt12*tt1&
&11*tt249*tt271*tt253-tt12*tt75*tt106*tt271*tt253
hes(1,4) = (tt235*tt337*tt111*tt249*tt268)/2.0E+0+tt343-tt12*tt11&
&1*tt327*tt268-(tt235*tt337*tt75*tt106*tt268)/2.0E+0+tt12*tt276*tt&
&106*tt268+tt12*tt336*tt75*tt268+tt12*tt111*tt249*tt329*tt253-tt12&
&*tt75*tt106*tt329*tt253
hes(1,5) = (tt347*tt337*tt111*tt249*tt268)/2.0E+0+tt388-(tt347*tt&
&337*tt75*tt106*tt268)/2.0E+0+tt12*tt350*tt106*tt268-tt12*tt111*tt&
&376*tt268+tt12*tt75*tt387*tt268+tt12*tt111*tt249*tt377*tt253-tt12&
&*tt75*tt106*tt377*tt253
hes(1,6) = (tt396*tt337*tt111*tt249*tt268)/2.0E+0+tt399-(tt396*tt&
&337*tt75*tt106*tt268)/2.0E+0+tt12*tt389*tt106*tt268+tt12*tt398*tt&
&75*tt268-tt12*tt111*tt393*tt268+tt12*tt111*tt249*tt397*tt253-tt12&
&*tt75*tt106*tt397*tt253
hes(1,7) = (tt141*tt337*tt111*tt249*tt268)/2.0E+0+tt442-(tt141*tt&
&337*tt75*tt106*tt268)/2.0E+0+tt12*tt406*tt106*tt268-tt12*tt111*tt&
&436*tt268+tt12*tt439*tt75*tt268+tt12*tt111*tt249*tt437*tt253-tt12&
&*tt75*tt106*tt437*tt253
hes(1,8) = (tt104*tt337*tt111*tt249*tt268)/2.0E+0+tt473-(tt104*tt&
&337*tt75*tt106*tt268)/2.0E+0+tt12*tt447*tt106*tt268-tt12*tt111*tt&
&464*tt268+tt12*tt75*tt472*tt268+tt12*tt111*tt249*tt465*tt253-tt12&
&*tt75*tt106*tt465*tt253
hes(1,9) = (tt478*tt337*tt111*tt249*tt268)/2.0E+0+tt481-(tt478*tt&
&337*tt75*tt106*tt268)/2.0E+0+tt12*tt474*tt106*tt268+tt12*tt480*tt&
&75*tt268-tt12*tt111*tt476*tt268+tt12*tt111*tt249*tt479*tt253-tt12&
&*tt75*tt106*tt479*tt253
hes(1,10) = tt499+tt12*tt483*tt106*tt268-tt12*tt497*tt111*tt268+t&
&t12*tt111*tt249*tt498*tt253-tt12*tt75*tt106*tt498*tt253
hes(1,11) = tt508+tt12*tt500*tt106*tt268+tt507-tt12*tt504*tt111*t&
&t268+tt12*tt111*tt249*tt505*tt253-tt12*tt75*tt106*tt505*tt253
hes(1,12) = tt512+tt509*tt12*tt106*tt268+tt511-tt12*tt248*tt111*t&
&t268+tt12*tt111*tt249*tt510*tt253-tt12*tt75*tt106*tt510*tt253
hes(2,1) = (-tt12*tt257*tt106*tt268)+tt12*tt111*tt266*tt268+tt12*&
&tt257*tt249*tt250*tt253-tt12*tt75*tt266*tt250*tt253
hes(2,2) = tt12*tt257*tt249*tt267*tt253-tt12*tt75*tt266*tt267*tt2&
&53
hes(2,3) = tt12*tt269*tt266*tt268-tt12*tt270*tt257*tt268+tt12*tt2&
&57*tt249*tt271*tt253-tt12*tt75*tt266*tt271*tt253
hes(2,4) = (tt235*tt337*tt257*tt249*tt268)/2.0E+0+tt514-tt12*tt25&
&7*tt327*tt268-(tt235*tt337*tt75*tt266*tt268)/2.0E+0+tt12*tt276*tt&
&266*tt268+tt12*tt75*tt513*tt268+tt12*tt257*tt249*tt329*tt253-tt12&
&*tt75*tt266*tt329*tt253
hes(2,5) = (tt347*tt337*tt257*tt249*tt268)/2.0E+0+tt519-(tt347*tt&
&337*tt75*tt266*tt268)/2.0E+0+tt12*tt350*tt266*tt268-tt12*tt257*tt&
&376*tt268+tt12*tt516*tt75*tt268+tt12*tt257*tt249*tt377*tt253-tt12&
&*tt75*tt266*tt377*tt253
hes(2,6) = (tt396*tt337*tt257*tt249*tt268)/2.0E+0+tt521-(tt396*tt&
&337*tt75*tt266*tt268)/2.0E+0+tt12*tt389*tt266*tt268+tt12*tt520*tt&
&75*tt268-tt12*tt257*tt393*tt268+tt12*tt257*tt249*tt397*tt253-tt12&
&*tt75*tt266*tt397*tt253
hes(2,7) = (tt141*tt337*tt257*tt249*tt268)/2.0E+0+tt524-tt12*tt25&
&7*tt436*tt268-(tt141*tt337*tt75*tt266*tt268)/2.0E+0+tt12*tt406*tt&
&266*tt268+tt12*tt75*tt523*tt268+tt12*tt257*tt249*tt437*tt253-tt12&
&*tt75*tt266*tt437*tt253
hes(2,8) = (tt104*tt337*tt257*tt249*tt268)/2.0E+0+tt528-(tt104*tt&
&337*tt75*tt266*tt268)/2.0E+0+tt12*tt447*tt266*tt268-tt12*tt257*tt&
&464*tt268+tt12*tt526*tt75*tt268+tt12*tt257*tt249*tt465*tt253-tt12&
&*tt75*tt266*tt465*tt253
hes(2,9) = (tt478*tt337*tt257*tt249*tt268)/2.0E+0+tt530-(tt478*tt&
&337*tt75*tt266*tt268)/2.0E+0+tt12*tt474*tt266*tt268+tt12*tt529*tt&
&75*tt268-tt12*tt257*tt476*tt268+tt12*tt257*tt249*tt479*tt253-tt12&
&*tt75*tt266*tt479*tt253
hes(2,10) = tt508+tt12*tt483*tt266*tt268+tt531-tt12*tt497*tt257*t&
&t268+tt12*tt257*tt249*tt498*tt253-tt12*tt75*tt266*tt498*tt253
hes(2,11) = tt532+tt12*tt500*tt266*tt268-tt12*tt504*tt257*tt268+t&
&t12*tt257*tt249*tt505*tt253-tt12*tt75*tt266*tt505*tt253
hes(2,12) = tt534+tt509*tt12*tt266*tt268+tt533-tt12*tt248*tt257*t&
&t268+tt12*tt257*tt249*tt510*tt253-tt12*tt75*tt266*tt510*tt253
hes(3,1) = (-tt12*tt269*tt106*tt268)+tt12*tt270*tt111*tt268+tt12*&
&tt269*tt249*tt250*tt253-tt12*tt270*tt75*tt250*tt253
hes(3,2) = (-tt12*tt269*tt266*tt268)+tt12*tt270*tt257*tt268+tt12*&
&tt269*tt249*tt267*tt253-tt12*tt270*tt75*tt267*tt253
hes(3,3) = tt12*tt269*tt249*tt271*tt253-tt12*tt270*tt75*tt271*tt2&
&53
hes(3,4) = (tt235*tt337*tt269*tt249*tt268)/2.0E+0+tt536-tt12*tt26&
&9*tt327*tt268-(tt235*tt337*tt270*tt75*tt268)/2.0E+0+tt12*tt535*tt&
&75*tt268+tt12*tt270*tt276*tt268+tt12*tt269*tt249*tt329*tt253-tt12&
&*tt270*tt75*tt329*tt253
hes(3,5) = (tt347*tt337*tt269*tt249*tt268)/2.0E+0+tt538-tt12*tt26&
&9*tt376*tt268-(tt347*tt337*tt270*tt75*tt268)/2.0E+0+tt12*tt537*tt&
&75*tt268+tt12*tt270*tt350*tt268+tt12*tt269*tt249*tt377*tt253-tt12&
&*tt270*tt75*tt377*tt253
hes(3,6) = (tt396*tt337*tt269*tt249*tt268)/2.0E+0+tt540+tt539*tt1&
&2*tt75*tt268-(tt396*tt337*tt270*tt75*tt268)/2.0E+0-tt12*tt269*tt3&
&93*tt268+tt12*tt270*tt389*tt268+tt12*tt269*tt249*tt397*tt253-tt12&
&*tt270*tt75*tt397*tt253
hes(3,7) = (tt141*tt337*tt269*tt249*tt268)/2.0E+0+tt542-tt12*tt26&
&9*tt436*tt268-(tt141*tt337*tt270*tt75*tt268)/2.0E+0+tt12*tt541*tt&
&75*tt268+tt12*tt270*tt406*tt268+tt12*tt269*tt249*tt437*tt253-tt12&
&*tt270*tt75*tt437*tt253
hes(3,8) = (tt104*tt337*tt269*tt249*tt268)/2.0E+0+tt544-tt12*tt26&
&9*tt464*tt268-(tt104*tt337*tt270*tt75*tt268)/2.0E+0+tt12*tt543*tt&
&75*tt268+tt12*tt270*tt447*tt268+tt12*tt269*tt249*tt465*tt253-tt12&
&*tt270*tt75*tt465*tt253
hes(3,9) = (tt478*tt337*tt269*tt249*tt268)/2.0E+0+tt546+tt545*tt1&
&2*tt75*tt268-(tt478*tt337*tt270*tt75*tt268)/2.0E+0-tt12*tt269*tt4&
&76*tt268+tt12*tt270*tt474*tt268+tt12*tt269*tt249*tt479*tt253-tt12&
&*tt270*tt75*tt479*tt253
hes(3,10) = tt512+tt547-tt12*tt497*tt269*tt268+tt12*tt483*tt270*t&
&t268+tt12*tt269*tt249*tt498*tt253-tt12*tt270*tt75*tt498*tt253
hes(3,11) = tt534+tt548-tt12*tt504*tt269*tt268+tt12*tt500*tt270*t&
&t268+tt12*tt269*tt249*tt505*tt253-tt12*tt270*tt75*tt505*tt253
hes(3,12) = tt549-tt12*tt248*tt269*tt268+tt509*tt12*tt270*tt268+t&
&t12*tt269*tt249*tt510*tt253-tt12*tt270*tt75*tt510*tt253
hes(4,1) = tt111*tt550*tt268+tt75*(tt12*tt336-(tt235*tt337*tt106)&
&/2.0E+0)*tt268+tt343-tt12*tt276*tt106*tt268-tt75*tt550*tt250*tt25&
&3+tt12*tt276*tt249*tt250*tt253
hes(4,2) = tt257*tt550*tt268+tt75*(tt12*tt513-(tt235*tt337*tt266)&
&/2.0E+0)*tt268+tt514-tt12*tt276*tt266*tt268-tt75*tt550*tt267*tt25&
&3+tt12*tt276*tt249*tt267*tt253
hes(4,3) = tt269*tt550*tt268+tt536+(tt12*tt535-(tt235*tt337*tt270&
&)/2.0E+0)*tt75*tt268-tt12*tt270*tt276*tt268-tt75*tt271*tt550*tt25&
&3+tt12*tt276*tt249*tt271*tt253
hes(4,4) = tt276*tt550*tt268+tt75*(tt568+(3.0E+0*tt235**2*tt567*t&
&t249)/4.0E+0-tt235*tt337*tt327+tt12*(((tt563+tt558+tt32)*X(2,3)+t&
&t454+tt566)*X(3,4)+((tt557+tt556+tt172)*X(2,4)+tt161+(tt108+tt563&
&+tt558)*X(2,1))*X(3,3)+(tt245*X(2,4)+tt562+tt209*X(2,1))*X(3,2)+(&
&tt561+(tt159+tt557+tt556)*X(2,3)+tt525)*X(3,1)))*tt268+(tt235*tt3&
&37*tt276*tt249*tt268)/2.0E+0-tt12*(tt555+tt553+tt552-tt325*X(2,4)&
&+tt378+tt551)*tt249*tt268-tt12*tt276*tt327*tt268-tt75*tt550*tt329&
&*tt253+tt12*tt276*tt249*tt329*tt253
hes(4,5) = tt350*tt550*tt268+tt572+(tt347*tt337*tt276*tt249*tt268&
&)/2.0E+0+tt569-tt12*tt276*tt376*tt268-tt75*tt550*tt377*tt253+tt12&
&*tt276*tt249*tt377*tt253
hes(4,6) = tt389*tt550*tt268+tt574+(tt396*tt337*tt276*tt249*tt268&
&)/2.0E+0+tt573-tt12*tt276*tt393*tt268-tt75*tt550*tt397*tt253+tt12&
&*tt276*tt249*tt397*tt253
hes(4,7) = tt406*tt550*tt268+tt582+(tt141*tt337*tt276*tt249*tt268&
&)/2.0E+0+tt580-tt12*tt276*tt436*tt268-tt75*tt550*tt437*tt253+tt12&
&*tt276*tt249*tt437*tt253
hes(4,8) = tt447*tt550*tt268+tt586+(tt104*tt337*tt276*tt249*tt268&
&)/2.0E+0+tt583-tt12*tt276*tt464*tt268-tt75*tt550*tt465*tt253+tt12&
&*tt276*tt249*tt465*tt253
hes(4,9) = tt474*tt550*tt268+tt588+(tt478*tt337*tt276*tt249*tt268&
&)/2.0E+0+tt587-tt12*tt276*tt476*tt268-tt75*tt550*tt479*tt253+tt12&
&*tt276*tt249*tt479*tt253
hes(4,10) = tt483*tt550*tt268+tt593+(tt589*tt12-(tt235*tt337*tt49&
&7)/2.0E+0)*tt75*tt268-tt12*tt497*tt276*tt268-tt75*tt498*tt550*tt2&
&53+tt12*tt276*tt249*tt498*tt253
hes(4,11) = tt500*tt550*tt268+tt597+(tt12*tt596-(tt235*tt337*tt50&
&4)/2.0E+0)*tt75*tt268-tt12*tt504*tt276*tt268-tt75*tt505*tt550*tt2&
&53+tt12*tt276*tt249*tt505*tt253
hes(4,12) = tt509*tt550*tt268+tt598+(tt12*tt326-(tt235*tt337*tt24&
&8)/2.0E+0)*tt75*tt268-tt12*tt248*tt276*tt268-tt75*tt510*tt550*tt2&
&53+tt12*tt276*tt249*tt510*tt253
hes(5,1) = tt111*tt599*tt268+tt75*(tt12*tt387-(tt347*tt337*tt106)&
&/2.0E+0)*tt268+tt388-tt12*tt350*tt106*tt268-tt75*tt599*tt250*tt25&
&3+tt12*tt350*tt249*tt250*tt253
hes(5,2) = tt257*tt599*tt268+tt75*(tt12*tt516-(tt347*tt337*tt266)&
&/2.0E+0)*tt268+tt519-tt12*tt350*tt266*tt268-tt75*tt599*tt267*tt25&
&3+tt12*tt350*tt249*tt267*tt253
hes(5,3) = tt269*tt599*tt268+tt538+(tt12*tt537-(tt347*tt337*tt270&
&)/2.0E+0)*tt75*tt268-tt12*tt270*tt350*tt268-tt75*tt271*tt599*tt25&
&3+tt12*tt350*tt249*tt271*tt253
hes(5,4) = tt276*tt599*tt268+tt572+(tt235*tt337*tt350*tt249*tt268&
&)/2.0E+0+tt569-tt12*tt350*tt327*tt268-tt75*tt599*tt329*tt253+tt12&
&*tt350*tt249*tt329*tt253
hes(5,5) = tt350*tt599*tt268+tt75*(tt568+(3.0E+0*tt347**2*tt567*t&
&t249)/4.0E+0-tt347*tt337*tt376+tt12*((2*tt237*X(2,3)+6*tt21*X(2,2&
&)+2*tt16*X(2,1))*X(3,4)+(tt602+2*tt193*X(2,3)+6*tt189*X(2,2)+2*tt&
&187*X(2,1))*X(3,3)+(tt601+2*tt155*X(2,3)+2*tt107*X(2,1))*X(3,2)+(&
&tt600+2*tt128*X(2,3)+6*tt119*X(2,2))*X(3,1)))*tt268+(tt347*tt337*&
&tt350*tt249*tt268)/2.0E+0-tt12*(tt555+tt553+tt552-2*tt23)*tt249*t&
&t268-tt12*tt350*tt376*tt268-tt75*tt599*tt377*tt253+tt12*tt350*tt2&
&49*tt377*tt253
hes(5,6) = tt389*tt599*tt268+tt604+(tt396*tt337*tt350*tt249*tt268&
&)/2.0E+0+tt603-tt12*tt350*tt393*tt268-tt75*tt599*tt397*tt253+tt12&
&*tt350*tt249*tt397*tt253
hes(5,7) = tt406*tt599*tt268+tt606+(tt141*tt337*tt350*tt249*tt268&
&)/2.0E+0+tt605-tt12*tt350*tt436*tt268-tt75*tt599*tt437*tt253+tt12&
&*tt350*tt249*tt437*tt253
hes(5,8) = tt447*tt599*tt268+tt610+(tt104*tt337*tt350*tt249*tt268&
&)/2.0E+0+tt609-tt12*tt350*tt464*tt268-tt75*tt599*tt465*tt253+tt12&
&*tt350*tt249*tt465*tt253
hes(5,9) = tt474*tt599*tt268+tt612+(tt478*tt337*tt350*tt249*tt268&
&)/2.0E+0+tt611-tt12*tt350*tt476*tt268-tt75*tt599*tt479*tt253+tt12&
&*tt350*tt249*tt479*tt253
hes(5,10) = tt483*tt599*tt268+tt615+(tt12*tt614-(tt347*tt337*tt49&
&7)/2.0E+0)*tt75*tt268-tt12*tt497*tt350*tt268-tt75*tt498*tt599*tt2&
&53+tt12*tt350*tt249*tt498*tt253
hes(5,11) = tt500*tt599*tt268+tt618+(tt616*tt12-(tt347*tt337*tt50&
&4)/2.0E+0)*tt75*tt268-tt12*tt504*tt350*tt268-tt75*tt505*tt599*tt2&
&53+tt12*tt350*tt249*tt505*tt253
hes(5,12) = tt509*tt599*tt268+tt619+(tt12*tt375-(tt347*tt337*tt24&
&8)/2.0E+0)*tt75*tt268-tt12*tt248*tt350*tt268-tt75*tt510*tt599*tt2&
&53+tt12*tt350*tt249*tt510*tt253
hes(6,1) = tt111*tt620*tt268+tt75*(tt12*tt398-(tt396*tt337*tt106)&
&/2.0E+0)*tt268+tt399-tt12*tt389*tt106*tt268-tt75*tt620*tt250*tt25&
&3+tt12*tt389*tt249*tt250*tt253
hes(6,2) = tt257*tt620*tt268+tt75*(tt12*tt520-(tt396*tt337*tt266)&
&/2.0E+0)*tt268+tt521-tt12*tt389*tt266*tt268-tt75*tt620*tt267*tt25&
&3+tt12*tt389*tt249*tt267*tt253
hes(6,3) = tt269*tt620*tt268+tt540+(tt539*tt12-(tt396*tt337*tt270&
&)/2.0E+0)*tt75*tt268-tt12*tt270*tt389*tt268-tt75*tt271*tt620*tt25&
&3+tt12*tt389*tt249*tt271*tt253
hes(6,4) = tt276*tt620*tt268+tt574+(tt235*tt337*tt389*tt249*tt268&
&)/2.0E+0+tt573-tt12*tt389*tt327*tt268-tt75*tt620*tt329*tt253+tt12&
&*tt389*tt249*tt329*tt253
hes(6,5) = tt350*tt620*tt268+tt604+(tt347*tt337*tt389*tt249*tt268&
&)/2.0E+0+tt603-tt12*tt389*tt376*tt268-tt75*tt620*tt377*tt253+tt12&
&*tt389*tt249*tt377*tt253
hes(6,6) = tt389*tt620*tt268+tt75*(tt568+(3.0E+0*tt396**2*tt567*t&
&t249)/4.0E+0-tt396*tt337*tt393+tt12*(tt621+2*tt218*X(3,3)+6*tt180&
&*X(3,2)+2*tt178*X(3,1)))*tt268+(tt396*tt337*tt389*tt249*tt268)/2.&
&0E+0+2*tt56*tt12*tt249*tt268-tt12*tt389*tt393*tt268-tt75*tt620*tt&
&397*tt253+tt12*tt389*tt249*tt397*tt253
hes(6,7) = tt406*tt620*tt268+tt623+(tt141*tt337*tt389*tt249*tt268&
&)/2.0E+0+tt622-tt12*tt389*tt436*tt268-tt75*tt620*tt437*tt253+tt12&
&*tt389*tt249*tt437*tt253
hes(6,8) = tt447*tt620*tt268+tt625+(tt104*tt337*tt389*tt249*tt268&
&)/2.0E+0+tt624-tt12*tt389*tt464*tt268-tt75*tt620*tt465*tt253+tt12&
&*tt389*tt249*tt465*tt253
hes(6,9) = tt474*tt620*tt268+tt627+(tt478*tt337*tt389*tt249*tt268&
&)/2.0E+0+tt626-tt12*tt389*tt476*tt268-tt75*tt620*tt479*tt253+tt12&
&*tt389*tt249*tt479*tt253
hes(6,10) = tt483*tt620*tt268+tt629+(tt12*tt628-(tt396*tt337*tt49&
&7)/2.0E+0)*tt75*tt268-tt12*tt497*tt389*tt268-tt75*tt498*tt620*tt2&
&53+tt12*tt389*tt249*tt498*tt253
hes(6,11) = tt500*tt620*tt268+tt631+(tt12*tt630-(tt396*tt337*tt50&
&4)/2.0E+0)*tt75*tt268-tt12*tt504*tt389*tt268-tt75*tt505*tt620*tt2&
&53+tt12*tt389*tt249*tt505*tt253
hes(6,12) = tt509*tt620*tt268+tt632+(tt392*tt12-(tt396*tt337*tt24&
&8)/2.0E+0)*tt75*tt268-tt12*tt248*tt389*tt268-tt75*tt510*tt620*tt2&
&53+tt12*tt389*tt249*tt510*tt253
hes(7,1) = tt111*tt633*tt268+tt75*(tt12*tt439-(tt141*tt337*tt106)&
&/2.0E+0)*tt268+tt442-tt12*tt406*tt106*tt268-tt75*tt633*tt250*tt25&
&3+tt12*tt406*tt249*tt250*tt253
hes(7,2) = tt257*tt633*tt268+tt75*(tt12*tt523-(tt141*tt337*tt266)&
&/2.0E+0)*tt268+tt524-tt12*tt406*tt266*tt268-tt75*tt633*tt267*tt25&
&3+tt12*tt406*tt249*tt267*tt253
hes(7,3) = tt269*tt633*tt268+tt542+(tt12*tt541-(tt141*tt337*tt270&
&)/2.0E+0)*tt75*tt268-tt12*tt270*tt406*tt268-tt75*tt271*tt633*tt25&
&3+tt12*tt406*tt249*tt271*tt253
hes(7,4) = tt276*tt633*tt268+tt582+(tt235*tt337*tt406*tt249*tt268&
&)/2.0E+0+tt580-tt12*tt406*tt327*tt268-tt75*tt633*tt329*tt253+tt12&
&*tt406*tt249*tt329*tt253
hes(7,5) = tt350*tt633*tt268+tt606+(tt347*tt337*tt406*tt249*tt268&
&)/2.0E+0+tt605-tt12*tt406*tt376*tt268-tt75*tt633*tt377*tt253+tt12&
&*tt406*tt249*tt377*tt253
hes(7,6) = tt389*tt633*tt268+tt623+(tt396*tt337*tt406*tt249*tt268&
&)/2.0E+0+tt622-tt12*tt406*tt393*tt268-tt75*tt633*tt397*tt253+tt12&
&*tt406*tt249*tt397*tt253
hes(7,7) = tt406*tt633*tt268+tt75*(tt568+(3.0E+0*tt141**2*tt567*t&
&t249)/4.0E+0-tt141*tt337*tt436+tt12*((2*tt46*X(2,3)+(tt559+2*tt23&
&0)*X(2,2)+tt566)*X(3,4)+(tt602+tt161+tt210*X(2,1))*X(3,3)+((tt564&
&+2*tt168)*X(2,4)+tt562+(tt159+tt559+tt305)*X(2,1))*X(3,2)+(tt561+&
&tt515+(tt108+tt564+tt319)*X(2,2))*X(3,1)))*tt268+(tt141*tt337*tt4&
&06*tt249*tt268)/2.0E+0-tt12*(tt634+tt467+tt585-(tt61+tt324)*X(2,4&
&)+tt408+tt321)*tt249*tt268-tt12*tt406*tt436*tt268-tt75*tt633*tt43&
&7*tt253+tt12*tt406*tt249*tt437*tt253
hes(7,8) = tt447*tt633*tt268+tt636+(tt104*tt337*tt406*tt249*tt268&
&)/2.0E+0+tt635-tt12*tt406*tt464*tt268-tt75*tt633*tt465*tt253+tt12&
&*tt406*tt249*tt465*tt253
hes(7,9) = tt474*tt633*tt268+tt638+(tt478*tt337*tt406*tt249*tt268&
&)/2.0E+0+tt637-tt12*tt406*tt476*tt268-tt75*tt633*tt479*tt253+tt12&
&*tt406*tt249*tt479*tt253
hes(7,10) = tt483*tt633*tt268+tt643+(tt639*tt12-(tt141*tt337*tt49&
&7)/2.0E+0)*tt75*tt268-tt12*tt497*tt406*tt268-tt75*tt498*tt633*tt2&
&53+tt12*tt406*tt249*tt498*tt253
hes(7,11) = tt500*tt633*tt268+tt645+(tt12*tt644-(tt141*tt337*tt50&
&4)/2.0E+0)*tt75*tt268-tt12*tt504*tt406*tt268-tt75*tt505*tt633*tt2&
&53+tt12*tt406*tt249*tt505*tt253
hes(7,12) = tt509*tt633*tt268+tt646+(tt12*tt435-(tt141*tt337*tt24&
&8)/2.0E+0)*tt75*tt268-tt12*tt248*tt406*tt268-tt75*tt510*tt633*tt2&
&53+tt12*tt406*tt249*tt510*tt253
hes(8,1) = tt111*tt647*tt268+tt75*(tt12*tt472-(tt104*tt337*tt106)&
&/2.0E+0)*tt268+tt473-tt12*tt447*tt106*tt268-tt75*tt647*tt250*tt25&
&3+tt12*tt447*tt249*tt250*tt253
hes(8,2) = tt257*tt647*tt268+tt75*(tt12*tt526-(tt104*tt337*tt266)&
&/2.0E+0)*tt268+tt528-tt12*tt447*tt266*tt268-tt75*tt647*tt267*tt25&
&3+tt12*tt447*tt249*tt267*tt253
hes(8,3) = tt269*tt647*tt268+tt544+(tt12*tt543-(tt104*tt337*tt270&
&)/2.0E+0)*tt75*tt268-tt12*tt270*tt447*tt268-tt75*tt271*tt647*tt25&
&3+tt12*tt447*tt249*tt271*tt253
hes(8,4) = tt276*tt647*tt268+tt586+(tt235*tt337*tt447*tt249*tt268&
&)/2.0E+0+tt583-tt12*tt447*tt327*tt268-tt75*tt647*tt329*tt253+tt12&
&*tt447*tt249*tt329*tt253
hes(8,5) = tt350*tt647*tt268+tt610+(tt347*tt337*tt447*tt249*tt268&
&)/2.0E+0+tt609-tt12*tt447*tt376*tt268-tt75*tt647*tt377*tt253+tt12&
&*tt447*tt249*tt377*tt253
hes(8,6) = tt389*tt647*tt268+tt625+(tt396*tt337*tt447*tt249*tt268&
&)/2.0E+0+tt624-tt12*tt447*tt393*tt268-tt75*tt647*tt397*tt253+tt12&
&*tt447*tt249*tt397*tt253
hes(8,7) = tt406*tt647*tt268+tt636+(tt141*tt337*tt447*tt249*tt268&
&)/2.0E+0+tt635-tt12*tt447*tt436*tt268-tt75*tt647*tt437*tt253+tt12&
&*tt447*tt249*tt437*tt253
hes(8,8) = tt447*tt647*tt268+tt75*(tt568+(3.0E+0*tt104**2*tt567*t&
&t249)/4.0E+0-tt104*tt337*tt464+tt12*((6*tt46*X(2,3)+2*tt241)*X(3,&
&4)+(tt602+2*tt198)*X(3,3)+(tt601+6*tt155*X(2,3)+2*tt162)*X(3,2)+(&
&tt600+6*tt110*X(2,3)+2*tt133*X(2,2))*X(3,1)))*tt268+(tt104*tt337*&
&tt447*tt249*tt268)/2.0E+0-tt12*(tt634+tt467+tt585-2*tt39)*tt249*t&
&t268-tt12*tt447*tt464*tt268-tt75*tt647*tt465*tt253+tt12*tt447*tt2&
&49*tt465*tt253
hes(8,9) = tt474*tt647*tt268+tt649+(tt478*tt337*tt447*tt249*tt268&
&)/2.0E+0+tt648-tt12*tt447*tt476*tt268-tt75*tt647*tt479*tt253+tt12&
&*tt447*tt249*tt479*tt253
hes(8,10) = tt483*tt647*tt268+tt651+(tt12*tt650-(tt104*tt337*tt49&
&7)/2.0E+0)*tt75*tt268-tt12*tt497*tt447*tt268-tt75*tt498*tt647*tt2&
&53+tt12*tt447*tt249*tt498*tt253
hes(8,11) = tt500*tt647*tt268+tt654+(tt652*tt12-(tt104*tt337*tt50&
&4)/2.0E+0)*tt75*tt268-tt12*tt504*tt447*tt268-tt75*tt505*tt647*tt2&
&53+tt12*tt447*tt249*tt505*tt253
hes(8,12) = tt509*tt647*tt268+tt655+(tt12*tt463-(tt104*tt337*tt24&
&8)/2.0E+0)*tt75*tt268-tt12*tt248*tt447*tt268-tt75*tt510*tt647*tt2&
&53+tt12*tt447*tt249*tt510*tt253
hes(9,1) = tt111*tt656*tt268+tt75*(tt12*tt480-(tt478*tt337*tt106)&
&/2.0E+0)*tt268+tt481-tt12*tt474*tt106*tt268-tt75*tt656*tt250*tt25&
&3+tt12*tt474*tt249*tt250*tt253
hes(9,2) = tt257*tt656*tt268+tt75*(tt12*tt529-(tt478*tt337*tt266)&
&/2.0E+0)*tt268+tt530-tt12*tt474*tt266*tt268-tt75*tt656*tt267*tt25&
&3+tt12*tt474*tt249*tt267*tt253
hes(9,3) = tt269*tt656*tt268+tt546+(tt545*tt12-(tt478*tt337*tt270&
&)/2.0E+0)*tt75*tt268-tt12*tt270*tt474*tt268-tt75*tt271*tt656*tt25&
&3+tt12*tt474*tt249*tt271*tt253
hes(9,4) = tt276*tt656*tt268+tt588+(tt235*tt337*tt474*tt249*tt268&
&)/2.0E+0+tt587-tt12*tt474*tt327*tt268-tt75*tt656*tt329*tt253+tt12&
&*tt474*tt249*tt329*tt253
hes(9,5) = tt350*tt656*tt268+tt612+(tt347*tt337*tt474*tt249*tt268&
&)/2.0E+0+tt611-tt12*tt474*tt376*tt268-tt75*tt656*tt377*tt253+tt12&
&*tt474*tt249*tt377*tt253
hes(9,6) = tt389*tt656*tt268+tt627+(tt396*tt337*tt474*tt249*tt268&
&)/2.0E+0+tt626-tt12*tt474*tt393*tt268-tt75*tt656*tt397*tt253+tt12&
&*tt474*tt249*tt397*tt253
hes(9,7) = tt406*tt656*tt268+tt638+(tt141*tt337*tt474*tt249*tt268&
&)/2.0E+0+tt637-tt12*tt474*tt436*tt268-tt75*tt656*tt437*tt253+tt12&
&*tt474*tt249*tt437*tt253
hes(9,8) = tt447*tt656*tt268+tt649+(tt104*tt337*tt474*tt249*tt268&
&)/2.0E+0+tt648-tt12*tt474*tt464*tt268-tt75*tt656*tt465*tt253+tt12&
&*tt474*tt249*tt465*tt253
hes(9,9) = tt474*tt656*tt268+tt75*(tt568+(3.0E+0*tt478**2*tt567*t&
&t249)/4.0E+0-tt478*tt337*tt476+tt12*(tt621+6*tt226*X(3,3)+2*tt225&
&))*tt268+(tt478*tt337*tt474*tt249*tt268)/2.0E+0+2*tt68*tt12*tt249&
&*tt268-tt12*tt474*tt476*tt268-tt75*tt656*tt479*tt253+tt12*tt474*t&
&t249*tt479*tt253
hes(9,10) = tt483*tt656*tt268+tt658+(tt12*tt657-(tt478*tt337*tt49&
&7)/2.0E+0)*tt75*tt268-tt12*tt497*tt474*tt268-tt75*tt498*tt656*tt2&
&53+tt12*tt474*tt249*tt498*tt253
hes(9,11) = tt500*tt656*tt268+tt660+(tt12*tt659-(tt478*tt337*tt50&
&4)/2.0E+0)*tt75*tt268-tt12*tt504*tt474*tt268-tt75*tt505*tt656*tt2&
&53+tt12*tt474*tt249*tt505*tt253
hes(9,12) = tt509*tt656*tt268+tt661+(tt475*tt12-(tt478*tt337*tt24&
&8)/2.0E+0)*tt75*tt268-tt12*tt248*tt474*tt268-tt75*tt510*tt656*tt2&
&53+tt12*tt474*tt249*tt510*tt253
hes(10,1) = tt499-tt12*tt483*tt106*tt268+tt12*tt497*tt111*tt268+t&
&t12*tt483*tt249*tt250*tt253-tt12*tt497*tt75*tt250*tt253
hes(10,2) = tt508-tt12*tt483*tt266*tt268+tt531+tt12*tt497*tt257*t&
&t268+tt12*tt483*tt249*tt267*tt253-tt12*tt497*tt75*tt267*tt253
hes(10,3) = tt512+tt547+tt12*tt497*tt269*tt268-tt12*tt483*tt270*t&
&t268+tt12*tt483*tt249*tt271*tt253-tt12*tt497*tt75*tt271*tt253
hes(10,4) = (tt235*tt337*tt483*tt249*tt268)/2.0E+0+tt593-tt12*tt4&
&83*tt327*tt268+tt589*tt12*tt75*tt268-(tt235*tt337*tt497*tt75*tt26&
&8)/2.0E+0+tt12*tt497*tt276*tt268+tt12*tt483*tt249*tt329*tt253-tt1&
&2*tt497*tt75*tt329*tt253
hes(10,5) = tt615+(tt347*tt337*tt483*tt249*tt268)/2.0E+0-tt12*tt4&
&83*tt376*tt268-(tt347*tt337*tt497*tt75*tt268)/2.0E+0+tt12*tt614*t&
&t75*tt268+tt12*tt497*tt350*tt268+tt12*tt483*tt249*tt377*tt253-tt1&
&2*tt497*tt75*tt377*tt253
hes(10,6) = tt629+(tt396*tt337*tt483*tt249*tt268)/2.0E+0-(tt396*t&
&t337*tt497*tt75*tt268)/2.0E+0+tt12*tt628*tt75*tt268-tt12*tt483*tt&
&393*tt268+tt12*tt497*tt389*tt268+tt12*tt483*tt249*tt397*tt253-tt1&
&2*tt497*tt75*tt397*tt253
hes(10,7) = tt643+(tt141*tt337*tt483*tt249*tt268)/2.0E+0-tt12*tt4&
&83*tt436*tt268+tt639*tt12*tt75*tt268-(tt141*tt337*tt497*tt75*tt26&
&8)/2.0E+0+tt12*tt497*tt406*tt268+tt12*tt483*tt249*tt437*tt253-tt1&
&2*tt497*tt75*tt437*tt253
hes(10,8) = tt651+(tt104*tt337*tt483*tt249*tt268)/2.0E+0-tt12*tt4&
&83*tt464*tt268-(tt104*tt337*tt497*tt75*tt268)/2.0E+0+tt12*tt650*t&
&t75*tt268+tt12*tt497*tt447*tt268+tt12*tt483*tt249*tt465*tt253-tt1&
&2*tt497*tt75*tt465*tt253
hes(10,9) = tt658+(tt478*tt337*tt483*tt249*tt268)/2.0E+0-(tt478*t&
&t337*tt497*tt75*tt268)/2.0E+0+tt12*tt657*tt75*tt268-tt12*tt483*tt&
&476*tt268+tt12*tt497*tt474*tt268+tt12*tt483*tt249*tt479*tt253-tt1&
&2*tt497*tt75*tt479*tt253
hes(10,10) = tt12*tt483*tt249*tt498*tt253-tt12*tt497*tt75*tt498*t&
&t253
hes(10,11) = tt12*tt500*tt497*tt268-tt12*tt483*tt504*tt268+tt12*t&
&t483*tt249*tt505*tt253-tt12*tt497*tt75*tt505*tt253
hes(10,12) = tt509*tt12*tt497*tt268-tt12*tt483*tt248*tt268+tt12*t&
&t483*tt249*tt510*tt253-tt12*tt497*tt75*tt510*tt253
hes(11,1) = tt508-tt12*tt500*tt106*tt268+tt507+tt12*tt504*tt111*t&
&t268+tt12*tt500*tt249*tt250*tt253-tt12*tt504*tt75*tt250*tt253
hes(11,2) = tt532-tt12*tt500*tt266*tt268+tt12*tt504*tt257*tt268+t&
&t12*tt500*tt249*tt267*tt253-tt12*tt504*tt75*tt267*tt253
hes(11,3) = tt534+tt548+tt12*tt504*tt269*tt268-tt12*tt500*tt270*t&
&t268+tt12*tt500*tt249*tt271*tt253-tt12*tt504*tt75*tt271*tt253
hes(11,4) = tt597+(tt235*tt337*tt500*tt249*tt268)/2.0E+0-tt12*tt5&
&00*tt327*tt268-(tt235*tt337*tt504*tt75*tt268)/2.0E+0+tt12*tt596*t&
&t75*tt268+tt12*tt504*tt276*tt268+tt12*tt500*tt249*tt329*tt253-tt1&
&2*tt504*tt75*tt329*tt253
hes(11,5) = (tt347*tt337*tt500*tt249*tt268)/2.0E+0+tt618-tt12*tt5&
&00*tt376*tt268+tt616*tt12*tt75*tt268-(tt347*tt337*tt504*tt75*tt26&
&8)/2.0E+0+tt12*tt504*tt350*tt268+tt12*tt500*tt249*tt377*tt253-tt1&
&2*tt504*tt75*tt377*tt253
hes(11,6) = tt631+(tt396*tt337*tt500*tt249*tt268)/2.0E+0-(tt396*t&
&t337*tt504*tt75*tt268)/2.0E+0+tt12*tt630*tt75*tt268-tt12*tt500*tt&
&393*tt268+tt12*tt504*tt389*tt268+tt12*tt500*tt249*tt397*tt253-tt1&
&2*tt504*tt75*tt397*tt253
hes(11,7) = tt645+(tt141*tt337*tt500*tt249*tt268)/2.0E+0-tt12*tt5&
&00*tt436*tt268-(tt141*tt337*tt504*tt75*tt268)/2.0E+0+tt12*tt644*t&
&t75*tt268+tt12*tt504*tt406*tt268+tt12*tt500*tt249*tt437*tt253-tt1&
&2*tt504*tt75*tt437*tt253
hes(11,8) = tt654+(tt104*tt337*tt500*tt249*tt268)/2.0E+0-tt12*tt5&
&00*tt464*tt268+tt652*tt12*tt75*tt268-(tt104*tt337*tt504*tt75*tt26&
&8)/2.0E+0+tt12*tt504*tt447*tt268+tt12*tt500*tt249*tt465*tt253-tt1&
&2*tt504*tt75*tt465*tt253
hes(11,9) = tt660+(tt478*tt337*tt500*tt249*tt268)/2.0E+0-(tt478*t&
&t337*tt504*tt75*tt268)/2.0E+0+tt12*tt659*tt75*tt268-tt12*tt500*tt&
&476*tt268+tt12*tt504*tt474*tt268+tt12*tt500*tt249*tt479*tt253-tt1&
&2*tt504*tt75*tt479*tt253
hes(11,10) = (-tt12*tt500*tt497*tt268)+tt12*tt483*tt504*tt268+tt1&
&2*tt500*tt249*tt498*tt253-tt12*tt504*tt75*tt498*tt253
hes(11,11) = tt12*tt500*tt249*tt505*tt253-tt12*tt504*tt75*tt505*t&
&t253
hes(11,12) = tt509*tt12*tt504*tt268-tt12*tt500*tt248*tt268+tt12*t&
&t500*tt249*tt510*tt253-tt12*tt504*tt75*tt510*tt253
hes(12,1) = tt512-tt509*tt12*tt106*tt268+tt511+tt12*tt248*tt111*t&
&t268+tt509*tt12*tt249*tt250*tt253-tt12*tt248*tt75*tt250*tt253
hes(12,2) = tt534-tt509*tt12*tt266*tt268+tt533+tt12*tt248*tt257*t&
&t268+tt509*tt12*tt249*tt267*tt253-tt12*tt248*tt75*tt267*tt253
hes(12,3) = tt549+tt12*tt248*tt269*tt268-tt509*tt12*tt270*tt268+t&
&t509*tt12*tt249*tt271*tt253-tt12*tt248*tt75*tt271*tt253
hes(12,4) = tt598+(tt235*tt509*tt337*tt249*tt268)/2.0E+0-tt509*tt&
&12*tt327*tt268-(tt235*tt337*tt248*tt75*tt268)/2.0E+0+tt12*tt326*t&
&t75*tt268+tt12*tt248*tt276*tt268+tt509*tt12*tt249*tt329*tt253-tt1&
&2*tt248*tt75*tt329*tt253
hes(12,5) = tt619+(tt347*tt509*tt337*tt249*tt268)/2.0E+0-tt509*tt&
&12*tt376*tt268-(tt347*tt337*tt248*tt75*tt268)/2.0E+0+tt12*tt375*t&
&t75*tt268+tt12*tt248*tt350*tt268+tt509*tt12*tt249*tt377*tt253-tt1&
&2*tt248*tt75*tt377*tt253
hes(12,6) = tt632+(tt396*tt509*tt337*tt249*tt268)/2.0E+0+tt392*tt&
&12*tt75*tt268-(tt396*tt337*tt248*tt75*tt268)/2.0E+0-tt509*tt12*tt&
&393*tt268+tt12*tt248*tt389*tt268+tt509*tt12*tt249*tt397*tt253-tt1&
&2*tt248*tt75*tt397*tt253
hes(12,7) = tt646+(tt141*tt509*tt337*tt249*tt268)/2.0E+0-tt509*tt&
&12*tt436*tt268-(tt141*tt337*tt248*tt75*tt268)/2.0E+0+tt12*tt435*t&
&t75*tt268+tt12*tt248*tt406*tt268+tt509*tt12*tt249*tt437*tt253-tt1&
&2*tt248*tt75*tt437*tt253
hes(12,8) = tt655+(tt104*tt509*tt337*tt249*tt268)/2.0E+0-tt509*tt&
&12*tt464*tt268-(tt104*tt337*tt248*tt75*tt268)/2.0E+0+tt12*tt463*t&
&t75*tt268+tt12*tt248*tt447*tt268+tt509*tt12*tt249*tt465*tt253-tt1&
&2*tt248*tt75*tt465*tt253
hes(12,9) = tt661+(tt478*tt509*tt337*tt249*tt268)/2.0E+0+tt475*tt&
&12*tt75*tt268-(tt478*tt337*tt248*tt75*tt268)/2.0E+0-tt509*tt12*tt&
&476*tt268+tt12*tt248*tt474*tt268+tt509*tt12*tt249*tt479*tt253-tt1&
&2*tt248*tt75*tt479*tt253
hes(12,10) = (-tt509*tt12*tt497*tt268)+tt12*tt483*tt248*tt268+tt5&
&09*tt12*tt249*tt498*tt253-tt12*tt248*tt75*tt498*tt253
hes(12,11) = (-tt509*tt12*tt504*tt268)+tt12*tt500*tt248*tt268+tt5&
&09*tt12*tt249*tt505*tt253-tt12*tt248*tt75*tt505*tt253
hes(12,12) = tt509*tt12*tt249*tt510*tt253-tt12*tt248*tt75*tt510*t&
&t253
END 
SUBROUTINE surf_bending(val, X, d, l, area) 
IMPLICIT NONE 
REAL(KIND=8) val(1, 1) 
REAL(KIND=8) X(3, 4) 
REAL(KIND=8) d(1, 1) 
REAL(KIND=8) l(1, 1) 
REAL(KIND=8) area(1, 1) 
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
tt1 = X(1,2)**2
tt2 = -2*X(1,2)*X(1,3)
tt3 = X(1,3)**2
tt4 = X(2,2)**2
tt5 = -2*X(2,2)*X(2,3)
tt6 = X(2,3)**2
tt7 = X(3,2)**2
tt8 = X(3,3)**2
tt9 = X(1,3)**3
tt10 = -tt1
tt11 = -tt3
tt12 = (tt11+2*X(1,2)*X(1,3)+tt10)*X(1,4)
tt13 = -X(1,4)
tt14 = tt13+X(1,3)
tt15 = X(2,2)**3
tt16 = X(1,2)**3
tt17 = -tt16
tt18 = tt3+tt2+tt1
tt19 = tt18*X(1,4)
tt20 = -X(1,2)
tt21 = -2*X(1,3)
tt22 = 3*X(1,4)
tt23 = 2*X(1,2)
tt24 = -3*X(1,4)
tt25 = X(1,4)+tt20
tt26 = X(2,3)**3
tt27 = -tt9
tt28 = -X(1,3)
tt29 = tt28+X(1,2)
tt30 = -2*X(1,2)
tt31 = 2*X(1,3)
tt32 = tt31+tt30
tt33 = X(1,4)+tt28
tt34 = X(1,1)*tt1
tt35 = -2*X(1,4)
tt36 = tt13+X(1,1)
tt37 = tt33*X(2,1)
tt38 = -2*X(1,1)
tt39 = 2*X(1,4)
tt40 = tt39+tt38
tt41 = -X(1,1)*tt1
tt42 = 2*X(1,1)*X(1,2)
tt43 = -X(1,1)
tt44 = X(1,3)+tt43
tt45 = 2*X(1,1)
tt46 = tt29*X(2,4)+tt25*X(2,3)+tt14*X(2,2)
tt47 = tt36*X(2,3)
tt48 = tt13+X(1,2)
tt49 = X(1,4)+tt43
tt50 = tt35+tt45
tt51 = tt48*X(2,1)
tt52 = tt49*X(2,2)
tt53 = -2*X(1,1)*X(1,2)
tt54 = tt20+X(1,1)
tt55 = 3*X(1,1)
tt56 = -3*X(1,1)
tt57 = X(1,3)+tt20
tt58 = tt28+X(1,1)
tt59 = X(1,2)+tt43
tt60 = tt21+tt23
tt61 = tt57*X(2,1)
tt62 = tt59*X(2,3)+tt58*X(2,2)+tt61
tt63 = X(1,2)*X(1,3)
tt64 = tt57*X(1,4)
tt65 = -X(1,1)*X(1,3)
tt66 = tt58*X(1,4)
tt67 = tt29*X(1,4)
tt68 = X(1,1)*X(1,2)
tt69 = (tt30+X(1,1))*X(1,3)
tt70 = (X(1,3)+X(1,2)+tt38)*X(1,4)
tt71 = -X(1,1)*X(1,2)
tt72 = tt54*X(1,4)
tt73 = (X(1,2)+X(1,1))*X(1,3)
tt74 = tt59*X(1,3)
tt75 = X(2,2)*X(2,3)
tt76 = -tt6
tt77 = -X(2,2)
tt78 = -X(2,3)
tt79 = -tt4
tt80 = X(2,1)*X(2,2)
tt81 = -X(2,1)*X(2,2)
val(1,1) = (3*l(1,1)**2*(atan2(((tt62*tt8+((tt30+tt45)*X(2,3)+(tt&
&31+tt38)*X(2,2)+tt60*X(2,1))*X(3,2)*X(3,3)+tt62*tt7+tt59*tt26+((t&
&t28+tt30+tt55)*X(2,2)+tt61)*tt6+((tt31+X(1,2)+tt56)*tt4+tt60*X(2,&
&1)*X(2,2)+tt59*tt3+(tt42-2*tt1)*X(1,3)+tt16+tt41)*X(2,3)+tt58*tt1&
&5+tt57*X(2,1)*tt4+(tt27+(tt23+X(1,1))*tt3+(tt10+tt53)*X(1,3)+tt34&
&)*X(2,2)+(tt9-3*X(1,2)*tt3+3*tt1*X(1,3)+tt17)*X(2,1))*X(3,4)+(tt5&
&4*X(2,4)+tt52+tt51)*X(3,3)**3+(((X(1,3)+tt23+tt56)*X(2,4)+tt47+tt&
&50*X(2,2)+(tt22+tt28+tt30)*X(2,1))*X(3,2)+tt46*X(3,1))*tt8+(((tt2&
&1+tt20+tt55)*X(2,4)+tt40*X(2,3)+tt52+(tt24+tt31+X(1,2))*X(2,1))*t&
&t7+(tt32*X(2,4)+(tt35+tt23)*X(2,3)+(tt39+tt21)*X(2,2))*X(3,1)*X(3&
&,2)+(tt54*tt6+(tt23+tt38)*X(2,2)*X(2,3)+tt54*tt4+tt54*tt3+(2*tt1+&
&tt53)*X(1,3)+tt17+tt34)*X(2,4)+(tt52+tt51)*tt6+(tt50*tt4+(tt39+tt&
&30)*X(2,1)*X(2,2))*X(2,3)+tt49*tt15+tt48*X(2,1)*tt4+(tt19-X(1,1)*&
&tt3+2*X(1,1)*X(1,2)*X(1,3)+tt41)*X(2,2)+(tt12+X(1,2)*tt3-2*tt1*X(&
&1,3)+tt16)*X(2,1))*X(3,3)+(tt44*X(2,4)+tt47+tt37)*X(3,2)**3+tt46*&
&X(3,1)*tt7+((tt44*tt6+(tt21+tt45)*X(2,2)*X(2,3)+tt44*tt4+tt9+(tt3&
&0+tt43)*tt3+(tt1+tt42)*X(1,3)+tt41)*X(2,4)+tt36*tt26+(tt40*X(2,2)&
&+tt37)*tt6+(tt36*tt4+(tt35+tt31)*X(2,1)*X(2,2)+tt12+X(1,1)*tt3-2*&
&X(1,1)*X(1,2)*X(1,3)+tt34)*X(2,3)+tt33*X(2,1)*tt4+(tt19+tt27+2*X(&
&1,2)*tt3-tt1*X(1,3))*X(2,1))*X(3,2)+((tt29*tt6+tt32*X(2,2)*X(2,3)&
&+tt29*tt4+tt27+3*X(1,2)*tt3-3*tt1*X(1,3)+tt16)*X(2,4)+tt25*tt26+(&
&tt24+X(1,3)+tt23)*X(2,2)*tt6+((tt22+tt21+tt20)*tt4+tt19-X(1,2)*tt&
&3+2*tt1*X(1,3)+tt17)*X(2,3)+tt14*tt15+(tt12+tt9-2*X(1,2)*tt3+tt1*&
&X(1,3))*X(2,2))*X(3,1))/sqrt(tt8-2*X(3,2)*X(3,3)+tt7+tt6+tt5+tt4+&
&tt3+tt2+tt1),(-(((X(2,2)-X(2,1))*X(2,3)+tt79+tt80+tt74+tt10+tt68)&
&*X(3,3)+(tt76+(X(2,2)+X(2,1))*X(2,3)+tt81+tt11+tt73+tt71)*X(3,2)+&
&(tt6+tt5+tt4+tt3+tt2+tt1)*X(3,1))*X(3,4))-((tt77+X(2,1))*X(2,4)+t&
&t4+tt81+tt72+tt1+tt71)*tt8-(((X(2,3)+X(2,2)-2*X(2,1))*X(2,4)+(X(2&
&,1)-2*X(2,2))*X(2,3)+tt80+tt70+tt69+tt68)*X(3,2)+((tt78+X(2,2))*X&
&(2,4)+tt75+tt79+tt67+tt63+tt10)*X(3,1))*X(3,3)-((tt78+X(2,1))*X(2&
&,4)+tt6-X(2,1)*X(2,3)+tt66+tt3+tt65)*tt7-((X(2,3)+tt77)*X(2,4)+tt&
&76+tt75+tt64+tt11+tt63)*X(3,1)*X(3,2)-((tt74+tt10+tt68)*X(2,3)+(t&
&t11+tt73+tt71)*X(2,2)+tt18*X(2,1))*X(2,4)-(tt72+tt1+tt71)*tt6-((t&
&t70+tt69+tt68)*X(2,2)+(tt67+tt63+tt10)*X(2,1))*X(2,3)-(tt66+tt3+t&
&t65)*tt4-(tt64+tt11+tt63)*X(2,1)*X(2,2))-d(1,1))**2)/area(1,1)
END 
SUBROUTINE surf_bending_jac(jac, X, d, l, area) 
IMPLICIT NONE 
REAL(KIND=8) jac(1, 12) 
REAL(KIND=8) X(3, 4) 
REAL(KIND=8) d(1, 1) 
REAL(KIND=8) l(1, 1) 
REAL(KIND=8) area(1, 1) 
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
tt1 = 1/area(1,1)
tt2 = l(1,1)**2
tt3 = X(1,2)**2
tt4 = -2*X(1,2)*X(1,3)
tt5 = X(1,3)**2
tt6 = X(2,2)**2
tt7 = -2*X(2,2)*X(2,3)
tt8 = X(2,3)**2
tt9 = X(3,2)**2
tt10 = X(3,3)**2
tt11 = tt10-2*X(3,2)*X(3,3)+tt9+tt8+tt7+tt6+tt5+tt4+tt3
tt12 = sqrt(tt11)
tt13 = 1/tt12
tt14 = X(1,2)*X(1,3)
tt15 = -tt5
tt16 = -X(1,2)
tt17 = X(1,3)+tt16
tt18 = tt17*X(1,4)
tt19 = tt18+tt15+tt14
tt20 = -X(1,1)*X(1,3)
tt21 = -X(1,3)
tt22 = tt21+X(1,1)
tt23 = tt22*X(1,4)
tt24 = tt23+tt5+tt20
tt25 = -tt3
tt26 = tt21+X(1,2)
tt27 = tt26*X(1,4)
tt28 = tt27+tt14+tt25
tt29 = X(1,1)*X(1,2)
tt30 = -2*X(1,2)
tt31 = (tt30+X(1,1))*X(1,3)
tt32 = -2*X(1,1)
tt33 = X(1,3)+X(1,2)+tt32
tt34 = tt33*X(1,4)
tt35 = tt34+tt31+tt29
tt36 = -X(1,1)*X(1,2)
tt37 = tt16+X(1,1)
tt38 = tt37*X(1,4)
tt39 = tt38+tt3+tt36
tt40 = tt5+tt4+tt3
tt41 = tt40*X(2,1)
tt42 = (X(1,2)+X(1,1))*X(1,3)
tt43 = tt15+tt42+tt36
tt44 = -X(1,1)
tt45 = X(1,2)+tt44
tt46 = tt45*X(1,3)
tt47 = tt46+tt25+tt29
tt48 = X(2,2)*X(2,3)
tt49 = -tt8
tt50 = -X(2,2)
tt51 = X(2,3)+tt50
tt52 = tt51*X(2,4)+tt49+tt48+tt18+tt15+tt14
tt53 = -X(2,3)
tt54 = tt53+X(2,1)
tt55 = tt54*X(2,4)+tt8-X(2,1)*X(2,3)+tt23+tt5+tt20
tt56 = -tt6
tt57 = tt53+X(2,2)
tt58 = tt57*X(2,4)+tt48+tt56+tt27+tt14+tt25
tt59 = X(2,1)*X(2,2)
tt60 = -2*X(2,2)
tt61 = -2*X(2,1)
tt62 = X(2,3)+X(2,2)+tt61
tt63 = tt62*X(2,4)+(tt60+X(2,1))*X(2,3)+tt59+tt34+tt31+tt29
tt64 = -X(2,1)*X(2,2)
tt65 = tt50+X(2,1)
tt66 = tt65*X(2,4)+tt6+tt64+tt38+tt3+tt36
tt67 = tt8+tt7+tt6+tt5+tt4+tt3
tt68 = tt49+(X(2,2)+X(2,1))*X(2,3)+tt64+tt15+tt42+tt36
tt69 = -X(2,1)
tt70 = X(2,2)+tt69
tt71 = tt70*X(2,3)+tt56+tt59+tt46+tt25+tt29
tt72 = (-(tt71*X(3,3)+tt68*X(3,2)+tt67*X(3,1))*X(3,4))-tt66*tt10-&
&(tt63*X(3,2)+tt58*X(3,1))*X(3,3)-tt55*tt9-tt52*X(3,1)*X(3,2)-(tt4&
&7*X(2,3)+tt43*X(2,2)+tt41)*X(2,4)-tt39*tt8-(tt35*X(2,2)+tt28*X(2,&
&1))*X(2,3)-tt24*tt6-tt19*X(2,1)*X(2,2)
tt73 = X(2,3)**3
tt74 = 2*X(1,2)*X(1,3)
tt75 = 2*X(2,2)*X(2,3)
tt76 = -X(2,4)
tt77 = tt76+X(2,3)
tt78 = X(3,2)**3
tt79 = tt15+tt74+tt25
tt80 = tt79*X(2,2)
tt81 = X(2,2)**3
tt82 = -tt81
tt83 = -2*X(2,3)
tt84 = 2*X(2,2)
tt85 = X(2,4)+tt50
tt86 = X(3,3)**3
tt87 = tt40*X(2,2)
tt88 = -tt73
tt89 = 2*X(2,3)
tt90 = tt89+tt60
tt91 = tt3*X(1,3)
tt92 = -2*X(1,2)*tt5
tt93 = X(1,3)**3
tt94 = tt79*X(1,4)
tt95 = (tt94+tt93+tt92+tt91)*X(2,2)
tt96 = -X(1,4)
tt97 = tt96+X(1,3)
tt98 = tt97*tt81
tt99 = X(1,2)**3
tt100 = -tt99
tt101 = 2*tt3*X(1,3)
tt102 = -X(1,2)*tt5
tt103 = tt40*X(1,4)
tt104 = -2*X(1,3)
tt105 = 3*X(1,4)
tt106 = tt105+tt104+tt16
tt107 = tt106*tt6
tt108 = (tt107+tt103+tt102+tt101+tt100)*X(2,3)
tt109 = 2*X(1,2)
tt110 = -3*X(1,4)
tt111 = tt110+X(1,3)+tt109
tt112 = tt111*X(2,2)*tt8
tt113 = X(1,4)+tt16
tt114 = tt113*tt73
tt115 = -tt93
tt116 = 2*X(1,3)
tt117 = tt116+tt30
tt118 = tt26*tt8+tt117*X(2,2)*X(2,3)+tt26*tt6+tt115+3*X(1,2)*tt5-&
&3*tt3*X(1,3)+tt99
tt119 = tt118*X(2,4)
tt120 = -tt3*X(1,3)
tt121 = 2*X(1,2)*tt5
tt122 = (tt103+tt115+tt121+tt120)*X(2,1)
tt123 = X(1,4)+tt21
tt124 = tt123*X(2,1)*tt6
tt125 = X(1,1)*tt3
tt126 = -2*X(1,1)*X(1,2)*X(1,3)
tt127 = X(1,1)*tt5
tt128 = -2*X(1,4)
tt129 = tt128+tt116
tt130 = tt129*X(2,1)*X(2,2)
tt131 = tt96+X(1,1)
tt132 = tt131*tt6
tt133 = (tt132+tt130+tt94+tt127+tt126+tt125)*X(2,3)
tt134 = tt123*X(2,1)
tt135 = 2*X(1,4)
tt136 = tt135+tt32
tt137 = tt136*X(2,2)+tt134
tt138 = tt137*tt8
tt139 = tt131*tt73
tt140 = -X(1,1)*tt3
tt141 = 2*X(1,1)*X(1,2)
tt142 = tt30+tt44
tt143 = X(1,3)+tt44
tt144 = 2*X(1,1)
tt145 = tt104+tt144
tt146 = tt143*tt8+tt145*X(2,2)*X(2,3)+tt143*tt6+tt93+tt142*tt5+(t&
&t3+tt141)*X(1,3)+tt140
tt147 = tt146*X(2,4)
tt148 = tt26*X(2,4)+tt113*X(2,3)+tt97*X(2,2)
tt149 = tt131*X(2,3)
tt150 = tt143*X(2,4)+tt149+tt134
tt151 = -2*tt3*X(1,3)
tt152 = X(1,2)*tt5
tt153 = (tt94+tt152+tt151+tt99)*X(2,1)
tt154 = 2*X(1,1)*X(1,2)*X(1,3)
tt155 = -X(1,1)*tt5
tt156 = (tt103+tt155+tt154+tt140)*X(2,2)
tt157 = tt96+X(1,2)
tt158 = tt157*X(2,1)*tt6
tt159 = X(1,4)+tt44
tt160 = tt159*tt81
tt161 = tt135+tt30
tt162 = tt161*X(2,1)*X(2,2)
tt163 = tt128+tt144
tt164 = tt163*tt6
tt165 = (tt164+tt162)*X(2,3)
tt166 = tt157*X(2,1)
tt167 = tt159*X(2,2)
tt168 = tt167+tt166
tt169 = tt168*tt8
tt170 = -2*X(1,1)*X(1,2)
tt171 = 2*tt3
tt172 = (tt171+tt170)*X(1,3)
tt173 = tt37*tt5
tt174 = tt37*tt6
tt175 = tt109+tt32
tt176 = tt175*X(2,2)*X(2,3)
tt177 = tt37*tt8
tt178 = (tt177+tt176+tt174+tt173+tt172+tt100+tt125)*X(2,4)
tt179 = tt135+tt104
tt180 = tt128+tt109
tt181 = tt117*X(2,4)+tt180*X(2,3)+tt179*X(2,2)
tt182 = tt181*X(3,1)*X(3,2)
tt183 = tt110+tt116+X(1,2)
tt184 = 3*X(1,1)
tt185 = tt104+tt16+tt184
tt186 = tt185*X(2,4)+tt136*X(2,3)+tt167+tt183*X(2,1)
tt187 = tt186*tt9
tt188 = tt105+tt21+tt30
tt189 = -3*X(1,1)
tt190 = X(1,3)+tt109+tt189
tt191 = tt190*X(2,4)+tt149+tt163*X(2,2)+tt188*X(2,1)
tt192 = tt191*X(3,2)+tt148*X(3,1)
tt193 = tt37*X(2,4)+tt167+tt166
tt194 = 3*tt3*X(1,3)
tt195 = -3*X(1,2)*tt5
tt196 = (tt25+tt170)*X(1,3)
tt197 = tt109+X(1,1)
tt198 = tt197*tt5
tt199 = -2*tt3
tt200 = (tt199+tt141)*X(1,3)
tt201 = tt45*tt5
tt202 = tt104+tt109
tt203 = tt202*X(2,1)*X(2,2)
tt204 = tt116+X(1,2)+tt189
tt205 = tt204*tt6
tt206 = tt17*X(2,1)
tt207 = tt21+tt30+tt184
tt208 = tt207*X(2,2)+tt206
tt209 = tt45*X(2,3)
tt210 = tt209+tt22*X(2,2)+tt206
tt211 = tt202*X(2,1)
tt212 = tt116+tt32
tt213 = tt30+tt144
tt214 = tt213*X(2,3)+tt212*X(2,2)+tt211
tt215 = tt210*tt10+tt214*X(3,2)*X(3,3)+tt210*tt9+tt45*tt73+tt208*&
&tt8+(tt205+tt203+tt201+tt200+tt99+tt140)*X(2,3)+tt22*tt81+tt17*X(&
&2,1)*tt6+(tt115+tt198+tt196+tt125)*X(2,2)+(tt93+tt195+tt194+tt100&
&)*X(2,1)
tt216 = tt215*X(3,4)+tt193*tt86+tt192*tt10+(tt187+tt182+tt178+tt1&
&69+tt165+tt160+tt158+tt156+tt153)*X(3,3)+tt150*tt78+tt148*X(3,1)*&
&tt9+(tt147+tt139+tt138+tt133+tt124+tt122)*X(3,2)+(tt119+tt114+tt1&
&12+tt108+tt98+tt95)*X(3,1)
tt217 = 1/(tt216**2/tt11+tt72**2)
tt218 = tt128+X(1,3)+X(1,2)
tt219 = atan2(tt13*tt216,tt72)-d(1,1)
tt220 = X(2,4)+tt53
tt221 = -2*X(2,4)
tt222 = X(1,4)+X(1,3)+tt30
tt223 = X(1,4)+tt104+X(1,1)
tt224 = tt96+tt109+tt44
tt225 = X(1,3)+tt30+X(1,1)
tt226 = -2*tt5
tt227 = tt117*X(1,4)
tt228 = -3*tt3
tt229 = 4*X(1,2)*X(1,3)
tt230 = tt202*X(1,4)
tt231 = 3*tt3
tt232 = -6*X(1,2)*X(1,3)
tt233 = 3*tt5
tt234 = 2*tt5
tt235 = tt227-2*X(1,1)*X(1,3)+tt141
tt236 = -4*X(1,2)*X(1,3)
tt237 = 2*X(1,1)*X(1,3)
tt238 = X(2,1)*tt6
tt239 = tt76+X(2,1)
tt240 = 2*X(2,4)
tt241 = 6*X(1,2)*X(1,3)
tt242 = -3*tt5
tt243 = -X(2,1)*tt6
tt244 = 2*X(2,1)*X(2,2)
tt245 = X(2,3)+tt69
tt246 = 2*X(2,1)
tt247 = 1/tt12**3
tt248 = tt83+tt84
tt249 = X(1,4)+tt104+X(1,2)
tt250 = tt96+tt116+tt44
tt251 = X(1,4)+tt30+X(1,1)
tt252 = tt104+X(1,2)+X(1,1)
tt253 = -2*tt6
tt254 = tt76+X(2,2)
tt255 = X(2,4)+tt69
tt256 = tt254*X(3,1)
tt257 = -2*X(2,1)*X(2,2)
tt258 = tt26*X(3,1)
jac(1,1) = 6*tt1*tt2*(tt13*tt72*((tt57*tt10+tt90*X(3,2)*X(3,3)+tt&
&57*tt9+tt88+3*X(2,2)*tt8+((-3*tt6)+tt15+tt74+tt25)*X(2,3)+tt81+tt&
&87)*X(3,4)+tt85*tt86+((-3*X(2,4))+X(2,3)+tt84)*X(3,2)*tt10+((3*X(&
&2,4)+tt83+tt50)*tt9+tt67*X(2,4)-X(2,2)*tt8+2*tt6*X(2,3)+tt82+tt80&
&)*X(3,3)+tt77*tt78+((tt49+tt75+tt56+tt15+tt74+tt25)*X(2,4)+tt73-2&
&*X(2,2)*tt8+(tt6+tt5+tt4+tt3)*X(2,3))*X(3,2))*tt217-tt13*((-(tt26&
&*X(3,3)+tt17*X(3,2))*X(3,4))-tt113*tt10-tt218*X(3,2)*X(3,3)-tt123&
&*tt9-(tt26*X(2,3)+tt17*X(2,2))*X(2,4)-tt113*tt8-tt218*X(2,2)*X(2,&
&3)-tt123*tt6)*tt216*tt217)*tt219
jac(1,2) = 6*tt1*tt2*(tt13*tt72*((tt17*tt10+tt202*X(3,2)*X(3,3)+t&
&t17*tt9+tt17*tt8+tt202*X(2,2)*X(2,3)+tt17*tt6+tt93+tt195+tt194+tt&
&100)*X(3,4)+tt157*tt86+tt188*X(3,2)*tt10+(tt183*tt9+tt157*tt8+tt1&
&61*X(2,2)*X(2,3)+tt157*tt6+tt94+tt152+tt151+tt99)*X(3,3)+tt123*tt&
&78+(tt123*tt8+tt129*X(2,2)*X(2,3)+tt123*tt6+tt103+tt115+tt121+tt1&
&20)*X(3,2))*tt217-tt13*((-(tt57*X(3,3)+tt51*X(3,2))*X(3,4))-tt85*&
&tt10-(tt221+X(2,3)+X(2,2))*X(3,2)*X(3,3)-tt220*tt9-tt40*X(2,4)-tt&
&28*X(2,3)-tt19*X(2,2))*tt216*tt217)*tt219
jac(1,3) = 6*tt1*tt2*(tt13*(tt148*tt10+tt181*X(3,2)*X(3,3)+tt148*&
&tt9+tt119+tt114+tt112+tt108+tt98+tt95)*tt72*tt217-tt13*((-tt67*X(&
&3,4))-tt58*X(3,3)-tt52*X(3,2))*tt216*tt217)*tt219
jac(1,4) = 6*tt1*tt2*(tt72*(tt13*((tt245*tt10+(tt83+tt246)*X(3,2)&
&*X(3,3)+tt245*tt9+tt73+(tt60+tt69)*tt8+(tt6+tt244+tt5+(tt144-4*X(&
&1,2))*X(1,3)+tt231+tt170)*X(2,3)+tt243+(tt234+(tt30+tt32)*X(1,3)+&
&tt141)*X(2,2)+(tt242+tt241+tt228)*X(2,1))*X(3,4)+tt239*tt86+((tt2&
&40+tt61)*X(3,2)+tt220*X(3,1))*tt10+(tt239*tt9+(tt221+tt89)*X(3,1)&
&*X(3,2)+(tt49+tt75+tt56+tt15+(4*X(1,2)+tt32)*X(1,3)+tt228+tt141)*&
&X(2,4)+X(2,1)*tt8-2*X(2,1)*X(2,2)*X(2,3)+tt238+(tt230+tt237+tt170&
&)*X(2,2)+(tt227+tt5+tt236+tt231)*X(2,1))*X(3,3)+tt220*X(3,1)*tt9+&
&((tt226+(tt109+tt144)*X(1,3)+tt170)*X(2,4)+tt235*X(2,3)+(tt230+tt&
&234+tt4)*X(2,1))*X(3,2)+((tt8+tt7+tt6+tt233+tt232+tt231)*X(2,4)+t&
&t88+2*X(2,2)*tt8+(tt56+tt230+tt15+tt229+tt228)*X(2,3)+(tt227+tt22&
&6+tt74)*X(2,2))*X(3,1))-(tt202*tt247*tt216)/2.0E+0)*tt217-tt13*((&
&-(tt225*X(3,3)+tt143*X(3,2)+tt202*X(3,1))*X(3,4))-tt224*tt10-(tt2&
&23*X(3,2)+tt222*X(3,1))*X(3,3)-tt97*X(3,1)*X(3,2)-(tt225*X(2,3)+t&
&t143*X(2,2)+tt211)*X(2,4)-tt224*tt8-(tt223*X(2,2)+tt222*X(2,1))*X&
&(2,3)-tt97*X(2,1)*X(2,2))*tt216*tt217)*tt219
jac(1,5) = 6*tt1*tt2*(tt72*(tt13*((tt22*tt10+tt212*X(3,2)*X(3,3)+&
&tt22*tt9+tt207*tt8+(2*tt204*X(2,2)+tt211)*X(2,3)+3*tt22*tt6+2*tt1&
&7*X(2,1)*X(2,2)+tt115+tt198+tt196+tt125)*X(3,4)+tt159*tt86+(tt163&
&*X(3,2)+tt97*X(3,1))*tt10+(tt159*tt9+tt179*X(3,1)*X(3,2)+(tt175*X&
&(2,3)+2*tt37*X(2,2))*X(2,4)+tt159*tt8+(2*tt163*X(2,2)+tt161*X(2,1&
&))*X(2,3)+3*tt159*tt6+2*tt157*X(2,1)*X(2,2)+tt103+tt155+tt154+tt1&
&40)*X(3,3)+tt97*X(3,1)*tt9+((tt145*X(2,3)+2*tt143*X(2,2))*X(2,4)+&
&tt136*tt8+(2*tt131*X(2,2)+tt129*X(2,1))*X(2,3)+2*tt123*X(2,1)*X(2&
&,2))*X(3,2)+((tt117*X(2,3)+2*tt26*X(2,2))*X(2,4)+tt111*tt8+2*tt10&
&6*X(2,2)*X(2,3)+3*tt97*tt6+tt94+tt93+tt92+tt91)*X(3,1))-(tt248*tt&
&247*tt216)/2.0E+0)*tt217-tt13*((-((X(2,3)+tt60+X(2,1))*X(3,3)+tt2&
&45*X(3,2)+tt248*X(3,1))*X(3,4))-(tt76+tt84+tt69)*tt10-((X(2,4)+tt&
&83+X(2,1))*X(3,2)+(X(2,4)+X(2,3)+tt60)*X(3,1))*X(3,3)-tt77*X(3,1)&
&*X(3,2)-tt43*X(2,4)-tt35*X(2,3)-2*tt24*X(2,2)-tt19*X(2,1))*tt216*&
&tt217)*tt219
jac(1,6) = 6*tt1*tt2*(tt72*(tt13*((tt214*X(3,3)+2*tt210*X(3,2))*X&
&(3,4)+tt191*tt10+(2*tt186*X(3,2)+tt181*X(3,1))*X(3,3)+3*tt150*tt9&
&+2*tt148*X(3,1)*X(3,2)+tt147+tt139+tt138+tt133+tt124+tt122)-((2*X&
&(3,2)-2*X(3,3))*tt247*tt216)/2.0E+0)*tt217-tt13*((-tt68*X(3,4))-t&
&t63*X(3,3)-2*tt55*X(3,2)-tt52*X(3,1))*tt216*tt217)*tt219
jac(1,7) = 6*tt1*tt2*(tt72*(tt13*((tt65*tt10+(tt84+tt61)*X(3,2)*X&
&(3,3)+tt65*tt9+tt65*tt8+(2*tt6+tt257+2*tt45*X(1,3)+tt199+tt141)*X&
&(2,3)+tt82+tt238+(tt242+2*tt197*X(1,3)+tt25+tt170)*X(2,2)+(tt233+&
&tt232+tt231)*X(2,1))*X(3,4)+(tt255*X(3,2)+tt256)*tt10+((tt221+tt2&
&46)*tt9+(tt240+tt60)*X(3,1)*X(3,2)+(2*tt37*X(1,3)+tt171+tt170)*X(&
&2,4)+tt235*X(2,2)+(tt230+tt74+tt199)*X(2,1))*X(3,3)+tt255*tt78+tt&
&254*X(3,1)*tt9+((tt8+tt7+tt6+tt233+2*tt142*X(1,3)+tt3+tt141)*X(2,&
&4)-X(2,1)*tt8+(tt244+tt230+tt237+tt170)*X(2,3)+tt243+(tt227+tt242&
&+tt229+tt25)*X(2,1))*X(3,2)+((tt49+tt75+tt56+tt242+tt241+tt228)*X&
&(2,4)+X(2,2)*tt8+(tt253+tt227+tt4+tt171)*X(2,3)+tt81+(tt230+tt233&
&+tt236+tt3)*X(2,2))*X(3,1))-(tt117*tt247*tt216)/2.0E+0)*tt217-tt1&
&3*((-(tt45*X(3,3)+tt252*X(3,2)+tt117*X(3,1))*X(3,4))-(tt251*X(3,2&
&)+tt157*X(3,1))*X(3,3)-tt250*tt9-tt249*X(3,1)*X(3,2)-(tt209+tt252&
&*X(2,2)+tt117*X(2,1))*X(2,4)-(tt251*X(2,2)+tt166)*X(2,3)-tt250*tt&
&6-tt249*X(2,1)*X(2,2))*tt216*tt217)*tt219
jac(1,8) = 6*tt1*tt2*(tt72*(tt13*((tt45*tt10+tt213*X(3,2)*X(3,3)+&
&tt45*tt9+3*tt45*tt8+2*tt208*X(2,3)+tt205+tt203+tt201+tt200+tt99+t&
&t140)*X(3,4)+(tt131*X(3,2)+tt113*X(3,1))*tt10+(tt136*tt9+tt180*X(&
&3,1)*X(3,2)+(2*tt37*X(2,3)+tt175*X(2,2))*X(2,4)+2*tt168*X(2,3)+tt&
&164+tt162)*X(3,3)+tt131*tt78+tt113*X(3,1)*tt9+((2*tt143*X(2,3)+tt&
&145*X(2,2))*X(2,4)+3*tt131*tt8+2*tt137*X(2,3)+tt132+tt130+tt94+tt&
&127+tt126+tt125)*X(3,2)+((2*tt26*X(2,3)+tt117*X(2,2))*X(2,4)+3*tt&
&113*tt8+2*tt111*X(2,2)*X(2,3)+tt107+tt103+tt102+tt101+tt100)*X(3,&
&1))-(tt90*tt247*tt216)/2.0E+0)*tt217-tt13*((-(tt70*X(3,3)+(tt83+X&
&(2,2)+X(2,1))*X(3,2)+tt90*X(3,1))*X(3,4))-((X(2,4)+tt60+X(2,1))*X&
&(3,2)+tt256)*X(3,3)-(tt76+tt89+tt69)*tt9-(X(2,4)+tt83+X(2,2))*X(3&
&,1)*X(3,2)-tt47*X(2,4)-2*tt39*X(2,3)-tt35*X(2,2)-tt28*X(2,1))*tt2&
&16*tt217)*tt219
jac(1,9) = 6*tt1*tt2*(tt72*(tt13*((2*tt210*X(3,3)+tt214*X(3,2))*X&
&(3,4)+3*tt193*tt10+2*tt192*X(3,3)+tt187+tt182+tt178+tt169+tt165+t&
&t160+tt158+tt156+tt153)-((2*X(3,3)-2*X(3,2))*tt247*tt216)/2.0E+0)&
&*tt217-tt13*((-tt71*X(3,4))-2*tt66*X(3,3)-tt63*X(3,2)-tt58*X(3,1)&
&)*tt216*tt217)*tt219
jac(1,10) = 6*tt1*tt2*(tt13*(tt70*tt86+((tt53+tt60+3*X(2,1))*X(3,&
&2)+tt51*X(3,1))*tt10+((tt89+X(2,2)-3*X(2,1))*tt9+tt248*X(3,1)*X(3&
&,2)+tt70*tt8+(tt253+tt244)*X(2,3)+tt81+tt243+tt87+tt79*X(2,1))*X(&
&3,3)+tt54*tt78+tt51*X(3,1)*tt9+(tt88+(tt84+X(2,1))*tt8+(tt56+tt25&
&7+tt15+tt74+tt25)*X(2,3)+tt238+tt41)*X(3,2)+(tt73-3*X(2,2)*tt8+(3&
&*tt6+tt5+tt4+tt3)*X(2,3)+tt82+tt80)*X(3,1))*tt72*tt217-tt13*((-tt&
&37*tt10)-(tt33*X(3,2)+tt258)*X(3,3)-tt22*tt9-tt17*X(3,1)*X(3,2)-t&
&t37*tt8-(tt33*X(2,2)+tt26*X(2,1))*X(2,3)-tt22*tt6-tt17*X(2,1)*X(2&
&,2))*tt216*tt217)*tt219
jac(1,11) = 6*tt1*tt2*(tt13*(tt37*tt86+(tt190*X(3,2)+tt258)*tt10+&
&(tt185*tt9+tt117*X(3,1)*X(3,2)+tt177+tt176+tt174+tt173+tt172+tt10&
&0+tt125)*X(3,3)+tt143*tt78+tt26*X(3,1)*tt9+tt146*X(3,2)+tt118*X(3&
&,1))*tt72*tt217-tt13*((-tt65*tt10)-(tt62*X(3,2)+tt57*X(3,1))*X(3,&
&3)-tt54*tt9-tt51*X(3,1)*X(3,2)-tt47*X(2,3)-tt43*X(2,2)-tt40*X(2,1&
&))*tt216*tt217)*tt219
jac(1,12) = 6*tt1*tt2*(tt13*tt215*tt72*tt217-((-tt71*X(3,3))-tt68&
&*X(3,2)-tt67*X(3,1))*tt13*tt216*tt217)*tt219
END 
SUBROUTINE surf_bending_hes(hes, X, d, l, area) 
IMPLICIT NONE 
REAL(KIND=8) hes(12, 12) 
REAL(KIND=8) X(3, 4) 
REAL(KIND=8) d(1, 1) 
REAL(KIND=8) l(1, 1) 
REAL(KIND=8) area(1, 1) 
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
tt1 = 1/area(1,1)
tt2 = l(1,1)**2
tt3 = X(1,2)**2
tt4 = -2*X(1,2)*X(1,3)
tt5 = X(1,3)**2
tt6 = X(2,2)**2
tt7 = -2*X(2,2)*X(2,3)
tt8 = X(2,3)**2
tt9 = X(3,2)**2
tt10 = -2*X(3,2)*X(3,3)
tt11 = X(3,3)**2
tt12 = tt11+tt10+tt9+tt8+tt7+tt6+tt5+tt4+tt3
tt13 = sqrt(tt12)
tt14 = 1/tt13
tt15 = X(1,2)*X(1,3)
tt16 = -tt5
tt17 = -X(1,2)
tt18 = X(1,3)+tt17
tt19 = tt18*X(1,4)
tt20 = tt19+tt16+tt15
tt21 = -X(1,1)*X(1,3)
tt22 = -X(1,3)
tt23 = tt22+X(1,1)
tt24 = tt23*X(1,4)
tt25 = tt24+tt5+tt21
tt26 = -tt3
tt27 = tt22+X(1,2)
tt28 = tt27*X(1,4)
tt29 = tt28+tt15+tt26
tt30 = X(1,1)*X(1,2)
tt31 = -2*X(1,2)
tt32 = tt31+X(1,1)
tt33 = tt32*X(1,3)
tt34 = -2*X(1,1)
tt35 = X(1,3)+X(1,2)+tt34
tt36 = tt35*X(1,4)
tt37 = tt36+tt33+tt30
tt38 = -X(1,1)*X(1,2)
tt39 = tt17+X(1,1)
tt40 = tt39*X(1,4)
tt41 = tt40+tt3+tt38
tt42 = tt5+tt4+tt3
tt43 = tt42*X(2,1)
tt44 = X(1,2)+X(1,1)
tt45 = tt44*X(1,3)
tt46 = tt16+tt45+tt38
tt47 = -X(1,1)
tt48 = X(1,2)+tt47
tt49 = tt48*X(1,3)
tt50 = tt49+tt26+tt30
tt51 = X(2,2)*X(2,3)
tt52 = -tt8
tt53 = -X(2,2)
tt54 = X(2,3)+tt53
tt55 = tt54*X(2,4)+tt52+tt51+tt19+tt16+tt15
tt56 = -X(2,3)
tt57 = tt56+X(2,1)
tt58 = tt57*X(2,4)+tt8-X(2,1)*X(2,3)+tt24+tt5+tt21
tt59 = -tt6
tt60 = tt56+X(2,2)
tt61 = tt60*X(2,4)+tt51+tt59+tt28+tt15+tt26
tt62 = X(2,1)*X(2,2)
tt63 = -2*X(2,2)
tt64 = tt63+X(2,1)
tt65 = -2*X(2,1)
tt66 = X(2,3)+X(2,2)+tt65
tt67 = tt66*X(2,4)+tt64*X(2,3)+tt62+tt36+tt33+tt30
tt68 = -X(2,1)*X(2,2)
tt69 = tt53+X(2,1)
tt70 = tt69*X(2,4)+tt6+tt68+tt40+tt3+tt38
tt71 = tt8+tt7+tt6+tt5+tt4+tt3
tt72 = X(2,2)+X(2,1)
tt73 = tt52+tt72*X(2,3)+tt68+tt16+tt45+tt38
tt74 = -X(2,1)
tt75 = X(2,2)+tt74
tt76 = tt75*X(2,3)+tt59+tt62+tt49+tt26+tt30
tt77 = (-(tt76*X(3,3)+tt73*X(3,2)+tt71*X(3,1))*X(3,4))-tt70*tt11-&
&(tt67*X(3,2)+tt61*X(3,1))*X(3,3)-tt58*tt9-tt55*X(3,1)*X(3,2)-(tt5&
&0*X(2,3)+tt46*X(2,2)+tt43)*X(2,4)-tt41*tt8-(tt37*X(2,2)+tt29*X(2,&
&1))*X(2,3)-tt25*tt6-tt20*X(2,1)*X(2,2)
tt78 = (tt6+tt5+tt4+tt3)*X(2,3)
tt79 = -2*X(2,2)*tt8
tt80 = X(2,3)**3
tt81 = 2*X(1,2)*X(1,3)
tt82 = 2*X(2,2)*X(2,3)
tt83 = tt52+tt82+tt59+tt16+tt81+tt26
tt84 = tt83*X(2,4)
tt85 = -X(2,4)
tt86 = tt85+X(2,3)
tt87 = X(3,2)**3
tt88 = tt16+tt81+tt26
tt89 = tt88*X(2,2)
tt90 = X(2,2)**3
tt91 = -tt90
tt92 = 2*tt6*X(2,3)
tt93 = -X(2,2)*tt8
tt94 = tt71*X(2,4)
tt95 = -2*X(2,3)
tt96 = 3*X(2,4)+tt95+tt53
tt97 = tt96*tt9
tt98 = 2*X(2,2)
tt99 = (-3*X(2,4))+X(2,3)+tt98
tt100 = X(2,4)+tt53
tt101 = X(3,3)**3
tt102 = tt42*X(2,2)
tt103 = -3*tt6
tt104 = -tt80
tt105 = 2*X(2,3)
tt106 = tt105+tt63
tt107 = tt60*tt11+tt106*X(3,2)*X(3,3)+tt60*tt9+tt104+3*X(2,2)*tt8&
&+(tt103+tt16+tt81+tt26)*X(2,3)+tt90+tt102
tt108 = tt107*X(3,4)+tt100*tt101+tt99*X(3,2)*tt11+(tt97+tt94+tt93&
&+tt92+tt91+tt89)*X(3,3)+tt86*tt87+(tt84+tt80+tt79+tt78)*X(3,2)
tt109 = 1/tt12
tt110 = tt3*X(1,3)
tt111 = -2*X(1,2)*tt5
tt112 = X(1,3)**3
tt113 = tt88*X(1,4)
tt114 = (tt113+tt112+tt111+tt110)*X(2,2)
tt115 = -X(1,4)
tt116 = tt115+X(1,3)
tt117 = tt116*tt90
tt118 = X(1,2)**3
tt119 = -tt118
tt120 = 2*tt3*X(1,3)
tt121 = -X(1,2)*tt5
tt122 = tt42*X(1,4)
tt123 = -2*X(1,3)
tt124 = 3*X(1,4)
tt125 = tt124+tt123+tt17
tt126 = tt125*tt6
tt127 = (tt126+tt122+tt121+tt120+tt119)*X(2,3)
tt128 = 2*X(1,2)
tt129 = -3*X(1,4)
tt130 = tt129+X(1,3)+tt128
tt131 = tt130*X(2,2)*tt8
tt132 = X(1,4)+tt17
tt133 = tt132*tt80
tt134 = -3*tt3*X(1,3)
tt135 = 3*X(1,2)*tt5
tt136 = -tt112
tt137 = tt27*tt6
tt138 = 2*X(1,3)
tt139 = tt138+tt31
tt140 = tt139*X(2,2)*X(2,3)
tt141 = tt27*tt8
tt142 = tt141+tt140+tt137+tt136+tt135+tt134+tt118
tt143 = tt142*X(2,4)
tt144 = -tt3*X(1,3)
tt145 = 2*X(1,2)*tt5
tt146 = (tt122+tt136+tt145+tt144)*X(2,1)
tt147 = X(1,4)+tt22
tt148 = tt147*X(2,1)*tt6
tt149 = X(1,1)*tt3
tt150 = -2*X(1,1)*X(1,2)*X(1,3)
tt151 = X(1,1)*tt5
tt152 = -2*X(1,4)
tt153 = tt152+tt138
tt154 = tt153*X(2,1)*X(2,2)
tt155 = tt115+X(1,1)
tt156 = tt155*tt6
tt157 = (tt156+tt154+tt113+tt151+tt150+tt149)*X(2,3)
tt158 = tt147*X(2,1)
tt159 = 2*X(1,4)
tt160 = tt159+tt34
tt161 = tt160*X(2,2)
tt162 = tt161+tt158
tt163 = tt162*tt8
tt164 = tt155*tt80
tt165 = -X(1,1)*tt3
tt166 = 2*X(1,1)*X(1,2)
tt167 = (tt3+tt166)*X(1,3)
tt168 = tt31+tt47
tt169 = tt168*tt5
tt170 = X(1,3)+tt47
tt171 = tt170*tt6
tt172 = 2*X(1,1)
tt173 = tt123+tt172
tt174 = tt173*X(2,2)*X(2,3)
tt175 = tt170*tt8
tt176 = tt175+tt174+tt171+tt112+tt169+tt167+tt165
tt177 = tt176*X(2,4)
tt178 = tt27*X(2,4)+tt132*X(2,3)+tt116*X(2,2)
tt179 = tt155*X(2,3)
tt180 = tt170*X(2,4)+tt179+tt158
tt181 = -2*tt3*X(1,3)
tt182 = X(1,2)*tt5
tt183 = (tt113+tt182+tt181+tt118)*X(2,1)
tt184 = 2*X(1,1)*X(1,2)*X(1,3)
tt185 = -X(1,1)*tt5
tt186 = (tt122+tt185+tt184+tt165)*X(2,2)
tt187 = tt115+X(1,2)
tt188 = tt187*X(2,1)*tt6
tt189 = X(1,4)+tt47
tt190 = tt189*tt90
tt191 = tt159+tt31
tt192 = tt191*X(2,1)*X(2,2)
tt193 = tt152+tt172
tt194 = tt193*tt6
tt195 = (tt194+tt192)*X(2,3)
tt196 = tt187*X(2,1)
tt197 = tt189*X(2,2)
tt198 = tt197+tt196
tt199 = tt198*tt8
tt200 = -2*X(1,1)*X(1,2)
tt201 = 2*tt3
tt202 = (tt201+tt200)*X(1,3)
tt203 = tt39*tt5
tt204 = tt39*tt6
tt205 = tt128+tt34
tt206 = tt205*X(2,2)*X(2,3)
tt207 = tt39*tt8
tt208 = (tt207+tt206+tt204+tt203+tt202+tt119+tt149)*X(2,4)
tt209 = tt159+tt123
tt210 = tt152+tt128
tt211 = tt139*X(2,4)
tt212 = tt211+tt210*X(2,3)+tt209*X(2,2)
tt213 = tt212*X(3,1)*X(3,2)
tt214 = tt129+tt138+X(1,2)
tt215 = tt160*X(2,3)
tt216 = 3*X(1,1)
tt217 = tt123+tt17+tt216
tt218 = tt217*X(2,4)+tt215+tt197+tt214*X(2,1)
tt219 = tt218*tt9
tt220 = tt124+tt22+tt31
tt221 = tt193*X(2,2)
tt222 = -3*X(1,1)
tt223 = X(1,3)+tt128+tt222
tt224 = tt223*X(2,4)+tt179+tt221+tt220*X(2,1)
tt225 = tt224*X(3,2)+tt178*X(3,1)
tt226 = tt39*X(2,4)+tt197+tt196
tt227 = 3*tt3*X(1,3)
tt228 = -3*X(1,2)*tt5
tt229 = (tt26+tt200)*X(1,3)
tt230 = tt128+X(1,1)
tt231 = tt230*tt5
tt232 = -2*tt3
tt233 = (tt232+tt166)*X(1,3)
tt234 = tt48*tt5
tt235 = tt123+tt128
tt236 = tt235*X(2,1)*X(2,2)
tt237 = tt138+X(1,2)+tt222
tt238 = tt237*tt6
tt239 = tt18*X(2,1)
tt240 = tt22+tt31+tt216
tt241 = tt240*X(2,2)+tt239
tt242 = tt48*X(2,3)
tt243 = tt242+tt23*X(2,2)+tt239
tt244 = tt235*X(2,1)
tt245 = tt138+tt34
tt246 = tt31+tt172
tt247 = tt246*X(2,3)+tt245*X(2,2)+tt244
tt248 = tt243*tt11+tt247*X(3,2)*X(3,3)+tt243*tt9+tt48*tt80+tt241*&
&tt8+(tt238+tt236+tt234+tt233+tt118+tt165)*X(2,3)+tt23*tt90+tt18*X&
&(2,1)*tt6+(tt136+tt231+tt229+tt149)*X(2,2)+(tt112+tt228+tt227+tt1&
&19)*X(2,1)
tt249 = tt248*X(3,4)+tt226*tt101+tt225*tt11+(tt219+tt213+tt208+tt&
&199+tt195+tt190+tt188+tt186+tt183)*X(3,3)+tt180*tt87+tt178*X(3,1)&
&*tt9+(tt177+tt164+tt163+tt157+tt148+tt146)*X(3,2)+(tt143+tt133+tt&
&131+tt127+tt117+tt114)*X(3,1)
tt250 = tt249**2
tt251 = tt109*tt250+tt77**2
tt252 = 1/tt251
tt253 = tt152+X(1,3)+X(1,2)
tt254 = (-(tt27*X(3,3)+tt18*X(3,2))*X(3,4))-tt132*tt11-tt253*X(3,&
&2)*X(3,3)-tt147*tt9-(tt27*X(2,3)+tt18*X(2,2))*X(2,4)-tt132*tt8-tt&
&253*X(2,2)*X(2,3)-tt147*tt6
tt255 = tt14*tt77*tt108*tt252-tt14*tt254*tt249*tt252
tt256 = 2*tt109*tt108*tt249+2*tt254*tt77
tt257 = 1/tt251**2
tt258 = atan2(tt14*tt249,tt77)-d(1,1)
tt259 = tt147*tt6
tt260 = tt153*X(2,2)*X(2,3)
tt261 = tt147*tt8
tt262 = tt187*tt6
tt263 = tt191*X(2,2)*X(2,3)
tt264 = tt187*tt8
tt265 = tt214*tt9
tt266 = tt18*tt11+tt235*X(3,2)*X(3,3)+tt18*tt9+tt18*tt8+tt235*X(2&
&,2)*X(2,3)+tt18*tt6+tt112+tt228+tt227+tt119
tt267 = tt266*X(3,4)+tt187*tt101+tt220*X(3,2)*tt11+(tt265+tt264+t&
&t263+tt262+tt113+tt182+tt181+tt118)*X(3,3)+tt147*tt87+(tt261+tt26&
&0+tt259+tt122+tt136+tt145+tt144)*X(3,2)
tt268 = X(2,4)+tt56
tt269 = -2*X(2,4)
tt270 = tt269+X(2,3)+X(2,2)
tt271 = (-(tt60*X(3,3)+tt54*X(3,2))*X(3,4))-tt100*tt11-tt270*X(3,&
&2)*X(3,3)-tt268*tt9-tt42*X(2,4)-tt29*X(2,3)-tt20*X(2,2)
tt272 = tt14*tt77*tt267*tt252-tt14*tt271*tt249*tt252
tt273 = 6*tt1*tt2*tt255*tt272
tt274 = 2*tt109*tt267*tt249+2*tt271*tt77
tt275 = tt178*tt11+tt212*X(3,2)*X(3,3)+tt178*tt9+tt143+tt133+tt13&
&1+tt127+tt117+tt114
tt276 = (-tt71*X(3,4))-tt61*X(3,3)-tt55*X(3,2)
tt277 = tt14*tt275*tt77*tt252-tt14*tt276*tt249*tt252
tt278 = 6*tt1*tt2*tt277*tt255
tt279 = 2*tt109*tt275*tt249+2*tt276*tt77
tt280 = X(1,4)+X(1,3)+tt31
tt281 = X(1,4)+tt123+X(1,1)
tt282 = tt115+tt128+tt47
tt283 = X(1,3)+tt31+X(1,1)
tt284 = (-(tt283*X(3,3)+tt170*X(3,2)+tt235*X(3,1))*X(3,4))-tt282*&
&tt11-(tt281*X(3,2)+tt280*X(3,1))*X(3,3)-tt116*X(3,1)*X(3,2)-(tt28&
&3*X(2,3)+tt170*X(2,2)+tt244)*X(2,4)-tt282*tt8-(tt281*X(2,2)+tt280&
&*X(2,1))*X(2,3)-tt116*X(2,1)*X(2,2)
tt285 = -2*tt5
tt286 = tt139*X(1,4)
tt287 = (tt286+tt285+tt81)*X(2,2)
tt288 = -3*tt3
tt289 = 4*X(1,2)*X(1,3)
tt290 = tt235*X(1,4)
tt291 = (tt59+tt290+tt16+tt289+tt288)*X(2,3)
tt292 = 2*X(2,2)*tt8
tt293 = 3*tt3
tt294 = -6*X(1,2)*X(1,3)
tt295 = 3*tt5
tt296 = tt8+tt7+tt6+tt295+tt294+tt293
tt297 = tt296*X(2,4)
tt298 = 2*tt5
tt299 = tt290+tt298+tt4
tt300 = tt299*X(2,1)
tt301 = -2*X(1,1)*X(1,3)
tt302 = tt286+tt301+tt166
tt303 = tt302*X(2,3)
tt304 = tt285+(tt128+tt172)*X(1,3)+tt200
tt305 = tt304*X(2,4)
tt306 = -4*X(1,2)*X(1,3)
tt307 = (tt286+tt5+tt306+tt293)*X(2,1)
tt308 = 2*X(1,1)*X(1,3)
tt309 = (tt290+tt308+tt200)*X(2,2)
tt310 = X(2,1)*tt6
tt311 = -2*X(2,1)*X(2,2)*X(2,3)
tt312 = X(2,1)*tt8
tt313 = 4*X(1,2)
tt314 = (tt313+tt34)*X(1,3)
tt315 = (tt52+tt82+tt59+tt16+tt314+tt288+tt166)*X(2,4)
tt316 = tt269+tt105
tt317 = tt316*X(3,1)*X(3,2)
tt318 = tt85+X(2,1)
tt319 = tt318*tt9
tt320 = 2*X(2,4)
tt321 = tt320+tt65
tt322 = tt321*X(3,2)+tt268*X(3,1)
tt323 = 6*X(1,2)*X(1,3)
tt324 = -3*tt5
tt325 = (tt31+tt34)*X(1,3)
tt326 = -X(2,1)*tt6
tt327 = -4*X(1,2)
tt328 = (tt327+tt172)*X(1,3)
tt329 = 2*X(2,1)*X(2,2)
tt330 = tt63+tt74
tt331 = X(2,3)+tt74
tt332 = 2*X(2,1)
tt333 = tt95+tt332
tt334 = tt331*tt11+tt333*X(3,2)*X(3,3)+tt331*tt9+tt80+tt330*tt8+(&
&tt6+tt329+tt5+tt328+tt293+tt200)*X(2,3)+tt326+(tt298+tt325+tt166)&
&*X(2,2)+(tt324+tt323+tt288)*X(2,1)
tt335 = tt334*X(3,4)+tt318*tt101+tt322*tt11+(tt319+tt317+tt315+tt&
&312+tt311+tt310+tt309+tt307)*X(3,3)+tt268*X(3,1)*tt9+(tt305+tt303&
&+tt300)*X(3,2)+(tt297+tt104+tt292+tt291+tt287)*X(3,1)
tt336 = 1/tt13**3
tt337 = tt14*tt335-(tt235*tt336*tt249)/2.0E+0
tt338 = tt77*tt337*tt252-tt14*tt284*tt249*tt252
tt339 = 6*tt1*tt2*tt255*tt338
tt340 = 1/tt12**2
tt341 = (-tt235*tt340*tt250)+2*tt109*tt335*tt249+2*tt284*tt77
tt342 = tt235*X(2,3)
tt343 = tt139*X(2,2)
tt344 = tt235*X(2,4)
tt345 = tt235*X(2,2)
tt346 = tt139*X(2,3)
tt347 = tt346+tt345
tt348 = tt347*X(3,4)+(tt344+tt343)*X(3,3)+(tt211+tt342)*X(3,2)
tt349 = -X(2,2)*X(2,3)
tt350 = -tt54*X(2,4)
tt351 = -X(3,2)*X(3,3)
tt352 = -X(3,2)
tt353 = -(X(3,3)+tt352)*X(3,4)
tt354 = -tt14*(tt353+tt11+tt351+tt350+tt8+tt349)*tt249*tt252
tt355 = X(2,4)+X(2,3)+tt63
tt356 = X(2,4)+tt95+X(2,1)
tt357 = tt85+tt98+tt74
tt358 = tt95+tt98
tt359 = tt358*X(3,1)
tt360 = X(2,3)+tt63+X(2,1)
tt361 = (-(tt360*X(3,3)+tt331*X(3,2)+tt359)*X(3,4))-tt357*tt11-(t&
&t356*X(3,2)+tt355*X(3,1))*X(3,3)-tt86*X(3,1)*X(3,2)-tt46*X(2,4)-t&
&t37*X(2,3)-2*tt25*X(2,2)-tt20*X(2,1)
tt362 = 3*tt116*tt6
tt363 = 2*tt125*X(2,2)*X(2,3)
tt364 = tt130*tt8
tt365 = tt346+2*tt27*X(2,2)
tt366 = tt365*X(2,4)
tt367 = 2*tt147*X(2,1)*X(2,2)
tt368 = tt153*X(2,1)
tt369 = 2*tt155*X(2,2)
tt370 = (tt369+tt368)*X(2,3)
tt371 = tt160*tt8
tt372 = tt173*X(2,3)+2*tt170*X(2,2)
tt373 = tt372*X(2,4)
tt374 = 2*tt187*X(2,1)*X(2,2)
tt375 = 3*tt189*tt6
tt376 = tt191*X(2,1)
tt377 = 2*tt193*X(2,2)
tt378 = (tt377+tt376)*X(2,3)
tt379 = tt189*tt8
tt380 = tt205*X(2,3)+2*tt39*X(2,2)
tt381 = tt380*X(2,4)
tt382 = tt209*X(3,1)*X(3,2)
tt383 = tt189*tt9
tt384 = tt193*X(3,2)+tt116*X(3,1)
tt385 = 2*tt237*X(2,2)
tt386 = tt23*tt11+tt245*X(3,2)*X(3,3)+tt23*tt9+tt240*tt8+(tt385+t&
&t244)*X(2,3)+3*tt23*tt6+2*tt18*X(2,1)*X(2,2)+tt136+tt231+tt229+tt&
&149
tt387 = tt386*X(3,4)+tt189*tt101+tt384*tt11+(tt383+tt382+tt381+tt&
&379+tt378+tt375+tt374+tt122+tt185+tt184+tt165)*X(3,3)+tt116*X(3,1&
&)*tt9+(tt373+tt371+tt370+tt367)*X(3,2)+(tt366+tt364+tt363+tt362+t&
&t113+tt112+tt111+tt110)*X(3,1)
tt388 = tt14*tt387-(tt358*tt336*tt249)/2.0E+0
tt389 = tt77*tt388*tt252-tt14*tt361*tt249*tt252
tt390 = 6*tt1*tt2*tt255*tt389
tt391 = (-tt358*tt340*tt250)+2*tt109*tt387*tt249+2*tt361*tt77
tt392 = -2*tt8
tt393 = tt106*X(2,4)
tt394 = 4*X(2,2)*X(2,3)
tt395 = tt358*X(2,4)
tt396 = -tt9
tt397 = -tt101
tt398 = 3*tt6
tt399 = -6*X(2,2)*X(2,3)
tt400 = 3*tt8
tt401 = (tt11+tt10+tt9+tt400+tt399+tt398+tt5+tt4+tt3)*X(3,4)+tt39&
&7+2*X(3,2)*tt11+(tt396+tt395+tt52+tt394+tt103+tt16+tt81+tt26)*X(3&
&,3)+(tt393+tt392+tt82)*X(3,2)
tt402 = -((-tt18*X(2,4))-tt253*X(2,3)-2*tt147*X(2,2))*tt14*tt249*&
&tt252
tt403 = (-tt73*X(3,4))-tt67*X(3,3)-2*tt58*X(3,2)-tt55*X(3,1)
tt404 = tt212*X(3,1)
tt405 = 2*tt218*X(3,2)
tt406 = tt247*X(3,3)+2*tt243*X(3,2)
tt407 = tt406*X(3,4)+tt224*tt11+(tt405+tt404)*X(3,3)+3*tt180*tt9+&
&2*tt178*X(3,1)*X(3,2)+tt177+tt164+tt163+tt157+tt148+tt146
tt408 = 2*X(3,2)
tt409 = -2*X(3,3)
tt410 = tt409+tt408
tt411 = tt14*tt407-(tt410*tt336*tt249)/2.0E+0
tt412 = tt77*tt411*tt252-tt14*tt403*tt249*tt252
tt413 = 6*tt1*tt2*tt255*tt412
tt414 = (-tt410*tt340*tt250)+2*tt109*tt407*tt249+2*tt403*tt77
tt415 = (tt106*X(3,3)+2*tt60*X(3,2))*X(3,4)+tt99*tt11+2*tt96*X(3,&
&2)*X(3,3)+3*tt86*tt9+tt84+tt80+tt79+tt78
tt416 = -tt14*((-tt18*X(3,4))-tt253*X(3,3)-2*tt147*X(3,2))*tt249*&
&tt252
tt417 = X(1,4)+tt123+X(1,2)
tt418 = tt115+tt138+tt47
tt419 = X(1,4)+tt31+X(1,1)
tt420 = tt139*X(2,1)
tt421 = tt123+X(1,2)+X(1,1)
tt422 = tt139*X(3,1)
tt423 = (-(tt48*X(3,3)+tt421*X(3,2)+tt422)*X(3,4))-(tt419*X(3,2)+&
&tt187*X(3,1))*X(3,3)-tt418*tt9-tt417*X(3,1)*X(3,2)-(tt242+tt421*X&
&(2,2)+tt420)*X(2,4)-(tt419*X(2,2)+tt196)*X(2,3)-tt418*tt6-tt417*X&
&(2,1)*X(2,2)
tt424 = (tt290+tt295+tt306+tt3)*X(2,2)
tt425 = -2*tt6
tt426 = (tt425+tt286+tt4+tt201)*X(2,3)
tt427 = X(2,2)*tt8
tt428 = tt52+tt82+tt59+tt324+tt323+tt288
tt429 = tt428*X(2,4)
tt430 = (tt286+tt324+tt289+tt26)*X(2,1)
tt431 = (tt329+tt290+tt308+tt200)*X(2,3)
tt432 = -X(2,1)*tt8
tt433 = tt8+tt7+tt6+tt295+2*tt168*X(1,3)+tt3+tt166
tt434 = tt433*X(2,4)
tt435 = tt85+X(2,2)
tt436 = X(2,4)+tt74
tt437 = (tt290+tt81+tt232)*X(2,1)
tt438 = tt302*X(2,2)
tt439 = 2*tt39*X(1,3)
tt440 = (tt439+tt201+tt200)*X(2,4)
tt441 = tt320+tt63
tt442 = tt441*X(3,1)*X(3,2)
tt443 = tt269+tt332
tt444 = tt443*tt9
tt445 = tt435*X(3,1)
tt446 = tt436*X(3,2)+tt445
tt447 = 2*tt230*X(1,3)
tt448 = 2*tt48*X(1,3)
tt449 = -2*X(2,1)*X(2,2)
tt450 = 2*tt6
tt451 = tt98+tt65
tt452 = tt69*tt11+tt451*X(3,2)*X(3,3)+tt69*tt9+tt69*tt8+(tt450+tt&
&449+tt448+tt232+tt166)*X(2,3)+tt91+tt310+(tt324+tt447+tt26+tt200)&
&*X(2,2)+(tt295+tt294+tt293)*X(2,1)
tt453 = tt452*X(3,4)+tt446*tt11+(tt444+tt442+tt440+tt438+tt437)*X&
&(3,3)+tt436*tt87+tt435*X(3,1)*tt9+(tt434+tt432+tt431+tt326+tt430)&
&*X(3,2)+(tt429+tt427+tt426+tt90+tt424)*X(3,1)
tt454 = tt14*tt453-(tt139*tt336*tt249)/2.0E+0
tt455 = tt77*tt454*tt252-tt14*tt423*tt249*tt252
tt456 = 6*tt1*tt2*tt255*tt455
tt457 = (-tt139*tt340*tt250)+2*tt109*tt453*tt249+2*tt423*tt77
tt458 = tt342+tt343
tt459 = tt458*X(3,4)+(tt211+tt345)*X(3,3)+(tt344+tt346)*X(3,2)
tt460 = -tt60*X(2,4)
tt461 = -(X(3,2)-X(3,3))*X(3,4)
tt462 = -tt14*(tt461+tt351+tt9+tt460+tt349+tt6)*tt249*tt252
tt463 = X(2,4)+tt95+X(2,2)
tt464 = tt85+tt105+tt74
tt465 = X(2,4)+tt63+X(2,1)
tt466 = tt95+X(2,2)+X(2,1)
tt467 = (-(tt75*X(3,3)+tt466*X(3,2)+tt106*X(3,1))*X(3,4))-(tt465*&
&X(3,2)+tt445)*X(3,3)-tt464*tt9-tt463*X(3,1)*X(3,2)-tt50*X(2,4)-2*&
&tt41*X(2,3)-tt37*X(2,2)-tt29*X(2,1)
tt468 = 2*tt130*X(2,2)*X(2,3)
tt469 = 3*tt132*tt8
tt470 = 2*tt27*X(2,3)+tt343
tt471 = tt470*X(2,4)
tt472 = 2*tt162*X(2,3)
tt473 = 3*tt155*tt8
tt474 = tt173*X(2,2)
tt475 = 2*tt170*X(2,3)+tt474
tt476 = tt475*X(2,4)
tt477 = 2*tt198*X(2,3)
tt478 = 2*tt39*X(2,3)+tt205*X(2,2)
tt479 = tt478*X(2,4)
tt480 = tt210*X(3,1)*X(3,2)
tt481 = tt160*tt9
tt482 = tt155*X(3,2)+tt132*X(3,1)
tt483 = tt48*tt11+tt246*X(3,2)*X(3,3)+tt48*tt9+3*tt48*tt8+2*tt241&
&*X(2,3)+tt238+tt236+tt234+tt233+tt118+tt165
tt484 = tt483*X(3,4)+tt482*tt11+(tt481+tt480+tt479+tt477+tt194+tt&
&192)*X(3,3)+tt155*tt87+tt132*X(3,1)*tt9+(tt476+tt473+tt472+tt156+&
&tt154+tt113+tt151+tt150+tt149)*X(3,2)+(tt471+tt469+tt468+tt126+tt&
&122+tt121+tt120+tt119)*X(3,1)
tt485 = tt14*tt484-(tt106*tt336*tt249)/2.0E+0
tt486 = tt77*tt485*tt252-tt14*tt467*tt249*tt252
tt487 = 6*tt1*tt2*tt255*tt486
tt488 = (-tt106*tt340*tt250)+2*tt109*tt484*tt249+2*tt467*tt77
tt489 = -4*X(2,2)*X(2,3)
tt490 = -2*tt9
tt491 = 6*X(2,2)*X(2,3)
tt492 = -3*tt8
tt493 = 2*X(3,2)*X(3,3)
tt494 = -tt11
tt495 = (tt494+tt493+tt396+tt492+tt491+tt103+tt16+tt81+tt26)*X(3,&
&4)+X(3,2)*tt11+(tt490+tt393+tt7+tt450)*X(3,3)+tt87+(tt395+tt400+t&
&t489+tt6+tt5+tt4+tt3)*X(3,2)
tt496 = -((-tt27*X(2,4))-2*tt132*X(2,3)-tt253*X(2,2))*tt14*tt249*&
&tt252
tt497 = (-tt76*X(3,4))-2*tt70*X(3,3)-tt67*X(3,2)-tt61*X(3,1)
tt498 = 2*tt243*X(3,3)+tt247*X(3,2)
tt499 = tt498*X(3,4)+3*tt226*tt11+2*tt225*X(3,3)+tt219+tt213+tt20&
&8+tt199+tt195+tt190+tt188+tt186+tt183
tt500 = -2*X(3,2)
tt501 = 2*X(3,3)+tt500
tt502 = tt14*tt499-(tt501*tt336*tt249)/2.0E+0
tt503 = tt77*tt502*tt252-tt14*tt497*tt249*tt252
tt504 = 6*tt1*tt2*tt255*tt503
tt505 = (-tt501*tt340*tt250)+2*tt109*tt499*tt249+2*tt497*tt77
tt506 = (2*tt60*X(3,3)+tt106*X(3,2))*X(3,4)+3*tt100*tt11+2*tt99*X&
&(3,2)*X(3,3)+tt97+tt94+tt93+tt92+tt91+tt89
tt507 = -tt14*((-tt27*X(3,4))-2*tt132*X(3,3)-tt253*X(3,2))*tt249*&
&tt252
tt508 = (tt398+tt5+tt4+tt3)*X(2,3)
tt509 = -3*X(2,2)*tt8
tt510 = (tt59+tt449+tt16+tt81+tt26)*X(2,3)
tt511 = tt98+X(2,1)
tt512 = tt511*tt8
tt513 = tt88*X(2,1)
tt514 = (tt425+tt329)*X(2,3)
tt515 = tt75*tt8
tt516 = tt358*X(3,1)*X(3,2)
tt517 = tt105+X(2,2)-3*X(2,1)
tt518 = tt517*tt9
tt519 = tt56+tt63+3*X(2,1)
tt520 = tt519*X(3,2)+tt54*X(3,1)
tt521 = tt75*tt101+tt520*tt11+(tt518+tt516+tt515+tt514+tt90+tt326&
&+tt102+tt513)*X(3,3)+tt57*tt87+tt54*X(3,1)*tt9+(tt104+tt512+tt510&
&+tt310+tt43)*X(3,2)+(tt80+tt509+tt508+tt91+tt89)*X(3,1)
tt522 = tt27*X(3,1)
tt523 = (-tt39*tt11)-(tt35*X(3,2)+tt522)*X(3,3)-tt23*tt9-tt18*X(3&
&,1)*X(3,2)-tt39*tt8-(tt35*X(2,2)+tt27*X(2,1))*X(2,3)-tt23*tt6-tt1&
&8*X(2,1)*X(2,2)
tt524 = tt14*tt521*tt77*tt252-tt14*tt523*tt249*tt252
tt525 = 6*tt1*tt2*tt524*tt255
tt526 = 2*tt109*tt521*tt249+2*tt523*tt77
tt527 = -(tt494+tt493+tt396+tt52+tt82+tt59)*tt14*tt249*tt252
tt528 = tt139*X(3,1)*X(3,2)
tt529 = tt217*tt9
tt530 = tt223*X(3,2)+tt522
tt531 = tt39*tt101+tt530*tt11+(tt529+tt528+tt207+tt206+tt204+tt20&
&3+tt202+tt119+tt149)*X(3,3)+tt170*tt87+tt27*X(3,1)*tt9+tt176*X(3,&
&2)+tt142*X(3,1)
tt532 = (-tt69*tt11)-(tt66*X(3,2)+tt60*X(3,1))*X(3,3)-tt57*tt9-tt&
&54*X(3,1)*X(3,2)-tt50*X(2,3)-tt46*X(2,2)-tt42*X(2,1)
tt533 = tt14*tt531*tt77*tt252-tt14*tt532*tt249*tt252
tt534 = 6*tt1*tt2*tt533*tt255
tt535 = 2*tt109*tt531*tt249+2*tt532*tt77
tt536 = -tt87
tt537 = tt14*(tt101-3*X(3,2)*tt11+(3*tt9+tt8+tt7+tt6+tt5+tt4+tt3)&
&*X(3,3)+tt536+tt83*X(3,2))*tt77*tt252
tt538 = -((-tt27*X(2,3))-tt18*X(2,2))*tt14*tt249*tt252
tt539 = (-tt76*X(3,3))-tt73*X(3,2)-tt71*X(3,1)
tt540 = tt14*tt248*tt77*tt252-tt539*tt14*tt249*tt252
tt541 = 6*tt1*tt2*tt540*tt255
tt542 = 2*tt109*tt248*tt249+2*tt539*tt77
tt543 = tt14*tt107*tt77*tt252
tt544 = -((-tt27*X(3,3))-tt18*X(3,2))*tt14*tt249*tt252
tt545 = 6*tt1*tt2*tt277*tt272
tt546 = 6*tt1*tt2*tt272*tt338
tt547 = (tt494+tt493+tt396+tt52+tt82+tt59+tt324+tt323+tt288)*X(3,&
&4)+tt101-2*X(3,2)*tt11+(tt9+tt8+tt7+tt6+tt286+tt5+tt306+tt293)*X(&
&3,3)+tt299*X(3,2)
tt548 = -((-tt235*X(2,4))-tt280*X(2,3)-tt116*X(2,2))*tt14*tt249*t&
&t252
tt549 = 6*tt1*tt2*tt272*tt389
tt550 = tt191*X(2,3)
tt551 = (tt342+2*tt18*X(2,2))*X(3,4)+(tt550+2*tt187*X(2,2))*X(3,3&
&)+(tt153*X(2,3)+2*tt147*X(2,2))*X(3,2)
tt552 = -X(1,2)*X(1,3)
tt553 = -tt18*X(1,4)
tt554 = -tt14*(tt353+tt11+tt351+tt553+tt5+tt552)*tt249*tt252
tt555 = 6*tt1*tt2*tt272*tt412
tt556 = (tt235*X(3,3)+2*tt18*X(3,2))*X(3,4)+tt220*tt11+2*tt214*X(&
&3,2)*X(3,3)+3*tt147*tt9+tt261+tt260+tt259+tt122+tt136+tt145+tt144
tt557 = -tt14*((-tt54*X(3,4))-tt270*X(3,3)-2*tt268*X(3,2))*tt249*&
&tt252
tt558 = 6*tt1*tt2*tt272*tt455
tt559 = 2*tt9
tt560 = (tt11+tt10+tt9+tt8+tt7+tt6+tt295+tt294+tt293)*X(3,4)-X(3,&
&2)*tt11+(tt559+tt290+tt81+tt232)*X(3,3)+tt536+(tt52+tt82+tt59+tt2&
&86+tt324+tt289+tt26)*X(3,2)
tt561 = -((-tt139*X(2,4))-tt187*X(2,3)-tt417*X(2,2))*tt14*tt249*t&
&t252
tt562 = 6*tt1*tt2*tt272*tt486
tt563 = tt153*X(2,2)
tt564 = (2*tt18*X(2,3)+tt345)*X(3,4)+(2*tt187*X(2,3)+tt191*X(2,2)&
&)*X(3,3)+(2*tt147*X(2,3)+tt563)*X(3,2)
tt565 = -tt27*X(1,4)
tt566 = -tt14*(tt461+tt351+tt9+tt565+tt552+tt3)*tt249*tt252
tt567 = 6*tt1*tt2*tt272*tt503
tt568 = (2*tt18*X(3,3)+tt235*X(3,2))*X(3,4)+3*tt187*tt11+2*tt220*&
&X(3,2)*X(3,3)+tt265+tt264+tt263+tt262+tt113+tt182+tt181+tt118
tt569 = -tt14*((-tt60*X(3,4))-2*tt100*X(3,3)-tt270*X(3,2))*tt249*&
&tt252
tt570 = 6*tt1*tt2*tt524*tt272
tt571 = tt14*(tt397+3*X(3,2)*tt11+((-3*tt9)+tt52+tt82+tt59+tt16+t&
&t81+tt26)*X(3,3)+tt87+tt71*X(3,2))*tt77*tt252
tt572 = 6*tt1*tt2*tt533*tt272
tt573 = -(tt494+tt493+tt396+tt16+tt81+tt26)*tt14*tt249*tt252
tt574 = 6*tt1*tt2*tt540*tt272
tt575 = tt14*tt266*tt77*tt252
tt576 = -((-tt60*X(3,3))-tt54*X(3,2))*tt14*tt249*tt252
tt577 = 6*tt1*tt2*tt277*tt338
tt578 = tt268*tt11+tt316*X(3,2)*X(3,3)+tt268*tt9+tt297+tt104+tt29&
&2+tt291+tt287
tt579 = -tt14*((-tt235*X(3,4))-tt280*X(3,3)-tt116*X(3,2))*tt249*t&
&t252
tt580 = 6*tt1*tt2*tt277*tt389
tt581 = tt116*tt11+tt209*X(3,2)*X(3,3)+tt116*tt9+tt366+tt364+tt36&
&3+tt362+tt113+tt112+tt111+tt110
tt582 = -tt14*((-tt358*X(3,4))-tt355*X(3,3)-tt86*X(3,2))*tt249*tt&
&252
tt583 = 6*tt1*tt2*tt277*tt412
tt584 = tt212*X(3,3)+2*tt178*X(3,2)
tt585 = -(tt350+tt8+tt349+tt553+tt5+tt552)*tt14*tt249*tt252
tt586 = 6*tt1*tt2*tt277*tt455
tt587 = tt435*tt11+tt441*X(3,2)*X(3,3)+tt435*tt9+tt429+tt427+tt42&
&6+tt90+tt424
tt588 = -tt14*((-tt139*X(3,4))-tt187*X(3,3)-tt417*X(3,2))*tt249*t&
&t252
tt589 = 6*tt1*tt2*tt277*tt486
tt590 = tt132*tt11+tt210*X(3,2)*X(3,3)+tt132*tt9+tt471+tt469+tt46&
&8+tt126+tt122+tt121+tt120+tt119
tt591 = -tt14*((-tt106*X(3,4))-tt435*X(3,3)-tt463*X(3,2))*tt249*t&
&t252
tt592 = 6*tt1*tt2*tt277*tt503
tt593 = 2*tt178*X(3,3)+tt212*X(3,2)
tt594 = -(tt460+tt349+tt6+tt565+tt552+tt3)*tt14*tt249*tt252
tt595 = 6*tt1*tt2*tt524*tt277
tt596 = tt14*(tt54*tt11+tt358*X(3,2)*X(3,3)+tt54*tt9+tt80+tt509+t&
&t508+tt91+tt89)*tt77*tt252
tt597 = 6*tt1*tt2*tt533*tt277
tt598 = tt14*(tt27*tt11+tt139*X(3,2)*X(3,3)+tt27*tt9+tt141+tt140+&
&tt137+tt136+tt135+tt134+tt118)*tt77*tt252
tt599 = 6*tt1*tt2*tt540*tt277
tt600 = -tt83*tt14*tt249*tt252
tt601 = 2*X(2,1)*X(2,3)
tt602 = 2*X(3,1)*X(3,3)
tt603 = -2*tt11
tt604 = 2*X(3,1)
tt605 = -(tt409+tt604)*X(3,4)
tt606 = -6*X(1,2)
tt607 = 4*X(1,3)
tt608 = 6*X(1,2)
tt609 = -6*X(1,3)
tt610 = tt609+tt608
tt611 = tt610*X(2,4)
tt612 = tt193*X(2,3)
tt613 = -4*X(1,3)
tt614 = 6*X(1,3)
tt615 = tt614+tt606
tt616 = tt615*X(2,1)
tt617 = 1/tt13**5
tt618 = -tt336*tt249
tt619 = 6*tt1*tt2*tt338*tt389
tt620 = -((-tt170*X(2,4))-tt281*X(2,3)-tt116*X(2,1))*tt14*tt249*t&
&t252
tt621 = 2*tt8
tt622 = tt393-2*X(2,1)*X(2,3)+tt329+tt290+tt308+tt200
tt623 = tt77*((3.0E+0*tt235*tt358*tt617*tt249)/4.0E+0-(tt358*tt33&
&6*tt335)/2.0E+0-(tt235*tt336*tt387)/2.0E+0+tt14*((tt392+(tt98+tt3&
&32)*X(2,3)+tt449+tt298+tt325+tt166)*X(3,4)+tt622*X(3,3)+(tt395+tt&
&621+tt7+tt286+tt285+tt81)*X(3,1)))*tt252
tt624 = 6*tt1*tt2*tt338*tt412
tt625 = -tt14*((-tt170*X(3,4))-tt281*X(3,3)-tt116*X(3,1))*tt249*t&
&t252
tt626 = tt77*((3.0E+0*tt235*tt410*tt617*tt249)/4.0E+0-(tt410*tt33&
&6*tt335)/2.0E+0-(tt235*tt336*tt407)/2.0E+0+tt14*((tt333*X(3,3)+2*&
&tt331*X(3,2))*X(3,4)+tt321*tt11+(2*tt318*X(3,2)+tt316*X(3,1))*X(3&
&,3)+2*tt268*X(3,1)*X(3,2)+tt305+tt303+tt300))*tt252
tt627 = 6*tt1*tt2*tt338*tt455
tt628 = -tt64*X(2,3)
tt629 = -tt66*X(2,4)
tt630 = -X(3,1)*X(3,2)
tt631 = -(tt500+X(3,1))*X(3,3)
tt632 = -(X(3,3)+X(3,2)-2*X(3,1))*X(3,4)
tt633 = -tt14*(tt632+tt631+tt630+tt629+tt628+tt68)*tt249*tt252
tt634 = tt336*tt249
tt635 = tt77*(tt634+(3.0E+0*tt235*tt139*tt617*tt249)/4.0E+0-(tt13&
&9*tt336*tt335)/2.0E+0-(tt235*tt336*tt453)/2.0E+0+tt14*(((tt138+tt&
&327+tt172)*X(2,3)+(tt607+tt31+tt34)*X(2,2)+tt610*X(2,1))*X(3,4)+(&
&(tt123+tt313+tt34)*X(2,4)+tt221+(tt159+tt138+tt327)*X(2,1))*X(3,3&
&)+((tt613+tt128+tt172)*X(2,4)+tt215+(tt152+tt607+tt31)*X(2,1))*X(&
&3,2)+(tt615*X(2,4)+(tt152+tt123+tt313)*X(2,3)+(tt159+tt613+tt128)&
&*X(2,2))*X(3,1)))*tt252
tt636 = 6*tt1*tt2*tt338*tt486
tt637 = -((-tt283*X(2,4))-2*tt282*X(2,3)-tt281*X(2,2)-tt280*X(2,1&
&))*tt14*tt249*tt252
tt638 = -X(3,1)*tt9
tt639 = 2*X(3,1)*X(3,2)
tt640 = tt77*((3.0E+0*tt235*tt106*tt617*tt249)/4.0E+0-(tt106*tt33&
&6*tt335)/2.0E+0-(tt235*tt336*tt484)/2.0E+0+tt14*((tt11+tt10+tt9+t&
&t400+2*tt330*X(2,3)+tt6+tt329+tt5+tt328+tt293+tt200)*X(3,4)-X(3,1&
&)*tt11+(tt639+tt395+tt601+tt449)*X(3,3)+tt638+tt302*X(3,2)+(tt393&
&+tt492+tt394+tt59+tt290+tt16+tt289+tt288)*X(3,1)))*tt252
tt641 = 6*tt1*tt2*tt338*tt503
tt642 = -tt14*((-tt283*X(3,4))-2*tt282*X(3,3)-tt281*X(3,2)-tt280*&
&X(3,1))*tt249*tt252
tt643 = tt77*((3.0E+0*tt235*tt501*tt617*tt249)/4.0E+0-(tt501*tt33&
&6*tt335)/2.0E+0-(tt235*tt336*tt499)/2.0E+0+tt14*((2*tt331*X(3,3)+&
&tt333*X(3,2))*X(3,4)+3*tt318*tt11+2*tt322*X(3,3)+tt319+tt317+tt31&
&5+tt312+tt311+tt310+tt309+tt307))*tt252
tt644 = 6*tt1*tt2*tt524*tt338
tt645 = (tt345+tt420)*X(3,3)+(tt346+tt244)*X(3,2)+tt458*X(3,1)
tt646 = -tt72*X(2,3)
tt647 = X(3,1)*X(3,2)
tt648 = -(X(3,2)+X(3,1))*X(3,3)
tt649 = -tt14*(tt11+tt648+tt647+tt8+tt646+tt62)*tt249*tt252
tt650 = 6*tt1*tt2*tt533*tt338
tt651 = X(3,1)*tt9
tt652 = -2*X(3,1)*X(3,2)
tt653 = tt397+(tt408+X(3,1))*tt11+(tt396+tt652+tt52+tt82+tt59+tt1&
&6+tt314+tt288+tt166)*X(3,3)+tt651+tt304*X(3,2)+tt296*X(3,1)
tt654 = -((-tt283*X(2,3))-tt170*X(2,2)-tt235*X(2,1))*tt14*tt249*t&
&t252
tt655 = 6*tt1*tt2*tt540*tt338
tt656 = -((-tt283*X(3,3))-tt170*X(3,2)-tt235*X(3,1))*tt14*tt249*t&
&t252
tt657 = 2*tt27*X(2,4)
tt658 = 2*tt170*X(2,4)
tt659 = 2*tt39*X(2,4)
tt660 = 6*tt1*tt2*tt389*tt412
tt661 = -tt14*((-tt331*X(3,4))-tt356*X(3,3)-tt86*X(3,1))*tt249*tt&
&252
tt662 = tt77*((3.0E+0*tt358*tt410*tt617*tt249)/4.0E+0-(tt410*tt33&
&6*tt387)/2.0E+0-(tt358*tt336*tt407)/2.0E+0+tt14*((tt245*X(3,3)+2*&
&tt23*X(3,2))*X(3,4)+tt193*tt11+(2*tt189*X(3,2)+tt209*X(3,1))*X(3,&
&3)+2*tt116*X(3,1)*X(3,2)+tt373+tt371+tt370+tt367))*tt252
tt663 = 6*tt1*tt2*tt455*tt389
tt664 = -((-tt421*X(2,4))-tt419*X(2,3)-2*tt418*X(2,2)-tt417*X(2,1&
&))*tt14*tt249*tt252
tt665 = tt77*((3.0E+0*tt139*tt358*tt617*tt249)/4.0E+0-(tt358*tt33&
&6*tt453)/2.0E+0-(tt139*tt336*tt387)/2.0E+0+tt14*((tt494+tt493+tt3&
&96+tt52+(4*X(2,2)+tt65)*X(2,3)+tt103+tt329+tt324+tt447+tt26+tt200&
&)*X(3,4)+X(3,1)*tt11+(tt652+tt286+tt301+tt166)*X(3,3)+tt651+(tt39&
&5+tt601+tt449)*X(3,2)+(tt393+tt8+tt489+tt398+tt290+tt295+tt306+tt&
&3)*X(3,1)))*tt252
tt666 = 6*tt1*tt2*tt389*tt486
tt667 = -tt32*X(1,3)
tt668 = -tt35*X(1,4)
tt669 = -tt14*(tt632+tt631+tt630+tt668+tt667+tt38)*tt249*tt252
tt670 = tt77*(tt634+(3.0E+0*tt358*tt106*tt617*tt249)/4.0E+0-(tt10&
&6*tt336*tt387)/2.0E+0-(tt358*tt336*tt484)/2.0E+0+tt14*((2*tt240*X&
&(2,3)+tt385+tt244)*X(3,4)+(tt205*X(2,4)+2*tt189*X(2,3)+tt377+tt37&
&6)*X(3,3)+(tt173*X(2,4)+2*tt160*X(2,3)+tt369+tt368)*X(3,2)+(tt211&
&+2*tt130*X(2,3)+2*tt125*X(2,2))*X(3,1)))*tt252
tt671 = 6*tt1*tt2*tt389*tt503
tt672 = -tt14*((-tt360*X(3,4))-2*tt357*X(3,3)-tt356*X(3,2)-tt355*&
&X(3,1))*tt249*tt252
tt673 = tt77*((3.0E+0*tt358*tt501*tt617*tt249)/4.0E+0-(tt501*tt33&
&6*tt387)/2.0E+0-(tt358*tt336*tt499)/2.0E+0+tt14*((2*tt23*X(3,3)+t&
&t245*X(3,2))*X(3,4)+3*tt189*tt11+2*tt384*X(3,3)+tt383+tt382+tt381&
&+tt379+tt378+tt375+tt374+tt122+tt185+tt184+tt165))*tt252
tt674 = 6*tt1*tt2*tt524*tt389
tt675 = -X(3,1)
tt676 = tt101+(tt500+tt675)*tt11+(tt9+tt639+tt8+(tt332-4*X(2,2))*&
&X(2,3)+tt398+tt449+tt5+tt4+tt3)*X(3,3)+tt638+(tt621+(tt63+tt65)*X&
&(2,3)+tt329)*X(3,2)+(tt492+tt491+tt103+tt16+tt81+tt26)*X(3,1)
tt677 = -((-tt35*X(2,3))-2*tt23*X(2,2)-tt18*X(2,1))*tt14*tt249*tt&
&252
tt678 = 6*tt1*tt2*tt533*tt389
tt679 = tt380*X(3,3)+tt372*X(3,2)+tt365*X(3,1)
tt680 = -tt44*X(1,3)
tt681 = -tt14*(tt11+tt648+tt647+tt5+tt680+tt30)*tt249*tt252
tt682 = 6*tt1*tt2*tt540*tt389
tt683 = -((-tt360*X(3,3))-tt331*X(3,2)-tt358*X(3,1))*tt14*tt249*t&
&t252
tt684 = 2*tt243*X(3,4)
tt685 = 6*tt1*tt2*tt455*tt412
tt686 = -tt14*((-tt421*X(3,4))-tt419*X(3,3)-2*tt418*X(3,2)-tt417*&
&X(3,1))*tt249*tt252
tt687 = tt77*((3.0E+0*tt139*tt410*tt617*tt249)/4.0E+0-(tt410*tt33&
&6*tt453)/2.0E+0-(tt139*tt336*tt407)/2.0E+0+tt14*((tt451*X(3,3)+2*&
&tt69*X(3,2))*X(3,4)+tt436*tt11+(2*tt443*X(3,2)+tt441*X(3,1))*X(3,&
&3)+3*tt436*tt9+2*tt435*X(3,1)*X(3,2)+tt434+tt432+tt431+tt326+tt43&
&0))*tt252
tt688 = 6*tt1*tt2*tt486*tt412
tt689 = -tt14*((-tt466*X(3,4))-tt465*X(3,3)-2*tt464*X(3,2)-tt463*&
&X(3,1))*tt249*tt252
tt690 = tt77*((3.0E+0*tt106*tt410*tt617*tt249)/4.0E+0-(tt410*tt33&
&6*tt484)/2.0E+0-(tt106*tt336*tt407)/2.0E+0+tt14*((tt246*X(3,3)+2*&
&tt48*X(3,2))*X(3,4)+tt155*tt11+(2*tt160*X(3,2)+tt210*X(3,1))*X(3,&
&3)+3*tt155*tt9+2*tt132*X(3,1)*X(3,2)+tt476+tt473+tt472+tt156+tt15&
&4+tt113+tt151+tt150+tt149))*tt252
tt691 = 6*tt1*tt2*tt412*tt503
tt692 = -(tt629+tt628+tt68+tt668+tt667+tt38)*tt14*tt249*tt252
tt693 = tt77*(tt634+(3.0E+0*tt410*tt501*tt617*tt249)/4.0E+0-(tt41&
&0*tt336*tt499)/2.0E+0-(tt501*tt336*tt407)/2.0E+0+tt14*(tt247*X(3,&
&4)+2*tt224*X(3,3)+tt405+tt404))*tt252
tt694 = 6*tt1*tt2*tt524*tt412
tt695 = tt519*tt11+(2*tt517*X(3,2)+tt359)*X(3,3)+3*tt57*tt9+2*tt5&
&4*X(3,1)*X(3,2)+tt104+tt512+tt510+tt310+tt43
tt696 = -((-tt35*X(3,3))-2*tt23*X(3,2)-tt18*X(3,1))*tt14*tt249*tt&
&252
tt697 = 6*tt1*tt2*tt533*tt412
tt698 = tt223*tt11+(2*tt217*X(3,2)+tt422)*X(3,3)+3*tt170*tt9+2*tt&
&27*X(3,1)*X(3,2)+tt175+tt174+tt171+tt112+tt169+tt167+tt165
tt699 = -((-tt66*X(3,3))-2*tt57*X(3,2)-tt54*X(3,1))*tt14*tt249*tt&
&252
tt700 = 6*tt1*tt2*tt540*tt412
tt701 = -(tt8+tt646+tt62+tt5+tt680+tt30)*tt14*tt249*tt252
tt702 = -(tt500+tt604)*X(3,4)
tt703 = 6*tt1*tt2*tt455*tt486
tt704 = -((-tt48*X(2,4))-tt419*X(2,2)-tt187*X(2,1))*tt14*tt249*tt&
&252
tt705 = tt77*((3.0E+0*tt139*tt106*tt617*tt249)/4.0E+0-(tt106*tt33&
&6*tt453)/2.0E+0-(tt139*tt336*tt484)/2.0E+0+tt14*((2*tt69*X(2,3)+t&
&t450+tt449+tt448+tt232+tt166)*X(3,4)+tt622*X(3,2)+(tt395+tt82+tt4&
&25+tt286+tt4+tt201)*X(3,1)))*tt252
tt706 = 6*tt1*tt2*tt455*tt503
tt707 = -tt14*((-tt48*X(3,4))-tt419*X(3,2)-tt187*X(3,1))*tt249*tt&
&252
tt708 = tt77*((3.0E+0*tt139*tt501*tt617*tt249)/4.0E+0-(tt501*tt33&
&6*tt453)/2.0E+0-(tt139*tt336*tt499)/2.0E+0+tt14*((2*tt69*X(3,3)+t&
&t451*X(3,2))*X(3,4)+2*tt446*X(3,3)+tt444+tt442+tt440+tt438+tt437)&
&)*tt252
tt709 = 6*tt1*tt2*tt524*tt455
tt710 = (tt343+tt244)*X(3,3)+(tt342+tt420)*X(3,2)+tt347*X(3,1)
tt711 = -tt75*X(2,3)
tt712 = X(3,2)+tt675
tt713 = -tt712*X(3,3)
tt714 = -(tt713+tt9+tt630+tt711+tt6+tt68)*tt14*tt249*tt252
tt715 = 6*tt1*tt2*tt533*tt455
tt716 = tt712*tt11+(tt490+tt639+tt439+tt201+tt200)*X(3,3)+tt87+tt&
&638+tt433*X(3,2)+tt428*X(3,1)
tt717 = -((-tt48*X(2,3))-tt421*X(2,2)-tt139*X(2,1))*tt14*tt249*tt&
&252
tt718 = 6*tt1*tt2*tt540*tt455
tt719 = -((-tt48*X(3,3))-tt421*X(3,2)-tt139*X(3,1))*tt14*tt249*tt&
&252
tt720 = 6*tt1*tt2*tt486*tt503
tt721 = -tt14*((-tt75*X(3,4))-tt465*X(3,2)-tt435*X(3,1))*tt249*tt&
&252
tt722 = tt77*((3.0E+0*tt106*tt501*tt617*tt249)/4.0E+0-(tt501*tt33&
&6*tt484)/2.0E+0-(tt106*tt336*tt499)/2.0E+0+tt14*((2*tt48*X(3,3)+t&
&t246*X(3,2))*X(3,4)+2*tt482*X(3,3)+tt481+tt480+tt479+tt477+tt194+&
&tt192))*tt252
tt723 = 6*tt1*tt2*tt524*tt486
tt724 = (tt352+X(3,1))*tt11+(tt559+tt652+2*tt75*X(2,3)+tt425+tt32&
&9)*X(3,3)+tt536+tt651+(tt492+2*tt511*X(2,3)+tt59+tt449+tt16+tt81+&
&tt26)*X(3,2)+(tt400+tt399+tt398+tt5+tt4+tt3)*X(3,1)
tt725 = -((-2*tt39*X(2,3))-tt35*X(2,2)-tt27*X(2,1))*tt14*tt249*tt&
&252
tt726 = 6*tt1*tt2*tt533*tt486
tt727 = tt478*X(3,3)+tt475*X(3,2)+tt470*X(3,1)
tt728 = -tt48*X(1,3)
tt729 = -(tt713+tt9+tt630+tt728+tt3+tt38)*tt14*tt249*tt252
tt730 = 6*tt1*tt2*tt540*tt486
tt731 = -((-tt75*X(3,3))-tt466*X(3,2)-tt106*X(3,1))*tt14*tt249*tt&
&252
tt732 = 6*tt1*tt2*tt524*tt503
tt733 = 3*tt75*tt11+2*tt520*X(3,3)+tt518+tt516+tt515+tt514+tt90+t&
&t326+tt102+tt513
tt734 = -((-2*tt39*X(3,3))-tt35*X(3,2)-tt27*X(3,1))*tt14*tt249*tt&
&252
tt735 = 6*tt1*tt2*tt533*tt503
tt736 = 3*tt39*tt11+2*tt530*X(3,3)+tt529+tt528+tt207+tt206+tt204+&
&tt203+tt202+tt119+tt149
tt737 = -((-2*tt69*X(3,3))-tt66*X(3,2)-tt60*X(3,1))*tt14*tt249*tt&
&252
tt738 = 6*tt1*tt2*tt540*tt503
tt739 = -(tt711+tt6+tt68+tt728+tt3+tt38)*tt14*tt249*tt252
tt740 = 6*tt1*tt2*tt524*tt533
tt741 = 6*tt1*tt2*tt524*tt540
tt742 = 6*tt1*tt2*tt533*tt540
hes(1,1) = 6*tt1*tt2*(tt14*tt254*tt249*tt256*tt257-tt14*tt77*tt10&
&8*tt256*tt257)*tt258+6*tt1*tt2*tt255**2
hes(1,2) = 6*tt1*tt2*(tt14*tt271*tt108*tt252-tt14*tt254*tt267*tt2&
&52+tt14*tt254*tt249*tt274*tt257-tt14*tt77*tt108*tt274*tt257)*tt25&
&8+tt273
hes(1,3) = 6*tt1*tt2*(tt14*tt276*tt108*tt252-tt14*tt275*tt254*tt2&
&52+tt14*tt254*tt249*tt279*tt257-tt14*tt77*tt108*tt279*tt257)*tt25&
&8+tt278
hes(1,4) = 6*tt1*tt2*((tt235*tt336*tt254*tt249*tt252)/2.0E+0+tt35&
&4-tt14*tt254*tt335*tt252-(tt235*tt336*tt77*tt108*tt252)/2.0E+0+tt&
&14*tt284*tt108*tt252+tt14*tt348*tt77*tt252+tt14*tt254*tt249*tt341&
&*tt257-tt14*tt77*tt108*tt341*tt257)*tt258+tt339
hes(1,5) = 6*tt1*tt2*((tt358*tt336*tt254*tt249*tt252)/2.0E+0+tt40&
&2-(tt358*tt336*tt77*tt108*tt252)/2.0E+0+tt14*tt361*tt108*tt252-tt&
&14*tt254*tt387*tt252+tt14*tt77*tt401*tt252+tt14*tt254*tt249*tt391&
&*tt257-tt14*tt77*tt108*tt391*tt257)*tt258+tt390
hes(1,6) = 6*tt1*tt2*((tt410*tt336*tt254*tt249*tt252)/2.0E+0+tt41&
&6-(tt410*tt336*tt77*tt108*tt252)/2.0E+0+tt14*tt403*tt108*tt252+tt&
&14*tt415*tt77*tt252-tt14*tt254*tt407*tt252+tt14*tt254*tt249*tt414&
&*tt257-tt14*tt77*tt108*tt414*tt257)*tt258+tt413
hes(1,7) = 6*tt1*tt2*((tt139*tt336*tt254*tt249*tt252)/2.0E+0+tt46&
&2-(tt139*tt336*tt77*tt108*tt252)/2.0E+0+tt14*tt423*tt108*tt252-tt&
&14*tt254*tt453*tt252+tt14*tt459*tt77*tt252+tt14*tt254*tt249*tt457&
&*tt257-tt14*tt77*tt108*tt457*tt257)*tt258+tt456
hes(1,8) = 6*tt1*tt2*((tt106*tt336*tt254*tt249*tt252)/2.0E+0+tt49&
&6-(tt106*tt336*tt77*tt108*tt252)/2.0E+0+tt14*tt467*tt108*tt252-tt&
&14*tt254*tt484*tt252+tt14*tt77*tt495*tt252+tt14*tt254*tt249*tt488&
&*tt257-tt14*tt77*tt108*tt488*tt257)*tt258+tt487
hes(1,9) = 6*tt1*tt2*((tt501*tt336*tt254*tt249*tt252)/2.0E+0+tt50&
&7-(tt501*tt336*tt77*tt108*tt252)/2.0E+0+tt14*tt497*tt108*tt252+tt&
&14*tt506*tt77*tt252-tt14*tt254*tt499*tt252+tt14*tt254*tt249*tt505&
&*tt257-tt14*tt77*tt108*tt505*tt257)*tt258+tt504
hes(1,10) = 6*tt1*tt2*(tt527+tt14*tt523*tt108*tt252-tt14*tt521*tt&
&254*tt252+tt14*tt254*tt249*tt526*tt257-tt14*tt77*tt108*tt526*tt25&
&7)*tt258+tt525
hes(1,11) = 6*tt1*tt2*(tt538+tt14*tt532*tt108*tt252+tt537-tt14*tt&
&531*tt254*tt252+tt14*tt254*tt249*tt535*tt257-tt14*tt77*tt108*tt53&
&5*tt257)*tt258+tt534
hes(1,12) = 6*tt1*tt2*(tt544+tt539*tt14*tt108*tt252+tt543-tt14*tt&
&248*tt254*tt252+tt14*tt254*tt249*tt542*tt257-tt14*tt77*tt108*tt54&
&2*tt257)*tt258+tt541
hes(2,1) = 6*tt1*tt2*((-tt14*tt271*tt108*tt252)+tt14*tt254*tt267*&
&tt252+tt14*tt271*tt249*tt256*tt257-tt14*tt77*tt267*tt256*tt257)*t&
&t258+tt273
hes(2,2) = 6*tt1*tt2*(tt14*tt271*tt249*tt274*tt257-tt14*tt77*tt26&
&7*tt274*tt257)*tt258+6*tt1*tt2*tt272**2
hes(2,3) = 6*tt1*tt2*(tt14*tt276*tt267*tt252-tt14*tt275*tt271*tt2&
&52+tt14*tt271*tt249*tt279*tt257-tt14*tt77*tt267*tt279*tt257)*tt25&
&8+tt545
hes(2,4) = 6*tt1*tt2*((tt235*tt336*tt271*tt249*tt252)/2.0E+0+tt54&
&8-tt14*tt271*tt335*tt252-(tt235*tt336*tt77*tt267*tt252)/2.0E+0+tt&
&14*tt284*tt267*tt252+tt14*tt77*tt547*tt252+tt14*tt271*tt249*tt341&
&*tt257-tt14*tt77*tt267*tt341*tt257)*tt258+tt546
hes(2,5) = 6*tt1*tt2*((tt358*tt336*tt271*tt249*tt252)/2.0E+0+tt55&
&4-(tt358*tt336*tt77*tt267*tt252)/2.0E+0+tt14*tt361*tt267*tt252-tt&
&14*tt271*tt387*tt252+tt14*tt551*tt77*tt252+tt14*tt271*tt249*tt391&
&*tt257-tt14*tt77*tt267*tt391*tt257)*tt258+tt549
hes(2,6) = 6*tt1*tt2*((tt410*tt336*tt271*tt249*tt252)/2.0E+0+tt55&
&7-(tt410*tt336*tt77*tt267*tt252)/2.0E+0+tt14*tt403*tt267*tt252+tt&
&14*tt556*tt77*tt252-tt14*tt271*tt407*tt252+tt14*tt271*tt249*tt414&
&*tt257-tt14*tt77*tt267*tt414*tt257)*tt258+tt555
hes(2,7) = 6*tt1*tt2*((tt139*tt336*tt271*tt249*tt252)/2.0E+0+tt56&
&1-tt14*tt271*tt453*tt252-(tt139*tt336*tt77*tt267*tt252)/2.0E+0+tt&
&14*tt423*tt267*tt252+tt14*tt77*tt560*tt252+tt14*tt271*tt249*tt457&
&*tt257-tt14*tt77*tt267*tt457*tt257)*tt258+tt558
hes(2,8) = 6*tt1*tt2*((tt106*tt336*tt271*tt249*tt252)/2.0E+0+tt56&
&6-(tt106*tt336*tt77*tt267*tt252)/2.0E+0+tt14*tt467*tt267*tt252-tt&
&14*tt271*tt484*tt252+tt14*tt564*tt77*tt252+tt14*tt271*tt249*tt488&
&*tt257-tt14*tt77*tt267*tt488*tt257)*tt258+tt562
hes(2,9) = 6*tt1*tt2*((tt501*tt336*tt271*tt249*tt252)/2.0E+0+tt56&
&9-(tt501*tt336*tt77*tt267*tt252)/2.0E+0+tt14*tt497*tt267*tt252+tt&
&14*tt568*tt77*tt252-tt14*tt271*tt499*tt252+tt14*tt271*tt249*tt505&
&*tt257-tt14*tt77*tt267*tt505*tt257)*tt258+tt567
hes(2,10) = 6*tt1*tt2*(tt538+tt14*tt523*tt267*tt252+tt571-tt14*tt&
&521*tt271*tt252+tt14*tt271*tt249*tt526*tt257-tt14*tt77*tt267*tt52&
&6*tt257)*tt258+tt570
hes(2,11) = 6*tt1*tt2*(tt573+tt14*tt532*tt267*tt252-tt14*tt531*tt&
&271*tt252+tt14*tt271*tt249*tt535*tt257-tt14*tt77*tt267*tt535*tt25&
&7)*tt258+tt572
hes(2,12) = 6*tt1*tt2*(tt576+tt539*tt14*tt267*tt252+tt575-tt14*tt&
&248*tt271*tt252+tt14*tt271*tt249*tt542*tt257-tt14*tt77*tt267*tt54&
&2*tt257)*tt258+tt574
hes(3,1) = 6*tt1*tt2*((-tt14*tt276*tt108*tt252)+tt14*tt275*tt254*&
&tt252+tt14*tt276*tt249*tt256*tt257-tt14*tt275*tt77*tt256*tt257)*t&
&t258+tt278
hes(3,2) = 6*tt1*tt2*((-tt14*tt276*tt267*tt252)+tt14*tt275*tt271*&
&tt252+tt14*tt276*tt249*tt274*tt257-tt14*tt275*tt77*tt274*tt257)*t&
&t258+tt545
hes(3,3) = 6*tt1*tt2*(tt14*tt276*tt249*tt279*tt257-tt14*tt275*tt7&
&7*tt279*tt257)*tt258+6*tt1*tt2*tt277**2
hes(3,4) = 6*tt1*tt2*((tt235*tt336*tt276*tt249*tt252)/2.0E+0+tt57&
&9-tt14*tt276*tt335*tt252-(tt235*tt336*tt275*tt77*tt252)/2.0E+0+tt&
&14*tt578*tt77*tt252+tt14*tt275*tt284*tt252+tt14*tt276*tt249*tt341&
&*tt257-tt14*tt275*tt77*tt341*tt257)*tt258+tt577
hes(3,5) = 6*tt1*tt2*((tt358*tt336*tt276*tt249*tt252)/2.0E+0+tt58&
&2-tt14*tt276*tt387*tt252-(tt358*tt336*tt275*tt77*tt252)/2.0E+0+tt&
&14*tt581*tt77*tt252+tt14*tt275*tt361*tt252+tt14*tt276*tt249*tt391&
&*tt257-tt14*tt275*tt77*tt391*tt257)*tt258+tt580
hes(3,6) = 6*tt1*tt2*((tt410*tt336*tt276*tt249*tt252)/2.0E+0+tt58&
&5+tt584*tt14*tt77*tt252-(tt410*tt336*tt275*tt77*tt252)/2.0E+0-tt1&
&4*tt276*tt407*tt252+tt14*tt275*tt403*tt252+tt14*tt276*tt249*tt414&
&*tt257-tt14*tt275*tt77*tt414*tt257)*tt258+tt583
hes(3,7) = 6*tt1*tt2*((tt139*tt336*tt276*tt249*tt252)/2.0E+0+tt58&
&8-tt14*tt276*tt453*tt252-(tt139*tt336*tt275*tt77*tt252)/2.0E+0+tt&
&14*tt587*tt77*tt252+tt14*tt275*tt423*tt252+tt14*tt276*tt249*tt457&
&*tt257-tt14*tt275*tt77*tt457*tt257)*tt258+tt586
hes(3,8) = 6*tt1*tt2*((tt106*tt336*tt276*tt249*tt252)/2.0E+0+tt59&
&1-tt14*tt276*tt484*tt252-(tt106*tt336*tt275*tt77*tt252)/2.0E+0+tt&
&14*tt590*tt77*tt252+tt14*tt275*tt467*tt252+tt14*tt276*tt249*tt488&
&*tt257-tt14*tt275*tt77*tt488*tt257)*tt258+tt589
hes(3,9) = 6*tt1*tt2*((tt501*tt336*tt276*tt249*tt252)/2.0E+0+tt59&
&4+tt593*tt14*tt77*tt252-(tt501*tt336*tt275*tt77*tt252)/2.0E+0-tt1&
&4*tt276*tt499*tt252+tt14*tt275*tt497*tt252+tt14*tt276*tt249*tt505&
&*tt257-tt14*tt275*tt77*tt505*tt257)*tt258+tt592
hes(3,10) = 6*tt1*tt2*(tt544+tt596-tt14*tt521*tt276*tt252+tt14*tt&
&523*tt275*tt252+tt14*tt276*tt249*tt526*tt257-tt14*tt275*tt77*tt52&
&6*tt257)*tt258+tt595
hes(3,11) = 6*tt1*tt2*(tt576+tt598-tt14*tt531*tt276*tt252+tt14*tt&
&532*tt275*tt252+tt14*tt276*tt249*tt535*tt257-tt14*tt275*tt77*tt53&
&5*tt257)*tt258+tt597
hes(3,12) = 6*tt1*tt2*(tt600-tt14*tt248*tt276*tt252+tt539*tt14*tt&
&275*tt252+tt14*tt276*tt249*tt542*tt257-tt14*tt275*tt77*tt542*tt25&
&7)*tt258+tt599
hes(4,1) = 6*tt1*tt2*(tt254*tt337*tt252+tt77*(tt14*tt348-(tt235*t&
&t336*tt108)/2.0E+0)*tt252+tt354-tt14*tt284*tt108*tt252-tt77*tt337&
&*tt256*tt257+tt14*tt284*tt249*tt256*tt257)*tt258+tt339
hes(4,2) = 6*tt1*tt2*(tt271*tt337*tt252+tt77*(tt14*tt547-(tt235*t&
&t336*tt267)/2.0E+0)*tt252+tt548-tt14*tt284*tt267*tt252-tt77*tt337&
&*tt274*tt257+tt14*tt284*tt249*tt274*tt257)*tt258+tt546
hes(4,3) = 6*tt1*tt2*(tt276*tt337*tt252+tt579+(tt14*tt578-(tt235*&
&tt336*tt275)/2.0E+0)*tt77*tt252-tt14*tt275*tt284*tt252-tt77*tt279&
&*tt337*tt257+tt14*tt284*tt249*tt279*tt257)*tt258+tt577
hes(4,4) = 6*tt1*tt2*(tt284*tt337*tt252+tt77*(tt618+(3.0E+0*tt235&
&**2*tt617*tt249)/4.0E+0-tt235*tt336*tt335+tt14*(((tt613+tt608+tt3&
&4)*X(2,3)+tt474+tt616)*X(3,4)+((tt607+tt606+tt172)*X(2,4)+tt161+(&
&tt152+tt613+tt608)*X(2,1))*X(3,3)+(tt245*X(2,4)+tt612+tt209*X(2,1&
&))*X(3,2)+(tt611+(tt159+tt607+tt606)*X(2,3)+tt563)*X(3,1)))*tt252&
&+(tt235*tt336*tt284*tt249*tt252)/2.0E+0-tt14*(tt605+tt603+tt602-t&
&t333*X(2,4)+tt392+tt601)*tt249*tt252-tt14*tt284*tt335*tt252-tt77*&
&tt337*tt341*tt257+tt14*tt284*tt249*tt341*tt257)*tt258+6*tt1*tt2*t&
&t338**2
hes(4,5) = 6*tt1*tt2*(tt361*tt337*tt252+tt623+(tt358*tt336*tt284*&
&tt249*tt252)/2.0E+0+tt620-tt14*tt284*tt387*tt252-tt77*tt337*tt391&
&*tt257+tt14*tt284*tt249*tt391*tt257)*tt258+tt619
hes(4,6) = 6*tt1*tt2*(tt403*tt337*tt252+tt626+(tt410*tt336*tt284*&
&tt249*tt252)/2.0E+0+tt625-tt14*tt284*tt407*tt252-tt77*tt337*tt414&
&*tt257+tt14*tt284*tt249*tt414*tt257)*tt258+tt624
hes(4,7) = 6*tt1*tt2*(tt423*tt337*tt252+tt635+(tt139*tt336*tt284*&
&tt249*tt252)/2.0E+0+tt633-tt14*tt284*tt453*tt252-tt77*tt337*tt457&
&*tt257+tt14*tt284*tt249*tt457*tt257)*tt258+tt627
hes(4,8) = 6*tt1*tt2*(tt467*tt337*tt252+tt640+(tt106*tt336*tt284*&
&tt249*tt252)/2.0E+0+tt637-tt14*tt284*tt484*tt252-tt77*tt337*tt488&
&*tt257+tt14*tt284*tt249*tt488*tt257)*tt258+tt636
hes(4,9) = 6*tt1*tt2*(tt497*tt337*tt252+tt643+(tt501*tt336*tt284*&
&tt249*tt252)/2.0E+0+tt642-tt14*tt284*tt499*tt252-tt77*tt337*tt505&
&*tt257+tt14*tt284*tt249*tt505*tt257)*tt258+tt641
hes(4,10) = 6*tt1*tt2*(tt523*tt337*tt252+tt649+(tt645*tt14-(tt235&
&*tt336*tt521)/2.0E+0)*tt77*tt252-tt14*tt521*tt284*tt252-tt77*tt52&
&6*tt337*tt257+tt14*tt284*tt249*tt526*tt257)*tt258+tt644
hes(4,11) = 6*tt1*tt2*(tt532*tt337*tt252+tt654+(tt14*tt653-(tt235&
&*tt336*tt531)/2.0E+0)*tt77*tt252-tt14*tt531*tt284*tt252-tt77*tt53&
&5*tt337*tt257+tt14*tt284*tt249*tt535*tt257)*tt258+tt650
hes(4,12) = 6*tt1*tt2*(tt539*tt337*tt252+tt656+(tt14*tt334-(tt235&
&*tt336*tt248)/2.0E+0)*tt77*tt252-tt14*tt248*tt284*tt252-tt77*tt54&
&2*tt337*tt257+tt14*tt284*tt249*tt542*tt257)*tt258+tt655
hes(5,1) = 6*tt1*tt2*(tt254*tt388*tt252+tt77*(tt14*tt401-(tt358*t&
&t336*tt108)/2.0E+0)*tt252+tt402-tt14*tt361*tt108*tt252-tt77*tt388&
&*tt256*tt257+tt14*tt361*tt249*tt256*tt257)*tt258+tt390
hes(5,2) = 6*tt1*tt2*(tt271*tt388*tt252+tt77*(tt14*tt551-(tt358*t&
&t336*tt267)/2.0E+0)*tt252+tt554-tt14*tt361*tt267*tt252-tt77*tt388&
&*tt274*tt257+tt14*tt361*tt249*tt274*tt257)*tt258+tt549
hes(5,3) = 6*tt1*tt2*(tt276*tt388*tt252+tt582+(tt14*tt581-(tt358*&
&tt336*tt275)/2.0E+0)*tt77*tt252-tt14*tt275*tt361*tt252-tt77*tt279&
&*tt388*tt257+tt14*tt361*tt249*tt279*tt257)*tt258+tt580
hes(5,4) = 6*tt1*tt2*(tt284*tt388*tt252+tt623+(tt235*tt336*tt361*&
&tt249*tt252)/2.0E+0+tt620-tt14*tt361*tt335*tt252-tt77*tt388*tt341&
&*tt257+tt14*tt361*tt249*tt341*tt257)*tt258+tt619
hes(5,5) = 6*tt1*tt2*(tt361*tt388*tt252+tt77*(tt618+(3.0E+0*tt358&
&**2*tt617*tt249)/4.0E+0-tt358*tt336*tt387+tt14*((2*tt237*X(2,3)+6&
&*tt23*X(2,2)+2*tt18*X(2,1))*X(3,4)+(tt659+2*tt193*X(2,3)+6*tt189*&
&X(2,2)+2*tt187*X(2,1))*X(3,3)+(tt658+2*tt155*X(2,3)+2*tt147*X(2,1&
&))*X(3,2)+(tt657+2*tt125*X(2,3)+6*tt116*X(2,2))*X(3,1)))*tt252+(t&
&t358*tt336*tt361*tt249*tt252)/2.0E+0-tt14*(tt605+tt603+tt602-2*tt&
&25)*tt249*tt252-tt14*tt361*tt387*tt252-tt77*tt388*tt391*tt257+tt1&
&4*tt361*tt249*tt391*tt257)*tt258+6*tt1*tt2*tt389**2
hes(5,6) = 6*tt1*tt2*(tt403*tt388*tt252+tt662+(tt410*tt336*tt361*&
&tt249*tt252)/2.0E+0+tt661-tt14*tt361*tt407*tt252-tt77*tt388*tt414&
&*tt257+tt14*tt361*tt249*tt414*tt257)*tt258+tt660
hes(5,7) = 6*tt1*tt2*(tt423*tt388*tt252+tt665+(tt139*tt336*tt361*&
&tt249*tt252)/2.0E+0+tt664-tt14*tt361*tt453*tt252-tt77*tt388*tt457&
&*tt257+tt14*tt361*tt249*tt457*tt257)*tt258+tt663
hes(5,8) = 6*tt1*tt2*(tt467*tt388*tt252+tt670+(tt106*tt336*tt361*&
&tt249*tt252)/2.0E+0+tt669-tt14*tt361*tt484*tt252-tt77*tt388*tt488&
&*tt257+tt14*tt361*tt249*tt488*tt257)*tt258+tt666
hes(5,9) = 6*tt1*tt2*(tt497*tt388*tt252+tt673+(tt501*tt336*tt361*&
&tt249*tt252)/2.0E+0+tt672-tt14*tt361*tt499*tt252-tt77*tt388*tt505&
&*tt257+tt14*tt361*tt249*tt505*tt257)*tt258+tt671
hes(5,10) = 6*tt1*tt2*(tt523*tt388*tt252+tt677+(tt14*tt676-(tt358&
&*tt336*tt521)/2.0E+0)*tt77*tt252-tt14*tt521*tt361*tt252-tt77*tt52&
&6*tt388*tt257+tt14*tt361*tt249*tt526*tt257)*tt258+tt674
hes(5,11) = 6*tt1*tt2*(tt532*tt388*tt252+tt681+(tt679*tt14-(tt358&
&*tt336*tt531)/2.0E+0)*tt77*tt252-tt14*tt531*tt361*tt252-tt77*tt53&
&5*tt388*tt257+tt14*tt361*tt249*tt535*tt257)*tt258+tt678
hes(5,12) = 6*tt1*tt2*(tt539*tt388*tt252+tt683+(tt14*tt386-(tt358&
&*tt336*tt248)/2.0E+0)*tt77*tt252-tt14*tt248*tt361*tt252-tt77*tt54&
&2*tt388*tt257+tt14*tt361*tt249*tt542*tt257)*tt258+tt682
hes(6,1) = 6*tt1*tt2*(tt254*tt411*tt252+tt77*(tt14*tt415-(tt410*t&
&t336*tt108)/2.0E+0)*tt252+tt416-tt14*tt403*tt108*tt252-tt77*tt411&
&*tt256*tt257+tt14*tt403*tt249*tt256*tt257)*tt258+tt413
hes(6,2) = 6*tt1*tt2*(tt271*tt411*tt252+tt77*(tt14*tt556-(tt410*t&
&t336*tt267)/2.0E+0)*tt252+tt557-tt14*tt403*tt267*tt252-tt77*tt411&
&*tt274*tt257+tt14*tt403*tt249*tt274*tt257)*tt258+tt555
hes(6,3) = 6*tt1*tt2*(tt276*tt411*tt252+tt585+(tt584*tt14-(tt410*&
&tt336*tt275)/2.0E+0)*tt77*tt252-tt14*tt275*tt403*tt252-tt77*tt279&
&*tt411*tt257+tt14*tt403*tt249*tt279*tt257)*tt258+tt583
hes(6,4) = 6*tt1*tt2*(tt284*tt411*tt252+tt626+(tt235*tt336*tt403*&
&tt249*tt252)/2.0E+0+tt625-tt14*tt403*tt335*tt252-tt77*tt411*tt341&
&*tt257+tt14*tt403*tt249*tt341*tt257)*tt258+tt624
hes(6,5) = 6*tt1*tt2*(tt361*tt411*tt252+tt662+(tt358*tt336*tt403*&
&tt249*tt252)/2.0E+0+tt661-tt14*tt403*tt387*tt252-tt77*tt411*tt391&
&*tt257+tt14*tt403*tt249*tt391*tt257)*tt258+tt660
hes(6,6) = 6*tt1*tt2*(tt403*tt411*tt252+tt77*(tt618+(3.0E+0*tt410&
&**2*tt617*tt249)/4.0E+0-tt410*tt336*tt407+tt14*(tt684+2*tt218*X(3&
&,3)+6*tt180*X(3,2)+2*tt178*X(3,1)))*tt252+(tt410*tt336*tt403*tt24&
&9*tt252)/2.0E+0+2*tt58*tt14*tt249*tt252-tt14*tt403*tt407*tt252-tt&
&77*tt411*tt414*tt257+tt14*tt403*tt249*tt414*tt257)*tt258+6*tt1*tt&
&2*tt412**2
hes(6,7) = 6*tt1*tt2*(tt423*tt411*tt252+tt687+(tt139*tt336*tt403*&
&tt249*tt252)/2.0E+0+tt686-tt14*tt403*tt453*tt252-tt77*tt411*tt457&
&*tt257+tt14*tt403*tt249*tt457*tt257)*tt258+tt685
hes(6,8) = 6*tt1*tt2*(tt467*tt411*tt252+tt690+(tt106*tt336*tt403*&
&tt249*tt252)/2.0E+0+tt689-tt14*tt403*tt484*tt252-tt77*tt411*tt488&
&*tt257+tt14*tt403*tt249*tt488*tt257)*tt258+tt688
hes(6,9) = 6*tt1*tt2*(tt497*tt411*tt252+tt693+(tt501*tt336*tt403*&
&tt249*tt252)/2.0E+0+tt692-tt14*tt403*tt499*tt252-tt77*tt411*tt505&
&*tt257+tt14*tt403*tt249*tt505*tt257)*tt258+tt691
hes(6,10) = 6*tt1*tt2*(tt523*tt411*tt252+tt696+(tt14*tt695-(tt410&
&*tt336*tt521)/2.0E+0)*tt77*tt252-tt14*tt521*tt403*tt252-tt77*tt52&
&6*tt411*tt257+tt14*tt403*tt249*tt526*tt257)*tt258+tt694
hes(6,11) = 6*tt1*tt2*(tt532*tt411*tt252+tt699+(tt14*tt698-(tt410&
&*tt336*tt531)/2.0E+0)*tt77*tt252-tt14*tt531*tt403*tt252-tt77*tt53&
&5*tt411*tt257+tt14*tt403*tt249*tt535*tt257)*tt258+tt697
hes(6,12) = 6*tt1*tt2*(tt539*tt411*tt252+tt701+(tt406*tt14-(tt410&
&*tt336*tt248)/2.0E+0)*tt77*tt252-tt14*tt248*tt403*tt252-tt77*tt54&
&2*tt411*tt257+tt14*tt403*tt249*tt542*tt257)*tt258+tt700
hes(7,1) = 6*tt1*tt2*(tt254*tt454*tt252+tt77*(tt14*tt459-(tt139*t&
&t336*tt108)/2.0E+0)*tt252+tt462-tt14*tt423*tt108*tt252-tt77*tt454&
&*tt256*tt257+tt14*tt423*tt249*tt256*tt257)*tt258+tt456
hes(7,2) = 6*tt1*tt2*(tt271*tt454*tt252+tt77*(tt14*tt560-(tt139*t&
&t336*tt267)/2.0E+0)*tt252+tt561-tt14*tt423*tt267*tt252-tt77*tt454&
&*tt274*tt257+tt14*tt423*tt249*tt274*tt257)*tt258+tt558
hes(7,3) = 6*tt1*tt2*(tt276*tt454*tt252+tt588+(tt14*tt587-(tt139*&
&tt336*tt275)/2.0E+0)*tt77*tt252-tt14*tt275*tt423*tt252-tt77*tt279&
&*tt454*tt257+tt14*tt423*tt249*tt279*tt257)*tt258+tt586
hes(7,4) = 6*tt1*tt2*(tt284*tt454*tt252+tt635+(tt235*tt336*tt423*&
&tt249*tt252)/2.0E+0+tt633-tt14*tt423*tt335*tt252-tt77*tt454*tt341&
&*tt257+tt14*tt423*tt249*tt341*tt257)*tt258+tt627
hes(7,5) = 6*tt1*tt2*(tt361*tt454*tt252+tt665+(tt358*tt336*tt423*&
&tt249*tt252)/2.0E+0+tt664-tt14*tt423*tt387*tt252-tt77*tt454*tt391&
&*tt257+tt14*tt423*tt249*tt391*tt257)*tt258+tt663
hes(7,6) = 6*tt1*tt2*(tt403*tt454*tt252+tt687+(tt410*tt336*tt423*&
&tt249*tt252)/2.0E+0+tt686-tt14*tt423*tt407*tt252-tt77*tt454*tt414&
&*tt257+tt14*tt423*tt249*tt414*tt257)*tt258+tt685
hes(7,7) = 6*tt1*tt2*(tt423*tt454*tt252+tt77*(tt618+(3.0E+0*tt139&
&**2*tt617*tt249)/4.0E+0-tt139*tt336*tt453+tt14*((2*tt48*X(2,3)+(t&
&t609+2*tt230)*X(2,2)+tt616)*X(3,4)+(tt659+tt161+tt210*X(2,1))*X(3&
&,3)+((tt614+2*tt168)*X(2,4)+tt612+(tt159+tt609+tt313)*X(2,1))*X(3&
&,2)+(tt611+tt550+(tt152+tt614+tt327)*X(2,2))*X(3,1)))*tt252+(tt13&
&9*tt336*tt423*tt249*tt252)/2.0E+0-tt14*(tt702+tt490+tt639-(tt63+t&
&t332)*X(2,4)+tt425+tt329)*tt249*tt252-tt14*tt423*tt453*tt252-tt77&
&*tt454*tt457*tt257+tt14*tt423*tt249*tt457*tt257)*tt258+6*tt1*tt2*&
&tt455**2
hes(7,8) = 6*tt1*tt2*(tt467*tt454*tt252+tt705+(tt106*tt336*tt423*&
&tt249*tt252)/2.0E+0+tt704-tt14*tt423*tt484*tt252-tt77*tt454*tt488&
&*tt257+tt14*tt423*tt249*tt488*tt257)*tt258+tt703
hes(7,9) = 6*tt1*tt2*(tt497*tt454*tt252+tt708+(tt501*tt336*tt423*&
&tt249*tt252)/2.0E+0+tt707-tt14*tt423*tt499*tt252-tt77*tt454*tt505&
&*tt257+tt14*tt423*tt249*tt505*tt257)*tt258+tt706
hes(7,10) = 6*tt1*tt2*(tt523*tt454*tt252+tt714+(tt710*tt14-(tt139&
&*tt336*tt521)/2.0E+0)*tt77*tt252-tt14*tt521*tt423*tt252-tt77*tt52&
&6*tt454*tt257+tt14*tt423*tt249*tt526*tt257)*tt258+tt709
hes(7,11) = 6*tt1*tt2*(tt532*tt454*tt252+tt717+(tt14*tt716-(tt139&
&*tt336*tt531)/2.0E+0)*tt77*tt252-tt14*tt531*tt423*tt252-tt77*tt53&
&5*tt454*tt257+tt14*tt423*tt249*tt535*tt257)*tt258+tt715
hes(7,12) = 6*tt1*tt2*(tt539*tt454*tt252+tt719+(tt14*tt452-(tt139&
&*tt336*tt248)/2.0E+0)*tt77*tt252-tt14*tt248*tt423*tt252-tt77*tt54&
&2*tt454*tt257+tt14*tt423*tt249*tt542*tt257)*tt258+tt718
hes(8,1) = 6*tt1*tt2*(tt254*tt485*tt252+tt77*(tt14*tt495-(tt106*t&
&t336*tt108)/2.0E+0)*tt252+tt496-tt14*tt467*tt108*tt252-tt77*tt485&
&*tt256*tt257+tt14*tt467*tt249*tt256*tt257)*tt258+tt487
hes(8,2) = 6*tt1*tt2*(tt271*tt485*tt252+tt77*(tt14*tt564-(tt106*t&
&t336*tt267)/2.0E+0)*tt252+tt566-tt14*tt467*tt267*tt252-tt77*tt485&
&*tt274*tt257+tt14*tt467*tt249*tt274*tt257)*tt258+tt562
hes(8,3) = 6*tt1*tt2*(tt276*tt485*tt252+tt591+(tt14*tt590-(tt106*&
&tt336*tt275)/2.0E+0)*tt77*tt252-tt14*tt275*tt467*tt252-tt77*tt279&
&*tt485*tt257+tt14*tt467*tt249*tt279*tt257)*tt258+tt589
hes(8,4) = 6*tt1*tt2*(tt284*tt485*tt252+tt640+(tt235*tt336*tt467*&
&tt249*tt252)/2.0E+0+tt637-tt14*tt467*tt335*tt252-tt77*tt485*tt341&
&*tt257+tt14*tt467*tt249*tt341*tt257)*tt258+tt636
hes(8,5) = 6*tt1*tt2*(tt361*tt485*tt252+tt670+(tt358*tt336*tt467*&
&tt249*tt252)/2.0E+0+tt669-tt14*tt467*tt387*tt252-tt77*tt485*tt391&
&*tt257+tt14*tt467*tt249*tt391*tt257)*tt258+tt666
hes(8,6) = 6*tt1*tt2*(tt403*tt485*tt252+tt690+(tt410*tt336*tt467*&
&tt249*tt252)/2.0E+0+tt689-tt14*tt467*tt407*tt252-tt77*tt485*tt414&
&*tt257+tt14*tt467*tt249*tt414*tt257)*tt258+tt688
hes(8,7) = 6*tt1*tt2*(tt423*tt485*tt252+tt705+(tt139*tt336*tt467*&
&tt249*tt252)/2.0E+0+tt704-tt14*tt467*tt453*tt252-tt77*tt485*tt457&
&*tt257+tt14*tt467*tt249*tt457*tt257)*tt258+tt703
hes(8,8) = 6*tt1*tt2*(tt467*tt485*tt252+tt77*(tt618+(3.0E+0*tt106&
&**2*tt617*tt249)/4.0E+0-tt106*tt336*tt484+tt14*((6*tt48*X(2,3)+2*&
&tt241)*X(3,4)+(tt659+2*tt198)*X(3,3)+(tt658+6*tt155*X(2,3)+2*tt16&
&2)*X(3,2)+(tt657+6*tt132*X(2,3)+2*tt130*X(2,2))*X(3,1)))*tt252+(t&
&t106*tt336*tt467*tt249*tt252)/2.0E+0-tt14*(tt702+tt490+tt639-2*tt&
&41)*tt249*tt252-tt14*tt467*tt484*tt252-tt77*tt485*tt488*tt257+tt1&
&4*tt467*tt249*tt488*tt257)*tt258+6*tt1*tt2*tt486**2
hes(8,9) = 6*tt1*tt2*(tt497*tt485*tt252+tt722+(tt501*tt336*tt467*&
&tt249*tt252)/2.0E+0+tt721-tt14*tt467*tt499*tt252-tt77*tt485*tt505&
&*tt257+tt14*tt467*tt249*tt505*tt257)*tt258+tt720
hes(8,10) = 6*tt1*tt2*(tt523*tt485*tt252+tt725+(tt14*tt724-(tt106&
&*tt336*tt521)/2.0E+0)*tt77*tt252-tt14*tt521*tt467*tt252-tt77*tt52&
&6*tt485*tt257+tt14*tt467*tt249*tt526*tt257)*tt258+tt723
hes(8,11) = 6*tt1*tt2*(tt532*tt485*tt252+tt729+(tt727*tt14-(tt106&
&*tt336*tt531)/2.0E+0)*tt77*tt252-tt14*tt531*tt467*tt252-tt77*tt53&
&5*tt485*tt257+tt14*tt467*tt249*tt535*tt257)*tt258+tt726
hes(8,12) = 6*tt1*tt2*(tt539*tt485*tt252+tt731+(tt14*tt483-(tt106&
&*tt336*tt248)/2.0E+0)*tt77*tt252-tt14*tt248*tt467*tt252-tt77*tt54&
&2*tt485*tt257+tt14*tt467*tt249*tt542*tt257)*tt258+tt730
hes(9,1) = 6*tt1*tt2*(tt254*tt502*tt252+tt77*(tt14*tt506-(tt501*t&
&t336*tt108)/2.0E+0)*tt252+tt507-tt14*tt497*tt108*tt252-tt77*tt502&
&*tt256*tt257+tt14*tt497*tt249*tt256*tt257)*tt258+tt504
hes(9,2) = 6*tt1*tt2*(tt271*tt502*tt252+tt77*(tt14*tt568-(tt501*t&
&t336*tt267)/2.0E+0)*tt252+tt569-tt14*tt497*tt267*tt252-tt77*tt502&
&*tt274*tt257+tt14*tt497*tt249*tt274*tt257)*tt258+tt567
hes(9,3) = 6*tt1*tt2*(tt276*tt502*tt252+tt594+(tt593*tt14-(tt501*&
&tt336*tt275)/2.0E+0)*tt77*tt252-tt14*tt275*tt497*tt252-tt77*tt279&
&*tt502*tt257+tt14*tt497*tt249*tt279*tt257)*tt258+tt592
hes(9,4) = 6*tt1*tt2*(tt284*tt502*tt252+tt643+(tt235*tt336*tt497*&
&tt249*tt252)/2.0E+0+tt642-tt14*tt497*tt335*tt252-tt77*tt502*tt341&
&*tt257+tt14*tt497*tt249*tt341*tt257)*tt258+tt641
hes(9,5) = 6*tt1*tt2*(tt361*tt502*tt252+tt673+(tt358*tt336*tt497*&
&tt249*tt252)/2.0E+0+tt672-tt14*tt497*tt387*tt252-tt77*tt502*tt391&
&*tt257+tt14*tt497*tt249*tt391*tt257)*tt258+tt671
hes(9,6) = 6*tt1*tt2*(tt403*tt502*tt252+tt693+(tt410*tt336*tt497*&
&tt249*tt252)/2.0E+0+tt692-tt14*tt497*tt407*tt252-tt77*tt502*tt414&
&*tt257+tt14*tt497*tt249*tt414*tt257)*tt258+tt691
hes(9,7) = 6*tt1*tt2*(tt423*tt502*tt252+tt708+(tt139*tt336*tt497*&
&tt249*tt252)/2.0E+0+tt707-tt14*tt497*tt453*tt252-tt77*tt502*tt457&
&*tt257+tt14*tt497*tt249*tt457*tt257)*tt258+tt706
hes(9,8) = 6*tt1*tt2*(tt467*tt502*tt252+tt722+(tt106*tt336*tt497*&
&tt249*tt252)/2.0E+0+tt721-tt14*tt497*tt484*tt252-tt77*tt502*tt488&
&*tt257+tt14*tt497*tt249*tt488*tt257)*tt258+tt720
hes(9,9) = 6*tt1*tt2*(tt497*tt502*tt252+tt77*(tt618+(3.0E+0*tt501&
&**2*tt617*tt249)/4.0E+0-tt501*tt336*tt499+tt14*(tt684+6*tt226*X(3&
&,3)+2*tt225))*tt252+(tt501*tt336*tt497*tt249*tt252)/2.0E+0+2*tt70&
&*tt14*tt249*tt252-tt14*tt497*tt499*tt252-tt77*tt502*tt505*tt257+t&
&t14*tt497*tt249*tt505*tt257)*tt258+6*tt1*tt2*tt503**2
hes(9,10) = 6*tt1*tt2*(tt523*tt502*tt252+tt734+(tt14*tt733-(tt501&
&*tt336*tt521)/2.0E+0)*tt77*tt252-tt14*tt521*tt497*tt252-tt77*tt52&
&6*tt502*tt257+tt14*tt497*tt249*tt526*tt257)*tt258+tt732
hes(9,11) = 6*tt1*tt2*(tt532*tt502*tt252+tt737+(tt14*tt736-(tt501&
&*tt336*tt531)/2.0E+0)*tt77*tt252-tt14*tt531*tt497*tt252-tt77*tt53&
&5*tt502*tt257+tt14*tt497*tt249*tt535*tt257)*tt258+tt735
hes(9,12) = 6*tt1*tt2*(tt539*tt502*tt252+tt739+(tt498*tt14-(tt501&
&*tt336*tt248)/2.0E+0)*tt77*tt252-tt14*tt248*tt497*tt252-tt77*tt54&
&2*tt502*tt257+tt14*tt497*tt249*tt542*tt257)*tt258+tt738
hes(10,1) = 6*tt1*tt2*(tt527-tt14*tt523*tt108*tt252+tt14*tt521*tt&
&254*tt252+tt14*tt523*tt249*tt256*tt257-tt14*tt521*tt77*tt256*tt25&
&7)*tt258+tt525
hes(10,2) = 6*tt1*tt2*(tt538-tt14*tt523*tt267*tt252+tt571+tt14*tt&
&521*tt271*tt252+tt14*tt523*tt249*tt274*tt257-tt14*tt521*tt77*tt27&
&4*tt257)*tt258+tt570
hes(10,3) = 6*tt1*tt2*(tt544+tt596+tt14*tt521*tt276*tt252-tt14*tt&
&523*tt275*tt252+tt14*tt523*tt249*tt279*tt257-tt14*tt521*tt77*tt27&
&9*tt257)*tt258+tt595
hes(10,4) = 6*tt1*tt2*((tt235*tt336*tt523*tt249*tt252)/2.0E+0+tt6&
&49-tt14*tt523*tt335*tt252+tt645*tt14*tt77*tt252-(tt235*tt336*tt52&
&1*tt77*tt252)/2.0E+0+tt14*tt521*tt284*tt252+tt14*tt523*tt249*tt34&
&1*tt257-tt14*tt521*tt77*tt341*tt257)*tt258+tt644
hes(10,5) = 6*tt1*tt2*(tt677+(tt358*tt336*tt523*tt249*tt252)/2.0E&
&+0-tt14*tt523*tt387*tt252-(tt358*tt336*tt521*tt77*tt252)/2.0E+0+t&
&t14*tt676*tt77*tt252+tt14*tt521*tt361*tt252+tt14*tt523*tt249*tt39&
&1*tt257-tt14*tt521*tt77*tt391*tt257)*tt258+tt674
hes(10,6) = 6*tt1*tt2*(tt696+(tt410*tt336*tt523*tt249*tt252)/2.0E&
&+0-(tt410*tt336*tt521*tt77*tt252)/2.0E+0+tt14*tt695*tt77*tt252-tt&
&14*tt523*tt407*tt252+tt14*tt521*tt403*tt252+tt14*tt523*tt249*tt41&
&4*tt257-tt14*tt521*tt77*tt414*tt257)*tt258+tt694
hes(10,7) = 6*tt1*tt2*(tt714+(tt139*tt336*tt523*tt249*tt252)/2.0E&
&+0-tt14*tt523*tt453*tt252+tt710*tt14*tt77*tt252-(tt139*tt336*tt52&
&1*tt77*tt252)/2.0E+0+tt14*tt521*tt423*tt252+tt14*tt523*tt249*tt45&
&7*tt257-tt14*tt521*tt77*tt457*tt257)*tt258+tt709
hes(10,8) = 6*tt1*tt2*(tt725+(tt106*tt336*tt523*tt249*tt252)/2.0E&
&+0-tt14*tt523*tt484*tt252-(tt106*tt336*tt521*tt77*tt252)/2.0E+0+t&
&t14*tt724*tt77*tt252+tt14*tt521*tt467*tt252+tt14*tt523*tt249*tt48&
&8*tt257-tt14*tt521*tt77*tt488*tt257)*tt258+tt723
hes(10,9) = 6*tt1*tt2*(tt734+(tt501*tt336*tt523*tt249*tt252)/2.0E&
&+0-(tt501*tt336*tt521*tt77*tt252)/2.0E+0+tt14*tt733*tt77*tt252-tt&
&14*tt523*tt499*tt252+tt14*tt521*tt497*tt252+tt14*tt523*tt249*tt50&
&5*tt257-tt14*tt521*tt77*tt505*tt257)*tt258+tt732
hes(10,10) = 6*tt1*tt2*(tt14*tt523*tt249*tt526*tt257-tt14*tt521*t&
&t77*tt526*tt257)*tt258+6*tt1*tt2*tt524**2
hes(10,11) = 6*tt1*tt2*(tt14*tt532*tt521*tt252-tt14*tt523*tt531*t&
&t252+tt14*tt523*tt249*tt535*tt257-tt14*tt521*tt77*tt535*tt257)*tt&
&258+tt740
hes(10,12) = 6*tt1*tt2*(tt539*tt14*tt521*tt252-tt14*tt523*tt248*t&
&t252+tt14*tt523*tt249*tt542*tt257-tt14*tt521*tt77*tt542*tt257)*tt&
&258+tt741
hes(11,1) = 6*tt1*tt2*(tt538-tt14*tt532*tt108*tt252+tt537+tt14*tt&
&531*tt254*tt252+tt14*tt532*tt249*tt256*tt257-tt14*tt531*tt77*tt25&
&6*tt257)*tt258+tt534
hes(11,2) = 6*tt1*tt2*(tt573-tt14*tt532*tt267*tt252+tt14*tt531*tt&
&271*tt252+tt14*tt532*tt249*tt274*tt257-tt14*tt531*tt77*tt274*tt25&
&7)*tt258+tt572
hes(11,3) = 6*tt1*tt2*(tt576+tt598+tt14*tt531*tt276*tt252-tt14*tt&
&532*tt275*tt252+tt14*tt532*tt249*tt279*tt257-tt14*tt531*tt77*tt27&
&9*tt257)*tt258+tt597
hes(11,4) = 6*tt1*tt2*(tt654+(tt235*tt336*tt532*tt249*tt252)/2.0E&
&+0-tt14*tt532*tt335*tt252-(tt235*tt336*tt531*tt77*tt252)/2.0E+0+t&
&t14*tt653*tt77*tt252+tt14*tt531*tt284*tt252+tt14*tt532*tt249*tt34&
&1*tt257-tt14*tt531*tt77*tt341*tt257)*tt258+tt650
hes(11,5) = 6*tt1*tt2*((tt358*tt336*tt532*tt249*tt252)/2.0E+0+tt6&
&81-tt14*tt532*tt387*tt252+tt679*tt14*tt77*tt252-(tt358*tt336*tt53&
&1*tt77*tt252)/2.0E+0+tt14*tt531*tt361*tt252+tt14*tt532*tt249*tt39&
&1*tt257-tt14*tt531*tt77*tt391*tt257)*tt258+tt678
hes(11,6) = 6*tt1*tt2*(tt699+(tt410*tt336*tt532*tt249*tt252)/2.0E&
&+0-(tt410*tt336*tt531*tt77*tt252)/2.0E+0+tt14*tt698*tt77*tt252-tt&
&14*tt532*tt407*tt252+tt14*tt531*tt403*tt252+tt14*tt532*tt249*tt41&
&4*tt257-tt14*tt531*tt77*tt414*tt257)*tt258+tt697
hes(11,7) = 6*tt1*tt2*(tt717+(tt139*tt336*tt532*tt249*tt252)/2.0E&
&+0-tt14*tt532*tt453*tt252-(tt139*tt336*tt531*tt77*tt252)/2.0E+0+t&
&t14*tt716*tt77*tt252+tt14*tt531*tt423*tt252+tt14*tt532*tt249*tt45&
&7*tt257-tt14*tt531*tt77*tt457*tt257)*tt258+tt715
hes(11,8) = 6*tt1*tt2*(tt729+(tt106*tt336*tt532*tt249*tt252)/2.0E&
&+0-tt14*tt532*tt484*tt252+tt727*tt14*tt77*tt252-(tt106*tt336*tt53&
&1*tt77*tt252)/2.0E+0+tt14*tt531*tt467*tt252+tt14*tt532*tt249*tt48&
&8*tt257-tt14*tt531*tt77*tt488*tt257)*tt258+tt726
hes(11,9) = 6*tt1*tt2*(tt737+(tt501*tt336*tt532*tt249*tt252)/2.0E&
&+0-(tt501*tt336*tt531*tt77*tt252)/2.0E+0+tt14*tt736*tt77*tt252-tt&
&14*tt532*tt499*tt252+tt14*tt531*tt497*tt252+tt14*tt532*tt249*tt50&
&5*tt257-tt14*tt531*tt77*tt505*tt257)*tt258+tt735
hes(11,10) = 6*tt1*tt2*((-tt14*tt532*tt521*tt252)+tt14*tt523*tt53&
&1*tt252+tt14*tt532*tt249*tt526*tt257-tt14*tt531*tt77*tt526*tt257)&
&*tt258+tt740
hes(11,11) = 6*tt1*tt2*(tt14*tt532*tt249*tt535*tt257-tt14*tt531*t&
&t77*tt535*tt257)*tt258+6*tt1*tt2*tt533**2
hes(11,12) = 6*tt1*tt2*(tt539*tt14*tt531*tt252-tt14*tt532*tt248*t&
&t252+tt14*tt532*tt249*tt542*tt257-tt14*tt531*tt77*tt542*tt257)*tt&
&258+tt742
hes(12,1) = 6*tt1*tt2*(tt544-tt539*tt14*tt108*tt252+tt543+tt14*tt&
&248*tt254*tt252+tt539*tt14*tt249*tt256*tt257-tt14*tt248*tt77*tt25&
&6*tt257)*tt258+tt541
hes(12,2) = 6*tt1*tt2*(tt576-tt539*tt14*tt267*tt252+tt575+tt14*tt&
&248*tt271*tt252+tt539*tt14*tt249*tt274*tt257-tt14*tt248*tt77*tt27&
&4*tt257)*tt258+tt574
hes(12,3) = 6*tt1*tt2*(tt600+tt14*tt248*tt276*tt252-tt539*tt14*tt&
&275*tt252+tt539*tt14*tt249*tt279*tt257-tt14*tt248*tt77*tt279*tt25&
&7)*tt258+tt599
hes(12,4) = 6*tt1*tt2*(tt656+(tt235*tt539*tt336*tt249*tt252)/2.0E&
&+0-tt539*tt14*tt335*tt252-(tt235*tt336*tt248*tt77*tt252)/2.0E+0+t&
&t14*tt334*tt77*tt252+tt14*tt248*tt284*tt252+tt539*tt14*tt249*tt34&
&1*tt257-tt14*tt248*tt77*tt341*tt257)*tt258+tt655
hes(12,5) = 6*tt1*tt2*(tt683+(tt358*tt539*tt336*tt249*tt252)/2.0E&
&+0-tt539*tt14*tt387*tt252-(tt358*tt336*tt248*tt77*tt252)/2.0E+0+t&
&t14*tt386*tt77*tt252+tt14*tt248*tt361*tt252+tt539*tt14*tt249*tt39&
&1*tt257-tt14*tt248*tt77*tt391*tt257)*tt258+tt682
hes(12,6) = 6*tt1*tt2*(tt701+(tt410*tt539*tt336*tt249*tt252)/2.0E&
&+0+tt406*tt14*tt77*tt252-(tt410*tt336*tt248*tt77*tt252)/2.0E+0-tt&
&539*tt14*tt407*tt252+tt14*tt248*tt403*tt252+tt539*tt14*tt249*tt41&
&4*tt257-tt14*tt248*tt77*tt414*tt257)*tt258+tt700
hes(12,7) = 6*tt1*tt2*(tt719+(tt139*tt539*tt336*tt249*tt252)/2.0E&
&+0-tt539*tt14*tt453*tt252-(tt139*tt336*tt248*tt77*tt252)/2.0E+0+t&
&t14*tt452*tt77*tt252+tt14*tt248*tt423*tt252+tt539*tt14*tt249*tt45&
&7*tt257-tt14*tt248*tt77*tt457*tt257)*tt258+tt718
hes(12,8) = 6*tt1*tt2*(tt731+(tt106*tt539*tt336*tt249*tt252)/2.0E&
&+0-tt539*tt14*tt484*tt252-(tt106*tt336*tt248*tt77*tt252)/2.0E+0+t&
&t14*tt483*tt77*tt252+tt14*tt248*tt467*tt252+tt539*tt14*tt249*tt48&
&8*tt257-tt14*tt248*tt77*tt488*tt257)*tt258+tt730
hes(12,9) = 6*tt1*tt2*(tt739+(tt501*tt539*tt336*tt249*tt252)/2.0E&
&+0+tt498*tt14*tt77*tt252-(tt501*tt336*tt248*tt77*tt252)/2.0E+0-tt&
&539*tt14*tt499*tt252+tt14*tt248*tt497*tt252+tt539*tt14*tt249*tt50&
&5*tt257-tt14*tt248*tt77*tt505*tt257)*tt258+tt738
hes(12,10) = 6*tt1*tt2*((-tt539*tt14*tt521*tt252)+tt14*tt523*tt24&
&8*tt252+tt539*tt14*tt249*tt526*tt257-tt14*tt248*tt77*tt526*tt257)&
&*tt258+tt741
hes(12,11) = 6*tt1*tt2*((-tt539*tt14*tt531*tt252)+tt14*tt532*tt24&
&8*tt252+tt539*tt14*tt249*tt535*tt257-tt14*tt248*tt77*tt535*tt257)&
&*tt258+tt742
hes(12,12) = 6*tt1*tt2*(tt539*tt14*tt249*tt542*tt257-tt14*tt248*t&
&t77*tt542*tt257)*tt258+6*tt1*tt2*tt540**2
END 
SUBROUTINE line_bending(val, X, d1, d2) 
IMPLICIT NONE 
REAL(KIND=8) val(1, 1) 
REAL(KIND=8) X(3, 4) 
REAL(KIND=8) d1(1, 1) 
REAL(KIND=8) d2(1, 1) 
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
tt1 = X(1,3)**2
tt2 = 4*tt1
tt3 = -8*X(1,3)*X(1,4)
tt4 = X(1,4)**2
tt5 = 4*tt4
tt6 = tt5+tt3+tt2
tt7 = X(2,1)**2
tt8 = -8*tt1
tt9 = 16*X(1,3)*X(1,4)
tt10 = -8*tt4
tt11 = X(2,2)**2
tt12 = 8*X(1,2)-8*X(1,1)
tt13 = tt12*X(1,3)
tt14 = 8*X(1,1)-8*X(1,2)
tt15 = tt14*X(1,4)
tt16 = tt15+tt13
tt17 = tt14*X(1,3)
tt18 = tt12*X(1,4)
tt19 = tt18+tt17
tt20 = X(1,1)**2
tt21 = 4*tt20
tt22 = -8*X(1,1)*X(1,2)
tt23 = X(1,2)**2
tt24 = 4*tt23
tt25 = tt24+tt22+tt21
tt26 = X(2,3)**2
tt27 = -8*tt20
tt28 = 16*X(1,1)*X(1,2)
tt29 = -8*tt23
tt30 = X(2,4)**2
tt31 = 4*tt30-8*X(2,3)*X(2,4)+4*tt26+tt5+tt3+tt2
tt32 = X(3,1)**2
tt33 = X(3,2)**2
tt34 = 8*X(2,2)-8*X(2,1)
tt35 = 8*X(2,1)-8*X(2,2)
tt36 = tt35*X(2,4)+tt34*X(2,3)+tt15+tt13
tt37 = tt34*X(2,4)+tt35*X(2,3)+tt18+tt17
tt38 = 4*tt11-8*X(2,1)*X(2,2)+4*tt7+tt24+tt22+tt21
tt39 = X(3,3)**2
tt40 = X(3,4)**2
tt41 = tt23-2*X(1,1)*X(1,2)+tt20
tt42 = 2*d1(1,1)*d2(1,1)
tt43 = 2*X(1,1)-2*X(1,2)
tt44 = tt43*X(1,3)
tt45 = 2*X(1,2)-2*X(1,1)
tt46 = tt45*X(1,4)
tt47 = tt46+tt44+tt42
tt48 = -2*d1(1,1)*d2(1,1)
tt49 = tt45*X(1,3)
tt50 = tt43*X(1,4)
tt51 = tt50+tt49+tt48
tt52 = tt11-2*X(2,1)*X(2,2)+tt7
tt53 = 2*X(2,1)-2*X(2,2)
tt54 = 2*X(2,2)-2*X(2,1)
tt55 = tt54*X(2,4)+tt53*X(2,3)+tt46+tt44+tt42
tt56 = tt53*X(2,4)+tt54*X(2,3)+tt50+tt49+tt48
tt57 = tt33-2*X(3,1)*X(3,2)+tt32
val(1,1) = (tt38*tt40+(((-8*tt11)+16*X(2,1)*X(2,2)-8*tt7+tt29+tt2&
&8+tt27)*X(3,3)+tt36*X(3,2)+tt37*X(3,1))*X(3,4)+tt38*tt39+(tt37*X(&
&3,2)+tt36*X(3,1))*X(3,3)+tt31*tt33+((-8*tt30)+16*X(2,3)*X(2,4)-8*&
&tt26+tt10+tt9+tt8)*X(3,1)*X(3,2)+tt31*tt32+tt25*tt30+((tt29+tt28+&
&tt27)*X(2,3)+tt16*X(2,2)+tt19*X(2,1))*X(2,4)+tt25*tt26+(tt19*X(2,&
&2)+tt16*X(2,1))*X(2,3)+tt6*tt11+(tt10+tt9+tt8)*X(2,1)*X(2,2)+tt6*&
&tt7)/((d2(1,1)+d1(1,1))*(tt57*tt40+(((-2*tt33)+4*X(3,1)*X(3,2)-2*&
&tt32)*X(3,3)+tt55*X(3,2)+tt56*X(3,1))*X(3,4)+tt57*tt39+(tt56*X(3,&
&2)+tt55*X(3,1))*X(3,3)+tt52*tt30+(((-2*tt11)+4*X(2,1)*X(2,2)-2*tt&
&7)*X(2,3)+tt47*X(2,2)+tt51*X(2,1))*X(2,4)+tt52*tt26+(tt51*X(2,2)+&
&tt47*X(2,1))*X(2,3)+tt41*tt4+(((-2*tt23)+4*X(1,1)*X(1,2)-2*tt20)*&
&X(1,3)+2*d1(1,1)*d2(1,1)*X(1,2)-2*X(1,1)*d1(1,1)*d2(1,1))*X(1,4)+&
&tt41*tt1+(2*X(1,1)*d1(1,1)*d2(1,1)-2*d1(1,1)*d2(1,1)*X(1,2))*X(1,&
&3)+d1(1,1)**2*d2(1,1)**2))
END 
SUBROUTINE line_bending_jac(jac, X, d1, d2) 
IMPLICIT NONE 
REAL(KIND=8) jac(1, 12) 
REAL(KIND=8) X(3, 4) 
REAL(KIND=8) d1(1, 1) 
REAL(KIND=8) d2(1, 1) 
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
tt1 = 1/(d2(1,1)+d1(1,1))
tt2 = 2*X(1,1)-2*X(1,2)
tt3 = X(1,3)**2
tt4 = -2*d1(1,1)*d2(1,1)
tt5 = X(1,4)**2
tt6 = 2*X(1,3)-2*X(1,4)
tt7 = 2*X(1,4)-2*X(1,3)
tt8 = tt7*X(2,2)+tt6*X(2,1)
tt9 = tt6*X(2,2)+tt7*X(2,1)
tt10 = tt7*X(3,2)+tt6*X(3,1)
tt11 = tt6*X(3,2)+tt7*X(3,1)
tt12 = 4*tt3
tt13 = -8*X(1,3)*X(1,4)
tt14 = 4*tt5
tt15 = tt14+tt13+tt12
tt16 = X(2,1)**2
tt17 = -8*tt3
tt18 = 16*X(1,3)*X(1,4)
tt19 = -8*tt5
tt20 = tt19+tt18+tt17
tt21 = X(2,2)**2
tt22 = 8*X(1,2)-8*X(1,1)
tt23 = tt22*X(1,3)
tt24 = 8*X(1,1)-8*X(1,2)
tt25 = tt24*X(1,4)
tt26 = tt25+tt23
tt27 = tt26*X(2,1)
tt28 = tt24*X(1,3)
tt29 = tt22*X(1,4)
tt30 = tt29+tt28
tt31 = tt30*X(2,2)
tt32 = X(1,1)**2
tt33 = 4*tt32
tt34 = -8*X(1,1)*X(1,2)
tt35 = X(1,2)**2
tt36 = 4*tt35
tt37 = tt36+tt34+tt33
tt38 = X(2,3)**2
tt39 = tt30*X(2,1)
tt40 = tt26*X(2,2)
tt41 = -8*tt32
tt42 = 16*X(1,1)*X(1,2)
tt43 = -8*tt35
tt44 = tt43+tt42+tt41
tt45 = tt44*X(2,3)
tt46 = X(2,4)**2
tt47 = 4*tt46-8*X(2,3)*X(2,4)+4*tt38+tt14+tt13+tt12
tt48 = X(3,1)**2
tt49 = (-8*tt46)+16*X(2,3)*X(2,4)-8*tt38+tt19+tt18+tt17
tt50 = X(3,2)**2
tt51 = 8*X(2,2)-8*X(2,1)
tt52 = 8*X(2,1)-8*X(2,2)
tt53 = tt52*X(2,4)+tt51*X(2,3)+tt25+tt23
tt54 = tt53*X(3,1)
tt55 = tt51*X(2,4)+tt52*X(2,3)+tt29+tt28
tt56 = tt55*X(3,2)
tt57 = 4*tt21-8*X(2,1)*X(2,2)+4*tt16+tt36+tt34+tt33
tt58 = X(3,3)**2
tt59 = tt55*X(3,1)
tt60 = tt53*X(3,2)
tt61 = (-8*tt21)+16*X(2,1)*X(2,2)-8*tt16+tt43+tt42+tt41
tt62 = tt61*X(3,3)
tt63 = X(3,4)**2
tt64 = tt57*tt63+(tt62+tt60+tt59)*X(3,4)+tt57*tt58+(tt56+tt54)*X(&
&3,3)+tt47*tt50+tt49*X(3,1)*X(3,2)+tt47*tt48+tt37*tt46+(tt45+tt40+&
&tt39)*X(2,4)+tt37*tt38+(tt31+tt27)*X(2,3)+tt15*tt21+tt20*X(2,1)*X&
&(2,2)+tt15*tt16
tt65 = 2*X(1,1)*d1(1,1)*d2(1,1)
tt66 = -2*d1(1,1)*d2(1,1)*X(1,2)
tt67 = tt35-2*X(1,1)*X(1,2)+tt32
tt68 = -2*X(1,1)*d1(1,1)*d2(1,1)
tt69 = 2*d1(1,1)*d2(1,1)*X(1,2)
tt70 = (-2*tt35)+4*X(1,1)*X(1,2)-2*tt32
tt71 = tt70*X(1,3)
tt72 = 2*d1(1,1)*d2(1,1)
tt73 = tt2*X(1,3)
tt74 = 2*X(1,2)-2*X(1,1)
tt75 = tt74*X(1,4)
tt76 = tt75+tt73+tt72
tt77 = tt76*X(2,1)
tt78 = tt74*X(1,3)
tt79 = tt2*X(1,4)
tt80 = tt79+tt78+tt4
tt81 = tt80*X(2,2)
tt82 = tt21-2*X(2,1)*X(2,2)+tt16
tt83 = tt80*X(2,1)
tt84 = tt76*X(2,2)
tt85 = (-2*tt21)+4*X(2,1)*X(2,2)-2*tt16
tt86 = tt85*X(2,3)
tt87 = 2*X(2,1)-2*X(2,2)
tt88 = tt87*X(2,3)
tt89 = 2*X(2,2)-2*X(2,1)
tt90 = tt89*X(2,4)
tt91 = tt90+tt88+tt75+tt73+tt72
tt92 = tt91*X(3,1)
tt93 = tt89*X(2,3)
tt94 = tt87*X(2,4)
tt95 = tt94+tt93+tt79+tt78+tt4
tt96 = tt95*X(3,2)
tt97 = tt50-2*X(3,1)*X(3,2)+tt48
tt98 = tt95*X(3,1)
tt99 = tt91*X(3,2)
tt100 = (-2*tt50)+4*X(3,1)*X(3,2)-2*tt48
tt101 = tt100*X(3,3)
tt102 = tt97*tt63+(tt101+tt99+tt98)*X(3,4)+tt97*tt58+(tt96+tt92)*&
&X(3,3)+tt82*tt46+(tt86+tt84+tt83)*X(2,4)+tt82*tt38+(tt81+tt77)*X(&
&2,3)+tt67*tt5+(tt71+tt69+tt68)*X(1,4)+tt67*tt3+(tt66+tt65)*X(1,3)&
&+d1(1,1)**2*d2(1,1)**2
tt103 = 1/tt102**2
tt104 = 8*X(1,4)-8*X(1,3)
tt105 = tt104*X(2,1)
tt106 = 8*X(1,3)-8*X(1,4)
tt107 = tt106*X(2,2)
tt108 = tt106*X(2,1)
tt109 = tt104*X(2,2)
tt110 = 16*X(1,2)-16*X(1,1)
tt111 = tt104*X(3,1)
tt112 = tt106*X(3,2)
tt113 = tt106*X(3,1)
tt114 = tt104*X(3,2)
tt115 = 1/tt102
tt116 = 2*X(2,3)-2*X(2,4)
tt117 = 2*X(2,4)-2*X(2,3)
tt118 = tt117*X(3,2)+tt116*X(3,1)
tt119 = tt116*X(3,2)+tt117*X(3,1)
tt120 = 8*X(2,4)-8*X(2,3)
tt121 = tt120*X(3,1)
tt122 = 8*X(2,3)-8*X(2,4)
tt123 = tt122*X(3,2)
tt124 = tt122*X(3,1)
tt125 = tt120*X(3,2)
tt126 = 2*X(3,1)-2*X(3,2)
tt127 = 16*X(1,1)-16*X(1,2)
tt128 = 2*X(3,2)-2*X(3,1)
tt129 = tt74*X(2,2)+tt2*X(2,1)
tt130 = tt2*X(2,2)+tt74*X(2,1)
tt131 = tt74*X(3,2)+tt2*X(3,1)
tt132 = tt2*X(3,2)+tt74*X(3,1)
tt133 = 16*X(1,4)-16*X(1,3)
tt134 = tt24*X(2,2)+tt22*X(2,1)
tt135 = tt22*X(2,2)+tt24*X(2,1)
tt136 = tt24*X(3,2)+tt22*X(3,1)
tt137 = tt22*X(3,2)+tt24*X(3,1)
tt138 = tt89*X(3,2)+tt87*X(3,1)
tt139 = tt87*X(3,2)+tt89*X(3,1)
tt140 = tt52*X(3,2)+tt51*X(3,1)
tt141 = tt51*X(3,2)+tt52*X(3,1)
tt142 = 16*X(1,3)-16*X(1,4)
jac(1,1) = tt1*(tt24*tt63+(tt110*X(3,3)+tt114+tt113)*X(3,4)+tt24*&
&tt58+(tt112+tt111)*X(3,3)+tt24*tt46+(tt110*X(2,3)+tt109+tt108)*X(&
&2,4)+tt24*tt38+(tt107+tt105)*X(2,3))*tt115-tt1*(tt11*X(3,4)+tt10*&
&X(3,3)+tt9*X(2,4)+tt8*X(2,3)+tt2*tt5+((4*X(1,2)-4*X(1,1))*X(1,3)+&
&tt4)*X(1,4)+tt2*tt3+2*d1(1,1)*d2(1,1)*X(1,3))*tt64*tt103
jac(1,2) = tt1*(tt52*tt63+((16*X(2,2)-16*X(2,1))*X(3,3)+tt125+tt1&
&24)*X(3,4)+tt52*tt58+(tt123+tt121)*X(3,3)+tt30*X(2,4)+tt26*X(2,3)&
&+tt20*X(2,2)+2*tt15*X(2,1))*tt115-tt1*(tt119*X(3,4)+tt118*X(3,3)+&
&tt87*tt46+((4*X(2,2)-4*X(2,1))*X(2,3)+tt79+tt78+tt4)*X(2,4)+tt87*&
&tt38+tt76*X(2,3))*tt64*tt103
jac(1,3) = tt1*(tt55*X(3,4)+tt53*X(3,3)+tt49*X(3,2)+2*tt47*X(3,1)&
&)*tt115-tt1*tt64*(tt126*tt63+((4*X(3,2)-4*X(3,1))*X(3,3)+tt94+tt9&
&3+tt79+tt78+tt4)*X(3,4)+tt126*tt58+tt91*X(3,3))*tt103
jac(1,4) = tt1*(tt22*tt63+(tt127*X(3,3)+tt112+tt111)*X(3,4)+tt22*&
&tt58+(tt114+tt113)*X(3,3)+tt22*tt46+(tt127*X(2,3)+tt107+tt105)*X(&
&2,4)+tt22*tt38+(tt109+tt108)*X(2,3))*tt115-tt1*(tt10*X(3,4)+tt11*&
&X(3,3)+tt8*X(2,4)+tt9*X(2,3)+tt74*tt5+((4*X(1,1)-4*X(1,2))*X(1,3)&
&+tt72)*X(1,4)+tt74*tt3-2*d1(1,1)*d2(1,1)*X(1,3))*tt64*tt103
jac(1,5) = tt1*(tt51*tt63+((16*X(2,1)-16*X(2,2))*X(3,3)+tt123+tt1&
&21)*X(3,4)+tt51*tt58+(tt125+tt124)*X(3,3)+tt26*X(2,4)+tt30*X(2,3)&
&+2*tt15*X(2,2)+tt20*X(2,1))*tt115-tt1*(tt118*X(3,4)+tt119*X(3,3)+&
&tt89*tt46+((4*X(2,1)-4*X(2,2))*X(2,3)+tt75+tt73+tt72)*X(2,4)+tt89&
&*tt38+tt80*X(2,3))*tt64*tt103
jac(1,6) = tt1*(tt53*X(3,4)+tt55*X(3,3)+2*tt47*X(3,2)+tt49*X(3,1)&
&)*tt115-tt1*tt64*(tt128*tt63+((4*X(3,1)-4*X(3,2))*X(3,3)+tt90+tt8&
&8+tt75+tt73+tt72)*X(3,4)+tt128*tt58+tt95*X(3,3))*tt103
jac(1,7) = tt1*(tt137*X(3,4)+tt136*X(3,3)+tt106*tt50+tt133*X(3,1)&
&*X(3,2)+tt106*tt48+tt135*X(2,4)+tt134*X(2,3)+tt106*tt21+tt133*X(2&
&,1)*X(2,2)+tt106*tt16)*tt115-tt1*(tt132*X(3,4)+tt131*X(3,3)+tt130&
&*X(2,4)+tt129*X(2,3)+tt70*X(1,4)+2*tt67*X(1,3)+tt66+tt65)*tt64*tt&
&103
jac(1,8) = tt1*(tt141*X(3,4)+tt140*X(3,3)+tt122*tt50+(16*X(2,4)-1&
&6*X(2,3))*X(3,1)*X(3,2)+tt122*tt48+tt44*X(2,4)+2*tt37*X(2,3)+tt31&
&+tt27)*tt115-tt1*(tt139*X(3,4)+tt138*X(3,3)+tt85*X(2,4)+2*tt82*X(&
&2,3)+tt81+tt77)*tt64*tt103
jac(1,9) = tt1*(tt61*X(3,4)+2*tt57*X(3,3)+tt56+tt54)*tt115-tt1*(t&
&t100*X(3,4)+2*tt97*X(3,3)+tt96+tt92)*tt64*tt103
jac(1,10) = tt1*(tt136*X(3,4)+tt137*X(3,3)+tt104*tt50+tt142*X(3,1&
&)*X(3,2)+tt104*tt48+tt134*X(2,4)+tt135*X(2,3)+tt104*tt21+tt142*X(&
&2,1)*X(2,2)+tt104*tt16)*tt115-tt1*(tt131*X(3,4)+tt132*X(3,3)+tt12&
&9*X(2,4)+tt130*X(2,3)+2*tt67*X(1,4)+tt71+tt69+tt68)*tt64*tt103
jac(1,11) = tt1*(tt140*X(3,4)+tt141*X(3,3)+tt120*tt50+(16*X(2,3)-&
&16*X(2,4))*X(3,1)*X(3,2)+tt120*tt48+2*tt37*X(2,4)+tt45+tt40+tt39)&
&*tt115-tt1*(tt138*X(3,4)+tt139*X(3,3)+2*tt82*X(2,4)+tt86+tt84+tt8&
&3)*tt64*tt103
jac(1,12) = tt1*(2*tt57*X(3,4)+tt62+tt60+tt59)*tt115-tt1*(2*tt97*&
&X(3,4)+tt101+tt99+tt98)*tt64*tt103
END 
SUBROUTINE line_bending_hes(hes, X, d1, d2) 
IMPLICIT NONE 
REAL(KIND=8) hes(12, 12) 
REAL(KIND=8) X(3, 4) 
REAL(KIND=8) d1(1, 1) 
REAL(KIND=8) d2(1, 1) 
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
tt1 = 1/(d2(1,1)+d1(1,1))
tt2 = 2*X(1,1)-2*X(1,2)
tt3 = X(1,3)**2
tt4 = -2*d1(1,1)*d2(1,1)
tt5 = 4*X(1,2)-4*X(1,1)
tt6 = tt5*X(1,3)
tt7 = X(1,4)**2
tt8 = 2*X(1,3)-2*X(1,4)
tt9 = 2*X(1,4)-2*X(1,3)
tt10 = tt9*X(2,2)+tt8*X(2,1)
tt11 = tt8*X(2,2)+tt9*X(2,1)
tt12 = tt9*X(3,2)+tt8*X(3,1)
tt13 = tt8*X(3,2)+tt9*X(3,1)
tt14 = tt13*X(3,4)+tt12*X(3,3)+tt11*X(2,4)+tt10*X(2,3)+tt2*tt7+(t&
&t6+tt4)*X(1,4)+tt2*tt3+2*d1(1,1)*d2(1,1)*X(1,3)
tt15 = 4*tt3
tt16 = -8*X(1,3)*X(1,4)
tt17 = 4*tt7
tt18 = tt17+tt16+tt15
tt19 = X(2,1)**2
tt20 = -8*tt3
tt21 = 16*X(1,3)*X(1,4)
tt22 = -8*tt7
tt23 = tt22+tt21+tt20
tt24 = X(2,2)**2
tt25 = 8*X(1,2)-8*X(1,1)
tt26 = tt25*X(1,3)
tt27 = 8*X(1,1)-8*X(1,2)
tt28 = tt27*X(1,4)
tt29 = tt28+tt26
tt30 = tt29*X(2,1)
tt31 = tt27*X(1,3)
tt32 = tt25*X(1,4)
tt33 = tt32+tt31
tt34 = tt33*X(2,2)
tt35 = X(1,1)**2
tt36 = 4*tt35
tt37 = -8*X(1,1)*X(1,2)
tt38 = X(1,2)**2
tt39 = 4*tt38
tt40 = tt39+tt37+tt36
tt41 = X(2,3)**2
tt42 = tt33*X(2,1)
tt43 = tt29*X(2,2)
tt44 = -8*tt35
tt45 = 16*X(1,1)*X(1,2)
tt46 = -8*tt38
tt47 = tt46+tt45+tt44
tt48 = tt47*X(2,3)
tt49 = X(2,4)**2
tt50 = 4*tt49-8*X(2,3)*X(2,4)+4*tt41+tt17+tt16+tt15
tt51 = X(3,1)**2
tt52 = -8*tt41
tt53 = 16*X(2,3)*X(2,4)
tt54 = -8*tt49
tt55 = tt54+tt53+tt52+tt22+tt21+tt20
tt56 = X(3,2)**2
tt57 = 8*X(2,2)-8*X(2,1)
tt58 = tt57*X(2,3)
tt59 = 8*X(2,1)-8*X(2,2)
tt60 = tt59*X(2,4)
tt61 = tt60+tt58+tt28+tt26
tt62 = tt61*X(3,1)
tt63 = tt59*X(2,3)
tt64 = tt57*X(2,4)
tt65 = tt64+tt63+tt32+tt31
tt66 = tt65*X(3,2)
tt67 = 4*tt24-8*X(2,1)*X(2,2)+4*tt19+tt39+tt37+tt36
tt68 = X(3,3)**2
tt69 = tt65*X(3,1)
tt70 = tt61*X(3,2)
tt71 = -8*tt19
tt72 = 16*X(2,1)*X(2,2)
tt73 = -8*tt24
tt74 = tt73+tt72+tt71+tt46+tt45+tt44
tt75 = tt74*X(3,3)
tt76 = X(3,4)**2
tt77 = tt67*tt76+(tt75+tt70+tt69)*X(3,4)+tt67*tt68+(tt66+tt62)*X(&
&3,3)+tt50*tt56+tt55*X(3,1)*X(3,2)+tt50*tt51+tt40*tt49+(tt48+tt43+&
&tt42)*X(2,4)+tt40*tt41+(tt34+tt30)*X(2,3)+tt18*tt24+tt23*X(2,1)*X&
&(2,2)+tt18*tt19
tt78 = 2*X(1,1)*d1(1,1)*d2(1,1)
tt79 = -2*d1(1,1)*d2(1,1)*X(1,2)
tt80 = tt38-2*X(1,1)*X(1,2)+tt35
tt81 = -2*X(1,1)*d1(1,1)*d2(1,1)
tt82 = 2*d1(1,1)*d2(1,1)*X(1,2)
tt83 = (-2*tt38)+4*X(1,1)*X(1,2)-2*tt35
tt84 = tt83*X(1,3)
tt85 = 2*d1(1,1)*d2(1,1)
tt86 = tt2*X(1,3)
tt87 = 2*X(1,2)-2*X(1,1)
tt88 = tt87*X(1,4)
tt89 = tt88+tt86+tt85
tt90 = tt89*X(2,1)
tt91 = tt87*X(1,3)
tt92 = tt2*X(1,4)
tt93 = tt92+tt91+tt4
tt94 = tt93*X(2,2)
tt95 = tt24-2*X(2,1)*X(2,2)+tt19
tt96 = tt93*X(2,1)
tt97 = tt89*X(2,2)
tt98 = (-2*tt24)+4*X(2,1)*X(2,2)-2*tt19
tt99 = tt98*X(2,3)
tt100 = 2*X(2,1)-2*X(2,2)
tt101 = tt100*X(2,3)
tt102 = 2*X(2,2)-2*X(2,1)
tt103 = tt102*X(2,4)
tt104 = tt103+tt101+tt88+tt86+tt85
tt105 = tt104*X(3,1)
tt106 = tt102*X(2,3)
tt107 = tt100*X(2,4)
tt108 = tt107+tt106+tt92+tt91+tt4
tt109 = tt108*X(3,2)
tt110 = tt56-2*X(3,1)*X(3,2)+tt51
tt111 = tt108*X(3,1)
tt112 = tt104*X(3,2)
tt113 = (-2*tt56)+4*X(3,1)*X(3,2)-2*tt51
tt114 = tt113*X(3,3)
tt115 = tt110*tt76+(tt114+tt112+tt111)*X(3,4)+tt110*tt68+(tt109+t&
&t105)*X(3,3)+tt95*tt49+(tt99+tt97+tt96)*X(2,4)+tt95*tt41+(tt94+tt&
&90)*X(2,3)+tt80*tt7+(tt84+tt82+tt81)*X(1,4)+tt80*tt3+(tt79+tt78)*&
&X(1,3)+d1(1,1)**2*d2(1,1)**2
tt116 = 1/tt115**3
tt117 = 8*X(1,4)-8*X(1,3)
tt118 = tt117*X(2,1)
tt119 = 8*X(1,3)-8*X(1,4)
tt120 = tt119*X(2,2)
tt121 = tt119*X(2,1)
tt122 = tt117*X(2,2)
tt123 = 16*X(1,2)-16*X(1,1)
tt124 = tt123*X(2,3)
tt125 = tt117*X(3,1)
tt126 = tt119*X(3,2)
tt127 = tt119*X(3,1)
tt128 = tt117*X(3,2)
tt129 = tt123*X(3,3)
tt130 = tt27*tt76+(tt129+tt128+tt127)*X(3,4)+tt27*tt68+(tt126+tt1&
&25)*X(3,3)+tt27*tt49+(tt124+tt122+tt121)*X(2,4)+tt27*tt41+(tt120+&
&tt118)*X(2,3)
tt131 = 1/tt115**2
tt132 = -tt1*(2*tt7-4*X(1,3)*X(1,4)+2*tt3)*tt77*tt131
tt133 = 8*tt68
tt134 = -16*X(3,3)*X(3,4)
tt135 = 8*tt76
tt136 = 1/tt115
tt137 = tt1*(tt135+tt134+tt133+8*tt49-16*X(2,3)*X(2,4)+8*tt41)*tt&
&136
tt138 = 4*X(2,2)-4*X(2,1)
tt139 = tt138*X(2,3)
tt140 = 2*X(2,3)-2*X(2,4)
tt141 = 2*X(2,4)-2*X(2,3)
tt142 = tt141*X(3,2)+tt140*X(3,1)
tt143 = tt140*X(3,2)+tt141*X(3,1)
tt144 = tt143*X(3,4)+tt142*X(3,3)+tt100*tt49+(tt139+tt92+tt91+tt4&
&)*X(2,4)+tt100*tt41+tt89*X(2,3)
tt145 = 8*X(2,4)-8*X(2,3)
tt146 = tt145*X(3,1)
tt147 = 8*X(2,3)-8*X(2,4)
tt148 = tt147*X(3,2)
tt149 = tt147*X(3,1)
tt150 = tt145*X(3,2)
tt151 = 16*X(2,2)-16*X(2,1)
tt152 = tt151*X(3,3)
tt153 = tt59*tt76+(tt152+tt150+tt149)*X(3,4)+tt59*tt68+(tt148+tt1&
&46)*X(3,3)+tt33*X(2,4)+tt29*X(2,3)+tt23*X(2,2)+2*tt18*X(2,1)
tt154 = -tt1*(tt9*X(2,4)+tt8*X(2,3))*tt77*tt131
tt155 = tt1*(tt119*X(2,4)+tt117*X(2,3))*tt136
tt156 = tt155+tt154-tt1*tt14*tt153*tt131-tt1*tt144*tt130*tt131+2*&
&tt1*tt14*tt144*tt77*tt116
tt157 = 2*X(3,1)-2*X(3,2)
tt158 = 4*X(3,2)-4*X(3,1)
tt159 = tt158*X(3,3)
tt160 = tt157*tt76+(tt159+tt107+tt106+tt92+tt91+tt4)*X(3,4)+tt157&
&*tt68+tt104*X(3,3)
tt161 = tt65*X(3,4)+tt61*X(3,3)+tt55*X(3,2)+2*tt50*X(3,1)
tt162 = -tt1*(tt9*X(3,4)+tt8*X(3,3))*tt77*tt131
tt163 = tt1*(tt119*X(3,4)+tt117*X(3,3))*tt136
tt164 = tt163-tt1*tt130*tt160*tt131+tt162-tt1*tt161*tt14*tt131+2*&
&tt1*tt14*tt77*tt160*tt116
tt165 = 4*X(1,1)-4*X(1,2)
tt166 = tt165*X(1,3)
tt167 = tt12*X(3,4)+tt13*X(3,3)+tt10*X(2,4)+tt11*X(2,3)+tt87*tt7+&
&(tt166+tt85)*X(1,4)+tt87*tt3-2*d1(1,1)*d2(1,1)*X(1,3)
tt168 = 16*X(1,1)-16*X(1,2)
tt169 = tt168*X(2,3)
tt170 = tt168*X(3,3)
tt171 = tt25*tt76+(tt170+tt126+tt125)*X(3,4)+tt25*tt68+(tt128+tt1&
&27)*X(3,3)+tt25*tt49+(tt169+tt120+tt118)*X(2,4)+tt25*tt41+(tt122+&
&tt121)*X(2,3)
tt172 = -8*tt68
tt173 = 16*X(3,3)*X(3,4)
tt174 = -8*tt76
tt175 = tt1*(tt174+tt173+tt172+tt54+tt53+tt52)*tt136-tt1*((-2*tt7&
&)+4*X(1,3)*X(1,4)-2*tt3)*tt77*tt131-tt1*tt14*tt171*tt131-tt1*tt16&
&7*tt130*tt131+2*tt1*tt14*tt167*tt77*tt116
tt176 = 4*X(2,1)-4*X(2,2)
tt177 = tt176*X(2,3)
tt178 = tt142*X(3,4)+tt143*X(3,3)+tt102*tt49+(tt177+tt88+tt86+tt8&
&5)*X(2,4)+tt102*tt41+tt93*X(2,3)
tt179 = 16*X(2,1)-16*X(2,2)
tt180 = tt179*X(3,3)
tt181 = tt57*tt76+(tt180+tt148+tt146)*X(3,4)+tt57*tt68+(tt150+tt1&
&49)*X(3,3)+tt29*X(2,4)+tt33*X(2,3)+2*tt18*X(2,2)+tt23*X(2,1)
tt182 = -tt1*(tt8*X(2,4)+tt9*X(2,3))*tt77*tt131
tt183 = tt1*(tt117*X(2,4)+tt119*X(2,3))*tt136
tt184 = tt183+tt182-tt1*tt14*tt181*tt131-tt1*tt178*tt130*tt131+2*&
&tt1*tt14*tt178*tt77*tt116
tt185 = 2*X(3,2)-2*X(3,1)
tt186 = 4*X(3,1)-4*X(3,2)
tt187 = tt186*X(3,3)
tt188 = tt185*tt76+(tt187+tt103+tt101+tt88+tt86+tt85)*X(3,4)+tt18&
&5*tt68+tt108*X(3,3)
tt189 = tt61*X(3,4)+tt65*X(3,3)+2*tt50*X(3,2)+tt55*X(3,1)
tt190 = -tt1*(tt8*X(3,4)+tt9*X(3,3))*tt77*tt131
tt191 = tt1*(tt117*X(3,4)+tt119*X(3,3))*tt136
tt192 = tt191-tt1*tt130*tt188*tt131+tt190-tt1*tt189*tt14*tt131+2*&
&tt1*tt14*tt77*tt188*tt116
tt193 = tt87*X(2,2)+tt2*X(2,1)
tt194 = tt2*X(2,2)+tt87*X(2,1)
tt195 = tt87*X(3,2)+tt2*X(3,1)
tt196 = tt2*X(3,2)+tt87*X(3,1)
tt197 = tt196*X(3,4)+tt195*X(3,3)+tt194*X(2,4)+tt193*X(2,3)+tt83*&
&X(1,4)+2*tt80*X(1,3)+tt79+tt78
tt198 = 16*X(1,4)-16*X(1,3)
tt199 = tt27*X(2,2)+tt25*X(2,1)
tt200 = tt25*X(2,2)+tt27*X(2,1)
tt201 = tt27*X(3,2)+tt25*X(3,1)
tt202 = tt25*X(3,2)+tt27*X(3,1)
tt203 = tt202*X(3,4)+tt201*X(3,3)+tt119*tt56+tt198*X(3,1)*X(3,2)+&
&tt119*tt51+tt200*X(2,4)+tt199*X(2,3)+tt119*tt24+tt198*X(2,1)*X(2,&
&2)+tt119*tt19
tt204 = tt157*X(3,3)
tt205 = tt185*X(3,4)
tt206 = 8*X(3,2)-8*X(3,1)
tt207 = tt206*X(3,3)
tt208 = 8*X(3,1)-8*X(3,2)
tt209 = tt208*X(3,4)
tt210 = tt1*(tt209+tt207+tt60+tt58)*tt136
tt211 = tt210-tt1*(tt205+tt204+tt103+tt101+tt5*X(1,4)+2*tt2*X(1,3&
&)+tt85)*tt77*tt131-tt1*tt197*tt130*tt131-tt1*tt203*tt14*tt131+2*t&
&t1*tt197*tt14*tt77*tt116
tt212 = tt102*X(3,2)+tt100*X(3,1)
tt213 = tt100*X(3,2)+tt102*X(3,1)
tt214 = tt213*X(3,4)+tt212*X(3,3)+tt98*X(2,4)+2*tt95*X(2,3)+tt94+&
&tt90
tt215 = 16*X(2,4)-16*X(2,3)
tt216 = tt59*X(3,2)+tt57*X(3,1)
tt217 = tt57*X(3,2)+tt59*X(3,1)
tt218 = tt217*X(3,4)+tt216*X(3,3)+tt147*tt56+tt215*X(3,1)*X(3,2)+&
&tt147*tt51+tt47*X(2,4)+2*tt40*X(2,3)+tt34+tt30
tt219 = -tt1*tt10*tt77*tt131
tt220 = tt1*(tt123*X(2,4)+2*tt27*X(2,3)+tt120+tt118)*tt136+tt219-&
&tt1*tt214*tt130*tt131-tt1*tt14*tt218*tt131+2*tt1*tt14*tt214*tt77*&
&tt116
tt221 = tt113*X(3,4)+2*tt110*X(3,3)+tt109+tt105
tt222 = tt74*X(3,4)+2*tt67*X(3,3)+tt66+tt62
tt223 = -tt1*tt12*tt77*tt131
tt224 = tt1*(tt123*X(3,4)+2*tt27*X(3,3)+tt126+tt125)*tt136+tt223-&
&tt1*tt221*tt130*tt131-tt1*tt222*tt14*tt131+2*tt1*tt14*tt221*tt77*&
&tt116
tt225 = tt195*X(3,4)+tt196*X(3,3)+tt193*X(2,4)+tt194*X(2,3)+2*tt8&
&0*X(1,4)+tt84+tt82+tt81
tt226 = 16*X(1,3)-16*X(1,4)
tt227 = tt201*X(3,4)+tt202*X(3,3)+tt117*tt56+tt226*X(3,1)*X(3,2)+&
&tt117*tt51+tt199*X(2,4)+tt200*X(2,3)+tt117*tt24+tt226*X(2,1)*X(2,&
&2)+tt117*tt19
tt228 = tt185*X(3,3)
tt229 = tt157*X(3,4)
tt230 = tt208*X(3,3)
tt231 = tt206*X(3,4)
tt232 = tt1*(tt231+tt230+tt64+tt63)*tt136
tt233 = tt232-tt1*(tt229+tt228+tt107+tt106+2*tt2*X(1,4)+tt6+tt4)*&
&tt77*tt131-tt1*tt225*tt130*tt131-tt1*tt227*tt14*tt131+2*tt1*tt225&
&*tt14*tt77*tt116
tt234 = tt212*X(3,4)+tt213*X(3,3)+2*tt95*X(2,4)+tt99+tt97+tt96
tt235 = 16*X(2,3)-16*X(2,4)
tt236 = tt216*X(3,4)+tt217*X(3,3)+tt145*tt56+tt235*X(3,1)*X(3,2)+&
&tt145*tt51+2*tt40*X(2,4)+tt48+tt43+tt42
tt237 = -tt1*tt11*tt77*tt131
tt238 = tt1*(2*tt27*X(2,4)+tt124+tt122+tt121)*tt136+tt237-tt1*tt2&
&34*tt130*tt131-tt1*tt14*tt236*tt131+2*tt1*tt14*tt234*tt77*tt116
tt239 = 2*tt110*X(3,4)+tt114+tt112+tt111
tt240 = 2*tt67*X(3,4)+tt75+tt70+tt69
tt241 = -tt1*tt13*tt77*tt131
tt242 = tt1*(2*tt27*X(3,4)+tt129+tt128+tt127)*tt136+tt241-tt1*tt2&
&39*tt130*tt131-tt1*tt240*tt14*tt131+2*tt1*tt14*tt239*tt77*tt116
tt243 = -tt1*(2*tt49-4*X(2,3)*X(2,4)+2*tt41)*tt77*tt131
tt244 = tt1*(tt135+tt134+tt133+2*tt18)*tt136
tt245 = -tt1*(tt141*X(3,4)+tt140*X(3,3))*tt77*tt131
tt246 = tt1*(tt147*X(3,4)+tt145*X(3,3))*tt136
tt247 = tt246-tt1*tt153*tt160*tt131+tt245-tt1*tt161*tt144*tt131+2&
&*tt1*tt144*tt77*tt160*tt116
tt248 = tt183+tt182-tt1*tt167*tt153*tt131-tt1*tt144*tt171*tt131+2&
&*tt1*tt167*tt144*tt77*tt116
tt249 = tt1*(tt174+tt173+tt172+tt22+tt21+tt20)*tt136-tt1*((-2*tt4&
&9)+4*X(2,3)*X(2,4)-2*tt41)*tt77*tt131-tt1*tt144*tt181*tt131-tt1*t&
&t178*tt153*tt131+2*tt1*tt144*tt178*tt77*tt116
tt250 = -tt1*(tt140*X(3,4)+tt141*X(3,3))*tt77*tt131
tt251 = tt1*(tt145*X(3,4)+tt147*X(3,3))*tt136
tt252 = tt251-tt1*tt153*tt188*tt131+tt250-tt1*tt189*tt144*tt131+2&
&*tt1*tt144*tt77*tt188*tt116
tt253 = -tt1*(tt87*X(2,4)+tt2*X(2,3))*tt77*tt131
tt254 = tt25*X(2,3)
tt255 = tt27*X(2,4)
tt256 = tt1*(tt255+tt254+tt198*X(2,2)+2*tt119*X(2,1))*tt136+tt253&
&-tt1*tt197*tt153*tt131-tt1*tt203*tt144*tt131+2*tt1*tt197*tt144*tt&
&77*tt116
tt257 = tt1*(tt209+tt207+tt28+tt26)*tt136
tt258 = tt257-tt1*(tt205+tt204+tt138*X(2,4)+2*tt100*X(2,3)+tt88+t&
&t86+tt85)*tt77*tt131-tt1*tt214*tt153*tt131-tt1*tt218*tt144*tt131+&
&2*tt1*tt214*tt144*tt77*tt116
tt259 = -tt1*tt142*tt77*tt131
tt260 = tt1*(tt151*X(3,4)+2*tt59*X(3,3)+tt148+tt146)*tt136+tt259-&
&tt1*tt221*tt153*tt131-tt1*tt222*tt144*tt131+2*tt1*tt144*tt221*tt7&
&7*tt116
tt261 = -tt1*(tt2*X(2,4)+tt87*X(2,3))*tt77*tt131
tt262 = tt27*X(2,3)
tt263 = tt25*X(2,4)
tt264 = tt1*(tt263+tt262+tt226*X(2,2)+2*tt117*X(2,1))*tt136+tt261&
&-tt1*tt225*tt153*tt131-tt1*tt227*tt144*tt131+2*tt1*tt225*tt144*tt&
&77*tt116
tt265 = tt1*(tt231+tt230+tt32+tt31)*tt136
tt266 = tt265-tt1*(tt229+tt228+2*tt100*X(2,4)+tt139+tt92+tt91+tt4&
&)*tt77*tt131-tt1*tt234*tt153*tt131-tt1*tt236*tt144*tt131+2*tt1*tt&
&234*tt144*tt77*tt116
tt267 = -tt1*tt143*tt77*tt131
tt268 = tt1*(2*tt59*X(3,4)+tt152+tt150+tt149)*tt136+tt267-tt1*tt2&
&39*tt153*tt131-tt1*tt240*tt144*tt131+2*tt1*tt144*tt239*tt77*tt116
tt269 = -tt1*(2*tt76-4*X(3,3)*X(3,4)+2*tt68)*tt77*tt131
tt270 = 2*tt1*tt50*tt136
tt271 = tt191-tt1*tt171*tt160*tt131+tt190-tt1*tt161*tt167*tt131+2&
&*tt1*tt167*tt77*tt160*tt116
tt272 = tt251-tt1*tt181*tt160*tt131+tt250-tt1*tt161*tt178*tt131+2&
&*tt1*tt178*tt77*tt160*tt116
tt273 = tt1*tt55*tt136-tt1*tt161*tt188*tt131-tt1*tt189*tt160*tt13&
&1-tt1*((-2*tt76)+4*X(3,3)*X(3,4)-2*tt68)*tt77*tt131+2*tt1*tt77*tt&
&160*tt188*tt116
tt274 = -tt1*(tt87*X(3,4)+tt2*X(3,3))*tt77*tt131
tt275 = tt25*X(3,3)
tt276 = tt27*X(3,4)
tt277 = tt1*(tt276+tt275+tt198*X(3,2)+2*tt119*X(3,1))*tt136-tt1*t&
&t203*tt160*tt131+tt274-tt1*tt161*tt197*tt131+2*tt1*tt197*tt77*tt1&
&60*tt116
tt278 = -tt1*(tt102*X(3,4)+tt100*X(3,3))*tt77*tt131
tt279 = tt57*X(3,3)
tt280 = tt59*X(3,4)
tt281 = tt1*(tt280+tt279+tt215*X(3,2)+2*tt147*X(3,1))*tt136-tt1*t&
&t218*tt160*tt131+tt278-tt1*tt161*tt214*tt131+2*tt1*tt214*tt77*tt1&
&60*tt116
tt282 = tt1*tt61*tt136
tt283 = tt282-tt1*tt222*tt160*tt131-tt1*(tt158*X(3,4)+2*tt157*X(3&
&,3)+tt103+tt101+tt88+tt86+tt85)*tt77*tt131-tt1*tt161*tt221*tt131+&
&2*tt1*tt221*tt77*tt160*tt116
tt284 = -tt1*(tt2*X(3,4)+tt87*X(3,3))*tt77*tt131
tt285 = tt27*X(3,3)
tt286 = tt25*X(3,4)
tt287 = tt1*(tt286+tt285+tt226*X(3,2)+2*tt117*X(3,1))*tt136-tt1*t&
&t227*tt160*tt131+tt284-tt1*tt161*tt225*tt131+2*tt1*tt225*tt77*tt1&
&60*tt116
tt288 = -tt1*(tt100*X(3,4)+tt102*X(3,3))*tt77*tt131
tt289 = tt59*X(3,3)
tt290 = tt57*X(3,4)
tt291 = tt1*(tt290+tt289+tt235*X(3,2)+2*tt145*X(3,1))*tt136-tt1*t&
&t236*tt160*tt131+tt288-tt1*tt161*tt234*tt131+2*tt1*tt234*tt77*tt1&
&60*tt116
tt292 = tt1*tt65*tt136
tt293 = tt292-tt1*tt240*tt160*tt131-tt1*(2*tt157*X(3,4)+tt159+tt1&
&07+tt106+tt92+tt91+tt4)*tt77*tt131-tt1*tt161*tt239*tt131+2*tt1*tt&
&239*tt77*tt160*tt116
tt294 = tt155+tt154-tt1*tt167*tt181*tt131-tt1*tt178*tt171*tt131+2&
&*tt1*tt167*tt178*tt77*tt116
tt295 = tt163-tt1*tt171*tt188*tt131+tt162-tt1*tt189*tt167*tt131+2&
&*tt1*tt167*tt77*tt188*tt116
tt296 = tt232-tt1*(tt229+tt228+tt107+tt106+tt165*X(1,4)+2*tt87*X(&
&1,3)+tt4)*tt77*tt131-tt1*tt197*tt171*tt131-tt1*tt203*tt167*tt131+&
&2*tt1*tt197*tt167*tt77*tt116
tt297 = tt1*(tt168*X(2,4)+2*tt25*X(2,3)+tt122+tt121)*tt136+tt237-&
&tt1*tt214*tt171*tt131-tt1*tt167*tt218*tt131+2*tt1*tt167*tt214*tt7&
&7*tt116
tt298 = tt1*(tt168*X(3,4)+2*tt25*X(3,3)+tt128+tt127)*tt136+tt241-&
&tt1*tt221*tt171*tt131-tt1*tt222*tt167*tt131+2*tt1*tt167*tt221*tt7&
&7*tt116
tt299 = tt210-tt1*(tt205+tt204+tt103+tt101+2*tt87*X(1,4)+tt166+tt&
&85)*tt77*tt131-tt1*tt225*tt171*tt131-tt1*tt227*tt167*tt131+2*tt1*&
&tt225*tt167*tt77*tt116
tt300 = tt1*(2*tt25*X(2,4)+tt169+tt120+tt118)*tt136+tt219-tt1*tt2&
&34*tt171*tt131-tt1*tt167*tt236*tt131+2*tt1*tt167*tt234*tt77*tt116
tt301 = tt1*(2*tt25*X(3,4)+tt170+tt126+tt125)*tt136+tt223-tt1*tt2&
&39*tt171*tt131-tt1*tt240*tt167*tt131+2*tt1*tt167*tt239*tt77*tt116
tt302 = tt246-tt1*tt181*tt188*tt131+tt245-tt1*tt189*tt178*tt131+2&
&*tt1*tt178*tt77*tt188*tt116
tt303 = tt1*(tt263+tt262+2*tt119*X(2,2)+tt198*X(2,1))*tt136+tt261&
&-tt1*tt197*tt181*tt131-tt1*tt203*tt178*tt131+2*tt1*tt197*tt178*tt&
&77*tt116
tt304 = tt265-tt1*(tt229+tt228+tt176*X(2,4)+2*tt102*X(2,3)+tt92+t&
&t91+tt4)*tt77*tt131-tt1*tt214*tt181*tt131-tt1*tt218*tt178*tt131+2&
&*tt1*tt214*tt178*tt77*tt116
tt305 = tt1*(tt179*X(3,4)+2*tt57*X(3,3)+tt150+tt149)*tt136+tt267-&
&tt1*tt221*tt181*tt131-tt1*tt222*tt178*tt131+2*tt1*tt178*tt221*tt7&
&7*tt116
tt306 = tt1*(tt255+tt254+2*tt117*X(2,2)+tt226*X(2,1))*tt136+tt253&
&-tt1*tt225*tt181*tt131-tt1*tt227*tt178*tt131+2*tt1*tt225*tt178*tt&
&77*tt116
tt307 = tt257-tt1*(tt205+tt204+2*tt102*X(2,4)+tt177+tt88+tt86+tt8&
&5)*tt77*tt131-tt1*tt234*tt181*tt131-tt1*tt236*tt178*tt131+2*tt1*t&
&t234*tt178*tt77*tt116
tt308 = tt1*(2*tt57*X(3,4)+tt180+tt148+tt146)*tt136+tt259-tt1*tt2&
&39*tt181*tt131-tt1*tt240*tt178*tt131+2*tt1*tt178*tt239*tt77*tt116
tt309 = tt1*(tt286+tt285+2*tt119*X(3,2)+tt198*X(3,1))*tt136-tt1*t&
&t203*tt188*tt131+tt284-tt1*tt189*tt197*tt131+2*tt1*tt197*tt77*tt1&
&88*tt116
tt310 = tt1*(tt290+tt289+2*tt147*X(3,2)+tt215*X(3,1))*tt136-tt1*t&
&t218*tt188*tt131+tt288-tt1*tt189*tt214*tt131+2*tt1*tt214*tt77*tt1&
&88*tt116
tt311 = tt292-tt1*tt222*tt188*tt131-tt1*(tt186*X(3,4)+2*tt185*X(3&
&,3)+tt107+tt106+tt92+tt91+tt4)*tt77*tt131-tt1*tt189*tt221*tt131+2&
&*tt1*tt221*tt77*tt188*tt116
tt312 = tt1*(tt276+tt275+2*tt117*X(3,2)+tt226*X(3,1))*tt136-tt1*t&
&t227*tt188*tt131+tt274-tt1*tt189*tt225*tt131+2*tt1*tt225*tt77*tt1&
&88*tt116
tt313 = tt1*(tt280+tt279+2*tt145*X(3,2)+tt235*X(3,1))*tt136-tt1*t&
&t236*tt188*tt131+tt278-tt1*tt189*tt234*tt131+2*tt1*tt234*tt77*tt1&
&88*tt116
tt314 = tt282-tt1*tt240*tt188*tt131-tt1*(2*tt185*X(3,4)+tt187+tt1&
&03+tt101+tt88+tt86+tt85)*tt77*tt131-tt1*tt189*tt239*tt131+2*tt1*t&
&t239*tt77*tt188*tt116
tt315 = -2*tt1*tt80*tt77*tt131
tt316 = 8*tt51
tt317 = -16*X(3,1)*X(3,2)
tt318 = 8*tt56
tt319 = tt1*(tt318+tt317+tt316+8*tt24-16*X(2,1)*X(2,2)+8*tt19)*tt&
&136
tt320 = -tt1*tt193*tt77*tt131
tt321 = tt1*tt199*tt136
tt322 = tt321+tt320-tt1*tt197*tt218*tt131-tt1*tt203*tt214*tt131+2&
&*tt1*tt197*tt214*tt77*tt116
tt323 = -tt1*tt195*tt77*tt131
tt324 = tt1*tt201*tt136
tt325 = tt324+tt323-tt1*tt203*tt221*tt131-tt1*tt222*tt197*tt131+2&
&*tt1*tt197*tt221*tt77*tt116
tt326 = -8*tt51
tt327 = 16*X(3,1)*X(3,2)
tt328 = -8*tt56
tt329 = tt1*(tt328+tt327+tt326+tt73+tt72+tt71)*tt136-tt1*tt83*tt7&
&7*tt131-tt1*tt225*tt203*tt131-tt1*tt227*tt197*tt131+2*tt1*tt197*t&
&t225*tt77*tt116
tt330 = -tt1*tt194*tt77*tt131
tt331 = tt1*tt200*tt136
tt332 = tt331+tt330-tt1*tt203*tt234*tt131-tt1*tt197*tt236*tt131+2&
&*tt1*tt197*tt234*tt77*tt116
tt333 = -tt1*tt196*tt77*tt131
tt334 = tt1*tt202*tt136
tt335 = tt334+tt333-tt1*tt203*tt239*tt131-tt1*tt240*tt197*tt131+2&
&*tt1*tt197*tt239*tt77*tt116
tt336 = -2*tt1*tt95*tt77*tt131
tt337 = tt1*(tt318+tt317+tt316+2*tt40)*tt136
tt338 = -tt1*tt212*tt77*tt131
tt339 = tt1*tt216*tt136
tt340 = tt339+tt338-tt1*tt218*tt221*tt131-tt1*tt222*tt214*tt131+2&
&*tt1*tt214*tt221*tt77*tt116
tt341 = tt331+tt330-tt1*tt225*tt218*tt131-tt1*tt227*tt214*tt131+2&
&*tt1*tt225*tt214*tt77*tt116
tt342 = tt1*(tt328+tt327+tt326+tt46+tt45+tt44)*tt136-tt1*tt98*tt7&
&7*tt131-tt1*tt234*tt218*tt131-tt1*tt236*tt214*tt131+2*tt1*tt214*t&
&t234*tt77*tt116
tt343 = -tt1*tt213*tt77*tt131
tt344 = tt1*tt217*tt136
tt345 = tt344+tt343-tt1*tt218*tt239*tt131-tt1*tt240*tt214*tt131+2&
&*tt1*tt214*tt239*tt77*tt116
tt346 = -2*tt1*tt110*tt77*tt131
tt347 = 2*tt1*tt67*tt136
tt348 = tt334+tt333-tt1*tt227*tt221*tt131-tt1*tt222*tt225*tt131+2&
&*tt1*tt225*tt221*tt77*tt116
tt349 = tt344+tt343-tt1*tt236*tt221*tt131-tt1*tt222*tt234*tt131+2&
&*tt1*tt234*tt221*tt77*tt116
tt350 = tt1*tt74*tt136-tt1*tt113*tt77*tt131-tt1*tt222*tt239*tt131&
&-tt1*tt240*tt221*tt131+2*tt1*tt221*tt239*tt77*tt116
tt351 = tt321+tt320-tt1*tt227*tt234*tt131-tt1*tt225*tt236*tt131+2&
&*tt1*tt225*tt234*tt77*tt116
tt352 = tt324+tt323-tt1*tt227*tt239*tt131-tt1*tt240*tt225*tt131+2&
&*tt1*tt225*tt239*tt77*tt116
tt353 = tt339+tt338-tt1*tt236*tt239*tt131-tt1*tt240*tt234*tt131+2&
&*tt1*tt234*tt239*tt77*tt116
hes(1,1) = tt137+tt132-2*tt1*tt14*tt130*tt131+2*tt1*tt14**2*tt77*&
&tt116
hes(1,2) = tt156
hes(1,3) = tt164
hes(1,4) = tt175
hes(1,5) = tt184
hes(1,6) = tt192
hes(1,7) = tt211
hes(1,8) = tt220
hes(1,9) = tt224
hes(1,10) = tt233
hes(1,11) = tt238
hes(1,12) = tt242
hes(2,1) = tt156
hes(2,2) = tt244+tt243-2*tt1*tt144*tt153*tt131+2*tt1*tt144**2*tt7&
&7*tt116
hes(2,3) = tt247
hes(2,4) = tt248
hes(2,5) = tt249
hes(2,6) = tt252
hes(2,7) = tt256
hes(2,8) = tt258
hes(2,9) = tt260
hes(2,10) = tt264
hes(2,11) = tt266
hes(2,12) = tt268
hes(3,1) = tt164
hes(3,2) = tt247
hes(3,3) = tt270-2*tt1*tt161*tt160*tt131+tt269+2*tt1*tt77*tt160**&
&2*tt116
hes(3,4) = tt271
hes(3,5) = tt272
hes(3,6) = tt273
hes(3,7) = tt277
hes(3,8) = tt281
hes(3,9) = tt283
hes(3,10) = tt287
hes(3,11) = tt291
hes(3,12) = tt293
hes(4,1) = tt175
hes(4,2) = tt248
hes(4,3) = tt271
hes(4,4) = tt137+tt132-2*tt1*tt167*tt171*tt131+2*tt1*tt167**2*tt7&
&7*tt116
hes(4,5) = tt294
hes(4,6) = tt295
hes(4,7) = tt296
hes(4,8) = tt297
hes(4,9) = tt298
hes(4,10) = tt299
hes(4,11) = tt300
hes(4,12) = tt301
hes(5,1) = tt184
hes(5,2) = tt249
hes(5,3) = tt272
hes(5,4) = tt294
hes(5,5) = tt244+tt243-2*tt1*tt178*tt181*tt131+2*tt1*tt178**2*tt7&
&7*tt116
hes(5,6) = tt302
hes(5,7) = tt303
hes(5,8) = tt304
hes(5,9) = tt305
hes(5,10) = tt306
hes(5,11) = tt307
hes(5,12) = tt308
hes(6,1) = tt192
hes(6,2) = tt252
hes(6,3) = tt273
hes(6,4) = tt295
hes(6,5) = tt302
hes(6,6) = tt270-2*tt1*tt189*tt188*tt131+tt269+2*tt1*tt77*tt188**&
&2*tt116
hes(6,7) = tt309
hes(6,8) = tt310
hes(6,9) = tt311
hes(6,10) = tt312
hes(6,11) = tt313
hes(6,12) = tt314
hes(7,1) = tt211
hes(7,2) = tt256
hes(7,3) = tt277
hes(7,4) = tt296
hes(7,5) = tt303
hes(7,6) = tt309
hes(7,7) = tt319+tt315-2*tt1*tt197*tt203*tt131+2*tt1*tt197**2*tt7&
&7*tt116
hes(7,8) = tt322
hes(7,9) = tt325
hes(7,10) = tt329
hes(7,11) = tt332
hes(7,12) = tt335
hes(8,1) = tt220
hes(8,2) = tt258
hes(8,3) = tt281
hes(8,4) = tt297
hes(8,5) = tt304
hes(8,6) = tt310
hes(8,7) = tt322
hes(8,8) = tt337+tt336-2*tt1*tt214*tt218*tt131+2*tt1*tt214**2*tt7&
&7*tt116
hes(8,9) = tt340
hes(8,10) = tt341
hes(8,11) = tt342
hes(8,12) = tt345
hes(9,1) = tt224
hes(9,2) = tt260
hes(9,3) = tt283
hes(9,4) = tt298
hes(9,5) = tt305
hes(9,6) = tt311
hes(9,7) = tt325
hes(9,8) = tt340
hes(9,9) = tt347+tt346-2*tt1*tt222*tt221*tt131+2*tt1*tt221**2*tt7&
&7*tt116
hes(9,10) = tt348
hes(9,11) = tt349
hes(9,12) = tt350
hes(10,1) = tt233
hes(10,2) = tt264
hes(10,3) = tt287
hes(10,4) = tt299
hes(10,5) = tt306
hes(10,6) = tt312
hes(10,7) = tt329
hes(10,8) = tt341
hes(10,9) = tt348
hes(10,10) = tt319+tt315-2*tt1*tt227*tt225*tt131+2*tt1*tt225**2*t&
&t77*tt116
hes(10,11) = tt351
hes(10,12) = tt352
hes(11,1) = tt238
hes(11,2) = tt266
hes(11,3) = tt291
hes(11,4) = tt300
hes(11,5) = tt307
hes(11,6) = tt313
hes(11,7) = tt332
hes(11,8) = tt342
hes(11,9) = tt349
hes(11,10) = tt351
hes(11,11) = tt337+tt336-2*tt1*tt236*tt234*tt131+2*tt1*tt234**2*t&
&t77*tt116
hes(11,12) = tt353
hes(12,1) = tt242
hes(12,2) = tt268
hes(12,3) = tt293
hes(12,4) = tt301
hes(12,5) = tt308
hes(12,6) = tt314
hes(12,7) = tt335
hes(12,8) = tt345
hes(12,9) = tt350
hes(12,10) = tt352
hes(12,11) = tt353
hes(12,12) = tt347+tt346-2*tt1*tt240*tt239*tt131+2*tt1*tt239**2*t&
&t77*tt116
END 
SUBROUTINE curve_bending(val, X, d1, d2) 
IMPLICIT NONE 
REAL(KIND=8) val(1, 1) 
REAL(KIND=8) X(3, 3) 
REAL(KIND=8) d1(1, 1) 
REAL(KIND=8) d2(1, 1) 
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
tt1 = X(1,2)**2
tt2 = 4*tt1
tt3 = -8*X(1,2)*X(1,3)
tt4 = X(1,3)**2
tt5 = 4*tt4
tt6 = X(2,1)**2
tt7 = -8*X(1,1)*X(1,2)
tt8 = 8*X(1,1)
tt9 = 8*X(1,2)
tt10 = (tt9+tt8)*X(1,3)
tt11 = -8*tt4
tt12 = X(1,1)**2
tt13 = 4*tt12
tt14 = -8*X(1,1)*X(1,3)
tt15 = X(2,2)**2
tt16 = 8*X(1,1)*X(1,2)
tt17 = -8*tt1
tt18 = (tt9-8*X(1,1))*X(1,3)
tt19 = -8*tt12
tt20 = (tt8-8*X(1,2))*X(1,3)
tt21 = X(2,3)**2
tt22 = 4*tt15
tt23 = 4*tt21
tt24 = X(3,1)**2
tt25 = -8*X(2,1)*X(2,2)
tt26 = 8*X(2,1)
tt27 = 8*X(2,2)
tt28 = 4*tt6
tt29 = X(3,2)**2
tt30 = 8*X(2,1)*X(2,2)
tt31 = X(3,3)**2
tt32 = -2*d1(1,1)*d2(1,1)
tt33 = X(1,2)**3
tt34 = 2*d1(1,1)*d2(1,1)
tt35 = -2*X(1,1)*X(1,2)
tt36 = 2*X(1,1)*X(1,2)
tt37 = -2*tt1
tt38 = (2*X(1,2)-2*X(1,1))*X(1,3)
tt39 = 2*tt1
tt40 = (2*X(1,1)-2*X(1,2))*X(1,3)
tt41 = X(2,2)**3
tt42 = -2*X(2,1)*X(2,2)
tt43 = 2*X(2,1)*X(2,2)
tt44 = -2*tt15
tt45 = (2*X(2,2)-2*X(2,1))*X(2,3)
tt46 = 2*tt15
tt47 = (2*X(2,1)-2*X(2,2))*X(2,3)
tt48 = X(3,2)**3
val(1,1) = ((tt22+tt25+tt28+tt2+tt7+tt13)*tt31+(((tt26-8*X(2,2))*&
&X(2,3)+tt30-8*tt6+tt20+tt16+tt19)*X(3,2)+((tt27-8*X(2,1))*X(2,3)-&
&8*tt15+tt30+tt18+tt17+tt16)*X(3,1))*X(3,3)+(tt23-8*X(2,1)*X(2,3)+&
&tt28+tt5+tt14+tt13)*tt29+((-8*tt21)+(tt27+tt26)*X(2,3)+tt25+tt11+&
&tt10+tt7)*X(3,1)*X(3,2)+(tt23-8*X(2,2)*X(2,3)+tt22+tt5+tt3+tt2)*t&
&t24+(tt2+tt7+tt13)*tt21+((tt20+tt16+tt19)*X(2,2)+(tt18+tt17+tt16)&
&*X(2,1))*X(2,3)+(tt5+tt14+tt13)*tt15+(tt11+tt10+tt7)*X(2,1)*X(2,2&
&)+(tt5+tt3+tt2)*tt6)/((d2(1,1)+d1(1,1))*((tt29-2*X(3,1)*X(3,2)+tt&
&24)*tt31+((-2*tt48)+4*X(3,1)*tt29+((-2*tt24)+tt45+tt44+tt43+tt38+&
&tt37+tt36+tt34)*X(3,2)+(tt47+tt46+tt42+tt40+tt39+tt35+tt32)*X(3,1&
&))*X(3,3)+X(3,2)**4-2*X(3,1)*tt48+(tt24+tt47+tt46+tt42+tt40+tt39+&
&tt35+tt32)*tt29+(tt45+tt44+tt43+tt38+tt37+tt36+tt34)*X(3,1)*X(3,2&
&)+(tt15+tt42+tt6)*tt21+((-2*tt41)+4*X(2,1)*tt15+((-2*tt6)+tt38+tt&
&37+tt36+tt34)*X(2,2)+(tt40+tt39+tt35+tt32)*X(2,1))*X(2,3)+X(2,2)*&
&*4-2*X(2,1)*tt41+(tt6+tt40+tt39+tt35+tt32)*tt15+(tt38+tt37+tt36+t&
&t34)*X(2,1)*X(2,2)+(tt1+tt35+tt12)*tt4+((-2*tt33)+4*X(1,1)*tt1+(t&
&t34-2*tt12)*X(1,2)-2*X(1,1)*d1(1,1)*d2(1,1))*X(1,3)+X(1,2)**4-2*X&
&(1,1)*tt33+(tt32+tt12)*tt1+2*X(1,1)*d1(1,1)*d2(1,1)*X(1,2)+d1(1,1&
&)**2*d2(1,1)**2))
END 
SUBROUTINE curve_bending_jac(jac, X, d1, d2) 
IMPLICIT NONE 
REAL(KIND=8) jac(1, 9) 
REAL(KIND=8) X(3, 3) 
REAL(KIND=8) d1(1, 1) 
REAL(KIND=8) d2(1, 1) 
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
tt1 = 1/(d2(1,1)+d1(1,1))
tt2 = X(1,2)**2
tt3 = X(1,2)**3
tt4 = -2*tt3
tt5 = -2*d1(1,1)*d2(1,1)
tt6 = 4*tt2
tt7 = 2*X(1,1)
tt8 = -2*X(1,2)
tt9 = tt8+tt7
tt10 = X(1,3)**2
tt11 = 2*X(1,2)
tt12 = -2*X(1,3)
tt13 = tt12+tt11
tt14 = 2*X(1,3)
tt15 = tt14+tt8
tt16 = X(2,2)**2
tt17 = X(3,2)**2
tt18 = -8*X(1,2)*X(1,3)
tt19 = 4*tt10
tt20 = tt19+tt18+tt6
tt21 = X(2,1)**2
tt22 = -8*X(1,1)*X(1,2)
tt23 = 8*X(1,1)
tt24 = 8*X(1,2)
tt25 = (tt24+tt23)*X(1,3)
tt26 = -8*tt10
tt27 = tt26+tt25+tt22
tt28 = X(1,1)**2
tt29 = 4*tt28
tt30 = -8*X(1,1)*X(1,3)
tt31 = tt19+tt30+tt29
tt32 = 8*X(1,1)*X(1,2)
tt33 = -8*tt2
tt34 = -8*X(1,1)
tt35 = tt24+tt34
tt36 = tt35*X(1,3)
tt37 = tt36+tt33+tt32
tt38 = tt37*X(2,1)
tt39 = -8*tt28
tt40 = -8*X(1,2)
tt41 = tt40+tt23
tt42 = tt41*X(1,3)
tt43 = tt42+tt32+tt39
tt44 = tt43*X(2,2)
tt45 = tt6+tt22+tt29
tt46 = X(2,3)**2
tt47 = 4*tt16
tt48 = 4*tt46
tt49 = tt48-8*X(2,2)*X(2,3)+tt47+tt19+tt18+tt6
tt50 = X(3,1)**2
tt51 = -8*X(2,1)*X(2,2)
tt52 = 8*X(2,1)
tt53 = 8*X(2,2)
tt54 = (-8*tt46)+(tt53+tt52)*X(2,3)+tt51+tt26+tt25+tt22
tt55 = 4*tt21
tt56 = tt48-8*X(2,1)*X(2,3)+tt55+tt19+tt30+tt29
tt57 = 8*X(2,1)*X(2,2)
tt58 = -8*X(2,1)
tt59 = tt53+tt58
tt60 = tt59*X(2,3)-8*tt16+tt57+tt36+tt33+tt32
tt61 = tt60*X(3,1)
tt62 = -8*X(2,2)
tt63 = tt62+tt52
tt64 = tt63*X(2,3)+tt57-8*tt21+tt42+tt32+tt39
tt65 = tt64*X(3,2)
tt66 = tt47+tt51+tt55+tt6+tt22+tt29
tt67 = X(3,3)**2
tt68 = tt66*tt67+(tt65+tt61)*X(3,3)+tt56*tt17+tt54*X(3,1)*X(3,2)+&
&tt49*tt50+tt45*tt46+(tt44+tt38)*X(2,3)+tt31*tt16+tt27*X(2,1)*X(2,&
&2)+tt20*tt21
tt69 = tt5+tt28
tt70 = -2*X(1,1)*d1(1,1)*d2(1,1)
tt71 = -2*tt28
tt72 = 2*d1(1,1)*d2(1,1)
tt73 = (tt72+tt71)*X(1,2)
tt74 = 4*X(1,1)*tt2
tt75 = -2*X(1,1)*X(1,2)
tt76 = tt2+tt75+tt28
tt77 = 2*X(1,1)*X(1,2)
tt78 = -2*tt2
tt79 = -2*X(1,1)
tt80 = tt11+tt79
tt81 = tt80*X(1,3)
tt82 = tt81+tt78+tt77+tt72
tt83 = 2*tt2
tt84 = tt9*X(1,3)
tt85 = tt21+tt84+tt83+tt75+tt5
tt86 = X(2,2)**3
tt87 = (tt84+tt83+tt75+tt5)*X(2,1)
tt88 = -2*tt21
tt89 = (tt88+tt81+tt78+tt77+tt72)*X(2,2)
tt90 = 4*X(2,1)*tt16
tt91 = -2*tt86
tt92 = -2*X(2,1)*X(2,2)
tt93 = tt16+tt92+tt21
tt94 = 2*X(2,1)*X(2,2)
tt95 = -2*tt16
tt96 = -2*X(2,1)
tt97 = 2*X(2,2)
tt98 = tt97+tt96
tt99 = tt98*X(2,3)
tt100 = tt99+tt95+tt94+tt81+tt78+tt77+tt72
tt101 = 2*tt16
tt102 = 2*X(2,1)
tt103 = -2*X(2,2)
tt104 = tt103+tt102
tt105 = tt104*X(2,3)
tt106 = tt50+tt105+tt101+tt92+tt84+tt83+tt75+tt5
tt107 = X(3,2)**3
tt108 = (tt105+tt101+tt92+tt84+tt83+tt75+tt5)*X(3,1)
tt109 = -2*tt50
tt110 = (tt109+tt99+tt95+tt94+tt81+tt78+tt77+tt72)*X(3,2)
tt111 = 4*X(3,1)*tt17
tt112 = -2*tt107
tt113 = tt17-2*X(3,1)*X(3,2)+tt50
tt114 = tt113*tt67+(tt112+tt111+tt110+tt108)*X(3,3)+X(3,2)**4-2*X&
&(3,1)*tt107+tt106*tt17+tt100*X(3,1)*X(3,2)+tt93*tt46+(tt91+tt90+t&
&t89+tt87)*X(2,3)+X(2,2)**4-2*X(2,1)*tt86+tt85*tt16+tt82*X(2,1)*X(&
&2,2)+tt76*tt10+(tt4+tt74+tt73+tt70)*X(1,3)+X(1,2)**4-2*X(1,1)*tt3&
&+tt69*tt2+2*X(1,1)*d1(1,1)*d2(1,1)*X(1,2)+d1(1,1)**2*d2(1,1)**2
tt115 = 1/tt114**2
tt116 = 8*X(1,3)
tt117 = tt116+tt40
tt118 = -8*X(1,3)
tt119 = tt118+tt23
tt120 = tt118+tt24
tt121 = tt116+tt24-16*X(1,1)
tt122 = 1/tt114
tt123 = -2*X(2,3)
tt124 = tt123+tt97
tt125 = 2*X(2,3)
tt126 = tt125+tt103
tt127 = 8*X(2,3)
tt128 = tt127+tt62
tt129 = -8*X(2,3)
tt130 = tt129+tt52
tt131 = tt129+tt53
tt132 = tt14-4*X(1,2)+tt7
tt133 = tt12+4*X(1,2)+tt79
tt134 = tt116+tt34
tt135 = tt116-16*X(1,2)+tt23
tt136 = tt125-4*X(2,2)+tt102
tt137 = tt123+4*X(2,2)+tt96
tt138 = tt127+tt58
tt139 = (-16*X(1,3))+tt24+tt23
jac(1,1) = tt1*(tt41*tt67+(tt121*X(3,2)+tt120*X(3,1))*X(3,3)+tt11&
&9*tt17+tt117*X(3,1)*X(3,2)+tt41*tt46+(tt121*X(2,2)+tt120*X(2,1))*&
&X(2,3)+tt119*tt16+tt117*X(2,1)*X(2,2))*tt122-tt1*((tt13*X(3,2)+tt&
&15*X(3,1))*X(3,3)+tt15*tt17+tt13*X(3,1)*X(3,2)+(tt13*X(2,2)+tt15*&
&X(2,1))*X(2,3)+tt15*tt16+tt13*X(2,1)*X(2,2)+tt9*tt10+(tt6-4*X(1,1&
&)*X(1,2)+tt5)*X(1,3)+tt4+2*X(1,1)*tt2+2*d1(1,1)*d2(1,1)*X(1,2))*t&
&t68*tt115
jac(1,2) = tt1*(tt63*tt67+((tt127+tt53-16*X(2,1))*X(3,2)+tt131*X(&
&3,1))*X(3,3)+tt130*tt17+tt128*X(3,1)*X(3,2)+tt37*X(2,3)+tt27*X(2,&
&2)+2*tt20*X(2,1))*tt122-tt1*((tt124*X(3,2)+tt126*X(3,1))*X(3,3)+t&
&t126*tt17+tt124*X(3,1)*X(3,2)+tt104*tt46+(tt47-4*X(2,1)*X(2,2)+tt&
&84+tt83+tt75+tt5)*X(2,3)+tt91+2*X(2,1)*tt16+tt82*X(2,2))*tt68*tt1&
&15
jac(1,3) = tt1*(tt60*X(3,3)+tt54*X(3,2)+2*tt49*X(3,1))*tt122-tt1*&
&tt68*((2*X(3,1)-2*X(3,2))*tt67+(4*tt17-4*X(3,1)*X(3,2)+tt105+tt10&
&1+tt92+tt84+tt83+tt75+tt5)*X(3,3)+tt112+2*X(3,1)*tt17+tt100*X(3,2&
&))*tt115
jac(1,4) = tt1*(tt35*tt67+(tt119*X(3,2)+tt135*X(3,1))*X(3,3)+tt13&
&4*X(3,1)*X(3,2)+tt120*tt50+tt35*tt46+(tt119*X(2,2)+tt135*X(2,1))*&
&X(2,3)+tt134*X(2,1)*X(2,2)+tt120*tt21)*tt122-tt1*((tt132*X(3,2)+t&
&t133*X(3,1))*X(3,3)+tt133*tt17+tt132*X(3,1)*X(3,2)+(tt132*X(2,2)+&
&tt133*X(2,1))*X(2,3)+tt133*tt16+tt132*X(2,1)*X(2,2)+tt80*tt10+((-&
&6*tt2)+tt32+tt72+tt71)*X(1,3)+4*tt3-6*X(1,1)*tt2+2*tt69*X(1,2)+2*&
&X(1,1)*d1(1,1)*d2(1,1))*tt68*tt115
jac(1,5) = tt1*(tt59*tt67+(tt130*X(3,2)+(tt127-16*X(2,2)+tt52)*X(&
&3,1))*X(3,3)+tt138*X(3,1)*X(3,2)+tt131*tt50+tt43*X(2,3)+2*tt31*X(&
&2,2)+tt27*X(2,1))*tt122-tt1*((tt136*X(3,2)+tt137*X(3,1))*X(3,3)+t&
&t137*tt17+tt136*X(3,1)*X(3,2)+tt98*tt46+((-6*tt16)+tt57+tt88+tt81&
&+tt78+tt77+tt72)*X(2,3)+4*tt86-6*X(2,1)*tt16+2*tt85*X(2,2)+tt82*X&
&(2,1))*tt68*tt115
jac(1,6) = tt1*(tt64*X(3,3)+2*tt56*X(3,2)+tt54*X(3,1))*tt122-tt1*&
&tt68*((2*X(3,2)-2*X(3,1))*tt67+((-6*tt17)+8*X(3,1)*X(3,2)+tt109+t&
&t99+tt95+tt94+tt81+tt78+tt77+tt72)*X(3,3)+4*tt107-6*X(3,1)*tt17+2&
&*tt106*X(3,2)+tt100*X(3,1))*tt115
jac(1,7) = tt1*((tt41*X(3,2)+tt35*X(3,1))*X(3,3)+tt134*tt17+tt139&
&*X(3,1)*X(3,2)+tt117*tt50+(tt41*X(2,2)+tt35*X(2,1))*X(2,3)+tt134*&
&tt16+tt139*X(2,1)*X(2,2)+tt117*tt21)*tt122-tt1*((tt80*X(3,2)+tt9*&
&X(3,1))*X(3,3)+tt9*tt17+tt80*X(3,1)*X(3,2)+(tt80*X(2,2)+tt9*X(2,1&
&))*X(2,3)+tt9*tt16+tt80*X(2,1)*X(2,2)+2*tt76*X(1,3)+tt4+tt74+tt73&
&+tt70)*tt68*tt115
jac(1,8) = tt1*((tt63*X(3,2)+tt59*X(3,1))*X(3,3)+tt138*tt17+((-16&
&*X(2,3))+tt53+tt52)*X(3,1)*X(3,2)+tt128*tt50+2*tt45*X(2,3)+tt44+t&
&t38)*tt122-tt1*((tt98*X(3,2)+tt104*X(3,1))*X(3,3)+tt104*tt17+tt98&
&*X(3,1)*X(3,2)+2*tt93*X(2,3)+tt91+tt90+tt89+tt87)*tt68*tt115
jac(1,9) = tt1*(2*tt66*X(3,3)+tt65+tt61)*tt122-tt1*(2*tt113*X(3,3&
&)+tt112+tt111+tt110+tt108)*tt68*tt115
END 
SUBROUTINE curve_bending_hes(hes, X, d1, d2) 
IMPLICIT NONE 
REAL(KIND=8) hes(9, 9) 
REAL(KIND=8) X(3, 3) 
REAL(KIND=8) d1(1, 1) 
REAL(KIND=8) d2(1, 1) 
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
tt1 = 1/(d2(1,1)+d1(1,1))
tt2 = X(1,2)**2
tt3 = X(1,2)**3
tt4 = -2*tt3
tt5 = -2*d1(1,1)*d2(1,1)
tt6 = -4*X(1,1)*X(1,2)
tt7 = 4*tt2
tt8 = 2*X(1,1)
tt9 = -2*X(1,2)
tt10 = tt9+tt8
tt11 = X(1,3)**2
tt12 = 2*X(1,2)
tt13 = -2*X(1,3)
tt14 = tt13+tt12
tt15 = 2*X(1,3)
tt16 = tt15+tt9
tt17 = X(2,2)**2
tt18 = tt14*X(2,2)
tt19 = tt18+tt16*X(2,1)
tt20 = X(3,2)**2
tt21 = tt14*X(3,2)
tt22 = tt21+tt16*X(3,1)
tt23 = tt22*X(3,3)+tt16*tt20+tt14*X(3,1)*X(3,2)+tt19*X(2,3)+tt16*&
&tt17+tt14*X(2,1)*X(2,2)+tt10*tt11+(tt7+tt6+tt5)*X(1,3)+tt4+2*X(1,&
&1)*tt2+2*d1(1,1)*d2(1,1)*X(1,2)
tt24 = -8*X(1,2)*X(1,3)
tt25 = 4*tt11
tt26 = tt25+tt24+tt7
tt27 = X(2,1)**2
tt28 = -8*X(1,1)*X(1,2)
tt29 = 8*X(1,1)
tt30 = 8*X(1,2)
tt31 = (tt30+tt29)*X(1,3)
tt32 = -8*tt11
tt33 = tt32+tt31+tt28
tt34 = X(1,1)**2
tt35 = 4*tt34
tt36 = -8*X(1,1)*X(1,3)
tt37 = tt25+tt36+tt35
tt38 = 8*X(1,1)*X(1,2)
tt39 = -8*tt2
tt40 = -8*X(1,1)
tt41 = tt30+tt40
tt42 = tt41*X(1,3)
tt43 = tt42+tt39+tt38
tt44 = tt43*X(2,1)
tt45 = -8*tt34
tt46 = -8*X(1,2)
tt47 = tt46+tt29
tt48 = tt47*X(1,3)
tt49 = tt48+tt38+tt45
tt50 = tt49*X(2,2)
tt51 = tt7+tt28+tt35
tt52 = X(2,3)**2
tt53 = 4*tt17
tt54 = 4*tt52
tt55 = tt54-8*X(2,2)*X(2,3)+tt53+tt25+tt24+tt7
tt56 = X(3,1)**2
tt57 = -8*X(2,1)*X(2,2)
tt58 = 8*X(2,1)
tt59 = 8*X(2,2)
tt60 = (tt59+tt58)*X(2,3)
tt61 = -8*tt52
tt62 = tt61+tt60+tt57+tt32+tt31+tt28
tt63 = 4*tt27
tt64 = tt54-8*X(2,1)*X(2,3)+tt63+tt25+tt36+tt35
tt65 = 8*X(2,1)*X(2,2)
tt66 = -8*tt17
tt67 = -8*X(2,1)
tt68 = tt59+tt67
tt69 = tt68*X(2,3)
tt70 = tt69+tt66+tt65+tt42+tt39+tt38
tt71 = tt70*X(3,1)
tt72 = -8*tt27
tt73 = -8*X(2,2)
tt74 = tt73+tt58
tt75 = tt74*X(2,3)
tt76 = tt75+tt65+tt72+tt48+tt38+tt45
tt77 = tt76*X(3,2)
tt78 = tt53+tt57+tt63+tt7+tt28+tt35
tt79 = X(3,3)**2
tt80 = tt78*tt79+(tt77+tt71)*X(3,3)+tt64*tt20+tt62*X(3,1)*X(3,2)+&
&tt55*tt56+tt51*tt52+(tt50+tt44)*X(2,3)+tt37*tt17+tt33*X(2,1)*X(2,&
&2)+tt26*tt27
tt81 = tt5+tt34
tt82 = -2*X(1,1)*d1(1,1)*d2(1,1)
tt83 = -2*tt34
tt84 = 2*d1(1,1)*d2(1,1)
tt85 = (tt84+tt83)*X(1,2)
tt86 = 4*X(1,1)*tt2
tt87 = -2*X(1,1)*X(1,2)
tt88 = tt2+tt87+tt34
tt89 = 2*X(1,1)*X(1,2)
tt90 = -2*tt2
tt91 = -2*X(1,1)
tt92 = tt12+tt91
tt93 = tt92*X(1,3)
tt94 = tt93+tt90+tt89+tt84
tt95 = 2*tt2
tt96 = tt10*X(1,3)
tt97 = tt27+tt96+tt95+tt87+tt5
tt98 = X(2,2)**3
tt99 = (tt96+tt95+tt87+tt5)*X(2,1)
tt100 = -2*tt27
tt101 = (tt100+tt93+tt90+tt89+tt84)*X(2,2)
tt102 = 4*X(2,1)*tt17
tt103 = -2*tt98
tt104 = -2*X(2,1)*X(2,2)
tt105 = tt17+tt104+tt27
tt106 = 2*X(2,1)*X(2,2)
tt107 = -2*tt17
tt108 = -2*X(2,1)
tt109 = 2*X(2,2)
tt110 = tt109+tt108
tt111 = tt110*X(2,3)
tt112 = tt111+tt107+tt106+tt93+tt90+tt89+tt84
tt113 = 2*tt17
tt114 = 2*X(2,1)
tt115 = -2*X(2,2)
tt116 = tt115+tt114
tt117 = tt116*X(2,3)
tt118 = tt56+tt117+tt113+tt104+tt96+tt95+tt87+tt5
tt119 = X(3,2)**3
tt120 = (tt117+tt113+tt104+tt96+tt95+tt87+tt5)*X(3,1)
tt121 = -2*tt56
tt122 = (tt121+tt111+tt107+tt106+tt93+tt90+tt89+tt84)*X(3,2)
tt123 = 4*X(3,1)*tt20
tt124 = -2*tt119
tt125 = -2*X(3,1)*X(3,2)
tt126 = tt20+tt125+tt56
tt127 = tt126*tt79+(tt124+tt123+tt122+tt120)*X(3,3)+X(3,2)**4-2*X&
&(3,1)*tt119+tt118*tt20+tt112*X(3,1)*X(3,2)+tt105*tt52+(tt103+tt10&
&2+tt101+tt99)*X(2,3)+X(2,2)**4-2*X(2,1)*tt98+tt97*tt17+tt94*X(2,1&
&)*X(2,2)+tt88*tt11+(tt4+tt86+tt85+tt82)*X(1,3)+X(1,2)**4-2*X(1,1)&
&*tt3+tt81*tt2+2*X(1,1)*d1(1,1)*d2(1,1)*X(1,2)+d1(1,1)**2*d2(1,1)*&
&*2
tt128 = 1/tt127**3
tt129 = 8*X(1,3)
tt130 = tt129+tt46
tt131 = -8*X(1,3)
tt132 = tt131+tt29
tt133 = tt131+tt30
tt134 = tt133*X(2,1)
tt135 = tt129+tt30-16*X(1,1)
tt136 = tt135*X(2,2)
tt137 = tt133*X(3,1)
tt138 = tt135*X(3,2)
tt139 = tt47*tt79+(tt138+tt137)*X(3,3)+tt132*tt20+tt130*X(3,1)*X(&
&3,2)+tt47*tt52+(tt136+tt134)*X(2,3)+tt132*tt17+tt130*X(2,1)*X(2,2&
&)
tt140 = 1/tt127**2
tt141 = 2*tt11
tt142 = 8*tt17
tt143 = 8*tt52
tt144 = 8*tt20
tt145 = -16*X(3,2)*X(3,3)
tt146 = 8*tt79
tt147 = 1/tt127
tt148 = -4*X(2,1)*X(2,2)
tt149 = -2*X(2,3)
tt150 = tt149+tt109
tt151 = 2*X(2,3)
tt152 = tt151+tt115
tt153 = tt150*X(3,2)
tt154 = tt153+tt152*X(3,1)
tt155 = tt154*X(3,3)+tt152*tt20+tt150*X(3,1)*X(3,2)+tt116*tt52+(t&
&t53+tt148+tt96+tt95+tt87+tt5)*X(2,3)+tt103+2*X(2,1)*tt17+tt94*X(2&
&,2)
tt156 = 8*X(2,3)
tt157 = tt156+tt73
tt158 = -8*X(2,3)
tt159 = tt158+tt58
tt160 = tt158+tt59
tt161 = tt160*X(3,1)
tt162 = tt156+tt59-16*X(2,1)
tt163 = tt162*X(3,2)
tt164 = tt74*tt79+(tt163+tt161)*X(3,3)+tt159*tt20+tt157*X(3,1)*X(&
&3,2)+tt43*X(2,3)+tt33*X(2,2)+2*tt26*X(2,1)
tt165 = tt1*(tt133*X(2,3)+tt130*X(2,2))*tt147-tt1*(tt16*X(2,3)+tt&
&18)*tt80*tt140-tt1*tt23*tt164*tt140-tt1*tt155*tt139*tt140+2*tt1*t&
&t23*tt155*tt80*tt128
tt166 = -4*X(3,1)*X(3,2)
tt167 = 4*tt20
tt168 = 2*X(3,1)-2*X(3,2)
tt169 = tt168*tt79+(tt167+tt166+tt117+tt113+tt104+tt96+tt95+tt87+&
&tt5)*X(3,3)+tt124+2*X(3,1)*tt20+tt112*X(3,2)
tt170 = tt70*X(3,3)+tt62*X(3,2)+2*tt55*X(3,1)
tt171 = tt1*(tt133*X(3,3)+tt130*X(3,2))*tt147-tt1*tt139*tt169*tt1&
&40-tt1*(tt16*X(3,3)+tt21)*tt80*tt140-tt1*tt170*tt23*tt140+2*tt1*t&
&t23*tt80*tt169*tt128
tt172 = -6*tt2
tt173 = tt15-4*X(1,2)+tt8
tt174 = tt13+4*X(1,2)+tt91
tt175 = tt173*X(2,2)
tt176 = tt175+tt174*X(2,1)
tt177 = tt173*X(3,2)
tt178 = tt177+tt174*X(3,1)
tt179 = tt178*X(3,3)+tt174*tt20+tt173*X(3,1)*X(3,2)+tt176*X(2,3)+&
&tt174*tt17+tt173*X(2,1)*X(2,2)+tt92*tt11+(tt172+tt38+tt84+tt83)*X&
&(1,3)+4*tt3-6*X(1,1)*tt2+2*tt81*X(1,2)+2*X(1,1)*d1(1,1)*d2(1,1)
tt180 = tt129+tt40
tt181 = tt129-16*X(1,2)+tt29
tt182 = tt181*X(2,1)
tt183 = tt132*X(2,2)
tt184 = tt181*X(3,1)
tt185 = tt132*X(3,2)
tt186 = tt41*tt79+(tt185+tt184)*X(3,3)+tt180*X(3,1)*X(3,2)+tt133*&
&tt56+tt41*tt52+(tt183+tt182)*X(2,3)+tt180*X(2,1)*X(2,2)+tt133*tt2&
&7
tt187 = 2*X(3,1)*X(3,2)
tt188 = -2*tt20
tt189 = 2*X(3,2)-2*X(3,1)
tt190 = tt189*X(3,3)
tt191 = -8*X(3,1)*X(3,2)
tt192 = 8*X(3,1)
tt193 = 8*X(3,2)
tt194 = (tt193+tt192)*X(3,3)
tt195 = -8*tt79
tt196 = tt1*(tt195+tt194+tt191+tt61+tt60+tt57)*tt147-tt1*(tt190+t&
&t188+tt187+tt111+tt107+tt106-2*tt11+(tt30-4*X(1,1))*X(1,3)+tt172+&
&4*X(1,1)*X(1,2)+tt84)*tt80*tt140-tt1*tt23*tt186*tt140-tt1*tt179*t&
&t139*tt140+2*tt1*tt23*tt179*tt80*tt128
tt197 = -6*tt17
tt198 = -4*X(2,2)
tt199 = tt151+tt198+tt114
tt200 = tt149+4*X(2,2)+tt108
tt201 = tt199*X(3,2)
tt202 = tt201+tt200*X(3,1)
tt203 = tt202*X(3,3)+tt200*tt20+tt199*X(3,1)*X(3,2)+tt110*tt52+(t&
&t197+tt65+tt100+tt93+tt90+tt89+tt84)*X(2,3)+4*tt98-6*X(2,1)*tt17+&
&2*tt97*X(2,2)+tt94*X(2,1)
tt204 = tt156+tt67
tt205 = tt156-16*X(2,2)+tt58
tt206 = tt205*X(3,1)
tt207 = tt159*X(3,2)
tt208 = tt68*tt79+(tt207+tt206)*X(3,3)+tt204*X(3,1)*X(3,2)+tt160*&
&tt56+tt49*X(2,3)+2*tt37*X(2,2)+tt33*X(2,1)
tt209 = tt1*(tt135*X(2,3)+2*tt132*X(2,2)+tt130*X(2,1))*tt147-tt1*&
&(tt14*X(2,3)+2*tt16*X(2,2)+tt14*X(2,1))*tt80*tt140-tt1*tt23*tt208&
&*tt140-tt1*tt203*tt139*tt140+2*tt1*tt23*tt203*tt80*tt128
tt210 = 8*X(3,1)*X(3,2)
tt211 = -6*tt20
tt212 = tt189*tt79+(tt211+tt210+tt121+tt111+tt107+tt106+tt93+tt90&
&+tt89+tt84)*X(3,3)+4*tt119-6*X(3,1)*tt20+2*tt118*X(3,2)+tt112*X(3&
&,1)
tt213 = tt76*X(3,3)+2*tt64*X(3,2)+tt62*X(3,1)
tt214 = tt1*(tt135*X(3,3)+2*tt132*X(3,2)+tt130*X(3,1))*tt147-tt1*&
&tt139*tt212*tt140-tt1*(tt14*X(3,3)+2*tt16*X(3,2)+tt14*X(3,1))*tt8&
&0*tt140-tt1*tt213*tt23*tt140+2*tt1*tt23*tt80*tt212*tt128
tt215 = tt92*X(2,2)
tt216 = tt215+tt10*X(2,1)
tt217 = tt92*X(3,2)
tt218 = tt217+tt10*X(3,1)
tt219 = tt218*X(3,3)+tt10*tt20+tt92*X(3,1)*X(3,2)+tt216*X(2,3)+tt&
&10*tt17+tt92*X(2,1)*X(2,2)+2*tt88*X(1,3)+tt4+tt86+tt85+tt82
tt220 = (-16*X(1,3))+tt30+tt29
tt221 = tt47*X(2,2)+tt41*X(2,1)
tt222 = tt47*X(3,2)+tt41*X(3,1)
tt223 = tt222*X(3,3)+tt180*tt20+tt220*X(3,1)*X(3,2)+tt130*tt56+tt&
&221*X(2,3)+tt180*tt17+tt220*X(2,1)*X(2,2)+tt130*tt27
tt224 = 2*tt20
tt225 = tt168*X(3,3)
tt226 = -8*tt20
tt227 = (tt193-8*X(3,1))*X(3,3)
tt228 = tt1*(tt227+tt226+tt210+tt69+tt66+tt65)*tt147-tt1*(tt225+t&
&t224+tt125+tt117+tt113+tt104+2*tt10*X(1,3)+tt7+tt6+tt5)*tt80*tt14&
&0-tt1*tt219*tt139*tt140-tt1*tt223*tt23*tt140+2*tt1*tt219*tt23*tt8&
&0*tt128
tt229 = tt110*X(3,2)
tt230 = tt229+tt116*X(3,1)
tt231 = tt230*X(3,3)+tt116*tt20+tt110*X(3,1)*X(3,2)+2*tt105*X(2,3&
&)+tt103+tt102+tt101+tt99
tt232 = (-16*X(2,3))+tt59+tt58
tt233 = tt74*X(3,2)+tt68*X(3,1)
tt234 = tt233*X(3,3)+tt204*tt20+tt232*X(3,1)*X(3,2)+tt157*tt56+2*&
&tt51*X(2,3)+tt50+tt44
tt235 = tt1*(2*tt47*X(2,3)+tt136+tt134)*tt147-tt1*tt19*tt80*tt140&
&-tt1*tt231*tt139*tt140-tt1*tt23*tt234*tt140+2*tt1*tt23*tt231*tt80&
&*tt128
tt236 = 2*tt126*X(3,3)+tt124+tt123+tt122+tt120
tt237 = 2*tt78*X(3,3)+tt77+tt71
tt238 = tt1*(2*tt47*X(3,3)+tt138+tt137)*tt147-tt1*tt22*tt80*tt140&
&-tt1*tt236*tt139*tt140-tt1*tt237*tt23*tt140+2*tt1*tt23*tt236*tt80&
&*tt128
tt239 = 2*tt52
tt240 = tt1*(tt160*X(3,3)+tt157*X(3,2))*tt147-tt1*tt164*tt169*tt1&
&40-tt1*(tt152*X(3,3)+tt153)*tt80*tt140-tt1*tt170*tt155*tt140+2*tt&
&1*tt155*tt80*tt169*tt128
tt241 = tt1*(tt181*X(2,3)+tt180*X(2,2)+2*tt133*X(2,1))*tt147-tt1*&
&(tt174*X(2,3)+tt175)*tt80*tt140-tt1*tt179*tt164*tt140-tt1*tt155*t&
&t186*tt140+2*tt1*tt179*tt155*tt80*tt128
tt242 = tt1*(tt195+tt194+tt191+tt32+tt31+tt28)*tt147-tt1*(tt190+t&
&t188+tt187-2*tt52+(tt59-4*X(2,1))*X(2,3)+tt197+4*X(2,1)*X(2,2)+tt&
&93+tt90+tt89+tt84)*tt80*tt140-tt1*tt155*tt208*tt140-tt1*tt203*tt1&
&64*tt140+2*tt1*tt155*tt203*tt80*tt128
tt243 = tt1*(tt162*X(3,3)+2*tt159*X(3,2)+tt157*X(3,1))*tt147-tt1*&
&tt164*tt212*tt140-tt1*(tt150*X(3,3)+2*tt152*X(3,2)+tt150*X(3,1))*&
&tt80*tt140-tt1*tt213*tt155*tt140+2*tt1*tt155*tt80*tt212*tt128
tt244 = tt1*(tt41*X(2,3)+tt220*X(2,2)+2*tt130*X(2,1))*tt147-tt1*(&
&tt10*X(2,3)+tt215)*tt80*tt140-tt1*tt219*tt164*tt140-tt1*tt223*tt1&
&55*tt140+2*tt1*tt219*tt155*tt80*tt128
tt245 = tt1*(tt227+tt226+tt210+tt42+tt39+tt38)*tt147-tt1*(tt225+t&
&t224+tt125+2*tt116*X(2,3)+tt53+tt148+tt96+tt95+tt87+tt5)*tt80*tt1&
&40-tt1*tt231*tt164*tt140-tt1*tt234*tt155*tt140+2*tt1*tt231*tt155*&
&tt80*tt128
tt246 = tt1*(2*tt74*X(3,3)+tt163+tt161)*tt147-tt1*tt154*tt80*tt14&
&0-tt1*tt236*tt164*tt140-tt1*tt237*tt155*tt140+2*tt1*tt155*tt236*t&
&t80*tt128
tt247 = 2*tt79
tt248 = tt1*(tt181*X(3,3)+tt180*X(3,2)+2*tt133*X(3,1))*tt147-tt1*&
&tt186*tt169*tt140-tt1*(tt174*X(3,3)+tt177)*tt80*tt140-tt1*tt170*t&
&t179*tt140+2*tt1*tt179*tt80*tt169*tt128
tt249 = tt1*(tt205*X(3,3)+tt204*X(3,2)+2*tt160*X(3,1))*tt147-tt1*&
&tt208*tt169*tt140-tt1*(tt200*X(3,3)+tt201)*tt80*tt140-tt1*tt170*t&
&t203*tt140+2*tt1*tt203*tt80*tt169*tt128
tt250 = tt1*tt62*tt147-tt1*tt170*tt212*tt140-tt1*tt213*tt169*tt14&
&0-tt1*((-2*tt79)+(tt193-4*X(3,1))*X(3,3)+tt211+4*X(3,1)*X(3,2)+tt&
&111+tt107+tt106+tt93+tt90+tt89+tt84)*tt80*tt140+2*tt1*tt80*tt169*&
&tt212*tt128
tt251 = tt1*(tt41*X(3,3)+tt220*X(3,2)+2*tt130*X(3,1))*tt147-tt1*t&
&t223*tt169*tt140-tt1*(tt10*X(3,3)+tt217)*tt80*tt140-tt1*tt170*tt2&
&19*tt140+2*tt1*tt219*tt80*tt169*tt128
tt252 = tt1*(tt68*X(3,3)+tt232*X(3,2)+2*tt157*X(3,1))*tt147-tt1*t&
&t234*tt169*tt140-tt1*(tt116*X(3,3)+tt229)*tt80*tt140-tt1*tt170*tt&
&231*tt140+2*tt1*tt231*tt80*tt169*tt128
tt253 = tt1*tt70*tt147-tt1*tt237*tt169*tt140-tt1*(2*tt168*X(3,3)+&
&tt167+tt166+tt117+tt113+tt104+tt96+tt95+tt87+tt5)*tt80*tt140-tt1*&
&tt170*tt236*tt140+2*tt1*tt236*tt80*tt169*tt128
tt254 = (4*X(3,1)-4*X(3,2))*X(3,3)
tt255 = 8*tt27
tt256 = 8*tt56
tt257 = -16*X(3,1)*X(3,3)
tt258 = tt1*(tt132*X(2,3)+tt180*X(2,1))*tt147-tt1*(tt173*X(2,3)+2&
&*tt174*X(2,2)+tt173*X(2,1))*tt80*tt140-tt1*tt179*tt208*tt140-tt1*&
&tt203*tt186*tt140+2*tt1*tt179*tt203*tt80*tt128
tt259 = tt1*(tt132*X(3,3)+tt180*X(3,1))*tt147-tt1*tt186*tt212*tt1&
&40-tt1*(tt173*X(3,3)+2*tt174*X(3,2)+tt173*X(3,1))*tt80*tt140-tt1*&
&tt213*tt179*tt140+2*tt1*tt179*tt80*tt212*tt128
tt260 = -8*tt56
tt261 = (tt192-8*X(3,2))*X(3,3)
tt262 = tt1*(tt261+tt210+tt260+tt75+tt65+tt72)*tt147-tt1*(tt190+t&
&t188+tt187+tt111+tt107+tt106+2*tt92*X(1,3)+tt172+tt38+tt84+tt83)*&
&tt80*tt140-tt1*tt219*tt186*tt140-tt1*tt223*tt179*tt140+2*tt1*tt21&
&9*tt179*tt80*tt128
tt263 = tt1*(2*tt41*X(2,3)+tt183+tt182)*tt147-tt1*tt176*tt80*tt14&
&0-tt1*tt231*tt186*tt140-tt1*tt179*tt234*tt140+2*tt1*tt179*tt231*t&
&t80*tt128
tt264 = tt1*(2*tt41*X(3,3)+tt185+tt184)*tt147-tt1*tt178*tt80*tt14&
&0-tt1*tt236*tt186*tt140-tt1*tt237*tt179*tt140+2*tt1*tt179*tt236*t&
&t80*tt128
tt265 = tt1*(tt159*X(3,3)+tt204*X(3,1))*tt147-tt1*tt208*tt212*tt1&
&40-tt1*(tt199*X(3,3)+2*tt200*X(3,2)+tt199*X(3,1))*tt80*tt140-tt1*&
&tt213*tt203*tt140+2*tt1*tt203*tt80*tt212*tt128
tt266 = tt1*(tt47*X(2,3)+2*tt180*X(2,2)+tt220*X(2,1))*tt147-tt1*(&
&tt92*X(2,3)+2*tt10*X(2,2)+tt92*X(2,1))*tt80*tt140-tt1*tt219*tt208&
&*tt140-tt1*tt223*tt203*tt140+2*tt1*tt219*tt203*tt80*tt128
tt267 = tt1*(tt261+tt210+tt260+tt48+tt38+tt45)*tt147-tt1*(tt190+t&
&t188+tt187+2*tt110*X(2,3)+tt197+tt65+tt100+tt93+tt90+tt89+tt84)*t&
&t80*tt140-tt1*tt231*tt208*tt140-tt1*tt234*tt203*tt140+2*tt1*tt231&
&*tt203*tt80*tt128
tt268 = tt1*(2*tt68*X(3,3)+tt207+tt206)*tt147-tt1*tt202*tt80*tt14&
&0-tt1*tt236*tt208*tt140-tt1*tt237*tt203*tt140+2*tt1*tt203*tt236*t&
&t80*tt128
tt269 = tt1*(tt47*X(3,3)+2*tt180*X(3,2)+tt220*X(3,1))*tt147-tt1*t&
&t223*tt212*tt140-tt1*(tt92*X(3,3)+2*tt10*X(3,2)+tt92*X(3,1))*tt80&
&*tt140-tt1*tt213*tt219*tt140+2*tt1*tt219*tt80*tt212*tt128
tt270 = tt1*(tt74*X(3,3)+2*tt204*X(3,2)+tt232*X(3,1))*tt147-tt1*t&
&t234*tt212*tt140-tt1*(tt110*X(3,3)+2*tt116*X(3,2)+tt110*X(3,1))*t&
&t80*tt140-tt1*tt213*tt231*tt140+2*tt1*tt231*tt80*tt212*tt128
tt271 = tt1*tt76*tt147-tt1*tt237*tt212*tt140-tt1*(2*tt189*X(3,3)+&
&tt211+tt210+tt121+tt111+tt107+tt106+tt93+tt90+tt89+tt84)*tt80*tt1&
&40-tt1*tt213*tt236*tt140+2*tt1*tt236*tt80*tt212*tt128
tt272 = -16*X(3,1)*X(3,2)
tt273 = tt1*tt221*tt147-tt1*tt216*tt80*tt140-tt1*tt223*tt231*tt14&
&0-tt1*tt219*tt234*tt140+2*tt1*tt219*tt231*tt80*tt128
tt274 = tt1*tt222*tt147-tt1*tt218*tt80*tt140-tt1*tt223*tt236*tt14&
&0-tt1*tt237*tt219*tt140+2*tt1*tt219*tt236*tt80*tt128
tt275 = tt1*tt233*tt147-tt1*tt230*tt80*tt140-tt1*tt234*tt236*tt14&
&0-tt1*tt237*tt231*tt140+2*tt1*tt231*tt236*tt80*tt128
hes(1,1) = tt1*(tt146+tt145+tt144+tt143-16*X(2,2)*X(2,3)+tt142)*t&
&t147-tt1*(tt141-4*X(1,2)*X(1,3)+tt95)*tt80*tt140-2*tt1*tt23*tt139&
&*tt140+2*tt1*tt23**2*tt80*tt128
hes(1,2) = tt165
hes(1,3) = tt171
hes(1,4) = tt196
hes(1,5) = tt209
hes(1,6) = tt214
hes(1,7) = tt228
hes(1,8) = tt235
hes(1,9) = tt238
hes(2,1) = tt165
hes(2,2) = tt1*(tt146+tt145+tt144+2*tt26)*tt147-tt1*(tt239-4*X(2,&
&2)*X(2,3)+tt113)*tt80*tt140-2*tt1*tt155*tt164*tt140+2*tt1*tt155**&
&2*tt80*tt128
hes(2,3) = tt240
hes(2,4) = tt241
hes(2,5) = tt242
hes(2,6) = tt243
hes(2,7) = tt244
hes(2,8) = tt245
hes(2,9) = tt246
hes(3,1) = tt171
hes(3,2) = tt240
hes(3,3) = 2*tt1*tt55*tt147-2*tt1*tt170*tt169*tt140-tt1*(tt247-4*&
&X(3,2)*X(3,3)+tt224)*tt80*tt140+2*tt1*tt80*tt169**2*tt128
hes(3,4) = tt248
hes(3,5) = tt249
hes(3,6) = tt250
hes(3,7) = tt251
hes(3,8) = tt252
hes(3,9) = tt253
hes(4,1) = tt196
hes(4,2) = tt241
hes(4,3) = tt248
hes(4,4) = tt1*(tt146+tt257+tt256+tt143-16*X(2,1)*X(2,3)+tt255)*t&
&t147-tt1*(tt254+tt167+tt166+(tt198+4*X(2,1))*X(2,3)+tt53+tt148+tt&
&141+(tt29-12*X(1,2))*X(1,3)+12*tt2-12*X(1,1)*X(1,2)+2*tt81)*tt80*&
&tt140-2*tt1*tt179*tt186*tt140+2*tt1*tt179**2*tt80*tt128
hes(4,5) = tt258
hes(4,6) = tt259
hes(4,7) = tt262
hes(4,8) = tt263
hes(4,9) = tt264
hes(5,1) = tt209
hes(5,2) = tt242
hes(5,3) = tt249
hes(5,4) = tt258
hes(5,5) = tt1*(tt146+tt257+tt256+2*tt37)*tt147-tt1*(tt254+tt167+&
&tt166+tt239+(tt58-12*X(2,2))*X(2,3)+12*tt17-12*X(2,1)*X(2,2)+2*tt&
&97)*tt80*tt140-2*tt1*tt203*tt208*tt140+2*tt1*tt203**2*tt80*tt128
hes(5,6) = tt265
hes(5,7) = tt266
hes(5,8) = tt267
hes(5,9) = tt268
hes(6,1) = tt214
hes(6,2) = tt243
hes(6,3) = tt250
hes(6,4) = tt259
hes(6,5) = tt265
hes(6,6) = 2*tt1*tt64*tt147-2*tt1*tt213*tt212*tt140-tt1*(tt247+(t&
&t192-12*X(3,2))*X(3,3)+12*tt20-12*X(3,1)*X(3,2)+2*tt118)*tt80*tt1&
&40+2*tt1*tt80*tt212**2*tt128
hes(6,7) = tt269
hes(6,8) = tt270
hes(6,9) = tt271
hes(7,1) = tt228
hes(7,2) = tt244
hes(7,3) = tt251
hes(7,4) = tt262
hes(7,5) = tt266
hes(7,6) = tt269
hes(7,7) = tt1*(tt144+tt272+tt256+tt142-16*X(2,1)*X(2,2)+tt255)*t&
&t147-2*tt1*tt88*tt80*tt140-2*tt1*tt223*tt219*tt140+2*tt1*tt219**2&
&*tt80*tt128
hes(7,8) = tt273
hes(7,9) = tt274
hes(8,1) = tt235
hes(8,2) = tt245
hes(8,3) = tt252
hes(8,4) = tt263
hes(8,5) = tt267
hes(8,6) = tt270
hes(8,7) = tt273
hes(8,8) = tt1*(tt144+tt272+tt256+2*tt51)*tt147-2*tt1*tt105*tt80*&
&tt140-2*tt1*tt234*tt231*tt140+2*tt1*tt231**2*tt80*tt128
hes(8,9) = tt275
hes(9,1) = tt238
hes(9,2) = tt246
hes(9,3) = tt253
hes(9,4) = tt264
hes(9,5) = tt268
hes(9,6) = tt271
hes(9,7) = tt274
hes(9,8) = tt275
hes(9,9) = 2*tt1*tt78*tt147-2*tt1*tt126*tt80*tt140-2*tt1*tt237*tt&
&236*tt140+2*tt1*tt236**2*tt80*tt128
END 
SUBROUTINE tri_linear_elas(val, X, D, area, lam, miu) 
IMPLICIT NONE 
REAL(KIND=8) val(1, 1) 
REAL(KIND=8) X(2, 3) 
REAL(KIND=8) D(2, 2) 
REAL(KIND=8) area(1, 1) 
REAL(KIND=8) lam(1, 1) 
REAL(KIND=8) miu(1, 1) 
REAL(KIND=8)  tt1 
REAL(KIND=8)  tt2 
tt1 = -X(1,1)
tt2 = -X(2,1)
val(1,1) = area(1,1)*(5.0E-1*lam(1,1)*(5.0E-1*(2*D(2,2)*(X(2,3)+t&
&t2)+2*D(1,2)*(X(2,2)+tt2))+5.0E-1*(2*(X(1,3)+tt1)*D(2,1)+2*D(1,1)&
&*(X(1,2)+tt1))-2)**2+miu(1,1)*((1.0E+0*D(2,2)*X(2,3)+1.0E+0*D(1,2&
&)*X(2,2)-1.0E+0*X(2,1)*D(2,2)-1.0E+0*D(1,2)*X(2,1)-1)**2+5.0E-1*(&
&D(2,1)*X(2,3)+D(1,1)*X(2,2)+X(1,3)*D(2,2)-X(1,1)*D(2,2)-D(2,1)*X(&
&2,1)-D(1,1)*X(2,1)+D(1,2)*X(1,2)-X(1,1)*D(1,2))**2+(1.0E+0*X(1,3)&
&*D(2,1)-1.0E+0*X(1,1)*D(2,1)+1.0E+0*D(1,1)*X(1,2)-1.0E+0*D(1,1)*X&
&(1,1)-1)**2))
END 
SUBROUTINE tri_linear_elas_jac(jac, X, D, area, lam, miu) 
IMPLICIT NONE 
REAL(KIND=8) jac(1, 6) 
REAL(KIND=8) X(2, 3) 
REAL(KIND=8) D(2, 2) 
REAL(KIND=8) area(1, 1) 
REAL(KIND=8) lam(1, 1) 
REAL(KIND=8) miu(1, 1) 
REAL(KIND=8)  tt1 
REAL(KIND=8)  tt2 
REAL(KIND=8)  tt3 
REAL(KIND=8)  tt4 
REAL(KIND=8)  tt5 
REAL(KIND=8)  tt6 
tt1 = 1.0E+0*X(1,3)*D(2,1)-1.0E+0*X(1,1)*D(2,1)+1.0E+0*D(1,1)*X(1&
&,2)-1.0E+0*D(1,1)*X(1,1)-1
tt2 = D(2,1)*X(2,3)+D(1,1)*X(2,2)+X(1,3)*D(2,2)-X(1,1)*D(2,2)-D(2&
&,1)*X(2,1)-D(1,1)*X(2,1)+D(1,2)*X(1,2)-X(1,1)*D(1,2)
tt3 = -X(1,1)
tt4 = -X(2,1)
tt5 = 5.0E-1*(2*D(2,2)*(X(2,3)+tt4)+2*D(1,2)*(X(2,2)+tt4))+5.0E-1&
&*(2*(X(1,3)+tt3)*D(2,1)+2*D(1,1)*(X(1,2)+tt3))-2
tt6 = 1.0E+0*D(2,2)*X(2,3)+1.0E+0*D(1,2)*X(2,2)-1.0E+0*X(2,1)*D(2&
&,2)-1.0E+0*D(1,2)*X(2,1)-1
jac(1,1) = area(1,1)*(5.0E-1*lam(1,1)*((-2*D(2,1))-2*D(1,1))*tt5+&
&miu(1,1)*(1.0E+0*((-D(2,2))-D(1,2))*tt2+2*((-1.0E+0*D(2,1))-1.0E+&
&0*D(1,1))*tt1))
jac(1,2) = area(1,1)*(5.0E-1*lam(1,1)*((-2*D(2,2))-2*D(1,2))*tt5+&
&miu(1,1)*(2*((-1.0E+0*D(2,2))-1.0E+0*D(1,2))*tt6+1.0E+0*((-D(2,1)&
&)-D(1,1))*tt2))
jac(1,3) = area(1,1)*(1.0E+0*D(1,1)*lam(1,1)*tt5+miu(1,1)*(1.0E+0&
&*D(1,2)*tt2+2.0E+0*D(1,1)*tt1))
jac(1,4) = area(1,1)*(1.0E+0*lam(1,1)*D(1,2)*tt5+miu(1,1)*(2.0E+0&
&*D(1,2)*tt6+1.0E+0*D(1,1)*tt2))
jac(1,5) = area(1,1)*(1.0E+0*lam(1,1)*D(2,1)*tt5+miu(1,1)*(1.0E+0&
&*D(2,2)*tt2+2.0E+0*D(2,1)*tt1))
jac(1,6) = area(1,1)*(1.0E+0*lam(1,1)*D(2,2)*tt5+miu(1,1)*(2.0E+0&
&*D(2,2)*tt6+1.0E+0*D(2,1)*tt2))
END 
SUBROUTINE tri_linear_elas_hes(hes, X, D, area, lam, miu) 
IMPLICIT NONE 
REAL(KIND=8) hes(6, 6) 
REAL(KIND=8) X(2, 3) 
REAL(KIND=8) D(2, 2) 
REAL(KIND=8) area(1, 1) 
REAL(KIND=8) lam(1, 1) 
REAL(KIND=8) miu(1, 1) 
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
tt1 = (-2*D(2,1))-2*D(1,1)
tt2 = (-1.0E+0*D(2,1))-1.0E+0*D(1,1)
tt3 = (-D(2,2))-D(1,2)
tt4 = (-2*D(2,2))-2*D(1,2)
tt5 = (-D(2,1))-D(1,1)
tt6 = area(1,1)*(1.0E+0*miu(1,1)*tt5*tt3+2.5E-1*lam(1,1)*tt1*tt4)&
&
tt7 = area(1,1)*(miu(1,1)*(1.0E+0*D(1,2)*tt3+2.0E+0*D(1,1)*tt2)+5&
&.0E-1*D(1,1)*lam(1,1)*tt1)
tt8 = area(1,1)*(1.0E+0*D(1,1)*miu(1,1)*tt3+5.0E-1*lam(1,1)*D(1,2&
&)*tt1)
tt9 = area(1,1)*(miu(1,1)*(1.0E+0*tt3*D(2,2)+2.0E+0*tt2*D(2,1))+5&
&.0E-1*lam(1,1)*tt1*D(2,1))
tt10 = area(1,1)*(5.0E-1*lam(1,1)*tt1*D(2,2)+1.0E+0*miu(1,1)*D(2,&
&1)*tt3)
tt11 = (-1.0E+0*D(2,2))-1.0E+0*D(1,2)
tt12 = area(1,1)*(5.0E-1*D(1,1)*lam(1,1)*tt4+1.0E+0*miu(1,1)*D(1,&
&2)*tt5)
tt13 = area(1,1)*(5.0E-1*lam(1,1)*D(1,2)*tt4+miu(1,1)*(2.0E+0*D(1&
&,2)*tt11+1.0E+0*D(1,1)*tt5))
tt14 = area(1,1)*(1.0E+0*miu(1,1)*tt5*D(2,2)+5.0E-1*lam(1,1)*D(2,&
&1)*tt4)
tt15 = area(1,1)*(miu(1,1)*(2.0E+0*tt11*D(2,2)+1.0E+0*tt5*D(2,1))&
&+5.0E-1*lam(1,1)*tt4*D(2,2))
tt16 = D(1,1)**2
tt17 = D(1,2)**2
tt18 = area(1,1)*(1.0E+0*D(1,1)*miu(1,1)*D(1,2)+1.0E+0*D(1,1)*lam&
&(1,1)*D(1,2))
tt19 = area(1,1)*(miu(1,1)*(1.0E+0*D(1,2)*D(2,2)+2.0E+0*D(1,1)*D(&
&2,1))+1.0E+0*D(1,1)*lam(1,1)*D(2,1))
tt20 = area(1,1)*(1.0E+0*D(1,1)*lam(1,1)*D(2,2)+1.0E+0*miu(1,1)*D&
&(1,2)*D(2,1))
tt21 = area(1,1)*(1.0E+0*D(1,1)*miu(1,1)*D(2,2)+1.0E+0*lam(1,1)*D&
&(1,2)*D(2,1))
tt22 = area(1,1)*(miu(1,1)*(2.0E+0*D(1,2)*D(2,2)+1.0E+0*D(1,1)*D(&
&2,1))+1.0E+0*lam(1,1)*D(1,2)*D(2,2))
tt23 = D(2,1)**2
tt24 = D(2,2)**2
tt25 = area(1,1)*(1.0E+0*miu(1,1)*D(2,1)*D(2,2)+1.0E+0*lam(1,1)*D&
&(2,1)*D(2,2))
hes(1,1) = area(1,1)*(miu(1,1)*(1.0E+0*tt3**2+2*tt2**2)+2.5E-1*la&
&m(1,1)*tt1**2)
hes(1,2) = tt6
hes(1,3) = tt7
hes(1,4) = tt8
hes(1,5) = tt9
hes(1,6) = tt10
hes(2,1) = tt6
hes(2,2) = area(1,1)*(miu(1,1)*(2*tt11**2+1.0E+0*tt5**2)+2.5E-1*l&
&am(1,1)*tt4**2)
hes(2,3) = tt12
hes(2,4) = tt13
hes(2,5) = tt14
hes(2,6) = tt15
hes(3,1) = tt7
hes(3,2) = tt12
hes(3,3) = area(1,1)*(miu(1,1)*(1.0E+0*tt17+2.0E+0*tt16)+1.0E+0*t&
&t16*lam(1,1))
hes(3,4) = tt18
hes(3,5) = tt19
hes(3,6) = tt20
hes(4,1) = tt8
hes(4,2) = tt13
hes(4,3) = tt18
hes(4,4) = area(1,1)*(miu(1,1)*(2.0E+0*tt17+1.0E+0*tt16)+1.0E+0*l&
&am(1,1)*tt17)
hes(4,5) = tt21
hes(4,6) = tt22
hes(5,1) = tt9
hes(5,2) = tt14
hes(5,3) = tt19
hes(5,4) = tt21
hes(5,5) = area(1,1)*(miu(1,1)*(1.0E+0*tt24+2.0E+0*tt23)+1.0E+0*l&
&am(1,1)*tt23)
hes(5,6) = tt25
hes(6,1) = tt10
hes(6,2) = tt15
hes(6,3) = tt20
hes(6,4) = tt22
hes(6,5) = tt25
hes(6,6) = area(1,1)*(miu(1,1)*(2.0E+0*tt24+1.0E+0*tt23)+1.0E+0*l&
&am(1,1)*tt24)
END 
SUBROUTINE tri_stvk_elas(val, X, D, area, lam, miu) 
IMPLICIT NONE 
REAL(KIND=8) val(1, 1) 
REAL(KIND=8) X(2, 3) 
REAL(KIND=8) D(2, 2) 
REAL(KIND=8) area(1, 1) 
REAL(KIND=8) lam(1, 1) 
REAL(KIND=8) miu(1, 1) 
REAL(KIND=8)  tt1 
REAL(KIND=8)  tt2 
REAL(KIND=8)  tt3 
REAL(KIND=8)  tt4 
REAL(KIND=8)  tt5 
REAL(KIND=8)  tt6 
REAL(KIND=8)  tt7 
REAL(KIND=8)  tt8 
tt1 = X(1,1)**2
tt2 = X(2,1)**2
tt3 = -X(1,1)
tt4 = X(1,2)+tt3
tt5 = X(1,3)+tt3
tt6 = -X(2,1)
tt7 = X(2,2)+tt6
tt8 = X(2,3)+tt6
val(1,1) = area(1,1)*(5.0E-1*lam(1,1)*(((D(2,2)*tt8+D(1,2)*tt7)**&
&2+(tt5*D(2,2)+D(1,2)*tt4)**2-1)/2.0E+0+((D(2,1)*tt8+D(1,1)*tt7)**&
&2+(tt5*D(2,1)+D(1,1)*tt4)**2-1)/2.0E+0)**2+miu(1,1)*(((D(2,2)*X(2&
&,3)+D(1,2)*X(2,2)-X(2,1)*D(2,2)-D(1,2)*X(2,1))**2+(X(1,3)*D(2,2)-&
&X(1,1)*D(2,2)+D(1,2)*X(1,2)-X(1,1)*D(1,2))**2-1)**2/4.0E+0+((D(2,&
&1)*X(2,3)+D(1,1)*X(2,2)-D(2,1)*X(2,1)-D(1,1)*X(2,1))**2+(X(1,3)*D&
&(2,1)-X(1,1)*D(2,1)+D(1,1)*X(1,2)-D(1,1)*X(1,1))**2-1)**2/4.0E+0+&
&(D(2,1)*D(2,2)*X(2,3)**2+D(1,1)*D(2,2)*X(2,2)*X(2,3)+D(1,2)*D(2,1&
&)*X(2,2)*X(2,3)-2*D(2,1)*X(2,1)*D(2,2)*X(2,3)-D(1,1)*X(2,1)*D(2,2&
&)*X(2,3)-D(1,2)*D(2,1)*X(2,1)*X(2,3)+D(1,1)*D(1,2)*X(2,2)**2-D(1,&
&1)*X(2,1)*D(2,2)*X(2,2)-D(1,2)*D(2,1)*X(2,1)*X(2,2)-2*D(1,1)*D(1,&
&2)*X(2,1)*X(2,2)+D(2,1)*tt2*D(2,2)+D(1,1)*tt2*D(2,2)+X(1,3)**2*D(&
&2,1)*D(2,2)-2*X(1,1)*X(1,3)*D(2,1)*D(2,2)+tt1*D(2,1)*D(2,2)+D(1,1&
&)*X(1,2)*X(1,3)*D(2,2)-D(1,1)*X(1,1)*X(1,3)*D(2,2)-D(1,1)*X(1,1)*&
&X(1,2)*D(2,2)+D(1,1)*tt1*D(2,2)+D(1,2)*D(2,1)*tt2+D(1,1)*D(1,2)*t&
&t2+D(1,2)*X(1,2)*X(1,3)*D(2,1)-X(1,1)*D(1,2)*X(1,3)*D(2,1)-X(1,1)&
&*D(1,2)*X(1,2)*D(2,1)+tt1*D(1,2)*D(2,1)+D(1,1)*D(1,2)*X(1,2)**2-2&
&*D(1,1)*X(1,1)*D(1,2)*X(1,2)+D(1,1)*tt1*D(1,2))**2/2.0E+0))
END 
SUBROUTINE tri_stvk_elas_jac(jac, X, D, area, lam, miu) 
IMPLICIT NONE 
REAL(KIND=8) jac(1, 6) 
REAL(KIND=8) X(2, 3) 
REAL(KIND=8) D(2, 2) 
REAL(KIND=8) area(1, 1) 
REAL(KIND=8) lam(1, 1) 
REAL(KIND=8) miu(1, 1) 
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
tt1 = X(1,1)**2
tt2 = X(2,1)**2
tt3 = D(2,1)*D(2,2)*X(2,3)**2+D(1,1)*D(2,2)*X(2,2)*X(2,3)+D(1,2)*&
&D(2,1)*X(2,2)*X(2,3)-2*D(2,1)*X(2,1)*D(2,2)*X(2,3)-D(1,1)*X(2,1)*&
&D(2,2)*X(2,3)-D(1,2)*D(2,1)*X(2,1)*X(2,3)+D(1,1)*D(1,2)*X(2,2)**2&
&-D(1,1)*X(2,1)*D(2,2)*X(2,2)-D(1,2)*D(2,1)*X(2,1)*X(2,2)-2*D(1,1)&
&*D(1,2)*X(2,1)*X(2,2)+D(2,1)*tt2*D(2,2)+D(1,1)*tt2*D(2,2)+X(1,3)*&
&*2*D(2,1)*D(2,2)-2*X(1,1)*X(1,3)*D(2,1)*D(2,2)+tt1*D(2,1)*D(2,2)+&
&D(1,1)*X(1,2)*X(1,3)*D(2,2)-D(1,1)*X(1,1)*X(1,3)*D(2,2)-D(1,1)*X(&
&1,1)*X(1,2)*D(2,2)+D(1,1)*tt1*D(2,2)+D(1,2)*D(2,1)*tt2+D(1,1)*D(1&
&,2)*tt2+D(1,2)*X(1,2)*X(1,3)*D(2,1)-X(1,1)*D(1,2)*X(1,3)*D(2,1)-X&
&(1,1)*D(1,2)*X(1,2)*D(2,1)+tt1*D(1,2)*D(2,1)+D(1,1)*D(1,2)*X(1,2)&
&**2-2*D(1,1)*X(1,1)*D(1,2)*X(1,2)+D(1,1)*tt1*D(1,2)
tt4 = (-D(2,1))-D(1,1)
tt5 = X(1,3)*D(2,1)-X(1,1)*D(2,1)+D(1,1)*X(1,2)-D(1,1)*X(1,1)
tt6 = D(2,1)*X(2,3)+D(1,1)*X(2,2)-D(2,1)*X(2,1)-D(1,1)*X(2,1)
tt7 = tt6**2+tt5**2-1
tt8 = (-D(2,2))-D(1,2)
tt9 = X(1,3)*D(2,2)-X(1,1)*D(2,2)+D(1,2)*X(1,2)-X(1,1)*D(1,2)
tt10 = D(2,2)*X(2,3)+D(1,2)*X(2,2)-X(2,1)*D(2,2)-D(1,2)*X(2,1)
tt11 = tt10**2+tt9**2-1
tt12 = -X(1,1)
tt13 = X(1,2)+tt12
tt14 = X(1,3)+tt12
tt15 = tt14*D(2,1)+D(1,1)*tt13
tt16 = tt14*D(2,2)+D(1,2)*tt13
tt17 = -X(2,1)
tt18 = X(2,2)+tt17
tt19 = X(2,3)+tt17
tt20 = D(2,1)*tt19+D(1,1)*tt18
tt21 = D(2,2)*tt19+D(1,2)*tt18
tt22 = (tt21**2+tt16**2-1)/2.0E+0+(tt20**2+tt15**2-1)/2.0E+0
tt23 = -X(1,1)*D(1,2)*D(2,1)
tt24 = -D(1,1)*X(1,1)*D(2,2)
tt25 = -D(1,2)*D(2,1)*X(2,1)
tt26 = -D(1,1)*X(2,1)*D(2,2)
jac(1,1) = area(1,1)*(1.0E+0*lam(1,1)*(tt8*tt16+tt4*tt15)*tt22+mi&
&u(1,1)*(tt8*tt9*tt11+tt4*tt5*tt7+((-2*X(1,3)*D(2,1)*D(2,2))+2*X(1&
&,1)*D(2,1)*D(2,2)-D(1,1)*X(1,3)*D(2,2)-D(1,1)*X(1,2)*D(2,2)+2*D(1&
&,1)*X(1,1)*D(2,2)-D(1,2)*X(1,3)*D(2,1)-D(1,2)*X(1,2)*D(2,1)+2*X(1&
&,1)*D(1,2)*D(2,1)-2*D(1,1)*D(1,2)*X(1,2)+2*D(1,1)*X(1,1)*D(1,2))*&
&tt3))
jac(1,2) = area(1,1)*(1.0E+0*lam(1,1)*(tt8*tt21+tt4*tt20)*tt22+mi&
&u(1,1)*(tt8*tt10*tt11+tt4*tt6*tt7+((-2*D(2,1)*D(2,2)*X(2,3))-D(1,&
&1)*D(2,2)*X(2,3)-D(1,2)*D(2,1)*X(2,3)-D(1,1)*D(2,2)*X(2,2)-D(1,2)&
&*D(2,1)*X(2,2)-2*D(1,1)*D(1,2)*X(2,2)+2*D(2,1)*X(2,1)*D(2,2)+2*D(&
&1,1)*X(2,1)*D(2,2)+2*D(1,2)*D(2,1)*X(2,1)+2*D(1,1)*D(1,2)*X(2,1))&
&*tt3))
jac(1,3) = area(1,1)*(1.0E+0*lam(1,1)*(D(1,2)*tt16+D(1,1)*tt15)*t&
&t22+miu(1,1)*(D(1,2)*tt9*tt11+D(1,1)*tt5*tt7+(D(1,1)*X(1,3)*D(2,2&
&)+tt24+D(1,2)*X(1,3)*D(2,1)+tt23+2*D(1,1)*D(1,2)*X(1,2)-2*D(1,1)*&
&X(1,1)*D(1,2))*tt3))
jac(1,4) = area(1,1)*(1.0E+0*lam(1,1)*(D(1,2)*tt21+D(1,1)*tt20)*t&
&t22+miu(1,1)*(D(1,2)*tt10*tt11+D(1,1)*tt6*tt7+(D(1,1)*D(2,2)*X(2,&
&3)+D(1,2)*D(2,1)*X(2,3)+2*D(1,1)*D(1,2)*X(2,2)+tt26+tt25-2*D(1,1)&
&*D(1,2)*X(2,1))*tt3))
jac(1,5) = area(1,1)*(1.0E+0*lam(1,1)*(D(2,2)*tt16+D(2,1)*tt15)*t&
&t22+miu(1,1)*(D(2,2)*tt9*tt11+D(2,1)*tt5*tt7+(2*X(1,3)*D(2,1)*D(2&
&,2)-2*X(1,1)*D(2,1)*D(2,2)+D(1,1)*X(1,2)*D(2,2)+tt24+D(1,2)*X(1,2&
&)*D(2,1)+tt23)*tt3))
jac(1,6) = area(1,1)*(1.0E+0*lam(1,1)*(D(2,2)*tt21+D(2,1)*tt20)*t&
&t22+miu(1,1)*(D(2,2)*tt10*tt11+D(2,1)*tt6*tt7+(2*D(2,1)*D(2,2)*X(&
&2,3)+D(1,1)*D(2,2)*X(2,2)+D(1,2)*D(2,1)*X(2,2)-2*D(2,1)*X(2,1)*D(&
&2,2)+tt26+tt25)*tt3))
END 
SUBROUTINE tri_stvk_elas_hes(hes, X, D, area, lam, miu) 
IMPLICIT NONE 
REAL(KIND=8) hes(6, 6) 
REAL(KIND=8) X(2, 3) 
REAL(KIND=8) D(2, 2) 
REAL(KIND=8) area(1, 1) 
REAL(KIND=8) lam(1, 1) 
REAL(KIND=8) miu(1, 1) 
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
tt1 = (-D(2,1))-D(1,1)
tt2 = -X(1,1)
tt3 = X(1,2)+tt2
tt4 = X(1,3)+tt2
tt5 = tt4*D(2,1)+D(1,1)*tt3
tt6 = (-D(2,2))-D(1,2)
tt7 = tt4*D(2,2)+D(1,2)*tt3
tt8 = tt6*tt7+tt1*tt5
tt9 = tt1**2
tt10 = X(1,3)*D(2,1)-X(1,1)*D(2,1)+D(1,1)*X(1,2)-D(1,1)*X(1,1)
tt11 = tt10**2
tt12 = tt6**2
tt13 = X(1,3)*D(2,2)-X(1,1)*D(2,2)+D(1,2)*X(1,2)-X(1,1)*D(1,2)
tt14 = tt13**2
tt15 = (-2*X(1,3)*D(2,1)*D(2,2))+2*X(1,1)*D(2,1)*D(2,2)-D(1,1)*X(&
&1,3)*D(2,2)-D(1,1)*X(1,2)*D(2,2)+2*D(1,1)*X(1,1)*D(2,2)-D(1,2)*X(&
&1,3)*D(2,1)-D(1,2)*X(1,2)*D(2,1)+2*X(1,1)*D(1,2)*D(2,1)-2*D(1,1)*&
&D(1,2)*X(1,2)+2*D(1,1)*X(1,1)*D(1,2)
tt16 = X(1,1)**2
tt17 = X(2,1)**2
tt18 = D(2,1)*D(2,2)*X(2,3)**2+D(1,1)*D(2,2)*X(2,2)*X(2,3)+D(1,2)&
&*D(2,1)*X(2,2)*X(2,3)-2*D(2,1)*X(2,1)*D(2,2)*X(2,3)-D(1,1)*X(2,1)&
&*D(2,2)*X(2,3)-D(1,2)*D(2,1)*X(2,1)*X(2,3)+D(1,1)*D(1,2)*X(2,2)**&
&2-D(1,1)*X(2,1)*D(2,2)*X(2,2)-D(1,2)*D(2,1)*X(2,1)*X(2,2)-2*D(1,1&
&)*D(1,2)*X(2,1)*X(2,2)+D(2,1)*tt17*D(2,2)+D(1,1)*tt17*D(2,2)+X(1,&
&3)**2*D(2,1)*D(2,2)-2*X(1,1)*X(1,3)*D(2,1)*D(2,2)+tt16*D(2,1)*D(2&
&,2)+D(1,1)*X(1,2)*X(1,3)*D(2,2)-D(1,1)*X(1,1)*X(1,3)*D(2,2)-D(1,1&
&)*X(1,1)*X(1,2)*D(2,2)+D(1,1)*tt16*D(2,2)+D(1,2)*D(2,1)*tt17+D(1,&
&1)*D(1,2)*tt17+D(1,2)*X(1,2)*X(1,3)*D(2,1)-X(1,1)*D(1,2)*X(1,3)*D&
&(2,1)-X(1,1)*D(1,2)*X(1,2)*D(2,1)+tt16*D(1,2)*D(2,1)+D(1,1)*D(1,2&
&)*X(1,2)**2-2*D(1,1)*X(1,1)*D(1,2)*X(1,2)+D(1,1)*tt16*D(1,2)
tt19 = (2*D(2,1)*D(2,2)+2*D(1,1)*D(2,2)+2*D(1,2)*D(2,1)+2*D(1,1)*&
&D(1,2))*tt18
tt20 = D(2,1)*X(2,3)+D(1,1)*X(2,2)-D(2,1)*X(2,1)-D(1,1)*X(2,1)
tt21 = tt20**2
tt22 = tt21+tt11-1
tt23 = tt9*tt22
tt24 = D(2,2)*X(2,3)+D(1,2)*X(2,2)-X(2,1)*D(2,2)-D(1,2)*X(2,1)
tt25 = tt24**2
tt26 = tt25+tt14-1
tt27 = tt12*tt26
tt28 = -X(2,1)
tt29 = X(2,2)+tt28
tt30 = X(2,3)+tt28
tt31 = D(2,1)*tt30+D(1,1)*tt29
tt32 = D(2,2)*tt30+D(1,2)*tt29
tt33 = (tt32**2+tt7**2-1)/2.0E+0+(tt31**2+tt5**2-1)/2.0E+0
tt34 = 1.0E+0*lam(1,1)*(tt12+tt9)*tt33
tt35 = (-2*D(2,1)*D(2,2)*X(2,3))-D(1,1)*D(2,2)*X(2,3)-D(1,2)*D(2,&
&1)*X(2,3)-D(1,1)*D(2,2)*X(2,2)-D(1,2)*D(2,1)*X(2,2)-2*D(1,1)*D(1,&
&2)*X(2,2)+2*D(2,1)*X(2,1)*D(2,2)+2*D(1,1)*X(2,1)*D(2,2)+2*D(1,2)*&
&D(2,1)*X(2,1)+2*D(1,1)*D(1,2)*X(2,1)
tt36 = tt6*tt32+tt1*tt31
tt37 = area(1,1)*(1.0E+0*lam(1,1)*tt8*tt36+miu(1,1)*(tt15*tt35+2*&
&tt12*tt13*tt24+2*tt9*tt10*tt20))
tt38 = D(1,2)*tt7+D(1,1)*tt5
tt39 = -X(1,1)*D(1,2)*D(2,1)
tt40 = -D(1,1)*X(1,1)*D(2,2)
tt41 = D(1,1)*X(1,3)*D(2,2)+tt40+D(1,2)*X(1,3)*D(2,1)+tt39+2*D(1,&
&1)*D(1,2)*X(1,2)-2*D(1,1)*X(1,1)*D(1,2)
tt42 = -D(1,2)*D(2,1)
tt43 = -D(1,1)*D(2,2)
tt44 = (tt43+tt42-2*D(1,1)*D(1,2))*tt18
tt45 = D(1,1)*tt1*tt22
tt46 = D(1,2)*tt6*tt26
tt47 = 1.0E+0*lam(1,1)*(D(1,2)*tt6+D(1,1)*tt1)*tt33
tt48 = area(1,1)*(tt47+miu(1,1)*(tt46+tt45+tt44+2*D(1,2)*tt6*tt14&
&+tt41*tt15+2*D(1,1)*tt1*tt11)+1.0E+0*lam(1,1)*tt38*tt8)
tt49 = 2*D(1,1)*tt1*tt10*tt20
tt50 = 2*D(1,2)*tt6*tt13*tt24
tt51 = -D(1,2)*D(2,1)*X(2,1)
tt52 = -D(1,1)*X(2,1)*D(2,2)
tt53 = D(1,1)*D(2,2)*X(2,3)+D(1,2)*D(2,1)*X(2,3)+2*D(1,1)*D(1,2)*&
&X(2,2)+tt52+tt51-2*D(1,1)*D(1,2)*X(2,1)
tt54 = D(1,2)*tt32+D(1,1)*tt31
tt55 = area(1,1)*(1.0E+0*lam(1,1)*tt8*tt54+miu(1,1)*(tt15*tt53+tt&
&50+tt49))
tt56 = D(2,2)*tt7+D(2,1)*tt5
tt57 = 2*X(1,3)*D(2,1)*D(2,2)-2*X(1,1)*D(2,1)*D(2,2)+D(1,1)*X(1,2&
&)*D(2,2)+tt40+D(1,2)*X(1,2)*D(2,1)+tt39
tt58 = ((-2*D(2,1)*D(2,2))+tt43+tt42)*tt18
tt59 = tt1*D(2,1)*tt22
tt60 = tt6*D(2,2)*tt26
tt61 = 1.0E+0*lam(1,1)*(tt6*D(2,2)+tt1*D(2,1))*tt33
tt62 = area(1,1)*(tt61+miu(1,1)*(tt60+tt59+tt58+2*tt6*D(2,2)*tt14&
&+tt15*tt57+2*tt1*D(2,1)*tt11)+1.0E+0*lam(1,1)*tt8*tt56)
tt63 = 2*tt1*D(2,1)*tt10*tt20
tt64 = 2*tt6*D(2,2)*tt13*tt24
tt65 = 2*D(2,1)*D(2,2)*X(2,3)+D(1,1)*D(2,2)*X(2,2)+D(1,2)*D(2,1)*&
&X(2,2)-2*D(2,1)*X(2,1)*D(2,2)+tt52+tt51
tt66 = D(2,2)*tt32+D(2,1)*tt31
tt67 = area(1,1)*(1.0E+0*lam(1,1)*tt8*tt66+miu(1,1)*(tt15*tt65+tt&
&64+tt63))
tt68 = area(1,1)*(1.0E+0*lam(1,1)*tt38*tt36+miu(1,1)*(tt41*tt35+t&
&t50+tt49))
tt69 = area(1,1)*(tt47+miu(1,1)*(tt46+2*D(1,2)*tt6*tt25+tt45+2*D(&
&1,1)*tt1*tt21+tt44+tt53*tt35)+1.0E+0*lam(1,1)*tt54*tt36)
tt70 = area(1,1)*(1.0E+0*lam(1,1)*tt56*tt36+miu(1,1)*(tt57*tt35+t&
&t64+tt63))
tt71 = area(1,1)*(tt61+miu(1,1)*(tt60+2*tt6*D(2,2)*tt25+tt59+2*tt&
&1*D(2,1)*tt21+tt58+tt35*tt65)+1.0E+0*lam(1,1)*tt36*tt66)
tt72 = D(1,1)**2
tt73 = D(1,2)**2
tt74 = 2*D(1,1)*D(1,2)*tt18
tt75 = tt72*tt22
tt76 = tt73*tt26
tt77 = 1.0E+0*lam(1,1)*(tt73+tt72)*tt33
tt78 = area(1,1)*(1.0E+0*lam(1,1)*tt38*tt54+miu(1,1)*(tt41*tt53+2&
&*tt73*tt13*tt24+2*tt72*tt10*tt20))
tt79 = (D(1,1)*D(2,2)+D(1,2)*D(2,1))*tt18
tt80 = D(1,1)*D(2,1)*tt22
tt81 = D(1,2)*D(2,2)*tt26
tt82 = 1.0E+0*lam(1,1)*(D(1,2)*D(2,2)+D(1,1)*D(2,1))*tt33
tt83 = area(1,1)*(tt82+miu(1,1)*(tt81+tt80+tt79+2*D(1,2)*D(2,2)*t&
&t14+tt41*tt57+2*D(1,1)*D(2,1)*tt11)+1.0E+0*lam(1,1)*tt38*tt56)
tt84 = 2*D(1,1)*D(2,1)*tt10*tt20
tt85 = 2*D(1,2)*D(2,2)*tt13*tt24
tt86 = area(1,1)*(1.0E+0*lam(1,1)*tt38*tt66+miu(1,1)*(tt41*tt65+t&
&t85+tt84))
tt87 = area(1,1)*(1.0E+0*lam(1,1)*tt56*tt54+miu(1,1)*(tt57*tt53+t&
&t85+tt84))
tt88 = area(1,1)*(tt82+miu(1,1)*(tt81+2*D(1,2)*D(2,2)*tt25+tt80+2&
&*D(1,1)*D(2,1)*tt21+tt79+tt53*tt65)+1.0E+0*lam(1,1)*tt54*tt66)
tt89 = D(2,1)**2
tt90 = D(2,2)**2
tt91 = 2*D(2,1)*D(2,2)*tt18
tt92 = tt89*tt22
tt93 = tt90*tt26
tt94 = 1.0E+0*lam(1,1)*(tt90+tt89)*tt33
tt95 = area(1,1)*(1.0E+0*lam(1,1)*tt56*tt66+miu(1,1)*(tt57*tt65+2&
&*tt90*tt13*tt24+2*tt89*tt10*tt20))
hes(1,1) = area(1,1)*(tt34+miu(1,1)*(tt27+tt23+tt19+tt15**2+2*tt1&
&2*tt14+2*tt9*tt11)+1.0E+0*lam(1,1)*tt8**2)
hes(1,2) = tt37
hes(1,3) = tt48
hes(1,4) = tt55
hes(1,5) = tt62
hes(1,6) = tt67
hes(2,1) = tt37
hes(2,2) = area(1,1)*(1.0E+0*lam(1,1)*tt36**2+tt34+miu(1,1)*(tt35&
&**2+tt27+2*tt12*tt25+tt23+2*tt9*tt21+tt19))
hes(2,3) = tt68
hes(2,4) = tt69
hes(2,5) = tt70
hes(2,6) = tt71
hes(3,1) = tt48
hes(3,2) = tt68
hes(3,3) = area(1,1)*(tt77+miu(1,1)*(tt76+tt75+tt74+tt41**2+2*tt7&
&3*tt14+2*tt72*tt11)+1.0E+0*lam(1,1)*tt38**2)
hes(3,4) = tt78
hes(3,5) = tt83
hes(3,6) = tt86
hes(4,1) = tt55
hes(4,2) = tt69
hes(4,3) = tt78
hes(4,4) = area(1,1)*(1.0E+0*lam(1,1)*tt54**2+tt77+miu(1,1)*(tt53&
&**2+tt76+2*tt73*tt25+tt75+2*tt72*tt21+tt74))
hes(4,5) = tt87
hes(4,6) = tt88
hes(5,1) = tt62
hes(5,2) = tt70
hes(5,3) = tt83
hes(5,4) = tt87
hes(5,5) = area(1,1)*(tt94+miu(1,1)*(tt93+tt92+tt91+tt57**2+2*tt9&
&0*tt14+2*tt89*tt11)+1.0E+0*lam(1,1)*tt56**2)
hes(5,6) = tt95
hes(6,1) = tt67
hes(6,2) = tt71
hes(6,3) = tt86
hes(6,4) = tt88
hes(6,5) = tt95
hes(6,6) = area(1,1)*(1.0E+0*lam(1,1)*tt66**2+tt94+miu(1,1)*(tt65&
&**2+tt93+2*tt90*tt25+tt92+2*tt89*tt21+tt91))
END 
SUBROUTINE tet_stvk(val, X, D, volume, lam, miu) 
IMPLICIT NONE 
REAL(KIND=8) val(1, 1) 
REAL(KIND=8) X(3, 4) 
REAL(KIND=8) D(3, 3) 
REAL(KIND=8) volume(1, 1) 
REAL(KIND=8) lam(1, 1) 
REAL(KIND=8) miu(1, 1) 
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
tt1 = X(1,1)**2
tt2 = X(1,2)**2
tt3 = X(2,1)**2
tt4 = X(1,3)**2
tt5 = X(2,2)**2
tt6 = X(2,3)**2
tt7 = X(3,1)**2
tt8 = X(1,4)**2
tt9 = X(2,4)**2
tt10 = X(3,2)**2
tt11 = X(3,3)**2
tt12 = X(3,4)**2
tt13 = -X(1,1)
tt14 = X(1,2)+tt13
tt15 = X(1,3)+tt13
tt16 = X(1,4)+tt13
tt17 = -X(2,1)
tt18 = X(2,2)+tt17
tt19 = X(2,3)+tt17
tt20 = X(2,4)+tt17
tt21 = -X(3,1)
tt22 = X(3,2)+tt21
tt23 = X(3,3)+tt21
tt24 = X(3,4)+tt21
val(1,1) = volume(1,1)*(5.0E-1*lam(1,1)*(5.0E-1*((D(3,3)*tt24+D(2&
&,3)*tt23+D(1,3)*tt22)**2+(tt20*D(3,3)+D(2,3)*tt19+D(1,3)*tt18)**2&
&+(tt16*D(3,3)+tt15*D(2,3)+tt14*D(1,3))**2-1)+5.0E-1*((D(3,2)*tt24&
&+D(2,2)*tt23+D(1,2)*tt22)**2+(tt20*D(3,2)+D(2,2)*tt19+D(1,2)*tt18&
&)**2+(tt16*D(3,2)+tt15*D(2,2)+D(1,2)*tt14)**2-1)+5.0E-1*((D(3,1)*&
&tt24+D(2,1)*tt23+D(1,1)*tt22)**2+(tt20*D(3,1)+D(2,1)*tt19+D(1,1)*&
&tt18)**2+(tt16*D(3,1)+tt15*D(2,1)+D(1,1)*tt14)**2-1))**2+miu(1,1)&
&*(2.5E-1*((D(3,3)*X(3,4)+D(2,3)*X(3,3)-X(3,1)*D(3,3)+D(1,3)*X(3,2&
&)-D(2,3)*X(3,1)-D(1,3)*X(3,1))**2+(X(2,4)*D(3,3)-X(2,1)*D(3,3)+D(&
&2,3)*X(2,3)-X(2,1)*D(2,3)+D(1,3)*X(2,2)-D(1,3)*X(2,1))**2+(X(1,4)&
&*D(3,3)-X(1,1)*D(3,3)+X(1,3)*D(2,3)-X(1,1)*D(2,3)+X(1,2)*D(1,3)-X&
&(1,1)*D(1,3))**2-1)**2+2.5E-1*((D(3,2)*X(3,4)+D(2,2)*X(3,3)+D(1,2&
&)*X(3,2)-X(3,1)*D(3,2)-D(2,2)*X(3,1)-D(1,2)*X(3,1))**2+(X(2,4)*D(&
&3,2)-X(2,1)*D(3,2)+D(2,2)*X(2,3)+D(1,2)*X(2,2)-X(2,1)*D(2,2)-D(1,&
&2)*X(2,1))**2+(X(1,4)*D(3,2)-X(1,1)*D(3,2)+X(1,3)*D(2,2)-X(1,1)*D&
&(2,2)+D(1,2)*X(1,2)-X(1,1)*D(1,2))**2-1)**2+2.5E-1*((D(3,1)*X(3,4&
&)+D(2,1)*X(3,3)+D(1,1)*X(3,2)-D(3,1)*X(3,1)-D(2,1)*X(3,1)-D(1,1)*&
&X(3,1))**2+(X(2,4)*D(3,1)-X(2,1)*D(3,1)+D(2,1)*X(2,3)+D(1,1)*X(2,&
&2)-D(2,1)*X(2,1)-D(1,1)*X(2,1))**2+(X(1,4)*D(3,1)-X(1,1)*D(3,1)+X&
&(1,3)*D(2,1)-X(1,1)*D(2,1)+D(1,1)*X(1,2)-D(1,1)*X(1,1))**2-1)**2+&
&5.0E-1*(D(3,2)*D(3,3)*tt12+D(2,2)*D(3,3)*X(3,3)*X(3,4)+D(2,3)*D(3&
&,2)*X(3,3)*X(3,4)+D(1,2)*X(3,2)*D(3,3)*X(3,4)-2*X(3,1)*D(3,2)*D(3&
&,3)*X(3,4)-D(2,2)*X(3,1)*D(3,3)*X(3,4)-D(1,2)*X(3,1)*D(3,3)*X(3,4&
&)+D(1,3)*D(3,2)*X(3,2)*X(3,4)-D(2,3)*X(3,1)*D(3,2)*X(3,4)-D(1,3)*&
&X(3,1)*D(3,2)*X(3,4)+D(2,2)*D(2,3)*tt11-D(2,2)*X(3,1)*D(3,3)*X(3,&
&3)+D(1,2)*D(2,3)*X(3,2)*X(3,3)+D(1,3)*D(2,2)*X(3,2)*X(3,3)-D(2,3)&
&*X(3,1)*D(3,2)*X(3,3)-2*D(2,2)*D(2,3)*X(3,1)*X(3,3)-D(1,2)*D(2,3)&
&*X(3,1)*X(3,3)-D(1,3)*D(2,2)*X(3,1)*X(3,3)-D(1,2)*X(3,1)*X(3,2)*D&
&(3,3)+tt7*D(3,2)*D(3,3)+tt9*D(3,2)*D(3,3)-2*X(2,1)*X(2,4)*D(3,2)*&
&D(3,3)+tt3*D(3,2)*D(3,3)+tt8*D(3,2)*D(3,3)-2*X(1,1)*X(1,4)*D(3,2)&
&*D(3,3)+tt1*D(3,2)*D(3,3)+D(2,2)*tt7*D(3,3)+D(1,2)*tt7*D(3,3)+D(2&
&,2)*X(2,3)*X(2,4)*D(3,3)+D(1,2)*X(2,2)*X(2,4)*D(3,3)-X(2,1)*D(2,2&
&)*X(2,4)*D(3,3)-D(1,2)*X(2,1)*X(2,4)*D(3,3)-X(2,1)*D(2,2)*X(2,3)*&
&D(3,3)-D(1,2)*X(2,1)*X(2,2)*D(3,3)+tt3*D(2,2)*D(3,3)+X(1,3)*X(1,4&
&)*D(2,2)*D(3,3)-X(1,1)*X(1,4)*D(2,2)*D(3,3)-X(1,1)*X(1,3)*D(2,2)*&
&D(3,3)+tt1*D(2,2)*D(3,3)+D(1,2)*tt3*D(3,3)+D(1,2)*X(1,2)*X(1,4)*D&
&(3,3)-X(1,1)*D(1,2)*X(1,4)*D(3,3)-X(1,1)*D(1,2)*X(1,2)*D(3,3)+tt1&
&*D(1,2)*D(3,3)+D(1,2)*D(1,3)*tt10-D(1,3)*X(3,1)*D(3,2)*X(3,2)-D(1&
&,2)*D(2,3)*X(3,1)*X(3,2)-D(1,3)*D(2,2)*X(3,1)*X(3,2)-2*D(1,2)*D(1&
&,3)*X(3,1)*X(3,2)+D(2,3)*tt7*D(3,2)+D(1,3)*tt7*D(3,2)+D(2,3)*X(2,&
&3)*X(2,4)*D(3,2)-X(2,1)*D(2,3)*X(2,4)*D(3,2)+D(1,3)*X(2,2)*X(2,4)&
&*D(3,2)-D(1,3)*X(2,1)*X(2,4)*D(3,2)-X(2,1)*D(2,3)*X(2,3)*D(3,2)+t&
&t3*D(2,3)*D(3,2)+X(1,3)*X(1,4)*D(2,3)*D(3,2)-X(1,1)*X(1,4)*D(2,3)&
&*D(3,2)-X(1,1)*X(1,3)*D(2,3)*D(3,2)+tt1*D(2,3)*D(3,2)-D(1,3)*X(2,&
&1)*X(2,2)*D(3,2)+D(1,3)*tt3*D(3,2)+X(1,2)*D(1,3)*X(1,4)*D(3,2)-X(&
&1,1)*D(1,3)*X(1,4)*D(3,2)-X(1,1)*X(1,2)*D(1,3)*D(3,2)+tt1*D(1,3)*&
&D(3,2)+D(2,2)*D(2,3)*tt7+D(1,2)*D(2,3)*tt7+D(1,3)*D(2,2)*tt7+D(1,&
&2)*D(1,3)*tt7+D(2,2)*D(2,3)*tt6+D(1,2)*X(2,2)*D(2,3)*X(2,3)-2*X(2&
&,1)*D(2,2)*D(2,3)*X(2,3)-D(1,2)*X(2,1)*D(2,3)*X(2,3)+D(1,3)*D(2,2&
&)*X(2,2)*X(2,3)-D(1,3)*X(2,1)*D(2,2)*X(2,3)-D(1,2)*X(2,1)*X(2,2)*&
&D(2,3)+tt3*D(2,2)*D(2,3)+tt4*D(2,2)*D(2,3)-2*X(1,1)*X(1,3)*D(2,2)&
&*D(2,3)+tt1*D(2,2)*D(2,3)+D(1,2)*tt3*D(2,3)+D(1,2)*X(1,2)*X(1,3)*&
&D(2,3)-X(1,1)*D(1,2)*X(1,3)*D(2,3)-X(1,1)*D(1,2)*X(1,2)*D(2,3)+tt&
&1*D(1,2)*D(2,3)+D(1,2)*D(1,3)*tt5-D(1,3)*X(2,1)*D(2,2)*X(2,2)-2*D&
&(1,2)*D(1,3)*X(2,1)*X(2,2)+D(1,3)*tt3*D(2,2)+X(1,2)*D(1,3)*X(1,3)&
&*D(2,2)-X(1,1)*D(1,3)*X(1,3)*D(2,2)-X(1,1)*X(1,2)*D(1,3)*D(2,2)+t&
&t1*D(1,3)*D(2,2)+D(1,2)*D(1,3)*tt3+D(1,2)*tt2*D(1,3)-2*X(1,1)*D(1&
&,2)*X(1,2)*D(1,3)+tt1*D(1,2)*D(1,3))**2+5.0E-1*(D(3,1)*D(3,3)*tt1&
&2+D(2,1)*D(3,3)*X(3,3)*X(3,4)+D(2,3)*D(3,1)*X(3,3)*X(3,4)+D(1,1)*&
&X(3,2)*D(3,3)*X(3,4)-2*D(3,1)*X(3,1)*D(3,3)*X(3,4)-D(2,1)*X(3,1)*&
&D(3,3)*X(3,4)-D(1,1)*X(3,1)*D(3,3)*X(3,4)+D(1,3)*D(3,1)*X(3,2)*X(&
&3,4)-D(2,3)*D(3,1)*X(3,1)*X(3,4)-D(1,3)*D(3,1)*X(3,1)*X(3,4)+D(2,&
&1)*D(2,3)*tt11-D(2,1)*X(3,1)*D(3,3)*X(3,3)+D(1,1)*D(2,3)*X(3,2)*X&
&(3,3)+D(1,3)*D(2,1)*X(3,2)*X(3,3)-D(2,3)*D(3,1)*X(3,1)*X(3,3)-2*D&
&(2,1)*D(2,3)*X(3,1)*X(3,3)-D(1,1)*D(2,3)*X(3,1)*X(3,3)-D(1,3)*D(2&
&,1)*X(3,1)*X(3,3)-D(1,1)*X(3,1)*X(3,2)*D(3,3)+D(3,1)*tt7*D(3,3)+D&
&(2,1)*tt7*D(3,3)+D(1,1)*tt7*D(3,3)+tt9*D(3,1)*D(3,3)-2*X(2,1)*X(2&
&,4)*D(3,1)*D(3,3)+tt3*D(3,1)*D(3,3)+tt8*D(3,1)*D(3,3)-2*X(1,1)*X(&
&1,4)*D(3,1)*D(3,3)+tt1*D(3,1)*D(3,3)+D(2,1)*X(2,3)*X(2,4)*D(3,3)+&
&D(1,1)*X(2,2)*X(2,4)*D(3,3)-D(2,1)*X(2,1)*X(2,4)*D(3,3)-D(1,1)*X(&
&2,1)*X(2,4)*D(3,3)-D(2,1)*X(2,1)*X(2,3)*D(3,3)-D(1,1)*X(2,1)*X(2,&
&2)*D(3,3)+D(2,1)*tt3*D(3,3)+D(1,1)*tt3*D(3,3)+X(1,3)*X(1,4)*D(2,1&
&)*D(3,3)-X(1,1)*X(1,4)*D(2,1)*D(3,3)-X(1,1)*X(1,3)*D(2,1)*D(3,3)+&
&tt1*D(2,1)*D(3,3)+D(1,1)*X(1,2)*X(1,4)*D(3,3)-D(1,1)*X(1,1)*X(1,4&
&)*D(3,3)-D(1,1)*X(1,1)*X(1,2)*D(3,3)+D(1,1)*tt1*D(3,3)+D(1,1)*D(1&
&,3)*tt10-D(1,3)*D(3,1)*X(3,1)*X(3,2)-D(1,1)*D(2,3)*X(3,1)*X(3,2)-&
&D(1,3)*D(2,1)*X(3,1)*X(3,2)-2*D(1,1)*D(1,3)*X(3,1)*X(3,2)+D(2,3)*&
&D(3,1)*tt7+D(1,3)*D(3,1)*tt7+D(2,1)*D(2,3)*tt7+D(1,1)*D(2,3)*tt7+&
&D(1,3)*D(2,1)*tt7+D(1,1)*D(1,3)*tt7+D(2,3)*X(2,3)*X(2,4)*D(3,1)-X&
&(2,1)*D(2,3)*X(2,4)*D(3,1)+D(1,3)*X(2,2)*X(2,4)*D(3,1)-D(1,3)*X(2&
&,1)*X(2,4)*D(3,1)-X(2,1)*D(2,3)*X(2,3)*D(3,1)+tt3*D(2,3)*D(3,1)+X&
&(1,3)*X(1,4)*D(2,3)*D(3,1)-X(1,1)*X(1,4)*D(2,3)*D(3,1)-X(1,1)*X(1&
&,3)*D(2,3)*D(3,1)+tt1*D(2,3)*D(3,1)-D(1,3)*X(2,1)*X(2,2)*D(3,1)+D&
&(1,3)*tt3*D(3,1)+X(1,2)*D(1,3)*X(1,4)*D(3,1)-X(1,1)*D(1,3)*X(1,4)&
&*D(3,1)-X(1,1)*X(1,2)*D(1,3)*D(3,1)+tt1*D(1,3)*D(3,1)+D(2,1)*D(2,&
&3)*tt6+D(1,1)*X(2,2)*D(2,3)*X(2,3)-2*D(2,1)*X(2,1)*D(2,3)*X(2,3)-&
&D(1,1)*X(2,1)*D(2,3)*X(2,3)+D(1,3)*D(2,1)*X(2,2)*X(2,3)-D(1,3)*D(&
&2,1)*X(2,1)*X(2,3)-D(1,1)*X(2,1)*X(2,2)*D(2,3)+D(2,1)*tt3*D(2,3)+&
&D(1,1)*tt3*D(2,3)+tt4*D(2,1)*D(2,3)-2*X(1,1)*X(1,3)*D(2,1)*D(2,3)&
&+tt1*D(2,1)*D(2,3)+D(1,1)*X(1,2)*X(1,3)*D(2,3)-D(1,1)*X(1,1)*X(1,&
&3)*D(2,3)-D(1,1)*X(1,1)*X(1,2)*D(2,3)+D(1,1)*tt1*D(2,3)+D(1,1)*D(&
&1,3)*tt5-D(1,3)*D(2,1)*X(2,1)*X(2,2)-2*D(1,1)*D(1,3)*X(2,1)*X(2,2&
&)+D(1,3)*D(2,1)*tt3+D(1,1)*D(1,3)*tt3+X(1,2)*D(1,3)*X(1,3)*D(2,1)&
&-X(1,1)*D(1,3)*X(1,3)*D(2,1)-X(1,1)*X(1,2)*D(1,3)*D(2,1)+tt1*D(1,&
&3)*D(2,1)+D(1,1)*tt2*D(1,3)-2*D(1,1)*X(1,1)*X(1,2)*D(1,3)+D(1,1)*&
&tt1*D(1,3))**2+5.0E-1*(D(3,1)*D(3,2)*tt12+D(2,1)*D(3,2)*X(3,3)*X(&
&3,4)+D(2,2)*D(3,1)*X(3,3)*X(3,4)+D(1,1)*D(3,2)*X(3,2)*X(3,4)+D(1,&
&2)*D(3,1)*X(3,2)*X(3,4)-2*D(3,1)*X(3,1)*D(3,2)*X(3,4)-D(2,1)*X(3,&
&1)*D(3,2)*X(3,4)-D(1,1)*X(3,1)*D(3,2)*X(3,4)-D(2,2)*D(3,1)*X(3,1)&
&*X(3,4)-D(1,2)*D(3,1)*X(3,1)*X(3,4)+D(2,1)*D(2,2)*tt11+D(1,1)*D(2&
&,2)*X(3,2)*X(3,3)+D(1,2)*D(2,1)*X(3,2)*X(3,3)-D(2,1)*X(3,1)*D(3,2&
&)*X(3,3)-D(2,2)*D(3,1)*X(3,1)*X(3,3)-2*D(2,1)*D(2,2)*X(3,1)*X(3,3&
&)-D(1,1)*D(2,2)*X(3,1)*X(3,3)-D(1,2)*D(2,1)*X(3,1)*X(3,3)+D(1,1)*&
&D(1,2)*tt10-D(1,1)*X(3,1)*D(3,2)*X(3,2)-D(1,2)*D(3,1)*X(3,1)*X(3,&
&2)-D(1,1)*D(2,2)*X(3,1)*X(3,2)-D(1,2)*D(2,1)*X(3,1)*X(3,2)-2*D(1,&
&1)*D(1,2)*X(3,1)*X(3,2)+D(3,1)*tt7*D(3,2)+D(2,1)*tt7*D(3,2)+D(1,1&
&)*tt7*D(3,2)+tt9*D(3,1)*D(3,2)-2*X(2,1)*X(2,4)*D(3,1)*D(3,2)+tt3*&
&D(3,1)*D(3,2)+tt8*D(3,1)*D(3,2)-2*X(1,1)*X(1,4)*D(3,1)*D(3,2)+tt1&
&*D(3,1)*D(3,2)+D(2,1)*X(2,3)*X(2,4)*D(3,2)+D(1,1)*X(2,2)*X(2,4)*D&
&(3,2)-D(2,1)*X(2,1)*X(2,4)*D(3,2)-D(1,1)*X(2,1)*X(2,4)*D(3,2)-D(2&
&,1)*X(2,1)*X(2,3)*D(3,2)-D(1,1)*X(2,1)*X(2,2)*D(3,2)+D(2,1)*tt3*D&
&(3,2)+D(1,1)*tt3*D(3,2)+X(1,3)*X(1,4)*D(2,1)*D(3,2)-X(1,1)*X(1,4)&
&*D(2,1)*D(3,2)-X(1,1)*X(1,3)*D(2,1)*D(3,2)+tt1*D(2,1)*D(3,2)+D(1,&
&1)*X(1,2)*X(1,4)*D(3,2)-D(1,1)*X(1,1)*X(1,4)*D(3,2)-D(1,1)*X(1,1)&
&*X(1,2)*D(3,2)+D(1,1)*tt1*D(3,2)+D(2,2)*D(3,1)*tt7+D(1,2)*D(3,1)*&
&tt7+D(2,1)*D(2,2)*tt7+D(1,1)*D(2,2)*tt7+D(1,2)*D(2,1)*tt7+D(1,1)*&
&D(1,2)*tt7+D(2,2)*X(2,3)*X(2,4)*D(3,1)+D(1,2)*X(2,2)*X(2,4)*D(3,1&
&)-X(2,1)*D(2,2)*X(2,4)*D(3,1)-D(1,2)*X(2,1)*X(2,4)*D(3,1)-X(2,1)*&
&D(2,2)*X(2,3)*D(3,1)-D(1,2)*X(2,1)*X(2,2)*D(3,1)+tt3*D(2,2)*D(3,1&
&)+X(1,3)*X(1,4)*D(2,2)*D(3,1)-X(1,1)*X(1,4)*D(2,2)*D(3,1)-X(1,1)*&
&X(1,3)*D(2,2)*D(3,1)+tt1*D(2,2)*D(3,1)+D(1,2)*tt3*D(3,1)+D(1,2)*X&
&(1,2)*X(1,4)*D(3,1)-X(1,1)*D(1,2)*X(1,4)*D(3,1)-X(1,1)*D(1,2)*X(1&
&,2)*D(3,1)+tt1*D(1,2)*D(3,1)+D(2,1)*D(2,2)*tt6+D(1,1)*D(2,2)*X(2,&
&2)*X(2,3)+D(1,2)*D(2,1)*X(2,2)*X(2,3)-2*D(2,1)*X(2,1)*D(2,2)*X(2,&
&3)-D(1,1)*X(2,1)*D(2,2)*X(2,3)-D(1,2)*D(2,1)*X(2,1)*X(2,3)+D(1,1)&
&*D(1,2)*tt5-D(1,1)*X(2,1)*D(2,2)*X(2,2)-D(1,2)*D(2,1)*X(2,1)*X(2,&
&2)-2*D(1,1)*D(1,2)*X(2,1)*X(2,2)+D(2,1)*tt3*D(2,2)+D(1,1)*tt3*D(2&
&,2)+tt4*D(2,1)*D(2,2)-2*X(1,1)*X(1,3)*D(2,1)*D(2,2)+tt1*D(2,1)*D(&
&2,2)+D(1,1)*X(1,2)*X(1,3)*D(2,2)-D(1,1)*X(1,1)*X(1,3)*D(2,2)-D(1,&
&1)*X(1,1)*X(1,2)*D(2,2)+D(1,1)*tt1*D(2,2)+D(1,2)*D(2,1)*tt3+D(1,1&
&)*D(1,2)*tt3+D(1,2)*X(1,2)*X(1,3)*D(2,1)-X(1,1)*D(1,2)*X(1,3)*D(2&
&,1)-X(1,1)*D(1,2)*X(1,2)*D(2,1)+tt1*D(1,2)*D(2,1)+D(1,1)*D(1,2)*t&
&t2-2*D(1,1)*X(1,1)*D(1,2)*X(1,2)+D(1,1)*tt1*D(1,2))**2))
END 
SUBROUTINE tet_stvk_jac(jac, X, D, volume, lam, miu) 
IMPLICIT NONE 
REAL(KIND=8) jac(1, 12) 
REAL(KIND=8) X(3, 4) 
REAL(KIND=8) D(3, 3) 
REAL(KIND=8) volume(1, 1) 
REAL(KIND=8) lam(1, 1) 
REAL(KIND=8) miu(1, 1) 
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
tt1 = X(1,1)**2
tt2 = X(1,2)**2
tt3 = X(2,1)**2
tt4 = X(1,3)**2
tt5 = X(2,2)**2
tt6 = X(2,3)**2
tt7 = X(3,1)**2
tt8 = X(1,4)**2
tt9 = X(2,4)**2
tt10 = X(3,2)**2
tt11 = X(3,3)**2
tt12 = X(3,4)**2
tt13 = D(3,1)*D(3,2)*tt12+D(2,1)*D(3,2)*X(3,3)*X(3,4)+D(2,2)*D(3,&
&1)*X(3,3)*X(3,4)+D(1,1)*D(3,2)*X(3,2)*X(3,4)+D(1,2)*D(3,1)*X(3,2)&
&*X(3,4)-2*D(3,1)*X(3,1)*D(3,2)*X(3,4)-D(2,1)*X(3,1)*D(3,2)*X(3,4)&
&-D(1,1)*X(3,1)*D(3,2)*X(3,4)-D(2,2)*D(3,1)*X(3,1)*X(3,4)-D(1,2)*D&
&(3,1)*X(3,1)*X(3,4)+D(2,1)*D(2,2)*tt11+D(1,1)*D(2,2)*X(3,2)*X(3,3&
&)+D(1,2)*D(2,1)*X(3,2)*X(3,3)-D(2,1)*X(3,1)*D(3,2)*X(3,3)-D(2,2)*&
&D(3,1)*X(3,1)*X(3,3)-2*D(2,1)*D(2,2)*X(3,1)*X(3,3)-D(1,1)*D(2,2)*&
&X(3,1)*X(3,3)-D(1,2)*D(2,1)*X(3,1)*X(3,3)+D(1,1)*D(1,2)*tt10-D(1,&
&1)*X(3,1)*D(3,2)*X(3,2)-D(1,2)*D(3,1)*X(3,1)*X(3,2)-D(1,1)*D(2,2)&
&*X(3,1)*X(3,2)-D(1,2)*D(2,1)*X(3,1)*X(3,2)-2*D(1,1)*D(1,2)*X(3,1)&
&*X(3,2)+D(3,1)*tt7*D(3,2)+D(2,1)*tt7*D(3,2)+D(1,1)*tt7*D(3,2)+tt9&
&*D(3,1)*D(3,2)-2*X(2,1)*X(2,4)*D(3,1)*D(3,2)+tt3*D(3,1)*D(3,2)+tt&
&8*D(3,1)*D(3,2)-2*X(1,1)*X(1,4)*D(3,1)*D(3,2)+tt1*D(3,1)*D(3,2)+D&
&(2,1)*X(2,3)*X(2,4)*D(3,2)+D(1,1)*X(2,2)*X(2,4)*D(3,2)-D(2,1)*X(2&
&,1)*X(2,4)*D(3,2)-D(1,1)*X(2,1)*X(2,4)*D(3,2)-D(2,1)*X(2,1)*X(2,3&
&)*D(3,2)-D(1,1)*X(2,1)*X(2,2)*D(3,2)+D(2,1)*tt3*D(3,2)+D(1,1)*tt3&
&*D(3,2)+X(1,3)*X(1,4)*D(2,1)*D(3,2)-X(1,1)*X(1,4)*D(2,1)*D(3,2)-X&
&(1,1)*X(1,3)*D(2,1)*D(3,2)+tt1*D(2,1)*D(3,2)+D(1,1)*X(1,2)*X(1,4)&
&*D(3,2)-D(1,1)*X(1,1)*X(1,4)*D(3,2)-D(1,1)*X(1,1)*X(1,2)*D(3,2)+D&
&(1,1)*tt1*D(3,2)+D(2,2)*D(3,1)*tt7+D(1,2)*D(3,1)*tt7+D(2,1)*D(2,2&
&)*tt7+D(1,1)*D(2,2)*tt7+D(1,2)*D(2,1)*tt7+D(1,1)*D(1,2)*tt7+D(2,2&
&)*X(2,3)*X(2,4)*D(3,1)+D(1,2)*X(2,2)*X(2,4)*D(3,1)-X(2,1)*D(2,2)*&
&X(2,4)*D(3,1)-D(1,2)*X(2,1)*X(2,4)*D(3,1)-X(2,1)*D(2,2)*X(2,3)*D(&
&3,1)-D(1,2)*X(2,1)*X(2,2)*D(3,1)+tt3*D(2,2)*D(3,1)+X(1,3)*X(1,4)*&
&D(2,2)*D(3,1)-X(1,1)*X(1,4)*D(2,2)*D(3,1)-X(1,1)*X(1,3)*D(2,2)*D(&
&3,1)+tt1*D(2,2)*D(3,1)+D(1,2)*tt3*D(3,1)+D(1,2)*X(1,2)*X(1,4)*D(3&
&,1)-X(1,1)*D(1,2)*X(1,4)*D(3,1)-X(1,1)*D(1,2)*X(1,2)*D(3,1)+tt1*D&
&(1,2)*D(3,1)+D(2,1)*D(2,2)*tt6+D(1,1)*D(2,2)*X(2,2)*X(2,3)+D(1,2)&
&*D(2,1)*X(2,2)*X(2,3)-2*D(2,1)*X(2,1)*D(2,2)*X(2,3)-D(1,1)*X(2,1)&
&*D(2,2)*X(2,3)-D(1,2)*D(2,1)*X(2,1)*X(2,3)+D(1,1)*D(1,2)*tt5-D(1,&
&1)*X(2,1)*D(2,2)*X(2,2)-D(1,2)*D(2,1)*X(2,1)*X(2,2)-2*D(1,1)*D(1,&
&2)*X(2,1)*X(2,2)+D(2,1)*tt3*D(2,2)+D(1,1)*tt3*D(2,2)+tt4*D(2,1)*D&
&(2,2)-2*X(1,1)*X(1,3)*D(2,1)*D(2,2)+tt1*D(2,1)*D(2,2)+D(1,1)*X(1,&
&2)*X(1,3)*D(2,2)-D(1,1)*X(1,1)*X(1,3)*D(2,2)-D(1,1)*X(1,1)*X(1,2)&
&*D(2,2)+D(1,1)*tt1*D(2,2)+D(1,2)*D(2,1)*tt3+D(1,1)*D(1,2)*tt3+D(1&
&,2)*X(1,2)*X(1,3)*D(2,1)-X(1,1)*D(1,2)*X(1,3)*D(2,1)-X(1,1)*D(1,2&
&)*X(1,2)*D(2,1)+tt1*D(1,2)*D(2,1)+D(1,1)*D(1,2)*tt2-2*D(1,1)*X(1,&
&1)*D(1,2)*X(1,2)+D(1,1)*tt1*D(1,2)
tt14 = D(3,1)*D(3,3)*tt12+D(2,1)*D(3,3)*X(3,3)*X(3,4)+D(2,3)*D(3,&
&1)*X(3,3)*X(3,4)+D(1,1)*X(3,2)*D(3,3)*X(3,4)-2*D(3,1)*X(3,1)*D(3,&
&3)*X(3,4)-D(2,1)*X(3,1)*D(3,3)*X(3,4)-D(1,1)*X(3,1)*D(3,3)*X(3,4)&
&+D(1,3)*D(3,1)*X(3,2)*X(3,4)-D(2,3)*D(3,1)*X(3,1)*X(3,4)-D(1,3)*D&
&(3,1)*X(3,1)*X(3,4)+D(2,1)*D(2,3)*tt11-D(2,1)*X(3,1)*D(3,3)*X(3,3&
&)+D(1,1)*D(2,3)*X(3,2)*X(3,3)+D(1,3)*D(2,1)*X(3,2)*X(3,3)-D(2,3)*&
&D(3,1)*X(3,1)*X(3,3)-2*D(2,1)*D(2,3)*X(3,1)*X(3,3)-D(1,1)*D(2,3)*&
&X(3,1)*X(3,3)-D(1,3)*D(2,1)*X(3,1)*X(3,3)-D(1,1)*X(3,1)*X(3,2)*D(&
&3,3)+D(3,1)*tt7*D(3,3)+D(2,1)*tt7*D(3,3)+D(1,1)*tt7*D(3,3)+tt9*D(&
&3,1)*D(3,3)-2*X(2,1)*X(2,4)*D(3,1)*D(3,3)+tt3*D(3,1)*D(3,3)+tt8*D&
&(3,1)*D(3,3)-2*X(1,1)*X(1,4)*D(3,1)*D(3,3)+tt1*D(3,1)*D(3,3)+D(2,&
&1)*X(2,3)*X(2,4)*D(3,3)+D(1,1)*X(2,2)*X(2,4)*D(3,3)-D(2,1)*X(2,1)&
&*X(2,4)*D(3,3)-D(1,1)*X(2,1)*X(2,4)*D(3,3)-D(2,1)*X(2,1)*X(2,3)*D&
&(3,3)-D(1,1)*X(2,1)*X(2,2)*D(3,3)+D(2,1)*tt3*D(3,3)+D(1,1)*tt3*D(&
&3,3)+X(1,3)*X(1,4)*D(2,1)*D(3,3)-X(1,1)*X(1,4)*D(2,1)*D(3,3)-X(1,&
&1)*X(1,3)*D(2,1)*D(3,3)+tt1*D(2,1)*D(3,3)+D(1,1)*X(1,2)*X(1,4)*D(&
&3,3)-D(1,1)*X(1,1)*X(1,4)*D(3,3)-D(1,1)*X(1,1)*X(1,2)*D(3,3)+D(1,&
&1)*tt1*D(3,3)+D(1,1)*D(1,3)*tt10-D(1,3)*D(3,1)*X(3,1)*X(3,2)-D(1,&
&1)*D(2,3)*X(3,1)*X(3,2)-D(1,3)*D(2,1)*X(3,1)*X(3,2)-2*D(1,1)*D(1,&
&3)*X(3,1)*X(3,2)+D(2,3)*D(3,1)*tt7+D(1,3)*D(3,1)*tt7+D(2,1)*D(2,3&
&)*tt7+D(1,1)*D(2,3)*tt7+D(1,3)*D(2,1)*tt7+D(1,1)*D(1,3)*tt7+D(2,3&
&)*X(2,3)*X(2,4)*D(3,1)-X(2,1)*D(2,3)*X(2,4)*D(3,1)+D(1,3)*X(2,2)*&
&X(2,4)*D(3,1)-D(1,3)*X(2,1)*X(2,4)*D(3,1)-X(2,1)*D(2,3)*X(2,3)*D(&
&3,1)+tt3*D(2,3)*D(3,1)+X(1,3)*X(1,4)*D(2,3)*D(3,1)-X(1,1)*X(1,4)*&
&D(2,3)*D(3,1)-X(1,1)*X(1,3)*D(2,3)*D(3,1)+tt1*D(2,3)*D(3,1)-D(1,3&
&)*X(2,1)*X(2,2)*D(3,1)+D(1,3)*tt3*D(3,1)+X(1,2)*D(1,3)*X(1,4)*D(3&
&,1)-X(1,1)*D(1,3)*X(1,4)*D(3,1)-X(1,1)*X(1,2)*D(1,3)*D(3,1)+tt1*D&
&(1,3)*D(3,1)+D(2,1)*D(2,3)*tt6+D(1,1)*X(2,2)*D(2,3)*X(2,3)-2*D(2,&
&1)*X(2,1)*D(2,3)*X(2,3)-D(1,1)*X(2,1)*D(2,3)*X(2,3)+D(1,3)*D(2,1)&
&*X(2,2)*X(2,3)-D(1,3)*D(2,1)*X(2,1)*X(2,3)-D(1,1)*X(2,1)*X(2,2)*D&
&(2,3)+D(2,1)*tt3*D(2,3)+D(1,1)*tt3*D(2,3)+tt4*D(2,1)*D(2,3)-2*X(1&
&,1)*X(1,3)*D(2,1)*D(2,3)+tt1*D(2,1)*D(2,3)+D(1,1)*X(1,2)*X(1,3)*D&
&(2,3)-D(1,1)*X(1,1)*X(1,3)*D(2,3)-D(1,1)*X(1,1)*X(1,2)*D(2,3)+D(1&
&,1)*tt1*D(2,3)+D(1,1)*D(1,3)*tt5-D(1,3)*D(2,1)*X(2,1)*X(2,2)-2*D(&
&1,1)*D(1,3)*X(2,1)*X(2,2)+D(1,3)*D(2,1)*tt3+D(1,1)*D(1,3)*tt3+X(1&
&,2)*D(1,3)*X(1,3)*D(2,1)-X(1,1)*D(1,3)*X(1,3)*D(2,1)-X(1,1)*X(1,2&
&)*D(1,3)*D(2,1)+tt1*D(1,3)*D(2,1)+D(1,1)*tt2*D(1,3)-2*D(1,1)*X(1,&
&1)*X(1,2)*D(1,3)+D(1,1)*tt1*D(1,3)
tt15 = D(3,2)*D(3,3)*tt12+D(2,2)*D(3,3)*X(3,3)*X(3,4)+D(2,3)*D(3,&
&2)*X(3,3)*X(3,4)+D(1,2)*X(3,2)*D(3,3)*X(3,4)-2*X(3,1)*D(3,2)*D(3,&
&3)*X(3,4)-D(2,2)*X(3,1)*D(3,3)*X(3,4)-D(1,2)*X(3,1)*D(3,3)*X(3,4)&
&+D(1,3)*D(3,2)*X(3,2)*X(3,4)-D(2,3)*X(3,1)*D(3,2)*X(3,4)-D(1,3)*X&
&(3,1)*D(3,2)*X(3,4)+D(2,2)*D(2,3)*tt11-D(2,2)*X(3,1)*D(3,3)*X(3,3&
&)+D(1,2)*D(2,3)*X(3,2)*X(3,3)+D(1,3)*D(2,2)*X(3,2)*X(3,3)-D(2,3)*&
&X(3,1)*D(3,2)*X(3,3)-2*D(2,2)*D(2,3)*X(3,1)*X(3,3)-D(1,2)*D(2,3)*&
&X(3,1)*X(3,3)-D(1,3)*D(2,2)*X(3,1)*X(3,3)-D(1,2)*X(3,1)*X(3,2)*D(&
&3,3)+tt7*D(3,2)*D(3,3)+tt9*D(3,2)*D(3,3)-2*X(2,1)*X(2,4)*D(3,2)*D&
&(3,3)+tt3*D(3,2)*D(3,3)+tt8*D(3,2)*D(3,3)-2*X(1,1)*X(1,4)*D(3,2)*&
&D(3,3)+tt1*D(3,2)*D(3,3)+D(2,2)*tt7*D(3,3)+D(1,2)*tt7*D(3,3)+D(2,&
&2)*X(2,3)*X(2,4)*D(3,3)+D(1,2)*X(2,2)*X(2,4)*D(3,3)-X(2,1)*D(2,2)&
&*X(2,4)*D(3,3)-D(1,2)*X(2,1)*X(2,4)*D(3,3)-X(2,1)*D(2,2)*X(2,3)*D&
&(3,3)-D(1,2)*X(2,1)*X(2,2)*D(3,3)+tt3*D(2,2)*D(3,3)+X(1,3)*X(1,4)&
&*D(2,2)*D(3,3)-X(1,1)*X(1,4)*D(2,2)*D(3,3)-X(1,1)*X(1,3)*D(2,2)*D&
&(3,3)+tt1*D(2,2)*D(3,3)+D(1,2)*tt3*D(3,3)+D(1,2)*X(1,2)*X(1,4)*D(&
&3,3)-X(1,1)*D(1,2)*X(1,4)*D(3,3)-X(1,1)*D(1,2)*X(1,2)*D(3,3)+tt1*&
&D(1,2)*D(3,3)+D(1,2)*D(1,3)*tt10-D(1,3)*X(3,1)*D(3,2)*X(3,2)-D(1,&
&2)*D(2,3)*X(3,1)*X(3,2)-D(1,3)*D(2,2)*X(3,1)*X(3,2)-2*D(1,2)*D(1,&
&3)*X(3,1)*X(3,2)+D(2,3)*tt7*D(3,2)+D(1,3)*tt7*D(3,2)+D(2,3)*X(2,3&
&)*X(2,4)*D(3,2)-X(2,1)*D(2,3)*X(2,4)*D(3,2)+D(1,3)*X(2,2)*X(2,4)*&
&D(3,2)-D(1,3)*X(2,1)*X(2,4)*D(3,2)-X(2,1)*D(2,3)*X(2,3)*D(3,2)+tt&
&3*D(2,3)*D(3,2)+X(1,3)*X(1,4)*D(2,3)*D(3,2)-X(1,1)*X(1,4)*D(2,3)*&
&D(3,2)-X(1,1)*X(1,3)*D(2,3)*D(3,2)+tt1*D(2,3)*D(3,2)-D(1,3)*X(2,1&
&)*X(2,2)*D(3,2)+D(1,3)*tt3*D(3,2)+X(1,2)*D(1,3)*X(1,4)*D(3,2)-X(1&
&,1)*D(1,3)*X(1,4)*D(3,2)-X(1,1)*X(1,2)*D(1,3)*D(3,2)+tt1*D(1,3)*D&
&(3,2)+D(2,2)*D(2,3)*tt7+D(1,2)*D(2,3)*tt7+D(1,3)*D(2,2)*tt7+D(1,2&
&)*D(1,3)*tt7+D(2,2)*D(2,3)*tt6+D(1,2)*X(2,2)*D(2,3)*X(2,3)-2*X(2,&
&1)*D(2,2)*D(2,3)*X(2,3)-D(1,2)*X(2,1)*D(2,3)*X(2,3)+D(1,3)*D(2,2)&
&*X(2,2)*X(2,3)-D(1,3)*X(2,1)*D(2,2)*X(2,3)-D(1,2)*X(2,1)*X(2,2)*D&
&(2,3)+tt3*D(2,2)*D(2,3)+tt4*D(2,2)*D(2,3)-2*X(1,1)*X(1,3)*D(2,2)*&
&D(2,3)+tt1*D(2,2)*D(2,3)+D(1,2)*tt3*D(2,3)+D(1,2)*X(1,2)*X(1,3)*D&
&(2,3)-X(1,1)*D(1,2)*X(1,3)*D(2,3)-X(1,1)*D(1,2)*X(1,2)*D(2,3)+tt1&
&*D(1,2)*D(2,3)+D(1,2)*D(1,3)*tt5-D(1,3)*X(2,1)*D(2,2)*X(2,2)-2*D(&
&1,2)*D(1,3)*X(2,1)*X(2,2)+D(1,3)*tt3*D(2,2)+X(1,2)*D(1,3)*X(1,3)*&
&D(2,2)-X(1,1)*D(1,3)*X(1,3)*D(2,2)-X(1,1)*X(1,2)*D(1,3)*D(2,2)+tt&
&1*D(1,3)*D(2,2)+D(1,2)*D(1,3)*tt3+D(1,2)*tt2*D(1,3)-2*X(1,1)*D(1,&
&2)*X(1,2)*D(1,3)+tt1*D(1,2)*D(1,3)
tt16 = (-D(3,1))-D(2,1)-D(1,1)
tt17 = X(1,4)*D(3,1)-X(1,1)*D(3,1)+X(1,3)*D(2,1)-X(1,1)*D(2,1)+D(&
&1,1)*X(1,2)-D(1,1)*X(1,1)
tt18 = X(2,4)*D(3,1)-X(2,1)*D(3,1)+D(2,1)*X(2,3)+D(1,1)*X(2,2)-D(&
&2,1)*X(2,1)-D(1,1)*X(2,1)
tt19 = D(3,1)*X(3,4)+D(2,1)*X(3,3)+D(1,1)*X(3,2)-D(3,1)*X(3,1)-D(&
&2,1)*X(3,1)-D(1,1)*X(3,1)
tt20 = tt19**2+tt18**2+tt17**2-1
tt21 = (-D(3,2))-D(2,2)-D(1,2)
tt22 = X(1,4)*D(3,2)-X(1,1)*D(3,2)+X(1,3)*D(2,2)-X(1,1)*D(2,2)+D(&
&1,2)*X(1,2)-X(1,1)*D(1,2)
tt23 = X(2,4)*D(3,2)-X(2,1)*D(3,2)+D(2,2)*X(2,3)+D(1,2)*X(2,2)-X(&
&2,1)*D(2,2)-D(1,2)*X(2,1)
tt24 = D(3,2)*X(3,4)+D(2,2)*X(3,3)+D(1,2)*X(3,2)-X(3,1)*D(3,2)-D(&
&2,2)*X(3,1)-D(1,2)*X(3,1)
tt25 = tt24**2+tt23**2+tt22**2-1
tt26 = (-D(3,3))-D(2,3)-D(1,3)
tt27 = X(1,4)*D(3,3)-X(1,1)*D(3,3)+X(1,3)*D(2,3)-X(1,1)*D(2,3)+X(&
&1,2)*D(1,3)-X(1,1)*D(1,3)
tt28 = X(2,4)*D(3,3)-X(2,1)*D(3,3)+D(2,3)*X(2,3)-X(2,1)*D(2,3)+D(&
&1,3)*X(2,2)-D(1,3)*X(2,1)
tt29 = D(3,3)*X(3,4)+D(2,3)*X(3,3)-X(3,1)*D(3,3)+D(1,3)*X(3,2)-D(&
&2,3)*X(3,1)-D(1,3)*X(3,1)
tt30 = tt29**2+tt28**2+tt27**2-1
tt31 = -X(1,1)
tt32 = X(1,2)+tt31
tt33 = X(1,3)+tt31
tt34 = X(1,4)+tt31
tt35 = tt34*D(3,1)+tt33*D(2,1)+D(1,1)*tt32
tt36 = tt34*D(3,2)+tt33*D(2,2)+D(1,2)*tt32
tt37 = tt34*D(3,3)+tt33*D(2,3)+tt32*D(1,3)
tt38 = -X(2,1)
tt39 = X(2,2)+tt38
tt40 = X(2,3)+tt38
tt41 = X(2,4)+tt38
tt42 = tt41*D(3,1)+D(2,1)*tt40+D(1,1)*tt39
tt43 = -X(3,1)
tt44 = X(3,2)+tt43
tt45 = X(3,3)+tt43
tt46 = X(3,4)+tt43
tt47 = D(3,1)*tt46+D(2,1)*tt45+D(1,1)*tt44
tt48 = tt41*D(3,2)+D(2,2)*tt40+D(1,2)*tt39
tt49 = D(3,2)*tt46+D(2,2)*tt45+D(1,2)*tt44
tt50 = tt41*D(3,3)+D(2,3)*tt40+D(1,3)*tt39
tt51 = D(3,3)*tt46+D(2,3)*tt45+D(1,3)*tt44
tt52 = 5.0E-1*(tt51**2+tt50**2+tt37**2-1)+5.0E-1*(tt49**2+tt48**2&
&+tt36**2-1)+5.0E-1*(tt47**2+tt42**2+tt35**2-1)
tt53 = -X(1,1)*D(1,2)*D(2,1)
tt54 = -D(1,1)*X(1,1)*D(2,2)
tt55 = -X(1,1)*D(1,2)*D(3,1)
tt56 = -D(1,1)*X(1,1)*D(3,2)
tt57 = -X(1,1)*D(1,3)*D(2,1)
tt58 = -D(1,1)*X(1,1)*D(2,3)
tt59 = -X(1,1)*D(1,3)*D(3,1)
tt60 = -D(1,1)*X(1,1)*D(3,3)
tt61 = -X(1,1)*D(1,3)*D(2,2)
tt62 = -X(1,1)*D(1,2)*D(2,3)
tt63 = -X(1,1)*D(1,3)*D(3,2)
tt64 = -X(1,1)*D(1,2)*D(3,3)
tt65 = -D(1,2)*D(2,1)*X(2,1)
tt66 = -D(1,1)*X(2,1)*D(2,2)
tt67 = -D(1,2)*X(2,1)*D(3,1)
tt68 = -D(1,1)*X(2,1)*D(3,2)
tt69 = -D(1,3)*D(2,1)*X(2,1)
tt70 = -D(1,1)*X(2,1)*D(2,3)
tt71 = -D(1,3)*X(2,1)*D(3,1)
tt72 = -D(1,1)*X(2,1)*D(3,3)
tt73 = -D(1,3)*X(2,1)*D(2,2)
tt74 = -D(1,2)*X(2,1)*D(2,3)
tt75 = -D(1,3)*X(2,1)*D(3,2)
tt76 = -D(1,2)*X(2,1)*D(3,3)
tt77 = -D(1,2)*D(2,1)*X(3,1)
tt78 = -D(1,1)*D(2,2)*X(3,1)
tt79 = -D(1,2)*D(3,1)*X(3,1)
tt80 = -D(1,1)*X(3,1)*D(3,2)
tt81 = -D(1,3)*D(2,1)*X(3,1)
tt82 = -D(1,1)*D(2,3)*X(3,1)
tt83 = -D(1,3)*D(3,1)*X(3,1)
tt84 = -D(1,1)*X(3,1)*D(3,3)
tt85 = -D(1,3)*D(2,2)*X(3,1)
tt86 = -D(1,2)*D(2,3)*X(3,1)
tt87 = -D(1,3)*X(3,1)*D(3,2)
tt88 = -D(1,2)*X(3,1)*D(3,3)
tt89 = -X(1,1)*D(2,2)*D(3,1)
tt90 = -X(1,1)*D(2,1)*D(3,2)
tt91 = -X(1,1)*D(2,3)*D(3,1)
tt92 = -X(1,1)*D(2,1)*D(3,3)
tt93 = -X(1,1)*D(2,3)*D(3,2)
tt94 = -X(1,1)*D(2,2)*D(3,3)
tt95 = -X(2,1)*D(2,2)*D(3,1)
tt96 = -D(2,1)*X(2,1)*D(3,2)
tt97 = -X(2,1)*D(2,3)*D(3,1)
tt98 = -D(2,1)*X(2,1)*D(3,3)
tt99 = -X(2,1)*D(2,3)*D(3,2)
tt100 = -X(2,1)*D(2,2)*D(3,3)
tt101 = -D(2,2)*D(3,1)*X(3,1)
tt102 = -D(2,1)*X(3,1)*D(3,2)
tt103 = -D(2,3)*D(3,1)*X(3,1)
tt104 = -D(2,1)*X(3,1)*D(3,3)
tt105 = -D(2,3)*X(3,1)*D(3,2)
tt106 = -D(2,2)*X(3,1)*D(3,3)
jac(1,1) = volume(1,1)*(1.0E+0*lam(1,1)*(1.0E+0*tt26*tt37+1.0E+0*&
&tt21*tt36+1.0E+0*tt16*tt35)*tt52+miu(1,1)*(1.0E+0*tt26*tt27*tt30+&
&1.0E+0*tt21*tt22*tt25+1.0E+0*tt16*tt17*tt20+1.0E+0*((-2*X(1,4)*D(&
&3,2)*D(3,3))+2*X(1,1)*D(3,2)*D(3,3)-X(1,4)*D(2,2)*D(3,3)-X(1,3)*D&
&(2,2)*D(3,3)+2*X(1,1)*D(2,2)*D(3,3)-D(1,2)*X(1,4)*D(3,3)-D(1,2)*X&
&(1,2)*D(3,3)+2*X(1,1)*D(1,2)*D(3,3)-X(1,4)*D(2,3)*D(3,2)-X(1,3)*D&
&(2,3)*D(3,2)+2*X(1,1)*D(2,3)*D(3,2)-D(1,3)*X(1,4)*D(3,2)-X(1,2)*D&
&(1,3)*D(3,2)+2*X(1,1)*D(1,3)*D(3,2)-2*X(1,3)*D(2,2)*D(2,3)+2*X(1,&
&1)*D(2,2)*D(2,3)-D(1,2)*X(1,3)*D(2,3)-D(1,2)*X(1,2)*D(2,3)+2*X(1,&
&1)*D(1,2)*D(2,3)-D(1,3)*X(1,3)*D(2,2)-X(1,2)*D(1,3)*D(2,2)+2*X(1,&
&1)*D(1,3)*D(2,2)-2*D(1,2)*X(1,2)*D(1,3)+2*X(1,1)*D(1,2)*D(1,3))*t&
&t15+1.0E+0*((-2*X(1,4)*D(3,1)*D(3,3))+2*X(1,1)*D(3,1)*D(3,3)-X(1,&
&4)*D(2,1)*D(3,3)-X(1,3)*D(2,1)*D(3,3)+2*X(1,1)*D(2,1)*D(3,3)-D(1,&
&1)*X(1,4)*D(3,3)-D(1,1)*X(1,2)*D(3,3)+2*D(1,1)*X(1,1)*D(3,3)-X(1,&
&4)*D(2,3)*D(3,1)-X(1,3)*D(2,3)*D(3,1)+2*X(1,1)*D(2,3)*D(3,1)-D(1,&
&3)*X(1,4)*D(3,1)-X(1,2)*D(1,3)*D(3,1)+2*X(1,1)*D(1,3)*D(3,1)-2*X(&
&1,3)*D(2,1)*D(2,3)+2*X(1,1)*D(2,1)*D(2,3)-D(1,1)*X(1,3)*D(2,3)-D(&
&1,1)*X(1,2)*D(2,3)+2*D(1,1)*X(1,1)*D(2,3)-D(1,3)*X(1,3)*D(2,1)-X(&
&1,2)*D(1,3)*D(2,1)+2*X(1,1)*D(1,3)*D(2,1)-2*D(1,1)*X(1,2)*D(1,3)+&
&2*D(1,1)*X(1,1)*D(1,3))*tt14+1.0E+0*((-2*X(1,4)*D(3,1)*D(3,2))+2*&
&X(1,1)*D(3,1)*D(3,2)-X(1,4)*D(2,1)*D(3,2)-X(1,3)*D(2,1)*D(3,2)+2*&
&X(1,1)*D(2,1)*D(3,2)-D(1,1)*X(1,4)*D(3,2)-D(1,1)*X(1,2)*D(3,2)+2*&
&D(1,1)*X(1,1)*D(3,2)-X(1,4)*D(2,2)*D(3,1)-X(1,3)*D(2,2)*D(3,1)+2*&
&X(1,1)*D(2,2)*D(3,1)-D(1,2)*X(1,4)*D(3,1)-D(1,2)*X(1,2)*D(3,1)+2*&
&X(1,1)*D(1,2)*D(3,1)-2*X(1,3)*D(2,1)*D(2,2)+2*X(1,1)*D(2,1)*D(2,2&
&)-D(1,1)*X(1,3)*D(2,2)-D(1,1)*X(1,2)*D(2,2)+2*D(1,1)*X(1,1)*D(2,2&
&)-D(1,2)*X(1,3)*D(2,1)-D(1,2)*X(1,2)*D(2,1)+2*X(1,1)*D(1,2)*D(2,1&
&)-2*D(1,1)*D(1,2)*X(1,2)+2*D(1,1)*X(1,1)*D(1,2))*tt13))
jac(1,2) = volume(1,1)*(1.0E+0*lam(1,1)*(1.0E+0*tt26*tt50+1.0E+0*&
&tt21*tt48+1.0E+0*tt16*tt42)*tt52+miu(1,1)*(1.0E+0*tt26*tt28*tt30+&
&1.0E+0*tt21*tt23*tt25+1.0E+0*tt16*tt18*tt20+1.0E+0*((-2*X(2,4)*D(&
&3,2)*D(3,3))+2*X(2,1)*D(3,2)*D(3,3)-D(2,2)*X(2,4)*D(3,3)-D(1,2)*X&
&(2,4)*D(3,3)-D(2,2)*X(2,3)*D(3,3)-D(1,2)*X(2,2)*D(3,3)+2*X(2,1)*D&
&(2,2)*D(3,3)+2*D(1,2)*X(2,1)*D(3,3)-D(2,3)*X(2,4)*D(3,2)-D(1,3)*X&
&(2,4)*D(3,2)-D(2,3)*X(2,3)*D(3,2)+2*X(2,1)*D(2,3)*D(3,2)-D(1,3)*X&
&(2,2)*D(3,2)+2*D(1,3)*X(2,1)*D(3,2)-2*D(2,2)*D(2,3)*X(2,3)-D(1,2)&
&*D(2,3)*X(2,3)-D(1,3)*D(2,2)*X(2,3)-D(1,2)*X(2,2)*D(2,3)+2*X(2,1)&
&*D(2,2)*D(2,3)+2*D(1,2)*X(2,1)*D(2,3)-D(1,3)*D(2,2)*X(2,2)-2*D(1,&
&2)*D(1,3)*X(2,2)+2*D(1,3)*X(2,1)*D(2,2)+2*D(1,2)*D(1,3)*X(2,1))*t&
&t15+1.0E+0*((-2*X(2,4)*D(3,1)*D(3,3))+2*X(2,1)*D(3,1)*D(3,3)-D(2,&
&1)*X(2,4)*D(3,3)-D(1,1)*X(2,4)*D(3,3)-D(2,1)*X(2,3)*D(3,3)-D(1,1)&
&*X(2,2)*D(3,3)+2*D(2,1)*X(2,1)*D(3,3)+2*D(1,1)*X(2,1)*D(3,3)-D(2,&
&3)*X(2,4)*D(3,1)-D(1,3)*X(2,4)*D(3,1)-D(2,3)*X(2,3)*D(3,1)+2*X(2,&
&1)*D(2,3)*D(3,1)-D(1,3)*X(2,2)*D(3,1)+2*D(1,3)*X(2,1)*D(3,1)-2*D(&
&2,1)*D(2,3)*X(2,3)-D(1,1)*D(2,3)*X(2,3)-D(1,3)*D(2,1)*X(2,3)-D(1,&
&1)*X(2,2)*D(2,3)+2*D(2,1)*X(2,1)*D(2,3)+2*D(1,1)*X(2,1)*D(2,3)-D(&
&1,3)*D(2,1)*X(2,2)-2*D(1,1)*D(1,3)*X(2,2)+2*D(1,3)*D(2,1)*X(2,1)+&
&2*D(1,1)*D(1,3)*X(2,1))*tt14+1.0E+0*((-2*X(2,4)*D(3,1)*D(3,2))+2*&
&X(2,1)*D(3,1)*D(3,2)-D(2,1)*X(2,4)*D(3,2)-D(1,1)*X(2,4)*D(3,2)-D(&
&2,1)*X(2,3)*D(3,2)-D(1,1)*X(2,2)*D(3,2)+2*D(2,1)*X(2,1)*D(3,2)+2*&
&D(1,1)*X(2,1)*D(3,2)-D(2,2)*X(2,4)*D(3,1)-D(1,2)*X(2,4)*D(3,1)-D(&
&2,2)*X(2,3)*D(3,1)-D(1,2)*X(2,2)*D(3,1)+2*X(2,1)*D(2,2)*D(3,1)+2*&
&D(1,2)*X(2,1)*D(3,1)-2*D(2,1)*D(2,2)*X(2,3)-D(1,1)*D(2,2)*X(2,3)-&
&D(1,2)*D(2,1)*X(2,3)-D(1,1)*D(2,2)*X(2,2)-D(1,2)*D(2,1)*X(2,2)-2*&
&D(1,1)*D(1,2)*X(2,2)+2*D(2,1)*X(2,1)*D(2,2)+2*D(1,1)*X(2,1)*D(2,2&
&)+2*D(1,2)*D(2,1)*X(2,1)+2*D(1,1)*D(1,2)*X(2,1))*tt13))
jac(1,3) = volume(1,1)*(1.0E+0*lam(1,1)*(1.0E+0*tt26*tt51+1.0E+0*&
&tt21*tt49+1.0E+0*tt16*tt47)*tt52+miu(1,1)*(1.0E+0*tt26*tt29*tt30+&
&1.0E+0*tt21*tt24*tt25+1.0E+0*tt16*tt19*tt20+1.0E+0*((-2*D(3,2)*D(&
&3,3)*X(3,4))-D(2,2)*D(3,3)*X(3,4)-D(1,2)*D(3,3)*X(3,4)-D(2,3)*D(3&
&,2)*X(3,4)-D(1,3)*D(3,2)*X(3,4)-D(2,2)*D(3,3)*X(3,3)-D(2,3)*D(3,2&
&)*X(3,3)-2*D(2,2)*D(2,3)*X(3,3)-D(1,2)*D(2,3)*X(3,3)-D(1,3)*D(2,2&
&)*X(3,3)-D(1,2)*X(3,2)*D(3,3)+2*X(3,1)*D(3,2)*D(3,3)+2*D(2,2)*X(3&
&,1)*D(3,3)+2*D(1,2)*X(3,1)*D(3,3)-D(1,3)*D(3,2)*X(3,2)-D(1,2)*D(2&
&,3)*X(3,2)-D(1,3)*D(2,2)*X(3,2)-2*D(1,2)*D(1,3)*X(3,2)+2*D(2,3)*X&
&(3,1)*D(3,2)+2*D(1,3)*X(3,1)*D(3,2)+2*D(2,2)*D(2,3)*X(3,1)+2*D(1,&
&2)*D(2,3)*X(3,1)+2*D(1,3)*D(2,2)*X(3,1)+2*D(1,2)*D(1,3)*X(3,1))*t&
&t15+1.0E+0*((-2*D(3,1)*D(3,3)*X(3,4))-D(2,1)*D(3,3)*X(3,4)-D(1,1)&
&*D(3,3)*X(3,4)-D(2,3)*D(3,1)*X(3,4)-D(1,3)*D(3,1)*X(3,4)-D(2,1)*D&
&(3,3)*X(3,3)-D(2,3)*D(3,1)*X(3,3)-2*D(2,1)*D(2,3)*X(3,3)-D(1,1)*D&
&(2,3)*X(3,3)-D(1,3)*D(2,1)*X(3,3)-D(1,1)*X(3,2)*D(3,3)+2*D(3,1)*X&
&(3,1)*D(3,3)+2*D(2,1)*X(3,1)*D(3,3)+2*D(1,1)*X(3,1)*D(3,3)-D(1,3)&
&*D(3,1)*X(3,2)-D(1,1)*D(2,3)*X(3,2)-D(1,3)*D(2,1)*X(3,2)-2*D(1,1)&
&*D(1,3)*X(3,2)+2*D(2,3)*D(3,1)*X(3,1)+2*D(1,3)*D(3,1)*X(3,1)+2*D(&
&2,1)*D(2,3)*X(3,1)+2*D(1,1)*D(2,3)*X(3,1)+2*D(1,3)*D(2,1)*X(3,1)+&
&2*D(1,1)*D(1,3)*X(3,1))*tt14+1.0E+0*((-2*D(3,1)*D(3,2)*X(3,4))-D(&
&2,1)*D(3,2)*X(3,4)-D(1,1)*D(3,2)*X(3,4)-D(2,2)*D(3,1)*X(3,4)-D(1,&
&2)*D(3,1)*X(3,4)-D(2,1)*D(3,2)*X(3,3)-D(2,2)*D(3,1)*X(3,3)-2*D(2,&
&1)*D(2,2)*X(3,3)-D(1,1)*D(2,2)*X(3,3)-D(1,2)*D(2,1)*X(3,3)-D(1,1)&
&*D(3,2)*X(3,2)-D(1,2)*D(3,1)*X(3,2)-D(1,1)*D(2,2)*X(3,2)-D(1,2)*D&
&(2,1)*X(3,2)-2*D(1,1)*D(1,2)*X(3,2)+2*D(3,1)*X(3,1)*D(3,2)+2*D(2,&
&1)*X(3,1)*D(3,2)+2*D(1,1)*X(3,1)*D(3,2)+2*D(2,2)*D(3,1)*X(3,1)+2*&
&D(1,2)*D(3,1)*X(3,1)+2*D(2,1)*D(2,2)*X(3,1)+2*D(1,1)*D(2,2)*X(3,1&
&)+2*D(1,2)*D(2,1)*X(3,1)+2*D(1,1)*D(1,2)*X(3,1))*tt13))
jac(1,4) = volume(1,1)*(1.0E+0*lam(1,1)*(1.0E+0*D(1,3)*tt37+1.0E+&
&0*D(1,2)*tt36+1.0E+0*D(1,1)*tt35)*tt52+miu(1,1)*(1.0E+0*D(1,3)*tt&
&27*tt30+1.0E+0*D(1,2)*tt22*tt25+1.0E+0*D(1,1)*tt17*tt20+1.0E+0*(D&
&(1,2)*X(1,4)*D(3,3)+tt64+D(1,3)*X(1,4)*D(3,2)+tt63+D(1,2)*X(1,3)*&
&D(2,3)+tt62+D(1,3)*X(1,3)*D(2,2)+tt61+2*D(1,2)*X(1,2)*D(1,3)-2*X(&
&1,1)*D(1,2)*D(1,3))*tt15+1.0E+0*(D(1,1)*X(1,4)*D(3,3)+tt60+D(1,3)&
&*X(1,4)*D(3,1)+tt59+D(1,1)*X(1,3)*D(2,3)+tt58+D(1,3)*X(1,3)*D(2,1&
&)+tt57+2*D(1,1)*X(1,2)*D(1,3)-2*D(1,1)*X(1,1)*D(1,3))*tt14+1.0E+0&
&*(D(1,1)*X(1,4)*D(3,2)+tt56+D(1,2)*X(1,4)*D(3,1)+tt55+D(1,1)*X(1,&
&3)*D(2,2)+tt54+D(1,2)*X(1,3)*D(2,1)+tt53+2*D(1,1)*D(1,2)*X(1,2)-2&
&*D(1,1)*X(1,1)*D(1,2))*tt13))
jac(1,5) = volume(1,1)*(1.0E+0*lam(1,1)*(1.0E+0*D(1,3)*tt50+1.0E+&
&0*D(1,2)*tt48+1.0E+0*D(1,1)*tt42)*tt52+miu(1,1)*(1.0E+0*D(1,3)*tt&
&28*tt30+1.0E+0*D(1,2)*tt23*tt25+1.0E+0*D(1,1)*tt18*tt20+1.0E+0*(D&
&(1,2)*X(2,4)*D(3,3)+tt76+D(1,3)*X(2,4)*D(3,2)+tt75+D(1,2)*D(2,3)*&
&X(2,3)+D(1,3)*D(2,2)*X(2,3)+tt74+2*D(1,2)*D(1,3)*X(2,2)+tt73-2*D(&
&1,2)*D(1,3)*X(2,1))*tt15+1.0E+0*(D(1,1)*X(2,4)*D(3,3)+tt72+D(1,3)&
&*X(2,4)*D(3,1)+tt71+D(1,1)*D(2,3)*X(2,3)+D(1,3)*D(2,1)*X(2,3)+tt7&
&0+2*D(1,1)*D(1,3)*X(2,2)+tt69-2*D(1,1)*D(1,3)*X(2,1))*tt14+1.0E+0&
&*(D(1,1)*X(2,4)*D(3,2)+tt68+D(1,2)*X(2,4)*D(3,1)+tt67+D(1,1)*D(2,&
&2)*X(2,3)+D(1,2)*D(2,1)*X(2,3)+2*D(1,1)*D(1,2)*X(2,2)+tt66+tt65-2&
&*D(1,1)*D(1,2)*X(2,1))*tt13))
jac(1,6) = volume(1,1)*(1.0E+0*lam(1,1)*(1.0E+0*D(1,3)*tt51+1.0E+&
&0*D(1,2)*tt49+1.0E+0*D(1,1)*tt47)*tt52+miu(1,1)*(1.0E+0*D(1,3)*tt&
&29*tt30+1.0E+0*D(1,2)*tt24*tt25+1.0E+0*D(1,1)*tt19*tt20+1.0E+0*(D&
&(1,2)*D(3,3)*X(3,4)+D(1,3)*D(3,2)*X(3,4)+D(1,2)*D(2,3)*X(3,3)+D(1&
&,3)*D(2,2)*X(3,3)+tt88+2*D(1,2)*D(1,3)*X(3,2)+tt87+tt86+tt85-2*D(&
&1,2)*D(1,3)*X(3,1))*tt15+1.0E+0*(D(1,1)*D(3,3)*X(3,4)+D(1,3)*D(3,&
&1)*X(3,4)+D(1,1)*D(2,3)*X(3,3)+D(1,3)*D(2,1)*X(3,3)+tt84+2*D(1,1)&
&*D(1,3)*X(3,2)+tt83+tt82+tt81-2*D(1,1)*D(1,3)*X(3,1))*tt14+1.0E+0&
&*(D(1,1)*D(3,2)*X(3,4)+D(1,2)*D(3,1)*X(3,4)+D(1,1)*D(2,2)*X(3,3)+&
&D(1,2)*D(2,1)*X(3,3)+2*D(1,1)*D(1,2)*X(3,2)+tt80+tt79+tt78+tt77-2&
&*D(1,1)*D(1,2)*X(3,1))*tt13))
jac(1,7) = volume(1,1)*(1.0E+0*lam(1,1)*(1.0E+0*D(2,3)*tt37+1.0E+&
&0*D(2,2)*tt36+1.0E+0*D(2,1)*tt35)*tt52+miu(1,1)*(1.0E+0*D(2,3)*tt&
&27*tt30+1.0E+0*D(2,2)*tt22*tt25+1.0E+0*D(2,1)*tt17*tt20+1.0E+0*(X&
&(1,4)*D(2,2)*D(3,3)+tt94+X(1,4)*D(2,3)*D(3,2)+tt93+2*X(1,3)*D(2,2&
&)*D(2,3)-2*X(1,1)*D(2,2)*D(2,3)+D(1,2)*X(1,2)*D(2,3)+tt62+X(1,2)*&
&D(1,3)*D(2,2)+tt61)*tt15+1.0E+0*(X(1,4)*D(2,1)*D(3,3)+tt92+X(1,4)&
&*D(2,3)*D(3,1)+tt91+2*X(1,3)*D(2,1)*D(2,3)-2*X(1,1)*D(2,1)*D(2,3)&
&+D(1,1)*X(1,2)*D(2,3)+tt58+X(1,2)*D(1,3)*D(2,1)+tt57)*tt14+1.0E+0&
&*(X(1,4)*D(2,1)*D(3,2)+tt90+X(1,4)*D(2,2)*D(3,1)+tt89+2*X(1,3)*D(&
&2,1)*D(2,2)-2*X(1,1)*D(2,1)*D(2,2)+D(1,1)*X(1,2)*D(2,2)+tt54+D(1,&
&2)*X(1,2)*D(2,1)+tt53)*tt13))
jac(1,8) = volume(1,1)*(1.0E+0*lam(1,1)*(1.0E+0*D(2,3)*tt50+1.0E+&
&0*D(2,2)*tt48+1.0E+0*D(2,1)*tt42)*tt52+miu(1,1)*(1.0E+0*D(2,3)*tt&
&28*tt30+1.0E+0*D(2,2)*tt23*tt25+1.0E+0*D(2,1)*tt18*tt20+1.0E+0*(D&
&(2,2)*X(2,4)*D(3,3)+tt100+D(2,3)*X(2,4)*D(3,2)+tt99+2*D(2,2)*D(2,&
&3)*X(2,3)+D(1,2)*X(2,2)*D(2,3)-2*X(2,1)*D(2,2)*D(2,3)+tt74+D(1,3)&
&*D(2,2)*X(2,2)+tt73)*tt15+1.0E+0*(D(2,1)*X(2,4)*D(3,3)+tt98+D(2,3&
&)*X(2,4)*D(3,1)+tt97+2*D(2,1)*D(2,3)*X(2,3)+D(1,1)*X(2,2)*D(2,3)-&
&2*D(2,1)*X(2,1)*D(2,3)+tt70+D(1,3)*D(2,1)*X(2,2)+tt69)*tt14+1.0E+&
&0*(D(2,1)*X(2,4)*D(3,2)+tt96+D(2,2)*X(2,4)*D(3,1)+tt95+2*D(2,1)*D&
&(2,2)*X(2,3)+D(1,1)*D(2,2)*X(2,2)+D(1,2)*D(2,1)*X(2,2)-2*D(2,1)*X&
&(2,1)*D(2,2)+tt66+tt65)*tt13))
jac(1,9) = volume(1,1)*(1.0E+0*lam(1,1)*(1.0E+0*D(2,3)*tt51+1.0E+&
&0*D(2,2)*tt49+1.0E+0*D(2,1)*tt47)*tt52+miu(1,1)*(1.0E+0*D(2,3)*tt&
&29*tt30+1.0E+0*D(2,2)*tt24*tt25+1.0E+0*D(2,1)*tt19*tt20+1.0E+0*(D&
&(2,2)*D(3,3)*X(3,4)+D(2,3)*D(3,2)*X(3,4)+2*D(2,2)*D(2,3)*X(3,3)+t&
&t106+D(1,2)*D(2,3)*X(3,2)+D(1,3)*D(2,2)*X(3,2)+tt105-2*D(2,2)*D(2&
&,3)*X(3,1)+tt86+tt85)*tt15+1.0E+0*(D(2,1)*D(3,3)*X(3,4)+D(2,3)*D(&
&3,1)*X(3,4)+2*D(2,1)*D(2,3)*X(3,3)+tt104+D(1,1)*D(2,3)*X(3,2)+D(1&
&,3)*D(2,1)*X(3,2)+tt103-2*D(2,1)*D(2,3)*X(3,1)+tt82+tt81)*tt14+1.&
&0E+0*(D(2,1)*D(3,2)*X(3,4)+D(2,2)*D(3,1)*X(3,4)+2*D(2,1)*D(2,2)*X&
&(3,3)+D(1,1)*D(2,2)*X(3,2)+D(1,2)*D(2,1)*X(3,2)+tt102+tt101-2*D(2&
&,1)*D(2,2)*X(3,1)+tt78+tt77)*tt13))
jac(1,10) = volume(1,1)*(1.0E+0*lam(1,1)*(1.0E+0*D(3,3)*tt37+1.0E&
&+0*D(3,2)*tt36+1.0E+0*D(3,1)*tt35)*tt52+miu(1,1)*(1.0E+0*D(3,3)*t&
&t27*tt30+1.0E+0*D(3,2)*tt22*tt25+1.0E+0*D(3,1)*tt17*tt20+1.0E+0*(&
&2*X(1,4)*D(3,2)*D(3,3)-2*X(1,1)*D(3,2)*D(3,3)+X(1,3)*D(2,2)*D(3,3&
&)+tt94+D(1,2)*X(1,2)*D(3,3)+tt64+X(1,3)*D(2,3)*D(3,2)+tt93+X(1,2)&
&*D(1,3)*D(3,2)+tt63)*tt15+1.0E+0*(2*X(1,4)*D(3,1)*D(3,3)-2*X(1,1)&
&*D(3,1)*D(3,3)+X(1,3)*D(2,1)*D(3,3)+tt92+D(1,1)*X(1,2)*D(3,3)+tt6&
&0+X(1,3)*D(2,3)*D(3,1)+tt91+X(1,2)*D(1,3)*D(3,1)+tt59)*tt14+1.0E+&
&0*(2*X(1,4)*D(3,1)*D(3,2)-2*X(1,1)*D(3,1)*D(3,2)+X(1,3)*D(2,1)*D(&
&3,2)+tt90+D(1,1)*X(1,2)*D(3,2)+tt56+X(1,3)*D(2,2)*D(3,1)+tt89+D(1&
&,2)*X(1,2)*D(3,1)+tt55)*tt13))
jac(1,11) = volume(1,1)*(1.0E+0*lam(1,1)*(1.0E+0*D(3,3)*tt50+1.0E&
&+0*D(3,2)*tt48+1.0E+0*D(3,1)*tt42)*tt52+miu(1,1)*(1.0E+0*D(3,3)*t&
&t28*tt30+1.0E+0*D(3,2)*tt23*tt25+1.0E+0*D(3,1)*tt18*tt20+1.0E+0*(&
&2*X(2,4)*D(3,2)*D(3,3)-2*X(2,1)*D(3,2)*D(3,3)+D(2,2)*X(2,3)*D(3,3&
&)+D(1,2)*X(2,2)*D(3,3)+tt100+tt76+D(2,3)*X(2,3)*D(3,2)+tt99+D(1,3&
&)*X(2,2)*D(3,2)+tt75)*tt15+1.0E+0*(2*X(2,4)*D(3,1)*D(3,3)-2*X(2,1&
&)*D(3,1)*D(3,3)+D(2,1)*X(2,3)*D(3,3)+D(1,1)*X(2,2)*D(3,3)+tt98+tt&
&72+D(2,3)*X(2,3)*D(3,1)+tt97+D(1,3)*X(2,2)*D(3,1)+tt71)*tt14+1.0E&
&+0*(2*X(2,4)*D(3,1)*D(3,2)-2*X(2,1)*D(3,1)*D(3,2)+D(2,1)*X(2,3)*D&
&(3,2)+D(1,1)*X(2,2)*D(3,2)+tt96+tt68+D(2,2)*X(2,3)*D(3,1)+D(1,2)*&
&X(2,2)*D(3,1)+tt95+tt67)*tt13))
jac(1,12) = volume(1,1)*(1.0E+0*lam(1,1)*(1.0E+0*D(3,3)*tt51+1.0E&
&+0*D(3,2)*tt49+1.0E+0*D(3,1)*tt47)*tt52+miu(1,1)*(1.0E+0*D(3,3)*t&
&t29*tt30+1.0E+0*D(3,2)*tt24*tt25+1.0E+0*D(3,1)*tt19*tt20+1.0E+0*(&
&2*D(3,2)*D(3,3)*X(3,4)+D(2,2)*D(3,3)*X(3,3)+D(2,3)*D(3,2)*X(3,3)+&
&D(1,2)*X(3,2)*D(3,3)-2*X(3,1)*D(3,2)*D(3,3)+tt106+tt88+D(1,3)*D(3&
&,2)*X(3,2)+tt105+tt87)*tt15+1.0E+0*(2*D(3,1)*D(3,3)*X(3,4)+D(2,1)&
&*D(3,3)*X(3,3)+D(2,3)*D(3,1)*X(3,3)+D(1,1)*X(3,2)*D(3,3)-2*D(3,1)&
&*X(3,1)*D(3,3)+tt104+tt84+D(1,3)*D(3,1)*X(3,2)+tt103+tt83)*tt14+1&
&.0E+0*(2*D(3,1)*D(3,2)*X(3,4)+D(2,1)*D(3,2)*X(3,3)+D(2,2)*D(3,1)*&
&X(3,3)+D(1,1)*D(3,2)*X(3,2)+D(1,2)*D(3,1)*X(3,2)-2*D(3,1)*X(3,1)*&
&D(3,2)+tt102+tt80+tt101+tt79)*tt13))
END 
SUBROUTINE tet_stvk_hes(hes, X, D, volume, lam, miu) 
IMPLICIT NONE 
REAL(KIND=8) hes(12, 12) 
REAL(KIND=8) X(3, 4) 
REAL(KIND=8) D(3, 3) 
REAL(KIND=8) volume(1, 1) 
REAL(KIND=8) lam(1, 1) 
REAL(KIND=8) miu(1, 1) 
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
tt1 = (-D(3,1))-D(2,1)-D(1,1)
tt2 = -X(1,1)
tt3 = X(1,2)+tt2
tt4 = X(1,3)+tt2
tt5 = X(1,4)+tt2
tt6 = tt5*D(3,1)+tt4*D(2,1)+D(1,1)*tt3
tt7 = (-D(3,2))-D(2,2)-D(1,2)
tt8 = tt5*D(3,2)+tt4*D(2,2)+D(1,2)*tt3
tt9 = (-D(3,3))-D(2,3)-D(1,3)
tt10 = tt5*D(3,3)+tt4*D(2,3)+tt3*D(1,3)
tt11 = 1.0E+0*tt9*tt10+1.0E+0*tt7*tt8+1.0E+0*tt1*tt6
tt12 = tt1**2
tt13 = X(1,4)*D(3,1)-X(1,1)*D(3,1)+X(1,3)*D(2,1)-X(1,1)*D(2,1)+D(&
&1,1)*X(1,2)-D(1,1)*X(1,1)
tt14 = tt13**2
tt15 = tt7**2
tt16 = X(1,4)*D(3,2)-X(1,1)*D(3,2)+X(1,3)*D(2,2)-X(1,1)*D(2,2)+D(&
&1,2)*X(1,2)-X(1,1)*D(1,2)
tt17 = tt16**2
tt18 = (-2*X(1,4)*D(3,1)*D(3,2))+2*X(1,1)*D(3,1)*D(3,2)-X(1,4)*D(&
&2,1)*D(3,2)-X(1,3)*D(2,1)*D(3,2)+2*X(1,1)*D(2,1)*D(3,2)-D(1,1)*X(&
&1,4)*D(3,2)-D(1,1)*X(1,2)*D(3,2)+2*D(1,1)*X(1,1)*D(3,2)-X(1,4)*D(&
&2,2)*D(3,1)-X(1,3)*D(2,2)*D(3,1)+2*X(1,1)*D(2,2)*D(3,1)-D(1,2)*X(&
&1,4)*D(3,1)-D(1,2)*X(1,2)*D(3,1)+2*X(1,1)*D(1,2)*D(3,1)-2*X(1,3)*&
&D(2,1)*D(2,2)+2*X(1,1)*D(2,1)*D(2,2)-D(1,1)*X(1,3)*D(2,2)-D(1,1)*&
&X(1,2)*D(2,2)+2*D(1,1)*X(1,1)*D(2,2)-D(1,2)*X(1,3)*D(2,1)-D(1,2)*&
&X(1,2)*D(2,1)+2*X(1,1)*D(1,2)*D(2,1)-2*D(1,1)*D(1,2)*X(1,2)+2*D(1&
&,1)*X(1,1)*D(1,2)
tt19 = tt9**2
tt20 = X(1,4)*D(3,3)-X(1,1)*D(3,3)+X(1,3)*D(2,3)-X(1,1)*D(2,3)+X(&
&1,2)*D(1,3)-X(1,1)*D(1,3)
tt21 = tt20**2
tt22 = (-2*X(1,4)*D(3,1)*D(3,3))+2*X(1,1)*D(3,1)*D(3,3)-X(1,4)*D(&
&2,1)*D(3,3)-X(1,3)*D(2,1)*D(3,3)+2*X(1,1)*D(2,1)*D(3,3)-D(1,1)*X(&
&1,4)*D(3,3)-D(1,1)*X(1,2)*D(3,3)+2*D(1,1)*X(1,1)*D(3,3)-X(1,4)*D(&
&2,3)*D(3,1)-X(1,3)*D(2,3)*D(3,1)+2*X(1,1)*D(2,3)*D(3,1)-D(1,3)*X(&
&1,4)*D(3,1)-X(1,2)*D(1,3)*D(3,1)+2*X(1,1)*D(1,3)*D(3,1)-2*X(1,3)*&
&D(2,1)*D(2,3)+2*X(1,1)*D(2,1)*D(2,3)-D(1,1)*X(1,3)*D(2,3)-D(1,1)*&
&X(1,2)*D(2,3)+2*D(1,1)*X(1,1)*D(2,3)-D(1,3)*X(1,3)*D(2,1)-X(1,2)*&
&D(1,3)*D(2,1)+2*X(1,1)*D(1,3)*D(2,1)-2*D(1,1)*X(1,2)*D(1,3)+2*D(1&
&,1)*X(1,1)*D(1,3)
tt23 = (-2*X(1,4)*D(3,2)*D(3,3))+2*X(1,1)*D(3,2)*D(3,3)-X(1,4)*D(&
&2,2)*D(3,3)-X(1,3)*D(2,2)*D(3,3)+2*X(1,1)*D(2,2)*D(3,3)-D(1,2)*X(&
&1,4)*D(3,3)-D(1,2)*X(1,2)*D(3,3)+2*X(1,1)*D(1,2)*D(3,3)-X(1,4)*D(&
&2,3)*D(3,2)-X(1,3)*D(2,3)*D(3,2)+2*X(1,1)*D(2,3)*D(3,2)-D(1,3)*X(&
&1,4)*D(3,2)-X(1,2)*D(1,3)*D(3,2)+2*X(1,1)*D(1,3)*D(3,2)-2*X(1,3)*&
&D(2,2)*D(2,3)+2*X(1,1)*D(2,2)*D(2,3)-D(1,2)*X(1,3)*D(2,3)-D(1,2)*&
&X(1,2)*D(2,3)+2*X(1,1)*D(1,2)*D(2,3)-D(1,3)*X(1,3)*D(2,2)-X(1,2)*&
&D(1,3)*D(2,2)+2*X(1,1)*D(1,3)*D(2,2)-2*D(1,2)*X(1,2)*D(1,3)+2*X(1&
&,1)*D(1,2)*D(1,3)
tt24 = X(1,1)**2
tt25 = X(1,2)**2
tt26 = X(2,1)**2
tt27 = X(1,3)**2
tt28 = X(2,2)**2
tt29 = X(2,3)**2
tt30 = X(3,1)**2
tt31 = X(1,4)**2
tt32 = X(2,4)**2
tt33 = X(3,2)**2
tt34 = X(3,3)**2
tt35 = X(3,4)**2
tt36 = D(3,1)*D(3,2)*tt35+D(2,1)*D(3,2)*X(3,3)*X(3,4)+D(2,2)*D(3,&
&1)*X(3,3)*X(3,4)+D(1,1)*D(3,2)*X(3,2)*X(3,4)+D(1,2)*D(3,1)*X(3,2)&
&*X(3,4)-2*D(3,1)*X(3,1)*D(3,2)*X(3,4)-D(2,1)*X(3,1)*D(3,2)*X(3,4)&
&-D(1,1)*X(3,1)*D(3,2)*X(3,4)-D(2,2)*D(3,1)*X(3,1)*X(3,4)-D(1,2)*D&
&(3,1)*X(3,1)*X(3,4)+D(2,1)*D(2,2)*tt34+D(1,1)*D(2,2)*X(3,2)*X(3,3&
&)+D(1,2)*D(2,1)*X(3,2)*X(3,3)-D(2,1)*X(3,1)*D(3,2)*X(3,3)-D(2,2)*&
&D(3,1)*X(3,1)*X(3,3)-2*D(2,1)*D(2,2)*X(3,1)*X(3,3)-D(1,1)*D(2,2)*&
&X(3,1)*X(3,3)-D(1,2)*D(2,1)*X(3,1)*X(3,3)+D(1,1)*D(1,2)*tt33-D(1,&
&1)*X(3,1)*D(3,2)*X(3,2)-D(1,2)*D(3,1)*X(3,1)*X(3,2)-D(1,1)*D(2,2)&
&*X(3,1)*X(3,2)-D(1,2)*D(2,1)*X(3,1)*X(3,2)-2*D(1,1)*D(1,2)*X(3,1)&
&*X(3,2)+D(3,1)*tt30*D(3,2)+D(2,1)*tt30*D(3,2)+D(1,1)*tt30*D(3,2)+&
&tt32*D(3,1)*D(3,2)-2*X(2,1)*X(2,4)*D(3,1)*D(3,2)+tt26*D(3,1)*D(3,&
&2)+tt31*D(3,1)*D(3,2)-2*X(1,1)*X(1,4)*D(3,1)*D(3,2)+tt24*D(3,1)*D&
&(3,2)+D(2,1)*X(2,3)*X(2,4)*D(3,2)+D(1,1)*X(2,2)*X(2,4)*D(3,2)-D(2&
&,1)*X(2,1)*X(2,4)*D(3,2)-D(1,1)*X(2,1)*X(2,4)*D(3,2)-D(2,1)*X(2,1&
&)*X(2,3)*D(3,2)-D(1,1)*X(2,1)*X(2,2)*D(3,2)+D(2,1)*tt26*D(3,2)+D(&
&1,1)*tt26*D(3,2)+X(1,3)*X(1,4)*D(2,1)*D(3,2)-X(1,1)*X(1,4)*D(2,1)&
&*D(3,2)-X(1,1)*X(1,3)*D(2,1)*D(3,2)+tt24*D(2,1)*D(3,2)+D(1,1)*X(1&
&,2)*X(1,4)*D(3,2)-D(1,1)*X(1,1)*X(1,4)*D(3,2)-D(1,1)*X(1,1)*X(1,2&
&)*D(3,2)+D(1,1)*tt24*D(3,2)+D(2,2)*D(3,1)*tt30+D(1,2)*D(3,1)*tt30&
&+D(2,1)*D(2,2)*tt30+D(1,1)*D(2,2)*tt30+D(1,2)*D(2,1)*tt30+D(1,1)*&
&D(1,2)*tt30+D(2,2)*X(2,3)*X(2,4)*D(3,1)+D(1,2)*X(2,2)*X(2,4)*D(3,&
&1)-X(2,1)*D(2,2)*X(2,4)*D(3,1)-D(1,2)*X(2,1)*X(2,4)*D(3,1)-X(2,1)&
&*D(2,2)*X(2,3)*D(3,1)-D(1,2)*X(2,1)*X(2,2)*D(3,1)+tt26*D(2,2)*D(3&
&,1)+X(1,3)*X(1,4)*D(2,2)*D(3,1)-X(1,1)*X(1,4)*D(2,2)*D(3,1)-X(1,1&
&)*X(1,3)*D(2,2)*D(3,1)+tt24*D(2,2)*D(3,1)+D(1,2)*tt26*D(3,1)+D(1,&
&2)*X(1,2)*X(1,4)*D(3,1)-X(1,1)*D(1,2)*X(1,4)*D(3,1)-X(1,1)*D(1,2)&
&*X(1,2)*D(3,1)+tt24*D(1,2)*D(3,1)+D(2,1)*D(2,2)*tt29+D(1,1)*D(2,2&
&)*X(2,2)*X(2,3)+D(1,2)*D(2,1)*X(2,2)*X(2,3)-2*D(2,1)*X(2,1)*D(2,2&
&)*X(2,3)-D(1,1)*X(2,1)*D(2,2)*X(2,3)-D(1,2)*D(2,1)*X(2,1)*X(2,3)+&
&D(1,1)*D(1,2)*tt28-D(1,1)*X(2,1)*D(2,2)*X(2,2)-D(1,2)*D(2,1)*X(2,&
&1)*X(2,2)-2*D(1,1)*D(1,2)*X(2,1)*X(2,2)+D(2,1)*tt26*D(2,2)+D(1,1)&
&*tt26*D(2,2)+tt27*D(2,1)*D(2,2)-2*X(1,1)*X(1,3)*D(2,1)*D(2,2)+tt2&
&4*D(2,1)*D(2,2)+D(1,1)*X(1,2)*X(1,3)*D(2,2)-D(1,1)*X(1,1)*X(1,3)*&
&D(2,2)-D(1,1)*X(1,1)*X(1,2)*D(2,2)+D(1,1)*tt24*D(2,2)+D(1,2)*D(2,&
&1)*tt26+D(1,1)*D(1,2)*tt26+D(1,2)*X(1,2)*X(1,3)*D(2,1)-X(1,1)*D(1&
&,2)*X(1,3)*D(2,1)-X(1,1)*D(1,2)*X(1,2)*D(2,1)+tt24*D(1,2)*D(2,1)+&
&D(1,1)*D(1,2)*tt25-2*D(1,1)*X(1,1)*D(1,2)*X(1,2)+D(1,1)*tt24*D(1,&
&2)
tt37 = 1.0E+0*(2*D(3,1)*D(3,2)+2*D(2,1)*D(3,2)+2*D(1,1)*D(3,2)+2*&
&D(2,2)*D(3,1)+2*D(1,2)*D(3,1)+2*D(2,1)*D(2,2)+2*D(1,1)*D(2,2)+2*D&
&(1,2)*D(2,1)+2*D(1,1)*D(1,2))*tt36
tt38 = D(3,1)*D(3,3)*tt35+D(2,1)*D(3,3)*X(3,3)*X(3,4)+D(2,3)*D(3,&
&1)*X(3,3)*X(3,4)+D(1,1)*X(3,2)*D(3,3)*X(3,4)-2*D(3,1)*X(3,1)*D(3,&
&3)*X(3,4)-D(2,1)*X(3,1)*D(3,3)*X(3,4)-D(1,1)*X(3,1)*D(3,3)*X(3,4)&
&+D(1,3)*D(3,1)*X(3,2)*X(3,4)-D(2,3)*D(3,1)*X(3,1)*X(3,4)-D(1,3)*D&
&(3,1)*X(3,1)*X(3,4)+D(2,1)*D(2,3)*tt34-D(2,1)*X(3,1)*D(3,3)*X(3,3&
&)+D(1,1)*D(2,3)*X(3,2)*X(3,3)+D(1,3)*D(2,1)*X(3,2)*X(3,3)-D(2,3)*&
&D(3,1)*X(3,1)*X(3,3)-2*D(2,1)*D(2,3)*X(3,1)*X(3,3)-D(1,1)*D(2,3)*&
&X(3,1)*X(3,3)-D(1,3)*D(2,1)*X(3,1)*X(3,3)-D(1,1)*X(3,1)*X(3,2)*D(&
&3,3)+D(3,1)*tt30*D(3,3)+D(2,1)*tt30*D(3,3)+D(1,1)*tt30*D(3,3)+tt3&
&2*D(3,1)*D(3,3)-2*X(2,1)*X(2,4)*D(3,1)*D(3,3)+tt26*D(3,1)*D(3,3)+&
&tt31*D(3,1)*D(3,3)-2*X(1,1)*X(1,4)*D(3,1)*D(3,3)+tt24*D(3,1)*D(3,&
&3)+D(2,1)*X(2,3)*X(2,4)*D(3,3)+D(1,1)*X(2,2)*X(2,4)*D(3,3)-D(2,1)&
&*X(2,1)*X(2,4)*D(3,3)-D(1,1)*X(2,1)*X(2,4)*D(3,3)-D(2,1)*X(2,1)*X&
&(2,3)*D(3,3)-D(1,1)*X(2,1)*X(2,2)*D(3,3)+D(2,1)*tt26*D(3,3)+D(1,1&
&)*tt26*D(3,3)+X(1,3)*X(1,4)*D(2,1)*D(3,3)-X(1,1)*X(1,4)*D(2,1)*D(&
&3,3)-X(1,1)*X(1,3)*D(2,1)*D(3,3)+tt24*D(2,1)*D(3,3)+D(1,1)*X(1,2)&
&*X(1,4)*D(3,3)-D(1,1)*X(1,1)*X(1,4)*D(3,3)-D(1,1)*X(1,1)*X(1,2)*D&
&(3,3)+D(1,1)*tt24*D(3,3)+D(1,1)*D(1,3)*tt33-D(1,3)*D(3,1)*X(3,1)*&
&X(3,2)-D(1,1)*D(2,3)*X(3,1)*X(3,2)-D(1,3)*D(2,1)*X(3,1)*X(3,2)-2*&
&D(1,1)*D(1,3)*X(3,1)*X(3,2)+D(2,3)*D(3,1)*tt30+D(1,3)*D(3,1)*tt30&
&+D(2,1)*D(2,3)*tt30+D(1,1)*D(2,3)*tt30+D(1,3)*D(2,1)*tt30+D(1,1)*&
&D(1,3)*tt30+D(2,3)*X(2,3)*X(2,4)*D(3,1)-X(2,1)*D(2,3)*X(2,4)*D(3,&
&1)+D(1,3)*X(2,2)*X(2,4)*D(3,1)-D(1,3)*X(2,1)*X(2,4)*D(3,1)-X(2,1)&
&*D(2,3)*X(2,3)*D(3,1)+tt26*D(2,3)*D(3,1)+X(1,3)*X(1,4)*D(2,3)*D(3&
&,1)-X(1,1)*X(1,4)*D(2,3)*D(3,1)-X(1,1)*X(1,3)*D(2,3)*D(3,1)+tt24*&
&D(2,3)*D(3,1)-D(1,3)*X(2,1)*X(2,2)*D(3,1)+D(1,3)*tt26*D(3,1)+X(1,&
&2)*D(1,3)*X(1,4)*D(3,1)-X(1,1)*D(1,3)*X(1,4)*D(3,1)-X(1,1)*X(1,2)&
&*D(1,3)*D(3,1)+tt24*D(1,3)*D(3,1)+D(2,1)*D(2,3)*tt29+D(1,1)*X(2,2&
&)*D(2,3)*X(2,3)-2*D(2,1)*X(2,1)*D(2,3)*X(2,3)-D(1,1)*X(2,1)*D(2,3&
&)*X(2,3)+D(1,3)*D(2,1)*X(2,2)*X(2,3)-D(1,3)*D(2,1)*X(2,1)*X(2,3)-&
&D(1,1)*X(2,1)*X(2,2)*D(2,3)+D(2,1)*tt26*D(2,3)+D(1,1)*tt26*D(2,3)&
&+tt27*D(2,1)*D(2,3)-2*X(1,1)*X(1,3)*D(2,1)*D(2,3)+tt24*D(2,1)*D(2&
&,3)+D(1,1)*X(1,2)*X(1,3)*D(2,3)-D(1,1)*X(1,1)*X(1,3)*D(2,3)-D(1,1&
&)*X(1,1)*X(1,2)*D(2,3)+D(1,1)*tt24*D(2,3)+D(1,1)*D(1,3)*tt28-D(1,&
&3)*D(2,1)*X(2,1)*X(2,2)-2*D(1,1)*D(1,3)*X(2,1)*X(2,2)+D(1,3)*D(2,&
&1)*tt26+D(1,1)*D(1,3)*tt26+X(1,2)*D(1,3)*X(1,3)*D(2,1)-X(1,1)*D(1&
&,3)*X(1,3)*D(2,1)-X(1,1)*X(1,2)*D(1,3)*D(2,1)+tt24*D(1,3)*D(2,1)+&
&D(1,1)*tt25*D(1,3)-2*D(1,1)*X(1,1)*X(1,2)*D(1,3)+D(1,1)*tt24*D(1,&
&3)
tt39 = 1.0E+0*(2*D(3,1)*D(3,3)+2*D(2,1)*D(3,3)+2*D(1,1)*D(3,3)+2*&
&D(2,3)*D(3,1)+2*D(1,3)*D(3,1)+2*D(2,1)*D(2,3)+2*D(1,1)*D(2,3)+2*D&
&(1,3)*D(2,1)+2*D(1,1)*D(1,3))*tt38
tt40 = D(3,2)*D(3,3)*tt35+D(2,2)*D(3,3)*X(3,3)*X(3,4)+D(2,3)*D(3,&
&2)*X(3,3)*X(3,4)+D(1,2)*X(3,2)*D(3,3)*X(3,4)-2*X(3,1)*D(3,2)*D(3,&
&3)*X(3,4)-D(2,2)*X(3,1)*D(3,3)*X(3,4)-D(1,2)*X(3,1)*D(3,3)*X(3,4)&
&+D(1,3)*D(3,2)*X(3,2)*X(3,4)-D(2,3)*X(3,1)*D(3,2)*X(3,4)-D(1,3)*X&
&(3,1)*D(3,2)*X(3,4)+D(2,2)*D(2,3)*tt34-D(2,2)*X(3,1)*D(3,3)*X(3,3&
&)+D(1,2)*D(2,3)*X(3,2)*X(3,3)+D(1,3)*D(2,2)*X(3,2)*X(3,3)-D(2,3)*&
&X(3,1)*D(3,2)*X(3,3)-2*D(2,2)*D(2,3)*X(3,1)*X(3,3)-D(1,2)*D(2,3)*&
&X(3,1)*X(3,3)-D(1,3)*D(2,2)*X(3,1)*X(3,3)-D(1,2)*X(3,1)*X(3,2)*D(&
&3,3)+tt30*D(3,2)*D(3,3)+tt32*D(3,2)*D(3,3)-2*X(2,1)*X(2,4)*D(3,2)&
&*D(3,3)+tt26*D(3,2)*D(3,3)+tt31*D(3,2)*D(3,3)-2*X(1,1)*X(1,4)*D(3&
&,2)*D(3,3)+tt24*D(3,2)*D(3,3)+D(2,2)*tt30*D(3,3)+D(1,2)*tt30*D(3,&
&3)+D(2,2)*X(2,3)*X(2,4)*D(3,3)+D(1,2)*X(2,2)*X(2,4)*D(3,3)-X(2,1)&
&*D(2,2)*X(2,4)*D(3,3)-D(1,2)*X(2,1)*X(2,4)*D(3,3)-X(2,1)*D(2,2)*X&
&(2,3)*D(3,3)-D(1,2)*X(2,1)*X(2,2)*D(3,3)+tt26*D(2,2)*D(3,3)+X(1,3&
&)*X(1,4)*D(2,2)*D(3,3)-X(1,1)*X(1,4)*D(2,2)*D(3,3)-X(1,1)*X(1,3)*&
&D(2,2)*D(3,3)+tt24*D(2,2)*D(3,3)+D(1,2)*tt26*D(3,3)+D(1,2)*X(1,2)&
&*X(1,4)*D(3,3)-X(1,1)*D(1,2)*X(1,4)*D(3,3)-X(1,1)*D(1,2)*X(1,2)*D&
&(3,3)+tt24*D(1,2)*D(3,3)+D(1,2)*D(1,3)*tt33-D(1,3)*X(3,1)*D(3,2)*&
&X(3,2)-D(1,2)*D(2,3)*X(3,1)*X(3,2)-D(1,3)*D(2,2)*X(3,1)*X(3,2)-2*&
&D(1,2)*D(1,3)*X(3,1)*X(3,2)+D(2,3)*tt30*D(3,2)+D(1,3)*tt30*D(3,2)&
&+D(2,3)*X(2,3)*X(2,4)*D(3,2)-X(2,1)*D(2,3)*X(2,4)*D(3,2)+D(1,3)*X&
&(2,2)*X(2,4)*D(3,2)-D(1,3)*X(2,1)*X(2,4)*D(3,2)-X(2,1)*D(2,3)*X(2&
&,3)*D(3,2)+tt26*D(2,3)*D(3,2)+X(1,3)*X(1,4)*D(2,3)*D(3,2)-X(1,1)*&
&X(1,4)*D(2,3)*D(3,2)-X(1,1)*X(1,3)*D(2,3)*D(3,2)+tt24*D(2,3)*D(3,&
&2)-D(1,3)*X(2,1)*X(2,2)*D(3,2)+D(1,3)*tt26*D(3,2)+X(1,2)*D(1,3)*X&
&(1,4)*D(3,2)-X(1,1)*D(1,3)*X(1,4)*D(3,2)-X(1,1)*X(1,2)*D(1,3)*D(3&
&,2)+tt24*D(1,3)*D(3,2)+D(2,2)*D(2,3)*tt30+D(1,2)*D(2,3)*tt30+D(1,&
&3)*D(2,2)*tt30+D(1,2)*D(1,3)*tt30+D(2,2)*D(2,3)*tt29+D(1,2)*X(2,2&
&)*D(2,3)*X(2,3)-2*X(2,1)*D(2,2)*D(2,3)*X(2,3)-D(1,2)*X(2,1)*D(2,3&
&)*X(2,3)+D(1,3)*D(2,2)*X(2,2)*X(2,3)-D(1,3)*X(2,1)*D(2,2)*X(2,3)-&
&D(1,2)*X(2,1)*X(2,2)*D(2,3)+tt26*D(2,2)*D(2,3)+tt27*D(2,2)*D(2,3)&
&-2*X(1,1)*X(1,3)*D(2,2)*D(2,3)+tt24*D(2,2)*D(2,3)+D(1,2)*tt26*D(2&
&,3)+D(1,2)*X(1,2)*X(1,3)*D(2,3)-X(1,1)*D(1,2)*X(1,3)*D(2,3)-X(1,1&
&)*D(1,2)*X(1,2)*D(2,3)+tt24*D(1,2)*D(2,3)+D(1,2)*D(1,3)*tt28-D(1,&
&3)*X(2,1)*D(2,2)*X(2,2)-2*D(1,2)*D(1,3)*X(2,1)*X(2,2)+D(1,3)*tt26&
&*D(2,2)+X(1,2)*D(1,3)*X(1,3)*D(2,2)-X(1,1)*D(1,3)*X(1,3)*D(2,2)-X&
&(1,1)*X(1,2)*D(1,3)*D(2,2)+tt24*D(1,3)*D(2,2)+D(1,2)*D(1,3)*tt26+&
&D(1,2)*tt25*D(1,3)-2*X(1,1)*D(1,2)*X(1,2)*D(1,3)+tt24*D(1,2)*D(1,&
&3)
tt41 = 1.0E+0*(2*D(3,2)*D(3,3)+2*D(2,2)*D(3,3)+2*D(1,2)*D(3,3)+2*&
&D(2,3)*D(3,2)+2*D(1,3)*D(3,2)+2*D(2,2)*D(2,3)+2*D(1,2)*D(2,3)+2*D&
&(1,3)*D(2,2)+2*D(1,2)*D(1,3))*tt40
tt42 = X(2,4)*D(3,1)-X(2,1)*D(3,1)+D(2,1)*X(2,3)+D(1,1)*X(2,2)-D(&
&2,1)*X(2,1)-D(1,1)*X(2,1)
tt43 = tt42**2
tt44 = D(3,1)*X(3,4)+D(2,1)*X(3,3)+D(1,1)*X(3,2)-D(3,1)*X(3,1)-D(&
&2,1)*X(3,1)-D(1,1)*X(3,1)
tt45 = tt44**2
tt46 = tt45+tt43+tt14-1
tt47 = 1.0E+0*tt12*tt46
tt48 = X(2,4)*D(3,2)-X(2,1)*D(3,2)+D(2,2)*X(2,3)+D(1,2)*X(2,2)-X(&
&2,1)*D(2,2)-D(1,2)*X(2,1)
tt49 = tt48**2
tt50 = D(3,2)*X(3,4)+D(2,2)*X(3,3)+D(1,2)*X(3,2)-X(3,1)*D(3,2)-D(&
&2,2)*X(3,1)-D(1,2)*X(3,1)
tt51 = tt50**2
tt52 = tt51+tt49+tt17-1
tt53 = 1.0E+0*tt15*tt52
tt54 = X(2,4)*D(3,3)-X(2,1)*D(3,3)+D(2,3)*X(2,3)-X(2,1)*D(2,3)+D(&
&1,3)*X(2,2)-D(1,3)*X(2,1)
tt55 = tt54**2
tt56 = D(3,3)*X(3,4)+D(2,3)*X(3,3)-X(3,1)*D(3,3)+D(1,3)*X(3,2)-D(&
&2,3)*X(3,1)-D(1,3)*X(3,1)
tt57 = tt56**2
tt58 = tt57+tt55+tt21-1
tt59 = 1.0E+0*tt19*tt58
tt60 = -X(2,1)
tt61 = X(2,2)+tt60
tt62 = X(2,3)+tt60
tt63 = X(2,4)+tt60
tt64 = tt63*D(3,1)+D(2,1)*tt62+D(1,1)*tt61
tt65 = -X(3,1)
tt66 = X(3,2)+tt65
tt67 = X(3,3)+tt65
tt68 = X(3,4)+tt65
tt69 = D(3,1)*tt68+D(2,1)*tt67+D(1,1)*tt66
tt70 = tt63*D(3,2)+D(2,2)*tt62+D(1,2)*tt61
tt71 = D(3,2)*tt68+D(2,2)*tt67+D(1,2)*tt66
tt72 = tt63*D(3,3)+D(2,3)*tt62+D(1,3)*tt61
tt73 = D(3,3)*tt68+D(2,3)*tt67+D(1,3)*tt66
tt74 = 5.0E-1*(tt73**2+tt72**2+tt10**2-1)+5.0E-1*(tt71**2+tt70**2&
&+tt8**2-1)+5.0E-1*(tt69**2+tt64**2+tt6**2-1)
tt75 = 1.0E+0*lam(1,1)*(1.0E+0*tt19+1.0E+0*tt15+1.0E+0*tt12)*tt74&
&
tt76 = 1.0E+0*tt9*tt72+1.0E+0*tt7*tt70+1.0E+0*tt1*tt64
tt77 = (-2*X(2,4)*D(3,1)*D(3,2))+2*X(2,1)*D(3,1)*D(3,2)-D(2,1)*X(&
&2,4)*D(3,2)-D(1,1)*X(2,4)*D(3,2)-D(2,1)*X(2,3)*D(3,2)-D(1,1)*X(2,&
&2)*D(3,2)+2*D(2,1)*X(2,1)*D(3,2)+2*D(1,1)*X(2,1)*D(3,2)-D(2,2)*X(&
&2,4)*D(3,1)-D(1,2)*X(2,4)*D(3,1)-D(2,2)*X(2,3)*D(3,1)-D(1,2)*X(2,&
&2)*D(3,1)+2*X(2,1)*D(2,2)*D(3,1)+2*D(1,2)*X(2,1)*D(3,1)-2*D(2,1)*&
&D(2,2)*X(2,3)-D(1,1)*D(2,2)*X(2,3)-D(1,2)*D(2,1)*X(2,3)-D(1,1)*D(&
&2,2)*X(2,2)-D(1,2)*D(2,1)*X(2,2)-2*D(1,1)*D(1,2)*X(2,2)+2*D(2,1)*&
&X(2,1)*D(2,2)+2*D(1,1)*X(2,1)*D(2,2)+2*D(1,2)*D(2,1)*X(2,1)+2*D(1&
&,1)*D(1,2)*X(2,1)
tt78 = (-2*X(2,4)*D(3,1)*D(3,3))+2*X(2,1)*D(3,1)*D(3,3)-D(2,1)*X(&
&2,4)*D(3,3)-D(1,1)*X(2,4)*D(3,3)-D(2,1)*X(2,3)*D(3,3)-D(1,1)*X(2,&
&2)*D(3,3)+2*D(2,1)*X(2,1)*D(3,3)+2*D(1,1)*X(2,1)*D(3,3)-D(2,3)*X(&
&2,4)*D(3,1)-D(1,3)*X(2,4)*D(3,1)-D(2,3)*X(2,3)*D(3,1)+2*X(2,1)*D(&
&2,3)*D(3,1)-D(1,3)*X(2,2)*D(3,1)+2*D(1,3)*X(2,1)*D(3,1)-2*D(2,1)*&
&D(2,3)*X(2,3)-D(1,1)*D(2,3)*X(2,3)-D(1,3)*D(2,1)*X(2,3)-D(1,1)*X(&
&2,2)*D(2,3)+2*D(2,1)*X(2,1)*D(2,3)+2*D(1,1)*X(2,1)*D(2,3)-D(1,3)*&
&D(2,1)*X(2,2)-2*D(1,1)*D(1,3)*X(2,2)+2*D(1,3)*D(2,1)*X(2,1)+2*D(1&
&,1)*D(1,3)*X(2,1)
tt79 = (-2*X(2,4)*D(3,2)*D(3,3))+2*X(2,1)*D(3,2)*D(3,3)-D(2,2)*X(&
&2,4)*D(3,3)-D(1,2)*X(2,4)*D(3,3)-D(2,2)*X(2,3)*D(3,3)-D(1,2)*X(2,&
&2)*D(3,3)+2*X(2,1)*D(2,2)*D(3,3)+2*D(1,2)*X(2,1)*D(3,3)-D(2,3)*X(&
&2,4)*D(3,2)-D(1,3)*X(2,4)*D(3,2)-D(2,3)*X(2,3)*D(3,2)+2*X(2,1)*D(&
&2,3)*D(3,2)-D(1,3)*X(2,2)*D(3,2)+2*D(1,3)*X(2,1)*D(3,2)-2*D(2,2)*&
&D(2,3)*X(2,3)-D(1,2)*D(2,3)*X(2,3)-D(1,3)*D(2,2)*X(2,3)-D(1,2)*X(&
&2,2)*D(2,3)+2*X(2,1)*D(2,2)*D(2,3)+2*D(1,2)*X(2,1)*D(2,3)-D(1,3)*&
&D(2,2)*X(2,2)-2*D(1,2)*D(1,3)*X(2,2)+2*D(1,3)*X(2,1)*D(2,2)+2*D(1&
&,2)*D(1,3)*X(2,1)
tt80 = volume(1,1)*(miu(1,1)*(1.0E+0*tt23*tt79+1.0E+0*tt22*tt78+2&
&.0E+0*tt19*tt20*tt54+1.0E+0*tt18*tt77+2.0E+0*tt15*tt16*tt48+2.0E+&
&0*tt12*tt13*tt42)+1.0E+0*lam(1,1)*tt11*tt76)
tt81 = (-2*D(3,1)*D(3,2)*X(3,4))-D(2,1)*D(3,2)*X(3,4)-D(1,1)*D(3,&
&2)*X(3,4)-D(2,2)*D(3,1)*X(3,4)-D(1,2)*D(3,1)*X(3,4)-D(2,1)*D(3,2)&
&*X(3,3)-D(2,2)*D(3,1)*X(3,3)-2*D(2,1)*D(2,2)*X(3,3)-D(1,1)*D(2,2)&
&*X(3,3)-D(1,2)*D(2,1)*X(3,3)-D(1,1)*D(3,2)*X(3,2)-D(1,2)*D(3,1)*X&
&(3,2)-D(1,1)*D(2,2)*X(3,2)-D(1,2)*D(2,1)*X(3,2)-2*D(1,1)*D(1,2)*X&
&(3,2)+2*D(3,1)*X(3,1)*D(3,2)+2*D(2,1)*X(3,1)*D(3,2)+2*D(1,1)*X(3,&
&1)*D(3,2)+2*D(2,2)*D(3,1)*X(3,1)+2*D(1,2)*D(3,1)*X(3,1)+2*D(2,1)*&
&D(2,2)*X(3,1)+2*D(1,1)*D(2,2)*X(3,1)+2*D(1,2)*D(2,1)*X(3,1)+2*D(1&
&,1)*D(1,2)*X(3,1)
tt82 = (-2*D(3,1)*D(3,3)*X(3,4))-D(2,1)*D(3,3)*X(3,4)-D(1,1)*D(3,&
&3)*X(3,4)-D(2,3)*D(3,1)*X(3,4)-D(1,3)*D(3,1)*X(3,4)-D(2,1)*D(3,3)&
&*X(3,3)-D(2,3)*D(3,1)*X(3,3)-2*D(2,1)*D(2,3)*X(3,3)-D(1,1)*D(2,3)&
&*X(3,3)-D(1,3)*D(2,1)*X(3,3)-D(1,1)*X(3,2)*D(3,3)+2*D(3,1)*X(3,1)&
&*D(3,3)+2*D(2,1)*X(3,1)*D(3,3)+2*D(1,1)*X(3,1)*D(3,3)-D(1,3)*D(3,&
&1)*X(3,2)-D(1,1)*D(2,3)*X(3,2)-D(1,3)*D(2,1)*X(3,2)-2*D(1,1)*D(1,&
&3)*X(3,2)+2*D(2,3)*D(3,1)*X(3,1)+2*D(1,3)*D(3,1)*X(3,1)+2*D(2,1)*&
&D(2,3)*X(3,1)+2*D(1,1)*D(2,3)*X(3,1)+2*D(1,3)*D(2,1)*X(3,1)+2*D(1&
&,1)*D(1,3)*X(3,1)
tt83 = (-2*D(3,2)*D(3,3)*X(3,4))-D(2,2)*D(3,3)*X(3,4)-D(1,2)*D(3,&
&3)*X(3,4)-D(2,3)*D(3,2)*X(3,4)-D(1,3)*D(3,2)*X(3,4)-D(2,2)*D(3,3)&
&*X(3,3)-D(2,3)*D(3,2)*X(3,3)-2*D(2,2)*D(2,3)*X(3,3)-D(1,2)*D(2,3)&
&*X(3,3)-D(1,3)*D(2,2)*X(3,3)-D(1,2)*X(3,2)*D(3,3)+2*X(3,1)*D(3,2)&
&*D(3,3)+2*D(2,2)*X(3,1)*D(3,3)+2*D(1,2)*X(3,1)*D(3,3)-D(1,3)*D(3,&
&2)*X(3,2)-D(1,2)*D(2,3)*X(3,2)-D(1,3)*D(2,2)*X(3,2)-2*D(1,2)*D(1,&
&3)*X(3,2)+2*D(2,3)*X(3,1)*D(3,2)+2*D(1,3)*X(3,1)*D(3,2)+2*D(2,2)*&
&D(2,3)*X(3,1)+2*D(1,2)*D(2,3)*X(3,1)+2*D(1,3)*D(2,2)*X(3,1)+2*D(1&
&,2)*D(1,3)*X(3,1)
tt84 = 1.0E+0*tt9*tt73+1.0E+0*tt7*tt71+1.0E+0*tt1*tt69
tt85 = volume(1,1)*(1.0E+0*lam(1,1)*tt11*tt84+miu(1,1)*(1.0E+0*tt&
&23*tt83+1.0E+0*tt22*tt82+2.0E+0*tt19*tt20*tt56+1.0E+0*tt18*tt81+2&
&.0E+0*tt15*tt16*tt50+2.0E+0*tt12*tt13*tt44))
tt86 = 1.0E+0*D(1,3)*tt10+1.0E+0*D(1,2)*tt8+1.0E+0*D(1,1)*tt6
tt87 = -X(1,1)*D(1,2)*D(2,1)
tt88 = -D(1,1)*X(1,1)*D(2,2)
tt89 = -X(1,1)*D(1,2)*D(3,1)
tt90 = -D(1,1)*X(1,1)*D(3,2)
tt91 = D(1,1)*X(1,4)*D(3,2)+tt90+D(1,2)*X(1,4)*D(3,1)+tt89+D(1,1)&
&*X(1,3)*D(2,2)+tt88+D(1,2)*X(1,3)*D(2,1)+tt87+2*D(1,1)*D(1,2)*X(1&
&,2)-2*D(1,1)*X(1,1)*D(1,2)
tt92 = -X(1,1)*D(1,3)*D(2,1)
tt93 = -D(1,1)*X(1,1)*D(2,3)
tt94 = -X(1,1)*D(1,3)*D(3,1)
tt95 = -D(1,1)*X(1,1)*D(3,3)
tt96 = D(1,1)*X(1,4)*D(3,3)+tt95+D(1,3)*X(1,4)*D(3,1)+tt94+D(1,1)&
&*X(1,3)*D(2,3)+tt93+D(1,3)*X(1,3)*D(2,1)+tt92+2*D(1,1)*X(1,2)*D(1&
&,3)-2*D(1,1)*X(1,1)*D(1,3)
tt97 = -X(1,1)*D(1,3)*D(2,2)
tt98 = -X(1,1)*D(1,2)*D(2,3)
tt99 = -X(1,1)*D(1,3)*D(3,2)
tt100 = -X(1,1)*D(1,2)*D(3,3)
tt101 = D(1,2)*X(1,4)*D(3,3)+tt100+D(1,3)*X(1,4)*D(3,2)+tt99+D(1,&
&2)*X(1,3)*D(2,3)+tt98+D(1,3)*X(1,3)*D(2,2)+tt97+2*D(1,2)*X(1,2)*D&
&(1,3)-2*X(1,1)*D(1,2)*D(1,3)
tt102 = -D(1,2)*D(2,1)
tt103 = -D(1,1)*D(2,2)
tt104 = -D(1,2)*D(3,1)
tt105 = -D(1,1)*D(3,2)
tt106 = 1.0E+0*(tt105+tt104+tt103+tt102-2*D(1,1)*D(1,2))*tt36
tt107 = -D(1,3)*D(2,1)
tt108 = -D(1,1)*D(2,3)
tt109 = -D(1,3)*D(3,1)
tt110 = -D(1,1)*D(3,3)
tt111 = 1.0E+0*(tt110+tt109+tt108+tt107-2*D(1,1)*D(1,3))*tt38
tt112 = -D(1,3)*D(2,2)
tt113 = -D(1,2)*D(2,3)
tt114 = -D(1,3)*D(3,2)
tt115 = -D(1,2)*D(3,3)
tt116 = 1.0E+0*(tt115+tt114+tt113+tt112-2*D(1,2)*D(1,3))*tt40
tt117 = 1.0E+0*D(1,1)*tt1*tt46
tt118 = 1.0E+0*D(1,2)*tt7*tt52
tt119 = 1.0E+0*D(1,3)*tt9*tt58
tt120 = 1.0E+0*lam(1,1)*(1.0E+0*D(1,3)*tt9+1.0E+0*D(1,2)*tt7+1.0E&
&+0*D(1,1)*tt1)*tt74
tt121 = volume(1,1)*(tt120+miu(1,1)*(tt119+tt118+tt117+tt116+tt11&
&1+tt106+2.0E+0*D(1,3)*tt9*tt21+1.0E+0*tt101*tt23+1.0E+0*tt96*tt22&
&+2.0E+0*D(1,2)*tt7*tt17+1.0E+0*tt91*tt18+2.0E+0*D(1,1)*tt1*tt14)+&
&1.0E+0*lam(1,1)*tt86*tt11)
tt122 = 1.0E+0*D(1,3)*tt72+1.0E+0*D(1,2)*tt70+1.0E+0*D(1,1)*tt64
tt123 = 2.0E+0*D(1,1)*tt1*tt13*tt42
tt124 = 2.0E+0*D(1,2)*tt7*tt16*tt48
tt125 = -D(1,2)*D(2,1)*X(2,1)
tt126 = -D(1,1)*X(2,1)*D(2,2)
tt127 = -D(1,2)*X(2,1)*D(3,1)
tt128 = -D(1,1)*X(2,1)*D(3,2)
tt129 = D(1,1)*X(2,4)*D(3,2)+tt128+D(1,2)*X(2,4)*D(3,1)+tt127+D(1&
&,1)*D(2,2)*X(2,3)+D(1,2)*D(2,1)*X(2,3)+2*D(1,1)*D(1,2)*X(2,2)+tt1&
&26+tt125-2*D(1,1)*D(1,2)*X(2,1)
tt130 = 2.0E+0*D(1,3)*tt9*tt20*tt54
tt131 = -D(1,3)*D(2,1)*X(2,1)
tt132 = -D(1,1)*X(2,1)*D(2,3)
tt133 = -D(1,3)*X(2,1)*D(3,1)
tt134 = -D(1,1)*X(2,1)*D(3,3)
tt135 = D(1,1)*X(2,4)*D(3,3)+tt134+D(1,3)*X(2,4)*D(3,1)+tt133+D(1&
&,1)*D(2,3)*X(2,3)+D(1,3)*D(2,1)*X(2,3)+tt132+2*D(1,1)*D(1,3)*X(2,&
&2)+tt131-2*D(1,1)*D(1,3)*X(2,1)
tt136 = -D(1,3)*X(2,1)*D(2,2)
tt137 = -D(1,2)*X(2,1)*D(2,3)
tt138 = -D(1,3)*X(2,1)*D(3,2)
tt139 = -D(1,2)*X(2,1)*D(3,3)
tt140 = D(1,2)*X(2,4)*D(3,3)+tt139+D(1,3)*X(2,4)*D(3,2)+tt138+D(1&
&,2)*D(2,3)*X(2,3)+D(1,3)*D(2,2)*X(2,3)+tt137+2*D(1,2)*D(1,3)*X(2,&
&2)+tt136-2*D(1,2)*D(1,3)*X(2,1)
tt141 = volume(1,1)*(miu(1,1)*(1.0E+0*tt140*tt23+1.0E+0*tt135*tt2&
&2+tt130+1.0E+0*tt129*tt18+tt124+tt123)+1.0E+0*lam(1,1)*tt11*tt122&
&)
tt142 = 2.0E+0*D(1,1)*tt1*tt13*tt44
tt143 = 2.0E+0*D(1,2)*tt7*tt16*tt50
tt144 = -D(1,2)*D(2,1)*X(3,1)
tt145 = -D(1,1)*D(2,2)*X(3,1)
tt146 = -D(1,2)*D(3,1)*X(3,1)
tt147 = -D(1,1)*X(3,1)*D(3,2)
tt148 = D(1,1)*D(3,2)*X(3,4)+D(1,2)*D(3,1)*X(3,4)+D(1,1)*D(2,2)*X&
&(3,3)+D(1,2)*D(2,1)*X(3,3)+2*D(1,1)*D(1,2)*X(3,2)+tt147+tt146+tt1&
&45+tt144-2*D(1,1)*D(1,2)*X(3,1)
tt149 = 2.0E+0*D(1,3)*tt9*tt20*tt56
tt150 = -D(1,3)*D(2,1)*X(3,1)
tt151 = -D(1,1)*D(2,3)*X(3,1)
tt152 = -D(1,3)*D(3,1)*X(3,1)
tt153 = -D(1,1)*X(3,1)*D(3,3)
tt154 = D(1,1)*D(3,3)*X(3,4)+D(1,3)*D(3,1)*X(3,4)+D(1,1)*D(2,3)*X&
&(3,3)+D(1,3)*D(2,1)*X(3,3)+tt153+2*D(1,1)*D(1,3)*X(3,2)+tt152+tt1&
&51+tt150-2*D(1,1)*D(1,3)*X(3,1)
tt155 = -D(1,3)*D(2,2)*X(3,1)
tt156 = -D(1,2)*D(2,3)*X(3,1)
tt157 = -D(1,3)*X(3,1)*D(3,2)
tt158 = -D(1,2)*X(3,1)*D(3,3)
tt159 = D(1,2)*D(3,3)*X(3,4)+D(1,3)*D(3,2)*X(3,4)+D(1,2)*D(2,3)*X&
&(3,3)+D(1,3)*D(2,2)*X(3,3)+tt158+2*D(1,2)*D(1,3)*X(3,2)+tt157+tt1&
&56+tt155-2*D(1,2)*D(1,3)*X(3,1)
tt160 = 1.0E+0*D(1,3)*tt73+1.0E+0*D(1,2)*tt71+1.0E+0*D(1,1)*tt69
tt161 = volume(1,1)*(1.0E+0*lam(1,1)*tt11*tt160+miu(1,1)*(1.0E+0*&
&tt23*tt159+1.0E+0*tt22*tt154+tt149+1.0E+0*tt18*tt148+tt143+tt142)&
&)
tt162 = 1.0E+0*D(2,3)*tt10+1.0E+0*D(2,2)*tt8+1.0E+0*D(2,1)*tt6
tt163 = -X(1,1)*D(2,2)*D(3,1)
tt164 = -X(1,1)*D(2,1)*D(3,2)
tt165 = X(1,4)*D(2,1)*D(3,2)+tt164+X(1,4)*D(2,2)*D(3,1)+tt163+2*X&
&(1,3)*D(2,1)*D(2,2)-2*X(1,1)*D(2,1)*D(2,2)+D(1,1)*X(1,2)*D(2,2)+t&
&t88+D(1,2)*X(1,2)*D(2,1)+tt87
tt166 = -X(1,1)*D(2,3)*D(3,1)
tt167 = -X(1,1)*D(2,1)*D(3,3)
tt168 = X(1,4)*D(2,1)*D(3,3)+tt167+X(1,4)*D(2,3)*D(3,1)+tt166+2*X&
&(1,3)*D(2,1)*D(2,3)-2*X(1,1)*D(2,1)*D(2,3)+D(1,1)*X(1,2)*D(2,3)+t&
&t93+X(1,2)*D(1,3)*D(2,1)+tt92
tt169 = -X(1,1)*D(2,3)*D(3,2)
tt170 = -X(1,1)*D(2,2)*D(3,3)
tt171 = X(1,4)*D(2,2)*D(3,3)+tt170+X(1,4)*D(2,3)*D(3,2)+tt169+2*X&
&(1,3)*D(2,2)*D(2,3)-2*X(1,1)*D(2,2)*D(2,3)+D(1,2)*X(1,2)*D(2,3)+t&
&t98+X(1,2)*D(1,3)*D(2,2)+tt97
tt172 = -D(2,2)*D(3,1)
tt173 = -D(2,1)*D(3,2)
tt174 = 1.0E+0*(tt173+tt172-2*D(2,1)*D(2,2)+tt103+tt102)*tt36
tt175 = -D(2,3)*D(3,1)
tt176 = -D(2,1)*D(3,3)
tt177 = 1.0E+0*(tt176+tt175-2*D(2,1)*D(2,3)+tt108+tt107)*tt38
tt178 = -D(2,3)*D(3,2)
tt179 = -D(2,2)*D(3,3)
tt180 = 1.0E+0*(tt179+tt178-2*D(2,2)*D(2,3)+tt113+tt112)*tt40
tt181 = 1.0E+0*D(2,1)*tt1*tt46
tt182 = 1.0E+0*D(2,2)*tt7*tt52
tt183 = 1.0E+0*D(2,3)*tt9*tt58
tt184 = 1.0E+0*lam(1,1)*(1.0E+0*D(2,3)*tt9+1.0E+0*D(2,2)*tt7+1.0E&
&+0*D(2,1)*tt1)*tt74
tt185 = volume(1,1)*(tt184+miu(1,1)*(tt183+tt182+tt181+tt180+tt17&
&7+tt174+2.0E+0*D(2,3)*tt9*tt21+1.0E+0*tt171*tt23+1.0E+0*tt168*tt2&
&2+2.0E+0*D(2,2)*tt7*tt17+1.0E+0*tt165*tt18+2.0E+0*D(2,1)*tt1*tt14&
&)+1.0E+0*lam(1,1)*tt162*tt11)
tt186 = 1.0E+0*D(2,3)*tt72+1.0E+0*D(2,2)*tt70+1.0E+0*D(2,1)*tt64
tt187 = 2.0E+0*D(2,1)*tt1*tt13*tt42
tt188 = 2.0E+0*D(2,2)*tt7*tt16*tt48
tt189 = -X(2,1)*D(2,2)*D(3,1)
tt190 = -D(2,1)*X(2,1)*D(3,2)
tt191 = D(2,1)*X(2,4)*D(3,2)+tt190+D(2,2)*X(2,4)*D(3,1)+tt189+2*D&
&(2,1)*D(2,2)*X(2,3)+D(1,1)*D(2,2)*X(2,2)+D(1,2)*D(2,1)*X(2,2)-2*D&
&(2,1)*X(2,1)*D(2,2)+tt126+tt125
tt192 = 2.0E+0*D(2,3)*tt9*tt20*tt54
tt193 = -X(2,1)*D(2,3)*D(3,1)
tt194 = -D(2,1)*X(2,1)*D(3,3)
tt195 = D(2,1)*X(2,4)*D(3,3)+tt194+D(2,3)*X(2,4)*D(3,1)+tt193+2*D&
&(2,1)*D(2,3)*X(2,3)+D(1,1)*X(2,2)*D(2,3)-2*D(2,1)*X(2,1)*D(2,3)+t&
&t132+D(1,3)*D(2,1)*X(2,2)+tt131
tt196 = -X(2,1)*D(2,3)*D(3,2)
tt197 = -X(2,1)*D(2,2)*D(3,3)
tt198 = D(2,2)*X(2,4)*D(3,3)+tt197+D(2,3)*X(2,4)*D(3,2)+tt196+2*D&
&(2,2)*D(2,3)*X(2,3)+D(1,2)*X(2,2)*D(2,3)-2*X(2,1)*D(2,2)*D(2,3)+t&
&t137+D(1,3)*D(2,2)*X(2,2)+tt136
tt199 = volume(1,1)*(miu(1,1)*(1.0E+0*tt198*tt23+1.0E+0*tt195*tt2&
&2+tt192+1.0E+0*tt191*tt18+tt188+tt187)+1.0E+0*lam(1,1)*tt11*tt186&
&)
tt200 = 2.0E+0*D(2,1)*tt1*tt13*tt44
tt201 = 2.0E+0*D(2,2)*tt7*tt16*tt50
tt202 = -D(2,2)*D(3,1)*X(3,1)
tt203 = -D(2,1)*X(3,1)*D(3,2)
tt204 = D(2,1)*D(3,2)*X(3,4)+D(2,2)*D(3,1)*X(3,4)+2*D(2,1)*D(2,2)&
&*X(3,3)+D(1,1)*D(2,2)*X(3,2)+D(1,2)*D(2,1)*X(3,2)+tt203+tt202-2*D&
&(2,1)*D(2,2)*X(3,1)+tt145+tt144
tt205 = 2.0E+0*D(2,3)*tt9*tt20*tt56
tt206 = -D(2,3)*D(3,1)*X(3,1)
tt207 = -D(2,1)*X(3,1)*D(3,3)
tt208 = D(2,1)*D(3,3)*X(3,4)+D(2,3)*D(3,1)*X(3,4)+2*D(2,1)*D(2,3)&
&*X(3,3)+tt207+D(1,1)*D(2,3)*X(3,2)+D(1,3)*D(2,1)*X(3,2)+tt206-2*D&
&(2,1)*D(2,3)*X(3,1)+tt151+tt150
tt209 = -D(2,3)*X(3,1)*D(3,2)
tt210 = -D(2,2)*X(3,1)*D(3,3)
tt211 = D(2,2)*D(3,3)*X(3,4)+D(2,3)*D(3,2)*X(3,4)+2*D(2,2)*D(2,3)&
&*X(3,3)+tt210+D(1,2)*D(2,3)*X(3,2)+D(1,3)*D(2,2)*X(3,2)+tt209-2*D&
&(2,2)*D(2,3)*X(3,1)+tt156+tt155
tt212 = 1.0E+0*D(2,3)*tt73+1.0E+0*D(2,2)*tt71+1.0E+0*D(2,1)*tt69
tt213 = volume(1,1)*(1.0E+0*lam(1,1)*tt11*tt212+miu(1,1)*(1.0E+0*&
&tt23*tt211+1.0E+0*tt22*tt208+tt205+1.0E+0*tt18*tt204+tt201+tt200)&
&)
tt214 = 1.0E+0*D(3,3)*tt10+1.0E+0*D(3,2)*tt8+1.0E+0*D(3,1)*tt6
tt215 = 2*X(1,4)*D(3,1)*D(3,2)-2*X(1,1)*D(3,1)*D(3,2)+X(1,3)*D(2,&
&1)*D(3,2)+tt164+D(1,1)*X(1,2)*D(3,2)+tt90+X(1,3)*D(2,2)*D(3,1)+tt&
&163+D(1,2)*X(1,2)*D(3,1)+tt89
tt216 = 2*X(1,4)*D(3,1)*D(3,3)-2*X(1,1)*D(3,1)*D(3,3)+X(1,3)*D(2,&
&1)*D(3,3)+tt167+D(1,1)*X(1,2)*D(3,3)+tt95+X(1,3)*D(2,3)*D(3,1)+tt&
&166+X(1,2)*D(1,3)*D(3,1)+tt94
tt217 = 2*X(1,4)*D(3,2)*D(3,3)-2*X(1,1)*D(3,2)*D(3,3)+X(1,3)*D(2,&
&2)*D(3,3)+tt170+D(1,2)*X(1,2)*D(3,3)+tt100+X(1,3)*D(2,3)*D(3,2)+t&
&t169+X(1,2)*D(1,3)*D(3,2)+tt99
tt218 = 1.0E+0*((-2*D(3,1)*D(3,2))+tt173+tt105+tt172+tt104)*tt36
tt219 = 1.0E+0*((-2*D(3,1)*D(3,3))+tt176+tt110+tt175+tt109)*tt38
tt220 = 1.0E+0*((-2*D(3,2)*D(3,3))+tt179+tt115+tt178+tt114)*tt40
tt221 = 1.0E+0*tt1*D(3,1)*tt46
tt222 = 1.0E+0*tt7*D(3,2)*tt52
tt223 = 1.0E+0*tt9*D(3,3)*tt58
tt224 = 1.0E+0*lam(1,1)*(1.0E+0*tt9*D(3,3)+1.0E+0*tt7*D(3,2)+1.0E&
&+0*tt1*D(3,1))*tt74
tt225 = volume(1,1)*(tt224+miu(1,1)*(tt223+tt222+tt221+tt220+tt21&
&9+tt218+2.0E+0*tt9*D(3,3)*tt21+1.0E+0*tt23*tt217+1.0E+0*tt22*tt21&
&6+2.0E+0*tt7*D(3,2)*tt17+1.0E+0*tt18*tt215+2.0E+0*tt1*D(3,1)*tt14&
&)+1.0E+0*lam(1,1)*tt11*tt214)
tt226 = 1.0E+0*D(3,3)*tt72+1.0E+0*D(3,2)*tt70+1.0E+0*D(3,1)*tt64
tt227 = 2.0E+0*tt1*D(3,1)*tt13*tt42
tt228 = 2.0E+0*tt7*D(3,2)*tt16*tt48
tt229 = 2*X(2,4)*D(3,1)*D(3,2)-2*X(2,1)*D(3,1)*D(3,2)+D(2,1)*X(2,&
&3)*D(3,2)+D(1,1)*X(2,2)*D(3,2)+tt190+tt128+D(2,2)*X(2,3)*D(3,1)+D&
&(1,2)*X(2,2)*D(3,1)+tt189+tt127
tt230 = 2.0E+0*tt9*D(3,3)*tt20*tt54
tt231 = 2*X(2,4)*D(3,1)*D(3,3)-2*X(2,1)*D(3,1)*D(3,3)+D(2,1)*X(2,&
&3)*D(3,3)+D(1,1)*X(2,2)*D(3,3)+tt194+tt134+D(2,3)*X(2,3)*D(3,1)+t&
&t193+D(1,3)*X(2,2)*D(3,1)+tt133
tt232 = 2*X(2,4)*D(3,2)*D(3,3)-2*X(2,1)*D(3,2)*D(3,3)+D(2,2)*X(2,&
&3)*D(3,3)+D(1,2)*X(2,2)*D(3,3)+tt197+tt139+D(2,3)*X(2,3)*D(3,2)+t&
&t196+D(1,3)*X(2,2)*D(3,2)+tt138
tt233 = volume(1,1)*(miu(1,1)*(1.0E+0*tt23*tt232+1.0E+0*tt22*tt23&
&1+tt230+1.0E+0*tt18*tt229+tt228+tt227)+1.0E+0*lam(1,1)*tt11*tt226&
&)
tt234 = 2.0E+0*tt1*D(3,1)*tt13*tt44
tt235 = 2.0E+0*tt7*D(3,2)*tt16*tt50
tt236 = 2*D(3,1)*D(3,2)*X(3,4)+D(2,1)*D(3,2)*X(3,3)+D(2,2)*D(3,1)&
&*X(3,3)+D(1,1)*D(3,2)*X(3,2)+D(1,2)*D(3,1)*X(3,2)-2*D(3,1)*X(3,1)&
&*D(3,2)+tt203+tt147+tt202+tt146
tt237 = 2.0E+0*tt9*D(3,3)*tt20*tt56
tt238 = 2*D(3,1)*D(3,3)*X(3,4)+D(2,1)*D(3,3)*X(3,3)+D(2,3)*D(3,1)&
&*X(3,3)+D(1,1)*X(3,2)*D(3,3)-2*D(3,1)*X(3,1)*D(3,3)+tt207+tt153+D&
&(1,3)*D(3,1)*X(3,2)+tt206+tt152
tt239 = 2*D(3,2)*D(3,3)*X(3,4)+D(2,2)*D(3,3)*X(3,3)+D(2,3)*D(3,2)&
&*X(3,3)+D(1,2)*X(3,2)*D(3,3)-2*X(3,1)*D(3,2)*D(3,3)+tt210+tt158+D&
&(1,3)*D(3,2)*X(3,2)+tt209+tt157
tt240 = 1.0E+0*D(3,3)*tt73+1.0E+0*D(3,2)*tt71+1.0E+0*D(3,1)*tt69
tt241 = volume(1,1)*(1.0E+0*lam(1,1)*tt11*tt240+miu(1,1)*(1.0E+0*&
&tt23*tt239+1.0E+0*tt22*tt238+tt237+1.0E+0*tt18*tt236+tt235+tt234)&
&)
tt242 = volume(1,1)*(1.0E+0*lam(1,1)*tt76*tt84+miu(1,1)*(1.0E+0*t&
&t79*tt83+1.0E+0*tt78*tt82+2.0E+0*tt19*tt54*tt56+1.0E+0*tt77*tt81+&
&2.0E+0*tt15*tt48*tt50+2.0E+0*tt12*tt42*tt44))
tt243 = volume(1,1)*(miu(1,1)*(1.0E+0*tt101*tt79+1.0E+0*tt96*tt78&
&+tt130+1.0E+0*tt91*tt77+tt124+tt123)+1.0E+0*lam(1,1)*tt86*tt76)
tt244 = volume(1,1)*(tt120+miu(1,1)*(tt119+tt118+tt117+tt116+tt11&
&1+tt106+2.0E+0*D(1,3)*tt9*tt55+1.0E+0*tt140*tt79+1.0E+0*tt135*tt7&
&8+2.0E+0*D(1,2)*tt7*tt49+1.0E+0*tt129*tt77+2.0E+0*D(1,1)*tt1*tt43&
&)+1.0E+0*lam(1,1)*tt122*tt76)
tt245 = 2.0E+0*D(1,1)*tt1*tt42*tt44
tt246 = 2.0E+0*D(1,2)*tt7*tt48*tt50
tt247 = 2.0E+0*D(1,3)*tt9*tt54*tt56
tt248 = volume(1,1)*(1.0E+0*lam(1,1)*tt76*tt160+miu(1,1)*(1.0E+0*&
&tt79*tt159+1.0E+0*tt78*tt154+tt247+1.0E+0*tt77*tt148+tt246+tt245)&
&)
tt249 = volume(1,1)*(miu(1,1)*(1.0E+0*tt171*tt79+1.0E+0*tt168*tt7&
&8+tt192+1.0E+0*tt165*tt77+tt188+tt187)+1.0E+0*lam(1,1)*tt162*tt76&
&)
tt250 = volume(1,1)*(tt184+miu(1,1)*(tt183+tt182+tt181+tt180+tt17&
&7+tt174+2.0E+0*D(2,3)*tt9*tt55+1.0E+0*tt198*tt79+1.0E+0*tt195*tt7&
&8+2.0E+0*D(2,2)*tt7*tt49+1.0E+0*tt191*tt77+2.0E+0*D(2,1)*tt1*tt43&
&)+1.0E+0*lam(1,1)*tt186*tt76)
tt251 = 2.0E+0*D(2,1)*tt1*tt42*tt44
tt252 = 2.0E+0*D(2,2)*tt7*tt48*tt50
tt253 = 2.0E+0*D(2,3)*tt9*tt54*tt56
tt254 = volume(1,1)*(1.0E+0*lam(1,1)*tt76*tt212+miu(1,1)*(1.0E+0*&
&tt79*tt211+1.0E+0*tt78*tt208+tt253+1.0E+0*tt77*tt204+tt252+tt251)&
&)
tt255 = volume(1,1)*(miu(1,1)*(1.0E+0*tt217*tt79+1.0E+0*tt216*tt7&
&8+tt230+1.0E+0*tt215*tt77+tt228+tt227)+1.0E+0*lam(1,1)*tt214*tt76&
&)
tt256 = volume(1,1)*(tt224+miu(1,1)*(tt223+tt222+tt221+tt220+tt21&
&9+tt218+2.0E+0*tt9*D(3,3)*tt55+1.0E+0*tt79*tt232+1.0E+0*tt78*tt23&
&1+2.0E+0*tt7*D(3,2)*tt49+1.0E+0*tt77*tt229+2.0E+0*tt1*D(3,1)*tt43&
&)+1.0E+0*lam(1,1)*tt76*tt226)
tt257 = 2.0E+0*tt1*D(3,1)*tt42*tt44
tt258 = 2.0E+0*tt7*D(3,2)*tt48*tt50
tt259 = 2.0E+0*tt9*D(3,3)*tt54*tt56
tt260 = volume(1,1)*(1.0E+0*lam(1,1)*tt76*tt240+miu(1,1)*(1.0E+0*&
&tt79*tt239+1.0E+0*tt78*tt238+tt259+1.0E+0*tt77*tt236+tt258+tt257)&
&)
tt261 = volume(1,1)*(1.0E+0*lam(1,1)*tt86*tt84+miu(1,1)*(1.0E+0*t&
&t101*tt83+1.0E+0*tt96*tt82+tt149+1.0E+0*tt91*tt81+tt143+tt142))
tt262 = volume(1,1)*(1.0E+0*lam(1,1)*tt122*tt84+miu(1,1)*(1.0E+0*&
&tt140*tt83+1.0E+0*tt135*tt82+tt247+1.0E+0*tt129*tt81+tt246+tt245)&
&)
tt263 = volume(1,1)*(tt120+miu(1,1)*(tt119+2.0E+0*D(1,3)*tt9*tt57&
&+tt118+2.0E+0*D(1,2)*tt7*tt51+tt117+2.0E+0*D(1,1)*tt1*tt45+tt116+&
&tt111+tt106+1.0E+0*tt159*tt83+1.0E+0*tt154*tt82+1.0E+0*tt148*tt81&
&)+1.0E+0*lam(1,1)*tt160*tt84)
tt264 = volume(1,1)*(1.0E+0*lam(1,1)*tt162*tt84+miu(1,1)*(1.0E+0*&
&tt171*tt83+1.0E+0*tt168*tt82+tt205+1.0E+0*tt165*tt81+tt201+tt200)&
&)
tt265 = volume(1,1)*(1.0E+0*lam(1,1)*tt186*tt84+miu(1,1)*(1.0E+0*&
&tt198*tt83+1.0E+0*tt195*tt82+tt253+1.0E+0*tt191*tt81+tt252+tt251)&
&)
tt266 = volume(1,1)*(tt184+miu(1,1)*(tt183+2.0E+0*D(2,3)*tt9*tt57&
&+tt182+2.0E+0*D(2,2)*tt7*tt51+tt181+2.0E+0*D(2,1)*tt1*tt45+tt180+&
&tt177+tt174+1.0E+0*tt211*tt83+1.0E+0*tt208*tt82+1.0E+0*tt204*tt81&
&)+1.0E+0*lam(1,1)*tt212*tt84)
tt267 = volume(1,1)*(1.0E+0*lam(1,1)*tt214*tt84+miu(1,1)*(1.0E+0*&
&tt217*tt83+1.0E+0*tt216*tt82+tt237+1.0E+0*tt215*tt81+tt235+tt234)&
&)
tt268 = volume(1,1)*(1.0E+0*lam(1,1)*tt226*tt84+miu(1,1)*(1.0E+0*&
&tt232*tt83+1.0E+0*tt231*tt82+tt259+1.0E+0*tt229*tt81+tt258+tt257)&
&)
tt269 = volume(1,1)*(tt224+miu(1,1)*(tt223+2.0E+0*tt9*D(3,3)*tt57&
&+tt222+2.0E+0*tt7*D(3,2)*tt51+tt221+2.0E+0*tt1*D(3,1)*tt45+tt220+&
&tt219+tt218+1.0E+0*tt83*tt239+1.0E+0*tt82*tt238+1.0E+0*tt81*tt236&
&)+1.0E+0*lam(1,1)*tt84*tt240)
tt270 = D(1,1)**2
tt271 = D(1,2)**2
tt272 = D(1,3)**2
tt273 = 2.0E+0*D(1,1)*D(1,2)*tt36
tt274 = 2.0E+0*D(1,1)*D(1,3)*tt38
tt275 = 2.0E+0*D(1,2)*D(1,3)*tt40
tt276 = 1.0E+0*tt270*tt46
tt277 = 1.0E+0*tt271*tt52
tt278 = 1.0E+0*tt272*tt58
tt279 = 1.0E+0*lam(1,1)*(1.0E+0*tt272+1.0E+0*tt271+1.0E+0*tt270)*&
&tt74
tt280 = volume(1,1)*(1.0E+0*lam(1,1)*tt86*tt122+miu(1,1)*(1.0E+0*&
&tt101*tt140+1.0E+0*tt96*tt135+2.0E+0*tt272*tt20*tt54+1.0E+0*tt91*&
&tt129+2.0E+0*tt271*tt16*tt48+2.0E+0*tt270*tt13*tt42))
tt281 = volume(1,1)*(1.0E+0*lam(1,1)*tt86*tt160+miu(1,1)*(1.0E+0*&
&tt101*tt159+1.0E+0*tt96*tt154+2.0E+0*tt272*tt20*tt56+1.0E+0*tt91*&
&tt148+2.0E+0*tt271*tt16*tt50+2.0E+0*tt270*tt13*tt44))
tt282 = 1.0E+0*(D(1,1)*D(2,2)+D(1,2)*D(2,1))*tt36
tt283 = 1.0E+0*(D(1,1)*D(2,3)+D(1,3)*D(2,1))*tt38
tt284 = 1.0E+0*(D(1,2)*D(2,3)+D(1,3)*D(2,2))*tt40
tt285 = 1.0E+0*D(1,1)*D(2,1)*tt46
tt286 = 1.0E+0*D(1,2)*D(2,2)*tt52
tt287 = 1.0E+0*D(1,3)*D(2,3)*tt58
tt288 = 1.0E+0*lam(1,1)*(1.0E+0*D(1,3)*D(2,3)+1.0E+0*D(1,2)*D(2,2&
&)+1.0E+0*D(1,1)*D(2,1))*tt74
tt289 = volume(1,1)*(tt288+miu(1,1)*(tt287+tt286+tt285+tt284+tt28&
&3+tt282+2.0E+0*D(1,3)*D(2,3)*tt21+1.0E+0*tt101*tt171+1.0E+0*tt96*&
&tt168+2.0E+0*D(1,2)*D(2,2)*tt17+1.0E+0*tt91*tt165+2.0E+0*D(1,1)*D&
&(2,1)*tt14)+1.0E+0*lam(1,1)*tt86*tt162)
tt290 = 2.0E+0*D(1,1)*D(2,1)*tt13*tt42
tt291 = 2.0E+0*D(1,2)*D(2,2)*tt16*tt48
tt292 = 2.0E+0*D(1,3)*D(2,3)*tt20*tt54
tt293 = volume(1,1)*(1.0E+0*lam(1,1)*tt86*tt186+miu(1,1)*(1.0E+0*&
&tt101*tt198+1.0E+0*tt96*tt195+tt292+1.0E+0*tt91*tt191+tt291+tt290&
&))
tt294 = 2.0E+0*D(1,1)*D(2,1)*tt13*tt44
tt295 = 2.0E+0*D(1,2)*D(2,2)*tt16*tt50
tt296 = 2.0E+0*D(1,3)*D(2,3)*tt20*tt56
tt297 = volume(1,1)*(1.0E+0*lam(1,1)*tt86*tt212+miu(1,1)*(1.0E+0*&
&tt101*tt211+1.0E+0*tt96*tt208+tt296+1.0E+0*tt91*tt204+tt295+tt294&
&))
tt298 = 1.0E+0*(D(1,1)*D(3,2)+D(1,2)*D(3,1))*tt36
tt299 = 1.0E+0*(D(1,1)*D(3,3)+D(1,3)*D(3,1))*tt38
tt300 = 1.0E+0*(D(1,2)*D(3,3)+D(1,3)*D(3,2))*tt40
tt301 = 1.0E+0*D(1,1)*D(3,1)*tt46
tt302 = 1.0E+0*D(1,2)*D(3,2)*tt52
tt303 = 1.0E+0*D(1,3)*D(3,3)*tt58
tt304 = 1.0E+0*lam(1,1)*(1.0E+0*D(1,3)*D(3,3)+1.0E+0*D(1,2)*D(3,2&
&)+1.0E+0*D(1,1)*D(3,1))*tt74
tt305 = volume(1,1)*(tt304+miu(1,1)*(tt303+tt302+tt301+tt300+tt29&
&9+tt298+2.0E+0*D(1,3)*D(3,3)*tt21+1.0E+0*tt101*tt217+1.0E+0*tt96*&
&tt216+2.0E+0*D(1,2)*D(3,2)*tt17+1.0E+0*tt91*tt215+2.0E+0*D(1,1)*D&
&(3,1)*tt14)+1.0E+0*lam(1,1)*tt86*tt214)
tt306 = 2.0E+0*D(1,1)*D(3,1)*tt13*tt42
tt307 = 2.0E+0*D(1,2)*D(3,2)*tt16*tt48
tt308 = 2.0E+0*D(1,3)*D(3,3)*tt20*tt54
tt309 = volume(1,1)*(miu(1,1)*(1.0E+0*tt101*tt232+1.0E+0*tt96*tt2&
&31+tt308+1.0E+0*tt91*tt229+tt307+tt306)+1.0E+0*lam(1,1)*tt86*tt22&
&6)
tt310 = 2.0E+0*D(1,1)*D(3,1)*tt13*tt44
tt311 = 2.0E+0*D(1,2)*D(3,2)*tt16*tt50
tt312 = 2.0E+0*D(1,3)*D(3,3)*tt20*tt56
tt313 = volume(1,1)*(1.0E+0*lam(1,1)*tt86*tt240+miu(1,1)*(1.0E+0*&
&tt101*tt239+1.0E+0*tt96*tt238+tt312+1.0E+0*tt91*tt236+tt311+tt310&
&))
tt314 = volume(1,1)*(1.0E+0*lam(1,1)*tt122*tt160+miu(1,1)*(1.0E+0&
&*tt140*tt159+1.0E+0*tt135*tt154+2.0E+0*tt272*tt54*tt56+1.0E+0*tt1&
&29*tt148+2.0E+0*tt271*tt48*tt50+2.0E+0*tt270*tt42*tt44))
tt315 = volume(1,1)*(1.0E+0*lam(1,1)*tt162*tt122+miu(1,1)*(1.0E+0&
&*tt171*tt140+1.0E+0*tt168*tt135+tt292+1.0E+0*tt165*tt129+tt291+tt&
&290))
tt316 = volume(1,1)*(tt288+miu(1,1)*(tt287+tt286+tt285+tt284+tt28&
&3+tt282+2.0E+0*D(1,3)*D(2,3)*tt55+1.0E+0*tt140*tt198+1.0E+0*tt135&
&*tt195+2.0E+0*D(1,2)*D(2,2)*tt49+1.0E+0*tt129*tt191+2.0E+0*D(1,1)&
&*D(2,1)*tt43)+1.0E+0*lam(1,1)*tt122*tt186)
tt317 = 2.0E+0*D(1,1)*D(2,1)*tt42*tt44
tt318 = 2.0E+0*D(1,2)*D(2,2)*tt48*tt50
tt319 = 2.0E+0*D(1,3)*D(2,3)*tt54*tt56
tt320 = volume(1,1)*(1.0E+0*lam(1,1)*tt122*tt212+miu(1,1)*(1.0E+0&
&*tt140*tt211+1.0E+0*tt135*tt208+tt319+1.0E+0*tt129*tt204+tt318+tt&
&317))
tt321 = volume(1,1)*(miu(1,1)*(1.0E+0*tt140*tt217+1.0E+0*tt135*tt&
&216+tt308+1.0E+0*tt129*tt215+tt307+tt306)+1.0E+0*lam(1,1)*tt214*t&
&t122)
tt322 = volume(1,1)*(tt304+miu(1,1)*(tt303+tt302+tt301+tt300+tt29&
&9+tt298+2.0E+0*D(1,3)*D(3,3)*tt55+1.0E+0*tt140*tt232+1.0E+0*tt135&
&*tt231+2.0E+0*D(1,2)*D(3,2)*tt49+1.0E+0*tt129*tt229+2.0E+0*D(1,1)&
&*D(3,1)*tt43)+1.0E+0*lam(1,1)*tt122*tt226)
tt323 = 2.0E+0*D(1,1)*D(3,1)*tt42*tt44
tt324 = 2.0E+0*D(1,2)*D(3,2)*tt48*tt50
tt325 = 2.0E+0*D(1,3)*D(3,3)*tt54*tt56
tt326 = volume(1,1)*(1.0E+0*lam(1,1)*tt122*tt240+miu(1,1)*(1.0E+0&
&*tt140*tt239+1.0E+0*tt135*tt238+tt325+1.0E+0*tt129*tt236+tt324+tt&
&323))
tt327 = volume(1,1)*(1.0E+0*lam(1,1)*tt162*tt160+miu(1,1)*(1.0E+0&
&*tt171*tt159+1.0E+0*tt168*tt154+tt296+1.0E+0*tt165*tt148+tt295+tt&
&294))
tt328 = volume(1,1)*(1.0E+0*lam(1,1)*tt186*tt160+miu(1,1)*(1.0E+0&
&*tt198*tt159+1.0E+0*tt195*tt154+tt319+1.0E+0*tt191*tt148+tt318+tt&
&317))
tt329 = volume(1,1)*(tt288+miu(1,1)*(tt287+2.0E+0*D(1,3)*D(2,3)*t&
&t57+tt286+2.0E+0*D(1,2)*D(2,2)*tt51+tt285+2.0E+0*D(1,1)*D(2,1)*tt&
&45+tt284+tt283+tt282+1.0E+0*tt159*tt211+1.0E+0*tt154*tt208+1.0E+0&
&*tt148*tt204)+1.0E+0*lam(1,1)*tt160*tt212)
tt330 = volume(1,1)*(1.0E+0*lam(1,1)*tt214*tt160+miu(1,1)*(1.0E+0&
&*tt217*tt159+1.0E+0*tt216*tt154+tt312+1.0E+0*tt215*tt148+tt311+tt&
&310))
tt331 = volume(1,1)*(1.0E+0*lam(1,1)*tt226*tt160+miu(1,1)*(1.0E+0&
&*tt232*tt159+1.0E+0*tt231*tt154+tt325+1.0E+0*tt229*tt148+tt324+tt&
&323))
tt332 = volume(1,1)*(tt304+miu(1,1)*(tt303+2.0E+0*D(1,3)*D(3,3)*t&
&t57+tt302+2.0E+0*D(1,2)*D(3,2)*tt51+tt301+2.0E+0*D(1,1)*D(3,1)*tt&
&45+tt300+tt299+tt298+1.0E+0*tt159*tt239+1.0E+0*tt154*tt238+1.0E+0&
&*tt148*tt236)+1.0E+0*lam(1,1)*tt160*tt240)
tt333 = D(2,1)**2
tt334 = D(2,2)**2
tt335 = D(2,3)**2
tt336 = 2.0E+0*D(2,1)*D(2,2)*tt36
tt337 = 2.0E+0*D(2,1)*D(2,3)*tt38
tt338 = 2.0E+0*D(2,2)*D(2,3)*tt40
tt339 = 1.0E+0*tt333*tt46
tt340 = 1.0E+0*tt334*tt52
tt341 = 1.0E+0*tt335*tt58
tt342 = 1.0E+0*lam(1,1)*(1.0E+0*tt335+1.0E+0*tt334+1.0E+0*tt333)*&
&tt74
tt343 = volume(1,1)*(1.0E+0*lam(1,1)*tt162*tt186+miu(1,1)*(1.0E+0&
&*tt171*tt198+1.0E+0*tt168*tt195+2.0E+0*tt335*tt20*tt54+1.0E+0*tt1&
&65*tt191+2.0E+0*tt334*tt16*tt48+2.0E+0*tt333*tt13*tt42))
tt344 = volume(1,1)*(1.0E+0*lam(1,1)*tt162*tt212+miu(1,1)*(1.0E+0&
&*tt171*tt211+1.0E+0*tt168*tt208+2.0E+0*tt335*tt20*tt56+1.0E+0*tt1&
&65*tt204+2.0E+0*tt334*tt16*tt50+2.0E+0*tt333*tt13*tt44))
tt345 = 1.0E+0*(D(2,1)*D(3,2)+D(2,2)*D(3,1))*tt36
tt346 = 1.0E+0*(D(2,1)*D(3,3)+D(2,3)*D(3,1))*tt38
tt347 = 1.0E+0*(D(2,2)*D(3,3)+D(2,3)*D(3,2))*tt40
tt348 = 1.0E+0*D(2,1)*D(3,1)*tt46
tt349 = 1.0E+0*D(2,2)*D(3,2)*tt52
tt350 = 1.0E+0*D(2,3)*D(3,3)*tt58
tt351 = 1.0E+0*lam(1,1)*(1.0E+0*D(2,3)*D(3,3)+1.0E+0*D(2,2)*D(3,2&
&)+1.0E+0*D(2,1)*D(3,1))*tt74
tt352 = volume(1,1)*(tt351+miu(1,1)*(tt350+tt349+tt348+tt347+tt34&
&6+tt345+2.0E+0*D(2,3)*D(3,3)*tt21+1.0E+0*tt171*tt217+1.0E+0*tt168&
&*tt216+2.0E+0*D(2,2)*D(3,2)*tt17+1.0E+0*tt165*tt215+2.0E+0*D(2,1)&
&*D(3,1)*tt14)+1.0E+0*lam(1,1)*tt162*tt214)
tt353 = 2.0E+0*D(2,1)*D(3,1)*tt13*tt42
tt354 = 2.0E+0*D(2,2)*D(3,2)*tt16*tt48
tt355 = 2.0E+0*D(2,3)*D(3,3)*tt20*tt54
tt356 = volume(1,1)*(miu(1,1)*(1.0E+0*tt171*tt232+1.0E+0*tt168*tt&
&231+tt355+1.0E+0*tt165*tt229+tt354+tt353)+1.0E+0*lam(1,1)*tt162*t&
&t226)
tt357 = 2.0E+0*D(2,1)*D(3,1)*tt13*tt44
tt358 = 2.0E+0*D(2,2)*D(3,2)*tt16*tt50
tt359 = 2.0E+0*D(2,3)*D(3,3)*tt20*tt56
tt360 = volume(1,1)*(1.0E+0*lam(1,1)*tt162*tt240+miu(1,1)*(1.0E+0&
&*tt171*tt239+1.0E+0*tt168*tt238+tt359+1.0E+0*tt165*tt236+tt358+tt&
&357))
tt361 = volume(1,1)*(1.0E+0*lam(1,1)*tt186*tt212+miu(1,1)*(1.0E+0&
&*tt198*tt211+1.0E+0*tt195*tt208+2.0E+0*tt335*tt54*tt56+1.0E+0*tt1&
&91*tt204+2.0E+0*tt334*tt48*tt50+2.0E+0*tt333*tt42*tt44))
tt362 = volume(1,1)*(miu(1,1)*(1.0E+0*tt198*tt217+1.0E+0*tt195*tt&
&216+tt355+1.0E+0*tt191*tt215+tt354+tt353)+1.0E+0*lam(1,1)*tt214*t&
&t186)
tt363 = volume(1,1)*(tt351+miu(1,1)*(tt350+tt349+tt348+tt347+tt34&
&6+tt345+2.0E+0*D(2,3)*D(3,3)*tt55+1.0E+0*tt198*tt232+1.0E+0*tt195&
&*tt231+2.0E+0*D(2,2)*D(3,2)*tt49+1.0E+0*tt191*tt229+2.0E+0*D(2,1)&
&*D(3,1)*tt43)+1.0E+0*lam(1,1)*tt186*tt226)
tt364 = 2.0E+0*D(2,1)*D(3,1)*tt42*tt44
tt365 = 2.0E+0*D(2,2)*D(3,2)*tt48*tt50
tt366 = 2.0E+0*D(2,3)*D(3,3)*tt54*tt56
tt367 = volume(1,1)*(1.0E+0*lam(1,1)*tt186*tt240+miu(1,1)*(1.0E+0&
&*tt198*tt239+1.0E+0*tt195*tt238+tt366+1.0E+0*tt191*tt236+tt365+tt&
&364))
tt368 = volume(1,1)*(1.0E+0*lam(1,1)*tt214*tt212+miu(1,1)*(1.0E+0&
&*tt217*tt211+1.0E+0*tt216*tt208+tt359+1.0E+0*tt215*tt204+tt358+tt&
&357))
tt369 = volume(1,1)*(1.0E+0*lam(1,1)*tt226*tt212+miu(1,1)*(1.0E+0&
&*tt232*tt211+1.0E+0*tt231*tt208+tt366+1.0E+0*tt229*tt204+tt365+tt&
&364))
tt370 = volume(1,1)*(tt351+miu(1,1)*(tt350+2.0E+0*D(2,3)*D(3,3)*t&
&t57+tt349+2.0E+0*D(2,2)*D(3,2)*tt51+tt348+2.0E+0*D(2,1)*D(3,1)*tt&
&45+tt347+tt346+tt345+1.0E+0*tt211*tt239+1.0E+0*tt208*tt238+1.0E+0&
&*tt204*tt236)+1.0E+0*lam(1,1)*tt212*tt240)
tt371 = D(3,1)**2
tt372 = D(3,2)**2
tt373 = D(3,3)**2
tt374 = 2.0E+0*D(3,1)*D(3,2)*tt36
tt375 = 2.0E+0*D(3,1)*D(3,3)*tt38
tt376 = 2.0E+0*D(3,2)*D(3,3)*tt40
tt377 = 1.0E+0*tt371*tt46
tt378 = 1.0E+0*tt372*tt52
tt379 = 1.0E+0*tt373*tt58
tt380 = 1.0E+0*lam(1,1)*(1.0E+0*tt373+1.0E+0*tt372+1.0E+0*tt371)*&
&tt74
tt381 = volume(1,1)*(miu(1,1)*(1.0E+0*tt217*tt232+1.0E+0*tt216*tt&
&231+2.0E+0*tt373*tt20*tt54+1.0E+0*tt215*tt229+2.0E+0*tt372*tt16*t&
&t48+2.0E+0*tt371*tt13*tt42)+1.0E+0*lam(1,1)*tt214*tt226)
tt382 = volume(1,1)*(1.0E+0*lam(1,1)*tt214*tt240+miu(1,1)*(1.0E+0&
&*tt217*tt239+1.0E+0*tt216*tt238+2.0E+0*tt373*tt20*tt56+1.0E+0*tt2&
&15*tt236+2.0E+0*tt372*tt16*tt50+2.0E+0*tt371*tt13*tt44))
tt383 = volume(1,1)*(1.0E+0*lam(1,1)*tt226*tt240+miu(1,1)*(1.0E+0&
&*tt232*tt239+1.0E+0*tt231*tt238+2.0E+0*tt373*tt54*tt56+1.0E+0*tt2&
&29*tt236+2.0E+0*tt372*tt48*tt50+2.0E+0*tt371*tt42*tt44))
hes(1,1) = volume(1,1)*(tt75+miu(1,1)*(tt59+tt53+tt47+tt41+tt39+t&
&t37+1.0E+0*tt23**2+1.0E+0*tt22**2+2.0E+0*tt19*tt21+1.0E+0*tt18**2&
&+2.0E+0*tt15*tt17+2.0E+0*tt12*tt14)+1.0E+0*lam(1,1)*tt11**2)
hes(1,2) = tt80
hes(1,3) = tt85
hes(1,4) = tt121
hes(1,5) = tt141
hes(1,6) = tt161
hes(1,7) = tt185
hes(1,8) = tt199
hes(1,9) = tt213
hes(1,10) = tt225
hes(1,11) = tt233
hes(1,12) = tt241
hes(2,1) = tt80
hes(2,2) = volume(1,1)*(tt75+miu(1,1)*(tt59+tt53+tt47+tt41+tt39+t&
&t37+1.0E+0*tt79**2+1.0E+0*tt78**2+2.0E+0*tt19*tt55+1.0E+0*tt77**2&
&+2.0E+0*tt15*tt49+2.0E+0*tt12*tt43)+1.0E+0*lam(1,1)*tt76**2)
hes(2,3) = tt242
hes(2,4) = tt243
hes(2,5) = tt244
hes(2,6) = tt248
hes(2,7) = tt249
hes(2,8) = tt250
hes(2,9) = tt254
hes(2,10) = tt255
hes(2,11) = tt256
hes(2,12) = tt260
hes(3,1) = tt85
hes(3,2) = tt242
hes(3,3) = volume(1,1)*(1.0E+0*lam(1,1)*tt84**2+tt75+miu(1,1)*(1.&
&0E+0*tt83**2+1.0E+0*tt82**2+tt59+2.0E+0*tt19*tt57+1.0E+0*tt81**2+&
&tt53+2.0E+0*tt15*tt51+tt47+2.0E+0*tt12*tt45+tt41+tt39+tt37))
hes(3,4) = tt261
hes(3,5) = tt262
hes(3,6) = tt263
hes(3,7) = tt264
hes(3,8) = tt265
hes(3,9) = tt266
hes(3,10) = tt267
hes(3,11) = tt268
hes(3,12) = tt269
hes(4,1) = tt121
hes(4,2) = tt243
hes(4,3) = tt261
hes(4,4) = volume(1,1)*(tt279+miu(1,1)*(tt278+tt277+tt276+tt275+t&
&t274+tt273+1.0E+0*tt101**2+1.0E+0*tt96**2+2.0E+0*tt272*tt21+1.0E+&
&0*tt91**2+2.0E+0*tt271*tt17+2.0E+0*tt270*tt14)+1.0E+0*lam(1,1)*tt&
&86**2)
hes(4,5) = tt280
hes(4,6) = tt281
hes(4,7) = tt289
hes(4,8) = tt293
hes(4,9) = tt297
hes(4,10) = tt305
hes(4,11) = tt309
hes(4,12) = tt313
hes(5,1) = tt141
hes(5,2) = tt244
hes(5,3) = tt262
hes(5,4) = tt280
hes(5,5) = volume(1,1)*(tt279+miu(1,1)*(tt278+tt277+tt276+tt275+t&
&t274+tt273+1.0E+0*tt140**2+1.0E+0*tt135**2+2.0E+0*tt272*tt55+1.0E&
&+0*tt129**2+2.0E+0*tt271*tt49+2.0E+0*tt270*tt43)+1.0E+0*lam(1,1)*&
&tt122**2)
hes(5,6) = tt314
hes(5,7) = tt315
hes(5,8) = tt316
hes(5,9) = tt320
hes(5,10) = tt321
hes(5,11) = tt322
hes(5,12) = tt326
hes(6,1) = tt161
hes(6,2) = tt248
hes(6,3) = tt263
hes(6,4) = tt281
hes(6,5) = tt314
hes(6,6) = volume(1,1)*(1.0E+0*lam(1,1)*tt160**2+tt279+miu(1,1)*(&
&1.0E+0*tt159**2+1.0E+0*tt154**2+tt278+2.0E+0*tt272*tt57+1.0E+0*tt&
&148**2+tt277+2.0E+0*tt271*tt51+tt276+2.0E+0*tt270*tt45+tt275+tt27&
&4+tt273))
hes(6,7) = tt327
hes(6,8) = tt328
hes(6,9) = tt329
hes(6,10) = tt330
hes(6,11) = tt331
hes(6,12) = tt332
hes(7,1) = tt185
hes(7,2) = tt249
hes(7,3) = tt264
hes(7,4) = tt289
hes(7,5) = tt315
hes(7,6) = tt327
hes(7,7) = volume(1,1)*(tt342+miu(1,1)*(tt341+tt340+tt339+tt338+t&
&t337+tt336+1.0E+0*tt171**2+1.0E+0*tt168**2+2.0E+0*tt335*tt21+1.0E&
&+0*tt165**2+2.0E+0*tt334*tt17+2.0E+0*tt333*tt14)+1.0E+0*lam(1,1)*&
&tt162**2)
hes(7,8) = tt343
hes(7,9) = tt344
hes(7,10) = tt352
hes(7,11) = tt356
hes(7,12) = tt360
hes(8,1) = tt199
hes(8,2) = tt250
hes(8,3) = tt265
hes(8,4) = tt293
hes(8,5) = tt316
hes(8,6) = tt328
hes(8,7) = tt343
hes(8,8) = volume(1,1)*(tt342+miu(1,1)*(tt341+tt340+tt339+tt338+t&
&t337+tt336+1.0E+0*tt198**2+1.0E+0*tt195**2+2.0E+0*tt335*tt55+1.0E&
&+0*tt191**2+2.0E+0*tt334*tt49+2.0E+0*tt333*tt43)+1.0E+0*lam(1,1)*&
&tt186**2)
hes(8,9) = tt361
hes(8,10) = tt362
hes(8,11) = tt363
hes(8,12) = tt367
hes(9,1) = tt213
hes(9,2) = tt254
hes(9,3) = tt266
hes(9,4) = tt297
hes(9,5) = tt320
hes(9,6) = tt329
hes(9,7) = tt344
hes(9,8) = tt361
hes(9,9) = volume(1,1)*(1.0E+0*lam(1,1)*tt212**2+tt342+miu(1,1)*(&
&1.0E+0*tt211**2+1.0E+0*tt208**2+tt341+2.0E+0*tt335*tt57+1.0E+0*tt&
&204**2+tt340+2.0E+0*tt334*tt51+tt339+2.0E+0*tt333*tt45+tt338+tt33&
&7+tt336))
hes(9,10) = tt368
hes(9,11) = tt369
hes(9,12) = tt370
hes(10,1) = tt225
hes(10,2) = tt255
hes(10,3) = tt267
hes(10,4) = tt305
hes(10,5) = tt321
hes(10,6) = tt330
hes(10,7) = tt352
hes(10,8) = tt362
hes(10,9) = tt368
hes(10,10) = volume(1,1)*(tt380+miu(1,1)*(tt379+tt378+tt377+tt376&
&+tt375+tt374+1.0E+0*tt217**2+1.0E+0*tt216**2+2.0E+0*tt373*tt21+1.&
&0E+0*tt215**2+2.0E+0*tt372*tt17+2.0E+0*tt371*tt14)+1.0E+0*lam(1,1&
&)*tt214**2)
hes(10,11) = tt381
hes(10,12) = tt382
hes(11,1) = tt233
hes(11,2) = tt256
hes(11,3) = tt268
hes(11,4) = tt309
hes(11,5) = tt322
hes(11,6) = tt331
hes(11,7) = tt356
hes(11,8) = tt363
hes(11,9) = tt369
hes(11,10) = tt381
hes(11,11) = volume(1,1)*(tt380+miu(1,1)*(tt379+tt378+tt377+tt376&
&+tt375+tt374+1.0E+0*tt232**2+1.0E+0*tt231**2+2.0E+0*tt373*tt55+1.&
&0E+0*tt229**2+2.0E+0*tt372*tt49+2.0E+0*tt371*tt43)+1.0E+0*lam(1,1&
&)*tt226**2)
hes(11,12) = tt383
hes(12,1) = tt241
hes(12,2) = tt260
hes(12,3) = tt269
hes(12,4) = tt313
hes(12,5) = tt326
hes(12,6) = tt332
hes(12,7) = tt360
hes(12,8) = tt367
hes(12,9) = tt370
hes(12,10) = tt382
hes(12,11) = tt383
hes(12,12) = volume(1,1)*(1.0E+0*lam(1,1)*tt240**2+tt380+miu(1,1)&
&*(1.0E+0*tt239**2+1.0E+0*tt238**2+tt379+2.0E+0*tt373*tt57+1.0E+0*&
&tt236**2+tt378+2.0E+0*tt372*tt51+tt377+2.0E+0*tt371*tt45+tt376+tt&
&375+tt374))
END 
SUBROUTINE tet_linear(val, X, D, volume, lam, miu) 
IMPLICIT NONE 
REAL(KIND=8) val(1, 1) 
REAL(KIND=8) X(3, 4) 
REAL(KIND=8) D(3, 3) 
REAL(KIND=8) volume(1, 1) 
REAL(KIND=8) lam(1, 1) 
REAL(KIND=8) miu(1, 1) 
REAL(KIND=8)  tt1 
REAL(KIND=8)  tt2 
REAL(KIND=8)  tt3 
tt1 = -X(1,1)
tt2 = -X(2,1)
tt3 = -X(3,1)
val(1,1) = volume(1,1)*(5.0E-1*lam(1,1)*(5.0E-1*(2*D(3,3)*(X(3,4)&
&+tt3)+2*D(2,3)*(X(3,3)+tt3)+2*D(1,3)*(X(3,2)+tt3))+5.0E-1*(2*(X(2&
&,4)+tt2)*D(3,2)+2*D(2,2)*(X(2,3)+tt2)+2*D(1,2)*(X(2,2)+tt2))+5.0E&
&-1*(2*(X(1,4)+tt1)*D(3,1)+2*(X(1,3)+tt1)*D(2,1)+2*D(1,1)*(X(1,2)+&
&tt1))-3)**2+miu(1,1)*((1.0E+0*D(3,3)*X(3,4)+1.0E+0*D(2,3)*X(3,3)-&
&1.0E+0*X(3,1)*D(3,3)+1.0E+0*D(1,3)*X(3,2)-1.0E+0*D(2,3)*X(3,1)-1.&
&0E+0*D(1,3)*X(3,1)-1)**2+5.0E-1*(D(3,2)*X(3,4)+D(2,2)*X(3,3)+X(2,&
&4)*D(3,3)-X(2,1)*D(3,3)+D(1,2)*X(3,2)-X(3,1)*D(3,2)-D(2,2)*X(3,1)&
&-D(1,2)*X(3,1)+D(2,3)*X(2,3)-X(2,1)*D(2,3)+D(1,3)*X(2,2)-D(1,3)*X&
&(2,1))**2+5.0E-1*(D(3,1)*X(3,4)+D(2,1)*X(3,3)+X(1,4)*D(3,3)-X(1,1&
&)*D(3,3)+D(1,1)*X(3,2)-D(3,1)*X(3,1)-D(2,1)*X(3,1)-D(1,1)*X(3,1)+&
&X(1,3)*D(2,3)-X(1,1)*D(2,3)+X(1,2)*D(1,3)-X(1,1)*D(1,3))**2+(1.0E&
&+0*X(2,4)*D(3,2)-1.0E+0*X(2,1)*D(3,2)+1.0E+0*D(2,2)*X(2,3)+1.0E+0&
&*D(1,2)*X(2,2)-1.0E+0*X(2,1)*D(2,2)-1.0E+0*D(1,2)*X(2,1)-1)**2+5.&
&0E-1*(X(1,4)*D(3,2)-X(1,1)*D(3,2)+X(2,4)*D(3,1)-X(2,1)*D(3,1)+D(2&
&,1)*X(2,3)+D(1,1)*X(2,2)+X(1,3)*D(2,2)-X(1,1)*D(2,2)-D(2,1)*X(2,1&
&)-D(1,1)*X(2,1)+D(1,2)*X(1,2)-X(1,1)*D(1,2))**2+(1.0E+0*X(1,4)*D(&
&3,1)-1.0E+0*X(1,1)*D(3,1)+1.0E+0*X(1,3)*D(2,1)-1.0E+0*X(1,1)*D(2,&
&1)+1.0E+0*D(1,1)*X(1,2)-1.0E+0*D(1,1)*X(1,1)-1)**2))
END 
SUBROUTINE tet_linear_jac(jac, X, D, volume, lam, miu) 
IMPLICIT NONE 
REAL(KIND=8) jac(1, 12) 
REAL(KIND=8) X(3, 4) 
REAL(KIND=8) D(3, 3) 
REAL(KIND=8) volume(1, 1) 
REAL(KIND=8) lam(1, 1) 
REAL(KIND=8) miu(1, 1) 
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
tt1 = 1.0E+0*X(1,4)*D(3,1)-1.0E+0*X(1,1)*D(3,1)+1.0E+0*X(1,3)*D(2&
&,1)-1.0E+0*X(1,1)*D(2,1)+1.0E+0*D(1,1)*X(1,2)-1.0E+0*D(1,1)*X(1,1&
&)-1
tt2 = (-D(3,2))-D(2,2)-D(1,2)
tt3 = X(1,4)*D(3,2)-X(1,1)*D(3,2)+X(2,4)*D(3,1)-X(2,1)*D(3,1)+D(2&
&,1)*X(2,3)+D(1,1)*X(2,2)+X(1,3)*D(2,2)-X(1,1)*D(2,2)-D(2,1)*X(2,1&
&)-D(1,1)*X(2,1)+D(1,2)*X(1,2)-X(1,1)*D(1,2)
tt4 = (-D(3,3))-D(2,3)-D(1,3)
tt5 = D(3,1)*X(3,4)+D(2,1)*X(3,3)+X(1,4)*D(3,3)-X(1,1)*D(3,3)+D(1&
&,1)*X(3,2)-D(3,1)*X(3,1)-D(2,1)*X(3,1)-D(1,1)*X(3,1)+X(1,3)*D(2,3&
&)-X(1,1)*D(2,3)+X(1,2)*D(1,3)-X(1,1)*D(1,3)
tt6 = -X(1,1)
tt7 = -X(2,1)
tt8 = -X(3,1)
tt9 = 5.0E-1*(2*D(3,3)*(X(3,4)+tt8)+2*D(2,3)*(X(3,3)+tt8)+2*D(1,3&
&)*(X(3,2)+tt8))+5.0E-1*(2*(X(2,4)+tt7)*D(3,2)+2*D(2,2)*(X(2,3)+tt&
&7)+2*D(1,2)*(X(2,2)+tt7))+5.0E-1*(2*(X(1,4)+tt6)*D(3,1)+2*(X(1,3)&
&+tt6)*D(2,1)+2*D(1,1)*(X(1,2)+tt6))-3
tt10 = (-D(3,1))-D(2,1)-D(1,1)
tt11 = 1.0E+0*X(2,4)*D(3,2)-1.0E+0*X(2,1)*D(3,2)+1.0E+0*D(2,2)*X(&
&2,3)+1.0E+0*D(1,2)*X(2,2)-1.0E+0*X(2,1)*D(2,2)-1.0E+0*D(1,2)*X(2,&
&1)-1
tt12 = D(3,2)*X(3,4)+D(2,2)*X(3,3)+X(2,4)*D(3,3)-X(2,1)*D(3,3)+D(&
&1,2)*X(3,2)-X(3,1)*D(3,2)-D(2,2)*X(3,1)-D(1,2)*X(3,1)+D(2,3)*X(2,&
&3)-X(2,1)*D(2,3)+D(1,3)*X(2,2)-D(1,3)*X(2,1)
tt13 = 1.0E+0*D(3,3)*X(3,4)+1.0E+0*D(2,3)*X(3,3)-1.0E+0*X(3,1)*D(&
&3,3)+1.0E+0*D(1,3)*X(3,2)-1.0E+0*D(2,3)*X(3,1)-1.0E+0*D(1,3)*X(3,&
&1)-1
jac(1,1) = volume(1,1)*(5.0E-1*lam(1,1)*((-2*D(3,1))-2*D(2,1)-2*D&
&(1,1))*tt9+miu(1,1)*(1.0E+0*tt4*tt5+1.0E+0*tt2*tt3+2*((-1.0E+0*D(&
&3,1))-1.0E+0*D(2,1)-1.0E+0*D(1,1))*tt1))
jac(1,2) = volume(1,1)*(5.0E-1*lam(1,1)*((-2*D(3,2))-2*D(2,2)-2*D&
&(1,2))*tt9+miu(1,1)*(1.0E+0*tt4*tt12+2*((-1.0E+0*D(3,2))-1.0E+0*D&
&(2,2)-1.0E+0*D(1,2))*tt11+1.0E+0*tt10*tt3))
jac(1,3) = volume(1,1)*(5.0E-1*lam(1,1)*((-2*D(3,3))-2*D(2,3)-2*D&
&(1,3))*tt9+miu(1,1)*(2*((-1.0E+0*D(3,3))-1.0E+0*D(2,3)-1.0E+0*D(1&
&,3))*tt13+1.0E+0*tt2*tt12+1.0E+0*tt10*tt5))
jac(1,4) = volume(1,1)*(1.0E+0*D(1,1)*lam(1,1)*tt9+miu(1,1)*(1.0E&
&+0*D(1,3)*tt5+1.0E+0*D(1,2)*tt3+2.0E+0*D(1,1)*tt1))
jac(1,5) = volume(1,1)*(1.0E+0*lam(1,1)*D(1,2)*tt9+miu(1,1)*(1.0E&
&+0*D(1,3)*tt12+2.0E+0*D(1,2)*tt11+1.0E+0*D(1,1)*tt3))
jac(1,6) = volume(1,1)*(1.0E+0*lam(1,1)*D(1,3)*tt9+miu(1,1)*(2.0E&
&+0*D(1,3)*tt13+1.0E+0*D(1,2)*tt12+1.0E+0*D(1,1)*tt5))
jac(1,7) = volume(1,1)*(1.0E+0*lam(1,1)*D(2,1)*tt9+miu(1,1)*(1.0E&
&+0*D(2,3)*tt5+1.0E+0*D(2,2)*tt3+2.0E+0*D(2,1)*tt1))
jac(1,8) = volume(1,1)*(1.0E+0*lam(1,1)*D(2,2)*tt9+miu(1,1)*(1.0E&
&+0*D(2,3)*tt12+2.0E+0*D(2,2)*tt11+1.0E+0*D(2,1)*tt3))
jac(1,9) = volume(1,1)*(1.0E+0*lam(1,1)*D(2,3)*tt9+miu(1,1)*(2.0E&
&+0*D(2,3)*tt13+1.0E+0*D(2,2)*tt12+1.0E+0*D(2,1)*tt5))
jac(1,10) = volume(1,1)*(1.0E+0*lam(1,1)*D(3,1)*tt9+miu(1,1)*(1.0&
&E+0*D(3,3)*tt5+1.0E+0*D(3,2)*tt3+2.0E+0*D(3,1)*tt1))
jac(1,11) = volume(1,1)*(1.0E+0*lam(1,1)*D(3,2)*tt9+miu(1,1)*(1.0&
&E+0*D(3,3)*tt12+2.0E+0*D(3,2)*tt11+1.0E+0*D(3,1)*tt3))
jac(1,12) = volume(1,1)*(1.0E+0*lam(1,1)*D(3,3)*tt9+miu(1,1)*(2.0&
&E+0*D(3,3)*tt13+1.0E+0*D(3,2)*tt12+1.0E+0*D(3,1)*tt5))
END 
SUBROUTINE tet_linear_hes(hes, X, D, volume, lam, miu) 
IMPLICIT NONE 
REAL(KIND=8) hes(12, 12) 
REAL(KIND=8) X(3, 4) 
REAL(KIND=8) D(3, 3) 
REAL(KIND=8) volume(1, 1) 
REAL(KIND=8) lam(1, 1) 
REAL(KIND=8) miu(1, 1) 
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
tt1 = (-2*D(3,1))-2*D(2,1)-2*D(1,1)
tt2 = (-1.0E+0*D(3,1))-1.0E+0*D(2,1)-1.0E+0*D(1,1)
tt3 = (-D(3,2))-D(2,2)-D(1,2)
tt4 = 1.0E+0*tt3**2
tt5 = (-D(3,3))-D(2,3)-D(1,3)
tt6 = 1.0E+0*tt5**2
tt7 = (-2*D(3,2))-2*D(2,2)-2*D(1,2)
tt8 = (-D(3,1))-D(2,1)-D(1,1)
tt9 = volume(1,1)*(1.0E+0*miu(1,1)*tt8*tt3+2.5E-1*lam(1,1)*tt1*tt&
&7)
tt10 = (-2*D(3,3))-2*D(2,3)-2*D(1,3)
tt11 = volume(1,1)*(1.0E+0*miu(1,1)*tt8*tt5+2.5E-1*lam(1,1)*tt1*t&
&t10)
tt12 = 1.0E+0*D(1,2)*tt3
tt13 = 1.0E+0*D(1,3)*tt5
tt14 = volume(1,1)*(miu(1,1)*(tt13+tt12+2.0E+0*D(1,1)*tt2)+5.0E-1&
&*D(1,1)*lam(1,1)*tt1)
tt15 = volume(1,1)*(1.0E+0*D(1,1)*miu(1,1)*tt3+5.0E-1*lam(1,1)*D(&
&1,2)*tt1)
tt16 = volume(1,1)*(1.0E+0*D(1,1)*miu(1,1)*tt5+5.0E-1*lam(1,1)*D(&
&1,3)*tt1)
tt17 = 1.0E+0*D(2,2)*tt3
tt18 = 1.0E+0*D(2,3)*tt5
tt19 = volume(1,1)*(miu(1,1)*(tt18+tt17+2.0E+0*D(2,1)*tt2)+5.0E-1&
&*lam(1,1)*D(2,1)*tt1)
tt20 = volume(1,1)*(1.0E+0*miu(1,1)*D(2,1)*tt3+5.0E-1*lam(1,1)*D(&
&2,2)*tt1)
tt21 = volume(1,1)*(1.0E+0*miu(1,1)*D(2,1)*tt5+5.0E-1*lam(1,1)*D(&
&2,3)*tt1)
tt22 = 1.0E+0*tt3*D(3,2)
tt23 = 1.0E+0*tt5*D(3,3)
tt24 = volume(1,1)*(miu(1,1)*(tt23+tt22+2.0E+0*tt2*D(3,1))+5.0E-1&
&*lam(1,1)*tt1*D(3,1))
tt25 = volume(1,1)*(5.0E-1*lam(1,1)*tt1*D(3,2)+1.0E+0*miu(1,1)*D(&
&3,1)*tt3)
tt26 = volume(1,1)*(5.0E-1*lam(1,1)*tt1*D(3,3)+1.0E+0*miu(1,1)*D(&
&3,1)*tt5)
tt27 = 1.0E+0*tt8**2
tt28 = (-1.0E+0*D(3,2))-1.0E+0*D(2,2)-1.0E+0*D(1,2)
tt29 = volume(1,1)*(1.0E+0*miu(1,1)*tt3*tt5+2.5E-1*lam(1,1)*tt7*t&
&t10)
tt30 = volume(1,1)*(5.0E-1*D(1,1)*lam(1,1)*tt7+1.0E+0*miu(1,1)*D(&
&1,2)*tt8)
tt31 = 1.0E+0*D(1,1)*tt8
tt32 = volume(1,1)*(miu(1,1)*(tt13+2.0E+0*D(1,2)*tt28+tt31)+5.0E-&
&1*lam(1,1)*D(1,2)*tt7)
tt33 = volume(1,1)*(1.0E+0*miu(1,1)*D(1,2)*tt5+5.0E-1*lam(1,1)*D(&
&1,3)*tt7)
tt34 = volume(1,1)*(5.0E-1*lam(1,1)*D(2,1)*tt7+1.0E+0*miu(1,1)*D(&
&2,2)*tt8)
tt35 = 1.0E+0*D(2,1)*tt8
tt36 = volume(1,1)*(miu(1,1)*(tt18+2.0E+0*D(2,2)*tt28+tt35)+5.0E-&
&1*lam(1,1)*D(2,2)*tt7)
tt37 = volume(1,1)*(1.0E+0*miu(1,1)*D(2,2)*tt5+5.0E-1*lam(1,1)*D(&
&2,3)*tt7)
tt38 = volume(1,1)*(1.0E+0*miu(1,1)*tt8*D(3,2)+5.0E-1*lam(1,1)*D(&
&3,1)*tt7)
tt39 = 1.0E+0*tt8*D(3,1)
tt40 = volume(1,1)*(miu(1,1)*(tt23+2.0E+0*tt28*D(3,2)+tt39)+5.0E-&
&1*lam(1,1)*tt7*D(3,2))
tt41 = volume(1,1)*(5.0E-1*lam(1,1)*tt7*D(3,3)+1.0E+0*miu(1,1)*D(&
&3,2)*tt5)
tt42 = (-1.0E+0*D(3,3))-1.0E+0*D(2,3)-1.0E+0*D(1,3)
tt43 = volume(1,1)*(5.0E-1*D(1,1)*lam(1,1)*tt10+1.0E+0*miu(1,1)*D&
&(1,3)*tt8)
tt44 = volume(1,1)*(5.0E-1*lam(1,1)*D(1,2)*tt10+1.0E+0*miu(1,1)*D&
&(1,3)*tt3)
tt45 = volume(1,1)*(5.0E-1*lam(1,1)*D(1,3)*tt10+miu(1,1)*(2.0E+0*&
&D(1,3)*tt42+tt12+tt31))
tt46 = volume(1,1)*(5.0E-1*lam(1,1)*D(2,1)*tt10+1.0E+0*miu(1,1)*D&
&(2,3)*tt8)
tt47 = volume(1,1)*(5.0E-1*lam(1,1)*D(2,2)*tt10+1.0E+0*miu(1,1)*D&
&(2,3)*tt3)
tt48 = volume(1,1)*(5.0E-1*lam(1,1)*D(2,3)*tt10+miu(1,1)*(2.0E+0*&
&D(2,3)*tt42+tt17+tt35))
tt49 = volume(1,1)*(1.0E+0*miu(1,1)*tt8*D(3,3)+5.0E-1*lam(1,1)*D(&
&3,1)*tt10)
tt50 = volume(1,1)*(1.0E+0*miu(1,1)*tt3*D(3,3)+5.0E-1*lam(1,1)*D(&
&3,2)*tt10)
tt51 = volume(1,1)*(miu(1,1)*(2.0E+0*tt42*D(3,3)+tt22+tt39)+5.0E-&
&1*lam(1,1)*tt10*D(3,3))
tt52 = D(1,1)**2
tt53 = D(1,2)**2
tt54 = 1.0E+0*tt53
tt55 = D(1,3)**2
tt56 = 1.0E+0*tt55
tt57 = volume(1,1)*(1.0E+0*D(1,1)*miu(1,1)*D(1,2)+1.0E+0*D(1,1)*l&
&am(1,1)*D(1,2))
tt58 = volume(1,1)*(1.0E+0*D(1,1)*miu(1,1)*D(1,3)+1.0E+0*D(1,1)*l&
&am(1,1)*D(1,3))
tt59 = 1.0E+0*D(1,2)*D(2,2)
tt60 = 1.0E+0*D(1,3)*D(2,3)
tt61 = volume(1,1)*(miu(1,1)*(tt60+tt59+2.0E+0*D(1,1)*D(2,1))+1.0&
&E+0*D(1,1)*lam(1,1)*D(2,1))
tt62 = volume(1,1)*(1.0E+0*D(1,1)*lam(1,1)*D(2,2)+1.0E+0*miu(1,1)&
&*D(1,2)*D(2,1))
tt63 = volume(1,1)*(1.0E+0*D(1,1)*lam(1,1)*D(2,3)+1.0E+0*miu(1,1)&
&*D(1,3)*D(2,1))
tt64 = 1.0E+0*D(1,2)*D(3,2)
tt65 = 1.0E+0*D(1,3)*D(3,3)
tt66 = volume(1,1)*(miu(1,1)*(tt65+tt64+2.0E+0*D(1,1)*D(3,1))+1.0&
&E+0*D(1,1)*lam(1,1)*D(3,1))
tt67 = volume(1,1)*(1.0E+0*D(1,1)*lam(1,1)*D(3,2)+1.0E+0*miu(1,1)&
&*D(1,2)*D(3,1))
tt68 = volume(1,1)*(1.0E+0*D(1,1)*lam(1,1)*D(3,3)+1.0E+0*miu(1,1)&
&*D(1,3)*D(3,1))
tt69 = 1.0E+0*tt52
tt70 = volume(1,1)*(1.0E+0*miu(1,1)*D(1,2)*D(1,3)+1.0E+0*lam(1,1)&
&*D(1,2)*D(1,3))
tt71 = volume(1,1)*(1.0E+0*D(1,1)*miu(1,1)*D(2,2)+1.0E+0*lam(1,1)&
&*D(1,2)*D(2,1))
tt72 = 1.0E+0*D(1,1)*D(2,1)
tt73 = volume(1,1)*(miu(1,1)*(tt60+2.0E+0*D(1,2)*D(2,2)+tt72)+1.0&
&E+0*lam(1,1)*D(1,2)*D(2,2))
tt74 = volume(1,1)*(1.0E+0*lam(1,1)*D(1,2)*D(2,3)+1.0E+0*miu(1,1)&
&*D(1,3)*D(2,2))
tt75 = volume(1,1)*(1.0E+0*D(1,1)*miu(1,1)*D(3,2)+1.0E+0*lam(1,1)&
&*D(1,2)*D(3,1))
tt76 = 1.0E+0*D(1,1)*D(3,1)
tt77 = volume(1,1)*(miu(1,1)*(tt65+2.0E+0*D(1,2)*D(3,2)+tt76)+1.0&
&E+0*lam(1,1)*D(1,2)*D(3,2))
tt78 = volume(1,1)*(1.0E+0*lam(1,1)*D(1,2)*D(3,3)+1.0E+0*miu(1,1)&
&*D(1,3)*D(3,2))
tt79 = volume(1,1)*(1.0E+0*D(1,1)*miu(1,1)*D(2,3)+1.0E+0*lam(1,1)&
&*D(1,3)*D(2,1))
tt80 = volume(1,1)*(1.0E+0*miu(1,1)*D(1,2)*D(2,3)+1.0E+0*lam(1,1)&
&*D(1,3)*D(2,2))
tt81 = volume(1,1)*(miu(1,1)*(2.0E+0*D(1,3)*D(2,3)+tt59+tt72)+1.0&
&E+0*lam(1,1)*D(1,3)*D(2,3))
tt82 = volume(1,1)*(1.0E+0*D(1,1)*miu(1,1)*D(3,3)+1.0E+0*lam(1,1)&
&*D(1,3)*D(3,1))
tt83 = volume(1,1)*(1.0E+0*miu(1,1)*D(1,2)*D(3,3)+1.0E+0*lam(1,1)&
&*D(1,3)*D(3,2))
tt84 = volume(1,1)*(miu(1,1)*(2.0E+0*D(1,3)*D(3,3)+tt64+tt76)+1.0&
&E+0*lam(1,1)*D(1,3)*D(3,3))
tt85 = D(2,1)**2
tt86 = D(2,2)**2
tt87 = 1.0E+0*tt86
tt88 = D(2,3)**2
tt89 = 1.0E+0*tt88
tt90 = volume(1,1)*(1.0E+0*miu(1,1)*D(2,1)*D(2,2)+1.0E+0*lam(1,1)&
&*D(2,1)*D(2,2))
tt91 = volume(1,1)*(1.0E+0*miu(1,1)*D(2,1)*D(2,3)+1.0E+0*lam(1,1)&
&*D(2,1)*D(2,3))
tt92 = 1.0E+0*D(2,2)*D(3,2)
tt93 = 1.0E+0*D(2,3)*D(3,3)
tt94 = volume(1,1)*(miu(1,1)*(tt93+tt92+2.0E+0*D(2,1)*D(3,1))+1.0&
&E+0*lam(1,1)*D(2,1)*D(3,1))
tt95 = volume(1,1)*(1.0E+0*lam(1,1)*D(2,1)*D(3,2)+1.0E+0*miu(1,1)&
&*D(2,2)*D(3,1))
tt96 = volume(1,1)*(1.0E+0*lam(1,1)*D(2,1)*D(3,3)+1.0E+0*miu(1,1)&
&*D(2,3)*D(3,1))
tt97 = 1.0E+0*tt85
tt98 = volume(1,1)*(1.0E+0*miu(1,1)*D(2,2)*D(2,3)+1.0E+0*lam(1,1)&
&*D(2,2)*D(2,3))
tt99 = volume(1,1)*(1.0E+0*miu(1,1)*D(2,1)*D(3,2)+1.0E+0*lam(1,1)&
&*D(2,2)*D(3,1))
tt100 = 1.0E+0*D(2,1)*D(3,1)
tt101 = volume(1,1)*(miu(1,1)*(tt93+2.0E+0*D(2,2)*D(3,2)+tt100)+1&
&.0E+0*lam(1,1)*D(2,2)*D(3,2))
tt102 = volume(1,1)*(1.0E+0*lam(1,1)*D(2,2)*D(3,3)+1.0E+0*miu(1,1&
&)*D(2,3)*D(3,2))
tt103 = volume(1,1)*(1.0E+0*miu(1,1)*D(2,1)*D(3,3)+1.0E+0*lam(1,1&
&)*D(2,3)*D(3,1))
tt104 = volume(1,1)*(1.0E+0*miu(1,1)*D(2,2)*D(3,3)+1.0E+0*lam(1,1&
&)*D(2,3)*D(3,2))
tt105 = volume(1,1)*(miu(1,1)*(2.0E+0*D(2,3)*D(3,3)+tt92+tt100)+1&
&.0E+0*lam(1,1)*D(2,3)*D(3,3))
tt106 = D(3,1)**2
tt107 = D(3,2)**2
tt108 = 1.0E+0*tt107
tt109 = D(3,3)**2
tt110 = 1.0E+0*tt109
tt111 = volume(1,1)*(1.0E+0*miu(1,1)*D(3,1)*D(3,2)+1.0E+0*lam(1,1&
&)*D(3,1)*D(3,2))
tt112 = volume(1,1)*(1.0E+0*miu(1,1)*D(3,1)*D(3,3)+1.0E+0*lam(1,1&
&)*D(3,1)*D(3,3))
tt113 = 1.0E+0*tt106
tt114 = volume(1,1)*(1.0E+0*miu(1,1)*D(3,2)*D(3,3)+1.0E+0*lam(1,1&
&)*D(3,2)*D(3,3))
hes(1,1) = volume(1,1)*(miu(1,1)*(tt6+tt4+2*tt2**2)+2.5E-1*lam(1,&
&1)*tt1**2)
hes(1,2) = tt9
hes(1,3) = tt11
hes(1,4) = tt14
hes(1,5) = tt15
hes(1,6) = tt16
hes(1,7) = tt19
hes(1,8) = tt20
hes(1,9) = tt21
hes(1,10) = tt24
hes(1,11) = tt25
hes(1,12) = tt26
hes(2,1) = tt9
hes(2,2) = volume(1,1)*(miu(1,1)*(tt6+2*tt28**2+tt27)+2.5E-1*lam(&
&1,1)*tt7**2)
hes(2,3) = tt29
hes(2,4) = tt30
hes(2,5) = tt32
hes(2,6) = tt33
hes(2,7) = tt34
hes(2,8) = tt36
hes(2,9) = tt37
hes(2,10) = tt38
hes(2,11) = tt40
hes(2,12) = tt41
hes(3,1) = tt11
hes(3,2) = tt29
hes(3,3) = volume(1,1)*(miu(1,1)*(2*tt42**2+tt4+tt27)+2.5E-1*lam(&
&1,1)*tt10**2)
hes(3,4) = tt43
hes(3,5) = tt44
hes(3,6) = tt45
hes(3,7) = tt46
hes(3,8) = tt47
hes(3,9) = tt48
hes(3,10) = tt49
hes(3,11) = tt50
hes(3,12) = tt51
hes(4,1) = tt14
hes(4,2) = tt30
hes(4,3) = tt43
hes(4,4) = volume(1,1)*(miu(1,1)*(tt56+tt54+2.0E+0*tt52)+1.0E+0*t&
&t52*lam(1,1))
hes(4,5) = tt57
hes(4,6) = tt58
hes(4,7) = tt61
hes(4,8) = tt62
hes(4,9) = tt63
hes(4,10) = tt66
hes(4,11) = tt67
hes(4,12) = tt68
hes(5,1) = tt15
hes(5,2) = tt32
hes(5,3) = tt44
hes(5,4) = tt57
hes(5,5) = volume(1,1)*(miu(1,1)*(tt56+2.0E+0*tt53+tt69)+1.0E+0*l&
&am(1,1)*tt53)
hes(5,6) = tt70
hes(5,7) = tt71
hes(5,8) = tt73
hes(5,9) = tt74
hes(5,10) = tt75
hes(5,11) = tt77
hes(5,12) = tt78
hes(6,1) = tt16
hes(6,2) = tt33
hes(6,3) = tt45
hes(6,4) = tt58
hes(6,5) = tt70
hes(6,6) = volume(1,1)*(miu(1,1)*(2.0E+0*tt55+tt54+tt69)+1.0E+0*l&
&am(1,1)*tt55)
hes(6,7) = tt79
hes(6,8) = tt80
hes(6,9) = tt81
hes(6,10) = tt82
hes(6,11) = tt83
hes(6,12) = tt84
hes(7,1) = tt19
hes(7,2) = tt34
hes(7,3) = tt46
hes(7,4) = tt61
hes(7,5) = tt71
hes(7,6) = tt79
hes(7,7) = volume(1,1)*(miu(1,1)*(tt89+tt87+2.0E+0*tt85)+1.0E+0*l&
&am(1,1)*tt85)
hes(7,8) = tt90
hes(7,9) = tt91
hes(7,10) = tt94
hes(7,11) = tt95
hes(7,12) = tt96
hes(8,1) = tt20
hes(8,2) = tt36
hes(8,3) = tt47
hes(8,4) = tt62
hes(8,5) = tt73
hes(8,6) = tt80
hes(8,7) = tt90
hes(8,8) = volume(1,1)*(miu(1,1)*(tt89+2.0E+0*tt86+tt97)+1.0E+0*l&
&am(1,1)*tt86)
hes(8,9) = tt98
hes(8,10) = tt99
hes(8,11) = tt101
hes(8,12) = tt102
hes(9,1) = tt21
hes(9,2) = tt37
hes(9,3) = tt48
hes(9,4) = tt63
hes(9,5) = tt74
hes(9,6) = tt81
hes(9,7) = tt91
hes(9,8) = tt98
hes(9,9) = volume(1,1)*(miu(1,1)*(2.0E+0*tt88+tt87+tt97)+1.0E+0*l&
&am(1,1)*tt88)
hes(9,10) = tt103
hes(9,11) = tt104
hes(9,12) = tt105
hes(10,1) = tt24
hes(10,2) = tt38
hes(10,3) = tt49
hes(10,4) = tt66
hes(10,5) = tt75
hes(10,6) = tt82
hes(10,7) = tt94
hes(10,8) = tt99
hes(10,9) = tt103
hes(10,10) = volume(1,1)*(miu(1,1)*(tt110+tt108+2.0E+0*tt106)+1.0&
&E+0*lam(1,1)*tt106)
hes(10,11) = tt111
hes(10,12) = tt112
hes(11,1) = tt25
hes(11,2) = tt40
hes(11,3) = tt50
hes(11,4) = tt67
hes(11,5) = tt77
hes(11,6) = tt83
hes(11,7) = tt95
hes(11,8) = tt101
hes(11,9) = tt104
hes(11,10) = tt111
hes(11,11) = volume(1,1)*(miu(1,1)*(tt110+2.0E+0*tt107+tt113)+1.0&
&E+0*lam(1,1)*tt107)
hes(11,12) = tt114
hes(12,1) = tt26
hes(12,2) = tt41
hes(12,3) = tt51
hes(12,4) = tt68
hes(12,5) = tt78
hes(12,6) = tt84
hes(12,7) = tt96
hes(12,8) = tt102
hes(12,9) = tt105
hes(12,10) = tt112
hes(12,11) = tt114
hes(12,12) = volume(1,1)*(miu(1,1)*(2.0E+0*tt109+tt108+tt113)+1.0&
&E+0*lam(1,1)*tt109)
END 
SUBROUTINE tet_coro(val, X, D, R, volume, lam, miu) 
IMPLICIT NONE 
REAL(KIND=8) val(1, 1) 
REAL(KIND=8) X(3, 4) 
REAL(KIND=8) D(3, 3) 
REAL(KIND=8) R(3, 3) 
REAL(KIND=8) volume(1, 1) 
REAL(KIND=8) lam(1, 1) 
REAL(KIND=8) miu(1, 1) 
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
tt1 = -X(1,1)
tt2 = X(1,2)+tt1
tt3 = X(1,3)+tt1
tt4 = X(1,4)+tt1
tt5 = -X(2,1)
tt6 = X(2,2)+tt5
tt7 = X(2,3)+tt5
tt8 = X(2,4)+tt5
tt9 = -X(3,1)
tt10 = X(3,2)+tt9
tt11 = X(3,3)+tt9
tt12 = X(3,4)+tt9
val(1,1) = volume(1,1)*((lam(1,1)*((2*R(3,3)*(D(3,3)*tt12+D(2,3)*&
&tt11+D(1,3)*tt10)+2*R(2,3)*(tt8*D(3,3)+D(2,3)*tt7+D(1,3)*tt6)+2*R&
&(1,3)*(tt4*D(3,3)+tt3*D(2,3)+tt2*D(1,3)))/2.0E+0+(2*R(3,2)*(D(3,2&
&)*tt12+D(2,2)*tt11+D(1,2)*tt10)+2*R(2,2)*(tt8*D(3,2)+D(2,2)*tt7+D&
&(1,2)*tt6)+2*R(1,2)*(tt4*D(3,2)+tt3*D(2,2)+D(1,2)*tt2))/2.0E+0+(2&
&*R(3,1)*(D(3,1)*tt12+D(2,1)*tt11+D(1,1)*tt10)+2*R(2,1)*(tt8*D(3,1&
&)+D(2,1)*tt7+D(1,1)*tt6)+2*R(1,1)*(tt4*D(3,1)+tt3*D(2,1)+D(1,1)*t&
&t2))/2.0E+0-3)**2)/2.0E+0+miu(1,1)*((D(3,3)*R(3,3)*X(3,4)+D(2,3)*&
&R(3,3)*X(3,3)-X(3,1)*D(3,3)*R(3,3)+D(1,3)*X(3,2)*R(3,3)-D(2,3)*X(&
&3,1)*R(3,3)-D(1,3)*X(3,1)*R(3,3)+R(2,3)*X(2,4)*D(3,3)-X(2,1)*R(2,&
&3)*D(3,3)+R(1,3)*X(1,4)*D(3,3)-X(1,1)*R(1,3)*D(3,3)+D(2,3)*R(2,3)&
&*X(2,3)-X(2,1)*D(2,3)*R(2,3)+D(1,3)*X(2,2)*R(2,3)-D(1,3)*X(2,1)*R&
&(2,3)+R(1,3)*X(1,3)*D(2,3)-X(1,1)*R(1,3)*D(2,3)+X(1,2)*D(1,3)*R(1&
&,3)-X(1,1)*D(1,3)*R(1,3)-1)**2+(D(3,2)*R(3,3)*X(3,4)+R(3,2)*D(3,3&
&)*X(3,4)+D(2,2)*R(3,3)*X(3,3)+D(2,3)*R(3,2)*X(3,3)+D(1,2)*X(3,2)*&
&R(3,3)-X(3,1)*D(3,2)*R(3,3)-D(2,2)*X(3,1)*R(3,3)-D(1,2)*X(3,1)*R(&
&3,3)-X(3,1)*R(3,2)*D(3,3)+R(2,2)*X(2,4)*D(3,3)-X(2,1)*R(2,2)*D(3,&
&3)+R(1,2)*X(1,4)*D(3,3)-X(1,1)*R(1,2)*D(3,3)+D(1,3)*R(3,2)*X(3,2)&
&-D(2,3)*X(3,1)*R(3,2)-D(1,3)*X(3,1)*R(3,2)+R(2,3)*X(2,4)*D(3,2)-X&
&(2,1)*R(2,3)*D(3,2)+R(1,3)*X(1,4)*D(3,2)-X(1,1)*R(1,3)*D(3,2)+D(2&
&,2)*R(2,3)*X(2,3)+R(2,2)*D(2,3)*X(2,3)+D(1,2)*X(2,2)*R(2,3)-X(2,1&
&)*D(2,2)*R(2,3)-D(1,2)*X(2,1)*R(2,3)-X(2,1)*R(2,2)*D(2,3)+R(1,2)*&
&X(1,3)*D(2,3)-X(1,1)*R(1,2)*D(2,3)+D(1,3)*R(2,2)*X(2,2)-D(1,3)*X(&
&2,1)*R(2,2)+R(1,3)*X(1,3)*D(2,2)-X(1,1)*R(1,3)*D(2,2)+D(1,2)*X(1,&
&2)*R(1,3)-X(1,1)*D(1,2)*R(1,3)+R(1,2)*X(1,2)*D(1,3)-X(1,1)*R(1,2)&
&*D(1,3))**2/2.0E+0+(D(3,1)*R(3,3)*X(3,4)+R(3,1)*D(3,3)*X(3,4)+D(2&
&,1)*R(3,3)*X(3,3)+D(2,3)*R(3,1)*X(3,3)+D(1,1)*X(3,2)*R(3,3)-D(3,1&
&)*X(3,1)*R(3,3)-D(2,1)*X(3,1)*R(3,3)-D(1,1)*X(3,1)*R(3,3)-R(3,1)*&
&X(3,1)*D(3,3)+R(2,1)*X(2,4)*D(3,3)-R(2,1)*X(2,1)*D(3,3)+R(1,1)*X(&
&1,4)*D(3,3)-R(1,1)*X(1,1)*D(3,3)+D(1,3)*R(3,1)*X(3,2)-D(2,3)*R(3,&
&1)*X(3,1)-D(1,3)*R(3,1)*X(3,1)+R(2,3)*X(2,4)*D(3,1)-X(2,1)*R(2,3)&
&*D(3,1)+R(1,3)*X(1,4)*D(3,1)-X(1,1)*R(1,3)*D(3,1)+D(2,1)*R(2,3)*X&
&(2,3)+R(2,1)*D(2,3)*X(2,3)+D(1,1)*X(2,2)*R(2,3)-D(2,1)*X(2,1)*R(2&
&,3)-D(1,1)*X(2,1)*R(2,3)-R(2,1)*X(2,1)*D(2,3)+R(1,1)*X(1,3)*D(2,3&
&)-R(1,1)*X(1,1)*D(2,3)+D(1,3)*R(2,1)*X(2,2)-D(1,3)*R(2,1)*X(2,1)+&
&R(1,3)*X(1,3)*D(2,1)-X(1,1)*R(1,3)*D(2,1)+D(1,1)*X(1,2)*R(1,3)-D(&
&1,1)*X(1,1)*R(1,3)+R(1,1)*X(1,2)*D(1,3)-R(1,1)*X(1,1)*D(1,3))**2/&
&2.0E+0+(D(3,2)*R(3,2)*X(3,4)+D(2,2)*R(3,2)*X(3,3)+D(1,2)*R(3,2)*X&
&(3,2)-X(3,1)*D(3,2)*R(3,2)-D(2,2)*X(3,1)*R(3,2)-D(1,2)*X(3,1)*R(3&
&,2)+R(2,2)*X(2,4)*D(3,2)-X(2,1)*R(2,2)*D(3,2)+R(1,2)*X(1,4)*D(3,2&
&)-X(1,1)*R(1,2)*D(3,2)+D(2,2)*R(2,2)*X(2,3)+D(1,2)*R(2,2)*X(2,2)-&
&X(2,1)*D(2,2)*R(2,2)-D(1,2)*X(2,1)*R(2,2)+R(1,2)*X(1,3)*D(2,2)-X(&
&1,1)*R(1,2)*D(2,2)+D(1,2)*R(1,2)*X(1,2)-X(1,1)*D(1,2)*R(1,2)-1)**&
&2+(D(3,1)*R(3,2)*X(3,4)+R(3,1)*D(3,2)*X(3,4)+D(2,1)*R(3,2)*X(3,3)&
&+D(2,2)*R(3,1)*X(3,3)+D(1,1)*R(3,2)*X(3,2)+D(1,2)*R(3,1)*X(3,2)-D&
&(3,1)*X(3,1)*R(3,2)-D(2,1)*X(3,1)*R(3,2)-D(1,1)*X(3,1)*R(3,2)-R(3&
&,1)*X(3,1)*D(3,2)+R(2,1)*X(2,4)*D(3,2)-R(2,1)*X(2,1)*D(3,2)+R(1,1&
&)*X(1,4)*D(3,2)-R(1,1)*X(1,1)*D(3,2)-D(2,2)*R(3,1)*X(3,1)-D(1,2)*&
&R(3,1)*X(3,1)+R(2,2)*X(2,4)*D(3,1)-X(2,1)*R(2,2)*D(3,1)+R(1,2)*X(&
&1,4)*D(3,1)-X(1,1)*R(1,2)*D(3,1)+D(2,1)*R(2,2)*X(2,3)+R(2,1)*D(2,&
&2)*X(2,3)+D(1,1)*R(2,2)*X(2,2)+D(1,2)*R(2,1)*X(2,2)-D(2,1)*X(2,1)&
&*R(2,2)-D(1,1)*X(2,1)*R(2,2)-R(2,1)*X(2,1)*D(2,2)+R(1,1)*X(1,3)*D&
&(2,2)-R(1,1)*X(1,1)*D(2,2)-D(1,2)*R(2,1)*X(2,1)+R(1,2)*X(1,3)*D(2&
&,1)-X(1,1)*R(1,2)*D(2,1)+D(1,1)*R(1,2)*X(1,2)+R(1,1)*D(1,2)*X(1,2&
&)-D(1,1)*X(1,1)*R(1,2)-R(1,1)*X(1,1)*D(1,2))**2/2.0E+0+(D(3,1)*R(&
&3,1)*X(3,4)+D(2,1)*R(3,1)*X(3,3)+D(1,1)*R(3,1)*X(3,2)-D(3,1)*R(3,&
&1)*X(3,1)-D(2,1)*R(3,1)*X(3,1)-D(1,1)*R(3,1)*X(3,1)+R(2,1)*X(2,4)&
&*D(3,1)-R(2,1)*X(2,1)*D(3,1)+R(1,1)*X(1,4)*D(3,1)-R(1,1)*X(1,1)*D&
&(3,1)+D(2,1)*R(2,1)*X(2,3)+D(1,1)*R(2,1)*X(2,2)-D(2,1)*R(2,1)*X(2&
&,1)-D(1,1)*R(2,1)*X(2,1)+R(1,1)*X(1,3)*D(2,1)-R(1,1)*X(1,1)*D(2,1&
&)+D(1,1)*R(1,1)*X(1,2)-D(1,1)*R(1,1)*X(1,1)-1)**2))
END 
SUBROUTINE tet_coro_jac(jac, X, D, R, volume, lam, miu) 
IMPLICIT NONE 
REAL(KIND=8) jac(1, 12) 
REAL(KIND=8) X(3, 4) 
REAL(KIND=8) D(3, 3) 
REAL(KIND=8) R(3, 3) 
REAL(KIND=8) volume(1, 1) 
REAL(KIND=8) lam(1, 1) 
REAL(KIND=8) miu(1, 1) 
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
tt1 = D(3,1)*R(3,1)*X(3,4)+D(2,1)*R(3,1)*X(3,3)+D(1,1)*R(3,1)*X(3&
&,2)-D(3,1)*R(3,1)*X(3,1)-D(2,1)*R(3,1)*X(3,1)-D(1,1)*R(3,1)*X(3,1&
&)+R(2,1)*X(2,4)*D(3,1)-R(2,1)*X(2,1)*D(3,1)+R(1,1)*X(1,4)*D(3,1)-&
&R(1,1)*X(1,1)*D(3,1)+D(2,1)*R(2,1)*X(2,3)+D(1,1)*R(2,1)*X(2,2)-D(&
&2,1)*R(2,1)*X(2,1)-D(1,1)*R(2,1)*X(2,1)+R(1,1)*X(1,3)*D(2,1)-R(1,&
&1)*X(1,1)*D(2,1)+D(1,1)*R(1,1)*X(1,2)-D(1,1)*R(1,1)*X(1,1)-1
tt2 = D(3,1)*R(3,2)*X(3,4)+R(3,1)*D(3,2)*X(3,4)+D(2,1)*R(3,2)*X(3&
&,3)+D(2,2)*R(3,1)*X(3,3)+D(1,1)*R(3,2)*X(3,2)+D(1,2)*R(3,1)*X(3,2&
&)-D(3,1)*X(3,1)*R(3,2)-D(2,1)*X(3,1)*R(3,2)-D(1,1)*X(3,1)*R(3,2)-&
&R(3,1)*X(3,1)*D(3,2)+R(2,1)*X(2,4)*D(3,2)-R(2,1)*X(2,1)*D(3,2)+R(&
&1,1)*X(1,4)*D(3,2)-R(1,1)*X(1,1)*D(3,2)-D(2,2)*R(3,1)*X(3,1)-D(1,&
&2)*R(3,1)*X(3,1)+R(2,2)*X(2,4)*D(3,1)-X(2,1)*R(2,2)*D(3,1)+R(1,2)&
&*X(1,4)*D(3,1)-X(1,1)*R(1,2)*D(3,1)+D(2,1)*R(2,2)*X(2,3)+R(2,1)*D&
&(2,2)*X(2,3)+D(1,1)*R(2,2)*X(2,2)+D(1,2)*R(2,1)*X(2,2)-D(2,1)*X(2&
&,1)*R(2,2)-D(1,1)*X(2,1)*R(2,2)-R(2,1)*X(2,1)*D(2,2)+R(1,1)*X(1,3&
&)*D(2,2)-R(1,1)*X(1,1)*D(2,2)-D(1,2)*R(2,1)*X(2,1)+R(1,2)*X(1,3)*&
&D(2,1)-X(1,1)*R(1,2)*D(2,1)+D(1,1)*R(1,2)*X(1,2)+R(1,1)*D(1,2)*X(&
&1,2)-D(1,1)*X(1,1)*R(1,2)-R(1,1)*X(1,1)*D(1,2)
tt3 = D(3,2)*R(3,2)*X(3,4)+D(2,2)*R(3,2)*X(3,3)+D(1,2)*R(3,2)*X(3&
&,2)-X(3,1)*D(3,2)*R(3,2)-D(2,2)*X(3,1)*R(3,2)-D(1,2)*X(3,1)*R(3,2&
&)+R(2,2)*X(2,4)*D(3,2)-X(2,1)*R(2,2)*D(3,2)+R(1,2)*X(1,4)*D(3,2)-&
&X(1,1)*R(1,2)*D(3,2)+D(2,2)*R(2,2)*X(2,3)+D(1,2)*R(2,2)*X(2,2)-X(&
&2,1)*D(2,2)*R(2,2)-D(1,2)*X(2,1)*R(2,2)+R(1,2)*X(1,3)*D(2,2)-X(1,&
&1)*R(1,2)*D(2,2)+D(1,2)*R(1,2)*X(1,2)-X(1,1)*D(1,2)*R(1,2)-1
tt4 = D(3,1)*R(3,3)*X(3,4)+R(3,1)*D(3,3)*X(3,4)+D(2,1)*R(3,3)*X(3&
&,3)+D(2,3)*R(3,1)*X(3,3)+D(1,1)*X(3,2)*R(3,3)-D(3,1)*X(3,1)*R(3,3&
&)-D(2,1)*X(3,1)*R(3,3)-D(1,1)*X(3,1)*R(3,3)-R(3,1)*X(3,1)*D(3,3)+&
&R(2,1)*X(2,4)*D(3,3)-R(2,1)*X(2,1)*D(3,3)+R(1,1)*X(1,4)*D(3,3)-R(&
&1,1)*X(1,1)*D(3,3)+D(1,3)*R(3,1)*X(3,2)-D(2,3)*R(3,1)*X(3,1)-D(1,&
&3)*R(3,1)*X(3,1)+R(2,3)*X(2,4)*D(3,1)-X(2,1)*R(2,3)*D(3,1)+R(1,3)&
&*X(1,4)*D(3,1)-X(1,1)*R(1,3)*D(3,1)+D(2,1)*R(2,3)*X(2,3)+R(2,1)*D&
&(2,3)*X(2,3)+D(1,1)*X(2,2)*R(2,3)-D(2,1)*X(2,1)*R(2,3)-D(1,1)*X(2&
&,1)*R(2,3)-R(2,1)*X(2,1)*D(2,3)+R(1,1)*X(1,3)*D(2,3)-R(1,1)*X(1,1&
&)*D(2,3)+D(1,3)*R(2,1)*X(2,2)-D(1,3)*R(2,1)*X(2,1)+R(1,3)*X(1,3)*&
&D(2,1)-X(1,1)*R(1,3)*D(2,1)+D(1,1)*X(1,2)*R(1,3)-D(1,1)*X(1,1)*R(&
&1,3)+R(1,1)*X(1,2)*D(1,3)-R(1,1)*X(1,1)*D(1,3)
tt5 = D(3,2)*R(3,3)*X(3,4)+R(3,2)*D(3,3)*X(3,4)+D(2,2)*R(3,3)*X(3&
&,3)+D(2,3)*R(3,2)*X(3,3)+D(1,2)*X(3,2)*R(3,3)-X(3,1)*D(3,2)*R(3,3&
&)-D(2,2)*X(3,1)*R(3,3)-D(1,2)*X(3,1)*R(3,3)-X(3,1)*R(3,2)*D(3,3)+&
&R(2,2)*X(2,4)*D(3,3)-X(2,1)*R(2,2)*D(3,3)+R(1,2)*X(1,4)*D(3,3)-X(&
&1,1)*R(1,2)*D(3,3)+D(1,3)*R(3,2)*X(3,2)-D(2,3)*X(3,1)*R(3,2)-D(1,&
&3)*X(3,1)*R(3,2)+R(2,3)*X(2,4)*D(3,2)-X(2,1)*R(2,3)*D(3,2)+R(1,3)&
&*X(1,4)*D(3,2)-X(1,1)*R(1,3)*D(3,2)+D(2,2)*R(2,3)*X(2,3)+R(2,2)*D&
&(2,3)*X(2,3)+D(1,2)*X(2,2)*R(2,3)-X(2,1)*D(2,2)*R(2,3)-D(1,2)*X(2&
&,1)*R(2,3)-X(2,1)*R(2,2)*D(2,3)+R(1,2)*X(1,3)*D(2,3)-X(1,1)*R(1,2&
&)*D(2,3)+D(1,3)*R(2,2)*X(2,2)-D(1,3)*X(2,1)*R(2,2)+R(1,3)*X(1,3)*&
&D(2,2)-X(1,1)*R(1,3)*D(2,2)+D(1,2)*X(1,2)*R(1,3)-X(1,1)*D(1,2)*R(&
&1,3)+R(1,2)*X(1,2)*D(1,3)-X(1,1)*R(1,2)*D(1,3)
tt6 = D(3,3)*R(3,3)*X(3,4)+D(2,3)*R(3,3)*X(3,3)-X(3,1)*D(3,3)*R(3&
&,3)+D(1,3)*X(3,2)*R(3,3)-D(2,3)*X(3,1)*R(3,3)-D(1,3)*X(3,1)*R(3,3&
&)+R(2,3)*X(2,4)*D(3,3)-X(2,1)*R(2,3)*D(3,3)+R(1,3)*X(1,4)*D(3,3)-&
&X(1,1)*R(1,3)*D(3,3)+D(2,3)*R(2,3)*X(2,3)-X(2,1)*D(2,3)*R(2,3)+D(&
&1,3)*X(2,2)*R(2,3)-D(1,3)*X(2,1)*R(2,3)+R(1,3)*X(1,3)*D(2,3)-X(1,&
&1)*R(1,3)*D(2,3)+X(1,2)*D(1,3)*R(1,3)-X(1,1)*D(1,3)*R(1,3)-1
tt7 = (-D(3,1))-D(2,1)-D(1,1)
tt8 = (-D(3,2))-D(2,2)-D(1,2)
tt9 = (-D(3,3))-D(2,3)-D(1,3)
tt10 = -X(1,1)
tt11 = X(1,2)+tt10
tt12 = X(1,3)+tt10
tt13 = X(1,4)+tt10
tt14 = -X(2,1)
tt15 = X(2,2)+tt14
tt16 = X(2,3)+tt14
tt17 = X(2,4)+tt14
tt18 = -X(3,1)
tt19 = X(3,2)+tt18
tt20 = X(3,3)+tt18
tt21 = X(3,4)+tt18
tt22 = (2*R(3,3)*(D(3,3)*tt21+D(2,3)*tt20+D(1,3)*tt19)+2*R(2,3)*(&
&tt17*D(3,3)+D(2,3)*tt16+D(1,3)*tt15)+2*R(1,3)*(tt13*D(3,3)+tt12*D&
&(2,3)+tt11*D(1,3)))/2.0E+0+(2*R(3,2)*(D(3,2)*tt21+D(2,2)*tt20+D(1&
&,2)*tt19)+2*R(2,2)*(tt17*D(3,2)+D(2,2)*tt16+D(1,2)*tt15)+2*R(1,2)&
&*(tt13*D(3,2)+tt12*D(2,2)+D(1,2)*tt11))/2.0E+0+(2*R(3,1)*(D(3,1)*&
&tt21+D(2,1)*tt20+D(1,1)*tt19)+2*R(2,1)*(tt17*D(3,1)+D(2,1)*tt16+D&
&(1,1)*tt15)+2*R(1,1)*(tt13*D(3,1)+tt12*D(2,1)+D(1,1)*tt11))/2.0E+&
&0-3
jac(1,1) = volume(1,1)*(lam(1,1)*(R(1,3)*tt9+R(1,2)*tt8+R(1,1)*tt&
&7)*tt22+miu(1,1)*(2*((-R(1,3)*D(3,3))-R(1,3)*D(2,3)-D(1,3)*R(1,3)&
&)*tt6+((-R(1,2)*D(3,3))-R(1,3)*D(3,2)-R(1,2)*D(2,3)-R(1,3)*D(2,2)&
&-D(1,2)*R(1,3)-R(1,2)*D(1,3))*tt5+((-R(1,1)*D(3,3))-R(1,3)*D(3,1)&
&-R(1,1)*D(2,3)-R(1,3)*D(2,1)-D(1,1)*R(1,3)-R(1,1)*D(1,3))*tt4+2*(&
&(-R(1,2)*D(3,2))-R(1,2)*D(2,2)-D(1,2)*R(1,2))*tt3+((-R(1,1)*D(3,2&
&))-R(1,2)*D(3,1)-R(1,1)*D(2,2)-R(1,2)*D(2,1)-D(1,1)*R(1,2)-R(1,1)&
&*D(1,2))*tt2+2*((-R(1,1)*D(3,1))-R(1,1)*D(2,1)-D(1,1)*R(1,1))*tt1&
&))
jac(1,2) = volume(1,1)*(lam(1,1)*(R(2,3)*tt9+R(2,2)*tt8+R(2,1)*tt&
&7)*tt22+miu(1,1)*(2*((-R(2,3)*D(3,3))-D(2,3)*R(2,3)-D(1,3)*R(2,3)&
&)*tt6+((-R(2,2)*D(3,3))-R(2,3)*D(3,2)-D(2,2)*R(2,3)-D(1,2)*R(2,3)&
&-R(2,2)*D(2,3)-D(1,3)*R(2,2))*tt5+((-R(2,1)*D(3,3))-R(2,3)*D(3,1)&
&-D(2,1)*R(2,3)-D(1,1)*R(2,3)-R(2,1)*D(2,3)-D(1,3)*R(2,1))*tt4+2*(&
&(-R(2,2)*D(3,2))-D(2,2)*R(2,2)-D(1,2)*R(2,2))*tt3+((-R(2,1)*D(3,2&
&))-R(2,2)*D(3,1)-D(2,1)*R(2,2)-D(1,1)*R(2,2)-R(2,1)*D(2,2)-D(1,2)&
&*R(2,1))*tt2+2*((-R(2,1)*D(3,1))-D(2,1)*R(2,1)-D(1,1)*R(2,1))*tt1&
&))
jac(1,3) = volume(1,1)*(lam(1,1)*(tt9*R(3,3)+tt8*R(3,2)+tt7*R(3,1&
&))*tt22+miu(1,1)*(2*((-D(3,3)*R(3,3))-D(2,3)*R(3,3)-D(1,3)*R(3,3)&
&)*tt6+((-D(3,2)*R(3,3))-D(2,2)*R(3,3)-D(1,2)*R(3,3)-R(3,2)*D(3,3)&
&-D(2,3)*R(3,2)-D(1,3)*R(3,2))*tt5+((-D(3,1)*R(3,3))-D(2,1)*R(3,3)&
&-D(1,1)*R(3,3)-R(3,1)*D(3,3)-D(2,3)*R(3,1)-D(1,3)*R(3,1))*tt4+2*(&
&(-D(3,2)*R(3,2))-D(2,2)*R(3,2)-D(1,2)*R(3,2))*tt3+((-D(3,1)*R(3,2&
&))-D(2,1)*R(3,2)-D(1,1)*R(3,2)-R(3,1)*D(3,2)-D(2,2)*R(3,1)-D(1,2)&
&*R(3,1))*tt2+2*((-D(3,1)*R(3,1))-D(2,1)*R(3,1)-D(1,1)*R(3,1))*tt1&
&))
jac(1,4) = volume(1,1)*(lam(1,1)*(D(1,3)*R(1,3)+D(1,2)*R(1,2)+D(1&
&,1)*R(1,1))*tt22+miu(1,1)*(2*D(1,3)*R(1,3)*tt6+(D(1,2)*R(1,3)+R(1&
&,2)*D(1,3))*tt5+(D(1,1)*R(1,3)+R(1,1)*D(1,3))*tt4+2*D(1,2)*R(1,2)&
&*tt3+(D(1,1)*R(1,2)+R(1,1)*D(1,2))*tt2+2*D(1,1)*R(1,1)*tt1))
jac(1,5) = volume(1,1)*(lam(1,1)*(D(1,3)*R(2,3)+D(1,2)*R(2,2)+D(1&
&,1)*R(2,1))*tt22+miu(1,1)*(2*D(1,3)*R(2,3)*tt6+(D(1,2)*R(2,3)+D(1&
&,3)*R(2,2))*tt5+(D(1,1)*R(2,3)+D(1,3)*R(2,1))*tt4+2*D(1,2)*R(2,2)&
&*tt3+(D(1,1)*R(2,2)+D(1,2)*R(2,1))*tt2+2*D(1,1)*R(2,1)*tt1))
jac(1,6) = volume(1,1)*(lam(1,1)*(D(1,3)*R(3,3)+D(1,2)*R(3,2)+D(1&
&,1)*R(3,1))*tt22+miu(1,1)*(2*D(1,3)*R(3,3)*tt6+(D(1,2)*R(3,3)+D(1&
&,3)*R(3,2))*tt5+(D(1,1)*R(3,3)+D(1,3)*R(3,1))*tt4+2*D(1,2)*R(3,2)&
&*tt3+(D(1,1)*R(3,2)+D(1,2)*R(3,1))*tt2+2*D(1,1)*R(3,1)*tt1))
jac(1,7) = volume(1,1)*(lam(1,1)*(R(1,3)*D(2,3)+R(1,2)*D(2,2)+R(1&
&,1)*D(2,1))*tt22+miu(1,1)*(2*R(1,3)*D(2,3)*tt6+(R(1,2)*D(2,3)+R(1&
&,3)*D(2,2))*tt5+(R(1,1)*D(2,3)+R(1,3)*D(2,1))*tt4+2*R(1,2)*D(2,2)&
&*tt3+(R(1,1)*D(2,2)+R(1,2)*D(2,1))*tt2+2*R(1,1)*D(2,1)*tt1))
jac(1,8) = volume(1,1)*(lam(1,1)*(D(2,3)*R(2,3)+D(2,2)*R(2,2)+D(2&
&,1)*R(2,1))*tt22+miu(1,1)*(2*D(2,3)*R(2,3)*tt6+(D(2,2)*R(2,3)+R(2&
&,2)*D(2,3))*tt5+(D(2,1)*R(2,3)+R(2,1)*D(2,3))*tt4+2*D(2,2)*R(2,2)&
&*tt3+(D(2,1)*R(2,2)+R(2,1)*D(2,2))*tt2+2*D(2,1)*R(2,1)*tt1))
jac(1,9) = volume(1,1)*(lam(1,1)*(D(2,3)*R(3,3)+D(2,2)*R(3,2)+D(2&
&,1)*R(3,1))*tt22+miu(1,1)*(2*D(2,3)*R(3,3)*tt6+(D(2,2)*R(3,3)+D(2&
&,3)*R(3,2))*tt5+(D(2,1)*R(3,3)+D(2,3)*R(3,1))*tt4+2*D(2,2)*R(3,2)&
&*tt3+(D(2,1)*R(3,2)+D(2,2)*R(3,1))*tt2+2*D(2,1)*R(3,1)*tt1))
jac(1,10) = volume(1,1)*(lam(1,1)*(R(1,3)*D(3,3)+R(1,2)*D(3,2)+R(&
&1,1)*D(3,1))*tt22+miu(1,1)*(2*R(1,3)*D(3,3)*tt6+(R(1,2)*D(3,3)+R(&
&1,3)*D(3,2))*tt5+(R(1,1)*D(3,3)+R(1,3)*D(3,1))*tt4+2*R(1,2)*D(3,2&
&)*tt3+(R(1,1)*D(3,2)+R(1,2)*D(3,1))*tt2+2*R(1,1)*D(3,1)*tt1))
jac(1,11) = volume(1,1)*(lam(1,1)*(R(2,3)*D(3,3)+R(2,2)*D(3,2)+R(&
&2,1)*D(3,1))*tt22+miu(1,1)*(2*R(2,3)*D(3,3)*tt6+(R(2,2)*D(3,3)+R(&
&2,3)*D(3,2))*tt5+(R(2,1)*D(3,3)+R(2,3)*D(3,1))*tt4+2*R(2,2)*D(3,2&
&)*tt3+(R(2,1)*D(3,2)+R(2,2)*D(3,1))*tt2+2*R(2,1)*D(3,1)*tt1))
jac(1,12) = volume(1,1)*(lam(1,1)*(D(3,3)*R(3,3)+D(3,2)*R(3,2)+D(&
&3,1)*R(3,1))*tt22+miu(1,1)*(2*D(3,3)*R(3,3)*tt6+(D(3,2)*R(3,3)+R(&
&3,2)*D(3,3))*tt5+(D(3,1)*R(3,3)+R(3,1)*D(3,3))*tt4+2*D(3,2)*R(3,2&
&)*tt3+(D(3,1)*R(3,2)+R(3,1)*D(3,2))*tt2+2*D(3,1)*R(3,1)*tt1))
END 
SUBROUTINE tet_coro_hes(hes, X, D, R, volume, lam, miu) 
IMPLICIT NONE 
REAL(KIND=8) hes(12, 12) 
REAL(KIND=8) X(3, 4) 
REAL(KIND=8) D(3, 3) 
REAL(KIND=8) R(3, 3) 
REAL(KIND=8) volume(1, 1) 
REAL(KIND=8) lam(1, 1) 
REAL(KIND=8) miu(1, 1) 
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
tt1 = (-D(3,1))-D(2,1)-D(1,1)
tt2 = (-D(3,2))-D(2,2)-D(1,2)
tt3 = (-D(3,3))-D(2,3)-D(1,3)
tt4 = R(1,3)*tt3+R(1,2)*tt2+R(1,1)*tt1
tt5 = (-R(1,1)*D(3,1))-R(1,1)*D(2,1)-D(1,1)*R(1,1)
tt6 = (-R(1,1)*D(3,2))-R(1,2)*D(3,1)-R(1,1)*D(2,2)-R(1,2)*D(2,1)-&
&D(1,1)*R(1,2)-R(1,1)*D(1,2)
tt7 = (-R(1,2)*D(3,2))-R(1,2)*D(2,2)-D(1,2)*R(1,2)
tt8 = (-R(1,1)*D(3,3))-R(1,3)*D(3,1)-R(1,1)*D(2,3)-R(1,3)*D(2,1)-&
&D(1,1)*R(1,3)-R(1,1)*D(1,3)
tt9 = (-R(1,2)*D(3,3))-R(1,3)*D(3,2)-R(1,2)*D(2,3)-R(1,3)*D(2,2)-&
&D(1,2)*R(1,3)-R(1,2)*D(1,3)
tt10 = (-R(1,3)*D(3,3))-R(1,3)*D(2,3)-D(1,3)*R(1,3)
tt11 = R(2,3)*tt3+R(2,2)*tt2+R(2,1)*tt1
tt12 = (-R(2,1)*D(3,1))-D(2,1)*R(2,1)-D(1,1)*R(2,1)
tt13 = (-R(2,1)*D(3,2))-R(2,2)*D(3,1)-D(2,1)*R(2,2)-D(1,1)*R(2,2)&
&-R(2,1)*D(2,2)-D(1,2)*R(2,1)
tt14 = (-R(2,2)*D(3,2))-D(2,2)*R(2,2)-D(1,2)*R(2,2)
tt15 = (-R(2,1)*D(3,3))-R(2,3)*D(3,1)-D(2,1)*R(2,3)-D(1,1)*R(2,3)&
&-R(2,1)*D(2,3)-D(1,3)*R(2,1)
tt16 = (-R(2,2)*D(3,3))-R(2,3)*D(3,2)-D(2,2)*R(2,3)-D(1,2)*R(2,3)&
&-R(2,2)*D(2,3)-D(1,3)*R(2,2)
tt17 = (-R(2,3)*D(3,3))-D(2,3)*R(2,3)-D(1,3)*R(2,3)
tt18 = volume(1,1)*(miu(1,1)*(2*tt10*tt17+tt9*tt16+tt8*tt15+2*tt7&
&*tt14+tt6*tt13+2*tt5*tt12)+lam(1,1)*tt4*tt11)
tt19 = tt3*R(3,3)+tt2*R(3,2)+tt1*R(3,1)
tt20 = (-D(3,1)*R(3,1))-D(2,1)*R(3,1)-D(1,1)*R(3,1)
tt21 = (-D(3,1)*R(3,2))-D(2,1)*R(3,2)-D(1,1)*R(3,2)-R(3,1)*D(3,2)&
&-D(2,2)*R(3,1)-D(1,2)*R(3,1)
tt22 = (-D(3,2)*R(3,2))-D(2,2)*R(3,2)-D(1,2)*R(3,2)
tt23 = (-D(3,1)*R(3,3))-D(2,1)*R(3,3)-D(1,1)*R(3,3)-R(3,1)*D(3,3)&
&-D(2,3)*R(3,1)-D(1,3)*R(3,1)
tt24 = (-D(3,2)*R(3,3))-D(2,2)*R(3,3)-D(1,2)*R(3,3)-R(3,2)*D(3,3)&
&-D(2,3)*R(3,2)-D(1,3)*R(3,2)
tt25 = (-D(3,3)*R(3,3))-D(2,3)*R(3,3)-D(1,3)*R(3,3)
tt26 = volume(1,1)*(miu(1,1)*(2*tt10*tt25+tt9*tt24+tt8*tt23+2*tt7&
&*tt22+tt6*tt21+2*tt5*tt20)+lam(1,1)*tt4*tt19)
tt27 = D(1,3)*R(1,3)+D(1,2)*R(1,2)+D(1,1)*R(1,1)
tt28 = D(1,1)*R(1,2)+R(1,1)*D(1,2)
tt29 = D(1,1)*R(1,3)+R(1,1)*D(1,3)
tt30 = D(1,2)*R(1,3)+R(1,2)*D(1,3)
tt31 = volume(1,1)*(miu(1,1)*(2*D(1,3)*R(1,3)*tt10+tt30*tt9+tt29*&
&tt8+2*D(1,2)*R(1,2)*tt7+tt28*tt6+2*D(1,1)*R(1,1)*tt5)+lam(1,1)*tt&
&27*tt4)
tt32 = D(1,3)*R(2,3)+D(1,2)*R(2,2)+D(1,1)*R(2,1)
tt33 = D(1,1)*R(2,2)+D(1,2)*R(2,1)
tt34 = D(1,1)*R(2,3)+D(1,3)*R(2,1)
tt35 = D(1,2)*R(2,3)+D(1,3)*R(2,2)
tt36 = volume(1,1)*(miu(1,1)*(2*D(1,3)*R(2,3)*tt10+tt35*tt9+tt34*&
&tt8+2*D(1,2)*R(2,2)*tt7+tt33*tt6+2*D(1,1)*R(2,1)*tt5)+lam(1,1)*tt&
&32*tt4)
tt37 = D(1,3)*R(3,3)+D(1,2)*R(3,2)+D(1,1)*R(3,1)
tt38 = D(1,1)*R(3,2)+D(1,2)*R(3,1)
tt39 = D(1,1)*R(3,3)+D(1,3)*R(3,1)
tt40 = D(1,2)*R(3,3)+D(1,3)*R(3,2)
tt41 = volume(1,1)*(miu(1,1)*(tt9*tt40+tt8*tt39+2*D(1,3)*tt10*R(3&
&,3)+tt6*tt38+2*D(1,2)*tt7*R(3,2)+2*D(1,1)*tt5*R(3,1))+lam(1,1)*tt&
&4*tt37)
tt42 = R(1,3)*D(2,3)+R(1,2)*D(2,2)+R(1,1)*D(2,1)
tt43 = R(1,1)*D(2,2)+R(1,2)*D(2,1)
tt44 = R(1,1)*D(2,3)+R(1,3)*D(2,1)
tt45 = R(1,2)*D(2,3)+R(1,3)*D(2,2)
tt46 = volume(1,1)*(miu(1,1)*(2*R(1,3)*D(2,3)*tt10+tt45*tt9+tt44*&
&tt8+2*R(1,2)*D(2,2)*tt7+tt43*tt6+2*R(1,1)*D(2,1)*tt5)+lam(1,1)*tt&
&42*tt4)
tt47 = D(2,3)*R(2,3)+D(2,2)*R(2,2)+D(2,1)*R(2,1)
tt48 = D(2,1)*R(2,2)+R(2,1)*D(2,2)
tt49 = D(2,1)*R(2,3)+R(2,1)*D(2,3)
tt50 = D(2,2)*R(2,3)+R(2,2)*D(2,3)
tt51 = volume(1,1)*(miu(1,1)*(2*D(2,3)*R(2,3)*tt10+tt50*tt9+tt49*&
&tt8+2*D(2,2)*R(2,2)*tt7+tt48*tt6+2*D(2,1)*R(2,1)*tt5)+lam(1,1)*tt&
&47*tt4)
tt52 = D(2,3)*R(3,3)+D(2,2)*R(3,2)+D(2,1)*R(3,1)
tt53 = D(2,1)*R(3,2)+D(2,2)*R(3,1)
tt54 = D(2,1)*R(3,3)+D(2,3)*R(3,1)
tt55 = D(2,2)*R(3,3)+D(2,3)*R(3,2)
tt56 = volume(1,1)*(miu(1,1)*(tt9*tt55+tt8*tt54+2*D(2,3)*tt10*R(3&
&,3)+tt6*tt53+2*D(2,2)*tt7*R(3,2)+2*D(2,1)*tt5*R(3,1))+lam(1,1)*tt&
&4*tt52)
tt57 = R(1,3)*D(3,3)+R(1,2)*D(3,2)+R(1,1)*D(3,1)
tt58 = R(1,1)*D(3,2)+R(1,2)*D(3,1)
tt59 = R(1,1)*D(3,3)+R(1,3)*D(3,1)
tt60 = R(1,2)*D(3,3)+R(1,3)*D(3,2)
tt61 = volume(1,1)*(miu(1,1)*(2*R(1,3)*D(3,3)*tt10+tt9*tt60+tt8*t&
&t59+2*R(1,2)*D(3,2)*tt7+tt6*tt58+2*R(1,1)*D(3,1)*tt5)+lam(1,1)*tt&
&4*tt57)
tt62 = R(2,3)*D(3,3)+R(2,2)*D(3,2)+R(2,1)*D(3,1)
tt63 = R(2,1)*D(3,2)+R(2,2)*D(3,1)
tt64 = R(2,1)*D(3,3)+R(2,3)*D(3,1)
tt65 = R(2,2)*D(3,3)+R(2,3)*D(3,2)
tt66 = volume(1,1)*(miu(1,1)*(tt9*tt65+tt8*tt64+2*R(2,3)*D(3,3)*t&
&t10+tt6*tt63+2*R(2,2)*D(3,2)*tt7+2*R(2,1)*D(3,1)*tt5)+lam(1,1)*tt&
&4*tt62)
tt67 = D(3,3)*R(3,3)+D(3,2)*R(3,2)+D(3,1)*R(3,1)
tt68 = D(3,1)*R(3,2)+R(3,1)*D(3,2)
tt69 = D(3,1)*R(3,3)+R(3,1)*D(3,3)
tt70 = D(3,2)*R(3,3)+R(3,2)*D(3,3)
tt71 = volume(1,1)*(miu(1,1)*(tt9*tt70+tt8*tt69+2*D(3,3)*tt10*R(3&
&,3)+tt6*tt68+2*D(3,2)*tt7*R(3,2)+2*D(3,1)*tt5*R(3,1))+lam(1,1)*tt&
&4*tt67)
tt72 = volume(1,1)*(miu(1,1)*(2*tt17*tt25+tt16*tt24+tt15*tt23+2*t&
&t14*tt22+tt13*tt21+2*tt12*tt20)+lam(1,1)*tt11*tt19)
tt73 = volume(1,1)*(miu(1,1)*(2*D(1,3)*R(1,3)*tt17+tt30*tt16+tt29&
&*tt15+2*D(1,2)*R(1,2)*tt14+tt28*tt13+2*D(1,1)*R(1,1)*tt12)+lam(1,&
&1)*tt27*tt11)
tt74 = volume(1,1)*(miu(1,1)*(2*D(1,3)*R(2,3)*tt17+tt35*tt16+tt34&
&*tt15+2*D(1,2)*R(2,2)*tt14+tt33*tt13+2*D(1,1)*R(2,1)*tt12)+lam(1,&
&1)*tt32*tt11)
tt75 = volume(1,1)*(miu(1,1)*(tt16*tt40+tt15*tt39+2*D(1,3)*tt17*R&
&(3,3)+tt13*tt38+2*D(1,2)*tt14*R(3,2)+2*D(1,1)*tt12*R(3,1))+lam(1,&
&1)*tt11*tt37)
tt76 = volume(1,1)*(miu(1,1)*(2*R(1,3)*D(2,3)*tt17+tt45*tt16+tt44&
&*tt15+2*R(1,2)*D(2,2)*tt14+tt43*tt13+2*R(1,1)*D(2,1)*tt12)+lam(1,&
&1)*tt42*tt11)
tt77 = volume(1,1)*(miu(1,1)*(2*D(2,3)*R(2,3)*tt17+tt50*tt16+tt49&
&*tt15+2*D(2,2)*R(2,2)*tt14+tt48*tt13+2*D(2,1)*R(2,1)*tt12)+lam(1,&
&1)*tt47*tt11)
tt78 = volume(1,1)*(miu(1,1)*(tt16*tt55+tt15*tt54+2*D(2,3)*tt17*R&
&(3,3)+tt13*tt53+2*D(2,2)*tt14*R(3,2)+2*D(2,1)*tt12*R(3,1))+lam(1,&
&1)*tt11*tt52)
tt79 = volume(1,1)*(miu(1,1)*(2*R(1,3)*D(3,3)*tt17+tt60*tt16+tt59&
&*tt15+2*R(1,2)*D(3,2)*tt14+tt58*tt13+2*R(1,1)*D(3,1)*tt12)+lam(1,&
&1)*tt11*tt57)
tt80 = volume(1,1)*(miu(1,1)*(2*R(2,3)*D(3,3)*tt17+tt16*tt65+tt15&
&*tt64+2*R(2,2)*D(3,2)*tt14+tt13*tt63+2*R(2,1)*D(3,1)*tt12)+lam(1,&
&1)*tt11*tt62)
tt81 = volume(1,1)*(miu(1,1)*(tt16*tt70+tt15*tt69+2*D(3,3)*tt17*R&
&(3,3)+tt13*tt68+2*D(3,2)*tt14*R(3,2)+2*D(3,1)*tt12*R(3,1))+lam(1,&
&1)*tt11*tt67)
tt82 = volume(1,1)*(miu(1,1)*(2*D(1,3)*R(1,3)*tt25+tt30*tt24+tt29&
&*tt23+2*D(1,2)*R(1,2)*tt22+tt28*tt21+2*D(1,1)*R(1,1)*tt20)+lam(1,&
&1)*tt27*tt19)
tt83 = volume(1,1)*(miu(1,1)*(2*D(1,3)*R(2,3)*tt25+tt35*tt24+tt34&
&*tt23+2*D(1,2)*R(2,2)*tt22+tt33*tt21+2*D(1,1)*R(2,1)*tt20)+lam(1,&
&1)*tt32*tt19)
tt84 = volume(1,1)*(miu(1,1)*(2*D(1,3)*R(3,3)*tt25+tt40*tt24+tt39&
&*tt23+2*D(1,2)*R(3,2)*tt22+tt38*tt21+2*D(1,1)*R(3,1)*tt20)+lam(1,&
&1)*tt37*tt19)
tt85 = volume(1,1)*(miu(1,1)*(2*R(1,3)*D(2,3)*tt25+tt45*tt24+tt44&
&*tt23+2*R(1,2)*D(2,2)*tt22+tt43*tt21+2*R(1,1)*D(2,1)*tt20)+lam(1,&
&1)*tt42*tt19)
tt86 = volume(1,1)*(miu(1,1)*(2*D(2,3)*R(2,3)*tt25+tt50*tt24+tt49&
&*tt23+2*D(2,2)*R(2,2)*tt22+tt48*tt21+2*D(2,1)*R(2,1)*tt20)+lam(1,&
&1)*tt47*tt19)
tt87 = volume(1,1)*(miu(1,1)*(2*D(2,3)*R(3,3)*tt25+tt55*tt24+tt54&
&*tt23+2*D(2,2)*R(3,2)*tt22+tt53*tt21+2*D(2,1)*R(3,1)*tt20)+lam(1,&
&1)*tt52*tt19)
tt88 = volume(1,1)*(miu(1,1)*(2*R(1,3)*D(3,3)*tt25+tt60*tt24+tt59&
&*tt23+2*R(1,2)*D(3,2)*tt22+tt58*tt21+2*R(1,1)*D(3,1)*tt20)+lam(1,&
&1)*tt57*tt19)
tt89 = volume(1,1)*(miu(1,1)*(2*R(2,3)*D(3,3)*tt25+tt65*tt24+tt64&
&*tt23+2*R(2,2)*D(3,2)*tt22+tt63*tt21+2*R(2,1)*D(3,1)*tt20)+lam(1,&
&1)*tt62*tt19)
tt90 = volume(1,1)*(miu(1,1)*(2*D(3,3)*R(3,3)*tt25+tt24*tt70+tt23&
&*tt69+2*D(3,2)*R(3,2)*tt22+tt21*tt68+2*D(3,1)*R(3,1)*tt20)+lam(1,&
&1)*tt19*tt67)
tt91 = D(1,1)**2
tt92 = R(1,1)**2
tt93 = D(1,2)**2
tt94 = R(1,2)**2
tt95 = D(1,3)**2
tt96 = R(1,3)**2
tt97 = volume(1,1)*(miu(1,1)*(tt30*tt35+tt29*tt34+2*tt95*R(1,3)*R&
&(2,3)+tt28*tt33+2*tt93*R(1,2)*R(2,2)+2*tt91*R(1,1)*R(2,1))+lam(1,&
&1)*tt27*tt32)
tt98 = volume(1,1)*(miu(1,1)*(tt30*tt40+tt29*tt39+2*tt95*R(1,3)*R&
&(3,3)+tt28*tt38+2*tt93*R(1,2)*R(3,2)+2*tt91*R(1,1)*R(3,1))+lam(1,&
&1)*tt27*tt37)
tt99 = volume(1,1)*(miu(1,1)*(tt30*tt45+tt29*tt44+2*D(1,3)*tt96*D&
&(2,3)+tt28*tt43+2*D(1,2)*tt94*D(2,2)+2*D(1,1)*tt92*D(2,1))+lam(1,&
&1)*tt27*tt42)
tt100 = 2*D(1,1)*R(1,1)*D(2,1)*R(2,1)
tt101 = 2*D(1,2)*R(1,2)*D(2,2)*R(2,2)
tt102 = 2*D(1,3)*R(1,3)*D(2,3)*R(2,3)
tt103 = volume(1,1)*(miu(1,1)*(tt30*tt50+tt29*tt49+tt102+tt28*tt4&
&8+tt101+tt100)+lam(1,1)*tt27*tt47)
tt104 = 2*D(1,1)*R(1,1)*D(2,1)*R(3,1)
tt105 = 2*D(1,2)*R(1,2)*D(2,2)*R(3,2)
tt106 = 2*D(1,3)*R(1,3)*D(2,3)*R(3,3)
tt107 = volume(1,1)*(miu(1,1)*(tt30*tt55+tt29*tt54+tt106+tt28*tt5&
&3+tt105+tt104)+lam(1,1)*tt27*tt52)
tt108 = volume(1,1)*(miu(1,1)*(tt30*tt60+tt29*tt59+2*D(1,3)*tt96*&
&D(3,3)+tt28*tt58+2*D(1,2)*tt94*D(3,2)+2*D(1,1)*tt92*D(3,1))+lam(1&
&,1)*tt27*tt57)
tt109 = 2*D(1,1)*R(1,1)*R(2,1)*D(3,1)
tt110 = 2*D(1,2)*R(1,2)*R(2,2)*D(3,2)
tt111 = 2*D(1,3)*R(1,3)*R(2,3)*D(3,3)
tt112 = volume(1,1)*(miu(1,1)*(tt30*tt65+tt29*tt64+tt111+tt28*tt6&
&3+tt110+tt109)+lam(1,1)*tt27*tt62)
tt113 = 2*D(1,1)*R(1,1)*D(3,1)*R(3,1)
tt114 = 2*D(1,2)*R(1,2)*D(3,2)*R(3,2)
tt115 = 2*D(1,3)*R(1,3)*D(3,3)*R(3,3)
tt116 = volume(1,1)*(miu(1,1)*(tt30*tt70+tt29*tt69+tt115+tt28*tt6&
&8+tt114+tt113)+lam(1,1)*tt27*tt67)
tt117 = R(2,1)**2
tt118 = R(2,2)**2
tt119 = R(2,3)**2
tt120 = volume(1,1)*(miu(1,1)*(tt35*tt40+tt34*tt39+2*tt95*R(2,3)*&
&R(3,3)+tt33*tt38+2*tt93*R(2,2)*R(3,2)+2*tt91*R(2,1)*R(3,1))+lam(1&
&,1)*tt32*tt37)
tt121 = volume(1,1)*(miu(1,1)*(tt45*tt35+tt44*tt34+tt102+tt43*tt3&
&3+tt101+tt100)+lam(1,1)*tt42*tt32)
tt122 = volume(1,1)*(miu(1,1)*(2*D(1,3)*D(2,3)*tt119+tt35*tt50+tt&
&34*tt49+2*D(1,2)*D(2,2)*tt118+tt33*tt48+2*D(1,1)*D(2,1)*tt117)+la&
&m(1,1)*tt32*tt47)
tt123 = 2*D(1,1)*D(2,1)*R(2,1)*R(3,1)
tt124 = 2*D(1,2)*D(2,2)*R(2,2)*R(3,2)
tt125 = 2*D(1,3)*D(2,3)*R(2,3)*R(3,3)
tt126 = volume(1,1)*(miu(1,1)*(tt35*tt55+tt34*tt54+tt125+tt33*tt5&
&3+tt124+tt123)+lam(1,1)*tt32*tt52)
tt127 = volume(1,1)*(miu(1,1)*(tt35*tt60+tt34*tt59+tt111+tt33*tt5&
&8+tt110+tt109)+lam(1,1)*tt32*tt57)
tt128 = volume(1,1)*(miu(1,1)*(tt35*tt65+tt34*tt64+2*D(1,3)*tt119&
&*D(3,3)+tt33*tt63+2*D(1,2)*tt118*D(3,2)+2*D(1,1)*tt117*D(3,1))+la&
&m(1,1)*tt32*tt62)
tt129 = 2*D(1,1)*R(2,1)*D(3,1)*R(3,1)
tt130 = 2*D(1,2)*R(2,2)*D(3,2)*R(3,2)
tt131 = 2*D(1,3)*R(2,3)*D(3,3)*R(3,3)
tt132 = volume(1,1)*(miu(1,1)*(tt35*tt70+tt34*tt69+tt131+tt33*tt6&
&8+tt130+tt129)+lam(1,1)*tt32*tt67)
tt133 = R(3,1)**2
tt134 = R(3,2)**2
tt135 = R(3,3)**2
tt136 = volume(1,1)*(miu(1,1)*(tt45*tt40+tt44*tt39+tt106+tt43*tt3&
&8+tt105+tt104)+lam(1,1)*tt42*tt37)
tt137 = volume(1,1)*(miu(1,1)*(tt50*tt40+tt49*tt39+tt125+tt48*tt3&
&8+tt124+tt123)+lam(1,1)*tt47*tt37)
tt138 = volume(1,1)*(miu(1,1)*(2*D(1,3)*D(2,3)*tt135+tt40*tt55+tt&
&39*tt54+2*D(1,2)*D(2,2)*tt134+tt38*tt53+2*D(1,1)*D(2,1)*tt133)+la&
&m(1,1)*tt37*tt52)
tt139 = volume(1,1)*(miu(1,1)*(tt60*tt40+tt59*tt39+tt115+tt58*tt3&
&8+tt114+tt113)+lam(1,1)*tt57*tt37)
tt140 = volume(1,1)*(miu(1,1)*(tt65*tt40+tt64*tt39+tt131+tt63*tt3&
&8+tt130+tt129)+lam(1,1)*tt62*tt37)
tt141 = volume(1,1)*(miu(1,1)*(2*D(1,3)*D(3,3)*tt135+tt40*tt70+tt&
&39*tt69+2*D(1,2)*D(3,2)*tt134+tt38*tt68+2*D(1,1)*D(3,1)*tt133)+la&
&m(1,1)*tt37*tt67)
tt142 = D(2,1)**2
tt143 = D(2,2)**2
tt144 = D(2,3)**2
tt145 = volume(1,1)*(miu(1,1)*(tt45*tt50+tt44*tt49+2*R(1,3)*tt144&
&*R(2,3)+tt43*tt48+2*R(1,2)*tt143*R(2,2)+2*R(1,1)*tt142*R(2,1))+la&
&m(1,1)*tt42*tt47)
tt146 = volume(1,1)*(miu(1,1)*(tt45*tt55+tt44*tt54+2*R(1,3)*tt144&
&*R(3,3)+tt43*tt53+2*R(1,2)*tt143*R(3,2)+2*R(1,1)*tt142*R(3,1))+la&
&m(1,1)*tt42*tt52)
tt147 = volume(1,1)*(miu(1,1)*(tt45*tt60+tt44*tt59+2*tt96*D(2,3)*&
&D(3,3)+tt43*tt58+2*tt94*D(2,2)*D(3,2)+2*tt92*D(2,1)*D(3,1))+lam(1&
&,1)*tt42*tt57)
tt148 = 2*R(1,1)*D(2,1)*R(2,1)*D(3,1)
tt149 = 2*R(1,2)*D(2,2)*R(2,2)*D(3,2)
tt150 = 2*R(1,3)*D(2,3)*R(2,3)*D(3,3)
tt151 = volume(1,1)*(miu(1,1)*(tt45*tt65+tt44*tt64+tt150+tt43*tt6&
&3+tt149+tt148)+lam(1,1)*tt42*tt62)
tt152 = 2*R(1,1)*D(2,1)*D(3,1)*R(3,1)
tt153 = 2*R(1,2)*D(2,2)*D(3,2)*R(3,2)
tt154 = 2*R(1,3)*D(2,3)*D(3,3)*R(3,3)
tt155 = volume(1,1)*(miu(1,1)*(tt45*tt70+tt44*tt69+tt154+tt43*tt6&
&8+tt153+tt152)+lam(1,1)*tt42*tt67)
tt156 = volume(1,1)*(miu(1,1)*(tt50*tt55+tt49*tt54+2*tt144*R(2,3)&
&*R(3,3)+tt48*tt53+2*tt143*R(2,2)*R(3,2)+2*tt142*R(2,1)*R(3,1))+la&
&m(1,1)*tt47*tt52)
tt157 = volume(1,1)*(miu(1,1)*(tt50*tt60+tt49*tt59+tt150+tt48*tt5&
&8+tt149+tt148)+lam(1,1)*tt47*tt57)
tt158 = volume(1,1)*(miu(1,1)*(tt50*tt65+tt49*tt64+2*D(2,3)*tt119&
&*D(3,3)+tt48*tt63+2*D(2,2)*tt118*D(3,2)+2*D(2,1)*tt117*D(3,1))+la&
&m(1,1)*tt47*tt62)
tt159 = 2*D(2,1)*R(2,1)*D(3,1)*R(3,1)
tt160 = 2*D(2,2)*R(2,2)*D(3,2)*R(3,2)
tt161 = 2*D(2,3)*R(2,3)*D(3,3)*R(3,3)
tt162 = volume(1,1)*(miu(1,1)*(tt50*tt70+tt49*tt69+tt161+tt48*tt6&
&8+tt160+tt159)+lam(1,1)*tt47*tt67)
tt163 = volume(1,1)*(miu(1,1)*(tt60*tt55+tt59*tt54+tt154+tt58*tt5&
&3+tt153+tt152)+lam(1,1)*tt57*tt52)
tt164 = volume(1,1)*(miu(1,1)*(tt65*tt55+tt64*tt54+tt161+tt63*tt5&
&3+tt160+tt159)+lam(1,1)*tt62*tt52)
tt165 = volume(1,1)*(miu(1,1)*(2*D(2,3)*D(3,3)*tt135+tt55*tt70+tt&
&54*tt69+2*D(2,2)*D(3,2)*tt134+tt53*tt68+2*D(2,1)*D(3,1)*tt133)+la&
&m(1,1)*tt52*tt67)
tt166 = D(3,1)**2
tt167 = D(3,2)**2
tt168 = D(3,3)**2
tt169 = volume(1,1)*(miu(1,1)*(2*R(1,3)*R(2,3)*tt168+tt60*tt65+tt&
&59*tt64+2*R(1,2)*R(2,2)*tt167+tt58*tt63+2*R(1,1)*R(2,1)*tt166)+la&
&m(1,1)*tt57*tt62)
tt170 = volume(1,1)*(miu(1,1)*(tt60*tt70+tt59*tt69+2*R(1,3)*tt168&
&*R(3,3)+tt58*tt68+2*R(1,2)*tt167*R(3,2)+2*R(1,1)*tt166*R(3,1))+la&
&m(1,1)*tt57*tt67)
tt171 = volume(1,1)*(miu(1,1)*(tt65*tt70+tt64*tt69+2*R(2,3)*tt168&
&*R(3,3)+tt63*tt68+2*R(2,2)*tt167*R(3,2)+2*R(2,1)*tt166*R(3,1))+la&
&m(1,1)*tt62*tt67)
hes(1,1) = volume(1,1)*(miu(1,1)*(2*tt10**2+tt9**2+tt8**2+2*tt7**&
&2+tt6**2+2*tt5**2)+lam(1,1)*tt4**2)
hes(1,2) = tt18
hes(1,3) = tt26
hes(1,4) = tt31
hes(1,5) = tt36
hes(1,6) = tt41
hes(1,7) = tt46
hes(1,8) = tt51
hes(1,9) = tt56
hes(1,10) = tt61
hes(1,11) = tt66
hes(1,12) = tt71
hes(2,1) = tt18
hes(2,2) = volume(1,1)*(miu(1,1)*(2*tt17**2+tt16**2+tt15**2+2*tt1&
&4**2+tt13**2+2*tt12**2)+lam(1,1)*tt11**2)
hes(2,3) = tt72
hes(2,4) = tt73
hes(2,5) = tt74
hes(2,6) = tt75
hes(2,7) = tt76
hes(2,8) = tt77
hes(2,9) = tt78
hes(2,10) = tt79
hes(2,11) = tt80
hes(2,12) = tt81
hes(3,1) = tt26
hes(3,2) = tt72
hes(3,3) = volume(1,1)*(miu(1,1)*(2*tt25**2+tt24**2+tt23**2+2*tt2&
&2**2+tt21**2+2*tt20**2)+lam(1,1)*tt19**2)
hes(3,4) = tt82
hes(3,5) = tt83
hes(3,6) = tt84
hes(3,7) = tt85
hes(3,8) = tt86
hes(3,9) = tt87
hes(3,10) = tt88
hes(3,11) = tt89
hes(3,12) = tt90
hes(4,1) = tt31
hes(4,2) = tt73
hes(4,3) = tt82
hes(4,4) = volume(1,1)*(lam(1,1)*tt27**2+miu(1,1)*(tt30**2+tt29**&
&2+2*tt95*tt96+tt28**2+2*tt93*tt94+2*tt91*tt92))
hes(4,5) = tt97
hes(4,6) = tt98
hes(4,7) = tt99
hes(4,8) = tt103
hes(4,9) = tt107
hes(4,10) = tt108
hes(4,11) = tt112
hes(4,12) = tt116
hes(5,1) = tt36
hes(5,2) = tt74
hes(5,3) = tt83
hes(5,4) = tt97
hes(5,5) = volume(1,1)*(lam(1,1)*tt32**2+miu(1,1)*(tt35**2+tt34**&
&2+2*tt95*tt119+tt33**2+2*tt93*tt118+2*tt91*tt117))
hes(5,6) = tt120
hes(5,7) = tt121
hes(5,8) = tt122
hes(5,9) = tt126
hes(5,10) = tt127
hes(5,11) = tt128
hes(5,12) = tt132
hes(6,1) = tt41
hes(6,2) = tt75
hes(6,3) = tt84
hes(6,4) = tt98
hes(6,5) = tt120
hes(6,6) = volume(1,1)*(lam(1,1)*tt37**2+miu(1,1)*(tt40**2+tt39**&
&2+2*tt95*tt135+tt38**2+2*tt93*tt134+2*tt91*tt133))
hes(6,7) = tt136
hes(6,8) = tt137
hes(6,9) = tt138
hes(6,10) = tt139
hes(6,11) = tt140
hes(6,12) = tt141
hes(7,1) = tt46
hes(7,2) = tt76
hes(7,3) = tt85
hes(7,4) = tt99
hes(7,5) = tt121
hes(7,6) = tt136
hes(7,7) = volume(1,1)*(lam(1,1)*tt42**2+miu(1,1)*(tt45**2+tt44**&
&2+2*tt96*tt144+tt43**2+2*tt94*tt143+2*tt92*tt142))
hes(7,8) = tt145
hes(7,9) = tt146
hes(7,10) = tt147
hes(7,11) = tt151
hes(7,12) = tt155
hes(8,1) = tt51
hes(8,2) = tt77
hes(8,3) = tt86
hes(8,4) = tt103
hes(8,5) = tt122
hes(8,6) = tt137
hes(8,7) = tt145
hes(8,8) = volume(1,1)*(lam(1,1)*tt47**2+miu(1,1)*(tt50**2+tt49**&
&2+2*tt144*tt119+tt48**2+2*tt143*tt118+2*tt142*tt117))
hes(8,9) = tt156
hes(8,10) = tt157
hes(8,11) = tt158
hes(8,12) = tt162
hes(9,1) = tt56
hes(9,2) = tt78
hes(9,3) = tt87
hes(9,4) = tt107
hes(9,5) = tt126
hes(9,6) = tt138
hes(9,7) = tt146
hes(9,8) = tt156
hes(9,9) = volume(1,1)*(lam(1,1)*tt52**2+miu(1,1)*(tt55**2+tt54**&
&2+2*tt144*tt135+tt53**2+2*tt143*tt134+2*tt142*tt133))
hes(9,10) = tt163
hes(9,11) = tt164
hes(9,12) = tt165
hes(10,1) = tt61
hes(10,2) = tt79
hes(10,3) = tt88
hes(10,4) = tt108
hes(10,5) = tt127
hes(10,6) = tt139
hes(10,7) = tt147
hes(10,8) = tt157
hes(10,9) = tt163
hes(10,10) = volume(1,1)*(lam(1,1)*tt57**2+miu(1,1)*(tt60**2+tt59&
&**2+2*tt96*tt168+tt58**2+2*tt94*tt167+2*tt92*tt166))
hes(10,11) = tt169
hes(10,12) = tt170
hes(11,1) = tt66
hes(11,2) = tt80
hes(11,3) = tt89
hes(11,4) = tt112
hes(11,5) = tt128
hes(11,6) = tt140
hes(11,7) = tt151
hes(11,8) = tt158
hes(11,9) = tt164
hes(11,10) = tt169
hes(11,11) = volume(1,1)*(lam(1,1)*tt62**2+miu(1,1)*(tt65**2+tt64&
&**2+2*tt119*tt168+tt63**2+2*tt118*tt167+2*tt117*tt166))
hes(11,12) = tt171
hes(12,1) = tt71
hes(12,2) = tt81
hes(12,3) = tt90
hes(12,4) = tt116
hes(12,5) = tt132
hes(12,6) = tt141
hes(12,7) = tt155
hes(12,8) = tt162
hes(12,9) = tt165
hes(12,10) = tt170
hes(12,11) = tt171
hes(12,12) = volume(1,1)*(lam(1,1)*tt67**2+miu(1,1)*(tt70**2+tt69&
&**2+2*tt168*tt135+tt68**2+2*tt167*tt134+2*tt166*tt133))
END 
SUBROUTINE tet_neohookean(val, X, D, volume, lam, miu) 
IMPLICIT NONE 
REAL(KIND=8) val(1, 1) 
REAL(KIND=8) X(3, 4) 
REAL(KIND=8) D(3, 3) 
REAL(KIND=8) volume(1, 1) 
REAL(KIND=8) lam(1, 1) 
REAL(KIND=8) miu(1, 1) 
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
tt1 = -X(1,1)
tt2 = X(1,2)+tt1
tt3 = X(1,3)+tt1
tt4 = X(1,4)+tt1
tt5 = tt4*D(3,1)+tt3*D(2,1)+D(1,1)*tt2
tt6 = tt5**2
tt7 = -X(2,1)
tt8 = X(2,2)+tt7
tt9 = X(2,3)+tt7
tt10 = X(2,4)+tt7
tt11 = tt10*D(3,1)+D(2,1)*tt9+D(1,1)*tt8
tt12 = tt11**2
tt13 = tt4*D(3,2)+tt3*D(2,2)+D(1,2)*tt2
tt14 = tt13**2
tt15 = tt10*D(3,2)+D(2,2)*tt9+D(1,2)*tt8
tt16 = tt15**2
tt17 = tt4*D(3,3)+tt3*D(2,3)+tt2*D(1,3)
tt18 = tt17**2
tt19 = tt10*D(3,3)+D(2,3)*tt9+D(1,3)*tt8
tt20 = tt19**2
tt21 = -X(3,1)
tt22 = X(3,2)+tt21
tt23 = X(3,3)+tt21
tt24 = X(3,4)+tt21
tt25 = D(3,1)*tt24+D(2,1)*tt23+D(1,1)*tt22
tt26 = tt25**2
tt27 = D(3,2)*tt24+D(2,2)*tt23+D(1,2)*tt22
tt28 = tt27**2
tt29 = D(3,3)*tt24+D(2,3)*tt23+D(1,3)*tt22
tt30 = tt29**2
tt31 = tt25*tt29+tt11*tt19+tt5*tt17
tt32 = tt25*tt27+tt11*tt15+tt5*tt13
tt33 = tt27*tt29+tt15*tt19+tt13*tt17
tt34 = tt28+tt16+tt14
tt35 = tt30+tt20+tt18
tt36 = log((tt26+tt12+tt6)*(tt34*tt35-tt33**2)-tt32*(tt32*tt35-tt&
&31*tt33)+tt31*(tt32*tt33-tt31*tt34))
val(1,1) = volume(1,1)*(1.25E-1*lam(1,1)*tt36**2+5.0E-1*miu(1,1)*&
&((-tt36)+tt30+tt28+tt26+tt20+tt18+tt16+tt14+tt12+tt6-3))
END 
SUBROUTINE tet_neohookean_jac(jac, X, D, volume, lam, miu) 
IMPLICIT NONE 
REAL(KIND=8) jac(1, 12) 
REAL(KIND=8) X(3, 4) 
REAL(KIND=8) D(3, 3) 
REAL(KIND=8) volume(1, 1) 
REAL(KIND=8) lam(1, 1) 
REAL(KIND=8) miu(1, 1) 
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
tt1 = (-D(3,1))-D(2,1)-D(1,1)
tt2 = -X(1,1)
tt3 = X(1,2)+tt2
tt4 = X(1,3)+tt2
tt5 = X(1,4)+tt2
tt6 = tt5*D(3,1)+tt4*D(2,1)+D(1,1)*tt3
tt7 = (-D(3,2))-D(2,2)-D(1,2)
tt8 = tt5*D(3,2)+tt4*D(2,2)+D(1,2)*tt3
tt9 = (-D(3,3))-D(2,3)-D(1,3)
tt10 = tt5*D(3,3)+tt4*D(2,3)+tt3*D(1,3)
tt11 = -X(2,1)
tt12 = X(2,2)+tt11
tt13 = X(2,3)+tt11
tt14 = X(2,4)+tt11
tt15 = tt14*D(3,1)+D(2,1)*tt13+D(1,1)*tt12
tt16 = tt14*D(3,3)+D(2,3)*tt13+D(1,3)*tt12
tt17 = -X(3,1)
tt18 = X(3,2)+tt17
tt19 = X(3,3)+tt17
tt20 = X(3,4)+tt17
tt21 = D(3,1)*tt20+D(2,1)*tt19+D(1,1)*tt18
tt22 = D(3,3)*tt20+D(2,3)*tt19+D(1,3)*tt18
tt23 = tt21*tt22+tt15*tt16+tt6*tt10
tt24 = tt7*tt10+tt8*tt9
tt25 = tt14*D(3,2)+D(2,2)*tt13+D(1,2)*tt12
tt26 = D(3,2)*tt20+D(2,2)*tt19+D(1,2)*tt18
tt27 = tt21*tt26+tt15*tt25+tt6*tt8
tt28 = tt1*tt8+tt6*tt7
tt29 = tt26*tt22+tt25*tt16+tt8*tt10
tt30 = tt1*tt10+tt6*tt9
tt31 = tt26**2+tt25**2+tt8**2
tt32 = tt27*tt29-tt23*tt31
tt33 = tt21**2+tt15**2+tt6**2
tt34 = tt22**2+tt16**2+tt10**2
tt35 = tt27*tt34-tt23*tt29
tt36 = tt31*tt34-tt29**2
tt37 = 2*tt1*tt6*tt36-tt28*tt35-tt27*(tt28*tt34-tt30*tt29-tt24*tt&
&23+2*tt9*tt10*tt27)+tt33*(2*tt7*tt8*tt34+2*tt9*tt10*tt31-2*tt24*t&
&t29)+tt30*tt32+tt23*((-tt30*tt31)+tt28*tt29-2*tt7*tt8*tt23+tt24*t&
&t27)
tt38 = tt33*tt36-tt27*tt35+tt23*tt32
tt39 = 1/tt38
tt40 = log(tt38)
tt41 = tt7*tt16+tt25*tt9
tt42 = tt1*tt25+tt15*tt7
tt43 = tt1*tt16+tt15*tt9
tt44 = 2*tt1*tt15*tt36-tt42*tt35-tt27*(tt42*tt34-tt43*tt29-tt41*t&
&t23+2*tt9*tt16*tt27)+tt33*(2*tt7*tt25*tt34+2*tt9*tt16*tt31-2*tt41&
&*tt29)+tt43*tt32+tt23*((-tt43*tt31)+tt42*tt29-2*tt7*tt25*tt23+tt4&
&1*tt27)
tt45 = tt7*tt22+tt9*tt26
tt46 = tt1*tt26+tt7*tt21
tt47 = tt1*tt22+tt9*tt21
tt48 = 2*tt1*tt21*tt36-tt46*tt35-tt27*(tt46*tt34-tt47*tt29-tt45*t&
&t23+2*tt9*tt22*tt27)+tt33*(2*tt7*tt26*tt34+2*tt9*tt22*tt31-2*tt45&
&*tt29)+tt47*tt32+tt23*((-tt47*tt31)+tt46*tt29-2*tt7*tt26*tt23+tt2&
&7*tt45)
tt49 = D(1,2)*tt10+D(1,3)*tt8
tt50 = D(1,1)*tt8+D(1,2)*tt6
tt51 = D(1,1)*tt10+D(1,3)*tt6
tt52 = 2*D(1,1)*tt6*tt36-tt50*tt35-tt27*(tt50*tt34-tt51*tt29-tt49&
&*tt23+2*D(1,3)*tt10*tt27)+tt33*(2*D(1,2)*tt8*tt34+2*D(1,3)*tt10*t&
&t31-2*tt49*tt29)+tt51*tt32+tt23*((-tt51*tt31)+tt50*tt29-2*D(1,2)*&
&tt8*tt23+tt49*tt27)
tt53 = D(1,2)*tt16+D(1,3)*tt25
tt54 = D(1,1)*tt25+D(1,2)*tt15
tt55 = D(1,1)*tt16+D(1,3)*tt15
tt56 = 2*D(1,1)*tt15*tt36-tt54*tt35-tt27*(tt54*tt34-tt55*tt29-tt5&
&3*tt23+2*D(1,3)*tt16*tt27)+tt33*(2*D(1,2)*tt25*tt34+2*D(1,3)*tt16&
&*tt31-2*tt53*tt29)+tt55*tt32+tt23*((-tt55*tt31)+tt54*tt29-2*D(1,2&
&)*tt25*tt23+tt53*tt27)
tt57 = D(1,2)*tt22+D(1,3)*tt26
tt58 = D(1,1)*tt26+D(1,2)*tt21
tt59 = D(1,1)*tt22+D(1,3)*tt21
tt60 = 2*D(1,1)*tt21*tt36-tt58*tt35-tt27*(tt58*tt34-tt59*tt29-tt5&
&7*tt23+2*D(1,3)*tt22*tt27)+tt33*(2*D(1,2)*tt26*tt34+2*D(1,3)*tt22&
&*tt31-2*tt57*tt29)+tt59*tt32+tt23*((-tt59*tt31)+tt58*tt29-2*D(1,2&
&)*tt26*tt23+tt27*tt57)
tt61 = D(2,2)*tt10+D(2,3)*tt8
tt62 = D(2,1)*tt8+D(2,2)*tt6
tt63 = D(2,1)*tt10+D(2,3)*tt6
tt64 = 2*D(2,1)*tt6*tt36-tt62*tt35-tt27*(tt62*tt34-tt63*tt29-tt61&
&*tt23+2*D(2,3)*tt10*tt27)+tt33*(2*D(2,2)*tt8*tt34+2*D(2,3)*tt10*t&
&t31-2*tt61*tt29)+tt63*tt32+tt23*((-tt63*tt31)+tt62*tt29-2*D(2,2)*&
&tt8*tt23+tt61*tt27)
tt65 = D(2,2)*tt16+D(2,3)*tt25
tt66 = D(2,1)*tt25+D(2,2)*tt15
tt67 = D(2,1)*tt16+D(2,3)*tt15
tt68 = 2*D(2,1)*tt15*tt36-tt66*tt35-tt27*(tt66*tt34-tt67*tt29-tt6&
&5*tt23+2*D(2,3)*tt16*tt27)+tt33*(2*D(2,2)*tt25*tt34+2*D(2,3)*tt16&
&*tt31-2*tt65*tt29)+tt67*tt32+tt23*((-tt67*tt31)+tt66*tt29-2*D(2,2&
&)*tt25*tt23+tt65*tt27)
tt69 = D(2,2)*tt22+D(2,3)*tt26
tt70 = D(2,1)*tt26+D(2,2)*tt21
tt71 = D(2,1)*tt22+D(2,3)*tt21
tt72 = 2*D(2,1)*tt21*tt36-tt70*tt35-tt27*(tt70*tt34-tt71*tt29-tt6&
&9*tt23+2*D(2,3)*tt22*tt27)+tt33*(2*D(2,2)*tt26*tt34+2*D(2,3)*tt22&
&*tt31-2*tt69*tt29)+tt71*tt32+tt23*((-tt71*tt31)+tt70*tt29-2*D(2,2&
&)*tt26*tt23+tt27*tt69)
tt73 = D(3,2)*tt10+tt8*D(3,3)
tt74 = D(3,1)*tt8+tt6*D(3,2)
tt75 = D(3,1)*tt10+tt6*D(3,3)
tt76 = 2*D(3,1)*tt6*tt36-tt74*tt35-tt27*(tt74*tt34-tt75*tt29-tt73&
&*tt23+2*D(3,3)*tt10*tt27)+tt33*(2*D(3,2)*tt8*tt34+2*D(3,3)*tt10*t&
&t31-2*tt73*tt29)+tt75*tt32+tt23*((-tt75*tt31)+tt74*tt29-2*D(3,2)*&
&tt8*tt23+tt73*tt27)
tt77 = D(3,2)*tt16+tt25*D(3,3)
tt78 = D(3,1)*tt25+tt15*D(3,2)
tt79 = D(3,1)*tt16+tt15*D(3,3)
tt80 = 2*D(3,1)*tt15*tt36-tt78*tt35-tt27*(tt78*tt34-tt79*tt29-tt7&
&7*tt23+2*D(3,3)*tt16*tt27)+tt33*(2*D(3,2)*tt25*tt34+2*D(3,3)*tt16&
&*tt31-2*tt77*tt29)+tt79*tt32+tt23*((-tt79*tt31)+tt78*tt29-2*D(3,2&
&)*tt25*tt23+tt77*tt27)
tt81 = D(3,2)*tt22+D(3,3)*tt26
tt82 = D(3,1)*tt26+D(3,2)*tt21
tt83 = D(3,1)*tt22+D(3,3)*tt21
tt84 = 2*D(3,1)*tt21*tt36-tt82*tt35-tt27*(tt82*tt34-tt83*tt29-tt8&
&1*tt23+2*D(3,3)*tt22*tt27)+tt33*(2*D(3,2)*tt26*tt34+2*D(3,3)*tt22&
&*tt31-2*tt81*tt29)+tt83*tt32+tt23*((-tt83*tt31)+tt82*tt29-2*D(3,2&
&)*tt26*tt23+tt27*tt81)
jac(1,1) = volume(1,1)*(2.5E-1*lam(1,1)*tt37*tt39*tt40+5.0E-1*miu&
&(1,1)*((-tt37*tt39)+2*tt9*tt10+2*tt7*tt8+2*tt1*tt6))
jac(1,2) = volume(1,1)*(2.5E-1*lam(1,1)*tt44*tt39*tt40+5.0E-1*miu&
&(1,1)*((-tt44*tt39)+2*tt9*tt16+2*tt7*tt25+2*tt1*tt15))
jac(1,3) = volume(1,1)*(2.5E-1*lam(1,1)*tt48*tt39*tt40+5.0E-1*miu&
&(1,1)*((-tt48*tt39)+2*tt9*tt22+2*tt7*tt26+2*tt1*tt21))
jac(1,4) = volume(1,1)*(2.5E-1*lam(1,1)*tt52*tt39*tt40+5.0E-1*miu&
&(1,1)*((-tt52*tt39)+2*D(1,3)*tt10+2*D(1,2)*tt8+2*D(1,1)*tt6))
jac(1,5) = volume(1,1)*(2.5E-1*lam(1,1)*tt56*tt39*tt40+5.0E-1*miu&
&(1,1)*((-tt56*tt39)+2*D(1,3)*tt16+2*D(1,2)*tt25+2*D(1,1)*tt15))
jac(1,6) = volume(1,1)*(2.5E-1*lam(1,1)*tt60*tt39*tt40+5.0E-1*miu&
&(1,1)*((-tt60*tt39)+2*D(1,3)*tt22+2*D(1,2)*tt26+2*D(1,1)*tt21))
jac(1,7) = volume(1,1)*(2.5E-1*lam(1,1)*tt64*tt39*tt40+5.0E-1*miu&
&(1,1)*((-tt64*tt39)+2*D(2,3)*tt10+2*D(2,2)*tt8+2*D(2,1)*tt6))
jac(1,8) = volume(1,1)*(2.5E-1*lam(1,1)*tt68*tt39*tt40+5.0E-1*miu&
&(1,1)*((-tt68*tt39)+2*D(2,3)*tt16+2*D(2,2)*tt25+2*D(2,1)*tt15))
jac(1,9) = volume(1,1)*(2.5E-1*lam(1,1)*tt72*tt39*tt40+5.0E-1*miu&
&(1,1)*((-tt72*tt39)+2*D(2,3)*tt22+2*D(2,2)*tt26+2*D(2,1)*tt21))
jac(1,10) = volume(1,1)*(2.5E-1*lam(1,1)*tt76*tt39*tt40+5.0E-1*mi&
&u(1,1)*((-tt76*tt39)+2*D(3,3)*tt10+2*D(3,2)*tt8+2*D(3,1)*tt6))
jac(1,11) = volume(1,1)*(2.5E-1*lam(1,1)*tt80*tt39*tt40+5.0E-1*mi&
&u(1,1)*((-tt80*tt39)+2*D(3,3)*tt16+2*D(3,2)*tt25+2*D(3,1)*tt15))
jac(1,12) = volume(1,1)*(2.5E-1*lam(1,1)*tt84*tt39*tt40+5.0E-1*mi&
&u(1,1)*((-tt84*tt39)+2*D(3,3)*tt22+2*D(3,2)*tt26+2*D(3,1)*tt21))
END 
SUBROUTINE tet_neohookean_hes(hes, X, D, volume, lam, miu) 
IMPLICIT NONE 
REAL(KIND=8) hes(12, 12) 
REAL(KIND=8) X(3, 4) 
REAL(KIND=8) D(3, 3) 
REAL(KIND=8) volume(1, 1) 
REAL(KIND=8) lam(1, 1) 
REAL(KIND=8) miu(1, 1) 
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
tt1 = -X(1,1)
tt2 = X(1,2)+tt1
tt3 = X(1,3)+tt1
tt4 = X(1,4)+tt1
tt5 = tt4*D(3,1)+tt3*D(2,1)+D(1,1)*tt2
tt6 = tt4*D(3,3)+tt3*D(2,3)+tt2*D(1,3)
tt7 = -X(2,1)
tt8 = X(2,2)+tt7
tt9 = X(2,3)+tt7
tt10 = X(2,4)+tt7
tt11 = tt10*D(3,1)+D(2,1)*tt9+D(1,1)*tt8
tt12 = tt10*D(3,3)+D(2,3)*tt9+D(1,3)*tt8
tt13 = -X(3,1)
tt14 = X(3,2)+tt13
tt15 = X(3,3)+tt13
tt16 = X(3,4)+tt13
tt17 = D(3,1)*tt16+D(2,1)*tt15+D(1,1)*tt14
tt18 = D(3,3)*tt16+D(2,3)*tt15+D(1,3)*tt14
tt19 = tt17*tt18+tt11*tt12+tt5*tt6
tt20 = tt4*D(3,2)+tt3*D(2,2)+D(1,2)*tt2
tt21 = (-D(3,3))-D(2,3)-D(1,3)
tt22 = (-D(3,2))-D(2,2)-D(1,2)
tt23 = tt22*tt6+tt20*tt21
tt24 = tt10*D(3,2)+D(2,2)*tt9+D(1,2)*tt8
tt25 = D(3,2)*tt16+D(2,2)*tt15+D(1,2)*tt14
tt26 = tt17*tt25+tt11*tt24+tt5*tt20
tt27 = (-D(3,1))-D(2,1)-D(1,1)
tt28 = tt27*tt20+tt5*tt22
tt29 = tt25*tt18+tt24*tt12+tt20*tt6
tt30 = tt27*tt6+tt5*tt21
tt31 = tt25**2+tt24**2+tt20**2
tt32 = (-tt30*tt31)+tt28*tt29-2*tt22*tt20*tt19+tt23*tt26
tt33 = tt26*tt29-tt19*tt31
tt34 = tt17**2+tt11**2+tt5**2
tt35 = tt18**2+tt12**2+tt6**2
tt36 = 2*tt22*tt20*tt35+2*tt21*tt6*tt31-2*tt23*tt29
tt37 = tt28*tt35-tt30*tt29-tt23*tt19+2*tt21*tt6*tt26
tt38 = tt26*tt35-tt19*tt29
tt39 = tt31*tt35-tt29**2
tt40 = 2*tt27*tt5*tt39-tt28*tt38-tt26*tt37+tt34*tt36+tt30*tt33+tt&
&19*tt32
tt41 = tt40**2
tt42 = tt34*tt39-tt26*tt38+tt19*tt33
tt43 = 1/tt42**2
tt44 = tt27**2
tt45 = 2*tt44
tt46 = tt22**2
tt47 = 2*tt46
tt48 = tt21**2
tt49 = 2*tt48
tt50 = 2*tt22*tt21*tt26
tt51 = -2*tt46*tt19
tt52 = 2*tt27*tt22*tt29
tt53 = -2*tt27*tt21*tt31
tt54 = 2*tt27*tt21*tt33
tt55 = -4*tt22*tt21*tt29
tt56 = 2*tt48*tt31
tt57 = 2*tt46*tt35
tt58 = 2*tt48*tt26
tt59 = -2*tt22*tt21*tt19
tt60 = -2*tt27*tt21*tt29
tt61 = 2*tt27*tt22*tt35
tt62 = -2*tt27*tt22*tt38
tt63 = 2*tt44*tt39
tt64 = tt63+tt62-2*tt28*tt37+4*tt27*tt5*tt36-tt26*(tt61+tt60+tt59&
&+tt58-2*tt30*tt23+4*tt28*tt21*tt6)+tt34*(tt57+tt56+tt55-2*tt23**2&
&+8*tt22*tt20*tt21*tt6)+tt54+2*tt30*tt32+tt19*(tt53+tt52+tt51+tt50&
&+2*tt28*tt23-4*tt22*tt20*tt30)
tt65 = 1/tt42
tt66 = log(tt42)
tt67 = tt22*tt12+tt24*tt21
tt68 = tt27*tt24+tt11*tt22
tt69 = tt27*tt12+tt11*tt21
tt70 = (-tt69*tt31)+tt68*tt29-2*tt22*tt24*tt19+tt67*tt26
tt71 = 2*tt22*tt24*tt35+2*tt21*tt12*tt31-2*tt67*tt29
tt72 = tt68*tt35-tt69*tt29-tt67*tt19+2*tt21*tt12*tt26
tt73 = 2*tt27*tt11*tt39-tt68*tt38-tt26*tt72+tt34*tt71+tt69*tt33+t&
&t19*tt70
tt74 = (-tt28*tt72)-tt68*tt37+2*tt27*tt5*tt71+2*tt27*tt11*tt36+tt&
&30*tt70+tt69*tt32+((-2*tt23*tt67)+4*tt22*tt20*tt21*tt12+4*tt22*tt&
&24*tt21*tt6)*tt34+(tt28*tt67-2*tt22*tt20*tt69+tt68*tt23-2*tt22*tt&
&24*tt30)*tt19-((-tt30*tt67)-tt23*tt69+2*tt28*tt21*tt12+2*tt68*tt2&
&1*tt6)*tt26
tt75 = volume(1,1)*(2.5E-1*lam(1,1)*tt74*tt65*tt66-2.5E-1*lam(1,1&
&)*tt40*tt73*tt43*tt66+5.0E-1*miu(1,1)*(tt40*tt73*tt43-tt74*tt65)+&
&2.5E-1*lam(1,1)*tt40*tt73*tt43)
tt76 = tt22*tt18+tt21*tt25
tt77 = tt27*tt25+tt22*tt17
tt78 = tt27*tt18+tt21*tt17
tt79 = (-tt78*tt31)+tt77*tt29-2*tt22*tt25*tt19+tt26*tt76
tt80 = 2*tt22*tt25*tt35+2*tt21*tt18*tt31-2*tt76*tt29
tt81 = tt77*tt35-tt78*tt29-tt76*tt19+2*tt21*tt18*tt26
tt82 = 2*tt27*tt17*tt39-tt77*tt38-tt26*tt81+tt34*tt80+tt78*tt33+t&
&t19*tt79
tt83 = (-tt28*tt81)+2*tt27*tt5*tt80-tt77*tt37+2*tt27*tt17*tt36+tt&
&30*tt79+tt78*tt32+((-2*tt23*tt76)+4*tt22*tt20*tt21*tt18+4*tt22*tt&
&21*tt6*tt25)*tt34-tt26*((-tt30*tt76)-tt23*tt78+2*tt21*tt6*tt77+2*&
&tt28*tt21*tt18)+tt19*(tt28*tt76-2*tt22*tt20*tt78+tt23*tt77-2*tt22&
&*tt30*tt25)
tt84 = volume(1,1)*(2.5E-1*lam(1,1)*tt83*tt65*tt66-2.5E-1*lam(1,1&
&)*tt40*tt82*tt43*tt66+5.0E-1*miu(1,1)*(tt40*tt82*tt43-tt83*tt65)+&
&2.5E-1*lam(1,1)*tt40*tt82*tt43)
tt85 = D(1,2)*tt6+D(1,3)*tt20
tt86 = D(1,1)*tt20+D(1,2)*tt5
tt87 = D(1,1)*tt6+D(1,3)*tt5
tt88 = (-tt87*tt31)+tt86*tt29-2*D(1,2)*tt20*tt19+tt85*tt26
tt89 = 2*D(1,2)*tt20*tt35+2*D(1,3)*tt6*tt31-2*tt85*tt29
tt90 = tt86*tt35-tt87*tt29-tt85*tt19+2*D(1,3)*tt6*tt26
tt91 = 2*D(1,1)*tt5*tt39-tt86*tt38-tt26*tt90+tt34*tt89+tt87*tt33+&
&tt19*tt88
tt92 = 2*D(1,1)*tt27
tt93 = 2*D(1,2)*tt22
tt94 = 2*D(1,3)*tt21
tt95 = D(1,2)*tt21+D(1,3)*tt22
tt96 = tt95*tt26
tt97 = -2*D(1,2)*tt22*tt19
tt98 = D(1,1)*tt22+D(1,2)*tt27
tt99 = tt98*tt29
tt100 = D(1,1)*tt21+D(1,3)*tt27
tt101 = -tt100*tt31
tt102 = tt100*tt33
tt103 = 2*D(1,3)*tt21*tt26
tt104 = -tt95*tt19
tt105 = -tt100*tt29
tt106 = tt98*tt35
tt107 = -2*tt95*tt29
tt108 = 2*D(1,3)*tt21*tt31
tt109 = 2*D(1,2)*tt22*tt35
tt110 = -tt98*tt38
tt111 = 2*D(1,1)*tt27*tt39
tt112 = tt111+tt110-tt86*tt37-tt28*tt90+2*D(1,1)*tt5*tt36+2*tt27*&
&tt5*tt89+tt34*(tt109+tt108+tt107-2*tt85*tt23+4*D(1,2)*tt20*tt21*t&
&t6+4*D(1,3)*tt22*tt20*tt6)-tt26*(tt106+tt105+tt104+tt103-tt87*tt2&
&3-tt85*tt30+2*tt86*tt21*tt6+2*D(1,3)*tt28*tt6)+tt102+tt87*tt32+tt&
&30*tt88+tt19*(tt101+tt99+tt97+tt96+tt86*tt23-2*D(1,2)*tt20*tt30+t&
&t28*tt85-2*tt22*tt20*tt87)
tt113 = volume(1,1)*(2.5E-1*lam(1,1)*tt112*tt65*tt66-2.5E-1*lam(1&
&,1)*tt91*tt40*tt43*tt66+5.0E-1*miu(1,1)*((-tt112*tt65)+tt91*tt40*&
&tt43+tt94+tt93+tt92)+2.5E-1*lam(1,1)*tt91*tt40*tt43)
tt114 = D(1,2)*tt12+D(1,3)*tt24
tt115 = D(1,1)*tt24+D(1,2)*tt11
tt116 = D(1,1)*tt12+D(1,3)*tt11
tt117 = (-tt116*tt31)+tt115*tt29-2*D(1,2)*tt24*tt19+tt114*tt26
tt118 = 2*D(1,2)*tt24*tt35+2*D(1,3)*tt12*tt31-2*tt114*tt29
tt119 = tt115*tt35-tt116*tt29-tt114*tt19+2*D(1,3)*tt12*tt26
tt120 = 2*D(1,1)*tt11*tt39-tt115*tt38-tt26*tt119+tt34*tt118+tt116&
&*tt33+tt19*tt117
tt121 = (-tt28*tt119)-tt115*tt37+2*tt27*tt5*tt118+2*D(1,1)*tt11*t&
&t36+tt30*tt117+tt116*tt32+((-2*tt23*tt114)+4*D(1,3)*tt22*tt20*tt1&
&2+4*D(1,2)*tt24*tt21*tt6)*tt34+(tt28*tt114-2*tt22*tt20*tt116+tt11&
&5*tt23-2*D(1,2)*tt24*tt30)*tt19-((-tt30*tt114)-tt23*tt116+2*D(1,3&
&)*tt28*tt12+2*tt115*tt21*tt6)*tt26
tt122 = volume(1,1)*(2.5E-1*lam(1,1)*tt121*tt65*tt66-2.5E-1*lam(1&
&,1)*tt40*tt120*tt43*tt66+5.0E-1*miu(1,1)*(tt40*tt120*tt43-tt121*t&
&t65)+2.5E-1*lam(1,1)*tt40*tt120*tt43)
tt123 = D(1,2)*tt18+D(1,3)*tt25
tt124 = D(1,1)*tt25+D(1,2)*tt17
tt125 = D(1,1)*tt18+D(1,3)*tt17
tt126 = (-tt125*tt31)+tt124*tt29-2*D(1,2)*tt25*tt19+tt26*tt123
tt127 = 2*D(1,2)*tt25*tt35+2*D(1,3)*tt18*tt31-2*tt123*tt29
tt128 = tt124*tt35-tt125*tt29-tt123*tt19+2*D(1,3)*tt18*tt26
tt129 = 2*D(1,1)*tt17*tt39-tt124*tt38-tt26*tt128+tt34*tt127+tt125&
&*tt33+tt19*tt126
tt130 = (-tt28*tt128)+2*tt27*tt5*tt127-tt124*tt37+2*D(1,1)*tt17*t&
&t36+tt30*tt126+tt125*tt32+((-2*tt23*tt123)+4*D(1,3)*tt22*tt20*tt1&
&8+4*D(1,2)*tt21*tt6*tt25)*tt34-tt26*((-tt30*tt123)-tt23*tt125+2*t&
&t21*tt6*tt124+2*D(1,3)*tt28*tt18)+tt19*(tt28*tt123-2*tt22*tt20*tt&
&125+tt23*tt124-2*D(1,2)*tt30*tt25)
tt131 = volume(1,1)*(2.5E-1*lam(1,1)*tt130*tt65*tt66-2.5E-1*lam(1&
&,1)*tt40*tt129*tt43*tt66+5.0E-1*miu(1,1)*(tt40*tt129*tt43-tt130*t&
&t65)+2.5E-1*lam(1,1)*tt40*tt129*tt43)
tt132 = D(2,2)*tt6+D(2,3)*tt20
tt133 = D(2,1)*tt20+D(2,2)*tt5
tt134 = D(2,1)*tt6+D(2,3)*tt5
tt135 = (-tt134*tt31)+tt133*tt29-2*D(2,2)*tt20*tt19+tt132*tt26
tt136 = 2*D(2,2)*tt20*tt35+2*D(2,3)*tt6*tt31-2*tt132*tt29
tt137 = tt133*tt35-tt134*tt29-tt132*tt19+2*D(2,3)*tt6*tt26
tt138 = 2*D(2,1)*tt5*tt39-tt133*tt38-tt26*tt137+tt34*tt136+tt134*&
&tt33+tt19*tt135
tt139 = 2*D(2,1)*tt27
tt140 = 2*D(2,2)*tt22
tt141 = 2*D(2,3)*tt21
tt142 = D(2,2)*tt21+D(2,3)*tt22
tt143 = tt142*tt26
tt144 = -2*D(2,2)*tt22*tt19
tt145 = D(2,1)*tt22+D(2,2)*tt27
tt146 = tt145*tt29
tt147 = D(2,1)*tt21+D(2,3)*tt27
tt148 = -tt147*tt31
tt149 = tt147*tt33
tt150 = 2*D(2,3)*tt21*tt26
tt151 = -tt142*tt19
tt152 = -tt147*tt29
tt153 = tt145*tt35
tt154 = -2*tt142*tt29
tt155 = 2*D(2,3)*tt21*tt31
tt156 = 2*D(2,2)*tt22*tt35
tt157 = -tt145*tt38
tt158 = 2*D(2,1)*tt27*tt39
tt159 = tt158+tt157-tt133*tt37-tt28*tt137+2*D(2,1)*tt5*tt36+2*tt2&
&7*tt5*tt136+tt34*(tt156+tt155+tt154-2*tt132*tt23+4*D(2,2)*tt20*tt&
&21*tt6+4*D(2,3)*tt22*tt20*tt6)-tt26*(tt153+tt152+tt151+tt150-tt13&
&4*tt23-tt132*tt30+2*tt133*tt21*tt6+2*D(2,3)*tt28*tt6)+tt149+tt134&
&*tt32+tt30*tt135+tt19*(tt148+tt146+tt144+tt143+tt133*tt23-2*D(2,2&
&)*tt20*tt30+tt28*tt132-2*tt22*tt20*tt134)
tt160 = volume(1,1)*(2.5E-1*lam(1,1)*tt159*tt65*tt66-2.5E-1*lam(1&
&,1)*tt138*tt40*tt43*tt66+5.0E-1*miu(1,1)*((-tt159*tt65)+tt138*tt4&
&0*tt43+tt141+tt140+tt139)+2.5E-1*lam(1,1)*tt138*tt40*tt43)
tt161 = D(2,2)*tt12+D(2,3)*tt24
tt162 = D(2,1)*tt24+D(2,2)*tt11
tt163 = D(2,1)*tt12+D(2,3)*tt11
tt164 = (-tt163*tt31)+tt162*tt29-2*D(2,2)*tt24*tt19+tt161*tt26
tt165 = 2*D(2,2)*tt24*tt35+2*D(2,3)*tt12*tt31-2*tt161*tt29
tt166 = tt162*tt35-tt163*tt29-tt161*tt19+2*D(2,3)*tt12*tt26
tt167 = 2*D(2,1)*tt11*tt39-tt162*tt38-tt26*tt166+tt34*tt165+tt163&
&*tt33+tt19*tt164
tt168 = (-tt28*tt166)-tt162*tt37+2*tt27*tt5*tt165+2*D(2,1)*tt11*t&
&t36+tt30*tt164+tt163*tt32+((-2*tt23*tt161)+4*D(2,3)*tt22*tt20*tt1&
&2+4*D(2,2)*tt24*tt21*tt6)*tt34+(tt28*tt161-2*tt22*tt20*tt163+tt16&
&2*tt23-2*D(2,2)*tt24*tt30)*tt19-((-tt30*tt161)-tt23*tt163+2*D(2,3&
&)*tt28*tt12+2*tt162*tt21*tt6)*tt26
tt169 = volume(1,1)*(2.5E-1*lam(1,1)*tt168*tt65*tt66-2.5E-1*lam(1&
&,1)*tt40*tt167*tt43*tt66+5.0E-1*miu(1,1)*(tt40*tt167*tt43-tt168*t&
&t65)+2.5E-1*lam(1,1)*tt40*tt167*tt43)
tt170 = D(2,2)*tt18+D(2,3)*tt25
tt171 = D(2,1)*tt25+D(2,2)*tt17
tt172 = D(2,1)*tt18+D(2,3)*tt17
tt173 = (-tt172*tt31)+tt171*tt29-2*D(2,2)*tt25*tt19+tt26*tt170
tt174 = 2*D(2,2)*tt25*tt35+2*D(2,3)*tt18*tt31-2*tt170*tt29
tt175 = tt171*tt35-tt172*tt29-tt170*tt19+2*D(2,3)*tt18*tt26
tt176 = 2*D(2,1)*tt17*tt39-tt171*tt38-tt26*tt175+tt34*tt174+tt172&
&*tt33+tt19*tt173
tt177 = (-tt28*tt175)+2*tt27*tt5*tt174-tt171*tt37+2*D(2,1)*tt17*t&
&t36+tt30*tt173+tt172*tt32+((-2*tt23*tt170)+4*D(2,3)*tt22*tt20*tt1&
&8+4*D(2,2)*tt21*tt6*tt25)*tt34-tt26*((-tt30*tt170)-tt23*tt172+2*t&
&t21*tt6*tt171+2*D(2,3)*tt28*tt18)+tt19*(tt28*tt170-2*tt22*tt20*tt&
&172+tt23*tt171-2*D(2,2)*tt30*tt25)
tt178 = volume(1,1)*(2.5E-1*lam(1,1)*tt177*tt65*tt66-2.5E-1*lam(1&
&,1)*tt40*tt176*tt43*tt66+5.0E-1*miu(1,1)*(tt40*tt176*tt43-tt177*t&
&t65)+2.5E-1*lam(1,1)*tt40*tt176*tt43)
tt179 = D(3,2)*tt6+tt20*D(3,3)
tt180 = D(3,1)*tt20+tt5*D(3,2)
tt181 = D(3,1)*tt6+tt5*D(3,3)
tt182 = (-tt181*tt31)+tt180*tt29-2*D(3,2)*tt20*tt19+tt179*tt26
tt183 = 2*D(3,2)*tt20*tt35+2*D(3,3)*tt6*tt31-2*tt179*tt29
tt184 = tt180*tt35-tt181*tt29-tt179*tt19+2*D(3,3)*tt6*tt26
tt185 = 2*D(3,1)*tt5*tt39-tt180*tt38-tt26*tt184+tt34*tt183+tt181*&
&tt33+tt19*tt182
tt186 = 2*tt27*D(3,1)
tt187 = 2*tt22*D(3,2)
tt188 = 2*tt21*D(3,3)
tt189 = tt22*D(3,3)+D(3,2)*tt21
tt190 = tt189*tt26
tt191 = -2*tt22*D(3,2)*tt19
tt192 = tt27*D(3,2)+D(3,1)*tt22
tt193 = tt192*tt29
tt194 = tt27*D(3,3)+D(3,1)*tt21
tt195 = -tt194*tt31
tt196 = tt194*tt33
tt197 = -2*tt189*tt29
tt198 = 2*tt21*D(3,3)*tt31
tt199 = 2*tt22*D(3,2)*tt35
tt200 = 2*tt21*D(3,3)*tt26
tt201 = -tt189*tt19
tt202 = -tt194*tt29
tt203 = tt192*tt35
tt204 = -tt192*tt38
tt205 = 2*tt27*D(3,1)*tt39
tt206 = tt205+tt204-tt28*tt184-tt180*tt37-tt26*(tt203+tt202+tt201&
&+tt200-tt30*tt179-tt181*tt23+2*tt28*D(3,3)*tt6+2*tt180*tt21*tt6)+&
&2*tt27*tt5*tt183+2*D(3,1)*tt5*tt36+tt34*(tt199+tt198+tt197-2*tt23&
&*tt179+4*tt22*tt20*D(3,3)*tt6+4*D(3,2)*tt20*tt21*tt6)+tt196+tt30*&
&tt182+tt181*tt32+tt19*(tt195+tt193+tt191+tt190+tt28*tt179+tt180*t&
&t23-2*tt22*tt20*tt181-2*D(3,2)*tt20*tt30)
tt207 = volume(1,1)*(2.5E-1*lam(1,1)*tt206*tt65*tt66-2.5E-1*lam(1&
&,1)*tt40*tt185*tt43*tt66+5.0E-1*miu(1,1)*((-tt206*tt65)+tt40*tt18&
&5*tt43+tt188+tt187+tt186)+2.5E-1*lam(1,1)*tt40*tt185*tt43)
tt208 = D(3,2)*tt12+tt24*D(3,3)
tt209 = D(3,1)*tt24+tt11*D(3,2)
tt210 = D(3,1)*tt12+tt11*D(3,3)
tt211 = (-tt210*tt31)+tt209*tt29-2*D(3,2)*tt24*tt19+tt208*tt26
tt212 = 2*D(3,2)*tt24*tt35+2*D(3,3)*tt12*tt31-2*tt208*tt29
tt213 = tt209*tt35-tt210*tt29-tt208*tt19+2*D(3,3)*tt12*tt26
tt214 = 2*D(3,1)*tt11*tt39-tt209*tt38-tt26*tt213+tt34*tt212+tt210&
&*tt33+tt19*tt211
tt215 = (-tt28*tt213)-tt209*tt37+2*tt27*tt5*tt212+2*D(3,1)*tt11*t&
&t36+tt30*tt211+tt210*tt32+((-2*tt23*tt208)+4*tt22*tt20*D(3,3)*tt1&
&2+4*D(3,2)*tt24*tt21*tt6)*tt34+(tt28*tt208-2*tt22*tt20*tt210+tt20&
&9*tt23-2*D(3,2)*tt24*tt30)*tt19-((-tt30*tt208)-tt23*tt210+2*tt28*&
&D(3,3)*tt12+2*tt209*tt21*tt6)*tt26
tt216 = volume(1,1)*(2.5E-1*lam(1,1)*tt215*tt65*tt66-2.5E-1*lam(1&
&,1)*tt40*tt214*tt43*tt66+5.0E-1*miu(1,1)*(tt40*tt214*tt43-tt215*t&
&t65)+2.5E-1*lam(1,1)*tt40*tt214*tt43)
tt217 = D(3,2)*tt18+D(3,3)*tt25
tt218 = D(3,1)*tt25+D(3,2)*tt17
tt219 = D(3,1)*tt18+D(3,3)*tt17
tt220 = (-tt219*tt31)+tt218*tt29-2*D(3,2)*tt25*tt19+tt26*tt217
tt221 = 2*D(3,2)*tt25*tt35+2*D(3,3)*tt18*tt31-2*tt217*tt29
tt222 = tt218*tt35-tt219*tt29-tt217*tt19+2*D(3,3)*tt18*tt26
tt223 = 2*D(3,1)*tt17*tt39-tt218*tt38-tt26*tt222+tt34*tt221+tt219&
&*tt33+tt19*tt220
tt224 = (-tt28*tt222)+2*tt27*tt5*tt221-tt218*tt37+2*D(3,1)*tt17*t&
&t36+tt30*tt220+tt219*tt32+((-2*tt23*tt217)+4*tt22*tt20*D(3,3)*tt1&
&8+4*D(3,2)*tt21*tt6*tt25)*tt34-tt26*((-tt30*tt217)-tt23*tt219+2*t&
&t21*tt6*tt218+2*tt28*D(3,3)*tt18)+tt19*(tt28*tt217-2*tt22*tt20*tt&
&219+tt23*tt218-2*D(3,2)*tt30*tt25)
tt225 = volume(1,1)*(2.5E-1*lam(1,1)*tt224*tt65*tt66-2.5E-1*lam(1&
&,1)*tt40*tt223*tt43*tt66+5.0E-1*miu(1,1)*(tt40*tt223*tt43-tt224*t&
&t65)+2.5E-1*lam(1,1)*tt40*tt223*tt43)
tt226 = tt73**2
tt227 = tt63+tt62-2*tt68*tt72+4*tt27*tt11*tt71-tt26*(tt61+tt60+tt&
&59+tt58-2*tt69*tt67+4*tt68*tt21*tt12)+tt34*(tt57+tt56+tt55-2*tt67&
&**2+8*tt22*tt24*tt21*tt12)+tt54+2*tt69*tt70+tt19*(tt53+tt52+tt51+&
&tt50+2*tt68*tt67-4*tt22*tt24*tt69)
tt228 = (-tt68*tt81)+2*tt27*tt11*tt80-tt77*tt72+2*tt27*tt17*tt71+&
&tt69*tt79+tt78*tt70+((-2*tt67*tt76)+4*tt22*tt24*tt21*tt18+4*tt22*&
&tt21*tt12*tt25)*tt34-tt26*((-tt69*tt76)-tt67*tt78+2*tt21*tt12*tt7&
&7+2*tt68*tt21*tt18)+tt19*(tt68*tt76-2*tt22*tt24*tt78+tt67*tt77-2*&
&tt22*tt69*tt25)
tt229 = volume(1,1)*(2.5E-1*lam(1,1)*tt228*tt65*tt66-2.5E-1*lam(1&
&,1)*tt73*tt82*tt43*tt66+5.0E-1*miu(1,1)*(tt73*tt82*tt43-tt228*tt6&
&5)+2.5E-1*lam(1,1)*tt73*tt82*tt43)
tt230 = (-tt86*tt72)-tt68*tt90+2*D(1,1)*tt5*tt71+2*tt27*tt11*tt89&
&+tt87*tt70+tt69*tt88+((-2*tt85*tt67)+4*D(1,2)*tt20*tt21*tt12+4*D(&
&1,3)*tt22*tt24*tt6)*tt34+(tt86*tt67-2*D(1,2)*tt20*tt69+tt68*tt85-&
&2*tt22*tt24*tt87)*tt19-((-tt87*tt67)-tt85*tt69+2*tt86*tt21*tt12+2&
&*D(1,3)*tt68*tt6)*tt26
tt231 = volume(1,1)*(2.5E-1*lam(1,1)*tt230*tt65*tt66-2.5E-1*lam(1&
&,1)*tt91*tt73*tt43*tt66+5.0E-1*miu(1,1)*(tt91*tt73*tt43-tt230*tt6&
&5)+2.5E-1*lam(1,1)*tt91*tt73*tt43)
tt232 = tt111+tt110-tt115*tt72-tt68*tt119+2*D(1,1)*tt11*tt71+2*tt&
&27*tt11*tt118+tt34*(tt109+tt108+tt107-2*tt114*tt67+4*D(1,2)*tt24*&
&tt21*tt12+4*D(1,3)*tt22*tt24*tt12)-tt26*(tt106+tt105+tt104+tt103-&
&tt116*tt67-tt114*tt69+2*tt115*tt21*tt12+2*D(1,3)*tt68*tt12)+tt102&
&+tt116*tt70+tt69*tt117+tt19*(tt101+tt99+tt97+tt96+tt115*tt67-2*D(&
&1,2)*tt24*tt69+tt68*tt114-2*tt22*tt24*tt116)
tt233 = volume(1,1)*(2.5E-1*lam(1,1)*tt232*tt65*tt66-2.5E-1*lam(1&
&,1)*tt120*tt73*tt43*tt66+5.0E-1*miu(1,1)*((-tt232*tt65)+tt120*tt7&
&3*tt43+tt94+tt93+tt92)+2.5E-1*lam(1,1)*tt120*tt73*tt43)
tt234 = (-tt68*tt128)+2*tt27*tt11*tt127-tt124*tt72+2*D(1,1)*tt17*&
&tt71+tt69*tt126+tt125*tt70+((-2*tt67*tt123)+4*D(1,3)*tt22*tt24*tt&
&18+4*D(1,2)*tt21*tt12*tt25)*tt34-tt26*((-tt69*tt123)-tt67*tt125+2&
&*tt21*tt12*tt124+2*D(1,3)*tt68*tt18)+tt19*(tt68*tt123-2*tt22*tt24&
&*tt125+tt67*tt124-2*D(1,2)*tt69*tt25)
tt235 = volume(1,1)*(2.5E-1*lam(1,1)*tt234*tt65*tt66-2.5E-1*lam(1&
&,1)*tt73*tt129*tt43*tt66+5.0E-1*miu(1,1)*(tt73*tt129*tt43-tt234*t&
&t65)+2.5E-1*lam(1,1)*tt73*tt129*tt43)
tt236 = (-tt133*tt72)-tt68*tt137+2*D(2,1)*tt5*tt71+2*tt27*tt11*tt&
&136+tt134*tt70+tt69*tt135+((-2*tt132*tt67)+4*D(2,2)*tt20*tt21*tt1&
&2+4*D(2,3)*tt22*tt24*tt6)*tt34+(tt133*tt67-2*D(2,2)*tt20*tt69+tt6&
&8*tt132-2*tt22*tt24*tt134)*tt19-((-tt134*tt67)-tt132*tt69+2*tt133&
&*tt21*tt12+2*D(2,3)*tt68*tt6)*tt26
tt237 = volume(1,1)*(2.5E-1*lam(1,1)*tt236*tt65*tt66-2.5E-1*lam(1&
&,1)*tt138*tt73*tt43*tt66+5.0E-1*miu(1,1)*(tt138*tt73*tt43-tt236*t&
&t65)+2.5E-1*lam(1,1)*tt138*tt73*tt43)
tt238 = tt158+tt157-tt162*tt72-tt68*tt166+2*D(2,1)*tt11*tt71+2*tt&
&27*tt11*tt165+tt34*(tt156+tt155+tt154-2*tt161*tt67+4*D(2,2)*tt24*&
&tt21*tt12+4*D(2,3)*tt22*tt24*tt12)-tt26*(tt153+tt152+tt151+tt150-&
&tt163*tt67-tt161*tt69+2*tt162*tt21*tt12+2*D(2,3)*tt68*tt12)+tt149&
&+tt163*tt70+tt69*tt164+tt19*(tt148+tt146+tt144+tt143+tt162*tt67-2&
&*D(2,2)*tt24*tt69+tt68*tt161-2*tt22*tt24*tt163)
tt239 = volume(1,1)*(2.5E-1*lam(1,1)*tt238*tt65*tt66-2.5E-1*lam(1&
&,1)*tt167*tt73*tt43*tt66+5.0E-1*miu(1,1)*((-tt238*tt65)+tt167*tt7&
&3*tt43+tt141+tt140+tt139)+2.5E-1*lam(1,1)*tt167*tt73*tt43)
tt240 = (-tt68*tt175)+2*tt27*tt11*tt174-tt171*tt72+2*D(2,1)*tt17*&
&tt71+tt69*tt173+tt172*tt70+((-2*tt67*tt170)+4*D(2,3)*tt22*tt24*tt&
&18+4*D(2,2)*tt21*tt12*tt25)*tt34-tt26*((-tt69*tt170)-tt67*tt172+2&
&*tt21*tt12*tt171+2*D(2,3)*tt68*tt18)+tt19*(tt68*tt170-2*tt22*tt24&
&*tt172+tt67*tt171-2*D(2,2)*tt69*tt25)
tt241 = volume(1,1)*(2.5E-1*lam(1,1)*tt240*tt65*tt66-2.5E-1*lam(1&
&,1)*tt73*tt176*tt43*tt66+5.0E-1*miu(1,1)*(tt73*tt176*tt43-tt240*t&
&t65)+2.5E-1*lam(1,1)*tt73*tt176*tt43)
tt242 = (-tt180*tt72)-tt68*tt184+2*D(3,1)*tt5*tt71+2*tt27*tt11*tt&
&183+tt181*tt70+tt69*tt182+((-2*tt179*tt67)+4*D(3,2)*tt20*tt21*tt1&
&2+4*tt22*tt24*D(3,3)*tt6)*tt34+(tt180*tt67-2*D(3,2)*tt20*tt69+tt6&
&8*tt179-2*tt22*tt24*tt181)*tt19-((-tt181*tt67)-tt179*tt69+2*tt180&
&*tt21*tt12+2*tt68*D(3,3)*tt6)*tt26
tt243 = volume(1,1)*(2.5E-1*lam(1,1)*tt242*tt65*tt66-2.5E-1*lam(1&
&,1)*tt185*tt73*tt43*tt66+5.0E-1*miu(1,1)*(tt185*tt73*tt43-tt242*t&
&t65)+2.5E-1*lam(1,1)*tt185*tt73*tt43)
tt244 = tt205+tt204-tt68*tt213-tt209*tt72-tt26*(tt203+tt202+tt201&
&+tt200-tt69*tt208-tt210*tt67+2*tt68*D(3,3)*tt12+2*tt209*tt21*tt12&
&)+2*tt27*tt11*tt212+2*D(3,1)*tt11*tt71+tt34*(tt199+tt198+tt197-2*&
&tt67*tt208+4*tt22*tt24*D(3,3)*tt12+4*D(3,2)*tt24*tt21*tt12)+tt196&
&+tt69*tt211+tt210*tt70+tt19*(tt195+tt193+tt191+tt190+tt68*tt208+t&
&t209*tt67-2*tt22*tt24*tt210-2*D(3,2)*tt24*tt69)
tt245 = volume(1,1)*(2.5E-1*lam(1,1)*tt244*tt65*tt66-2.5E-1*lam(1&
&,1)*tt73*tt214*tt43*tt66+5.0E-1*miu(1,1)*((-tt244*tt65)+tt73*tt21&
&4*tt43+tt188+tt187+tt186)+2.5E-1*lam(1,1)*tt73*tt214*tt43)
tt246 = (-tt68*tt222)+2*tt27*tt11*tt221-tt218*tt72+2*D(3,1)*tt17*&
&tt71+tt69*tt220+tt219*tt70+((-2*tt67*tt217)+4*tt22*tt24*D(3,3)*tt&
&18+4*D(3,2)*tt21*tt12*tt25)*tt34-tt26*((-tt69*tt217)-tt67*tt219+2&
&*tt21*tt12*tt218+2*tt68*D(3,3)*tt18)+tt19*(tt68*tt217-2*tt22*tt24&
&*tt219+tt67*tt218-2*D(3,2)*tt69*tt25)
tt247 = volume(1,1)*(2.5E-1*lam(1,1)*tt246*tt65*tt66-2.5E-1*lam(1&
&,1)*tt73*tt223*tt43*tt66+5.0E-1*miu(1,1)*(tt73*tt223*tt43-tt246*t&
&t65)+2.5E-1*lam(1,1)*tt73*tt223*tt43)
tt248 = tt82**2
tt249 = tt63+tt34*((-2*tt76**2)+tt57+tt56+tt55+8*tt22*tt21*tt25*t&
&t18)+tt62-2*tt77*tt81+4*tt27*tt17*tt80-tt26*(tt61+tt60+tt59-2*tt7&
&8*tt76+tt58+4*tt21*tt18*tt77)+tt54+2*tt78*tt79+tt19*(tt53+tt52+tt&
&51+2*tt77*tt76-4*tt22*tt25*tt78+tt50)
tt250 = (-tt86*tt81)+2*D(1,1)*tt5*tt80-tt77*tt90+2*tt27*tt17*tt89&
&+tt87*tt79+tt78*tt88+((-2*tt85*tt76)+4*D(1,2)*tt20*tt21*tt18+4*D(&
&1,3)*tt22*tt6*tt25)*tt34-tt26*((-tt87*tt76)-tt85*tt78+2*D(1,3)*tt&
&6*tt77+2*tt86*tt21*tt18)+tt19*(tt86*tt76-2*D(1,2)*tt20*tt78+tt85*&
&tt77-2*tt22*tt87*tt25)
tt251 = volume(1,1)*(2.5E-1*lam(1,1)*tt250*tt65*tt66-2.5E-1*lam(1&
&,1)*tt91*tt82*tt43*tt66+5.0E-1*miu(1,1)*(tt91*tt82*tt43-tt250*tt6&
&5)+2.5E-1*lam(1,1)*tt91*tt82*tt43)
tt252 = (-tt115*tt81)+2*D(1,1)*tt11*tt80-tt77*tt119+2*tt27*tt17*t&
&t118+tt116*tt79+tt78*tt117+((-2*tt114*tt76)+4*D(1,2)*tt24*tt21*tt&
&18+4*D(1,3)*tt22*tt12*tt25)*tt34-tt26*((-tt116*tt76)-tt114*tt78+2&
&*D(1,3)*tt12*tt77+2*tt115*tt21*tt18)+tt19*(tt115*tt76-2*D(1,2)*tt&
&24*tt78+tt114*tt77-2*tt22*tt116*tt25)
tt253 = volume(1,1)*(2.5E-1*lam(1,1)*tt252*tt65*tt66-2.5E-1*lam(1&
&,1)*tt120*tt82*tt43*tt66+5.0E-1*miu(1,1)*(tt120*tt82*tt43-tt252*t&
&t65)+2.5E-1*lam(1,1)*tt120*tt82*tt43)
tt254 = tt111+tt110-tt124*tt81-tt77*tt128+2*D(1,1)*tt17*tt80+2*tt&
&27*tt17*tt127+tt34*(tt109+tt108+tt107-2*tt123*tt76+4*D(1,2)*tt21*&
&tt25*tt18+4*D(1,3)*tt22*tt25*tt18)-tt26*(tt106+tt105+tt104-tt125*&
&tt76-tt123*tt78+tt103+2*D(1,3)*tt18*tt77+2*tt21*tt18*tt124)+tt102&
&+tt125*tt79+tt78*tt126+tt19*(tt101+tt99+tt97+tt124*tt76-2*D(1,2)*&
&tt25*tt78+tt77*tt123-2*tt22*tt25*tt125+tt96)
tt255 = volume(1,1)*(2.5E-1*lam(1,1)*tt254*tt65*tt66-2.5E-1*lam(1&
&,1)*tt129*tt82*tt43*tt66+5.0E-1*miu(1,1)*((-tt254*tt65)+tt129*tt8&
&2*tt43+tt94+tt93+tt92)+2.5E-1*lam(1,1)*tt129*tt82*tt43)
tt256 = (-tt133*tt81)+2*D(2,1)*tt5*tt80-tt77*tt137+2*tt27*tt17*tt&
&136+tt134*tt79+tt78*tt135+((-2*tt132*tt76)+4*D(2,2)*tt20*tt21*tt1&
&8+4*D(2,3)*tt22*tt6*tt25)*tt34-tt26*((-tt134*tt76)-tt132*tt78+2*D&
&(2,3)*tt6*tt77+2*tt133*tt21*tt18)+tt19*(tt133*tt76-2*D(2,2)*tt20*&
&tt78+tt132*tt77-2*tt22*tt134*tt25)
tt257 = volume(1,1)*(2.5E-1*lam(1,1)*tt256*tt65*tt66-2.5E-1*lam(1&
&,1)*tt138*tt82*tt43*tt66+5.0E-1*miu(1,1)*(tt138*tt82*tt43-tt256*t&
&t65)+2.5E-1*lam(1,1)*tt138*tt82*tt43)
tt258 = (-tt162*tt81)+2*D(2,1)*tt11*tt80-tt77*tt166+2*tt27*tt17*t&
&t165+tt163*tt79+tt78*tt164+((-2*tt161*tt76)+4*D(2,2)*tt24*tt21*tt&
&18+4*D(2,3)*tt22*tt12*tt25)*tt34-tt26*((-tt163*tt76)-tt161*tt78+2&
&*D(2,3)*tt12*tt77+2*tt162*tt21*tt18)+tt19*(tt162*tt76-2*D(2,2)*tt&
&24*tt78+tt161*tt77-2*tt22*tt163*tt25)
tt259 = volume(1,1)*(2.5E-1*lam(1,1)*tt258*tt65*tt66-2.5E-1*lam(1&
&,1)*tt167*tt82*tt43*tt66+5.0E-1*miu(1,1)*(tt167*tt82*tt43-tt258*t&
&t65)+2.5E-1*lam(1,1)*tt167*tt82*tt43)
tt260 = tt158+tt157-tt171*tt81-tt77*tt175+2*D(2,1)*tt17*tt80+2*tt&
&27*tt17*tt174+tt34*(tt156+tt155+tt154-2*tt170*tt76+4*D(2,2)*tt21*&
&tt25*tt18+4*D(2,3)*tt22*tt25*tt18)-tt26*(tt153+tt152+tt151-tt172*&
&tt76-tt170*tt78+tt150+2*D(2,3)*tt18*tt77+2*tt21*tt18*tt171)+tt149&
&+tt172*tt79+tt78*tt173+tt19*(tt148+tt146+tt144+tt171*tt76-2*D(2,2&
&)*tt25*tt78+tt77*tt170-2*tt22*tt25*tt172+tt143)
tt261 = volume(1,1)*(2.5E-1*lam(1,1)*tt260*tt65*tt66-2.5E-1*lam(1&
&,1)*tt176*tt82*tt43*tt66+5.0E-1*miu(1,1)*((-tt260*tt65)+tt176*tt8&
&2*tt43+tt141+tt140+tt139)+2.5E-1*lam(1,1)*tt176*tt82*tt43)
tt262 = (-tt180*tt81)+2*D(3,1)*tt5*tt80-tt77*tt184+2*tt27*tt17*tt&
&183+tt181*tt79+tt78*tt182+((-2*tt179*tt76)+4*D(3,2)*tt20*tt21*tt1&
&8+4*tt22*D(3,3)*tt6*tt25)*tt34-tt26*((-tt181*tt76)-tt179*tt78+2*D&
&(3,3)*tt6*tt77+2*tt180*tt21*tt18)+tt19*(tt180*tt76-2*D(3,2)*tt20*&
&tt78+tt179*tt77-2*tt22*tt181*tt25)
tt263 = volume(1,1)*(2.5E-1*lam(1,1)*tt262*tt65*tt66-2.5E-1*lam(1&
&,1)*tt185*tt82*tt43*tt66+5.0E-1*miu(1,1)*(tt185*tt82*tt43-tt262*t&
&t65)+2.5E-1*lam(1,1)*tt185*tt82*tt43)
tt264 = (-tt209*tt81)+2*D(3,1)*tt11*tt80-tt77*tt213+2*tt27*tt17*t&
&t212+tt210*tt79+tt78*tt211+((-2*tt208*tt76)+4*D(3,2)*tt24*tt21*tt&
&18+4*tt22*D(3,3)*tt12*tt25)*tt34-tt26*((-tt210*tt76)-tt208*tt78+2&
&*D(3,3)*tt12*tt77+2*tt209*tt21*tt18)+tt19*(tt209*tt76-2*D(3,2)*tt&
&24*tt78+tt208*tt77-2*tt22*tt210*tt25)
tt265 = volume(1,1)*(2.5E-1*lam(1,1)*tt264*tt65*tt66-2.5E-1*lam(1&
&,1)*tt214*tt82*tt43*tt66+5.0E-1*miu(1,1)*(tt214*tt82*tt43-tt264*t&
&t65)+2.5E-1*lam(1,1)*tt214*tt82*tt43)
tt266 = tt205+tt204-tt77*tt222-tt218*tt81+2*tt27*tt17*tt221+2*D(3&
&,1)*tt17*tt80-tt26*(tt203+tt202+tt201-tt78*tt217-tt219*tt76+tt200&
&+2*tt21*tt18*tt218+2*D(3,3)*tt18*tt77)+tt34*(tt199+tt198+tt197-2*&
&tt76*tt217+4*tt22*D(3,3)*tt25*tt18+4*D(3,2)*tt21*tt25*tt18)+tt196&
&+tt78*tt220+tt219*tt79+tt19*(tt195+tt193+tt191+tt77*tt217+tt218*t&
&t76-2*tt22*tt25*tt219-2*D(3,2)*tt25*tt78+tt190)
tt267 = volume(1,1)*(2.5E-1*lam(1,1)*tt266*tt65*tt66-2.5E-1*lam(1&
&,1)*tt82*tt223*tt43*tt66+5.0E-1*miu(1,1)*((-tt266*tt65)+tt82*tt22&
&3*tt43+tt188+tt187+tt186)+2.5E-1*lam(1,1)*tt82*tt223*tt43)
tt268 = tt91**2
tt269 = D(1,1)**2
tt270 = 2*tt269
tt271 = D(1,2)**2
tt272 = 2*tt271
tt273 = D(1,3)**2
tt274 = 2*tt273
tt275 = 2*D(1,2)*D(1,3)*tt26
tt276 = -2*tt271*tt19
tt277 = 2*D(1,1)*D(1,2)*tt29
tt278 = -2*D(1,1)*D(1,3)*tt31
tt279 = 2*D(1,1)*D(1,3)*tt33
tt280 = 2*tt273*tt26
tt281 = -2*D(1,2)*D(1,3)*tt19
tt282 = -2*D(1,1)*D(1,3)*tt29
tt283 = 2*D(1,1)*D(1,2)*tt35
tt284 = -4*D(1,2)*D(1,3)*tt29
tt285 = 2*tt273*tt31
tt286 = 2*tt271*tt35
tt287 = -2*D(1,1)*D(1,2)*tt38
tt288 = 2*tt269*tt39
tt289 = tt288+tt287-2*tt86*tt90+4*D(1,1)*tt5*tt89+tt34*(tt286+tt2&
&85+tt284-2*tt85**2+8*D(1,2)*D(1,3)*tt20*tt6)-tt26*(tt283+tt282+tt&
&281+tt280-2*tt87*tt85+4*D(1,3)*tt86*tt6)+tt279+2*tt87*tt88+tt19*(&
&tt278+tt277+tt276+tt275+2*tt86*tt85-4*D(1,2)*tt20*tt87)
tt290 = (-tt86*tt119)-tt115*tt90+2*D(1,1)*tt5*tt118+2*D(1,1)*tt11&
&*tt89+tt87*tt117+tt116*tt88+((-2*tt85*tt114)+4*D(1,2)*D(1,3)*tt20&
&*tt12+4*D(1,2)*D(1,3)*tt24*tt6)*tt34+(tt86*tt114-2*D(1,2)*tt20*tt&
&116+tt115*tt85-2*D(1,2)*tt24*tt87)*tt19-((-tt87*tt114)-tt85*tt116&
&+2*D(1,3)*tt86*tt12+2*D(1,3)*tt115*tt6)*tt26
tt291 = volume(1,1)*(2.5E-1*lam(1,1)*tt290*tt65*tt66-2.5E-1*lam(1&
&,1)*tt91*tt120*tt43*tt66+5.0E-1*miu(1,1)*(tt91*tt120*tt43-tt290*t&
&t65)+2.5E-1*lam(1,1)*tt91*tt120*tt43)
tt292 = (-tt86*tt128)+2*D(1,1)*tt5*tt127-tt124*tt90+2*D(1,1)*tt17&
&*tt89+tt87*tt126+tt125*tt88+((-2*tt85*tt123)+4*D(1,2)*D(1,3)*tt20&
&*tt18+4*D(1,2)*D(1,3)*tt6*tt25)*tt34-tt26*((-tt87*tt123)-tt85*tt1&
&25+2*D(1,3)*tt6*tt124+2*D(1,3)*tt86*tt18)+tt19*(tt86*tt123-2*D(1,&
&2)*tt20*tt125+tt85*tt124-2*D(1,2)*tt87*tt25)
tt293 = volume(1,1)*(2.5E-1*lam(1,1)*tt292*tt65*tt66-2.5E-1*lam(1&
&,1)*tt91*tt129*tt43*tt66+5.0E-1*miu(1,1)*(tt91*tt129*tt43-tt292*t&
&t65)+2.5E-1*lam(1,1)*tt91*tt129*tt43)
tt294 = 2*D(1,1)*D(2,1)
tt295 = 2*D(1,2)*D(2,2)
tt296 = 2*D(1,3)*D(2,3)
tt297 = D(1,2)*D(2,3)+D(1,3)*D(2,2)
tt298 = tt297*tt26
tt299 = -2*D(1,2)*D(2,2)*tt19
tt300 = D(1,1)*D(2,2)+D(1,2)*D(2,1)
tt301 = tt300*tt29
tt302 = D(1,1)*D(2,3)+D(1,3)*D(2,1)
tt303 = -tt302*tt31
tt304 = tt302*tt33
tt305 = -2*tt297*tt29
tt306 = 2*D(1,3)*D(2,3)*tt31
tt307 = 2*D(1,2)*D(2,2)*tt35
tt308 = 2*D(1,3)*D(2,3)*tt26
tt309 = -tt297*tt19
tt310 = -tt302*tt29
tt311 = tt300*tt35
tt312 = -tt300*tt38
tt313 = 2*D(1,1)*D(2,1)*tt39
tt314 = tt313+tt312-tt86*tt137-tt133*tt90+2*D(1,1)*tt5*tt136+2*D(&
&2,1)*tt5*tt89-tt26*(tt311+tt310+tt309+tt308-tt87*tt132-tt85*tt134&
&+2*D(1,3)*tt133*tt6+2*D(2,3)*tt86*tt6)+tt34*(tt307+tt306+tt305-2*&
&tt85*tt132+4*D(1,2)*D(2,3)*tt20*tt6+4*D(1,3)*D(2,2)*tt20*tt6)+tt3&
&04+tt87*tt135+tt134*tt88+tt19*(tt303+tt301+tt299+tt298+tt86*tt132&
&-2*D(1,2)*tt20*tt134+tt133*tt85-2*D(2,2)*tt20*tt87)
tt315 = volume(1,1)*(2.5E-1*lam(1,1)*tt314*tt65*tt66-2.5E-1*lam(1&
&,1)*tt91*tt138*tt43*tt66+5.0E-1*miu(1,1)*((-tt314*tt65)+tt91*tt13&
&8*tt43+tt296+tt295+tt294)+2.5E-1*lam(1,1)*tt91*tt138*tt43)
tt316 = (-tt86*tt166)-tt162*tt90+2*D(1,1)*tt5*tt165+2*D(2,1)*tt11&
&*tt89+tt87*tt164+tt163*tt88+((-2*tt85*tt161)+4*D(1,2)*D(2,3)*tt20&
&*tt12+4*D(1,3)*D(2,2)*tt24*tt6)*tt34+(tt86*tt161-2*D(1,2)*tt20*tt&
&163+tt162*tt85-2*D(2,2)*tt24*tt87)*tt19-((-tt87*tt161)-tt85*tt163&
&+2*D(2,3)*tt86*tt12+2*D(1,3)*tt162*tt6)*tt26
tt317 = volume(1,1)*(2.5E-1*lam(1,1)*tt316*tt65*tt66-2.5E-1*lam(1&
&,1)*tt91*tt167*tt43*tt66+5.0E-1*miu(1,1)*(tt91*tt167*tt43-tt316*t&
&t65)+2.5E-1*lam(1,1)*tt91*tt167*tt43)
tt318 = (-tt86*tt175)+2*D(1,1)*tt5*tt174-tt171*tt90+2*D(2,1)*tt17&
&*tt89+tt87*tt173+tt172*tt88+((-2*tt85*tt170)+4*D(1,2)*D(2,3)*tt20&
&*tt18+4*D(1,3)*D(2,2)*tt6*tt25)*tt34-tt26*((-tt87*tt170)-tt85*tt1&
&72+2*D(1,3)*tt6*tt171+2*D(2,3)*tt86*tt18)+tt19*(tt86*tt170-2*D(1,&
&2)*tt20*tt172+tt85*tt171-2*D(2,2)*tt87*tt25)
tt319 = volume(1,1)*(2.5E-1*lam(1,1)*tt318*tt65*tt66-2.5E-1*lam(1&
&,1)*tt91*tt176*tt43*tt66+5.0E-1*miu(1,1)*(tt91*tt176*tt43-tt318*t&
&t65)+2.5E-1*lam(1,1)*tt91*tt176*tt43)
tt320 = 2*D(1,1)*D(3,1)
tt321 = 2*D(1,2)*D(3,2)
tt322 = 2*D(1,3)*D(3,3)
tt323 = D(1,2)*D(3,3)+D(1,3)*D(3,2)
tt324 = tt323*tt26
tt325 = -2*D(1,2)*D(3,2)*tt19
tt326 = D(1,1)*D(3,2)+D(1,2)*D(3,1)
tt327 = tt326*tt29
tt328 = D(1,1)*D(3,3)+D(1,3)*D(3,1)
tt329 = -tt328*tt31
tt330 = tt328*tt33
tt331 = -2*tt323*tt29
tt332 = 2*D(1,3)*D(3,3)*tt31
tt333 = 2*D(1,2)*D(3,2)*tt35
tt334 = 2*D(1,3)*D(3,3)*tt26
tt335 = -tt323*tt19
tt336 = -tt328*tt29
tt337 = tt326*tt35
tt338 = -tt326*tt38
tt339 = 2*D(1,1)*D(3,1)*tt39
tt340 = tt339+tt338-tt86*tt184-tt180*tt90+2*D(1,1)*tt5*tt183+2*D(&
&3,1)*tt5*tt89-tt26*(tt337+tt336+tt335+tt334-tt87*tt179-tt85*tt181&
&+2*tt86*D(3,3)*tt6+2*D(1,3)*tt180*tt6)+tt34*(tt333+tt332+tt331-2*&
&tt85*tt179+4*D(1,2)*tt20*D(3,3)*tt6+4*D(1,3)*D(3,2)*tt20*tt6)+tt3&
&30+tt87*tt182+tt181*tt88+tt19*(tt329+tt327+tt325+tt324+tt86*tt179&
&-2*D(1,2)*tt20*tt181+tt180*tt85-2*D(3,2)*tt20*tt87)
tt341 = volume(1,1)*(2.5E-1*lam(1,1)*tt340*tt65*tt66-2.5E-1*lam(1&
&,1)*tt91*tt185*tt43*tt66+5.0E-1*miu(1,1)*((-tt340*tt65)+tt91*tt18&
&5*tt43+tt322+tt321+tt320)+2.5E-1*lam(1,1)*tt91*tt185*tt43)
tt342 = (-tt86*tt213)-tt209*tt90+2*D(1,1)*tt5*tt212+2*D(3,1)*tt11&
&*tt89+tt87*tt211+tt210*tt88+((-2*tt85*tt208)+4*D(1,2)*tt20*D(3,3)&
&*tt12+4*D(1,3)*D(3,2)*tt24*tt6)*tt34+(tt86*tt208-2*D(1,2)*tt20*tt&
&210+tt209*tt85-2*D(3,2)*tt24*tt87)*tt19-((-tt87*tt208)-tt85*tt210&
&+2*tt86*D(3,3)*tt12+2*D(1,3)*tt209*tt6)*tt26
tt343 = volume(1,1)*(2.5E-1*lam(1,1)*tt342*tt65*tt66-2.5E-1*lam(1&
&,1)*tt91*tt214*tt43*tt66+5.0E-1*miu(1,1)*(tt91*tt214*tt43-tt342*t&
&t65)+2.5E-1*lam(1,1)*tt91*tt214*tt43)
tt344 = (-tt86*tt222)+2*D(1,1)*tt5*tt221-tt218*tt90+2*D(3,1)*tt17&
&*tt89+tt87*tt220+tt219*tt88+((-2*tt85*tt217)+4*D(1,2)*tt20*D(3,3)&
&*tt18+4*D(1,3)*D(3,2)*tt6*tt25)*tt34-tt26*((-tt87*tt217)-tt85*tt2&
&19+2*D(1,3)*tt6*tt218+2*tt86*D(3,3)*tt18)+tt19*(tt86*tt217-2*D(1,&
&2)*tt20*tt219+tt85*tt218-2*D(3,2)*tt87*tt25)
tt345 = volume(1,1)*(2.5E-1*lam(1,1)*tt344*tt65*tt66-2.5E-1*lam(1&
&,1)*tt91*tt223*tt43*tt66+5.0E-1*miu(1,1)*(tt91*tt223*tt43-tt344*t&
&t65)+2.5E-1*lam(1,1)*tt91*tt223*tt43)
tt346 = tt120**2
tt347 = tt288+tt287-2*tt115*tt119+4*D(1,1)*tt11*tt118+tt34*(tt286&
&+tt285+tt284-2*tt114**2+8*D(1,2)*D(1,3)*tt24*tt12)-tt26*(tt283+tt&
&282+tt281+tt280-2*tt116*tt114+4*D(1,3)*tt115*tt12)+tt279+2*tt116*&
&tt117+tt19*(tt278+tt277+tt276+tt275+2*tt115*tt114-4*D(1,2)*tt24*t&
&t116)
tt348 = (-tt115*tt128)+2*D(1,1)*tt11*tt127-tt124*tt119+2*D(1,1)*t&
&t17*tt118+tt116*tt126+tt125*tt117+((-2*tt114*tt123)+4*D(1,2)*D(1,&
&3)*tt24*tt18+4*D(1,2)*D(1,3)*tt12*tt25)*tt34-tt26*((-tt116*tt123)&
&-tt114*tt125+2*D(1,3)*tt12*tt124+2*D(1,3)*tt115*tt18)+tt19*(tt115&
&*tt123-2*D(1,2)*tt24*tt125+tt114*tt124-2*D(1,2)*tt116*tt25)
tt349 = volume(1,1)*(2.5E-1*lam(1,1)*tt348*tt65*tt66-2.5E-1*lam(1&
&,1)*tt120*tt129*tt43*tt66+5.0E-1*miu(1,1)*(tt120*tt129*tt43-tt348&
&*tt65)+2.5E-1*lam(1,1)*tt120*tt129*tt43)
tt350 = (-tt133*tt119)-tt115*tt137+2*D(2,1)*tt5*tt118+2*D(1,1)*tt&
&11*tt136+tt134*tt117+tt116*tt135+((-2*tt132*tt114)+4*D(1,3)*D(2,2&
&)*tt20*tt12+4*D(1,2)*D(2,3)*tt24*tt6)*tt34+(tt133*tt114-2*D(2,2)*&
&tt20*tt116+tt115*tt132-2*D(1,2)*tt24*tt134)*tt19-((-tt134*tt114)-&
&tt132*tt116+2*D(1,3)*tt133*tt12+2*D(2,3)*tt115*tt6)*tt26
tt351 = volume(1,1)*(2.5E-1*lam(1,1)*tt350*tt65*tt66-2.5E-1*lam(1&
&,1)*tt138*tt120*tt43*tt66+5.0E-1*miu(1,1)*(tt138*tt120*tt43-tt350&
&*tt65)+2.5E-1*lam(1,1)*tt138*tt120*tt43)
tt352 = tt313+tt312-tt115*tt166-tt162*tt119+2*D(1,1)*tt11*tt165+2&
&*D(2,1)*tt11*tt118-tt26*(tt311+tt310+tt309+tt308-tt116*tt161-tt11&
&4*tt163+2*D(1,3)*tt162*tt12+2*D(2,3)*tt115*tt12)+tt34*(tt307+tt30&
&6+tt305-2*tt114*tt161+4*D(1,2)*D(2,3)*tt24*tt12+4*D(1,3)*D(2,2)*t&
&t24*tt12)+tt304+tt116*tt164+tt163*tt117+tt19*(tt303+tt301+tt299+t&
&t298+tt115*tt161-2*D(1,2)*tt24*tt163+tt162*tt114-2*D(2,2)*tt24*tt&
&116)
tt353 = volume(1,1)*(2.5E-1*lam(1,1)*tt352*tt65*tt66-2.5E-1*lam(1&
&,1)*tt120*tt167*tt43*tt66+5.0E-1*miu(1,1)*((-tt352*tt65)+tt120*tt&
&167*tt43+tt296+tt295+tt294)+2.5E-1*lam(1,1)*tt120*tt167*tt43)
tt354 = (-tt115*tt175)+2*D(1,1)*tt11*tt174-tt171*tt119+2*D(2,1)*t&
&t17*tt118+tt116*tt173+tt172*tt117+((-2*tt114*tt170)+4*D(1,2)*D(2,&
&3)*tt24*tt18+4*D(1,3)*D(2,2)*tt12*tt25)*tt34-tt26*((-tt116*tt170)&
&-tt114*tt172+2*D(1,3)*tt12*tt171+2*D(2,3)*tt115*tt18)+tt19*(tt115&
&*tt170-2*D(1,2)*tt24*tt172+tt114*tt171-2*D(2,2)*tt116*tt25)
tt355 = volume(1,1)*(2.5E-1*lam(1,1)*tt354*tt65*tt66-2.5E-1*lam(1&
&,1)*tt120*tt176*tt43*tt66+5.0E-1*miu(1,1)*(tt120*tt176*tt43-tt354&
&*tt65)+2.5E-1*lam(1,1)*tt120*tt176*tt43)
tt356 = (-tt180*tt119)-tt115*tt184+2*D(3,1)*tt5*tt118+2*D(1,1)*tt&
&11*tt183+tt181*tt117+tt116*tt182+((-2*tt179*tt114)+4*D(1,3)*D(3,2&
&)*tt20*tt12+4*D(1,2)*tt24*D(3,3)*tt6)*tt34+(tt180*tt114-2*D(3,2)*&
&tt20*tt116+tt115*tt179-2*D(1,2)*tt24*tt181)*tt19-((-tt181*tt114)-&
&tt179*tt116+2*D(1,3)*tt180*tt12+2*tt115*D(3,3)*tt6)*tt26
tt357 = volume(1,1)*(2.5E-1*lam(1,1)*tt356*tt65*tt66-2.5E-1*lam(1&
&,1)*tt185*tt120*tt43*tt66+5.0E-1*miu(1,1)*(tt185*tt120*tt43-tt356&
&*tt65)+2.5E-1*lam(1,1)*tt185*tt120*tt43)
tt358 = tt339+tt338-tt115*tt213-tt209*tt119+2*D(1,1)*tt11*tt212+2&
&*D(3,1)*tt11*tt118-tt26*(tt337+tt336+tt335+tt334-tt116*tt208-tt11&
&4*tt210+2*tt115*D(3,3)*tt12+2*D(1,3)*tt209*tt12)+tt34*(tt333+tt33&
&2+tt331-2*tt114*tt208+4*D(1,2)*tt24*D(3,3)*tt12+4*D(1,3)*D(3,2)*t&
&t24*tt12)+tt330+tt116*tt211+tt210*tt117+tt19*(tt329+tt327+tt325+t&
&t324+tt115*tt208-2*D(1,2)*tt24*tt210+tt209*tt114-2*D(3,2)*tt24*tt&
&116)
tt359 = volume(1,1)*(2.5E-1*lam(1,1)*tt358*tt65*tt66-2.5E-1*lam(1&
&,1)*tt120*tt214*tt43*tt66+5.0E-1*miu(1,1)*((-tt358*tt65)+tt120*tt&
&214*tt43+tt322+tt321+tt320)+2.5E-1*lam(1,1)*tt120*tt214*tt43)
tt360 = (-tt115*tt222)+2*D(1,1)*tt11*tt221-tt218*tt119+2*D(3,1)*t&
&t17*tt118+tt116*tt220+tt219*tt117+((-2*tt114*tt217)+4*D(1,2)*tt24&
&*D(3,3)*tt18+4*D(1,3)*D(3,2)*tt12*tt25)*tt34-tt26*((-tt116*tt217)&
&-tt114*tt219+2*D(1,3)*tt12*tt218+2*tt115*D(3,3)*tt18)+tt19*(tt115&
&*tt217-2*D(1,2)*tt24*tt219+tt114*tt218-2*D(3,2)*tt116*tt25)
tt361 = volume(1,1)*(2.5E-1*lam(1,1)*tt360*tt65*tt66-2.5E-1*lam(1&
&,1)*tt120*tt223*tt43*tt66+5.0E-1*miu(1,1)*(tt120*tt223*tt43-tt360&
&*tt65)+2.5E-1*lam(1,1)*tt120*tt223*tt43)
tt362 = tt129**2
tt363 = tt288+tt34*((-2*tt123**2)+tt286+tt285+tt284+8*D(1,2)*D(1,&
&3)*tt25*tt18)+tt287-2*tt124*tt128+4*D(1,1)*tt17*tt127-tt26*(tt283&
&+tt282+tt281-2*tt125*tt123+tt280+4*D(1,3)*tt18*tt124)+tt279+2*tt1&
&25*tt126+tt19*(tt278+tt277+tt276+2*tt124*tt123-4*D(1,2)*tt25*tt12&
&5+tt275)
tt364 = (-tt133*tt128)+2*D(2,1)*tt5*tt127-tt124*tt137+2*D(1,1)*tt&
&17*tt136+tt134*tt126+tt125*tt135+((-2*tt132*tt123)+4*D(1,3)*D(2,2&
&)*tt20*tt18+4*D(1,2)*D(2,3)*tt6*tt25)*tt34-tt26*((-tt134*tt123)-t&
&t132*tt125+2*D(2,3)*tt6*tt124+2*D(1,3)*tt133*tt18)+tt19*(tt133*tt&
&123-2*D(2,2)*tt20*tt125+tt132*tt124-2*D(1,2)*tt134*tt25)
tt365 = volume(1,1)*(2.5E-1*lam(1,1)*tt364*tt65*tt66-2.5E-1*lam(1&
&,1)*tt138*tt129*tt43*tt66+5.0E-1*miu(1,1)*(tt138*tt129*tt43-tt364&
&*tt65)+2.5E-1*lam(1,1)*tt138*tt129*tt43)
tt366 = (-tt162*tt128)+2*D(2,1)*tt11*tt127-tt124*tt166+2*D(1,1)*t&
&t17*tt165+tt163*tt126+tt125*tt164+((-2*tt161*tt123)+4*D(1,3)*D(2,&
&2)*tt24*tt18+4*D(1,2)*D(2,3)*tt12*tt25)*tt34-tt26*((-tt163*tt123)&
&-tt161*tt125+2*D(2,3)*tt12*tt124+2*D(1,3)*tt162*tt18)+tt19*(tt162&
&*tt123-2*D(2,2)*tt24*tt125+tt161*tt124-2*D(1,2)*tt163*tt25)
tt367 = volume(1,1)*(2.5E-1*lam(1,1)*tt366*tt65*tt66-2.5E-1*lam(1&
&,1)*tt167*tt129*tt43*tt66+5.0E-1*miu(1,1)*(tt167*tt129*tt43-tt366&
&*tt65)+2.5E-1*lam(1,1)*tt167*tt129*tt43)
tt368 = tt313+tt312-tt124*tt175-tt171*tt128+2*D(1,1)*tt17*tt174+2&
&*D(2,1)*tt17*tt127-tt26*(tt311+tt310+tt309-tt125*tt170-tt123*tt17&
&2+tt308+2*D(1,3)*tt18*tt171+2*D(2,3)*tt18*tt124)+tt34*(tt307+tt30&
&6+tt305-2*tt123*tt170+4*D(1,2)*D(2,3)*tt25*tt18+4*D(1,3)*D(2,2)*t&
&t25*tt18)+tt304+tt125*tt173+tt172*tt126+tt19*(tt303+tt301+tt299+t&
&t124*tt170-2*D(1,2)*tt25*tt172+tt171*tt123-2*D(2,2)*tt25*tt125+tt&
&298)
tt369 = volume(1,1)*(2.5E-1*lam(1,1)*tt368*tt65*tt66-2.5E-1*lam(1&
&,1)*tt129*tt176*tt43*tt66+5.0E-1*miu(1,1)*((-tt368*tt65)+tt129*tt&
&176*tt43+tt296+tt295+tt294)+2.5E-1*lam(1,1)*tt129*tt176*tt43)
tt370 = (-tt180*tt128)+2*D(3,1)*tt5*tt127-tt124*tt184+2*D(1,1)*tt&
&17*tt183+tt181*tt126+tt125*tt182+((-2*tt179*tt123)+4*D(1,3)*D(3,2&
&)*tt20*tt18+4*D(1,2)*D(3,3)*tt6*tt25)*tt34-tt26*((-tt181*tt123)-t&
&t179*tt125+2*D(3,3)*tt6*tt124+2*D(1,3)*tt180*tt18)+tt19*(tt180*tt&
&123-2*D(3,2)*tt20*tt125+tt179*tt124-2*D(1,2)*tt181*tt25)
tt371 = volume(1,1)*(2.5E-1*lam(1,1)*tt370*tt65*tt66-2.5E-1*lam(1&
&,1)*tt185*tt129*tt43*tt66+5.0E-1*miu(1,1)*(tt185*tt129*tt43-tt370&
&*tt65)+2.5E-1*lam(1,1)*tt185*tt129*tt43)
tt372 = (-tt209*tt128)+2*D(3,1)*tt11*tt127-tt124*tt213+2*D(1,1)*t&
&t17*tt212+tt210*tt126+tt125*tt211+((-2*tt208*tt123)+4*D(1,3)*D(3,&
&2)*tt24*tt18+4*D(1,2)*D(3,3)*tt12*tt25)*tt34-tt26*((-tt210*tt123)&
&-tt208*tt125+2*D(3,3)*tt12*tt124+2*D(1,3)*tt209*tt18)+tt19*(tt209&
&*tt123-2*D(3,2)*tt24*tt125+tt208*tt124-2*D(1,2)*tt210*tt25)
tt373 = volume(1,1)*(2.5E-1*lam(1,1)*tt372*tt65*tt66-2.5E-1*lam(1&
&,1)*tt214*tt129*tt43*tt66+5.0E-1*miu(1,1)*(tt214*tt129*tt43-tt372&
&*tt65)+2.5E-1*lam(1,1)*tt214*tt129*tt43)
tt374 = tt339+tt338-tt124*tt222-tt218*tt128+2*D(1,1)*tt17*tt221+2&
&*D(3,1)*tt17*tt127-tt26*(tt337+tt336+tt335-tt125*tt217-tt123*tt21&
&9+tt334+2*D(1,3)*tt18*tt218+2*D(3,3)*tt18*tt124)+tt34*(tt333+tt33&
&2+tt331-2*tt123*tt217+4*D(1,2)*D(3,3)*tt25*tt18+4*D(1,3)*D(3,2)*t&
&t25*tt18)+tt330+tt125*tt220+tt219*tt126+tt19*(tt329+tt327+tt325+t&
&t124*tt217-2*D(1,2)*tt25*tt219+tt218*tt123-2*D(3,2)*tt25*tt125+tt&
&324)
tt375 = volume(1,1)*(2.5E-1*lam(1,1)*tt374*tt65*tt66-2.5E-1*lam(1&
&,1)*tt129*tt223*tt43*tt66+5.0E-1*miu(1,1)*((-tt374*tt65)+tt129*tt&
&223*tt43+tt322+tt321+tt320)+2.5E-1*lam(1,1)*tt129*tt223*tt43)
tt376 = tt138**2
tt377 = D(2,1)**2
tt378 = 2*tt377
tt379 = D(2,2)**2
tt380 = 2*tt379
tt381 = D(2,3)**2
tt382 = 2*tt381
tt383 = 2*D(2,2)*D(2,3)*tt26
tt384 = -2*tt379*tt19
tt385 = 2*D(2,1)*D(2,2)*tt29
tt386 = -2*D(2,1)*D(2,3)*tt31
tt387 = 2*D(2,1)*D(2,3)*tt33
tt388 = 2*tt381*tt26
tt389 = -2*D(2,2)*D(2,3)*tt19
tt390 = -2*D(2,1)*D(2,3)*tt29
tt391 = 2*D(2,1)*D(2,2)*tt35
tt392 = -4*D(2,2)*D(2,3)*tt29
tt393 = 2*tt381*tt31
tt394 = 2*tt379*tt35
tt395 = -2*D(2,1)*D(2,2)*tt38
tt396 = 2*tt377*tt39
tt397 = tt396+tt395-2*tt133*tt137+4*D(2,1)*tt5*tt136+tt34*(tt394+&
&tt393+tt392-2*tt132**2+8*D(2,2)*D(2,3)*tt20*tt6)-tt26*(tt391+tt39&
&0+tt389+tt388-2*tt134*tt132+4*D(2,3)*tt133*tt6)+tt387+2*tt134*tt1&
&35+tt19*(tt386+tt385+tt384+tt383+2*tt133*tt132-4*D(2,2)*tt20*tt13&
&4)
tt398 = (-tt133*tt166)-tt162*tt137+2*D(2,1)*tt5*tt165+2*D(2,1)*tt&
&11*tt136+tt134*tt164+tt163*tt135+((-2*tt132*tt161)+4*D(2,2)*D(2,3&
&)*tt20*tt12+4*D(2,2)*D(2,3)*tt24*tt6)*tt34+(tt133*tt161-2*D(2,2)*&
&tt20*tt163+tt162*tt132-2*D(2,2)*tt24*tt134)*tt19-((-tt134*tt161)-&
&tt132*tt163+2*D(2,3)*tt133*tt12+2*D(2,3)*tt162*tt6)*tt26
tt399 = volume(1,1)*(2.5E-1*lam(1,1)*tt398*tt65*tt66-2.5E-1*lam(1&
&,1)*tt138*tt167*tt43*tt66+5.0E-1*miu(1,1)*(tt138*tt167*tt43-tt398&
&*tt65)+2.5E-1*lam(1,1)*tt138*tt167*tt43)
tt400 = (-tt133*tt175)+2*D(2,1)*tt5*tt174-tt171*tt137+2*D(2,1)*tt&
&17*tt136+tt134*tt173+tt172*tt135+((-2*tt132*tt170)+4*D(2,2)*D(2,3&
&)*tt20*tt18+4*D(2,2)*D(2,3)*tt6*tt25)*tt34-tt26*((-tt134*tt170)-t&
&t132*tt172+2*D(2,3)*tt6*tt171+2*D(2,3)*tt133*tt18)+tt19*(tt133*tt&
&170-2*D(2,2)*tt20*tt172+tt132*tt171-2*D(2,2)*tt134*tt25)
tt401 = volume(1,1)*(2.5E-1*lam(1,1)*tt400*tt65*tt66-2.5E-1*lam(1&
&,1)*tt138*tt176*tt43*tt66+5.0E-1*miu(1,1)*(tt138*tt176*tt43-tt400&
&*tt65)+2.5E-1*lam(1,1)*tt138*tt176*tt43)
tt402 = 2*D(2,1)*D(3,1)
tt403 = 2*D(2,2)*D(3,2)
tt404 = 2*D(2,3)*D(3,3)
tt405 = D(2,2)*D(3,3)+D(2,3)*D(3,2)
tt406 = tt405*tt26
tt407 = -2*D(2,2)*D(3,2)*tt19
tt408 = D(2,1)*D(3,2)+D(2,2)*D(3,1)
tt409 = tt408*tt29
tt410 = D(2,1)*D(3,3)+D(2,3)*D(3,1)
tt411 = -tt410*tt31
tt412 = tt410*tt33
tt413 = -2*tt405*tt29
tt414 = 2*D(2,3)*D(3,3)*tt31
tt415 = 2*D(2,2)*D(3,2)*tt35
tt416 = 2*D(2,3)*D(3,3)*tt26
tt417 = -tt405*tt19
tt418 = -tt410*tt29
tt419 = tt408*tt35
tt420 = -tt408*tt38
tt421 = 2*D(2,1)*D(3,1)*tt39
tt422 = tt421+tt420-tt133*tt184-tt180*tt137-tt26*(tt419+tt418+tt4&
&17+tt416-tt134*tt179-tt132*tt181+2*tt133*D(3,3)*tt6+2*D(2,3)*tt18&
&0*tt6)+2*D(2,1)*tt5*tt183+2*D(3,1)*tt5*tt136+tt34*(tt415+tt414+tt&
&413-2*tt132*tt179+4*D(2,2)*tt20*D(3,3)*tt6+4*D(2,3)*D(3,2)*tt20*t&
&t6)+tt412+tt134*tt182+tt181*tt135+tt19*(tt411+tt409+tt407+tt406+t&
&t133*tt179-2*D(2,2)*tt20*tt181+tt180*tt132-2*D(3,2)*tt20*tt134)
tt423 = volume(1,1)*(2.5E-1*lam(1,1)*tt422*tt65*tt66-2.5E-1*lam(1&
&,1)*tt138*tt185*tt43*tt66+5.0E-1*miu(1,1)*((-tt422*tt65)+tt138*tt&
&185*tt43+tt404+tt403+tt402)+2.5E-1*lam(1,1)*tt138*tt185*tt43)
tt424 = (-tt133*tt213)-tt209*tt137+2*D(2,1)*tt5*tt212+2*D(3,1)*tt&
&11*tt136+tt134*tt211+tt210*tt135+((-2*tt132*tt208)+4*D(2,2)*tt20*&
&D(3,3)*tt12+4*D(2,3)*D(3,2)*tt24*tt6)*tt34+(tt133*tt208-2*D(2,2)*&
&tt20*tt210+tt209*tt132-2*D(3,2)*tt24*tt134)*tt19-((-tt134*tt208)-&
&tt132*tt210+2*tt133*D(3,3)*tt12+2*D(2,3)*tt209*tt6)*tt26
tt425 = volume(1,1)*(2.5E-1*lam(1,1)*tt424*tt65*tt66-2.5E-1*lam(1&
&,1)*tt138*tt214*tt43*tt66+5.0E-1*miu(1,1)*(tt138*tt214*tt43-tt424&
&*tt65)+2.5E-1*lam(1,1)*tt138*tt214*tt43)
tt426 = (-tt133*tt222)+2*D(2,1)*tt5*tt221-tt218*tt137+2*D(3,1)*tt&
&17*tt136+tt134*tt220+tt219*tt135+((-2*tt132*tt217)+4*D(2,2)*tt20*&
&D(3,3)*tt18+4*D(2,3)*D(3,2)*tt6*tt25)*tt34-tt26*((-tt134*tt217)-t&
&t132*tt219+2*D(2,3)*tt6*tt218+2*tt133*D(3,3)*tt18)+tt19*(tt133*tt&
&217-2*D(2,2)*tt20*tt219+tt132*tt218-2*D(3,2)*tt134*tt25)
tt427 = volume(1,1)*(2.5E-1*lam(1,1)*tt426*tt65*tt66-2.5E-1*lam(1&
&,1)*tt138*tt223*tt43*tt66+5.0E-1*miu(1,1)*(tt138*tt223*tt43-tt426&
&*tt65)+2.5E-1*lam(1,1)*tt138*tt223*tt43)
tt428 = tt167**2
tt429 = tt396+tt395-2*tt162*tt166+4*D(2,1)*tt11*tt165+tt34*(tt394&
&+tt393+tt392-2*tt161**2+8*D(2,2)*D(2,3)*tt24*tt12)-tt26*(tt391+tt&
&390+tt389+tt388-2*tt163*tt161+4*D(2,3)*tt162*tt12)+tt387+2*tt163*&
&tt164+tt19*(tt386+tt385+tt384+tt383+2*tt162*tt161-4*D(2,2)*tt24*t&
&t163)
tt430 = (-tt162*tt175)+2*D(2,1)*tt11*tt174-tt171*tt166+2*D(2,1)*t&
&t17*tt165+tt163*tt173+tt172*tt164+((-2*tt161*tt170)+4*D(2,2)*D(2,&
&3)*tt24*tt18+4*D(2,2)*D(2,3)*tt12*tt25)*tt34-tt26*((-tt163*tt170)&
&-tt161*tt172+2*D(2,3)*tt12*tt171+2*D(2,3)*tt162*tt18)+tt19*(tt162&
&*tt170-2*D(2,2)*tt24*tt172+tt161*tt171-2*D(2,2)*tt163*tt25)
tt431 = volume(1,1)*(2.5E-1*lam(1,1)*tt430*tt65*tt66-2.5E-1*lam(1&
&,1)*tt167*tt176*tt43*tt66+5.0E-1*miu(1,1)*(tt167*tt176*tt43-tt430&
&*tt65)+2.5E-1*lam(1,1)*tt167*tt176*tt43)
tt432 = (-tt180*tt166)-tt162*tt184+2*D(3,1)*tt5*tt165+2*D(2,1)*tt&
&11*tt183+tt181*tt164+tt163*tt182+((-2*tt179*tt161)+4*D(2,3)*D(3,2&
&)*tt20*tt12+4*D(2,2)*tt24*D(3,3)*tt6)*tt34+(tt180*tt161-2*D(3,2)*&
&tt20*tt163+tt162*tt179-2*D(2,2)*tt24*tt181)*tt19-((-tt181*tt161)-&
&tt179*tt163+2*D(2,3)*tt180*tt12+2*tt162*D(3,3)*tt6)*tt26
tt433 = volume(1,1)*(2.5E-1*lam(1,1)*tt432*tt65*tt66-2.5E-1*lam(1&
&,1)*tt185*tt167*tt43*tt66+5.0E-1*miu(1,1)*(tt185*tt167*tt43-tt432&
&*tt65)+2.5E-1*lam(1,1)*tt185*tt167*tt43)
tt434 = tt421+tt420-tt162*tt213-tt209*tt166+2*D(2,1)*tt11*tt212+2&
&*D(3,1)*tt11*tt165-tt26*(tt419+tt418+tt417+tt416-tt163*tt208-tt16&
&1*tt210+2*tt162*D(3,3)*tt12+2*D(2,3)*tt209*tt12)+tt34*(tt415+tt41&
&4+tt413-2*tt161*tt208+4*D(2,2)*tt24*D(3,3)*tt12+4*D(2,3)*D(3,2)*t&
&t24*tt12)+tt412+tt163*tt211+tt210*tt164+tt19*(tt411+tt409+tt407+t&
&t406+tt162*tt208-2*D(2,2)*tt24*tt210+tt209*tt161-2*D(3,2)*tt24*tt&
&163)
tt435 = volume(1,1)*(2.5E-1*lam(1,1)*tt434*tt65*tt66-2.5E-1*lam(1&
&,1)*tt167*tt214*tt43*tt66+5.0E-1*miu(1,1)*((-tt434*tt65)+tt167*tt&
&214*tt43+tt404+tt403+tt402)+2.5E-1*lam(1,1)*tt167*tt214*tt43)
tt436 = (-tt162*tt222)+2*D(2,1)*tt11*tt221-tt218*tt166+2*D(3,1)*t&
&t17*tt165+tt163*tt220+tt219*tt164+((-2*tt161*tt217)+4*D(2,2)*tt24&
&*D(3,3)*tt18+4*D(2,3)*D(3,2)*tt12*tt25)*tt34-tt26*((-tt163*tt217)&
&-tt161*tt219+2*D(2,3)*tt12*tt218+2*tt162*D(3,3)*tt18)+tt19*(tt162&
&*tt217-2*D(2,2)*tt24*tt219+tt161*tt218-2*D(3,2)*tt163*tt25)
tt437 = volume(1,1)*(2.5E-1*lam(1,1)*tt436*tt65*tt66-2.5E-1*lam(1&
&,1)*tt167*tt223*tt43*tt66+5.0E-1*miu(1,1)*(tt167*tt223*tt43-tt436&
&*tt65)+2.5E-1*lam(1,1)*tt167*tt223*tt43)
tt438 = tt176**2
tt439 = tt396+tt34*((-2*tt170**2)+tt394+tt393+tt392+8*D(2,2)*D(2,&
&3)*tt25*tt18)+tt395-2*tt171*tt175+4*D(2,1)*tt17*tt174-tt26*(tt391&
&+tt390+tt389-2*tt172*tt170+tt388+4*D(2,3)*tt18*tt171)+tt387+2*tt1&
&72*tt173+tt19*(tt386+tt385+tt384+2*tt171*tt170-4*D(2,2)*tt25*tt17&
&2+tt383)
tt440 = (-tt180*tt175)+2*D(3,1)*tt5*tt174-tt171*tt184+2*D(2,1)*tt&
&17*tt183+tt181*tt173+tt172*tt182+((-2*tt179*tt170)+4*D(2,3)*D(3,2&
&)*tt20*tt18+4*D(2,2)*D(3,3)*tt6*tt25)*tt34-tt26*((-tt181*tt170)-t&
&t179*tt172+2*D(3,3)*tt6*tt171+2*D(2,3)*tt180*tt18)+tt19*(tt180*tt&
&170-2*D(3,2)*tt20*tt172+tt179*tt171-2*D(2,2)*tt181*tt25)
tt441 = volume(1,1)*(2.5E-1*lam(1,1)*tt440*tt65*tt66-2.5E-1*lam(1&
&,1)*tt185*tt176*tt43*tt66+5.0E-1*miu(1,1)*(tt185*tt176*tt43-tt440&
&*tt65)+2.5E-1*lam(1,1)*tt185*tt176*tt43)
tt442 = (-tt209*tt175)+2*D(3,1)*tt11*tt174-tt171*tt213+2*D(2,1)*t&
&t17*tt212+tt210*tt173+tt172*tt211+((-2*tt208*tt170)+4*D(2,3)*D(3,&
&2)*tt24*tt18+4*D(2,2)*D(3,3)*tt12*tt25)*tt34-tt26*((-tt210*tt170)&
&-tt208*tt172+2*D(3,3)*tt12*tt171+2*D(2,3)*tt209*tt18)+tt19*(tt209&
&*tt170-2*D(3,2)*tt24*tt172+tt208*tt171-2*D(2,2)*tt210*tt25)
tt443 = volume(1,1)*(2.5E-1*lam(1,1)*tt442*tt65*tt66-2.5E-1*lam(1&
&,1)*tt214*tt176*tt43*tt66+5.0E-1*miu(1,1)*(tt214*tt176*tt43-tt442&
&*tt65)+2.5E-1*lam(1,1)*tt214*tt176*tt43)
tt444 = tt421+tt420-tt171*tt222-tt218*tt175+2*D(2,1)*tt17*tt221+2&
&*D(3,1)*tt17*tt174-tt26*(tt419+tt418+tt417-tt172*tt217-tt170*tt21&
&9+tt416+2*D(2,3)*tt18*tt218+2*D(3,3)*tt18*tt171)+tt34*(tt415+tt41&
&4+tt413-2*tt170*tt217+4*D(2,2)*D(3,3)*tt25*tt18+4*D(2,3)*D(3,2)*t&
&t25*tt18)+tt412+tt172*tt220+tt219*tt173+tt19*(tt411+tt409+tt407+t&
&t171*tt217-2*D(2,2)*tt25*tt219+tt218*tt170-2*D(3,2)*tt25*tt172+tt&
&406)
tt445 = volume(1,1)*(2.5E-1*lam(1,1)*tt444*tt65*tt66-2.5E-1*lam(1&
&,1)*tt176*tt223*tt43*tt66+5.0E-1*miu(1,1)*((-tt444*tt65)+tt176*tt&
&223*tt43+tt404+tt403+tt402)+2.5E-1*lam(1,1)*tt176*tt223*tt43)
tt446 = tt185**2
tt447 = D(3,1)**2
tt448 = 2*tt447
tt449 = D(3,2)**2
tt450 = 2*tt449
tt451 = D(3,3)**2
tt452 = 2*tt451
tt453 = 2*D(3,2)*D(3,3)*tt26
tt454 = -2*tt449*tt19
tt455 = 2*D(3,1)*D(3,2)*tt29
tt456 = -2*D(3,1)*D(3,3)*tt31
tt457 = 2*D(3,1)*D(3,3)*tt33
tt458 = 2*tt451*tt26
tt459 = -2*D(3,2)*D(3,3)*tt19
tt460 = -2*D(3,1)*D(3,3)*tt29
tt461 = 2*D(3,1)*D(3,2)*tt35
tt462 = -4*D(3,2)*D(3,3)*tt29
tt463 = 2*tt451*tt31
tt464 = 2*tt449*tt35
tt465 = -2*D(3,1)*D(3,2)*tt38
tt466 = 2*tt447*tt39
tt467 = tt466+tt465+tt34*(tt464+tt463+tt462-2*tt179**2+8*D(3,2)*t&
&t20*D(3,3)*tt6)-2*tt180*tt184+4*D(3,1)*tt5*tt183-tt26*(tt461+tt46&
&0+tt459+tt458-2*tt181*tt179+4*tt180*D(3,3)*tt6)+tt457+2*tt181*tt1&
&82+tt19*(tt456+tt455+tt454+tt453+2*tt180*tt179-4*D(3,2)*tt20*tt18&
&1)
tt468 = (-tt180*tt213)-tt209*tt184+2*D(3,1)*tt5*tt212+2*D(3,1)*tt&
&11*tt183+tt181*tt211+tt210*tt182+((-2*tt179*tt208)+4*D(3,2)*tt20*&
&D(3,3)*tt12+4*D(3,2)*tt24*D(3,3)*tt6)*tt34+(tt180*tt208-2*D(3,2)*&
&tt20*tt210+tt209*tt179-2*D(3,2)*tt24*tt181)*tt19-((-tt181*tt208)-&
&tt179*tt210+2*tt180*D(3,3)*tt12+2*tt209*D(3,3)*tt6)*tt26
tt469 = volume(1,1)*(2.5E-1*lam(1,1)*tt468*tt65*tt66-2.5E-1*lam(1&
&,1)*tt185*tt214*tt43*tt66+5.0E-1*miu(1,1)*(tt185*tt214*tt43-tt468&
&*tt65)+2.5E-1*lam(1,1)*tt185*tt214*tt43)
tt470 = (-tt180*tt222)+2*D(3,1)*tt5*tt221-tt218*tt184+2*D(3,1)*tt&
&17*tt183+tt181*tt220+tt219*tt182+((-2*tt179*tt217)+4*D(3,2)*tt20*&
&D(3,3)*tt18+4*D(3,2)*D(3,3)*tt6*tt25)*tt34-tt26*((-tt181*tt217)-t&
&t179*tt219+2*D(3,3)*tt6*tt218+2*tt180*D(3,3)*tt18)+tt19*(tt180*tt&
&217-2*D(3,2)*tt20*tt219+tt179*tt218-2*D(3,2)*tt181*tt25)
tt471 = volume(1,1)*(2.5E-1*lam(1,1)*tt470*tt65*tt66-2.5E-1*lam(1&
&,1)*tt185*tt223*tt43*tt66+5.0E-1*miu(1,1)*(tt185*tt223*tt43-tt470&
&*tt65)+2.5E-1*lam(1,1)*tt185*tt223*tt43)
tt472 = tt214**2
tt473 = tt466+tt465+tt34*(tt464+tt463+tt462-2*tt208**2+8*D(3,2)*t&
&t24*D(3,3)*tt12)-2*tt209*tt213+4*D(3,1)*tt11*tt212-tt26*(tt461+tt&
&460+tt459+tt458-2*tt210*tt208+4*tt209*D(3,3)*tt12)+tt457+2*tt210*&
&tt211+tt19*(tt456+tt455+tt454+tt453+2*tt209*tt208-4*D(3,2)*tt24*t&
&t210)
tt474 = (-tt209*tt222)+2*D(3,1)*tt11*tt221-tt218*tt213+2*D(3,1)*t&
&t17*tt212+tt210*tt220+tt219*tt211+((-2*tt208*tt217)+4*D(3,2)*tt24&
&*D(3,3)*tt18+4*D(3,2)*D(3,3)*tt12*tt25)*tt34-tt26*((-tt210*tt217)&
&-tt208*tt219+2*D(3,3)*tt12*tt218+2*tt209*D(3,3)*tt18)+tt19*(tt209&
&*tt217-2*D(3,2)*tt24*tt219+tt208*tt218-2*D(3,2)*tt210*tt25)
tt475 = volume(1,1)*(2.5E-1*lam(1,1)*tt474*tt65*tt66-2.5E-1*lam(1&
&,1)*tt214*tt223*tt43*tt66+5.0E-1*miu(1,1)*(tt214*tt223*tt43-tt474&
&*tt65)+2.5E-1*lam(1,1)*tt214*tt223*tt43)
tt476 = tt223**2
tt477 = tt466+tt34*((-2*tt217**2)+tt464+tt463+tt462+8*D(3,2)*D(3,&
&3)*tt25*tt18)+tt465-2*tt218*tt222+4*D(3,1)*tt17*tt221-tt26*(tt461&
&+tt460+tt459-2*tt219*tt217+tt458+4*D(3,3)*tt18*tt218)+tt457+2*tt2&
&19*tt220+tt19*(tt456+tt455+tt454+2*tt218*tt217-4*D(3,2)*tt25*tt21&
&9+tt453)
hes(1,1) = volume(1,1)*(2.5E-1*lam(1,1)*tt64*tt65*tt66-2.5E-1*lam&
&(1,1)*tt41*tt43*tt66+5.0E-1*miu(1,1)*((-tt64*tt65)+tt41*tt43+tt49&
&+tt47+tt45)+2.5E-1*lam(1,1)*tt41*tt43)
hes(1,2) = tt75
hes(1,3) = tt84
hes(1,4) = tt113
hes(1,5) = tt122
hes(1,6) = tt131
hes(1,7) = tt160
hes(1,8) = tt169
hes(1,9) = tt178
hes(1,10) = tt207
hes(1,11) = tt216
hes(1,12) = tt225
hes(2,1) = tt75
hes(2,2) = volume(1,1)*(2.5E-1*lam(1,1)*tt227*tt65*tt66-2.5E-1*la&
&m(1,1)*tt226*tt43*tt66+5.0E-1*miu(1,1)*((-tt227*tt65)+tt226*tt43+&
&tt49+tt47+tt45)+2.5E-1*lam(1,1)*tt226*tt43)
hes(2,3) = tt229
hes(2,4) = tt231
hes(2,5) = tt233
hes(2,6) = tt235
hes(2,7) = tt237
hes(2,8) = tt239
hes(2,9) = tt241
hes(2,10) = tt243
hes(2,11) = tt245
hes(2,12) = tt247
hes(3,1) = tt84
hes(3,2) = tt229
hes(3,3) = volume(1,1)*(2.5E-1*lam(1,1)*tt249*tt65*tt66-2.5E-1*la&
&m(1,1)*tt248*tt43*tt66+5.0E-1*miu(1,1)*((-tt249*tt65)+tt248*tt43+&
&tt49+tt47+tt45)+2.5E-1*lam(1,1)*tt248*tt43)
hes(3,4) = tt251
hes(3,5) = tt253
hes(3,6) = tt255
hes(3,7) = tt257
hes(3,8) = tt259
hes(3,9) = tt261
hes(3,10) = tt263
hes(3,11) = tt265
hes(3,12) = tt267
hes(4,1) = tt113
hes(4,2) = tt231
hes(4,3) = tt251
hes(4,4) = volume(1,1)*(2.5E-1*lam(1,1)*tt289*tt65*tt66-2.5E-1*la&
&m(1,1)*tt268*tt43*tt66+5.0E-1*miu(1,1)*((-tt289*tt65)+tt268*tt43+&
&tt274+tt272+tt270)+2.5E-1*lam(1,1)*tt268*tt43)
hes(4,5) = tt291
hes(4,6) = tt293
hes(4,7) = tt315
hes(4,8) = tt317
hes(4,9) = tt319
hes(4,10) = tt341
hes(4,11) = tt343
hes(4,12) = tt345
hes(5,1) = tt122
hes(5,2) = tt233
hes(5,3) = tt253
hes(5,4) = tt291
hes(5,5) = volume(1,1)*(2.5E-1*lam(1,1)*tt347*tt65*tt66-2.5E-1*la&
&m(1,1)*tt346*tt43*tt66+5.0E-1*miu(1,1)*((-tt347*tt65)+tt346*tt43+&
&tt274+tt272+tt270)+2.5E-1*lam(1,1)*tt346*tt43)
hes(5,6) = tt349
hes(5,7) = tt351
hes(5,8) = tt353
hes(5,9) = tt355
hes(5,10) = tt357
hes(5,11) = tt359
hes(5,12) = tt361
hes(6,1) = tt131
hes(6,2) = tt235
hes(6,3) = tt255
hes(6,4) = tt293
hes(6,5) = tt349
hes(6,6) = volume(1,1)*(2.5E-1*lam(1,1)*tt363*tt65*tt66-2.5E-1*la&
&m(1,1)*tt362*tt43*tt66+5.0E-1*miu(1,1)*((-tt363*tt65)+tt362*tt43+&
&tt274+tt272+tt270)+2.5E-1*lam(1,1)*tt362*tt43)
hes(6,7) = tt365
hes(6,8) = tt367
hes(6,9) = tt369
hes(6,10) = tt371
hes(6,11) = tt373
hes(6,12) = tt375
hes(7,1) = tt160
hes(7,2) = tt237
hes(7,3) = tt257
hes(7,4) = tt315
hes(7,5) = tt351
hes(7,6) = tt365
hes(7,7) = volume(1,1)*(2.5E-1*lam(1,1)*tt397*tt65*tt66-2.5E-1*la&
&m(1,1)*tt376*tt43*tt66+5.0E-1*miu(1,1)*((-tt397*tt65)+tt376*tt43+&
&tt382+tt380+tt378)+2.5E-1*lam(1,1)*tt376*tt43)
hes(7,8) = tt399
hes(7,9) = tt401
hes(7,10) = tt423
hes(7,11) = tt425
hes(7,12) = tt427
hes(8,1) = tt169
hes(8,2) = tt239
hes(8,3) = tt259
hes(8,4) = tt317
hes(8,5) = tt353
hes(8,6) = tt367
hes(8,7) = tt399
hes(8,8) = volume(1,1)*(2.5E-1*lam(1,1)*tt429*tt65*tt66-2.5E-1*la&
&m(1,1)*tt428*tt43*tt66+5.0E-1*miu(1,1)*((-tt429*tt65)+tt428*tt43+&
&tt382+tt380+tt378)+2.5E-1*lam(1,1)*tt428*tt43)
hes(8,9) = tt431
hes(8,10) = tt433
hes(8,11) = tt435
hes(8,12) = tt437
hes(9,1) = tt178
hes(9,2) = tt241
hes(9,3) = tt261
hes(9,4) = tt319
hes(9,5) = tt355
hes(9,6) = tt369
hes(9,7) = tt401
hes(9,8) = tt431
hes(9,9) = volume(1,1)*(2.5E-1*lam(1,1)*tt439*tt65*tt66-2.5E-1*la&
&m(1,1)*tt438*tt43*tt66+5.0E-1*miu(1,1)*((-tt439*tt65)+tt438*tt43+&
&tt382+tt380+tt378)+2.5E-1*lam(1,1)*tt438*tt43)
hes(9,10) = tt441
hes(9,11) = tt443
hes(9,12) = tt445
hes(10,1) = tt207
hes(10,2) = tt243
hes(10,3) = tt263
hes(10,4) = tt341
hes(10,5) = tt357
hes(10,6) = tt371
hes(10,7) = tt423
hes(10,8) = tt433
hes(10,9) = tt441
hes(10,10) = volume(1,1)*(2.5E-1*lam(1,1)*tt467*tt65*tt66-2.5E-1*&
&lam(1,1)*tt446*tt43*tt66+5.0E-1*miu(1,1)*((-tt467*tt65)+tt446*tt4&
&3+tt452+tt450+tt448)+2.5E-1*lam(1,1)*tt446*tt43)
hes(10,11) = tt469
hes(10,12) = tt471
hes(11,1) = tt216
hes(11,2) = tt245
hes(11,3) = tt265
hes(11,4) = tt343
hes(11,5) = tt359
hes(11,6) = tt373
hes(11,7) = tt425
hes(11,8) = tt435
hes(11,9) = tt443
hes(11,10) = tt469
hes(11,11) = volume(1,1)*(2.5E-1*lam(1,1)*tt473*tt65*tt66-2.5E-1*&
&lam(1,1)*tt472*tt43*tt66+5.0E-1*miu(1,1)*((-tt473*tt65)+tt472*tt4&
&3+tt452+tt450+tt448)+2.5E-1*lam(1,1)*tt472*tt43)
hes(11,12) = tt475
hes(12,1) = tt225
hes(12,2) = tt247
hes(12,3) = tt267
hes(12,4) = tt345
hes(12,5) = tt361
hes(12,6) = tt375
hes(12,7) = tt427
hes(12,8) = tt437
hes(12,9) = tt445
hes(12,10) = tt471
hes(12,11) = tt475
hes(12,12) = volume(1,1)*(2.5E-1*lam(1,1)*tt477*tt65*tt66-2.5E-1*&
&lam(1,1)*tt476*tt43*tt66+5.0E-1*miu(1,1)*((-tt477*tt65)+tt476*tt4&
&3+tt452+tt450+tt448)+2.5E-1*lam(1,1)*tt476*tt43)
END 
SUBROUTINE tet_arap(val, X, D, R, volume) 
IMPLICIT NONE 
REAL(KIND=8) val(1, 1) 
REAL(KIND=8) X(3, 4) 
REAL(KIND=8) D(3, 3) 
REAL(KIND=8) R(3, 3) 
REAL(KIND=8) volume(1, 1) 
val(1,1) = volume(1,1)*((D(3,3)*X(3,4)+D(2,3)*X(3,3)-R(3,3)-X(3,1&
&)*D(3,3)+D(1,3)*X(3,2)-D(2,3)*X(3,1)-D(1,3)*X(3,1))**2+(D(3,2)*X(&
&3,4)+D(2,2)*X(3,3)+D(1,2)*X(3,2)-R(3,2)-X(3,1)*D(3,2)-D(2,2)*X(3,&
&1)-D(1,2)*X(3,1))**2+(D(3,1)*X(3,4)+D(2,1)*X(3,3)+D(1,1)*X(3,2)-D&
&(3,1)*X(3,1)-D(2,1)*X(3,1)-D(1,1)*X(3,1)-R(3,1))**2+(X(2,4)*D(3,3&
&)-X(2,1)*D(3,3)+D(2,3)*X(2,3)-R(2,3)-X(2,1)*D(2,3)+D(1,3)*X(2,2)-&
&D(1,3)*X(2,1))**2+(X(1,4)*D(3,3)-X(1,1)*D(3,3)+X(1,3)*D(2,3)-X(1,&
&1)*D(2,3)-R(1,3)+X(1,2)*D(1,3)-X(1,1)*D(1,3))**2+(X(2,4)*D(3,2)-X&
&(2,1)*D(3,2)+D(2,2)*X(2,3)+D(1,2)*X(2,2)-R(2,2)-X(2,1)*D(2,2)-D(1&
&,2)*X(2,1))**2+(X(1,4)*D(3,2)-X(1,1)*D(3,2)+X(1,3)*D(2,2)-X(1,1)*&
&D(2,2)+D(1,2)*X(1,2)-R(1,2)-X(1,1)*D(1,2))**2+(X(2,4)*D(3,1)-X(2,&
&1)*D(3,1)+D(2,1)*X(2,3)+D(1,1)*X(2,2)-D(2,1)*X(2,1)-D(1,1)*X(2,1)&
&-R(2,1))**2+(X(1,4)*D(3,1)-X(1,1)*D(3,1)+X(1,3)*D(2,1)-X(1,1)*D(2&
&,1)+D(1,1)*X(1,2)-D(1,1)*X(1,1)-R(1,1))**2)
END 
SUBROUTINE tet_arap_jac(jac, X, D, R, volume) 
IMPLICIT NONE 
REAL(KIND=8) jac(1, 12) 
REAL(KIND=8) X(3, 4) 
REAL(KIND=8) D(3, 3) 
REAL(KIND=8) R(3, 3) 
REAL(KIND=8) volume(1, 1) 
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
tt1 = (-D(3,1))-D(2,1)-D(1,1)
tt2 = X(1,4)*D(3,1)-X(1,1)*D(3,1)+X(1,3)*D(2,1)-X(1,1)*D(2,1)+D(1&
&,1)*X(1,2)-D(1,1)*X(1,1)-R(1,1)
tt3 = (-D(3,2))-D(2,2)-D(1,2)
tt4 = X(1,4)*D(3,2)-X(1,1)*D(3,2)+X(1,3)*D(2,2)-X(1,1)*D(2,2)+D(1&
&,2)*X(1,2)-R(1,2)-X(1,1)*D(1,2)
tt5 = (-D(3,3))-D(2,3)-D(1,3)
tt6 = X(1,4)*D(3,3)-X(1,1)*D(3,3)+X(1,3)*D(2,3)-X(1,1)*D(2,3)-R(1&
&,3)+X(1,2)*D(1,3)-X(1,1)*D(1,3)
tt7 = X(2,4)*D(3,1)-X(2,1)*D(3,1)+D(2,1)*X(2,3)+D(1,1)*X(2,2)-D(2&
&,1)*X(2,1)-D(1,1)*X(2,1)-R(2,1)
tt8 = X(2,4)*D(3,2)-X(2,1)*D(3,2)+D(2,2)*X(2,3)+D(1,2)*X(2,2)-R(2&
&,2)-X(2,1)*D(2,2)-D(1,2)*X(2,1)
tt9 = X(2,4)*D(3,3)-X(2,1)*D(3,3)+D(2,3)*X(2,3)-R(2,3)-X(2,1)*D(2&
&,3)+D(1,3)*X(2,2)-D(1,3)*X(2,1)
tt10 = D(3,1)*X(3,4)+D(2,1)*X(3,3)+D(1,1)*X(3,2)-D(3,1)*X(3,1)-D(&
&2,1)*X(3,1)-D(1,1)*X(3,1)-R(3,1)
tt11 = D(3,2)*X(3,4)+D(2,2)*X(3,3)+D(1,2)*X(3,2)-R(3,2)-X(3,1)*D(&
&3,2)-D(2,2)*X(3,1)-D(1,2)*X(3,1)
tt12 = D(3,3)*X(3,4)+D(2,3)*X(3,3)-R(3,3)-X(3,1)*D(3,3)+D(1,3)*X(&
&3,2)-D(2,3)*X(3,1)-D(1,3)*X(3,1)
jac(1,1) = volume(1,1)*(2*tt5*tt6+2*tt3*tt4+2*tt1*tt2)
jac(1,2) = volume(1,1)*(2*tt5*tt9+2*tt3*tt8+2*tt1*tt7)
jac(1,3) = volume(1,1)*(2*tt5*tt12+2*tt3*tt11+2*tt1*tt10)
jac(1,4) = volume(1,1)*(2*D(1,3)*tt6+2*D(1,2)*tt4+2*D(1,1)*tt2)
jac(1,5) = volume(1,1)*(2*D(1,3)*tt9+2*D(1,2)*tt8+2*D(1,1)*tt7)
jac(1,6) = volume(1,1)*(2*D(1,3)*tt12+2*D(1,2)*tt11+2*D(1,1)*tt10&
&)
jac(1,7) = volume(1,1)*(2*D(2,3)*tt6+2*D(2,2)*tt4+2*D(2,1)*tt2)
jac(1,8) = volume(1,1)*(2*D(2,3)*tt9+2*D(2,2)*tt8+2*D(2,1)*tt7)
jac(1,9) = volume(1,1)*(2*D(2,3)*tt12+2*D(2,2)*tt11+2*D(2,1)*tt10&
&)
jac(1,10) = volume(1,1)*(2*D(3,3)*tt6+2*D(3,2)*tt4+2*D(3,1)*tt2)
jac(1,11) = volume(1,1)*(2*D(3,3)*tt9+2*D(3,2)*tt8+2*D(3,1)*tt7)
jac(1,12) = volume(1,1)*(2*D(3,3)*tt12+2*D(3,2)*tt11+2*D(3,1)*tt1&
&0)
END 
SUBROUTINE tet_arap_hes(hes, X, D, R, volume) 
IMPLICIT NONE 
REAL(KIND=8) hes(12, 12) 
REAL(KIND=8) X(3, 4) 
REAL(KIND=8) D(3, 3) 
REAL(KIND=8) R(3, 3) 
REAL(KIND=8) volume(1, 1) 
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
tt1 = (-D(3,1))-D(2,1)-D(1,1)
tt2 = (-D(3,2))-D(2,2)-D(1,2)
tt3 = (-D(3,3))-D(2,3)-D(1,3)
tt4 = volume(1,1)*(2*tt3**2+2*tt2**2+2*tt1**2)
tt5 = volume(1,1)*(2*D(1,3)*tt3+2*D(1,2)*tt2+2*D(1,1)*tt1)
tt6 = volume(1,1)*(2*D(2,3)*tt3+2*D(2,2)*tt2+2*D(2,1)*tt1)
tt7 = volume(1,1)*(2*tt3*D(3,3)+2*tt2*D(3,2)+2*tt1*D(3,1))
tt8 = volume(1,1)*(2*D(1,3)**2+2*D(1,2)**2+2*D(1,1)**2)
tt9 = volume(1,1)*(2*D(1,3)*D(2,3)+2*D(1,2)*D(2,2)+2*D(1,1)*D(2,1&
&))
tt10 = volume(1,1)*(2*D(1,3)*D(3,3)+2*D(1,2)*D(3,2)+2*D(1,1)*D(3,&
&1))
tt11 = volume(1,1)*(2*D(2,3)**2+2*D(2,2)**2+2*D(2,1)**2)
tt12 = volume(1,1)*(2*D(2,3)*D(3,3)+2*D(2,2)*D(3,2)+2*D(2,1)*D(3,&
&1))
tt13 = volume(1,1)*(2*D(3,3)**2+2*D(3,2)**2+2*D(3,1)**2)
hes(1,1) = tt4
hes(1,2) = 0
hes(1,3) = 0
hes(1,4) = tt5
hes(1,5) = 0
hes(1,6) = 0
hes(1,7) = tt6
hes(1,8) = 0
hes(1,9) = 0
hes(1,10) = tt7
hes(1,11) = 0
hes(1,12) = 0
hes(2,1) = 0
hes(2,2) = tt4
hes(2,3) = 0
hes(2,4) = 0
hes(2,5) = tt5
hes(2,6) = 0
hes(2,7) = 0
hes(2,8) = tt6
hes(2,9) = 0
hes(2,10) = 0
hes(2,11) = tt7
hes(2,12) = 0
hes(3,1) = 0
hes(3,2) = 0
hes(3,3) = tt4
hes(3,4) = 0
hes(3,5) = 0
hes(3,6) = tt5
hes(3,7) = 0
hes(3,8) = 0
hes(3,9) = tt6
hes(3,10) = 0
hes(3,11) = 0
hes(3,12) = tt7
hes(4,1) = tt5
hes(4,2) = 0
hes(4,3) = 0
hes(4,4) = tt8
hes(4,5) = 0
hes(4,6) = 0
hes(4,7) = tt9
hes(4,8) = 0
hes(4,9) = 0
hes(4,10) = tt10
hes(4,11) = 0
hes(4,12) = 0
hes(5,1) = 0
hes(5,2) = tt5
hes(5,3) = 0
hes(5,4) = 0
hes(5,5) = tt8
hes(5,6) = 0
hes(5,7) = 0
hes(5,8) = tt9
hes(5,9) = 0
hes(5,10) = 0
hes(5,11) = tt10
hes(5,12) = 0
hes(6,1) = 0
hes(6,2) = 0
hes(6,3) = tt5
hes(6,4) = 0
hes(6,5) = 0
hes(6,6) = tt8
hes(6,7) = 0
hes(6,8) = 0
hes(6,9) = tt9
hes(6,10) = 0
hes(6,11) = 0
hes(6,12) = tt10
hes(7,1) = tt6
hes(7,2) = 0
hes(7,3) = 0
hes(7,4) = tt9
hes(7,5) = 0
hes(7,6) = 0
hes(7,7) = tt11
hes(7,8) = 0
hes(7,9) = 0
hes(7,10) = tt12
hes(7,11) = 0
hes(7,12) = 0
hes(8,1) = 0
hes(8,2) = tt6
hes(8,3) = 0
hes(8,4) = 0
hes(8,5) = tt9
hes(8,6) = 0
hes(8,7) = 0
hes(8,8) = tt11
hes(8,9) = 0
hes(8,10) = 0
hes(8,11) = tt12
hes(8,12) = 0
hes(9,1) = 0
hes(9,2) = 0
hes(9,3) = tt6
hes(9,4) = 0
hes(9,5) = 0
hes(9,6) = tt9
hes(9,7) = 0
hes(9,8) = 0
hes(9,9) = tt11
hes(9,10) = 0
hes(9,11) = 0
hes(9,12) = tt12
hes(10,1) = tt7
hes(10,2) = 0
hes(10,3) = 0
hes(10,4) = tt10
hes(10,5) = 0
hes(10,6) = 0
hes(10,7) = tt12
hes(10,8) = 0
hes(10,9) = 0
hes(10,10) = tt13
hes(10,11) = 0
hes(10,12) = 0
hes(11,1) = 0
hes(11,2) = tt7
hes(11,3) = 0
hes(11,4) = 0
hes(11,5) = tt10
hes(11,6) = 0
hes(11,7) = 0
hes(11,8) = tt12
hes(11,9) = 0
hes(11,10) = 0
hes(11,11) = tt13
hes(11,12) = 0
hes(12,1) = 0
hes(12,2) = 0
hes(12,3) = tt7
hes(12,4) = 0
hes(12,5) = 0
hes(12,6) = tt10
hes(12,7) = 0
hes(12,8) = 0
hes(12,9) = tt12
hes(12,10) = 0
hes(12,11) = 0
hes(12,12) = tt13
END 
SUBROUTINE hex_linear(val, X, h, lam, miu) 
IMPLICIT NONE 
REAL(KIND=8) val(1, 1) 
REAL(KIND=8) X(3, 8) 
REAL(KIND=8) h(1, 1) 
REAL(KIND=8) lam(1, 1) 
REAL(KIND=8) miu(1, 1) 
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
tt1 = 1/h(1,1)
tt2 = -2.5E-1*X(1,1)*tt1
tt3 = -2.5E-1*tt1*X(1,3)
tt4 = 2.5E-1*tt1*X(1,4)
tt5 = -2.5E-1*tt1*X(1,5)
tt6 = 2.5E-1*tt1*X(1,6)
tt7 = 2.5E-1*tt1*X(1,8)
tt8 = -2.5E-1*tt1*X(1,2)
tt9 = 2.5E-1*tt1*X(1,7)
tt10 = -2.5E-1*tt1*X(2,1)
tt11 = -2.5E-1*tt1*X(2,3)
tt12 = 2.5E-1*tt1*X(2,4)
tt13 = -2.5E-1*tt1*X(2,5)
tt14 = 2.5E-1*tt1*X(2,6)
tt15 = 2.5E-1*tt1*X(2,8)
tt16 = -2.5E-1*tt1*X(2,2)
tt17 = 2.5E-1*tt1*X(2,7)
tt18 = -2.5E-1*tt1*X(3,1)
tt19 = -2.5E-1*tt1*X(3,3)
tt20 = 2.5E-1*tt1*X(3,4)
tt21 = -2.5E-1*tt1*X(3,5)
tt22 = 2.5E-1*tt1*X(3,6)
tt23 = 2.5E-1*tt1*X(3,8)
tt24 = -2.5E-1*tt1*X(3,2)
tt25 = 2.5E-1*tt1*X(3,7)
val(1,1) = h(1,1)**3*(5.0E-1*lam(1,1)*(5.0E-1*(5.0E-1*tt1*X(3,8)+&
&5.0E-1*tt1*X(3,7)+5.0E-1*tt1*X(3,6)+5.0E-1*tt1*X(3,5)-5.0E-1*tt1*&
&X(3,4)-5.0E-1*tt1*X(3,3)-5.0E-1*tt1*X(3,2)-5.0E-1*tt1*X(3,1))+5.0&
&E-1*(5.0E-1*tt1*X(2,8)+5.0E-1*tt1*X(2,7)-5.0E-1*tt1*X(2,6)-5.0E-1&
&*tt1*X(2,5)+5.0E-1*tt1*X(2,4)+5.0E-1*tt1*X(2,3)-5.0E-1*tt1*X(2,2)&
&-5.0E-1*tt1*X(2,1))+5.0E-1*(5.0E-1*tt1*X(1,8)-5.0E-1*tt1*X(1,7)+5&
&.0E-1*tt1*X(1,6)-5.0E-1*tt1*X(1,5)+5.0E-1*tt1*X(1,4)-5.0E-1*tt1*X&
&(1,3)+5.0E-1*tt1*X(1,2)-5.0E-1*X(1,1)*tt1)-3)**2+miu(1,1)*((tt23+&
&tt25+tt22+2.5E-1*tt1*X(3,5)-2.5E-1*tt1*X(3,4)+tt19+tt24+tt18-1)**&
&2+5.0E-1*(tt23+tt25-2.5E-1*tt1*X(3,6)+tt21+tt20+2.5E-1*tt1*X(3,3)&
&+tt24+tt18+tt15+tt17+tt14+2.5E-1*tt1*X(2,5)-2.5E-1*tt1*X(2,4)+tt1&
&1+tt16+tt10)**2+5.0E-1*(tt23-2.5E-1*tt1*X(3,7)+tt22+tt21+tt20+tt1&
&9+2.5E-1*tt1*X(3,2)+tt18+tt7+tt9+tt6+2.5E-1*tt1*X(1,5)-2.5E-1*tt1&
&*X(1,4)+tt3+tt8+tt2)**2+(tt15+tt17-2.5E-1*tt1*X(2,6)+tt13+tt12+2.&
&5E-1*tt1*X(2,3)+tt16+tt10-1)**2+5.0E-1*(tt15-2.5E-1*tt1*X(2,7)+tt&
&14+tt13+tt12+tt11+2.5E-1*tt1*X(2,2)+tt10+tt7+tt9-2.5E-1*tt1*X(1,6&
&)+tt5+tt4+2.5E-1*tt1*X(1,3)+tt8+tt2)**2+(tt7-2.5E-1*tt1*X(1,7)+tt&
&6+tt5+tt4+tt3+2.5E-1*tt1*X(1,2)+tt2-1)**2))
END 
SUBROUTINE hex_linear_jac(jac, X, h, lam, miu) 
IMPLICIT NONE 
REAL(KIND=8) jac(1, 24) 
REAL(KIND=8) X(3, 8) 
REAL(KIND=8) h(1, 1) 
REAL(KIND=8) lam(1, 1) 
REAL(KIND=8) miu(1, 1) 
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
tt1 = h(1,1)**3
tt2 = 1/h(1,1)
tt3 = -2.5E-1*X(1,1)*tt2
tt4 = -2.5E-1*tt2*X(1,3)
tt5 = 2.5E-1*tt2*X(1,4)
tt6 = -2.5E-1*tt2*X(1,5)
tt7 = 2.5E-1*tt2*X(1,6)
tt8 = 2.5E-1*tt2*X(1,8)
tt9 = tt8-2.5E-1*tt2*X(1,7)+tt7+tt6+tt5+tt4+2.5E-1*tt2*X(1,2)+tt3&
&-1
tt10 = -5.0E-1*tt2*tt9
tt11 = -2.5E-1*tt2*X(1,2)
tt12 = 2.5E-1*tt2*X(1,7)
tt13 = -2.5E-1*tt2*X(2,1)
tt14 = -2.5E-1*tt2*X(2,3)
tt15 = 2.5E-1*tt2*X(2,4)
tt16 = -2.5E-1*tt2*X(2,5)
tt17 = 2.5E-1*tt2*X(2,6)
tt18 = 2.5E-1*tt2*X(2,8)
tt19 = tt18-2.5E-1*tt2*X(2,7)+tt17+tt16+tt15+tt14+2.5E-1*tt2*X(2,&
&2)+tt13+tt8+tt12-2.5E-1*tt2*X(1,6)+tt6+tt5+2.5E-1*tt2*X(1,3)+tt11&
&+tt3
tt20 = -2.5E-1*tt2*tt19
tt21 = -2.5E-1*tt2*X(3,1)
tt22 = -2.5E-1*tt2*X(3,3)
tt23 = 2.5E-1*tt2*X(3,4)
tt24 = -2.5E-1*tt2*X(3,5)
tt25 = 2.5E-1*tt2*X(3,6)
tt26 = 2.5E-1*tt2*X(3,8)
tt27 = tt26-2.5E-1*tt2*X(3,7)+tt25+tt24+tt23+tt22+2.5E-1*tt2*X(3,&
&2)+tt21+tt8+tt12+tt7+2.5E-1*tt2*X(1,5)-2.5E-1*tt2*X(1,4)+tt4+tt11&
&+tt3
tt28 = -2.5E-1*tt2*tt27
tt29 = 5.0E-1*(5.0E-1*tt2*X(3,8)+5.0E-1*tt2*X(3,7)+5.0E-1*tt2*X(3&
&,6)+5.0E-1*tt2*X(3,5)-5.0E-1*tt2*X(3,4)-5.0E-1*tt2*X(3,3)-5.0E-1*&
&tt2*X(3,2)-5.0E-1*tt2*X(3,1))+5.0E-1*(5.0E-1*tt2*X(2,8)+5.0E-1*tt&
&2*X(2,7)-5.0E-1*tt2*X(2,6)-5.0E-1*tt2*X(2,5)+5.0E-1*tt2*X(2,4)+5.&
&0E-1*tt2*X(2,3)-5.0E-1*tt2*X(2,2)-5.0E-1*tt2*X(2,1))+5.0E-1*(5.0E&
&-1*tt2*X(1,8)-5.0E-1*tt2*X(1,7)+5.0E-1*tt2*X(1,6)-5.0E-1*tt2*X(1,&
&5)+5.0E-1*tt2*X(1,4)-5.0E-1*tt2*X(1,3)+5.0E-1*tt2*X(1,2)-5.0E-1*X&
&(1,1)*tt2)-3
tt30 = -2.5E-1*tt2*lam(1,1)*tt29
tt31 = -2.5E-1*tt2*X(2,2)
tt32 = 2.5E-1*tt2*X(2,7)
tt33 = tt18+tt32-2.5E-1*tt2*X(2,6)+tt16+tt15+2.5E-1*tt2*X(2,3)+tt&
&31+tt13-1
tt34 = -5.0E-1*tt2*tt33
tt35 = -2.5E-1*tt2*X(3,2)
tt36 = 2.5E-1*tt2*X(3,7)
tt37 = tt26+tt36-2.5E-1*tt2*X(3,6)+tt24+tt23+2.5E-1*tt2*X(3,3)+tt&
&35+tt21+tt18+tt32+tt17+2.5E-1*tt2*X(2,5)-2.5E-1*tt2*X(2,4)+tt14+t&
&t31+tt13
tt38 = -2.5E-1*tt2*tt37
tt39 = tt26+tt36+tt25+2.5E-1*tt2*X(3,5)-2.5E-1*tt2*X(3,4)+tt22+tt&
&35+tt21-1
tt40 = -5.0E-1*tt2*tt39
tt41 = 5.0E-1*tt2*tt9
tt42 = 2.5E-1*tt2*lam(1,1)*tt29
tt43 = 2.5E-1*tt2*tt19
tt44 = 2.5E-1*tt2*tt27
tt45 = 5.0E-1*tt2*tt33
tt46 = 2.5E-1*tt2*tt37
tt47 = 5.0E-1*tt2*tt39
jac(1,1) = tt1*(tt30+miu(1,1)*(tt28+tt20+tt10))
jac(1,2) = tt1*(tt30+miu(1,1)*(tt38+tt34+tt20))
jac(1,3) = tt1*(tt30+miu(1,1)*(tt40+tt38+tt28))
jac(1,4) = tt1*(tt42+miu(1,1)*(tt28+tt20+tt41))
jac(1,5) = tt1*(tt30+miu(1,1)*(tt38+tt34+tt43))
jac(1,6) = tt1*(tt30+miu(1,1)*(tt40+tt38+tt44))
jac(1,7) = tt1*(tt30+miu(1,1)*(tt28+tt43+tt10))
jac(1,8) = tt1*(tt42+miu(1,1)*(tt38+tt45+tt20))
jac(1,9) = tt1*(tt30+miu(1,1)*(tt40+tt46+tt28))
jac(1,10) = tt1*(tt42+miu(1,1)*(tt28+tt43+tt41))
jac(1,11) = tt1*(tt42+miu(1,1)*(tt38+tt45+tt43))
jac(1,12) = tt1*(tt30+miu(1,1)*(tt40+tt46+tt44))
jac(1,13) = tt1*(tt30+miu(1,1)*(tt44+tt20+tt10))
jac(1,14) = tt1*(tt30+miu(1,1)*(tt46+tt34+tt20))
jac(1,15) = tt1*(tt42+miu(1,1)*(tt47+tt38+tt28))
jac(1,16) = tt1*(tt42+miu(1,1)*(tt44+tt20+tt41))
jac(1,17) = tt1*(tt30+miu(1,1)*(tt46+tt34+tt43))
jac(1,18) = tt1*(tt42+miu(1,1)*(tt47+tt38+tt44))
jac(1,19) = tt1*(tt30+miu(1,1)*(tt44+tt43+tt10))
jac(1,20) = tt1*(tt42+miu(1,1)*(tt46+tt45+tt20))
jac(1,21) = tt1*(tt42+miu(1,1)*(tt47+tt46+tt28))
jac(1,22) = tt1*(tt42+miu(1,1)*(tt44+tt43+tt41))
jac(1,23) = tt1*(tt42+miu(1,1)*(tt46+tt45+tt43))
jac(1,24) = tt1*(tt42+miu(1,1)*(tt47+tt46+tt44))
END 
SUBROUTINE hex_linear_hes(hes, X, h, lam, miu) 
IMPLICIT NONE 
REAL(KIND=8) hes(24, 24) 
REAL(KIND=8) X(3, 8) 
REAL(KIND=8) h(1, 1) 
REAL(KIND=8) lam(1, 1) 
REAL(KIND=8) miu(1, 1) 
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
tt1 = h(1,1)**3
tt2 = 1/h(1,1)**2
tt3 = 6.25E-2*tt2*lam(1,1)
tt4 = tt1*(2.5E-1*tt2*miu(1,1)+tt3)
tt5 = 6.25E-2*tt2*miu(1,1)
tt6 = tt1*(tt5+tt3)
tt7 = -6.25E-2*h(1,1)*lam(1,1)
tt8 = -6.25E-2*tt2*miu(1,1)
tt9 = tt1*(tt8+tt3)
tt10 = tt1*(1.25E-1*tt2*miu(1,1)+tt3)
tt11 = -6.25E-2*tt2*lam(1,1)
tt12 = tt1*(tt5+tt11)
tt13 = tt1*(tt11-1.25E-1*tt2*miu(1,1))
tt14 = tt1*(tt8+tt11)
tt15 = 6.25E-2*h(1,1)*lam(1,1)
tt16 = tt1*(tt11-2.5E-1*tt2*miu(1,1))
hes(1,1) = tt4
hes(1,2) = tt6
hes(1,3) = tt6
hes(1,4) = tt7
hes(1,5) = tt9
hes(1,6) = tt9
hes(1,7) = tt10
hes(1,8) = tt12
hes(1,9) = tt6
hes(1,10) = tt13
hes(1,11) = tt14
hes(1,12) = tt9
hes(1,13) = tt10
hes(1,14) = tt6
hes(1,15) = tt12
hes(1,16) = tt13
hes(1,17) = tt9
hes(1,18) = tt14
hes(1,19) = tt15
hes(1,20) = tt12
hes(1,21) = tt12
hes(1,22) = tt16
hes(1,23) = tt14
hes(1,24) = tt14
hes(2,1) = tt6
hes(2,2) = tt4
hes(2,3) = tt6
hes(2,4) = tt12
hes(2,5) = tt10
hes(2,6) = tt6
hes(2,7) = tt9
hes(2,8) = tt7
hes(2,9) = tt9
hes(2,10) = tt14
hes(2,11) = tt13
hes(2,12) = tt9
hes(2,13) = tt6
hes(2,14) = tt10
hes(2,15) = tt12
hes(2,16) = tt12
hes(2,17) = tt15
hes(2,18) = tt12
hes(2,19) = tt9
hes(2,20) = tt13
hes(2,21) = tt14
hes(2,22) = tt14
hes(2,23) = tt16
hes(2,24) = tt14
hes(3,1) = tt6
hes(3,2) = tt6
hes(3,3) = tt4
hes(3,4) = tt12
hes(3,5) = tt6
hes(3,6) = tt10
hes(3,7) = tt6
hes(3,8) = tt12
hes(3,9) = tt10
hes(3,10) = tt12
hes(3,11) = tt12
hes(3,12) = tt15
hes(3,13) = tt9
hes(3,14) = tt9
hes(3,15) = tt7
hes(3,16) = tt14
hes(3,17) = tt9
hes(3,18) = tt13
hes(3,19) = tt9
hes(3,20) = tt14
hes(3,21) = tt13
hes(3,22) = tt14
hes(3,23) = tt14
hes(3,24) = tt16
hes(4,1) = tt7
hes(4,2) = tt12
hes(4,3) = tt12
hes(4,4) = tt4
hes(4,5) = tt14
hes(4,6) = tt14
hes(4,7) = tt13
hes(4,8) = tt6
hes(4,9) = tt12
hes(4,10) = tt10
hes(4,11) = tt9
hes(4,12) = tt14
hes(4,13) = tt13
hes(4,14) = tt12
hes(4,15) = tt6
hes(4,16) = tt10
hes(4,17) = tt14
hes(4,18) = tt9
hes(4,19) = tt16
hes(4,20) = tt6
hes(4,21) = tt6
hes(4,22) = tt15
hes(4,23) = tt9
hes(4,24) = tt9
hes(5,1) = tt9
hes(5,2) = tt10
hes(5,3) = tt6
hes(5,4) = tt14
hes(5,5) = tt4
hes(5,6) = tt6
hes(5,7) = tt6
hes(5,8) = tt13
hes(5,9) = tt9
hes(5,10) = tt12
hes(5,11) = tt7
hes(5,12) = tt9
hes(5,13) = tt9
hes(5,14) = tt15
hes(5,15) = tt12
hes(5,16) = tt14
hes(5,17) = tt10
hes(5,18) = tt12
hes(5,19) = tt6
hes(5,20) = tt16
hes(5,21) = tt14
hes(5,22) = tt12
hes(5,23) = tt13
hes(5,24) = tt14
hes(6,1) = tt9
hes(6,2) = tt6
hes(6,3) = tt10
hes(6,4) = tt14
hes(6,5) = tt6
hes(6,6) = tt4
hes(6,7) = tt9
hes(6,8) = tt12
hes(6,9) = tt15
hes(6,10) = tt14
hes(6,11) = tt12
hes(6,12) = tt10
hes(6,13) = tt6
hes(6,14) = tt9
hes(6,15) = tt13
hes(6,16) = tt12
hes(6,17) = tt9
hes(6,18) = tt7
hes(6,19) = tt6
hes(6,20) = tt14
hes(6,21) = tt16
hes(6,22) = tt12
hes(6,23) = tt14
hes(6,24) = tt13
hes(7,1) = tt10
hes(7,2) = tt9
hes(7,3) = tt6
hes(7,4) = tt13
hes(7,5) = tt6
hes(7,6) = tt9
hes(7,7) = tt4
hes(7,8) = tt14
hes(7,9) = tt6
hes(7,10) = tt7
hes(7,11) = tt12
hes(7,12) = tt9
hes(7,13) = tt15
hes(7,14) = tt9
hes(7,15) = tt12
hes(7,16) = tt16
hes(7,17) = tt6
hes(7,18) = tt14
hes(7,19) = tt10
hes(7,20) = tt14
hes(7,21) = tt12
hes(7,22) = tt13
hes(7,23) = tt12
hes(7,24) = tt14
hes(8,1) = tt12
hes(8,2) = tt7
hes(8,3) = tt12
hes(8,4) = tt6
hes(8,5) = tt13
hes(8,6) = tt12
hes(8,7) = tt14
hes(8,8) = tt4
hes(8,9) = tt14
hes(8,10) = tt9
hes(8,11) = tt10
hes(8,12) = tt14
hes(8,13) = tt12
hes(8,14) = tt13
hes(8,15) = tt6
hes(8,16) = tt6
hes(8,17) = tt16
hes(8,18) = tt6
hes(8,19) = tt14
hes(8,20) = tt10
hes(8,21) = tt9
hes(8,22) = tt9
hes(8,23) = tt15
hes(8,24) = tt9
hes(9,1) = tt6
hes(9,2) = tt9
hes(9,3) = tt10
hes(9,4) = tt12
hes(9,5) = tt9
hes(9,6) = tt15
hes(9,7) = tt6
hes(9,8) = tt14
hes(9,9) = tt4
hes(9,10) = tt12
hes(9,11) = tt14
hes(9,12) = tt10
hes(9,13) = tt9
hes(9,14) = tt6
hes(9,15) = tt13
hes(9,16) = tt14
hes(9,17) = tt6
hes(9,18) = tt16
hes(9,19) = tt9
hes(9,20) = tt12
hes(9,21) = tt7
hes(9,22) = tt14
hes(9,23) = tt12
hes(9,24) = tt13
hes(10,1) = tt13
hes(10,2) = tt14
hes(10,3) = tt12
hes(10,4) = tt10
hes(10,5) = tt12
hes(10,6) = tt14
hes(10,7) = tt7
hes(10,8) = tt9
hes(10,9) = tt12
hes(10,10) = tt4
hes(10,11) = tt6
hes(10,12) = tt14
hes(10,13) = tt16
hes(10,14) = tt14
hes(10,15) = tt6
hes(10,16) = tt15
hes(10,17) = tt12
hes(10,18) = tt9
hes(10,19) = tt13
hes(10,20) = tt9
hes(10,21) = tt6
hes(10,22) = tt10
hes(10,23) = tt6
hes(10,24) = tt9
hes(11,1) = tt14
hes(11,2) = tt13
hes(11,3) = tt12
hes(11,4) = tt9
hes(11,5) = tt7
hes(11,6) = tt12
hes(11,7) = tt12
hes(11,8) = tt10
hes(11,9) = tt14
hes(11,10) = tt6
hes(11,11) = tt4
hes(11,12) = tt14
hes(11,13) = tt14
hes(11,14) = tt16
hes(11,15) = tt6
hes(11,16) = tt9
hes(11,17) = tt13
hes(11,18) = tt6
hes(11,19) = tt12
hes(11,20) = tt15
hes(11,21) = tt9
hes(11,22) = tt6
hes(11,23) = tt10
hes(11,24) = tt9
hes(12,1) = tt9
hes(12,2) = tt9
hes(12,3) = tt15
hes(12,4) = tt14
hes(12,5) = tt9
hes(12,6) = tt10
hes(12,7) = tt9
hes(12,8) = tt14
hes(12,9) = tt10
hes(12,10) = tt14
hes(12,11) = tt14
hes(12,12) = tt4
hes(12,13) = tt6
hes(12,14) = tt6
hes(12,15) = tt16
hes(12,16) = tt12
hes(12,17) = tt6
hes(12,18) = tt13
hes(12,19) = tt6
hes(12,20) = tt12
hes(12,21) = tt13
hes(12,22) = tt12
hes(12,23) = tt12
hes(12,24) = tt7
hes(13,1) = tt10
hes(13,2) = tt6
hes(13,3) = tt9
hes(13,4) = tt13
hes(13,5) = tt9
hes(13,6) = tt6
hes(13,7) = tt15
hes(13,8) = tt12
hes(13,9) = tt9
hes(13,10) = tt16
hes(13,11) = tt14
hes(13,12) = tt6
hes(13,13) = tt4
hes(13,14) = tt6
hes(13,15) = tt14
hes(13,16) = tt7
hes(13,17) = tt9
hes(13,18) = tt12
hes(13,19) = tt10
hes(13,20) = tt12
hes(13,21) = tt14
hes(13,22) = tt13
hes(13,23) = tt14
hes(13,24) = tt12
hes(14,1) = tt6
hes(14,2) = tt10
hes(14,3) = tt9
hes(14,4) = tt12
hes(14,5) = tt15
hes(14,6) = tt9
hes(14,7) = tt9
hes(14,8) = tt13
hes(14,9) = tt6
hes(14,10) = tt14
hes(14,11) = tt16
hes(14,12) = tt6
hes(14,13) = tt6
hes(14,14) = tt4
hes(14,15) = tt14
hes(14,16) = tt12
hes(14,17) = tt10
hes(14,18) = tt14
hes(14,19) = tt9
hes(14,20) = tt7
hes(14,21) = tt12
hes(14,22) = tt14
hes(14,23) = tt13
hes(14,24) = tt12
hes(15,1) = tt12
hes(15,2) = tt12
hes(15,3) = tt7
hes(15,4) = tt6
hes(15,5) = tt12
hes(15,6) = tt13
hes(15,7) = tt12
hes(15,8) = tt6
hes(15,9) = tt13
hes(15,10) = tt6
hes(15,11) = tt6
hes(15,12) = tt16
hes(15,13) = tt14
hes(15,14) = tt14
hes(15,15) = tt4
hes(15,16) = tt9
hes(15,17) = tt14
hes(15,18) = tt10
hes(15,19) = tt14
hes(15,20) = tt9
hes(15,21) = tt10
hes(15,22) = tt9
hes(15,23) = tt9
hes(15,24) = tt15
hes(16,1) = tt13
hes(16,2) = tt12
hes(16,3) = tt14
hes(16,4) = tt10
hes(16,5) = tt14
hes(16,6) = tt12
hes(16,7) = tt16
hes(16,8) = tt6
hes(16,9) = tt14
hes(16,10) = tt15
hes(16,11) = tt9
hes(16,12) = tt12
hes(16,13) = tt7
hes(16,14) = tt12
hes(16,15) = tt9
hes(16,16) = tt4
hes(16,17) = tt14
hes(16,18) = tt6
hes(16,19) = tt13
hes(16,20) = tt6
hes(16,21) = tt9
hes(16,22) = tt10
hes(16,23) = tt9
hes(16,24) = tt6
hes(17,1) = tt9
hes(17,2) = tt15
hes(17,3) = tt9
hes(17,4) = tt14
hes(17,5) = tt10
hes(17,6) = tt9
hes(17,7) = tt6
hes(17,8) = tt16
hes(17,9) = tt6
hes(17,10) = tt12
hes(17,11) = tt13
hes(17,12) = tt6
hes(17,13) = tt9
hes(17,14) = tt10
hes(17,15) = tt14
hes(17,16) = tt14
hes(17,17) = tt4
hes(17,18) = tt14
hes(17,19) = tt6
hes(17,20) = tt13
hes(17,21) = tt12
hes(17,22) = tt12
hes(17,23) = tt7
hes(17,24) = tt12
hes(18,1) = tt14
hes(18,2) = tt12
hes(18,3) = tt13
hes(18,4) = tt9
hes(18,5) = tt12
hes(18,6) = tt7
hes(18,7) = tt14
hes(18,8) = tt6
hes(18,9) = tt16
hes(18,10) = tt9
hes(18,11) = tt6
hes(18,12) = tt13
hes(18,13) = tt12
hes(18,14) = tt14
hes(18,15) = tt10
hes(18,16) = tt6
hes(18,17) = tt14
hes(18,18) = tt4
hes(18,19) = tt12
hes(18,20) = tt9
hes(18,21) = tt15
hes(18,22) = tt6
hes(18,23) = tt9
hes(18,24) = tt10
hes(19,1) = tt15
hes(19,2) = tt9
hes(19,3) = tt9
hes(19,4) = tt16
hes(19,5) = tt6
hes(19,6) = tt6
hes(19,7) = tt10
hes(19,8) = tt14
hes(19,9) = tt9
hes(19,10) = tt13
hes(19,11) = tt12
hes(19,12) = tt6
hes(19,13) = tt10
hes(19,14) = tt9
hes(19,15) = tt14
hes(19,16) = tt13
hes(19,17) = tt6
hes(19,18) = tt12
hes(19,19) = tt4
hes(19,20) = tt14
hes(19,21) = tt14
hes(19,22) = tt7
hes(19,23) = tt12
hes(19,24) = tt12
hes(20,1) = tt12
hes(20,2) = tt13
hes(20,3) = tt14
hes(20,4) = tt6
hes(20,5) = tt16
hes(20,6) = tt14
hes(20,7) = tt14
hes(20,8) = tt10
hes(20,9) = tt12
hes(20,10) = tt9
hes(20,11) = tt15
hes(20,12) = tt12
hes(20,13) = tt12
hes(20,14) = tt7
hes(20,15) = tt9
hes(20,16) = tt6
hes(20,17) = tt13
hes(20,18) = tt9
hes(20,19) = tt14
hes(20,20) = tt4
hes(20,21) = tt6
hes(20,22) = tt9
hes(20,23) = tt10
hes(20,24) = tt6
hes(21,1) = tt12
hes(21,2) = tt14
hes(21,3) = tt13
hes(21,4) = tt6
hes(21,5) = tt14
hes(21,6) = tt16
hes(21,7) = tt12
hes(21,8) = tt9
hes(21,9) = tt7
hes(21,10) = tt6
hes(21,11) = tt9
hes(21,12) = tt13
hes(21,13) = tt14
hes(21,14) = tt12
hes(21,15) = tt10
hes(21,16) = tt9
hes(21,17) = tt12
hes(21,18) = tt15
hes(21,19) = tt14
hes(21,20) = tt6
hes(21,21) = tt4
hes(21,22) = tt9
hes(21,23) = tt6
hes(21,24) = tt10
hes(22,1) = tt16
hes(22,2) = tt14
hes(22,3) = tt14
hes(22,4) = tt15
hes(22,5) = tt12
hes(22,6) = tt12
hes(22,7) = tt13
hes(22,8) = tt9
hes(22,9) = tt14
hes(22,10) = tt10
hes(22,11) = tt6
hes(22,12) = tt12
hes(22,13) = tt13
hes(22,14) = tt14
hes(22,15) = tt9
hes(22,16) = tt10
hes(22,17) = tt12
hes(22,18) = tt6
hes(22,19) = tt7
hes(22,20) = tt9
hes(22,21) = tt9
hes(22,22) = tt4
hes(22,23) = tt6
hes(22,24) = tt6
hes(23,1) = tt14
hes(23,2) = tt16
hes(23,3) = tt14
hes(23,4) = tt9
hes(23,5) = tt13
hes(23,6) = tt14
hes(23,7) = tt12
hes(23,8) = tt15
hes(23,9) = tt12
hes(23,10) = tt6
hes(23,11) = tt10
hes(23,12) = tt12
hes(23,13) = tt14
hes(23,14) = tt13
hes(23,15) = tt9
hes(23,16) = tt9
hes(23,17) = tt7
hes(23,18) = tt9
hes(23,19) = tt12
hes(23,20) = tt10
hes(23,21) = tt6
hes(23,22) = tt6
hes(23,23) = tt4
hes(23,24) = tt6
hes(24,1) = tt14
hes(24,2) = tt14
hes(24,3) = tt16
hes(24,4) = tt9
hes(24,5) = tt14
hes(24,6) = tt13
hes(24,7) = tt14
hes(24,8) = tt9
hes(24,9) = tt13
hes(24,10) = tt9
hes(24,11) = tt9
hes(24,12) = tt7
hes(24,13) = tt12
hes(24,14) = tt12
hes(24,15) = tt15
hes(24,16) = tt6
hes(24,17) = tt12
hes(24,18) = tt10
hes(24,19) = tt12
hes(24,20) = tt6
hes(24,21) = tt10
hes(24,22) = tt6
hes(24,23) = tt6
hes(24,24) = tt4
END 
SUBROUTINE hex_stvk(val, X, h, lam, miu) 
IMPLICIT NONE 
REAL(KIND=8) val(1, 1) 
REAL(KIND=8) X(3, 8) 
REAL(KIND=8) h(1, 1) 
REAL(KIND=8) lam(1, 1) 
REAL(KIND=8) miu(1, 1) 
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
tt1 = 1/h(1,1)**2
tt2 = 6.25E-2*X(1,1)**2*tt1
tt3 = X(1,2)**2
tt4 = -6.25E-2*tt1*tt3
tt5 = X(1,3)**2
tt6 = -6.25E-2*tt1*tt5
tt7 = X(1,4)**2
tt8 = X(1,5)**2
tt9 = 1.25E-1*tt1*X(1,3)*X(1,6)
tt10 = X(1,6)**2
tt11 = -6.25E-2*tt1*tt10
tt12 = 1.25E-1*tt1*X(1,2)*X(1,7)
tt13 = X(1,7)**2
tt14 = -6.25E-2*tt1*tt13
tt15 = -1.25E-1*X(1,1)*tt1*X(1,8)
tt16 = 6.25E-2*tt1*X(1,8)**2
tt17 = 6.25E-2*tt1*X(2,1)**2
tt18 = X(2,2)**2
tt19 = -6.25E-2*tt1*tt18
tt20 = X(2,3)**2
tt21 = -6.25E-2*tt1*tt20
tt22 = X(2,4)**2
tt23 = X(2,5)**2
tt24 = 1.25E-1*tt1*X(2,3)*X(2,6)
tt25 = X(2,6)**2
tt26 = -6.25E-2*tt1*tt25
tt27 = 1.25E-1*tt1*X(2,2)*X(2,7)
tt28 = X(2,7)**2
tt29 = -6.25E-2*tt1*tt28
tt30 = -1.25E-1*tt1*X(2,1)*X(2,8)
tt31 = 6.25E-2*tt1*X(2,8)**2
tt32 = 6.25E-2*tt1*X(3,1)**2
tt33 = X(3,2)**2
tt34 = -6.25E-2*tt1*tt33
tt35 = X(3,3)**2
tt36 = -6.25E-2*tt1*tt35
tt37 = X(3,4)**2
tt38 = X(3,5)**2
tt39 = 1.25E-1*tt1*X(3,3)*X(3,6)
tt40 = X(3,6)**2
tt41 = -6.25E-2*tt1*tt40
tt42 = 1.25E-1*tt1*X(3,2)*X(3,7)
tt43 = X(3,7)**2
tt44 = -6.25E-2*tt1*tt43
tt45 = -1.25E-1*tt1*X(3,1)*X(3,8)
tt46 = 6.25E-2*tt1*X(3,8)**2
tt47 = -6.25E-2*tt1*tt7
tt48 = 1.25E-1*tt1*X(1,4)*X(1,5)
tt49 = -6.25E-2*tt1*tt8
tt50 = -6.25E-2*tt1*tt22
tt51 = 1.25E-1*tt1*X(2,4)*X(2,5)
tt52 = -6.25E-2*tt1*tt23
tt53 = -6.25E-2*tt1*tt37
tt54 = 1.25E-1*tt1*X(3,4)*X(3,5)
tt55 = -6.25E-2*tt1*tt38
tt56 = 1/h(1,1)
tt57 = -2.5E-1*X(1,1)*tt56
tt58 = -2.5E-1*tt56*X(1,3)
tt59 = 2.5E-1*tt56*X(1,4)
tt60 = -2.5E-1*tt56*X(1,5)
tt61 = 2.5E-1*tt56*X(1,6)
tt62 = 2.5E-1*tt56*X(1,8)
tt63 = -2.5E-1*tt56*X(2,1)
tt64 = -2.5E-1*tt56*X(2,3)
tt65 = 2.5E-1*tt56*X(2,4)
tt66 = -2.5E-1*tt56*X(2,5)
tt67 = 2.5E-1*tt56*X(2,6)
tt68 = 2.5E-1*tt56*X(2,8)
tt69 = -2.5E-1*tt56*X(3,1)
tt70 = -2.5E-1*tt56*X(3,3)
tt71 = 2.5E-1*tt56*X(3,4)
tt72 = -2.5E-1*tt56*X(3,5)
tt73 = 2.5E-1*tt56*X(3,6)
tt74 = 2.5E-1*tt56*X(3,8)
tt75 = (tt74-2.5E-1*tt56*X(3,7)+tt73+tt72+tt71+tt70+2.5E-1*tt56*X&
&(3,2)+tt69)**2+(tt68-2.5E-1*tt56*X(2,7)+tt67+tt66+tt65+tt64+2.5E-&
&1*tt56*X(2,2)+tt63)**2+(tt62-2.5E-1*tt56*X(1,7)+tt61+tt60+tt59+tt&
&58+2.5E-1*tt56*X(1,2)+tt57)**2-1
tt76 = -2.5E-1*tt56*X(1,2)
tt77 = 2.5E-1*tt56*X(1,7)
tt78 = -2.5E-1*tt56*X(2,2)
tt79 = 2.5E-1*tt56*X(2,7)
tt80 = -2.5E-1*tt56*X(3,2)
tt81 = 2.5E-1*tt56*X(3,7)
tt82 = (tt74+tt81-2.5E-1*tt56*X(3,6)+tt72+tt71+2.5E-1*tt56*X(3,3)&
&+tt80+tt69)**2+(tt68+tt79-2.5E-1*tt56*X(2,6)+tt66+tt65+2.5E-1*tt5&
&6*X(2,3)+tt78+tt63)**2+(tt62+tt77-2.5E-1*tt56*X(1,6)+tt60+tt59+2.&
&5E-1*tt56*X(1,3)+tt76+tt57)**2-1
tt83 = (tt74+tt81+tt73+2.5E-1*tt56*X(3,5)-2.5E-1*tt56*X(3,4)+tt70&
&+tt80+tt69)**2+(tt68+tt79+tt67+2.5E-1*tt56*X(2,5)-2.5E-1*tt56*X(2&
&,4)+tt64+tt78+tt63)**2+(tt62+tt77+tt61+2.5E-1*tt56*X(1,5)-2.5E-1*&
&tt56*X(1,4)+tt58+tt76+tt57)**2-1
val(1,1) = h(1,1)**3*(5.0E-1*lam(1,1)*(5.0E-1*tt83+5.0E-1*tt82+5.&
&0E-1*tt75)**2+miu(1,1)*(2.5E-1*tt83**2+2.5E-1*tt82**2+2.5E-1*tt75&
&**2+5.0E-1*(tt46+1.25E-1*tt1*X(3,7)*X(3,8)-1.25E-1*tt1*X(3,2)*X(3&
&,8)+tt45+6.25E-2*tt1*tt43-1.25E-1*tt1*X(3,2)*X(3,7)-1.25E-1*tt1*X&
&(3,1)*X(3,7)+tt41-1.25E-1*tt1*X(3,5)*X(3,6)+1.25E-1*tt1*X(3,4)*X(&
&3,6)+tt39+tt55+tt54+1.25E-1*tt1*X(3,3)*X(3,5)+tt53-1.25E-1*tt1*X(&
&3,3)*X(3,4)+tt36+6.25E-2*tt1*tt33+1.25E-1*tt1*X(3,1)*X(3,2)+tt32+&
&tt31+1.25E-1*tt1*X(2,7)*X(2,8)-1.25E-1*tt1*X(2,2)*X(2,8)+tt30+6.2&
&5E-2*tt1*tt28-1.25E-1*tt1*X(2,2)*X(2,7)-1.25E-1*tt1*X(2,1)*X(2,7)&
&+tt26-1.25E-1*tt1*X(2,5)*X(2,6)+1.25E-1*tt1*X(2,4)*X(2,6)+tt24+tt&
&52+tt51+1.25E-1*tt1*X(2,3)*X(2,5)+tt50-1.25E-1*tt1*X(2,3)*X(2,4)+&
&tt21+6.25E-2*tt1*tt18+1.25E-1*tt1*X(2,1)*X(2,2)+tt17+tt16+1.25E-1&
&*tt1*X(1,7)*X(1,8)-1.25E-1*tt1*X(1,2)*X(1,8)+tt15+6.25E-2*tt1*tt1&
&3-1.25E-1*tt1*X(1,2)*X(1,7)-1.25E-1*X(1,1)*tt1*X(1,7)+tt11-1.25E-&
&1*tt1*X(1,5)*X(1,6)+1.25E-1*tt1*X(1,4)*X(1,6)+tt9+tt49+tt48+1.25E&
&-1*tt1*X(1,3)*X(1,5)+tt47-1.25E-1*tt1*X(1,3)*X(1,4)+tt6+6.25E-2*t&
&t1*tt3+1.25E-1*X(1,1)*tt1*X(1,2)+tt2)**2+5.0E-1*(tt46+1.25E-1*tt1&
&*X(3,6)*X(3,8)-1.25E-1*tt1*X(3,3)*X(3,8)+tt45+tt44-1.25E-1*tt1*X(&
&3,5)*X(3,7)+1.25E-1*tt1*X(3,4)*X(3,7)+tt42+6.25E-2*tt1*tt40-1.25E&
&-1*tt1*X(3,3)*X(3,6)-1.25E-1*tt1*X(3,1)*X(3,6)+tt55+tt54+1.25E-1*&
&tt1*X(3,2)*X(3,5)+tt53-1.25E-1*tt1*X(3,2)*X(3,4)+6.25E-2*tt1*tt35&
&+1.25E-1*tt1*X(3,1)*X(3,3)+tt34+tt32+tt31+1.25E-1*tt1*X(2,6)*X(2,&
&8)-1.25E-1*tt1*X(2,3)*X(2,8)+tt30+tt29-1.25E-1*tt1*X(2,5)*X(2,7)+&
&1.25E-1*tt1*X(2,4)*X(2,7)+tt27+6.25E-2*tt1*tt25-1.25E-1*tt1*X(2,3&
&)*X(2,6)-1.25E-1*tt1*X(2,1)*X(2,6)+tt52+tt51+1.25E-1*tt1*X(2,2)*X&
&(2,5)+tt50-1.25E-1*tt1*X(2,2)*X(2,4)+6.25E-2*tt1*tt20+1.25E-1*tt1&
&*X(2,1)*X(2,3)+tt19+tt17+tt16+1.25E-1*tt1*X(1,6)*X(1,8)-1.25E-1*t&
&t1*X(1,3)*X(1,8)+tt15+tt14-1.25E-1*tt1*X(1,5)*X(1,7)+1.25E-1*tt1*&
&X(1,4)*X(1,7)+tt12+6.25E-2*tt1*tt10-1.25E-1*tt1*X(1,3)*X(1,6)-1.2&
&5E-1*X(1,1)*tt1*X(1,6)+tt49+tt48+1.25E-1*tt1*X(1,2)*X(1,5)+tt47-1&
&.25E-1*tt1*X(1,2)*X(1,4)+6.25E-2*tt1*tt5+1.25E-1*X(1,1)*tt1*X(1,3&
&)+tt4+tt2)**2+5.0E-1*(tt46-1.25E-1*tt1*X(3,5)*X(3,8)+1.25E-1*tt1*&
&X(3,4)*X(3,8)+tt45+tt44+1.25E-1*tt1*X(3,6)*X(3,7)-1.25E-1*tt1*X(3&
&,3)*X(3,7)+tt42+tt41+tt39-1.25E-1*tt1*X(3,2)*X(3,6)+6.25E-2*tt1*t&
&t38-1.25E-1*tt1*X(3,4)*X(3,5)+1.25E-1*tt1*X(3,1)*X(3,5)+6.25E-2*t&
&t1*tt37-1.25E-1*tt1*X(3,1)*X(3,4)+tt36+1.25E-1*tt1*X(3,2)*X(3,3)+&
&tt34+tt32+tt31-1.25E-1*tt1*X(2,5)*X(2,8)+1.25E-1*tt1*X(2,4)*X(2,8&
&)+tt30+tt29+1.25E-1*tt1*X(2,6)*X(2,7)-1.25E-1*tt1*X(2,3)*X(2,7)+t&
&t27+tt26+tt24-1.25E-1*tt1*X(2,2)*X(2,6)+6.25E-2*tt1*tt23-1.25E-1*&
&tt1*X(2,4)*X(2,5)+1.25E-1*tt1*X(2,1)*X(2,5)+6.25E-2*tt1*tt22-1.25&
&E-1*tt1*X(2,1)*X(2,4)+tt21+1.25E-1*tt1*X(2,2)*X(2,3)+tt19+tt17+tt&
&16-1.25E-1*tt1*X(1,5)*X(1,8)+1.25E-1*tt1*X(1,4)*X(1,8)+tt15+tt14+&
&1.25E-1*tt1*X(1,6)*X(1,7)-1.25E-1*tt1*X(1,3)*X(1,7)+tt12+tt11+tt9&
&-1.25E-1*tt1*X(1,2)*X(1,6)+6.25E-2*tt1*tt8-1.25E-1*tt1*X(1,4)*X(1&
&,5)+1.25E-1*X(1,1)*tt1*X(1,5)+6.25E-2*tt1*tt7-1.25E-1*X(1,1)*tt1*&
&X(1,4)+tt6+1.25E-1*tt1*X(1,2)*X(1,3)+tt4+tt2)**2))
END 
SUBROUTINE hex_stvk_jac(jac, X, h, lam, miu) 
IMPLICIT NONE 
REAL(KIND=8) jac(1, 24) 
REAL(KIND=8) X(3, 8) 
REAL(KIND=8) h(1, 1) 
REAL(KIND=8) lam(1, 1) 
REAL(KIND=8) miu(1, 1) 
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
tt1 = h(1,1)**3
tt2 = 1/h(1,1)
tt3 = -2.5E-1*X(1,1)*tt2
tt4 = -2.5E-1*tt2*X(1,3)
tt5 = 2.5E-1*tt2*X(1,4)
tt6 = -2.5E-1*tt2*X(1,5)
tt7 = 2.5E-1*tt2*X(1,6)
tt8 = 2.5E-1*tt2*X(1,8)
tt9 = tt8-2.5E-1*tt2*X(1,7)+tt7+tt6+tt5+tt4+2.5E-1*tt2*X(1,2)+tt3&
&
tt10 = -2.5E-1*tt2*tt9
tt11 = -2.5E-1*tt2*X(1,2)
tt12 = 2.5E-1*tt2*X(1,7)
tt13 = tt8+tt12-2.5E-1*tt2*X(1,6)+tt6+tt5+2.5E-1*tt2*X(1,3)+tt11+&
&tt3
tt14 = -2.5E-1*tt2*tt13
tt15 = tt8+tt12+tt7+2.5E-1*tt2*X(1,5)-2.5E-1*tt2*X(1,4)+tt4+tt11+&
&tt3
tt16 = -2.5E-1*tt2*tt15
tt17 = -2.5E-1*tt2*X(2,1)
tt18 = -2.5E-1*tt2*X(2,3)
tt19 = 2.5E-1*tt2*X(2,4)
tt20 = -2.5E-1*tt2*X(2,5)
tt21 = 2.5E-1*tt2*X(2,6)
tt22 = 2.5E-1*tt2*X(2,8)
tt23 = tt22-2.5E-1*tt2*X(2,7)+tt21+tt20+tt19+tt18+2.5E-1*tt2*X(2,&
&2)+tt17
tt24 = -2.5E-1*tt2*X(3,1)
tt25 = -2.5E-1*tt2*X(3,3)
tt26 = 2.5E-1*tt2*X(3,4)
tt27 = -2.5E-1*tt2*X(3,5)
tt28 = 2.5E-1*tt2*X(3,6)
tt29 = 2.5E-1*tt2*X(3,8)
tt30 = tt29-2.5E-1*tt2*X(3,7)+tt28+tt27+tt26+tt25+2.5E-1*tt2*X(3,&
&2)+tt24
tt31 = tt30**2+tt23**2+tt9**2-1
tt32 = -2.5E-1*tt2*X(2,2)
tt33 = 2.5E-1*tt2*X(2,7)
tt34 = tt22+tt33-2.5E-1*tt2*X(2,6)+tt20+tt19+2.5E-1*tt2*X(2,3)+tt&
&32+tt17
tt35 = -2.5E-1*tt2*X(3,2)
tt36 = 2.5E-1*tt2*X(3,7)
tt37 = tt29+tt36-2.5E-1*tt2*X(3,6)+tt27+tt26+2.5E-1*tt2*X(3,3)+tt&
&35+tt24
tt38 = tt37**2+tt34**2+tt13**2-1
tt39 = tt22+tt33+tt21+2.5E-1*tt2*X(2,5)-2.5E-1*tt2*X(2,4)+tt18+tt&
&32+tt17
tt40 = tt29+tt36+tt28+2.5E-1*tt2*X(3,5)-2.5E-1*tt2*X(3,4)+tt25+tt&
&35+tt24
tt41 = tt40**2+tt39**2+tt15**2-1
tt42 = 5.0E-1*tt41+5.0E-1*tt38+5.0E-1*tt31
tt43 = 1/h(1,1)**2
tt44 = 1.25E-1*X(1,1)*tt43
tt45 = -1.25E-1*tt43*X(1,4)
tt46 = 1.25E-1*tt43*X(1,5)
tt47 = -1.25E-1*tt43*X(1,8)
tt48 = 6.25E-2*X(1,1)**2*tt43
tt49 = X(1,2)**2
tt50 = -6.25E-2*tt43*tt49
tt51 = X(1,3)**2
tt52 = -6.25E-2*tt43*tt51
tt53 = X(1,4)**2
tt54 = X(1,5)**2
tt55 = 1.25E-1*tt43*X(1,3)*X(1,6)
tt56 = X(1,6)**2
tt57 = -6.25E-2*tt43*tt56
tt58 = 1.25E-1*tt43*X(1,2)*X(1,7)
tt59 = X(1,7)**2
tt60 = -6.25E-2*tt43*tt59
tt61 = -1.25E-1*X(1,1)*tt43*X(1,8)
tt62 = 6.25E-2*tt43*X(1,8)**2
tt63 = 6.25E-2*tt43*X(2,1)**2
tt64 = X(2,2)**2
tt65 = -6.25E-2*tt43*tt64
tt66 = X(2,3)**2
tt67 = -6.25E-2*tt43*tt66
tt68 = X(2,4)**2
tt69 = X(2,5)**2
tt70 = 1.25E-1*tt43*X(2,3)*X(2,6)
tt71 = X(2,6)**2
tt72 = -6.25E-2*tt43*tt71
tt73 = 1.25E-1*tt43*X(2,2)*X(2,7)
tt74 = X(2,7)**2
tt75 = -6.25E-2*tt43*tt74
tt76 = -1.25E-1*tt43*X(2,1)*X(2,8)
tt77 = 6.25E-2*tt43*X(2,8)**2
tt78 = 6.25E-2*tt43*X(3,1)**2
tt79 = X(3,2)**2
tt80 = -6.25E-2*tt43*tt79
tt81 = X(3,3)**2
tt82 = -6.25E-2*tt43*tt81
tt83 = X(3,4)**2
tt84 = X(3,5)**2
tt85 = 1.25E-1*tt43*X(3,3)*X(3,6)
tt86 = X(3,6)**2
tt87 = -6.25E-2*tt43*tt86
tt88 = 1.25E-1*tt43*X(3,2)*X(3,7)
tt89 = X(3,7)**2
tt90 = -6.25E-2*tt43*tt89
tt91 = -1.25E-1*tt43*X(3,1)*X(3,8)
tt92 = 6.25E-2*tt43*X(3,8)**2
tt93 = tt92-1.25E-1*tt43*X(3,5)*X(3,8)+1.25E-1*tt43*X(3,4)*X(3,8)&
&+tt91+tt90+1.25E-1*tt43*X(3,6)*X(3,7)-1.25E-1*tt43*X(3,3)*X(3,7)+&
&tt88+tt87+tt85-1.25E-1*tt43*X(3,2)*X(3,6)+6.25E-2*tt43*tt84-1.25E&
&-1*tt43*X(3,4)*X(3,5)+1.25E-1*tt43*X(3,1)*X(3,5)+6.25E-2*tt43*tt8&
&3-1.25E-1*tt43*X(3,1)*X(3,4)+tt82+1.25E-1*tt43*X(3,2)*X(3,3)+tt80&
&+tt78+tt77-1.25E-1*tt43*X(2,5)*X(2,8)+1.25E-1*tt43*X(2,4)*X(2,8)+&
&tt76+tt75+1.25E-1*tt43*X(2,6)*X(2,7)-1.25E-1*tt43*X(2,3)*X(2,7)+t&
&t73+tt72+tt70-1.25E-1*tt43*X(2,2)*X(2,6)+6.25E-2*tt43*tt69-1.25E-&
&1*tt43*X(2,4)*X(2,5)+1.25E-1*tt43*X(2,1)*X(2,5)+6.25E-2*tt43*tt68&
&-1.25E-1*tt43*X(2,1)*X(2,4)+tt67+1.25E-1*tt43*X(2,2)*X(2,3)+tt65+&
&tt63+tt62-1.25E-1*tt43*X(1,5)*X(1,8)+1.25E-1*tt43*X(1,4)*X(1,8)+t&
&t61+tt60+1.25E-1*tt43*X(1,6)*X(1,7)-1.25E-1*tt43*X(1,3)*X(1,7)+tt&
&58+tt57+tt55-1.25E-1*tt43*X(1,2)*X(1,6)+6.25E-2*tt43*tt54-1.25E-1&
&*tt43*X(1,4)*X(1,5)+1.25E-1*X(1,1)*tt43*X(1,5)+6.25E-2*tt43*tt53-&
&1.25E-1*X(1,1)*tt43*X(1,4)+tt52+1.25E-1*tt43*X(1,2)*X(1,3)+tt50+t&
&t48
tt94 = 1.0E+0*(tt47+tt46+tt45+tt44)*tt93
tt95 = 1.25E-1*tt43*X(1,3)
tt96 = -1.25E-1*tt43*X(1,6)
tt97 = -6.25E-2*tt43*tt53
tt98 = 1.25E-1*tt43*X(1,4)*X(1,5)
tt99 = -6.25E-2*tt43*tt54
tt100 = -6.25E-2*tt43*tt68
tt101 = 1.25E-1*tt43*X(2,4)*X(2,5)
tt102 = -6.25E-2*tt43*tt69
tt103 = -6.25E-2*tt43*tt83
tt104 = 1.25E-1*tt43*X(3,4)*X(3,5)
tt105 = -6.25E-2*tt43*tt84
tt106 = tt92+1.25E-1*tt43*X(3,6)*X(3,8)-1.25E-1*tt43*X(3,3)*X(3,8&
&)+tt91+tt90-1.25E-1*tt43*X(3,5)*X(3,7)+1.25E-1*tt43*X(3,4)*X(3,7)&
&+tt88+6.25E-2*tt43*tt86-1.25E-1*tt43*X(3,3)*X(3,6)-1.25E-1*tt43*X&
&(3,1)*X(3,6)+tt105+tt104+1.25E-1*tt43*X(3,2)*X(3,5)+tt103-1.25E-1&
&*tt43*X(3,2)*X(3,4)+6.25E-2*tt43*tt81+1.25E-1*tt43*X(3,1)*X(3,3)+&
&tt80+tt78+tt77+1.25E-1*tt43*X(2,6)*X(2,8)-1.25E-1*tt43*X(2,3)*X(2&
&,8)+tt76+tt75-1.25E-1*tt43*X(2,5)*X(2,7)+1.25E-1*tt43*X(2,4)*X(2,&
&7)+tt73+6.25E-2*tt43*tt71-1.25E-1*tt43*X(2,3)*X(2,6)-1.25E-1*tt43&
&*X(2,1)*X(2,6)+tt102+tt101+1.25E-1*tt43*X(2,2)*X(2,5)+tt100-1.25E&
&-1*tt43*X(2,2)*X(2,4)+6.25E-2*tt43*tt66+1.25E-1*tt43*X(2,1)*X(2,3&
&)+tt65+tt63+tt62+1.25E-1*tt43*X(1,6)*X(1,8)-1.25E-1*tt43*X(1,3)*X&
&(1,8)+tt61+tt60-1.25E-1*tt43*X(1,5)*X(1,7)+1.25E-1*tt43*X(1,4)*X(&
&1,7)+tt58+6.25E-2*tt43*tt56-1.25E-1*tt43*X(1,3)*X(1,6)-1.25E-1*X(&
&1,1)*tt43*X(1,6)+tt99+tt98+1.25E-1*tt43*X(1,2)*X(1,5)+tt97-1.25E-&
&1*tt43*X(1,2)*X(1,4)+6.25E-2*tt43*tt51+1.25E-1*X(1,1)*tt43*X(1,3)&
&+tt50+tt48
tt107 = 1.0E+0*(tt47+tt96+tt95+tt44)*tt106
tt108 = 1.25E-1*tt43*X(1,2)
tt109 = -1.25E-1*tt43*X(1,7)
tt110 = tt92+1.25E-1*tt43*X(3,7)*X(3,8)-1.25E-1*tt43*X(3,2)*X(3,8&
&)+tt91+6.25E-2*tt43*tt89-1.25E-1*tt43*X(3,2)*X(3,7)-1.25E-1*tt43*&
&X(3,1)*X(3,7)+tt87-1.25E-1*tt43*X(3,5)*X(3,6)+1.25E-1*tt43*X(3,4)&
&*X(3,6)+tt85+tt105+tt104+1.25E-1*tt43*X(3,3)*X(3,5)+tt103-1.25E-1&
&*tt43*X(3,3)*X(3,4)+tt82+6.25E-2*tt43*tt79+1.25E-1*tt43*X(3,1)*X(&
&3,2)+tt78+tt77+1.25E-1*tt43*X(2,7)*X(2,8)-1.25E-1*tt43*X(2,2)*X(2&
&,8)+tt76+6.25E-2*tt43*tt74-1.25E-1*tt43*X(2,2)*X(2,7)-1.25E-1*tt4&
&3*X(2,1)*X(2,7)+tt72-1.25E-1*tt43*X(2,5)*X(2,6)+1.25E-1*tt43*X(2,&
&4)*X(2,6)+tt70+tt102+tt101+1.25E-1*tt43*X(2,3)*X(2,5)+tt100-1.25E&
&-1*tt43*X(2,3)*X(2,4)+tt67+6.25E-2*tt43*tt64+1.25E-1*tt43*X(2,1)*&
&X(2,2)+tt63+tt62+1.25E-1*tt43*X(1,7)*X(1,8)-1.25E-1*tt43*X(1,2)*X&
&(1,8)+tt61+6.25E-2*tt43*tt59-1.25E-1*tt43*X(1,2)*X(1,7)-1.25E-1*X&
&(1,1)*tt43*X(1,7)+tt57-1.25E-1*tt43*X(1,5)*X(1,6)+1.25E-1*tt43*X(&
&1,4)*X(1,6)+tt55+tt99+tt98+1.25E-1*tt43*X(1,3)*X(1,5)+tt97-1.25E-&
&1*tt43*X(1,3)*X(1,4)+tt52+6.25E-2*tt43*tt49+1.25E-1*X(1,1)*tt43*X&
&(1,2)+tt48
tt111 = 1.0E+0*(tt47+tt109+tt108+tt44)*tt110
tt112 = -2.5E-1*tt2*tt9*tt31
tt113 = -2.5E-1*tt2*tt13*tt38
tt114 = -2.5E-1*tt2*tt15*tt41
tt115 = -2.5E-1*tt2*tt23
tt116 = -2.5E-1*tt2*tt34
tt117 = -2.5E-1*tt2*tt39
tt118 = 1.25E-1*tt43*X(2,1)
tt119 = -1.25E-1*tt43*X(2,4)
tt120 = 1.25E-1*tt43*X(2,5)
tt121 = -1.25E-1*tt43*X(2,8)
tt122 = 1.0E+0*(tt121+tt120+tt119+tt118)*tt93
tt123 = 1.25E-1*tt43*X(2,3)
tt124 = -1.25E-1*tt43*X(2,6)
tt125 = 1.0E+0*(tt121+tt124+tt123+tt118)*tt106
tt126 = 1.25E-1*tt43*X(2,2)
tt127 = -1.25E-1*tt43*X(2,7)
tt128 = 1.0E+0*(tt121+tt127+tt126+tt118)*tt110
tt129 = -2.5E-1*tt2*tt23*tt31
tt130 = -2.5E-1*tt2*tt34*tt38
tt131 = -2.5E-1*tt2*tt39*tt41
tt132 = -2.5E-1*tt2*tt30
tt133 = -2.5E-1*tt2*tt37
tt134 = -2.5E-1*tt2*tt40
tt135 = 1.25E-1*tt43*X(3,1)
tt136 = -1.25E-1*tt43*X(3,4)
tt137 = 1.25E-1*tt43*X(3,5)
tt138 = -1.25E-1*tt43*X(3,8)
tt139 = 1.0E+0*(tt138+tt137+tt136+tt135)*tt93
tt140 = 1.25E-1*tt43*X(3,3)
tt141 = -1.25E-1*tt43*X(3,6)
tt142 = 1.0E+0*(tt138+tt141+tt140+tt135)*tt106
tt143 = 1.25E-1*tt43*X(3,2)
tt144 = -1.25E-1*tt43*X(3,7)
tt145 = 1.0E+0*(tt138+tt144+tt143+tt135)*tt110
tt146 = -2.5E-1*tt2*tt30*tt31
tt147 = -2.5E-1*tt2*tt37*tt38
tt148 = -2.5E-1*tt2*tt40*tt41
tt149 = 2.5E-1*tt2*tt9
tt150 = -1.25E-1*tt43*X(1,2)
tt151 = 1.25E-1*tt43*X(1,7)
tt152 = 1.0E+0*(tt151+tt96+tt95+tt150)*tt93
tt153 = 1.0E+0*(tt151+tt46+tt45+tt150)*tt106
tt154 = 2.5E-1*tt2*tt9*tt31
tt155 = 2.5E-1*tt2*tt23
tt156 = -1.25E-1*tt43*X(2,2)
tt157 = 1.25E-1*tt43*X(2,7)
tt158 = 1.0E+0*(tt157+tt124+tt123+tt156)*tt93
tt159 = 1.0E+0*(tt157+tt120+tt119+tt156)*tt106
tt160 = 2.5E-1*tt2*tt23*tt31
tt161 = 2.5E-1*tt2*tt30
tt162 = -1.25E-1*tt43*X(3,2)
tt163 = 1.25E-1*tt43*X(3,7)
tt164 = 1.0E+0*(tt163+tt141+tt140+tt162)*tt93
tt165 = 1.0E+0*(tt163+tt137+tt136+tt162)*tt106
tt166 = 2.5E-1*tt2*tt30*tt31
tt167 = 2.5E-1*tt2*tt13
tt168 = -1.25E-1*tt43*X(1,3)
tt169 = 1.25E-1*tt43*X(1,6)
tt170 = 1.0E+0*(tt109+tt169+tt168+tt108)*tt93
tt171 = 1.0E+0*(tt169+tt46+tt45+tt168)*tt110
tt172 = 2.5E-1*tt2*tt13*tt38
tt173 = 2.5E-1*tt2*tt34
tt174 = -1.25E-1*tt43*X(2,3)
tt175 = 1.25E-1*tt43*X(2,6)
tt176 = 1.0E+0*(tt127+tt175+tt174+tt126)*tt93
tt177 = 1.0E+0*(tt175+tt120+tt119+tt174)*tt110
tt178 = 2.5E-1*tt2*tt34*tt38
tt179 = 2.5E-1*tt2*tt37
tt180 = -1.25E-1*tt43*X(3,3)
tt181 = 1.25E-1*tt43*X(3,6)
tt182 = 1.0E+0*(tt144+tt181+tt180+tt143)*tt93
tt183 = 1.0E+0*(tt181+tt137+tt136+tt180)*tt110
tt184 = 2.5E-1*tt2*tt37*tt38
tt185 = -1.25E-1*X(1,1)*tt43
tt186 = 1.25E-1*tt43*X(1,4)
tt187 = -1.25E-1*tt43*X(1,5)
tt188 = 1.25E-1*tt43*X(1,8)
tt189 = 1.0E+0*(tt188+tt187+tt186+tt185)*tt93
tt190 = -1.25E-1*tt43*X(2,1)
tt191 = 1.25E-1*tt43*X(2,4)
tt192 = -1.25E-1*tt43*X(2,5)
tt193 = 1.25E-1*tt43*X(2,8)
tt194 = 1.0E+0*(tt193+tt192+tt191+tt190)*tt93
tt195 = -1.25E-1*tt43*X(3,1)
tt196 = 1.25E-1*tt43*X(3,4)
tt197 = -1.25E-1*tt43*X(3,5)
tt198 = 1.25E-1*tt43*X(3,8)
tt199 = 1.0E+0*(tt198+tt197+tt196+tt195)*tt93
tt200 = 2.5E-1*tt2*tt15
tt201 = 1.0E+0*(tt109+tt187+tt186+tt108)*tt106
tt202 = 1.0E+0*(tt96+tt187+tt186+tt95)*tt110
tt203 = 2.5E-1*tt2*tt15*tt41
tt204 = 2.5E-1*tt2*tt39
tt205 = 1.0E+0*(tt127+tt192+tt191+tt126)*tt106
tt206 = 1.0E+0*(tt124+tt192+tt191+tt123)*tt110
tt207 = 2.5E-1*tt2*tt39*tt41
tt208 = 2.5E-1*tt2*tt40
tt209 = 1.0E+0*(tt144+tt197+tt196+tt143)*tt106
tt210 = 1.0E+0*(tt141+tt197+tt196+tt140)*tt110
tt211 = 2.5E-1*tt2*tt40*tt41
tt212 = 1.0E+0*(tt188+tt169+tt168+tt185)*tt106
tt213 = 1.0E+0*(tt193+tt175+tt174+tt190)*tt106
tt214 = 1.0E+0*(tt198+tt181+tt180+tt195)*tt106
tt215 = 1.0E+0*(tt188+tt151+tt150+tt185)*tt110
tt216 = 1.0E+0*(tt193+tt157+tt156+tt190)*tt110
tt217 = 1.0E+0*(tt198+tt163+tt162+tt195)*tt110
jac(1,1) = tt1*(miu(1,1)*(tt114+tt113+tt112+tt111+tt107+tt94)+1.0&
&E+0*lam(1,1)*(tt16+tt14+tt10)*tt42)
jac(1,2) = tt1*(miu(1,1)*(tt131+tt130+tt129+tt128+tt125+tt122)+1.&
&0E+0*lam(1,1)*(tt117+tt116+tt115)*tt42)
jac(1,3) = tt1*(miu(1,1)*(tt148+tt147+tt146+tt145+tt142+tt139)+1.&
&0E+0*lam(1,1)*(tt134+tt133+tt132)*tt42)
jac(1,4) = tt1*(miu(1,1)*(tt114+tt113+tt154+tt111+tt153+tt152)+1.&
&0E+0*lam(1,1)*(tt16+tt14+tt149)*tt42)
jac(1,5) = tt1*(miu(1,1)*(tt131+tt130+tt160+tt128+tt159+tt158)+1.&
&0E+0*lam(1,1)*(tt117+tt116+tt155)*tt42)
jac(1,6) = tt1*(miu(1,1)*(tt148+tt147+tt166+tt145+tt165+tt164)+1.&
&0E+0*lam(1,1)*(tt134+tt133+tt161)*tt42)
jac(1,7) = tt1*(miu(1,1)*(tt114+tt172+tt112+tt171+tt107+tt170)+1.&
&0E+0*lam(1,1)*(tt16+tt167+tt10)*tt42)
jac(1,8) = tt1*(miu(1,1)*(tt131+tt178+tt129+tt177+tt125+tt176)+1.&
&0E+0*lam(1,1)*(tt117+tt173+tt115)*tt42)
jac(1,9) = tt1*(miu(1,1)*(tt148+tt184+tt146+tt183+tt142+tt182)+1.&
&0E+0*lam(1,1)*(tt134+tt179+tt132)*tt42)
jac(1,10) = tt1*(miu(1,1)*(tt114+tt172+tt154+tt171+tt153+tt189)+1&
&.0E+0*lam(1,1)*(tt16+tt167+tt149)*tt42)
jac(1,11) = tt1*(miu(1,1)*(tt131+tt178+tt160+tt177+tt159+tt194)+1&
&.0E+0*lam(1,1)*(tt117+tt173+tt155)*tt42)
jac(1,12) = tt1*(miu(1,1)*(tt148+tt184+tt166+tt183+tt165+tt199)+1&
&.0E+0*lam(1,1)*(tt134+tt179+tt161)*tt42)
jac(1,13) = tt1*(miu(1,1)*(tt203+tt113+tt112+tt202+tt201+tt94)+1.&
&0E+0*lam(1,1)*(tt200+tt14+tt10)*tt42)
jac(1,14) = tt1*(miu(1,1)*(tt207+tt130+tt129+tt206+tt205+tt122)+1&
&.0E+0*lam(1,1)*(tt204+tt116+tt115)*tt42)
jac(1,15) = tt1*(miu(1,1)*(tt211+tt147+tt146+tt210+tt209+tt139)+1&
&.0E+0*lam(1,1)*(tt208+tt133+tt132)*tt42)
jac(1,16) = tt1*(miu(1,1)*(tt203+tt113+tt154+tt202+tt212+tt152)+1&
&.0E+0*lam(1,1)*(tt200+tt14+tt149)*tt42)
jac(1,17) = tt1*(miu(1,1)*(tt207+tt130+tt160+tt206+tt213+tt158)+1&
&.0E+0*lam(1,1)*(tt204+tt116+tt155)*tt42)
jac(1,18) = tt1*(miu(1,1)*(tt211+tt147+tt166+tt210+tt214+tt164)+1&
&.0E+0*lam(1,1)*(tt208+tt133+tt161)*tt42)
jac(1,19) = tt1*(miu(1,1)*(tt203+tt172+tt112+tt215+tt201+tt170)+1&
&.0E+0*lam(1,1)*(tt200+tt167+tt10)*tt42)
jac(1,20) = tt1*(miu(1,1)*(tt207+tt178+tt129+tt216+tt205+tt176)+1&
&.0E+0*lam(1,1)*(tt204+tt173+tt115)*tt42)
jac(1,21) = tt1*(miu(1,1)*(tt211+tt184+tt146+tt217+tt209+tt182)+1&
&.0E+0*lam(1,1)*(tt208+tt179+tt132)*tt42)
jac(1,22) = tt1*(miu(1,1)*(tt203+tt172+tt154+tt215+tt212+tt189)+1&
&.0E+0*lam(1,1)*(tt200+tt167+tt149)*tt42)
jac(1,23) = tt1*(miu(1,1)*(tt207+tt178+tt160+tt216+tt213+tt194)+1&
&.0E+0*lam(1,1)*(tt204+tt173+tt155)*tt42)
jac(1,24) = tt1*(miu(1,1)*(tt211+tt184+tt166+tt217+tt214+tt199)+1&
&.0E+0*lam(1,1)*(tt208+tt179+tt161)*tt42)
END 
SUBROUTINE hex_stvk_hes(hes, X, h, lam, miu) 
IMPLICIT NONE 
REAL(KIND=8) hes(24, 24) 
REAL(KIND=8) X(3, 8) 
REAL(KIND=8) h(1, 1) 
REAL(KIND=8) lam(1, 1) 
REAL(KIND=8) miu(1, 1) 
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
tt1 = h(1,1)**3
tt2 = 1/h(1,1)
tt3 = -2.5E-1*X(1,1)*tt2
tt4 = -2.5E-1*tt2*X(1,3)
tt5 = 2.5E-1*tt2*X(1,4)
tt6 = -2.5E-1*tt2*X(1,5)
tt7 = 2.5E-1*tt2*X(1,6)
tt8 = 2.5E-1*tt2*X(1,8)
tt9 = tt8-2.5E-1*tt2*X(1,7)+tt7+tt6+tt5+tt4+2.5E-1*tt2*X(1,2)+tt3&
&
tt10 = -2.5E-1*tt2*tt9
tt11 = -2.5E-1*tt2*X(1,2)
tt12 = 2.5E-1*tt2*X(1,7)
tt13 = tt8+tt12-2.5E-1*tt2*X(1,6)+tt6+tt5+2.5E-1*tt2*X(1,3)+tt11+&
&tt3
tt14 = -2.5E-1*tt2*tt13
tt15 = tt8+tt12+tt7+2.5E-1*tt2*X(1,5)-2.5E-1*tt2*X(1,4)+tt4+tt11+&
&tt3
tt16 = -2.5E-1*tt2*tt15
tt17 = tt16+tt14+tt10
tt18 = 1/h(1,1)**2
tt19 = tt9**2
tt20 = -2.5E-1*tt2*X(2,1)
tt21 = -2.5E-1*tt2*X(2,3)
tt22 = 2.5E-1*tt2*X(2,4)
tt23 = -2.5E-1*tt2*X(2,5)
tt24 = 2.5E-1*tt2*X(2,6)
tt25 = 2.5E-1*tt2*X(2,8)
tt26 = tt25-2.5E-1*tt2*X(2,7)+tt24+tt23+tt22+tt21+2.5E-1*tt2*X(2,&
&2)+tt20
tt27 = tt26**2
tt28 = -2.5E-1*tt2*X(3,1)
tt29 = -2.5E-1*tt2*X(3,3)
tt30 = 2.5E-1*tt2*X(3,4)
tt31 = -2.5E-1*tt2*X(3,5)
tt32 = 2.5E-1*tt2*X(3,6)
tt33 = 2.5E-1*tt2*X(3,8)
tt34 = tt33-2.5E-1*tt2*X(3,7)+tt32+tt31+tt30+tt29+2.5E-1*tt2*X(3,&
&2)+tt28
tt35 = tt34**2
tt36 = tt35+tt27+tt19-1
tt37 = tt13**2
tt38 = -2.5E-1*tt2*X(2,2)
tt39 = 2.5E-1*tt2*X(2,7)
tt40 = tt25+tt39-2.5E-1*tt2*X(2,6)+tt23+tt22+2.5E-1*tt2*X(2,3)+tt&
&38+tt20
tt41 = tt40**2
tt42 = -2.5E-1*tt2*X(3,2)
tt43 = 2.5E-1*tt2*X(3,7)
tt44 = tt33+tt43-2.5E-1*tt2*X(3,6)+tt31+tt30+2.5E-1*tt2*X(3,3)+tt&
&42+tt28
tt45 = tt44**2
tt46 = tt45+tt41+tt37-1
tt47 = tt15**2
tt48 = tt25+tt39+tt24+2.5E-1*tt2*X(2,5)-2.5E-1*tt2*X(2,4)+tt21+tt&
&38+tt20
tt49 = tt48**2
tt50 = tt33+tt43+tt32+2.5E-1*tt2*X(3,5)-2.5E-1*tt2*X(3,4)+tt29+tt&
&42+tt28
tt51 = tt50**2
tt52 = tt51+tt49+tt47-1
tt53 = 5.0E-1*tt52+5.0E-1*tt46+5.0E-1*tt36
tt54 = 1.875E-1*tt18*lam(1,1)*tt53
tt55 = 1.25E-1*X(1,1)*tt18
tt56 = -1.25E-1*tt18*X(1,4)
tt57 = 1.25E-1*tt18*X(1,5)
tt58 = -1.25E-1*tt18*X(1,8)
tt59 = tt58+tt57+tt56+tt55
tt60 = 1.0E+0*tt59**2
tt61 = 1.25E-1*tt18*X(1,3)
tt62 = -1.25E-1*tt18*X(1,6)
tt63 = tt58+tt62+tt61+tt55
tt64 = 1.0E+0*tt63**2
tt65 = 1.25E-1*tt18*X(1,2)
tt66 = -1.25E-1*tt18*X(1,7)
tt67 = tt58+tt66+tt65+tt55
tt68 = 1.0E+0*tt67**2
tt69 = 1.25E-1*tt18*tt19
tt70 = 1.25E-1*tt18*tt37
tt71 = 1.25E-1*tt18*tt47
tt72 = 6.25E-2*X(1,1)**2*tt18
tt73 = X(1,2)**2
tt74 = -6.25E-2*tt18*tt73
tt75 = X(1,3)**2
tt76 = -6.25E-2*tt18*tt75
tt77 = X(1,4)**2
tt78 = X(1,5)**2
tt79 = 1.25E-1*tt18*X(1,3)*X(1,6)
tt80 = X(1,6)**2
tt81 = -6.25E-2*tt18*tt80
tt82 = 1.25E-1*tt18*X(1,2)*X(1,7)
tt83 = X(1,7)**2
tt84 = -6.25E-2*tt18*tt83
tt85 = -1.25E-1*X(1,1)*tt18*X(1,8)
tt86 = 6.25E-2*tt18*X(1,8)**2
tt87 = 6.25E-2*tt18*X(2,1)**2
tt88 = X(2,2)**2
tt89 = -6.25E-2*tt18*tt88
tt90 = X(2,3)**2
tt91 = -6.25E-2*tt18*tt90
tt92 = X(2,4)**2
tt93 = X(2,5)**2
tt94 = 1.25E-1*tt18*X(2,3)*X(2,6)
tt95 = X(2,6)**2
tt96 = -6.25E-2*tt18*tt95
tt97 = 1.25E-1*tt18*X(2,2)*X(2,7)
tt98 = X(2,7)**2
tt99 = -6.25E-2*tt18*tt98
tt100 = -1.25E-1*tt18*X(2,1)*X(2,8)
tt101 = 6.25E-2*tt18*X(2,8)**2
tt102 = 6.25E-2*tt18*X(3,1)**2
tt103 = X(3,2)**2
tt104 = -6.25E-2*tt18*tt103
tt105 = X(3,3)**2
tt106 = -6.25E-2*tt18*tt105
tt107 = X(3,4)**2
tt108 = X(3,5)**2
tt109 = 1.25E-1*tt18*X(3,3)*X(3,6)
tt110 = X(3,6)**2
tt111 = -6.25E-2*tt18*tt110
tt112 = 1.25E-1*tt18*X(3,2)*X(3,7)
tt113 = X(3,7)**2
tt114 = -6.25E-2*tt18*tt113
tt115 = -1.25E-1*tt18*X(3,1)*X(3,8)
tt116 = 6.25E-2*tt18*X(3,8)**2
tt117 = tt116-1.25E-1*tt18*X(3,5)*X(3,8)+1.25E-1*tt18*X(3,4)*X(3,&
&8)+tt115+tt114+1.25E-1*tt18*X(3,6)*X(3,7)-1.25E-1*tt18*X(3,3)*X(3&
&,7)+tt112+tt111+tt109-1.25E-1*tt18*X(3,2)*X(3,6)+6.25E-2*tt18*tt1&
&08-1.25E-1*tt18*X(3,4)*X(3,5)+1.25E-1*tt18*X(3,1)*X(3,5)+6.25E-2*&
&tt18*tt107-1.25E-1*tt18*X(3,1)*X(3,4)+tt106+1.25E-1*tt18*X(3,2)*X&
&(3,3)+tt104+tt102+tt101-1.25E-1*tt18*X(2,5)*X(2,8)+1.25E-1*tt18*X&
&(2,4)*X(2,8)+tt100+tt99+1.25E-1*tt18*X(2,6)*X(2,7)-1.25E-1*tt18*X&
&(2,3)*X(2,7)+tt97+tt96+tt94-1.25E-1*tt18*X(2,2)*X(2,6)+6.25E-2*tt&
&18*tt93-1.25E-1*tt18*X(2,4)*X(2,5)+1.25E-1*tt18*X(2,1)*X(2,5)+6.2&
&5E-2*tt18*tt92-1.25E-1*tt18*X(2,1)*X(2,4)+tt91+1.25E-1*tt18*X(2,2&
&)*X(2,3)+tt89+tt87+tt86-1.25E-1*tt18*X(1,5)*X(1,8)+1.25E-1*tt18*X&
&(1,4)*X(1,8)+tt85+tt84+1.25E-1*tt18*X(1,6)*X(1,7)-1.25E-1*tt18*X(&
&1,3)*X(1,7)+tt82+tt81+tt79-1.25E-1*tt18*X(1,2)*X(1,6)+6.25E-2*tt1&
&8*tt78-1.25E-1*tt18*X(1,4)*X(1,5)+1.25E-1*X(1,1)*tt18*X(1,5)+6.25&
&E-2*tt18*tt77-1.25E-1*X(1,1)*tt18*X(1,4)+tt76+1.25E-1*tt18*X(1,2)&
&*X(1,3)+tt74+tt72
tt118 = 1.25E-1*tt18*tt117
tt119 = -6.25E-2*tt18*tt77
tt120 = 1.25E-1*tt18*X(1,4)*X(1,5)
tt121 = -6.25E-2*tt18*tt78
tt122 = -6.25E-2*tt18*tt92
tt123 = 1.25E-1*tt18*X(2,4)*X(2,5)
tt124 = -6.25E-2*tt18*tt93
tt125 = -6.25E-2*tt18*tt107
tt126 = 1.25E-1*tt18*X(3,4)*X(3,5)
tt127 = -6.25E-2*tt18*tt108
tt128 = tt116+1.25E-1*tt18*X(3,6)*X(3,8)-1.25E-1*tt18*X(3,3)*X(3,&
&8)+tt115+tt114-1.25E-1*tt18*X(3,5)*X(3,7)+1.25E-1*tt18*X(3,4)*X(3&
&,7)+tt112+6.25E-2*tt18*tt110-1.25E-1*tt18*X(3,3)*X(3,6)-1.25E-1*t&
&t18*X(3,1)*X(3,6)+tt127+tt126+1.25E-1*tt18*X(3,2)*X(3,5)+tt125-1.&
&25E-1*tt18*X(3,2)*X(3,4)+6.25E-2*tt18*tt105+1.25E-1*tt18*X(3,1)*X&
&(3,3)+tt104+tt102+tt101+1.25E-1*tt18*X(2,6)*X(2,8)-1.25E-1*tt18*X&
&(2,3)*X(2,8)+tt100+tt99-1.25E-1*tt18*X(2,5)*X(2,7)+1.25E-1*tt18*X&
&(2,4)*X(2,7)+tt97+6.25E-2*tt18*tt95-1.25E-1*tt18*X(2,3)*X(2,6)-1.&
&25E-1*tt18*X(2,1)*X(2,6)+tt124+tt123+1.25E-1*tt18*X(2,2)*X(2,5)+t&
&t122-1.25E-1*tt18*X(2,2)*X(2,4)+6.25E-2*tt18*tt90+1.25E-1*tt18*X(&
&2,1)*X(2,3)+tt89+tt87+tt86+1.25E-1*tt18*X(1,6)*X(1,8)-1.25E-1*tt1&
&8*X(1,3)*X(1,8)+tt85+tt84-1.25E-1*tt18*X(1,5)*X(1,7)+1.25E-1*tt18&
&*X(1,4)*X(1,7)+tt82+6.25E-2*tt18*tt80-1.25E-1*tt18*X(1,3)*X(1,6)-&
&1.25E-1*X(1,1)*tt18*X(1,6)+tt121+tt120+1.25E-1*tt18*X(1,2)*X(1,5)&
&+tt119-1.25E-1*tt18*X(1,2)*X(1,4)+6.25E-2*tt18*tt75+1.25E-1*X(1,1&
&)*tt18*X(1,3)+tt74+tt72
tt129 = 1.25E-1*tt18*tt128
tt130 = tt116+1.25E-1*tt18*X(3,7)*X(3,8)-1.25E-1*tt18*X(3,2)*X(3,&
&8)+tt115+6.25E-2*tt18*tt113-1.25E-1*tt18*X(3,2)*X(3,7)-1.25E-1*tt&
&18*X(3,1)*X(3,7)+tt111-1.25E-1*tt18*X(3,5)*X(3,6)+1.25E-1*tt18*X(&
&3,4)*X(3,6)+tt109+tt127+tt126+1.25E-1*tt18*X(3,3)*X(3,5)+tt125-1.&
&25E-1*tt18*X(3,3)*X(3,4)+tt106+6.25E-2*tt18*tt103+1.25E-1*tt18*X(&
&3,1)*X(3,2)+tt102+tt101+1.25E-1*tt18*X(2,7)*X(2,8)-1.25E-1*tt18*X&
&(2,2)*X(2,8)+tt100+6.25E-2*tt18*tt98-1.25E-1*tt18*X(2,2)*X(2,7)-1&
&.25E-1*tt18*X(2,1)*X(2,7)+tt96-1.25E-1*tt18*X(2,5)*X(2,6)+1.25E-1&
&*tt18*X(2,4)*X(2,6)+tt94+tt124+tt123+1.25E-1*tt18*X(2,3)*X(2,5)+t&
&t122-1.25E-1*tt18*X(2,3)*X(2,4)+tt91+6.25E-2*tt18*tt88+1.25E-1*tt&
&18*X(2,1)*X(2,2)+tt87+tt86+1.25E-1*tt18*X(1,7)*X(1,8)-1.25E-1*tt1&
&8*X(1,2)*X(1,8)+tt85+6.25E-2*tt18*tt83-1.25E-1*tt18*X(1,2)*X(1,7)&
&-1.25E-1*X(1,1)*tt18*X(1,7)+tt81-1.25E-1*tt18*X(1,5)*X(1,6)+1.25E&
&-1*tt18*X(1,4)*X(1,6)+tt79+tt121+tt120+1.25E-1*tt18*X(1,3)*X(1,5)&
&+tt119-1.25E-1*tt18*X(1,3)*X(1,4)+tt76+6.25E-2*tt18*tt73+1.25E-1*&
&X(1,1)*tt18*X(1,2)+tt72
tt131 = 1.25E-1*tt18*tt130
tt132 = 6.25E-2*tt18*tt36
tt133 = 6.25E-2*tt18*tt46
tt134 = 6.25E-2*tt18*tt52
tt135 = -2.5E-1*tt2*tt26
tt136 = -2.5E-1*tt2*tt40
tt137 = -2.5E-1*tt2*tt48
tt138 = tt137+tt136+tt135
tt139 = 1.25E-1*tt18*X(2,1)
tt140 = -1.25E-1*tt18*X(2,4)
tt141 = 1.25E-1*tt18*X(2,5)
tt142 = -1.25E-1*tt18*X(2,8)
tt143 = tt142+tt141+tt140+tt139
tt144 = 1.0E+0*tt59*tt143
tt145 = 1.25E-1*tt18*X(2,3)
tt146 = -1.25E-1*tt18*X(2,6)
tt147 = tt142+tt146+tt145+tt139
tt148 = 1.0E+0*tt63*tt147
tt149 = 1.25E-1*tt18*X(2,2)
tt150 = -1.25E-1*tt18*X(2,7)
tt151 = tt142+tt150+tt149+tt139
tt152 = 1.0E+0*tt67*tt151
tt153 = 1.25E-1*tt18*tt9*tt26
tt154 = 1.25E-1*tt18*tt13*tt40
tt155 = 1.25E-1*tt18*tt15*tt48
tt156 = tt1*(miu(1,1)*(tt155+tt154+tt153+tt152+tt148+tt144)+1.0E+&
&0*lam(1,1)*tt17*tt138)
tt157 = -2.5E-1*tt2*tt34
tt158 = -2.5E-1*tt2*tt44
tt159 = -2.5E-1*tt2*tt50
tt160 = tt159+tt158+tt157
tt161 = 1.25E-1*tt18*X(3,1)
tt162 = -1.25E-1*tt18*X(3,4)
tt163 = 1.25E-1*tt18*X(3,5)
tt164 = -1.25E-1*tt18*X(3,8)
tt165 = tt164+tt163+tt162+tt161
tt166 = 1.0E+0*tt59*tt165
tt167 = 1.25E-1*tt18*X(3,3)
tt168 = -1.25E-1*tt18*X(3,6)
tt169 = tt164+tt168+tt167+tt161
tt170 = 1.0E+0*tt63*tt169
tt171 = 1.25E-1*tt18*X(3,2)
tt172 = -1.25E-1*tt18*X(3,7)
tt173 = tt164+tt172+tt171+tt161
tt174 = 1.0E+0*tt67*tt173
tt175 = 1.25E-1*tt18*tt9*tt34
tt176 = 1.25E-1*tt18*tt13*tt44
tt177 = 1.25E-1*tt18*tt15*tt50
tt178 = tt1*(miu(1,1)*(tt177+tt176+tt175+tt174+tt170+tt166)+1.0E+&
&0*lam(1,1)*tt17*tt160)
tt179 = 2.5E-1*tt2*tt9
tt180 = tt16+tt14+tt179
tt181 = 6.25E-2*tt18*lam(1,1)*tt53
tt182 = -1.25E-1*tt18*X(1,2)
tt183 = 1.25E-1*tt18*X(1,7)
tt184 = tt183+tt62+tt61+tt182
tt185 = 1.0E+0*tt184*tt59
tt186 = tt183+tt57+tt56+tt182
tt187 = 1.0E+0*tt186*tt63
tt188 = -1.25E-1*tt18*tt19
tt189 = -6.25E-2*tt18*tt36
tt190 = tt1*(miu(1,1)*(tt134+tt133+tt189+tt131+tt71+tt70+tt188+tt&
&68+tt187+tt185)+tt181+1.0E+0*lam(1,1)*tt17*tt180)
tt191 = 2.5E-1*tt2*tt26
tt192 = tt137+tt136+tt191
tt193 = -1.25E-1*tt18*X(2,2)
tt194 = 1.25E-1*tt18*X(2,7)
tt195 = tt194+tt141+tt140+tt193
tt196 = 1.0E+0*tt63*tt195
tt197 = tt194+tt146+tt145+tt193
tt198 = 1.0E+0*tt59*tt197
tt199 = -1.25E-1*tt18*tt9*tt26
tt200 = tt1*(miu(1,1)*(tt155+tt154+tt199+tt152+tt198+tt196)+1.0E+&
&0*lam(1,1)*tt17*tt192)
tt201 = 2.5E-1*tt2*tt34
tt202 = tt159+tt158+tt201
tt203 = -1.25E-1*tt18*X(3,2)
tt204 = 1.25E-1*tt18*X(3,7)
tt205 = tt204+tt163+tt162+tt203
tt206 = 1.0E+0*tt63*tt205
tt207 = tt204+tt168+tt167+tt203
tt208 = 1.0E+0*tt59*tt207
tt209 = -1.25E-1*tt18*tt9*tt34
tt210 = tt1*(miu(1,1)*(tt177+tt176+tt209+tt174+tt208+tt206)+1.0E+&
&0*lam(1,1)*tt17*tt202)
tt211 = 2.5E-1*tt2*tt13
tt212 = tt16+tt211+tt10
tt213 = -1.25E-1*tt18*X(1,3)
tt214 = 1.25E-1*tt18*X(1,6)
tt215 = tt66+tt214+tt213+tt65
tt216 = 1.0E+0*tt215*tt59
tt217 = tt214+tt57+tt56+tt213
tt218 = 1.0E+0*tt217*tt67
tt219 = -1.25E-1*tt18*tt37
tt220 = -6.25E-2*tt18*tt46
tt221 = tt1*(miu(1,1)*(tt134+tt220+tt132+tt129+tt71+tt219+tt69+tt&
&64+tt218+tt216)+tt181+1.0E+0*lam(1,1)*tt17*tt212)
tt222 = 2.5E-1*tt2*tt40
tt223 = tt137+tt222+tt135
tt224 = -1.25E-1*tt18*X(2,3)
tt225 = 1.25E-1*tt18*X(2,6)
tt226 = tt225+tt141+tt140+tt224
tt227 = 1.0E+0*tt67*tt226
tt228 = tt150+tt225+tt224+tt149
tt229 = 1.0E+0*tt59*tt228
tt230 = -1.25E-1*tt18*tt13*tt40
tt231 = tt1*(miu(1,1)*(tt155+tt230+tt153+tt148+tt229+tt227)+1.0E+&
&0*lam(1,1)*tt17*tt223)
tt232 = 2.5E-1*tt2*tt44
tt233 = tt159+tt232+tt157
tt234 = -1.25E-1*tt18*X(3,3)
tt235 = 1.25E-1*tt18*X(3,6)
tt236 = tt235+tt163+tt162+tt234
tt237 = 1.0E+0*tt67*tt236
tt238 = tt172+tt235+tt234+tt171
tt239 = 1.0E+0*tt59*tt238
tt240 = -1.25E-1*tt18*tt13*tt44
tt241 = tt1*(miu(1,1)*(tt177+tt240+tt175+tt170+tt239+tt237)+1.0E+&
&0*lam(1,1)*tt17*tt233)
tt242 = tt16+tt211+tt179
tt243 = -6.25E-2*tt18*lam(1,1)*tt53
tt244 = -1.25E-1*X(1,1)*tt18
tt245 = 1.25E-1*tt18*X(1,4)
tt246 = -1.25E-1*tt18*X(1,5)
tt247 = 1.25E-1*tt18*X(1,8)
tt248 = tt247+tt246+tt245+tt244
tt249 = 1.0E+0*tt59*tt248
tt250 = -1.25E-1*tt18*tt117
tt251 = tt1*(miu(1,1)*(tt134+tt220+tt189+tt250+tt71+tt219+tt188+t&
&t249+tt218+tt187)+tt243+1.0E+0*lam(1,1)*tt17*tt242)
tt252 = tt137+tt222+tt191
tt253 = -1.25E-1*tt18*X(2,1)
tt254 = 1.25E-1*tt18*X(2,4)
tt255 = -1.25E-1*tt18*X(2,5)
tt256 = 1.25E-1*tt18*X(2,8)
tt257 = tt256+tt255+tt254+tt253
tt258 = 1.0E+0*tt59*tt257
tt259 = tt1*(miu(1,1)*(tt155+tt230+tt199+tt258+tt196+tt227)+1.0E+&
&0*lam(1,1)*tt17*tt252)
tt260 = tt159+tt232+tt201
tt261 = -1.25E-1*tt18*X(3,1)
tt262 = 1.25E-1*tt18*X(3,4)
tt263 = -1.25E-1*tt18*X(3,5)
tt264 = 1.25E-1*tt18*X(3,8)
tt265 = tt264+tt263+tt262+tt261
tt266 = 1.0E+0*tt59*tt265
tt267 = tt1*(miu(1,1)*(tt177+tt240+tt209+tt266+tt206+tt237)+1.0E+&
&0*lam(1,1)*tt17*tt260)
tt268 = 2.5E-1*tt2*tt15
tt269 = tt268+tt14+tt10
tt270 = tt66+tt246+tt245+tt65
tt271 = 1.0E+0*tt270*tt63
tt272 = tt62+tt246+tt245+tt61
tt273 = 1.0E+0*tt272*tt67
tt274 = -1.25E-1*tt18*tt47
tt275 = -6.25E-2*tt18*tt52
tt276 = tt1*(miu(1,1)*(tt275+tt133+tt132+tt118+tt274+tt70+tt69+tt&
&60+tt273+tt271)+tt181+1.0E+0*lam(1,1)*tt17*tt269)
tt277 = 2.5E-1*tt2*tt48
tt278 = tt277+tt136+tt135
tt279 = tt146+tt255+tt254+tt145
tt280 = 1.0E+0*tt67*tt279
tt281 = tt150+tt255+tt254+tt149
tt282 = 1.0E+0*tt63*tt281
tt283 = -1.25E-1*tt18*tt15*tt48
tt284 = tt1*(miu(1,1)*(tt283+tt154+tt153+tt144+tt282+tt280)+1.0E+&
&0*lam(1,1)*tt17*tt278)
tt285 = 2.5E-1*tt2*tt50
tt286 = tt285+tt158+tt157
tt287 = tt168+tt263+tt262+tt167
tt288 = 1.0E+0*tt67*tt287
tt289 = tt172+tt263+tt262+tt171
tt290 = 1.0E+0*tt63*tt289
tt291 = -1.25E-1*tt18*tt15*tt50
tt292 = tt1*(miu(1,1)*(tt291+tt176+tt175+tt166+tt290+tt288)+1.0E+&
&0*lam(1,1)*tt17*tt286)
tt293 = tt268+tt14+tt179
tt294 = tt247+tt214+tt213+tt244
tt295 = 1.0E+0*tt63*tt294
tt296 = -1.25E-1*tt18*tt128
tt297 = tt1*(miu(1,1)*(tt275+tt133+tt189+tt296+tt274+tt70+tt188+t&
&t295+tt273+tt185)+tt243+1.0E+0*lam(1,1)*tt17*tt293)
tt298 = tt277+tt136+tt191
tt299 = tt256+tt225+tt224+tt253
tt300 = 1.0E+0*tt63*tt299
tt301 = tt1*(miu(1,1)*(tt283+tt154+tt199+tt300+tt198+tt280)+1.0E+&
&0*lam(1,1)*tt17*tt298)
tt302 = tt285+tt158+tt201
tt303 = tt264+tt235+tt234+tt261
tt304 = 1.0E+0*tt63*tt303
tt305 = tt1*(miu(1,1)*(tt291+tt176+tt209+tt304+tt208+tt288)+1.0E+&
&0*lam(1,1)*tt17*tt302)
tt306 = tt268+tt211+tt10
tt307 = tt247+tt183+tt182+tt244
tt308 = 1.0E+0*tt67*tt307
tt309 = -1.25E-1*tt18*tt130
tt310 = tt1*(miu(1,1)*(tt275+tt220+tt132+tt309+tt274+tt219+tt69+t&
&t308+tt271+tt216)+tt243+1.0E+0*lam(1,1)*tt17*tt306)
tt311 = tt277+tt222+tt135
tt312 = tt256+tt194+tt193+tt253
tt313 = 1.0E+0*tt67*tt312
tt314 = tt1*(miu(1,1)*(tt283+tt230+tt153+tt313+tt229+tt282)+1.0E+&
&0*lam(1,1)*tt17*tt311)
tt315 = tt285+tt232+tt157
tt316 = tt264+tt204+tt203+tt261
tt317 = 1.0E+0*tt67*tt316
tt318 = tt1*(miu(1,1)*(tt291+tt240+tt175+tt317+tt239+tt290)+1.0E+&
&0*lam(1,1)*tt17*tt315)
tt319 = tt268+tt211+tt179
tt320 = -1.875E-1*tt18*lam(1,1)*tt53
tt321 = tt1*(miu(1,1)*(tt275+tt220+tt189+tt309+tt296+tt250+tt274+&
&tt219+tt188+tt308+tt295+tt249)+tt320+1.0E+0*lam(1,1)*tt17*tt319)
tt322 = tt277+tt222+tt191
tt323 = tt1*(miu(1,1)*(tt283+tt230+tt199+tt313+tt300+tt258)+1.0E+&
&0*lam(1,1)*tt17*tt322)
tt324 = tt285+tt232+tt201
tt325 = tt1*(miu(1,1)*(tt291+tt240+tt209+tt317+tt304+tt266)+1.0E+&
&0*lam(1,1)*tt17*tt324)
tt326 = 1.0E+0*tt143**2
tt327 = 1.0E+0*tt147**2
tt328 = 1.0E+0*tt151**2
tt329 = 1.25E-1*tt18*tt27
tt330 = 1.25E-1*tt18*tt41
tt331 = 1.25E-1*tt18*tt49
tt332 = 1.0E+0*tt143*tt165
tt333 = 1.0E+0*tt147*tt169
tt334 = 1.0E+0*tt151*tt173
tt335 = 1.25E-1*tt18*tt26*tt34
tt336 = 1.25E-1*tt18*tt40*tt44
tt337 = 1.25E-1*tt18*tt48*tt50
tt338 = tt1*(miu(1,1)*(tt337+tt336+tt335+tt334+tt333+tt332)+1.0E+&
&0*lam(1,1)*tt138*tt160)
tt339 = 1.0E+0*tt184*tt143
tt340 = 1.0E+0*tt186*tt147
tt341 = tt1*(miu(1,1)*(tt155+tt154+tt199+tt152+tt340+tt339)+1.0E+&
&0*lam(1,1)*tt180*tt138)
tt342 = 1.0E+0*tt197*tt143
tt343 = 1.0E+0*tt195*tt147
tt344 = -1.25E-1*tt18*tt27
tt345 = tt1*(miu(1,1)*(tt134+tt133+tt189+tt131+tt331+tt330+tt344+&
&tt328+tt343+tt342)+tt181+1.0E+0*lam(1,1)*tt138*tt192)
tt346 = 1.0E+0*tt147*tt205
tt347 = 1.0E+0*tt143*tt207
tt348 = -1.25E-1*tt18*tt26*tt34
tt349 = tt1*(miu(1,1)*(tt337+tt336+tt348+tt334+tt347+tt346)+1.0E+&
&0*lam(1,1)*tt138*tt202)
tt350 = 1.0E+0*tt215*tt143
tt351 = 1.0E+0*tt217*tt151
tt352 = tt1*(miu(1,1)*(tt155+tt230+tt153+tt351+tt148+tt350)+1.0E+&
&0*lam(1,1)*tt212*tt138)
tt353 = 1.0E+0*tt228*tt143
tt354 = 1.0E+0*tt226*tt151
tt355 = -1.25E-1*tt18*tt41
tt356 = tt1*(miu(1,1)*(tt134+tt220+tt132+tt129+tt331+tt355+tt329+&
&tt327+tt354+tt353)+tt181+1.0E+0*lam(1,1)*tt138*tt223)
tt357 = 1.0E+0*tt151*tt236
tt358 = 1.0E+0*tt143*tt238
tt359 = -1.25E-1*tt18*tt40*tt44
tt360 = tt1*(miu(1,1)*(tt337+tt359+tt335+tt333+tt358+tt357)+1.0E+&
&0*lam(1,1)*tt138*tt233)
tt361 = 1.0E+0*tt248*tt143
tt362 = tt1*(miu(1,1)*(tt155+tt230+tt199+tt351+tt340+tt361)+1.0E+&
&0*lam(1,1)*tt242*tt138)
tt363 = 1.0E+0*tt143*tt257
tt364 = tt1*(miu(1,1)*(tt134+tt220+tt189+tt250+tt331+tt355+tt344+&
&tt363+tt354+tt343)+tt243+1.0E+0*lam(1,1)*tt138*tt252)
tt365 = 1.0E+0*tt143*tt265
tt366 = tt1*(miu(1,1)*(tt337+tt359+tt348+tt365+tt346+tt357)+1.0E+&
&0*lam(1,1)*tt138*tt260)
tt367 = 1.0E+0*tt270*tt147
tt368 = 1.0E+0*tt272*tt151
tt369 = tt1*(miu(1,1)*(tt283+tt154+tt153+tt368+tt367+tt144)+1.0E+&
&0*lam(1,1)*tt269*tt138)
tt370 = 1.0E+0*tt281*tt147
tt371 = 1.0E+0*tt279*tt151
tt372 = -1.25E-1*tt18*tt49
tt373 = tt1*(miu(1,1)*(tt275+tt133+tt132+tt118+tt372+tt330+tt329+&
&tt326+tt371+tt370)+tt181+1.0E+0*lam(1,1)*tt138*tt278)
tt374 = 1.0E+0*tt151*tt287
tt375 = 1.0E+0*tt147*tt289
tt376 = -1.25E-1*tt18*tt48*tt50
tt377 = tt1*(miu(1,1)*(tt376+tt336+tt335+tt332+tt375+tt374)+1.0E+&
&0*lam(1,1)*tt138*tt286)
tt378 = 1.0E+0*tt294*tt147
tt379 = tt1*(miu(1,1)*(tt283+tt154+tt199+tt368+tt378+tt339)+1.0E+&
&0*lam(1,1)*tt293*tt138)
tt380 = 1.0E+0*tt147*tt299
tt381 = tt1*(miu(1,1)*(tt275+tt133+tt189+tt296+tt372+tt330+tt344+&
&tt380+tt371+tt342)+tt243+1.0E+0*lam(1,1)*tt138*tt298)
tt382 = 1.0E+0*tt147*tt303
tt383 = tt1*(miu(1,1)*(tt376+tt336+tt348+tt382+tt347+tt374)+1.0E+&
&0*lam(1,1)*tt138*tt302)
tt384 = 1.0E+0*tt307*tt151
tt385 = tt1*(miu(1,1)*(tt283+tt230+tt153+tt384+tt367+tt350)+1.0E+&
&0*lam(1,1)*tt306*tt138)
tt386 = 1.0E+0*tt151*tt312
tt387 = tt1*(miu(1,1)*(tt275+tt220+tt132+tt309+tt372+tt355+tt329+&
&tt386+tt370+tt353)+tt243+1.0E+0*lam(1,1)*tt138*tt311)
tt388 = 1.0E+0*tt151*tt316
tt389 = tt1*(miu(1,1)*(tt376+tt359+tt335+tt388+tt358+tt375)+1.0E+&
&0*lam(1,1)*tt138*tt315)
tt390 = tt1*(miu(1,1)*(tt283+tt230+tt199+tt384+tt378+tt361)+1.0E+&
&0*lam(1,1)*tt319*tt138)
tt391 = tt1*(miu(1,1)*(tt275+tt220+tt189+tt309+tt296+tt250+tt372+&
&tt355+tt344+tt386+tt380+tt363)+tt320+1.0E+0*lam(1,1)*tt138*tt322)
tt392 = tt1*(miu(1,1)*(tt376+tt359+tt348+tt388+tt382+tt365)+1.0E+&
&0*lam(1,1)*tt138*tt324)
tt393 = 1.0E+0*tt165**2
tt394 = 1.0E+0*tt169**2
tt395 = 1.0E+0*tt173**2
tt396 = 1.25E-1*tt18*tt35
tt397 = 1.25E-1*tt18*tt45
tt398 = 1.25E-1*tt18*tt51
tt399 = 1.0E+0*tt184*tt165
tt400 = 1.0E+0*tt186*tt169
tt401 = tt1*(miu(1,1)*(tt177+tt176+tt209+tt174+tt400+tt399)+1.0E+&
&0*lam(1,1)*tt180*tt160)
tt402 = 1.0E+0*tt197*tt165
tt403 = 1.0E+0*tt195*tt169
tt404 = tt1*(miu(1,1)*(tt337+tt336+tt348+tt334+tt403+tt402)+1.0E+&
&0*lam(1,1)*tt192*tt160)
tt405 = 1.0E+0*tt207*tt165
tt406 = 1.0E+0*tt205*tt169
tt407 = -1.25E-1*tt18*tt35
tt408 = tt1*(miu(1,1)*(tt134+tt398+tt133+tt397+tt189+tt407+tt395+&
&tt131+tt406+tt405)+tt181+1.0E+0*lam(1,1)*tt160*tt202)
tt409 = 1.0E+0*tt215*tt165
tt410 = 1.0E+0*tt217*tt173
tt411 = tt1*(miu(1,1)*(tt177+tt240+tt175+tt410+tt170+tt409)+1.0E+&
&0*lam(1,1)*tt212*tt160)
tt412 = 1.0E+0*tt228*tt165
tt413 = 1.0E+0*tt226*tt173
tt414 = tt1*(miu(1,1)*(tt337+tt359+tt335+tt413+tt333+tt412)+1.0E+&
&0*lam(1,1)*tt223*tt160)
tt415 = 1.0E+0*tt238*tt165
tt416 = 1.0E+0*tt236*tt173
tt417 = -1.25E-1*tt18*tt45
tt418 = tt1*(miu(1,1)*(tt134+tt398+tt220+tt417+tt132+tt396+tt394+&
&tt129+tt416+tt415)+tt181+1.0E+0*lam(1,1)*tt160*tt233)
tt419 = 1.0E+0*tt248*tt165
tt420 = tt1*(miu(1,1)*(tt177+tt240+tt209+tt410+tt400+tt419)+1.0E+&
&0*lam(1,1)*tt242*tt160)
tt421 = 1.0E+0*tt257*tt165
tt422 = tt1*(miu(1,1)*(tt337+tt359+tt348+tt413+tt403+tt421)+1.0E+&
&0*lam(1,1)*tt252*tt160)
tt423 = 1.0E+0*tt165*tt265
tt424 = tt1*(miu(1,1)*(tt134+tt398+tt220+tt417+tt189+tt407+tt250+&
&tt423+tt416+tt406)+tt243+1.0E+0*lam(1,1)*tt160*tt260)
tt425 = 1.0E+0*tt270*tt169
tt426 = 1.0E+0*tt272*tt173
tt427 = tt1*(miu(1,1)*(tt291+tt176+tt175+tt426+tt425+tt166)+1.0E+&
&0*lam(1,1)*tt269*tt160)
tt428 = 1.0E+0*tt281*tt169
tt429 = 1.0E+0*tt279*tt173
tt430 = tt1*(miu(1,1)*(tt376+tt336+tt335+tt429+tt428+tt332)+1.0E+&
&0*lam(1,1)*tt278*tt160)
tt431 = 1.0E+0*tt289*tt169
tt432 = 1.0E+0*tt287*tt173
tt433 = -1.25E-1*tt18*tt51
tt434 = tt1*(miu(1,1)*(tt275+tt433+tt133+tt397+tt132+tt396+tt393+&
&tt118+tt432+tt431)+tt181+1.0E+0*lam(1,1)*tt160*tt286)
tt435 = 1.0E+0*tt294*tt169
tt436 = tt1*(miu(1,1)*(tt291+tt176+tt209+tt426+tt435+tt399)+1.0E+&
&0*lam(1,1)*tt293*tt160)
tt437 = 1.0E+0*tt299*tt169
tt438 = tt1*(miu(1,1)*(tt376+tt336+tt348+tt429+tt437+tt402)+1.0E+&
&0*lam(1,1)*tt298*tt160)
tt439 = 1.0E+0*tt169*tt303
tt440 = tt1*(miu(1,1)*(tt275+tt433+tt133+tt397+tt189+tt407+tt296+&
&tt439+tt432+tt405)+tt243+1.0E+0*lam(1,1)*tt160*tt302)
tt441 = 1.0E+0*tt307*tt173
tt442 = tt1*(miu(1,1)*(tt291+tt240+tt175+tt441+tt425+tt409)+1.0E+&
&0*lam(1,1)*tt306*tt160)
tt443 = 1.0E+0*tt312*tt173
tt444 = tt1*(miu(1,1)*(tt376+tt359+tt335+tt443+tt428+tt412)+1.0E+&
&0*lam(1,1)*tt311*tt160)
tt445 = 1.0E+0*tt173*tt316
tt446 = tt1*(miu(1,1)*(tt275+tt433+tt220+tt417+tt132+tt396+tt309+&
&tt445+tt431+tt415)+tt243+1.0E+0*lam(1,1)*tt160*tt315)
tt447 = tt1*(miu(1,1)*(tt291+tt240+tt209+tt441+tt435+tt419)+1.0E+&
&0*lam(1,1)*tt319*tt160)
tt448 = tt1*(miu(1,1)*(tt376+tt359+tt348+tt443+tt437+tt421)+1.0E+&
&0*lam(1,1)*tt322*tt160)
tt449 = tt1*(miu(1,1)*(tt275+tt433+tt220+tt417+tt189+tt407+tt309+&
&tt296+tt250+tt445+tt439+tt423)+tt320+1.0E+0*lam(1,1)*tt160*tt324)
tt450 = 1.0E+0*tt186**2
tt451 = 1.0E+0*tt184**2
tt452 = 1.0E+0*tt186*tt195
tt453 = 1.0E+0*tt184*tt197
tt454 = tt1*(miu(1,1)*(tt155+tt154+tt153+tt152+tt453+tt452)+1.0E+&
&0*lam(1,1)*tt180*tt192)
tt455 = 1.0E+0*tt186*tt205
tt456 = 1.0E+0*tt184*tt207
tt457 = tt1*(miu(1,1)*(tt177+tt176+tt175+tt174+tt456+tt455)+1.0E+&
&0*lam(1,1)*tt180*tt202)
tt458 = 1.0E+0*tt215*tt184
tt459 = tt1*(miu(1,1)*(tt134+tt220+tt189+tt118+tt71+tt219+tt188+t&
&t218+tt187+tt458)+tt243+1.0E+0*lam(1,1)*tt180*tt212)
tt460 = 1.0E+0*tt184*tt228
tt461 = tt1*(miu(1,1)*(tt155+tt230+tt199+tt340+tt460+tt227)+1.0E+&
&0*lam(1,1)*tt180*tt223)
tt462 = 1.0E+0*tt184*tt238
tt463 = tt1*(miu(1,1)*(tt177+tt240+tt209+tt400+tt462+tt237)+1.0E+&
&0*lam(1,1)*tt180*tt233)
tt464 = 1.0E+0*tt184*tt248
tt465 = tt1*(miu(1,1)*(tt134+tt220+tt132+tt296+tt71+tt219+tt69+tt&
&464+tt218+tt450)+tt181+1.0E+0*lam(1,1)*tt180*tt242)
tt466 = 1.0E+0*tt184*tt257
tt467 = tt1*(miu(1,1)*(tt155+tt230+tt153+tt466+tt452+tt227)+1.0E+&
&0*lam(1,1)*tt180*tt252)
tt468 = 1.0E+0*tt184*tt265
tt469 = tt1*(miu(1,1)*(tt177+tt240+tt175+tt468+tt455+tt237)+1.0E+&
&0*lam(1,1)*tt180*tt260)
tt470 = 1.0E+0*tt270*tt186
tt471 = tt1*(miu(1,1)*(tt275+tt133+tt189+tt129+tt274+tt70+tt188+t&
&t273+tt185+tt470)+tt243+1.0E+0*lam(1,1)*tt180*tt269)
tt472 = 1.0E+0*tt186*tt281
tt473 = tt1*(miu(1,1)*(tt283+tt154+tt199+tt339+tt472+tt280)+1.0E+&
&0*lam(1,1)*tt180*tt278)
tt474 = 1.0E+0*tt186*tt289
tt475 = tt1*(miu(1,1)*(tt291+tt176+tt209+tt399+tt474+tt288)+1.0E+&
&0*lam(1,1)*tt180*tt286)
tt476 = 1.0E+0*tt186*tt294
tt477 = tt1*(miu(1,1)*(tt275+tt133+tt132+tt250+tt274+tt70+tt69+tt&
&476+tt273+tt451)+tt181+1.0E+0*lam(1,1)*tt180*tt293)
tt478 = 1.0E+0*tt186*tt299
tt479 = tt1*(miu(1,1)*(tt283+tt154+tt153+tt478+tt453+tt280)+1.0E+&
&0*lam(1,1)*tt180*tt298)
tt480 = 1.0E+0*tt186*tt303
tt481 = tt1*(miu(1,1)*(tt291+tt176+tt175+tt480+tt456+tt288)+1.0E+&
&0*lam(1,1)*tt180*tt302)
tt482 = tt1*(miu(1,1)*(tt275+tt220+tt189+tt309+tt129+tt118+tt274+&
&tt219+tt188+tt308+tt458+tt470)+tt320+1.0E+0*lam(1,1)*tt180*tt306)
tt483 = tt1*(miu(1,1)*(tt283+tt230+tt199+tt313+tt460+tt472)+1.0E+&
&0*lam(1,1)*tt180*tt311)
tt484 = tt1*(miu(1,1)*(tt291+tt240+tt209+tt317+tt462+tt474)+1.0E+&
&0*lam(1,1)*tt180*tt315)
tt485 = tt1*(miu(1,1)*(tt275+tt220+tt132+tt309+tt274+tt219+tt69+t&
&t308+tt476+tt464)+tt243+1.0E+0*lam(1,1)*tt180*tt319)
tt486 = tt1*(miu(1,1)*(tt283+tt230+tt153+tt313+tt478+tt466)+1.0E+&
&0*lam(1,1)*tt180*tt322)
tt487 = tt1*(miu(1,1)*(tt291+tt240+tt175+tt317+tt480+tt468)+1.0E+&
&0*lam(1,1)*tt180*tt324)
tt488 = 1.0E+0*tt195**2
tt489 = 1.0E+0*tt197**2
tt490 = 1.0E+0*tt195*tt205
tt491 = 1.0E+0*tt197*tt207
tt492 = tt1*(miu(1,1)*(tt337+tt336+tt335+tt334+tt491+tt490)+1.0E+&
&0*lam(1,1)*tt192*tt202)
tt493 = 1.0E+0*tt215*tt197
tt494 = tt1*(miu(1,1)*(tt155+tt230+tt199+tt351+tt493+tt196)+1.0E+&
&0*lam(1,1)*tt212*tt192)
tt495 = 1.0E+0*tt228*tt197
tt496 = tt1*(miu(1,1)*(tt134+tt220+tt189+tt118+tt331+tt355+tt344+&
&tt354+tt343+tt495)+tt243+1.0E+0*lam(1,1)*tt192*tt223)
tt497 = 1.0E+0*tt197*tt238
tt498 = tt1*(miu(1,1)*(tt337+tt359+tt348+tt403+tt497+tt357)+1.0E+&
&0*lam(1,1)*tt192*tt233)
tt499 = 1.0E+0*tt248*tt197
tt500 = tt1*(miu(1,1)*(tt155+tt230+tt153+tt351+tt499+tt452)+1.0E+&
&0*lam(1,1)*tt242*tt192)
tt501 = 1.0E+0*tt197*tt257
tt502 = tt1*(miu(1,1)*(tt134+tt220+tt132+tt296+tt331+tt355+tt329+&
&tt501+tt354+tt488)+tt181+1.0E+0*lam(1,1)*tt192*tt252)
tt503 = 1.0E+0*tt197*tt265
tt504 = tt1*(miu(1,1)*(tt337+tt359+tt335+tt503+tt490+tt357)+1.0E+&
&0*lam(1,1)*tt192*tt260)
tt505 = 1.0E+0*tt270*tt195
tt506 = tt1*(miu(1,1)*(tt283+tt154+tt199+tt368+tt198+tt505)+1.0E+&
&0*lam(1,1)*tt269*tt192)
tt507 = 1.0E+0*tt281*tt195
tt508 = tt1*(miu(1,1)*(tt275+tt133+tt189+tt129+tt372+tt330+tt344+&
&tt371+tt342+tt507)+tt243+1.0E+0*lam(1,1)*tt192*tt278)
tt509 = 1.0E+0*tt195*tt289
tt510 = tt1*(miu(1,1)*(tt376+tt336+tt348+tt402+tt509+tt374)+1.0E+&
&0*lam(1,1)*tt192*tt286)
tt511 = 1.0E+0*tt294*tt195
tt512 = tt1*(miu(1,1)*(tt283+tt154+tt153+tt368+tt453+tt511)+1.0E+&
&0*lam(1,1)*tt293*tt192)
tt513 = 1.0E+0*tt195*tt299
tt514 = tt1*(miu(1,1)*(tt275+tt133+tt132+tt250+tt372+tt330+tt329+&
&tt513+tt371+tt489)+tt181+1.0E+0*lam(1,1)*tt192*tt298)
tt515 = 1.0E+0*tt195*tt303
tt516 = tt1*(miu(1,1)*(tt376+tt336+tt335+tt515+tt491+tt374)+1.0E+&
&0*lam(1,1)*tt192*tt302)
tt517 = tt1*(miu(1,1)*(tt283+tt230+tt199+tt384+tt493+tt505)+1.0E+&
&0*lam(1,1)*tt306*tt192)
tt518 = tt1*(miu(1,1)*(tt275+tt220+tt189+tt309+tt129+tt118+tt372+&
&tt355+tt344+tt386+tt495+tt507)+tt320+1.0E+0*lam(1,1)*tt192*tt311)
tt519 = tt1*(miu(1,1)*(tt376+tt359+tt348+tt388+tt497+tt509)+1.0E+&
&0*lam(1,1)*tt192*tt315)
tt520 = tt1*(miu(1,1)*(tt283+tt230+tt153+tt384+tt499+tt511)+1.0E+&
&0*lam(1,1)*tt319*tt192)
tt521 = tt1*(miu(1,1)*(tt275+tt220+tt132+tt309+tt372+tt355+tt329+&
&tt386+tt513+tt501)+tt243+1.0E+0*lam(1,1)*tt192*tt322)
tt522 = tt1*(miu(1,1)*(tt376+tt359+tt335+tt388+tt515+tt503)+1.0E+&
&0*lam(1,1)*tt192*tt324)
tt523 = 1.0E+0*tt205**2
tt524 = 1.0E+0*tt207**2
tt525 = 1.0E+0*tt215*tt207
tt526 = tt1*(miu(1,1)*(tt177+tt240+tt209+tt410+tt525+tt206)+1.0E+&
&0*lam(1,1)*tt212*tt202)
tt527 = 1.0E+0*tt228*tt207
tt528 = tt1*(miu(1,1)*(tt337+tt359+tt348+tt413+tt527+tt346)+1.0E+&
&0*lam(1,1)*tt223*tt202)
tt529 = 1.0E+0*tt238*tt207
tt530 = tt1*(miu(1,1)*(tt134+tt398+tt220+tt417+tt189+tt407+tt118+&
&tt416+tt406+tt529)+tt243+1.0E+0*lam(1,1)*tt202*tt233)
tt531 = 1.0E+0*tt248*tt207
tt532 = tt1*(miu(1,1)*(tt177+tt240+tt175+tt410+tt531+tt455)+1.0E+&
&0*lam(1,1)*tt242*tt202)
tt533 = 1.0E+0*tt257*tt207
tt534 = tt1*(miu(1,1)*(tt337+tt359+tt335+tt413+tt533+tt490)+1.0E+&
&0*lam(1,1)*tt252*tt202)
tt535 = 1.0E+0*tt207*tt265
tt536 = tt1*(miu(1,1)*(tt134+tt398+tt220+tt417+tt132+tt396+tt296+&
&tt535+tt416+tt523)+tt181+1.0E+0*lam(1,1)*tt202*tt260)
tt537 = 1.0E+0*tt270*tt205
tt538 = tt1*(miu(1,1)*(tt291+tt176+tt209+tt426+tt208+tt537)+1.0E+&
&0*lam(1,1)*tt269*tt202)
tt539 = 1.0E+0*tt281*tt205
tt540 = tt1*(miu(1,1)*(tt376+tt336+tt348+tt429+tt347+tt539)+1.0E+&
&0*lam(1,1)*tt278*tt202)
tt541 = 1.0E+0*tt289*tt205
tt542 = tt1*(miu(1,1)*(tt275+tt433+tt133+tt397+tt189+tt407+tt129+&
&tt432+tt405+tt541)+tt243+1.0E+0*lam(1,1)*tt202*tt286)
tt543 = 1.0E+0*tt294*tt205
tt544 = tt1*(miu(1,1)*(tt291+tt176+tt175+tt426+tt456+tt543)+1.0E+&
&0*lam(1,1)*tt293*tt202)
tt545 = 1.0E+0*tt299*tt205
tt546 = tt1*(miu(1,1)*(tt376+tt336+tt335+tt429+tt491+tt545)+1.0E+&
&0*lam(1,1)*tt298*tt202)
tt547 = 1.0E+0*tt205*tt303
tt548 = tt1*(miu(1,1)*(tt275+tt433+tt133+tt397+tt132+tt396+tt250+&
&tt547+tt432+tt524)+tt181+1.0E+0*lam(1,1)*tt202*tt302)
tt549 = tt1*(miu(1,1)*(tt291+tt240+tt209+tt441+tt525+tt537)+1.0E+&
&0*lam(1,1)*tt306*tt202)
tt550 = tt1*(miu(1,1)*(tt376+tt359+tt348+tt443+tt527+tt539)+1.0E+&
&0*lam(1,1)*tt311*tt202)
tt551 = tt1*(miu(1,1)*(tt275+tt433+tt220+tt417+tt189+tt407+tt309+&
&tt129+tt118+tt445+tt529+tt541)+tt320+1.0E+0*lam(1,1)*tt202*tt315)
tt552 = tt1*(miu(1,1)*(tt291+tt240+tt175+tt441+tt531+tt543)+1.0E+&
&0*lam(1,1)*tt319*tt202)
tt553 = tt1*(miu(1,1)*(tt376+tt359+tt335+tt443+tt533+tt545)+1.0E+&
&0*lam(1,1)*tt322*tt202)
tt554 = tt1*(miu(1,1)*(tt275+tt433+tt220+tt417+tt132+tt396+tt309+&
&tt445+tt547+tt535)+tt243+1.0E+0*lam(1,1)*tt202*tt324)
tt555 = 1.0E+0*tt217**2
tt556 = 1.0E+0*tt215**2
tt557 = 1.0E+0*tt217*tt226
tt558 = 1.0E+0*tt215*tt228
tt559 = tt1*(miu(1,1)*(tt155+tt154+tt153+tt148+tt558+tt557)+1.0E+&
&0*lam(1,1)*tt212*tt223)
tt560 = 1.0E+0*tt217*tt236
tt561 = 1.0E+0*tt215*tt238
tt562 = tt1*(miu(1,1)*(tt177+tt176+tt175+tt170+tt561+tt560)+1.0E+&
&0*lam(1,1)*tt212*tt233)
tt563 = 1.0E+0*tt215*tt248
tt564 = tt1*(miu(1,1)*(tt134+tt133+tt189+tt309+tt71+tt70+tt188+tt&
&563+tt187+tt555)+tt181+1.0E+0*lam(1,1)*tt212*tt242)
tt565 = 1.0E+0*tt215*tt257
tt566 = tt1*(miu(1,1)*(tt155+tt154+tt199+tt565+tt196+tt557)+1.0E+&
&0*lam(1,1)*tt212*tt252)
tt567 = 1.0E+0*tt215*tt265
tt568 = tt1*(miu(1,1)*(tt177+tt176+tt209+tt567+tt206+tt560)+1.0E+&
&0*lam(1,1)*tt212*tt260)
tt569 = 1.0E+0*tt272*tt217
tt570 = tt1*(miu(1,1)*(tt275+tt220+tt132+tt131+tt274+tt219+tt69+t&
&t271+tt216+tt569)+tt243+1.0E+0*lam(1,1)*tt212*tt269)
tt571 = 1.0E+0*tt217*tt279
tt572 = tt1*(miu(1,1)*(tt283+tt230+tt153+tt350+tt282+tt571)+1.0E+&
&0*lam(1,1)*tt212*tt278)
tt573 = 1.0E+0*tt217*tt287
tt574 = tt1*(miu(1,1)*(tt291+tt240+tt175+tt409+tt290+tt573)+1.0E+&
&0*lam(1,1)*tt212*tt286)
tt575 = tt1*(miu(1,1)*(tt275+tt220+tt189+tt131+tt296+tt118+tt274+&
&tt219+tt188+tt295+tt458+tt569)+tt320+1.0E+0*lam(1,1)*tt212*tt293)
tt576 = tt1*(miu(1,1)*(tt283+tt230+tt199+tt300+tt493+tt571)+1.0E+&
&0*lam(1,1)*tt212*tt298)
tt577 = tt1*(miu(1,1)*(tt291+tt240+tt209+tt304+tt525+tt573)+1.0E+&
&0*lam(1,1)*tt212*tt302)
tt578 = 1.0E+0*tt217*tt307
tt579 = tt1*(miu(1,1)*(tt275+tt133+tt132+tt250+tt274+tt70+tt69+tt&
&578+tt271+tt556)+tt181+1.0E+0*lam(1,1)*tt212*tt306)
tt580 = 1.0E+0*tt217*tt312
tt581 = tt1*(miu(1,1)*(tt283+tt154+tt153+tt580+tt558+tt282)+1.0E+&
&0*lam(1,1)*tt212*tt311)
tt582 = 1.0E+0*tt217*tt316
tt583 = tt1*(miu(1,1)*(tt291+tt176+tt175+tt582+tt561+tt290)+1.0E+&
&0*lam(1,1)*tt212*tt315)
tt584 = tt1*(miu(1,1)*(tt275+tt133+tt189+tt296+tt274+tt70+tt188+t&
&t578+tt295+tt563)+tt243+1.0E+0*lam(1,1)*tt212*tt319)
tt585 = tt1*(miu(1,1)*(tt283+tt154+tt199+tt580+tt300+tt565)+1.0E+&
&0*lam(1,1)*tt212*tt322)
tt586 = tt1*(miu(1,1)*(tt291+tt176+tt209+tt582+tt304+tt567)+1.0E+&
&0*lam(1,1)*tt212*tt324)
tt587 = 1.0E+0*tt226**2
tt588 = 1.0E+0*tt228**2
tt589 = 1.0E+0*tt226*tt236
tt590 = 1.0E+0*tt228*tt238
tt591 = tt1*(miu(1,1)*(tt337+tt336+tt335+tt333+tt590+tt589)+1.0E+&
&0*lam(1,1)*tt223*tt233)
tt592 = 1.0E+0*tt248*tt228
tt593 = tt1*(miu(1,1)*(tt155+tt154+tt199+tt340+tt592+tt557)+1.0E+&
&0*lam(1,1)*tt242*tt223)
tt594 = 1.0E+0*tt228*tt257
tt595 = tt1*(miu(1,1)*(tt134+tt133+tt189+tt309+tt331+tt330+tt344+&
&tt594+tt343+tt587)+tt181+1.0E+0*lam(1,1)*tt223*tt252)
tt596 = 1.0E+0*tt228*tt265
tt597 = tt1*(miu(1,1)*(tt337+tt336+tt348+tt596+tt346+tt589)+1.0E+&
&0*lam(1,1)*tt223*tt260)
tt598 = 1.0E+0*tt272*tt226
tt599 = tt1*(miu(1,1)*(tt283+tt230+tt153+tt367+tt229+tt598)+1.0E+&
&0*lam(1,1)*tt269*tt223)
tt600 = 1.0E+0*tt279*tt226
tt601 = tt1*(miu(1,1)*(tt275+tt220+tt132+tt131+tt372+tt355+tt329+&
&tt370+tt353+tt600)+tt243+1.0E+0*lam(1,1)*tt223*tt278)
tt602 = 1.0E+0*tt226*tt287
tt603 = tt1*(miu(1,1)*(tt376+tt359+tt335+tt412+tt375+tt602)+1.0E+&
&0*lam(1,1)*tt223*tt286)
tt604 = tt1*(miu(1,1)*(tt283+tt230+tt199+tt378+tt460+tt598)+1.0E+&
&0*lam(1,1)*tt293*tt223)
tt605 = tt1*(miu(1,1)*(tt275+tt220+tt189+tt131+tt296+tt118+tt372+&
&tt355+tt344+tt380+tt495+tt600)+tt320+1.0E+0*lam(1,1)*tt223*tt298)
tt606 = tt1*(miu(1,1)*(tt376+tt359+tt348+tt382+tt527+tt602)+1.0E+&
&0*lam(1,1)*tt223*tt302)
tt607 = 1.0E+0*tt307*tt226
tt608 = tt1*(miu(1,1)*(tt283+tt154+tt153+tt367+tt558+tt607)+1.0E+&
&0*lam(1,1)*tt306*tt223)
tt609 = 1.0E+0*tt226*tt312
tt610 = tt1*(miu(1,1)*(tt275+tt133+tt132+tt250+tt372+tt330+tt329+&
&tt609+tt370+tt588)+tt181+1.0E+0*lam(1,1)*tt223*tt311)
tt611 = 1.0E+0*tt226*tt316
tt612 = tt1*(miu(1,1)*(tt376+tt336+tt335+tt611+tt590+tt375)+1.0E+&
&0*lam(1,1)*tt223*tt315)
tt613 = tt1*(miu(1,1)*(tt283+tt154+tt199+tt378+tt592+tt607)+1.0E+&
&0*lam(1,1)*tt319*tt223)
tt614 = tt1*(miu(1,1)*(tt275+tt133+tt189+tt296+tt372+tt330+tt344+&
&tt609+tt380+tt594)+tt243+1.0E+0*lam(1,1)*tt223*tt322)
tt615 = tt1*(miu(1,1)*(tt376+tt336+tt348+tt611+tt382+tt596)+1.0E+&
&0*lam(1,1)*tt223*tt324)
tt616 = 1.0E+0*tt236**2
tt617 = 1.0E+0*tt238**2
tt618 = 1.0E+0*tt248*tt238
tt619 = tt1*(miu(1,1)*(tt177+tt176+tt209+tt400+tt618+tt560)+1.0E+&
&0*lam(1,1)*tt242*tt233)
tt620 = 1.0E+0*tt257*tt238
tt621 = tt1*(miu(1,1)*(tt337+tt336+tt348+tt403+tt620+tt589)+1.0E+&
&0*lam(1,1)*tt252*tt233)
tt622 = 1.0E+0*tt238*tt265
tt623 = tt1*(miu(1,1)*(tt134+tt398+tt133+tt397+tt189+tt407+tt309+&
&tt622+tt406+tt616)+tt181+1.0E+0*lam(1,1)*tt233*tt260)
tt624 = 1.0E+0*tt272*tt236
tt625 = tt1*(miu(1,1)*(tt291+tt240+tt175+tt425+tt239+tt624)+1.0E+&
&0*lam(1,1)*tt269*tt233)
tt626 = 1.0E+0*tt279*tt236
tt627 = tt1*(miu(1,1)*(tt376+tt359+tt335+tt428+tt358+tt626)+1.0E+&
&0*lam(1,1)*tt278*tt233)
tt628 = 1.0E+0*tt287*tt236
tt629 = tt1*(miu(1,1)*(tt275+tt433+tt220+tt417+tt132+tt396+tt131+&
&tt431+tt415+tt628)+tt243+1.0E+0*lam(1,1)*tt233*tt286)
tt630 = tt1*(miu(1,1)*(tt291+tt240+tt209+tt435+tt462+tt624)+1.0E+&
&0*lam(1,1)*tt293*tt233)
tt631 = tt1*(miu(1,1)*(tt376+tt359+tt348+tt437+tt497+tt626)+1.0E+&
&0*lam(1,1)*tt298*tt233)
tt632 = tt1*(miu(1,1)*(tt275+tt433+tt220+tt417+tt189+tt407+tt131+&
&tt296+tt118+tt439+tt529+tt628)+tt320+1.0E+0*lam(1,1)*tt233*tt302)
tt633 = 1.0E+0*tt307*tt236
tt634 = tt1*(miu(1,1)*(tt291+tt176+tt175+tt425+tt561+tt633)+1.0E+&
&0*lam(1,1)*tt306*tt233)
tt635 = 1.0E+0*tt312*tt236
tt636 = tt1*(miu(1,1)*(tt376+tt336+tt335+tt428+tt590+tt635)+1.0E+&
&0*lam(1,1)*tt311*tt233)
tt637 = 1.0E+0*tt236*tt316
tt638 = tt1*(miu(1,1)*(tt275+tt433+tt133+tt397+tt132+tt396+tt250+&
&tt637+tt431+tt617)+tt181+1.0E+0*lam(1,1)*tt233*tt315)
tt639 = tt1*(miu(1,1)*(tt291+tt176+tt209+tt435+tt618+tt633)+1.0E+&
&0*lam(1,1)*tt319*tt233)
tt640 = tt1*(miu(1,1)*(tt376+tt336+tt348+tt437+tt620+tt635)+1.0E+&
&0*lam(1,1)*tt322*tt233)
tt641 = tt1*(miu(1,1)*(tt275+tt433+tt133+tt397+tt189+tt407+tt296+&
&tt637+tt439+tt622)+tt243+1.0E+0*lam(1,1)*tt233*tt324)
tt642 = 1.0E+0*tt248**2
tt643 = 1.0E+0*tt248*tt257
tt644 = tt1*(miu(1,1)*(tt155+tt154+tt153+tt643+tt452+tt557)+1.0E+&
&0*lam(1,1)*tt242*tt252)
tt645 = 1.0E+0*tt248*tt265
tt646 = tt1*(miu(1,1)*(tt177+tt176+tt175+tt645+tt455+tt560)+1.0E+&
&0*lam(1,1)*tt242*tt260)
tt647 = tt1*(miu(1,1)*(tt275+tt220+tt189+tt131+tt129+tt250+tt274+&
&tt219+tt188+tt249+tt470+tt569)+tt320+1.0E+0*lam(1,1)*tt242*tt269)
tt648 = tt1*(miu(1,1)*(tt283+tt230+tt199+tt361+tt472+tt571)+1.0E+&
&0*lam(1,1)*tt242*tt278)
tt649 = tt1*(miu(1,1)*(tt291+tt240+tt209+tt419+tt474+tt573)+1.0E+&
&0*lam(1,1)*tt242*tt286)
tt650 = tt1*(miu(1,1)*(tt275+tt220+tt132+tt131+tt274+tt219+tt69+t&
&t476+tt464+tt569)+tt243+1.0E+0*lam(1,1)*tt242*tt293)
tt651 = tt1*(miu(1,1)*(tt283+tt230+tt153+tt478+tt499+tt571)+1.0E+&
&0*lam(1,1)*tt242*tt298)
tt652 = tt1*(miu(1,1)*(tt291+tt240+tt175+tt480+tt531+tt573)+1.0E+&
&0*lam(1,1)*tt242*tt302)
tt653 = tt1*(miu(1,1)*(tt275+tt133+tt189+tt129+tt274+tt70+tt188+t&
&t578+tt563+tt470)+tt243+1.0E+0*lam(1,1)*tt242*tt306)
tt654 = tt1*(miu(1,1)*(tt283+tt154+tt199+tt580+tt592+tt472)+1.0E+&
&0*lam(1,1)*tt242*tt311)
tt655 = tt1*(miu(1,1)*(tt291+tt176+tt209+tt582+tt618+tt474)+1.0E+&
&0*lam(1,1)*tt242*tt315)
tt656 = tt1*(miu(1,1)*(tt275+tt133+tt132+tt118+tt274+tt70+tt69+tt&
&642+tt578+tt476)+tt181+1.0E+0*lam(1,1)*tt242*tt319)
tt657 = tt1*(miu(1,1)*(tt283+tt154+tt153+tt580+tt478+tt643)+1.0E+&
&0*lam(1,1)*tt242*tt322)
tt658 = tt1*(miu(1,1)*(tt291+tt176+tt175+tt582+tt480+tt645)+1.0E+&
&0*lam(1,1)*tt242*tt324)
tt659 = 1.0E+0*tt257**2
tt660 = 1.0E+0*tt257*tt265
tt661 = tt1*(miu(1,1)*(tt337+tt336+tt335+tt660+tt490+tt589)+1.0E+&
&0*lam(1,1)*tt252*tt260)
tt662 = tt1*(miu(1,1)*(tt283+tt230+tt199+tt258+tt505+tt598)+1.0E+&
&0*lam(1,1)*tt269*tt252)
tt663 = tt1*(miu(1,1)*(tt275+tt220+tt189+tt131+tt129+tt250+tt372+&
&tt355+tt344+tt363+tt507+tt600)+tt320+1.0E+0*lam(1,1)*tt252*tt278)
tt664 = tt1*(miu(1,1)*(tt376+tt359+tt348+tt421+tt509+tt602)+1.0E+&
&0*lam(1,1)*tt252*tt286)
tt665 = tt1*(miu(1,1)*(tt283+tt230+tt153+tt466+tt511+tt598)+1.0E+&
&0*lam(1,1)*tt293*tt252)
tt666 = tt1*(miu(1,1)*(tt275+tt220+tt132+tt131+tt372+tt355+tt329+&
&tt513+tt501+tt600)+tt243+1.0E+0*lam(1,1)*tt252*tt298)
tt667 = tt1*(miu(1,1)*(tt376+tt359+tt335+tt515+tt533+tt602)+1.0E+&
&0*lam(1,1)*tt252*tt302)
tt668 = tt1*(miu(1,1)*(tt283+tt154+tt199+tt565+tt505+tt607)+1.0E+&
&0*lam(1,1)*tt306*tt252)
tt669 = tt1*(miu(1,1)*(tt275+tt133+tt189+tt129+tt372+tt330+tt344+&
&tt609+tt594+tt507)+tt243+1.0E+0*lam(1,1)*tt252*tt311)
tt670 = tt1*(miu(1,1)*(tt376+tt336+tt348+tt611+tt620+tt509)+1.0E+&
&0*lam(1,1)*tt252*tt315)
tt671 = tt1*(miu(1,1)*(tt283+tt154+tt153+tt643+tt511+tt607)+1.0E+&
&0*lam(1,1)*tt319*tt252)
tt672 = tt1*(miu(1,1)*(tt275+tt133+tt132+tt118+tt372+tt330+tt329+&
&tt659+tt609+tt513)+tt181+1.0E+0*lam(1,1)*tt252*tt322)
tt673 = tt1*(miu(1,1)*(tt376+tt336+tt335+tt611+tt515+tt660)+1.0E+&
&0*lam(1,1)*tt252*tt324)
tt674 = 1.0E+0*tt265**2
tt675 = tt1*(miu(1,1)*(tt291+tt240+tt209+tt266+tt537+tt624)+1.0E+&
&0*lam(1,1)*tt269*tt260)
tt676 = tt1*(miu(1,1)*(tt376+tt359+tt348+tt365+tt539+tt626)+1.0E+&
&0*lam(1,1)*tt278*tt260)
tt677 = tt1*(miu(1,1)*(tt275+tt433+tt220+tt417+tt189+tt407+tt131+&
&tt129+tt250+tt423+tt541+tt628)+tt320+1.0E+0*lam(1,1)*tt260*tt286)
tt678 = tt1*(miu(1,1)*(tt291+tt240+tt175+tt468+tt543+tt624)+1.0E+&
&0*lam(1,1)*tt293*tt260)
tt679 = tt1*(miu(1,1)*(tt376+tt359+tt335+tt503+tt545+tt626)+1.0E+&
&0*lam(1,1)*tt298*tt260)
tt680 = tt1*(miu(1,1)*(tt275+tt433+tt220+tt417+tt132+tt396+tt131+&
&tt547+tt535+tt628)+tt243+1.0E+0*lam(1,1)*tt260*tt302)
tt681 = tt1*(miu(1,1)*(tt291+tt176+tt209+tt567+tt537+tt633)+1.0E+&
&0*lam(1,1)*tt306*tt260)
tt682 = tt1*(miu(1,1)*(tt376+tt336+tt348+tt596+tt539+tt635)+1.0E+&
&0*lam(1,1)*tt311*tt260)
tt683 = tt1*(miu(1,1)*(tt275+tt433+tt133+tt397+tt189+tt407+tt129+&
&tt637+tt622+tt541)+tt243+1.0E+0*lam(1,1)*tt260*tt315)
tt684 = tt1*(miu(1,1)*(tt291+tt176+tt175+tt645+tt543+tt633)+1.0E+&
&0*lam(1,1)*tt319*tt260)
tt685 = tt1*(miu(1,1)*(tt376+tt336+tt335+tt660+tt545+tt635)+1.0E+&
&0*lam(1,1)*tt322*tt260)
tt686 = tt1*(miu(1,1)*(tt275+tt433+tt133+tt397+tt132+tt396+tt674+&
&tt118+tt637+tt547)+tt181+1.0E+0*lam(1,1)*tt260*tt324)
tt687 = 1.0E+0*tt272**2
tt688 = 1.0E+0*tt270**2
tt689 = 1.0E+0*tt272*tt279
tt690 = 1.0E+0*tt270*tt281
tt691 = tt1*(miu(1,1)*(tt155+tt154+tt153+tt144+tt690+tt689)+1.0E+&
&0*lam(1,1)*tt269*tt278)
tt692 = 1.0E+0*tt272*tt287
tt693 = 1.0E+0*tt270*tt289
tt694 = tt1*(miu(1,1)*(tt177+tt176+tt175+tt166+tt693+tt692)+1.0E+&
&0*lam(1,1)*tt269*tt286)
tt695 = 1.0E+0*tt270*tt294
tt696 = tt1*(miu(1,1)*(tt134+tt133+tt189+tt309+tt71+tt70+tt188+tt&
&695+tt185+tt687)+tt181+1.0E+0*lam(1,1)*tt269*tt293)
tt697 = 1.0E+0*tt270*tt299
tt698 = tt1*(miu(1,1)*(tt155+tt154+tt199+tt697+tt198+tt689)+1.0E+&
&0*lam(1,1)*tt269*tt298)
tt699 = 1.0E+0*tt270*tt303
tt700 = tt1*(miu(1,1)*(tt177+tt176+tt209+tt699+tt208+tt692)+1.0E+&
&0*lam(1,1)*tt269*tt302)
tt701 = 1.0E+0*tt272*tt307
tt702 = tt1*(miu(1,1)*(tt134+tt220+tt132+tt296+tt71+tt219+tt69+tt&
&701+tt216+tt688)+tt181+1.0E+0*lam(1,1)*tt269*tt306)
tt703 = 1.0E+0*tt272*tt312
tt704 = tt1*(miu(1,1)*(tt155+tt230+tt153+tt703+tt229+tt690)+1.0E+&
&0*lam(1,1)*tt269*tt311)
tt705 = 1.0E+0*tt272*tt316
tt706 = tt1*(miu(1,1)*(tt177+tt240+tt175+tt705+tt239+tt693)+1.0E+&
&0*lam(1,1)*tt269*tt315)
tt707 = tt1*(miu(1,1)*(tt134+tt220+tt189+tt250+tt71+tt219+tt188+t&
&t701+tt695+tt249)+tt243+1.0E+0*lam(1,1)*tt269*tt319)
tt708 = tt1*(miu(1,1)*(tt155+tt230+tt199+tt703+tt697+tt258)+1.0E+&
&0*lam(1,1)*tt269*tt322)
tt709 = tt1*(miu(1,1)*(tt177+tt240+tt209+tt705+tt699+tt266)+1.0E+&
&0*lam(1,1)*tt269*tt324)
tt710 = 1.0E+0*tt279**2
tt711 = 1.0E+0*tt281**2
tt712 = 1.0E+0*tt279*tt287
tt713 = 1.0E+0*tt281*tt289
tt714 = tt1*(miu(1,1)*(tt337+tt336+tt335+tt332+tt713+tt712)+1.0E+&
&0*lam(1,1)*tt278*tt286)
tt715 = 1.0E+0*tt294*tt281
tt716 = tt1*(miu(1,1)*(tt155+tt154+tt199+tt339+tt715+tt689)+1.0E+&
&0*lam(1,1)*tt293*tt278)
tt717 = 1.0E+0*tt281*tt299
tt718 = tt1*(miu(1,1)*(tt134+tt133+tt189+tt309+tt331+tt330+tt344+&
&tt717+tt342+tt710)+tt181+1.0E+0*lam(1,1)*tt278*tt298)
tt719 = 1.0E+0*tt281*tt303
tt720 = tt1*(miu(1,1)*(tt337+tt336+tt348+tt719+tt347+tt712)+1.0E+&
&0*lam(1,1)*tt278*tt302)
tt721 = 1.0E+0*tt307*tt279
tt722 = tt1*(miu(1,1)*(tt155+tt230+tt153+tt350+tt690+tt721)+1.0E+&
&0*lam(1,1)*tt306*tt278)
tt723 = 1.0E+0*tt279*tt312
tt724 = tt1*(miu(1,1)*(tt134+tt220+tt132+tt296+tt331+tt355+tt329+&
&tt723+tt353+tt711)+tt181+1.0E+0*lam(1,1)*tt278*tt311)
tt725 = 1.0E+0*tt279*tt316
tt726 = tt1*(miu(1,1)*(tt337+tt359+tt335+tt725+tt358+tt713)+1.0E+&
&0*lam(1,1)*tt278*tt315)
tt727 = tt1*(miu(1,1)*(tt155+tt230+tt199+tt361+tt715+tt721)+1.0E+&
&0*lam(1,1)*tt319*tt278)
tt728 = tt1*(miu(1,1)*(tt134+tt220+tt189+tt250+tt331+tt355+tt344+&
&tt723+tt717+tt363)+tt243+1.0E+0*lam(1,1)*tt278*tt322)
tt729 = tt1*(miu(1,1)*(tt337+tt359+tt348+tt725+tt719+tt365)+1.0E+&
&0*lam(1,1)*tt278*tt324)
tt730 = 1.0E+0*tt287**2
tt731 = 1.0E+0*tt289**2
tt732 = 1.0E+0*tt294*tt289
tt733 = tt1*(miu(1,1)*(tt177+tt176+tt209+tt399+tt732+tt692)+1.0E+&
&0*lam(1,1)*tt293*tt286)
tt734 = 1.0E+0*tt299*tt289
tt735 = tt1*(miu(1,1)*(tt337+tt336+tt348+tt402+tt734+tt712)+1.0E+&
&0*lam(1,1)*tt298*tt286)
tt736 = 1.0E+0*tt289*tt303
tt737 = tt1*(miu(1,1)*(tt134+tt398+tt133+tt397+tt189+tt407+tt309+&
&tt736+tt405+tt730)+tt181+1.0E+0*lam(1,1)*tt286*tt302)
tt738 = 1.0E+0*tt307*tt287
tt739 = tt1*(miu(1,1)*(tt177+tt240+tt175+tt409+tt693+tt738)+1.0E+&
&0*lam(1,1)*tt306*tt286)
tt740 = 1.0E+0*tt312*tt287
tt741 = tt1*(miu(1,1)*(tt337+tt359+tt335+tt412+tt713+tt740)+1.0E+&
&0*lam(1,1)*tt311*tt286)
tt742 = 1.0E+0*tt287*tt316
tt743 = tt1*(miu(1,1)*(tt134+tt398+tt220+tt417+tt132+tt396+tt296+&
&tt742+tt415+tt731)+tt181+1.0E+0*lam(1,1)*tt286*tt315)
tt744 = tt1*(miu(1,1)*(tt177+tt240+tt209+tt419+tt732+tt738)+1.0E+&
&0*lam(1,1)*tt319*tt286)
tt745 = tt1*(miu(1,1)*(tt337+tt359+tt348+tt421+tt734+tt740)+1.0E+&
&0*lam(1,1)*tt322*tt286)
tt746 = tt1*(miu(1,1)*(tt134+tt398+tt220+tt417+tt189+tt407+tt250+&
&tt742+tt736+tt423)+tt243+1.0E+0*lam(1,1)*tt286*tt324)
tt747 = 1.0E+0*tt294**2
tt748 = 1.0E+0*tt294*tt299
tt749 = tt1*(miu(1,1)*(tt155+tt154+tt153+tt748+tt453+tt689)+1.0E+&
&0*lam(1,1)*tt293*tt298)
tt750 = 1.0E+0*tt294*tt303
tt751 = tt1*(miu(1,1)*(tt177+tt176+tt175+tt750+tt456+tt692)+1.0E+&
&0*lam(1,1)*tt293*tt302)
tt752 = tt1*(miu(1,1)*(tt134+tt220+tt189+tt118+tt71+tt219+tt188+t&
&t701+tt695+tt458)+tt243+1.0E+0*lam(1,1)*tt293*tt306)
tt753 = tt1*(miu(1,1)*(tt155+tt230+tt199+tt703+tt460+tt715)+1.0E+&
&0*lam(1,1)*tt293*tt311)
tt754 = tt1*(miu(1,1)*(tt177+tt240+tt209+tt705+tt462+tt732)+1.0E+&
&0*lam(1,1)*tt293*tt315)
tt755 = tt1*(miu(1,1)*(tt134+tt220+tt132+tt129+tt71+tt219+tt69+tt&
&747+tt701+tt464)+tt181+1.0E+0*lam(1,1)*tt293*tt319)
tt756 = tt1*(miu(1,1)*(tt155+tt230+tt153+tt703+tt748+tt466)+1.0E+&
&0*lam(1,1)*tt293*tt322)
tt757 = tt1*(miu(1,1)*(tt177+tt240+tt175+tt705+tt750+tt468)+1.0E+&
&0*lam(1,1)*tt293*tt324)
tt758 = 1.0E+0*tt299**2
tt759 = 1.0E+0*tt299*tt303
tt760 = tt1*(miu(1,1)*(tt337+tt336+tt335+tt759+tt491+tt712)+1.0E+&
&0*lam(1,1)*tt298*tt302)
tt761 = tt1*(miu(1,1)*(tt155+tt230+tt199+tt697+tt493+tt721)+1.0E+&
&0*lam(1,1)*tt306*tt298)
tt762 = tt1*(miu(1,1)*(tt134+tt220+tt189+tt118+tt331+tt355+tt344+&
&tt723+tt717+tt495)+tt243+1.0E+0*lam(1,1)*tt298*tt311)
tt763 = tt1*(miu(1,1)*(tt337+tt359+tt348+tt725+tt497+tt734)+1.0E+&
&0*lam(1,1)*tt298*tt315)
tt764 = tt1*(miu(1,1)*(tt155+tt230+tt153+tt748+tt499+tt721)+1.0E+&
&0*lam(1,1)*tt319*tt298)
tt765 = tt1*(miu(1,1)*(tt134+tt220+tt132+tt129+tt331+tt355+tt329+&
&tt758+tt723+tt501)+tt181+1.0E+0*lam(1,1)*tt298*tt322)
tt766 = tt1*(miu(1,1)*(tt337+tt359+tt335+tt725+tt759+tt503)+1.0E+&
&0*lam(1,1)*tt298*tt324)
tt767 = 1.0E+0*tt303**2
tt768 = tt1*(miu(1,1)*(tt177+tt240+tt209+tt699+tt525+tt738)+1.0E+&
&0*lam(1,1)*tt306*tt302)
tt769 = tt1*(miu(1,1)*(tt337+tt359+tt348+tt719+tt527+tt740)+1.0E+&
&0*lam(1,1)*tt311*tt302)
tt770 = tt1*(miu(1,1)*(tt134+tt398+tt220+tt417+tt189+tt407+tt118+&
&tt742+tt736+tt529)+tt243+1.0E+0*lam(1,1)*tt302*tt315)
tt771 = tt1*(miu(1,1)*(tt177+tt240+tt175+tt750+tt531+tt738)+1.0E+&
&0*lam(1,1)*tt319*tt302)
tt772 = tt1*(miu(1,1)*(tt337+tt359+tt335+tt759+tt533+tt740)+1.0E+&
&0*lam(1,1)*tt322*tt302)
tt773 = tt1*(miu(1,1)*(tt134+tt398+tt220+tt417+tt132+tt396+tt767+&
&tt129+tt742+tt535)+tt181+1.0E+0*lam(1,1)*tt302*tt324)
tt774 = 1.0E+0*tt307**2
tt775 = 1.0E+0*tt307*tt312
tt776 = tt1*(miu(1,1)*(tt155+tt154+tt153+tt775+tt558+tt690)+1.0E+&
&0*lam(1,1)*tt306*tt311)
tt777 = 1.0E+0*tt307*tt316
tt778 = tt1*(miu(1,1)*(tt177+tt176+tt175+tt777+tt561+tt693)+1.0E+&
&0*lam(1,1)*tt306*tt315)
tt779 = tt1*(miu(1,1)*(tt134+tt133+tt189+tt131+tt71+tt70+tt188+tt&
&774+tt695+tt563)+tt181+1.0E+0*lam(1,1)*tt306*tt319)
tt780 = tt1*(miu(1,1)*(tt155+tt154+tt199+tt775+tt697+tt565)+1.0E+&
&0*lam(1,1)*tt306*tt322)
tt781 = tt1*(miu(1,1)*(tt177+tt176+tt209+tt777+tt699+tt567)+1.0E+&
&0*lam(1,1)*tt306*tt324)
tt782 = 1.0E+0*tt312**2
tt783 = 1.0E+0*tt312*tt316
tt784 = tt1*(miu(1,1)*(tt337+tt336+tt335+tt783+tt590+tt713)+1.0E+&
&0*lam(1,1)*tt311*tt315)
tt785 = tt1*(miu(1,1)*(tt155+tt154+tt199+tt775+tt592+tt715)+1.0E+&
&0*lam(1,1)*tt319*tt311)
tt786 = tt1*(miu(1,1)*(tt134+tt133+tt189+tt131+tt331+tt330+tt344+&
&tt782+tt717+tt594)+tt181+1.0E+0*lam(1,1)*tt311*tt322)
tt787 = tt1*(miu(1,1)*(tt337+tt336+tt348+tt783+tt719+tt596)+1.0E+&
&0*lam(1,1)*tt311*tt324)
tt788 = 1.0E+0*tt316**2
tt789 = tt1*(miu(1,1)*(tt177+tt176+tt209+tt777+tt618+tt732)+1.0E+&
&0*lam(1,1)*tt319*tt315)
tt790 = tt1*(miu(1,1)*(tt337+tt336+tt348+tt783+tt620+tt734)+1.0E+&
&0*lam(1,1)*tt322*tt315)
tt791 = tt1*(miu(1,1)*(tt134+tt398+tt133+tt397+tt189+tt407+tt788+&
&tt131+tt736+tt622)+tt181+1.0E+0*lam(1,1)*tt315*tt324)
tt792 = tt1*(miu(1,1)*(tt155+tt154+tt153+tt775+tt748+tt643)+1.0E+&
&0*lam(1,1)*tt319*tt322)
tt793 = tt1*(miu(1,1)*(tt177+tt176+tt175+tt777+tt750+tt645)+1.0E+&
&0*lam(1,1)*tt319*tt324)
tt794 = tt1*(miu(1,1)*(tt337+tt336+tt335+tt783+tt759+tt660)+1.0E+&
&0*lam(1,1)*tt322*tt324)
hes(1,1) = tt1*(miu(1,1)*(tt134+tt133+tt132+tt131+tt129+tt118+tt7&
&1+tt70+tt69+tt68+tt64+tt60)+tt54+1.0E+0*lam(1,1)*tt17**2)
hes(1,2) = tt156
hes(1,3) = tt178
hes(1,4) = tt190
hes(1,5) = tt200
hes(1,6) = tt210
hes(1,7) = tt221
hes(1,8) = tt231
hes(1,9) = tt241
hes(1,10) = tt251
hes(1,11) = tt259
hes(1,12) = tt267
hes(1,13) = tt276
hes(1,14) = tt284
hes(1,15) = tt292
hes(1,16) = tt297
hes(1,17) = tt301
hes(1,18) = tt305
hes(1,19) = tt310
hes(1,20) = tt314
hes(1,21) = tt318
hes(1,22) = tt321
hes(1,23) = tt323
hes(1,24) = tt325
hes(2,1) = tt156
hes(2,2) = tt1*(miu(1,1)*(tt134+tt133+tt132+tt131+tt129+tt118+tt3&
&31+tt330+tt329+tt328+tt327+tt326)+tt54+1.0E+0*lam(1,1)*tt138**2)
hes(2,3) = tt338
hes(2,4) = tt341
hes(2,5) = tt345
hes(2,6) = tt349
hes(2,7) = tt352
hes(2,8) = tt356
hes(2,9) = tt360
hes(2,10) = tt362
hes(2,11) = tt364
hes(2,12) = tt366
hes(2,13) = tt369
hes(2,14) = tt373
hes(2,15) = tt377
hes(2,16) = tt379
hes(2,17) = tt381
hes(2,18) = tt383
hes(2,19) = tt385
hes(2,20) = tt387
hes(2,21) = tt389
hes(2,22) = tt390
hes(2,23) = tt391
hes(2,24) = tt392
hes(3,1) = tt178
hes(3,2) = tt338
hes(3,3) = tt1*(1.0E+0*lam(1,1)*tt160**2+miu(1,1)*(tt134+tt398+tt&
&133+tt397+tt132+tt396+tt395+tt394+tt393+tt131+tt129+tt118)+tt54)
hes(3,4) = tt401
hes(3,5) = tt404
hes(3,6) = tt408
hes(3,7) = tt411
hes(3,8) = tt414
hes(3,9) = tt418
hes(3,10) = tt420
hes(3,11) = tt422
hes(3,12) = tt424
hes(3,13) = tt427
hes(3,14) = tt430
hes(3,15) = tt434
hes(3,16) = tt436
hes(3,17) = tt438
hes(3,18) = tt440
hes(3,19) = tt442
hes(3,20) = tt444
hes(3,21) = tt446
hes(3,22) = tt447
hes(3,23) = tt448
hes(3,24) = tt449
hes(4,1) = tt190
hes(4,2) = tt341
hes(4,3) = tt401
hes(4,4) = tt1*(miu(1,1)*(tt134+tt133+tt132+tt131+tt296+tt250+tt7&
&1+tt70+tt69+tt68+tt451+tt450)+tt54+1.0E+0*lam(1,1)*tt180**2)
hes(4,5) = tt454
hes(4,6) = tt457
hes(4,7) = tt459
hes(4,8) = tt461
hes(4,9) = tt463
hes(4,10) = tt465
hes(4,11) = tt467
hes(4,12) = tt469
hes(4,13) = tt471
hes(4,14) = tt473
hes(4,15) = tt475
hes(4,16) = tt477
hes(4,17) = tt479
hes(4,18) = tt481
hes(4,19) = tt482
hes(4,20) = tt483
hes(4,21) = tt484
hes(4,22) = tt485
hes(4,23) = tt486
hes(4,24) = tt487
hes(5,1) = tt200
hes(5,2) = tt345
hes(5,3) = tt404
hes(5,4) = tt454
hes(5,5) = tt1*(miu(1,1)*(tt134+tt133+tt132+tt131+tt296+tt250+tt3&
&31+tt330+tt329+tt328+tt489+tt488)+tt54+1.0E+0*lam(1,1)*tt192**2)
hes(5,6) = tt492
hes(5,7) = tt494
hes(5,8) = tt496
hes(5,9) = tt498
hes(5,10) = tt500
hes(5,11) = tt502
hes(5,12) = tt504
hes(5,13) = tt506
hes(5,14) = tt508
hes(5,15) = tt510
hes(5,16) = tt512
hes(5,17) = tt514
hes(5,18) = tt516
hes(5,19) = tt517
hes(5,20) = tt518
hes(5,21) = tt519
hes(5,22) = tt520
hes(5,23) = tt521
hes(5,24) = tt522
hes(6,1) = tt210
hes(6,2) = tt349
hes(6,3) = tt408
hes(6,4) = tt457
hes(6,5) = tt492
hes(6,6) = tt1*(1.0E+0*lam(1,1)*tt202**2+miu(1,1)*(tt134+tt398+tt&
&133+tt397+tt132+tt396+tt395+tt131+tt296+tt250+tt524+tt523)+tt54)
hes(6,7) = tt526
hes(6,8) = tt528
hes(6,9) = tt530
hes(6,10) = tt532
hes(6,11) = tt534
hes(6,12) = tt536
hes(6,13) = tt538
hes(6,14) = tt540
hes(6,15) = tt542
hes(6,16) = tt544
hes(6,17) = tt546
hes(6,18) = tt548
hes(6,19) = tt549
hes(6,20) = tt550
hes(6,21) = tt551
hes(6,22) = tt552
hes(6,23) = tt553
hes(6,24) = tt554
hes(7,1) = tt221
hes(7,2) = tt352
hes(7,3) = tt411
hes(7,4) = tt459
hes(7,5) = tt494
hes(7,6) = tt526
hes(7,7) = tt1*(miu(1,1)*(tt134+tt133+tt132+tt309+tt129+tt250+tt7&
&1+tt70+tt69+tt64+tt556+tt555)+tt54+1.0E+0*lam(1,1)*tt212**2)
hes(7,8) = tt559
hes(7,9) = tt562
hes(7,10) = tt564
hes(7,11) = tt566
hes(7,12) = tt568
hes(7,13) = tt570
hes(7,14) = tt572
hes(7,15) = tt574
hes(7,16) = tt575
hes(7,17) = tt576
hes(7,18) = tt577
hes(7,19) = tt579
hes(7,20) = tt581
hes(7,21) = tt583
hes(7,22) = tt584
hes(7,23) = tt585
hes(7,24) = tt586
hes(8,1) = tt231
hes(8,2) = tt356
hes(8,3) = tt414
hes(8,4) = tt461
hes(8,5) = tt496
hes(8,6) = tt528
hes(8,7) = tt559
hes(8,8) = tt1*(miu(1,1)*(tt134+tt133+tt132+tt309+tt129+tt250+tt3&
&31+tt330+tt329+tt327+tt588+tt587)+tt54+1.0E+0*lam(1,1)*tt223**2)
hes(8,9) = tt591
hes(8,10) = tt593
hes(8,11) = tt595
hes(8,12) = tt597
hes(8,13) = tt599
hes(8,14) = tt601
hes(8,15) = tt603
hes(8,16) = tt604
hes(8,17) = tt605
hes(8,18) = tt606
hes(8,19) = tt608
hes(8,20) = tt610
hes(8,21) = tt612
hes(8,22) = tt613
hes(8,23) = tt614
hes(8,24) = tt615
hes(9,1) = tt241
hes(9,2) = tt360
hes(9,3) = tt418
hes(9,4) = tt463
hes(9,5) = tt498
hes(9,6) = tt530
hes(9,7) = tt562
hes(9,8) = tt591
hes(9,9) = tt1*(1.0E+0*lam(1,1)*tt233**2+miu(1,1)*(tt134+tt398+tt&
&133+tt397+tt132+tt396+tt394+tt309+tt129+tt250+tt617+tt616)+tt54)
hes(9,10) = tt619
hes(9,11) = tt621
hes(9,12) = tt623
hes(9,13) = tt625
hes(9,14) = tt627
hes(9,15) = tt629
hes(9,16) = tt630
hes(9,17) = tt631
hes(9,18) = tt632
hes(9,19) = tt634
hes(9,20) = tt636
hes(9,21) = tt638
hes(9,22) = tt639
hes(9,23) = tt640
hes(9,24) = tt641
hes(10,1) = tt251
hes(10,2) = tt362
hes(10,3) = tt420
hes(10,4) = tt465
hes(10,5) = tt500
hes(10,6) = tt532
hes(10,7) = tt564
hes(10,8) = tt593
hes(10,9) = tt619
hes(10,10) = tt1*(miu(1,1)*(tt134+tt133+tt132+tt309+tt296+tt118+t&
&t71+tt70+tt69+tt642+tt450+tt555)+tt54+1.0E+0*lam(1,1)*tt242**2)
hes(10,11) = tt644
hes(10,12) = tt646
hes(10,13) = tt647
hes(10,14) = tt648
hes(10,15) = tt649
hes(10,16) = tt650
hes(10,17) = tt651
hes(10,18) = tt652
hes(10,19) = tt653
hes(10,20) = tt654
hes(10,21) = tt655
hes(10,22) = tt656
hes(10,23) = tt657
hes(10,24) = tt658
hes(11,1) = tt259
hes(11,2) = tt364
hes(11,3) = tt422
hes(11,4) = tt467
hes(11,5) = tt502
hes(11,6) = tt534
hes(11,7) = tt566
hes(11,8) = tt595
hes(11,9) = tt621
hes(11,10) = tt644
hes(11,11) = tt1*(miu(1,1)*(tt134+tt133+tt132+tt309+tt296+tt118+t&
&t331+tt330+tt329+tt659+tt488+tt587)+tt54+1.0E+0*lam(1,1)*tt252**2&
&)
hes(11,12) = tt661
hes(11,13) = tt662
hes(11,14) = tt663
hes(11,15) = tt664
hes(11,16) = tt665
hes(11,17) = tt666
hes(11,18) = tt667
hes(11,19) = tt668
hes(11,20) = tt669
hes(11,21) = tt670
hes(11,22) = tt671
hes(11,23) = tt672
hes(11,24) = tt673
hes(12,1) = tt267
hes(12,2) = tt366
hes(12,3) = tt424
hes(12,4) = tt469
hes(12,5) = tt504
hes(12,6) = tt536
hes(12,7) = tt568
hes(12,8) = tt597
hes(12,9) = tt623
hes(12,10) = tt646
hes(12,11) = tt661
hes(12,12) = tt1*(1.0E+0*lam(1,1)*tt260**2+miu(1,1)*(tt134+tt398+&
&tt133+tt397+tt132+tt396+tt674+tt309+tt296+tt118+tt523+tt616)+tt54&
&)
hes(12,13) = tt675
hes(12,14) = tt676
hes(12,15) = tt677
hes(12,16) = tt678
hes(12,17) = tt679
hes(12,18) = tt680
hes(12,19) = tt681
hes(12,20) = tt682
hes(12,21) = tt683
hes(12,22) = tt684
hes(12,23) = tt685
hes(12,24) = tt686
hes(13,1) = tt276
hes(13,2) = tt369
hes(13,3) = tt427
hes(13,4) = tt471
hes(13,5) = tt506
hes(13,6) = tt538
hes(13,7) = tt570
hes(13,8) = tt599
hes(13,9) = tt625
hes(13,10) = tt647
hes(13,11) = tt662
hes(13,12) = tt675
hes(13,13) = tt1*(miu(1,1)*(tt134+tt133+tt132+tt309+tt296+tt118+t&
&t71+tt70+tt69+tt60+tt688+tt687)+tt54+1.0E+0*lam(1,1)*tt269**2)
hes(13,14) = tt691
hes(13,15) = tt694
hes(13,16) = tt696
hes(13,17) = tt698
hes(13,18) = tt700
hes(13,19) = tt702
hes(13,20) = tt704
hes(13,21) = tt706
hes(13,22) = tt707
hes(13,23) = tt708
hes(13,24) = tt709
hes(14,1) = tt284
hes(14,2) = tt373
hes(14,3) = tt430
hes(14,4) = tt473
hes(14,5) = tt508
hes(14,6) = tt540
hes(14,7) = tt572
hes(14,8) = tt601
hes(14,9) = tt627
hes(14,10) = tt648
hes(14,11) = tt663
hes(14,12) = tt676
hes(14,13) = tt691
hes(14,14) = tt1*(miu(1,1)*(tt134+tt133+tt132+tt309+tt296+tt118+t&
&t331+tt330+tt329+tt326+tt711+tt710)+tt54+1.0E+0*lam(1,1)*tt278**2&
&)
hes(14,15) = tt714
hes(14,16) = tt716
hes(14,17) = tt718
hes(14,18) = tt720
hes(14,19) = tt722
hes(14,20) = tt724
hes(14,21) = tt726
hes(14,22) = tt727
hes(14,23) = tt728
hes(14,24) = tt729
hes(15,1) = tt292
hes(15,2) = tt377
hes(15,3) = tt434
hes(15,4) = tt475
hes(15,5) = tt510
hes(15,6) = tt542
hes(15,7) = tt574
hes(15,8) = tt603
hes(15,9) = tt629
hes(15,10) = tt649
hes(15,11) = tt664
hes(15,12) = tt677
hes(15,13) = tt694
hes(15,14) = tt714
hes(15,15) = tt1*(1.0E+0*lam(1,1)*tt286**2+miu(1,1)*(tt134+tt398+&
&tt133+tt397+tt132+tt396+tt393+tt309+tt296+tt118+tt731+tt730)+tt54&
&)
hes(15,16) = tt733
hes(15,17) = tt735
hes(15,18) = tt737
hes(15,19) = tt739
hes(15,20) = tt741
hes(15,21) = tt743
hes(15,22) = tt744
hes(15,23) = tt745
hes(15,24) = tt746
hes(16,1) = tt297
hes(16,2) = tt379
hes(16,3) = tt436
hes(16,4) = tt477
hes(16,5) = tt512
hes(16,6) = tt544
hes(16,7) = tt575
hes(16,8) = tt604
hes(16,9) = tt630
hes(16,10) = tt650
hes(16,11) = tt665
hes(16,12) = tt678
hes(16,13) = tt696
hes(16,14) = tt716
hes(16,15) = tt733
hes(16,16) = tt1*(miu(1,1)*(tt134+tt133+tt132+tt309+tt129+tt250+t&
&t71+tt70+tt69+tt747+tt451+tt687)+tt54+1.0E+0*lam(1,1)*tt293**2)
hes(16,17) = tt749
hes(16,18) = tt751
hes(16,19) = tt752
hes(16,20) = tt753
hes(16,21) = tt754
hes(16,22) = tt755
hes(16,23) = tt756
hes(16,24) = tt757
hes(17,1) = tt301
hes(17,2) = tt381
hes(17,3) = tt438
hes(17,4) = tt479
hes(17,5) = tt514
hes(17,6) = tt546
hes(17,7) = tt576
hes(17,8) = tt605
hes(17,9) = tt631
hes(17,10) = tt651
hes(17,11) = tt666
hes(17,12) = tt679
hes(17,13) = tt698
hes(17,14) = tt718
hes(17,15) = tt735
hes(17,16) = tt749
hes(17,17) = tt1*(miu(1,1)*(tt134+tt133+tt132+tt309+tt129+tt250+t&
&t331+tt330+tt329+tt758+tt489+tt710)+tt54+1.0E+0*lam(1,1)*tt298**2&
&)
hes(17,18) = tt760
hes(17,19) = tt761
hes(17,20) = tt762
hes(17,21) = tt763
hes(17,22) = tt764
hes(17,23) = tt765
hes(17,24) = tt766
hes(18,1) = tt305
hes(18,2) = tt383
hes(18,3) = tt440
hes(18,4) = tt481
hes(18,5) = tt516
hes(18,6) = tt548
hes(18,7) = tt577
hes(18,8) = tt606
hes(18,9) = tt632
hes(18,10) = tt652
hes(18,11) = tt667
hes(18,12) = tt680
hes(18,13) = tt700
hes(18,14) = tt720
hes(18,15) = tt737
hes(18,16) = tt751
hes(18,17) = tt760
hes(18,18) = tt1*(1.0E+0*lam(1,1)*tt302**2+miu(1,1)*(tt134+tt398+&
&tt133+tt397+tt132+tt396+tt767+tt309+tt129+tt250+tt524+tt730)+tt54&
&)
hes(18,19) = tt768
hes(18,20) = tt769
hes(18,21) = tt770
hes(18,22) = tt771
hes(18,23) = tt772
hes(18,24) = tt773
hes(19,1) = tt310
hes(19,2) = tt385
hes(19,3) = tt442
hes(19,4) = tt482
hes(19,5) = tt517
hes(19,6) = tt549
hes(19,7) = tt579
hes(19,8) = tt608
hes(19,9) = tt634
hes(19,10) = tt653
hes(19,11) = tt668
hes(19,12) = tt681
hes(19,13) = tt702
hes(19,14) = tt722
hes(19,15) = tt739
hes(19,16) = tt752
hes(19,17) = tt761
hes(19,18) = tt768
hes(19,19) = tt1*(miu(1,1)*(tt134+tt133+tt132+tt131+tt296+tt250+t&
&t71+tt70+tt69+tt774+tt556+tt688)+tt54+1.0E+0*lam(1,1)*tt306**2)
hes(19,20) = tt776
hes(19,21) = tt778
hes(19,22) = tt779
hes(19,23) = tt780
hes(19,24) = tt781
hes(20,1) = tt314
hes(20,2) = tt387
hes(20,3) = tt444
hes(20,4) = tt483
hes(20,5) = tt518
hes(20,6) = tt550
hes(20,7) = tt581
hes(20,8) = tt610
hes(20,9) = tt636
hes(20,10) = tt654
hes(20,11) = tt669
hes(20,12) = tt682
hes(20,13) = tt704
hes(20,14) = tt724
hes(20,15) = tt741
hes(20,16) = tt753
hes(20,17) = tt762
hes(20,18) = tt769
hes(20,19) = tt776
hes(20,20) = tt1*(miu(1,1)*(tt134+tt133+tt132+tt131+tt296+tt250+t&
&t331+tt330+tt329+tt782+tt588+tt711)+tt54+1.0E+0*lam(1,1)*tt311**2&
&)
hes(20,21) = tt784
hes(20,22) = tt785
hes(20,23) = tt786
hes(20,24) = tt787
hes(21,1) = tt318
hes(21,2) = tt389
hes(21,3) = tt446
hes(21,4) = tt484
hes(21,5) = tt519
hes(21,6) = tt551
hes(21,7) = tt583
hes(21,8) = tt612
hes(21,9) = tt638
hes(21,10) = tt655
hes(21,11) = tt670
hes(21,12) = tt683
hes(21,13) = tt706
hes(21,14) = tt726
hes(21,15) = tt743
hes(21,16) = tt754
hes(21,17) = tt763
hes(21,18) = tt770
hes(21,19) = tt778
hes(21,20) = tt784
hes(21,21) = tt1*(1.0E+0*lam(1,1)*tt315**2+miu(1,1)*(tt134+tt398+&
&tt133+tt397+tt132+tt396+tt788+tt131+tt296+tt250+tt617+tt731)+tt54&
&)
hes(21,22) = tt789
hes(21,23) = tt790
hes(21,24) = tt791
hes(22,1) = tt321
hes(22,2) = tt390
hes(22,3) = tt447
hes(22,4) = tt485
hes(22,5) = tt520
hes(22,6) = tt552
hes(22,7) = tt584
hes(22,8) = tt613
hes(22,9) = tt639
hes(22,10) = tt656
hes(22,11) = tt671
hes(22,12) = tt684
hes(22,13) = tt707
hes(22,14) = tt727
hes(22,15) = tt744
hes(22,16) = tt755
hes(22,17) = tt764
hes(22,18) = tt771
hes(22,19) = tt779
hes(22,20) = tt785
hes(22,21) = tt789
hes(22,22) = tt1*(miu(1,1)*(tt134+tt133+tt132+tt131+tt129+tt118+t&
&t71+tt70+tt69+tt774+tt747+tt642)+tt54+1.0E+0*lam(1,1)*tt319**2)
hes(22,23) = tt792
hes(22,24) = tt793
hes(23,1) = tt323
hes(23,2) = tt391
hes(23,3) = tt448
hes(23,4) = tt486
hes(23,5) = tt521
hes(23,6) = tt553
hes(23,7) = tt585
hes(23,8) = tt614
hes(23,9) = tt640
hes(23,10) = tt657
hes(23,11) = tt672
hes(23,12) = tt685
hes(23,13) = tt708
hes(23,14) = tt728
hes(23,15) = tt745
hes(23,16) = tt756
hes(23,17) = tt765
hes(23,18) = tt772
hes(23,19) = tt780
hes(23,20) = tt786
hes(23,21) = tt790
hes(23,22) = tt792
hes(23,23) = tt1*(miu(1,1)*(tt134+tt133+tt132+tt131+tt129+tt118+t&
&t331+tt330+tt329+tt782+tt758+tt659)+tt54+1.0E+0*lam(1,1)*tt322**2&
&)
hes(23,24) = tt794
hes(24,1) = tt325
hes(24,2) = tt392
hes(24,3) = tt449
hes(24,4) = tt487
hes(24,5) = tt522
hes(24,6) = tt554
hes(24,7) = tt586
hes(24,8) = tt615
hes(24,9) = tt641
hes(24,10) = tt658
hes(24,11) = tt673
hes(24,12) = tt686
hes(24,13) = tt709
hes(24,14) = tt729
hes(24,15) = tt746
hes(24,16) = tt757
hes(24,17) = tt766
hes(24,18) = tt773
hes(24,19) = tt781
hes(24,20) = tt787
hes(24,21) = tt791
hes(24,22) = tt793
hes(24,23) = tt794
hes(24,24) = tt1*(1.0E+0*lam(1,1)*tt324**2+miu(1,1)*(tt134+tt398+&
&tt133+tt397+tt132+tt396+tt788+tt767+tt674+tt131+tt129+tt118)+tt54&
&)
END 
