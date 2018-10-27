SUBROUTINE calc_stretch(val, X) 
IMPLICIT NONE 
REAL(KIND=8) val(1, 1) 
REAL(KIND=8) X(3, 2) 
val(1,1) = sqrt((X(3,1)-X(3,2))**2+(X(2,1)-X(2,2))**2+(X(1,1)-X(1&
&,2))**2)
END 
SUBROUTINE calc_stretch_jac(jac, X) 
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
SUBROUTINE calc_stretch_hes(hes, X) 
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
SUBROUTINE calc_bend(val, X) 
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
SUBROUTINE calc_bend_jac(jac, X) 
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
SUBROUTINE calc_bend_hes(hes, X) 
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
