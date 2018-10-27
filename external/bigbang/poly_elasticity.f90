SUBROUTINE polynomial_elas(val, X, Dm, coef) 
IMPLICIT NONE 
REAL(KIND=8) val(1, 1) 
REAL(KIND=8) X(2, 3) 
REAL(KIND=8) Dm(2, 2) 
REAL(KIND=8) coef(20, 1) 
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
tt1 = -X(1,1)
tt2 = X(1,2)+tt1
tt3 = X(1,3)+tt1
tt4 = tt3*Dm(2,1)+Dm(1,1)*tt2
tt5 = -X(2,1)
tt6 = X(2,2)+tt5
tt7 = X(2,3)+tt5
tt8 = Dm(2,1)*tt7+Dm(1,1)*tt6
tt9 = tt8**2+tt4**2
tt10 = tt3*Dm(2,2)+Dm(1,2)*tt2
tt11 = Dm(2,2)*tt7+Dm(1,2)*tt6
tt12 = tt11**2+tt10**2
tt13 = tt8*tt11+tt4*tt10
tt14 = tt9**2
tt15 = tt12**2
tt16 = tt13**2
val(1,1) = tt13*tt9*tt12*coef(20,1)+tt16*tt12*coef(19,1)+tt16*tt9&
&*coef(18,1)+tt13*tt15*coef(17,1)+tt9*tt15*coef(16,1)+tt13*tt14*co&
&ef(15,1)+tt14*tt12*coef(14,1)+tt13**3*coef(13,1)+tt12**3*coef(12,&
&1)+tt9**3*coef(11,1)+tt13*tt9*coef(10,1)+tt13*tt12*coef(9,1)+tt9*&
&tt12*coef(8,1)+tt16*coef(7,1)+tt15*coef(6,1)+tt14*coef(5,1)+tt13*&
&coef(4,1)+tt12*coef(3,1)+coef(2,1)*tt9+coef(1,1)
END 
SUBROUTINE polynomial_elas_jac(jac, X, Dm, coef) 
IMPLICIT NONE 
REAL(KIND=8) jac(1, 6) 
REAL(KIND=8) X(2, 3) 
REAL(KIND=8) Dm(2, 2) 
REAL(KIND=8) coef(20, 1) 
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
tt1 = (-Dm(2,1))-Dm(1,1)
tt2 = -X(1,1)
tt3 = X(1,2)+tt2
tt4 = X(1,3)+tt2
tt5 = tt4*Dm(2,1)+Dm(1,1)*tt3
tt6 = (-Dm(2,2))-Dm(1,2)
tt7 = tt4*Dm(2,2)+Dm(1,2)*tt3
tt8 = tt1*tt7+tt5*tt6
tt9 = -X(2,1)
tt10 = X(2,2)+tt9
tt11 = X(2,3)+tt9
tt12 = Dm(2,1)*tt11+Dm(1,1)*tt10
tt13 = tt12**2+tt5**2
tt14 = Dm(2,2)*tt11+Dm(1,2)*tt10
tt15 = tt14**2+tt7**2
tt16 = tt12*tt14+tt5*tt7
tt17 = tt13**2
tt18 = tt15**2
tt19 = tt16**2
tt20 = tt1*tt14+tt6*tt12
tt21 = Dm(1,1)*tt7+Dm(1,2)*tt5
tt22 = Dm(1,1)*tt14+Dm(1,2)*tt12
tt23 = Dm(2,1)*tt7+tt5*Dm(2,2)
tt24 = Dm(2,1)*tt14+Dm(2,2)*tt12
jac(1,1) = tt8*tt13*tt15*coef(20,1)+2*tt1*tt5*tt16*tt15*coef(20,1&
&)+2*tt6*tt7*tt16*tt13*coef(20,1)+2*tt6*tt7*tt19*coef(19,1)+2*tt8*&
&tt16*tt15*coef(19,1)+2*tt1*tt5*tt19*coef(18,1)+2*tt8*tt16*tt13*co&
&ef(18,1)+tt8*tt18*coef(17,1)+4*tt6*tt7*tt16*tt15*coef(17,1)+2*tt1&
&*tt5*tt18*coef(16,1)+4*tt6*tt7*tt13*tt15*coef(16,1)+tt8*tt17*coef&
&(15,1)+4*tt1*tt5*tt16*tt13*coef(15,1)+2*tt6*tt7*tt17*coef(14,1)+4&
&*tt1*tt5*tt13*tt15*coef(14,1)+3*tt8*tt19*coef(13,1)+6*tt6*tt7*tt1&
&8*coef(12,1)+6*tt1*tt5*tt17*coef(11,1)+tt8*tt13*coef(10,1)+2*tt1*&
&tt5*tt16*coef(10,1)+tt8*tt15*coef(9,1)+2*tt6*tt7*tt16*coef(9,1)+2&
&*tt1*tt5*tt15*coef(8,1)+2*tt6*tt7*tt13*coef(8,1)+2*tt8*tt16*coef(&
&7,1)+4*tt6*tt7*tt15*coef(6,1)+4*tt1*tt5*tt13*coef(5,1)+tt8*coef(4&
&,1)+2*tt6*tt7*coef(3,1)+2*tt1*tt5*coef(2,1)
jac(1,2) = tt20*tt13*tt15*coef(20,1)+2*tt1*tt12*tt16*tt15*coef(20&
&,1)+2*tt6*tt14*tt16*tt13*coef(20,1)+2*tt6*tt14*tt19*coef(19,1)+2*&
&tt20*tt16*tt15*coef(19,1)+2*tt1*tt12*tt19*coef(18,1)+2*tt20*tt16*&
&tt13*coef(18,1)+tt20*tt18*coef(17,1)+4*tt6*tt14*tt16*tt15*coef(17&
&,1)+2*tt1*tt12*tt18*coef(16,1)+4*tt6*tt14*tt13*tt15*coef(16,1)+tt&
&20*tt17*coef(15,1)+4*tt1*tt12*tt16*tt13*coef(15,1)+2*tt6*tt14*tt1&
&7*coef(14,1)+4*tt1*tt12*tt13*tt15*coef(14,1)+3*tt20*tt19*coef(13,&
&1)+6*tt6*tt14*tt18*coef(12,1)+6*tt1*tt12*tt17*coef(11,1)+tt20*tt1&
&3*coef(10,1)+2*tt1*tt12*tt16*coef(10,1)+tt20*tt15*coef(9,1)+2*tt6&
&*tt14*tt16*coef(9,1)+2*tt1*tt12*tt15*coef(8,1)+2*tt6*tt14*tt13*co&
&ef(8,1)+2*tt20*tt16*coef(7,1)+4*tt6*tt14*tt15*coef(6,1)+4*tt1*tt1&
&2*tt13*coef(5,1)+tt20*coef(4,1)+2*tt6*tt14*coef(3,1)+2*tt1*coef(2&
&,1)*tt12
jac(1,3) = tt21*tt13*tt15*coef(20,1)+2*Dm(1,1)*tt5*tt16*tt15*coef&
&(20,1)+2*Dm(1,2)*tt7*tt16*tt13*coef(20,1)+2*Dm(1,2)*tt7*tt19*coef&
&(19,1)+2*tt21*tt16*tt15*coef(19,1)+2*Dm(1,1)*tt5*tt19*coef(18,1)+&
&2*tt21*tt16*tt13*coef(18,1)+tt21*tt18*coef(17,1)+4*Dm(1,2)*tt7*tt&
&16*tt15*coef(17,1)+2*Dm(1,1)*tt5*tt18*coef(16,1)+4*Dm(1,2)*tt7*tt&
&13*tt15*coef(16,1)+tt21*tt17*coef(15,1)+4*Dm(1,1)*tt5*tt16*tt13*c&
&oef(15,1)+2*Dm(1,2)*tt7*tt17*coef(14,1)+4*Dm(1,1)*tt5*tt13*tt15*c&
&oef(14,1)+3*tt21*tt19*coef(13,1)+6*Dm(1,2)*tt7*tt18*coef(12,1)+6*&
&Dm(1,1)*tt5*tt17*coef(11,1)+tt21*tt13*coef(10,1)+2*Dm(1,1)*tt5*tt&
&16*coef(10,1)+tt21*tt15*coef(9,1)+2*Dm(1,2)*tt7*tt16*coef(9,1)+2*&
&Dm(1,1)*tt5*tt15*coef(8,1)+2*Dm(1,2)*tt7*tt13*coef(8,1)+2*tt21*tt&
&16*coef(7,1)+4*Dm(1,2)*tt7*tt15*coef(6,1)+4*Dm(1,1)*tt5*tt13*coef&
&(5,1)+tt21*coef(4,1)+2*Dm(1,2)*tt7*coef(3,1)+2*Dm(1,1)*tt5*coef(2&
&,1)
jac(1,4) = tt22*tt13*tt15*coef(20,1)+2*Dm(1,1)*tt12*tt16*tt15*coe&
&f(20,1)+2*Dm(1,2)*tt14*tt16*tt13*coef(20,1)+2*Dm(1,2)*tt14*tt19*c&
&oef(19,1)+2*tt22*tt16*tt15*coef(19,1)+2*Dm(1,1)*tt12*tt19*coef(18&
&,1)+2*tt22*tt16*tt13*coef(18,1)+tt22*tt18*coef(17,1)+4*Dm(1,2)*tt&
&14*tt16*tt15*coef(17,1)+2*Dm(1,1)*tt12*tt18*coef(16,1)+4*Dm(1,2)*&
&tt14*tt13*tt15*coef(16,1)+tt22*tt17*coef(15,1)+4*Dm(1,1)*tt12*tt1&
&6*tt13*coef(15,1)+2*Dm(1,2)*tt14*tt17*coef(14,1)+4*Dm(1,1)*tt12*t&
&t13*tt15*coef(14,1)+3*tt22*tt19*coef(13,1)+6*Dm(1,2)*tt14*tt18*co&
&ef(12,1)+6*Dm(1,1)*tt12*tt17*coef(11,1)+tt22*tt13*coef(10,1)+2*Dm&
&(1,1)*tt12*tt16*coef(10,1)+tt22*tt15*coef(9,1)+2*Dm(1,2)*tt14*tt1&
&6*coef(9,1)+2*Dm(1,1)*tt12*tt15*coef(8,1)+2*Dm(1,2)*tt14*tt13*coe&
&f(8,1)+2*tt22*tt16*coef(7,1)+4*Dm(1,2)*tt14*tt15*coef(6,1)+4*Dm(1&
&,1)*tt12*tt13*coef(5,1)+tt22*coef(4,1)+2*Dm(1,2)*tt14*coef(3,1)+2&
&*Dm(1,1)*coef(2,1)*tt12
jac(1,5) = tt23*tt13*tt15*coef(20,1)+2*Dm(2,1)*tt5*tt16*tt15*coef&
&(20,1)+2*Dm(2,2)*tt7*tt16*tt13*coef(20,1)+2*Dm(2,2)*tt7*tt19*coef&
&(19,1)+2*tt23*tt16*tt15*coef(19,1)+2*Dm(2,1)*tt5*tt19*coef(18,1)+&
&2*tt23*tt16*tt13*coef(18,1)+tt23*tt18*coef(17,1)+4*Dm(2,2)*tt7*tt&
&16*tt15*coef(17,1)+2*Dm(2,1)*tt5*tt18*coef(16,1)+4*Dm(2,2)*tt7*tt&
&13*tt15*coef(16,1)+tt23*tt17*coef(15,1)+4*Dm(2,1)*tt5*tt16*tt13*c&
&oef(15,1)+2*Dm(2,2)*tt7*tt17*coef(14,1)+4*Dm(2,1)*tt5*tt13*tt15*c&
&oef(14,1)+3*tt23*tt19*coef(13,1)+6*Dm(2,2)*tt7*tt18*coef(12,1)+6*&
&Dm(2,1)*tt5*tt17*coef(11,1)+tt23*tt13*coef(10,1)+2*Dm(2,1)*tt5*tt&
&16*coef(10,1)+tt23*tt15*coef(9,1)+2*Dm(2,2)*tt7*tt16*coef(9,1)+2*&
&Dm(2,1)*tt5*tt15*coef(8,1)+2*Dm(2,2)*tt7*tt13*coef(8,1)+2*tt23*tt&
&16*coef(7,1)+4*Dm(2,2)*tt7*tt15*coef(6,1)+4*Dm(2,1)*tt5*tt13*coef&
&(5,1)+tt23*coef(4,1)+2*Dm(2,2)*tt7*coef(3,1)+2*Dm(2,1)*tt5*coef(2&
&,1)
jac(1,6) = tt24*tt13*tt15*coef(20,1)+2*Dm(2,1)*tt12*tt16*tt15*coe&
&f(20,1)+2*Dm(2,2)*tt14*tt16*tt13*coef(20,1)+2*Dm(2,2)*tt14*tt19*c&
&oef(19,1)+2*tt24*tt16*tt15*coef(19,1)+2*Dm(2,1)*tt12*tt19*coef(18&
&,1)+2*tt24*tt16*tt13*coef(18,1)+tt24*tt18*coef(17,1)+4*Dm(2,2)*tt&
&14*tt16*tt15*coef(17,1)+2*Dm(2,1)*tt12*tt18*coef(16,1)+4*Dm(2,2)*&
&tt14*tt13*tt15*coef(16,1)+tt24*tt17*coef(15,1)+4*Dm(2,1)*tt12*tt1&
&6*tt13*coef(15,1)+2*Dm(2,2)*tt14*tt17*coef(14,1)+4*Dm(2,1)*tt12*t&
&t13*tt15*coef(14,1)+3*tt24*tt19*coef(13,1)+6*Dm(2,2)*tt14*tt18*co&
&ef(12,1)+6*Dm(2,1)*tt12*tt17*coef(11,1)+tt24*tt13*coef(10,1)+2*Dm&
&(2,1)*tt12*tt16*coef(10,1)+tt24*tt15*coef(9,1)+2*Dm(2,2)*tt14*tt1&
&6*coef(9,1)+2*Dm(2,1)*tt12*tt15*coef(8,1)+2*Dm(2,2)*tt14*tt13*coe&
&f(8,1)+2*tt24*tt16*coef(7,1)+4*Dm(2,2)*tt14*tt15*coef(6,1)+4*Dm(2&
&,1)*tt12*tt13*coef(5,1)+tt24*coef(4,1)+2*Dm(2,2)*tt14*coef(3,1)+2&
&*Dm(2,1)*coef(2,1)*tt12
END 
SUBROUTINE polynomial_elas_hes(hes, X, Dm, coef) 
IMPLICIT NONE 
REAL(KIND=8) hes(6, 6) 
REAL(KIND=8) X(2, 3) 
REAL(KIND=8) Dm(2, 2) 
REAL(KIND=8) coef(20, 1) 
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
tt1 = (-Dm(2,1))-Dm(1,1)
tt2 = tt1**2
tt3 = 2*tt2*coef(2,1)
tt4 = (-Dm(2,2))-Dm(1,2)
tt5 = tt4**2
tt6 = 2*tt5*coef(3,1)
tt7 = 2*tt1*tt4*coef(4,1)
tt8 = -X(1,1)
tt9 = X(1,2)+tt8
tt10 = X(1,3)+tt8
tt11 = tt10*Dm(2,1)+Dm(1,1)*tt9
tt12 = tt11**2
tt13 = -X(2,1)
tt14 = X(2,2)+tt13
tt15 = X(2,3)+tt13
tt16 = Dm(2,1)*tt15+Dm(1,1)*tt14
tt17 = tt16**2
tt18 = tt17+tt12
tt19 = 4*tt2*tt18*coef(5,1)
tt20 = tt10*Dm(2,2)+Dm(1,2)*tt9
tt21 = tt20**2
tt22 = Dm(2,2)*tt15+Dm(1,2)*tt14
tt23 = tt22**2
tt24 = tt23+tt21
tt25 = 4*tt5*tt24*coef(6,1)
tt26 = tt1*tt20+tt11*tt4
tt27 = tt26**2
tt28 = tt16*tt22+tt11*tt20
tt29 = 4*tt1*tt4*tt28*coef(7,1)
tt30 = 2*tt5*tt18*coef(8,1)
tt31 = 2*tt2*tt24*coef(8,1)
tt32 = 2*tt5*tt28*coef(9,1)
tt33 = 2*tt1*tt4*tt24*coef(9,1)
tt34 = 2*tt2*tt28*coef(10,1)
tt35 = 2*tt1*tt4*tt18*coef(10,1)
tt36 = tt18**2
tt37 = 6*tt2*tt36*coef(11,1)
tt38 = tt24**2
tt39 = 6*tt5*tt38*coef(12,1)
tt40 = tt28**2
tt41 = 6*tt1*tt4*tt40*coef(13,1)
tt42 = 4*tt2*tt18*tt24*coef(14,1)
tt43 = 2*tt5*tt36*coef(14,1)
tt44 = 4*tt2*tt28*tt18*coef(15,1)
tt45 = 2*tt1*tt4*tt36*coef(15,1)
tt46 = 4*tt5*tt18*tt24*coef(16,1)
tt47 = 2*tt2*tt38*coef(16,1)
tt48 = 4*tt5*tt28*tt24*coef(17,1)
tt49 = 2*tt1*tt4*tt38*coef(17,1)
tt50 = 4*tt1*tt4*tt28*tt18*coef(18,1)
tt51 = 2*tt2*tt40*coef(18,1)
tt52 = 4*tt1*tt4*tt28*tt24*coef(19,1)
tt53 = 2*tt5*tt40*coef(19,1)
tt54 = 2*tt5*tt28*tt18*coef(20,1)
tt55 = 2*tt2*tt28*tt24*coef(20,1)
tt56 = 2*tt1*tt4*tt18*tt24*coef(20,1)
tt57 = tt1*tt22+tt4*tt16
tt58 = 2*tt1*tt11*tt57*tt24*coef(20,1)+2*tt1*tt26*tt16*tt24*coef(&
&20,1)+2*tt4*tt20*tt57*tt18*coef(20,1)+2*tt4*tt26*tt22*tt18*coef(2&
&0,1)+4*tt1*tt11*tt4*tt22*tt28*coef(20,1)+4*tt1*tt4*tt20*tt16*tt28&
&*coef(20,1)+2*tt26*tt57*tt24*coef(19,1)+4*tt4*tt20*tt57*tt28*coef&
&(19,1)+4*tt4*tt26*tt22*tt28*coef(19,1)+2*tt26*tt57*tt18*coef(18,1&
&)+4*tt1*tt11*tt57*tt28*coef(18,1)+4*tt1*tt26*tt16*tt28*coef(18,1)&
&+4*tt4*tt20*tt57*tt24*coef(17,1)+4*tt4*tt26*tt22*tt24*coef(17,1)+&
&8*tt5*tt20*tt22*tt28*coef(17,1)+8*tt1*tt11*tt4*tt22*tt24*coef(16,&
&1)+8*tt1*tt4*tt20*tt16*tt24*coef(16,1)+8*tt5*tt20*tt22*tt18*coef(&
&16,1)+4*tt1*tt11*tt57*tt18*coef(15,1)+4*tt1*tt26*tt16*tt18*coef(1&
&5,1)+8*tt2*tt11*tt16*tt28*coef(15,1)+8*tt2*tt11*tt16*tt24*coef(14&
&,1)+8*tt1*tt11*tt4*tt22*tt18*coef(14,1)+8*tt1*tt4*tt20*tt16*tt18*&
&coef(14,1)+6*tt26*tt57*tt28*coef(13,1)+24*tt5*tt20*tt22*tt24*coef&
&(12,1)+24*tt2*tt11*tt16*tt18*coef(11,1)+2*tt1*tt11*tt57*coef(10,1&
&)+2*tt1*tt26*tt16*coef(10,1)+2*tt4*tt20*tt57*coef(9,1)+2*tt4*tt26&
&*tt22*coef(9,1)+4*tt1*tt11*tt4*tt22*coef(8,1)+4*tt1*tt4*tt20*tt16&
&*coef(8,1)+2*tt26*tt57*coef(7,1)+8*tt5*tt20*tt22*coef(6,1)+8*tt2*&
&tt11*tt16*coef(5,1)
tt59 = 2*Dm(1,1)*tt1*coef(2,1)
tt60 = 2*Dm(1,2)*tt4*coef(3,1)
tt61 = Dm(1,1)*tt4+Dm(1,2)*tt1
tt62 = tt61*coef(4,1)
tt63 = 4*Dm(1,1)*tt1*tt18*coef(5,1)
tt64 = 4*Dm(1,2)*tt4*tt24*coef(6,1)
tt65 = Dm(1,1)*tt20+Dm(1,2)*tt11
tt66 = 2*tt61*tt28*coef(7,1)
tt67 = 2*Dm(1,2)*tt4*tt18*coef(8,1)
tt68 = 2*Dm(1,1)*tt1*tt24*coef(8,1)
tt69 = 2*Dm(1,2)*tt4*tt28*coef(9,1)
tt70 = tt61*tt24*coef(9,1)
tt71 = 2*Dm(1,1)*tt1*tt28*coef(10,1)
tt72 = tt61*tt18*coef(10,1)
tt73 = 6*Dm(1,1)*tt1*tt36*coef(11,1)
tt74 = 6*Dm(1,2)*tt4*tt38*coef(12,1)
tt75 = 3*tt61*tt40*coef(13,1)
tt76 = 4*Dm(1,1)*tt1*tt18*tt24*coef(14,1)
tt77 = 2*Dm(1,2)*tt4*tt36*coef(14,1)
tt78 = 4*Dm(1,1)*tt1*tt28*tt18*coef(15,1)
tt79 = tt61*tt36*coef(15,1)
tt80 = 4*Dm(1,2)*tt4*tt18*tt24*coef(16,1)
tt81 = 2*Dm(1,1)*tt1*tt38*coef(16,1)
tt82 = 4*Dm(1,2)*tt4*tt28*tt24*coef(17,1)
tt83 = tt61*tt38*coef(17,1)
tt84 = 2*tt61*tt28*tt18*coef(18,1)
tt85 = 2*Dm(1,1)*tt1*tt40*coef(18,1)
tt86 = 2*tt61*tt28*tt24*coef(19,1)
tt87 = 2*Dm(1,2)*tt4*tt40*coef(19,1)
tt88 = 2*Dm(1,2)*tt4*tt28*tt18*coef(20,1)
tt89 = 2*Dm(1,1)*tt1*tt28*tt24*coef(20,1)
tt90 = tt61*tt18*tt24*coef(20,1)
tt91 = tt90+tt89+2*Dm(1,1)*tt11*tt26*tt24*coef(20,1)+2*tt1*tt11*t&
&t65*tt24*coef(20,1)+tt88+2*Dm(1,2)*tt20*tt26*tt18*coef(20,1)+2*tt&
&4*tt20*tt65*tt18*coef(20,1)+4*Dm(1,1)*tt11*tt4*tt20*tt28*coef(20,&
&1)+4*Dm(1,2)*tt1*tt11*tt20*tt28*coef(20,1)+tt87+tt86+2*tt65*tt26*&
&tt24*coef(19,1)+4*Dm(1,2)*tt20*tt26*tt28*coef(19,1)+4*tt4*tt20*tt&
&65*tt28*coef(19,1)+tt85+tt84+2*tt65*tt26*tt18*coef(18,1)+4*Dm(1,1&
&)*tt11*tt26*tt28*coef(18,1)+4*tt1*tt11*tt65*tt28*coef(18,1)+tt83+&
&tt82+4*Dm(1,2)*tt20*tt26*tt24*coef(17,1)+4*tt4*tt20*tt65*tt24*coe&
&f(17,1)+8*Dm(1,2)*tt4*tt21*tt28*coef(17,1)+tt81+tt80+8*Dm(1,1)*tt&
&11*tt4*tt20*tt24*coef(16,1)+8*Dm(1,2)*tt1*tt11*tt20*tt24*coef(16,&
&1)+8*Dm(1,2)*tt4*tt21*tt18*coef(16,1)+tt79+tt78+4*Dm(1,1)*tt11*tt&
&26*tt18*coef(15,1)+4*tt1*tt11*tt65*tt18*coef(15,1)+8*Dm(1,1)*tt1*&
&tt12*tt28*coef(15,1)+tt77+tt76+8*Dm(1,1)*tt1*tt12*tt24*coef(14,1)&
&+8*Dm(1,1)*tt11*tt4*tt20*tt18*coef(14,1)+8*Dm(1,2)*tt1*tt11*tt20*&
&tt18*coef(14,1)+tt75+6*tt65*tt26*tt28*coef(13,1)+tt74+24*Dm(1,2)*&
&tt4*tt21*tt24*coef(12,1)+tt73+24*Dm(1,1)*tt1*tt12*tt18*coef(11,1)&
&+tt72+tt71+2*Dm(1,1)*tt11*tt26*coef(10,1)+2*tt1*tt11*tt65*coef(10&
&,1)+tt70+tt69+2*Dm(1,2)*tt20*tt26*coef(9,1)+2*tt4*tt20*tt65*coef(&
&9,1)+tt68+tt67+4*Dm(1,1)*tt11*tt4*tt20*coef(8,1)+4*Dm(1,2)*tt1*tt&
&11*tt20*coef(8,1)+tt66+2*tt65*tt26*coef(7,1)+tt64+8*Dm(1,2)*tt4*t&
&t21*coef(6,1)+tt63+8*Dm(1,1)*tt1*tt12*coef(5,1)+tt62+tt60+tt59
tt92 = 8*Dm(1,1)*tt1*tt11*tt16*coef(5,1)
tt93 = 8*Dm(1,2)*tt4*tt20*tt22*coef(6,1)
tt94 = Dm(1,1)*tt22+Dm(1,2)*tt16
tt95 = 24*Dm(1,1)*tt1*tt11*tt16*tt18*coef(11,1)
tt96 = 24*Dm(1,2)*tt4*tt20*tt22*tt24*coef(12,1)
tt97 = 8*Dm(1,1)*tt1*tt11*tt16*tt24*coef(14,1)
tt98 = 8*Dm(1,1)*tt1*tt11*tt16*tt28*coef(15,1)
tt99 = 8*Dm(1,2)*tt4*tt20*tt22*tt18*coef(16,1)
tt100 = 8*Dm(1,2)*tt4*tt20*tt22*tt28*coef(17,1)
tt101 = 2*tt1*tt11*tt94*tt24*coef(20,1)+2*Dm(1,1)*tt26*tt16*tt24*&
&coef(20,1)+2*tt4*tt20*tt94*tt18*coef(20,1)+2*Dm(1,2)*tt26*tt22*tt&
&18*coef(20,1)+4*Dm(1,2)*tt1*tt11*tt22*tt28*coef(20,1)+4*Dm(1,1)*t&
&t4*tt20*tt16*tt28*coef(20,1)+2*tt26*tt94*tt24*coef(19,1)+4*tt4*tt&
&20*tt94*tt28*coef(19,1)+4*Dm(1,2)*tt26*tt22*tt28*coef(19,1)+2*tt2&
&6*tt94*tt18*coef(18,1)+4*tt1*tt11*tt94*tt28*coef(18,1)+4*Dm(1,1)*&
&tt26*tt16*tt28*coef(18,1)+4*tt4*tt20*tt94*tt24*coef(17,1)+4*Dm(1,&
&2)*tt26*tt22*tt24*coef(17,1)+tt100+8*Dm(1,2)*tt1*tt11*tt22*tt24*c&
&oef(16,1)+8*Dm(1,1)*tt4*tt20*tt16*tt24*coef(16,1)+tt99+4*tt1*tt11&
&*tt94*tt18*coef(15,1)+4*Dm(1,1)*tt26*tt16*tt18*coef(15,1)+tt98+tt&
&97+8*Dm(1,2)*tt1*tt11*tt22*tt18*coef(14,1)+8*Dm(1,1)*tt4*tt20*tt1&
&6*tt18*coef(14,1)+6*tt26*tt94*tt28*coef(13,1)+tt96+tt95+2*tt1*tt1&
&1*tt94*coef(10,1)+2*Dm(1,1)*tt26*tt16*coef(10,1)+2*tt4*tt20*tt94*&
&coef(9,1)+2*Dm(1,2)*tt26*tt22*coef(9,1)+4*Dm(1,2)*tt1*tt11*tt22*c&
&oef(8,1)+4*Dm(1,1)*tt4*tt20*tt16*coef(8,1)+2*tt26*tt94*coef(7,1)+&
&tt93+tt92
tt102 = 2*tt1*Dm(2,1)*coef(2,1)
tt103 = 2*tt4*Dm(2,2)*coef(3,1)
tt104 = tt1*Dm(2,2)+Dm(2,1)*tt4
tt105 = tt104*coef(4,1)
tt106 = 4*tt1*Dm(2,1)*tt18*coef(5,1)
tt107 = 4*tt4*Dm(2,2)*tt24*coef(6,1)
tt108 = Dm(2,1)*tt20+tt11*Dm(2,2)
tt109 = 2*tt104*tt28*coef(7,1)
tt110 = 2*tt4*Dm(2,2)*tt18*coef(8,1)
tt111 = 2*tt1*Dm(2,1)*tt24*coef(8,1)
tt112 = 2*tt4*Dm(2,2)*tt28*coef(9,1)
tt113 = tt104*tt24*coef(9,1)
tt114 = 2*tt1*Dm(2,1)*tt28*coef(10,1)
tt115 = tt104*tt18*coef(10,1)
tt116 = 6*tt1*Dm(2,1)*tt36*coef(11,1)
tt117 = 6*tt4*Dm(2,2)*tt38*coef(12,1)
tt118 = 3*tt104*tt40*coef(13,1)
tt119 = 4*tt1*Dm(2,1)*tt18*tt24*coef(14,1)
tt120 = 2*tt4*Dm(2,2)*tt36*coef(14,1)
tt121 = 4*tt1*Dm(2,1)*tt28*tt18*coef(15,1)
tt122 = tt104*tt36*coef(15,1)
tt123 = 4*tt4*Dm(2,2)*tt18*tt24*coef(16,1)
tt124 = 2*tt1*Dm(2,1)*tt38*coef(16,1)
tt125 = 4*tt4*Dm(2,2)*tt28*tt24*coef(17,1)
tt126 = tt104*tt38*coef(17,1)
tt127 = 2*tt104*tt28*tt18*coef(18,1)
tt128 = 2*tt1*Dm(2,1)*tt40*coef(18,1)
tt129 = 2*tt104*tt28*tt24*coef(19,1)
tt130 = 2*tt4*Dm(2,2)*tt40*coef(19,1)
tt131 = 2*tt4*Dm(2,2)*tt28*tt18*coef(20,1)
tt132 = 2*tt1*Dm(2,1)*tt28*tt24*coef(20,1)
tt133 = tt104*tt18*tt24*coef(20,1)
tt134 = tt133+tt132+2*tt1*tt11*tt108*tt24*coef(20,1)+2*Dm(2,1)*tt&
&11*tt26*tt24*coef(20,1)+tt131+2*tt4*tt20*tt108*tt18*coef(20,1)+2*&
&Dm(2,2)*tt20*tt26*tt18*coef(20,1)+4*tt1*tt11*Dm(2,2)*tt20*tt28*co&
&ef(20,1)+4*Dm(2,1)*tt11*tt4*tt20*tt28*coef(20,1)+tt130+tt129+2*tt&
&26*tt108*tt24*coef(19,1)+4*tt4*tt20*tt108*tt28*coef(19,1)+4*Dm(2,&
&2)*tt20*tt26*tt28*coef(19,1)+tt128+tt127+2*tt26*tt108*tt18*coef(1&
&8,1)+4*tt1*tt11*tt108*tt28*coef(18,1)+4*Dm(2,1)*tt11*tt26*tt28*co&
&ef(18,1)+tt126+tt125+4*tt4*tt20*tt108*tt24*coef(17,1)+4*Dm(2,2)*t&
&t20*tt26*tt24*coef(17,1)+8*tt4*Dm(2,2)*tt21*tt28*coef(17,1)+tt124&
&+tt123+8*tt1*tt11*Dm(2,2)*tt20*tt24*coef(16,1)+8*Dm(2,1)*tt11*tt4&
&*tt20*tt24*coef(16,1)+8*tt4*Dm(2,2)*tt21*tt18*coef(16,1)+tt122+tt&
&121+4*tt1*tt11*tt108*tt18*coef(15,1)+4*Dm(2,1)*tt11*tt26*tt18*coe&
&f(15,1)+8*tt1*Dm(2,1)*tt12*tt28*coef(15,1)+tt120+tt119+8*tt1*Dm(2&
&,1)*tt12*tt24*coef(14,1)+8*tt1*tt11*Dm(2,2)*tt20*tt18*coef(14,1)+&
&8*Dm(2,1)*tt11*tt4*tt20*tt18*coef(14,1)+tt118+6*tt26*tt108*tt28*c&
&oef(13,1)+tt117+24*tt4*Dm(2,2)*tt21*tt24*coef(12,1)+tt116+24*tt1*&
&Dm(2,1)*tt12*tt18*coef(11,1)+tt115+tt114+2*tt1*tt11*tt108*coef(10&
&,1)+2*Dm(2,1)*tt11*tt26*coef(10,1)+tt113+tt112+2*tt4*tt20*tt108*c&
&oef(9,1)+2*Dm(2,2)*tt20*tt26*coef(9,1)+tt111+tt110+4*tt1*tt11*Dm(&
&2,2)*tt20*coef(8,1)+4*Dm(2,1)*tt11*tt4*tt20*coef(8,1)+tt109+2*tt2&
&6*tt108*coef(7,1)+tt107+8*tt4*Dm(2,2)*tt21*coef(6,1)+tt106+8*tt1*&
&Dm(2,1)*tt12*coef(5,1)+tt105+tt103+tt102
tt135 = 8*tt1*Dm(2,1)*tt11*tt16*coef(5,1)
tt136 = 8*tt4*Dm(2,2)*tt20*tt22*coef(6,1)
tt137 = Dm(2,1)*tt22+Dm(2,2)*tt16
tt138 = 24*tt1*Dm(2,1)*tt11*tt16*tt18*coef(11,1)
tt139 = 24*tt4*Dm(2,2)*tt20*tt22*tt24*coef(12,1)
tt140 = 8*tt1*Dm(2,1)*tt11*tt16*tt24*coef(14,1)
tt141 = 8*tt1*Dm(2,1)*tt11*tt16*tt28*coef(15,1)
tt142 = 8*tt4*Dm(2,2)*tt20*tt22*tt18*coef(16,1)
tt143 = 8*tt4*Dm(2,2)*tt20*tt22*tt28*coef(17,1)
tt144 = 2*tt1*tt11*tt137*tt24*coef(20,1)+2*Dm(2,1)*tt26*tt16*tt24&
&*coef(20,1)+2*tt4*tt20*tt137*tt18*coef(20,1)+2*Dm(2,2)*tt26*tt22*&
&tt18*coef(20,1)+4*tt1*tt11*Dm(2,2)*tt22*tt28*coef(20,1)+4*Dm(2,1)&
&*tt4*tt20*tt16*tt28*coef(20,1)+2*tt26*tt137*tt24*coef(19,1)+4*tt4&
&*tt20*tt137*tt28*coef(19,1)+4*Dm(2,2)*tt26*tt22*tt28*coef(19,1)+2&
&*tt26*tt137*tt18*coef(18,1)+4*tt1*tt11*tt137*tt28*coef(18,1)+4*Dm&
&(2,1)*tt26*tt16*tt28*coef(18,1)+4*tt4*tt20*tt137*tt24*coef(17,1)+&
&4*Dm(2,2)*tt26*tt22*tt24*coef(17,1)+tt143+8*tt1*tt11*Dm(2,2)*tt22&
&*tt24*coef(16,1)+8*Dm(2,1)*tt4*tt20*tt16*tt24*coef(16,1)+tt142+4*&
&tt1*tt11*tt137*tt18*coef(15,1)+4*Dm(2,1)*tt26*tt16*tt18*coef(15,1&
&)+tt141+tt140+8*tt1*tt11*Dm(2,2)*tt22*tt18*coef(14,1)+8*Dm(2,1)*t&
&t4*tt20*tt16*tt18*coef(14,1)+6*tt26*tt137*tt28*coef(13,1)+tt139+t&
&t138+2*tt1*tt11*tt137*coef(10,1)+2*Dm(2,1)*tt26*tt16*coef(10,1)+2&
&*tt4*tt20*tt137*coef(9,1)+2*Dm(2,2)*tt26*tt22*coef(9,1)+4*tt1*tt1&
&1*Dm(2,2)*tt22*coef(8,1)+4*Dm(2,1)*tt4*tt20*tt16*coef(8,1)+2*tt26&
&*tt137*coef(7,1)+tt136+tt135
tt145 = tt57**2
tt146 = 2*Dm(1,1)*tt11*tt57*tt24*coef(20,1)+2*tt1*tt65*tt16*tt24*&
&coef(20,1)+2*Dm(1,2)*tt20*tt57*tt18*coef(20,1)+2*tt4*tt65*tt22*tt&
&18*coef(20,1)+4*Dm(1,1)*tt11*tt4*tt22*tt28*coef(20,1)+4*Dm(1,2)*t&
&t1*tt20*tt16*tt28*coef(20,1)+2*tt65*tt57*tt24*coef(19,1)+4*Dm(1,2&
&)*tt20*tt57*tt28*coef(19,1)+4*tt4*tt65*tt22*tt28*coef(19,1)+2*tt6&
&5*tt57*tt18*coef(18,1)+4*Dm(1,1)*tt11*tt57*tt28*coef(18,1)+4*tt1*&
&tt65*tt16*tt28*coef(18,1)+4*Dm(1,2)*tt20*tt57*tt24*coef(17,1)+4*t&
&t4*tt65*tt22*tt24*coef(17,1)+tt100+8*Dm(1,1)*tt11*tt4*tt22*tt24*c&
&oef(16,1)+8*Dm(1,2)*tt1*tt20*tt16*tt24*coef(16,1)+tt99+4*Dm(1,1)*&
&tt11*tt57*tt18*coef(15,1)+4*tt1*tt65*tt16*tt18*coef(15,1)+tt98+tt&
&97+8*Dm(1,1)*tt11*tt4*tt22*tt18*coef(14,1)+8*Dm(1,2)*tt1*tt20*tt1&
&6*tt18*coef(14,1)+6*tt65*tt57*tt28*coef(13,1)+tt96+tt95+2*Dm(1,1)&
&*tt11*tt57*coef(10,1)+2*tt1*tt65*tt16*coef(10,1)+2*Dm(1,2)*tt20*t&
&t57*coef(9,1)+2*tt4*tt65*tt22*coef(9,1)+4*Dm(1,1)*tt11*tt4*tt22*c&
&oef(8,1)+4*Dm(1,2)*tt1*tt20*tt16*coef(8,1)+2*tt65*tt57*coef(7,1)+&
&tt93+tt92
tt147 = tt90+tt89+2*Dm(1,1)*tt16*tt57*tt24*coef(20,1)+2*tt1*tt16*&
&tt94*tt24*coef(20,1)+tt88+2*Dm(1,2)*tt22*tt57*tt18*coef(20,1)+2*t&
&t4*tt22*tt94*tt18*coef(20,1)+4*Dm(1,1)*tt4*tt16*tt22*tt28*coef(20&
&,1)+4*Dm(1,2)*tt1*tt16*tt22*tt28*coef(20,1)+tt87+tt86+2*tt94*tt57&
&*tt24*coef(19,1)+4*Dm(1,2)*tt22*tt57*tt28*coef(19,1)+4*tt4*tt22*t&
&t94*tt28*coef(19,1)+tt85+tt84+2*tt94*tt57*tt18*coef(18,1)+4*Dm(1,&
&1)*tt16*tt57*tt28*coef(18,1)+4*tt1*tt16*tt94*tt28*coef(18,1)+tt83&
&+tt82+4*Dm(1,2)*tt22*tt57*tt24*coef(17,1)+4*tt4*tt22*tt94*tt24*co&
&ef(17,1)+8*Dm(1,2)*tt4*tt23*tt28*coef(17,1)+tt81+tt80+8*Dm(1,1)*t&
&t4*tt16*tt22*tt24*coef(16,1)+8*Dm(1,2)*tt1*tt16*tt22*tt24*coef(16&
&,1)+8*Dm(1,2)*tt4*tt23*tt18*coef(16,1)+tt79+tt78+4*Dm(1,1)*tt16*t&
&t57*tt18*coef(15,1)+4*tt1*tt16*tt94*tt18*coef(15,1)+8*Dm(1,1)*tt1&
&*tt17*tt28*coef(15,1)+tt77+tt76+8*Dm(1,1)*tt1*tt17*tt24*coef(14,1&
&)+8*Dm(1,1)*tt4*tt16*tt22*tt18*coef(14,1)+8*Dm(1,2)*tt1*tt16*tt22&
&*tt18*coef(14,1)+tt75+6*tt94*tt57*tt28*coef(13,1)+tt74+24*Dm(1,2)&
&*tt4*tt23*tt24*coef(12,1)+tt73+24*Dm(1,1)*tt1*tt17*tt18*coef(11,1&
&)+tt72+tt71+2*Dm(1,1)*tt16*tt57*coef(10,1)+2*tt1*tt16*tt94*coef(1&
&0,1)+tt70+tt69+2*Dm(1,2)*tt22*tt57*coef(9,1)+2*tt4*tt22*tt94*coef&
&(9,1)+tt68+tt67+4*Dm(1,1)*tt4*tt16*tt22*coef(8,1)+4*Dm(1,2)*tt1*t&
&t16*tt22*coef(8,1)+tt66+2*tt94*tt57*coef(7,1)+tt64+8*Dm(1,2)*tt4*&
&tt23*coef(6,1)+tt63+8*Dm(1,1)*tt1*tt17*coef(5,1)+tt62+tt60+tt59
tt148 = 2*Dm(2,1)*tt11*tt57*tt24*coef(20,1)+2*tt1*tt108*tt16*tt24&
&*coef(20,1)+2*Dm(2,2)*tt20*tt57*tt18*coef(20,1)+2*tt4*tt108*tt22*&
&tt18*coef(20,1)+4*Dm(2,1)*tt11*tt4*tt22*tt28*coef(20,1)+4*tt1*Dm(&
&2,2)*tt20*tt16*tt28*coef(20,1)+2*tt108*tt57*tt24*coef(19,1)+4*Dm(&
&2,2)*tt20*tt57*tt28*coef(19,1)+4*tt4*tt108*tt22*tt28*coef(19,1)+2&
&*tt108*tt57*tt18*coef(18,1)+4*Dm(2,1)*tt11*tt57*tt28*coef(18,1)+4&
&*tt1*tt108*tt16*tt28*coef(18,1)+4*Dm(2,2)*tt20*tt57*tt24*coef(17,&
&1)+4*tt4*tt108*tt22*tt24*coef(17,1)+tt143+8*Dm(2,1)*tt11*tt4*tt22&
&*tt24*coef(16,1)+8*tt1*Dm(2,2)*tt20*tt16*tt24*coef(16,1)+tt142+4*&
&Dm(2,1)*tt11*tt57*tt18*coef(15,1)+4*tt1*tt108*tt16*tt18*coef(15,1&
&)+tt141+tt140+8*Dm(2,1)*tt11*tt4*tt22*tt18*coef(14,1)+8*tt1*Dm(2,&
&2)*tt20*tt16*tt18*coef(14,1)+6*tt108*tt57*tt28*coef(13,1)+tt139+t&
&t138+2*Dm(2,1)*tt11*tt57*coef(10,1)+2*tt1*tt108*tt16*coef(10,1)+2&
&*Dm(2,2)*tt20*tt57*coef(9,1)+2*tt4*tt108*tt22*coef(9,1)+4*Dm(2,1)&
&*tt11*tt4*tt22*coef(8,1)+4*tt1*Dm(2,2)*tt20*tt16*coef(8,1)+2*tt10&
&8*tt57*coef(7,1)+tt136+tt135
tt149 = tt133+tt132+2*tt1*tt16*tt137*tt24*coef(20,1)+2*Dm(2,1)*tt&
&16*tt57*tt24*coef(20,1)+tt131+2*tt4*tt22*tt137*tt18*coef(20,1)+2*&
&Dm(2,2)*tt22*tt57*tt18*coef(20,1)+4*tt1*Dm(2,2)*tt16*tt22*tt28*co&
&ef(20,1)+4*Dm(2,1)*tt4*tt16*tt22*tt28*coef(20,1)+tt130+tt129+2*tt&
&57*tt137*tt24*coef(19,1)+4*tt4*tt22*tt137*tt28*coef(19,1)+4*Dm(2,&
&2)*tt22*tt57*tt28*coef(19,1)+tt128+tt127+2*tt57*tt137*tt18*coef(1&
&8,1)+4*tt1*tt16*tt137*tt28*coef(18,1)+4*Dm(2,1)*tt16*tt57*tt28*co&
&ef(18,1)+tt126+tt125+4*tt4*tt22*tt137*tt24*coef(17,1)+4*Dm(2,2)*t&
&t22*tt57*tt24*coef(17,1)+8*tt4*Dm(2,2)*tt23*tt28*coef(17,1)+tt124&
&+tt123+8*tt1*Dm(2,2)*tt16*tt22*tt24*coef(16,1)+8*Dm(2,1)*tt4*tt16&
&*tt22*tt24*coef(16,1)+8*tt4*Dm(2,2)*tt23*tt18*coef(16,1)+tt122+tt&
&121+4*tt1*tt16*tt137*tt18*coef(15,1)+4*Dm(2,1)*tt16*tt57*tt18*coe&
&f(15,1)+8*tt1*Dm(2,1)*tt17*tt28*coef(15,1)+tt120+tt119+8*tt1*Dm(2&
&,1)*tt17*tt24*coef(14,1)+8*tt1*Dm(2,2)*tt16*tt22*tt18*coef(14,1)+&
&8*Dm(2,1)*tt4*tt16*tt22*tt18*coef(14,1)+tt118+6*tt57*tt137*tt28*c&
&oef(13,1)+tt117+24*tt4*Dm(2,2)*tt23*tt24*coef(12,1)+tt116+24*tt1*&
&Dm(2,1)*tt17*tt18*coef(11,1)+tt115+tt114+2*tt1*tt16*tt137*coef(10&
&,1)+2*Dm(2,1)*tt16*tt57*coef(10,1)+tt113+tt112+2*tt4*tt22*tt137*c&
&oef(9,1)+2*Dm(2,2)*tt22*tt57*coef(9,1)+tt111+tt110+4*tt1*Dm(2,2)*&
&tt16*tt22*coef(8,1)+4*Dm(2,1)*tt4*tt16*tt22*coef(8,1)+tt109+2*tt5&
&7*tt137*coef(7,1)+tt107+8*tt4*Dm(2,2)*tt23*coef(6,1)+tt106+8*tt1*&
&Dm(2,1)*tt17*coef(5,1)+tt105+tt103+tt102
tt150 = Dm(1,1)**2
tt151 = 2*tt150*coef(2,1)
tt152 = Dm(1,2)**2
tt153 = 2*tt152*coef(3,1)
tt154 = 2*Dm(1,1)*Dm(1,2)*coef(4,1)
tt155 = 4*tt150*tt18*coef(5,1)
tt156 = 4*tt152*tt24*coef(6,1)
tt157 = tt65**2
tt158 = 4*Dm(1,1)*Dm(1,2)*tt28*coef(7,1)
tt159 = 2*tt152*tt18*coef(8,1)
tt160 = 2*tt150*tt24*coef(8,1)
tt161 = 2*tt152*tt28*coef(9,1)
tt162 = 2*Dm(1,1)*Dm(1,2)*tt24*coef(9,1)
tt163 = 2*tt150*tt28*coef(10,1)
tt164 = 2*Dm(1,1)*Dm(1,2)*tt18*coef(10,1)
tt165 = 6*tt150*tt36*coef(11,1)
tt166 = 6*tt152*tt38*coef(12,1)
tt167 = 6*Dm(1,1)*Dm(1,2)*tt40*coef(13,1)
tt168 = 4*tt150*tt18*tt24*coef(14,1)
tt169 = 2*tt152*tt36*coef(14,1)
tt170 = 4*tt150*tt28*tt18*coef(15,1)
tt171 = 2*Dm(1,1)*Dm(1,2)*tt36*coef(15,1)
tt172 = 4*tt152*tt18*tt24*coef(16,1)
tt173 = 2*tt150*tt38*coef(16,1)
tt174 = 4*tt152*tt28*tt24*coef(17,1)
tt175 = 2*Dm(1,1)*Dm(1,2)*tt38*coef(17,1)
tt176 = 4*Dm(1,1)*Dm(1,2)*tt28*tt18*coef(18,1)
tt177 = 2*tt150*tt40*coef(18,1)
tt178 = 4*Dm(1,1)*Dm(1,2)*tt28*tt24*coef(19,1)
tt179 = 2*tt152*tt40*coef(19,1)
tt180 = 2*tt152*tt28*tt18*coef(20,1)
tt181 = 2*tt150*tt28*tt24*coef(20,1)
tt182 = 2*Dm(1,1)*Dm(1,2)*tt18*tt24*coef(20,1)
tt183 = 2*Dm(1,1)*tt11*tt94*tt24*coef(20,1)+2*Dm(1,1)*tt65*tt16*t&
&t24*coef(20,1)+2*Dm(1,2)*tt20*tt94*tt18*coef(20,1)+2*Dm(1,2)*tt65&
&*tt22*tt18*coef(20,1)+4*Dm(1,1)*Dm(1,2)*tt11*tt22*tt28*coef(20,1)&
&+4*Dm(1,1)*Dm(1,2)*tt20*tt16*tt28*coef(20,1)+2*tt65*tt94*tt24*coe&
&f(19,1)+4*Dm(1,2)*tt20*tt94*tt28*coef(19,1)+4*Dm(1,2)*tt65*tt22*t&
&t28*coef(19,1)+2*tt65*tt94*tt18*coef(18,1)+4*Dm(1,1)*tt11*tt94*tt&
&28*coef(18,1)+4*Dm(1,1)*tt65*tt16*tt28*coef(18,1)+4*Dm(1,2)*tt20*&
&tt94*tt24*coef(17,1)+4*Dm(1,2)*tt65*tt22*tt24*coef(17,1)+8*tt152*&
&tt20*tt22*tt28*coef(17,1)+8*Dm(1,1)*Dm(1,2)*tt11*tt22*tt24*coef(1&
&6,1)+8*Dm(1,1)*Dm(1,2)*tt20*tt16*tt24*coef(16,1)+8*tt152*tt20*tt2&
&2*tt18*coef(16,1)+4*Dm(1,1)*tt11*tt94*tt18*coef(15,1)+4*Dm(1,1)*t&
&t65*tt16*tt18*coef(15,1)+8*tt150*tt11*tt16*tt28*coef(15,1)+8*tt15&
&0*tt11*tt16*tt24*coef(14,1)+8*Dm(1,1)*Dm(1,2)*tt11*tt22*tt18*coef&
&(14,1)+8*Dm(1,1)*Dm(1,2)*tt20*tt16*tt18*coef(14,1)+6*tt65*tt94*tt&
&28*coef(13,1)+24*tt152*tt20*tt22*tt24*coef(12,1)+24*tt150*tt11*tt&
&16*tt18*coef(11,1)+2*Dm(1,1)*tt11*tt94*coef(10,1)+2*Dm(1,1)*tt65*&
&tt16*coef(10,1)+2*Dm(1,2)*tt20*tt94*coef(9,1)+2*Dm(1,2)*tt65*tt22&
&*coef(9,1)+4*Dm(1,1)*Dm(1,2)*tt11*tt22*coef(8,1)+4*Dm(1,1)*Dm(1,2&
&)*tt20*tt16*coef(8,1)+2*tt65*tt94*coef(7,1)+8*tt152*tt20*tt22*coe&
&f(6,1)+8*tt150*tt11*tt16*coef(5,1)
tt184 = 2*Dm(1,1)*Dm(2,1)*coef(2,1)
tt185 = 2*Dm(1,2)*Dm(2,2)*coef(3,1)
tt186 = Dm(1,1)*Dm(2,2)+Dm(1,2)*Dm(2,1)
tt187 = tt186*coef(4,1)
tt188 = 4*Dm(1,1)*Dm(2,1)*tt18*coef(5,1)
tt189 = 4*Dm(1,2)*Dm(2,2)*tt24*coef(6,1)
tt190 = 2*tt186*tt28*coef(7,1)
tt191 = 2*Dm(1,2)*Dm(2,2)*tt18*coef(8,1)
tt192 = 2*Dm(1,1)*Dm(2,1)*tt24*coef(8,1)
tt193 = 2*Dm(1,2)*Dm(2,2)*tt28*coef(9,1)
tt194 = tt186*tt24*coef(9,1)
tt195 = 2*Dm(1,1)*Dm(2,1)*tt28*coef(10,1)
tt196 = tt186*tt18*coef(10,1)
tt197 = 6*Dm(1,1)*Dm(2,1)*tt36*coef(11,1)
tt198 = 6*Dm(1,2)*Dm(2,2)*tt38*coef(12,1)
tt199 = 3*tt186*tt40*coef(13,1)
tt200 = 4*Dm(1,1)*Dm(2,1)*tt18*tt24*coef(14,1)
tt201 = 2*Dm(1,2)*Dm(2,2)*tt36*coef(14,1)
tt202 = 4*Dm(1,1)*Dm(2,1)*tt28*tt18*coef(15,1)
tt203 = tt186*tt36*coef(15,1)
tt204 = 4*Dm(1,2)*Dm(2,2)*tt18*tt24*coef(16,1)
tt205 = 2*Dm(1,1)*Dm(2,1)*tt38*coef(16,1)
tt206 = 4*Dm(1,2)*Dm(2,2)*tt28*tt24*coef(17,1)
tt207 = tt186*tt38*coef(17,1)
tt208 = 2*tt186*tt28*tt18*coef(18,1)
tt209 = 2*Dm(1,1)*Dm(2,1)*tt40*coef(18,1)
tt210 = 2*tt186*tt28*tt24*coef(19,1)
tt211 = 2*Dm(1,2)*Dm(2,2)*tt40*coef(19,1)
tt212 = 2*Dm(1,2)*Dm(2,2)*tt28*tt18*coef(20,1)
tt213 = 2*Dm(1,1)*Dm(2,1)*tt28*tt24*coef(20,1)
tt214 = tt186*tt18*tt24*coef(20,1)
tt215 = tt214+tt213+2*Dm(1,1)*tt11*tt108*tt24*coef(20,1)+2*Dm(2,1&
&)*tt11*tt65*tt24*coef(20,1)+tt212+2*Dm(1,2)*tt20*tt108*tt18*coef(&
&20,1)+2*Dm(2,2)*tt20*tt65*tt18*coef(20,1)+4*Dm(1,1)*tt11*Dm(2,2)*&
&tt20*tt28*coef(20,1)+4*Dm(1,2)*Dm(2,1)*tt11*tt20*tt28*coef(20,1)+&
&tt211+tt210+2*tt65*tt108*tt24*coef(19,1)+4*Dm(1,2)*tt20*tt108*tt2&
&8*coef(19,1)+4*Dm(2,2)*tt20*tt65*tt28*coef(19,1)+tt209+tt208+2*tt&
&65*tt108*tt18*coef(18,1)+4*Dm(1,1)*tt11*tt108*tt28*coef(18,1)+4*D&
&m(2,1)*tt11*tt65*tt28*coef(18,1)+tt207+tt206+4*Dm(1,2)*tt20*tt108&
&*tt24*coef(17,1)+4*Dm(2,2)*tt20*tt65*tt24*coef(17,1)+8*Dm(1,2)*Dm&
&(2,2)*tt21*tt28*coef(17,1)+tt205+tt204+8*Dm(1,1)*tt11*Dm(2,2)*tt2&
&0*tt24*coef(16,1)+8*Dm(1,2)*Dm(2,1)*tt11*tt20*tt24*coef(16,1)+8*D&
&m(1,2)*Dm(2,2)*tt21*tt18*coef(16,1)+tt203+tt202+4*Dm(1,1)*tt11*tt&
&108*tt18*coef(15,1)+4*Dm(2,1)*tt11*tt65*tt18*coef(15,1)+8*Dm(1,1)&
&*Dm(2,1)*tt12*tt28*coef(15,1)+tt201+tt200+8*Dm(1,1)*Dm(2,1)*tt12*&
&tt24*coef(14,1)+8*Dm(1,1)*tt11*Dm(2,2)*tt20*tt18*coef(14,1)+8*Dm(&
&1,2)*Dm(2,1)*tt11*tt20*tt18*coef(14,1)+tt199+6*tt65*tt108*tt28*co&
&ef(13,1)+tt198+24*Dm(1,2)*Dm(2,2)*tt21*tt24*coef(12,1)+tt197+24*D&
&m(1,1)*Dm(2,1)*tt12*tt18*coef(11,1)+tt196+tt195+2*Dm(1,1)*tt11*tt&
&108*coef(10,1)+2*Dm(2,1)*tt11*tt65*coef(10,1)+tt194+tt193+2*Dm(1,&
&2)*tt20*tt108*coef(9,1)+2*Dm(2,2)*tt20*tt65*coef(9,1)+tt192+tt191&
&+4*Dm(1,1)*tt11*Dm(2,2)*tt20*coef(8,1)+4*Dm(1,2)*Dm(2,1)*tt11*tt2&
&0*coef(8,1)+tt190+2*tt65*tt108*coef(7,1)+tt189+8*Dm(1,2)*Dm(2,2)*&
&tt21*coef(6,1)+tt188+8*Dm(1,1)*Dm(2,1)*tt12*coef(5,1)+tt187+tt185&
&+tt184
tt216 = 8*Dm(1,1)*Dm(2,1)*tt11*tt16*coef(5,1)
tt217 = 8*Dm(1,2)*Dm(2,2)*tt20*tt22*coef(6,1)
tt218 = 24*Dm(1,1)*Dm(2,1)*tt11*tt16*tt18*coef(11,1)
tt219 = 24*Dm(1,2)*Dm(2,2)*tt20*tt22*tt24*coef(12,1)
tt220 = 8*Dm(1,1)*Dm(2,1)*tt11*tt16*tt24*coef(14,1)
tt221 = 8*Dm(1,1)*Dm(2,1)*tt11*tt16*tt28*coef(15,1)
tt222 = 8*Dm(1,2)*Dm(2,2)*tt20*tt22*tt18*coef(16,1)
tt223 = 8*Dm(1,2)*Dm(2,2)*tt20*tt22*tt28*coef(17,1)
tt224 = 2*Dm(1,1)*tt11*tt137*tt24*coef(20,1)+2*Dm(2,1)*tt65*tt16*&
&tt24*coef(20,1)+2*Dm(1,2)*tt20*tt137*tt18*coef(20,1)+2*Dm(2,2)*tt&
&65*tt22*tt18*coef(20,1)+4*Dm(1,1)*tt11*Dm(2,2)*tt22*tt28*coef(20,&
&1)+4*Dm(1,2)*Dm(2,1)*tt20*tt16*tt28*coef(20,1)+2*tt65*tt137*tt24*&
&coef(19,1)+4*Dm(1,2)*tt20*tt137*tt28*coef(19,1)+4*Dm(2,2)*tt65*tt&
&22*tt28*coef(19,1)+2*tt65*tt137*tt18*coef(18,1)+4*Dm(1,1)*tt11*tt&
&137*tt28*coef(18,1)+4*Dm(2,1)*tt65*tt16*tt28*coef(18,1)+4*Dm(1,2)&
&*tt20*tt137*tt24*coef(17,1)+4*Dm(2,2)*tt65*tt22*tt24*coef(17,1)+t&
&t223+8*Dm(1,1)*tt11*Dm(2,2)*tt22*tt24*coef(16,1)+8*Dm(1,2)*Dm(2,1&
&)*tt20*tt16*tt24*coef(16,1)+tt222+4*Dm(1,1)*tt11*tt137*tt18*coef(&
&15,1)+4*Dm(2,1)*tt65*tt16*tt18*coef(15,1)+tt221+tt220+8*Dm(1,1)*t&
&t11*Dm(2,2)*tt22*tt18*coef(14,1)+8*Dm(1,2)*Dm(2,1)*tt20*tt16*tt18&
&*coef(14,1)+6*tt65*tt137*tt28*coef(13,1)+tt219+tt218+2*Dm(1,1)*tt&
&11*tt137*coef(10,1)+2*Dm(2,1)*tt65*tt16*coef(10,1)+2*Dm(1,2)*tt20&
&*tt137*coef(9,1)+2*Dm(2,2)*tt65*tt22*coef(9,1)+4*Dm(1,1)*tt11*Dm(&
&2,2)*tt22*coef(8,1)+4*Dm(1,2)*Dm(2,1)*tt20*tt16*coef(8,1)+2*tt65*&
&tt137*coef(7,1)+tt217+tt216
tt225 = tt94**2
tt226 = 2*Dm(2,1)*tt11*tt94*tt24*coef(20,1)+2*Dm(1,1)*tt108*tt16*&
&tt24*coef(20,1)+2*Dm(2,2)*tt20*tt94*tt18*coef(20,1)+2*Dm(1,2)*tt1&
&08*tt22*tt18*coef(20,1)+4*Dm(1,2)*Dm(2,1)*tt11*tt22*tt28*coef(20,&
&1)+4*Dm(1,1)*Dm(2,2)*tt20*tt16*tt28*coef(20,1)+2*tt108*tt94*tt24*&
&coef(19,1)+4*Dm(2,2)*tt20*tt94*tt28*coef(19,1)+4*Dm(1,2)*tt108*tt&
&22*tt28*coef(19,1)+2*tt108*tt94*tt18*coef(18,1)+4*Dm(2,1)*tt11*tt&
&94*tt28*coef(18,1)+4*Dm(1,1)*tt108*tt16*tt28*coef(18,1)+4*Dm(2,2)&
&*tt20*tt94*tt24*coef(17,1)+4*Dm(1,2)*tt108*tt22*tt24*coef(17,1)+t&
&t223+8*Dm(1,2)*Dm(2,1)*tt11*tt22*tt24*coef(16,1)+8*Dm(1,1)*Dm(2,2&
&)*tt20*tt16*tt24*coef(16,1)+tt222+4*Dm(2,1)*tt11*tt94*tt18*coef(1&
&5,1)+4*Dm(1,1)*tt108*tt16*tt18*coef(15,1)+tt221+tt220+8*Dm(1,2)*D&
&m(2,1)*tt11*tt22*tt18*coef(14,1)+8*Dm(1,1)*Dm(2,2)*tt20*tt16*tt18&
&*coef(14,1)+6*tt108*tt94*tt28*coef(13,1)+tt219+tt218+2*Dm(2,1)*tt&
&11*tt94*coef(10,1)+2*Dm(1,1)*tt108*tt16*coef(10,1)+2*Dm(2,2)*tt20&
&*tt94*coef(9,1)+2*Dm(1,2)*tt108*tt22*coef(9,1)+4*Dm(1,2)*Dm(2,1)*&
&tt11*tt22*coef(8,1)+4*Dm(1,1)*Dm(2,2)*tt20*tt16*coef(8,1)+2*tt108&
&*tt94*coef(7,1)+tt217+tt216
tt227 = tt214+tt213+2*Dm(1,1)*tt16*tt137*tt24*coef(20,1)+2*Dm(2,1&
&)*tt16*tt94*tt24*coef(20,1)+tt212+2*Dm(1,2)*tt22*tt137*tt18*coef(&
&20,1)+2*Dm(2,2)*tt22*tt94*tt18*coef(20,1)+4*Dm(1,1)*Dm(2,2)*tt16*&
&tt22*tt28*coef(20,1)+4*Dm(1,2)*Dm(2,1)*tt16*tt22*tt28*coef(20,1)+&
&tt211+tt210+2*tt94*tt137*tt24*coef(19,1)+4*Dm(1,2)*tt22*tt137*tt2&
&8*coef(19,1)+4*Dm(2,2)*tt22*tt94*tt28*coef(19,1)+tt209+tt208+2*tt&
&94*tt137*tt18*coef(18,1)+4*Dm(1,1)*tt16*tt137*tt28*coef(18,1)+4*D&
&m(2,1)*tt16*tt94*tt28*coef(18,1)+tt207+tt206+4*Dm(1,2)*tt22*tt137&
&*tt24*coef(17,1)+4*Dm(2,2)*tt22*tt94*tt24*coef(17,1)+8*Dm(1,2)*Dm&
&(2,2)*tt23*tt28*coef(17,1)+tt205+tt204+8*Dm(1,1)*Dm(2,2)*tt16*tt2&
&2*tt24*coef(16,1)+8*Dm(1,2)*Dm(2,1)*tt16*tt22*tt24*coef(16,1)+8*D&
&m(1,2)*Dm(2,2)*tt23*tt18*coef(16,1)+tt203+tt202+4*Dm(1,1)*tt16*tt&
&137*tt18*coef(15,1)+4*Dm(2,1)*tt16*tt94*tt18*coef(15,1)+8*Dm(1,1)&
&*Dm(2,1)*tt17*tt28*coef(15,1)+tt201+tt200+8*Dm(1,1)*Dm(2,1)*tt17*&
&tt24*coef(14,1)+8*Dm(1,1)*Dm(2,2)*tt16*tt22*tt18*coef(14,1)+8*Dm(&
&1,2)*Dm(2,1)*tt16*tt22*tt18*coef(14,1)+tt199+6*tt94*tt137*tt28*co&
&ef(13,1)+tt198+24*Dm(1,2)*Dm(2,2)*tt23*tt24*coef(12,1)+tt197+24*D&
&m(1,1)*Dm(2,1)*tt17*tt18*coef(11,1)+tt196+tt195+2*Dm(1,1)*tt16*tt&
&137*coef(10,1)+2*Dm(2,1)*tt16*tt94*coef(10,1)+tt194+tt193+2*Dm(1,&
&2)*tt22*tt137*coef(9,1)+2*Dm(2,2)*tt22*tt94*coef(9,1)+tt192+tt191&
&+4*Dm(1,1)*Dm(2,2)*tt16*tt22*coef(8,1)+4*Dm(1,2)*Dm(2,1)*tt16*tt2&
&2*coef(8,1)+tt190+2*tt94*tt137*coef(7,1)+tt189+8*Dm(1,2)*Dm(2,2)*&
&tt23*coef(6,1)+tt188+8*Dm(1,1)*Dm(2,1)*tt17*coef(5,1)+tt187+tt185&
&+tt184
tt228 = Dm(2,1)**2
tt229 = 2*tt228*coef(2,1)
tt230 = Dm(2,2)**2
tt231 = 2*tt230*coef(3,1)
tt232 = 2*Dm(2,1)*Dm(2,2)*coef(4,1)
tt233 = 4*tt228*tt18*coef(5,1)
tt234 = 4*tt230*tt24*coef(6,1)
tt235 = tt108**2
tt236 = 4*Dm(2,1)*Dm(2,2)*tt28*coef(7,1)
tt237 = 2*tt230*tt18*coef(8,1)
tt238 = 2*tt228*tt24*coef(8,1)
tt239 = 2*tt230*tt28*coef(9,1)
tt240 = 2*Dm(2,1)*Dm(2,2)*tt24*coef(9,1)
tt241 = 2*tt228*tt28*coef(10,1)
tt242 = 2*Dm(2,1)*Dm(2,2)*tt18*coef(10,1)
tt243 = 6*tt228*tt36*coef(11,1)
tt244 = 6*tt230*tt38*coef(12,1)
tt245 = 6*Dm(2,1)*Dm(2,2)*tt40*coef(13,1)
tt246 = 4*tt228*tt18*tt24*coef(14,1)
tt247 = 2*tt230*tt36*coef(14,1)
tt248 = 4*tt228*tt28*tt18*coef(15,1)
tt249 = 2*Dm(2,1)*Dm(2,2)*tt36*coef(15,1)
tt250 = 4*tt230*tt18*tt24*coef(16,1)
tt251 = 2*tt228*tt38*coef(16,1)
tt252 = 4*tt230*tt28*tt24*coef(17,1)
tt253 = 2*Dm(2,1)*Dm(2,2)*tt38*coef(17,1)
tt254 = 4*Dm(2,1)*Dm(2,2)*tt28*tt18*coef(18,1)
tt255 = 2*tt228*tt40*coef(18,1)
tt256 = 4*Dm(2,1)*Dm(2,2)*tt28*tt24*coef(19,1)
tt257 = 2*tt230*tt40*coef(19,1)
tt258 = 2*tt230*tt28*tt18*coef(20,1)
tt259 = 2*tt228*tt28*tt24*coef(20,1)
tt260 = 2*Dm(2,1)*Dm(2,2)*tt18*tt24*coef(20,1)
tt261 = 2*Dm(2,1)*tt11*tt137*tt24*coef(20,1)+2*Dm(2,1)*tt108*tt16&
&*tt24*coef(20,1)+2*Dm(2,2)*tt20*tt137*tt18*coef(20,1)+2*Dm(2,2)*t&
&t108*tt22*tt18*coef(20,1)+4*Dm(2,1)*tt11*Dm(2,2)*tt22*tt28*coef(2&
&0,1)+4*Dm(2,1)*Dm(2,2)*tt20*tt16*tt28*coef(20,1)+2*tt108*tt137*tt&
&24*coef(19,1)+4*Dm(2,2)*tt20*tt137*tt28*coef(19,1)+4*Dm(2,2)*tt10&
&8*tt22*tt28*coef(19,1)+2*tt108*tt137*tt18*coef(18,1)+4*Dm(2,1)*tt&
&11*tt137*tt28*coef(18,1)+4*Dm(2,1)*tt108*tt16*tt28*coef(18,1)+4*D&
&m(2,2)*tt20*tt137*tt24*coef(17,1)+4*Dm(2,2)*tt108*tt22*tt24*coef(&
&17,1)+8*tt230*tt20*tt22*tt28*coef(17,1)+8*Dm(2,1)*tt11*Dm(2,2)*tt&
&22*tt24*coef(16,1)+8*Dm(2,1)*Dm(2,2)*tt20*tt16*tt24*coef(16,1)+8*&
&tt230*tt20*tt22*tt18*coef(16,1)+4*Dm(2,1)*tt11*tt137*tt18*coef(15&
&,1)+4*Dm(2,1)*tt108*tt16*tt18*coef(15,1)+8*tt228*tt11*tt16*tt28*c&
&oef(15,1)+8*tt228*tt11*tt16*tt24*coef(14,1)+8*Dm(2,1)*tt11*Dm(2,2&
&)*tt22*tt18*coef(14,1)+8*Dm(2,1)*Dm(2,2)*tt20*tt16*tt18*coef(14,1&
&)+6*tt108*tt137*tt28*coef(13,1)+24*tt230*tt20*tt22*tt24*coef(12,1&
&)+24*tt228*tt11*tt16*tt18*coef(11,1)+2*Dm(2,1)*tt11*tt137*coef(10&
&,1)+2*Dm(2,1)*tt108*tt16*coef(10,1)+2*Dm(2,2)*tt20*tt137*coef(9,1&
&)+2*Dm(2,2)*tt108*tt22*coef(9,1)+4*Dm(2,1)*tt11*Dm(2,2)*tt22*coef&
&(8,1)+4*Dm(2,1)*Dm(2,2)*tt20*tt16*coef(8,1)+2*tt108*tt137*coef(7,&
&1)+8*tt230*tt20*tt22*coef(6,1)+8*tt228*tt11*tt16*coef(5,1)
tt262 = tt137**2
hes(1,1) = tt56+tt55+4*tt1*tt11*tt26*tt24*coef(20,1)+tt54+4*tt4*t&
&t20*tt26*tt18*coef(20,1)+8*tt1*tt11*tt4*tt20*tt28*coef(20,1)+tt53&
&+tt52+2*tt27*tt24*coef(19,1)+8*tt4*tt20*tt26*tt28*coef(19,1)+tt51&
&+tt50+2*tt27*tt18*coef(18,1)+8*tt1*tt11*tt26*tt28*coef(18,1)+tt49&
&+tt48+8*tt4*tt20*tt26*tt24*coef(17,1)+8*tt5*tt21*tt28*coef(17,1)+&
&tt47+tt46+16*tt1*tt11*tt4*tt20*tt24*coef(16,1)+8*tt5*tt21*tt18*co&
&ef(16,1)+tt45+tt44+8*tt1*tt11*tt26*tt18*coef(15,1)+8*tt2*tt12*tt2&
&8*coef(15,1)+tt43+tt42+8*tt2*tt12*tt24*coef(14,1)+16*tt1*tt11*tt4&
&*tt20*tt18*coef(14,1)+tt41+6*tt27*tt28*coef(13,1)+tt39+24*tt5*tt2&
&1*tt24*coef(12,1)+tt37+24*tt2*tt12*tt18*coef(11,1)+tt35+tt34+4*tt&
&1*tt11*tt26*coef(10,1)+tt33+tt32+4*tt4*tt20*tt26*coef(9,1)+tt31+t&
&t30+8*tt1*tt11*tt4*tt20*coef(8,1)+tt29+2*tt27*coef(7,1)+tt25+8*tt&
&5*tt21*coef(6,1)+tt19+8*tt2*tt12*coef(5,1)+tt7+tt6+tt3
hes(1,2) = tt58
hes(1,3) = tt91
hes(1,4) = tt101
hes(1,5) = tt134
hes(1,6) = tt144
hes(2,1) = tt58
hes(2,2) = tt56+tt55+4*tt1*tt16*tt57*tt24*coef(20,1)+tt54+4*tt4*t&
&t22*tt57*tt18*coef(20,1)+8*tt1*tt4*tt16*tt22*tt28*coef(20,1)+tt53&
&+2*tt145*tt24*coef(19,1)+tt52+8*tt4*tt22*tt57*tt28*coef(19,1)+tt5&
&1+2*tt145*tt18*coef(18,1)+tt50+8*tt1*tt16*tt57*tt28*coef(18,1)+tt&
&49+tt48+8*tt4*tt22*tt57*tt24*coef(17,1)+8*tt5*tt23*tt28*coef(17,1&
&)+tt47+tt46+16*tt1*tt4*tt16*tt22*tt24*coef(16,1)+8*tt5*tt23*tt18*&
&coef(16,1)+tt45+tt44+8*tt1*tt16*tt57*tt18*coef(15,1)+8*tt2*tt17*t&
&t28*coef(15,1)+tt43+tt42+8*tt2*tt17*tt24*coef(14,1)+16*tt1*tt4*tt&
&16*tt22*tt18*coef(14,1)+tt41+6*tt145*tt28*coef(13,1)+tt39+24*tt5*&
&tt23*tt24*coef(12,1)+tt37+24*tt2*tt17*tt18*coef(11,1)+tt35+tt34+4&
&*tt1*tt16*tt57*coef(10,1)+tt33+tt32+4*tt4*tt22*tt57*coef(9,1)+tt3&
&1+tt30+8*tt1*tt4*tt16*tt22*coef(8,1)+2*tt145*coef(7,1)+tt29+tt25+&
&8*tt5*tt23*coef(6,1)+tt19+8*tt2*tt17*coef(5,1)+tt7+tt6+tt3
hes(2,3) = tt146
hes(2,4) = tt147
hes(2,5) = tt148
hes(2,6) = tt149
hes(3,1) = tt91
hes(3,2) = tt146
hes(3,3) = tt182+tt181+4*Dm(1,1)*tt11*tt65*tt24*coef(20,1)+tt180+&
&4*Dm(1,2)*tt20*tt65*tt18*coef(20,1)+8*Dm(1,1)*Dm(1,2)*tt11*tt20*t&
&t28*coef(20,1)+tt179+tt178+2*tt157*tt24*coef(19,1)+8*Dm(1,2)*tt20&
&*tt65*tt28*coef(19,1)+tt177+tt176+2*tt157*tt18*coef(18,1)+8*Dm(1,&
&1)*tt11*tt65*tt28*coef(18,1)+tt175+tt174+8*Dm(1,2)*tt20*tt65*tt24&
&*coef(17,1)+8*tt152*tt21*tt28*coef(17,1)+tt173+tt172+16*Dm(1,1)*D&
&m(1,2)*tt11*tt20*tt24*coef(16,1)+8*tt152*tt21*tt18*coef(16,1)+tt1&
&71+tt170+8*Dm(1,1)*tt11*tt65*tt18*coef(15,1)+8*tt150*tt12*tt28*co&
&ef(15,1)+tt169+tt168+8*tt150*tt12*tt24*coef(14,1)+16*Dm(1,1)*Dm(1&
&,2)*tt11*tt20*tt18*coef(14,1)+tt167+6*tt157*tt28*coef(13,1)+tt166&
&+24*tt152*tt21*tt24*coef(12,1)+tt165+24*tt150*tt12*tt18*coef(11,1&
&)+tt164+tt163+4*Dm(1,1)*tt11*tt65*coef(10,1)+tt162+tt161+4*Dm(1,2&
&)*tt20*tt65*coef(9,1)+tt160+tt159+8*Dm(1,1)*Dm(1,2)*tt11*tt20*coe&
&f(8,1)+tt158+2*tt157*coef(7,1)+tt156+8*tt152*tt21*coef(6,1)+tt155&
&+8*tt150*tt12*coef(5,1)+tt154+tt153+tt151
hes(3,4) = tt183
hes(3,5) = tt215
hes(3,6) = tt224
hes(4,1) = tt101
hes(4,2) = tt147
hes(4,3) = tt183
hes(4,4) = tt182+tt181+4*Dm(1,1)*tt16*tt94*tt24*coef(20,1)+tt180+&
&4*Dm(1,2)*tt22*tt94*tt18*coef(20,1)+8*Dm(1,1)*Dm(1,2)*tt16*tt22*t&
&t28*coef(20,1)+tt179+2*tt225*tt24*coef(19,1)+tt178+8*Dm(1,2)*tt22&
&*tt94*tt28*coef(19,1)+tt177+2*tt225*tt18*coef(18,1)+tt176+8*Dm(1,&
&1)*tt16*tt94*tt28*coef(18,1)+tt175+tt174+8*Dm(1,2)*tt22*tt94*tt24&
&*coef(17,1)+8*tt152*tt23*tt28*coef(17,1)+tt173+tt172+16*Dm(1,1)*D&
&m(1,2)*tt16*tt22*tt24*coef(16,1)+8*tt152*tt23*tt18*coef(16,1)+tt1&
&71+tt170+8*Dm(1,1)*tt16*tt94*tt18*coef(15,1)+8*tt150*tt17*tt28*co&
&ef(15,1)+tt169+tt168+8*tt150*tt17*tt24*coef(14,1)+16*Dm(1,1)*Dm(1&
&,2)*tt16*tt22*tt18*coef(14,1)+tt167+6*tt225*tt28*coef(13,1)+tt166&
&+24*tt152*tt23*tt24*coef(12,1)+tt165+24*tt150*tt17*tt18*coef(11,1&
&)+tt164+tt163+4*Dm(1,1)*tt16*tt94*coef(10,1)+tt162+tt161+4*Dm(1,2&
&)*tt22*tt94*coef(9,1)+tt160+tt159+8*Dm(1,1)*Dm(1,2)*tt16*tt22*coe&
&f(8,1)+2*tt225*coef(7,1)+tt158+tt156+8*tt152*tt23*coef(6,1)+tt155&
&+8*tt150*tt17*coef(5,1)+tt154+tt153+tt151
hes(4,5) = tt226
hes(4,6) = tt227
hes(5,1) = tt134
hes(5,2) = tt148
hes(5,3) = tt215
hes(5,4) = tt226
hes(5,5) = tt260+tt259+4*Dm(2,1)*tt11*tt108*tt24*coef(20,1)+tt258&
&+4*Dm(2,2)*tt20*tt108*tt18*coef(20,1)+8*Dm(2,1)*tt11*Dm(2,2)*tt20&
&*tt28*coef(20,1)+tt257+tt256+2*tt235*tt24*coef(19,1)+8*Dm(2,2)*tt&
&20*tt108*tt28*coef(19,1)+tt255+tt254+2*tt235*tt18*coef(18,1)+8*Dm&
&(2,1)*tt11*tt108*tt28*coef(18,1)+tt253+tt252+8*Dm(2,2)*tt20*tt108&
&*tt24*coef(17,1)+8*tt230*tt21*tt28*coef(17,1)+tt251+tt250+16*Dm(2&
&,1)*tt11*Dm(2,2)*tt20*tt24*coef(16,1)+8*tt230*tt21*tt18*coef(16,1&
&)+tt249+tt248+8*Dm(2,1)*tt11*tt108*tt18*coef(15,1)+8*tt228*tt12*t&
&t28*coef(15,1)+tt247+tt246+8*tt228*tt12*tt24*coef(14,1)+16*Dm(2,1&
&)*tt11*Dm(2,2)*tt20*tt18*coef(14,1)+tt245+6*tt235*tt28*coef(13,1)&
&+tt244+24*tt230*tt21*tt24*coef(12,1)+tt243+24*tt228*tt12*tt18*coe&
&f(11,1)+tt242+tt241+4*Dm(2,1)*tt11*tt108*coef(10,1)+tt240+tt239+4&
&*Dm(2,2)*tt20*tt108*coef(9,1)+tt238+tt237+8*Dm(2,1)*tt11*Dm(2,2)*&
&tt20*coef(8,1)+tt236+2*tt235*coef(7,1)+tt234+8*tt230*tt21*coef(6,&
&1)+tt233+8*tt228*tt12*coef(5,1)+tt232+tt231+tt229
hes(5,6) = tt261
hes(6,1) = tt144
hes(6,2) = tt149
hes(6,3) = tt224
hes(6,4) = tt227
hes(6,5) = tt261
hes(6,6) = tt260+tt259+4*Dm(2,1)*tt16*tt137*tt24*coef(20,1)+tt258&
&+4*Dm(2,2)*tt22*tt137*tt18*coef(20,1)+8*Dm(2,1)*Dm(2,2)*tt16*tt22&
&*tt28*coef(20,1)+tt257+2*tt262*tt24*coef(19,1)+tt256+8*Dm(2,2)*tt&
&22*tt137*tt28*coef(19,1)+tt255+2*tt262*tt18*coef(18,1)+tt254+8*Dm&
&(2,1)*tt16*tt137*tt28*coef(18,1)+tt253+tt252+8*Dm(2,2)*tt22*tt137&
&*tt24*coef(17,1)+8*tt230*tt23*tt28*coef(17,1)+tt251+tt250+16*Dm(2&
&,1)*Dm(2,2)*tt16*tt22*tt24*coef(16,1)+8*tt230*tt23*tt18*coef(16,1&
&)+tt249+tt248+8*Dm(2,1)*tt16*tt137*tt18*coef(15,1)+8*tt228*tt17*t&
&t28*coef(15,1)+tt247+tt246+8*tt228*tt17*tt24*coef(14,1)+16*Dm(2,1&
&)*Dm(2,2)*tt16*tt22*tt18*coef(14,1)+tt245+6*tt262*tt28*coef(13,1)&
&+tt244+24*tt230*tt23*tt24*coef(12,1)+tt243+24*tt228*tt17*tt18*coe&
&f(11,1)+tt242+tt241+4*Dm(2,1)*tt16*tt137*coef(10,1)+tt240+tt239+4&
&*Dm(2,2)*tt22*tt137*coef(9,1)+tt238+tt237+8*Dm(2,1)*Dm(2,2)*tt16*&
&tt22*coef(8,1)+2*tt262*coef(7,1)+tt236+tt234+8*tt230*tt23*coef(6,&
&1)+tt233+8*tt228*tt17*coef(5,1)+tt232+tt231+tt229
END 
SUBROUTINE polynomial_elas_invar(val, X, Dm, coef) 
IMPLICIT NONE 
REAL(KIND=8) val(1, 1) 
REAL(KIND=8) X(2, 3) 
REAL(KIND=8) Dm(2, 2) 
REAL(KIND=8) coef(3, 1) 
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
tt5 = tt4**2
tt6 = tt3*Dm(2,2)+Dm(1,2)*tt2
tt7 = tt6**2
tt8 = -X(2,1)
tt9 = X(2,2)+tt8
tt10 = X(2,3)+tt8
tt11 = Dm(2,1)*tt10+Dm(1,1)*tt9
tt12 = tt11**2
tt13 = Dm(2,2)*tt10+Dm(1,2)*tt9
tt14 = tt13**2
tt15 = sqrt((tt14+tt12+tt7+tt5)**2-4*((tt12+tt5)*(tt14+tt7)-(tt11&
&*tt13+tt4*tt6)**2)+1.0E-10)
tt16 = ((-tt15)+tt14+tt12+tt7+tt5)/2.0E+0
tt17 = (tt15+tt14+tt12+tt7+tt5)/2.0E+0
val(1,1) = (tt17+tt16-2)**2*coef(3,1)+coef(1,1)*(tt17-1)**2+coef(&
&2,1)*(tt16-1)**2
END 
SUBROUTINE polynomial_elas_invar_jac(jac, X, Dm, coef) 
IMPLICIT NONE 
REAL(KIND=8) jac(1, 6) 
REAL(KIND=8) X(2, 3) 
REAL(KIND=8) Dm(2, 2) 
REAL(KIND=8) coef(3, 1) 
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
tt1 = (-Dm(2,1))-Dm(1,1)
tt2 = -X(1,1)
tt3 = X(1,2)+tt2
tt4 = X(1,3)+tt2
tt5 = tt4*Dm(2,1)+Dm(1,1)*tt3
tt6 = 2*tt1*tt5
tt7 = (-Dm(2,2))-Dm(1,2)
tt8 = tt4*Dm(2,2)+Dm(1,2)*tt3
tt9 = 2*tt7*tt8
tt10 = tt5**2
tt11 = tt8**2
tt12 = -X(2,1)
tt13 = X(2,2)+tt12
tt14 = X(2,3)+tt12
tt15 = Dm(2,1)*tt14+Dm(1,1)*tt13
tt16 = tt15**2
tt17 = Dm(2,2)*tt14+Dm(1,2)*tt13
tt18 = tt17**2
tt19 = tt18+tt16+tt11+tt10
tt20 = tt15*tt17+tt5*tt8
tt21 = tt16+tt10
tt22 = tt18+tt11
tt23 = 2*(tt9+tt6)*tt19-4*(2*tt1*tt5*tt22+2*tt7*tt8*tt21-2*(tt1*t&
&t8+tt5*tt7)*tt20)
tt24 = sqrt(tt19**2-4*(tt21*tt22-tt20**2)+1.0E-10)
tt25 = 1/tt24
tt26 = (-(tt23*tt25)/2.0E+0)+tt9+tt6
tt27 = ((-tt24)+tt18+tt16+tt11+tt10)/2.0E+0
tt28 = tt27-1
tt29 = (tt23*tt25)/2.0E+0+tt9+tt6
tt30 = (tt24+tt18+tt16+tt11+tt10)/2.0E+0
tt31 = tt30-1
tt32 = tt30+tt27-2
tt33 = 2*tt1*tt15
tt34 = 2*tt7*tt17
tt35 = 2*(tt34+tt33)*tt19-4*(2*tt1*tt15*tt22+2*tt7*tt17*tt21-2*(t&
&t1*tt17+tt7*tt15)*tt20)
tt36 = (-(tt35*tt25)/2.0E+0)+tt34+tt33
tt37 = (tt35*tt25)/2.0E+0+tt34+tt33
tt38 = 2*Dm(1,1)*tt5
tt39 = 2*Dm(1,2)*tt8
tt40 = 2*(tt39+tt38)*tt19-4*(2*Dm(1,1)*tt5*tt22+2*Dm(1,2)*tt8*tt2&
&1-2*(Dm(1,1)*tt8+Dm(1,2)*tt5)*tt20)
tt41 = (-(tt40*tt25)/2.0E+0)+tt39+tt38
tt42 = (tt40*tt25)/2.0E+0+tt39+tt38
tt43 = 2*Dm(1,1)*tt15
tt44 = 2*Dm(1,2)*tt17
tt45 = 2*(tt44+tt43)*tt19-4*(2*Dm(1,1)*tt15*tt22+2*Dm(1,2)*tt17*t&
&t21-2*(Dm(1,1)*tt17+Dm(1,2)*tt15)*tt20)
tt46 = (-(tt45*tt25)/2.0E+0)+tt44+tt43
tt47 = (tt45*tt25)/2.0E+0+tt44+tt43
tt48 = 2*Dm(2,1)*tt5
tt49 = 2*Dm(2,2)*tt8
tt50 = 2*(tt49+tt48)*tt19-4*(2*Dm(2,1)*tt5*tt22+2*Dm(2,2)*tt8*tt2&
&1-2*(Dm(2,1)*tt8+tt5*Dm(2,2))*tt20)
tt51 = (-(tt50*tt25)/2.0E+0)+tt49+tt48
tt52 = (tt50*tt25)/2.0E+0+tt49+tt48
tt53 = 2*Dm(2,1)*tt15
tt54 = 2*Dm(2,2)*tt17
tt55 = 2*(tt54+tt53)*tt19-4*(2*Dm(2,1)*tt15*tt22+2*Dm(2,2)*tt17*t&
&t21-2*(Dm(2,1)*tt17+Dm(2,2)*tt15)*tt20)
tt56 = (-(tt55*tt25)/2.0E+0)+tt54+tt53
tt57 = (tt55*tt25)/2.0E+0+tt54+tt53
jac(1,1) = 2*(tt29/2.0E+0+tt26/2.0E+0)*tt32*coef(3,1)+coef(1,1)*t&
&t29*tt31+coef(2,1)*tt26*tt28
jac(1,2) = 2*(tt37/2.0E+0+tt36/2.0E+0)*tt32*coef(3,1)+coef(1,1)*t&
&t37*tt31+coef(2,1)*tt36*tt28
jac(1,3) = 2*(tt42/2.0E+0+tt41/2.0E+0)*tt32*coef(3,1)+coef(1,1)*t&
&t42*tt31+coef(2,1)*tt41*tt28
jac(1,4) = 2*(tt47/2.0E+0+tt46/2.0E+0)*tt32*coef(3,1)+coef(1,1)*t&
&t47*tt31+coef(2,1)*tt46*tt28
jac(1,5) = 2*(tt52/2.0E+0+tt51/2.0E+0)*tt32*coef(3,1)+coef(1,1)*t&
&t52*tt31+coef(2,1)*tt51*tt28
jac(1,6) = 2*(tt57/2.0E+0+tt56/2.0E+0)*tt32*coef(3,1)+coef(1,1)*t&
&t57*tt31+coef(2,1)*tt56*tt28
END 
SUBROUTINE polynomial_elas_invar_hes(hes, X, Dm, coef) 
IMPLICIT NONE 
REAL(KIND=8) hes(6, 6) 
REAL(KIND=8) X(2, 3) 
REAL(KIND=8) Dm(2, 2) 
REAL(KIND=8) coef(3, 1) 
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
tt1 = (-Dm(2,1))-Dm(1,1)
tt2 = tt1**2
tt3 = 2*tt2
tt4 = (-Dm(2,2))-Dm(1,2)
tt5 = tt4**2
tt6 = 2*tt5
tt7 = -X(1,1)
tt8 = X(1,2)+tt7
tt9 = X(1,3)+tt7
tt10 = tt9*Dm(2,1)+Dm(1,1)*tt8
tt11 = 2*tt1*tt10
tt12 = tt9*Dm(2,2)+Dm(1,2)*tt8
tt13 = 2*tt4*tt12
tt14 = tt13+tt11
tt15 = tt10**2
tt16 = tt12**2
tt17 = -X(2,1)
tt18 = X(2,2)+tt17
tt19 = X(2,3)+tt17
tt20 = Dm(2,1)*tt19+Dm(1,1)*tt18
tt21 = tt20**2
tt22 = Dm(2,2)*tt19+Dm(1,2)*tt18
tt23 = tt22**2
tt24 = tt23+tt21+tt16+tt15
tt25 = tt1*tt12+tt10*tt4
tt26 = tt20*tt22+tt10*tt12
tt27 = tt21+tt15
tt28 = tt23+tt16
tt29 = 2*tt14*tt24-4*(2*tt1*tt10*tt28+2*tt4*tt12*tt27-2*tt25*tt26&
&)
tt30 = tt29**2
tt31 = sqrt(tt24**2-4*(tt27*tt28-tt26**2)+1.0E-10)
tt32 = 1/tt31**3
tt33 = 2*(tt6+tt3)*tt24
tt34 = -4*tt1*tt4*tt26
tt35 = 2*tt5*tt27
tt36 = 2*tt2*tt28
tt37 = (-4*(tt36+tt35+tt34-2*tt25**2+8*tt1*tt10*tt4*tt12))+tt33+2&
&*tt14**2
tt38 = 1/tt31
tt39 = (-(tt37*tt38)/2.0E+0)+(tt30*tt32)/4.0E+0+tt6+tt3
tt40 = ((-tt31)+tt23+tt21+tt16+tt15)/2.0E+0
tt41 = tt40-1
tt42 = (tt37*tt38)/2.0E+0-(tt30*tt32)/4.0E+0+tt6+tt3
tt43 = (tt31+tt23+tt21+tt16+tt15)/2.0E+0
tt44 = tt43-1
tt45 = (-(tt29*tt38)/2.0E+0)+tt13+tt11
tt46 = (tt29*tt38)/2.0E+0+tt13+tt11
tt47 = tt43+tt40-2
tt48 = tt46/2.0E+0+tt45/2.0E+0
tt49 = 2*tt1*tt20
tt50 = 2*tt4*tt22
tt51 = tt50+tt49
tt52 = tt1*tt22+tt4*tt20
tt53 = 2*tt51*tt24-4*(2*tt1*tt20*tt28+2*tt4*tt22*tt27-2*tt52*tt26&
&)
tt54 = (-(tt53*tt38)/2.0E+0)+tt50+tt49
tt55 = (tt53*tt38)/2.0E+0+tt50+tt49
tt56 = 2*tt14*tt51-4*((-2*tt25*tt52)+4*tt1*tt10*tt4*tt22+4*tt1*tt&
&4*tt12*tt20)
tt57 = (tt29*tt53*tt32)/4.0E+0-(tt56*tt38)/2.0E+0
tt58 = (tt56*tt38)/2.0E+0-(tt29*tt53*tt32)/4.0E+0
tt59 = tt55/2.0E+0+tt54/2.0E+0
tt60 = 2*(tt58/2.0E+0+tt57/2.0E+0)*tt47*coef(3,1)+2*tt48*tt59*coe&
&f(3,1)+coef(1,1)*tt58*tt44+coef(2,1)*tt57*tt41+(coef(1,1)*tt46*tt&
&55)/2.0E+0+(coef(2,1)*tt45*tt54)/2.0E+0
tt61 = 2*Dm(1,1)*tt10
tt62 = 2*Dm(1,2)*tt12
tt63 = tt62+tt61
tt64 = Dm(1,1)*tt12+Dm(1,2)*tt10
tt65 = 2*tt63*tt24-4*(2*Dm(1,1)*tt10*tt28+2*Dm(1,2)*tt12*tt27-2*t&
&t64*tt26)
tt66 = (-(tt65*tt38)/2.0E+0)+tt62+tt61
tt67 = (tt65*tt38)/2.0E+0+tt62+tt61
tt68 = 2*Dm(1,1)*tt1
tt69 = 2*Dm(1,2)*tt4
tt70 = 2*(tt69+tt68)*tt24
tt71 = -2*(Dm(1,1)*tt4+Dm(1,2)*tt1)*tt26
tt72 = 2*Dm(1,2)*tt4*tt27
tt73 = 2*Dm(1,1)*tt1*tt28
tt74 = (-4*(tt73+tt72+tt71-2*tt64*tt25+4*Dm(1,1)*tt10*tt4*tt12+4*&
&Dm(1,2)*tt1*tt10*tt12))+tt70+2*tt63*tt14
tt75 = (-(tt74*tt38)/2.0E+0)+(tt65*tt29*tt32)/4.0E+0+tt69+tt68
tt76 = (tt74*tt38)/2.0E+0-(tt65*tt29*tt32)/4.0E+0+tt69+tt68
tt77 = tt67/2.0E+0+tt66/2.0E+0
tt78 = 2*(tt76/2.0E+0+tt75/2.0E+0)*tt47*coef(3,1)+2*tt77*tt48*coe&
&f(3,1)+coef(1,1)*tt76*tt44+coef(2,1)*tt75*tt41+(coef(1,1)*tt67*tt&
&46)/2.0E+0+(coef(2,1)*tt66*tt45)/2.0E+0
tt79 = 2*Dm(1,1)*tt20
tt80 = 2*Dm(1,2)*tt22
tt81 = tt80+tt79
tt82 = Dm(1,1)*tt22+Dm(1,2)*tt20
tt83 = 2*tt81*tt24-4*(2*Dm(1,1)*tt20*tt28+2*Dm(1,2)*tt22*tt27-2*t&
&t82*tt26)
tt84 = (-(tt83*tt38)/2.0E+0)+tt80+tt79
tt85 = (tt83*tt38)/2.0E+0+tt80+tt79
tt86 = 2*tt14*tt81-4*((-2*tt25*tt82)+4*Dm(1,2)*tt1*tt10*tt22+4*Dm&
&(1,1)*tt4*tt12*tt20)
tt87 = (tt29*tt83*tt32)/4.0E+0-(tt86*tt38)/2.0E+0
tt88 = (tt86*tt38)/2.0E+0-(tt29*tt83*tt32)/4.0E+0
tt89 = tt85/2.0E+0+tt84/2.0E+0
tt90 = 2*(tt88/2.0E+0+tt87/2.0E+0)*tt47*coef(3,1)+2*tt48*tt89*coe&
&f(3,1)+coef(1,1)*tt88*tt44+coef(2,1)*tt87*tt41+(coef(1,1)*tt46*tt&
&85)/2.0E+0+(coef(2,1)*tt45*tt84)/2.0E+0
tt91 = 2*Dm(2,1)*tt10
tt92 = 2*Dm(2,2)*tt12
tt93 = tt92+tt91
tt94 = Dm(2,1)*tt12+tt10*Dm(2,2)
tt95 = 2*tt93*tt24-4*(2*Dm(2,1)*tt10*tt28+2*Dm(2,2)*tt12*tt27-2*t&
&t94*tt26)
tt96 = (-(tt95*tt38)/2.0E+0)+tt92+tt91
tt97 = (tt95*tt38)/2.0E+0+tt92+tt91
tt98 = 2*tt1*Dm(2,1)
tt99 = 2*tt4*Dm(2,2)
tt100 = 2*(tt99+tt98)*tt24
tt101 = -2*(tt1*Dm(2,2)+Dm(2,1)*tt4)*tt26
tt102 = 2*tt4*Dm(2,2)*tt27
tt103 = 2*tt1*Dm(2,1)*tt28
tt104 = (-4*(tt103+tt102+tt101-2*tt25*tt94+4*tt1*tt10*Dm(2,2)*tt1&
&2+4*Dm(2,1)*tt10*tt4*tt12))+tt100+2*tt14*tt93
tt105 = (-(tt104*tt38)/2.0E+0)+(tt29*tt95*tt32)/4.0E+0+tt99+tt98
tt106 = (tt104*tt38)/2.0E+0-(tt29*tt95*tt32)/4.0E+0+tt99+tt98
tt107 = tt97/2.0E+0+tt96/2.0E+0
tt108 = 2*(tt106/2.0E+0+tt105/2.0E+0)*tt47*coef(3,1)+2*tt48*tt107&
&*coef(3,1)+coef(1,1)*tt106*tt44+coef(2,1)*tt105*tt41+(coef(1,1)*t&
&t46*tt97)/2.0E+0+(coef(2,1)*tt45*tt96)/2.0E+0
tt109 = 2*Dm(2,1)*tt20
tt110 = 2*Dm(2,2)*tt22
tt111 = tt110+tt109
tt112 = Dm(2,1)*tt22+Dm(2,2)*tt20
tt113 = 2*tt111*tt24-4*(2*Dm(2,1)*tt20*tt28+2*Dm(2,2)*tt22*tt27-2&
&*tt112*tt26)
tt114 = (-(tt113*tt38)/2.0E+0)+tt110+tt109
tt115 = (tt113*tt38)/2.0E+0+tt110+tt109
tt116 = 2*tt14*tt111-4*((-2*tt25*tt112)+4*tt1*tt10*Dm(2,2)*tt22+4&
&*Dm(2,1)*tt4*tt12*tt20)
tt117 = (tt29*tt113*tt32)/4.0E+0-(tt116*tt38)/2.0E+0
tt118 = (tt116*tt38)/2.0E+0-(tt29*tt113*tt32)/4.0E+0
tt119 = tt115/2.0E+0+tt114/2.0E+0
tt120 = 2*(tt118/2.0E+0+tt117/2.0E+0)*tt47*coef(3,1)+2*tt48*tt119&
&*coef(3,1)+coef(1,1)*tt118*tt44+coef(2,1)*tt117*tt41+(coef(1,1)*t&
&t46*tt115)/2.0E+0+(coef(2,1)*tt45*tt114)/2.0E+0
tt121 = tt53**2
tt122 = 2*tt51**2-4*((-2*tt52**2)+tt36+tt35+tt34+8*tt1*tt4*tt20*t&
&t22)+tt33
tt123 = (-(tt122*tt38)/2.0E+0)+(tt121*tt32)/4.0E+0+tt6+tt3
tt124 = (tt122*tt38)/2.0E+0-(tt121*tt32)/4.0E+0+tt6+tt3
tt125 = 2*tt63*tt51-4*((-2*tt64*tt52)+4*Dm(1,1)*tt10*tt4*tt22+4*D&
&m(1,2)*tt1*tt12*tt20)
tt126 = (tt65*tt53*tt32)/4.0E+0-(tt125*tt38)/2.0E+0
tt127 = (tt125*tt38)/2.0E+0-(tt65*tt53*tt32)/4.0E+0
tt128 = 2*(tt127/2.0E+0+tt126/2.0E+0)*tt47*coef(3,1)+2*tt77*tt59*&
&coef(3,1)+coef(1,1)*tt127*tt44+coef(2,1)*tt126*tt41+(coef(1,1)*tt&
&67*tt55)/2.0E+0+(coef(2,1)*tt66*tt54)/2.0E+0
tt129 = (-4*(tt73+tt72+tt71-2*tt82*tt52+4*Dm(1,1)*tt4*tt20*tt22+4&
&*Dm(1,2)*tt1*tt20*tt22))+tt70+2*tt81*tt51
tt130 = (-(tt129*tt38)/2.0E+0)+(tt83*tt53*tt32)/4.0E+0+tt69+tt68
tt131 = (tt129*tt38)/2.0E+0-(tt83*tt53*tt32)/4.0E+0+tt69+tt68
tt132 = 2*(tt131/2.0E+0+tt130/2.0E+0)*tt47*coef(3,1)+2*tt89*tt59*&
&coef(3,1)+coef(1,1)*tt131*tt44+coef(2,1)*tt130*tt41+(coef(1,1)*tt&
&85*tt55)/2.0E+0+(coef(2,1)*tt84*tt54)/2.0E+0
tt133 = 2*tt93*tt51-4*((-2*tt94*tt52)+4*Dm(2,1)*tt10*tt4*tt22+4*t&
&t1*Dm(2,2)*tt12*tt20)
tt134 = (tt95*tt53*tt32)/4.0E+0-(tt133*tt38)/2.0E+0
tt135 = (tt133*tt38)/2.0E+0-(tt95*tt53*tt32)/4.0E+0
tt136 = 2*(tt135/2.0E+0+tt134/2.0E+0)*tt47*coef(3,1)+2*tt107*tt59&
&*coef(3,1)+coef(1,1)*tt135*tt44+coef(2,1)*tt134*tt41+(coef(1,1)*t&
&t97*tt55)/2.0E+0+(coef(2,1)*tt96*tt54)/2.0E+0
tt137 = (-4*(tt103+tt102+tt101-2*tt52*tt112+4*tt1*Dm(2,2)*tt20*tt&
&22+4*Dm(2,1)*tt4*tt20*tt22))+tt100+2*tt51*tt111
tt138 = (-(tt137*tt38)/2.0E+0)+(tt53*tt113*tt32)/4.0E+0+tt99+tt98&
&
tt139 = (tt137*tt38)/2.0E+0-(tt53*tt113*tt32)/4.0E+0+tt99+tt98
tt140 = 2*(tt139/2.0E+0+tt138/2.0E+0)*tt47*coef(3,1)+2*tt59*tt119&
&*coef(3,1)+coef(1,1)*tt139*tt44+coef(2,1)*tt138*tt41+(coef(1,1)*t&
&t55*tt115)/2.0E+0+(coef(2,1)*tt54*tt114)/2.0E+0
tt141 = Dm(1,1)**2
tt142 = 2*tt141
tt143 = Dm(1,2)**2
tt144 = 2*tt143
tt145 = tt65**2
tt146 = 2*(tt144+tt142)*tt24
tt147 = -4*Dm(1,1)*Dm(1,2)*tt26
tt148 = 2*tt143*tt27
tt149 = 2*tt141*tt28
tt150 = (-4*(tt149+tt148+tt147-2*tt64**2+8*Dm(1,1)*Dm(1,2)*tt10*t&
&t12))+tt146+2*tt63**2
tt151 = (-(tt150*tt38)/2.0E+0)+(tt145*tt32)/4.0E+0+tt144+tt142
tt152 = (tt150*tt38)/2.0E+0-(tt145*tt32)/4.0E+0+tt144+tt142
tt153 = 2*tt63*tt81-4*((-2*tt64*tt82)+4*Dm(1,1)*Dm(1,2)*tt10*tt22&
&+4*Dm(1,1)*Dm(1,2)*tt12*tt20)
tt154 = (tt65*tt83*tt32)/4.0E+0-(tt153*tt38)/2.0E+0
tt155 = (tt153*tt38)/2.0E+0-(tt65*tt83*tt32)/4.0E+0
tt156 = 2*(tt155/2.0E+0+tt154/2.0E+0)*tt47*coef(3,1)+2*tt77*tt89*&
&coef(3,1)+coef(1,1)*tt155*tt44+coef(2,1)*tt154*tt41+(coef(1,1)*tt&
&67*tt85)/2.0E+0+(coef(2,1)*tt66*tt84)/2.0E+0
tt157 = 2*Dm(1,1)*Dm(2,1)
tt158 = 2*Dm(1,2)*Dm(2,2)
tt159 = 2*(tt158+tt157)*tt24
tt160 = -2*(Dm(1,1)*Dm(2,2)+Dm(1,2)*Dm(2,1))*tt26
tt161 = 2*Dm(1,2)*Dm(2,2)*tt27
tt162 = 2*Dm(1,1)*Dm(2,1)*tt28
tt163 = (-4*(tt162+tt161+tt160-2*tt64*tt94+4*Dm(1,1)*tt10*Dm(2,2)&
&*tt12+4*Dm(1,2)*Dm(2,1)*tt10*tt12))+tt159+2*tt63*tt93
tt164 = (-(tt163*tt38)/2.0E+0)+(tt65*tt95*tt32)/4.0E+0+tt158+tt15&
&7
tt165 = (tt163*tt38)/2.0E+0-(tt65*tt95*tt32)/4.0E+0+tt158+tt157
tt166 = 2*(tt165/2.0E+0+tt164/2.0E+0)*tt47*coef(3,1)+2*tt77*tt107&
&*coef(3,1)+coef(1,1)*tt165*tt44+coef(2,1)*tt164*tt41+(coef(1,1)*t&
&t67*tt97)/2.0E+0+(coef(2,1)*tt66*tt96)/2.0E+0
tt167 = 2*tt63*tt111-4*((-2*tt64*tt112)+4*Dm(1,1)*tt10*Dm(2,2)*tt&
&22+4*Dm(1,2)*Dm(2,1)*tt12*tt20)
tt168 = (tt65*tt113*tt32)/4.0E+0-(tt167*tt38)/2.0E+0
tt169 = (tt167*tt38)/2.0E+0-(tt65*tt113*tt32)/4.0E+0
tt170 = 2*(tt169/2.0E+0+tt168/2.0E+0)*tt47*coef(3,1)+2*tt77*tt119&
&*coef(3,1)+coef(1,1)*tt169*tt44+coef(2,1)*tt168*tt41+(coef(1,1)*t&
&t67*tt115)/2.0E+0+(coef(2,1)*tt66*tt114)/2.0E+0
tt171 = tt83**2
tt172 = 2*tt81**2-4*((-2*tt82**2)+tt149+tt148+tt147+8*Dm(1,1)*Dm(&
&1,2)*tt20*tt22)+tt146
tt173 = (-(tt172*tt38)/2.0E+0)+(tt171*tt32)/4.0E+0+tt144+tt142
tt174 = (tt172*tt38)/2.0E+0-(tt171*tt32)/4.0E+0+tt144+tt142
tt175 = 2*tt93*tt81-4*((-2*tt94*tt82)+4*Dm(1,2)*Dm(2,1)*tt10*tt22&
&+4*Dm(1,1)*Dm(2,2)*tt12*tt20)
tt176 = (tt95*tt83*tt32)/4.0E+0-(tt175*tt38)/2.0E+0
tt177 = (tt175*tt38)/2.0E+0-(tt95*tt83*tt32)/4.0E+0
tt178 = 2*(tt177/2.0E+0+tt176/2.0E+0)*tt47*coef(3,1)+2*tt107*tt89&
&*coef(3,1)+coef(1,1)*tt177*tt44+coef(2,1)*tt176*tt41+(coef(1,1)*t&
&t97*tt85)/2.0E+0+(coef(2,1)*tt96*tt84)/2.0E+0
tt179 = (-4*(tt162+tt161+tt160-2*tt82*tt112+4*Dm(1,1)*Dm(2,2)*tt2&
&0*tt22+4*Dm(1,2)*Dm(2,1)*tt20*tt22))+tt159+2*tt81*tt111
tt180 = (-(tt179*tt38)/2.0E+0)+(tt83*tt113*tt32)/4.0E+0+tt158+tt1&
&57
tt181 = (tt179*tt38)/2.0E+0-(tt83*tt113*tt32)/4.0E+0+tt158+tt157
tt182 = 2*(tt181/2.0E+0+tt180/2.0E+0)*tt47*coef(3,1)+2*tt89*tt119&
&*coef(3,1)+coef(1,1)*tt181*tt44+coef(2,1)*tt180*tt41+(coef(1,1)*t&
&t85*tt115)/2.0E+0+(coef(2,1)*tt84*tt114)/2.0E+0
tt183 = Dm(2,1)**2
tt184 = 2*tt183
tt185 = Dm(2,2)**2
tt186 = 2*tt185
tt187 = tt95**2
tt188 = 2*(tt186+tt184)*tt24
tt189 = -4*Dm(2,1)*Dm(2,2)*tt26
tt190 = 2*tt185*tt27
tt191 = 2*tt183*tt28
tt192 = (-4*(tt191+tt190+tt189-2*tt94**2+8*Dm(2,1)*tt10*Dm(2,2)*t&
&t12))+tt188+2*tt93**2
tt193 = (-(tt192*tt38)/2.0E+0)+(tt187*tt32)/4.0E+0+tt186+tt184
tt194 = (tt192*tt38)/2.0E+0-(tt187*tt32)/4.0E+0+tt186+tt184
tt195 = 2*tt93*tt111-4*((-2*tt94*tt112)+4*Dm(2,1)*tt10*Dm(2,2)*tt&
&22+4*Dm(2,1)*Dm(2,2)*tt12*tt20)
tt196 = (tt95*tt113*tt32)/4.0E+0-(tt195*tt38)/2.0E+0
tt197 = (tt195*tt38)/2.0E+0-(tt95*tt113*tt32)/4.0E+0
tt198 = 2*(tt197/2.0E+0+tt196/2.0E+0)*tt47*coef(3,1)+2*tt107*tt11&
&9*coef(3,1)+coef(1,1)*tt197*tt44+coef(2,1)*tt196*tt41+(coef(1,1)*&
&tt97*tt115)/2.0E+0+(coef(2,1)*tt96*tt114)/2.0E+0
tt199 = tt113**2
tt200 = 2*tt111**2-4*((-2*tt112**2)+tt191+tt190+tt189+8*Dm(2,1)*D&
&m(2,2)*tt20*tt22)+tt188
tt201 = (-(tt200*tt38)/2.0E+0)+(tt199*tt32)/4.0E+0+tt186+tt184
tt202 = (tt200*tt38)/2.0E+0-(tt199*tt32)/4.0E+0+tt186+tt184
hes(1,1) = 2*tt48**2*coef(3,1)+2*(tt42/2.0E+0+tt39/2.0E+0)*tt47*c&
&oef(3,1)+(coef(1,1)*tt46**2)/2.0E+0+(coef(2,1)*tt45**2)/2.0E+0+co&
&ef(1,1)*tt42*tt44+coef(2,1)*tt39*tt41
hes(1,2) = tt60
hes(1,3) = tt78
hes(1,4) = tt90
hes(1,5) = tt108
hes(1,6) = tt120
hes(2,1) = tt60
hes(2,2) = 2*tt59**2*coef(3,1)+2*(tt124/2.0E+0+tt123/2.0E+0)*tt47&
&*coef(3,1)+(coef(1,1)*tt55**2)/2.0E+0+(coef(2,1)*tt54**2)/2.0E+0+&
&coef(1,1)*tt124*tt44+coef(2,1)*tt123*tt41
hes(2,3) = tt128
hes(2,4) = tt132
hes(2,5) = tt136
hes(2,6) = tt140
hes(3,1) = tt78
hes(3,2) = tt128
hes(3,3) = 2*tt77**2*coef(3,1)+2*(tt152/2.0E+0+tt151/2.0E+0)*tt47&
&*coef(3,1)+(coef(1,1)*tt67**2)/2.0E+0+(coef(2,1)*tt66**2)/2.0E+0+&
&coef(1,1)*tt152*tt44+coef(2,1)*tt151*tt41
hes(3,4) = tt156
hes(3,5) = tt166
hes(3,6) = tt170
hes(4,1) = tt90
hes(4,2) = tt132
hes(4,3) = tt156
hes(4,4) = 2*tt89**2*coef(3,1)+2*(tt174/2.0E+0+tt173/2.0E+0)*tt47&
&*coef(3,1)+(coef(1,1)*tt85**2)/2.0E+0+(coef(2,1)*tt84**2)/2.0E+0+&
&coef(1,1)*tt174*tt44+coef(2,1)*tt173*tt41
hes(4,5) = tt178
hes(4,6) = tt182
hes(5,1) = tt108
hes(5,2) = tt136
hes(5,3) = tt166
hes(5,4) = tt178
hes(5,5) = 2*tt107**2*coef(3,1)+2*(tt194/2.0E+0+tt193/2.0E+0)*tt4&
&7*coef(3,1)+(coef(1,1)*tt97**2)/2.0E+0+(coef(2,1)*tt96**2)/2.0E+0&
&+coef(1,1)*tt194*tt44+coef(2,1)*tt193*tt41
hes(5,6) = tt198
hes(6,1) = tt120
hes(6,2) = tt140
hes(6,3) = tt170
hes(6,4) = tt182
hes(6,5) = tt198
hes(6,6) = 2*tt119**2*coef(3,1)+2*(tt202/2.0E+0+tt201/2.0E+0)*tt4&
&7*coef(3,1)+(coef(1,1)*tt115**2)/2.0E+0+(coef(2,1)*tt114**2)/2.0E&
&+0+coef(1,1)*tt202*tt44+coef(2,1)*tt201*tt41
END 
SUBROUTINE polynomial_elas_invar1(val, X, Dm, coef, area) 
IMPLICIT NONE 
REAL(KIND=8) val(1, 1) 
REAL(KIND=8) X(2, 3) 
REAL(KIND=8) Dm(2, 2) 
REAL(KIND=8) coef(5, 1) 
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
tt1 = X(1,1)**2
tt2 = X(2,1)**2
tt3 = -X(1,1)
tt4 = X(1,2)+tt3
tt5 = X(1,3)+tt3
tt6 = tt5*Dm(2,1)+Dm(1,1)*tt4
tt7 = tt6**2
tt8 = -X(2,1)
tt9 = X(2,2)+tt8
tt10 = X(2,3)+tt8
tt11 = Dm(2,1)*tt10+Dm(1,1)*tt9
tt12 = tt11**2
tt13 = tt5*Dm(2,2)+Dm(1,2)*tt4
tt14 = tt13**2
tt15 = Dm(2,2)*tt10+Dm(1,2)*tt9
tt16 = tt15**2
tt17 = -tt13*tt11
tt18 = tt6*tt15
tt19 = tt18+tt17
val(1,1) = area(1,1)*(coef(5,1)*log(tt19)**2+(tt18+tt17-1)**2*coe&
&f(4,1)+((tt16+tt12+tt14+tt7)/tt19-2)*coef(3,1)+coef(2,1)*((tt16+t&
&t14-1)/2.0E+0+(tt12+tt7-1)/2.0E+0)**2+coef(1,1)*(((Dm(2,2)*X(2,3)&
&+Dm(1,2)*X(2,2)-X(2,1)*Dm(2,2)-Dm(1,2)*X(2,1))**2+(X(1,3)*Dm(2,2)&
&-X(1,1)*Dm(2,2)+Dm(1,2)*X(1,2)-X(1,1)*Dm(1,2))**2-1)**2/4.0E+0+((&
&Dm(2,1)*X(2,3)+Dm(1,1)*X(2,2)-Dm(2,1)*X(2,1)-Dm(1,1)*X(2,1))**2+(&
&X(1,3)*Dm(2,1)-X(1,1)*Dm(2,1)+Dm(1,1)*X(1,2)-Dm(1,1)*X(1,1))**2-1&
&)**2/4.0E+0+(Dm(2,1)*Dm(2,2)*X(2,3)**2+Dm(1,1)*Dm(2,2)*X(2,2)*X(2&
&,3)+Dm(1,2)*Dm(2,1)*X(2,2)*X(2,3)-2*Dm(2,1)*X(2,1)*Dm(2,2)*X(2,3)&
&-Dm(1,1)*X(2,1)*Dm(2,2)*X(2,3)-Dm(1,2)*Dm(2,1)*X(2,1)*X(2,3)+Dm(1&
&,1)*Dm(1,2)*X(2,2)**2-Dm(1,1)*X(2,1)*Dm(2,2)*X(2,2)-Dm(1,2)*Dm(2,&
&1)*X(2,1)*X(2,2)-2*Dm(1,1)*Dm(1,2)*X(2,1)*X(2,2)+Dm(2,1)*tt2*Dm(2&
&,2)+Dm(1,1)*tt2*Dm(2,2)+X(1,3)**2*Dm(2,1)*Dm(2,2)-2*X(1,1)*X(1,3)&
&*Dm(2,1)*Dm(2,2)+tt1*Dm(2,1)*Dm(2,2)+Dm(1,1)*X(1,2)*X(1,3)*Dm(2,2&
&)-Dm(1,1)*X(1,1)*X(1,3)*Dm(2,2)-Dm(1,1)*X(1,1)*X(1,2)*Dm(2,2)+Dm(&
&1,1)*tt1*Dm(2,2)+Dm(1,2)*Dm(2,1)*tt2+Dm(1,1)*Dm(1,2)*tt2+Dm(1,2)*&
&X(1,2)*X(1,3)*Dm(2,1)-X(1,1)*Dm(1,2)*X(1,3)*Dm(2,1)-X(1,1)*Dm(1,2&
&)*X(1,2)*Dm(2,1)+tt1*Dm(1,2)*Dm(2,1)+Dm(1,1)*Dm(1,2)*X(1,2)**2-2*&
&Dm(1,1)*X(1,1)*Dm(1,2)*X(1,2)+Dm(1,1)*tt1*Dm(1,2))**2/2.0E+0))
END 
SUBROUTINE polynomial_elas_invar1_jac(jac, X, Dm, coef, area) 
IMPLICIT NONE 
REAL(KIND=8) jac(1, 6) 
REAL(KIND=8) X(2, 3) 
REAL(KIND=8) Dm(2, 2) 
REAL(KIND=8) coef(5, 1) 
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
tt1 = X(1,1)**2
tt2 = X(2,1)**2
tt3 = Dm(2,1)*Dm(2,2)*X(2,3)**2+Dm(1,1)*Dm(2,2)*X(2,2)*X(2,3)+Dm(&
&1,2)*Dm(2,1)*X(2,2)*X(2,3)-2*Dm(2,1)*X(2,1)*Dm(2,2)*X(2,3)-Dm(1,1&
&)*X(2,1)*Dm(2,2)*X(2,3)-Dm(1,2)*Dm(2,1)*X(2,1)*X(2,3)+Dm(1,1)*Dm(&
&1,2)*X(2,2)**2-Dm(1,1)*X(2,1)*Dm(2,2)*X(2,2)-Dm(1,2)*Dm(2,1)*X(2,&
&1)*X(2,2)-2*Dm(1,1)*Dm(1,2)*X(2,1)*X(2,2)+Dm(2,1)*tt2*Dm(2,2)+Dm(&
&1,1)*tt2*Dm(2,2)+X(1,3)**2*Dm(2,1)*Dm(2,2)-2*X(1,1)*X(1,3)*Dm(2,1&
&)*Dm(2,2)+tt1*Dm(2,1)*Dm(2,2)+Dm(1,1)*X(1,2)*X(1,3)*Dm(2,2)-Dm(1,&
&1)*X(1,1)*X(1,3)*Dm(2,2)-Dm(1,1)*X(1,1)*X(1,2)*Dm(2,2)+Dm(1,1)*tt&
&1*Dm(2,2)+Dm(1,2)*Dm(2,1)*tt2+Dm(1,1)*Dm(1,2)*tt2+Dm(1,2)*X(1,2)*&
&X(1,3)*Dm(2,1)-X(1,1)*Dm(1,2)*X(1,3)*Dm(2,1)-X(1,1)*Dm(1,2)*X(1,2&
&)*Dm(2,1)+tt1*Dm(1,2)*Dm(2,1)+Dm(1,1)*Dm(1,2)*X(1,2)**2-2*Dm(1,1)&
&*X(1,1)*Dm(1,2)*X(1,2)+Dm(1,1)*tt1*Dm(1,2)
tt4 = (-Dm(2,1))-Dm(1,1)
tt5 = X(1,3)*Dm(2,1)-X(1,1)*Dm(2,1)+Dm(1,1)*X(1,2)-Dm(1,1)*X(1,1)&
&
tt6 = Dm(2,1)*X(2,3)+Dm(1,1)*X(2,2)-Dm(2,1)*X(2,1)-Dm(1,1)*X(2,1)&
&
tt7 = tt6**2+tt5**2-1
tt8 = (-Dm(2,2))-Dm(1,2)
tt9 = X(1,3)*Dm(2,2)-X(1,1)*Dm(2,2)+Dm(1,2)*X(1,2)-X(1,1)*Dm(1,2)&
&
tt10 = Dm(2,2)*X(2,3)+Dm(1,2)*X(2,2)-X(2,1)*Dm(2,2)-Dm(1,2)*X(2,1&
&)
tt11 = tt10**2+tt9**2-1
tt12 = -X(1,1)
tt13 = X(1,2)+tt12
tt14 = X(1,3)+tt12
tt15 = tt14*Dm(2,1)+Dm(1,1)*tt13
tt16 = tt14*Dm(2,2)+Dm(1,2)*tt13
tt17 = tt15**2
tt18 = -X(2,1)
tt19 = X(2,2)+tt18
tt20 = X(2,3)+tt18
tt21 = Dm(2,1)*tt20+Dm(1,1)*tt19
tt22 = tt21**2
tt23 = tt16**2
tt24 = Dm(2,2)*tt20+Dm(1,2)*tt19
tt25 = tt24**2
tt26 = (tt25+tt23-1)/2.0E+0+(tt22+tt17-1)/2.0E+0
tt27 = tt4*tt24-tt8*tt21
tt28 = -tt16*tt21
tt29 = tt15*tt24
tt30 = tt29+tt28
tt31 = 1/tt30**2
tt32 = tt25+tt22+tt23+tt17
tt33 = 1/tt30
tt34 = tt29+tt28-1
tt35 = log(tt30)
tt36 = tt15*tt8-tt4*tt16
tt37 = -X(1,1)*Dm(1,2)*Dm(2,1)
tt38 = -Dm(1,1)*X(1,1)*Dm(2,2)
tt39 = Dm(1,1)*tt24-Dm(1,2)*tt21
tt40 = -Dm(1,2)*Dm(2,1)*X(2,1)
tt41 = -Dm(1,1)*X(2,1)*Dm(2,2)
tt42 = Dm(1,2)*tt15-Dm(1,1)*tt16
tt43 = Dm(2,1)*tt24-Dm(2,2)*tt21
tt44 = tt15*Dm(2,2)-Dm(2,1)*tt16
jac(1,1) = area(1,1)*(2*tt27*tt33*coef(5,1)*tt35+2*tt27*tt34*coef&
&(4,1)+((2*tt8*tt16+2*tt4*tt15)*tt33-tt27*tt31*tt32)*coef(3,1)+2*c&
&oef(2,1)*(tt8*tt16+tt4*tt15)*tt26+coef(1,1)*(tt8*tt9*tt11+tt4*tt5&
&*tt7+((-2*X(1,3)*Dm(2,1)*Dm(2,2))+2*X(1,1)*Dm(2,1)*Dm(2,2)-Dm(1,1&
&)*X(1,3)*Dm(2,2)-Dm(1,1)*X(1,2)*Dm(2,2)+2*Dm(1,1)*X(1,1)*Dm(2,2)-&
&Dm(1,2)*X(1,3)*Dm(2,1)-Dm(1,2)*X(1,2)*Dm(2,1)+2*X(1,1)*Dm(1,2)*Dm&
&(2,1)-2*Dm(1,1)*Dm(1,2)*X(1,2)+2*Dm(1,1)*X(1,1)*Dm(1,2))*tt3))
jac(1,2) = area(1,1)*(2*tt36*tt33*coef(5,1)*tt35+2*tt36*tt34*coef&
&(4,1)+(tt33*(2*tt8*tt24+2*tt4*tt21)-tt36*tt31*tt32)*coef(3,1)+2*c&
&oef(2,1)*(tt8*tt24+tt4*tt21)*tt26+coef(1,1)*(tt8*tt10*tt11+tt4*tt&
&6*tt7+((-2*Dm(2,1)*Dm(2,2)*X(2,3))-Dm(1,1)*Dm(2,2)*X(2,3)-Dm(1,2)&
&*Dm(2,1)*X(2,3)-Dm(1,1)*Dm(2,2)*X(2,2)-Dm(1,2)*Dm(2,1)*X(2,2)-2*D&
&m(1,1)*Dm(1,2)*X(2,2)+2*Dm(2,1)*X(2,1)*Dm(2,2)+2*Dm(1,1)*X(2,1)*D&
&m(2,2)+2*Dm(1,2)*Dm(2,1)*X(2,1)+2*Dm(1,1)*Dm(1,2)*X(2,1))*tt3))
jac(1,3) = area(1,1)*(2*tt39*tt33*coef(5,1)*tt35+2*tt39*tt34*coef&
&(4,1)+((2*Dm(1,2)*tt16+2*Dm(1,1)*tt15)*tt33-tt39*tt31*tt32)*coef(&
&3,1)+2*coef(2,1)*(Dm(1,2)*tt16+Dm(1,1)*tt15)*tt26+coef(1,1)*(Dm(1&
&,2)*tt9*tt11+Dm(1,1)*tt5*tt7+(Dm(1,1)*X(1,3)*Dm(2,2)+tt38+Dm(1,2)&
&*X(1,3)*Dm(2,1)+tt37+2*Dm(1,1)*Dm(1,2)*X(1,2)-2*Dm(1,1)*X(1,1)*Dm&
&(1,2))*tt3))
jac(1,4) = area(1,1)*(2*tt42*tt33*coef(5,1)*tt35+2*tt42*tt34*coef&
&(4,1)+((2*Dm(1,2)*tt24+2*Dm(1,1)*tt21)*tt33-tt42*tt31*tt32)*coef(&
&3,1)+2*coef(2,1)*(Dm(1,2)*tt24+Dm(1,1)*tt21)*tt26+coef(1,1)*(Dm(1&
&,2)*tt10*tt11+Dm(1,1)*tt6*tt7+(Dm(1,1)*Dm(2,2)*X(2,3)+Dm(1,2)*Dm(&
&2,1)*X(2,3)+2*Dm(1,1)*Dm(1,2)*X(2,2)+tt41+tt40-2*Dm(1,1)*Dm(1,2)*&
&X(2,1))*tt3))
jac(1,5) = area(1,1)*(2*tt43*tt33*coef(5,1)*tt35+2*tt43*tt34*coef&
&(4,1)+((2*Dm(2,2)*tt16+2*Dm(2,1)*tt15)*tt33-tt43*tt31*tt32)*coef(&
&3,1)+2*coef(2,1)*(Dm(2,2)*tt16+Dm(2,1)*tt15)*tt26+coef(1,1)*(Dm(2&
&,2)*tt9*tt11+Dm(2,1)*tt5*tt7+(2*X(1,3)*Dm(2,1)*Dm(2,2)-2*X(1,1)*D&
&m(2,1)*Dm(2,2)+Dm(1,1)*X(1,2)*Dm(2,2)+tt38+Dm(1,2)*X(1,2)*Dm(2,1)&
&+tt37)*tt3))
jac(1,6) = area(1,1)*(2*tt44*tt33*coef(5,1)*tt35+2*tt44*tt34*coef&
&(4,1)+(tt33*(2*Dm(2,2)*tt24+2*Dm(2,1)*tt21)-tt44*tt31*tt32)*coef(&
&3,1)+2*coef(2,1)*(Dm(2,2)*tt24+Dm(2,1)*tt21)*tt26+coef(1,1)*(Dm(2&
&,2)*tt10*tt11+Dm(2,1)*tt6*tt7+(2*Dm(2,1)*Dm(2,2)*X(2,3)+Dm(1,1)*D&
&m(2,2)*X(2,2)+Dm(1,2)*Dm(2,1)*X(2,2)-2*Dm(2,1)*X(2,1)*Dm(2,2)+tt4&
&1+tt40)*tt3))
END 
SUBROUTINE polynomial_elas_invar1_hes(hes, X, Dm, coef, area) 
IMPLICIT NONE 
REAL(KIND=8) hes(6, 6) 
REAL(KIND=8) X(2, 3) 
REAL(KIND=8) Dm(2, 2) 
REAL(KIND=8) coef(5, 1) 
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
tt1 = (-Dm(2,1))-Dm(1,1)
tt2 = -X(1,1)
tt3 = X(1,2)+tt2
tt4 = X(1,3)+tt2
tt5 = tt4*Dm(2,1)+Dm(1,1)*tt3
tt6 = (-Dm(2,2))-Dm(1,2)
tt7 = tt4*Dm(2,2)+Dm(1,2)*tt3
tt8 = tt6*tt7+tt1*tt5
tt9 = tt1**2
tt10 = X(1,3)*Dm(2,1)-X(1,1)*Dm(2,1)+Dm(1,1)*X(1,2)-Dm(1,1)*X(1,1&
&)
tt11 = tt10**2
tt12 = tt6**2
tt13 = X(1,3)*Dm(2,2)-X(1,1)*Dm(2,2)+Dm(1,2)*X(1,2)-X(1,1)*Dm(1,2&
&)
tt14 = tt13**2
tt15 = (-2*X(1,3)*Dm(2,1)*Dm(2,2))+2*X(1,1)*Dm(2,1)*Dm(2,2)-Dm(1,&
&1)*X(1,3)*Dm(2,2)-Dm(1,1)*X(1,2)*Dm(2,2)+2*Dm(1,1)*X(1,1)*Dm(2,2)&
&-Dm(1,2)*X(1,3)*Dm(2,1)-Dm(1,2)*X(1,2)*Dm(2,1)+2*X(1,1)*Dm(1,2)*D&
&m(2,1)-2*Dm(1,1)*Dm(1,2)*X(1,2)+2*Dm(1,1)*X(1,1)*Dm(1,2)
tt16 = X(1,1)**2
tt17 = X(2,1)**2
tt18 = Dm(2,1)*Dm(2,2)*X(2,3)**2+Dm(1,1)*Dm(2,2)*X(2,2)*X(2,3)+Dm&
&(1,2)*Dm(2,1)*X(2,2)*X(2,3)-2*Dm(2,1)*X(2,1)*Dm(2,2)*X(2,3)-Dm(1,&
&1)*X(2,1)*Dm(2,2)*X(2,3)-Dm(1,2)*Dm(2,1)*X(2,1)*X(2,3)+Dm(1,1)*Dm&
&(1,2)*X(2,2)**2-Dm(1,1)*X(2,1)*Dm(2,2)*X(2,2)-Dm(1,2)*Dm(2,1)*X(2&
&,1)*X(2,2)-2*Dm(1,1)*Dm(1,2)*X(2,1)*X(2,2)+Dm(2,1)*tt17*Dm(2,2)+D&
&m(1,1)*tt17*Dm(2,2)+X(1,3)**2*Dm(2,1)*Dm(2,2)-2*X(1,1)*X(1,3)*Dm(&
&2,1)*Dm(2,2)+tt16*Dm(2,1)*Dm(2,2)+Dm(1,1)*X(1,2)*X(1,3)*Dm(2,2)-D&
&m(1,1)*X(1,1)*X(1,3)*Dm(2,2)-Dm(1,1)*X(1,1)*X(1,2)*Dm(2,2)+Dm(1,1&
&)*tt16*Dm(2,2)+Dm(1,2)*Dm(2,1)*tt17+Dm(1,1)*Dm(1,2)*tt17+Dm(1,2)*&
&X(1,2)*X(1,3)*Dm(2,1)-X(1,1)*Dm(1,2)*X(1,3)*Dm(2,1)-X(1,1)*Dm(1,2&
&)*X(1,2)*Dm(2,1)+tt16*Dm(1,2)*Dm(2,1)+Dm(1,1)*Dm(1,2)*X(1,2)**2-2&
&*Dm(1,1)*X(1,1)*Dm(1,2)*X(1,2)+Dm(1,1)*tt16*Dm(1,2)
tt19 = (2*Dm(2,1)*Dm(2,2)+2*Dm(1,1)*Dm(2,2)+2*Dm(1,2)*Dm(2,1)+2*D&
&m(1,1)*Dm(1,2))*tt18
tt20 = Dm(2,1)*X(2,3)+Dm(1,1)*X(2,2)-Dm(2,1)*X(2,1)-Dm(1,1)*X(2,1&
&)
tt21 = tt20**2
tt22 = tt21+tt11-1
tt23 = tt9*tt22
tt24 = Dm(2,2)*X(2,3)+Dm(1,2)*X(2,2)-X(2,1)*Dm(2,2)-Dm(1,2)*X(2,1&
&)
tt25 = tt24**2
tt26 = tt25+tt14-1
tt27 = tt12*tt26
tt28 = tt5**2
tt29 = -X(2,1)
tt30 = X(2,2)+tt29
tt31 = X(2,3)+tt29
tt32 = Dm(2,1)*tt31+Dm(1,1)*tt30
tt33 = tt32**2
tt34 = tt7**2
tt35 = Dm(2,2)*tt31+Dm(1,2)*tt30
tt36 = tt35**2
tt37 = (tt36+tt34-1)/2.0E+0+(tt33+tt28-1)/2.0E+0
tt38 = 2*coef(2,1)*(tt12+tt9)*tt37
tt39 = tt1*tt35-tt6*tt32
tt40 = tt39**2
tt41 = -tt7*tt32
tt42 = tt5*tt35
tt43 = tt42+tt41
tt44 = 1/tt43**3
tt45 = tt36+tt33+tt34+tt28
tt46 = 2*tt6*tt7+2*tt1*tt5
tt47 = 1/tt43**2
tt48 = 1/tt43
tt49 = (2*tt12+2*tt9)*tt48
tt50 = log(tt43)
tt51 = (-2*Dm(2,1)*Dm(2,2)*X(2,3))-Dm(1,1)*Dm(2,2)*X(2,3)-Dm(1,2)&
&*Dm(2,1)*X(2,3)-Dm(1,1)*Dm(2,2)*X(2,2)-Dm(1,2)*Dm(2,1)*X(2,2)-2*D&
&m(1,1)*Dm(1,2)*X(2,2)+2*Dm(2,1)*X(2,1)*Dm(2,2)+2*Dm(1,1)*X(2,1)*D&
&m(2,2)+2*Dm(1,2)*Dm(2,1)*X(2,1)+2*Dm(1,1)*Dm(1,2)*X(2,1)
tt52 = tt6*tt35+tt1*tt32
tt53 = 2*tt6*tt35+2*tt1*tt32
tt54 = tt5*tt6-tt1*tt7
tt55 = area(1,1)*((-2*tt54*tt39*tt47*coef(5,1)*tt50)+2*tt54*tt39*&
&tt47*coef(5,1)+2*tt54*tt39*coef(4,1)+((-tt54*tt46*tt47)+2*tt54*tt&
&39*tt44*tt45-tt39*tt47*tt53)*coef(3,1)+2*coef(2,1)*tt8*tt52+coef(&
&1,1)*(tt15*tt51+2*tt12*tt13*tt24+2*tt9*tt10*tt20))
tt56 = Dm(1,2)*tt7+Dm(1,1)*tt5
tt57 = -X(1,1)*Dm(1,2)*Dm(2,1)
tt58 = -Dm(1,1)*X(1,1)*Dm(2,2)
tt59 = Dm(1,1)*X(1,3)*Dm(2,2)+tt58+Dm(1,2)*X(1,3)*Dm(2,1)+tt57+2*&
&Dm(1,1)*Dm(1,2)*X(1,2)-2*Dm(1,1)*X(1,1)*Dm(1,2)
tt60 = -Dm(1,2)*Dm(2,1)
tt61 = -Dm(1,1)*Dm(2,2)
tt62 = (tt61+tt60-2*Dm(1,1)*Dm(1,2))*tt18
tt63 = Dm(1,1)*tt1*tt22
tt64 = Dm(1,2)*tt6*tt26
tt65 = 2*coef(2,1)*(Dm(1,2)*tt6+Dm(1,1)*tt1)*tt37
tt66 = Dm(1,1)*tt35-Dm(1,2)*tt32
tt67 = 2*Dm(1,2)*tt7+2*Dm(1,1)*tt5
tt68 = (2*Dm(1,2)*tt6+2*Dm(1,1)*tt1)*tt48
tt69 = area(1,1)*((-2*tt66*tt39*tt47*coef(5,1)*tt50)+2*tt66*tt39*&
&tt47*coef(5,1)+2*tt66*tt39*coef(4,1)+(tt68-tt67*tt39*tt47-tt46*tt&
&66*tt47+2*tt66*tt39*tt44*tt45)*coef(3,1)+tt65+coef(1,1)*(tt64+tt6&
&3+tt62+2*Dm(1,2)*tt6*tt14+tt59*tt15+2*Dm(1,1)*tt1*tt11)+2*coef(2,&
&1)*tt56*tt8)
tt70 = 2*Dm(1,1)*tt1*tt10*tt20
tt71 = 2*Dm(1,2)*tt6*tt13*tt24
tt72 = -Dm(1,2)*Dm(2,1)*X(2,1)
tt73 = -Dm(1,1)*X(2,1)*Dm(2,2)
tt74 = Dm(1,1)*Dm(2,2)*X(2,3)+Dm(1,2)*Dm(2,1)*X(2,3)+2*Dm(1,1)*Dm&
&(1,2)*X(2,2)+tt73+tt72-2*Dm(1,1)*Dm(1,2)*X(2,1)
tt75 = Dm(1,2)*tt35+Dm(1,1)*tt32
tt76 = Dm(1,2)*tt5-Dm(1,1)*tt7
tt77 = Dm(1,2)*tt1-Dm(1,1)*tt6
tt78 = 2*Dm(1,2)*tt35+2*Dm(1,1)*tt32
tt79 = tt42+tt41-1
tt80 = area(1,1)*(2*tt77*tt48*coef(5,1)*tt50-2*tt76*tt39*tt47*coe&
&f(5,1)*tt50+2*tt76*tt39*tt47*coef(5,1)+2*tt77*tt79*coef(4,1)+2*tt&
&76*tt39*coef(4,1)+((-tt78*tt39*tt47)-tt76*tt46*tt47-tt77*tt47*tt4&
&5+2*tt76*tt39*tt44*tt45)*coef(3,1)+2*coef(2,1)*tt8*tt75+coef(1,1)&
&*(tt15*tt74+tt71+tt70))
tt81 = Dm(2,2)*tt7+Dm(2,1)*tt5
tt82 = 2*X(1,3)*Dm(2,1)*Dm(2,2)-2*X(1,1)*Dm(2,1)*Dm(2,2)+Dm(1,1)*&
&X(1,2)*Dm(2,2)+tt58+Dm(1,2)*X(1,2)*Dm(2,1)+tt57
tt83 = ((-2*Dm(2,1)*Dm(2,2))+tt61+tt60)*tt18
tt84 = tt1*Dm(2,1)*tt22
tt85 = tt6*Dm(2,2)*tt26
tt86 = 2*coef(2,1)*(tt6*Dm(2,2)+tt1*Dm(2,1))*tt37
tt87 = Dm(2,1)*tt35-Dm(2,2)*tt32
tt88 = 2*Dm(2,2)*tt7+2*Dm(2,1)*tt5
tt89 = (2*tt6*Dm(2,2)+2*tt1*Dm(2,1))*tt48
tt90 = area(1,1)*((-2*tt39*tt87*tt47*coef(5,1)*tt50)+2*tt39*tt87*&
&tt47*coef(5,1)+2*tt39*tt87*coef(4,1)+(tt89-tt46*tt87*tt47-tt88*tt&
&39*tt47+2*tt39*tt87*tt44*tt45)*coef(3,1)+tt86+coef(1,1)*(tt85+tt8&
&4+tt83+2*tt6*Dm(2,2)*tt14+tt15*tt82+2*tt1*Dm(2,1)*tt11)+2*coef(2,&
&1)*tt8*tt81)
tt91 = 2*tt1*Dm(2,1)*tt10*tt20
tt92 = 2*tt6*Dm(2,2)*tt13*tt24
tt93 = 2*Dm(2,1)*Dm(2,2)*X(2,3)+Dm(1,1)*Dm(2,2)*X(2,2)+Dm(1,2)*Dm&
&(2,1)*X(2,2)-2*Dm(2,1)*X(2,1)*Dm(2,2)+tt73+tt72
tt94 = Dm(2,2)*tt35+Dm(2,1)*tt32
tt95 = 2*Dm(2,2)*tt35+2*Dm(2,1)*tt32
tt96 = tt5*Dm(2,2)-Dm(2,1)*tt7
tt97 = tt1*Dm(2,2)-Dm(2,1)*tt6
tt98 = area(1,1)*(2*tt97*tt48*coef(5,1)*tt50-2*tt96*tt39*tt47*coe&
&f(5,1)*tt50+2*tt96*tt39*tt47*coef(5,1)+2*tt97*tt79*coef(4,1)+2*tt&
&96*tt39*coef(4,1)+((-tt96*tt46*tt47)-tt97*tt47*tt45+2*tt96*tt39*t&
&t44*tt45-tt39*tt47*tt95)*coef(3,1)+2*coef(2,1)*tt8*tt94+coef(1,1)&
&*(tt15*tt93+tt92+tt91))
tt99 = tt54**2
tt100 = Dm(1,1)*tt6-Dm(1,2)*tt1
tt101 = area(1,1)*(2*tt100*tt48*coef(5,1)*tt50-2*tt54*tt66*tt47*c&
&oef(5,1)*tt50+2*tt54*tt66*tt47*coef(5,1)+2*tt100*tt79*coef(4,1)+2&
&*tt54*tt66*coef(4,1)+((-tt67*tt54*tt47)-tt100*tt47*tt45+2*tt54*tt&
&66*tt44*tt45-tt66*tt47*tt53)*coef(3,1)+2*coef(2,1)*tt56*tt52+coef&
&(1,1)*(tt59*tt51+tt71+tt70))
tt102 = area(1,1)*((-2*tt76*tt54*tt47*coef(5,1)*tt50)+2*tt76*tt54&
&*tt47*coef(5,1)+2*tt76*tt54*coef(4,1)+(tt68-tt54*tt78*tt47+2*tt76&
&*tt54*tt44*tt45-tt76*tt47*tt53)*coef(3,1)+tt65+coef(1,1)*(tt64+2*&
&Dm(1,2)*tt6*tt25+tt63+2*Dm(1,1)*tt1*tt21+tt62+tt74*tt51)+2*coef(2&
&,1)*tt75*tt52)
tt103 = Dm(2,1)*tt6-tt1*Dm(2,2)
tt104 = area(1,1)*(2*tt103*tt48*coef(5,1)*tt50-2*tt54*tt87*tt47*c&
&oef(5,1)*tt50+2*tt54*tt87*tt47*coef(5,1)+2*tt103*tt79*coef(4,1)+2&
&*tt54*tt87*coef(4,1)+((-tt54*tt88*tt47)-tt103*tt47*tt45+2*tt54*tt&
&87*tt44*tt45-tt87*tt47*tt53)*coef(3,1)+2*coef(2,1)*tt81*tt52+coef&
&(1,1)*(tt82*tt51+tt92+tt91))
tt105 = area(1,1)*((-2*tt54*tt96*tt47*coef(5,1)*tt50)+2*tt54*tt96&
&*tt47*coef(5,1)+2*tt54*tt96*coef(4,1)+(tt89+2*tt54*tt96*tt44*tt45&
&-tt54*tt47*tt95-tt96*tt47*tt53)*coef(3,1)+tt86+coef(1,1)*(tt85+2*&
&tt6*Dm(2,2)*tt25+tt84+2*tt1*Dm(2,1)*tt21+tt83+tt51*tt93)+2*coef(2&
&,1)*tt52*tt94)
tt106 = Dm(1,1)**2
tt107 = Dm(1,2)**2
tt108 = 2*Dm(1,1)*Dm(1,2)*tt18
tt109 = tt106*tt22
tt110 = tt107*tt26
tt111 = 2*(tt107+tt106)*coef(2,1)*tt37
tt112 = tt66**2
tt113 = (2*tt107+2*tt106)*tt48
tt114 = area(1,1)*((-2*tt76*tt66*tt47*coef(5,1)*tt50)+2*tt76*tt66&
&*tt47*coef(5,1)+2*tt76*tt66*coef(4,1)+((-tt66*tt78*tt47)-tt76*tt6&
&7*tt47+2*tt76*tt66*tt44*tt45)*coef(3,1)+2*coef(2,1)*tt56*tt75+coe&
&f(1,1)*(tt59*tt74+2*tt107*tt13*tt24+2*tt106*tt10*tt20))
tt115 = Dm(1,2)*Dm(2,1)
tt116 = Dm(1,1)*Dm(2,2)
tt117 = (tt116+tt115)*tt18
tt118 = Dm(1,1)*Dm(2,1)*tt22
tt119 = Dm(1,2)*Dm(2,2)*tt26
tt120 = 2*coef(2,1)*(Dm(1,2)*Dm(2,2)+Dm(1,1)*Dm(2,1))*tt37
tt121 = (2*Dm(1,2)*Dm(2,2)+2*Dm(1,1)*Dm(2,1))*tt48
tt122 = area(1,1)*((-2*tt66*tt87*tt47*coef(5,1)*tt50)+2*tt66*tt87&
&*tt47*coef(5,1)+2*tt66*tt87*coef(4,1)+(tt121-tt67*tt87*tt47-tt88*&
&tt66*tt47+2*tt66*tt87*tt44*tt45)*coef(3,1)+tt120+coef(1,1)*(tt119&
&+tt118+tt117+2*Dm(1,2)*Dm(2,2)*tt14+tt59*tt82+2*Dm(1,1)*Dm(2,1)*t&
&t11)+2*coef(2,1)*tt56*tt81)
tt123 = 2*Dm(1,1)*Dm(2,1)*tt10*tt20
tt124 = 2*Dm(1,2)*Dm(2,2)*tt13*tt24
tt125 = tt116+tt60
tt126 = area(1,1)*(2*tt125*tt48*coef(5,1)*tt50-2*tt96*tt66*tt47*c&
&oef(5,1)*tt50+2*tt96*tt66*tt47*coef(5,1)+2*tt125*tt79*coef(4,1)+2&
&*tt96*tt66*coef(4,1)+((-tt67*tt96*tt47)-tt125*tt47*tt45+2*tt96*tt&
&66*tt44*tt45-tt66*tt47*tt95)*coef(3,1)+2*coef(2,1)*tt56*tt94+coef&
&(1,1)*(tt59*tt93+tt124+tt123))
tt127 = tt76**2
tt128 = tt61+tt115
tt129 = area(1,1)*(2*tt128*tt48*coef(5,1)*tt50-2*tt76*tt87*tt47*c&
&oef(5,1)*tt50+2*tt76*tt87*tt47*coef(5,1)+2*tt128*tt79*coef(4,1)+2&
&*tt76*tt87*coef(4,1)+((-tt78*tt87*tt47)-tt76*tt88*tt47-tt128*tt47&
&*tt45+2*tt76*tt87*tt44*tt45)*coef(3,1)+2*coef(2,1)*tt81*tt75+coef&
&(1,1)*(tt82*tt74+tt124+tt123))
tt130 = area(1,1)*((-2*tt76*tt96*tt47*coef(5,1)*tt50)+2*tt76*tt96&
&*tt47*coef(5,1)+2*tt76*tt96*coef(4,1)+(tt121-tt96*tt78*tt47+2*tt7&
&6*tt96*tt44*tt45-tt76*tt47*tt95)*coef(3,1)+tt120+coef(1,1)*(tt119&
&+2*Dm(1,2)*Dm(2,2)*tt25+tt118+2*Dm(1,1)*Dm(2,1)*tt21+tt117+tt74*t&
&t93)+2*coef(2,1)*tt75*tt94)
tt131 = Dm(2,1)**2
tt132 = Dm(2,2)**2
tt133 = 2*Dm(2,1)*Dm(2,2)*tt18
tt134 = tt131*tt22
tt135 = tt132*tt26
tt136 = 2*coef(2,1)*(tt132+tt131)*tt37
tt137 = tt87**2
tt138 = (2*tt132+2*tt131)*tt48
tt139 = area(1,1)*((-2*tt96*tt87*tt47*coef(5,1)*tt50)+2*tt96*tt87&
&*tt47*coef(5,1)+2*tt96*tt87*coef(4,1)+((-tt96*tt88*tt47)+2*tt96*t&
&t87*tt44*tt45-tt87*tt47*tt95)*coef(3,1)+2*coef(2,1)*tt81*tt94+coe&
&f(1,1)*(tt82*tt93+2*tt132*tt13*tt24+2*tt131*tt10*tt20))
tt140 = tt96**2
hes(1,1) = area(1,1)*((-2*tt40*tt47*coef(5,1)*tt50)+2*tt40*tt47*c&
&oef(5,1)+2*tt40*coef(4,1)+(tt49-2*tt46*tt39*tt47+2*tt40*tt44*tt45&
&)*coef(3,1)+tt38+coef(1,1)*(tt27+tt23+tt19+tt15**2+2*tt12*tt14+2*&
&tt9*tt11)+2*coef(2,1)*tt8**2)
hes(1,2) = tt55
hes(1,3) = tt69
hes(1,4) = tt80
hes(1,5) = tt90
hes(1,6) = tt98
hes(2,1) = tt55
hes(2,2) = area(1,1)*((-2*tt99*tt47*coef(5,1)*tt50)+2*tt99*tt47*c&
&oef(5,1)+2*tt99*coef(4,1)+(tt49+2*tt99*tt44*tt45-2*tt54*tt47*tt53&
&)*coef(3,1)+2*coef(2,1)*tt52**2+tt38+coef(1,1)*(tt51**2+tt27+2*tt&
&12*tt25+tt23+2*tt9*tt21+tt19))
hes(2,3) = tt101
hes(2,4) = tt102
hes(2,5) = tt104
hes(2,6) = tt105
hes(3,1) = tt69
hes(3,2) = tt101
hes(3,3) = area(1,1)*((-2*tt112*tt47*coef(5,1)*tt50)+2*tt112*tt47&
&*coef(5,1)+2*tt112*coef(4,1)+(tt113-2*tt67*tt66*tt47+2*tt112*tt44&
&*tt45)*coef(3,1)+tt111+coef(1,1)*(tt110+tt109+tt108+tt59**2+2*tt1&
&07*tt14+2*tt106*tt11)+2*coef(2,1)*tt56**2)
hes(3,4) = tt114
hes(3,5) = tt122
hes(3,6) = tt126
hes(4,1) = tt80
hes(4,2) = tt102
hes(4,3) = tt114
hes(4,4) = area(1,1)*((-2*tt127*tt47*coef(5,1)*tt50)+2*tt127*tt47&
&*coef(5,1)+2*tt127*coef(4,1)+(tt113-2*tt76*tt78*tt47+2*tt127*tt44&
&*tt45)*coef(3,1)+2*coef(2,1)*tt75**2+tt111+coef(1,1)*(tt74**2+tt1&
&10+2*tt107*tt25+tt109+2*tt106*tt21+tt108))
hes(4,5) = tt129
hes(4,6) = tt130
hes(5,1) = tt90
hes(5,2) = tt104
hes(5,3) = tt122
hes(5,4) = tt129
hes(5,5) = area(1,1)*((-2*tt137*tt47*coef(5,1)*tt50)+2*tt137*tt47&
&*coef(5,1)+2*tt137*coef(4,1)+(tt138-2*tt88*tt87*tt47+2*tt137*tt44&
&*tt45)*coef(3,1)+tt136+coef(1,1)*(tt135+tt134+tt133+tt82**2+2*tt1&
&32*tt14+2*tt131*tt11)+2*coef(2,1)*tt81**2)
hes(5,6) = tt139
hes(6,1) = tt98
hes(6,2) = tt105
hes(6,3) = tt126
hes(6,4) = tt130
hes(6,5) = tt139
hes(6,6) = area(1,1)*((-2*tt140*tt47*coef(5,1)*tt50)+2*tt140*tt47&
&*coef(5,1)+2*tt140*coef(4,1)+(tt138+2*tt140*tt44*tt45-2*tt96*tt47&
&*tt95)*coef(3,1)+2*coef(2,1)*tt94**2+tt136+coef(1,1)*(tt93**2+tt1&
&35+2*tt132*tt25+tt134+2*tt131*tt21+tt133))
END 
