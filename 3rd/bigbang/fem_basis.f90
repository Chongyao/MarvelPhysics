SUBROUTINE quad4_shift_shape_func(val, eps, midx, midy) 
IMPLICIT NONE 
REAL(KIND=8) val(4, 1) 
REAL(KIND=8) eps(2, 1) 
REAL(KIND=8) midx(1, 2) 
REAL(KIND=8) midy(1, 2) 
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
tt1 = eps(1,1)-1
tt2 = 1/(midx(1,1)+1)
tt3 = 1/midx(1,2)
tt4 = 2*midx(1,2)
tt5 = -midx(1,1)
tt6 = 2*midx(1,1)*midx(1,2)
tt7 = 1/(eps(1,1)*tt2*tt3*(tt4+midx(1,1)-1)-tt2*tt3*(tt6+tt5+1))
tt8 = eps(2,1)-1
tt9 = 1/(midy(1,1)+1)
tt10 = 1/midy(1,2)
tt11 = -midy(1,1)
tt12 = 2*midy(1,1)*midy(1,2)
tt13 = 2*midy(1,2)
tt14 = 1/(tt9*tt10*(tt13+midy(1,1)-1)*eps(2,1)-tt9*tt10*(tt12+tt1&
&1+1))
tt15 = eps(2,1)+1
tt16 = 1/(midy(1,1)-1)
tt17 = 1/(tt16*tt10*(tt12+tt11-1)-tt16*tt10*(tt13+tt11-1)*eps(2,1&
&))
tt18 = eps(1,1)+1
tt19 = 1/(midx(1,1)-1)
tt20 = 1/(tt19*tt3*(tt6+tt5-1)-eps(1,1)*tt19*tt3*(tt4+tt5-1))
val(1,1) = tt1*tt7*tt8*tt14
val(2,1) = tt1*tt7*tt15*tt17
val(3,1) = tt18*tt20*tt8*tt14
val(4,1) = tt18*tt20*tt15*tt17
END 
SUBROUTINE quad4_shift_shape_func_jac(jac, eps, midx, midy) 
IMPLICIT NONE 
REAL(KIND=8) jac(4, 2) 
REAL(KIND=8) eps(2, 1) 
REAL(KIND=8) midx(1, 2) 
REAL(KIND=8) midy(1, 2) 
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
tt1 = eps(1,1)-1
tt2 = 1/(midx(1,1)+1)
tt3 = 1/midx(1,2)
tt4 = 2*midx(1,2)
tt5 = tt4+midx(1,1)-1
tt6 = -midx(1,1)
tt7 = 2*midx(1,1)*midx(1,2)
tt8 = eps(1,1)*tt2*tt3*tt5-tt2*tt3*(tt7+tt6+1)
tt9 = 1/tt8**2
tt10 = eps(2,1)-1
tt11 = 1/(midy(1,1)+1)
tt12 = 1/midy(1,2)
tt13 = -midy(1,1)
tt14 = 2*midy(1,1)*midy(1,2)
tt15 = 2*midy(1,2)
tt16 = tt15+midy(1,1)-1
tt17 = tt11*tt12*tt16*eps(2,1)-tt11*tt12*(tt14+tt13+1)
tt18 = 1/tt17
tt19 = 1/tt8
tt20 = 1/tt17**2
tt21 = eps(2,1)+1
tt22 = 1/(midy(1,1)-1)
tt23 = tt15+tt13-1
tt24 = tt22*tt12*(tt14+tt13-1)-tt22*tt12*tt23*eps(2,1)
tt25 = 1/tt24
tt26 = 1/tt24**2
tt27 = eps(1,1)+1
tt28 = 1/(midx(1,1)-1)
tt29 = tt4+tt6-1
tt30 = tt28*tt3*(tt7+tt6-1)-eps(1,1)*tt28*tt3*tt29
tt31 = 1/tt30**2
tt32 = 1/tt30
jac(1,1) = tt19*tt10*tt18-tt1*tt2*tt3*tt5*tt9*tt10*tt18
jac(1,2) = tt1*tt19*tt18-tt1*tt11*tt19*tt12*tt16*tt10*tt20
jac(2,1) = tt19*tt21*tt25-tt1*tt2*tt3*tt5*tt9*tt21*tt25
jac(2,2) = tt1*tt19*tt25+tt1*tt22*tt19*tt12*tt23*tt21*tt26
jac(3,1) = tt32*tt10*tt18+tt27*tt28*tt3*tt29*tt31*tt10*tt18
jac(3,2) = tt27*tt32*tt18-tt27*tt11*tt32*tt12*tt16*tt10*tt20
jac(4,1) = tt32*tt21*tt25+tt27*tt28*tt3*tt29*tt31*tt21*tt25
jac(4,2) = tt27*tt32*tt25+tt27*tt22*tt32*tt12*tt23*tt21*tt26
END 
SUBROUTINE quad4_shape_function(val, eps, ptx2, pty2) 
IMPLICIT NONE 
REAL(KIND=8) val(4, 1) 
REAL(KIND=8) eps(2, 1) 
REAL(KIND=8) ptx2(2, 1) 
REAL(KIND=8) pty2(2, 1) 
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
tt1 = -ptx2(2,1)
tt2 = tt1+eps(1,1)
tt3 = 1/(tt1+ptx2(1,1))
tt4 = -pty2(2,1)
tt5 = 1/(tt4+pty2(1,1))
tt6 = tt4+eps(2,1)
tt7 = -pty2(1,1)
tt8 = eps(2,1)+tt7
tt9 = 1/(pty2(2,1)+tt7)
tt10 = -ptx2(1,1)
tt11 = tt10+eps(1,1)
tt12 = 1/(ptx2(2,1)+tt10)
val(1,1) = tt2*tt3*tt5*tt6
val(2,1) = tt8*tt2*tt3*tt9
val(3,1) = tt11*tt12*tt5*tt6
val(4,1) = tt11*tt8*tt12*tt9
END 
SUBROUTINE quad4_shape_function_jac(jac, eps, ptx2, pty2) 
IMPLICIT NONE 
REAL(KIND=8) jac(4, 2) 
REAL(KIND=8) eps(2, 1) 
REAL(KIND=8) ptx2(2, 1) 
REAL(KIND=8) pty2(2, 1) 
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
tt1 = -ptx2(2,1)
tt2 = 1/(tt1+ptx2(1,1))
tt3 = -pty2(2,1)
tt4 = 1/(tt3+pty2(1,1))
tt5 = tt3+eps(2,1)
tt6 = tt1+eps(1,1)
tt7 = -pty2(1,1)
tt8 = eps(2,1)+tt7
tt9 = 1/(pty2(2,1)+tt7)
tt10 = -ptx2(1,1)
tt11 = 1/(ptx2(2,1)+tt10)
tt12 = tt10+eps(1,1)
jac(1,1) = tt2*tt4*tt5
jac(1,2) = tt6*tt2*tt4
jac(2,1) = tt8*tt2*tt9
jac(2,2) = tt6*tt2*tt9
jac(3,1) = tt11*tt4*tt5
jac(3,2) = tt12*tt11*tt4
jac(4,1) = tt8*tt11*tt9
jac(4,2) = tt12*tt11*tt9
END 
SUBROUTINE quad4_shape_function_hes(hes, eps, ptx2, pty2) 
IMPLICIT NONE 
REAL(KIND=8) hes(8, 2) 
REAL(KIND=8) eps(2, 1) 
REAL(KIND=8) ptx2(2, 1) 
REAL(KIND=8) pty2(2, 1) 
REAL(KIND=8)  tt1 
REAL(KIND=8)  tt2 
REAL(KIND=8)  tt3 
REAL(KIND=8)  tt4 
REAL(KIND=8)  tt5 
REAL(KIND=8)  tt6 
REAL(KIND=8)  tt7 
REAL(KIND=8)  tt8 
tt1 = 1/(ptx2(1,1)-ptx2(2,1))
tt2 = 1/(pty2(1,1)-pty2(2,1))
tt3 = tt1*tt2
tt4 = 1/(pty2(2,1)-pty2(1,1))
tt5 = tt1*tt4
tt6 = 1/(ptx2(2,1)-ptx2(1,1))
tt7 = tt6*tt2
tt8 = tt6*tt4
hes(1,1) = 0
hes(1,2) = tt3
hes(2,1) = 0
hes(2,2) = tt5
hes(3,1) = 0
hes(3,2) = tt7
hes(4,1) = 0
hes(4,2) = tt8
hes(5,1) = tt3
hes(5,2) = 0
hes(6,1) = tt5
hes(6,2) = 0
hes(7,1) = tt7
hes(7,2) = 0
hes(8,1) = tt8
hes(8,2) = 0
END 
SUBROUTINE quad9_shape_function(val, eps, ptx3, pty3) 
IMPLICIT NONE 
REAL(KIND=8) val(9, 1) 
REAL(KIND=8) eps(2, 1) 
REAL(KIND=8) ptx3(3, 1) 
REAL(KIND=8) pty3(3, 1) 
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
tt1 = -ptx3(2,1)
tt2 = tt1+eps(1,1)
tt3 = 1/(tt1+ptx3(1,1))
tt4 = -pty3(2,1)
tt5 = 1/(tt4+pty3(1,1))
tt6 = tt4+eps(2,1)
tt7 = -ptx3(3,1)
tt8 = tt7+eps(1,1)
tt9 = 1/(tt7+ptx3(1,1))
tt10 = -pty3(3,1)
tt11 = 1/(tt10+pty3(1,1))
tt12 = tt10+eps(2,1)
tt13 = -pty3(1,1)
tt14 = eps(2,1)+tt13
tt15 = 1/(pty3(2,1)+tt13)
tt16 = 1/(tt10+pty3(2,1))
tt17 = 1/(pty3(3,1)+tt13)
tt18 = 1/(pty3(3,1)+tt4)
tt19 = -ptx3(1,1)
tt20 = tt19+eps(1,1)
tt21 = 1/(ptx3(2,1)+tt19)
tt22 = 1/(tt7+ptx3(2,1))
tt23 = 1/(ptx3(3,1)+tt19)
tt24 = 1/(ptx3(3,1)+tt1)
val(1,1) = tt2*tt3*tt5*tt6*tt8*tt9*tt11*tt12
val(2,1) = tt14*tt2*tt3*tt15*tt8*tt9*tt12*tt16
val(3,1) = tt14*tt2*tt3*tt6*tt8*tt9*tt17*tt18
val(4,1) = tt20*tt21*tt5*tt6*tt8*tt22*tt11*tt12
val(5,1) = tt20*tt14*tt21*tt15*tt8*tt22*tt12*tt16
val(6,1) = tt20*tt14*tt21*tt6*tt8*tt22*tt17*tt18
val(7,1) = tt20*tt2*tt5*tt6*tt23*tt24*tt11*tt12
val(8,1) = tt20*tt14*tt2*tt15*tt23*tt24*tt12*tt16
val(9,1) = tt20*tt14*tt2*tt6*tt23*tt24*tt17*tt18
END 
SUBROUTINE quad9_shape_function_jac(jac, eps, ptx3, pty3) 
IMPLICIT NONE 
REAL(KIND=8) jac(9, 2) 
REAL(KIND=8) eps(2, 1) 
REAL(KIND=8) ptx3(3, 1) 
REAL(KIND=8) pty3(3, 1) 
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
tt1 = -ptx3(2,1)
tt2 = tt1+eps(1,1)
tt3 = 1/(tt1+ptx3(1,1))
tt4 = -pty3(2,1)
tt5 = 1/(tt4+pty3(1,1))
tt6 = tt4+eps(2,1)
tt7 = -ptx3(3,1)
tt8 = 1/(tt7+ptx3(1,1))
tt9 = -pty3(3,1)
tt10 = 1/(tt9+pty3(1,1))
tt11 = tt9+eps(2,1)
tt12 = tt7+eps(1,1)
tt13 = -pty3(1,1)
tt14 = eps(2,1)+tt13
tt15 = 1/(pty3(2,1)+tt13)
tt16 = 1/(tt9+pty3(2,1))
tt17 = 1/(pty3(3,1)+tt13)
tt18 = 1/(pty3(3,1)+tt4)
tt19 = -ptx3(1,1)
tt20 = tt19+eps(1,1)
tt21 = 1/(ptx3(2,1)+tt19)
tt22 = 1/(tt7+ptx3(2,1))
tt23 = 1/(ptx3(3,1)+tt19)
tt24 = 1/(ptx3(3,1)+tt1)
jac(1,1) = tt3*tt5*tt6*tt12*tt8*tt10*tt11+tt2*tt3*tt5*tt6*tt8*tt1&
&0*tt11
jac(1,2) = tt2*tt3*tt5*tt12*tt8*tt10*tt11+tt2*tt3*tt5*tt6*tt12*tt&
&8*tt10
jac(2,1) = tt14*tt3*tt15*tt12*tt8*tt11*tt16+tt14*tt2*tt3*tt15*tt8&
&*tt11*tt16
jac(2,2) = tt2*tt3*tt15*tt12*tt8*tt11*tt16+tt14*tt2*tt3*tt15*tt12&
&*tt8*tt16
jac(3,1) = tt14*tt3*tt6*tt12*tt8*tt17*tt18+tt14*tt2*tt3*tt6*tt8*t&
&t17*tt18
jac(3,2) = tt2*tt3*tt6*tt12*tt8*tt17*tt18+tt14*tt2*tt3*tt12*tt8*t&
&t17*tt18
jac(4,1) = tt21*tt5*tt6*tt12*tt22*tt10*tt11+tt20*tt21*tt5*tt6*tt2&
&2*tt10*tt11
jac(4,2) = tt20*tt21*tt5*tt12*tt22*tt10*tt11+tt20*tt21*tt5*tt6*tt&
&12*tt22*tt10
jac(5,1) = tt14*tt21*tt15*tt12*tt22*tt11*tt16+tt20*tt14*tt21*tt15&
&*tt22*tt11*tt16
jac(5,2) = tt20*tt21*tt15*tt12*tt22*tt11*tt16+tt20*tt14*tt21*tt15&
&*tt12*tt22*tt16
jac(6,1) = tt14*tt21*tt6*tt12*tt22*tt17*tt18+tt20*tt14*tt21*tt6*t&
&t22*tt17*tt18
jac(6,2) = tt20*tt21*tt6*tt12*tt22*tt17*tt18+tt20*tt14*tt21*tt12*&
&tt22*tt17*tt18
jac(7,1) = tt2*tt5*tt6*tt23*tt24*tt10*tt11+tt20*tt5*tt6*tt23*tt24&
&*tt10*tt11
jac(7,2) = tt20*tt2*tt5*tt23*tt24*tt10*tt11+tt20*tt2*tt5*tt6*tt23&
&*tt24*tt10
jac(8,1) = tt14*tt2*tt15*tt23*tt24*tt11*tt16+tt20*tt14*tt15*tt23*&
&tt24*tt11*tt16
jac(8,2) = tt20*tt2*tt15*tt23*tt24*tt11*tt16+tt20*tt14*tt2*tt15*t&
&t23*tt24*tt16
jac(9,1) = tt14*tt2*tt6*tt23*tt24*tt17*tt18+tt20*tt14*tt6*tt23*tt&
&24*tt17*tt18
jac(9,2) = tt20*tt2*tt6*tt23*tt24*tt17*tt18+tt20*tt14*tt2*tt23*tt&
&24*tt17*tt18
END 
SUBROUTINE quad9_shape_function_hes(hes, eps, ptx3, pty3) 
IMPLICIT NONE 
REAL(KIND=8) hes(18, 2) 
REAL(KIND=8) eps(2, 1) 
REAL(KIND=8) ptx3(3, 1) 
REAL(KIND=8) pty3(3, 1) 
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
tt1 = -ptx3(2,1)
tt2 = 1/(tt1+ptx3(1,1))
tt3 = -pty3(2,1)
tt4 = 1/(tt3+pty3(1,1))
tt5 = tt3+eps(2,1)
tt6 = -ptx3(3,1)
tt7 = 1/(tt6+ptx3(1,1))
tt8 = -pty3(3,1)
tt9 = 1/(tt8+pty3(1,1))
tt10 = tt8+eps(2,1)
tt11 = tt1+eps(1,1)
tt12 = tt6+eps(1,1)
tt13 = tt2*tt4*tt12*tt7*tt9*tt10+tt11*tt2*tt4*tt7*tt9*tt10+tt2*tt&
&4*tt5*tt12*tt7*tt9+tt11*tt2*tt4*tt5*tt7*tt9
tt14 = -pty3(1,1)
tt15 = eps(2,1)+tt14
tt16 = 1/(pty3(2,1)+tt14)
tt17 = 1/(tt8+pty3(2,1))
tt18 = tt2*tt16*tt12*tt7*tt10*tt17+tt11*tt2*tt16*tt7*tt10*tt17+tt&
&15*tt2*tt16*tt12*tt7*tt17+tt15*tt11*tt2*tt16*tt7*tt17
tt19 = 1/(pty3(3,1)+tt14)
tt20 = 1/(pty3(3,1)+tt3)
tt21 = tt2*tt5*tt12*tt7*tt19*tt20+tt15*tt2*tt12*tt7*tt19*tt20+tt1&
&1*tt2*tt5*tt7*tt19*tt20+tt15*tt11*tt2*tt7*tt19*tt20
tt22 = -ptx3(1,1)
tt23 = 1/(ptx3(2,1)+tt22)
tt24 = 1/(tt6+ptx3(2,1))
tt25 = tt22+eps(1,1)
tt26 = tt23*tt4*tt12*tt24*tt9*tt10+tt25*tt23*tt4*tt24*tt9*tt10+tt&
&23*tt4*tt5*tt12*tt24*tt9+tt25*tt23*tt4*tt5*tt24*tt9
tt27 = tt23*tt16*tt12*tt24*tt10*tt17+tt25*tt23*tt16*tt24*tt10*tt1&
&7+tt15*tt23*tt16*tt12*tt24*tt17+tt25*tt15*tt23*tt16*tt24*tt17
tt28 = tt23*tt5*tt12*tt24*tt19*tt20+tt15*tt23*tt12*tt24*tt19*tt20&
&+tt25*tt23*tt5*tt24*tt19*tt20+tt25*tt15*tt23*tt24*tt19*tt20
tt29 = 1/(ptx3(3,1)+tt22)
tt30 = 1/(ptx3(3,1)+tt1)
tt31 = tt11*tt4*tt29*tt30*tt9*tt10+tt25*tt4*tt29*tt30*tt9*tt10+tt&
&11*tt4*tt5*tt29*tt30*tt9+tt25*tt4*tt5*tt29*tt30*tt9
tt32 = tt11*tt16*tt29*tt30*tt10*tt17+tt25*tt16*tt29*tt30*tt10*tt1&
&7+tt15*tt11*tt16*tt29*tt30*tt17+tt25*tt15*tt16*tt29*tt30*tt17
tt33 = tt11*tt5*tt29*tt30*tt19*tt20+tt25*tt5*tt29*tt30*tt19*tt20+&
&tt15*tt11*tt29*tt30*tt19*tt20+tt25*tt15*tt29*tt30*tt19*tt20
hes(1,1) = 2*tt2*tt4*tt5*tt7*tt9*tt10
hes(1,2) = tt13
hes(2,1) = 2*tt15*tt2*tt16*tt7*tt10*tt17
hes(2,2) = tt18
hes(3,1) = 2*tt15*tt2*tt5*tt7*tt19*tt20
hes(3,2) = tt21
hes(4,1) = 2*tt23*tt4*tt5*tt24*tt9*tt10
hes(4,2) = tt26
hes(5,1) = 2*tt15*tt23*tt16*tt24*tt10*tt17
hes(5,2) = tt27
hes(6,1) = 2*tt15*tt23*tt5*tt24*tt19*tt20
hes(6,2) = tt28
hes(7,1) = 2*tt4*tt5*tt29*tt30*tt9*tt10
hes(7,2) = tt31
hes(8,1) = 2*tt15*tt16*tt29*tt30*tt10*tt17
hes(8,2) = tt32
hes(9,1) = 2*tt15*tt5*tt29*tt30*tt19*tt20
hes(9,2) = tt33
hes(10,1) = tt13
hes(10,2) = 2*tt11*tt2*tt4*tt12*tt7*tt9
hes(11,1) = tt18
hes(11,2) = 2*tt11*tt2*tt16*tt12*tt7*tt17
hes(12,1) = tt21
hes(12,2) = 2*tt11*tt2*tt12*tt7*tt19*tt20
hes(13,1) = tt26
hes(13,2) = 2*tt25*tt23*tt4*tt12*tt24*tt9
hes(14,1) = tt27
hes(14,2) = 2*tt25*tt23*tt16*tt12*tt24*tt17
hes(15,1) = tt28
hes(15,2) = 2*tt25*tt23*tt12*tt24*tt19*tt20
hes(16,1) = tt31
hes(16,2) = 2*tt25*tt11*tt4*tt29*tt30*tt9
hes(17,1) = tt32
hes(17,2) = 2*tt25*tt11*tt16*tt29*tt30*tt17
hes(18,1) = tt33
hes(18,2) = 2*tt25*tt11*tt29*tt30*tt19*tt20
END 
SUBROUTINE quad16_shape_function(val, eps, ptx4, pty4) 
IMPLICIT NONE 
REAL(KIND=8) val(16, 1) 
REAL(KIND=8) eps(2, 1) 
REAL(KIND=8) ptx4(4, 1) 
REAL(KIND=8) pty4(4, 1) 
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
tt1 = -ptx4(2,1)
tt2 = tt1+eps(1,1)
tt3 = 1/(tt1+ptx4(1,1))
tt4 = -pty4(2,1)
tt5 = 1/(tt4+pty4(1,1))
tt6 = tt4+eps(2,1)
tt7 = -ptx4(3,1)
tt8 = tt7+eps(1,1)
tt9 = 1/(tt7+ptx4(1,1))
tt10 = -pty4(3,1)
tt11 = 1/(tt10+pty4(1,1))
tt12 = tt10+eps(2,1)
tt13 = -ptx4(4,1)
tt14 = tt13+eps(1,1)
tt15 = 1/(tt13+ptx4(1,1))
tt16 = -pty4(4,1)
tt17 = 1/(tt16+pty4(1,1))
tt18 = tt16+eps(2,1)
tt19 = -pty4(1,1)
tt20 = eps(2,1)+tt19
tt21 = 1/(pty4(2,1)+tt19)
tt22 = 1/(tt10+pty4(2,1))
tt23 = 1/(tt16+pty4(2,1))
tt24 = 1/(pty4(3,1)+tt19)
tt25 = 1/(pty4(3,1)+tt4)
tt26 = 1/(tt16+pty4(3,1))
tt27 = 1/(pty4(4,1)+tt19)
tt28 = 1/(pty4(4,1)+tt4)
tt29 = 1/(pty4(4,1)+tt10)
tt30 = -ptx4(1,1)
tt31 = tt30+eps(1,1)
tt32 = 1/(ptx4(2,1)+tt30)
tt33 = 1/(tt7+ptx4(2,1))
tt34 = 1/(tt13+ptx4(2,1))
tt35 = 1/(ptx4(3,1)+tt30)
tt36 = 1/(ptx4(3,1)+tt1)
tt37 = 1/(tt13+ptx4(3,1))
tt38 = 1/(ptx4(4,1)+tt30)
tt39 = 1/(ptx4(4,1)+tt1)
tt40 = 1/(ptx4(4,1)+tt7)
val(1,1) = tt2*tt3*tt5*tt6*tt8*tt9*tt11*tt12*tt14*tt15*tt17*tt18
val(2,1) = tt20*tt2*tt3*tt21*tt8*tt9*tt12*tt22*tt14*tt15*tt18*tt2&
&3
val(3,1) = tt20*tt2*tt3*tt6*tt8*tt9*tt24*tt25*tt14*tt15*tt18*tt26&
&
val(4,1) = tt20*tt2*tt3*tt6*tt8*tt9*tt12*tt14*tt15*tt27*tt28*tt29&
&
val(5,1) = tt31*tt32*tt5*tt6*tt8*tt33*tt11*tt12*tt14*tt34*tt17*tt&
&18
val(6,1) = tt31*tt20*tt32*tt21*tt8*tt33*tt12*tt22*tt14*tt34*tt18*&
&tt23
val(7,1) = tt31*tt20*tt32*tt6*tt8*tt33*tt24*tt25*tt14*tt34*tt18*t&
&t26
val(8,1) = tt31*tt20*tt32*tt6*tt8*tt33*tt12*tt14*tt34*tt27*tt28*t&
&t29
val(9,1) = tt31*tt2*tt5*tt6*tt35*tt36*tt11*tt12*tt14*tt37*tt17*tt&
&18
val(10,1) = tt31*tt20*tt2*tt21*tt35*tt36*tt12*tt22*tt14*tt37*tt18&
&*tt23
val(11,1) = tt31*tt20*tt2*tt6*tt35*tt36*tt24*tt25*tt14*tt37*tt18*&
&tt26
val(12,1) = tt31*tt20*tt2*tt6*tt35*tt36*tt12*tt14*tt37*tt27*tt28*&
&tt29
val(13,1) = tt31*tt2*tt5*tt6*tt8*tt11*tt12*tt38*tt39*tt40*tt17*tt&
&18
val(14,1) = tt31*tt20*tt2*tt21*tt8*tt12*tt22*tt38*tt39*tt40*tt18*&
&tt23
val(15,1) = tt31*tt20*tt2*tt6*tt8*tt24*tt25*tt38*tt39*tt40*tt18*t&
&t26
val(16,1) = tt31*tt20*tt2*tt6*tt8*tt12*tt38*tt39*tt40*tt27*tt28*t&
&t29
END 
SUBROUTINE quad16_shape_function_jac(jac, eps, ptx4, pty4) 
IMPLICIT NONE 
REAL(KIND=8) jac(16, 2) 
REAL(KIND=8) eps(2, 1) 
REAL(KIND=8) ptx4(4, 1) 
REAL(KIND=8) pty4(4, 1) 
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
tt1 = -ptx4(2,1)
tt2 = tt1+eps(1,1)
tt3 = 1/(tt1+ptx4(1,1))
tt4 = -pty4(2,1)
tt5 = 1/(tt4+pty4(1,1))
tt6 = tt4+eps(2,1)
tt7 = -ptx4(3,1)
tt8 = tt7+eps(1,1)
tt9 = 1/(tt7+ptx4(1,1))
tt10 = -pty4(3,1)
tt11 = 1/(tt10+pty4(1,1))
tt12 = tt10+eps(2,1)
tt13 = -ptx4(4,1)
tt14 = 1/(tt13+ptx4(1,1))
tt15 = -pty4(4,1)
tt16 = 1/(tt15+pty4(1,1))
tt17 = tt15+eps(2,1)
tt18 = tt13+eps(1,1)
tt19 = -pty4(1,1)
tt20 = eps(2,1)+tt19
tt21 = 1/(pty4(2,1)+tt19)
tt22 = 1/(tt10+pty4(2,1))
tt23 = 1/(tt15+pty4(2,1))
tt24 = 1/(pty4(3,1)+tt19)
tt25 = 1/(pty4(3,1)+tt4)
tt26 = 1/(tt15+pty4(3,1))
tt27 = 1/(pty4(4,1)+tt19)
tt28 = 1/(pty4(4,1)+tt4)
tt29 = 1/(pty4(4,1)+tt10)
tt30 = -ptx4(1,1)
tt31 = tt30+eps(1,1)
tt32 = 1/(ptx4(2,1)+tt30)
tt33 = 1/(tt7+ptx4(2,1))
tt34 = 1/(tt13+ptx4(2,1))
tt35 = 1/(ptx4(3,1)+tt30)
tt36 = 1/(ptx4(3,1)+tt1)
tt37 = 1/(tt13+ptx4(3,1))
tt38 = 1/(ptx4(4,1)+tt30)
tt39 = 1/(ptx4(4,1)+tt1)
tt40 = 1/(ptx4(4,1)+tt7)
jac(1,1) = tt3*tt5*tt6*tt8*tt9*tt11*tt12*tt18*tt14*tt16*tt17+tt2*&
&tt3*tt5*tt6*tt9*tt11*tt12*tt18*tt14*tt16*tt17+tt2*tt3*tt5*tt6*tt8&
&*tt9*tt11*tt12*tt14*tt16*tt17
jac(1,2) = tt2*tt3*tt5*tt8*tt9*tt11*tt12*tt18*tt14*tt16*tt17+tt2*&
&tt3*tt5*tt6*tt8*tt9*tt11*tt18*tt14*tt16*tt17+tt2*tt3*tt5*tt6*tt8*&
&tt9*tt11*tt12*tt18*tt14*tt16
jac(2,1) = tt20*tt3*tt21*tt8*tt9*tt12*tt22*tt18*tt14*tt17*tt23+tt&
&20*tt2*tt3*tt21*tt9*tt12*tt22*tt18*tt14*tt17*tt23+tt20*tt2*tt3*tt&
&21*tt8*tt9*tt12*tt22*tt14*tt17*tt23
jac(2,2) = tt2*tt3*tt21*tt8*tt9*tt12*tt22*tt18*tt14*tt17*tt23+tt2&
&0*tt2*tt3*tt21*tt8*tt9*tt22*tt18*tt14*tt17*tt23+tt20*tt2*tt3*tt21&
&*tt8*tt9*tt12*tt22*tt18*tt14*tt23
jac(3,1) = tt20*tt3*tt6*tt8*tt9*tt24*tt25*tt18*tt14*tt17*tt26+tt2&
&0*tt2*tt3*tt6*tt9*tt24*tt25*tt18*tt14*tt17*tt26+tt20*tt2*tt3*tt6*&
&tt8*tt9*tt24*tt25*tt14*tt17*tt26
jac(3,2) = tt2*tt3*tt6*tt8*tt9*tt24*tt25*tt18*tt14*tt17*tt26+tt20&
&*tt2*tt3*tt8*tt9*tt24*tt25*tt18*tt14*tt17*tt26+tt20*tt2*tt3*tt6*t&
&t8*tt9*tt24*tt25*tt18*tt14*tt26
jac(4,1) = tt20*tt3*tt6*tt8*tt9*tt12*tt18*tt14*tt27*tt28*tt29+tt2&
&0*tt2*tt3*tt6*tt9*tt12*tt18*tt14*tt27*tt28*tt29+tt20*tt2*tt3*tt6*&
&tt8*tt9*tt12*tt14*tt27*tt28*tt29
jac(4,2) = tt2*tt3*tt6*tt8*tt9*tt12*tt18*tt14*tt27*tt28*tt29+tt20&
&*tt2*tt3*tt8*tt9*tt12*tt18*tt14*tt27*tt28*tt29+tt20*tt2*tt3*tt6*t&
&t8*tt9*tt18*tt14*tt27*tt28*tt29
jac(5,1) = tt32*tt5*tt6*tt8*tt33*tt11*tt12*tt18*tt34*tt16*tt17+tt&
&31*tt32*tt5*tt6*tt33*tt11*tt12*tt18*tt34*tt16*tt17+tt31*tt32*tt5*&
&tt6*tt8*tt33*tt11*tt12*tt34*tt16*tt17
jac(5,2) = tt31*tt32*tt5*tt8*tt33*tt11*tt12*tt18*tt34*tt16*tt17+t&
&t31*tt32*tt5*tt6*tt8*tt33*tt11*tt18*tt34*tt16*tt17+tt31*tt32*tt5*&
&tt6*tt8*tt33*tt11*tt12*tt18*tt34*tt16
jac(6,1) = tt20*tt32*tt21*tt8*tt33*tt12*tt22*tt18*tt34*tt17*tt23+&
&tt31*tt20*tt32*tt21*tt33*tt12*tt22*tt18*tt34*tt17*tt23+tt31*tt20*&
&tt32*tt21*tt8*tt33*tt12*tt22*tt34*tt17*tt23
jac(6,2) = tt31*tt32*tt21*tt8*tt33*tt12*tt22*tt18*tt34*tt17*tt23+&
&tt31*tt20*tt32*tt21*tt8*tt33*tt22*tt18*tt34*tt17*tt23+tt31*tt20*t&
&t32*tt21*tt8*tt33*tt12*tt22*tt18*tt34*tt23
jac(7,1) = tt20*tt32*tt6*tt8*tt33*tt24*tt25*tt18*tt34*tt17*tt26+t&
&t31*tt20*tt32*tt6*tt33*tt24*tt25*tt18*tt34*tt17*tt26+tt31*tt20*tt&
&32*tt6*tt8*tt33*tt24*tt25*tt34*tt17*tt26
jac(7,2) = tt31*tt32*tt6*tt8*tt33*tt24*tt25*tt18*tt34*tt17*tt26+t&
&t31*tt20*tt32*tt8*tt33*tt24*tt25*tt18*tt34*tt17*tt26+tt31*tt20*tt&
&32*tt6*tt8*tt33*tt24*tt25*tt18*tt34*tt26
jac(8,1) = tt20*tt32*tt6*tt8*tt33*tt12*tt18*tt34*tt27*tt28*tt29+t&
&t31*tt20*tt32*tt6*tt33*tt12*tt18*tt34*tt27*tt28*tt29+tt31*tt20*tt&
&32*tt6*tt8*tt33*tt12*tt34*tt27*tt28*tt29
jac(8,2) = tt31*tt32*tt6*tt8*tt33*tt12*tt18*tt34*tt27*tt28*tt29+t&
&t31*tt20*tt32*tt8*tt33*tt12*tt18*tt34*tt27*tt28*tt29+tt31*tt20*tt&
&32*tt6*tt8*tt33*tt18*tt34*tt27*tt28*tt29
jac(9,1) = tt2*tt5*tt6*tt35*tt36*tt11*tt12*tt18*tt37*tt16*tt17+tt&
&31*tt5*tt6*tt35*tt36*tt11*tt12*tt18*tt37*tt16*tt17+tt31*tt2*tt5*t&
&t6*tt35*tt36*tt11*tt12*tt37*tt16*tt17
jac(9,2) = tt31*tt2*tt5*tt35*tt36*tt11*tt12*tt18*tt37*tt16*tt17+t&
&t31*tt2*tt5*tt6*tt35*tt36*tt11*tt18*tt37*tt16*tt17+tt31*tt2*tt5*t&
&t6*tt35*tt36*tt11*tt12*tt18*tt37*tt16
jac(10,1) = tt20*tt2*tt21*tt35*tt36*tt12*tt22*tt18*tt37*tt17*tt23&
&+tt31*tt20*tt21*tt35*tt36*tt12*tt22*tt18*tt37*tt17*tt23+tt31*tt20&
&*tt2*tt21*tt35*tt36*tt12*tt22*tt37*tt17*tt23
jac(10,2) = tt31*tt2*tt21*tt35*tt36*tt12*tt22*tt18*tt37*tt17*tt23&
&+tt31*tt20*tt2*tt21*tt35*tt36*tt22*tt18*tt37*tt17*tt23+tt31*tt20*&
&tt2*tt21*tt35*tt36*tt12*tt22*tt18*tt37*tt23
jac(11,1) = tt20*tt2*tt6*tt35*tt36*tt24*tt25*tt18*tt37*tt17*tt26+&
&tt31*tt20*tt6*tt35*tt36*tt24*tt25*tt18*tt37*tt17*tt26+tt31*tt20*t&
&t2*tt6*tt35*tt36*tt24*tt25*tt37*tt17*tt26
jac(11,2) = tt31*tt2*tt6*tt35*tt36*tt24*tt25*tt18*tt37*tt17*tt26+&
&tt31*tt20*tt2*tt35*tt36*tt24*tt25*tt18*tt37*tt17*tt26+tt31*tt20*t&
&t2*tt6*tt35*tt36*tt24*tt25*tt18*tt37*tt26
jac(12,1) = tt20*tt2*tt6*tt35*tt36*tt12*tt18*tt37*tt27*tt28*tt29+&
&tt31*tt20*tt6*tt35*tt36*tt12*tt18*tt37*tt27*tt28*tt29+tt31*tt20*t&
&t2*tt6*tt35*tt36*tt12*tt37*tt27*tt28*tt29
jac(12,2) = tt31*tt2*tt6*tt35*tt36*tt12*tt18*tt37*tt27*tt28*tt29+&
&tt31*tt20*tt2*tt35*tt36*tt12*tt18*tt37*tt27*tt28*tt29+tt31*tt20*t&
&t2*tt6*tt35*tt36*tt18*tt37*tt27*tt28*tt29
jac(13,1) = tt2*tt5*tt6*tt8*tt11*tt12*tt38*tt39*tt40*tt16*tt17+tt&
&31*tt5*tt6*tt8*tt11*tt12*tt38*tt39*tt40*tt16*tt17+tt31*tt2*tt5*tt&
&6*tt11*tt12*tt38*tt39*tt40*tt16*tt17
jac(13,2) = tt31*tt2*tt5*tt8*tt11*tt12*tt38*tt39*tt40*tt16*tt17+t&
&t31*tt2*tt5*tt6*tt8*tt11*tt38*tt39*tt40*tt16*tt17+tt31*tt2*tt5*tt&
&6*tt8*tt11*tt12*tt38*tt39*tt40*tt16
jac(14,1) = tt20*tt2*tt21*tt8*tt12*tt22*tt38*tt39*tt40*tt17*tt23+&
&tt31*tt20*tt21*tt8*tt12*tt22*tt38*tt39*tt40*tt17*tt23+tt31*tt20*t&
&t2*tt21*tt12*tt22*tt38*tt39*tt40*tt17*tt23
jac(14,2) = tt31*tt2*tt21*tt8*tt12*tt22*tt38*tt39*tt40*tt17*tt23+&
&tt31*tt20*tt2*tt21*tt8*tt22*tt38*tt39*tt40*tt17*tt23+tt31*tt20*tt&
&2*tt21*tt8*tt12*tt22*tt38*tt39*tt40*tt23
jac(15,1) = tt20*tt2*tt6*tt8*tt24*tt25*tt38*tt39*tt40*tt17*tt26+t&
&t31*tt20*tt6*tt8*tt24*tt25*tt38*tt39*tt40*tt17*tt26+tt31*tt20*tt2&
&*tt6*tt24*tt25*tt38*tt39*tt40*tt17*tt26
jac(15,2) = tt31*tt2*tt6*tt8*tt24*tt25*tt38*tt39*tt40*tt17*tt26+t&
&t31*tt20*tt2*tt8*tt24*tt25*tt38*tt39*tt40*tt17*tt26+tt31*tt20*tt2&
&*tt6*tt8*tt24*tt25*tt38*tt39*tt40*tt26
jac(16,1) = tt20*tt2*tt6*tt8*tt12*tt38*tt39*tt40*tt27*tt28*tt29+t&
&t31*tt20*tt6*tt8*tt12*tt38*tt39*tt40*tt27*tt28*tt29+tt31*tt20*tt2&
&*tt6*tt12*tt38*tt39*tt40*tt27*tt28*tt29
jac(16,2) = tt31*tt2*tt6*tt8*tt12*tt38*tt39*tt40*tt27*tt28*tt29+t&
&t31*tt20*tt2*tt8*tt12*tt38*tt39*tt40*tt27*tt28*tt29+tt31*tt20*tt2&
&*tt6*tt8*tt38*tt39*tt40*tt27*tt28*tt29
END 
SUBROUTINE quad16_shape_function_hes(hes, eps, ptx4, pty4) 
IMPLICIT NONE 
REAL(KIND=8) hes(32, 2) 
REAL(KIND=8) eps(2, 1) 
REAL(KIND=8) ptx4(4, 1) 
REAL(KIND=8) pty4(4, 1) 
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
tt1 = -ptx4(2,1)
tt2 = tt1+eps(1,1)
tt3 = 1/(tt1+ptx4(1,1))
tt4 = -pty4(2,1)
tt5 = 1/(tt4+pty4(1,1))
tt6 = tt4+eps(2,1)
tt7 = -ptx4(3,1)
tt8 = 1/(tt7+ptx4(1,1))
tt9 = -pty4(3,1)
tt10 = 1/(tt9+pty4(1,1))
tt11 = tt9+eps(2,1)
tt12 = -ptx4(4,1)
tt13 = 1/(tt12+ptx4(1,1))
tt14 = -pty4(4,1)
tt15 = 1/(tt14+pty4(1,1))
tt16 = tt14+eps(2,1)
tt17 = tt7+eps(1,1)
tt18 = tt12+eps(1,1)
tt19 = tt3*tt5*tt17*tt8*tt10*tt11*tt18*tt13*tt15*tt16+tt2*tt3*tt5&
&*tt8*tt10*tt11*tt18*tt13*tt15*tt16+tt3*tt5*tt6*tt17*tt8*tt10*tt18&
&*tt13*tt15*tt16+tt2*tt3*tt5*tt6*tt8*tt10*tt18*tt13*tt15*tt16+tt2*&
&tt3*tt5*tt17*tt8*tt10*tt11*tt13*tt15*tt16+tt2*tt3*tt5*tt6*tt17*tt&
&8*tt10*tt13*tt15*tt16+tt3*tt5*tt6*tt17*tt8*tt10*tt11*tt18*tt13*tt&
&15+tt2*tt3*tt5*tt6*tt8*tt10*tt11*tt18*tt13*tt15+tt2*tt3*tt5*tt6*t&
&t17*tt8*tt10*tt11*tt13*tt15
tt20 = -pty4(1,1)
tt21 = eps(2,1)+tt20
tt22 = 1/(pty4(2,1)+tt20)
tt23 = 1/(tt9+pty4(2,1))
tt24 = 1/(tt14+pty4(2,1))
tt25 = tt3*tt22*tt17*tt8*tt11*tt23*tt18*tt13*tt16*tt24+tt2*tt3*tt&
&22*tt8*tt11*tt23*tt18*tt13*tt16*tt24+tt21*tt3*tt22*tt17*tt8*tt23*&
&tt18*tt13*tt16*tt24+tt21*tt2*tt3*tt22*tt8*tt23*tt18*tt13*tt16*tt2&
&4+tt2*tt3*tt22*tt17*tt8*tt11*tt23*tt13*tt16*tt24+tt21*tt2*tt3*tt2&
&2*tt17*tt8*tt23*tt13*tt16*tt24+tt21*tt3*tt22*tt17*tt8*tt11*tt23*t&
&t18*tt13*tt24+tt21*tt2*tt3*tt22*tt8*tt11*tt23*tt18*tt13*tt24+tt21&
&*tt2*tt3*tt22*tt17*tt8*tt11*tt23*tt13*tt24
tt26 = 1/(pty4(3,1)+tt20)
tt27 = 1/(pty4(3,1)+tt4)
tt28 = 1/(tt14+pty4(3,1))
tt29 = tt3*tt6*tt17*tt8*tt26*tt27*tt18*tt13*tt16*tt28+tt21*tt3*tt&
&17*tt8*tt26*tt27*tt18*tt13*tt16*tt28+tt2*tt3*tt6*tt8*tt26*tt27*tt&
&18*tt13*tt16*tt28+tt21*tt2*tt3*tt8*tt26*tt27*tt18*tt13*tt16*tt28+&
&tt2*tt3*tt6*tt17*tt8*tt26*tt27*tt13*tt16*tt28+tt21*tt2*tt3*tt17*t&
&t8*tt26*tt27*tt13*tt16*tt28+tt21*tt3*tt6*tt17*tt8*tt26*tt27*tt18*&
&tt13*tt28+tt21*tt2*tt3*tt6*tt8*tt26*tt27*tt18*tt13*tt28+tt21*tt2*&
&tt3*tt6*tt17*tt8*tt26*tt27*tt13*tt28
tt30 = 1/(pty4(4,1)+tt20)
tt31 = 1/(pty4(4,1)+tt4)
tt32 = 1/(pty4(4,1)+tt9)
tt33 = tt3*tt6*tt17*tt8*tt11*tt18*tt13*tt30*tt31*tt32+tt21*tt3*tt&
&17*tt8*tt11*tt18*tt13*tt30*tt31*tt32+tt2*tt3*tt6*tt8*tt11*tt18*tt&
&13*tt30*tt31*tt32+tt21*tt2*tt3*tt8*tt11*tt18*tt13*tt30*tt31*tt32+&
&tt21*tt3*tt6*tt17*tt8*tt18*tt13*tt30*tt31*tt32+tt21*tt2*tt3*tt6*t&
&t8*tt18*tt13*tt30*tt31*tt32+tt2*tt3*tt6*tt17*tt8*tt11*tt13*tt30*t&
&t31*tt32+tt21*tt2*tt3*tt17*tt8*tt11*tt13*tt30*tt31*tt32+tt21*tt2*&
&tt3*tt6*tt17*tt8*tt13*tt30*tt31*tt32
tt34 = -ptx4(1,1)
tt35 = tt34+eps(1,1)
tt36 = 1/(ptx4(2,1)+tt34)
tt37 = 1/(tt7+ptx4(2,1))
tt38 = 1/(tt12+ptx4(2,1))
tt39 = tt36*tt5*tt17*tt37*tt10*tt11*tt18*tt38*tt15*tt16+tt35*tt36&
&*tt5*tt37*tt10*tt11*tt18*tt38*tt15*tt16+tt36*tt5*tt6*tt17*tt37*tt&
&10*tt18*tt38*tt15*tt16+tt35*tt36*tt5*tt6*tt37*tt10*tt18*tt38*tt15&
&*tt16+tt35*tt36*tt5*tt17*tt37*tt10*tt11*tt38*tt15*tt16+tt35*tt36*&
&tt5*tt6*tt17*tt37*tt10*tt38*tt15*tt16+tt36*tt5*tt6*tt17*tt37*tt10&
&*tt11*tt18*tt38*tt15+tt35*tt36*tt5*tt6*tt37*tt10*tt11*tt18*tt38*t&
&t15+tt35*tt36*tt5*tt6*tt17*tt37*tt10*tt11*tt38*tt15
tt40 = tt36*tt22*tt17*tt37*tt11*tt23*tt18*tt38*tt16*tt24+tt35*tt3&
&6*tt22*tt37*tt11*tt23*tt18*tt38*tt16*tt24+tt21*tt36*tt22*tt17*tt3&
&7*tt23*tt18*tt38*tt16*tt24+tt35*tt21*tt36*tt22*tt37*tt23*tt18*tt3&
&8*tt16*tt24+tt35*tt36*tt22*tt17*tt37*tt11*tt23*tt38*tt16*tt24+tt3&
&5*tt21*tt36*tt22*tt17*tt37*tt23*tt38*tt16*tt24+tt21*tt36*tt22*tt1&
&7*tt37*tt11*tt23*tt18*tt38*tt24+tt35*tt21*tt36*tt22*tt37*tt11*tt2&
&3*tt18*tt38*tt24+tt35*tt21*tt36*tt22*tt17*tt37*tt11*tt23*tt38*tt2&
&4
tt41 = tt36*tt6*tt17*tt37*tt26*tt27*tt18*tt38*tt16*tt28+tt21*tt36&
&*tt17*tt37*tt26*tt27*tt18*tt38*tt16*tt28+tt35*tt36*tt6*tt37*tt26*&
&tt27*tt18*tt38*tt16*tt28+tt35*tt21*tt36*tt37*tt26*tt27*tt18*tt38*&
&tt16*tt28+tt35*tt36*tt6*tt17*tt37*tt26*tt27*tt38*tt16*tt28+tt35*t&
&t21*tt36*tt17*tt37*tt26*tt27*tt38*tt16*tt28+tt21*tt36*tt6*tt17*tt&
&37*tt26*tt27*tt18*tt38*tt28+tt35*tt21*tt36*tt6*tt37*tt26*tt27*tt1&
&8*tt38*tt28+tt35*tt21*tt36*tt6*tt17*tt37*tt26*tt27*tt38*tt28
tt42 = tt36*tt6*tt17*tt37*tt11*tt18*tt38*tt30*tt31*tt32+tt21*tt36&
&*tt17*tt37*tt11*tt18*tt38*tt30*tt31*tt32+tt35*tt36*tt6*tt37*tt11*&
&tt18*tt38*tt30*tt31*tt32+tt35*tt21*tt36*tt37*tt11*tt18*tt38*tt30*&
&tt31*tt32+tt21*tt36*tt6*tt17*tt37*tt18*tt38*tt30*tt31*tt32+tt35*t&
&t21*tt36*tt6*tt37*tt18*tt38*tt30*tt31*tt32+tt35*tt36*tt6*tt17*tt3&
&7*tt11*tt38*tt30*tt31*tt32+tt35*tt21*tt36*tt17*tt37*tt11*tt38*tt3&
&0*tt31*tt32+tt35*tt21*tt36*tt6*tt17*tt37*tt38*tt30*tt31*tt32
tt43 = 1/(ptx4(3,1)+tt34)
tt44 = 1/(ptx4(3,1)+tt1)
tt45 = 1/(tt12+ptx4(3,1))
tt46 = tt2*tt5*tt43*tt44*tt10*tt11*tt18*tt45*tt15*tt16+tt35*tt5*t&
&t43*tt44*tt10*tt11*tt18*tt45*tt15*tt16+tt2*tt5*tt6*tt43*tt44*tt10&
&*tt18*tt45*tt15*tt16+tt35*tt5*tt6*tt43*tt44*tt10*tt18*tt45*tt15*t&
&t16+tt35*tt2*tt5*tt43*tt44*tt10*tt11*tt45*tt15*tt16+tt35*tt2*tt5*&
&tt6*tt43*tt44*tt10*tt45*tt15*tt16+tt2*tt5*tt6*tt43*tt44*tt10*tt11&
&*tt18*tt45*tt15+tt35*tt5*tt6*tt43*tt44*tt10*tt11*tt18*tt45*tt15+t&
&t35*tt2*tt5*tt6*tt43*tt44*tt10*tt11*tt45*tt15
tt47 = tt2*tt22*tt43*tt44*tt11*tt23*tt18*tt45*tt16*tt24+tt35*tt22&
&*tt43*tt44*tt11*tt23*tt18*tt45*tt16*tt24+tt21*tt2*tt22*tt43*tt44*&
&tt23*tt18*tt45*tt16*tt24+tt35*tt21*tt22*tt43*tt44*tt23*tt18*tt45*&
&tt16*tt24+tt35*tt2*tt22*tt43*tt44*tt11*tt23*tt45*tt16*tt24+tt35*t&
&t21*tt2*tt22*tt43*tt44*tt23*tt45*tt16*tt24+tt21*tt2*tt22*tt43*tt4&
&4*tt11*tt23*tt18*tt45*tt24+tt35*tt21*tt22*tt43*tt44*tt11*tt23*tt1&
&8*tt45*tt24+tt35*tt21*tt2*tt22*tt43*tt44*tt11*tt23*tt45*tt24
tt48 = tt2*tt6*tt43*tt44*tt26*tt27*tt18*tt45*tt16*tt28+tt35*tt6*t&
&t43*tt44*tt26*tt27*tt18*tt45*tt16*tt28+tt21*tt2*tt43*tt44*tt26*tt&
&27*tt18*tt45*tt16*tt28+tt35*tt21*tt43*tt44*tt26*tt27*tt18*tt45*tt&
&16*tt28+tt35*tt2*tt6*tt43*tt44*tt26*tt27*tt45*tt16*tt28+tt35*tt21&
&*tt2*tt43*tt44*tt26*tt27*tt45*tt16*tt28+tt21*tt2*tt6*tt43*tt44*tt&
&26*tt27*tt18*tt45*tt28+tt35*tt21*tt6*tt43*tt44*tt26*tt27*tt18*tt4&
&5*tt28+tt35*tt21*tt2*tt6*tt43*tt44*tt26*tt27*tt45*tt28
tt49 = tt2*tt6*tt43*tt44*tt11*tt18*tt45*tt30*tt31*tt32+tt35*tt6*t&
&t43*tt44*tt11*tt18*tt45*tt30*tt31*tt32+tt21*tt2*tt43*tt44*tt11*tt&
&18*tt45*tt30*tt31*tt32+tt35*tt21*tt43*tt44*tt11*tt18*tt45*tt30*tt&
&31*tt32+tt21*tt2*tt6*tt43*tt44*tt18*tt45*tt30*tt31*tt32+tt35*tt21&
&*tt6*tt43*tt44*tt18*tt45*tt30*tt31*tt32+tt35*tt2*tt6*tt43*tt44*tt&
&11*tt45*tt30*tt31*tt32+tt35*tt21*tt2*tt43*tt44*tt11*tt45*tt30*tt3&
&1*tt32+tt35*tt21*tt2*tt6*tt43*tt44*tt45*tt30*tt31*tt32
tt50 = 1/(ptx4(4,1)+tt34)
tt51 = 1/(ptx4(4,1)+tt1)
tt52 = 1/(ptx4(4,1)+tt7)
tt53 = tt2*tt5*tt17*tt10*tt11*tt50*tt51*tt52*tt15*tt16+tt35*tt5*t&
&t17*tt10*tt11*tt50*tt51*tt52*tt15*tt16+tt35*tt2*tt5*tt10*tt11*tt5&
&0*tt51*tt52*tt15*tt16+tt2*tt5*tt6*tt17*tt10*tt50*tt51*tt52*tt15*t&
&t16+tt35*tt5*tt6*tt17*tt10*tt50*tt51*tt52*tt15*tt16+tt35*tt2*tt5*&
&tt6*tt10*tt50*tt51*tt52*tt15*tt16+tt2*tt5*tt6*tt17*tt10*tt11*tt50&
&*tt51*tt52*tt15+tt35*tt5*tt6*tt17*tt10*tt11*tt50*tt51*tt52*tt15+t&
&t35*tt2*tt5*tt6*tt10*tt11*tt50*tt51*tt52*tt15
tt54 = tt2*tt22*tt17*tt11*tt23*tt50*tt51*tt52*tt16*tt24+tt35*tt22&
&*tt17*tt11*tt23*tt50*tt51*tt52*tt16*tt24+tt35*tt2*tt22*tt11*tt23*&
&tt50*tt51*tt52*tt16*tt24+tt21*tt2*tt22*tt17*tt23*tt50*tt51*tt52*t&
&t16*tt24+tt35*tt21*tt22*tt17*tt23*tt50*tt51*tt52*tt16*tt24+tt35*t&
&t21*tt2*tt22*tt23*tt50*tt51*tt52*tt16*tt24+tt21*tt2*tt22*tt17*tt1&
&1*tt23*tt50*tt51*tt52*tt24+tt35*tt21*tt22*tt17*tt11*tt23*tt50*tt5&
&1*tt52*tt24+tt35*tt21*tt2*tt22*tt11*tt23*tt50*tt51*tt52*tt24
tt55 = tt2*tt6*tt17*tt26*tt27*tt50*tt51*tt52*tt16*tt28+tt35*tt6*t&
&t17*tt26*tt27*tt50*tt51*tt52*tt16*tt28+tt21*tt2*tt17*tt26*tt27*tt&
&50*tt51*tt52*tt16*tt28+tt35*tt21*tt17*tt26*tt27*tt50*tt51*tt52*tt&
&16*tt28+tt35*tt2*tt6*tt26*tt27*tt50*tt51*tt52*tt16*tt28+tt35*tt21&
&*tt2*tt26*tt27*tt50*tt51*tt52*tt16*tt28+tt21*tt2*tt6*tt17*tt26*tt&
&27*tt50*tt51*tt52*tt28+tt35*tt21*tt6*tt17*tt26*tt27*tt50*tt51*tt5&
&2*tt28+tt35*tt21*tt2*tt6*tt26*tt27*tt50*tt51*tt52*tt28
tt56 = tt2*tt6*tt17*tt11*tt50*tt51*tt52*tt30*tt31*tt32+tt35*tt6*t&
&t17*tt11*tt50*tt51*tt52*tt30*tt31*tt32+tt21*tt2*tt17*tt11*tt50*tt&
&51*tt52*tt30*tt31*tt32+tt35*tt21*tt17*tt11*tt50*tt51*tt52*tt30*tt&
&31*tt32+tt35*tt2*tt6*tt11*tt50*tt51*tt52*tt30*tt31*tt32+tt35*tt21&
&*tt2*tt11*tt50*tt51*tt52*tt30*tt31*tt32+tt21*tt2*tt6*tt17*tt50*tt&
&51*tt52*tt30*tt31*tt32+tt35*tt21*tt6*tt17*tt50*tt51*tt52*tt30*tt3&
&1*tt32+tt35*tt21*tt2*tt6*tt50*tt51*tt52*tt30*tt31*tt32
hes(1,1) = 2*tt3*tt5*tt6*tt8*tt10*tt11*tt18*tt13*tt15*tt16+2*tt3*&
&tt5*tt6*tt17*tt8*tt10*tt11*tt13*tt15*tt16+2*tt2*tt3*tt5*tt6*tt8*t&
&t10*tt11*tt13*tt15*tt16
hes(1,2) = tt19
hes(2,1) = 2*tt21*tt3*tt22*tt8*tt11*tt23*tt18*tt13*tt16*tt24+2*tt&
&21*tt3*tt22*tt17*tt8*tt11*tt23*tt13*tt16*tt24+2*tt21*tt2*tt3*tt22&
&*tt8*tt11*tt23*tt13*tt16*tt24
hes(2,2) = tt25
hes(3,1) = 2*tt21*tt3*tt6*tt8*tt26*tt27*tt18*tt13*tt16*tt28+2*tt2&
&1*tt3*tt6*tt17*tt8*tt26*tt27*tt13*tt16*tt28+2*tt21*tt2*tt3*tt6*tt&
&8*tt26*tt27*tt13*tt16*tt28
hes(3,2) = tt29
hes(4,1) = 2*tt21*tt3*tt6*tt8*tt11*tt18*tt13*tt30*tt31*tt32+2*tt2&
&1*tt3*tt6*tt17*tt8*tt11*tt13*tt30*tt31*tt32+2*tt21*tt2*tt3*tt6*tt&
&8*tt11*tt13*tt30*tt31*tt32
hes(4,2) = tt33
hes(5,1) = 2*tt36*tt5*tt6*tt37*tt10*tt11*tt18*tt38*tt15*tt16+2*tt&
&36*tt5*tt6*tt17*tt37*tt10*tt11*tt38*tt15*tt16+2*tt35*tt36*tt5*tt6&
&*tt37*tt10*tt11*tt38*tt15*tt16
hes(5,2) = tt39
hes(6,1) = 2*tt21*tt36*tt22*tt37*tt11*tt23*tt18*tt38*tt16*tt24+2*&
&tt21*tt36*tt22*tt17*tt37*tt11*tt23*tt38*tt16*tt24+2*tt35*tt21*tt3&
&6*tt22*tt37*tt11*tt23*tt38*tt16*tt24
hes(6,2) = tt40
hes(7,1) = 2*tt21*tt36*tt6*tt37*tt26*tt27*tt18*tt38*tt16*tt28+2*t&
&t21*tt36*tt6*tt17*tt37*tt26*tt27*tt38*tt16*tt28+2*tt35*tt21*tt36*&
&tt6*tt37*tt26*tt27*tt38*tt16*tt28
hes(7,2) = tt41
hes(8,1) = 2*tt21*tt36*tt6*tt37*tt11*tt18*tt38*tt30*tt31*tt32+2*t&
&t21*tt36*tt6*tt17*tt37*tt11*tt38*tt30*tt31*tt32+2*tt35*tt21*tt36*&
&tt6*tt37*tt11*tt38*tt30*tt31*tt32
hes(8,2) = tt42
hes(9,1) = 2*tt5*tt6*tt43*tt44*tt10*tt11*tt18*tt45*tt15*tt16+2*tt&
&2*tt5*tt6*tt43*tt44*tt10*tt11*tt45*tt15*tt16+2*tt35*tt5*tt6*tt43*&
&tt44*tt10*tt11*tt45*tt15*tt16
hes(9,2) = tt46
hes(10,1) = 2*tt21*tt22*tt43*tt44*tt11*tt23*tt18*tt45*tt16*tt24+2&
&*tt21*tt2*tt22*tt43*tt44*tt11*tt23*tt45*tt16*tt24+2*tt35*tt21*tt2&
&2*tt43*tt44*tt11*tt23*tt45*tt16*tt24
hes(10,2) = tt47
hes(11,1) = 2*tt21*tt6*tt43*tt44*tt26*tt27*tt18*tt45*tt16*tt28+2*&
&tt21*tt2*tt6*tt43*tt44*tt26*tt27*tt45*tt16*tt28+2*tt35*tt21*tt6*t&
&t43*tt44*tt26*tt27*tt45*tt16*tt28
hes(11,2) = tt48
hes(12,1) = 2*tt21*tt6*tt43*tt44*tt11*tt18*tt45*tt30*tt31*tt32+2*&
&tt21*tt2*tt6*tt43*tt44*tt11*tt45*tt30*tt31*tt32+2*tt35*tt21*tt6*t&
&t43*tt44*tt11*tt45*tt30*tt31*tt32
hes(12,2) = tt49
hes(13,1) = 2*tt5*tt6*tt17*tt10*tt11*tt50*tt51*tt52*tt15*tt16+2*t&
&t2*tt5*tt6*tt10*tt11*tt50*tt51*tt52*tt15*tt16+2*tt35*tt5*tt6*tt10&
&*tt11*tt50*tt51*tt52*tt15*tt16
hes(13,2) = tt53
hes(14,1) = 2*tt21*tt22*tt17*tt11*tt23*tt50*tt51*tt52*tt16*tt24+2&
&*tt21*tt2*tt22*tt11*tt23*tt50*tt51*tt52*tt16*tt24+2*tt35*tt21*tt2&
&2*tt11*tt23*tt50*tt51*tt52*tt16*tt24
hes(14,2) = tt54
hes(15,1) = 2*tt21*tt6*tt17*tt26*tt27*tt50*tt51*tt52*tt16*tt28+2*&
&tt21*tt2*tt6*tt26*tt27*tt50*tt51*tt52*tt16*tt28+2*tt35*tt21*tt6*t&
&t26*tt27*tt50*tt51*tt52*tt16*tt28
hes(15,2) = tt55
hes(16,1) = 2*tt21*tt6*tt17*tt11*tt50*tt51*tt52*tt30*tt31*tt32+2*&
&tt21*tt2*tt6*tt11*tt50*tt51*tt52*tt30*tt31*tt32+2*tt35*tt21*tt6*t&
&t11*tt50*tt51*tt52*tt30*tt31*tt32
hes(16,2) = tt56
hes(17,1) = tt19
hes(17,2) = 2*tt2*tt3*tt5*tt17*tt8*tt10*tt18*tt13*tt15*tt16+2*tt2&
&*tt3*tt5*tt17*tt8*tt10*tt11*tt18*tt13*tt15+2*tt2*tt3*tt5*tt6*tt17&
&*tt8*tt10*tt18*tt13*tt15
hes(18,1) = tt25
hes(18,2) = 2*tt2*tt3*tt22*tt17*tt8*tt23*tt18*tt13*tt16*tt24+2*tt&
&2*tt3*tt22*tt17*tt8*tt11*tt23*tt18*tt13*tt24+2*tt21*tt2*tt3*tt22*&
&tt17*tt8*tt23*tt18*tt13*tt24
hes(19,1) = tt29
hes(19,2) = 2*tt2*tt3*tt17*tt8*tt26*tt27*tt18*tt13*tt16*tt28+2*tt&
&2*tt3*tt6*tt17*tt8*tt26*tt27*tt18*tt13*tt28+2*tt21*tt2*tt3*tt17*t&
&t8*tt26*tt27*tt18*tt13*tt28
hes(20,1) = tt33
hes(20,2) = 2*tt2*tt3*tt17*tt8*tt11*tt18*tt13*tt30*tt31*tt32+2*tt&
&2*tt3*tt6*tt17*tt8*tt18*tt13*tt30*tt31*tt32+2*tt21*tt2*tt3*tt17*t&
&t8*tt18*tt13*tt30*tt31*tt32
hes(21,1) = tt39
hes(21,2) = 2*tt35*tt36*tt5*tt17*tt37*tt10*tt18*tt38*tt15*tt16+2*&
&tt35*tt36*tt5*tt17*tt37*tt10*tt11*tt18*tt38*tt15+2*tt35*tt36*tt5*&
&tt6*tt17*tt37*tt10*tt18*tt38*tt15
hes(22,1) = tt40
hes(22,2) = 2*tt35*tt36*tt22*tt17*tt37*tt23*tt18*tt38*tt16*tt24+2&
&*tt35*tt36*tt22*tt17*tt37*tt11*tt23*tt18*tt38*tt24+2*tt35*tt21*tt&
&36*tt22*tt17*tt37*tt23*tt18*tt38*tt24
hes(23,1) = tt41
hes(23,2) = 2*tt35*tt36*tt17*tt37*tt26*tt27*tt18*tt38*tt16*tt28+2&
&*tt35*tt36*tt6*tt17*tt37*tt26*tt27*tt18*tt38*tt28+2*tt35*tt21*tt3&
&6*tt17*tt37*tt26*tt27*tt18*tt38*tt28
hes(24,1) = tt42
hes(24,2) = 2*tt35*tt36*tt17*tt37*tt11*tt18*tt38*tt30*tt31*tt32+2&
&*tt35*tt36*tt6*tt17*tt37*tt18*tt38*tt30*tt31*tt32+2*tt35*tt21*tt3&
&6*tt17*tt37*tt18*tt38*tt30*tt31*tt32
hes(25,1) = tt46
hes(25,2) = 2*tt35*tt2*tt5*tt43*tt44*tt10*tt18*tt45*tt15*tt16+2*t&
&t35*tt2*tt5*tt43*tt44*tt10*tt11*tt18*tt45*tt15+2*tt35*tt2*tt5*tt6&
&*tt43*tt44*tt10*tt18*tt45*tt15
hes(26,1) = tt47
hes(26,2) = 2*tt35*tt2*tt22*tt43*tt44*tt23*tt18*tt45*tt16*tt24+2*&
&tt35*tt2*tt22*tt43*tt44*tt11*tt23*tt18*tt45*tt24+2*tt35*tt21*tt2*&
&tt22*tt43*tt44*tt23*tt18*tt45*tt24
hes(27,1) = tt48
hes(27,2) = 2*tt35*tt2*tt43*tt44*tt26*tt27*tt18*tt45*tt16*tt28+2*&
&tt35*tt2*tt6*tt43*tt44*tt26*tt27*tt18*tt45*tt28+2*tt35*tt21*tt2*t&
&t43*tt44*tt26*tt27*tt18*tt45*tt28
hes(28,1) = tt49
hes(28,2) = 2*tt35*tt2*tt43*tt44*tt11*tt18*tt45*tt30*tt31*tt32+2*&
&tt35*tt2*tt6*tt43*tt44*tt18*tt45*tt30*tt31*tt32+2*tt35*tt21*tt2*t&
&t43*tt44*tt18*tt45*tt30*tt31*tt32
hes(29,1) = tt53
hes(29,2) = 2*tt35*tt2*tt5*tt17*tt10*tt50*tt51*tt52*tt15*tt16+2*t&
&t35*tt2*tt5*tt17*tt10*tt11*tt50*tt51*tt52*tt15+2*tt35*tt2*tt5*tt6&
&*tt17*tt10*tt50*tt51*tt52*tt15
hes(30,1) = tt54
hes(30,2) = 2*tt35*tt2*tt22*tt17*tt23*tt50*tt51*tt52*tt16*tt24+2*&
&tt35*tt2*tt22*tt17*tt11*tt23*tt50*tt51*tt52*tt24+2*tt35*tt21*tt2*&
&tt22*tt17*tt23*tt50*tt51*tt52*tt24
hes(31,1) = tt55
hes(31,2) = 2*tt35*tt2*tt17*tt26*tt27*tt50*tt51*tt52*tt16*tt28+2*&
&tt35*tt2*tt6*tt17*tt26*tt27*tt50*tt51*tt52*tt28+2*tt35*tt21*tt2*t&
&t17*tt26*tt27*tt50*tt51*tt52*tt28
hes(32,1) = tt56
hes(32,2) = 2*tt35*tt2*tt17*tt11*tt50*tt51*tt52*tt30*tt31*tt32+2*&
&tt35*tt2*tt6*tt17*tt50*tt51*tt52*tt30*tt31*tt32+2*tt35*tt21*tt2*t&
&t17*tt50*tt51*tt52*tt30*tt31*tt32
END 
SUBROUTINE quad4_shape_func_jac(val, eps) 
IMPLICIT NONE 
REAL(KIND=8) val(1, 8) 
REAL(KIND=8) eps(2, 1) 
REAL(KIND=8)  tt1 
REAL(KIND=8)  tt2 
REAL(KIND=8)  tt3 
REAL(KIND=8)  tt4 
tt1 = eps(2,1)-1
tt2 = eps(2,1)+1
tt3 = eps(1,1)-1
tt4 = eps(1,1)+1
val(1,1) = tt1/4.0E+0
val(1,2) = -tt2/4.0E+0
val(1,3) = -tt1/4.0E+0
val(1,4) = tt2/4.0E+0
val(1,5) = tt3/4.0E+0
val(1,6) = -tt3/4.0E+0
val(1,7) = -tt4/4.0E+0
val(1,8) = tt4/4.0E+0
END 
SUBROUTINE quad9_shape_func_jac(val, eps) 
IMPLICIT NONE 
REAL(KIND=8) val(1, 18) 
REAL(KIND=8) eps(2, 1) 
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
tt1 = eps(1,1)-1
tt2 = eps(2,1)-1
tt3 = (eps(1,1)*tt2*eps(2,1))/4.0E+0
tt4 = eps(2,1)+1
tt5 = -(eps(1,1)*tt2*tt4)/2.0E+0
tt6 = (eps(1,1)*eps(2,1)*tt4)/4.0E+0
tt7 = eps(1,1)+1
tt8 = (tt1*eps(1,1)*eps(2,1))/4.0E+0
tt9 = -(tt1*tt7*eps(2,1))/2.0E+0
tt10 = (eps(1,1)*tt7*eps(2,1))/4.0E+0
val(1,1) = tt3+(tt1*tt2*eps(2,1))/4.0E+0
val(1,2) = tt5-(tt1*tt2*tt4)/2.0E+0
val(1,3) = tt6+(tt1*eps(2,1)*tt4)/4.0E+0
val(1,4) = (-(tt7*tt2*eps(2,1))/2.0E+0)-(tt1*tt2*eps(2,1))/2.0E+0&
&
val(1,5) = tt7*tt2*tt4+tt1*tt2*tt4
val(1,6) = (-(tt7*eps(2,1)*tt4)/2.0E+0)-(tt1*eps(2,1)*tt4)/2.0E+0&
&
val(1,7) = (tt7*tt2*eps(2,1))/4.0E+0+tt3
val(1,8) = tt5-(tt7*tt2*tt4)/2.0E+0
val(1,9) = (tt7*eps(2,1)*tt4)/4.0E+0+tt6
val(1,10) = tt8+(tt1*eps(1,1)*tt2)/4.0E+0
val(1,11) = (-(tt1*eps(1,1)*tt4)/2.0E+0)-(tt1*eps(1,1)*tt2)/2.0E+&
&0
val(1,12) = (tt1*eps(1,1)*tt4)/4.0E+0+tt8
val(1,13) = tt9-(tt1*tt7*tt2)/2.0E+0
val(1,14) = tt1*tt7*tt4+tt1*tt7*tt2
val(1,15) = tt9-(tt1*tt7*tt4)/2.0E+0
val(1,16) = tt10+(eps(1,1)*tt7*tt2)/4.0E+0
val(1,17) = (-(eps(1,1)*tt7*tt4)/2.0E+0)-(eps(1,1)*tt7*tt2)/2.0E+&
&0
val(1,18) = (eps(1,1)*tt7*tt4)/4.0E+0+tt10
END 
SUBROUTINE quad16_shape_func_jac(val, eps) 
IMPLICIT NONE 
REAL(KIND=8) val(1, 32) 
REAL(KIND=8) eps(2, 1) 
REAL(KIND=8)  tt1 
REAL(KIND=8)  tt2 
REAL(KIND=8)  tt3 
REAL(KIND=8)  tt4 
REAL(KIND=8)  tt5 
REAL(KIND=8)  tt6 
REAL(KIND=8)  tt7 
REAL(KIND=8)  tt8 
tt1 = eps(1,1)-1
tt2 = eps(1,1)+(-1.0E+0)/3.0E+0
tt3 = eps(2,1)-1
tt4 = eps(2,1)+(-1.0E+0)/3.0E+0
tt5 = eps(2,1)+1.0E+0/3.0E+0
tt6 = eps(1,1)+1.0E+0/3.0E+0
tt7 = eps(2,1)+1
tt8 = eps(1,1)+1
val(1,1) = (8.1E+1*tt2*tt6*tt3*tt4*tt5)/2.56E+2+(8.1E+1*tt1*tt6*t&
&t3*tt4*tt5)/2.56E+2+(8.1E+1*tt1*tt2*tt3*tt4*tt5)/2.56E+2
val(1,2) = ((-2.43E+2)*tt2*tt6*tt3*tt4*tt7)/2.56E+2+((-2.43E+2)*t&
&t1*tt6*tt3*tt4*tt7)/2.56E+2+((-2.43E+2)*tt1*tt2*tt3*tt4*tt7)/2.56&
&E+2
val(1,3) = (2.43E+2*tt2*tt6*tt3*tt5*tt7)/2.56E+2+(2.43E+2*tt1*tt6&
&*tt3*tt5*tt7)/2.56E+2+(2.43E+2*tt1*tt2*tt3*tt5*tt7)/2.56E+2
val(1,4) = ((-8.1E+1)*tt2*tt6*tt4*tt5*tt7)/2.56E+2+((-8.1E+1)*tt1&
&*tt6*tt4*tt5*tt7)/2.56E+2+((-8.1E+1)*tt1*tt2*tt4*tt5*tt7)/2.56E+2
val(1,5) = ((-2.43E+2)*tt2*tt8*tt3*tt4*tt5)/2.56E+2+((-2.43E+2)*t&
&t1*tt8*tt3*tt4*tt5)/2.56E+2+((-2.43E+2)*tt1*tt2*tt3*tt4*tt5)/2.56&
&E+2
val(1,6) = (7.29E+2*tt2*tt8*tt3*tt4*tt7)/2.56E+2+(7.29E+2*tt1*tt8&
&*tt3*tt4*tt7)/2.56E+2+(7.29E+2*tt1*tt2*tt3*tt4*tt7)/2.56E+2
val(1,7) = ((-7.29E+2)*tt2*tt8*tt3*tt5*tt7)/2.56E+2+((-7.29E+2)*t&
&t1*tt8*tt3*tt5*tt7)/2.56E+2+((-7.29E+2)*tt1*tt2*tt3*tt5*tt7)/2.56&
&E+2
val(1,8) = (2.43E+2*tt2*tt8*tt4*tt5*tt7)/2.56E+2+(2.43E+2*tt1*tt8&
&*tt4*tt5*tt7)/2.56E+2+(2.43E+2*tt1*tt2*tt4*tt5*tt7)/2.56E+2
val(1,9) = (2.43E+2*tt6*tt8*tt3*tt4*tt5)/2.56E+2+(2.43E+2*tt1*tt8&
&*tt3*tt4*tt5)/2.56E+2+(2.43E+2*tt1*tt6*tt3*tt4*tt5)/2.56E+2
val(1,10) = ((-7.29E+2)*tt6*tt8*tt3*tt4*tt7)/2.56E+2+((-7.29E+2)*&
&tt1*tt8*tt3*tt4*tt7)/2.56E+2+((-7.29E+2)*tt1*tt6*tt3*tt4*tt7)/2.5&
&6E+2
val(1,11) = (7.29E+2*tt6*tt8*tt3*tt5*tt7)/2.56E+2+(7.29E+2*tt1*tt&
&8*tt3*tt5*tt7)/2.56E+2+(7.29E+2*tt1*tt6*tt3*tt5*tt7)/2.56E+2
val(1,12) = ((-2.43E+2)*tt6*tt8*tt4*tt5*tt7)/2.56E+2+((-2.43E+2)*&
&tt1*tt8*tt4*tt5*tt7)/2.56E+2+((-2.43E+2)*tt1*tt6*tt4*tt5*tt7)/2.5&
&6E+2
val(1,13) = ((-8.1E+1)*tt6*tt8*tt3*tt4*tt5)/2.56E+2+((-8.1E+1)*tt&
&2*tt8*tt3*tt4*tt5)/2.56E+2+((-8.1E+1)*tt2*tt6*tt3*tt4*tt5)/2.56E+&
&2
val(1,14) = (2.43E+2*tt6*tt8*tt3*tt4*tt7)/2.56E+2+(2.43E+2*tt2*tt&
&8*tt3*tt4*tt7)/2.56E+2+(2.43E+2*tt2*tt6*tt3*tt4*tt7)/2.56E+2
val(1,15) = ((-2.43E+2)*tt6*tt8*tt3*tt5*tt7)/2.56E+2+((-2.43E+2)*&
&tt2*tt8*tt3*tt5*tt7)/2.56E+2+((-2.43E+2)*tt2*tt6*tt3*tt5*tt7)/2.5&
&6E+2
val(1,16) = (8.1E+1*tt6*tt8*tt4*tt5*tt7)/2.56E+2+(8.1E+1*tt2*tt8*&
&tt4*tt5*tt7)/2.56E+2+(8.1E+1*tt2*tt6*tt4*tt5*tt7)/2.56E+2
val(1,17) = (8.1E+1*tt1*tt2*tt6*tt4*tt5)/2.56E+2+(8.1E+1*tt1*tt2*&
&tt6*tt3*tt5)/2.56E+2+(8.1E+1*tt1*tt2*tt6*tt3*tt4)/2.56E+2
val(1,18) = ((-2.43E+2)*tt1*tt2*tt6*tt4*tt7)/2.56E+2+((-2.43E+2)*&
&tt1*tt2*tt6*tt3*tt7)/2.56E+2+((-2.43E+2)*tt1*tt2*tt6*tt3*tt4)/2.5&
&6E+2
val(1,19) = (2.43E+2*tt1*tt2*tt6*tt5*tt7)/2.56E+2+(2.43E+2*tt1*tt&
&2*tt6*tt3*tt7)/2.56E+2+(2.43E+2*tt1*tt2*tt6*tt3*tt5)/2.56E+2
val(1,20) = ((-8.1E+1)*tt1*tt2*tt6*tt5*tt7)/2.56E+2+((-8.1E+1)*tt&
&1*tt2*tt6*tt4*tt7)/2.56E+2+((-8.1E+1)*tt1*tt2*tt6*tt4*tt5)/2.56E+&
&2
val(1,21) = ((-2.43E+2)*tt1*tt2*tt8*tt4*tt5)/2.56E+2+((-2.43E+2)*&
&tt1*tt2*tt8*tt3*tt5)/2.56E+2+((-2.43E+2)*tt1*tt2*tt8*tt3*tt4)/2.5&
&6E+2
val(1,22) = (7.29E+2*tt1*tt2*tt8*tt4*tt7)/2.56E+2+(7.29E+2*tt1*tt&
&2*tt8*tt3*tt7)/2.56E+2+(7.29E+2*tt1*tt2*tt8*tt3*tt4)/2.56E+2
val(1,23) = ((-7.29E+2)*tt1*tt2*tt8*tt5*tt7)/2.56E+2+((-7.29E+2)*&
&tt1*tt2*tt8*tt3*tt7)/2.56E+2+((-7.29E+2)*tt1*tt2*tt8*tt3*tt5)/2.5&
&6E+2
val(1,24) = (2.43E+2*tt1*tt2*tt8*tt5*tt7)/2.56E+2+(2.43E+2*tt1*tt&
&2*tt8*tt4*tt7)/2.56E+2+(2.43E+2*tt1*tt2*tt8*tt4*tt5)/2.56E+2
val(1,25) = (2.43E+2*tt1*tt6*tt8*tt4*tt5)/2.56E+2+(2.43E+2*tt1*tt&
&6*tt8*tt3*tt5)/2.56E+2+(2.43E+2*tt1*tt6*tt8*tt3*tt4)/2.56E+2
val(1,26) = ((-7.29E+2)*tt1*tt6*tt8*tt4*tt7)/2.56E+2+((-7.29E+2)*&
&tt1*tt6*tt8*tt3*tt7)/2.56E+2+((-7.29E+2)*tt1*tt6*tt8*tt3*tt4)/2.5&
&6E+2
val(1,27) = (7.29E+2*tt1*tt6*tt8*tt5*tt7)/2.56E+2+(7.29E+2*tt1*tt&
&6*tt8*tt3*tt7)/2.56E+2+(7.29E+2*tt1*tt6*tt8*tt3*tt5)/2.56E+2
val(1,28) = ((-2.43E+2)*tt1*tt6*tt8*tt5*tt7)/2.56E+2+((-2.43E+2)*&
&tt1*tt6*tt8*tt4*tt7)/2.56E+2+((-2.43E+2)*tt1*tt6*tt8*tt4*tt5)/2.5&
&6E+2
val(1,29) = ((-8.1E+1)*tt2*tt6*tt8*tt4*tt5)/2.56E+2+((-8.1E+1)*tt&
&2*tt6*tt8*tt3*tt5)/2.56E+2+((-8.1E+1)*tt2*tt6*tt8*tt3*tt4)/2.56E+&
&2
val(1,30) = (2.43E+2*tt2*tt6*tt8*tt4*tt7)/2.56E+2+(2.43E+2*tt2*tt&
&6*tt8*tt3*tt7)/2.56E+2+(2.43E+2*tt2*tt6*tt8*tt3*tt4)/2.56E+2
val(1,31) = ((-2.43E+2)*tt2*tt6*tt8*tt5*tt7)/2.56E+2+((-2.43E+2)*&
&tt2*tt6*tt8*tt3*tt7)/2.56E+2+((-2.43E+2)*tt2*tt6*tt8*tt3*tt5)/2.5&
&6E+2
val(1,32) = (8.1E+1*tt2*tt6*tt8*tt5*tt7)/2.56E+2+(8.1E+1*tt2*tt6*&
&tt8*tt4*tt7)/2.56E+2+(8.1E+1*tt2*tt6*tt8*tt4*tt5)/2.56E+2
END 
SUBROUTINE quad4_shape_func_val(val, eps) 
IMPLICIT NONE 
REAL(KIND=8) val(4, 1) 
REAL(KIND=8) eps(2, 1) 
REAL(KIND=8)  tt1 
REAL(KIND=8)  tt2 
REAL(KIND=8)  tt3 
REAL(KIND=8)  tt4 
tt1 = eps(1,1)-1
tt2 = eps(2,1)-1
tt3 = eps(2,1)+1
tt4 = eps(1,1)+1
val(1,1) = (tt1*tt2)/4.0E+0
val(2,1) = -(tt1*tt3)/4.0E+0
val(3,1) = -(tt4*tt2)/4.0E+0
val(4,1) = (tt4*tt3)/4.0E+0
END 
SUBROUTINE quad9_shape_func_val(val, eps) 
IMPLICIT NONE 
REAL(KIND=8) val(9, 1) 
REAL(KIND=8) eps(2, 1) 
REAL(KIND=8)  tt1 
REAL(KIND=8)  tt2 
REAL(KIND=8)  tt3 
REAL(KIND=8)  tt4 
tt1 = eps(1,1)-1
tt2 = eps(2,1)-1
tt3 = eps(2,1)+1
tt4 = eps(1,1)+1
val(1,1) = (tt1*eps(1,1)*tt2*eps(2,1))/4.0E+0
val(2,1) = -(tt1*eps(1,1)*tt2*tt3)/2.0E+0
val(3,1) = (tt1*eps(1,1)*eps(2,1)*tt3)/4.0E+0
val(4,1) = -(tt1*tt4*tt2*eps(2,1))/2.0E+0
val(5,1) = tt1*tt4*tt2*tt3
val(6,1) = -(tt1*tt4*eps(2,1)*tt3)/2.0E+0
val(7,1) = (eps(1,1)*tt4*tt2*eps(2,1))/4.0E+0
val(8,1) = -(eps(1,1)*tt4*tt2*tt3)/2.0E+0
val(9,1) = (eps(1,1)*tt4*eps(2,1)*tt3)/4.0E+0
END 
SUBROUTINE quad16_shape_func_val(val, eps) 
IMPLICIT NONE 
REAL(KIND=8) val(16, 1) 
REAL(KIND=8) eps(2, 1) 
REAL(KIND=8)  tt1 
REAL(KIND=8)  tt2 
REAL(KIND=8)  tt3 
REAL(KIND=8)  tt4 
REAL(KIND=8)  tt5 
REAL(KIND=8)  tt6 
REAL(KIND=8)  tt7 
REAL(KIND=8)  tt8 
tt1 = eps(1,1)-1
tt2 = eps(1,1)+(-1.0E+0)/3.0E+0
tt3 = eps(1,1)+1.0E+0/3.0E+0
tt4 = eps(2,1)-1
tt5 = eps(2,1)+(-1.0E+0)/3.0E+0
tt6 = eps(2,1)+1.0E+0/3.0E+0
tt7 = eps(2,1)+1
tt8 = eps(1,1)+1
val(1,1) = (8.1E+1*tt1*tt2*tt3*tt4*tt5*tt6)/2.56E+2
val(2,1) = ((-2.43E+2)*tt1*tt2*tt3*tt4*tt5*tt7)/2.56E+2
val(3,1) = (2.43E+2*tt1*tt2*tt3*tt4*tt6*tt7)/2.56E+2
val(4,1) = ((-8.1E+1)*tt1*tt2*tt3*tt5*tt6*tt7)/2.56E+2
val(5,1) = ((-2.43E+2)*tt1*tt2*tt8*tt4*tt5*tt6)/2.56E+2
val(6,1) = (7.29E+2*tt1*tt2*tt8*tt4*tt5*tt7)/2.56E+2
val(7,1) = ((-7.29E+2)*tt1*tt2*tt8*tt4*tt6*tt7)/2.56E+2
val(8,1) = (2.43E+2*tt1*tt2*tt8*tt5*tt6*tt7)/2.56E+2
val(9,1) = (2.43E+2*tt1*tt3*tt8*tt4*tt5*tt6)/2.56E+2
val(10,1) = ((-7.29E+2)*tt1*tt3*tt8*tt4*tt5*tt7)/2.56E+2
val(11,1) = (7.29E+2*tt1*tt3*tt8*tt4*tt6*tt7)/2.56E+2
val(12,1) = ((-2.43E+2)*tt1*tt3*tt8*tt5*tt6*tt7)/2.56E+2
val(13,1) = ((-8.1E+1)*tt2*tt3*tt8*tt4*tt5*tt6)/2.56E+2
val(14,1) = (2.43E+2*tt2*tt3*tt8*tt4*tt5*tt7)/2.56E+2
val(15,1) = ((-2.43E+2)*tt2*tt3*tt8*tt4*tt6*tt7)/2.56E+2
val(16,1) = (8.1E+1*tt2*tt3*tt8*tt5*tt6*tt7)/2.56E+2
END 
SUBROUTINE hex8_shape_func(val, eps) 
IMPLICIT NONE 
REAL(KIND=8) val(8, 1) 
REAL(KIND=8) eps(3, 1) 
REAL(KIND=8)  tt1 
REAL(KIND=8)  tt2 
REAL(KIND=8)  tt3 
REAL(KIND=8)  tt4 
REAL(KIND=8)  tt5 
REAL(KIND=8)  tt6 
tt1 = eps(1,1)-1
tt2 = eps(2,1)-1
tt3 = eps(3,1)-1
tt4 = eps(3,1)+1
tt5 = eps(2,1)+1
tt6 = eps(1,1)+1
val(1,1) = -(tt1*tt2*tt3)/8.0E+0
val(2,1) = (tt1*tt2*tt4)/8.0E+0
val(3,1) = (tt1*tt5*tt3)/8.0E+0
val(4,1) = -(tt1*tt5*tt4)/8.0E+0
val(5,1) = (tt6*tt2*tt3)/8.0E+0
val(6,1) = -(tt6*tt2*tt4)/8.0E+0
val(7,1) = -(tt6*tt5*tt3)/8.0E+0
val(8,1) = (tt6*tt5*tt4)/8.0E+0
END 
SUBROUTINE hex8_shape_func_jac(jac, eps) 
IMPLICIT NONE 
REAL(KIND=8) jac(8, 3) 
REAL(KIND=8) eps(3, 1) 
REAL(KIND=8)  tt1 
REAL(KIND=8)  tt2 
REAL(KIND=8)  tt3 
REAL(KIND=8)  tt4 
REAL(KIND=8)  tt5 
REAL(KIND=8)  tt6 
tt1 = eps(2,1)-1
tt2 = eps(3,1)-1
tt3 = eps(1,1)-1
tt4 = eps(3,1)+1
tt5 = eps(2,1)+1
tt6 = eps(1,1)+1
jac(1,1) = -(tt1*tt2)/8.0E+0
jac(1,2) = -(tt3*tt2)/8.0E+0
jac(1,3) = -(tt3*tt1)/8.0E+0
jac(2,1) = (tt1*tt4)/8.0E+0
jac(2,2) = (tt3*tt4)/8.0E+0
jac(2,3) = (tt3*tt1)/8.0E+0
jac(3,1) = (tt5*tt2)/8.0E+0
jac(3,2) = (tt3*tt2)/8.0E+0
jac(3,3) = (tt3*tt5)/8.0E+0
jac(4,1) = -(tt5*tt4)/8.0E+0
jac(4,2) = -(tt3*tt4)/8.0E+0
jac(4,3) = -(tt3*tt5)/8.0E+0
jac(5,1) = (tt1*tt2)/8.0E+0
jac(5,2) = (tt6*tt2)/8.0E+0
jac(5,3) = (tt6*tt1)/8.0E+0
jac(6,1) = -(tt1*tt4)/8.0E+0
jac(6,2) = -(tt6*tt4)/8.0E+0
jac(6,3) = -(tt6*tt1)/8.0E+0
jac(7,1) = -(tt5*tt2)/8.0E+0
jac(7,2) = -(tt6*tt2)/8.0E+0
jac(7,3) = -(tt6*tt5)/8.0E+0
jac(8,1) = (tt5*tt4)/8.0E+0
jac(8,2) = (tt6*tt4)/8.0E+0
jac(8,3) = (tt6*tt5)/8.0E+0
END 
SUBROUTINE hex27_shape_func(val, eps) 
IMPLICIT NONE 
REAL(KIND=8) val(27, 1) 
REAL(KIND=8) eps(3, 1) 
REAL(KIND=8)  tt1 
REAL(KIND=8)  tt2 
REAL(KIND=8)  tt3 
REAL(KIND=8)  tt4 
REAL(KIND=8)  tt5 
REAL(KIND=8)  tt6 
tt1 = eps(1,1)-1
tt2 = eps(2,1)-1
tt3 = eps(3,1)-1
tt4 = eps(3,1)+1
tt5 = eps(2,1)+1
tt6 = eps(1,1)+1
val(1,1) = (tt1*eps(1,1)*tt2*eps(2,1)*tt3*eps(3,1))/8.0E+0
val(2,1) = -(tt1*eps(1,1)*tt2*eps(2,1)*tt3*tt4)/4.0E+0
val(3,1) = (tt1*eps(1,1)*tt2*eps(2,1)*eps(3,1)*tt4)/8.0E+0
val(4,1) = -(tt1*eps(1,1)*tt2*tt5*tt3*eps(3,1))/4.0E+0
val(5,1) = (tt1*eps(1,1)*tt2*tt5*tt3*tt4)/2.0E+0
val(6,1) = -(tt1*eps(1,1)*tt2*tt5*eps(3,1)*tt4)/4.0E+0
val(7,1) = (tt1*eps(1,1)*eps(2,1)*tt5*tt3*eps(3,1))/8.0E+0
val(8,1) = -(tt1*eps(1,1)*eps(2,1)*tt5*tt3*tt4)/4.0E+0
val(9,1) = (tt1*eps(1,1)*eps(2,1)*tt5*eps(3,1)*tt4)/8.0E+0
val(10,1) = -(tt1*tt6*tt2*eps(2,1)*tt3*eps(3,1))/4.0E+0
val(11,1) = (tt1*tt6*tt2*eps(2,1)*tt3*tt4)/2.0E+0
val(12,1) = -(tt1*tt6*tt2*eps(2,1)*eps(3,1)*tt4)/4.0E+0
val(13,1) = (tt1*tt6*tt2*tt5*tt3*eps(3,1))/2.0E+0
val(14,1) = -tt1*tt6*tt2*tt5*tt3*tt4
val(15,1) = (tt1*tt6*tt2*tt5*eps(3,1)*tt4)/2.0E+0
val(16,1) = -(tt1*tt6*eps(2,1)*tt5*tt3*eps(3,1))/4.0E+0
val(17,1) = (tt1*tt6*eps(2,1)*tt5*tt3*tt4)/2.0E+0
val(18,1) = -(tt1*tt6*eps(2,1)*tt5*eps(3,1)*tt4)/4.0E+0
val(19,1) = (eps(1,1)*tt6*tt2*eps(2,1)*tt3*eps(3,1))/8.0E+0
val(20,1) = -(eps(1,1)*tt6*tt2*eps(2,1)*tt3*tt4)/4.0E+0
val(21,1) = (eps(1,1)*tt6*tt2*eps(2,1)*eps(3,1)*tt4)/8.0E+0
val(22,1) = -(eps(1,1)*tt6*tt2*tt5*tt3*eps(3,1))/4.0E+0
val(23,1) = (eps(1,1)*tt6*tt2*tt5*tt3*tt4)/2.0E+0
val(24,1) = -(eps(1,1)*tt6*tt2*tt5*eps(3,1)*tt4)/4.0E+0
val(25,1) = (eps(1,1)*tt6*eps(2,1)*tt5*tt3*eps(3,1))/8.0E+0
val(26,1) = -(eps(1,1)*tt6*eps(2,1)*tt5*tt3*tt4)/4.0E+0
val(27,1) = (eps(1,1)*tt6*eps(2,1)*tt5*eps(3,1)*tt4)/8.0E+0
END 
SUBROUTINE hex27_shape_func_jac(jac, eps) 
IMPLICIT NONE 
REAL(KIND=8) jac(27, 3) 
REAL(KIND=8) eps(3, 1) 
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
tt1 = eps(1,1)-1
tt2 = eps(2,1)-1
tt3 = eps(3,1)-1
tt4 = (eps(1,1)*tt2*eps(2,1)*tt3*eps(3,1))/8.0E+0
tt5 = (tt1*eps(1,1)*eps(2,1)*tt3*eps(3,1))/8.0E+0
tt6 = (tt1*eps(1,1)*tt2*eps(2,1)*eps(3,1))/8.0E+0
tt7 = eps(3,1)+1
tt8 = -(eps(1,1)*tt2*eps(2,1)*tt3*tt7)/4.0E+0
tt9 = -(tt1*eps(1,1)*eps(2,1)*tt3*tt7)/4.0E+0
tt10 = (eps(1,1)*tt2*eps(2,1)*eps(3,1)*tt7)/8.0E+0
tt11 = (tt1*eps(1,1)*eps(2,1)*eps(3,1)*tt7)/8.0E+0
tt12 = eps(2,1)+1
tt13 = -(eps(1,1)*tt2*tt12*tt3*eps(3,1))/4.0E+0
tt14 = -(tt1*eps(1,1)*tt2*tt12*eps(3,1))/4.0E+0
tt15 = (eps(1,1)*tt2*tt12*tt3*tt7)/2.0E+0
tt16 = -(eps(1,1)*tt2*tt12*eps(3,1)*tt7)/4.0E+0
tt17 = (eps(1,1)*eps(2,1)*tt12*tt3*eps(3,1))/8.0E+0
tt18 = (tt1*eps(1,1)*eps(2,1)*tt12*eps(3,1))/8.0E+0
tt19 = -(eps(1,1)*eps(2,1)*tt12*tt3*tt7)/4.0E+0
tt20 = (eps(1,1)*eps(2,1)*tt12*eps(3,1)*tt7)/8.0E+0
tt21 = eps(1,1)+1
tt22 = -(tt1*tt21*eps(2,1)*tt3*eps(3,1))/4.0E+0
tt23 = -(tt1*tt21*tt2*eps(2,1)*eps(3,1))/4.0E+0
tt24 = (tt1*tt21*eps(2,1)*tt3*tt7)/2.0E+0
tt25 = -(tt1*tt21*eps(2,1)*eps(3,1)*tt7)/4.0E+0
tt26 = (tt1*tt21*tt2*tt12*eps(3,1))/2.0E+0
tt27 = -(tt1*tt21*eps(2,1)*tt12*eps(3,1))/4.0E+0
tt28 = (eps(1,1)*tt21*eps(2,1)*tt3*eps(3,1))/8.0E+0
tt29 = (eps(1,1)*tt21*tt2*eps(2,1)*eps(3,1))/8.0E+0
tt30 = -(eps(1,1)*tt21*eps(2,1)*tt3*tt7)/4.0E+0
tt31 = (eps(1,1)*tt21*eps(2,1)*eps(3,1)*tt7)/8.0E+0
tt32 = -(eps(1,1)*tt21*tt2*tt12*eps(3,1))/4.0E+0
tt33 = (eps(1,1)*tt21*eps(2,1)*tt12*eps(3,1))/8.0E+0
jac(1,1) = tt4+(tt1*tt2*eps(2,1)*tt3*eps(3,1))/8.0E+0
jac(1,2) = tt5+(tt1*eps(1,1)*tt2*tt3*eps(3,1))/8.0E+0
jac(1,3) = tt6+(tt1*eps(1,1)*tt2*eps(2,1)*tt3)/8.0E+0
jac(2,1) = tt8-(tt1*tt2*eps(2,1)*tt3*tt7)/4.0E+0
jac(2,2) = tt9-(tt1*eps(1,1)*tt2*tt3*tt7)/4.0E+0
jac(2,3) = (-(tt1*eps(1,1)*tt2*eps(2,1)*tt7)/4.0E+0)-(tt1*eps(1,1&
&)*tt2*eps(2,1)*tt3)/4.0E+0
jac(3,1) = tt10+(tt1*tt2*eps(2,1)*eps(3,1)*tt7)/8.0E+0
jac(3,2) = tt11+(tt1*eps(1,1)*tt2*eps(3,1)*tt7)/8.0E+0
jac(3,3) = (tt1*eps(1,1)*tt2*eps(2,1)*tt7)/8.0E+0+tt6
jac(4,1) = tt13-(tt1*tt2*tt12*tt3*eps(3,1))/4.0E+0
jac(4,2) = (-(tt1*eps(1,1)*tt12*tt3*eps(3,1))/4.0E+0)-(tt1*eps(1,&
&1)*tt2*tt3*eps(3,1))/4.0E+0
jac(4,3) = tt14-(tt1*eps(1,1)*tt2*tt12*tt3)/4.0E+0
jac(5,1) = tt15+(tt1*tt2*tt12*tt3*tt7)/2.0E+0
jac(5,2) = (tt1*eps(1,1)*tt12*tt3*tt7)/2.0E+0+(tt1*eps(1,1)*tt2*t&
&t3*tt7)/2.0E+0
jac(5,3) = (tt1*eps(1,1)*tt2*tt12*tt7)/2.0E+0+(tt1*eps(1,1)*tt2*t&
&t12*tt3)/2.0E+0
jac(6,1) = tt16-(tt1*tt2*tt12*eps(3,1)*tt7)/4.0E+0
jac(6,2) = (-(tt1*eps(1,1)*tt12*eps(3,1)*tt7)/4.0E+0)-(tt1*eps(1,&
&1)*tt2*eps(3,1)*tt7)/4.0E+0
jac(6,3) = tt14-(tt1*eps(1,1)*tt2*tt12*tt7)/4.0E+0
jac(7,1) = tt17+(tt1*eps(2,1)*tt12*tt3*eps(3,1))/8.0E+0
jac(7,2) = (tt1*eps(1,1)*tt12*tt3*eps(3,1))/8.0E+0+tt5
jac(7,3) = tt18+(tt1*eps(1,1)*eps(2,1)*tt12*tt3)/8.0E+0
jac(8,1) = tt19-(tt1*eps(2,1)*tt12*tt3*tt7)/4.0E+0
jac(8,2) = tt9-(tt1*eps(1,1)*tt12*tt3*tt7)/4.0E+0
jac(8,3) = (-(tt1*eps(1,1)*eps(2,1)*tt12*tt7)/4.0E+0)-(tt1*eps(1,&
&1)*eps(2,1)*tt12*tt3)/4.0E+0
jac(9,1) = tt20+(tt1*eps(2,1)*tt12*eps(3,1)*tt7)/8.0E+0
jac(9,2) = (tt1*eps(1,1)*tt12*eps(3,1)*tt7)/8.0E+0+tt11
jac(9,3) = (tt1*eps(1,1)*eps(2,1)*tt12*tt7)/8.0E+0+tt18
jac(10,1) = (-(tt21*tt2*eps(2,1)*tt3*eps(3,1))/4.0E+0)-(tt1*tt2*e&
&ps(2,1)*tt3*eps(3,1))/4.0E+0
jac(10,2) = tt22-(tt1*tt21*tt2*tt3*eps(3,1))/4.0E+0
jac(10,3) = tt23-(tt1*tt21*tt2*eps(2,1)*tt3)/4.0E+0
jac(11,1) = (tt21*tt2*eps(2,1)*tt3*tt7)/2.0E+0+(tt1*tt2*eps(2,1)*&
&tt3*tt7)/2.0E+0
jac(11,2) = tt24+(tt1*tt21*tt2*tt3*tt7)/2.0E+0
jac(11,3) = (tt1*tt21*tt2*eps(2,1)*tt7)/2.0E+0+(tt1*tt21*tt2*eps(&
&2,1)*tt3)/2.0E+0
jac(12,1) = (-(tt21*tt2*eps(2,1)*eps(3,1)*tt7)/4.0E+0)-(tt1*tt2*e&
&ps(2,1)*eps(3,1)*tt7)/4.0E+0
jac(12,2) = tt25-(tt1*tt21*tt2*eps(3,1)*tt7)/4.0E+0
jac(12,3) = tt23-(tt1*tt21*tt2*eps(2,1)*tt7)/4.0E+0
jac(13,1) = (tt21*tt2*tt12*tt3*eps(3,1))/2.0E+0+(tt1*tt2*tt12*tt3&
&*eps(3,1))/2.0E+0
jac(13,2) = (tt1*tt21*tt12*tt3*eps(3,1))/2.0E+0+(tt1*tt21*tt2*tt3&
&*eps(3,1))/2.0E+0
jac(13,3) = tt26+(tt1*tt21*tt2*tt12*tt3)/2.0E+0
jac(14,1) = (-tt21*tt2*tt12*tt3*tt7)-tt1*tt2*tt12*tt3*tt7
jac(14,2) = (-tt1*tt21*tt12*tt3*tt7)-tt1*tt21*tt2*tt3*tt7
jac(14,3) = (-tt1*tt21*tt2*tt12*tt7)-tt1*tt21*tt2*tt12*tt3
jac(15,1) = (tt21*tt2*tt12*eps(3,1)*tt7)/2.0E+0+(tt1*tt2*tt12*eps&
&(3,1)*tt7)/2.0E+0
jac(15,2) = (tt1*tt21*tt12*eps(3,1)*tt7)/2.0E+0+(tt1*tt21*tt2*eps&
&(3,1)*tt7)/2.0E+0
jac(15,3) = (tt1*tt21*tt2*tt12*tt7)/2.0E+0+tt26
jac(16,1) = (-(tt21*eps(2,1)*tt12*tt3*eps(3,1))/4.0E+0)-(tt1*eps(&
&2,1)*tt12*tt3*eps(3,1))/4.0E+0
jac(16,2) = tt22-(tt1*tt21*tt12*tt3*eps(3,1))/4.0E+0
jac(16,3) = tt27-(tt1*tt21*eps(2,1)*tt12*tt3)/4.0E+0
jac(17,1) = (tt21*eps(2,1)*tt12*tt3*tt7)/2.0E+0+(tt1*eps(2,1)*tt1&
&2*tt3*tt7)/2.0E+0
jac(17,2) = (tt1*tt21*tt12*tt3*tt7)/2.0E+0+tt24
jac(17,3) = (tt1*tt21*eps(2,1)*tt12*tt7)/2.0E+0+(tt1*tt21*eps(2,1&
&)*tt12*tt3)/2.0E+0
jac(18,1) = (-(tt21*eps(2,1)*tt12*eps(3,1)*tt7)/4.0E+0)-(tt1*eps(&
&2,1)*tt12*eps(3,1)*tt7)/4.0E+0
jac(18,2) = tt25-(tt1*tt21*tt12*eps(3,1)*tt7)/4.0E+0
jac(18,3) = tt27-(tt1*tt21*eps(2,1)*tt12*tt7)/4.0E+0
jac(19,1) = (tt21*tt2*eps(2,1)*tt3*eps(3,1))/8.0E+0+tt4
jac(19,2) = tt28+(eps(1,1)*tt21*tt2*tt3*eps(3,1))/8.0E+0
jac(19,3) = tt29+(eps(1,1)*tt21*tt2*eps(2,1)*tt3)/8.0E+0
jac(20,1) = tt8-(tt21*tt2*eps(2,1)*tt3*tt7)/4.0E+0
jac(20,2) = tt30-(eps(1,1)*tt21*tt2*tt3*tt7)/4.0E+0
jac(20,3) = (-(eps(1,1)*tt21*tt2*eps(2,1)*tt7)/4.0E+0)-(eps(1,1)*&
&tt21*tt2*eps(2,1)*tt3)/4.0E+0
jac(21,1) = (tt21*tt2*eps(2,1)*eps(3,1)*tt7)/8.0E+0+tt10
jac(21,2) = tt31+(eps(1,1)*tt21*tt2*eps(3,1)*tt7)/8.0E+0
jac(21,3) = (eps(1,1)*tt21*tt2*eps(2,1)*tt7)/8.0E+0+tt29
jac(22,1) = tt13-(tt21*tt2*tt12*tt3*eps(3,1))/4.0E+0
jac(22,2) = (-(eps(1,1)*tt21*tt12*tt3*eps(3,1))/4.0E+0)-(eps(1,1)&
&*tt21*tt2*tt3*eps(3,1))/4.0E+0
jac(22,3) = tt32-(eps(1,1)*tt21*tt2*tt12*tt3)/4.0E+0
jac(23,1) = (tt21*tt2*tt12*tt3*tt7)/2.0E+0+tt15
jac(23,2) = (eps(1,1)*tt21*tt12*tt3*tt7)/2.0E+0+(eps(1,1)*tt21*tt&
&2*tt3*tt7)/2.0E+0
jac(23,3) = (eps(1,1)*tt21*tt2*tt12*tt7)/2.0E+0+(eps(1,1)*tt21*tt&
&2*tt12*tt3)/2.0E+0
jac(24,1) = tt16-(tt21*tt2*tt12*eps(3,1)*tt7)/4.0E+0
jac(24,2) = (-(eps(1,1)*tt21*tt12*eps(3,1)*tt7)/4.0E+0)-(eps(1,1)&
&*tt21*tt2*eps(3,1)*tt7)/4.0E+0
jac(24,3) = tt32-(eps(1,1)*tt21*tt2*tt12*tt7)/4.0E+0
jac(25,1) = (tt21*eps(2,1)*tt12*tt3*eps(3,1))/8.0E+0+tt17
jac(25,2) = (eps(1,1)*tt21*tt12*tt3*eps(3,1))/8.0E+0+tt28
jac(25,3) = tt33+(eps(1,1)*tt21*eps(2,1)*tt12*tt3)/8.0E+0
jac(26,1) = tt19-(tt21*eps(2,1)*tt12*tt3*tt7)/4.0E+0
jac(26,2) = tt30-(eps(1,1)*tt21*tt12*tt3*tt7)/4.0E+0
jac(26,3) = (-(eps(1,1)*tt21*eps(2,1)*tt12*tt7)/4.0E+0)-(eps(1,1)&
&*tt21*eps(2,1)*tt12*tt3)/4.0E+0
jac(27,1) = (tt21*eps(2,1)*tt12*eps(3,1)*tt7)/8.0E+0+tt20
jac(27,2) = (eps(1,1)*tt21*tt12*eps(3,1)*tt7)/8.0E+0+tt31
jac(27,3) = (eps(1,1)*tt21*eps(2,1)*tt12*tt7)/8.0E+0+tt33
END 
SUBROUTINE quad4_rbf_shape_func_jac(val, eps) 
IMPLICIT NONE 
REAL(KIND=8) val(1, 8) 
REAL(KIND=8) eps(2, 1) 
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
tt1 = eps(1,1)+1
tt2 = -tt1**2
tt3 = eps(2,1)+1
tt4 = -tt3**2
tt5 = exp(tt4+tt2)
tt6 = eps(1,1)-1
tt7 = -tt6**2
tt8 = eps(2,1)-1
tt9 = -tt8**2
tt10 = exp(tt9+tt7)
tt11 = exp(tt9+tt2)
tt12 = exp(tt4+tt7)
tt13 = tt5+tt12+tt11+tt10
tt14 = 1/tt13
tt15 = 1/tt13**2
tt16 = (-2*tt5*tt1)-2*tt11*tt1-2*tt12*tt6-2*tt10*tt6
tt17 = (-2*tt5*tt3)-2*tt12*tt3-2*tt11*tt8-2*tt10*tt8
val(1,1) = (-tt5*tt15*tt16)-2*tt5*tt14*tt1
val(1,2) = (-tt11*tt15*tt16)-2*tt11*tt14*tt1
val(1,3) = (-tt12*tt15*tt16)-2*tt12*tt14*tt6
val(1,4) = (-tt10*tt15*tt16)-2*tt10*tt14*tt6
val(1,5) = (-tt5*tt15*tt17)-2*tt5*tt14*tt3
val(1,6) = (-tt11*tt15*tt17)-2*tt11*tt14*tt8
val(1,7) = (-tt12*tt15*tt17)-2*tt12*tt14*tt3
val(1,8) = (-tt10*tt15*tt17)-2*tt10*tt14*tt8
END 
SUBROUTINE quad9_rbf_shape_func_jac(val, eps) 
IMPLICIT NONE 
REAL(KIND=8) val(1, 18) 
REAL(KIND=8) eps(2, 1) 
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
tt1 = eps(1,1)+1
tt2 = -tt1**2
tt3 = eps(2,1)+1
tt4 = -tt3**2
tt5 = exp(tt4+tt2)
tt6 = eps(1,1)-1
tt7 = -tt6**2
tt8 = eps(2,1)-1
tt9 = -tt8**2
tt10 = exp(tt9+tt7)
tt11 = -eps(1,1)**2
tt12 = exp(tt9+tt11)
tt13 = exp(tt9+tt2)
tt14 = -eps(2,1)**2
tt15 = exp(tt14+tt7)
tt16 = exp(tt14+tt11)
tt17 = exp(tt14+tt2)
tt18 = exp(tt4+tt7)
tt19 = exp(tt4+tt11)
tt20 = tt5+tt19+tt18+tt17+tt16+tt15+tt13+tt12+tt10
tt21 = 1/tt20
tt22 = 1/tt20**2
tt23 = (-2*tt5*tt1)-2*tt17*tt1-2*tt13*tt1-2*tt19*eps(1,1)-2*tt16*&
&eps(1,1)-2*tt12*eps(1,1)-2*tt18*tt6-2*tt15*tt6-2*tt10*tt6
tt24 = (-2*tt5*tt3)-2*tt19*tt3-2*tt18*tt3-2*tt17*eps(2,1)-2*tt16*&
&eps(2,1)-2*tt15*eps(2,1)-2*tt13*tt8-2*tt12*tt8-2*tt10*tt8
val(1,1) = (-tt5*tt22*tt23)-2*tt5*tt21*tt1
val(1,2) = (-tt17*tt22*tt23)-2*tt17*tt21*tt1
val(1,3) = (-tt13*tt22*tt23)-2*tt13*tt21*tt1
val(1,4) = (-tt19*tt22*tt23)-2*tt19*tt21*eps(1,1)
val(1,5) = (-tt16*tt22*tt23)-2*tt16*tt21*eps(1,1)
val(1,6) = (-tt12*tt22*tt23)-2*tt12*tt21*eps(1,1)
val(1,7) = (-tt18*tt22*tt23)-2*tt18*tt21*tt6
val(1,8) = (-tt15*tt22*tt23)-2*tt15*tt21*tt6
val(1,9) = (-tt10*tt22*tt23)-2*tt10*tt21*tt6
val(1,10) = (-tt5*tt22*tt24)-2*tt5*tt21*tt3
val(1,11) = (-tt17*tt22*tt24)-2*tt17*tt21*eps(2,1)
val(1,12) = (-tt13*tt22*tt24)-2*tt13*tt21*tt8
val(1,13) = (-tt19*tt22*tt24)-2*tt19*tt21*tt3
val(1,14) = (-tt16*tt22*tt24)-2*tt16*tt21*eps(2,1)
val(1,15) = (-tt12*tt22*tt24)-2*tt12*tt21*tt8
val(1,16) = (-tt18*tt22*tt24)-2*tt18*tt21*tt3
val(1,17) = (-tt15*tt22*tt24)-2*tt15*tt21*eps(2,1)
val(1,18) = (-tt10*tt22*tt24)-2*tt10*tt21*tt8
END 
SUBROUTINE quad16_rbf_shape_func_jac(val, eps) 
IMPLICIT NONE 
REAL(KIND=8) val(1, 32) 
REAL(KIND=8) eps(2, 1) 
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
tt1 = eps(1,1)+1
tt2 = -tt1**2
tt3 = eps(2,1)+1
tt4 = -tt3**2
tt5 = exp(tt4+tt2)
tt6 = eps(1,1)-1
tt7 = -tt6**2
tt8 = eps(2,1)-1
tt9 = -tt8**2
tt10 = exp(tt9+tt7)
tt11 = eps(1,1)+(-1.0E+0)/3.0E+0
tt12 = -tt11**2
tt13 = exp(tt9+tt12)
tt14 = eps(1,1)+1.0E+0/3.0E+0
tt15 = -tt14**2
tt16 = exp(tt9+tt15)
tt17 = exp(tt9+tt2)
tt18 = eps(2,1)+(-1.0E+0)/3.0E+0
tt19 = -tt18**2
tt20 = exp(tt19+tt7)
tt21 = exp(tt19+tt12)
tt22 = exp(tt19+tt15)
tt23 = exp(tt19+tt2)
tt24 = eps(2,1)+1.0E+0/3.0E+0
tt25 = -tt24**2
tt26 = exp(tt25+tt7)
tt27 = exp(tt25+tt12)
tt28 = exp(tt25+tt15)
tt29 = exp(tt25+tt2)
tt30 = exp(tt4+tt7)
tt31 = exp(tt4+tt12)
tt32 = exp(tt4+tt15)
tt33 = tt5+tt32+tt31+tt30+tt29+tt28+tt27+tt26+tt23+tt22+tt21+tt20&
&+tt17+tt16+tt13+tt10
tt34 = 1/tt33
tt35 = 1/tt33**2
tt36 = (-2*tt5*tt1)-2*tt29*tt1-2*tt23*tt1-2*tt17*tt1-2*tt32*tt14-&
&2*tt28*tt14-2*tt22*tt14-2*tt16*tt14-2*tt31*tt11-2*tt27*tt11-2*tt2&
&1*tt11-2*tt13*tt11-2*tt30*tt6-2*tt26*tt6-2*tt20*tt6-2*tt10*tt6
tt37 = (-2*tt5*tt3)-2*tt32*tt3-2*tt31*tt3-2*tt30*tt3-2*tt29*tt24-&
&2*tt28*tt24-2*tt27*tt24-2*tt26*tt24-2*tt23*tt18-2*tt22*tt18-2*tt2&
&1*tt18-2*tt20*tt18-2*tt17*tt8-2*tt16*tt8-2*tt13*tt8-2*tt10*tt8
val(1,1) = (-tt5*tt35*tt36)-2*tt5*tt34*tt1
val(1,2) = (-tt29*tt35*tt36)-2*tt29*tt34*tt1
val(1,3) = (-tt23*tt35*tt36)-2*tt23*tt34*tt1
val(1,4) = (-tt17*tt35*tt36)-2*tt17*tt34*tt1
val(1,5) = (-tt32*tt35*tt36)-2*tt32*tt34*tt14
val(1,6) = (-tt28*tt35*tt36)-2*tt28*tt34*tt14
val(1,7) = (-tt22*tt35*tt36)-2*tt22*tt34*tt14
val(1,8) = (-tt16*tt35*tt36)-2*tt16*tt34*tt14
val(1,9) = (-tt31*tt35*tt36)-2*tt31*tt34*tt11
val(1,10) = (-tt27*tt35*tt36)-2*tt27*tt34*tt11
val(1,11) = (-tt21*tt35*tt36)-2*tt21*tt34*tt11
val(1,12) = (-tt13*tt35*tt36)-2*tt13*tt34*tt11
val(1,13) = (-tt30*tt35*tt36)-2*tt30*tt34*tt6
val(1,14) = (-tt26*tt35*tt36)-2*tt26*tt34*tt6
val(1,15) = (-tt20*tt35*tt36)-2*tt20*tt34*tt6
val(1,16) = (-tt10*tt35*tt36)-2*tt10*tt34*tt6
val(1,17) = (-tt5*tt35*tt37)-2*tt5*tt34*tt3
val(1,18) = (-tt29*tt35*tt37)-2*tt29*tt34*tt24
val(1,19) = (-tt23*tt35*tt37)-2*tt23*tt34*tt18
val(1,20) = (-tt17*tt35*tt37)-2*tt17*tt34*tt8
val(1,21) = (-tt32*tt35*tt37)-2*tt32*tt34*tt3
val(1,22) = (-tt28*tt35*tt37)-2*tt28*tt34*tt24
val(1,23) = (-tt22*tt35*tt37)-2*tt22*tt34*tt18
val(1,24) = (-tt16*tt35*tt37)-2*tt16*tt34*tt8
val(1,25) = (-tt31*tt35*tt37)-2*tt31*tt34*tt3
val(1,26) = (-tt27*tt35*tt37)-2*tt27*tt34*tt24
val(1,27) = (-tt21*tt35*tt37)-2*tt21*tt34*tt18
val(1,28) = (-tt13*tt35*tt37)-2*tt13*tt34*tt8
val(1,29) = (-tt30*tt35*tt37)-2*tt30*tt34*tt3
val(1,30) = (-tt26*tt35*tt37)-2*tt26*tt34*tt24
val(1,31) = (-tt20*tt35*tt37)-2*tt20*tt34*tt18
val(1,32) = (-tt10*tt35*tt37)-2*tt10*tt34*tt8
END 
SUBROUTINE quad4_rbf_shape_func_val(val, eps) 
IMPLICIT NONE 
REAL(KIND=8) val(4, 1) 
REAL(KIND=8) eps(2, 1) 
REAL(KIND=8)  tt1 
REAL(KIND=8)  tt2 
REAL(KIND=8)  tt3 
REAL(KIND=8)  tt4 
REAL(KIND=8)  tt5 
REAL(KIND=8)  tt6 
REAL(KIND=8)  tt7 
REAL(KIND=8)  tt8 
REAL(KIND=8)  tt9 
tt1 = -(eps(1,1)+1)**2
tt2 = -(eps(2,1)+1)**2
tt3 = exp(tt2+tt1)
tt4 = -(eps(1,1)-1)**2
tt5 = -(eps(2,1)-1)**2
tt6 = exp(tt5+tt4)
tt7 = exp(tt5+tt1)
tt8 = exp(tt2+tt4)
tt9 = 1/(tt3+tt8+tt7+tt6)
val(1,1) = tt3*tt9
val(2,1) = tt7*tt9
val(3,1) = tt8*tt9
val(4,1) = tt6*tt9
END 
SUBROUTINE quad9_rbf_shape_func_val(val, eps) 
IMPLICIT NONE 
REAL(KIND=8) val(9, 1) 
REAL(KIND=8) eps(2, 1) 
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
tt1 = -(eps(1,1)+1)**2
tt2 = -(eps(2,1)+1)**2
tt3 = exp(tt2+tt1)
tt4 = -(eps(1,1)-1)**2
tt5 = -(eps(2,1)-1)**2
tt6 = exp(tt5+tt4)
tt7 = -eps(1,1)**2
tt8 = exp(tt5+tt7)
tt9 = exp(tt5+tt1)
tt10 = -eps(2,1)**2
tt11 = exp(tt10+tt4)
tt12 = exp(tt10+tt7)
tt13 = exp(tt10+tt1)
tt14 = exp(tt2+tt4)
tt15 = exp(tt2+tt7)
tt16 = 1/(tt3+tt15+tt14+tt13+tt12+tt11+tt9+tt8+tt6)
val(1,1) = tt3*tt16
val(2,1) = tt13*tt16
val(3,1) = tt9*tt16
val(4,1) = tt15*tt16
val(5,1) = tt12*tt16
val(6,1) = tt8*tt16
val(7,1) = tt14*tt16
val(8,1) = tt11*tt16
val(9,1) = tt6*tt16
END 
SUBROUTINE quad16_rbf_shape_func_val(val, eps) 
IMPLICIT NONE 
REAL(KIND=8) val(16, 1) 
REAL(KIND=8) eps(2, 1) 
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
tt1 = -(eps(1,1)+1)**2
tt2 = -(eps(2,1)+1)**2
tt3 = exp(tt2+tt1)
tt4 = -(eps(1,1)-1)**2
tt5 = -(eps(2,1)-1)**2
tt6 = exp(tt5+tt4)
tt7 = -(eps(1,1)+(-1.0E+0)/3.0E+0)**2
tt8 = exp(tt5+tt7)
tt9 = -(eps(1,1)+1.0E+0/3.0E+0)**2
tt10 = exp(tt5+tt9)
tt11 = exp(tt5+tt1)
tt12 = -(eps(2,1)+(-1.0E+0)/3.0E+0)**2
tt13 = exp(tt12+tt4)
tt14 = exp(tt12+tt7)
tt15 = exp(tt12+tt9)
tt16 = exp(tt12+tt1)
tt17 = -(eps(2,1)+1.0E+0/3.0E+0)**2
tt18 = exp(tt17+tt4)
tt19 = exp(tt17+tt7)
tt20 = exp(tt17+tt9)
tt21 = exp(tt17+tt1)
tt22 = exp(tt2+tt4)
tt23 = exp(tt2+tt7)
tt24 = exp(tt2+tt9)
tt25 = 1/(tt3+tt24+tt23+tt22+tt21+tt20+tt19+tt18+tt16+tt15+tt14+t&
&t13+tt11+tt10+tt8+tt6)
val(1,1) = tt3*tt25
val(2,1) = tt21*tt25
val(3,1) = tt16*tt25
val(4,1) = tt11*tt25
val(5,1) = tt24*tt25
val(6,1) = tt20*tt25
val(7,1) = tt15*tt25
val(8,1) = tt10*tt25
val(9,1) = tt23*tt25
val(10,1) = tt19*tt25
val(11,1) = tt14*tt25
val(12,1) = tt8*tt25
val(13,1) = tt22*tt25
val(14,1) = tt18*tt25
val(15,1) = tt13*tt25
val(16,1) = tt6*tt25
END 
