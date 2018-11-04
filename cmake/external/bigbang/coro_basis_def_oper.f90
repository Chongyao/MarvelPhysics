SUBROUTINE coro_ani_basis_df_2d(val, X, X0, R, Nu, Nx, dNude, dNxde, Dm) 
IMPLICIT NONE 
REAL(KIND=8) val(1, 4) 
REAL(KIND=8) X(2, 1) 
REAL(KIND=8) X0(2, 1) 
REAL(KIND=8) R(2, 2) 
REAL(KIND=8) Nu(2, 1) 
REAL(KIND=8) Nx(1, 1) 
REAL(KIND=8) dNude(2, 2) 
REAL(KIND=8) dNxde(1, 2) 
REAL(KIND=8) Dm(2, 2) 
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
tt1 = -X0(1,1)
tt2 = R(2,1)*X(2,1)+tt1+R(1,1)*X(1,1)
tt3 = R(1,2)*X0(2,1)+R(1,1)*X0(1,1)+tt1
tt4 = -X0(2,1)
tt5 = X(2,1)*R(2,2)+tt4+X(1,1)*R(1,2)
tt6 = R(1,2)*dNude(2,1)*tt5+dNxde(1,1)*tt3+R(1,1)*dNude(1,1)*tt2
tt7 = R(1,2)*tt5*dNude(2,2)+dNxde(1,2)*tt3+R(1,1)*dNude(1,2)*tt2
tt8 = X0(2,1)*R(2,2)+tt4+X0(1,1)*R(2,1)
tt9 = dNxde(1,1)*tt8+dNude(2,1)*R(2,2)*tt5+dNude(1,1)*R(2,1)*tt2
tt10 = R(2,2)*tt5*dNude(2,2)+dNxde(1,2)*tt8+dNude(1,2)*R(2,1)*tt2&
&
val(1,1) = Dm(2,1)*tt7+Dm(1,1)*tt6
val(1,2) = Dm(2,1)*tt10+Dm(1,1)*tt9
val(1,3) = Dm(2,2)*tt7+Dm(1,2)*tt6
val(1,4) = Dm(2,2)*tt10+Dm(1,2)*tt9
END 
SUBROUTINE coro_ani_basis_df_2d_jac(jac, X, X0, R, Nu, Nx, dNude, dNxde, Dm) 
IMPLICIT NONE 
REAL(KIND=8) jac(4, 2) 
REAL(KIND=8) X(2, 1) 
REAL(KIND=8) X0(2, 1) 
REAL(KIND=8) R(2, 2) 
REAL(KIND=8) Nu(2, 1) 
REAL(KIND=8) Nx(1, 1) 
REAL(KIND=8) dNude(2, 2) 
REAL(KIND=8) dNxde(1, 2) 
REAL(KIND=8) Dm(2, 2) 
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
tt1 = R(1,1)**2
tt2 = R(1,2)**2
tt3 = tt2*dNude(2,1)+tt1*dNude(1,1)
tt4 = tt2*dNude(2,2)+tt1*dNude(1,2)
tt5 = R(1,2)*dNude(2,1)*R(2,2)+R(1,1)*dNude(1,1)*R(2,1)
tt6 = R(1,2)*R(2,2)*dNude(2,2)+R(1,1)*dNude(1,2)*R(2,1)
tt7 = Dm(2,1)*tt6+Dm(1,1)*tt5
tt8 = R(2,1)**2
tt9 = R(2,2)**2
tt10 = dNude(2,1)*tt9+dNude(1,1)*tt8
tt11 = tt9*dNude(2,2)+dNude(1,2)*tt8
tt12 = Dm(2,2)*tt6+Dm(1,2)*tt5
jac(1,1) = Dm(2,1)*tt4+Dm(1,1)*tt3
jac(1,2) = tt7
jac(2,1) = tt7
jac(2,2) = Dm(2,1)*tt11+Dm(1,1)*tt10
jac(3,1) = Dm(2,2)*tt4+Dm(1,2)*tt3
jac(3,2) = tt12
jac(4,1) = tt12
jac(4,2) = Dm(2,2)*tt11+Dm(1,2)*tt10
END 
SUBROUTINE coro_aug_ani_df_2d(val, X, X0, R, Nu, Nx, dNude, dNxde, Dm) 
IMPLICIT NONE 
REAL(KIND=8) val(1, 4) 
REAL(KIND=8) X(2, 1) 
REAL(KIND=8) X0(2, 1) 
REAL(KIND=8) R(2, 2) 
REAL(KIND=8) Nu(2, 2) 
REAL(KIND=8) Nx(1, 1) 
REAL(KIND=8) dNude(4, 2) 
REAL(KIND=8) dNxde(1, 2) 
REAL(KIND=8) Dm(2, 2) 
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
tt1 = -X0(1,1)
tt2 = R(2,1)*X(2,1)+tt1+R(1,1)*X(1,1)
tt3 = R(1,2)*X0(2,1)+R(1,1)*X0(1,1)+tt1
tt4 = -X0(2,1)
tt5 = X(2,1)*R(2,2)+tt4+X(1,1)*R(1,2)
tt6 = R(1,2)*tt5*dNude(4,1)+R(1,1)*tt5*dNude(3,1)+R(1,2)*tt2*dNud&
&e(2,1)+dNxde(1,1)*tt3+R(1,1)*dNude(1,1)*tt2
tt7 = R(1,2)*tt5*dNude(4,2)+R(1,1)*tt5*dNude(3,2)+R(1,2)*tt2*dNud&
&e(2,2)+dNxde(1,2)*tt3+R(1,1)*dNude(1,2)*tt2
tt8 = X0(2,1)*R(2,2)+tt4+X0(1,1)*R(2,1)
tt9 = R(2,2)*tt5*dNude(4,1)+R(2,1)*tt5*dNude(3,1)+dNxde(1,1)*tt8+&
&tt2*dNude(2,1)*R(2,2)+dNude(1,1)*R(2,1)*tt2
tt10 = R(2,2)*tt5*dNude(4,2)+R(2,1)*tt5*dNude(3,2)+tt2*R(2,2)*dNu&
&de(2,2)+dNxde(1,2)*tt8+dNude(1,2)*R(2,1)*tt2
val(1,1) = Dm(2,1)*tt7+Dm(1,1)*tt6
val(1,2) = Dm(2,1)*tt10+Dm(1,1)*tt9
val(1,3) = Dm(2,2)*tt7+Dm(1,2)*tt6
val(1,4) = Dm(2,2)*tt10+Dm(1,2)*tt9
END 
SUBROUTINE coro_aug_ani_df_2d_jac(jac, X, X0, R, Nu, Nx, dNude, dNxde, Dm) 
IMPLICIT NONE 
REAL(KIND=8) jac(4, 2) 
REAL(KIND=8) X(2, 1) 
REAL(KIND=8) X0(2, 1) 
REAL(KIND=8) R(2, 2) 
REAL(KIND=8) Nu(2, 2) 
REAL(KIND=8) Nx(1, 1) 
REAL(KIND=8) dNude(4, 2) 
REAL(KIND=8) dNxde(1, 2) 
REAL(KIND=8) Dm(2, 2) 
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
tt1 = R(1,1)**2
tt2 = R(1,2)**2
tt3 = tt2*dNude(4,1)+R(1,1)*R(1,2)*dNude(3,1)+R(1,1)*R(1,2)*dNude&
&(2,1)+tt1*dNude(1,1)
tt4 = tt2*dNude(4,2)+R(1,1)*R(1,2)*dNude(3,2)+R(1,1)*R(1,2)*dNude&
&(2,2)+tt1*dNude(1,2)
tt5 = R(1,1)*dNude(1,1)*R(2,1)
tt6 = R(1,2)*R(2,2)*dNude(4,1)
tt7 = tt6+R(1,1)*R(2,2)*dNude(3,1)+R(1,2)*R(2,1)*dNude(2,1)+tt5
tt8 = R(1,1)*dNude(1,2)*R(2,1)
tt9 = R(1,2)*R(2,2)*dNude(4,2)
tt10 = tt9+R(1,1)*R(2,2)*dNude(3,2)+R(1,2)*R(2,1)*dNude(2,2)+tt8
tt11 = tt6+R(1,2)*R(2,1)*dNude(3,1)+R(1,1)*dNude(2,1)*R(2,2)+tt5
tt12 = tt9+R(1,2)*R(2,1)*dNude(3,2)+R(1,1)*R(2,2)*dNude(2,2)+tt8
tt13 = R(2,1)**2
tt14 = R(2,2)**2
tt15 = tt14*dNude(4,1)+R(2,1)*R(2,2)*dNude(3,1)+R(2,1)*dNude(2,1)&
&*R(2,2)+dNude(1,1)*tt13
tt16 = tt14*dNude(4,2)+R(2,1)*R(2,2)*dNude(3,2)+R(2,1)*R(2,2)*dNu&
&de(2,2)+dNude(1,2)*tt13
jac(1,1) = Dm(2,1)*tt4+Dm(1,1)*tt3
jac(1,2) = Dm(2,1)*tt10+Dm(1,1)*tt7
jac(2,1) = Dm(2,1)*tt12+Dm(1,1)*tt11
jac(2,2) = Dm(2,1)*tt16+Dm(1,1)*tt15
jac(3,1) = Dm(2,2)*tt4+Dm(1,2)*tt3
jac(3,2) = Dm(2,2)*tt10+Dm(1,2)*tt7
jac(4,1) = Dm(2,2)*tt12+Dm(1,2)*tt11
jac(4,2) = Dm(2,2)*tt16+Dm(1,2)*tt15
END 
SUBROUTINE coro_iso_basis_df_2d(val, X, X0, R, Nu, Nx, dNude, dNxde, Dm) 
IMPLICIT NONE 
REAL(KIND=8) val(1, 4) 
REAL(KIND=8) X(2, 1) 
REAL(KIND=8) X0(2, 1) 
REAL(KIND=8) R(2, 2) 
REAL(KIND=8) Nu(1, 1) 
REAL(KIND=8) Nx(1, 1) 
REAL(KIND=8) dNude(1, 2) 
REAL(KIND=8) dNxde(1, 2) 
REAL(KIND=8) Dm(2, 2) 
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
tt1 = -X0(1,1)
tt2 = R(1,2)*X0(2,1)+R(1,1)*X0(1,1)+tt1
tt3 = R(2,1)*X(2,1)+tt1+R(1,1)*X(1,1)
tt4 = -X0(2,1)
tt5 = X(2,1)*R(2,2)+tt4+X(1,1)*R(1,2)
tt6 = R(1,2)*tt5+R(1,1)*tt3
tt7 = dNude(1,1)*tt6+dNxde(1,1)*tt2
tt8 = dNude(1,2)*tt6+dNxde(1,2)*tt2
tt9 = X0(2,1)*R(2,2)+tt4+X0(1,1)*R(2,1)
tt10 = R(2,2)*tt5+R(2,1)*tt3
tt11 = dNude(1,1)*tt10+dNxde(1,1)*tt9
tt12 = dNude(1,2)*tt10+dNxde(1,2)*tt9
val(1,1) = Dm(2,1)*tt8+Dm(1,1)*tt7
val(1,2) = Dm(2,1)*tt12+Dm(1,1)*tt11
val(1,3) = Dm(2,2)*tt8+Dm(1,2)*tt7
val(1,4) = Dm(2,2)*tt12+Dm(1,2)*tt11
END 
SUBROUTINE coro_iso_basis_df_2d_jac(jac, X, X0, R, Nu, Nx, dNude, dNxde, Dm) 
IMPLICIT NONE 
REAL(KIND=8) jac(4, 2) 
REAL(KIND=8) X(2, 1) 
REAL(KIND=8) X0(2, 1) 
REAL(KIND=8) R(2, 2) 
REAL(KIND=8) Nu(1, 1) 
REAL(KIND=8) Nx(1, 1) 
REAL(KIND=8) dNude(1, 2) 
REAL(KIND=8) dNxde(1, 2) 
REAL(KIND=8) Dm(2, 2) 
REAL(KIND=8)  tt1 
REAL(KIND=8)  tt2 
REAL(KIND=8)  tt3 
REAL(KIND=8)  tt4 
REAL(KIND=8)  tt5 
tt1 = R(1,2)**2+R(1,1)**2
tt2 = R(1,2)*R(2,2)+R(1,1)*R(2,1)
tt3 = dNude(1,2)*Dm(2,1)*tt2+Dm(1,1)*dNude(1,1)*tt2
tt4 = R(2,2)**2+R(2,1)**2
tt5 = dNude(1,2)*Dm(2,2)*tt2+dNude(1,1)*Dm(1,2)*tt2
jac(1,1) = tt1*dNude(1,2)*Dm(2,1)+Dm(1,1)*dNude(1,1)*tt1
jac(1,2) = tt3
jac(2,1) = tt3
jac(2,2) = dNude(1,2)*Dm(2,1)*tt4+Dm(1,1)*dNude(1,1)*tt4
jac(3,1) = tt1*dNude(1,2)*Dm(2,2)+dNude(1,1)*Dm(1,2)*tt1
jac(3,2) = tt5
jac(4,1) = tt5
jac(4,2) = dNude(1,2)*Dm(2,2)*tt4+dNude(1,1)*Dm(1,2)*tt4
END 
SUBROUTINE coro_ani_basis_df_3d(val, X, X0, R, Nu, Nx, dNude, dNxde, Dm) 
IMPLICIT NONE 
REAL(KIND=8) val(1, 9) 
REAL(KIND=8) X(3, 1) 
REAL(KIND=8) X0(3, 1) 
REAL(KIND=8) R(3, 3) 
REAL(KIND=8) Nu(3, 1) 
REAL(KIND=8) Nx(1, 1) 
REAL(KIND=8) dNude(3, 3) 
REAL(KIND=8) dNxde(1, 3) 
REAL(KIND=8) Dm(3, 3) 
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
tt1 = -X0(1,1)
tt2 = R(3,1)*X(3,1)+R(2,1)*X(2,1)+tt1+R(1,1)*X(1,1)
tt3 = R(1,3)*X0(3,1)+R(1,2)*X0(2,1)+R(1,1)*X0(1,1)+tt1
tt4 = -X0(2,1)
tt5 = X(3,1)*R(3,2)+X(2,1)*R(2,2)+tt4+X(1,1)*R(1,2)
tt6 = -X0(3,1)
tt7 = X(3,1)*R(3,3)+tt6+X(2,1)*R(2,3)+X(1,1)*R(1,3)
tt8 = R(1,3)*dNude(3,1)*tt7+R(1,2)*dNude(2,1)*tt5+dNxde(1,1)*tt3+&
&R(1,1)*dNude(1,1)*tt2
tt9 = R(1,3)*dNude(3,2)*tt7+R(1,2)*dNude(2,2)*tt5+dNxde(1,2)*tt3+&
&R(1,1)*dNude(1,2)*tt2
tt10 = R(1,3)*tt7*dNude(3,3)+R(1,2)*dNude(2,3)*tt5+dNxde(1,3)*tt3&
&+R(1,1)*dNude(1,3)*tt2
tt11 = R(2,3)*X0(3,1)+X0(2,1)*R(2,2)+tt4+X0(1,1)*R(2,1)
tt12 = R(2,3)*dNude(3,1)*tt7+dNude(2,1)*R(2,2)*tt5+dNxde(1,1)*tt1&
&1+dNude(1,1)*R(2,1)*tt2
tt13 = R(2,3)*dNude(3,2)*tt7+R(2,2)*dNude(2,2)*tt5+dNxde(1,2)*tt1&
&1+dNude(1,2)*R(2,1)*tt2
tt14 = R(2,3)*tt7*dNude(3,3)+R(2,2)*dNude(2,3)*tt5+dNxde(1,3)*tt1&
&1+dNude(1,3)*R(2,1)*tt2
tt15 = X0(3,1)*R(3,3)+X0(2,1)*R(3,2)+tt6+X0(1,1)*R(3,1)
tt16 = dNxde(1,1)*tt15+dNude(3,1)*R(3,3)*tt7+dNude(2,1)*R(3,2)*tt&
&5+dNude(1,1)*R(3,1)*tt2
tt17 = dNxde(1,2)*tt15+dNude(3,2)*R(3,3)*tt7+dNude(2,2)*R(3,2)*tt&
&5+dNude(1,2)*R(3,1)*tt2
tt18 = R(3,3)*tt7*dNude(3,3)+dNxde(1,3)*tt15+dNude(2,3)*R(3,2)*tt&
&5+dNude(1,3)*R(3,1)*tt2
val(1,1) = Dm(3,1)*tt10+Dm(2,1)*tt9+Dm(1,1)*tt8
val(1,2) = Dm(3,1)*tt14+Dm(2,1)*tt13+Dm(1,1)*tt12
val(1,3) = Dm(3,1)*tt18+Dm(2,1)*tt17+Dm(1,1)*tt16
val(1,4) = Dm(3,2)*tt10+Dm(2,2)*tt9+Dm(1,2)*tt8
val(1,5) = Dm(3,2)*tt14+Dm(2,2)*tt13+Dm(1,2)*tt12
val(1,6) = Dm(3,2)*tt18+Dm(2,2)*tt17+Dm(1,2)*tt16
val(1,7) = Dm(3,3)*tt10+Dm(2,3)*tt9+Dm(1,3)*tt8
val(1,8) = Dm(3,3)*tt14+Dm(2,3)*tt13+Dm(1,3)*tt12
val(1,9) = Dm(3,3)*tt18+Dm(2,3)*tt17+Dm(1,3)*tt16
END 
SUBROUTINE coro_ani_basis_df_3d_jac(jac, X, X0, R, Nu, Nx, dNude, dNxde, Dm) 
IMPLICIT NONE 
REAL(KIND=8) jac(9, 3) 
REAL(KIND=8) X(3, 1) 
REAL(KIND=8) X0(3, 1) 
REAL(KIND=8) R(3, 3) 
REAL(KIND=8) Nu(3, 1) 
REAL(KIND=8) Nx(1, 1) 
REAL(KIND=8) dNude(3, 3) 
REAL(KIND=8) dNxde(1, 3) 
REAL(KIND=8) Dm(3, 3) 
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
tt1 = R(1,1)**2
tt2 = R(1,2)**2
tt3 = R(1,3)**2
tt4 = tt3*dNude(3,1)+tt2*dNude(2,1)+tt1*dNude(1,1)
tt5 = tt3*dNude(3,2)+tt2*dNude(2,2)+tt1*dNude(1,2)
tt6 = tt3*dNude(3,3)+tt2*dNude(2,3)+tt1*dNude(1,3)
tt7 = R(1,3)*R(2,3)*dNude(3,1)+R(1,2)*dNude(2,1)*R(2,2)+R(1,1)*dN&
&ude(1,1)*R(2,1)
tt8 = R(1,3)*R(2,3)*dNude(3,2)+R(1,2)*R(2,2)*dNude(2,2)+R(1,1)*dN&
&ude(1,2)*R(2,1)
tt9 = R(1,3)*R(2,3)*dNude(3,3)+R(1,2)*R(2,2)*dNude(2,3)+R(1,1)*dN&
&ude(1,3)*R(2,1)
tt10 = Dm(3,1)*tt9+Dm(2,1)*tt8+Dm(1,1)*tt7
tt11 = R(1,3)*dNude(3,1)*R(3,3)+R(1,2)*dNude(2,1)*R(3,2)+R(1,1)*d&
&Nude(1,1)*R(3,1)
tt12 = R(1,3)*dNude(3,2)*R(3,3)+R(1,2)*dNude(2,2)*R(3,2)+R(1,1)*d&
&Nude(1,2)*R(3,1)
tt13 = R(1,3)*R(3,3)*dNude(3,3)+R(1,2)*dNude(2,3)*R(3,2)+R(1,1)*d&
&Nude(1,3)*R(3,1)
tt14 = Dm(3,1)*tt13+Dm(2,1)*tt12+Dm(1,1)*tt11
tt15 = R(2,1)**2
tt16 = R(2,2)**2
tt17 = R(2,3)**2
tt18 = tt17*dNude(3,1)+dNude(2,1)*tt16+dNude(1,1)*tt15
tt19 = tt17*dNude(3,2)+tt16*dNude(2,2)+dNude(1,2)*tt15
tt20 = tt17*dNude(3,3)+tt16*dNude(2,3)+dNude(1,3)*tt15
tt21 = R(2,3)*dNude(3,1)*R(3,3)+dNude(2,1)*R(2,2)*R(3,2)+dNude(1,&
&1)*R(2,1)*R(3,1)
tt22 = R(2,3)*dNude(3,2)*R(3,3)+R(2,2)*dNude(2,2)*R(3,2)+dNude(1,&
&2)*R(2,1)*R(3,1)
tt23 = R(2,3)*R(3,3)*dNude(3,3)+R(2,2)*dNude(2,3)*R(3,2)+dNude(1,&
&3)*R(2,1)*R(3,1)
tt24 = Dm(3,1)*tt23+Dm(2,1)*tt22+Dm(1,1)*tt21
tt25 = R(3,1)**2
tt26 = R(3,2)**2
tt27 = R(3,3)**2
tt28 = dNude(3,1)*tt27+dNude(2,1)*tt26+dNude(1,1)*tt25
tt29 = dNude(3,2)*tt27+dNude(2,2)*tt26+dNude(1,2)*tt25
tt30 = tt27*dNude(3,3)+dNude(2,3)*tt26+dNude(1,3)*tt25
tt31 = Dm(3,2)*tt9+Dm(2,2)*tt8+Dm(1,2)*tt7
tt32 = Dm(3,2)*tt13+Dm(2,2)*tt12+Dm(1,2)*tt11
tt33 = Dm(3,2)*tt23+Dm(2,2)*tt22+Dm(1,2)*tt21
tt34 = Dm(3,3)*tt9+Dm(2,3)*tt8+Dm(1,3)*tt7
tt35 = Dm(3,3)*tt13+Dm(2,3)*tt12+Dm(1,3)*tt11
tt36 = Dm(3,3)*tt23+Dm(2,3)*tt22+Dm(1,3)*tt21
jac(1,1) = Dm(3,1)*tt6+Dm(2,1)*tt5+Dm(1,1)*tt4
jac(1,2) = tt10
jac(1,3) = tt14
jac(2,1) = tt10
jac(2,2) = Dm(3,1)*tt20+Dm(2,1)*tt19+Dm(1,1)*tt18
jac(2,3) = tt24
jac(3,1) = tt14
jac(3,2) = tt24
jac(3,3) = Dm(3,1)*tt30+Dm(2,1)*tt29+Dm(1,1)*tt28
jac(4,1) = Dm(3,2)*tt6+Dm(2,2)*tt5+Dm(1,2)*tt4
jac(4,2) = tt31
jac(4,3) = tt32
jac(5,1) = tt31
jac(5,2) = Dm(3,2)*tt20+Dm(2,2)*tt19+Dm(1,2)*tt18
jac(5,3) = tt33
jac(6,1) = tt32
jac(6,2) = tt33
jac(6,3) = Dm(3,2)*tt30+Dm(2,2)*tt29+Dm(1,2)*tt28
jac(7,1) = Dm(3,3)*tt6+Dm(2,3)*tt5+Dm(1,3)*tt4
jac(7,2) = tt34
jac(7,3) = tt35
jac(8,1) = tt34
jac(8,2) = Dm(3,3)*tt20+Dm(2,3)*tt19+Dm(1,3)*tt18
jac(8,3) = tt36
jac(9,1) = tt35
jac(9,2) = tt36
jac(9,3) = Dm(3,3)*tt30+Dm(2,3)*tt29+Dm(1,3)*tt28
END 
SUBROUTINE coro_aug_ani_df_3d(val, X, X0, R, Nu, Nx, dNude, dNxde, Dm) 
IMPLICIT NONE 
REAL(KIND=8) val(1, 9) 
REAL(KIND=8) X(3, 1) 
REAL(KIND=8) X0(3, 1) 
REAL(KIND=8) R(3, 3) 
REAL(KIND=8) Nu(3, 3) 
REAL(KIND=8) Nx(1, 1) 
REAL(KIND=8) dNude(9, 3) 
REAL(KIND=8) dNxde(1, 3) 
REAL(KIND=8) Dm(3, 3) 
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
tt1 = -X0(1,1)
tt2 = R(3,1)*X(3,1)+R(2,1)*X(2,1)+tt1+R(1,1)*X(1,1)
tt3 = R(1,3)*X0(3,1)+R(1,2)*X0(2,1)+R(1,1)*X0(1,1)+tt1
tt4 = -X0(2,1)
tt5 = X(3,1)*R(3,2)+X(2,1)*R(2,2)+tt4+X(1,1)*R(1,2)
tt6 = -X0(3,1)
tt7 = X(3,1)*R(3,3)+tt6+X(2,1)*R(2,3)+X(1,1)*R(1,3)
tt8 = R(1,3)*tt7*dNude(9,1)+R(1,2)*tt7*dNude(8,1)+R(1,1)*tt7*dNud&
&e(7,1)+R(1,3)*tt5*dNude(6,1)+R(1,2)*tt5*dNude(5,1)+R(1,1)*tt5*dNu&
&de(4,1)+R(1,3)*tt2*dNude(3,1)+dNxde(1,1)*tt3+R(1,2)*dNude(2,1)*tt&
&2+R(1,1)*dNude(1,1)*tt2
tt9 = R(1,3)*tt7*dNude(9,2)+R(1,2)*tt7*dNude(8,2)+R(1,1)*tt7*dNud&
&e(7,2)+R(1,3)*tt5*dNude(6,2)+R(1,2)*tt5*dNude(5,2)+R(1,1)*tt5*dNu&
&de(4,2)+R(1,3)*tt2*dNude(3,2)+dNxde(1,2)*tt3+R(1,2)*dNude(2,2)*tt&
&2+R(1,1)*dNude(1,2)*tt2
tt10 = R(1,3)*tt7*dNude(9,3)+R(1,2)*tt7*dNude(8,3)+R(1,1)*tt7*dNu&
&de(7,3)+R(1,3)*tt5*dNude(6,3)+R(1,2)*tt5*dNude(5,3)+R(1,1)*tt5*dN&
&ude(4,3)+R(1,3)*tt2*dNude(3,3)+dNxde(1,3)*tt3+R(1,2)*dNude(2,3)*t&
&t2+R(1,1)*dNude(1,3)*tt2
tt11 = R(2,3)*X0(3,1)+X0(2,1)*R(2,2)+tt4+X0(1,1)*R(2,1)
tt12 = R(2,3)*tt7*dNude(9,1)+R(2,2)*tt7*dNude(8,1)+R(2,1)*tt7*dNu&
&de(7,1)+R(2,3)*tt5*dNude(6,1)+R(2,2)*tt5*dNude(5,1)+R(2,1)*tt5*dN&
&ude(4,1)+R(2,3)*tt2*dNude(3,1)+dNxde(1,1)*tt11+dNude(2,1)*R(2,2)*&
&tt2+dNude(1,1)*R(2,1)*tt2
tt13 = R(2,3)*tt7*dNude(9,2)+R(2,2)*tt7*dNude(8,2)+R(2,1)*tt7*dNu&
&de(7,2)+R(2,3)*tt5*dNude(6,2)+R(2,2)*tt5*dNude(5,2)+R(2,1)*tt5*dN&
&ude(4,2)+R(2,3)*tt2*dNude(3,2)+dNxde(1,2)*tt11+R(2,2)*dNude(2,2)*&
&tt2+dNude(1,2)*R(2,1)*tt2
tt14 = R(2,3)*tt7*dNude(9,3)+R(2,2)*tt7*dNude(8,3)+R(2,1)*tt7*dNu&
&de(7,3)+R(2,3)*tt5*dNude(6,3)+R(2,2)*tt5*dNude(5,3)+R(2,1)*tt5*dN&
&ude(4,3)+R(2,3)*tt2*dNude(3,3)+dNxde(1,3)*tt11+R(2,2)*dNude(2,3)*&
&tt2+dNude(1,3)*R(2,1)*tt2
tt15 = X0(3,1)*R(3,3)+X0(2,1)*R(3,2)+tt6+X0(1,1)*R(3,1)
tt16 = R(3,3)*tt7*dNude(9,1)+R(3,2)*tt7*dNude(8,1)+R(3,1)*tt7*dNu&
&de(7,1)+tt5*R(3,3)*dNude(6,1)+R(3,2)*tt5*dNude(5,1)+R(3,1)*tt5*dN&
&ude(4,1)+dNxde(1,1)*tt15+tt2*dNude(3,1)*R(3,3)+dNude(2,1)*tt2*R(3&
&,2)+dNude(1,1)*R(3,1)*tt2
tt17 = R(3,3)*tt7*dNude(9,2)+R(3,2)*tt7*dNude(8,2)+R(3,1)*tt7*dNu&
&de(7,2)+tt5*R(3,3)*dNude(6,2)+R(3,2)*tt5*dNude(5,2)+R(3,1)*tt5*dN&
&ude(4,2)+dNxde(1,2)*tt15+tt2*dNude(3,2)*R(3,3)+dNude(2,2)*tt2*R(3&
&,2)+dNude(1,2)*R(3,1)*tt2
tt18 = R(3,3)*tt7*dNude(9,3)+R(3,2)*tt7*dNude(8,3)+R(3,1)*tt7*dNu&
&de(7,3)+tt5*R(3,3)*dNude(6,3)+R(3,2)*tt5*dNude(5,3)+R(3,1)*tt5*dN&
&ude(4,3)+tt2*R(3,3)*dNude(3,3)+dNxde(1,3)*tt15+dNude(2,3)*tt2*R(3&
&,2)+dNude(1,3)*R(3,1)*tt2
val(1,1) = Dm(3,1)*tt10+Dm(2,1)*tt9+Dm(1,1)*tt8
val(1,2) = Dm(3,1)*tt14+Dm(2,1)*tt13+Dm(1,1)*tt12
val(1,3) = Dm(3,1)*tt18+Dm(2,1)*tt17+Dm(1,1)*tt16
val(1,4) = Dm(3,2)*tt10+Dm(2,2)*tt9+Dm(1,2)*tt8
val(1,5) = Dm(3,2)*tt14+Dm(2,2)*tt13+Dm(1,2)*tt12
val(1,6) = Dm(3,2)*tt18+Dm(2,2)*tt17+Dm(1,2)*tt16
val(1,7) = Dm(3,3)*tt10+Dm(2,3)*tt9+Dm(1,3)*tt8
val(1,8) = Dm(3,3)*tt14+Dm(2,3)*tt13+Dm(1,3)*tt12
val(1,9) = Dm(3,3)*tt18+Dm(2,3)*tt17+Dm(1,3)*tt16
END 
SUBROUTINE coro_aug_ani_df_3d_jac(jac, X, X0, R, Nu, Nx, dNude, dNxde, Dm) 
IMPLICIT NONE 
REAL(KIND=8) jac(9, 3) 
REAL(KIND=8) X(3, 1) 
REAL(KIND=8) X0(3, 1) 
REAL(KIND=8) R(3, 3) 
REAL(KIND=8) Nu(3, 3) 
REAL(KIND=8) Nx(1, 1) 
REAL(KIND=8) dNude(9, 3) 
REAL(KIND=8) dNxde(1, 3) 
REAL(KIND=8) Dm(3, 3) 
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
tt1 = R(1,1)**2
tt2 = R(1,2)**2
tt3 = R(1,3)**2
tt4 = tt3*dNude(9,1)+R(1,2)*R(1,3)*dNude(8,1)+R(1,1)*R(1,3)*dNude&
&(7,1)+R(1,2)*R(1,3)*dNude(6,1)+tt2*dNude(5,1)+R(1,1)*R(1,2)*dNude&
&(4,1)+R(1,1)*R(1,3)*dNude(3,1)+R(1,1)*R(1,2)*dNude(2,1)+tt1*dNude&
&(1,1)
tt5 = tt3*dNude(9,2)+R(1,2)*R(1,3)*dNude(8,2)+R(1,1)*R(1,3)*dNude&
&(7,2)+R(1,2)*R(1,3)*dNude(6,2)+tt2*dNude(5,2)+R(1,1)*R(1,2)*dNude&
&(4,2)+R(1,1)*R(1,3)*dNude(3,2)+R(1,1)*R(1,2)*dNude(2,2)+tt1*dNude&
&(1,2)
tt6 = tt3*dNude(9,3)+R(1,2)*R(1,3)*dNude(8,3)+R(1,1)*R(1,3)*dNude&
&(7,3)+R(1,2)*R(1,3)*dNude(6,3)+tt2*dNude(5,3)+R(1,1)*R(1,2)*dNude&
&(4,3)+R(1,1)*R(1,3)*dNude(3,3)+R(1,1)*R(1,2)*dNude(2,3)+tt1*dNude&
&(1,3)
tt7 = R(1,1)*dNude(1,1)*R(2,1)
tt8 = R(1,2)*R(2,2)*dNude(5,1)
tt9 = R(1,3)*R(2,3)*dNude(9,1)
tt10 = tt9+R(1,2)*R(2,3)*dNude(8,1)+R(1,1)*R(2,3)*dNude(7,1)+R(1,&
&3)*R(2,2)*dNude(6,1)+tt8+R(1,1)*R(2,2)*dNude(4,1)+R(1,3)*R(2,1)*d&
&Nude(3,1)+R(1,2)*R(2,1)*dNude(2,1)+tt7
tt11 = R(1,1)*dNude(1,2)*R(2,1)
tt12 = R(1,2)*R(2,2)*dNude(5,2)
tt13 = R(1,3)*R(2,3)*dNude(9,2)
tt14 = tt13+R(1,2)*R(2,3)*dNude(8,2)+R(1,1)*R(2,3)*dNude(7,2)+R(1&
&,3)*R(2,2)*dNude(6,2)+tt12+R(1,1)*R(2,2)*dNude(4,2)+R(1,3)*R(2,1)&
&*dNude(3,2)+R(1,2)*R(2,1)*dNude(2,2)+tt11
tt15 = R(1,1)*dNude(1,3)*R(2,1)
tt16 = R(1,2)*R(2,2)*dNude(5,3)
tt17 = R(1,3)*R(2,3)*dNude(9,3)
tt18 = tt17+R(1,2)*R(2,3)*dNude(8,3)+R(1,1)*R(2,3)*dNude(7,3)+R(1&
&,3)*R(2,2)*dNude(6,3)+tt16+R(1,1)*R(2,2)*dNude(4,3)+R(1,3)*R(2,1)&
&*dNude(3,3)+R(1,2)*R(2,1)*dNude(2,3)+tt15
tt19 = R(1,1)*dNude(1,1)*R(3,1)
tt20 = R(1,2)*R(3,2)*dNude(5,1)
tt21 = R(1,3)*R(3,3)*dNude(9,1)
tt22 = tt21+R(1,2)*R(3,3)*dNude(8,1)+R(1,1)*R(3,3)*dNude(7,1)+R(1&
&,3)*R(3,2)*dNude(6,1)+tt20+R(1,1)*R(3,2)*dNude(4,1)+R(1,3)*R(3,1)&
&*dNude(3,1)+R(1,2)*dNude(2,1)*R(3,1)+tt19
tt23 = R(1,1)*dNude(1,2)*R(3,1)
tt24 = R(1,2)*R(3,2)*dNude(5,2)
tt25 = R(1,3)*R(3,3)*dNude(9,2)
tt26 = tt25+R(1,2)*R(3,3)*dNude(8,2)+R(1,1)*R(3,3)*dNude(7,2)+R(1&
&,3)*R(3,2)*dNude(6,2)+tt24+R(1,1)*R(3,2)*dNude(4,2)+R(1,3)*R(3,1)&
&*dNude(3,2)+R(1,2)*dNude(2,2)*R(3,1)+tt23
tt27 = R(1,1)*dNude(1,3)*R(3,1)
tt28 = R(1,2)*R(3,2)*dNude(5,3)
tt29 = R(1,3)*R(3,3)*dNude(9,3)
tt30 = tt29+R(1,2)*R(3,3)*dNude(8,3)+R(1,1)*R(3,3)*dNude(7,3)+R(1&
&,3)*R(3,2)*dNude(6,3)+tt28+R(1,1)*R(3,2)*dNude(4,3)+R(1,3)*R(3,1)&
&*dNude(3,3)+R(1,2)*dNude(2,3)*R(3,1)+tt27
tt31 = tt9+R(1,3)*R(2,2)*dNude(8,1)+R(1,3)*R(2,1)*dNude(7,1)+R(1,&
&2)*R(2,3)*dNude(6,1)+tt8+R(1,2)*R(2,1)*dNude(4,1)+R(1,1)*R(2,3)*d&
&Nude(3,1)+R(1,1)*dNude(2,1)*R(2,2)+tt7
tt32 = tt13+R(1,3)*R(2,2)*dNude(8,2)+R(1,3)*R(2,1)*dNude(7,2)+R(1&
&,2)*R(2,3)*dNude(6,2)+tt12+R(1,2)*R(2,1)*dNude(4,2)+R(1,1)*R(2,3)&
&*dNude(3,2)+R(1,1)*R(2,2)*dNude(2,2)+tt11
tt33 = tt17+R(1,3)*R(2,2)*dNude(8,3)+R(1,3)*R(2,1)*dNude(7,3)+R(1&
&,2)*R(2,3)*dNude(6,3)+tt16+R(1,2)*R(2,1)*dNude(4,3)+R(1,1)*R(2,3)&
&*dNude(3,3)+R(1,1)*R(2,2)*dNude(2,3)+tt15
tt34 = R(2,1)**2
tt35 = R(2,2)**2
tt36 = R(2,3)**2
tt37 = tt36*dNude(9,1)+R(2,2)*R(2,3)*dNude(8,1)+R(2,1)*R(2,3)*dNu&
&de(7,1)+R(2,2)*R(2,3)*dNude(6,1)+tt35*dNude(5,1)+R(2,1)*R(2,2)*dN&
&ude(4,1)+R(2,1)*R(2,3)*dNude(3,1)+R(2,1)*dNude(2,1)*R(2,2)+dNude(&
&1,1)*tt34
tt38 = tt36*dNude(9,2)+R(2,2)*R(2,3)*dNude(8,2)+R(2,1)*R(2,3)*dNu&
&de(7,2)+R(2,2)*R(2,3)*dNude(6,2)+tt35*dNude(5,2)+R(2,1)*R(2,2)*dN&
&ude(4,2)+R(2,1)*R(2,3)*dNude(3,2)+R(2,1)*R(2,2)*dNude(2,2)+dNude(&
&1,2)*tt34
tt39 = tt36*dNude(9,3)+R(2,2)*R(2,3)*dNude(8,3)+R(2,1)*R(2,3)*dNu&
&de(7,3)+R(2,2)*R(2,3)*dNude(6,3)+tt35*dNude(5,3)+R(2,1)*R(2,2)*dN&
&ude(4,3)+R(2,1)*R(2,3)*dNude(3,3)+R(2,1)*R(2,2)*dNude(2,3)+dNude(&
&1,3)*tt34
tt40 = dNude(1,1)*R(2,1)*R(3,1)
tt41 = R(2,2)*R(3,2)*dNude(5,1)
tt42 = R(2,3)*R(3,3)*dNude(9,1)
tt43 = tt42+R(2,2)*R(3,3)*dNude(8,1)+R(2,1)*R(3,3)*dNude(7,1)+R(2&
&,3)*R(3,2)*dNude(6,1)+tt41+R(2,1)*R(3,2)*dNude(4,1)+R(2,3)*R(3,1)&
&*dNude(3,1)+dNude(2,1)*R(2,2)*R(3,1)+tt40
tt44 = dNude(1,2)*R(2,1)*R(3,1)
tt45 = R(2,2)*R(3,2)*dNude(5,2)
tt46 = R(2,3)*R(3,3)*dNude(9,2)
tt47 = tt46+R(2,2)*R(3,3)*dNude(8,2)+R(2,1)*R(3,3)*dNude(7,2)+R(2&
&,3)*R(3,2)*dNude(6,2)+tt45+R(2,1)*R(3,2)*dNude(4,2)+R(2,3)*R(3,1)&
&*dNude(3,2)+R(2,2)*dNude(2,2)*R(3,1)+tt44
tt48 = dNude(1,3)*R(2,1)*R(3,1)
tt49 = R(2,2)*R(3,2)*dNude(5,3)
tt50 = R(2,3)*R(3,3)*dNude(9,3)
tt51 = tt50+R(2,2)*R(3,3)*dNude(8,3)+R(2,1)*R(3,3)*dNude(7,3)+R(2&
&,3)*R(3,2)*dNude(6,3)+tt49+R(2,1)*R(3,2)*dNude(4,3)+R(2,3)*R(3,1)&
&*dNude(3,3)+R(2,2)*dNude(2,3)*R(3,1)+tt48
tt52 = tt21+R(1,3)*R(3,2)*dNude(8,1)+R(1,3)*R(3,1)*dNude(7,1)+R(1&
&,2)*R(3,3)*dNude(6,1)+tt20+R(1,2)*R(3,1)*dNude(4,1)+R(1,1)*dNude(&
&3,1)*R(3,3)+R(1,1)*dNude(2,1)*R(3,2)+tt19
tt53 = tt25+R(1,3)*R(3,2)*dNude(8,2)+R(1,3)*R(3,1)*dNude(7,2)+R(1&
&,2)*R(3,3)*dNude(6,2)+tt24+R(1,2)*R(3,1)*dNude(4,2)+R(1,1)*dNude(&
&3,2)*R(3,3)+R(1,1)*dNude(2,2)*R(3,2)+tt23
tt54 = tt29+R(1,3)*R(3,2)*dNude(8,3)+R(1,3)*R(3,1)*dNude(7,3)+R(1&
&,2)*R(3,3)*dNude(6,3)+tt28+R(1,2)*R(3,1)*dNude(4,3)+R(1,1)*R(3,3)&
&*dNude(3,3)+R(1,1)*dNude(2,3)*R(3,2)+tt27
tt55 = tt42+R(2,3)*R(3,2)*dNude(8,1)+R(2,3)*R(3,1)*dNude(7,1)+R(2&
&,2)*R(3,3)*dNude(6,1)+tt41+R(2,2)*R(3,1)*dNude(4,1)+R(2,1)*dNude(&
&3,1)*R(3,3)+R(2,1)*dNude(2,1)*R(3,2)+tt40
tt56 = tt46+R(2,3)*R(3,2)*dNude(8,2)+R(2,3)*R(3,1)*dNude(7,2)+R(2&
&,2)*R(3,3)*dNude(6,2)+tt45+R(2,2)*R(3,1)*dNude(4,2)+R(2,1)*dNude(&
&3,2)*R(3,3)+R(2,1)*dNude(2,2)*R(3,2)+tt44
tt57 = tt50+R(2,3)*R(3,2)*dNude(8,3)+R(2,3)*R(3,1)*dNude(7,3)+R(2&
&,2)*R(3,3)*dNude(6,3)+tt49+R(2,2)*R(3,1)*dNude(4,3)+R(2,1)*R(3,3)&
&*dNude(3,3)+R(2,1)*dNude(2,3)*R(3,2)+tt48
tt58 = R(3,1)**2
tt59 = R(3,2)**2
tt60 = R(3,3)**2
tt61 = tt60*dNude(9,1)+R(3,2)*R(3,3)*dNude(8,1)+R(3,1)*R(3,3)*dNu&
&de(7,1)+R(3,2)*R(3,3)*dNude(6,1)+tt59*dNude(5,1)+R(3,1)*R(3,2)*dN&
&ude(4,1)+R(3,1)*dNude(3,1)*R(3,3)+dNude(2,1)*R(3,1)*R(3,2)+dNude(&
&1,1)*tt58
tt62 = tt60*dNude(9,2)+R(3,2)*R(3,3)*dNude(8,2)+R(3,1)*R(3,3)*dNu&
&de(7,2)+R(3,2)*R(3,3)*dNude(6,2)+tt59*dNude(5,2)+R(3,1)*R(3,2)*dN&
&ude(4,2)+R(3,1)*dNude(3,2)*R(3,3)+dNude(2,2)*R(3,1)*R(3,2)+dNude(&
&1,2)*tt58
tt63 = tt60*dNude(9,3)+R(3,2)*R(3,3)*dNude(8,3)+R(3,1)*R(3,3)*dNu&
&de(7,3)+R(3,2)*R(3,3)*dNude(6,3)+tt59*dNude(5,3)+R(3,1)*R(3,2)*dN&
&ude(4,3)+R(3,1)*R(3,3)*dNude(3,3)+dNude(2,3)*R(3,1)*R(3,2)+dNude(&
&1,3)*tt58
jac(1,1) = Dm(3,1)*tt6+Dm(2,1)*tt5+Dm(1,1)*tt4
jac(1,2) = Dm(3,1)*tt18+Dm(2,1)*tt14+Dm(1,1)*tt10
jac(1,3) = Dm(3,1)*tt30+Dm(2,1)*tt26+Dm(1,1)*tt22
jac(2,1) = Dm(3,1)*tt33+Dm(2,1)*tt32+Dm(1,1)*tt31
jac(2,2) = Dm(3,1)*tt39+Dm(2,1)*tt38+Dm(1,1)*tt37
jac(2,3) = Dm(3,1)*tt51+Dm(2,1)*tt47+Dm(1,1)*tt43
jac(3,1) = Dm(3,1)*tt54+Dm(2,1)*tt53+Dm(1,1)*tt52
jac(3,2) = Dm(3,1)*tt57+Dm(2,1)*tt56+Dm(1,1)*tt55
jac(3,3) = Dm(3,1)*tt63+Dm(2,1)*tt62+Dm(1,1)*tt61
jac(4,1) = Dm(3,2)*tt6+Dm(2,2)*tt5+Dm(1,2)*tt4
jac(4,2) = Dm(3,2)*tt18+Dm(2,2)*tt14+Dm(1,2)*tt10
jac(4,3) = Dm(3,2)*tt30+Dm(2,2)*tt26+Dm(1,2)*tt22
jac(5,1) = Dm(3,2)*tt33+Dm(2,2)*tt32+Dm(1,2)*tt31
jac(5,2) = Dm(3,2)*tt39+Dm(2,2)*tt38+Dm(1,2)*tt37
jac(5,3) = Dm(3,2)*tt51+Dm(2,2)*tt47+Dm(1,2)*tt43
jac(6,1) = Dm(3,2)*tt54+Dm(2,2)*tt53+Dm(1,2)*tt52
jac(6,2) = Dm(3,2)*tt57+Dm(2,2)*tt56+Dm(1,2)*tt55
jac(6,3) = Dm(3,2)*tt63+Dm(2,2)*tt62+Dm(1,2)*tt61
jac(7,1) = Dm(3,3)*tt6+Dm(2,3)*tt5+Dm(1,3)*tt4
jac(7,2) = Dm(3,3)*tt18+Dm(2,3)*tt14+Dm(1,3)*tt10
jac(7,3) = Dm(3,3)*tt30+Dm(2,3)*tt26+Dm(1,3)*tt22
jac(8,1) = Dm(3,3)*tt33+Dm(2,3)*tt32+Dm(1,3)*tt31
jac(8,2) = Dm(3,3)*tt39+Dm(2,3)*tt38+Dm(1,3)*tt37
jac(8,3) = Dm(3,3)*tt51+Dm(2,3)*tt47+Dm(1,3)*tt43
jac(9,1) = Dm(3,3)*tt54+Dm(2,3)*tt53+Dm(1,3)*tt52
jac(9,2) = Dm(3,3)*tt57+Dm(2,3)*tt56+Dm(1,3)*tt55
jac(9,3) = Dm(3,3)*tt63+Dm(2,3)*tt62+Dm(1,3)*tt61
END 
SUBROUTINE coro_iso_basis_df_3d(val, X, X0, R, Nu, Nx, dNude, dNxde, Dm) 
IMPLICIT NONE 
REAL(KIND=8) val(1, 9) 
REAL(KIND=8) X(3, 1) 
REAL(KIND=8) X0(3, 1) 
REAL(KIND=8) R(3, 3) 
REAL(KIND=8) Nu(1, 1) 
REAL(KIND=8) Nx(1, 1) 
REAL(KIND=8) dNude(1, 3) 
REAL(KIND=8) dNxde(1, 3) 
REAL(KIND=8) Dm(3, 3) 
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
tt1 = -X0(1,1)
tt2 = R(1,3)*X0(3,1)+R(1,2)*X0(2,1)+R(1,1)*X0(1,1)+tt1
tt3 = R(3,1)*X(3,1)+R(2,1)*X(2,1)+tt1+R(1,1)*X(1,1)
tt4 = -X0(2,1)
tt5 = X(3,1)*R(3,2)+X(2,1)*R(2,2)+tt4+X(1,1)*R(1,2)
tt6 = -X0(3,1)
tt7 = X(3,1)*R(3,3)+tt6+X(2,1)*R(2,3)+X(1,1)*R(1,3)
tt8 = R(1,3)*tt7+R(1,2)*tt5+R(1,1)*tt3
tt9 = dNude(1,1)*tt8+dNxde(1,1)*tt2
tt10 = dNude(1,2)*tt8+dNxde(1,2)*tt2
tt11 = dNude(1,3)*tt8+dNxde(1,3)*tt2
tt12 = R(2,3)*X0(3,1)+X0(2,1)*R(2,2)+tt4+X0(1,1)*R(2,1)
tt13 = R(2,3)*tt7+R(2,2)*tt5+R(2,1)*tt3
tt14 = dNude(1,1)*tt13+dNxde(1,1)*tt12
tt15 = dNude(1,2)*tt13+dNxde(1,2)*tt12
tt16 = dNude(1,3)*tt13+dNxde(1,3)*tt12
tt17 = X0(3,1)*R(3,3)+X0(2,1)*R(3,2)+tt6+X0(1,1)*R(3,1)
tt18 = R(3,3)*tt7+R(3,2)*tt5+R(3,1)*tt3
tt19 = dNude(1,1)*tt18+dNxde(1,1)*tt17
tt20 = dNude(1,2)*tt18+dNxde(1,2)*tt17
tt21 = dNude(1,3)*tt18+dNxde(1,3)*tt17
val(1,1) = Dm(3,1)*tt11+Dm(2,1)*tt10+Dm(1,1)*tt9
val(1,2) = Dm(3,1)*tt16+Dm(2,1)*tt15+Dm(1,1)*tt14
val(1,3) = Dm(3,1)*tt21+Dm(2,1)*tt20+Dm(1,1)*tt19
val(1,4) = Dm(3,2)*tt11+Dm(2,2)*tt10+Dm(1,2)*tt9
val(1,5) = Dm(3,2)*tt16+Dm(2,2)*tt15+Dm(1,2)*tt14
val(1,6) = Dm(3,2)*tt21+Dm(2,2)*tt20+Dm(1,2)*tt19
val(1,7) = Dm(3,3)*tt11+Dm(2,3)*tt10+Dm(1,3)*tt9
val(1,8) = Dm(3,3)*tt16+Dm(2,3)*tt15+Dm(1,3)*tt14
val(1,9) = Dm(3,3)*tt21+Dm(2,3)*tt20+Dm(1,3)*tt19
END 
SUBROUTINE coro_iso_basis_df_3d_jac(jac, X, X0, R, Nu, Nx, dNude, dNxde, Dm) 
IMPLICIT NONE 
REAL(KIND=8) jac(9, 3) 
REAL(KIND=8) X(3, 1) 
REAL(KIND=8) X0(3, 1) 
REAL(KIND=8) R(3, 3) 
REAL(KIND=8) Nu(1, 1) 
REAL(KIND=8) Nx(1, 1) 
REAL(KIND=8) dNude(1, 3) 
REAL(KIND=8) dNxde(1, 3) 
REAL(KIND=8) Dm(3, 3) 
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
tt1 = R(1,3)**2+R(1,2)**2+R(1,1)**2
tt2 = R(1,3)*R(2,3)+R(1,2)*R(2,2)+R(1,1)*R(2,1)
tt3 = dNude(1,3)*tt2*Dm(3,1)+dNude(1,2)*Dm(2,1)*tt2+Dm(1,1)*dNude&
&(1,1)*tt2
tt4 = R(1,3)*R(3,3)+R(1,2)*R(3,2)+R(1,1)*R(3,1)
tt5 = dNude(1,3)*Dm(3,1)*tt4+dNude(1,2)*Dm(2,1)*tt4+Dm(1,1)*dNude&
&(1,1)*tt4
tt6 = R(2,3)**2+R(2,2)**2+R(2,1)**2
tt7 = R(2,3)*R(3,3)+R(2,2)*R(3,2)+R(2,1)*R(3,1)
tt8 = dNude(1,3)*Dm(3,1)*tt7+dNude(1,2)*Dm(2,1)*tt7+Dm(1,1)*dNude&
&(1,1)*tt7
tt9 = R(3,3)**2+R(3,2)**2+R(3,1)**2
tt10 = dNude(1,3)*tt2*Dm(3,2)+dNude(1,2)*Dm(2,2)*tt2+dNude(1,1)*D&
&m(1,2)*tt2
tt11 = dNude(1,3)*Dm(3,2)*tt4+dNude(1,2)*Dm(2,2)*tt4+dNude(1,1)*D&
&m(1,2)*tt4
tt12 = dNude(1,3)*Dm(3,2)*tt7+dNude(1,2)*Dm(2,2)*tt7+dNude(1,1)*D&
&m(1,2)*tt7
tt13 = dNude(1,3)*tt2*Dm(3,3)+dNude(1,2)*Dm(2,3)*tt2+dNude(1,1)*D&
&m(1,3)*tt2
tt14 = dNude(1,3)*Dm(3,3)*tt4+dNude(1,2)*Dm(2,3)*tt4+dNude(1,1)*D&
&m(1,3)*tt4
tt15 = dNude(1,3)*Dm(3,3)*tt7+dNude(1,2)*Dm(2,3)*tt7+dNude(1,1)*D&
&m(1,3)*tt7
jac(1,1) = tt1*dNude(1,3)*Dm(3,1)+dNude(1,2)*tt1*Dm(2,1)+Dm(1,1)*&
&dNude(1,1)*tt1
jac(1,2) = tt3
jac(1,3) = tt5
jac(2,1) = tt3
jac(2,2) = dNude(1,3)*tt6*Dm(3,1)+dNude(1,2)*Dm(2,1)*tt6+Dm(1,1)*&
&dNude(1,1)*tt6
jac(2,3) = tt8
jac(3,1) = tt5
jac(3,2) = tt8
jac(3,3) = dNude(1,3)*Dm(3,1)*tt9+dNude(1,2)*Dm(2,1)*tt9+Dm(1,1)*&
&dNude(1,1)*tt9
jac(4,1) = tt1*dNude(1,3)*Dm(3,2)+dNude(1,2)*tt1*Dm(2,2)+dNude(1,&
&1)*Dm(1,2)*tt1
jac(4,2) = tt10
jac(4,3) = tt11
jac(5,1) = tt10
jac(5,2) = dNude(1,3)*tt6*Dm(3,2)+dNude(1,2)*Dm(2,2)*tt6+dNude(1,&
&1)*Dm(1,2)*tt6
jac(5,3) = tt12
jac(6,1) = tt11
jac(6,2) = tt12
jac(6,3) = dNude(1,3)*Dm(3,2)*tt9+dNude(1,2)*Dm(2,2)*tt9+dNude(1,&
&1)*Dm(1,2)*tt9
jac(7,1) = tt1*dNude(1,3)*Dm(3,3)+dNude(1,2)*tt1*Dm(2,3)+dNude(1,&
&1)*Dm(1,3)*tt1
jac(7,2) = tt13
jac(7,3) = tt14
jac(8,1) = tt13
jac(8,2) = dNude(1,3)*tt6*Dm(3,3)+dNude(1,2)*Dm(2,3)*tt6+dNude(1,&
&1)*Dm(1,3)*tt6
jac(8,3) = tt15
jac(9,1) = tt14
jac(9,2) = tt15
jac(9,3) = dNude(1,3)*Dm(3,3)*tt9+dNude(1,2)*Dm(2,3)*tt9+dNude(1,&
&1)*Dm(1,3)*tt9
END 
