SUBROUTINE go_thru_nnet(val, W, X) 
IMPLICIT NONE 
REAL(KIND=8) val(1, 1) 
REAL(KIND=8) W(31, 1) 
REAL(KIND=8) X(3, 1) 
val(1,1) = W(31,1)+W(30,1)/(exp((-W(24,1))-X(3,1)*W(18,1)-X(2,1)*&
&W(12,1)-X(1,1)*W(6,1))+1)+W(29,1)/(exp((-W(23,1))-X(3,1)*W(17,1)-&
&X(2,1)*W(11,1)-X(1,1)*W(5,1))+1)+W(28,1)/(exp((-W(22,1))-X(3,1)*W&
&(16,1)-X(2,1)*W(10,1)-X(1,1)*W(4,1))+1)+W(27,1)/(exp((-W(21,1))-X&
&(3,1)*W(15,1)-X(2,1)*W(9,1)-X(1,1)*W(3,1))+1)+W(26,1)/(exp((-W(20&
&,1))-X(3,1)*W(14,1)-X(2,1)*W(8,1)-X(1,1)*W(2,1))+1)+W(25,1)/(exp(&
&(-W(19,1))-X(3,1)*W(13,1)-X(2,1)*W(7,1)-W(1,1)*X(1,1))+1)
END 
SUBROUTINE nnet_squared_error(val, W, X, Y) 
IMPLICIT NONE 
REAL(KIND=8) val(1, 1) 
REAL(KIND=8) W(31, 1) 
REAL(KIND=8) X(3, 1) 
REAL(KIND=8) Y(1, 1) 
val(1,1) = (W(31,1)+W(30,1)/(exp((-W(24,1))-X(3,1)*W(18,1)-X(2,1)&
&*W(12,1)-X(1,1)*W(6,1))+1)+W(29,1)/(exp((-W(23,1))-X(3,1)*W(17,1)&
&-X(2,1)*W(11,1)-X(1,1)*W(5,1))+1)+W(28,1)/(exp((-W(22,1))-X(3,1)*&
&W(16,1)-X(2,1)*W(10,1)-X(1,1)*W(4,1))+1)+W(27,1)/(exp((-W(21,1))-&
&X(3,1)*W(15,1)-X(2,1)*W(9,1)-X(1,1)*W(3,1))+1)+W(26,1)/(exp((-W(2&
&0,1))-X(3,1)*W(14,1)-X(2,1)*W(8,1)-X(1,1)*W(2,1))+1)+W(25,1)/(exp&
&((-W(19,1))-X(3,1)*W(13,1)-X(2,1)*W(7,1)-W(1,1)*X(1,1))+1)-Y(1,1)&
&)**2
END 
SUBROUTINE nnet_squared_error_jac(jac, W, X, Y) 
IMPLICIT NONE 
REAL(KIND=8) jac(1, 31) 
REAL(KIND=8) W(31, 1) 
REAL(KIND=8) X(3, 1) 
REAL(KIND=8) Y(1, 1) 
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
tt1 = exp((-W(19,1))-X(3,1)*W(13,1)-X(2,1)*W(7,1)-W(1,1)*X(1,1))
tt2 = tt1+1
tt3 = 1/tt2**2
tt4 = 1/tt2
tt5 = exp((-W(20,1))-X(3,1)*W(14,1)-X(2,1)*W(8,1)-X(1,1)*W(2,1))
tt6 = tt5+1
tt7 = 1/tt6
tt8 = exp((-W(21,1))-X(3,1)*W(15,1)-X(2,1)*W(9,1)-X(1,1)*W(3,1))
tt9 = tt8+1
tt10 = 1/tt9
tt11 = exp((-W(22,1))-X(3,1)*W(16,1)-X(2,1)*W(10,1)-X(1,1)*W(4,1)&
&)
tt12 = tt11+1
tt13 = 1/tt12
tt14 = exp((-W(23,1))-X(3,1)*W(17,1)-X(2,1)*W(11,1)-X(1,1)*W(5,1)&
&)
tt15 = tt14+1
tt16 = 1/tt15
tt17 = exp((-W(24,1))-X(3,1)*W(18,1)-X(2,1)*W(12,1)-X(1,1)*W(6,1)&
&)
tt18 = tt17+1
tt19 = 1/tt18
tt20 = W(31,1)+tt19*W(30,1)+tt16*W(29,1)+tt13*W(28,1)+tt10*W(27,1&
&)+tt7*W(26,1)+tt4*W(25,1)-Y(1,1)
tt21 = 1/tt6**2
tt22 = 1/tt9**2
tt23 = 1/tt12**2
tt24 = 1/tt15**2
tt25 = 1/tt18**2
jac(1,1) = 2*tt1*tt3*X(1,1)*W(25,1)*tt20
jac(1,2) = 2*tt5*tt21*X(1,1)*W(26,1)*tt20
jac(1,3) = 2*tt8*tt22*X(1,1)*W(27,1)*tt20
jac(1,4) = 2*tt11*tt23*X(1,1)*W(28,1)*tt20
jac(1,5) = 2*tt14*tt24*X(1,1)*W(29,1)*tt20
jac(1,6) = 2*tt17*tt25*X(1,1)*W(30,1)*tt20
jac(1,7) = 2*tt1*tt3*X(2,1)*W(25,1)*tt20
jac(1,8) = 2*tt5*tt21*X(2,1)*W(26,1)*tt20
jac(1,9) = 2*tt8*tt22*X(2,1)*W(27,1)*tt20
jac(1,10) = 2*tt11*tt23*X(2,1)*W(28,1)*tt20
jac(1,11) = 2*tt14*tt24*X(2,1)*W(29,1)*tt20
jac(1,12) = 2*tt17*tt25*X(2,1)*W(30,1)*tt20
jac(1,13) = 2*tt1*tt3*X(3,1)*W(25,1)*tt20
jac(1,14) = 2*tt5*tt21*X(3,1)*W(26,1)*tt20
jac(1,15) = 2*tt8*tt22*X(3,1)*W(27,1)*tt20
jac(1,16) = 2*tt11*tt23*X(3,1)*W(28,1)*tt20
jac(1,17) = 2*tt14*tt24*X(3,1)*W(29,1)*tt20
jac(1,18) = 2*tt17*tt25*X(3,1)*W(30,1)*tt20
jac(1,19) = 2*tt1*tt3*W(25,1)*tt20
jac(1,20) = 2*tt5*tt21*W(26,1)*tt20
jac(1,21) = 2*tt8*tt22*W(27,1)*tt20
jac(1,22) = 2*tt11*tt23*W(28,1)*tt20
jac(1,23) = 2*tt14*tt24*W(29,1)*tt20
jac(1,24) = 2*tt17*tt25*W(30,1)*tt20
jac(1,25) = 2*tt4*tt20
jac(1,26) = 2*tt7*tt20
jac(1,27) = 2*tt10*tt20
jac(1,28) = 2*tt13*tt20
jac(1,29) = 2*tt16*tt20
jac(1,30) = 2*tt19*tt20
jac(1,31) = 2*tt20
END 
SUBROUTINE nnet_squared_error_hes(hes, W, X, Y) 
IMPLICIT NONE 
REAL(KIND=8) hes(31, 31) 
REAL(KIND=8) W(31, 1) 
REAL(KIND=8) X(3, 1) 
REAL(KIND=8) Y(1, 1) 
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
tt1 = exp((-2*W(19,1))-2*X(3,1)*W(13,1)-2*X(2,1)*W(7,1)-2*W(1,1)*&
&X(1,1))
tt2 = -W(1,1)*X(1,1)
tt3 = -X(2,1)*W(7,1)
tt4 = -X(3,1)*W(13,1)
tt5 = -W(19,1)
tt6 = exp(tt5+tt4+tt3+tt2)
tt7 = tt6+1
tt8 = 1/tt7**4
tt9 = X(1,1)**2
tt10 = W(25,1)**2
tt11 = 1/tt7**3
tt12 = 1/tt7
tt13 = -X(1,1)*W(2,1)
tt14 = -X(2,1)*W(8,1)
tt15 = -X(3,1)*W(14,1)
tt16 = -W(20,1)
tt17 = exp(tt16+tt15+tt14+tt13)
tt18 = tt17+1
tt19 = 1/tt18
tt20 = -X(1,1)*W(3,1)
tt21 = -X(2,1)*W(9,1)
tt22 = -X(3,1)*W(15,1)
tt23 = -W(21,1)
tt24 = exp(tt23+tt22+tt21+tt20)
tt25 = tt24+1
tt26 = 1/tt25
tt27 = -X(1,1)*W(4,1)
tt28 = -X(2,1)*W(10,1)
tt29 = -X(3,1)*W(16,1)
tt30 = -W(22,1)
tt31 = exp(tt30+tt29+tt28+tt27)
tt32 = tt31+1
tt33 = 1/tt32
tt34 = -X(1,1)*W(5,1)
tt35 = -X(2,1)*W(11,1)
tt36 = -X(3,1)*W(17,1)
tt37 = -W(23,1)
tt38 = exp(tt37+tt36+tt35+tt34)
tt39 = tt38+1
tt40 = 1/tt39
tt41 = -X(1,1)*W(6,1)
tt42 = -X(2,1)*W(12,1)
tt43 = -X(3,1)*W(18,1)
tt44 = -W(24,1)
tt45 = exp(tt44+tt43+tt42+tt41)
tt46 = tt45+1
tt47 = 1/tt46
tt48 = W(31,1)+tt47*W(30,1)+tt40*W(29,1)+tt33*W(28,1)+tt26*W(27,1&
&)+tt19*W(26,1)+tt12*W(25,1)-Y(1,1)
tt49 = 1/tt7**2
tt50 = 1/tt18**2
tt51 = exp(tt16+tt5+tt15+tt4+tt14+tt3+tt13+tt2)
tt52 = 2*tt49*tt50*tt51*tt9*W(25,1)*W(26,1)
tt53 = 1/tt25**2
tt54 = exp(tt23+tt5+tt22+tt4+tt21+tt3+tt20+tt2)
tt55 = 2*tt49*tt53*tt54*tt9*W(25,1)*W(27,1)
tt56 = 1/tt32**2
tt57 = exp(tt30+tt5+tt29+tt4+tt28+tt3+tt27+tt2)
tt58 = 2*tt49*tt56*tt57*tt9*W(25,1)*W(28,1)
tt59 = 1/tt39**2
tt60 = exp(tt37+tt5+tt36+tt4+tt35+tt3+tt34+tt2)
tt61 = 2*tt49*tt59*tt60*tt9*W(25,1)*W(29,1)
tt62 = 1/tt46**2
tt63 = exp(tt44+tt5+tt43+tt4+tt42+tt3+tt41+tt2)
tt64 = 2*tt49*tt62*tt63*tt9*W(25,1)*W(30,1)
tt65 = (-2*tt6*tt49*X(1,1)*X(2,1)*W(25,1)*tt48)+4*tt1*tt11*X(1,1)&
&*X(2,1)*W(25,1)*tt48+2*tt1*tt8*X(1,1)*X(2,1)*tt10
tt66 = 2*tt49*tt50*tt51*X(1,1)*X(2,1)*W(25,1)*W(26,1)
tt67 = 2*tt49*tt53*tt54*X(1,1)*X(2,1)*W(25,1)*W(27,1)
tt68 = 2*tt49*tt56*tt57*X(1,1)*X(2,1)*W(25,1)*W(28,1)
tt69 = 2*tt49*tt59*tt60*X(1,1)*X(2,1)*W(25,1)*W(29,1)
tt70 = 2*tt49*tt62*tt63*X(1,1)*X(2,1)*W(25,1)*W(30,1)
tt71 = (-2*tt6*tt49*X(1,1)*X(3,1)*W(25,1)*tt48)+4*tt1*tt11*X(1,1)&
&*X(3,1)*W(25,1)*tt48+2*tt1*tt8*X(1,1)*X(3,1)*tt10
tt72 = 2*tt49*tt50*tt51*X(1,1)*X(3,1)*W(25,1)*W(26,1)
tt73 = 2*tt49*tt53*tt54*X(1,1)*X(3,1)*W(25,1)*W(27,1)
tt74 = 2*tt49*tt56*tt57*X(1,1)*X(3,1)*W(25,1)*W(28,1)
tt75 = 2*tt49*tt59*tt60*X(1,1)*X(3,1)*W(25,1)*W(29,1)
tt76 = 2*tt49*tt62*tt63*X(1,1)*X(3,1)*W(25,1)*W(30,1)
tt77 = (-2*tt6*tt49*X(1,1)*W(25,1)*tt48)+4*tt1*tt11*X(1,1)*W(25,1&
&)*tt48+2*tt1*tt8*X(1,1)*tt10
tt78 = 2*tt49*tt50*tt51*X(1,1)*W(25,1)*W(26,1)
tt79 = 2*tt49*tt53*tt54*X(1,1)*W(25,1)*W(27,1)
tt80 = 2*tt49*tt56*tt57*X(1,1)*W(25,1)*W(28,1)
tt81 = 2*tt49*tt59*tt60*X(1,1)*W(25,1)*W(29,1)
tt82 = 2*tt49*tt62*tt63*X(1,1)*W(25,1)*W(30,1)
tt83 = 2*tt6*tt49*X(1,1)*tt48+2*tt6*tt11*X(1,1)*W(25,1)
tt84 = 2*tt6*tt49*tt19*X(1,1)*W(25,1)
tt85 = 2*tt6*tt49*tt26*X(1,1)*W(25,1)
tt86 = 2*tt6*tt49*tt33*X(1,1)*W(25,1)
tt87 = 2*tt6*tt49*tt40*X(1,1)*W(25,1)
tt88 = 2*tt6*tt49*tt47*X(1,1)*W(25,1)
tt89 = 2*tt6*tt49*X(1,1)*W(25,1)
tt90 = exp((-2*W(20,1))-2*X(3,1)*W(14,1)-2*X(2,1)*W(8,1)-2*X(1,1)&
&*W(2,1))
tt91 = 1/tt18**4
tt92 = W(26,1)**2
tt93 = 1/tt18**3
tt94 = exp(tt23+tt16+tt22+tt15+tt21+tt14+tt20+tt13)
tt95 = 2*tt50*tt53*tt94*tt9*W(26,1)*W(27,1)
tt96 = exp(tt30+tt16+tt29+tt15+tt28+tt14+tt27+tt13)
tt97 = 2*tt50*tt56*tt96*tt9*W(26,1)*W(28,1)
tt98 = exp(tt37+tt16+tt36+tt15+tt35+tt14+tt34+tt13)
tt99 = 2*tt50*tt59*tt98*tt9*W(26,1)*W(29,1)
tt100 = exp(tt44+tt16+tt43+tt15+tt42+tt14+tt41+tt13)
tt101 = 2*tt50*tt62*tt100*tt9*W(26,1)*W(30,1)
tt102 = (-2*tt17*tt50*X(1,1)*X(2,1)*W(26,1)*tt48)+4*tt90*tt93*X(1&
&,1)*X(2,1)*W(26,1)*tt48+2*tt90*tt91*X(1,1)*X(2,1)*tt92
tt103 = 2*tt50*tt53*tt94*X(1,1)*X(2,1)*W(26,1)*W(27,1)
tt104 = 2*tt50*tt56*tt96*X(1,1)*X(2,1)*W(26,1)*W(28,1)
tt105 = 2*tt50*tt59*tt98*X(1,1)*X(2,1)*W(26,1)*W(29,1)
tt106 = 2*tt50*tt62*tt100*X(1,1)*X(2,1)*W(26,1)*W(30,1)
tt107 = (-2*tt17*tt50*X(1,1)*X(3,1)*W(26,1)*tt48)+4*tt90*tt93*X(1&
&,1)*X(3,1)*W(26,1)*tt48+2*tt90*tt91*X(1,1)*X(3,1)*tt92
tt108 = 2*tt50*tt53*tt94*X(1,1)*X(3,1)*W(26,1)*W(27,1)
tt109 = 2*tt50*tt56*tt96*X(1,1)*X(3,1)*W(26,1)*W(28,1)
tt110 = 2*tt50*tt59*tt98*X(1,1)*X(3,1)*W(26,1)*W(29,1)
tt111 = 2*tt50*tt62*tt100*X(1,1)*X(3,1)*W(26,1)*W(30,1)
tt112 = (-2*tt17*tt50*X(1,1)*W(26,1)*tt48)+4*tt90*tt93*X(1,1)*W(2&
&6,1)*tt48+2*tt90*tt91*X(1,1)*tt92
tt113 = 2*tt50*tt53*tt94*X(1,1)*W(26,1)*W(27,1)
tt114 = 2*tt50*tt56*tt96*X(1,1)*W(26,1)*W(28,1)
tt115 = 2*tt50*tt59*tt98*X(1,1)*W(26,1)*W(29,1)
tt116 = 2*tt50*tt62*tt100*X(1,1)*W(26,1)*W(30,1)
tt117 = 2*tt12*tt17*tt50*X(1,1)*W(26,1)
tt118 = 2*tt17*tt50*X(1,1)*tt48+2*tt17*tt93*X(1,1)*W(26,1)
tt119 = 2*tt17*tt50*tt26*X(1,1)*W(26,1)
tt120 = 2*tt17*tt50*tt33*X(1,1)*W(26,1)
tt121 = 2*tt17*tt50*tt40*X(1,1)*W(26,1)
tt122 = 2*tt17*tt50*tt47*X(1,1)*W(26,1)
tt123 = 2*tt17*tt50*X(1,1)*W(26,1)
tt124 = exp((-2*W(21,1))-2*X(3,1)*W(15,1)-2*X(2,1)*W(9,1)-2*X(1,1&
&)*W(3,1))
tt125 = 1/tt25**4
tt126 = W(27,1)**2
tt127 = 1/tt25**3
tt128 = exp(tt30+tt23+tt29+tt22+tt28+tt21+tt27+tt20)
tt129 = 2*tt53*tt56*tt128*tt9*W(27,1)*W(28,1)
tt130 = exp(tt37+tt23+tt36+tt22+tt35+tt21+tt34+tt20)
tt131 = 2*tt53*tt59*tt130*tt9*W(27,1)*W(29,1)
tt132 = exp(tt44+tt23+tt43+tt22+tt42+tt21+tt41+tt20)
tt133 = 2*tt53*tt62*tt132*tt9*W(27,1)*W(30,1)
tt134 = (-2*tt24*tt53*X(1,1)*X(2,1)*W(27,1)*tt48)+4*tt124*tt127*X&
&(1,1)*X(2,1)*W(27,1)*tt48+2*tt124*tt125*X(1,1)*X(2,1)*tt126
tt135 = 2*tt53*tt56*tt128*X(1,1)*X(2,1)*W(27,1)*W(28,1)
tt136 = 2*tt53*tt59*tt130*X(1,1)*X(2,1)*W(27,1)*W(29,1)
tt137 = 2*tt53*tt62*tt132*X(1,1)*X(2,1)*W(27,1)*W(30,1)
tt138 = (-2*tt24*tt53*X(1,1)*X(3,1)*W(27,1)*tt48)+4*tt124*tt127*X&
&(1,1)*X(3,1)*W(27,1)*tt48+2*tt124*tt125*X(1,1)*X(3,1)*tt126
tt139 = 2*tt53*tt56*tt128*X(1,1)*X(3,1)*W(27,1)*W(28,1)
tt140 = 2*tt53*tt59*tt130*X(1,1)*X(3,1)*W(27,1)*W(29,1)
tt141 = 2*tt53*tt62*tt132*X(1,1)*X(3,1)*W(27,1)*W(30,1)
tt142 = (-2*tt24*tt53*X(1,1)*W(27,1)*tt48)+4*tt124*tt127*X(1,1)*W&
&(27,1)*tt48+2*tt124*tt125*X(1,1)*tt126
tt143 = 2*tt53*tt56*tt128*X(1,1)*W(27,1)*W(28,1)
tt144 = 2*tt53*tt59*tt130*X(1,1)*W(27,1)*W(29,1)
tt145 = 2*tt53*tt62*tt132*X(1,1)*W(27,1)*W(30,1)
tt146 = 2*tt12*tt24*tt53*X(1,1)*W(27,1)
tt147 = 2*tt19*tt24*tt53*X(1,1)*W(27,1)
tt148 = 2*tt24*tt53*X(1,1)*tt48+2*tt24*tt127*X(1,1)*W(27,1)
tt149 = 2*tt24*tt53*tt33*X(1,1)*W(27,1)
tt150 = 2*tt24*tt53*tt40*X(1,1)*W(27,1)
tt151 = 2*tt24*tt53*tt47*X(1,1)*W(27,1)
tt152 = 2*tt24*tt53*X(1,1)*W(27,1)
tt153 = exp((-2*W(22,1))-2*X(3,1)*W(16,1)-2*X(2,1)*W(10,1)-2*X(1,&
&1)*W(4,1))
tt154 = 1/tt32**4
tt155 = W(28,1)**2
tt156 = 1/tt32**3
tt157 = exp(tt37+tt30+tt36+tt29+tt35+tt28+tt34+tt27)
tt158 = 2*tt56*tt59*tt157*tt9*W(28,1)*W(29,1)
tt159 = exp(tt44+tt30+tt43+tt29+tt42+tt28+tt41+tt27)
tt160 = 2*tt56*tt62*tt159*tt9*W(28,1)*W(30,1)
tt161 = (-2*tt31*tt56*X(1,1)*X(2,1)*W(28,1)*tt48)+4*tt153*tt156*X&
&(1,1)*X(2,1)*W(28,1)*tt48+2*tt153*tt154*X(1,1)*X(2,1)*tt155
tt162 = 2*tt56*tt59*tt157*X(1,1)*X(2,1)*W(28,1)*W(29,1)
tt163 = 2*tt56*tt62*tt159*X(1,1)*X(2,1)*W(28,1)*W(30,1)
tt164 = (-2*tt31*tt56*X(1,1)*X(3,1)*W(28,1)*tt48)+4*tt153*tt156*X&
&(1,1)*X(3,1)*W(28,1)*tt48+2*tt153*tt154*X(1,1)*X(3,1)*tt155
tt165 = 2*tt56*tt59*tt157*X(1,1)*X(3,1)*W(28,1)*W(29,1)
tt166 = 2*tt56*tt62*tt159*X(1,1)*X(3,1)*W(28,1)*W(30,1)
tt167 = (-2*tt31*tt56*X(1,1)*W(28,1)*tt48)+4*tt153*tt156*X(1,1)*W&
&(28,1)*tt48+2*tt153*tt154*X(1,1)*tt155
tt168 = 2*tt56*tt59*tt157*X(1,1)*W(28,1)*W(29,1)
tt169 = 2*tt56*tt62*tt159*X(1,1)*W(28,1)*W(30,1)
tt170 = 2*tt12*tt31*tt56*X(1,1)*W(28,1)
tt171 = 2*tt19*tt31*tt56*X(1,1)*W(28,1)
tt172 = 2*tt26*tt31*tt56*X(1,1)*W(28,1)
tt173 = 2*tt31*tt56*X(1,1)*tt48+2*tt31*tt156*X(1,1)*W(28,1)
tt174 = 2*tt31*tt56*tt40*X(1,1)*W(28,1)
tt175 = 2*tt31*tt56*tt47*X(1,1)*W(28,1)
tt176 = 2*tt31*tt56*X(1,1)*W(28,1)
tt177 = exp((-2*W(23,1))-2*X(3,1)*W(17,1)-2*X(2,1)*W(11,1)-2*X(1,&
&1)*W(5,1))
tt178 = 1/tt39**4
tt179 = W(29,1)**2
tt180 = 1/tt39**3
tt181 = exp(tt44+tt37+tt43+tt36+tt42+tt35+tt41+tt34)
tt182 = 2*tt59*tt62*tt181*tt9*W(29,1)*W(30,1)
tt183 = (-2*tt38*tt59*X(1,1)*X(2,1)*W(29,1)*tt48)+4*tt177*tt180*X&
&(1,1)*X(2,1)*W(29,1)*tt48+2*tt177*tt178*X(1,1)*X(2,1)*tt179
tt184 = 2*tt59*tt62*tt181*X(1,1)*X(2,1)*W(29,1)*W(30,1)
tt185 = (-2*tt38*tt59*X(1,1)*X(3,1)*W(29,1)*tt48)+4*tt177*tt180*X&
&(1,1)*X(3,1)*W(29,1)*tt48+2*tt177*tt178*X(1,1)*X(3,1)*tt179
tt186 = 2*tt59*tt62*tt181*X(1,1)*X(3,1)*W(29,1)*W(30,1)
tt187 = (-2*tt38*tt59*X(1,1)*W(29,1)*tt48)+4*tt177*tt180*X(1,1)*W&
&(29,1)*tt48+2*tt177*tt178*X(1,1)*tt179
tt188 = 2*tt59*tt62*tt181*X(1,1)*W(29,1)*W(30,1)
tt189 = 2*tt12*tt38*tt59*X(1,1)*W(29,1)
tt190 = 2*tt19*tt38*tt59*X(1,1)*W(29,1)
tt191 = 2*tt26*tt38*tt59*X(1,1)*W(29,1)
tt192 = 2*tt33*tt38*tt59*X(1,1)*W(29,1)
tt193 = 2*tt38*tt59*X(1,1)*tt48+2*tt38*tt180*X(1,1)*W(29,1)
tt194 = 2*tt38*tt59*tt47*X(1,1)*W(29,1)
tt195 = 2*tt38*tt59*X(1,1)*W(29,1)
tt196 = exp((-2*W(24,1))-2*X(3,1)*W(18,1)-2*X(2,1)*W(12,1)-2*X(1,&
&1)*W(6,1))
tt197 = 1/tt46**4
tt198 = W(30,1)**2
tt199 = 1/tt46**3
tt200 = (-2*tt45*tt62*X(1,1)*X(2,1)*W(30,1)*tt48)+4*tt196*tt199*X&
&(1,1)*X(2,1)*W(30,1)*tt48+2*tt196*tt197*X(1,1)*X(2,1)*tt198
tt201 = (-2*tt45*tt62*X(1,1)*X(3,1)*W(30,1)*tt48)+4*tt196*tt199*X&
&(1,1)*X(3,1)*W(30,1)*tt48+2*tt196*tt197*X(1,1)*X(3,1)*tt198
tt202 = (-2*tt45*tt62*X(1,1)*W(30,1)*tt48)+4*tt196*tt199*X(1,1)*W&
&(30,1)*tt48+2*tt196*tt197*X(1,1)*tt198
tt203 = 2*tt12*tt45*tt62*X(1,1)*W(30,1)
tt204 = 2*tt19*tt45*tt62*X(1,1)*W(30,1)
tt205 = 2*tt26*tt45*tt62*X(1,1)*W(30,1)
tt206 = 2*tt33*tt45*tt62*X(1,1)*W(30,1)
tt207 = 2*tt40*tt45*tt62*X(1,1)*W(30,1)
tt208 = 2*tt45*tt62*X(1,1)*tt48+2*tt45*tt199*X(1,1)*W(30,1)
tt209 = 2*tt45*tt62*X(1,1)*W(30,1)
tt210 = X(2,1)**2
tt211 = 2*tt49*tt50*tt51*tt210*W(25,1)*W(26,1)
tt212 = 2*tt49*tt53*tt54*tt210*W(25,1)*W(27,1)
tt213 = 2*tt49*tt56*tt57*tt210*W(25,1)*W(28,1)
tt214 = 2*tt49*tt59*tt60*tt210*W(25,1)*W(29,1)
tt215 = 2*tt49*tt62*tt63*tt210*W(25,1)*W(30,1)
tt216 = (-2*tt6*tt49*X(2,1)*X(3,1)*W(25,1)*tt48)+4*tt1*tt11*X(2,1&
&)*X(3,1)*W(25,1)*tt48+2*tt1*tt8*X(2,1)*X(3,1)*tt10
tt217 = 2*tt49*tt50*tt51*X(2,1)*X(3,1)*W(25,1)*W(26,1)
tt218 = 2*tt49*tt53*tt54*X(2,1)*X(3,1)*W(25,1)*W(27,1)
tt219 = 2*tt49*tt56*tt57*X(2,1)*X(3,1)*W(25,1)*W(28,1)
tt220 = 2*tt49*tt59*tt60*X(2,1)*X(3,1)*W(25,1)*W(29,1)
tt221 = 2*tt49*tt62*tt63*X(2,1)*X(3,1)*W(25,1)*W(30,1)
tt222 = (-2*tt6*tt49*X(2,1)*W(25,1)*tt48)+4*tt1*tt11*X(2,1)*W(25,&
&1)*tt48+2*tt1*tt8*X(2,1)*tt10
tt223 = 2*tt49*tt50*tt51*X(2,1)*W(25,1)*W(26,1)
tt224 = 2*tt49*tt53*tt54*X(2,1)*W(25,1)*W(27,1)
tt225 = 2*tt49*tt56*tt57*X(2,1)*W(25,1)*W(28,1)
tt226 = 2*tt49*tt59*tt60*X(2,1)*W(25,1)*W(29,1)
tt227 = 2*tt49*tt62*tt63*X(2,1)*W(25,1)*W(30,1)
tt228 = 2*tt6*tt49*X(2,1)*tt48+2*tt6*tt11*X(2,1)*W(25,1)
tt229 = 2*tt6*tt49*tt19*X(2,1)*W(25,1)
tt230 = 2*tt6*tt49*tt26*X(2,1)*W(25,1)
tt231 = 2*tt6*tt49*tt33*X(2,1)*W(25,1)
tt232 = 2*tt6*tt49*tt40*X(2,1)*W(25,1)
tt233 = 2*tt6*tt49*tt47*X(2,1)*W(25,1)
tt234 = 2*tt6*tt49*X(2,1)*W(25,1)
tt235 = 2*tt50*tt53*tt94*tt210*W(26,1)*W(27,1)
tt236 = 2*tt50*tt56*tt96*tt210*W(26,1)*W(28,1)
tt237 = 2*tt50*tt59*tt98*tt210*W(26,1)*W(29,1)
tt238 = 2*tt50*tt62*tt100*tt210*W(26,1)*W(30,1)
tt239 = (-2*tt17*tt50*X(2,1)*X(3,1)*W(26,1)*tt48)+4*tt90*tt93*X(2&
&,1)*X(3,1)*W(26,1)*tt48+2*tt90*tt91*X(2,1)*X(3,1)*tt92
tt240 = 2*tt50*tt53*tt94*X(2,1)*X(3,1)*W(26,1)*W(27,1)
tt241 = 2*tt50*tt56*tt96*X(2,1)*X(3,1)*W(26,1)*W(28,1)
tt242 = 2*tt50*tt59*tt98*X(2,1)*X(3,1)*W(26,1)*W(29,1)
tt243 = 2*tt50*tt62*tt100*X(2,1)*X(3,1)*W(26,1)*W(30,1)
tt244 = (-2*tt17*tt50*X(2,1)*W(26,1)*tt48)+4*tt90*tt93*X(2,1)*W(2&
&6,1)*tt48+2*tt90*tt91*X(2,1)*tt92
tt245 = 2*tt50*tt53*tt94*X(2,1)*W(26,1)*W(27,1)
tt246 = 2*tt50*tt56*tt96*X(2,1)*W(26,1)*W(28,1)
tt247 = 2*tt50*tt59*tt98*X(2,1)*W(26,1)*W(29,1)
tt248 = 2*tt50*tt62*tt100*X(2,1)*W(26,1)*W(30,1)
tt249 = 2*tt12*tt17*tt50*X(2,1)*W(26,1)
tt250 = 2*tt17*tt50*X(2,1)*tt48+2*tt17*tt93*X(2,1)*W(26,1)
tt251 = 2*tt17*tt50*tt26*X(2,1)*W(26,1)
tt252 = 2*tt17*tt50*tt33*X(2,1)*W(26,1)
tt253 = 2*tt17*tt50*tt40*X(2,1)*W(26,1)
tt254 = 2*tt17*tt50*tt47*X(2,1)*W(26,1)
tt255 = 2*tt17*tt50*X(2,1)*W(26,1)
tt256 = 2*tt53*tt56*tt128*tt210*W(27,1)*W(28,1)
tt257 = 2*tt53*tt59*tt130*tt210*W(27,1)*W(29,1)
tt258 = 2*tt53*tt62*tt132*tt210*W(27,1)*W(30,1)
tt259 = (-2*tt24*tt53*X(2,1)*X(3,1)*W(27,1)*tt48)+4*tt124*tt127*X&
&(2,1)*X(3,1)*W(27,1)*tt48+2*tt124*tt125*X(2,1)*X(3,1)*tt126
tt260 = 2*tt53*tt56*tt128*X(2,1)*X(3,1)*W(27,1)*W(28,1)
tt261 = 2*tt53*tt59*tt130*X(2,1)*X(3,1)*W(27,1)*W(29,1)
tt262 = 2*tt53*tt62*tt132*X(2,1)*X(3,1)*W(27,1)*W(30,1)
tt263 = (-2*tt24*tt53*X(2,1)*W(27,1)*tt48)+4*tt124*tt127*X(2,1)*W&
&(27,1)*tt48+2*tt124*tt125*X(2,1)*tt126
tt264 = 2*tt53*tt56*tt128*X(2,1)*W(27,1)*W(28,1)
tt265 = 2*tt53*tt59*tt130*X(2,1)*W(27,1)*W(29,1)
tt266 = 2*tt53*tt62*tt132*X(2,1)*W(27,1)*W(30,1)
tt267 = 2*tt12*tt24*tt53*X(2,1)*W(27,1)
tt268 = 2*tt19*tt24*tt53*X(2,1)*W(27,1)
tt269 = 2*tt24*tt53*X(2,1)*tt48+2*tt24*tt127*X(2,1)*W(27,1)
tt270 = 2*tt24*tt53*tt33*X(2,1)*W(27,1)
tt271 = 2*tt24*tt53*tt40*X(2,1)*W(27,1)
tt272 = 2*tt24*tt53*tt47*X(2,1)*W(27,1)
tt273 = 2*tt24*tt53*X(2,1)*W(27,1)
tt274 = 2*tt56*tt59*tt157*tt210*W(28,1)*W(29,1)
tt275 = 2*tt56*tt62*tt159*tt210*W(28,1)*W(30,1)
tt276 = (-2*tt31*tt56*X(2,1)*X(3,1)*W(28,1)*tt48)+4*tt153*tt156*X&
&(2,1)*X(3,1)*W(28,1)*tt48+2*tt153*tt154*X(2,1)*X(3,1)*tt155
tt277 = 2*tt56*tt59*tt157*X(2,1)*X(3,1)*W(28,1)*W(29,1)
tt278 = 2*tt56*tt62*tt159*X(2,1)*X(3,1)*W(28,1)*W(30,1)
tt279 = (-2*tt31*tt56*X(2,1)*W(28,1)*tt48)+4*tt153*tt156*X(2,1)*W&
&(28,1)*tt48+2*tt153*tt154*X(2,1)*tt155
tt280 = 2*tt56*tt59*tt157*X(2,1)*W(28,1)*W(29,1)
tt281 = 2*tt56*tt62*tt159*X(2,1)*W(28,1)*W(30,1)
tt282 = 2*tt12*tt31*tt56*X(2,1)*W(28,1)
tt283 = 2*tt19*tt31*tt56*X(2,1)*W(28,1)
tt284 = 2*tt26*tt31*tt56*X(2,1)*W(28,1)
tt285 = 2*tt31*tt56*X(2,1)*tt48+2*tt31*tt156*X(2,1)*W(28,1)
tt286 = 2*tt31*tt56*tt40*X(2,1)*W(28,1)
tt287 = 2*tt31*tt56*tt47*X(2,1)*W(28,1)
tt288 = 2*tt31*tt56*X(2,1)*W(28,1)
tt289 = 2*tt59*tt62*tt181*tt210*W(29,1)*W(30,1)
tt290 = (-2*tt38*tt59*X(2,1)*X(3,1)*W(29,1)*tt48)+4*tt177*tt180*X&
&(2,1)*X(3,1)*W(29,1)*tt48+2*tt177*tt178*X(2,1)*X(3,1)*tt179
tt291 = 2*tt59*tt62*tt181*X(2,1)*X(3,1)*W(29,1)*W(30,1)
tt292 = (-2*tt38*tt59*X(2,1)*W(29,1)*tt48)+4*tt177*tt180*X(2,1)*W&
&(29,1)*tt48+2*tt177*tt178*X(2,1)*tt179
tt293 = 2*tt59*tt62*tt181*X(2,1)*W(29,1)*W(30,1)
tt294 = 2*tt12*tt38*tt59*X(2,1)*W(29,1)
tt295 = 2*tt19*tt38*tt59*X(2,1)*W(29,1)
tt296 = 2*tt26*tt38*tt59*X(2,1)*W(29,1)
tt297 = 2*tt33*tt38*tt59*X(2,1)*W(29,1)
tt298 = 2*tt38*tt59*X(2,1)*tt48+2*tt38*tt180*X(2,1)*W(29,1)
tt299 = 2*tt38*tt59*tt47*X(2,1)*W(29,1)
tt300 = 2*tt38*tt59*X(2,1)*W(29,1)
tt301 = (-2*tt45*tt62*X(2,1)*X(3,1)*W(30,1)*tt48)+4*tt196*tt199*X&
&(2,1)*X(3,1)*W(30,1)*tt48+2*tt196*tt197*X(2,1)*X(3,1)*tt198
tt302 = (-2*tt45*tt62*X(2,1)*W(30,1)*tt48)+4*tt196*tt199*X(2,1)*W&
&(30,1)*tt48+2*tt196*tt197*X(2,1)*tt198
tt303 = 2*tt12*tt45*tt62*X(2,1)*W(30,1)
tt304 = 2*tt19*tt45*tt62*X(2,1)*W(30,1)
tt305 = 2*tt26*tt45*tt62*X(2,1)*W(30,1)
tt306 = 2*tt33*tt45*tt62*X(2,1)*W(30,1)
tt307 = 2*tt40*tt45*tt62*X(2,1)*W(30,1)
tt308 = 2*tt45*tt62*X(2,1)*tt48+2*tt45*tt199*X(2,1)*W(30,1)
tt309 = 2*tt45*tt62*X(2,1)*W(30,1)
tt310 = X(3,1)**2
tt311 = 2*tt49*tt50*tt51*tt310*W(25,1)*W(26,1)
tt312 = 2*tt49*tt53*tt54*tt310*W(25,1)*W(27,1)
tt313 = 2*tt49*tt56*tt57*tt310*W(25,1)*W(28,1)
tt314 = 2*tt49*tt59*tt60*tt310*W(25,1)*W(29,1)
tt315 = 2*tt49*tt62*tt63*tt310*W(25,1)*W(30,1)
tt316 = (-2*tt6*tt49*X(3,1)*W(25,1)*tt48)+4*tt1*tt11*X(3,1)*W(25,&
&1)*tt48+2*tt1*tt8*X(3,1)*tt10
tt317 = 2*tt49*tt50*tt51*X(3,1)*W(25,1)*W(26,1)
tt318 = 2*tt49*tt53*tt54*X(3,1)*W(25,1)*W(27,1)
tt319 = 2*tt49*tt56*tt57*X(3,1)*W(25,1)*W(28,1)
tt320 = 2*tt49*tt59*tt60*X(3,1)*W(25,1)*W(29,1)
tt321 = 2*tt49*tt62*tt63*X(3,1)*W(25,1)*W(30,1)
tt322 = 2*tt6*tt49*X(3,1)*tt48+2*tt6*tt11*X(3,1)*W(25,1)
tt323 = 2*tt6*tt49*tt19*X(3,1)*W(25,1)
tt324 = 2*tt6*tt49*tt26*X(3,1)*W(25,1)
tt325 = 2*tt6*tt49*tt33*X(3,1)*W(25,1)
tt326 = 2*tt6*tt49*tt40*X(3,1)*W(25,1)
tt327 = 2*tt6*tt49*tt47*X(3,1)*W(25,1)
tt328 = 2*tt6*tt49*X(3,1)*W(25,1)
tt329 = 2*tt50*tt53*tt94*tt310*W(26,1)*W(27,1)
tt330 = 2*tt50*tt56*tt96*tt310*W(26,1)*W(28,1)
tt331 = 2*tt50*tt59*tt98*tt310*W(26,1)*W(29,1)
tt332 = 2*tt50*tt62*tt100*tt310*W(26,1)*W(30,1)
tt333 = (-2*tt17*tt50*X(3,1)*W(26,1)*tt48)+4*tt90*tt93*X(3,1)*W(2&
&6,1)*tt48+2*tt90*tt91*X(3,1)*tt92
tt334 = 2*tt50*tt53*tt94*X(3,1)*W(26,1)*W(27,1)
tt335 = 2*tt50*tt56*tt96*X(3,1)*W(26,1)*W(28,1)
tt336 = 2*tt50*tt59*tt98*X(3,1)*W(26,1)*W(29,1)
tt337 = 2*tt50*tt62*tt100*X(3,1)*W(26,1)*W(30,1)
tt338 = 2*tt12*tt17*tt50*X(3,1)*W(26,1)
tt339 = 2*tt17*tt50*X(3,1)*tt48+2*tt17*tt93*X(3,1)*W(26,1)
tt340 = 2*tt17*tt50*tt26*X(3,1)*W(26,1)
tt341 = 2*tt17*tt50*tt33*X(3,1)*W(26,1)
tt342 = 2*tt17*tt50*tt40*X(3,1)*W(26,1)
tt343 = 2*tt17*tt50*tt47*X(3,1)*W(26,1)
tt344 = 2*tt17*tt50*X(3,1)*W(26,1)
tt345 = 2*tt53*tt56*tt128*tt310*W(27,1)*W(28,1)
tt346 = 2*tt53*tt59*tt130*tt310*W(27,1)*W(29,1)
tt347 = 2*tt53*tt62*tt132*tt310*W(27,1)*W(30,1)
tt348 = (-2*tt24*tt53*X(3,1)*W(27,1)*tt48)+4*tt124*tt127*X(3,1)*W&
&(27,1)*tt48+2*tt124*tt125*X(3,1)*tt126
tt349 = 2*tt53*tt56*tt128*X(3,1)*W(27,1)*W(28,1)
tt350 = 2*tt53*tt59*tt130*X(3,1)*W(27,1)*W(29,1)
tt351 = 2*tt53*tt62*tt132*X(3,1)*W(27,1)*W(30,1)
tt352 = 2*tt12*tt24*tt53*X(3,1)*W(27,1)
tt353 = 2*tt19*tt24*tt53*X(3,1)*W(27,1)
tt354 = 2*tt24*tt53*X(3,1)*tt48+2*tt24*tt127*X(3,1)*W(27,1)
tt355 = 2*tt24*tt53*tt33*X(3,1)*W(27,1)
tt356 = 2*tt24*tt53*tt40*X(3,1)*W(27,1)
tt357 = 2*tt24*tt53*tt47*X(3,1)*W(27,1)
tt358 = 2*tt24*tt53*X(3,1)*W(27,1)
tt359 = 2*tt56*tt59*tt157*tt310*W(28,1)*W(29,1)
tt360 = 2*tt56*tt62*tt159*tt310*W(28,1)*W(30,1)
tt361 = (-2*tt31*tt56*X(3,1)*W(28,1)*tt48)+4*tt153*tt156*X(3,1)*W&
&(28,1)*tt48+2*tt153*tt154*X(3,1)*tt155
tt362 = 2*tt56*tt59*tt157*X(3,1)*W(28,1)*W(29,1)
tt363 = 2*tt56*tt62*tt159*X(3,1)*W(28,1)*W(30,1)
tt364 = 2*tt12*tt31*tt56*X(3,1)*W(28,1)
tt365 = 2*tt19*tt31*tt56*X(3,1)*W(28,1)
tt366 = 2*tt26*tt31*tt56*X(3,1)*W(28,1)
tt367 = 2*tt31*tt56*X(3,1)*tt48+2*tt31*tt156*X(3,1)*W(28,1)
tt368 = 2*tt31*tt56*tt40*X(3,1)*W(28,1)
tt369 = 2*tt31*tt56*tt47*X(3,1)*W(28,1)
tt370 = 2*tt31*tt56*X(3,1)*W(28,1)
tt371 = 2*tt59*tt62*tt181*tt310*W(29,1)*W(30,1)
tt372 = (-2*tt38*tt59*X(3,1)*W(29,1)*tt48)+4*tt177*tt180*X(3,1)*W&
&(29,1)*tt48+2*tt177*tt178*X(3,1)*tt179
tt373 = 2*tt59*tt62*tt181*X(3,1)*W(29,1)*W(30,1)
tt374 = 2*tt12*tt38*tt59*X(3,1)*W(29,1)
tt375 = 2*tt19*tt38*tt59*X(3,1)*W(29,1)
tt376 = 2*tt26*tt38*tt59*X(3,1)*W(29,1)
tt377 = 2*tt33*tt38*tt59*X(3,1)*W(29,1)
tt378 = 2*tt38*tt59*X(3,1)*tt48+2*tt38*tt180*X(3,1)*W(29,1)
tt379 = 2*tt38*tt59*tt47*X(3,1)*W(29,1)
tt380 = 2*tt38*tt59*X(3,1)*W(29,1)
tt381 = (-2*tt45*tt62*X(3,1)*W(30,1)*tt48)+4*tt196*tt199*X(3,1)*W&
&(30,1)*tt48+2*tt196*tt197*X(3,1)*tt198
tt382 = 2*tt12*tt45*tt62*X(3,1)*W(30,1)
tt383 = 2*tt19*tt45*tt62*X(3,1)*W(30,1)
tt384 = 2*tt26*tt45*tt62*X(3,1)*W(30,1)
tt385 = 2*tt33*tt45*tt62*X(3,1)*W(30,1)
tt386 = 2*tt40*tt45*tt62*X(3,1)*W(30,1)
tt387 = 2*tt45*tt62*X(3,1)*tt48+2*tt45*tt199*X(3,1)*W(30,1)
tt388 = 2*tt45*tt62*X(3,1)*W(30,1)
tt389 = 2*tt49*tt50*tt51*W(25,1)*W(26,1)
tt390 = 2*tt49*tt53*tt54*W(25,1)*W(27,1)
tt391 = 2*tt49*tt56*tt57*W(25,1)*W(28,1)
tt392 = 2*tt49*tt59*tt60*W(25,1)*W(29,1)
tt393 = 2*tt49*tt62*tt63*W(25,1)*W(30,1)
tt394 = 2*tt6*tt49*tt48+2*tt6*tt11*W(25,1)
tt395 = 2*tt6*tt49*tt19*W(25,1)
tt396 = 2*tt6*tt49*tt26*W(25,1)
tt397 = 2*tt6*tt49*tt33*W(25,1)
tt398 = 2*tt6*tt49*tt40*W(25,1)
tt399 = 2*tt6*tt49*tt47*W(25,1)
tt400 = 2*tt6*tt49*W(25,1)
tt401 = 2*tt50*tt53*tt94*W(26,1)*W(27,1)
tt402 = 2*tt50*tt56*tt96*W(26,1)*W(28,1)
tt403 = 2*tt50*tt59*tt98*W(26,1)*W(29,1)
tt404 = 2*tt50*tt62*tt100*W(26,1)*W(30,1)
tt405 = 2*tt12*tt17*tt50*W(26,1)
tt406 = 2*tt17*tt50*tt48+2*tt17*tt93*W(26,1)
tt407 = 2*tt17*tt50*tt26*W(26,1)
tt408 = 2*tt17*tt50*tt33*W(26,1)
tt409 = 2*tt17*tt50*tt40*W(26,1)
tt410 = 2*tt17*tt50*tt47*W(26,1)
tt411 = 2*tt17*tt50*W(26,1)
tt412 = 2*tt53*tt56*tt128*W(27,1)*W(28,1)
tt413 = 2*tt53*tt59*tt130*W(27,1)*W(29,1)
tt414 = 2*tt53*tt62*tt132*W(27,1)*W(30,1)
tt415 = 2*tt12*tt24*tt53*W(27,1)
tt416 = 2*tt19*tt24*tt53*W(27,1)
tt417 = 2*tt24*tt53*tt48+2*tt24*tt127*W(27,1)
tt418 = 2*tt24*tt53*tt33*W(27,1)
tt419 = 2*tt24*tt53*tt40*W(27,1)
tt420 = 2*tt24*tt53*tt47*W(27,1)
tt421 = 2*tt24*tt53*W(27,1)
tt422 = 2*tt56*tt59*tt157*W(28,1)*W(29,1)
tt423 = 2*tt56*tt62*tt159*W(28,1)*W(30,1)
tt424 = 2*tt12*tt31*tt56*W(28,1)
tt425 = 2*tt19*tt31*tt56*W(28,1)
tt426 = 2*tt26*tt31*tt56*W(28,1)
tt427 = 2*tt31*tt56*tt48+2*tt31*tt156*W(28,1)
tt428 = 2*tt31*tt56*tt40*W(28,1)
tt429 = 2*tt31*tt56*tt47*W(28,1)
tt430 = 2*tt31*tt56*W(28,1)
tt431 = 2*tt59*tt62*tt181*W(29,1)*W(30,1)
tt432 = 2*tt12*tt38*tt59*W(29,1)
tt433 = 2*tt19*tt38*tt59*W(29,1)
tt434 = 2*tt26*tt38*tt59*W(29,1)
tt435 = 2*tt33*tt38*tt59*W(29,1)
tt436 = 2*tt38*tt59*tt48+2*tt38*tt180*W(29,1)
tt437 = 2*tt38*tt59*tt47*W(29,1)
tt438 = 2*tt38*tt59*W(29,1)
tt439 = 2*tt12*tt45*tt62*W(30,1)
tt440 = 2*tt19*tt45*tt62*W(30,1)
tt441 = 2*tt26*tt45*tt62*W(30,1)
tt442 = 2*tt33*tt45*tt62*W(30,1)
tt443 = 2*tt40*tt45*tt62*W(30,1)
tt444 = 2*tt45*tt62*tt48+2*tt45*tt199*W(30,1)
tt445 = 2*tt45*tt62*W(30,1)
tt446 = 2*tt12*tt19
tt447 = 2*tt12*tt26
tt448 = 2*tt12*tt33
tt449 = 2*tt12*tt40
tt450 = 2*tt12*tt47
tt451 = 2*tt12
tt452 = 2*tt19*tt26
tt453 = 2*tt19*tt33
tt454 = 2*tt19*tt40
tt455 = 2*tt19*tt47
tt456 = 2*tt19
tt457 = 2*tt26*tt33
tt458 = 2*tt26*tt40
tt459 = 2*tt26*tt47
tt460 = 2*tt26
tt461 = 2*tt33*tt40
tt462 = 2*tt33*tt47
tt463 = 2*tt33
tt464 = 2*tt40*tt47
tt465 = 2*tt40
tt466 = 2*tt47
hes(1,1) = (-2*tt6*tt49*tt9*W(25,1)*tt48)+4*tt1*tt11*tt9*W(25,1)*&
&tt48+2*tt1*tt8*tt9*tt10
hes(1,2) = tt52
hes(1,3) = tt55
hes(1,4) = tt58
hes(1,5) = tt61
hes(1,6) = tt64
hes(1,7) = tt65
hes(1,8) = tt66
hes(1,9) = tt67
hes(1,10) = tt68
hes(1,11) = tt69
hes(1,12) = tt70
hes(1,13) = tt71
hes(1,14) = tt72
hes(1,15) = tt73
hes(1,16) = tt74
hes(1,17) = tt75
hes(1,18) = tt76
hes(1,19) = tt77
hes(1,20) = tt78
hes(1,21) = tt79
hes(1,22) = tt80
hes(1,23) = tt81
hes(1,24) = tt82
hes(1,25) = tt83
hes(1,26) = tt84
hes(1,27) = tt85
hes(1,28) = tt86
hes(1,29) = tt87
hes(1,30) = tt88
hes(1,31) = tt89
hes(2,1) = tt52
hes(2,2) = (-2*tt17*tt50*tt9*W(26,1)*tt48)+4*tt90*tt93*tt9*W(26,1&
&)*tt48+2*tt90*tt91*tt9*tt92
hes(2,3) = tt95
hes(2,4) = tt97
hes(2,5) = tt99
hes(2,6) = tt101
hes(2,7) = tt66
hes(2,8) = tt102
hes(2,9) = tt103
hes(2,10) = tt104
hes(2,11) = tt105
hes(2,12) = tt106
hes(2,13) = tt72
hes(2,14) = tt107
hes(2,15) = tt108
hes(2,16) = tt109
hes(2,17) = tt110
hes(2,18) = tt111
hes(2,19) = tt78
hes(2,20) = tt112
hes(2,21) = tt113
hes(2,22) = tt114
hes(2,23) = tt115
hes(2,24) = tt116
hes(2,25) = tt117
hes(2,26) = tt118
hes(2,27) = tt119
hes(2,28) = tt120
hes(2,29) = tt121
hes(2,30) = tt122
hes(2,31) = tt123
hes(3,1) = tt55
hes(3,2) = tt95
hes(3,3) = (-2*tt24*tt53*tt9*W(27,1)*tt48)+4*tt124*tt127*tt9*W(27&
&,1)*tt48+2*tt124*tt125*tt9*tt126
hes(3,4) = tt129
hes(3,5) = tt131
hes(3,6) = tt133
hes(3,7) = tt67
hes(3,8) = tt103
hes(3,9) = tt134
hes(3,10) = tt135
hes(3,11) = tt136
hes(3,12) = tt137
hes(3,13) = tt73
hes(3,14) = tt108
hes(3,15) = tt138
hes(3,16) = tt139
hes(3,17) = tt140
hes(3,18) = tt141
hes(3,19) = tt79
hes(3,20) = tt113
hes(3,21) = tt142
hes(3,22) = tt143
hes(3,23) = tt144
hes(3,24) = tt145
hes(3,25) = tt146
hes(3,26) = tt147
hes(3,27) = tt148
hes(3,28) = tt149
hes(3,29) = tt150
hes(3,30) = tt151
hes(3,31) = tt152
hes(4,1) = tt58
hes(4,2) = tt97
hes(4,3) = tt129
hes(4,4) = (-2*tt31*tt56*tt9*W(28,1)*tt48)+4*tt153*tt156*tt9*W(28&
&,1)*tt48+2*tt153*tt154*tt9*tt155
hes(4,5) = tt158
hes(4,6) = tt160
hes(4,7) = tt68
hes(4,8) = tt104
hes(4,9) = tt135
hes(4,10) = tt161
hes(4,11) = tt162
hes(4,12) = tt163
hes(4,13) = tt74
hes(4,14) = tt109
hes(4,15) = tt139
hes(4,16) = tt164
hes(4,17) = tt165
hes(4,18) = tt166
hes(4,19) = tt80
hes(4,20) = tt114
hes(4,21) = tt143
hes(4,22) = tt167
hes(4,23) = tt168
hes(4,24) = tt169
hes(4,25) = tt170
hes(4,26) = tt171
hes(4,27) = tt172
hes(4,28) = tt173
hes(4,29) = tt174
hes(4,30) = tt175
hes(4,31) = tt176
hes(5,1) = tt61
hes(5,2) = tt99
hes(5,3) = tt131
hes(5,4) = tt158
hes(5,5) = (-2*tt38*tt59*tt9*W(29,1)*tt48)+4*tt177*tt180*tt9*W(29&
&,1)*tt48+2*tt177*tt178*tt9*tt179
hes(5,6) = tt182
hes(5,7) = tt69
hes(5,8) = tt105
hes(5,9) = tt136
hes(5,10) = tt162
hes(5,11) = tt183
hes(5,12) = tt184
hes(5,13) = tt75
hes(5,14) = tt110
hes(5,15) = tt140
hes(5,16) = tt165
hes(5,17) = tt185
hes(5,18) = tt186
hes(5,19) = tt81
hes(5,20) = tt115
hes(5,21) = tt144
hes(5,22) = tt168
hes(5,23) = tt187
hes(5,24) = tt188
hes(5,25) = tt189
hes(5,26) = tt190
hes(5,27) = tt191
hes(5,28) = tt192
hes(5,29) = tt193
hes(5,30) = tt194
hes(5,31) = tt195
hes(6,1) = tt64
hes(6,2) = tt101
hes(6,3) = tt133
hes(6,4) = tt160
hes(6,5) = tt182
hes(6,6) = (-2*tt45*tt62*tt9*W(30,1)*tt48)+4*tt196*tt199*tt9*W(30&
&,1)*tt48+2*tt196*tt197*tt9*tt198
hes(6,7) = tt70
hes(6,8) = tt106
hes(6,9) = tt137
hes(6,10) = tt163
hes(6,11) = tt184
hes(6,12) = tt200
hes(6,13) = tt76
hes(6,14) = tt111
hes(6,15) = tt141
hes(6,16) = tt166
hes(6,17) = tt186
hes(6,18) = tt201
hes(6,19) = tt82
hes(6,20) = tt116
hes(6,21) = tt145
hes(6,22) = tt169
hes(6,23) = tt188
hes(6,24) = tt202
hes(6,25) = tt203
hes(6,26) = tt204
hes(6,27) = tt205
hes(6,28) = tt206
hes(6,29) = tt207
hes(6,30) = tt208
hes(6,31) = tt209
hes(7,1) = tt65
hes(7,2) = tt66
hes(7,3) = tt67
hes(7,4) = tt68
hes(7,5) = tt69
hes(7,6) = tt70
hes(7,7) = (-2*tt6*tt49*tt210*W(25,1)*tt48)+4*tt1*tt11*tt210*W(25&
&,1)*tt48+2*tt1*tt8*tt210*tt10
hes(7,8) = tt211
hes(7,9) = tt212
hes(7,10) = tt213
hes(7,11) = tt214
hes(7,12) = tt215
hes(7,13) = tt216
hes(7,14) = tt217
hes(7,15) = tt218
hes(7,16) = tt219
hes(7,17) = tt220
hes(7,18) = tt221
hes(7,19) = tt222
hes(7,20) = tt223
hes(7,21) = tt224
hes(7,22) = tt225
hes(7,23) = tt226
hes(7,24) = tt227
hes(7,25) = tt228
hes(7,26) = tt229
hes(7,27) = tt230
hes(7,28) = tt231
hes(7,29) = tt232
hes(7,30) = tt233
hes(7,31) = tt234
hes(8,1) = tt66
hes(8,2) = tt102
hes(8,3) = tt103
hes(8,4) = tt104
hes(8,5) = tt105
hes(8,6) = tt106
hes(8,7) = tt211
hes(8,8) = (-2*tt17*tt50*tt210*W(26,1)*tt48)+4*tt90*tt93*tt210*W(&
&26,1)*tt48+2*tt90*tt91*tt210*tt92
hes(8,9) = tt235
hes(8,10) = tt236
hes(8,11) = tt237
hes(8,12) = tt238
hes(8,13) = tt217
hes(8,14) = tt239
hes(8,15) = tt240
hes(8,16) = tt241
hes(8,17) = tt242
hes(8,18) = tt243
hes(8,19) = tt223
hes(8,20) = tt244
hes(8,21) = tt245
hes(8,22) = tt246
hes(8,23) = tt247
hes(8,24) = tt248
hes(8,25) = tt249
hes(8,26) = tt250
hes(8,27) = tt251
hes(8,28) = tt252
hes(8,29) = tt253
hes(8,30) = tt254
hes(8,31) = tt255
hes(9,1) = tt67
hes(9,2) = tt103
hes(9,3) = tt134
hes(9,4) = tt135
hes(9,5) = tt136
hes(9,6) = tt137
hes(9,7) = tt212
hes(9,8) = tt235
hes(9,9) = (-2*tt24*tt53*tt210*W(27,1)*tt48)+4*tt124*tt127*tt210*&
&W(27,1)*tt48+2*tt124*tt125*tt210*tt126
hes(9,10) = tt256
hes(9,11) = tt257
hes(9,12) = tt258
hes(9,13) = tt218
hes(9,14) = tt240
hes(9,15) = tt259
hes(9,16) = tt260
hes(9,17) = tt261
hes(9,18) = tt262
hes(9,19) = tt224
hes(9,20) = tt245
hes(9,21) = tt263
hes(9,22) = tt264
hes(9,23) = tt265
hes(9,24) = tt266
hes(9,25) = tt267
hes(9,26) = tt268
hes(9,27) = tt269
hes(9,28) = tt270
hes(9,29) = tt271
hes(9,30) = tt272
hes(9,31) = tt273
hes(10,1) = tt68
hes(10,2) = tt104
hes(10,3) = tt135
hes(10,4) = tt161
hes(10,5) = tt162
hes(10,6) = tt163
hes(10,7) = tt213
hes(10,8) = tt236
hes(10,9) = tt256
hes(10,10) = (-2*tt31*tt56*tt210*W(28,1)*tt48)+4*tt153*tt156*tt21&
&0*W(28,1)*tt48+2*tt153*tt154*tt210*tt155
hes(10,11) = tt274
hes(10,12) = tt275
hes(10,13) = tt219
hes(10,14) = tt241
hes(10,15) = tt260
hes(10,16) = tt276
hes(10,17) = tt277
hes(10,18) = tt278
hes(10,19) = tt225
hes(10,20) = tt246
hes(10,21) = tt264
hes(10,22) = tt279
hes(10,23) = tt280
hes(10,24) = tt281
hes(10,25) = tt282
hes(10,26) = tt283
hes(10,27) = tt284
hes(10,28) = tt285
hes(10,29) = tt286
hes(10,30) = tt287
hes(10,31) = tt288
hes(11,1) = tt69
hes(11,2) = tt105
hes(11,3) = tt136
hes(11,4) = tt162
hes(11,5) = tt183
hes(11,6) = tt184
hes(11,7) = tt214
hes(11,8) = tt237
hes(11,9) = tt257
hes(11,10) = tt274
hes(11,11) = (-2*tt38*tt59*tt210*W(29,1)*tt48)+4*tt177*tt180*tt21&
&0*W(29,1)*tt48+2*tt177*tt178*tt210*tt179
hes(11,12) = tt289
hes(11,13) = tt220
hes(11,14) = tt242
hes(11,15) = tt261
hes(11,16) = tt277
hes(11,17) = tt290
hes(11,18) = tt291
hes(11,19) = tt226
hes(11,20) = tt247
hes(11,21) = tt265
hes(11,22) = tt280
hes(11,23) = tt292
hes(11,24) = tt293
hes(11,25) = tt294
hes(11,26) = tt295
hes(11,27) = tt296
hes(11,28) = tt297
hes(11,29) = tt298
hes(11,30) = tt299
hes(11,31) = tt300
hes(12,1) = tt70
hes(12,2) = tt106
hes(12,3) = tt137
hes(12,4) = tt163
hes(12,5) = tt184
hes(12,6) = tt200
hes(12,7) = tt215
hes(12,8) = tt238
hes(12,9) = tt258
hes(12,10) = tt275
hes(12,11) = tt289
hes(12,12) = (-2*tt45*tt62*tt210*W(30,1)*tt48)+4*tt196*tt199*tt21&
&0*W(30,1)*tt48+2*tt196*tt197*tt210*tt198
hes(12,13) = tt221
hes(12,14) = tt243
hes(12,15) = tt262
hes(12,16) = tt278
hes(12,17) = tt291
hes(12,18) = tt301
hes(12,19) = tt227
hes(12,20) = tt248
hes(12,21) = tt266
hes(12,22) = tt281
hes(12,23) = tt293
hes(12,24) = tt302
hes(12,25) = tt303
hes(12,26) = tt304
hes(12,27) = tt305
hes(12,28) = tt306
hes(12,29) = tt307
hes(12,30) = tt308
hes(12,31) = tt309
hes(13,1) = tt71
hes(13,2) = tt72
hes(13,3) = tt73
hes(13,4) = tt74
hes(13,5) = tt75
hes(13,6) = tt76
hes(13,7) = tt216
hes(13,8) = tt217
hes(13,9) = tt218
hes(13,10) = tt219
hes(13,11) = tt220
hes(13,12) = tt221
hes(13,13) = (-2*tt6*tt49*tt310*W(25,1)*tt48)+4*tt1*tt11*tt310*W(&
&25,1)*tt48+2*tt1*tt8*tt310*tt10
hes(13,14) = tt311
hes(13,15) = tt312
hes(13,16) = tt313
hes(13,17) = tt314
hes(13,18) = tt315
hes(13,19) = tt316
hes(13,20) = tt317
hes(13,21) = tt318
hes(13,22) = tt319
hes(13,23) = tt320
hes(13,24) = tt321
hes(13,25) = tt322
hes(13,26) = tt323
hes(13,27) = tt324
hes(13,28) = tt325
hes(13,29) = tt326
hes(13,30) = tt327
hes(13,31) = tt328
hes(14,1) = tt72
hes(14,2) = tt107
hes(14,3) = tt108
hes(14,4) = tt109
hes(14,5) = tt110
hes(14,6) = tt111
hes(14,7) = tt217
hes(14,8) = tt239
hes(14,9) = tt240
hes(14,10) = tt241
hes(14,11) = tt242
hes(14,12) = tt243
hes(14,13) = tt311
hes(14,14) = (-2*tt17*tt50*tt310*W(26,1)*tt48)+4*tt90*tt93*tt310*&
&W(26,1)*tt48+2*tt90*tt91*tt310*tt92
hes(14,15) = tt329
hes(14,16) = tt330
hes(14,17) = tt331
hes(14,18) = tt332
hes(14,19) = tt317
hes(14,20) = tt333
hes(14,21) = tt334
hes(14,22) = tt335
hes(14,23) = tt336
hes(14,24) = tt337
hes(14,25) = tt338
hes(14,26) = tt339
hes(14,27) = tt340
hes(14,28) = tt341
hes(14,29) = tt342
hes(14,30) = tt343
hes(14,31) = tt344
hes(15,1) = tt73
hes(15,2) = tt108
hes(15,3) = tt138
hes(15,4) = tt139
hes(15,5) = tt140
hes(15,6) = tt141
hes(15,7) = tt218
hes(15,8) = tt240
hes(15,9) = tt259
hes(15,10) = tt260
hes(15,11) = tt261
hes(15,12) = tt262
hes(15,13) = tt312
hes(15,14) = tt329
hes(15,15) = (-2*tt24*tt53*tt310*W(27,1)*tt48)+4*tt124*tt127*tt31&
&0*W(27,1)*tt48+2*tt124*tt125*tt310*tt126
hes(15,16) = tt345
hes(15,17) = tt346
hes(15,18) = tt347
hes(15,19) = tt318
hes(15,20) = tt334
hes(15,21) = tt348
hes(15,22) = tt349
hes(15,23) = tt350
hes(15,24) = tt351
hes(15,25) = tt352
hes(15,26) = tt353
hes(15,27) = tt354
hes(15,28) = tt355
hes(15,29) = tt356
hes(15,30) = tt357
hes(15,31) = tt358
hes(16,1) = tt74
hes(16,2) = tt109
hes(16,3) = tt139
hes(16,4) = tt164
hes(16,5) = tt165
hes(16,6) = tt166
hes(16,7) = tt219
hes(16,8) = tt241
hes(16,9) = tt260
hes(16,10) = tt276
hes(16,11) = tt277
hes(16,12) = tt278
hes(16,13) = tt313
hes(16,14) = tt330
hes(16,15) = tt345
hes(16,16) = (-2*tt31*tt56*tt310*W(28,1)*tt48)+4*tt153*tt156*tt31&
&0*W(28,1)*tt48+2*tt153*tt154*tt310*tt155
hes(16,17) = tt359
hes(16,18) = tt360
hes(16,19) = tt319
hes(16,20) = tt335
hes(16,21) = tt349
hes(16,22) = tt361
hes(16,23) = tt362
hes(16,24) = tt363
hes(16,25) = tt364
hes(16,26) = tt365
hes(16,27) = tt366
hes(16,28) = tt367
hes(16,29) = tt368
hes(16,30) = tt369
hes(16,31) = tt370
hes(17,1) = tt75
hes(17,2) = tt110
hes(17,3) = tt140
hes(17,4) = tt165
hes(17,5) = tt185
hes(17,6) = tt186
hes(17,7) = tt220
hes(17,8) = tt242
hes(17,9) = tt261
hes(17,10) = tt277
hes(17,11) = tt290
hes(17,12) = tt291
hes(17,13) = tt314
hes(17,14) = tt331
hes(17,15) = tt346
hes(17,16) = tt359
hes(17,17) = (-2*tt38*tt59*tt310*W(29,1)*tt48)+4*tt177*tt180*tt31&
&0*W(29,1)*tt48+2*tt177*tt178*tt310*tt179
hes(17,18) = tt371
hes(17,19) = tt320
hes(17,20) = tt336
hes(17,21) = tt350
hes(17,22) = tt362
hes(17,23) = tt372
hes(17,24) = tt373
hes(17,25) = tt374
hes(17,26) = tt375
hes(17,27) = tt376
hes(17,28) = tt377
hes(17,29) = tt378
hes(17,30) = tt379
hes(17,31) = tt380
hes(18,1) = tt76
hes(18,2) = tt111
hes(18,3) = tt141
hes(18,4) = tt166
hes(18,5) = tt186
hes(18,6) = tt201
hes(18,7) = tt221
hes(18,8) = tt243
hes(18,9) = tt262
hes(18,10) = tt278
hes(18,11) = tt291
hes(18,12) = tt301
hes(18,13) = tt315
hes(18,14) = tt332
hes(18,15) = tt347
hes(18,16) = tt360
hes(18,17) = tt371
hes(18,18) = (-2*tt45*tt62*tt310*W(30,1)*tt48)+4*tt196*tt199*tt31&
&0*W(30,1)*tt48+2*tt196*tt197*tt310*tt198
hes(18,19) = tt321
hes(18,20) = tt337
hes(18,21) = tt351
hes(18,22) = tt363
hes(18,23) = tt373
hes(18,24) = tt381
hes(18,25) = tt382
hes(18,26) = tt383
hes(18,27) = tt384
hes(18,28) = tt385
hes(18,29) = tt386
hes(18,30) = tt387
hes(18,31) = tt388
hes(19,1) = tt77
hes(19,2) = tt78
hes(19,3) = tt79
hes(19,4) = tt80
hes(19,5) = tt81
hes(19,6) = tt82
hes(19,7) = tt222
hes(19,8) = tt223
hes(19,9) = tt224
hes(19,10) = tt225
hes(19,11) = tt226
hes(19,12) = tt227
hes(19,13) = tt316
hes(19,14) = tt317
hes(19,15) = tt318
hes(19,16) = tt319
hes(19,17) = tt320
hes(19,18) = tt321
hes(19,19) = (-2*tt6*tt49*W(25,1)*tt48)+4*tt1*tt11*W(25,1)*tt48+2&
&*tt1*tt8*tt10
hes(19,20) = tt389
hes(19,21) = tt390
hes(19,22) = tt391
hes(19,23) = tt392
hes(19,24) = tt393
hes(19,25) = tt394
hes(19,26) = tt395
hes(19,27) = tt396
hes(19,28) = tt397
hes(19,29) = tt398
hes(19,30) = tt399
hes(19,31) = tt400
hes(20,1) = tt78
hes(20,2) = tt112
hes(20,3) = tt113
hes(20,4) = tt114
hes(20,5) = tt115
hes(20,6) = tt116
hes(20,7) = tt223
hes(20,8) = tt244
hes(20,9) = tt245
hes(20,10) = tt246
hes(20,11) = tt247
hes(20,12) = tt248
hes(20,13) = tt317
hes(20,14) = tt333
hes(20,15) = tt334
hes(20,16) = tt335
hes(20,17) = tt336
hes(20,18) = tt337
hes(20,19) = tt389
hes(20,20) = (-2*tt17*tt50*W(26,1)*tt48)+4*tt90*tt93*W(26,1)*tt48&
&+2*tt90*tt91*tt92
hes(20,21) = tt401
hes(20,22) = tt402
hes(20,23) = tt403
hes(20,24) = tt404
hes(20,25) = tt405
hes(20,26) = tt406
hes(20,27) = tt407
hes(20,28) = tt408
hes(20,29) = tt409
hes(20,30) = tt410
hes(20,31) = tt411
hes(21,1) = tt79
hes(21,2) = tt113
hes(21,3) = tt142
hes(21,4) = tt143
hes(21,5) = tt144
hes(21,6) = tt145
hes(21,7) = tt224
hes(21,8) = tt245
hes(21,9) = tt263
hes(21,10) = tt264
hes(21,11) = tt265
hes(21,12) = tt266
hes(21,13) = tt318
hes(21,14) = tt334
hes(21,15) = tt348
hes(21,16) = tt349
hes(21,17) = tt350
hes(21,18) = tt351
hes(21,19) = tt390
hes(21,20) = tt401
hes(21,21) = (-2*tt24*tt53*W(27,1)*tt48)+4*tt124*tt127*W(27,1)*tt&
&48+2*tt124*tt125*tt126
hes(21,22) = tt412
hes(21,23) = tt413
hes(21,24) = tt414
hes(21,25) = tt415
hes(21,26) = tt416
hes(21,27) = tt417
hes(21,28) = tt418
hes(21,29) = tt419
hes(21,30) = tt420
hes(21,31) = tt421
hes(22,1) = tt80
hes(22,2) = tt114
hes(22,3) = tt143
hes(22,4) = tt167
hes(22,5) = tt168
hes(22,6) = tt169
hes(22,7) = tt225
hes(22,8) = tt246
hes(22,9) = tt264
hes(22,10) = tt279
hes(22,11) = tt280
hes(22,12) = tt281
hes(22,13) = tt319
hes(22,14) = tt335
hes(22,15) = tt349
hes(22,16) = tt361
hes(22,17) = tt362
hes(22,18) = tt363
hes(22,19) = tt391
hes(22,20) = tt402
hes(22,21) = tt412
hes(22,22) = (-2*tt31*tt56*W(28,1)*tt48)+4*tt153*tt156*W(28,1)*tt&
&48+2*tt153*tt154*tt155
hes(22,23) = tt422
hes(22,24) = tt423
hes(22,25) = tt424
hes(22,26) = tt425
hes(22,27) = tt426
hes(22,28) = tt427
hes(22,29) = tt428
hes(22,30) = tt429
hes(22,31) = tt430
hes(23,1) = tt81
hes(23,2) = tt115
hes(23,3) = tt144
hes(23,4) = tt168
hes(23,5) = tt187
hes(23,6) = tt188
hes(23,7) = tt226
hes(23,8) = tt247
hes(23,9) = tt265
hes(23,10) = tt280
hes(23,11) = tt292
hes(23,12) = tt293
hes(23,13) = tt320
hes(23,14) = tt336
hes(23,15) = tt350
hes(23,16) = tt362
hes(23,17) = tt372
hes(23,18) = tt373
hes(23,19) = tt392
hes(23,20) = tt403
hes(23,21) = tt413
hes(23,22) = tt422
hes(23,23) = (-2*tt38*tt59*W(29,1)*tt48)+4*tt177*tt180*W(29,1)*tt&
&48+2*tt177*tt178*tt179
hes(23,24) = tt431
hes(23,25) = tt432
hes(23,26) = tt433
hes(23,27) = tt434
hes(23,28) = tt435
hes(23,29) = tt436
hes(23,30) = tt437
hes(23,31) = tt438
hes(24,1) = tt82
hes(24,2) = tt116
hes(24,3) = tt145
hes(24,4) = tt169
hes(24,5) = tt188
hes(24,6) = tt202
hes(24,7) = tt227
hes(24,8) = tt248
hes(24,9) = tt266
hes(24,10) = tt281
hes(24,11) = tt293
hes(24,12) = tt302
hes(24,13) = tt321
hes(24,14) = tt337
hes(24,15) = tt351
hes(24,16) = tt363
hes(24,17) = tt373
hes(24,18) = tt381
hes(24,19) = tt393
hes(24,20) = tt404
hes(24,21) = tt414
hes(24,22) = tt423
hes(24,23) = tt431
hes(24,24) = (-2*tt45*tt62*W(30,1)*tt48)+4*tt196*tt199*W(30,1)*tt&
&48+2*tt196*tt197*tt198
hes(24,25) = tt439
hes(24,26) = tt440
hes(24,27) = tt441
hes(24,28) = tt442
hes(24,29) = tt443
hes(24,30) = tt444
hes(24,31) = tt445
hes(25,1) = tt83
hes(25,2) = tt117
hes(25,3) = tt146
hes(25,4) = tt170
hes(25,5) = tt189
hes(25,6) = tt203
hes(25,7) = tt228
hes(25,8) = tt249
hes(25,9) = tt267
hes(25,10) = tt282
hes(25,11) = tt294
hes(25,12) = tt303
hes(25,13) = tt322
hes(25,14) = tt338
hes(25,15) = tt352
hes(25,16) = tt364
hes(25,17) = tt374
hes(25,18) = tt382
hes(25,19) = tt394
hes(25,20) = tt405
hes(25,21) = tt415
hes(25,22) = tt424
hes(25,23) = tt432
hes(25,24) = tt439
hes(25,25) = 2*tt49
hes(25,26) = tt446
hes(25,27) = tt447
hes(25,28) = tt448
hes(25,29) = tt449
hes(25,30) = tt450
hes(25,31) = tt451
hes(26,1) = tt84
hes(26,2) = tt118
hes(26,3) = tt147
hes(26,4) = tt171
hes(26,5) = tt190
hes(26,6) = tt204
hes(26,7) = tt229
hes(26,8) = tt250
hes(26,9) = tt268
hes(26,10) = tt283
hes(26,11) = tt295
hes(26,12) = tt304
hes(26,13) = tt323
hes(26,14) = tt339
hes(26,15) = tt353
hes(26,16) = tt365
hes(26,17) = tt375
hes(26,18) = tt383
hes(26,19) = tt395
hes(26,20) = tt406
hes(26,21) = tt416
hes(26,22) = tt425
hes(26,23) = tt433
hes(26,24) = tt440
hes(26,25) = tt446
hes(26,26) = 2*tt50
hes(26,27) = tt452
hes(26,28) = tt453
hes(26,29) = tt454
hes(26,30) = tt455
hes(26,31) = tt456
hes(27,1) = tt85
hes(27,2) = tt119
hes(27,3) = tt148
hes(27,4) = tt172
hes(27,5) = tt191
hes(27,6) = tt205
hes(27,7) = tt230
hes(27,8) = tt251
hes(27,9) = tt269
hes(27,10) = tt284
hes(27,11) = tt296
hes(27,12) = tt305
hes(27,13) = tt324
hes(27,14) = tt340
hes(27,15) = tt354
hes(27,16) = tt366
hes(27,17) = tt376
hes(27,18) = tt384
hes(27,19) = tt396
hes(27,20) = tt407
hes(27,21) = tt417
hes(27,22) = tt426
hes(27,23) = tt434
hes(27,24) = tt441
hes(27,25) = tt447
hes(27,26) = tt452
hes(27,27) = 2*tt53
hes(27,28) = tt457
hes(27,29) = tt458
hes(27,30) = tt459
hes(27,31) = tt460
hes(28,1) = tt86
hes(28,2) = tt120
hes(28,3) = tt149
hes(28,4) = tt173
hes(28,5) = tt192
hes(28,6) = tt206
hes(28,7) = tt231
hes(28,8) = tt252
hes(28,9) = tt270
hes(28,10) = tt285
hes(28,11) = tt297
hes(28,12) = tt306
hes(28,13) = tt325
hes(28,14) = tt341
hes(28,15) = tt355
hes(28,16) = tt367
hes(28,17) = tt377
hes(28,18) = tt385
hes(28,19) = tt397
hes(28,20) = tt408
hes(28,21) = tt418
hes(28,22) = tt427
hes(28,23) = tt435
hes(28,24) = tt442
hes(28,25) = tt448
hes(28,26) = tt453
hes(28,27) = tt457
hes(28,28) = 2*tt56
hes(28,29) = tt461
hes(28,30) = tt462
hes(28,31) = tt463
hes(29,1) = tt87
hes(29,2) = tt121
hes(29,3) = tt150
hes(29,4) = tt174
hes(29,5) = tt193
hes(29,6) = tt207
hes(29,7) = tt232
hes(29,8) = tt253
hes(29,9) = tt271
hes(29,10) = tt286
hes(29,11) = tt298
hes(29,12) = tt307
hes(29,13) = tt326
hes(29,14) = tt342
hes(29,15) = tt356
hes(29,16) = tt368
hes(29,17) = tt378
hes(29,18) = tt386
hes(29,19) = tt398
hes(29,20) = tt409
hes(29,21) = tt419
hes(29,22) = tt428
hes(29,23) = tt436
hes(29,24) = tt443
hes(29,25) = tt449
hes(29,26) = tt454
hes(29,27) = tt458
hes(29,28) = tt461
hes(29,29) = 2*tt59
hes(29,30) = tt464
hes(29,31) = tt465
hes(30,1) = tt88
hes(30,2) = tt122
hes(30,3) = tt151
hes(30,4) = tt175
hes(30,5) = tt194
hes(30,6) = tt208
hes(30,7) = tt233
hes(30,8) = tt254
hes(30,9) = tt272
hes(30,10) = tt287
hes(30,11) = tt299
hes(30,12) = tt308
hes(30,13) = tt327
hes(30,14) = tt343
hes(30,15) = tt357
hes(30,16) = tt369
hes(30,17) = tt379
hes(30,18) = tt387
hes(30,19) = tt399
hes(30,20) = tt410
hes(30,21) = tt420
hes(30,22) = tt429
hes(30,23) = tt437
hes(30,24) = tt444
hes(30,25) = tt450
hes(30,26) = tt455
hes(30,27) = tt459
hes(30,28) = tt462
hes(30,29) = tt464
hes(30,30) = 2*tt62
hes(30,31) = tt466
hes(31,1) = tt89
hes(31,2) = tt123
hes(31,3) = tt152
hes(31,4) = tt176
hes(31,5) = tt195
hes(31,6) = tt209
hes(31,7) = tt234
hes(31,8) = tt255
hes(31,9) = tt273
hes(31,10) = tt288
hes(31,11) = tt300
hes(31,12) = tt309
hes(31,13) = tt328
hes(31,14) = tt344
hes(31,15) = tt358
hes(31,16) = tt370
hes(31,17) = tt380
hes(31,18) = tt388
hes(31,19) = tt400
hes(31,20) = tt411
hes(31,21) = tt421
hes(31,22) = tt430
hes(31,23) = tt438
hes(31,24) = tt445
hes(31,25) = tt451
hes(31,26) = tt456
hes(31,27) = tt460
hes(31,28) = tt463
hes(31,29) = tt465
hes(31,30) = tt466
hes(31,31) = 2
END 
SUBROUTINE nnet_tri_elas(val, NODS, Dm, W) 
IMPLICIT NONE 
REAL(KIND=8) val(1, 1) 
REAL(KIND=8) NODS(2, 3) 
REAL(KIND=8) Dm(2, 2) 
REAL(KIND=8) W(31, 1) 
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
tt1 = -NODS(1,1)
tt2 = NODS(1,2)+tt1
tt3 = NODS(1,3)+tt1
tt4 = tt3*Dm(2,1)+Dm(1,1)*tt2
tt5 = -NODS(2,1)
tt6 = NODS(2,2)+tt5
tt7 = NODS(2,3)+tt5
tt8 = Dm(2,1)*tt7+Dm(1,1)*tt6
tt9 = tt8**2+tt4**2
tt10 = tt3*Dm(2,2)+Dm(1,2)*tt2
tt11 = Dm(2,2)*tt7+Dm(1,2)*tt6
tt12 = tt11**2+tt10**2
tt13 = tt8*tt11+tt4*tt10
val(1,1) = W(31,1)+W(30,1)/(exp((-W(24,1))-tt13*W(18,1)-tt12*W(12&
&,1)-tt9*W(6,1))+1)+W(29,1)/(exp((-W(23,1))-tt13*W(17,1)-tt12*W(11&
&,1)-tt9*W(5,1))+1)+W(28,1)/(exp((-W(22,1))-tt13*W(16,1)-tt12*W(10&
&,1)-tt9*W(4,1))+1)+W(27,1)/(exp((-W(21,1))-tt13*W(15,1)-tt12*W(9,&
&1)-tt9*W(3,1))+1)+W(26,1)/(exp((-W(20,1))-tt13*W(14,1)-tt12*W(8,1&
&)-W(2,1)*tt9)+1)+W(25,1)/(exp((-W(19,1))-tt13*W(13,1)-tt12*W(7,1)&
&-W(1,1)*tt9)+1)
END 
SUBROUTINE nnet_tri_elas_jac(jac, NODS, Dm, W) 
IMPLICIT NONE 
REAL(KIND=8) jac(1, 6) 
REAL(KIND=8) NODS(2, 3) 
REAL(KIND=8) Dm(2, 2) 
REAL(KIND=8) W(31, 1) 
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
tt1 = -NODS(1,1)
tt2 = NODS(1,2)+tt1
tt3 = NODS(1,3)+tt1
tt4 = tt3*Dm(2,1)+Dm(1,1)*tt2
tt5 = -NODS(2,1)
tt6 = NODS(2,2)+tt5
tt7 = NODS(2,3)+tt5
tt8 = Dm(2,1)*tt7+Dm(1,1)*tt6
tt9 = tt8**2+tt4**2
tt10 = tt3*Dm(2,2)+Dm(1,2)*tt2
tt11 = Dm(2,2)*tt7+Dm(1,2)*tt6
tt12 = tt11**2+tt10**2
tt13 = tt8*tt11+tt4*tt10
tt14 = exp((-W(19,1))-tt13*W(13,1)-tt12*W(7,1)-W(1,1)*tt9)
tt15 = 1/(tt14+1)**2
tt16 = (-Dm(2,1))-Dm(1,1)
tt17 = (-Dm(2,2))-Dm(1,2)
tt18 = tt16*tt10+tt4*tt17
tt19 = exp((-W(20,1))-tt13*W(14,1)-tt12*W(8,1)-W(2,1)*tt9)
tt20 = 1/(tt19+1)**2
tt21 = exp((-W(21,1))-tt13*W(15,1)-tt12*W(9,1)-tt9*W(3,1))
tt22 = 1/(tt21+1)**2
tt23 = exp((-W(22,1))-tt13*W(16,1)-tt12*W(10,1)-tt9*W(4,1))
tt24 = 1/(tt23+1)**2
tt25 = exp((-W(23,1))-tt13*W(17,1)-tt12*W(11,1)-tt9*W(5,1))
tt26 = 1/(tt25+1)**2
tt27 = exp((-W(24,1))-tt13*W(18,1)-tt12*W(12,1)-tt9*W(6,1))
tt28 = 1/(tt27+1)**2
tt29 = tt16*tt11+tt17*tt8
tt30 = Dm(1,1)*tt10+Dm(1,2)*tt4
tt31 = Dm(1,1)*tt11+Dm(1,2)*tt8
tt32 = Dm(2,1)*tt10+tt4*Dm(2,2)
tt33 = Dm(2,1)*tt11+Dm(2,2)*tt8
jac(1,1) = (-tt27*tt28*((-tt18*W(18,1))-2*tt17*tt10*W(12,1)-2*tt1&
&6*tt4*W(6,1))*W(30,1))-tt25*tt26*((-tt18*W(17,1))-2*tt17*tt10*W(1&
&1,1)-2*tt16*tt4*W(5,1))*W(29,1)-tt23*tt24*((-tt18*W(16,1))-2*tt17&
&*tt10*W(10,1)-2*tt16*tt4*W(4,1))*W(28,1)-tt21*tt22*((-tt18*W(15,1&
&))-2*tt17*tt10*W(9,1)-2*tt16*tt4*W(3,1))*W(27,1)-tt19*tt20*((-tt1&
&8*W(14,1))-2*tt17*tt10*W(8,1)-2*tt16*tt4*W(2,1))*W(26,1)-tt14*tt1&
&5*((-tt18*W(13,1))-2*tt17*tt10*W(7,1)-2*W(1,1)*tt16*tt4)*W(25,1)
jac(1,2) = (-tt27*tt28*((-tt29*W(18,1))-2*tt17*tt11*W(12,1)-2*tt1&
&6*tt8*W(6,1))*W(30,1))-tt25*tt26*((-tt29*W(17,1))-2*tt17*tt11*W(1&
&1,1)-2*tt16*tt8*W(5,1))*W(29,1)-tt23*tt24*((-tt29*W(16,1))-2*tt17&
&*tt11*W(10,1)-2*tt16*tt8*W(4,1))*W(28,1)-tt21*tt22*((-tt29*W(15,1&
&))-2*tt17*tt11*W(9,1)-2*tt16*tt8*W(3,1))*W(27,1)-tt19*tt20*((-tt2&
&9*W(14,1))-2*tt17*tt11*W(8,1)-2*tt16*W(2,1)*tt8)*W(26,1)-tt14*tt1&
&5*((-tt29*W(13,1))-2*tt17*tt11*W(7,1)-2*W(1,1)*tt16*tt8)*W(25,1)
jac(1,3) = (-tt27*tt28*((-tt30*W(18,1))-2*Dm(1,2)*tt10*W(12,1)-2*&
&Dm(1,1)*tt4*W(6,1))*W(30,1))-tt25*tt26*((-tt30*W(17,1))-2*Dm(1,2)&
&*tt10*W(11,1)-2*Dm(1,1)*tt4*W(5,1))*W(29,1)-tt23*tt24*((-tt30*W(1&
&6,1))-2*Dm(1,2)*tt10*W(10,1)-2*Dm(1,1)*tt4*W(4,1))*W(28,1)-tt21*t&
&t22*((-tt30*W(15,1))-2*Dm(1,2)*tt10*W(9,1)-2*Dm(1,1)*tt4*W(3,1))*&
&W(27,1)-tt19*tt20*((-tt30*W(14,1))-2*Dm(1,2)*tt10*W(8,1)-2*Dm(1,1&
&)*tt4*W(2,1))*W(26,1)-tt14*tt15*((-tt30*W(13,1))-2*Dm(1,2)*tt10*W&
&(7,1)-2*Dm(1,1)*W(1,1)*tt4)*W(25,1)
jac(1,4) = (-tt27*tt28*((-tt31*W(18,1))-2*Dm(1,2)*tt11*W(12,1)-2*&
&Dm(1,1)*tt8*W(6,1))*W(30,1))-tt25*tt26*((-tt31*W(17,1))-2*Dm(1,2)&
&*tt11*W(11,1)-2*Dm(1,1)*tt8*W(5,1))*W(29,1)-tt23*tt24*((-tt31*W(1&
&6,1))-2*Dm(1,2)*tt11*W(10,1)-2*Dm(1,1)*tt8*W(4,1))*W(28,1)-tt21*t&
&t22*((-tt31*W(15,1))-2*Dm(1,2)*tt11*W(9,1)-2*Dm(1,1)*tt8*W(3,1))*&
&W(27,1)-tt19*tt20*((-tt31*W(14,1))-2*Dm(1,2)*tt11*W(8,1)-2*Dm(1,1&
&)*W(2,1)*tt8)*W(26,1)-tt14*tt15*((-tt31*W(13,1))-2*Dm(1,2)*tt11*W&
&(7,1)-2*Dm(1,1)*W(1,1)*tt8)*W(25,1)
jac(1,5) = (-tt27*tt28*((-tt32*W(18,1))-2*Dm(2,2)*tt10*W(12,1)-2*&
&Dm(2,1)*tt4*W(6,1))*W(30,1))-tt25*tt26*((-tt32*W(17,1))-2*Dm(2,2)&
&*tt10*W(11,1)-2*Dm(2,1)*tt4*W(5,1))*W(29,1)-tt23*tt24*((-tt32*W(1&
&6,1))-2*Dm(2,2)*tt10*W(10,1)-2*Dm(2,1)*tt4*W(4,1))*W(28,1)-tt21*t&
&t22*((-tt32*W(15,1))-2*Dm(2,2)*tt10*W(9,1)-2*Dm(2,1)*tt4*W(3,1))*&
&W(27,1)-tt19*tt20*((-tt32*W(14,1))-2*Dm(2,2)*tt10*W(8,1)-2*Dm(2,1&
&)*tt4*W(2,1))*W(26,1)-tt14*tt15*((-tt32*W(13,1))-2*Dm(2,2)*tt10*W&
&(7,1)-2*W(1,1)*Dm(2,1)*tt4)*W(25,1)
jac(1,6) = (-tt27*tt28*((-tt33*W(18,1))-2*Dm(2,2)*tt11*W(12,1)-2*&
&Dm(2,1)*tt8*W(6,1))*W(30,1))-tt25*tt26*((-tt33*W(17,1))-2*Dm(2,2)&
&*tt11*W(11,1)-2*Dm(2,1)*tt8*W(5,1))*W(29,1)-tt23*tt24*((-tt33*W(1&
&6,1))-2*Dm(2,2)*tt11*W(10,1)-2*Dm(2,1)*tt8*W(4,1))*W(28,1)-tt21*t&
&t22*((-tt33*W(15,1))-2*Dm(2,2)*tt11*W(9,1)-2*Dm(2,1)*tt8*W(3,1))*&
&W(27,1)-tt19*tt20*((-tt33*W(14,1))-2*Dm(2,2)*tt11*W(8,1)-2*Dm(2,1&
&)*W(2,1)*tt8)*W(26,1)-tt14*tt15*((-tt33*W(13,1))-2*Dm(2,2)*tt11*W&
&(7,1)-2*W(1,1)*Dm(2,1)*tt8)*W(25,1)
END 
SUBROUTINE nnet_tri_elas_hes(hes, NODS, Dm, W) 
IMPLICIT NONE 
REAL(KIND=8) hes(6, 6) 
REAL(KIND=8) NODS(2, 3) 
REAL(KIND=8) Dm(2, 2) 
REAL(KIND=8) W(31, 1) 
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
tt1 = -NODS(1,1)
tt2 = NODS(1,2)+tt1
tt3 = NODS(1,3)+tt1
tt4 = tt3*Dm(2,1)+Dm(1,1)*tt2
tt5 = -NODS(2,1)
tt6 = NODS(2,2)+tt5
tt7 = NODS(2,3)+tt5
tt8 = Dm(2,1)*tt7+Dm(1,1)*tt6
tt9 = tt8**2+tt4**2
tt10 = tt3*Dm(2,2)+Dm(1,2)*tt2
tt11 = Dm(2,2)*tt7+Dm(1,2)*tt6
tt12 = tt11**2+tt10**2
tt13 = tt8*tt11+tt4*tt10
tt14 = exp((-W(19,1))-tt13*W(13,1)-tt12*W(7,1)-W(1,1)*tt9)
tt15 = tt14+1
tt16 = 1/tt15**2
tt17 = (-Dm(2,1))-Dm(1,1)
tt18 = tt17**2
tt19 = (-Dm(2,2))-Dm(1,2)
tt20 = tt19**2
tt21 = -tt14*tt16*((-2*tt17*tt19*W(13,1))-2*tt20*W(7,1)-2*W(1,1)*&
&tt18)*W(25,1)
tt22 = exp((-2*W(19,1))-2*tt13*W(13,1)-2*tt12*W(7,1)-2*W(1,1)*tt9&
&)
tt23 = 1/tt15**3
tt24 = tt17*tt10+tt4*tt19
tt25 = (-tt24*W(13,1))-2*tt19*tt10*W(7,1)-2*W(1,1)*tt17*tt4
tt26 = tt25**2
tt27 = exp((-W(20,1))-tt13*W(14,1)-tt12*W(8,1)-W(2,1)*tt9)
tt28 = tt27+1
tt29 = 1/tt28**2
tt30 = -tt27*tt29*((-2*tt17*tt19*W(14,1))-2*tt20*W(8,1)-2*tt18*W(&
&2,1))*W(26,1)
tt31 = exp((-2*W(20,1))-2*tt13*W(14,1)-2*tt12*W(8,1)-2*W(2,1)*tt9&
&)
tt32 = 1/tt28**3
tt33 = (-tt24*W(14,1))-2*tt19*tt10*W(8,1)-2*tt17*tt4*W(2,1)
tt34 = tt33**2
tt35 = exp((-W(21,1))-tt13*W(15,1)-tt12*W(9,1)-tt9*W(3,1))
tt36 = tt35+1
tt37 = 1/tt36**2
tt38 = -tt35*tt37*((-2*tt17*tt19*W(15,1))-2*tt20*W(9,1)-2*tt18*W(&
&3,1))*W(27,1)
tt39 = exp((-2*W(21,1))-2*tt13*W(15,1)-2*tt12*W(9,1)-2*tt9*W(3,1)&
&)
tt40 = 1/tt36**3
tt41 = (-tt24*W(15,1))-2*tt19*tt10*W(9,1)-2*tt17*tt4*W(3,1)
tt42 = tt41**2
tt43 = exp((-W(22,1))-tt13*W(16,1)-tt12*W(10,1)-tt9*W(4,1))
tt44 = tt43+1
tt45 = 1/tt44**2
tt46 = -tt43*tt45*((-2*tt17*tt19*W(16,1))-2*tt20*W(10,1)-2*tt18*W&
&(4,1))*W(28,1)
tt47 = exp((-2*W(22,1))-2*tt13*W(16,1)-2*tt12*W(10,1)-2*tt9*W(4,1&
&))
tt48 = 1/tt44**3
tt49 = (-tt24*W(16,1))-2*tt19*tt10*W(10,1)-2*tt17*tt4*W(4,1)
tt50 = tt49**2
tt51 = exp((-W(23,1))-tt13*W(17,1)-tt12*W(11,1)-tt9*W(5,1))
tt52 = tt51+1
tt53 = 1/tt52**2
tt54 = -tt51*tt53*((-2*tt17*tt19*W(17,1))-2*tt20*W(11,1)-2*tt18*W&
&(5,1))*W(29,1)
tt55 = exp((-2*W(23,1))-2*tt13*W(17,1)-2*tt12*W(11,1)-2*tt9*W(5,1&
&))
tt56 = 1/tt52**3
tt57 = (-tt24*W(17,1))-2*tt19*tt10*W(11,1)-2*tt17*tt4*W(5,1)
tt58 = tt57**2
tt59 = exp((-W(24,1))-tt13*W(18,1)-tt12*W(12,1)-tt9*W(6,1))
tt60 = tt59+1
tt61 = 1/tt60**2
tt62 = -tt59*tt61*((-2*tt17*tt19*W(18,1))-2*tt20*W(12,1)-2*tt18*W&
&(6,1))*W(30,1)
tt63 = exp((-2*W(24,1))-2*tt13*W(18,1)-2*tt12*W(12,1)-2*tt9*W(6,1&
&))
tt64 = 1/tt60**3
tt65 = (-tt24*W(18,1))-2*tt19*tt10*W(12,1)-2*tt17*tt4*W(6,1)
tt66 = tt65**2
tt67 = tt17*tt11+tt19*tt8
tt68 = (-tt67*W(13,1))-2*tt19*tt11*W(7,1)-2*W(1,1)*tt17*tt8
tt69 = (-tt67*W(14,1))-2*tt19*tt11*W(8,1)-2*tt17*W(2,1)*tt8
tt70 = (-tt67*W(15,1))-2*tt19*tt11*W(9,1)-2*tt17*tt8*W(3,1)
tt71 = (-tt67*W(16,1))-2*tt19*tt11*W(10,1)-2*tt17*tt8*W(4,1)
tt72 = (-tt67*W(17,1))-2*tt19*tt11*W(11,1)-2*tt17*tt8*W(5,1)
tt73 = (-tt67*W(18,1))-2*tt19*tt11*W(12,1)-2*tt17*tt8*W(6,1)
tt74 = (-tt59*tt61*tt65*tt73*W(30,1))+2*tt63*tt64*tt65*tt73*W(30,&
&1)-tt51*tt53*tt57*tt72*W(29,1)+2*tt55*tt56*tt57*tt72*W(29,1)-tt43&
&*tt45*tt49*tt71*W(28,1)+2*tt47*tt48*tt49*tt71*W(28,1)-tt35*tt37*t&
&t41*tt70*W(27,1)+2*tt39*tt40*tt41*tt70*W(27,1)-tt27*tt29*tt33*tt6&
&9*W(26,1)+2*tt31*tt32*tt33*tt69*W(26,1)-tt14*tt16*tt25*tt68*W(25,&
&1)+2*tt22*tt23*tt25*tt68*W(25,1)
tt75 = Dm(1,1)*tt19+Dm(1,2)*tt17
tt76 = -tt14*tt16*((-tt75*W(13,1))-2*Dm(1,2)*tt19*W(7,1)-2*Dm(1,1&
&)*W(1,1)*tt17)*W(25,1)
tt77 = Dm(1,1)*tt10+Dm(1,2)*tt4
tt78 = (-tt77*W(13,1))-2*Dm(1,2)*tt10*W(7,1)-2*Dm(1,1)*W(1,1)*tt4&
&
tt79 = -tt27*tt29*((-tt75*W(14,1))-2*Dm(1,2)*tt19*W(8,1)-2*Dm(1,1&
&)*tt17*W(2,1))*W(26,1)
tt80 = (-tt77*W(14,1))-2*Dm(1,2)*tt10*W(8,1)-2*Dm(1,1)*tt4*W(2,1)&
&
tt81 = -tt35*tt37*((-tt75*W(15,1))-2*Dm(1,2)*tt19*W(9,1)-2*Dm(1,1&
&)*tt17*W(3,1))*W(27,1)
tt82 = (-tt77*W(15,1))-2*Dm(1,2)*tt10*W(9,1)-2*Dm(1,1)*tt4*W(3,1)&
&
tt83 = -tt43*tt45*((-tt75*W(16,1))-2*Dm(1,2)*tt19*W(10,1)-2*Dm(1,&
&1)*tt17*W(4,1))*W(28,1)
tt84 = (-tt77*W(16,1))-2*Dm(1,2)*tt10*W(10,1)-2*Dm(1,1)*tt4*W(4,1&
&)
tt85 = -tt51*tt53*((-tt75*W(17,1))-2*Dm(1,2)*tt19*W(11,1)-2*Dm(1,&
&1)*tt17*W(5,1))*W(29,1)
tt86 = (-tt77*W(17,1))-2*Dm(1,2)*tt10*W(11,1)-2*Dm(1,1)*tt4*W(5,1&
&)
tt87 = -tt59*tt61*((-tt75*W(18,1))-2*Dm(1,2)*tt19*W(12,1)-2*Dm(1,&
&1)*tt17*W(6,1))*W(30,1)
tt88 = (-tt77*W(18,1))-2*Dm(1,2)*tt10*W(12,1)-2*Dm(1,1)*tt4*W(6,1&
&)
tt89 = (-tt59*tt61*tt88*tt65*W(30,1))+2*tt63*tt64*tt88*tt65*W(30,&
&1)+tt87-tt51*tt53*tt86*tt57*W(29,1)+2*tt55*tt56*tt86*tt57*W(29,1)&
&+tt85-tt43*tt45*tt84*tt49*W(28,1)+2*tt47*tt48*tt84*tt49*W(28,1)+t&
&t83-tt35*tt37*tt82*tt41*W(27,1)+2*tt39*tt40*tt82*tt41*W(27,1)+tt8&
&1-tt27*tt29*tt80*tt33*W(26,1)+2*tt31*tt32*tt80*tt33*W(26,1)+tt79-&
&tt14*tt16*tt78*tt25*W(25,1)+2*tt22*tt23*tt78*tt25*W(25,1)+tt76
tt90 = Dm(1,1)*tt11+Dm(1,2)*tt8
tt91 = (-tt90*W(13,1))-2*Dm(1,2)*tt11*W(7,1)-2*Dm(1,1)*W(1,1)*tt8&
&
tt92 = (-tt90*W(14,1))-2*Dm(1,2)*tt11*W(8,1)-2*Dm(1,1)*W(2,1)*tt8&
&
tt93 = (-tt90*W(15,1))-2*Dm(1,2)*tt11*W(9,1)-2*Dm(1,1)*tt8*W(3,1)&
&
tt94 = (-tt90*W(16,1))-2*Dm(1,2)*tt11*W(10,1)-2*Dm(1,1)*tt8*W(4,1&
&)
tt95 = (-tt90*W(17,1))-2*Dm(1,2)*tt11*W(11,1)-2*Dm(1,1)*tt8*W(5,1&
&)
tt96 = (-tt90*W(18,1))-2*Dm(1,2)*tt11*W(12,1)-2*Dm(1,1)*tt8*W(6,1&
&)
tt97 = (-tt59*tt61*tt65*tt96*W(30,1))+2*tt63*tt64*tt65*tt96*W(30,&
&1)-tt51*tt53*tt57*tt95*W(29,1)+2*tt55*tt56*tt57*tt95*W(29,1)-tt43&
&*tt45*tt49*tt94*W(28,1)+2*tt47*tt48*tt49*tt94*W(28,1)-tt35*tt37*t&
&t41*tt93*W(27,1)+2*tt39*tt40*tt41*tt93*W(27,1)-tt27*tt29*tt33*tt9&
&2*W(26,1)+2*tt31*tt32*tt33*tt92*W(26,1)-tt14*tt16*tt25*tt91*W(25,&
&1)+2*tt22*tt23*tt25*tt91*W(25,1)
tt98 = tt17*Dm(2,2)+Dm(2,1)*tt19
tt99 = -tt14*tt16*((-tt98*W(13,1))-2*tt19*Dm(2,2)*W(7,1)-2*W(1,1)&
&*tt17*Dm(2,1))*W(25,1)
tt100 = Dm(2,1)*tt10+tt4*Dm(2,2)
tt101 = (-tt100*W(13,1))-2*Dm(2,2)*tt10*W(7,1)-2*W(1,1)*Dm(2,1)*t&
&t4
tt102 = -tt27*tt29*((-tt98*W(14,1))-2*tt19*Dm(2,2)*W(8,1)-2*tt17*&
&Dm(2,1)*W(2,1))*W(26,1)
tt103 = (-tt100*W(14,1))-2*Dm(2,2)*tt10*W(8,1)-2*Dm(2,1)*tt4*W(2,&
&1)
tt104 = -tt35*tt37*((-tt98*W(15,1))-2*tt19*Dm(2,2)*W(9,1)-2*tt17*&
&Dm(2,1)*W(3,1))*W(27,1)
tt105 = (-tt100*W(15,1))-2*Dm(2,2)*tt10*W(9,1)-2*Dm(2,1)*tt4*W(3,&
&1)
tt106 = -tt43*tt45*((-tt98*W(16,1))-2*tt19*Dm(2,2)*W(10,1)-2*tt17&
&*Dm(2,1)*W(4,1))*W(28,1)
tt107 = (-tt100*W(16,1))-2*Dm(2,2)*tt10*W(10,1)-2*Dm(2,1)*tt4*W(4&
&,1)
tt108 = -tt51*tt53*((-tt98*W(17,1))-2*tt19*Dm(2,2)*W(11,1)-2*tt17&
&*Dm(2,1)*W(5,1))*W(29,1)
tt109 = (-tt100*W(17,1))-2*Dm(2,2)*tt10*W(11,1)-2*Dm(2,1)*tt4*W(5&
&,1)
tt110 = -tt59*tt61*((-tt98*W(18,1))-2*tt19*Dm(2,2)*W(12,1)-2*tt17&
&*Dm(2,1)*W(6,1))*W(30,1)
tt111 = (-tt100*W(18,1))-2*Dm(2,2)*tt10*W(12,1)-2*Dm(2,1)*tt4*W(6&
&,1)
tt112 = (-tt59*tt61*tt65*tt111*W(30,1))+2*tt63*tt64*tt65*tt111*W(&
&30,1)+tt110-tt51*tt53*tt57*tt109*W(29,1)+2*tt55*tt56*tt57*tt109*W&
&(29,1)+tt108-tt43*tt45*tt49*tt107*W(28,1)+2*tt47*tt48*tt49*tt107*&
&W(28,1)+tt106-tt35*tt37*tt41*tt105*W(27,1)+2*tt39*tt40*tt41*tt105&
&*W(27,1)+tt104-tt27*tt29*tt33*tt103*W(26,1)+2*tt31*tt32*tt33*tt10&
&3*W(26,1)+tt102-tt14*tt16*tt25*tt101*W(25,1)+2*tt22*tt23*tt25*tt1&
&01*W(25,1)+tt99
tt113 = Dm(2,1)*tt11+Dm(2,2)*tt8
tt114 = (-tt113*W(13,1))-2*Dm(2,2)*tt11*W(7,1)-2*W(1,1)*Dm(2,1)*t&
&t8
tt115 = (-tt113*W(14,1))-2*Dm(2,2)*tt11*W(8,1)-2*Dm(2,1)*W(2,1)*t&
&t8
tt116 = (-tt113*W(15,1))-2*Dm(2,2)*tt11*W(9,1)-2*Dm(2,1)*tt8*W(3,&
&1)
tt117 = (-tt113*W(16,1))-2*Dm(2,2)*tt11*W(10,1)-2*Dm(2,1)*tt8*W(4&
&,1)
tt118 = (-tt113*W(17,1))-2*Dm(2,2)*tt11*W(11,1)-2*Dm(2,1)*tt8*W(5&
&,1)
tt119 = (-tt113*W(18,1))-2*Dm(2,2)*tt11*W(12,1)-2*Dm(2,1)*tt8*W(6&
&,1)
tt120 = (-tt59*tt61*tt65*tt119*W(30,1))+2*tt63*tt64*tt65*tt119*W(&
&30,1)-tt51*tt53*tt57*tt118*W(29,1)+2*tt55*tt56*tt57*tt118*W(29,1)&
&-tt43*tt45*tt49*tt117*W(28,1)+2*tt47*tt48*tt49*tt117*W(28,1)-tt35&
&*tt37*tt41*tt116*W(27,1)+2*tt39*tt40*tt41*tt116*W(27,1)-tt27*tt29&
&*tt33*tt115*W(26,1)+2*tt31*tt32*tt33*tt115*W(26,1)-tt14*tt16*tt25&
&*tt114*W(25,1)+2*tt22*tt23*tt25*tt114*W(25,1)
tt121 = tt68**2
tt122 = tt69**2
tt123 = tt70**2
tt124 = tt71**2
tt125 = tt72**2
tt126 = tt73**2
tt127 = (-tt59*tt61*tt88*tt73*W(30,1))+2*tt63*tt64*tt88*tt73*W(30&
&,1)-tt51*tt53*tt86*tt72*W(29,1)+2*tt55*tt56*tt86*tt72*W(29,1)-tt4&
&3*tt45*tt84*tt71*W(28,1)+2*tt47*tt48*tt84*tt71*W(28,1)-tt35*tt37*&
&tt82*tt70*W(27,1)+2*tt39*tt40*tt82*tt70*W(27,1)-tt27*tt29*tt80*tt&
&69*W(26,1)+2*tt31*tt32*tt80*tt69*W(26,1)-tt14*tt16*tt78*tt68*W(25&
&,1)+2*tt22*tt23*tt78*tt68*W(25,1)
tt128 = (-tt59*tt61*tt96*tt73*W(30,1))+2*tt63*tt64*tt96*tt73*W(30&
&,1)+tt87-tt51*tt53*tt95*tt72*W(29,1)+2*tt55*tt56*tt95*tt72*W(29,1&
&)+tt85-tt43*tt45*tt94*tt71*W(28,1)+2*tt47*tt48*tt94*tt71*W(28,1)+&
&tt83-tt35*tt37*tt93*tt70*W(27,1)+2*tt39*tt40*tt93*tt70*W(27,1)+tt&
&81-tt27*tt29*tt92*tt69*W(26,1)+2*tt31*tt32*tt92*tt69*W(26,1)+tt79&
&-tt14*tt16*tt91*tt68*W(25,1)+2*tt22*tt23*tt91*tt68*W(25,1)+tt76
tt129 = (-tt59*tt61*tt111*tt73*W(30,1))+2*tt63*tt64*tt111*tt73*W(&
&30,1)-tt51*tt53*tt109*tt72*W(29,1)+2*tt55*tt56*tt109*tt72*W(29,1)&
&-tt43*tt45*tt107*tt71*W(28,1)+2*tt47*tt48*tt107*tt71*W(28,1)-tt35&
&*tt37*tt105*tt70*W(27,1)+2*tt39*tt40*tt105*tt70*W(27,1)-tt27*tt29&
&*tt103*tt69*W(26,1)+2*tt31*tt32*tt103*tt69*W(26,1)-tt14*tt16*tt10&
&1*tt68*W(25,1)+2*tt22*tt23*tt101*tt68*W(25,1)
tt130 = (-tt59*tt61*tt73*tt119*W(30,1))+2*tt63*tt64*tt73*tt119*W(&
&30,1)+tt110-tt51*tt53*tt72*tt118*W(29,1)+2*tt55*tt56*tt72*tt118*W&
&(29,1)+tt108-tt43*tt45*tt71*tt117*W(28,1)+2*tt47*tt48*tt71*tt117*&
&W(28,1)+tt106-tt35*tt37*tt70*tt116*W(27,1)+2*tt39*tt40*tt70*tt116&
&*W(27,1)+tt104-tt27*tt29*tt69*tt115*W(26,1)+2*tt31*tt32*tt69*tt11&
&5*W(26,1)+tt102-tt14*tt16*tt68*tt114*W(25,1)+2*tt22*tt23*tt68*tt1&
&14*W(25,1)+tt99
tt131 = Dm(1,1)**2
tt132 = Dm(1,2)**2
tt133 = -tt14*tt16*((-2*Dm(1,1)*Dm(1,2)*W(13,1))-2*tt132*W(7,1)-2&
&*tt131*W(1,1))*W(25,1)
tt134 = tt78**2
tt135 = -tt27*tt29*((-2*Dm(1,1)*Dm(1,2)*W(14,1))-2*tt132*W(8,1)-2&
&*tt131*W(2,1))*W(26,1)
tt136 = tt80**2
tt137 = -tt35*tt37*((-2*Dm(1,1)*Dm(1,2)*W(15,1))-2*tt132*W(9,1)-2&
&*tt131*W(3,1))*W(27,1)
tt138 = tt82**2
tt139 = -tt43*tt45*((-2*Dm(1,1)*Dm(1,2)*W(16,1))-2*tt132*W(10,1)-&
&2*tt131*W(4,1))*W(28,1)
tt140 = tt84**2
tt141 = -tt51*tt53*((-2*Dm(1,1)*Dm(1,2)*W(17,1))-2*tt132*W(11,1)-&
&2*tt131*W(5,1))*W(29,1)
tt142 = tt86**2
tt143 = -tt59*tt61*((-2*Dm(1,1)*Dm(1,2)*W(18,1))-2*tt132*W(12,1)-&
&2*tt131*W(6,1))*W(30,1)
tt144 = tt88**2
tt145 = (-tt59*tt61*tt88*tt96*W(30,1))+2*tt63*tt64*tt88*tt96*W(30&
&,1)-tt51*tt53*tt86*tt95*W(29,1)+2*tt55*tt56*tt86*tt95*W(29,1)-tt4&
&3*tt45*tt84*tt94*W(28,1)+2*tt47*tt48*tt84*tt94*W(28,1)-tt35*tt37*&
&tt82*tt93*W(27,1)+2*tt39*tt40*tt82*tt93*W(27,1)-tt27*tt29*tt80*tt&
&92*W(26,1)+2*tt31*tt32*tt80*tt92*W(26,1)-tt14*tt16*tt78*tt91*W(25&
&,1)+2*tt22*tt23*tt78*tt91*W(25,1)
tt146 = Dm(1,1)*Dm(2,2)+Dm(1,2)*Dm(2,1)
tt147 = -tt14*tt16*((-tt146*W(13,1))-2*Dm(1,2)*Dm(2,2)*W(7,1)-2*D&
&m(1,1)*W(1,1)*Dm(2,1))*W(25,1)
tt148 = -tt27*tt29*((-tt146*W(14,1))-2*Dm(1,2)*Dm(2,2)*W(8,1)-2*D&
&m(1,1)*Dm(2,1)*W(2,1))*W(26,1)
tt149 = -tt35*tt37*((-tt146*W(15,1))-2*Dm(1,2)*Dm(2,2)*W(9,1)-2*D&
&m(1,1)*Dm(2,1)*W(3,1))*W(27,1)
tt150 = -tt43*tt45*((-tt146*W(16,1))-2*Dm(1,2)*Dm(2,2)*W(10,1)-2*&
&Dm(1,1)*Dm(2,1)*W(4,1))*W(28,1)
tt151 = -tt51*tt53*((-tt146*W(17,1))-2*Dm(1,2)*Dm(2,2)*W(11,1)-2*&
&Dm(1,1)*Dm(2,1)*W(5,1))*W(29,1)
tt152 = -tt59*tt61*((-tt146*W(18,1))-2*Dm(1,2)*Dm(2,2)*W(12,1)-2*&
&Dm(1,1)*Dm(2,1)*W(6,1))*W(30,1)
tt153 = (-tt59*tt61*tt88*tt111*W(30,1))+2*tt63*tt64*tt88*tt111*W(&
&30,1)+tt152-tt51*tt53*tt86*tt109*W(29,1)+2*tt55*tt56*tt86*tt109*W&
&(29,1)+tt151-tt43*tt45*tt84*tt107*W(28,1)+2*tt47*tt48*tt84*tt107*&
&W(28,1)+tt150-tt35*tt37*tt82*tt105*W(27,1)+2*tt39*tt40*tt82*tt105&
&*W(27,1)+tt149-tt27*tt29*tt80*tt103*W(26,1)+2*tt31*tt32*tt80*tt10&
&3*W(26,1)+tt148-tt14*tt16*tt78*tt101*W(25,1)+2*tt22*tt23*tt78*tt1&
&01*W(25,1)+tt147
tt154 = (-tt59*tt61*tt88*tt119*W(30,1))+2*tt63*tt64*tt88*tt119*W(&
&30,1)-tt51*tt53*tt86*tt118*W(29,1)+2*tt55*tt56*tt86*tt118*W(29,1)&
&-tt43*tt45*tt84*tt117*W(28,1)+2*tt47*tt48*tt84*tt117*W(28,1)-tt35&
&*tt37*tt82*tt116*W(27,1)+2*tt39*tt40*tt82*tt116*W(27,1)-tt27*tt29&
&*tt80*tt115*W(26,1)+2*tt31*tt32*tt80*tt115*W(26,1)-tt14*tt16*tt78&
&*tt114*W(25,1)+2*tt22*tt23*tt78*tt114*W(25,1)
tt155 = tt91**2
tt156 = tt92**2
tt157 = tt93**2
tt158 = tt94**2
tt159 = tt95**2
tt160 = tt96**2
tt161 = (-tt59*tt61*tt111*tt96*W(30,1))+2*tt63*tt64*tt111*tt96*W(&
&30,1)-tt51*tt53*tt109*tt95*W(29,1)+2*tt55*tt56*tt109*tt95*W(29,1)&
&-tt43*tt45*tt107*tt94*W(28,1)+2*tt47*tt48*tt107*tt94*W(28,1)-tt35&
&*tt37*tt105*tt93*W(27,1)+2*tt39*tt40*tt105*tt93*W(27,1)-tt27*tt29&
&*tt103*tt92*W(26,1)+2*tt31*tt32*tt103*tt92*W(26,1)-tt14*tt16*tt10&
&1*tt91*W(25,1)+2*tt22*tt23*tt101*tt91*W(25,1)
tt162 = (-tt59*tt61*tt96*tt119*W(30,1))+2*tt63*tt64*tt96*tt119*W(&
&30,1)+tt152-tt51*tt53*tt95*tt118*W(29,1)+2*tt55*tt56*tt95*tt118*W&
&(29,1)+tt151-tt43*tt45*tt94*tt117*W(28,1)+2*tt47*tt48*tt94*tt117*&
&W(28,1)+tt150-tt35*tt37*tt93*tt116*W(27,1)+2*tt39*tt40*tt93*tt116&
&*W(27,1)+tt149-tt27*tt29*tt92*tt115*W(26,1)+2*tt31*tt32*tt92*tt11&
&5*W(26,1)+tt148-tt14*tt16*tt91*tt114*W(25,1)+2*tt22*tt23*tt91*tt1&
&14*W(25,1)+tt147
tt163 = Dm(2,1)**2
tt164 = Dm(2,2)**2
tt165 = -tt14*tt16*((-2*Dm(2,1)*Dm(2,2)*W(13,1))-2*tt164*W(7,1)-2&
&*W(1,1)*tt163)*W(25,1)
tt166 = tt101**2
tt167 = -tt27*tt29*((-2*Dm(2,1)*Dm(2,2)*W(14,1))-2*tt164*W(8,1)-2&
&*tt163*W(2,1))*W(26,1)
tt168 = tt103**2
tt169 = -tt35*tt37*((-2*Dm(2,1)*Dm(2,2)*W(15,1))-2*tt164*W(9,1)-2&
&*tt163*W(3,1))*W(27,1)
tt170 = tt105**2
tt171 = -tt43*tt45*((-2*Dm(2,1)*Dm(2,2)*W(16,1))-2*tt164*W(10,1)-&
&2*tt163*W(4,1))*W(28,1)
tt172 = tt107**2
tt173 = -tt51*tt53*((-2*Dm(2,1)*Dm(2,2)*W(17,1))-2*tt164*W(11,1)-&
&2*tt163*W(5,1))*W(29,1)
tt174 = tt109**2
tt175 = -tt59*tt61*((-2*Dm(2,1)*Dm(2,2)*W(18,1))-2*tt164*W(12,1)-&
&2*tt163*W(6,1))*W(30,1)
tt176 = tt111**2
tt177 = (-tt59*tt61*tt111*tt119*W(30,1))+2*tt63*tt64*tt111*tt119*&
&W(30,1)-tt51*tt53*tt109*tt118*W(29,1)+2*tt55*tt56*tt109*tt118*W(2&
&9,1)-tt43*tt45*tt107*tt117*W(28,1)+2*tt47*tt48*tt107*tt117*W(28,1&
&)-tt35*tt37*tt105*tt116*W(27,1)+2*tt39*tt40*tt105*tt116*W(27,1)-t&
&t27*tt29*tt103*tt115*W(26,1)+2*tt31*tt32*tt103*tt115*W(26,1)-tt14&
&*tt16*tt101*tt114*W(25,1)+2*tt22*tt23*tt101*tt114*W(25,1)
tt178 = tt114**2
tt179 = tt115**2
tt180 = tt116**2
tt181 = tt117**2
tt182 = tt118**2
tt183 = tt119**2
hes(1,1) = (-tt59*tt61*tt66*W(30,1))+2*tt63*tt64*tt66*W(30,1)+tt6&
&2-tt51*tt53*tt58*W(29,1)+2*tt55*tt56*tt58*W(29,1)+tt54-tt43*tt45*&
&tt50*W(28,1)+2*tt47*tt48*tt50*W(28,1)+tt46-tt35*tt37*tt42*W(27,1)&
&+2*tt39*tt40*tt42*W(27,1)+tt38-tt27*tt29*tt34*W(26,1)+2*tt31*tt32&
&*tt34*W(26,1)+tt30-tt14*tt16*tt26*W(25,1)+2*tt22*tt23*tt26*W(25,1&
&)+tt21
hes(1,2) = tt74
hes(1,3) = tt89
hes(1,4) = tt97
hes(1,5) = tt112
hes(1,6) = tt120
hes(2,1) = tt74
hes(2,2) = (-tt59*tt61*tt126*W(30,1))+2*tt63*tt64*tt126*W(30,1)+t&
&t62-tt51*tt53*tt125*W(29,1)+2*tt55*tt56*tt125*W(29,1)+tt54-tt43*t&
&t45*tt124*W(28,1)+2*tt47*tt48*tt124*W(28,1)+tt46-tt35*tt37*tt123*&
&W(27,1)+2*tt39*tt40*tt123*W(27,1)+tt38-tt27*tt29*tt122*W(26,1)+2*&
&tt31*tt32*tt122*W(26,1)+tt30-tt14*tt16*tt121*W(25,1)+2*tt22*tt23*&
&tt121*W(25,1)+tt21
hes(2,3) = tt127
hes(2,4) = tt128
hes(2,5) = tt129
hes(2,6) = tt130
hes(3,1) = tt89
hes(3,2) = tt127
hes(3,3) = (-tt59*tt61*tt144*W(30,1))+2*tt63*tt64*tt144*W(30,1)+t&
&t143-tt51*tt53*tt142*W(29,1)+2*tt55*tt56*tt142*W(29,1)+tt141-tt43&
&*tt45*tt140*W(28,1)+2*tt47*tt48*tt140*W(28,1)+tt139-tt35*tt37*tt1&
&38*W(27,1)+2*tt39*tt40*tt138*W(27,1)+tt137-tt27*tt29*tt136*W(26,1&
&)+2*tt31*tt32*tt136*W(26,1)+tt135-tt14*tt16*tt134*W(25,1)+2*tt22*&
&tt23*tt134*W(25,1)+tt133
hes(3,4) = tt145
hes(3,5) = tt153
hes(3,6) = tt154
hes(4,1) = tt97
hes(4,2) = tt128
hes(4,3) = tt145
hes(4,4) = (-tt59*tt61*tt160*W(30,1))+2*tt63*tt64*tt160*W(30,1)+t&
&t143-tt51*tt53*tt159*W(29,1)+2*tt55*tt56*tt159*W(29,1)+tt141-tt43&
&*tt45*tt158*W(28,1)+2*tt47*tt48*tt158*W(28,1)+tt139-tt35*tt37*tt1&
&57*W(27,1)+2*tt39*tt40*tt157*W(27,1)+tt137-tt27*tt29*tt156*W(26,1&
&)+2*tt31*tt32*tt156*W(26,1)+tt135-tt14*tt16*tt155*W(25,1)+2*tt22*&
&tt23*tt155*W(25,1)+tt133
hes(4,5) = tt161
hes(4,6) = tt162
hes(5,1) = tt112
hes(5,2) = tt129
hes(5,3) = tt153
hes(5,4) = tt161
hes(5,5) = (-tt59*tt61*tt176*W(30,1))+2*tt63*tt64*tt176*W(30,1)+t&
&t175-tt51*tt53*tt174*W(29,1)+2*tt55*tt56*tt174*W(29,1)+tt173-tt43&
&*tt45*tt172*W(28,1)+2*tt47*tt48*tt172*W(28,1)+tt171-tt35*tt37*tt1&
&70*W(27,1)+2*tt39*tt40*tt170*W(27,1)+tt169-tt27*tt29*tt168*W(26,1&
&)+2*tt31*tt32*tt168*W(26,1)+tt167-tt14*tt16*tt166*W(25,1)+2*tt22*&
&tt23*tt166*W(25,1)+tt165
hes(5,6) = tt177
hes(6,1) = tt120
hes(6,2) = tt130
hes(6,3) = tt154
hes(6,4) = tt162
hes(6,5) = tt177
hes(6,6) = (-tt59*tt61*tt183*W(30,1))+2*tt63*tt64*tt183*W(30,1)+t&
&t175-tt51*tt53*tt182*W(29,1)+2*tt55*tt56*tt182*W(29,1)+tt173-tt43&
&*tt45*tt181*W(28,1)+2*tt47*tt48*tt181*W(28,1)+tt171-tt35*tt37*tt1&
&80*W(27,1)+2*tt39*tt40*tt180*W(27,1)+tt169-tt27*tt29*tt179*W(26,1&
&)+2*tt31*tt32*tt179*W(26,1)+tt167-tt14*tt16*tt178*W(25,1)+2*tt22*&
&tt23*tt178*W(25,1)+tt165
END 
