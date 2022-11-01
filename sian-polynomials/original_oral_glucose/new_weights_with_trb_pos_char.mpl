read "../custom_f4.mpl":
kernelopts(printbytes=false, assertlevel=1):
interface(echo=0, prettyprint=0):
infolevel[Groebner]:=10:
et_hat:=[435620177-G_0, G_0*X_0+G_0*p1-Gb_0*p1-R_0*v+G_1, 838422005-Ib_0, Ib_1, 1129443801-Gb_0, Gb_1, -G_1+252426677469780452, (p1+X_0)*G_1+G_2-p1*Gb_1-v*R_1+G_0*X_1, R_1-k, Ib_0*p3+X_0*p2+X_1-1254509621*p3, -G_2-207594809671662063111718668, 2*X_1*G_1+(p1+X_0)*G_2+G_3-p1*Gb_2-v*R_2+G_0*X_2, Gb_2, R_2, (-781502031+Ib_1)*p3+p2*X_1+X_2, -G_3+63144106005154002667250965271305402, 3*X_2*G_1+3*X_1*G_2+(p1+X_0)*G_3+G_4-p1*Gb_3-v*R_3+G_0*X_3, Gb_3, R_3, (-465240455+Ib_2)*p3+p2*X_2+X_3, Ib_2, -G_4+249364193432466386530856940979326343803899122, 4*X_3*G_1+6*X_2*G_2+4*X_1*G_3+(p1+X_0)*G_4+G_5-p1*Gb_4-v*R_4+G_0*X_4, Gb_4, R_4, (-870002259+Ib_3)*p3+p2*X_3+X_4, Ib_3, -G_5-762386449299869026948005654193954086404152648336121238, 5*X_4*G_1+10*X_3*G_2+10*X_2*G_3+5*X_1*G_4+(p1+X_0)*G_5+G_6-p1*Gb_5-v*R_5+G_0*X_5, Gb_5, R_5, (-147479033+Ib_4)*p3+p2*X_4+X_5, Ib_4, -G_6+1203485281435409151508860425025341998356778487256904453498443742, 6*X_5*G_1+15*X_4*G_2+20*X_3*G_3+15*X_2*G_4+6*X_1*G_5+(p1+X_0)*G_6+G_7-p1*Gb_6-v*R_6+G_0*X_6, Gb_6, R_6, (-40742143+Ib_5)*p3+p2*X_5+X_6, Ib_5, -G_7-116657440864305234326612779992542403516432802935711067526639227011834618, -Ib_1, -Ib_2, -Ib_3, -Ib_4, -Ib_5, -Gb_1, -Gb_2, -Gb_3, -Gb_4, -Gb_5, -Gb_6, z_aux-1]:
vars:=[G_7, Gb_6, X_6, R_6, G_6, Ib_5, Gb_5, X_5, R_5, G_5, Ib_4, Gb_4, X_4, R_4, G_4, Ib_3, Gb_3, X_3, R_3, G_3, Ib_2, Gb_2, X_2, R_2, G_2, Ib_1, Gb_1, X_1, R_1, G_1, Ib_0, Gb_0, X_0, R_0, G_0, z_aux, w_aux, k, p1, p2, p3, v]:

out := [entries(BasicF4Alg(et_hat, vars, weights=""), `pairs`)]: 
_mean := add([seq(rhs(x), x in out)])/numelems(out):
w:= {seq(lhs(v)=lhs(v)^2, v in select(x->rhs(x)<_mean, out))}:
printf("\n\nWeights:    %a\n\n", w):

et_hat := subs(w, et_hat):

gb:=CodeTools[Usage](Groebner[Basis](et_hat, tdeg(op(vars)), characteristic=11863279),output='all'):

printf("\n%a\n", gb):
# {}
quit: