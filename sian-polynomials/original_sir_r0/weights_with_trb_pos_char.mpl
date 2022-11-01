read "../custom_f4.mpl":
kernelopts(printbytes=false, assertlevel=1):
interface(echo=0, prettyprint=0):
infolevel[Groebner]:=10:
et_hat:=[1400629169-In_0, -In_0*S_0^2*b+In_0*g+In_1, -AUX_0*g-b+187546228204439091/407617433*g, AUX_1, -In_1+412768363916792920674516045, (-S_0^2*b+g)*In_1+In_2-In_0*b*S_1^2, In_0*S_0^2*b+S_1^2, -In_2-293316353548240230270086789959835921303809887, (-In_0*S_2^2-2*In_1*S_1^2-In_2*S_0^2)*b+In_2*g+In_3, In_0*S_1^2*b+In_1*S_0^2*b+S_2^2, -In_3-36146365061743518144543713531607467909486763301490878486680771, (-In_0*S_3^2-3*In_1*S_2^2-3*In_2*S_1^2-In_3*S_0^2)*b+In_3*g+In_4, (In_0*S_2^2+2*In_1*S_1^2+In_2*S_0^2)*b+S_3^2, -In_4+547085439059617083593357132947341838062656832140012518755013109390842278905182041, -g*AUX_1, g*z_aux-1]:
vars:=[In_4, In_3, S_3, In_2, S_2, AUX_1, In_1, S_1, AUX_0, In_0, S_0, z_aux, w_aux, b, g]:

# out := [entries(BasicF4Alg(et_hat, vars, weights=""), `pairs`)]: 
# _mean := add([seq(rhs(x), x in out)])/numelems(out):
# w:= {seq(lhs(v)=lhs(v)^2, v in select(x->rhs(x)<_mean, out))}:
# printf("\n\nWeights:    %a\n\n", w):
gb:=CodeTools[Usage](Groebner[Basis](et_hat, tdeg(op(vars)), characteristic=11863279),output='all'):

printf("\n%a\n", gb):
# {In_0 = In_0, In_1 = In_1, In_2 = In_2, In_3 = In_3, In_4 = In_4, S_0 = S_0^2, S_1 = S_1^2, S_2 = S_2^2, S_3 = S_3^2, z_aux = z_aux}
quit: