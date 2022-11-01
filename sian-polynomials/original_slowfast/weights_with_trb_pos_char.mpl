read "../custom_f4.mpl":
kernelopts(printbytes=false, assertlevel=1):
interface(echo=0, prettyprint=0):
infolevel[Groebner]:=10:
et_hat:=[420458555-xC_0, -k2^2*xB_0+xC_1, -eA_0*xA_0-eB*xB_0-eC_0*xC_0+192480162076840014, eA_1, eC_1, k1^2*xA_0+xA_1, -k1^2*xA_0+k2^2*xB_0+xB_1, 394290306-eA_0, 322196644-eC_0, -xC_1+55568366666740760, -k2^2*xB_1+xC_2, -eA_0*xA_1-eA_1*xA_0-eB*xB_1-eC_0*xC_1-eC_1*xC_0+511631736710344140210292, eA_2, eC_2, k1^2*xA_1+xA_2, -k1^2*xA_1+k2^2*xB_1+xB_2, -xC_2-16832259261323345418019280, -k2^2*xB_2+xC_3, -eA_0*xA_2-2*eA_1*xA_1-eA_2*xA_0-eB*xB_2-eC_0*xC_2-2*eC_1*xC_1-eC_2*xC_0-136513472938525199750710316207368, eA_3, eC_3, k1^2*xA_2+xA_3, -k1^2*xA_2+k2^2*xB_2+xB_3, -xC_3+5029875339627846373577924695212320, -eA_0*xA_3-3*eA_1*xA_2-3*eA_2*xA_1-eA_3*xA_0-eB*xB_3-eC_0*xC_3-3*eC_1*xC_2-3*eC_2*xC_1-eC_3*xC_0+27394786562064691023337422523806344893072, -eA_1, -eA_2, -eA_3, -eC_1, -eC_2, -eC_3, z_aux-1]:
vars:=[xC_3, xB_3, xA_3, eC_3, eA_3, xC_2, xB_2, xA_2, eC_2, eA_2, xC_1, xB_1, xA_1, eC_1, eA_1, xC_0, xB_0, xA_0, eC_0, eA_0, z_aux, w_aux, eB, k1, k2]:

# out := [entries(BasicF4Alg(et_hat, vars, weights=""), `pairs`)]: 
# _mean := add([seq(rhs(x), x in out)])/numelems(out):
# w:= {seq(lhs(v)=lhs(v)^2, v in select(x->rhs(x)<_mean, out))}:
# printf("\n\nWeights:    %a\n\n", w):
gb:=CodeTools[Usage](Groebner[Basis](et_hat, tdeg(op(vars)), characteristic=11863279),output='all'):

printf("\n%a\n", gb):
# {k1 = k1^2, k2 = k2^2, xA_0 = xA_0, xA_1 = xA_1, xA_2 = xA_2, xA_3 = xA_3, xB_0 = xB_0, xB_1 = xB_1, xB_2 = xB_2, xB_3 = xB_3, z_aux = z_aux}
quit: