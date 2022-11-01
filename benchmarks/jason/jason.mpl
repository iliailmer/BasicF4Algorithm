kernelopts(printbytes=false, assertlevel=1):
interface(echo=0, prettyprint=0):

read "../custom_f4.mpl";
infolevel[Groebner]:= 5;
var := [x1, x2, x3, x4, x5, x6, x7, x8]: 

sys := [x1^2*x3^4+x1*x2*x3^2*x5^2+x1*x2*x3*x4*x5*x7+x1*x2*x3*x4*x6*x8 + x1*x2*x4^2*x6^2+x2^2*x4^4, x2^6, x1^6]:
# pairs_table := GetSelectionTable(sys, tdeg(op(var))):
out := [entries(BasicF4Alg(sys, var, weights=""), `pairs`)]; # [x4 = 5, x5 = 9, x6 = 5, x7 = 9, x8 = 0, x1 = 13, x2 = 13, x3 = 5]

# print(Groebner:-IsBasis(out, tdeg(op(var))));
_mean := add([seq(rhs(x), x in out)])/numelems(out):
w:= {seq(lhs(v)=lhs(v)^2, v in select(x->rhs(x)<_mean, out))};
print(w):

gb := Groebner:-Basis(subs(w, sys), tdeg(op(var)), characteristic=11863279):

gb := Groebner:-Basis(sys, tdeg(op(var)), characteristic=11863279):
# pairs_table_w := GetSelectionTable(subs({x7=x7^2}, sys), tdeg(op(var))):
# printf(`%a`, [entries(pairs_table, `pairs`)]);
# printf(`%a`, [entries(pairs_table_w, `pairs`)]);