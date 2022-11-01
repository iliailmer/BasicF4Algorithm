# bayes148
# sparse with many monomial divisions
kernelopts(printbytes=false, assertlevel=1):
interface(echo=0, prettyprint=0):

read "../custom_f4.mpl";
infolevel[Groebner]:= 5;
var := [x1,x2,x3,x4,x5,x6,x7]:
sys := [
2*x1^2 - 2*x2^2 + 2*x3^2 - 2*x4^2 + 2*x5^2 - 2*x6^2 + 2*x7^2 - 1,
2*x1^3 - 2*x2^3 + 2*x3^3 - 2*x4^3 + 2*x5^3 - 2*x6^3 + 2*x7^3 - 1,
2*x1^4 - 2*x2^4 + 2*x3^4 - 2*x4^4 + 2*x5^4 - 2*x6^4 + 2*x7^4 - 1,
2*x1^5 - 2*x2^5 + 2*x3^5 - 2*x4^5 + 2*x5^5 - 2*x6^5 + 2*x7^5 - 1,
2*x1^6 - 2*x2^6 + 2*x3^6 - 2*x4^6 + 2*x5^6 - 2*x6^6 + 2*x7^6 - 1,
2*x1^7 - 2*x2^7 + 2*x3^7 - 2*x4^7 + 2*x5^7 - 2*x6^7 + 2*x7^7 - 1,
2*x1^8 - 2*x2^8 + 2*x3^8 - 2*x4^8 + 2*x5^8 - 2*x6^8 + 2*x7^8 - 1
]:

out := [entries(BasicF4Alg(sys, var, weights=""), `pairs`)]; # [x4 = 5, x5 = 9, x6 = 5, x7 = 9, x8 = 0, x1 = 13, x2 = 13, x3 = 5]

# print(Groebner:-IsBasis(out, tdeg(op(var))));
_mean := add([seq(rhs(x), x in out)])/numelems(out):
w:= {seq(lhs(v)=lhs(v)^2, v in select(x->rhs(x)<_mean, out))};
print(w):

gb := Groebner:-Basis(subs(w, sys), tdeg(op(var)), characteristic=11863279):

gb := Groebner:-Basis(sys, tdeg(op(var)), characteristic=11863279):