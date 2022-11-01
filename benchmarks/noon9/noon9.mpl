kernelopts(printbytes=false, assertlevel=1):
interface(echo=0, prettyprint=0):

read "../custom_f4.mpl";
infolevel[Groebner]:= 5;
# noon9
# sparse all-around challenge
var := [x1,x2,x3,x4,x5,x6,x7,x8,x9]:
sys := [
10*x1^2*x9+10*x2^2*x9+10*x3^2*x9+10*x4^2*x9+10*x5^2*x9+10*x6^2*x9+10*x7^2*x9+10*x8^2*x9-11*x9+10,
10*x1^2*x8+10*x2^2*x8+10*x3^2*x8+10*x4^2*x8+10*x5^2*x8+10*x6^2*x8+10*x7^2*x8+10*x8*x9^2-11*x8+10,
10*x1^2*x7+10*x2^2*x7+10*x3^2*x7+10*x4^2*x7+10*x5^2*x7+10*x6^2*x7+10*x7*x8^2+10*x7*x9^2-11*x7+10,
10*x1^2*x6+10*x2^2*x6+10*x3^2*x6+10*x4^2*x6+10*x5^2*x6+10*x6*x7^2+10*x6*x8^2+10*x6*x9^2-11*x6+10,
10*x1^2*x5+10*x2^2*x5+10*x3^2*x5+10*x4^2*x5+10*x5*x6^2+10*x5*x7^2+10*x5*x8^2+10*x5*x9^2-11*x5+10,
10*x1^2*x4+10*x2^2*x4+10*x3^2*x4+10*x4*x5^2+10*x4*x6^2+10*x4*x7^2+10*x4*x8^2+10*x4*x9^2-11*x4+10,
10*x1^2*x3+10*x2^2*x3+10*x3*x4^2+10*x3*x5^2+10*x3*x6^2+10*x3*x7^2+10*x3*x8^2+10*x3*x9^2-11*x3+10,
10*x1*x2^2+10*x1*x3^2+10*x1*x4^2+10*x1*x5^2+10*x1*x6^2+10*x1*x7^2+10*x1*x8^2+10*x1*x9^2-11*x1+10,
10*x1^2*x2+10*x2*x3^2+10*x2*x4^2+10*x2*x5^2+10*x2*x6^2+10*x2*x7^2+10*x2*x8^2+10*x2*x9^2-11*x2+10
]:

out := [entries(BasicF4Alg(sys, var, weights=""), `pairs`)]; # [x4 = 5, x5 = 9, x6 = 5, x7 = 9, x8 = 0, x1 = 13, x2 = 13, x3 = 5]

# print(Groebner:-IsBasis(out, tdeg(op(var))));
_mean := add([seq(rhs(x), x in out)])/numelems(out):
w:= {seq(lhs(v)=lhs(v)^2, v in select(x->rhs(x)<_mean, out))};
print(w):

gb := Groebner:-Basis(subs(w, sys), tdeg(op(var)), characteristic=11863279):

gb := Groebner:-Basis(sys, tdeg(op(var)), characteristic=11863279):