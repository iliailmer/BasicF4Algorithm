# bayes148
# sparse with many monomial divisions
kernelopts(printbytes=false, assertlevel=1):
interface(echo=0, prettyprint=0):

read "../custom_f4.mpl";
infolevel[Groebner]:= 5;
var := [x0,x1,x2,x3,x4,x5,x6,x7]:
sys := [x0+x1+x2+x3+x4+x5+x6+x7,
x0*x1+x0*x7+x1*x2+x2*x3+x3*x4+x4*x5+x5*x6+x6*x7,
x0*x1*x2+x0*x1*x7+x0*x6*x7+x1*x2*x3+x2*x3*x4+x3*x4*x5+x4*x5*x6+x5*x6*x7,
x0*x1*x2*x3+x0*x1*x2*x7+x0*x1*x6*x7+x0*x5*x6*x7+x1*x2*x3*x4+x2*x3*x4*x5+x3*x4*x5*x6+x4*x5*x6*x7,
x0*x1*x2*x3*x4+x0*x1*x2*x3*x7+x0*x1*x2*x6*x7+x0*x1*x5*x6*x7+x0*x4*x5*x6*x7+x1*x2*x3*x4*x5+x2*x3*x4*x5*x6+x3*x4*x5*x6*x7,
x0*x1*x2*x3*x4*x5+x0*x1*x2*x3*x4*x7+x0*x1*x2*x3*x6*x7+x0*x1*x2*x5*x6*x7+x0*x1*x4*x5*x6*x7+x0*x3*x4*x5*x6*x7+x1*x2*x3*x4*x5*x6+x2*x3*x4*x5*x6*x7,
x0*x1*x2*x3*x4*x5*x6+x0*x1*x2*x3*x4*x5*x7+x0*x1*x2*x3*x4*x6*x7+x0*x1*x2*x3*x5*x6*x7+x0*x1*x2*x4*x5*x6*x7+x0*x1*x3*x4*x5*x6*x7+x0*x2*x3*x4*x5*x6*x7+x1*x2*x3*x4*x5*x6*x7,
x0*x1*x2*x3*x4*x5*x6*x7-1]:


out := [entries(BasicF4Alg(sys, var, weights=""), `pairs`)]; # [x4 = 5, x5 = 9, x6 = 5, x7 = 9, x8 = 0, x1 = 13, x2 = 13, x3 = 5]

# print(Groebner:-IsBasis(out, tdeg(op(var))));
_mean := add([seq(rhs(x), x in out)])/numelems(out):
w:= {seq(lhs(v)=lhs(v)^2, v in select(x->rhs(x)=0, out))}; #{seq(lhs(v)=lhs(v)^2, v in select(x->rhs(x)<_mean, out))};
print(w):

gb := Groebner:-Basis(subs(w, sys), tdeg(op(var)), characteristic=11863279):

gb := Groebner:-Basis(sys, tdeg(op(var)), characteristic=11863279):