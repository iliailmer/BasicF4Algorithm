read "../custom_f4.mpl";
infolevel[Groebner]:= 5;
var := [x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12]:
sys := [
x1*x12*x2 + x1*x12 + x10*x11*x12 + x10*x12*x9 + x12*x2*x3 + x12*x3*x4 + x12* x4*x5 + x12*x5*x6 + x12*x6*x7 + x12*x7*x8 + x12*x8*x9 - 1,
x1*x12*x3 + x10*x12*x8 + x11*x12*x9 + x12*x2*x4 + x12*x2 + x12*x3*x5 + x12* x4*x6 + x12*x5*x7 + x12*x6*x8 + x12*x7*x9 - 2,
x1*x12*x4 + x10*x12*x7 + x11*x12*x8 + x12*x2*x5 + x12*x3*x6 + x12*x3 + x12* x4*x7 + x12*x5*x8 + x12*x6*x9 - 3,
x1*x12*x5 + x10*x12*x6 + x11*x12*x7 + x12*x2*x6 + x12*x3*x7 + x12*x4*x8 + x12*x4 + x12*x5*x9 - 4,
x1*x12*x6 + x10*x12*x5 + x11*x12*x6 + x12*x2*x7 + x12*x3*x8 + x12*x4*x9 + x12*x5 - 5,
x1*x12*x7 + x10*x12*x4 + x11*x12*x5 + x12*x2*x8 + x12*x3*x9 + x12*x6 - 6,
x1*x12*x8 + x10*x12*x3 + x11*x12*x4 + x12*x2*x9 + x12*x7 - 7,
x1*x12*x9 + x10*x12*x2 + x11*x12*x3 + x12*x8 - 8,
x1*x10*x12 + x11*x12*x2 + x12*x9 - 9,
x1*x11*x12 + x10*x12 - 10,
x11*x12 - 11,
x1 + x10 + x11 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + 1
]:

out := [entries(BasicF4Alg(sys, var, weights=""), `pairs`)]; # [x4 = 5, x5 = 9, x6 = 5, x7 = 9, x8 = 0, x1 = 13, x2 = 13, x3 = 5]

# print(Groebner:-IsBasis(out, tdeg(op(var))));
_mean := add([seq(rhs(x), x in out)])/numelems(out):
w:= {seq(lhs(v)=lhs(v)^2, v in select(x->rhs(x)=0, out))};
print(w):

gb := Groebner:-Basis(subs(w, sys), tdeg(op(var))):