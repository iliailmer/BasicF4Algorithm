read "../custom_f4.mpl";
infolevel[Groebner]:= 5;
# katsura9
# dense with many zero-reductions
var := [x0, x1, x2, x3, x4, x5, x6, x7, x8, x9]:
sys := [2*x0*x8+2*x1*x7+2*x1*x9+2*x2*x6+2*x3*x5+x4^2-x8,
2*x0*x7+2*x1*x6+2*x1*x8+2*x2*x5+2*x2*x9+2*x3*x4-x7,
2*x0*x6+2*x1*x5+2*x1*x7+2*x2*x4+2*x2*x8+x3^2+2*x3*x9-x6,
2*x0*x5+2*x1*x4+2*x1*x6+2*x2*x3+2*x2*x7+2*x3*x8+2*x4*x9-x5,
2*x0*x4+2*x1*x3+2*x1*x5+x2^2+2*x2*x6+2*x3*x7+2*x4*x8+2*x5*x9-x4,
2*x0*x3+2*x1*x2+2*x1*x4+2*x2*x5+2*x3*x6+2*x4*x7+2*x5*x8+2*x6*x9-x3,
2*x0*x2+x1^2+2*x1*x3+2*x2*x4+2*x3*x5+2*x4*x6+2*x5*x7+2*x6*x8+2*x7*x9-x2,
2*x0*x1+2*x1*x2+2*x2*x3+2*x3*x4+2*x4*x5+2*x5*x6+2*x6*x7+2*x7*x8+2*x8*x9-x1,
x0^2+2*x1^2+2*x2^2+2*x3^2+2*x4^2+2*x5^2+2*x6^2+2*x7^2+2*x8^2+2*x9^2-x0,
2*x9+2*x8+2*x7+2*x6+2*x5+2*x4+2*x3+2*x2+2*x1+x0-1]:

# pairs_table := GetSelectionTable(sys, tdeg(op(var))):
out := [entries(BasicF4Alg(sys, var, weights=""), `pairs`)]; # [x4 = 5, x5 = 9, x6 = 5, x7 = 9, x8 = 0, x1 = 13, x2 = 13, x3 = 5]

# print(Groebner:-IsBasis(out, tdeg(op(var))));
_mean := add([seq(rhs(x), x in out)])/numelems(out):
w:= {seq(lhs(v)=lhs(v)^2, v in select(x->rhs(x)<_mean, out))};
print(w):

gb := Groebner:-Basis(subs(w, sys), tdeg(op(var))):