read "../custom_f4.mpl";
# infolevel[Groebner]:= 5;
var := [x0,x1,x2,x3,x4,x5]:
sys := [5*x0^2*x3+37*x1*x3*x4+32*x1*x3*x5+21*x3*x5+55*x4*x5,
39*x0*x1*x5+23*x1^2*x4+57*x1*x2*x4+56*x1*x4^2+10*x2^2+52*x3*x4*x5,
33*x0^2*x3+51*x0^2+42*x0*x3*x5+51*x1^2*x4+32*x1*x3^2+x5^3,
44*x0*x3^2+42*x1*x3+47*x1*x4^2+12*x2*x3+2*x2*x4*x5+43*x3*x4^2,
49*x0^2*x2+11*x0*x1*x2+39*x0*x3*x4+44*x0*x3*x4+54*x0*x3+45*x1^2*x4,
48*x0*x2*x3+2*x2^2*x3+59*x2^2*x5+17*x2+36*x3^3+45*x4]:
# pairs_table := GetSelectionTable(sys, tdeg(op(var))):
out := [entries(BasicF4Alg(sys, var, weights=""), `pairs`)]; # [x4 = 5, x5 = 9, x6 = 5, x7 = 9, x8 = 0, x1 = 13, x2 = 13, x3 = 5]

# print(Groebner:-IsBasis(out, tdeg(op(var))));
_mean := add([seq(rhs(x), x in out)])/numelems(out):
w:= {seq(lhs(v)=lhs(v)^2, v in select(x->rhs(x)<_mean, out))};
print(w):

gb := Groebner:-Basis(subs(w, sys), tdeg(op(var))):
# pairs_table_w := GetSelectionTable(subs({x7=x7^2}, sys), tdeg(op(var))):
# printf(`%a`, [entries(pairs_table, `pairs`)]);
# printf(`%a`, [entries(pairs_table_w, `pairs`)]);