# yang1

kernelopts(printbytes=false, assertlevel=1):
interface(echo=0, prettyprint=0):

read "../custom_f4.mpl";
infolevel[Groebner]:= 5;

# sparse with many pair updates
var := [x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15, x16, x17, x18, x19, x20, x21, x22, x23, x24, x25, x26, x27, x28, x29, x30, x31, x32, x33, x34, x35, x36, x37, x38, x39, x40, x41, x42, x43, x44, x45, x46, x47, x48]:
sys := [
x21*x45+x22*x46+x23*x47+x24*x48, x17*x45+x18*x46+x19*x47+x20*x48, x13*x45+x14*x46+x15*x47+x16*x48, x10*x46+x11*x47+x12*x48+x45*x9, x45*x5+x46*x6+x47*x7+x48*x8, x1*x45+x2*x46+x3*x47+x4*x48, x21*x41+x22*x42+x23*x43+x24*x44,
x17*x41+x18*x42+x19*x43+x20*x44, x13*x41+x14*x42+x15*x43+x16*x44, x10*x42+x11*x43+x12*x44+x41*x9, x41*x5+x42*x6+x43*x7+x44*x8, x1*x41+x2*x42+x3*x43+x4*x44, x21*x37+x22*x38+x23*x39+x24*x40, x17*x37+x18*x38+x19*x39+x20*x40,
x13*x37+x14*x38+x15*x39+x16*x40, x10*x38+x11*x39+x12*x40+x37*x9, x37*x5+x38*x6+x39*x7+x40*x8, x1*x37+x2*x38+x3*x39+x4*x40, x21*x33+x22*x34+x23*x35+x24*x36, x17*x33+x18*x34+x19*x35+x20*x36, x13*x33+x14*x34+x15*x35+x16*x36,
x10*x34+x11*x35+x12*x36+x33*x9, x33*x5+x34*x6+x35*x7+x36*x8, x1*x33+x2*x34+x3*x35+x36*x4, x21*x29+x22*x30+x23*x31+x24*x32, x17*x29+x18*x30+x19*x31+x20*x32, x13*x29+x14*x30+x15*x31+x16*x32, x10*x30+x11*x31+x12*x32+x29*x9,
x29*x5+x30*x6+x31*x7+x32*x8, x1*x29+x2*x30+x3*x31+x32*x4, x21*x25+x22*x26+x23*x27+x24*x28, x17*x25+x18*x26+x19*x27+x20*x28, x13*x25+x14*x26+x15*x27+x16*x28, x10*x26+x11*x27+x12*x28+x25*x9, x25*x5+x26*x6+x27*x7+x28*x8, x1*x25+x2*x26+x27*x3+x28*x4,
x33*x38*x43*x48-x33*x38*x44*x47-x33*x39*x42*x48+x33*x39*x44*x46+x33*x40*x42*x47-x33*x40*x43*x46-x34*x37*x43*x48+x34*x37*x44*x47+x34*x39*x41*x48-x34*x39*x44*x45-x34*x40*x41*x47+x34*x40*x43*x45+x35*x37*x42*x48-x35*x37*x44*x46-x35*x38*x41*x48+x35*x38*x44*x45+x35*x40*x41*x46-x35*x40*x42*x45-x36*x37*x42*x47+x36*x37*x43*x46+x36*x38*x41*x47-x36*x38*x43*x45-x36*x39*x41*x46+x36*x39*x42*x45,
x29*x38*x43*x48-x29*x38*x44*x47-x29*x39*x42*x48+x29*x39*x44*x46+x29*x40*x42*x47-x29*x40*x43*x46-x30*x37*x43*x48+x30*x37*x44*x47+x30*x39*x41*x48-x30*x39*x44*x45-x30*x40*x41*x47+x30*x40*x43*x45+x31*x37*x42*x48-x31*x37*x44*x46-x31*x38*x41*x48+x31*x38*x44*x45+x31*x40*x41*x46-x31*x40*x42*x45-x32*x37*x42*x47+x32*x37*x43*x46+x32*x38*x41*x47-x32*x38*x43*x45-x32*x39*x41*x46+x32*x39*x42*x45,
x25*x38*x43*x48-x25*x38*x44*x47-x25*x39*x42*x48+x25*x39*x44*x46+x25*x40*x42*x47-x25*x40*x43*x46-x26*x37*x43*x48+x26*x37*x44*x47+x26*x39*x41*x48-x26*x39*x44*x45-x26*x40*x41*x47+x26*x40*x43*x45+x27*x37*x42*x48-x27*x37*x44*x46-x27*x38*x41*x48+x27*x38*x44*x45+x27*x40*x41*x46-x27*x40*x42*x45-x28*x37*x42*x47+x28*x37*x43*x46+x28*x38*x41*x47-x28*x38*x43*x45-x28*x39*x41*x46+x28*x39*x42*x45,
x29*x34*x43*x48-x29*x34*x44*x47-x29*x35*x42*x48+x29*x35*x44*x46+x29*x36*x42*x47-x29*x36*x43*x46-x30*x33*x43*x48+x30*x33*x44*x47+x30*x35*x41*x48-x30*x35*x44*x45-x30*x36*x41*x47+x30*x36*x43*x45+x31*x33*x42*x48-x31*x33*x44*x46-x31*x34*x41*x48+x31*x34*x44*x45+x31*x36*x41*x46-x31*x36*x42*x45-x32*x33*x42*x47+x32*x33*x43*x46+x32*x34*x41*x47-x32*x34*x43*x45-x32*x35*x41*x46+x32*x35*x42*x45,
x25*x34*x43*x48-x25*x34*x44*x47-x25*x35*x42*x48+x25*x35*x44*x46+x25*x36*x42*x47-x25*x36*x43*x46-x26*x33*x43*x48+x26*x33*x44*x47+x26*x35*x41*x48-x26*x35*x44*x45-x26*x36*x41*x47+x26*x36*x43*x45+x27*x33*x42*x48-x27*x33*x44*x46-x27*x34*x41*x48+x27*x34*x44*x45+x27*x36*x41*x46-x27*x36*x42*x45-x28*x33*x42*x47+x28*x33*x43*x46+x28*x34*x41*x47-x28*x34*x43*x45-x28*x35*x41*x46+x28*x35*x42*x45,
x25*x30*x43*x48-x25*x30*x44*x47-x25*x31*x42*x48+x25*x31*x44*x46+x25*x32*x42*x47-x25*x32*x43*x46-x26*x29*x43*x48+x26*x29*x44*x47+x26*x31*x41*x48-x26*x31*x44*x45-x26*x32*x41*x47+x26*x32*x43*x45+x27*x29*x42*x48-x27*x29*x44*x46-x27*x30*x41*x48+x27*x30*x44*x45+x27*x32*x41*x46-x27*x32*x42*x45-x28*x29*x42*x47+x28*x29*x43*x46+x28*x30*x41*x47-x28*x30*x43*x45-x28*x31*x41*x46+x28*x31*x42*x45,
x29*x34*x39*x48-x29*x34*x40*x47-x29*x35*x38*x48+x29*x35*x40*x46+x29*x36*x38*x47-x29*x36*x39*x46-x30*x33*x39*x48+x30*x33*x40*x47+x30*x35*x37*x48-x30*x35*x40*x45-x30*x36*x37*x47+x30*x36*x39*x45+x31*x33*x38*x48-x31*x33*x40*x46-x31*x34*x37*x48+x31*x34*x40*x45+x31*x36*x37*x46-x31*x36*x38*x45-x32*x33*x38*x47+x32*x33*x39*x46+x32*x34*x37*x47-x32*x34*x39*x45-x32*x35*x37*x46+x32*x35*x38*x45,
x25*x34*x39*x48-x25*x34*x40*x47-x25*x35*x38*x48+x25*x35*x40*x46+x25*x36*x38*x47-x25*x36*x39*x46-x26*x33*x39*x48+x26*x33*x40*x47+x26*x35*x37*x48-x26*x35*x40*x45-x26*x36*x37*x47+x26*x36*x39*x45+x27*x33*x38*x48-x27*x33*x40*x46-x27*x34*x37*x48+x27*x34*x40*x45+x27*x36*x37*x46-x27*x36*x38*x45-x28*x33*x38*x47+x28*x33*x39*x46+x28*x34*x37*x47-x28*x34*x39*x45-x28*x35*x37*x46+x28*x35*x38*x45,
x25*x30*x39*x48-x25*x30*x40*x47-x25*x31*x38*x48+x25*x31*x40*x46+x25*x32*x38*x47-x25*x32*x39*x46-x26*x29*x39*x48+x26*x29*x40*x47+x26*x31*x37*x48-x26*x31*x40*x45-x26*x32*x37*x47+x26*x32*x39*x45+x27*x29*x38*x48-x27*x29*x40*x46-x27*x30*x37*x48+x27*x30*x40*x45+x27*x32*x37*x46-x27*x32*x38*x45-x28*x29*x38*x47+x28*x29*x39*x46+x28*x30*x37*x47-x28*x30*x39*x45-x28*x31*x37*x46+x28*x31*x38*x45,
x25*x30*x35*x48-x25*x30*x36*x47-x25*x31*x34*x48+x25*x31*x36*x46+x25*x32*x34*x47-x25*x32*x35*x46-x26*x29*x35*x48+x26*x29*x36*x47+x26*x31*x33*x48-x26*x31*x36*x45-x26*x32*x33*x47+x26*x32*x35*x45+x27*x29*x34*x48-x27*x29*x36*x46-x27*x30*x33*x48+x27*x30*x36*x45+x27*x32*x33*x46-x27*x32*x34*x45-x28*x29*x34*x47+x28*x29*x35*x46+x28*x30*x33*x47-x28*x30*x35*x45-x28*x31*x33*x46+x28*x31*x34*x45,
x29*x34*x39*x44-x29*x34*x40*x43-x29*x35*x38*x44+x29*x35*x40*x42+x29*x36*x38*x43-x29*x36*x39*x42-x30*x33*x39*x44+x30*x33*x40*x43+x30*x35*x37*x44-x30*x35*x40*x41-x30*x36*x37*x43+x30*x36*x39*x41+x31*x33*x38*x44-x31*x33*x40*x42-x31*x34*x37*x44+x31*x34*x40*x41+x31*x36*x37*x42-x31*x36*x38*x41-x32*x33*x38*x43+x32*x33*x39*x42+x32*x34*x37*x43-x32*x34*x39*x41-x32*x35*x37*x42+x32*x35*x38*x41,
x25*x34*x39*x44-x25*x34*x40*x43-x25*x35*x38*x44+x25*x35*x40*x42+x25*x36*x38*x43-x25*x36*x39*x42-x26*x33*x39*x44+x26*x33*x40*x43+x26*x35*x37*x44-x26*x35*x40*x41-x26*x36*x37*x43+x26*x36*x39*x41+x27*x33*x38*x44-x27*x33*x40*x42-x27*x34*x37*x44+x27*x34*x40*x41+x27*x36*x37*x42-x27*x36*x38*x41-x28*x33*x38*x43+x28*x33*x39*x42+x28*x34*x37*x43-x28*x34*x39*x41-x28*x35*x37*x42+x28*x35*x38*x41,
x25*x30*x39*x44-x25*x30*x40*x43-x25*x31*x38*x44+x25*x31*x40*x42+x25*x32*x38*x43-x25*x32*x39*x42-x26*x29*x39*x44+x26*x29*x40*x43+x26*x31*x37*x44-x26*x31*x40*x41-x26*x32*x37*x43+x26*x32*x39*x41+x27*x29*x38*x44-x27*x29*x40*x42-x27*x30*x37*x44+x27*x30*x40*x41+x27*x32*x37*x42-x27*x32*x38*x41-x28*x29*x38*x43+x28*x29*x39*x42+x28*x30*x37*x43-x28*x30*x39*x41-x28*x31*x37*x42+x28*x31*x38*x41,
x25*x30*x35*x44-x25*x30*x36*x43-x25*x31*x34*x44+x25*x31*x36*x42+x25*x32*x34*x43-x25*x32*x35*x42-x26*x29*x35*x44+x26*x29*x36*x43+x26*x31*x33*x44-x26*x31*x36*x41-x26*x32*x33*x43+x26*x32*x35*x41+x27*x29*x34*x44-x27*x29*x36*x42-x27*x30*x33*x44+x27*x30*x36*x41+x27*x32*x33*x42-x27*x32*x34*x41-x28*x29*x34*x43+x28*x29*x35*x42+x28*x30*x33*x43-x28*x30*x35*x41-x28*x31*x33*x42+x28*x31*x34*x41,
x25*x30*x35*x40-x25*x30*x36*x39-x25*x31*x34*x40+x25*x31*x36*x38+x25*x32*x34*x39-x25*x32*x35*x38-x26*x29*x35*x40+x26*x29*x36*x39+x26*x31*x33*x40-x26*x31*x36*x37-x26*x32*x33*x39+x26*x32*x35*x37+x27*x29*x34*x40-x27*x29*x36*x38-x27*x30*x33*x40+x27*x30*x36*x37+x27*x32*x33*x38-x27*x32*x34*x37-x28*x29*x34*x39+x28*x29*x35*x38+x28*x30*x33*x39-x28*x30*x35*x37-x28*x31*x33*x38+x28*x31*x34*x37,
-x10*x13*x19*x24+x10*x13*x20*x23+x10*x15*x17*x24-x10*x15*x20*x21-x10*x16*x17*x23+x10*x16*x19*x21+x11*x13*x18*x24-x11*x13*x20*x22-x11*x14*x17*x24+x11*x14*x20*x21+x11*x16*x17*x22-x11*x16*x18*x21-x12*x13*x18*x23+x12*x13*x19*x22+x12*x14*x17*x23-x12*x14*x19*x21-x12*x15*x17*x22+x12*x15*x18*x21+x14*x19*x24*x9-x14*x20*x23*x9-x15*x18*x24*x9+x15*x20*x22*x9+x16*x18*x23*x9-x16*x19*x22*x9,
-x13*x18*x23*x8+x13*x18*x24*x7+x13*x19*x22*x8-x13*x19*x24*x6-x13*x20*x22*x7+x13*x20*x23*x6+x14*x17*x23*x8-x14*x17*x24*x7-x14*x19*x21*x8+x14*x19*x24*x5+x14*x20*x21*x7-x14*x20*x23*x5-x15*x17*x22*x8+x15*x17*x24*x6+x15*x18*x21*x8-x15*x18*x24*x5-x15*x20*x21*x6+x15*x20*x22*x5+x16*x17*x22*x7-x16*x17*x23*x6-x16*x18*x21*x7+x16*x18*x23*x5+x16*x19*x21*x6-x16*x19*x22*x5,
x1*x14*x19*x24-x1*x14*x20*x23-x1*x15*x18*x24+x1*x15*x20*x22+x1*x16*x18*x23-x1*x16*x19*x22-x13*x18*x23*x4+x13*x18*x24*x3-x13*x19*x2*x24+x13*x19*x22*x4+x13*x2*x20*x23-x13*x20*x22*x3+x14*x17*x23*x4-x14*x17*x24*x3-x14*x19*x21*x4+x14*x20*x21*x3+x15*x17*x2*x24-x15*x17*x22*x4+x15*x18*x21*x4-x15*x2*x20*x21-x16*x17*x2*x23+x16*x17*x22*x3-x16*x18*x21*x3+x16*x19*x2*x21,
x10*x17*x23*x8-x10*x17*x24*x7-x10*x19*x21*x8+x10*x19*x24*x5+x10*x20*x21*x7-x10*x20*x23*x5-x11*x17*x22*x8+x11*x17*x24*x6+x11*x18*x21*x8-x11*x18*x24*x5-x11*x20*x21*x6+x11*x20*x22*x5+x12*x17*x22*x7-x12*x17*x23*x6-x12*x18*x21*x7+x12*x18*x23*x5+x12*x19*x21*x6-x12*x19*x22*x5-x18*x23*x8*x9+x18*x24*x7*x9+x19*x22*x8*x9-x19*x24*x6*x9-x20*x22*x7*x9+x20*x23*x6*x9,
x1*x10*x19*x24-x1*x10*x20*x23-x1*x11*x18*x24+x1*x11*x20*x22+x1*x12*x18*x23-x1*x12*x19*x22+x10*x17*x23*x4-x10*x17*x24*x3-x10*x19*x21*x4+x10*x20*x21*x3+x11*x17*x2*x24-x11*x17*x22*x4+x11*x18*x21*x4-x11*x2*x20*x21-x12*x17*x2*x23+x12*x17*x22*x3-x12*x18*x21*x3+x12*x19*x2*x21-x18*x23*x4*x9+x18*x24*x3*x9-x19*x2*x24*x9+x19*x22*x4*x9+x2*x20*x23*x9-x20*x22*x3*x9,
x1*x18*x23*x8-x1*x18*x24*x7-x1*x19*x22*x8+x1*x19*x24*x6+x1*x20*x22*x7-x1*x20*x23*x6-x17*x2*x23*x8+x17*x2*x24*x7+x17*x22*x3*x8-x17*x22*x4*x7+x17*x23*x4*x6-x17*x24*x3*x6-x18*x21*x3*x8+x18*x21*x4*x7-x18*x23*x4*x5+x18*x24*x3*x5+x19*x2*x21*x8-x19*x2*x24*x5-x19*x21*x4*x6+x19*x22*x4*x5-x2*x20*x21*x7+x2*x20*x23*x5+x20*x21*x3*x6-x20*x22*x3*x5,
x10*x13*x23*x8-x10*x13*x24*x7-x10*x15*x21*x8+x10*x15*x24*x5+x10*x16*x21*x7-x10*x16*x23*x5-x11*x13*x22*x8+x11*x13*x24*x6+x11*x14*x21*x8-x11*x14*x24*x5-x11*x16*x21*x6+x11*x16*x22*x5+x12*x13*x22*x7-x12*x13*x23*x6-x12*x14*x21*x7+x12*x14*x23*x5+x12*x15*x21*x6-x12*x15*x22*x5-x14*x23*x8*x9+x14*x24*x7*x9+x15*x22*x8*x9-x15*x24*x6*x9-x16*x22*x7*x9+x16*x23*x6*x9,
x1*x10*x15*x24-x1*x10*x16*x23-x1*x11*x14*x24+x1*x11*x16*x22+x1*x12*x14*x23-x1*x12*x15*x22+x10*x13*x23*x4-x10*x13*x24*x3-x10*x15*x21*x4+x10*x16*x21*x3+x11*x13*x2*x24-x11*x13*x22*x4+x11*x14*x21*x4-x11*x16*x2*x21-x12*x13*x2*x23+x12*x13*x22*x3-x12*x14*x21*x3+x12*x15*x2*x21-x14*x23*x4*x9+x14*x24*x3*x9-x15*x2*x24*x9+x15*x22*x4*x9+x16*x2*x23*x9-x16*x22*x3*x9,
x1*x14*x23*x8-x1*x14*x24*x7-x1*x15*x22*x8+x1*x15*x24*x6+x1*x16*x22*x7-x1*x16*x23*x6-x13*x2*x23*x8+x13*x2*x24*x7+x13*x22*x3*x8-x13*x22*x4*x7+x13*x23*x4*x6-x13*x24*x3*x6-x14*x21*x3*x8+x14*x21*x4*x7-x14*x23*x4*x5+x14*x24*x3*x5+x15*x2*x21*x8-x15*x2*x24*x5-x15*x21*x4*x6+x15*x22*x4*x5-x16*x2*x21*x7+x16*x2*x23*x5+x16*x21*x3*x6-x16*x22*x3*x5,
x1*x10*x23*x8-x1*x10*x24*x7-x1*x11*x22*x8+x1*x11*x24*x6+x1*x12*x22*x7-x1*x12*x23*x6-x10*x21*x3*x8+x10*x21*x4*x7-x10*x23*x4*x5+x10*x24*x3*x5+x11*x2*x21*x8-x11*x2*x24*x5-x11*x21*x4*x6+x11*x22*x4*x5-x12*x2*x21*x7+x12*x2*x23*x5+x12*x21*x3*x6-x12*x22*x3*x5-x2*x23*x8*x9+x2*x24*x7*x9+x22*x3*x8*x9-x22*x4*x7*x9+x23*x4*x6*x9-x24*x3*x6*x9,
x10*x13*x19*x8-x10*x13*x20*x7-x10*x15*x17*x8+x10*x15*x20*x5+x10*x16*x17*x7-x10*x16*x19*x5-x11*x13*x18*x8+x11*x13*x20*x6+x11*x14*x17*x8-x11*x14*x20*x5-x11*x16*x17*x6+x11*x16*x18*x5+x12*x13*x18*x7-x12*x13*x19*x6-x12*x14*x17*x7+x12*x14*x19*x5+x12*x15*x17*x6-x12*x15*x18*x5-x14*x19*x8*x9+x14*x20*x7*x9+x15*x18*x8*x9-x15*x20*x6*x9-x16*x18*x7*x9+x16*x19*x6*x9,
x1*x10*x15*x20-x1*x10*x16*x19-x1*x11*x14*x20+x1*x11*x16*x18+x1*x12*x14*x19-x1*x12*x15*x18+x10*x13*x19*x4-x10*x13*x20*x3-x10*x15*x17*x4+x10*x16*x17*x3-x11*x13*x18*x4+x11*x13*x2*x20+x11*x14*x17*x4-x11*x16*x17*x2+x12*x13*x18*x3-x12*x13*x19*x2-x12*x14*x17*x3+x12*x15*x17*x2-x14*x19*x4*x9+x14*x20*x3*x9+x15*x18*x4*x9-x15*x2*x20*x9-x16*x18*x3*x9+x16*x19*x2*x9,
x1*x14*x19*x8-x1*x14*x20*x7-x1*x15*x18*x8+x1*x15*x20*x6+x1*x16*x18*x7-x1*x16*x19*x6+x13*x18*x3*x8-x13*x18*x4*x7-x13*x19*x2*x8+x13*x19*x4*x6+x13*x2*x20*x7-x13*x20*x3*x6-x14*x17*x3*x8+x14*x17*x4*x7-x14*x19*x4*x5+x14*x20*x3*x5+x15*x17*x2*x8-x15*x17*x4*x6+x15*x18*x4*x5-x15*x2*x20*x5-x16*x17*x2*x7+x16*x17*x3*x6-x16*x18*x3*x5+x16*x19*x2*x5,
x1*x10*x19*x8-x1*x10*x20*x7-x1*x11*x18*x8+x1*x11*x20*x6+x1*x12*x18*x7-x1*x12*x19*x6-x10*x17*x3*x8+x10*x17*x4*x7-x10*x19*x4*x5+x10*x20*x3*x5+x11*x17*x2*x8-x11*x17*x4*x6+x11*x18*x4*x5-x11*x2*x20*x5-x12*x17*x2*x7+x12*x17*x3*x6-x12*x18*x3*x5+x12*x19*x2*x5+x18*x3*x8*x9-x18*x4*x7*x9-x19*x2*x8*x9+x19*x4*x6*x9+x2*x20*x7*x9-x20*x3*x6*x9,
x1*x10*x15*x8-x1*x10*x16*x7-x1*x11*x14*x8+x1*x11*x16*x6+x1*x12*x14*x7-x1*x12*x15*x6-x10*x13*x3*x8+x10*x13*x4*x7-x10*x15*x4*x5+x10*x16*x3*x5+x11*x13*x2*x8-x11*x13*x4*x6+x11*x14*x4*x5-x11*x16*x2*x5-x12*x13*x2*x7+x12*x13*x3*x6-x12*x14*x3*x5+x12*x15*x2*x5+x14*x3*x8*x9-x14*x4*x7*x9-x15*x2*x8*x9+x15*x4*x6*x9+x16*x2*x7*x9-x16*x3*x6*x9
]:

out := [entries(BasicF4Alg(sys, var, weights=""), `pairs`)]; # [x4 = 5, x5 = 9, x6 = 5, x7 = 9, x8 = 0, x1 = 13, x2 = 13, x3 = 5]

# print(Groebner:-IsBasis(out, tdeg(op(var))));
_mean := add([seq(rhs(x), x in out)])/numelems(out):
w:= {seq(lhs(v)=lhs(v)^2, v in select(x->rhs(x)<_mean, out))};
print(w):

gb := Groebner:-Basis(subs(w, sys), tdeg(op(var)), characteristic=11863279):

gb := Groebner:-Basis(sys, tdeg(op(var)), characteristic=11863279):