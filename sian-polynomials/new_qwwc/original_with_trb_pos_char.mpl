read "../custom_f4.mpl":
kernelopts(printbytes=false, assertlevel=1):
interface(echo=0, prettyprint=0):
infolevel[Groebner]:=10:
et_hat:=[676953570707823-x_0, a*x_0-a*y_0-y_0*z_0+x_1, -434409430990738661355525264262-x_1, a*x_1+x_2+(-a-z_0)*y_1-y_0*z_1, -b*x_0-b*y_0+x_0*z_0+y_1, c*z_0+d*w_0-x_0*y_0+z_1, 210239701336440377032089319314007229897659840-x_2, a*x_2-2*z_1*y_1+x_3+(-a-z_0)*y_2-y_0*z_2, (-b+z_0)*x_1-b*y_1+y_2+x_0*z_1, c*z_1+d*w_1-x_0*y_1-x_1*y_0+z_2, -e*z_0+f*w_0-x_0*y_0+w_1, 211703468362889950294305214225245785486202320793544145859412-x_3, a*x_3-3*z_2*y_1-3*z_1*y_2+x_4+(-a-z_0)*y_3-y_0*z_3, 2*z_1*x_1+(-b+z_0)*x_2-b*y_2+y_3+x_0*z_2, c*z_2+d*w_2-x_0*y_2-2*x_1*y_1-x_2*y_0+z_3, -e*z_1+f*w_1-x_0*y_1-x_1*y_0+w_2, -568383971752183866675102812797549366934411173665173691368255557913427319628-x_4, a*x_4-4*z_3*y_1-6*z_2*y_2-4*z_1*y_3+x_5+(-a-z_0)*y_4-y_0*z_4, 3*z_2*x_1+3*z_1*x_2+(-b+z_0)*x_3-b*y_3+y_4+x_0*z_3, c*z_3+d*w_3-x_0*y_3-3*x_1*y_2-3*x_2*y_1-x_3*y_0+z_4, -e*z_2+f*w_2-x_0*y_2-2*x_1*y_1-x_2*y_0+w_3, 160868338866749046230862486035742522786773796821438329932115766569787191261931082320487944-x_5, a*x_5-5*z_4*y_1-10*z_3*y_2-10*z_2*y_3-5*z_1*y_4+x_6+(-a-z_0)*y_5-y_0*z_5, 4*z_3*x_1+6*z_2*x_2+4*z_1*x_3+(-b+z_0)*x_4-b*y_4+y_5+x_0*z_4, c*z_4+d*w_4-x_0*y_4-4*x_1*y_3-6*x_2*y_2-4*x_3*y_1-x_4*y_0+z_5, -e*z_3+f*w_3-x_0*y_3-3*x_1*y_2-3*x_2*y_1-x_3*y_0+w_4, 3392046785251851215549104703131688074856484564217145958816588544822730211251347038892591529599202084939320-x_6, a*x_6-6*z_5*y_1-15*z_4*y_2-20*z_3*y_3-15*z_2*y_4-6*z_1*y_5+x_7+(-a-z_0)*y_6-y_0*z_6, 5*z_4*x_1+10*z_3*x_2+10*z_2*x_3+5*z_1*x_4+(-b+z_0)*x_5-b*y_5+y_6+x_0*z_5, c*z_5+d*w_5-x_0*y_5-5*x_1*y_4-10*x_2*y_3-10*x_3*y_2-5*x_4*y_1-x_5*y_0+z_6, -e*z_4+f*w_4-x_0*y_4-4*x_1*y_3-6*x_2*y_2-4*x_3*y_1-x_4*y_0+w_5, -13373940991802130713166019922013147408203489945156929476920246909968175540248108056230159910648050739178682322258571415984-x_7, a*x_7-7*z_6*y_1-21*z_5*y_2-35*z_4*y_3-35*z_3*y_4-21*z_2*y_5-7*z_1*y_6+x_8+(-a-z_0)*y_7-y_0*z_7, 6*z_5*x_1+15*z_4*x_2+20*z_3*x_3+15*z_2*x_4+6*z_1*x_5+(-b+z_0)*x_6-b*y_6+y_7+x_0*z_6, c*z_6+d*w_6-x_0*y_6-6*x_1*y_5-15*x_2*y_4-20*x_3*y_3-15*x_4*y_2-6*x_5*y_1-x_6*y_0+z_7, -e*z_5+f*w_5-x_0*y_5-5*x_1*y_4-10*x_2*y_3-10*x_3*y_2-5*x_4*y_1-x_5*y_0+w_6, 19520838128380269785572344475603895460903526784292974755255287150119483748786449054447903538595948875676434928930784060360497742107008080-x_8, a*x_8-8*z_7*y_1-28*z_6*y_2-56*z_5*y_3-70*z_4*y_4-56*z_3*y_5-28*z_2*y_6-8*z_1*y_7+x_9+(-a-z_0)*y_8-y_0*z_8, 7*z_6*x_1+21*z_5*x_2+35*z_4*x_3+35*z_3*x_4+21*z_2*x_5+7*z_1*x_6+(-b+z_0)*x_7-b*y_7+y_8+x_0*z_7, c*z_7+d*w_7-x_0*y_7-7*x_1*y_6-21*x_2*y_5-35*x_3*y_4-35*x_4*y_3-21*x_5*y_2-7*x_6*y_1-x_7*y_0+z_8, -e*z_6+f*w_6-x_0*y_6-6*x_1*y_5-15*x_2*y_4-20*x_3*y_3-15*x_4*y_2-6*x_5*y_1-x_6*y_0+w_7, 82361162483862775772977807477189272053885742423513091050345472071115236609582002844850896373957052441602075410922994433468222500460860402618849911155024-x_9, z_aux-1]:
vars:=[x_9, z_8, y_8, x_8, z_7, y_7, x_7, w_7, z_6, y_6, x_6, w_6, z_5, y_5, x_5, w_5, z_4, y_4, x_4, w_4, z_3, y_3, x_3, w_3, z_2, y_2, x_2, w_2, z_1, y_1, x_1, w_1, z_0, y_0, x_0, w_0, z_aux, w_aux, a, b, c, d, e, f]:

# out := [entries(BasicF4Alg(et_hat, vars, weights=""), `pairs`)]: 
# _mean := add([seq(rhs(x), x in out)])/numelems(out):
# w:= {seq(lhs(v)=lhs(v)^2, v in select(x->rhs(x)<_mean, out))}:
printf("\n\nWeights:    %a\n\n", w):
gb:=CodeTools[Usage](Groebner[Basis](et_hat, tdeg(op(vars)), characteristic=11863279),output='all'):

printf("\n%a\n", gb):
# {}
quit: