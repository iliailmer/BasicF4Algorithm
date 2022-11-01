read "../custom_f4.mpl":
kernelopts(printbytes=false, assertlevel=1):
interface(echo=0, prettyprint=0):
infolevel[Groebner]:=10:
et_hat:=[170132672640658494-N_0, N_0*P_0*delta_NE+N_0*mu_N+N_1, 77456066668774370-M_0, M_0*mu_M-S_0*delta_LM+M_1, 291461906629641491-S_0-E_0, E_0^2*mu_EE-E_0*P_0*rho_E-N_0*P_0*delta_NE+E_0*delta_EL+E_1, E_0*S_0*mu_LE+S_0^2*mu_LL-S_0*delta_EL+S_0*delta_LM+S_1, -N_1-9317619436030301701513870719338318497041674693045250, (P_0*delta_NE+mu_N)*N_1+N_2+N_0*delta_NE*P_1, E_0*P_0*mu_PE-P_0^2*rho_P+P_0*S_0*mu_PL+P_0*mu_P+P_1, -M_1+9210120666295968615854855753142287, M_1*mu_M-S_1*delta_LM+M_2, -E_1-S_1-7467219314311807835581553993175930554891257251310846, (2*E_0*mu_EE-P_0*rho_E+delta_EL)*E_1+(-N_0*P_1-N_1*P_0)*delta_NE-E_0*P_1*rho_E+E_2, S_0*mu_LE*E_1+(E_0*mu_LE+2*S_0*mu_LL-delta_EL+delta_LM)*S_1+S_2, -N_2+859247441044979630772970431876145056129249442863597393116741724354254791496063463547314, (N_0*P_2+2*N_1*P_1+N_2*P_0)*delta_NE+N_2*mu_N+N_3, (E_0*mu_PE-2*P_0*rho_P+S_0*mu_PL+mu_P)*P_1+(E_1*mu_PE+S_1*mu_PL)*P_0+P_2, -M_2-808890194356936048219052867435122419642096273203373874948249416522532, M_2*mu_M-S_2*delta_LM+M_3, -E_2-S_2+1076721872917539485491561849194099250808339637881991284702347727261355113180676179575222, (-N_0*P_2-2*N_1*P_1-N_2*P_0)*delta_NE+(-E_0*P_2-2*E_1*P_1-E_2*P_0)*rho_E+(2*E_0*mu_EE+delta_EL)*E_2+2*E_1^2*mu_EE+E_3, (E_0*mu_LE+2*S_0*mu_LL-delta_EL+delta_LM)*S_2+(2*E_1*S_1+E_2*S_0)*mu_LE+2*S_1^2*mu_LL+S_3, -N_3-128098605617606621921664242614215750973015972073596482041203971158063942972500260671686628034919157715681723110961788713266, (N_0*P_3+3*N_1*P_2+3*N_2*P_1+N_3*P_0)*delta_NE+N_3*mu_N+N_4, (E_0*mu_PE-2*P_0*rho_P+S_0*mu_PL+mu_P)*P_2+(2*E_1*P_1+E_2*P_0)*mu_PE+(P_0*S_2+2*P_1*S_1)*mu_PL-2*P_1^2*rho_P+P_3, -M_3+67968559215936684748522818787925263745371422991771366220968687280017014518652748642636519489032020170886, M_3*mu_M-S_3*delta_LM+M_4, -E_3-S_3-154975409025696318555096188738313814276566651200240376713150232949975424974004590287932120658234484970425422536614841901334, (-N_0*P_3-3*N_1*P_2-3*N_2*P_1-N_3*P_0)*delta_NE+(-E_0*P_3-3*E_1*P_2-3*E_2*P_1-E_3*P_0)*rho_E+(2*E_0*mu_EE+delta_EL)*E_3+6*E_1*E_2*mu_EE+E_4, (E_0*S_3+3*E_1*S_2+3*E_2*S_1+E_3*S_0)*mu_LE+(2*S_0*mu_LL-delta_EL+delta_LM)*S_3+6*S_1*S_2*mu_LL+S_4, -N_4+25264405343318853141132300874085487089988997569019438380461725384046414668118164314336271532088720651935960350288903593643188050947937236794369593032594946082, (N_0*P_4+4*N_1*P_3+6*N_2*P_2+4*N_3*P_1+N_4*P_0)*delta_NE+N_4*mu_N+N_5, (E_0*P_3+3*E_1*P_2+3*E_2*P_1+E_3*P_0)*mu_PE+(P_0*S_3+3*P_1*S_2+3*P_2*S_1+P_3*S_0)*mu_PL+(-2*P_0*rho_P+mu_P)*P_3-6*P_1*P_2*rho_P+P_4, -M_4-9148391323629461694329671106555469728526104738499334240560421332462356151299855735064164918988949870195453221214120342895487058272095952178, M_4*mu_M-S_4*delta_LM+M_5, -E_4-S_4+28508657388986726093008722637599421387701297974603797662259226357820235203519608094109639013953395561279276952562110186236796128082300993936020634417944461998, (-N_0*P_4-4*N_1*P_3-6*N_2*P_2-4*N_3*P_1-N_4*P_0)*delta_NE+(-E_0*P_4-4*E_1*P_3-6*E_2*P_2-4*E_3*P_1-E_4*P_0)*rho_E+(2*E_0*E_4+8*E_1*E_3+6*E_2^2)*mu_EE+E_4*delta_EL+E_5, (E_0*S_4+4*E_1*S_3+6*E_2*S_2+4*E_3*S_1+E_4*S_0)*mu_LE+(2*S_0*mu_LL-delta_EL+delta_LM)*S_4+(8*S_1*S_3+6*S_2^2)*mu_LL+S_5, -N_5-6176344998261860822436545240670683070784283720074202206013608231456276168528848901610886740940228540258898954588002196686360069953280491119910385634611729356708056045560767427755966375753023122, (N_0*P_5+5*N_1*P_4+10*N_2*P_3+10*N_3*P_2+5*N_4*P_1+N_5*P_0)*delta_NE+N_5*mu_N+N_6, (E_0*P_4+4*E_1*P_3+6*E_2*P_2+4*E_3*P_1+E_4*P_0)*mu_PE+(P_0*S_4+4*P_1*S_3+6*P_2*S_2+4*P_3*S_1+P_4*S_0)*mu_PL+(-2*P_0*rho_P+mu_P)*P_4+(-8*P_1*P_3-6*P_2^2)*rho_P+P_5, -M_5+1783313973395224265365692167152113457208092575975255543003423126589357670541884421695801699567619224482461696584433796096087142555718680031397675405278064715900647181172614882, M_5*mu_M-S_5*delta_LM+M_6, -N_6+1805070239948704112595159752570323152453237130313187397589295687763608463870640374940985898117032904040694100386341009164249009072312065283437624957481664848330049980868140519819752166230991830021653361708915005980500661254520578, -M_6-428601119803332388261897613415252716452740280369446885867164323898867016193656735104015109138158041485349301449992105933450774581171122375624462012330638685833453143004124684538972127165741656526748893881911378, -E_5-S_5-6543615311317796969532470979393948005104897153823784999163718172132697808640725474314121383318494806405857759927162261426783112533631698708564387600559454330520627592445085739073736643064382030, z_aux-1]:
vars:=[N_6, M_6, S_5, P_5, N_5, M_5, E_5, S_4, P_4, N_4, M_4, E_4, S_3, P_3, N_3, M_3, E_3, S_2, P_2, N_2, M_2, E_2, S_1, P_1, N_1, M_1, E_1, S_0, P_0, N_0, M_0, E_0, z_aux, w_aux, delta_EL, delta_LM, delta_NE, mu_EE, mu_LE, mu_LL, mu_M, mu_N, mu_P, mu_PE, mu_PL, rho_E, rho_P]:

# out := [entries(BasicF4Alg(et_hat, vars, weights=""), `pairs`)]: 
# _mean := add([seq(rhs(x), x in out)])/numelems(out):
# w:= {seq(lhs(v)=lhs(v)^2, v in select(x->rhs(x)<_mean, out))}:
printf("\n\nWeights:    %a\n\n", w):
gb:=CodeTools[Usage](Groebner[Basis](et_hat, tdeg(op(vars)), characteristic=11863279),output='all'):

printf("\n%a\n", gb):
# {}
quit: