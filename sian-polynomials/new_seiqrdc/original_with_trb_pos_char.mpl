read "../custom_f4.mpl":
kernelopts(printbytes=false, assertlevel=1):
interface(echo=0, prettyprint=0):
infolevel[Groebner]:=10:
et_hat:=[23648947516418358715858379-c_0, -a*s_0+c_0*mu+c_0*tau0+c_1, -c_1-598134519197371292432711266626076760235795679011941, (mu+tau0)*c_1+c_2-a*s_1, b*i_0*n^2*s_0+a*s_0-mu*n+mu*s_0+s_1, -c_2-3557402871899362507766674025675725024345601093773415407663188868470900459744542125526292815970271840842434572098611444127746021937191986833422759982941, (mu+tau0)*c_2+c_3-a*s_2, b*n^2*s_0*i_1+(b*i_0*n^2+a+mu)*s_1+s_2, i_0*mu*s_0+dlt*i_0-e_0*g+i_1, -c_3+152754802824929938942193523463131543640536603967140049776301357497861411913081145051586605294481570187553223381879177226090678862545603773486673735275625947800204065520264060057448150299616606925033006881219469296040945242291936518042597744737869781179, (mu+tau0)*c_3+c_4-a*s_3, (b*i_0*n^2+a+mu)*s_2+n^2*(2*i_1*s_1+i_2*s_0)*b+s_3, -e_1*g+(mu*s_0+dlt)*i_1+i_2+i_0*mu*s_1, -b*i_0*n*s_0+e_0*g+e_0*mu+e_1, -c_4-6559287948633368603315690386330455173983097288020269582607145261893869713806931305387556086367018918166477682215574044361178959815638016600584283605316228280641907369985473052432728664651892577802778087018131399155239587994481443845056706621511403698621890961877997324622110682112940791169582876272763243507529848471610655207582675107956019895402182541, (mu+tau0)*c_4+c_5-a*s_4, 3*n^2*(s_1*i_2+s_2*i_1+1/3*s_3*i_0+1/3*i_3*s_0)*b+(a+mu)*s_3+s_4, (i_0*s_2+2*i_1*s_1+i_2*s_0)*mu-e_2*g+dlt*i_2+i_3, (g+mu)*e_1-n*(i_0*s_1+i_1*s_0)*b+e_2, -c_5+281655683470695337151751650087079958205414511927839593776283817146149851592428634289215499303458074913469607993687901630457682613129160471168796906250076993895939671759349654420437329225015929459564234028516402721382081442460725271703687668470576239143584612642693232117183096330939164179237444849156112255981020355438580079603116689384777607657487727288428976286748986683054515029353852659264968000710568433515617710186644073183248696739194355735381419, (mu+tau0)*c_5+c_6-a*s_5, n^2*(i_0*s_4+4*i_1*s_3+6*i_2*s_2+4*i_3*s_1+i_4*s_0)*b+(a+mu)*s_4+s_5, (i_0*s_3+3*i_1*s_2+3*i_2*s_1+i_3*s_0)*mu-e_3*g+dlt*i_3+i_4, -n*(i_0*s_2+2*i_1*s_1+i_2*s_0)*b+(g+mu)*e_2+e_3, -c_6-12094288991821582953853195672323493586340647908205245225349318797049853119376010966931554952008378581074023447701290277843587339529132079579155976641240730965825943881278513188748275799415748591730406455792670453215745758578971003606175622776370145845994282239586546838436844741043539668179617249218229665585943766183004044131102264968839322260755338563397096767354554582014897401459310674493614932720193846840364233789033695868126706456830347124000625827084040911944487033358303993996650312683748062017551128764512659918060597203931937212661931276985181, (mu+tau0)*c_6+c_7-a*s_6, n^2*(i_0*s_5+5*i_1*s_4+10*i_2*s_3+10*i_3*s_2+5*i_4*s_1+i_5*s_0)*b+(a+mu)*s_5+s_6, (i_0*s_4+4*i_1*s_3+6*i_2*s_2+4*i_3*s_1+i_4*s_0)*mu-e_4*g+dlt*i_4+i_5, -n*(i_0*s_3+3*i_1*s_2+3*i_2*s_1+i_3*s_0)*b+(g+mu)*e_3+e_4, -c_7+519328509246700389415692050142170706117325637376831766579578815512290603470761279987211763245974389510706893295030317619375512349407369543894754496285354939901294015247764148174597725340570069835619075150916486911427181026904558356697422070513587070870299689021777454107989937942716935045600883363662599016253682421111502506838701444752266543106105369812058491761726188834507342969420531137284327164322982026913882290146304725100348568366057234862238150486178802302854055051879383544069552176615956892288421343522067800090111468759021288395830589563925442572671464183367563249681120245286200995227031986915393650566770618124793938946608830192529333355739, (mu+tau0)*c_7+c_8-a*s_7, 6*n^2*(i_5*s_1+1/6*i_6*s_0+5/2*s_2*i_4+10/3*s_3*i_3+5/2*s_4*i_2+s_5*i_1+1/6*s_6*i_0)*b+(a+mu)*s_6+s_7, (i_0*s_5+5*i_1*s_4+10*i_2*s_3+10*i_3*s_2+5*i_4*s_1+i_5*s_0)*mu-e_5*g+dlt*i_5+i_6, -4*(i_3*s_1+1/4*i_4*s_0+3/2*s_2*i_2+s_3*i_1+1/4*s_4*i_0)*n*b+(g+mu)*e_4+e_5, -c_8-22299955019991543457793045051248421271384011044064764363861962859224580340459326119674200686993911271205039600928230443230028009025208476216470312766934992043780842908903339369366500898199385709686013807131379267561996747371076160101499851212211345303925172399485249914926253193606083644507962141013352714551827064226346732920026273429941065219073384819273141706906598393988175742105816049000361029157628158455639474057282500738614746496530498415239842033887829726327062231659250156624324811288084082904311039713750522843130294180694319324071913226220460872971458259496468311116357695781791111720609135866572235368645992078918015274660894505433558605843926244256399714493516283072248266555436695999521328367232866228214255674255903372110537863927411281581, (mu+tau0)*c_8+c_9-a*s_8, 21*(i_5*s_2+1/3*i_6*s_1+1/21*i_7*s_0+5/3*s_3*i_4+5/3*s_4*i_3+s_5*i_2+1/3*s_6*i_1+1/21*s_7*i_0)*n^2*b+(a+mu)*s_7+s_8, (i_0*s_6+6*i_1*s_5+15*i_2*s_4+20*i_3*s_3+15*i_4*s_2+6*i_5*s_1+i_6*s_0)*mu-e_6*g+dlt*i_6+i_7, -10*n*(i_3*s_2+1/2*i_4*s_1+1/10*i_5*s_0+s_3*i_2+1/2*s_4*i_1+1/10*s_5*i_0)*b+(g+mu)*e_5+e_6, -c_9+957559589045044542070725953505951061133193259183115607704809204960475948849603740161066956872041362761595211896401914241738840522514703930554407016273570098400135362093892903742250488046573066709170268664345639466324202403720341109902643049075805638927654186530940804692306568694383912428587228150624614685318143575865868760242998625324232727508934265308290916087621366715634775214306449421688030227389033105014402621996271410054841509834054281277566160772956039012888433340538329310694814935925089464656214740932684463700365841106585047193037042545285247218419637435449757416383673978036211949129931359817281699382967169353938794918620397162061418598384578240058591806837958231904812116884354876819077024111100726119790221139092429965980340952904407121325351070172752737997885531779927579311788046791671322732567781214254339101148059236025416706103303179, (mu+tau0)*c_9+c_10-a*s_9, 56*n^2*(i_5*s_3+1/2*i_6*s_2+1/7*i_7*s_1+1/56*i_8*s_0+5/4*s_4*i_4+s_5*i_3+1/2*s_6*i_2+1/7*s_7*i_1+1/56*s_8*i_0)*b+(a+mu)*s_8+s_9, (i_0*s_7+7*i_1*s_6+21*i_2*s_5+35*i_3*s_4+35*i_4*s_3+21*i_5*s_2+7*i_6*s_1+i_7*s_0)*mu-e_7*g+dlt*i_7+i_8, -20*n*(s_3*i_3+3/4*s_2*i_4+3/10*i_5*s_1+1/20*i_6*s_0+3/4*s_4*i_2+3/10*s_5*i_1+1/20*s_6*i_0)*b+(g+mu)*e_6+e_7, -c_10-41117588163299456736108486095840334041131739315959109440521341489283902099365627221053922281749364406721191283825836289810405192961150890114602988904218192479953200830931371831290323044521337897060541197632198883132664364190055123301876879250039177077181054787410574901577374581474689375184615910095976284204433467283180763899823520664350509334373893279706227779526389748422128443857810838424316484322990130349758898241627055791568432792711534392400146389049185631232778573891720027225518352114327006354830020609013774448716384582169889966124472341732689374198649835026647393072910883158455693729180531337812975674205903382489275709214438398763696801980237700055535780310845219524314790463407916228545267700804336361439547232628413880198021474680136311007143989751209634566076714900351716005727907381987279940655002426407615642023994335579643991233126493585972288663715549161157005496936939030150060215512544954020936942182978413913816148850957120349319421, c_11+(mu+tau0)*c_10-a*s_10, 126*n^2*(i_5*s_4+2/3*i_6*s_3+2/7*i_7*s_2+1/14*i_8*s_1+1/126*i_9*s_0+s_5*i_4+2/3*s_6*i_3+2/7*s_7*i_2+1/14*s_8*i_1+1/126*s_9*i_0)*b+(a+mu)*s_9+s_10, (i_0*s_8+8*i_1*s_7+28*i_2*s_6+56*i_3*s_5+70*i_4*s_4+56*i_5*s_3+28*i_6*s_2+8*i_7*s_1+i_8*s_0)*mu-e_8*g+dlt*i_8+i_9, -35*n*(s_4*i_3+s_3*i_4+3/5*i_5*s_2+1/5*i_6*s_1+1/35*i_7*s_0+3/5*s_5*i_2+1/5*s_6*i_1+1/35*s_7*i_0)*b+(g+mu)*e_7+e_8, -c_11+1765588351585264604647401574674355307557191600289598468758495342988778304809045258169051299259800725503310319637204552541141626963898101383368722632471329216785229193441846121040268073873970023411570478459601441848835720824686224720584039338355160221405782751104341247319803376993630589118492655465456652294017191989046337613834079958745171201065234013512030046277645640580477487767276494208000057253375181513223325075977865405149691944747515374024284609689905110570805820554484635725544632261746357499651511954498461971897380778383491299478757152856763491352112936399245964316063543724955658987557338425936824365140782889153466925165913339486952794973971611924081673412659210844513140757002572495108587598737175849107619440455210073622238458043179185245741960842987713788025384568722674640488673109130962616111814321348393828818103541601459191561561590634449784252950517531200730332492389561931723647336450710991633382566193742124770950027225208637415046576274006943880662011909816993575565892642114465325470379528764016399074350777973100627008733569507899, z_aux-1]:
vars:=[c_11, s_10, c_10, s_9, i_9, c_9, s_8, i_8, e_8, c_8, s_7, i_7, e_7, c_7, s_6, i_6, e_6, c_6, s_5, i_5, e_5, c_5, s_4, i_4, e_4, c_4, s_3, i_3, e_3, c_3, s_2, i_2, e_2, c_2, s_1, i_1, e_1, c_1, s_0, i_0, e_0, c_0, z_aux, w_aux, a, b, dlt, g, mu, n, tau0]:

# out := [entries(BasicF4Alg(et_hat, vars, weights=""), `pairs`)]: 
# _mean := add([seq(rhs(x), x in out)])/numelems(out):
# w:= {seq(lhs(v)=lhs(v)^2, v in select(x->rhs(x)<_mean, out))}:
# printf("\n\nWeights:    %a\n\n", w):
gb:=CodeTools[Usage](Groebner[Basis](et_hat, tdeg(op(vars)), characteristic=11863279),output='all'):

printf("\n%a\n", gb):
# {}
quit: