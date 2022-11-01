read "../custom_f4.mpl":
kernelopts(printbytes=false, assertlevel=1):
interface(echo=0, prettyprint=0):
infolevel[Groebner]:=10:
et_hat:=[876909521211110465952102503707379663570321-Dd_0-R_0-Tt_0, Dd_0*eta+Dd_0*rho-In_0*eps+Dd_1, -A_0*th-Dd_0*eta+R_0*ksi+R_0*nu+R_1, -A_0*mu-R_0*nu+Tt_0*sgm+Tt_0*ta+Tt_1, -Dd_1-R_1-Tt_1-147466189730264241810337336256757823822034831154709753585280397205259075474601814372, (eta+rho)*Dd_1+Dd_2-eps*In_1, -A_1*th-Dd_1*eta+(nu+ksi)*R_1+R_2, -A_1*mu-R_1*nu+(sgm+ta)*Tt_1+Tt_2, A_0*kappa+A_0*mu+A_0*th-In_0*zeta+A_1, -A_0*S_0*g-Dd_0*S_0*b-In_0*S_0*a-R_0*S_0*dlt+In_0*eps+In_0*lam+In_0*zeta+In_1, -Dd_2-R_2-Tt_2+223229342151986639502132434321298248073817809309131748715417891823349053461003995465797552184882303367447359484941969859334360813164151101636781825403379295797263896328, (eta+rho)*Dd_2+Dd_3-eps*In_2, -A_2*th-Dd_2*eta+(nu+ksi)*R_2+R_3, -A_2*mu-R_2*nu+(sgm+ta)*Tt_2+Tt_3, (th+mu+kappa)*A_1+A_2-zeta*In_1, (-S_0*a+eps+lam+zeta)*In_1+(-A_1*g-Dd_1*b-R_1*dlt)*S_0+(-A_0*g-Dd_0*b-In_0*a-R_0*dlt)*S_1+In_2, A_0*S_0*g+Dd_0*S_0*b+In_0*S_0*a+R_0*S_0*dlt+S_1, -Dd_3-R_3-Tt_3-233289504866366358487958519265030096330550189197221181937105074938087171656831890504182810304488165394563624695574288992384270734497493198106831601097475099085076204949543620898948139570061908575443590011900577830901768814742254936842871651866174029881, (eta+rho)*Dd_3+Dd_4-eps*In_3, -A_3*th-Dd_3*eta+(nu+ksi)*R_3+R_4, -A_3*mu-R_3*nu+(sgm+ta)*Tt_3+Tt_4, (th+mu+kappa)*A_2+A_3-zeta*In_2, (-S_0*a+eps+lam+zeta)*In_2+(-A_2*g-Dd_2*b-R_2*dlt)*S_0+(-2*A_1*g-2*Dd_1*b-2*In_1*a-2*R_1*dlt)*S_1+(-A_0*g-Dd_0*b-In_0*a-R_0*dlt)*S_2+In_3, (A_1*g+Dd_1*b+In_1*a+R_1*dlt)*S_0+(A_0*g+Dd_0*b+In_0*a+R_0*dlt)*S_1+S_2, -Dd_4-R_4-Tt_4+80464397702305863568840117689783826114503368152812058934403787937068860888706960606572880696562100995289691479931913288080501049209237020821689053519563359223418460575483907397527796098702161932794258759919974670638885478265006327977639798737835045776944353953401959337408743259431690198284697847627048788802492928531436130870848759335, (eta+rho)*Dd_4+Dd_5-eps*In_4, -A_4*th-Dd_4*eta+(nu+ksi)*R_4+R_5, -A_4*mu-R_4*nu+(sgm+ta)*Tt_4+Tt_5, (th+mu+kappa)*A_3+A_4-zeta*In_3, (-S_0*a+eps+lam+zeta)*In_3+(-A_3*g-Dd_3*b-R_3*dlt)*S_0+(-3*A_2*g-3*Dd_2*b-3*In_2*a-3*R_2*dlt)*S_1+(-3*A_1*g-3*Dd_1*b-3*In_1*a-3*R_1*dlt)*S_2+(-A_0*g-Dd_0*b-In_0*a-R_0*dlt)*S_3+In_4, (A_2*g+Dd_2*b+In_2*a+R_2*dlt)*S_0+(2*A_1*g+2*Dd_1*b+2*In_1*a+2*R_1*dlt)*S_1+(A_0*g+Dd_0*b+In_0*a+R_0*dlt)*S_2+S_3, -Dd_5-R_5-Tt_5+428008590611055724686339302827110631245449124219457111643495243731171888697080375465689949016299066246672479508495682170670414015438533195792236219406190477522183647739451707153227091667200111191598567930314415173509769625534277894444859028236660955997537424087224793586818461605885790630325042822356681264135653316153773623925569683256835072650079036193599484938967472266498786752147268695154616624112639215211772240570, (eta+rho)*Dd_5+Dd_6-eps*In_5, -A_5*th-Dd_5*eta+(nu+ksi)*R_5+R_6, -A_5*mu-R_5*nu+(sgm+ta)*Tt_5+Tt_6, (th+mu+kappa)*A_4+A_5-zeta*In_4, (-R_0*S_4-4*R_1*S_3-6*R_2*S_2-4*R_3*S_1-R_4*S_0)*dlt+(-In_0*S_4-4*In_1*S_3-6*In_2*S_2-4*In_3*S_1-In_4*S_0)*a+(-Dd_0*S_4-4*Dd_1*S_3-6*Dd_2*S_2-4*Dd_3*S_1-Dd_4*S_0)*b+(-A_0*S_4-4*A_1*S_3-6*A_2*S_2-4*A_3*S_1-A_4*S_0)*g+(eps+zeta+lam)*In_4+In_5, (A_3*g+Dd_3*b+In_3*a+R_3*dlt)*S_0+(3*A_2*g+3*Dd_2*b+3*In_2*a+3*R_2*dlt)*S_1+(3*A_1*g+3*Dd_1*b+3*In_1*a+3*R_1*dlt)*S_2+(A_0*g+Dd_0*b+In_0*a+R_0*dlt)*S_3+S_4, -Dd_6-R_6-Tt_6-1217980888790561946720622138008389345542771727421917528260291911611637719825696616379634187266872876772138363070095635595343830984373395963484669750003309004313941411589157598723809906851489039214369215174632234022906160783227497704852391469209038679582286598646603837000190744578045508143829401399876449722198175813931138752953156054311084165513034221212204361112521735237410746509601419390617039121584050293729283291566696361332690212588020041727554392856801996768530179854755484101906359932891570160210, (eta+rho)*Dd_6+Dd_7-eps*In_6, -A_6*th-Dd_6*eta+(nu+ksi)*R_6+R_7, -A_6*mu-R_6*nu+(sgm+ta)*Tt_6+Tt_7, (th+mu+kappa)*A_5+A_6-zeta*In_5, (-R_0*S_5-5*R_1*S_4-10*R_2*S_3-10*R_3*S_2-5*R_4*S_1-R_5*S_0)*dlt+(-In_0*S_5-5*In_1*S_4-10*In_2*S_3-10*In_3*S_2-5*In_4*S_1-In_5*S_0)*a+(-Dd_0*S_5-5*Dd_1*S_4-10*Dd_2*S_3-10*Dd_3*S_2-5*Dd_4*S_1-Dd_5*S_0)*b+(-A_0*S_5-5*A_1*S_4-10*A_2*S_3-10*A_3*S_2-5*A_4*S_1-A_5*S_0)*g+(eps+zeta+lam)*In_5+In_6, (R_0*S_4+4*R_1*S_3+6*R_2*S_2+4*R_3*S_1+R_4*S_0)*dlt+(In_0*S_4+4*In_1*S_3+6*In_2*S_2+4*In_3*S_1+In_4*S_0)*a+(Dd_0*S_4+4*Dd_1*S_3+6*Dd_2*S_2+4*Dd_3*S_1+Dd_4*S_0)*b+(A_0*S_4+4*A_1*S_3+6*A_2*S_2+4*A_3*S_1+A_4*S_0)*g+S_5, -Dd_7-R_7-Tt_7+322282730748796397662642875822476230737205835297666604752727752210280328575088283277665236848347205165614047365496409021948581941627895549443780336075148398190169455517261047321737512468225089081856658874099295329475184988071874674729667845959695305173310728422931657090140129809152353854858135526708631482610819251877771580812598330401941586928964489906302005647222275800403850813340242371850982774165157803555032775302174818045140908432174964969091118762884572015954025597070788433682537889400314469726478662894584463665274778371149091773476955656713303990748113569480968850772429356085, (eta+rho)*Dd_7+Dd_8-eps*In_7, -A_7*th-Dd_7*eta+(nu+ksi)*R_7+R_8, -A_7*mu-R_7*nu+(sgm+ta)*Tt_7+Tt_8, (th+mu+kappa)*A_6+A_7-zeta*In_6, (-R_0*S_6-6*R_1*S_5-15*R_2*S_4-20*R_3*S_3-15*R_4*S_2-6*R_5*S_1-R_6*S_0)*dlt+(-In_0*S_6-6*In_1*S_5-15*In_2*S_4-20*In_3*S_3-15*In_4*S_2-6*In_5*S_1-In_6*S_0)*a+(-Dd_0*S_6-6*Dd_1*S_5-15*Dd_2*S_4-20*Dd_3*S_3-15*Dd_4*S_2-6*Dd_5*S_1-Dd_6*S_0)*b+(-A_0*S_6-6*A_1*S_5-15*A_2*S_4-20*A_3*S_3-15*A_4*S_2-6*A_5*S_1-A_6*S_0)*g+(eps+zeta+lam)*In_6+In_7, (R_0*S_5+5*R_1*S_4+10*R_2*S_3+10*R_3*S_2+5*R_4*S_1+R_5*S_0)*dlt+(In_0*S_5+5*In_1*S_4+10*In_2*S_3+10*In_3*S_2+5*In_4*S_1+In_5*S_0)*a+(Dd_0*S_5+5*Dd_1*S_4+10*Dd_2*S_3+10*Dd_3*S_2+5*Dd_4*S_1+Dd_5*S_0)*b+(A_0*S_5+5*A_1*S_4+10*A_2*S_3+10*A_3*S_2+5*A_4*S_1+A_5*S_0)*g+S_6, -Dd_8-R_8-Tt_8+9707571132279639936169756241617242391085981334236330507502332852386894077084067753661672133636965979189965577580882269997575605629176997522961782498707829659095594013456705278836144412401587566354564768125639734881638235866448432123806606663498821459398981032618952969842186157095985771876361252644983080111059736538277328368464841816764889148906841307450119795135876747320398652884642913211665083377769848927118401732700917385613111344472936444441677553075093511734585085311278633542083757662534027621673298354713646551350150169905545612138837766929804743584233464274953202706794890639746465339735053836152014353209582590617814726631682245075167198545912023185728064877405, (eta+rho)*Dd_8+Dd_9-eps*In_8, -A_8*th-Dd_8*eta+(nu+ksi)*R_8+R_9, -A_8*mu-R_8*nu+(sgm+ta)*Tt_8+Tt_9, (th+mu+kappa)*A_7+A_8-zeta*In_7, (-R_0*S_7-7*R_1*S_6-21*R_2*S_5-35*R_3*S_4-35*R_4*S_3-21*R_5*S_2-7*R_6*S_1-R_7*S_0)*dlt+(-In_0*S_7-7*In_1*S_6-21*In_2*S_5-35*In_3*S_4-35*In_4*S_3-21*In_5*S_2-7*In_6*S_1-In_7*S_0)*a+(-Dd_0*S_7-7*Dd_1*S_6-21*Dd_2*S_5-35*Dd_3*S_4-35*Dd_4*S_3-21*Dd_5*S_2-7*Dd_6*S_1-Dd_7*S_0)*b+(-A_0*S_7-7*A_1*S_6-21*A_2*S_5-35*A_3*S_4-35*A_4*S_3-21*A_5*S_2-7*A_6*S_1-A_7*S_0)*g+(eps+zeta+lam)*In_7+In_8, (R_0*S_6+6*R_1*S_5+15*R_2*S_4+20*R_3*S_3+15*R_4*S_2+6*R_5*S_1+R_6*S_0)*dlt+(In_0*S_6+6*In_1*S_5+15*In_2*S_4+20*In_3*S_3+15*In_4*S_2+6*In_5*S_1+In_6*S_0)*a+(Dd_0*S_6+6*Dd_1*S_5+15*Dd_2*S_4+20*Dd_3*S_3+15*Dd_4*S_2+6*Dd_5*S_1+Dd_6*S_0)*b+(A_0*S_6+6*A_1*S_5+15*A_2*S_4+20*A_3*S_3+15*A_4*S_2+6*A_5*S_1+A_6*S_0)*g+S_7, -Dd_9-R_9-Tt_9-35305570051498872901899187810586935726870953102876149564834068957871605234037458699815168339663800411783991062573328140082597007991758397027961428706606090480060672634398380559890823783313337143081350516446855501831485022626154687258721308034889715423591772341060145605853411982742719713004578882252325401453414858184324482712142881551630104960551764160044151449986145685775901820877418406999645239568549291226575711657419964230680388615484767451548391429397572053442321145767338714611917586196788981624453029464594146634412647864964768764815867017233468434776358294332844868696056367556534575220335619159163100843009995277635334992121073792682735991707418508730759204953384383915733634769136267392186957466486945167324666107486836915208271746702992226781960, (eta+rho)*Dd_9+Dd_10-eps*In_9, -A_9*th-Dd_9*eta+(nu+ksi)*R_9+R_10, -A_9*mu-R_9*nu+(sgm+ta)*Tt_9+Tt_10, (th+mu+kappa)*A_8+A_9-zeta*In_8, (-R_0*S_8-8*R_1*S_7-28*R_2*S_6-56*R_3*S_5-70*R_4*S_4-56*R_5*S_3-28*R_6*S_2-8*R_7*S_1-R_8*S_0)*dlt+(-In_0*S_8-8*In_1*S_7-28*In_2*S_6-56*In_3*S_5-70*In_4*S_4-56*In_5*S_3-28*In_6*S_2-8*In_7*S_1-In_8*S_0)*a+(-Dd_0*S_8-8*Dd_1*S_7-28*Dd_2*S_6-56*Dd_3*S_5-70*Dd_4*S_4-56*Dd_5*S_3-28*Dd_6*S_2-8*Dd_7*S_1-Dd_8*S_0)*b+(-A_0*S_8-8*A_1*S_7-28*A_2*S_6-56*A_3*S_5-70*A_4*S_4-56*A_5*S_3-28*A_6*S_2-8*A_7*S_1-A_8*S_0)*g+(eps+zeta+lam)*In_8+In_9, (R_0*S_7+7*R_1*S_6+21*R_2*S_5+35*R_3*S_4+35*R_4*S_3+21*R_5*S_2+7*R_6*S_1+R_7*S_0)*dlt+(In_0*S_7+7*In_1*S_6+21*In_2*S_5+35*In_3*S_4+35*In_4*S_3+21*In_5*S_2+7*In_6*S_1+In_7*S_0)*a+(Dd_0*S_7+7*Dd_1*S_6+21*Dd_2*S_5+35*Dd_3*S_4+35*Dd_4*S_3+21*Dd_5*S_2+7*Dd_6*S_1+Dd_7*S_0)*b+(A_0*S_7+7*A_1*S_6+21*A_2*S_5+35*A_3*S_4+35*A_4*S_3+21*A_5*S_2+7*A_6*S_1+A_7*S_0)*g+S_8, -Dd_10-R_10-Tt_10-16054706713533615494195321297196184193172449393748102134286566825225422714260983758030395377955541926476373456458548085271884014795044975450513067290249789531684644892723430604879949936092148909001453994873479892691646887550141824691451255560073564842523919424877139520537214337270679607973488050443341023922878566577341351700501215798753467190727050851652320540942214847227873685480920225458511682951495852209784021498016041501842575330921463352000878897563065660880659937129252438096841394865903561066031260918826549098492565465638150227397537918215561565947456282985151870449174943389759696074347356064936014053406858334287420348889153082550013095200577896727778010504142002266238249616333308176126463857335057341002847628076774780896179381039247537154485952215284249755192112843215001278413490717759340155256149389251459010720443494700300, Dd_11+(eta+rho)*Dd_10-eps*In_10, -A_10*th-Dd_10*eta+R_11+(nu+ksi)*R_10, -A_10*mu-R_10*nu+Tt_11+(sgm+ta)*Tt_10, (th+mu+kappa)*A_9+A_10-zeta*In_9, (-R_0*S_9-9*R_1*S_8-36*R_2*S_7-84*R_3*S_6-126*R_4*S_5-126*R_5*S_4-84*R_6*S_3-36*R_7*S_2-9*R_8*S_1-R_9*S_0)*dlt+(-In_0*S_9-9*In_1*S_8-36*In_2*S_7-84*In_3*S_6-126*In_4*S_5-126*In_5*S_4-84*In_6*S_3-36*In_7*S_2-9*In_8*S_1-In_9*S_0)*a+(-Dd_0*S_9-9*Dd_1*S_8-36*Dd_2*S_7-84*Dd_3*S_6-126*Dd_4*S_5-126*Dd_5*S_4-84*Dd_6*S_3-36*Dd_7*S_2-9*Dd_8*S_1-Dd_9*S_0)*b+(-A_0*S_9-9*A_1*S_8-36*A_2*S_7-84*A_3*S_6-126*A_4*S_5-126*A_5*S_4-84*A_6*S_3-36*A_7*S_2-9*A_8*S_1-A_9*S_0)*g+(eps+zeta+lam)*In_9+In_10, (R_0*S_8+8*R_1*S_7+28*R_2*S_6+56*R_3*S_5+70*R_4*S_4+56*R_5*S_3+28*R_6*S_2+8*R_7*S_1+R_8*S_0)*dlt+(In_0*S_8+8*In_1*S_7+28*In_2*S_6+56*In_3*S_5+70*In_4*S_4+56*In_5*S_3+28*In_6*S_2+8*In_7*S_1+In_8*S_0)*a+(Dd_0*S_8+8*Dd_1*S_7+28*Dd_2*S_6+56*Dd_3*S_5+70*Dd_4*S_4+56*Dd_5*S_3+28*Dd_6*S_2+8*Dd_7*S_1+Dd_8*S_0)*b+(A_0*S_8+8*A_1*S_7+28*A_2*S_6+56*A_3*S_5+70*A_4*S_4+56*A_5*S_3+28*A_6*S_2+8*A_7*S_1+A_8*S_0)*g+S_9, -Dd_11-R_11-Tt_11+724677074915669440456808252243006906289889471798383809324837722545384199779033755490958812363152517391774847069824284909677263295705347174091734873513763963826635287648047565161640685082363895232181460864150617035490155432301702805416708194520452558680659163623987845536443183373123250500264542060700079829055808333349882715778438617558381087652269397874164864388740352123054384978736183944817465677695980543024529754034432349237021586812976128703681413404742199931024777636975340001960389387092270486002590272062500490978934929118052850488179268156430537032723021615021136802283567208632861977741072455676270133576926424412391415354582382814187510750496253690078556773314535372895148163008702279043908202770364588225700851236007391034146947598937367800039572401937568625525714178650202251905126487097355730198169000211901065367955576688019200831502016371097337140851579222722472457009544665589928946759235311900051516545621955, (eta+rho)*Dd_11+Dd_12-eps*In_11, -A_11*th-Dd_11*eta+(nu+ksi)*R_11+R_12, -A_11*mu-R_11*nu+(sgm+ta)*Tt_11+Tt_12, A_11+(th+mu+kappa)*A_10-zeta*In_10, (-R_0*S_10-10*R_1*S_9-R_10*S_0-45*R_2*S_8-120*R_3*S_7-210*R_4*S_6-252*R_5*S_5-210*R_6*S_4-120*R_7*S_3-45*R_8*S_2-10*R_9*S_1)*dlt+(-In_0*S_10-10*In_1*S_9-In_10*S_0-45*In_2*S_8-120*In_3*S_7-210*In_4*S_6-252*In_5*S_5-210*In_6*S_4-120*In_7*S_3-45*In_8*S_2-10*In_9*S_1)*a+(-Dd_0*S_10-10*Dd_1*S_9-Dd_10*S_0-45*Dd_2*S_8-120*Dd_3*S_7-210*Dd_4*S_6-252*Dd_5*S_5-210*Dd_6*S_4-120*Dd_7*S_3-45*Dd_8*S_2-10*Dd_9*S_1)*b+(-A_0*S_10-10*A_1*S_9-A_10*S_0-45*A_2*S_8-120*A_3*S_7-210*A_4*S_6-252*A_5*S_5-210*A_6*S_4-120*A_7*S_3-45*A_8*S_2-10*A_9*S_1)*g+(eps+zeta+lam)*In_10+In_11, (R_0*S_9+9*R_1*S_8+36*R_2*S_7+84*R_3*S_6+126*R_4*S_5+126*R_5*S_4+84*R_6*S_3+36*R_7*S_2+9*R_8*S_1+R_9*S_0)*dlt+(In_0*S_9+9*In_1*S_8+36*In_2*S_7+84*In_3*S_6+126*In_4*S_5+126*In_5*S_4+84*In_6*S_3+36*In_7*S_2+9*In_8*S_1+In_9*S_0)*a+(Dd_0*S_9+9*Dd_1*S_8+36*Dd_2*S_7+84*Dd_3*S_6+126*Dd_4*S_5+126*Dd_5*S_4+84*Dd_6*S_3+36*Dd_7*S_2+9*Dd_8*S_1+Dd_9*S_0)*b+(A_0*S_9+9*A_1*S_8+36*A_2*S_7+84*A_3*S_6+126*A_4*S_5+126*A_5*S_4+84*A_6*S_3+36*A_7*S_2+9*A_8*S_1+A_9*S_0)*g+S_10, -Dd_12-R_12-Tt_12-2869632190638622693691380752771829842910388600855724963782767876945605328072761516619011144566671411843380556907343572254065249166652775551097423923422862561967093483487653172435888255212760269411887966653985839565808122688505040296691996718773192522706707663150697492914610542531769592119126107703904558514396743365214112032593541192463331288175814295419359118523101784298293126022633532859773072334288130648915564912934666450003119178084010234626687951785494122974587518045120035686482458352260465943650219457951336536716925854667099733337798621857623045933862998541956511845730176483585819901625255210053451215593912010905377963498023064421186057296302427615024094700569200228641779418438706004013945930382251341934211469460377001385509384601131911083996828977324467761856914795422490769302941928859911071156964547931397869749464467118709607085538346640372217143760002580303977505331435743814346513986326538769541197465142950418893266314608822771798861546107118282722280520208543645633096991897960236440280525, (eta+rho)*Dd_12+Dd_13-eps*In_12, -A_12*th-Dd_12*eta+(nu+ksi)*R_12+R_13, -A_12*mu-R_12*nu+(sgm+ta)*Tt_12+Tt_13, (th+mu+kappa)*A_11+A_12-zeta*In_11, (-R_0*S_11-11*R_1*S_10-11*R_10*S_1-R_11*S_0-55*R_2*S_9-165*R_3*S_8-330*R_4*S_7-462*R_5*S_6-462*R_6*S_5-330*R_7*S_4-165*R_8*S_3-55*R_9*S_2)*dlt+(-In_0*S_11-11*In_1*S_10-11*In_10*S_1-In_11*S_0-55*In_2*S_9-165*In_3*S_8-330*In_4*S_7-462*In_5*S_6-462*In_6*S_5-330*In_7*S_4-165*In_8*S_3-55*In_9*S_2)*a+(-Dd_0*S_11-11*Dd_1*S_10-11*Dd_10*S_1-Dd_11*S_0-55*Dd_2*S_9-165*Dd_3*S_8-330*Dd_4*S_7-462*Dd_5*S_6-462*Dd_6*S_5-330*Dd_7*S_4-165*Dd_8*S_3-55*Dd_9*S_2)*b+(-A_0*S_11-11*A_1*S_10-11*A_10*S_1-A_11*S_0-55*A_2*S_9-165*A_3*S_8-330*A_4*S_7-462*A_5*S_6-462*A_6*S_5-330*A_7*S_4-165*A_8*S_3-55*A_9*S_2)*g+(eps+zeta+lam)*In_11+In_12, (R_0*S_10+10*R_1*S_9+R_10*S_0+45*R_2*S_8+120*R_3*S_7+210*R_4*S_6+252*R_5*S_5+210*R_6*S_4+120*R_7*S_3+45*R_8*S_2+10*R_9*S_1)*dlt+(In_0*S_10+10*In_1*S_9+In_10*S_0+45*In_2*S_8+120*In_3*S_7+210*In_4*S_6+252*In_5*S_5+210*In_6*S_4+120*In_7*S_3+45*In_8*S_2+10*In_9*S_1)*a+(Dd_0*S_10+10*Dd_1*S_9+Dd_10*S_0+45*Dd_2*S_8+120*Dd_3*S_7+210*Dd_4*S_6+252*Dd_5*S_5+210*Dd_6*S_4+120*Dd_7*S_3+45*Dd_8*S_2+10*Dd_9*S_1)*b+(A_0*S_10+10*A_1*S_9+A_10*S_0+45*A_2*S_8+120*A_3*S_7+210*A_4*S_6+252*A_5*S_5+210*A_6*S_4+120*A_7*S_3+45*A_8*S_2+10*A_9*S_1)*g+S_11, -Dd_13-R_13-Tt_13-5872466778469181254773551560630882503495443378356278628972965969584741669943669860862884666131887208365937874260739649027719486673544268277896854432964590375332899108092802781511754487511882178202710348449406081872722129275045314078655792667815829514211723149042898745252208878019067143801175305211060286211305516165997896133575793097768803139946324570061385321774839192487159412089659510794291937997145773595392625124764355705833786939802500111429924168796980471434460343537840991619474971973041753896363319106384616607444714170888253361934730652458173114192642173475563326867269694182319656211571148129720607711293122877141326637369831898003460184604513226622125658720349204790377178512834183110325066733330995783913696469870105793220441610519812435587597426654033686321515484267544739839519792620087329486736308441102374057734339887826235590545235709992608315995930382071154743575487791823172560620035013918235294751591648448792705786341502072602894045670436445469425442537848447243552147154438436079401393348374171900404476193999532364669770763413498181560441858508263453439014172497912102570, (eta+rho)*Dd_13+Dd_14-eps*In_13, -A_13*th-Dd_13*eta+(nu+ksi)*R_13+R_14, -A_13*mu-R_13*nu+(sgm+ta)*Tt_13+Tt_14, (th+mu+kappa)*A_12+A_13-zeta*In_12, (-R_0*S_12-12*R_1*S_11-66*R_10*S_2-12*R_11*S_1-R_12*S_0-66*R_2*S_10-220*R_3*S_9-495*R_4*S_8-792*R_5*S_7-924*R_6*S_6-792*R_7*S_5-495*R_8*S_4-220*R_9*S_3)*dlt+(-In_0*S_12-12*In_1*S_11-66*In_10*S_2-12*In_11*S_1-In_12*S_0-66*In_2*S_10-220*In_3*S_9-495*In_4*S_8-792*In_5*S_7-924*In_6*S_6-792*In_7*S_5-495*In_8*S_4-220*In_9*S_3)*a+(-Dd_0*S_12-12*Dd_1*S_11-66*Dd_10*S_2-12*Dd_11*S_1-Dd_12*S_0-66*Dd_2*S_10-220*Dd_3*S_9-495*Dd_4*S_8-792*Dd_5*S_7-924*Dd_6*S_6-792*Dd_7*S_5-495*Dd_8*S_4-220*Dd_9*S_3)*b+(-A_0*S_12-12*A_1*S_11-66*A_10*S_2-12*A_11*S_1-A_12*S_0-66*A_2*S_10-220*A_3*S_9-495*A_4*S_8-792*A_5*S_7-924*A_6*S_6-792*A_7*S_5-495*A_8*S_4-220*A_9*S_3)*g+(eps+zeta+lam)*In_12+In_13, (R_0*S_11+11*R_1*S_10+11*R_10*S_1+R_11*S_0+55*R_2*S_9+165*R_3*S_8+330*R_4*S_7+462*R_5*S_6+462*R_6*S_5+330*R_7*S_4+165*R_8*S_3+55*R_9*S_2)*dlt+(In_0*S_11+11*In_1*S_10+11*In_10*S_1+In_11*S_0+55*In_2*S_9+165*In_3*S_8+330*In_4*S_7+462*In_5*S_6+462*In_6*S_5+330*In_7*S_4+165*In_8*S_3+55*In_9*S_2)*a+(Dd_0*S_11+11*Dd_1*S_10+11*Dd_10*S_1+Dd_11*S_0+55*Dd_2*S_9+165*Dd_3*S_8+330*Dd_4*S_7+462*Dd_5*S_6+462*Dd_6*S_5+330*Dd_7*S_4+165*Dd_8*S_3+55*Dd_9*S_2)*b+(A_0*S_11+11*A_1*S_10+11*A_10*S_1+A_11*S_0+55*A_2*S_9+165*A_3*S_8+330*A_4*S_7+462*A_5*S_6+462*A_6*S_5+330*A_7*S_4+165*A_8*S_3+55*A_9*S_2)*g+S_12, -Dd_14-R_14-Tt_14+123894754587633728197628488389686869572634940204179894835856239208864830607735708276530929979557259183713657090713745136020945697623096338475463853866376829381409370533582244705888680295814690732530344008656596605544195877225813107562709406820796662569130198466211260953691317097710338069598697326209129310547165332276243237325060753175924672541160831762938268730117985193237470471462669022124277602236058660115869933662908691577622158061023256839180653892023814849332535750945157160658140485010879222682702777294106784711701691609504498421619878003396324409432425550205928236946589248609774301982062469326182895214309180221244874295070024089846085347108217178436589420676746549964643487142129501013212548974275691141307873816708185962289921419864074273654283986743325672516959249253445357492344627793387022006119962782113930057908636020822445994545065724771929399410156581704044581480094592322187116479364591669595289906126294529367515095764274433165137258324964890990297642320064739216740602098939851525259629048397392198475237208276942956101788947046538546936770659850530810734321358301961213591092056609758285415082778820516934621770210086021871383502010778926797379666532954330, (eta+rho)*Dd_14+Dd_15-eps*In_14, -A_14*th-Dd_14*eta+(nu+ksi)*R_14+R_15, -A_14*mu-R_14*nu+(sgm+ta)*Tt_14+Tt_15, (th+mu+kappa)*A_13+A_14-zeta*In_13, (-R_0*S_13-13*R_1*S_12-286*R_10*S_3-78*R_11*S_2-13*R_12*S_1-R_13*S_0-78*R_2*S_11-286*R_3*S_10-715*R_4*S_9-1287*R_5*S_8-1716*R_6*S_7-1716*R_7*S_6-1287*R_8*S_5-715*R_9*S_4)*dlt+(-In_0*S_13-13*In_1*S_12-286*In_10*S_3-78*In_11*S_2-13*In_12*S_1-In_13*S_0-78*In_2*S_11-286*In_3*S_10-715*In_4*S_9-1287*In_5*S_8-1716*In_6*S_7-1716*In_7*S_6-1287*In_8*S_5-715*In_9*S_4)*a+(-Dd_0*S_13-13*Dd_1*S_12-286*Dd_10*S_3-78*Dd_11*S_2-13*Dd_12*S_1-Dd_13*S_0-78*Dd_2*S_11-286*Dd_3*S_10-715*Dd_4*S_9-1287*Dd_5*S_8-1716*Dd_6*S_7-1716*Dd_7*S_6-1287*Dd_8*S_5-715*Dd_9*S_4)*b+(-A_0*S_13-13*A_1*S_12-286*A_10*S_3-78*A_11*S_2-13*A_12*S_1-A_13*S_0-78*A_2*S_11-286*A_3*S_10-715*A_4*S_9-1287*A_5*S_8-1716*A_6*S_7-1716*A_7*S_6-1287*A_8*S_5-715*A_9*S_4)*g+(eps+zeta+lam)*In_13+In_14, (R_0*S_12+12*R_1*S_11+66*R_10*S_2+12*R_11*S_1+R_12*S_0+66*R_2*S_10+220*R_3*S_9+495*R_4*S_8+792*R_5*S_7+924*R_6*S_6+792*R_7*S_5+495*R_8*S_4+220*R_9*S_3)*dlt+(In_0*S_12+12*In_1*S_11+66*In_10*S_2+12*In_11*S_1+In_12*S_0+66*In_2*S_10+220*In_3*S_9+495*In_4*S_8+792*In_5*S_7+924*In_6*S_6+792*In_7*S_5+495*In_8*S_4+220*In_9*S_3)*a+(Dd_0*S_12+12*Dd_1*S_11+66*Dd_10*S_2+12*Dd_11*S_1+Dd_12*S_0+66*Dd_2*S_10+220*Dd_3*S_9+495*Dd_4*S_8+792*Dd_5*S_7+924*Dd_6*S_6+792*Dd_7*S_5+495*Dd_8*S_4+220*Dd_9*S_3)*b+(A_0*S_12+12*A_1*S_11+66*A_10*S_2+12*A_11*S_1+A_12*S_0+66*A_2*S_10+220*A_3*S_9+495*A_4*S_8+792*A_5*S_7+924*A_6*S_6+792*A_7*S_5+495*A_8*S_4+220*A_9*S_3)*g+S_13, -Dd_15-R_15-Tt_15-482339900023417484519101253668914869373981148491633840267304811512920732184099188958731525502470718232653381735111685352037136376709793034230517207143204568529456067418545169069258831530934182063745772219522416047445985757443209448762302039730630906734172730148051068639258051325033967036390557787037597308626443669823999339666832632460913402590878110661454676564463741913140343654931633304429215351207349392328703919500852848831338899400252156557600225554716704582541415247934957080539750737229802481388546356373352792583827925606483661081746575631625084768501569200506276263230646858082625493762679789304586792561874650533322329790973421345350020676775391313628918223157706837893123027559207566322173410485146382940256511526395813908520123751722973047107367406370661312230186271130172073467597843056617406669301257131684609415663269239251817219678386407401524521468932247749531918956785101934570132441051516080413126627921410093910147954015617635645166673638260442884650366167919961714400269405156137121776074349350757121191050928858274011359552565652578845480753385079636043140084992188606556880553353112230770145173454260355131940518664901729457803960692562807002391516086802862504047153632017542131386280376487006178708840945544021327728300213895709954993020815, (eta+rho)*Dd_15+Dd_16-eps*In_15, -A_15*th-Dd_15*eta+(nu+ksi)*R_15+R_16, -A_15*mu-R_15*nu+(sgm+ta)*Tt_15+Tt_16, (th+mu+kappa)*A_14+A_15-zeta*In_14, (-R_0*S_14-14*R_1*S_13-1001*R_10*S_4-364*R_11*S_3-91*R_12*S_2-14*R_13*S_1-R_14*S_0-91*R_2*S_12-364*R_3*S_11-1001*R_4*S_10-2002*R_5*S_9-3003*R_6*S_8-3432*R_7*S_7-3003*R_8*S_6-2002*R_9*S_5)*dlt+(-In_0*S_14-14*In_1*S_13-1001*In_10*S_4-364*In_11*S_3-91*In_12*S_2-14*In_13*S_1-In_14*S_0-91*In_2*S_12-364*In_3*S_11-1001*In_4*S_10-2002*In_5*S_9-3003*In_6*S_8-3432*In_7*S_7-3003*In_8*S_6-2002*In_9*S_5)*a+(-Dd_0*S_14-14*Dd_1*S_13-1001*Dd_10*S_4-364*Dd_11*S_3-91*Dd_12*S_2-14*Dd_13*S_1-Dd_14*S_0-91*Dd_2*S_12-364*Dd_3*S_11-1001*Dd_4*S_10-2002*Dd_5*S_9-3003*Dd_6*S_8-3432*Dd_7*S_7-3003*Dd_8*S_6-2002*Dd_9*S_5)*b+(-A_0*S_14-14*A_1*S_13-1001*A_10*S_4-364*A_11*S_3-91*A_12*S_2-14*A_13*S_1-A_14*S_0-91*A_2*S_12-364*A_3*S_11-1001*A_4*S_10-2002*A_5*S_9-3003*A_6*S_8-3432*A_7*S_7-3003*A_8*S_6-2002*A_9*S_5)*g+(eps+zeta+lam)*In_14+In_15, (R_0*S_13+13*R_1*S_12+286*R_10*S_3+78*R_11*S_2+13*R_12*S_1+R_13*S_0+78*R_2*S_11+286*R_3*S_10+715*R_4*S_9+1287*R_5*S_8+1716*R_6*S_7+1716*R_7*S_6+1287*R_8*S_5+715*R_9*S_4)*dlt+(In_0*S_13+13*In_1*S_12+286*In_10*S_3+78*In_11*S_2+13*In_12*S_1+In_13*S_0+78*In_2*S_11+286*In_3*S_10+715*In_4*S_9+1287*In_5*S_8+1716*In_6*S_7+1716*In_7*S_6+1287*In_8*S_5+715*In_9*S_4)*a+(Dd_0*S_13+13*Dd_1*S_12+286*Dd_10*S_3+78*Dd_11*S_2+13*Dd_12*S_1+Dd_13*S_0+78*Dd_2*S_11+286*Dd_3*S_10+715*Dd_4*S_9+1287*Dd_5*S_8+1716*Dd_6*S_7+1716*Dd_7*S_6+1287*Dd_8*S_5+715*Dd_9*S_4)*b+(A_0*S_13+13*A_1*S_12+286*A_10*S_3+78*A_11*S_2+13*A_12*S_1+A_13*S_0+78*A_2*S_11+286*A_3*S_10+715*A_4*S_9+1287*A_5*S_8+1716*A_6*S_7+1716*A_7*S_6+1287*A_8*S_5+715*A_9*S_4)*g+S_14, -Dd_16-R_16-Tt_16-2462307318011844058284498429424258712899024929439971428104704968275174951103091069066712652668448719394868452954859974013407714425189445304645491549732638841876797732868612923733504042819498968408765756028862223178223471566839735179609018044190206285642041103867096157685814347155438831263275717040316866856007545213919493760256020376071821590147713552229288262758676136120759694512505929507235004630118237148415125251504874672077880385792515805609304555305527846435530180956417266747476927365960600449659267979835302002719761214361584681644774267022228726078696229614794233353921396749919786223724796416027849955792949710791695206427208486135661369379303017196647343467242981008889922049375962809416421855728473583526019644367350455190297031839039099880191127381371593358397968031619011060241092886587669896120352652989273410723509395450281986812595882115355452538187000990467152551187441476600822690671572309861771442390012352937104379805374044239616920418656431361508518619261085166945012961954923664254728175004408799641836293984513580539685807484442552419914936657668596756086048440291546365570369959675334737262618776937249783460643537366155341505449950441733625287808178650081123682443061053362228178988670685516009318460722205294456955748086586175554568497784610158075438212501004510437563387166165328297418326408512451936699173184116913096535, (eta+rho)*Dd_16+Dd_17-eps*In_16, -A_16*th-Dd_16*eta+(nu+ksi)*R_16+R_17, -A_16*mu-R_16*nu+(sgm+ta)*Tt_16+Tt_17, (th+mu+kappa)*A_15+A_16-zeta*In_15, (-15*S_14*R_1-105*S_13*R_2-455*S_12*R_3-1365*S_11*R_4-3003*S_10*R_5-5005*S_9*R_6-6435*S_8*R_7-6435*S_7*R_8-5005*S_6*R_9-R_0*S_15-3003*S_5*R_10-1365*S_4*R_11-455*S_3*R_12-105*S_2*R_13-15*S_1*R_14-S_0*R_15)*dlt+(-1365*S_4*In_11-455*S_3*In_12-105*S_2*In_13-15*S_1*In_14-5005*S_6*In_9-3003*S_5*In_10-S_0*In_15-6435*S_7*In_8-In_0*S_15-15*S_14*In_1-105*S_13*In_2-455*S_12*In_3-1365*S_11*In_4-3003*S_10*In_5-5005*S_9*In_6-6435*S_8*In_7)*a+(-1365*S_4*Dd_11-455*S_3*Dd_12-105*S_2*Dd_13-15*S_1*Dd_14-S_0*Dd_15-5005*S_6*Dd_9-3003*S_5*Dd_10-Dd_0*S_15-15*S_14*Dd_1-105*S_13*Dd_2-455*S_12*Dd_3-1365*S_11*Dd_4-3003*S_10*Dd_5-5005*S_9*Dd_6-6435*S_8*Dd_7-6435*S_7*Dd_8)*b+(-3003*A_10*S_5-1365*A_11*S_4-455*A_12*S_3-105*A_13*S_2-15*A_14*S_1-A_15*S_0-1365*A_4*S_11-3003*A_5*S_10-S_15*A_0-15*A_1*S_14-105*A_2*S_13-455*A_3*S_12-5005*S_9*A_6-6435*A_7*S_8-6435*A_8*S_7-5005*A_9*S_6)*g+(eps+zeta+lam)*In_15+In_16, (R_0*S_14+14*R_1*S_13+1001*R_10*S_4+364*R_11*S_3+91*R_12*S_2+14*R_13*S_1+R_14*S_0+91*R_2*S_12+364*R_3*S_11+1001*R_4*S_10+2002*R_5*S_9+3003*R_6*S_8+3432*R_7*S_7+3003*R_8*S_6+2002*R_9*S_5)*dlt+(In_0*S_14+14*In_1*S_13+1001*In_10*S_4+364*In_11*S_3+91*In_12*S_2+14*In_13*S_1+In_14*S_0+91*In_2*S_12+364*In_3*S_11+1001*In_4*S_10+2002*In_5*S_9+3003*In_6*S_8+3432*In_7*S_7+3003*In_8*S_6+2002*In_9*S_5)*a+(Dd_0*S_14+14*Dd_1*S_13+1001*Dd_10*S_4+364*Dd_11*S_3+91*Dd_12*S_2+14*Dd_13*S_1+Dd_14*S_0+91*Dd_2*S_12+364*Dd_3*S_11+1001*Dd_4*S_10+2002*Dd_5*S_9+3003*Dd_6*S_8+3432*Dd_7*S_7+3003*Dd_8*S_6+2002*Dd_9*S_5)*b+(A_0*S_14+14*A_1*S_13+1001*A_10*S_4+364*A_11*S_3+91*A_12*S_2+14*A_13*S_1+A_14*S_0+91*A_2*S_12+364*A_3*S_11+1001*A_4*S_10+2002*A_5*S_9+3003*A_6*S_8+3432*A_7*S_7+3003*A_8*S_6+2002*A_9*S_5)*g+S_15, -Dd_17-R_17-Tt_17+40017286368992248698211175148289777191084640692171392119197744350077745594324845284889014377568904370943096786996879700317033416270223111610963765204809768738060305934352089539510802499088014929558086888031821201830278891132183748940758653240006886640005362196584108506803775002525826589566066190710061983886295644147071225751911941668589042495853442395006105962971966735213849801562401087785351640945520871475737449648108839886650960405153150753257075126648582438482138167819609005958824150093121875363455386401921193252330426837976277164122933900051562044486274622906131618392421880458213234917208996283318594245697144608431167111479408146202511726804434989126848127425596031218001169186017388639727856817030811603869848609634584240520540666438055424221960850882398251788630637791619489620848569265246299857364889337081755442062469392756384907461221011875157763714929795946415169488592832330352420559923578011021755450374617775107271270661879393067429425561074476948946840652961909487187607223264716132057944904838165843300786910747472639897500826521214748494661707315446207703456684440181642230312757861238214730859966967378278072781388134653626438771414379239240245331594354578371039758847793320087268517459085132346707024521343122259430336337054859994718450228524351243683349384758191661228183725932214670135250033904841844108461081281416866018755086146474236845799532927210556326059571231564286635284633515408829561170412134823380, (eta+rho)*Dd_17+Dd_18-eps*In_17, -A_17*th-Dd_17*eta+(nu+ksi)*R_17+R_18, -A_17*mu-R_17*nu+(sgm+ta)*Tt_17+Tt_18, (th+mu+kappa)*A_16+A_17-zeta*In_16, (-16*S_15*R_1-120*S_14*R_2-4368*S_5*R_11-1820*S_4*R_12-560*S_3*R_13-120*S_2*R_14-16*S_1*R_15-S_0*R_16-560*S_13*R_3-1820*S_12*R_4-4368*S_11*R_5-8008*S_10*R_6-11440*S_9*R_7-12870*S_8*R_8-11440*S_7*R_9-8008*S_6*R_10-R_0*S_16)*dlt+(-16*S_15*In_1-120*S_14*In_2-4368*S_5*In_11-1820*S_4*In_12-560*S_3*In_13-120*S_2*In_14-16*S_1*In_15-560*S_13*In_3-1820*S_12*In_4-4368*S_11*In_5-8008*S_10*In_6-11440*S_9*In_7-12870*S_8*In_8-11440*S_7*In_9-8008*S_6*In_10-In_0*S_16-S_0*In_16)*a+(-16*S_15*Dd_1-120*S_14*Dd_2-4368*S_5*Dd_11-1820*S_4*Dd_12-560*S_3*Dd_13-120*S_2*Dd_14-16*S_1*Dd_15-S_0*Dd_16-560*S_13*Dd_3-1820*S_12*Dd_4-4368*S_11*Dd_5-8008*S_10*Dd_6-11440*S_9*Dd_7-12870*S_8*Dd_8-11440*S_7*Dd_9-8008*S_6*Dd_10-Dd_0*S_16)*b+(-8008*A_10*S_6-4368*A_11*S_5-1820*A_12*S_4-560*A_13*S_3-120*A_14*S_2-16*A_15*S_1-A_16*S_0-16*A_1*S_15-120*A_2*S_14-560*A_3*S_13-1820*A_4*S_12-4368*A_5*S_11-8008*S_10*A_6-11440*A_7*S_9-12870*A_8*S_8-11440*A_9*S_7-S_16*A_0)*g+(eps+zeta+lam)*In_16+In_17, (5005*S_9*R_6+6435*S_8*R_7+6435*S_7*R_8+5005*S_6*R_9+3003*S_5*R_10+15*S_14*R_1+105*S_13*R_2+1365*S_4*R_11+455*S_3*R_12+105*S_2*R_13+15*S_1*R_14+S_0*R_15+455*S_12*R_3+R_0*S_15+1365*S_11*R_4+3003*S_10*R_5)*dlt+(15*S_14*In_1+105*S_13*In_2+1365*S_4*In_11+455*S_3*In_12+105*S_2*In_13+15*S_1*In_14+S_0*In_15+455*S_12*In_3+1365*S_11*In_4+3003*S_10*In_5+5005*S_9*In_6+6435*S_8*In_7+6435*S_7*In_8+5005*S_6*In_9+3003*S_5*In_10+In_0*S_15)*a+(3003*S_5*Dd_10+15*S_14*Dd_1+105*S_13*Dd_2+1365*S_4*Dd_11+455*S_3*Dd_12+105*S_2*Dd_13+15*S_1*Dd_14+S_0*Dd_15+455*S_12*Dd_3+1365*S_11*Dd_4+3003*S_10*Dd_5+5005*S_9*Dd_6+6435*S_8*Dd_7+6435*S_7*Dd_8+5005*S_6*Dd_9+Dd_0*S_15)*b+(15*A_1*S_14+105*A_2*S_13+1365*A_11*S_4+455*A_12*S_3+105*A_13*S_2+15*A_14*S_1+A_15*S_0+455*A_3*S_12+1365*A_4*S_11+3003*A_5*S_10+5005*S_9*A_6+6435*A_7*S_8+6435*A_8*S_7+5005*A_9*S_6+3003*A_10*S_5+S_15*A_0)*g+S_16, -Dd_18-R_18-Tt_18-137179092896284468016611879564460088609231070703481637589576014374253733678813714008080870275511044353235534697983763709664037317308341407013593330717727638300625839002524429661576388829094045701058838442379628074303002392506147220013848909338857213442217970159588847388191268363799185995672019875948634876693032036480996371999889683414846817105132188678600369816382999118194389381311906423164833532871610076482524192935479182768625544071387556684259071452089328656609811014039803822471000267648984077108817702699594856859425809035497633857115375268121000949751683460899749884995210938913271548938848971982644458517008366392960639937531519995089062227325954609626167019925883029445384141760438056366948235588153751397936551509696663110604747675774115714795915155257121670653945702621442508966513354469484795919131977206856628399452987106390612538131555971946689824205795766511851632208384373315667812176985895396651088695648193654164505240987241529977108872662017389729843261499130878866405122980567831360520671466194620699078943339408022946952950627351705705304347971481033884246117764058729640028997845207222350156384239554334468396822820979027555630830417548725627985942664345033761850820490910191801302792295394827929115901210750915795511227660903911710640025452644871840913745426691324834639073426277047308103805608498705524902810618877139879099008722861172651463142004119060762802158518559528114628980890488939943630716171708600563112181510669236502838226739903868488648328946558395221467816086158311783599086185120, z_aux-1]:
vars:=[Tt_18, Dd_18, R_18, Tt_17, In_17, Dd_17, R_17, A_17, Tt_16, In_16, Dd_16, S_16, R_16, A_16, Tt_15, In_15, Dd_15, S_15, R_15, A_15, Tt_14, In_14, Dd_14, S_14, R_14, A_14, Tt_13, In_13, Dd_13, S_13, R_13, A_13, Tt_12, In_12, Dd_12, S_12, R_12, A_12, Tt_11, In_11, Dd_11, S_11, R_11, A_11, Tt_10, In_10, Dd_10, S_10, R_10, A_10, Tt_9, In_9, Dd_9, S_9, R_9, A_9, Tt_8, In_8, Dd_8, S_8, R_8, A_8, Tt_7, In_7, Dd_7, S_7, R_7, A_7, Tt_6, In_6, Dd_6, S_6, R_6, A_6, Tt_5, In_5, Dd_5, S_5, R_5, A_5, Tt_4, In_4, Dd_4, S_4, R_4, A_4, Tt_3, In_3, Dd_3, S_3, R_3, A_3, Tt_2, In_2, Dd_2, S_2, R_2, A_2, Tt_1, In_1, Dd_1, S_1, R_1, A_1, Tt_0, In_0, Dd_0, S_0, R_0, A_0, z_aux, w_aux, a, b, dlt, eps, eta, g, kappa, ksi, lam, mu, nu, rho, sgm, ta, th, zeta]:

out := [entries(BasicF4Alg(et_hat, vars, weights=""), `pairs`)]: 
_mean := add([seq(rhs(x), x in out)])/numelems(out):
w:= {seq(lhs(v)=lhs(v)^2, v in select(x->rhs(x)<_mean, out))}:
printf("\n\nWeights:    %a\n\n", w):

et_hat := subs(w, et_hat):

gb:=CodeTools[Usage](Groebner[Basis](et_hat, tdeg(op(vars)), characteristic=11863279),output='all'):

printf("\n%a\n", gb):
# {}
quit: