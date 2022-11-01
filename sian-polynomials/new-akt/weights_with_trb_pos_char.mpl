read "../custom_f4.mpl":
kernelopts(printbytes=false, assertlevel=1):
interface(echo=0, prettyprint=0):
infolevel[Groebner]:=10:
et_hat:=[-a1*pEGFR_0-a1*pEGFR_Akt_0+15638839515288668066461703843094128594169605064330793746074578986386550, Akt_0^3*pEGFR_0*reaction_2_k1-pEGFR_Akt_0*reaction_2_k2^3-EGF_EGFR_0^2*reaction_9_k1+pEGFR_0*reaction_4_k1-pEGFR_Akt_0*reaction_3_k1+pEGFR_1, -Akt_0^3*pEGFR_0*reaction_2_k1+pEGFR_Akt_0*reaction_2_k2^3+pEGFR_Akt_0*reaction_3_k1+pEGFR_Akt_1, -a2*pAkt_0-a2*pAkt_S6_0+384166097839744341671193321313930419843203255918094641551877910543875, S6_0^3*pAkt_0*reaction_5_k1-pAkt_S6_0*reaction_5_k2^3+pAkt_0*reaction_7_k1-pAkt_S6_0*reaction_6_k1-pEGFR_Akt_0*reaction_3_k1+pAkt_1, -S6_0^3*pAkt_0*reaction_5_k1+pAkt_S6_0*reaction_5_k2^3+pAkt_S6_0*reaction_6_k1+pAkt_S6_1, -a3*pS6_0+2519858979947415320848169343976054814743658599487587553141397033318930, pS6_0*reaction_8_k1-pAkt_S6_0*reaction_6_k1+pS6_1, (-pEGFR_1-pEGFR_Akt_1)*a1-480911319748682666468098641237663046808436094779706213443412542453316625896808903594249810663509963292625, (Akt_0^3*pEGFR_1+Akt_1^3*pEGFR_0)*reaction_2_k1+(-reaction_2_k2^3-reaction_3_k1)*pEGFR_Akt_1+pEGFR_1*reaction_4_k1-EGF_EGFR_1^2*reaction_9_k1+pEGFR_2, (-Akt_0^3*pEGFR_1-Akt_1^3*pEGFR_0)*reaction_2_k1+(reaction_2_k2^3+reaction_3_k1)*pEGFR_Akt_1+pEGFR_Akt_2, Akt_0^3*pEGFR_0*reaction_2_k1-pEGFR_Akt_0*reaction_2_k2^3+Akt_1^3-pAkt_0*reaction_7_k1, -EGF_EGFR_0^2*reaction_1_k1+EGF_EGFR_0^2*reaction_1_k2+EGF_EGFR_0^2*reaction_9_k1+EGF_EGFR_1^2, (-pAkt_1-pAkt_S6_1)*a2+402085642442990767691257299703180310619176845351520507656612300398428119488331326128364070591679725303825, (S6_0^3*pAkt_1+S6_1^3*pAkt_0)*reaction_5_k1+(-reaction_5_k2^3-reaction_6_k1)*pAkt_S6_1+pAkt_1*reaction_7_k1-reaction_3_k1*pEGFR_Akt_1+pAkt_2, (-S6_0^3*pAkt_1-S6_1^3*pAkt_0)*reaction_5_k1+(reaction_5_k2^3+reaction_6_k1)*pAkt_S6_1+pAkt_S6_2, S6_0^3*pAkt_0*reaction_5_k1-pAkt_S6_0*reaction_5_k2^3+S6_1^3-pS6_0*reaction_8_k1, -a3*pS6_1-52843194490783290813024844599023929906423174205566644075684708924299122239036147606208820316530639796978, pS6_1*reaction_8_k1-pAkt_S6_1*reaction_6_k1+pS6_2, (-pEGFR_2-pEGFR_Akt_2)*a1+3135978861257847544528885722073674866134311317336263605948058872857919220241796315693893983402769162695849303368977218474970760494107538286435497065063676161807162282403831325, (Akt_0^3*pEGFR_2+2*Akt_1^3*pEGFR_1+Akt_2^3*pEGFR_0)*reaction_2_k1+(-reaction_2_k2^3-reaction_3_k1)*pEGFR_Akt_2+pEGFR_2*reaction_4_k1-EGF_EGFR_2^2*reaction_9_k1+pEGFR_3, (-Akt_0^3*pEGFR_2-2*Akt_1^3*pEGFR_1-Akt_2^3*pEGFR_0)*reaction_2_k1+(reaction_2_k2^3+reaction_3_k1)*pEGFR_Akt_2+pEGFR_Akt_3, Akt_0^3*pEGFR_1*reaction_2_k1+Akt_1^3*pEGFR_0*reaction_2_k1-pEGFR_Akt_1*reaction_2_k2^3+Akt_2^3-pAkt_1*reaction_7_k1, (-reaction_1_k1+reaction_1_k2+reaction_9_k1)*EGF_EGFR_1^2+EGF_EGFR_2^2, (-pAkt_2-pAkt_S6_2)*a2+1285680678564151681267033250726662864601697694128975856122238159773839295251373200028222736855514556594979049641120575679627704490789109350822136483521168327974943416240649450, (S6_0^3*pAkt_2+2*S6_1^3*pAkt_1+S6_2^3*pAkt_0)*reaction_5_k1+(-reaction_5_k2^3-reaction_6_k1)*pAkt_S6_2-reaction_3_k1*pEGFR_Akt_2+pAkt_2*reaction_7_k1+pAkt_3, (-S6_0^3*pAkt_2-2*S6_1^3*pAkt_1-S6_2^3*pAkt_0)*reaction_5_k1+(reaction_5_k2^3+reaction_6_k1)*pAkt_S6_2+pAkt_S6_3, S6_0^3*pAkt_1*reaction_5_k1+S6_1^3*pAkt_0*reaction_5_k1-pAkt_S6_1*reaction_5_k2^3+S6_2^3-pS6_1*reaction_8_k1, -a3*pS6_2+291210367933518719734414835955534083877260904197672508116372572109191302626095057445084588977962600269495711526783236129139712591689111643999382922609661272753840307203331092, pS6_2*reaction_8_k1-pAkt_S6_2*reaction_6_k1+pS6_3, (-pEGFR_3-pEGFR_Akt_3)*a1-34678872847570915196393207900353738650719698275306275461965337688487058134393218150843643415760640304733715365580009682647981296763498180125465523165499486434906245693801689429246977798070252168377035068566232298309952280327255073472570365398725, (Akt_0^3*pEGFR_3+3*Akt_1^3*pEGFR_2+3*Akt_2^3*pEGFR_1+Akt_3^3*pEGFR_0)*reaction_2_k1+(-reaction_2_k2^3-reaction_3_k1)*pEGFR_Akt_3-EGF_EGFR_3^2*reaction_9_k1+pEGFR_3*reaction_4_k1+pEGFR_4, (-Akt_0^3*pEGFR_3-3*Akt_1^3*pEGFR_2-3*Akt_2^3*pEGFR_1-Akt_3^3*pEGFR_0)*reaction_2_k1+(reaction_2_k2^3+reaction_3_k1)*pEGFR_Akt_3+pEGFR_Akt_4, (Akt_0^3*pEGFR_2+2*Akt_1^3*pEGFR_1+Akt_2^3*pEGFR_0)*reaction_2_k1-reaction_2_k2^3*pEGFR_Akt_2-pAkt_2*reaction_7_k1+Akt_3^3, (-reaction_1_k1+reaction_1_k2+reaction_9_k1)*EGF_EGFR_2^2+EGF_EGFR_3^2, (-pAkt_3-pAkt_S6_3)*a2-14196822285209630592979825772998839261974726800415934144794140540987517828058337737293523929188690988453591142143755752689438459819701715926503691730819727802254897427709577072054615607243396516517702636954418816353613513472407945527651724180800, (S6_0^3*pAkt_3+3*S6_1^3*pAkt_2+3*S6_2^3*pAkt_1+S6_3^3*pAkt_0)*reaction_5_k1+(-reaction_5_k2^3-reaction_6_k1)*pAkt_S6_3-reaction_3_k1*pEGFR_Akt_3+pAkt_3*reaction_7_k1+pAkt_4, (-S6_0^3*pAkt_3-3*S6_1^3*pAkt_2-3*S6_2^3*pAkt_1-S6_3^3*pAkt_0)*reaction_5_k1+(reaction_5_k2^3+reaction_6_k1)*pAkt_S6_3+pAkt_S6_4, (S6_0^3*pAkt_2+2*S6_1^3*pAkt_1+S6_2^3*pAkt_0)*reaction_5_k1-reaction_5_k2^3*pAkt_S6_2-reaction_8_k1*pS6_2+S6_3^3, -a3*pS6_3-3007110911406685798471679224550787263948564236142090546083459169721573968219847848804050761445832679678073304574110518339784440523901066574555555659711927315983964490788285659843324454687856860701426438142131556140668071305586202884562487484552, pS6_3*reaction_8_k1-pAkt_S6_3*reaction_6_k1+pS6_4, (-pEGFR_4-pEGFR_Akt_4)*a1+574763317005539309830345010573484610499064412520820496232583839502036062897633372745092941765525185257298518298215551506099477325308940897647787840233572695729798899282159695388087058385174842380894752730655630041741429832870324263733667616201102433734075701703614624225260530520247747320030456252093180920621208725, (Akt_0^3*pEGFR_4+4*Akt_1^3*pEGFR_3+6*Akt_2^3*pEGFR_2+4*Akt_3^3*pEGFR_1+Akt_4^3*pEGFR_0)*reaction_2_k1+(-reaction_2_k2^3-reaction_3_k1)*pEGFR_Akt_4+pEGFR_4*reaction_4_k1-EGF_EGFR_4^2*reaction_9_k1+pEGFR_5, (-Akt_0^3*pEGFR_4-4*Akt_1^3*pEGFR_3-6*Akt_2^3*pEGFR_2-4*Akt_3^3*pEGFR_1-Akt_4^3*pEGFR_0)*reaction_2_k1+(reaction_2_k2^3+reaction_3_k1)*pEGFR_Akt_4+pEGFR_Akt_5, (Akt_0^3*pEGFR_3+3*Akt_1^3*pEGFR_2+3*Akt_2^3*pEGFR_1+Akt_3^3*pEGFR_0)*reaction_2_k1-reaction_2_k2^3*pEGFR_Akt_3-pAkt_3*reaction_7_k1+Akt_4^3, (-reaction_1_k1+reaction_1_k2+reaction_9_k1)*EGF_EGFR_3^2+EGF_EGFR_4^2, (-pAkt_4-pAkt_S6_4)*a2+233907316091622314784054325962672231633607662799813379566362280308576916332690643331306022173833673165448102794078980768723888089747460679620487525905985625121059539565429594749322530515325611253136753266421370062770038185618706499378829194539264364257723896348441149833802631656559648849715297841292359205231546700, (S6_0^3*pAkt_4+4*S6_1^3*pAkt_3+6*S6_2^3*pAkt_2+4*S6_3^3*pAkt_1+S6_4^3*pAkt_0)*reaction_5_k1+(-reaction_5_k2^3-reaction_6_k1)*pAkt_S6_4+pAkt_4*reaction_7_k1-reaction_3_k1*pEGFR_Akt_4+pAkt_5, (-S6_0^3*pAkt_4-4*S6_1^3*pAkt_3-6*S6_2^3*pAkt_2-4*S6_3^3*pAkt_1-S6_4^3*pAkt_0)*reaction_5_k1+(reaction_5_k2^3+reaction_6_k1)*pAkt_S6_4+pAkt_S6_5, (S6_0^3*pAkt_3+3*S6_1^3*pAkt_2+3*S6_2^3*pAkt_1+S6_3^3*pAkt_0)*reaction_5_k1-reaction_5_k2^3*pAkt_S6_3-reaction_8_k1*pS6_3+S6_4^3, -a3*pS6_4+35556167045841694247541498067734590493768751933319854042260320696287876855657827679549441017764664820856194098375464449529083600920473610788072300025144673365559672523721428004916190619259500160237626922816446315804629217608232282567556114497143378456068642114708026219037005903253900275378831967096925785891880448, pS6_4*reaction_8_k1-pAkt_S6_4*reaction_6_k1+pS6_5, (-pEGFR_5-pEGFR_Akt_5)*a1-12701398766509785545694685950522926977597015255789443453687934393048788344232077671299111357640800543010793256324414131860241714875203173229462064532301203465678777498798916155586402918611295956219861607072451262825801212345627337770506966574120351426853675605562945416871905114025647616718430995358868708479868258332252918735256110963821413875222644649946817055105045002984783763667325, (Akt_0^3*pEGFR_5+5*Akt_1^3*pEGFR_4+10*Akt_2^3*pEGFR_3+10*Akt_3^3*pEGFR_2+5*Akt_4^3*pEGFR_1+Akt_5^3*pEGFR_0)*reaction_2_k1+(-reaction_2_k2^3-reaction_3_k1)*pEGFR_Akt_5+pEGFR_5*reaction_4_k1-EGF_EGFR_5^2*reaction_9_k1+pEGFR_6, (-Akt_0^3*pEGFR_5-5*Akt_1^3*pEGFR_4-10*Akt_2^3*pEGFR_3-10*Akt_3^3*pEGFR_2-5*Akt_4^3*pEGFR_1-Akt_5^3*pEGFR_0)*reaction_2_k1+(reaction_2_k2^3+reaction_3_k1)*pEGFR_Akt_5+pEGFR_Akt_6, (Akt_0^3*pEGFR_4+4*Akt_1^3*pEGFR_3+6*Akt_2^3*pEGFR_2+4*Akt_3^3*pEGFR_1+Akt_4^3*pEGFR_0)*reaction_2_k1-reaction_2_k2^3*pEGFR_Akt_4-pAkt_4*reaction_7_k1+Akt_5^3, (-reaction_1_k1+reaction_1_k2+reaction_9_k1)*EGF_EGFR_4^2+EGF_EGFR_5^2, (-pAkt_5-pAkt_S6_5)*a2-5141860124997004660888002851890983848423280754859370619024313272651799263541928224324630903534954274881513465772432530912755279383217255224847960502736750662096154569509594201042590777846610564232689727213622074204630003917923930983285196791013391271955173030165583216883966585334242876197549355969271434954948884571010444404312079312181169393616894420588646548343782005605408712379300, (S6_0^3*pAkt_5+5*S6_1^3*pAkt_4+10*S6_2^3*pAkt_3+10*S6_3^3*pAkt_2+5*S6_4^3*pAkt_1+S6_5^3*pAkt_0)*reaction_5_k1+(-reaction_5_k2^3-reaction_6_k1)*pAkt_S6_5+pAkt_5*reaction_7_k1-reaction_3_k1*pEGFR_Akt_5+pAkt_6, (-S6_0^3*pAkt_5-5*S6_1^3*pAkt_4-10*S6_2^3*pAkt_3-10*S6_3^3*pAkt_2-5*S6_4^3*pAkt_1-S6_5^3*pAkt_0)*reaction_5_k1+(reaction_5_k2^3+reaction_6_k1)*pAkt_S6_5+pAkt_S6_6, (S6_0^3*pAkt_4+4*S6_1^3*pAkt_3+6*S6_2^3*pAkt_2+4*S6_3^3*pAkt_1+S6_4^3*pAkt_0)*reaction_5_k1-reaction_5_k2^3*pAkt_S6_4-reaction_8_k1*pS6_4+S6_5^3, -a3*pS6_5-506689810897205651671567159857092041746480296424452601031213908698731386306034915142039874574695213430217179379116123654599864164748377980522439215960154798977967023773265263975323284962703525545024618079138149206794649940112534803274553243485954032712952521932373491369073331632887508858203862411382516374869517519759005618924415378878711896301891108138420195171553733886816184634496, pS6_5*reaction_8_k1-pAkt_S6_5*reaction_6_k1+pS6_6, (-pEGFR_6-pEGFR_Akt_6)*a1+350852141763489645716592097419448451347796065705638239672120947587632521399091879500543292042718016071351067015923505762473725545260437592963767814091121765383874035426899240484794358413901089895736154130248530479442708905430643275266221217575921105847594760581042083353568930403028590898344754106469053703674664577369501517262212689205851356508361338477860946637402825616241994978732995857980633262866774851650903162379797007573862121409136602272272595125, (Akt_0^3*pEGFR_6+6*Akt_1^3*pEGFR_5+15*Akt_2^3*pEGFR_4+20*Akt_3^3*pEGFR_3+15*Akt_4^3*pEGFR_2+6*Akt_5^3*pEGFR_1+Akt_6^3*pEGFR_0)*reaction_2_k1+(-reaction_2_k2^3-reaction_3_k1)*pEGFR_Akt_6+pEGFR_6*reaction_4_k1-EGF_EGFR_6^2*reaction_9_k1+pEGFR_7, (-Akt_0^3*pEGFR_6-6*Akt_1^3*pEGFR_5-15*Akt_2^3*pEGFR_4-20*Akt_3^3*pEGFR_3-15*Akt_4^3*pEGFR_2-6*Akt_5^3*pEGFR_1-Akt_6^3*pEGFR_0)*reaction_2_k1+(reaction_2_k2^3+reaction_3_k1)*pEGFR_Akt_6+pEGFR_Akt_7, (Akt_0^3*pEGFR_5+5*Akt_1^3*pEGFR_4+10*Akt_2^3*pEGFR_3+10*Akt_3^3*pEGFR_2+5*Akt_4^3*pEGFR_1+Akt_5^3*pEGFR_0)*reaction_2_k1-reaction_2_k2^3*pEGFR_Akt_5-pAkt_5*reaction_7_k1+Akt_6^3, (-reaction_1_k1+reaction_1_k2+reaction_9_k1)*EGF_EGFR_5^2+EGF_EGFR_6^2, (-pAkt_6-pAkt_S6_6)*a2+141535899004004659629999287428974570816134179688673256599907709797775688032849778364906377275293825149122140113222186464536777977556379971075345745852368477600937607169489612823455961859764910902602114158909578348387979386447215097971231166864923184753331699856255311584854166741684085747801935486658819697136772633297867559111448631548400633349819509099260733543441846145466111931359521443716714858551825644882104411534260261034724761043881191230479039200, (S6_0^3*pAkt_6+6*S6_1^3*pAkt_5+15*S6_2^3*pAkt_4+20*S6_3^3*pAkt_3+15*S6_4^3*pAkt_2+6*S6_5^3*pAkt_1+S6_6^3*pAkt_0)*reaction_5_k1+(-reaction_5_k2^3-reaction_6_k1)*pAkt_S6_6+pAkt_6*reaction_7_k1-reaction_3_k1*pEGFR_Akt_6+pAkt_7, (-S6_0^3*pAkt_6-6*S6_1^3*pAkt_5-15*S6_2^3*pAkt_4-20*S6_3^3*pAkt_3-15*S6_4^3*pAkt_2-6*S6_5^3*pAkt_1-S6_6^3*pAkt_0)*reaction_5_k1+(reaction_5_k2^3+reaction_6_k1)*pAkt_S6_6+pAkt_S6_7, (S6_0^3*pAkt_5+5*S6_1^3*pAkt_4+10*S6_2^3*pAkt_3+10*S6_3^3*pAkt_2+5*S6_4^3*pAkt_1+S6_5^3*pAkt_0)*reaction_5_k1-reaction_5_k2^3*pAkt_S6_5-reaction_8_k1*pS6_5+S6_6^3, (-pEGFR_7-pEGFR_Akt_7)*a1-11629952984306971506006899696162124950692223983838682351750959344786203463106319458053544746422626632045071848476378422537863519206022703798656364039714128913914807855730207485402562647733441062210664058703996278089078625911225018871363796000605427202838627830776014134189045885202792490533621062910731986524850611950417702863383046836250460326489578616229666332263987682623392137034001858502942747295287631836455196517818301640991726972475461714383367800278307991162333752858846808890808677515778531151465199253354495057184925, (-pAkt_7-pAkt_S6_7)*a2-4681235279227668711493556187542614377621560929037703117660807697642045590541418212223149115150859894269687101083607804268585585886455185647860837588479443402864153238325170496505935514482150650701663615531547516201080041684697215986539360463653107927930571846302082786566742568296529137897527226304051203577904874620681771160338102663875647960727464602210457338417243183424975294118194425572497904743977140018540035366493386728862103531943490683917194914650992177848046626377871426698196784972221871750216057607311205344353300, -a3*pS6_6+8872716124852577803797777794125596103852543679211894202562854815225720281220636797536874415236977907978433788540004261503453210567654297600252432577522514403478865856677716784448716802377726204314837946924060583186041964826133346025338187387173908056257833883051582634964341302405471819566424828094504627681592781409532275365674918349380443086810624445987722632313151016788944382662600351506229583966191655148087388767813541120118715783153257613991755512, z_aux^2-1]:
vars:=[pEGFR_Akt_7, pAkt_S6_7, pEGFR_7, pAkt_7, pEGFR_Akt_6, EGF_EGFR_6, pAkt_S6_6, pEGFR_6, pAkt_6, pS6_6, Akt_6, S6_6, pEGFR_Akt_5, EGF_EGFR_5, pAkt_S6_5, pEGFR_5, pAkt_5, pS6_5, Akt_5, S6_5, pEGFR_Akt_4, EGF_EGFR_4, pAkt_S6_4, pEGFR_4, pAkt_4, pS6_4, Akt_4, S6_4, pEGFR_Akt_3, EGF_EGFR_3, pAkt_S6_3, pEGFR_3, pAkt_3, pS6_3, Akt_3, S6_3, pEGFR_Akt_2, EGF_EGFR_2, pAkt_S6_2, pEGFR_2, pAkt_2, pS6_2, Akt_2, S6_2, pEGFR_Akt_1, EGF_EGFR_1, pAkt_S6_1, pEGFR_1, pAkt_1, pS6_1, Akt_1, S6_1, pEGFR_Akt_0, EGF_EGFR_0, pAkt_S6_0, pEGFR_0, pAkt_0, pS6_0, Akt_0, S6_0, z_aux, w_aux, EGFR_turnover, a1, a2, a3, reaction_1_k1, reaction_1_k2, reaction_2_k1, reaction_2_k2, reaction_3_k1, reaction_4_k1, reaction_5_k1, reaction_5_k2, reaction_6_k1, reaction_7_k1, reaction_8_k1, reaction_9_k1]:

# out := [entries(BasicF4Alg(et_hat, vars, weights=""), `pairs`)]: 
# _mean := add([seq(rhs(x), x in out)])/numelems(out):
# w:= {seq(lhs(v)=lhs(v)^2, v in select(x->rhs(x)<_mean, out))}:
printf("\n\nWeights:    %a\n\n", w):
gb:=CodeTools[Usage](Groebner[Basis](et_hat, tdeg(op(vars)), characteristic=11863279),output='all'):

printf("\n%a\n", gb):
# {Akt_0 = Akt_0^3, Akt_1 = Akt_1^3, Akt_2 = Akt_2^3, Akt_3 = Akt_3^3, Akt_4 = Akt_4^3, Akt_5 = Akt_5^3, Akt_6 = Akt_6^3, S6_0 = S6_0^3, S6_1 = S6_1^3, S6_2 = S6_2^3, S6_3 = S6_3^3, S6_4 = S6_4^3, S6_5 = S6_5^3, S6_6 = S6_6^3, pAkt_0 = pAkt_0, pAkt_1 = pAkt_1, pAkt_2 = pAkt_2, pAkt_3 = pAkt_3, pAkt_4 = pAkt_4, pAkt_5 = pAkt_5, pAkt_6 = pAkt_6, pAkt_7 = pAkt_7, pEGFR_0 = pEGFR_0, pEGFR_1 = pEGFR_1, pEGFR_2 = pEGFR_2, pEGFR_3 = pEGFR_3, pEGFR_4 = pEGFR_4, pEGFR_5 = pEGFR_5, pEGFR_6 = pEGFR_6, pEGFR_7 = pEGFR_7, pS6_0 = pS6_0, pS6_1 = pS6_1, pS6_2 = pS6_2, pS6_3 = pS6_3, pS6_4 = pS6_4, pS6_5 = pS6_5, pS6_6 = pS6_6, z_aux = z_aux^2, EGF_EGFR_0 = EGF_EGFR_0^2, EGF_EGFR_1 = EGF_EGFR_1^2, EGF_EGFR_2 = EGF_EGFR_2^2, EGF_EGFR_3 = EGF_EGFR_3^2, EGF_EGFR_4 = EGF_EGFR_4^2, EGF_EGFR_5 = EGF_EGFR_5^2, EGF_EGFR_6 = EGF_EGFR_6^2, pAkt_S6_0 = pAkt_S6_0, pAkt_S6_1 = pAkt_S6_1, pAkt_S6_2 = pAkt_S6_2, pAkt_S6_3 = pAkt_S6_3, pAkt_S6_4 = pAkt_S6_4, pAkt_S6_5 = pAkt_S6_5, pAkt_S6_6 = pAkt_S6_6, pAkt_S6_7 = pAkt_S6_7, pEGFR_Akt_0 = pEGFR_Akt_0, pEGFR_Akt_1 = pEGFR_Akt_1, pEGFR_Akt_2 = pEGFR_Akt_2, pEGFR_Akt_3 = pEGFR_Akt_3, pEGFR_Akt_4 = pEGFR_Akt_4, pEGFR_Akt_5 = pEGFR_Akt_5, pEGFR_Akt_6 = pEGFR_Akt_6, pEGFR_Akt_7 = pEGFR_Akt_7, reaction_2_k2 = reaction_2_k2^3, reaction_5_k2 = reaction_5_k2^3}
quit: