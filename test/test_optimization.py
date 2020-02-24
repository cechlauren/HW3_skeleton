import sys
from hw3align.sequences import *
from hw3align.optimalmatrix import * 
from hw3align.optimalmatrix import true_pos_align
from hw3align.optimalmatrix import true_neg_align
from hw3align.getmatrix import blosum50
#######################################################################################################
#To test the optimization of these matrices, need to get true positive and true negative alignments

TP_Alignments, TN_Alignments = getAlignments(scorematrix=blosum50, gap_start=-10, gap_extend= -1)

def test_get_alignments():
	assert TP_Alignments[7] == ("NSNQIKILGNQGSFLTKG-PSKLNDRADSRRSLW--------DQGNFPLIIK------NLKI", "NCSTFYVVKEDGTIVYTGTATSMFD-NDTKETVYIADFSSVNEEGTYYLAVPGVGKSVNFKI")
#making sure the alignment example in my sw testing looks as expected.   

#######################################################################################################
#Now the scoring from sw will be redone

def test_scoring():
	
    testA, testB = TP_Alignments[7]
    assert scoreAlignment(testA, testB, blosum50, -10, -1) == 35.0
#making sure that the scoring example in my sw testing looks as expected
#######################################################################################################
#BRACE YOURSELF!!! Messy but needed to be done. Pulled from jupyter notebook for : getAlignments(blosum50, -10, -1)
true_pos_align = [("SLEAAQKSNVTSSWAKASAAWGTAGPEFFMALFDAHDDVFAKFSGLFSGAAKGTVKNTPEMAAQAQSFKGLVSNWV--DNLDNAGALEGQCKTFAANHKARGISAGQLEAAFKVLSGFMK---SYGG-----DE---GAWTAVAGALMGEI",
   "SLEHAKVD--TSNEARQD------GIDLYKHMFENYPPLRKYFKSREEYTAE-DVQNDPFFAKQGQ--KILLACHVLCATYDDRETFNAYTRELLDRH-AR----DHVHMPPEVWTDFWKLFEEYLGKKTTLDEPTKQAWHEIGREFAKEI"),
  ("AQVKDALTKMRAAALDAQKATPPKLE--DKSP----DSPEMKDFRHGFDILVGQIDDALKL---ANEGKVKEAQAAAEQLKTTRNAYHQKYR",
   "AQIKANVEVLKTLTALPWAAFGPGTEGGDARPEIWSDAASFKQKQQAFQ------DNIVKLSAAADAGDLDKLRAAFGDVGASCKACHDAYR"),
  ("LQEIIKTLNSLTEQKTLCTELTVTDIFAASKNTTEKETFCRAATVLRQFYSHHEKDTRCLG-ATAQQFH-RHKQLIRFLKRLDRNLWGLAGLNSCPVKEANQSTLENFLER",
   "LQMILNGINNYKNPKL--TRMLTFKFYMPKKATELKHLQC-----LEEELKPLEE---VLNLAQSKNFHLRPRDLISNINVIVLEL-------KCEYADETATIVE-FLNR"),
  ("CNACH--GTGLLNAPKVGDSAAWKTRADAKGG--LDGLLAQSLSGLNAMPPKGTCAD-CSDDELKA",
   "CAACHMGGRNSVMPEKTLDKAALEQYLD--GGFKVESIIYQVENGKGAMPAW---ADRLSEEEIQA"),
  ("LTEEQIAEFKEAFALF-DKDGDGTITTKELGTVMR--------SLG-------QNPTEAELQDMINEVDADGNGTIDFPEFLSLMARKMKEQDSEEELIE--------AFKVFDRDGNGLISAAELRHVMTNLGEKLTD-DEVDEMIREADIDGDGHINY------EEFVRMMVS",
   "LNDFQKQKIKFTFDFFLDMNHDGSIQDNDFEDMMTRYKEVNKGSLSDADYKSMQASLEDEWRDLKGRADINKDDVVSWEEYLAMWEKTIATCKSVADLPAWCQNRIPFLFKGMDVSGDGIVDLEEFQNYCKNFQLQCADVPAVYNVI----TDG-GKVTFDLNRYKELYYRLLTS"),
  ("RAECIQR-GVSPSQAQGLGSNLVTE", "RKRKIDRDAVLNMWQQGLGASHISK"),
  ("YEDFQKVYNAIALKLREDDEYDNYIG--YGPVLVRLAWHTSGTWDKHDNTGGSYGGTYRFKKEFNDPSNAG-LQNGFKFLEPIHKEFPWISSGDLFSLGGVTAVQEMQ-GPKIP--------WRC---GRVDTP----EDTTPDNGRLPDADKDA-DYVRTFFQRLNMNDREVVALMGAH--ALGKTHLKNSGYEG--PW------------GAANNVFTNEFYLNLLNEDWKLEKNDANNEQWDSKSGYMMLPTDYSLIQDPKYLSIV-KEYAN--DQDKFFKD----FSKAFEKLLENGITFPK",
   "YQEARKIVGAMVQII----TYRDYLPLVLGPTAMR-----------------KYLPTYRSYNDSVDPRIANVFTNAFRYGHTLIQPF-------MFRLD--NRYQPMEPNPRVPLSRVFFASWRVVLEGGID-PILRGLMATP--AKLNRQNQIAVDEIR---ERL----FEQVMRIGLDLPALNMQRSRDHGLPGYNAWRRFCGLPQPETVGQLGTVLRN---LKLARKLMEQYGTPNNIDIW---MGGVSEPLKRKGRVGPLLACIIGTQFRKLRDGDRFWWENEGVFSMQQRQALAQ-ISLPR"),
  ("NSNQIKILGNQGSFLTKG-PSKLNDRADSRRSLW--------DQGNFPLIIK------NLKI",
   "NCSTFYVVKEDGTIVYTGTATSMFD-NDTKETVYIADFSSVNEEGTYYLAVPGVGKSVNFKI"),
  ("GKPITVKCSVADV-YPF--DRLEID", "GDSVSLTCSTTGCESPFFSWRTQID"),
  ("PTVTWLRKGQV------LSTSARHQVTTTKYKSTFEIS---SVQASDEGNYSVVVENSEGKQ",
   "PTIAW---GNTKFAIVEVDQAATAYNNLVKVKNAADVSVSWNLWNGDTGTTAKILLN--GKE"),
  ("PNAP----------KLTGITCQADKAEIHWE--QQGDNRSP--ILHYTIQFNTSFTPASW",
   "PDPPIALNWTLLNVSLTGI--HAD-IQVRWEAPRNADIQKGWMVLEYELQYK-EVNETKW"),
  ("QTITELCSEYRNTQIYTINDKILSYTESMAGKREMVIITFKSGETFQVEVPGSQHIDSQKKAIERMKDTLRITYLTETKI----DKLC",
   "QPVIGACTSPYDGKYWSMYSRL----------RKMLYLIYVAGISVRV------HVSKEEQYYDYEDATFETYALTGISICNPGSSLC"),
  ("LTAGHCTN--ISASWSIGT---RTGTSFPNNDYGIIRHS--NPAAADGRVYL-------YNGSYQDI-TTAG----NAFVGQAVQRSGSTTGLRSGSVTGLNATVN-----------YGSSGIVYGMIQTNVCA------QP--GDSGGSLFAGSTAL--GLTS--GGSGNCRTGGTTFYQPVT",
   "LTAAHCVHDAVSVVVYLGSAVQYEGEAVVNSER-IISHSMFNPDTYLNDVALIKIPHVEYTDNIQPIRLPSGEELNNKF--ENIWATVSGWGQSNTDTVILQYTYNLVIDNDRCAQEY-PPGII---VESTICGDTSDGKSPCFGDSGGPFVLSDKNLLIGVVSFVSGAG-CESGKPVGFSRVT"),
  ("ISKFLGFWYEI---AFASKMGTPGLAHKEEKMGAM-----VVELKENLLALTTTYYSEDHCVLEKV-----TATEGDGPAKFQVTRLSGKKEVVVEATD----YLTYAIID---ITSLVAG-AV-HRTMK",
   "VDAFLGTWKLVDSKNFDDYMKSLGVGFATRQVASMTKPTTIIEKNGDILTLKTHSTFKNTEISFKLGVEFDETTADDRKVKSIVT-LDGGKLVHLQKWDGQETTLVRELIDGKLILTLTHGTAVCTRTYE"),
  ("IDVLLGADDGSLAFVPSEFSISPGEKIVF-KNNAGFPHNIVFDEDSIPSGVDASKISMSEEDLLNAKGETFEVALSNKGEYSFYCSPHQGAGMVGKVTV",
   "VKMLNSGPGGMMVFDPALVRLKPGDSIKFLPTDKG--HNVETIKGMAPDGADYVKTTVGQEAV---------VKFDKEGVYGFKCAPHYMMGMVALVVV"),
  ("FISRHNSNFFSDK--LVLTSVTPASSAP-------VLQTPKATSSTLYFDSLTVNAGNG-GFLHCIQMDTSVNAANQVVSVGA-DIAFDA",
   "FMSQYNPITLQQKRFTILKDVTLNCSLTGESIKDRIINLPGQLVN--YNGATAVAASNGPGAIFMLQIGDSL--------VGLWDSSYEA"),
  ("WTKS-EPHSWELIFPIEVCGPNNGFEMWSSEWANQTS------WHLSFLVD---------NPKQST---TFDVLLGISQNFEIAGNTLMP-AF",
   "WTASCAAAEAKVTSAITISLPN---EL-SSERNKQLKVGRVLLW-LGLLPSVSGTVKSCVTETQTTAAASFQVALAVADNSKDVVAAMYPEAF"),
  ("PGKGIPDRFEGKVVTRKDVLNQSINFTANRDTFILIAPTPGVAYWVADVPAGTFPISTTTFNA-VNFPGFNS-----MFGNAAASRSDQVSSFRYASM---NVGIY-PTSNLMQFAGSITVWKCPVKLSNVQFPVATTPATSALVHTLVGLDGVLAVGPDNFS--ESFIK----GV---FSQSVCNEP----DFEFSDILEGIQTLP------PANVTVATSGQPFNLAAGAEAVSGIVGWGNMDTIVIRVSAPTGAVNSAI--LKTWACLEYRPNPNAMLYQFGHDSPPCDEVALQEYRTV----ARSLPVAVIAAQ",
   "PGKG--DE-NSKAMLGQQSMPNRPNYIAFRDNFI------GLMYYNSTGNMGVLAGQASQLNAVVDLQDRNTELSYQLLLDSIGDRT------RYFSMWNQAVDSYDPDVRIIENHGT------EDELPNYCFPLGGIGVTDTYQAIKANGNGSGDNGDTTWTKDETFATRNEIGVGNNFAMEINLNANLWRNFLYSNI---ALYLPDKLKYNPTNVEI--SDNP----------------NTYDYMNKRVVAP-GLVDCYINLGARWS-LDYMDNVN----PFNH-----HRNAGLRYRSMLLGNGRYVPFHIQVPQ"),
  ("PFVYYKCDLEVTLSPHTSGAHGLLVR--------W-----C-----PTGT------PTKPTTQVLH-----EVSSLSEGRTPQVYSAGPGTSNQISFVVPYNSPLSVLPAVWYNGHKRFDNTGDL-----GIAPNSDFG-----TLFFAGTKP",
   "PFEYIRIPLPHVLSGEDGGVFGATLRRHYLVKTGWRVQVQCNASQFHAGSLLVFMAPEYPTLDVFAMDNRWSKDNLPNGTRTQTNRKGPFAMDHQNFW-----QWTLYPHQFLN--LRTNTTVDLEVPYVNIAPTSSWTQHASWTLVIAVVAP"),
  ("DAVINHMCGSGAAAGTGTTCGSY-----CNPG-----SREFPAV------------P-YSAWDFNDGKCKTASGGIESYNDPYQVRDCQLV-----------GLLD----LALEKDYVRSMIADYLNKLIDIGVAGFRIDASKHMWPGDIKAVLDKLHNLNTN-WFPAGSRPFIFQEVID----LGGEAIKSSEYFGNGRVTEFKYGAKLGTVVRKWSGEKMSYLKNWG---EGWGFM------PSDRALVFVDNHDNQRGHGAGGSSILTF--WDARL----------------YKVAVGFMLAHP-------YGFTRVMSSYRWARNFVNGEDVNDWIGPPNNN-----GVIKEVTI",
   "NALRDYAEARGIKIG---TCVNYPFYNNSDPTYNSILQREFSMVVCENEMKFDALQPRQNVFDFSKGDQLLA----FAERNGMQMRGHTLIWHNQNPSWLTNGNWNRDSLLAVMKNHITTVMTHYKGKIVEWDVANECMD--------------DSGNGLRSSIW-----RNVIGQDYLDYAFRYAREADPDALLFYNDYNIE-DLGPKSNAVF-----NMIKSMKERGVPIDGVGFQCHFINGMSPEYLASID--QNIKRYAEIG-VIVSFTEIDIRIPQSENPATAFQVQANNYKELMKICLANPNCNTFVMWGFT---DKYTWIPGTFPGYG-NPLIYDSNYNPKPAYNAIKEALM"),
  ("YAEARGIKIGT-CVNYPFYNNS---D--PTYN-SILQREFSMVVC----------ENEMKFDAL----QPRQNVFDFSKGDQLLAFAERNGMQMRGHTLIWHNQNPSWLTNGNWNRDSLLAVMKNHI-TTVMTHYKGKIVEWDVANE--CMDDSGNGLR--SSI------WRNVIGQDYLDYAFRYAREADPDALLFYNDYNIEDLG-PKSNAVFNMIKSMKERGVPIDGVGFQCHFINGMSPEYLASIDQNIKRYAEIGVIVSFTEIDIRIPQSENPATAFQVQANNYKEL-MKICLANPNCNTFVMWGFTDKYTWIPG",
   "YKQNSGKVVGSYFVEWGVYGRNFTVDKIPAQNLTHLLYGF-IPICGGNGINDSLKEIEGSFQALQRSCQGRED-FKISIHDPFAALQK----AQKGVT-AW--DDP---YKGNFGQ--LMALKQAHPDLKILPSIGG----WTLSDPFFFMGDKVKRDRFVGSVKEFLQTWKFFDGVD-IDWEFPGGKGANP------------NLGSPQDGETYVLL--MKELRAMLDQLSTE----TGRKYELTSAISAGKDKIDKVAYNVAQNSMDHIFLMSYDFYGAF-----DLKNLGHQTALNAP------AWKPDTAYTTVNG"),
  ("KHYASAISVLTDEKYFQGSFNFLPIVSQIAP---QPILCKDFIID--PYQIYLARYYQADACLLMLSVLDDDQYRQLAA---VAHSLEMGVLTEVSNEEEQERAIALG----AKVVGINNRDLRDLSIDLNRTRELAPKLGHNVTVISESGINTYAQVRE--LSH-FANGFLIGSALMAHDDLHAAVRRV",
   "ERYENLFAQLNDRR--EGA--FVPFVTLGDPGIEQSLKIIDTLIDAGADALELGVPFSADG-----PTIQNANLRAFAAGVTPAQCFEMLALI-----REKHPTIPIGLLMYANLV-FNN------GIDAFYAR--CEQVGVDSVLVADVPVEESAPFRQAALRHNIAPIFICPPN--ADDDL---LRQV"),
  ("EITLALLRN-AISDNVKANK--------------HK-------F-------LIDGFPRKMDQAISFERDI-VESKFILFFDC-----PEDIMLERLLERGKT-----SGRSDDN-------IESIKKRFNTFKETSMPVIEYFETKSKVVRVRCDRSVEDVYKDVQDAIR",
   "EYDLKRLRNIGIAAHIDAGKTTTTERILYYTGRIHKTAAVTTCFWKDHRINIID-TPGHVDFTIEVERSMRVLDGAIVVFDSSQGVEPQSETVWRQAEKYKVPRIAFANKMDKTGADLWLVIRTMQERLG-----ARPV---------VMQLPIGR--EDTFSGIIDVLR"),
  ("KGLF--GKKEMRILMVGLDAAGKTTILYKLKL--------------GEIVTTIPT-----IGFNVETVEYKNIS--FTVWDVGGQ-DKIRPLWRHYFQNTQGLIFVVDSNDRERVNEAREELMRMLAEDELRDAVLLVFANKQDL---PNAMNAAEI",
   "KGEFIRTKPHVNVGTIGHVDHGKTTLTAALTFVTAAENPNVEVKDYGDI-DKAPEERARGITINTAHVEYETAKRHYSHVDCPGHADYIKNMITGAAQ-MDGAILVVSAADGP-MPQTREHI--LLAR-QVGVPYIVVFMNKVDMVDDPELLDLVEM"),
  ("PMILGYWNVRGLTHPIRLLLEYTDSSYEEKRYAMGDAPDYDRSQWLNEKFKLGLDFPNLPYLIDGSRKITQSNAIMRYLAR",
   "PYTIVYFPVRGRCEAMRMLLADQGQSWKEEVVTI--------DTWMQGLLKPTCLYGQLPKFEDGDLTLYQSNAILRHLGR"),
  ("GWGQAKPVERA---FGVGGLEEYVDKHPA-LAGAYTG---------P--------YGDRAVDAELFLKTLAEGVAM------------FNARTGRETEMCGGKLSFDDVFEREYARTIVRKP",
   "GYGVIRQVGRQLSYLGSGCIRTKVDDLPSRLKLIYAGVTEIITQFQPDYFAIEQVFMAKNADSALKLGQ-ARGVAIVAAVNQELPVFEYAARQVKQTVVGIGSAEKSQV--QHMVRTLLKLP"),
  ("DTYAATRYPVILVHGLAGTDK-----FANVVDYWYGIQS------DLQSHGAKVYVANLSGFQSDDGPNGR--GEQLLAYVKQVLAAT---GATKVNLIGHSQG----GLTSRYVAAVAPQLV----ASVTTIGTPH----RGSEFADFVQDVLKTDPTGLSSTVIAAFVNVFGTLVSSSHNTDQDALAALRTLTTAQTATYNRNFPSAGLGAPGSCQTGAATETVGGSQHLLYSWGGT",
   "NTGRKTRF---IIHGF--IDKGEESWLSTMCQNMFKVESVNCICVDWKSGSRTAY--------SQASQNVRIVGAE-VAYLVGVLQSSFDYSPSNVHIIGHSLGSHAAGEAGRRTNGAVGRITGLDPAEPCFQGTPELVRLDPSD-AQFV-DVIHTD--------IAPFIPNLGFGMS-------------------QTAGHLDFFPNGGKEMPG-CQKNVLSQIV----DIDGIWQGT"),
  ("KDKLMASHPGLVVELVPMVTRGDVIGKGLFVKELEVALLENRADIAVHSMKDVPVEFPQGLGLVTICEREDPRDAFVSNNYDSLDALP",
   "KDKFEIVTPSESILAEPTVS---VVDK----------VVEKKDTKAVAEAYLKYLYSPEG-QEIAAKNFYRPRDADVAKKYD--DAFP"),
  ("AAGCVRHTVEDAIAKGFRPIIPRETIGDRVPGVVQWNLYDIDNKFGDVESTDSVVQYLDALPQFEDTVPKTLSDPQPEVEAPADP",
   "AAGQIGYSLLFRIAAG-------EMLGKDQPVILQ--LLEIPQA---MKALEGVVM------ELEDCAFPLLAG----LEATDDP"),
  ("AGPSG---SEAARVLMESGYTVHLTDTAEKIGGHLNQVAALPGLGEWSYHRDYRETQITKLLKKNKESQLALGQKPMTADDVLQYGAD",
   "AGGVGLIACQWAKAL--GAKLIGTVGTAQKAQSALKA-------GAWQV-INYREEDLVERLK-----EITGGKKVRVVYD--SVGRD"),
  ("GCKAAGAARIIGVDINKDKFAKAKEVGATECVN----PQDYKKPIQEVLTEMSNGGVDFSFE-VIGRL--DTMV---TAL-----SCCQEAYG---VSVIVGVP--PDSQNLSMN",
   "GIASINRAAMYG---QKCALIEAKELGGT-CVNVGCVPKKVMWHAAQIREAIHMYGPDYGFDTTINKFNWETLIASRTAYIDRIHTSYENVLGKNNVDVIKGFARFVDAKTLEVN"),
  ("GATLGAAFAAEEIDPKKRVILFIGDGSLQLTVQEISTMIRWGLKPYLFVLNND-------------GYTIEKLIHGPKAQYNEIQGWD-HLSLLPTFGAKDYETHR--VATTGE--WDKLTQDKSFNDNSKIRMI",
   "GGGTGGATAAKYIK--------LADPSIEVTLIEPNT------DYYTCYLSNEVIGGDRKLESIKHGYDGLR-AHGIQVVHDSATGIDPDKKLVKTAGGAEFGYDRCVVAPGIELIYDKIEXQRA----GKIAQI"),
  ("GTVGT-AQKAQSALKAGAWQVI-NYREEDLVERLKEITGGKKVRVVY-DSVGRDTWERSLDCLQRRGLMVSFGNSSGAV--TGVNLGIL",
   "GARGLGAEAARQAVAAGARVVLADVLDEEGAATAREL--GDAARYQHLDVTIEEDWQRVV-AYARE----EFGSVDGLVNNAGISTGMF"),
  ("DGDRVPALVIGSGYGGAV---AALRLTQAG------IPTQIVEMG---RSWDTPGSDGKIFCGMLNPDKRSMWLADKTDQP",
   "DG----VVIVNCSRGRLVDTDAVIRGLDSGKIFGFVMDTYEDEVGVFNKDW-----EGKEF-----PDKR---LADLIDRP"),
  ("IIAGG-GLTG--LTTAARLTENPNISVLVIESGSYESDRGPIIEDLNAYGDIFGSSVDHAYETVELATNNQTALIRSGNGLGGSTL-----VNGGTWTRPHKAQVDSWETVFGNEG--WN---WDNVAAYSLQAERARAPNAKQIAAGHYFNASCHG---------VNGTVH------AGPRDTGDDYSPIVKALMSAVED-----RGVPTKKDFGCGDPHGVSMFPNTLHEDQVRSDAAREWLLPNYQRPN----LQVLTGQYVGKVLLS",
   "LVYGGRGALGSRCVQAFRARNWWVASIDVVEN---EEASASVIVKMT---DSFTEQADQV--TAEV-----------GKLLGDQKVDAILCVAGG-WAGGNAKS----KSLFKNCDLMWKQSIWTSTISSHLATKHLKEGGLLTLAGA---KAALDGTPGMIGYGMAKGAVHQLCQSLAG-KNSGMPSGAAAIAVLPVTLDTPMNRKSMP-EADFSSWTP--LEFLVETFH----------DWITGN-KRPNSGSLIQVVTTD--GKTELT"),
  ("VIITGGARGLGAEAARQAVAAGARVVLA-----DVLDEEGAA----TARELGDAARYQHLDVTIEEDWQRVVAYAREEFGSVDG-LVNNAGISTGMFLETESV",
   "VVVGGGYIGLELGIAYRKLGAQVSVVEARERILPTYDSELTAPVAESLKKLGIA---LHLGHSVE--------------GYENGCLLANDGKGGQLRLEADRV"),
  ("VTG-GSGYIGSHTCVQLLQNGHDVIILDNLCNSKR-SVLPVIERLGGKHPTFVEGDIRN-----EALMTEILHDHAIDTVIHFAGLKAVG-----ESVQKP-----LEYYDNNVNGTLRLISAMRAANVKNFIFSSSATVYGDQPKIPYVESFPTGTPQSPYGKSKLMVEQILTDLQKAQPDW--SIALLRYFNPVGAHP-SGDMGEDPQGIPNNLMPYIAQVAVGRRDSLAIFGND--YPTEDGTGVRDYIHVMDLA--DGHVVAMEKLANKP",
   "VIGAGSG--G-------LEAGWNAASL----HKKRVAVIDLQKHHGPPHYAALGGTCVNVGCVPKKLM--VTGANYMDTIRESAGF---GWELDRESV-RPNWKALIAAKNKAVSGINDSYEGM-FADTEGLTFHQGFGALQDNHTVLVRES---ADPNSAVLET-LDTEYILL----ATGSWPQHLGIXRV-------PRSQTLQLEKAGVE------VAKNGAIKVDAYSKTNVDNIYAIGD---VTDRVMLTPVAINEGAAFVDTVFANKP"),
  ("DCIKDFTD---DQAQAEAFIEHFSY-RAHD---VTDAASY-AVLKEAIEEAADKFDIDGNRIFYMSVAPRFF--GTIAKYL----KSEGLL",
   "DAVKDDFDVFIDFTRPEGTLNHLAFCRQHGKGMVIGTTGFDEAGKQAIRDAAADIAI----VFAANFSXMTFANGAVRSALWLSGKESGLF"),
  ("VGINGFGRIGRNVFRAALKNPDIEVVAVNDLTDA--NTLAHLLKYDSVHG--RLDAEVS-VNGNNLVVNGKEIIVKAERDPENLAWG-----EIGVDIVVESTGRFTKREDAAKHLEAGAKKVII---SAPAKNED-----ITIVMGVNQDKYDPKAHHVISN-ASXNETG",
   "IALIGLAVMGQNLI-LNMNDHGFVVCAFNRTVSKVDDFLANEAKGTKVLGAHSLEEMVSKLKKPRRII----LLVKAGQAVDNFIEKLVPLLDIG-DIIID--GGNSEYRDTMRRCRDLKDKGILFVGSGVSGGEDGARYGPSLMPGGNKEAW-PHIKAIFQGIAAKVGTG"),
  ("VAAGRIGLAVLRRLAPFDVHLHYTDRHRLPESVEKELNLTWHATREDMYPVCDVVTLNCPLHPETEHMI",
   "VAATGATASVIYVPAPF-------CKDSILEAIDAGIKLI--ITITEGIPTLDMLTVKVKLDEAGVRMI"),
  ("RRAG-EGEKMIRTRSWPGWEPLELVGEKLDNKTL-GIYG--FGSIGQALAKRAQGFDMDIDYFD",
   "RRDGFRAEKILIKRLMDKVENGNIILH--TNRTLEEVTGDQMGVTGVRL-RDTQNSD-NIESLD"),
  ("GGY-------ISIEFAGIFNAYKARGGQ--VDLAYRGDM--ILRGFDSELRKQLTEQLRANGINVRTHENPAKVTKNADGTRH",
   "GGYEDTAGSDVVVITAGI----PRQPGQTRIDLA--GDNAPIMEDIQSSL-----DEHNDDYISLTT-SNPVDLL-----NRH"),
  ("KVAVIGGGNTAVEEAL", "KVLVV--GNPANTNAL"),
  ("IIGGGPGGYVAA--IRAGQLGIPTVLVEGQALGGTCLNIGCIPSKALIHVAEQFHQASRFTEPSPLGISVASPRLDIGQSVAWKDGIVDRLTTGVAALLKKHGVKVVHGWAKVL--DGKQVE-VDGQRIQCEHLLLATGSSSVELPXRRPRTK--GFNLECLDLKMNGAAIAID-ERCQTSMH-NVWAIGDVA-GEPM--LAHRAMAQGEMVAEII",
   "VVGGGTGGATAAKYIKLADPSIEVTLIEPNTDYYTC------------YLSNEVIGGDRKLESIKHGY----------------DG------------LRAHGIQVVHDSATGIDPDKKLVKTAGGAEFGYDRCVVAPG---IELIYDKIEXQRAGKIAQIAGLTNDAGWCPVDIKTFESSIHKGIHVIGDASIANPMPKSGYSANSQGKVAAAAV"),
  ("KVAVLGAAGGIGQALALLLKTQLPSGSELSLYDIA--------P----VTPGVAVDLSH---IPTAVKIKGFSGEDATPALEGADVVLISAGVRRKPGMDRSDLFNVNAGIVKNLVQQVAKTCPKACIGIIT-NPVN",
   "KVSVVGAAGTVGAA----------AGYNIALRDIADEVVFVDIPDKEDDTVGQAADTNHGIAYDSNTRVRQGGYEDTA----GSDVVVITAGIPRQPGQTRIDLAGDNAPIMED-IQSSLDEHNDDYISLTTSNPVD"),
  ("LDNYRGY--SLGNWVCA--------AKFESNFNTQATNRNTDGSTDYGILQINSRWWCNDGRTPGSRN-LCNIPCSALLSSDITASVNCAKKIVS---D---GNGMNAWVAWRNRCKGTDVQAWIR",
   "MDRYKTIIKKVGEKLCVEPAVIAGIISRESHAGKVLKNGWGDRGNGFGLMQVDKR----SHKPQGTWNGEVHITQGTTILTDFIKRIQ--KKFPSWTKDQQLKGGISAY----NAGAG-NVRSYAR"),
  ("FMARRGFL----PM--TLPSYAREKAFLGTGHFPAYRDQVWAIAETDLYLTGTAEVVLNALHSGEILPYEALPLRYAGYAPAFRSEAGSFGKDVRGLMRVHQFHKVEQYVLTEASLEASDRAFQELLENAEEILRLL-ELPYRLVEVATG----DMG-PGKWRQVDIEVYLPSEGRYRETHSCSALLDWQARRANLRYRDPEG-RVRYAYTLNNTALATPRILAMLLENHQLQDGRV-RVPQALIPYMGKEVLEP",
   "FMVARGFMEVETPMMQVIPGGASARPFI-THH---------NALDLDMYLRIAPELYLKRLVVG---GFE----RVFEINRNFRNE----GISVR--------HNPE-FTMMELYMAYAD--YHDLIELTESLFRTLAQEVLGTTKVTYGEHVFDFGKP--FEKLTMREAIK---KYRPETDMADLDNFDAAKA---LAESIGITVEKSWGLGR---IVTEIFDEVAEAHLIQPTFITEYPAEVSPLARRNDVNP"),
  ("DSARFRYVLSEEFDAPVQNVEGTILGEHGDAQVPVFSKVRVDGTDPEFSGDEKEQLLGDLQESAMDVIERKGATEWGPARGVAHMVEAILHDTGEVLPASVKLEGEFGHEDTAFGVPVRLGSNGVEEIVEWDLDDYEQDLMADAAEKL",
   "DTARFRFLLGEYFSVAPQNVHAYIIGEHGDTELPVWSQAYI-GVMP----IDLERIFVNVRDAAYQIIEKKGATYYGIAMGLARVTRAILHNENAILTVSAYLDGLYGERDVYIGVPAVINRNGIREVIEIELNDDEKNRFHHSAATL"),
  ("CVNGGECFMVKDLSNPSRYLCKCQPGFTGARC", "CLNQGHC---KD--GIGDYTCTCAEGFEGKNC"),
  ("ICCTKCHK--GTYLYNDCPGPGQDTDCRECESGSFTA--SENHLRHC",
   "LSCSKCRKEMGQVEISSCT-VDRDTVC-GCRKNQYRHYWSEN-LFQC")]
 
true_neg_align = [("PADEEMLFIYS", "PQDRESLFYFN"),
  ("CTGKHFLNEQQLMQASQYAGYAEHKKAHDD---FIHKLDTW-DGDV-TYAKNWLVNHIKTIDF-KYRG",
   "CHGTH----------GNSVGPASPSIAQMDPMVFVEVMEGFKSGEIASTIMGRIAKGYSTADFEKMAG"),
  ("NMC", "NKC"),
  ("KKIHENEKRLE----AGDHPVELLARDFEKNYNMYIFPVHWQFGQLDQHPID",
   "KKAGEGAKVIELQGIAGTSA----ARERGEGFQQAV-AAH-KFNVLASQPAD"),
  ("CLKGIARGIHNLNEDNAR--SIPPKCGVNLPYTISLNIDCSRV",
   "CYDGAYRNDDETCEERAARIKIGPKC-------VKAFKDCCYI"),
  ("DPDEVARRWGERKSKPNMNYDKL", "DADMVVITAGPRQ-KPGQSRLEL"),
  ("NGIEYGDMQLICEAYHLMK---DVLGLGHKEMAKA-FEEWNKTELDSF---LIEITASILKFQDADGKHLLPKIRDSAGQKGTGKWTAISALEYGVPVTLI-----GEAVFARCLS",
   "NDLSRDDLNLVLATAAKLKANPQPELLKHKVIASCFFEASTRTRL-SFETSMHRLGASVVGFSDSANTSL--------GKKGETLADTISVISTYVDAIVMRHPQEGAARLATEFS"),
  ("RMKQLEDKVEELLSKNYHLENE", "RLTEIWGNGNEETSEVFPLKTK"),
  ("GKKGDTVELTCTASQKKSIQFHWKNSNQIKILGNQGSF-----LTKGPSKLNDRADSRRSLWDQGNF",
   "GKK-------CQSWSSMTPHRHQKTPENYPNAGLTMNYCRNPDADKGPWCFTTDPSVR---WEYCNL"),
  ("DAPSQIEVKDV", "DACEQAAIQCV"),
  ("DTFW-H--------TFYGANGDPKPPH---------------------------TYT-----IDM-----KTTQNVNGLSMLPRQDG-NQNGW--IGRH----EVYLSSDGTN-WGSPVASGSWFA---DSTT--KYSNFE",
   "DTFFMHCRDGSLAKVAFGTEGTPEITHTEVFHPEDEFLINHPAYSQKAGRLVWPTYTGKIHQIDLSSGDAKFLPAVEALTEAERADGWRPGGWQQVAYHRALDRIYLLVDQRDEWRHKTA--SRFVVVLDAKTGERLAKFE"),
  ("IIG", "IIG"),
  ("ALKSRIALTVEDS", "ALENQHTIDLTDS"),
  ("GTYYLLPQVWAQGGGV----QLAKTGEETCPLTVVQSPNELSDGKPIRIESRLRSAFIPDDDKVRIGF---AYAPKCAPSPWWTVVEGLSVKLSEDESTQFDYPF--KFEQVSDQLHSYKLLYCE----GKH-EKCASIGI",
   "GSQDLANQYKSEGKNVVSALQLDMTNYKGSAQDVV-FITDYTDSNFTQYLTQLMDEYLP---SLTYGFDTCGYA--CSDHASWHNA-GYPAAM----------PFESKFNDYNPRIHTTQDTLANSDPTGSHAKKFTQLGL"),
  ("SEVDGQTIYTPSKSTTAKLLSGATWSISYGDGSSSSGDV",
   "TQIHG--LYRSSDKT------GGYWKITMNDGSTYQSDL"),
  ("KCSLTGKWTNDLGSNMTIGAVNSRGEFTGTYTTAVTATSNEIKESPLHGTENTINKRTQPTFGFTVNWKFSESTTVFTGQCFID-RN-GKEVLKTMWLLRSSVND",
   "RSSLPGFYR----TSLTLAAPEAAGE-----VERLIGHPLPLRLDAITGPEEE-GGRLETILG----WPLAERTVVIPSAIPTDPRNVGGD------LDPSSIPD"),
  ("GSLAFVPSEFSISPGEKIVFKNNAGFPHNIVFDEDSIPSGVDAS",
   "GSSKYPNCAYKTTQANKHIIVACEGNPY--------VPVHFDAS"),
  ("MTPIVFTNLAEGETVSIKASGSVNR", "MKPVTLYDVAEYAGVSYQ---TVSR"),
  ("GDALVVL", "GDAEIVL"),
  ("LCTRSGNENEFRDMVTRCNNVGVRIYVDAVINHMCGSGAAAGTGTTCGSYCNPGSREFPAVPYSAWDFNDG-KCKTASGG--IESYNDP-YQVRDCQLVGLLDLALEKDYVRSMIADYLNKLIDIGVAG-FRIDASKHMWPGDIKAV",
   "LIVRIAAANGIGRLVIGQNGILSTPAVSCIIRKIKAIGGIILTASH-----NPGG------PNG--DF--GIKFNISNGGPAPEAITDKIFQISKT----IEEYAICPDL-----------KVDLGVLGKQQFDLENKFKPFTVEIV"),
  ("GGPAQFVNATSKEVVEWAAKLGLP---LVFRWDDSNAQRKEYAGLYHEVAQAHGVDVSQVRHKLTLLVNQNVDGEAAR--AEARVYLEEFVRESYSNTDFEQKMGELLSENAIGTYEESTQAARVAIECCGA",
   "GKTCAFIDAEHALDPIYARKLGVDIDNLLCSQPDTGEQALE---ICDALARSGAVDVIVVDSVAALTPKAEIEGLAARMMSQAMRKLAGNLKQSNTLLIFINQTG----GNALKFY----ASVRLDIRRIGA"),
  ("CGIC-----GEHGGDPSSVEFCHKV------GLNYVSC",
   "CAVCNDYASGYHYG-VWSCEGCKAFFKRSIQGHNDYMC"),
  ("VQNEPDEAEQDC", "VQGGPGPSGQCC"),
  ("W--TPAKIIDDAVA--RIREQVGDDKVILGLSGGVDSSVTAMLLHRAIGK--NLTCVFVDNGLLRLNEAEQVLDMFGDHFGLNIV------HV----PAEDRFLSALA--GENDPEAKRKIIGRVFVEV-FDEEALKLEDVKWLAQGTIYPDVIESAAKMGLVEPLKEL",
   "WQTQPA-LLTASVALYRVWQQ----------QGG---KAPAMMAGHSLGEYSALVCA----GVIDFADAVRLVEMRG-KFMQEAVPEXVPSHCALMKPAADKLAVELAKITFNAPTV--PVVNNVDVKCETNGDAIRDALVRQLYNPVQWTKSVEYMAAQG-VEHLYEV"),
  ("GGIGKSTTTQNLVAALAEMGKKVMIVGCDPKADSTRLILHSKAQ",
   "GQLGIPTV---LVEGQA-LGGTCLNIGCIP----SKALIHVAEQ"),
  ("RWLLAAAGVEFEEKF", "RWINNDVLI-YERKF"),
  ("RPYALVVTEDDVISISANIDGGQPWRRTV--------GTD--NIVYTDWQRDNYFAAIQQALPKARRIGIEHDHLNLQNRDKLAARYPDAELVDVAAACM",
   "RIHYIQVRHEEVGAMAAAADAKLTGKIGVCFGSAGPGGTHLMNGLY-DAREDH--VPVLALIGQFGTTGMNMDTFQEMNENPIYADVADYNVTAVNAATL"),
  ("VYKD-ETGKTPVLTS----VKKAEQYLLENETTKNYLGIDGIPEFGRCTQELLFGKGSALINDKRARTAQTPGGTGALRVAADFLAKNTSVKRV--WVSNPSWPNHKSVFNSAG",
   "IYTNYENGKNDYVKALPGHLKPFETLLSQNQGGKAFIVGDQISFADYNLLDLL------LIHQVLA-----PGCLDNFPLLSAYVARLSARPKIKAFLSSPEHVNRP--INGNG"),
  ("PTSQRLQLLEPF---DKWDGKDLEDLQILIKVK--GKCTTDHISA------AG-P------WLKF--RGHLDNI-----SNNLLIGA---INSENRK-ANSVRNAVTQEFGPVPDTARYYKQHGIRWVVIGDENYGEGSSREHSAL-EP-RFLGGRAIITKSFARIHETNLKKQGLLPLTFADPADYNKIHPVDKLTIQGLKDFAPGKPLTCIIKHPNGTQETILL",
   "PTADR-----PFVLGLPTGGTPMTTYKALVEMHKAGQVSFKHVVTFNMDEYVGLPKEHPESYYSFMHRNFFDHVDIPAENINLLNGNAPDIDAECRQYEEKIRS--------------YGKIH----LFMG----GVGND-GHIAFNEPASSLASRTRI-KTLT--HDTRVANS-----RFFD-NDVNQV-PKYALTV-GVGTLLDAEEVMILVL---GSQKALAL"),
  ("INGEQI--RSISDLHQTLKKELA--------LPEYYGE-----NLDALWDALTGWVE",
   "VTGEALVQRS-DDNADALKKRLAAYHAQTEPIVDFYKKTGIWAGVDASQPPATVWAD"),
  ("WVA-------SDFDALIPSLKAKKIDAIISSLSITDKRQQEIAFSDKLY-AADSRLIAAKGSPIQPTLESL----KGKHVGVLQGSTQEAYANDNWRTKGVDVVAYANQ---------DLIYSDLTAGRLDAALQDEVAASEGFLK---QPAGK----------EYAFA------GPSV--KDKK------YFGDGTGVGLR-KDDTELKAAF--DKALTELRQDGTYDKMAK-KYFDFNVYGD",
   "WLAMEKPESFTDKD-----IRAAK-NAAFNALKIYDEAEVNKYFVRAQYGAGDS-------ADFKPYLEELDVPADSKNNTFIAGELQ--FDLPRW--EGVPFYVRSGKRLAAKQTRVDIVFK---AGTFNFGSEQE--AQEAVLSIIIDPKGAIELKLNAKSVEDAFNTRTIDLGWTVSDEDKKNTPXGSNFADWNGVSIAWKFVDAISAVYTADKAPLETYKSGSMGPEASDKLLAAN--GD"),
  ("AALKNPDIEVVAVNDLTDANTLAHLLKYDSVHGRLDAEVSVNGNNLVVNGKEII-VKAERDPENLAWGE-IGVDIVVESTGRFTKREDAAKHLEAGAKKVIISAPAKNEDITIVMGVNQDKYDPKAHHVISNASXNETGY-SHRVVDLAA",
   "AHLAAENIDAAIFTSYHNINYYSDFL-YCS-FGRPYA--------LVVTEDDVISISANIDG-GQPWRRTVGTDNIVYTDWQRDNYFAAIQQALPKARRI---------------GIEHD------HLNLQNRDKLAARYPDAELVDVAA"),
  ("TKVSVVGAAGTVGAAAGYNIALRDIADEVVFVDIP",
   "SNIFLTGAFGTI--ASGPFXXXXXXXXXXXXLDIP"),
  ("VREIGYVHSDGFDANSCAVLSAIGKQ", "VRDGAYGYQNLFDTTVDAFYTAMGKH"),
  ("TDAKKVITTPLTVGYF---DGKDGLKQDARLLKVISYLDVGDGNYW",
   "SDKPK---RPLSAYMLWLNSARESIKRENPGIKVTEVAKRG-GELW"),
  ("FGKLHVLSTPNQDN-PVMEGFTPIVGIDVW--EHAYYLKYQNR",
   "FKQQKYLSAPEREHLASMIHLTP-TQVKIWFQNHRYKMKRQAK"),
  ("CESCVEIAPGAFAMDPEIEKAYVKDVEGASQEEVEEAMDTCPVQC",
   "C-SCTK------SMPP---KCRCSDIRN----------DFCYEPC"),
  ("AAHFAAESSTGTR--------GVDALVYEVDEARELTKIAY",
   "AAILSNEKFRGVRLYVSENQLKITANNPEQEEAEEILDVTY"),
  ("APSNFANGVAEWISSNSRSQAYKVTCSV--RQSSAQNRKYTIKVEVPKV---ATQTVGGV",
   "APSGTALAMGEAI-AHALDKDLK-DCAVYSREGHTGER-------VPGTIGFATVRAGDI"),
  ("PETGTLRILQEEEWRGLGIT-QSLG", "PDGMQIKITRQEIGQIVGCSRETVG"),
  ("DPEGNK", "DPICNK"),
  ("VCANKGVTAERP----KNDADMKPGQGYG", "ICRRRSAGFKGPCMSNKNCAQVCQQEGWG"),
  ("VYYFGQEGLHNVLVIDLLGPSLEDLLDLCG-RKFSVKTVAMAAKQMLARVQSIHEKSLVYRDIKPDNFLIGRPNSKNANMIYV-VDFGMVK--FYRDPVTKQHIPYREKKNLSGTARYMSINTHLGREQSRRDDLEA",
   "LMYKGQPMTFRLLLVDT--P--ETKHPKKGVEKYGPEASAF-TKKMVENAKKIEVE--FDKGQRTDKY--GR------GLAYIYADGKMVNEALVRQGLAK--VAYVYKPN----------NTH---EQHLRGKSEA"),
  ("QTIATKNDDGTYTLNGSKIFITNGGAADI--YIVFAMTDK",
   "QRWVFKNDGSIYSLYDDMVMDVKGSDPSLKQIILWPYTGK"),
  ("GGNGDMTYARLGFKGETQINSDLTGYGQWEYNFQGNNSEGADAQTGNKTRLAFAGLKYADVGSFDYGRNYGVVYDALGYTDMLPEFGGDTAYSDDFFVGRVGGVATYRNSNFFGLVDGL------------NFAVQYL--GKNERDTARRSNGDGVGGSISYEYEGFGIVGAYGAADRTNLQEAQPLGNGKKAEQWATGLKYDANNIYLAANYGETRNATPITNKFTNTSGFANKTQDVLLVAQYQFDF---GLRPSIAYTKSKAKDVEGIGDVDLV-NYF",
   "GGRGAL--------GSRCVQA-FRARNWWVASIDVVENEEASASVIVKMTDSFT--EQADQVTAEVGKLLG---DQ--KVDAILCVAG----------GWAGGNA--KSKSLFKNCDLMWKQSIWTSTISSHLATKHLKEGGLLTLAGAKAALDGTPGMI-----GYGM--AKGAVH----QLCQSLA-GKN-----SGMPSGAAAI---AVLPVTLD-TPMNRKSMPEADFSSWTPLEFLVETFH-DWITGNKRPN---SGSLIQVVTTDGKTELTPAYF"),
  ("ARIENKQDIQIVK", "AKLSNAQVIDVTK"),
  ("MDSVCP---QG----KYIHP", "IDASSPFSQKGDERYKYVDP"),
  ("CKG", "CYG"),
  ("QACDICRLKKLKCSKEKPKCAKCLKNNWECRYSP", "KSVSLC-VKRLIYTNDAGETIKGVCSNFLCDLKP"),
  ("DQRWQSVLARDPNADGEFVFAVRTTGI--FCRPSCR-ARHALR---EN--VSFYANAS-EALAAGF",
   "NQAW------DEQVASVYIPELEALGVELIEQPVGRENTQALRRLSDNNRVAIMADESLSTLASAF")]
#Utilize running output to verify that actual function is running correctly.
#Confirm data stays the right length:

def test_geneticAlg():
	
	pop, scores, library, objectiveMeans = optimizeMatrix_geneticAlg(
		blosum50, 1, 1, 1, 10, 5, 5, 3, -10, -1, true_pos_align, true_neg_align)
	return None
	#assert len(pop) == 10
	#assert len(scores) == 10
	#assert len(library) == 3
	#assert len(objectiveMeans) == 5 + 1 # there's an extra -inf at the beginning so that it starts off increasing

#######################################################################################################
#If time, will test objective function:
#true_pos_aligns, true_neg_aligns = getAlignments(scorematrix=blosum50, gap_start=-10, gap_extend=-1)
#def test_objective_function():
	#recall Obj. fxn: def scoreMatrix(true_pos_align, true_neg_align, scorematrix, gap_start, gap_extend):
	#true_pos_align = true_pos_aligns[7, ?, ?]
	#true_neg_align = true_neg_aligns[?, ?, ?]
	#assert scoreMatrix(true_pos_align, true_neg_align, blosum50, -10, -1) == XXX

#######################################################################################################	
#Also important to check that the matrix changes but will not make excessive # of copies
#def test_mutate_matrix():
	#return None
    #M = blosum.50.copy()
    #M2 = mutateMatrix(M, 1, 1)
    #assert M2 = M
    #assert M2 != blosum50
    

#def test_selection():
	#return None
	# Very random and hard to test, but I want to make sure that my weights
	# are working and that I'm sampling with replacement: 
	#L = ["a", "b", "c"]
	#w = [0, 0, 1]
	#assert selection(L, w) == ["c", "c", "c"]

#def test_scale_scores():
	#return None
	
	# just make sure I'm doing what I said: largest/smallest is 10^selectivePressure
	# and the sum is 1.
	#a,b = scaleScores([1, 2], 1)
	#assert a+b == 1
	#assert b/a == 10
	


		



















