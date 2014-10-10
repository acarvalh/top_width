/////////////////////////////////
// cuts
double const  weight =1.;///10000;//0.001;//
double Mjj =0;//400; 
double PTjj = 0;//400; 
double DeltayVBF = 0;//3;
double DeltaRVBF = 0;//3;
bool shower=true;
// To be applied only to hadron level events
double const jet_ptmin=30.0; // for jet reconstruction
double const rapmax=5.0; // for jet reconstruction
double const jet_ptminvbf=0.0; // We cut on all jets below 50 GeV
double const jet1_ptminvbf=0.0; // We cut on all jets below 50 GeV
double const jet2_ptminvbf=0.0; // We cut on all jets below 50 Ge
double const bjetpt = 30;
double const etab = 2.5;
double const etaj=4.7;
double const higgs_mass = 125.0;
int const cat =1; // minimum number of btag
double const RR =0.4;
////////////////////////////////////////
// weights b-tag
double const subjet2b=1;
double const fatjet2b=1;
double const normalb=1;
double const normalc=1;
double const normall=1;
double const misb=1;
/////////////////////////////////////////
double const H1_ptmin=0.0; // We cut on all jets below 50 GeV
double const H2_ptmin=0.0; // We cut on all jets below 50 GeV
double const HH_ptmin=0.0; // We cut on all jets below 50 GeV
double const MHH = 0; // minimum
///////////////////////////////////
// for substructure
// mass drop
double const Rsb = 1.1; // CA recluster
double const mu = 0.67;
double const ycut = 0.09;
double const Mfat =100;
//
double const Rfilt = 0.1;
int const n_subjet =3;
///////////////////////////////////
// on the 4b's analysis
double const minMH=0;
double const toleranceH2=150;
double const toleranceH1=150;
double const toleranceX=200;
double const HThiggses = 0;
double const DetaHH = 20;//1.3;
double const DetaH = 150;//1.5;
//////////////////////////////////
// on the wwbb analysis
double const wmass = 80.4;
double const bmass = 4.7;
double const tmass = 173.0;
double const ptlepton = 20.0;
double const lepiso = 0.3;
double const MeeMax = 3000.0;
double const MetMin = 20.0;

/*
PYTHIA Abort from Pythia::next: reached end of Les Houches Events File  

 *-------  PYTHIA Event and Cross Section Statistics  -------------------------------------------------------------*
 |                                                                                                                 |
 | Subprocess                                    Code |            Number of events       |      sigma +- delta    |
 |                                                    |       Tried   Selected   Accepted |     (estimated) (mb)   |
 |                                                    |                                   |                        |
 |-----------------------------------------------------------------------------------------------------------------|
 |                                                    |                                   |                        |
 | Les Houches User Process(es)                  9999 |      100000     100000     100000 |   4.481e-08  0.000e+00 |
 |    ... whereof user classification code          0 |                            100000 |                        | 
 |                                                    |                                   |                        |
 | sum                                                |      100000     100000     100000 |   4.481e-08  0.000e+00 |
 |                                                                                                                 |
 *-------  End PYTHIA Event and Cross Section Statistics ----------------------------------------------------------*

 *-------  PYTHIA Error and Warning Messages Statistics  ----------------------------------------------------------* 
 |                                                                                                                 | 
 |  times   message                                                                                                | 
 |                                                                                                                 | 
 |      1   Abort from Pythia::next: reached end of Les Houches Events File                                        | 
 |      1   Error in Pythia::check: energy-momentum not conserved                                                  | 
 |      1   Error in Pythia::next: check of event revealed problems                                                | 
 |    944   Error in Pythia::next: hadronLevel failed; try again                                                   | 
 |     57   Error in SpaceShower::pT2nearQCDthreshold: stuck in loop                                               | 
 |    772   Error in StringFragmentation::fragment: stuck in joining                                               | 
 |    172   Error in StringFragmentation::fragmentToJunction: caught in junction flavour loop                      | 
 |    365   Warning in MultipleInteractions::pTnext: weight above unity                                            | 
 |      3   Warning in ParticleDataEntry::initBWmass: switching off width                                          | 
 |      8   Warning in Pythia::check: energy-momentum not quite conserved                                          | 
 |      1   Warning in Pythia::initSLHA: No MODSEL found, keeping internal SUSY switched off                       | 
 |     21   Warning in SpaceShower::pT2nextQCD: small daughter PDF                                                 | 
 |    104   Warning in SpaceShower::pT2nextQCD: weight above unity                                                 | 
 |    256   Warning in StringFragmentation::fragmentToJunction: bad convergence junction rest frame                | 
 |    274   Warning in TauDecays::decay: unknown tau production, assuming unpolarized and uncorrelated             | 
 |    442   Warning in TimeShower::findMEcorr: ME weight above PS one                                              | 
 |                                                                                                                 | 
 *-------  End PYTHIA Error and Warning Messages Statistics  ------------------------------------------------------* 

*/
