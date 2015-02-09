/////////////////////////////////
// cuts
//double weight =1.;///10000;//0.001;//
double Mjj =0;//400; 
double PTjj = 0;//400; 
double DeltayVBF = 0;//3;
double DeltaRVBF = 0;//3;
bool shower=true;
// To be applied only to hadron level events
//
double const genmasshad=2000; // genmass
double const genmasshadmin=0; // genmass
//
double const genmasslep=2000; // genmass
double const genmasslepmin=0; // genmass
// the gen-level cuts
double const bjetpt = 30.0; 
double const mbblow = 0.0; 
// basline
double const jet_ptmin=35.0; // parton for jet reconstruction
double const jet_ptminfinal=20.0; // in final jet reconstruction
double const rapmax=5.0; // for jet reconstruction
double const etab = 2.5;
double const etal = 2.5;
double const etaj=15;
double const RR =0.4;
double const ptlepton = 20.;
double const lepiso = 0.3;
// analysis cuts
double const mblcut = 190;
int const cat =2; // minimum number of btag
////////////////////////////////////////
// weights b-tag
double const subjet2b=1;
double const fatjet2b=1;
double const normalb=1;
double const normalc=1;
double const normall=1;
double const misb=1;
/////////////////////////////////////////
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
//////////////////////////////////
// on the wwbb analysis
double const wmass = 80.4;
double const bmass = 4.7;
double const tmass = 173.0;
double const higgs_mass = 125.0;
double const MeeMax = 30000.0;
double const MetMin = 0.0;
double const  wbtransmassMax = 3000;
double const  wbtransmassMin = 0;

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