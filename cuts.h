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
double const tmass = 175.0;
double const ptlepton = 0.0;
double const lepiso = 0.0;
double const MeeMax = 3000.0;
double const MetMin = 0.0;
