//////////////////
// to run:
// make HH_VBF
// ./HH_VBF
/////////////////
/*
mg5>define bb = b b~
Defined multiparticle bb = b b~
mg5>define ww = w+ w-
Defined multiparticle ww = w+ w-
mg5>define ll = e+ e- mu+ mu-
Defined multiparticle ll = e- mu- e+ mu+
mg5>define nu = vl vl~
Defined multiparticle nu = ve vm vt ve~ vm~ vt~
mg5>define j = g u c d s u~ c~ d~ s~ b b~
Defined multiparticle j = g u c d s u~ c~ d~ s~ b b~
mg5>generate p p > t t~, (tt > ww bb, w+ > ll nu, w- > j j )
Command "generate p p > t t~, (tt > ww bb, w+ > ll nu, w- > j j )" interrupted with error:
InvalidCmd : No particle tt in model
mg5>define tt = t t~
Defined multiparticle tt = t t~
mg5>generate p p > t t~, (tt > ww bb, w+ > ll nu, w- > j j ) on-on: 13 diagrams ?
generate p p > t w- b~, w- > ll nu, (t > w+ b, w+ > j j) on-off:  21 diagrams
generate p p > t~ w+ b, w+ > j j , (t~ > w- b~, w- > ll nu) off-on: 21 diagrams
*/
#include "fastjet/ClusterSequence.hh"
#include "fastjet/Selector.hh"
#include "fastjet/tools/MassDropTagger.hh"
#include <fstream>
using namespace fastjet;
using namespace std;
#include "Functions.h"
#include "choices.h"
int main() {
  srand( time(NULL) );
  hello();
  ////////////////////////////////
  // input
  vector<string> filename;
  string file, path,data;
  path="/home/xanda/Documents/ggAnalysis/Parton/code/top_width/";
  string sample[5] = {"onon","onoff","offon","offoff","wbwb_13tev_100k"};// the last one have 50k
  data = ".lhe.decayed";
  for(unsigned int isample=4; isample<5;isample++)
   for(unsigned int type=0; type<6;type++){
  decla(0);
  file = path + sample[isample]+ data;
  //////////////////////////////////
  ifstream in1;
  cout<<"\n\n reading file = "<<file<<endl;
  //return 0;
  in1.open(file.c_str());
  for(unsigned int ievent=0;ievent<100000;ievent++){ // for each event  // 
     string c;
     in1>>c;
     double Px, Py , Pz, E;
     int pID;
     unsigned int nparticles;
     vector<PseudoJet> particles;//jets 
     vector<PseudoJet> neutrinos;
     vector<PseudoJet> leptons; 
     vector<PseudoJet> tops;                    
     int nb = 0;
     in1>>nparticles; unsigned int counter=0,countert=0,counterl=0,countern=0;
       for(unsigned int ipart=0;ipart<nparticles;ipart++){ // loop on particles
          in1 >> pID >> Px >> Py >> Pz >> E ;//>> idup;
          if (abs(pID) < 6 || pID==21){  // if a quark/gluon -- neglect hadrons
		particles.push_back(fastjet::PseudoJet(Px,Py,Pz,E)); 
                particles.at(counter).set_user_index(pID); 
		if(abs(pID) == 5) nb++;  // count b's and no-b's
                //cout<<"particle flavour "<< particles.at(counter).user_index()<<endl;
		counter++;
	  } else if (abs(pID)==6) {
		tops.push_back(fastjet::PseudoJet(Px,Py,Pz,E));
		countert++;
          } else if (abs(pID)==11 || abs(pID)==13) {
		leptons.push_back(fastjet::PseudoJet(Px,Py,Pz,E));
		counterl++;
          } else if (abs(pID)==12 || abs(pID)==14) {
		neutrinos.push_back(fastjet::PseudoJet(Px,Py,Pz,E));
		countern++;
          }
       /////////////////////////////////////////////////
       } // close for each particle
       vector<PseudoJet> jets; 
       vector<int> btag, bmistag, fattag, btrue; int bh,bl; double met;
       int njets = recojets(particles, jets,btag,bmistag,fattag,btrue);
       //cout<<"here"<<endl;
       // lepton isolation and basic cuts
       bool lepcuts=false, hadtopreco=false, lepwreco=false;
       int count=0; for(unsigned int i = 0; i < njets; i++) if(btag[i]>0)count++;
       if(counterl==1 && njets>3 // only two jets to reco the hadronic W and make enought balance to reco met 
          //&& count>1
         ){ //cout<<"njets "<<njets<<" nleptons "<<counterl<<" pzl "<< neutrinos.at(0).pz()<<endl;
         lepcuts=recol(met,jets,leptons,neutrinos); // returned the leptons
         int true_tops = truetops(jets,leptons,neutrinos,btag,btrue);  
              // by now it only works at parton level
         if(lepcuts && reco == 0 && true_tops==type ) hadtopreco=recohadt(bh,bl,jets,leptons,neutrinos,btag,btrue,met); 
         if(lepcuts && reco == 1 && true_tops==type ) lepwreco = recolept2step(bh,bl,jets,leptons,neutrinos,btag,btrue,met); 
         //if(lepcuts && true_tops==type )lepwreco = recotlepeq(bh,bl,jets,leptons,neutrinos,btag,btrue,met);
       }
  } // close for each event
  cout<<"\n\n closing file = "<<file<<endl;
  save_hist(isample,reco,type);
  in1.close();
  } // close fo sample
}
