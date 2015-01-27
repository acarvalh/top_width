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
    //vector<string> filename;
    string file, data;
    // we are going to have only four samples, with 180 GeV mt selection, merge them and separate
    string path[1]={
        "/afs/cern.ch/work/a/acarvalh/tt_fulylep/ttbar_gen180/"
        //"/data2/TopWidth/MG5_aMC_v2_1_0/wwbbfullep/top_Wvary/",
        //"/data2/TopWidth/MG5_aMC_v2_1_0/wwbbfullep_OnOn/top_Wvary/",
        //"/data2/TopWidth/MG5_aMC_v2_1_0/wwbbfullep_OnOff/top_Wvary/",
        //"/data2/TopWidth/MG5_aMC_v2_1_0/wwbbfullep_OffOn/top_Wvary/",
        //"/data2/TopWidth/MG5_aMC_v2_1_0/wwbbfullep_OffOff/top_Wvary/"
    };
    //string sample[7] = {"Wt_0","Wt_1","Wt_2","Wt_3","Wt_4","Wt_5","Wt_6"};// the last one have 50k
    string sample[5] = {"ttOnOnLep180","ttOnOffLep180","ttOffOnLep180","ttOffOffLep180","ttFullLep180"};// the last one have 50k    
    string label[5] = {"OnOn","OnOff","OffOn","OffOff","Full"};// the last one have 50k  
    data = ".lhe.decayed";
    //data = ".lhe.shower";
    double cut[10] = {180,185, 190, 195, 200, 210, 220, 230, 240, 250};
    /////////////////////////////////////////////////////////
    // information I want to make a table
    //
    // mtdef sample type CX net_eff
    // just need net_eff / sample / mtdeff ===> vector
    /////////////////////////////////////////////////////////
    // gen ifo deffinitions
    //
    // had/lep | plus/minus
    // fill if: 0 = (< mt1,mt2) | 1 = (< m1 , >m2) | 2 = (>m1 , m2<) | 3 = (m1,m2 >) 
    /////////////////////////////////////////////////////////
    double CX[5] = {22.67,1.103,1.103,100, 24.92};// 1,1,1,1,1};//
    double nevents =100000, lumi = 50;// /fb
    vector< vector< double > > finalevents; // save nevents / file / mtdeff
    vector< vector< double > > finaleventsfrom; // trace / file / mtdeff
    // to each cut deffinition and each one of the four region deffintions I pass by the four files and then save 
    for(unsigned int mtdef=2; mtdef<3; mtdef++) // for gen cut deffinition
      for(unsigned int isample=0; isample<1;isample++)
        for(unsigned int type=0; type<1;type++) { 
             decla(0);
             double finalevents=0; // counter for net eff
             ///////////////////////////////////////////////
             double weight; 
             for(unsigned i=0; i<5; i++ ) if(isample==i && nicepic) weight = CX[i]*lumi/nevents;   
               else if(isample==i && !nicepic) weight = CX[i]*lumi/nevents; // this one makes raw efficiencies
             //////////////////////////////////////////////
             file = path[0] + sample[isample]+ data;
             cout<<"\n\n reading file = "<<file<<endl;
             ifstream in1; in1.open(file.c_str());
             for(unsigned int ievent=0;ievent<100000;ievent++){ // read and process for each event 
                double Px, Py , Pz, E; int pID, mother; unsigned int nparticles;
                unsigned int counter=0,countert=0,counterl=0,countern=0, nb = 0;
                vector<PseudoJet> particles;//jets 
                vector<PseudoJet> neutrinos;
                vector<PseudoJet> leptons; 
                vector<PseudoJet> tops;                    
                string c; in1>>c; in1>>nparticles; 
                /////////////////////////////////////////////////////////////////////
                // read and understand
                for(unsigned int ipart=0;ipart<nparticles;ipart++){ // loop on particles
                  in1 >> pID >> Px >> Py >> Pz >> E ;//>> idup;
                  if(abs(pID) < 6 || pID==21){ // quark / gluon
                      particles.push_back(fastjet::PseudoJet(Px,Py,Pz,E)); 
                      particles.at(counter).set_user_index(pID); //cout<<"particle flavour "<< particles.at(counter).user_index()<<endl;
                    if(abs(pID) == 5) nb++; counter++; // count b's and no-b's
                    } else if (abs(pID)==6) {tops.push_back(fastjet::PseudoJet(Px,Py,Pz,E)); countert++; // b--quarks
                    } else if (abs(pID)==11 || abs(pID)==13) {
                      leptons.push_back(fastjet::PseudoJet(Px,Py,Pz,E)); // counterl++;
                      leptons.at(counterl).set_user_index(pID); // save charge for gen deffinition
                    } else if (abs(pID)==12 || abs(pID)==14) {
                      neutrinos.push_back(fastjet::PseudoJet(Px,Py,Pz,E)); countern++;
                      neutrinos.at(countern).set_user_index(pID);
                    } // close if 
                  } // close for each particle
                  ////////////////////////////////////////////////////////////////////
                  // Analyse
                  vector<PseudoJet> jets; vector<int> btag, bmistag, fattag, btrue; int bh,bl; double met=0;
                  bool lepcuts=false, hadtopreco=false, lepwreco=false;
                  int numbb=0; for(unsigned i =0; i< btrue.size() ; i++ ) if(abs(btrue[i])==5) numbb++; // "b--taggable gen-b's
                  //
                  int nlep =recol(jets,leptons,neutrinos); // lepton isolation and basic cuts  
                  int njets = recojets(particles, jets,btag,bmistag,fattag,btrue); // jet deffinition and basic cuts
                  if(semilep && nlep>0 && njets>3 && numbb >1){ 
                    // NEED TO RE-IMPLEMENT THE MET cut!!!!
                    if(lepcuts && reco == 0) hadtopreco=recohadt(bh,bl,jets,leptons,neutrinos,btag,btrue,met,weight,cut[mtdef],type); // reco by had
                    if(lepcuts && reco == 1) lepwreco = recolept2step(bh,bl,jets,leptons,neutrinos,btag,btrue,met,weight,cut[mtdef],type); // by lep
                    }else if(!semilep && nlep>1 && njets>1 && numbb >1){ // close if semilep
                      hadtopreco = fullylep(bh,bl,jets,leptons,neutrinos,btag,btrue,met,weight,cut[mtdef],type); if(hadtopreco)finalevents++; 
                    } // close if !semilep
              } in1.close(); // close for each event
              save_hist(1,imtdef,type); // hitogram / type / mtdef
              ///////////////////////////////////////////////////////////////////// 
              // save relevant info
              finalevents0[isample] = finalevents;
              double neventslumi = (finalevents0[0]*CX[0]+finalevents0[1]*CX[1]+finalevents0[2]*CX[2]+finalevents0[3]*CX[3])*1000*lumi/nevents;
              double neventsall = (CX[0]+CX[1]+CX[2])*1000*lumi;
              ///////////////////////////////////////////////////////////////////
              // intermediate check
              if(!nicepic)cout<<" raw nevents passed = "<< finalevents0[isample]<<" sample "<<isample<<endl;               
              cout<<"\n\n closing file = "<<file<<endl; cout<<" "<<endl;
              cout<<"mt cut ="<< cut[mtdef]<<endl;
        } // close fo type
        //////////////////////////////////////////////////////////////////////////
        // make the table 
    
    }