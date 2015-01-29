//////////////////
// to run:
// make HH_VBF
// ./HH_VBF
/////////////////
/*
 To macht gen numbering scheeme
 
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
//////////////////////////////////////////////////////////
#include "fastjet/ClusterSequence.hh"
#include "fastjet/Selector.hh"
#include "fastjet/tools/MassDropTagger.hh"
#include <fstream>
#include <iostream>
#include <vector>
using namespace fastjet;
using namespace std;
#include "Functions.h"
#include "choices.h"
int main() {
    srand( time(NULL) );
    hello();
    string file, data;
    //////////////////////////////////////////////////
    // input
    // we are going to have only four samples, with 180 GeV mt selection, merge them and separate
    string path[1]={"/afs/cern.ch/work/a/acarvalh/tt_fulylep/ttbar_gen180/"};
    string sample[5] = {"ttOnOnLep180","ttOnOffLep180","ttOffOnLep180","ttOffOffLep180","ttFullLep180"};// the last one have 50k    
    string label[5] = {"OnOn","OnOff","OffOn","OffOff","Full"}; 
    if(!showering)data = ".lhe.decayed"; else data = ".lhe.shower";
    double cut[10] = {180,185, 190, 195, 200, 210, 220, 230, 240, 250};
    /////////////////////////////////////////////////////////
    // information I want to make a table
    //
    // mtdef sample type CX net_eff
    // just need net_eff / sample / mtdeff ===> vector
    //vector< vector< vector< double > > > finaleventsN; // save nevents / sample / type / mtdeff
    //vector<vector< double > > finaleventsfrom[4]; // trace / sample / type / mtdeff
    vector< vector<  vector<double> > > finaleventsN(10, vector< vector<double> >(5, vector<double>(4))); // trace / sample / type / mtdeff
    vector< vector<  vector<int> > > finaleventsfrom(10, vector< vector<int> >(5, vector<int>(4))); // trace / sample / type / mtdeff
    /////////////////////////////////////////////////////////
    // gen ifo deffinitions
    //
    // had/lep | plus/minus
    // fill if: 0 = (< mt1,mt2) | 1 = (< m1 , >m2) | 2 = (>m1 , m2<) | 3 = (m1,m2 >) 
    /////////////////////////////////////////////////////////
    double CX[5] = {22.67,1.103,1.103,100, 24.92};// 1,1,1,1,1};//
    double nevents =100000, lumi = 50;// /fb
    //////////////////////////////////////////////////////////////////////////////////
    // to each cut deffinition and each one of the four region deffintions I pass by the four files and then save 
    for(unsigned int mtdef=0; mtdef<10; mtdef++) // for gen cut deffinition
      for(unsigned int type=0; type<4;type++) { 
          double finalevents0[4]; // to be obsolete
          for(unsigned int isample=0; isample<5;isample++) {
             decla(0);
             double finalevents=0; // counter for net eff
             ///////////////////////////////////////////////
             double weight; 
             for(unsigned i=0; i<5; i++ ) 
               {if(isample==i && nicepic) weight = CX[i]*lumi/nevents;   
                   else if(isample==i && !nicepic) weight = lumi/nevents;} // this one makes raw efficiencies
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
                      leptons.push_back(fastjet::PseudoJet(Px,Py,Pz,E)); // 
                      leptons.at(counterl).set_user_index(pID); counterl++; // save charge for gen deffinition
                    } else if (abs(pID)==12 || abs(pID)==14) {
                      neutrinos.push_back(fastjet::PseudoJet(Px,Py,Pz,E)); 
                      neutrinos.at(countern).set_user_index(pID); countern++;
                    } // close if 
                  } // close for each particle
                  ////////////////////////////////////////////////////////////////////
                  // Analyse
                  vector<PseudoJet> jets; vector<int> btag, bmistag, fattag, btrue; int bh,bl; double met=0;
                  bool lepcuts=false, hadtopreco=false, lepwreco=false;
                  int nlep =recol(jets,leptons,neutrinos); // lepton isolation and basic cuts  
                  int njets = recojets(particles, jets,btag,bmistag,fattag,btrue); // jet deffinition and basic cuts
                  int numbb=0; for(unsigned i =0; i< btrue.size() ; i++ ) if(abs(btrue[i])==5) numbb++; // "b--taggable gen-b's
                  if(semilep && nlep>0 && njets>3 && numbb >1){ 
                    // NEED TO RE-IMPLEMENT THE MET cut!!!!
                    if(lepcuts && reco == 0) hadtopreco=recohadt(bh,bl,jets,leptons,neutrinos,btag,btrue,met,weight,cut[mtdef],type); // reco by had
                    if(lepcuts && reco == 1) lepwreco = recolept2step(bh,bl,jets,leptons,neutrinos,btag,btrue,met,weight,cut[mtdef],type); // by lep
                    }else if(!semilep && nlep>1 && njets>1 && numbb >1){ // close if semilep
                      hadtopreco = fullylep(bh,bl,jets,leptons,neutrinos,btag,btrue,met,weight,cut[mtdef],type); 
                        if(hadtopreco)finalevents++;// else if(type ==0 && isample==0) cout<<"ops"<<endl; 
                        //cout<<"here"<<endl;
                    } // close if !semilep
                } in1.close(); // close for each event
              /////////////////////////////////////////////////////////////////////
              // intermediate check
              finalevents0[isample] = finalevents;
              double neventslumi = (finalevents0[0]*CX[0]+finalevents0[1]*CX[1]+finalevents0[2]*CX[2]+finalevents0[3]*CX[3])*1000*lumi/nevents;
              double neventsa9ll = (CX[0]+CX[1]+CX[2])*1000*lumi;
              if(!nicepic)cout<<" raw nevents passed = "<< finalevents0[isample]<<" sample "<<isample<<endl;               
              cout<<"\n\n closing file = "<<file<<endl; //cout<<" "<<endl;
              cout<<"mt cut ="<< cut[mtdef]<<endl;   
              // save to table
              int tablecounter=0;
              finaleventsN[mtdef][isample][type] = finalevents; // raw events always! 
              //finaleventsfrom[mtdef][isample][type] = isample + type; //finaleventsfrom.at(tablecounter).push_back(type);
              } // close for sample
              save_hist(1,mtdef,type); // hitogram / type / mtdef ==> after pass by the 4 samples
        } // close fo type and mtcut
        //////////////////////////////////////////////////////////////////////////
        // make the table 
        ofstream NetEffMtGen;
        NetEffMtGen.open("NetEffMtGen.txt"); // file to save
        cout<<"mtcut type sample NEtEv"<<endl;
        NetEffMtGen<<"mtcut type sample NEtEv"<<endl;
        for(unsigned int mtdef=0; mtdef<10; mtdef++) for(unsigned int type=0; type<2;type++) for(unsigned int isample=0; isample<5;isample++) {
        cout<<cut[mtdef]<<" " <<label[type]<<" "<<sample[isample]<<" "<<finaleventsN[mtdef][isample][type]<<endl;
        NetEffMtGen<< cut[mtdef]<<" " << type <<" "<< isample<<" "<<finaleventsN[mtdef][isample][type]<<endl;
        //cout<<" "<<endl;
        //cout<<" "<<mtdef<<" "<< type<<" "<<isample<<" "<<finaleventsfrom[mtdef][isample][type]<<endl;
        } //
        NetEffMtGen.close();
    }