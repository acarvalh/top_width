//////////////////
// to run:
// source /afs/cern.ch/sw/lcg/app/releases/ROOT/5.32.00/x86_64-slc5-gcc43-opt/root/bin/thisroot.sh
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
    string fileparton , fileshower , dataparton = ".lhe.decayed" , datashower;
    if(showering) datashower = ".lhe.shower"; else datashower = ".lhe.decayed"; 
    //////////////////////////////////////////////////
    // input
    // we are going to have only four samples, with 180 GeV mt selection, merge them and separate
//    string path[1]={"/afs/cern.ch/work/a/acarvalh/tt_fulylep/"};
    string path[1]={"/Users/Xanda/Documents/codes/git/top_width/mtdef/"};
    string sample[5] = {"OnOnVary","OnOffVary","OffOnVary","OffOffVary","FullVary"};// folder
    string fileGam[13] = {"/Wt_0","/Wt_1","/Wt_2","/Wt_3","/Wt_4","/Wt_5","/Wt_6","/Wt_7","/Wt_8","/Wt_9","/Wt_10","/Wt_11","/Wt_12"}; // width
    double Gamma[13] ={0.6, 0.8, 1.01, 1.20, 1.40, 1.608, 1.8, 2, 2.2, 2.4, 2.6, 2.8, 3};
    string label[5] = {"OnOn","OnOff","OffOn","OffOff","Full"};
    // our files are now all showered --- if there are less than 20 particles it is genlevel + EW radiation
    // the list is ordered by status, so it is enought to take the first 6 particles and neglect FSR kiks as well IS remnants
    string file , data = ".lhe.shower"; 
    int maxnpart = -1 , minnpart = -1; if(showering) {maxnpart = 10000; minnpart = 9;}  else {maxnpart = 15; minnpart = 5;}  
    cout<<" maxnpart "<<maxnpart<<endl;
    double cut[4] = {190,250,300,350};//180,185, 190, 195, 200, 210, 220, 230, 240, 250};
    // If shower I read the two files
    /////////////////////////////////////////////////////////
    // information I want to make a table
    //
    // mtdef sample type CX net_eff
    // just need net_eff / sample / mtdeff ===> vector
    //vector< vector< vector< double > > > finaleventsN; // save nevents / sample / type / mtdeff
    //vector<vector< double > > finaleventsfrom[4]; // trace / sample / type / mtdeff
    vector< vector< vector<  vector<double> > > > finaleventsN(13, 
                                                        vector< vector< vector<double> > > (10, 
                                                                vector< vector<double> >(5, 
                                                                        vector<double>(4) ) ) ); // trace / sample / type / mtdeff
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    vector< vector< vector<  vector<double> > > > finaleventsN0(13, 
                                                               vector< vector< vector<double> > > (10, 
                                                                                                   vector< vector<double> >(5, 
                                                                                                                            vector<double>(4) ) ) ); // trace / sample / type / mtdeff
    vector< vector< vector<  vector<double> > > > finaleventsN1(13, 
                                                               vector< vector< vector<double> > > (10, 
                                                                                                   vector< vector<double> >(5, 
                                                                                                                            vector<double>(4) ) ) ); // trace / sample / type / mtdeff
    vector< vector< vector<  vector<double> > > > finaleventsN2(13, 
                                                               vector< vector< vector<double> > > (10, 
                                                                                                   vector< vector<double> >(5, 
                                                                                                                            vector<double>(4) ) ) ); // trace / sample / type / mtdeff
    vector< vector< vector<  vector<double> > > > finaleventsN3(13, 
                                                               vector< vector< vector<double> > > (10, 
                                                                                                   vector< vector<double> >(5, 
                                                                                                                            vector<double>(4) ) ) ); // trace / sample / type / mtdeff
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    vector< vector< vector<  vector<double> > > > finaleventsN1D0(13, 
                                                                vector< vector< vector<double> > > (10, 
                                                                                                    vector< vector<double> >(5, 
                                                                                                                             vector<double>(4) ) ) ); // trace / sample / type / mtdeff
    vector< vector< vector<  vector<double> > > > finaleventsN1D1(13, 
                                                                vector< vector< vector<double> > > (10, 
                                                                                                    vector< vector<double> >(5, 
                                                                                                                             vector<double>(4) ) ) ); // trace / sample / type / mtdeff
    vector< vector< vector<  vector<double> > > > finaleventsN1D2(13, 
                                                                vector< vector< vector<double> > > (10, 
                                                                                                    vector< vector<double> >(5, 
                                                                                                                             vector<double>(4) ) ) ); // trace / sample / type / mtdeff
    vector< vector< vector<  vector<double> > > > finaleventsN1D3(13, 
                                                                vector< vector< vector<double> > > (10, 
                                                                                                    vector< vector<double> >(5, 
                                                                                                                             vector<double>(4) ) ) ); // trace / sample / type / mtdeff
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    vector< vector< vector<  vector<int> > > > finaleventsfrom(13,
                                                               vector< vector< vector<int> > > (10, 
                                                                        vector< vector<int> >(5, 
                                                                                vector<int>(4) ) ) ); // trace / sample / type / mtdeff                                                               
    /////////////////////////////////////////////////////////
    // gen ifo deffinitions
    //
    // had/lep | plus/minus
    // fill if: 0 = (< mt1,mt2) | 1 = (< m1 , >m2) | 2 = (>m1 , m2<) | 3 = (m1,m2 >) 
    /////////////////////////////////////////////////////////
    double CX[5] = {25.8,1.177,1.176,0.0572, 24.92};// 1,1,1,1,1};// Fix CX by Gamma
    double nevents =100000, lumi = 50;// /fb
    //////////////////////////////////////////////////////////////////////////////////
    // to each cut deffinition and each one of the four region deffintions I pass by the four files and then save 
    for(unsigned int files=0; files<13; files++) // 13 width
    for(unsigned int mtdef=0; mtdef<2; mtdef++) // 10 for gen cut deffinition
      for(unsigned int type=0; type<4;type++) { // OnOn ... 
          double finaleventscheck[4]; // to be obsolete
          decla(0);
           for(unsigned int isample=0; isample<4;isample++) { // ononFile ....              
             double finalevents=0; 
             double finalevents0=0; double finalevents1=0; double finalevents2=0; double finalevents3=0; // counter for net eff
             double finalevents1D0=0; double finalevents1D1=0; double finalevents1D2=0; double finalevents1D3=0; // counter for net eff
             ///////////////////////////////////////////////
             double weight; 
             for(unsigned i=0; i<5; i++ ) 
               {if(isample==i && nicepic) weight = CX[i]*lumi/nevents; else if(isample==i && !nicepic) weight = lumi/nevents;} // this one makes raw efficiencies
             //////////////////////////////////////////////
             fileparton = path[0] + sample[isample] + fileGam[files] + dataparton; cout<<"\n\n reading file = "<<fileparton<<endl;
             fileshower = path[0] + sample[isample] + fileGam[files] + datashower; cout<<"\n\n reading file = "<<fileshower<<endl;
             ifstream in1; in1.open(fileparton.c_str());
             ifstream in1shower; in1shower.open(fileshower.c_str());
             //
             for(unsigned int ievent=0;ievent<100000;ievent++){ // read and process for each event 
                double Pxgen, Pygen , Pzgen, Egen; int pIDgen, mothergen; unsigned int nparticlesgen;
                unsigned int countergen=0,countertgen=0,counterlgen=0,counterngen=0, nbgen = 0;
                vector<PseudoJet> particlesgen; //jets 
                vector<PseudoJet> neutrinosgen;
                vector<PseudoJet> leptonsgen; 
                vector<PseudoJet> topsgen;   
                //
                double Px, Py , Pz, E; int pID, mother; unsigned int nparticles;
                unsigned int counter=0,countert=0,counterl=0,countern=0, nb = 0;
                vector<PseudoJet> particles; //jets 
                vector<PseudoJet> neutrinos;
                vector<PseudoJet> leptons; 
                vector<PseudoJet> tops;   
                //
                /////////////////////////////////////////////////////////////////////
                // read and understand
                string cgen; in1>>cgen; in1>>nparticlesgen; 
                for(unsigned int ipart=0;ipart<nparticlesgen;ipart++){ // loop on particles
                   in1 >> pIDgen >> Pxgen >> Pygen >> Pzgen >> Egen ;//>> idup;
                   if(abs(pIDgen)<6||pIDgen==21){particlesgen.push_back(fastjet::PseudoJet(Pxgen,Pygen,Pzgen,Egen));particlesgen.at(countergen).set_user_index(pIDgen);if(abs(pIDgen)==5)nbgen++;countergen++;} 
                   else if (abs(pIDgen)==6) {topsgen.push_back(fastjet::PseudoJet(Pxgen,Pygen,Pzgen,Egen)); countertgen++;} 
                   else if (abs(pIDgen)==11 || abs(pIDgen)==13) {leptonsgen.push_back(fastjet::PseudoJet(Pxgen,Pygen,Pzgen,Egen));leptonsgen.at(counterlgen).set_user_index(pIDgen); counterlgen++;} 
                   else if (abs(pIDgen)==12 || abs(pIDgen)==14) {neutrinosgen.push_back(fastjet::PseudoJet(Pxgen,Pygen,Pzgen,Egen));neutrinosgen.at(counterngen).set_user_index(pIDgen);counterngen++;} 
                 } // close for each particle
                string cshower; in1shower>>cshower; in1shower>>nparticles; 
                 // read showered
                 for(unsigned int ipart=0;ipart<nparticles;ipart++){ // loop on particles
                     in1shower >> pID >> Px >> Py >> Pz >> E ;//>> idup;
                     if(abs(pID) < 6 || pID==21){particles.push_back(fastjet::PseudoJet(Px,Py,Pz,E));particles.at(counter).set_user_index(pID); if(abs(pID) == 5) nb++; counter++;} 
                     else if (abs(pID)==6) {tops.push_back(fastjet::PseudoJet(Px,Py,Pz,E)); countert++;} 
                     else if (abs(pID)==11 || abs(pID)==13) {leptons.push_back(fastjet::PseudoJet(Px,Py,Pz,E)); leptons.at(counterl).set_user_index(pID); counterl++;} 
                     else if (abs(pID)==12 || abs(pID)==14) {neutrinos.push_back(fastjet::PseudoJet(Px,Py,Pz,E)); neutrinos.at(countern).set_user_index(pID); countern++;} // close if 
                 } // close for each particle
                 ////////////////////////////////////////////////////////////////////
                 // Parton level cuts 
                 bool GenLevel = GenLevelDilep(particlesgen, leptonsgen, neutrinosgen, cut[mtdef], weight, type); //cout<<"GenLevel "<<GenLevel<<endl;
                 ////////////////////////////////////////////////////////////////////
                 // Analyse
                 int reco1D1=-1 , reco1D2=-1 , reco2D=-1;
                 if(GenLevel && nparticles > minnpart ){
                    //cout<<ievent<<" nparticles "<<nparticles<<endl; 
                    //cout<<ievent<<" nparticlesgen "<<nparticlesgen<<endl;
                    vector<PseudoJet> jets; vector<int> btag, bmistag, fattag, btrue; int bh,bl; double met=0;
                    bool lepcuts=false, lepwreco=false; int hadtopreco;
                    int njets = recojets(particles, jets,btag,bmistag,fattag,btrue,weight); // jet deffinition and basic cuts
                    int nlep =recol(jets,leptons,neutrinos,weight); // lepton isolation and basic cuts  
                    //int numbb=0;
                    //if(!showering) {for(unsigned i =0; i< btrue.size() ; i++ ) if(abs(btrue[i])==5) numbb++;} else numbb = 5; // "b--taggable gen-b's
                    //cout<<"njets "<<njets<<" nlep "<<nlep<<endl;
                    if(semilep && nlep>0 && njets>3){ 
                    // NEED TO RE-IMPLEMENT THE MET cut!!!!
                    //if(lepcuts && reco == 0) hadtopreco=recohadt(bh,bl,jets,leptons,neutrinos,btag,btrue,met,weight,cut[mtdef],type); // reco by had
                    //if(lepcuts && reco == 1) lepwreco = recolept2step(bh,bl,jets,leptons,neutrinos,btag,btrue,met,weight,cut[mtdef],type); // by lep
                      cout<<"ops not semileptonic "<<endl;
                     } else if(!semilep && nlep>1 && njets>1){ 
                         hadtopreco = fullylep(reco1D1,reco1D2,reco2D,jets,leptons,neutrinos,btag,btrue,met,weight); 
                         if(hadtopreco>-1)finalevents++;
                         if(reco2D==0)finalevents0++; else if(reco2D==1) finalevents1++; else if(reco2D==2) finalevents2++; else if(reco2D==3) finalevents3++; 
                         if(reco1D1==0)finalevents1D0++; else if(reco1D1==1) finalevents1D1++; 
                         if(reco1D2==2) finalevents1D2++; else if(reco1D2==3) finalevents1D3++;
                     } // close fullylep 
                 } // close GenLevel
                } in1.close(); // close for each event
              /////////////////////////////////////////////////////////////////////
              // intermediate check
              finaleventscheck[isample] = finalevents; 
              double neventslumi = (finaleventscheck[0]*CX[0]+finaleventscheck[1]*CX[1]+finaleventscheck[2]*CX[2]+finaleventscheck[3]*CX[3])*1000*lumi/nevents;
              //double neventsa9ll = (CX[0]+CX[1]+CX[2])*1000*lumi;
              cout<<" raw nevents passed = "<< finaleventscheck[isample]<<" sample "<<isample<<endl;               
              cout<<"\n\n closing file = "<<fileshower<<endl; //cout<<" "<<endl;
              cout<<"mt cut ="<< cut[mtdef]<<endl;   
              // save to table
              //int tablecounter=0;
              finaleventsN[files][mtdef][isample][type] = finalevents; // raw events always! 
              // for 2D couting
              finaleventsN0[files][mtdef][isample][type] = finalevents0; // raw events always!                
              finaleventsN1[files][mtdef][isample][type] = finalevents1; // raw events always! 
              finaleventsN2[files][mtdef][isample][type] = finalevents2; // raw events always!           
              finaleventsN3[files][mtdef][isample][type] = finalevents3; // raw events always!           
              // for 1D couting
              finaleventsN1D0[files][mtdef][isample][type] = finalevents1D0; // raw events always!                
              finaleventsN1D1[files][mtdef][isample][type] = finalevents1D1; // raw events always! 
              finaleventsN1D2[files][mtdef][isample][type] = finalevents1D2; // raw events always!           
              finaleventsN1D3[files][mtdef][isample][type] = finalevents1D3; // raw events always!           
               //finaleventsfrom[mtdef][isample][type] = isample + type; //finaleventsfrom.at(tablecounter).push_back(type);
              } // close for sample
              save_hist(files,mtdef,type,showering); // hitogram / type / mtdef ==> after pass by the 4 samples
        } // close fo type and mtcut and gam
        //////////////////////////////////////////////////////////////////////////
        // make the table 
        ofstream NetEffMtGen;
        NetEffMtGen.open("NetEffMtGen.txt"); // file to save
        cout<<" Gamm mtcut type sample NEtEv2D0 NEtEv2D1 NEtEv2D2 NEtEv2D3 NEtEv1D0 NEtEv1D1 NEtEv1D2 NEtEv1D3 total"<<endl;
        NetEffMtGen<<"Gamm mtcut type sample NEtEv2D0 NEtEv2D1 NEtEv2D2 NEtEv2D3 NEtEv1D0 NEtEv1D1 NEtEv1D2 NEtEv1D3 total"<<endl;
        int nev; if(showering) nev=50000; else nev=100000; 
        for(unsigned int files=0; files<13; files++)
            for(unsigned int mtdef=0; mtdef<2; mtdef++) for(unsigned int type=0; type<4;type++) for(unsigned int isample=0; isample<4;isample++) {
             cout<<        Gamma[files]<<" "<< cut[mtdef]<<" " <<label[type]<<" "<<sample[isample]
                           <<" "<<finaleventsN0[files][mtdef][isample][type]<<" "<<finaleventsN1[files][mtdef][isample][type]
                           <<" "<<finaleventsN2[files][mtdef][isample][type]<<" "<<finaleventsN3[files][mtdef][isample][type] 
                           <<" "<<finaleventsN1D0[files][mtdef][isample][type]<<" "<<finaleventsN1D1[files][mtdef][isample][type]
                           <<" "<<finaleventsN1D2[files][mtdef][isample][type]<<" "<<finaleventsN1D3[files][mtdef][isample][type] 
                           <<" "<<finaleventsN[files][mtdef][isample][type]<<endl;
             NetEffMtGen<< Gamma[files]<<" "<< cut[mtdef]<<" " << type      <<" "<<        isample
                          <<" "<<finaleventsN0[files][mtdef][isample][type]<<" "<<finaleventsN1[files][mtdef][isample][type]
                          <<" "<<finaleventsN2[files][mtdef][isample][type]<<" "<<finaleventsN3[files][mtdef][isample][type] 
                          <<" "<<finaleventsN1D0[files][mtdef][isample][type]<<" "<<finaleventsN1D1[files][mtdef][isample][type]
                          <<" "<<finaleventsN1D2[files][mtdef][isample][type]<<" "<<finaleventsN1D3[files][mtdef][isample][type] 
                          <<" "<<finaleventsN[files][mtdef][isample][type]<<endl;
        //cout<<" "<<endl;
        //cout<<" "<<mtdef<<" "<< type<<" "<<isample<<" "<<finaleventsfrom[mtdef][isample][type]<<endl;
        } //
        NetEffMtGen.close();
    }