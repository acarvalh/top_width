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
    data = ".lhe.decayed";
    //data = ".lhe.shower";
    double cut[10] = {180,185, 190, 195, 200, 210, 220, 230, 240, 250};
    for(unsigned int ifolder=0; ifolder<1;ifolder++) // { four samples by folder
        // to each cut deffinition and each one of the four region deffintions I pass by the four files and then save
        for(unsigned int mtdef=2; mtdef<3; mtdef++)
        for(unsigned int type=0; type<1;type++) { // for cut deffinition
            decla(0);
            // had/lep | plus/minus
            // fill if: 0 = (< mt1,mt2) | 1 = (< m1 , >m2) | 2 = (>m1 , m2<) | 3 = (m1,m2 >) 
            for(unsigned int isample=0; isample<3;isample++){ 
            // samples with weights 
                double weight; double nevents =100000, lumi = 20;// /fb
                double finalevents; // counter 
                double CX[5] = {22.67,1.103,1.103,100, 24.92};// 1,1,1,1,1};//
            for (unsigned i=0; i<5; i++ ) if(isample==i) weight = CX[i]*lumi/nevents; // decide the mt cut    
                file = path[ifolder] + sample[isample]+ data;
                //////////////////////////////////
                cout<<"\n\n reading file = "<<file<<endl;
                ifstream in1;
                in1.open(file.c_str());
                //return 0;
                for(unsigned int ievent=0;ievent<100000;ievent++){ // for each event  //
                    string c;
                    in1>>c;
                    double Px, Py , Pz, E;
                    int pID, mother;
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
                            leptons.at(counterl).set_user_index(pID); 
                            counterl++;
                        } else if (abs(pID)==12 || abs(pID)==14) {
                            neutrinos.push_back(fastjet::PseudoJet(Px,Py,Pz,E));
                            neutrinos.at(countern).set_user_index(pID);
                            countern++;
                        }
                        /////////////////////////////////////////////////
                    } // close for each particle
                    //if(ievent== 32634) {
                    //    cout<<"jet1 "<<particles.at(0).px()<<" "<<particles.at(0).py()<<" "<<particles.at(0).pz()<<" "<<particles.at(0).e()<<endl;
                    //    cout<<"jet2 "<<particles.at(1).px()<<" "<<particles.at(1).py()<<" "<<particles.at(1).pz()<<" "<<particles.at(1).e()<<endl;
                    //    cout<<"jet3 "<<particles.at(2).px()<<" "<<particles.at(2).py()<<" "<<particles.at(2).pz()<<" "<<particles.at(2).e()<<endl;
                    //    cout<<"jet4 "<<particles.at(3).px()<<" "<<particles.at(3).py()<<" "<<particles.at(3).pz()<<" "<<particles.at(3).e()<<endl;
                    //    cout<<"lepton "<<leptons.at(0).px()<<" "<<leptons.at(0).py()<<" "<<leptons.at(0).pz()<<" "<<leptons.at(0).e()<<endl;
                    //    cout<<"neutrino "<<particles.at(0).px()<<" "<<particles.at(0).py()<<" "<<particles.at(0).pz()<<" "<<particles.at(0).e()<<endl;
                    //}
                    vector<PseudoJet> jets;       // cout<<"here "<< counterl <<endl;
                    vector<int> btag, bmistag, fattag, btrue; int bh,bl; double met=0;
                    int njets = recojets(particles, jets,btag,bmistag,fattag,btrue);
                    int numbb=0; for(unsigned i =0; i< btrue.size() ; i++ ) if(abs(btrue[i])==5) numbb++; // "b--taggable gen-b's
                    //cout<<"here "<< numbb <<endl;
                    // lepton isolation and basic cuts
                    int nlep =recol(jets,leptons,neutrinos); // returned the leptons 
                    bool lepcuts=false, hadtopreco=false, lepwreco=false;
                    int count=0; for(unsigned int i = 0; i < njets; i++) if(btag[i]>0)count++;
                    if( semilep && nlep>0 && njets>3 && numbb >1// only two jets to reco the hadronic W and make enought balance to reco met
                       //&& count>1
                       ){ cout<<"njets "<<njets<<" nleptons "<<counterl<<" pzl "<< neutrinos.at(0).pz()<<endl;
                        //int true_tops = truetops(jets,leptons,neutrinos,btag,btrue);  
                        // NEED TO RE-IMPLEMENT THE MET cut!!!!
                        if(lepcuts && reco == 0) hadtopreco=recohadt(bh,bl,jets,leptons,neutrinos,btag,btrue,met,weight,cut[mtdef],type); // && true_tops==type )
                        if(lepcuts && reco == 1) lepwreco = recolept2step(bh,bl,jets,leptons,neutrinos,btag,btrue,met,weight,cut[mtdef],type); // && true_tops==type )
                        //if(lepcuts && true_tops==type )lepwreco = recotlepeq(bh,bl,jets,leptons,neutrinos,btag,btrue,met, weight);
                    } else if ( !semilep && nlep>1 && njets>1 && numbb >1){
                        //cout<<"njets "<<njets<<" nleptons "<<counterl<<" pzl "<< neutrinos.size()<<endl;
                        hadtopreco = fullylep(bh,bl,jets,leptons,neutrinos,btag,btrue,met,weight,cut[mtdef],type); 
                        if(hadtopreco) finalevents=+ 
                        //cout<<"enterednjets "<<njets<<" nleptons "<<counterl<<" pzl "<< neutrinos.at(0).pz()<<endl;
                    } // close if !semilep
                    // To sum all and divide accordingly with a higher mt deffinition 
                    // we need to do a ttree insted                    
                } // close for each event
               in1.close();  
               } // close for sample
               cout<<"\n\n closing file = "<<file<<endl; //cout<<ifolder<<endl;
               save_hist(1,ifolder,type);
               if(type==0) cout<<"OnOn, mtcut ="<< cut<<endl;
               //cout<<"\n\n closing file = "<<endl; 
        } // close fo type
    }