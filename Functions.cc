#include <fstream>
#include "fastjet/ClusterSequence.hh"
#include "fastjet/tools/Filter.hh"
#include "fastjet/Selector.hh"
#include "fastjet/tools/MassDropTagger.hh"
#include <TCanvas.h>
#include <TFile.h>
#include <TArray.h>
#include <TH1D.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TArray.h>
#include <TVector.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <TString.h>
#include <TObject.h>
#include <TMath.h>
#include <TVector2.h>
#include <stdio.h>      /* Standard Library of Input and Output */
#include <complex.h>    /* Standard c++ Library of Complex Numbers */
#include "rootHistos.h" //declare the histos
#include "cuts.h" // basic cuts 
using namespace fastjet;
using namespace std;
/////////////////////////////////////////////////////////////////
void hello(){cout<<"\n\n\n HELLO!!!! \n\n"<<endl;}
/////////////////////////////////////////////////////////////////
bool GenLevelDilep(vector<PseudoJet> particles,  vector<PseudoJet> leptons, vector<PseudoJet> neutrinos , double cut, double weight, int type){
    ///////////////////////////////////////////////////////////////////
    // gen level info // had == plus
    bool passGencut = false;
    int blll=-1, bhhh=-1, ell=-1, ehh=-1, null=-1, nuhh=-1;
    for(unsigned int nj1=0; nj1< particles.size(); nj1++) if(particles.at(nj1).user_index()==5) bhhh=nj1; else if(particles.at(nj1).user_index()==-5) blll=nj1; 
    for(unsigned int nj1=0; nj1< leptons.size(); nj1++) {if(leptons.at(nj1).user_index() > 0) ell=nj1; else if(leptons.at(nj1).user_index() < 0) ehh=nj1;}
    for(unsigned int nj1=0; nj1< neutrinos.size(); nj1++) {if(neutrinos.at(nj1).user_index() > 0) nuhh=nj1; else  null=nj1;}
    PseudoJet lepTtrue, hadTtrue; int counttruth =-1;
    if(blll!=-1 && bhhh!=-1 && ell!=-1 && ehh!=-1 && null!=-1 && nuhh!=-1){ 
        lepTtrue = leptons.at(ell) + neutrinos.at(null) + particles.at(blll);
        hadTtrue = leptons.at(ehh) + neutrinos.at(nuhh) + particles.at(bhhh);
        //////////////////////////////////////////////////////////////////////////////
        // fill if: 0 = (< mt1,mt2) | 1 = (< m1 , >m2) | 2 = (>m1 , m2<) | 3 = (m1,m2 >)
        if(hadTtrue.m()<cut && lepTtrue.m() < cut) counttruth=0;
        else if(hadTtrue.m()<cut && lepTtrue.m()>=cut) counttruth=1;
        else if(hadTtrue.m()>=cut && lepTtrue.m()<cut) counttruth=2;            
        else if(hadTtrue.m()>=cut && lepTtrue.m()>=cut) counttruth=3;
        ////////////////////////////////////////////////////////////////////////////
        if(counttruth==type && particles.at(blll).pt() > bjetpt && particles.at(bhhh).pt() > bjetpt && (particles.at(blll)+particles.at(bhhh)).m() > mbblow ) {    
            genmbb->Fill((particles.at(bhhh)+particles.at(blll)).m(),weight);        
            leptop->Fill(lepTtrue.m(),weight);
            hadtop->Fill(hadTtrue.m(),weight);
            genbhad->Fill(particles.at(bhhh).pt(),weight);
            genblep->Fill(particles.at(blll).pt(),weight);
            genbhadeta->Fill(particles.at(bhhh).eta(),weight);
            genblepeta->Fill(particles.at(blll).eta(),weight);
            passGencut = true;
        } // close if cut
    } // close if have enough particles
    /////////////////////////////////////////////////////////////////////////
    return passGencut;
} // close gen level cuts
/////////////////////////////////////////////////////////////////
int fullylep(int & reco1D1, int & reco1D2, int & reco2D,vector<PseudoJet> jets, vector<PseudoJet> leptons_pt,vector<PseudoJet> neutrinos, vector<int> btag, vector<int> btrue, double met, double weight){
      bool passcut = false; int typeReco = -1;  // 0 = onon ; 1 = onoff ; 2 = offon ; 3 = offoff
      // I have at least 2 isolated leptons, I take the two leading: sort by pt
      int ell=-1, ehh=-1;  int blll =0, bhhh=1;
      vector<PseudoJet> leptons = sorted_by_pt(leptons_pt);
      unsigned int jsize = jets.size() , lsize = leptons.size();
      for(unsigned int nj1=0; nj1< lsize; nj1++) {if(leptons.at(nj1).user_index()>0 && ell==-1) ell=nj1; else if(leptons.at(nj1).user_index()<0 && ehh==-1) ehh=nj1;}
      /////////////////////////////////////////////////////////////////////////
      // start the analysis --- if passed the gencut
      PseudoJet MET = jets.at(0); for(unsigned n=1; n<jsize; n++) MET += jets.at(n); for(unsigned n=1; n<lsize; n++) MET += leptons.at(n);   
      met = MET.pt();// (neutrinos[0]+neutrinos[1]).pt(); // fix to MET
      // basic control plots + invariant mass of pairs
      double mll =  (leptons.at(0) + leptons.at(1)).m();
      double mjj = (jets.at(0) + jets.at(1)).m();
      ///////////////////////////////////////////////////
      // to use b--tag
      vector<PseudoJet> btaggedjets; 
      //for(unsigned n=0; n<jsize; n++) if(btrue[n] > 0) btaggedjets.push_back(jets.at(n));
      // not to use bt--tag
      btaggedjets = jets;
      if(btaggedjets.size() >1){
        // //////////////////////////////////////////////// works for shower
        // minimize to all jets
        vector<double> minplus , minminus; 
        for(unsigned int nj1=0; nj1< jsize; nj1++)  { //if(btag[nj1]>0) 
              //if(btag[nj1]>0) 
              // minimize invariant masses
              minplus.push_back((leptons.at(ehh) + btaggedjets.at(nj1)).m()); //minP.push_back(nj1); // I do not need to keep the jet, but anyway
              minminus.push_back((leptons.at(ell) + btaggedjets.at(nj1)).m()); //minM.push_back(nj1); // I do not need to keep the jet, but anyway 
        } // loop on jets
        int minM = TMath::LocMin(minminus.size(), &minminus[0]);
        int minP = TMath::LocMin(minplus.size(), &minplus[0]);
        //////////////////////////////////////////////////////////////          
        // minimize the balance instead
        vector<double> minBal; vector<int> j1v , j2v; int j1 , j2; //int minP , minM , minB; 
        for(unsigned int nj1=0; nj1< jsize; nj1++) for(unsigned int nj2=nj1+1; nj2< jsize; nj2++) { // if(btag[nj1]>0 || btag[nj2]>0) 
              //if(btag[nj1]>0) 
              minBal.push_back(abs((leptons.at(ell)+btaggedjets.at(nj1)).m()-(leptons.at(ehh)+btaggedjets.at(nj2)).m())); 
              j1v.push_back(nj1); j2v.push_back(nj2); // I need to keep the jet!! 
        } // loop on jets
        int minB = TMath::LocMin(minBal.size(), &minBal[0]);    
        ///////////////////////////////////////////////////
        // true instead
        //double mblLeadT= (leptons.at(ell) + jets.at(blll)).m(); double mblSubT= (leptons.at(ehh) + jets.at(bhhh)).m(); 
        ////////////////////////////////////////////////////
        if(minplus.size()>1 ){
          j1=j1v[minB]; j2=j2v[minB]; 
          bool minlj=true; double mbl1 = -10 , mbl2=-10;
          double mblBal1 = (leptons.at(ell)+btaggedjets.at(j1)).m() , mblBal2 = (leptons.at(ehh)+btaggedjets.at(j2)).m();
          if(minlj) { mbl1 = minplus[minP]; mbl2 = minplus[minM]; } else {mbl1 = mblBal1 ; mbl2 = mblBal2;}
          ///////////////////////////////////////////
          // 2D
          if( mbl1 > mblcut && mbl2 > mblcut ) {reco2D=3;} //if(truthMB ==1) {typeRecoTruth =6;} else typeRecoTruth =7;
          else if( mbl1 > mblcut && mbl2 < mblcut) {reco2D = 2;}
          else if( mbl1 < mblcut && mbl2 > mblcut) {reco2D = 1;} //if(truthMB ==1) {typeRecoTruth =3;} else typeRecoTruth =4;  
          else if( mbl1 < mblcut && mbl2 < mblcut) {reco2D = 0;} //if(truthMB ==1) {typeRecoTruth =0;} else typeRecoTruth =1; 
          else cout << "oups !!" <<endl;
          ///////////////////////////////////////////
          // 1D
          // balance
          if(mbl1 > mblcut) {reco1D2=3;} else if(mbl1 < mblcut) {reco1D2 = 2;} else cout << "oups !!" <<endl;
          if(mbl2 > mblcut) {reco1D1=1;} else if(mbl2 < mblcut) {reco1D1 = 0;} else cout << "oups !!" <<endl;
          ////////////////////////////////////////////
          //if(typeReco > -1 && jsize == 2){
          double vectorLep[14] = {leptons[0].pt(), leptons[0].eta(), leptons[1].pt(), leptons[1].eta(), met, mll, mjj, minplus[minP], minminus[minM], mblBal1, mblBal2, reco1D1, reco1D2, reco2D};
          for(unsigned i=0;i<14;i++) basicLeptons[i]->Fill(vectorLep[i],weight);
          // tranverse mass --- total transverse mass
          ///////////////////////////////////////////////////
          // we do not reco the tops --- fill with other vectors with zero
          double leptop[13]={170,0,0,0, //lept.m(),lept.pt(),lept.eta(),lept.phi(),
            90,0,0,0, //lepW.m(),lepW.pt(),lepW.eta(),lepW.phi(),
            0,0,0,0,0 //pnuzerror,truth,mterror,wmt,tmt};
          };
          for(unsigned i=1;i<13;i++) basicLeptop[i]->Fill(leptop[i],weight);
          double hadtop[10]={170,0,0,0, //hadt.m(),hadt.pt(),hadt.eta(),hadt.phi(),
            0,0,0,0,//hadW.m(),hadW.pt(),hadW.eta(),hadW.phi(),
            0,0};//truth,0}; ---> do we need truth?
          for(unsigned i=0;i<10;i++) basicHadtop[i]->Fill(hadtop[i],weight);  
          passcut=true; typeReco=1;
        } // 2btag 
    } // close analysis /
    //return passcut;
    return typeReco;
} //close fullylep    
/////////////////////////////////////////////////////////////////
int recol( vector<PseudoJet> jets,vector<PseudoJet> leptons,vector<PseudoJet> neutrinos, double weight){
    // from the jet collection find the hadronic W
    // construct met from all the rest
    // do a isolation vector
    vector<double> LepIso; int nlep=0; int nlepsurvive=0; //cout<<"njet "<<jets.size()<<endl;
    unsigned int jsize = jets.size();
    for(unsigned int j = 0;j<leptons.size();j++) {
        for(unsigned int i = 0;i<jsize;i++) LepIso.push_back(leptons.at(j).delta_R(jets.at(i)));
        if(leptons.size()>1 && jsize >0) for(unsigned int i = 0;i<leptons.size();i++) if (i!=j) LepIso.push_back(leptons.at(j).delta_R(leptons.at(i)));
        double MinDRLep = TMath::LocMin(LepIso.size(), &LepIso[0]);
        if(1>0 && leptons.size()>1 && jsize >0
           && LepIso[MinDRLep] > lepiso  
           && leptons.at(j).pt()> ptlepton 
           && abs(leptons.at(j).eta())< etal
           ) nlep++; // close basic cuts
    } // close for nlep
    //cout<<nlep<<endl;
    Nlep_passing_kLooseID->Fill(nlep,weight);
    return nlep;
}// close recow
/////////////////////////////////////////////////////////////////
int isbtagged(vector<PseudoJet> jets, vector<int> & btag, vector<int> & bmistag, vector<int> & btrue);
/////////////////////////////////////////////////////////////////
int recojets(vector<PseudoJet> particles,vector<PseudoJet> & jets_akt, vector<int> & btag, vector<int> & bmistag, vector<int> & fattag, vector<int> & btrue, double weight){
    JetDefinition akt(antikt_algorithm, RR);
    ClusterSequence cs_akt(particles, akt);
    //vector<PseudoJet> jets_akt;
    Selector jet_selector = SelectorPtMin(jet_ptmin) && SelectorAbsRapMax(rapmax);
    if(shower){jets_akt = sorted_by_pt(jet_selector(cs_akt.inclusive_jets()));} // first we do akt jets from particles
    else{double const ptmin=0.0; Selector jet_selector_parton = SelectorPtMin(ptmin);jets_akt = sorted_by_pt(jet_selector_parton(cs_akt.inclusive_jets()));}//if parton jet pt min zero
    //////////////////// do a jet cut
    vector<PseudoJet> jets_final;     unsigned int njets;
    for (unsigned int i = 0; i < jets_akt.size(); i++) if(jets_akt.at(i).pt()>jet_ptminfinal) jets_final.push_back(jets_akt.at(i));
    //for (unsigned int i = 0; i < particles.size(); i++) if(particles.at(i).pt()<jet_ptminfinal) {njets =0; return njets; } //jets_final.push_back(particles.at(i));
    njets = jets_final.size();
    int nbtag = isbtagged(jets_akt, btag, bmistag,btrue); // check wheather the b(c)jet is b--(mis)tagable
    // fill btags
    int count=0; for(unsigned int i = 0; i < njets; i++) if(btag[i]>0)count++; //btagselected->Fill(count,weight);
    ///////////////////// check tag
    JetDefinition CA10(cambridge_algorithm, Rsb);
    // Filter definition to improve mass resolution
    Filter filter(JetDefinition(cambridge_algorithm, Rfilt), SelectorNHardest(n_subjet));
    PseudoJet tagged_jet;
    for (unsigned int i = 0; i < njets; i++) { // to each akt jet
        // first recluster with some large CA (needed for mass-drop)
        ClusterSequence cs_tmp(jets_akt[i].constituents(), CA10);
        // next get hardest jet
        PseudoJet ca_jet = sorted_by_pt(cs_tmp.inclusive_jets())[0]; // find the cores
        // now run mass drop tagger
        MassDropTagger md_tagger(mu, ycut); // define the cut on mass drop
        // mu: ratio in between mass of cores, symetric splitting
        tagged_jet = md_tagger(ca_jet);
        if(tagged_jet.m()>10) {
            PseudoJet filtered_jet = filter(jets_akt.at(i)); // filter to tag
            if(filtered_jet.m()>Mfat) fattag.push_back(i); //see = 1; else see = 0; // no fat tag 
        } //else see = 0; 
    } // close find mass drop
    //jets = jets_akt;---
    Njets_passing_kLooseID->Fill(njets,weight);
    btagselected->Fill(nbtag,weight); 
    return njets;
} // close cluster jets
////////////////////////////////////////////////////////////////////////////////////////////////
int isbtagged(vector<PseudoJet> jets, vector<int> & btag, vector<int> & bmistag, vector<int> & btrue){ 
    unsigned int jsize = jets.size();
    int nbtag=0;
    for (unsigned int i=0; i<jsize; i++) { // check wheter jet have inside a b's are taggable 
        vector<PseudoJet> constitu=jets.at(i).constituents();
        unsigned int csize = constitu.size();
        int see=0,see2=0,seetruth=0,id;
        for (unsigned int j=0; j<csize; j++) {
            //cout<<"constituents flavour "<<constitu.at(j).user_index()<<endl;
            if(constitu.at(j).user_index() == 5 || constitu.at(j).user_index() == -5) {
                seetruth++; id=constitu.at(j).user_index();
                if(constitu.at(j).pt() > bjetpt && constitu.at(j).eta() < etab) {see++; nbtag++;}// to reco btag 
            }  
            //if( abs(constitu.at(j).user_index()) == 4  // work !!
            //     && constitu.at(j).pt() > bjetpt
            // && constitu.at(j).eta() < etab
            //  ) {see2++;}// bmistag.push_back(1);} bmistag.push_back(0);
        } // close constituents
        //bmistag.push_back(see2);
        //cout<<see<<endl;
        btag.push_back(see); //else btag.push_back(0); // count all tag/jet
        //if(see==0 && see2>0) bmistag.push_back(1); else bmistag.push_back(0); // count only one tag/jet
        if(seetruth>0 && see>0) btrue.push_back(id); else btrue.push_back(0); // to mctruth
        //cout<<"b-quarks mistagged = " <<bmistag[i] <<" b-quark = " <<btag[i] <<endl;
    } // close for each jet
    //int numbb=0; for(unsigned i=0 ; i< btrue.size() ; i++ ) numbb += btrue[i];
    //cout<<" total "<<numbb<<endl;
    return nbtag;
} // close isbtagged
/////////////////////////////////////////////////////////////////////////
// save the histos
int save_hist(int isample,int reco,int sample, bool shower){
    const char* Mass;
    Mass = Form("Control_mtdef_%d_type_%d_Gamm_%d_shower_%d_minlj.root",reco,sample,isample,shower); 
    cout<<" saved "<< Form("Control_mtdef_%d_type_%d_Gamm_%d_shower_%d_minlj.root",reco,sample,isample,shower)<<endl;
    TFile f1(Mass, "recreate");
    f1.cd();
    Njets_passing_kLooseID->Write();
    Nlep_passing_kLooseID->Write();
    btagselected->Write();
    genmbb->Write();
    leptop->Write();
    hadtop->Write();
    genbhad->Write();
    genblep->Write();
    genbhadeta->Write();
    genblepeta->Write();
    for(unsigned i=0;i<14;i++) basicLeptons[i]->Write();
    for(unsigned i=0;i<10;i++) basicHadtop[i]->Write();
    for(unsigned i=0;i<19;i++) basicLeptop[i]->Write();
    f1.Close();
    Njets_passing_kLooseID->Reset();
    Nlep_passing_kLooseID->Reset();
    btagselected->Reset();
    genmbb->Reset();
    leptop->Reset();
    hadtop->Reset();
    genbhad->Reset();
    genblep->Reset();
    genbhadeta->Reset();
    genblepeta->Reset();
    basicLeptons.clear();
    basicHadtop.clear();
    basicLeptop.clear();
    //  basicLeptons[0]->Reset();
    //  basicLeptons[1]->Reset();
    //  basicLeptons[2]->Reset();
    //  for(unsigned i=0;i<10;i++) basicHadtop[i]->Reset();
    //  for(unsigned i=0;i<13;i++) basicLeptop[i]->Reset();

    return 0;
}
/////////////////////////////////////////////////////////////////////////
int decla(int mass){
    
    delete gDirectory->FindObject("leptop1");
    delete gDirectory->FindObject("hadtop1");
    delete gDirectory->FindObject("njets_passing_kLooseID_ct4");
    delete gDirectory->FindObject("btagselected");
    delete gDirectory->FindObject("E1histpt");
    delete gDirectory->FindObject("E1histeta");
    delete gDirectory->FindObject("MetMass_ct4");
    delete gDirectory->FindObject("H1hist");
    delete gDirectory->FindObject("H1histpt");
    delete gDirectory->FindObject("H1histeta");
    delete gDirectory->FindObject("H1histphi");
    delete gDirectory->FindObject("HW1hist");
    delete gDirectory->FindObject("HW1histpt");
    delete gDirectory->FindObject("HW1histeta");
    delete gDirectory->FindObject("HW1histphi");
    delete gDirectory->FindObject("recotruth");
    delete gDirectory->FindObject("detabb");
    delete gDirectory->FindObject("H1LepThist");
    delete gDirectory->FindObject("H1LepThistpt");
    delete gDirectory->FindObject("H1LepThisteta");
    delete gDirectory->FindObject("H1LepThistphi");
    delete gDirectory->FindObject("HW1LepThist");
    delete gDirectory->FindObject("HW1LepThistpt");
    delete gDirectory->FindObject("HW1LepThisteta");
    delete gDirectory->FindObject("HW1LepThistphi");
    delete gDirectory->FindObject("pnuzerror");
    delete gDirectory->FindObject("recotruthlept");
    delete gDirectory->FindObject("mterror");
    delete gDirectory->FindObject("wmt");
    delete gDirectory->FindObject("tmt");
    delete gDirectory->FindObject("detalb");
    delete gDirectory->FindObject("mraz");
    delete gDirectory->FindObject("mrazt");
    delete gDirectory->FindObject("razratio");
    
    
    const char* label="without btag im reco";
    
    Njets_passing_kLooseID = new TH1D("njets_passing_kLooseID_ct4",  
                                      label, 
                                      13, -0.5, 12.5);
    Njets_passing_kLooseID->GetYaxis()->SetTitle("");
    Njets_passing_kLooseID->GetXaxis()->SetTitle("Njets after showering"); 

    Nlep_passing_kLooseID = new TH1D("nlep_passing_kLooseID_ct4",  
                                      label, 
                                      13, -0.5, 12.5);
    Nlep_passing_kLooseID->GetYaxis()->SetTitle("");
    Nlep_passing_kLooseID->GetXaxis()->SetTitle("Nleptons after showering"); 
    
    btagselected = new TH1D("btagselected",  
                            label, 
                            13, -0.5, 12.5);
    btagselected->GetYaxis()->SetTitle("");
    btagselected->GetXaxis()->SetTitle("b-tagable b's on selected events");
    
    genmbb = new TH1D("genmbb",  
                      label, 
                      70, 0, 1000);
    genmbb->GetYaxis()->SetTitle("");
    genmbb->GetXaxis()->SetTitle("true Mbb"); 
    
    leptop = new TH1D("leptop1",  
                      label, 
                      170, 0, 500);
    //leptop->SetLogY(1);    
    leptop->GetYaxis()->SetTitle("");
    leptop->GetXaxis()->SetTitle("true M lep top"); 
    
    hadtop = new TH1D("hadtop1",  
                      label, 
                      140, 0, 500);
    hadtop->GetYaxis()->SetTitle("");
    hadtop->GetXaxis()->SetTitle("true M had top"); 
    
    genbhad = new TH1D("genbhad",
                       label,
                       70, 0, 200);
    genbhad->GetYaxis()->SetTitle("");
    genbhad->GetXaxis()->SetTitle("true b pt had top");
    
    genblep = new TH1D("genblep",
                       label,
                       70, 0, 200);
    genblep->GetYaxis()->SetTitle("");
    genblep->GetXaxis()->SetTitle("true b pt lep top");
    
    
    genbhadeta = new TH1D("genbhadeta",
                          label,
                          20, -5, 5);
    genbhadeta->GetYaxis()->SetTitle("");
    genbhadeta->GetXaxis()->SetTitle("true b eta had top");
    
    genblepeta = new TH1D("genblepeta",
                          label,
                          20, -5, 5);
    genblepeta->GetYaxis()->SetTitle("");
    genblepeta->GetXaxis()->SetTitle("true b eta lep top");
    // gen 
    /////////////////////////////////////////////////////////////////////////////
    // for leptons
    TH1D *E1histpt = new TH1D("E1histpt",  
                              label, 
                              50, 0, 200);
    E1histpt->GetYaxis()->SetTitle("Events/ 2 GeV");
    E1histpt->GetXaxis()->SetTitle("lepton 1 P_T (GeV)");
    basicLeptons.push_back (E1histpt); 
    
    TH1D *E1histeta = new TH1D("E1histeta",  
                               label, 
                               30, -6, 6);
    E1histeta->GetYaxis()->SetTitle("Events/ 2 GeV");
    E1histeta->GetXaxis()->SetTitle("#eta_{lepton1} (GeV)");
    basicLeptons.push_back (E1histeta); 

    TH1D *E2histpt = new TH1D("E2histpt",  
                              label, 
                              50, 0, 200);
    E2histpt->GetYaxis()->SetTitle("Events/ 2 GeV");
    E2histpt->GetXaxis()->SetTitle("lepton 2 P_T (GeV)");
    basicLeptons.push_back (E2histpt); 
    
    TH1D *E2histeta = new TH1D("E2histeta",  
                               label, 
                               30, -6, 6);
    E2histeta->GetYaxis()->SetTitle("Events/ 2 GeV");
    E2histeta->GetXaxis()->SetTitle("#eta_{lepton2} (GeV)");
    basicLeptons.push_back (E2histeta); 
    
    TH1D *MetMass = new TH1D("MetMass_ct4",  
                             label, 
                             50, 0, 700);
    MetMass->GetXaxis()->SetTitle("MET (GeV)");
    basicLeptons.push_back (MetMass);

    //TH1D *DetaLep = new TH1D("DetaLep_ct4",  
    //                         label, 
    //                         50, -6, 6);
    //DetaLep->GetXaxis()->SetTitle("#Delta #eta (ll)");
    //basicLeptons.push_back (DetaLep);

    TH1D *llMass = new TH1D("llMass",  
                             label, 
                             50, 0, 700);
    llMass->GetXaxis()->SetTitle("M_{ll} (GeV)");
    basicLeptons.push_back (llMass);    

    TH1D *jjMass = new TH1D("jjMass",  
                            label, 
                            50, 0, 700);
    jjMass->GetXaxis()->SetTitle("M_{jj} (GeV)");
    basicLeptons.push_back (jjMass);    
    
    TH1D *jbleadMass = new TH1D("jbleadMass",  
                            label, 
                            50, 150, 370);
    jbleadMass->GetXaxis()->SetTitle("M_{l-j} (GeV) m(j l+) min");
    basicLeptons.push_back (jbleadMass);    
    
    TH1D *jbsubleadMass = new TH1D("jbsubleadMass",  
                                label, 
                                50, 150, 370);    
    jbsubleadMass->GetXaxis()->SetTitle("M_{l-j} (GeV) m(j l-) min");
    basicLeptons.push_back (jbsubleadMass);    
    
    TH1D *jbleadXMass = new TH1D("jbleadXMass",  
                                 label, 
                                 50, 330, 200);
    jbleadXMass->GetXaxis()->SetTitle("M_{l-j} (GeV) leading Balance");
    basicLeptons.push_back (jbleadXMass);    
    
    TH1D *jbsubleadXMass = new TH1D("jbsubleadXMass",  
                                    label, 
                                    50, 0, 700);
    jbsubleadXMass->GetXaxis()->SetTitle("M_{l-j} (GeV) subleading Balance");
    basicLeptons.push_back (jbsubleadXMass);      

    TH1D *truthMBplot1 = new TH1D("truthMBplot1",  
                                 label, 
                                 7, -1.5, 5.5);
    truthMBplot1->GetXaxis()->SetTitle("1 lepton-jet category 1D");
    basicLeptons.push_back (truthMBplot1);      

    
    TH1D *truthMBplot = new TH1D("truthMBplot",  
                                    label, 
                                 7, -1.5, 5.5);
    truthMBplot->GetXaxis()->SetTitle("1 lepton-jet category 1D");
    basicLeptons.push_back (truthMBplot);      

    TH1D *typetruthMBplot = new TH1D("typetruthMBplot",  
                                 label, 
                                 7, -1.5, 5.5);
    typetruthMBplot->GetYaxis()->SetTitle("% from selected events");
    typetruthMBplot->GetXaxis()->SetTitle("lepton-jet category 2D");
    basicLeptons.push_back (typetruthMBplot);      
    
    ///////////////////////////////////////////////////////////////////////////
    // for hadronic tops
    TH1D *H1hist = new TH1D("H1hist",  
                            label, 
                            70, 0, 1000);
    H1hist->GetYaxis()->SetTitle("Events/ 2 GeV");
    H1hist->GetXaxis()->SetTitle("mass t_{had} (GeV)");
    basicHadtop.push_back (H1hist); 
    
    TH1D *H1histpt = new TH1D("H1histpt",  
                              label, 
                              100, 0, 1000);
    H1histpt->GetYaxis()->SetTitle("Events/ 2 GeV");
    H1histpt->GetXaxis()->SetTitle("t_{had} P_T (GeV)");
    basicHadtop.push_back (H1histpt); 
    
    TH1D *H1histeta = new TH1D("H1histeta",  
                               label, 
                               30, -6, 6);
    H1histeta->GetYaxis()->SetTitle("Events/ 2 GeV");
    H1histeta->GetXaxis()->SetTitle("#eta_{t_{had}} (GeV)");
    basicHadtop.push_back (H1histeta); 
    
    TH1D *H1histphi = new TH1D("H1histphi",  
                               label, 
                               30, 0, 5);
    H1histphi->GetYaxis()->SetTitle("Events/ 2 GeV");
    H1histphi->GetXaxis()->SetTitle("#phi_{t_{had}} (GeV)");
    basicHadtop.push_back (H1histphi); 
    // had w
    TH1D *HW1hist = new TH1D("HW1hist",  
                             label, 
                             80, 70, 90);
    HW1hist->GetYaxis()->SetTitle("Events/ 2 GeV");
    HW1hist->GetXaxis()->SetTitle("mass W_{had} (GeV)");
    basicHadtop.push_back (HW1hist); 
    
    TH1D *HW1histpt = new TH1D("HW1histpt",  
                               label, 
                               100, 0, 600);
    HW1histpt->GetYaxis()->SetTitle("Events/ 2 GeV");
    HW1histpt->GetXaxis()->SetTitle("W_{had} P_T (GeV)");
    basicHadtop.push_back (HW1histpt); 
    
    TH1D *HW1histeta = new TH1D("HW1histeta",  
                                label, 
                                30, -6, 6);
    HW1histeta->GetYaxis()->SetTitle("Events/ 2 GeV");
    HW1histeta->GetXaxis()->SetTitle("#eta_{W_{had}} (GeV)");
    basicHadtop.push_back (HW1histeta); 
    
    TH1D *HW1histphi = new TH1D("HW1histphi",  
                                label, 
                                30, 0, 5);
    HW1histphi->GetYaxis()->SetTitle("Events/ 2 GeV");
    HW1histphi->GetXaxis()->SetTitle("#phi_{W_{had}} (GeV)");
    basicHadtop.push_back (HW1histphi); 
    
    TH1D *recotruth = new TH1D("recotruth",  
                               label, 
                               3, -0.5, 2.5);
    recotruth->GetYaxis()->SetTitle("Events/ 2 GeV");
    recotruth->GetXaxis()->SetTitle("reco truth");
    basicHadtop.push_back (recotruth); 
    
    TH1D *detabb = new TH1D("detabb",  
                            label, 
                            30, 0., 5);
    detabb->GetYaxis()->SetTitle("Events/ 2 GeV");
    detabb->GetXaxis()->SetTitle("#Delta#eta bb");
    basicHadtop.push_back (detabb); 
    ///////////////////////////////////////////////////////////////////////////
    // for leptonic tops
    TH1D *H1LepThist = new TH1D("H1LepThist",  
                                label, 
                                70, 0, 1000);
    H1LepThist->GetYaxis()->SetTitle("Events/ 2 GeV");
    H1LepThist->GetXaxis()->SetTitle("mass t_{lep} (GeV)");
    basicLeptop.push_back (H1LepThist); 
    
    TH1D *H1LepThistpt = new TH1D("H1LepThistpt",  
                                  label, 
                                  20, 0, 300);
    H1LepThistpt->GetYaxis()->SetTitle("Events/ 2 GeV");
    H1LepThistpt->GetXaxis()->SetTitle("t_{lep} P_T (GeV)");
    basicLeptop.push_back (H1LepThistpt); 
    
    TH1D *H1LepThisteta = new TH1D("H1LepThisteta",  
                                   label, 
                                   30, -6, 6);
    H1LepThisteta->GetYaxis()->SetTitle("Events/ 2 GeV");
    H1LepThisteta->GetXaxis()->SetTitle("#eta_{t_{lep}} (GeV)");
    basicLeptop.push_back (H1LepThisteta); 
    
    TH1D *H1LepThistphi = new TH1D("H1LepThistphi",  
                                   label, 
                                   30, 0, 5);
    H1LepThistphi->GetYaxis()->SetTitle("Events/ 2 GeV");
    H1LepThistphi->GetXaxis()->SetTitle("#phi_{t_{lep}} (GeV)");
    basicLeptop.push_back (H1LepThistphi); 
    // had w
    TH1D *HW1LepThist = new TH1D("HW1LepThist",  
                                 label, 
                                 100, 80, 81);
    HW1LepThist->GetYaxis()->SetTitle("Events/ 2 GeV");
    HW1LepThist->GetXaxis()->SetTitle("mass W_{lep} (GeV)");
    basicLeptop.push_back (HW1LepThist); 
    
    TH1D *HW1LepThistpt = new TH1D("HW1LepThistpt",  
                                   label, 
                                   20, 0, 300);
    HW1LepThistpt->GetYaxis()->SetTitle("Events/ 2 GeV");
    HW1LepThistpt->GetXaxis()->SetTitle("W_{lep} P_{T} (GeV)"); 
    basicLeptop.push_back (HW1LepThistpt); 
    
    TH1D *HW1LepThisteta = new TH1D("HW1LepThisteta",  
                                    label, 
                                    30, -6, 6);
    HW1LepThisteta->GetYaxis()->SetTitle("Events/ 2 GeV");
    HW1LepThisteta->GetXaxis()->SetTitle("#eta_{W_{lep}} (GeV)");
    basicLeptop.push_back (HW1LepThisteta); 
    
    TH1D *HW1LepThistphi = new TH1D("HW1LepThistphi",  
                                    label, 
                                    30, 0, 5);
    HW1LepThistphi->GetYaxis()->SetTitle("Events/ 2 GeV");
    HW1LepThistphi->GetXaxis()->SetTitle("#phi_{W_{lep}} (GeV)");
    basicLeptop.push_back (HW1LepThistphi); 
    
    TH1D *pnuzerror = new TH1D("pnuzerror",  
                               label, 
                               120, -60, 60);
    pnuzerror->GetYaxis()->SetTitle("Events/ 2 GeV");
    pnuzerror->GetXaxis()->SetTitle("mW(reco) - mW(truth)");
    basicLeptop.push_back (pnuzerror); 
    
    TH1D *recotruthlept = new TH1D("recotruthlept",  
                                   label, 
                                   3, -0.5, 2.5);
    recotruthlept->GetYaxis()->SetTitle("Events/ 2 GeV");
    recotruthlept->GetXaxis()->SetTitle("reco truth lept");
    basicLeptop.push_back (recotruthlept); 
    
    TH1D *mterror = new TH1D("mterror",  
                             label, 
                             120, -60, 60);
    mterror->GetYaxis()->SetTitle("Events/ 2 GeV");
    mterror->GetXaxis()->SetTitle("mt(reco) - mt(truth)");
    basicLeptop.push_back (mterror); 
    
    TH1D *wmt = new TH1D("wmt",  
                         label, 
                         75, 0, 150);
    wmt->GetYaxis()->SetTitle("Events/ 2 GeV");
    wmt->GetXaxis()->SetTitle("W transverse mass");
    basicLeptop.push_back (wmt); 
    
    TH1D *tmt = new TH1D("tmt",  
                         label, 
                         50, 0, 600);
    tmt->GetYaxis()->SetTitle("Events/ 2 GeV");
    tmt->GetXaxis()->SetTitle("Top transverse mass (Wb)");
    basicLeptop.push_back (tmt); 
    
    TH1D *detalb = new TH1D("detalb",  
                            label, 
                            30, 0, 6);
    detalb->GetYaxis()->SetTitle("Events/ 2 GeV");
    detalb->GetXaxis()->SetTitle("#Delta#eta b l");
    basicLeptop.push_back (detalb); 
    
    TH1D *mraz = new TH1D("mraz",  
                          label, 
                          150, 0, 1000);
    mraz->GetYaxis()->SetTitle("Events/ 2 GeV");
    mraz->GetXaxis()->SetTitle("MR");
    basicLeptop.push_back (mraz); 
    
    TH1D *mrazt = new TH1D("mrazt",  
                           label, 
                           150, 0, 300);
    mrazt->GetYaxis()->SetTitle("Events/ 2 GeV");
    mrazt->GetXaxis()->SetTitle("MRt");
    basicLeptop.push_back (mrazt); 
    
    TH1D *razratio = new TH1D("razratio",  
                              label, 
                              100, 0, 1);
    razratio->GetYaxis()->SetTitle("Events/ 2 GeV");
    razratio->GetXaxis()->SetTitle("MRt/MR");
    basicLeptop.push_back (razratio); 
    
    TH1D *tmtal = new TH1D("tmtal",  
                           label, 
                           50, 0, 600);
    tmtal->GetYaxis()->SetTitle("Events/ 2 GeV");
    tmtal->GetXaxis()->SetTitle("top transverse mass (lb - #nu)");
    basicLeptop.push_back (tmtal); 
    
    TH1D *mbl = new TH1D("massbl",  
                         label, 
                         30, 0, 1000);
    mbl->GetYaxis()->SetTitle("Events/ 2 GeV");
    mbl->GetXaxis()->SetTitle("m_{lb}");
    basicLeptop.push_back (mbl); 
    ///////////////////////////////////////////////////////////////////////////////////
    return 0;
}


/*
 */ 
/*  /////////////////////////////////////////////////////////////////////
 // find w -- closest hadronic mass
 vector<int> j1,j2; vector<double> a3; int minM;
 for(unsigned int nj1=0; nj1< jsize; nj1++) 
 for(unsigned int nj2=nj1+1; nj2< jsize; nj2++) 
 if(btag[nj1]==0 && btag[nj2]==0)
 { 
 //std::cout<<nj1<<nj2<<" "<<nj3<<nj4<<std::endl; 
 double invmassA =  (jets.at(nj1)+jets.at(nj2)).m();
 a3.push_back((invmassA-wmass)*(invmassA-wmass)); 
 j1.push_back(nj1); j2.push_back(nj2); 
 // we also what to keep the nj...           
 } // loop on jets  
 minM = TMath::LocMin(a3.size(), &a3[0]);
 PseudoJet hadW = jets[j1[minM]] + jets[j2[minM]];
 //cout<<" the hadronic w's are: "<<j1[minM]<<" "<<j2[minM]<<" "<<btrue[j1[minM]]<<" "<<btrue[j2[minM]]<<endl;
 // find t -- closest jet
 unsigned wjet1=j1[minM], wjet2=j2[minM]; 
 vector<int> j3; vector<double> a4; int minMt;
 for(unsigned int nj1=0; nj1< jsize; nj1++) 
 if(nj1!=wjet1 && nj1!=wjet2) 
 if(btag[nj1]>0) 
 { 
 double invmassA =  (jets.at(nj1)+hadW).m();
 a4.push_back((invmassA-tmass)*(invmassA-tmass)); j3.push_back(nj1); 
 //double dr=  jets.at(nj1).delta_R(hadW);
 //a4.push_back(dr); j3.push_back(nj1);
 //cout<<jsize<<" "<<btag.size()<<" "<<btag[nj1]<<endl;
 } // loop on jets
 //if(a4.size()>0){
 minMt = TMath::LocMin(a4.size(), &a4[0]); bh = j3[minMt];
 //cout<<" the hadronic t's are: "<<j1[minM]<<" "<<j2[minM]<<" "<<j3[minMt]<<endl;
 PseudoJet hadt = jets[bh] + hadW;
 */
/////////////////////////////////////////////////
//TVector2 teste = jets.at(j1[minM]).pt()+jets.at(j2[minM]).pt();
//cout<<" teste "<<
//(jets.at(j1[minM]).px()+jets.at(j2[minM]).px())*(jets.at(j1[minM]).px()+jets.at(j2[minM]).px())+
//(jets.at(j1[minM]).py()+jets.at(j2[minM]).py())*(jets.at(j1[minM]).py()+jets.at(j2[minM]).py())<<" "<<
//(jets.at(j1[minM]).pt()+jets.at(j2[minM]).pt())*(jets.at(j1[minM]).pt()+jets.at(j2[minM]).pt())<<" "<<
// teste*teste<<" "<<endl;
