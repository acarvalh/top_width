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
int truetops(vector<PseudoJet> jets, vector<PseudoJet> leptons,vector<PseudoJet> neutrinos, vector<int> btag, vector<int> btrue){
    int bl=-1,bh=-1; vector<int> wj;
    unsigned int jsize = jets.size();
    for(unsigned int nj1=0; nj1< jsize; nj1++) // only works at parton level
        if(btrue[nj1]==-5) bh=nj1; else if(btrue[nj1]==5) bl=nj1; else wj.push_back(nj1);   
    //cout<<wj.size()<<" "<<jsize<<" "<<bl<<" "<<bh<<" "<<endl; 
    int sample;
    if(bh!=-1 && bl!=-1 && wj.size()>1){
        //   cout<<jsize<<endl;
        //cout<<wj[0]<<" "<<wj[1]<<" "<<bl<<" "<<bh<<" "<<endl;  
        PseudoJet lepTtrue = leptons.at(0) + neutrinos.at(0) + jets.at(bl);
        PseudoJet hadTtrue = jets.at(wj[0]) + jets.at(wj[1]) + jets.at(bh);
        // separate samples
        // leptonic 0 = on shell (mass within mtop +- 0.5 GeV), 1 = off down, 2 = off up
        // hadronic 3 = on shell (mass within mtop +- 0.5 GeV), 4 = off down, 5 = off up
        //if( hadTtrue.m()<tmass+0.5 && hadTtrue.m()>tmass-0.5) sample =0;
        //if( hadTtrue.m()<tmass-0.5 ) sample =1;
        //if( hadTtrue.m()>tmass+0.5 ) sample =2;
        //
        if( lepTtrue.m()<tmass+27 && lepTtrue.m()>=tmass-17) sample =3;
        else if( lepTtrue.m()<tmass-17 ) sample =4;
        else if( lepTtrue.m()>=tmass+27 ) sample =5;
    }
    //return 0;
    return sample;
} // close truetops 
/////////////////////////////////////////////////////////////////
bool fullylep(int & bh, int & bl,vector<PseudoJet> jets, vector<PseudoJet> leptons,vector<PseudoJet> neutrinos, vector<int> btag, vector<int> btrue, double met, double weight, double cut, int type){
    ///////////////////////////////////////////////////////////////////
    // gen level info // had == plus
    int blll=-1, bhhh=-1, lep1=-1, lep2=-1, nu1=-1, nu2=-1;
    unsigned int jsize = jets.size();
    for(unsigned int nj1=0; nj1< jsize; nj1++) if(btrue[nj1]==5) bhhh=nj1; else if(btrue[nj1]==-5) blll=nj1; 
    int counttruth =-1; 
    // I have 2 isolated leptons, I assume are the two leading
    for(unsigned int nj1=0; nj1< 2; nj1++) {if(leptons.at(nj1).user_index() > 0) lep1=nj1; else lep2=nj1;  }
    for(unsigned int nj1=0; nj1< 2; nj1++) {if(neutrinos.at(nj1).user_index() > 0) nu1=nj1; else  nu2=nj1;  }
    //cout<<wj.size()<<" "<<jsize<<" "<<blll<<" "<<bhhh<<" "<<endl;  
    PseudoJet lepTtrue, hadTtrue;
    if(blll!=-1 && bhhh!=-1 && lep1!=-1 && lep2!=-1 && nu1!=-1 && nu2!=-1){ 
        // cout<<jsize<<endl;
        //cout<<wj[0]<<" "<<wj[1]<<" "<<blll<<" "<<bhhh<<" "<<endl;  
        lepTtrue = leptons.at(lep1) + neutrinos.at(nu2) + jets.at(blll);
        hadTtrue = leptons.at(lep2) + neutrinos.at(nu1) + jets.at(bhhh);
        //
        //if(ifolder==0) counttruth==1; // full
        //else if(hadTtrue.m() < genmasshad && hadTtrue.m() > genmasshadmin && lepTtrue.m() < genmasslep && lepTtrue.m() > genmasslepmin) counttruth=1; // onon
        //else if(ifolder==2 && hadTtrue.m() < genmass && lepTtrue.m() > genmass) counttruth=1; // onoff
        //else if(ifolder==3 && hadTtrue.m() > genmass && lepTtrue.m() < genmass) counttruth=1; // offon
        //else if(ifolder==4 && hadTtrue.m() > genmass && lepTtrue.m() > genmass) counttruth=1; // offoff
        
        // }
        //int truth=-10;
        // fill if: 0 = (< mt1,mt2) | 1 = (< m1 , >m2) | 2 = (>m1 , m2<) | 3 = (m1,m2 >)
        if(hadTtrue.m()<cut && lepTtrue.m() < cut) counttruth=0;
          else if(hadTtrue.m()<cut && lepTtrue.m()>=cut) counttruth=1;
          else if(hadTtrue.m()>=cut && lepTtrue.m()<cut) counttruth=2;            
          else if(hadTtrue.m()>=cut && lepTtrue.m()>=cut) counttruth=3;
        //cout<<"here 3 "<<cut<<endl;
        //if (counttruth==2 || counttruth==3)cout<<"here 3 "<<counttruth<<endl;
        //
        if(1>0
           && counttruth==type
           ) {    
        //truth=1;
            //cout<<jsize<<endl;
            //cout<<(jets.at(bhhh)+jets.at(blll)).m()<<endl;
            //cout<<counttruth<<endl;
            //
        Njets_passing_kLooseID->Fill(jsize,weight);
        genmbb->Fill((jets.at(bhhh)+jets.at(blll)).m(),weight);        
        leptop->Fill(lepTtrue.m(),weight);
        hadtop->Fill(hadTtrue.m(),weight);
        genbhad->Fill(jets.at(bhhh).pt(),weight);
        genblep->Fill(jets.at(blll).pt(),weight);
        genbhadeta->Fill(jets.at(bhhh).eta(),weight);
        genblepeta->Fill(jets.at(blll).eta(),weight);
            return true;
        } else return false; // close if cut
    } else return false;
    //////////////////////////////  
    //cout<<cut<<" "<<type<<endl;
    //double hadtop[10]={1,1,1,1, //hadt.m(),hadt.pt(),hadt.eta(),hadt.phi(),
    //    1,1,1,1, //hadW.m(),hadW.pt(),hadW.eta(),hadW.phi(),
    //    1,1}; //truth,detabb};
    //for(unsigned i=0;i<10;i++) basicHadtop[i]->Fill(hadtop[i],weight);
    //cout<<cut<<" "<<type<<endl;
    
} //close fullylep    
/////////////////////////////////////////////////////////////////
bool recolept2step(int & bh, int & bl,vector<PseudoJet> jets, vector<PseudoJet> leptons,vector<PseudoJet> neutrinos, vector<int> btag, vector<int> btrue, double met, double weight, double cut, int type){
    // I did not reco the hadronic --- do not know which b to take
    // first try in the 2-step way 
    //double mw = (leptons.at(0)+neutrinos.at(0)).m();//teste
    //cout<<"met "<<met<<" pzl "<<neutrinos.at(0).pt()<<endl;
    int blll=-1,bhhh=-1; vector<int> wj;
    unsigned int jsize = jets.size();
    for(unsigned int nj1=0; nj1< jsize; nj1++) // only works at parton level
        if(btrue[nj1]==-5) bhhh=nj1; else if(btrue[nj1]==5) blll=nj1; else wj.push_back(nj1);   
    //cout<<wj.size()<<" "<<jsize<<" "<<bl<<" "<<bh<<" "<<endl; 
    if(bhhh!=-1 && blll!=-1 && wj.size()>1){
        //   cout<<jsize<<endl;
        //cout<<wj[0]<<" "<<wj[1]<<" "<<bl<<" "<<bh<<" "<<endl;  
        PseudoJet lepTtrue = leptons.at(0) + neutrinos.at(0) + jets.at(blll);
        PseudoJet hadTtrue = jets.at(wj[0]) + jets.at(wj[1]) + jets.at(bhhh);
        leptop->Fill(lepTtrue.m(),weight);
        hadtop->Fill(hadTtrue.m(),weight);
    }
    /////////////////////////////////////////////////////////////
    double wt = (leptons.at(0).px()*neutrinos.at(0).px()) + (leptons.at(0).py()*neutrinos.at(0).py());
    double mu = (pow(wmass,2)/2 + wt);
    double aw = (pow(leptons.at(0).pz(),2)-pow(leptons.at(0).e(),2));
    double bw = 2*mu*leptons.at(0).pz();
    double cw = pow(mu,2)-pow(leptons.at(0).e(),2)*pow(met,2);
    double discriminant = pow(bw,2) - 4*aw*cw ; 
    double pnuzerror,pnuz;//,recowm,recowpt,recoweta,recowphi;  
    PseudoJet lepW; int recotruth;
    if(discriminant>=0){
        //recotruth=1; 
        double pznu1 = (-bw - sqrt(discriminant))/(2*aw);
        double pznu2 = (-bw + sqrt(discriminant))/(2*aw);
        //cout<<"pz1 "<<pznu1 << " pz2 " <<pznu2 <<" pz "<<neutrinos.at(0).pz()<<" pzl "<<leptons.at(0).pz()<<endl;
        //cout<<"dumb "<< (neutrinos.at(0)+leptons.at(0)).m() << " calculated " <<
        //sqrt(pow(leptons.at(0).e()+neutrinos.at(0).e(),2)-wt-pow(leptons.at(0).pz()+neutrinos.at(0).pz(),2))
        //  <<endl;
        pnuz=TMath::Min(pznu1,pznu2); 
        double enu = sqrt(pow(neutrinos.at(0).px(),2)+pow(neutrinos.at(0).py(),2)+pow(pnuz,2));
        double pxnu=neutrinos.at(0).px(),pynu=neutrinos.at(0).py();
        lepW = leptons.at(0)+ fastjet::PseudoJet(pxnu,pynu,pnuz,enu);
        pnuzerror=(neutrinos.at(0)+leptons.at(0)).m()-lepW.m();
        ///////////////////////////////////////////////////////////////
        // choose the other jet with the top mass  
        unsigned int jsize = jets.size();
        vector<int> j3; vector<double> a4; int minMt;
        for(unsigned int nj1=0; nj1< jsize; nj1++) {
            //if(btag[nj1]>0) 
            double invmassA =  (jets.at(nj1)+lepW).m();
            a4.push_back((invmassA-tmass)*(invmassA-tmass)); j3.push_back(nj1); 
        } // loop on jets
        minMt = TMath::LocMin(a4.size(), &a4[0]); bl = j3[minMt];   
        PseudoJet lept=lepW+jets[bl]; 
        //////////////////////////////////////////////////////////////
        PseudoJet lepWtransv = leptons.at(0) + jets.at(bl);
        double wtransM = // w transverse mass
        sqrt(
             2*neutrinos.at(0).e()*leptons.at(0).e()-
             2*(neutrinos.at(0).px()*leptons.at(0).px()+neutrinos.at(0).py()*leptons.at(0).py())); // fix neutrino
        double ttransM = // t transverse mass
        sqrt(
             2*neutrinos.at(0).e()*lepWtransv.e()-
             2*(neutrinos.at(0).px()*lepWtransv.px()+neutrinos.at(0).py()*lepWtransv.py())); // fix neutrino
        /////////////////////////////////////////////////////////////////// 
        // find hadronic w -- closest hadronic mass -- but not the hadronic top
        vector<int> j1,j2; vector<double> a3; int minM;
        for(unsigned int nj1=0; nj1< jsize; nj1++) 
            for(unsigned int nj2=nj1+1; nj2< jsize; nj2++) 
                if (nj1!=nj2 && nj1!=bl && nj2!=bl )
                    //if(btag[nj1]==0 && btag[nj2]==0)
                {  //std::cout<<nj1<<" "<<nj2<<" "<<bl<<" jsize "<<jsize<<std::endl; 
                    double invmassA =  (jets.at(nj1)+jets.at(nj2)).m();
                    a3.push_back((invmassA-wmass)*(invmassA-wmass)); 
                    j1.push_back(nj1); j2.push_back(nj2); 
                    // we also what to keep the nj...           
                } // loop on jets  
        minM = TMath::LocMin(a3.size(), &a3[0]);
        PseudoJet hadW = jets[j1[minM]] + jets[j2[minM]];
        ////////////////////////////////////////////////////////////
        // do the hadronic top with the hardest pt jet
        vector<double> ptj; vector<int> nptj;
        for(unsigned int nj1=0; nj1< jsize; nj1++) 
            if (nj1!=j1[minM] && nj1!=j2[minM] && nj1!=bl)
            { ptj.push_back(jets.at(nj1).pt()); nptj.push_back(nj1);  }
        int maxPt = TMath::LocMax(ptj.size(), &ptj[0]); 
        PseudoJet hadt = hadW + jets.at(nptj[maxPt]);
        ////////////////////////////////////////////////////////////////////
        // test truth 
        int truth; 
        if(btrue[j1[minM]]==0 && btrue[j2[minM]]==0  && btrue[bl]==-5 && btrue[nptj[maxPt]]==5 ) truth=1; else truth=0; 
        ///////////////////////////////////////////////////////////////////
        double detabb  = abs(jets[nptj[maxPt]].eta() - jets[bl].eta());
        /////////////////////////////////////////////////////////////////// 
        // fill diagrams
        double leptop[13]={lept.m(),lept.pt(),lept.eta(),lept.phi(),
            lepW.m(),lepW.pt(),lepW.eta(),lepW.phi(),
            pnuzerror,0,0,wtransM,ttransM}; // pnuzerror,truth,mterror,wmt,tmt
        for(unsigned i=0;i<13;i++) basicLeptop[i]->Fill(leptop[i],weight);
        basicLeptons[0]->Fill(leptons[0].pt(),weight);
        basicLeptons[1]->Fill(leptons[0].eta(),weight);
        basicLeptons[2]->Fill(met,weight);
        double hadtop[10]={hadt.m(),hadt.pt(),hadt.eta(),hadt.phi(),
            hadW.m(),hadW.pt(),hadW.eta(),hadW.phi(),
            truth,detabb};
        for(unsigned i=0;i<10;i++) basicHadtop[i]->Fill(hadtop[i],weight);
        //cout<<" the hadronic w's are: "<<lept.m()<<" "<<hadt.m()<<" "<<lepW.phi()<<" "<<lepW.pt()<<" "<<lepW.phi()<<endl;
        /*  double mtop=lept.m(); basicLeptop[0]->Fill(mtop,weight);
         double leptop[13]={170,0,0,0, //lept.m(),lept.pt(),lept.eta(),lept.phi(),
         90,0,0,0, //lepW.m(),lepW.pt(),lepW.eta(),lepW.phi(),
         0,0,0,wtransM,ttransM //pnuzerror,truth,mterror,wmt,tmt};
         };
         for(unsigned i=1;i<13;i++) basicLeptop[i]->Fill(leptop[i],weight);
         basicLeptons[0]->Fill(leptons[0].pt(),weight);
         basicLeptons[1]->Fill(leptons[0].eta(),weight);
         basicLeptons[2]->Fill(neutrinos[0].pt(),weight); // fix to met
         double hadtop[10]={hadt.m(),hadt.pt(),hadt.eta(),hadt.phi(),
         hadW.m(),hadW.pt(),hadW.eta(),hadW.phi(),
         truth,0};
         for(unsigned i=0;i<10;i++) basicHadtop[i]->Fill(hadtop[i],weight);
         */
        return true;
        /////////////////////////////////////////////////////////////
    } else return false;// recotruth=0;
    //basicLeptop[9]->Fill(recotruth,weight);
    //return true;
}// close recolept
//double complex discriminant = pow(bw,2) -4*aw*cw;
//cout<< creal(discriminant)<<" "<<cimag(discriminant) <<endl;
/////////////////////////////////////////////////////////////////
bool recotlepeq(int & bh, int & bl,vector<PseudoJet> jets, vector<PseudoJet> leptons,vector<PseudoJet> neutrinos, vector<int> btag, vector<int> btrue, double met, double weight, double cut, int type){
    double mw = wmass;// (leptons.at(0)+neutrinos.at(0)).m();//teste
    double wt = (leptons.at(0).px()*neutrinos.at(0).px()) + (leptons.at(0).py()*neutrinos.at(0).py());
    double mu = (pow(mw,2)/2 + wt);
    double aw = (pow(leptons.at(0).pz(),2)-pow(leptons.at(0).e(),2));
    double bw = 2*mu*leptons.at(0).pz();
    double cw = pow(mu,2)-pow(leptons.at(0).e(),2)*pow(neutrinos.at(0).pt(),2); // change to met
    //double mw = (neutrinos.at(0)+leptons.at(0)).m(); // teste
    double discriminant = pow(bw,2) - 4*aw*cw ;
    //
    unsigned int jsize = jets.size(); vector<double> a4,sol;
    for(unsigned int i=0; i< jsize; i++) {
        double mb=4.7;//jets.at(i).m();//teste
        PseudoJet plb = leptons.at(0)+jets.at(i);
        double wtt = (plb.px()*neutrinos.at(0).px()) + (plb.py()*neutrinos.at(0).py());
        double plpj = leptons.at(0).e()*jets.at(i).e() - leptons.at(0).px()*jets.at(i).px() 
        -leptons.at(0).py()*jets.at(i).py()-leptons.at(0).pz()*jets.at(i).pz();//p_l\,.\,p_b
        double mt = 173;// (neutrinos.at(0)+plb).m(); // teste
        double mut = ((pow(mt,2)-pow(mb,2))/2 + wtt -plpj);
        double at = (pow(plb.pz(),2)-pow(plb.e(),2));
        double bt = 2*mut*plb.pz();
        double ct = pow(mut,2)-pow(plb.e(),2)*pow(neutrinos.at(0).pt(),2); // change to met
        /////////////////////////////////////////////////////////////////////////////
        // solution of the quadratic equation
        //double discriminant2 = pow(bt,2) - 4*at*ct ; 
        //double pznu1 = (-bt - sqrt(discriminant2))/(2*at);
        //double pznu2 = (-bt + sqrt(discriminant2))/(2*at);
        //cout<<"quadratic pz1 "<<pznu1 << " pz2 " <<pznu2 <<" pz "<<neutrinos.at(0).pz()<<" btrue "<<btrue[i]<<endl;
        //////////////////////////////////////////////////////////////////////////////
        double X=-9*pow((at*bt + aw*bw),2) + 6*(pow(at,2) + pow(aw,2))*(pow(bt,2) + pow(bw,2) + 2*at*ct + 2*aw*cw);
        double Y = 54*at*pow(aw,2)*pow(bt,3)-108*pow(at,2)*aw*pow(bt,2)*bw+54*pow(aw,3)*pow(bt,2)*bw+
        54*pow(at,3)*bt*pow(bw,2)-108*at*pow(aw,2)*bt*pow(bw,2)+54*pow(at,2)*aw*pow(bw,3)-
        108*pow(at,2)*pow(aw,2)*bt*ct-108*pow(aw,4)*bt*ct+108*pow(at,3)*aw*bw*ct+
        108*at*pow(aw,3)*bw*ct+108*pow(at,3)*aw*bt*cw+108*at*pow(aw,3)*bt*cw-108*pow(at,4)*bw*cw-
        108*pow(at,2)*pow(aw,2)*bw*cw;
        //double complex det = Y + csqrt(cpow(Y,2) + 4*cpow(X,3));
        double det = 0;
        double AA = X/(3*pow(2.,0.666667)*(pow(at,2) + pow(aw,2)));
        /*double complex pnuzsol =-(at*bt + aw*bw)/(2*(pow(at,2) + pow(aw,2))) 
         -cpow(det,-0.333333)*AA
         +cpow(det,0.333333)/(6*pow(2.,0.333333)*(pow(at,2) + pow(aw,2)));
         if(abs(cimag(pnuzsol))<0.1){ 
         a4.push_back( pow(ct+creal(pnuzsol)*(bt+at*creal(pnuzsol)),2) +
         pow(cw+creal(pnuzsol)*(bw+aw*creal(pnuzsol)),2) ); 
         sol.push_back(creal(pnuzsol));
         //cout<<i<<" pz3 "<<creal(pnuzsol)<<" "<< a4[i]<<" "<< (neutrinos.at(0)+plb).m()<<" "<<btrue[i]<<endl;
         }// close if no imaginary    */
    } //close loop on jets
    // make the neutrino vector 
    int minM = TMath::LocMin(a4.size(), &a4[0]); //cout<<minM<<endl;
    double enu = sqrt(pow(neutrinos.at(0).px(),2)+pow(neutrinos.at(0).py(),2)+pow(sol[minM],2));
    double pxnu=neutrinos.at(0).px(),pynu=neutrinos.at(0).py();
    // the solution jet is the minimum 
    PseudoJet lepW = leptons.at(0)+ fastjet::PseudoJet(pxnu,pynu,sol[minM],enu);
    PseudoJet lept = jets[minM] + lepW; bl = minM;  
    double pnuzerror=(neutrinos.at(0)+leptons.at(0)).m()-lepW.m();
    int bll; for(unsigned int i=0; i< jsize; i++) if(btrue[i]==5) bll=i; // to know true mass
    double mterror=(neutrinos.at(0)+leptons.at(0)+jets[bll]).m()-lept.m();
    //////////////////////////////////////////////////////////////
    PseudoJet lepWtransv = leptons.at(0) + jets.at(bl);
    double wtransM = // w transverse mass
    sqrt(
         2*neutrinos.at(0).e()*leptons.at(0).e()-
         2*(neutrinos.at(0).px()*leptons.at(0).px()+neutrinos.at(0).py()*leptons.at(0).py())); // fix neutrino
    double ttransM = // t transverse mass
    sqrt(
         2*neutrinos.at(0).e()*lepWtransv.e()-
         2*(neutrinos.at(0).px()*lepWtransv.px()+neutrinos.at(0).py()*lepWtransv.py())); // fix neutrino
    /////////////////////////////////////////////////////////////////// 
    // find hadronic w -- closest hadronic mass -- but not the hadronic top
    vector<int> j1,j2; vector<double> a3; int minMw;
    for(unsigned int nj1=0; nj1< jsize; nj1++) 
        for(unsigned int nj2=nj1+1; nj2< jsize; nj2++) 
            if (nj1!=nj2 && nj1!=bl && nj2!=bl )
                //if(btag[nj1]==0 && btag[nj2]==0)
            {  //std::cout<<nj1<<" "<<nj2<<" "<<bl<<" jsize "<<jsize<<std::endl; 
                double invmassA =  (jets.at(nj1)+jets.at(nj2)).m();
                a3.push_back((invmassA-wmass)*(invmassA-wmass)); 
                j1.push_back(nj1); j2.push_back(nj2); 
                // we also what to keep the nj...           
            } // loop on jets  
    minMw = TMath::LocMin(a3.size(), &a3[0]);
    PseudoJet hadW = jets[j1[minMw]] + jets[j2[minMw]];
    ////////////////////////////////////////////////////////////
    // do the hadronic top with the hardest pt jet
    vector<double> ptj; vector<int> nptj;
    for(unsigned int nj1=0; nj1< jsize; nj1++) 
        if (nj1!=j1[minMw] && nj1!=j2[minMw] && nj1!=bl)
        { ptj.push_back(jets.at(nj1).pt()); nptj.push_back(nj1);  }
    int maxPt = TMath::LocMax(ptj.size(), &ptj[0]); 
    PseudoJet hadt = hadW + jets.at(nptj[maxPt]);
    ////////////////////////////////////////////////////////////////////
    // test truth 
    int truth=0; 
    if(btrue[j1[minMw]]==0 && btrue[j2[minMw]]==0  && btrue[bl]==-5 && btrue[nptj[maxPt]]==5 ) truth=1; else truth=0; 
    //cout<<j1[minMw]<<" "<<j2[minMw]<<" "<<bl<<" "<<nptj[maxPt]<<" "<<endl; 
    ///////////////////////////////////////////////////////////////////
    double detabb  = abs(jets[nptj[maxPt]].eta() - jets[bl].eta());
    /////////////////////////////////////////////////////////////////// 
    // fill diagrams
    double leptop[13]={lept.m(),lept.pt(),lept.eta(),lept.phi(),
        lepW.m(),lepW.pt(),lepW.eta(),lepW.phi(),
        pnuzerror,0,0,wtransM,ttransM}; // pnuzerror,truth,mterror,wmt,tmt
    for(unsigned i=0;i<13;i++) basicLeptop[i]->Fill(leptop[i],weight);
    basicLeptons[0]->Fill(leptons[0].pt(),weight);
    basicLeptons[1]->Fill(leptons[0].eta(),weight);
    basicLeptons[2]->Fill(met,weight);
    double hadtop[10]={hadt.m(),hadt.pt(),hadt.eta(),hadt.phi(),
        hadW.m(),hadW.pt(),hadW.eta(),hadW.phi(),
        truth,detabb};
    for(unsigned i=0;i<10;i++) basicHadtop[i]->Fill(hadtop[i],weight);
    return true;
} // close recotlepeq
/////////////////////////////////////////////////////////////////
bool recohadt(int & bh, int & bl, vector<PseudoJet> jets, vector<PseudoJet> leptons,vector<PseudoJet> neutrinos, vector<int> btag, vector<int> btrue, double met, double weight, double cut, int type){
    ///////////////////////////////////////////////////////////////////
    // gen level info
    int blll=-1,bhhh=-1; vector<int> wj;
    unsigned int jsize = jets.size();
    for(unsigned int nj1=0; nj1< jsize; nj1++) // only works at parton level
        if(btrue[nj1]==5) bhhh=nj1; else if(btrue[nj1]==-5) blll=nj1; else wj.push_back(nj1);   
    //cout<<wj.size()<<" "<<jsize<<" "<<blll<<" "<<bhhh<<" "<<endl;  
    if(bhhh!=-1 && blll!=-1 && wj.size()>1){
        // cout<<jsize<<endl;
        //cout<<wj[0]<<" "<<wj[1]<<" "<<blll<<" "<<bhhh<<" "<<endl;  
        PseudoJet lepTtrue = leptons.at(0) + neutrinos.at(0) + jets.at(blll);
        PseudoJet hadTtrue = jets.at(wj[0]) + jets.at(wj[1]) + jets.at(bhhh);
        leptop->Fill(lepTtrue.m(),weight);
        hadtop->Fill(hadTtrue.m(),weight);
        genbhad->Fill(jets.at(bhhh).pt(),weight);
        genblep->Fill(jets.at(blll).pt(),weight);
        genbhadeta->Fill(jets.at(bhhh).eta(),weight);
        genblepeta->Fill(jets.at(blll).eta(),weight);
    }
    /////////////////////////////////////////////////////////////
    // find the three jets that are hadronic t
    //unsigned int jsize = jets.size();
    vector<int> j1,j2,j3; vector<double> a3; int minM;
    for(unsigned int nj1=0; nj1< jsize; nj1++) 
        for(unsigned int nj2=0; nj2< jsize; nj2++) if (nj2!=nj1)
            for(unsigned int nj3=0; nj3< jsize; nj3++) if (nj3!=nj2 && nj3!=nj1){
                // closest from the quadratic sum of w and t 
                PseudoJet wteste = jets.at(nj1)+ jets.at(nj2);
                PseudoJet tteste = wteste+ jets.at(nj3);
                double invmass = (wteste.m()-wmass)*(wteste.m()-wmass) + (tteste.m()-tmass)*(tteste.m()-tmass);
                a3.push_back(invmass); 
                j1.push_back(nj1); j2.push_back(nj2);j3.push_back(nj3); 
            } // loop on jets 
    minM = TMath::LocMin(a3.size(), &a3[0]);
    PseudoJet hadW = jets[j1[minM]] + jets[j2[minM]];
    PseudoJet hadt = jets[j3[minM]] + hadW;
    ///////////////////////////////////////////////////////////////////
    // take the other hardest jet for the lep t transverse mass
    vector<double> ptj; vector<int> nptj;
    for(unsigned int nj1=0; nj1< jsize; nj1++) 
        if (nj1!=j1[minM] && nj1!=j2[minM] && nj1!=j3[minM])
        { ptj.push_back(jets.at(nj1).pt()); nptj.push_back(nj1);  }
    int maxPt = TMath::LocMax(ptj.size(), &ptj[0]);
    PseudoJet lepWtransv = leptons.at(0) + jets.at(nptj[maxPt]);
    //////////////////////////////////////////////////////////////////
    // chose the neutrino solving to mw2
    double wt = (leptons.at(0).px()*neutrinos.at(0).px()) + (leptons.at(0).py()*neutrinos.at(0).py());
    double mu = (pow(wmass,2)/2 + wt);
    double aw = (pow(leptons.at(0).pz(),2)-pow(leptons.at(0).e(),2));
    double bw = 2*mu*leptons.at(0).pz();
    double cw = pow(mu,2)-pow(leptons.at(0).e(),2)*pow(met,2);
    double discriminant = pow(bw,2) - 4*aw*cw ; 
    double pnuzerror,pnuz=-10;//,recowm,recowpt,recoweta,recowphi;  
    PseudoJet lepW; int recotruth=0; double mterror, recoWmass=-10, recotopmass=-10;
    if(discriminant>=0){
        //recotruth=1; 
        double pznu1 = (-bw - sqrt(discriminant))/(2*aw);
        double pznu2 = (-bw + sqrt(discriminant))/(2*aw);
        //cout<<"pz1 "<<pznu1 << " pz2 " <<pznu2 <<" pz "<<neutrinos.at(0).pz()<<" pzl "<<leptons.at(0).pz()<<endl;
        //cout<<"dumb "<< (neutrinos.at(0)+leptons.at(0)).m() << " calculated " <<
        //sqrt(pow(leptons.at(0).e()+neutrinos.at(0).e(),2)-wt-pow(leptons.at(0).pz()+neutrinos.at(0).pz(),2))
        //  <<endl;
        pnuz=TMath::Min(pznu1,pznu2); 
        double enu = sqrt(pow(neutrinos.at(0).px(),2)+pow(neutrinos.at(0).py(),2)+pow(pnuz,2));
        double pxnu=neutrinos.at(0).px(),pynu=neutrinos.at(0).py();
        lepW = leptons.at(0)+ fastjet::PseudoJet(pxnu,pynu,pnuz,enu);
        recoWmass=lepW.m(); //cout<<recoWmass<<endl;
        recotopmass = (lepW+ jets.at(nptj[maxPt])).m();
        mterror=(neutrinos.at(0)+leptons.at(0)+jets.at(nptj[maxPt])).m()-recotopmass;
        pnuzerror=(neutrinos.at(0)+leptons.at(0)).m()-lepW.m();
        if (abs(pnuz-neutrinos.at(0).pz())<0.1) recotruth=1;
    } 
    ////////////////////////////////////////////////////////////////////
    // define tranverse masses
    double wwmass = pow(neutrinos.at(0).e()+leptons.at(0).e(),2)
    - pow(neutrinos.at(0).px()+leptons.at(0).px(),2)
    - pow(neutrinos.at(0).py()+leptons.at(0).py(),2)
    - pow(neutrinos.at(0).pz()+leptons.at(0).pz(),2);
    double wtransE = sqrt(pow(wmass,2) + pow(neutrinos.at(0).px()+leptons.at(0).px(),2)
                          + pow(neutrinos.at(0).py()+leptons.at(0).py(),2));// w transverse mass                   
    double wtransM = pow(
                         pow(neutrinos.at(0).pt()+leptons.at(0).pt(),2)-
                         pow(neutrinos.at(0).px()+leptons.at(0).px(),2)-pow(neutrinos.at(0).py()+leptons.at(0).py(),2)
                         ,0.5); // fix neutrino
    //cout<<neutrinos.at(0).pt()<<" "<<leptons.at(0).pt()<<" "<<wtransM<<endl;
    double btransE = sqrt( pow(bmass,2) +jets.at(nptj[maxPt]).px()*jets.at(nptj[maxPt]).px()+
                          jets.at(nptj[maxPt]).py()*jets.at(nptj[maxPt]).py());
    double ttransM = sqrt(pow(wtransE + btransE ,2)-
                          pow(jets.at(nptj[maxPt]).px()+leptons.at(0).px()+neutrinos.at(0).px(),2)- // fix neutrino
                          pow(jets.at(nptj[maxPt]).py()+leptons.at(0).py()+neutrinos.at(0).py(),2));
    ////////////////////////////////////////////////////////////////////
    // define alternative transverse quantities
    double lbmass = sqrt( pow(jets.at(nptj[maxPt]).e()+leptons.at(0).e(),2)
                         - pow(jets.at(nptj[maxPt]).px()+leptons.at(0).px(),2)
                         - pow(jets.at(nptj[maxPt]).py()+leptons.at(0).py(),2)
                         - pow(jets.at(nptj[maxPt]).pz()+leptons.at(0).pz(),2));
    double lbtransE = sqrt(pow(lbmass,2) + pow(jets.at(nptj[maxPt]).px()+leptons.at(0).px(),2)
                           + pow(jets.at(nptj[maxPt]).py()+leptons.at(0).py(),2));// w transverse mass    
    double ttransMal = sqrt(pow(lbtransE + neutrinos.at(0).pt(),2)-
                            pow(jets.at(nptj[maxPt]).px()+leptons.at(0).px()+neutrinos.at(0).px(),2)- // fix neutrino
                            pow(jets.at(nptj[maxPt]).py()+leptons.at(0).py()+neutrinos.at(0).py(),2));
    ////////////////////////////////////////////////////////////////////
    // define razor variables
    //cout<<wtransM<<" "<<btransE<<" "<<ttransM<<endl;
    //cout<<ttransM<<endl;
    double mrazor = sqrt(
                         pow( 
                             sqrt(pow(leptons.at(0).px(),2)+pow(leptons.at(0).py(),2)+pow(leptons.at(0).pz(),2)) +
                             sqrt(pow(jets.at(nptj[maxPt]).px(),2)+pow(jets.at(nptj[maxPt]).py(),2)+pow(jets.at(nptj[maxPt]).pz(),2))
                             ,2)-
                         pow(jets.at(nptj[maxPt]).pz()+leptons.at(0).pz(),2)
                         );
    double mrazortrans = sqrt(
                              neutrinos.at(0).pt()*(jets.at(nptj[maxPt]).pt()+leptons.at(0).pt())-
                              neutrinos.at(0).px()*(jets.at(nptj[maxPt]).px()+leptons.at(0).px())-
                              neutrinos.at(0).py()*(jets.at(nptj[maxPt]).py()+leptons.at(0).py()) 
                              )*sqrt(0.5);
    double ratiorazor = mrazortrans/mrazor; 
    /////////////////////////////////////////////////////////////////// 
    // test trueth
    int truth;
    if(btrue[j1[minM]]==0 && btrue[j2[minM]]==0 && btrue[j3[minM]]==5 && btrue[nptj[maxPt]]==-5 ) truth=1; else truth=0; 
    //std::cout<<j1[minM]<<" "<<j2[minM]<<" "<<j3[minM]<<" "<<nptj[maxPt]<<" "<<std::endl; 
    //std::cout<<btrue[j1[minM]]<<" "<<btrue[j2[minM]]<<" "<<btrue[j3[minM]]<<" "<<btrue[nptj[maxPt]]<<" "<<std::endl; 
    ///////////////////////////////////////////////////////////////////
    double detabb  = abs(jets[nptj[maxPt]].eta() - jets[j3[minM]].eta());
    double detalb  = abs(leptons.at(0).eta() - jets[j3[minM]].eta());
    /////////////////////////////////////////////////////////
    // fill diagrams
    /////////////////////////////////////////////////////////
    if(1>0
       && ttransM > wbtransmassMin
       && ttransM < wbtransmassMax
       //     && mrazor>150
       //&& ratiorazor<0.8 
       //&& wtransM>80
       ){
        double leptop[19]={recotopmass,0,0,0, //lept.m(),lept.pt(),lept.eta(),lept.phi(),
            recoWmass,0,0,0, //lepW.m(),lepW.pt(),lepW.eta(),lepW.phi(),
            pnuzerror,recotruth,0,wtransM,ttransM,detalb,mrazor,mrazortrans, ratiorazor,ttransMal,lbmass
            //pnuzerror,truth,mterror,wmt,tmt};
        };
        for(unsigned i=0;i<19;i++) basicLeptop[i]->Fill(leptop[i],weight);
        //
        basicLeptons[0]->Fill(leptons[0].pt(),weight);
        basicLeptons[1]->Fill(leptons[0].eta(),weight);
        basicLeptons[2]->Fill(neutrinos[0].pt(),weight); // fix to met
        double hadtop[10]={hadt.m(),hadt.pt(),hadt.eta(),hadt.phi(),
            hadW.m(),hadW.pt(),hadW.eta(),hadW.phi(),
            truth,detabb};
        for(unsigned i=0;i<10;i++) basicHadtop[i]->Fill(hadtop[i],weight);
    }
    // reco the lep top
    return true;//} else return false;
}// close top reco
/////////////////////////////////////////////////////////////////
int recol( vector<PseudoJet> jets,vector<PseudoJet> leptons,vector<PseudoJet> neutrinos){
    // from the jet collection find the hadronic W
    // construct met from all the rest
    // do a isolation vector
    vector<double> LepIso; int nlep=0; int nlepsurvive=0;
    //unsigned int MinDRLep;
    unsigned int jsize = jets.size();
    for(unsigned int j = 0;j<leptons.size();j++) {
        for(unsigned int i = 0;i<jsize;i++) LepIso.push_back(leptons.at(j).delta_R(jets.at(i)));
        if(leptons.size()>1) for(unsigned int i = 0;i<leptons.size();i++) if (i!=j) LepIso.push_back(leptons.at(j).delta_R(leptons.at(i)));
        double MinDRLep = TMath::LocMin(LepIso.size(), &LepIso[0]);
        if(1>0
           && LepIso[MinDRLep] > lepiso  
           && leptons.at(j).pt()> ptlepton 
           && abs(leptons.at(j).eta())< etal
           ) nlep++; // close basic cuts
    } // close for nlep
    //cout<<nlep<<endl;
    return nlep;
}// close recow
/////////////////////////////////////////////////////////////////
void isbtagged(vector<PseudoJet> jets, vector<int> & btag, vector<int> & bmistag, vector<int> & btrue);
/////////////////////////////////////////////////////////////////
int recojets(vector<PseudoJet> particles,vector<PseudoJet> & jets_akt, vector<int> & btag, vector<int> & bmistag, vector<int> & fattag, vector<int> & btrue){
    JetDefinition akt(antikt_algorithm, RR);
    ClusterSequence cs_akt(particles, akt);
    //vector<PseudoJet> jets_akt;
    
    Selector jet_selector = SelectorPtMin(jet_ptmin) && SelectorAbsRapMax(rapmax);
    if(shower){
        // first we do akt jets from particles
        jets_akt = sorted_by_pt(jet_selector(cs_akt.inclusive_jets()));
    }
    else{
        double const ptmin=0.0; // if parton jet pt min zero
        Selector jet_selector_parton = SelectorPtMin(ptmin);
        jets_akt = sorted_by_pt(jet_selector_parton(cs_akt.inclusive_jets()));
    }
    //////////////////// do a jet cut
    vector<PseudoJet> jets_final;     unsigned int njets;
    for (unsigned int i = 0; i < jets_akt.size(); i++) if(jets_akt.at(i).pt()>jet_ptminfinal) jets_final.push_back(jets_akt.at(i));
    //for (unsigned int i = 0; i < particles.size(); i++) if(particles.at(i).pt()<jet_ptminfinal) {njets =0; return njets; } //jets_final.push_back(particles.at(i));
    njets = jets_final.size();
    //cout<<njets<<endl;
    isbtagged(jets_akt, btag, bmistag,btrue); // check wheather the b(c)jet is b--(mis)tagable
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
    //jets = jets_akt;
    return njets;
} // close cluster jets
////////////////////////////////////////////////////////////////////////////////////////////////
void isbtagged(vector<PseudoJet> jets, vector<int> & btag, vector<int> & bmistag, vector<int> & btrue){ 
    unsigned int jsize = jets.size();
    for (unsigned int i=0; i<jsize; i++) { // check wheter jet have inside a b's are taggable 
        vector<PseudoJet> constitu=jets.at(i).constituents();
        unsigned int csize = constitu.size();
        int see=0,see2=0,seetruth=0,id;
        for (unsigned int j=0; j<csize; j++) {
            //cout<<"constituents flavour "<<constitu.at(j).user_index()<<endl;
            if(constitu.at(j).user_index() == 5 || constitu.at(j).user_index() == -5) {
                seetruth++; id=constitu.at(j).user_index();
                if(constitu.at(j).pt() > bjetpt && constitu.at(j).eta() < etab) see++;// to reco btag 
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
} // close isbtagged
/////////////////////////////////////////////////////////////////////////
// save the histos
int save_hist(int isample,int reco,int sample){
    const char* Mass;
    Mass = Form("Control_reco_%d_place_%d_.root",reco,isample); cout<<sample<<endl;
    TFile f1(Mass, "recreate");
    f1.cd();
    Njets_passing_kLooseID->Write();
    btagselected->Write();
    genmbb->Write();
    leptop->Write();
    hadtop->Write();
    genbhad->Write();
    genblep->Write();
    genbhadeta->Write();
    genblepeta->Write();
    //basicLeptons[0]->Write();
    //basicLeptons[1]->Write();
    //basicLeptons[2]->Write();
    //for(unsigned i=0;i<10;i++) basicHadtop[i]->Write();
    //for(unsigned i=0;i<19;i++) basicLeptop[i]->Write();
    f1.Close();
    //
    Njets_passing_kLooseID->Reset();
    btagselected->Reset();
    genmbb->Reset();
    leptop->Reset();
    hadtop->Reset();
    genbhad->Reset();
    genblep->Reset();
    genbhadeta->Reset();
    genblepeta->Reset();
    //basicLeptons.clear();
    //basicHadtop.clear();
    //basicLeptop.clear();
    //  basicLeptons[0]->Reset();
    //  basicLeptons[1]->Reset();
    //  basicLeptons[2]->Reset();
    //  for(unsigned i=0;i<10;i++) basicHadtop[i]->Reset();
    //  for(unsigned i=0;i<13;i++) basicLeptop[i]->Reset();
    cout<<sample<<endl;
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
    
    TH1D *MetMass = new TH1D("MetMass_ct4",  
                             label, 
                             50, 0, 300);
    MetMass->GetXaxis()->SetTitle("MET (GeV)");
    basicLeptons.push_back (MetMass);
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
