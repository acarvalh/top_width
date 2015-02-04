organization of files

    string path[1]={
        "/afs/cern.ch/work/a/acarvalh/tt_fulylep/ttbar_gen180/"
        //"/data2/TopWidth/MG5_aMC_v2_1_0/wwbbfullep/top_Wvary/",
        //"/data2/TopWidth/MG5_aMC_v2_1_0/wwbbfullep_OnOn/top_Wvary/",
        //"/data2/TopWidth/MG5_aMC_v2_1_0/wwbbfullep_OnOff/top_Wvary/",
        //"/data2/TopWidth/MG5_aMC_v2_1_0/wwbbfullep_OffOn/top_Wvary/",
        //"/data2/TopWidth/MG5_aMC_v2_1_0/wwbbfullep_OffOff/top_Wvary/"
    };

    //string sample[7] = {"Wt_0","Wt_1","Wt_2","Wt_3","Wt_4","Wt_5","Wt_6"};// the last one have 50k

see one particular event

                    //if(ievent== 32634) {
                    //    cout<<"jet1 "<<particles.at(0).px()<<" "<<particles.at(0).py()<<" "<<particles.at(0).pz()<<" "<<particles.at(0).e()<<endl;
                    //    cout<<"jet2 "<<particles.at(1).px()<<" "<<particles.at(1).py()<<" "<<particles.at(1).pz()<<" "<<particles.at(1).e()<<endl;
                    //    cout<<"jet3 "<<particles.at(2).px()<<" "<<particles.at(2).py()<<" "<<particles.at(2).pz()<<" "<<particles.at(2).e()<<endl;
                    //    cout<<"jet4 "<<particles.at(3).px()<<" "<<particles.at(3).py()<<" "<<particles.at(3).pz()<<" "<<particles.at(3).e()<<endl;
                    //    cout<<"lepton "<<leptons.at(0).px()<<" "<<leptons.at(0).py()<<" "<<leptons.at(0).pz()<<" "<<leptons.at(0).e()<<endl;
                    //    cout<<"neutrino "<<particles.at(0).px()<<" "<<particles.at(0).py()<<" "<<particles.at(0).pz()<<" "<<particles.at(0).e()<<endl;
                    //}


checks


cout<<"njets "<<njets<<" nleptons "<<counterl<<" pzl "<< neutrinos.at(0).pz()<<endl;

                  //
                    //int count=0; for(unsigned int i = 0; i < njets; i++) if(btag[i]>0)count++;

functions to semileptonic

                        //if(lepcuts && true_tops==type )lepwreco = recotlepeq(bh,bl,jets,leptons,neutrinos,btag,btrue,met, weight);
                        //int true_tops = truetops(jets,leptons,neutrinos,btag,btrue);  

print results


               cout<<"\n\n total nevents 50/fb "<< neventsall<<" "
               //<< (finalevents0[0]*CX[0])*1000*lumi/nevents << " "
               //<< (finalevents0[1]*CX[1]+finalevents0[2]*CX[2]+finalevents0[3]*CX[3])*1000*lumi/nevents << " "
               //<< (finalevents0[0]/nevents)*(CX[0]+CX[1]+CX[2])*1000*lumi
               <<endl; 


                cout<< label[type] <<" nevents 50/fb = "<< neventslumi << " net efficiency " << neventslumi/neventsall <<endl;
               //





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



