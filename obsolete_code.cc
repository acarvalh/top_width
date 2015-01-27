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

