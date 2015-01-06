{
    /////////////////////////////////////////////
    // here we put the options for ploting
    TStyle *defaultStyle = new TStyle("defaultStyle","Default Style");
    defaultStyle->SetOptStat(0000);
    defaultStyle->SetOptFit(000); 
    defaultStyle->SetPalette(1);
    /////// pad ////////////
    defaultStyle->SetPadBorderMode(1);
    defaultStyle->SetPadBorderSize(1);
    defaultStyle->SetPadColor(0);
    defaultStyle->SetPadTopMargin(0.05);
    defaultStyle->SetPadBottomMargin(0.13);
    defaultStyle->SetPadLeftMargin(0.13);
    defaultStyle->SetPadRightMargin(0.02);
    /////// canvas /////////
    defaultStyle->SetCanvasBorderMode(0);
    defaultStyle->SetCanvasColor(0);
    defaultStyle->SetCanvasDefH(600);
    defaultStyle->SetCanvasDefW(600);
    /////// frame //////////
    defaultStyle->SetFrameBorderMode(0);
    defaultStyle->SetFrameBorderSize(1);
    defaultStyle->SetFrameFillColor(0); 
    defaultStyle->SetFrameLineColor(1);
    /////// label //////////
    defaultStyle->SetLabelOffset(0.005,"XY");
    defaultStyle->SetLabelSize(0.07,"XY");
    defaultStyle->SetLabelFont(46,"XY");
    /////// title //////////
    defaultStyle->SetTitleOffset(1.1,"X");
    defaultStyle->SetTitleSize(0.01,"X");
    defaultStyle->SetTitleOffset(1.25,"Y");
    defaultStyle->SetTitleSize(0.07,"Y");
    defaultStyle->SetTitleFont(44, "XYZ");
    /////// various ////////
    defaultStyle->SetNdivisions(505,"Y");
    defaultStyle->SetLegendBorderSize(0);  // For the axis titles:
    
    defaultStyle->SetTitleColor(1, "XYZ");
    defaultStyle->SetTitleFont(42, "XYZ");
    defaultStyle->SetTitleSize(0.06, "XYZ");
    
    // defaultStyle->SetTitleYSize(Float_t size = 0.02);
    defaultStyle->SetTitleXOffset(0.9);
    defaultStyle->SetTitleYOffset(1.05);
    // defaultStyle->SetTitleOffset(1.1, "Y"); // Another way to set the Offset
    
    // For the axis labels:
    defaultStyle->SetLabelColor(1, "XYZ");
    defaultStyle->SetLabelFont(42, "XYZ");
    defaultStyle->SetLabelOffset(0.007, "XYZ");
    defaultStyle->SetLabelSize(0.04, "XYZ");
    
    // For the axis:
    defaultStyle->SetAxisColor(1, "XYZ");
    defaultStyle->SetStripDecimals(kTRUE);
    defaultStyle->SetTickLength(0.03, "XYZ");
    defaultStyle->SetNdivisions(510, "XYZ");
    defaultStyle->SetPadTickX(1);   // To get tick marks on the opposite side of the frame
    defaultStyle->SetPadTickY(1);
    defaultStyle->cd();
    /////////////////////////////////////////////////////
    int nmass = 30;
    const char* channel[nmass]={
        "Control_reco_0_place_0_.root",
        "Control_reco_1_place_0_.root",
        "Control_reco_2_place_0_.root",
        "Control_reco_3_place_0_.root",
        "Control_reco_4_place_0_.root",
        "Control_reco_5_place_0_.root",
        "Control_reco_6_place_0_.root",
        "Control_reco_7_place_0_.root",
        "Control_reco_8_place_0_.root",
        "Control_reco_9_place_0_.root",
        //
        "Control_reco_0_place_2_.root",
        "Control_reco_1_place_2_.root",
        "Control_reco_2_place_2_.root",
        "Control_reco_3_place_2_.root",
        "Control_reco_4_place_2_.root",
        "Control_reco_5_place_2_.root",
        "Control_reco_6_place_2_.root",
        "Control_reco_7_place_2_.root",
        "Control_reco_8_place_2_.root",
        "Control_reco_9_place_2_.root",
        //
        "Control_reco_0_place_4_.root",
        "Control_reco_1_place_4_.root",
        "Control_reco_2_place_4_.root",
        "Control_reco_3_place_4_.root",
        "Control_reco_4_place_4_.root",
        "Control_reco_5_place_4_.root",
        "Control_reco_6_place_4_.root",
        "Control_reco_7_place_4_.root",
        "Control_reco_8_place_4_.root",
        "Control_reco_9_place_4_.root"
    };
    const char* lege[nmass]={"UU 0.85","DD 0.85","DU 0.85","UD 0.85","UOn 0.85","DOn 0.85","OnD 0.85","OnU 0.85","OnOn 0.85","full 0.85",
        "UU","DD","DU","UD","UOn","DOn","OnD","OnU","OnOn","full",
        "UU 1.15","DD 1.15","DU 1.15","UD 1.15","UOn 1.15","DOn 1.15","OnD 1.15","OnU 1.15","OnOn 1.15","full 1.15"};
    int maxtodo=8;
    int todo[maxtodo]={18,17,16,14,15,11,13,12};//14,17,18};//13 ,2 ,1 };//27,28,34, // 0,21,16,13, 9 ,5 ,
    double masses[maxtodo] = {8,17,16};//,23,1,1,1,1,1}; //8,0,1,2,3,4,5,6,7,9};//,6,7,8,9};//,
    //
    TLegend *leg = new TLegend(0.65,0.60,0.99,0.99);
    leg->SetTextSize(0.04146853);
    leg->SetLineColor(1);
    leg->SetLineStyle(1);
    leg->SetLineWidth(1);
    leg->SetFillColor(0);
    //
    int nplots =39;
    const char* namplots[nplots] = {
        "njets_passing_kLooseID_ct4.png",
        "btagselected.png",
        "E1histpt.png",  
        "E1histeta.png",  
        "MetMass_ct4.png",
        "H1hist.png",
        "H1histpt.png",
        "H1histeta.png",
        "H1histphi.png",
        "HW1hist.png",
        "HW1histpt.png",
        "HW1histeta.png",
        "HW1histphi.png",
        "truth.png",
        "detabb.png",
        //
        "H1LepThist.png",
        "H1LepThistpt.png",
        "H1LepThisteta.png",
        "H1LepThistphi.png",
        "HW1LepThist.png",
        "HW1LepThistpt.png",
        "HW1LepThisteta.png",
        "HW1LepThistphi.png",
        "pnuzerror.png",
        "recotruthlept.png",
        //
        "hadtop.png",
        "leptop.png",
        //
        "wmt.png",
        "tmt.png",
        "detalb.png",
        //
        "mraz.png",
        "mrazt.png",
        "razratio.png",
        "tmtal.png",
        "mlb.png",
        "genbhad.png",
        "genblep.png",
        "genbhadeta.png",
        "genblepeta.png"
    };
    
    TH1D* plots[maxtodo][nplots];//[file][plot]
    TFile *file[maxtodo];
    for(int i=0;i<maxtodo;i++){
        TFile *file[i] = TFile::Open(channel[todo[i]]);
        cout<<channel[todo[i]]<<endl;
        TH1D* plots[i][0] = (TH1D* ) file[i]->Get("njets_passing_kLooseID_ct4;1"); 
        TH1D* plots[i][1] = (TH1D* ) file[i]->Get("btagselected;1"); 
        TH1D* plots[i][2] = (TH1D* ) file[i]->Get("E1histpt;1");   
        TH1D* plots[i][3] = (TH1D* ) file[i]->Get("E1histeta;1");   
        TH1D* plots[i][4] = (TH1D* ) file[i]->Get("MetMass_ct4;1"); 
        TH1D* plots[i][5] = (TH1D* ) file[i]->Get("H1hist;1"); 
        TH1D* plots[i][6] = (TH1D* ) file[i]->Get("H1histpt;1"); 
        TH1D* plots[i][7] = (TH1D* ) file[i]->Get("H1histeta;1"); 
        TH1D* plots[i][8] = (TH1D* ) file[i]->Get("H1histphi;1");
        TH1D* plots[i][9] = (TH1D* ) file[i]->Get("HW1hist;1"); 
        TH1D* plots[i][10] = (TH1D* ) file[i]->Get("HW1histpt;1"); 
        TH1D* plots[i][11] = (TH1D* ) file[i]->Get("HW1histeta;1"); 
        TH1D* plots[i][12] = (TH1D* ) file[i]->Get("HW1histphi;1");
        TH1D* plots[i][13] = (TH1D* ) file[i]->Get("recotruth;1");
        TH1D* plots[i][14] = (TH1D* ) file[i]->Get("detabb;1");
        //
        TH1D* plots[i][15] = (TH1D* ) file[i]->Get("H1LepThist;1");
        TH1D* plots[i][16] = (TH1D* ) file[i]->Get("H1LepThistpt;1");
        TH1D* plots[i][17] = (TH1D* ) file[i]->Get("H1LepThisteta;1");
        TH1D* plots[i][18] = (TH1D* ) file[i]->Get("H1LepThistphi;1");
        TH1D* plots[i][19] = (TH1D* ) file[i]->Get("HW1LepThist;1");
        TH1D* plots[i][20] = (TH1D* ) file[i]->Get("HW1LepThistpt;1");
        TH1D* plots[i][21] = (TH1D* ) file[i]->Get("HW1LepThisteta;1");
        TH1D* plots[i][22] = (TH1D* ) file[i]->Get("HW1LepThistphi;1");
        TH1D* plots[i][23] = (TH1D* ) file[i]->Get("pnuzerror;1");
        TH1D* plots[i][24] = (TH1D* ) file[i]->Get("recotruthlept;1");
        //
        TH1D* plots[i][25] = (TH1D* ) file[i]->Get("hadtop1;1");
        TH1D* plots[i][26] = (TH1D* ) file[i]->Get("leptop1;1");
        //
        TH1D* plots[i][27] = (TH1D* ) file[i]->Get("wmt;1");
        TH1D* plots[i][28] = (TH1D* ) file[i]->Get("tmt;1");
        TH1D* plots[i][29] = (TH1D* ) file[i]->Get("detalb;1");
        //
        TH1D* plots[i][30] = (TH1D* ) file[i]->Get("mraz;1");
        TH1D* plots[i][31] = (TH1D* ) file[i]->Get("mrazt;1");
        TH1D* plots[i][32] = (TH1D* ) file[i]->Get("razratio;1");
        TH1D* plots[i][33] = (TH1D* ) file[i]->Get("tmtal;1");
        TH1D* plots[i][34] = (TH1D* ) file[i]->Get("massbl;1");
        TH1D* plots[i][35] = (TH1D* ) file[i]->Get("genbhad;1");
        TH1D* plots[i][36] = (TH1D* ) file[i]->Get("genblep;1");
        TH1D* plots[i][37] = (TH1D* ) file[i]->Get("genbhadeta;1");
        TH1D* plots[i][38] = (TH1D* ) file[i]->Get("genblepeta;1");
    }
    const int sigcolor[nmass]={1,28,90,8,93,6,7,8,9,50,
        1,28,90,8,93,6,7,8,9,50,
        1,28,90,8,93,6,7,8,9,50}; // 1,2,3,4,5,6 };
    //const int sigline[nmass]={5,5,5,5,1,1,1,1,1,10}; // 1,2,3,4,5,6 };
    /*const int sigline[nmass]={5,5,5,5,5,5,5,5,5,5,
     1,1,1,1,1,1,1,1,1,1,
     10,10,10,10,10,10,10,10,10,10};
     */
    const int sigline[nmass]={5,5,5,5,1,1,1,1,1,10,
        5,5,5,5,1,1,1,1,1,10,
        5,5,5,5,1,1,1,1,1,10};
    for(int k=0;k<nplots;k++) for(int l=0;l<maxtodo;l++){
        plots[l][k]->SetLineColor(sigcolor[todo[l]]);
        plots[l][k]->SetLineStyle(sigline[todo[l]]);
        plots[l][k]->SetLineWidth(3);
        //cout<<"here "<<k<<" "<<l<<endl;
    }
    TCanvas* PT_HAT = new TCanvas();
    
    int max=nmass;
    double high[nplots]={2,1.5,1.2,1.1,1.2,
        1.5,1.2,1.7,1.7,1.5,
        1.2,1.2,1.2,1.5,1.5,
        1.2,1.2,1.2,1.5,1.5,
        1.5,1.2,1.7,1.7,1.5,
        1.2,1.2,1.2,1.5,1.5,
        1.2,1.2,1.2,1.5,1.5}; 
    PT_HAT->cd();
    vector<double> norm; for(int j=1;j<maxtodo;j++) norm.push_back(1./plots[0][13].Integral());
        double fixnorm = 1./100000.;
        for(int i=0;i<nplots;i++) {
            //if(i==16 || i==4 || i==5 || i==6 || i==7  || i==12) PT_HAT->SetLogy(1); else 
            PT_HAT->SetLogy(0);
            for(int j=0;j<maxtodo;j++) {
                leg->AddEntry(plots[j][i],lege[todo[j]],"l");
            }
            //double normalize0 = 1./plots[0][13].Integral();
            //plots[0][i].Scale(norm[0]);
            //plots[0][i].Scale(fixnorm);
            //plots[0][i].SetMaximum(high[i]*plots[0][i].GetMaximum());
            plots[0][i].Draw("Hist");
            for(int j=0;j<maxtodo;j++) {
                //plots[j][i].Scale(1./10000.);
                //else
                double norma; if(plots[j][i].Integral()>0) norma = 1./plots[j][i].Integral(); else norma=1;
                if (i!=13) plots[j][i].Scale(norma);
                //plots[j][i].SetMaximum(high[i]*plots[j][i].GetMaximum());
                plots[j][i].Draw("Hist,same");
            }
            if(i==0) {TLine li(4,0.00001,4,0.25); li->Draw("same");}
            if(i==5) {TLine li3(173,0.0002,173,0.6); li3->Draw("same");}
            if(i==9) {TLine li2(80,0.00001,80,0.70); li2->Draw("same");}
            if(i==25 || i==26) {TLine li3(203,0.2,203,10000); li3->Draw("same");}
            if(i==25 || i==26) {TLine li2(143,0.2,143,10000); li2->Draw("same");}
            leg->Draw("same");
            //	if( i==1 || i==5) {
            //          PT_HAT->SetLogy(1); 
            //TLine li(500,0.00001,500,0.001);
            //li->Draw("same");
            //       } else PT_HAT->SetLogy(0); 	
            PT_HAT->SaveAs(namplots[i]);
            //	if(i==28 || i==22 || i==9 || i==1 || i==0 || i==29) PT_HAT->SaveAs(namplotspdf[i]);
            PT_HAT->Clear();
            leg->Clear();
        }
    PT_HAT->Close(); 
    ////////////////////////////////////////////////////////////////////////////
    TMultiGraph *mg1 = new TMultiGraph();
    //mg1->SetMaximum(0.6);
    //mg->GetXaxis()->SetRangeUser(490,520);
    
    int nevents=10000;//20000;
    vector<double> cat0,cat1,cat2,ctot;
    for(int j=0;j<maxtodo;j++) {
        cout<<lege[todo[j]]<<endl;
        // get number of njets
        int nbins = plots[j][13].GetNbinsX();
        //cout<<nbins<<endl;
        for(int k=0;k<nbins;k++) {
            double njets = plots[j][13].GetBinContent(k); 
            if (k == 0)      {cat0.push_back(njets); cout<<"zero "<<cat0[j]<<endl;}
            else if (k == 1) {cat1.push_back(njets); cout<<"mistag "<<cat1[j]<<endl;}
            else if (k == 2) {cat2.push_back(njets); cout<<"true "<<cat2[j]<<endl;}
        } // close for bins
        //
        ctot.push_back(cat0[j]+cat1[j]+cat2[j]);
        cout<<"total "<<ctot[j]/(ctot[0]+ctot[1]+ctot[2])<<endl;
        cout<<" "<<endl;
    } // close for masses
    for(int j=0;j<maxtodo;j++) cout<<lege[todo[j]]<<" "<<ctot[j]/nevents<<" "<<cat1[j]/nevents<<" "<<cat2[j]/nevents<<endl;
    
    
}// end file


