// main11.cc is a part of the PYTHIA event generator.
// Copyright (C) 2014 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.
// This is a simple test program.
// It illustrates how Les Houches Event File input can be used in Pythia8.
// It uses the ttsample.lhe input file, the latter only with 100 events.
#include "Pythia8/Pythia.h"
using namespace Pythia8;
using namespace std;

int main() {

    // Generator. We here stick with default values, but changes
  // could be inserted with readString or readFile.
  for(unsigned int fol = 0; fol<4;fol++)
   for(unsigned int gam = 7; gam<13; gam++){  
      Pythia pythia;
       //unsigned int fol = 0; unsigned int gam = 0;
       // Initialize Les Houches Event File run. List initialization information.
       pythia.readString("Beams:frameType = 4");
       string LHEInput = "Beams:LHEF = ";
       string path = "/afs/cern.ch/work/a/acarvalh/tt_fulylep/";
       string folder, files;
      if(fol==0) folder = "OnOnVary/";
      if(fol==1) folder = "OnOffVary/";
      if(fol==2) folder = "OffOnVary/";
      if(fol==3) folder = "OffOffVary/";
       if(fol==4) folder = "FullVary/";      
      ////////////////////////////////
      if(gam==0) files = "Wt_0.lhe"; 
      if(gam==1) files = "Wt_1.lhe"; 
      if(gam==2) files = "Wt_2.lhe"; 
      if(gam==3) files = "Wt_3.lhe"; 
      if(gam==4) files = "Wt_4.lhe"; 
      if(gam==5) files = "Wt_5.lhe"; 
      if(gam==6) files = "Wt_6.lhe"; 
      if(gam==7) files = "Wt_7.lhe"; 
      if(gam==8) files = "Wt_8.lhe"; 
      if(gam==9) files = "Wt_9.lhe"; 
      if(gam==10) files = "Wt_10.lhe"; 
      if(gam==11) files = "Wt_11.lhe"; 
      if(gam==12) files = "Wt_12.lhe"; 
      ///////////////////////////////
      cout<<LHEInput+folder+files<<endl;
      pythia.readString(LHEInput+path+folder+files);
      pythia.readString("PartonLevel:all = on"); // Of hadronization
      pythia.readString("HadronLevel:all = off"); // Of hadronization
      //
      pythia.readString("5:mayDecay = no");
      pythia.readString("-5:mayDecay = no");  
       //cout<<namefile_out<<endl;  
    string namefile_out=path+folder+files+".shower";
    ofstream out_pythia;
    out_pythia.precision(5); 
    out_pythia.open(namefile_out.c_str());
       //
       string namefile_outME=path+folder+files+".decayed";       
       ofstream out_pythiaME;
       out_pythiaME.precision(5); 
       out_pythiaME.open(namefile_outME.c_str());
       //
      // Allow for possibility of a few faulty events.
      int nAbort = 100;
      int iAbort = 0;
      // Begin event loop; generate until none left in input file.
      pythia.init();
    for (int iEvent = 0; ; ++iEvent) {// iEvent<10000
    // Generate events, and check whether generation failed.
    if (!pythia.next()) {

      // If failure because reached end of file then exit event loop.
      if (pythia.info.atEndOfFile()) break;

      // First few failures write off as "acceptable" errors, then quit.
      if (++iAbort < nAbort) continue;
      break;
     }
      //cout<<"Number of particles = "<<pythia.process.size()<<endl;
      vector<int> pID , pIDgen;
      vector<double> px , pxgen;
      vector<double> py , pygen;
      vector<double> pz , pzgen;
      vector<double> E , Egen;
      vector<int> mother , mothergen;
      vector<int> code , codegen;
      // Some checks on the event record
      // Check for example that at least we have two bs and two bbars
      for (int i = 0; i < pythia.event.size(); i++){
          int particle_id = pythia.event[i].id();
          int particle_status = pythia.event[i].status();
          int particle_mother = pythia.event[i].mother1();
          // save only final state particles
          if(particle_status>0){
              //cout<<i<<" "<<particle_id<<" "<<particle_mother<<endl;
              double ppx= pythia.event[i].px();
              double ppy= pythia.event[i].py();
              double ppz= pythia.event[i].pz();
              double EE= pythia.event[i].e();
              //cout<<ppx<<" "<<ppy<<" "<<ppz<<" "<<EE<<endl;
              pID.push_back(particle_id);
              px.push_back(ppx);
              py.push_back(ppy);
              pz.push_back(ppz);
              E.push_back(EE);
              mother.push_back(particle_mother);
              code.push_back(particle_id);
         } // close particle status
        } // close event size 
        ////////////////////////////////////////////////////////////////////////
      } // close process size
      //////////////////////////////////////////////////////////////////
      // Save into file
        //cout<<"Number of final state particles = "<<E.size()<<"\n"<<endl;
        out_pythia.flush();
        /////////////////////////////////////////////////////////////////////
        out_pythia<<"#"<<endl;
        out_pythia<<E.size()<<endl;
        for(unsigned i=0;i<E.size();i++) {out_pythia<<pID.at(i)<<" "<<px.at(i)<<" "<<py.at(i)<<" "<<pz.at(i)<<" "<<E.at(i)<<" "<<endl;}
        ///////////////////////////////////////////////////////////////////////
     	/////////////////////////////////////////////////////////////////////// <<mothergen.at(i)<<" "
  // End of event loop.
  }
  out_pythia.close();
  cout<<namefile_out<<endl;
  // Give statistics. Print histogram.
  pythia.stat();
  //cout << nCharged;
  } // close for folder
  // Done.
  return 0;
}
