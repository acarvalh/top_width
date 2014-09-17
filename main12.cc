// main12.cc is a part of the PYTHIA event generator.
// Copyright (C) 2011 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This is a simple test program. 
// It illustrates how Les Houches Event File input can be used in Pythia8.
// It uses the ttsample.lhe input file, the latter only with 100 events.
// Other samples could be generated as illustrated by main53.f.

#include "Pythia.h"
//#include <fstream>
using namespace Pythia8;
int main() {
  Pythia pythia;                            

  // Stick with default values, so do not bother with a separate file
  // for changes. However, do one change, to show readString in action.
  pythia.readString("Beams:frameType = 4");
  // the analysis program
  //pythia.readString("SLHA:readFrom = 2");
  //pythia.readString("SLHA:file = Susy.txt "); // input the decay table
  pythia.readString("PartonLevel:MI = off"); // Off multiple interactions
  pythia.readString("PartonLevel:ISR = off"); // Shower on
  pythia.readString("PartonLevel:FSR = off"); // Shower on
  pythia.readString("HadronLevel:all = off"); // Of hadronization 
  //
  // select semileptonic events t>w+b 
  pythia.readString("6:onMode = off");
  pythia.readString("6:onIfMatch = 24 5"); // mu numu
  pythia.readString("24:onMode = off");
  pythia.readString("24:onIfMatch = 12 11"); // e ve
  pythia.readString("24:onIfMatch = 14 13"); // mu numu
  pythia.readString("-24:onMode = off");
  pythia.readString("-24:onIfMatch = 12 -11"); // e ve
  pythia.readString("-24:onIfMatch = 14 -13"); // mu numu 
  string path;
  //path = "/afs/cern.ch/user/a/acarvalh/ttbar_lo/Events/14tev/unweighted_events.lhe";
  path = "/afs/cern.ch/work/a/acarvalh/ttbar_end_lo/Events/14tev/unweighted_events.lhe";
  //string namefile_in=path;
  //string sfile = "Beams:LHEF ="+namefile_in;
  //pythia.readString(sfile.c_str());
  // Initialize Les Houches Event File run. List initialization information.
  pythia.init(path);      
  string namefile_out;
  namefile_out = "/afs/cern.ch/user/a/acarvalh/ttbar_lo/Events/14tev/unweighted_events.lhe.decayed";
  ofstream out_pythia;
  out_pythia.precision(3);
  out_pythia.open(namefile_out.c_str());
  // Allow for possibility of a few faulty events.
  int nAbort = 10;
  int iAbort = 0;
  ////////////////
  // Begin event loop; generate until none left in input file.     
  for (int iEvent = 0;iEvent<10 ; ++iEvent) {
    cout<<"\n ievent = "<<iEvent<<"\n"<<endl;
    // Generate events, and check whether generation failed.
    if (!pythia.next()) {
      if (pythia.info.atEndOfFile()) break; 
      if (++iAbort < nAbort) continue;
      break;
    }
    // List first few events: Les Houches, hard process and complete.
    //if (iEvent < nPrint) {     
    //  pythia.LHAeventList();               
    //  pythia.info.list();          
    //  pythia.process.list();          
    //  pythia.event.list();           
    //}                           
    cout<<"Number of particles showered = "<<pythia.process.size()<<endl;
    vector<int> pID;
    vector<double> px;
    vector<double> py;
    vector<double> pz;
    vector<double> E;
    vector<int> mother;
    vector<int> code;
    int counter=0;
    for (int i = 0; i < pythia.process.size(); i++){
      int particle_id = pythia.process[i].id();
      int particle_status = pythia.process[i].status(); // change for "event" in hadron level 
      int particle_mother = pythia.process[i].mother1();
      // save only final state particles
      if( 1>0
         && (particle_status>0 )//|| abs(particle_id)==6 || abs(particle_id)==24)
         //&& particle_id!=2101 && particle_id!=2103 // remove remnants 
        ) { 
        cout<<i<<" "<<particle_id<<" "<<particle_mother<<" "<<particle_status<<endl;
        double ppx= pythia.process[i].px();
        double ppy= pythia.process[i].py();
        double ppz= pythia.process[i].pz();
        double EE= pythia.process[i].e();
        //cout<<px<<" "<<py<<" "<<pz<<" "<<E<<endl;
        pID.push_back(particle_id);
        px.push_back(ppx);
        py.push_back(ppy);
        pz.push_back(ppz);
        E.push_back(EE);
        mother.push_back(particle_mother);
        code.push_back(particle_id);
	if(particle_status>0) counter++;
      }
    }
    // Save into file
    out_pythia<<"#"<<endl;
    cout<<"Number of final state particles = "<<counter<<"\n"<<endl;
    out_pythia<<E.size()<<endl;
     for(unsigned i=0;i<E.size();i++){
       out_pythia<<pID.at(i)<<" "<<px.at(i)<<" "<<py.at(i)<<" "<<pz.at(i)<<" "<<E.at(i)<<" "<<endl;
     }    
  // End of event loop.
  }                                         
  out_pythia.close();
  // Give statistics. Print histogram.
  pythia.statistics();
  // Done.                           
  return 0;
}
