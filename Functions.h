// Test function //
void hello();
//////////////////////////////////////////////////////////
// histos
int decla(int);
int save_hist(int,int,int);
//////////////////////////////////////////////////////////
int truetops(vector<PseudoJet> jets, vector<PseudoJet> leptons,vector<PseudoJet> neutrinos, vector<int> btag, vector<int> btrue);
//////////////////////////////////////////////////////////
// tags
int recojets(vector<PseudoJet> particles,vector<PseudoJet> & jets, vector<int> & btag, 
             vector<int> & bmistag, vector<int> & fattag, vector<int> & btrue);
int recol(vector<PseudoJet> jets,vector<PseudoJet> leptons,vector<PseudoJet> neutrinos);
bool recohadt(int & bh, int & bl, vector<PseudoJet> jets,vector<PseudoJet> leptons,vector<PseudoJet> neutrinos, 
              vector<int> btag, vector<int> btrue, double met);
bool recolept2step(int & bh, int & bl,vector<PseudoJet> jets, vector<PseudoJet> leptons,vector<PseudoJet> neutrinos, 
              vector<int> btag, vector<int> btrue, double met);
bool recotlepeq(int & bh, int & bl,vector<PseudoJet> jets, vector<PseudoJet> leptons,vector<PseudoJet> neutrinos, 
              vector<int> btag, vector<int> btrue, double met);

bool fullylep(int & bh, int & bl, vector<PseudoJet> jets,vector<PseudoJet> leptons,vector<PseudoJet> neutrinos, 
              vector<int> btag, vector<int> btrue, double met);
