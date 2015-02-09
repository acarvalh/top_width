// Test function
void hello();
//////////////////////////////////////////////////////////
// histos
int decla(int);
int save_hist(int,int,int,bool);
//////////////////////////////////////////////////////////
int truetops(vector<PseudoJet> jets, vector<PseudoJet> leptons,vector<PseudoJet> neutrinos, vector<int> btag, vector<int> btrue);
//////////////////////////////////////////////////////////
// tags
bool GenLevelDilep(vector<PseudoJet> particles,  vector<PseudoJet> leptons, vector<PseudoJet> neutrinos , double cut, double weight, int type);
int recojets(vector<PseudoJet> particles,vector<PseudoJet> & jets, vector<int> & btag, 
             vector<int> & bmistag, vector<int> & fattag, vector<int> & btrue, double weight);
int recol(vector<PseudoJet> jets,vector<PseudoJet> leptons,vector<PseudoJet> neutrinos, double weight);
bool recohadt(int & bh, int & bl, vector<PseudoJet> jets,vector<PseudoJet> leptons,vector<PseudoJet> neutrinos, 
              vector<int> btag, vector<int> btrue, double met, double weight, double cut, int type);
bool recolept2step(int & bh, int & bl,vector<PseudoJet> jets, vector<PseudoJet> leptons,vector<PseudoJet> neutrinos, 
                   vector<int> btag, vector<int> btrue, double met, double weight, double cut, int type);
bool recotlepeq(int & bh, int & bl,vector<PseudoJet> jets, vector<PseudoJet> leptons,vector<PseudoJet> neutrinos, 
                vector<int> btag, vector<int> btrue, double met, double weight, double cut, int type);

bool fullylep(int & bh, int & bl, vector<PseudoJet> jets,vector<PseudoJet> leptons,vector<PseudoJet> neutrinos, 
              vector<int> btag, vector<int> btrue, double met,double weight);
