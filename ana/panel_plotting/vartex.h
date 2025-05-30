// Using var pair

std::map<string, string> varLatex;
std::map<string, string> varTexX;
//dalphat dphit dpt dptx dpty pn tp ptheta
void InitVarLatex() {
string oname[15] = {
  "dalphat",
  "dphit",
  "dPPi", 
  "dPP",
  "dPRi",
  "dPR",
  "dpt",
  "dptx",
  "dpl",
  "dpty",
  "dthetaP", 
  "pn", 
  "ptmu",
  "tp",
  "ptheta"
};


string llname[15] = {
  "d#alpha_{t}",
  "d#phi_{t}",
  "dP_{p,i}",
  "dP_{p}",
  "dP_{r,i}",
  "dP_{r}",
  "dP_{t}",
  "dP_{t,x}",
  "dP_{l}",
  "dP_{t,y}",
  "d#theta_{p}",
  "p_{n}",
  "p_{t#mu}",
  "T_{p}",
  "p_{#theta}"
};


string llnameX[15] = {
  "d#alpha_{t}",
  "d#phi_{t}",
  "dP_{p,i}",
  "dP_{p}",
  "dP_{r,i}",
  "dP_{r}",
  "dP_{t}",
  "dP_{t,x}",
  "dP_{l}",
  "dP_{t,y}",
  "d#theta_{p}",
  "p_{n} (GeV)",
  "p_{t#mu} (GeV)",
  "T_{p}",
  "p_{#theta}"
};


 for(int i = 0; i < 14; i++) {
   varLatex.insert(pair<string, string>(oname[i],llname[i]));
   varTexX.insert(pair<string, string>(oname[i],llnameX[i]));
 }
}

