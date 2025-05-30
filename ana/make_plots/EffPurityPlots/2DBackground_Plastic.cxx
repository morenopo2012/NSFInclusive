//#include "include/CCQENuPlotUtils.h"
#include "include/NukeCCUtilsNSF.h"
#include "TParameter.h"
#include "../drawUtils.h"
#include <string>

void makeplot(string histsuffix, TFile* myfile, bool doprelimlabel=false){
//  string varnames;

         //varnames = "Enu_Ehad";

 vector<string> vars;
   vars.clear();
   vars.push_back("Emu_Ehad");
   vars.push_back("Enu_Ehad");
  //PlotUtils::MnvH2D *dem = (PlotUtils::MnvH2D*)myfile->Get(Form("%s_truth_%s",varnames.c_str(),histsuffix.c_str()));
  //TFile *f1 = TFile:T:Open("Hists_Efficiency_t1_z26_Nu_v1_.root");
 for( unsigned int i = 0; i != vars.size(); ++i ){
  const string& varnames = vars[i];
 
  PlotUtils::MnvH2D *USFe = (PlotUtils::MnvH2D*)myfile->Get(Form("selected_mc_US_Fe_%s", varnames.c_str()));
  PlotUtils::MnvH2D *USPb = (PlotUtils::MnvH2D*)myfile->Get(Form("selected_mc_US_Pb_%s", varnames.c_str()));
  PlotUtils::MnvH2D *USC = (PlotUtils::MnvH2D*)myfile->Get(Form("selected_mc_US_C_%s", varnames.c_str()));
  PlotUtils::MnvH2D *USOther = (PlotUtils::MnvH2D*)myfile->Get(Form("selected_mc_US_other_%s", varnames.c_str()));
  PlotUtils::MnvH2D *USregUS = (PlotUtils::MnvH2D*)myfile->Get(Form("selected_mc_US_regUS_%s", varnames.c_str()));
  PlotUtils::MnvH2D *USregDS = (PlotUtils::MnvH2D*)myfile->Get(Form("selected_mc_US_regDS_%s", varnames.c_str()));
                                  
  PlotUtils::MnvH2D *DSFe = (PlotUtils::MnvH2D*)myfile->Get(Form("selected_mc_DS_Fe_%s", varnames.c_str()));
  PlotUtils::MnvH2D *DSPb = (PlotUtils::MnvH2D*)myfile->Get(Form("selected_mc_DS_Pb_%s", varnames.c_str()));
  PlotUtils::MnvH2D *DSC = (PlotUtils::MnvH2D*)myfile->Get(Form("selected_mc_DS_C_%s", varnames.c_str()));
  PlotUtils::MnvH2D *DSOther = (PlotUtils::MnvH2D*)myfile->Get(Form("selected_mc_DS_other_%s", varnames.c_str()));
  PlotUtils::MnvH2D *DSregDS = (PlotUtils::MnvH2D*)myfile->Get(Form("selected_mc_DS_regDS_%s", varnames.c_str()));
  PlotUtils::MnvH2D *DSregUS = (PlotUtils::MnvH2D*)myfile->Get(Form("selected_mc_DS_regUS_%s", varnames.c_str()));

  DSFe->GetXaxis()->SetTitle("Reconstructed Muon Energy");
  DSFe->GetYaxis()->SetTitle("Reconstructed Hadronic Energy");
  
   cout<<"Print my name"<<endl;

  TCanvas c2("test","test");
  c2.SetRightMargin(0.1788009);
  DSFe->Draw("colz");

  if (doprelimlabel) {  
    TLatex* labl=new TLatex(7, 20.6,"MINER#nuA Preliminary  \n  POT: 2.52 x10^{21}" );
    labl->SetTextSize(0.03);
    labl->SetTextFont(112);//22);                                                                            
    labl->SetTextColor(kRed);
    labl->Draw();
  }
  
  //c2.Print(Form("%s_%s_Eff.C",varnames.c_str(),histsuffix.c_str()));
  //c2.Print(Form("%s_%s_Eff.eps",varnames.c_str(),histsuffix.c_str()));
  //c1.Print(Form("%s_%s_ContinInTrans.png",varnames.c_str(),histsuffix.c_str()));
  c2.Print(Form("%s_%s_USFe.png",varnames.c_str(),histsuffix.c_str()));
 
}
}

int main( int argc, char *argv[]){

  if(argc==1){
    std::cout<<"-----------------------------------------------------------------------------------------\
------"<<std::endl;
    std::cout<<"MACROS HELP:\n\n"<<
      "\t-./EffPlots2D Path_to_Output_file Target_number Material_atomic_number doPreminaryLabel\n\n"<<
       "\t-Path_to_Output_file\t =\t Path to the directory where the output ROOT file will be created \n"<\
<
      "\t-Target_number\t \t = \t Number of target you want to run over eg. 1 \n"<<
       "\t-Material_atomic_number\t =\t Atomic number of material, eg. 26 to run iron, 82 to run lead  \n"\
	     <<
      "\t-doPreliminaryLabel= Add MINERvA Preliminary to plot?"<< std::endl;
    std::cout<<"-----------------------------------------------------------------------------------------\
------"<<std::endl;
    return 0;
  }


  PlotUtils::MnvPlotter *plotter = new PlotUtils::MnvPlotter;
  plotter->SetROOT6Palette(87);
  gStyle->SetNumberContours(500);
  gStyle->SetTitleSize(1,"x");
  gStyle->SetTitleSize(1,"y");
  gStyle->SetTitleSize(1,"z");
  gStyle->SetOptStat(0);
  TH1::AddDirectory(false);


  string location=argv[1];
  int targetID=atoi(argv[2]);
  int targetZ=atoi(argv[3]);
  bool doprelimlabel=argv[4];
  TFile *infile = new TFile(Form("%s/Hists_PlasticSideband_t%d_z%d_Nu.root",location.c_str(),targetID,targetZ));
  makeplot("cc",infile, doprelimlabel);
  //makeplot("qe",infile, doprelimlabel);
  //makeplot("res",infile, doprelimlabel);
  //makeplot("dis",infile, doprelimlabel);
  //makeplot("true_dis",infile, doprelimlabel);
  //makeplot("sis",infile, doprelimlabel);
  //makeplot("2p2h",infile, doprelimlabel);
  //makeplot("oth",infile, doprelimlabel);

  return 0;
}
