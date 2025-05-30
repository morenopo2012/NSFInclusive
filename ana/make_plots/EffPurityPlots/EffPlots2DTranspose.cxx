#include "include/NukeCCUtilsNSF.h"
#include "TParameter.h"
#include "../drawUtils.h"
#include "TColor.h"

void makeplot(string histsuffix, TFile* myfile, bool doprelimlabel=false) {
  string varnames = "Eav_vs_q3";

  // Load histograms
  PlotUtils::MnvH2D* num = (PlotUtils::MnvH2D*)myfile->Get("h_mc_Eavailable_q3");
  PlotUtils::MnvH2D* dem = (PlotUtils::MnvH2D*)myfile->Get("h_truth_Eavailable_q3");

  // Divide to get efficiency
  num->Divide(num, dem);

  // Transpose helper
  auto TransposeHist = [](PlotUtils::MnvH2D* hist) -> PlotUtils::MnvH2D* {
    int nx = hist->GetNbinsX();
    int ny = hist->GetNbinsY();
    const TArrayD* xbins = hist->GetXaxis()->GetXbins();
    const TArrayD* ybins = hist->GetYaxis()->GetXbins();

    PlotUtils::MnvH2D* transposed;
    if (xbins->fN && ybins->fN)
      transposed = new PlotUtils::MnvH2D("transposed", "transposed", ny, ybins->GetArray(), nx, xbins->GetArray());
    else
      transposed = new PlotUtils::MnvH2D("transposed", "transposed", ny, hist->GetYaxis()->GetXmin(), hist->GetYaxis()->GetXmax(),
                                        nx, hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmax());

    for (int ix = 1; ix <= nx; ++ix) {
      for (int iy = 1; iy <= ny; ++iy) {
        double content = hist->GetBinContent(ix, iy);
        double error = hist->GetBinError(ix, iy);
        transposed->SetBinContent(iy, ix, content);
        transposed->SetBinError(iy, ix, error);
      }
    }

    return transposed;
  };

  PlotUtils::MnvH2D* ratio = TransposeHist(num);

  // Style and labels
  ratio->GetXaxis()->SetTitle("q_{3} (GeV)");
  ratio->GetYaxis()->SetTitle("E_{av} (GeV)");
  ratio->GetZaxis()->SetTitle("Efficiency");

  applyStyle(ratio);
  gStyle->SetStripDecimals(0);
  ratio->GetXaxis()->SetTitleOffset(1.3);
  ratio->GetZaxis()->SetTitleOffset(1.5);
  ratio->GetZaxis()->SetRangeUser(0, 1);

  // Canvas and drawing
  TCanvas c2("test", "test");
  c2.SetRightMargin(0.1788009);
  gStyle->SetPalette(112);
  ratio->Draw("COLZ");

  if (doprelimlabel) {
    TLatex* labl = new TLatex(2.1, 1.3, "MINER#nuA Preliminary");
    labl->SetTextSize(0.03);
    labl->SetTextFont(112);
    labl->SetTextColor(2);
    labl->Draw();
  }

  c2.Print(Form("%s_%s_Eff_Transposed.C", varnames.c_str(), histsuffix.c_str()));
  c2.Print(Form("%s_%s_Eff_Transposed.eps", varnames.c_str(), histsuffix.c_str()));
  c2.Print(Form("%s_%s_Eff_Transposed.png", varnames.c_str(), histsuffix.c_str()));
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
  //ROOT::Cintex::Cintex::Enable();
  TH1::AddDirectory(false);


  string location=argv[1];
  int targetID=atoi(argv[2]);
  int targetZ=atoi(argv[3]);
  bool doprelimlabel=argv[4];
  //TFile *infile = new TFile("/exp/minerva/app/users/omorenop/cmtuser/git-Mat/NSFNukeCCInclusive/ana/make_hists/Migration/Efficiency_1D_minervame1L_t99_z99_nosys.root");
  TFile *infile = new TFile("/exp/minerva/app/users/omorenop/cmtuser/git-Mat/NSFNukeCCInclusive/ana/make_hists/Migration/Efficiency_1D_minervame1G_t95_z95_nosys.root");

  //makeplot("cc",infile, doprelimlabel);
  //makeplot("qe",infile, doprelimlabel);
  //makeplot("res",infile, doprelimlabel);
  //makeplot("dis",infile, doprelimlabel);
  //makeplot("true_dis",infile, doprelimlabel);
  //makeplot("sis",infile, doprelimlabel);
  //makeplot("2p2h",infile, doprelimlabel);
  //makeplot("oth",infile, doprelimlabel);
  makeplot("new",infile,doprelimlabel);

  return 0;
}

