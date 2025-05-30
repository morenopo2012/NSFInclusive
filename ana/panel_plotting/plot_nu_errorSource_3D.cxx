
//#include "myPlotStyle.h"
#include "GridCanvas.h"
#include "PlotUtils/MnvH2D.h"
#include "PlotUtils/HyperDimLinearizer.h"//THIS HAS TO CHANGE TO BE INCLUDED IN THE MAKE FILE EVENTUALLY.
// #include "../util/plot/plot.h"


#include "TCanvas.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TStopwatch.h"
#include "TEnv.h"
#include "TChain.h"
#include "TF2.h"
#include "Math/DistFunc.h"
#include "TGraphErrors.h"
#include "TStyle.h"
#include "TFile.h"
//#include "localColor.h"
#include "PlotUtils/MnvColors.h"
#include "Cintex/Cintex.h"

#include "myPlotStyle.h"

#include "plot.h"

using namespace PlotUtils;
void makePlots(string location,string errorsource, bool doRatio)
{

  ROOT::Cintex::Cintex::Enable();
  myPlotStyle();
  TH1::AddDirectory(false);
  TH1::SetDefaultSumw2();
  gStyle->SetErrorX(0);
  gStyle->SetEndErrorSize(2);
  
  TFile f1(location.c_str());

  MnvH2D* dataMnv=(MnvH2D*)f1.Get("h_pzptrec_data_nobck_unfold_effcor_cross_section");
  vector<string> vertnames = dataMnv->GetVertErrorBandNames();
  vector<string> latnames = dataMnv->GetLatErrorBandNames();
  MnvH2D* uup,*udown;
  if(find(vertnames.begin(),vertnames.end(),errorsource)!=vertnames.end()){
    MnvVertErrorBand2D *err = dataMnv->GetVertErrorBand(errorsource);
    udown = new MnvH2D(*err->GetHist(0));
    uup = new MnvH2D(*err->GetHist(1));
  }
  else if(find(latnames.begin(),latnames.end(),errorsource)!=latnames.end()){
    MnvLatErrorBand2D *err = dataMnv->GetLatErrorBand(errorsource);
    udown = new MnvH2D(*err->GetHist(0));
    uup = new MnvH2D(*err->GetHist(1));
  }
  else{
    cout << "You gave an error band that doesn't exist in this container? Try again." << endl;
    return;
  }
  vector<double> recoil3Dbins;
  vector<double> pt3Dbins;
  vector<double> pz3Dbins;

  //Finer version
  pt3Dbins.push_back(0.0);
  pt3Dbins.push_back(0.075);//added
  pt3Dbins.push_back(0.15);
  pt3Dbins.push_back(0.25);
  pt3Dbins.push_back(0.325);//added
  pt3Dbins.push_back(0.4);
  pt3Dbins.push_back(0.475);//added
  pt3Dbins.push_back(0.55);
  pt3Dbins.push_back(0.7);
  pt3Dbins.push_back(0.85);
  pt3Dbins.push_back(1.0);
  //  pt3Dbins.push_back(1.25);
  //  pt3Dbins.push_back(1.5);
  pt3Dbins.push_back(2.5);


  pz3Dbins.push_back(1.5);
  pz3Dbins.push_back(3.5);//added ME
  //These were added to fix unfolding
  pz3Dbins.push_back(4.5);//fix unf
  pz3Dbins.push_back(7.0);//fix unf
  //
  pz3Dbins.push_back(8.0);
  //These are added to fix unfolding
  pz3Dbins.push_back(10.0);//fix unf
  pz3Dbins.push_back(20.0);

  recoil3Dbins.push_back(0.0);
  recoil3Dbins.push_back(20.0);
  for(int i=0;i<4;i++)recoil3Dbins.push_back(i*40+40);//40,80,120,160
  for(int i=0;i<2;i++)recoil3Dbins.push_back(i*80+240);//240,320
  for(int i=0;i<2;i++)recoil3Dbins.push_back(i*200+400);
  recoil3Dbins.push_back(799.0);
  for(int i=0;i<recoil3Dbins.size();i++) recoil3Dbins[i]=recoil3Dbins[i]/1000.;

  std::vector<std::vector<double> > full3D;
  full3D.push_back(recoil3Dbins);
  full3D.push_back(pt3Dbins);
  full3D.push_back(pz3Dbins);
  
  HyperDimLinearizer *my3d = new HyperDimLinearizer(full3D,0);

  std::cout << "Starting up getting the projections" << std::endl;
  std::vector<TH2D*> dataresults = my3d->Get2DHistos(dataMnv,false);
  std::vector<TH2D*> universeUp = my3d->Get2DHistos(uup,false);
  std::vector<TH2D*> universeDown = my3d->Get2DHistos(udown,false);
  



  for(unsigned int i=1;i<pz3Dbins.size();i++){
    double pzwidth = pz3Dbins[i]-pz3Dbins[i-1];
    dataresults[i]->Scale(1e39/pzwidth,"width");
    universeUp[i]->Scale(1e39/pzwidth,"width");
    universeDown[i]->Scale(1e39/pzwidth,"width");
  }

    


  vector<int> mycolors = MnvColors::GetColors(9);
  vector<int> mycolors2 = MnvColors::GetColors(6);
  for(unsigned int i=1;i<pz3Dbins.size();i++){//skip underflow and overflow pz bins
    cout << "DOING BIN " << i << endl;

    // Make a list of the histograms we want to draw, along with the
    // draw options we want to use for them. You can add "graph" to the
    // draw options if you want the histogram to be converted to a graph
    // and then drawn. In that case the draw options are interpreted as
    // options to TGraphErrors::Draw().
    //
    // I don't know what happens if you put a "graph" first in the list,
    // so don't do that. Make sure the first item doesn't have "graph"
    // in its options

    dataresults[i]->SetLineColor(1);
    universeUp[i]->SetLineColor(2);
    universeDown[i]->SetLineColor(4);
    universeUp[i]->SetLineStyle(1);
    universeDown[i]->SetLineStyle(1);
    
    if(doRatio){
      universeUp[i]->Divide(dataresults[i]);
      universeDown[i]->Divide(dataresults[i]);
      dataresults[i]->Divide(dataresults[i]);
    }
    
    
    std::vector<std::pair<TH2*, const char*> > histAndOpts;
    if(doRatio)histAndOpts.push_back(std::make_pair(dataresults[i], "hist"));
    else histAndOpts.push_back(std::make_pair(dataresults[i], "histPE"));
    histAndOpts.push_back(std::make_pair(universeUp[i], "hist"));
    histAndOpts.push_back(std::make_pair(universeDown[i], "hist"));
    TLegend* leg=new TLegend(0.8, 0.1, 1, 0.3);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->SetTextSize(0.03);
    leg->AddEntry(dataresults[i], "MINERvA data CV", "l");
    leg->AddEntry(universeUp[i], Form("+1#sigma %s",errorsource.c_str()),"l");
    leg->AddEntry(universeDown[i], Form("-1#sigma %s",errorsource.c_str()),"l");
    // // ----------------------------------------------------------------------------------
    // //
    // // Now make pz in bins of pt. It's all the same

    // // Values to multiply each bin by to get them on a similar range
    vector<double> multiplierspz1 = GetScales(histAndOpts,true,5.49,0.75,false);
    
    // // plotXAxis1D fiddles the x axis values to squash up the tail so it
    // // doesn't take up all the horizontal space.
    GridCanvas* gc2= NULL;

    gc2=plotXAxis1D(histAndOpts, "#Sigma T_{p} (GeV)", "P_{t,muon}", doRatio?NULL:&multiplierspz1[0]);
    if(doRatio){
      gc2->SetYLimits(0.81,1.19);
      gc2->SetYTitle("Ratio to CV value");
    }
    else{
      gc2->SetYLimits(0.01,5.49);
      gc2->SetYTitle("d^{3}#sigma/dp_{T}dp_{||}d#SigmaT_{p} (x10^{-39} cm^{2}/GeV^{3}/c^{2}/Nucleon)");
    }
    gc2->Modified();
    TLatex mytex;
    mytex.SetTextSize(0.05);
    string mystring =     Form("%.2f < P_{||} [GeV/c] < %.2f",pz3Dbins[i-1],pz3Dbins[i]);

    mytex.DrawLatex(0.35,0.96,mystring.c_str());
    leg->SetX1(0.79);
    leg->SetY1(0.1);
    leg->SetX2(1.0);
    leg->SetY2(0.4);
    leg->Draw("SAME");
    if(doRatio){
      gc2->Print(Form("nu-3d-errorsource-plot-pz_bin_%d_Ratio.eps",i));
      gc2->Print(Form("nu-3d-errorsource-plot-pz_bin_%d_Ratio.png",i));
      gc2->Print(Form("nu-3d-errorsource-plot-pz_bin_%d_Ratio.C",i));
    }
    else{
      gc2->Print(Form("nu-3d-errorsource-plot-pz_bin_%d.eps",i));
      gc2->Print(Form("nu-3d-errorsource-plot-pz_bin_%d.png",i));
      gc2->Print(Form("nu-3d-errorsource-plot-pz_bin_%d.C",i));
    }
  }
}

int main(int argc, char* argv[])
{
  makePlots(argv[1],argv[2],false);
  makePlots(argv[1],argv[2],true);
  
  return 0;
}
