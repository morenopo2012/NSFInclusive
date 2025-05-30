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

void makePlots(string location)
{

  ROOT::Cintex::Cintex::Enable();
  myPlotStyle();
  TH1::AddDirectory(false);
  TH1::SetDefaultSumw2();
  gStyle->SetErrorX(0);
  gStyle->SetEndErrorSize(2);


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

  /*
  recoil3Dbins.push_back(0.0);
  recoil3Dbins.push_back(0.0200);
  for(int i=0;i<9;i++)recoil3Dbins.push_back(i*0.040+0.040);
  for(int i=0;i<2;i++)recoil3Dbins.push_back(i*0.2+0.4);
  recoil3Dbins.push_back(0.7999999);
  */
  recoil3Dbins.push_back(0.0);
  recoil3Dbins.push_back(20.0);
  //for(int i=0;i<10;i++)recoil3Dbins.push_back(i*40); 
  //for(int i=0;i<9;i++)recoil3Dbins.push_back(i*40+40);//40,80,120,160,200,240,280,320,360
  for(int i=0;i<4;i++)recoil3Dbins.push_back(i*40+40);//40,80,120,160
  for(int i=0;i<2;i++)recoil3Dbins.push_back(i*80+240);//240,320
  for(int i=0;i<2;i++)recoil3Dbins.push_back(i*200+400);
  recoil3Dbins.push_back(799.0);
  for(int i=0;i<recoil3Dbins.size();i++) recoil3Dbins[i]=recoil3Dbins[i]/1000.;//make GeV



  std::vector<std::vector<double> > full3D;
  full3D.push_back(recoil3Dbins);
  full3D.push_back(pt3Dbins);
  full3D.push_back(pz3Dbins);

  
  HyperDimLinearizer *my3d = new HyperDimLinearizer(full3D,0);




  //three files 1-track, 2+track, N-track
  //CV
  TFile f1(Form("%s_CV_Zexp_Bubble/EffPurity_MakeFlux-1_CombinedPlaylists.root",location.c_str()));//1track
  
  //numerator
  MnvH2D *mnvqelike = (MnvH2D*)f1.Get("h_pzptrec_qelike");
  MnvH2D *mnvqelike_qe = (MnvH2D*)f1.Get("h_pzptrec_qelike_qe");
  MnvH2D *mnvqelike_res = (MnvH2D*)f1.Get("h_pzptrec_qelike_res");
  MnvH2D *mnvqelike_2p2h = (MnvH2D*)f1.Get("h_pzptrec_qelike_2p2h");
  MnvH2D *mnvqelike_dis = (MnvH2D*)f1.Get("h_pzptrec_qelike_dis");
  //denominator
  MnvH2D *mnv_den_qelike = (MnvH2D*)f1.Get("h_pzptrec_truth_qelike");
  MnvH2D *mnv_den_qelike_qe = (MnvH2D*)f1.Get("h_pzptrec_truth_qelike_qe");
  MnvH2D *mnv_den_qelike_res = (MnvH2D*)f1.Get("h_pzptrec_truth_qelike_res");
  MnvH2D *mnv_den_qelike_2p2h = (MnvH2D*)f1.Get("h_pzptrec_truth_qelike_2p2h");
  MnvH2D *mnv_den_qelike_dis = (MnvH2D*)f1.Get("h_pzptrec_truth_qelike_dis");



  mnvqelike->Divide(mnvqelike,mnv_den_qelike);
  mnvqelike_qe->Divide(mnvqelike_qe,mnv_den_qelike_qe);
  mnvqelike_res->Divide(mnvqelike_res,mnv_den_qelike_res);
  mnvqelike_2p2h->Divide(mnvqelike_2p2h,mnv_den_qelike_2p2h);
  mnvqelike_dis->Divide(mnvqelike_dis,mnv_den_qelike_dis);


  std::vector<TH2D*> mc_qelike = my3d->Get2DHistos(mnvqelike,false);
  std::vector<TH2D*> mc_qelike_qe = my3d->Get2DHistos(mnvqelike_qe,false);
  std::vector<TH2D*> mc_qelike_res = my3d->Get2DHistos(mnvqelike_res,false);
  std::vector<TH2D*> mc_qelike_2p2h = my3d->Get2DHistos(mnvqelike_2p2h,false);
  std::vector<TH2D*> mc_qelike_dis = my3d->Get2DHistos(mnvqelike_dis,false);

  //  vector<int> mycolors = getColors(2);
  vector<int> mycolors = MnvColors::GetColors(9);
  vector<int> mycolors2 = MnvColors::GetColors(6);

  for(int i=1;i<pz3Dbins.size();i++){

    //need to add signal and bkg colors
    mc_qelike[i]->SetLineColor(1);
    mc_qelike[i]->SetLineWidth(4);
    mc_qelike_qe[i]->SetLineColor(mycolors[2]);
    mc_qelike_res[i]->SetLineColor(mycolors[1]);
    mc_qelike_dis[i]->SetLineColor(mycolors[7]);
    mc_qelike_2p2h[i]->SetLineColor(mycolors[0]);
    // Make a list of the histograms we want to draw, along with the
    // draw options we want to use for them. You can add "graph" to the
    // draw options if you want the histogram to be converted to a graph
    // and then drawn. In that case the draw options are interpreted as
    // options to TGraphErrors::Draw().
    //
    // I don't know what happens if you put a "graph" first in the list,
    // so don't do that. Make sure the first item doesn't have "graph"
    // in its options
    std::vector<std::pair<TH2*, const char*> > histAndOpts;
    histAndOpts.push_back(std::make_pair(mc_qelike[i],       "hist"));
    histAndOpts.push_back(std::make_pair(mc_qelike_qe[i],       "hist"));
    histAndOpts.push_back(std::make_pair(mc_qelike_res[i],       "hist"));
    histAndOpts.push_back(std::make_pair(mc_qelike_2p2h[i],       "hist"));
    //    histAndOpts.push_back(std::make_pair(mc_qelike_dis[i],       "hist"));

    
    
    
    // ----------------------------------------------------------------------------------
    //
    // First make pt in bins of pz

    bool doMultipliers = false;
    // Values to multiply each bin by to get them on a similar range
    double multipliers[] ={1,1,1,1,1,1,2,2,2,2,2,5,10,20,1};
    double multipliers2[]={0.5,0.5,0.5,0.5,0.5,0.5,50,1,1,1,1,2,5,25,1};
    double multipliers3[]={10,10,10,10,10,10,10,20,20,20,20,20,100,300};
    double multipliers4[]={10,10,10,50,100,1,1,1,1,1,1,1,1,1,1};
    GridCanvas* gc=NULL;
    gc = plotYAxis1D(histAndOpts,"P_{t} (GeV)","Vis. E (GeV)", doMultipliers ? multipliers : NULL);

    // Set the y range manually. Can also use gc->Remax() to guess automatically
    gc->SetYLimits(0, 0.99);
    gc->SetYTitle("Selection Efficiency");
    gc->Modified();
    // Example of adding a legend. The co-ordinate system is NDC on the
    // entire canvas, ie (0,0) in the bottom left corner of the canvas
    // (not the individual pad), and (1,1) in the top right
    TLegend* leg=new TLegend(0.79, 0.1, 1.0, 0.4);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->SetTextSize(0.03);
    leg->AddEntry(mc_qelike[i],"QELike","l");
    leg->AddEntry(mc_qelike_qe[i],"QELike and QE","l");
    leg->AddEntry(mc_qelike_res[i],"QELike and Resonant","l");
    leg->AddEntry(mc_qelike_2p2h[i],"QELike and 2p2h","l");
  

    TLatex mytex;
    mytex.SetTextSize(0.05);
    string mystring =     Form("%.2f < P_{||} [GeV] < %.2f",pz3Dbins[i-1],pz3Dbins[i]);
    mytex.DrawLatex(0.35,0.96,mystring.c_str());    
    leg->Draw("SAME");
    gc->Print(doMultipliers ? Form("nu-2d-efficiency-pt-multiplier-bin-%d.eps",i) : Form("nu-2d-efficiency-pt-bin-%d.eps",i));
    gc->Print(doMultipliers ? Form("nu-2d-efficiency-pt-multiplier-bin-%d.png",i) : Form("nu-2d-efficiency-pt-bin-%d.png",i));
    gc->Print(doMultipliers ? Form("nu-2d-efficiency-pt-multiplier-bin-%d.C",i) : Form("nu-2d-efficiency-pt-bin-%d.C",i));
    
    // ----------------------------------------------------------------------------------
    //
    // Now make pz in bins of pt. It's all the same
    
    // Values to multiply each bin by to get them on a similar range
    double  multipliers5[]={4,1,1,1,4,200,75,1,1,1,1,1,1,1};
    double  multipliers6[]={2,0.7,0.5,0.5,1,10,5,500,1,1,1,1,1,1};
    double  multipliers7[]={25,10,10,10,20,100,5,100,1,1,1,1,1,1};
    double  multipliers8[]={100,50,25,25,25,50,100,500,1,1,1,1,1,1};
			    
    
    // plotXAxis1D fiddles the x axis values to squash up the tail so it
    // doesn't take up all the horizontal space.
    GridCanvas* gc2=NULL;
    gc2 = plotXAxis1D(histAndOpts,"Visible Energy (GeV)","P_{t} (GeV)", doMultipliers ? multipliers5 : NULL);
    
    mytex.DrawLatex(0.35,0.96,mystring.c_str());
    gc2->SetYLimits(0, 0.99);
    gc2->SetYTitle("Selection Efficiency");
    leg->Draw("SAME");

    gc2->Modified();
    gc2->Print(doMultipliers ? Form("nu-2d-efficiency-pz-multiplier-bin-%d.eps",i) : Form("nu-2d-efficiency-pz-bin-%d.eps",i));
    gc2->Print(doMultipliers ? Form("nu-2d-efficiency-pz-multiplier-bin-%d.png",i) : Form("nu-2d-efficiency-pz-bin-%d.png",i));
    gc2->Print(doMultipliers ? Form("nu-2d-efficiency-pz-multiplier-bin-%d.C",i) : Form("nu-2d-efficiency-pz-bin-%d.C",i));
  }
}

int main(int argc, char* argv[])
{

  string location = argv[1];
  makePlots(location);

  return 0;
}
