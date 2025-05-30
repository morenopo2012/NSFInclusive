//#include "myPlotStyle.h"
#include "GridCanvas.h"
#include "PlotUtils/MnvH2D.h"
#include "PlotUtils/HyperDimLinearizer.h"//THIS HAS TO CHANGE TO BE INCLUDED IN THE MAKE FILE EVENTUALLY.
#include "PlotUtils/MnvColors.h"
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
#include "TLine.h"
#include "TArrow.h"
//#include "localColor.h"
#include "Cintex/Cintex.h"

#include "myPlotStyle.h"

#include "plot.h"

using namespace PlotUtils;

vector<double> GetScales(std::vector<std::pair<TH2*, const char*> >histopts,double ymax, bool pzProj){
  vector<double> tmpvect;
  int nbins = histopts[0].first->GetNbinsX()+2;
  if(pzProj) nbins = histopts[0].first->GetNbinsY()+2;
  for(int i=1;i<nbins;i++){
    double maxval = 0;
    for(uint j=0;j<histopts.size();j++){
      TH1D *tmp = pzProj? histopts[j].first->ProjectionX("tmp",i,i): histopts[j].first->ProjectionY("tmp",i,i);
      int maxbin = tmp->GetMaximumBin();
      double content = tmp->GetBinContent(maxbin);
      if(content>maxval) maxval=content;
    }
    //we want abaout 75% of the 1.49, so 1.15
    double scale = ymax*(0.9)/maxval;
    if(scale>1){
      int tmpscale = floor(scale*10);
      scale = tmpscale/10.0;
    }
    else{
      int tmpscale = ceil(scale*10);
      scale = tmpscale/10.0;
    }
    cout << scale << endl;
    tmpvect.push_back(scale);
  }
    return tmpvect;
}


void makePlots(bool doMultipliers,bool doRatio, string location, bool pzptrec)
{

  cout << "I'm running (multipliers,ratio,pzptrec)= " << doMultipliers<<","<< doRatio << "," << pzptrec << endl;
  ROOT::Cintex::Cintex::Enable();
  myPlotStyle();
  TH1::AddDirectory(false);
  TH1::SetDefaultSumw2();
  gStyle->SetErrorX(0);
  gStyle->SetEndErrorSize(2);

  string variableset = "pzptrec";
  string xaxislabel = "#Sigma T_{p} (GeV)";
  string yaxislabel = "P_{t} (GeV/c)";
  string zaxislabel = "P_{||} (GeV/c)";
  string yaxislong = "Muon Transverse Momentum (GeV/c)";
  string crosssectionlabel = "d^{3}#sigma/dp_{T}dp_{||}d#SigmaT_{p} (x10^{-39} cm^{2}/GeV^{3}/c^{2}/Nucleon)";
  
  if(!pzptrec){
    variableset = "enuproxyE";
    xaxislabel = "#Sigma T_{p} (GeV)";
    yaxislabel = "q_{0}^{QE} (GeV)";
    zaxislabel = "E_{#mu} (GeV)";
    yaxislong = yaxislabel;
    crosssectionlabel = "d^{3}#sigma/dq_{0}^{QE}dE_{#mu}d#SigmaT_{p} (x10^{-39} cm^{2}/GeV^{3}/Nucleon)";
  }


  
  TFile f1(Form("%s",location.c_str()));//final

  MnvH2D* dataMnv=(MnvH2D*)f1.Get(Form("h_%s_data_nobck_unfold_effcor_cross_section",variableset.c_str()));
  MnvH2D* mcMnv=(MnvH2D*)f1.Get(Form("h_%s_mc_nobck_unfold_effcor_cross_section",variableset.c_str()));

  dataMnv->PopVertErrorBand("GENIE_FrElas_N");
  dataMnv->PopVertErrorBand("GENIE_Theta_Delta2Npi");
  dataMnv->PopVertErrorBand("GENIE_NormCCQE");
  dataMnv->PopVertErrorBand("GENIE_MaCCQEshape");

 
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
  for(int i=0;i<recoil3Dbins.size();i++) recoil3Dbins[i]=recoil3Dbins[i]/1000.;//make GeV

  std::vector<std::vector<double> > temp_full3D;
  temp_full3D.push_back(recoil3Dbins);
  temp_full3D.push_back(pt3Dbins);
  temp_full3D.push_back(pz3Dbins);

  std::vector<std::vector<double> > temp_enuproxyE;
  temp_enuproxyE.push_back(recoil3Dbins);
  temp_enuproxyE.push_back(recoil3Dbins);
  temp_enuproxyE.push_back(pz3Dbins);

  std::vector<std::vector<double> > full3D;
  if(pzptrec) full3D = temp_full3D;
  else full3D = temp_enuproxyE;


  
  HyperDimLinearizer *my3d = new HyperDimLinearizer(full3D,0);

  std::cout << "Starting up getting the projections" << std::endl;
  std::vector<TH2D*> dataresults = my3d->Get2DHistos(dataMnv,true);
  std::vector<TH2D*> dataresults_statonly = my3d->Get2DHistos(dataMnv,false);
  std::vector<TH2D*> dataresults_statonly_simple = my3d->Get2DHistos(dataMnv,false);
  std::vector<TH2D*> mcresults = my3d->Get2DHistos(mcMnv,false);

  bool dooverall = true;
  vector<double> overallscale;
  if(!doMultipliers){
    overallscale.push_back(180);
    overallscale.push_back(10);
    overallscale.push_back(1);
    overallscale.push_back(0.2);
    overallscale.push_back(0.06);
    overallscale.push_back(0.05);
  }
  else{
    for(unsigned int i=1;i<pz3Dbins.size();i++) overallscale.push_back(1);
  }

  cout << "I found " << mcresults.size() << " histograms to work with " << endl;
  for(unsigned int i=1;i<pz3Dbins.size();i++){
      double pzwidth = pz3Dbins[i]-pz3Dbins[i-1];
      dataresults[i]->Scale(1e39/pzwidth*overallscale[i-1],"width");
      dataresults_statonly[i]->Scale(1e39/pzwidth*overallscale[i-1],"width");
      dataresults_statonly_simple[i]->Scale(1e39/pzwidth*overallscale[i-1],"width");
      mcresults[i]->Scale(1e39/pzwidth*overallscale[i-1],"width");
      cout << i << "\t" << pz3Dbins[i-1] << "\t" << pz3Dbins[i] << "\t" << dataresults[i]->Integral(1,dataresults[i]->GetNbinsX(),1,dataresults[i]->GetNbinsY(),"width") << endl;
  }

    


  vector<int> mycolors = MnvColors::GetColors(10);//Light
  vector<int> mycolors2 = MnvColors::GetColors(11);//Dark
  vector<int> index;
  index.push_back(0);
  index.push_back(1);
  index.push_back(2);
  index.push_back(3);
  index.push_back(4);
  index.push_back(5);
  index.push_back(6);
  index.push_back(7);
  index.push_back(8);
  index.push_back(9);

  std::vector<std::pair<TH2*, const char*> > histAndOpts;
  std::vector<std::pair<TH2*, const char*> > histAndOptsDataPoints;
  //  TLegend *leg = new TLegend(0.15,0.95,0.39,1);//if I do 3 separate legends
  TLegend *leg = new TLegend(0.30,0.95,0.95,1);
  TLegend *leg2 = new TLegend(0.44,0.95,0.68,1);
  TLegend *leg3 = new TLegend(0.74,0.95,0.98,1);
  TLegend *leg4 = new TLegend(0.78,0.05,0.95,0.37);

  leg->SetNColumns(6);
  leg->SetLineColor(kWhite);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.03);

  leg2->SetNColumns(1);
  leg2->SetLineColor(kWhite);
  leg2->SetFillStyle(0);
  leg2->SetTextSize(0.05);

  leg3->SetNColumns(1);
  leg3->SetLineColor(kWhite);
  leg3->SetFillStyle(0);
  leg3->SetTextSize(0.05);

  leg4->SetTextSize(0.03);
  //  leg4->SetNColumns(1);
  leg4->SetLineColor(kWhite);
  //  leg4->SetFillStyle(0);
  leg4->SetHeader("P_{||} Range, Scale Factors");
  
  for(unsigned int i=1;i<pz3Dbins.size();i++){//skip underflow and overflow pz bins
    cout << "DOING BIN " << i << endl;
    // Get the data histogram with stat error and with total error
    // separately so we can plot them both for inner and outer ticks
    // These line and marker styles will be propagated to the 1D plots
      mcresults[i]->SetLineColor(mycolors2[index[i]]);
      mcresults[i]->SetLineWidth(2);

      //      if(i==2)leg2->AddEntry(mcresults[i],Form("[ %0.1f , %0.1f )",pz3Dbins[i-1],pz3Dbins[i]),"L");
      //      if(i==3)leg3->AddEntry(mcresults[i],Form("[ %0.1f , %0.1f )",pz3Dbins[i-1],pz3Dbins[i]),"L");
      // These line and marker styles will be propagated to the 1D plots
      dataresults[i]->SetMarkerStyle(1);
      //      dataresults[i]->SetMarkerSize(0.7);
      dataresults[i]->SetLineColor(mycolors[index[i]]);
      dataresults[i]->SetLineWidth(2);

      dataresults[i]->SetFillColor(mycolors[index[i]]);
      dataresults[i]->SetFillStyle(3003);
      

      dataresults_statonly[i]->SetMarkerStyle(1);
      dataresults_statonly[i]->SetMarkerSize(1);
      dataresults_statonly[i]->SetLineColor(mycolors2[index[i]]);
      dataresults_statonly[i]->SetLineWidth(10);
      dataresults_statonly[i]->SetFillColor(mycolors2[index[i]]);
      dataresults_statonly[i]->SetFillStyle(3001);
      dataresults_statonly_simple[i]->SetMarkerStyle(8);
      dataresults_statonly_simple[i]->SetMarkerSize(0.4);
      //      dataresults_statonly_simple[i]->SetMarkerColor(mycolors[index[i]]);

      leg4->AddEntry(mcresults[i],Form("[ %0.1f , %0.1f ), #times %0.2f",pz3Dbins[i-1],pz3Dbins[i],overallscale[i-1]),"l");
      //      leg->AddEntry(mcresults[i],Form("[ %0.1f , %0.1f )",pz3Dbins[i-1],pz3Dbins[i]),"L");
      if(doRatio){
	TF1 *myconstant = new TF1("myfunc","1",0,0.8);
	dataresults[i]->Divide(mcresults[i]);
	dataresults_statonly[i]->Divide(mcresults[i]);
	dataresults_statonly_simple[i]->Divide(mcresults[i]);
	mcresults[i]->Divide(mcresults[i]);

	double sepsize=2;
	dataresults[i]->Add(myconstant,(i-1)*sepsize);
	dataresults_statonly[i]->Add(myconstant,(i-1)*sepsize);
	dataresults_statonly_simple[i]->Add(myconstant,(i-1)*sepsize);
	mcresults[i]->Add(myconstant,(i-1)*sepsize);



      }
      //      dataresults_statonly_simple[i]->SetMarkerColor(mycolors2[index[i]]);
    // Make a list of the histograms we want to draw, along with the
    // draw options we want to use for them. You can add "graph" to the
    // draw options if you want the histogram to be converted to a graph
    // and then drawn. In that case the draw options are interpreted as
    // options to TGraphErrors::Draw().
    //
    // I don't know what happens if you put a "graph" first in the list,
    // so don't do that. Make sure the first item doesn't have "graph"
    // in its options
  

      histAndOpts.push_back(std::make_pair(dataresults[i], "histp"));
      histAndOpts.push_back(std::make_pair(mcresults[i],       "histl"));
      if(!doRatio){
	histAndOpts.push_back(std::make_pair(dataresults[i], "histle5"));
	histAndOpts.push_back(std::make_pair(dataresults_statonly[i], "histpe5"));
      }
      else{
       	histAndOpts.push_back(std::make_pair(dataresults[i], "histle5"));
	histAndOpts.push_back(std::make_pair(dataresults_statonly[i], "histpe5"));
      }
      histAndOptsDataPoints.push_back(std::make_pair(dataresults_statonly_simple[i], "histp"));


      //      histAndOpts.push_back(std::make_pair(dataresults[i],     "graph0 ep"));
  }
  cout << "Mults" << endl;
  
  histAndOpts.insert(histAndOpts.end(),histAndOptsDataPoints.begin(),histAndOptsDataPoints.end());
  
  // // ----------------------------------------------------------------------------------
  // //
  // // Now make pz in bins of pt. It's all the same
  
  // // Values to multiply each bin by to get them on a similar range
  
  // // plotXAxis1D fiddles the x axis values to squash up the tail so it
  // // doesn't take up all the horizontal space.

  TLatex *pzlabel = new TLatex(0.13,0.97, Form("%s :",zaxislabel.c_str()));
  pzlabel->SetTextSize(0.04);
  pzlabel->SetTextFont(42);
  
  TLatex *RegionA = new TLatex(0.25,0.86, "A");
  RegionA->SetTextSize(0.03);
  RegionA->SetTextFont(42);

  TArrow *LineA = new TArrow(0.15,0.85,0.37,0.85,0.01,"|-|");
  //  LineA->SetLineColor(kBlack);
  //  LineA->SetAngle(90);
  //  LineA->SetLineWidth(2);


  TLatex *RegionB = new TLatex(0.435,0.7, "B");
  RegionB->SetTextSize(0.03);
  RegionB->SetTextFont(42);
  TArrow *LineB = new TArrow(0.42,0.69,0.46,0.69,0.01,"|-|");

  TLatex *RegionB2 = new TLatex(0.135,0.3, "B");
  RegionB2->SetTextSize(0.03);
  RegionB2->SetTextFont(42);
  TArrow *LineB2 = new TArrow(0.12,0.29,0.16,0.29,0.01,"|-|");

  TLatex *RegionC1 = new TLatex(0.88,0.65, "C");
  RegionC1->SetTextSize(0.03);
  RegionC1->SetTextFont(42);
  TArrow *LineC = new TArrow(0.8,0.64,0.96,0.64,0.01,"|-|");

  TLatex *RegionD1 = new TLatex(0.405,0.3, "D");
  RegionD1->SetTextSize(0.03);
  RegionD1->SetTextFont(42);
  TArrow *LineD = new TArrow(0.4,0.29,0.42,0.29,0.01,"|-|");

  TLatex *RegionD2 = new TLatex(0.7,0.4, "D");
  RegionD2->SetTextSize(0.03);
  RegionD2->SetTextFont(42);
  TArrow *LineD2 = new TArrow(0.695,0.39,0.715,0.39,0.01,"|-|");
    

  const int nLabels=18;
  const double positions[nLabels]={0.5,1,1.5,2.5,3,3.5,4.5,5,5.5,6.5,7,7.5,8.5,9,9.5,10.5,11,11.5};
  const char* valueStrings[nLabels]={"0.5","1","1.5","0.5","1","1.5","0.5","1","1.5","0.5","1","1.5","0.5","1","1.5","0.5","1","1.5"};
  
  GridCanvas* gc= NULL;
  gc=plotXAxis1D_ReducedXRange(histAndOpts, xaxislabel, yaxislabel, 4,3,0,0.79,800,500, NULL);
  gc->SetYTitle(crosssectionlabel.c_str());
  if(!doMultipliers){
    if(!doRatio) gc->SetYLimits(0.7e-6, 12000);
    else{
      gc->SetYTitle("Ratio to Minerva Tune v1");
      gc->SetYLimits(0,13.99);
      //      gc->SetManualYLabels(nLabels, positions, valueStrings,0.03);
    }
  }
  else{
    gc->SetYLimits(0.00003, 10.49);
  }

  if(!doRatio)gc->SetLogy(true);
  //  leg->Draw("SAME");
  //  pzlabel->Draw("SAME");
  leg4->Draw("SAME");
  /*
  RegionA->Draw("SAME");
  LineA->Draw("");
  LineB->Draw("");
  LineB2->Draw("");
  LineC->Draw("");
  LineD->Draw("");
  LineD2->Draw("");
  RegionB->Draw("SAME");
  RegionB2->Draw("SAME");
  RegionC1->Draw("SAME");
  RegionD1->Draw("SAME");
  RegionD2->Draw("SAME");
  //  leg2->Draw("SAME");
  //  leg3->Draw("SAME");
  */
  gc->Modified();
  if(!doRatio){
    gc->Print(doMultipliers ? Form("nu-3d-xsec-comps-%s-pz-multiplier_bin_log_all3D.eps",variableset.c_str()) : Form("nu-3d-xsec-comps-%s-pz_bin_log_all3D.eps",variableset.c_str()));
    gc->Print(doMultipliers ? Form("nu-3d-xsec-comps-%s-pz-multiplier_bin_log_all3D.png",variableset.c_str()) : Form("nu-3d-xsec-comps-%s-pz_bin_log_all3D.png",variableset.c_str()));
    gc->Print(doMultipliers ? Form("nu-3d-xsec-comps-%s-pz-multiplier_bin_log_all3D.C",variableset.c_str()) : Form("nu-3d-xsec-comps-%s-pz_bin_log_all3D.C",variableset.c_str()));
  }
  else{
    gc->Print(doMultipliers ? Form("nu-3d-xsec-comps-%s-pz-multiplier_bin_log_all3D_ratio.eps",variableset.c_str()) : Form("nu-3d-xsec-comps-%s-pz_bin_log_all3D_ratio.eps",variableset.c_str()));
    gc->Print(doMultipliers ? Form("nu-3d-xsec-comps-%s-pz-multiplier_bin_log_all3D_ratio.png",variableset.c_str()) : Form("nu-3d-xsec-comps-%s-pz_bin_log_all3D_ratio.png",variableset.c_str()));
    gc->Print(doMultipliers ? Form("nu-3d-xsec-comps-%s-pz-multiplier_bin_log_all3D_ratio.C",variableset.c_str()) : Form("nu-3d-xsec-comps-%s-pz_bin_log_all3D_ratio.C",variableset.c_str()));

  }


  //  if(doRatio) return ;

  GridCanvas* gc2= NULL;
  vector<double> multi_1 = GetScales(histAndOpts,10.5, true);
  gc2=plotXAxis1D_ReducedXRange(histAndOpts, xaxislabel, yaxislabel, 3,2,0,0.79,800,500, &multi_1[0]);
  //gc2=plotXAxis1D(histAndOpts, xaxislabel, "P_{t,muon}", 3,2,0,1000, NULL);
  gc2->SetYLimits(0.00003, 10.49);
  gc2->SetYTitle(crosssectionlabel.c_str());
  gc2->SetLogy(true);
  leg->Draw("SAME");
  pzlabel->Draw("SAME");

  
  gc2->SetLogy(false);
  gc2->Modified();
  if(doMultipliers){//Overall scale doesn't make sense on these!
    gc2->Print(doMultipliers ? Form("nu-3d-xsec-comps-%s-pz-multiplier_bin_all3D.eps",variableset.c_str()) : Form("nu-3d-xsec-comps-%s-pz_bin_all3D.eps",variableset.c_str()));
    gc2->Print(doMultipliers ? Form("nu-3d-xsec-comps-%s-pz-multiplier_bin_all3D.png",variableset.c_str()) : Form("nu-3d-xsec-comps-%s-pz_bin_all3D.png",variableset.c_str()));
    gc2->Print(doMultipliers ? Form("nu-3d-xsec-comps-%s-pz-multiplier_bin_all3D.C",variableset.c_str()) : Form("nu-3d-xsec-comps-%s-pz_bin_all3D.C",variableset.c_str()));
  }

  GridCanvas* gc3= NULL;
  vector<double> multi_2 = GetScales(histAndOpts,10.5, true);
  gc3=plotXAxis1D_ReducedXRange(histAndOpts, xaxislabel, yaxislabel, 3,2,0,0.4,800,500, &multi_2[0]);
  //gc3=plotXAxis1D(histAndOpts, xaxislabel, "P_{t,muon}", 3,2,0,1000, NULL);
  gc3->SetYLimits(0.00003, 10.49);
  gc3->SetYTitle(crosssectionlabel.c_str());
  gc3->SetLogy(true);
  leg->Draw("SAME");
  pzlabel->Draw("SAME");

  
  gc3->SetLogy(false);
  gc3->Modified();
  if(doMultipliers){
    gc3->Print(doMultipliers ? Form("nu-3d-xsec-comps-%s-pz-multiplier_bin_zoom_low_all3D.eps",variableset.c_str()) : Form("nu-3d-xsec-comps-%s-pz_bin_zoom_low_all3D.eps",variableset.c_str()));
    gc3->Print(doMultipliers ? Form("nu-3d-xsec-comps-%s-pz-multiplier_bin_zoom_low_all3D.png",variableset.c_str()) : Form("nu-3d-xsec-comps-%s-pz_bin_zoom_low_all3D.png",variableset.c_str()));
    gc3->Print(doMultipliers ? Form("nu-3d-xsec-comps-%s-pz-multiplier_bin_zoom_low_all3D.C",variableset.c_str()) : Form("nu-3d-xsec-comps-%s-pz_bin_zoom_low_all3D.C",variableset.c_str()));
  }
  GridCanvas* gc4= NULL;
  vector<double> multi_3 = GetScales(histAndOpts,10.5, true);
  gc4=plotXAxis1D_ReducedXRange(histAndOpts, xaxislabel, yaxislabel, 3,2,0.4,0.95,800,500, &multi_3[0]);
  //gc4=plotXAxis1D(histAndOpts, xaxislabel, "P_{t,muon}", 3,2,0,1000, NULL);
  gc4->SetYLimits(0.00003, 10.49);
  gc4->SetYTitle(crosssectionlabel.c_str());
  gc4->SetLogy(true);
  leg->Draw("SAME");
  pzlabel->Draw("SAME");

  
  gc4->SetLogy(false);
  gc4->Modified();
  if(doMultipliers){ 
    gc4->Print(doMultipliers ? Form("nu-3d-xsec-comps-%s-pz-multiplier_bin_zoom_high_all3D.eps",variableset.c_str()) : Form("nu-3d-xsec-comps-%s-pz_bin_zoom_high_all3D.eps",variableset.c_str()));
    gc4->Print(doMultipliers ? Form("nu-3d-xsec-comps-%s-pz-multiplier_bin_zoom_high_all3D.png",variableset.c_str()) : Form("nu-3d-xsec-comps-%s-pz_binz_zoom_high_all3D.png",variableset.c_str()));
    gc4->Print(doMultipliers ? Form("nu-3d-xsec-comps-%s-pz-multiplier_bin_zoom_high_all3D.C",variableset.c_str()) : Form("nu-3d-xsec-comps-%s-pz_bin_zoom_high_all3D.C",variableset.c_str()));
  }

  GridCanvas* gc5= NULL;
  vector<double> multi_4 = GetScales(histAndOpts,10.5, true);
  gc5=plotXAxis1D_ReducedXRange(histAndOpts, xaxislabel, yaxislabel, 3,2,0.4,0.79,800,500, &multi_4[0]);
  //gc5=plotXAxis1D(histAndOpts, xaxislabel, "P_{t,muon}", 3,2,0,1000, NULL);
  gc5->SetYLimits(0.00003, 10.49);
  gc5->SetYTitle(crosssectionlabel.c_str());
  gc5->SetLogy(true);
  leg->Draw("SAME");
  pzlabel->Draw("SAME");

  
  gc5->SetLogy(false);
  gc5->Modified();
  if(doMultipliers){
    gc5->Print(doMultipliers ? Form("nu-3d-xsec-comps-%s-pz-multiplier_bin_zoom_allhigh_all3D.eps",variableset.c_str()) : Form("nu-3d-xsec-comps-%s-pz_bin_zoom_allhigh_all3D.eps",variableset.c_str()));
    gc5->Print(doMultipliers ? Form("nu-3d-xsec-comps-%s-pz-multiplier_bin_zoom_allhigh_all3D.png",variableset.c_str()) : Form("nu-3d-xsec-comps-%s-pz_bin_zoom_allhigh_all3D.png",variableset.c_str()));
    gc5->Print(doMultipliers ? Form("nu-3d-xsec-comps-%s-pz-multiplier_bin_zoom_allhigh_all3D.C",variableset.c_str()) : Form("nu-3d-xsec-comps-%s-pz_bin_zoom_allhigh_all3D.C",variableset.c_str()));
  }
  GridCanvas* gc6= NULL;
  vector<double> multi_5 = GetScales(histAndOpts,10.5, true);
  gc6=plotXAxis1D_ReducedXRange(histAndOpts, xaxislabel, yaxislabel, 3,2,0.,0.79,800,500, &multi_5[0]);
  //gc6=plotXAxis1D(histAndOpts, xaxislabel, "P_{t,muon}", 3,2,0,1000, NULL);
  gc6->SetYLimits(0.00003, 10.49);
  gc6->SetYTitle(crosssectionlabel.c_str());
  gc6->SetLogy(true);
  leg->Draw("SAME");
  pzlabel->Draw("SAME");

  
  gc6->SetLogy(false);
  gc6->Modified();
  if(doMultipliers){
    gc6->Print(doMultipliers ? Form("nu-3d-xsec-comps-%s-pz-multiplier_bin_zoom_allbins_all3D.eps",variableset.c_str()) : Form("nu-3d-xsec-comps-%s-pz_bin_allbins_all3D.eps",variableset.c_str()));
    gc6->Print(doMultipliers ? Form("nu-3d-xsec-comps-%s-pz-multiplier_bin_zoom_allbins_all3D.png",variableset.c_str()) : Form("nu-3d-xsec-comps-%s-pz_bin_allbins_all3D.png",variableset.c_str()));
    gc6->Print(doMultipliers ? Form("nu-3d-xsec-comps-%s-pz-multiplier_bin_zoom_allbins_all3D.C",variableset.c_str()) : Form("nu-3d-xsec-comps-%s-pz_bin_allbins_all3D.C",variableset.c_str()));
  }



  //Other direction


  GridCanvas* gc299= NULL;
  gc299=plotYAxis1D(histAndOpts, yaxislabel, xaxislabel, 4,4,800,500, NULL);
  //gc299=plotYAxis1D(histAndOpts, "P_{t,muon} (GeV)", "#Sigma T_{p}", 3,2,0,1000, NULL);
  gc299->SetYTitle(crosssectionlabel.c_str());
  if(!doMultipliers){
    if(!doRatio) gc299->SetYLimits(5e-6, 12000);
    else{
      gc299->SetYTitle("Ratio to Minerva Tune v1");
      gc299->SetYLimits(0,13.99);
      //      gc299->SetManualYLabels(nLabels, positions, valueStrings,0.03);
    }
  }
  else{
    gc299->SetYLimits(0.00003, 10.49);
  }

  if(!doRatio)gc299->SetLogy(true);
  leg->Draw("SAME");
  pzlabel->Draw("SAME");

  /*
  RegionA->Draw("SAME");
  LineA->Draw("");
  LineB->Draw("");
  LineB2->Draw("");
  LineC->Draw("");
  LineD->Draw("");
  LineD2->Draw("");
  RegionB->Draw("SAME");
  RegionB2->Draw("SAME");
  RegionC1->Draw("SAME");
  RegionD1->Draw("SAME");
  RegionD2->Draw("SAME");
  //  leg2->Draw("SAME");
  //  leg3->Draw("SAME");
  */
  gc299->Modified();
  if(!doRatio){
    gc299->Print(doMultipliers ? Form("nu-3d-xsec-comps-%s-pt-multiplier_bin_log_all3D.eps",variableset.c_str()) : Form("nu-3d-xsec-comps-%s-pt_bin_log_all3D.eps",variableset.c_str()));
    gc299->Print(doMultipliers ? Form("nu-3d-xsec-comps-%s-pt-multiplier_bin_log_all3D.png",variableset.c_str()) : Form("nu-3d-xsec-comps-%s-pt_bin_log_all3D.png",variableset.c_str()));
    gc299->Print(doMultipliers ? Form("nu-3d-xsec-comps-%s-pt-multiplier_bin_log_all3D.C",variableset.c_str()) : Form("nu-3d-xsec-comps-%s-pt_bin_log_all3D.C",variableset.c_str()));
  }
  else{
    gc299->Print(doMultipliers ? Form("nu-3d-xsec-comps-%s-pt-multiplier_bin_log_all3D_ratio.eps",variableset.c_str()) : Form("nu-3d-xsec-comps-%s-pt_bin_log_all3D_ratio.eps",variableset.c_str()));
    gc299->Print(doMultipliers ? Form("nu-3d-xsec-comps-%s-pt-multiplier_bin_log_all3D_ratio.png",variableset.c_str()) : Form("nu-3d-xsec-comps-%s-pt_bin_log_all3D_ratio.png",variableset.c_str()));
    gc299->Print(doMultipliers ? Form("nu-3d-xsec-comps-%s-pt-multiplier_bin_log_all3D_ratio.C",variableset.c_str()) : Form("nu-3d-xsec-comps-%s-pt_bin_log_all3D_ratio.C",variableset.c_str()));

  }


  //  if(doRatio) return ;

  GridCanvas* gc22= NULL;
  vector<double> multi_21 = GetScales(histAndOpts,10.5, true);
  gc22=plotYAxis1D_ReducedXRange(histAndOpts, yaxislabel, xaxislabel, 4,4,0,2.5,800,500, &multi_21[0]);
  //gc22=plotYAxis1D(histAndOpts, "P_{t,muon} (GeV)", "#Sigma T_{p}", 4,4,0,1000, NULL);
  gc22->SetYLimits(0.00003, 10.49);
  gc22->SetYTitle(crosssectionlabel.c_str());
  gc22->SetLogy(true);
  leg->Draw("SAME");
  pzlabel->Draw("SAME");

  
  gc22->SetLogy(false);
  gc22->Modified();
  if(doMultipliers){//Overall scale doesn't make sense on these!
    gc22->Print(doMultipliers ? Form("nu-3d-xsec-comps-%s-pt-multiplier_bin_all3D.eps",variableset.c_str()) : Form("nu-3d-xsec-comps-%s-pt_bin_all3D.eps",variableset.c_str()));
    gc22->Print(doMultipliers ? Form("nu-3d-xsec-comps-%s-pt-multiplier_bin_all3D.png",variableset.c_str()) : Form("nu-3d-xsec-comps-%s-pt_bin_all3D.png",variableset.c_str()));
    gc22->Print(doMultipliers ? Form("nu-3d-xsec-comps-%s-pt-multiplier_bin_all3D.C",variableset.c_str()) : Form("nu-3d-xsec-comps-%s-pt_bin_all3D.C",variableset.c_str()));
  }

  GridCanvas* gc23= NULL;
  vector<double> multi_22 = GetScales(histAndOpts,10.5, true);
  gc23=plotYAxis1D_ReducedXRange(histAndOpts, yaxislabel, xaxislabel, 4,4,0,0.4,800,500, &multi_22[0]);
  //gc23=plotYAxis1D(histAndOpts, "P_{t,muon} (GeV)", "#Sigma T_{p}", 4,4,0,1000, NULL);
  gc23->SetYLimits(0.00003, 10.49);
  gc23->SetYTitle(crosssectionlabel.c_str());
  gc23->SetLogy(true);
  leg->Draw("SAME");
  pzlabel->Draw("SAME");

  
  gc23->SetLogy(false);
  gc23->Modified();
  if(doMultipliers){
    gc23->Print(doMultipliers ? Form("nu-3d-xsec-comps-%s-pt-multiplier_bin_zoom_low_all3D.eps",variableset.c_str()) : Form("nu-3d-xsec-comps-%s-pt_bin_zoom_low_all3D.eps",variableset.c_str()));
    gc23->Print(doMultipliers ? Form("nu-3d-xsec-comps-%s-pt-multiplier_bin_zoom_low_all3D.png",variableset.c_str()) : Form("nu-3d-xsec-comps-%s-pt_bin_zoom_low_all3D.png",variableset.c_str()));
    gc23->Print(doMultipliers ? Form("nu-3d-xsec-comps-%s-pt-multiplier_bin_zoom_low_all3D.C",variableset.c_str()) : Form("nu-3d-xsec-comps-%s-pt_bin_zoom_low_all3D.C",variableset.c_str()));
  }
  GridCanvas* gc24= NULL;
  vector<double> multi_23 = GetScales(histAndOpts,10.5, true);
  gc24=plotYAxis1D_ReducedXRange(histAndOpts, yaxislabel, xaxislabel, 4,4,0.4,0.95,800,500, &multi_23[0]);
  //gc24=plotYAxis1D(histAndOpts, "P_{t,muon} (GeV)", "#Sigma T_{p}", 4,4,0,1000, NULL);
  gc24->SetYLimits(0.00003, 10.49);
  gc24->SetYTitle(crosssectionlabel.c_str());
  gc24->SetLogy(true);
  leg->Draw("SAME");
  pzlabel->Draw("SAME");

  
  gc24->SetLogy(false);
  gc24->Modified();
  if(doMultipliers){ 
    gc24->Print(doMultipliers ? Form("nu-3d-xsec-comps-%s-pt-multiplier_bin_zoom_high_all3D.eps",variableset.c_str()) : Form("nu-3d-xsec-comps-%s-pt_bin_zoom_high_all3D.eps",variableset.c_str()));
    gc24->Print(doMultipliers ? Form("nu-3d-xsec-comps-%s-pt-multiplier_bin_zoom_high_all3D.png",variableset.c_str()) : Form("nu-3d-xsec-comps-%s-pt_binz_zoom_high_all3D.png",variableset.c_str()));
    gc24->Print(doMultipliers ? Form("nu-3d-xsec-comps-%s-pt-multiplier_bin_zoom_high_all3D.C",variableset.c_str()) : Form("nu-3d-xsec-comps-%s-pt_bin_zoom_high_all3D.C",variableset.c_str()));
  }

  GridCanvas* gc25= NULL;
  vector<double> multi_24 = GetScales(histAndOpts,10.5, true);
  gc25=plotYAxis1D_ReducedXRange(histAndOpts, yaxislabel, xaxislabel, 4,4,0.4,2.5,800,500, &multi_24[0]);
  //gc25=plotYAxis1D(histAndOpts, "P_{t,muon} (GeV)", "#Sigma T_{p}", 4,4,0,1000, NULL);
  gc25->SetYLimits(0.00003, 10.49);
  gc25->SetYTitle(crosssectionlabel.c_str());
  gc25->SetLogy(true);
  leg->Draw("SAME");
  pzlabel->Draw("SAME");

  
  gc25->SetLogy(false);
  gc25->Modified();
  if(doMultipliers){
    gc25->Print(doMultipliers ? Form("nu-3d-xsec-comps-%s-pt-multiplier_bin_zoom_allhigh_all3D.eps",variableset.c_str()) : Form("nu-3d-xsec-comps-%s-pt_bin_zoom_allhigh_all3D.eps",variableset.c_str()));
    gc25->Print(doMultipliers ? Form("nu-3d-xsec-comps-%s-pt-multiplier_bin_zoom_allhigh_all3D.png",variableset.c_str()) : Form("nu-3d-xsec-comps-%s-pt_bin_zoom_allhigh_all3D.png",variableset.c_str()));
    gc25->Print(doMultipliers ? Form("nu-3d-xsec-comps-%s-pt-multiplier_bin_zoom_allhigh_all3D.C",variableset.c_str()) : Form("nu-3d-xsec-comps-%s-pt_bin_zoom_allhigh_all3D.C",variableset.c_str()));
  }
  GridCanvas* gc26= NULL;
  vector<double> multi_25 = GetScales(histAndOpts,10.5, true);
  gc26=plotYAxis1D_ReducedXRange(histAndOpts, yaxislabel, xaxislabel, 4,4,0.,2.5,800,500, &multi_25[0]);
  //gc26=plotYAxisy1D(histAndOpts, "P_{t,muon} (GeV)", "#Sigma T_{p}", 4,4,0,1000, NULL);
  gc26->SetYLimits(0.00003, 10.49);
  gc26->SetYTitle(crosssectionlabel.c_str());
  gc26->SetLogy(true);
  leg->Draw("SAME");
  pzlabel->Draw("SAME");

  
  gc26->SetLogy(false);
  gc26->Modified();
  if(doMultipliers){
    gc26->Print(doMultipliers ? Form("nu-3d-xsec-comps-%s-pt-multiplier_bin_zoom_allbins_all3D.eps",variableset.c_str()) : Form("nu-3d-xsec-comps-%s-pt_bin_allbins_all3D.eps",variableset.c_str()));
    gc26->Print(doMultipliers ? Form("nu-3d-xsec-comps-%s-pt-multiplier_bin_zoom_allbins_all3D.png",variableset.c_str()) : Form("nu-3d-xsec-comps-%s-pt_bin_allbins_all3D.png",variableset.c_str()));
    gc26->Print(doMultipliers ? Form("nu-3d-xsec-comps-%s-pt-multiplier_bin_zoom_allbins_all3D.C",variableset.c_str()) : Form("nu-3d-xsec-comps-%s-pt_bin_allbins_all3D.C",variableset.c_str()));
  }





}

int main(int argc, char* argv[])
{

  string s_pzptrec = argv[2];
  bool doPzPtRec = s_pzptrec=="1"? true:false;

  //no ratio
  makePlots(true,false,argv[1],doPzPtRec);
  makePlots(false,false,argv[1],doPzPtRec);
  
  //ratio
  makePlots(false,true,argv[1],doPzPtRec);

  return 0;
}
