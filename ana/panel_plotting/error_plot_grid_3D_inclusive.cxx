//  #include "../util/plot/myPlotStyle.h"
#include "GridCanvas.h"
#include "PlotUtils/MnvH2D.h"
#include "PlotUtils/MnvPlotter.h"
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

#include "Cintex/Cintex.h"

#include "myPlotStyle.h"

#include "plot.h"

#include <iostream>

using namespace std;
using namespace PlotUtils;

// =====================================================================
// axis==1: x axis. axis==2: y axis

// Take a vector of 1D histograms which all have the same binning and
// stack them together into a 2D histogram. The axis argument says
// which axis to stack them on
TH2* concatenateHists(vector<TH1*>& hists1D, int axis)
{
  assert(hists1D.size());
  
  //These are eavail vs pt vs pz
  
   int nyBins_pzptrec=11;
   double yBins_pzptrec[nyBins_pzptrec+1]={0,0.075,0.15,0.25,0.325,0.4,0.475,0.55,0.7,0.85,1.0,2.5};
  
  
   int nxBins_pzptrec=14;
   double xBins_pzptrec[nxBins_pzptrec+1]={0,0.06,0.17,0.34,0.6,0.98,1.56,2.34,3.3,4.49,7.79,12.63,19.55,29.63,50.0};


  //These are eavail vs q0qe vs emu
   int nxBins_enuproxy=10;
   double xBins_enuproxy[nxBins_enuproxy+1]={0,0.020,0.040,0.080,0.120,0.160,0.240,0.320,0.400,0.600,0.7999};

   int nyBins_enuproxy=10;
   double yBins_enuproxy[nyBins_enuproxy+1]={0,0.020,0.040,0.080,0.120,0.160,0.240,0.320,0.400,0.600,0.7999};
  
  


   string histname = hists1D[0]->GetName();

   int nxBins=histname.find("pzptrec")==string::npos? nxBins_enuproxy:nxBins_pzptrec;
   double xBins[nxBins+1];
   if(histname.find("pzptrec")==string::npos){
     for(int i=0;i<nxBins+1;i++)xBins[i]=xBins_enuproxy[i];
   }
   else{
     for(int i=0;i<nxBins+1;i++)xBins[i]=xBins_pzptrec[i];
   }

   int nyBins=histname.find("pzptrec")==string::npos? nyBins_enuproxy:nyBins_pzptrec;
   double yBins[nyBins+1];
   if(histname.find("pzptrec")==string::npos){
     for(int i=0;i<nyBins+1;i++) yBins[i]=yBins_enuproxy[i];
   }
   else{
     for(int i=0;i<nyBins+1;i++) yBins[i]=yBins_pzptrec[i];
   }
   

  


  TH2* ret=0;
  if(axis==1){
    ret=new TH2D(uniq(), TString::Format(";%s", hists1D[0]->GetXaxis()->GetTitle()),
                 hists1D[0]->GetXaxis()->GetNbins(), hists1D[0]->GetXaxis()->GetXbins()->GetArray(),
                 nyBins, yBins);
  }
  else{
    ret=new TH2D(uniq(), TString::Format(";;%s", hists1D[0]->GetXaxis()->GetTitle()),
                 nxBins, xBins,
                 hists1D[0]->GetXaxis()->GetNbins(), hists1D[0]->GetXaxis()->GetXbins()->GetArray());
  }

  ret->SetLineColor(hists1D[0]->GetLineColor());
  ret->SetLineStyle(hists1D[0]->GetLineStyle());

  for(unsigned int iHist=0; iHist<hists1D.size(); ++iHist){
    for(int j=0; j<hists1D[0]->GetXaxis()->GetNbins()+1; ++j){
      int ixBin=axis==1 ? j       : iHist+1;
      int iyBin=axis==1 ? iHist+1 : j;
      double content=hists1D[iHist]->GetBinContent(j);
      ret->SetBinContent(ixBin, iyBin, content);
    }
  }

  return ret;
}

// =====================================================================
vector<std::pair<TH2*, const char*> > getSystHistsAndOpts(MnvH2D* data, bool pt,TLegend *&leg, string group = "",bool isFrac=true)
{
  MnvPlotter plotter;
  plotter.ApplyStyle(kCCQENuStyle);

  vector<string> vertnames = data->GetVertErrorBandNames();
  vector<string> latnames = data->GetLatErrorBandNames();
  // For each bin in the other variable, make a vector of the
  // systematic histograms
  const int nBins=pt ? 10 : 11;
  vector<vector<TH1*> > histsPT;
  histsPT.resize(nBins);

  // Get MnvPlotter to plot all the histograms, and slurp them into histsPT
  //  data->SaveAs(Form("hist_%s.root",data->GetName()));
    
  for(int i=0; i<nBins; ++i){
    // First plot the histograms in the dummy canvas...
    TCanvas c;
    MnvH1D* proj=pt ? data->ProjectionY(uniq(), i+1, i+1) : data->ProjectionX(uniq(), i+1, i+1);
    if(group==""){
      plotter.DrawErrorSummary(proj, "TR", true, true, -1, false,"Others", isFrac);
      //plotter.DrawErrorSummary(proj, "TR", true, true, -1, false,"Cross Section Models", true);
      if(i==0){
	leg = new TLegend(*getPadLegend(&c));
	leg->SetName("MyLegend");
	
      }
    }
    else{
      /*
      if(plotter.error_summary_group_map.find(group)!=plotter.error_summary_group_map.end()){
	plotter.DrawErrorSummary(proj, "TR", true, true, -1, false,group, isFrac);
	//plotter.DrawErrorSummary(proj, "TR", true, true, -1, false,"Cross Section Models", true);
	TCanvas c2;
	c2.cd();
	leg = new TLegend(*getPadLegend(&c));
	leg->SetName("MyLegend");
	leg->SetNColumns(3);
	leg->SetX1(0);
	leg->SetX2(1);
	leg->SetY1(0);
	leg->SetY2(1);
	leg->Draw();
	c2.Print(Form("Legend_%s_%d.png",group.c_str(),i));
	c2.Print(Form("Legend_%s_%d.C",group.c_str(),i));
      }
      else{
      */
	if(std::count(vertnames.begin(),vertnames.end(),group)){
	  TH1D err = proj->GetVertErrorBand(group)->GetErrorBand(true);
	  err.DrawClone();
	}
	else{
	  TH1D err = proj->GetLatErrorBand(group)->GetErrorBand(true);
	  err.DrawClone();
	}
	
	//      }
    }
    c.ls();
    std::vector<TH1*> padHists=getPadHists(&c);
    histsPT[i]=padHists;
  }

  // concatenateHists wants a vector of hists for each of the bins of
  // a given systematic. But histsPT is the other way round (the inner
  // vector loops over systematics).  So we have this fiddly loop to
  // make a transposed version of the vector-of-vector

  // It would have been easier to just pass the original
  // vector-of-vector into concatenateHists, and tell it which
  // systematic we wanted, but I've written and debugged this now, so
  // not changing it

  //  First index is systematic, second
  // index is bin
  vector<vector<TH1*> > histsPT_transpose;
  int nSyst=histsPT[0].size();
  cout << "There are " << nSyst << " systematics" << endl;
  histsPT_transpose.resize(nSyst);

  for(int iSyst=0; iSyst<nSyst; ++iSyst){
    for(unsigned int iBin=0; iBin<histsPT.size(); ++iBin){  
      histsPT_transpose[iSyst].push_back(histsPT[iBin][iSyst]);
    }
  }

  vector<std::pair<TH2*, const char*> > histsPT2D;
  // TODO: Figure out why the last systematic is crashing
  for(int iSyst=0; iSyst<histsPT_transpose.size(); ++iSyst){
    TH2* h2d=concatenateHists(histsPT_transpose[iSyst], pt ? 2 : 1);
    // We want to draw all of these histograms as graphs, and exclude
    // the zero bins, to get rid of ROOT artifacts. The "graph0" draw
    // option does that (and I made it safe to pass all graphs)
    histsPT2D.push_back(std::make_pair(h2d, "graph0 l"));
  }

  return histsPT2D;
}

// =====================================================================
int makePlots(bool pt, string location, bool drawGroups, string histoname, bool isFrac)
{
  // This turns out to be complicated (well, I could have made it less
  // bad, but this way is complicated and generalizable):
  //
  // MnvPlotter knows how to make the histograms we need, including
  // fancy stuff like grouping systematics together and sticking to
  // colour schemes. But it doesn't know about GridCanvas, and it
  // goes off and draws the histograms in its own way
  //
  // So we're going to take projections, and ask MnvPlotter to plot
  // the projection in a dummy canvas. Then we'll grab all the
  // histograms from that canvas and hold onto them.
  //
  // We could just grab all those 1D histograms and plot them straight
  // into the appropriate panel of our GridCanvas, but I want to reuse
  // the plotpz1D and plotpT1D functions in plot.h, because they know
  // how to do fanciness like squashing the tail in pz. Those
  // functions take a vector of 2D histograms, so we have to take our
  // 1D histograms from MnvPlotter, and stack them back together into
  // 2D histograms. That's what the concatenateHists function does

  ROOT::Cintex::Cintex::Enable();
  myPlotStyle();
  TH1::AddDirectory(false);
  TH1::SetDefaultSumw2();
  gStyle->SetErrorX(0);
  gStyle->SetEndErrorSize(2);

  TFile f(Form("%s",location.c_str()));
  MnvH2D* dataMnv=(MnvH2D*)f.Get(histoname.c_str());
  //  MnvH2D* dataMnv=(MnvH2D*)f.Get("h_pzptrec_data_nobck_unfold_effcor_cross_section");
  //  MnvH2D* dataMnv=(MnvH2D*)f.Get("h_pzptrec_mc_nobck_unfold_effcor");
  //MnvH2D* dataMnv=(MnvH2D*)f.Get("h_pzptrec_mc_nobck_unfold");

  //  dataMnv->PopSysErrorMatrix("unfoldingCov");

  if(!dataMnv){
    cout << "Failed to get the histogram" << endl;
    return 1;
  }
  cout << dataMnv << endl;
  //  dataMnv->PopVertErrorBand("Reweight_Neutron");
  dataMnv->PopVertErrorBand("GENIE_Theta_Delta2Npi");
  dataMnv->PopVertErrorBand("GENIE_FrElas_N");
  dataMnv->PopVertErrorBand("GENIE_NormCCQE");
  dataMnv->PopVertErrorBand("GENIE_MaCCQEshape");
  dataMnv->PopVertErrorBand("HybridCut_Efficiency");

  vector<double> recoil3Dbins;
  vector<double> pt3Dbins;
  vector<double> pz3Dbins;

  //Nominal W&C
  /*
  pt3Dbins.push_back(0.0);
  pt3Dbins.push_back(0.15);
  pt3Dbins.push_back(0.25);//added ME
  pt3Dbins.push_back(0.4);
  pt3Dbins.push_back(0.7);//added ME
  pt3Dbins.push_back(1.0);
  pt3Dbins.push_back(2.5);
  */

  
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


  recoil3Dbins.push_back(0);
  recoil3Dbins.push_back(0.06);
  recoil3Dbins.push_back(0.17);
  recoil3Dbins.push_back(0.34);
  recoil3Dbins.push_back(0.6);
  recoil3Dbins.push_back(0.98);
  recoil3Dbins.push_back(1.56);
  recoil3Dbins.push_back(2.34);
  recoil3Dbins.push_back(3.3);
  recoil3Dbins.push_back(4.49);
  //  recoil3Dbins.push_back(5.95);
  recoil3Dbins.push_back(7.79);
  //  recoil3Dbins.push_back(10.01);
  recoil3Dbins.push_back(12.63);
  //  recoil3Dbins.push_back(15.75);
  recoil3Dbins.push_back(19.55);
  //  recoil3Dbins.push_back(24.11);
  recoil3Dbins.push_back(29.63);
  //  recoil3Dbins.push_back(36.22);
  recoil3Dbins.push_back(50.0);

  for(int i=0;i<recoil3Dbins.size();i++) cout << recoil3Dbins[i] << endl;


  std::vector<std::vector<double> > full3D_pzptrec;
  full3D_pzptrec.push_back(recoil3Dbins);
  full3D_pzptrec.push_back(pt3Dbins);
  full3D_pzptrec.push_back(pz3Dbins);

  std::vector<std::vector<double> > full3D_enuproxy;
  full3D_enuproxy.push_back(recoil3Dbins);
  full3D_enuproxy.push_back(recoil3Dbins);
  full3D_enuproxy.push_back(pz3Dbins);
  
  std::vector<std::vector<double> > full3D;
  if(histoname.find("enuproxyE")!=string::npos) full3D=full3D_enuproxy;
  else full3D = full3D_pzptrec;

  HyperDimLinearizer *my3d = new HyperDimLinearizer(full3D,0);
  PlotUtils::MnvPlotter *plotter = new PlotUtils::MnvPlotter();
  plotter->ApplyStyle(kCCQENuInclusiveStyle);
  vector<string> names;
  for(std::map<string,vector<string> >::iterator it=plotter->error_summary_group_map.begin(); it!=plotter->error_summary_group_map.end();++it) names.push_back(it->first);

  string labelx = "#Sigma T_{P} (GeV)";
  string labely = "P_{t} (GeV/c)";
  string labelz = "P_{||} (GeV/c)";

  if(histoname.find("enuproxyE")!=string::npos){
    labely = "q^{QE}_{0} (GeV)";
    labelz = "E_{#mu} (GeV)";
  }

  std::cout << "Starting up getting the projections" << std::endl;
  std::vector<MnvH2D*> dataresults = my3d->Get2DMnvHistos(dataMnv,true);
  TLegend *leg = NULL;  
  for(int i=1;i<dataresults.size()-1;i++){
    MnvH2D *data = dataresults[i];
    vector<string> names2 = data->GetErrorBandNames();

    TLatex mytex;
    mytex.SetTextSize(0.05);
    string mystring =     Form("%.2f < %s < %.2f",pz3Dbins[i-1],labelz.c_str(),pz3Dbins[i]);
    if(!drawGroups){

      vector<std::pair<TH2*, const char*> > histsPT2D=getSystHistsAndOpts(data,  pt, *&leg, "", isFrac);   

      leg->SetLineColor(kWhite);
      leg->SetFillColor(kWhite);

      if(pt){
	vector<double> mults = GetScales(histsPT2D,false,1e-41,0.9,false);
	for(int i=0;i<mults.size();i++){
	  if(isFrac){
	    mults[i]=1;
	  }
	  else{
	    mults[i]<0?mults[i]=1.0:mults[i]=mults[i];
	  }
	}
	GridCanvas* gcPT=plotYAxis1D(histsPT2D,labely,labelx,4,3,800,500,&mults[0]);
	if(isFrac) gcPT->SetYTitle("Fractional uncertainty");
	else  gcPT->SetYTitle("Uncertainty");
	if(isFrac)gcPT->SetYLimits(0, 0.49);
	else gcPT->SetYLimits(0, 1e-41);
	//gcPT->SetYLimits(0, 1e12);
	//	gcPT->SetLogy(true);
	leg->SetX1(0.58);
	leg->SetY1(0.05);
	leg->SetX2(0.95);
	leg->SetY2(0.32);
	leg->Draw("SAME");
	mytex.DrawLatex(0.35,0.96,mystring.c_str());
	gcPT->Modified();
	if(isFrac){
	  gcPT->Print(Form("errors-pt_%d_fractional.eps",i));
	  gcPT->Print(Form("errors-pt_%d_fractional.png",i));
	  gcPT->Print(Form("errors-pt_%d_fractional.C",i));
	}
	else{
	  gcPT->Print(Form("errors-pt_%d.eps",i));
	  gcPT->Print(Form("errors-pt_%d.png",i));
	  gcPT->Print(Form("errors-pt_%d.C",i));
	}
      }
      else{
	vector<double> mults = GetScales(histsPT2D,true,1e-41,0.9,false);
	for(int i=0;i<mults.size();i++){
	  if(isFrac){
	    mults[i]=1;
	  }
	  else{
	    mults[i]<0?mults[i]=1.0:mults[i]=mults[i];
	  }
	}
	GridCanvas* gcPT=plotXAxis1D(histsPT2D,labelx,labely,4,3,800,500,&mults[0]);
	if(isFrac) gcPT->SetYTitle("Fractional uncertainty");
	else  gcPT->SetYTitle("Uncertainty");
	if(isFrac) gcPT->SetYLimits(0, 0.49);
	else gcPT->SetYLimits(0, 1e-41);
	//gcPT->SetYLimits(0, 1e3);
	//	gcPT->SetLogy(true);
	mytex.DrawLatex(0.35,0.96,mystring.c_str());
	leg->SetX1(0.77);
	leg->SetY1(0.1);
	leg->SetX2(1.0);
	leg->SetY2(0.4);
	leg->Draw("SAME");

	gcPT->Modified();
	if(isFrac){
	  gcPT->Print(Form("errors-pz_%d_fractional.eps",i));
	  gcPT->Print(Form("errors-pz_%d_fractional.png",i));
	  gcPT->Print(Form("errors-pz_%d_fractional.C",i));
	}
	else{
	  gcPT->Print(Form("errors-pz_%d.eps",i));
	  gcPT->Print(Form("errors-pz_%d.png",i));
	  gcPT->Print(Form("errors-pz_%d.C",i));
	}
      }
    }
    else{
      /*
      for(int n=0;n<names.size();n++){
	string group = names[n];
	
	vector<std::pair<TH2*, const char*> > histsPT2D=getSystHistsAndOpts(data,  pt, *&leg, group);   
	TLegend *leg = new TLegend(0.17, 0.7, 0.31, 0.9);
	leg->SetLineColor(kWhite);
	leg->SetFillColor(kWhite);
      	for(int j=0;j<histsPT2D.size();j++)leg->AddEntry(histsPT2D[j].first,group.c_str(),"l");
	if(pt){
	  GridCanvas* gcPT=plotYAxis1D(histsPT2D,"P_{t,#mu} (GeV/c)","#Sigma T_{p}",4,4,800,500);
	  gcPT->SetYTitle("Fractional uncertainty");
	  gcPT->SetYLimits(0, 0.249);
	  leg->Draw("SAME");
	  gcPT->Modified();
	  gcPT->Print(Form("errors-%s-pt_%d.eps",group.c_str(),i));
	  gcPT->Print(Form("errors-%s-pt_%d.png",group.c_str(),i));
	  gcPT->Print(Form("errors-%s-pt_%d.C",group.c_str(),i));
	}
	else{
	  GridCanvas* gcPT=plotXAxis1D(histsPT2D,"#Sigma T_{p} (GeV)","P_{t,#mu}",4,3,800,500);
	  gcPT->SetYTitle("Fractional uncertainty");
	  gcPT->SetYLimits(0, 0.249);
	  gcPT->Modified();
	  gcPT->Print(Form("errors-%s-pz_%d.eps",group.c_str(),i));
	  gcPT->Print(Form("errors-%s-pz_%d.png",group.c_str(),i));
	  gcPT->Print(Form("errors-%s-pz_%d.C",group.c_str(),i));
	}
      }
      */
      for(int n=0;n<names2.size();n++){
	string group = names2[n];
	
	vector<std::pair<TH2*, const char*> > histsPT2D=getSystHistsAndOpts(data,  pt, *&leg, group);   
	TLegend *leg = new TLegend(0.17, 0.7, 0.31, 0.9);
	leg->SetLineColor(kWhite);
	leg->SetFillColor(kWhite);
      	for(int j=0;j<histsPT2D.size();j++)leg->AddEntry(histsPT2D[j].first,group.c_str(),"l");
	if(pt){
	  GridCanvas* gcPT=plotYAxis1D(histsPT2D,"P_{t,#mu} (GeV/c)","#Sigma T_{p}",4,4,800,500);
	  gcPT->SetYTitle("Fractional uncertainty");
	  gcPT->SetYLimits(0, 0.249);
	  leg->Draw("SAME");
	  gcPT->Modified();
	  gcPT->Print(Form("errors-%s-pt_%d.eps",group.c_str(),i));
	  gcPT->Print(Form("errors-%s-pt_%d.png",group.c_str(),i));
	  gcPT->Print(Form("errors-%s-pt_%d.C",group.c_str(),i));
	}
	else{
	  GridCanvas* gcPT=plotXAxis1D(histsPT2D,"#Sigma T_{p} (GeV)","P_{t,#mu}",4,3,800,500);
	  gcPT->SetYTitle("Fractional uncertainty");
	  gcPT->SetYLimits(0, 0.249);
	  gcPT->Modified();
	  gcPT->Print(Form("errors-%s-pz_%d.eps",group.c_str(),i));
	  gcPT->Print(Form("errors-%s-pz_%d.png",group.c_str(),i));
	  gcPT->Print(Form("errors-%s-pz_%d.C",group.c_str(),i));
	}
      }
    }
  }
}
int main(int argc, char* argv[])
{
  //  makePlots(true,argv[1],false, argv[2], false);
  makePlots(false,argv[1],false, argv[2], false);

    //makePlots(true,argv[1],false, argv[2], true);
  makePlots(false,argv[1],false, argv[2], true);
  //Draw individual sources
  //makePlots(true,argv[1],true, argv[2],true);
   makePlots(false,argv[1],true, argv[2],true);
  
  return 0;
}
