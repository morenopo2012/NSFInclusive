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
#include "TMath.h"
#include "TRandom3.h"
#include "TDecompSVD.h"
//#include "localColor.h"
#include "Cintex/Cintex.h"
#include "PlotUtils/MnvColors.h"
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

void drawLabel(int xcells, int ycells, vector<string> labels, GridCanvas *c, int linenumber){
  c->cd();
  TLatex *mystring = new TLatex();
  mystring->SetTextFont(42);
  mystring->SetTextSize(0.03);
  
  //These are the top
  double xoff = (1-(c->GetLeftMargin()+c->GetRightMargin()))/xcells;
  double yoff = (1-(c->GetTopMargin()+c->GetBottomMargin()))/ycells;
  double x_shift = -0.005;
  double y_shift = -yoff*0.55;
  double additionaloff = -0.03*(1-linenumber)*3.5;

  //These are at the bottom
  double xoff_b = (1-(c->GetLeftMargin()+c->GetRightMargin()))/xcells;
  double yoff_b = (1-(c->GetBottomMargin()))/ycells;
  double x_shift_b = -0.005;
  double y_shift_b = -yoff*0.08;
  double additionaloff_b = -0.03*(1-linenumber)*3.5;



  //double additionaloff = 0;
  cout << xoff << "\t" << yoff << endl;
  int counter = 0;
  for(int j=0;j<xcells;j++){//reading order is loop over y then x
    for(int i=0;i<xcells;i++){
      if(counter>=(int)labels.size()) break;
      if(counter>=8){
	mystring->DrawLatexNDC(xoff_b*(i+1)+x_shift_b-additionaloff_b,1-(yoff_b*(j+1)+y_shift_b),labels[counter].c_str());
	counter+=1;
      }
      else{
	mystring->DrawLatexNDC(xoff*(i+1)+x_shift-additionaloff_b,1-(yoff*(j+1)+y_shift),labels[counter].c_str());
	counter+=1;
      }
    }
  }
}

void makePlots(string location, string location2="", string modelhistname="",string modelprettyname="", string location3="", string modelhistname3="",string modelprettyname3="")
{

  ROOT::Cintex::Cintex::Enable();
  myPlotStyle();
  TH1::AddDirectory(false);
  TH1::SetDefaultSumw2();
  gStyle->SetErrorX(0);
  gStyle->SetEndErrorSize(2);
  
  TFile f1(Form("%s",location.c_str()));//final

  MnvH2D* dataMnv=(MnvH2D*)f1.Get("h_pzptrec_data_nobck_unfold_effcor_cross_section");
  MnvH2D* mcMnv=(MnvH2D*)f1.Get("h_pzptrec_mc_nobck_unfold_effcor_cross_section");

  TFile *f2;
  MnvH2D* modelMnv;
  if(location2!=""){
    f2 = new TFile(location2.c_str());
    modelMnv = (MnvH2D*)f2->Get(modelhistname.c_str());
  }

  TFile *f3;
  MnvH2D* modelMnv2;
  if(location3!=""){
    f3 = new TFile(location3.c_str());
    modelMnv2 = (MnvH2D*)f3->Get(modelhistname3.c_str());
  }


  //  dataMnv->PopVertErrorBand("Reweight_Neutron");
 
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
  //for(int i=0;i<10;i++)recoil3Dbins.push_back(i*40); 
  //for(int i=0;i<9;i++)recoil3Dbins.push_back(i*40+40);//40,80,120,160,200,240,280,320,360
  for(int i=0;i<4;i++)recoil3Dbins.push_back(i*40+40);//40,80,120,160
  for(int i=0;i<2;i++)recoil3Dbins.push_back(i*80+240);//240,320
  for(int i=0;i<2;i++)recoil3Dbins.push_back(i*200+400);
  recoil3Dbins.push_back(799.0);
  for(int i=0;i<recoil3Dbins.size();i++) recoil3Dbins[i]=recoil3Dbins[i]/1000.;

  MnvH2D *Mnv_averageRecoil_Data = new MnvH2D("averageRecoil_Data","averageRecoil_Data",pz3Dbins.size()-1,&pz3Dbins[0],pt3Dbins.size()-1,&pt3Dbins[0]);
  MnvH2D *Mnv_averageRecoil_DataStat = new MnvH2D("averageRecoil_DataStat","averageRecoil_DataStat",pz3Dbins.size()-1,&pz3Dbins[0],pt3Dbins.size()-1,&pt3Dbins[0]);
  MnvH2D *Mnv_averageRecoil_MC = new MnvH2D("averageRecoil_MC","averageRecoil_MC",pz3Dbins.size()-1,&pz3Dbins[0],pt3Dbins.size()-1,&pt3Dbins[0]);
  MnvH2D *Mnv_averageRecoil_Model = new MnvH2D("averageRecoil_MC_Model","averageRecoil_MC_Model",pz3Dbins.size()-1,&pz3Dbins[0],pt3Dbins.size()-1,&pt3Dbins[0]);
  MnvH2D *Mnv_averageRecoil_Model2 = new MnvH2D("averageRecoil_MC_Model2","averageRecoil_MC_Model2",pz3Dbins.size()-1,&pz3Dbins[0],pt3Dbins.size()-1,&pt3Dbins[0]);


  std::vector<std::vector<double> > full3D;
  full3D.push_back(recoil3Dbins);
  full3D.push_back(pt3Dbins);
  full3D.push_back(pz3Dbins);
  
  HyperDimLinearizer *my3d = new HyperDimLinearizer(full3D,0);

  std::cout << "Starting up getting the projections" << std::endl;
  std::vector<MnvH2D*> dataresults = my3d->Get2DMnvHistos(dataMnv,true);
  std::vector<MnvH2D*> mcresults = my3d->Get2DMnvHistos(mcMnv,false);
  std::vector<MnvH2D*> modelresults;
  if(location2!="") modelresults = my3d->Get2DMnvHistos(modelMnv,false);

  std::vector<MnvH2D*> model2results;
  if(location3!="") model2results = my3d->Get2DMnvHistos(modelMnv2,false);


  TMatrixD mymat(dataMnv->GetSysErrorMatrix("unfoldingCov"));


  std::vector<std::pair<TH2*, const char*> > histAndOpts;
  std::vector<std::pair<TH2*, const char*> > histAndOptsRatio;
  std::vector<std::pair<TH2*, const char*> > histAndOptsDiff;
  std::vector<std::pair<TH2*, const char*> > histAndOptsBoth;

  double dataXSecTotal = 0;
  vector<double> dataPtSum;
  for(unsigned int i=0;i<pt3Dbins.size();i++) dataPtSum.push_back(0);

  double mcXSecTotal = 0;
  vector<double> mcPtSum;
  for(unsigned int i=0;i<pt3Dbins.size();i++) mcPtSum.push_back(0);

  double modelXSecTotal = 0;
  vector<double> modelPtSum;
  for(unsigned int i=0;i<pt3Dbins.size();i++) modelPtSum.push_back(0);

  double model2XSecTotal = 0;
  vector<double> model2PtSum;
  for(unsigned int i=0;i<pt3Dbins.size();i++) model2PtSum.push_back(0);
  
  for(unsigned int i=1;i<pz3Dbins.size();i++){//skip underflow and overflow pz bins
    cout << "DOING BIN " << i << endl;
    // Get the data histogram with stat error and with total error
    // separately so we can plot them both for inner and outer ticks
    // These line and marker styles will be propagated to the 1D plots
    //CV
    //dataresults[i] is a ptmu vs recoil distribution for a particular pzmubin
    for(int j=1;j<dataresults[i]->GetNbinsY()+1;j++){//y-bins, ptmu
      double data_recoil = 0;
      double data_recoil_err = 0;
      double mc_recoil = 0;
      double mc_recoil_err = 0;
      double model_recoil = 0;
      double model_recoil_err = 0;
      double model2_recoil = 0;
      double model2_recoil_err = 0;

      for(int k=1;k<dataresults[i]->GetNbinsX()+1;k++){//x-bins, recoil
	
	dataXSecTotal+=dataresults[i]->GetBinContent(k,j);
	dataPtSum[j-1]+=dataresults[i]->GetBinContent(k,j);

	mcXSecTotal+=mcresults[i]->GetBinContent(k,j);
	mcPtSum[j-1]+=mcresults[i]->GetBinContent(k,j);

	if(location2!=""){
	  modelXSecTotal+=modelresults[i]->GetBinContent(k,j);
	  modelPtSum[j-1]+=modelresults[i]->GetBinContent(k,j);
	}
	if(location3!=""){
	  model2XSecTotal+=model2results[i]->GetBinContent(k,j);
	  model2PtSum[j-1]+=model2results[i]->GetBinContent(k,j);
	}
	


	data_recoil += dataresults[i]->GetBinContent(k,j)*dataresults[i]->GetXaxis()->GetBinCenter(k);
	mc_recoil += mcresults[i]->GetBinContent(k,j)*mcresults[i]->GetXaxis()->GetBinCenter(k);
	if(location2!="")	model_recoil += modelresults[i]->GetBinContent(k,j)*modelresults[i]->GetXaxis()->GetBinCenter(k);
	if(location3!="")	model2_recoil += model2results[i]->GetBinContent(k,j)*model2results[i]->GetXaxis()->GetBinCenter(k);

	data_recoil_err += dataresults[i]->GetBinError(k,j)*dataresults[i]->GetBinError(k,j);
	mc_recoil_err += mcresults[i]->GetBinError(k,j)*mcresults[i]->GetBinError(k,j);
	if(location2!="")	model_recoil_err += modelresults[i]->GetBinError(k,j)*modelresults[i]->GetBinError(k,j);
	if(location3!="")	model2_recoil_err += model2results[i]->GetBinError(k,j)*model2results[i]->GetBinError(k,j);

	vector<double> myvals_i;
	myvals_i.push_back(recoil3Dbins[k]+0.00005);
	myvals_i.push_back(pt3Dbins[j]+0.00005);
	myvals_i.push_back(pz3Dbins[i]+0.00005);

	pair<int,int> my3Di = my3d->GetBin(myvals_i);
	int global3Di = dataMnv->GetBin(my3Di.first,j);
	//	cout << global3Di << endl;;
	for(int ll=k;ll<dataresults[i]->GetNbinsX()+1;ll++){
	  vector<double> myvals_j;
	  myvals_j.push_back(recoil3Dbins[ll]+0.00005);//recoil
	  myvals_j.push_back(pt3Dbins[j]+0.00005);//pt
	  myvals_j.push_back(pz3Dbins[i]+0.00005);//pz
	  pair<int,int> my3Dj = my3d->GetBin(myvals_j);
	  int global3Dj = dataMnv->GetBin(my3Dj.first,j);
	  //	  cout << "\t" << global3Dj << "\t" << my3Dj.first <<"\t"  << my3Dj.second << "\t" << 2*mymat[global3Di][global3Dj] << endl;;
	  data_recoil_err += 2*mymat[global3Di][global3Dj];
	}

      }

      
      
      
      double data_int = dataresults[i]->ProjectionX("tmp",j,j)->Integral(1,dataresults[i]->GetNbinsX()+1);
      double mc_int = mcresults[i]->ProjectionX("tmpmc",j,j)->Integral(1,dataresults[i]->GetNbinsX()+1);
      double model_int = modelresults[i]->ProjectionX("tmpmodel",j,j)->Integral(1,dataresults[i]->GetNbinsX()+1);
      double model2_int = model2results[i]->ProjectionX("tmpmodel2",j,j)->Integral(1,dataresults[i]->GetNbinsX()+1);
      if(i==2) cout << j << "\tData: " << data_recoil/data_int << "\tMC: " << mc_recoil/mc_int  << "\t" << TMath::Sqrt(data_recoil_err)/data_int<< endl;
      Mnv_averageRecoil_Data->SetBinContent(i,j,data_recoil/data_int);
      Mnv_averageRecoil_MC->SetBinContent(i,j,mc_recoil/mc_int);
      Mnv_averageRecoil_Model->SetBinContent(i,j,model_recoil/model_int);
      Mnv_averageRecoil_Model2->SetBinContent(i,j,model2_recoil/model2_int);
      Mnv_averageRecoil_DataStat->SetBinContent(i,j,data_recoil/data_int);
      Mnv_averageRecoil_DataStat->SetBinError(i,j,TMath::Sqrt(data_recoil_err)/data_int);
      Mnv_averageRecoil_Data->SetBinError(i,j,TMath::Sqrt(data_recoil_err)/data_int);
      Mnv_averageRecoil_MC->SetBinError(i,j,TMath::Sqrt(mc_recoil_err)/mc_int);
      Mnv_averageRecoil_Model->SetBinError(i,j,TMath::Sqrt(model_recoil_err)/model_int);
      Mnv_averageRecoil_Model2->SetBinError(i,j,TMath::Sqrt(model2_recoil_err)/model2_int);
    }    
    Mnv_averageRecoil_Data->AddMissingErrorBandsAndFillWithCV(*dataresults[i]);
    Mnv_averageRecoil_MC->AddMissingErrorBandsAndFillWithCV(*dataresults[i]);
    Mnv_averageRecoil_Model->AddMissingErrorBandsAndFillWithCV(*dataresults[i]);
    Mnv_averageRecoil_Model2->AddMissingErrorBandsAndFillWithCV(*dataresults[i]);


    
    vector<string> verts = Mnv_averageRecoil_Data->GetVertErrorBandNames();
    vector<string> lats = Mnv_averageRecoil_Data->GetLatErrorBandNames();
    
    for(unsigned int v=0;v<verts.size();v++){
      for(unsigned int u=0;u<dataresults[i]->GetVertErrorBand(verts[v])->GetNHists();u++){
	for(int j=1;j<dataresults[i]->GetNbinsY()+1;j++){
	  double data_recoil_vert = 0;
	  for(int k=1;k<dataresults[i]->GetNbinsX()+1;k++){
	    data_recoil_vert += dataresults[i]->GetVertErrorBand(verts[v])->GetHist(u)->GetBinContent(k,j)*dataresults[i]->GetXaxis()->GetBinCenter(k);
	  }
	  double data_int = dataresults[i]->GetVertErrorBand(verts[v])->GetHist(u)->ProjectionX("tmp",j,j)->Integral(1,dataresults[i]->GetNbinsX());
	  Mnv_averageRecoil_Data->GetVertErrorBand(verts[v])->GetHist(u)->SetBinContent(i,j,data_recoil_vert/data_int);
	}    
      }//loop over universes
    }//loop over vertical error bands


    for(unsigned int v=0;v<lats.size();v++){
      for(unsigned int u=0;u<dataresults[i]->GetLatErrorBand(lats[v])->GetNHists();u++){

	for(int j=1;j<dataresults[i]->GetNbinsY()+1;j++){
	  double data_recoil_lat = 0;
	  for(int k=1;k<dataresults[i]->GetNbinsX()+1;k++){
	    data_recoil_lat += dataresults[i]->GetLatErrorBand(lats[v])->GetHist(u)->GetBinContent(k,j)*dataresults[i]->GetXaxis()->GetBinCenter(k);
	  }
	  double data_int = dataresults[i]->GetLatErrorBand(lats[v])->GetHist(u)->ProjectionX("tmp",j,j)->Integral(1,dataresults[i]->GetNbinsX());
	  Mnv_averageRecoil_Data->GetLatErrorBand(lats[v])->GetHist(u)->SetBinContent(i,j,data_recoil_lat/data_int);
	}    
      }//loop over universes
    }//loop over lateral error bands
  }//loop over pz bins
  
  TH2D *averageRecoil_Data = new TH2D(Mnv_averageRecoil_Data->GetCVHistoWithError());
  TH2D *averageRecoil_Data2 = new TH2D(Mnv_averageRecoil_Data->GetCVHistoWithError());
  TH2D *averageRecoil_DataStat = new TH2D(Mnv_averageRecoil_DataStat->GetCVHistoWithError());
  TH2D *averageRecoil_MC = new TH2D(Mnv_averageRecoil_MC->GetCVHistoWithError());
  TH2D *averageRecoil_Model = new TH2D(Mnv_averageRecoil_Model->GetCVHistoWithError());
  TH2D *averageRecoil_Model2 = new TH2D(Mnv_averageRecoil_Model2->GetCVHistoWithError());
  TH2D totalFrac = Mnv_averageRecoil_Data->GetTotalError(true,true);
  TH2D statFrac = Mnv_averageRecoil_Data->GetStatError(true);
  averageRecoil_Data->SaveAs("test.root");
  totalFrac.SaveAs("FractionalError.root");
  statFrac.SaveAs("StatFracError.root");


  vector<int> mycolors = MnvColors::GetColors(9);
  averageRecoil_Data->SetMarkerStyle(8);
  averageRecoil_Data->SetMarkerSize(0.7);
  averageRecoil_Data2->SetMarkerStyle(1);
  averageRecoil_DataStat->SetMarkerStyle(1);
  averageRecoil_DataStat->SetLineWidth(10);
  averageRecoil_DataStat->SetLineColor(kMagenta-9);
  averageRecoil_MC->SetLineColor(2);
  averageRecoil_Model->SetLineColor(mycolors[4]);
  averageRecoil_Model2->SetLineColor(mycolors[2]);

  TH2D *ratioData = (TH2D*)averageRecoil_Data->Clone("RatioData");
  TH2D *ratioData2 = (TH2D*)averageRecoil_Data2->Clone("RatioData2");
  TH2D *ratioDataStat = (TH2D*)averageRecoil_DataStat->Clone("RatioDataStat");
  TH2D *ratioMC = (TH2D*)averageRecoil_MC->Clone("RatioMC");
  TH2D *ratioModel = (TH2D*)averageRecoil_Model->Clone("RatioModel");
  TH2D *ratioModel2 = (TH2D*)averageRecoil_Model2->Clone("RatioModel2");
  ratioData->Divide(averageRecoil_Data,averageRecoil_MC);
  ratioData2->Divide(averageRecoil_Data2,averageRecoil_MC);
  ratioDataStat->Divide(averageRecoil_DataStat,averageRecoil_MC);
  ratioModel->Divide(averageRecoil_Model,averageRecoil_MC);
  ratioModel2->Divide(averageRecoil_Model2,averageRecoil_MC);
  ratioMC->Divide(averageRecoil_MC,averageRecoil_MC);

  TH2D *diffData = (TH2D*)averageRecoil_Data->Clone("DiffData");
  TH2D *diffData2 = (TH2D*)averageRecoil_Data2->Clone("DiffData2");
  TH2D *diffDataStat = (TH2D*)averageRecoil_DataStat->Clone("DiffDataStat");
  TH2D *diffMC = (TH2D*)averageRecoil_MC->Clone("DiffMC");
  TH2D *diffModel = (TH2D*)averageRecoil_Model->Clone("DiffModel");
  TH2D *diffModel2 = (TH2D*)averageRecoil_Model2->Clone("DiffModel2");
  diffData->Add(averageRecoil_MC,-1);
  diffData2->Add(averageRecoil_MC,-1);
  diffDataStat->Add(averageRecoil_MC,-1);
  diffModel->Add(averageRecoil_MC,-1);
  diffModel2->Add(averageRecoil_MC,-1);
  diffMC->Add(averageRecoil_MC,-1);



  TH2D *upperMC = (TH2D*)averageRecoil_MC->Clone("MCUpper");
  TH2D *lowerMC = (TH2D*)averageRecoil_MC->Clone("MCLower");

  for(int i=1;i<upperMC->GetNbinsX()+1;i++){
    for(int j=1;j<upperMC->GetNbinsY()+1;j++){
      upperMC->SetBinContent(i,j,averageRecoil_MC->GetBinContent(i,j)+0.050);
      lowerMC->SetBinContent(i,j,averageRecoil_MC->GetBinContent(i,j)-0.050);
    }
  }
  upperMC->SetLineStyle(2);
  lowerMC->SetLineStyle(2);


  histAndOpts.push_back(std::make_pair(averageRecoil_Data, "histpe"));
  histAndOpts.push_back(std::make_pair(averageRecoil_DataStat, "histpe"));
  //  histAndOpts.push_back(std::make_pair(averageRecoil_Data2, "histpe"));
  histAndOpts.push_back(std::make_pair(averageRecoil_MC, "hist"));
  histAndOpts.push_back(std::make_pair(upperMC,"hist"));
  histAndOpts.push_back(std::make_pair(lowerMC,"hist"));
  if(location2!="")   histAndOpts.push_back(std::make_pair(averageRecoil_Model,"hist"));
  if(location3!="")   histAndOpts.push_back(std::make_pair(averageRecoil_Model2,"hist"));
  histAndOpts.push_back(std::make_pair(averageRecoil_MC, "hist"));
  histAndOpts.push_back(std::make_pair(averageRecoil_Data, "histpe1"));
  
  histAndOptsRatio.push_back(std::make_pair(ratioData, "histpe"));
  histAndOptsRatio.push_back(std::make_pair(ratioDataStat, "histpe"));
  histAndOptsRatio.push_back(std::make_pair(ratioData2, "histpe"));
  histAndOptsRatio.push_back(std::make_pair(ratioMC, "hist"));
  if(location2!="")  histAndOptsRatio.push_back(std::make_pair(ratioModel, "hist"));
  if(location3!="")  histAndOptsRatio.push_back(std::make_pair(ratioModel2, "hist"));
  histAndOptsRatio.push_back(std::make_pair(ratioData, "histpe"));

  histAndOptsDiff.push_back(std::make_pair(diffData, "histpe"));
  histAndOptsDiff.push_back(std::make_pair(diffDataStat, "histpe"));
  histAndOptsDiff.push_back(std::make_pair(diffData2, "histpe"));
  histAndOptsDiff.push_back(std::make_pair(diffMC, "hist"));
  if(location2!="")histAndOptsDiff.push_back(std::make_pair(diffModel, "hist"));
  if(location3!="")histAndOptsDiff.push_back(std::make_pair(diffModel2, "hist"));
  histAndOptsDiff.push_back(std::make_pair(diffData, "histpe"));


  TH2D *diffMCClone = (TH2D*)diffMC->Clone("mcdiffclone");
  diffMCClone->SetLineStyle(2);
  
  histAndOptsBoth.push_back(std::make_pair(averageRecoil_Data, "histpe"));
  histAndOptsBoth.push_back(std::make_pair(averageRecoil_DataStat, "histpe"));
  histAndOptsBoth.push_back(std::make_pair(averageRecoil_Data2, "histpe"));
  histAndOptsBoth.push_back(std::make_pair(averageRecoil_MC, "hist"));
  histAndOptsBoth.push_back(std::make_pair(averageRecoil_Data, "histpe"));
  histAndOptsBoth.push_back(std::make_pair(diffData, "histe"));
  histAndOptsBoth.push_back(std::make_pair(diffData2, "histe"));
  histAndOptsBoth.push_back(std::make_pair(diffMCClone, "hist"));
  histAndOptsBoth.push_back(std::make_pair(diffData, "histe"));

  TLegend *leg = new TLegend(0.79,0.1,1.0,0.4);
  leg->SetLineColor(0);

  leg->AddEntry(averageRecoil_Data,"MINERvA Data","P");
  leg->AddEntry(averageRecoil_MC, "MINERvA Tune v1","L");
  leg->AddEntry(upperMC,"#splitline{#pm 50 MeV of}{recoil energy}","L");
  if(location2!="")leg->AddEntry(averageRecoil_Model,modelprettyname.c_str(),"L");
  if(location3!="")leg->AddEntry(averageRecoil_Model2,modelprettyname3.c_str(),"L");


  vector<string> datafractions;

  for(int i=0;i<11;i++){
    string a = Form("Data: %0.2f",dataPtSum[i]/dataXSecTotal*100);
    datafractions.push_back(a+"%");
  }
  vector<string> mcfractions;
  for(int i=0;i<11;i++){
    string b = Form("MC: %0.2f",mcPtSum[i]/mcXSecTotal*100);
    mcfractions.push_back(b+"%");
  }
  
  
  GridCanvas *gc=plotXAxis1D(histAndOpts, "P_{||} (GeV/c)", "P_{t} (GeV/c)", NULL);
  gc->SetYLimits(0,0.59);
  drawLabel(4,3,datafractions,gc,1);
  drawLabel(4,3,mcfractions,gc,2);
  gc->SetYTitle("Average #Sigma T_{p} (GeV)");
  leg->Draw("SAME");
  gc->Print("AverageRecoil.png");
  gc->Print("AverageRecoil.eps");
  gc->Print("AverageRecoil.C");

  GridCanvas *gc2=plotXAxis1D(histAndOptsRatio, "P_{||} (GeV/c)", "P_{t} (GeV/c)", NULL);
  gc2->SetYLimits(0,1.99);
  drawLabel(4,3,datafractions,gc2,1);
  drawLabel(4,3,mcfractions,gc2,2);
  gc2->SetYTitle("Ratio Data To Minerva Tune v1");
  leg->Draw("SAME");
  gc2->Print("AverageRecoilRatio.png");
  gc2->Print("AverageRecoilRatio.eps");
  gc2->Print("AverageRecoilRatio.C");


  GridCanvas *gc3=plotXAxis1D(histAndOptsDiff, "P_{||} (GeV/c)", "P_{t} (GeV/c)", NULL);
  gc3->SetYLimits(-0.24,0.24);
  drawLabel(4,3,datafractions,gc3,1);
  drawLabel(4,3,mcfractions,gc3,2);
  gc3->SetYTitle("Difference Data - Minerva Tune v1");
  leg->Draw("SAME");
  gc3->Print("AverageRecoilDiff.png");
  gc3->Print("AverageRecoilDiff.eps");
  gc3->Print("AverageRecoilDiff.C");


  GridCanvas *gc4=plotXAxis1D(histAndOptsBoth, "P_{||} (GeV/c)", "P_{t} (GeV/c)", NULL);
  gc4->SetYLimits(-0.11,0.59);
  drawLabel(4,3,datafractions,gc4,1);
  drawLabel(4,3,mcfractions,gc4,2);
  gc4->SetYTitle("Average #Sigma T_{p} (GeV) And Difference (GeV)");
  leg->Draw("SAME");
  gc4->Print("AverageRecoilBoth.png");
  gc4->Print("AverageRecoilBoth.eps");
  gc4->Print("AverageRecoilBoth.C");





  
}

int main(int argc, char* argv[])
{
  //no ratio
  cout << argc << endl;
  for(int i=0;i<argc;i++) cout << i<<"\t"<<argv[i] << endl;
  if(argc==2)makePlots(argv[1],"","","","","","");
  if(argc==6)makePlots(argv[1],argv[2],argv[3],"","","");
  if(argc==8)makePlots(argv[1],argv[2],argv[3],argv[4],argv[5],argv[6],argv[7]);


  return 0;
}
