#include <bits/stdc++.h>
#include "TClass.h"
#include "TList.h"
#include "TH1.h"
#include "TH2.h"
#include "TPad.h"
#include "TKey.h"


#include <iostream>
#include <vector>
#include <math.h>
#include <iostream>

#include "GridCanvas.h"
#include "PlotUtils/MnvH2D.h"
#include "PlotUtils/MnvPlotter.h"

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
#include "localColor.h"
#include "Cintex/Cintex.h"

#include "myPlotStyle.h"

// Deepika 
// Changed to fit the non-grid plots
// plot.h dependency is removed. 

using namespace std;
#include "vartex.h"
std::map<string, string> varname;

using namespace PlotUtils;

TString uniq() {
  static int i=0;
  return TString::Format("uniq%d", i++);
}

std::vector<TH1*> getPadHists(TPad* pad)
{
  std::vector<TH1*> ret;
  std::vector<string> names;
  TIter next(pad->GetListOfPrimitives());
  TObject *obj;
  while (( obj=next() )) {
    //    cout << obj->GetName() << endl;
    if ( obj->IsA()->InheritsFrom(TH1::Class()) ) {
      string name = obj->GetName();
      //      cout << "getPadAxisHist returning hist with name " << obj->GetName() <<"t"<< obj->GetTitle() << endl;
      bool isFound = false;
      for(unsigned int i=0;i<names.size();i++){
        if(names[i]==name){
          isFound = true;
          break;
        }
      }
      if(!isFound){
        ret.push_back((TH1*)obj);
        names.push_back(name);
      }
      
    }
  }
  return ret;
}




vector<vector<double>> GetKinValue(int binsx,int binsy,TH1 *rebinnedhisto,TH2 *templatehist) {
  vector<vector <double>> tmpv;
  for (int i = 1; i < binsy-1; i++) {
    vector<double>binsptbin;
    for (int j = 0; j < binsx-1; j++) {
      int binid      = rebinnedhisto->FindBin(i*binsx+j+1);
      int mapped_low = rebinnedhisto->GetBinLowEdge(binid);
      int xindex     = mapped_low%binsx;
      double kin_val = templatehist->GetXaxis()->GetBinLowEdge(int(xindex));
      binsptbin.push_back(kin_val);    
      cout << i << " " << j << "  "<< binid << " " << mapped_low << " " << xindex  << " " << kin_val << endl; 
    } 
    //it = unique(binsptbin.begin(), binsptbin.end());
    //binsptbin.resize(distance(binsptbin.begin(), it));
    
    std::sort(binsptbin.begin(), binsptbin.end()); 
    auto last = std::unique(binsptbin.begin(), binsptbin.end());
    binsptbin.erase(last, binsptbin.end());
    
    cout << " Bin " << i << " : ";
    for(auto it:binsptbin) {
      cout <<it << " ";
    }
    cout << endl;
    
    tmpv.push_back(binsptbin);
  }
  return tmpv;
}

MnvH1D* makeTH1D(int pt_bin,double pt_width, vector<double> boundaries,MnvH1D* rebinnedhisto,int startindex,const char *prefix) {
  
  int n = boundaries.size();
    
  double arr[n];
  cout << " Ptbin " << pt_bin << " " << startindex << " ";
  
  for (unsigned int i = 0; i < boundaries.size(); i++) {
    arr[i] = boundaries[i];
    cout  << arr[i] << " ";
    
  }
  
  MnvH1D* retTH1D = new MnvH1D(Form("%s_pt_bin_%d",prefix,pt_bin),Form("%s_pt_bin_%d", prefix,pt_bin),n-1,arr);
  cout << " ------ Lat ------" << endl;
  vector <string> var_leb = rebinnedhisto->GetLatErrorBandNames();
  for(unsigned int i = 0; i < var_leb.size(); i++) {
    cout << var_leb[i].c_str() <<endl;
    retTH1D->AddLatErrorBand(var_leb[i].c_str());
  }
  cout << " ------ Vet ------" << endl;
  vector <string> var_veb = rebinnedhisto->GetVertErrorBandNames();
  for(unsigned int i = 0; i < var_veb.size(); i++) {
    cout << var_veb[i].c_str() <<endl;
    retTH1D->AddVertErrorBand(var_veb[i].c_str());
  }
  
  /*vector <string> var_eb = dataMnv->GetErrorBandNames();
    for(unsigned int i = 0; i < var_eb.size(); i++) {
    cout << var_eb[i].c_str() <<endl;
    
    }*/
  
  
  for (unsigned int i = 0; i < boundaries.size(); i++) {
    
    double bv      = rebinnedhisto->GetBinContent(startindex+i)/(pt_width);
    double bv_err  = rebinnedhisto->GetBinError(startindex+i)/(pt_width);
    retTH1D->SetBinContent(i+1,bv);
    retTH1D->SetBinError(i+1,bv_err);
    
    for(unsigned int ii = 0; ii < var_veb.size(); ii++) {
      int nhists =  rebinnedhisto->GetVertErrorBand(var_veb[ii].c_str())->GetNHists();
      cout << "=========> " << i << "  " << ii << "  " << nhists << endl;
      for (int nh = 0; nh < nhists; nh++) {
	double ssa = rebinnedhisto->GetVertErrorBand(var_veb[ii].c_str())->GetHists().at(nh)->GetBinContent(startindex+i);
	cout <<  startindex+i << " " << i << " " << ii << " " << nh << " " << ssa << endl;
	retTH1D->GetVertErrorBand(var_veb[ii].c_str())->GetHists().at(nh)->SetBinContent(i+1, (ssa/pt_width));
      }
    }
    
    for(unsigned int ii = 0; ii < var_leb.size(); ii++) {
      int nhists =  rebinnedhisto->GetLatErrorBand(var_leb[ii].c_str())->GetNHists();
      cout << "=========> " << i << "  " << ii << "  " << nhists << endl;
      for (int nh = 0; nh < nhists; nh++) {
         double ssa = rebinnedhisto->GetLatErrorBand(var_leb[ii].c_str())->GetHists().at(nh)->GetBinContent(startindex+i);
	 cout <<  startindex+i << " " << i << " " << ii << " " << nh << " " << ssa << endl;  
	 retTH1D->GetLatErrorBand(var_leb[ii].c_str())->GetHists().at(nh)->SetBinContent(i+1, (ssa/pt_width));
      }
    }
  }
  cout << " " <<retTH1D->GetNbinsX() << endl;

  return retTH1D;
}

//_________________________________
vector <MnvH1D*> MakeTH1Ds(vector<vector<double>> kinVs, MnvH1D *rebinnedhisto,TH2 *templatehist, const char *prefix) {
  vector <MnvH1D*> myTH1Ds;
  int counter = templatehist->GetNbinsX()+4;
  int i = 0;
  cout <<"--------------------------------------------"<<endl;
  for (auto& it : kinVs) {
    
    MnvH1D* tmp_th1D = makeTH1D(i,templatehist->GetYaxis()->GetBinWidth(i+1),it,rebinnedhisto,counter,prefix);
    counter+=it.size()+1;
    myTH1Ds.push_back(tmp_th1D);

    cout << " Bin " << i <<" : ";
    for (auto& itt : it) {
      cout << itt << " ";
    }
    cout <<endl;
    i++;
  }
  
//vector <MnvH1D*> myTH1DsT;


  return myTH1Ds;
}

vector<double> GetScales(vector <vector <TH1*>> histopts){
  vector<double> tmpvect;
  cout  << " <<<<<<< Calculating the Scale Factor >>>>>>>>>> " << endl;
  int nh = histopts.size();
  int nb = histopts[0].size();
  
  double mxv[nb][nh];
  
  cout << nh  << "  " << nb << endl;
  int i = 0;
  for (auto& it : histopts) {
    int bin = 0;
    for (auto& itt : it) {
      int maxbin = itt->GetMaximumBin();
      mxv[bin][i] = itt->GetBinContent(maxbin);
      cout << itt->GetName() << " " << i << "  " << bin << " " << mxv[bin][i] << " " << maxbin  << endl;
      bin++;
    }
    i++;
  }
  

  for(int ii = 0; ii < nb; ii++){
    double maxval = 0;
    for(int jj = 0; jj < nh; jj++){
      cout  << ii << " " << " " << jj << " " << mxv[ii][jj] << endl;
      double content =  mxv[ii][jj];
      if (content > maxval) maxval = content;
    }
    cout << maxval << endl;
    double scale = 3.6/maxval;
    if(scale>1){
       int tmpscale = floor(scale*10);
       scale = (tmpscale/10.0)-0.5;
       if (scale > 10) scale = scale - 4;
    }
    else{
      int tmpscale = ceil(scale*1000);
      scale = tmpscale/1000.0;
      }
    //scale = scale;
    cout << "Bin " << ii << " Scale Factor " << scale << endl;
    tmpvect.push_back(scale);
  }


  cout << " <<<<<<<<<<< ------------- end scale value calculation  --------- >>>>>>>>> " << endl;
  return tmpvect;
}

//___________________________________________________________________________
vector <MnvH1D*> PrepareHist(TH2* hisy_temp, MnvH1D *dataMnv, const char *prefix){
  int n_binsx_original = hisy_temp->GetNbinsX()+2;
  int n_binsy_original = hisy_temp->GetNbinsY()+2;
  cout << "==================== " << n_binsx_original << "  " << n_binsy_original <<endl;
  TH1* data           = new TH1D(dataMnv->GetCVHistoWithError());
  vector<vector<double>> kin_value = GetKinValue(n_binsx_original, n_binsy_original, data, hisy_temp);
  vector <MnvH1D*> mk1Ds = MakeTH1Ds(kin_value, dataMnv, hisy_temp, prefix);


  MnvPlotter plotter;
  plotter.ApplyStyle(kCCQENuStyle);
  plotter.axis_maximum = 5;
  for (auto& it : mk1Ds) {
    TCanvas c;
    //    plotter.DrawErrorSummary(it, "TL", false, true, 0.00001,0,"",1, "Uncertainities",false);
    plotter.DrawErrorSummary(it, "TR", true, true, -1, false, "", true);
    c.SaveAs(Form("fig_bin_wise_%s.png", it->GetName()));
  }



//vector <MnvH1D*> mk1Ds;

/*
  for (unsigned int i = 0; i < mk1Ds.size(); i++){
    mk1Ds[i]->SetLineColor(dataMnv->GetLineColor());
    mk1Ds[i]->SetLineStyle(dataMnv->GetLineStyle());
    mk1Ds[i]->SetLineWidth(dataMnv->GetLineWidth());
    mk1Ds[i]->SetMarkerStyle(dataMnv->GetMarkerStyle());
    mk1Ds[i]->SetMarkerSize(dataMnv->GetMarkerSize());
    mk1Ds[i]->GetXaxis()->SetNdivisions(504);
    mk1Ds[i]->GetYaxis()->SetNdivisions(504);
  }*/
  return mk1Ds;

}

void DrawBinRange(TH2* h, int axis, int bin, const char* varName, const char* numFormatStr=".2f", bool left=false, int gridx = 0, int gridy = 0) {

  double varmin = h->GetYaxis()->GetBinLowEdge(bin);
  double varmax = h->GetYaxis()->GetBinUpEdge(bin);
  //double binwidth = h->GetYaxis()->GetBinWidth(bin);

  //  cout << " Ptmu Bin " << bin << " : " << varmin << " " << binwidth << " " << varmin+binwidth << " " << varmax << endl;

  TString formatStr(TString::Format("%%%s < %%s < %%%s", numFormatStr, numFormatStr));

  TLatex* la=0;
  TString text(TString::Format(formatStr.Data(), varmin, varName, varmax));

  if(left){
    la=new TLatex(gPad->GetLeftMargin()+0.02,
                  1-gPad->GetTopMargin()-0.01,
                  text);
    la->SetTextAlign(13); // top left
  }
  else{
    la=new TLatex(1-gPad->GetRightMargin()-0.01,
                  1-gPad->GetTopMargin()-0.01,
                  text);
    la->SetTextAlign(33); // top right
  }

  la->SetNDC();
  la->SetTextFont(42);
  if(gridx==0||gridy==1) la->SetTextSize(0.03);
  else if(gridx<gridy) la->SetTextSize(28/(150.0*gridy));//Tall
  else if(gridx==gridy) la->SetTextSize(10/(150.0*gridx));//Wide
  else la->SetTextSize(18/(150.0*gridx));//Wide
  la->Draw();
  
}

//_____________________________________________________________
void makePlots(bool doMultipliers,bool doGenies,string location,string varx, string vary)
{
  ROOT::Cintex::Cintex::Enable();
  myPlotStyle();
  TH1::AddDirectory(false);
  TH1::SetDefaultSumw2();
  gStyle->SetErrorX(0);
  gStyle->SetEndErrorSize(2);

  string var = varx+"_"+vary;
  string varspecfilename = var;
  string varspec = var;
  if(vary=="ptmu") varspecfilename=varx;
  cout << varspecfilename << endl;
  string filename = Form("%s/CrossSection_CombinedPlaylists.root_%s.root",location.c_str(),varspecfilename.c_str());//Final result

  cout << filename << endl;


  if(vary!="ptmu")filename = Form("%s/CrossSection_%s_CombinedPlaylists.root",location.c_str(),varspecfilename.c_str());//Final result
  TFile f1(filename.c_str());
  //f1.ls();

  MnvH2D* hisy_temp_mnv   = (MnvH2D*)f1.Get(Form("h_%s_qelike_templatebinning",varspec.c_str()));//Get from N track
  MnvH1D* dataMnv         = (MnvH1D*)f1.Get(Form("h_%s_data_nobck_unfold_effcor_cross_section",varspec.c_str()));
 
  TH2* hisy_temp          = new TH2D(hisy_temp_mnv->GetCVHistoWithStatError());
  vector <MnvH1D*> data1D = PrepareHist(hisy_temp, dataMnv, "data1D");



  TFile *file = new TFile(Form("MnvHistOutputFiles_%s_%s.root",varx.c_str(),vary.c_str()),"RECREATE");
  file->cd();
  for (auto& it : data1D) {

   it->Write();
  }



  //name =  dataMnv->GetErrorBandNames();

//  cout << dataMnv->GetErrorBandNames() << endl;

vector <string> var_eb = dataMnv->GetErrorBandNames();   
 for(unsigned int i = 0; i < var_eb.size(); i++) {
 cout << var_eb[i].c_str() <<endl;

 }
cout << " ------ Lat ------" << endl;
vector <string> var_leb = dataMnv->GetLatErrorBandNames();   
 for(unsigned int i = 0; i < var_leb.size(); i++) {
 cout << var_leb[i].c_str() <<endl;

 }
cout << " ------ Vet ------" << endl;
vector <string> var_veb = dataMnv->GetVertErrorBandNames();   
 for(unsigned int i = 0; i < var_veb.size(); i++) {
 cout << var_veb[i].c_str() <<endl;

 }





  //dataMnv->Scale(1e-6, "width");
  //dataMnv->SetMaximum(1);
  MnvPlotter plotter;
  plotter.ApplyStyle(kCCQENuStyle);
  plotter.axis_maximum = 5;
  TCanvas c;
  plotter.DrawErrorSummary(dataMnv, "TL", false, true, 0.00001,0,"",1, "Uncertainities",false);
  //plotter.DrawErrorSummary(dataMnv, "TR", true, true, -1, false, "", true);
  c.SaveAs("test.png");

  /*
  MnvH1D* mcMnv             = (MnvH1D*)f1.Get(Form("h_%s_mc_nobck_unfold_effcor_cross_section",varspec.c_str()));
  MnvH1D* nomGenieMnv       = (MnvH1D*)f1.Get(Form("h_%s_mc_nobck_unfold_effcor_cross_section",varspec.c_str()));//2
  MnvH1D* bestGenieMnv      = (MnvH1D*)f1.Get(Form("h_%s_mc_nobck_unfold_effcor_cross_section",varspec.c_str()));//3
  MnvH1D* mcMnv_qelike_qe   = (MnvH1D*)f1.Get(Form("h_%s_cross_section_qelike_qe",varspec.c_str()));//Get from N track
  MnvH1D* mcMnv_qelike_res  = (MnvH1D*)f1.Get(Form("h_%s_cross_section_qelike_res",varspec.c_str()));//Get from N track
  MnvH1D* mcMnv_qelike_dis  = (MnvH1D*)f1.Get(Form("h_%s_cross_section_qelike_dis",varspec.c_str()));//Get from N track
  MnvH1D* mcMnv_qelike_2p2h = (MnvH1D*)f1.Get(Form("h_%s_cross_section_qelike_2p2h",varspec.c_str()));//Get from N track
  MnvH2D* hisy_temp_mnv     = (MnvH2D*)f1.Get(Form("h_%s_qelike_templatebinning",varspec.c_str()));//Get from N track
  */

  /*
  dataMnv->Scale(1e39, "width");
  mcMnv->Scale(1e39, "width");
  nomGenieMnv->Scale(1e39, "width");
  bestGenieMnv->Scale(1e39, "width");

  mcMnv_qelike_qe->Scale(1e39,"width");
  mcMnv_qelike_res->Scale(1e39,"width");
  mcMnv_qelike_dis->Scale(1e39,"width");
  mcMnv_qelike_2p2h->Scale(1e39,"width");
  hisy_temp_mnv->Scale(1e39,"width");
  */
 /*
  TH2* hisy_temp      = new TH2D(hisy_temp_mnv->GetCVHistoWithStatError());
 
  TH1* dataStat       = new TH1D(dataMnv->GetCVHistoWithStatError());
  
  TH1* data           = new TH1D(dataMnv->GetCVHistoWithError());
  TH1* mc             = new TH1D(mcMnv->GetCVHistoWithStatError());
  TH1* nomGenie       = new TH1D(nomGenieMnv->GetCVHistoWithStatError());
  TH1* bestGenie      = new TH1D(bestGenieMnv->GetCVHistoWithStatError());
  TH1* mc_qelike_qe   = new TH1D(mcMnv_qelike_qe->GetCVHistoWithStatError());
  TH1* mc_qelike_res  = new TH1D(mcMnv_qelike_res->GetCVHistoWithStatError());
  TH1* mc_qelike_dis  = new TH1D(mcMnv_qelike_dis->GetCVHistoWithStatError());
  TH1* mc_qelike_2p2h = new TH1D(mcMnv_qelike_2p2h->GetCVHistoWithStatError());

  vector<int> mycolors = getColors(2);
  mc->SetLineColor(kRed);
  mc->SetLineWidth(2);
  nomGenie->SetLineColor(kBlue);
  nomGenie->SetLineWidth(2);
  bestGenie->SetLineColor(mycolors[10]);
  bestGenie->SetLineWidth(2);
  mc_qelike_qe->SetLineColor(mycolors[3]);
  mc_qelike_res->SetLineColor(mycolors[4]);
  mc_qelike_dis->SetLineColor(mycolors[5]);
  mc_qelike_2p2h->SetLineColor(mycolors[16]);
  data->SetMarkerStyle(kFullCircle);
  data->SetMarkerSize(0.6);
  //data->SetMarkerStyle(1);
  data->SetLineColor(kBlack);
  data->SetLineWidth(2);
  
  dataStat->SetLineColor(kBlack);
  dataStat->SetMarkerStyle(1);
  
  vector <TH1*> dataStat1D       = PrepareHist(hisy_temp, dataStat, "dataStat1D");
  vector <TH1*> data1D           = PrepareHist(hisy_temp, data, "data1D");
  vector <TH1*> mc1D             = PrepareHist(hisy_temp, mc, "mc1D");
  vector <TH1*> nomGenie1D       = PrepareHist(hisy_temp, nomGenie, "nomGenie1D");
  vector <TH1*> bestGenie1D      = PrepareHist(hisy_temp, bestGenie, "bestGenie1D");
  vector <TH1*> mc_qelike_qe1D   = PrepareHist(hisy_temp, mc_qelike_qe, "mc_qelike_qe1D");
  vector <TH1*> mc_qelike_res1D  = PrepareHist(hisy_temp, mc_qelike_res, "mc_qelike_res1D");
  vector <TH1*> mc_qelike_dis1D  = PrepareHist(hisy_temp, mc_qelike_dis, "mc_qelike_dis1D");
  vector <TH1*> mc_qelike_2p2h1D = PrepareHist(hisy_temp, mc_qelike_2p2h, "mc_qelike_2p2h1D");
  
  vector <vector <TH1*>> allHists;
  allHists.push_back(dataStat1D);
  allHists.push_back(mc1D);
  if(doGenies){  
    allHists.push_back(nomGenie1D);
    allHists.push_back(bestGenie1D);
  } else {
    allHists.push_back(mc_qelike_qe1D);
    allHists.push_back(mc_qelike_res1D);
    allHists.push_back(mc_qelike_dis1D);
    allHists.push_back(mc_qelike_2p2h1D);
  }
  allHists.push_back(data1D);

  vector <double>vptmubw;
  for(unsigned int i = 0; i < dataStat1D.size(); i++) {
    int bin = i + 1;
    double varmin = hisy_temp->GetYaxis()->GetBinLowEdge(bin);
    double varmax = hisy_temp->GetYaxis()->GetBinUpEdge(bin);
    double ptmu_binwidth = hisy_temp->GetYaxis()->GetBinWidth(bin);
    vptmubw.push_back(ptmu_binwidth);
    cout << " Ptmu Bin " << i << " : " << varmin << " " << ptmu_binwidth << " " << varmin+ptmu_binwidth << " " << varmax << endl;
  }
  

  //Normalization, each bin by (1e39/ptmu,width)
  for (auto& it : allHists) {
    int j = 0;
    for (auto& itt : it) {
      //itt->Scale((1e39/vptmubw[j]),"width");
      itt->Scale((1e39),"width");
      j++;
    }
  }

  vector <double> multi = GetScales(allHists);


  TFile *file = new TFile(Form("HistOutputFiles_%s_%s.root",varx.c_str(),vary.c_str()),"RECREATE");
  file->cd();
  for (auto& it : allHists) {
    for (auto& itt : it) {
      
   itt->Write();   
 }
  }
  
  hisy_temp->Write();
  dataStat->Write();
  data->Write();
  mc->Write();
  nomGenie->Write();
  bestGenie->Write();
  mc_qelike_qe->Write();
  mc_qelike_res->Write();
  mc_qelike_dis->Write();
  mc_qelike_2p2h->Write();
  file->Close();
  

  int nhist = dataStat1D.size();
  int grid_x = sqrt(nhist)+1;
  int grid_y = nhist/(grid_x-1);
  if(grid_x*grid_y-nhist==grid_x) grid_y-=1;

  GridCanvas* gc=new GridCanvas(uniq(), grid_x, grid_y, 1000, 700);
  gc->SetRightMargin(0.01);
  gc->SetLeftMargin(0.1);
  //  gc->SetBottomMargin(0.1);
  gc->ResetPads();

  string cell_title = varLatex.at(vary);

  string ytitle = "d^{2}#sigma/"+varLatex.at(varx)+varLatex.at(vary)+" (x10^{-39} cm^{2}/GeV^{2}/c^{2}/Nucleon)";
  gc->SetYTitle(ytitle.c_str());
  gc->SetXTitle(varLatex.at(varx).c_str());

  gc->Modified();

  for(int i=0; i<nhist; ++i){
    TPad* pad=(TPad*)gc->cd(i+1);
    //dataStat1D[i]->GetYaxis()->SetRangeUser(0.0, 4);

    if (doMultipliers) {
      dataStat1D[i]->Scale(multi[i]);
      mc1D[i]->Scale(multi[i]);
      
      nomGenie1D[i]->Scale(multi[i]);
      bestGenie1D[i]->Scale(multi[i]);
      
      mc_qelike_qe1D[i]->Scale(multi[i]);
      mc_qelike_res1D[i]->Scale(multi[i]);
      mc_qelike_dis1D[i]->Scale(multi[i]);
      mc_qelike_2p2h1D[i]->Scale(multi[i]);
      
      data1D[i]->Scale(multi[i]);
    }

    dataStat1D[i]->Draw("histpe1");
    mc1D[i]->Draw("hist l same");
  
    if(doGenies){
      nomGenie1D[i]->Draw("hist l same");
      bestGenie1D[i]->Draw("hist l same");
    }
    else {
      mc_qelike_qe1D[i]->Draw("hist l same");
      mc_qelike_res1D[i]->Draw("hist l same");
      mc_qelike_dis1D[i]->Draw("hist l same");
      mc_qelike_2p2h1D[i]->Draw("hist l same");
    }
    data1D[i]->Draw("HIST PE1 same");
        
    DrawBinRange(hisy_temp, 1, i+1, cell_title.c_str(), ".2f", false);
    
    if(doMultipliers){
      TLatex* la2=new TLatex(1-pad->GetRightMargin()-0.01,
			     1-pad->GetTopMargin()-0.06,
			     TString::Format("#times %.2f", multi[i]));
      la2->SetTextAlign(33); // top right
      la2->SetNDC();
      la2->SetTextFont(42);
      la2->SetTextSize(0.035);
      la2->Draw();
    } 
    
  }
  
  TLegend* leg=new TLegend(0.82, 0.1, 0.9, 0.35);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.03);
  leg->AddEntry(data, "MINERvA data", "lpe");
  leg->AddEntry(mc, "MnvGENIE", "l");
  if(doGenies){
    leg->AddEntry(nomGenie, "GENIE 2.8.4", "l");
    leg->AddEntry(bestGenie, "2p2h and RPA", "l");
  }
  else{
    leg->AddEntry(mc_qelike_qe,"QE","l");
    leg->AddEntry(mc_qelike_res,"Resonant","l");
    leg->AddEntry(mc_qelike_dis,"DIS","l");
    leg->AddEntry(mc_qelike_2p2h,"2p2h","l");
    
  }
  leg->Draw();

  gc->SetYLimits(0, 6.99);
  gc->Modified();


  //  gc->Print("test.png");

  if(doGenies){
    gc->Print(doMultipliers ? Form("nu-2d-xsec-genies-nongrid-%s-multiplier.png",var.c_str()) : Form("nu-2d-xsec-genies-nongrid-%s.png",var.c_str()));
  }
  else {
    gc->Print(doMultipliers ? Form("nu-2d-xsec-comps-nongrid-%s-multiplier.png",var.c_str()) : Form("nu-2d-xsec-comps-nongrid-%s.png",var.c_str()));

  }
*/

}

int main(int argc, char* argv[])
{

  InitVarLatex();
  makePlots(true,false,argv[1],argv[2],argv[3]);
  return 0;
}
