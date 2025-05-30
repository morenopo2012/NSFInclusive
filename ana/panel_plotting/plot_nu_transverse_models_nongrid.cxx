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
#include "PlotUtils/MnvColors.h"
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
      //      cout << i << " " << j << "  "<< binid << " " << mapped_low << " " << xindex  << " " << kin_val << endl; 
    } 
    //it = unique(binsptbin.begin(), binsptbin.end());
    //binsptbin.resize(distance(binsptbin.begin(), it));
    
    std::sort(binsptbin.begin(), binsptbin.end()); 
    auto last = std::unique(binsptbin.begin(), binsptbin.end());
    binsptbin.erase(last, binsptbin.end());

    /*
    cout << " Bin " << i << " : ";
    for(auto it:binsptbin) {
      
      cout <<it << " ";
    }
    cout << endl;
    */
    tmpv.push_back(binsptbin);
  }
  return tmpv;
}

TH1* makeTH1D(int pt_bin,double pt_width, vector<double> boundaries,TH1* rebinnedhisto,int startindex,const char *prefix) {
  
  int n = boundaries.size();
  


  double arr[n];
  //cout << " Ptbin " << pt_bin << " " << startindex << " ";

  for (unsigned int i = 0; i < boundaries.size(); i++) {
    arr[i] = boundaries[i];
    //cout  << arr[i] << " ";

  }
  
  TH1* retTH1D = new TH1D(Form("%s_pt_bin_%d",prefix,pt_bin),Form("%s_pt_bin_%d", prefix,pt_bin),n-1,arr);
  
  for (unsigned int i = 0; i < boundaries.size(); i++) {
    
    double bv      = rebinnedhisto->GetBinContent(startindex+i)/(pt_width);
    double bv_err  = rebinnedhisto->GetBinError(startindex+i)/(pt_width);
//    double bv      = rebinnedhisto->GetBinContent(startindex+i);
//    double bv_err  = rebinnedhisto->GetBinError(startindex+i);
    retTH1D->SetBinContent(i+1,bv);
    retTH1D->SetBinError(i+1,bv_err);
  }
  //cout << " " <<retTH1D->GetNbinsX() << endl;

  return retTH1D;
}

//_________________________________
vector <TH1*> MakeTH1Ds(vector<vector<double>> kinVs, TH1 *rebinnedhisto,TH2 *templatehist, const char *prefix) {
  vector <TH1*> myTH1Ds;
  int counter = templatehist->GetNbinsX()+4;
  int i = 0;
  //  cout <<"--------------------------------------------"<<endl;
  for (auto& it : kinVs) {
    
    TH1* tmp_th1D = makeTH1D(i,templatehist->GetYaxis()->GetBinWidth(i+1),it,rebinnedhisto,counter,prefix);
    counter+=it.size()+1;
    myTH1Ds.push_back(tmp_th1D);

    /*
      cout << " Bin " << i <<" : ";
      for (auto& itt : it) {
      cout << itt << " ";
      }
      cout <<endl;
    */
    
    i++;
  }
  
  return myTH1Ds;
}

vector<double> GetScalesNew(vector<vector<TH1*>>histopts, bool pxProj, double plotMax, double fracMax, bool limitRange = false, double rmin=0, double rmax=100){
  vector<double> tmpvect;
  
  int nh = histopts.size();
  int nb = histopts[0].size();
  
  double mxv[nb][nh];
  
  int i = 0;
  for (auto& it : histopts) {
    int bin = 0;
    for (auto& itt : it) {
      int maxbin = itt->GetMaximumBin();
      mxv[bin][i] = itt->GetBinContent(maxbin);
      //double content = itt->GetBinContent(maxbin);
      bin++;
    }
    i++;
  }
  
  for(int ii = 0; ii < nb; ii++){
    double maxval = 0;
    for(int jj = 0; jj < nh; jj++){
      double content =  mxv[ii][jj];
      if (content > maxval) maxval = content;
    } 
    double scale = (plotMax*fracMax)/maxval;
    if(scale>1){
      int tmpscale = floor(scale*10);
      scale = tmpscale/10.0;
    }
    else{
      int tmpscale = ceil(scale*10);
      scale = tmpscale/10.0;
    }
    //scale = scale;
    //cout << "Bin " << ii << " Scale Factor " << scale << endl;
    tmpvect.push_back(scale);
  }
  
  return tmpvect;
}


vector<double> GetScalesRatio(vector<vector<TH1*>> histopts, bool pxProj, double plotMax){
  vector<double> tmpvect;
  
  int nh = histopts.size();
  int nb = histopts[0].size();

  cout << nh << " " << nb << endl;
  
  double mxv[nb][nh];
  
  int i = 0;
  for (auto& it : histopts) {
    int bin = 0;
    for (auto& itt : it) {
      int maxbin = itt->GetMaximumBin();
      mxv[bin][i] = itt->GetBinContent(maxbin);
      //double content = itt->GetBinContent(maxbin);
      bin++;
    }
    i++;
  }
  
  for(int ii = 0; ii < nb; ii++){
    double maxval = 0;
    for(int jj = 0; jj < nh; jj++){
      double content =  mxv[ii][jj];
      if (content > maxval) maxval = content;
    } 
    //we want about fracMax of the plotMax. 
    //This is a ratio -> we only want easy numbers
    double scale = (plotMax*0.75)/maxval;
    if(scale>1) scale = 1;
    if(scale<1 && scale >0.5) scale = 0.5;
    if(scale<0.5) scale = 0.25;
    tmpvect.push_back(scale);
  }
  
  return tmpvect;
}

vector<double> GetScales(vector <vector <TH1*>> histopts){
  vector<double> tmpvect;
  //cout  << " <<<<<<< Calculating the Scale Factor >>>>>>>>>> " << endl;
  int nh = histopts.size();
  int nb = histopts[0].size();
  
  double mxv[nb][nh];
  
  // cout << nh  << "  " << nb << endl;
  int i = 0;
  for (auto& it : histopts) {
    int bin = 0;
    for (auto& itt : it) {
      int maxbin = itt->GetMaximumBin();
      mxv[bin][i] = itt->GetBinContent(maxbin);
      //cout << itt->GetName() << " " << i << "  " << bin << " " << mxv[bin][i] << " " << maxbin  << endl;
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
    //    cout << "Bin " << ii << " Scale Factor " << scale << endl;
    tmpvect.push_back(scale);
  }


  // cout << " <<<<<<<<<<< ------------- end scale value calculation  --------- >>>>>>>>> " << endl;
  return tmpvect;
}

//___________________________________________________________________________
vector <TH1*> PrepareHist(TH2* hisy_temp, TH1* dataMnv, const char *prefix){
  int n_binsx_original = hisy_temp->GetNbinsX()+2;
  int n_binsy_original = hisy_temp->GetNbinsY()+2;
  // cout << "==================== " << n_binsx_original << "  " << n_binsy_original <<endl;
  vector<vector<double>> kin_value = GetKinValue(n_binsx_original, n_binsy_original, dataMnv, hisy_temp);
  vector <TH1*> mk1Ds = MakeTH1Ds(kin_value, dataMnv, hisy_temp, prefix);


  for (unsigned int i = 0; i < mk1Ds.size(); i++){
   mk1Ds[i]->SetLineColor(dataMnv->GetLineColor());
   mk1Ds[i]->SetLineStyle(dataMnv->GetLineStyle());
   mk1Ds[i]->SetLineWidth(dataMnv->GetLineWidth());
   mk1Ds[i]->SetMarkerStyle(dataMnv->GetMarkerStyle());
   mk1Ds[i]->SetMarkerSize(dataMnv->GetMarkerSize());
    mk1Ds[i]->GetXaxis()->SetNdivisions(504);
    mk1Ds[i]->GetYaxis()->SetNdivisions(504);
  }

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
void makePlots(bool doMultipliers,bool doGenies,string location,string varx, string vary, string iiter)
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
  
  TString model[] = {"CV", "CV_ISI_NuWro_LFG", "CV_ISI_NuWro_SF", "CV_minerva_joint_lowq2", "CV_MK", "CV_Zexp_Bubble", "default", "pion_2p2h", "pion_rpa", "piontune"};
  TString model_name[] = {"MINERvA Tune V1", "ISI NuWro LFG", "ISI NuWro SF", "MINERvA joint lowq2", "Minoo Model", "V4", "GENIE 2.12.6 + Valencia 2p2h", "Pion + 2p2h", "Pion + RPA", "Pion Tune"};
  MnvH1D* dataMnv[10];
  MnvH1D* mcMnv[10];
  
  /*MnvH1D* mcMnv_qelike_qe[10];
  MnvH1D* mcMnv_qelike_res[10];
  MnvH1D* mcMnv_qelike_dis[10];
  MnvH1D* mcMnv_qelike_2p2h[10];*/
  
  MnvH2D* hisy_temp_mnv[10];
  
  TH2* hisy_temp[10];
  TH1* dataStat[10];      
  TH1* data[10];          
  TH1* mc[10];            
  /*TH1* nomGenie[10];      
    TH1* bestGenie[10];     
    TH1* mc_qelike_qe[10];  
    TH1* mc_qelike_res[10]; 
    TH1* mc_qelike_dis[10]; 
    TH1* mc_qelike_2p2h[10];*/
  
  vector <TH1*> dataStat1D[10];
  vector <TH1*> data1D[10];
  vector <TH1*> mc1D[10];
  
  for (int i = 0; i < 10; i++) {
    
    string filename = Form("%s/Transverse_%s/CrossSection_iter_%s_CombinedPlaylists.root_%s.root",location.c_str(),model[i].Data(),iiter.c_str(), varspecfilename.c_str());//Final result
    cout << filename << endl;
    if(vary!="ptmu")filename = Form("%s/Transverse_%s/CrossSection_iter_%s_%s_CombinedPlaylists.root",location.c_str(),model[i].Data(),iiter.c_str(),varspecfilename.c_str());//Final result
    TFile f1(filename.c_str());
    
    dataMnv[i]           = (MnvH1D*)f1.Get(Form("h_%s_data_nobck_unfold_effcor_cross_section",varspec.c_str()));
    mcMnv[i]             = (MnvH1D*)f1.Get(Form("h_%s_mc_nobck_unfold_effcor_cross_section",varspec.c_str()));
    /*mcMnv_qelike_qe[i]   = (MnvH1D*)f1.Get(Form("h_%s_cross_section_qelike_qe",varspec.c_str()));//Get from N track
      mcMnv_qelike_res[i]  = (MnvH1D*)f1.Get(Form("h_%s_cross_section_qelike_res",varspec.c_str()));//Get from N track
      mcMnv_qelike_dis[i]  = (MnvH1D*)f1.Get(Form("h_%s_cross_section_qelike_dis",varspec.c_str()));//Get from N track
      mcMnv_qelike_2p2h[i] = (MnvH1D*)f1.Get(Form("h_%s_cross_section_qelike_2p2h",varspec.c_str()));//Get from N track
    */      
    hisy_temp_mnv[i]     = (MnvH2D*)f1.Get(Form("h_%s_qelike_templatebinning",varspec.c_str()));//Get from N track
    
    hisy_temp[i]         = new TH2D(hisy_temp_mnv[i]->GetCVHistoWithStatError());
    dataStat[i]          = new TH1D(dataMnv[i]->GetCVHistoWithStatError());
    data[i]              = new TH1D(dataMnv[i]->GetCVHistoWithError());
    mc[i]                = new TH1D(mcMnv[i]->GetCVHistoWithStatError());
    /* 
    mc_qelike_qe[i]      = new TH1D(mcMnv_qelike_qe->GetCVHistoWithStatError());
    mc_qelike_res[i]     = new TH1D(mcMnv_qelike_res->GetCVHistoWithStatError());
    mc_qelike_dis[i]     = new TH1D(mcMnv_qelike_dis->GetCVHistoWithStatError());
    mc_qelike_2p2h[i]    = new TH1D(mcMnv_qelike_2p2h->GetCVHistoWithStatError());
    */
    
    
    dataStat1D[i]        = PrepareHist(hisy_temp[i], dataStat[i], Form("dataStat1D_%d",i));
    data1D[i]            = PrepareHist(hisy_temp[i], data[i], Form("data1D_%d",i));
    mc1D[i]              = PrepareHist(hisy_temp[i], mc[i], Form("mc1D_%d",i));

   
    for (int k = 0; k < dataStat1D[i].size(); k++) {
      dataStat1D[i][k]->Scale(1e39,"width");
    }

    for (int k = 0; k < data1D[i].size();k++) {
      data1D[i][k]->Scale(1e39,"width");
      data1D[i][k]->SetMarkerSize(0.8);
    }

    for (int k = 0; k < mc1D[i].size();k++) {
      mc1D[i][k]->Scale(1e39,"width");
    }
    
    

    
    /*vector <TH1*> nomGenie1D       = PrepareHist(hisy_temp, nomGenie, "nomGenie1D");
      vector <TH1*> bestGenie1D      = PrepareHist(hisy_temp, bestGenie, "bestGenie1D");
      vector <TH1*> mc_qelike_qe1D   = PrepareHist(hisy_temp, mc_qelike_qe, "mc_qelike_qe1D");
      vector <TH1*> mc_qelike_res1D  = PrepareHist(hisy_temp, mc_qelike_res, "mc_qelike_res1D");
      vector <TH1*> mc_qelike_dis1D  = PrepareHist(hisy_temp, mc_qelike_dis, "mc_qelike_dis1D");
      
      vector <TH1*> mc_qelike_2p2h1D = PrepareHist(hisy_temp, mc_qelike_2p2h, "mc_qelike_2p2h1D");
    */
    
    
  }
  vector <vector <TH1*>> allHists;
  //allHists.push_back(dataStat1D[0]);
  allHists.push_back(data1D[0]);
  
  for (int i = 0; i < 10; i++) {
    allHists.push_back(mc1D[i]);  
  }
  vector <TH1*> mc_denominator;
  for (auto& it : mc1D[0]) {
    const char *name_tmp = Form("%s_denominator",it->GetName());
    TH1* h = (TH1*)it->Clone(name_tmp);
    mc_denominator.push_back(h);
  }
  
  
  
  
  //Normalization, each bin by (1e39/ptmu,width)
  
  
  
  vector <double> multi = GetScalesRatio(allHists,1,2.5);

  for (auto& it : multi){
    cout << it <<endl;
  }

  
 // vector<int> mycolors = MnvColors::GetColors(9);

  /*
    mycolors.push_back(1052);
    mycolors.push_back(1053);
    mycolors.push_back(1053);
    mycolors.push_back(1054);
    mycolors.push_back(219);
    mycolors.push_back(2);
    mycolors.push_back(4);
    mycolors.push_back(kGreen+2);
    mycolors.push_back(kOrange+2);
    mycolors.push_back(kViolet+1);
    mycolors.push_back(3);
    mycolors.push_back(kRed+2);
  */
  int mycolors[] = {2, 4, 3, 219, 810, 418, 616, 810, 860, 880, 862, 882, 800, 810};
  
  //for (auto& it : mycolors){
  //  cout << it <<endl;
  //}
  
  
  int nhist = dataStat1D[0].size();
  int grid_x = sqrt(nhist)+1;
  int grid_y = nhist/(grid_x-1);
  if(grid_x*grid_y-nhist==grid_x) grid_y-=1;
  
  if (nhist < 16) {
    grid_x = 5;
    grid_y = 3;
  }


  GridCanvas* gc=new GridCanvas(uniq(), grid_x, grid_y, 1000, 700);
  gc->SetRightMargin(0.01);
  gc->SetLeftMargin(0.1);
  //  gc->SetBottomMargin(0.1);
  gc->ResetPads();
  
  string cell_title = varLatex.at(vary);
  string ytitle = "Ratio to MINERvA Tune - v1";
  //string ytitle = "d^{2}#sigma/"+varLatex.at(varx)+varLatex.at(vary)+" (x10^{-39} cm^{2}/GeV^{2}/c^{2}/Nucleon)";
  gc->SetYTitle(ytitle.c_str());
  gc->SetXTitle(varLatex.at(varx).c_str());
  
  gc->Modified();
  
  for(int i=0; i<nhist; ++i){
    //TPad* pad=(TPad*)
    gc->cd(i+1);
    dataStat1D[0][i]->Divide(mc_denominator[i]);
    for(int j = 0; j < 10; j++) {
      mc1D[j][i]->Divide(mc_denominator[i]);
    }
    data1D[0][i]->Divide(mc_denominator[i]);
    
    //dataStat1D[i]->GetYaxis()->SetRangeUser(0.0, 4);
    
    /*    dataStat1D[0][i]->Scale(multi[i]);
	  for(int j = 0; j < 10; j++) 
	  mc1D[j][i]->Scale(multi[i]);      
	  
	  data1D[0][i]->Scale(multi[i]);
    */
    //    dataStat1D[0][i]->Draw("histpe1");
    
    data1D[0][i]->Draw("histpe");
    for(int j = 0; j < 10; j++) {
      mc1D[j][i]->SetLineColor(mycolors[j]);
      mc1D[j][i]->Draw("hist l same");
    }
    data1D[0][i]->Draw("histpe same");
    
    
    DrawBinRange(hisy_temp[0], 1, i+1, cell_title.c_str(), ".2f", false);
    
    // if(doMultipliers){
    /*
      TLatex* la2=new TLatex(1-pad->GetRightMargin()-0.01,1-pad->GetTopMargin()-0.06,
      TString::Format("#times %.2f", multi[i]));
      la2->SetTextAlign(33); // top right
      la2->SetNDC();
      la2->SetTextFont(42);
      la2->SetTextSize(0.035);
      la2->Draw();
      //  } 
      */
  }
  
  
    TLegend* leg=new TLegend(0.32, 0.1, 0.9, 0.35);
    leg->SetNColumns(2);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->SetTextSize(0.03);
    leg->AddEntry(data1D[0][0], "MINERvA data", "lpe");
    for(int j = 0; j < 10; j++) {
      leg->AddEntry(mc1D[j][0], model_name[j].Data(), "l");
    }
    leg->Draw();
    
    //gc->SetYLimits(0.49,1.49);
    gc->SetYLimits(-0.99,2.99);
   
  /*  if(varx=="dpt") gc->SetXLimits(0,1.99);
    if(varx=="dpl") gc->SetXLimits(-0.49,0.49);
    if(varx=="dpty") gc->SetXLimits(-1.69,0.89);
    if(varx=="dpn") gc->SetXLimits(0,1.99);
    if(varx=="tp") gc->SetXLimits(0.09,0.49);
*/
    gc->Modified();
    
    
    //  gc->Print("test.png");
    
    if(doGenies){
      gc->Print(doMultipliers ? Form("nu-ratio-nongrid-iter-%s-%s-multiplier.png",iiter.c_str(), var.c_str()) : Form("nu-ratio-nongrid-iter-%s-%s.png",iiter.c_str(), var.c_str()));
    }
    else {
      gc->Print(doMultipliers ? Form("nu-ratio-nongrid-iter-%s-%s-multiplier.png",iiter.c_str(), var.c_str()) : Form("nu-ratio-nongrid-iter-%s-%s.png",iiter.c_str(), var.c_str()));
      
    }
    
}

int main(int argc, char* argv[]){
  InitVarLatex();
  string a = argv[4];
  makePlots(true,false,argv[1],argv[2],argv[3],a);
  return 0;
}
