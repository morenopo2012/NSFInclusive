//#include "myPlotStyle.h"
#include "GridCanvas.h"
#include "PlotUtils/MnvH2D.h"
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
#include "localColor.h"
#include "PlotUtils/MnvColors.h"
//#include "Cintex/Cintex.h"
#include <THStack.h>
#include "myPlotStyle.h"

#include "plot.h"
#include <algorithm>
using namespace PlotUtils;
using namespace std;
void makePlots(bool doMultipliers,bool doRatio,string location, const std::string& varName,  int targetID, int targetZ,bool ispT=false,bool IsTargetCombined=false)
{
 // ROOT::Cintex::Cintex::Enable();
  myPlotStyle();
  TH1::AddDirectory(false);
  TH1::SetDefaultSumw2();
  //gStyle->SetErrorX(0);
//  gStyle->SetEndErrorSize(2);
//  gStyle->SetLabelSize(0.04);
  TFile* f1=nullptr;
  if(targetID==99)
	  f1=TFile::Open(Form("EventSelection_fullp4_ns_thesis_tuneV4_t%d_z%d_sys.root",targetID,targetZ));
 if(targetZ==06)
	    f1=TFile::Open(Form("EventSelection_fullp4_ns_thesis_tuneV4_t3_z06_sys.root"));
  if(IsTargetCombined==false && targetZ!=99 && targetZ!=06) 
  f1=TFile::Open(Form("EventSelection_single_fullp4_ns_thesis_tuneV4_t%d_z%02d_sys.root",targetID,targetZ));
 // f1=TFile::Open(Form("EventSelection_fullp4_ns_thesis_tuneV4_t%d_z%02d_sys.root",targetID,targetZ)); 
if(IsTargetCombined==true)
	  f1=TFile::Open(Form("EventSelection_fullp4_ns_targetCombined_not_target1_thesis_tuneV4_t%d_z%02d_sys.root",targetID,targetZ));
// f1=TFile::Open(Form("EventSelection_fullp4_ns_targetCombined_thesis_tuneV4_t%d_z%02d_sys.root",targetID,targetZ));
//if(targetID==99)
//	f1=TFile::Open(Form("EventSelection_fullp4_ns_thesis_tuneV4_t%d_z%d_sys.root",targetID,targetZ));
//TFile f1(Form("EventSelection_fullp4_ns_targetCombined_thesis_tuneV4_t%d_z%02d_sys.root",targetID,targetZ));	
  MnvH2D* mcMnv=(MnvH2D*)f1->Get(Form("selected_mc_reco2d_%s",varName.c_str()));
  MnvH2D* dataMnv=(MnvH2D*)f1->Get(Form("selected_data_reco2d_%s",varName.c_str()));
  MnvH2D* mcMnv_signal=(MnvH2D*)f1->Get(Form("selected_mc_reco2d_signal_%s",varName.c_str()));
  MnvH2D* QE=(MnvH2D*)f1->Get(Form("selected_mc_reco2d_QE_%s",varName.c_str()));
  MnvH2D* RES=(MnvH2D*)f1->Get(Form("selected_mc_reco2d_RES_%s",varName.c_str()));
  MnvH2D* DIS =(MnvH2D*)f1->Get(Form("selected_mc_reco2d_DIS_%s",varName.c_str()));
  MnvH2D* mc_2p2h1=(MnvH2D*)f1->Get(Form("selected_mc_reco2d_2p2h_%s",varName.c_str()));
  MnvH2D* other=(MnvH2D*)f1->Get(Form("selected_mc_reco2d_OtherIT_%s",varName.c_str()));
  MnvH2D* bkg = (MnvH2D*)f1->Get(Form("selected_mc_reco2d_bkg_%s",varName.c_str()));


  double DataPOT = 1.11707110844e+21;
  double MCPOT = 4.96197407132e+21;
//  double DataPOT =5.24737e+19;
//  double MCPOT = 3.13434e+20;

  double scale = DataPOT/MCPOT;
  
   mcMnv->Scale(scale*1e-4, "width");
    mcMnv_signal->Scale(scale*1e-4, "width");
    QE->Scale(scale*1e-4, "width");
    RES->Scale(scale*1e-4, "width");
    DIS->Scale(scale*1e-4, "width");
    other->Scale(scale*1e-4, "width");
    mc_2p2h1->Scale(scale*1e-4, "width");
    bkg->Scale(scale*1e-4, "width");
   dataMnv->Scale(1e-4, "width");

 // if(targetID==99)
 // {  mcMnv->Scale(scale*1e-4, "width");
  //   dataMnv->Scale(1e-4, "width");
 // } 
  // Get the data histogram with stat error and with total error
  // separately so we can plot them both for inner and outer ticks
  TH2* dataStat=new TH2D(dataMnv->GetCVHistoWithStatError());
  TH2* data=new TH2D(dataMnv->GetCVHistoWithError());
  TH2* mc=new TH2D(mcMnv->GetCVHistoWithStatError());
  TH2* mc_line=new TH2D(mcMnv->GetCVHistoWithStatError());
  TH2* mcTotalError=new TH2D(mcMnv->GetCVHistoWithError());
    TH2* signal=new TH2D(mcMnv_signal->GetCVHistoWithStatError());
  TH2* mc_QE=new TH2D(QE->GetCVHistoWithError());
  TH2* mc_QE_line=new TH2D(QE->GetCVHistoWithError());
   TH2* mc_RES=new TH2D(RES->GetCVHistoWithError());
   TH2* mc_RES_line=new TH2D(RES->GetCVHistoWithError());
  TH2* mc_DIS=new TH2D(DIS->GetCVHistoWithError());
  TH2* mc_DIS_line=new TH2D(DIS->GetCVHistoWithError());
  TH2* mc_2p2h=new TH2D(mc_2p2h1->GetCVHistoWithError());
  TH2* mc_other=new TH2D(other->GetCVHistoWithError());
  TH2* mc_bkg = new TH2D(bkg->GetCVHistoWithError());
 
  THStack* hs = new THStack("hs","Stacked histograms");


/*  total_bkg->Reset();
  total_bkg->Add(bkg_DS);
  total_bkg->Add(bkg_US);
  total_bkg->Add(bkg_other);
  total_bkg->Add(bkg_WS);
  total_bkg->Add(bkg_NC);
*/  
  // These line and marker styles will be propagated to the 1D plots
  vector<int> mycolors = MnvColors::GetColors(2);
  if(!doRatio){
   mc->SetLineColor(kBlack);
   mc_line->SetLineColor(kBlack);
   mc_line->SetLineWidth(2.2);
   mcTotalError->SetFillColor(kBlue-10);
   mcTotalError->SetFillStyle(3144);
   mc->SetLineWidth(2.2);
   mc->SetFillColor(kBlue-7);
   mc_QE->SetLineColor(kGray+3);
   mc_QE->SetLineWidth(1.9);
  
  mc_other->SetLineColor(kGray+3);
  mc_other->SetLineWidth(1.9);
 mc_other->SetFillColor(kGreen-9);
  
  mc_QE_line->SetLineColor(kGray+3);
   mc_QE_line->SetLineWidth(1.9);

   //mc_QE->SetFillStyle(3002);
   mc_QE->SetFillColor(kBlue-7);

   mc_RES->SetLineColor(kGray+3);
   mc_RES->SetLineWidth(1.9);
//mc_RES->SetFillStyle(3002);
   mc_RES->SetFillColor(kYellow-9);
  
   mc_RES_line->SetLineColor(kGray+3);
   mc_RES_line->SetLineWidth(1.95);

   mc_DIS->SetLineColor(kGray+3);
   mc_DIS->SetLineWidth(1.9);
//mc_DIS->SetFillStyle(3002);
   mc_DIS->SetFillColor(kAzure+10);
   
mc_DIS_line->SetLineColor(kGray+3);
   mc_DIS_line->SetLineWidth(1.95);

   mc_2p2h->SetLineColor(kGray+3);
   mc_2p2h->SetLineWidth(1.9);
//mc_2p2ph->SetFillStyle(3002);
   mc_2p2h->SetFillColor(kMagenta-9);
   mc_QE_line->Add(mc_2p2h);
  
  //mc_QE_line->Add(2p2h);
 mc_DIS_line->Add(mc_other);
 mc_DIS->Add(mc_other);
 mc_2p2h->Add(mc_DIS); 
 mc_QE->Add(mc_2p2h);
  mc_RES->Add(mc_QE);
 mc_RES_line->Reset();
 mc_RES_line->Add(mc_RES); 

   mc_bkg->SetLineColor(kGray+3);
   mc_bkg->SetLineWidth(1.9);
//mc_bkg->SetFillStyle(3002);
   mc_bkg->SetFillColor(kYellow-7);


  signal->SetLineColor(kBlack);
  signal->SetLineWidth(3);
  signal->SetFillColor(kBlue-9); 

  data->SetMarkerStyle(kFullCircle);
  data->SetMarkerSize(0.9); 
  data->SetMarkerColor(kPink+5);
  data->SetLineColor(kBlack);
  data->SetLineWidth(2); 

  THStack* hs = new THStack("hs","Stacked histograms");
    hs->Add(mc);
    hs->Add(signal);
    hs->Add(mc_QE);
    hs->Add(mc_RES);
    hs->Add(mc_DIS);
    hs->Add(mc_2p2h);
   // hs->Add(h_other);


  }
else {
  mc->SetLineColor(kBlack);
  mcTotalError->SetFillColor(kBlue-10);
  mc->SetLineWidth(2);
  mc->SetFillColor(kAzure+6);
  signal->SetLineColor(kBlack);
  signal->SetLineWidth(3);
  signal->SetFillColor(kBlue-9);
  
    data->SetMarkerStyle(kFullCircle);
  data->SetMarkerSize(0.92);
  data->SetMarkerColor(kPink+5);
  data->SetLineColor(kBlack);
  data->SetLineWidth(2.1);
}
  //Add 2p2h with qe to reduce number of categories
//  mc_qe->Add(mc_2p2h);

 
  if(doRatio)
  {
   data->Divide(mc);
   dataStat->Divide(mc);   
   mc_line->Divide(mc_line);
   mcTotalError->Divide(mc);
   mc->Divide(mc);
  //  mcTotalError->Divide(mc);
   //mc_line->Divide(mc);
   signal->Divide(signal);
  // mcTotalError->Divide(mc);
  // dataStat->Divide(mc);
  
  }
    

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
  std::vector<std::pair<THStack*, const char*> > histAndOptsStack;
if(!doRatio){
  histAndOpts.push_back(std::make_pair(mcTotalError,       "E2"));	
//histAndOptsStack.push_back(std::make_pair(hs,       "hist"));
 

//  histAndOpts.push_back(std::make_pair(mc,       "hist"));
//  histAndOpts.push_back(std::make_pair(signal,       "hist"));
  histAndOpts.push_back(std::make_pair(mc_RES,       "hist"));
  histAndOpts.push_back(std::make_pair(mc_QE,       "hist"));
  histAndOpts.push_back(std::make_pair(mc_2p2h,       "hist"));
  histAndOpts.push_back(std::make_pair(mc_DIS,       "hist"));
  histAndOpts.push_back(std::make_pair(mc_other,       "hist"));
  histAndOpts.push_back(std::make_pair(data,  "X0"     "histpe"));
 histAndOpts.push_back(std::make_pair(mcTotalError,       "E2"));
   histAndOpts.push_back(std::make_pair(mc_line,       "hist")); 
    histAndOpts.push_back(std::make_pair(mc_RES_line,       "hist"));
histAndOpts.push_back(std::make_pair(mc_DIS_line,       "hist"));
    histAndOpts.push_back(std::make_pair(data,  "X0"     "histpe1"));  

   }
else{
histAndOpts.push_back(std::make_pair(mc,       "hist"));
//histAndOpts.push_back(std::make_pair(data,     "histpe1"));
// histAndOpts.push_back(std::make_pair(dataStat, "X0" "histpe1"));
histAndOpts.push_back(std::make_pair(mcTotalError,       "E2"));
 histAndOpts.push_back(std::make_pair(mc_line,       "hist"));
histAndOpts.push_back(std::make_pair(data, "X0"    "histpe1"));
}

  // ------------------------------------------------------
  // ----------------------------
  //
  // First make pt in bins of pz

  // Values to multiply each bin by to get them on a similar range
//  vector<double> multipliers = GetScales(histAndOpts, true, 5,0.75);
 // double  multipliersxQ21[]={1,3, 20, 30,1,
   //                  8, 1, 1, 1, 1, 1, 1, 1};
//  double  multipliersxQ21[]={1, 5, 15, 40, 2,
//                      4, 1, 1, 1, 1, 1, 1, 1};
//  double  multipliersxQ21[]={1,5, 20, 80, 700,50000};
double  multipliersxQ21[9]={1, 5, 15,65,500,15000};
double  multipliersWQ21[]={2,3,10,50,700,25000};
double multiplierspZpT[]={15,4,2,2,2,2,4,25,800};
double multiplierspTpZ[]={20,12,10,10,12,25,150,800};
 std::string var1,var2,var3,var4;
 int ncol,nrow; 

  if (varName=="x_Q2")
{ var1 = "Bjorken x";
   var2 = "Q^{2}(GeV^{2})"; 
   var3 = "x";
   var4 = "Q^{2}";
  ncol=3; nrow=2;}
  else if(varName=="W_Q2") 
{  var1 = "W (GeV)";
   var2 =  "Q^{2}(GeV^{2})";
   var3 = "W";
   var4 = "Q^{2}";
   size_t size = sizeof(multipliersWQ21) / sizeof(multipliersWQ21[0]);
   for (size_t i = 0; i < size; ++i) {
   multipliersxQ21[i]=multipliersWQ21[i];}
   ncol=3; nrow=2;}
 else if(varName=="pZmu_pTmu") 
{  var1 = "p_{z}(GeV)";
   var2 = "p_{t}(GeV)"; 
   var3 ="p_{z}";
   var4 = "p_{t}";
  size_t size = sizeof(multiplierspZpT) / sizeof(multiplierspZpT[0]);
   for (size_t i = 0; i < size; ++i) {
   multipliersxQ21[i]=multiplierspZpT[i];
   }
   ncol=3; nrow=3; }
   GridCanvas* gc=NULL;     

if(varName=="pZmu_pTmu" && ispT==true)
gc=plotYAxis1D( histAndOpts, var2, var1,ncol,nrow,850,520, doMultipliers ? multiplierspTpZ :NULL);
else
  gc=plotXAxis1D(histAndOpts, var1, var2,ncol,nrow,850,520, doMultipliers ? multipliersxQ21 :NULL);


   //if(varName=="pZmu_pTmu")
//  gc=plotYAxis1D( histAndOpts, var2, var1,ncol,nrow,850,520, doMultipliers ? multiplierspTpZ :NULL); 

 //GridCanvas* gc=plotXAxis1D(histAndOpts, "p_{z}(GeV)", "p_t(GeV)",3,2,850,520, doMultipliers ? multipliersxQ21 :NULL);
  // Set the y range manually. Can also use gc->Remax() to guess automatically
   gc->Remax();
 //  gc->SetXLimits(0,9.99);
//  gc->SetYLimits(0, 6000);
//  if(doRatio) gc->SetYLimits(0,1.99);
//  else  gc->SetYLimits(0, 1.99);
  if(doRatio){ gc->SetYTitle("Data/MC ratio");
  //   if(varName=="W_Q2")
	  gc->SetYLimits(0,1.99);
  }
// else gc->SetYTitle("Events  per (GeV)^{2}");
   // else gc->SetYTitle("Event Rate x 10^{-4}");
else{   if(varName!="W_Q2")
            gc->SetYTitle(Form("d^{2}N/d%sd%s  x(10^{-4}) GeV^{-2}",var3.c_str(),var4.c_str()));
            else
             gc->SetYTitle(Form("d^{2}N/d%sd%s  x(10^{-4}) GeV^{-3}",var3.c_str(),var4.c_str()));
                     }

  gc->Modified();
  // Example of adding a legend. The co-ordinate system is NDC on the
  // entire canvas, ie (0,0) in the bottom left corner of the canvas
  // (not the individual pad), and (1,1) in the top right
  TLegend* leg=new TLegend(0.0, 0.0, 0.4, 0.1);
TLegend* leg0=new TLegend(0.68, 0.0, 0.95, 0.1);
 TLatex* latex = new TLatex(0, 0, "#color[2]{#font[72]{MINERvA Work in Progress}}");
    leg0->AddEntry(latex, "", "");
    leg0->SetBorderSize(0);
    leg0->SetTextSize(0.035);  

  leg->SetNColumns(3);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.026);
  leg->SetTextFont(62);
 // leg->AddEntry(data, "Data", "lpe1");
if(!doRatio){ }// leg->AddEntry(mc, "MINERvA tuneV4", "f");}
else{ leg->SetTextSize(0.035);
//    leg->AddEntry(data, "MINERvA data ", "lpe1");
    leg->AddEntry(mc, "Total MC","f");
    leg->AddEntry(data, "MINERvA data ", "lpe1");
}
if(!doRatio){
//	leg->AddEntry(mc, "Total MC", "f");
 // leg->AddEntry(signal,"Signal","f");
  leg->AddEntry(mc_2p2h,"2p2h","f");
  leg->AddEntry(mc_QE,"QuasiElastic","f");
  leg->AddEntry(mc_RES, "Resonance","f");
  leg->AddEntry(mc_DIS, "DIS","f");
//  leg->AddEntry(mc_2p2h, "2p2h","f");
  leg->AddEntry(mc_other, "Other","f");
  leg->AddEntry(data, "MINERvA Data","lpe1");
  
 // leg->AddEntry(data_bkg,"Data constrained bkg prediction","lpe1");
}

  TLegend* leg2 = new TLegend(0.7, 0.0, 1, 0.1);
  leg2->SetNColumns(2);
  leg2->SetFillStyle(0);
  leg2->SetBorderSize(0);
  leg2->SetTextSize(0.025);
 

// leg2->AddEntry(mc_Trans,"Trans","l");
// leg2->AddEntry(mc_Contin, "Contin", "l");
// leg2->AddEntry(mc_US,"US","l");
 //leg2->AddEntry(mc_DS,"DS","l");


//leg2->AddEntry(mc_dis_sis,"Wrong sign","l");
//leg2->AddEntry(mc_2p2h,"Other","l");
  TLegend* leg3=new TLegend(0.13, 0.81, 0.4, 0.91);
  leg3->SetNColumns(1);
  leg3->SetFillStyle(0);
  leg3->SetBorderSize(0);
  leg3->SetTextSize(0.02);
//  leg3->AddEntry((TObject*)0, "MINERvA Preliminary", "");
  leg3->SetTextColor(4);
//  leg3->AddEntry((TObject*)0,"DataPOT:1.61e+20","");
//  leg3->AddEntry((TObject*)0,"MCPOT:6.4e+20", "");
   leg3->AddEntry((TObject*)0,"DataPOT:1.118e+21","");
  leg3->AddEntry((TObject*)0,"MCPOT:4.962e+21", "");

  std::string material="Goofy";
  if(targetZ==6)
     material="Carbon";
  else if(targetZ==26)
	material="Iron";
  else if(targetZ==82)
	material="Lead";
//TLatex* latex1 = new TLatex(0, 0, "#color[2]{#font[72]{MINERvA Work in Progress}}");
//title->AddEntry(latex1, "", "");
  TLegend* title = new TLegend(0.29, 0.95, 0.95, 1);
  if(targetID==99){
   TLatex* latex2 = new TLatex(0, 0, Form("#color[1]{#font[72]{Tracker}}"));
   title->AddEntry(latex2, "", "");	
//	  title->SetHeader(Form("Tracker","C")); // option "C" allows to center the header
  }
   else {
    if(IsTargetCombined==false)	   
    { TLatex* latex1 = new TLatex(0, 0, Form("#color[1]{#font[72]{Target %d %s}}",targetID,material.c_str()));
   // title->SetHeader(Form("Target %d %s", targetID, material.c_str(),"C"));
      title->AddEntry(latex1, "", "");}
    else
    {
	    if(targetZ==26){
		    TLatex* latex1 = new TLatex(0, 0, "#color[1]{#font[72]{Combined Iron (t2+t3+t5)}}");
   title->AddEntry(latex1, "", "");}
	    else if(targetZ==82){
		    TLatex* latex1 = new TLatex(0, 0, "#color[1]{#font[72]{Combined Lead (t2+t3+t4+t5)}}");
		    title->AddEntry(latex1, "", "");}
   //title->SetHeader(Form("Combined %s", material.c_str(),"C"));	    
   }
   };
  title->SetBorderSize(0);
  title->SetFillStyle(0);
  title->SetTextSize(0.04);
  title-> SetTextFont(72);
  title->Draw();



  
  TLatex *mytex = new TLatex();
  TLine *myline = new TLine();
  mytex->SetTextFont(42);
  mytex->SetTextSize(0.050);
  leg->Draw("SAME");
//  leg0->Draw("SAME");//prints MINERvA work in Progress
  leg2->Draw("SAME");
  // leg3->Draw("SAME");
  //  mytex->DrawLatex(0.67,0.23,"GENIE Components:");
  //  myline->DrawLine(0.67,0.225,0.847,0.225);
    
  if(doRatio){
   if(IsTargetCombined==false){	 
	  if(ispT==true)
		  gc->Print(doMultipliers ? Form("Event_ptpz_reverse_%d_%02d_channel_breakdown_ratio.pdf",targetID,targetZ) : Form("Event_ptpz_reverse_%d_%02d_mult_channel_breakdown_ratio.pdf",targetID,targetZ));
	  else 
    gc->Print(doMultipliers ? Form("Event_%s_%d_%02d_channel_breakdown_ratio.pdf",varName.c_str(),targetID,targetZ) : Form("Event_%s_%d_%02d_mult_channel_breakdown_ratio.pdf",varName.c_str(),targetID,targetZ));
   } else{
	    if(ispT==true)
		   gc->Print(doMultipliers ? Form("Event_ptpz_reverse_%d_%02d_channel_breakdown_ratio_TargetCombined.pdf",targetID,targetZ) : Form("Event_ptpz_reverse_%d_%02d_mult_channel_breakdown_ratio_TargetCombined.pdf",targetID,targetZ));
	    else 
	   gc->Print(doMultipliers ? Form("Event_%s_%d_%02d_channel_breakdown_ratio_TargetCombined.pdf",varName.c_str(),targetID,targetZ) : Form("Event_%s_%d_%02d_mult_channel_breakdown_ratio_TargetCombined.pdf",varName.c_str(),targetID,targetZ));
   // gc->Print(doMultipliers ? Form("Event_%s_%d_%02d_channel_breakdown_mul_ratio.png",varName.c_str(),targetID,targetZ) : Form("Event_%s_%d_%02d_channel_breakdown_ratio.png",varName.c_str(),targetID,targetZ));
//    gc->Print(doMultipliers ? "nu-2d-xsec-comps-pt-multiplier_ratio.C" : "nu-2d-xsec-comps-pt_ratio.C");
   // gc->SetYLimits(0,2);

 }
  }

  else{
   // gc->Print(doMultipliers ? Form("Event_%s_%d_%02d_mult_channel_breakdown.pdf",varName.c_str(),targetID,targetZ): Form("Event_%s_%d_%02d_channel_breakdown.pdf",varName.c_str(),targetID,targetZ));
   // gc->Print(doMultipliers ? Form("Event_%s_%d_%02d_channel_breakdown_mul.png",varName.c_str(),targetID,targetZ) : Form("Event_%s_%d_%02d_channel_breakdown.png",varName.c_str(),targetID,targetZ));
   // gc->Print(doMultipliers ? "Event_x_Q2_target_multiplier_bkg_Sub.C" : "Event_x_Q2_target_bkg_Sub.C");
 if(IsTargetCombined==true){
	 if(ispT==true)
   gc->Print(doMultipliers ? Form("Event_ptpz_reverse_%d_%02d_mult_channel_breakdown_TargetCombined.pdf",targetID,targetZ): Form("Event_%s_%d_%02d_channel_breakdown_TargetCombined.pdf",varName.c_str(),targetID,targetZ));
	 else
  gc->Print(doMultipliers ? Form("Event_%s_%d_%02d_mult_channel_breakdown_TargetCombined.pdf",varName.c_str(),targetID,targetZ): Form("Event_%s_%d_%02d_channel_breakdown_TargetCombined.pdf",varName.c_str(),targetID,targetZ));
 }
 	 
 else{
	   if(ispT==true)
		    gc->Print(doMultipliers ? Form("Event_ptpz_reverse_%d_%02d_mult_channel_breakdown.pdf",targetID,targetZ): Form("Event_%s_%d_%02d_channel_breakdown.pdf",varName.c_str(),targetID,targetZ));
	   else
 gc->Print(doMultipliers ? Form("Event_%s_%d_%02d_mult_channel_breakdown.pdf",varName.c_str(),targetID,targetZ): Form("Event_%s_%d_%02d_channel_breakdown.pdf",varName.c_str(),targetID,targetZ));
 }

  
}

  // ----------------------------------------------------------------------------------
  //
  // Now make pz in bins of pt. It's all the same

  // Values to multiply each bin by to get them on a similar range
//  vector<double> multipliers2 = GetScales(histAndOpts, false, 2.5,0.75);
 
 /* double  multipliersxQ22[]={3, 1, 1, 2, 10,
                           100, 1, 1, 1, 1, 1, 1, 1};

  

  // plotpz1D fiddles the x axis values to squash up the tail so it
  // doesn't take up all the horizontal space.
  //GridCanvas* gc2=plotYAxis1D(histAndOpts, "Muon Transverse Momentum (GeV/c)", "p_{||}",4,4,800,500,doMultipliers ? &multipliers2[0] : NULL);
  GridCanvas* gc2=plotYAxis1D(histAndOpts, "Q^2", "x",4,4,800,500,doMultipliers ? multipliersxQ22 : NULL);
	gc->Remax();
//  if(doRatio) gc2->SetYLimits(0,1.99);
//  else  gc2->SetYLimits(0, 2.49);
  if(doRatio) gc2->SetYTitle("Ratio data/MINERvA Tune v1");
  else gc2->SetYTitle("Events  per (GeV/c)^{2}");
  gc2->Modified();
  if(doRatio){
    gc2->Print(doMultipliers ? "nu-2d-xsec-comps-pz-multiplier_ratio.eps" : "nu-2d-xsec-comps-pz_ratio.eps");
    gc2->Print(doMultipliers ? "nu-2d-xsec-comps-pz-multiplier_ratio.png" : "nu-2d-xsec-comps-pz_ratio.png");
//    gc2->Print(doMultipliers ? "nu-2d-xsec-comps-pz-multiplier_ratio.C" : "nu-2d-xsec-comps-pz_ratio.C");
  }
  else{
    gc2->Print(doMultipliers ? "nu-2d-xsec-comps-pz-multiplier.eps" : "nu-2d-xsec-comps-pz.eps");
    gc2->Print(doMultipliers ? "nu-2d-xsec-comps-pz-multiplier.png" : "nu-2d-xsec-comps-pz.png");
//    gc2->Print(doMultipliers ? "nu-2d-xsec-comps-pz-multiplier.C" : "nu-2d-xsec-comps-pz.C");
  }*/

}

int main(int argc, char* argv[])
{
	std::string varName = argv[2];
	int targetID = std::atoi(argv[3]);
	int targetZ = std::atoi(argv[4]);
	//bool ispT= argv[5];
//	bool IsTargetCombined = argv[6];

	 bool ispT = std::string(argv[5]) == "true" || std::string(argv[5]) == "1";
    bool IsTargetCombined = std::string(argv[6]) == "true" || std::string(argv[6]) == "1";
  //  makePlots(true,true,argv[1]); //multipliers on ratio don't make sense
  makePlots(true,false,argv[1],argv[2],targetID,targetZ,ispT,IsTargetCombined);
  makePlots(false,true,argv[1],argv[2],targetID, targetZ,ispT,IsTargetCombined);
 // makePlots(false,false,argv[1],argv[2],targetID, targetZ,IsTargetCombined);
  return 0;
}
