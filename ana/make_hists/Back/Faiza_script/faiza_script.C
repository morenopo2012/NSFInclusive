#include "TFile.h"
#include "TFrame.h"
#include "TArrow.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TH1F.h"
#include "MnvH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TGaxis.h"
#include "TH1.h"
#include "TMath.h"
#include <stdio.h>
#include "TStyle.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TRandom.h"
#include "TColor.h"
#include <vector>

#include "PlotUtils/MacroUtil.h" 
#include "PlotUtils/ChainWrapper.h"
#include "PlotUtils/makeChainWrapper.h"
#include "PlotUtils/HistWrapper.h"
#include "PlotUtils/Hist2DWrapper.h"
#include "PlotUtils/GenieSystematics.h"
#include "PlotUtils/FluxSystematics.h"
#include "PlotUtils/MnvTuneSystematics.h"
#include "PlotUtils/MinosEfficiencySystematics.h"
#include "PlotUtils/MnvHadronReweight.h" 
#include "PlotUtils/FluxReweighter.h"
#include "PlotUtils/MinosMuonEfficiencyCorrection.h"
#include "PlotUtils/AngleSystematics.h"
#include "PlotUtils/MuonSystematics.h"
#include "PlotUtils/ResponseSystematics.h"
#include "PlotUtils/MuonResolutionSystematics.h"

using namespace std;

void faiza_script(){
   gStyle->SetOptStat(0);
   gStyle->SetLineWidth(3);
   gStyle->SetFrameLineWidth(3);
   gStyle->SetLegendFont(22);
   gStyle->SetLegendBorderSize(0);
   TGaxis::SetMaxDigits(6);

//   TFile *caffe   = new TFile("/minerva/data/users/afilkins/NukeHists/faizaML_v1_CaffeMENNDL/minervame6A/Hists_PlaneStacks_t2_z26_AntiNu_faizaML_v1_CaffeMENNDL.root");
   //TFile *menndl  = new TFile("/minerva/data/users/afilkins/NukeHists/faizaML_v1_CaffeArtisinal/minervame6A/Hists_PlaneStacks_t2_z26_AntiNu_faizaML_v1_CaffeArtisinal.root");
   TFile *menndl = new TFile("/minerva/data/users/fakbar/NukeHists/v1_CaffeMENNDL/minervame6A/Hists_PlaneStacksDIS_t2_z26_AntiNu_v1_CaffeMENNDL.root");
   TFile *caffe  = new TFile("/minerva/data/users/fakbar/NukeHists/v1_CaffeArtisinal/minervame6A/Hists_PlaneStacksDIS_t2_z26_AntiNu_v1_CaffeArtisinal.root");
   double potMC_menndl   = 3.184474e+20;
   double potData_menndl = 1.622756e+20; 
   double potData_caffe  = 1.622756e+20;
   double potMC_caffe    = 3.191510e+20; 
   double binwidth       = 1.7;
   double scale_menndl   = potData_menndl/(potMC_menndl*binwidth);
   double scale_caffe    = potData_caffe/(potMC_caffe*binwidth);
   cout<<scale_menndl<<"  "<<scale_caffe<<endl;

   caffe->cd();
   MnvH1D* caffe_Fe      =  (MnvH1D*)caffe->Get("hFe_vtxzMC_t2_z26");
   MnvH1D* caffe_Pb      =  (MnvH1D*)caffe->Get("hPb_vtxzMC_t2_z26");
   MnvH1D* caffe_C       =  (MnvH1D*)caffe->Get("hC_vtxzMC_t2_z26");
   MnvH1D* caffe_US      =  (MnvH1D*)caffe->Get("hUS_vtxzMC_t2_z26");
   MnvH1D* caffe_DS      =  (MnvH1D*)caffe->Get("hDS_vtxzMC_t2_z26");
   MnvH1D* caffe_Other   =  (MnvH1D*)caffe->Get("hOther_vtxzMC_t2_z26");

   menndl->cd();
   MnvH1D* menndl_Fe     =  (MnvH1D*)menndl->Get("hFe_vtxzMC_t2_z26");
   MnvH1D* menndl_Pb     =  (MnvH1D*)menndl->Get("hPb_vtxzMC_t2_z26");
   MnvH1D* menndl_C      =  (MnvH1D*)menndl->Get("hC_vtxzMC_t2_z26");
   MnvH1D* menndl_US     =  (MnvH1D*)menndl->Get("hUS_vtxzMC_t2_z26");
   MnvH1D* menndl_DS     =  (MnvH1D*)menndl->Get("hDS_vtxzMC_t2_z26");
   MnvH1D* menndl_Other  =  (MnvH1D*)menndl->Get("hOther_vtxzMC_t2_z26");

   caffe_Fe->Scale(scale_caffe);
   caffe_Pb->Scale(scale_caffe);
   caffe_C->Scale(scale_caffe);
   caffe_US->Scale(scale_caffe);
   caffe_DS->Scale(scale_caffe);
   caffe_Other->Scale(scale_caffe);

   menndl_Fe->Scale(scale_menndl);
   menndl_Pb->Scale(scale_menndl);
   menndl_C->Scale(scale_menndl);
   menndl_US->Scale(scale_menndl);
   menndl_DS->Scale(scale_menndl);
   menndl_Other->Scale(scale_menndl);

   MnvH1D* ratio_Fe      = (MnvH1D*)menndl_Fe->Clone("ratio_");
   MnvH1D* ratio_Pb      = (MnvH1D*)menndl_Pb->Clone("ratio_");
   MnvH1D* ratio_C       = (MnvH1D*)menndl_C ->Clone("ratio_");
   MnvH1D* ratio_US      = (MnvH1D*)menndl_US->Clone("ratio_");
   MnvH1D* ratio_DS      = (MnvH1D*)menndl_DS->Clone("ratio_");
   MnvH1D* ratio_Other   = (MnvH1D*)menndl_Other->Clone("ratio_Other");

   ratio_Fe              -> Divide(caffe_Fe,menndl_Fe,1.0,1.0,"B");
   ratio_Pb              -> Divide(caffe_Pb,menndl_Pb,1.0,1.0,"B");
   ratio_C               -> Divide(caffe_C ,menndl_C,1.0,1.0,"B");
   ratio_US              -> Divide(caffe_US,menndl_US,1.0,1.0,"B");
   ratio_DS              -> Divide(caffe_DS,menndl_DS,1.0,1.0,"B");
   ratio_Other           -> Divide(caffe_Other,menndl_Other,1.0,1.0,"B");


   TCanvas *c1 = new TCanvas("c1","A Simple Graph Example",200,10,1200,800);
    caffe_Other->Draw("hist");
    caffe_US->Draw("histsame");
    caffe_DS->Draw("histsame");
    caffe_Fe->Draw("histsame");

   TCanvas *c2 = new TCanvas("c2","A Simple Graph Example",200,10,1200,800);
    menndl_Other->Draw("hist");
    menndl_US->Draw("histsame");
    menndl_DS->Draw("histsame");
    menndl_Fe->Draw("histsame");

   TCanvas *c3 = new TCanvas("c3","A Simple Graph Example",200,10,1200,800);
    ratio_Other->Draw("l");
    ratio_US->Draw("lsame");
    ratio_DS->Draw("lsame");
    ratio_Fe->Draw("lsame");


    //
    /*ratio_Fe->Draw("c");
    ratio_Pb->Draw("lpcsame");
    ratio_C->Draw("lpcsame");
    ratio_US->Draw("lpcsame");
    ratio_DS->Draw("lpcsame");
    ratio_Other->Draw("lpcsame");
*/
    ratio_Other->SetLineWidth(4);
    ratio_Fe->SetLineWidth(4);
    ratio_US->SetLineWidth(4);
    ratio_DS->SetLineWidth(4);
    ratio_Other->GetXaxis()->SetRangeUser(455,485);
    ratio_Other->GetYaxis()->SetRangeUser(0.0, 2.0);
    caffe_Other->GetXaxis()->SetRangeUser(455,485);
    caffe_Other->GetYaxis()->SetRangeUser(0,800);
    menndl_Other->GetXaxis()->SetRangeUser(455,485);
    menndl_Other->GetYaxis()->SetRangeUser(0,800);
    caffe_Other->SetTitle("DED-Caffe");
    menndl_Other->SetTitle("MENNDL");
    ratio_Other->SetTitle("DED-Caffe/MENNDL");
    caffe_Other->GetXaxis()->SetTitle("Vertex Z (cm)");
    menndl_Other->GetXaxis()->SetTitle("Vertex Z (cm)");
    ratio_Other->GetXaxis()->SetTitle("Vertex Z (cm)");

    auto legend = new TLegend(0.5,0.65,0.85,0.85);
    legend->AddEntry(ratio_Fe,"Fe","l"); 
    legend->AddEntry(ratio_US,"US","l"); 
    legend->AddEntry(ratio_DS,"DS","l"); 
    legend->AddEntry(ratio_Other,"Other","l"); 
    legend->Draw();
}

