#include "PlotUtils/MnvH1D.h"
#include "TFile.h"
#include "TArrow.h"
#include "TLegend.h"
#include "TCanvas.h"
using namespace PlotUtils;

void plot_phy(){ 

  TString hname1 = Form("/minerva/data/users/fakbar/NukeHists/v1_MAD/Iron_combined_phy.root");
  TString hname1 = Form("/minerva/data/users/fakbar/NukeHists/v1_MAD/Hists_EventSelection_phys_without_sys_t2_z26_Nu_v1_MAD_.root");
  cout<<hname1<<endl;
  TFile *f1 = new TFile( hname1,"read" );

  TString hname2 = Form("/minerva/data/users/fakbar/NukeHists/v1_nukeCC_validation/minervame1A/Hists_SideBands_t2_z26_Nu_v1_nukeCC_validation.root");
  cout<<hname2<<endl;
  TFile *f2 = new TFile( hname2,"read" );

   f1->cd();
   MnvH1D *Trans_data       = (MnvH1D*)f1->Get("selected_data_reco_trans_Emu");
   MnvH1D *Trans_mc_TinT    = (MnvH1D*)f1->Get("hists_trans_in_trans_Emu");
   MnvH1D *Trans_mc_CinT    = (MnvH1D*)f1->Get("hists_contin_in_trans_Emu");
   MnvH1D *Trans_mc_SinT    = (MnvH1D*)f1->Get("hists_signal_in_trans_Emu");

   MnvH1D *Contin_data      = (MnvH1D*)f1->Get("selected_data_reco_contin_Emu");
   MnvH1D *Contin_mc_TinC   = (MnvH1D*)f1->Get("hists_trans_in_contin_Emu");
   MnvH1D *Contin_mc_CinC   = (MnvH1D*)f1->Get("hists_contin_in_contin_Emu");
   MnvH1D *Contin_mc_SinC   = (MnvH1D*)f1->Get("hists_signal_in_contin_Emu");

   MnvH1D *Trans_mc_CH_US   = (MnvH1D*)f1->Get("hists_trans_mc_CH_US_Emu");
   MnvH1D *Trans_mc_CH_DS   = (MnvH1D*)f1->Get("hists_trans_mc_CH_DS_Emu");
   MnvH1D *Contin_mc_CH_US  = (MnvH1D*)f1->Get("hists_contin_mc_CH_US_Emu");
   MnvH1D *Contin_mc_CH_DS  = (MnvH1D*)f1->Get("hists_contin_mc_CH_DS_Emu");

   f2->cd();
   MnvH1D *Trans_data_nuke       = (MnvH1D*)f2->Get("Data_trans_Emu_t2_z26");
   MnvH1D *Trans_mc_TinT_nuke    = (MnvH1D*)f2->Get("MC_true_trans_in_trans_Emu_t2_z26");
   MnvH1D *Trans_mc_CinT_nuke    = (MnvH1D*)f2->Get("MC_true_contin_in_trans_Emu_t2_z26");
   MnvH1D *Trans_mc_SinT_nuke    = (MnvH1D*)f2->Get("MC_true_signal_in_trans_Emu_t2_z26");

   MnvH1D *Contin_data_nuke      = (MnvH1D*)f2->Get("Data_contin_Emu_t2_z26");
   MnvH1D *Contin_mc_TinC_nuke   = (MnvH1D*)f2->Get("MC_true_trans_in_contin_Emu_t2_z26");
   MnvH1D *Contin_mc_CinC_nuke   = (MnvH1D*)f2->Get("MC_true_contin_in_contin_Emu_t2_z26");
   MnvH1D *Contin_mc_SinC_nuke   = (MnvH1D*)f2->Get("MC_true_signal_in_contin_Emu_t2_z26");

   MnvH1D *Trans_mc_CH_US_nuke   = (MnvH1D*)f2->Get("MC_trans_CH_US_Emu_t2_z26");
   MnvH1D *Trans_mc_CH_DS_nuke   = (MnvH1D*)f2->Get("MC_trans_CH_DS_Emu_t2_z26");
   MnvH1D *Contin_mc_CH_US_nuke  = (MnvH1D*)f2->Get("MC_contin_CH_US_Emu_t2_z26");
   MnvH1D *Contin_mc_CH_DS_nuke  = (MnvH1D*)f2->Get("MC_contin_CH_DS_Emu_t2_z26");
//*******************************************************
//			TRANSITION
//*******************************************************
{
   Trans_data->SetMarkerColor(kBlack);
   Trans_mc_TinT->SetMarkerColor(kOrange+2);
   Trans_mc_CinT->SetMarkerColor(kBlue);
   Trans_mc_SinT->SetMarkerColor(kGreen+2);

   Trans_data->SetLineColor(kBlack);
   Trans_mc_TinT->SetLineColor(kOrange+2);
   Trans_mc_CinT->SetLineColor(kBlue);
   Trans_mc_SinT->SetLineColor(kGreen+2);

   Trans_data_nuke->SetMarkerColor(kBlack);
   Trans_mc_TinT_nuke->SetMarkerColor(kOrange+2);
   Trans_mc_CinT_nuke->SetMarkerColor(kBlue);
   Trans_mc_SinT_nuke->SetMarkerColor(kGreen+2);

   Trans_data_nuke->SetLineColor(kBlack);
   Trans_mc_TinT_nuke->SetLineColor(kOrange+2);
   Trans_mc_CinT_nuke->SetLineColor(kBlue);
   Trans_mc_SinT_nuke->SetLineColor(kGreen+2);

//Ratio = MAD/NukeCC
  MnvH1D* Trans_data_histo_nuke     = (MnvH1D*)Trans_data_nuke->Clone("trans_histo_data");
  MnvH1D* Trans_mc_TinT_histo_nuke  = (MnvH1D*)Trans_mc_TinT_nuke->Clone("trans_histo_mc_tint");
  MnvH1D* Trans_mc_CinT_histo_nuke  = (MnvH1D*)Trans_mc_CinT_nuke->Clone("trans_histo_mc_cint");
  MnvH1D* Trans_mc_SinT_histo_nuke  = (MnvH1D*)Trans_mc_SinT_nuke->Clone("trans_histo_mc_sint");

   Trans_data_histo_nuke->Divide(Trans_data, Trans_data_nuke, 1,1,"B");
   Trans_mc_TinT_histo_nuke->Divide(Trans_mc_TinT, Trans_mc_TinT_nuke, 1,1,"B");
   Trans_mc_CinT_histo_nuke->Divide(Trans_mc_CinT, Trans_mc_CinT_nuke, 1,1,"B");
   Trans_mc_SinT_histo_nuke->Divide(Trans_mc_SinT, Trans_mc_SinT_nuke, 1,1,"B");

   TCanvas *c1 = new TCanvas("c1","Trans Data Ratio:NSF/OSF",600,500);
   c1->SetGrid(1);
   Trans_data_histo_nuke->Draw();
//   TCanvas *c2 = new TCanvas("c2","Trans in Trans Ratio:NSF/OSF",600,500);
//   c2->SetGrid(1);
   Trans_mc_TinT_histo_nuke->Draw("Psame][");
//   TCanvas *c3 = new TCanvas("c3","Contin in Trans Ratio:NSF/OSF",600,500);
//   c3->SetGrid(1);
   Trans_mc_CinT_histo_nuke->Draw("Psame][");
//   TCanvas *c4 = new TCanvas("c4","Signal in Trans Ratio:NSF/OSF",600,500);
//   c4->SetGrid(1);
   Trans_mc_SinT_histo_nuke->Draw("Psame][");

   Trans_data_histo_nuke   ->GetYaxis()->SetRangeUser(0.5,1.5);
   Trans_mc_TinT_histo_nuke->GetYaxis()->SetRangeUser(0.5,1.5);
   Trans_mc_CinT_histo_nuke->GetYaxis()->SetRangeUser(0.5,1.5);
   Trans_mc_SinT_histo_nuke->GetYaxis()->SetRangeUser(0.5,1.5);

         TLegend* leg = new TLegend(0.5,0.6,0.65,0.85);
         leg->AddEntry(Trans_data_histo_nuke,"Data.", "lep");
         leg->AddEntry(Trans_mc_TinT_histo_nuke,"Trans.", "lep");
         leg->AddEntry(Trans_mc_CinT_histo_nuke,"Contin.", "lep");
         leg->AddEntry(Trans_mc_SinT_histo_nuke,"Signal.", "lep");
         leg->Draw();
   gPad->SetGrid(1);
   gPad->Print("FluxValidation/Ratio_NSF_OSF_Trans_Emu_Iron.pdf","pdf");
   gPad->Print("FluxValidation/Ratio_NSF_OSF_Trans_Emu_Iron.png","png");
   gPad->Update();

}
//*******************************************************
//			CONTINUUM
//*******************************************************
{
   Contin_data->SetMarkerColor(kBlack);
   Contin_mc_TinC->SetMarkerColor(kOrange+2);
   Contin_mc_CinC->SetMarkerColor(kBlue);
   Contin_mc_SinC->SetMarkerColor(kGreen+2);

   Contin_data->SetLineColor(kBlack);
   Contin_mc_TinC->SetLineColor(kOrange+2);
   Contin_mc_CinC->SetLineColor(kBlue);
   Contin_mc_SinC->SetLineColor(kGreen+2);

   Contin_data_nuke->SetMarkerColor(kBlack);
   Contin_mc_TinC_nuke->SetMarkerColor(kOrange+2);
   Contin_mc_CinC_nuke->SetMarkerColor(kBlue);
   Contin_mc_SinC_nuke->SetMarkerColor(kGreen+2);

   Contin_data_nuke->SetLineColor(kBlack);
   Contin_mc_TinC_nuke->SetLineColor(kOrange+2);
   Contin_mc_CinC_nuke->SetLineColor(kBlue);
   Contin_mc_SinC_nuke->SetLineColor(kGreen+2);

//Ratio = MAD/NukeCC
  MnvH1D* Contin_data_histo_nuke     = (MnvH1D*)Contin_data_nuke->Clone("contin_histo_data");
  MnvH1D* Contin_mc_TinC_histo_nuke  = (MnvH1D*)Contin_mc_TinC_nuke->Clone("contin_histo_mc_tinc");
  MnvH1D* Contin_mc_CinC_histo_nuke  = (MnvH1D*)Contin_mc_CinC_nuke->Clone("contin_histo_mc_cinc");
  MnvH1D* Contin_mc_SinC_histo_nuke  = (MnvH1D*)Contin_mc_SinC_nuke->Clone("contin_histo_mc_sinc");

   Contin_data_histo_nuke->Divide(Contin_data, Contin_data_nuke, 1,1,"B");
   Contin_mc_TinC_histo_nuke->Divide(Contin_mc_TinC, Contin_mc_TinC_nuke, 1,1,"B");
   Contin_mc_CinC_histo_nuke->Divide(Contin_mc_CinC, Contin_mc_CinC_nuke, 1,1,"B");
   Contin_mc_SinC_histo_nuke->Divide(Contin_mc_SinC, Contin_mc_SinC_nuke, 1,1,"B");

   TCanvas *c11 = new TCanvas("c11","Contin Data Ratio:NSF/OSF",600,500);
   c11->SetGrid(1);
   Contin_data_histo_nuke->Draw();
//   TCanvas *c22 = new TCanvas("c22","Trans in Contin Ratio:NSF/OSF",600,500);
//   c22->SetGrid(1);
   Contin_mc_TinC_histo_nuke->Draw("Psame][");
//   TCanvas *c33 = new TCanvas("c33","Contin in Contin Ratio:NSF/OSF",600,500);
//   c33->SetGrid(1);
   Contin_mc_CinC_histo_nuke->Draw("Psame][");
//   TCanvas *c44 = new TCanvas("c44","Signal in Contin Ratio:NSF/OSF",600,500);
//   c44->SetGrid(1);
   Contin_mc_SinC_histo_nuke->Draw("Psame][");

   Contin_data_histo_nuke   ->GetYaxis()->SetRangeUser(0.5,1.5);
   Contin_mc_TinC_histo_nuke->GetYaxis()->SetRangeUser(0.5,1.5);
   Contin_mc_CinC_histo_nuke->GetYaxis()->SetRangeUser(0.5,1.5);
   Contin_mc_SinC_histo_nuke->GetYaxis()->SetRangeUser(0.5,1.5);

         TLegend* leg = new TLegend(0.5,0.6,0.65,0.85);
         leg->AddEntry(Contin_data_histo_nuke,"Data.", "lep");
         leg->AddEntry(Contin_mc_TinC_histo_nuke,"Trans.", "lep");
         leg->AddEntry(Contin_mc_CinC_histo_nuke,"Contin.", "lep");
         leg->AddEntry(Contin_mc_SinC_histo_nuke,"Signal.", "lep");
         leg->Draw();
   gPad->SetGrid(1);
   gPad->Print("FluxValidation/Ratio_NSF_OSF_Contin_Emu_Iron.pdf","pdf");
   gPad->Print("FluxValidation/Ratio_NSF_OSF_Contin_Emu_Iron.png","png");
   gPad->Update();

}



}

