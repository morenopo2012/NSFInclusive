#include "PlotUtils/MnvH1D.h"
#include "TFile.h"
#include "TArrow.h"
#include "TLegend.h"
#include "TCanvas.h"
using namespace PlotUtils;

void plot(){ 

  TString hname1 = Form("/minerva/data/users/fakbar/NukeHists/v1_NuMAD/Hists_PlasticBackgd_with_SYS_FullDet_q2WfromBranch_ME1A_t1_z26_Nu_v1_NuMAD.root");
  cout<<hname1<<endl;
  TFile *f1 = new TFile( hname1,"read" );

  TString hname2 = Form("/pnfs/minerva/persistent/users/afilkins/NukeHists/Mnv080321_NuEOff_fulldetectorSBs/minervame1A/Hists_PlasticSideband_t1_z26_Nu_Mnv080321_NuEOff_fulldetectorSBs.root");
  cout<<hname2<<endl;
  TFile *f2 = new TFile( hname2,"read" );

   f1->cd();
   MnvH1D *US_data_var   = (MnvH1D*)f1->Get("selected_data_reco_US_planeDNN");
   MnvH1D *US_mat_var    = (MnvH1D*)f1->Get("US_Fe_Iron_planeDNN");
   MnvH1D *US_other_var  = (MnvH1D*)f1->Get("US_other_Iron_planeDNN");
   MnvH1D *US_regUS_var  = (MnvH1D*)f1->Get("US_regUS_Iron_planeDNN");
   MnvH1D *US_regDS_var  = (MnvH1D*)f1->Get("US_regDS_Iron_planeDNN");

   MnvH1D *DS_data_var   = (MnvH1D*)f1->Get("selected_data_reco_DS_planeDNN");
   MnvH1D *DS_mat_var    = (MnvH1D*)f1->Get("DS_Fe_Iron_planeDNN");
   MnvH1D *DS_other_var  = (MnvH1D*)f1->Get("DS_other_Iron_planeDNN");
   MnvH1D *DS_regUS_var  = (MnvH1D*)f1->Get("DS_regUS_Iron_planeDNN");
   MnvH1D *DS_regDS_var  = (MnvH1D*)f1->Get("DS_regDS_Iron_planeDNN");

   f2->cd();
   MnvH1D *US_data_var_nuke   = (MnvH1D*)f2->Get("DataUSSideband_planeDNN_t1_z26");
   MnvH1D *US_mat_var_nuke    = (MnvH1D*)f2->Get("hFe_planeDNNMCUSScintSideband_t1_z26");
   MnvH1D *US_other_var_nuke  = (MnvH1D*)f2->Get("hOther_planeDNNMCUSScintSideband_t1_z26");
   MnvH1D *US_regUS_var_nuke  = (MnvH1D*)f2->Get("hUS_planeDNNMCUSScintSideband_t1_z26");
   MnvH1D *US_regDS_var_nuke  = (MnvH1D*)f2->Get("hDS_planeDNNMCUSScintSideband_t1_z26");

   MnvH1D *DS_data_var_nuke   = (MnvH1D*)f2->Get("DataDSSideband_planeDNN_t1_z26");
   MnvH1D *DS_mat_var_nuke    = (MnvH1D*)f2->Get("hFe_planeDNNMCDSScintSideband_t1_z26");
   MnvH1D *DS_other_var_nuke  = (MnvH1D*)f2->Get("hOther_planeDNNMCDSScintSideband_t1_z26");
   MnvH1D *DS_regUS_var_nuke  = (MnvH1D*)f2->Get("hUS_planeDNNMCDSScintSideband_t1_z26");
   MnvH1D *DS_regDS_var_nuke  = (MnvH1D*)f2->Get("hDS_planeDNNMCDSScintSideband_t1_z26");

//*******************************************************
//			UPSTREAM
//*******************************************************
{
//MasterAnaDev
/*   TCanvas *c1 = new TCanvas("c1","US NSF:MasterAnaDev",600,500);
   c1->SetGrid(1);
   US_data_var->Draw();
   US_mat_var->Draw("Psame][");
   US_other_var->Draw("Psame][");
   US_regUS_var->Draw("Psame][");
   US_regDS_var->Draw("Psame][");
*/
   US_data_var->SetMarkerColor(kBlack);
   US_mat_var->SetMarkerColor(kRed);
   US_other_var->SetMarkerColor(kBlue);
   US_regUS_var->SetMarkerColor(kGreen+2);
   US_regDS_var->SetMarkerColor(kOrange+2);

   US_data_var->SetLineColor(kBlack);
   US_mat_var->SetLineColor(kRed);
   US_other_var->SetLineColor(kBlue);
   US_regUS_var->SetLineColor(kGreen+2);
   US_regDS_var->SetLineColor(kOrange+2);

//NukeCCInclusive
/*   TCanvas *c2 = new TCanvas("c2","US OSF:NukeCCInclusive",600,500);
   c2->SetGrid(1);
   US_data_var_nuke->Draw();
   US_mat_var_nuke->Draw("Psame][");
   US_other_var_nuke->Draw("Psame][");
   US_regUS_var_nuke->Draw("Psame][");
   US_regDS_var_nuke->Draw("Psame][");
*/
   US_data_var_nuke->SetMarkerColor(kBlack);
   US_mat_var_nuke->SetMarkerColor(kRed);
   US_other_var_nuke->SetMarkerColor(kBlue);
   US_regUS_var_nuke->SetMarkerColor(kGreen+2);
   US_regDS_var_nuke->SetMarkerColor(kOrange+2);

   US_data_var_nuke->SetLineColor(kBlack);
   US_mat_var_nuke->SetLineColor(kRed);
   US_other_var_nuke->SetLineColor(kBlue);
   US_regUS_var_nuke->SetLineColor(kGreen+2);
   US_regDS_var_nuke->SetLineColor(kOrange+2);

//Ratio = MAD/NukeCC
  MnvH1D* US_data_histo_nuke       = (MnvH1D*)US_data_var_nuke->Clone("US_histo_data");
  MnvH1D* US_signal_histo_nuke     = (MnvH1D*)US_mat_var_nuke->Clone("US_histo_mc_signal");
  MnvH1D* US_other_histo_nuke      = (MnvH1D*)US_other_var_nuke->Clone("US_histo_mc_other_targ");
  MnvH1D* US_plastic_us_histo_nuke = (MnvH1D*)US_regUS_var_nuke->Clone("US_histo_mc_plasticUS");
  MnvH1D* US_plastic_ds_histo_nuke = (MnvH1D*)US_regDS_var_nuke->Clone("US_histo_mc_plasticDS");

   US_data_histo_nuke->Divide(US_data_var, US_data_var_nuke, 1,1,"B");
   US_signal_histo_nuke->Divide(US_mat_var, US_mat_var_nuke, 1,1,"B");
   US_other_histo_nuke->Divide(US_other_var, US_other_var_nuke, 1,1,"B");
   US_plastic_us_histo_nuke->Divide(US_regUS_var, US_regUS_var_nuke, 1,1,"B");
   US_plastic_ds_histo_nuke->Divide(US_regDS_var, US_regDS_var_nuke, 1,1,"B");

   TCanvas *c1 = new TCanvas("c1","US Data Ratio:NSF/OSF",600,500);
   c1->SetGrid(1);
   US_data_histo_nuke->Draw();
//   TCanvas *c2 = new TCanvas("c2","US Material Ratio:NSF/OSF",600,500);
//   c2->SetGrid(1);
   US_signal_histo_nuke->Draw("Psame][");
//  TCanvas *c3 = new TCanvas("c3","US Other Ratio:NSF/OSF",600,500);
//   c3->SetGrid(1);
   US_other_histo_nuke->Draw("Psame][");
//   TCanvas *c4 = new TCanvas("c4","US regUS Ratio:NSF/OSF",600,500);
//   c4->SetGrid(1);
   US_plastic_us_histo_nuke->Draw("Psame][");
//   TCanvas *c5 = new TCanvas("c5","US regDS Ratio:NSF/OSF",600,500);
//   c5->SetGrid(1);
   US_plastic_ds_histo_nuke->Draw("Psame][");
 
   US_data_histo_nuke->GetYaxis()->SetRangeUser(0.99,1.01);
   US_signal_histo_nuke->GetYaxis()->SetRangeUser(0.99,1.01);
   US_other_histo_nuke->GetYaxis()->SetRangeUser(0.99,1.01);
   US_plastic_us_histo_nuke->GetYaxis()->SetRangeUser(0.99,1.01);
   US_plastic_ds_histo_nuke->GetYaxis()->SetRangeUser(0.99,1.01);
   US_data_histo_nuke->GetXaxis()->SetRangeUser(0,70);
   US_data_histo_nuke->GetXaxis()->SetTitle("planeDNN"); 

         TLegend* leg = new TLegend(0.5,0.6,0.65,0.85);
         leg->AddEntry(US_data_histo_nuke,"Data.", "lep");
         leg->AddEntry(US_signal_histo_nuke,"Material.", "lep");
         leg->AddEntry(US_other_histo_nuke,"Other", "lep");
         leg->AddEntry(US_plastic_us_histo_nuke,"US", "lep");
         leg->AddEntry(US_plastic_ds_histo_nuke,"DS", "lep");
         leg->Draw();
   gPad->SetGrid(1);
   gPad->Print("FluxValidation/Ratio_NSF_OSF_US_planeDNN_Iron.pdf","pdf");
   gPad->Print("FluxValidation/Ratio_NSF_OSF_US_planeDNN_Iron.png","png");
   gPad->Update();

}
//*******************************************************
//			DOWNSTREAM
//*******************************************************
{
//MasterAnaDev
//   TCanvas *c11 = new TCanvas("c11","DS NSF:MasterAnaDev",600,500);
//   c11->SetGrid(1);
//   DS_data_var->Draw();
//   DS_mat_var->Draw("Psame][");
//   DS_other_var->Draw("Psame][");
//   DS_regUS_var->Draw("Psame][");
//   DS_regDS_var->Draw("Psame][");

   DS_data_var->SetMarkerColor(kBlack);
   DS_mat_var->SetMarkerColor(kRed);
   DS_other_var->SetMarkerColor(kBlue);
   DS_regUS_var->SetMarkerColor(kGreen+2);
   DS_regDS_var->SetMarkerColor(kOrange+2);

   DS_data_var->SetLineColor(kBlack);
   DS_mat_var->SetLineColor(kRed);
   DS_other_var->SetLineColor(kBlue);
   DS_regUS_var->SetLineColor(kGreen+2);
   DS_regDS_var->SetLineColor(kOrange+2);

//NukeCCInclusive
//   TCanvas *c22 = new TCanvas("c22","DS OSF:NukeCCInclusive",600,500);
//   c22->SetGrid(1);
//   DS_data_var_nuke->Draw();
//   DS_mat_var_nuke->Draw("Psame][");
//   DS_other_var_nuke->Draw("Psame][");
//   DS_regUS_var_nuke->Draw("Psame][");
//   DS_regDS_var_nuke->Draw("Psame][");

   DS_data_var_nuke->SetMarkerColor(kBlack);
   DS_mat_var_nuke->SetMarkerColor(kRed);
   DS_other_var_nuke->SetMarkerColor(kBlue);
   DS_regUS_var_nuke->SetMarkerColor(kGreen+2);
   DS_regDS_var_nuke->SetMarkerColor(kOrange+2);

   DS_data_var_nuke->SetLineColor(kBlack);
   DS_mat_var_nuke->SetLineColor(kRed);
   DS_other_var_nuke->SetLineColor(kBlue);
   DS_regUS_var_nuke->SetLineColor(kGreen+2);
   DS_regDS_var_nuke->SetLineColor(kOrange+2);

//Ratio = MAD/NukeCC
  MnvH1D* DS_data_histo_nuke       = (MnvH1D*)DS_data_var_nuke->Clone("DS_histo_data");
  MnvH1D* DS_signal_histo_nuke     = (MnvH1D*)DS_mat_var_nuke->Clone("DS_histo_mc_signal");
  MnvH1D* DS_other_histo_nuke      = (MnvH1D*)DS_other_var_nuke->Clone("DS_histo_mc_other_targ");
  MnvH1D* DS_plastic_us_histo_nuke = (MnvH1D*)DS_regUS_var_nuke->Clone("DS_histo_mc_plasticUS");
  MnvH1D* DS_plastic_ds_histo_nuke = (MnvH1D*)DS_regDS_var_nuke->Clone("DS_histo_mc_plasticDS");

   DS_data_histo_nuke->Divide(DS_data_var, DS_data_var_nuke, 1,1,"B");
   DS_signal_histo_nuke->Divide(DS_mat_var, DS_mat_var_nuke, 1,1,"B");
   DS_other_histo_nuke->Divide(DS_other_var, DS_other_var_nuke, 1,1,"B");
   DS_plastic_us_histo_nuke->Divide(DS_regUS_var, DS_regUS_var_nuke, 1,1,"B");
   DS_plastic_ds_histo_nuke->Divide(DS_regDS_var, DS_regDS_var_nuke, 1,1,"B");

   TCanvas *c11 = new TCanvas("c11","DS Data Ratio:NSF/OSF",600,500);
   c11->SetGrid(1);
   DS_data_histo_nuke->Draw();
//   TCanvas *c22 = new TCanvas("c22","DS Material Ratio:NSF/OSF",600,500);
//   c22->SetGrid(1);
   DS_signal_histo_nuke->Draw("Psame][");
//   TCanvas *c33 = new TCanvas("c33","DS Other Ratio:NSF/OSF",600,500);
//   c33->SetGrid(1);
   DS_other_histo_nuke->Draw("Psame][");
//   TCanvas *c44 = new TCanvas("c44","DS regUS Ratio:NSF/OSF",600,500);
//   c44->SetGrid(1);
   DS_plastic_us_histo_nuke->Draw("Psame][");
//   TCanvas *c55 = new TCanvas("c55","DS regDS Ratio:NSF/OSF",600,500);
//   c55->SetGrid(1);
   DS_plastic_ds_histo_nuke->Draw("Psame][");
 
   DS_data_histo_nuke->GetYaxis()->SetRangeUser(0.99,1.01);
   DS_signal_histo_nuke->GetYaxis()->SetRangeUser(0.99,1.01);
   DS_other_histo_nuke->GetYaxis()->SetRangeUser(0.99,1.01);
   DS_plastic_us_histo_nuke->GetYaxis()->SetRangeUser(0.99,1.01);
   DS_plastic_ds_histo_nuke->GetYaxis()->SetRangeUser(0.99,1.01);
   DS_data_histo_nuke->GetXaxis()->SetRangeUser(0,70);
   DS_data_histo_nuke->GetXaxis()->SetTitle("planeDNN"); 
//*************************************************************************

//   TCanvas *c4 = new TCanvas("c4","Legend",600,500);
         TLegend* leg = new TLegend(0.5,0.6,0.65,0.85);
         leg->AddEntry(US_data_histo_nuke,"Data.", "lep");
         leg->AddEntry(US_signal_histo_nuke,"Material.", "lep");
         leg->AddEntry(US_other_histo_nuke,"Other", "lep");
         leg->AddEntry(US_plastic_us_histo_nuke,"US", "lep");
         leg->AddEntry(US_plastic_ds_histo_nuke,"DS", "lep");
         leg->Draw();

   gPad->SetGrid(1);
   gPad->Print("FluxValidation/Ratio_NSF_OSF_DS_planeDNN_Iron.pdf","pdf");
   gPad->Print("FluxValidation/Ratio_NSF_OSF_DS_planeDNN_Iron.png","png");
   gPad->Update();

  }
}

