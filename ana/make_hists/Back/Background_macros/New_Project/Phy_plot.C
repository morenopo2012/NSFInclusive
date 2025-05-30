#include "PlotUtils/MnvH1D.h"
#include "TFile.h"
#include "TArrow.h"
#include "TLegend.h"
#include "TCanvas.h"
using namespace PlotUtils;

void Phy_plot(){ 

   string target = "t5_z26";
   string var = "Emu";

  TString hname1 = Form("/minerva/data/users/fakbar/NukeHists/v1_MAD/Hists_EventSelection_phys_without_sys_%s_Nu_v1_MAD_.root", target.c_str());
  cout<<hname1<<endl;
  TFile *f1 = new TFile( hname1,"read" );

  TString hname2 = Form("/minerva/data/users/fakbar/NukeHists/v1_nukeCC_validation/minervame1A/Hists_SideBands_%s_Nu_v1_nukeCC_validation.root", target.c_str());
  cout<<hname2<<endl;
  TFile *f2 = new TFile( hname2,"read" );

   f1->cd();
   MnvH1D *TinT   = (MnvH1D*)f1->Get(Form("hists_trans_in_trans_%s", var.c_str() ));
   MnvH1D *CinT   = (MnvH1D*)f1->Get(Form("hists_contin_in_trans_%s", var.c_str()));
   MnvH1D *SinT   = (MnvH1D*)f1->Get(Form("hists_signal_in_trans_%s", var.c_str()));
   MnvH1D *TinC   = (MnvH1D*)f1->Get(Form("hists_trans_in_contin_%s", var.c_str()));
   MnvH1D *CinC   = (MnvH1D*)f1->Get(Form("hists_contin_in_contin_%s", var.c_str()));
   MnvH1D *SinC   = (MnvH1D*)f1->Get(Form("hists_signal_in_contin_%s", var.c_str()));

   MnvH1D *DatainT   = (MnvH1D*)f1->Get(Form("selected_data_reco_trans_%s", var.c_str()));
   MnvH1D *DatainC   = (MnvH1D*)f1->Get(Form("selected_data_reco_contin_%s", var.c_str()));

   MnvH1D *USinT   = (MnvH1D*)f1->Get(Form("hists_trans_mc_CH_US_%s", var.c_str()));
   MnvH1D *USinC   = (MnvH1D*)f1->Get(Form("hists_contin_mc_CH_US_%s", var.c_str()));
   MnvH1D *DSinT   = (MnvH1D*)f1->Get(Form("hists_trans_mc_CH_DS_%s", var.c_str()));
   MnvH1D *DSinC   = (MnvH1D*)f1->Get(Form("hists_contin_mc_CH_DS_%s", var.c_str()));

   f2->cd();
   MnvH1D *TinT_nuke   = (MnvH1D*)f2->Get(Form("MC_true_trans_in_trans_%s_%s",   var.c_str(), target.c_str()));
   MnvH1D *CinT_nuke   = (MnvH1D*)f2->Get(Form("MC_true_contin_in_trans_%s_%s",  var.c_str(), target.c_str()));
   MnvH1D *SinT_nuke   = (MnvH1D*)f2->Get(Form("MC_true_signal_in_trans_%s_%s",  var.c_str(), target.c_str()));
   MnvH1D *TinC_nuke   = (MnvH1D*)f2->Get(Form("MC_true_trans_in_contin_%s_%s",  var.c_str(), target.c_str()));
   MnvH1D *CinC_nuke   = (MnvH1D*)f2->Get(Form("MC_true_contin_in_contin_%s_%s", var.c_str(), target.c_str()));
   MnvH1D *SinC_nuke   = (MnvH1D*)f2->Get(Form("MC_true_signal_in_contin_%s_%s", var.c_str(), target.c_str()));

   MnvH1D *DatainT_nuke   = (MnvH1D*)f2->Get(Form("Data_trans_%s_%s", var.c_str(), target.c_str()));
   MnvH1D *DatainC_nuke   = (MnvH1D*)f2->Get(Form("Data_contin_%s_%s", var.c_str(), target.c_str()));

   MnvH1D *USinT_nuke   = (MnvH1D*)f2->Get(Form("MC_trans_CH_US_%s_%s", var.c_str(), target.c_str()));
   MnvH1D *USinC_nuke   = (MnvH1D*)f2->Get(Form("MC_contin_CH_US_%s_%s", var.c_str(), target.c_str()));
   MnvH1D *DSinT_nuke   = (MnvH1D*)f2->Get(Form("MC_trans_CH_DS_%s_%s",  var.c_str(), target.c_str()));
   MnvH1D *DSinC_nuke   = (MnvH1D*)f2->Get(Form("MC_contin_CH_DS_%s_%s", var.c_str(), target.c_str()));


//*******************************************************
//			TRANSITION
//*******************************************************
{
//MasterAnaDev
/*   TCanvas *c1 = new TCanvas("c1","Trans in Trans NSF:MasterAnaDev",600,500);
   c1->SetGrid(1);
   TinT ->Draw();
   TinT_nuke ->Draw("same");
   TinT->GetXaxis()->SetTitle("Emu (GeV)"); 
   TinT ->GetXaxis()->SetRangeUser(0,50);

   TCanvas *c2 = new TCanvas("c2","Contin in contin NSF:MasterAnaDev",600,500);
   c2->SetGrid(1);
   CinC ->Draw();
   CinC_nuke ->Draw("same");
   CinC->GetXaxis()->SetTitle("Emu (GeV)"); 
   CinC ->GetXaxis()->SetRangeUser(0,50);
*/
   TinT->SetLineColor(kBlack);       
   CinT->SetLineColor(kBlack);          
   SinT->SetLineColor(kBlack);       
 
   TinC->SetLineColor(kBlack);       
   CinC->SetLineColor(kBlack);          
   SinC->SetLineColor(kBlack);       
 
   TinT->SetMarkerColor(kBlack);       
   CinT->SetMarkerColor(kBlack);          
   SinT->SetMarkerColor(kBlack);       
 
   TinC->SetMarkerColor(kBlack);       
   CinC->SetMarkerColor(kBlack);          
   SinC->SetMarkerColor(kBlack);       
 
/*
   DatainT->SetLineColor(kBlack);
   TinT->SetLineColor(kOrange);       
   CinT->SetLineColor(kBlue);          
   SinT->SetLineColor(kGreen+1);       
 
   DatainC->SetLineColor(kBlack);
   TinC->SetLineColor(kOrange);       
   CinC->SetLineColor(kBlue);          
   SinC->SetLineColor(kGreen+1);       
 
   DatainT->SetMarkerColor(kBlack);
   TinT->SetMarkerColor(kOrange);       
   CinT->SetMarkerColor(kBlue);          
   SinT->SetMarkerColor(kGreen+1);       
 
   DatainC->SetMarkerColor(kBlack);
   TinC->SetMarkerColor(kOrange);       
   CinC->SetMarkerColor(kBlue);          
   SinC->SetMarkerColor(kGreen+1);       
 */
/*
   TCanvas *c3 = new TCanvas("c3","Contin in Trans OSF:NukeCC",600,500);
   c3->SetGrid(1);
   CinT->Draw();
   CinT_nuke->Draw("same");
   CinT->GetXaxis()->SetTitle("Emu (GeV)"); 
   CinT ->GetXaxis()->SetRangeUser(0,50);

   TCanvas *c4 = new TCanvas("c4","Trans in Contin OSF:NukeCC",600,500);
   c4->SetGrid(1);
   TinC->Draw();
   TinC_nuke->Draw("same");
   TinC->GetXaxis()->SetTitle("Emu (GeV)"); 
   TinC ->GetXaxis()->SetRangeUser(0,50);

   TCanvas *c5 = new TCanvas("c5","Signal in Contin OSF:NukeCC",600,500);
   c5->SetGrid(1);
   SinC->Draw();
   SinC_nuke->Draw("same");
   SinC->GetXaxis()->SetTitle("Emu (GeV)"); 
   SinC ->GetXaxis()->SetRangeUser(0,50);

   TCanvas *c6 = new TCanvas("c6","Signal in Trans OSF:NukeCC",600,500);
   c6->SetGrid(1);
   SinT->Draw();
   SinT_nuke->Draw("same");
   SinT->GetXaxis()->SetTitle("Emu (GeV)"); 
   SinT ->GetXaxis()->SetRangeUser(0,50);
*/
   DatainT_nuke->SetLineColor(kRed);
   TinT_nuke->SetLineColor(kRed);       
   CinT_nuke->SetLineColor(kRed);          
   SinT_nuke->SetLineColor(kRed);       
 
   DatainC_nuke->SetLineColor(kRed);
   TinC_nuke->SetLineColor(kRed);       
   CinC_nuke->SetLineColor(kRed);          
   SinC_nuke->SetLineColor(kRed);       
 
   DatainT_nuke->SetMarkerColor(kRed);
   TinT_nuke->SetMarkerColor(kRed);       
   CinT_nuke->SetMarkerColor(kRed);          
   SinT_nuke->SetMarkerColor(kRed);       
 
   DatainC_nuke->SetMarkerColor(kRed);
   TinC_nuke->SetMarkerColor(kRed);       
   CinC_nuke->SetMarkerColor(kRed);          
   SinC_nuke->SetMarkerColor(kRed);       
                                                 
/*
   DatainT_nuke->SetLineColor(kBlack);
   TinT_nuke->SetLineColor(kOrange);       
   CinT_nuke->SetLineColor(kBlue);          
   SinT_nuke->SetLineColor(kGreen+1);       
 
   DatainC_nuke->SetLineColor(kBlack);
   TinC_nuke->SetLineColor(kOrange);       
   CinC_nuke->SetLineColor(kBlue);          
   SinC_nuke->SetLineColor(kGreen+1);       
 
   DatainT_nuke->SetMarkerColor(kBlack);
   TinT_nuke->SetMarkerColor(kOrange);       
   CinT_nuke->SetMarkerColor(kBlue);          
   SinT_nuke->SetMarkerColor(kGreen+1);       
 
   DatainC_nuke->SetMarkerColor(kBlack);
   TinC_nuke->SetMarkerColor(kOrange);       
   CinC_nuke->SetMarkerColor(kBlue);          
   SinC_nuke->SetMarkerColor(kGreen+1);       
*/              
//Ratio = MAD/NukeCC

  MnvH1D* DataInTrans_nuke       = (MnvH1D*)DatainT_nuke->Clone("DatainTrans_nuke");
  DataInTrans_nuke->Divide(DatainT, DatainT_nuke, 1,1,"B");
  TCanvas *c11 = new TCanvas("c11","Trans Data Ratio:NSF/OSF",600,500);
  c11->SetGrid(1);
  DataInTrans_nuke->Draw();
  DataInTrans_nuke->GetYaxis()->SetRangeUser(0.99,1.01);
  DataInTrans_nuke->GetXaxis()->SetTitle("Emu (GeV)"); 
  DataInTrans_nuke ->GetXaxis()->SetRangeUser(0,50);
   gPad->SetGrid(1);
   gPad->Print(Form("Ratio_NSF_OSF_Data_in_Trans_5_26.pdf","pdf",target.c_str() ));
   gPad->Update();


  MnvH1D* DataInContin_nuke       = (MnvH1D*)DatainC_nuke->Clone("DatainContin_nuke");
  DataInContin_nuke->Divide(DatainC, DatainC_nuke, 1,1,"B");
  TCanvas *c12 = new TCanvas("c12","Contin Data Ratio:NSF/OSF",600,500);
  c12->SetGrid(1);
  DataInContin_nuke->Draw();
  DataInContin_nuke->GetYaxis()->SetRangeUser(0.99,1.01);
  DataInContin_nuke->GetXaxis()->SetTitle("Emu (GeV)"); 
  DataInContin_nuke->GetXaxis()->SetRangeUser(0,50);
   gPad->SetGrid(1);
   gPad->Print(Form("Ratio_NSF_OSF_Data_in_Contin_5_26.pdf","pdf",target.c_str() ));
   gPad->Update();


  MnvH1D* TransInContin_nuke       = (MnvH1D*)TinC_nuke->Clone("TransinContin_nuke");
  TransInContin_nuke->Divide(TinC, TinC_nuke, 1,1,"B");
  TCanvas *c13 = new TCanvas("c13","Trans in contin Ratio:NSF/OSF",600,500);
  c13->SetGrid(1);
  TransInContin_nuke->Draw();
  TransInContin_nuke->GetYaxis()->SetRangeUser(0.99,1.01);
  TransInContin_nuke->GetXaxis()->SetTitle("Emu (GeV)"); 
  TransInContin_nuke->GetXaxis()->SetRangeUser(0,50);
   gPad->SetGrid(1);
   gPad->Print(Form("Ratio_NSF_OSF_Trans_in_Contin_5_26.pdf","pdf",target.c_str() ));
   gPad->Update();

  MnvH1D* ContinInContin_nuke       = (MnvH1D*)CinC_nuke->Clone("ContininContin_nuke");
  ContinInContin_nuke->Divide(CinC, CinC_nuke, 1,1,"B");
  TCanvas *c14 = new TCanvas("c14","Contin in contin Ratio:NSF/OSF",600,500);
  c14->SetGrid(1);
  ContinInContin_nuke->Draw();
  ContinInContin_nuke->GetYaxis()->SetRangeUser(0.99,1.01);
  ContinInContin_nuke->GetXaxis()->SetTitle("Emu (GeV)"); 
  ContinInContin_nuke->GetXaxis()->SetRangeUser(0,50);
   gPad->SetGrid(1);
   gPad->Print(Form("Ratio_NSF_OSF_Contin_in_Contin_5_26.pdf","pdf",target.c_str() ));
   gPad->Update();


  MnvH1D* TransInTrans_nuke       = (MnvH1D*)TinT_nuke->Clone("TransinTrans_nuke");
  TransInTrans_nuke->Divide(TinT, TinT_nuke, 1,1,"B");
  TCanvas *c15 = new TCanvas("c15","Trans in trans Ratio:NSF/OSF",600,500);
  c15->SetGrid(1);
  TransInTrans_nuke->Draw();
  TransInTrans_nuke->GetYaxis()->SetRangeUser(0.99,1.01);
  TransInTrans_nuke->GetXaxis()->SetTitle("Emu (GeV)"); 
  TransInTrans_nuke->GetXaxis()->SetRangeUser(0,50);
   gPad->SetGrid(1);
   gPad->Print(Form("Ratio_NSF_OSF_Trans_in_Trans_5_26.pdf","pdf",target.c_str() ));
   gPad->Update();


  MnvH1D* ContinInTrans_nuke       = (MnvH1D*)CinT_nuke->Clone("ContininTrans_nuke");
  ContinInTrans_nuke->Divide(CinT, CinT_nuke, 1,1,"B");
  TCanvas *c16 = new TCanvas("c16","Contin in trans Ratio:NSF/OSF",600,500);
  c16->SetGrid(1);
  ContinInTrans_nuke->Draw();
  ContinInTrans_nuke->GetYaxis()->SetRangeUser(0.99,1.01);
  ContinInTrans_nuke->GetXaxis()->SetTitle("Emu (GeV)"); 
  ContinInTrans_nuke->GetXaxis()->SetRangeUser(0,50);
   gPad->SetGrid(1);
   gPad->Print(Form("Ratio_NSF_OSF_Contin_in_Trans_5_26.pdf","pdf",target.c_str() ));
   gPad->Update();

  MnvH1D* SigInTrans_nuke       = (MnvH1D*)SinT_nuke->Clone("SignalinTrans_nuke");
  SigInTrans_nuke->Divide(SinT, SinT_nuke, 1,1,"B");
  TCanvas *c17 = new TCanvas("c17","Signal in trans Ratio:NSF/OSF",600,500);
  c17->SetGrid(1);
  SigInTrans_nuke->Draw();
  SigInTrans_nuke->GetYaxis()->SetRangeUser(0.99,1.01);
  SigInTrans_nuke->GetXaxis()->SetTitle("Emu (GeV)"); 
  SigInTrans_nuke->GetXaxis()->SetRangeUser(0,50);
   gPad->SetGrid(1);
   gPad->Print(Form("Ratio_NSF_OSF_Signal_in_Trans_5_26.pdf","pdf",target.c_str() ));
   gPad->Update();

  MnvH1D* SigInContin_nuke       = (MnvH1D*)SinC_nuke->Clone("SignalinContin_nuke");
  SigInContin_nuke->Divide(SinC, SinC_nuke, 1,1,"B");
  TCanvas *c18 = new TCanvas("c18","Signal in Contin Ratio:NSF/OSF",600,500);
  c18->SetGrid(1);
  SigInContin_nuke->Draw();
  SigInContin_nuke->GetYaxis()->SetRangeUser(0.99,1.01);
  SigInContin_nuke->GetXaxis()->SetTitle("Emu (GeV)"); 
  SigInContin_nuke->GetXaxis()->SetRangeUser(0,50);
   gPad->SetGrid(1);
   gPad->Print(Form("Ratio_NSF_OSF_Signal_in_Contin_5_26.pdf","pdf",target.c_str() ));
   gPad->Update();

cout<<target.c_str()<<endl;



/*
         TLegend* leg = new TLegend(0.5,0.6,0.65,0.85);
         leg->AddEntry(US_data_histo_nuke,"Data.", "lep");
         leg->AddEntry(US_signal_histo_nuke,"Material.", "lep");
         leg->AddEntry(US_other_histo_nuke,"Other", "lep");
         leg->AddEntry(US_plastic_us_histo_nuke,"US", "lep");
         leg->AddEntry(US_plastic_ds_histo_nuke,"DS", "lep");
         leg->Draw();
   gPad->SetGrid(1);
   gPad->Print("Flu%sValidation/Ratio_NSF_OSF_US_planeDNN_Iron.pdf","pdf");
   gPad->Print("Flu%sValidation/Ratio_NSF_OSF_US_planeDNN_Iron.pdf","pdf");
   gPad->Update();
*/
}
}

