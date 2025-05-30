//******************************data MC***************************

{
    TFile *f1 = TFile::Open("Hists_EventSelection_t1_z26_Nu.root");
    //TFile *f1 = TFile::Open("/minerva/app/users/zdar/cmtuser/Minerva_v22r1p1_MADNew/Ana/NSFNukeCCInclusive/ana/make_hists/Background_macros/HistsFullTargetIronME1A/Hists_PhysicsBackgd_with_SYS_FullDet_q2WfromBranch_ME1A_targetsCombined_t2_z26_Nu.root");
  
  TFile *f2 = TFile::Open("/pnfs/minerva/persistent/users/afilkins/NukeHists/Mnv210924_git/minervame1L/Hists_Energy_Data_t1_z26_Nu_Mnv210924_git.root");
  //TFile *f1 = TFile::Open("/minerva/data/users/afilkins/NukeHists/MuonKludged/minervame1A/Hists_Energy_MC_t1_z26_Nu_MuonKludged.root");
  //  TFile *f1 = TFile::Open("$HISTS/v1/Hists_Energy_t3_z06_Nu_v1_.root");

  //**********************************************Z Mass************************
  TCanvas *c1 =  new TCanvas("c1","",323,144,534,523);
  c1->Draw();
  c1->cd();
  
  c1_2 = new TPad("c1_2","newpad",0.008032129,0.1866667,0.9879518,0.9911111);
  c1_2->SetBottomMargin(0.1);
  c1_2->Draw();
  c1_2->cd();
  c1_2->SetTickx(1);
  c1_2->SetTicky(1);
  TH1F *h1_mc = (TH1F*)f1->Get("h_data_Enu");
  //TH1F *h1_data = (TH1F*)f2->Get("hTransition_EnuMC_trans_t1_z26");
  TH1F *h1_data = (TH1F*)f2->Get("sample_dis_Data_Enu_t1_z26");

  //TH1F *h1_mc = (TH1F*)f1->Get("DS_regDS_Iron_Enu");
  //TH1F *h1_data = (TH1F*)f2->Get("hDS_EnuMCDSScintSideband_t1_z26");
  
  //TH1F *h1_mc = (TH1F*)f1->Get("selected_data_reco_DS_Enu");
  //TH1F *h1_data = (TH1F*)f2->Get("DataDSSideband_Enu_t1_z26"); 

  //TH1F *h1_mc = (TH1F*)f1->Get("selected_data_reco_contin_Ehad");
  //TH1F *h1_data = (TH1F*)f2->Get("Data_contin_Ehad_t1_z26"); 
  
  h1_mc->SetLineStyle(1);
  h1_mc->SetLineColor(1);
 

  //h1_data->GetXaxis()->SetLimits(1,50);
  //h1_mc->GetXaxis()->SetLimits(1,50);
  //h1_data->GetXaxis()->SetRangeUser(2,50);
  //h1_data->SetLineColor(kRed-2);
  h1_mc->SetTitle("Chor"); 
  h1_data->SetLineColor(2);
  h1_data->SetFillColor(kBlue-2);
  //h1_mc->SetLineColor(3);
  //h1_mc->SetFillColor(kBlue-2);
  h1_mc->SetStats(0);
  h1_data->SetStats(0);
  //h1_mc->GetXaxis()->SetRangeUser(0.,10.);
  //h1_mc->GetYaxis()->SetRangeUser(0.,200);
  h1_mc->GetYaxis()->SetRangeUser(0.5, 1.5);
  //h1_mc->GetXaxis()->SetTitle("Reconstructed Q^{2}");
  //h1_mc->GetYaxis()->SetTitle("Data/ MnvGENIE");
  //h1_mc->Draw("hist");
  h1_data->Divide(h1_mc);
  h1_data->GetXaxis()->SetTitle("Reconstructed Ehad");
  h1_data->GetYaxis()->SetTitle("MC Contin NSF/NukeCC");
  h1_data->Draw();
  //TLegend *leg1 = new TLegend(0.6905444,0.5759494,0.7825215,0.8417722,NULL,"brNDC");
  TLegend *leg1 = new TLegend(0.5905444,0.5759494,0.6825215,0.8417722,NULL,"brNDC");
  leg1->AddEntry(h1_mc,"MC", "l");
  leg1->AddEntry(h1_data,"Data");
//  leg1->AddEntry("" ,"POT Normalized");
//  leg1->AddEntry("" ,"Iron of Target 1");
  leg1->SetFillColor(0);
  leg1->SetLineColor(1);
  leg1->SetTextSize(0.02888889);
  leg1->SetTextColor(4);
  //leg1->Draw();

  c1->cd();
  c1->Draw();
  c1->Update();
  c1->SaveAs("Ehad_mc.png");
}
