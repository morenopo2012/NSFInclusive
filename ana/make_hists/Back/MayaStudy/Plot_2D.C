//******************************data MC***************************

{
    //TFile *f1 = TFile::Open("Hists_EventSelection_t1_z26_Nu_v1_.root");
    TFile *f1 = TFile::Open("../Hists_EventSelection_t1_z26_Nu_v1_.root");
  
  //TFile *f1 = TFile::Open("/minerva/data/users/afilkins/NukeHists/MuonKludged/minervame1A/Hists_Energy_Data_t1_z26_Nu_MuonKludged.root");
  //TFile *f1 = TFile::Open("/minerva/data/users/afilkins/NukeHists/MuonKludged/minervame1A/Hists_Energy_MC_t1_z26_Nu_MuonKludged.root");

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
  
  //TH2F *h1_mc = (TH2F*)f1->Get("selected_mc_reco2d_Emu_Ehad");
  //TH2F *h1_data = (TH2F*)f1->Get("selected_data2d_reco_Emu_Ehad");
  TH2F *h1_mc = (TH2F*)f1->Get("selected_mc_reco2d_W_Q2");
  TH2F *h1_data = (TH2F*)f1->Get("selected_data2d_reco_W_Q2");
  

  h1_mc->SetLineStyle(1);
  h1_mc->SetLineColor(1);

  //h1_data->SetLineColor(kRed-2);
  h1_mc->SetTitle("Chor"); 
  h1_data->SetLineColor(2);
  h1_data->SetFillColor(kBlue-2);
  h1_mc->SetStats(0);
  h1_data->SetStats(0);
  h1_data->Scale(1e-3, "width");
  h1_mc->Scale(1e-3, "width");
  //h1_mc->Scale(0.448);
  h1_mc->Scale(0.421);
 ///// h1_mc->Scale(0.433);
  //h1_mc->GetXaxis()->SetRangeUser(0.,10.);
  //h1_mc->GetYaxis()->SetRangeUser(0.,200);
  h1_mc->GetYaxis()->SetRangeUser(0.,15);
  //h1_mc->GetXaxis()->SetTitle("Reconstructed Q^{2}");
  //h1_mc->GetYaxis()->SetTitle("Data/ MnvGENIE");
  //h1_mc->SetTitle("Chor"); 
  //h1_mc->Draw("hist");
  //h1_data->Draw("same");
  h1_data->Divide(h1_mc);
  //h1_data->GetXaxis()->SetTitle("Muon Energy (GeV)");
  //h1_data->GetYaxis()->SetTitle("Hadronic Energy (GeV)");
  h1_data->GetXaxis()->SetTitle("Invariant Mass W (GeV)");
  h1_data->GetYaxis()->SetTitle("Q^{2} (GeV^{2})");
  h1_data->Draw("colz");
  //h1_data->cd(); 
  //c1->cd();
  //c1->Draw();
  //c1->Update();
  TFile f("WQ2_Iron.root","recreate");
  f->cd();
  h1_data->Write();
  c1->cd();
  c1->Draw();
  c1->Update();
  c1->SaveAs("Ehad_mc.png");
  //c1->SaveAs("Ehad.root");
}
