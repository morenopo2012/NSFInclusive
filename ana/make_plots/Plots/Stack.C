//******************************data MC***************************

{

  TFile *f1 = TFile::Open("$HISTS/v1/Hists_Energy2D_t1_z26_Nu_v1_noW.root");
  TFile *f2 = TFile::Open("$HISTS/v1/Hists_Energy2D_t1_z26_Nu_v1_wG.root");
  TFile *f3 = TFile::Open("$HISTS/v1/Hists_Energy2D_t1_z26_Nu_v1_wGF.root");
  TFile *f4 = TFile::Open("$HISTS/v1/Hists_Energy2D_t1_z26_Nu_v1_wGF2p2h.root");
  TFile *f5=  TFile::Open("$HISTS/v1/Hists_Energy2D_t1_z26_Nu_v1_wGF2p2hRPA.root");
  TFile *f6 = TFile::Open("$HISTS/v1/Hists_Energy2D_t1_z26_Nu_v1_wGF2p2hRPAMi.root");
  TFile *f7 = TFile::Open("$HISTS/v1/Hists_Energy2D_t1_z26_Nu_v1_wGF2p2hRPAMiNon.root");
  TFile *f8 = TFile::Open("$HISTS/v1/Hists_Energy2D_t1_z26_Nu_v1_wGF2p2hRPAMiNonNuE.root");
  TFile *f9 = TFile::Open("$HISTS/v1/Hists_Energy2D_t1_z26_Nu_v1_wGF2p2hRPAMiNonNuE.root");


  //**********************************************Z Mass************************
  TCanvas *c1 =  new TCanvas("c1","",323,144,534,523);
  c1->Draw();
  c1->cd();
  

  TH1F *h8_mc = (TH1F*)f8->Get("Ehad_data");
  TH1F *h1_mc = (TH1F*)f8->Get("Ehad_dis");
  TH1F *h2_mc = (TH1F*)f8->Get("Ehad_lowW");
  TH1F *h3_mc = (TH1F*)f8->Get("Ehad_lowQ2");


  h1_mc->SetLineColor(1);
  h1_mc->SetFillColor(kBlue-2);

  //h2_mc->SetLineColor(2);
  //h2_mc->SetLineColor(kOrange-2);
  //h2_mc->SetFillColor(kOrange-2);
  h2_mc->SetLineColor(2);  //Purple
  h2_mc->SetFillColor(kGreen);
  h3_mc->SetLineColor(51);
  h3_mc->SetFillColor(kBlue-2);
  h8_mc->SetLineColor(kRed-2);
  h8_mc->SetFillColor(kBlue-2);

  h1_mc->SetStats(0);
  h2_mc->SetStats(0);
  h3_mc->SetStats(0);
  h8_mc->SetStats(0);


  h8_mc->Scale(1e-3, "width");
  h1_mc->Scale(1e-3, "width");
  h2_mc->Scale(1e-3, "width");
  h3_mc->Scale(1e-3, "width");
   
  h1_mc->Scale(0.449);
  h2_mc->Scale(0.449);
  h3_mc->Scale(0.449);
 
  THStack *hStack1 = new THStack("hs","");
  hStack1->Add(h1_mc);
  hStack1->Add(h2_mc);
  hStack1->Add(h3_mc);
  
  hStack1->Draw("hist");
  h8_mc->Draw("same p");

  hStack1->GetYaxis()->SetTitle("Number of Events");
  hStack1->GetYaxis()->SetTitleOffset(0.90);
  hStack1->GetYaxis()->SetTitleSize(0.04);


  hStack1->GetYaxis()->SetRangeUser(0.,2.);
  hStack1->GetYaxis()->SetTitle("Events/10^3");
  hStack1->GetXaxis()->SetTitle("Ehad_dis");
  hStack1->SetTitle("Iron of Target 1"); 
  

  TLegend *leg1 = new TLegend(0.6905444,0.5759494,0.8825215,0.8417722,NULL,"brNDC");
  leg1->AddEntry(h8_mc,"Data","lp");
  leg1->AddEntry(h1_mc,"DIS","l");
  leg1->AddEntry(h2_mc,"Low W","l");
  leg1->AddEntry(h3_mc,"Low Q2","l");
  
  leg1->AddEntry(" ", "POT Normalized");
  leg1->SetFillColor(0);
  leg1->SetLineColor(1);
  leg1->SetTextSize(0.025);
  leg1->SetTextColor(kRed);
  leg1->Draw();
  c1->cd();

  c1->Draw();
  c1->Update();
  c1->SaveAs("Ehad.png");

}
