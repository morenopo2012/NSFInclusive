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

  c1->SetTicky(5);
  c1->SetGridx();
  c1->SetGridy();
  c1->cd();
  c1->Draw();

  TCanvas *c2 =  new TCanvas("c2"," ",323,144,534,523);
  c2->SetGridx(2);
  c2->SetGridy();
  c2->Draw();
  c2->cd();
  
  TH1F *h9_mc = (TH1F*)f9->Get("Q2_data");
  TH1F *h1_mc = (TH1F*)f1->Get("Q2_dis");
  TH1F *h2_mc = (TH1F*)f2->Get("Q2_dis");
  TH1F *h3_mc = (TH1F*)f3->Get("Q2_dis");
  TH1F *h4_mc = (TH1F*)f4->Get("Q2_dis");
  TH1F *h5_mc = (TH1F*)f5->Get("Q2_dis");
  TH1F *h6_mc = (TH1F*)f6->Get("Q2_dis");
  TH1F *h7_mc = (TH1F*)f7->Get("Q2_dis");
  TH1F *h8_mc = (TH1F*)f8->Get("Q2_dis");


  h4_mc->SetLineColor(1);
  h4_mc->SetFillColor(kBlue-2);

  h1_mc->SetStats(0);
  h2_mc->SetStats(0);
  h3_mc->SetStats(0);
  h4_mc->SetStats(0);
  h5_mc->SetStats(0);
  h6_mc->SetStats(0);
  h7_mc->SetStats(0);
  h8_mc->SetStats(0);
  h9_mc->SetStats(0);

  h9_mc->Scale(1e-3, "width");
  h1_mc->Scale(1e-3, "width");
  h2_mc->Scale(1e-3, "width");
  h3_mc->Scale(1e-3, "width");
  h4_mc->Scale(1e-3, "width");
  h5_mc->Scale(1e-3, "width");
  h6_mc->Scale(1e-3, "width");
  h7_mc->Scale(1e-3, "width");
  h8_mc->Scale(1e-3, "width");
  //h1_mc->Scale(0.45);
  h1_mc->Scale(0.449);
  h2_mc->Scale(0.449);
  h3_mc->Scale(0.449);
  h4_mc->Scale(0.449);
  h5_mc->Scale(0.449);
  h6_mc->Scale(0.449);
  h7_mc->Scale(0.449);
  h8_mc->Scale(0.449);
  
  h4_mc->GetYaxis()->SetRangeUser(0.,10.);
  h4_mc->GetYaxis()->SetTitle("Events/10^3");
  h4_mc->GetXaxis()->SetTitle("Q2_dis");
  h4_mc->SetTitle("Iron of Target 1"); 
  
  h4_mc->Draw("hist");
  h4_mc->Draw("same");


  TLegend *leg1 = new TLegend(0.6905444,0.5759494,0.8825215,0.8417722,NULL,"brNDC");
  leg1->AddEntry(h9_mc,"Data","lp");
  //leg1->AddEntry(h1_mc,"No Weight","l");
  //leg1->AddEntry(h2_mc,"Genie only","l");
  //leg1->AddEntry(h3_mc,"Genie+Flux","l");
  leg1->AddEntry(h4_mc,"Ge+Fl+2p2h","l");
  //leg1->AddEntry(h5_mc,"Ge+Fl+2p2h+RPA","l");
  //leg1->AddEntry(h6_mc,"Ge+Fl+2p2h+RPA+Minos","l");
  //leg1->AddEntry(h7_mc,"G+F+2p2h+RPA+Minos+NonRes","l");
  //leg1->AddEntry(h8_mc,"G+F+2p2h+RPA+Minos+NonRes+NuE","l");
  leg1->AddEntry(" ", "POT Normalized", "l");
  leg1->SetFillColor(0);
  leg1->SetLineColor(1);
  leg1->SetTextSize(0.025);
  leg1->SetTextColor(kRed);
  leg1->Draw();
  c1->cd();

  TH1F *h1ratio = (TH1F *)h9_mc->Clone();
  h1ratio->Divide(h4_mc);
     
  h1ratio->SetMarkerColor(kBlack);
  h1ratio->SetMarkerStyle(20);
  h1ratio->SetLineColor(kBlack);
  h1ratio->SetMarkerSize(0.8);

  h1ratio->SetTitle("Iron of target 1 ");  
  h1ratio->GetXaxis()->SetTitle("Q2_dis");
  h1ratio->GetXaxis()->SetLabelFont(42);
  h1ratio->GetYaxis()->SetTitle("Data/MC");
  h1ratio->GetYaxis()->SetRangeUser(0.,2.);
  h1ratio->GetYaxis()->SetNdivisions(10);
  h1ratio->Draw();
  

  TLine l1(1.,1.0,30.,1.0);
  l1->SetLineWidth(2);
  l1->SetLineColor(kRed);
  l1->SetLineStyle(1);
  l1->Draw();
  c1->Draw();
  c1->Update();
  c2->Draw();
  c2->Update();
  c1->SaveAs("Q2_dis_dist.png");
  c2->SaveAs("Q2_dis_ratio.png");

}
