//******************************data MC***************************

{

  TFile *f1 = TFile::Open("$HISTS/v1/Hists_Energy_t1_z26_Nu_v1_G.root");
  TFile *f2 = TFile::Open("$HISTS/v1/Hists_Energy_t1_z26_Nu_v1_GF.root");
  TFile *f3 = TFile::Open("$HISTS/v1/Hists_Energy_t1_z26_Nu_v1_GF2.root");
  TFile *f4 = TFile::Open("$HISTS/v1/Hists_Energy_t1_z26_Nu_v1_GF2RP.root");
  TFile *f5 = TFile::Open("$HISTS/v1/Hists_Energy_t1_z26_Nu_v1_GF2RPM.root");
  TFile *f6 = TFile::Open("$HISTS/v1/Hists_Energy_t1_z26_Nu_v1_GF2RPMNu.root");
  TFile *f7 = TFile::Open("$HISTS/v1/Hists_Energy_t1_z26_Nu_v1_GF2RPMNuNon.root");
  TFile *f8 = TFile::Open("$HISTS/v1/Hists_Energy_t1_z26_Nu_v1_GF2RPMNuNon.root");

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

  TH1F *h8_mc = (TH1F*)f8->Get("Ehad_data");
  TH1F *h1_mc = (TH1F*)f1->Get("Ehad_mc");
  TH1F *h2_mc = (TH1F*)f2->Get("Ehad_mc");
  TH1F *h3_mc = (TH1F*)f3->Get("Ehad_mc");
  TH1F *h4_mc = (TH1F*)f4->Get("Ehad_mc");
  TH1F *h5_mc = (TH1F*)f5->Get("Ehad_mc");
  TH1F *h6_mc = (TH1F*)f6->Get("Ehad_mc");
  TH1F *h7_mc = (TH1F*)f7->Get("Ehad_mc");


  //h1_mc->SetMarkerColor(1);
  //h1_mc->SetMarkerStyle(20);
  //h1_mc->SetMarkerSize(0.8); 
  h4_mc->SetLineColor(1);
  h4_mc->SetFillColor(kBlue-2);
/*
  //h2_mc->SetLineColor(2);
  //h2_mc->SetLineColor(kOrange-2);
  //h2_mc->SetFillColor(kOrange-2);
  h3_mc->SetLineColor(2);  //Purple
  //h3_mc->SetFillColor(2);
  h4_mc->SetLineColor(kOrange-8);
  //h4_mc->SetFillColor(kOrange-8);
  h5_mc->SetLineColor(51);
  //h5_mc->SetFillColor(51);
  h6_mc->SetLineColor(kGreen-2);
 // h6_mc->SetFillColor(kGreen-2);
  h7_mc->SetLineColor(kOrange-8);
  //h7_mc->SetFillColor(kOrange-8);
  h8_mc->SetLineColor(kRed-2);
  h8_mc->SetFillColor(kBlue-2);

*/

  h1_mc->SetStats(0);
  h2_mc->SetStats(0);
  h3_mc->SetStats(0);
  h4_mc->SetStats(0);
  h5_mc->SetStats(0);
  h6_mc->SetStats(0);
  h7_mc->SetStats(0);
  h8_mc->SetStats(0);


  h8_mc->Scale(1e-3, "width");
  h1_mc->Scale(1e-3, "width");
  h2_mc->Scale(1e-3, "width");
  h3_mc->Scale(1e-3, "width");
  h4_mc->Scale(1e-3, "width");
  h5_mc->Scale(1e-3, "width");
  h6_mc->Scale(1e-3, "width");
  h7_mc->Scale(1e-3, "width");
  //h1_mc->Scale(0.45);
  h1_mc->Scale(0.449);
  h2_mc->Scale(0.449);
  h3_mc->Scale(0.449);
  h4_mc->Scale(0.449);
  h5_mc->Scale(0.449);
  h6_mc->Scale(0.449);
  h7_mc->Scale(0.449);
  
  h4_mc->GetYaxis()->SetRangeUser(0.,8.);
  h4_mc->GetYaxis()->SetTitle("Events/10^3");
  h4_mc->GetXaxis()->SetTitle("Ehad_mc");
  h4_mc->SetTitle("Iron of Target 1"); 
  
  h4_mc->Draw("hist");
  h8_mc->Draw("same");


  TLegend *leg1 = new TLegend(0.6905444,0.5759494,0.8825215,0.8417722,NULL,"brNDC");
  leg1->AddEntry(h8_mc,"Data","lp");
  //leg1->AddEntry(h1_mc,"Genie only","l");
  //leg1->AddEntry(h2_mc,"Genie+Flux","l");
  //leg1->AddEntry(h3_mc,"Ge+Fl+2p2h","l");
  leg1->AddEntry(h4_mc,"Ge+Fl+2p2h+RPA","l");
  //leg1->AddEntry(h5_mc,"Ge+Fl+2p2h+RPA+Minos","l");
  //leg1->AddEntry(h6_mc,"G+F+2p2h+RPA+Minos+NuE","l");
  //leg1->AddEntry(h7_mc,"G+F+2p2h+RPA+Minos+NuE+NoN","l");
  leg1->AddEntry(" ", "POT Normalized");
  leg1->SetFillColor(0);
  leg1->SetLineColor(1);
  leg1->SetTextSize(0.025);
  leg1->SetTextColor(kRed);
  leg1->Draw();
  c1->cd();

  TH1F *h1ratio = (TH1F *)h8_mc->Clone();
  h1ratio->Divide(h4_mc);
     
  h1ratio->SetMarkerColor(kBlack);
  h1ratio->SetMarkerStyle(20);
  h1ratio->SetLineColor(kBlack);
  h1ratio->SetMarkerSize(0.8);

  h1ratio->SetTitle(" ");  
  h1ratio->GetXaxis()->SetTitle("Ehad_mc");
  h1ratio->GetXaxis()->SetLabelFont(42);
  h1ratio->GetXaxis()->SetLabelSize(0.15);
  h1ratio->GetXaxis()->SetTitleSize(0.14);
  h1ratio->GetXaxis()->SetTitleOffset(1);
  h1ratio->GetXaxis()->SetRangeUser(2.5,125.);
  h1ratio->GetXaxis()->SetMoreLogLabels();
  h1ratio->GetXaxis()->SetNoExponent();

  h1ratio->GetYaxis()->SetTitle("Data/MC");
  h1ratio->GetYaxis()->SetTitleSize(0.15);
  h1ratio->GetYaxis()->SetTitleOffset(0.3);
  h1ratio->GetYaxis()->SetLabelFont(42);
  h1ratio->GetYaxis()->SetLabelSize(0.1);
  h1ratio->GetYaxis()->SetRangeUser(0.,2.);
  h1ratio->GetYaxis()->SetNdivisions(5);

  c1_2 = new TPad("c1_2", "newpad",0.008064516,0.0116071,0.9899194,0.2299107);
  c1_2->Draw();
  c1_2->cd();
  
  c1_2->Range(-85.9335,-19.83656,785.9335,21.48034);
  c1_2->SetFillColor(0);
  c1_2->SetBorderMode(0);
  c1_2->SetBorderSize(1);
  c1_2->SetTopMargin(0.03067478);
  c1_2->SetBottomMargin(0.3047036);
  c1_2->SetFrameBorderMode(0);
  c1_2->SetFrameBorderMode(0);

  h1ratio->Draw();
  TLine l1(0.,1.0,55.,1.0);
  l1->SetLineWidth(2);
  l1->SetLineColor(kRed);
  l1->SetLineStyle(1);
  l1->Draw();
  c1->Draw();
  c1->Update();
  c1->SaveAs("Ehad.png");

}
