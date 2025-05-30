//******************************data MC***************************

{

  //TFile *f1 = TFile::Open("$HISTS/v1/Hists_Energy_t1_z26_Nu_v1_woGenie.root");
  TFile *f1 = TFile::Open("$HISTS/v1/Hists_Energy_t1_z26_Nu_v1_.root");
  TFile *f2 = TFile::Open("$HISTS/v1/Hists_Energy_t1_z26_Nu_v1_.root");

  //**********************************************Ratio************************
  TCanvas *c1 =  new TCanvas("c1","",323,144,534,523);
  c1->Draw();
  c1->cd();
  
  c1_2 = new TPad("c1_2","newpad",0.008032129,0.1866667,0.9879518,0.9911111);
  c1_2->SetBottomMargin(0.1);
  c1_2->Draw();
  c1_2->cd();
  c1_2->SetTickx(1);
  c1_2->SetTicky(1);

  TH1F *h1_mc = (TH1F*)f1->Get("Ehad_dis");
  TH1F *h2_mc = (TH1F*)f2->Get("Ehad_dis");


  h1_mc->SetMarkerColor(1);
  h1_mc->SetMarkerStyle(20);
  h1_mc->SetMarkerSize(0.8); 
  h1_mc->SetLineColor(1);

  h2_mc->SetMarkerColor(1);
  h2_mc->SetMarkerStyle(20);
  h2_mc->SetMarkerSize(0.8); 
  h2_mc->SetLineColor(1);
  
  //h2_mc->SetLineColor(kRed-2);
  //h2_mc->SetFillColor(kBlue-2);

  h1_mc->SetStats(0);
  h2_mc->SetStats(0);

  h1_mc->Scale(1e-3, "width");
  h2_mc->Scale(1e-3, "width");
  h1_mc->Scale(0.449);
  h1_mc->GetYaxis()->SetRangeUser(0.,3);
  h1_mc->GetYaxis()->SetTitle("Events/10^3");
  h1_mc->SetTitle("Tracker"); 
  h1_mc->Draw("hist");
  h2_mc->Draw("SAME");


  TLegend *leg1 = new TLegend(0.6905444,0.5759494,0.8825215,0.8417722,NULL,"brNDC");
  leg1->AddEntry(h1_mc,"WoGenie","lp");
  leg1->AddEntry(h2_mc,"WGenie","f");
  leg1->SetFillColor(0);
  leg1->SetLineColor(0);
  leg1->SetTextSize(0.02888889);
  leg1->Draw();

  c1->cd();


  TH1F *h1ratio = (TH1F *)h1_mc->Clone();
  h1ratio->Divide(h2_mc);
     
  h1ratio->SetMarkerColor(kBlack);
  h1ratio->SetMarkerStyle(20);
  h1ratio->SetLineColor(kBlack);
  h1ratio->SetMarkerSize(0.8);

  h1ratio->SetTitle("  ");  
  h1ratio->GetXaxis()->SetTitle("Ehad_dis");
  h1ratio->GetXaxis()->SetLabelFont(42);
  h1ratio->GetXaxis()->SetLabelSize(0.15);
  h1ratio->GetXaxis()->SetTitleSize(0.14);
  h1ratio->GetXaxis()->SetTitleOffset(1);
  h1ratio->GetXaxis()->SetRangeUser(15.,3000.);
  h1ratio->GetXaxis()->SetMoreLogLabels();
  h1ratio->GetXaxis()->SetNoExponent();

  h1ratio->GetYaxis()->SetTitle("woG/wG");
  h1ratio->GetYaxis()->SetTitleSize(0.15);
  h1ratio->GetYaxis()->SetTitleOffset(0.3);
  h1ratio->GetYaxis()->SetLabelFont(42);
  h1ratio->GetYaxis()->SetLabelSize(0.1);
  h1ratio->GetYaxis()->SetRangeUser(0.5,1.5);
  h1ratio->GetYaxis()->SetNdivisions(5);

  c1_1 = new TPad("c1_1", "newpad",0.008064516,0.0116071,0.9899194,0.2299107);
  c1_1->Draw();
  c1_1->cd();
  c1_1->SetLogx();
  c1_1->SetGridx();
  c1_1->SetGridy();
  
  c1_1->Range(-85.9335,-19.83656,785.9335,21.48034);
  c1_1->SetFillColor(0);
  c1_1->SetBorderMode(0);
  c1_1->SetBorderSize(1);
  c1_1->SetTopMargin(0.03067478);
  c1_1->SetBottomMargin(0.3047036);
  c1_1->SetFrameBorderMode(0);
  c1_1->SetFrameBorderMode(0);

  h1ratio->Draw();
  TLine l1(15.,1.0,3000.,1.0);
  l1->SetLineWidth(1);
  l1->SetLineColor(kRed);
  l1->SetLineStyle(5);
  l1->Draw("same");
  c1->Draw();
  c1->Update();

  c1->SaveAs("Ehad_dis_ratio.png");

}
