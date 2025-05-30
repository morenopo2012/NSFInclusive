//******************************data MC***************************

{
  TFile *f1 = TFile::Open("Hists_EventSelection_t1_z26_Nu_v1_.root");

  //**********************************************Z Mass************************
  TCanvas *c1 =  new TCanvas("c1"," ",323,144,534,523);
  //c1->SetTickx(2);
  c1->SetTicky(10);
  c1->SetGridx();
  c1->SetGridy();
  c1->Draw();
  c1->cd();
  
  TH1F *h1_mc = (TH1F*)f1->Get("selected_mc_reco_Emu");
  

  h1_mc->SetLineStyle(1);
  h1_mc->SetLineColor(1);

  h1_mc->Scale(1e-3, "width");
  //h1_mc->Scale(0.45);
  h1_mc->Scale(0.449);
  h1_mc->GetYaxis()->SetRangeUser(0.,20.);
  h1_mc->GetYaxis()->SetTitle("Events/10^3");
  h1_mc->GetXaxis()->SetTitle("Enu_dis");
  h1_mc->SetTitle("Tracker"); 
  h1_mc->Draw("hist");
  

  TLegend *leg = new TLegend(0.4661654,0.5955734,0.8195489,0.861167,NULL,"brNDC");
  leg->AddEntry(h1_mc,"MC", "l");
  leg->AddEntry("" ,"POT Normalized", "l");
  leg->SetBorderSize(1);
  leg->SetTextColor(4);
  leg->SetTextFont(62);
  leg->SetTextSize(0.03888889);
  leg->SetLineColor(1);
  leg->SetLineStyle(1);
  leg->SetLineWidth(2);
  leg->SetFillColor(0);
  leg->SetFillStyle(1001);
  leg->Draw();

  c1->cd();
 
  c1->SaveAs("Emu_dis_ratio.png");
}
