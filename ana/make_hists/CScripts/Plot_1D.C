//******************************data MC***************************

{
    TFile *f1 = TFile::Open("/pnfs/minerva/persistent/users/zdar/data_ME1A_ML_tuples_v2/grid/minerva/ana/numibeam/v22r1p1/00/00/60/94/MV_00006094_Subruns_0001_MasterAnaDev_AnaData_Tuple_v22r1p1.root");
  
    TFile *f2 = TFile::Open("/pnfs/minerva/persistent/users/vansari/data_ana_minervame1A_TBV/grid/minerva/ana/numibeam/v22r1p1/00/00/60/94/MV_00006094_Subruns_0001_MasterAnaDev_AnaData_Tuple_v22r1p1.root");


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
  TTree *t1=(TTree*)f1->Get("MasterAnaDev");
  TTree *t2=(TTree*)f2->Get("MasterAnaDev");


  TH1F *h1_mc = (TH1F*)t1->Get("eventID");
  TH1F *h1_data = (TH1F*)t2->Get("eventID");
  

  h1_mc->SetLineStyle(1);
  h1_mc->SetLineColor(1);

  //h1_data->SetLineColor(kRed-2);
  h1_mc->SetTitle("Chor"); 
  h1_data->SetLineColor(2);
  h1_data->SetFillColor(kBlue-2);
  h1_mc->SetStats(0);
  h1_data->SetStats(0);
  //h1_data->Scale(1e-3, "width");
  //h1_mc->Scale(1e-3, "width");
  //h1_mc->Scale(0.448);
  //h1_mc->Scale(0.421);
 ///// h1_mc->Scale(0.433);
  //h1_mc->GetXaxis()->SetRangeUser(0.,10.);
  //h1_mc->GetYaxis()->SetRangeUser(0.,200);
  //h1_mc->GetYaxis()->SetRangeUser(0.,15);
  //h1_mc->GetXaxis()->SetTitle("Reconstructed Q^{2}");
  //h1_mc->GetYaxis()->SetTitle("Data/ MnvGENIE");
  //h1_mc->SetTitle("Chor"); 
  //h1_mc->Draw("hist");
  //h1_data->Draw("same");
  h1_data->Divide(h1_mc);
  //h1_data->GetXaxis()->SetTitle("Reconstructed Q^{2}");
  //h1_data->GetYaxis()->SetTitle("Data/ MnvGENIE");
  h1_data->Draw("");
  //TLegend *leg1 = new TLegend(0.6905444,0.5759494,0.7825215,0.8417722,NULL,"brNDC");
  TLegend *leg1 = new TLegend(0.5905444,0.5759494,0.6825215,0.8417722,NULL,"brNDC");
  leg1->AddEntry(h1_mc,"MC", "l");
  leg1->AddEntry(h1_data,"Data");
  leg1->AddEntry("" ,"POT Normalized");
  leg1->AddEntry("" ,"Iron of Target 1");
  leg1->SetFillColor(0);
  leg1->SetLineColor(1);
  leg1->SetTextSize(0.02888889);
  leg1->SetTextColor(4);
  leg1->Draw();

  c1->cd();
  c1->Draw();
  c1->Update();
  c1->SaveAs("Ehad_mc.png");
}
