void test_2dfits() {
  TFile *f = new TFile("xy_Tracker.root","READ");
   
  TH2D *h2  = (TH2D*)f->Get("selected_data2d_reco_x_y");
  //TF2 *f2 = new TF2("f2","1.5*([2]*(x*x) + [1]*(y*y)+ [0])",0,50,0,25);
  TF2 *f2 = new TF2("f2","1.5*([2]*(x*x) + [1]*(y*y)+ [0])",0,0.8,0,1);
  //TF2 *f2 = new TF2("f2","1.5*([2]*(x*x) + [1]*(y*y)+ [0])",2,16,0,30);
  f2->SetParameters(2,2);
  TGraph2D *gr = new TGraph2D(h2);
  gr->Fit(f2);
  gr->SetMarkerColor(1);
  std::cout<<" evaluating at 3,4 "<<f2->Eval(3,2)<<std::endl;
  TCanvas *cgr = new TCanvas();
  cgr->cd();
  gr->Draw("colz");
  cgr->Print("fitted_graph.png");
  TCanvas *cfr = new TCanvas();
  cfr->cd(); 
  f2->Draw("colz");
  cfr->Print("fit_function.png");
  
}
