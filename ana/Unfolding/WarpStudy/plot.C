{
TFile *f = new TFile("yWy_TrackerP.root");
TH1F *h = (TH1F*)f->Get("Chi2_Iteration_Dists/m_avg_chi2_modelData_trueData_iter_chi2");
h->Scale(0.166);
h->GetYaxis()->SetTitle("Average #Chi^{2} / NDF");
h->GetYaxis()->SetRange(0,1);
h->Draw();
}
