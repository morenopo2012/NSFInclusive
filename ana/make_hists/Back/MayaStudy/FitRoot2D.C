
#include "TF2.h"
#include "TH2.h"
#include "TMath.h"
Double_t g2(Double_t *x, Double_t *par) {
   Double_t r1 = Double_t((x[0]-par[1])/par[2]);
   Double_t r2 = Double_t((x[1]-par[3])/par[4]);
   return par[0]*TMath::Exp(-0.5*(r1*r1+r2*r2));
}
Double_t fun2(Double_t *x, Double_t *par) {
   Double_t *p1 = &par[0];
   Double_t *p2 = &par[5];
   Double_t *p3 = &par[10];
   Double_t result = g2(x,p1) + g2(x,p2) + g2(x,p3);
   return result;
}
void FitRoot2D() {
   const Int_t npar = 15;
   Double_t f2params[npar] =
      {100,-3,3,-3,3,160,0,0.8,0,0.9,40,4,0.7,4,0.7};
   TF2 *f2 = new TF2("f2",fun2,-10,10,-10,10, npar);
   f2->SetParameters(f2params);
   //Create an histogram and fill it randomly with f2
     TH2F *h2 = new TH2F("h2","from f2",40,-10,10,40,-10,10);
     //TFile *f1 = TFile::Open("Ehad.root");
     //TH2F *h2 = new TH2F("c1");
     Int_t nentries = 100000;
     h2->FillRandom("f2",nentries);
   //Fit h2 with original function f2
     Float_t ratio = 4*nentries/100000;
     f2params[ 0] *= ratio;
     f2params[ 5] *= ratio;
     f2params[10] *= ratio;
     f2->SetParameters(f2params);
     h2->Fit("f2");
     f2->Draw("cont1 same");
     //f2->Draw("colz");
     }
