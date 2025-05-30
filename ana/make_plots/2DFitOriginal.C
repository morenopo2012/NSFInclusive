#include "TH2D.h"
#include "PlotUtils/MnvH2D.h"
#include "TF2.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TRandom3.h"
#include "TVirtualFitter.h"
#include "TList.h"
#include "TMath.h"
#include "TFile.h"  
#include <iostream>



double gauss2D(double *x, double *par) {
   double angle = 180*par[5]/TMath::Pi();
   double X = double(x[0]*cos(angle)+x[1]*sin(angle));
   double Y = double(-x[0]*sin(angle)+x[1]*cos(angle));
   
   double z1 = double((X-par[1])/par[2]);
   double z2 = double((Y-par[3])/par[4]);

   return - par[0]*exp(-0.5*(z1*z1+z2*z2));
}

/*double my2Dfunc(double *x, double *par) {
   return gauss2D(x,&par[0]) + 0.1* gauss2D(x,&par[4]);
}*/

// data need to be globals to be visible by fcn
    TRandom3 rndm;
    TH2D *h1, *h2;
    Int_t npfits;
 
//wrapping   
void myFcn(Int_t &npar, Double_t *gin, Double_t &fval, Double_t *p, Int_t iflag) {      
        
  double chi2 = 0;
  double x[2];
  double tmp;
  npfits = 0;
  
  TAxis *ax = h1->GetXaxis();
  TAxis *ay = h1->GetYaxis();  
  int nbinX1 = h1->GetNbinsX();
  int nbinY1 = h1->GetNbinsY();  
    
  for (int ix = 1; ix <= nbinX1; ++ix) {     
    for (int iy = 1; iy <= nbinY1; ++iy) {
         x[0] = ax->GetBinCenter(ix);
          
        if ( h1->GetBinError(ix,iy) > 0 ) {
             x[1] = ay->GetBinCenter(iy);
             tmp = (h1->GetBinContent(ix,iy) - gauss2D(x,p))/h1->GetBinError(ix,iy);
            
            chi2 += tmp*tmp;
            npfits++;
      }
    }
  }
  fval = chi2;
}
  
  double *p;
  const double mx1 = p[1];		//moyenne en x
  const double my1 = p[3];		//moyenne en y
  const double sx1 = p[2];		//sigma en x
  const double sy1 = p[4];		//sigma en y
  double x, y;


int fit2dHist(int option=1) {
  // create two histograms
  int nbx1 = 10;
  int nby1 = 8;
  double xlow1 = 0.;
  double ylow1 = 0.;
  double xup1 = 50.;
  double yup1 = 25.;
  
  double iniParams[6] = {100, 100., 6., 40., 10., 55.};
  
  
  // create fit function
  TF2 * func = new TF2("func",gauss2D,xlow1,xup1,ylow1,yup1, 6);
  func->SetParameters(iniParams);
   cout<<"Hi I am here"<<endl; 
  //Histo
   TFile * f = new TFile("test.root");
   TCanvas *c = (TCanvas *)f->Get("selected_data2d_reco_Emu_Ehad");
   h1 = (TH2D *)c->GetPrimitive("selected_data2d_reco_Emu_Ehad");
   
   cout<<"Hi I am here"<<endl; 

 
  bool global = false;
  if (option > 10) global = true;
  if (global) {
    // fill data structure for fit (coordinates + values + errors)
    std::cout << "Do global fit" << std::endl;
    
    // fit now all the function together
   TVirtualFitter::SetDefaultFitter("Minuit");
   TVirtualFitter * minuit = TVirtualFitter::Fitter(0,6);
    
    for (int i = 0; i < 6; ++i) {
      minuit->SetParameter(i, func->GetParName(i), func->GetParameter(i), 0.01, 0,0);
    }
    
    minuit->SetFCN(myFcn);
    double arglist[100];
    arglist[0] = 0;
    // set print level
    minuit->ExecuteCommand("SET PRINT",arglist,2);
    
    // minimize
    arglist[0] = 5000; // number of function calls
    arglist[1] = 0.01; // tolerance
    minuit->ExecuteCommand("MIGRAD",arglist,2);
    
    //get result
    double minParams[6];
    double parErrors[6];
    for (int i = 0; i < 6; ++i) {
      minParams[i] = minuit->GetParameter(i);
      parErrors[i] = minuit->GetParError(i);
    }
    
    double chi2, edm, errdef;
    int nvpar, nparx;
    minuit->GetStats(chi2,edm,errdef,nvpar,nparx);
    func->SetParameters(minParams);
    func->SetParErrors(parErrors);
    func->SetChisquare(chi2);
    int ndf = npfits-nvpar;
    func->SetNDF(ndf);
    
    // add to list of functions
    h1->GetListOfFunctions()->Add(func);
  }
  
  else {
    // fit independently
    h1->Fit(func);
  }
  
  
  // Create a new canvas.
  TCanvas * c1 = new TCanvas("c1","Mass and Width Fit",100,10,800,500);
  c1->Divide(2,1);
  gStyle->SetOptFit();
  gStyle->SetStatY(0.5);
  c1->cd(1);
  h1->Draw();
  func->SetRange(xlow1,ylow1,xup1,yup1);
  func->DrawCopy("cont1 same");
  c1->cd(2);
  h1->Draw("lego");
  func->DrawCopy("surf1 same");
  return 0;
}
