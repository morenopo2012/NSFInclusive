#include "PlotUtils/MnvH2D.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TPad.h"
#include "TStyle.h"

int main() {
  // Open your ROOT file
  TFile *f = new TFile("XSec_thesis_nuclearflux_minervame1L_t99_z99.root", "READ");
  if (!f || f->IsZombie()) {
    std::cerr << "Could not open file!" << std::endl;
    return 0;
  }

  // Get the MnvH2D histogram
  PlotUtils::MnvH2D* hist2D = (PlotUtils::MnvH2D*)f->Get("simulatedCrossSection_minervame1L");
  if (!hist2D) {
    std::cerr << "Could not find histogram!" << std::endl;
    return 0;
  }

  // Find out how many bins in q3 (Y axis)
  int nbinsY = hist2D->GetYaxis()->GetNbins();
  
  // Create a canvas
  int nx = 4;  // Number of columns
  int ny = (nbinsY+nx-1)/nx;  // Enough rows
  TCanvas* c = new TCanvas("c", "Slices", 3000, 2000); // Bigger size
  c->Divide(nx, ny);

  // Loop over each q3 bin
  for (int iy = 1; iy <= nbinsY; ++iy) { // bin 0 is underflow
    c->cd(iy);

    // Project Eavailable (X axis) for this fixed q3 bin
    TH1D* hProj = hist2D->ProjectionX(Form("proj_q3bin%d", iy), iy, iy);

    hProj->SetLineWidth(2);
    hProj->SetLineColor(kBlue+1);
    hProj->SetTitle(Form("Slice at q3 bin %d", iy));
    hProj->GetXaxis()->SetTitle("Eavailable (GeV)");
    hProj->GetYaxis()->SetTitle("SimulatedCrossSection");

    hProj->Draw("hist");
  }

  // Save the panel
  c->SaveAs("Slices_Eavailable_per_q3bin.png");
  c->SaveAs("Slices_Eavailable_per_q3bin.pdf");

  return 0;
}

