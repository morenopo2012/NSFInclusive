//#include "myPlotStyle.h"
#include "GridCanvas.h"
#include "PlotUtils/MnvH2D.h"
// #include "../util/plot/plot.h"
#include <fstream>

#include "TCanvas.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TStopwatch.h"
#include "TEnv.h"
#include "TChain.h"
#include "TF2.h"
#include "Math/DistFunc.h"
#include "TGraphErrors.h"
#include "TStyle.h"
#include "TFile.h"
#include "localColor.h"
#include "PlotUtils/MnvColors.h"
//#include "Cintex/Cintex.h"

#include "myPlotStyle.h"
#include <THStack.h>

#include "plot.h"
#include <algorithm>
using namespace PlotUtils;
using namespace std;
void makePlots()
{
 // ROOT::Cintex::Cintex::Enable();
 myPlotStyle();
 TH1::AddDirectory(false);
 TH1::SetDefaultSumw2();
 gStyle->SetErrorX(0);
 gStyle->SetEndErrorSize(2);
 gStyle->SetLabelSize(0.04);

 std::cout<<"Good"<<std::endl;
 // Open the ROOT file
 TFile* file = TFile::Open("../KE_C_Carbon_ResonantPion_FullStatistics_Weighted.root", "READ");// KE_PB_Lead_QE_FullStatistics.root //KE_PB_QE_Weighted.root KE_PB_ResonantOnlyPions_Weighted.root KE_PB_Inclusive_Weighted.root
 if (!file || file->IsZombie()) {
    std::cerr << "Error opening file KE_QE.root\n";
    return;
 }

 // Retrieve histograms
 const char* hist_names[] = {
    "hist_fate_NC", "hist_fate_Elas", "hist_fate_ChEx", "hist_fate_KnockOut", 
    "hist_fate_MultN", "hist_fate_Abs",  "hist_fate_nn_50",
    "hist_fate_pn_51", "hist_fate_pp_52" ,"hist_fate_PionP"
 };
 const char* legend_entries[] = {"No FSI", "Elas", "ChargEx", "Inelastic", "MultN", "Abs", "nn","pn","pp" , "PionP"};
 const int num_hists = 10;

 TH1* histograms[num_hists];
 for (int i = 0; i < num_hists; ++i) {
    histograms[i] = (TH1*)file->Get(hist_names[i]);
    if (!histograms[i]) {
        std::cerr << "Error: Histogram " << hist_names[i] << " not found in file.\n";
        return;
    }
 }


// Create output file for the table
    std::ofstream tableFile("FateTable.txt");
    if (!tableFile.is_open()) {
        std::cerr << "Error: Unable to open output file.\n";
        return;
    }

    // Write header
    tableFile << "Bin\tTotal Entries\t";
    for (int i = 0; i < num_hists; ++i) {
        tableFile << legend_entries[i] << "\t" << legend_entries[i] << " frac\t";
    }
    tableFile << "\n";

    // Loop over bins and calculate totals and fractions
    int nBins = histograms[0]->GetNbinsX();
    for (int bin = 1; bin <= nBins; ++bin) {
        double totalEntries = 0.0;

        // Calculate total entries for the bin
        for (int i = 0; i < num_hists; ++i) {
            totalEntries += histograms[i]->GetBinContent(bin);
        }

    // Write bin number and total entries
        tableFile << bin << "\t" << totalEntries << "\t";

        // Write entries and fractions for each fate
        for (int i = 0; i < num_hists; ++i) {
            double binContent = histograms[i]->GetBinContent(bin);
            double fraction = (totalEntries > 0) ? binContent / totalEntries : 0.0;
            tableFile << binContent << "\t" << fraction << "\t";
        }
        tableFile << "\n";
    }

    // Close the file
    tableFile.close();
    std::cout << "Table saved to FateTable.txt\n";





 // Create a stacked histogram
 THStack* hs = new THStack("hs", "Stacked Histograms;Hadron KE (MeV);Events/ 40 MeV");

 // Define colors for each histogram
// int colors[] = {kRed, kRed-3, kGreen-3 , kGreen ,kMagenta-3, kMagenta, kMagenta+1, kMagenta+2, kMagenta+3 , kCyan, }; //kOrange, kYellow+1
 int colors[] = {TColor::GetColor("#A2C7E8"), TColor::GetColor("#8FBDE5"),TColor::GetColor("#A2E8C7"),TColor::GetColor("#8FE5BD"),TColor::GetColor("#D1C0B0"),TColor::GetColor("#C849A9"),TColor::GetColor("#CABBE9"),TColor::GetColor("#D4A5E6"),TColor::GetColor("#E0B7E4"),TColor::GetColor("#F4A582") };

 for (int i = 0; i < num_hists; ++i) {
     histograms[i]->SetFillColor(colors[i]);
     histograms[i]->SetLineColor(kBlack);
     hs->Add(histograms[i]);
 }
 
 // Create a canvas and draw the stack
 TCanvas* c = new TCanvas("c", "Stacked Histograms", 800, 600);
 hs->Draw("HIST");

 // Add a legend
 TLegend* legend = new TLegend(0.7, 0.7, 0.9, 0.9);
 legend->AddEntry(histograms[9], legend_entries[9], "f"); //Pion
 legend->AddEntry(histograms[8], legend_entries[8], "f"); //pp
 legend->AddEntry(histograms[7], legend_entries[7], "f"); //pn
 legend->AddEntry(histograms[6], legend_entries[6], "f"); //nn
 legend->AddEntry(histograms[5], legend_entries[5], "f"); //Abs
 //legend->AddEntry(histograms[4], legend_entries[4], "f"); //MultNucl
 legend->AddEntry(histograms[3], legend_entries[3], "f"); //Inelastic
 legend->AddEntry(histograms[2], legend_entries[2], "f"); //ChEx
 legend->AddEntry(histograms[1], legend_entries[1], "f"); //Elast
 legend->AddEntry(histograms[0], legend_entries[0], "f"); //No Scat
 //for (int i = 0; i < num_hists; ++i) {
 //    legend->AddEntry(histograms[i], legend_entries[i], "f");
 //}
 legend->Draw();

 // Save the plot
 c->SaveAs("KE.png");
 c->SaveAs("KE.C");

 // Clean up
 file->Close();
 delete file;

}

int main(int argc, char* argv[])
{
  makePlots();
  return 0;
}
