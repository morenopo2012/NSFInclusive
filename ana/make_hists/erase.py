import ROOT

# Open the ROOT files
file_1 = ROOT.TFile("/exp/minerva/app/users/omorenop/cmtuser/git-Mat/MINERvA101/opt/etc/MParamFiles/data/Reweight/outNievesRPAratio-anu12C-20GeV-20170202.root")
file_2 = ROOT.TFile("/exp/minerva/app/users/omorenop/cmtuser/git-Mat/MINERvA101/opt/etc/MParamFiles/data/Reweight/outNievesRPAratio-nu12C-20GeV-20170202.root")

# Access the 2D histograms from each file
hist_1 = file_1.Get("hnonrelratio")
hist_2 = file_2.Get("hnonrelratio")

# Check if the histograms were properly loaded
if not hist_1 or not hist_2:
    print("Error: Could not retrieve histograms.")
    exit()

# Create a clone of the first histogram for the ratio
hist_ratio = hist_1.Clone("hist_ratio")
hist_ratio.SetTitle("Ratio of hnonrelratio ( anu 12C / nu 12C)")

# Perform the division of the two histograms
hist_ratio.Divide(hist_2)

# Create a canvas to draw the histograms
canvas = ROOT.TCanvas("canvas", "Canvas", 800, 600)

# Draw the ratio histogram
hist_ratio.Draw("COLZ")

# Save the ratio histogram to a file (optional)
output_file = ROOT.TFile("ratio_histogram.root", "RECREATE")
hist_ratio.Write()
output_file.Close()

# Keep the canvas open
canvas.SaveAs("ratio_histogram_carb.png")  # Save as PNG file
canvas.Update()
canvas.Draw()

