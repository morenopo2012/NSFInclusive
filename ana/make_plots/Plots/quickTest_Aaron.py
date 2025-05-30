import ROOT,PlotUtils
from plottingClasses import *
from PlotUtils import MnvH1D,MnvH2D,MnvPlotter
#Load and implement Phil's plot style header file
#ROOT.gROOT.ProcessLine(".L ../style/myPlotStyle.h")
#ROOT.myPlotStyle()

plotter = PlotUtils.MnvPlotter()

filePath = 'Hists_EventSelection_t14_z82_Nu_v1_.root'
#filePath = '/pnfs/minerva/persistent/users/zdar/NukeHists/v1/Hists_Energy2D_v1_t1_z26_Nu_v1_.root'
histFile = ROOT.TFile(filePath)
# Zubair
#hist = histFile.Get("Enu_dis")
#hist = histFile.Get("Enu_mc")
hist = histFile.Get("selected_mc_reco_Emu")

ebs = hist.GetVertErrorBandNames()
for eb in ebs:
  print eb

branch = "Flux"
#branch = "BeamAngleX"
#branch = "genie_NormCCRES"
#branch = "genie_MFP_pi"
#branch = "Genie_InteractionModel"
#branch = "Genie_FSI"

EB = hist.GetVertErrorBand(branch)
x = EB.GetNHists()

print 'nHists in {0} error band: {1}'.format(EB,x)

with makeEnv_TCanvas('aaronTest.png'):
  plotter.DrawErrorSummary(hist)#,"TR",True,True,0.00001,False,branch)
  input("Press enter to end program")

#exec("hist.GetVertErrorBand(\"{0}\").DrawAll()".format(branch))

input("Press enter to end program")

