import ROOT,os

class makeEnv_TCanvas(object):

  def __init__(self,plotName,logy=False):
    self.plotName = plotName
    self.logy = logy

  def __enter__(self):
    self.canvas = ROOT.TCanvas( "canvas" , "canvas" , 10 , 10 , 1000, 750 )
    if self.logy: self.canvas.SetLogy(0)
    return self

  def __exit__(self,*exc):
    # If directory for plotName doesn't exist yet, make it
    plotNameComponents = self.plotName.split('/')
    plotDir = '/'.join(plotNameComponents[:-1]) 
    if not os.path.isdir(plotDir):
      print "Making plot directory {0}".format(plotDir)
      os.system( "mkdir %s" % plotDir )
    self.canvas.Print(self.plotName)
    del self.canvas

class makeEnv_TCanvas_nuMigrationMatrix(makeEnv_TCanvas):

  def __init__(self,plotName,nu,bigBins):
    makeEnv_TCanvas.__init__(self,plotName)
    self.nu = nu
    self.bigBins = bigBins

  def __exit__(self,*exc):
   
    scale = 10 if not self.bigBins else 1
    horizontalLine = ROOT.TLine(0.0,self.nu*scale,5.0*scale,self.nu*scale)
    horizontalLine.SetLineColor(ROOT.kRed)
    horizontalLine.Draw()
    verticalLine = ROOT.TLine(self.nu*scale,0.0,self.nu*scale,5.0*scale)
    verticalLine.SetLineColor(ROOT.kRed)
    verticalLine.Draw()

    makeEnv_TCanvas.__exit__(self,*exc)

def setPlotSpecs_nuMigrationMatrix(MM):
  MM.GetXaxis().SetTitle("Reco #nu")
  MM.GetYaxis().SetTitle("True #nu")

def setPlotSpecs_ENu(hist,reco=False,fullBinning=False):
  if reco: hist.GetXaxis().SetTitle('reconstructed E_{#nu}(GeV)')
  else: hist.GetXaxis().SetTitle('E_{#nu}(GeV)')
  if not fullBinning:
    #hist.GetXaxis().SetRangeUser(0,120)
    hist.GetXaxis().SetRangeUser(2,22)
    #hist.GetXaxis().SetRangeUser(2,50)

def setPlotSpecs_eff(hist,hor='ENu'):
  if hor == 'ENu': setPlotSpecs_ENu(hist)
  hist.GetYaxis().SetTitle('Efficiency, #eta')

def setPlotSpecs_effNumerator(hist,hor='ENu'):
  if hor == 'ENu': setPlotSpecs_ENu(hist)
  hist.GetYaxis().SetTitleSize(0.04)
  hist.GetYaxis().SetTitle('#Reco MC Events x 10^{4}/GeV (#eta Numerator)')

def setPlotSpecs_effDenominator(hist,hor='ENu'):
  if hor == 'ENu': setPlotSpecs_ENu(hist)
  hist.GetYaxis().SetTitleSize(0.04)
  hist.GetYaxis().SetTitle('#True MC Events x 10^{4}/GeV (#eta Denominator)')

def setPlotSpecs_dataRate(hist,hor='ENu'):
  if hor == 'ENu': setPlotSpecs_ENu(hist)
  hist.GetYaxis().SetTitle('#Reco Data Events x 10^{4}/GeV')

def setPlotSpecs_flux(hist,hor='ENu'):
  if hor == 'ENu': setPlotSpecs_ENu(hist)
  hist.GetYaxis().SetTitle('Neutrinos/m^{2}/10^{6} POT/GeV')

def setPlotSpecs_fluxRatio(hist,hor='ENu'):
  if hor == 'ENu': setPlotSpecs_ENu(hist)
  hist.GetYaxis().SetTitle('Flux Ratio')
  hist.GetYaxis().SetRangeUser(0.75,1.25)

def setPlotSpecs_dataMCRatio(hist,hor='ENu'):
  if hor == 'ENu': setPlotSpecs_ENu(hist)
  hist.GetYaxis().SetTitle('data/MC Ratio')

def setPlotSpecs_xSection(hist,hor='ENu'):
  if hor == 'ENu': setPlotSpecs_ENu(hist)
  hist.GetYaxis().SetTitle('#sigma (10^{-38}cm^{2}/nucleon)')

def setPlotSpecs_xSectionPerE(hist,hor='ENu'):
  if hor == 'ENu': setPlotSpecs_ENu(hist)
  hist.GetYaxis().SetTitle('#sigma/E (10^{-38}cm^{2}/nucleon/GeV)')

def setPlotSpecs_2D(hist,hor='ENu'):
  if hor == 'ENu': setPlotSpecs_ENu(hist,reco=True)
  hist.GetYaxis().SetTitle('#nu (GeV)')
  global line
  line = ROOT.TLine(3,0,3,5)
  line.SetLineStyle(2)
  line.SetLineWidth(2)
  line.Draw()
  
  return

def setTitles( hist , title ):
  
  hist.SetXTitle('E_{#nu} (GeV)')
  hist.SetYTitle('#nu_{#mu}s/m^{2}/POT/GeV')
  hist.SetTitle('%s'%title)

def declareLegend( size , n_columns , pos ):

  n_lines = size / n_columns

  if pos == "UL":
    leg = ROOT.TLegend(0.17,0.9-(0.07*n_lines),0.5,0.9)
  if pos == "LL":
    leg = ROOT.TLegend(0.17,0.55-(0.07*n_lines),0.5,0.55)
  if pos == "UR":
    leg = ROOT.TLegend(0.7,0.85-(0.07*n_lines),0.9,0.85)
  if pos == "UR-Flux":
    leg = ROOT.TLegend(0.45,0.85-(0.07*n_lines),0.9,0.85)
  if pos == "UR-Flux2":
    leg = ROOT.TLegend(0.58,0.90-(0.07*n_lines),0.9,0.90)
  if pos == "LR":
    leg = ROOT.TLegend(0.55,0.15,0.95,0.15+(0.07*n_lines))
  leg.SetNColumns(n_columns)
  leg.SetBorderSize(0)
  leg.SetFillStyle(0)
  leg.SetTextSize(0.04)
  leg.SetTextFont(42)

  return leg

def declareChi2TextLine( inChi2 ):

  texbox = ROOT.TLatex( 0.78 , 0.25 , "#splitline{{#chi^{{2}}: {0:.1f}}}{{#chi^{{2}}/dof: {1:.2f}}}".format(inChi2,inChi2/16.) )
  texbox.SetNDC() # set position to be in NDC
  texbox.SetTextAlign(23)
  texbox.SetTextSize(0.03)

  return texbox

def setColors( systematicUniverse , color ):

  exec("flux_minus{0}.SetLineColor({1})".format(systematicUniverse,color))
  exec("flux_minus{0}.SetMarkerColor({1})".format(systematicUniverse,color))
  exec("flux_plus{0}.SetLineColor({1})".format(systematicUniverse,color))
  exec("flux_plus{0}.SetMarkerColor({1})".format(systematicUniverse,color))

def writeHist(hist,outFile):
  #exec('localHist = {0}'.format(histName))
  #localHist.SetName(histName)
  outFile.cd()
  print 'Writing {0} to output file'.format(hist.GetName())
  hist.Write()


