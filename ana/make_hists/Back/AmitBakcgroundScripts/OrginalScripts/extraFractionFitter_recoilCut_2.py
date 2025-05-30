###################################################
#!/usr/bin/env python

from optparse import OptionParser
import sys, os, subprocess, datetime
import gc


usage = "usage: %prog [options]"
parser = OptionParser(usage=usage)

###################################################
##
## Helper Functions
##
###################################################


###################################################
#
# Setup parser that reads in options before importing ROOT
#Otherwise --help reads the some weird pyROOT python
###################################################

parser.add_option("-r","--recoil_file",dest="recoil_file",
                  help= "recoil file for bkg subtraction",default="NA")

parser.add_option("-o","--output_file",dest="output_file",
                  help="name of the output file for bkg subtraction",default="NA")

parser.add_option("-s","--run_systematics",dest="run_systematics",action="store_true",
                  help="Do the fit over all the systematic universes",default=False)

parser.add_option("-c","--cheryl_method",dest="cheryl_method",action="store_true",
                  help="Use Cheryl Method for Fit",default=False)

parser.add_option("-p","--partial_fitrange",dest="partial_fitrange",action="store_true",
                  help="Do the partial fitting over recoil energy range",default=False)

parser.add_option("-b","--rebin",dest="rebin",action="store_true",
                  help="Rebin the Recoil Energy Bins",default=False)

parser.add_option("-d","--e_low",dest="e_low",
                  help="The starting fit bin of Recoil Energy",default=0.1)

parser.add_option("-u","--e_high",dest="e_high",
                  help="The ending fit bin of the Recoil Energy",default=0.5)

parser.add_option("-f","--fit_bin",dest="fit_bin",
                  help = "Particular Super PTPZ bin over which you want to run the fit instead of over all fitting bins (Use for debugging purpose can go from 0 to 13)",default=-1)

parser.add_option("-x","--constrain_par",dest="constrain_par",action="store_true",
                  help = "Constrain or Unconstrain parameter",default=False)

parser.add_option("-i","--initial_frac",dest="initial_frac",action="store_true",
                  help = "Setting True means initial fraction will default to 0.5 for fitting",default=False)

parser.add_option("-D","--debug_mode",dest="debug_mode",action="store_true",
                  help="Setting True means the code stops at multiple places, creates plots etc for debugging",default=False)

parser.add_option("-H","--hist_file",dest="hist_file",
                  help="Create MnvH2D's of scale factor to make bkg subtraction easier",default="NA")

parser.add_option("-q","--q2_enu_hist",dest="q2_enu_hist",
                  help="Template Q2 Enu to fill the Scale Factors",default="NA")                 

parser.add_option("-P","--pzpt_hist",dest="pzpt_hist",
                  help="Template PZPT to fill the Scale Factors",default="NA")

parser.add_option("-V","--verbose_mode",dest="verbose_mode",action="store_true",
                  help= "Set True for Verbose or Runs Quietly ",default=False)

parser.add_option("-F","--flux",dest="flux",action="store_true",
                  help= "Run on Flux only ",default=False)
		  
		  
parser.add_option("-S","--start_univ",dest="start_univ",
		  help="Start Number of Flux Universe ",default=0)
		  
parser.add_option("-E","--end_univ",dest="end_univ",
		  help="Start Number of Flux Universe ",default=100)
		  
parser.add_option("-G","--genie1",dest="genie1",action="store_true",
                  help= "Run on GENIE only ",default=False)
		  
parser.add_option("-Z","--genie2",dest="genie2",action="store_true",
                  help= "Run on GENIE only ",default=False)
(options, args) = parser.parse_args()

###################################################
##
## NOW IMPORT ROOT FUNCTIONS
##
###################################################

from ROOT import *
from PlotUtils import *
import array

###################################################
##
## HELPER FUNCTIONS
##
###################################################


genie_list1=['GENIE_AGKYxF1pi', 'GENIE_BhtBY', 'GENIE_CV1uBY', 'GENIE_EtaNCEL', 'GENIE_FrAbs_pi', 'GENIE_FrCEx_pi', 'GENIE_FrElas_pi', 'GENIE_FrPiProd_N', 'GENIE_MFP_N', 'GENIE_MaCCQE', 'GENIE_MaRES', 'GENIE_NormDISCC', 'GENIE_RDecBR1gamma', 'GENIE_Rvn2pi', 'GENIE_Rvp2pi', 'GENIE_VecFFCCQEshape']


genie_list2 = ['GENIE_AhtBY', 'GENIE_CCQEPauliSupViaKF', 'GENIE_CV2uBY', 'GENIE_FrAbs_N', 'GENIE_FrCEx_N', 'GENIE_FrElas_N', 'GENIE_FrInel_N', 'GENIE_FrPiProd_pi', 'GENIE_MFP_pi', 'GENIE_MaNCEL', 'GENIE_MvRES', 'GENIE_NormNCRES', 'GENIE_Rvn1pi', 'GENIE_Rvp1pi', 'GENIE_Theta_Delta2Npi']



def CreateHistDict(f,i):
    newname = "ptpz_bin_"+str(i)+"_"
    pass_name = "passed_bin_"+str(i)+"_"
    fail_name = "failed_bin_"+str(i)+"_"
    
    h_passqelike = f.Get(pass_name+"qelike")
    h_passqelikenot = f.Get(pass_name+"qelikenot")
    h_passdata = f.Get(pass_name+"data")
    h_qelike = f.Get(newname+"qelike")
    h_qelikenot = f.Get(newname+"qelikenot")
    h_data = f.Get(newname+"data")
    pot = f.Get("pot")
    hist_dict = {"qelike":h_qelike,
                 "qelikenot":h_qelikenot,
                 "data":h_data,
                 "qelike_pass":h_passqelike,
                 "qelikenot_pass":h_passqelikenot,
                 "data_pass":h_passdata,
                 
                 "pot":pot
    }
    return hist_dict

def ClearHistDict(hist_dict):
    del hist_dict["qelike"]
    del hist_dict["qelikenot"]
    del hist_dict["data"]
    del hist_dict["qelike_pass"]
    del hist_dict["qelikenot_pass"]
    del hist_dict["data_pass"]
    hist_dict.clear()


def CreateScaleFactorHist(_mnvhist,vertnames,name="bla"):
    scale_hist = MnvH1D("scale_hist_"+name,"scale_hist_"+name,15,0,15)
    scale_hist.ClearAllErrorBands()
      
    return scale_hist

def CreateCVScaleFactorWithRootFile(rfile,name="bla"):
    prefix_dummy = "ptpz_bin_0_"
    dummy_mnvh1d = rfile.Get(prefix_dummy+"mc")
    vertnames = dummy_mnvh1d.GetVertErrorBandNames()
    scale_hist = CreateScaleFactorHist(dummy_mnvh1d,vertnames,name)
    return scale_hist


#Create the Dictionary to get the pot Fit informations from each Fit bin
#bkg frac,signal frac (before/after fit)
#scale factor for bkg
#chi2 from the FIT
#NDF from the FIT
#NDF from the FIT
def CreateScaleHistDictONE(f,name):
    scale_hist_sig = CreateCVScaleFactorWithRootFile(f,name)
    scale_hist_bkg = scale_hist_sig.Clone("scale_hist_res")

    frac_sig_before = scale_hist_sig.Clone("frac_sig_before")
    frac_sig_after = scale_hist_sig.Clone("frac_sig_after")
    frac_bkg_before = scale_hist_sig.Clone("frac_bkg_before")
    frac_bkg_after = scale_hist_sig.Clone("frac_bkg_after")

    chi2_hist = scale_hist_sig.Clone("chi2_hist")
    ndf_hist = scale_hist_sig.Clone("ndf_hist")

    scale_hist_dict ={"sig":scale_hist_sig,"bkg":scale_hist_bkg,
                      "frac sig before":frac_sig_before,"frac sig after":frac_sig_after,
                      "frac bkg before":frac_bkg_before,"frac bkg after":frac_bkg_after,
                      "chi2 hist":chi2_hist,
                      "ndf hist":ndf_hist
    }
    
    return scale_hist_dict

def CreateScaleHistDict(o):
    scale_hist_bkg = o.Get("scale_hist_res")
    scale_hist_sig = o.Get("scale_hist_sig")
    frac_sig_before = o.Get("frac_sig_before")
    frac_sig_after =  o.Get("frac_sig_after")
    frac_bkg_before = o.Get("frac_bkg_before")
    frac_bkg_after = o.Get("frac_bkg_after")

    chi2_hist = o.Get("chi2_hist")
    ndf_hist = o.Get("ndf_hist")

    scale_hist_dict ={"sig":scale_hist_sig,"bkg":scale_hist_bkg,
                      "frac sig before":frac_sig_before,"frac sig after":frac_sig_after,
                      "frac bkg before":frac_bkg_before,"frac bkg after":frac_bkg_after,
                      "chi2 hist":chi2_hist,
                      "ndf hist":ndf_hist
                      }
    return scale_hist_dict

def GetTotalBestComponent(nabove,xcom,frac_com_after): #xcom = signal or bkg depending on which you want
    return frac_com_after*nabove/xcom

def extraFitter(indata,insignal,inbackground,indata_pass,insignal_pass,inbackground_pass,options,loval=0,hival=0):
    mc = TObjArray(2)
    mc.Add(insignal)
    mc.Add(inbackground)
    tot_event = insignal.Integral()+inbackground.Integral()
    sig_frac = insignal.Integral()/tot_event
    bkg_frac = inbackground.Integral()/tot_event

    MODE="Q"
    if options.verbose_mode:
        MODE="V"
    fit = TFractionFitter(indata,mc,MODE) #V for VERBOSE RUNNING
    version = gROOT.GetVersion()#VERSION 5 AND 6 HANDLES THINGS DIFFERENTLY
    if "6." in version:
        print "ROOT 6"
        fitter = fit.GetFitter()#GET THE FITTER OBJECT
        fitconfig = fitter.Config()#CONFIG OBJECT TO MANIPULATE THE FIT PARAMETERS.
        if not options.initial_frac:
            fitconfig.ParSettings(0).Set("QELike",sig_frac,0.0001)#PUT THE OPTION FOR INITIAL FRACTION HERE
            fitconfig.ParSettings(1).Set("QELikeNot",bkg_frac,0.0001)#PUT THE OPTION FOR INITIAL FRACTION HERE
    else:
        print "ROOT 5"
        if not options.initial_frac:
            fitter= fit.GetFitter()
            fitter.SetParameter(0,"QELIKE",sig_frac,0.0001,0.0,1.0)
            fitter.SetParameter(1,"QELIKENOT",bkg_frac,0.0001,0.0,1.0)

    #PUT THE OPTION TO CONSTRAIN OR NOT CONSTRAIN FIT HERE
    if options.constrain_par:
        print "YOU HAVE CHOSEN TO CONSTRAIN THE PARAMETER. IT WILL ONLY FIT FROM 0.0001 TO 1.0 RANGE"
        #okay ALSO IF I AM WORKING IN GPVM....IT NEEDS TO CALL VIRTUAL FITTER
        fit.Constrain(0,0.0001,1.0)
        fit.Constrain(1,0.0001,1.0)

    lo = indata.GetXaxis().FindBin(loval)
    hi = indata.GetXaxis().FindBin(hival)

    print "FITTING BIN RANGES ARE ",lo," TO ",hi 

    ndof = 0

    #SET THE FIT RANGE

    fit.SetRangeX(lo,hi)

    status = fit.Fit()

    print "STATUS OF THE FIT IS ",status

    #get various integrals for inside and outside locut area
    lo_signal = insignal.Integral(0,lo-1) #signal outside fit range
    hi_signal = insignal.Integral(lo,hi) #signal  inside fit range
    lo_bkg = inbackground.Integral(0,lo-1) #bkg outside fit range
    hi_bkg = inbackground.Integral(lo,hi) #signal inside fit range
    sig_tot = insignal.Integral(0,hi)
    bkg_tot = inbackground.Integral(0,hi)
    nabove = hi_signal+hi_bkg #total events in fitting region
    ntot = sig_tot+bkg_tot

    
    #Now the best fit fraction from the fit...
    fitf_sig, fite_sig = Double(0), Double()
    fit.GetResult(0,fitf_sig,fite_sig)
    fitf_bkg, fite_bkg = Double(0), Double()
    fit.GetResult(1,fitf_bkg,fite_bkg)

    #now the fractions....
    #AFTER FIT
    frac_sig_after = fitf_sig
    frac_bkg_after = fitf_bkg

    #result = TH1F()
    result = fit.GetPlot().Clone()
    chi2f = fit.GetChisquare()
    if options.debug_mode:
        c = TCanvas()
        c.cd()
        result.Print("ALL")
        indata.Draw("hist")
        result.Draw("histSame")
        raw_input()
        del c
    ndf = fit.GetNDF()

    #BEFORE FIT
    frac_sig_before = hi_signal/nabove
    frac_bkg_before = hi_bkg/nabove

    #FRACTION OF SIGNAL/BKG INSIDE THE FIT REGION
    xS = hi_signal/(lo_signal+hi_signal)
    xB = hi_bkg/(lo_bkg+hi_bkg)
    
    tot_best_sigfrac = GetTotalSignalFraction(xS,xB,frac_sig_after) #FOR THE DATA
    tot_best_bkgfrac = GetTotalBkgFraction(xS,xB,frac_bkg_after) #FOR THE DATA


    
    nabove_best = result.Integral(lo,hi)
    tot_best_bkg = GetTotalBestComponent(nabove_best,xB,frac_bkg_after)#xB
    tot_best_sig = GetTotalBestComponent(nabove_best,xS,frac_sig_after)#xS

    eff_sig = insignal_pass.Integral(0,hi)/sig_tot #EFFICIENCY OF SIGNAL INSIDE RECOIL REGION
    eff_bkg = inbackground_pass.Integral(0,hi)/bkg_tot #EFFICIENCY OF SIGNAL OUTSIDE RECOIL REGION

    print "BEST SIGNAL AND BACKGROUND ",tot_best_sig,tot_best_bkg
    print "EFFICIENCY CORRECTED VAL SIGNAL AND BACKGROUND ",tot_best_sig*eff_sig,tot_best_bkg*eff_bkg
    effcor_best_sigfrac = eff_sig*tot_best_sig/(eff_sig*tot_best_sig+eff_bkg*tot_best_bkg)
    effcor_best_bkgfrac = eff_bkg*tot_best_bkg/(eff_sig*tot_best_sig+eff_bkg*tot_best_bkg)
    print "BEST BKG AND SIGNAL INTEGRALS AND DATA ",tot_best_bkg,tot_best_sig,indata.Integral(lo,hi)
    ntot_best = tot_best_bkg+tot_best_sig

    #NOW WE NEED TO CORRECT FOR THE EFFICIENCY FOR THE SIGNAL AND BACKGROUDN
    #effcor_best_sigfrac = eff_sig*tot_best_sigfrac/(eff_sig*tot_best_sigfrac+eff_bkg*tot_best_bkgfrac)
    #effcor_best_bkgfrac = eff_sig*tot_best_bkgfrac/(eff_sig*tot_best_sigfrac+eff_bkg*tot_best_bkgfrac)

    print "EFFICIENCY SIGNAL AND BACKGROUND ",eff_sig,eff_bkg
    print "EFFICIENCY CORRECTED SIGNAL AND BACKGROUND FRACTION ",effcor_best_sigfrac,effcor_best_bkgfrac
    print "BEST SIGNAL AND BACKGROUND FRACTION ",tot_best_sig/(tot_best_sig+tot_best_bkg),tot_best_bkg/(tot_best_sig+tot_best_bkg)
    
    #sys.exit()
    tot_sigfrac_before = sig_tot/ntot
    tot_bkgfrac_before = bkg_tot/ntot

    print "FIT RANGE SIGNAL,BKG FRAC",fitf_sig,fite_sig
    print "XS and XB ",xS,xB
    print tot_best_sigfrac,tot_best_bkgfrac," BEST SIG AND BKG FRACTION"
    print tot_sigfrac_before,tot_bkgfrac_before," INITIAL SIG AND BKG FRACTION"
    print "SIG CORRECTION FACTOR ",tot_best_sig/(hi_signal+lo_signal)
    print "BKG CORRECTION FACTOR ",tot_best_bkg/(hi_bkg+lo_bkg)

    del fitter
    del mc,fitf_sig,fitf_bkg
    del fit,indata,insignal,inbackground,indata_pass,insignal_pass,inbackground_pass,options
    return {
        "chi2":chi2f,
        "ndof":ndf,
        "frac sig before":tot_sigfrac_before,
        "frac sig after":effcor_best_sigfrac,
        "frac bkg before":tot_bkgfrac_before,
        "frac bkg after":effcor_best_bkgfrac,
        "frac sig after tot":effcor_best_sigfrac,
        "frac bkg after tot":effcor_best_bkgfrac,
        #"corr sig":tot_best_sigfrac/tot_sigfrac_before,
        #"corr bkg":tot_best_bkgfrac/tot_bkgfrac_before,
	"tot sig after":tot_best_sigfrac,
	"tot bkg after":tot_best_bkgfrac,
        "scale sig":tot_best_sig/sig_tot,
        "scale bkg":tot_best_bkg/bkg_tot,
        "err sig":fite_sig,
        "err bkg":fite_bkg,
	"frac sig fit":frac_sig_after,
	"frac bkg fit":frac_bkg_after,
        "best mc":result
        }

def data(f,hist_dict,scale_hist_list,options,bin,isSys=False,sys_name="bla",univ=0):
    
    data = hist_dict["data"]
    #data.Reset()
    qelike = hist_dict["qelike"]
    #data.Add(qelike)
    qelikenot = hist_dict["qelikenot"]
    data_pass = hist_dict["data_pass"]
    qelike_pass = hist_dict["qelike_pass"]
    qelikenot_pass = hist_dict["qelikenot_pass"]
    
    pot = hist_dict["pot"]
    #data.Scale(1.2)
    #data.Add(qelikenot,0.5)
    h_data = TH1D() #DONT FORGET TO DELETE THEM AT THE END OF FUNCTION
    h_qelike = TH1D() #DONT FORGET TO DELETE THEM AT THE END OF FUNCTION
    h_qelikenot = TH1D() #DONT FORGET TO DELETE THEM AT THE END OF FUNCTION
    h_data_pass = TH1D()
    h_qelike_pass = TH1D()
    h_qelikenot_pass = TH1D()
    if not isSys:
        h_data = TH1D(data)
        h_qelike = TH1D(qelike)
        h_qelikenot = TH1D(qelikenot)
        h_data_pass = TH1D(data_pass)
        h_qelike_pass = TH1D(qelike_pass)
        h_qelikenot_pass = TH1D(qelikenot_pass)
        h_data_pass = TH1D(data_pass)
    if isSys:
        h_data = TH1D(data)
        print type(qelike)
        h_qelike = TH1D(qelike.GetVertErrorBand(sys_name).GetHist(univ))
        h_qelikenot =  TH1D(qelikenot.GetVertErrorBand(sys_name).GetHist(univ))
        h_qelike_pass =  TH1D(qelike_pass.GetVertErrorBand(sys_name).GetHist(univ))
        h_qelikenot_pass =  TH1D(qelikenot_pass.GetVertErrorBand(sys_name).GetHist(univ))
        h_data_pass = TH1D(data_pass)

    if options.rebin:
        h_data.Rebin()
        h_qelikenot.Rebin()
        h_qelike.Rebin()
        h_qelike_pass.Rebin()
        h_qelikenot_pass.Rebin()
        h_data_pass.Rebin()
        
    pot_wgt = pot.X()/pot.Y()    
    #h_qelike.Scale(pot_wgt)
    #h_qelikenot.Scale(pot_wgt)
    #h_qelike_pass.Scale(pot_wgt)
    #h_qelikenot_pass.Scale(pot_wgt)
    area_norm =  h_data.Integral()/(h_qelike.Integral()+h_qelikenot.Integral())

    print "AREA NORM FACTOR ",area_norm
    print "POT NORM FACTOR ",pot_wgt
    #fit doesnt need the area normalization....
    #h_qelike.Scale(area_norm)
    #h_qelikenot.Scale(area_norm)
    #h_qelike_pass.Scale(area_norm)
    #h_qelikenot_pass.Scale(area_norm)

    e_high = options.e_high #Recoil Energy in GeV
    e_low = options.e_low  #Recoil Energy in GeV

    if options.cheryl_method:
        e_low = 0.2 #Recoil Energy in GeV
        e_high = 0.5 #Recoil Energy in GeV

    print "FITTING RANGE IS ",e_low," TO ",e_high

    result = extraFitter(h_data,h_qelike,h_qelikenot,h_data_pass,h_qelike_pass,h_qelikenot_pass,options,e_low,e_high)

    print "MEMORY SIZE OF RESULT ",sys.getsizeof(result)
    #NOW GET THE RESULTS...
    scale_hist_sig = scale_hist_list["sig"]
    scale_hist_bkg = scale_hist_list["bkg"]
    
    frac_sig_before = scale_hist_list["frac sig before"]
    frac_sig_after = scale_hist_list["frac sig after"]
    frac_bkg_before = scale_hist_list["frac bkg before"]
    frac_bkg_after = scale_hist_list["frac bkg after"]
    chi2_hist = scale_hist_list["chi2 hist"]
    ndf_hist = scale_hist_list["ndf hist"]
    if isSys==False:
        #I WANT TO CREATE THE BEST MC, BEST SIGNAL AND BEST BACKGROUND AS WELL....
        _data = TH1D(data.Clone())
        _qelike = TH1D(hist_dict["qelike"].Clone())
        _qelikenot = TH1D(hist_dict["qelikenot"].Clone())
        if options.rebin:
            _data.Rebin()
            _qelike.Rebin()
            _qelikenot.Rebin()
        best_qelike = _qelike.Clone("best_"+h_qelike.GetName())
        best_qelikenot = _qelikenot.Clone("best_"+h_qelikenot.GetName())

        #best_qelike.Reset()
        #best_qelikenot.Reset()

        best_qelike.SetDirectory(0)
        best_qelikenot.SetDirectory(0)
        #SCALE THESE WITH THE CORRECTION FACTOR FROM THE FIT

        #I THINK THE EASIEST WAY TO DO IS CREATE TH1DS FOR THIS MULTIPLICATION AND DO
        #BINOMIAL MULTIPLICATIONS.....

        h_scale_sig = best_qelike.Clone()
        h_scale_bkg = best_qelike.Clone()

        h_scale_sig.Reset()
        h_scale_bkg.Reset()

        
        nom_stack = THStack(newname+"nomstack","")
        h_qelikenot.SetTitle("QELIKENOT")
        h_qelike.SetLineColor(4)
        h_qelike.SetTitle("QELIKE")
        h_qelikenot.SetLineColor(2)
        h_qelikenot.SetLineWidth(2)
        h_qelikenot.SetMarkerSize(0)
        h_qelike.SetLineWidth(3)
        h_qelike.SetMarkerSize(0)
        lo = h_qelike.GetXaxis().FindBin(float(options.e_low))
        hi = h_qelike.GetXaxis().FindBin(float(options.e_high))
        best_mc = result["best mc"].Clone(newname+"bestmc")
        best_tot = best_mc.Integral()
        best_sig = best_tot*result["frac sig fit"]
        best_bkg = best_tot*result["frac bkg fit"]
        print best_sig,best_bkg
        #print lo,hi
        #sys.exit()
        tot_mc = _qelike.Integral(lo,hi)+_qelikenot.Integral(lo,hi)
        best_qelike.Scale(best_sig/_qelike.Integral(lo,hi))
        best_qelikenot.Scale(best_bkg/_qelikenot.Integral(lo,hi))
        
        nom_stack.Add(h_qelikenot)
        nom_stack.Add(h_qelike)
  
        best_stack = THStack(newname+"beststack","")
        best_qelikenot.SetTitle("QELIKENOT")
        best_qelike.SetTitle("QELIKE")
        best_qelike.SetLineColor(4)
        best_qelikenot.SetLineColor(2)
        best_qelikenot.SetLineWidth(2)
        best_qelikenot.SetMarkerSize(0)
        best_qelike.SetLineWidth(3)
        best_qelike.SetMarkerSize(0)

        #best_area_norm = h_data.Integral()/(best_qelike.Integral()+best_qelikenot.Integral())
        #best_qelike.Scale(best_area_norm)
        #best_qelikenot.Scale(best_area_norm)
        #best_qelike.Scale(pot_wgt)
        #best_qelikenot.Scale(pot_wgt)
        h_mc = _qelike.Clone(newname+"mc")
        h_mc.Add(_qelikenot)
        best_stack.Add(best_qelikenot)
        best_stack.Add(best_qelike)
        _data.Write()
        h_mc.Write()
        _qelike.Write()
        _qelikenot.Write()
        best_qelike.Write()
        best_qelikenot.Write()
        best_stack.Write()
        nom_stack.Write()
        
        #NOW GET THE STUFF........
        scale_hist_sig.SetBinContent(bin+1,result["scale sig"])
        scale_hist_sig.SetBinError(bin+1,result["err sig"])
        scale_hist_bkg.SetBinContent(bin+1,result["scale bkg"])
        scale_hist_bkg.SetBinError(bin+1,result["err bkg"])
        frac_sig_before.SetBinContent(bin+1,result["frac sig before"])
        frac_sig_after.SetBinContent(bin+1,result["frac sig after tot"])
        frac_sig_after.SetBinError(bin+1,result["err sig"])
        frac_bkg_before.SetBinContent(bin+1,result["frac bkg before"])
        frac_bkg_after.SetBinContent(bin+1,result["frac bkg after tot"])
        frac_bkg_after.SetBinError(bin+1,result["err bkg"])
        chi2_hist.SetBinContent(bin+1,result["chi2"])
        ndf_hist.SetBinContent(bin+1,result["ndof"])
        del _qelike,_qelikenot,_data
        del best_qelike,best_qelikenot
        del nom_stack,best_stack,best_mc,h_mc
    if isSys==True:
        univ_scale_sig = scale_hist_sig.GetVertErrorBand(sys_name).GetHist(univ)
        univ_scale_bkg  = scale_hist_bkg.GetVertErrorBand(sys_name).GetHist(univ)
        univ_frac_sig_before = frac_sig_before.GetVertErrorBand(sys_name).GetHist(univ)
        univ_frac_sig_after = frac_sig_after.GetVertErrorBand(sys_name).GetHist(univ)
        univ_frac_bkg_before = frac_bkg_before.GetVertErrorBand(sys_name).GetHist(univ)
        univ_frac_bkg_after = frac_bkg_after.GetVertErrorBand(sys_name).GetHist(univ)
        univ_chi2 = chi2_hist.GetVertErrorBand(sys_name).GetHist(univ)
        univ_ndf = ndf_hist.GetVertErrorBand(sys_name).GetHist(univ)
        univ_scale_sig.SetBinContent(bin+1,result["scale sig"])
        univ_scale_bkg.SetBinContent(bin+1,result["scale bkg"])
        univ_frac_sig_before.SetBinContent(bin+1,result["frac sig before"])
        univ_frac_sig_after.SetBinContent(bin+1,result["frac sig after tot"])
        univ_frac_bkg_before.SetBinContent(bin+1,result["frac bkg before"])
        univ_frac_bkg_after.SetBinContent(bin+1,result["frac bkg after tot"])
        univ_chi2.SetBinContent(bin+1,result["chi2"])
        univ_ndf.SetBinContent(bin+1,result["ndof"])
        del univ_scale_sig,univ_scale_bkg,univ_frac_sig_before,univ_frac_sig_after,univ_frac_bkg_before,univ_frac_bkg_after,univ_chi2,univ_ndf
    del h_qelike,h_qelikenot,h_data
    del h_qelike_pass,h_qelikenot_pass,h_data_pass
    del data,qelike,qelikenot,data_pass,qelike_pass,qelikenot_pass
    ClearHistDict(hist_dict)
    del result["best mc"]
    result.clear()
    del result

        
  
  
    
def GetTotalSignalFraction(xS,xB,fSignal):
    #num = (1.0 - (xS-1)/xS)*fSignal
    #den = 1.0+(fSignal-1)/xB*(xB-1) - fSignal*(xS-1)/xS
    num = fSignal*xB
    den = xS-fSignal*xS+fSignal*xB
    tot_frac= num/den
    if tot_frac<0 or tot_frac>1:
        print "Unphysical best signal fraction....exiting",tot_frac
        sys.exit()
    return num/den

def GetTotalBkgFraction(xS,xB,fBkg):
    num  = fBkg*xS
    den = xB-fBkg*xB+fBkg*xS
    tot_frac = num/den
    if tot_frac<0 or tot_frac>1:
        print "Unphysical best background fraction....exiting",tot_frac
        sys.exit()
    return num/den

###################################################
##
## END OF HELPER FUNCTIONS........WRITE STUFF
##BELOW THIS. WRITE FUNCTIONS ABOVE THIS
###################################################



#THE INPUT FILE
if str(options.recoil_file)=="NA":
    print "Need to Give proper recoil_file"
    sys.exit()


q2 = MnvH2D()
pz = MnvH2D()
if str(options.hist_file)!="NA":
    hist_filename = str(options.hist_file)
    if str(options.q2_enu_hist)=="NA":
        print "Need to give the name of Q2 Enu histograms template "
        sys.exit()
    if str(options.pzpt_hist)=="NA":
        print "Need to give the name of PZPT histograms template "
        sys.exit()
    q2_hist = str(options.q2_enu_hist)
    pz_hist = str(options.pzpt_hist)
    hist_file = TFile.Open(hist_filename,"READ")
    q2 = hist_file.Get(q2_hist)
    pz = hist_file.Get(pz_hist)

filename = str(options.recoil_file)
output = str(options.output_file)
f = TFile.Open(filename,"READONLY")

exist_output=False
if os.path.exists(output):
    exist_output=True

print exist_output
o = TFile.Open(output,"RECREATE")

pot = f.Get("pot")



if options.fit_bin>-1:
    scale_hist_list = CreateScaleHistDict(f,"sig")
    bin = str(options.fit_bin)
    newname= "ptpz_bin_"+bin+"_"
    pass_name = "passed_bin_"+str(bin)+"_"
    fail_name = "failed_bin_"+str(bin)+"_"
    h_qelike = f.Get(newname+"qelike")
    h_qelikenot = f.Get(newname+"qelikenot")
    h_data = f.Get(newname+"data")
    h_passqelike = f.Get(pass_name+"qelike")
    h_passqelikenot = f.Get(pass_name+"qelikenot")
    h_failqelike = f.Get(fail_name+"qelike")
    h_failqelikenot = f.Get(fail_name+"qelikenot")
    h_passdata = f.Get(pass_name+"data")
    h_faildata = f.Get(fail_name+"data")
    pot = f.Get("pot")
    hist_dict = {"qelike":h_qelike,
                 "qelikenot":h_qelikenot,
                 "data":h_data,
                 "qelike_pass":h_passqelike,
                 "qelikenot_pass":h_passqelikenot,
                 "data_pass":h_passdata,
                 
                 "pot":pot
    }
    
    data(f,hist_dict,scale_hist_list,options,False)


else:
 
    scale_hist_list = CreateScaleHistDictONE(f,"sig")
    for i in range(0,14):
        newname = "ptpz_bin_"+str(i)+"_"
        pass_name = "passed_bin_"+str(i)+"_"
        fail_name = "failed_bin_"+str(i)+"_"
        
        h_qelike = f.Get(newname+"qelike")
        h_qelikenot = f.Get(newname+"qelikenot")
        h_data = f.Get(newname+"data")
        h_passqelike = f.Get(pass_name+"qelike")
        h_passqelikenot = f.Get(pass_name+"qelikenot")
        h_passdata = f.Get(pass_name+"data")
        pot = f.Get("pot")
        hist_dict = {"qelike":h_qelike,
                     "qelikenot":h_qelikenot,
                     "data":h_data,
                     "qelike_pass":h_passqelike,
                     "qelikenot_pass":h_passqelikenot,
                     "data_pass":h_passdata,

                     "pot":pot
                     }
        """
        if options.rebin:
            hist_dict = {"qelike":h_qelike.Rebin(),
                         "qelikenot":h_qelikenot.Rebin(),
                         "data":h_data.Rebin(),
                         "pot":pot
                         }
        """
        print "************START RUNNING OVER CV OF FIT SUPER BIN***************",i
        data(f,hist_dict,scale_hist_list,options,i,False)
        print "******************END RUNNING OVER CV OF FIT SUPER BIN***********",i
        del h_qelike,h_qelikenot,h_data
        del h_passqelike,h_passqelikenot,h_passdata
        print len(hist_dict)
        #sys.exit()
        #del hist_dict["qelike"]
        #del hist_dict["qelikenot"]
        #del hist_dict["data"]
        #del hist_dict["qelike_pass"]
        #del hist_dict["qelikenot_pass"]
        #del hist_dict["data_pass"]
        
        hist_dict.clear()
        del hist_dict
        gc.collect()
    gROOT.Reset()
    #NOW THE SYSTEMATICS
    counter = 0
    if options.run_systematics:
        for i in range(0,14):
            hist_dict = CreateHistDict(f,i)
            if i==0:
                for stuff in scale_hist_list:
                    scale_hist_list[stuff].AddMissingErrorBandsAndFillWithCV(hist_dict["qelike"])

            vertnames = hist_dict["qelike"].GetVertErrorBandNames()
            ClearHistDict(hist_dict)
            for name in vertnames:
                hist_dict = CreateHistDict(f,i)
                if options.flux:
                    if "Flux" not in name:
                        continue
                elif options.genie1:
		    if name not in genie_list1:
		        continue
		elif options.genie2:
		     if name not in genie_list2:
		         continue
		    #if counter<=500:
		    #    continue

		else:
		    if "Flux" in name:
		        continue
		    if "GENIE" in name:
		        continue
                nuniv = hist_dict["qelike"].GetVertErrorBand(name).GetNHists()
                start_uni = 0
                end_uni = nuniv
                if options.flux:
                    start_uni = int(options.start_univ)
                    end_uni = int(options.end_univ)
                ClearHistDict(hist_dict)
                del hist_dict
                for u in range(start_uni,end_uni):
                    hist_dict = CreateHistDict(f,i)
                    print name
                    print "###START RUNNING OVER SYSTEMATICS ",name," OF UNIVERSE ",u," ###"
                    data(f,hist_dict,scale_hist_list,options,i,True,name,u)
                    print "MEMORY SIZE OF DATA ",sys.getsizeof(data)
                    print "###END RUNNING OVER SYSTEMATICS ",name," OF UNIVERSE ",u," ###"
                    counter +=1
                    print "################DONE WITH ",counter," SYSTEMATICS############# "
                    print "MEMORY SIZE OF HIST_DICT ",sys.getsizeof(hist_dict)
                    del hist_dict
                    gc.collect()
                    print gc.get_count()," GARBAGE COLLECTOR"


#NOW WRITE DOWN THE MEMBERS OF SCALE_HIST_LIST
print "NOW WRITE DOWN THE STUFFS"
if str(options.hist_file)!="NA":
    scale_hist = scale_hist_list["frac sig after"]
    for i in range(0,14):
        print "WRITING FRACTION INFO FOR ",i
        scale = scale_hist.GetBinContent(i+1)
        scale_err = scale_hist.GetBinError(i+1)
        temp_q2 = q2.Clone("h_enu_q2_scale_"+str(i))
        temp_q2.Reset()
        xbins = temp_q2.GetNbinsX()
        ybins = temp_q2.GetNbinsY()
        for x in range(0,xbins+1):
            for y in range(0,ybins+1):
                temp_q2.SetBinContent(x,y,scale)
                temp_q2.SetBinError(x,y,scale_err)
                temp_q2.AddMissingErrorBandsAndFillWithCV(q2)
        if options.run_systematics:
            vertnames = temp_q2.GetVertErrorBandNames()
            for name in vertnames:
                nhist = temp_q2.GetVertErrorBand(name).GetNHists()
                for univ in range(0,nhist):
                    temp_hist = temp_q2.GetVertErrorBand(name).GetHist(univ)
                    temp_scale = scale_hist.GetVertErrorBand(name).GetHist(univ)
                    for x in range(0,xbins+1):
                        for y in range(0,ybins+1):
                            temp_hist.SetBinContent(x,y,temp_scale.GetBinContent(i+1))
        temp_pz = pz.Clone("h_pzmu_ptmu_scale_"+str(i))
        temp_pz.Reset()
        xbins = temp_pz.GetNbinsX()
        ybins = temp_pz.GetNbinsY()
        for x in range(0,xbins+1):
            for y in range(0,ybins+1):
                temp_pz.SetBinContent(x,y,scale)
                temp_pz.SetBinError(x,y,scale_err)
                temp_pz.AddMissingErrorBandsAndFillWithCV(pz)
        if options.run_systematics:
            vertnames = temp_pz.GetVertErrorBandNames()
            for name in vertnames:
                nhist = temp_pz.GetVertErrorBand(name).GetNHists()
                for univ in range(0,nhist):
                    temp_hist = temp_pz.GetVertErrorBand(name).GetHist(univ)
                    temp_scale = scale_hist.GetVertErrorBand(name).GetHist(univ)
                    for x in range(0,xbins+1):
                        for y in range(0,ybins+1):
                            temp_hist.SetBinContent(x,y,temp_scale.GetBinContent(i+1))
        temp_pz.Write()
        temp_q2.Write()
        del temp_pz,temp_q2


    
for stuff in scale_hist_list:
    scale_hist_list[stuff].Write()

    

#AND FINALLY THE POT
pot = f.Get("pot")
pot.Write("pot")

del pz,q2

#o.Close()

print "******************************END OF CODE********************"
