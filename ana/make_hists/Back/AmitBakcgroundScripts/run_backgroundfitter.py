import ROOT
from ROOT import PlotUtils
from PlotUtils import MnvH2D,MnvH1D,MnvVertErrorBand
import os,sys
#RUN THE FITTER TWO TIMES ONCE WITH THE OPTION OF FLUX WITHOUT FLUX SYSTEMATICS
#ANOTHER WITH RUNNING ONLY IN THE FLUX SYSTEMATICS
#COMBINE THE TWO OUTPUT AND GET THE FINAL OVERALL OUTPUT
#HAD TO DO THIS BECAUSE OF SOME MEMORY LEAK WHICH I CANNOT FIGURE WHERE IT IS COMING FROM....

if len(sys.argv)<3:
    print "python run_backgroundfitter.py recoil_file.root template_file.root output_file.root"
    sys.exit()

rfile = str(sys.argv[1])
tfile = str(sys.argv[2])
ofile = str(sys.argv[3])

o1 = ofile.replace(".root","_1.root")
o2 = ofile.replace(".root","_2.root")
o3 =  ofile.replace(".root","_3.root")
command1 = "python extraFractionFitter_recoilCut_2.py -r "+rfile+" -x -H "+tfile+" -q h_enu_q2_bin_0_mc -P h_pzmu_ptmu_bin_0_mc -o "+o1+" -b -s"

os.system(command1)
#print command1
command2 = "python extraFractionFitter_recoilCut_2.py -r "+rfile+" -x -H "+tfile+" -q h_enu_q2_bin_0_mc -P h_pzmu_ptmu_bin_0_mc -o "+o2+" -b -F -s"
#print command2
os.system(command2)
#sys.exit()

command3 = "python extraFractionFitter_recoilCut_2.py -r "+rfile+" -x -H "+tfile+" -q h_enu_q2_bin_0_mc -P h_pzmu_ptmu_bin_0_mc -o "+o3+" -b -G -s"
#print command2
os.system(command3)
#sys.exit()

ofile1 = ROOT.TFile(o1,"READ")
ofile2 = ROOT.TFile(o2,"READ")
ofile3 = ROOT.TFile(o3,"READ")
final_ofile = ROOT.TFile(ofile,"RECREATE")

keys = ofile1.GetListOfKeys()
for key in keys:
    stuff1 = ofile1.Get(key.GetName())
    stuff2 = ofile2.Get(key.GetName())
    stuff3 = ofile3.Get(key.GetName())
    #print stuff1.ClassName()
    if stuff1.ClassName()=="PlotUtils::MnvH1D" or stuff1.ClassName()=="PlotUtils::MnvH2D":
    	#print "Hey you!!"
    	stuff1.PopVertErrorBand("Flux")
	stuff2.TransferVertErrorBand(stuff1,"Flux")
	vnames = stuff3.GetVertErrorBandNames()
	for name in vnames:
	    if "GENIE" in name:
	        stuff1.PopVertErrorBand(name)
	        stuff3.TransferVertErrorBand(stuff1,name)
        stuff1.Write()
    elif stuff1.ClassName()=="TVector2":
        stuff1.Write("pot")
    else:
        stuff1.Write()

final_ofile.Close()
        
