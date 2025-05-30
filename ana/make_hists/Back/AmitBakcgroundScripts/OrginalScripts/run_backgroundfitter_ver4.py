import ROOT
import os,sys
from ROOT import PlotUtils
#from PlotUtils import MnvH2D,MnvH1D,MnvVertErrorBand

#RUN THE FITTER TWO TIMES ONCE WITH THE OPTION OF FLUX WITHOUT FLUX SYSTEMATICS
#ANOTHER WITH RUNNING ONLY IN THE FLUX SYSTEMATICS
#COMBINE THE TWO OUTPUT AND GET THE FINAL OVERALL OUTPUT
#HAD TO DO THIS BECAUSE OF SOME MEMORY LEAK WHICH I CANNOT FIGURE WHERE IT IS COMING FROM....


genie_list1=['GENIE_AGKYxF1pi', 'GENIE_BhtBY', 'GENIE_CV1uBY', 'GENIE_EtaNCEL', 'GENIE_FrAbs_pi', 'GENIE_FrCEx_pi', 'GENIE_FrElas_pi', 'GENIE_FrPiProd_N', 'GENIE_MFP_N', 'GENIE_MaCCQE', 'GENIE_MaRES', 'GENIE_NormDISCC', 'GENIE_RDecBR1gamma', 'GENIE_Rvn2pi', 'GENIE_Rvp2pi', 'GENIE_VecFFCCQEshape']


genie_list2 = ['GENIE_AhtBY', 'GENIE_CCQEPauliSupViaKF', 'GENIE_CV2uBY', 'GENIE_FrAbs_N', 'GENIE_FrCEx_N', 'GENIE_FrElas_N', 'GENIE_FrInel_N', 'GENIE_FrPiProd_pi', 'GENIE_MFP_pi', 'GENIE_MaNCEL', 'GENIE_MvRES', 'GENIE_NormNCRES', 'GENIE_Rvn1pi', 'GENIE_Rvp1pi', 'GENIE_Theta_Delta2Npi']



if len(sys.argv)<3:
    print "python run_backgroundfitter.py recoil_file.root template_file.root output_file.root"
    sys.exit()

rfile = str(sys.argv[1])
tfile = str(sys.argv[2])
ofile = str(sys.argv[3])

o1 = ofile.replace(".root","_1.root")
o2 = ofile.replace(".root","_2.root")
o2_1 = ofile.replace(".root","_2_1.root")
o2_2 = ofile.replace(".root","_2_2.root")
o2_3 = ofile.replace(".root","_2_3.root")
o2_4 = ofile.replace(".root","_2_4.root")
o2_5 = ofile.replace(".root","_2_5.root")
o2_6 = ofile.replace(".root","_2_6.root")
o2_7 = ofile.replace(".root","_2_7.root")
o2_8 = ofile.replace(".root","_2_8.root")
o2_9 = ofile.replace(".root","_2_9.root")
o2_10 = ofile.replace(".root","_2_10.root")

o2_11 = ofile.replace(".root","_2_11.root")
o2_12 = ofile.replace(".root","_2_12.root")
o2_13 = ofile.replace(".root","_2_13.root")
o2_14 = ofile.replace(".root","_2_14.root")
o2_15 = ofile.replace(".root","_2_15.root")
o2_16 = ofile.replace(".root","_2_16.root")
o2_17 = ofile.replace(".root","_2_17.root")
o2_18 = ofile.replace(".root","_2_18.root")
o2_19 = ofile.replace(".root","_2_19.root")
o2_20 = ofile.replace(".root","_2_20.root")


o3_1 =  ofile.replace(".root","_3_1.root")
o3_2  = ofile.replace(".root","_3_2.root")
f_start = 0
f_end = 25
f_end2 = 50
f_end3 = 75
f_end4 = 100
f_end5 = 125
f_end6 = 150
f_end7 = 175
f_end8 = 200
f_end9 = 225
f_end10 = 250
f_end11 = 275
f_end12 = 300
f_end13 = 325
f_end14 = 350
f_end15 = 375
f_end16 = 400
f_end17 = 425
f_end18 = 450
f_end19 = 475
f_end20 = 500


command1 = "python extraFractionFitter_recoilCut_2.py -r "+rfile+" -x -H "+tfile+" -q h_enu_q2_bin_0_mc -P h_pzmu_ptmu_bin_0_mc -o "+o1+" -b -s"

os.system(command1)
#print command1
"""
command2_1 = "python extraFractionFitter_recoilCut_2.py -r "+rfile+" -x -H "+tfile+" -q h_enu_q2_bin_0_mc -P h_pzmu_ptmu_bin_0_mc -o "+o2_1+" -b -F -s -S "+str(f_start)+" -E "+str(f_end)
#print command2
os.system(command2_1)
#sys.exit()

command2_2 = "python extraFractionFitter_recoilCut_2.py -r "+rfile+" -x -H "+tfile+" -q h_enu_q2_bin_0_mc -P h_pzmu_ptmu_bin_0_mc -o "+o2_2+" -b -F -s -S "+str(f_end)+" -E "+str(f_end2)
#print command2
os.system(command2_2)
#sys.exit()

command2_3 = "python extraFractionFitter_recoilCut_2.py -r "+rfile+" -x -H "+tfile+" -q h_enu_q2_bin_0_mc -P h_pzmu_ptmu_bin_0_mc -o "+o2_3+" -b -F -s -S "+str(f_end2)+" -E "+str(f_end3)
#print command2
os.system(command2_3)
#sys.exit()

command2_4 = "python extraFractionFitter_recoilCut_2.py -r "+rfile+" -x -H "+tfile+" -q h_enu_q2_bin_0_mc -P h_pzmu_ptmu_bin_0_mc -o "+o2_4+" -b -F -s -S "+str(f_end3)+" -E "+str(f_end4)
#print command2
os.system(command2_4)
#sys.exit()

command2_5 = "python extraFractionFitter_recoilCut_2.py -r "+rfile+" -x -H "+tfile+" -q h_enu_q2_bin_0_mc -P h_pzmu_ptmu_bin_0_mc -o "+o2_5+" -b -F -s -S "+str(f_end4)+" -E "+str(f_end5)
#print command2
os.system(command2_5)
#sys.exit()

command2_6 = "python extraFractionFitter_recoilCut_2.py -r "+rfile+" -x -H "+tfile+" -q h_enu_q2_bin_0_mc -P h_pzmu_ptmu_bin_0_mc -o "+o2_6+" -b -F -s -S "+str(f_end5)+" -E "+str(f_end6)
#print command2
os.system(command2_6)
#sys.exit()

command2_7 = "python extraFractionFitter_recoilCut_2.py -r "+rfile+" -x -H "+tfile+" -q h_enu_q2_bin_0_mc -P h_pzmu_ptmu_bin_0_mc -o "+o2_7+" -b -F -s -S "+str(f_end6)+" -E "+str(f_end7)
#print command2
os.system(command2_7)
#sys.exit()

command2_8 = "python extraFractionFitter_recoilCut_2.py -r "+rfile+" -x -H "+tfile+" -q h_enu_q2_bin_0_mc -P h_pzmu_ptmu_bin_0_mc -o "+o2_8+" -b -F -s -S "+str(f_end7)+" -E "+str(f_end8)
#print command2
os.system(command2_8)
#sys.exit()

command2_9 = "python extraFractionFitter_recoilCut_2.py -r "+rfile+" -x -H "+tfile+" -q h_enu_q2_bin_0_mc -P h_pzmu_ptmu_bin_0_mc -o "+o2_9+" -b -F -s -S "+str(f_end8)+" -E "+str(f_end9)
#print command2
os.system(command2_9)
#sys.exit()

command2_10 = "python extraFractionFitter_recoilCut_2.py -r "+rfile+" -x -H "+tfile+" -q h_enu_q2_bin_0_mc -P h_pzmu_ptmu_bin_0_mc -o "+o2_10+" -b -F -s -S "+str(f_end9)+" -E "+str(f_end10)
#print command2
os.system(command2_10)
#sys.exit()

command2_11 = "python extraFractionFitter_recoilCut_2.py -r "+rfile+" -x -H "+tfile+" -q h_enu_q2_bin_0_mc -P h_pzmu_ptmu_bin_0_mc -o "+o2_11+" -b -F -s -S "+str(f_end10)+" -E "+str(f_end11)
#print command2
os.system(command2_11)
#sys.exit()

command2_12 = "python extraFractionFitter_recoilCut_2.py -r "+rfile+" -x -H "+tfile+" -q h_enu_q2_bin_0_mc -P h_pzmu_ptmu_bin_0_mc -o "+o2_12+" -b -F -s -S "+str(f_end11)+" -E "+str(f_end12)
#print command2
os.system(command2_12)
#sys.exit()

command2_13 = "python extraFractionFitter_recoilCut_2.py -r "+rfile+" -x -H "+tfile+" -q h_enu_q2_bin_0_mc -P h_pzmu_ptmu_bin_0_mc -o "+o2_13+" -b -F -s -S "+str(f_end12)+" -E "+str(f_end13)
#print command2
os.system(command2_13)
#sys.exit()

command2_14 = "python extraFractionFitter_recoilCut_2.py -r "+rfile+" -x -H "+tfile+" -q h_enu_q2_bin_0_mc -P h_pzmu_ptmu_bin_0_mc -o "+o2_14+" -b -F -s -S "+str(f_end13)+" -E "+str(f_end14)
#print command2
os.system(command2_14)
#sys.exit()

command2_15 = "python extraFractionFitter_recoilCut_2.py -r "+rfile+" -x -H "+tfile+" -q h_enu_q2_bin_0_mc -P h_pzmu_ptmu_bin_0_mc -o "+o2_15+" -b -F -s -S "+str(f_end14)+" -E "+str(f_end15)
#print command2
os.system(command2_15)
#sys.exit()

command2_16 = "python extraFractionFitter_recoilCut_2.py -r "+rfile+" -x -H "+tfile+" -q h_enu_q2_bin_0_mc -P h_pzmu_ptmu_bin_0_mc -o "+o2_16+" -b -F -s -S "+str(f_end15)+" -E "+str(f_end16)
#print command2
os.system(command2_16)
#sys.exit()


command2_17 = "python extraFractionFitter_recoilCut_2.py -r "+rfile+" -x -H "+tfile+" -q h_enu_q2_bin_0_mc -P h_pzmu_ptmu_bin_0_mc -o "+o2_17+" -b -F -s -S "+str(f_end16)+" -E "+str(f_end17)
#print command2
os.system(command2_17)
#sys.exit()

command2_18 = "python extraFractionFitter_recoilCut_2.py -r "+rfile+" -x -H "+tfile+" -q h_enu_q2_bin_0_mc -P h_pzmu_ptmu_bin_0_mc -o "+o2_18+" -b -F -s -S "+str(f_end17)+" -E "+str(f_end18)
#print command2
os.system(command2_18)
#sys.exit()

command2_19 = "python extraFractionFitter_recoilCut_2.py -r "+rfile+" -x -H "+tfile+" -q h_enu_q2_bin_0_mc -P h_pzmu_ptmu_bin_0_mc -o "+o2_19+" -b -F -s -S "+str(f_end18)+" -E "+str(f_end19)
#print command2
os.system(command2_19)
#sys.exit()

command2_20 = "python extraFractionFitter_recoilCut_2.py -r "+rfile+" -x -H "+tfile+" -q h_enu_q2_bin_0_mc -P h_pzmu_ptmu_bin_0_mc -o "+o2_20+" -b -F -s -S "+str(f_end19)+" -E "+str(f_end20)
#print command2
os.system(command2_20)
#sys.exit()
"""
command3_1 = "python extraFractionFitter_recoilCut_2.py -r "+rfile+" -x -H "+tfile+" -q h_enu_q2_bin_0_mc -P h_pzmu_ptmu_bin_0_mc -o "+o3_1+" -b -G -s"
#print command2
#os.system(command3_1)

command3_2 = "python extraFractionFitter_recoilCut_2.py -r "+rfile+" -x -H "+tfile+" -q h_enu_q2_bin_0_mc -P h_pzmu_ptmu_bin_0_mc -o "+o3_2+" -b -Z -s"
#print command2
#os.system(command3_2)


#sys.exit()


ofile1 = ROOT.TFile(o1,"READ")
#ofile2 = ROOT.TFile(o2,"READ")
ofile2_1 = ROOT.TFile(o2_1,"READ")
ofile2_2 = ROOT.TFile(o2_2,"READ")
ofile2_3 = ROOT.TFile(o2_3,"READ")
ofile2_4 = ROOT.TFile(o2_4,"READ")
ofile2_5 = ROOT.TFile(o2_5,"READ")

ofile2_6 = ROOT.TFile(o2_6,"READ")
ofile2_7 = ROOT.TFile(o2_7,"READ")
ofile2_8 = ROOT.TFile(o2_8,"READ")
ofile2_9 = ROOT.TFile(o2_9,"READ")
ofile2_10 = ROOT.TFile(o2_10,"READ")
ofile2_11 = ROOT.TFile(o2_11,"READ")
ofile2_12 = ROOT.TFile(o2_12,"READ")
ofile2_13 = ROOT.TFile(o2_13,"READ")
ofile2_14 = ROOT.TFile(o2_14,"READ")
ofile2_15 = ROOT.TFile(o2_15,"READ")
ofile2_16 = ROOT.TFile(o2_16,"READ")
ofile2_17 = ROOT.TFile(o2_17,"READ")
ofile2_18 = ROOT.TFile(o2_18,"READ")
ofile2_19 = ROOT.TFile(o2_19,"READ")
ofile2_20 = ROOT.TFile(o2_20,"READ")


ofile3_1 = ROOT.TFile(o3_1,"READ")
ofile3_2 = ROOT.TFile(o3_2,"READ")
tfile = ROOT.TFile(tfile,"READ")
final_ofile = ROOT.TFile(ofile,"RECREATE")

hfile = tfile.Get("h_enu_q2_bin_0_mc")

keys = ofile1.GetListOfKeys()
for key in keys:
    stuff1 = ofile1.Get(key.GetName())
    #stuff2 = ofile2.Get(key.GetName())
    stuff2_1 = ofile2_1.Get(key.GetName())
    stuff2_2 = ofile2_2.Get(key.GetName())
    stuff2_3 = ofile2_3.Get(key.GetName())
    stuff2_4 = ofile2_4.Get(key.GetName())
    stuff2_5 = ofile2_5.Get(key.GetName())
    stuff2_6 = ofile2_6.Get(key.GetName())
    stuff2_7 = ofile2_7.Get(key.GetName())
    stuff2_8 = ofile2_8.Get(key.GetName())
    stuff2_9 = ofile2_9.Get(key.GetName())
    stuff2_10 = ofile2_10.Get(key.GetName())
    stuff2_11 = ofile2_11.Get(key.GetName())
    stuff2_12 = ofile2_12.Get(key.GetName())
    stuff2_13 = ofile2_13.Get(key.GetName())
    stuff2_14 = ofile2_14.Get(key.GetName())
    stuff2_15 = ofile2_15.Get(key.GetName())
    stuff2_16 = ofile2_16.Get(key.GetName())
    stuff2_17 = ofile2_17.Get(key.GetName())
    stuff2_18 = ofile2_18.Get(key.GetName())
    stuff2_19 = ofile2_19.Get(key.GetName())
    stuff2_20 = ofile2_20.Get(key.GetName())                
    stuff3_1 = ofile3_1.Get(key.GetName())
    stuff3_2  = ofile3_2.Get(key.GetName())
    #print stuff1.ClassName()
    if stuff1.ClassName()=="PlotUtils::MnvH1D" or stuff1.ClassName()=="PlotUtils::MnvH2D":
    	#print "Hey you!!"
    	#stuff1.PopVertErrorBand("Flux")
	#stuff1.AddVertErrorBandAndFillWithCV("Flux",200)
	#stuff2.TransferVertErrorBand(stuff1,"Flux")
	vnames = stuff3_1.GetVertErrorBandNames()
	for name in vnames:
	    if name in genie_list1:
	        stuff1.PopVertErrorBand(name)
	        stuff3_1.TransferVertErrorBand(stuff1,name)
	    if name in genie_list2:
	        stuff1.PopVertErrorBand(name)
		stuff3_2.TransferVertErrorBand(stuff1,name)
            if "Flux" in name:
                stuff1.PopVertErrorBand("Flux")
                stuff1.AddVertErrorBand("Flux",f_end20)
                for i in range(f_start,f_end):
                    stuff1.GetVertErrorBand("Flux").GetHist(i).Add(stuff2_1.GetVertErrorBand("Flux").GetHist(i))
                    stuff1.GetVertErrorBand("Flux").SetUnivWgt(i,stuff2_1.GetVertErrorBand("Flux").GetUnivWgt(i))
                for i in range(f_end,f_end2):
                    stuff1.GetVertErrorBand("Flux").GetHist(i).Add(stuff2_2.GetVertErrorBand("Flux").GetHist(i))
                    stuff1.GetVertErrorBand("Flux").SetUnivWgt(i,stuff2_2.GetVertErrorBand("Flux").GetUnivWgt(i))
                for i in range(f_end2,f_end3):
                    stuff1.GetVertErrorBand("Flux").GetHist(i).Add(stuff2_3.GetVertErrorBand("Flux").GetHist(i))
                    stuff1.GetVertErrorBand("Flux").SetUnivWgt(i,stuff2_3.GetVertErrorBand("Flux").GetUnivWgt(i))
                for i in range(f_end3,f_end4):
                    stuff1.GetVertErrorBand("Flux").GetHist(i).Add(stuff2_4.GetVertErrorBand("Flux").GetHist(i))
                    stuff1.GetVertErrorBand("Flux").SetUnivWgt(i,stuff2_4.GetVertErrorBand("Flux").GetUnivWgt(i))
                for i in range(f_end4,f_end5):
                    stuff1.GetVertErrorBand("Flux").GetHist(i).Add(stuff2_5.GetVertErrorBand("Flux").GetHist(i))
                    stuff1.GetVertErrorBand("Flux").SetUnivWgt(i,stuff2_5.GetVertErrorBand("Flux").GetUnivWgt(i))
		   
                for i in range(f_end5,f_end6):
                    stuff1.GetVertErrorBand("Flux").GetHist(i).Add(stuff2_6.GetVertErrorBand("Flux").GetHist(i))
                    stuff1.GetVertErrorBand("Flux").SetUnivWgt(i,stuff2_6.GetVertErrorBand("Flux").GetUnivWgt(i))
                for i in range(f_end6,f_end7):
                    stuff1.GetVertErrorBand("Flux").GetHist(i).Add(stuff2_7.GetVertErrorBand("Flux").GetHist(i))
                    stuff1.GetVertErrorBand("Flux").SetUnivWgt(i,stuff2_7.GetVertErrorBand("Flux").GetUnivWgt(i))
                for i in range(f_end7,f_end8):
                    stuff1.GetVertErrorBand("Flux").GetHist(i).Add(stuff2_8.GetVertErrorBand("Flux").GetHist(i))
                    stuff1.GetVertErrorBand("Flux").SetUnivWgt(i,stuff2_8.GetVertErrorBand("Flux").GetUnivWgt(i))
                for i in range(f_end8,f_end9):
                    stuff1.GetVertErrorBand("Flux").GetHist(i).Add(stuff2_9.GetVertErrorBand("Flux").GetHist(i))
                    stuff1.GetVertErrorBand("Flux").SetUnivWgt(i,stuff2_9.GetVertErrorBand("Flux").GetUnivWgt(i))
                for i in range(f_end9,f_end10):
                    stuff1.GetVertErrorBand("Flux").GetHist(i).Add(stuff2_10.GetVertErrorBand("Flux").GetHist(i))
                    stuff1.GetVertErrorBand("Flux").SetUnivWgt(i,stuff2_10.GetVertErrorBand("Flux").GetUnivWgt(i))
                for i in range(f_end10,f_end11):
                    stuff1.GetVertErrorBand("Flux").GetHist(i).Add(stuff2_11.GetVertErrorBand("Flux").GetHist(i))
                    stuff1.GetVertErrorBand("Flux").SetUnivWgt(i,stuff2_11.GetVertErrorBand("Flux").GetUnivWgt(i))
                for i in range(f_end11,f_end12):
                    stuff1.GetVertErrorBand("Flux").GetHist(i).Add(stuff2_12.GetVertErrorBand("Flux").GetHist(i))
                    stuff1.GetVertErrorBand("Flux").SetUnivWgt(i,stuff2_12.GetVertErrorBand("Flux").GetUnivWgt(i))
                for i in range(f_end12,f_end13):
                    stuff1.GetVertErrorBand("Flux").GetHist(i).Add(stuff2_13.GetVertErrorBand("Flux").GetHist(i))
                    stuff1.GetVertErrorBand("Flux").SetUnivWgt(i,stuff2_13.GetVertErrorBand("Flux").GetUnivWgt(i))
                for i in range(f_end13,f_end14):
                    stuff1.GetVertErrorBand("Flux").GetHist(i).Add(stuff2_14.GetVertErrorBand("Flux").GetHist(i))
                    stuff1.GetVertErrorBand("Flux").SetUnivWgt(i,stuff2_14.GetVertErrorBand("Flux").GetUnivWgt(i))
                for i in range(f_end14,f_end15):
                    stuff1.GetVertErrorBand("Flux").GetHist(i).Add(stuff2_15.GetVertErrorBand("Flux").GetHist(i))
                    stuff1.GetVertErrorBand("Flux").SetUnivWgt(i,stuff2_15.GetVertErrorBand("Flux").GetUnivWgt(i))
                for i in range(f_end15,f_end16):
                    stuff1.GetVertErrorBand("Flux").GetHist(i).Add(stuff2_16.GetVertErrorBand("Flux").GetHist(i))
                    stuff1.GetVertErrorBand("Flux").SetUnivWgt(i,stuff2_16.GetVertErrorBand("Flux").GetUnivWgt(i))
                for i in range(f_end16,f_end17):
                    stuff1.GetVertErrorBand("Flux").GetHist(i).Add(stuff2_17.GetVertErrorBand("Flux").GetHist(i))
                    stuff1.GetVertErrorBand("Flux").SetUnivWgt(i,stuff2_17.GetVertErrorBand("Flux").GetUnivWgt(i))
                for i in range(f_end17,f_end18):
                    stuff1.GetVertErrorBand("Flux").GetHist(i).Add(stuff2_18.GetVertErrorBand("Flux").GetHist(i))
                    stuff1.GetVertErrorBand("Flux").SetUnivWgt(i,stuff2_18.GetVertErrorBand("Flux").GetUnivWgt(i))
                for i in range(f_end18,f_end19):
                    stuff1.GetVertErrorBand("Flux").GetHist(i).Add(stuff2_19.GetVertErrorBand("Flux").GetHist(i))
                    stuff1.GetVertErrorBand("Flux").SetUnivWgt(i,stuff2_19.GetVertErrorBand("Flux").GetUnivWgt(i))
                for i in range(f_end19,f_end20):
                    stuff1.GetVertErrorBand("Flux").GetHist(i).Add(stuff2_20.GetVertErrorBand("Flux").GetHist(i))
                    stuff1.GetVertErrorBand("Flux").SetUnivWgt(i,stuff2_20.GetVertErrorBand("Flux").GetUnivWgt(i))		    		    		    
                
                
                
        stuff1.Write()
    elif stuff1.ClassName()=="TVector2":
        stuff1.Write("pot")
    else:
        stuff1.Write()

final_ofile.Close()
       

