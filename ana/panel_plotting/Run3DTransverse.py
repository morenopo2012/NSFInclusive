import os,sys,multiprocessing

def runCmd(mycmd):
    os.system(mycmd)

variable_dict =[]
variable_dict.append("muonPz")#0
variable_dict.append("dalphat")
variable_dict.append("dphit")
variable_dict.append("pn")
variable_dict.append("dpt")
variable_dict.append("dptx")#5
variable_dict.append("dpty")
variable_dict.append("tp")
#variable_dict.append("ptheta")#8
#variable_dict.append("signed")#9
#variable_dict.append("signeddalphat")#10
#variable_dict.append("signeddphit")#11
#variable_dict.append("dthetaR")
#variable_dict.append("dthetaP")
#variable_dict.append("dPR")
#variable_dict.append("dPP")
#variable_dict.append("dPRi")
#variable_dict.append("dPPi")#17
variable_dict.append("recoil")#18
variable_dict.append("dpl")#19
variable_dict.append("thetapn")#20


commands = []

for i in range(len(variable_dict)):
    for j in range(len(variable_dict)):
        if(j>=i): continue
        if(variable_dict[i]=="recoil" or variable_dict[j]=="recoil"):continue
        cmd = " rm *.png *.eps *.C;mkdir -p TKI3D/%s_%s/ratio; rm TKI3D/%s_%s/ratio/*;./plot_nu_eventrate_3D_transverse /pnfs/minerva/persistent/users/drut1186/Test_TKI_TKI/Test3_Renamed/Transverse_CV/selection_Signal/MuonEventSelection_MakeFlux-1_Multiplicity-2_Sample-Signal_%s_%s_CombinedPlaylists.root 1 %s %s;mv *.png TKI3D/%s_%s/ratio"%(variable_dict[i],variable_dict[j],variable_dict[i],variable_dict[j],variable_dict[i],variable_dict[j],variable_dict[i],variable_dict[j],variable_dict[i],variable_dict[j])
        cmd2 = " rm *.png *.eps *.C;mkdir -p TKI3D/%s_%s/rate; rm TKI3D/%s_%s/rate/*;./plot_nu_eventrate_3D_transverse /pnfs/minerva/persistent/users/drut1186/Test_TKI_TKI/Test3_Renamed/Transverse_CV/selection_Signal/MuonEventSelection_MakeFlux-1_Multiplicity-2_Sample-Signal_%s_%s_CombinedPlaylists.root 0 %s %s;mv *.png TKI3D/%s_%s/rate"%(variable_dict[i],variable_dict[j],variable_dict[i],variable_dict[j],variable_dict[i],variable_dict[j],variable_dict[i],variable_dict[j],variable_dict[i],variable_dict[j])
        #os.system(cmd)
        #os.system(cmd2)
        commands.append(cmd)
        commands.append(cmd2)


pool = multiprocessing.Pool(processes=3)
output = pool.map(runCmd,commands)

