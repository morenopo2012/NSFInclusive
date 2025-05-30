import os,sys

variable_dict =[]
variable_dict.append("muonPz")#0
variable_dict.append("dalphat")
variable_dict.append("dphit")
variable_dict.append("pn")
variable_dict.append("dpt")
variable_dict.append("dptx")#5
variable_dict.append("dpty")
variable_dict.append("tp")
variable_dict.append("ptheta")
variable_dict.append("signed")
variable_dict.append("signeddalphat")
variable_dict.append("signeddphit")#11
variable_dict.append("dthetaR")
variable_dict.append("dthetaP")
variable_dict.append("dPR")
variable_dict.append("dPP")
variable_dict.append("dPRi")
variable_dict.append("dPPi")#17


for xvar in range(0,len(variable_dict)):
    for yvar in range(0,len(variable_dict)):
        if yvar>=xvar: continue;
#        cmd = ./plot_nu_transverse_efficiency /minerva/data/users/drut1186/ME_Transverse/Signal_TKITKI_06172020_rerun/Transverse dpty dalphat
        cmd = "./plot_nu_transverse_efficiency %s %s %s"%(sys.argv[1],variable_dict[xvar],variable_dict[yvar])
        cmd_mkdir = "mkdir -p TKITKI/%s_%s"%(variable_dict[xvar],variable_dict[yvar])
        cmd_save = "mv *.png *.eps *.C TKITKI/%s_%s"%(variable_dict[xvar],variable_dict[yvar])
        os.system(cmd_mkdir)
        os.system(cmd)
        os.system(cmd_save)
