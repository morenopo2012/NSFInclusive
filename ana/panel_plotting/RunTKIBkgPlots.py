import os,sys

#./plot_nu_transverse_bkg_constrained /minerva/data/users/drut1186/ME_Transverse/Processing_07032020_SameUpdates_OlderTuples/Transverse_CV_Elastic_FSI/ 0 pn ptmu
indir = sys.argv[1]
variable_dict =[]
#variable_dict.append("muonPz")#0
variable_dict.append("dalphat")
variable_dict.append("dphit")
variable_dict.append("pn")
variable_dict.append("dpt")
variable_dict.append("dptx")#5
variable_dict.append("dpty")
variable_dict.append("tp")
variable_dict.append("ptheta")
variable_dict.append("signed")
#variable_dict.append("signeddalphat")
#variable_dict.append("signeddphit")#11
#variable_dict.append("dthetaR")
#variable_dict.append("dthetaP")
#variable_dict.append("dPR")
#variable_dict.append("dPP")
#variable_dict.append("dPRi")
#variable_dict.append("dPPi")#17
variable_dict.append("dpl")


for i in range(2):
    for v in variable_dict:
        cmd = "./plot_nu_transverse_bkg_constrained %s %d %s ptmu"%(indir,i,v)
        os.system(cmd)
