import os,sys

variables = ["pn","dalphat","dphit","dpt","dptx","dpty","dpl","ptheta","tp"]

for v in variables:
    cmd = "./plot_nu_xsec2_transverse %s %s %s Standard"%(sys.argv[1],v,"ptmu")
    os.system(cmd)
