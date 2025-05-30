import os,sys


indir = "/minerva/data/users/kleykamp/makehistoutput/xsec/XSec_Final_2021-04-25_test/muon_p_2d/flux_corrected/"


cmd = "./plot_nuclearTargetRatios_2D "
cmd = cmd + "%stracker/out.root "%(indir)
cmd = cmd + "%slead/out.root "%(indir)
cmd = cmd + "%siron/out.root "%(indir)
cmd = cmd + "%scarbon/out.root "%(indir)
cmd = cmd + "%swater/out.root "%(indir)

os.system(cmd)
