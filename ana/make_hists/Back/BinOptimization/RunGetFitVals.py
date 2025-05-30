import os,sys

variables = ["Emu_Ehad","W_Q2","x_Q2","x_y"]
#variables = ["Emu","Ehad", "x", "y"]

for variable in variables:
    cmd = "python GetFitVals.py %s %s"%(sys.argv[1],variable)
    os.system(cmd)
