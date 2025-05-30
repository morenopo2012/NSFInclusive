import sys,os

loc = sys.argv[1]
target = "EventSelection/background2D"
ratio = ["0","1"]
scale_qelike = ["0","1"]

for r in ratio:
    for s2 in scale_qelike:
        cmdclean = "rm *.png *.C *.eps"
        cmdmake = "./plot_nu_bkg_samples_constrained %s %s %s"%(loc,r,s2)
        cmdmv = "mv *.png *.C *.eps %s"%(target)
        os.system("mkdir -p %s"%(target))
        os.system(cmdclean)
        os.system(cmdmake)
        os.system(cmdmv)

loc = sys.argv[1]
target = "EventSelection/background2D_Samples"
ratio = ["0","1"]

for r in ratio:
    cmdclean = "rm *.png *.C *.eps"
    cmdmake = "./plot_nu_bkg_samples %s %s"%(loc,r)
    cmdmv = "mv *.png *.C *.eps %s"%(target)
    os.system("mkdir -p %s"%(target))
    os.system(cmdclean)
    os.system(cmdmake)
    os.system(cmdmv)

#sys.exit()
loc = sys.argv[1]
target = "EventSelection/background3D"
ratio = ["0","1"]
scale_qelike = ["0","1"]

for r in ratio:
    for s2 in scale_qelike:
        cmdclean = "rm *.png *.C *.eps"
        cmdmake = "./plot_nu_bkg %s %s %s"%(loc,r,s2)
        cmdmv = "mv *.png *.C *.eps %s"%(target)
        os.system("mkdir -p %s"%(target))
        os.system(cmdclean)
        os.system(cmdmake)
        os.system(cmdmv)
#        sys.exit()
