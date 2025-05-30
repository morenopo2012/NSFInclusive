import sys,os

loc = sys.argv[1]
low_recoil_removed = ["0","1"]#unique
resfsi = ["0","1"]#unique
qefsi = ["0","1"]#unique
recoil_fit = ["0","1"]#unique
ratio = ["0","1"]
zoom = ["0","1"]
constraint = ["0","1"]

for r in ratio:

    for l in low_recoil_removed:
        target = "EventSelection/standard"
        os.system("mkdir -p %s"%(target))
        cmdclean = "rm *.png *.C *.eps"
        cmdmake = "./plot_nu_eventrate_3D %s %s %s %s %s %s %s 0 0"%(loc,l,resfsi[0],qefsi[0],recoil_fit[0],r,zoom[0])
        cmdmv = "mv *.png *.C *.eps %s"%(target)
        os.system(cmdclean)
        os.system(cmdmake)
        os.system(cmdmv)

    for res in resfsi:
        target = "EventSelection/resfsi"
        os.system("mkdir -p %s"%(target))
        cmdclean = "rm *.png *.C *.eps"
        cmdmake = "./plot_nu_eventrate_3D %s %s %s %s %s %s %s 0 0"%(loc,low_recoil_removed[0],res,qefsi[0],recoil_fit[0],r,zoom[0])
        cmdmv = "mv *.png *.C *.eps %s"%(target)
        os.system(cmdclean)
        os.system(cmdmake)
        os.system(cmdmv)

    for qe in qefsi:
        target = "EventSelection/qefsi"
        os.system("mkdir -p %s"%(target))
        cmdclean = "rm *.png *.C *.eps"
        cmdmake = "./plot_nu_eventrate_3D %s %s %s %s %s %s %s 0 0"%(loc,low_recoil_removed[0],resfsi[0],qe,recoil_fit[0],r,zoom[0])
        cmdmv = "mv *.png *.C *.eps %s"%(target)
        os.system(cmdclean)
        os.system(cmdmake)
        os.system(cmdmv)

    for rec in recoil_fit:
        target = "EventSelection/recoil_fit"
        os.system("mkdir -p %s"%(target))
        cmdclean = "rm *.png *.C *.eps"
        cmdmake = "./plot_nu_eventrate_3D %s %s %s %s %s %s %s 0 0"%(loc,low_recoil_removed[0],resfsi[0],qefsi[0],rec,r,zoom[0])
        cmdmv = "mv *.png *.C *.eps %s"%(target)
        os.system(cmdclean)
        os.system(cmdmake)
        os.system(cmdmv)
        
    target = "EventSelection/no_constraint"
    os.system("mkdir -p %s"%(target))
    cmdclean = "rm *.png *.C *.eps"
    cmdmake = "./plot_nu_eventrate_3D %s %s %s %s %s %s %s 1 0"%(loc,l,resfsi[0],qefsi[0],recoil_fit[0],r,zoom[0])
    cmdmv = "mv *.png *.C *.eps %s"%(target)
    os.system(cmdclean)
    os.system(cmdmake)
    os.system(cmdmv)

    target = "EventSelection/constraint"
    os.system("mkdir -p %s"%(target))
    cmdclean = "rm *.png *.C *.eps"
    cmdmake = "./plot_nu_eventrate_3D %s %s %s %s %s %s %s 1 1"%(loc,l,resfsi[0],qefsi[0],recoil_fit[0],r,zoom[0])
    cmdmv = "mv *.png *.C *.eps %s"%(target)
    os.system(cmdclean)
    os.system(cmdmake)
    os.system(cmdmv)
