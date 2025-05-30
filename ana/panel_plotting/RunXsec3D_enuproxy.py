import sys,os

loc = sys.argv[1]
outdir = sys.argv[2]
target = ["CrossSection_Enuproxy/standard","CrossSection_Enuproxy/resfsi","CrossSection_Enuproxy/qefsi","CrossSection_Enuproxy/resisi","CrossSection_Enuproxy/2p2htunes","CrossSection_Enuproxy/2p2h_nn_np"]
target2 = "CrossSection_Enuproxy/models"
firstset = ["1 0 0 0 0 0","0 1 0 0 0","0 0 1 0 0 0","0 0 0 1 0 0","0 0 0 0 1 0","0 0 0 0 0 1"]
secondset = ["0","1"]
thirdset = ["0", "1"]

for i,s1 in enumerate(firstset):
    for s2 in secondset:
        for s3 in thirdset:
            for j in range(1,11):
                if(j!=10) : continue
                cmdclean = "rm *.png *.C *.eps"
                cmdmake = "./plot_nu_xsec3 %s %s %s %s %s 0"%(loc,s1,s2,s3,j)
                cmdmv = "mv *.png *.C *.eps %s/%s/iteration_%d"%(outdir,target[i],j)
                os.system("mkdir -p %s/%s/iteration_%d"%(outdir,target[i],j))
                os.system(cmdclean)
                os.system(cmdmake)
                os.system(cmdmv)

sys.exit("I'm not doing the model comparisons")
for s2 in secondset:
    for s3 in thirdset:
        for j in range(1,11):
            if(j!=10) : continue
            cmdclean = "rm *.png *.C *.eps"
            cmdmake = "./plot_nu_xsec3_models %s %s %s %d"%(loc,s2,s3,j)
            cmdmv = "mv *.png *.C *.eps %s/%s/iteraction_%d"%(outdir,target2,j)
            os.system("mkdir -p %s/%s/iteration_%d"%(outdir,target2,j))
            os.system(cmdclean)
            os.system(cmdmake)
            os.system(cmdmv)



