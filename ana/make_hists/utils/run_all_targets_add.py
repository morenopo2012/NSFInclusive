import os
import sys


#*************irrelevant variables that don't get used********************
mats = []
mats.append("z06")
mats.append("z26")
mats.append("z82")

targets = []
targets.append("t1")
targets.append("t2")
targets.append("t3")
targets.append("t4")
targets.append("t5")

target_ints = []
target_ints.append("1")
target_ints.append("2")
target_ints.append("3")
target_ints.append("4")
target_ints.append("5")
target_ints.append("14")
target_ints.append("24")
target_ints.append("34")
target_ints.append("44")
target_ints.append("54")
target_ints.append("64")
target_ints.append("74")
target_ints.append("84")
target_ints.append("94")

#*********
carbon = [3]
iron   = [1,2,3,5]
lead   = [1,2,3,4,5]

strs = []
strs_int = []

for i in carbon:
    tgt = "t"+str(i)+"_z06"
    strs.append(tgt)
    strs_int.append("3")

for i in iron:
    tgt= "t"+str(i)+"_z26"
    strs.append(tgt)
    strs_int.append("1")

for i in lead:
    tgt= "t"+str(i)+"_z82"
    strs.append(tgt)
    strs_int.append("1")
#Trackers
#for j in range(1,10):
#    tgt = "t"+str(j)+"4_z82"
#tmp since trackers needed to be rerun
tgt = "t14_z82"
strs.append(tgt)
strs_int.append("14")

if "SideBand" in str(sys.argv[1]):
    del strs[:]
    del strs_int[:]
    strs.append("t3_z06")
    strs.append("t1_z26")
    strs.append("t1_z82")
    strs.append("t14_z82")
    strs_int.append("3")
    strs_int.append("1")
    strs_int.append("1")
    strs_int.append("14")

if "Plastic" in str(sys.argv[1]):
    del strs[:]
    del strs_int[:]
    strs.append("t3_z06")
    strs.append("t1_z26")
    strs.append("t1_z82")
    strs_int.append("3")
    strs_int.append("1")
    strs_int.append("1")

    

print strs

harg=str(sys.argv[1])
nchar=harg.find("NukeOnlyMC")
for targ, targetID in zip(strs,strs_int):
    hname=harg
    if( nchar!=-1  and int(targetID)>10 ): 
        hname= harg[0:nchar] + harg[nchar+8:]
    plot = hname+"_"+str(targ)+"_Nu_"+os.getenv("NUKECC_TAG")
    print plot
    command = "./PlaylistAdder "+plot+" "+targetID+" minervame1A minervame1B minervame1C minervame1D minervame1E minervame1F minervame1G minervame1L minervame1M minervame1N minervame1O minervame1P"
    os.system(command)
    print command
