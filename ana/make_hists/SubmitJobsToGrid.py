#sample command to run: python SubmitJobsToGrid.py --stage eventLoop --playlist minervame1A --basedir /minerva/app/users/mmehmood/MAT_AL9/ --pdg 14 --filename test --outdir /pnfs/minerva/persistent/users/mmehmood/default_analysis_loc

import os,sys,time
from optparse import OptionParser
import datetime

###########################################################
#  Script: JobSubmission python script for submitting
#          analysis jobs to the grid
###########################################################

# Write the command you used to run your analysis
def writeEventLoop(mywrapper,outdir,targetID, targetZ, playlist):
    mywrapper.write("./runEventLoop_MEC_FateWeightCRESWeighted . " " "+str(targetID)+" "+str(targetZ)+" "+playlist)
    #mywrapper.write("./VertexZRecoVars . " " "+str(targetID)+" "+str(targetZ)+" "+playlist)
    mywrapper.write("\n")
    #mywrapper.write("./ES . 6A\n") #run code (I modified)
    mywrapper.write("ifdh cp ./*.root "+outdir)

def writeMigration(mywrapper,outdir,targetID, targetZ, playlist):
    mywrapper.write("./MigrationHists . " " "+str(targetID)+" "+str(targetZ)+" " +playlist)
    #mywrapper.write("./MigrationHistsDaisy . " " "+str(targetID)+" "+str(targetZ)+" " +playlist)
    mywrapper.write("\n")
    #mywrapper.write("./M . 6A\n")
    mywrapper.write("ifdh cp ./*.root "+outdir)

def writeEfficiency(mywrapper,outdir, targetID, targetZ, playlist):
    mywrapper.write("./EfficiencyHists . " " "+str(targetID)+" "+str(targetZ)+" " +playlist)
    #mywrapper.write("./EfficiencyHistsDaisy . " " "+str(targetID)+" "+str(targetZ)+" " +playlist)
    mywrapper.write("\n")
    #mywrapper.write("./Eff . 6A\n")
    mywrapper.write("ifdh cp ./*.root "+outdir)

def writePlasticSideband(mywrapper,outdir, targetID, targetZ, playlist):
    mywrapper.write("./TargetCombined_RunEventLoop_plasticSB_2D . " " "+str(targetID)+" "+str(targetZ)+" " +playlist)
    mywrapper.write("\n")
    #mywrapper.write("./Eff . 6A\n")
    mywrapper.write("ifdh cp ./*.root "+outdir)

def writePhysicsSideband(mywrapper,outdir, targetID, targetZ, playlist):
    mywrapper.write("./TargetCombined_RunEventLoop_physSB_2D . " " "+str(targetID)+" "+str(targetZ)+" " +playlist)
    mywrapper.write("\n")
    #mywrapper.write("./Eff . 6A\n")
    mywrapper.write("ifdh cp ./*.root "+outdir)

def writeAll(mywrapper,outdir):
    mywrapper.write("./ES . 6J\n")
    mywrapper.write("ifdh cp ./*.root "+outdir)
    mywrapper.write("./M . 6J\n")
    mywrapper.write("ifdh cp ./*.root "+outdir)
    mywrapper.write("./Eff . 6J\n")
    mywrapper.write("ifdh cp ./*.root "+outdir)
 
def closureTest(mywrapper,outdir):
    mywrapper.write("./c\n")
    mywrapper.write("ifdh cp ./*.root "+outdir)

def writeSetups(mywrapper,basedir):
    # create topdir + source minerva products
    topdir = "$CONDOR_DIR_INPUT" + basedir
    mywrapper.write("cd "+topdir+"\n")
    #mywrapper.write("source /cvmfs/minerva.opensciencegrid.org/minerva/setup/setup_minerva_products.sh\n")
    # set environment variables
    #mywrapper.write("cd "+topdir+"\n")

    mywrapper.write("export IFDH_BASE_URI=https://samminerva.fnal.gov:8483/sam/minerva/api\n")
    mywrapper.write("export SAM_WEB_HOST=samminerva.fnal.gov\n") 
    mywrapper.write("source /cvmfs/larsoft.opensciencegrid.org/spack-packages/setup-env.sh\n")
    mywrapper.write("spack load cmake@3.27.7\n")
    mywrapper.write("spack load fife-utils@3.7.4\n")

    mywrapper.write("source /cvmfs/larsoft.opensciencegrid.org/spack-v0.22.0-fermi/setup-env.sh\n")
    mywrapper.write("spack load root@6.28.12 arch=linux-almalinux9-x86_64_v3\n")
    mywrapper.write("spack load gcc@12.2.0 arch=linux-almalinux9-x86_64_v3\n")

    mywrapper.write("export LD_LIBRARY_PATH=${ROOTSYS}/lib/root:${LD_LIBRARY_PATH}\n")
    mywrapper.write("cd "+topdir+"\n")
    mywrapper.write("cd opt/bin/\n")
    mywrapper.write("source setup.sh\n")
   
    mywrapper.write("export NUKECCSRCROOT=${CONDOR_DIR_INPUT}/exp/minerva/app/users/omorenop/cmtuser/git-Mat/NSFNukeCCInclusive/NUKECCSRC\n")
    mywrapper.write("export MY_NSFNUKECC=${CONDOR_DIR_INPUT}/exp/minerva/app/users/omorenop/cmtuser/git-Mat/NSFNukeCCInclusive/ana/\n")
    #mywrapper.write("export /cvmfs/larsoft.opensciencegrid.org/spack-packages/opt/spack/linux-almalinux9-x86_64_v2/gcc-12.2.0/root-6.28.12-sfwfmqorvxttrxgfrfhoq5kplou2pddd/bin/thisroot.sh\n")
    mywrapper.write("export NUKECCSRC_ANA=${NUKECCSRCROOT}/ana_common\n")
    mywrapper.write("export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${CONDOR_DIR_INPUT}/exp/minerva/app/users/omorenop/cmtuser/git-Mat/NSFNukeCCInclusive/NUKECCSRC/ana_common/src\n")

    # source bash script
    mywrapper.write("cd "+topdir+"\n")
    mywrapper.write("source setup_GRID.sh\n")

    mywrapper.write("export PYTHONPATH=/exp/minerva/app/users/omorenop/cmtuser/git-Mat/MAT-MINERvA/python:/exp/minerva/app/users/omorenop/cmtuser/git-Mat/MAT-MINERvA/python/PlotUtils\n")
    #mywrapper.write("source /cvmfs/larsoft.opensciencegrid.org/spack-packages/opt/spack/linux-almalinux9-x86_64_v2/gcc-12.2.0/root-6.28.12-sfwfmqorvxttrxgfrfhoq5kplou2pddd/bin/thisroot.sh\n")
    # get to directory to run code
    mywrapper.write("cd "+topdir+"\n")
    mywrapper.write("cd NSFNukeCCInclusive/ana/gridSubmission/\n")
    mywrapper.write("echo Rint.Logon: ./rootlogon_grid.C > ./.rootrc\n")

# Define function to create tarball
def createTarball(outdir,tag,basedir):
    found = os.path.isfile("%s/myareatar_%s.tar.gz"%(outdir,tag))
    if(not found):
        cmd = "tar -czf /exp/minerva/app/users/$USER/myareatar_%s.tar.gz %s"%(tag,basedir)
        print("Making tar",cmd)
        os.system(cmd)
        cmd2 = "cp /exp/minerva/app/users/$USER/myareatar_%s.tar.gz %s/"%(tag,outdir)
        print("Copying tar",cmd2)
        os.system(cmd2)
# Define function to unpack tarball
def writeTarballProceedure(mywrapper,tag,basedir):
    print("I will be making the tarball upacking with this version") 
    print("Path is",basedir)
    #mywrapper.write("source /cvmfs/minerva.opensciencegrid.org/minerva/software_releases/v22r1p1/setup.sh\n")
    mywrapper.write("export XRD_NETWORKSTACK=IPv4\n")
    #mywrapper.write("export UPS_OVERRIDE='-H Linux64bit+3.10-2.17'\n")
    mywrapper.write("cd $CONDOR_DIR_INPUT\n")
    mywrapper.write("tar -xvzf myareatar_%s.tar.gz\n"%tag)
    mywrapper.write("cd $CONDOR_DIR_INPUT\n")
    mywrapper.write("export MINERVA_PREFIX=${CONDOR_DIR_INPUT}/exp/minerva/app/users/omorenop/cmtuser/git-Mat/opt\n")
    writeSetups(mywrapper,basedir)

def writeOptions(parser):
    print("Now write options to the parser") 
    # Directory to write output
    parser.add_option('--outdir', dest='outdir', help='Directory to write output to', default = "/pnfs/minerva/persistent/users/"+_user_+"/default_analysis_loc/AnaMay2025")
    parser.add_option('--basedir', dest='basedir', help='Base directory for making tarball', default = "NONE")
    parser.add_option('--stage', dest='stage', help='Process type', default="NONE")
    parser.add_option('--sample', dest='sample', help='Sample type', default="NONE")
    parser.add_option('--playlist', dest='playlist', help='Playlist type', default="NONE")
    parser.add_option('--targetID',dest='targetID',help='Which target to analyze',action="store")
    parser.add_option('--targetZ',dest='targetZ',help='Which material to analyze within chosen target',action="store")
    ##########################################################################
    #  Options for making tar files....Basically you can make tarfiles 
    #######################################################################
    parser.add_option('--tag', dest='tag', help="Tag your release",default="tag_")
    parser.add_option('--mail',dest='mail',help="Want mail after job completion or failure",default=False,action="store_true")
    parser.add_option('--sametar',dest='sametar',help="Recycle the same tar file for jobsubmission",default=False,action="store_true")
    parser.add_option('--tarfilename',dest='tarfilename',help='Name of the tarfile you want to use',default="NA")
    parser.add_option('--notimestamp',dest='notimestamp',help='Flag to TURN OFF time stamp in the tag',default=False,action="store_true")
    parser.add_option('--filename',dest='filename',help='Name we want to use for output files',default="NONE")
    parser.add_option('--pdg',dest='pdg',help='The mode that we want to run the analysis in',default="NONE")

# Valid stages for the neutrinos (you can add more for your own analysis)
valid_stages=["eventLoop", "migration", "efficiency", "transwarp", "all", "closure", "plasticSideband", "physicsSideband", ]
avail_playlists=["minervame6A", "minervame6B", "minervame6C", "minervame6D", "minervame6E", "minervame6F", "minervame6G", "minervame6H", "minervame6I", "minervame6J", "minervame5A", "minervame1A", "minervame1B", "minervame1C", "minervame1D", "minervame1E", "minervame1F", "minervame1G", "minervame1L", "minervame1M", "minervame1N", "minervame1O", "minervame1P"]
# Get user's name
_user_ = os.getenv("USER")
usage = "usage: %prog[opts]" 
parser = OptionParser(usage=usage)
writeOptions(parser)
(opts,args) = parser.parse_args()

# Now print some information for user's selected options
if opts.mail:
    print("############################################") 
    print("You have opted to Get mails once the jobs are finished of failed.") 
    print("#############################################") 

##############################################
#  TODO: Here I want to put the tag scheme...
##############################################
tag_name = str(opts.tag)
time_stamp = int(time.time())
if tag_name=="tag_":
    print("YOU DIDNT SPECIFY ANY ANY TAG...SO I WILL USE MY OWN TAGGING SCHEME") 
    # You can add in your own tagging scheme
    tag_name += "default_" 

print("Is everything fine till here*********") 
# Add the time stamp to tag
if not opts.notimestamp:
    tag_name += str(time_stamp)

print("Tag for this version is ",tag_name)
print("******************************************************") 

# This is the output directory after the job is finished
output_dir = "$CONDOR_DIR_HISTS/" 

# Make outdir if not exist
if(not os.path.isdir(opts.outdir)):
    os.makedirs(opts.outdir)

memory = 60000
#memory = 10000
#memory = 3500
disk = 60 #added in May20 b/c jobs going to 'held'

# Check the stage is valid or not
if opts.stage not in valid_stages:
    print(opts.stage,"Selected stage is not valid. Here are the valid options",valid_stages)
    sys.exit()

##############################################
#  Create wrapper
##############################################

wrapper_name = "%s_%s_wrapper_%s.sh"%(opts.stage,opts.playlist,tag_name)

mywrapper = open(wrapper_name,"w")
mywrapper.write("#!/bin/sh\n")
# Now create tarball
if (opts.sametar==False):
    createTarball(opts.outdir,tag_name,opts.basedir)
else:
    tarname = str(opts.tarfilename)
    if not os.path.exists(opts.outdir+opts.tarfilename):
        print("Tar File "+opts.outdir+opts.tarfilename+" doesn't Exist!") 
        sys.exit()
    #change the tag to the current one...
    cmd="cp "+opts.outdir+opts.tarfilename+" "+opts.outdir+"myareatar_"+tag_name+".tar.gz" 
    os.system(cmd)

# This will unpack the tarball we just made above
writeTarballProceedure(mywrapper,tag_name,opts.basedir)

# Now the add the command to run event loop 
if(opts.stage=="eventLoop"):
    writeEventLoop(mywrapper,opts.outdir,opts.targetID, opts.targetZ, opts.playlist)
elif (opts.stage=="migration"):
    writeMigration(mywrapper,opts.outdir,opts.targetID, opts.targetZ, opts.playlist)
elif (opts.stage=="efficiency"): 
    writeEfficiency(mywrapper,opts.outdir,opts.targetID, opts.targetZ, opts.playlist)
elif (opts.stage=="transwarp"):
    writeTransWarp(mywrapper,opts.outdir)
elif (opts.stage=="all"):
    writeAll(mywrapper,opts.outdir)
elif (opts.stage=="plasticSideband"):
    writePlasticSideband(mywrapper,opts.outdir, opts.targetID, opts.targetZ, opts.playlist)
elif (opts.stage=="physicsSideband"):
    writePhysicsSideband(mywrapper,opts.outdir, opts.targetID, opts.targetZ, opts.playlist)
    closureTest(mywrapper,opts.outdir)
mywrapper.close()

# Making the wrapper we just created readable, writable and executable by everyone
os.system("chmod 777 %s"%(wrapper_name))

# TODO: not sure if this is needed
configstring = "" 
# Get gcc release
#gccstring = "x86_64-slc7-gcc49-opt" 
# Now add the execute command
cmd = "" 
cmd += "jobsub_submit --group minerva " #Group of experiment
#cmd += "--cmtconfig "+gccstring+" " #Setup minerva soft release built with minerva configuration

#cmd += "--OS sl7 " #Operating system #Not needed in SL7
#cmd += "--singularity-image /cvmfs/singularity.opensciencegrid.org/fermilab/fnal-wn-sl7:latest "
cmd += "--singularity-image /cvmfs/singularity.opensciencegrid.org/fermilab/fnal-wn-el9:latest "

if opts.mail:
    cmd += "-M " #this option to make decide if you want the mail or not
#cmd += "--subgroup=Nightly " #This is only for high priority jobs
#cmd += "-d "+tag_name+opts.outdir+" " #Nov7, command takes care of copying files back to dCache, so that don't have to worry about writing your own dCache commands, -David
cmd += "--resource-provides=usage_model=DEDICATED,OPPORTUNISTIC " 
cmd += "--role=Analysis " 
cmd += "--expected-lifetime 72h " 
cmd += "--memory "+str(memory)+"MB " 
cmd += "--disk "+str(disk)+"GB " #added in May20 b/c jobs going to 'held'
cmd += configstring+" " #the environments for the tunes to bee applied
cmd += "-f "+opts.outdir+"/myareatar_"+tag_name+".tar.gz " 
#cmd += "-i /cvmfs/minerva.opensciencegrid.org/minerva/software_releases/v22r1p1"+" " 
cmd += "file://"+os.environ["PWD"]+"/"+wrapper_name
print(cmd)

os.system(cmd)

print("Deleting the app area tar..... ") 
os.system("rm -rf /exp/minerva/app/users/$USER/myareatar_"+tag_name+".tar.gz")
print("Sleeping") 

time.sleep(1)
