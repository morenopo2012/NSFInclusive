import os,sys,time
from optparse import OptionParser
import datetime

###########################################################
#Author: Amit Bashyal
#Script: JobSubmission wrapper creater for CCQENUNSF 
#
#
#
####################################################

def createTarball(outdir,tag):
    #tooldir = os.environ["CCQENUNSFROOT"]
    tooldir = os.environ["NSFNUKECCINCLUSIVEROOT"]
    basedir = "%s/../../"%(tooldir)
    found = os.path.isfile("%s/myareatar_%s.tar.gz"%(outdir,tag))
    if(not found):
        cmd = "tar -czf /minerva/app/users/$USER/myareatar_%s.tar.gz %s"%(tag,basedir)
        print "Making tar",cmd
        os.system(cmd)
        cmd2 = "cp /minerva/app/users/$USER/myareatar_%s.tar.gz %s/"%(tag,outdir)
        print "Copying tar",cmd2
        os.system(cmd2)

def writeTarballProceedure(mywrapper,tag):
    #some uh magic to get the my_ccqenu setup using the submission setup... otherwise we don't know which one to setup on the grid node...
    myccqenu = os.environ["NSFNUKECCINCLUSIVEROOT"]
    #myversion = myccqenu.split("/")[-2]#Count from the end and skip the ana.
    print "I will be making the tarball upacking with this version"
    print "Path is",myccqenu
    mywrapper.write("source /cvmfs/minerva.opensciencegrid.org/minerva/software_releases/v22r1p1/setup.sh\n")
    mywrapper.write("cd $CONDOR_DIR_INPUT\n")
    mywrapper.write("tar -xvzf myareatar_%s.tar.gz\n"%tag)
    mywrapper.write("cd Tools/ProductionScriptsLite/cmt\n")
    mywrapper.write("cmt config\n")
    mywrapper.write("source setup.sh\n")
    mywrapper.write("cd $CONDOR_DIR_INPUT\n")
    mywrapper.write("cd Ana/NSFNukeCCInclusive/cmt/ \n")
    mywrapper.write("cmt config\n")
    mywrapper.write("source setup.sh\n")
    mywrapper.write("export LD_LIBRARY_PATH=$NSFNUKECCINCLUSIVEROOT/NUKECCSRC/ana_common/src:$LD_LIBRARY_PATH\n")
    #need to kill and relink our main so
    mywrapper.write("cd $_CONDOR_SCRATCH_DIR\n")#go back home...
    mywrapper.write("mkdir playlists\n")#go back home...
    mywrapper.write("cp $CONDOR_DIR_INPUT/Ana/NSFNukeCCInclusive/ana/include/playlists/*.txt $_CONDOR_SCRATCH_DIR/playlists\n")#go back home...
    mywrapper.write("echo $LD_LIBRARY_PATH\n")
    mywrapper.write("ls /usr/lib64\n")


#valid stages for the anti neutrinos
valid_stages=["PTPZ"]
avail_playlists=["minervame1A"]

_user_ = os.getenv("USER") #make sure the files goes to the user's persistent working place......
usage = "usage: %prog[opts]"
parser = OptionParser(usage=usage)
parser.add_option('--outdir', dest='outdir', help='Directory to write output to', default = "/pnfs/minerva/persistent/users/"+_user_+"/default_analysis_loc/")
parser.add_option('--stage', dest='stage', help='Process type', default="NONE")
parser.add_option('--sample', dest='sample', help='Sample type', default="NONE")
parser.add_option('--playlist', dest='playlist', help='Playlist type', default="NONE")
##########################################################################
#Options for making tar files....Basically you can make tarfiles 
#######################################################################
parser.add_option('--tag', dest='tag', help="Tag your release",default="tag_")
parser.add_option('--mail',dest='mail',help="Want mail after job completion or failure",default=False,action="store_true")
parser.add_option('--sametar',dest='sametar',help="Recycle the same tar file for jobsubmission",default=False,action="store_true")
parser.add_option('--tarfilename',dest='tarfilename',help='Name of the tarfile you want to use',default="NA")
parser.add_option('--notimestamp',dest='notimestamp',help='Flag to TURN OFF time stamp in the tag',default=False,action="store_true")
parser.add_option('--targetID',dest='targetID',help='Which target to analyze',action="store")
parser.add_option('--targetZ',dest='targetZ',help='Which material to analyze within chosen target',action="store")
parser.add_option('--CV',dest='CV',help="Do piontune+CCQERPA+2p2hTune-->MnvGENIEv1",default=False,action="store_true")
parser.add_option('--CV2',dest='CV2',help="Do piontune+CCQERPA+2p2hTune+MnvLowQ2-->MnvGENIEv2",default=False,action="store_true")
############################################################################
#A note on tag...
#Always put tunes_timestamp 
#Anything else should be put by invoking "--tag"
############################################################################

#need to put the option of interactive as well...
parser.add_option('--interactive',dest='interactive',
                  help='Run it interactively or on grid',default=False,action="store_true")
#parser.add_option('--stage', dest='stage', help='Process type', default="NONE")

(opts,args) = parser.parse_args()



################MAIL#############################
if opts.mail:
    print "############################################"
    print "You have opted to Get mails once the jobs are finished of failed."
    print "#############################################"
#################MAIL##############################

###############INTERACTIVE######################
if opts.interactive:
    print "###################################################"
    print "You have chosen to run this interactively ."
    print "###################################################"

##############################################
#here I want to put the tag scheme...
##############################################
tag_name = str(opts.tag)
time_stamp = int(time.time())
if tag_name=="tag_":
    print "YOU DIDNT SPECIFY ANY ANY TAG...SO I WILL USE MY OWN TAGGING SCHEME"
    if opts.CV:
        tag_name.replace("piontune_","")
        tag_name.replace("rpa_","")
        tag_name.replace("resrpanieves_","")
        tag_name.replace("lowq2_","")
        tag_name += "CV_"

    elif opts.CV2:
        tag_name.replace("piontune_","")
        tag_name.replace("rpa_","")
        tag_name.replace("resrpanieves_","")
        tag_name.replace("lowq2_","")
        tag_name += "CV2_"
    else:
        tag_name += "default_"

print "Is everything fine till here*********"
if not opts.notimestamp:
    tag_name += str(time_stamp)

print "Tag for this version is ",tag_name
print "******************************************************"

##############################################
#here I want to put the tag scheme...
##############################################

print "************************************************"
temp_dir = " "
if not opts.interactive:
    temp_dir = "$_CONDOR_SCRATCH_DIR/"


output_dir=" "
if not opts.interactive:
    output_dir = "$CONDOR_DIR_HISTS/"

input_dir=" "
if not opts.interactive:
    input_dir = "$CONDOR_DIR_INPUT/"

if(not os.path.isdir(opts.outdir)):
    os.makedirs(opts.outdir)

#TODO: Figure this out using /usr/bin/time -v interactively
memory = 2000


if opts.stage not in valid_stages:
    print opts.stage,"Selected stage is not valid. Here are the valid options",valid_stages
    sys.exit()

#Create wrapper
wrapper_name = "%s_%s_wrapper_%s.sh"%(opts.stage,opts.playlist,tag_name)

my_wrapper = open(wrapper_name,"w")
my_wrapper.write("#!/bin/sh\n")
if (opts.sametar==False):
    createTarball(opts.outdir,tag_name)
else:
    tarname = str(opts.tarfilename)
    if not os.path.exists(opts.outdir+opts.tarfilename):
        print "Tar File "+opts.outdir+opts.tarfilename+" doesn't Exist!"
        sys.exit()
    #change the tag to the current one...
    cmd="cp "+opts.outdir+opts.tarfilename+" "+opts.outdir+"myareatar_"+tag_name+".tar.gz"
    os.system(cmd)

#this will unpack the tarball we just made above
writeTarballProceedure(my_wrapper,tag_name)

PTPZ_output = "PTPZ_"+tag_name+tag_name+str(opts.playlist)+".root"

#Now the PTPZ distributions
if(opts.stage=="PTPZ"):                              
    my_wrapper.write("$NSFNUKECCINCLUSIVEROOT/ana/make_hists/test/runEventLoop "+str(output_dir)+" "+str(opts.targetID)+" "+str(opts.targetZ)+" "+str(opts.playlist))
my_wrapper.close()


os.system("chmod 777 %s"%(wrapper_name))

configstring = ""
if(opts.CV):
    configstring +=" -e DOPIONTUNE "
    os.environ["DOPIONTUNE"]="1"

    configstring +=" -e DOCCQERPA "
    os.environ["DOCCQERPA"]="1"

    configstring += " -e DO2P2HTUNE "
    os.environ["DO2P2HTUNE"]="1"

if(opts.CV2):
    configstring +=" -e DOPIONTUNE "
    os.environ["DOPIONTUNE"]="1"

    configstring +=" -e DOCCQERPA "
    os.environ["DOCCQERPA"]="1"

    configstring += " -e DO2P2HTUNE "
    os.environ["DO2P2HTUNE"]="1"

    configstring += "-e DOMINERVALOWQ2"
    os.environ["DOMINERVALOWQ2"]="1"


mnvRelease = os.environ["MINERVA_RELEASE"]
mnvRelease = mnvRelease[ :3]
print mnvRelease
if(mnvRelease == "v10" ):
    gccstring= "x86_64-slc6-gcc44-opt" 
if(mnvRelease == "v20" or mnvRelease=="v21"):
    gccstring= "x86_64-slc6-gcc49-opt"
if(mnvRelease == "v22"):
    gccstring = "x86_64-slc7-gcc49-opt"

cmd = ""
cmd += "jobsub_submit --group minerva " #Group of experiment
cmd += "--cmtconfig "+gccstring+" " #Setup minerva soft release built with minerva configuration
#if self.osversion == "SL6" or self.osversion == "sl6":
#cmd += " --append_condor_requirements='(TARGET.HAS_SINGULARITY=?=true)' "
#cmd += " --lines '+SingularityImage=\\\"/cvmfs/singularity.opensciencegrid.org/fermilab/fnal-wn-sl6:latest\\\"' "
#else:
#self.jobsub_flags += " --OS=%s" % self.osversion

cmd += "--OS sl7 " #Operating system #Not needed in SL7

if opts.mail:
    cmd += "-M " #this option to make decide if you want the mail or not
#cmd += "--OS=SL6 " #Operating system
#cmd += "--subgroup=Nightly " #This is only for high priority jobs
cmd += "--resource-provides=usage_model=DEDICATED,OPPORTUNISTIC "
cmd += "--role=Analysis "
cmd += "--expected-lifetime 12h "
#cmd += "--expected-lifetime 8h "
cmd += "--memory "+str(memory)+"MB "
cmd += configstring+" " #the environments for the tunes to bee applied
cmd += "-f "+opts.outdir+"/myareatar_"+tag_name+".tar.gz "
cmd += "-d HISTS "+str(opts.outdir)+" " #writable directory and path were outputs will be copied after the completion of the job
cmd += "-i /cvmfs/minerva.opensciencegrid.org/minerva/software_releases/"+os.environ["MINERVA_RELEASE"]+" "
cmd += "file://"+os.environ["PWD"]+"/"+wrapper_name
print cmd
if opts.interactive:
    print "Execute the command to run interactively "
    print "Command is: ", cmd
else:
    os.system(cmd) #TODO: Put me back
    #print "Ooops, you would have submitted the test script to the grid!  Edit line 183 if you really want to do that."

print "Deleting the app area tar..... "
os.system("rm -rf /minerva/app/users/$USER/myareatar_"+tag_name+".tar.gz")
print "Sleeping"

time.sleep(10)
