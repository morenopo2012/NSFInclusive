import os,sys

version="ver24"
type="_pzmu"
tune = "_2p2hrpa"
#path = "/Users/schellma/data/Nominal_"+version+"/"
path1 = "/pnfs/minerva/persistent/users/bashyal8/postNSFMigration/Nominal_"+version+"/"
path = "/minerva/app/users/zdar/cmtuser/Minerva_v22r1p1_MADNew/Ana/NSFNukeCCInclusive/ana/make_hists/XSTest/"
if len(sys.argv) < 2:
  print " need to enter pz or enu"
  sys.exit(1)
  
if sys.argv[1] == "pz":
  type = "_pzmu"
  code = "0"
  _iter = "4"
else:
  type = "_enu"
  code = "1"
  _iter = "8"

targetID = "14"
targetZ  = "82"
file_suffix = tune+"_minervame5A6A6B6C6D6E6F6G_Nominal_ver7.root"
file_suffix = file_suffix.replace("ver7",version)
#file_suffix = "_CV_minervame6G_Nominal_ver3.root"

eff_file =  path+"Hists_Efficiency_t"+targetID+"_z"+targetZ+"_Nu_minervame1D.root"

scale_file = path1+"BackgroundScale/merged_BkgScaleFactor"+tune+"_minervame5A6A6B6C6D6E6F6G_Nominal_ver7.root"
scale_file = scale_file.replace("ver7",version)

event_file = path+"Hists_EventSelection_t"+targetID+"_z"+targetZ+"_Nu_minervame1D.root"

migration_file = path+"Hists_Migration_t"+targetID+"_z"+targetZ+"_Nu_minervame1D.root"
#migration_file = path+"MigrationMatrix_enu/merged_files/merged_MigrationMatrix"+file_suffix.replace(".root","_enu.root")
#migration_file = migration_file.replace("_enu",type)
print "migration is ", migration_file
output_file = "XS"+type+"_t"+targetID+"_z"+targetZ+".root"
num_iter = "4" #8 for enu q2 and 4 for pzpt
ana_type = "0" #0 for pz pt and 1 for enu q2
ana_type = code
num_iter = _iter

#command = "./XS "+event_file+" "+scale_file+" "+migration_file+" "+eff_file+" "+output_file+ " "+num_iter+" "+ana_type
command = "./XS "+event_file+" "+scale_file+" "+migration_file+" "+eff_file+" "+output_file+ " "+num_iter+" "+ana_type+" "+targetID+" "+targetZ
print command
os.system(command)
