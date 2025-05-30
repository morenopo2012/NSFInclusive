import os,sys

version="ver24"
type="_pzmu"
tune = "_2p2hrpa"
#path = "/Users/schellma/data/Nominal_"+version+"/"
path = "/pnfs/minerva/persistent/users/bashyal8/postNSFMigration/Nominal_"+version+"/"
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

file_suffix = tune+"_minervame5A6A6B6C6D6E6F6G_Nominal_ver7.root"
file_suffix = file_suffix.replace("ver7",version)
#file_suffix = "_CV_minervame6G_Nominal_ver3.root"

eff_file =  path+"EffXPurity/merged_files/merged_EffXPurity"+file_suffix

scale_file = path+"BackgroundScale/merged_BkgScaleFactor"+tune+"_minervame5A6A6B6C6D6E6F6G_Nominal_ver7.root"
scale_file = scale_file.replace("ver7",version)

event_file = path+"PTPZ/merged_files/merged_PTPZ"+file_suffix

migration_file = path+"MigrationMatrix_enu/merged_files/merged_MigrationMatrix"+file_suffix.replace(".root","_enu.root")
migration_file = migration_file.replace("_enu",type)
print "migration is ", migration_file
output_file = "XS"+type+".root"

num_iter = "4" #8 for enu q2 and 4 for pzpt
ana_type = "0" #0 for pz pt and 1 for enu q2
ana_type = code
num_iter = _iter

command = "./XS "+event_file+" "+scale_file+" "+migration_file+" "+eff_file+" "+output_file+ " "+num_iter+" "+ana_type
print command
os.system(command)
