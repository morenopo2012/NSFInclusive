import os,sys
import time

playlists=["minervame1A"]
#playlists=["minervame6B"]
#stages=["RecoilEnergy","PTPZ","MigrationMatrix","EffXPurity"]
#stages=["PTPZ","RecoilEnergy","EffXPurity"]
#stages = ["RecoilEnergy"]
stages=["PTPZ"]
migration = "enu" #"pzmu" or "enu"
#stages = ["PTPZ_SideBand","MigrationMatrix"]
#stages = ["MigrationMatrix"]
#tunes = ["--piontune","--2p2h","--rpa","--2p2h --rpa","--rpa --piontune","--CV","--CV2"]
#tunes = ["--CV"]
#tunes = ["--2p2h --piontune","--2p2h --rpa","--NA"]
tunes = ["NA"]
#outdir = "/pnfs/minerva/persistent/users/bashyal8/postJulyCollab/postMinerba/Nominal/"
#outdir = "/pnfs/minerva/persistent/users/bashyal8/postJulyCollab/post_Bug_Fix/newpZbinsSystematics/Optimized/"
outdir = "/pnfs/minerva/persistent/users/zdar/histsLoc/"
#Also use the same tarfile for >1 loop thingy....
time_stamp = int(time.time())
flag_tag = "Nominal_"+str(time_stamp)
#tarfilename="myareatar_optimized_1573147512piontunepiontune_.tar.gz"
tarfilename=""
for i in range(0,len(playlists)):
    for j in range(0,len(tunes)):
        for k in range(0,len(stages)):
            cmd=""
            if i==0 and \
               j==0 and \
               k==0:
                tune__tag=tunes[j].replace("--","")
                #in case of multiple tags, space needs to be removed.
                tune_tag = tune__tag.replace(" ","")
                tag_name=tune_tag
                tarfilename = "myareatar_"+tag_name+".tar.gz"
                #cmd = "python SubmitPlaylistToGrid_antiNu.py --outdir "+outdir+" --playlist "+playlists[i]+" "+tunes[j]+" --tag "+tag_name+" --stage"+" "+stages[k]+" "+"--sample Signal --notimestamp  --samplemigration "+migration
                cmd = "python SubmitPlaylistToGrid_antiNu.py --outdir "+outdir+" --playlist "+playlists[i]+" "+tunes[j]+" --tag "+tag_name+" --stage"+" "+stages[k]+" "+"--sample Signal --notimestamp  --samplemigration "+migration
		os.system(cmd)
		#once the tarfile is copied...go to the directory and find the tar file...
		_list = os.listdir(outdir)
		for things in _list:
		    if tag_name in things:
		        tarfilename=things
			break
			
		
            else:
                tune__tag = tunes[j].replace("--","")
                tune_tag = tune__tag.replace(" ","")
                tag_name = tune_tag
                cmd = "python SubmitPlaylistToGrid_antiNu.py --outdir "+ outdir+" --playlist "+playlists[i]+" "+tunes[j]+" --tag "+tag_name+" --stage"+" "+stages[k]+" "+"--sample Signal --samplemigration "+migration+" --sametar --tarfilename "+tarfilename
            	os.system(cmd)
        #print tag_name
        #print cmd
	#print tarfilename
