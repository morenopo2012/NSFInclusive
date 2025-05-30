# Sample command: . gridSubmit.sh false -14 MnvTunev1_RHC_neg14

neutrinoMode=$1
echo $neutrinoMode

echo "Let's get started!"
#outdir=/pnfs/minerva/persistent/users/zdar/default_analysis_loc/$filename

if $neutrinoMode; then 
	echo "Neutrino Mode"
	#for STAGE in "eventLoop" "migration" "efficiency" "plasticSideband" "physicsSideband" ;do
	#for STAGE in "migration" "efficiency" "eventLoop" "physicsSideband"; do
	for STAGE in "migration";do
           for PLAYLIST in  "minervame1B" "minervame1C" "minervame1D" "minervame1E" "minervame1F" "minervame1G" "minervame1M" "minervame1N" "minervame1P";do
           #for PLAYLIST in "minervame1A" "minervame1D" "minervame1E" "minervame1F" "minervame1G" "minervame1M" "minervame1N"; do
           #for PLAYLIST in "minervame1E" "minervame1F"; do
              #python SubmitJobsToGrid.py --stage $STAGE --playlist $PLAYLIST --basedir /exp/minerva/app/users/zdar/MAT_AL9/ --targetID 3  --targetZ 82 
              #echo "***JUST FINISHED SUBMITTING $PLAYLIST for $STAGE"
              python SubmitJobsToGrid.py --stage $STAGE --playlist $PLAYLIST --basedir /exp/minerva/app/users/zdar/MAT_AL9/ --targetID 3  --targetZ 6 
              echo "***JUST FINISHED SUBMITTING $PLAYLIST for $STAGE"
              #python SubmitJobsToGrid.py --stage $STAGE --playlist $PLAYLIST --basedir /exp/minerva/app/users/zdar/MAT_AL9/ --targetID 3  --targetZ 26 
              #echo "***JUST FINISHED SUBMITTING $PLAYLIST for $STAGE"
              python SubmitJobsToGrid.py --stage $STAGE --playlist $PLAYLIST --basedir /exp/minerva/app/users/zdar/MAT_AL9/ --targetID 14  --targetZ 82 
              echo "***JUST FINISHED SUBMITTING $PLAYLIST for $STAGE"
           done
        done

   # RUNNING OVER RHC PLAYLISTS
else
	echo "Anti Neutrino mode"
        for STAGE in "eventLoop";do
           for PLAYLIST in "minervame6B" "minervame5A" "minervame6A" "minervame6C" "minervame6D" "minervame6E" "minervame6F" "minervame6G" "minervame6H" "minervame6I" "minervame6J";do
           #for PLAYLIST in "minervame6F"; do
              #python SubmitJobsToGrid.py --stage $STAGE --playlist $PLAYLIST --basedir /exp/minerva/app/users/zdar/MAT_AL9/ --targetID 3  --targetZ 82 
              #echo "***JUST FINISHED SUBMITTING $PLAYLIST for $STAGE"
              python SubmitJobsToGrid.py --stage $STAGE --playlist $PLAYLIST --basedir /exp/minerva/app/users/zdar/MAT_AL9/ --targetID 3  --targetZ 26 
              echo "***JUST FINISHED SUBMITTING $PLAYLIST for $STAGE"
              python SubmitJobsToGrid.py --stage $STAGE --playlist $PLAYLIST --basedir /exp/minerva/app/users/zdar/MAT_AL9/ --targetID 14  --targetZ 82 
              echo "***JUST FINISHED SUBMITTING $PLAYLIST for $STAGE"
              #python SubmitJobsToGrid.py --stage $STAGE --playlist $PLAYLIST --basedir /exp/minerva/app/users/zdar/MAT_AL9/ --targetID 3  --targetZ 6 
              #echo "***JUST FINISHED SUBMITTING $PLAYLIST for $STAGE"
           done
        done

fi

rm -rf eventLoop_minervame*sh
rm -rf migration_minervame*sh
rm -rf efficiency_minervame*sh 
rm -rf plastic*sh 
rm -rf physics*sh 
