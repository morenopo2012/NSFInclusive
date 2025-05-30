#!/bin/sh
source /cvmfs/minerva.opensciencegrid.org/minerva/software_releases/v22r1p1/setup.sh
cd $CONDOR_DIR_INPUT
tar -xvzf myareatar_tag_default_1602902521.tar.gz
cd Tools/ProductionScriptsLite/cmt
cmt config
source setup.sh
cd $CONDOR_DIR_INPUT
cd Ana/NSFNukeCCInclusive/cmt/ 
cmt config
source setup.sh
export LD_LIBRARY_PATH=$NSFNUKECCINCLUSIVEROOT/NUKECCSRC/ana_common/src:$LD_LIBRARY_PATH
cd $_CONDOR_SCRATCH_DIR
mkdir playlists
cp $CONDOR_DIR_INPUT/Ana/NSFNukeCCInclusive/ana/include/playlists/*.txt $_CONDOR_SCRATCH_DIR/playlists
echo $LD_LIBRARY_PATH
ls /usr/lib64
$NSFNUKECCINCLUSIVEROOT/ana/make_hists/test/runEventLoop $CONDOR_DIR_HISTS/ 14 82 minervame1A