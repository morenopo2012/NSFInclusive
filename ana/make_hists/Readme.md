All the Files  generated for histograms, no plotting were stored at this location:

/exp/minerva/data/users/omorenop/FateStudiesPlots

Feb 2025
Changes you need to make to switch to different modes in the MEC_FateWeight script.

   ECAL mode:
     Fiducial cut: 
        NukeCC_Cuts.cxx
        NukeCC_Cuts::isFiducialXYZ 
     playlist:
        NSFNukeCCInclusive/ana/include/playlists (change the files to use)



How to run Tracker/DSCAL Samples for the ingredients to extract the cross section.


runEventLoop
	File required: runEventLoop_MEC_FateWeight.cxx



