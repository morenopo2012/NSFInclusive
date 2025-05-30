import ROOT

myFile = ROOT.TFile("/pnfs/minerva/persistent/users/vansari/LatestMergedImages/Merged_LatticeImages_ME1LData/MasterAnaDev_data_AnaTuple_run00019168_Playlist.root")
#myFile = ROOT.TFile("")
#myTree = myFile.Get("MasterAnaDev")
myTree = myFile.Get("MasterAnaDev")
branches = myTree.GetListOfBranches()

totalSize = 0
for branch in branches:
  print branch.GetName() + "'s size: " + str(branch.GetTotBytes()/1e6) + "MB"
  totalSize += branch.GetTotBytes()/1e6

print "total size of file should be: " + str(totalSize) + "MB" 
