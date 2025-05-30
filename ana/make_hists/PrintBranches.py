import	ROOT

myFile = ROOT.TFile("/pnfs/minerva/persistent/users/zdar/Merged_mc_ana_me1B_v1.61/MasterAnaDev_mc_AnaTuple_run00111000_Playlist.root")
myTree = myFile.Get("MasterAnaDev")
branches = myTree.GetListOfBranches()
leaves   = type(branches).__name__

totalSize = 0
for branch in branches:
  #print branch.GetName() + "'s size: " + str(branch.GetTotBytes()/1e6) + "MB"
  print branch.GetName()
  #print (type(branch).__name__) 
  totalSize += branch.GetTotBytes()/1e6

print "total size of file should be: " + str(totalSize) + "MB"
