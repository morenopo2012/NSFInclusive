import os,sys

variables = ["nu_hadQ2","W_hadQ2","hadQ2_x","hadQ2_y","x_y"]
grid = ["5 4","5 4","4 3","5 4","4 3"]

for i,v in enumerate(variables):
    os.system("mkdir -p Inclusive_Gen/%s"%(v))
    #cmd = "make;./plot_nu_xsec2_inclusive_general /minerva/data/users/drut1186/3D_ME_Inclusive/ProtoType_v6/Inclusive_CV/CrossSection_per_nucleon_iterations_10_CombinedPlaylists.root_%s.root %s"%(v,v)
    cmd = "make;./error_plot_grid_generic /minerva/data/users/drut1186/3D_ME_Inclusive/ProtoType_v6/Inclusive_CV/CrossSection_per_nucleon_iterations_10_CombinedPlaylists.root_%s.root %s %s %s 800 500 %s 800 500 CrossSection"%(v,v.split("_")[0],v.split("_")[1],grid[i],grid[i])
    os.system(cmd)
    os.system("mv *.png Inclusive_Gen/%s"%(v))
