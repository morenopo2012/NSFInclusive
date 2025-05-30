import os,sys

histograms = ["h_pzptrec_data_nobck_unfold_effcor_cross_section",
              "h_pzptrec_data_nobck_unfold_effcor",
              "h_pzptrec_data_nobck_unfold",
              "h_pzptrec_data_nobck",
              "h_pzptrec_data_Alltrack",
              "h_pzptrec_mc_nobck_unfold_effcor_cross_section",
              "h_pzptrec_mc_nobck_unfold_effcor",
              "h_pzptrec_mc_nobck_unfold",
              "h_pzptrec_mc_nobck",
              "h_pzptrec_mc_Alltrack"]

target = ["3DError/xsec",
          "3DError/eff",
          "3DError/unf",
          "3DError/nobck",
          "3DError/sel",
          "3DError_MC/xsec",
          "3DError_MC/eff",
          "3DError_MC/unf",
          "3DError_MC/nobck",
          "3DError_MC/sel"]


indir = sys.argv[1]
outdir = sys.argv[2]
for i,h in enumerate(histograms):
    for j in range(1,11):
        if(j!=10): continue
        os.system("mkdir -p %s/%s/iteration_%d"%(outdir,target[i],j))
        os.system("rm *.png")
        os.system("rm *.eps")
        os.system("rm *.C")
        cmd = "./error_plot_grid_3D %s/CrossSection_per_nucleon_3D_pzptreco_iterations_%s_CombinedPlaylists.root_big3d.root %s"%(sys.argv[1],j,h)
        os.system(cmd)
        os.system("mv *.png %s/%s/iteration_%d"%(outdir,target[i],j))
        os.system("mv *.eps %s/%s/iteration_%d"%(outdir,target[i],j))
        os.system("mv *.C %s/%s/iteration_%d"%(outdir,target[i],j))
    
              
              
