//#include "include/CCQENuPlotUtils.h"
//#include "MinervaUnfold/MnvUnfold.h"
//#include "include/NukeApp.h"
//#include "include/NukeCC.h"
//#include "include/NukeUtils.h"
#include "PlotUtils/TargetUtils.h"
#include "PlotUtils/HistogramUtils.h"
#include "TParameter.h"
#include "TList.h"

bool useMCOnly = true; 

using namespace NukeCC_Ana;
double GetPOTScale(std::string playlist, int targetID ){
  FileType::t_FileType mc_type = FileType::kMC;
  if(targetID < 10 && useMCOnly ) {mc_type = FileType::kNukeOnlyMC; cout<<"I'm getting a NukeOnly POT "<<endl;}
  const double dataMCScale = NukeUtils::Get().GetDataPOT(playlist)/NukeUtils::Get().GetPlaylistPOT(playlist, mc_type);
    if(dataMCScale > 1) cout << "I expected more POT in MC than data. Make sure you ran your jobs correctly or know this is true" << endl;
    cout <<"Scale " << playlist <<"  "<< dataMCScale << endl;
    return dataMCScale;
}

double GetPOTTotalScale(std::vector<string> playlists){
    double MC_sum = 0.0;
    double Data_sum = 0.0;
    
    for(int i=0;i<playlists.size();i++){
        MC_sum += NukeUtils::Get().GetMCPOT(playlists[i]);
        Data_sum += NukeUtils::Get().GetDataPOT(playlists[i]);
    }
    cout<<"The Total MC POT for your samples is "<<MC_sum<<endl;
    cout<<"The Total Data POT for your samples is "<<Data_sum<<endl;
    cout<<"The Ratio is "<<Data_sum/MC_sum<<endl;
    
//    return Data_sum/MC_sum;
    
//    return MC_sum/Data_sum;
    
    return 1.0;
}

int PlaylistAdder(string basefilename, string target, std::vector<string> playlists)
{
  

  int targetID = stoi(target);
    //Need to open files
    const int num_playlists = playlists.size();//Full par input vector less 3 header values
    cout<<"I'm gonna add your POT and get a scale factor!"<<endl;
    double sum_scale = 1.0; //GetPOTTotalScale(playlists);
    const std::string HISTS( getenv("HISTS") );
    cout<<HISTS<<endl;
    const std::string NUKECC_TAG( getenv("NUKECC_TAG") );
    cout<<NUKECC_TAG<<endl;
    const std::string USER( getenv("USER") );
    TFile *f_input[num_playlists];
    std::string histogram_path = HISTS+"/"+NUKECC_TAG;
    cout<<"The path for your histograms is "<<histogram_path<<endl;
    cout<<"You are combining "<<num_playlists<<" playlists"<<endl;
    // read files to get histogram
    for(int i=0;i<num_playlists;i++){
        cout<<basefilename<<"  "<<playlists[i]<<endl;
        string filename = histogram_path+"/"+playlists[i]+"/"+basefilename+".root";
        cout << "Opening " << filename << endl;
        f_input[i] = new TFile(filename.c_str());
        
        if (f_input[i]->IsZombie() || f_input[i]->GetListOfKeys()->IsEmpty()){
            Error("Playlistadder","Could not get histogram ROOT file or it was empty.");
            return 1;
        }
    }
  
    string outname = "";
    outname = "/minerva/data/users/"+USER+"/NukeHists/"+NUKECC_TAG+"/AllNuME/"+basefilename+".root";
   
    TFile *f_output = new TFile(outname.c_str(),"RECREATE");
   

    //Add hist
  
    cout<<"file name "<<f_input[0]->GetName()<<endl;
    TList *file0List = f_input[0]->GetListOfKeys();
  
    TIter next(file0List);
    TObject *ob = 0;
 
    while((ob=next())){
        bool dopotnorm = true;
        string name = ob->GetName();
        string classname = ((TKey*)ob)->GetClassName();
        cout<<name<<endl;
        if(name.find("Data") != std::string::npos || name.find("data") != std::string::npos) {
	  //            cout<<"The object name contains 'data', I'm not going to scale it by POT "<<endl;
            dopotnorm = false;
            
        }
        cout<<basefilename<<endl;
        if(basefilename.find("Data") != std::string::npos || basefilename.find("data") != std::string::npos){
	  //            cout<<"The base name contains 'data', I'm not going to scale it by POT "<<endl;
            dopotnorm = false;
        }
	if(basefilename.find("Efficiency") != std::string::npos || name.find("Efficiency") != std::string::npos){
	  //	  cout<<"The base name contains 'data', I'm not going to scale it by POT "<<endl;
	  dopotnorm = false;
	}
	if(basefilename.find("Migration") != std::string::npos){
	  dopotnorm = false;
	}
        
        cout<<"for "<<name<<"  in "<<basefilename<<"  am I going to normalize by pot? "<<dopotnorm<<endl;
        
        if(classname=="TH1D"){
            cout<<"I'm summing a set of TH1Ds for "<<name<<endl;
            TH1D *hist1D = (TH1D*)f_input[0]->Get(name.c_str());
            if(dopotnorm) hist1D->Scale(GetPOTScale(playlists[0],targetID));
            for(int i=1;i<num_playlists;i++){
                TH1D *tmp1D = (TH1D*)f_input[i]->Get(name.c_str());
                if(dopotnorm) tmp1D->Scale(GetPOTScale(playlists[i],targetID));
                hist1D->Add(tmp1D);
                delete tmp1D;
            }
            if(dopotnorm) hist1D->Scale(1.0/sum_scale);
            f_output->cd();
            hist1D->Write(name.c_str());
            delete hist1D;
        }
        else if(classname=="PlotUtils::MnvH1D"){
            cout<<"I'm summing a set of PlotUtils::MnvH1D for "<<name<<endl;
            PlotUtils::MnvH1D *hist1D = (PlotUtils::MnvH1D*)f_input[0]->Get(name.c_str());
            if(dopotnorm) hist1D->Scale(GetPOTScale(playlists[0],targetID));
            for(int i=1;i<num_playlists;i++){
                PlotUtils::MnvH1D *tmp1D = (PlotUtils::MnvH1D*)f_input[i]->Get(name.c_str());
                if(dopotnorm) tmp1D->Scale(GetPOTScale(playlists[i],targetID));

		if(hist1D->HasUncorrError("Plastic_SB") || hist1D->HasUncorrError("Phys_SB") ){ //errors are getting added as binomial errors, don't want that for this. 
		  if(!hist1D->HasUncorrError("Plastic_SB")){ 
		    hist1D->AddUncorrErrorAndFillWithCV("Plastic_SB");
		  }
		  if(!hist1D->HasUncorrError("Phys_SB")){ 
		    hist1D->AddUncorrErrorAndFillWithCV("Phys_SB");
		  }
		  tmp1D->AddMissingErrorBandsAndFillWithCV(*hist1D);
		  TH1D* histCHtmp=(TH1D*) hist1D->GetUncorrError("Plastic_SB")->Clone();
		  TH1D* tmpCH=(TH1D*) tmp1D->GetUncorrError("Plastic_SB");
		  TH1D* histPhystmp=(TH1D*) hist1D->GetUncorrError("Phys_SB")->Clone();
		  TH1D* tmpPhys=(TH1D*) tmp1D->GetUncorrError("Phys_SB");

		  hist1D->Add(tmp1D);

		  TH1D* histCH= hist1D->GetUncorrError("Plastic_SB");
		  TH1D* histPhys= hist1D->GetUncorrError("Phys_SB");
		  for(int ibin=0; ibin<= hist1D->GetNbinsX()+1; ++ibin){
		    double binerr= histCHtmp->GetBinError(ibin)+tmpCH->GetBinError(ibin);
		    histCH->SetBinError(ibin, binerr);

		    binerr= histPhystmp->GetBinError(ibin)+tmpPhys->GetBinError(ibin);
		    histPhys->SetBinError(ibin, binerr);

		  }//end bin loop
		}//end Plastic or Phys SB
		else  hist1D->Add(tmp1D);
                delete tmp1D;
            }
            f_output->cd();
          if(dopotnorm)  hist1D->Scale(1.0/sum_scale);
            hist1D->Write(name.c_str());
            delete hist1D;
        }
        else if(classname=="TH2D"){
            cout<<"I'm summing a set of TH2D for "<<name<<endl;
            TH2D *hist2D = (TH2D*)f_input[0]->Get(name.c_str());
            if(dopotnorm) hist2D->Scale(GetPOTScale(playlists[0],targetID));
            for(int i=1;i<num_playlists;i++){
                TH2D *tmp2D = (TH2D*)f_input[i]->Get(name.c_str());
                if(dopotnorm) tmp2D->Scale(GetPOTScale(playlists[i],targetID));
                hist2D->Add(tmp2D);
                delete tmp2D;
            }
            f_output->cd();
          if(dopotnorm) hist2D->Scale(1.0/sum_scale);
            hist2D->Write(name.c_str());
            delete hist2D;
        }
        else if(classname=="PlotUtils::MnvH2D"){
            cout<<"I'm summing a set of PlotUtils::MnvH2D for "<<name<<endl;
            PlotUtils::MnvH2D *hist2D = (PlotUtils::MnvH2D*)f_input[0]->Get(name.c_str());
            if(dopotnorm) hist2D->Scale(GetPOTScale(playlists[0],targetID));
            for(int i=1;i<num_playlists;i++){
                PlotUtils::MnvH2D *tmp2D = (PlotUtils::MnvH2D*)f_input[i]->Get(name.c_str());
                if(dopotnorm) tmp2D->Scale(GetPOTScale(playlists[i],targetID));
                hist2D->Add(tmp2D);
                delete tmp2D;
            }
            f_output->cd();
            if(dopotnorm) hist2D->Scale(1.0/sum_scale);
            hist2D->Write(name.c_str());
            delete hist2D;
        }
        else{
            cout << "This is a type of class I don't know how to handle: "<< classname << "\tskipping. Get your life right. " << endl;
        }
    }//obj next
    ob = NULL;
    file0List = NULL;
    f_output->Close();
    return 0;
};

int main( int argc, char *argv[])
{
    //ROOT::Cintex::Cintex::Enable();
    TH1::AddDirectory(false);
    if (argc==1){
        std::cout<<"-----------------------------------------------------------------------------------------------"<<std::endl;
        std::cout<<"MACROS HELP:\n\n"<<
        "PlaylistAdder outfile_root  file_genie_variation_xsection\n\n"<<
        "\t-Base string of files. Assumes files are of the form XXXXXX_minervaY.root where Y is any playlist(s).\n"<<
        "\t-Type of file we are adding. Valid choices: MuonEventSelection, EffPurity, BackgroundTemples, Migration\n"<<
        "\t-MCOnly sample. Valid choices: 0 = no special mconly sample (2p2h for instance), 1 = yes\n"<<
        "\t-List of playlists. Example would be minerva1 minerva7 minerva 9. This will add 1,7,9 together.\n"
        "\t*********************************************\n"<<
        
        std::cout<<"-----------------------------------------------------------------------------------------------"<<std::endl;
        
        return 0;
    }

    
    //! Default parameters
    std::vector<std::string> par;
    std::vector<std::string> playlists;
    par.push_back("PlaylistAdder");
    par.push_back(argv[1]);
    par.push_back(argv[2]);
    //! Set user parameters
    for( int i=3; i<argc; ++i){
        par.push_back(argv[i]);
        playlists.push_back(argv[i]);
    }

    for( unsigned int i=0; i<par.size(); ++i)
        std::cout<<"Parameter "<< i << ": " << par[i] << std::endl;
    std::cout<<"I'm about to try and run this"<<std::endl;
    return PlaylistAdder(par[1],par[2], playlists);
}
