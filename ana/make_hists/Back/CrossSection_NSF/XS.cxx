

#include "include/NukeCCUtilsNSF.h"
#include "include/NukeCC_Cuts.h"
#include "MinervaUnfold/MnvUnfold.h"
#include "PlotUtils/TargetUtils.h"
#include "PlotUtils/FluxReweighter.h"
#include "include/GlobalIncludes.h"
#include "TROOT.h"
#include "TFile.h"
#include "TParameter.h"

#define NEWCOV 1
#define DEBUG 1

using namespace NUKECC_ANA;

//Declare the functions here.....


void RebinFluxHist(TH1D* h_flux, TH1D *&h_rebinned_flux){
   cout << "TRACE: RebinFluxHist " << h_flux->GetName() << endl;
  
  //strategy is to recale orig by bin width (undo bin width normalization) then combine bins and then rescale by binwidth again
  TH1D* scaler = (TH1D*)h_flux->Clone("fluxscaler");
  TH1D* flux_cv = (TH1D*)h_flux->Clone("fluxcvtomod");
  for(int i=1;i<scaler->GetNbinsX();i++) scaler->SetBinContent(i,scaler->GetBinWidth(i));
  flux_cv->Multiply(scaler);//undid bin width normalization
  vector<double>rebinned_flux_bin_edges;
  for(int i=1;i<h_rebinned_flux->GetNbinsX()+2;i++){//need low edge of overflow (high edge of last bin)
    rebinned_flux_bin_edges.push_back(h_rebinned_flux->GetBinLowEdge(i));
  }
  h_rebinned_flux = (TH1D*)flux_cv->Rebin(rebinned_flux_bin_edges.size()-1,"Fluxrebinned",&rebinned_flux_bin_edges[0]);
  h_rebinned_flux->Scale(1.0,"width");//And redo bin width norm
  //h_rebinned_flux->SaveAs("FluxRebinned.root");
  
}
MnvH2D* GetTotalComponent(TFile *eventfile,std::string hist_type,std::string category){
  MnvH2D *_mc = (MnvH2D*)eventfile->Get((hist_type+"bin_0_"+category).c_str());
  //std::cout<<hist_type<<"bin_0"<<category<<std::endl;
  //std::cout<<_mc->GetNbinsX()<<std::endl;
  MnvH2D *tot_mc = _mc->Clone((hist_type+category).c_str());
  tot_mc->Reset();
  tot_mc->SetDirectory(0);
  for(int i=0;i<14;i++){
    MnvH2D *temp_mc = (MnvH2D*)eventfile->Get((hist_type+"bin_"+std::to_string(i)+"_"+category).c_str());
    tot_mc->Add(temp_mc);
  
  }
  return tot_mc;
}


std::map<std::string,MnvH2D*>GetTotalSignalComponent(TFile *eventfile,std::string hist_type){
   const char *sig_list[5] = {"qelike_2p2h","qelike_res","qelike_dis","qelike_qe","qelike"};
   std::map<std::string,MnvH2D*> sig_dict;
   for(int i=0;i<5;i++){
     sig_dict[sig_list[i]] = GetTotalComponent(eventfile,hist_type,sig_list[i]);
     
   }
   return sig_dict;
}

std::map<std::string,MnvH2D*>GetTotalBackgroundComponent(TFile *eventfile,std::string hist_type){
   const char *sig_list[6] = {"qelikenot_2p2h","qelikenot_res","qelikenot_dis","qelikenot_qe","qelikenot_coh","qelikenot"};
   std::map<std::string,MnvH2D*> sig_dict;
   for(int i=0;i<6;i++){
     sig_dict[sig_list[i]] = GetTotalComponent(eventfile,hist_type,sig_list[i]);
     
   }
   return sig_dict;
}

//Get The flux rebinned in the bins of enu.....
MnvH2D* GetEnuFluxHists(MnvH2D *enu_hist,MnvH1D *h_flux){
    MnvH2D *h_flux_normalization = (MnvH2D*)enu_hist->Clone("h_flux_normalization");
    h_flux_normalization->ClearAllErrorBands();
    h_flux_normalization->Reset();
    
    //Calculate the flux integrated value for each universe and set those values in bins of pt-pz
    //Integrate from 0.0 to 100.0 GeV (JO:we're removing the neutrino energy cut)
    double e_min = 0.0;//GeV
    double e_max = 100.;//GeV
    int b_min = h_flux->FindBin( e_min );
    int b_max = h_flux->FindBin( e_max );
    
   // double flux_cv = h_flux->Integral( b_min, b_max, "width" );
    
    //const int lowBin = 0;
    const int highBinX = h_flux_normalization->GetNbinsX()+1;
    const int highBinY = h_flux_normalization->GetNbinsY()+1;
    //const int highBin = h_flux_normalization->GetBin( highBinX, highBinY );


    cout << "Doing this by energy" << endl;
    //strategy is get 1D th1d and fill up h_flux_normalization
    TH1D* h_rebinned_flux = new TH1D(h_flux_normalization->ProjectionX("proj_x",0,-1)->GetCVHistoWithStatError());
    TH1D* tmphist = new TH1D(h_flux->GetCVHistoWithStatError());
    RebinFluxHist(tmphist, h_rebinned_flux);
    for(int i=0;i<highBinX;i++){
      for(int j=0;j<highBinY;j++){
	h_flux_normalization->SetBinContent(i,j,h_rebinned_flux->GetBinContent(i));
      }
    }
    cout << "Now for verts" << endl;
    //Do the same with the vertical error bands
    std::vector<std::string> vertNames = h_flux->GetVertErrorBandNames();
    for(unsigned int k=0; k<vertNames.size(); ++k ){
      MnvVertErrorBand *errBand = h_flux->GetVertErrorBand( vertNames[k] );
      //const int universes = errBand->GetNHists();
      int universes = errBand->GetNHists();
      //const int nuniverses = errBand->GetNHists();
      //int universes = 0;
      //if(nuniverses > 100) universes = 100;
      //else universes = nuniverses;  
      if(vertNames[k]=="Flux") universes = 100;//cout << "Flux Universes = " << universes << endl;
      std::vector<TH2D*> vert_hists;
      if(vertNames[k].find("_BeamFocus")!=-1||vertNames[k].find("ppfx1_Total")!=-1) continue;
      cout << "Working on " << vertNames[k] << endl;
      for(int u=0; u< universes; ++u)
	{
	  TH2D *h_vert = new TH2D( (TH2D)*h_flux_normalization);
	  h_vert->SetName( Form("h_vert_%s_%i", vertNames[k].c_str(), u) );
	  
	  //strategy is get 1D th1d and fill up h_flux_normalization
	  TH1D* h_rebinned_flux = new TH1D(h_flux_normalization->ProjectionX("proj_x",0,-1)->GetCVHistoWithStatError());
	  TH1D* tmphist = new TH1D(*errBand->GetHist( u ));
	  RebinFluxHist(tmphist, h_rebinned_flux);
	  for(int i=0;i<highBinX;i++){
	    for(int j=0;j<highBinY;j++){
	      h_flux_normalization->SetBinContent(i,j,h_rebinned_flux->GetBinContent(i));
	    }
	  }
	  vert_hists.push_back( h_vert );
	}
      h_flux_normalization->AddVertErrorBand( vertNames[k], vert_hists );
      
      
      
      //cleaning
      for(std::vector<TH2D*>::iterator itHist = vert_hists.begin(); itHist != vert_hists.end(); ++itHist)
	delete *itHist;
  }     

    return h_flux_normalization;
}

std::map<std::string,MnvH2D*>GetFitBinBkgSubtraction(TFile *f,TFile *scale_file,std::string hist_type){

   std::map<std::string,MnvH2D*> subtracted_components;
   MnvH2D *_template = (MnvH2D*) f->Get((hist_type+"bin_0_mc").c_str());
   MnvH2D *scaled_databkgsub = _template->Clone((hist_type+"databkgsub").c_str());
   scaled_databkgsub->Reset();
   MnvH2D *scaled_bkg =scaled_databkgsub->Clone((hist_type+"qelikenotscaled").c_str());
   std::vector<std::string>vertnames = scaled_bkg->GetVertErrorBandNames();
   //TVector2 *pot = (TVector2*)f->Get("pot");
   //double weight = pot->X()/pot->Y();
   //for(size_t i = 0;i!=vertnames.size();i++)std::cout<<vertnames[i]<<std::endl;
   for(int i =0;i<14;i++){
      MnvH2D *qelikenot = (MnvH2D*) f->Get((hist_type+"bin_"+std::to_string(i)+"_qelikenot").c_str());
      MnvH2D *data = (MnvH2D*) f->Get((hist_type+"bin_"+std::to_string(i)+"_data").c_str());
      MnvH2D *scale = (MnvH2D*)scale_file->Get((hist_type+"scale_"+std::to_string(i)).c_str());
      
      MnvH2D *scaled_qelikenot = (MnvH2D*)qelikenot->Clone((hist_type+"bin_"+to_string(i)+"scaledqelikenot").c_str());
      scaled_qelikenot->Multiply(qelikenot,scale,1.0,1.0);
      scaled_bkg->Add(scaled_qelikenot);
      MnvH2D *temp_data = data->Clone();
      
      temp_data->AddMissingErrorBandsAndFillWithCV(*qelikenot);
      MnvH2D *scaled_data = temp_data->Clone();
      //scaled_data->Reset();
      
      scaled_data->Multiply(temp_data,scale,1.0,1.0);
      
      scaled_databkgsub->Add(scaled_data);
      
      delete temp_data; delete scaled_data; //delete bkgsub_data;
      delete qelikenot; delete data; delete scale;
   }
   
   subtracted_components["scaled_bkg"] = scaled_bkg;
   subtracted_components["bkgsub_data"] = scaled_databkgsub;
   return subtracted_components;
}


void CorrectByEfficiency( MnvH2D *&h_effcor, MnvH2D* h_unfolded, MnvH2D* h_effhist_num, MnvH2D* h_effhist_dem )
{
  //create efficiency corrected histo
  h_effcor = (MnvH2D*)h_unfolded->Clone( Form( "%s_effcor", h_unfolded->GetName() ) );


  MnvH2D *h_eff_cor_temp;
    h_effhist_dem->AddMissingErrorBandsAndFillWithCV( *h_unfolded );
    h_eff_cor_temp = (MnvH2D*)h_effhist_num->Clone(Form("temp_eff_cor_%s",h_unfolded->GetName()));
    h_eff_cor_temp->Divide(h_effhist_num,h_effhist_dem,1.0,1.0,"B");
  
  //efficiency doesn't have lateral error bands
  //so add to efficiency the lateral error bands found in unfolded data 
  //and use the efficiency central value to fill the universes 
  //This is for the efficiency denominator which is NOT affected by this... So CV is fine.
  //h_eff_cor_temp[kQELike]->AddMissingErrorBandsAndFillWithCV( *h_unfolded );
    
  //  h_effhist[kQELike]->Scale(0.2);

  //divide by efficiency 
  h_effcor->Divide( h_unfolded, h_eff_cor_temp );


}

double getBinArea(TH2* h, int globalBin)
{
  int binx, biny, binz;
  h->GetBinXYZ(globalBin, binx, biny, binz);
  return h->GetXaxis()->GetBinWidth(binx)*h->GetYaxis()->GetBinWidth(biny);
}

#ifndef NEWCOV
TMatrixD divideCovByHist(TMatrixD& m, TH2D* h)
{
  TMatrixD ret(m);
  for(int i=0; i<h->fN; ++i){
    for(int j=0; j<h->fN; ++j){
      ret(i,j)=m(i,j)/(h->fArray[i]*h->fArray[j]);
    }
  }
  return ret;
}

TMatrixD divideCovByHists(TMatrixD m, TH2D* num, TH2D* dem)
{ 

  TH2D *tmp = new TH2D(*num);
  tmp->Divide(dem);
  TMatrixD ret(m);
  cout << "Dividing cov by hists" << num->fN << endl;
  for(int i=0; i<num->fN; ++i){
    for(int j=0; j<num->fN; ++j){
      int x,y,z,usex1,usey1,usex2,usey2;
      num->GetBinXYZ(i,x,y,z);
      usex1=x;
      usey1=y;
      num->GetBinXYZ(j,x,y,z);
      usex2=x;
      usey2=y;
      double eff1 = tmp->GetBinContent(usex1,usey1);
      double eff2 = tmp->GetBinContent(usex2,usey2);
      double val = eff1*eff2;
      if(val!=0)ret(i,j)=m(i,j)/(eff1*eff2);
    }
  }
  return ret;
}

TMatrixD divideCovByBinWidth(TMatrixD& m, TH2D* h)
{
  TMatrixD ret(m);
  for(int i=0; i<h->fN; ++i){
    for(int j=0; j<h->fN; ++j){
      ret(i,j)=m(i,j)/(getBinArea(h, i)*getBinArea(h, j));
    }
  }
  return ret;
}
#endif

//need to clean out result for unreportable bins
void SetBinZero(MnvH2D *&xs, int x, int y){

  xs->SetBinContent(x,y,0);
  vector<string> verterrnames = xs->GetVertErrorBandNames();
  vector<string> laterrnames = xs->GetLatErrorBandNames();
  //Vertical Errors
  for(int i=0;i<verterrnames.size();i++){
    MnvVertErrorBand2D *tmperr = xs->GetVertErrorBand(verterrnames[i]);
    int nhist = tmperr->GetNHists();
    for(int j=0;j<nhist;j++){
      tmperr->GetHist(j)->SetBinContent(x,y,0);
    }
  }
  //Lateral Errors
  for(int i=0;i<laterrnames.size();i++){
    MnvLatErrorBand2D *tmperr = xs->GetLatErrorBand(laterrnames[i]);
    int nhist = tmperr->GetNHists();
    for(int j=0;j<nhist;j++){
      tmperr->GetHist(j)->SetBinContent(x,y,0);
    }
  }
}

void SetBinValue(MnvH2D *&xs, int x, int y, double v){

  xs->SetBinContent(x,y,v);
  vector<string> verterrnames = xs->GetVertErrorBandNames();
  vector<string> laterrnames = xs->GetLatErrorBandNames();
  //Vertical Errors
  for(int i=0;i<verterrnames.size();i++){
    MnvVertErrorBand2D *tmperr = xs->GetVertErrorBand(verterrnames[i]);
    int nhist = tmperr->GetNHists();
    for(int j=0;j<nhist;j++){
      tmperr->GetHist(j)->SetBinContent(x,y,v);
    }
  }
  //Lateral Errors
  for(int i=0;i<laterrnames.size();i++){
    MnvLatErrorBand2D *tmperr = xs->GetLatErrorBand(laterrnames[i]);
    int nhist = tmperr->GetNHists();
    for(int j=0;j<nhist;j++){
      tmperr->GetHist(j)->SetBinContent(x,y,v);
    }
  }
}

void ZeroUnreported(MnvH2D *&xs){
  int nBinsX = xs->GetNbinsX()+2;//Include over/under flows
  int nBinsY = xs->GetNbinsY()+2;//Include over/under flows
  for(int i=0;i<nBinsX;i++){
    for(int j=0;j<nBinsY;j++){
  
      if(i==0 || j==0 || i==nBinsX-1 || j==nBinsY-1) SetBinZero(xs,i,j);
      //if(xs->GetBinContent(i,j) == 0.) cout<<"Zero!"<<endl;//SetBinValue(xs,i,j,2.);
      //if(!(std::isnan(xs->GetBinContent(i,j)))) SetBinValue(xs,i,j,3.);
      //double binvalue = xs->GetBinContent(i,j);
      //if(binvalue == binvalue) SetBinValue(xs,i,j,3.);
      //if(!(std::isnan(xs->GetBinContent(i,j)))) cout<<"NaN!"<<endl;//SetBinValue(xs,i,j,3.);

      //if(i==1 && j>=10) SetBinZero(xs,i,j);
      //else if(i==1 && j>=10) SetBinZero(xs,i,j);
      //else if(i==2 && j>=11) SetBinZero(xs,i,j);
      //else if(i==3 && j>=12) SetBinZero(xs,i,j);
      //else if(i==4 && j>=12) SetBinZero(xs,i,j);
      //else if(i==5 && j>=13) SetBinZero(xs,i,j);
      //else if(i==6 && j>=14) SetBinZero(xs,i,j);
      //if(i==1 && j>=9) SetBinZero(xs,i,j);
      //else if(i==2 && j>=10) SetBinZero(xs,i,j);
      //else if(i==3 && j>=11) SetBinZero(xs,i,j);
      //else if(i==4 && j>=12) SetBinZero(xs,i,j);
      //else if(i==5 && j>=13) SetBinZero(xs,i,j);
      //else if(i==6 && j>=13) SetBinZero(xs,i,j);
      //else if(j==1) SetBinZero(xs,i,j);
    }
  }
}

void ZeroUnreported(TMatrixD &matd,MnvH2D*xs){
  //Zero the stupid stat covariance matrix we have to cart around
  //mapping of cov matrix to x,y bin
  //Big bins = y axis bins
  //small bins = x axis bins
  int nBinsX = xs->GetNbinsX()+2;//Include over/under flows
  int nBinsY = xs->GetNbinsY()+2;//Include over/under flows
  //20 pt and 30 pz. So bin 35 is pt bin 1 and the 35-30 or bin 5 of pz
  // 35/30 = 1 (with 5 left over) 35-1*30 = 5. 0 to 29 is the first and 30-49 is the second so you subtract 1 off the little bin
  int covMatrixN = nBinsX*nBinsY;
  for(int i=0;i<covMatrixN;i++){
    int bigBin_i = i/nBinsX;
    int littleBin_i = i%nBinsX;
    for(int j=0;j<covMatrixN;j++){
      int bigBin_j = j/nBinsX;
      int littleBin_j = j%nBinsX;
      
      if(littleBin_i==0 || bigBin_i==0 || littleBin_i==nBinsX-1 || bigBin_i==nBinsY-1)matd[i][j]=0;


      //if(littleBin_i==1 && bigBin_i>=10) matd[i][j]=0;
      //else if(littleBin_i==1 && bigBin_i>=10) matd[i][j]=0;
      //else if(littleBin_i==2 && bigBin_i>=11) matd[i][j]=0;
      //else if(littleBin_i==3 && bigBin_i>=12) matd[i][j]=0;
      //else if(littleBin_i==4 && bigBin_i>=12) matd[i][j]=0;
      //else if(littleBin_i==5 && bigBin_i>=13) matd[i][j]=0;
      //else if(littleBin_i==6 && bigBin_i>=14) matd[i][j]=0;
      
      if(littleBin_j==0 || bigBin_j==0 || littleBin_j==nBinsX-1 || bigBin_j==nBinsY-1)matd[i][j]=0;

      //if(littleBin_j==1 && bigBin_j>=10) matd[i][j]=0;
      //else if(littleBin_j==1 && bigBin_j>=10) matd[i][j]=0;
      //else if(littleBin_j==2 && bigBin_j>=11) matd[i][j]=0;
      //else if(littleBin_j==3 && bigBin_j>=12) matd[i][j]=0;
      //else if(littleBin_j==4 && bigBin_j>=12) matd[i][j]=0;
      //else if(littleBin_j==5 && bigBin_j>=13) matd[i][j]=0;
      //else if(littleBin_j==6 && bigBin_j>=14) matd[i][j]=0;

      //if(xs->GetBinContent(i,j) == 0.)  cout<<"Zero!"<<endl;//matd[i][j]=2;
      if(std::isnan(xs->GetBinContent(i,j)))  cout<<"NaN!"<<endl;//matd[i][j]=3;
      //if(matd[i][j] == 0)  cout<<"Zero! Matrix"<<endl;//matd[i][j]=2;
      if(std::isnan(matd[i][j]))  cout<<"NaN! Matrix"<<endl;//matd[i][j]=3;
      //if(matd[i][j]==matd[i][j])  matd[i][j]=3;
      //if(matd[i][j] == 0) matd[i][j]=2;
      //if(std::isnan(matd[i][j])) matd[i][j]=3;

    }
  }
}



void NormalizeByFluxAndTargets( MnvH2D *&h_normalized, MnvH2D *h_effcor, MnvH1D* h_flux, double pot_scale, bool applyFluxConstraint, NukeCCUtilsNSF* utils, TMatrixD unfoldCov, int targetID = 0, int targetZ = 0, MnvH2D *flux_norm = NULL)
{

  //-------------------------------------------------
  //! Normalize by Flux
  //--------------------------------------------------
  h_flux->AddMissingErrorBandsAndFillWithCV( *h_effcor->ProjectionX() );
 
  MnvH2D *h_flux_normalization = (MnvH2D*)h_effcor->Clone("h_flux_normalization");
  h_flux_normalization->ClearAllErrorBands();
  h_flux_normalization->Reset();

  //Calculate the flux integrated value for each universe and set those values in bins of pt-pz
  //Integrate from 0.0 to 100.0 GeV (JO:we're removing the neutrino energy cut)
  double e_min = 0.0;//GeV
  double e_max = 100.;//GeV
  int b_min = h_flux->FindBin( e_min );
  int b_max = h_flux->FindBin( e_max );

  double flux_cv = h_flux->Integral( b_min, b_max, "width" );

  const int lowBin = 0;
  const int highBinX = h_flux_normalization->GetNbinsX()+1;
  const int highBinY = h_flux_normalization->GetNbinsY()+1;
  const int highBin = h_flux_normalization->GetBin( highBinX, highBinY );
  if(flux_norm==NULL){
    cout <<"cv\t" << flux_cv << endl;
    for( int i=lowBin; i <= highBin; ++i ){
      h_flux_normalization->SetBinContent( i, flux_cv );
    }

    cout << "Now for verts" << endl;
    //Do the same with the vertical error bands
    std::vector<std::string> vertNames = h_flux->GetVertErrorBandNames();
    for(unsigned int k=0; k<vertNames.size(); ++k )
      {
        cout << k << endl;
	MnvVertErrorBand *errBand = h_flux->GetVertErrorBand( vertNames[k] );
	//const int universes = errBand->GetNHists();
	int universes = errBand->GetNHists();
	//const int nuniverses = errBand->GetNHists();
        //int universes = 0;
        //if(nuniverses > 100) universes = 100;
        //else universes = nuniverses;  
        if(vertNames[k]=="Flux") universes = 100;//cout << "Flux Universes = " << universes << endl;
        cout << "N UNIVERSES VERT BANDS = " << universes << endl;
	std::vector<TH2D*> vert_hists;
	for(int u=0; u< universes; ++u)
	  {
	    
	    TH2D *h_vert = new TH2D( (TH2D)*h_flux_normalization);
	    h_vert->SetName( Form("h_vert_%s_%i", vertNames[k].c_str(), u) );
	    if(flux_norm==NULL){
	      double flux_vert = errBand->GetHist( u )->Integral( b_min, b_max, "width");
	      //	cout <<"vert\t" << u<<"\t" << flux_vert << endl;
	      for( int i=0; i <= highBin; ++i ){
		h_vert->SetBinContent( i, flux_vert );
	      }
	    }
	    vert_hists.push_back( h_vert );
	  }
	h_flux_normalization->AddVertErrorBand( vertNames[k], vert_hists );
      
	
	//cleaning
	for(std::vector<TH2D*>::iterator itHist = vert_hists.begin(); itHist != vert_hists.end(); ++itHist)
	  delete *itHist;
      }
    cout << "Now for lateral" << endl;
    //Do the same with the lateral error bands
    std::vector<std::string> latNames = h_flux->GetLatErrorBandNames();
    for(unsigned int k=0; k<latNames.size(); ++k )
      {
	cout << k << "\t" << latNames[k] << endl;
	MnvLatErrorBand *errBand = h_flux->GetLatErrorBand( latNames[k] );
	const int universes = errBand->GetNHists();
	//const int nuniverses = errBand->GetNHists();
        //int universes = 0;
        //if(nuniverses > 100) universes = 100;
        //else universes = nuniverses;  
        cout << "N UNIVERSES LAT BANDS = " << universes << endl;
	std::vector<TH2D*> lat_hists;
	cout << k << "\t" << latNames[k] << "\t" <<universes << endl;
	for(int u=0; u< universes; ++u)
	  {
	    TH2D *h_lat = new TH2D( (TH2D)*h_flux_normalization);
	    h_lat->SetName( Form("h_lat_%s_%i", latNames[k].c_str(), u) );
	    if(flux_norm==NULL){
	      double flux_lat = errBand->GetHist( u )->Integral( b_min, b_max, "width");
	      for( int i=0; i <= highBin; ++i ){
		h_lat->SetBinContent( i, flux_lat );
	      }
	    }
	    lat_hists.push_back( h_lat );
	  }
	h_flux_normalization->AddLatErrorBand( latNames[k], lat_hists );
	
	//cleaning
	for(std::vector<TH2D*>::iterator itHist = lat_hists.begin(); itHist != lat_hists.end(); ++itHist)
	  delete *itHist;
      }
  }
  else{
    h_flux_normalization=(MnvH2D*)flux_norm->Clone("enuflux");
    //h_flux_normalization->PopVertErrorBand( "Flux" );
    std::cout << "TRACE: Add missing error bands to " << h_flux_normalization->GetName() << std::endl;
    h_flux_normalization->AddMissingErrorBandsAndFillWithCV( *h_effcor->ProjectionX() );
  }


  cout << "Now the scales" << endl;
  //Convert flux units from m^2/POT to cm^2/POT
  h_flux_normalization->Scale( 1.0e-4 );
  h_normalized = (MnvH2D*)h_effcor->Clone( Form("%s_cross_section", h_effcor->GetName() ) );
   
/*
  std::vector<std::string> ErrorBandNames = h_normalized->GetVertErrorBandNames();
  for( std::vector<std::string>::iterator itName = ErrorBandNames.begin(); itName != ErrorBandNames.end(); ++itName ) {
    cout << "h_normalized Vert error band name is " << *itName << endl;
  }
  std::vector<std::string> ErrorBandNames2 = h_flux_normalization->GetVertErrorBandNames();
  for( std::vector<std::string>::iterator itName = ErrorBandNames2.begin(); itName != ErrorBandNames2.end(); ++itName ) {
    cout << "h_flux_normalization Vert error band name is " << *itName << endl;
  }
*/
  std::cout << "TRACE: Divide h_normalized by h_flux_normalization " << h_normalized->GetName() << std::endl;
  h_normalized->Divide( h_normalized, h_flux_normalization );

  //flux correct cov matrix
#ifndef NEWCOV
  TMatrixD unfoldingMatrixToScale = divideCovByHist(unfoldCov,h_flux_normalization);
#endif
  //-------------------------------------------------
  //! Normalize by number of CarbonAtoms on targets
  //--------------------------------------------------
  const int nplanes = 2 * ( 80 - 27 + 1 );//fiducial volume -> modules 27-80
  TString name = h_normalized->GetName();

//////*********From NukeCC package about nucleon correction**********//////////
  double Nucleons = -999.;
         
     if(targetZ==6)
       Nucleons =  TargetUtils::Get().GetPassiveTargetNNucleons( 3, targetZ, 0 );
     if(targetZ==26)
       Nucleons =  TargetUtils::Get().GetPassiveTargetNNucleons( targetID, targetZ, 0 );
     if(targetZ==82)
       Nucleons =  TargetUtils::Get().GetPassiveTargetNNucleons( targetID, targetZ, 0 );
     if(targetID>10 )
       Nucleons =  TargetUtils::Get().GetTrackerNNucleons( 12.0, 0);



  cout << " Number of Nucleons: " << Nucleons << endl;

  double targets = 0.0;
  if( name.Contains( "qelike" ) ){ //My MC is background subtracted so take the qelike component
    //targets = TargetUtils::Get().GetTrackerNCarbonAtoms( nplanes, /*isMC =*/ true, 850 );
    targets = TargetUtils::Get().GetTrackerNNucleons( nplanes, /*isMC =*/ true, 850 );
    //targets = TargetUtils::Get().GetTrackerNNucleons( nplanes, /*isMC =*/ false, 850 );
  }
  else {
    //targets = TargetUtils::Get().GetTrackerNCarbonAtoms( nplanes, /*isMC =*/ false, 850 );
    targets = TargetUtils::Get().GetTrackerNNucleons( nplanes, /*isMC =*/ false, 850 );
  }

  cout << TargetUtils::Get().GetTrackerNProtons(5980.,8422.,true,850.) << "\t" << TargetUtils::Get().GetNPlanes(5980,8422) <<"'\t This code." << endl;
  cout << TargetUtils::Get().GetTrackerNProtons(5990.,8340.,true,850.) << "\t" << TargetUtils::Get().GetNPlanes(5990,8340)<< "\t extractor. " << endl;
  cout << TargetUtils::Get().GetTrackerNCarbonAtoms( nplanes, /*isMC =*/ true, 850 ) << endl;
  cout << TargetUtils::Get().GetTrackerNProtons( nplanes, /*isMC =*/ true, 850 ) << endl;
  cout <<"MC\t" <<TargetUtils::Get().GetTrackerNNucleons( nplanes, /*isMC =*/ true, 850 ) << endl;//MC
  cout <<"Data\t" <<TargetUtils::Get().GetTrackerNNucleons( nplanes, /*isMC =*/ false, 850 ) << endl;//Data
  cout <<"MC Protons\t" <<TargetUtils::Get().GetTrackerNProtons( nplanes, /*isMC =*/ true, 850 ) << endl;//MC
  cout <<"Data Protons\t" <<TargetUtils::Get().GetTrackerNProtons( nplanes, /*isMC =*/ false, 850 ) << endl;//Data
  //Print out some useful information...
  double flux = h_flux_normalization->GetBinContent(1);
  double nC12 = TargetUtils::Get().GetTrackerNCarbonAtoms( nplanes, true, 850 );
  cout << " number of targets: " << targets << " neutrons" << endl;
  cout << " number of C12(MC): " << nC12 << endl;
  cout << " ratio of nC12/neutrons: " << nC12/targets << endl;
  cout << " flux integrated: " << flux << endl;
  
  // Scale by targets and POT scale
  //double scale = 1.0 / ( targets * pot_scale );
  double scale = 1.0 / ( Nucleons * pot_scale );
  std::cout << "TRACE: Scale h_normalized " << h_normalized->GetName() << std::endl;
  h_normalized->Scale( scale );
  
  //Change units in Z axis
  h_normalized->GetZaxis()->SetTitle( "d#sigma_{DIS}/d{#mu}_{P_Z}d{#mu}_{P_T} (cm^{2}/GeV^{2}/C^{12})" );

  //scale,bin width norm, pushback
#ifndef NEWCOV
    cout << "Scale the matrix" << endl;
  unfoldingMatrixToScale*= scale*scale;
  TMatrixD finalCovMatrix = divideCovByBinWidth(unfoldingMatrixToScale,h_normalized);
  // Jeremy tells me that the covariance matrix has the diagonal
  // errors on it, which are already included elsewhere, so we have to
  // subtract them off before adding the unfolding covariance matrix
  // back on
  for(int i=0; i<finalCovMatrix.GetNrows(); ++i) finalCovMatrix(i,i)=0;


  //OKAY ZERO OUT BINS (ONLY IN pzmu_ptmu
  string mytitle = h_normalized->GetName();
  cout << "@@@@@@@@@@@@@@@@@@@@@@@@@@" << endl;
  cout << mytitle << endl;
  cout << mytitle.find("pzmu") << endl;
  //if(mytitle.find("pzmu")!=string::npos){
    //cout << "I'm going zero out bins" << endl;
  ZeroUnreported(h_normalized);
  ZeroUnreported(finalCovMatrix,h_normalized);
  //}
  h_normalized->PushCovMatrix("unfoldingCov",finalCovMatrix);
  
#endif
  
#ifdef NEWCOV
   string mytitle = h_normalized->GetName();
    cout << "TRACE: new way, get the unfolding cov matrix,  remove unreported and diagonal errors and put it back" << endl;
    // remove the overflow bins and put it back.
    if (h_normalized->HasErrorMatrix("unfoldingCov")){
      std::cout << "TRACE there is an unfoldingCov here" << std::endl;
      TMatrixD* CovMatrix = (h_normalized->PopSysErrorMatrix("unfoldingCov"));
      TMatrixD test = TMatrixD(*CovMatrix);
      test.ResizeTo(20,20);
      test.Print();
      ZeroUnreported(h_normalized);
      
      ZeroUnreported(*CovMatrix,h_normalized);
      
      h_normalized->PushSysErrorMatrix("unfoldingCov",CovMatrix);
    }
    else{
      cout << "ERROR? << expect to have an unfoldingCov at this point " << endl;
    }
  std::cout << "TRACE: put matrix into finalCovMatrix " << h_normalized->GetName() << std::endl;
    TMatrixD finalCovMatrix = h_normalized->GetSysErrorMatrix("unfoldingCov");
  TMatrixD test = finalCovMatrix;
  test.ResizeTo(20,20);
  test.Print();
    
#endif
  
  finalCovMatrix.SaveAs("CovMat_tocheck.root");

  if(mytitle.find("Emu_Ehad_data")!=string::npos) finalCovMatrix.SaveAs("CovMat_tocheck_pzmu.root");
  if(mytitle.find("x_y_data")!=string::npos) finalCovMatrix.SaveAs("CovMat_tocheck_enu.root");
  if(mytitle.find("q2_ptmu_data")!=string::npos) finalCovMatrix.SaveAs("CovMat_tocheck_q2.root");

  TCanvas c1("can","can");
  finalCovMatrix.Draw("COLZ");
  //if(mytitle.find("Emu_Ehad_data")!=string::npos) c1.Print("CovMat_tocheck_pzmu_noNanto3.png");
  c1.Print("CovMat_tocheck_pzmu_noNanto3.png");
  if(mytitle.find("x_y_data")!=string::npos) c1.Print("CovMat_tocheck_enu_noNanto3.png");
  if(mytitle.find("q2_ptmu_data")!=string::npos) c1.Print("CovMat_tocheck_q2_noNanto3.png");

  delete h_flux_normalization; 
}

void ZeroDiagonal(TMatrixD &m){
  std::cout << " TRACE: enter ZeroDiagonal  "   << std::endl;
  int n = m.GetNrows();
  for (int i = 0; i < n; i++){
    m[i][i] = 0;
  }
}


int CrossSectionHists(std::string event_filename,std::string scale_filename,std::string migration_filename,std::string efficiency_filename,std::string output_filename,int num_iter,int ana_type, int targetID, int targetZ ){
   TFile *scale_file = new TFile(scale_filename.c_str(),"READ");
   TFile *event_file = new TFile(event_filename.c_str(),"READ");
   TFile *migration_file = new TFile(migration_filename.c_str(),"READ");
   TFile *eff_file = new TFile(efficiency_filename.c_str(),"READ");
   
///////////*****Modern way of POT counting************/////////////////////
   TParameter<double> *mcPOT = (TParameter<double>*)event_file->Get("MCPOT");
   TParameter<double> *dataPOT = (TParameter<double>*)event_file->Get("DataPOT");
   double pot_mc_event = mcPOT->GetVal();
   double pot_data_event = dataPOT->GetVal();
//   double dataMCScale = datapot/mcpot;   
  
   TParameter<double> *pot_mig = (TParameter<double>*)migration_file->Get("MCPOT");
   TParameter<double> *pot_efficiency = (TParameter<double>*)eff_file->Get("MCPOT");
   double pot_migration = pot_mig->GetVal();
   double pot_eff = pot_efficiency->GetVal();

   std::cout<<"***Event File*** "<<event_file->GetName()<<std::endl; 
   std::cout<<"***Efficiency File*** "<<eff_file->GetName()<<std::endl; 
   std::cout<<"***Migration File*** "<<migration_file->GetName()<<std::endl; 
   
   std::cout<<"***POT: MC event "<<pot_mc_event<<std::endl; 
   std::cout<<"***POT: Data event "<<pot_data_event<<std::endl; 
   std::cout<<"***POT: migration "<<pot_migration<<std::endl;
   std::cout<<"***POT: efficiency "<<pot_eff<<std::endl;  
 
   //FIRST THING FIRST.....CHECK THAT THE POT IS CONSISTENT OVER ALL INPUT FILES.....
   std::cout<<"CHECKING POT INFO IN ALL INPUT ROOT FILES "<<std::endl;
   TVector2 *pot_sf = (TVector2*)scale_file->Get("pot");
   //TVector2 *pot_event = (TVector2*)event_file->Get("pot");
   //TVector2 *pot_migration = (TVector2*)migration_file->Get("pot");
   //TVector2 *pot_eff = (TVector2*)eff_file->Get("pot");
   /*
   std::cout<<"***POT: scale "<<pot_sf->X()<<" "<<pot_sf->Y()<<std::endl;
   std::cout<<"***POT: event "<<pot_event->X()<<" "<<pot_event->Y()<<std::endl; 
   std::cout<<"***POT: migration "<<pot_migration->X()<<" "<<pot_migration->Y()<<std::endl;
   std::cout<<"***POT: efficiency "<<pot_eff->X()<<" "<<pot_eff->Y()<<std::endl;  
   */

   
   bool applyFluxConstraint = false;
   
   
   

  //NukeCCUtilsNSF *utils = new NukeCCUtilsNSF(false,"minervame5A"); //you can give any name you want to give here by the way....
  NukeCCUtilsNSF *utils = new NukeCCUtilsNSF("minervame1A"); //you can give any name you want to give here by the way....
  MinervaUnfold::MnvUnfold unfold;
  

  
   
 //FluxReweighter *frw = new FluxReweighter(-14,false,FluxReweighter::minervame5A,FluxReweighter::gen2thin, FluxReweighter::g4numiv6, m_fluxUniverses); 
 //FluxReweighter *frw = new FluxReweighter(-14,false,FluxReweighter::minervame5A,FluxReweighter::gen2thin, FluxReweighter::g4numiv6, 100); 
 FluxReweighter *frw = new FluxReweighter(14,false,FluxReweighter::minervame1A,FluxReweighter::gen2thin, FluxReweighter::g4numiv6, 100); 

   
   //const std::string hist_type[2]  = {"h_mc_Emu_Ehad","h_mc_x_y"};
   //const std::string truth_type[2]  = {"h_truth_Emu_Ehad","h_truth_Emu_Ehad"};
   
   const std::string migration_type[2] = {"selected_mc_response2d_Emu_Emu_migration","selected_mc_response2d_x_y_migration"};
 
//#ifndef DEBUG

   std::cout<<" I Am HERE 1"<<std::endl;
//#endif   
   
   
 //  std::map<std::string,MnvH2D*>bkgsub_hist_type = GetFitBinBkgSubtraction(event_file,scale_file,hist_type[ana_type]);   

    MnvH2D *mc_component = (MnvH2D*)event_file->Get("h_mc_Emu_Ehad");  
    MnvH2D *data_component = (MnvH2D*)event_file->Get("h_data_Emu_Ehad");  
  

  std::cout << " data component " << data_component->GetName() << std::endl;

/*
        for( int bin1=1; bin1 < data_component->GetNbinsX()+1; bin1++ ){
        for( int bin2=1; bin2 < data_component->GetNbinsY()+1; bin2++ ){
        cout << "data bin content, bin error? type? " << (bin1,bin2) << ", " << data_component->GetBinContent(bin1, bin2) << ", " << data_component->GetBinError(bin1,bin2) << endl;
       } }
        for( int bin1=1; bin1 < mc_component->GetNbinsX()+1; bin1++ ){
        for( int bin2=1; bin2 < mc_component->GetNbinsY()+1; bin2++ ){
        cout << "mc bin content, bin error? type? " << (bin1, bin2) << ", " << mc_component->GetBinContent(bin1, bin2) << ", " << mc_component->GetBinError(bin1, bin2) << endl;
       }}
   */
 
    //MnvH2D *mc_component = GetTotalComponent(event_file,hist_type[ana_type],"mc");
    //MnvH2D *data_component = GetTotalComponent(event_file,hist_type[ana_type],"data");
    //MnvH2D *bkg_component = GetTotalComponent(event_file,hist_type[ana_type],"qelikenot"); 
    //std::map<std::string,MnvH2D*>signal_dict = GetTotalSignalComponent(event_file,hist_type[ana_type]);
    MnvH2D *signal_dict = (MnvH2D*)event_file->Get("h_mc_Emu_Ehad");
//    std::map<std::string,MnvH2D*>bkg_dict = GetTotalBackgroundComponent(event_file,hist_type[ana_type]);

//#ifndef DEBUG
std::cout<<" I Am HERE 2"<<std::endl;
//#endif    
    
  migration_file->ls();
  //MnvH2D* h_hist_type_reco = (MnvH2D*)migration_file->Get((migration_type[ana_type]+"qelike_reco").c_str());
  MnvH2D* h_hist_type_reco = (MnvH2D*)migration_file->Get("h_mc_reco_Emu_Ehad");
  MnvH2D* h_hist_type_generated = (MnvH2D*)migration_file->Get("h_mc_true_Emu_Ehad");
  MnvH2D *h_hist_type_migration = (MnvH2D*)migration_file->Get("selected_mc_response2d_Emu_Ehad_migration");
  
  //I dont want this CV universe in my list of systematics...need to figure out how to get rid of it in MigrationMatrix....
  std::cout << " hist_type_reco " << h_hist_type_reco->GetName() << std::endl;
  
  std::vector<std::string>error_names = h_hist_type_reco->GetErrorBandNames();
  for(int i=0; i <error_names.size();i++){
    if(error_names[i]=="cv"){
      h_hist_type_reco->PopVertErrorBand("cv");
      h_hist_type_generated->PopVertErrorBand("cv");
      h_hist_type_migration->PopVertErrorBand("cv");
      
    }
    
  }
  
    
  TH2D *migstatDummy = new TH2D(h_hist_type_migration->GetCVHistoWithStatError());
  migstatDummy->SaveAs("h_Emu_Ehad_migration.root");    
    
     
   //now the Efficiency correction.....
   
   //MnvH1D *h_flux = frw->GetFluxReweighted(-14);
   MnvH1D *h_flux = frw->GetFluxReweighted(14);
   
   h_flux->SaveAs("h_flux.root");
   
   
    double pot_scale = pot_data_event/pot_mc_event;
   //double pot_migration_scale = pot_migration->X()/pot_migration->Y();
   //double pot_eff_scale = pot_eff->X()/pot_eff->Y();
   
 
   if(pot_scale==0)pot_scale = 1.0;
  cout<< "__________"<<endl;
  cout<< "POT Information: "<<endl;
  cout<< "POT Data = " << pot_data_event << endl;
  cout<< "POT MC   = " << pot_mc_event << endl;
  cout<< "POT Scale factor = " << pot_scale << endl;
  cout<< "---------------------------------------"  << endl;
  
  
   //signal_dict["qelike"]->Scale(pot_scale);
  // h_hist_type_reco->Scale(pot_migration_scale);
  // h_hist_type_generated->Scale(pot_migration_scale);
  // h_hist_type_migration->Scale(pot_migration_scale);
   

   
   //std::map<std::string,MnvH2D*>bkgsub_hist = GetFitBinBkgSubtraction(event_file,scale_file,hist_type[ana_type]);
 

   //MnvH2D *h_mc_no_bck = h_hist_type_reco->Clone((hist_type[ana_type]+"qelikebkgsub").c_str());
   MnvH2D *h_mc_no_bck = signal_dict->Clone("h_mc_Emu_Ehad");
   h_mc_no_bck->AddMissingErrorBandsAndFillWithCV(*mc_component);
  // bkgsub_hist["scaled_bkg"]->AddMissingErrorBandsAndFillWithCV(*mc_component);
  // h_mc_no_bck->Add(bkgsub_hist["scaled_bkg"],-1.0);
//   MnvH2D *h_data_no_bck = bkgsub_hist["bkgsub_data"];

   std::cout<<"***No BKG hist name*** "<<h_mc_no_bck->GetName()<<std::endl; 
   //MnvH2D *h_mc_tuned_bck = bkgsub_hist["scaled_bkg"];

   //now the Unfolding...........
  
   //MnvH2D *h_mc_unfolded = new MnvH2D(*h_mc_no_bck);
   MnvH2D *h_mc_unfolded = (MnvH2D*)event_file->Get("h_mc_Emu_Ehad");
   h_mc_unfolded->SetName("h_mc_Emu_Ehad_unfolded");
   std::cout<<"***MC Unfold  name*** "<<h_mc_unfolded->GetName()<<std::endl; 

   MnvH2D *h_data_unfolded = new MnvH2D(*h_mc_unfolded);
   //MnvH2D *h_data_unfolded = (MnvH2D*)event_file->Get("h_data_Emu_Ehad");
   h_data_unfolded->SetName("h_data_Emu_Ehad_unfolded");
   h_data_unfolded->AddMissingErrorBandsAndFillWithCV(*h_hist_type_migration);

   MnvH2D *h_data_no_bck = (MnvH2D*)event_file->Get("h_data_Emu_Ehad");
   h_data_no_bck->AddMissingErrorBandsAndFillWithCV(*h_hist_type_migration);
   TMatrixD unfoldingCovMatrixOrig_hist_type;
   std::cout<<"***Data Unfold  name*** "<<h_data_unfolded->GetName()<<std::endl; 
  /* 
        for( int bin1=1; bin1 < h_data_unfolded->GetNbinsX()+1; bin1++ ){
        for( int bin2=1; bin2 < h_data_unfolded->GetNbinsY()+1; bin2++ ){
        cout << "data bin content, bin error? type? " << (bin1,bin2) << ", " << h_data_unfolded->GetBinContent(bin1, bin2) << ", " << h_data_unfolded->GetBinError(bin1,bin2) << endl;
       } }
        for( int bin1=1; bin1 < h_mc_unfolded->GetNbinsX()+1; bin1++ ){
        for( int bin2=1; bin2 < h_mc_unfolded->GetNbinsY()+1; bin2++ ){
        cout << "mc Before bin content, bin error? type? " << (bin1,bin2) << ", " << h_mc_unfolded->GetBinContent(bin1, bin2) << ", " << h_mc_unfolded->GetBinError(bin1,bin2) << endl;
       } }
*/
   
   bool data_unfolded = unfold.UnfoldHisto2D(h_data_unfolded,h_hist_type_migration,h_hist_type_reco,h_hist_type_generated,h_data_no_bck,num_iter,true,true);
   bool mc_unfolded = unfold.UnfoldHisto2D(h_mc_unfolded,h_hist_type_migration,h_hist_type_reco,h_hist_type_generated,h_mc_no_bck,num_iter,true,true);
     
  /*      for( int bin1=1; bin1 < h_mc_unfolded->GetNbinsX()+1; bin1++ ){
        for( int bin2=1; bin2 < h_mc_unfolded->GetNbinsY()+1; bin2++ ){
        cout << "mc after unfolding bin content, bin error? type? " << (bin1,bin2) << ", " << h_mc_unfolded->GetBinContent(bin1, bin2) << ", " << h_mc_unfolded->GetBinError(bin1,bin2) << endl;
       } }
*/
   
   std::cout<<" Status of the unfolding is DATA/MC"<<data_unfolded<<"/"<<mc_unfolded<<std::endl;
   
   
   //now we need to get the covariance matrix out of the unfolded distribution.....
   
   
  TH2D* hUnfoldedDummy=new TH2D(h_data_unfolded->GetCVHistoWithStatError());
  TH2D* hMigrationDummy=new TH2D(h_hist_type_migration->GetCVHistoWithStatError());
  TH2D* hRecoDummy=new TH2D(h_hist_type_reco->GetCVHistoWithStatError());
  TH2D* hTruthDummy=new TH2D(h_hist_type_generated->GetCVHistoWithStatError());
  TH2D* hBGSubDataDummy=new TH2D(h_data_no_bck->GetCVHistoWithStatError());



 unfold.UnfoldHisto2D(hUnfoldedDummy, unfoldingCovMatrixOrig_hist_type, hMigrationDummy, hRecoDummy, hTruthDummy, hBGSubDataDummy, num_iter);
 // There's a bug in RooUnfold that's making it return covariance
  std::cout << "TRACE : make the cov matrix " << std::endl;
  TMatrixD test = TMatrixD(unfoldingCovMatrixOrig_hist_type);
  
  test.ResizeTo(20,20);
  test.Print();
   
   int correctNbins = hUnfoldedDummy->fN;
   int matrixRows = unfoldingCovMatrixOrig_hist_type.GetNrows();
   
   if(correctNbins!=matrixRows){
   
 	cout << "****************************************************************************" << endl;
	cout << "*  Fixing unfolding matrix size because of RooUnfold bug. From " << matrixRows << " to " << correctNbins << endl;
	cout << "****************************************************************************" << endl;
	// It looks like this DTRT, since the extra last two bins don't have any content
	unfoldingCovMatrixOrig_hist_type.ResizeTo(correctNbins, correctNbins);
	
	

        //unfoldingCovMatrixOrig_pzmu_ptmu.Print();
    
   }
  ZeroDiagonal(unfoldingCovMatrixOrig_hist_type);
  #ifdef NEWCOV
          cout << "TRACE: NEWCOV put the unfolding matrix into the MNvH2D"  << endl;
          //ZeroDiagonal(unfoldingCovMatrixOrig_pzmu_ptmu); // don;'t double count stat error
          h_data_unfolded->PushCovMatrix("unfoldingCov",unfoldingCovMatrixOrig_hist_type);
  #endif
   
   //NOW THE EFFICIENCY CORRECTION......

   //I want to do this for all the components in fact.....
   
   //create maps for efficiency numerator/denominator efficiency and efficiency corrected histograms...
   
   
   MnvH2D *h_hist_type_data_effcor = NULL;
   MnvH2D *h_hist_type_mc_effcor = NULL;
   
   std::string components[6] = {"qelike","qelike_2p2h","qelike_res","qelike_2p2h","qelike_dis","qelike_qe"};
      
 
  
  cout<< "---------------------------------------"  << endl;
  cout<< "Correcting By Efficiency"  << endl;
  cout<< "---------------------------------------"  << endl;
  //Let's calculate the efficiency on the fly now
  
  MnvH2D *h_hist_type_effhist_num = (MnvH2D*)eff_file->Get("h_mc_Emu_Ehad");
  MnvH2D *h_hist_type_effhist_den = (MnvH2D*)eff_file->Get("h_truth_Emu_Ehad");
   
  std::cout<<"***Efficiency corrected Num*** "<<h_hist_type_effhist_num->GetName()<<std::endl; 
  std::cout<<"***Efficiency corrected Denom*** "<<h_hist_type_effhist_den->GetName()<<std::endl; 
  
  h_hist_type_effhist_num->AddMissingErrorBandsAndFillWithCV(*h_data_unfolded);  // make certain you get all the bands
  h_hist_type_effhist_den->AddMissingErrorBandsAndFillWithCV(*h_data_unfolded);
  
  
  
  
  //do the pot normalization of the efficiency components as well...
 // h_hist_type_effhist_num->Scale(pot_eff_scale);
//  h_hist_type_effhist_den->Scale(pot_eff_scale);
  
  CorrectByEfficiency(h_hist_type_data_effcor,h_data_unfolded,h_hist_type_effhist_num,h_hist_type_effhist_den);
  CorrectByEfficiency(h_hist_type_mc_effcor,h_mc_unfolded,h_hist_type_effhist_num,h_hist_type_effhist_den);
  
/*        for( int binX=1; binX < h_hist_type_effhist_num->GetNbinsX()+1; binX++ ){
        for( int binY=1; binY < h_hist_type_effhist_num->GetNbinsY()+1; binY++ ){
        cout << "binX content, binX error? type? " << (binX, binY) << ", " << h_hist_type_effhist_num->GetBinContent(binX, binY) << ", " << h_hist_type_effhist_num->GetBinError(binX, binY) << endl;
       }}
        for( int binX=1; binX < h_hist_type_effhist_den->GetNbinsX()+1; binX++ ){
        for( int binY=1; binY < h_hist_type_effhist_den->GetNbinsY()+1; binY++ ){
        cout << "binX content, binX error? type? " << (binX, binY) << ", " << h_hist_type_effhist_den->GetBinContent(binX, binY) << ", " << h_hist_type_effhist_den->GetBinError(binX, binY) << endl;
       }}


        for( int bin1=1; bin1 < h_data_unfolded->GetNbinsX()+1; bin1++ ){
        for( int bin2=1; bin2 < h_data_unfolded->GetNbinsY()+1; bin2++ ){
        cout << "data bin content, bin error? type? " << (bin1,bin2) << ", " << h_data_unfolded->GetBinContent(bin1, bin2) << ", " << h_data_unfolded->GetBinError(bin1,bin2) << endl;
       } }
        for( int bin1=1; bin1 < h_mc_unfolded->GetNbinsX()+1; bin1++ ){
        for( int bin2=1; bin2 < h_mc_unfolded->GetNbinsY()+1; bin2++ ){
        cout << "mc Before bin content, bin error? type? " << (bin1,bin2) << ", " << h_mc_unfolded->GetBinContent(bin1, bin2) << ", " << h_mc_unfolded->GetBinError(bin1,bin2) << endl;
       } }
*/
   
   std::string effcorname = ana_type==0?"h_effcor_Emu_Ehad_":"h_effcor_x_y_";
      
   
   std::cout<<"Efficiency correction by the unfolding matrix "<<std::endl;
#ifndef NEWCOV   
   TMatrixD unfoldingCovMatrixEff_hist_type;
   std::cout<<" Here the covariance matrix is divided by efficiency "<<std::endl;
   TMatrixD tmpMat = divideCovByHists(unfoldingCovMatrixOrig_hist_type,h_hist_type_effhist_num,h_hist_type_effhist_den);
   int ncols = tmpMat.GetNcols();
   unfoldingCovMatrixEff_hist_type.ResizeTo(ncols,ncols);
   unfoldingCovMatrixEff_hist_type = tmpMat;
   
#else
   TMatrixD unfoldingCovMatrixEff_hist_type = h_hist_type_data_effcor->GetSysErrorMatrix("unfoldingCov");
  std::cout << "TRACE: check Eff " <<  unfoldingCovMatrixEff_hist_type[40][41] << std::endl;
   h_hist_type_data_effcor->MnvH2DToCSV("effcorNew","./",1.0,false,true,false);
   std::cout<<"TRACE::Got UnfoldingCov "<<unfoldingCovMatrixEff_hist_type.GetNrows()<<std::endl;
   test = unfoldingCovMatrixEff_hist_type;
  test.ResizeTo(20,20);
  test.Print();
#endif
   
   //Now the Efficiency correction of the various compoments of the MC (true compoments) to get the signal component cross-sections....
      
   std::cout<<"--------NORMALIZATION BY FLUX AND TARGET---------"<<std::endl;
   //Now the Normalization by flux and number of targets
   MnvH2D *h_flux_normalization=NULL; //needs special enu binned flux if enu analysis
   
   if(ana_type==1)h_flux_normalization = GetEnuFluxHists(h_hist_type_effhist_den,h_flux);
      
    std::map<std::string,MnvH2D*>comp_crosslist;  
    /* 
    for(int i =0;i<6;i++){
    //Get the Truth distribution (the efficiency denominator)......
      MnvH2D *h_hist_type_truthcomp = (MnvH2D*)eff_file->Get((truth_type[ana_type]).c_str());
      h_hist_type_truthcomp->AddMissingErrorBandsAndFillWithCV(*h_data_unfolded);
      //normalize by POT
     // h_hist_type_truthcomp->Scale(pot_eff_scale);
      
      MnvH2D *h_hist_type_cross_section_comp = new MnvH2D();
      //Normalize these components by Flux and Number of targets.....
      
      
    std::cout << "doing specials Basically events divided by flux and pot" << std::endl;
    std::cout<<" SIGNAL COMPONENT: "<< components[i]<<std::endl;
    if(ana_type==0)NormalizeByFluxAndTargets( h_hist_type_cross_section_comp, h_hist_type_truthcomp, h_flux, pot_data_event, applyFluxConstraint, utils, unfoldingCovMatrixEff_hist_type);
    
     else{
     NormalizeByFluxAndTargets( h_hist_type_cross_section_comp, h_hist_type_truthcomp, h_flux, pot_data_event, applyFluxConstraint, utils, unfoldingCovMatrixEff_hist_type, h_flux_normalization);
    
      }
      comp_crosslist[hist_type[ana_type]+components[i]+"_crossSection"] = h_hist_type_cross_section_comp;
      

    }
    */
  std::cout << "TRACE: make and fill the cross section histograms " << std::endl;
    //now the cross-secton for unfolded distributions...needs to be efficiency corrected by the unfoldig matrix as well....
    MnvH2D *h_hist_type_cross_section_mc = (MnvH2D*) h_hist_type_mc_effcor->Clone("temp_xs_mc"); // clone to make certain bands get in.
    MnvH2D *h_hist_type_cross_section_data = (MnvH2D*) h_hist_type_data_effcor->Clone("temp_xs_data");
     
    
    if(ana_type==0){  
    NormalizeByFluxAndTargets( h_hist_type_cross_section_data, h_hist_type_data_effcor, h_flux, pot_data_event, applyFluxConstraint, utils, unfoldingCovMatrixEff_hist_type, targetID, targetZ);
    NormalizeByFluxAndTargets( h_hist_type_cross_section_mc, h_hist_type_mc_effcor, h_flux, pot_data_event, applyFluxConstraint, utils, unfoldingCovMatrixEff_hist_type, targetID, targetZ);
    
    }
    
    if(ana_type==1){
    
    NormalizeByFluxAndTargets( h_hist_type_cross_section_data, h_hist_type_data_effcor, h_flux, pot_data_event, applyFluxConstraint, utils, unfoldingCovMatrixEff_hist_type, targetID, targetZ, h_flux_normalization);
    NormalizeByFluxAndTargets( h_hist_type_cross_section_mc, h_hist_type_mc_effcor, h_flux, pot_data_event, applyFluxConstraint, utils, unfoldingCovMatrixEff_hist_type, targetID, targetZ, h_flux_normalization);
    
    }
  
    cout<<"*******After Normalization outputs"<<endl;
       /* 
        for( int binX=1; binX < h_hist_type_mc_effcor->GetNbinsX()+1; binX++ ){
        for( int binY=1; binY < h_hist_type_mc_effcor->GetNbinsY()+1; binY++ ){
        cout << "MC Eff corrected binX content, binX error? type? " << (binX, binY) << ", " << h_hist_type_mc_effcor->GetBinContent(binX, binY) << ", " << h_hist_type_mc_effcor->GetBinError(binX, binY) << endl;
       }}
        for( int binX=1; binX < h_hist_type_data_effcor->GetNbinsX()+1; binX++ ){
        for( int binY=1; binY < h_hist_type_data_effcor->GetNbinsY()+1; binY++ ){
        cout << "Data Eff corrected binX content, binX error? type? " << (binX, binY) << ", " << h_hist_type_data_effcor->GetBinContent(binX, binY) << ", " << h_hist_type_data_effcor->GetBinError(binX, binY) << endl;
       }}

        for( int binX=1; binX < h_hist_type_cross_section_mc->GetNbinsX()+1; binX++ ){
        for( int binY=1; binY < h_hist_type_cross_section_mc->GetNbinsY()+1; binY++ ){
        cout << "MC CrossSection binX content, binX error? type? " << (binX, binY) << ", " << h_hist_type_cross_section_mc->GetBinContent(binX, binY) << ", " << h_hist_type_cross_section_mc->GetBinError(binX, binY) << endl;
       }}
        for( int binX=1; binX < h_hist_type_cross_section_data->GetNbinsX()+1; binX++ ){
        for( int binY=1; binY < h_hist_type_cross_section_data->GetNbinsY()+1; binY++ ){
        cout << "Data CrossSection binX content, binX error? type? " << (binX, binY) << ", " << h_hist_type_cross_section_data->GetBinContent(binX, binY) << ", " << h_hist_type_cross_section_data->GetBinError(binX, binY) << endl;
       }}

    */
   
   
   
   
   //now check and fix the flux error bands
   
     // This is for the nue constraint only
//   frw->CheckAndFixFluxErrorBand(h_hist_type_cross_section_mc);
//   frw->CheckAndFixFluxErrorBand(h_hist_type_cross_section_data);

   TFile *outfile = TFile::Open(output_filename.c_str(),"RECREATE");
   outfile->cd();
   for(auto hists:comp_crosslist)(hists.second)->Write((hists.first).c_str());
   //for(auto hists:signal_dict)(hists.second)->Write();
//   for(auto hists:bkg_dict)(hists.second)->Write();
   mc_component->Write("h_mc_Emu_Ehad");
   data_component->Write("h_data_Emu_Ehad");
//   bkg_component->Write((hist_type[ana_type]+"qelikenot").c_str());
   h_data_unfolded->Write("h_data_Emu_Ehad_unfolded");
   h_mc_unfolded->Write("h_mc_Emu_Ehad_unfolded");
   h_hist_type_migration->Write();
   h_hist_type_generated->Write();
   //h_data_no_bck->Write();
   h_mc_no_bck->Write();
   h_hist_type_effhist_num->Write();
   h_hist_type_effhist_den->Write();
   h_hist_type_migration->Write();
   h_hist_type_data_effcor->Write((effcorname+"data").c_str());
   h_hist_type_mc_effcor->Write((effcorname+"mc").c_str());
   h_hist_type_cross_section_data->Write("h_dataCrossSection");
   h_hist_type_cross_section_mc->Write("h_recoCrossSection");
   //pot_mc->Write("MCPOT");
   //pot_data->Write("DataPOT");
#ifndef NEWCOV
   unfoldingCovMatrixOrig_hist_type.Write("Unfolding_CovMarix");
#endif
   
   
   outfile->Close();
   for(auto hists:comp_crosslist)delete hists.second;
   comp_crosslist.clear();
   delete utils;
   delete mc_component;
   delete data_component;
//   delete bkg_component;
   delete h_data_unfolded;
   delete h_mc_unfolded;
   delete h_hist_type_migration;
   delete h_hist_type_generated;
   //delete h_data_no_bck;
   delete h_mc_no_bck;
   delete h_hist_type_effhist_num;
   delete h_hist_type_effhist_den;
   delete h_hist_type_cross_section_data;
   delete h_hist_type_cross_section_mc;
   delete h_hist_type_data_effcor;
   delete h_hist_type_mc_effcor;
   
   
   delete scale_file;
   delete migration_file;
   delete eff_file;
   delete event_file;
   
   delete pot_sf;
   //delete pot_mc;
   //delete pot_data;
   //delete pot_migration;
   //delete pot_eff;
   
   


  return 1;
}





int main( int argc, char *argv[])
{
  #ifndef ROOT6
  //ROOT::Cintex::Cintex::Enable();
  #endif
  TH1::AddDirectory(false);

  if (argc==1){
    std::cout<<"-----------------------------------------------------------------------------------------------"<<std::endl;
    std::cout<<"MACROS HELP:\n\n"<<
      "\t-./CrossSectionHists  Name_and_path_to_EventRate_file Name_and_Path_to_BkgScaleFile Name_and_path_to_MigrationMatrixFile Name_and_Path_to_EfficiencyFile  Name_and_Path_to_Output_Filename num_iter ana_type"<<std::endl;
    std::cout<<"-----------------------------------------------------------------------------------------------"<<std::endl;
    return 0; 
  }
  
  //! Default parameters
  std::vector<std::string> par;
  par.push_back("CrossSectionHists");
  par.push_back( Form("%s/ana/rootfiles/EventRates.root",getenv("CCQENuNSFROOT") ) );
  par.push_back(Form("%s/ana/rootfiles/scale_file.root",getenv("CCQENUNSFROOT")));
  par.push_back( Form("%s/ana/rootfiles/Migration.root",getenv("CCQENuNSFROOT")));
  par.push_back(Form("%s/ana/rootfiles/Efficiency.root",getenv("CCQENuNSFROOT")));
  par.push_back(Form("%s/ana/rootfiles/output.root",getenv("CCQENuNSFROOT")));
  par.push_back(Form("4"));
  par.push_back(Form("1"));
  par.push_back(Form("1"));
  par.push_back(Form("26"));
  


  //! Set user parameters
  for( int i=0; i<argc; ++i){
    par.at(i) = argv[i];
  }


  for( unsigned int i=0; i<par.size(); ++i)
    std::cout<<"Parameter "<< i << ": " << par[i] << std::endl;

  return CrossSectionHists(par[1], par[2],par[3],par[4],par[5],atoi(par[6].c_str()),atoi(par[7].c_str()), atoi(par[8].c_str()), atoi(par[9].c_str()));
  
}
