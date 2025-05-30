#ifndef CVUNIVERSE_cxx
#define CVUNIVERSE_cxx 1

#include "../include/CVUniverse.h"
#include "PlotUtils/FluxReweighter.h"
#include "PlotUtils/HyperDimLinearizer.h"
#include "PlotUtils/MinosMuonEfficiencyCorrection.h"
//#include "PlotUtils/MinosMuonPlusEfficiencyCorrection.h"
#include "PlotUtils/PlotUtilsPhysicalConstants.h"
#include <cmath>
#include "PlotUtils/MinervaUniverse.h"
#include "PlotUtils/ChainWrapper.h"
//using namespace globalV;
using namespace NUKECC_ANA;

std::string pwd = "$MPARAMFILESROOT/data/Calibrations/energy_calib/CalorimetryTunings.txt";

util::CaloCorrection Nu_Tgt1_Fe(pwd.c_str(), "NukeCC_Nu_Tgt1_Fe");
util::CaloCorrection Nu_Tgt1_Pb(pwd.c_str(), "NukeCC_Nu_Tgt1_Pb");
util::CaloCorrection Nu_Tgt2_Fe(pwd.c_str(), "NukeCC_Nu_Tgt2_Fe");
util::CaloCorrection Nu_Tgt2_Pb(pwd.c_str(), "NukeCC_Nu_Tgt2_Pb");
util::CaloCorrection Nu_Tgt3_C(pwd.c_str(),  "NukeCC_Nu_Tgt3_C");
util::CaloCorrection Nu_Tgt3_Fe(pwd.c_str(), "NukeCC_Nu_Tgt3_Fe");
util::CaloCorrection Nu_Tgt3_Pb(pwd.c_str(), "NukeCC_Nu_Tgt3_Pb");
util::CaloCorrection Nu_Tgt4_Pb(pwd.c_str(), "NukeCC_Nu_Tgt4_Pb");
util::CaloCorrection Nu_Tgt5_Fe(pwd.c_str(), "NukeCC_Nu_Tgt5_Fe");
util::CaloCorrection Nu_Tgt5_Pb(pwd.c_str(), "NukeCC_Nu_Tgt5_Pb");
util::CaloCorrection Nu_Tracker(pwd.c_str(), "NukeCC_AntiNu_Tracker");
//util::CaloCorrection Nu_Tracker(pwd.c_str(), "CCInclusiveTrackerTunings");


//again the constructor.....
CVUniverse::CVUniverse(PlotUtils::ChainWrapper *chw,double nsigma):PlotUtils::MinervaUniverse(chw,nsigma)
{
  //just go with 2 universe world
 // recoShifter = new MnvRecoShifter(neutrinoMode,2);


}

  //destructor....
CVUniverse::~CVUniverse(){
 // delete recoShifter;
  }
 Long64_t CVUniverse::GetMaxEntries()
{

      TTree          *fChain;   //!pointer to the analyzed TTree or TChain
    //! If MAX_ENTRIES environmental is set to something positive, use it.  Otherwise use them all.
    if( getenv("MAX_ENTRIES") )
    {
        int maxEntries = atoi(getenv( "MAX_ENTRIES") );
        if( maxEntries > 0  )
             return maxEntries;
    }
      return fChain->GetEntries();
}

//latest way of getting recoil energy at macro level instead from tuples
double CVUniverse::GetCalRecoilEnergy() const {
    //return 1.0;
    return GetDouble("part_response_total_recoil_passive_allNonMuonClusters_id") + GetDouble("part_response_total_recoil_passive_allNonMuonClusters_od");
}
double CVUniverse::GetNonCalRecoilEnergy() const {return 0.0;} // added for New Particle Response Systematic



//TBV version
double CVUniverse::ApplyCaloTuning(double calRecoilE) const{

    int targetID = GetInt((GetAnaToolName() +"_targetID").c_str());
    int targetZ = GetInt((GetAnaToolName() +"_targetZ").c_str());
      if (targetID == 1){
         if(targetZ == 26) return Nu_Tgt1_Fe.eCorrection(calRecoilE*mev_to_gev)/mev_to_gev; //MeV
         else if(targetZ == 82) return Nu_Tgt1_Pb.eCorrection(calRecoilE*mev_to_gev)/mev_to_gev; //MeV
         else{
             std::cout << "No spline found for target 1." << std::endl;
             return calRecoilE; //MeV
         }
      }
      else if (targetID == 2){
        if(targetZ == 26) return Nu_Tgt2_Fe.eCorrection(calRecoilE*mev_to_gev)/mev_to_gev; //MeV
        else if(targetZ == 82) return Nu_Tgt2_Pb.eCorrection(calRecoilE*mev_to_gev)/mev_to_gev; //MeV
        else{
            std::cout << "No spline found for target 2." << std::endl;
            return calRecoilE; //MeV
        }
     }
     else if (targetID == 3){
        if(targetZ == 6) return Nu_Tgt3_C.eCorrection(calRecoilE*mev_to_gev)/mev_to_gev; //MeV
        else if(targetZ == 26) return Nu_Tgt3_Fe.eCorrection(calRecoilE*mev_to_gev)/mev_to_gev; //MeV
        else if(targetZ == 82) return Nu_Tgt3_Pb.eCorrection(calRecoilE*mev_to_gev)/mev_to_gev; //MeV
        else{
            std::cout << "No spline found for target 3." << std::endl;
            return calRecoilE; //MeV
        }
     }
     else if(targetID == 4){
        if(targetZ == 82) return Nu_Tgt4_Pb.eCorrection(calRecoilE*mev_to_gev)/mev_to_gev; //MeV
        else{
            std::cout << "No spline found for target 4." << std::endl;
            return calRecoilE; //MeV
        }
     }
     else if(targetID == 5){
        if(targetZ == 26) return Nu_Tgt5_Fe.eCorrection(calRecoilE*mev_to_gev)/mev_to_gev; //MeV
        else if(targetZ == 82) return Nu_Tgt5_Pb.eCorrection(calRecoilE*mev_to_gev)/mev_to_gev; //MeV
        else{
            std::cout << "No spline found for target 5." << std::endl;
            return calRecoilE; //MeV
        }
     }
     else{
        //std::cout<<"Applying caloTuning for calRecoilE: "<<calRecoilE<<std::endl;
        return Nu_Tracker.eCorrection(calRecoilE*mev_to_gev)/mev_to_gev; //MeV
     }

     //return 1.0;
}



//calling branches and variables
//if (useDNN) {
//double CVUniverse::GetRecoilEnergy()  const { return GetDouble((GetAnaToolName() + "_recoil_E").c_str()); }
//} else {
//double CVUniverse::GetRecoilEnergy()  const { return GetDouble((GetAnaToolName() + "_ANN_recoil_E").c_str());}
//}

//double CVUniverse::GetMuonCurve()  const { return 1/GetDouble("(GetAnaToolName() +_minos_trk_eqp_qp"); }
double CVUniverse::GetMuonCurve()  const { return 1/GetDouble((GetAnaToolName() + "_minos_trk_eqp_qp").c_str()); }
double CVUniverse::GetHelicity()  const { return GetInt((GetAnaToolName() +"_nuHelicity").c_str()); }
double CVUniverse::GetTrueHelicity()  const { return GetInt("mc_incomig"); }
double CVUniverse::GetFiducial()  const { return GetInt((GetAnaToolName() +"_in_fiducial_area").c_str()); }
double CVUniverse::GetTdead()  const { return GetInt("phys_n_dead_discr_pair_upstream_prim_track_proj"); }

double CVUniverse::GetTargetID()   const{return GetInt((GetAnaToolName() + "_targetID").c_str());}
double CVUniverse::GetANNTargetID()   const{return GetInt("ANN_targetID");}

std::tuple<std::vector<int>,std::vector<int>,std::vector<double>,std::vector<double> > CVUniverse::GetFateNum()   const{return calFateNumVector();}

//double CVUniverse::GetThetamu()       const {return GetDouble("muon_theta");}
double CVUniverse::GetEnu()           const {return GetEmu()+ GetRecoilEnergy();}
double CVUniverse::GetQ2Reco()        const { return   calcQ2(GetEnu(),  GetEmu(), GetThetamu() ); }//in MeV^2
double CVUniverse::GetWReco()         const{return (calcW(GetQ2Reco(), GetRecoilEnergy()));}//in MeV
double CVUniverse::GetxReco()         const{return calcX(GetQ2Reco(), GetEnu(),GetEmu());}
double CVUniverse::GetyReco()         const{return calcY(GetEnu(), GetRecoilEnergy());}

// for true variables

double CVUniverse::GetTruthNuPDG() const{return GetInt("mc_incoming");}
double CVUniverse::GetCurrent()    const{return GetInt("mc_current");}
double CVUniverse::GetEhadTrue()    const{return GetEnuTrue() - GetElepTrue();}
double CVUniverse::GetThetamuTrue( )const{return GetDouble("truth_muon_theta");}

double CVUniverse::GetQ2IncTrue()      const{return (calcQ2(GetEnuTrue(),  GetElepTrue(), GetThetamuTrue())); }
double CVUniverse::GetWTrue()       const{return (calcW(GetQ2IncTrue(), GetEhadTrue()));}

double CVUniverse::GetxTrue()       const{return calcX(GetQ2IncTrue(), GetEnuTrue(),GetElepTrue());}
double CVUniverse::GetyTrue()       const{ return calcY(GetEnuTrue(), GetEhadTrue());}

//double CVUniverse::Getq3True()      const{return calcq3(GetQ2True(),  GetEnuTrue(), GetElepTrue());}
//double CVUniverse::Getq0True()      const{return calcq0(GetEnuTrue(), GetElepTrue());}
double CVUniverse::Getq0Truth()     const{return calcq0(GetEnuTrue(), GetElepTrue());}
double CVUniverse::GetThetamuDeg()  const{return GetThetamu()*rad_to_deg;}
double CVUniverse::GetThetamuTrueDeg()  const{return GetThetamuTrue()*rad_to_deg;}
double CVUniverse::GetMuonP()       const{return GetPmu();}

// in GeV

double CVUniverse::GetMuonEGeV()           const {return GetEmu()*mev_to_gev;}
double CVUniverse::GetMuonETrueGeV()           const {return GetElepTrue()*mev_to_gev;}
double CVUniverse::GetEnuGeV()           const {return (GetEmu()+ GetRecoilEnergy())*mev_to_gev;}

double CVUniverse::GetEnuTrueGeV()           const {return (GetEnuTrue())*mev_to_gev;}
double CVUniverse::GetEhadGeV()  const { return GetRecoilEnergy()*mev_to_gev; }

double CVUniverse::GetEhadTrueGeV()    const{return (GetEnuTrue() - GetElepTrue())*mev_to_gev;}


double CVUniverse::GetQ2RecoGeV()        const { return (GetQ2Reco() * mev_to_gev *mev_to_gev); }
double CVUniverse::GetWRecoGeV()         const{return (GetWReco() *mev_to_gev);}

double CVUniverse::GetQ2TrueGeV()      const{return (GetQ2IncTrue()*mev_to_gev*mev_to_gev); }
double CVUniverse::GetWTrueGeV()       const{return (GetWTrue()*mev_to_gev);}

double CVUniverse::GetMuonPZ()      const{return GetPmu()/1000. * cos(GetThetamu()); }
double CVUniverse::GetMuonPZTrue()  const{return GetPlepTrue()/1000. * cos(GetThetalepTrue()); } // Gev/c

double CVUniverse::GetEavReco()	   const{ return GetEavailReco(); }
//double CVUniverse::GetEavRecoECal() const{return GetEavailRecoECal();}
double CVUniverse::GetCVWeight()   const{return calcCVWeight();}
double CVUniverse::Getq3Reco()     const{ return calcRecoq3();}
double CVUniverse::GetECALHCALAvEnergy()     const{ return calcECALHCALAvEn();}
double CVUniverse::GetTrackerECALAvEnergy()     const{ return  calcTrackerECALAvEn();}
double CVUniverse::GetEavailTrue()    const{ return calcEavailTrue();}
double CVUniverse::Getq3Truth()     const{ return calcTrueq3();}
//double CVUniverse::GetKEFSITruth()     const{ return calcKEFSI();}

//MECAna way to calculate the variables
double CVUniverse::MECAnaGetTrackerECALAvEnergy()     const{ return  MECAnacalcTrackerECALAvEn();}
double CVUniverse::MECAnaGetq3Reco()     const{ return MECAnacalcRecoq3();}

//Low recoil
double CVUniverse::GetEavalReco_Low()	   const{ return GetEavailRecoECal(); }

double CVUniverse::GetMuonPt()const{
   return GetMuonP()*sin(GetThetamu());
}
//=======================================================================
double CVUniverse::GetMuonPz()const{
  return GetMuon4V().Z();
}
//=======================================================================
double CVUniverse::GetlepPtTrue()const{
   return GetPlepTrue()*sin(GetThetalepTrue());
}
//=======================================================================
double CVUniverse::GetVertexZNew() const {
    double vertex = GetVecElem("ANN_vtx", 2)/10; // in this current ntuple, Z position in mm not cm, hence divide
    //double vertex = GetVecElem((GetAnaToolName() + "_vtx").c_str(), 2)/10; // in this current ntuple, Z position in mm not cm, hence divide
    return vertex;
}
//=======================================================================
double CVUniverse::GetVertexZTrueNew() const {
    return GetVecElem("mc_vtx", 2)/10;
}
//=======================================================================
double CVUniverse::GetVertexZMy() const {
    double vertex = GetVecElem((GetAnaToolName() + "_vtx").c_str(), 2)/10; // in this current ntuple, Z position in mm not cm, hence divide
    return vertex;
}
//=======================================================================
double CVUniverse::GetVertexZTrueMy() const {
    return GetVecElem("mc_vtx", 2)/10;
}
//=======================================================================
//For Background///////
double CVUniverse::GetplaneDNNReco()  const{
           int ANN_vtx_module, ANN_vtx_plane;
           int Segment = GetVecElem("ANN_segments",0);
           int targetID = GetTargetFromSegment(Segment, ANN_vtx_module, ANN_vtx_plane);
           double ANN_VTX_module = (double)ANN_vtx_module;
           double ANN_VTX_plane = (double)ANN_vtx_plane;
           double test=ANN_VTX_module*2.+ANN_VTX_plane+10.;
//           return  ANN_VTX_module*2.+ANN_VTX_plane+10.;}
           return test;
}
//=======================================================================
double CVUniverse::GetplaneDNNTrue() const{
  return static_cast<double>(GetInt("truth_vtx_module")*2+GetInt("truth_vtx_plane")+10);
}
//======================================================================
std::pair<std::vector<double>,std::vector<double> > CVUniverse::getFateWeight(double FateNum) const {
   double weight = -99.0;
   std::vector<double> preFSIEnergy;
   std::vector<double> fateWeightVal;

   int mc_er_nPart = GetInt("mc_er_nPart");
   int targetA = GetInt("mc_targetA");

   for(int i=0; i < mc_er_nPart; ++i){

        int mc_er_status = GetVecElem("mc_er_status",i);
        double piEnergy = GetVecElem("mc_er_E",i);
        if( mc_er_status == 14){
            // Only track fates for pizero, piplus, piminus.
            int mc_er_ID = GetVecElem("mc_er_ID",i);
            //cout << "mc_er_ID: " << mc_er_ID << endl;
            //if(FateNum == 2){cout << "Entered to mc_er_status"<< " mc_er_ID: "<<mc_er_ID<<endl;}
            if(mc_er_ID == 111 || mc_er_ID==211 || mc_er_ID ==-211){
                int mot= GetVecElem("mc_er_mother",i);
                //Here is where the weight is going to be applied
                if (targetA == 207){
                      //if(FateNum ==2){
                      //cout << " Energy of the pion is: "<<piEnergy<<" the particle is: "<<mc_er_ID<< " fateWeight: "<<getFateTable(piEnergy, FateNum)<<" FateNum: "<<FateNum<<endl;}
			                double fateW = getFateTable(piEnergy, FateNum);
                      //cout << " piEnergy: "<< piEnergy << "fateW: "<< fateW << " fate: "<<FateNum <<endl;
                      fateWeightVal.push_back(fateW);
                      preFSIEnergy.push_back(piEnergy);
		            }
 		            else {
                  fateWeightVal.push_back(1.0);
                  preFSIEnergy.push_back(piEnergy);}

            }//end track fates for pions
            else {
              fateWeightVal.push_back(1.0);
              preFSIEnergy.push_back(piEnergy); }
         }//end mc_er_status == 14
         //else {return 1.0;}
   }//end mc_er_nPart

   return std::make_pair(fateWeightVal, preFSIEnergy);
}
//=======================================================================
std::pair<double,double> CVUniverse::getFateWeightSingle(double FateNum) const{
   double weight = -99.0;

   int mc_er_nPart = GetInt("mc_er_nPart");
   int targetA = GetInt("mc_targetA");
   //cout << "mc_er_nPart" << endl;
   for(int i=0; i < mc_er_nPart; ++i){

        int mc_er_status = GetVecElem("mc_er_status",i);
        double piEnergy = GetVecElem("mc_er_E",i);
        //cout << "mc_er_status weight: " << mc_er_status << endl;
        if( mc_er_status == 14){

            // Only track fates for pizero, piplus, piminus.
            int mc_er_ID = GetVecElem("mc_er_ID",i);
            //cout << "mc_er_ID: " << mc_er_ID << endl;
            //if(FateNum == 2){cout << "Entered to mc_er_status"<< " mc_er_ID: "<<mc_er_ID<<endl;}
            if(mc_er_ID == 111 || mc_er_ID==211 || mc_er_ID ==-211){
                int mot= GetVecElem("mc_er_mother",i);

                //Here is where the weight is going to be applied
                if (targetA == 207){
                      //if(FateNum ==2){
                      //cout << " Energy of the pion is: "<<piEnergy<<" the particle is: "<<mc_er_ID<< " fateWeight: "<<getFateTable(piEnergy, FateNum)<<" FateNum: "<<FateNum<<endl;}
			                double fateW = getFateTable(piEnergy, FateNum);
                      //cout << " piEnergy: "<< piEnergy << "fateW: "<< fateW << " fate: "<<FateNum <<endl;
                      return std::make_pair(fateW,piEnergy);
		            }
 		            else { return std::make_pair(1.0,piEnergy);}

            }//end track fates for pions
            else { return std::make_pair(1.0,piEnergy); }
         }//end mc_er_status == 14
         //else {return 1.0;}
   }//end mc_er_nPart
   //return 1.0;
}
//=============This the the old version to calculate the weight============================================
double CVUniverse::getFateTable(double piEnergy, int FateNum) const{

    //Inelastic
    if(FateNum == 1){
        if(piEnergy > 0 && piEnergy <=50 ){ return 0.163; }
        if(piEnergy > 50 && piEnergy <=100 ){ return 0.204; }
        if(piEnergy > 100 && piEnergy <=150 ){ return 0.500; }
        if(piEnergy > 150 && piEnergy <=200 ){ return 0.490; }
        if(piEnergy > 200 && piEnergy <=250 ){ return 1.117; }
        if(piEnergy > 250 && piEnergy <=300 ){ return 0.364; }
        if(piEnergy > 300 && piEnergy <=350 ){ return 0.614; }
        if(piEnergy > 350 && piEnergy <=400 ){ return 0.342; }
        if(piEnergy > 400 && piEnergy <=450 ){ return 0.346; }
        if(piEnergy > 450 && piEnergy <=500 ){ return 0.394; }
        if(piEnergy > 500 && piEnergy <=550 ){ return 0.347; }
        if(piEnergy > 550 && piEnergy <=600 ){ return 0.329; }
        if(piEnergy > 600 && piEnergy <=650 ){ return 0.300; }
        if(piEnergy > 650 && piEnergy <=700 ){ return 0.266; }
        if(piEnergy > 700 && piEnergy <=750 ){ return 0.245; }
        if(piEnergy > 750 && piEnergy <=800 ){ return 0.238; }
        if(piEnergy > 800 && piEnergy <=850 ){ return 0.236; }
        if(piEnergy > 850 && piEnergy <=900 ){ return 0.232; }
        if(piEnergy > 900 && piEnergy <=950 ){ return 0.223; }
        if(piEnergy > 950 && piEnergy <=1000 ){ return 0.212; }
        if(piEnergy > 1000){ return 0.212; }
    }
    else if(FateNum == 2){ //CEX
        if(piEnergy > 0 && piEnergy <=50 ){ return 2.8; }
        if(piEnergy > 50 && piEnergy <=100 ){ return 1.359; }
        if(piEnergy > 100 && piEnergy <=150 ){ return 2.425; }
        if(piEnergy > 150 && piEnergy <=200 ){ return 2.484; }
        if(piEnergy > 200 && piEnergy <=250 ){ return 1.902; }
        if(piEnergy > 250 && piEnergy <=300 ){ return 4.104; }
        if(piEnergy > 300 && piEnergy <=350 ){ return 3.054; }
        if(piEnergy > 350 && piEnergy <=400 ){ return 3.328; }
        if(piEnergy > 400 && piEnergy <=450 ){ return 4.201; }
        if(piEnergy > 450 && piEnergy <=500 ){ return 0.0; }
        if(piEnergy > 500 && piEnergy <=550 ){ return 0.0; }
        if(piEnergy > 550 && piEnergy <=600 ){ return 0.0; }
        if(piEnergy > 600 && piEnergy <=650 ){ return 0.0; }
        if(piEnergy > 650 && piEnergy <=700 ){ return 0.0; }
        if(piEnergy > 700 && piEnergy <=750 ){ return 0.0; }
        if(piEnergy > 750 && piEnergy <=800 ){ return 0.0; }
        if(piEnergy > 800 && piEnergy <=850 ){ return 0.0; }
        if(piEnergy > 850 && piEnergy <=900 ){ return 0.0; }
        if(piEnergy > 900 && piEnergy <=950 ){ return 0.0; }
        if(piEnergy > 950 && piEnergy <=1000 ){ return 0.0; }
        if(piEnergy > 1000){ return 0.0; }
    }
   else if(FateNum == 8){ //Pion Prod
        if(piEnergy > 0 && piEnergy <=50 ){ return 0.0; }
        if(piEnergy > 50 && piEnergy <=100 ){ return 0.0; }
        if(piEnergy > 100 && piEnergy <=150 ){ return 0.0; }
        if(piEnergy > 150 && piEnergy <=200 ){ return 0.0; }
        if(piEnergy > 200 && piEnergy <=250 ){ return 0.0; }
        if(piEnergy > 250 && piEnergy <=300 ){ return 0.0; }
        if(piEnergy > 300 && piEnergy <=350 ){ return 0.0; }
        if(piEnergy > 350 && piEnergy <=400 ){ return 1.0; }
        if(piEnergy > 400 && piEnergy <=450 ){ return 4.351; }
        if(piEnergy > 450 && piEnergy <=500 ){ return 6.492; }
        if(piEnergy > 500 && piEnergy <=550 ){ return 5.418; }
        if(piEnergy > 550 && piEnergy <=600 ){ return 5.275; }
        if(piEnergy > 600 && piEnergy <=650 ){ return 4.974; }
        if(piEnergy > 650 && piEnergy <=700 ){ return 4.515; }
        if(piEnergy > 700 && piEnergy <=750 ){ return 4.150; }
        if(piEnergy > 750 && piEnergy <=800 ){ return 3.941; }
        if(piEnergy > 800 && piEnergy <=850 ){ return 3.867; }
        if(piEnergy > 850 && piEnergy <=900 ){ return 3.881; }
        if(piEnergy > 900 && piEnergy <=950 ){ return 3.846; }
        if(piEnergy > 950 && piEnergy <=1000 ){ return 3.762; }
        if(piEnergy > 1000){ return 3.762; }
    }
   else if(FateNum == 5 || FateNum == 50 || FateNum == 51 || FateNum == 52){ //Absorbtion
        if(piEnergy > 0 && piEnergy <=50 ){ return 2.113; }
        if(piEnergy > 50 && piEnergy <=100 ){ return 1.985; }
        if(piEnergy > 100 && piEnergy <=150 ){ return 1.475; }
        if(piEnergy > 150 && piEnergy <=200 ){ return 1.183; }
        if(piEnergy > 200 && piEnergy <=250 ){ return 0.779; }
        if(piEnergy > 250 && piEnergy <=300 ){ return 0.969; }
        if(piEnergy > 300 && piEnergy <=350 ){ return 0.936; }
        if(piEnergy > 350 && piEnergy <=400 ){ return 0.998; }
        if(piEnergy > 400 && piEnergy <=450 ){ return 0.969; }
        if(piEnergy > 450 && piEnergy <=500 ){ return 0.737; }
        if(piEnergy > 500 && piEnergy <=550 ){ return 0.865; }
        if(piEnergy > 550 && piEnergy <=600 ){ return 0.791; }
        if(piEnergy > 600 && piEnergy <=650 ){ return 0.713; }
        if(piEnergy > 650 && piEnergy <=700 ){ return 0.659; }
        if(piEnergy > 700 && piEnergy <=750 ){ return 0.620; }
        if(piEnergy > 750 && piEnergy <=800 ){ return 0.594; }
        if(piEnergy > 800 && piEnergy <=850 ){ return 0.572; }
        if(piEnergy > 850 && piEnergy <=900 ){ return 0.549; }
        if(piEnergy > 900 && piEnergy <=950 ){ return 0.529; }
        if(piEnergy > 950 && piEnergy <=1000 ){ return 0.512; }
        if(piEnergy > 1000){ return 0.512; }
    }
   else {return 1.0;}

}


//=======Fate Weight using functional fits=================================
double CVUniverse::getResonantPionWeightJan2024( double piEnergy, int FateNum) const {
//No Scattering or Elastic
 if(FateNum == 1 || FateNum == 3){
   if(piEnergy < 0) { return 1.0; }
   if(piEnergy >= 0  && piEnergy <= 1000) {
     double weight = (-5.21168359e-14)*pow(piEnergy,5)+(1.57328334e-10)*pow(piEnergy,4)+(-1.77472513e-07)*pow(piEnergy,3)+(9.02098170e-05)*pow(piEnergy,2)+(-1.87718127e-02)*(piEnergy)+(1.60234611);
     return weight;
   }
   if(piEnergy > 1000) { return 0.8345;}
 }
 else if(FateNum == 2){ //CEX
    if(piEnergy < 0) { return 1.0; }
    if(piEnergy >= 0 && piEnergy <= 440){
      double weight = (8.58451494e-14)*pow(piEnergy,6) + (-1.20718361e-10)*pow(piEnergy,5)+(6.28942994e-8)*pow(piEnergy,4)+ (-1.49996818e-05)*pow(piEnergy,3) + (1.58236369e-03)*pow(piEnergy,2) + (-3.92965591e-02)*(piEnergy) + (8.81257014e-01);
      return weight;
    }
    if(piEnergy >= 440){ return 1.0;} // choose 1.0 or 0!
 }
 else if(FateNum == 10){ //Inelastic
   if(piEnergy < 0) { return 1.0; }
   if(piEnergy >= 0 && piEnergy <= 1000) {
     double weight = (3.04418208e-14)*pow(piEnergy,5) + (-9.79384698e-11)*pow(piEnergy,4) + (1.18010433e-07)*pow(piEnergy,3) + (-6.38560139e-05)*pow(piEnergy,2) + (1.36510822e-02)*(piEnergy) + (-6.79747222e-02);
     return weight; }
   if (piEnergy > 1000){ return 0.2497;}
 }
 else if(FateNum == 5 || FateNum == 50 || FateNum == 51 || FateNum == 52 || FateNum == 9){ //Absorption
   if(piEnergy < 0) { return 1.0; }
   if(piEnergy >= 0 && piEnergy <= 1000){
     double weight = (7.70276638e-19)*pow(piEnergy,7) + (-2.89092444e-15)*pow(piEnergy,6) + (4.32625973e-12)*pow(piEnergy,5) + (-3.27910968e-09)*pow(piEnergy,4) + (1.31889237e-06)*pow(piEnergy,3) + (-2.6699991e-04)*pow(piEnergy,2) + (2.07448722e-02)*(piEnergy) + (1.42728628);
     return weight;
   }
   if(piEnergy > 1000) { return 0.5103;}
 }
 else if(FateNum == 8){//Pion Production
   if(piEnergy < 0) { return 1.0; }
   if(piEnergy > 0 && piEnergy <= 360){return 0.0;}
   if(piEnergy > 360 && piEnergy <= 1000){
     double weight = (7.9964078e-15)*pow(piEnergy,6) + (-3.2627557e-11)*pow(piEnergy,5) + (5.39194844e-08)*pow(piEnergy,4) + (-4.600247e-05)*pow(piEnergy,3) + (2.127698e-02)*pow(piEnergy,2) + (-5.0428784)*(piEnergy) + (4.8417941e+02);
     return weight;
   }
   if(piEnergy > 1000) { return 4.0499;}
 }
 else{return 1.0;}
}

double CVUniverse::getQEWeightJan2024( double piEnergy, int FateNum) const {
  //No Scattering and Elastic
  if(FateNum == 1 || FateNum == 3){
    if(piEnergy < 0) { return 1.0; }
    if(piEnergy > 0 && piEnergy <= 1000) {
	 double weight = (1.52e-14)*pow(piEnergy,5)+(-4.43e-11)*pow(piEnergy,4)+(4.79e-08)*pow(piEnergy,3)+(-2.31e-05)*pow(piEnergy,2)+(4.7e-03)*(piEnergy)+(4.13e-01);
     return weight;
   }
   if(piEnergy > 1000) { return 0.8216;}
 }
 else if(FateNum == 2 || FateNum == 10){ //CEX and Inelastic (KnockOut)
    if(piEnergy < 0) { return 1.0; }
    if(piEnergy > 0 && piEnergy <= 1000){
      double weight = (-4.01e-14)*pow(piEnergy,5)+(1.13e-10)*pow(piEnergy,4)+(-1.17e-07)*pow(piEnergy,3) + (5.33e-05)*pow(piEnergy,2) + (-1e-02)*(piEnergy) + (1.8);
      return weight;
    }
    if(piEnergy > 1000){ return 1.1292;}
 }
 else if(FateNum == 5 || FateNum == 50 || FateNum == 51 || FateNum == 52 || FateNum == 9){ //Absorption
   if(piEnergy < 0) { return 1.0; }
   if(piEnergy > 0 && piEnergy <= 1000){
     double weight = (-2.4e-16)*pow(piEnergy,5) + (1.51e-12)*pow(piEnergy,4) + (-2.02e-9)*pow(piEnergy,3) +  (4.9e-7)*pow(piEnergy,2) +  (3.58e-5)*(piEnergy) +(1.21);
     return weight;
   }
   if(piEnergy > 1000) { return 0.9502;}
 }
 else if(FateNum == 8){ //Pion Production
   if(piEnergy < 0) { return 1.0; }
   if(piEnergy > 0 && piEnergy <= 360){return 0.0;}
   if(piEnergy > 360 && piEnergy <= 1000){
     double weight = (-2.66e-13)*pow(piEnergy,5) + (8.52e-10)*pow(piEnergy,4) + (-1.04e-06)*pow(piEnergy,3) +  (5.99e-04)*pow(piEnergy,2) +  (-1.59e-01)*(piEnergy) + (1.64e01);
     return weight;
   }
   if(piEnergy > 1000) { return 1.1157;}
 }
 else{return 1.0;}
}

double CVUniverse::getCarbonResonantPionWeightJan2024( double piEnergy, int FateNum) const {
//No Scattering or Elastic
 if(FateNum == 1 || FateNum == 3){
   if(piEnergy < 0) { return 1.0; }
   if(piEnergy >= 0  && piEnergy <= 1000) {
     double weight = (-1.53944923e-21)*pow(piEnergy,8)+(6.55392717e-18)*pow(piEnergy,7)+(-1.14569751e-14)*pow(piEnergy,6)+ (1.05403844e-11)*pow(piEnergy,5)+(-5.41061752e-09)*pow(piEnergy,4)+(1.49586063e-06)*pow(piEnergy,3)+(-1.88846670e-04)*pow(piEnergy,2)+(5.61591486e-03)*(piEnergy)+(1.05900673);
     return weight;
   }
   if(piEnergy > 1000) { return 1.0169;}
 }
 else if(FateNum == 2){ //CEX
    if(piEnergy < 0) { return 1.0; }
    if(piEnergy >= 0 && piEnergy <= 440){
      double weight = (1.06314699e-13)*pow(piEnergy,6) + (-1.55971182e-10)*pow(piEnergy,5)+(8.50310740e-08)*pow(piEnergy,4)+ (-2.14587214e-05)*pow(piEnergy,3) + (2.50531351e-03)*pow(piEnergy,2) + (-1.05860176e-01)*(piEnergy) + (3.38068686);
      return weight;
    }
    if(piEnergy >= 440){ return 1.0;}
 }
 else if(FateNum == 10){ //Inelastic
   if(piEnergy < 0) { return 1.0; }
   if(piEnergy >= 0 && piEnergy <= 1000) {
     double weight = (1.10590464e-16)*pow(piEnergy,6) + (-2.86347992e-13)*pow(piEnergy,5) + (2.37045715e-10)*pow(piEnergy,4) + (-4.22114879e-08)*pow(piEnergy,3) + (-2.86684181e-05)*pow(piEnergy,2) + (9.56608344e-03)*(piEnergy) + (3.54050813e-01);
     return weight; }
   if (piEnergy > 1000){ return 0.2210;}
 }
 else if(FateNum == 5 || FateNum == 50 || FateNum == 51 || FateNum == 52 || FateNum == 9){ //Absorption
   if(piEnergy < 0) { return 1.0; }
   if(piEnergy >= 0 && piEnergy <= 1000){
     double weight = (1.45521601e-14)*pow(piEnergy,5) + (-4.35072911e-11)*pow(piEnergy,4) + (4.78128955e-08)*pow(piEnergy,3) + (-2.21406008e-05)*pow(piEnergy,2) + (2.59438410e-03)*(piEnergy) + (1.03895039);
     return weight;
   }
   if(piEnergy > 1000) { return 0.3442;}
 }
 else if(FateNum == 8){//Pion Production
   if(piEnergy < 0) { return 1.0; }
   if(piEnergy > 0 && piEnergy <= 360){return 0.0;}
   if(piEnergy > 360 && piEnergy <= 1000){
     double weight = (9.43677192e-13)*pow(piEnergy,5) + (-3.41232242e-09)*pow(piEnergy,4) + (4.88204991e-06)*pow(piEnergy,3) + (-3.43614130e-03)*pow(piEnergy,2) + (1.17678529)*(piEnergy) + (-1.49577542e+02);
     return weight;
   }
   if(piEnergy > 1000) { return 4.1299;}
 }
 else{return 1.0;}
}

double CVUniverse::getCarbonQEWeightJan2024( double piEnergy, int FateNum) const {
  //No Scattering and Elastic
  if(FateNum == 1 || FateNum == 3){
    if(piEnergy < 0) { return 1.0; }
    if(piEnergy > 0 && piEnergy <= 1000) {
	 double weight = (-2.12994606e-17)*pow(piEnergy,6)+ (7.29793661e-14)*pow(piEnergy,5)+(-9.84990949e-11)*pow(piEnergy,4)+(6.59301972e-08)*pow(piEnergy,3)+(-2.21221306e-05)*pow(piEnergy,2)+(3.16651535e-03)*(piEnergy)+(9.19404645e-01);
     return weight;
   }
   if(piEnergy > 1000) { return 1.0803;}
 }
 else if(FateNum == 2 || FateNum == 10){ //CEX and Inelastic (KnockOut)
    if(piEnergy < 0) { return 1.0; }
    if(piEnergy > 0 && piEnergy <= 1000){
      double weight = (9.27682029e-17)*pow(piEnergy,6)+ (-3.05673815e-13)*pow(piEnergy,5)+(3.91370282e-10)*pow(piEnergy,4)+(-2.45118247e-07)*pow(piEnergy,3) + ( 7.68746741e-05)*pow(piEnergy,2) + (-1.06603108e-02)*(piEnergy) + (1.38836281);
      return weight;
    }
    if(piEnergy > 1000){ return 0.9161;}
 }
 else if(FateNum == 5 || FateNum == 50 || FateNum == 51 || FateNum == 52 || FateNum == 9){ //Absorption
   if(piEnergy < 0) { return 1.0; }
   if(piEnergy > 0 && piEnergy <= 1000){
     double weight = (-1.93637157e-19)*pow(piEnergy,7) + (6.36592798e-16)*pow(piEnergy,6) + (-8.09055354e-13)*pow(piEnergy,5) + (4.96355058e-10)*pow(piEnergy,4) + (-1.48321018e-07)*pow(piEnergy,3) +  (1.80881639e-05)*pow(piEnergy,2) +  (3.35405156e-06)*(piEnergy) +(8.25891852e-01);
     return weight;
   }
   if(piEnergy > 1000) { return 0.8312;}
 }
 else if(FateNum == 8){ //Pion Production
   if(piEnergy < 0) { return 1.0; }
   if(piEnergy > 0 && piEnergy <= 360){return 0.0;}
   if(piEnergy > 360 && piEnergy <= 1000){
     double weight = (-7.02486249e-14)*pow(piEnergy,5) + (2.11430221e-10)*pow(piEnergy,4) + (-2.29180148e-07)*pow(piEnergy,3) +  (1.03245313e-04)*pow(piEnergy,2) +  (-1.43686684e-02)*(piEnergy) + (4.48736544e-02);
     return weight;
   }
   if(piEnergy > 1000) { return 0.9271;}
 }
 else{return 1.0;}
}


//======Fate weight using tables =====================================================
//Lead Resonant
double CVUniverse::getLeadResPionTable(double piEnergy, int FateNum) const{
    //Inelastic
    if(FateNum == 1 || FateNum == 3){
        if(piEnergy < 0) { return 1.0; }
	    if(piEnergy > 0 && piEnergy <= 40) { return 1.5302; }
        if(piEnergy > 40 && piEnergy <= 80) { return 1.1878; }
        if(piEnergy > 80 && piEnergy <= 120) { return 0.5296; }
        if(piEnergy > 120 && piEnergy <= 160) { return 0.2770; }
        if(piEnergy > 160 && piEnergy <= 200) { return 0.2139; }
        if(piEnergy > 200 && piEnergy <= 240) { return 0.2668; }
        if(piEnergy > 240 && piEnergy <= 280) { return 0.3606; }
        if(piEnergy > 280 && piEnergy <= 320) { return 0.4544; }
        if(piEnergy > 320 && piEnergy <= 360) { return 0.5363; }
        if(piEnergy > 360 && piEnergy <= 400) { return 0.6329; }
        if(piEnergy > 400 && piEnergy <= 440) { return 0.6646; }
        if(piEnergy > 440 && piEnergy <= 480) { return 0.6981; }
        if(piEnergy > 480 && piEnergy <= 520) { return 0.7715; }
        if(piEnergy > 520 && piEnergy <= 560) { return 0.7495; }
        if(piEnergy > 560 && piEnergy <= 600) { return 0.7392; }
        if(piEnergy > 600 && piEnergy <= 640) { return 0.7734; }
        if(piEnergy > 640 && piEnergy <= 680) { return 0.8445; }
        if(piEnergy > 680 && piEnergy <= 720) { return 0.8377; }
        if(piEnergy > 720 && piEnergy <= 760) { return 0.8541; }
        if(piEnergy > 760 && piEnergy <= 800) { return 0.8397; }
        if(piEnergy > 800 && piEnergy <= 840) { return 0.8336; }
        if(piEnergy > 840 && piEnergy <= 880) { return 0.8344; }
        if(piEnergy > 880 && piEnergy <= 920) { return 0.8817; }
        if(piEnergy > 920 && piEnergy <= 960) { return 0.7635; }
        if(piEnergy > 960 && piEnergy <= 1000) { return 0.7650; }
        if(piEnergy > 1000 && piEnergy <= 1040) { return 0.8510;}
	    if(piEnergy > 1040){ return 0.8345; }
    }
    else if(FateNum == 2){ //CEX
        if(piEnergy < 0) { return 1.0; } 
        if(piEnergy > 0 && piEnergy <= 40) { return 0.9527; }
        if(piEnergy > 40 && piEnergy <= 80) { return 0.7160; }
        if(piEnergy > 80 && piEnergy <= 120) { return 2.7969; }
        if(piEnergy > 120 && piEnergy <= 160) { return 3.2923; }
        if(piEnergy > 160 && piEnergy <= 200) { return 3.7592; }
        if(piEnergy > 200 && piEnergy <= 240) { return 2.6846; }
        if(piEnergy > 240 && piEnergy <= 280) { return 5.5597; }
        if(piEnergy > 280 && piEnergy <= 320) { return 4.8354; }
        if(piEnergy > 320 && piEnergy <= 360) { return 4.5724; }
        if(piEnergy > 360 && piEnergy <= 400) { return 5.6335; }
        if(piEnergy > 400 && piEnergy <= 440) { return 4.1540; }
        if(piEnergy > 440 && piEnergy <= 480) { return 1.4223; }
        if(piEnergy > 480 && piEnergy <= 520) { return 0.4668; }
        if(piEnergy > 520 && piEnergy <= 560) { return 1.000; }
        if(piEnergy > 560 && piEnergy <= 600) { return 1.000; }
        if(piEnergy > 600 && piEnergy <= 640) { return 1.000; }
        if(piEnergy > 640 && piEnergy <= 680) { return 1.000; }
        if(piEnergy > 680 && piEnergy <= 720) { return 1.000; }
        if(piEnergy > 720 && piEnergy <= 760) { return 1.000; }
        if(piEnergy > 760 && piEnergy <= 800) { return 1.000; }
        if(piEnergy > 800 && piEnergy <= 840) { return 1.000; }
        if(piEnergy > 840 && piEnergy <= 880) { return 1.000; }
        if(piEnergy > 880 && piEnergy <= 920) { return 1.000; }
        if(piEnergy > 920 && piEnergy <= 960) { return 1.000; }
        if(piEnergy > 960 && piEnergy <= 1000) { return 1.000; }
	    if(piEnergy > 1000){ return 1.0; }
    }
   else if(FateNum == 8){ //Pion Prod
        if(piEnergy < 0){ return 1.0;}
        if(piEnergy > 0 && piEnergy <= 40) { return 0.000; }
        if(piEnergy > 40 && piEnergy <= 80) { return 0.000; }
        if(piEnergy > 80 && piEnergy <= 120) { return 0.000; }
        if(piEnergy > 120 && piEnergy <= 160) { return 0.000; }
        if(piEnergy > 160 && piEnergy <= 200) { return 0.000; }
        if(piEnergy > 200 && piEnergy <= 240) { return 0.000; }
        if(piEnergy > 240 && piEnergy <= 280) { return 0.000; }
        if(piEnergy > 280 && piEnergy <= 320) { return 0.000; }
        if(piEnergy > 320 && piEnergy <= 360) { return 6.5642; }
        if(piEnergy > 360 && piEnergy <= 400) { return 5.1401; }
        if(piEnergy > 400 && piEnergy <= 440) { return 5.9052; }
        if(piEnergy > 440 && piEnergy <= 480) { return 7.8059; }
        if(piEnergy > 480 && piEnergy <= 520) { return 7.0154; }
        if(piEnergy > 520 && piEnergy <= 560) { return 6.5287; }
        if(piEnergy > 560 && piEnergy <= 600) { return 5.9629; }
        if(piEnergy > 600 && piEnergy <= 640) { return 5.3695; }
        if(piEnergy > 640 && piEnergy <= 680) { return 4.9282; }
        if(piEnergy > 680 && piEnergy <= 720) { return 4.5877; }
        if(piEnergy > 720 && piEnergy <= 760) { return 4.0836; }
        if(piEnergy > 760 && piEnergy <= 800) { return 4.2895; }
        if(piEnergy > 800 && piEnergy <= 840) { return 4.0582; }
        if(piEnergy > 840 && piEnergy <= 880) { return 4.0318; }
        if(piEnergy > 880 && piEnergy <= 920) { return 3.8852; }
        if(piEnergy > 920 && piEnergy <= 960) { return 3.9641; }
        if(piEnergy > 960 && piEnergy <= 1000) { return 4.2284; }
        if(piEnergy > 1000 && piEnergy <= 1040) { return 3.9596;}
	    if(piEnergy > 1040){ return 4.0499; }
    }
   else if(FateNum == 5 || FateNum == 50 || FateNum == 51 || FateNum == 52 || FateNum == 9 ){ //Absorbtion
        if(piEnergy < 0){ return 1.0;}
        if(piEnergy > 0 && piEnergy <= 40) { return 1.3881; }
        if(piEnergy > 40 && piEnergy <= 80) { return 1.9619; }
        if(piEnergy > 80 && piEnergy <= 120) { return 2.0162; }
        if(piEnergy > 120 && piEnergy <= 160) { return 1.6898; }
        if(piEnergy > 160 && piEnergy <= 200) { return 1.5069; }
        if(piEnergy > 200 && piEnergy <= 240) { return 1.2925; }
        if(piEnergy > 240 && piEnergy <= 280) { return 1.4535; }
        if(piEnergy > 280 && piEnergy <= 320) { return 1.3605; }
        if(piEnergy > 320 && piEnergy <= 360) { return 1.1812; }
        if(piEnergy > 360 && piEnergy <= 400) { return 1.1997; }
        if(piEnergy > 400 && piEnergy <= 440) { return 1.2080; }
        if(piEnergy > 440 && piEnergy <= 480) { return 1.1704; }
        if(piEnergy > 480 && piEnergy <= 520) { return 0.9966; }
        if(piEnergy > 520 && piEnergy <= 560) { return 0.8951; }
        if(piEnergy > 560 && piEnergy <= 600) { return 0.8472; }
        if(piEnergy > 600 && piEnergy <= 640) { return 0.7971; }
        if(piEnergy > 640 && piEnergy <= 680) { return 0.6743; }
        if(piEnergy > 680 && piEnergy <= 720) { return 0.6773; }
        if(piEnergy > 720 && piEnergy <= 760) { return 0.6616; }
        if(piEnergy > 760 && piEnergy <= 800) { return 0.6208; }
        if(piEnergy > 800 && piEnergy <= 840) { return 0.5678; }
        if(piEnergy > 840 && piEnergy <= 880) { return 0.5848; }
        if(piEnergy > 880 && piEnergy <= 920) { return 0.5472; }
        if(piEnergy > 920 && piEnergy <= 960) { return 0.5600; }
        if(piEnergy > 960 && piEnergy <= 1000) { return 0.5379; }
        if(piEnergy > 1000 && piEnergy <= 1040) { return 0.5367;}
	    if(piEnergy > 1040){ return 0.5103 ; }
    }
    else if(FateNum == 10 ){ //Inelastic (KnockOut)
        if(piEnergy < 0){ return 1.0; }
	    if(piEnergy > 0 && piEnergy <= 40) { return 0.1051; }
        if(piEnergy > 40 && piEnergy <= 80) { return 0.1737; }
        if(piEnergy > 80 && piEnergy <= 120) { return 0.5331; }
        if(piEnergy > 120 && piEnergy <= 160) { return 0.8109; }
        if(piEnergy > 160 && piEnergy <= 200) { return 0.8782; }
        if(piEnergy > 200 && piEnergy <= 240) { return 1.4738; }
        if(piEnergy > 240 && piEnergy <= 280) { return 0.6700; }
        if(piEnergy > 280 && piEnergy <= 320) { return 0.7411; }
        if(piEnergy > 320 && piEnergy <= 360) { return 0.7595; }
        if(piEnergy > 360 && piEnergy <= 400) { return 0.4457; }
        if(piEnergy > 400 && piEnergy <= 440) { return 0.4367; }
        if(piEnergy > 440 && piEnergy <= 480) { return 0.4283; }
        if(piEnergy > 480 && piEnergy <= 520) { return 0.4198; }
        if(piEnergy > 520 && piEnergy <= 560) { return 0.3875; }
        if(piEnergy > 560 && piEnergy <= 600) { return 0.3875; }
        if(piEnergy > 600 && piEnergy <= 640) { return 0.3391; }
        if(piEnergy > 640 && piEnergy <= 680) { return 0.3356; }
        if(piEnergy > 680 && piEnergy <= 720) { return 0.2848; }
        if(piEnergy > 720 && piEnergy <= 760) { return 0.2438; }
        if(piEnergy > 760 && piEnergy <= 800) { return 0.2156; }
        if(piEnergy > 800 && piEnergy <= 840) { return 0.2736; }
        if(piEnergy > 840 && piEnergy <= 880) { return 0.2240; }
        if(piEnergy > 880 && piEnergy <= 920) { return 0.2647; }
        if(piEnergy > 920 && piEnergy <= 960) { return 0.2660; }
        if(piEnergy > 960 && piEnergy <= 1000) { return 0.2208; }
        if(piEnergy > 1000 && piEnergy <= 1040) { return 0.2430;}
	    if(piEnergy > 1040){ return 0.2497; }
    }
   else {return 1.0;}
}

double CVUniverse::getLeadQETable(double piEnergy, int FateNum) const{
    //Inelastic and No Scattering
    if(FateNum == 1 || FateNum == 3){
        if(piEnergy < 0) { return 1.0; }
        if(piEnergy > 0 && piEnergy <= 40) { return 0.4109; }
        if(piEnergy > 40 && piEnergy <= 80) { return 0.5353; }
        if(piEnergy > 80 && piEnergy <= 120) { return 0.7266; }
        if(piEnergy > 120 && piEnergy <= 160) { return 0.7333; }
        if(piEnergy > 160 && piEnergy <= 200) { return 0.7184; }
        if(piEnergy > 200 && piEnergy <= 240) { return 0.7319; }
        if(piEnergy > 240 && piEnergy <= 280) { return 0.7067; }
        if(piEnergy > 280 && piEnergy <= 320) { return 0.6983; }
        if(piEnergy > 320 && piEnergy <= 360) { return 0.7384; }
        if(piEnergy > 360 && piEnergy <= 400) { return 0.6945; }
        if(piEnergy > 400 && piEnergy <= 440) { return 0.7070; }
        if(piEnergy > 440 && piEnergy <= 480) { return 0.7109; }
        if(piEnergy > 480 && piEnergy <= 520) { return 0.6623; }
        if(piEnergy > 520 && piEnergy <= 560) { return 0.6553; }
        if(piEnergy > 560 && piEnergy <= 600) { return 0.6821; }
        if(piEnergy > 600 && piEnergy <= 640) { return 0.7106; }
        if(piEnergy > 640 && piEnergy <= 680) { return 0.6889; }
        if(piEnergy > 680 && piEnergy <= 720) { return 0.7471; }
        if(piEnergy > 720 && piEnergy <= 760) { return 0.7519; }
        if(piEnergy > 760 && piEnergy <= 800) { return 0.7763; }
        if(piEnergy > 800 && piEnergy <= 840) { return 0.7046; }
        if(piEnergy > 840 && piEnergy <= 880) { return 0.7772; }
        if(piEnergy > 880 && piEnergy <= 920) { return 0.8352; }
        if(piEnergy > 920 && piEnergy <= 960) { return 0.7915; }
        if(piEnergy > 960 && piEnergy <= 1000) { return 0.7896; }
        if(piEnergy > 1000 && piEnergy <= 1040) { return 0.8825; }
        if(piEnergy > 1040){ return 0.8216; }
        }
    else if(FateNum == 2 || FateNum == 10){ //KnockOut (Inelastic and CEX)
        if(piEnergy < 0) { return 1.0; }
        if(piEnergy > 0 && piEnergy <= 40) { return 1.9333; }
        if(piEnergy > 40 && piEnergy <= 80) { return 1.3265; }
        if(piEnergy > 80 && piEnergy <= 120) { return 1.1852; }
        if(piEnergy > 120 && piEnergy <= 160) { return 1.1793; }
        if(piEnergy > 160 && piEnergy <= 200) { return 1.2160; }
        if(piEnergy > 200 && piEnergy <= 240) { return 1.1920; }
        if(piEnergy > 240 && piEnergy <= 280) { return 1.2850; }
        if(piEnergy > 280 && piEnergy <= 320) { return 1.2707; }
        if(piEnergy > 320 && piEnergy <= 360) { return 1.2726; }
        if(piEnergy > 360 && piEnergy <= 400) { return 1.2188; }
        if(piEnergy > 400 && piEnergy <= 440) { return 1.2153; }
        if(piEnergy > 440 && piEnergy <= 480) { return 1.2360; }
        if(piEnergy > 480 && piEnergy <= 520) { return 1.3037; }
        if(piEnergy > 520 && piEnergy <= 560) { return 1.2681; }
        if(piEnergy > 560 && piEnergy <= 600) { return 1.2016; }
        if(piEnergy > 600 && piEnergy <= 640) { return 1.2025; }
        if(piEnergy > 640 && piEnergy <= 680) { return 1.1988; }
        if(piEnergy > 680 && piEnergy <= 720) { return 1.1515; }
        if(piEnergy > 720 && piEnergy <= 760) { return 1.1505; }
        if(piEnergy > 760 && piEnergy <= 800) { return 1.0806; }
        if(piEnergy > 800 && piEnergy <= 840) { return 1.1138; }
        if(piEnergy > 840 && piEnergy <= 880) { return 1.2427; }
        if(piEnergy > 880 && piEnergy <= 920) { return 1.0180; }
        if(piEnergy > 920 && piEnergy <= 960) { return 1.0322; }
        if(piEnergy > 960 && piEnergy <= 1000) { return 1.1228; }
        if(piEnergy > 1000 && piEnergy <= 1040) { return 1.0682; }
        if(piEnergy > 1040){ return 1.1292; }
         }
    else if(FateNum == 5 || FateNum == 50 || FateNum == 51 || FateNum == 52 || FateNum == 9){ //Absorption
        if(piEnergy < 0) { return 1.0; }
        if(piEnergy > 0 && piEnergy <= 40) { return 1.1615; }
        if(piEnergy > 40 && piEnergy <= 80) { return 1.3456; }
        if(piEnergy > 80 && piEnergy <= 120) { return 1.1585; }
        if(piEnergy > 120 && piEnergy <= 160) { return 1.2078; }
        if(piEnergy > 160 && piEnergy <= 200) { return 1.2314; }
        if(piEnergy > 200 && piEnergy <= 240) { return 1.2516; }
        if(piEnergy > 240 && piEnergy <= 280) { return 1.1977; }
        if(piEnergy > 280 && piEnergy <= 320) { return 1.2319; }
        if(piEnergy > 320 && piEnergy <= 360) { return 1.1743; }
        if(piEnergy > 360 && piEnergy <= 400) { return 1.2825; }
        if(piEnergy > 400 && piEnergy <= 440) { return 1.2688; }
        if(piEnergy > 440 && piEnergy <= 480) { return 1.1895; }
        if(piEnergy > 480 && piEnergy <= 520) { return 1.1497; }
        if(piEnergy > 520 && piEnergy <= 560) { return 1.1860; }
        if(piEnergy > 560 && piEnergy <= 600) { return 1.1819; }
        if(piEnergy > 600 && piEnergy <= 640) { return 1.1401; }
        if(piEnergy > 640 && piEnergy <= 680) { return 1.2157; }
        if(piEnergy > 680 && piEnergy <= 720) { return 1.0525; }
        if(piEnergy > 720 && piEnergy <= 760) { return 1.0643; }
        if(piEnergy > 760 && piEnergy <= 800) { return 1.1194; }
        if(piEnergy > 800 && piEnergy <= 840) { return 1.1174; }
        if(piEnergy > 840 && piEnergy <= 880) { return 0.9660; }
        if(piEnergy > 880 && piEnergy <= 920) { return 1.0066; }
        if(piEnergy > 920 && piEnergy <= 960) { return 1.0866; }
        if(piEnergy > 960 && piEnergy <= 1000) { return 0.9786; }
        if(piEnergy > 1000 && piEnergy <= 1040) { return 0.9826; }
        if(piEnergy > 1040){ return 0.9502; }   
        }
    else if(FateNum == 8){ //Pion Prod
        if(piEnergy < 0) { return 1.0; }
        if(piEnergy > 0 && piEnergy <= 40) { return 0.0; }
        if(piEnergy > 40 && piEnergy <= 80) { return 0.0; }
        if(piEnergy > 80 && piEnergy <= 120) { return 0.0; }
        if(piEnergy > 120 && piEnergy <= 160) { return 0.0; }
        if(piEnergy > 160 && piEnergy <= 200) { return 0.0; }
        if(piEnergy > 200 && piEnergy <= 240) { return 0.0; }
        if(piEnergy > 240 && piEnergy <= 280) { return 0.0; }
        if(piEnergy > 280 && piEnergy <= 320) { return 0.8362; }
        if(piEnergy > 320 && piEnergy <= 360) { return 0.2129; }
        if(piEnergy > 360 && piEnergy <= 400) { return 0.8680; }
        if(piEnergy > 400 && piEnergy <= 440) { return 0.9886; }
        if(piEnergy > 440 && piEnergy <= 480) { return 1.1805; }
        if(piEnergy > 480 && piEnergy <= 520) { return 1.2430; }
        if(piEnergy > 520 && piEnergy <= 560) { return 1.2580; }
        if(piEnergy > 560 && piEnergy <= 600) { return 1.2242; }
        if(piEnergy > 600 && piEnergy <= 640) { return 1.1646; }
        if(piEnergy > 640 && piEnergy <= 680) { return 1.0499; }
        if(piEnergy > 680 && piEnergy <= 720) { return 1.2001; }
        if(piEnergy > 720 && piEnergy <= 760) { return 1.1765; }
        if(piEnergy > 760 && piEnergy <= 800) { return 1.0871; }
        if(piEnergy > 800 && piEnergy <= 840) { return 1.1818; }
        if(piEnergy > 840 && piEnergy <= 880) { return 1.0432; }
        if(piEnergy > 880 && piEnergy <= 920) { return 1.1787; }
        if(piEnergy > 920 && piEnergy <= 960) { return 1.1501; }
        if(piEnergy > 960 && piEnergy <= 1000) { return 1.1080; }
        if(piEnergy > 1000 && piEnergy <= 1040) { return 1.0464; }
        if(piEnergy > 1040){ return 1.1157; }
    }
    else {return 1.0;}
    
}
    
//Carbon
double CVUniverse::getCarbonResPionTable(double piEnergy, int FateNum) const{
    //Elastic or No Scattering
    if(FateNum == 1 || FateNum == 3){
        if(piEnergy < 0) { return 1.0; }
	    if(piEnergy > 0 && piEnergy <= 40) { return 1.0581; }
        if(piEnergy > 40 && piEnergy <= 80) { return 1.0508; }
        if(piEnergy > 80 && piEnergy <= 120) { return 0.9222; }
        if(piEnergy > 120 && piEnergy <= 160) { return 0.6904; }
        if(piEnergy > 160 && piEnergy <= 200) { return 0.5873; }
        if(piEnergy > 200 && piEnergy <= 240) { return 0.6534; }
        if(piEnergy > 240 && piEnergy <= 280) { return 0.7737; }
        if(piEnergy > 280 && piEnergy <= 320) { return 0.8858; }
        if(piEnergy > 320 && piEnergy <= 360) { return 0.9561; }
        if(piEnergy > 360 && piEnergy <= 400) { return 0.9974; }
        if(piEnergy > 400 && piEnergy <= 440) { return 1.0238; }
        if(piEnergy > 440 && piEnergy <= 480) { return 1.0319; }
        if(piEnergy > 480 && piEnergy <= 520) { return 1.0375; }
        if(piEnergy > 520 && piEnergy <= 560) { return 1.0575; }
        if(piEnergy > 560 && piEnergy <= 600) { return 1.0363; }
        if(piEnergy > 600 && piEnergy <= 640) { return 1.0274; }
        if(piEnergy > 640 && piEnergy <= 680) { return 1.0539; }
        if(piEnergy > 680 && piEnergy <= 720) { return 1.0555; }
        if(piEnergy > 720 && piEnergy <= 760) { return 1.0454; }
        if(piEnergy > 760 && piEnergy <= 800) { return 1.0665; }
        if(piEnergy > 800 && piEnergy <= 840) { return 1.0530; }
        if(piEnergy > 840 && piEnergy <= 880) { return 1.0400; }
        if(piEnergy > 880 && piEnergy <= 920) { return 1.1024; }
        if(piEnergy > 920 && piEnergy <= 960) { return 1.0211; }
        if(piEnergy > 960 && piEnergy <= 1000) { return 1.0509; }
        if(piEnergy > 1000 && piEnergy <= 1040) { return 0.9644; }
	    if(piEnergy > 1040){ return 1.0169; }
    }
    else if(FateNum == 2){ //CEX
        if(piEnergy < 0) { return 1.0; } 
        if(piEnergy > 0 && piEnergy <= 40) { return 3.4713; }
        if(piEnergy > 40 && piEnergy <= 80) { return 1.6215; }
        if(piEnergy > 80 && piEnergy <= 120) { return 3.3349; }
        if(piEnergy > 120 && piEnergy <= 160) { return 3.8336; }
        if(piEnergy > 160 && piEnergy <= 200) { return 4.0685; }
        if(piEnergy > 200 && piEnergy <= 240) { return 2.0837; }
        if(piEnergy > 240 && piEnergy <= 280) { return 5.5453; }
        if(piEnergy > 280 && piEnergy <= 320) { return 4.5974; }
        if(piEnergy > 320 && piEnergy <= 360) { return 4.2784; }
        if(piEnergy > 360 && piEnergy <= 400) { return 5.8548; }
        if(piEnergy > 400 && piEnergy <= 440) { return 3.5994; }
        if(piEnergy > 440 && piEnergy <= 480) { return 0.1504; }
        if(piEnergy > 480 && piEnergy <= 520) { return 1.000; }
        if(piEnergy > 520 && piEnergy <= 560) { return 1.000; }
        if(piEnergy > 560 && piEnergy <= 600) { return 1.000; }
        if(piEnergy > 600 && piEnergy <= 640) { return 1.000; }
        if(piEnergy > 640 && piEnergy <= 680) { return 1.000; }
        if(piEnergy > 680 && piEnergy <= 720) { return 1.000; }
        if(piEnergy > 720 && piEnergy <= 760) { return 1.000; }
        if(piEnergy > 760 && piEnergy <= 800) { return 1.000; }
        if(piEnergy > 800 && piEnergy <= 840) { return 1.000; }
        if(piEnergy > 840 && piEnergy <= 880) { return 1.000; }
        if(piEnergy > 880 && piEnergy <= 920) { return 1.000; }
        if(piEnergy > 920 && piEnergy <= 960) { return 1.000; }
        if(piEnergy > 960 && piEnergy <= 1000) { return 1.000; }
	    if(piEnergy > 1000){ return 1.0; }
    }
   else if(FateNum == 8){ //Pion Prod
        if(piEnergy < 0){ return 1.0;}
        if(piEnergy > 0 && piEnergy <= 40) { return 0.0; }
        if(piEnergy > 40 && piEnergy <= 80) { return 0.0; }
        if(piEnergy > 80 && piEnergy <= 120) { return 0.0; }
        if(piEnergy > 120 && piEnergy <= 160) { return 0.0; }
        if(piEnergy > 160 && piEnergy <= 200) { return 0.0; }
        if(piEnergy > 200 && piEnergy <= 240) { return 0.0; }
        if(piEnergy > 240 && piEnergy <= 280) { return 0.0; }
        if(piEnergy > 280 && piEnergy <= 320) { return 0.0; }
        if(piEnergy > 320 && piEnergy <= 360) { return 2.7737; }
        if(piEnergy > 360 && piEnergy <= 400) { return 4.4940; }
        if(piEnergy > 400 && piEnergy <= 440) { return 5.5121; }
        if(piEnergy > 440 && piEnergy <= 480) { return 7.5147; }
        if(piEnergy > 480 && piEnergy <= 520) { return 6.6067; }
        if(piEnergy > 520 && piEnergy <= 560) { return 5.4925; }
        if(piEnergy > 560 && piEnergy <= 600) { return 5.7969; }
        if(piEnergy > 600 && piEnergy <= 640) { return 5.1540; }
        if(piEnergy > 640 && piEnergy <= 680) { return 4.6858; }
        if(piEnergy > 680 && piEnergy <= 720) { return 4.2884; }
        if(piEnergy > 720 && piEnergy <= 760) { return 4.2028; }
        if(piEnergy > 760 && piEnergy <= 800) { return 3.8684; }
        if(piEnergy > 800 && piEnergy <= 840) { return 4.0217; }
        if(piEnergy > 840 && piEnergy <= 880) { return 4.1118; }
        if(piEnergy > 880 && piEnergy <= 920) { return 3.6247; }
        if(piEnergy > 920 && piEnergy <= 960) { return 3.5825; }
        if(piEnergy > 960 && piEnergy <= 1000) { return 3.6775; }
        if(piEnergy > 1000 && piEnergy <= 1040) { return 4.6379; }
	    if(piEnergy > 1040){ return 4.1299; }
    }
   else if(FateNum == 5 || FateNum == 50 || FateNum == 51 || FateNum == 52 || FateNum == 9 ){ //Absorbtion
        if(piEnergy < 0){ return 1.0;}
        if(piEnergy > 0 && piEnergy <= 40) { return 0.9566; }
        if(piEnergy > 40 && piEnergy <= 80) { return 1.1902; }
        if(piEnergy > 80 && piEnergy <= 120) { return 1.2035; }
        if(piEnergy > 120 && piEnergy <= 160) { return 1.1330; }
        if(piEnergy > 160 && piEnergy <= 200) { return 0.9957; }
        if(piEnergy > 200 && piEnergy <= 240) { return 0.9289; }
        if(piEnergy > 240 && piEnergy <= 280) { return 0.9207; }
        if(piEnergy > 280 && piEnergy <= 320) { return 0.8339; }
        if(piEnergy > 320 && piEnergy <= 360) { return 0.7017; }
        if(piEnergy > 360 && piEnergy <= 400) { return 0.6756; }
        if(piEnergy > 400 && piEnergy <= 440) { return 0.6497; }
        if(piEnergy > 440 && piEnergy <= 480) { return 0.6289; }
        if(piEnergy > 480 && piEnergy <= 520) { return 0.5924; }
        if(piEnergy > 520 && piEnergy <= 560) { return 0.5072; }
        if(piEnergy > 560 && piEnergy <= 600) { return 0.4994; }
        if(piEnergy > 600 && piEnergy <= 640) { return 0.4328; }
        if(piEnergy > 640 && piEnergy <= 680) { return 0.4104; }
        if(piEnergy > 680 && piEnergy <= 720) { return 0.3860; }
        if(piEnergy > 720 && piEnergy <= 760) { return 0.3148; }
        if(piEnergy > 760 && piEnergy <= 800) { return 0.3967; }
        if(piEnergy > 800 && piEnergy <= 840) { return 0.3724; }
        if(piEnergy > 840 && piEnergy <= 880) { return 0.3596; }
        if(piEnergy > 880 && piEnergy <= 920) { return 0.3954; }
        if(piEnergy > 920 && piEnergy <= 960) { return 0.3572; }
        if(piEnergy > 960 && piEnergy <= 1000) { return 0.3527; }
        if(piEnergy > 1000 && piEnergy <= 1040) { return 0.3190; }
	    if(piEnergy > 1040){ return 0.3442; }
    }
    else if(FateNum == 10 ){ //Inelastic (KnockOut)
        if(piEnergy < 0){ return 1.0; }
	    if(piEnergy > 0 && piEnergy <= 40) { return 0.4748; }
        if(piEnergy > 40 && piEnergy <= 80) { return 0.5983; }
        if(piEnergy > 80 && piEnergy <= 120) { return 0.7132; }
        if(piEnergy > 120 && piEnergy <= 160) { return 0.9291; }
        if(piEnergy > 160 && piEnergy <= 200) { return 1.2356; }
        if(piEnergy > 200 && piEnergy <= 240) { return 1.7500; }
        if(piEnergy > 240 && piEnergy <= 280) { return 0.7447; }
        if(piEnergy > 280 && piEnergy <= 320) { return 0.7164; }
        if(piEnergy > 320 && piEnergy <= 360) { return 0.7285; }
        if(piEnergy > 360 && piEnergy <= 400) { return 0.4188; }
        if(piEnergy > 400 && piEnergy <= 440) { return 0.3973; }
        if(piEnergy > 440 && piEnergy <= 480) { return 0.3841; }
        if(piEnergy > 480 && piEnergy <= 520) { return 0.3839; }
        if(piEnergy > 520 && piEnergy <= 560) { return 0.3560; }
        if(piEnergy > 560 && piEnergy <= 600) { return 0.3335; }
        if(piEnergy > 600 && piEnergy <= 640) { return 0.3575; }
        if(piEnergy > 640 && piEnergy <= 680) { return 0.3148; }
        if(piEnergy > 680 && piEnergy <= 720) { return 0.2496; }
        if(piEnergy > 720 && piEnergy <= 760) { return 0.3096; }
        if(piEnergy > 760 && piEnergy <= 800) { return 0.2433; }
        if(piEnergy > 800 && piEnergy <= 840) { return 0.2508; }
        if(piEnergy > 840 && piEnergy <= 880) { return 0.2458; }
        if(piEnergy > 880 && piEnergy <= 920) { return 0.2060; }
        if(piEnergy > 920 && piEnergy <= 960) { return 0.2780; }
        if(piEnergy > 960 && piEnergy <= 1000) { return 0.2118; }
        if(piEnergy > 1000 && piEnergy <= 1040) { return 0.2763; }
	    if(piEnergy > 1040){ return 0.2210; }
    }
   else {return 1.0;}
}

double CVUniverse::getCarbonQETable(double piEnergy, int FateNum) const{
        //Inelastic and No Scattering
    if(FateNum == 1 || FateNum == 3){
        if(piEnergy < 0) { return 1.0; }
        if(piEnergy > 0 && piEnergy <= 40) { return 0.9040; }
        if(piEnergy > 40 && piEnergy <= 80) { return 1.0149; }
        if(piEnergy > 80 && piEnergy <= 120) { return 1.1158; }
        if(piEnergy > 120 && piEnergy <= 160) { return 1.0734; }
        if(piEnergy > 160 && piEnergy <= 200) { return 1.0356; }
        if(piEnergy > 200 && piEnergy <= 240) { return 1.0443; }
        if(piEnergy > 240 && piEnergy <= 280) { return 1.0311; }
        if(piEnergy > 280 && piEnergy <= 320) { return 1.0329; }
        if(piEnergy > 320 && piEnergy <= 360) { return 1.0184; }
        if(piEnergy > 360 && piEnergy <= 400) { return 1.0159; }
        if(piEnergy > 400 && piEnergy <= 440) { return 1.0147; }
        if(piEnergy > 440 && piEnergy <= 480) { return 1.0377; }
        if(piEnergy > 480 && piEnergy <= 520) { return 1.0206; }
        if(piEnergy > 520 && piEnergy <= 560) { return 0.9981; }
        if(piEnergy > 560 && piEnergy <= 600) { return 1.0288; }
        if(piEnergy > 600 && piEnergy <= 640) { return 0.9405; }
        if(piEnergy > 640 && piEnergy <= 680) { return 0.9591; }
        if(piEnergy > 680 && piEnergy <= 720) { return 1.0267; }
        if(piEnergy > 720 && piEnergy <= 760) { return 1.0484; }
        if(piEnergy > 760 && piEnergy <= 800) { return 1.0925; }
        if(piEnergy > 800 && piEnergy <= 840) { return 1.0623; }
        if(piEnergy > 840 && piEnergy <= 880) { return 0.9773; }
        if(piEnergy > 880 && piEnergy <= 920) { return 1.0922; }
        if(piEnergy > 920 && piEnergy <= 960) { return 1.0436; }
        if(piEnergy > 960 && piEnergy <= 1000) { return 1.0627; }
        if(piEnergy > 1000 && piEnergy <= 1040) { return 1.0894; }
        if(piEnergy > 1040){ return 1.0803; }
    }
    else if(FateNum == 2 || FateNum == 10){ //KnockOut (Inelastic and CEX)
        if(piEnergy < 0) { return 1.0; }
        if(piEnergy > 0 && piEnergy <= 40) { return 1.4330; }
        if(piEnergy > 40 && piEnergy <= 80) { return 0.9991; }
        if(piEnergy > 80 && piEnergy <= 120) { return 0.8865; }
        if(piEnergy > 120 && piEnergy <= 160) { return 0.8826; }
        if(piEnergy > 160 && piEnergy <= 200) { return 0.9436; }
        if(piEnergy > 200 && piEnergy <= 240) { return 0.9102; }
        if(piEnergy > 240 && piEnergy <= 280) { return 0.9715; }
        if(piEnergy > 280 && piEnergy <= 320) { return 0.9083; }
        if(piEnergy > 320 && piEnergy <= 360) { return 1.0017; }
        if(piEnergy > 360 && piEnergy <= 400) { return 0.9815; }
        if(piEnergy > 400 && piEnergy <= 440) { return 1.0245; }
        if(piEnergy > 440 && piEnergy <= 480) { return 0.9881; }
        if(piEnergy > 480 && piEnergy <= 520) { return 0.9183; }
        if(piEnergy > 520 && piEnergy <= 560) { return 1.0856; }
        if(piEnergy > 560 && piEnergy <= 600) { return 0.9630; }
        if(piEnergy > 600 && piEnergy <= 640) { return 1.0503; }
        if(piEnergy > 640 && piEnergy <= 680) { return 1.0772; }
        if(piEnergy > 680 && piEnergy <= 720) { return 0.9547; }
        if(piEnergy > 720 && piEnergy <= 760) { return 1.0746; }
        if(piEnergy > 760 && piEnergy <= 800) { return 0.9076; }
        if(piEnergy > 800 && piEnergy <= 840) { return 1.0750; }
        if(piEnergy > 840 && piEnergy <= 880) { return 0.9869; }
        if(piEnergy > 880 && piEnergy <= 920) { return 0.9529; }
        if(piEnergy > 920 && piEnergy <= 960) { return 1.0122; }
        if(piEnergy > 960 && piEnergy <= 1000) { return 0.8850; }
        if(piEnergy > 1000 && piEnergy <= 1040) { return 0.9490; }
        if(piEnergy > 1040){ return 0.9161; }
    }
    else if(FateNum == 5 || FateNum == 50 || FateNum == 51 || FateNum == 52 || FateNum == 9){ //Absorption
        if(piEnergy < 0) { return 1.0; }
        if(piEnergy > 0 && piEnergy <= 40) { return 0.7848; }
        if(piEnergy > 40 && piEnergy <= 80) { return 0.9578; }
        if(piEnergy > 80 && piEnergy <= 120) { return 0.8241; }
        if(piEnergy > 120 && piEnergy <= 160) { return 0.8911; }
        if(piEnergy > 160 && piEnergy <= 200) { return 0.9302; }
        if(piEnergy > 200 && piEnergy <= 240) { return 0.9291; }
        if(piEnergy > 240 && piEnergy <= 280) { return 0.9106; }
        if(piEnergy > 280 && piEnergy <= 320) { return 0.9764; }
        if(piEnergy > 320 && piEnergy <= 360) { return 0.9468; }
        if(piEnergy > 360 && piEnergy <= 400) { return 0.9700; }
        if(piEnergy > 400 && piEnergy <= 440) { return 0.9106; }
        if(piEnergy > 440 && piEnergy <= 480) { return 0.9204; }
        if(piEnergy > 480 && piEnergy <= 520) { return 1.0143; }
        if(piEnergy > 520 && piEnergy <= 560) { return 0.8994; }
        if(piEnergy > 560 && piEnergy <= 600) { return 0.9346; }
        if(piEnergy > 600 && piEnergy <= 640) { return 1.1141; }
        if(piEnergy > 640 && piEnergy <= 680) { return 1.0308; }
        if(piEnergy > 680 && piEnergy <= 720) { return 1.0590; }
        if(piEnergy > 720 && piEnergy <= 760) { return 0.8078; }
        if(piEnergy > 760 && piEnergy <= 800) { return 0.8713; }
        if(piEnergy > 800 && piEnergy <= 840) { return 0.9448; }
        if(piEnergy > 840 && piEnergy <= 880) { return 1.1493; }
        if(piEnergy > 880 && piEnergy <= 920) { return 0.8560; }
        if(piEnergy > 920 && piEnergy <= 960) { return 0.7919; }
        if(piEnergy > 960 && piEnergy <= 1000) { return 1.1637; }
        if(piEnergy > 1000 && piEnergy <= 1040) { return 0.7933; }
        if(piEnergy > 1040){ return 0.8312; }   
    }
    else if(FateNum == 8){ //Pion Prod
        if(piEnergy < 0) { return 1.0; }
        if(piEnergy > 0 && piEnergy <= 40) { return 0.0; }
        if(piEnergy > 40 && piEnergy <= 80) { return 0.0; }
        if(piEnergy > 80 && piEnergy <= 120) { return 0.0; }
        if(piEnergy > 120 && piEnergy <= 160) { return 0.0; }
        if(piEnergy > 160 && piEnergy <= 200) { return 0.0; }
        if(piEnergy > 200 && piEnergy <= 240) { return 0.0; }
        if(piEnergy > 240 && piEnergy <= 280) { return 0.0; }
        if(piEnergy > 280 && piEnergy <= 320) { return 0.3030; }
        if(piEnergy > 320 && piEnergy <= 360) { return 0.3459; }
        if(piEnergy > 360 && piEnergy <= 400) { return 0.7433; }
        if(piEnergy > 400 && piEnergy <= 440) { return 1.1505; }
        if(piEnergy > 440 && piEnergy <= 480) { return 0.6535; }
        if(piEnergy > 480 && piEnergy <= 520) { return 0.9714; }
        if(piEnergy > 520 && piEnergy <= 560) { return 1.1058; }
        if(piEnergy > 560 && piEnergy <= 600) { return 1.0284; }
        if(piEnergy > 600 && piEnergy <= 640) { return 1.1271; }
        if(piEnergy > 640 && piEnergy <= 680) { return 1.0572; }
        if(piEnergy > 680 && piEnergy <= 720) { return 0.8285; }
        if(piEnergy > 720 && piEnergy <= 760) { return 0.9396; }
        if(piEnergy > 760 && piEnergy <= 800) { return 0.9158; }
        if(piEnergy > 800 && piEnergy <= 840) { return 0.7362; }
        if(piEnergy > 840 && piEnergy <= 880) { return 0.9422; }
        if(piEnergy > 880 && piEnergy <= 920) { return 0.8819; }
        if(piEnergy > 920 && piEnergy <= 960) { return 1.0331; }
        if(piEnergy > 960 && piEnergy <= 1000) { return 0.8242; }
        if(piEnergy > 1000 && piEnergy <= 1040) { return 0.9525; }
        if(piEnergy > 1040){ return 0.9271; }
    }
    else {return 1.0;}
    
}

//=======================================================================
std::tuple<std::vector<int>,std::vector<int>,std::vector<double>,std::vector<double> > CVUniverse::calFateNumVector() const {

 int mc_er_nPart = GetInt("mc_er_nPart");
  double tempfate = -999;
  std::vector<int> tempfateV;
  std::vector<int> prefsiV;
  std::vector<double> prefsiEnV;
  std::vector<double> FDEnV;
  int prefsi=0;

  //std::cout<<" Here "<<std::endl;

  for(int i=0; i < mc_er_nPart; ++i){
    int mc_er_status = GetVecElem("mc_er_status",i);
    double prefsiEn = GetVecElem("mc_er_E",i);

    //std::cout<<" mc_er_status: "<<mc_er_status<<std::endl;
    if( mc_er_status == 14){
      tempfate = -1;

      // Only track fates for proton, neutron, pizero, piplus, piminus.
      int mc_er_ID = GetVecElem("mc_er_ID",i);
 	  if(mc_er_ID == 2212 || mc_er_ID == 2112 || mc_er_ID == 111 || mc_er_ID==211 || mc_er_ID ==-211){
        prefsi = mc_er_ID;
 		// The distinguishing feature is first daughter last daughter
        // get an easy handle to these for use later
        int fd = GetVecElem("mc_er_FD",i);
 	    int ld = GetVecElem("mc_er_LD",i);
 		int mot= GetVecElem("mc_er_mother",i);

 		//is one of these a Pion
 		int isaPion = 0;
 		if(mc_er_ID == 111 || mc_er_ID == 211 || mc_er_ID == -211) isaPion = 1;

        //How many daughthers are pions
		int FSpions = 0;
		for(int jj = fd; jj <= ld; ++jj){
          if( GetVecElem("mc_er_ID",jj) == 211 || GetVecElem("mc_er_ID",jj) == -211 || GetVecElem("mc_er_ID",jj) == 111)FSpions++;
        }

        /*
        if(ld-fd == 1){
            if(mc_er_ID == 111 || mc_er_ID == 211 || mc_er_ID == -211){
                if( GetVecElem("mc_er_ID",fd) == 111 || GetVecElem("mc_er_ID",fd) == 211 || GetVecElem("mc_er_ID",fd) == -211  ){
                    std::cout<<"mc_er_ID: "<<mc_er_ID<< " fd: " << GetVecElem("mc_er_ID",fd)<<" ld: "<<GetVecElem("mc_er_ID",ld) <<std::endl;
                }

            }
        }*/

        if(fd == ld){  //only one daughter: fate can be no interaction,elastic or absorption

          int mc_targetA = GetInt("mc_targetA");
          double offset  = getGenieBEinMeV(mc_targetA);
          double tolerance = 0.0001;//withing machine precision
          int mc_intType = GetInt("mc_intType");
          // unknown nuclei need a larger tolerance
          if(TMath::Abs(offset - 8.0) < 0.1)tolerance = 1.2;
          if(offset < 6.0)tolerance = 1.2;
          if(mc_intType != 1)offset = 0.0;

          if(GetInt("mc_intType") != 1)offset = 0.0;
          if(GetVecElem("mc_er_ID",fd) > 2000000000){
            tempfate = 5;   // absorption on 3 nucleons
            tempfateV.push_back(tempfate);
            prefsiV.push_back(prefsi);
            prefsiEnV.push_back(prefsiEn);
            FDEnV.push_back(GetVecElem("mc_er_E",fd));

            //result = getGenieEan(chw,entry,i,fd,ld);
            //printDebugFun(chw,entry,i,fd,ld,mot,tempfate);
          }
          else if(TMath::Abs( GetVecElem("mc_er_E",i) - GetVecElem("mc_er_E",fd) - offset) < tolerance){

            tempfate = 1;   // no scattering
            tempfateV.push_back(tempfate);
            prefsiV.push_back(prefsi);
            prefsiEnV.push_back(prefsiEn);
            FDEnV.push_back(GetVecElem("mc_er_E",fd));

            //result = getGenieEan(chw,entry,i,fd,ld);
            //printDebugFun(chw,entry,i,fd,ld,mot,tempfate);
        }else if(GetVecElem("mc_er_ID",fd) != GetVecElem("mc_er_ID",i)){
            cout << "Now Here some day?" << endl;
            tempfate = 2;   // charge exchange ?
            tempfateV.push_back(tempfate);
            prefsiV.push_back(prefsi);
            prefsiEnV.push_back(prefsiEn);
            FDEnV.push_back(GetVecElem("mc_er_E",fd));
            //std::cout<<" hola 2"<<std::endl;
            //result = getGenieEan(chw,entry,i,fd,ld);
            //printDebugFun(chw,entry,i,fd,ld,mot,tempfate);
        }else {
            tempfate = 3;   // elastic fate
            tempfateV.push_back(tempfate);
            prefsiV.push_back(prefsi);
            prefsiEnV.push_back(prefsiEn);
            FDEnV.push_back(GetVecElem("mc_er_E",fd));
            //std::cout<<" hola 3"<<std::endl;
            //getGenieEan(chw,entry,i,fd,ld);
            //printDebugFun(chw,entry,i,fd,ld,mot,offset);
        }
      }
      else if(FSpions > 0 && !isaPion){

          tempfate = 8; // nucleon goes to pions
          tempfateV.push_back(tempfate);
          prefsiV.push_back(prefsi);
          prefsiEnV.push_back(prefsiEn);
          FDEnV.push_back(GetVecElem("mc_er_E",fd));
          //std::cout<<" hola 8"<<std::endl;
          //result = getGenieEan(chw,entry,i,fd,ld);
          //printDebugFun(chw,entry,i,fd,ld,mot,tempfate);
      }
      else if(isaPion && FSpions == 0){

  		tempfate = 5; //two body absorption?

  		if (ld-fd ==1){
  		  int p=0;
  		  int n=0;
  		  if(GetVecElem("mc_er_ID",fd) == 2212)p++;
          if(GetVecElem("mc_er_ID",fd) == 2112)n++;
          if(GetVecElem("mc_er_ID",ld) == 2212)p++;
          if(GetVecElem("mc_er_ID",ld) == 2112)n++;
          // use special codes for producing nn, pn, pp
          if(p==1 && n==1)tempfate = 51;
          if(p==2)tempfate = 52;
          if(n==2)tempfate = 50;
  		}
        tempfateV.push_back(tempfate);
        prefsiV.push_back(prefsi);
        prefsiEnV.push_back(prefsiEn);
        FDEnV.push_back(GetVecElem("mc_er_E",fd));

  	  } else if(isaPion && FSpions >= 2){

          tempfate = 8; //nucleon goes to pions
          tempfateV.push_back(tempfate);
          prefsiV.push_back(prefsi);
          prefsiEnV.push_back(prefsiEn);
          FDEnV.push_back(GetVecElem("mc_er_E",fd));
          //std::cout<<" hola 8-"<<std::endl;
          //result = getGenieEan(chw,entry,i,fd,ld);
          //printDebugFun(chw,entry,i,fd,ld,mot,tempfate);
  		} else if(isaPion && FSpions == 1){
          for(int jj = fd; jj <= ld; ++jj){
            if((GetVecElem("mc_er_ID",jj) == 111 || TMath::Abs(GetVecElem("mc_er_ID",jj)) == 211) && GetVecElem("mc_er_ID",jj) != GetVecElem("mc_er_ID",i)){
                tempfate = 2; //charge exchange?
                }
          }
          if(tempfate == 2){
            tempfateV.push_back(tempfate);
            prefsiV.push_back(prefsi);
            prefsiEnV.push_back(prefsiEn);
            FDEnV.push_back(GetVecElem("mc_er_E",fd));
            //std::cout<<" hola 2-"<<std::endl;
          }
          if(tempfate != 2) {
            //Here lands the Inelastic fate for hA
            tempfate = 10;
            tempfateV.push_back(tempfate);
            prefsiV.push_back(prefsi);
            prefsiEnV.push_back(prefsiEn);
            FDEnV.push_back(GetVecElem("mc_er_E",fd));
            //std::cout<<" hola 4"<<std::endl;
          }
        } else {
          tempfate = 4;  // or four, can't tell. Other fate
          //result = getGenieEan(chw,entry,i,fd,ld);
          if(GetVecElem("mc_er_ID",fd) ==  2000000300){
             tempfate = 9; //Multinucleons
             tempfateV.push_back(tempfate);
             prefsiV.push_back(prefsi);
             prefsiEnV.push_back(prefsiEn);
             FDEnV.push_back(GetVecElem("mc_er_E",fd));
             //std::cout<<" hola 9"<<std::endl;
             //result = getGenieEan(chw,entry,i,fd,ld);
             //printDebugFun(chw,entry,i,fd,ld,mot,tempfate);
          }
          else if(ld-fd == 1 && tempfate != 9){
             tempfate = 10; //knock out
             tempfateV.push_back(tempfate);
             prefsiV.push_back(prefsi);
             prefsiEnV.push_back(prefsiEn);
             FDEnV.push_back(GetVecElem("mc_er_E",fd));
             //std::cout<<" hola 10"<<std::endl;
             //printDebugFun(chw,entry,i,fd,ld,mot,tempfate);
             //getGenieEan2(chw,entry,i,fd,ld);
          }
          else if (ld-fd == 1){
            std::cout <<"here"<< std::endl;
            if (GetVecElem("mc_er_ID",fd) == 211) std::cout <<"woow"<< std::endl;
            if (GetVecElem("mc_er_ID",fd) == -211) std::cout <<"woow"<< std::endl;
            if (GetVecElem("mc_er_ID",fd) == 111) std::cout <<"woow"<< std::endl;
            if (GetVecElem("mc_er_ID",ld) == 211) std::cout <<"woow"<< std::endl;
            if (GetVecElem("mc_er_ID",ld) == -211) std::cout <<"woow"<< std::endl;
            if (GetVecElem("mc_er_ID",ld) == 111) std::cout <<"woow"<< std::endl;

          }
        }
        if(tempfate == 4) {
                tempfate = 10;
                //result = getGenieEan(chw,entry,i,fd,ld);
                //printDebugFun(chw,entry,i,fd,ld);
        }

        //printDebugFun(chw,entry,i,fd,ld,mot,tempfate);
        //getGenieEan2(chw,entry,i,fd,ld);
  		//std::cout <<"Reco q3: "<< getRecoq3(chw,entry) << std::endl;
        //std::cout<< "tempfate: "<<tempfate<<" prefsi: "<<prefsi <<" prefsiEn: "<<prefsiEn<<" fd: "<<GetVecElem("mc_er_ID",fd)<<" fd energy: "<< GetVecElem("mc_er_E",fd)<<" ld: "<<GetVecElem("mc_er_ID",ld)<<std::endl;
  		    } //end fates for p,n,pi
    //std::cout <<"tempfate: "<<tempfate <<std::endl;
    } //End mc_er_status == 14

  } //End for loop int i=0; i < mc_er_nPart; ++i
  //std::cout<<" size vector inside function:  "<< tempfateV.size() << std::endl;
  //std::cout <<"size vect: "<< tempfateV.size() <<std::endl;
  return std::make_tuple(tempfateV, prefsiV, prefsiEnV,FDEnV);

}
//=======================================================================
int CVUniverse::calFateNum() const {

  int mc_er_nPart = GetInt("mc_er_nPart");
  double tempfate = -999;
  std::vector<double> result;
  int prefsi=0;

  for(int i=0; i < mc_er_nPart; ++i){

    int mc_er_status = GetVecElem("mc_er_status",i);
    if( mc_er_status == 14){
         tempfate = -1;

         // Only track fates for proton, neutron, pizero, piplus, piminus.
 	      int mc_er_ID = GetVecElem("mc_er_ID",i);
 	      if(mc_er_ID == 2212 || mc_er_ID == 2112 || mc_er_ID == 111 || mc_er_ID==211 || mc_er_ID ==-211){
            prefsi = mc_er_ID;
 		        // The distinguishing feature is first daughter last daughter
            // get an easy handle to these for use later
            int fd = GetVecElem("mc_er_FD",i);
 	          int ld = GetVecElem("mc_er_LD",i);
 		        int mot= GetVecElem("mc_er_mother",i);

 		//is one of these a Pion
 		int isaPion = 0;
 		if(mc_er_ID == 111 || mc_er_ID == 211 || mc_er_ID == -211) isaPion = 1;

    //How many daughthers are pions
		int FSpions = 0;
		for(int jj = fd; jj <= ld; ++jj){
                    if( GetVecElem("mc_er_ID",jj) == 211 || GetVecElem("mc_er_ID",jj) == -211 || GetVecElem("mc_er_ID",jj) == 111)FSpions++;
                }

    if(fd == ld){  //only one daughter: fate can be no interaction,elastic or absorption

        int mc_targetA = GetInt("mc_targetA");
        double offset  = getGenieBEinMeV(mc_targetA);
        double tolerance = 0.0001;//withing machine precision
        int mc_intType = GetInt("mc_intType");
        // unknown nuclei need a larger tolerance
        if(TMath::Abs(offset - 8.0) < 0.1)tolerance = 1.2;
            if(offset < 6.0)tolerance = 1.2;
        if(mc_intType != 1)offset = 0.0;

        if(GetInt("mc_intType") != 1)offset = 0.0;
        if(GetVecElem("mc_er_ID",fd) > 2000000000){

        tempfate = 5;   // absorption on 3 nucleons
        //result = getGenieEan(chw,entry,i,fd,ld);
        //printDebugFun(chw,entry,i,fd,ld,mot,tempfate);
        }
        else if(TMath::Abs( GetVecElem("mc_er_E",i) - GetVecElem("mc_er_E",fd) - offset) < tolerance){

            tempfate = 1;   // no scattering
            //result = getGenieEan(chw,entry,i,fd,ld);
            //printDebugFun(chw,entry,i,fd,ld,mot,tempfate);
        }else if(GetVecElem("mc_er_ID",fd) != GetVecElem("mc_er_ID",i)){
            cout << "Now Here some day?" << endl;
            tempfate = 2;   // charge exchange ?
            //result = getGenieEan(chw,entry,i,fd,ld);
            //printDebugFun(chw,entry,i,fd,ld,mot,tempfate);
        }else {
            tempfate = 3;   // elastic fate
            //getGenieEan(chw,entry,i,fd,ld);
            //printDebugFun(chw,entry,i,fd,ld,mot,offset);
        }
      }
      else if(FSpions > 0 && !isaPion){

          tempfate = 8; // nucleon goes to pions
          //result = getGenieEan(chw,entry,i,fd,ld);
          //printDebugFun(chw,entry,i,fd,ld,mot,tempfate);
      }
      else if(isaPion && FSpions == 0){

  			tempfate = 5; //two body absorption?

  			if (ld-fd ==1){
  				int p=0;
  				int n=0;
  				if(GetVecElem("mc_er_ID",fd) == 2212)p++;
                          	if(GetVecElem("mc_er_ID",fd) == 2112)n++;
                          	if(GetVecElem("mc_er_ID",ld) == 2212)p++;
                          	if(GetVecElem("mc_er_ID",ld) == 2112)n++;
                              	// use special codes for producing nn, pn, pp
                                  if(p==1 && n==1)tempfate = 51;
                                  if(p==2)tempfate = 52;
                                  if(n==2)tempfate = 50;
  			}

  		} else if(isaPion && FSpions >= 2){

                      tempfate = 8; //nucleon goes to pions
                      //result = getGenieEan(chw,entry,i,fd,ld);
                      //printDebugFun(chw,entry,i,fd,ld,mot,tempfate);
  		} else if(isaPion && FSpions == 1){

                      tempfate = 4;
                      for(int jj = fd; jj <= ld; ++jj){
                          if(( GetVecElem("mc_er_ID",jj) == 111 || TMath::Abs(GetVecElem("mc_er_ID",jj)) == 211) && GetVecElem("mc_er_ID",jj) != GetVecElem("mc_er_ID",i)) tempfate = 2; //charge exchange?
                      }
                  } else {
                     //std::cout << " Got unknown fate " << tempfate << std::endl;
  		                tempfate = 4;  // or four, can't tell. Other fate
                     //result = getGenieEan(chw,entry,i,fd,ld);
  		                if(GetVecElem("mc_er_ID",fd) ==  2000000300){
                      tempfate = 9; //Multinucleons
                      //result = getGenieEan(chw,entry,i,fd,ld);
  		               //printDebugFun(chw,entry,i,fd,ld,mot,tempfate);
  		   }
                     if(ld-fd == 1 ){ tempfate = 10; //knock out
                          //printDebugFun(chw,entry,i,fd,ld,mot,tempfate);
          	        //getGenieEan2(chw,entry,i,fd,ld);
  			//std::cout <<"Reco q3: "<< getRecoq3(chw,entry) << std::endl;
  		    }

  		  }
                      if(tempfate == 4) {
  			                tempfate =10;
                        //result = getGenieEan(chw,entry,i,fd,ld);
  		                 //printDebugFun(chw,entry,i,fd,ld);
  		     }

  		} //end fates for p,n,pi
  	       //std::cout<<" fate is: "<< tempfate << " prefsi: "<<prefsi <<std::endl;


    } //End mc_er_status == 14



  } //End for loop int i=0; i < mc_er_nPart; ++i

  return tempfate;

}
//=======================================================================
double CVUniverse::getGenieBEinMeV(int A) const{
  // From UserPhysicsOptions.xml file
  // these are hard coded, so we can get an exact match
  // but beware if you try to port this code to any new GENIE.

  if (A == 1) return 0; // hydrogen
  if (A == 6) return 17.0; // lithium
  if (A == 12) return 25.0; // carbon
  if (A == 16) return 27.0; // oxygen
  if (A == 24) return 32.0; // magnesium
  if (A == 40) return 29.5; // argon
  if (A == 48) return 30.0; // Ti48
  if (A == 56) return 36.0; // 56 iron
  if (A == 58) return 36.0; // 58 nickel
  if (A >= 206 && A <= 208) return 44.0; // 208 lead

  // Specialty in MINERvA,picked off numbers by hand.
  if (A == 28) return 8.219751;  // silicon
  if (A == 27) return 8.115287;  // aluminum
  if (A == 14) return 7.185166;  // nitrogen
  if (A == 55) return 8.653063;  // iron55 or manganese55
  if (A == 35) return 8.347164;  // chlorine

  if (A == 4) return 5.0;  // this is rough,

  //else
  // this is a problem, all other numbers come from a semi-empirical binding energy formula
  // need to back off the QE precision for those.
  return 8.0;   // GENIE defaults to something like this.  Needs to be exact.  Check it.
}
//=======================================================================
double CVUniverse::GetEavailReco() const{

 double magicEnergyScale=1.17;

  return magicEnergyScale;
}
//=======================================================================
double CVUniverse::GetEavailRecoECal() const{
   double magicEnergyScale = 1.17;
   double recoilE_tracker_nom = 1e-3*GetDouble("blob_recoil_E_tracker");
   double recoilE_ecal_nom = 1e-3*GetDouble("blob_recoil_E_ecal");
   double recoilE_hcal_nom = 1e-3*GetDouble("lowrecoil_hcal_calE_1");

   std::vector<double> FuzzEvts = GetVecDouble("muon_fuzz_per_plane_r80_energies");
   double mufuzzTrackerEcal = 0.0;

   for(unsigned int l=0; l<FuzzEvts.size(); l++){
    double mufuzzTotal_E = 1e-3*GetVecElem("muon_fuzz_per_plane_r80_energies",l);
    if(mufuzzTotal_E<0) mufuzzTotal_E=0;
    int mufuzzTrackerEcal_PlnID = GetVecElem("muon_fuzz_per_plane_r80_planeIDs",l);
    if(mufuzzTrackerEcal_PlnID < 1.850E+09){ //6 module HCAL
          mufuzzTrackerEcal += mufuzzTotal_E;
    }
  }
  bool subtractMuFuzz = true;
  if(mufuzzTrackerEcal <0) mufuzzTrackerEcal=0; // Protect against -99999 default values. Dunno how they got in
  // Protect against mufuzz>actual energy. Don't know how that happens either
  mufuzzTrackerEcal = std::min(mufuzzTrackerEcal, recoilE_tracker_nom+recoilE_ecal_nom+recoilE_hcal_nom);

  double recoilE_trackerEcal=(recoilE_tracker_nom+recoilE_ecal_nom+recoilE_hcal_nom) - (subtractMuFuzz ? mufuzzTrackerEcal : 0);
  double ret_oscar = magicEnergyScale*(recoilE_trackerEcal);

  if(ret_oscar<0){
    cerr << "Tracker+ECAL+HCAL energy is " << ret_oscar << ". Bailing" << endl;
    exit(1);
  }
  return ret_oscar;
}
//===================================Reco Available Energy for ECAL HCAL===================================
double CVUniverse::calcECALHCALAvEn() const{
  bool subtractMuFuzz=true;
  double magicEnergyScale=1.17;
  double recoilE_tracker_nom=0;
  double recoilE_ecal_nom=0;
  double recoilE_hcal_nom=0;

  recoilE_tracker_nom = 1e-3*GetDouble("blob_recoil_E_tracker");
  recoilE_ecal_nom    = 1e-3*GetDouble("blob_recoil_E_ecal");
  recoilE_hcal_nom    = 1e-3*GetDouble("lowrecoil_hcal_calE_1");

  std::vector<double> FuzzEvts = GetVecDouble("muon_fuzz_per_plane_r80_energies");
  double mufuzzTrackerEcal = 0.0;

  for(unsigned int l=0; l<FuzzEvts.size(); l++){
    double mufuzzTotal_E = 1e-3*GetVecElem("muon_fuzz_per_plane_r80_energies",l);
    if(mufuzzTotal_E<0) mufuzzTotal_E=0;
    int mufuzzTrackerEcal_PlnID = GetVecElem("muon_fuzz_per_plane_r80_planeIDs",l);
    if(mufuzzTrackerEcal_PlnID < 1.850E+09){ //6 module HCAL
          mufuzzTrackerEcal += mufuzzTotal_E;
    }
  }
  // Protect against -99999 default values. Dunno how they got in
  if(mufuzzTrackerEcal <0) mufuzzTrackerEcal=0;
  // Protect against mufuzz>actual energy. Don't know how that happens either
  mufuzzTrackerEcal = std::min(mufuzzTrackerEcal, recoilE_tracker_nom+recoilE_ecal_nom+recoilE_hcal_nom);

  double recoilE_trackerEcal=(recoilE_tracker_nom+recoilE_ecal_nom+recoilE_hcal_nom) - (subtractMuFuzz ? mufuzzTrackerEcal : 0);
  double ret = magicEnergyScale*(recoilE_trackerEcal);

  if(ret<0){
    cerr << "Tracker+ECAL+HCAL energy is " << ret << ". Bailing" << endl;
    exit(1);
  }
  return ret;
}
//===================================Reco Available Energy for Tracker ECAL===================================
double CVUniverse::calcTrackerECALAvEn() const {

  bool subtractMuFuzz=true;
  double magicEnergyScale=1.17;
  double recoilE_tracker_nom=0;
  double recoilE_ecal_nom=0;

  recoilE_tracker_nom = 1e-3*GetDouble("blob_recoil_E_tracker");
  recoilE_ecal_nom    = 1e-3*GetDouble("blob_recoil_E_ecal");

  std::vector<double> FuzzEvts = GetVecDouble("muon_fuzz_per_plane_r80_energies");
  double mufuzzTrackerEcal = 0.0;

  //Implementation for the right distance downstream the vertex (300 mm)
  double vtx;
  
  //double mintrack;
  //mintrack=GetDouble("minervaTrackZPosition");
  //mintrack=mintrack + 300;
  //int startingPlaneID2= Convert2PlaneID(mintrack);

  //vtx = GetVecElem("vtx",2);
  //vtx = vtx + 300;   
  //int startingPlaneID = Convert2PlaneID(vtx);


  for(unsigned int l=0; l<FuzzEvts.size(); l++){
    double mufuzzTotal_E = 1e-3*GetVecElem("muon_fuzz_per_plane_r80_energies",l);
    if(mufuzzTotal_E<0) mufuzzTotal_E=0;
    int mufuzzTrackerEcal_PlnID = GetVecElem("muon_fuzz_per_plane_r80_planeIDs",l);
    //if(mufuzzTrackerEcal_PlnID < 1844969472 && mufuzzTrackerEcal_PlnID > startingPlaneID ){ //LastECal module
    //if(mufuzzTrackerEcal_PlnID < 1844969472 && mufuzzTrackerEcal_PlnID > startingPlaneID2 ){ //LastECal module
    if(mufuzzTrackerEcal_PlnID < 1844969472){ //Use 1846018048 for First Module HCAL  //Use 1844969472 for last ECAL module
          mufuzzTrackerEcal += mufuzzTotal_E;
    }
  }
  // Protect against -99999 default values. Dunno how they got in
  if(mufuzzTrackerEcal <0) mufuzzTrackerEcal=0;

  //if(mufuzzTrackerEcal >  recoilE_tracker_nom+recoilE_ecal_nom ) mufuzzTrackerEcal = 0; //Oscar added this to protect when muon fuzz is to high
  // Protect against mufuzz>actual energy. Don't know how that happens either
  mufuzzTrackerEcal = std::min(mufuzzTrackerEcal, recoilE_tracker_nom+recoilE_ecal_nom);

  double recoilE_trackerEcal=(recoilE_tracker_nom+recoilE_ecal_nom) - (subtractMuFuzz ? mufuzzTrackerEcal : 0);
  double ret = magicEnergyScale*magicEnergyScale*(recoilE_trackerEcal);

  if(ret<0){
    cerr << "Tracker+ECAL energy is " << ret << ". Bailing" << endl;
    exit(1);
  }
  return ret;
}

//Only for MEC Ana tuples
double CVUniverse::MECAnacalcTrackerECALAvEn() const {

    double recoilE_tracker_nom=0;
    double recoilE_ecal_nom=0;
    double magicEnergyScale =1.17;
    bool subtractMuFuzz = true;

    recoilE_tracker_nom = 1e-3*GetDouble("blob_recoil_E_tracker");
    recoilE_ecal_nom    = 1e-3*GetDouble("blob_recoil_E_ecal");

    double mufuzz_tracker      = 1e-3*GetDouble("lrMuonFuzzEnergyTracker");
    double mufuzz_ecal         = 1e-3*GetDouble("lrMuonFuzzEnergyEcal");

    // Protect against -99999 default values. Dunno how they got in
    if(mufuzz_tracker<0) mufuzz_tracker=0;
    if(mufuzz_ecal<0)    mufuzz_ecal=0;

    // Protect against mufuzz>actual energy. Don't know how that happens either
    mufuzz_tracker = std::min(mufuzz_tracker, recoilE_tracker_nom);
    mufuzz_ecal    = std::min(mufuzz_ecal,    recoilE_ecal_nom);

    double recoilE_tracker=recoilE_tracker_nom - (subtractMuFuzz ? mufuzz_tracker : 0);
    double recoilE_ecal=   recoilE_ecal_nom    - (subtractMuFuzz ? mufuzz_ecal    : 0);

    double ret = magicEnergyScale*(recoilE_tracker+recoilE_ecal);

    if(ret<0){
        cerr << "Tracker+ECAL energy is " << ret << ". Bailing" << endl;
        exit(1);
    }
      return ret;
      //return 1.0; //ONLY FOR SPECIAL SITUATION
}

double CVUniverse::MECAnacalcRecoq3() const{

    return 1e-3*GetDouble("MECAna_q3");
    //return 1.0; //ONLY FOR SPECIAL SITUATION
}



//==================================================================================================
double CVUniverse::calcRecoq3() const{
  double Q2Reco = GetQ2Reco(); //MeV^2
  double EhadREco=GetRecoilEnergy();//MeV
  return (sqrt(Q2Reco + EhadREco*EhadREco)/1000);
}
//=================================================================================
double CVUniverse::calcCVWeight() const{
  double weight=GetDouble("mc_cvweight_total");
  if(std::isnan(weight)) {weight=1;}

  return weight;
}
//=================================================================================
double CVUniverse::Convert2PlaneID(double vtx) const{
  std::vector<std::pair<double, int>> vtxToPlaneID;

  vtxToPlaneID = {
            {4293.04, 1336147968}, {4313.68, 1336410112}, {4337.25, 1337196544},
            {4357.9, 1337458688}, {4381.47, 1338245120}, {4402.11, 1338507264},
            {4425.68, 1339293696}, {4446.33, 1339555840}, {4514.11, 1208221696},
            {4534.76, 1208483840}, {4558.33, 1209270272}, {4578.97, 1209532416},
            {4602.54, 1210318848}, {4623.19, 1210580992}, {4646.76, 1211367424},
            {4667.4, 1211629568}, {4735.19, 1213464576}, {4755.83, 1213726720},
            {4779.4, 1214513152}, {4800.05, 1214775296}, {4823.62, 1215561728},
            {4844.26, 1215823872}, {4867.83, 1216610304}, {4888.48, 1216872448},
            {5000.48, 1219756032}, {5021.12, 1220018176}, {5044.69, 1220804608},
            {5065.34, 1221066752}, {5088.91, 1221853184}, {5109.55, 1222115328},
            {5133.12, 1222901760}, {5153.77, 1223163904}, {5456.74, 1223950336},
            {5477.38, 1224212480}, {5500.95, 1224998912}, {5521.6, 1225261056},
            {5545.17, 1226047488}, {5565.81, 1226309632}, {5589.38, 1227096064},
            {5610.02, 1227358208}, {5677.81, 1229193216}, {5698.45, 1229455360},
            {5722.03, 1230241792}, {5742.67, 1230503936}, {5810.45, 1500774400},
            {5831.1, 1501036544}, {5855.68, 1501822976}, {5876.33, 1502085120},
            {5900.91, 1502871552}, {5921.56, 1503133696}, {5946.14, 1503920128},
            {5966.79, 1504182272}, {5991.37, 1504968704}, {6012.01, 1505230848},
            {6036.6, 1506017280}, {6057.24, 1506279424}, {6081.83, 1507065856},
            {6102.47, 1507328000}, {6127.06, 1508114432}, {6147.7, 1508376576},
            {6172.29, 1509163008}, {6192.93, 1509425152}, {6217.52, 1510211584},
            {6238.16, 1510473728}, {6262.74, 1511260160}, {6283.39, 1511522304},
            {6307.97, 1512308736}, {6328.62, 1512570880}, {6353.2, 1513357312},
            {6373.85, 1513619456}, {6398.43, 1514405888}, {6419.08, 1514668032},
            {6443.66, 1515454464}, {6464.3, 1515716608}, {6488.89, 1516503040},
            {6509.53, 1516765184}, {6534.12, 1517551616}, {6554.76, 1517813760},
            {6579.35, 1518600192}, {6599.99, 1518862336}, {6624.58, 1519648768},
            {6645.22, 1519910912}, {6669.81, 1520697344}, {6690.45, 1520959488},
            {6715.03, 1521745920}, {6735.68, 1522008064}, {6760.26, 1522794496},
            {6780.91, 1523056640}, {6805.49, 1523843072}, {6826.14, 1524105216},
	    {6871.37, 1525153792}, {6895.95, 1525940224}, {6916.59, 1526202368},
            {6941.18, 1526988800}, {6961.82, 1527250944}, {6986.41, 1528037376}, {7007.05, 1528299520},
            {7031.64, 1529085952}, {7052.28, 1529348096}, {7076.87, 1530134528}, {7097.51, 1530396672},
            {7122.1, 1531183104}, {7142.74, 1531445248}, {7167.32, 1532231680}, {7187.97, 1532493824},
            {7212.55, 1533280256}, {7233.2, 1533542400}, {7257.78, 1534328832}, {7278.43, 1534590976},
            {7303.01, 1535377408}, {7323.66, 1535639552}, {7348.24, 1536425984}, {7368.88, 1536688128},
            {7393.47, 1537474560}, {7414.11, 1537736704}, {7438.7, 1538523136}, {7459.34, 1538785280},
            {7483.93, 1539571712}, {7504.57, 1539833856}, {7529.16, 1540620288}, {7549.8, 1540882432},
            {7574.39, 1541668864}, {7595.03, 1541931008}, {7619.61, 1542717440}, {7640.26, 1542979584},
            {7664.84, 1543766016}, {7685.49, 1544028160}, {7710.07, 1544814592}, {7730.72, 1545076736},
            {7755.3, 1545863168}, {7775.95, 1546125312}, {7800.53, 1546911744}, {7821.17, 1547173888},
            {7845.76, 1547960320}, {7866.4, 1548222464}, {7890.99, 1549008896}, {7911.63, 1549271040},
            {7936.22, 1550057472}, {7956.86, 1550319616}, {7981.45, 1551106048}, {8002.09, 1551368192},
            {8026.68, 1552154624}, {8047.32, 1552416768}, {8071.9, 1553203200}, {8092.55, 1553465344},
            {8117.13, 1554251776}, {8137.78, 1554513920}, {8162.36, 1555300352}, {8183.01, 1555562496},
            {8207.59, 1556348928}, {8228.24, 1556611072}, {8252.82, 1557397504}, {8273.46, 1557659648},
            {8298.05, 1558446080}, {8318.69, 1558708224}, {8343.28, 1559494656}, {8363.92, 1559756800},
            {8388.51, 1560543232}, {8409.15, 1560805376}, {8433.74, 1561591808}, {8454.38, 1561853952},
            {8478.97, 1562640384}, {8499.61, 1562902528}, {8524.19, 1563688960}, {8544.84, 1563951104},
            {8569.42, 1564737536}, {8590.07, 1564999680}, {8614.65, 1700003840}, {8635.3, 1700265984},
            {8659.46, 1701052416}, {8680.1, 1701314560}, {8704.26, 1702100992}, {8724.9, 1702363136},
            {8749.06, 1703149568}, {8769.71, 1703411712}, {8793.86, 1704198144}, {8814.51, 1704460288},
            {8838.67, 1705246720}, {8859.31, 1705508864}, {8883.47, 1706295296}, {8904.11, 1706557440},
            {8928.27, 1707343872}, {8948.92, 1707606016}, {8973.08, 1708392448}, {8993.72, 1708654592},
            {9017.88, 1709441024}, {9038.52, 1709703168}, {9088.08, 1844969472}, {9135.41, 1846018048},
            {9182.75, 1847066624}, {9230.08, 1848115200}, {9277.41, 1849163776}, {9324.74, 1850212352},
            {9372.08, 1851260928}, {9419.41, 1852309504}, {9466.74, 1853358080}, {9514.07, 1854406656},
            {9561.41, 1855455232}, {9608.74, 1856503808}, {9656.07, 1857552384}, {9703.4, 1858600960},
            {9750.74, 1859649536}, {9798.07, 1860698112}, {9845.4, 1861746688}, {9892.73, 1862795264},
            {9940.07, 1863843840}, {9987.4, 1864892416}   
     };

  for (size_t i = 0; i < vtxToPlaneID.size() - 1; ++i) {
      double lowerBound = vtxToPlaneID[i].first;
      double upperBound = vtxToPlaneID[i + 1].first;

      if (vtx >= lowerBound && vtx < upperBound) {
         return vtxToPlaneID[i].second;
      }
  }

  std::cerr << "Error: "<<vtx<<" vtx value out of range!\n";
  return -1; // Indicate an error
}

//=================================================================================
/*
double CVUniverse::calcRecoQ2() const{

   double M_muon = 105.6583; //MeV
   double pmu2_Reco = ( GetVecElem("MasterAnaDev_leptonE",3)*GetVecElem("MasterAnaDev_leptonE",3)) - M_muon*M_muon;
   double pmu_Reco= pmu2_Reco>0 ? sqrt(pmu2_Reco) : 0. ;

   double MasterAnaDev_ANN_recoil_E = GetDouble("MasterAnaDev_recoil_E")/1.0e3;
   double LepReco3 = GetDouble("MasterAnaDev_leptonE",3);

   double Emu_cc_Reco = LepReco3/1.e3; //Reco Muon Energy
   double MasterAnaDev_E_reco = Emu_cc_Reco + MasterAnaDev_ANN_recoil_E;

   double LepReco0 = GetDouble("MasterAnaDev_leptonE",0);//reco_muon_px
   double LepReco1 = GetDouble("MasterAnaDev_leptonE",1);//reco_muon_py
   double LepReco2 = GetDouble("MasterAnaDev_leptonE",2);//reco_muon_pz
   double numi_beam_angle_rad =  -0.05887;
   double pyprime = -1.0*sin(numi_beam_angle_rad)*LepReco2 + cos(numi_beam_angle_rad)*LepReco1;
   double pzprime =  1.0*cos(numi_beam_angle_rad)*LepReco2 + sin(numi_beam_angle_rad)*LepReco1;

   double pSquare = pow(LepReco0,2) + pow(pyprime,2) + pow(pzprime,2);
   double theta_Reco = acos( pzprime / sqrt(pSquare) );
   double q2_reco = ((2*MasterAnaDev_E_reco*1.0e3)*(LepReco3 - (pmu_Reco*cos((theta_Reco))))- M_muon*M_muon)/1.0e6;

   return q2_reco; //GeV^2
}
*/

double CVUniverse::calcEavailTrue() const {

  int mc_nFSPart = GetInt("mc_nFSPart"); //read the number of particles in the final state

  double protonE = 0.0;
  double neutronE = 0.0;
  double pionE = 0.0;
  double pizeroE = 0.0;
  double electronE = 0.0;
  double gammaE = 0.0;
  double otherE = 0.0;
  double missingE = 0.0;

  int protonN = 0;
  int neutronN = 0;
  int pionN = 0;
  int pizeroN = 0;
  int electronN = 0;
  int gammaN = 0;
  int otherN = 0;

  double M_pion = 135;
  double M_p = 938.27;
  double M_n= 940.6;
  double mass_sigma =1189.37; //Sigma zero 1192.642, sigma minus 1197.449
  double mass_kaon = 493.677;
  //double mass_lambda = 938
  double recoil = 0.0;

  for(int k=0; k< mc_nFSPart; k++){
      int pdg = GetVecElemInt("mc_FSPartPDG",k); //read the PDG code for the given particle
      double E = GetVecElem("mc_FSPartE",k); //read the Energy of the given FS particle

      if(pdg==2212){
        protonE += E - M_p;
        if(E - M_p > 5.0)protonN++;
      } else if(pdg==2112){
        neutronE += E - M_n;
        if(E - M_n > 5.0)neutronN++;
      } else if(abs(pdg)==211){
        pionE += E;
        pionN++;
      } else if(pdg==111){
        pizeroE += E;
        pizeroN++;
      } else if(abs(pdg)==11){
        electronE += E;
        electronN++;
      } else if(fabs(pdg)==22){
        gammaE += E;
        gammaN++;
      } else if(abs(pdg) != 13){
        // The Other Category.
        //cout << "other " <<  pdg << " energy " << E << endl;
        if(pdg>=2000000000){
          missingE += E;
        } else if(pdg>1000000000){
          // do nothing, assume negligible energy for nucleons
          // save me the trouble of asking what that nucleon's mass was.
        } else if(pdg>=2000){
          otherE += E - M_p;
          otherN++;
          //primarily strange baryons 3122s
        } else if(pdg<=-2000){
          otherE += E + M_p;
          otherN++;
          //primarily anti-baryons -2112 and -2212s,
        } else {
          otherE += E;
          otherN++;
          //primarily kaons and eta mesons.
        }
      }
  }//end loop over fs particles.

  double visibleE = protonE + pionE - pionN*M_pion + pizeroE + electronE + gammaE + otherE;

  return visibleE = visibleE/1.0e3;

  //return 1.0;
}


double CVUniverse::calcTrueq3() const{

  double LepTruth0 = GetVecElem("mc_primFSLepton",0);
  double LepTruth1 = GetVecElem("mc_primFSLepton",1);
  double LepTruth2 = GetVecElem("mc_primFSLepton",2);

  double inc_Px = GetVecElem("mc_incomingPartVec",0);
  double inc_Py = GetVecElem("mc_incomingPartVec",1);
  double inc_Pz = GetVecElem("mc_incomingPartVec",2);

  double px = LepTruth0 - inc_Px;
  double py = LepTruth1 - inc_Py;
  double pz = LepTruth2 - inc_Pz;

  double q3 = sqrt(px*px+py*py+pz*pz);
  return q3/1000.0;
}

double CVUniverse::calcQ2(const double Enu,const double Emu,const double Thetamu) const {
  double pmu2= Emu*Emu- MinervaUnits::M_mu*MinervaUnits::M_mu;
  double pmu= pmu2>0? sqrt(pmu2) : 0. ;
  double Q2 = 2* Enu * (Emu - pmu * cos(Thetamu) ) - pow(MinervaUnits::M_mu, 2);
    //use unapproximated version to match with branch value  //4*Enu*Emu*pow( sin( Thetamu/2), 2.);
  return Q2;
}


double CVUniverse::calcW(const double Q2,const double Ehad) const {
  double nuclMass = (M_neutron + M_proton)/2. ;//Tmp change for validation with (GetAnaToolName() + //M_nucleon;
  double W2 = ( pow(nuclMass, 2) +  2. * ( Ehad ) * nuclMass - Q2);
  W2 = W2 > 0 ? sqrt(W2) : 0.0;
  return W2;
}

double CVUniverse::calcX(const double Q2,const double Enu, const double Emu) const {
  double nuclMass = M_nucleon;
  //if( NEUTRON_PDG  == GetInt("mc_targetNucleon") )
  //  nuclMass = M_neutron;
  //else if( PROTON_PDG == GetInt("mc_targetNucleon") )
  //  nuclMass = M_proton;
  double x = Q2 / (2. * ( Enu -  Emu) * nuclMass);
  return x;
}

double CVUniverse::calcY(const double Enu,const double Ehad) const {
  return Ehad/Enu;
}


double CVUniverse::calcq3(const double 	Q2,const double Enu,const double Emu	 )const
{
	return sqrt(Q2 + pow(Enu - Emu,2.0));
}

double CVUniverse::calcq0(const double Enu,const double Emu)const
{
	return Enu-Emu;
}

double CVUniverse::GetWeightEmu() const {
   //const bool do_warping = true;
   double wgt_flux=1., wgt_2p2h=1.;
   double wgt_rpa=1.,   wgt_nrp=1.,  wgt_lowq2=1.;
   double wgt_genie=1., wgt_mueff=1.;
   double wgt_anisodd=1.;
   double wgt_emu=1.;

   // genie
   wgt_genie = GetGenieWeight();
   double Enu  = GetDouble("mc_incomingE")*1e-3;
   int nu_type = GetInt("mc_incoming");
   // flux
   wgt_flux = GetFluxAndCVWeight(Enu, nu_type);
   ///For Emu Iron
   //wgt_emu = 1.06991e-07 * GetMuonETrueGeV()*GetMuonETrueGeV()*GetMuonETrueGeV()*GetMuonETrueGeV()*GetMuonETrueGeV() - 1.02357e-05* GetMuonETrueGeV()*GetMuonETrueGeV()*GetMuonETrueGeV()*GetMuonETrueGeV()+ 0.000378062* GetMuonETrueGeV()*GetMuonETrueGeV()*GetMuonETrueGeV() - 0.00713933*GetMuonETrueGeV()*GetMuonETrueGeV() + 0.0670028*GetMuonETrueGeV() + 0.919948;
   //For Iron 2D
    //wgt_emu = 1.34*(5.4546e-05 * GetMuonETrueGeV()* GetMuonETrueGeV() + 9.65024e-06 * GetEhadTrueGeV() * GetEhadTrueGeV() + 0.760717);
    //For Tracker 2D
    wgt_emu = 1.34*(4.21486e-05 * GetMuonETrueGeV()* GetMuonETrueGeV() + 7.51935e-05 * GetEhadTrueGeV() * GetEhadTrueGeV() + 0.734119);
   //For Emu Lead
   //wgt_emu = 2.22052e-07 * GetMuonETrueGeV()*GetMuonETrueGeV()*GetMuonETrueGeV()*GetMuonETrueGeV()*GetMuonETrueGeV() - 2.39719e-05* GetMuonETrueGeV()*GetMuonETrueGeV()*GetMuonETrueGeV()*GetMuonETrueGeV()+ 0.000979981* GetMuonETrueGeV()*GetMuonETrueGeV()*GetMuonETrueGeV() - 0.0188648*GetMuonETrueGeV()*GetMuonETrueGeV() + 0.163442*GetMuonETrueGeV() + 0.668094;
   //For Emu Tracker
   //wgt_emu = 1.06991e-07 * GetMuonETrueGeV()*GetMuonETrueGeV()*GetMuonETrueGeV()*GetMuonETrueGeV()*GetMuonETrueGeV() - 1.02357e-05* GetMuonETrueGeV()*GetMuonETrueGeV()*GetMuonETrueGeV()*GetMuonETrueGeV()+ 0.000378062* GetMuonETrueGeV()*GetMuonETrueGeV()*GetMuonETrueGeV() - 0.00713933*GetMuonETrueGeV()*GetMuonETrueGeV() + 0.0670028*GetMuonETrueGeV() + 0.919948;
   // 2p2h
   wgt_2p2h = GetLowRecoil2p2hWeight();
   // rpa
   wgt_rpa = GetRPAWeight();
   //MINOS Efficiency
   wgt_mueff = GetMinosEfficiencyWeight();

   return wgt_flux*wgt_genie*wgt_rpa*wgt_nrp*wgt_lowq2*wgt_mueff*wgt_2p2h*wgt_emu;
   //cout<< "\t NonResWeight: "<<GetGenieWeight()<<"\t 2p2hWeight: "<<GetLowRecoil2p2hWeight()<<"\t RPAWeight: "<<GetRPAWeight()<<"\t MINOS Eff: "<<GetMinosEfficiencyWeight()<<"\t\t  wgt_emu"<<wgt_emu<<endl;

}
double CVUniverse::GetWeightQ2() const {
   //const bool do_warping = true;
   double wgt_flux=1., wgt_2p2h=1.;
   double wgt_rpa=1.,   wgt_nrp=1.,  wgt_lowq2=1.;
   double wgt_genie=1., wgt_mueff=1.;
   double wgt_q2=1.;

   // genie
   wgt_genie = GetGenieWeight();
   double Enu  = GetDouble("mc_incomingE")*1e-3;
   int nu_type = GetInt("mc_incoming");
   // flux
   wgt_flux = GetFluxAndCVWeight(Enu, nu_type);
   //For Iron
   wgt_q2 = 0.00143089*GetQ2TrueGeV()+1.16667;
   //For Lead
   //wgt_q2 = 0.0013294*GetQ2TrueGeV()+1.11518;
   //For Tracker
   //wgt_q2 = -0.000423594*GetQ2TrueGeV()+1.12566;
   // 2p2h
   wgt_2p2h = GetLowRecoil2p2hWeight();
   // rpa
   wgt_rpa = GetRPAWeight();
   //MINOS Efficiency
   wgt_mueff = GetMinosEfficiencyWeight();

   return wgt_flux*wgt_genie*wgt_rpa*wgt_nrp*wgt_q2*wgt_mueff*wgt_2p2h;

}

double CVUniverse::GetWeighty() const {
   //const bool do_warping = true;
   double wgt_flux=1., wgt_2p2h=1.;
   double wgt_rpa=1.,   wgt_nrp=1.,  wgt_lowq2=1.;
   double wgt_genie=1., wgt_mueff=1.;
   double wgt_anisodd=1.;
   double wgt_y=1.;

   // genie
   wgt_genie = GetGenieWeight();
   double Enu  = GetDouble("mc_incomingE")*1e-3;
   int nu_type = GetInt("mc_incoming");
   // flux
   wgt_flux = GetFluxAndCVWeight(Enu, nu_type);
   //For 2D Iron with xy weight function
   //wgt_y = 1.3*(0.00117*GetxTrue()*GetxTrue() + 0.00724* GetyTrue()*GetyTrue() + 0.00146);
   //For 2D Tracker with xy weight function
   wgt_y = 1.3*(-9.76e-05*GetxTrue()*GetxTrue() + 0.0180* GetyTrue()*GetyTrue() + 0.00362);
   //For Iron
   //wgt_y = 1.24836*GetyTrue()*GetyTrue()* GetyTrue() - 1.64309* GetyTrue()*GetyTrue() + 0.49131*GetyTrue() + 1.16857;
   //For Lead
   //wgt_y = 1.95476*GetyTrue()*GetyTrue()* GetyTrue() - 2.80867* GetyTrue()*GetyTrue() + 1.19048*GetyTrue() + 0.979472;
   //For Lead
   //wgt_y = 3.02641*GetyTrue()*GetyTrue()* GetyTrue() - 4.12468* GetyTrue()*GetyTrue() + 1.5317*GetyTrue() + 0.988857;
   // 2p2h
   wgt_2p2h = GetLowRecoil2p2hWeight();
   // rpa
   wgt_rpa = GetRPAWeight();
   //MINOS Efficiency
   wgt_mueff = GetMinosEfficiencyWeight();

   return wgt_flux*wgt_genie*wgt_rpa*wgt_nrp*wgt_lowq2*wgt_mueff*wgt_2p2h*wgt_y;

}
double CVUniverse::GetWeight() const {
   //const bool do_warping = true;
   double wgt_flux=1., wgt_2p2h=1.;
   double wgt_rpa=1.,   wgt_nrp=1.,  wgt_lowq2=1.;
   double wgt_genie=1., wgt_mueff=1.;
   double wgt_anisodd=1.;

   // genie
   wgt_genie = GetGenieWeight();
   double Enu  = GetDouble("mc_incomingE")*1e-3;
   int nu_type = GetInt("mc_incoming");
   // flux
   wgt_flux = GetFluxAndCVWeight(Enu, nu_type);
   // 2p2h
   wgt_2p2h = GetLowRecoil2p2hWeight();
   // rpa
   wgt_rpa = GetRPAWeight();
   //MINOS Efficiency
   wgt_mueff = GetMinosEfficiencyWeight();
   //LowQ2 weight
   //wgt_lowq2 *= GetLowQ2PiWeight("JOINT");
   //cout<< "\t NonResWeight: "<<GetGenieWeight()<<"\t 2p2hWeight: "<<GetLowRecoil2p2hWeight()<<"\t RPAWeight: "<<GetRPAWeight()<<"\t MINOS Eff: "<<GetMinosEfficiencyWeight()<<"\t\t Q0: "<<Getq0True()/1000.0<<"\t Q3: "<<Getq3True()/1000.0<<"LowQ2Weight"<< GetLowQ2PiWeight("JOINT")<<endl;


   return wgt_flux*wgt_genie*wgt_rpa*wgt_nrp*wgt_lowq2*wgt_mueff*wgt_2p2h;

}

double CVUniverse::GetTruthWeight()const{

   double wgt_flux=1., wgt_2p2h=1.;
   double wgt_rpa=1.,   wgt_nrp=1.,  wgt_lowq2=1.;
   double wgt_genie=1., wgt_mueff=1.;
   double wgt_anisodd=1.;

   //There need to be flags added to this to turn on and off different tunes. Same goes for GetWeight().  -- ANF 2020-3-18

  double Enu  = GetDouble("mc_incomingE")*1e-3;
   int nu_type = GetInt("mc_incoming");

    // genie
   wgt_genie = GetGenieWeight();
   // flux
   wgt_flux = GetFluxAndCVWeight(Enu,nu_type);
   // 2p2h
   wgt_2p2h = GetLowRecoil2p2hWeight();
   // rpa
   wgt_rpa = GetRPAWeight();
   //LowQ2 weight
   wgt_lowq2 *= GetLowQ2PiWeight("JOINT");
   return wgt_genie*wgt_flux*wgt_rpa*wgt_nrp*wgt_lowq2*wgt_mueff*wgt_2p2h;
}

//From Faiza Background package//

double CVUniverse::Var( const std::string& name, bool useTrue )
{
   /* if( name == "Emu" )
    {
        if( useTrue )      return GetElepTrue() * mev_to_gev;
        else               return GetEmu() * mev_to_gev ;
    }
    if ( name == "Enu" || name == "E" )
    {
        if( useTrue )      return GetEnuTrue() * mev_to_gev;
        else               return GetEnu() * mev_to_gev;
    }
    if ( name == "Ehad" )
    {
        if( useTrue )      return ( GetEhadTrue() ) * mev_to_gev;
        else               return GetRecoilEnergy()* mev_to_gev;
    }
    if ( name == "CCQE-Recoil")
    {
        //not sure what true would mean for this
        if( useTrue)   return ( GetEhadTrue() ) * mev_to_gev;
        else           return GetCCQERecoil() * mev_to_gev;
    }
    if( name == "VtxEnergy")
    {
        //true means nothing here
        if( useTrue ) return 0.;
        else          return GetVtxEnergy() * mev_to_gev;
    }
    if( name == "ETheta" )
    {
     if( useTrue ) return ( (GetElepTrue() * mev_to_gev ) * ( 1. - cos(GetThetamuTrue()) ));
        else  return ((GetEmu()  * mev_to_gev ) * ( 1. - cos( GetThetamu())) );
    }

    if( name == "x" ){

            if( useTrue ) return GetxTrue();
            else return GetxReco();
    }

    if( name == "y" )
    {
        if( useTrue )    return GetyTrue();
        else             return GetyReco();
    }
    if( name == "Q2" || name == "q2" )
    {
        if( useTrue ){

            return (GetQ2IncTrue());
        }

        else             return ( GetQ2Reco()  );
    }


    if( name == "W" )
    {
        if( useTrue ) return GetWTrue();
        else return GetWReco() ;
    }


    if( name == "PhiMu" )
    {
        if( useTrue )     return GetDouble("truth_muon_phi")*rad_to_deg;
        else              return       GetDouble("muon_phi")*rad_to_deg;
    }

    if( name == "ThetaMu" )
    {
        if( useTrue )      return GetThetamuTrue()*rad_to_deg;
        else               return     GetThetamu()*rad_to_deg;
    }

    if( name == "muonPt" )
    {
        if( useTrue )      return GetlepPtTrue()*rad_to_deg;
                       return     GetMuonPt()*rad_to_deg;
    }

    if( name == "CosThetaMu" )
    {
        if( useTrue )      return cos(GetThetamuTrue());
        else               return cos(GetThetamu());
    }


    if( name == "moduleDNN" )
    {
        if( useTrue ) return GetInt("truth_vtx_module");
        else          return GetInt("(GetAnaToolName() +_vtx_module");
    }
    */
    if( name == "planeDNN" )
    {
        if( useTrue) return GetInt("truth_vtx_module")*2+GetInt("truth_vtx_plane")+10;
        else{
           //return ANN_vtx_modules[0]*2+ANN_vtx_planes[0]+10;
           int ANN_vtx_module, ANN_vtx_plane;
           int targetID = GetTargetFromSegment( GetVecElem("ANN_segments",0), ANN_vtx_module, ANN_vtx_plane );
    return ANN_vtx_module*2+ANN_vtx_plane+10;
        }
    }


    cout << "Error [CVUniverse::Var] : I do not know how to interpret that variable \"" << name << "\", so I am throwing an exception." << endl;
    throw 1;

}


int CVUniverse::GetAtomicNumber(){
   return GetInt("mc_targetA");
}

//========================================================================================================
//    Adding Bunch of Mapping Function To Split The Nuclear Target Region Based on The Target Number
//========================================================================================================
int CVUniverse::AtomicNumberToMass(int targetZ ) const{
    if( targetZ == 6 ) return 12;
    if( targetZ == 26 ) return 56;
    if( targetZ == 207 ) return 207;
    return -999;
}

int CVUniverse::GetTargetMinBin( int targetID ) const{
    if( targetID == 1 ) return 1;
    else if( targetID == 2 ) return 11;
    else if( targetID == 3 ) return 21;
    else if( targetID == 4 ) return 41;
    else if( targetID == 5 ) return 51;
    //else if( targetID == 6 ) return 57;  //faiza
    else{
        std::cerr << "Bad targetID choice in CVUniverse::GetTargetMinBin! Try again. Hint: we only have 5 passive target" << std::endl;
    //}
        return -999;
   }
}

int CVUniverse::GetTargetMaxBin( int targetID ) const{
    if( targetID == 1 ) return 19;
    else if( targetID == 2 ) return 29;
    else if( targetID == 3 ) return 41;
    else if( targetID == 4 ) return 55;
    else if( targetID == 5 ) return 65;
    //else if( targetID == 0 ) return 9;   //faiza
    else{
        std::cerr << "Bad targetID choice in CVUniverse::GetTargetMaxBin! Try again. Hint: we only have 5 passive target" << std::endl;
   // }
        return -999;
    }
}

int CVUniverse::GetTargetPlane( int targetID ) const{
    if( targetID == 1 ) return 9;
    else if( targetID == 2 ) return 19;
    else if( targetID == 3 ) return 30;
    else if( targetID == 4 ) return 49;
    else if( targetID == 5 ) return 55;
    else{
        std::cerr << "Bad targetID choice in CVUniverse::GetTargetMaxBin! Try again. Hint: we only have 5 passive target" << std::endl;
     }
    return -999;
    //}
}

int CVUniverse::GetPlaneTargetID( int plane ) const{
    if( plane == 9 ) return 1;
    else if( plane == 19 ) return 2;
    else if( plane == 30 ) return 3;
    else if( plane == 49 ) return 4;
    else if( plane == 55 ) return 5;
    //else{
    //std::cerr << "Bad targetID choice in CVUniverse::GetPlaneTargetID! Try again." << std::endl;
    //}
    return -999;
}

int CVUniverse::GetTargetUSPlane( int targetID ) const{
    if( targetID == 1 ) return 8;
    else if( targetID == 2 ) return 18;
    else if( targetID == 3 ) return 28;
    else if( targetID == 4 ) return 48;
    else if( targetID == 5 ) return 54;
    else if( targetID == 6 ) return 64;
    //  else{
    //    std::cerr << "Bad targetID choice in CVUniverse::GetTargetUSPlane! Try again. Hint: we only have 5 passive target" << std::endl;
    //  }
    return -999;
}

int CVUniverse::GetTargetDSPlane( int targetID ) const{
    if( targetID == 1 ) return 11;
    else if( targetID == 2 ) return 21;
    else if( targetID == 3 ) return 33;
    else if( targetID == 4 ) return 51;
    else if( targetID == 5 ) return 57;
    else if( targetID == 0 ) return 1;
    //  else{
    //    std::cerr << "Bad targetID choice in CVUniverse::GetTargetUSPlane! Try again. Hint: we only have 5 passive target" << std::endl;
    //  }
    return -999;
}

double CVUniverse::GetTargetZStart( int targetID ) const
{
    if( targetID == 1 ) return 4469.46;
    else if( targetID == 2 ) return 4689.73;
    else if( targetID == 3 ) return 4907.82;

    else if( targetID == 5 ) return 5769.88;
    //  else{
    //    std::cerr << "Bad targetID choice in CVUniverse::GetTargetUSPlane! Try again. Hint: we only have 5 passive target" << std::endl;
    //  }
    return -999.9;
}

double CVUniverse::GetTargetZEnd( int targetID ) const
{
    if( targetID == 1 ) return 4495.2;
    else if( targetID == 2 ) return 4715.47;
    else if( targetID == 3 ) return 4984.96;
    else if( targetID == 4 ) return 5645.73;
    else if( targetID == 5 ) return 5782.87;
    //  else{
    //    std::cerr << "Bad targetID choice in CVUniverse::GetTargetUSPlane! Try again. Hint: we only have 5 passive target" << std::endl;
    //  }
    return -999.9;
}



int CVUniverse::GetTargetFromSegment( int segment, int& vtx_module, int& vtx_plane ) const
{
    //int targetID = -999;
    int targetID = -1;
    vtx_module = -999;
    vtx_plane = -999;

    if( segment == 0 ){
        vtx_module = -5;
        vtx_plane = 0;
        targetID = -1;
    }
    else if( segment == 1 ){
        vtx_module = -5;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 2 ){
        vtx_module = -5;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 3 ){
        vtx_module = -4;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 4 ){
        vtx_module = -4;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 5 ){
        vtx_module = -3;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 6 ){
        vtx_module = -3;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 7 ){
        vtx_module = -2;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 8 ){
        vtx_module = -2;
        vtx_plane = 2;
        targetID = -1;
        //targetID = 1;
    }
    else if( segment == 9 ){
        vtx_module = -1;
        vtx_plane = 1;
        //targetID = -1;
        targetID = 1;
    }
    else if( segment == 10 ){
        vtx_module = 0;
        vtx_plane = 1;
        targetID = -1;
        //targetID = 1;
    }
    else if( segment == 11 ){
        vtx_module = 0;
        vtx_plane = 2;
        targetID = -1;
        //targetID = 1;
    }
    else if( segment == 12 ){
        vtx_module = 1;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 13 ){
        vtx_module = 1;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 14 ){
        vtx_module = 2;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 15 ){
        vtx_module = 2;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 16 ){
        vtx_module = 3;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 17 ){
        vtx_module = 3;
        vtx_plane = 2;
        targetID = -1;
        //targetID = 2;
    }
    else if( segment == 18 ){
        vtx_module = 4;
        vtx_plane = 1;
        //targetID = -1;
        targetID = 2;
    }
    else if( segment == 19 ){
        vtx_module = 5;
        vtx_plane = 1;
        targetID = -1;
        //targetID = 2;
    }
    else if( segment == 20 ){
        vtx_module = 5;
        vtx_plane = 2;
        targetID = -1;
        //targetID = 2;
    }
    else if( segment == 21 ){
        vtx_module = 6;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 22 ){
        vtx_module = 6;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 23 ){
        vtx_module = 7;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 24 ){
        vtx_module = 7;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 25 ){
        vtx_module = 8;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 26 ){
        vtx_module = 8;
        vtx_plane = 2;
        targetID = -1;
        //targetID = 3;
    }
    else if( segment == 27 ){
        vtx_module = 9;
        vtx_plane = 2;
        //targetID = -1;
        targetID = 3;
    }
    else if( segment == 28 ){
        vtx_module = 11;
        vtx_plane = 1;
        targetID = -1;
        //targetID = 3;
    }
    else if( segment == 29 ){
        vtx_module = 11;
        vtx_plane = 2;
        targetID = -1;
        //targetID = 3;
    }
    else if( segment == 30 ){
        vtx_module = 12;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 31 ){
        vtx_module = 12;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 32 ){
        vtx_module = 13;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 33 ){
        vtx_module = 13;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 34 ){
        vtx_module = 14;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 35 ){
        vtx_module = 14;
        vtx_plane = 2;
        targetID = -1;
    }
   //segment 36 is water!!
    else if( segment == 36 ){
        vtx_module = 15;
        vtx_plane = 1;
        targetID = 6;
    }
    /*else if( segment == 36 ){
        vtx_module = -999;
        vtx_plane = -999;
        targetID = 6;
    }
    else if( segment == 37 ){
        vtx_module = 15;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 38 ){
        vtx_module = 15;
        vtx_plane = 2;
        targetID = -1;
    }*/
    else if( segment == 39 ){
        vtx_module = 16;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 40 ){
        vtx_module = 16;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 41 ){
        vtx_module = 17;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 42 ){
        vtx_module = 17;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 43 ){
        vtx_module = 18;
        vtx_plane = 1;
        targetID = -1;
        //targetID = 4;
    }
    else if( segment == 44 ){
        vtx_module = 18;
        vtx_plane = 2;
        targetID = -1;
        //targetID = 4;
    }
    else if( segment == 45 ){
        vtx_module = 19;
        vtx_plane = 1;
        //targetID = -1;
        targetID = 4;
    }
    else if( segment == 46 ){
        vtx_module = 20;
        vtx_plane = 1;
        targetID = -1;
        //targetID = 4;
    }
    else if( segment == 47 ){
        vtx_module = 20;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 48 ){
        vtx_module = 21;
        vtx_plane = 1;
        targetID = -1;
        //targetID = 5;
    }
    else if( segment == 49 ){
        vtx_module = 21;
        vtx_plane = 2;
        targetID = -1;
        //targetID = 5;
    }
    else if( segment == 50 ){
        vtx_module = 22;
        vtx_plane = 1;
        //targetID = -1;
        targetID = 5;
    }
    else if( segment == 51 ){
        vtx_module = 23;
        vtx_plane = 1;
        targetID = -1;
        //targetID = 5;
    }
    else if( segment == 52 ){
        vtx_module = 23;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 53 ){
        vtx_module = 24;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 54 ){
        vtx_module = 24;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 55 ){
        vtx_module = 25;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 56 ){
        vtx_module = 25;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 57 ){
        vtx_module = 26;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 58 ){
        vtx_module = 26;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 59 ){
        vtx_module = 27;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 60 ){
        vtx_module = 27;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 61 ){
        vtx_module = 28;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 62 ){
        vtx_module = 28;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 63 ){
        vtx_module = 29;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 64 ){
        vtx_module = 29;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 65 ){
        vtx_module = 30;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 66 ){
        vtx_module = 30;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 67 ){
        vtx_module = 31;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 68 ){
        vtx_module = 31;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 69 ){
        vtx_module = 32;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 70 ){
        vtx_module = 32;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 71 ){
        vtx_module = 33;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 72 ){
        vtx_module = 33;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 73 ){
        vtx_module = 34;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 74 ){
        vtx_module = 34;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 75 ){
        vtx_module = 35;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 76 ){
        vtx_module = 35;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 77 ){
        vtx_module = 36;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 78 ){
        vtx_module = 36;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 79 ){
        vtx_module = 37;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 80 ){
        vtx_module = 37;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 81 ){
        vtx_module = 38;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 82 ){
        vtx_module = 38;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 83 ){
        vtx_module = 39;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 84 ){
        vtx_module = 39;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 85 ){
        vtx_module = 40;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 86 ){
        vtx_module = 40;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 87 ){
        vtx_module = 41;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 88 ){
        vtx_module = 41;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 89 ){
        vtx_module = 42;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 90 ){
        vtx_module = 42;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 91 ){
        vtx_module = 43;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 92 ){
        vtx_module = 43;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 93 ){
        vtx_module = 44;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 94 ){
        vtx_module = 44;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 95 ){
        vtx_module = 45;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 96 ){
        vtx_module = 45;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 97 ){
        vtx_module = 46;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 98 ){
        vtx_module = 46;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 99 ){
        vtx_module = 47;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 100 ){
        vtx_module = 47;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 101 ){
        vtx_module = 48;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 102 ){
        vtx_module = 48;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 103 ){
        vtx_module = 49;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 104 ){
        vtx_module = 49;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 105 ){
        vtx_module = 50;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 106 ){
        vtx_module = 50;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 107 ){
        vtx_module = 51;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 108 ){
        vtx_module = 51;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 109 ){
        vtx_module = 52;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 110 ){
        vtx_module = 52;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 111 ){
        vtx_module = 53;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 112 ){
        vtx_module = 53;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 113 ){
        vtx_module = 54;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 114 ){
        vtx_module = 54;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 115 ){
        vtx_module = 55;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 116 ){
        vtx_module = 55;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 117 ){
        vtx_module = 56;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 118 ){
        vtx_module = 56;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 119 ){
        vtx_module = 57;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 120 ){
        vtx_module = 57;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 121 ){
        vtx_module = 58;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 122 ){
        vtx_module = 58;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 123 ){
        vtx_module = 59;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 124 ){
        vtx_module = 59;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 125 ){
        vtx_module = 60;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 126 ){
        vtx_module = 60;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 127 ){
        vtx_module = 61;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 128 ){
        vtx_module = 61;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 129 ){
        vtx_module = 62;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 130 ){
        vtx_module = 62;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 131 ){
        vtx_module = 63;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 132 ){
        vtx_module = 63;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 133 ){
        vtx_module = 64;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 134 ){
        vtx_module = 64;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 135 ){
        vtx_module = 65;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 136 ){
        vtx_module = 65;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 137 ){
        vtx_module = 66;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 138 ){
        vtx_module = 66;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 139 ){
        vtx_module = 67;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 140 ){
        vtx_module = 67;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 141 ){
        vtx_module = 68;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 142 ){
        vtx_module = 68;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 143 ){
        vtx_module = 69;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 144 ){
        vtx_module = 69;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 145 ){
        vtx_module = 70;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 146 ){
        vtx_module = 70;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 147 ){
        vtx_module = 71;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 148 ){
        vtx_module = 71;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 149 ){
        vtx_module = 72;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 150 ){
        vtx_module = 72;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 151 ){
        vtx_module = 73;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 152 ){
        vtx_module = 73;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 153 ){
        vtx_module = 74;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 154 ){
        vtx_module = 74;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 155 ){
        vtx_module = 75;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 156 ){
        vtx_module = 75;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 157 ){
        vtx_module = 76;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 158 ){
        vtx_module = 76;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 159 ){
        vtx_module = 77;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 160 ){
        vtx_module = 77;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 161 ){
        vtx_module = 78;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 162 ){
        vtx_module = 78;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 163 ){
        vtx_module = 79;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 164 ){
        vtx_module = 79;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 165 ){
        vtx_module = 80;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 166 ){
        vtx_module = 80;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 167 ){
        vtx_module = 81;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 168 ){
        vtx_module = 81;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 169 ){
        vtx_module = 82;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 170 ){
        vtx_module = 82;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 171 ){
        vtx_module = 83;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 172 ){
        vtx_module = 83;
        vtx_plane = 2;
        targetID = -1;
    }
    else if( segment == 173 ){
        vtx_module = 84;
        vtx_plane = 1;
        targetID = -1;
    }
    else if( segment == 174 ){
        vtx_module = 84;
        vtx_plane = 2;
        targetID = -1;
    }
  return targetID;

}

// ARACHNE EVENT DISPLAY
// data
int CVUniverse::GetRunN() const{return GetInt("ev_run");}
int CVUniverse::GetSubRunN() const{return GetInt("ev_subrun");}
int CVUniverse::GetGateN() const{return GetInt("ev_gate");}
int CVUniverse::GetSliceN() const{return GetVecElem("slice_numbers",0);}
int CVUniverse::Getvtx0N() const{return GetVecElem("vtx",0);}
int CVUniverse::Getvtx1N() const{return GetVecElem("vtx",1);}
int CVUniverse::Getvtx2N() const{return GetVecElem("vtx",2);}
int CVUniverse::Getvtx3N() const{return GetVecElem("vtx",3);}

// MC
int CVUniverse::GetMCRunN() const{return GetInt("mc_run");}
int CVUniverse::GetMCSubRunN() const{return GetInt("mc_subrun");}
int CVUniverse::GetMCGateN() const{return GetInt("mc_nthEvtInFile");}
int CVUniverse::GetMCSliceN() const{return GetVecElem("slice_numbers",0);}

//=====================
// Systematic shift gets
//=====================
//double CVUniverse::GetCCQERecoil( unsigned int shift )
double CVUniverse::GetCCQERecoil( unsigned int shift )
{
    // if( 0 == shift ) return blob_iso_E_nucl + blob_iso_E_tracker + blob_iso_E_ecal + blob_disp_E_nucl + blob_disp_E_tracker + blob_disp_E_ecal;
    // if( 1 == shift ) return blob_iso_alt01_E_nucl + blob_iso_alt01_E_tracker + blob_iso_alt01_E_ecal + blob_disp_alt01_E_nucl + blob_disp_alt01_E_tracker + blob_disp_alt01_E_ecal;
    // if( 2 == shift ) return blob_iso_alt02_E_nucl + blob_iso_alt02_E_tracker + blob_iso_alt02_E_ecal + blob_disp_alt02_E_nucl + blob_disp_alt02_E_tracker + blob_disp_alt02_E_ecal;
    // if( 3 == shift ) return blob_iso_alt03_E_nucl + blob_iso_alt03_E_tracker + blob_iso_alt03_E_ecal + blob_disp_alt03_E_nucl + blob_disp_alt03_E_tracker + blob_disp_alt03_E_ecal;
    // if( 4 == shift ) return blob_iso_alt04_E_nucl + blob_iso_alt04_E_tracker + blob_iso_alt04_E_ecal + blob_disp_alt04_E_nucl + blob_disp_alt04_E_tracker + blob_disp_alt04_E_ecal;

    Error("CVUniverse::GetCCQERecoil", "Valid values for the shift are 0-4");
    throw 1;
}

double CVUniverse::GetVtxEnergy( unsigned int shift )
{
    if( 0 == shift ) return GetDouble("blob_vtx_E");
    // if( 1 == shift ) return blob_vtx_alt01_E;
    // if( 2 == shift ) return blob_vtx_alt02_E;
    // if( 3 == shift ) return blob_vtx_alt03_E;
    // if( 4 == shift ) return blob_vtx_alt04_E;

    Error("CVUniverse::GetVtxEnergy", "Valid values for the shift are 0-4");
    throw 1;
}


#endif
