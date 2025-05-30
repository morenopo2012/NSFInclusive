#ifndef CVUNIVERSE_cxx
#define CVUNIVERSE_cxx 1

#include "../include/CVUniverse.h" 
//#include "../include/NukeCC_Cuts.h" 
#include "PlotUtils/FluxReweighter.h"
#include "PlotUtils/HyperDimLinearizer.h"
#include "PlotUtils/MinosMuonEfficiencyCorrection.h"
#include "PlotUtils/MinosMuonPlusEfficiencyCorrection.h"
#include "PlotUtils/PlotUtilsPhysicalConstants.h"

#include "PlotUtils/DefaultCVUniverse.h"
#include "PlotUtils/ChainWrapper.h"
//using namespace globalV;
using namespace NUKECC_ANA;
//again the constructor.....
CVUniverse::CVUniverse(PlotUtils::ChainWrapper *chw,double nsigma):PlotUtils::DefaultCVUniverse(chw,nsigma)
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



//calling branches and variables
//
double CVUniverse::GetRecoilEnergy()  const { return GetDouble("NukeCC_recoil_E"); }
double CVUniverse::GetThetamu()       const {return GetDouble("muon_theta");}
double CVUniverse::GetEnu()           const {return GetEmu()+ GetRecoilEnergy();} 
double CVUniverse::GetQ2Reco()        const { return calcRecoQ2(GetEnu(),  GetEmu(), GetThetamu()); }
double CVUniverse::GetWReco()         const{return calcWReco(GetQ2Reco(), GetRecoilEnergy());}
double CVUniverse::GetxReco()         const{return calcXReco(GetQ2Reco(), GetEnu(),GetEmu());}
double CVUniverse::GetyReco()         const{return calcYReco(GetEnu(), GetRecoilEnergy());}

// for true variables
double CVUniverse::GetEhadTrue()    const{return GetEnuTrue() - GetElepTrue();}
double CVUniverse::GetThetamuTrue( )const{return GetDouble("truth_muon_theta");}

double CVUniverse::GetQ2IncTrue()      const{return calcTrueQ2(GetEnuTrue(),  GetElepTrue(), GetThetamuTrue()); }
double CVUniverse::GetWTrue()       const{return calcWTrue(GetQ2IncTrue(), GetEhadTrue());}

double CVUniverse::GetxTrue()       const{return calcXTrue(GetQ2IncTrue(), GetEnuTrue(),GetElepTrue());}
double CVUniverse::GetyTrue()       const{return calcYTrue(GetEnuTrue(), GetEhadTrue());}
double CVUniverse::Getq3True()      const{return calcq3(GetQ2True(),  GetEnuTrue(), GetElepTrue());}
double CVUniverse::Getq0True()      const{return calcq0(GetEnuTrue(), GetElepTrue());}
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

double CVUniverse::GetQ2RecoGeV()        const { return (calcRecoQ2(GetEnu(),  GetEmu(), GetThetamu()))*mev_to_gev*mev_to_gev; }
double CVUniverse::GetWRecoGeV()         const{return (calcWReco(GetQ2Reco(), GetRecoilEnergy()))*mev_to_gev;}

double CVUniverse::GetQ2TrueGeV()      const{return (calcTrueQ2(GetEnuTrue(),  GetElepTrue(), GetThetamuTrue()))*mev_to_gev*mev_to_gev; }

double CVUniverse::GetWTrueGeV()       const{return (calcWTrue(GetQ2IncTrue(), GetEhadTrue()))*mev_to_gev;}

double CVUniverse::GetMuonPt()const{
   return GetMuonP()*sin(GetThetamu());
}

double CVUniverse::GetMuonPz()const{
  return GetMuon4V().Z();
}

double CVUniverse::GetlepPtTrue()const{
   return GetPlepTrue()*sin(GetThetalepTrue());

}
double CVUniverse::calcTrueQ2(const double EnuTrue,const double EmuTrue,const double ThetamuTrue)	
const {
double Q2 = 4*EnuTrue*EmuTrue*pow( sin( ThetamuTrue/2), 2.);
return Q2; 
}
  
double CVUniverse::calcWTrue(const double Q2True,const double EhadTrue)
const {
double nuclMass = M_nucleon;
 //if( NEUTRON_PDG  == GetInt("mc_targetNucleon") )
                //nuclMass = M_neutron;
            //else if( PROTON_PDG == GetInt("mc_targetNucleon") )
                nuclMass = M_proton;

double W2 = ( pow(nuclMass, 2) +  2. * ( EhadTrue ) * nuclMass - Q2True);
W2 = W2 > 0 ? sqrt(W2) : 0.0;
return W2;
}  

double CVUniverse::calcRecoQ2(const double Enu,const double Emu,const double Thetamu)	
const {
double Q2Reco = 4*Enu*Emu*pow( sin( Thetamu/2), 2.);
return Q2Reco; 
}
  
double CVUniverse::calcWReco(const double Q2,const double Ehad)
const {
double nuclMass = M_nucleon;
 //if( NEUTRON_PDG  == GetInt("mc_targetNucleon") )
                //nuclMass = M_neutron;
           // else if( PROTON_PDG == GetInt("mc_targetNucleon") )
                nuclMass = M_proton;

double W2Reco = ( pow(nuclMass, 2) +  2. * ( Ehad ) * nuclMass - Q2);
W2Reco = W2Reco > 0 ? sqrt(W2Reco) : 0.0;
return W2Reco;
}  


double CVUniverse::calcXTrue(const double Q2True,const double EnuTrue, const double EmuTrue)
const {
double nuclMass = M_nucleon;
 //if( NEUTRON_PDG  == GetInt("mc_targetNucleon") )
                //nuclMass = M_neutron;
            //else if( PROTON_PDG == GetInt("mc_targetNucleon") )
                nuclMass = M_proton;
double x = Q2True / (2. * ( EnuTrue -  EmuTrue) * nuclMass);
return x;
}  

double CVUniverse::calcXReco(const double Q2,const double Enu, const double Emu)
const {
double nuclMass = M_nucleon;
 //if( NEUTRON_PDG  == GetInt("mc_targetNucleon") )
                //nuclMass = M_neutron;
            //else if( PROTON_PDG == GetInt("mc_targetNucleon") )
                nuclMass = M_proton;
double x = Q2 / (2. * ( Enu -  Emu) * nuclMass);
return x;
}  


double CVUniverse::calcYTrue(const double EnuTrue,const double EhadTrue)	
const {
return EhadTrue/EnuTrue; 
}
	
double CVUniverse::calcYReco(const double Enu,const double Ehad)	
const {
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
   // wgt_emu = 1.34*(4.21486e-05 * GetMuonETrueGeV()* GetMuonETrueGeV() + 7.51935e-05 * GetEhadTrueGeV() * GetEhadTrueGeV() + 0.734119);
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
   //For Iron WQ2 weight
   // wgt_emu = (0.0001 * GetWTrueGeV()* GetWTrueGeV() - 3.22e-06 * GetQ2TrueGeV() * GetQ2TrueGeV() + 0.01047);
   //For Tracker WQ2 weight
    wgt_emu = (0.00013 * GetWTrueGeV()* GetWTrueGeV() - 8.97e-06 * GetQ2TrueGeV() * GetQ2TrueGeV() + 0.0247);

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
   double wgt_mk=1.;
   double wgt_dis=1.;
 
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
   //MK Weight
   //wgt_mk = GetMKWeight(); 
   //DIS Weight
   //wgt_dis = GetDISWeight(weightDIS::kAMU); 
 

   //cout<< "\t NonResWeight: "<<GetGenieWeight()<<"\t 2p2hWeight: "<<GetLowRecoil2p2hWeight()<<"\t RPAWeight: "<<GetRPAWeight()<<"\t MINOS Eff: "<<GetMinosEfficiencyWeight()<<"\t\t Q0: "<<Getq0True()/1000.0<<"\t Q3: "<<Getq3True()/1000.0<<endl; 
   

   return wgt_flux*wgt_genie*wgt_rpa*wgt_nrp*wgt_lowq2*wgt_mueff*wgt_2p2h*wgt_mk*wgt_dis;

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
   return wgt_genie*wgt_flux*wgt_rpa*wgt_nrp*wgt_lowq2*wgt_mueff*wgt_2p2h;
}



int CVUniverse::GetTargetFromSegment( int segment, int& vtx_module, int& vtx_plane )
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
    }
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
