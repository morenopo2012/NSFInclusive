//#ifndef MNV_NUKECC_cxx
//#define MNV_NUKECC_cxx 1

#ifndef MNV_NUKECC_BINNING_CXX
#define MNV_NUKECC_BINNING_CXX 1

#include "../include/NukeCC_Binning.h"
#include "PlotUtils/TargetUtils.h"
#include "../include/CVUniverse.h"
//#include "include/GlobalIncludes.h"
#include "../include/CommonIncludes.h"
#include <PlotUtils/MnvNormalization.h>
#include <PlotUtils/NuclModUtils.h>
#include <PlotUtils/FluxReweighter.h>
//#include <PlotUtils/FluxReweighterWithWiggleFit.h>
#include <PlotUtils/MnvNuclearModelWeight.h>
//#include <PlotUtils/MnvNormalizerME.h>
//=======
#include "PlotUtils/MinosMuonEfficiencyCorrection.h"

#include "TFileCollection.h"


#include "TVector3.h"

using namespace std;

using namespace NUKECC_ANA;

NukeCC_Binning::NukeCC_Binning()
{
}

NukeCC_Binning::~NukeCC_Binning()
{
}


//======================
// get singleton
//======================
NukeCC_Binning& NukeCC_Binning::Get()
{
    static NukeCC_Binning singleton;
    return singleton;
}


void NukeCC_Binning::AddBins( std::vector<double>& binsLowEdge, const int nBins, const double binWidth ) const
{
    //! If the bins have no defined low entry, assume it should be 0
    double x = 0.;
    if( binsLowEdge.size() > 0 )
        x = binsLowEdge[ binsLowEdge.size() - 1 ];
    else
        binsLowEdge.push_back( x );

    //! add nBins spaced by this width
    for( int iBin = 0; iBin < nBins; ++iBin )
    {
        x += binWidth;
        binsLowEdge.push_back( x );
    }
}

void NukeCC_Binning::AddBins( std::vector<double>& binsLowEdge, const double binWidth, const double lowBin, const double highBin ) const
{
    double x = lowBin;
    if( binsLowEdge.size() && lowBin == binsLowEdge.back() )
        x += binWidth;
    while( x <= highBin )
    {
        binsLowEdge.push_back( x );
        x += binWidth;
    }
}


std::vector<double>NukeCC_Binning::GetBins( const std::string& var, BinType::t_BinType type ) const
{
    vector<double> binsLowEdge;
    if(      type == BinType::kDIS )
        binsLowEdge = GetDISBins( var );
    else if( type == BinType::kEnergy )
        binsLowEdge = GetEnergyBins( var );
    else if( type == BinType::kResolution )
        binsLowEdge = GetResolutionBins( var );
    else if( type == BinType::kResidual )
        binsLowEdge = GetResidualBins( var );
    else if( type == BinType::kPosition )
        binsLowEdge = GetPosBins( var );
    else if( type == BinType::kSideband )
        binsLowEdge = GetSidebandBins( var );
    else if( type == BinType::kTruth )
        binsLowEdge = GetTruthBins( var );
    else if( type == BinType::kDeltaE )
        binsLowEdge = GetDeltaEBins( var );

    if( binsLowEdge.empty() )
        Warning( "NukeCC_Binning::GetBins", Form("No bins found for variable %s, bintype %d", var.c_str(), int(type) ) );

    return binsLowEdge;
}

std::vector<double> NukeCC_Binning::GetEnergyBins( const std::string& var, bool DISbins /*=true*/ ) const
{
    vector<double> binsLowEdge;
    //if( DISbins ) return binsLowEdge = GetDISBins( var );

    if( var == "Q2"){
      double tmpBins[] = { 0, 0.05, 0.1, 0.2, 0.3, 0.5, 0.75,1.0, 1.5, 3.0, 5.0, 8.0, 50. };
      //double tmpBins[] = { 0, 0.5, 1.0, 1.5, 5.0, 8.0};
      binsLowEdge.assign( tmpBins, tmpBins + sizeof(tmpBins) / sizeof(tmpBins[0]) );
    }
    else if( var == "W" ){
      double tmpBins[]  = {0, 1.1, 1.3, 1.7, 2., 2.5, 3.5, 5., 8., 10.};
      //double tmpBins[]  = {0, 2., 5., 8., 10.};
      binsLowEdge.assign( tmpBins, tmpBins + sizeof(tmpBins) / sizeof(tmpBins[0]) );
    }
    else if( var == "Eavail"){
     //double tmpBins[] ={0, 0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14, 0.16, 0.2, 0.25,0.299, 0.35, 0.399, 0.499, 0.599, 0.799, 0.999, 1.3, 1.9, 2.1};
     //double tmpBins[] = {0.0, 0.04, 0.08, 0.12, 0.16, 0.24, 0.32, 0.4, 0.6, 0.8, 1.0, 1.2}; //Marvin's Binning II for CrossSectionExtraction
     double tmpBins[] = {0.0, 0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14, 0.16, 0.2, 0.25,0.299, 0.35, 0.399,  0.45, 0.499,  0.55, 0.599, 0.65, 0.699, 0.75, 0.799, 0.85, 0.899, 0.95, 0.999, 1.3}; //Marvin's thesis
     //double tmpBins[] ={0, 0.05, 0.10 ,0.15, 0.20, 0.25, 0.3, 0.35, 0.4 ,0.45, 0.5 ,0.55,  0.6 , 0.65,  0.7 ,  0.75,  0.8 ,  0.85, 0.9 ,  0.95, 1.0, 1.3, 1.9, 2.1}; //MEC binning 
     //double tmpBins[] ={0.01,0.02, 0.04, 0.06, 0.08,0.1,0.12,0.14,0.16, 0.2, 0.25, 0.299, 0.35, 0.4 ,0.45, 0.5 ,0.55,  0.6 , 0.65,  0.7 ,  0.75,  0.8 ,  0.85, 0.9 ,  0.95, 1.0, 1.3, 1.9, 2.1}; //MEC binning
     //double tmpBins[] = {0, 0.03, 0.06, 0.09, 0.12 , 0.15, 0.18, 0.21,0.24,0.27,0.3,0.33,0.36,0.39,0.42,0.45,0.48,0.51,0.54,0.57,0.6,0.63,0.66,0.69,0.72,0.75,0.78,0.81,0.84,0.87,0.9,0.93,0.96,0.99,1.02,1.05,1.08,1.11,1.15 };
     //double tmpBins[] = {0.005,0.01,0.015,0.02,0.025,0.03,0.035,0.04,0.045,0.05,0.055,0.06,0.065,0.07,0.075,0.08,0.085,0.09,0.095,0.1,0.105,0.11,0.115,0.12,0.125,0.13,0.135,0.14,0.145,0.15,0.155,0.16,0.165,0.17,0.175,0.18,0.185,0.19,0.195,0.2,0.205,0.21,0.215,0.22,0.225,0.23,0.235,0.24,0.245,0.25,0.255,0.26,0.265,0.27,0.275,0.28,0.285,0.29,0.295,0.3,0.305,0.31,0.315,0.32,0.325,0.33,0.335,0.34,0.345,0.35,0.355,0.36,0.365,0.37,0.375,0.38,0.385,0.39,0.395,0.4,0.405,0.41,0.415,0.42,0.425,0.43,0.435,0.44,0.445,0.45,0.455,0.46,0.465,0.47,0.475,0.48,0.485,0.49,0.495,0.5,0.505,0.51,0.515,0.52,0.525,0.53,0.535,0.54,0.545,0.55,0.555,0.56,0.565,0.57,0.575,0.58,0.585,0.59,0.595,0.6,0.605,0.61,0.615,0.62,0.625,0.63,0.635,0.64,0.645,0.65,0.655,0.66,0.665,0.67,0.675,0.68,0.685,0.69,0.695,0.7,0.705,0.71,0.715,0.72,0.725,0.73,0.735,0.74,0.745,0.75,0.755,0.76,0.765,0.77,0.775,0.78,0.785,0.79,0.795,0.8,0.805,0.81,0.815,0.82,0.825,0.83,0.835,0.84,0.845,0.85,0.855,0.86,0.865,0.87,0.875,0.88,0.885,0.89,0.895,0.9,0.905,0.91,0.915,0.92,0.925,0.93,0.935,0.94,0.945,0.95,0.955,0.96,0.965,0.97,0.975,0.98,0.985,0.99,0.995,1,1.005,1.01,1.015,1.02,1.025,1.03,1.035,1.04,1.045,1.05,1.055,1.06,1.065,1.07,1.075,1.08,1.085,1.09,1.095,1.1,1.105,1.11,1.115,1.12,1.125,1.13,1.135,1.14,1.145,1.15,1.155,1.16,1.165,1.17,1.175,1.18,1.185,1.19,1.195,1.2,1.205,1.21,1.215,1.22,1.225,1.23,1.235,1.24,1.245,1.25,1.255,1.26,1.265,1.27,1.275,1.28,1.285,1.29,1.295,1.3,1.305,1.31,1.315,1.32,1.325,1.33,1.335,1.34,1.345,1.35,1.355,1.36,1.365,1.37,1.375,1.38,1.385,1.39,1.395,1.4,1.405,1.41,1.415,1.42,1.425,1.43,1.435,1.44,1.445,1.45,1.455,1.46,1.465,1.47,1.475,1.48,1.485,1.49,1.495,1.5,1.505,1.51,1.515,1.52,1.525,1.53,1.535,1.54,1.545,1.55,1.555,1.56,1.565,1.57,1.575,1.58,1.585,1.59,1.595,1.6}; //Rik Binning
     binsLowEdge.assign( tmpBins, tmpBins + sizeof(tmpBins) / sizeof(tmpBins[0]) );
    }
    else if( var == "q3"){
     //double tmpBins[] = {0.01, 0.2, 0.3, 0.4, 0.6, 0.9, 1.2, 1.5, 2.1, 3.0};
     double tmpBins[] = {0.0, 0.2, 0.3, 0.4, 0.6, 0.9, 1.2};
     //double tmpBins[] = {0.01, 0.12, 0.14, 0.18, 0.20, 0.24, 0.26, 0.30, 0.34, 0.38, 0.40,
     //                    0.41, 0.44, 0.46, 0.49, 0.50, 0.54, 0.58, 0.60, 0.64, 0.68, 0.70, 0.80,
     //                    0.81, 0.84, 0.86, 0.89, 0.92, 0.94, 0.98, 1.00, 1.10, 1.20,
     //                    1.29, 1.40, 1.50, 1.60,
     //                    1.70, 1.80, 1.90, 2.00,
     //                    2.10, 2.20, 2.30, 2.40,
      //                   2.50, 2.60, 2.70,
        //                 2.81, 2.89, 3.00};

     //double tmpBins[] = {0.0,0.2,0.3,0.4,0.5,0.6,0.8,0.9,1.0,1.1};
     binsLowEdge.assign( tmpBins, tmpBins + sizeof(tmpBins) / sizeof(tmpBins[0]) );
    }
    else if( var == "KE_preFSI"){
     double tmpBins[] = {0.0,40.0,80.0,120.0,160.0,200.0,240.0,280.0,320.0,360.0,400.0,440.0,480.0,520.0,560.0,600.0,640.0,680.0,720.0,760.0,800.0,840.0,880.0,920.0,960.0,1000.0,1040.0,1080.0,1120.0,1160.0,1200.0,1240.0,1280.0,1320.0,1360.0,1400.0,1440.0,1480.0,1520.0,1560.0,1600.0,1640.0,1680.0,1720.0,1760.0,1800.0,1840.0,1880.0,1920.0,1960.0};
    }
    else if( var == "Emu" || var == "Emum" ){
      //double tmpBins[]  = { 0. ,1. ,2. ,3. ,4. ,5. ,7. ,9. ,12. ,15. ,18. ,22. ,36. ,50. ,75. ,100. ,120. };
      double tmpBins[]  = {2.0, 3.5, 4.75, 5.5, 7.5, 10., 13., 16., 20., 25., 35., 50 };
      binsLowEdge.assign( tmpBins, tmpBins + sizeof(tmpBins) / sizeof(tmpBins[0]) );
    }
    else if(var == "pZmu"){
       double tmpBins[] = { 2.0, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 8.0, 10.0, 15.0, 20.0, 40.0 };
       binsLowEdge.assign(tmpBins, tmpBins+sizeof(tmpBins) / sizeof(tmpBins[0]));
    }
    else if( var == "Enu" ){
      double tmpBins[]  = { 2.0, 3.0, 4.0, 5.0, 6.25, 7.5, 8.75, 10., 12.5, 15., 20., 25.0, 30., 40., 50. };
      binsLowEdge.assign( tmpBins, tmpBins + sizeof(tmpBins) / sizeof(tmpBins[0]) );
    }
    else if( var == "Ehad" ){
      double tmpBins[]  = { 0., 0.25, 0.75, 1.5, 2.5, 3.5, 5.0, 8.0, 12., 18., 25. };
      binsLowEdge.assign( tmpBins, tmpBins + sizeof(tmpBins) / sizeof(tmpBins[0]) );
    }
    else if( var == "x" ){

      binsLowEdge.push_back(0.);
      binsLowEdge.push_back(0.05);
      binsLowEdge.push_back(0.1);
      binsLowEdge.push_back(0.2);
      binsLowEdge.push_back(0.3);
      binsLowEdge.push_back(0.5);
      binsLowEdge.push_back(1.75);

    }
    else if( var == "y" )
    {
       double tmpBins[] = { 0, 0.03, 0.15, 0.3, 0.45, 0.6, 1.};
       binsLowEdge.assign(tmpBins,tmpBins + sizeof(tmpBins) / sizeof(tmpBins[0]) );
    }


    else if( var == "CCQE-Recoil" || var == "VtxEnergy" )
    {
        //.5 * recoil for now
        binsLowEdge = GetBins( "Ehad" );
        for( vector<double>::iterator i = binsLowEdge.begin(); i != binsLowEdge.end(); ++i )
            *i = .5*(*i);
    }
    else if( var == "ETheta" )
    {
        AddBins( binsLowEdge, 4, .02 );
        AddBins( binsLowEdge, 2, .04 );
        AddBins( binsLowEdge, 2, .08 );
    }
    else if( var == "PhiMu" )
    {
        //4 deg bins from -180 to 180
        AddBins( binsLowEdge, 12., -180, 180 );
    }
    else if( var == "ThetaMu")
    {
        // 1 deg bins up to 8
        //AddBins( binsLowEdge, 10, 1.5 );
        AddBins( binsLowEdge, 8, 1.0);
        // 2 deg bins up to 14,
        //AddBins( binsLowEdge, 3, 2.0 );
        // 3 deg bins up to 26
        AddBins( binsLowEdge, 4, 3. );
        //--- overflow --//
        // 5 deg bins up to 51
        AddBins( binsLowEdge, 5, 5. );
        // 2 20 deg bins up to 91
        AddBins( binsLowEdge, 2, 20. );
    }
    else if( var == "ThetaMuX" || var == "ThetaMuY" )
    {
        //3 deg bins from -29 to -20
        AddBins( binsLowEdge, 3., -29., -20. );
        //---underflow---//
        //2 deg bins from -20 to -10
        AddBins( binsLowEdge, 2., -20., -10. );
        //1 deg bins from -10 to 10
        AddBins( binsLowEdge, 1., -10., 10. );
        //2 deg bins from 10 to 20
        AddBins( binsLowEdge, 2., 10., 20. );
        //--- overflow ---//
        //3 deg bins from 20 to 29
        AddBins( binsLowEdge, 3., 20., 29. );
    }
    else if( var == "CosThetaMu" )
    {
        binsLowEdge = GetBins( "ThetaMu" );
        for( vector<double>::iterator i = binsLowEdge.begin(); i != binsLowEdge.end(); ++i )
            *i = cos( deg_to_rad * (*i) );
        reverse( binsLowEdge.begin(), binsLowEdge.end() );
    }

    else if( var == "NTracks" )
    {
        // 1-10
        AddBins( binsLowEdge, 10, 1. );
    }

    else if( var == "plane" || var == "planeDNN"  )
    {
        // 0-66
        double tmpBins[] = { 0., 1., 2., 3., 4., 5., 6., 7., 8., 9., 11, 12, 13., 14., 15., 16., 17., 18., 19., 21., 22., 23., 24., 25., 26., 27., 28., 29., 33., 34., 35., 36., 37., 38., 39., 40., 41., 42., 43., 44., 45., 46., 47., 48., 49., 51., 52., 53., 54., 55., 57., 58., 59., 60., 61., 62., 63., 64., 65., 66. };
        //double tmpBins[] = { 0., 1., 2., 3., 4., 5., 6., 7., 8., 9., 11., 12., 16., 17., 18., 19., 21., 22., 28., 29., 33., 34., 35., 36., 37., 38., 39., 40., 41., 42., 43., 44., 45., 46., 47., 48., 49., 51., 52., 53., 54., 55., 57., 58., 59., 60., 61., 62., 63., 64., 65., 66. };
        binsLowEdge.assign( tmpBins, tmpBins + sizeof(tmpBins) / sizeof(tmpBins[0]) );
    }
    else if( var == "vtxz" || var == "vtxy" || var == "vtxx"){
        binsLowEdge = GetPosBins( var );

    }

    else if(var == "muonPt"){
       double tmpBins[] = { 0.0, 0.075, 0.15, 0.25, 0.325, 0.4, 0.475, 0.55, 0.7, 0.85, 1.0, 1.25, 2.5};
       binsLowEdge.assign(tmpBins,tmpBins + sizeof(tmpBins) / sizeof(tmpBins[0]) );
    }


    else if( var == "vtxz_all" ){
      // 0.1 GeV bins from 0 to 100 GeV
      //AddBins( binsLowEdge, 0.05, 0., 1000. );
      double tmpBins[]  =
        {4293.04/10,
        4313.68/10,
        4337.25/10,
        4357.9/10,
        4381.47/10,
        4402.11/10,
        4425.68/10,
        4446.33/10,
        4514.11/10,
        4534.76/10,
        4558.33/10,
        4578.97/10,
        4602.54/10,
        4623.19/10,
        4646.76/10,
        4667.4/10,
        4735.19/10,
        4755.83/10,
        4779.4/10,
        4800.05/10,
        4823.62/10,
        4844.26/10,
        4867.83/10,
        4888.48/10,
        5000.48/10,
        5021.12/10,
        5044.69/10,
        5065.34/10,
        5088.91/10,
        5109.55/10,
        5133.12/10,
        5153.77/10,
        5456.74/10,
        5477.38/10,
        5500.95/10,
        5521.6/10,
        5545.17/10,
        5565.81/10,
        5589.38/10,
        5610.02/10,
        5677.81/10,
        5698.45/10,
        5722.03/10,
        5742.67/10,
        5810.45/10,
        5831.1/10,
        5855.68/10,
        5876.33/10,
        5900.91/10,
        };
      binsLowEdge.assign( tmpBins, tmpBins + sizeof(tmpBins) / sizeof(tmpBins[0]) );
    }



    return binsLowEdge;

}


std::vector<double> NukeCC_Binning::GetDISBins(const std::string& var ) const
{


    vector<double> binsLowEdge;
    if( var == "Q2"){
        //double tmpBins[] = { 1.0, 1.25, 3.0, 5.0, 8.0, 30. };
        double tmpBins[] = { 4.0, 5.0, 8.0, 10.0, 30. };
        //double tmpBins[] = { 1.0, 1.5, 1.7, 2.0, 3.0, 4., 10.0, 20.0 };
        binsLowEdge.assign( tmpBins, tmpBins + sizeof(tmpBins) / sizeof(tmpBins[0]) );
    }
    else if( var == "W" )
    {
        //for DIS: one underflow bin (0 - 2 GeV) same otherwise
        //double tmpBins[]  = {2., 2.5, 3.5, 5., 8., 16.};
        double tmpBins[]  = {3.5, 5., 8., 12.0, 16.};
        //double tmpBins[]  = {1., 1.7, 2., 3.5, 5.0, 8., 16.};
        binsLowEdge.assign( tmpBins, tmpBins + sizeof(tmpBins) / sizeof(tmpBins[0]) );
    }
    else if( var == "Emu" || var == "Emum" )
    {
        //for DIS: one underflow bin (0 - 2 GeV) same otherwise
        double tmpBins[]  = { 2., 3., 4., 5., 6., 8., 10., 13., 16., 20., 25., 35., 50.};
        //double tmpBins[]  = { 3.5, 5, 6., 8., 10., 13., 16., 20., 25., 35., 50.};
        //double tmpBins[]  = {2.0, 3.5, 4.75, 5.5, 7.5, 10., 13., 16., 20., 25., 35., 50 };
        binsLowEdge.assign( tmpBins, tmpBins + sizeof(tmpBins) / sizeof(tmpBins[0]) );
    }
    else if(var == "pZmu"){
       double tmpBins[] = { 2.0, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 8.0, 10.0, 15.0, 20.0, 40.0 };
       binsLowEdge.assign(tmpBins, tmpBins+sizeof(tmpBins) / sizeof(tmpBins[0]));
    }
    else if( var == "Enu" )
    {
        double tmpBins[]  = { 5., 7.5, 10., 12.5, 15., 20., 30., 40., 50. };
        binsLowEdge.assign( tmpBins, tmpBins + sizeof(tmpBins) / sizeof(tmpBins[0]) );
    }
    else if( var == "Ehad" )
    {
        //for DIS: one underflow bin (0 - 2 GeV) same otherwise
        //double tmpBins[]  = { 3., 5., 8., 11., 15., 20., 25. };
        double tmpBins[]  = { 1.5, 3., 5., 8., 11., 15., 20., 25. };
        binsLowEdge.assign( tmpBins, tmpBins + sizeof(tmpBins) / sizeof(tmpBins[0]) );
    }
    else if( var == "x" )
    {

        //from Jorge
        //I guess this should be overflow since we don't have anything generated below 0.005
        //binsLowEdge.push_back(0.);

        // 0-.1 (shadowing)
        binsLowEdge.push_back(0.005);
        binsLowEdge.push_back(0.05);
        binsLowEdge.push_back(0.1);

        //.1-.3 (antishadowing)
        binsLowEdge.push_back(0.2);
        binsLowEdge.push_back(0.3);

        //.3-.8 (EMC)
        binsLowEdge.push_back(0.5);
        binsLowEdge.push_back(0.8);
        //binsLowEdge.push_back(1.5);
    }
    else if( var == "y" )
    {
       double tmpBins[] = { 0, 0.03, 0.15, 0.3, 0.45, 0.6, 1.};
       binsLowEdge.assign( tmpBins, tmpBins + sizeof(tmpBins) / sizeof(tmpBins[0]) );
    }

    else if( var == "planeDNN"  )
    {
        // 0-214
        double tmpBins[] = { 0., 1., 2., 3., 4., 5., 6., 7., 8., 9., 11, 12, 13., 14., 15., 16., 17., 18., 19., 21., 22., 23., 24., 25., 26., 27., 28., 29., 33., 34., 35., 36., 37., 38., 39., 40., 41., 42., 43., 44., 45., 46., 47., 48., 49., 51., 52., 53., 54., 55., 57., 58., 59., 60., 61., 62., 63., 64., 65., 66., 67., 68., 69., 70., 71., 72., 73., 74., 75., 76., 77., 78., 79., 80., 81., 82., 83., 84., 85., 86., 87., 88., 89., 90., 91., 92., 93., 94., 95., 96., 97., 98., 99., 100., 101., 102., 103., 104.,105.,106.,107.,108.,109.,110.,111.,112.,113.,114.,115.,116.,117.,118.,119.,120.,121.,122.,123.,124.,125.,126.,127.,128.,129.,130.,131.,132.,133.,134.,135.,136.,137.,138.,139.,140.,141.,142.,143.,144.,145.,146.,147.,148.,149.,150.,151.,152.,153.,154.,155.,156.,157.,158.,159.,160.,161.,162.,163.,164.,165.,166.,167.,168.,169.,170.,171.,172.,173.,174.,175.,176.,177.,178.,179.,180.,181.,182.,183.,184.,185.,186.,187.,188.,189.,190.,191.,192.,193.,194.,195.,196.,197.,198.,199.,200.,201.,202.,203.,204.,205.,206.,207.,208.,209.,210.,211.,212.,213.,214. };
        //double tmpBins[] = { 0., 1., 2., 3., 4., 5., 6., 7., 8., 9., 11., 12., 16., 17., 18., 19., 21., 22., 28., 29., 33., 34., 35., 36., 37., 38., 39., 40., 41., 42., 43., 44., 45., 46., 47., 48., 49., 51., 52., 53., 54., 55., 57., 58., 59., 60., 61., 62., 63., 64., 65., 66. };
        binsLowEdge.assign( tmpBins, tmpBins + sizeof(tmpBins) / sizeof(tmpBins[0]) );
    }


    else if( var == "vtxz" || var == "vtxy" || var == "vtxx"){
        binsLowEdge = GetPosBins( var );

    }

    else if( var == "ThetaMu")
    {
        // 1 deg bins up to 8
        //AddBins( binsLowEdge, 10, 1.5 );
        AddBins( binsLowEdge, 8, 1.0);
        // 2 deg bins up to 14,
        //AddBins( binsLowEdge, 3, 2.0 );
        // 3 deg bins up to 26
        AddBins( binsLowEdge, 4, 3. );
        //--- overflow --//
        // 5 deg bins up to 51
        AddBins( binsLowEdge, 5, 5. );
        // 2 20 deg bins up to 91
        AddBins( binsLowEdge, 2, 20. );
    }
    return binsLowEdge;

}


std::vector<double> NukeCC_Binning::GetSidebandBins(const std::string& var ) const
{
    vector<double> binsLowEdge;

    if( var == "W" )
    {
        AddBins( binsLowEdge, .05, 0., 25. );
        //for DIS: one underflow bin (0 - 2 GeV) same otherwise
        //double tmpBins[]  = {0, 1.3, 1.5, 1.7, 2.0, 2.5, 3., 3.5, 4., 4.5, 5., 6., 10., 25, 100.};
        //binsLowEdge.assign( tmpBins, tmpBins + sizeof(tmpBins) / sizeof(tmpBins[0]) );
    }
    else if( var == "Q2" )
    {
        AddBins( binsLowEdge, .05, 0., 25. );
        //for DIS: one underflow bin (0 - 2 GeV) same otherwise
        //double tmpBins[]  = {0, 1.3, 1.5, 1.7, 2.0, 2.5, 3., 3.5, 4., 4.5, 5., 6., 10., 25, 100.};
        //binsLowEdge.assign( tmpBins, tmpBins + sizeof(tmpBins) / sizeof(tmpBins[0]) );
    }
    else
        binsLowEdge =  GetEnergyBins( var );

    return binsLowEdge;

}

std::vector<double> NukeCC_Binning::GetResidualBins( const std::string& var ) const
{
    vector<double> binsLowEdge;
    if( var == "Emu" )
    {
        // .2 GeV bins from -2 to 2
        AddBins( binsLowEdge, 0.2, -6.1, 6.1 );
    }
    else if( var == "Enu" )
    {
        // .2 GeV bins from -5 to 5
        AddBins( binsLowEdge, 0.2, -4.1, 4.1 );
    }
    else if( var == "Emum")
    {
        // .2 GeV bins from -2 to 2
        AddBins( binsLowEdge, .2, -2.1, 2.1 );
    }
    else if( var == "Ehad" || var == "CCQE-Recoil" || var == "VtxEnergy" )
    {
        AddBins( binsLowEdge, 0.4, -10, 10 );
    }
    else if( var == "PhiMu" )
    {
        // -10 to 10 deg in bins of .25 deg
        AddBins( binsLowEdge, .5, -15.25, 15.15 );
    }
    else if( var == "ThetaMu" || var == "ThetaMuX" || var == "ThetaMuY" )
    {
        // -5 deg to 5 deg in bins of .2 deg
        AddBins( binsLowEdge, .2, -5.1, 5.1 );
    }
    else if( var == "x" )
    {
        AddBins( binsLowEdge, .05, -1.025, 1.025 );
    }
    else if( var == "y" )
    {
        AddBins( binsLowEdge, 0.02, -.31, .31 );
   }
    else if( var == "Q2" )
    {
        AddBins( binsLowEdge, .1, -1.05, 1.05 );
    }
    else if( var == "W" )
    {
        AddBins( binsLowEdge, .2, -2.1, 2.1 );
    }
    else if( var == "vtxt" )
    {
        // 10ns bins from -200 to 200
        AddBins( binsLowEdge, 10., -205., 205.);
    }
    else if( //all transverse directions get the same bins
            string::npos != var.find("vtxx") ||
            string::npos != var.find("vtxy") ||
            string::npos != var.find("vtxd") ||
            string::npos != var.find("vtxu") ||
            string::npos != var.find("vtxc")
            )
    {
        AddBins( binsLowEdge, .05, -2.525,2.525);
    }
    else if( string::npos != var.find("vtxz" ) )
    {
        //Try to make these ~1 strip wide
        AddBins( binsLowEdge, 17, -150, 420 );
    }
    else if ( string::npos != var.find("vtxr" ) )
    {
        AddBins( binsLowEdge, .05, -.075, 2.425);
    }


    return binsLowEdge;
}

std::vector<double> NukeCC_Binning::GetResolutionBins( const std::string& var ) const
{
    vector<double> binsLowEdge;

    AddBins( binsLowEdge, .05, -1.25, 2.25 );

    return binsLowEdge;
}

std::vector<double> NukeCC_Binning::GetDeltaEBins( const std::string& var ) const
{
    vector<double> binsLowEdge;

    AddBins( binsLowEdge, .01, -0.35, 1.25 );

    return binsLowEdge;
}


//==========================================
// start of Det bins
//==========================================
const unsigned int NukeCC_Binning::fNDetBins;

double NukeCC_Binning::fDetBins[126] = {
    424.11224, //-5
    428.22447, //-4
    432.64595, //-3
    437.06743, //-2
    441.48891, //tgt1
    445.91039, //0
    450.33188,
    454.75336,
    459.17484, //3
    463.59632, //tgt2
    468.01780, //5
    472.43928,
    476.86077,
    481.28225, //8
    485.70373, //tgt3a
    490.12521, //tgt3b
    494.54669, //11
    498.96817,
    503.38965,
    507.81114,
    512.23262,
    516.44336,
    544.59410,
    549.01558,
    553.43706,
    557.85854,
    562.28002,
    566.70151,
    571.12299,
    575.54447,
    579.91524,
    584.43814,
    588.96104,
    593.48394,
    598.00684,
    602.52974,
    607.05264,
    611.57554,
    616.09844,
    620.62134,
    625.14424,
    629.66714,
    634.19004,
    638.71294,
    643.23584,
    647.75874, //
    652.28164,
    656.80454,
    661.32744, //40
    665.85034,
    670.37324,
    674.89614,
    679.41904,
    683.94194,
    688.46484,
    692.98774,
    697.51064,
    702.03354,
    706.55644,
    711.07934,
    715.60224,
    720.12514,
    724.64804,
    729.17094,
    733.69384,
    738.21674,
    742.73964,
    747.26254,
    751.78544,
    756.30834,
    760.83124,
    765.35414,
    769.87704,
    774.39994,
    778.92284,
    783.44574,
    787.96864,
    792.49154,
    797.01444,
    801.53734,
    806.06024,
    810.58314,
    815.10604,
    819.62894,
    824.15184,
    828.67474,
    833.19764,
    837.72054,
    842.24344,
    846.76634,
    851.28924,
    855.81214,
    860.35635,
    864.83663, //end of mod 85
    869.38084,
    873.90374,
    878.42664,
    882.94954,
    887.47244,
    891.99534,
    896.51824,
    901.04114,
    905.56404,
    910.08694,
    914.60984,
    919.13274,
    923.65564,
    928.17854,
    932.70144,
    937.22434,
    941.74724,
    946.27014,
    950.79304,
    955.31594,
    959.83884,
    964.36174,
    968.88464,
    973.40754,
    977.93044,
    982.45334,
    986.97624,
    991.49914,
    996.02204,
    1000.54494,
    1005.06784,
    1009.59074
};
//==========================================
// end of Det bins
//==========================================


//distances here in cm, ns
std::vector<double> NukeCC_Binning::GetPosBins( const std::string& var, const int targetID /* = -1 */ ) const
{
    vector<double> binsLowEdge;

    const int xyspace = 2.5; //cm
    const double widthXUDC = 200.;
    const double widthY = 220.;
    if( string::npos != var.find("vtxt") )
    {
        // from 0 to 10000 in 100ns intervals
        static const double deltaT = 100;
        AddBins( binsLowEdge, 10000/deltaT, deltaT );
    }
    if(
       string::npos != var.find("vtxx" ) ||
       string::npos != var.find("vtxu" ) ||
       string::npos != var.find("vtxd" ) ||
       string::npos != var.find("vtxc" )
       )
    {
        binsLowEdge.push_back( -widthXUDC/2. );
        AddBins( binsLowEdge, widthXUDC/xyspace, xyspace );
    }
    else if( string::npos != var.find("vtxy" ) )
    {
        binsLowEdge.push_back( -widthY/2. );
        AddBins( binsLowEdge, widthY/xyspace, xyspace );
    }
    else if (string::npos != var.find("vtxr") )
    {
        binsLowEdge.push_back( 0 );
        AddBins( binsLowEdge, 1.17E2/xyspace, xyspace );
    }
    else if( string::npos != var.find("vtxz") )
    {
        int zspace = 1.;//cm
        if( targetID == -1 )
        {
            binsLowEdge.assign(fDetBins, fDetBins + fNDetBins); //return the detector bins
        }
        else if( targetID == 1 )
        {
            // from 440 to 460
            binsLowEdge.push_back( 440. );
            AddBins( binsLowEdge, (460-440)/zspace, zspace );
        }
        else if( targetID == 2 )
        {
            // from 462 to 482
            binsLowEdge.push_back( 462. );
            AddBins( binsLowEdge, (482-462)/zspace, zspace );
        }
        else if( targetID == 3 )
        {
            // from 485 to 505
            binsLowEdge.push_back( 485 );
            AddBins( binsLowEdge, (505-485)/zspace, zspace );
        }
        else if( targetID == 4 )
        {
            // from 555 to 575
            binsLowEdge.push_back( 555 );
            AddBins( binsLowEdge, (575-555)/zspace, zspace );
        }
        else if( targetID == 5 )
        {
            // from 575 to 585
            binsLowEdge.push_back( 575 );
            AddBins( binsLowEdge, (585-575)/zspace, zspace );
        }
        else
        {
            const int totalLowestMod = 27;
            const double totalLowestZ = 595.;
            const double deltaZPerMod = 4.5; //4.5cm per module

            int thisLowestMod = targetID % 1000;
            if( targetID < 100 )
            {
                int fauxSet = ( targetID - targetID%10 ) / 10;
                thisLowestMod = totalLowestMod + 6*(fauxSet-1);
            }
            double lowZ = totalLowestZ + deltaZPerMod * ( thisLowestMod - totalLowestMod );
            double highZ = lowZ + deltaZPerMod * 6;

            binsLowEdge.push_back( lowZ );
            AddBins( binsLowEdge, (int(highZ)-int(lowZ))/zspace, zspace );
        }
    }
    else if( "DistToTarg" == var )
    {
        const double targPos = GetBGZ( targetID );

        //skip this many modules before getting to the first (last) active faux module
        //mod 85 is the most DS module included in module z bins
        //mod 80 is the most DS module of the most DS faux target
        const int firstPlaneCount = 6;

        int counter = 0;
        int planesSinceLast = 0; //trick loop into making an active bin at the first chance

        std::vector<double> planeBins = GetPosBins( "vtxz", -1 ); //get bins by plane
        for( std::vector<double>::reverse_iterator i = planeBins.rbegin(); i != planeBins.rend(); ++i, ++counter )
        {
            if( counter < firstPlaneCount )
                continue;

            //faux targets every 6 modules
            if( 0 == planesSinceLast % 6 )
                binsLowEdge.push_back( targPos - *i );

            ++planesSinceLast;
        }
    }
    else if( var == "NTracks" )
    {
        // 1-10
        AddBins( binsLowEdge, 10, 1. );
    }
    else if( var == "plane" || var == "planeDNN"  )
    {
        // 0-66
        double tmpBins[] = { 0., 1., 2., 3., 4., 5., 6., 7., 8., 9., 11, 12, 13., 14., 15., 16., 17., 18., 19., 21., 22., 23., 24., 25., 26., 27., 28., 29., 33., 34., 35., 36., 37., 38., 39., 40., 41., 42., 43., 44., 45., 46., 47., 48., 49., 51., 52., 53., 54., 55., 57., 58., 59., 60., 61., 62., 63., 64., 65., 66. };
        binsLowEdge.assign( tmpBins, tmpBins + sizeof(tmpBins) / sizeof(tmpBins[0]) );
    }

    return binsLowEdge;
}

double NukeCC_Binning::GetBGZ( const int targetID ) const
{
    //all in cm
    if( targetID == 1 )  return 452.4;
    if( targetID == 2 )  return 474.5 ;
    if( targetID == 3 )  return 501.05;
    if( targetID == 4 )  return 568.75;
    if( targetID == 5 )  return 582.05;

    Error( "GetBGZ", Form("Cannot get background Z for targetID = %d", targetID ) );
    throw 1;
}

vector<double> NukeCC_Binning::GetTruthBins(const std::string& var) const{

    vector<double> binsLowEdge;

    if(var == "Enu" || var == "Emu" || var == "Emum" ){
        AddBins( binsLowEdge, 5, 1. );
        //2 GeV bin up to 7 GeV
        /*AddBins( binsLowEdge, 1, 2. );
         //3 GeV bin up to 10 GeV
         AddBins( binsLowEdge, 1, 3. );*/
        //5GeV bins up to 30GeV
        AddBins( binsLowEdge, 5, 5. );
        //10GeV Bins up to 50 GeV
        AddBins( binsLowEdge, 2, 10 );
        AddBins( binsLowEdge, 1, 50. ); //50-100
        //Overflow
        AddBins( binsLowEdge, 1, 100);
    }

    else if(var == "Ehad")
        NukeCC_Binning::Get().AddBins(binsLowEdge, 100, 0.5);

    else if(var == "Q2"){
        NukeCC_Binning::Get().AddBins(binsLowEdge, 200, 0.1);

    }

    else if(var == "x" ){
        NukeCC_Binning::Get().AddBins(binsLowEdge, 110, 0.1);

    }

    else if(var == "xGen" ){
        NukeCC_Binning::Get().AddBins(binsLowEdge, 110, 0.1);

    }

    else if(var == "y"){
        NukeCC_Binning::Get().AddBins(binsLowEdge, 110, 0.01);

    }

    else if(var == "W"){
        NukeCC_Binning::Get().AddBins(binsLowEdge, 200, 0.1);

    }

    else if(var == "ThetaMu"){
        NukeCC_Binning::Get().AddBins(binsLowEdge, 180, 1.0);

    }

    else if (var == "NTracks"){

        AddBins( binsLowEdge, 10, 1. );
    }

    else if (var == "plane" || var == "planeDNN" ){
        // 0-66
        double tmpBins[] = { 0., 1., 2., 3., 4., 5., 6., 7., 8., 9., 11, 12, 13., 14., 15., 16., 17., 18., 19., 21., 22., 23., 24., 25., 26., 27., 28., 29., 33., 34., 35., 36., 37., 38., 39., 40., 41., 42., 43., 44., 45., 46., 47., 48., 49., 51., 52., 53., 54., 55., 57., 58., 59., 60., 61., 62., 63., 64., 65., 66. };
        binsLowEdge.assign( tmpBins, tmpBins + sizeof(tmpBins) / sizeof(tmpBins[0]) );
    }

    else{
        Error( "NukeCC_Binning::GetTruthBins", Form("No bins found for variable %s", var.c_str() ) );
    }

    return binsLowEdge;


}//End of GetTruthBins

std::string NukeCC_Binning::GetXaxisTitle( const std::string& var ) const
{
    if( var == "Emu" )     return "Emu(GeV)";
    else if( var == "pZmu" )     return "Muon Longitudinal Momemtum(GeV/c)";
    else if( var == "Ehad" )    return "Ehad (GeV)";
    else if( var == "Eavail") return "Energy available (GeV)";
    else if( var == "q3") return "q3 (GeV)";
    else if( var == "CCQE-Recoil" )   return "CCQE-style Recoil (GeV)";
    else if( var == "VtxEnergy" )     return "Vertex Energy (GeV)";
    else if( var == "Enu" )     return "Enu (GeV)";
    else if( var == "Q2" )      return "Q^{2} (GeV/c)^{2}";
    else if( var == "W" )      return "W (GeV/c^{2})";
    else if( var == "y" )       return "y";
    else if( var == "x" )       return " x";
    else if( var == "PhiMu" )   return "Muon Phi wrt Beam (deg)";
    else if( var == "CosThetaMu" ) return "Cos(Muon Theta wrt Beam)";
    else if( var == "ThetaMu" ) return "Muon Theta wrt Beam (deg)";
    else if( var == "muonPt" ) return "Muon Pt(GeV)";
    else if( var == "ThetaMuX" )return "Muon Theta_{X} wrt Beam (deg)";
    else if( var == "ThetaMuY" )return "Muon Theta_{Y} wrt Beam (deg)";
    else if( var == "Emum" )    return "Muon Energy at MINOS (GeV)";
    else if( var == "ETheta" )  return "E_{\\mu}*(1-cos(\\theta_{\\mu}))";

    else if( string::npos != var.find("vtxt") )   return "Vertex Time (ns)";
    else if( string::npos != var.find("vtxx") )   return "Vertex X (cm)";
    else if( string::npos != var.find("vtxy") )   return "Vertex Y (cm)";
    else if( string::npos != var.find("vtxz") )   return "Vertex Z (cm)";
    else if( string::npos != var.find("vtxu") )   return "Vertex U (cm)";
    else if( string::npos != var.find("vtxd") )   return "Vertex D (cm)";
    else if( string::npos != var.find("vtxc") )   return "Vertex C (cm)";

    else if( var == "plane" ) return "Plane Number Track-Based";
    else if( var == "planeDNN" )   return "Plane Number DNN";
    else if( var == "NTracks" )   return "Number of Tracks";
    else if( var == "segment" )   return "Segment Number DNN";

    cout << "ERROR [CVUniverse::GetXaxisTitle] : Could not get title for variable: " << var << endl;
    throw 1;
}

std::string NukeCC_Binning::GetYaxisTitle( const std::string& var ) const
{
    double normWidth = GetVarNormWidth( var );
    char words[256];
    if( var == "Emu" )          sprintf(words, "N Events / %.1f GeV", normWidth );
    else if( var == "pZmu" )          sprintf(words, "N Events / %.1f GeV", normWidth );
    else if( var == "Ehad" )    sprintf(words, "N Events / %.1f GeV", normWidth );
    else if( var == "CCQE-Recoil" )    sprintf(words, "N Events / %.1f GeV", normWidth );
    else if( var == "VtxEnergy" )      sprintf(words, "N Events / %.1f GeV", normWidth );
    else if( var == "Enu" )     sprintf(words, "N Events / %.1f GeV", normWidth );
    else if( var == "muonPt" )     sprintf(words, "N Events / %.1f GeV", normWidth );
    else if( var == "Emum" )    sprintf(words, "N Events / %.1f GeV", normWidth );
    else if( var == "ETheta" )  sprintf(words, "N Events / %.01f MeV", normWidth*1000. );
    else if( var == "Q2" )      sprintf(words, "N Events / %.1f (GeV/c)^{2}", normWidth );
    else if( var == "W" )       sprintf(words, "N Events / %.1f (GeV/c^{2})", normWidth );
    //else if( var == "x" )       sprintf(words, "N Events / %.0f ", normWidth ); //norm 1?
    else if( var == "x" )       sprintf(words, "N Events / %.1f ", normWidth );
    else if( var == "y" )       sprintf(words, "N Events / %.2f ", normWidth );

    else if( var == "PhiMu" )   sprintf(words, "N Event / %.0f deg", normWidth );

    else if( var == "ThetaMu"
            ||   var == "CosThetaMu"
            ||   var == "ThetaMuX"
            ||   var == "ThetaMuY" ) sprintf(words, "N Events / %.1f deg", normWidth );

    else if( string::npos != var.find("vtxt") ) sprintf(words, "N Events / %.0f ns", normWidth );
    else if( string::npos != var.find("vtxx") ) sprintf(words, "N Events / %.0f cm", normWidth );
    else if( string::npos != var.find("vtxy") ) sprintf(words, "N Events / %.0f cm", normWidth );
    else if( string::npos != var.find("vtxz") ) sprintf(words, "N Events / %.1f cm", normWidth );
    else if( string::npos != var.find("vtxu") ) sprintf(words, "N Events / %.0f cm", normWidth );
    else if( string::npos != var.find("vtxd") ) sprintf(words, "N Events / %.0f cm", normWidth );
    else if( string::npos != var.find("vtxc") ) sprintf(words, "N Events / %.0f cm", normWidth );

    else if( var == "plane" || var == "planeDNN" ) sprintf(words, "N Events" );
    else if( var == "NTracks") sprintf(words, "N Events" );
    else if( var == "segment") sprintf(words, "N Events" );

    else
    {
        cout << "ERROR [CVUnivers::GetYaxisTitle] : Could not get title for variable: " << var << endl;
        throw 1;
    }

    return string(words);
}


double NukeCC_Binning::GetVarNormWidth( const std::string& var ) const
{
    if( var == "Emu" )    return 1.0;
    if( var == "pZmu" )    return 1.0;
    if( var == "Ehad" )   return .4;
    if( var == "CCQE-Recoil" )   return .25;
    if( var == "VtxEnergy" )     return .25;
    if( var == "Enu" )    return 1.0;
    if( var == "muonPt" )    return 1.0;
    if( var == "Emum" )   return 1.0;
    if( var == "ETheta" ) return .02;
    if( var == "Q2" )     return .1;
    if( var == "W" )     return 1.0;
    //if( var == "y" )      return .1;//?
    if( var == "y" )      return 1;//?
    //if( var == "x" )      return 1.;//?//if norm 1
    if( var == "x" )      return .1;//?

    if( var == "PhiMu" )  return 6.;
    if( var == "ThetaMu" || var == "ThetaMuX" || var == "ThetaMuY" ) return  1.;
    if( var == "CosThetaMu" )  return 1.;

    if( string::npos != var.find("vtxt") )      return 100.;
    if( string::npos != var.find("vtxx") )      return 10.;
    if( string::npos != var.find("vtxy") )      return 10.;
    if( string::npos != var.find("vtxz") )      return 1.7; //if it uses module/plane bins then you need to set the norm width to -1 yourself
    if( string::npos != var.find("vtxd") )      return 10.;
    if( string::npos != var.find("vtxu") )      return 10.;
    if( string::npos != var.find("vtxc") )      return 10.;

    if( var == "plane" || var == "planeDNN"  ) return 1.;
    if( var == "NTracks" || var == "segment" )   return 1.;


    Error("GetVarNormWidth", Form("Could not get the normalization width for variable: %s", var.c_str() ) );
    throw 1;
}

double NukeCC_Binning::GetDISVarMinVal( const std::string& var ) const
{
    if( var == "Emu" )    return 2.0;
    if( var == "pZmu" )    return 2.0;
    if( var == "Ehad" )   return 1.5;
    if( var == "CCQE-Recoil" )   return 0.;
    if( var == "VtxEnergy" )     return 0.;
    if( var == "Enu" )    return 5.0;
    if( var == "Emum" )   return 2.0;
    if( var == "ETheta" ) return 0.;
    if( var == "Q2" )     return 1.0;
    if( var == "W" )      return 2.0;
    if( var == "y" )      return 0.;
    if( var == "x" )      return 0.005;

    if( var == "PhiMu" ) return -180.;

    if( var == "ThetaMu" ) return 0.;
    if( var == "ThetaMuX" || var == "ThetaMuY" ) return  -20.;

    if( var == "plane" || var == "planeDNN"  ) return 0;
    if( var == "NTracks" ) return 0;
    if( var == "segment" ) return 0.;

    if( var == "vtxx" || var == "vtxy" ) return -100.;
    if( var == "vtxz" ) return 400.0;

    Error("GetVarMinVal", Form("Could not get the min value to plot for variable: %s", var.c_str() ) );
    throw 1;
}
double NukeCC_Binning::GetDISVarMaxVal( const std::string& var ) const
{
    const double e = EPSILON;
    if( var == "Emu" )    return 50.-e;
    if( var == "pZmu" )    return 40.-e;
    if( var == "Ehad" )   return 20.-e;
    if( var == "CCQE-Recoil" )   return 3.-e;
    if( var == "VtxEnergy" )     return 3.-e;
    if( var == "Enu" )    return 50.-e;
    if( var == "Emum" )   return 50.-e;
    if( var == "ETheta" ) return .5-e;
    if( var == "Q2" )     return 8.-e;
    if( var == "W" )      return 8.-e;
    if( var == "y" )      return 1.-e;
    if( var == "x" )      return 0.8 - e;

    if( var == "PhiMu" ) return 180.-e;
    if( var == "plane" || var == "planeDNN"  ) return 68-e;
    if( var == "NTracks" ) return 10-e;
    if( var == "segment" ) return 174.;

    if( var == "ThetaMu" ) return 17.-e; //was 24
    if( var == "ThetaMuX" || var == "ThetaMuY" ) return  20.-e;

    if( var == "vtxx" || var == "vtxy" ) return 100.-e;
    if( var == "vtxz" ) return 850.-e;

    Error("GetVarMaxVal", Form("Could not get the max value to plot for variable: %s", var.c_str() ) );
    throw 1;
}
#endif
