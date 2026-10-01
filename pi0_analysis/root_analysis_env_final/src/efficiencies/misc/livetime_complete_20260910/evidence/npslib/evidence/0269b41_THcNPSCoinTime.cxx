#include "THcNPSCoinTime.h"
#include "THcGlobals.h"
#include "THcParmList.h"
#include "THaEvData.h"
#include "THcNPSCluster.h"
#include <iostream>

using namespace std;

//________________________________________________________________________
THcNPSCoinTime::THcNPSCoinTime( const char* name, const char* description,
				const char* elecArmName, const char* hadArmName,
				const char* coinname) :
  THaPhysicsModule(name, description),
  fCoinDetName(coinname),
  fCoinDet(nullptr),
  felecArmName(elecArmName),
  fhadArmName(hadArmName),
  fHMSSpect(nullptr),
  fNPSSpect(nullptr)
{
  // Constructor  
}
  
//________________________________________________________________________
THcNPSCoinTime::~THcNPSCoinTime()
{
  // Destructor
  RemoveVariables();
}

//________________________________________________________________________
void THcNPSCoinTime::Clear( Option_t* opt )
{
  fROC1_RAW_CoinTime1 = fROC1_RAW_CoinTime2 = fNPS_RAW_CoinTime = fHMS_RAW_CoinTime = kBig;
  fROC1_RAW_CoinTime1_NoTrack = fROC1_RAW_CoinTime2_NoTrack = kBig;
  fROC1_epCoinTime1 = fROC1_epCoinTime2 =  kBig;
  fNPS_epCoinTime = fHMS_epCoinTime = kBig;
  elec_coinCorr = had_coinCorr_proton = kBig;
}

//________________________________________________________________________
void THcNPSCoinTime::Reset( Option_t* opt )
{
  Clear(opt);
}

//________________________________________________________________________
THaAnalysisObject::EStatus THcNPSCoinTime::Init( const TDatime& run_time )
{

  // Locate the spectrometer appratus and save pointer to it
  cout << "*************************************************" << endl;
  cout << "Initializing THcNPSCointTime Physics Modue" << endl;
  cout << "Hadron Arm   -------> " << fhadArmName << endl;
  cout << "Electron Arm -------> " << felecArmName << endl;
  cout << "TrigDet  -------> " << fCoinDetName << endl;
  cout << "**************************************************" << endl;

  fStatus = kOK;

  // Get spectrometers
  if (felecArmName == "H") {
    // HMS
    fHMSSpect = dynamic_cast<THcHallCSpectrometer*>
      ( FindModule( felecArmName.Data(), "THcHallCSpectrometer"));
    // NPS
    fNPSSpect = dynamic_cast<THcNPSApparatus*>
      ( FindModule( fhadArmName.Data(), "THcNPSApparatus"));
  }
  else {
    fHMSSpect = dynamic_cast<THcHallCSpectrometer*>
      ( FindModule( fhadArmName.Data(), "THcHallCSpectrometer"));
    
    fNPSSpect = dynamic_cast<THcNPSApparatus*>
      ( FindModule( felecArmName.Data(), "THcNPSApparatus"));
  }

  if( !fHMSSpect || !fNPSSpect ) {
    fStatus = kInitError;
    return fStatus;
  }

  fCoinDet = dynamic_cast<THcTrigDet*>
    ( FindModule( fCoinDetName.Data(), "THcTrigDet"));
  if( !fCoinDet ) {
    cout << "THcCoinTime module  Cannnot find TrigDet = " << fCoinDetName.Data() << endl;
    fStatus = kInitError;
    return fStatus;
  }
  
  // Get NPS calorimeter
  fNPSCalo = dynamic_cast<THcNPSCalorimeter*>(fNPSSpect->GetDetector("cal"));
  if( !fNPSCalo ) {
    fStatus = kInitError;
    return fStatus;
  }

  // NPS angle   
  fNPSAngle = TMath::DegToRad()*fNPSSpect->GetNPSAngle();
  
  if ((fStatus=THaPhysicsModule::Init( run_time )) != kOK )
    return fStatus;

  return fStatus;
}

//________________________________________________________________________
Int_t THcNPSCoinTime::ReadDatabase( const TDatime& date )
{
  
  DBRequest list[] = {
    {"eHadCoinTime_Offset", &feHad_CT_Offset,    kDouble, 0, 1}, // coin time offset for ep coincidence
    {"HMS_CentralPathLen",  &fHMScentralPathLen, kDouble, 0, 1},
    {"NPS_CentralPathLen",  &fNPScentralPathLen, kDouble, 0, 1},
    {0}
  };

  // Default values if not read from param file
  feHad_CT_Offset = 0.0;
  fHMScentralPathLen = 22.*100;
  fNPScentralPathLen = 9.5*100; 

  gHcParms->LoadParmValues((DBRequest*)&list, "");

  return kOK;
}

//________________________________________________________________________
Int_t THcNPSCoinTime::DefineVariables( EMode mode )
{

  if( mode == kDefine && fIsSetup ) return kOK;
  fIsSetup = ( mode == kDefine );

  RVarDef vars[] = {
    {"epCoinTime1_ROC1",   "ROC1 Corrected ep Coincidence Time, NPS & HMS 3/4", "fROC1_epCoinTime1"},
    {"epCoinTime2_ROC1",   "ROC1 Corrected ep Coincidence Time, NPS & El-Real", "fROC1_epCoinTime2"},
    {"epCoinTime_NPS",     "NPS Corrected ep Coincidence Time",  "fNPS_epCoinTime"},
    {"epCoinTime_HMS",     "HMS Corrected ep Coincidence Time",  "fHMS_epCoinTime"},
    
    {"CoinTime1_RAW_ROC1", "ROC1 RAW Coincidence Time, NPS & HMS 3/4",     "fROC1_RAW_CoinTime1"},
    {"CoinTime2_RAW_ROC1", "ROC1 RAW Coincidence Time, NPS & El-Real",     "fROC1_RAW_CoinTime2"},
    {"CoinTime1_RAW_ROC1_NoTrack", "ROC1 RAW Coincidence Time w/o Track Param, NPS & HMS 3/4",     "fROC1_RAW_CoinTime1_NoTrack"},
    {"CoinTime2_RAW_ROC1_NoTrack", "ROC1 RAW Coincidence Time w/o Track Paarm, NPS & El-Real",     "fROC1_RAW_CoinTime2_NoTrack"},
    {"CoinTime_RAW_NPS",   "NPS RAW Coincidence Time",                     "fNPS_RAW_CoinTime"},
    {"CoinTime_RAW_HMS",   "HMS RAW Coincidence Time",                     "fHMS_RAW_CoinTime"},

    {"DeltaHMSPathLength", "DeltaHMSpathLength (cm)",  "fDeltaHMSpathLength"},

    {"elec_coinCorr",      "",                         "elec_coinCorr"},
    {"had_coinCorr_proton","",                         "had_coinCorr_proton"},
    {0}
  };

  return DefineVarsFromList( vars, mode );
}

//________________________________________________________________________
Int_t THcNPSCoinTime::Process( const THaEvData& evdata )
{

  if( !IsOK() ) return -1;

  // Raw Tdctime
  // We only use ROC1
  pNPS_TdcTime_ROC1 = fCoinDet->Get_CT_Trigtime(0);   // pTRIG1_ROC1: NPS VTP
  pHMS_TdcTime_ROC1 = fCoinDet->Get_CT_Trigtime(1);   // pTRIG3_ROC1: HMS 3/4
  pELRE_TdcTime_ROC1 = fCoinDet->Get_CT_Trigtime(3);  // pTRIG4_ROC1: HMS ELREAL
  
  // Raw, Uncorrected Coincidence Time with out Track
  fROC1_RAW_CoinTime1_NoTrack =  pNPS_TdcTime_ROC1 - pHMS_TdcTime_ROC1;
  fROC1_RAW_CoinTime2_NoTrack =  pNPS_TdcTime_ROC1 - pELRE_TdcTime_ROC1;

  // HMS Track, NPS cluster
  // Get HMS tracks
  THaTrackInfo* hms_trkifo = fHMSSpect->GetTrackInfo();
  if( !hms_trkifo ) return 1;

  // Get NPS clusters
  if( fNPSCalo->GetNClusters() == 0 ) return 1;

  theHMSTrack = fHMSSpect->GetGoldenTrack();
  if( !theHMSTrack )
    return 1;

  // Loop over all clusters
  // We don't have a "Golden cluster", use highest energy cluster for now
  Double_t ClusterMaxE = 0.;
  Double_t NPS_FPtime = kBig;
  TVector3 pvect;
  TVector3 vertex = theHMSTrack->GetVertex();  // Use vertex from HMS
  for(auto& cluster : fNPSCalo->GetClusters()) {
    if( cluster.E() > ClusterMaxE ) {
      ClusterMaxE = cluster.E();

      // Get P vector in lab frame
      cluster.RotateToLab(fNPSAngle, vertex, pvect);      

      // Cluster time
      NPS_FPtime = cluster.GetTime(); 
    }// if found higher energy cluster
  }

  // HMS
  Double_t hms_xptar = theHMSTrack->GetTTheta();
  Double_t hms_dP = theHMSTrack->GetDp();
  Double_t hms_xfp = theHMSTrack->GetX();
  Double_t hms_xpfp = theHMSTrack->GetTheta();
  Double_t hms_ypfp = theHMSTrack->GetPhi();
  Double_t HMS_FPtime = theHMSTrack->GetFPTime();
  
  // Assume these are the values when the variable is initialized
  if (NPS_FPtime == -2000 || HMS_FPtime == -2000) return 1;
  if (NPS_FPtime == -1000 || HMS_FPtime == -1000) return 1;

  fDeltaHMSpathLength = -1.0*(12.462*hms_xpfp + 0.1138*hms_xpfp*hms_xfp - 0.0154*hms_xfp - 72.292*hms_xpfp*hms_xpfp - 0.0000544*hms_xfp*hms_xfp - 116.52*hms_ypfp*hms_ypfp);
  fDeltaHMSpathLength = (.12*hms_xptar*1000 +0.17*hms_dP/100.);
  
  // First assume NPS is the electron arm
  Double_t ElecPathLength = fNPScentralPathLen;
  Double_t HadPathLength  = fHMScentralPathLen + fDeltaHMSpathLength;
  Double_t elec_P = pvect.Mag();
  Double_t had_P = theHMSTrack->GetP();
  Int_t sign= -1;
  if(felecArmName == "H") {
    ElecPathLength = fHMScentralPathLen + fDeltaHMSpathLength;
    HadPathLength  = fNPScentralPathLen;    // NPS dP = 0
    elec_P = theHMSTrack->GetP();                     // HMS golden track momentum
    had_P = pvect.Mag();                             // NPS golden cluster momentum
    sign = 1;
  }
  
  // beta calculations beta = v/c = p/E
  Double_t elecArm_BetaCalc = elec_P / sqrt(elec_P*elec_P + MASS_ELECTRON*MASS_ELECTRON);
  Double_t hadArm_BetaCalc_proton = had_P / sqrt(had_P*had_P + MASS_PROTON*MASS_PROTON);

  //Coincidence Corrections
  elec_coinCorr = (ElecPathLength) / (LIGHTSPEED * elecArm_BetaCalc );
  had_coinCorr_proton = (HadPathLength) / (LIGHTSPEED * hadArm_BetaCalc_proton );

  // Raw, Uncorrected Coincidence Time
  // Since we have two coincidence trig, check both

  fROC1_RAW_CoinTime1 =  (pNPS_TdcTime_ROC1 + NPS_FPtime) - (pHMS_TdcTime_ROC1 + HMS_FPtime);
  fROC1_RAW_CoinTime2 =  (pNPS_TdcTime_ROC1 + NPS_FPtime) - (pELRE_TdcTime_ROC1 + HMS_FPtime);

  // Not sure if we want to keep those
  fNPS_RAW_CoinTime = NPS_FPtime - HMS_FPtime; // essentially same as fHMS_RAW_CoinTime
  fHMS_RAW_CoinTime = NPS_FPtime - HMS_FPtime;

  // Corrected Coincidence Time for ROC1
  fROC1_epCoinTime1 = fROC1_RAW_CoinTime1 + sign*( elec_coinCorr-had_coinCorr_proton) - feHad_CT_Offset; 
  fROC1_epCoinTime2 = fROC1_RAW_CoinTime2 + sign*( elec_coinCorr-had_coinCorr_proton) - feHad_CT_Offset; 
  fNPS_epCoinTime = fNPS_RAW_CoinTime + sign*( elec_coinCorr-had_coinCorr_proton) - feHad_CT_Offset; // essentially same as fHMS_epCoinTime
  fHMS_epCoinTime = fHMS_RAW_CoinTime + sign*( elec_coinCorr-had_coinCorr_proton) - feHad_CT_Offset; 

  return 0;
}

ClassImp(THcNPSCoinTime)
