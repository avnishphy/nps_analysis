#ifndef ROOT_THcNPSCoinTime
#define ROOT_THcNPSCoinTime

//////////////////////////////////////////////
//
// THcNPSCoinTime
//
// Calculate coincidence times for a pair of
// NPS + HMS spectrometer. 
// timing offsets for different paddles are handled
// in the detector class (e.g. THcNPSArray)
//
//////////////////////////////////////////////

#include "TString.h"
#include "THaPhysicsModule.h"
#include "THcTrigDet.h"
#include "THcHallCSpectrometer.h"
#include "THaTrack.h"
#include "THcNPSApparatus.h"
#include "THcNPSCalorimeter.h"

#define LIGHTSPEED 29.9792         // cm/ns
#define MASS_ELECTRON 0.000510998  // electron mass in GeV/c^2
#define MASS_PROTON 0.93827208     // proton mass in GeV/c^2

class THcNPSCoinTime : public THaPhysicsModule {
 public:
  THcNPSCoinTime( const char* name, const char* description,
		  const char* elecArmName="", const char* hadArmName="",
		  const char* coinname="");

  virtual ~THcNPSCoinTime();
    
  virtual EStatus Init( const TDatime& run_time ) ;
  virtual Int_t   Process( const THaEvData& );

  void            Clear( Option_t* opt="" );
  void            Reset( Option_t* opt = "" );

 protected:

  virtual Int_t ReadDatabase( const TDatime& date );
  virtual Int_t DefineVariables( EMode mode = kDefine );

  TString               fCoinDetName;
  THcTrigDet           *fCoinDet;
  TString               felecArmName;
  TString               fhadArmName;
  THcHallCSpectrometer *fHMSSpect;    // HMS
  THcNPSApparatus      *fNPSSpect;    // NPS
  THcNPSCalorimeter    *fNPSCalo;

  THaTrack             *theHMSTrack;

  // Read from Database
  Double_t fNPSAngle;
  Double_t feHad_CT_Offset;
  Double_t fHMScentralPathLen;
  Double_t fNPScentralPathLen;

  // Calculated using HMS track
  Double_t fDeltaHMSpathLength;

  Double_t fROC1_RAW_CoinTime1;
  Double_t fROC1_RAW_CoinTime2;
  Double_t fROC1_RAW_CoinTime1_NoTrack;
  Double_t fROC1_RAW_CoinTime2_NoTrack;

  Double_t fNPS_RAW_CoinTime;
  Double_t fHMS_RAW_CoinTime;

  Double_t fROC1_epCoinTime1;
  Double_t fROC1_epCoinTime2;
  Double_t fNPS_epCoinTime;
  Double_t fHMS_epCoinTime;

  // Trigger time pTrig1 T1 (NPS VTP||EDTM) and T3 (HMS 3/4)
  Double_t pNPS_TdcTime_ROC1;
  Double_t pHMS_TdcTime_ROC1;
  Double_t pELRE_TdcTime_ROC1;

  Double_t elec_coinCorr;
  Double_t had_coinCorr_proton;

  ClassDef(THcNPSCoinTime,0)
};

#endif
