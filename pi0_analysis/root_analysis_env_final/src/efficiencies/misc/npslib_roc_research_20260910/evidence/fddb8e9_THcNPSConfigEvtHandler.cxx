#include "THcNPSConfigEvtHandler.h"
#include "THaGlobals.h"
#include "THcGlobals.h"
#include "THaRunBase.h"
#include "THaEvData.h"
#include "DAQconfig.h"
#include "THcParmList.h"
#include <sstream>

using namespace std;

/////////////////////////////////////////////////////////////////////
//
// Parse NPS FADC250 and VTP config data (ev 137) -- regular text
// 
/////////////////////////////////////////////////////////////////////

//_______________________________________________________________
THcNPSConfigEvtHandler::THcNPSConfigEvtHandler( const char* name,
						const char* description ) :
  THaEvtTypeHandler(name, description)
{
}

//_______________________________________________________________
THcNPSConfigEvtHandler::~THcNPSConfigEvtHandler()
{

  for( auto& cfg : fNPSConfigList ) {
    gHcParms->RemoveString(cfg.par);
  }
}

//_______________________________________________________________
THaAnalysisObject::EStatus THcNPSConfigEvtHandler::Init( const TDatime& date )
{
  // default event type for nps config
  if( fEvtTypes.empty() ) {
    fEvtTypes.push_back(137);
  }

  // default parameter list
  fNPSConfigList.push_back( {Form("g%s_fadc250_sparsification", GetName()),"FADC250_SPARSIFICATION"} );
  fNPSConfigList.push_back( {Form("g%s_vtp_clus_trig_thr",      GetName()),"VTP_NPS_ECALCLUSTER_CLUSTER_TRIGGER_THR"} );
  fNPSConfigList.push_back( {Form("g%s_vtp_clus_ro_thr",        GetName()),"VTP_NPS_ECALCLUSTER_CLUSTER_READOUT_THR"} );
  fNPSConfigList.push_back( {Form("g%s_vtp_pair_trig_thr",      GetName()),"VTP_NPS_ECALCLUSTER_CLUSTER_PAIR_TRIGGER_THR"} );      

  return THaEvtTypeHandler::Init(date);
}

//_______________________________________________________________
void THcNPSConfigEvtHandler::AddEvtType( UInt_t evtype )
{
  // We don't want to add this event type to the evt type list of THaEvtTypeHandler
  // eventtypes from THaEvtTypeHandler is looked up by all inherited EvtTypeHandler classes
  // Instead, we set the event types only relevant for this class
  
  if( std::find(fEvtTypes.begin(), fEvtTypes.end(), evtype ) == fEvtTypes.end() )
    fEvtTypes.push_back(evtype);
}

//_______________________________________________________________
void THcNPSConfigEvtHandler::AddParameter(std::string parname, std::string keyname)
{
  fNPSConfigList.push_back( {Form("g%s_%s", GetName(), parname.data()), keyname} );
}
 
//_______________________________________________________________
void THcNPSConfigEvtHandler::MakeParms()
{
  // all data parsed as a single string varaible 

  for( auto& cfg : fNPSConfigList ) {
    if( !GetInfo(cfg.key.data()).empty() ) {
      gHcParms->RemoveString(cfg.par);
      gHcParms->AddString(cfg.par, GetInfo(cfg.key.data()) );
    }
  }      

}

//_______________________________________________________________
const string& THcNPSConfigEvtHandler::GetInfo(const char* parname )
{
  static const string nullstr;
  auto it = fNPSConfigData.find(parname);
  if( it != fNPSConfigData.end() )
    return it->second;
  else
    return nullstr;
}

//_______________________________________________________________
Int_t THcNPSConfigEvtHandler::Analyze( THaEvData* evdata )
{

  //std::cout << "********* THcNPSConfigEvtHandler::Analyze *********" << std::endl;

  UInt_t evtype = evdata->GetEvType();

  // Check event type
  if( std::find(fEvtTypes.begin(), fEvtTypes.end(), evtype) == fEvtTypes.end() )
    return -1;

  auto* ifo = DAQInfoExtra::GetFrom(evdata->GetExtra());
  if( !ifo ) return -1;

  // Parse string and store (key, value) pairs into a new container
  // VTP_* keys are missing when using THaRunBase::GetDAQInfo
  // So we parse them separately here

  for( auto& this_info : ifo->strings ) {
    istringstream ifstr(this_info);
    string line;
    while( getline(ifstr, line) ) {
      // skip blank lines
      if( line.find_first_not_of(" \t") == string::npos ) 
	continue;

      auto items = Podd::vsplit(line);
      if( !items.empty() ) {
	string& key = items[0];
	string val;
	val.reserve(line.size());
	for( size_t j = 1, e = items.size(); j < e; ++j ) {
	  val.append(items[j]);
	  if( j + 1 != e )
	    val.append(" ");
	}
	if( val != "end" )
	  fNPSConfigData.emplace(std::move(key), std::move(val));
      }
    }// getline
  }

  // Add to the gHcParm list
  MakeParms();

  /*
  for( auto& keyval : ifo->keyval )
    std::cout << keyval.first << "\t" << keyval.second <<  std::endl;
  std::cout << "Print FADC Mode: " << gHaRun->GetDAQInfo("FADC250_SPARSIFICATION") << std::endl;
  std::cout << "Print VTP Parameters: " << gHaRun->GetDAQInfo("VTP_FIRMWARETYPE") << std::endl;
  std::cout << "Print VTP Parameters: " << gHaRun->GetDAQInfo("VTP_NPS_ECALCLUSTER_CLUSTER_TRIGGER_THR") << std::endl;
  */

  return 0;
}

ClassImp(THcNPSConfigEvtHandler)
