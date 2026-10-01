#ifndef ROOT_THcNPSConfigEvtHandler
#define ROOT_THcNPSConfigEvtHandler

#include "THaEvtTypeHandler.h"
#include <vector>
#include <string>
#include <map>

class THcNPSConfigEvtHandler : public THaEvtTypeHandler {
 public:

  THcNPSConfigEvtHandler( const char* name, const char* description = "" );
  virtual ~THcNPSConfigEvtHandler();
  
  virtual Int_t   Analyze( THaEvData *evdata );
  virtual EStatus Init( const TDatime& date );
  virtual void    AddEvtType( UInt_t evtype );
  virtual void    AddParameter(std::string parname, std::string keyname);

  std::vector<UInt_t> GetEvtTypes() { return fEvtTypes; }
  const std::string&  GetInfo( const char* parname );

 private:

  std::vector<UInt_t> fEvtTypes;

  struct NPSConfig {
    std::string par;
    std::string key;
  };
  std::vector<NPSConfig> fNPSConfigList;

  std::map<std::string, std::string> fNPSConfigData; // container for parsed data

  void MakeParms();

  ClassDef(THcNPSConfigEvtHandler,0)
};

#endif
