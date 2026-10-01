#include "THaCodaFile.h"
#include "CodaDecoder.h"
#include <TSystem.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <map>
#include <vector>
#include <functional>
#include <stdexcept>
#include <cstring>
// Reads ONLY the single user-authorized file. No detector replay or staging.
// EVIO handles byte order; saved bank words are host-endian uint32_t.
void raw_extract(long long max_records=0) {
  const char* path="/cache/hallc/c-nps/raw/nps_coin_4305.dat.0";
  gSystem->mkdir("raw4305/nonphysics",true);gSystem->mkdir("raw4305/scaler_banks",true);
  Decoder::THaCodaFile file(path);if(!file.isOpen())throw std::runtime_error("Cannot open raw file");
  std::ofstream idx("raw4305/index.jsonl"),sc("raw4305/scaler_banks.jsonl"),events("raw4305/physics.bin",std::ios::binary);
  std::map<unsigned,unsigned long long> tags;unsigned long long record=0,physics=0,last=0,nbanks=0;int status=0;
  while((!max_records||record<(unsigned long long)max_records)&&(status=file.codaRead())==CODA_OK) {
    ++record;auto* v=file.getEvBuffer();unsigned len=v[0]+1,tag=v[1]>>16;tags[tag]++;
    bool phys=tag==0xff50||tag==0xff58||tag==0xff70||tag==0xff78;
    if(phys) {
      unsigned block=v[1]&255;if(block!=1)throw std::runtime_error("Unexpected multi-event block");
      Decoder::CodaDecoder::TBOBJ t;t.Fill(v+2,block,1);
      unsigned long long tm=0,bits=0;if(t.evTS)tm=t.evTS[0];else if(t.TSROC){std::memcpy(&tm,t.TSROC,8);tm&=0x0000ffffffffffffULL;}
      if(t.withTriggerBits())bits=t.TSROC[2]&63;
      if(t.evtNum!=last+1)throw std::runtime_error("Physics event-number gap");
      last=t.evtNum;++physics;unsigned long long row[4]={record,last,tm,bits};events.write((char*)row,sizeof(row));
    } else {
      std::string name="nonphysics/record"+std::to_string(record)+"_tag"+std::to_string(tag)+".bin";
      std::ofstream out("raw4305/"+name,std::ios::binary);out.write((char*)v,len*4);
      idx<<"{\"record\":"<<record<<",\"last_physics\":"<<last<<",\"tag\":"<<tag<<",\"words\":"<<len<<",\"file\":\""<<name<<"\"}\n";
      if(record<20)std::cout<<"RAW_USER "<<record<<" tag "<<tag<<" words "<<len<<std::endl;
    }
    std::function<void(const unsigned*,unsigned,std::string)> walk;
    walk=[&](const unsigned* p,unsigned available,std::string route) {
      if(available<2||p[0]+1>available||p[0]<1)throw std::runtime_error("Invalid bank length");
      unsigned n=p[0]+1,bt=p[1]>>16,dt=(p[1]>>8)&63;route+="/"+std::to_string(bt);
      if(dt==0x10||dt==0x0e) {
        unsigned pos=2;while(pos<n){if(p[pos]+1>n-pos)throw std::runtime_error("Child bank exceeds parent");walk(p+pos,n-pos,route);pos+=p[pos]+1;}
      } else if(dt==1 && n>2 && (bt==3801 || p[2]==0x00000620 || p[2]==0x00200720 || p[2]==0x00400820 || p[2]==0x00600920 || p[2]==0x00800a20 || p[2]==0x00a00b20 || p[2]==0x00c00c20)) {
        std::string name="scaler_banks/bank"+std::to_string(++nbanks)+".bin";std::ofstream out("raw4305/"+name,std::ios::binary);out.write((char*)p,n*4);
        sc<<"{\"record\":"<<record<<",\"last_physics\":"<<last<<",\"top_tag\":"<<tag<<",\"route\":\""<<route<<"\",\"words\":"<<n<<",\"file\":\""<<name<<"\"}\n";
      }
    };
    walk(v,len,"");
    if(physics&&physics%100000==0)std::cout<<"RAW_PROGRESS "<<record<<" physics "<<physics<<" banks "<<nbanks<<std::endl;
  }
  if(!max_records&&status!=CODA_EOF)throw std::runtime_error("Raw read did not end at EOF");
  std::ofstream summary("raw4305/extraction_summary.json");summary<<"{\"records\":"<<record<<",\"physics_events\":"<<physics<<",\"last_physics\":"<<last<<",\"scaler_banks\":"<<nbanks<<",\"EOF\":"<<(status==CODA_EOF?"true":"false")<<",\"tags\":{";
  bool first=true;for(auto x:tags){if(!first)summary<<",";first=false;summary<<"\""<<x.first<<"\":"<<x.second;}summary<<"}}\n";
  std::cout<<"RAW_COMPLETE records "<<record<<" physics "<<physics<<" scaler banks "<<nbanks<<std::endl;
}
