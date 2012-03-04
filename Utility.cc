#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <fstream>
#include "Option.h"
#include "Utility.h"
#include "FileReader.h"

// (C) Yasuda, Tomohiro: The university of Tokyo

void trap_point(){}

BITVECTOR_LIB_BEGIN

debug_values::debug_values(){
  for(uint32_t i=0;i<nelems(m_debug_values);++i) m_debug_values[i]=0;
}
void debug_values::Show(std::ostream& stream, int32_t indent){
  for(uint32_t i=0;i<nelems(m_debug_values);++i){
    if((i%10)==0){
      if(i!=0) stream<<std::endl;
      stream<<i<<":";
    }
    stream<<"\t"<<m_debug_values[i];
  }
  stream<<std::endl;
}
debug_values g_debug_values;


int64_t str2size(const char *nptr){
  char *p=0;
  int64_t v=strtoll(nptr,&p,10);
  if(nptr==p) Quit("Number expected: "<<nptr);
  switch(*p){
  case '\0':
    break;
  case 'g':
  case 'G':
    v *= 1024;
  case 'm':
  case 'M':
    v *= 1024;
  case 'k':
  case 'K':
    v *= 1024;
    break;
  default:
    Quit("Extra characters after an number: "<<nptr);
  }
  return v;
}

/*
void ShowIndent(std::ostream &stream, int32_t indent){
  if(indent==0) return;
  stream<<std::endl;
  for(int32_t i=0;i<indent;++i) stream<<' ';
}
*/
std::ostream& ShowIndent(std::ostream& stream, int32_t indent){
  if(indent<0) return stream;
  stream<<std::endl;
  for(int i=0;i<indent;++i) stream<<' ';
  return stream;
}

///////////////////////////////////////////////////////////////////////
// Write an integer in hexadecimal representation
Hex::Hex(int64_t v){
  std::sprintf(buffer, "0x%016llx", static_cast<long long unsigned>(v));
};
std::ostream& operator<<(std::ostream &stream, const Hex &h){
  stream<<h.buffer;
  return stream;
}

///////////////////////////////////////////////////////////////////////
// Write an integer in hexadecimal representation
Hex32::Hex32(int64_t v){
  std::sprintf(buffer, "0x%08llx", static_cast<long long unsigned>(v));
};
std::ostream& operator<<(std::ostream &stream, const Hex32 &h){
  stream<<h.buffer;
  return stream;
}

///////////////////////////////////////////////////////////////////////
// Process memory

int64_t CProcessMemory::get_size(const char *p){
  while(! isdigit(*p)) ++p;
  return atoll(p);
}

CProcessMemory::CProcessMemory(){
  Size=Peak=Data=0;
  std::ostringstream oss;
  oss<<"/proc/"<<getpid()<<"/status";
  std::ifstream file(oss.str().c_str());
  if(!file) Quit("Cannot open "<<oss.str());
  std::string line;
  while( getline(file,line) ){
    if( line.compare(0,6,"VmPeak")==0 ){Peak=get_size(line.c_str()); continue;}
    if( line.compare(0,6,"VmSize")==0 ){Size=get_size(line.c_str()); continue;}
    if( line.compare(0,6,"VmData")==0 ){Data=get_size(line.c_str()); continue;}
  }
  //std::cerr<<"Peak="<<Peak<<std::endl;
  //std::cerr<<"Size="<<Size<<std::endl;
  //std::cerr<<"Data="<<Data<<std::endl;
}

void CProcessMemory::Write(std::ostream &stream, int32_t indent) const {
  stream<<"{\"peak\":"<<Peak<<",\"size\":"<<Size<<",\"data\":"<<Data<<"}";
}

///////////////////////////////////////////////////////////////////////
// Mapping from strings to integers
void str2int_t::Write(std::ostream &stream, int32_t indent) const{
  stream<<'{';
  foreach_const(str2int_t,iter,*this){
    if(iter!=begin()) stream<<",";
    stream<<"\""<<(*iter).first<<"\":"<<(*iter).second;
  }
  stream<<'}';
}

// Mapping from strings to string
void str2str_t::Write(std::ostream &stream, int32_t indent) const{
  stream<<'{';
  foreach_const(str2str_t,iter,*this){
    if(iter!=begin()) stream<<",";
    stream<<"\""<<(*iter).first<<"\":\""<<(*iter).second<<"\"";
  }
  stream<<'}';
}


///////////////////////////////////////////////////////////////////////
CCountLines::CCountLines(const char *t){
  m_tag=t;
  m_n_all=m_n_chr=0;
  m_is_valid=Option().Find("verbose");
}

///////////////////////////////////////////////////////////////////////
bool check_prefix(const char *prefix, const char *str){
  Assert(prefix);
  Assert(*prefix);
  Assert(str);
  // Count the length of the prefix
  int32_t len=0;
  while(prefix[len]) ++len;
  // Compare
  for(int32_t i=0;i<len;++i){
    if(prefix[i]!=str[i]) return false;
  }
  return true;
}


BITVECTOR_LIB_END


////////////////////////////////////////////////////////////////////////////////
// Key-value store
/*
bool CKVStore::Read(const char *target_file){
  CFileReader fr(target_file);
  fr.SetProgressInterval(Option().RequireInteger("progress-interval"),"lines");
  for(;;){
    const char *p=fr.GetContentLine();
    if(!p) break;
    Parse(&p);
  }
  return true;
}

bool CKVStore::Add(const char *key, const char *value){
  mapping_t::const_iterator iter = m_map.find(key);
  if(iter!=m_map.end()){
    if((*iter).second==value) return true;
    Quit(key<<" is associated with "<<(*iter).second<<", conflicting with "<<value);
    //return false;
  }
  m_map[key]=value;
  return true;
}

bool CKVStore::Parse(const char **pp){
  CTokenizer td(*pp,"\t");
  //const char *p=*pp;
  std::string key,value;
  key   = td.NextString();
  value = td.NextString();
  //std::cerr<<"CKVStore::Parse: key='"<<key<<"', value='"<<value<<"', p='"<<p<<"'\n";
  //m_map[key]=value;
  Add(key.c_str(),value.c_str());
  *pp=td.Peep();
  return true;
}

void CKVStore::Write(std::ostream &stream, int32_t indent) const {
  foreach_const(mapping_t,iter,m_map){
    stream<<(iter==m_map.begin()? "{": ",\n ");
    stream<<"\""<<(*iter).first<<"\": \""<<(*iter).second<<"\"";
  }
  stream<<'}';
}

const char *CKVStore::Find(const char *key) const {
  mapping_t::const_iterator iter = m_map.find(key);
  if(iter==m_map.end()) return 0;
  return (*iter).second.c_str();
}
*/
