#ifndef _UTILITY_H_
#define _UTILITY_H_

// (C) Yasuda, Tomohiro: The university of Tokyo

#include <stdint.h>
#include <iostream>
#include <cstdlib> // for Quit
#include <sstream>
#include <map>
#include <string>

#ifndef BITVECTOR_LIB_BEGIN
#define BITVECTOR_LIB_BEGIN namespace BitVectorLib {
#define BITVECTOR_LIB_END   }
#endif

extern void trap_point();
typedef std::string CError;
#define Quit(x) if(1){std::ostringstream s; s<<"Error: "<<__FILE__<<':'<<__LINE__<<": "<<x; trap_point(); std::string ss(s.str()); throw ss;}
#define Assert(x) if(!(x)){std::cerr<<"Assertion: "<<__FILE__<<':'<<__LINE__<<": Assertion "<<__STRING(x)<<" failed."<<std::endl; trap_point(); std::exit(1);}
//#define Assert(x) Quit("Assertion "<<__STRING(x)<<" failed.")
#define Warning(x) if(1){std::cerr<<"Warning: "<<x<<std::endl;}

BITVECTOR_LIB_BEGIN

extern int64_t str2size(const char *nptr);


std::ostream& ShowIndent(std::ostream& stream, int32_t indent=0);
//void ShowIndent(std::ostream &stream, int32_t indent=0);

class Hex {
private:
public:
  char buffer[50];
  Hex(int64_t v);
  //inline Hex(int64_t v){std::sprintf(buffer, "0x%016llx", static_cast<long long unsigned>(v));}
};
std::ostream& operator<<(std::ostream &stream, const Hex &h);

class Hex32 {
private:
public:
  char buffer[50];
  Hex32(int64_t v);
  //inline Hex(int64_t v){std::sprintf(buffer, "0x%016llx", static_cast<long long unsigned>(v));}
};
std::ostream& operator<<(std::ostream &stream, const Hex32 &h);

#define foreach(t,v,c) for(t::iterator v=(c).begin(); v!=(c).end(); ++v)
#define foreach_const(t,v,c) for(t::const_iterator v=(c).begin(); v!=(c).end(); ++v)
//#define nelem(a) (sizeof(a)/sizeof(a[0]))
#define nelems(a) (sizeof(a)/sizeof(a[0]))
#define sign_cast static_cast

class debug_values {
 private:
  int64_t m_debug_values[50];
 public:
  //debug_values(){for(uint32_t i=0;i<nelem(m_debug_values);++i) m_debug_values[i]=0;}
  debug_values();
  void Show(std::ostream& stream, int32_t indent=0);
  inline int64_t& operator[](uint32_t i){return m_debug_values[i];}
};

extern debug_values g_debug_values;

class CProcessMemory {
 private:
  static int64_t get_size(const char *p);
 public:
  int64_t Size;
  int64_t Peak;
  int64_t Data;
  CProcessMemory();
  void Write(std::ostream &stream, int32_t indent=0) const;
};
inline std::ostream& operator<<(std::ostream &stream, const CProcessMemory &pm){pm.Write(stream); return stream;}

class str2int_t : public std::map<std::string,int64_t> {
 public:
  void Write(std::ostream &s, int32_t indent=0) const;
};
inline std::ostream& operator<<(std::ostream &stream, const str2int_t &si){si.Write(stream); return stream;}

class str2str_t : public std::map<std::string,std::string> {
 public:
  void Write(std::ostream &s, int32_t indent=0) const;
};
inline std::ostream& operator<<(std::ostream &stream, const str2str_t &ss){ss.Write(stream); return stream;}

////////////////////////////////////////////////////////////////////////////////
class CCountLines {
private:
  int64_t m_n_all;
  int64_t m_n_chr;
  std::string m_tag;
  bool m_is_valid;
public:
  //inline CCountLines(const char *t){m_tag=t; m_n_all=m_n_chr=0; m_is_valid=Option().Find("verbose");}
  CCountLines(const char *t);
  inline void IncrementAll(){++m_n_all;}
  inline void IncrementChr(){++m_n_chr;}
  inline ~CCountLines(){if(m_is_valid){std::cerr<<"# "<<m_tag<<": #all="<<m_n_all<<", #chr="<<m_n_chr<<std::endl;}}
};

extern bool check_prefix(const char *prefix, const char *str);

BITVECTOR_LIB_END

/*
////////////////////////////////////////////////////////////////////////////////
class CKVStore {
 private:
  typedef std::map<std::string,std::string> mapping_t;
  mapping_t m_map;
  inline void Initialize(){}
 public:
  CKVStore(){Initialize();}
  bool Parse(const char **p);
  bool Read(const char *target_file);
  void Write(std::ostream &stream, int32_t indent=0) const;
  void Show (std::ostream &stream, int32_t indent=0) const {Write(stream,indent); stream<<std::endl;}
  const char *Find(const char *key) const;
  bool Add(const char *key, const char *value);
};
*/

#endif // _UTILITY_H_
