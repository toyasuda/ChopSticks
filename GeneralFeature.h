/**
 * @file    General feature.h
 * @brief   Parse general features in GFF format, and represent it.
 *
 * @author  Tomohiro Yasuda
 * @date    2011-12-31
 *
 */

#ifndef _GENERAL_FEATURE_H_
#define _GENERAL_FEATURE_H_

#include <iostream>
#include <string>
#include <vector>
#include <stdint.h>

#ifdef BITVECTOR_LIB_BEGIN
#define _COVERAGE_ARRAY_BV_NS_ BitVectorLib::
#else
#define _COVERAGE_ARRAY_BV_NS_
#endif // BITVECTOR_LIB_BEGIN

class CGeneralFeature {
private:
  std::string m_name;
  std::string m_source;
  std::string m_type;
  int64_t m_start;
  int64_t m_end;
  int64_t m_score;
  int32_t m_strand;
  std::string m_frame;
  std::string m_group;
  //void Initialize(){m_start=m_end=m_score=m_strand=0;}
  void Clear();
  bool ParseBED(const char** pp);
  bool ParseGFF(const char** pp);
  bool ParseBED(const char* p){const char *q=p; return ParseBED(&q);}
  bool ParseGFF(const char* p){const char *q=p; return ParseGFF(&q);}
public:
  static const int32_t UNKNOWN=0;
  static const int32_t GFF=1;
  static const int32_t BED=2;
  CGeneralFeature(){Clear();}
  //CGeneralFeature(const char **pp){Clear(); ParseGFF(pp);}
  //CGeneralFeature(const char  *p ){Clear(); ParseGFF(p );}
  static int32_t Filetype(const char *filetype);
  bool Parse(const char* p, int32_t filetype)   __attribute__((warn_unused_result));
  void WriteGFF(std::ostream &stream, int32_t indent=0) const;
  void WriteBED(std::ostream &stream, int32_t indent=0) const;
  //void Show (std::ostream &stream, int32_t indent=0) const {Write(stream,indent); stream<<std::endl;}
  // Field data
  const char* Name()   const {return m_name.c_str();}
  const char* Source() const {return m_source.c_str();}
  const char* Type()   const {return m_type.c_str();}
  int64_t Start()  const {return m_start;}
  int64_t End()    const {return m_end;}
  int64_t Score()  const {return m_score;}
  int32_t Strand() const {return m_strand;}
  const char* Frame()  const {return m_frame.c_str();}
  const char* Group()  const {return m_group.c_str();}
  bool operator<(const CGeneralFeature &a) const;
  inline void Show(std::ostream &stream) const {WriteBED(stream,0); stream<<std::endl;}
};
inline std::ostream& operator<<(std::ostream& s, const CGeneralFeature& gf){gf.WriteBED(s); return s;}

class CChromosomeNormalizer {
 private:
  CKVStore m_kvs;
  _COVERAGE_ARRAY_BV_NS_ str2int_t m_unknown_refname;
 public:
  void Read(const char* table_filename){m_kvs.Read(table_filename);}
  int32_t Chr(const char *refname);
  //const char* Find(const char *refname){return m_kvs.Find(refname);}
  void ShowUnknown(std::ostream &stream) const {stream<<m_unknown_refname<<std::endl;}
};

class CGeneralFeatureVector: public std::vector<CGeneralFeature> {
 public:
  //void ReadGFF(int32_t chr, class CCoverageArray& ca, const char *gff_file);
  void ReadGFF(int32_t chr, class CChromosomeNormalizer& cn, const char *gff_file);
  void Sort();
};


#endif // _GENERAL_FEATURE_H_
