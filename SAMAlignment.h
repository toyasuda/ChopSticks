/**
 * @file    SAMAlignment.h
 * @brief   Parse alignments in SAM format, and represent it.
 *
 * @author  Tomohiro Yasuda
 * @date    2011-4-15
 *
 */

#ifndef _SAM_ALIGNMENT_H_
#define _SAM_ALIGNMENT_H_

#include <iostream>
#include <string>
#include <map>
#include <vector>
#include <stdint.h>
#include "Utility.h"

#ifdef BITVECTOR_LIB_BEGIN
#define _COVERAGE_ARRAY_BV_NS_ BitVectorLib::
#else
#define _COVERAGE_ARRAY_BV_NS_
#endif // BITVECTOR_LIB_BEGIN

class CCigarString {
private:
  struct operation_t {
    int32_t Len;
    char Sym;
    operation_t(){Len=0; Sym=0;}
    bool IsClip() const {return Sym=='S' && Len>=s_clip_length_threshold;}
    void Show(std::ostream &s, int32_t indent=0) const {s<<Sym<<Len;}
  };
  std::vector<operation_t> m_editscript;
  static int32_t s_clip_length_threshold;
  int32_t m_length;
  int32_t m_alignment_length;
  int32_t m_flag;
  int32_t __IsClipped();
public:
  CCigarString();
  CCigarString(const CCigarString& cs);
  bool Parse(const char *cigar_string);
  void Write(std::ostream &s, int32_t indent=0) const;
  int64_t ClipPositionLeft (int64_t p=0) const {return m_flag&1? p: -1;}
  int64_t ClipPositionRight(int64_t p=0) const;
  int32_t Length() const {return m_length;}
  int32_t IsClipped() const {return m_flag;}
  int32_t AlignmentLength() const {return m_alignment_length;}
};
inline std::ostream& operator<<(std::ostream& s, const CCigarString& x){x.Write(s,0); return s;}

class CGenomePosition {
private:
  static const uint64_t BOUNDARY=32;
  static std::map<std::string,uint32_t>* s_chromosome2int;
public:
  static uint64_t GenomePosition(const char *ch, int32_t pos) {return GenomePosition(Chromosome(ch),pos);}
  static uint64_t GenomePosition(uint64_t    ch, int32_t pos) {return (ch<<BOUNDARY)+pos;}
  static uint64_t Chromosome(uint64_t gpos) {return gpos>>BOUNDARY;}
  static uint64_t Position  (uint64_t gpos) {return gpos & ((1ULL<<BOUNDARY)-1);}
  static uint64_t Chromosome(const char *ch);
};

class CClipPositions {
private:
  std::vector<int64_t> m_positions;
  int64_t m_cluster_ID;
  static int32_t s_cluster_size_threshold;
  static int32_t s_cluster_gap_threshold;
  std::string m_prefix;
public:
  //CClipPositions(){m_cluster_ID=0;}
  CClipPositions(const char *p);
  void AddPosition(int64_t gpos){m_positions.push_back(gpos);}
  void FindClusters();
  int64_t Size() const {return m_cluster_ID;}
  void MakeOneCluster(std::vector<int64_t>::const_iterator b, std::vector<int64_t>::const_iterator e);
};

class CSAMAlignment {
private:
  int32_t m_chr;
  std::string m_qname;
  int32_t m_flag;
  std::string m_rname;
  int64_t m_pos;
  int32_t m_mapq;
  std::string m_cigar;
  std::string m_rnext;
  int64_t m_pnext;
  int32_t m_tlen;
  std::string m_seq;
  std::string m_qual;
  CCigarString m_cigar_info;
  _COVERAGE_ARRAY_BV_NS_ str2str_t m_options;

  void Initialize(){m_chr=m_flag=m_pos=m_mapq=m_pnext=m_tlen=0;}

protected:
  void Copy(const CSAMAlignment& sa);

public:
  static const int32_t MULTIPLE_FRAGMENTS=0x01;
  static const int32_t PROPERLY_ALIGNED=0x02;
  static const int32_t UNMAPPED=0x04;
  static const int32_t NEXT_UNMAPPED=0x08;
  static const int32_t REVERSE_STRAND=0x10;
  static const int32_t NEXT_REVERSE_STRAND=0x20;
  static const int32_t FIRST_FRAGMENT=0x40;
  static const int32_t LAST_FRAGMENT=0x80;
  static const int32_t SECONDARY_ALIGNMENT=0x100;
  static const int32_t DISQUALIFIED=0x200;
  static const int32_t DUPLICATE=0x400;

  inline bool IsMultipleFragments()  const {return Flag() & MULTIPLE_FRAGMENTS;}
  inline bool IsProperlyAligned()    const {return Flag() & PROPERLY_ALIGNED;}
  inline bool IsUnmapped()           const {return Flag() & UNMAPPED;}
  inline bool IsNextUnmapped()       const {return Flag() & NEXT_UNMAPPED;}
  inline bool IsReverseStrand()      const {return Flag() & REVERSE_STRAND;}
  inline bool IsNextReversestrand()  const {return Flag() & NEXT_REVERSE_STRAND;}
  inline bool IsFirstFragment()      const {return Flag() & FIRST_FRAGMENT;}
  inline bool IsLastFragment()       const {return Flag() & LAST_FRAGMENT;}
  inline bool IsSecondaryAlignment() const {return Flag() & SECONDARY_ALIGNMENT;}
  inline bool IsDisqualified()       const {return Flag() & DISQUALIFIED;}
  inline bool IsDuplicate()          const {return Flag() & DUPLICATE;}

  const char* Option(const std::string& k) const;

  inline CSAMAlignment(){Initialize();}
  inline CSAMAlignment(const char **pp){Initialize(); Parse(pp);}
  inline CSAMAlignment(const char  *p ){Initialize(); Parse(p );}
  inline CSAMAlignment(const CSAMAlignment& sa){Copy(sa);}
  inline CSAMAlignment& operator=(const CSAMAlignment& sa){Copy(sa); return *this;}
  bool operator<(const CSAMAlignment& sa) const;
  bool Parse(const char** pp);
  bool Parse(const char* p){const char *q=p; return Parse(&q);}
  void Write(std::ostream &stream, int32_t indent=0) const;
  void Show (std::ostream &stream, int32_t indent=0) const {Write(stream,indent); stream<<std::endl;}
  // Field data
  const char* QName() const {return m_qname.c_str();}
  int32_t Flag() const {return m_flag;}
  const char* RName() const {return m_rname.c_str();}
  int64_t Pos()  const {return m_pos;}
  int32_t Mapq() const {return m_mapq;}
  const char* Cigar() const {return m_cigar.c_str();}
  const char* RNext() const {return m_rnext.c_str();}
  int64_t PNext() const {return m_pnext;}
  int64_t TLen()  const {return m_tlen;}
  const char* Seq()  const {return m_seq.c_str();}
  const char* Qual() const {return m_qual.c_str();}

  inline bool Unmapped() const {return m_flag & 0x04;}
  inline int64_t Start() const {return Pos();}
  inline int64_t End  () const {return Pos()? Pos()+m_cigar_info.AlignmentLength()-1: -1;}
};
inline std::ostream& operator<<(std::ostream &s, const CSAMAlignment& a){a.Write(s); return s;}

class CQualityAnalyzer {
 private:
  int32_t m_min_quality_symbol;
  static const int32_t GREATER_THAN_ANY_QUALITY=1024*1024*1024;
 public:
  inline CQualityAnalyzer(int32_t min=0)  {m_min_quality_symbol=min?min:'!';}
  inline CQualityAnalyzer(const char *min){m_min_quality_symbol=(min&&min[0])?min[0]:'!';}
  double Average(const char *q) const;
  int32_t Min(const char *q) const;
  int32_t Max(const char *q) const;
  inline double  Average(const CSAMAlignment& a) const {return Average(a.Qual());}
  inline int32_t Min(const CSAMAlignment& a)     const {return Min(a.Qual());}
  inline int32_t Max(const CSAMAlignment& a)     const {return Max(a.Qual());}
};

/*
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
};
//inline std::ostream& operator<<(std::ostream &s, const CGeneralFeature& a){a.Write(s); return s;}
*/

#endif // _SAM_ALIGNMENT_H_
