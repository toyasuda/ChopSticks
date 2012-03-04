/**
 * @file    SequenceSet.h
 * @brief   Parse alignments in FASTQ format, and store it.
 *
 * @author  Tomohiro Yasuda
 * @date    2011-4-30
 *
 */

#ifndef _SEQUENCE_SET_H_
#define _SEQUENCE_SET_H_

#include <iostream>
//#include <string>
//#include <map>
#include <vector>
#include <stdint.h>

class CFileReader;

//////////////////////////////////////////////////////////////////////
class CSequenceAnchor {
 private:
  uint64_t m_position;
  uint8_t m_length;
 public:
  bool Read(CFileReader &file_reader);
  void WriteStat(std::ostream &s, int32_t indent=0) const;
};

//////////////////////////////////////////////////////////////////////
class CPositionStatistics {
 private:
  uint64_t m_length_count;
  uint64_t m_quality_sum;
  uint64_t m_base_count[5];
 public:
  CPositionStatistics();
  void IncrementLengthCount(){++m_length_count;}
  void AddQuality  (int32_t q){m_quality_sum+=q;}
  void AddBaseCount(int32_t basecode){m_base_count[basecode]+=1;}
  uint64_t LengthCount() const {return m_length_count;}
  uint64_t QualitySum() const {return m_quality_sum;}
  uint64_t BaseCount(int32_t basecode) const {return m_base_count[basecode];}
  void Write(std::ostream &s, int32_t indent=0, uint64_t n_sequences=1) const;
};

//////////////////////////////////////////////////////////////////////
class CSequenceSet {
private:
  //char **m_sequence_area;
  //uint64_t m_n_areas;
  std::vector<char*> m_sequence_areas;
  CPositionStatistics *m_position_statistics;
  uint32_t m_size_in_last_page;
  uint64_t m_n_sequences;
  uint64_t m_total_length;
  uint32_t m_length_limit;
  uint32_t m_min_quality;
  uint32_t m_max_quality;
  static const uint32_t AREA_SIZE =16*1024*1024;
  std::vector<CSequenceAnchor> m_sequence_anchor;
  void AddArea();
  char *MemoryToStoreSequence(int32_t len);
  static const int32_t BASES_PER_BYTE=4;
  static int32_t Char2BaseCode(int32_t c);
public:
  CSequenceSet();
  ~CSequenceSet();
  bool Read(CFileReader &file_reader);
  void Write(std::ostream &s, int32_t indent=0) const;
  void WriteStat(std::ostream &s, int32_t indent=0) const;
  void ShowStat(std::ostream &s) const {WriteStat(s,1); s<<std::endl;}
};
inline std::ostream& operator<<(std::ostream& s, const CSequenceSet& x){x.Write(s,0); return s;}


#endif // _SEQUENCE_SET_H_
