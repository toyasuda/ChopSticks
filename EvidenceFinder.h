/**
 * @file    EvidenceFinder.h
 * @brief   Find evidences for each deletions
 *
 * @author  Tomohiro Yasuda
 * @date    2012-01-01
 *
 */

#ifndef _EVIDENCE_FINDER_H_
#define _EVIDENCE_FINDER_H_

#include <stdint.h>
#include <vector>
#include <queue>
#include "Utility.h"
#include "FileReader.h"
#include "GeneralFeature.h"
#include "SAMReader.h"
#include "SAMAlignment.h"
#include "Tool.h"

#ifdef BITVECTOR_LIB_BEGIN
#define _COVERAGE_ARRAY_BV_NS_ BitVectorLib::
#else
#define _COVERAGE_ARRAY_BV_NS_
#endif // BITVECTOR_LIB_BEGIN

class CEvidenceFinderFeatures {
  //int32_t m_feature_paired;
  //int32_t m_feature_toofar;
  //int32_t m_feature_dangling_near_sv;
  //int32_t m_feature_independent;
 protected:
  static const int32_t PAIRED=1;
  static const int32_t TOOFAR=2;
  static const int32_t DANGLING_NEAR_SV=3;
  static const int32_t DANGLING_INDEPENDENT=4;
  static const int32_t MAX_QUEUE_SIZE=5;
  inline void OneMore(int32_t f){Assert(0<=f && f<=N_FEATURES); m_count[f]++;}
  inline void GetMax (int32_t f, int64_t x){Assert(0<=f && f<=N_FEATURES); if(m_count[f]<x) m_count[f]=x;}
 private:
  static const int32_t N_FEATURES=5;
  int32_t m_count[N_FEATURES+1];
  static const int32_t NIL=0;
 public:
  CEvidenceFinderFeatures();
  void Show(std::ostream& stream) const;
};


class CEvidenceFinder : public CSAMReader, public CEvidenceFinderFeatures {
private:
  /*
  class comp_t {
  private:
    const class CEvidenceFinder* m_finder;
  public:
    inline comp_t(const class CEvidenceFinder *f){m_finder=f;}
    bool operator()(int32_t a, int32_t b) const;
  };
  friend class comp_t;
  comp_t m_comparator;
  */
  class aln_t : public CSAMAlignment {
  private:
    std::string m_text;
    bool m_used;
  public:
    inline aln_t(){m_used=false;}
    //inline aln_t(const CSAMAlignment& a){Copy(a); Used=false;}
    //inline aln_t& operator=(const CSAMAlignment& a){Copy(a); Used=false; return *this;}
    inline void Import(const CSAMAlignment& a, const char *text){CSAMAlignment::Copy(a); m_used=false; m_text=text;}
    inline bool IsUsed() const {return m_used;}
    inline void SetUsed(){m_used=true;}
    inline const std::string& Text() const {return m_text;}
  };
  /*
  std::priority_queue<int32_t,std::vector<int32_t>, comp_t> m_position_queue;
  inline void QueuePush(int32_t slot){m_position_queue.push(slot);}
  inline void QueuePop(){m_position_queue.pop();}
  inline bool QueueEmpty() const {return m_position_queue.empty();}
  inline int32_t QueueTop(){return m_position_queue.top();}
  inline int32_t QueueSize(){return m_position_queue.size();}
  */
  static const int32_t MAX_GENOME_SIZE=300*1000*1000;
  //static const int32_t MAX_N_READS=1000*1000;
  static const int32_t MAX_N_READS=1000*1000;
  static const int32_t BIN_SIZE=1000000;
  //aln_t m_slots[MAX_N_READS];
  char m_output_stderr;
  aln_t* m_slots;
  std::vector<std::string>* m_overlapping_reads;
  int32_t m_n_bins;
  int32_t m_dangling_distance;
  std::deque<int32_t> m_position_queue;
  inline void QueuePush(int32_t slot){m_position_queue.push_back(slot);}
  inline void QueuePop(){m_position_queue.pop_front();}
  inline bool QueueEmpty() const {return m_position_queue.empty();}
  inline int32_t QueueTop(){return m_position_queue.front();}
  inline int32_t QueueSize(){return m_position_queue.size();}

  //std::vector<CGeneralFeature>::const_iterator *m_SVs;
  int32_t* m_SV_index;
  //std::vector<CGeneralFeature>::const_iterator m_end_iter;
  const std::vector<CGeneralFeature>* m_variants;
  _COVERAGE_ARRAY_BV_NS_ str2int_t m_known;
  CSlotManager m_slot_manager;
  void Flush(int64_t position);
  void Treat(const CSAMAlignment& aln, const char *text);
  void TreatHeader(const char *text);
public:
  CEvidenceFinder();
  ~CEvidenceFinder(){delete [] m_slots; delete [] m_overlapping_reads;}
  void MakeBin(const CGeneralFeatureVector& variants);
  void ReadSAM(const char *sam_file, CChromosomeNormalizer& cn);
};

#endif // _EVIDENCE_FINDER_H_
