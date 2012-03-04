/**
 * @file    CoverageBase.h
 * @brief   Base class for CCoverageArray
 *
 * @author  Tomohiro Yasuda
 * @date    2011-12-31
 *
 */

#ifndef _COVERAGE_BASE_H_
#define _COVERAGE_BASE_H_

#include <stdint.h>
#include <vector>
#include "Utility.h"
#include "FileReader.h"
#include "GeneralFeature.h"
#include "CoverageDistribution.h"

//class CGeneralFeature;
class CSAMAlignment;

#ifdef BITVECTOR_LIB_BEGIN
#define _COVERAGE_ARRAY_BV_NS_ BitVectorLib::
#else
#define _COVERAGE_ARRAY_BV_NS_
#endif // BITVECTOR_LIB_BEGIN

class CSAMReader {
 private:
  int32_t m_chr;
  int32_t m_margin_size;
  int64_t m_genome_size;
  int64_t m_max_position;
  int64_t m_min_position;
  int64_t m_n_total_bases;

  CKVStore m_library_threshold;
 protected:
  CSAMReader();
  inline bool LibraryAvailable() const {return !m_library_threshold.Empty();}
  inline const char *LibraryThreshold(const char *lib) const {return m_library_threshold.Find(lib);}
  inline void UpdateMaxPosition(int64_t p){if(m_max_position<p) m_max_position=p;}
  inline void UpdateMinPosition(int64_t p){
    if(m_min_position<0 || m_min_position>p){
      //std::cerr<<"UpdateMinPosition("<<p<<") -> "<<m_min_position<<std::endl;
      m_min_position=p;
    }
  }
  inline int64_t MaxPosition() const {return m_max_position;}
  inline int64_t MinPosition() const {return m_min_position;}
  inline int32_t MyChr() const {return m_chr;}
  inline int32_t MarginSize() const {return m_margin_size;}
  inline int64_t GenomeSize() const {return m_genome_size;}
  inline void AddTotalBases(int32_t nb){m_n_total_bases+=nb;}
  virtual void Treat(const CSAMAlignment& aln, const char *text)=0;
  virtual void TreatHeader(const char *text);
  void Initialize();
 public:
  void SetUp(int32_t chrNo);
  void ReadSAM(const char *sam_file, CChromosomeNormalizer& cn);
};

#endif // _COVERAGE_BASE_H_
