/**
 * @file    CoverageArray.h
 * @brief   Build coverage array to improve SV prediction
 *
 * @author  Tomohiro Yasuda
 * @date    2011-6-18
 *
 */

#ifndef _COVERAGE_ARRAY_H_
#define _COVERAGE_ARRAY_H_

#include <stdint.h>
#include <vector>
#include "Utility.h"
#include "FileReader.h"
#include "GeneralFeature.h"
#include "CoverageDistribution.h"
#include "SAMReader.h"

//class CGeneralFeature;
class CSAMAlignment;

#ifdef BITVECTOR_LIB_BEGIN
#define _COVERAGE_ARRAY_BV_NS_ BitVectorLib::
#else
#define _COVERAGE_ARRAY_BV_NS_
#endif // BITVECTOR_LIB_BEGIN

/*
class CCoverageDistribution {
private:
  std::vector<int32_t> m_distribution;
  double m_sum;
  double m_sum_squared;
  int32_t m_n_values;
  int32_t m_n_larger;
  static int32_t s_scale_factor;
  static int32_t s_index_limit;
  static int32_t s_window_size;
public:
  inline CCoverageDistribution(){Initialize();}
  void Initialize();
  void Increment(int32_t c);
  typedef std::vector<int32_t>::const_iterator const_iterator;
  typedef std::vector<int32_t>::iterator iterator;
  inline const_iterator begin() const {return m_distribution.begin();}
  inline const_iterator end  () const {return m_distribution.end();}
  inline double Average()  const {return m_sum/m_n_values;}
  inline double Variance() const {double a=Average(); return m_sum_squared/m_n_values-a*a;}
  void Write(std::ostream& stream) const;
};
*/

class CReadPositions {
public:
  class read_position_t {
  private:
    int32_t m_begin;
    int32_t m_end;
  public:
    inline read_position_t(int32_t b, int32_t e){m_begin=b; m_end=e;}
    bool operator<(const read_position_t& another);
    inline int32_t Begin() const {return m_begin;}
    inline int32_t End  () const {return m_end;}
    inline void Write(std::ostream &s) const {s<<'['<<m_begin<<','<<m_end<<']';}
  };
private:
  std::vector<read_position_t> m_read_positions;
  //protected:
  int32_t m_bin_size;
  //protected:
public:
  CReadPositions();
  bool TTest(int64_t begin_pos, int64_t end_pos, double mu, int32_t binsize) const;
  int64_t NAlignmentsIn(int64_t begin_pos, int64_t end_pos, int64_t *index) const;
  inline void AddRead(int64_t b,int64_t e){m_read_positions.push_back(read_position_t(b,e));}
  inline int64_t NReads() const {return m_read_positions.size();}
};
inline std::ostream& operator<<(std::ostream &s, const CReadPositions::read_position_t& rp){rp.Write(s); return s;}

class CCoverageArray : private CSAMReader {
  //class CCoverageArray {
public:
  typedef enum {REFINE_COVERAGE, REFINE_EDGE, REFINE_FRAGMENTED_EDGE} refine_type_t;
private:
  CReadPositions m_read_positions;

  int32_t m_max_coverage;
  double m_average_density;
  double m_fragment_threshold_rate;
  double m_margin_parameter;
  double m_refinement_coverage_threshold;
  uint16_t* m_coverage;
  static int32_t s_window_size;
  int32_t m_coverage_threshold;
  //int32_t m_ margin_size;
  int32_t m_max_chop_length;
  const char *m_output_format;

  void RefineByCoverage(std::vector<int64_t> *regions, int64_t front_pos, int64_t back_pos) const;
  int64_t RefineOneSide(int64_t front_pos, int64_t back_pos, refine_type_t rt) const;
  void RefineEdge(std::vector<int64_t> *regions, int64_t front_pos, int64_t back_pos, refine_type_t rt) const;
  void AnalyzeRegion(std::ostream &stream, int32_t no, int64_t front_pos, int64_t back_pos, refine_type_t rt) const;
  //int64_t MinPosition() const;
  double CoverageIn(int64_t front_pos, int64_t back_pos) const;
  mutable CCoverageDistribution m_covdist;
  bool CoverageDistribution(int64_t front_pos, int64_t back_pos, int64_t w, CCoverageDistribution *covdist, double my_coverage) const;
  virtual void Treat(const CSAMAlignment& aln, const char *text);
public:
  CCoverageArray();
  void SetUp(int32_t chrNo);
  void ReadSAM(const char *sam_file, CChromosomeNormalizer& cn);
  void Show(std::ostream &stream, int32_t indent=0) const;
  void RefineRegion(std::ostream &stream, refine_type_t rt, CGeneralFeatureVector& variants);
  void ShowCoverageDistribution(std::ostream &stream);
  void StatisticalAnalysis(int32_t coverage_unit) const;
};
//inline std::ostream& operator<<(std::ostream &s, const CCoverageArray::read_position_t& rp){rp.Write(s); return s;}

#endif // _COVERAGE_ARRAY_H_
