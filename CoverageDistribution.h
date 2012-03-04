/**
 * @file    CoverageArray.h
 * @brief   Build coverage array to improve SV prediction
 *
 * @author  Tomohiro Yasuda
 * @date    2011-6-18
 *
 */

#ifndef _COVERAGE_DISTRIBUTION_H_
#define _COVERAGE_DISTRIBUTION_H_

#include <stdint.h>
#include <vector>
#include "Utility.h"
#include "FileReader.h"
#include "GeneralFeature.h"

//class CGeneralFeature;

/*
#ifdef BITVECTOR_LIB_BEGIN
#define _COVERAGE_ARRAY_BV_NS_ BitVectorLib::
#else
#define _COVERAGE_ARRAY_BV_NS_
#endif // BITVECTOR_LIB_BEGIN
*/

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

#endif // _COVERAGE_DISTRIBUTION_H_
