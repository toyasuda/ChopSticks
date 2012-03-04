/**
 * @file    DescriptiveStatistics.cc
 * @brief   Calculate descriptive statistics
 *
 * @author  Tomohiro Yasuda
 * @date    2011-1-9
 */

#ifndef _DESCRIPTIVE_STATISTICS_H_
#define _DESCRIPTIVE_STATISTICS_H_

#include <iostream>
#include <vector>
#include <stdint.h>
#include "Utility.h"

class CDescriptiveStatistics {
private:
  double m_sum;
  double m_sum2;
  bool m_is_empty;
  double m_cutoff_min;
  double m_cutoff_max;
public:
  void Update(double value);
  int64_t Cardinality;
  double Average;
  double Variance;
  double StandardDeviation;
  double Max;
  double Min;
  double Upper3Sigma() const {return Average+3*StandardDeviation;}
  double Lower3Sigma() const {return Average-3*StandardDeviation;}
  void Clear();
  void Finalize();
  CDescriptiveStatistics(){Clear();}
  void Initialize(int32_t cutoff_min, int32_t cutoff_max);
  inline void Initialize(){Initialize(-1,-2);}
  inline void Calculate(const std::vector<int64_t>& data, int32_t cmin=-1, int32_t cmax=-2){
    Initialize(cmin,cmax);
    foreach_const(std::vector<int64_t>, iter, data) Update(*iter);
    Finalize();
  }
  void Show(std::ostream &stream) const;
  static void test(std::istream& is, std::ostream& os);
};

#endif // _DESCRIPTIVE_STATISTICS_H_
