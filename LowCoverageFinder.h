/**
 * @file    MappingReader.h
 * @brief   Read mapping into memory
 *
 * @author  Tomohiro Yasuda
 * @date    2010-10-2
 *
 */

#ifndef _LOW_COVERAGE_FINDER_H_
#define _LOW_COVERAGE_FINDER_H_

#include <iostream>
#include <deque>
#include <stdint.h>
#include "DescriptiveStatistics.h"

class CLowCoverageFinder {
protected:
  int64_t m_previous_changepoint;
  int32_t m_previous_coverage;
  int32_t *m_n_disappear;
  int32_t m_buffer_size;
  int32_t m_coverage;
  int64_t m_position;
  int64_t m_final_position;
  std::vector<int32_t> m_read_lengths_at_current_position;
  std::vector<int64_t> m_coverage_distribution;
  CDescriptiveStatistics m_coverage_statistics;
  //static const int32_t MAX_COVERAGE=10000;
  static const int32_t MAX_COVERAGE=500;
 public:
  class event_handler_t {
  public:
    virtual void EventChangeCoverage(const CLowCoverageFinder& lcf, int32_t diff){}
    virtual void EventEnterReads    (const CLowCoverageFinder& lcf, int32_t diff){}
    virtual void EventExitReads     (const CLowCoverageFinder& lcf, int32_t diff){}
    virtual void EventOneRead       (const CLowCoverageFinder& lcf, int32_t len){}
    virtual void EventFinal         (const CLowCoverageFinder& lcf){}
  };
private:
  std::vector<event_handler_t*> m_event_handlers;
public:
  void CommitUntil(int64_t pos);
  int32_t AddNewReads();
  void AddEventHandler(event_handler_t* eh){if(eh) m_event_handlers.push_back(eh);}
  CLowCoverageFinder(int sz, event_handler_t* eh=0);
  ~CLowCoverageFinder(){delete [] m_n_disappear;}
  void Add(int64_t pos, int64_t len);
  int Coverage() const {return m_coverage;}
  int64_t PreviousChangepoint() const {return m_previous_changepoint;}
  int32_t PreviousCoverage() const {return m_previous_coverage;}
  int64_t Position() const {return m_position;}
  void Finalize();
  void Analyze(const mapping_vec_t& mappings);
  void Analyze(CMappingReader *mreader);
  void ShowCoverageDistribution(std::ostream &stream, const char *prefix) const;
};



class CMinimalFinder : public CLowCoverageFinder::event_handler_t {
private:
  int64_t m_region_start;
  int m_previous_diff;
  int m_threshold;
  int m_lower_threshold;
  bool m_lower_threshold_observed;
  int64_t m_area;
public:
  struct region_t {
    int64_t Begin;
    int64_t End;
    int Coverage;
    region_t(){SetUp(0,0,0);}
    void SetUp(int64_t b, int64_t e, int c){Begin=b; End=e; Coverage=c;}
  };
  typedef std::vector<region_t> region_vec_t;
  region_vec_t Minimals;
  void AddMinimal(int64_t b, int64_t e, int32_t c);
  CMinimalFinder(int threshold, int lower_threshold=-1);
  void EventChangeCoverage(const CLowCoverageFinder &lcf, int32_t diff);
  void EventFinal         (const CLowCoverageFinder &lcf);
  void ShowMinimals(std::ostream &stream) const;
};
extern std::ostream& operator<<(std::ostream &stream, const CMinimalFinder::region_t& r);
//extern std::ostream& operator<<(std::ostream &stream, const region_t& r);

class CRegionRefiner : public CLowCoverageFinder::event_handler_t {
 public:
  class region_t {
  private:
    std::vector<int64_t> m_boundaries;
  public:
    int64_t Begin;
    int64_t End;
    double Coverage;
    region_t(){Begin=End=0; Coverage=0;}
    void SetUp(int64_t b, int64_t e){Begin=b; End=e; m_boundaries.push_back(b); m_boundaries.push_back(e);}
    void UpdateCoverage(const CLowCoverageFinder &lcf);
    void Refine(int64_t begin, int64_t end);
    static void Refine_test();
    bool IsValid() const {return !m_boundaries.empty();}
    void Show(std::ostream &stream) const;
  };
 private:
  std::deque<region_t> m_region_deque;
  std::vector<region_t> m_modified_regions;
  int32_t m_max_read_length;
  int64_t m_previous_begin;
  CFileReader m_prediction_reader;
  void FlushFinishedResults(int64_t pos);
  void EventOneRead(const CLowCoverageFinder& lcf, int32_t len);
  void EventChangeCoverage(const CLowCoverageFinder &lcf, int32_t diff);
  void EventFinal  (const CLowCoverageFinder& lcf);
  void ReadPrediction(int64_t examine_begin, int64_t examine_end);
 public:
  CRegionRefiner(const char *file);
  void ShowResults(std::ostream &stream, const char *prefix) const;
};
inline std::ostream& operator<<(std::ostream &s, const CRegionRefiner::region_t& r){r.Show(s); return s;}


#endif // _LOW_COVERAGE_FINDER_H_
