/**
 * @file    MappingReader.h
 * @brief   Read mapping into memory
 *
 * @author  Tomohiro Yasuda
 * @date    2010-10-2
 *
 */

#ifndef _MAPPING_READER_H_
#define _MAPPING_READER_H_

#include <iostream>
//#include <map>
//#include <string>
#include <algorithm>
#include <vector>
#include <stdint.h>

#include "FileReader.h"

/*
class CFileReader {
private:
  std::string m_filename;
  FILE *m_fp;
  FILE *m_my_fp;
  FILE *m_pipe_fp;
  char *m_buffer;
  static const int BUFFER_SIZE=65536;
  uint64_t m_line_number;
public:
  CFileReader(){m_fp=m_my_fp=m_pipe_fp=0; m_buffer=0;}
  ~CFileReader(){Close();}
  void Open(const char *filename);
  const char *GetLine();
  const char *GetContentLine();
  void Close();
  void Chomp();
  uint64_t Line() const {return m_line_number;}
};
*/

struct mapping_t {
  uint32_t Pos;
  uint16_t Len;
  char Strand;
  mapping_t(){Pos=Len=Strand=0;}
  bool IsValid() const {return Len>0;}
  void Purge(){Len=0;}
  bool operator<(const mapping_t &a) const {return Pos<a.Pos || (Pos==a.Pos && Len<a.Len);}
};

class mapping_vec_t : public std::vector<mapping_t> {
public:
  void Show(std::ostream &stream) const;
  void Sort(){sort(begin(),end());}
};

class CMappingReader {
private:
  mapping_vec_t *m_bins;
  int64_t m_expected_genome_length;
  int64_t m_n_mappings;
  static const int BIN_WIDTH=100; // must be no greater than 256
  int64_t NBins() const {return (m_expected_genome_length*2+BIN_WIDTH-1)/BIN_WIDTH;}
  int32_t m_begin_col;
  int32_t m_end_col;
  int32_t m_len_col;
  //void ParseFile(CFileReader* fr);
  void ParseFile();
  CFileReader m_file_reader;
  //bool ParseLine(const char *buffer, int64_t *pos, int* len, int* sign, int64_t linenumber);
  bool ParseLine(const char *buffer, mapping_t *mp, int64_t linenumber);
  //bool ParseAnother(CFileReader* fr);
public:
  bool ParseAnother(mapping_t* mp);
  CMappingReader(const char* gsize);
  ~CMappingReader();
  //void ParseFile(int begin_col, int end_col, int len_col, FILE* fp);
  //void ParseFile(int begin_col, int end_col, int len_col, CFileReader* fr);
  void SetUp(int begin_col, int end_col, int len_col, const char *file);
  void ParseFile(int begin_col, int end_col, int len_col, const char *file);
  void AddMapping(int64_t pos, int len, int sign);
  void SortBins();
  void Show(std::ostream &stream) const;
  //uint64_t Export(mapping_vec_t  *mv);
  void Export(mapping_vec_t  *mv);
};

//class CMappingReader;
class CMappingContainer {
public:
private:
  mapping_vec_t m_mappings;
  uint64_t m_genome_size;
  struct DataProperty_t {
    // Basic statistics
    uint64_t GenomeSize;
    uint64_t TotalReadLength;
    double Coverage;
    double AverageReadLength;
    int LongestReadLength;
    uint64_t NReads;
    void Show(std::ostream &stream) const;
    void SetUp(const mapping_vec_t &mappings);
    //inline double Significance(int gaplen) const {return exp(-Coverage*(gaplen-1)/AverageReadLength);}
    double Significance(int gaplen) const;
    //inline double P_of_falling_into(double x) const {return 1.0-exp(-x*NReads/GenomeSize);}
    inline double P_of_falling_into(double x) const {return x/GenomeSize;}
    void ShowSignificantLengths(std::ostream &stream) const;
    void SetUp(int G,int L,int N);
  };
  DataProperty_t m_data_property;
public:
  void Import(CMappingReader &mr);
  void Show(std::ostream &stream) const {m_mappings.Show(stream);}
  void FindGaps(mapping_vec_t *gp) const;
  void FindAberrentGaps(mapping_vec_t *gp) const;
  void FindBreak() const;
};

/*
class CLowCoverageFinder {
private:
  int *m_n_disappear;
  int m_buffer_size;
  int m_coverage;
  int64_t m_position;
  int64_t m_minpos;
  int m_direction;
  int m_threshold;
  //std::vector<int64_t> m_minimals;
public:
  struct region_t {
    int64_t Begin;
    int64_t End;
    int Coverage;
    region_t(){Begin=End=Coverage=0;}
    void SetUp(int64_t b, int64_t e, int c){Begin=b; End=e; Coverage=c;}
  };
  void RegisterMinimal();
public:
  //mapping_vec_t Minimals;
  typedef std::vector<region_t> region_vec_t;
  region_vec_t Minimals;
  CLowCoverageFinder(int sz, int threshold);
  ~CLowCoverageFinder(){delete [] m_n_disappear;}
  void CheckDisappearingReads(int64_t pos); //< call this before calling Add()
  void Add(int64_t pos, int64_t len);
  int Coverage() const {return m_coverage;}
};
*/

#endif // _MAPPING_READER_H_
