/**
 * @file    SequenceSet.cc
 * @brief   Parse alignments in FASTQ format, and store it.
 *
 * @author  Tomohiro Yasuda
 * @date    2011-4-30
 *
 */

#include <algorithm>
#include <limits>
#include <cstring>

#include "Option.h"
#include "Utility.h"
#include "FileReader.h"
//#include "SAMAlignment.h"
#include "SequenceSet.h"

using namespace BitVectorLib;

////////////////////////////////////////////////////////////////////////////////

CPositionStatistics::CPositionStatistics(){
  m_length_count=0;
  m_quality_sum=0;
  for(uint32_t i=0;i<nelems(m_base_count);++i) m_base_count[i]=0;
}

void CPositionStatistics::Write(std::ostream &stream, int32_t indent, uint64_t n_sequences) const{
  if(n_sequences==0) n_sequences=1;
  double average_quality = m_quality_sum;
  average_quality /= n_sequences;
  stream<<"{\"base_count\":[";
  for(uint32_t i=0;i<nelems(m_base_count);++i){
    if(i!=0) stream<<',';
    //stream<<m_base_count[i];
    double count = m_base_count[i];
    count /= n_sequences;
    stream<<count;
  }
  stream<<"], \"average_quality\":"<<average_quality<<"}";
}

////////////////////////////////////////////////////////////////////////////////

CSequenceSet::CSequenceSet(){
  m_n_sequences=0;
  m_total_length=0;
  m_min_quality=255;
  m_max_quality=0;
  // Sequence area
  m_size_in_last_page=AREA_SIZE;
  // Property area
  m_length_limit = Option().RequireInteger("max-read-length");
  m_position_statistics = new CPositionStatistics[m_length_limit];
}

CSequenceSet::~CSequenceSet(){
  foreach(std::vector<char*>,iter,m_sequence_areas) delete [] *iter;
  delete [] m_position_statistics;
}

int32_t CSequenceSet::Char2BaseCode(int32_t c){
  c &= 0xdf; // Capitalize
  //c>>=1;
  //c&=7;
  if(c=='A') return 0;
  if(c=='T') return 1;
  if(c=='G') return 2;
  if(c=='C') return 3;
  //if(nelems(m_base_count)<=4) Quit("Unexpected base: "<<static_cast<char>(c));
  if(c=='N') return 4;
  Warning("Unexpected symbol in sequence: "<<static_cast<char>(c));
  return 4;
}

void CSequenceSet::AddArea(){
  m_sequence_areas.push_back(new char[AREA_SIZE]);
  m_size_in_last_page=0;
}

char *CSequenceSet::MemoryToStoreSequence(int32_t len){
  len = (len+BASES_PER_BYTE-1)/BASES_PER_BYTE;
  if(m_size_in_last_page+len+1>=AREA_SIZE) AddArea();
  char *p = m_sequence_areas.back() + m_size_in_last_page;
  m_size_in_last_page += len;
  return p;
}

bool CSequenceSet::Read(CFileReader &file_reader){
  for(;;++m_n_sequences){
    // First header line
    const char *data = file_reader.GetLine();
    if(!data) return false;
    if(data[0]!='@') Quit("'@' expected at the head of the first line: "<<data);

    // Sequence line
    data = file_reader.GetLine();
    if(!data) Quit("Sequence expected at the second line.");
    uint32_t len=0;
    for(;data[len];++len){
      if(len>=m_length_limit) Quit("Longer sequence than predetermind: "<<m_n_sequences<<"-th sequence: "<<data);
      int32_t basecode = Char2BaseCode(data[len]);
      m_position_statistics[len].AddBaseCount(basecode);
    }
    m_position_statistics[len].IncrementLengthCount();
    m_total_length+=len;

    // Second header line
    data = file_reader.GetLine();
    if(!data) return false;
    if(data[0]!='+') Quit("'+' expected at the head of the third line.");

    // Quality line
    data = file_reader.GetLine();
    if(!data) Quit("Qualities expected at the fourth line.");
    for(uint32_t i=0;i<len;++i){
      if(! data[i]){
        Warning("Qualities of "<<m_n_sequences<<"-th sequence is shorter than sequence: "<<data);
        break;
      }
      uint32_t quality = data[i]-'!';
      if(m_min_quality>quality) m_min_quality=quality;
      if(m_max_quality<quality) m_max_quality=quality;
      m_position_statistics[i].AddQuality(quality);
    }
    if(data[len]) Warning("Qualities of "<<m_n_sequences<<"-th sequence is longer than sequence: "<<data);
  }
  return true; // never reached
}


void CSequenceSet::WriteStat(std::ostream &stream, int32_t indent) const {
  stream<<"{\"n_sequences\":"<<m_n_sequences
        <<", \"total_length\":"<<m_total_length
        <<", \"max_length\":"<<m_length_limit
        <<", \"max_quality\":"<<m_max_quality
        <<", \"min_quality\":"<<m_min_quality
        <<",";
  ShowIndent(stream,indent);
  stream<<"\"length_distribution\": [";
  bool is_first=true;
  uint32_t observed_max_length = 0;
  for(uint32_t i=0;i<m_length_limit;++i){
    if(m_position_statistics[i].LengthCount()==0) continue;
    if(is_first){is_first=false;}else{stream<<',';}
    stream<<'['<<i<<','<<m_position_statistics[i].LengthCount()<<']';
    observed_max_length=i;
  }
  stream<<"],";
  ShowIndent(stream,indent);
  stream<<"\"position_properties\": [";
  //for(uint32_t i=0;i<m_length_limit;++i){
  for(uint32_t i=0;i<=observed_max_length;++i){
    if(i!=0) stream<<',';
    ShowIndent(stream,indent+1);
    stream<<"\""<<i<<"\": ";
    m_position_statistics[i].Write(stream,0,m_n_sequences);
  }
  stream<<"]}";
}
