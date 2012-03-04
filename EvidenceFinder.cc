#include <iostream>
#include <cassert>
#include <cstdio>
#include <cstdlib>

#include "Option.h"
#include "Utility.h"
//#include "CoverageArray.h"
#include "SAMAlignment.h"
#include "EvidenceFinder.h"
#include <algorithm>
#include <limits>
#include <cmath>
#include <cstring>

#ifdef BITVECTOR_LIB_BEGIN
using namespace BitVectorLib;
#endif


////////////////////////////////////////////////////////////////////////////////
CEvidenceFinderFeatures::CEvidenceFinderFeatures(){
  for(int32_t i=0;i<=N_FEATURES;++i) m_count[i]=0;
}

void CEvidenceFinderFeatures::Show(std::ostream& stream) const {
  static const char *feature_names[] = {"?","paired","toofar","dangling_near_sv","dangling_independent","max_queue_size",0};
  stream<<"# EvidenceFinder: {";
  bool is_first=true;
  for(int32_t i=1;i<=N_FEATURES;++i){
    if(! feature_names[i]) break;
    if(is_first){is_first=false;}else{stream<<',';}
    stream<<'"'<<feature_names[i]<<"\":"<<m_count[i];
  }
  stream<<'}';
  stream<<std::endl;
}

////////////////////////////////////////////////////////////////////////////////
/*
bool CEvidenceFinder::comp_t::operator()(int32_t a, int32_t b) const {
  const CSAMAlignment& aref = m_finder->m_slots[a];
  const CSAMAlignment& bref = m_finder->m_slots[b];
  int32_t cmp = aref<bref;
  std::cout<<"aref="<<aref<<std::endl;
  std::cout<<"bref="<<bref<<std::endl;
  std::cout<<"cmp="<<cmp<<std::endl;
  return cmp;
}
*/

//CEvidenceFinder::CEvidenceFinder() : m_comparator(this) {
CEvidenceFinder::CEvidenceFinder(){
  m_n_bins = GenomeSize()/BIN_SIZE+1;
  std::cerr<<"# genome_size="<<GenomeSize()<<std::endl;
  if(m_n_bins<=1) m_n_bins=MAX_GENOME_SIZE/BIN_SIZE+1;
  std::cerr<<"# #bins="<<m_n_bins<<std::endl;
  //m_SVs = new std::vector<CGeneralFeature>::const_iterator[m_n_bins];
  m_SV_index = new int32_t[m_n_bins];
  m_slots = new aln_t[MAX_N_READS];
  m_slot_manager.SetUp(MAX_N_READS);
  m_dangling_distance = Option().RequireInteger("dangling-distance");
  std::cerr<<"# dangling reads will be reported if within "<<m_dangling_distance<<std::endl;
  const char *v = Option().Require("output-stderr");
  switch(*v){
  case 'D':
  case 'V':
    m_output_stderr=*v;
    break;
  default:
    Quit("Unknown output type: "<<v);
  }
}

void CEvidenceFinder::MakeBin(const CGeneralFeatureVector& variants){
  int32_t k=0;
  m_variants = &variants;
  for(uint32_t i=0;i<variants.size();++i){
    const CGeneralFeature& gf = variants[i];
    while(k<=gf.Start()/BIN_SIZE){
      if(k>=m_n_bins) Quit("Alignment beyond # of bins: m_n_bins="<<m_n_bins<<", k="<<k<<", gf="<<gf);
      m_SV_index[k]=i;
      ++k;
    }
  }
  while(k<m_n_bins){
    m_SV_index[k]=variants.size();
    ++k;
  }
  m_overlapping_reads = new std::vector<std::string>[variants.size()];
}

void CEvidenceFinder::TreatHeader(const char *text){
  std::cout<<text<<std::endl;
}

void CEvidenceFinder::Treat(const CSAMAlignment& downstream, const char *text){
  std::string qname(downstream.QName());

  _COVERAGE_ARRAY_BV_NS_ str2int_t::iterator iter=m_known.find(qname);
  if( m_known.find(qname)==m_known.end() ){
    int32_t slot=m_slot_manager.Allocate();
    if(slot<1) Quit("Too many reads to remember: "<<m_slot_manager.Size());
    m_slots[slot].Import(downstream,text);
    m_known[qname] = slot;
    QueuePush(slot);
    GetMax(MAX_QUEUE_SIZE,QueueSize());
  }else{
    aln_t& upstream=m_slots[(*iter).second];
    int32_t bin=upstream.Start()/BIN_SIZE;
    int64_t span_start = upstream.  End()  -MarginSize();
    int64_t span_end   = downstream.Start()+MarginSize();
    for(uint32_t i=m_SV_index[bin];i<m_variants->size();++i){
      const CGeneralFeature& sv = (*m_variants)[i];
      if(downstream.End() < sv.Start()) break;
      if(sv.Start()-m_dangling_distance<=span_start && // Mon Jan  9 23:17:34 2012
         span_start <= sv.Start() &&
         sv.End()<=span_end &&
         span_end<=sv.End()+m_dangling_distance){
        int32_t is_discordant=1;
        std::ostringstream oss;
        const char *lb=downstream.Option("LB");
        const char *thr=lb?LibraryThreshold(lb):0;
        if(!thr){
          lb=downstream.Option("RG");
          thr=lb?LibraryThreshold(lb):0;
        }
        //std::cerr<<"lb="<<(lb?lb:"N/A")<<", thr="<<(thr?thr:"N/A")<<std::endl;
        if(thr){
          int32_t w=std::atof(thr);
          int32_t len=downstream.End()-upstream.Start();
          if(len<w) is_discordant=0;
          oss<<"\tYT:f:"<<thr<<"\tYL:i:"<<len<<"\tYD:i:"<<is_discordant;
        }
        std::cout<<upstream.Text()<<oss.str()<<std::endl;
        std::cout<<text<<oss.str()<<std::endl;
        qname+=':';
        qname+=lb?lb:"NA";
        if(thr) qname+=(is_discordant?":d":":c");
        m_overlapping_reads[i].push_back(qname);
        break;
      }
    }
    upstream.SetUsed();
    m_known.erase(iter);
  }

  Flush(downstream.Start()-BIN_SIZE);
}

void CEvidenceFinder::Flush(int64_t position){
  while(! QueueEmpty() ){
    int32_t slot = QueueTop();
    aln_t& leftmost = m_slots[ slot ];
    if(position < leftmost.End()) break;
    int64_t d=leftmost.PNext()-leftmost.Pos();
    if(d<0) d=-d;
    if( leftmost.IsUsed() ) OneMore(PAIRED);
    else if(d>BIN_SIZE) OneMore(TOOFAR);
    else{
      // Look for the nearest
      int32_t bin=leftmost.Start()/BIN_SIZE-1;
      if(bin<0) bin=0;
      int32_t nearest_index=m_SV_index[bin];
      const CGeneralFeature& sv0=(*m_variants)[nearest_index];
      int64_t nearest_distance = interval_distance( sv0.Start(), sv0.End(), leftmost.Start(), leftmost.End() );
      for(uint32_t i=nearest_index;i<m_variants->size();++i){
        const CGeneralFeature& sv=(*m_variants)[i];
        int64_t distance = interval_distance( sv.Start(), sv.End(), leftmost.Start(), leftmost.End() );
        if(nearest_distance>distance){
          nearest_distance = distance;
          nearest_index=i;
        }
        if(distance>=BIN_SIZE) break;
      }

      if(nearest_distance>=m_dangling_distance) OneMore(DANGLING_INDEPENDENT);
      else{
        //std::cout<<"Dangling: "<<leftmost<<std::endl;
        //std::cout<<"nearest:  "<<(*nearest_iter)<<std::endl;
        //std::cout<<"nearest:  "<<(*m_variants)[nearest_index]<<std::endl;
        //m_feature_dangling_near_sv++;
        if(m_output_stderr=='D'){
          std::cerr<<"# nearest of the following:  "<<(*m_variants)[nearest_index]<<std::endl;
          std::cerr<<leftmost.Text()<<std::endl;
        }
        OneMore(DANGLING_NEAR_SV);
      }
    }
    QueuePop();
    m_slot_manager.Release(slot);
  }
}

void CEvidenceFinder::ReadSAM(const char *sam_file, CChromosomeNormalizer& cn){
  CSAMReader::ReadSAM(sam_file,cn);
  Flush(MaxPosition());
  if(m_output_stderr=='V'){
    for(uint32_t i=0;i<m_variants->size();++i){
      std::cerr<<(*m_variants)[i];
      const std::vector<std::string>& overlapping_reads = m_overlapping_reads[i];
      foreach_const(std::vector<std::string>, iter, overlapping_reads) std::cerr<<'\t'<<(*iter);
      std::cerr<<std::endl;
    }
  }
  CEvidenceFinderFeatures::Show(std::cerr);
}

