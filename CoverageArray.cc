#include <iostream>
#include <cassert>
#include <cstdio>
#include <cstdlib>

#include "Option.h"
#include "Utility.h"
#include "CoverageArray.h"
#include "SAMAlignment.h"
#include <algorithm>
#include <limits>
#include <cmath>
#include <cstring>

#ifdef BITVECTOR_LIB_BEGIN
using namespace BitVectorLib;
#endif

int32_t CCoverageArray::s_window_size =-1;

#define dT(x) 

////////////////////////////////////////////////////////////////////////////////
CReadPositions::CReadPositions(){
  m_bin_size = Option().RequireInteger("bin-size");
}

int64_t CReadPositions::NAlignmentsIn(int64_t begin_pos, int64_t end_pos, int64_t *index) const {
  //int64_t r_limit=m_read_positions.size();
  int64_t r_limit=NReads();
  int64_t r1=*index;
  for(; 0<=r1 && r1<r_limit && m_read_positions[r1].Begin()>=begin_pos; --r1);
  for(;          r1<r_limit && m_read_positions[r1].Begin()< begin_pos; ++r1);
  int64_t r2=r1;
  for(;          r2<r_limit && m_read_positions[r2].Begin()<=end_pos;   ++r2);
  *index=r2;
  return r2-r1;
}

bool CReadPositions::TTest(int64_t begin_pos, int64_t end_pos, double mu, int32_t binsize) const{
  std::cerr<<"CCoverageArray::TTest("<<begin_pos<<"-"<<end_pos<<":"<<end_pos-begin_pos+1<<", "<<mu<<", "<<binsize<<")"<<std::endl;
  mu *= binsize;
  int64_t r1 = 0;
  int64_t r2 = NReads()-1;
  int32_t n=1;
  for(;r2-r1>1;++n){
    int64_t m=(r1+r2)/2;
    if(begin_pos<=m_read_positions[m].Begin()){ r2=m; }else{ r1=m; }
  }
  int64_t r0=r1;
  int64_t n_aln = NAlignmentsIn(begin_pos,end_pos,&r1);
  std::cerr<<"Total "<<n_aln<<" reads."<<std::endl;

  for(int64_t pos=begin_pos; pos<end_pos; pos+=binsize){
    int64_t pos2 = pos+binsize-1;
    if(pos2>=end_pos) pos2=end_pos-1;
    int64_t n = NAlignmentsIn(pos,pos2,&r0);
    std::cerr<<((pos-begin_pos)/binsize+1)<<": pos="<<pos<<": n="<<n<<std::endl;
  }
  return true;
}


////////////////////////////////////////////////////////////////////////////////
CCoverageArray::CCoverageArray(){
  m_coverage=0;
  m_max_coverage=0;
  m_average_density=0;
  m_refinement_coverage_threshold = Option().RequireInteger("refinement-threshold");
  m_margin_parameter = Option().RequireInteger("margin-parameter");
  m_fragment_threshold_rate = Option().RequireDouble("fragment-threshold");
  if(s_window_size<0) s_window_size=Option().RequireInteger("coverage-window");
  m_output_format=Option().Require("output-format");
  std::cerr<<"# output format="<<m_output_format<<std::endl;
  m_max_chop_length=Option().RequireInteger("max-chop-length");
  std::cerr<<"# max chop length="<<m_max_chop_length<<std::endl;
}

void CCoverageArray::ShowCoverageDistribution(std::ostream &stream){
  std::vector<int64_t> coverage_distribution(m_max_coverage+1,0);
  for(int64_t position=MinPosition();position<=MaxPosition();++position){
    ++coverage_distribution[ m_coverage[position] ];
  }
  uint32_t coverage_upper_bound = Option().RequireInteger("coverage-upper-bound");
  for(uint32_t i=0;i<coverage_distribution.size();++i){
    if(i>=coverage_upper_bound){
      int64_t c=0;
      for(uint32_t j=i;j<coverage_distribution.size();++j) c+=m_coverage[j];
      stream<<"D\t"<<i<<"\t"<<c<<std::endl;
      break;
    }
    stream<<"D\t"<<i<<"\t"<<coverage_distribution[i]<<std::endl;
  }
}

void CCoverageArray::SetUp(int32_t chrNo){
  CSAMReader::SetUp(chrNo);

  // Get array
  m_coverage = new uint16_t[GenomeSize()+1];
  for(int64_t i=0;i<=GenomeSize();++i) m_coverage[i]=0;

  // Coverage threshold
  m_coverage_threshold = Option().RequireInteger("coverage-threshold");
  m_max_coverage=0;
  //m_ margin_size = Option().RequireInteger("margin-size");
  m_margin_parameter = Option().RequireInteger("margin-parameter");
  m_fragment_threshold_rate = Option().RequireDouble("fragment-threshold");
}

void CCoverageArray::RefineByCoverage(std::vector<int64_t> *regions, int64_t front_pos, int64_t back_pos) const{
  assert(regions);
  assert(front_pos<=back_pos);
  for(int64_t i=front_pos;i<=back_pos;++i){
    if((regions->size() & 1)==0 && m_coverage[i]<m_coverage_threshold){
      regions->push_back(i);
    }else if((regions->size() & 1)==1 && m_coverage[i]>=m_coverage_threshold){
      // Close the latest region
      regions->push_back(i);
    }
  }
  int64_t end_pos=back_pos+1;
  if((regions->size()&1)!=0) regions->push_back(end_pos);
}

bool CCoverageArray::CoverageDistribution(int64_t front_pos, int64_t back_pos, int64_t w, CCoverageDistribution *covdist, double my_coverage) const{
  Assert(covdist);
  if(s_window_size<=0) return true; // If window is not used, do nothing.
  //return true;
  std::cout<<"C\t"<<front_pos<<'\t'<<back_pos<<'\t'<<w<<'\t'<<my_coverage<<'\t';

  covdist->Initialize();
  if(front_pos>back_pos){int64_t x=front_pos; front_pos=back_pos; back_pos=x;}
  double c=0;
  int64_t len=back_pos-front_pos+1;
  // First w-1 bases
  int32_t i=0;
  for(;i<w;++i) c+=m_coverage[front_pos+i];
  // Rest len-w+1 bases
  for(;i<len;++i){
    covdist->Increment(c/w);
    c+=m_coverage[front_pos+i];
    c-=m_coverage[front_pos+i-w];
  }

  covdist->Write(std::cout);
  std::cout<<std::endl;
  return true;
}

int64_t CCoverageArray::RefineOneSide(int64_t front_pos, int64_t back_pos, refine_type_t rt) const{
  int32_t direction = 1;
  int64_t lower_limit=front_pos;
  int64_t upper_limit=back_pos;
  if(front_pos>back_pos){
    direction=-1;
    lower_limit=back_pos;
    upper_limit=front_pos;
  }
  double coverage_threshold=m_coverage_threshold;
  dT("1:coverage_threshold="<<coverage_threshold<<"("<<m_coverage_threshold<<")");
  if(MarginSize()>0){
    std::cerr
      <<"target: "
      <<front_pos<<'\t'<<back_pos<<'\t'<<front_pos-direction<<'\t'<<front_pos-direction*MarginSize()<<std::endl;
    coverage_threshold = CoverageIn(front_pos-direction,front_pos-direction*MarginSize());
    coverage_threshold *= m_margin_parameter; // CHECK, FIXIT
  }

  dT("2:coverage_threshold="<<coverage_threshold<<"("<<m_coverage_threshold<<")");
  // First, truncate contiguous high coverage region.
  for(;front_pos!=back_pos && m_coverage[front_pos]>=coverage_threshold;front_pos+=direction);
  int64_t x=front_pos;
  // Then, truncate gap-high coverage regions.
  dT("first_x="<<x);
  while((back_pos-x)*direction>=0){
    double c=0;
    double fragment_begin=x;
    dT("front="<<front_pos<<", back_pos="<<back_pos);
    if(x==back_pos) break;
    dT("Skipping low cov region");
    for(;lower_limit<=x && x<=upper_limit && m_coverage[x]< coverage_threshold;x+=direction) c+=m_coverage[x]; // skip gap
    if(x==back_pos) break;
    dT("Skipping high cov region");
    for(;lower_limit<=x && x<=upper_limit && m_coverage[x]>=coverage_threshold;x+=direction) c+=m_coverage[x]; // skip next region
    double w = (x-fragment_begin)*direction;
    dT("[fragment_begin,x]=["<<fragment_begin<<","<<x<<"]"<<back_pos);
    if(w==0) break;
    dT("w="<<w<<", not zero.");
    assert(w>1);
    //if(c/w<f_threshold){x=fragment_begin; break;} // Discard the last thin region
    CoverageDistribution(front_pos-direction, front_pos-direction*(s_window_size-1), w, &m_covdist, c/w);
    dT("c/w="<<c<<"/"<<w<<"="<<c/w<<", rate*cov_thr="<<m_fragment_threshold_rate<<"*"<<coverage_threshold<<"="<<m_fragment_threshold_rate*coverage_threshold);
    if(c/w<m_fragment_threshold_rate*coverage_threshold){x=fragment_begin; break;} // Discard the last thin region
  }
  dT("x="<<x);
  return x;
}

void CCoverageArray::RefineEdge(std::vector<int64_t> *regions, int64_t front_pos, int64_t back_pos, refine_type_t rt) const{
  assert(regions);
  if(front_pos>back_pos){
    // CHECK
    Warning("Degenerage region: "<<front_pos<<">"<<back_pos);
    return;
  }
  int64_t b=front_pos;
  int64_t e=back_pos;
  if(MarginSize()>0) std::cerr<<"Warning: margine-size("<<MarginSize()<<") was ignored in initial clipping.\n";
  while(b<=back_pos  && m_coverage[b]>=m_coverage_threshold) ++b;
  while(e>=front_pos && m_coverage[e]>=m_coverage_threshold) --e;
  if(b<e){
    b=RefineOneSide(b,e,rt);
    e=RefineOneSide(e,b,rt);
    if(m_max_chop_length>0){
      if(b-front_pos>m_max_chop_length){
        std::cerr<<"# Canceling ["<<b<<",*], return to ["<<front_pos<<",*]\n";
        b=front_pos;
      }
      if(back_pos-e>m_max_chop_length){
        std::cerr<<"# Canceling [*,"<<e<<"], return to [*,"<<back_pos<<"]\n";
        e=back_pos;
      }
    }
  }
  if(front_pos<=b && b<e && e<=back_pos){
    regions->push_back(b);
    regions->push_back(e);
  }else{
    regions->push_back(front_pos);
    regions->push_back(back_pos);
  }
}


double CCoverageArray::CoverageIn(int64_t front_pos, int64_t back_pos) const {
  Assert(front_pos!=0);
  if(front_pos>back_pos){
    int64_t x=front_pos; front_pos=back_pos; back_pos=x;
  }

  // Calculate coverage in this region
  int64_t coverage_integer=0;
  for(int64_t i=front_pos;i<=back_pos;++i) coverage_integer += m_coverage[i];
  double coverage_double = coverage_integer;
  coverage_double /= back_pos-front_pos+1;
  return coverage_double;
}

void WriteRefinedBED(std::ostream& stream, int32_t chr, int64_t front_pos, int64_t back_pos, int64_t refined_front_pos=-1, int64_t refined_back_pos=-1){
  if(refined_front_pos<0) refined_front_pos=front_pos;
  if(refined_back_pos <0) refined_back_pos =back_pos;
  stream<<"chr"<<chr
        <<'\t'<<front_pos-1
        <<'\t'<<back_pos
        <<'\t'<<"refined"
        <<'\t'<<1
        <<'\t'<<'+'
        <<'\t'<<refined_front_pos-1
        <<'\t'<<refined_back_pos
    //<<'\t'<<"255,0,0"
        <<std::endl;
}

void CCoverageArray::AnalyzeRegion(std::ostream &stream, int32_t no, int64_t front_pos, int64_t back_pos, refine_type_t rt) const{
  double coverage_double = CoverageIn(front_pos,back_pos);
  Assert(front_pos<=back_pos);
  //if(coverage_double==0) return; // remove this line.
  if(coverage_double==0){
    if(!m_output_format) return;
    if(std::strcmp("bed",m_output_format)!=0) return;
    WriteRefinedBED(stream,MyChr(),front_pos,back_pos);
    return;
  }
  std::vector<int64_t> region;
  if(m_refinement_coverage_threshold>0 && coverage_double>m_refinement_coverage_threshold){
    std::cerr<<"# Skipping ["<<front_pos<<","<<back_pos<<"]: threshold="<<m_refinement_coverage_threshold<<", coverage="<<coverage_double<<std::endl;
    // No change
    region.push_back(front_pos);
    region.push_back(back_pos);
  }else{
    // Calculate split regions
    switch(rt){
    case REFINE_COVERAGE:        RefineByCoverage(&region,front_pos,back_pos); break;
    case REFINE_EDGE:            RefineEdge(&region,front_pos,back_pos,rt); break;
    case REFINE_FRAGMENTED_EDGE: RefineEdge(&region,front_pos,back_pos,rt); break;
    default:
      Quit("Unexpected refine_type: "<<rt);
    }
  }
  // Write the result
  if(! *m_output_format){
    stream<<"V\t"<<no<<'\t'<<coverage_double<<"\t["<<front_pos<<","<<back_pos+1<<")\t";
    for(uint32_t i=0;i<region.size();i+=2){
      if(i>0) stream<<' ';
      stream<<region[i]<<'-'<<region[i+1];
    }
    stream<<std::endl;
  }else if( !std::strcmp("bed",m_output_format) ){
    double rate = 0;
    if(! region.empty()){
      int32_t len=back_pos-front_pos;
      Assert(len>=0); //if(len<0) rate *= -1;
      Assert(region.size()==2);
      len+=1;
      rate = static_cast<double>(region[1]-region[0]+1)/len;
    }
    if(region.empty()){
      WriteRefinedBED(stream,MyChr(),front_pos,back_pos);
    }else{
      WriteRefinedBED(stream,MyChr(),front_pos,back_pos,region[0],region[1]);
    }
  }else{
    Quit("Unknown output format: "<<m_output_format);
  }
}

void CCoverageArray::RefineRegion(std::ostream &stream, refine_type_t rt, CGeneralFeatureVector& variants){
  if(variants.empty()) Quit("Empty variation sets");
  std::sort(variants.begin(),variants.end());
  for(uint32_t i=0;i<variants.size();++i){
    const CGeneralFeature& var = variants[i];
    AnalyzeRegion(stream,i,var.Start(),var.End(),rt);
  }
}

/*
int64_t CCoverageArray::MinPosition() const {
    int64_t min_position=0;
    while(m_coverage[min_position]==0 && min_position<=MaxPosition()) ++min_position;
    return min_position;
}
*/

void CCoverageArray::ReadSAM(const char *sam_file, CChromosomeNormalizer& normalizer){
  CSAMReader::ReadSAM(sam_file,normalizer);
  m_average_density =  m_read_positions.NReads();
  m_average_density /= MaxPosition()-MinPosition();
  if(Option().Find("verbose")){
    std::cerr<<"# max_coverage="<<m_max_coverage<<std::endl;
  }
}

void CCoverageArray::Treat(const CSAMAlignment& aln, const char *text){
  m_read_positions.AddRead(aln.Start(),aln.End());
  for(int64_t i=aln.Start();i<=aln.End();++i){
    if(m_coverage[i]==std::numeric_limits<uint16_t>::max()) Quit("Coverage overflow: i="<<i<<": "<<text);
    m_coverage[i]++;
    if(m_max_coverage<m_coverage[i]) m_max_coverage=m_coverage[i];
  }
}

void CCoverageArray::Show(std::ostream &stream, int32_t indent) const {
  stream<<"{\"chr\":"<<MyChr()<<", \"genome_size\":"<<GenomeSize()<<", ["<<std::endl;
  for(int64_t i=0;i<=GenomeSize()+1;++i){ // !!! array index
    if(i>0) stream<<","<<std::endl;
    stream<<" ["<<i<<", "<<m_coverage[i]<<"]";
  }
  stream<<'}'<<std::endl;
}

////////////////////////////////////////////////////////////////////////////////
void CCoverageArray::StatisticalAnalysis(int32_t coverage_unit) const {
  int64_t binsize = Option().RequireInteger("bin-size");
  std::cerr<<"# binsize="<<binsize<<std::endl;
  std::cerr<<"# range=["<<MinPosition()<<','<<MaxPosition()-binsize<<"]"<<std::endl;
  int32_t n_bins = (MaxPosition()-binsize-MinPosition())/binsize;
  std::cerr<<"# total "<<n_bins<<" bins."<<std::endl;
  int32_t max_frequency = binsize*1024/coverage_unit;
  int32_t frequency[max_frequency];
  int64_t total_n_bases = 0;
  for(int32_t i=0;i<max_frequency;++i) frequency[i]=0;
  for(int64_t window=MinPosition();window<MaxPosition()-binsize;window+=binsize){
    int64_t b=0;
    for(int64_t i=window;i<window+binsize;++i) b += m_coverage[i];
    total_n_bases += b;
    if(b/coverage_unit>=max_frequency){
      std::cerr<<"Warning: Too much frequency: window="<<window<<", nbases="<<b<<std::endl;
      continue;
    }
    ++frequency[b/coverage_unit];
  }

  std::cerr<<"# average coverage = "<<static_cast<double>(total_n_bases)/(MaxPosition()-MinPosition())<<std::endl;

  double cumulative_n_bins=0;
  for(int32_t i=0;i<max_frequency;++i){
    if(frequency[i]==0) continue;
    std::cout<<i<<'\t'<<static_cast<double>(i)*coverage_unit/binsize<<'\t'
             <<frequency[i]<<'\t'<<cumulative_n_bins<<'\t'<<cumulative_n_bins/n_bins
             <<std::endl;
    cumulative_n_bins += frequency[i];
  }
}

