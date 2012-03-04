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

int32_t CCoverageDistribution::s_scale_factor=-1;
int32_t CCoverageDistribution::s_index_limit =-1;
int32_t CCoverageDistribution::s_window_size =-1;

//int32_t CCoverageArray::s_window_size =-1;


////////////////////////////////////////////////////////////////////////////////
void CCoverageDistribution::Initialize(){
  if(s_scale_factor<=0){
    char *opt = const_cast<char*>( Option().Require("coverage-window") );
    s_window_size=std::strtol(opt,&opt,10);
    if(*opt!=':' || !isdigit(*(opt+1)) ) Quit("Specify <window size>:<scale factor>:<index_limit> by --coverage-window option");
    ++opt;
    s_scale_factor=std::strtol(opt,&opt,10);
    if(*opt!=':' || !isdigit(*(opt+1)) ) Quit("Specify <window size>:<scale factor>:<index_limit> by --coverage-window option");
    ++opt;
    s_index_limit=std::strtol(opt,&opt,10);
    std::cerr<<"# CCoverageDistribution::Initialize: window size="<<s_window_size<<", scale factor="<<s_scale_factor<<", index limit="<<s_index_limit<<std::endl;
    if(s_scale_factor<=0) Quit("Scale factor must be greater than 0");
    if(s_index_limit <=0) Quit("Index limit must be greater than 0");
  }

  if(sign_cast<int32_t>(m_distribution.size())!=s_index_limit){
    m_distribution.assign(s_index_limit,0);
  }
  else if(m_n_values!=0){
    foreach(CCoverageDistribution, iter, m_distribution) *iter=0;
  }
  m_sum=m_sum_squared=0;
  m_n_values=0;
  m_n_larger=0;
}
void CCoverageDistribution::Increment(int32_t c){
  Assert(s_index_limit==sign_cast<int32_t>(m_distribution.size()));
  m_sum += c;
  m_sum_squared += c*c;
  ++m_n_values;
  int32_t i = c*s_scale_factor;
  Assert(i>=0);
  if(i>=s_index_limit){
    //std::cerr<<"# Too large coverage value: "<<c<<std::endl;
    ++m_n_larger;
    i=s_index_limit-1;
  }
  //if(! (0<=i && i<s_index_limit)) std::cerr<<"i="<<i<<", limit="<<s_index_limit<<", c="<<c<<", factor="<<s_scale_factor<<std::endl;
  ++( m_distribution[i] );
}

void CCoverageDistribution::Write(std::ostream& stream) const{
  stream<<Average()<<'\t';
  stream<<Variance()<<'\t';
  stream<<std::sqrt(Variance())<<'\t';
  stream<<m_n_values<<'\t';
  stream<<m_n_larger<<'\t';
  foreach_const(std::vector<int32_t>, iter, m_distribution){
    if(iter!=m_distribution.begin()) stream<<',';
    stream<<*iter;
  }
}

