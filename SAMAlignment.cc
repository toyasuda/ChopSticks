/**
 * @file    SAMAlignment.h
 * @brief   Parse alignments in SAM format, and represent it.
 *
 * @author  Tomohiro Yasuda
 * @date    2011-4-15
 *
 */

#include <algorithm>
#include <limits>
#include <cstring>

#include "Option.h"
#include "Utility.h"
#include "FileReader.h"
#include "SAMAlignment.h"

using namespace BitVectorLib;

////////////////////////////////////////////////////////////////////////////////

int32_t CCigarString::s_clip_length_threshold = -1;

CCigarString::CCigarString(const CCigarString& cs){
  m_editscript = cs.m_editscript;
  m_length = cs.m_length;
  m_alignment_length = cs.m_alignment_length;
  m_flag = cs.m_flag;
}

int64_t CCigarString::ClipPositionRight(int64_t p) const {
  if((m_flag&2)==0) return -1;
  p += m_length-m_editscript.back().Len;
  if( m_editscript.front().Sym=='S' ) p -= m_editscript.front().Len;
  return p;
}

CCigarString::CCigarString(){
  m_alignment_length=m_length=m_flag=0;
  if(s_clip_length_threshold<0){
    //const char *threshold_str = Option().Find("clipping-length-threshold");
    //s_clip_length_threshold = threshold_str? std::atoi( threshold_str ): 10;
    s_clip_length_threshold = Option().RequireInteger("clipping-length-threshold");
  }
}

int32_t CCigarString::__IsClipped(){
  if(m_editscript.size()<1) return 0;
  //int32_t flag=0;
  m_flag=0;
  if(m_editscript.front().IsClip()) m_flag+=1;
  if(m_editscript.back(). IsClip()) m_flag+=2;
  return m_flag;
}

bool CCigarString::Parse(const char *cigar_string){
  m_editscript.clear();
  if(!cigar_string) Quit("Empty CIGAR string");
  const char *c0=cigar_string;
  m_length=0;
  m_alignment_length=0;
  while(*cigar_string){
    int32_t len=0;
    for(;isdigit(*cigar_string);++cigar_string) len=len*10+*cigar_string-'0';
    if(*cigar_string<'A' || 'Z'<*cigar_string) Quit("Invalid CIGAR string: "<<c0);
    char sym=*cigar_string++;
    m_editscript.push_back(operation_t());
    m_editscript.back().Len=len;
    m_editscript.back().Sym=sym;
    if(sym!='D'){
      m_length+=len;
      if(sym!='S' && sym!='H') m_alignment_length+=len;
    }
  }
  __IsClipped();
  return true;
}

void CCigarString::Write(std::ostream &stream, int32_t indent) const {
  stream<<"[";
  bool is_first=true;
  foreach_const(std::vector<operation_t>,it,m_editscript){
    if(is_first){is_first=false;}else{stream<<' ';}
    (*it).Show(stream,indent+1);
  }
  stream<<"]";
  //stream<<std::endl;
}

////////////////////////////////////////////////////////////////////////////////
std::map<std::string,uint32_t>* CGenomePosition::s_chromosome2int=0;

uint64_t CGenomePosition::Chromosome(const char *ch){
  //const char *ch0=ch;
  int32_t order=0;
  if(!s_chromosome2int) s_chromosome2int = new std::map<std::string,uint32_t>;
  //while('a'<*ch && *ch<'u' && *ch!='m') ++ch;
  if(std::strncmp(ch,"chromosome",10)==0) ch+=10;
  if(std::strncmp(ch,"chr",       3) ==0) ch+=3;
  if(std::strncmp(ch,"ch",        2) ==0) ch+=2;
  if( isdigit(*ch) ){
    order = strtol(ch,0,10);
  }else{
    switch(*ch & 0xdf){
    case 'X': order=23; ++ch; break;
    case 'Y': order=24; ++ch; break;
    case 'U': order=26; ++ch; break;
    case 'M': if((ch[1]&0xdf)=='T'){ order=25; ch+=2; break; }
    default:{
      std::map<std::string,uint32_t>::iterator it=s_chromosome2int->find(ch);
      if(it==s_chromosome2int->end()){
        order=s_chromosome2int->size();
        (*s_chromosome2int)[ch]=order;
        std::cerr<<"Warning: Unknown chromosome name: "<<ch<<" -> "<<order+100<<std::endl;
      }else{
        order=(*it).second;
      }
      order+=100;
    }
    }
  }
  return order;
}


int32_t CClipPositions::s_cluster_size_threshold = -1;
int32_t CClipPositions::s_cluster_gap_threshold  = -1;

CClipPositions::CClipPositions(const char *p){
  if(!p) Quit("Empty prefix for output");
  m_cluster_ID=0;
  m_prefix=p;
  if(s_cluster_size_threshold<0){
    s_cluster_size_threshold = Option().RequireInteger("cluster-size-threshold");
    s_cluster_gap_threshold  = Option().RequireInteger("cluster-gap-threshold");
  }
}

void CClipPositions::MakeOneCluster(std::vector<int64_t>::const_iterator b, std::vector<int64_t>::const_iterator e){
  std::ostream& stream=std::cout;

  if(e-b<s_cluster_size_threshold) return;
  ++m_cluster_ID;
  stream<<m_prefix<<'\t'<<e-b<<'\t'
        <<CGenomePosition::Chromosome(*b)<<'\t'
        <<CGenomePosition::Position(*b)<<'\t'
        <<CGenomePosition::Position(*(e-1))<<"\t[";
  for(std::vector<int64_t>::const_iterator it=b;it!=e;++it){
    if(it!=b) stream<<',';
    stream<<CGenomePosition::Position(*it);
  }
  stream<<']'<<std::endl;
}

void CClipPositions::FindClusters(){
  if(m_positions.empty()) return;
  sort(m_positions.begin(),m_positions.end());
  std::vector<int64_t>::const_iterator cluster_start=m_positions.begin();
  int64_t prev=*cluster_start;
  for(std::vector<int64_t>::const_iterator iter=cluster_start+1;iter!=m_positions.end();++iter){
    if(*iter-prev>=s_cluster_gap_threshold){
      MakeOneCluster(cluster_start,iter);
      cluster_start=iter;
    }
    prev=*iter;
  }
  MakeOneCluster(cluster_start,m_positions.end());
}

////////////////////////////////////////////////////////////////////////////////
// SAM Alignment

void CSAMAlignment::Copy(const CSAMAlignment& sa){
  m_chr = sa.m_chr;
  m_qname = sa.m_qname;
  m_flag = sa.m_flag;
  m_rname = sa.m_rname;
  m_pos = sa.m_pos;
  m_mapq = sa.m_mapq;
  m_cigar = sa.m_cigar;
  m_rnext = sa.m_rnext;
  m_pnext = sa.m_pnext;
  m_tlen = sa.m_tlen;
  m_seq = sa.m_seq;
  m_qual = sa.m_qual;
  m_cigar_info = sa.m_cigar_info;
}

bool CSAMAlignment::operator<(const CSAMAlignment& sa) const {
  if(m_chr<sa.m_chr) return true;
  if(m_chr>sa.m_chr) return false;
  if(m_pos<sa.m_pos) return true;
  if(m_pos>sa.m_pos) return false;
  if(m_tlen<sa.m_tlen) return true;
  if(m_tlen>sa.m_tlen) return false;
  if(m_qname<sa.m_qname) return true;
  if(m_qname>sa.m_qname) return false;
  Quit("Cannot determine order: "<<*this<<", "<<sa);
}

void CSAMAlignment::Write(std::ostream &stream, int32_t indent) const {
  const char *delimiter = ",";
  stream<<m_qname<<delimiter
        <<m_flag <<delimiter
        <<m_rname<<delimiter
        <<m_pos  <<delimiter
        <<m_mapq <<delimiter
        <<m_cigar<<delimiter
        <<m_rnext<<delimiter
        <<m_pnext<<delimiter
        <<m_tlen <<delimiter
        <<"\t"
        <<m_seq  <<delimiter
        <<m_qual;
}

const char* CSAMAlignment::Option(const std::string& k) const {
  str2str_t::const_iterator iter=m_options.find(k);
  if(iter==m_options.end()) return 0;
  return (*iter).second.c_str();
}

bool CSAMAlignment::Parse(const char** pp){
  CTokenizer td(*pp,"\t");
  static const char *field_name[] =
    {"(no field)",
     "QNAME","FLAG","RNAME","POS","MAPQ","CIGAR","RNEXT","PNEXT","TLEN","SEQ","QUAL",0};

  //const char *p="";
  m_options.clear();
  std::string num;
  try {
    m_qname = td.NextString();
    m_flag  = td.NextInteger();
    m_rname = td.NextString();
    num=td.NextString();
    m_pos   = num=="*"? 0: std::atoll(num.c_str());
    num=td.NextString();
    m_mapq  = num=="*"? 0: std::atoll(num.c_str());
    m_cigar = td.NextString();
    m_rnext = td.NextString();
    num=td.NextString();
    m_pnext = num=="*"? 0: std::atoll(num.c_str());
    num=td.NextString();
    m_tlen  = num=="*"? 0: std::atoll(num.c_str());
    m_seq   = td.NextString();
    m_qual  = td.NextString();
    if((m_flag&0x04)!=0 || m_rname=="*" || m_pos==0){
      // Unmapped
      m_flag &= 0xffff-(0x0002+0x0010+0x0100);
      m_rname.erase();
      m_pos=0;
      m_cigar.erase();
      m_mapq=255;
    }
    if((m_flag&0x08)!=0 || m_rnext=="*" || m_pnext==0){
      m_flag &= 0xffff-(0x0020);
      m_rnext.erase();
      m_pnext=0;
    }
    if(m_qname=="*") m_qname.erase();
    if(m_seq  =="*") m_seq.erase();
    if(m_qual =="*") m_qual.erase();
    if(m_cigar!="*") m_cigar_info.Parse(m_cigar.c_str());

    while(td.hasNext()){
      std::string option=td.NextString();
      uint32_t i=0;
      while(i<option.length() && option[i]!=':') ++i;
      if(i==0) Quit("Unexpected option value: "<<option);
      if(i!=2) Warning("Length of option name is not 2: "<<option);
      if(i>=option.length()) Quit("Option does not have value: "<<option);
      uint32_t j=i+1;
      while(j<option.length() && option[j]!=':') ++j;
      if(j>=option.length()) Quit("Option does not have value: "<<option);
      std::string k(option,0,i);
      std::string v(option,j+1);
      m_options[k]=v;
      //std::cerr<<"option='"<<option<<"', k='"<<k<<"', v='"<<v<<"'"<<std::endl;
    }
  }catch(CError &e){
    uint32_t fID = td.FieldID();
    Assert(fID<nelems(field_name));
    Quit("Error in mandatory field"<<fID<<" "<<field_name[fID]<<": \""<<td.Prev()<<"\""<<", ("<<e<<")");
  }
  *pp=td.Peep();
  return true;
}

////////////////////////////////////////////////////////////////////////////////
// GFF: General feature format
/*
void CGeneralFeature::WriteGFF(std::ostream &stream, int32_t indent) const {
  const char *delimiter = ",";
  const char *strand_str[] = {"-","?","+"};
  stream<<m_name  <<delimiter
        <<m_source<<delimiter
        <<m_type  <<delimiter
        <<m_start <<delimiter
        <<m_end   <<delimiter
        <<m_score <<delimiter
        <<strand_str[m_strand+1]<<delimiter
        <<m_frame <<delimiter
        <<m_group;
}
void CGeneralFeature::WriteBED(std::ostream &stream, int32_t indent) const {
  const char *delimiter = ",";
  const char *strand_str[] = {"-","?","+"};
  stream<<m_name  <<delimiter
        <<m_start-1<<delimiter
        <<m_end   <<delimiter
        <<m_type  <<delimiter
        <<m_score <<delimiter
        <<strand_str[m_strand+1]<<delimiter
    //<<m_frame <<delimiter
    //<<m_group;
    ;
}

bool CGeneralFeature::ParseGFF(const char** pp){
  if(check_prefix("browser",*pp)) return false;
  if(check_prefix("track",  *pp)) return false;

  //std::cerr<<"ParseGFF: "<<*pp<<std::endl;
  CTokenizer td(*pp,"\t");
  Clear();
  static const char *field_name[] =
    {"(no field)",
     "NAME","SOURCE","TYPE","START","END","SCORE","STRAND","FRAME","GROUP",0};

  //const char *p="";
  try {
    std::string num;
    m_name  = td.NextString();
    if(m_name==".") m_name.erase();
    m_source= td.NextString();
    if(m_source==".") m_source.erase();
    m_type  = td.NextString();
    if(m_type==".") m_source.erase();
    num=td.NextString();
    m_start = num=="."? 0: std::atoll(num.c_str());
    num=td.NextString();
    m_end   = num=="."? 0: std::atoll(num.c_str());
    num=td.NextString();
    m_score = num=="."? 0: std::atoll(num.c_str());
    num=td.NextString();
    m_strand = 0;
    if(! num.empty()){
      if     (num[0]=='+' || num[0]=='f' || num[0]=='F') m_strand= 1;
      else if(num[0]=='-' || num[0]=='r' || num[0]=='R') m_strand=-1;
    }
    m_frame = td.NextString();
    if(m_frame==".") m_frame.erase();
    m_group = td.NextString();
    if(m_group==".") m_group.erase();
  }catch(CError &e){
    uint32_t fID = td.FieldID();
    Assert(fID<nelems(field_name));
    Quit("Error in mandatory field"<<fID<<" "<<field_name[fID]<<": \""<<td.Prev()<<"\""<<", ("<<e<<")");
  }
  *pp=td.Peep();
  return true;
}

void CGeneralFeature::Clear(){
  m_start=m_end=m_score=m_strand=0;
  m_name.erase();
  m_source.erase();
  m_type.erase();
  m_frame.erase();
  m_group.erase();
}

bool CGeneralFeature::ParseBED(const char** pp){
  if(check_prefix("browser",*pp)) return false;
  if(check_prefix("track",  *pp)) return false;

  //std::cerr<<"ParseBED: "<<*pp<<std::endl;
  Clear();
  CTokenizer td(*pp,"\t");
  static const char *field_name[] =
    {"(no field)",
     "CHROM","START","END","TYPE","SCORE","STRAND",0};

  //const char *p="";
  try {
    std::string num;
    m_name  = td.NextString();
    if(m_name==".") m_name.erase();

    num=td.NextString();
    m_start = num=="."? 0: std::atoll(num.c_str());

    num=td.NextString();
    m_end   = num=="."? 0: std::atoll(num.c_str());

    if(td.hasNext()){
      m_type = td.NextString();

      if(td.hasNext()){
          num=td.NextString();
          m_score = std::atof(num.c_str());

          if(td.hasNext()){
            num=td.NextString();
            m_strand = 0;
            if(! num.empty()){
              if     (num[0]=='+' || num[0]=='f' || num[0]=='F') m_strand= 1;
              else if(num[0]=='-' || num[0]=='r' || num[0]=='R') m_strand=-1;
            }
          }
      }
    }

    //if(m_strand>=0){ ++m_start; }
    if(m_strand<0) std::cerr<<"Warning: feature on reverse strand: "<<*pp<<std::endl;
    if(m_strand>=0) ++m_start;
    //std::cerr<<"ParseBED: ";
    //WriteGFF(std::cerr);
    //std::cerr<<std::endl;
  }catch(CError &e){
    uint32_t fID = td.FieldID();
    Assert(fID<nelems(field_name));
    Quit("Error in mandatory field"<<fID<<" "<<field_name[fID]<<": \""<<td.Prev()<<"\""<<", ("<<e<<")");
  }
  *pp=td.Peep();
  return true;
}

bool CGeneralFeature::Parse(const char* p, int32_t filetype){
  //std::cerr<<"CGeneralFeature::Parse(p, "<<filetype<<")\n";
  switch(filetype){
  case GFF:
    //std::cerr<<"CGeneralFeature::Parse:GFF: filetype="<<filetype<<std::endl;
    return ParseGFF(p);
  case BED:
    //std::cerr<<"CGeneralFeature::Parse:BED: filetype="<<filetype<<std::endl;
    return ParseBED(p);
  default:
    Quit("Unknown file type: "<<filetype);
  }
}

// Returns 0:unknown, 1:gff, 2:bed
int32_t CGeneralFeature::Filetype(const char *filename){
  const char *q = filename;
  while(*q) ++q;
  while(filename<q && !isalnum(*(q-1))) --q;
  const char *ext = q;
  while(filename<ext && *ext!='.') --ext;
  //std::cerr<<"1: file="<<filename<<", q="<<q<<", ext="<<ext<<std::endl;
  if(*ext!='.') return UNKNOWN;
  if(ext[1]=='g' && ext[2]=='z' && !isalnum(ext[3])){
    if(filename==ext) return UNKNOWN;
    --ext;
    while(filename<ext && *ext!='.') --ext;
    if(*ext!='.') return UNKNOWN;
  }
  //std::cerr<<"2: file="<<filename<<", q="<<q<<", ext="<<ext<<std::endl;
  if(ext[1]=='g' && ext[2]=='f' && ext[3]=='f' && !isalnum(ext[4])) return GFF;
  if(ext[1]=='b' && ext[2]=='e' && ext[3]=='d' && !isalnum(ext[4])) return BED;
  return UNKNOWN;
}

bool CGeneralFeature::operator<(const CGeneralFeature &a) const {
  if(m_name<a.m_name) return true;
  if(m_name>a.m_name) return false;
  if(m_start<a.m_start) return true;
  if(m_start>a.m_start) return false;
  if(m_end<a.m_end) return true;
  if(m_end>a.m_end) return false;
  if(m_strand>a.m_strand) return true;
  if(m_strand<a.m_strand) return false;
  if(m_type<a.m_type) return true;
  if(m_type>a.m_type) return false;
  if(m_source<a.m_source) return true;
  if(m_source>a.m_source) return false;
  if(m_frame<a.m_frame) return true;
  if(m_frame>a.m_frame) return false;
  return false; // m_group is ignored
}
*/

////////////////////////////////////////////////////////////////////////////////
// Quality values
double CQualityAnalyzer::Average(const char *q) const {
  Assert(q && *q);
  double sum=0;
  int32_t len=0;
  for(len=0;q[len];++len) sum+=q[len]-m_min_quality_symbol;
  return sum/len;
}
int32_t CQualityAnalyzer::Min(const char *q) const {
  Assert(q && *q);
  int32_t x=GREATER_THAN_ANY_QUALITY;
  for(;*q;++q){
    int32_t v=*q-m_min_quality_symbol;
    if(x>v) x=v;
  }
  return x;
}
int32_t CQualityAnalyzer::Max(const char *q) const {
  Assert(q && *q);
  int32_t x=0;
  for(;*q;++q){
    int32_t v=*q-m_min_quality_symbol;
    if(x<v) x=v;
  }
  return x;
}
