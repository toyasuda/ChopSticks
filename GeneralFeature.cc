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

/*
int32_t CCoverageDistribution::s_scale_factor=-1;
int32_t CCoverageDistribution::s_index_limit =-1;
int32_t CCoverageDistribution::s_window_size =-1;

int32_t CCoverageArray::s_window_size =-1;
*/

////////////////////////////////////////////////////////////////////////////////

int32_t CChromosomeNormalizer::Chr(const char *refname){
  //if(!refname) Quit("Unexpected reference name: "<<refname);
  if(!refname) Quit("Reference name cannot be null");
  const char *p=refname;
  while(isdigit(*p)) ++p;
  if(!*p) return std::atoi(refname);

  const char *normalized = m_kvs.Find(refname);
  if(!normalized){++ m_unknown_refname[refname]; return 0;}

  char *q=0;
  int32_t chr = std::strtol(normalized,&q,10);
  if(*q) Quit("Bad chromosome number "<<normalized<<" for "<<refname);
  return chr;
}



////////////////////////////////////////////////////////////////////////////////

//void CCoverageArray::ReadGFF(const char *gff_file){
//void CGeneralFeatureVector::ReadGFF(int32_t chr, CCoverageArray& ca, const char *gff_file){
void CGeneralFeatureVector::ReadGFF(int32_t chr, class CChromosomeNormalizer& cn, const char *gff_file){
  int32_t filetype = CGeneralFeature::Filetype(gff_file);
  std::cerr<<"# Filetype="<<filetype<<std::endl;
  CCountLines cl("CCoverageArray::ReadGFF");
  CFileReader fr(gff_file);
  fr.SetProgressInterval(Option().RequireInteger("progress-interval"),"lines");
  while(fr.GetContentLine("#")){
    cl.IncrementAll();
    CGeneralFeature gff;
    if(! gff.Parse(fr.CurrentLine(),filetype)) continue;
    if(gff.Start()<1 || gff.End()<1 || gff.Start()>gff.End()){
      gff.WriteGFF(std::cerr);
      std::cerr<<std::endl;
      Quit("Invalid GFF line: "<<fr.GetContentLine());
    }
    if(gff.Start()==gff.End()) std::cerr<<"Range size 1: "<<fr.GetContentLine()<<std::endl;

    //if(Chr(gff.Name()) != m_chr) continue;
    if(cn.Chr(gff.Name()) != chr) continue;
    //m_dbvar_variants.push_back(gff);
    push_back(gff);
    cl.IncrementChr();
  }
}

void CGeneralFeatureVector::Sort(){
  std::sort(begin(),end());
}

////////////////////////////////////////////////////////////////////////////////
// GFF: General feature format
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

