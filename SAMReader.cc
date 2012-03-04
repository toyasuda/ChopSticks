#include <iostream>
#include <cassert>
#include <cstdio>
#include <cstdlib>

#include "Option.h"
#include "Utility.h"
#include "SAMReader.h"
#include "SAMAlignment.h"

#ifdef BITVECTOR_LIB_BEGIN
using namespace BitVectorLib;
#endif

////////////////////////////////////////////////////////////////////////////////
void CSAMReader::Initialize(){
  m_chr=0;
  m_genome_size=0;
  m_margin_size=0;
  m_max_position=0;
  m_min_position=-1;
  m_n_total_bases=0;
}

CSAMReader::CSAMReader(){
  //m_genome_size=m_chr=0;
  //m_margin_size = Option().RequireInteger("margin-size");
  //m_max_position=0;
  //m_min_position=-1;
  //m_n_total_bases=0;
  Initialize();
  m_margin_size = Option().RequireInteger("margin-size");
}
void CSAMReader::SetUp(int32_t chrNo){
  // Get chromosome
  Initialize();
  m_chr = chrNo;
  m_genome_size = Option().RequireInteger("genome-size");
  m_margin_size = Option().RequireInteger("margin-size");
  const char *v = Option().Require("library-threshold");
  if(v && *v){
    m_library_threshold.Read(v);
  }else{
    v = "N/A";
  }
  if(Option().Find("verbose")){
    std::cerr<<"# chromosome="<<MyChr()<<std::endl;
    std::cerr<<"# genome-size="<<GenomeSize()<<std::endl;
    std::cerr<<"# margin-size="<<MarginSize()<<std::endl;
    std::cerr<<"# library-threshold="<<v<<std::endl;
  }
  //m_max_position=0;
  //m_min_position=-1;
  //m_n_total_bases=0;
}

void CSAMReader::TreatHeader(const char *text){}

//void CCoverageArray::ReadSAM(const char *sam_file, CChromosomeNormalizer& normalizer){
void CSAMReader::ReadSAM(const char *sam_file, CChromosomeNormalizer& normalizer){
  CCountLines cl("CCoverageArray::ReadSAM");
  CFileReader fr(sam_file);
  fr.SetProgressInterval(Option().RequireInteger("progress-interval"),"lines");
  int64_t n_invalid=0;
  int64_t prev_begin=0;
  int32_t prev_chr=0;
  //while(fr.GetContentLine("@")){
  while(fr.GetContentLine("")){
    if(fr.CurrentLine()[0]=='@'){
      TreatHeader(fr.CurrentLine());
      continue;
    }
    cl.IncrementAll();
    CSAMAlignment aln(fr.CurrentLine());

    int32_t chr = normalizer.Chr(aln.RName());
    if(chr==0){ std::cerr<<"Warning: Unknown reference name : "<<aln.QName()<<std::endl; continue;}
    if(prev_chr>chr) Quit("SAM alignment are not sorted by chromosome: "<<prev_chr<<">"<<chr);
    prev_chr=chr;
    if(chr<MyChr()) continue;
    if(chr>MyChr()) return;

    if(aln.Unmapped()) continue;
    if(aln.Start()==0 && aln.End()==-1){n_invalid++; continue;}
    if(aln.Start()<1 || aln.Start()>aln.End() || aln.End()>=GenomeSize()){
      std::cerr<<"Warning: Alignment start="<<aln.Start()<<", while end="<<aln.End()<<std::endl;
      continue;
    }

    if(prev_begin>aln.Start()) Quit("SAM alignments are not sorted by position: "<<prev_begin<<">"<<aln.Start());
    prev_begin=aln.Start();
    //m_read_positions.AddRead(aln.Start(),aln.End());

    if(normalizer.Chr(aln.RName()) != MyChr()) continue;
    if(aln.End()>=GenomeSize()) Quit("Going beyond the end of genome ("<<GenomeSize()<<"bp): "<<aln);
    //m_n_total_bases += aln.End()-aln.Start()+1;
    if(aln.End()>aln.Start()){
      //if(aln.Start()<10) std::cerr<<fr.CurrentLine()<<std::endl;
      AddTotalBases(aln.End()-aln.Start()+1);
      UpdateMinPosition(aln.Start());
      UpdateMaxPosition(aln.End());
    }
    cl.IncrementChr();

    Treat(aln, fr.CurrentLine());
  }
  if(n_invalid){
    std::cerr<<"Warning: Number of invalid alignments = "<<n_invalid<<std::endl;
  }

  //void CCoverageArray::ReadSAM(const char *sam_file, CChromosomeNormalizer& normalizer){
  if(Option().Find("verbose")){
    //int64_t min_position = MinPosition();
    //std::cerr<<"# max_coverage="<<m_max_coverage<<std::endl;
    std::cerr<<"# max_position="<<MaxPosition()<<std::endl;
    //std::cerr<<"# min_position="<<min_position<<std::endl;
    std::cerr<<"# min_position="<<MinPosition()<<std::endl;
  }
}
