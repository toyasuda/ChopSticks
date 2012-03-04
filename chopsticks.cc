////////////////////////////////////////////////////////////////////////////////
// sam2clippoisition
//
// 2011.4.2 T.Yasuda
////////////////////////////////////////////////////////////////////////////////

#include <iostream>
//#include <map>

#include "Option.h"
#include "Utility.h"
#include "FileReader.h"
#include "SAMAlignment.h"
#include "SequenceSet.h"
#include "CoverageArray.h"
#include "EvidenceFinder.h"
#include "Tool.h"

#include <limits>
#include <algorithm>
#include <cassert>

#include <cmath>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>

#ifdef BITVECTOR_LIB_BEGIN
using namespace BitVectorLib;
#endif

////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
COption::Option_t g_option_spec[] = {
  {"show-bam-line",            "A",0,"Flag to determine whether alignment in BAM file should be shown",0},
  {"margin-parameter",         "a",1,"The parameter for determining threshold based on coverage of margin region","0"},
  {"library-threshold",        "B",1,"File that contain threshold of discordant pairs for each library",""},
  {"bin-size",                 "b",1,"Size of bin for statistical test","1"},
  {"show-clipping",            "C",0,"Flag to determine whether all cordinates should be shown",0},
  {"coverage-upper-bound",     "d",1,"Threshold of coverage shown","500"},
  {"output-format",            "F",1,"Output format",""},
  {"fragment-threshold",       "f",1,"The threshold of joining fragmented edges","0"},
  {"genome-size",              "G",1,"Genome size","300000000"},
  {"cluster-gap-threshold",    "g",1,"Threshold of gaps between positions in clusters","3"},
  {"coverage-threshold",       "k",1,"Threshold of read coverage to fill SV region","1"},
  {"dangling-distance",        "L",1,"Threshold of dangling distance displayed","10000"},
  {"clipping-length-threshold","l",1,"Threshold of clipping length","10"},
  {"margin-size",              "M",1,"The size of margine region for calculating outside coverage","0"},
  {"show-memory-usage",        "m",0,"Show memory usage",0},
  {"output-stderr",            "O",1,"Type of output to stderr (D(angling),V(ariation)","D"},
  {"progress-interval",        "p",1,"Interval for progress report","0"},
  {"mininum-quality-symbol",   "q",1,"The symbol representing the minimum quality in SAM format","'!'"},
  {"refinement-threshold",     "R",1,"Highest coverage where refinement will be applied","-1"},
  {"max-read-length",          "r",1,"Maximum length of short reads","256"},
  {"cluster-size-threshold",   "s",1,"Threshold of cluster size","2"},
  {"verbose",                  "V",0,"Show extra messages",0},
  {"coverage-window",          "W",1,"Size and scale factor of coverage distribution","0:100:100"},
  {"max-chop-length",          "x",1,"Maximum allowed chop length","200"},
  {0,0,0,0}
};

int true_main(int argc, char *argv[]){
  // Command line arguments
  int skip = Option().SetUp(argc, argv, g_option_spec);
  int n_args = argc-skip;
  if(n_args<1){
    std::cerr<<"Usage: "<<argv[0]<<" <subcommand> ..."<<std::endl;
    COption::ShowHelp(std::cerr, g_option_spec);
    return 1;
  }

  std::string subcommand(argv[skip]);
  if(subcommand=="read_sam"){
    if(n_args<2) Quit("Usage: "<<argv[0]<<" read_feature <target file>");
    std::string target_file(argv[skip+1]);
    CFileReader fr(target_file.c_str());
    fr.SetProgressInterval(Option().RequireInteger("progress-interval"),"lines");
    while(fr.GetContentLine("@")) CSAMAlignment(fr.CurrentLine()).Show(std::cout);
  }
  else if(subcommand=="read_feature"){
    if(n_args<2) Quit("Usage: "<<argv[0]<<" read_feature <target file>");
    std::string target_file(argv[skip+1]);
    CFileReader fr(target_file.c_str());
    fr.SetProgressInterval(Option().RequireInteger("progress-interval"),"lines");
    int32_t filetype = CGeneralFeature::Filetype(target_file.c_str());
    std::cerr<<"# filetype = "<<filetype<<std::endl;
    while(fr.GetContentLine()){
      CGeneralFeature feature;
      if(! feature.Parse(fr.CurrentLine(),filetype)) continue;
      feature.WriteGFF(std::cout);
      std::cout<<std::endl;
    }
  }
  else if(subcommand=="trim"){
    if(n_args<5) Quit("Usage: "<<argv[0]<<" trim <chromosome no.> <accession table> <sam file> <gff file>");
    std::string chr_str(argv[skip+1]);
    CCoverageArray ca;
    CGeneralFeatureVector gfv;
    CChromosomeNormalizer cn;
    cn.Read(argv[skip+2]);
    ca.SetUp( cn.Chr(chr_str.c_str()) );
    //gfv.ReadGFF(std::atoi(chr_str.c_str()),cn,argv[skip+4]);
    gfv.ReadGFF( cn.Chr(chr_str.c_str()), cn, argv[skip+4]);
    ca.ReadSAM(argv[skip+3],cn);
    if(Option().Find("verbose")){
      std::cerr<<"# Unknown references: ";
      cn.ShowUnknown(std::cerr);
    }
    ca.RefineRegion(std::cout,CCoverageArray::REFINE_FRAGMENTED_EDGE,gfv);
  }
  else if(subcommand=="evidence"){
    if(n_args<5) Quit("Usage: "<<argv[0]<<" evidence <chromosome no.> <accession table> <sam file> <gff file>");
    std::string chr_str(argv[skip+1]);
    CEvidenceFinder ef;
    CGeneralFeatureVector gfv;
    CChromosomeNormalizer cn;
    cn.Read(argv[skip+2]);
    gfv.ReadGFF( cn.Chr(chr_str.c_str()), cn, argv[skip+4]);
    gfv.Sort();
    ef.SetUp( cn.Chr(chr_str.c_str()) );
    ef.MakeBin( gfv );
    ef.ReadSAM(argv[skip+3],cn);
    if(Option().Find("verbose")){
      std::cerr<<"Unknown references: ";
      cn.ShowUnknown(std::cerr);
    }
  }
  else if(subcommand=="sequence_set"){
    std::string target_file(argv[skip+1]);
    if(n_args<2) Quit("Usage: "<<argv[0]<<" sequence_set <FASTQ file>");
    CFileReader fr(target_file.c_str());
    CSequenceSet sset;
    sset.Read(fr);
    sset.ShowStat(std::cout);
  }
  else if(subcommand=="options"){
    Option().Show(std::cout);
  }
  else if(subcommand=="test"){
    CSlotManager::Test(std::cout);
  }
  else{
    std::cerr<<"# Unexpected subcommand: "<<subcommand<<std::endl;
  }

  if(Option().Find("show-memory-usage")){
    CProcessMemory pm;
    std::cerr<<"PeakMemory="<<pm.Peak<<std::endl;
  }
  return 0;
}

int main(int argc, char *argv[]){
  try{
    return true_main(argc,argv);
  }
  catch(CError &e){
    std::cerr<<e<<std::endl;
  }
}
