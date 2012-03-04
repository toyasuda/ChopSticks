/**
 * @file    FileReader.cc
 * @brief   Read tab, comma separated text file
 *
 * @author  Tomohiro Yasuda
 * @date    2010-11-20
 *
 */

//#include <algorithm>
#include <cstdio>
#include <cstring>
#include <cassert>

#include "Utility.h"
#include "Option.h"
#include "FileReader.h"

CProgressReport::CProgressReport(){
  m_previous_turn=0;
  m_progress_interval=-1;
}
void CProgressReport::SetInterval(int64_t interval, const char *unit){
  m_progress_interval=interval;
  m_unit=unit;
}
void CProgressReport::ShowProgress(int64_t progress){
  if(m_progress_interval>0){
    //uint64_t current_turn = m_read_amount/m_progress_interval;
    uint64_t current_turn = progress/m_progress_interval;
    if(m_previous_turn < current_turn){
      //std::cerr<<current_turn*m_progress_interval<<" bytes."<<std::endl;
      std::cerr<<current_turn*m_progress_interval<<' '<<m_unit<<'.'<<std::endl;
    }
    m_previous_turn=current_turn;
  }
}


void CFileReader::Initialize(){
  m_fp=m_my_fp=m_pipe_fp=0;
  m_buffer=0;
  m_read_amount=0;
  m_erase_newline=true;
}

FILE *CFileReader::OpenPipe(const char *filename, const char *suffix, const char *cmd){
  // Calculate length
  int len=0;
  while(filename[len]) ++len;
  while(len>0 && isspace(filename[len-1])) --len;

  // Check suffix
  int suffix_len=0;
  Assert(suffix);
  while(suffix[suffix_len]) ++suffix_len;
  if(suffix_len>=len) return 0;
  //std::cerr<<"Comparing '"<<filename+len-suffix_len<<"' and "<<suffix<<std::endl;
  if(std::strncmp(filename+len-suffix_len,suffix,suffix_len)) return 0;

  // Build command line
  char command_line[len+suffix_len+10];
  int j=0;
  for(;cmd[j];++j) command_line[j]=cmd[j];
  if(*cmd) command_line[j++]=' ';
  for(int i=0;i<len;++i) command_line[j++]=filename[i];
  command_line[j]=0;
  if(!*cmd && suffix[0]=='|' && !suffix[1]) command_line[j-1]=0;
  if(Option().Find("verbose")) std::cerr<<"# Opening pipe with '"<<command_line<<"'"<<std::endl;
  m_pipe_fp = popen(command_line,"r");
  return m_pipe_fp;
}

void CFileReader::Open(const char *filename){
  if(m_fp || m_my_fp || m_buffer) Quit("Duplicatedly opening '"<<filename<<"' using instance for '"<<m_filename);
  if(filename[0]=='-' && !filename[1]){
    m_fp = stdin;
    if(!m_fp) Quit("Cannot read stdard input");
  }else{
    //std::cerr<<"void CFileReader::Open("<<filename<<")"<<std::endl;
    /*
    int len=0;
    while(filename[len]) ++len;
    if(len==0) Quit("Opening filename is empty");
    if(len>1 && filename[len-1]=='|'){
      // Opening pipe
      char argv_body[len+1];
      for(int i=0;i<=len;++i) argv_body[i]=filename[i];
      argv_body[--len]=0;
      m_fp = m_pipe_fp = popen(argv_body,"r");
      if(!m_fp) Quit("Cannot open '"<<filename<<"'");
    }
    else if(len>3 && filename[len-3]=='.' && filename[len-2]=='g' && filename[len-1]=='z'){
      // Opening pipe using zcat
      static const char *zcat = "zcat ";
      char argv_body[len+10];
      int j=0;
      for(j=0;zcat[j];++j) argv_body[j]=zcat[j];
      for(int i=0;filename[i];++i,++j) argv_body[j]=filename[i];
      argv_body[j]=0;
      std::cerr<<"Opening '"<<argv_body<<"'"<<std::endl;
      m_fp = m_pipe_fp = popen(argv_body,"r");
      if(!m_fp) Quit("Cannot open '"<<filename<<"'");
    */
    if( (m_fp=OpenPipe(filename,"|",""))==0 &&
        (m_fp=OpenPipe(filename,".gz","zcat"))==0 &&
        (m_fp=OpenPipe(filename,".bam","samtools view -h"))==0){
      // Normal file
      m_fp = m_my_fp = std::fopen(filename,"r");
      if(!m_fp) Quit("Cannot open "<<filename);
    }
  }
  m_filename=filename;
  m_buffer = new char[BUFFER_SIZE];
  m_line_number=0;
  m_read_amount=0;
}

void CFileReader::Close(){
  if(m_my_fp)  { fclose(m_my_fp);   m_my_fp=0; }
  if(m_pipe_fp){ pclose(m_pipe_fp); m_pipe_fp=0; }
  m_fp=0;
  if(m_buffer){ delete [] m_buffer; m_buffer=0; }
}

const char *CFileReader::GetLine(){
  //assert(m_fp);
  //assert(m_buffer);
  if(!m_fp || !m_buffer) return 0;
  if(fgets(m_buffer,BUFFER_SIZE,m_fp)){
    ++m_line_number;
    //const char *p=m_buffer;
    //while(*p){++p; ++m_read_amount;}
    //m_progress_reporter.ShowProgress(m_read_amount);
    m_progress_reporter.ShowProgress(m_line_number);
    if(m_erase_newline){
      char *p=m_buffer;
      for(;*p;++p);
      while(m_buffer<p && 0<*(p-1) && *(p-1)<' ') --p;
      *p=0;
    }
    return m_buffer;
  }
  return 0;
}

const char *CFileReader::GetContentLine(const char *comment_sym){
  //if(!comment_sym) comment_sym="#";
  for(;;){
    const char *buffer = GetLine();
    if(!buffer) return 0;
    //if(buffer[0]=='#' || buffer[0]=='\n' || buffer[0]==0) continue;
    if(buffer[0]=='\n' || buffer[0]==0) continue;
    if(comment_sym && comment_sym[0]){
      if(strchr(comment_sym,buffer[0])) continue;
    }
    return buffer;
  }
  return 0; // never reached
}

bool CFileReader::GetContentLine(std::vector<std::string> *str_vec){
  assert(str_vec);
  str_vec->clear();
  const char *buffer = GetContentLine();
  if(!buffer) return false;
  while(*buffer){
    const char *start = buffer;
    //while(*buffer && *buffer!='\t') ++buffer;
    while(*buffer && *buffer!='\t' && *buffer!='\n') ++buffer;
    str_vec->push_back(std::string());
    str_vec->back().assign(start,buffer-start);
    if(*buffer=='\t' || *buffer=='\n') ++buffer;
  }
  return true;
}

void CFileReader::Chomp(){
  int len=0;
  while(m_buffer[len]) ++len;
  for(;;){
    if(len<=0) break;
    --len;
    if(m_buffer[len]!='\n' && m_buffer[len]!='\r') break;
    m_buffer[len]=0;
  }
  return;
}

const char *CTokenizer::Next(const char *d){
  if(!d) d=m_delimiter;
  // Skip delimiters
  /*
  for(;*m_source_ptr; ++m_source_ptr){
    if(! std::strchr(d,*m_csr)) break;
  }
  */
  m_token.erase();
  if(!m_source_ptr || ! *m_source_ptr) return 0;
  ++m_field_id;
  m_prev=m_source_ptr;
  for(;*m_source_ptr; ++m_source_ptr){
    if(std::strchr(d,*m_source_ptr)) break;
    m_token+=*m_source_ptr;
  }
  if(*m_source_ptr) ++m_source_ptr; // skip delimiter at the right to the token
  //return m_token.c_str();
  return CStr();
}

/*
int64_t CTokenizer::NextInteger(const char *d){
  //const char *p = Next(d);
  //if(!p) Quit("Missing integer data");
  Next(d);
  return Integer();
}
*/

int64_t CTokenizer::Integer() const {
  const char *p=CStr();
  if(!p || ! *p) Quit("Missing integer data");
  if(! isdigit(*p) && *p!='-') Quit("Integer expected: "<<*p);
  //return std::atoll(p);
  char *last=0;
  int64_t retval = std::strtoll(p,&last,10);
  if(*last || sign_cast<uint32_t>(last-CStr())>m_token.length()) Quit("Extra character(s) after number");
  return retval;
}

////////////////////////////////////////////////////////////////////////////////
// Key-value store

bool CKVStore::Read(const char *target_file){
  CFileReader fr(target_file);
  fr.SetProgressInterval(Option().RequireInteger("progress-interval"),"lines");
  for(;;){
    const char *p=fr.GetContentLine();
    if(!p) break;
    Parse(&p);
  }
  return true;
}

bool CKVStore::Add(const char *key, const char *value){
  mapping_t::const_iterator iter = m_map.find(key);
  if(iter!=m_map.end()){
    if((*iter).second==value) return true;
    Quit(key<<" is associated with "<<(*iter).second<<", conflicting with "<<value);
    //return false;
  }
  m_map[key]=value;
  return true;
}

bool CKVStore::Parse(const char **pp){
  CTokenizer td(*pp,"\t");
  //const char *p=*pp;
  std::string key,value;
  key   = td.NextString();
  value = td.NextString();
  //std::cerr<<"CKVStore::Parse: key='"<<key<<"', value='"<<value<<"', p='"<<p<<"'\n";
  //m_map[key]=value;
  Add(key.c_str(),value.c_str());
  *pp=td.Peep();
  return true;
}

void CKVStore::Write(std::ostream &stream, int32_t indent) const {
  foreach_const(mapping_t,iter,m_map){
    stream<<(iter==m_map.begin()? "{": ",\n ");
    stream<<"\""<<(*iter).first<<"\": \""<<(*iter).second<<"\"";
  }
  stream<<'}';
}

const char *CKVStore::Find(const char *key) const {
  mapping_t::const_iterator iter = m_map.find(key);
  if(iter==m_map.end()) return 0;
  return (*iter).second.c_str();
}
