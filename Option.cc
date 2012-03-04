/**
 * @file    Option.cc
 * @brief   Parse options with getopt(3) and keep the values in map<string,const char*>
 *
 * @author  Tomohiro Yasuda
 * @date    2010-9-20
 *
 */

// (C) Yasuda, Tomohiro: The university of Tokyo

//#include <iostream>
//#include <map>
#include <cstdlib>
//#include <cstdio>
//#include <string>
#include <getopt.h>
#include <cstring>

#include "Utility.h"
#include "Option.h"

COption *COption::s_option=0;
COption& Option(){
  if(! COption::s_option) COption::s_option = new COption;
  return *COption::s_option;
}

/**
 * @brief desctructor
 *
 * Delete each value of options
 */
COption::~COption(){
  for(std::map<std::string,const char*>::iterator iter= m_values.begin(); iter!=m_values.end(); ++iter){
    delete [] (*iter).second;
  }
}

/**
 * @brief Add a value of specified option
 *
 * Add (key,value) pair of an option.
 * @param k Long name of an option.
 * @param v Value of the option.
 */
void COption::Add(const char *k, const char *v){
  int len=0;
  while(v[len]) ++len;
  char *p=0;
  const char *prev = Find(k);
  m_values[k]=p=new char[len+1];
  if(!prev) delete [] prev;
  for(int i=0;i<=len;++i) p[i]=v[i];
}

/**
 * @brief Show values of options
 *
 * Write values of options in JSON format.
 * @param stream iostream to which the values is written.
 */
void COption::Show(std::ostream &stream) const{
  stream<<"{ ";
  for(std::map<std::string,const char*>::const_iterator iter= m_values.begin(); iter!=m_values.end(); ++iter){
    if(iter!=m_values.begin()) stream<<",\n  ";
    stream<<"\""<<(*iter).first<<"\": \""<<(*iter).second<<'"';
  }
  stream<<"\n}\n";
}

/**
 * @brief Show help message of options
 *
 * Show help message of options.
 * @param stream iostream to which help messages are written.
 * @param opt    the array of option specifications.
 */
void COption::ShowHelp(std::ostream &stream, const COption::Option_t* opt){
  stream<<" Options:"<<std::endl;
  for(int i=0;opt[i].Name;++i){
    stream<<"  -"<<opt[i].Short;
    if(opt[i].HasArg) stream<<"<value>";
    stream<<"\t--"<<opt[i].Name;
    if(opt[i].HasArg) stream<<"=<value>";
    if(opt[i].Default) stream<<"    [default: "<<opt[i].Default<<"]";
    stream<<std::endl;
    if(opt[i].Usage) stream<<"     "<<opt[i].Usage;
    stream<<std::endl;
  }
}

/*
 * @brief Find the value of an option
 *
 * Find the value of the specified option.
 * @param k long name of an option
 * @return
 * \li pointer to the value, if the argument takes an value and set.
 * \li pointer to "", if the argument does not take any argument and set.
 * \li 0 if the option is not set.
 */
/*
const char *COption::Find(const char *k) const {
  std::map<std::string,const char*>::const_iterator iter = m_values.find(k);
  if(iter==m_values.end()) return 0;
  return (*iter).second;
}
*/

/**
 * @brief Find the value of an option
 *
 * Find the value of the specified option.
 * @param k long name of an option
 * @param value_if_not_set return value used when no value was set
 * @return
 * \li pointer to the value, if the argument takes an value and set.
 * \li pointer to "", if the argument does not take any argument and set.
 * \li 0 if the option is not set.
 */
const char *COption::Find(const char *k, const char *value_if_not_set) const {
  //const char *COption::Find(const char *k) const {
  std::map<std::string,const char*>::const_iterator iter = m_values.find(k);
  if(iter==m_values.end()) return value_if_not_set;
  return (*iter).second;
}

/**
 * @brief Require the value of an option
 *
 * Find the value of the specified option. If not set, cause an exception.
 * @param k long name of an option
 * @return
 * \li pointer to the value, if the argument takes an value and set.
 * \li pointer to "", if the argument does not take any argument and set.
 * \li 0 if the option is not set.
 */
const char *COption::Require(const char *k) const {
  static const char *p0="not set";
  const char *p = Find(k,p0);
  if(p==p0 || !p) Quit("Option '"<<k<<"' required.");
  return p;
}

/**
 * @brief Require the value of an option
 *
 * Find the value of the specified option. If not set, cause an exception.
 * @param k long name of an option
 * @return
 * \li pointer to the value, if the argument takes an value and set.
 * \li pointer to "", if the argument does not take any argument and set.
 * \li 0 if the option is not set.
 */
int64_t COption::RequireInteger(const char *k) const {
  return std::atoll( Require(k) );
}

/**
 * @brief Require the value of an option
 *
 * Find the value of the specified option. If not set, cause an exception.
 * @param k long name of an option
 * @return
 * \li pointer to the value, if the argument takes an value and set.
 * \li pointer to "", if the argument does not take any argument and set.
 * \li 0 if the option is not set.
 */
double COption::RequireDouble(const char *k) const {
  return std::atof( Require(k) );
}

/**
 * @brief Parse command line arguments.
 *
 * Set up the instance of COption and parse command line arguments.
 * @param argc the argc value given to main() function.
 * @param argv the argv value given to main() function.
 * @param opt  the array of option specifications.
 * @return the index of the first command line argument other than options.
 */
int COption::SetUp(int argc, char **argv, COption::Option_t *opt){
  int n=0;
  m_usage=opt;
  // First, count # of options
  while(m_usage[n].Name) ++n;
  // Next, get option array
  struct option *long_options = 0;
  try{
    std::map<int,int> short2index;
    std::string short_options("D:");
    long_options = new struct option[n+1];
    long_options[n].name=0;
    long_options[n].has_arg=0;
    long_options[n].flag=0;
    long_options[n].val=0;
    for(int i=0;i<n;++i){
      // Setup long option
      long_options[i].name=m_usage[i].Name;
      long_options[i].has_arg=m_usage[i].HasArg;
      long_options[i].flag=0;
      long_options[i].val=1000+i;
      // Setup short option
      short2index[ m_usage[i].Short[0] ]=1000+i;
      if(! std::strcmp(m_usage[i].Short,"D") ) Quit("Option 'D' is reserved");
      short_options+=m_usage[i].Short;
      if(m_usage[i].HasArg) short_options+=':';
    }

    for(int32_t i=0;opt[i].Name;++i){
      if(!opt[i].Default) continue;
      //std::cerr<<"COption::SetUp: Adding default value ("<<opt[i].Name<<") -> "<<opt[i].Default<<std::endl;
      Add(opt[i].Name,opt[i].Default);
    }
    for(;;){
      int found = getopt_long(argc,argv,short_options.c_str(),long_options,0);
      if(found<0) break;
      if(found=='?') Quit("Ambigous option: "<<optind);
      // Temporary options
      if(found=='D'){
        std::string k;
        const char *p=optarg;
        while(*p && *p!='=') k+=*p++;
        if(*p=='=') ++p;
        Add(k.c_str(),p);
        continue;
      }
      // Short options
      if(found<1000){
        if(short2index.find(found)!=short2index.end()){
          found=short2index[found];
        }else{
          Quit("Unexpected short option: "<<static_cast<char>(found));
        }
      }
      if(found<1000 || 1000+n <= found ) Quit("COption error: found="<<found<<", index="<<optind);
      found-=1000;
      Add(m_usage[found].Name, m_usage[found].HasArg? optarg: "");
    }
  }catch(...){
    delete [] long_options;
    throw;
  }
  delete [] long_options;
  return optind;
}
