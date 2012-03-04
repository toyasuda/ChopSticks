/**
 * @file    Option.cc
 * @brief   Parse options with getopt(3) and keep the values in map<string,const char*>
 *
 * @author  Tomohiro Yasuda
 * @date    2010-9-20
 *
 */

// (C) Yasuda, Tomohiro: The university of Tokyo

#ifndef _OPTION_H_
#define _OPTION_H_

#include <iostream>
#include <map>
#include <string>
#include <stdint.h>

/**
 * @brief Parse options with getopt(3) and keep the values in map<string,const char*>
 *
 * This class takes argc, argv, and the array of specifications of options
 * by SetUp() method and keep the (key,value) pairs.
 * The values can be obtained by Find() method.
 * If the value for options without arguments, its value will be "".
 * If values of options not specified by command line arguments is 0.
 */
class COption {
public:
  struct Option_t {
    const char *Name;  ///< Long name of an option
    const char *Short; ///< Short name of an option
    int HasArg;        ///< 1 if this option takes a value, 0 otherwise
    const char *Usage; ///< Help message of this option
    const char *Default; ///< Help message of this option
  };
private:
  Option_t *m_usage; ///< Pointer to the array of option specifications
  std::map<std::string,const char *> m_values; ///< map instance that keep values of options
  /**
   * @brief Global variable that can keep an instance of COption class.
   */
  static COption *s_option;
  friend COption& Option();
public:
  //const char *Find(const char *k) const; ///< Get the value of an option
  const char *Find(const char *k, const char *value_if_not_set=0) const; ///< Get the value of an option
  const char *Require(const char *k) const; ///< Get the value of an required option
  int64_t RequireInteger(const char *k) const; ///< Get the value of an required option that is an integer
  double RequireDouble(const char *k) const; ///< Get the value of an required option that is an integer
  int SetUp(int argc, char **argv, COption::Option_t *opt); ///< Parse command line arguments
  void Add(const char *k, const char *v);
  void Show(std::ostream &stream) const; ///< Show all values of options
  static void ShowHelp(std::ostream &stream, const COption::Option_t* opt); ///< Show help messages for options
  COption(){m_usage=0;} ///< Constructor
  ~COption();           ///< Destructor
};

/**
 * @brief Instance of COption for this process
 *
 * This global function returns a reference to an instance of COption.
 * If the instance has not created yet, this function creats it.
 */
COption& Option();

#endif // _OPTION_H_
