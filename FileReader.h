/**
 * @file    FileReader.h
 * @brief   Read tab, comma separated test files
 *
 * @author  Tomohiro Yasuda
 * @date    2010-11-20
 *
 */

#ifndef _FILE_READER_H_
#define _FILE_READER_H_

#include <iostream>
#include <string>
#include <vector>
#include <stdint.h>
#include <map>

class CProgressReport {
 private:
  uint64_t m_previous_turn;
  int64_t m_progress_interval;
  std::string m_unit;
 public:
  CProgressReport();
  //inline void SetInterval(int64_t interval, const char *unit){m_progress_interval=interval; m_unit=unit;}
  void SetInterval(int64_t interval, const char *unit);
  void ShowProgress(int64_t progress);
};

class CFileReader {
private:
  std::string m_filename;
  FILE *m_fp;
  FILE *m_my_fp;
  FILE *m_pipe_fp;
  char *m_buffer;
  bool m_erase_newline;
  //static const int BUFFER_SIZE=65536;
  static const int BUFFER_SIZE=1024*1024*16;
  uint64_t m_line_number;
  uint64_t m_read_amount;
  //void Initialize(){m_fp=m_my_fp=m_pipe_fp=0; m_buffer=0; m_read_amount=0;}
  void Initialize();
  CProgressReport m_progress_reporter;
  FILE *OpenPipe(const char *filename, const char *suffix, const char *cmd);
public:
  CFileReader(){Initialize();}
  CFileReader(const char *fn){Initialize(); Open(fn);}
  //void SetProgressInterval(int64_t interval){m_progress_interval=interval;}
  void SetProgressInterval(int64_t t, const char *p){m_progress_reporter.SetInterval(t,p);}
  ~CFileReader(){Close();}
  void Open(const char *filename);
  const char *CurrentLine() const {return m_buffer;}
  const char *GetLine();
  const char *GetContentLine(const char *comment_sym=0);
  bool GetContentLine(std::vector<std::string> *str_vec);
  void Close();
  void Chomp();
  uint64_t Line() const {return m_line_number;}
  uint64_t ReadAmount() const {return m_read_amount;}
};

class CTokenizer {
private:
  std::string m_token;
  const char *m_delimiter;
  const char *m_source_ptr;
  const char *m_prev;
  int32_t m_field_id;
public:
  CTokenizer(const char *p, const char *d){m_source_ptr=p; m_delimiter=d; m_field_id=0; m_prev="";}
  int32_t FieldID() const {return m_field_id;}
  inline bool hasNext() const {return m_source_ptr && *m_source_ptr;}
  const char *Next(const char *d=0);
  const char *Prev(){return m_prev;}
  const char *CStr() const {return m_token.c_str();}
  const char *Peep() const {return m_source_ptr;}
  //int64_t NextInteger(const char *d=0);
  inline int64_t NextInteger(const char *d=0){Next(d); return Integer();}
  inline std::string NextString(const char *d=0){Next(d); return m_token;}
  int64_t Integer() const;
};

class CKVStore {
 private:
  typedef std::map<std::string,std::string> mapping_t;
  mapping_t m_map;
  inline void Initialize(){}
 public:
  CKVStore(){Initialize();}
  bool Parse(const char **p);
  bool Read(const char *target_file);
  void Write(std::ostream &stream, int32_t indent=0) const;
  void Show (std::ostream &stream, int32_t indent=0) const {Write(stream,indent); stream<<std::endl;}
  const char *Find(const char *key) const;
  bool Add(const char *key, const char *value);
  bool Empty() const {return m_map.empty();}
};


#endif // _FILE_READER_H_
