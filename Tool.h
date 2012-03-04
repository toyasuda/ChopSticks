////////////////////////////////////////////////////////////////////////////////
// Various useful classes
// (C) Yasuda, Tomohiro: The university of Tokyo

#ifndef _TOOL_H_
#define _TOOL_H_

#include <stdint.h>
#include <iostream>
#include <string>

class CSlotManager {
 private:
  int32_t* m_slots;
  int32_t m_free_top;
  int32_t m_n_slots;
  int32_t m_n_used;
  void Initialize(){m_slots=0; m_free_top=-1; m_n_slots=m_n_used=0;}
 public:
  static const int32_t INVALID=0;
  static const int32_t USED=-2;
  inline CSlotManager(){Initialize();}
  inline CSlotManager(int32_t sz){Initialize(); SetUp(sz);}
  inline int32_t Size() const {return m_n_used;}
  void SetUp(int32_t sz);
  int32_t Allocate();
  void Release(int32_t slot);
  void Write(std::ostream &stream, int32_t indent=0) const;
  static void Test(std::ostream &stream);
};
inline std::ostream& operator<<(std::ostream& s, const CSlotManager& m){m.Write(s); return s;}

int64_t interval_distance(int32_t start1, int32_t end1, int32_t start2, int32_t end2);

#endif // _TOOL_H_
