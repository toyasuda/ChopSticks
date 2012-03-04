////////////////////////////////////////////////////////////////////////////////
// Various useful classes
// (C) Yasuda, Tomohiro: The university of Tokyo

#include "Tool.h"
#include "Utility.h"

////////////////////////////////////////////////////////////////////////////////
// Slot manager
void CSlotManager::SetUp(int32_t sz){
  if(sz<1) Quit("Invalid size of slots: "<<sz);
  delete [] m_slots;
  m_slots=0;
  m_n_used=0;
  m_slots = new int32_t[m_n_slots=sz];
  m_free_top=1;
  for(int32_t i=0;i<m_n_slots-1;++i) m_slots[i]=i+1;
  m_slots[m_n_slots-1]=INVALID;
  m_slots[0]=INVALID;
}

int32_t CSlotManager::Allocate(){
  if(m_free_top<=0){
    Assert(1+m_n_used==m_n_slots);
    Assert(m_free_top==INVALID);
    return INVALID;
  }
  int32_t new_one=m_free_top;
  m_free_top = m_slots[m_free_top];
  m_slots[new_one]=USED;
  ++m_n_used;
  return new_one;
}

void CSlotManager::Release(int32_t slot){
  if(slot<=0 || m_n_slots<=slot) Quit("Invalid slot: "<<slot);
  if(m_slots[slot]!=USED) Quit("Slot "<<slot<<" is already free");
  m_slots[slot]=m_free_top;
  --m_n_used;
  m_free_top=slot;
}

void CSlotManager::Write(std::ostream& stream,int32_t indent) const {
  stream<<"{\"freetop\":"<<m_free_top<<",\"used\":"<<m_n_used<<",\"slots\":[";
  for(int32_t i=0;i<m_n_slots;++i){
    if(i) stream<<',';
    stream<<m_slots[i];
  }
  stream<<']';
}

void CSlotManager::Test(std::ostream& stream){
  stream<<"CSlotManager::Test: --------------------------------------------------"<<std::endl;
  CSlotManager sm;
  sm.SetUp(5);
  for(int32_t i=0;i<5+1;++i){
    int32_t a=sm.Allocate();
    stream<<"Allocate: i="<<i<<": "<<a<<": "<<sm<<std::endl;
  }
  for(int32_t i=0;i<5+1;++i){
    try{
      sm.Release(i);
      stream<<"Release: i="<<i<<": "<<sm<<std::endl;
    }catch(...){
      stream<<"Release: i="<<i<<": Exception occured."<<std::endl;
    }
  }
  for(int32_t i=0;i<5+1;++i){
    int32_t a=sm.Allocate();
    stream<<"Allocate: i="<<i<<": "<<a<<": "<<sm<<std::endl;
  }
  for(int32_t i=0;i<5+1;++i){
    try{
      sm.Release(5-i);
      stream<<"Release: i="<<i<<": "<<sm<<std::endl;
    }catch(...){
      stream<<"Release: i="<<i<<": Exception occured."<<std::endl;
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
// Mon Jan  2 11:30:45 2012
int64_t interval_distance(int32_t start1, int32_t end1, int32_t start2, int32_t end2){
  Assert(start1<=end1);
  Assert(start2<=end2);
  if(start1<=start2 && start2<end1) return 0; // overlap
  if(start1<end2    && end2 <=end1) return 0; // overlap
  if(end1<=start2) return start2-end1;
  if(end2<=start1) return start1-end2;
  Quit("Unexpected intervals: ["<<start1<<","<<end1<<") vs. ["<<start2<<","<<end2<<")");
}

