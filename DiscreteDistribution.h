/**
 * @file    DiscreteDistribution.h
 * @brief   Calculate binomial or poisson distribution
 *
 * @author  Tomohiro Yasuda
 * @date    2010-10-2
 *
 */

#ifndef  _DISCRETE_DISTRIBUTION_H_
#define  _DISCRETE_DISTRIBUTION_H_

#include <iostream>
#include <stdint.h>
//#include <cassert>

class CDiscreteDistribution {
private:
  double *m_P;
  double *m_logP;
  double *m_F;
  double m_p;
  int64_t m_N;
  int64_t m_max_k;
  enum {NONE, BINOMIAL, POISSON} m_type;
  void SetUp(int64_t maxk, double logP0, double Pnext);
public:
  CDiscreteDistribution(){m_P=m_logP=m_F=0; m_p=m_N=m_max_k=0;}
  ~CDiscreteDistribution();
  void SetUpBinomial(int64_t n, int64_t maxk, double p);
  void SetUpPoisson (int64_t n, int64_t maxk, double p){SetUpPoisson(maxk, n*p);}
  void SetUpPoisson (int64_t maxk, double lambda);
  inline double    P(int k) const {return (0<=k && k<=m_max_k)? m_P[k]:    0;}
  inline double logP(int k) const {return (0<=k && k<=m_max_k)? m_logP[k]: 0;}
  inline double    F(int k) const {return (0<=k && k<=m_max_k)? m_F[k]:    0;}
  double gslP(int k) const;
  inline double operator()(int k) const {return P(k);}
  void Show(std::ostream &stream) const;
};


#endif // _DISCRETE_DISTRIBUTION_H_
