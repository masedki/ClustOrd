#include "Results.h"

#ifndef XEM_H
#define XEM_H

class XEM{
  public:
  Results * p_results;
  Col<double> loglikeSmall, maxtmplogproba, rowsums;
  vector<Param> paramCand;
  Param * paramCurrent_p;
  Col<double>tmpval, tmpdep, tmpindep;
  vector <vector< Mat<double> > > m_probablock;
  Mat<double> m_probacompo;
  int iterCurrent;
  double loglikeoutput, bicoutput;
  vector< uvec > idxvbles;
  Mat<double> onDiag;
  
  XEM(){};
  XEM(Results * );
  ~XEM(){};
  
  void Redefinition(Results *);
  double ComputeLogLike();
  void Estep();
  void Mstep();
  void OneEM();
  void SwitchParamCurrent(int);
  void Run();  
  void ProbaComputation();
};
#endif
