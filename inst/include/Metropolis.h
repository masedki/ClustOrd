#include "XEM.h"

#ifndef Metropolis_H
#define Metropolis_H

class Metropolis{
  public:
  Results * p_best;
  Results * p_current;
  Results * p_candidate;
  XEM * p_xem;
//  Mat<double> Bestmodel, Currentmodel, Candidatemodel;
  Mat<double> allbic;
  
  Metropolis(){};
  Metropolis( const S4 *);
  ~Metropolis(){};
  void Run();
};
#endif
