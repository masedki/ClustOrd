
#ifndef Param_H
#define Param_H
  
#include "Model.h"
  
class Param {
  public:
    Col<double> m_pi;
    vector< vector< Mat<double> > > m_alpha;
    vector< Mat<double> > m_beta;
    Mat<double> m_epsilon;
  
  Param(){};
  Param(const Model *, const int );
  ~Param(){};
};
#endif