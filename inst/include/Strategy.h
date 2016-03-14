#ifndef Strategy_H
#define Strategy_H

#include <iostream>
#include <iomanip>
#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace std;
using namespace arma;
using namespace Rcpp;

class Strategy{
  public:
  int m_iterMH, m_nbSmall, m_iterSmall, m_nbKeep, m_iterKeep;
  double m_tolKeep;
  
  Strategy(){}
  Strategy(const S4 obj){
   this->m_iterMH = as<int>(obj.slot("iterMH"));
   this->m_nbSmall = as<int>(obj.slot("nbSmall"));
   this->m_iterSmall = as<int>(obj.slot("iterSmall"));
   this->m_nbKeep = as<int>(obj.slot("nbKeep"));
   this->m_iterKeep = as<int>(obj.slot("iterKeep"));
   this->m_tolKeep = as<double>(obj.slot("tolKeep"));
  };
  ~Strategy(){};

};
#endif
