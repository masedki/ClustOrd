/*
  Cette classe contient les éléments relatifs aux données 

Ces élements sont:
  m_nrows : nombre d'observations
  m_ncols : nombre de variables
*/
#ifndef Criteria_H
#define Criteria_H

#include <iostream>
#include <iomanip>
#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace std;
using namespace arma;
using namespace Rcpp;

class Criteria{
  public:
  double m_bic, m_loglike, m_icl, m_nbparam;

  Criteria(){m_bic= log(0); m_loglike= log(0); m_nbparam= 0; m_icl=log(0);};
  Criteria(const S4  obj){m_bic= as<double>(obj.slot("BIC")); m_icl= as<double>(obj.slot("ICL")); m_loglike= as<double>(obj.slot("loglikelihood")); m_nbparam= as<double>(obj.slot("nbparam"));};
  ~Criteria(){};
  
};
#endif
