/*
  Cette classe contient les éléments relatifs aux données 

Ces élements sont:
  m_nrows : nombre d'observations
  m_ncols : nombre de variables
*/
#ifndef Data_H
#define Data_H

#include <iostream>
#include <iomanip>
#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace std;
using namespace arma;
using namespace Rcpp;

class Data{
  public:
  int m_n, m_d, m_modalities;
  Mat<double> m_data;
  vector < vector < uvec > > m_whotake;
  Data(){};
  Data(const S4 obj){
   this->m_data = as<mat>(obj.slot("data"));
   this->m_n = as<int>(obj.slot("n"));
   this->m_d = as<int>(obj.slot("d"));
   this->m_modalities = as<int>(obj.slot("modalities"));
   m_whotake.resize(m_d);
   for (int j=0; j<m_d; j++){
     m_whotake[j].resize(m_modalities);
     for (int h=0; h<m_modalities; h++){
       m_whotake[j][h] = find(m_data.col(j)==h);
     }
   }
  }

  ~Data(){};
  
};
#endif
