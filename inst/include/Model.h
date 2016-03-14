#ifndef Model_H
#define Model_H

#include <iostream>
#include <iomanip>
#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace std;
using namespace arma;
using namespace Rcpp;

class Model{
  public:
  int m_g;
  Col<double> m_omega;
  Col<int> m_sizeblock;
  
  Model(){};
  Model(const S4  obj){
   this->m_g = as<int>(obj.slot("g"));
   this->m_omega = as<vec>(obj.slot("omega"));
   m_sizeblock.resize(m_omega.n_rows);
   for (int b=0; b<m_sizeblock.n_rows;b++){
     m_sizeblock(b)= sum(m_omega==b);
   }
  };
  ~Model(){};
  
  void NeighborModel(const Model * current_p){
    m_g = current_p->m_g;
    m_omega = current_p->m_omega;
    ivec who = randi<ivec>(1, distr_param(0, m_omega.n_rows-1));
    int bsup = max(m_omega);
    if (bsup<m_omega.n_rows)
      bsup++;
    ivec affect = randi<ivec>(1, distr_param(0, bsup));
    Col<double> tmpomega = m_omega;
    tmpomega(who(0)) = affect(0);
    int locmax = 0;
    uvec place=find(m_omega==0);
    m_omega(0) = 0;
    for (int u=1; u<m_omega.n_rows; u++){
      if (sum(tmpomega.head(u)==tmpomega(u) )>0){
        place=find(tmpomega.head(u)==tmpomega(u) );
        m_omega(u) = m_omega(place(0));
      }else{
        locmax++;
        m_omega(u)=locmax;
      }
    }
   // cout << "Converti " << trans(m_omega) << endl;
    
    
    m_sizeblock.resize(m_omega.n_rows);
    for (int b=0; b<m_sizeblock.n_rows;b++){
     m_sizeblock(b)= sum(m_omega==b);
    }
    
  }
    
};
#endif
