#include "Param.h"

Param::Param(const Model * p_model, const int modalities){
  this->m_alpha.resize(p_model->m_omega.n_rows);
  this->m_beta.resize(p_model->m_omega.n_rows);
  this->m_pi = ones<vec>(p_model->m_g)/p_model->m_g;  
  m_epsilon = zeros<mat>(p_model->m_g, p_model->m_omega.n_rows);
  for (int b=0; b<=max(p_model->m_omega); b++){
    if (p_model->m_sizeblock(b) == 0){
      m_alpha[b].resize(0);
      m_beta[b] = randu<mat>(0,0);
    }else if (p_model->m_sizeblock(b) == 1){
      m_alpha[b].resize(1);
      m_alpha[b][0] = randu(p_model->m_g, modalities);
      for (int k=0; k<p_model->m_g; k++) m_alpha[b][0].row(k) = m_alpha[b][0].row(k) / sum(m_alpha[b][0].row(k));
      m_beta[b] = randu<mat>(0,0);
    }else{
      m_alpha[b].resize(p_model->m_sizeblock(b));
      for (int j=0; j<p_model->m_sizeblock(b); j++){
        m_alpha[b][j] = randu(p_model->m_g, modalities);
        for (int k=0; k<p_model->m_g; k++) m_alpha[b][j].row(k) = m_alpha[b][j].row(k) / sum(m_alpha[b][j].row(k));
      }      
      m_beta[b] = randu(p_model->m_g, modalities);
      for (int k=0; k<p_model->m_g; k++) m_beta[b].row(k) = m_beta[b].row(k) / sum(m_beta[b].row(k));
      m_epsilon.col(b) = randu(p_model->m_g);
    }    
  }
};
