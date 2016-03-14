#include "XEM.h"

XEM::XEM(Results * p_input){
  p_results = p_input;
  loglikeSmall = ones<vec>(p_results->p_strategy->m_nbSmall) * log(0);
  for (int i=0;i<p_results->p_strategy->m_nbSmall;i++) paramCand.push_back( Param(p_results->p_model, p_results->p_data->m_modalities ) );
  m_probablock.resize(p_results->p_model->m_g);
  for (int k=0; k<p_results->p_model->m_g; k++){
    m_probablock[k].resize(max(p_results->p_model->m_omega)+1);
  }
  m_probacompo=ones<mat>(p_results->p_data->m_n, p_results->p_model->m_g);
  onDiag = ones<mat>(p_results->p_data->m_n, max(p_results->p_model->m_omega)+1);
  idxvbles.resize(p_results->p_data->m_d);
  for (int b=0; b<=max(p_results->p_model->m_omega); b++){
    if (p_results->p_model->m_sizeblock(b)>1){
      idxvbles[b] = find(p_results->p_model->m_omega == b);
      for (int j=0; j<idxvbles[b].n_rows; j++){
        onDiag.col(b) = onDiag.col(b) % (p_results->p_data->m_data.col(idxvbles[b](0)) == p_results->p_data->m_data.col(idxvbles[b](j)));
      }
      for (int k=0; k<p_results->p_model->m_g; k++) m_probablock[k][b] = ones<mat>(p_results->p_data->m_n, 2);
    }else if (p_results->p_model->m_sizeblock(b)==1){
      idxvbles[b] = find(p_results->p_model->m_omega == b);
      for (int k=0; k<p_results->p_model->m_g; k++) m_probablock[k][b] = ones<mat>(0,0);
    }else{
      idxvbles[b].resize(0);
      for (int k=0; k<p_results->p_model->m_g; k++) m_probablock[k][b] = ones<mat>(0,0);
    }
  }
  iterCurrent = p_results->p_strategy->m_iterSmall;
  loglikeoutput = log(0);
}

void XEM::Redefinition(Results * p_input){
  p_results = p_input;
  loglikeSmall = ones<vec>(p_results->p_strategy->m_nbSmall) * log(0);
  paramCand.resize(p_results->p_strategy->m_nbSmall);
  for (int i=0;i<p_results->p_strategy->m_nbSmall;i++) paramCand[i] = Param(p_results->p_model, p_results->p_data->m_modalities );
  m_probablock.resize(p_results->p_model->m_g);
  for (int k=0; k<p_results->p_model->m_g; k++){
    m_probablock[k].resize(max(p_results->p_model->m_omega)+1);
  }
  //m_probacompo=ones<mat>(p_results->p_data->m_n, p_results->p_model->m_g);
  onDiag = ones<mat>(p_results->p_data->m_n, max(p_results->p_model->m_omega)+1);
  idxvbles.resize(p_results->p_data->m_d);
  for (int b=0; b<=max(p_results->p_model->m_omega); b++){
    if (p_results->p_model->m_sizeblock(b)>1){
      idxvbles[b] = find(p_results->p_model->m_omega == b);
      for (int j=0; j<idxvbles[b].n_rows; j++){
        onDiag.col(b) = onDiag.col(b) % (p_results->p_data->m_data.col(idxvbles[b](0)) == p_results->p_data->m_data.col(idxvbles[b](j)));
      }
      for (int k=0; k<p_results->p_model->m_g; k++)      m_probablock[k][b] = ones<mat>(p_results->p_data->m_n, 2);
    }else if (p_results->p_model->m_sizeblock(b)==1){
      idxvbles[b] = find(p_results->p_model->m_omega == b);
      for (int k=0; k<p_results->p_model->m_g; k++) m_probablock[k][b] = ones<mat>(0,0);
    }else{
      idxvbles[b].resize(0);
      for (int k=0; k<p_results->p_model->m_g; k++) m_probablock[k][b] = ones<mat>(0,0);
    }
  }
  iterCurrent = p_results->p_strategy->m_iterSmall;
  loglikeoutput = log(0);
}

void XEM::Run(){
  // Partie Small EM
  for (int ini=0; ini<p_results->p_strategy->m_nbSmall; ini++){
    SwitchParamCurrent(ini);
    OneEM();
    loglikeSmall(ini) = ComputeLogLike();
  }
  uvec indices = sort_index(loglikeSmall);
  iterCurrent = p_results->p_strategy->m_iterKeep;
  if (p_results->p_strategy->m_nbSmall > p_results->p_strategy->m_nbKeep)   
  loglikeSmall( indices.head(p_results->p_strategy->m_nbSmall - p_results->p_strategy->m_nbKeep) ) = loglikeSmall( indices.head(p_results->p_strategy->m_nbSmall - p_results->p_strategy->m_nbKeep) ) + log(0);
  
  
  for (int tmp1=0; tmp1<p_results->p_strategy->m_nbKeep; tmp1++){
    SwitchParamCurrent(indices(p_results->p_strategy->m_nbSmall - tmp1 - 1));
    OneEM();
    loglikeSmall(indices(p_results->p_strategy->m_nbSmall - tmp1 - 1)) = ComputeLogLike();
  }
  uword  index;
  double indicebest = (loglikeSmall).max(index);
  SwitchParamCurrent(index);
  p_results->p_criteria->m_loglike = ComputeLogLike();
  p_results->p_criteria->m_nbparam = p_results->p_model->m_g * p_results->p_data->m_d*(p_results->p_data->m_modalities-1)   + (p_results->p_model->m_g -1) + p_results->p_model->m_g * (p_results->p_data->m_modalities-1) *sum(p_results->p_model->m_sizeblock>1);
  p_results->p_criteria->m_bic = p_results->p_criteria->m_loglike - 0.5*p_results->p_criteria->m_nbparam*log(p_results->p_data->m_n);
    for (int k=0; k<m_probacompo.n_cols; k++) m_probacompo.col(k) = m_probacompo.col(k)/rowsums;

  colvec ent = max(m_probacompo,1);
  double entropie = sum(log(ent));
  p_results->p_criteria->m_icl = p_results->p_criteria->m_bic + entropie;
  p_results->p_param = paramCurrent_p;
}


void XEM::ProbaComputation(){
  m_probacompo=ones<mat>(p_results->p_data->m_n, p_results->p_model->m_g);
  for (int k=0; k<p_results->p_model->m_g; k++){
    tmpval = ones<vec>(p_results->p_data->m_n) * paramCurrent_p->m_pi(k);
    for (int b=0; b<=max(p_results->p_model->m_omega);b++){
      if (p_results->p_model->m_sizeblock(b) == 1){
        for (int h=0; h<p_results->p_data->m_modalities; h++){
          tmpval(p_results->p_data->m_whotake[idxvbles[b](0)][h]) = tmpval(p_results->p_data->m_whotake[idxvbles[b](0)][h])  * paramCurrent_p->m_alpha[b][0](k, h);
        }          
      }else if (p_results->p_model->m_sizeblock(b) > 1){
        tmpindep = (1-paramCurrent_p->m_epsilon(k,b)) * ones<vec>(p_results->p_data->m_n);
        for (int j=0; j<p_results->p_model->m_sizeblock(b);j++){
          for (int h=0; h<p_results->p_data->m_modalities; h++){
            tmpindep(p_results->p_data->m_whotake[idxvbles[b](j)][h]) = tmpindep(p_results->p_data->m_whotake[idxvbles[b](j)][h])  * paramCurrent_p->m_alpha[b][j](k, h);
          }   
        }
        m_probablock[k][b].col(0) = tmpindep;
        tmpdep = (paramCurrent_p->m_epsilon(k,b)) * onDiag.col(b);
        for (int h=0; h<p_results->p_data->m_modalities; h++){
          tmpdep(p_results->p_data->m_whotake[idxvbles[b](0)][h]) = tmpdep(p_results->p_data->m_whotake[idxvbles[b](0)][h]) * paramCurrent_p->m_beta[b](k,h);
        }
        m_probablock[k][b].col(1) = tmpdep;
        tmpval = tmpval % (tmpdep + tmpindep);
      }   
    }
    m_probacompo.col(k)=tmpval;
  }
  m_probacompo = log(m_probacompo);
}


double XEM::ComputeLogLike(){
  ProbaComputation();
  maxtmplogproba = max(m_probacompo, 1);
  double output=0;
  if (min(maxtmplogproba) == 0){
    output = -999999999999;
  }else{
    for (int k=0; k<m_probacompo.n_cols; k++) m_probacompo.col(k)-=maxtmplogproba;
    m_probacompo = exp(m_probacompo);
    rowsums = sum(m_probacompo,1);
    output = sum(maxtmplogproba) + sum(log(rowsums));
  }
  return  output;
}

void XEM::SwitchParamCurrent(int ini){paramCurrent_p = &paramCand[ini];}

void XEM::Estep(){
  for (int k=0; k<m_probacompo.n_cols; k++) {
    m_probacompo.col(k) = m_probacompo.col(k)/rowsums;
    for (int b=0; b<=max(p_results->p_model->m_omega);b++){
      if (p_results->p_model->m_sizeblock(b) > 1){
        m_probablock[k][b].col(0) = m_probablock[k][b].col(0) / (m_probablock[k][b].col(0) + m_probablock[k][b].col(1));
        m_probablock[k][b].col(1) = ones<vec>(m_probablock[k][b].n_rows) - m_probablock[k][b].col(0);
      }
    }
  }
}


void XEM::Mstep(){
  for (int k=0; k<m_probacompo.n_cols; k++) {
    paramCurrent_p->m_pi(k) = sum(m_probacompo.col(k))/sum(sum(m_probacompo));
    rowsums=m_probacompo.col(k);    
    for (int b=0; b<=max(p_results->p_model->m_omega);b++){
      if (p_results->p_model->m_sizeblock(b) == 1){
        for (int h=0; h<p_results->p_data->m_modalities; h++){
          paramCurrent_p->m_alpha[b][0](k, h) = sum(rowsums(p_results->p_data->m_whotake[idxvbles[b](0)][h]));
        }
        paramCurrent_p->m_alpha[b][0].row(k) = paramCurrent_p->m_alpha[b][0].row(k)/sum(paramCurrent_p->m_alpha[b][0].row(k));
      }else if (p_results->p_model->m_sizeblock(b) > 1){
        if (any(rowsums==0)){
          uvec exi=find(rowsums==0);
          for (int loci=0; loci<exi.n_rows;loci++)
          m_probablock[k][b].row(exi(loci)) = trans(zeros<vec>(2));
        }
        tmpindep = m_probablock[k][b].col(0)%rowsums;
        for (int j=0; j<p_results->p_model->m_sizeblock(b);j++){
          for (int h=0; h<p_results->p_data->m_modalities; h++){
            paramCurrent_p->m_alpha[b][j](k, h) = sum(tmpindep(p_results->p_data->m_whotake[idxvbles[b](j)][h]));
          }
          paramCurrent_p->m_alpha[b][j].row(k) = paramCurrent_p->m_alpha[b][j].row(k)/sum(paramCurrent_p->m_alpha[b][j].row(k));
        }
        tmpdep = m_probablock[k][b].col(1)%rowsums;
        for (int h=0; h<p_results->p_data->m_modalities; h++){
          paramCurrent_p->m_beta[b](k,h) = sum(tmpdep(p_results->p_data->m_whotake[idxvbles[b](0)][h]));
        }
        paramCurrent_p->m_beta[b].row(k) = paramCurrent_p->m_beta[b].row(k)/sum(paramCurrent_p->m_beta[b].row(k));
        paramCurrent_p->m_epsilon(k,b) = sum(tmpdep) / (sum(tmpdep) + sum(tmpindep));
      }   
    }
  }
}

void XEM::OneEM(){
  double loglike = ComputeLogLike(), prec = log(0);
  int it=0;
  while ( (it<iterCurrent) && ((loglike-prec)>p_results->p_strategy->m_tolKeep) ){
    it ++;
    Estep();
    Mstep();
    prec = loglike;
    loglike = ComputeLogLike();
  }
  // Une verif
  if (prec>(loglike+p_results->p_strategy->m_tolKeep)) cout << "pb EM " << "prec" << prec << " loglike " << loglike << endl;;
}
