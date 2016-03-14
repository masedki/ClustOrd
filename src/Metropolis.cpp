#include "Metropolis.h"

Metropolis::Metropolis(const S4 * obj_p){
  p_best = new Results(obj_p);
  p_current = new Results(obj_p);
  p_candidate = new Results(obj_p);
  p_xem = new XEM(p_candidate);
  allbic = zeros<mat>(p_best->p_strategy->m_iterMH, 3);
};

void Metropolis::Run(){
  Col<double> proba_accept = randu( p_best->p_strategy->m_iterMH );
  double probaproposal = 0, proposalcand=0, proposalcurrent=0;
  proposalcurrent = max(p_current->p_model->m_omega);
  if (proposalcurrent < p_current->p_model->m_omega.n_rows)
  proposalcurrent ++;
  proposalcurrent = 1/proposalcurrent;
  
  for (int iter=0; iter<p_best->p_strategy->m_iterMH; iter++){
    p_candidate->SampleNeighbor(p_current);
    proposalcand = max(p_candidate->p_model->m_omega);
    if (proposalcand < p_candidate->p_model->m_omega.n_rows) proposalcand ++;
    proposalcand = 1/proposalcand;
    p_xem->Redefinition(p_candidate);
    p_xem->Run();
    p_candidate = p_xem->p_results;
    
    // sauvegarde
    /*allbic(iter,0) = p_best->p_criteria->m_bic ;
    allbic(iter,1) = p_current->p_criteria->m_bic ;
    allbic(iter,2) = p_candidate->p_criteria->m_bic;*/
    allbic(iter,0) = p_best->p_criteria->m_icl ;
    allbic(iter,1) = p_current->p_criteria->m_icl ;
    allbic(iter,2) = p_candidate->p_criteria->m_icl;
    
    
    
    probaproposal = proposalcurrent/proposalcand;
//    if (exp(p_candidate->p_criteria->m_bic - p_current->p_criteria->m_bic)*probaproposal > proba_accept(iter)){
   if (exp(p_candidate->p_criteria->m_icl - p_current->p_criteria->m_icl)*probaproposal > proba_accept(iter)){
      p_current->CopySparse(p_candidate);
      proposalcurrent = proposalcand;
 //     if (p_current->p_criteria->m_bic > p_best->p_criteria->m_bic)
      if (p_current->p_criteria->m_icl > p_best->p_criteria->m_icl)
      p_best->CopySparse(p_current);
    }
     
  }
};
