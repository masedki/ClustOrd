#include "Metropolis.h"

//[[Rcpp::export]]
S4  OptimizeBIC(S4 reference){
  S4 * reference_p=&reference;  
  Metropolis * metro_p = new Metropolis(reference_p);
  metro_p->Run();
  
  S4 model = reference.slot("model");
  model.slot("g") = metro_p->p_best->p_model->m_g;
  model.slot("omega") = trans(metro_p->p_best->p_model->m_omega);

  S4 criteria = reference.slot("criteria");
  criteria.slot("BIC") = metro_p->p_best->p_criteria->m_bic;
  criteria.slot("ICL") = metro_p->p_best->p_criteria->m_icl;
  criteria.slot("loglikelihood") = metro_p->p_best->p_criteria->m_loglike;
  criteria.slot("nbparam") = metro_p->p_best->p_criteria->m_nbparam;
  
  S4 param = reference.slot("param");
  param.slot("pi") = trans(metro_p->p_best->p_param->m_pi);
  param.slot("alpha")= metro_p->p_best->p_param->m_alpha ;
  param.slot("beta")= metro_p->p_best->p_param->m_beta ;
  param.slot("epsilon")= metro_p->p_best->p_param->m_epsilon ;
  
  S4 detailsMH = reference.slot("detailsMH");
  detailsMH.slot("allbic") = metro_p->allbic;
  
  delete metro_p, reference_p;
  
  return reference;
}

