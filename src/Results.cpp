#include "Results.h"

Results::Results(const S4 * obj_p){
  p_data = new Data(as<S4>(obj_p->slot("data")));
  p_model = new Model(as<S4>(obj_p->slot("model")));
  p_strategy = new Strategy(as<S4>(obj_p->slot("strategy")));
  p_criteria = new Criteria(as<S4>(obj_p->slot("criteria")));
  p_param = new Param(p_model, p_data->m_modalities);
};
 


void Results::SampleNeighbor(const Results * current_p){
  p_model->NeighborModel(current_p->p_model);
}

void Results::CopySparse(const Results * current_p){
  p_model->m_omega = current_p->p_model->m_omega;
  p_criteria->m_bic = current_p->p_criteria->m_bic;
  p_criteria->m_icl = current_p->p_criteria->m_icl;
  p_criteria->m_nbparam = current_p->p_criteria->m_nbparam;
  p_criteria->m_loglike = current_p->p_criteria->m_loglike;
  p_param->m_pi = current_p->p_param->m_pi;
  p_param->m_alpha = current_p->p_param->m_alpha;
  p_param->m_beta = current_p->p_param->m_beta;
  p_param->m_epsilon = current_p->p_param->m_epsilon;
}
