#include "Data.h"
#include "Model.h"
#include "Strategy.h"
#include "Criteria.h"
#include "Param.h"

#ifndef Results_H
#define Results_H

class Results{
  public:
  const Data* p_data;
  Model* p_model;
  const Strategy* p_strategy;
  Criteria* p_criteria;
  Param* p_param;

  Results(){};
  Results( const S4 *);
  ~Results(){};
  
  void UpdateCriteria();
  void SampleNeighbor(const Results *);
  void CopySparse(const Results *);
};
#endif
