#include "MutationJointOverdisp.hpp"


MutationJointOverdisp::MutationJointOverdisp(std::vector<MutationInd*>* mutIndVec, bool verbose, bool gradientDebug) : MutationListContainer(nullptr, 0, 1, verbose, gradientDebug) { // 0 mutCounts, 1 overdispVar
  this->mutIndVec = mutIndVec;
}

MutationJointOverdisp::~MutationJointOverdisp() {
  for(unsigned int mutIndIdx = 0; mutIndIdx < this->mutIndVec->size(); mutIndIdx++) {
    delete (*this->mutIndVec)[mutIndIdx];
  }
}

/*
 * function to calculate the joint likelihood of all the cells, given the current omega (overdispersion) parameter for the beta binomial
 */
double MutationJointOverdisp::getLogLikelihood() {
  double currOmega = this->getMutationVarEst(0);
  double totalLl = 0;
  double currLl = 0;

  // for each MutationInd object
  for(unsigned int mutIndIdx = 0; mutIndIdx < this->mutIndVec->size(); mutIndIdx++) {
    // set current omega value in internal MutationLists
    (*this->mutIndVec)[mutIndIdx]->setAllMutOverdispOmega(currOmega);

    // calc likelihood, using whatever value of X=mutation count is already stored
    currLl = (*this->mutIndVec)[mutIndIdx]->getLogLikelihood();
    totalLl += currLl;
  }

  return totalLl;
}

// no upper limit on param value omega
double MutationJointOverdisp::getValidOptimParamMax() const {
  return std::numeric_limits<double>::max();
}

