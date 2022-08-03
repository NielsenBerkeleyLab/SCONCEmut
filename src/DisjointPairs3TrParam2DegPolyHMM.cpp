#include "DisjointPairs3TrParam2DegPolyHMM.hpp"

/*
* only depths and samples to be used should be passed here. That is, if given 100 depths originally, but only 50 should be used
* here, then the depths passed here should have only the 50 to be used, since depths/sampleList is passed to the parent ctor
*/
DisjointPairs3TrParam2DegPolyHMM::DisjointPairs3TrParam2DegPolyHMM(std::vector<DepthPair*>* depths, std::vector<std::string>* sampleList, gsl_vector* fixedParams, int maxPloidy, int numFixedLibs, int numPairs, int numBranchesToEst) : AllPairs3TrParam2DegPolyHMM(depths, sampleList, fixedParams, maxPloidy, 2, 1, numFixedLibs, numPairs, numBranchesToEst) { // 2 shared transition params (beta/lambda), 1 fixed transition param (alpha), numFixedLibs fixed libs, numPairs pairs
  // everything is delegated to super ctor, there are no member variables
  //std::cout << "in Disjoint ctor; NUM_SHARED_TRANSITION_PARAMS_TO_EST = " << this->NUM_SHARED_TRANSITION_PARAMS_TO_EST << std::endl;
}
//DisjointPairs3TrParam2DegPolyHMM::DisjointPairs3TrParam2DegPolyHMM(const DisjointPairs3TrParam2DegPolyHMM& otherDisjointPairs3TrParam2DegPolyHMM) : AllPairs3TrParam2DegPolyHMM(otherDisjointPairs3TrParam2DegPolyHMM) {
//  // everything is delegated to super ctor, there are no member variables
//}
DisjointPairs3TrParam2DegPolyHMM* DisjointPairs3TrParam2DegPolyHMM::create(std::vector<DepthPair*>* depths, std::vector<std::string>* sampleList, int maxPloidy, int numPairs, gsl_vector* meanVarianceCoefVec, bool preallocIntermediates) {
  DisjointPairs3TrParam2DegPolyHMM* hmm = new DisjointPairs3TrParam2DegPolyHMM(depths, sampleList, nullptr, maxPloidy, 0, numPairs, numPairs * 3); // 0 fixedLibs, numPairs * 3 numBranchesToEst
  hmm->makeHMMPairs(meanVarianceCoefVec, preallocIntermediates);
  return hmm;
}
DisjointPairs3TrParam2DegPolyHMM* DisjointPairs3TrParam2DegPolyHMM::create(const DisjointPairs3TrParam2DegPolyHMM& otherDisjointPairs3TrParam2DegPolyHMM) {
  DisjointPairs3TrParam2DegPolyHMM* hmm = new DisjointPairs3TrParam2DegPolyHMM(otherDisjointPairs3TrParam2DegPolyHMM);
  //hmm->makeHMMPairs();
  hmm->hmmVec = otherDisjointPairs3TrParam2DegPolyHMM.hmmVec; // TODO not enough memory to copy?
  /*for(unsigned int i = 0; i < hmm->hmmVec->size(); i++) {
    (*hmm->hmmVec)[i] = hmm->copyHMM((*otherDisjointPairs3TrParam2DegPolyHMM.hmmVec)[i]);
  }*/
  return hmm;
}
DisjointPairs3TrParam2DegPolyHMM::~DisjointPairs3TrParam2DegPolyHMM() {
  // TODO
}
void DisjointPairs3TrParam2DegPolyHMM::makeHMMPairs(gsl_vector* meanVarianceCoefVec, gsl_vector* transitionParams, bool preallocIntermediates) {
  //std::cout << "in DisjointPairs3TrParam2DegPolyHMM::makeHMMPairs" << std::endl;
  std::vector<DepthPair*>* currDepths = nullptr;
  gsl_vector* currMeanVarCoefVec = nullptr; // make them all have their own copies of this vector
  int hmmIdx = 0;
  for(unsigned int i = 0; i < this->depthsVec->size() && hmmIdx < this->NUM_PAIRS; i+=2) {
      currDepths = new std::vector<DepthPair*>();
      currDepths->push_back((*this->depthsVec)[i]);
      currDepths->push_back((*this->depthsVec)[i+1]);
      TwoCell3TrParam2DegPolyHMM* hmm = new TwoCell3TrParam2DegPolyHMM(currDepths, this->getKploidy(), preallocIntermediates);
      (*this->hmmVec)[hmmIdx] = hmm;

      // do rest of HMM set up (usually happens in main.cpp)
      currMeanVarCoefVec = gsl_vector_alloc(meanVarianceCoefVec->size);
      gsl_vector_memcpy(currMeanVarCoefVec, meanVarianceCoefVec);
      hmm->setMeanVarianceFn(currMeanVarCoefVec);
      hmm->setTransition(transitionParams); // this vector isn't saved anywhere
      hmm->setLibScalingFactorsToTotalRatio();
      hmm->setAlpha(this->getAlpha());

      this->setLibScalingFactors(i, i+1, hmm->getLibScalingFactor(0), hmm->getLibScalingFactor(1));
      this->storeHMMIdxForCells(i, i+1, hmmIdx);
      //hmm->print(stdout);
      hmmIdx++;
  }
}

int DisjointPairs3TrParam2DegPolyHMM::getCell0IdxFromHMMIdx(int hmmIdx) {
  //std::cout <<" ###################### in DisjointPairs3TrParam2DegPolyHMM::getCell0IdxFromHMMIdx" << std::endl;
  return hmmIdx * 2;
}
int DisjointPairs3TrParam2DegPolyHMM::getCell1IdxFromHMMIdx(int hmmIdx) {
  return hmmIdx * 2 + 1;
}
int DisjointPairs3TrParam2DegPolyHMM::getHMMIdxFromCellPair(int cell0Idx, int cell1Idx) {
  // if cell0 and cell1 aren't adjacent to each other (ie if input (0,2))
  if(cell1Idx - cell0Idx != 1) {
    return -1;
  }
  // if cell0 isn't even (ie if input (1,2))
  if(cell0Idx % 2 != 0) {
    return -1;
  }
  return cell0Idx / 2;
}

DisjointPairs3TrParam2DegPolyHMM* DisjointPairs3TrParam2DegPolyHMM::bfgs(gsl_vector* initGuess, int maxIters, bool verbose, bool debug) {
  // create new HMM with the best guess parameters and return it
  //std::cout << "subclass::bfgs this paramsToEst: ";
  //printColVector(this->getParamsToEst());
  //DisjointPairs3TrParam2DegPolyHMM* bestGuessOptim = new DisjointPairs3TrParam2DegPolyHMM(*this);
  //DisjointPairs3TrParam2DegPolyHMM* bestGuessOptim = DisjointPairs3TrParam2DegPolyHMM::create(*this);
  DisjointPairs3TrParam2DegPolyHMM* bestGuessOptim = this;//DisjointPairs3TrParam2DegPolyHMM::create(*this);
  //bestGuessOptim->print(stdout);
  //std::cout << "subclass::bfgs bestGuessOptim paramsToEst: ";
  //printColVector(bestGuessOptim->getParamsToEst());
  gsl_vector* initGuessAsParams = gsl_vector_alloc(initGuess->size);
  //std::cout << "initGuess: " << std::endl;;
  //printColVector(initGuess);
  this->convertProbToParam(initGuessAsParams, initGuess);
  //std::cout << "initGuessAsParams: " << std::endl;
  //printColVector(initGuessAsParams);
  //exit(0);
  Optimizable::bfgs(initGuessAsParams, bestGuessOptim, maxIters, verbose, debug);
  gsl_vector_free(initGuessAsParams);
  return bestGuessOptim;
}

