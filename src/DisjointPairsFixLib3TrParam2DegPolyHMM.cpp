#include "DisjointPairsFixLib3TrParam2DegPolyHMM.hpp"

/*
* only depths and samples to be used should be passed here. That is, if given 100 depths originally, but only 50 should be used
* here, then the depths passed here should have only the 50 to be used, since depths/sampleList is passed to the parent ctor
*/
//DisjointPairs3TrParam2DegPolyHMM::DisjointPairs3TrParam2DegPolyHMM(std::vector<DepthPair*>* depths, std::vector<std::string>* sampleList, int maxPloidy, int numPairs) : AllPairs3TrParam2DegPolyHMM(depths, sampleList, maxPloidy, 2, 1, 0, numPairs) { // 2 shared transition params (beta/lambda), 1 fixed transition param (alpha), 0 fixed libs, numPairs pairs
DisjointPairsFixLib3TrParam2DegPolyHMM::DisjointPairsFixLib3TrParam2DegPolyHMM(std::vector<DepthPair*>* depths, std::vector<std::string>* sampleList, gsl_vector* fixedParams, int maxPloidy, int numPairs, int numBranchesToEst) : DisjointPairs3TrParam2DegPolyHMM(depths, sampleList, fixedParams, maxPloidy, depths->size(), numPairs, numBranchesToEst) { // all libs are fixed (ie number of cells)
  // everything is delegated to super ctor, there are no member variables
}
//DisjointPairsFixLib3TrParam2DegPolyHMM::DisjointPairsFixLib3TrParam2DegPolyHMM(const DisjointPairsFixLib3TrParam2DegPolyHMM& otherDisjointPairsFixLib3TrParam2DegPolyHMM) : DisjointPairs3TrParam2DegPolyHMM(otherDisjointPairsFixLib3TrParam2DegPolyHMM) {
//  // everything is delegated to super ctor, there are no member variables
//}
DisjointPairsFixLib3TrParam2DegPolyHMM* DisjointPairsFixLib3TrParam2DegPolyHMM::create(std::vector<DepthPair*>* depths, std::vector<std::string>* sampleList, gsl_vector* fixedParams, int maxPloidy, int numPairs, gsl_vector* meanVarianceCoefVec, bool preallocIntermediates) {
  DisjointPairsFixLib3TrParam2DegPolyHMM* hmm = new DisjointPairsFixLib3TrParam2DegPolyHMM(depths, sampleList, fixedParams, maxPloidy, numPairs, numPairs * 3); // numPairs * 3 numBranchesToEst
  hmm->makeHMMPairs(meanVarianceCoefVec, preallocIntermediates);
  return hmm;
}
DisjointPairsFixLib3TrParam2DegPolyHMM* DisjointPairsFixLib3TrParam2DegPolyHMM::create(const DisjointPairsFixLib3TrParam2DegPolyHMM& otherDisjointPairsFixLib3TrParam2DegPolyHMM) {
  DisjointPairsFixLib3TrParam2DegPolyHMM* hmm = new DisjointPairsFixLib3TrParam2DegPolyHMM(otherDisjointPairsFixLib3TrParam2DegPolyHMM);
  //hmm->makeHMMPairs();
  hmm->hmmVec = otherDisjointPairsFixLib3TrParam2DegPolyHMM.hmmVec; // TODO not enough memory to copy?
  return hmm;
}
DisjointPairsFixLib3TrParam2DegPolyHMM::~DisjointPairsFixLib3TrParam2DegPolyHMM() {
  // TODO
}
void DisjointPairsFixLib3TrParam2DegPolyHMM::makeHMMPairs(gsl_vector* meanVarianceCoefVec, gsl_vector* transitionParams, bool preallocIntermediates) {
  //std::cout << "in DisjointPairsFixLib3TrParam2DegPolyHMM::makeHMMPairs" << std::endl;
  std::vector<DepthPair*>* currDepths = nullptr;
  gsl_vector* currMeanVarCoefVec = nullptr; // make them all have their own copies of this vector
  int hmmIdx = 0;
  gsl_vector* currFixedParams = nullptr;
  for(unsigned int i = 0; i < this->depthsVec->size() && hmmIdx < this->NUM_PAIRS; i+=2) {
      currDepths = new std::vector<DepthPair*>();
      currDepths->push_back((*this->depthsVec)[i]);
      currDepths->push_back((*this->depthsVec)[i+1]);

      currFixedParams = gsl_vector_alloc(2 + 1); // 2 libs, alpha
      // set libs
      gsl_vector_set(currFixedParams, 0, gsl_vector_get(this->fixedParams, this->LIB_SIZE_SCALING_FACTOR_START_IDX + i));
      gsl_vector_set(currFixedParams, 1, gsl_vector_get(this->fixedParams, this->LIB_SIZE_SCALING_FACTOR_START_IDX + i+1));

      // set shared transition params
      gsl_vector_set(currFixedParams, 2, gsl_vector_get(this->fixedParams, this->FIXED_TRANSITION_PROB_START_IDX + 0)); // alpha
      //std::cout << "currFixedParams" << std::endl;
      //printColVector(currFixedParams);
      //printColVector(transitionParams);

      TwoCellFixLib3TrParam2DegPolyHMM* hmm = new TwoCellFixLib3TrParam2DegPolyHMM(currDepths, currFixedParams, this->getKploidy(), preallocIntermediates);
      (*this->hmmVec)[hmmIdx] = hmm;

      // do rest of HMM set up (usually happens in main.cpp)
      currMeanVarCoefVec = gsl_vector_alloc(meanVarianceCoefVec->size);
      gsl_vector_memcpy(currMeanVarCoefVec, meanVarianceCoefVec);
      hmm->setMeanVarianceFn(currMeanVarCoefVec);
      hmm->setTransition(transitionParams); // this vector isn't saved anywhere
      //hmm->setLibScalingFactorsToTotalRatio();
      hmm->setAlpha(this->getAlpha());

      this->setLibScalingFactors(i, i+1, hmm->getLibScalingFactor(0), hmm->getLibScalingFactor(1));
      this->storeHMMIdxForCells(i, i+1, hmmIdx);
      hmmIdx++;
  }
  //std::cout << "this->paramsToEst:" << std::endl;
  //printColVector(this->paramsToEst);
  //std::cout << "this->fixedParams:" << std::endl;
  //printColVector(this->fixedParams);
}

DisjointPairsFixLib3TrParam2DegPolyHMM* DisjointPairsFixLib3TrParam2DegPolyHMM::bfgs(gsl_vector* initGuess, int maxIters, bool verbose, bool debug) {
  // create new HMM with the best guess parameters and return it
  //DisjointPairsFixLib3TrParam2DegPolyHMM* bestGuessOptim = DisjointPairsFixLib3TrParam2DegPolyHMM::create(*this);
  DisjointPairsFixLib3TrParam2DegPolyHMM* bestGuessOptim = this;//DisjointPairsFixLib3TrParam2DegPolyHMM::create(*this);
  gsl_vector* initGuessAsParams = gsl_vector_alloc(initGuess->size);
  this->convertProbToParam(initGuessAsParams, initGuess);
  //printColVector(initGuessAsParams);
  Optimizable::bfgs(initGuessAsParams, bestGuessOptim, maxIters, verbose, debug);
  gsl_vector_free(initGuessAsParams);
  return bestGuessOptim;
}
//void DisjointPairsFixLib3TrParam2DegPolyHMM::setLibScalingFactor(int cellNum, double libScalingFactor) {
void DisjointPairsFixLib3TrParam2DegPolyHMM::setLibScalingFactors(int cell0Idx, int cell1Idx, double lib0, double lib1) {
  gsl_vector_set(this->fixedParams, this->LIB_SIZE_SCALING_FACTOR_START_IDX + cell0Idx, lib0);
  gsl_vector_set(this->fixedParams, this->LIB_SIZE_SCALING_FACTOR_START_IDX + cell1Idx, lib1);
  int hmmIdx = getHMMIdxFromCellPair(cell0Idx, cell1Idx);
  (*this->hmmVec)[hmmIdx]->setLibScalingFactor(0, lib0);
  (*this->hmmVec)[hmmIdx]->setLibScalingFactor(1, lib1);
}

double DisjointPairsFixLib3TrParam2DegPolyHMM::getLibScalingFactor(int cellNum) const {
  return gsl_vector_get(this->fixedParams, this->LIB_SIZE_SCALING_FACTOR_START_IDX + cellNum);
}

/*
 * sets all cells in indicated pair to passed libScalingFactor. Sets in this class, as well as each individual HMM
 * cellNumInPair should be either 0 or 1 (assumed, no checks)
 * overrides AllPairs3TrParam2DegPolyHMM since that one assumes libs go into paramsToEst
 */
 void DisjointPairsFixLib3TrParam2DegPolyHMM::setAllLibScalingFactors(int cellNumInPair, double libScalingFactor) {
  for(unsigned int i = 0; i < this->hmmVec->size(); i++) {
    gsl_vector_set(this->fixedParams, this->LIB_SIZE_SCALING_FACTOR_START_IDX + i * 2 + cellNumInPair, libScalingFactor);
    (*this->hmmVec)[i]->setLibScalingFactor(cellNumInPair, libScalingFactor);
  }
}

