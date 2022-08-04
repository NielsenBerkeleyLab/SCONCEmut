#include "AllPairs0TrParam2DegPolyHMMWithMuts.hpp"

AllPairs0TrParam2DegPolyHMMWithMuts::AllPairs0TrParam2DegPolyHMMWithMuts(std::vector<DepthPair*>* depths, std::vector<std::string>* sampleList, std::vector<MutationList*>* mutListVec, gsl_vector* allIndEstMutCounts, gsl_vector* fixedParams, int maxPloidy, int numFixedLibs, int numPairs, int numBranchesToEst) : AllPairs0TrParam2DegPolyHMM(depths, sampleList, fixedParams, maxPloidy, numFixedLibs, numPairs, numBranchesToEst) {
  this->mutListVec = mutListVec;
  this->allIndEstMutCounts = allIndEstMutCounts;
  this->mutPairVec = new std::vector<MutationPair*>(this->NUM_HMMS); // each TwoCell HMM has its own MutationPair
}
AllPairs0TrParam2DegPolyHMMWithMuts* AllPairs0TrParam2DegPolyHMMWithMuts::create(std::vector<DepthPair*>* depths, std::vector<std::string>* sampleList, std::vector<MutationList*>* mutListVec, gsl_vector* allIndEstMutCounts, gsl_vector* fixedParams, int maxPloidy, int numPairs, gsl_vector* meanVarianceCoefVec, bool preallocIntermediates) {
  AllPairs0TrParam2DegPolyHMMWithMuts* hmm = new AllPairs0TrParam2DegPolyHMMWithMuts(depths, sampleList, mutListVec, allIndEstMutCounts, fixedParams, maxPloidy, 0, numPairs, numPairs * 3); // 0 numFixedLibs, numPairs*3 numBranchesToEst
  hmm->makeHMMPairs(meanVarianceCoefVec, preallocIntermediates);
  return hmm;
}
void AllPairs0TrParam2DegPolyHMMWithMuts::makeOneHMMPair(int i, int j, bool preallocIntermediates) {
  std::vector<DepthPair*>* currDepths = new std::vector<DepthPair*>();
  currDepths->push_back((*this->depthsVec)[i]);
  currDepths->push_back((*this->depthsVec)[j]);

  std::vector<MutationList*>* currMutLists = new std::vector<MutationList*>();
  currMutLists->push_back((*this->mutListVec)[i]);
  currMutLists->push_back((*this->mutListVec)[j]);
  gsl_vector* currInitGuess = nullptr;
  if(this->allIndEstMutCounts != nullptr) {
    // similar to IndThenPairs2Stages3TrParam2DegPolyHMM::copyStage1EstsIntoStage2FixedParams tree branch estimate setting
    currInitGuess = gsl_vector_alloc(3);
    double xi = gsl_vector_get(this->allIndEstMutCounts, i);
    double xj = gsl_vector_get(this->allIndEstMutCounts, j);
    double x1Est = std::min(xi, xj) / 2;
    double x2Est = xi - x1Est;
    double x3Est = xj - x1Est;
    gsl_vector_set(currInitGuess, 0, x1Est);
    gsl_vector_set(currInitGuess, 1, x2Est);
    gsl_vector_set(currInitGuess, 2, x3Est);
  }
  MutationPair* mutPair = new MutationPair(currMutLists, currInitGuess, this->verbose, this->gradientDebug);
  gsl_vector_free(currInitGuess);

  gsl_vector* currFixedParams = gsl_vector_alloc(this->fixedParams->size);
  gsl_vector_memcpy(currFixedParams, this->fixedParams);
  TwoCell0TrParam2DegPolyHMMWithMuts* hmm = new TwoCell0TrParam2DegPolyHMMWithMuts(currDepths, mutPair, currFixedParams, this->getKploidy(), preallocIntermediates);
  int hmmIdx = this->getHMMIdxFromCellPair(i, j);
  (*this->hmmVec)[hmmIdx] = hmm;
  (*this->mutPairVec)[hmmIdx] = mutPair;

  // do rest of HMM set up (usually happens in main.cpp)
  gsl_vector* currMeanVarCoefVec = gsl_vector_alloc(this->meanVarianceCoefVec->size);
  gsl_vector_memcpy(currMeanVarCoefVec, this->meanVarianceCoefVec);
  hmm->setMeanVarianceFn(currMeanVarCoefVec);
  hmm->setLibScalingFactorsToTotalRatio();
  hmm->setAlpha(this->getAlpha());

  this->setLibScalingFactors(i, j, hmm->getLibScalingFactor(0), hmm->getLibScalingFactor(1));
}

AllPairs0TrParam2DegPolyHMMWithMuts::~AllPairs0TrParam2DegPolyHMMWithMuts() {
  // TODO
}

gsl_vector* AllPairs0TrParam2DegPolyHMMWithMuts::getAllPairedEstMutCounts() {
  gsl_vector* allPairedEstMutCounts = gsl_vector_alloc(this->NUM_BRANCH_LENGTHS_TO_EST);
  gsl_vector_set_all(allPairedEstMutCounts, -1);
  MutationPair* currMutPair = nullptr;
  int mutIdx = 0;
  for(unsigned int mutPairIdx = 0; mutPairIdx < this->mutPairVec->size(); mutPairIdx++) {
    currMutPair = (*this->mutPairVec)[mutPairIdx];
    if(currMutPair != nullptr) {
      gsl_vector_set(allPairedEstMutCounts, mutIdx + 0, currMutPair->getMutCountEst(0));
      gsl_vector_set(allPairedEstMutCounts, mutIdx + 1, currMutPair->getMutCountEst(1));
      gsl_vector_set(allPairedEstMutCounts, mutIdx + 2, currMutPair->getMutCountEst(2));
    }
    mutIdx +=3;
  }
  return allPairedEstMutCounts;
}

