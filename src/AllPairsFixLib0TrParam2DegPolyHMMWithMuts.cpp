#include "AllPairsFixLib0TrParam2DegPolyHMMWithMuts.hpp"

AllPairsFixLib0TrParam2DegPolyHMMWithMuts::AllPairsFixLib0TrParam2DegPolyHMMWithMuts(std::vector<DepthPair*>* depths, std::vector<std::string>* sampleList, std::vector<MutationList*>* mutListVec, gsl_vector* allIndEstMutCounts, gsl_vector* fixedParams, int maxPloidy, int numPairs, int numBranchesToEst) : AllPairsFixLib0TrParam2DegPolyHMM(depths, sampleList, fixedParams, maxPloidy, numPairs, numBranchesToEst) {
  this->mutListVec = mutListVec;
  this->mutPairVec = new std::vector<MutationPair*>(this->NUM_HMMS); // each TwoCell HMM has its own MutationPair
  this->allIndEstMutCounts = allIndEstMutCounts;
}
AllPairsFixLib0TrParam2DegPolyHMMWithMuts* AllPairsFixLib0TrParam2DegPolyHMMWithMuts::create(std::vector<DepthPair*>* depths, std::vector<std::string>* sampleList, std::vector<MutationList*>* mutListVec, gsl_vector* allIndEstMutCounts, gsl_vector* fixedParams, int maxPloidy, int numPairs, gsl_vector* meanVarianceCoefVec, bool preallocIntermediates) {
  AllPairsFixLib0TrParam2DegPolyHMMWithMuts* hmm = new AllPairsFixLib0TrParam2DegPolyHMMWithMuts(depths, sampleList, mutListVec, allIndEstMutCounts, fixedParams, maxPloidy, numPairs, numPairs * 3); // numPairs*3 numBranchesToEst
  hmm->makeHMMPairs(meanVarianceCoefVec, preallocIntermediates);
  return hmm;
}
void AllPairsFixLib0TrParam2DegPolyHMMWithMuts::makeOneHMMPair(int i, int j, bool preallocIntermediates) {
  int hmmIdx = this->getHMMIdxFromCellPair(i, j);
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
    //std::cout << "makeOneHMMPair i " << i << ", j " << j << std::endl;
    //printColVector(this->allIndEstMutCounts);
    double xi = gsl_vector_get(this->allIndEstMutCounts, i);
    double xj = gsl_vector_get(this->allIndEstMutCounts, j);
    double x1Est = std::min(xi, xj) / 2;
    double x2Est = xi - x1Est;
    double x3Est = xj - x1Est;
    gsl_vector_set(currInitGuess, 0, x1Est);
    gsl_vector_set(currInitGuess, 1, x2Est);
    gsl_vector_set(currInitGuess, 2, x3Est);
  }
  if(this->verbose) {
    std::lock_guard<std::mutex> lock(Optimizable::mtx); // from https://stackoverflow.com/a/18277334
    //std::cout << "\n#####\nCalling BFGS for MutationPair initialization (paired mutation counts) for HMM " << hmmIdx << "(" << (*this->sampleList)[i] << ", " << (*this->sampleList)[j] << "):" << std::endl;
    std::cout << "\n#####\nCalling SIMPLEX for MutationPair initialization (paired mutation counts) for HMM " << hmmIdx << "(" << (*this->sampleList)[i] << ", " << (*this->sampleList)[j] << "):" << std::endl;
    printColVector(currInitGuess);
  }
  MutationPair* mutPair = new MutationPair(currMutLists, currInitGuess, this->verbose, this->gradientDebug);
  gsl_vector_free(currInitGuess);

  gsl_vector* currFixedParams = gsl_vector_alloc(2 + 3 + 1); // 2 libs, alpha/beta/lambda, cnaToMutRateMu
  // set libs
  gsl_vector_set(currFixedParams, 0, gsl_vector_get(this->fixedParams, this->LIB_SIZE_SCALING_FACTOR_START_IDX + i));
  gsl_vector_set(currFixedParams, 1, gsl_vector_get(this->fixedParams, this->LIB_SIZE_SCALING_FACTOR_START_IDX + j));

  // set shared transition params
  gsl_vector_set(currFixedParams, 2, gsl_vector_get(this->fixedParams, this->FIXED_TRANSITION_PROB_START_IDX + 0)); // alpha
  gsl_vector_set(currFixedParams, 3, gsl_vector_get(this->fixedParams, this->FIXED_TRANSITION_PROB_START_IDX + 1)); // beta
  gsl_vector_set(currFixedParams, 4, gsl_vector_get(this->fixedParams, this->FIXED_TRANSITION_PROB_START_IDX + 2)); // lambda

  // set mutation params
  gsl_vector_set(currFixedParams, 5, gsl_vector_get(this->fixedParams, this->FIXED_MUTATION_PARAM_START_IDX + 0)); // cnaToMutRateMu

  TwoCellFixLib0TrParam2DegPolyHMMWithMuts* hmm = new TwoCellFixLib0TrParam2DegPolyHMMWithMuts(currDepths, mutPair, currFixedParams, this->getKploidy(), preallocIntermediates);

  (*this->hmmVec)[hmmIdx] = hmm;
  (*this->mutPairVec)[hmmIdx] = mutPair;

  // do rest of HMM set up
  gsl_vector* currMeanVarCoefVec = gsl_vector_alloc(this->meanVarianceCoefVec->size); // make them all have their own copies of this vector
  gsl_vector_memcpy(currMeanVarCoefVec, this->meanVarianceCoefVec);
  hmm->setMeanVarianceFn(currMeanVarCoefVec);
  hmm->setAlpha(this->getAlpha());

  this->setLibScalingFactors(i, j, hmm->getLibScalingFactor(0), hmm->getLibScalingFactor(1));
}

AllPairsFixLib0TrParam2DegPolyHMMWithMuts::~AllPairsFixLib0TrParam2DegPolyHMMWithMuts() {
  // TODO
}

///*
// * function to estimate mutation counts X1/X2/X3 for each MutationPair ==> actually, moved to MutationPair ctor
// */
//void AllPairsFixLib0TrParam2DegPolyHMMWithMuts::estimateAllMutationPairMutCounts(std::string filename) {
//  // for each MutationPair
//  gsl_vector* initGuess = gsl_vector_alloc(3);
//  for(unsigned int mutPairIdx = 0; mutPairIdx < this->mutPairVec->size; mutPairIdx++) {
//    // nullptr guard
//    if((*this->mutPairVec)[mutPairIdx] == nullptr) {
//      continue;
//    }
//    // call bfgs
//    (*this->mutPairVec)[mutPairIdx]->bfgs(
//  }
//  gsl_vector_free(initGuess);
//}
//

gsl_vector* AllPairsFixLib0TrParam2DegPolyHMMWithMuts::getAllPairedEstMutCounts() {
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

