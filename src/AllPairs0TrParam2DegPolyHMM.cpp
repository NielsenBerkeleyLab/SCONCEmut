#include "AllPairs0TrParam2DegPolyHMM.hpp"

// ctors and destructor
//}
AllPairs0TrParam2DegPolyHMM::AllPairs0TrParam2DegPolyHMM(std::vector<DepthPair*>* depths, std::vector<std::string>* sampleList, gsl_vector* fixedParams, int maxPloidy, int numFixedLibs, int numPairs, int numBranchesToEst) : AllPairs3TrParam2DegPolyHMM(depths, sampleList, fixedParams, maxPloidy, 0, 3, numFixedLibs, numPairs, numBranchesToEst) { // 0 numSharedTrParamsToEst, 3 numFixedTrParams (alpha/beta/lambda), passed numFixedLibs, passed numPairs, passed numBranchesToEst
  this->maxNumBFGSStarts = 3;
}
AllPairs0TrParam2DegPolyHMM* AllPairs0TrParam2DegPolyHMM::create(std::vector<DepthPair*>* depths, std::vector<std::string>* sampleList, gsl_vector* fixedParams, int maxPloidy, int numPairs, gsl_vector* meanVarianceCoefVec, bool preallocIntermediates) {
  AllPairs0TrParam2DegPolyHMM* hmm = new AllPairs0TrParam2DegPolyHMM(depths, sampleList, fixedParams, maxPloidy, 0, numPairs, numPairs * 3); // 0 numFixedLibs, numPairs*3 numBranchesToEst
  hmm->makeHMMPairs(meanVarianceCoefVec, preallocIntermediates);
  return hmm;
}

/*
 * helper method to prep for HMM set up
 */
void AllPairs0TrParam2DegPolyHMM::makeHMMPairs(gsl_vector* meanVarianceCoefVec, bool preallocIntermediates) {
  // prep for HMM set up
  if(meanVarianceCoefVec == nullptr) {
    meanVarianceCoefVec = HMM::createMeanVarianceCoefVec();
  }
  if(this->meanVarianceCoefVec == nullptr) {
    this->meanVarianceCoefVec = gsl_vector_alloc(meanVarianceCoefVec->size);
  }
  gsl_vector_memcpy(this->meanVarianceCoefVec, meanVarianceCoefVec);

  // prep transition mat
  gsl_vector* transitionParams = gsl_vector_alloc(3);
  gsl_vector_set(transitionParams, 0, 0.05);  // t1
  gsl_vector_set(transitionParams, 1, 0.05);  // t2
  gsl_vector_set(transitionParams, 2, 0.05);  // t3

  // save into paramsToEst
  for(int i = 0; i < this->NUM_PAIRS; i++) {
    gsl_vector_set(this->paramsToEst, this->BRANCH_LENGTH_START_IDX + 0 + i * 3, gsl_vector_get(transitionParams, 0)); // t1
    gsl_vector_set(this->paramsToEst, this->BRANCH_LENGTH_START_IDX + 1 + i * 3, gsl_vector_get(transitionParams, 1)); // t2
    gsl_vector_set(this->paramsToEst, this->BRANCH_LENGTH_START_IDX + 2 + i * 3, gsl_vector_get(transitionParams, 2)); // t3
  }

  // process all pairs
  makeHMMPairs(meanVarianceCoefVec, transitionParams, preallocIntermediates);
  this->setHMMNames();
  gsl_vector_free(transitionParams);
}
void AllPairs0TrParam2DegPolyHMM::makeHMMPairs(gsl_vector* meanVarianceCoefVec, gsl_vector* transitionParams, bool preallocIntermediates) {
  int hmmIdx = 0;
  for(unsigned int i = 0; i < this->depthsVec->size() && hmmIdx < this->NUM_PAIRS; i++) {
    for(unsigned int j = i+1; j < this->depthsVec->size() && hmmIdx < this->NUM_PAIRS; j++) {
      this->makeOneHMMPair(i, j, preallocIntermediates);
      this->storeHMMIdxForCells(i, j, hmmIdx);
      hmmIdx++;
    }
  }
}
void AllPairs0TrParam2DegPolyHMM::makeOneHMMPair(int i, int j, bool preallocIntermediates) {
  std::vector<DepthPair*>* currDepths = new std::vector<DepthPair*>();
  currDepths->push_back((*this->depthsVec)[i]);
  currDepths->push_back((*this->depthsVec)[j]);
  gsl_vector* currFixedParams = gsl_vector_alloc(this->fixedParams->size);
  gsl_vector_memcpy(currFixedParams, this->fixedParams);
  TwoCell0TrParam2DegPolyHMM* hmm = new TwoCell0TrParam2DegPolyHMM(currDepths, currFixedParams, this->getKploidy(), preallocIntermediates);
  int hmmIdx = this->getHMMIdxFromCellPair(i, j);
  (*this->hmmVec)[hmmIdx] = hmm;

  // do rest of HMM set up (usually happens in main.cpp)
  gsl_vector* currMeanVarCoefVec = gsl_vector_alloc(this->meanVarianceCoefVec->size);
  gsl_vector_memcpy(currMeanVarCoefVec, this->meanVarianceCoefVec);
  hmm->setMeanVarianceFn(currMeanVarCoefVec);
  hmm->setLibScalingFactorsToTotalRatio();
  hmm->setAlpha(this->getAlpha());

  this->setLibScalingFactors(i, j, hmm->getLibScalingFactor(0), hmm->getLibScalingFactor(1));
}

AllPairs0TrParam2DegPolyHMM* AllPairs0TrParam2DegPolyHMM::create(const AllPairs0TrParam2DegPolyHMM& otherAllPairsHMM) {
  AllPairs0TrParam2DegPolyHMM* hmm = new AllPairs0TrParam2DegPolyHMM(otherAllPairsHMM);
  hmm->hmmVec = otherAllPairsHMM.hmmVec;
  return hmm;
}
AllPairs0TrParam2DegPolyHMM::~AllPairs0TrParam2DegPolyHMM() {
  // TODO
}

void AllPairs0TrParam2DegPolyHMM::setFixedParams(gsl_vector* params) {
  if(params != nullptr) {
    gsl_vector_memcpy(this->fixedParams, params);
  }

  // nullptr guard
  HMM* hmm = this->getFirstNonNullHMM();
  int numFixedHMMParams = hmm->getNumFixedParams();
  int hmmFixedLibIdx = hmm->LIB_SIZE_SCALING_FACTOR_START_IDX;
  int hmmFixedTrIdx = hmm->FIXED_TRANSITION_PROB_START_IDX;
  int hmmNumFixedTrParams = hmm->NUM_FIXED_TRANSITION_PARAMS;

  // call each HMM's setFixedParams, which calls setTransition
  int cell0Idx = -1;
  int cell1Idx = -1;
  double currLib0 = 0;
  double currLib1 = 0;
  gsl_vector* currHMMParams = gsl_vector_alloc(numFixedHMMParams);;
  for(unsigned int hmmIdx = 0; hmmIdx < this->hmmVec->size(); hmmIdx++) {
    // nullptr guard
    if((*this->hmmVec)[hmmIdx] == nullptr) {
      continue;
    }

    // copy fixed libs, if any
    if(this->NUM_FIXED_LIBS > 0) {
      // get cell0 and cell1 indices for libs
      cell0Idx = getCell0IdxFromHMMIdx(hmmIdx);
      cell1Idx = getCell1IdxFromHMMIdx(hmmIdx);

      currLib0 = this->getLibScalingFactor(cell0Idx);
      currLib1 = this->getLibScalingFactor(cell1Idx);

      gsl_vector_set(currHMMParams, hmmFixedLibIdx + 0, currLib0);
      gsl_vector_set(currHMMParams, hmmFixedLibIdx + 1, currLib1);
    }
    // copy fixed transition params
    for(int i = 0; i < hmmNumFixedTrParams; i++) {
      gsl_vector_set(currHMMParams, hmmFixedTrIdx + i, gsl_vector_get(this->fixedParams, this->FIXED_TRANSITION_PROB_START_IDX + i));
    }

    (*this->hmmVec)[hmmIdx]->setFixedParams(currHMMParams);
  }
}

gsl_vector* AllPairs0TrParam2DegPolyHMM::getFixedParams() {
  return this->fixedParams;
}
/*
 * copies transition parameters only; potential memory leak TODO
 */
void AllPairs0TrParam2DegPolyHMM::setFixedTrParams(gsl_vector* params) {
  for(unsigned int i = 0; i < params->size; i++) {
    gsl_vector_set(this->fixedParams, this->FIXED_TRANSITION_PROB_START_IDX + 1 + i, gsl_vector_get(params, i)); // offset by 1 for alpha; hardcoding this seems like a potentially bad idea TODO
  }

  // call setFixedParams so these new transition params will be distributed to the HMMs
  this->setFixedParams(nullptr);
}

/*
 * same as AllPairs3TrParam2DegPolyHMM::setAllTransition but skips over beta/lambda terms, and only takes branches
 * for this->paramsToEst and for initializing all HMMs to same branches
 * branches shouldbe [t1, t2, t3] (ie shared for all HMMs)
 */
double AllPairs0TrParam2DegPolyHMM::setAllBranches(gsl_vector* branches) {
  double t1 = gsl_vector_get(branches, 0);
  double t2 = gsl_vector_get(branches, 1);
  double t3 = gsl_vector_get(branches, 2);
  double status = 0;
  for(unsigned int hmmIdx = 0; hmmIdx < this->hmmVec->size(); hmmIdx++) {
    // nullptr guard
    if((*this->hmmVec)[hmmIdx] == nullptr) {
      continue;
    }
    gsl_vector_set(this->paramsToEst, this->BRANCH_LENGTH_START_IDX + hmmIdx * 3 + 0, t1);
    gsl_vector_set(this->paramsToEst, this->BRANCH_LENGTH_START_IDX + hmmIdx * 3 + 1, t2);
    gsl_vector_set(this->paramsToEst, this->BRANCH_LENGTH_START_IDX + hmmIdx * 3 + 2, t3);

    status += (*this->hmmVec)[hmmIdx]->setTransition(branches);
  }
  return status;
}

/*
 * given params from hmmIdx'th HMM, save into this class's paramsToEst. Does not save into the hmm itself.
 * assumes params = [lib0, lib1, t1, t2, t3]
 */
void AllPairs0TrParam2DegPolyHMM::setParamsToEstFromIthHMM(gsl_vector* params, int hmmIdx) {
  gsl_vector_set(this->paramsToEst, this->LIB_SIZE_SCALING_FACTOR_START_IDX + hmmIdx + 0, gsl_vector_get(params, 0)); // lib0
  gsl_vector_set(this->paramsToEst, this->LIB_SIZE_SCALING_FACTOR_START_IDX + hmmIdx + 1, gsl_vector_get(params, 1)); // lib1
  gsl_vector_set(this->paramsToEst, this->BRANCH_LENGTH_START_IDX + hmmIdx + 0, gsl_vector_get(params, 2)); // t1
  gsl_vector_set(this->paramsToEst, this->BRANCH_LENGTH_START_IDX + hmmIdx + 1, gsl_vector_get(params, 3)); // t2
  gsl_vector_set(this->paramsToEst, this->BRANCH_LENGTH_START_IDX + hmmIdx + 2, gsl_vector_get(params, 4)); // t3
}
/*
 * given params from hmmIdx'th HMM, save into this class's fixedParams. Does not save into the hmm itself.
 * assumes params = [alpha, beta, lambda]
 */
void AllPairs0TrParam2DegPolyHMM::setFixedParamsFromIthHMM(gsl_vector* params, int hmmIdx) {
  gsl_vector_set(this->fixedParams, this->FIXED_TRANSITION_PROB_START_IDX + 0, gsl_vector_get(params, 0)); // alpha
  gsl_vector_set(this->fixedParams, this->FIXED_TRANSITION_PROB_START_IDX + 1, gsl_vector_get(params, 1)); // beta
  gsl_vector_set(this->fixedParams, this->FIXED_TRANSITION_PROB_START_IDX + 2, gsl_vector_get(params, 2)); // lambda
}

// needed to compile for subclasses
gsl_vector* AllPairs0TrParam2DegPolyHMM::getAllPairedEstMutCounts() {
  return nullptr;
}

/*
 * copied from AllPairs3TrParam2DegPolyHMM
 */
double AllPairs0TrParam2DegPolyHMM::setParamsToEst(gsl_vector* params) {
  gsl_vector_memcpy(this->paramsToEst, params);

  HMM* hmm = this->getFirstNonNullHMM();
  int hmmLibIdx = hmm->LIB_SIZE_SCALING_FACTOR_START_IDX;
  int hmmBranchIdx = hmm->BRANCH_LENGTH_START_IDX;
  int numHMMParams = hmm->getNumParamsToEst();
  gsl_vector* currHMMParams = gsl_vector_alloc(numHMMParams);

  int cell0Idx = -1;
  int cell1Idx = -1;

  double currLib0BFGS = 0;
  double currLib1BFGS = 0;
  double currT1BFGS = 0;
  double currT2BFGS = 0;
  double currT3BFGS = 0;

  double status = 0;
  // for each HMM, call that HMM's setParamsToEst method
  for(unsigned int hmmIdx = 0; hmmIdx < this->hmmVec->size(); hmmIdx++) {
    // nullptr guard
    if((*this->hmmVec)[hmmIdx] == nullptr) {
      continue;
    }

    // get cell0 and cell1 indices for libs
    cell0Idx = getCell0IdxFromHMMIdx(hmmIdx);
    cell1Idx = getCell1IdxFromHMMIdx(hmmIdx);

    // get appropriate libs
    currLib0BFGS = gsl_vector_get(params, this->LIB_SIZE_SCALING_FACTOR_START_IDX + cell0Idx);
    currLib1BFGS = gsl_vector_get(params, this->LIB_SIZE_SCALING_FACTOR_START_IDX + cell1Idx);

    // get apropriate branch lengths (there are 3 branch lengths stored per HMM
    currT1BFGS = gsl_vector_get(params, this->BRANCH_LENGTH_START_IDX + 3 * hmmIdx + 0);
    currT2BFGS = gsl_vector_get(params, this->BRANCH_LENGTH_START_IDX + 3 * hmmIdx + 1);
    currT3BFGS = gsl_vector_get(params, this->BRANCH_LENGTH_START_IDX + 3 * hmmIdx + 2);

    // set everything into currHMMParams
    gsl_vector_set(currHMMParams, hmmLibIdx + 0, currLib0BFGS);
    gsl_vector_set(currHMMParams, hmmLibIdx + 1, currLib1BFGS);
    gsl_vector_set(currHMMParams, hmmBranchIdx + 0, currT1BFGS);
    gsl_vector_set(currHMMParams, hmmBranchIdx + 1, currT2BFGS);
    gsl_vector_set(currHMMParams, hmmBranchIdx + 2, currT3BFGS);

    // call setParamsToEst on subclass, summing return status for each one (GSL_SUCCESS = 0, as returned by HMM::findSteadyStateDist)
    status += (*this->hmmVec)[hmmIdx]->setParamsToEst(currHMMParams);
  }
  return status;

}
/*
 * library sizes should be >= 0
 * pairwise branch lengths should be >=0, <=1
 * copied from AllPairs3TrParam2DegPolyHMM
 */
void AllPairs0TrParam2DegPolyHMM::convertProbToParam(gsl_vector* dest, const gsl_vector* src) const {
  // lib scaling factors
  double r = 0;
  for(int cellIdx = 0; cellIdx < this->NUM_LIBS_TO_EST; cellIdx++) {
    r = gsl_vector_get(src, this->LIB_SIZE_SCALING_FACTOR_START_IDX + cellIdx);
    gsl_vector_set(dest, this->LIB_SIZE_SCALING_FACTOR_START_IDX + cellIdx, log(r));
  }

  // pairwise branch lengths
  double t1 = 0;
  double t2 = 0;
  double t3 = 0;
  for(int pairIdx = 0; pairIdx < this->NUM_PAIRS; pairIdx++) {
    t1 = gsl_vector_get(src, this->BRANCH_LENGTH_START_IDX + 3 * pairIdx + 0);
    t2 = gsl_vector_get(src, this->BRANCH_LENGTH_START_IDX + 3 * pairIdx + 1);
    t3 = gsl_vector_get(src, this->BRANCH_LENGTH_START_IDX + 3 * pairIdx + 2);
    gsl_vector_set(dest, this->BRANCH_LENGTH_START_IDX + 3 * pairIdx + 0, log(t1)); // set T1
    gsl_vector_set(dest, this->BRANCH_LENGTH_START_IDX + 3 * pairIdx + 1, log(t2)); // set T2
    gsl_vector_set(dest, this->BRANCH_LENGTH_START_IDX + 3 * pairIdx + 2, log(t3)); // set T3
  }
}
void AllPairs0TrParam2DegPolyHMM::convertParamToProb(gsl_vector* dest, const gsl_vector* src) const {
  // lib scaling factors
  double w = 0;
  for(int cellIdx = 0; cellIdx < this->NUM_LIBS_TO_EST; cellIdx++) {
    w = gsl_vector_get(src, this->LIB_SIZE_SCALING_FACTOR_START_IDX + cellIdx);
    gsl_vector_set(dest, this->LIB_SIZE_SCALING_FACTOR_START_IDX + cellIdx, exp(w));
  }

  // pairwise branch lengths
  double T1 = 0;
  double T2 = 0;
  double T3 = 0;
  for(int pairIdx = 0; pairIdx < this->NUM_PAIRS; pairIdx++) {
    T1 = gsl_vector_get(src, this->BRANCH_LENGTH_START_IDX + 3 * pairIdx + 0);
    T2 = gsl_vector_get(src, this->BRANCH_LENGTH_START_IDX + 3 * pairIdx + 1);
    T3 = gsl_vector_get(src, this->BRANCH_LENGTH_START_IDX + 3 * pairIdx + 2);
    gsl_vector_set(dest, this->BRANCH_LENGTH_START_IDX + 3 * pairIdx + 0, exp(T1)); // set t1
    gsl_vector_set(dest, this->BRANCH_LENGTH_START_IDX + 3 * pairIdx + 1, exp(T2)); // set t2
    gsl_vector_set(dest, this->BRANCH_LENGTH_START_IDX + 3 * pairIdx + 2, exp(T3)); // set t3
  }
}
/*
 * same as HMM::checkOptimProbValidity
 */
double AllPairs0TrParam2DegPolyHMM::checkOptimProbValidity(gsl_vector* probs) const {
  return 0; // Thu 27 Aug 2020 08:27:51 PM PDT debugging no validity check
  // shortcut for bad library sizes (too large or too small)
  double currLibSizeScalingFactor = -1;
  for(int i = 0; i < this->NUM_LIBS_TO_EST; i++) {
    currLibSizeScalingFactor = gsl_vector_get(probs, this->LIB_SIZE_SCALING_FACTOR_START_IDX + i);
    if(currLibSizeScalingFactor < 1e-2 || currLibSizeScalingFactor > 1e2) {
      return GSL_NAN;
    }
  }
  // shortcut for any probabilities becoming too small
  double probMin = gsl_vector_min(probs);
  if(probMin < 1e-5 || gsl_isnan(probMin)) {
    return GSL_NAN;
  }
  return 0;
}

// BFGS
/*
 * copied from AllPairs3TrParam2DegPolyHMM, but not transition probs
 */
void AllPairs0TrParam2DegPolyHMM::setInitGuessNthTime(gsl_vector* initGuess, int iter, int numTotalRuns) const {
  if(iter == 0) {
    // lib scaling factors
    for(int cellIdx = 0; cellIdx < this->NUM_LIBS_TO_EST; cellIdx++) {
      gsl_vector_set(initGuess, this->LIB_SIZE_SCALING_FACTOR_START_IDX + cellIdx, 1);
    }

    // pairwise branch lengths
    for(int pairIdx = 0; pairIdx < this->NUM_PAIRS; pairIdx++) {
      gsl_vector_set(initGuess, this->BRANCH_LENGTH_START_IDX + 3 * pairIdx + 0, 0.013); // set t1
      gsl_vector_set(initGuess, this->BRANCH_LENGTH_START_IDX + 3 * pairIdx + 1, 0.017); // set t2
      gsl_vector_set(initGuess, this->BRANCH_LENGTH_START_IDX + 3 * pairIdx + 2, 0.019); // set t3
    }
  }
  else if(iter == 1) {
    // lib scaling factors
    for(int cellIdx = 0; cellIdx < this->NUM_LIBS_TO_EST; cellIdx++) {
      gsl_vector_set(initGuess, this->LIB_SIZE_SCALING_FACTOR_START_IDX + cellIdx, .75);
    }

    // pairwise branch lengths
    for(int pairIdx = 0; pairIdx < this->NUM_PAIRS; pairIdx++) {
      gsl_vector_set(initGuess, this->BRANCH_LENGTH_START_IDX + 3 * pairIdx + 0, 0.1); // set t1
      gsl_vector_set(initGuess, this->BRANCH_LENGTH_START_IDX + 3 * pairIdx + 1, 0.1); // set t2
      gsl_vector_set(initGuess, this->BRANCH_LENGTH_START_IDX + 3 * pairIdx + 2, 0.1); // set t3
    }
  }
  else if(iter == 2) {
    // lib scaling factors
    for(int cellIdx = 0; cellIdx < this->NUM_LIBS_TO_EST; cellIdx++) {
      gsl_vector_set(initGuess, this->LIB_SIZE_SCALING_FACTOR_START_IDX + cellIdx, 1.25);
    }

    // pairwise branch lengths
    for(int pairIdx = 0; pairIdx < this->NUM_PAIRS; pairIdx++) {
      gsl_vector_set(initGuess, this->BRANCH_LENGTH_START_IDX + 3 * pairIdx + 0, 0.02); // set t1
      gsl_vector_set(initGuess, this->BRANCH_LENGTH_START_IDX + 3 * pairIdx + 1, 0.02); // set t2
      gsl_vector_set(initGuess, this->BRANCH_LENGTH_START_IDX + 3 * pairIdx + 2, 0.02); // set t3
    }
  }
}

/*
 * copied from AllPairs3TrParam2DegPolyHMM
 */
AllPairs0TrParam2DegPolyHMM* AllPairs0TrParam2DegPolyHMM::bfgs(gsl_vector* initGuess, int maxIters, bool verbose, bool debug) {
  // create new HMM with the best guess parameters and return it
  AllPairs0TrParam2DegPolyHMM* bestGuessOptim = this;//AllPairs0TrParam2DegPolyHMM::create(*this);
  gsl_vector* initGuessAsParams = gsl_vector_alloc(initGuess->size);
  this->convertProbToParam(initGuessAsParams, initGuess);
  Optimizable::bfgs(initGuessAsParams, bestGuessOptim, maxIters, verbose, debug);
  gsl_vector_free(initGuessAsParams);
  return bestGuessOptim;
}

/*
 * Mon 10 Jan 2022 02:29:34 PM PST added to make compiler happy about adding a filename to AllPairsFixLib0TrParam2DegPolyHMM
 */
AllPairs0TrParam2DegPolyHMM* AllPairs0TrParam2DegPolyHMM::bfgs(gsl_vector* initGuess, std::string filename, int maxIters, bool verbose, bool debug) {
  return this->bfgs(initGuess, maxIters, verbose, debug);
}

