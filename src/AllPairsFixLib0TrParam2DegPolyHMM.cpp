#include "AllPairsFixLib0TrParam2DegPolyHMM.hpp"

// ctors and destructor
AllPairsFixLib0TrParam2DegPolyHMM::AllPairsFixLib0TrParam2DegPolyHMM(std::vector<DepthPair*>* depths, std::vector<std::string>* sampleList, gsl_vector* fixedParams, int maxPloidy, int numPairs, int numBranchesToEst) : AllPairs0TrParam2DegPolyHMM(depths, sampleList, fixedParams, maxPloidy, depths->size(), numPairs, numBranchesToEst) { // all libs are fixed (ie number of cells)
  // track which HMM's should have bfgs called on them, but don't necessarily prealloc all of them ahead of time to save memory (ie can't depend on hmmVec != nullptr)
  this->shouldCallBFGSOnHmmIdx = new std::vector<bool>(numPairs, true);
}
AllPairsFixLib0TrParam2DegPolyHMM* AllPairsFixLib0TrParam2DegPolyHMM::create(std::vector<DepthPair*>* depths, std::vector<std::string>* sampleList, gsl_vector* fixedParams, int maxPloidy, int numPairs, gsl_vector* meanVarianceCoefVec, bool preallocIntermediates) {
  AllPairsFixLib0TrParam2DegPolyHMM* hmm = new AllPairsFixLib0TrParam2DegPolyHMM(depths, sampleList, fixedParams, maxPloidy, numPairs, numPairs * 3); // numPairs*3 numBranchesToEst
  hmm->makeHMMPairs(meanVarianceCoefVec, preallocIntermediates);
  return hmm;
}

/*
 * helper method to prep for HMM set up
 */
void AllPairsFixLib0TrParam2DegPolyHMM::makeHMMPairs(gsl_vector* meanVarianceCoefVec, gsl_vector* transitionParams, bool preallocIntermediates) {
  int hmmIdx = 0;
  bool atLeastOneHmmCreated = false;
  for(unsigned int i = 0; i < this->depthsVec->size() && hmmIdx < this->NUM_PAIRS; i++) {
    for(unsigned int j = i+1; j < this->depthsVec->size() && hmmIdx < this->NUM_PAIRS; j++) {

      // always create at least one hmm (needed for indexing constants)
      if(!atLeastOneHmmCreated) {
        this->makeOneHMMPair(i, j, preallocIntermediates);
        atLeastOneHmmCreated = true;
      }
      else if(preallocIntermediates) {
        this->makeOneHMMPair(i, j, preallocIntermediates);
      }
      hmmIdx = getHMMIdxFromCellPair(i, j);
      (*this->shouldCallBFGSOnHmmIdx)[hmmIdx] = true;
      this->storeHMMIdxForCells(i, j, hmmIdx);
    }
  }
}

void AllPairsFixLib0TrParam2DegPolyHMM::makeOneHMMPair(int i, int j, bool preallocIntermediates) {
  std::vector<DepthPair*>* currDepths = new std::vector<DepthPair*>();
  currDepths->push_back((*this->depthsVec)[i]);
  currDepths->push_back((*this->depthsVec)[j]);

  gsl_vector* currFixedParams = gsl_vector_alloc(2 + 3); // 2 libs, alpha/beta/lambda
  // set libs
  gsl_vector_set(currFixedParams, 0, gsl_vector_get(this->fixedParams, this->LIB_SIZE_SCALING_FACTOR_START_IDX + i));
  gsl_vector_set(currFixedParams, 1, gsl_vector_get(this->fixedParams, this->LIB_SIZE_SCALING_FACTOR_START_IDX + j));

  // set shared transition params
  gsl_vector_set(currFixedParams, 2, gsl_vector_get(this->fixedParams, this->FIXED_TRANSITION_PROB_START_IDX + 0)); // alpha
  gsl_vector_set(currFixedParams, 3, gsl_vector_get(this->fixedParams, this->FIXED_TRANSITION_PROB_START_IDX + 1)); // beta
  gsl_vector_set(currFixedParams, 4, gsl_vector_get(this->fixedParams, this->FIXED_TRANSITION_PROB_START_IDX + 2)); // lambda

  TwoCellFixLib0TrParam2DegPolyHMM* hmm = new TwoCellFixLib0TrParam2DegPolyHMM(currDepths, currFixedParams, this->getKploidy(), preallocIntermediates);

  int hmmIdx = this->getHMMIdxFromCellPair(i, j);
  (*this->hmmVec)[hmmIdx] = hmm;

  // do rest of HMM set up
  gsl_vector* currMeanVarCoefVec = gsl_vector_alloc(this->meanVarianceCoefVec->size); // make them all have their own copies of this vector
  gsl_vector_memcpy(currMeanVarCoefVec, this->meanVarianceCoefVec);
  hmm->setMeanVarianceFn(currMeanVarCoefVec);
  hmm->setAlpha(this->getAlpha());

  this->setLibScalingFactors(i, j, hmm->getLibScalingFactor(0), hmm->getLibScalingFactor(1));
}

AllPairsFixLib0TrParam2DegPolyHMM* AllPairsFixLib0TrParam2DegPolyHMM::create(const AllPairsFixLib0TrParam2DegPolyHMM& otherAllPairsHMM) {
  AllPairsFixLib0TrParam2DegPolyHMM* hmm = new AllPairsFixLib0TrParam2DegPolyHMM(otherAllPairsHMM);
  hmm->hmmVec = otherAllPairsHMM.hmmVec;
  return hmm;
}

AllPairsFixLib0TrParam2DegPolyHMM::~AllPairsFixLib0TrParam2DegPolyHMM() {
  // TODO
}

/*
 * given params from hmmIdx'th HMM, save into this class's paramsToEst. Does not save into the hmm itself.
 * assumes params = [t1, t2, t3]
 */
void AllPairsFixLib0TrParam2DegPolyHMM::setParamsToEstFromIthHMM(gsl_vector* params, int hmmIdx) {
  gsl_vector_set(this->paramsToEst, this->BRANCH_LENGTH_START_IDX + hmmIdx + 0, gsl_vector_get(params, 0)); // t1
  gsl_vector_set(this->paramsToEst, this->BRANCH_LENGTH_START_IDX + hmmIdx + 1, gsl_vector_get(params, 1)); // t2
  gsl_vector_set(this->paramsToEst, this->BRANCH_LENGTH_START_IDX + hmmIdx + 2, gsl_vector_get(params, 2)); // t3
}
/*
 * given params from hmmIdx'th HMM, save into this class's fixedParams. Does not save into the hmm itself.
 * assumes params = [lib0, lib1, alpha, beta, lambda]
 */
void AllPairsFixLib0TrParam2DegPolyHMM::setFixedParamsFromIthHMM(gsl_vector* params, int hmmIdx) {
  gsl_vector_set(this->fixedParams, this->LIB_SIZE_SCALING_FACTOR_START_IDX + hmmIdx + 0, gsl_vector_get(params, 0)); // lib0
  gsl_vector_set(this->fixedParams, this->LIB_SIZE_SCALING_FACTOR_START_IDX + hmmIdx + 1, gsl_vector_get(params, 1)); // lib1
  gsl_vector_set(this->fixedParams, this->FIXED_TRANSITION_PROB_START_IDX + 0, gsl_vector_get(params, 2)); // alpha
  gsl_vector_set(this->fixedParams, this->FIXED_TRANSITION_PROB_START_IDX + 1, gsl_vector_get(params, 3)); // beta
  gsl_vector_set(this->fixedParams, this->FIXED_TRANSITION_PROB_START_IDX + 2, gsl_vector_get(params, 4)); // lambda
}

/*
 * copied from AllPairs3TrParam2DegPolyHMM
 */
double AllPairsFixLib0TrParam2DegPolyHMM::setParamsToEst(gsl_vector* params) {
  gsl_vector_memcpy(this->paramsToEst, params); // Thu 30 Jan 2020 10:47:35 AM PST changed to be a memcpy instead of free/alloc cycle

  // nullptr guard
  HMM* hmm = this->getFirstNonNullHMM();
  int hmmBranchIdx = hmm->BRANCH_LENGTH_START_IDX;
  int numHMMParams = hmm->getNumParamsToEst();
  gsl_vector* currHMMParams = gsl_vector_alloc(numHMMParams);

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
    // get apropriate branch lengths (there are 3 branch lengths stored per HMM
    currT1BFGS = gsl_vector_get(params, this->BRANCH_LENGTH_START_IDX + 3 * hmmIdx + 0);
    currT2BFGS = gsl_vector_get(params, this->BRANCH_LENGTH_START_IDX + 3 * hmmIdx + 1);
    currT3BFGS = gsl_vector_get(params, this->BRANCH_LENGTH_START_IDX + 3 * hmmIdx + 2);

    // set everything into currHMMParams
    gsl_vector_set(currHMMParams, hmmBranchIdx + 0, currT1BFGS);
    gsl_vector_set(currHMMParams, hmmBranchIdx + 1, currT2BFGS);
    gsl_vector_set(currHMMParams, hmmBranchIdx + 2, currT3BFGS);

    // call setParamsToEst on subclass, summing return status for each one (GSL_SUCCESS = 0, as returned by HMM::findSteadyStateDist)
    status += (*this->hmmVec)[hmmIdx]->setParamsToEst(currHMMParams);
  }
  return status;

}

void AllPairsFixLib0TrParam2DegPolyHMM::setLibScalingFactors(int cell0Idx, int cell1Idx, double lib0, double lib1) {
  gsl_vector_set(this->fixedParams, this->LIB_SIZE_SCALING_FACTOR_START_IDX + cell0Idx, lib0);
  gsl_vector_set(this->fixedParams, this->LIB_SIZE_SCALING_FACTOR_START_IDX + cell1Idx, lib1);
  int hmmIdx = getHMMIdxFromCellPair(cell0Idx, cell1Idx);
  (*this->hmmVec)[hmmIdx]->setLibScalingFactor(0, lib0);
  (*this->hmmVec)[hmmIdx]->setLibScalingFactor(1, lib1);
}

double AllPairsFixLib0TrParam2DegPolyHMM::getLibScalingFactor(int cellNum) const {
  return gsl_vector_get(this->fixedParams, this->LIB_SIZE_SCALING_FACTOR_START_IDX + cellNum);
}
/*
 * sets all cells in indicated pair to passed libScalingFactor. Sets in this class, as well as each individual HMM
 * cellNumInPair should be either 0 or 1 (assumed, no checks)
 * overrides AllPairs3TrParam2DegPolyHMM since that one assumes libs go into paramsToEst
 */
 void AllPairsFixLib0TrParam2DegPolyHMM::setAllLibScalingFactors(int cellNumInPair, double libScalingFactor) {
  for(unsigned int i = 0; i < this->hmmVec->size(); i++) {
    gsl_vector_set(this->fixedParams, this->LIB_SIZE_SCALING_FACTOR_START_IDX + i * 2 + cellNumInPair, libScalingFactor);
    (*this->hmmVec)[i]->setLibScalingFactor(cellNumInPair, libScalingFactor);
  }
}

// BFGS
/*
 * Because there are no shared parameters between each HMM, we can easily parallelize here during each
 * HMM's bfgs call. Each call is independent
 * Mon 10 Jan 2022 02:23:43 PM PST this used to be an override, changed to be simple implementation to add filename for intermediate saving
 */
AllPairsFixLib0TrParam2DegPolyHMM* AllPairsFixLib0TrParam2DegPolyHMM::bfgs(gsl_vector* initGuess, std::string filename, int maxIters, bool verbose, bool debug) {
  std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now(); // for entire BFGS function
  HMM* currHMM = nullptr;
  boost::asio::thread_pool pool(this->numThreads);

  for(unsigned int hmmIdx = 0; hmmIdx < this->hmmVec->size(); hmmIdx++) {
    // if shouldn't call bfgs on this hmm, just skip to the next hmmIdx
    if(!(*this->shouldCallBFGSOnHmmIdx)[hmmIdx]) {
      continue;
    }
    // Submit a lambda object to the pool
    boost::asio::post(pool,
      [this, hmmIdx, initGuess, filename, maxIters, verbose, debug]() { // see https://stackoverflow.com/a/7627218 for lambda function syntax explanation
        this->callIndvHMMBFGS(hmmIdx, initGuess, filename, maxIters, verbose, debug);
    });
  }
  // Wait for all tasks in the pool to complete.
  pool.join();

  std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
  double elapsedSec = std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() / 1000000.0;
  double totalLik = -std::numeric_limits<double>::infinity();
  for(unsigned int hmmIdx = 0; hmmIdx < this->hmmVec->size(); hmmIdx++) {
    currHMM = (*this->hmmVec)[hmmIdx];
    // nullptr guard
    if(currHMM  == nullptr) {
      continue;
    }
    if(currHMM->finalLl != -std::numeric_limits<double>::infinity()) {
      totalLik += currHMM->finalLl;
    }
    else {
      std::cout << "WARNING: hmmIdx " << hmmIdx << " had a final loglikelihood of -Inf. Skipping." << std::endl;
    }
  }
  this->finalLl = totalLik;

  if(verbose) {
    if(totalLik != -std::numeric_limits<double>::infinity()) {
      printf("OPTIMIZED STAGE 2 TOTAL LOGLIKELIHOOD: %.40f\n", totalLik);
    }
    printf("STAGE 2 TOTAL TIME (sec): %.10f\n", elapsedSec);
  }
  if(debug) {
    std::cerr << "##################################################################" << std::endl; // add a buffer line to the err file
  }

  return this;
}

void AllPairsFixLib0TrParam2DegPolyHMM::callIndvHMMBFGS(int hmmIdx, gsl_vector* initGuess, std::string filename, int maxIters, bool verbose, bool debug) {
  int cell0Idx = getCell0IdxFromHMMIdx(hmmIdx);
  int cell1Idx = getCell1IdxFromHMMIdx(hmmIdx);

  // if this hmm should exist but doesn't yet, create it
  if((*this->hmmVec)[hmmIdx] == nullptr) {
    this->makeOneHMMPair(cell0Idx, cell1Idx, false);
    this->getOnePairedOptimParamsToEstFromFile(hmmIdx, this->pairedEstimatesPath, 2 + 3 + 3); // 2 libs + alpha/beta/lambda + 3 branches per pair. includes a file.exists check
  }

  if(verbose) {
    std::lock_guard<std::mutex> lock(Optimizable::mtx); // from https://stackoverflow.com/a/18277334
    if(verbose) {
      std::cout << "\nCalling BFGS for branch lengths on HMM " << hmmIdx << "(" << (*this->sampleList)[cell0Idx] << ", " << (*this->sampleList)[cell1Idx] << "):" << std::endl;
    }
    if(debug) {
      std::cerr << "##################################################################" << std::endl; // add a buffer line to the err file before each bfgs
    }
  }

  HMM* currHMM = (*this->hmmVec)[hmmIdx];
  // if the current HMM has been read from file, don't need to reestimate
  if(currHMM->hasBeenReadFromFile) {
    if(verbose) {
      std::lock_guard<std::mutex> lock(Optimizable::mtx); // from https://stackoverflow.com/a/18277334
      std::cout << "HMM " << hmmIdx << "(" << (*this->sampleList)[cell0Idx] << ", " << (*this->sampleList)[cell1Idx] << ") has been read from file, skipping BFGS parameter reestimation. Read params are:" << std::endl;
      printRowVector(currHMM->getParamsToEst());
    }
  }
  else {
    // extract relevant branch lengths to est
    gsl_vector* currInitGuess = gsl_vector_alloc(currHMM->getNumParamsToEst());
    for(unsigned int i = 0; i < currInitGuess->size; i++) {
      gsl_vector_set(currInitGuess, i, gsl_vector_get(initGuess, this->BRANCH_LENGTH_START_IDX + 3 * hmmIdx + i));
    }

    // call bfgs
    currHMM->bfgs(currInitGuess, maxIters, verbose, debug);

    // save paramEstimates to file
    std::string hmmName = (*this->hmmNames)[hmmIdx];
    boost::replace_all(hmmName, ",", "__");
    currHMM->saveParamEstimates(filename + "__" + hmmName + "__k" + std::to_string(this->MAX_PLOIDY) + ".pairedParams");

    gsl_vector_free(currInitGuess);
  }

  // save paramsToEst from hmm into this->paramsToEst
  gsl_vector* savedBestEstParams = this->getParamsToEst();
  gsl_vector* currBestEst = currHMM->getParamsToEst();
  for(unsigned int i = 0; i < currBestEst->size; i++) {
    gsl_vector_set(savedBestEstParams, 3 * hmmIdx + i, gsl_vector_get(currBestEst, i));
  }

  // save viterbi decoding
  currHMM->allocDecodingIntermediates();
  currHMM->viterbiDecode();

  // only keep minimal info to save memory
  currHMM->freeToMinimum();
}

