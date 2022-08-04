#include "IndThenPairs2Stages3TrParam2DegPolyHMM.hpp"

// constructor and destructor
IndThenPairs2Stages3TrParam2DegPolyHMM::IndThenPairs2Stages3TrParam2DegPolyHMM(std::vector<DepthPair*>* depths, std::vector<std::string>* sampleList, int maxPloidy, bool runIndBFGS, int numPairsStage2, bool verbose, bool gradientDebug, bool centralDiff) :
                                        NUM_CELLS(depths->size()),
                                        NUM_SHARED_TRANSITION_PARAMS_TO_EST(2), // beta, lambda
                                        MAX_PLOIDY(maxPloidy),
                                        NUM_PAIRS_STAGE_2(numPairsStage2) { // create all possible pairs (n choose 2)
  this->depthsVec = depths;
  this->sampleList = sampleList;
  this->verbose = verbose;
  this->gradientDebug = gradientDebug;
  this->centralDiff = centralDiff;
  this->stage1PairwiseDistMat = nullptr;
  this->stage1NearestCellIdxMat = nullptr;

  this->runIndBFGS = runIndBFGS;
  this->ind2StagesVec = nullptr;
  this->stage2Pairs = nullptr;
  this->numThreads = 1;
}

IndThenPairs2Stages3TrParam2DegPolyHMM::~IndThenPairs2Stages3TrParam2DegPolyHMM() {
  // TODO
}

void IndThenPairs2Stages3TrParam2DegPolyHMM::cleanUpIndOptim() {
  // need to clean up the ind2StagesVec, but don't want to delete the DepthPairs (or any shared info). But, since I think they all hold pointers it should be ok
  if(this->ind2StagesVec != nullptr) {
    AllInd2Stages3TrParam2DegPolyHMM* currInd2Stage = nullptr;
    for(unsigned int i = 0; i < this->ind2StagesVec->size(); i++) {
      currInd2Stage = (*this->ind2StagesVec)[i];
      delete currInd2Stage;
    }
    this->ind2StagesVec = nullptr;
  }
}

std::vector<AllInd2Stages3TrParam2DegPolyHMM*>* IndThenPairs2Stages3TrParam2DegPolyHMM::getInd2StagesVec() {
  return this->ind2StagesVec;
}
AllPairs0TrParam2DegPolyHMM* IndThenPairs2Stages3TrParam2DegPolyHMM::getStage2HMM() {
  return this->stage2Pairs;
}

void IndThenPairs2Stages3TrParam2DegPolyHMM::print(FILE* stream) {
  if(this->ind2StagesVec != nullptr) {
    fprintf(stream, "\nIndThenPairs2Stages3TrParam2DegPolyHMM ind2StagesVec:\n");

    AllInd2Stages3TrParam2DegPolyHMM* currInd2Stage = nullptr;
    for(unsigned int i = 0; i < this->ind2StagesVec->size(); i++) {
      currInd2Stage = (*this->ind2StagesVec)[i];
      currInd2Stage->print(stream);
    }
  }

  if(this->stage2Pairs != nullptr) {
    fprintf(stream, "\nIndThenPairs2Stages3TrParam2DegPolyHMM stage2Pairs:\n");
    this->stage2Pairs->print(stream);
  }
}

std::vector<std::string>* IndThenPairs2Stages3TrParam2DegPolyHMM::getSampleList() {
  return this->sampleList;
}
void IndThenPairs2Stages3TrParam2DegPolyHMM::setSampleList(std::vector<std::string>* sampleList) {
  this->sampleList = sampleList;
}

int IndThenPairs2Stages3TrParam2DegPolyHMM::getIndInitGuessSize() const {
  // initGuess for the independent stages should be shared across all cells (ie one est for [beta, lambda, t])
  return (*this->ind2StagesVec)[0]->getBWInitGuessSize();
}

int IndThenPairs2Stages3TrParam2DegPolyHMM::getPairsInitGuessSize() const {
  // initGuess for the pairs stage is the number of params in the stage2 obj (includes all params that will be estimated, and already accounts for if libs are being estimated or not)
  return this->stage2Pairs->getNumParamsToEst();
}

/*
 * function to return a vector of the estimated libs from the independent stage.
 * entries are [libA, libB, ...]
 */
gsl_vector* IndThenPairs2Stages3TrParam2DegPolyHMM::getIndOptimAllLibEsts() const {
  gsl_vector* estLibs = gsl_vector_alloc(this->ind2StagesVec->size());
  double currLib = 0;
  AllInd2Stages3TrParam2DegPolyHMM* currInd2Stage = nullptr;
  AllInd3TrParam2DegPolyHMM* currIndHMM = nullptr;
  for(unsigned int i = 0; i < this->ind2StagesVec->size(); i++) {
    currInd2Stage = (*this->ind2StagesVec)[i];
    if(this->runIndBFGS) {
      currIndHMM = currInd2Stage->getBFGSAllInd();
    }
    else {
      currIndHMM = currInd2Stage->getBWAllInd();
    }
    currLib = currIndHMM->getLibScalingFactor(0); // always 0th cell
    gsl_vector_set(estLibs, i, currLib);
  }
  return estLibs;
}

/*
 * function to return a vector of the estimated ith transition param values from the independent stage.
 * beta is trParamIdx 0, lambda is trParamIdx1
 * so if trParamIdx == 0, entries in returned vector are [beta1, beta2, ...]
 */
gsl_vector* IndThenPairs2Stages3TrParam2DegPolyHMM::getIndOptimIthTrParams(int trParamIdx) const {
  gsl_vector* estTrParams = gsl_vector_alloc(this->ind2StagesVec->size());
  gsl_vector* currParamsToEst = nullptr;
  double currTrParam = 0;
  AllInd2Stages3TrParam2DegPolyHMM* currInd2Stage = nullptr;
  AllInd3TrParam2DegPolyHMM* currIndHMM = nullptr;
  for(unsigned int i = 0; i < this->ind2StagesVec->size(); i++) {
    currInd2Stage = (*this->ind2StagesVec)[i];
    if(this->runIndBFGS) {
      currIndHMM = currInd2Stage->getBFGSAllInd();
    }
    else {
      currIndHMM = currInd2Stage->getBWAllInd();
    }
    currParamsToEst = currIndHMM->getParamsToEst();
    currTrParam = gsl_vector_get(currParamsToEst, currIndHMM->SHARED_TRANSITION_PROB_START_IDX + trParamIdx);
    gsl_vector_set(estTrParams, i, currTrParam);
  }
  return estTrParams;
}

/*
 * function to return a vector of the summarized optimized transition params estimated from the independent stage.
 * entries are [median(beta), median(lambda)]
 */
gsl_vector* IndThenPairs2Stages3TrParam2DegPolyHMM::getIndOptimSummaryTrParams() const {
  gsl_vector* summaryEstTrParams = gsl_vector_alloc(this->NUM_SHARED_TRANSITION_PARAMS_TO_EST);
  double medianEst = 0;
  for(unsigned int trParamIdx = 0; trParamIdx < summaryEstTrParams->size; trParamIdx++) {
    gsl_vector* currTrParamEsts = this->getIndOptimIthTrParams(trParamIdx);
    gsl_sort(currTrParamEsts->data, 1, currTrParamEsts->size); // https://www.gnu.org/software/gsl/doc/html/sort.html#c.gsl_sort
    medianEst = gsl_stats_median_from_sorted_data(currTrParamEsts->data, 1, currTrParamEsts->size); // https://www.gnu.org/software/gsl/doc/html/statistics.html#c.gsl_stats_median_from_sorted_data
    gsl_vector_set(summaryEstTrParams, trParamIdx, medianEst);
    gsl_vector_free(currTrParamEsts);
  }
  return summaryEstTrParams;
}

/*
 * function to return a vector of the estimated branches from the independent stage.
 * entries are [tA, tB, ...]
 */
gsl_vector* IndThenPairs2Stages3TrParam2DegPolyHMM::getIndOptimAllBranchEsts() const {
  gsl_vector* estBranches = gsl_vector_alloc(this->ind2StagesVec->size());
  gsl_vector* currParamsToEst = nullptr;
  double currBranch = 0;
  AllInd2Stages3TrParam2DegPolyHMM* currInd2Stage = nullptr;
  AllInd3TrParam2DegPolyHMM* currIndHMM = nullptr;
  for(unsigned int i = 0; i < this->ind2StagesVec->size(); i++) {
    currInd2Stage = (*this->ind2StagesVec)[i];
    if(this->runIndBFGS) {
      currIndHMM = currInd2Stage->getBFGSAllInd();
    }
    else {
      currIndHMM = currInd2Stage->getBWAllInd();
    }
    currParamsToEst = currIndHMM->getParamsToEst();
    currBranch = gsl_vector_get(currParamsToEst, currIndHMM->BRANCH_LENGTH_START_IDX + 0); // always 0th branch
    gsl_vector_set(estBranches, i, currBranch);
  }
  return estBranches;
}

/*
 * function to return a vector of all the chrToViterbiPathMap's from sconce, in one vector
 * entries are [mapA, mapB, ...]
 */
std::vector<std::unordered_map<std::string, std::vector<int>*>*>* IndThenPairs2Stages3TrParam2DegPolyHMM::getIndOptimAllChrToViterbiPathMapVec() const {
  std::vector<std::unordered_map<std::string, std::vector<int>*>*>* allChrToViterbiPathMapVec = new std::vector<std::unordered_map<std::string, std::vector<int>*>*>(this->ind2StagesVec->size());
  std::vector<std::unordered_map<std::string, std::vector<int>*>*>* currChrToViterbiPathMapVec = nullptr;
  HMM* hmm = nullptr;
  for(unsigned int i = 0; i < this->ind2StagesVec->size(); i++) {
    if(this->runIndBFGS) {
      hmm = (*(*this->ind2StagesVec)[i]->getBFGSAllInd()->getHMMs())[0]; // assumes each AllInd obj only holds 1 cell (in the 0th position of hmmVec)
    } else {
      hmm = (*(*this->ind2StagesVec)[i]->getBWAllInd()->getHMMs())[0];
    }
    currChrToViterbiPathMapVec = hmm->getChrToViterbiPathMapVec();
    (*allChrToViterbiPathMapVec)[i] = (*currChrToViterbiPathMapVec)[0];
  }
  return allChrToViterbiPathMapVec;
}

/*
 * function to return a vector of the summarized optimized params estimated from the independent stage.
 * entries are [libA, libB, ..., alpha, median(beta), median(lambda), tA, tB, ...]
 */
gsl_vector* IndThenPairs2Stages3TrParam2DegPolyHMM::getIndOptimParamsToEstSummary() const {
  // get lib ests, summarized transition params, and branches
  gsl_vector* estLibs = this->getIndOptimAllLibEsts();
  gsl_vector* summaryEstTrParams = this->getIndOptimSummaryTrParams();
  gsl_vector* estBranches = this->getIndOptimAllBranchEsts();

  // copy into paramstoEstSummary
  gsl_vector* paramsToEstSummary = gsl_vector_alloc(estLibs->size + 1 + summaryEstTrParams->size + estBranches->size); // add one for alpha
  // libs
  unsigned int summaryIdx = 0;
  for(unsigned int libIdx = 0; libIdx < estLibs->size && summaryIdx < paramsToEstSummary->size; libIdx++, summaryIdx++) {
    gsl_vector_set(paramsToEstSummary, summaryIdx, gsl_vector_get(estLibs, libIdx));
  }

  // transition params
  gsl_vector_set(paramsToEstSummary, summaryIdx, (*this->ind2StagesVec)[0]->getBWAllInd()->getAlpha()); // alpha, from the first BW HMM
  summaryIdx++;
  for(unsigned int trParamIdx = 0; trParamIdx < summaryEstTrParams->size && summaryIdx < paramsToEstSummary->size; trParamIdx++, summaryIdx++) {
    gsl_vector_set(paramsToEstSummary, summaryIdx, gsl_vector_get(summaryEstTrParams, trParamIdx));
  }

  // branch lengths
  for(unsigned int branchIdx = 0; branchIdx < estBranches->size && summaryIdx < paramsToEstSummary->size; branchIdx++, summaryIdx++) {
    gsl_vector_set(paramsToEstSummary, summaryIdx, gsl_vector_get(estBranches, branchIdx));
  }

  // clean up
  gsl_vector_free(estLibs);
  gsl_vector_free(summaryEstTrParams);
  gsl_vector_free(estBranches);
  return paramsToEstSummary;
}

/*
 * function to set mean variance function for the independent cell estimation
 */
void IndThenPairs2Stages3TrParam2DegPolyHMM::setAllIndMeanVarianceFn(gsl_vector* meanVarianceCoefVec) {
  AllInd2Stages3TrParam2DegPolyHMM* currInd2Stage = nullptr;
  for(unsigned int i = 0; i < this->ind2StagesVec->size(); i++) {
    currInd2Stage = (*this->ind2StagesVec)[i];
    currInd2Stage->setBWAllMeanVarianceFn(meanVarianceCoefVec);
    currInd2Stage->setBFGSAllMeanVarianceFn(meanVarianceCoefVec);
  }
}

void IndThenPairs2Stages3TrParam2DegPolyHMM::setAllPairsMeanVarianceFn(gsl_vector* meanVarianceCoefVec) {
  this->stage2Pairs->setAllMeanVarianceFn(meanVarianceCoefVec);
}

double IndThenPairs2Stages3TrParam2DegPolyHMM::getTotalIndLogLikelihood() {
  double totalLogLik = 0;
  AllInd2Stages3TrParam2DegPolyHMM* currInd2Stage = nullptr;
  for(unsigned int i = 0; i < this->ind2StagesVec->size(); i++) {
    currInd2Stage = (*this->ind2StagesVec)[i];
    if(this->runIndBFGS) {
      totalLogLik += currInd2Stage->getBFGSAllIndLogLikelihood();
    }
    else {
      totalLogLik += currInd2Stage->getBWAllIndLogLikelihood();
    }
  }
  return totalLogLik;
}

void IndThenPairs2Stages3TrParam2DegPolyHMM::setNumThreads(int newNumThreads) {
  this->numThreads = newNumThreads;
  if(this->stage2Pairs != nullptr) {
    this->stage2Pairs->setNumThreads(newNumThreads);
  }
}

int IndThenPairs2Stages3TrParam2DegPolyHMM::getNumThreads() const {
  return this->numThreads;
}

double IndThenPairs2Stages3TrParam2DegPolyHMM::getStage2LogLikelihood() {
  return this->stage2Pairs->getLogLikelihood();
}

void IndThenPairs2Stages3TrParam2DegPolyHMM::setUpInd2Stages() {
  // this is similar to the set up in the sconce main function
  this->ind2StagesVec = new std::vector<AllInd2Stages3TrParam2DegPolyHMM*>();
  std::vector<DepthPair*>* currDepths = nullptr;
  std::vector<std::string>* currSample = nullptr;
  for(int cellIdx = 0; cellIdx < this->NUM_CELLS; cellIdx++) {
    currDepths = new std::vector<DepthPair*>(1);
    currSample = new std::vector<std::string>(1);

    // extract the current DepthPair and sample name
    (*currDepths)[0] = (*this->depthsVec)[cellIdx];
    (*currSample)[0] = (*this->sampleList)[cellIdx];

    // create and store AllInd2Stages3TrParam2DegPolyHMM obj
    AllInd2Stages3TrParam2DegPolyHMM* currInd2Stage = new AllInd2Stages3TrParam2DegPolyHMM(currDepths, currSample, this->MAX_PLOIDY, this->gradientDebug);
    currInd2Stage->setUpBaumWelch(); // TODO store meanVarianceCoefVec
    this->ind2StagesVec->push_back(currInd2Stage);
  }
}

/*
 * function to copy estimates from independent stage 1 into stage 2 vectors.
 * assumes all params are already allocated
 *
 * stage2FixedParams (fixed params for stage 2): [[libs,] alpha, beta, lambda], from indEstParams
 * stage2InitGuess (starting points for [libs and] branches): [[libs,] {t1,t2,t3}_AB, {t1,t2,t3}_AC,...]
 * indEstParams (estimated params from the independent stage; expected from getIndOptimParamsToEstSummary): [libA, libB, ..., alpha, median(beta), median(lambda), tA, tB, ...]
 * numStage2LibsToFix (how many libs should be stored into stage2FixedParams): as of Wed 25 Aug 2021 02:15:30 PM PDT, should be all, since we want all the stage 2 HMM's to be independent (ie parallelizable)
 */
void IndThenPairs2Stages3TrParam2DegPolyHMM::copyStage1EstsIntoStage2FixedParams(gsl_vector* stage2FixedParams, gsl_vector* stage2InitGuess, gsl_vector* indEstParams, int numStage2LibsToFix) {
  // store libs
  // if fixing libs
  unsigned int stage2FixedParamsIdx = 0;
  unsigned int stage2InitGuessIdx = 0;
  if(numStage2LibsToFix > 0) {
    for(int libIdx = 0; libIdx < numStage2LibsToFix && stage2FixedParamsIdx < stage2FixedParams->size; libIdx++, stage2FixedParamsIdx++) {
      gsl_vector_set(stage2FixedParams, stage2FixedParamsIdx, gsl_vector_get(indEstParams, libIdx));
    }
  }
  // else, estimating libs again, put into stage2InitGuess
  else {
    for(int libIdx = 0; libIdx < numStage2LibsToFix && stage2InitGuessIdx < stage2InitGuess->size; libIdx++, stage2InitGuessIdx++) {
      gsl_vector_set(stage2InitGuess, stage2InitGuessIdx, gsl_vector_get(indEstParams, libIdx));
    }
  }

  // store transition params
  //gsl_vector_set(stage2FixedParams, stage2FixedParamsIdx, (*this->ind2StagesVec)[0]->getBWAllInd()->getAlpha());// alpha, from the first BW HMM
  gsl_vector_set(stage2FixedParams, stage2FixedParamsIdx, gsl_vector_get(indEstParams, this->NUM_CELLS + 0));// alpha
  stage2FixedParamsIdx++;
  gsl_vector_set(stage2FixedParams, stage2FixedParamsIdx, gsl_vector_get(indEstParams, this->NUM_CELLS + 1)); // beta
  stage2FixedParamsIdx++;
  gsl_vector_set(stage2FixedParams, stage2FixedParamsIdx, gsl_vector_get(indEstParams, this->NUM_CELLS + 2));// lambda
  stage2FixedParamsIdx++;

  // store branches
  double tA = 0;
  double tB = 0;
  double t1Est = 0;
  double t2Est = 0;
  double t3Est = 0;
  for(unsigned int estBranchCounterA = 0; estBranchCounterA < (unsigned int) this->NUM_CELLS && this->NUM_CELLS + 1 + this->NUM_SHARED_TRANSITION_PARAMS_TO_EST + estBranchCounterA < indEstParams->size && stage2InitGuessIdx < stage2InitGuess->size; estBranchCounterA++) { // add 1 for alpha
    tA = gsl_vector_get(indEstParams, this->NUM_CELLS + 1 + this->NUM_SHARED_TRANSITION_PARAMS_TO_EST + estBranchCounterA);
    for(unsigned int estBranchCounterB = estBranchCounterA + 1; estBranchCounterB < (unsigned int) this->NUM_CELLS && this->NUM_CELLS + 1 + this->NUM_SHARED_TRANSITION_PARAMS_TO_EST + estBranchCounterB < indEstParams->size && stage2InitGuessIdx < stage2InitGuess->size; estBranchCounterB++) {
      tB = gsl_vector_get(indEstParams, this->NUM_CELLS + 1 + this->NUM_SHARED_TRANSITION_PARAMS_TO_EST + estBranchCounterB);
      t1Est = std::min(tA, tB) / 2;
      t2Est = tA - t1Est;
      t3Est = tB - t1Est;
      gsl_vector_set(stage2InitGuess, stage2InitGuessIdx + 0, t1Est);
      gsl_vector_set(stage2InitGuess, stage2InitGuessIdx + 1, t2Est);
      gsl_vector_set(stage2InitGuess, stage2InitGuessIdx + 2, t3Est);
      stage2InitGuessIdx += 3;
    }
  }
}

void IndThenPairs2Stages3TrParam2DegPolyHMM::setUpAllPairsBFGS(gsl_vector* fixedParams, gsl_vector* initGuess, gsl_vector* meanVarianceCoefVec, bool fixLib, bool preallocIntermediates) {
  if(fixLib) {
    this->stage2Pairs = AllPairsFixLib0TrParam2DegPolyHMM::create(this->depthsVec, this->sampleList, fixedParams, this->MAX_PLOIDY, this->NUM_PAIRS_STAGE_2, meanVarianceCoefVec, preallocIntermediates);
  }
  else {
    this->stage2Pairs = AllPairs0TrParam2DegPolyHMM::create(this->depthsVec, this->sampleList, fixedParams, this->MAX_PLOIDY, this->NUM_PAIRS_STAGE_2, meanVarianceCoefVec, preallocIntermediates);
  }
  this->stage2Pairs->gradientDebug = this->gradientDebug;
  this->stage2Pairs->setParamsToEst(initGuess);
  this->stage2Pairs->setNumThreads(this->numThreads);
  this->stage2Pairs->setCentralDiffFlag(centralDiff);
}

void IndThenPairs2Stages3TrParam2DegPolyHMM::setUpSelectPairsBFGS(gsl_vector* fixedParams, gsl_vector* initGuess, int numPairsToSummarize, int ordering, int seed, gsl_vector* meanVarianceCoefVec, bool fixLib, bool preallocIntermediates) {
  if(fixLib) {
    this->stage2Pairs = SelectPairsFixLib0TrParam2DegPolyHMM::create(this->depthsVec, this->sampleList, fixedParams, this->MAX_PLOIDY, this->NUM_PAIRS_STAGE_2, numPairsToSummarize, ordering, seed, this->stage1NearestCellIdxMat, meanVarianceCoefVec, preallocIntermediates);
  }
  else {
    this->stage2Pairs = SelectPairs0TrParam2DegPolyHMM::create(this->depthsVec, this->sampleList, fixedParams, this->MAX_PLOIDY, this->NUM_PAIRS_STAGE_2, numPairsToSummarize, ordering, seed, this->stage1NearestCellIdxMat, meanVarianceCoefVec, preallocIntermediates);
  }
  this->stage2Pairs->gradientDebug = this->gradientDebug;
  this->stage2Pairs->setParamsToEst(initGuess);
  this->stage2Pairs->setNumThreads(this->numThreads);
  this->stage2Pairs->setCentralDiffFlag(this->centralDiff);
}


// initGuess for the independent stages should be shared across all cells (ie one est for [beta, lambda, t])
void IndThenPairs2Stages3TrParam2DegPolyHMM::optimIndCells(gsl_vector* initGuess, std::string filename, int numBWIters, int numLibStarts, int libStartVal, int maxBFGSIters, bool verbose, bool debug) {
  // initGuess is used for doLeastSquares, and is the same size as paramsToEst for each bwIndHMM. But, the values in initGuess are overwritten and so the initial values don't matter. sconce's main fn saves the values of initGuess for bfgs

  std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

  // thread_pool code from https://www.boost.org/doc/libs/1_77_0/doc/html/boost_asio/reference/thread_pool.html
  boost::asio::thread_pool pool(this->numThreads);
  for(unsigned int hmmIdx = 0; hmmIdx < this->ind2StagesVec->size(); hmmIdx++) {
    // copy initGuess
    gsl_vector* currInitGuess = gsl_vector_alloc(initGuess->size);
    gsl_vector_memcpy(currInitGuess, initGuess);

    // Submit a lambda object to the pool
    boost::asio::post(pool,
      [this, hmmIdx, currInitGuess, numBWIters, numLibStarts, libStartVal, verbose, debug]() { // see https://stackoverflow.com/a/7627218 for lambda function syntax explanation
        this->callIndvBaumWelch(hmmIdx, currInitGuess, numBWIters, numLibStarts, libStartVal, verbose, debug);
    });
  }
  // Wait for all tasks in the pool to complete.
  pool.join();

  std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
  double bwTime = std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() / 1000000.0;
  printf("OPTIMINDCELLS BW+LS TOTAL TIME (sec): %.5f\n", bwTime);

  // for each cell, reest lib, transition params, and tree branch length using bfgs (the equivalent of the second half of sconce)
  if(this->runIndBFGS) {
    std::cout << "STARTING OPTIMINDCELLS BFGS SET UP" << std::endl;
    begin = std::chrono::steady_clock::now();

    // create a thread pool
    boost::asio::thread_pool pool(this->numThreads);
    gsl_vector* fixedParams = nullptr; // estimating everything with bfgs, keep fixedParams nullptr (alpha is set automatically in ctor)
    gsl_vector* meanVarianceCoefVec = (*(*this->ind2StagesVec)[0]->getBWAllInd()->getHMMs())[0]->getMeanVarianceCoefVec(); // get meanVar vec from 0'th HMM
    gsl_vector* bwResults = nullptr;
    gsl_vector* bfgsInitGuess = nullptr;
    for(unsigned int hmmIdx = 0; hmmIdx < this->ind2StagesVec->size(); hmmIdx++) {
      // first call set up if haven't already set up from reading file
      AllInd2Stages3TrParam2DegPolyHMM* currInd2Stage = (*this->ind2StagesVec)[hmmIdx];
      if(!currInd2Stage->hasHMMBeenReadFromFile()) {
        currInd2Stage->setUpAllIndBFGS(false, true, this->centralDiff, fixedParams); // fixLib=false (don't fix the libs, want to reest them); estTrParamsBFGS=true (reest tr params)

        // copy the mean/var relationship
        currInd2Stage->setBFGSAllMeanVarianceFn(meanVarianceCoefVec);

        // copy ests from bw+ls into bfgs initGuess
        bwResults = currInd2Stage->getBWAllInd()->getParamsToEst();
        bfgsInitGuess = gsl_vector_alloc(bwResults->size);
        gsl_vector_memcpy(bfgsInitGuess, bwResults);
        bool bwResultAdjusted = false;
        for(unsigned int i = 0; i < bfgsInitGuess->size; i++) {
          double bwResultValue = gsl_vector_get(bfgsInitGuess, i);
          if(bwResultValue < 1e-4) {
            bwResultAdjusted = true;
            bwResultValue = 1e-3;
            gsl_vector_set(bfgsInitGuess, i, bwResultValue);
          }
        }
        if(bwResultAdjusted) {
          std::cout << "adjusted bfgsInitGuess values < 1e-4 to be 1e-3" << std::endl;
          printColVector(bfgsInitGuess);
        }
        currInd2Stage->getBFGSAllInd()->setParamsToEst(bfgsInitGuess);
      }

      // call bfgs
      // Submit a lambda object to the pool
      boost::asio::post(pool,
        [this, hmmIdx, bfgsInitGuess, filename, maxBFGSIters, verbose, debug]() { // see https://stackoverflow.com/a/7627218 for lambda function syntax explanation
          this->callIndvBFGS(hmmIdx, bfgsInitGuess, filename, maxBFGSIters, verbose, debug);
      });

    }
    // Wait for all tasks in the pool to complete.
    pool.join();

    end = std::chrono::steady_clock::now();
    double bfgsTime = std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() / 1000000.0;
    printf("OPTIMINDCELLS BFGS TOTAL TIME (sec): %.5f\n", bfgsTime);
  }
}

/*
 * helper method for parallelizing bw+ls
 */
void IndThenPairs2Stages3TrParam2DegPolyHMM::callIndvBaumWelch(int hmmIdx, gsl_vector* initGuess, int numBWIters, int numLibStarts, double libStartVal, bool verbose, bool debug) {
  if(verbose) {
    std::lock_guard<std::mutex> lock(Optimizable::mtx); // from https://stackoverflow.com/a/18277334
    std::cout << "######## STARTING BAUM WELCH INITIALIZATION FOR IDX " << hmmIdx << " ########" << std::endl;
  }

  AllInd2Stages3TrParam2DegPolyHMM* currInd2Stage = (*this->ind2StagesVec)[hmmIdx];

  // if this hmm has already been read from file, don't run bw initialization
  if(currInd2Stage->hasHMMBeenReadFromFile()) {
    if(verbose) {
      std::lock_guard<std::mutex> lock(Optimizable::mtx); // from https://stackoverflow.com/a/18277334
      std::cout << "Cell " << hmmIdx << "(" << (*this->sampleList)[hmmIdx] << ") has been read from file, skipping Baum Welch initialization." << std::endl;
      printRowVector(currInd2Stage->getHMMParamsToEstReadFromFile());
    }
    return;
  }
  currInd2Stage->setBaumWelchInitGuess(initGuess, numBWIters, numLibStarts, libStartVal, verbose, debug);

  if(verbose) {
    std::lock_guard<std::mutex> lock(Optimizable::mtx); // from https://stackoverflow.com/a/18277334
    std::cout << "######## DONE WITH BAUM WELCH INITIALIZATION FOR IDX " << hmmIdx << " ########" << std::endl;
  }
}

/*
 * helper method for parallelizing indv bfgs
 */
void IndThenPairs2Stages3TrParam2DegPolyHMM::callIndvBFGS(int hmmIdx, gsl_vector* initGuess, std::string filename, int maxIters, bool verbose, bool debug) {
  if(verbose) {
    std::lock_guard<std::mutex> lock(Optimizable::mtx); // from https://stackoverflow.com/a/18277334
    std::cout << "######## STARTING INDV BFGS FOR IDX " << hmmIdx << " ########" << std::endl;
  }

  AllInd2Stages3TrParam2DegPolyHMM* currInd2Stage = (*this->ind2StagesVec)[hmmIdx];

  // if this hmm has already been read from file, don't run bfgs optimization
  if(currInd2Stage->hasHMMBeenReadFromFile()) {
    if(verbose) {
      std::lock_guard<std::mutex> lock(Optimizable::mtx); // from https://stackoverflow.com/a/18277334
      std::cout << "Cell " << hmmIdx << "(" << (*this->sampleList)[hmmIdx] << ") has been read from file, skipping BFGS parameter reestimation. Read params are:" << std::endl;
      printRowVector(currInd2Stage->getHMMParamsToEstReadFromFile());
    }
    return;
  }
  currInd2Stage->callBFGSNTimes(initGuess, 1, maxIters, verbose, debug);

  // save params
  std::string hmmName = (*this->sampleList)[hmmIdx];
  boost::replace_all(hmmName, ",", "__");
  currInd2Stage->saveParamEstimates(0, filename + "__" + hmmName + "__k" + std::to_string(this->MAX_PLOIDY) + ".sconceParams"); // each currInd2Stage holds 1 cell, so always working with 0'th cellIdx

  if(verbose) {
    std::lock_guard<std::mutex> lock(Optimizable::mtx); // from https://stackoverflow.com/a/18277334
    std::cout << "######## DONE WITH INDV BFGS FOR IDX " << hmmIdx << " ########" << std::endl;
  }
}
AllPairs0TrParam2DegPolyHMM* IndThenPairs2Stages3TrParam2DegPolyHMM::bfgsStage2(gsl_vector* initGuess, int maxIters, std::string filename, bool verbose, bool debug) {
  this->stage2Pairs = this->stage2Pairs->bfgs(initGuess, filename, maxIters, verbose, debug);
  return this->stage2Pairs;
}

/*
 * function to calc pariwise euclidean distanace between 2 cells.
 * cell0Idx and cell1Idx index into ind2StagesVec
 */
double IndThenPairs2Stages3TrParam2DegPolyHMM::calcPairwiseDist(int cell0Idx, int cell1Idx) {
  HMM* hmm0 = nullptr;
  HMM* hmm1 = nullptr;
  if(this->runIndBFGS) {
    hmm0 = (*(*this->ind2StagesVec)[cell0Idx]->getBFGSAllInd()->getHMMs())[0]; // assumes each AllInd obj only holds 1 cell (in the 0th position of hmmVec)
    hmm1 = (*(*this->ind2StagesVec)[cell1Idx]->getBFGSAllInd()->getHMMs())[0];
  }
  else {
    hmm0 = (*(*this->ind2StagesVec)[cell0Idx]->getBWAllInd()->getHMMs())[0];
    hmm1 = (*(*this->ind2StagesVec)[cell1Idx]->getBWAllInd()->getHMMs())[0];
  }

  double dist = 0;
  double chrDist = 0;
  std::vector<std::string>* chrVec = hmm0->getChrVec();
  std::string currChr;
  std::unordered_map<std::string, std::vector<int>*>* hmm0VitMap = (*hmm0->getChrToViterbiPathMapVec())[0];
  std::unordered_map<std::string, std::vector<int>*>* hmm1VitMap = (*hmm1->getChrToViterbiPathMapVec())[0];
  std::vector<int>* currPath0 = nullptr;
  std::vector<int>* currPath1 = nullptr;
  for(unsigned int chrIdx = 0; chrIdx < chrVec->size(); chrIdx++) {
    currChr = (*chrVec)[chrIdx];
    currPath0 = (*hmm0VitMap)[currChr];
    currPath1 = (*hmm1VitMap)[currChr];
    chrDist = calcEuclideanDistOfVectors(currPath0, currPath1);
    dist += pow(chrDist, 2.0); // want to resquare to sum up across chr, then sqrt genome wide
  }
  return sqrt(dist);
}

/*
 * function to alloc and create a symmetric pairwise distance matrix, stored in stage1PairwiseDistMat
 * dimensions are n x n
 *       | cell0 | cell1 | ... | celln
 * cell0
 * cell1
 * ...
 * celln
 */
void IndThenPairs2Stages3TrParam2DegPolyHMM::calcStage1PairwiseDistMat() {
  this->stage1PairwiseDistMat = gsl_matrix_alloc(this->NUM_CELLS, this->NUM_CELLS);
  double currDist = -1;
  for(int i = 0; i < this->NUM_CELLS; i++) {
    for(int j = i+1; j < this->NUM_CELLS; j++) {
      currDist = this->calcPairwiseDist(i, j);
      gsl_matrix_set(this->stage1PairwiseDistMat, i, j, currDist);
      gsl_matrix_set(this->stage1PairwiseDistMat, j, i, currDist); // set in both places to make symmetric
    }
  }
  std::cout << "this->stage1PairwiseDistMat" << std::endl;
  for(unsigned int row = 0; row < this->stage1PairwiseDistMat->size1; row++) {
    std::cout << "row: " << row << "\t";
    for(unsigned int col = 0; col < this->stage1PairwiseDistMat->size2; col++) {
      std::cout << gsl_matrix_get(this->stage1PairwiseDistMat, row, col) << "\t";
    }
    std::cout << std::endl;
  }
}

/*
 * function to alloc and create a matrix of which indices are closest to each cell, stored in stage1NearestCellIdxMat
 * assumes calcStage1PairwiseDistMat has already been called
 * dimensions are n x (n-1) (ie excludes self)
 *       | nearest cell idx to 0 | next nearest idx | ... | furthest idx
 * cell0
 * cell1
 * ...
 * celln
 */
void IndThenPairs2Stages3TrParam2DegPolyHMM::calcStage1NearestCellIdxMat() {
  gsl_permutation* perm = gsl_permutation_alloc(this->NUM_CELLS);
  int permVal = -1;
  this->stage1NearestCellIdxMat = gsl_matrix_alloc(this->NUM_CELLS, this->NUM_CELLS - 1);

  for(int cellRowIdx = 0; cellRowIdx < this->NUM_CELLS; cellRowIdx++) {
    // extract the vector to sort
    gsl_vector_view currRow = gsl_matrix_row(this->stage1PairwiseDistMat, cellRowIdx);

    // sort vector for indices
    // https://www.gnu.org/software/gsl/doc/html/sort.html#c.gsl_sort_vector_index
    // https://www.gnu.org/software/gsl/doc/html/sort.html#examples
    gsl_sort_vector_index(perm, &currRow.vector);

    // copy into stage1NearestCellIdxMat
    for(unsigned int permIdx = 0, cellColIdx = 0; permIdx < perm->size;) {
      permVal = gsl_permutation_get(perm, permIdx);

      // don't save self
      if(permVal == cellRowIdx) {
        permIdx++;
        continue;
      }

      gsl_matrix_set(this->stage1NearestCellIdxMat, cellRowIdx, cellColIdx, permVal);
      cellColIdx++;
      permIdx++;
    }
  }

  std::cout << "this->stage1NearestCellIdxMat" << std::endl;
  for(unsigned int row = 0; row < this->stage1NearestCellIdxMat->size1; row++) {
    std::cout << "row: " << row << "\t";
    for(unsigned int col = 0; col < this->stage1NearestCellIdxMat->size2; col++) {
      std::cout << gsl_matrix_get(this->stage1NearestCellIdxMat, row, col) << "\t";
    }
    std::cout << std::endl;
  }
}

void IndThenPairs2Stages3TrParam2DegPolyHMM::viterbiDecodeStage1() {
  for(unsigned int hmmIdx = 0; hmmIdx < this->ind2StagesVec->size(); hmmIdx++) {
    AllInd2Stages3TrParam2DegPolyHMM* currInd2Stage = (*this->ind2StagesVec)[hmmIdx];
    currInd2Stage->viterbiDecode();
  }
}

void IndThenPairs2Stages3TrParam2DegPolyHMM::viterbiDecodeStage2() {
  this->stage2Pairs->viterbiDecodeAll();
}

void IndThenPairs2Stages3TrParam2DegPolyHMM::saveStage1ViterbiDecodedCNA(std::string filename) {
  for(unsigned int hmmIdx = 0; hmmIdx < this->ind2StagesVec->size(); hmmIdx++) {
    AllInd2Stages3TrParam2DegPolyHMM* currInd2Stage = (*this->ind2StagesVec)[hmmIdx];
    // don't overwrite if read from file
    if(!currInd2Stage->hasHMMBeenReadFromFile()) {
      currInd2Stage->saveViterbiDecodedCNA(filename);
    }
  }
}

void IndThenPairs2Stages3TrParam2DegPolyHMM::saveStage2ViterbiDecodedCNA(std::string filename) {
  this->stage2Pairs->saveAllViterbiDecodedCNA(filename);
}

void IndThenPairs2Stages3TrParam2DegPolyHMM::saveStage1CNAToBed(std::string filename) {
  for(unsigned int hmmIdx = 0; hmmIdx < this->ind2StagesVec->size(); hmmIdx++) {
    AllInd2Stages3TrParam2DegPolyHMM* currInd2Stage = (*this->ind2StagesVec)[hmmIdx];
    currInd2Stage->saveCNAToBed(filename);
  }
}

void IndThenPairs2Stages3TrParam2DegPolyHMM::saveStage2CNAToBed(std::string filename) {
  this->stage2Pairs->saveAllCNAToBed(filename);
}

void IndThenPairs2Stages3TrParam2DegPolyHMM::saveStage1ParamEstimates(std::string filename) {
  for(unsigned int hmmIdx = 0; hmmIdx < this->ind2StagesVec->size(); hmmIdx++) {
    AllInd2Stages3TrParam2DegPolyHMM* currInd2Stage = (*this->ind2StagesVec)[hmmIdx];
    std::string hmmName = (*this->sampleList)[hmmIdx];
    boost::replace_all(hmmName, ",", "__");
    currInd2Stage->saveParamEstimates(0, filename + "__" + hmmName + "__k" + std::to_string(this->MAX_PLOIDY) + ".sconceParams"); // each currInd2Stage holds 1 cell, so always working with 0'th cellIdx
  }

  // also save joint file block
  gsl_vector* optimizedIndParamSummary = this->getIndOptimParamsToEstSummary();
  FILE* outFile = fopen((filename + "__k" + std::to_string(this->MAX_PLOIDY) + ".summarizedSconceParams").c_str(), "w");
  gsl_block_fprintf(outFile, optimizedIndParamSummary->block, "%.40f");
  fclose(outFile);
}

void IndThenPairs2Stages3TrParam2DegPolyHMM::summaryDecodeAllCellsAcrossPairs(int summaryMethod, bool runDecoding) {
  this->stage2Pairs->summaryDecodeAllCellsAcrossPairs(summaryMethod, runDecoding);
}

void IndThenPairs2Stages3TrParam2DegPolyHMM::saveAllSummaryDecodedCNAToBed(std::string filename) {
  this->stage2Pairs->saveAllSummaryDecodedCNAToBed(filename);
}

/*
 * function to read in individual cell .sconceParams file (ie from running just sconce alone) and combine them into one vector
 * assumes each file is in order of [lib, alpha, beta, lambda, t]
 * entries are [libA, libB, ..., alpha, median(beta), median(lambda), tA, tB, ...]
 */
void IndThenPairs2Stages3TrParam2DegPolyHMM::getIndOptimParamsToEstFromIndvFiles(std::string tumorFileList, std::string sconceEstimatesPath, int numExpectedLines, gsl_vector* meanVarianceCoefVec) const {
  for(unsigned int cellIdx = 0; cellIdx < this->depthsVec->size(); cellIdx++) {
    std::string hmmName = (*this->sampleList)[cellIdx];
    boost::replace_all(hmmName, ",", "__");
    std::string currFile = sconceEstimatesPath + "__" + hmmName + "__k" + std::to_string(this->MAX_PLOIDY) + ".sconceParams"; // based on saveStage1ParamEstimates() naming conventions

    AllInd2Stages3TrParam2DegPolyHMM* currInd2Stage = (*this->ind2StagesVec)[cellIdx];

    currInd2Stage->getParamsToEstFromFile(0, currFile, numExpectedLines, meanVarianceCoefVec); // each currInd2Stage holds 1 cell, so always working with 0'th cellIdx. subclass will set hmm's hasBeenReadFromFile
  }
}

/*
 * function to read in a summarized .summarizedSconceParams file (ie from running just sconce alone) and store them directly into a gsl_vector
 * Note: because this is a summarized file, cannot run any of the find nearest cell stuff that depends on having indv param estimates
 * entries are [libA, libB, ..., alpha, median(beta), median(lambda), tA, tB, ...]
 */
gsl_vector* IndThenPairs2Stages3TrParam2DegPolyHMM::getIndOptimParamsToEstSummaryFromJointFile(std::string sconceEstimatesJointFile, int numExpectedLines) const {
  gsl_vector* paramsToEstSummary = nullptr;

  // count how many lines are in sconceEstimatesJointFile
  std::ifstream jointFileStream(sconceEstimatesJointFile);
  jointFileStream.unsetf(std::ios_base::skipws);
  int numLines = std::count(std::istream_iterator<char>(jointFileStream), std::istream_iterator<char>(), '\n'); // line counting from https://stackoverflow.com/a/3482093
  jointFileStream.close();

  // check have correct number of lines
  if(numLines != numExpectedLines) {
    std::cerr << "Error: expected to read " << numExpectedLines << " lines from " << sconceEstimatesJointFile << ", but detected " << numLines << " lines. Exiting." << std::endl;
    exit(EXIT_FAILURE);
  }

  // read in sconceEstimatesJointFile into a vector
  paramsToEstSummary = gsl_vector_alloc(numLines);
  FILE* jointFile = fopen(sconceEstimatesJointFile.c_str(), "r");
  gsl_block_fscanf(jointFile, paramsToEstSummary->block);
  fclose(jointFile);

  return paramsToEstSummary;
}

/*
 * function for reading in paired param estimates from individual files (ie. ran part of the paired stuff but then had to kill it for whatever reason). If doesn't exist, skips
 */
void IndThenPairs2Stages3TrParam2DegPolyHMM::getPairedOptimParamsToEstFromFiles(std::string pairedEstimatesPath, int numExpectedLinesPerFile) {
  this->stage2Pairs->getPairedOptimParamsToEstFromFiles(pairedEstimatesPath, numExpectedLinesPerFile);
}

