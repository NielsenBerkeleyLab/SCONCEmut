#include "AllInd2Stages3TrParam2DegPolyHMM.hpp"

// ctors and destructors
AllInd2Stages3TrParam2DegPolyHMM::AllInd2Stages3TrParam2DegPolyHMM(std::vector<DepthPair*>* depths, std::vector<std::string>* sampleList, int maxPloidy, bool gradientDebug) :
                                        NUM_CELLS(depths->size()),
                                        NUM_SHARED_TRANSITION_PARAMS_TO_EST(2), // beta, lambda
                                        MAX_PLOIDY(maxPloidy) {

  this->depthsVec = depths;
  this->sampleList = sampleList;
  this->gradientDebug = gradientDebug;

  this->bwAllInd = nullptr;
  this->bfgsAllInd = nullptr;
}
AllInd2Stages3TrParam2DegPolyHMM::~AllInd2Stages3TrParam2DegPolyHMM() {
  delete this->bwAllInd;
  delete this->bfgsAllInd;
}

// accessors and mutators
int AllInd2Stages3TrParam2DegPolyHMM::getKploidy() const {
  return this->MAX_PLOIDY;
}

void AllInd2Stages3TrParam2DegPolyHMM::setUpBaumWelch() {
  this->bwAllInd = AllInd3TrParam2DegPolyHMM::create(this->depthsVec, this->sampleList, this->MAX_PLOIDY);
  this->bwAllInd->gradientDebug = this->gradientDebug;
}

void AllInd2Stages3TrParam2DegPolyHMM::setUpAllIndBFGS(bool fixLib, bool estTrParamsBFGS, bool centralDiff, gsl_vector* fixedParams) {
  if(estTrParamsBFGS) {
    if(fixLib) {
      this->bfgsAllInd = AllIndFixLib3TrParam2DegPolyHMM::create(this->depthsVec, this->sampleList, fixedParams, this->MAX_PLOIDY);
    }
    else {
      this->bfgsAllInd = AllInd3TrParam2DegPolyHMM::create(this->depthsVec, this->sampleList, this->MAX_PLOIDY);
    }
  }
  else {
    if(fixLib) {
      this->bfgsAllInd = AllIndFixLib0TrParam2DegPolyHMM::create(this->depthsVec, this->sampleList, fixedParams, this->MAX_PLOIDY);
    }
    else {
      this->bfgsAllInd = AllInd0TrParam2DegPolyHMM::create(this->depthsVec, this->sampleList, fixedParams, this->MAX_PLOIDY);
    }
  }
  this->bfgsAllInd->gradientDebug = this->gradientDebug;
  this->bfgsAllInd->setCentralDiffFlag(centralDiff);
}
AllInd3TrParam2DegPolyHMM* AllInd2Stages3TrParam2DegPolyHMM::getBWAllInd() {
  return this->bwAllInd;
}
AllInd3TrParam2DegPolyHMM* AllInd2Stages3TrParam2DegPolyHMM::getBFGSAllInd() {
  return this->bfgsAllInd;
}
void AllInd2Stages3TrParam2DegPolyHMM::print(FILE* stream) {
  // if read sconce from file, don't print bwAllInd, since params weren't initialized in bwAllInd
  if(this->bwAllInd != nullptr && !this->hasHMMBeenReadFromFile()) {
    fprintf(stream, "\nAllInd2Stages3TrParam2DegPolyHMM bwAllInd:\n");
    this->bwAllInd->print(stream);
  }
  if(this->bfgsAllInd != nullptr) {
    fprintf(stream, "\nAllInd2Stages3TrParam2DegPolyHMM bfgsAllInd:\n");
    this->bfgsAllInd->print(stream);
  }
  fprintf(stream, "\n\n################################################################################\n\n");
}
std::vector<std::string>* AllInd2Stages3TrParam2DegPolyHMM::getSampleList() {
  return this->sampleList;
}
void AllInd2Stages3TrParam2DegPolyHMM::setSampleList(std::vector<std::string>* sampleList) {
  this->sampleList = sampleList;
}

int AllInd2Stages3TrParam2DegPolyHMM::getInitGuessSize() const {
  return this->bfgsAllInd->getNumParamsToEst(); // Wed 01 Jul 2020 03:36:29 PM PDT every param to be estimated has an entry here
}

int AllInd2Stages3TrParam2DegPolyHMM::getBWInitGuessSize() const {
  if(this->bwAllInd == nullptr) {
    return 0;
  }
  return this->bwAllInd->getNumParamsToEst(); // Wed 01 Jul 2020 03:36:29 PM PDT every param to be estimated has an entry here
}

void AllInd2Stages3TrParam2DegPolyHMM::setBWAllMeanVarianceFn(gsl_vector* meanVarianceCoefVec) {
  if(this->bwAllInd != nullptr) {
    this->bwAllInd->setAllMeanVarianceFn(meanVarianceCoefVec);
  }
}

void AllInd2Stages3TrParam2DegPolyHMM::setBFGSAllMeanVarianceFn(gsl_vector* meanVarianceCoefVec) {
  if(this->bfgsAllInd != nullptr) {
    this->bfgsAllInd->setAllMeanVarianceFn(meanVarianceCoefVec);
  }
}

// likelihood methods
double AllInd2Stages3TrParam2DegPolyHMM::getBWAllIndLogLikelihood() {
  return this->bwAllInd->getLogLikelihood();
}
double AllInd2Stages3TrParam2DegPolyHMM::getBFGSAllIndLogLikelihood() {
  return this->bfgsAllInd->getLogLikelihood();
}

// bfgs methods. conversion between initGuess (in probability space) and initGuessAsParams (in BFGS space)
// is handled in subclasses
AllInd3TrParam2DegPolyHMM* AllInd2Stages3TrParam2DegPolyHMM::bfgs(gsl_vector* initGuess, int maxIters, bool verbose, bool debug) {
  this->bfgsAllInd = (AllInd3TrParam2DegPolyHMM*) this->bfgsAllInd->bfgs(initGuess, maxIters, verbose, debug);
  return this->bfgsAllInd;
}

/*
 * method to call bfgs numRuns times. uses initGuess if only doing one run. else, uses
 * baum welch results
 * returns the HMM with the best loglik
 */
AllInd3TrParam2DegPolyHMM* AllInd2Stages3TrParam2DegPolyHMM::callBFGSNTimes(gsl_vector* initGuess, int numRuns, int maxIters, bool verbose, bool debug) {
  std::vector<gsl_vector*>* baumWelchParamResults = this->bwAllInd->getBaumWelchParamResults();
  std::vector<gsl_vector*>* baumWelchParamResultsCopy = new std::vector<gsl_vector*>();
  for(unsigned int i = 0; i < baumWelchParamResults->size(); i++) {
    gsl_vector* currParams = (*baumWelchParamResults)[i];
    gsl_vector* currParamsCopy = gsl_vector_alloc(currParams->size);
    gsl_vector_memcpy(currParamsCopy, currParams);
    baumWelchParamResultsCopy->push_back(currParamsCopy);
  }
  this->bfgsAllInd->setBaumWelchParamResults(baumWelchParamResultsCopy);
  this->bfgsAllInd = (AllInd3TrParam2DegPolyHMM*) this->bfgsAllInd->callBFGSNTimes(numRuns, maxIters, verbose, debug);
  return this->bfgsAllInd;
}

void AllInd2Stages3TrParam2DegPolyHMM::setInitGuessNthTime(gsl_vector* initGuess, int iter) const {
  if(iter == 0) {
    this->setSharedInitGuess(initGuess, 1, 0.075, 12.5, 0.0225);
  }
  else if(iter == 1) {
    this->setSharedInitGuess(initGuess, 0.9, 0.025, 10, 0.13);
  }
  else if(iter == 2) {
    this->setSharedInitGuess(initGuess, 1.25, 0.12, 15, 0.02);
  }
}
/*
 * function to set elements of initGuess (which is the same size as paramsToEst) to shared
 * values of lib, beta, lambda, t
 */
void AllInd2Stages3TrParam2DegPolyHMM::setSharedInitGuess(gsl_vector* initGuess, double lib, double beta, double lambda, double t) const {
  int initGuessIdx = 0;
  // lib scaling factors
  for(int cellIdx = 0; cellIdx < this->bfgsAllInd->NUM_LIBS_TO_EST; cellIdx++) {
    gsl_vector_set(initGuess, initGuessIdx, lib);
    initGuessIdx++;
  }

  // shared transition parameters
  gsl_vector_set(initGuess, initGuessIdx, beta);
  initGuessIdx++;
  gsl_vector_set(initGuess, initGuessIdx, lambda);
  initGuessIdx++;

  // pairwise branch lengths
  for(int pairIdx = 0; pairIdx < this->NUM_CELLS; pairIdx++) {
    gsl_vector_set(initGuess, initGuessIdx + pairIdx, t);
  }
}
void AllInd2Stages3TrParam2DegPolyHMM::setBaumWelchInitGuess(gsl_vector* initGuess, int numBWIters, int numLibStarts, double libStartVal, bool verbose, bool debug) {
  if(this->bwAllInd != nullptr) {
    this->bwAllInd->setBaumWelchInitGuess(initGuess, numBWIters, numLibStarts, libStartVal, verbose, debug);
  }
}
void AllInd2Stages3TrParam2DegPolyHMM::copyBaumWelchIntoBFGSInitGuess(gsl_vector* bfgsInitGuess, gsl_vector* bwInitGuess) {
  int bfgsIdx = 0;
  int bwIdx = 0;
  // lib scaling factors
  // either lib scaling factors go into initGuess
  for(int cellIdx = 0; cellIdx < this->bfgsAllInd->NUM_LIBS_TO_EST; cellIdx++, bfgsIdx++, bwIdx++) {
    gsl_vector_set(bfgsInitGuess, bfgsIdx, gsl_vector_get(bwInitGuess, bwIdx));
  }
  // or they go into fixedParams
  double lib = 0;
  for(int cellIdx = 0; cellIdx < this->bfgsAllInd->NUM_FIXED_LIBS; cellIdx++, bwIdx++) {
    lib = gsl_vector_get(bwInitGuess, bwIdx);
    this->bfgsAllInd->setLibScalingFactor(cellIdx, lib);
  }

  // shared transition parameters
  for(; bwIdx < (int) bwInitGuess->size; bfgsIdx++, bwIdx++) {
    gsl_vector_set(bfgsInitGuess, bfgsIdx, gsl_vector_get(bwInitGuess, bwIdx));
  }
}

// methods to save results
void AllInd2Stages3TrParam2DegPolyHMM::viterbiDecode() {
  this->bfgsAllInd->viterbiDecodeAll();
}
void AllInd2Stages3TrParam2DegPolyHMM::saveViterbiDecodedCNA(std::string filename) {
  this->bfgsAllInd->saveAllViterbiDecodedCNA(filename);
}
void AllInd2Stages3TrParam2DegPolyHMM::saveCNAToBed(std::string filename) {
  this->bfgsAllInd->saveAllCNAToBed(filename);
}

void AllInd2Stages3TrParam2DegPolyHMM::saveParamEstimates(int cellIdx, std::string filename) {
  this->bfgsAllInd->saveParamEstimates(cellIdx, filename);
}

void AllInd2Stages3TrParam2DegPolyHMM::getParamsToEstFromFile(int cellIdx, std::string filename, int numExpectedLines, gsl_vector* meanVarianceCoefVec) {
  // first set up this->bfgsAllInd
  if(this->bfgsAllInd == nullptr) {
    this->bfgsAllInd = AllInd3TrParam2DegPolyHMM::create(this->depthsVec, this->sampleList, this->MAX_PLOIDY); // assumes alwasy resuming into this type of hmm, which is already done in IndThenPairs2Stages3TrParam2DegPolyHMM::optimIndCells
    this->bfgsAllInd->gradientDebug = this->gradientDebug;
  }
  // then read from file
  this->bfgsAllInd->getParamsToEstFromFile(cellIdx, filename, numExpectedLines, meanVarianceCoefVec);
}

/*
 * during getParamsToEstFromFile for getting sconce results from a file, bfgsAllInd is initialized. if it's nullptr, then haven't read from file
 */
bool AllInd2Stages3TrParam2DegPolyHMM::hasHMMBeenReadFromFile() {
  if(this->bfgsAllInd == nullptr) {
    return false;
  }
  return this->bfgsAllInd->hasIthHMMBeenReadFromFile(0); // always checking on 0'th, since bfgsAllInd only has 1 cell
}

gsl_vector* AllInd2Stages3TrParam2DegPolyHMM::getHMMParamsToEstReadFromFile() {
  if(this->bfgsAllInd == nullptr) {
    return nullptr;
  }
  return this->bfgsAllInd->getParamsToEstFromIthHMM(0); // always checking on 0'th, since bfgsAllInd only has 1 cell
}

