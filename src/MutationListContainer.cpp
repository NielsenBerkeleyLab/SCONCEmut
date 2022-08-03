#include "MutationListContainer.hpp"

MutationListContainer::MutationListContainer(std::vector<MutationList*>* mutListVec, int numMutCountsToEst, int numOverdispVarsToEst, bool verbose, bool gradientDebug) : NUM_MUTATION_COUNTS_TO_EST(numMutCountsToEst), NUM_MUTATION_VARS_TO_EST(numOverdispVarsToEst), MUTATION_COUNT_START_IDX(0), MUTATION_VARS_START_IDX(numMutCountsToEst) {
  this->mutListVec = mutListVec;
  this->internalInitGuess = nullptr;
  this->verbose = verbose;
  this->gradientDebug = gradientDebug;

  // Optimizable params
  this->paramsToEst = gsl_vector_alloc(this->NUM_MUTATION_COUNTS_TO_EST + NUM_MUTATION_VARS_TO_EST); // X1,X2,X3 (counts of muts on each branch) and/or mu+omega
  gsl_vector_set_zero(this->paramsToEst);
  this->fixedParams = nullptr;//gsl_vector_alloc(1); // seqErrorRate
  this->optimSuccess = false;
  this->maxNumBFGSStarts = 1; // maximum number of starting points for optimization, usually set by setInitGuessNthTime
  this->changeInBFGSLoglikelihood = 0; // total change in loglikelihood after BFGS
  this->initGuessCopyBeforeBFGS = nullptr; // copy of initGuess, before BFGS. useful for keeping variables from changing too much from initial values
  this->bestOptimLLInitGuess = nullptr; // initGuess with best optimized loglikelihood, found via calling BFGS N times
  this->BFGSParamResults = nullptr; // vector of paramsToEst from BFGS, starting from different parameter sets, stored in order discovered

  this->probParamConversionVec = gsl_vector_alloc(this->getNumParamsToEst()); // intermediate for evalLikelihoodAtPoint. Hold param conversion somewhere instead of constantly allocating new memory

  this->simParamsToEst = nullptr; // copies of paramsToEst and fixedParams, useful for comparing to simualation params if we have them
  this->simFixedParams = nullptr;
}

MutationListContainer::~MutationListContainer() {
  if(this->mutListVec != nullptr) {
    for(unsigned int mutListIdx = 0; mutListIdx < this->mutListVec->size(); mutListIdx++) {
      delete (*this->mutListVec)[mutListIdx];
    }
    delete this->mutListVec;
  }
  gsl_vector_free(this->internalInitGuess);
}

// chrToViterbiPathMapVec is a cellIdx:[chr:vector of viterbi decoded paths] map
MutationListContainer* MutationListContainer::estMutCountsPerBranch(std::vector<std::unordered_map<std::string, std::vector<int>*>*>* chrToViterbiPathMapVec, gsl_vector* initGuessParam, int maxIters, bool verbose) {
  //std::cout << "MutationListContainer::estMutCountsPerBranch" << std::endl;
  // update MutationLists with these viterbi decodings
  if(chrToViterbiPathMapVec != nullptr) {
    this->setAllCoordSconceCNMap(chrToViterbiPathMapVec);
  }

  // if given a new initGuess, use it. Otherwise, use whatever is currently the best guess
  gsl_vector* initGuess = nullptr;
  if(initGuessParam == nullptr) {
    initGuess = gsl_vector_alloc(this->getNumParamsToEst());
    gsl_vector_memcpy(initGuess, this->paramsToEst);
    //gsl_vector_set_all(initGuess, 150);
  }
  else {
    initGuess = initGuessParam;
  }
  // save somewhere to reset to if bfgs fails
  gsl_vector* initGuessCopy = gsl_vector_alloc(initGuess->size);
  gsl_vector_memcpy(initGuessCopy, initGuess);
  if(this->gradientDebug) {
    std::lock_guard<std::mutex> lock(Optimizable::mtx); // from https://stackoverflow.com/a/18277334
    printRowVector(initGuess);
  }

  // call bfgs
  this->bfgs(initGuess, maxIters, verbose, this->gradientDebug);
  // if no improvement in ll && was given a new initGuess, then new initGuess failed
  if(compareDoubles(0.0, this->changeInBFGSLoglikelihood) && initGuessParam != nullptr) {
  //if(!this->optimSuccess) {
    //printRowVector(initGuess);
    //gsl_vector_set_all(initGuess, 100);
    double maxMutCount = 0;
    double currMutCount = 0;
    for(int i = 0; i < this->NUM_MUTATION_COUNTS_TO_EST; i++) {
      currMutCount = gsl_vector_get(this->paramsToEst, this->MUTATION_COUNT_START_IDX + i);
      if(maxMutCount < currMutCount) {
        maxMutCount = currMutCount;
      }
    }
    gsl_vector_set_all(initGuess, maxMutCount * 2.0 / 3.0);
    //std::cout << "WARNING: MutationListContainer::estMutCountsPerBranch failed, retrying with initGuess:" << std::endl;
    //printRowVector(initGuess);
    if(verbose) {
      std::lock_guard<std::mutex> lock(Optimizable::mtx); // from https://stackoverflow.com/a/18277334
      std::cout << "WARNING: MutationListContainer::estMutCountsPerBranch failed, retrying with initGuess:" << std::endl;
      printRowVector(initGuess);
    }
    this->bfgs(initGuess, maxIters, verbose, this->gradientDebug);
  //gsl_vector_memcpy(initGuess, this->initGuessCopyBeforeBFGS);
    //std::cout << this->changeInBFGSLoglikelihood << std::endl;
  }
  //if(!this->optimSuccess) {
  // if no improvement in ll && was given a new initGuess, then new initGuess failed
  if(compareDoubles(0.0, this->changeInBFGSLoglikelihood) && initGuessParam != nullptr) {
    //gsl_vector_set_all(initGuess, 50);
    gsl_vector_set_all(initGuess, 3);
    //std::cout << "WARNING: MutationListContainer::estMutCountsPerBranch failed, retrying with initGuess:" << std::endl;
    //printRowVector(initGuess);
    if(verbose) {
      std::lock_guard<std::mutex> lock(Optimizable::mtx); // from https://stackoverflow.com/a/18277334
      std::cout << "WARNING: MutationListContainer::estMutCountsPerBranch failed, retrying with initGuess:" << std::endl;
      printRowVector(initGuess);
    }
    this->bfgs(initGuess, maxIters, verbose, this->gradientDebug);
    //std::cout << this->changeInBFGSLoglikelihood << std::endl;
  }
  // if no improvement in ll && was given a new initGuess, then new initGuess (and all retries) failed
  if(compareDoubles(0.0, this->changeInBFGSLoglikelihood) && initGuessParam != nullptr) {
    gsl_vector_memcpy(initGuess, initGuessCopy);
    gsl_vector_memcpy(this->paramsToEst, initGuessCopy);
    if(verbose) {
      std::lock_guard<std::mutex> lock(Optimizable::mtx); // from https://stackoverflow.com/a/18277334
      std::cout << "Reseting to: ";
      printRowVector(initGuess);
    }
    //gsl_vector_memcpy(initGuess, this->initGuessCopyBeforeBFGS);
  }
  gsl_vector_free(initGuessCopy);
  if(initGuess != initGuessParam) {
    gsl_vector_free(initGuess);
  }
  return this;
}

// chrToViterbiPathMapVec is a cellIdx:[chr:vector of viterbi decoded paths] map
void MutationListContainer::setAllCoordSconceCNMap(std::vector<std::unordered_map<std::string, std::vector<int>*>*>* chrToViterbiPathMapVec) {
  // update MutationLists with these viterbi decodings
  for(unsigned int mutListIdx = 0; mutListIdx < this->mutListVec->size(); mutListIdx++) {
    //std::cout << "MutationListContainer::setAllCoordSconceCNMap mutListIdx: " << mutListIdx << std::endl;
    (*this->mutListVec)[mutListIdx]->setCoordSconceCNMap((*chrToViterbiPathMapVec)[mutListIdx]);
    //(*this->mutListVec)[mutListIdx]->print(stdout);
  }
}

double MutationListContainer::getMutCountEst(int branchIdx) const {
  return gsl_vector_get(this->paramsToEst, this->MUTATION_COUNT_START_IDX + branchIdx);
}
double MutationListContainer::getMutationVarEst(int varIdx) {
  return gsl_vector_get(this->paramsToEst, this->MUTATION_VARS_START_IDX + varIdx);
}
void MutationListContainer::setAllMutOverdispOmega(double omega) {
  for(unsigned int mutListIdx = 0; mutListIdx < this->mutListVec->size(); mutListIdx++) {
    (*this->mutListVec)[mutListIdx]->setMutOverdispOmega(omega);
  }
}

/*
 * function to save mutation params for the given mutListIdx to the given file
 *
 * format is
 * MutationList[mutListIdx]::mutOverdispOmega
 * fixedParams from calling object
 * paramsToEst from calling object
 */
void MutationListContainer::saveMutParamsToFile(int mutListIdx, std::string filename) {
  this->saveMutParamsToFile(mutListIdx, filename, "w");
}
void MutationListContainer::saveMutParamsToFile(int mutListIdx, std::string filename, const char* mode) {
  gsl_vector* writingBuffer = nullptr;
  if(this->fixedParams != nullptr) {
    writingBuffer = gsl_vector_alloc(1 + this->fixedParams->size + this->paramsToEst->size); // MutationList::mutOverdispOmega + fixedParams + paramsToEst
  }
  else {
    writingBuffer = gsl_vector_alloc(1 + this->paramsToEst->size); // MutationList::mutOverdispOmega + paramsToEst
  }
  int bufferIdx = 0;

  // store MutationList::mutOverdispOmega
  gsl_vector_set(writingBuffer, bufferIdx, (*this->mutListVec)[mutListIdx]->getMutOverdispOmega());
  bufferIdx++;

  // store fixedParams
  if(this->fixedParams != nullptr) {
    for(unsigned int paramIdx = 0; paramIdx < this->fixedParams->size; paramIdx++, bufferIdx++) {
      gsl_vector_set(writingBuffer, bufferIdx, gsl_vector_get(this->fixedParams, paramIdx));
    }
  }

  // store paramsToEst
  for(unsigned int paramIdx = 0; paramIdx < this->paramsToEst->size; paramIdx++, bufferIdx++) {
    gsl_vector_set(writingBuffer, bufferIdx, gsl_vector_get(this->paramsToEst, paramIdx));
  }

  FILE* outFile = fopen(filename.c_str(), mode);
  gsl_block_fprintf(outFile, writingBuffer->block, "%.40f");
  fclose(outFile);

  gsl_vector_free(writingBuffer);
}

double MutationListContainer::getMutParamsFromFile(std::string filename, int numExpectedLinesPerFile) {
  // construct paramsToEst from the readingBuffer, then call set methods
  int numParamsToRead = this->paramsToEst->size;
  gsl_vector* readingBuffer = gsl_vector_alloc(numParamsToRead); // x's, mu/omega
  gsl_vector* readParamsToEst = gsl_vector_alloc(this->paramsToEst->size);
  int bufferIdx = 0;

  // read in file
  FILE* currFile = fopen(filename.c_str(), "r");
  gsl_block_fscanf(currFile, readingBuffer->block);
  fclose(currFile);

  // store mutation counts
  for(int mutIdx = 0; mutIdx < this->NUM_MUTATION_COUNTS_TO_EST; mutIdx++, bufferIdx++) {
    gsl_vector_set(readParamsToEst, this->MUTATION_COUNT_START_IDX + mutIdx, gsl_vector_get(readingBuffer, bufferIdx));
  }
  // store mutation parameters (mu/omega)
  for(int mutIdx = 0; mutIdx < this->NUM_MUTATION_VARS_TO_EST; mutIdx++, bufferIdx++) {
    gsl_vector_set(readParamsToEst, this->MUTATION_VARS_START_IDX + mutIdx, gsl_vector_get(readingBuffer, bufferIdx));
  }

  if(bufferIdx != numExpectedLinesPerFile) {
    std::lock_guard<std::mutex> lock(Optimizable::mtx); // from https://stackoverflow.com/a/18277334
    std::cerr << "WARNING: expected " << numExpectedLinesPerFile << " entries from " << filename << " but only saved " << bufferIdx << " entries. Ignoring " << filename << "." << std::endl;
    return GSL_NAN;
  }

  // save and update all variables
  this->setParamsToEst(readParamsToEst);

  gsl_vector_free(readingBuffer);
  gsl_vector_free(readParamsToEst);

  return 0;
}

double MutationListContainer::setParamsToEst(gsl_vector* params) {
  gsl_vector_memcpy(this->paramsToEst, params);
  return 0;
}

/*
 * probs = constrained space (probabilities in other classes, nonnegative here for Poisson-ish dist), X[1,X2,X3] or overdisp var
 * params = unconstrained bfgs space, y1,y2,y3
 * x1 = exp(y1) > 0
 * y1 = ln(x1) \in [-inf, inf]
 */
void MutationListContainer::convertProbToParam(gsl_vector* dest, const gsl_vector* src) const {
  gsl_vector_memcpy(dest, src);
  return;
  double x = 0;
  for(int mutIdx = 0; mutIdx < this->NUM_MUTATION_COUNTS_TO_EST; mutIdx++) {
    x = gsl_vector_get(src, this->MUTATION_COUNT_START_IDX + mutIdx);
    gsl_vector_set(dest, this->MUTATION_COUNT_START_IDX + mutIdx, log(x));
  }
  for(int varIdx = 0; varIdx < this->NUM_MUTATION_VARS_TO_EST; varIdx++) {
    x = gsl_vector_get(src, this->MUTATION_VARS_START_IDX + varIdx);
    gsl_vector_set(dest, this->MUTATION_VARS_START_IDX + varIdx, log(x));
  }
}
void MutationListContainer::convertParamToProb(gsl_vector* dest, const gsl_vector* src) const {
  gsl_vector_memcpy(dest, src);
  return;
  double y = 0;
  for(int mutIdx = 0; mutIdx < this->NUM_MUTATION_COUNTS_TO_EST; mutIdx++) {
    y = gsl_vector_get(src, this->MUTATION_COUNT_START_IDX + mutIdx);
    gsl_vector_set(dest, this->MUTATION_COUNT_START_IDX + mutIdx, exp(y));
  }
  for(int varIdx = 0; varIdx < this->NUM_MUTATION_VARS_TO_EST; varIdx++) {
    y = gsl_vector_get(src, this->MUTATION_VARS_START_IDX + varIdx);
    gsl_vector_set(dest, this->MUTATION_VARS_START_IDX + varIdx, exp(y));
  }
}
void MutationListContainer::setInitGuessNthTime(gsl_vector* initGuess, int iter, int numTotalRuns) const {
  switch(iter) {
    case 0:
      gsl_vector_set_all(initGuess, 10);
      break;
  }
}
Optimizable* MutationListContainer::bfgs(gsl_vector* initGuess, int maxIters, bool verbose, bool debug) { // this one should convert initGuess into BFGS space, and call Optimizable::bfgs
  MutationListContainer* bestGuessMutPair = this;
  gsl_vector* initGuessAsParams = gsl_vector_alloc(initGuess->size);
  bestGuessMutPair->convertProbToParam(initGuessAsParams, initGuess);
  //Optimizable::bfgs(initGuessAsParams, bestGuessMutPair, maxIters, verbose, debug);
  Optimizable::simplex(initGuessAsParams, bestGuessMutPair, maxIters, verbose, debug);
  gsl_vector_free(initGuessAsParams);
  return bestGuessMutPair;
}

void MutationListContainer::print(FILE* stream) {
  if(this->mutListVec != nullptr) {
    for(unsigned int mutListIdx = 0; mutListIdx < this->mutListVec->size(); mutListIdx++) {
      fprintf(stream, "MutationList[%i]:\n", mutListIdx);
      (*this->mutListVec)[mutListIdx]->print(stream);
    }
  }

  fprintf(stream, "MutationListContainer paramsToEst:\n");
  printColVector(stream, this->paramsToEst);

  fprintf(stream, "MutationListContainer loglikelihood: %.10f\n", this->getLogLikelihood());
}

double MutationListContainer::checkOptimProbValidity(gsl_vector* probs) const {
  // check if any prob (mutation count) is negative
  double curr = 0;
  double validMax = this->getValidOptimParamMax();
  for(unsigned int i = 0; i < probs->size; i++) {
    curr = gsl_vector_get(probs, i);
    //if(curr < 1e-2) {
    if(curr < 0) {
      //std::cout << "NAN: curr " << curr << " < 1e-2" << std::endl;
      return GSL_NAN;
    }
    if(curr > validMax) {
      //std::cout << "NAN: curr " << curr << " > " << validMax << std::endl;
      return GSL_NAN;
    }
    if(gsl_isnan(curr)) {
      //std::cout << "NAN: curr is nan" << std::endl;
      return GSL_NAN;
    }
  }
  double paramSum = gsl_blas_dasum(probs);
  if(paramSum > validMax) {
    //std::cout << "NAN: paramSum " << paramSum << " > " << validMax << std::endl;
    return GSL_NAN;
  }
  return 0;
}

double MutationListContainer::checkStateValidity(double epsilon) const {
  return 0;
}
void MutationListContainer::setSimParamsToEst(gsl_vector* params) {
  return;
}
void MutationListContainer::setSimFixedParams(gsl_vector* params) {
  return;
}
void MutationListContainer::miscFunctions() {
  return;
}

