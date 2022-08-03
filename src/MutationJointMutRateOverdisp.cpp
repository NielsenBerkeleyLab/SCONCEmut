#include "MutationJointMutRateOverdisp.hpp"

MutationJointMutRateOverdisp::MutationJointMutRateOverdisp(std::vector<MutationList*>* mutListVec, std::vector<std::unordered_map<std::string, std::vector<int>*>*>* chrToViterbiPathMapVec, gsl_vector* branchLengths, bool verbose, bool gradientDebug) : MutationListContainer(mutListVec, 0, 2, verbose, gradientDebug) {
//MutationJointMutRateOverdisp::MutationJointMutRateOverdisp(std::vector<MutationList*>* mutListVec, std::vector<std::unordered_map<std::string, std::vector<int>*>*>* chrToViterbiPathMapVec, gsl_vector* branchLengths, bool verbose, bool gradientDebug) : MutationListContainer(mutListVec, 0, 1, verbose, gradientDebug) {
  this->branchLengths = branchLengths;
  this->setAllCoordSconceCNMap(chrToViterbiPathMapVec);

  // create MutationInd's out of MutationLists
  this->mutIndVec = new std::vector<MutationInd*>(this->mutListVec->size());
  std::vector<MutationList*>* currMutListVec = nullptr;
  for(unsigned int mutListIdx = 0; mutListIdx < this->mutListVec->size(); mutListIdx++) {
    currMutListVec = new std::vector<MutationList*>(1);
    (*currMutListVec)[0] = (*this->mutListVec)[mutListIdx];
    //(*this->mutIndVec)[mutListIdx] = new MutationInd(currMutListVec, this->gradientDebug, this->gradientDebug); // make these very quiet by default
    (*this->mutIndVec)[mutListIdx] = new MutationInd(currMutListVec, this->verbose, this->gradientDebug);
  }

  this->internalInitGuess = gsl_vector_alloc(1);
}

MutationJointMutRateOverdisp::~MutationJointMutRateOverdisp() {
  // TODO
}

double MutationJointMutRateOverdisp::getCnaToMutRateMu() {
  return gsl_vector_get(this->paramsToEst, this->MUTATION_VARS_START_IDX + 0);
}
double MutationJointMutRateOverdisp::getMutOverdispOmega() {
  return gsl_vector_get(this->paramsToEst, this->MUTATION_VARS_START_IDX + 1);
}

/*
 * function to return a vector of all the estimated mutation counts, collated across mutIndVec
 */
gsl_vector* MutationJointMutRateOverdisp::getAllIndEstMutCounts() {
  gsl_vector* mutCounts = gsl_vector_alloc(this->mutIndVec->size());
  for(unsigned int mutIndIdx = 0; mutIndIdx < this->mutIndVec->size(); mutIndIdx++) {
    gsl_vector_set(mutCounts, mutIndIdx, (*this->mutIndVec)[mutIndIdx]->getMutCountEst());
  }
  return mutCounts;
}

/*
 * optim l(mu, omega): 
 * omega ==> l(X1) = sum_i p(D_i|S_i=1) X1/n + p(D_i|S_i=0) (N-X1)/N, repeated for all cells
 * mu ==> log P(mu*t1*X1) + log P(mu*t2*X2) + ..., across all cells
 */
double MutationJointMutRateOverdisp::getLogLikelihood() {
  //std::cout << "MutationJointMutRateOverdisp::getLogLikelihood" << std::endl;
  //double currMu = gsl_vector_get(this->paramsToEst, this->MUTATION_VARS_START_IDX + 0);
  double currMu = this->getMutationVarEst(0);
  //double currMu = 0.01; // debugging, use a constant to keep mu and omega at initGuess values
  //double currOmega = gsl_vector_get(this->paramsToEst, this->MUTATION_VARS_START_IDX + 1);
  double currOmega = this->getMutationVarEst(1);
  //double currOmega = gsl_vector_get(this->paramsToEst, this->MUTATION_VARS_START_IDX + 0);
  //currOmega = 10;
  if(gsl_isnan(currMu) || gsl_isnan(currOmega)) {
    return GSL_NAN;
  }
  double currLl = 0;
  double totalLl = 0;
  double currT = 0;
  double currX = 0;
  //gsl_vector_set_all(this->internalInitGuess, 300);

  // for each MutationInd object, reest mut count X, for currOmega (optim omega)
  boost::asio::thread_pool pool(this->numThreads);
  for(unsigned int mutIndIdx = 0; mutIndIdx < this->mutIndVec->size(); mutIndIdx++) {
    //std::cout << "optim mutIndIdx: " << mutIndIdx << std::endl;
    // set omega
    (*this->mutIndVec)[mutIndIdx]->setAllMutOverdispOmega(currOmega);
    //(*this->mutIndVec)[mutIndIdx]->print(stdout);

    // optimize X's for this omega value
    // Submit a lambda object to the pool
    boost::asio::post(pool,
      [this, mutIndIdx]() { // see https://stackoverflow.com/a/7627218 for lambda function syntax explanation
        this->callIndvMutCountEst(mutIndIdx);
    });
  }
  // Wait for all tasks in the pool to complete.
  pool.join();

  // then sum up log(P(mu*t*x) across all cells (optim mu)
  for(unsigned int mutIndIdx = 0; mutIndIdx < this->mutIndVec->size(); mutIndIdx++) {
    if(gsl_isnan((*this->mutIndVec)[mutIndIdx]->initLl) || gsl_isinf((*this->mutIndVec)[mutIndIdx]->initLl)) {
      return GSL_NAN;
    }
    currX = (*this->mutIndVec)[mutIndIdx]->getMutCountEst(0);
    //std::cout << currX << std::endl;
    // calc sum log(P(mu*t*x) (similar to TwoCellFixLib0TrParam2DegPolyHMMWithMuts::getMutLogLikelihood, where mutation counts have a Poisson-like dist, with mean = cnaToMutRateMu * branchLength), but here there's only one branch
    currT = gsl_vector_get(this->branchLengths, mutIndIdx);
    currLl = currX * log(currMu * currT) - currMu * currT - lgamma(currX+1); // poisson dist, with mean = cnaToMutRateMu * currT
    totalLl += currLl;
  }
  return totalLl;
}
void MutationJointMutRateOverdisp::callIndvMutCountEst(int mutIndIdx) {
  if(this->gradientDebug) {
  //if(this->verbose) {
    double currMu = gsl_vector_get(this->paramsToEst, this->MUTATION_VARS_START_IDX + 0);
    double currOmega = gsl_vector_get(this->paramsToEst, this->MUTATION_VARS_START_IDX + 1);
    //double currMu = 0.01;
    //double currOmega = gsl_vector_get(this->paramsToEst, this->MUTATION_VARS_START_IDX + 0);
    std::lock_guard<std::mutex> lock(Optimizable::mtx); // from https://stackoverflow.com/a/18277334
    //std::cout << "Calling BFGS for MutationInd[" << mutIndIdx << "] (mutation count) in MutationJointMutRateOverdisp::getLogLikelihood(); currMu(" << currMu << "), currOmega(" << currOmega << ")" << std::endl;
    std::cout << "Calling SIMPLEX for MutationInd[" << mutIndIdx << "] (mutation count) in MutationJointMutRateOverdisp::getLogLikelihood(); currMu(" << currMu << "), currOmega(" << currOmega << ")" << std::endl;
    //printRowVector(this->internalInitGuess);
    (*this->mutIndVec)[mutIndIdx]->print(stdout);
  }
  //(*this->mutIndVec)[mutIndIdx]->bfgs(this->internalInitGuess, 500, this->verbose, this->gradientDebug);
  //(*this->mutIndVec)[mutIndIdx]->bfgs(this->internalInitGuess, 500, this->gradientDebug, this->gradientDebug);
  //(*this->mutIndVec)[mutIndIdx]->estMutCountsPerBranch(nullptr, this->internalInitGuess, 500, this->gradientDebug); // mutation bfgs short circuit here
  gsl_vector* initGuess = nullptr;
  //if(compareDoubles(0.0, (*this->mutIndVec)[mutIndIdx]->getMutCountEst(0))) { // if already estimated mutation count, don't estimate again (saves time but might make optimization get stuck)/ reuseMutEsts
    initGuess = gsl_vector_alloc((*this->mutIndVec)[mutIndIdx]->NUM_MUTATION_COUNTS_TO_EST);
    //double validInitGuessMax = (*this->mutIndVec)[mutIndIdx]->getValidOptimParamMax();
    //double initGuessValToUse = validInitGuessMax * 0.5; // set to half so don't get stuck in upper or lower boundary
    //double initGuessValToUse = (*(*this->mutIndVec)[mutIndIdx]->mutListVec)[0]->coordVec->size() * 0.5; // set to half so don't get stuck in upper or lower boundary
    //double initGuessValToUse = (*(*this->mutIndVec)[mutIndIdx]->mutListVec)[0]->coordVec->size() * 0.8; // set to half so don't get stuck in upper or lower boundary
    double initGuessValToUse = (*(*this->mutIndVec)[mutIndIdx]->mutListVec)[0]->coordVec->size() * 0.1; // set to half so don't get stuck in upper or lower boundary
    //double initGuessValToUse = (*(*this->mutIndVec)[mutIndIdx]->mutListVec)[0]->coordVec->size() * 0.01; // set to half so don't get stuck in upper or lower boundary
    gsl_vector_set_all(initGuess, initGuessValToUse);
  //}
  (*this->mutIndVec)[mutIndIdx]->estMutCountsPerBranch(nullptr, initGuess, 500, this->gradientDebug); // mutation bfgs short circuit here
}

void MutationJointMutRateOverdisp::print(FILE* stream) {
  for(unsigned int mutIndIdx = 0; mutIndIdx < this->mutIndVec->size(); mutIndIdx++) {
    fprintf(stream, "MutationInd[%i] (mut count X_[%i])\n", mutIndIdx, mutIndIdx);
    (*this->mutIndVec)[mutIndIdx]->print(stream);
    fprintf(stream, "\n");
  }

  fprintf(stream, "MutationJointMutRateOverdisp branchLengths:\n");
  printColVector(stream, this->branchLengths);

  fprintf(stream, "MutationJointMutRateOverdisp paramsToEst ([mu, omega]):\n");
  printColVector(stream, this->paramsToEst);
}
Optimizable* MutationJointMutRateOverdisp::bfgs(gsl_vector* initGuess, int maxIters, bool verbose, bool debug) { // this one should convert initGuess into BFGS space, and call Optimizable::bfgs
  MutationListContainer* bestGuessMutPair = this;
  gsl_vector* initGuessAsParams = gsl_vector_alloc(initGuess->size);
  bestGuessMutPair->convertProbToParam(initGuessAsParams, initGuess);
  Optimizable::bfgs(initGuessAsParams, bestGuessMutPair, maxIters, verbose, debug);
  //Optimizable::simplex(initGuessAsParams, bestGuessMutPair, maxIters, verbose, debug);
  gsl_vector_free(initGuessAsParams);
  return bestGuessMutPair;
}
// no upper limit on param values mu or omega
double MutationJointMutRateOverdisp::getValidOptimParamMax() const {
  return std::numeric_limits<double>::max();
}

// If all sconceMutParams files exist, then can skip over mu/omega estimation (returns true). If missing some, must redo estimation (returns false)
double MutationJointMutRateOverdisp::getMutParamsFromFile(std::vector<std::string>* sampleList, int maxPloidy, std::string sconceEstimatesPath) {
  std::lock_guard<std::mutex> lock(Optimizable::mtx); // from https://stackoverflow.com/a/18277334
  std::cout << "######## ATTEMPTING TO READ INDV CELL .sconceMutParams FILES WITH PATH " << sconceEstimatesPath << " ######## " << std::endl;
  double status = 0;
  double runningCnaToMutRateMu = 0;
  double runningMutOverdispOmega = 0;
  MutationInd* mutInd = nullptr;
  for(unsigned int mutIndIdx = 0; mutIndIdx < this->mutIndVec->size(); mutIndIdx++) {
    std::string hmmName = (*sampleList)[mutIndIdx];
    boost::replace_all(hmmName, ",", "__");
    std::string currFile = sconceEstimatesPath + "__" + hmmName + "__k" + std::to_string(maxPloidy) + ".sconceMutParams"; // based on IndThenPairs2Stages3TrParam2DegPolyHMMWithMuts::estimateMutParamsFromIndCells() naming conventions
    if(!boost::filesystem::exists(currFile)) {
      std::cerr << currFile << " does not exist, will not read .sconceMutParams from file." << std::endl;
      return GSL_NAN;
    }
    mutInd = (*this->mutIndVec)[mutIndIdx];
    status = mutInd->getMutParamsFromFile(currFile, 3); // mu, omega, num mutations
    if(gsl_isnan(status)) {
      return GSL_NAN;
    }
    if(mutIndIdx == 0) {
      runningCnaToMutRateMu = mutInd->getCnaToMutRateMu();
      runningMutOverdispOmega = mutInd->getMutOverdispOmega();
    }
    else {
      // if any of the estimates that should be the same are actually different, return NAN
      if(!compareDoubles(mutInd->getCnaToMutRateMu(), runningCnaToMutRateMu) || !compareDoubles(mutInd->getMutOverdispOmega(), runningMutOverdispOmega)) {
        return GSL_NAN;
      }
    }
  }

  // save
  gsl_vector_set(this->paramsToEst, this->MUTATION_VARS_START_IDX + 0, runningCnaToMutRateMu); // mu
  gsl_vector_set(this->paramsToEst, this->MUTATION_VARS_START_IDX + 1, runningMutOverdispOmega); // omega
  return 0;
}

