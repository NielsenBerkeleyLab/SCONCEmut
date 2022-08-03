#include "IndThenPairs2Stages3TrParam2DegPolyHMMWithMuts.hpp"


IndThenPairs2Stages3TrParam2DegPolyHMMWithMuts::IndThenPairs2Stages3TrParam2DegPolyHMMWithMuts(std::vector<DepthPair*>* depths, std::vector<std::string>* sampleList, std::vector<MutationList*>* mutListVec, int maxPloidy, bool runIndBFGS, int numPairsStage2, bool verbose, bool gradientDebug, bool centralDiff) : IndThenPairs2Stages3TrParam2DegPolyHMM(depths, sampleList, maxPloidy, runIndBFGS, numPairsStage2, verbose, gradientDebug, centralDiff) {
  this->mutListVec = mutListVec;
  this->mutIndVec = nullptr;
  this->mutJointOverdisp = nullptr;
  this->mutJointMutRateOverdisp = nullptr;
  this->cnaToMutRateMu = 0;
  this->mutOverdispOmega = 0;
  this->allIndEstMutCounts = nullptr;
  this->allPairedEstMutCounts = nullptr;
}

IndThenPairs2Stages3TrParam2DegPolyHMMWithMuts::~IndThenPairs2Stages3TrParam2DegPolyHMMWithMuts() {
  // TODO
}
//void IndThenPairs2Stages3TrParam2DegPolyHMMWithMuts::cleanUpIndOptim(); // override to also delete mutIndVec and mutJointDisp? need to keep mutListVec TODO

void IndThenPairs2Stages3TrParam2DegPolyHMMWithMuts::print(FILE* stream) {
  if(this->ind2StagesVec != nullptr) {
    fprintf(stream, "\nIndThenPairs2Stages3TrParam2DegPolyHMM ind2StagesVec:\n");

    AllInd2Stages3TrParam2DegPolyHMM* currInd2Stage = nullptr;
    for(unsigned int i = 0; i < this->ind2StagesVec->size(); i++) {
      currInd2Stage = (*this->ind2StagesVec)[i];
      currInd2Stage->print(stream);
    }
  }

  if(this->mutJointMutRateOverdisp != nullptr) {
    this->mutJointMutRateOverdisp->print(stream);
  }
  else if(this->mutJointOverdisp != nullptr) {
    this->mutJointOverdisp->print(stream);
  }
  else if(this->mutIndVec != nullptr && this->mutIndVec->size() > 0) {
    for(unsigned int mutIndIdx = 0; mutIndIdx < this->mutIndVec->size(); mutIndIdx++) {
      fprintf(stream, "MutationInd[%i] (mut count X_[%i])\n", mutIndIdx, mutIndIdx);
      (*this->mutIndVec)[mutIndIdx]->print(stream);
      fprintf(stream, "\n");
    }
  }
  else {
    for(unsigned int mutListIdx = 0; mutListIdx < this->mutListVec->size(); mutListIdx++) {
      fprintf(stream, "MutationList[%i]\n", mutListIdx);
      (*this->mutListVec)[mutListIdx]->print(stream);
      fprintf(stream, "\n");
    }
  }

  if(this->stage2Pairs != nullptr) {
    fprintf(stream, "\nIndThenPairs2Stages3TrParam2DegPolyHMMWithMuts stage2Pairs:\n");
    this->stage2Pairs->print(stream);
  }

  if(this->allIndEstMutCounts != nullptr) {
    fprintf(stream, "allIndEstMutCounts:\n");
    printColVector(stream, this->allIndEstMutCounts);
  }
  if(this->allPairedEstMutCounts != nullptr) {
    fprintf(stream, "allPairedEstMutCounts:\n");
    printColVector(stream, this->allPairedEstMutCounts);
  }

  fprintf(stream, "IndThenPairs2Stages3TrParam2DegPolyHMMWithMuts cnaToMutRateMu, mutOverdispOmega:\n%.40f\n%.40f\n", this->cnaToMutRateMu, this->mutOverdispOmega);
}
double IndThenPairs2Stages3TrParam2DegPolyHMMWithMuts::getCnaToMutRateMu() {
  return this->cnaToMutRateMu;
}
double IndThenPairs2Stages3TrParam2DegPolyHMMWithMuts::getMutOverdispOmega() {
  return this->mutOverdispOmega;
}
gsl_vector* IndThenPairs2Stages3TrParam2DegPolyHMMWithMuts::getAllIndEstMutCounts() {
  if(this->allIndEstMutCounts != nullptr) {
    return this->allIndEstMutCounts;
  }
  else if(this->mutJointMutRateOverdisp != nullptr) {
    this->allIndEstMutCounts = this->mutJointMutRateOverdisp->getAllIndEstMutCounts();
    return this->allIndEstMutCounts;
  }
  return nullptr;
}
gsl_vector* IndThenPairs2Stages3TrParam2DegPolyHMMWithMuts::getAllPairedEstMutCounts() {
  if(this->allPairedEstMutCounts != nullptr) {
    return this->allPairedEstMutCounts;
  }
  else if(this->mutJointMutRateOverdisp != nullptr) {
    this->allPairedEstMutCounts = this->stage2Pairs->getAllPairedEstMutCounts();
    return this->allPairedEstMutCounts;
  }
  return nullptr;
}

//void IndThenPairs2Stages3TrParam2DegPolyHMMWithMuts::setUpInd2Stages() {
//  // don't think actually need to do this. definitely need to add a step to est mu, but I think that can be local to this class and don't think that has to be OneCell*WithMuts
//}

/*
 * function to estimate mutation count X (and rate mu, given t) for each individual cell, then get the median of all mu estimates and store in this->cnaToMutRateMu
 */
void IndThenPairs2Stages3TrParam2DegPolyHMMWithMuts::estCnaToMutRateMuIndependentlyForAllIndCells() {
  // create vector of MutationInd
  this->mutIndVec = new std::vector<MutationInd*>(this->mutListVec->size());
  std::vector<MutationList*>* currMutListVec = nullptr;
  for(unsigned int mutListIdx = 0; mutListIdx < this->mutListVec->size(); mutListIdx++) {
    currMutListVec = new std::vector<MutationList*>(1);
    (*currMutListVec)[0] = (*this->mutListVec)[mutListIdx];
    (*mutIndVec)[mutListIdx] = new MutationInd(currMutListVec, this->verbose, this->gradientDebug);
  }

  // optimize mutation count for each MutationInd. parallelize?
  HMM* hmm = nullptr;
  std::vector<std::unordered_map<std::string, std::vector<int>*>*>* chrToViterbiPathMapVec = nullptr;
  gsl_vector* initGuess = gsl_vector_alloc(1);
  for(unsigned int mutIndIdx = 0; mutIndIdx < this->mutIndVec->size(); mutIndIdx++) {
    gsl_vector_set(initGuess, 0, 10);
    if(this->runIndBFGS) {
      hmm = (*(*this->ind2StagesVec)[mutIndIdx]->getBFGSAllInd()->getHMMs())[0]; // assumes each AllInd obj only holds 1 cell (in the 0th position of hmmVec)
    } else {
      hmm = (*(*this->ind2StagesVec)[mutIndIdx]->getBWAllInd()->getHMMs())[0];
    }
    chrToViterbiPathMapVec = hmm->getChrToViterbiPathMapVec();
    (*this->mutIndVec)[mutIndIdx]->estMutCountsPerBranch(chrToViterbiPathMapVec, initGuess, 500, this->verbose);
  }
  gsl_vector_free(initGuess);

  // get t from each AllInd obj
  gsl_vector* indBranchEsts = this->getIndOptimAllBranchEsts();

  // calculate mu = X/t for each cell
  gsl_vector* medianMuEsts = gsl_vector_alloc(this->mutIndVec->size());
  double currX = 0;
  double currT = 0;
  double currMu = 0;
  for(unsigned int mutIndIdx = 0; mutIndIdx < this->mutIndVec->size(); mutIndIdx++) {
    currX = (*this->mutIndVec)[mutIndIdx]->getMutCountEst();
    currT = gsl_vector_get(indBranchEsts, mutIndIdx);
    currMu = currX / currT;
    gsl_vector_set(medianMuEsts, mutIndIdx, currMu);
  }

  // get median of mu estimates
  gsl_sort(medianMuEsts->data, 1, medianMuEsts->size); // https://www.gnu.org/software/gsl/doc/html/sort.html#c.gsl_sort

   // store in this->cnaToMutRateMu
  this->cnaToMutRateMu = gsl_stats_median_from_sorted_data(medianMuEsts->data, 1, medianMuEsts->size); // https://www.gnu.org/software/gsl/doc/html/statistics.html#c.gsl_stats_median_from_sorted_data
  for(unsigned int mutIndIdx = 0; mutIndIdx < this->mutIndVec->size(); mutIndIdx++) {
    (*this->mutIndVec)[mutIndIdx]->setCnaToMutRateMu(this->cnaToMutRateMu); // store summarized mu in each mutInd object so it gets saved to file later
  }
}

/*
 * function to estimate the overdispersion parameter for the beta binomial distribution jointly across all individual cells, then store in this->mutOverdispOmega
 */
void IndThenPairs2Stages3TrParam2DegPolyHMMWithMuts::estMutOverdispOmegaJointlyAcrossAllIndCells() {
  // create an MutationJointOverdisp (Optimizable) object that holds a list of MutationInd objects
  this->mutJointOverdisp = new MutationJointOverdisp(this->mutIndVec, this->verbose, this->gradientDebug);

  // optimize omega
  gsl_vector* initGuess = gsl_vector_alloc(1);
  gsl_vector_set_all(initGuess, 10);
  this->mutJointOverdisp->bfgs(initGuess, 500, this->verbose, this->gradientDebug);
  gsl_vector_free(initGuess);

  // store omega in calling object and in MutationLists directly
  this->mutOverdispOmega = this->mutJointOverdisp->getMutationVarEst(0);
  for(unsigned int mutListIdx = 0; mutListIdx < this->mutListVec->size(); mutListIdx++) {
    (*this->mutListVec)[mutListIdx]->setMutOverdispOmega(this->mutOverdispOmega);
  }
}

/*
 * function to jointly estimate cnaToMutRateMu and mutOverdispOmega jointly across all individual cells.
 * That is, estimate shared global estimates for mu and omega.
 */
void IndThenPairs2Stages3TrParam2DegPolyHMMWithMuts::estCnaToMutRateMuAndMutOverdispOmegaJointlyAcrossAllIndCells(std::string sconceEstimatesPath) {
  gsl_vector* indBranchEsts = this->getIndOptimAllBranchEsts();
  //printColVector(indBranchEsts);
  std::vector<std::unordered_map<std::string, std::vector<int>*>*>* allChrToViterbiPathMapVec = this->getIndOptimAllChrToViterbiPathMapVec();

  this->mutJointMutRateOverdisp = new MutationJointMutRateOverdisp(this->mutListVec, allChrToViterbiPathMapVec, indBranchEsts, this->verbose, this->gradientDebug);
  this->mutJointMutRateOverdisp->setNumThreads(this->numThreads);
  //this->mutJointMutRateOverdisp->print(stdout);

  // check if sconceMutParams files have already been created. If all of them exist, then can skip over mu/omega estimation. If missing some, must redo estimation
  double status = this->mutJointMutRateOverdisp->getMutParamsFromFile(this->sampleList, this->MAX_PLOIDY, sconceEstimatesPath); // TODO put this back in to read sconceMutParams
  //double status = GSL_NAN;
  //std::cout << "getMutParamsFromFile status: " << status << std::endl;

  // if didn't successfully read from files, do joint estimation
  if(gsl_isnan(status)) {
    // optimize mu and omega jointly
    gsl_vector* initGuess = gsl_vector_alloc(2);
    //gsl_vector_set_all(initGuess, 5000);
    //gsl_vector_set_all(initGuess, 1);
    //gsl_vector_set_all(initGuess, 10);
    gsl_vector_set_all(initGuess, 20000);
    //gsl_vector_set(initGuess, 0, 0.01); // cnaToMutRateMu
    //gsl_vector_set(initGuess, 0, 10); // cnaToMutRateMu
    //gsl_vector_set(initGuess, 1, 9.686635); // mutOverdispOmega # value learned from Navin 2011 data
    //gsl_vector* initGuess = gsl_vector_alloc(1);
    gsl_vector_set(initGuess, 1, 10); // mutOverdispOmega
    //gsl_vector_set(initGuess, 1, 0.1); // mutOverdispOmega
    //gsl_vector_set(initGuess, 1, 100); // mutOverdispOmega
    //gsl_vector_set(initGuess, 1, 1000); // mutOverdispOmega
    if(this->verbose) {
      std::lock_guard<std::mutex> lock(Optimizable::mtx); // from https://stackoverflow.com/a/18277334
      std::cout << "Calling BFGS for MutationJointMutRateOverdisp (mu and omega):" << std::endl;
    }
    this->mutJointMutRateOverdisp->bfgs(initGuess, 500, this->verbose, this->gradientDebug);
    if(!this->mutJointMutRateOverdisp->getOptimSuccess() || compareDoubles(0, this->mutJointMutRateOverdisp->getChangeInBFGSLoglikelihood())) {
      if(this->verbose) {
        std::lock_guard<std::mutex> lock(Optimizable::mtx); // from https://stackoverflow.com/a/18277334
        std::cout << "Calling BFGS for MutationJointMutRateOverdisp (mu and omega) FAILED, retrying:" << std::endl;
      }
      gsl_vector_set(initGuess, 0, 20);
      gsl_vector_set(initGuess, 1, 2);
      this->mutJointMutRateOverdisp->bfgs(initGuess, 500, this->verbose, this->gradientDebug);
    }
    if(!this->mutJointMutRateOverdisp->getOptimSuccess() || compareDoubles(0, this->mutJointMutRateOverdisp->getChangeInBFGSLoglikelihood())) {
      if(this->verbose) {
        std::lock_guard<std::mutex> lock(Optimizable::mtx); // from https://stackoverflow.com/a/18277334
        std::cout << "Calling BFGS for MutationJointMutRateOverdisp (mu and omega) FAILED, retrying:" << std::endl;
      }
      gsl_vector_set(initGuess, 0, 20000);
      gsl_vector_set(initGuess, 1, 20);
      this->mutJointMutRateOverdisp->bfgs(initGuess, 500, this->verbose, this->gradientDebug);
    }
    //this->mutJointMutRateOverdisp->bfgs(initGuess, 1, this->verbose, this->gradientDebug);
    gsl_vector_free(initGuess);
  }

  // store omega in calling object and in MutationLists directly (cnaToMutRate doesn't belong to MutationList class)
  this->cnaToMutRateMu = this->mutJointMutRateOverdisp->getCnaToMutRateMu();
  this->mutOverdispOmega = this->mutJointMutRateOverdisp->getMutOverdispOmega();
  //this->mutOverdispOmega = 20; // debugging
  for(unsigned int mutListIdx = 0; mutListIdx < this->mutListVec->size(); mutListIdx++) {
    (*this->mutListVec)[mutListIdx]->setMutOverdispOmega(this->mutOverdispOmega);
  }

  // copy this->mutJointMutRateOverdisp->mutIndVec into this->mutIndVec
  this->mutIndVec = this->mutJointMutRateOverdisp->mutIndVec;
  for(unsigned int mutIndIdx = 0; mutIndIdx < this->mutIndVec->size(); mutIndIdx++) {
    (*this->mutIndVec)[mutIndIdx]->setCnaToMutRateMu(this->cnaToMutRateMu);
  }

  this->allIndEstMutCounts = this->getAllIndEstMutCounts();

  //double jointLl = this->mutJointMutRateOverdisp->getLogLikelihood(); // no need to optim mut counts again, already set
  if(this->verbose) {
    std::lock_guard<std::mutex> lock(Optimizable::mtx); // from https://stackoverflow.com/a/18277334
    fprintf(stdout, "JOINT MUTATION PARAMETERS cnaToMutRateMu and mutOverdispOmega:\n%.40f\n%.40f\n", this->cnaToMutRateMu, this->mutOverdispOmega);
    fprintf(stdout, "\nallIndEstMutCounts:\n");
    printColVector(this->allIndEstMutCounts);
    //std::cout << "\nJOINT MUTATION PARAMETERS loglikelihood: " << jointLl << std::endl;
  }
}

/*
 * function to estimate, set, and save mu and omega parameters using independent cell results from sconce
 */
void IndThenPairs2Stages3TrParam2DegPolyHMMWithMuts::estimateMutParamsFromIndCells(std::string sconceEstimatesPath, std::string filename) {
  //// est cnaToMutRateMu and mutOverdispOmega separately
  //this->estCnaToMutRateMuIndependentlyForAllIndCells(); // mu = median(mu_1, mu_2, ...)
  //this->estMutOverdispOmegaJointlyAcrossAllIndCells(); // omega = one param optimized over all cells

  // reest t's with same median(beta) and median(lambda)
  this->reestIndBranchLengths();

  // est cnaToMutRateMu and mutOverdispOmega jointly across all cels (ie global ests for mu and omega)
  this->estCnaToMutRateMuAndMutOverdispOmegaJointlyAcrossAllIndCells(sconceEstimatesPath);

  //this->print(stdout);
  // save mu and omega params to file
  for(unsigned int mutIndIdx = 0; mutIndIdx < this->mutListVec->size(); mutIndIdx++) {
    std::string cellName = (*this->sampleList)[mutIndIdx];
    boost::replace_all(cellName, ",", "__");
    //std::cout << "about to saveMutParamsToFile:" << std::endl;
    //(*this->mutIndVec)[mutIndIdx]->print(stdout);
    (*this->mutIndVec)[mutIndIdx]->saveMutParamsToFile(0, filename + "__" + cellName + "__k" + std::to_string(this->MAX_PLOIDY) + ".sconceMutParams"); // each MutationInd holds 1 MutationList, so always working with 0'th mutListIdx
    //this->mutJointOverdisp->saveMutParamsToFile(mutIndIdx, filename + "__" + cellName + "__k" + std::to_string(this->MAX_PLOIDY) + ".sconceJointMutParams");
  }
}

void IndThenPairs2Stages3TrParam2DegPolyHMMWithMuts::reestIndBranchLengths() {
  // create OneCellFixLib0TrParam2DegPolyHMM objects
  std::vector<OneCellFixLib0TrParam2DegPolyHMM*>* reestOneCellTsVec = new std::vector<OneCellFixLib0TrParam2DegPolyHMM*>(this->ind2StagesVec->size());
  gsl_vector* indOptimSummaryTrParams = this->getIndOptimSummaryTrParams();
  double medianBeta = gsl_vector_get(indOptimSummaryTrParams, 0);
  double medianLambda = gsl_vector_get(indOptimSummaryTrParams, 1);
  gsl_vector* currFixedParams = nullptr;
  gsl_vector* branchEst = gsl_vector_alloc(1);
  for(unsigned int cellIdx = 0; cellIdx < this->ind2StagesVec->size(); cellIdx++) {
    AllInd3TrParam2DegPolyHMM* currStage1Results = (*this->ind2StagesVec)[cellIdx]->getBFGSAllInd();
    HMM* currHmm = (*currStage1Results->getHMMs())[0];
    std::vector<DepthPair*>* depths = currHmm->getDepths();
    currFixedParams = gsl_vector_alloc(4);
    gsl_vector_set(currFixedParams, 0, currHmm->getLibScalingFactor(0));
    gsl_vector_set(currFixedParams, 1, currHmm->getAlpha());
    gsl_vector_set(currFixedParams, 2, medianBeta);
    gsl_vector_set(currFixedParams, 3, medianLambda);
    gsl_vector_set(branchEst, 0, gsl_vector_get(currHmm->getParamsToEst(), currHmm->BRANCH_LENGTH_START_IDX));
    gsl_vector* meanVar = currHmm->getMeanVarianceCoefVec();
    OneCellFixLib0TrParam2DegPolyHMM* newHmm = new OneCellFixLib0TrParam2DegPolyHMM(depths, currFixedParams, currHmm->getKploidy());
    newHmm->setMeanVarianceFn(meanVar);
    newHmm->setParamsToEst(branchEst);
    (*reestOneCellTsVec)[cellIdx] = newHmm;
    //std::cout << "setting up branch length reest for cellIdx " << cellIdx << std::endl;
    //newHmm->print(stdout);
  }



  // for each cell, optimize t using median beta/lambda
  boost::asio::thread_pool pool(this->numThreads);
  for(unsigned int cellIdx = 0; cellIdx < reestOneCellTsVec->size(); cellIdx++) {

    // call bfgs, in parallel?
    HMM* currHmm = (*reestOneCellTsVec)[cellIdx];
    gsl_vector* initGuess = gsl_vector_alloc(1);
    gsl_vector_memcpy(initGuess, currHmm->getParamsToEst());
    //currHmm->bfgs(initGuess, 500, this->verbose, this->gradientDebug);
    //currHmm->print(stdout);

    // Submit a lambda object to the pool
    boost::asio::post(pool,
      [cellIdx, currHmm, initGuess, this]() { // see https://stackoverflow.com/a/7627218 for lambda function syntax explanation
        if(this->gradientDebug) {
          std::lock_guard<std::mutex> lock(Optimizable::mtx); // from https://stackoverflow.com/a/18277334
          std::cout << "Calling BFGS for reestimating independent branch lengths using common transition parameters for cellIdx " << cellIdx << std::endl;
        }
        currHmm->bfgs(initGuess, 500, this->verbose, this->gradientDebug);
    });
  }
  // Wait for all tasks in the pool to complete.
  pool.join();

  // save branch lengths
  gsl_vector* currParamsToEst = gsl_vector_alloc(4);
  currFixedParams = gsl_vector_alloc(1);
  for(unsigned int cellIdx = 0; cellIdx < reestOneCellTsVec->size(); cellIdx++) {
    AllInd3TrParam2DegPolyHMM* currStage1Results = (*this->ind2StagesVec)[cellIdx]->getBFGSAllInd();
    HMM* optimThmm = (*reestOneCellTsVec)[cellIdx];
    //std::cout << "optimThmm, cellIdx " << cellIdx << std::endl;
    //optimThmm->print(stdout);
    gsl_vector_set(currParamsToEst, 0, optimThmm->getLibScalingFactor(0));
    gsl_vector_set(currParamsToEst, 1, medianBeta);
    gsl_vector_set(currParamsToEst, 2, medianLambda);
    gsl_vector_set(currParamsToEst, 3, gsl_vector_get(optimThmm->getParamsToEst(), 0));
    gsl_vector_set(currFixedParams, 0, optimThmm->getAlpha());
    currStage1Results->setParamsToEstFromIthHMM(currParamsToEst, 0);
    currStage1Results->setFixedParamsFromIthHMM(currFixedParams, 0);
    (*currStage1Results->getHMMs())[0]->setParamsToEst(currParamsToEst);
    (*currStage1Results->getHMMs())[0]->setFixedParams(currFixedParams);
    //std::cout << "set to currParamsToEst, currFixedParams" << std::endl;
    //printRowVector(currParamsToEst);
    //printRowVector(currFixedParams);
    //std::cout << "currStage1Resutls" << std::endl;
    //currStage1Results->print(stdout);
  }
  gsl_vector_free(currParamsToEst);
  gsl_vector_free(branchEst);
  gsl_vector_free(currFixedParams);
  //std::cout << "saving branch lengths" << std::endl;
  //this->print(stdout);
}


/*
 * function to copy estimates from indv cell estimation into pairs. defers to parent class implementation, then additionally sets mu values
 */
void IndThenPairs2Stages3TrParam2DegPolyHMMWithMuts::copyStage1EstsIntoStage2FixedParams(gsl_vector* stage2FixedParams, gsl_vector* stage2InitGuess, gsl_vector* indEstParams, int numStage2LibsToFix) {
  // first call parent class
  IndThenPairs2Stages3TrParam2DegPolyHMM::copyStage1EstsIntoStage2FixedParams(stage2FixedParams, stage2InitGuess, indEstParams, numStage2LibsToFix);

  // then fill in cnaToMutRateMu as the last fixed param (assumes there's only one param). omega is already stored in each MutationList
  gsl_vector_set(stage2FixedParams, stage2FixedParams->size - 1, this->cnaToMutRateMu);
}

//void IndThenPairs2Stages3TrParam2DegPolyHMMWithMuts::estimateAllMutationPairMutCounts(std::string filename) {
//  this->stage2Pairs->estimateAllMutationPairMutCounts(filename);
//}

void IndThenPairs2Stages3TrParam2DegPolyHMMWithMuts::setUpAllPairsBFGS(gsl_vector* fixedParams, gsl_vector* initGuess, gsl_vector* meanVarianceCoefVec, bool fixLib, bool preallocIntermediates) {
  if(fixLib) {
    this->stage2Pairs = AllPairsFixLib0TrParam2DegPolyHMMWithMuts::create(this->depthsVec, this->sampleList, this->mutListVec, this->allIndEstMutCounts, fixedParams, this->MAX_PLOIDY, this->NUM_PAIRS_STAGE_2, meanVarianceCoefVec, preallocIntermediates);
  }
  else {
    this->stage2Pairs = AllPairs0TrParam2DegPolyHMMWithMuts::create(this->depthsVec, this->sampleList, this->mutListVec, this->allIndEstMutCounts, fixedParams, this->MAX_PLOIDY, this->NUM_PAIRS_STAGE_2, meanVarianceCoefVec, preallocIntermediates);
  }
  this->stage2Pairs->gradientDebug = this->gradientDebug;
  this->stage2Pairs->setParamsToEst(initGuess);
  this->stage2Pairs->setNumThreads(this->numThreads);
  this->stage2Pairs->setCentralDiffFlag(centralDiff);
  //std::cout << "end of IndThenPairs2Stages3TrParam2DegPolyHMMWithMuts::setUpAllPairsBFGS:" << std::endl;
  //this->stage2Pairs->print(stdout);

}
void IndThenPairs2Stages3TrParam2DegPolyHMMWithMuts::setUpSelectPairsBFGS(gsl_vector* fixedParams, gsl_vector* initGuess, int numPairsToSummarize, int ordering, int seed, gsl_vector* meanVarianceCoefVec, bool fixLib, bool preallocIntermediates) {
  if(fixLib) {
    this->stage2Pairs = SelectPairsFixLib0TrParam2DegPolyHMMWithMuts::create(this->depthsVec, this->sampleList, this->mutListVec, this->allIndEstMutCounts, fixedParams, this->MAX_PLOIDY, this->NUM_PAIRS_STAGE_2, numPairsToSummarize, ordering, seed, this->stage1NearestCellIdxMat, meanVarianceCoefVec, preallocIntermediates);
  }
  else {
    this->stage2Pairs = SelectPairs0TrParam2DegPolyHMMWithMuts::create(this->depthsVec, this->sampleList, this->mutListVec, this->allIndEstMutCounts, fixedParams, this->MAX_PLOIDY, this->NUM_PAIRS_STAGE_2, numPairsToSummarize, ordering, seed, this->stage1NearestCellIdxMat, meanVarianceCoefVec, preallocIntermediates);
  }
  this->stage2Pairs->gradientDebug = this->gradientDebug;
  this->stage2Pairs->setParamsToEst(initGuess);
  this->stage2Pairs->setNumThreads(this->numThreads);
  this->stage2Pairs->setCentralDiffFlag(this->centralDiff);
}
//void IndThenPairs2Stages3TrParam2DegPolyHMMWithMuts::optimIndCells(gsl_vector* initGuess, std::string filename, int numBWIters, int numLibStarts, int libStartVal, int maxBFGSIters, bool verbose, bool debug) {
//  // don't think need to override. can just tack on extra step at end of copyStage1EstsIntoStage2FixedParams
//}
AllPairs0TrParam2DegPolyHMM* IndThenPairs2Stages3TrParam2DegPolyHMMWithMuts::bfgsStage2(gsl_vector* initGuess, int maxIters, std::string filename, bool verbose, bool debug) {
  //std::cout << "HERE starting bfgsStage2" << std::endl;
  //this->print(stdout);
  IndThenPairs2Stages3TrParam2DegPolyHMM::bfgsStage2(initGuess, maxIters, filename, verbose, debug);
  this->allPairedEstMutCounts = this->getAllPairedEstMutCounts();
  return this->stage2Pairs;
}

