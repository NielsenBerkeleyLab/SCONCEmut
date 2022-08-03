#include "TwoCellFixLib0TrParam2DegPolyHMMWithMuts.hpp"

/*
 ********
 * constructors and destructor
 ********
 */
TwoCellFixLib0TrParam2DegPolyHMMWithMuts::TwoCellFixLib0TrParam2DegPolyHMMWithMuts(std::vector<DepthPair*>* depths, MutationPair* mutPair, gsl_vector* fixedParams, int maxPloidy, bool preallocIntermediates) : TwoCellFixLib0TrParam2DegPolyHMMWithMuts(depths, mutPair, fixedParams, maxPloidy, 0, 3, 2, 3, preallocIntermediates) { // fixedParams, 0 transition params to est, 3 fixedTrParams, 2 fixedLibs, 3 branches
}
TwoCellFixLib0TrParam2DegPolyHMMWithMuts::TwoCellFixLib0TrParam2DegPolyHMMWithMuts(std::vector<DepthPair*>* depths, MutationPair* mutPair, gsl_vector* fixedParams, int maxPloidy, int numTrParamsToEst, int numFixedTrParams, int numFixedLibs, int numBranches, bool preallocIntermediates) : TwoCellFixLib0TrParam2DegPolyHMM(depths, fixedParams, maxPloidy, numTrParamsToEst, numFixedTrParams, numFixedLibs, numBranches, preallocIntermediates) {
  //this->cnaToMutRateMu = 1;
  //std::cout << "TwoCellFixLib0TrParam2DegPolyHMMWithMuts::TwoCellFixLib0TrParam2DegPolyHMMWithMuts ctor" << std::endl;
  //mutPair->print(stdout);
  this->mutPair = mutPair;
}

TwoCellFixLib0TrParam2DegPolyHMMWithMuts::~TwoCellFixLib0TrParam2DegPolyHMMWithMuts() {
}

/*
 ********
 * accessors and mutators
 ********
 */
double TwoCellFixLib0TrParam2DegPolyHMMWithMuts::getCnaToMutRateMu() {
  //return this->cnaToMutRateMu;
  return gsl_vector_get(this->fixedParams, this->FIXED_MUTATION_PARAM_START_IDX + 0);
}

void TwoCellFixLib0TrParam2DegPolyHMMWithMuts::setCnaToMutRateMu(double rate) {
  //this->cnaToMutRateMu = rate;
  gsl_vector_set(this->fixedParams, this->FIXED_MUTATION_PARAM_START_IDX + 0, rate);
}

TwoCellFixLib0TrParam2DegPolyHMMWithMuts* TwoCellFixLib0TrParam2DegPolyHMMWithMuts::bfgs(gsl_vector* initGuess, int maxIters, bool verbose, bool debug) {
  TwoCellFixLib0TrParam2DegPolyHMMWithMuts* bestGuessHMM = this;
  gsl_vector* initGuessAsParams = gsl_vector_alloc(initGuess->size);
  bestGuessHMM->convertProbToParam(initGuessAsParams, initGuess);
  Optimizable::bfgs(initGuessAsParams, bestGuessHMM, maxIters, verbose, debug);
  gsl_vector_free(initGuessAsParams);
  return bestGuessHMM;
}

void TwoCellFixLib0TrParam2DegPolyHMMWithMuts::print(FILE* stream) {
  HMM::print(stream);
  fprintf(stream, "TwoCellFixLib0TrParam2DegPolyHMMWithMuts MutationPair:\n");
  this->mutPair->print(stream);

  //fprintf(stream, "before manual ll: %0.20f\n", this->getLogLikelihood());

  //gsl_vector* trueMutCounts = gsl_vector_alloc(3);
  //gsl_vector_set(trueMutCounts, 0, 1816);
  //gsl_vector_set(trueMutCounts, 1, 486);
  //gsl_vector_set(trueMutCounts, 2, 565);
  //
  //gsl_vector* truebranchLengths = gsl_vector_alloc(3);
  //gsl_vector_set(truebranchLengths, 0, 0.6093750);
  //gsl_vector_set(truebranchLengths, 1, 0.3828125);
  //gsl_vector_set(truebranchLengths, 2, 0.3828125);
  //
  //this->mutPair->setParamsToEst(trueMutCounts);
  //this->setParamsToEst(truebranchLengths);
  //std::cout << "MANUAL SET: MUTPAIR AND THIS->GETLIKELIHOOD" << std::endl;
  //this->mutPair->print(stream);
  //fprintf(stream, "manual ll: %0.20f\n", this->getLogLikelihood());
}

double TwoCellFixLib0TrParam2DegPolyHMMWithMuts::getLogLikelihood() {
  double cnaLl = this->runForwardAlg();
  //this->estMutCountsPerBranch(); // if want to reest X1/X2/X3 every BFGS iteration
  double mutLl = this->getMutLogLikelihood(); // different from MutationPair likelihood because this calcs prob of branch lengths with MutationPair mutations
  std::cout << "TwoCellFixLib0TrParam2DegPolyHMMWithMuts::getLogLikelihood: " << cnaLl << " + " << mutLl << " = " << cnaLl + mutLl << std::endl;
  return cnaLl + mutLl;
}

/*
 * function to calculate the likelihood of counts of somatic mutations on each branch, t1/t2/t3
 * counts have a Poisson-like dist, with mean = cnaToMutRateMu * branchLength
 * (Poisson-like instead of Poisson bc mean may be fractional)
 */
double TwoCellFixLib0TrParam2DegPolyHMMWithMuts::getMutLogLikelihood() {
  double mutLl = 0;
  double currT = 0;
  double currX = 0;
  double currLl = 0;
  double cnaToMutRateMu = this->getCnaToMutRateMu();
  // sum up loglik of prob of counts of muts
  for(int branchIdx = 0; branchIdx < this->NUM_BRANCH_LENGTHS_TO_EST; branchIdx++) {
    currT = gsl_vector_get(this->paramsToEst, this->BRANCH_LENGTH_START_IDX + branchIdx);
    currX = this->mutPair->getMutCountEst(branchIdx);
    //currLl = log(pow(this->cnaToMutRateMu * currT, currX) * exp(-this->cnaToMutRateMu * currT) / tgamma(currX+1));
    currLl = currX * log(cnaToMutRateMu * currT) - cnaToMutRateMu * currT - lgamma(currX+1); // poisson dist, with mean = cnaToMutRateMu * currT
    mutLl += currLl;
  }
  return mutLl;
}

/*
 * function to estimate X1/X2/X3. Assumes HMM is fully set up to run viterbiDecode successfully (ie transition params and meanVar need to be set; only a problem for creating these HMM's individually in testing)
 */
void TwoCellFixLib0TrParam2DegPolyHMMWithMuts::estMutCountsPerBranch() {
  if(this->gradientDebug) {
    std::lock_guard<std::mutex> lock(Optimizable::mtx); // from https://stackoverflow.com/a/18277334
    //std::cout << "Calling BFGS for MutationPair to estimate mutation counts" << std::endl;
    std::cout << "Calling SIMPLEX for MutationPair to estimate mutation counts" << std::endl;
  }
  //std::cout << "TwoCellFixLib0TrParam2DegPolyHMMWithMuts::estMutCountsPerBranch" << std::endl;
  // first run viterbi decoding so have current copy number ests
  this->allocDecodingIntermediates(); // if reest X1/X2/X3 on ever iter of BFGS, need viterbi decoding
  this->viterbiDecode();

  // then est counts on each branch (ie optimize ll(X1,X2,X3) to get ests of X1,X2,X3)
  //this->mutPair->setInitGuessNthTime(this->mutPair->internalInitGuess, 0, 0);
  gsl_vector_memcpy(this->mutPair->internalInitGuess, this->mutPair->getParamsToEst());
  //this->mutPair->estMutCountsPerBranch(this->chrToViterbiPathMapVec, this->mutPair->internalInitGuess, 500, this->gradientDebug);
}


/*
 * function to save param estimates (including mutation params) to file
 */
void TwoCellFixLib0TrParam2DegPolyHMMWithMuts::saveParamEstimates(std::string filename) const {
  HMM::saveParamEstimates(filename);
  mutPair->saveMutParamsToFile(0, filename, "a"); // the only index dependent variable is cnaToMutRateMu, which should be constant across all MutationList objects, so just pass 0. "a" for append
}

