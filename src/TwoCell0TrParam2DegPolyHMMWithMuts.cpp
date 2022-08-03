#include "TwoCell0TrParam2DegPolyHMMWithMuts.hpp"

TwoCell0TrParam2DegPolyHMMWithMuts::TwoCell0TrParam2DegPolyHMMWithMuts(std::vector<DepthPair*>* depths, MutationPair* mutPair, gsl_vector* fixedParams, int maxPloidy, bool preallocIntermediates) : TwoCell0TrParam2DegPolyHMMWithMuts(depths, mutPair, fixedParams, maxPloidy, 0, 3, 0, 3, preallocIntermediates) { // 0 transition params to est, 3 fixedTrParams, 0 fixedLibs, 3 branches
}

TwoCell0TrParam2DegPolyHMMWithMuts::TwoCell0TrParam2DegPolyHMMWithMuts(std::vector<DepthPair*>* depths, MutationPair* mutPair, gsl_vector* fixedParams, int maxPloidy, int numTrParamsToEst, int numFixedTrParams, int numFixedLibs, int numBranches, bool preallocIntermediates) : TwoCell0TrParam2DegPolyHMM(depths, fixedParams, maxPloidy, numTrParamsToEst, numFixedTrParams, numFixedLibs, numBranches, preallocIntermediates) {
  //this->cnaToMutRateMu = 1;
  this->mutPair = mutPair;
}
TwoCell0TrParam2DegPolyHMMWithMuts::~TwoCell0TrParam2DegPolyHMMWithMuts() {
  // TODO
}

double TwoCell0TrParam2DegPolyHMMWithMuts::getCnaToMutRateMu() {
  //return this->cnaToMutRateMu;
  return gsl_vector_get(this->fixedParams, this->FIXED_MUTATION_PARAM_START_IDX + 0);
}
void TwoCell0TrParam2DegPolyHMMWithMuts::setCnaToMutRateMu(double rate) {
  //this->cnaToMutRateMu = rate;
  gsl_vector_set(this->fixedParams, this->FIXED_MUTATION_PARAM_START_IDX + 0, rate);
}
TwoCell0TrParam2DegPolyHMMWithMuts* TwoCell0TrParam2DegPolyHMMWithMuts::bfgs(gsl_vector* initGuess, int maxIters, bool verbose, bool debug) {
  TwoCell0TrParam2DegPolyHMMWithMuts* bestGuessHMM = this;
  gsl_vector* initGuessAsParams = gsl_vector_alloc(initGuess->size);
  this->convertProbToParam(initGuessAsParams, initGuess);
  Optimizable::bfgs(initGuessAsParams, bestGuessHMM, maxIters, verbose, debug);
  gsl_vector_free(initGuessAsParams);
  return bestGuessHMM;
}

void TwoCell0TrParam2DegPolyHMMWithMuts::print(FILE* stream) {
  HMM::print(stream);
  fprintf(stream, "TwoCell0TrParam2DegPolyHMMWithMuts MutationPair:\n");
  this->mutPair->print(stream);
}


double TwoCell0TrParam2DegPolyHMMWithMuts::getLogLikelihood() {
  double cnaLl = this->runForwardAlg();
  double mutLl = this->getMutLogLikelihood();
  return cnaLl + mutLl;
}
/*
 * copied from TwoCellFixLib0TrParam2DegPolyHMMWithMuts
 * function to calculate the genotype likelihoods of somatic mutations
 */
double TwoCell0TrParam2DegPolyHMMWithMuts::getMutLogLikelihood() {
  double mutLl = 0;
  double currT = 0;
  double currX = 0;
  double currLl = 0;
  double cnaToMutRateMu = this->getCnaToMutRateMu();
  // sum up loglik of prob of counts of muts
  for(int branchIdx = 0; branchIdx < this->NUM_BRANCH_LENGTHS_TO_EST; branchIdx++) {
    currT = gsl_vector_get(this->paramsToEst, this->BRANCH_LENGTH_START_IDX + branchIdx);
    currX = this->mutPair->getMutCountEst(branchIdx);
    currLl = currX * log(cnaToMutRateMu * currT) - cnaToMutRateMu * currT - lgamma(currX+1); // poisson dist, with mean = cnaToMutRateMu * currT
    mutLl += currLl;
  }
  return mutLl;
}

/*
 * copied from TwoCellFixLib0TrParam2DegPolyHMMWithMuts
 * function to estimate X1/X2/X3
 */
void TwoCell0TrParam2DegPolyHMMWithMuts::estMutCountsPerBranch() {
  // first run viterbi decoding so have current copy number ests
  this->viterbiDecode();

  // then est counts on each branch (ie optimize ll(X1,X2,X3) to get ests of X1,X2,X3)
  gsl_vector* initGuess = gsl_vector_alloc(this->NUM_BRANCH_LENGTHS_TO_EST);
  this->mutPair->setInitGuessNthTime(initGuess, 0, 0);
  this->mutPair->estMutCountsPerBranch(this->chrToViterbiPathMapVec, initGuess, 500, this->gradientDebug);
}

