#include "MutationInd.hpp"

MutationInd::MutationInd(std::vector<MutationList*>* mutListVec, bool verbose, bool gradientDebug) : MutationListContainer(mutListVec, 1, 0, verbose, gradientDebug) { // 1 mutCount, 0 overdispVars
  this->cnaToMutRateMu = 0;
  this->hasMutEsts = false;
}

MutationInd::~MutationInd() {
}
MutationListContainer* MutationInd::estMutCountsPerBranch(std::vector<std::unordered_map<std::string, std::vector<int>*>*>* chrToViterbiPathMapVec, gsl_vector* initGuess, int maxIters, bool verbose) {
  if(this->hasMutEsts && compareDoubles(0, this->changeInBFGSLoglikelihood)) {
    return this;
  }
  this->hasMutEsts = true;
  return MutationListContainer::estMutCountsPerBranch(chrToViterbiPathMapVec, initGuess, maxIters, verbose);
}

double MutationInd::getMutCountEst() {
  return this->getMutCountEst(0);
}

void MutationInd::setCnaToMutRateMu(double mu) {
  this->cnaToMutRateMu = mu;
}
double MutationInd::getCnaToMutRateMu() {
  return this->cnaToMutRateMu;
}
double MutationInd::getMutOverdispOmega(int mutListIdx) {
  return (*this->mutListVec)[mutListIdx]->getMutOverdispOmega();
}

/*
 * saves ([mu, omega] ordering matches MutationJointMutRateOverdisp):
 * this->cnaToMutRateMu
 * MutationList[mutListIdx]::mutOverdispOmega
 * paramsToEst from calling object
 */
void MutationInd::saveMutParamsToFile(int mutListIdx, std::string filename) {
  gsl_vector* writingBuffer = gsl_vector_alloc(1 + this->paramsToEst->size + 1); // MutationList::mutOverdispOmega + paramsToEst + this->cnaToMutRateMu
  int bufferIdx = 0;

  // store this->cnaToMutRateMu
  gsl_vector_set(writingBuffer, bufferIdx, this->cnaToMutRateMu);
  bufferIdx++;

  // store MutationList::mutOverdispOmega
  gsl_vector_set(writingBuffer, bufferIdx, (*this->mutListVec)[mutListIdx]->getMutOverdispOmega());
  bufferIdx++;

  // store paramsToEst
  for(unsigned int paramIdx = 0; paramIdx < this->paramsToEst->size; paramIdx++, bufferIdx++) {
    gsl_vector_set(writingBuffer, bufferIdx, gsl_vector_get(this->paramsToEst, paramIdx));
  }

  FILE* outFile = fopen(filename.c_str(), "w");
  gsl_block_fprintf(outFile, writingBuffer->block, "%.40f");
  fclose(outFile);

  gsl_vector_free(writingBuffer);
}

/*
 * counterpart to MutationInd::saveMutParamsToFile
 * stores cnaToMutRateMu, MutationList::mutOverdispOmega, and paramsToEst
 */
double MutationInd::getMutParamsFromFile(std::string filename, int numExpectedLinesPerFile) {
  // construct paramsToEst from the readingBuffer, then call set methods
  int numParamsToRead = numExpectedLinesPerFile;
  gsl_vector* readingBuffer = gsl_vector_alloc(numParamsToRead); // mu/omega/x
  gsl_vector* readParamsToEst = gsl_vector_alloc(this->paramsToEst->size);
  int bufferIdx = 0;

  // read in file
  FILE* currFile = fopen(filename.c_str(), "r");
  gsl_block_fscanf(currFile, readingBuffer->block);
  fclose(currFile);

  // store this->cnaToMutRateMu
  this->setCnaToMutRateMu(gsl_vector_get(readingBuffer, bufferIdx));
  bufferIdx++;

  // store MutationList::mutOverdispOmega
  this->setAllMutOverdispOmega(gsl_vector_get(readingBuffer, bufferIdx));
  bufferIdx++;

  // store paramsToEst
  for(unsigned int paramIdx = 0; paramIdx < this->paramsToEst->size; paramIdx++, bufferIdx++) {
    gsl_vector_set(readParamsToEst, paramIdx, gsl_vector_get(readingBuffer, bufferIdx));
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

/*
 * calculates l(X) = sum_sites [log(p(data_i | S_i=1) * X/N + p(data_i | S_i=0) * (N-X)/N)]
 */
double MutationInd::getLogLikelihood() {
  MutationList* mutList = (*this->mutListVec)[0]; // for convenience
  double totalLl = 0;
  double currLl = 0;
  double x = gsl_vector_get(this->paramsToEst, 0);
  double N = mutList->coordVec->size();
  double prob_s1 = 0; // P(D_i | S_i = 1)
  double prob_s0 = 0;
  std::string currSite;

  // for each site
  for(unsigned int mutIdx = 0; mutIdx < mutList->coordVec->size(); mutIdx++) {
    currSite = (*mutList->coordVec)[mutIdx];

    // calc ll for each type of mutation
    prob_s1 = mutList->getLikelihood(currSite, true);
    prob_s0 = mutList->getLikelihood(currSite, false);
    double prob_s1_norm = prob_s1 / (prob_s1 + prob_s0);
    double prob_s0_norm = prob_s0 / (prob_s1 + prob_s0);

    // then sum up and store
    currLl = log(prob_s1_norm * x/N + prob_s0_norm * (N-x)/N); // P(D|S=1) * X/N + P(D|S=0) * (N-X)/N
    if(gsl_isinf(currLl) || gsl_isnan(currLl)) { // inf from prob_s1/prob_s0 log(0), nan from prob_s1_norm/prob_s1_norm 0/0
      currLl = 0;
    }
    totalLl += currLl;
  }
  return totalLl;
}

// ex estimated x shouldn't be greater than N
double MutationInd::getValidOptimParamMax() const {
  return (*this->mutListVec)[0]->coordVec->size();
}

