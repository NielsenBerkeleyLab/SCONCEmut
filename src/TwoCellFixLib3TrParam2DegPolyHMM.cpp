#include "TwoCellFixLib3TrParam2DegPolyHMM.hpp"


TwoCellFixLib3TrParam2DegPolyHMM::TwoCellFixLib3TrParam2DegPolyHMM(std::vector<DepthPair*>* depths, gsl_vector* fixedParams, int maxPloidy, bool preallocIntermediates) : TwoCellFixLib3TrParam2DegPolyHMM(depths, fixedParams, maxPloidy, 2, 1, 2, 3, preallocIntermediates) { // fixedParams, 2 numTrParamsToEst (beta, lambda), 1 numFixedTrParams (alpha), 2 numFixedLibs, 3 numBranches
}
TwoCellFixLib3TrParam2DegPolyHMM::TwoCellFixLib3TrParam2DegPolyHMM(std::vector<DepthPair*>* depths, gsl_vector* fixedParams, int maxPloidy, int numTrParamsToEst, int numFixedTrParams, int numFixedLibs, int numBranches, bool preallocIntermediates) : TwoCell3TrParam2DegPolyHMM(depths, fixedParams, maxPloidy, numTrParamsToEst, numFixedTrParams, numFixedLibs, numBranches, preallocIntermediates) {
  this->maxNumBFGSStarts = 3;

  // set up totalLogEmissionLookup
  std::string currChr;
  int currNumWindows = -1;
  std::vector<std::string>* chrVec = this->getChrVec();
  this->totalLogEmissionLookup = new std::vector<gsl_matrix*>();
  this->totalEmissionLookup = new std::vector<gsl_matrix*>();
  for(unsigned int chrIdx = 0; chrIdx < chrVec->size(); chrIdx++) {
    currChr = (*chrVec)[chrIdx];
    currNumWindows = (*(*this->depths)[0]->regions)[currChr]->size();
    //this->totalLogEmissionLookup->push_back(gsl_matrix_alloc(this->states->size(), currNumWindows));
    this->totalLogEmissionLookup->push_back(gsl_matrix_alloc(currNumWindows, this->states->size())); // transposed
    this->totalEmissionLookup->push_back(gsl_matrix_alloc(currNumWindows, this->states->size())); // transposed
    gsl_matrix_set_zero((*this->totalLogEmissionLookup)[chrIdx]);
    gsl_matrix_set_zero((*this->totalEmissionLookup)[chrIdx]);
    //printMatrix((*this->totalLogEmissionLookup)[chrIdx]);
  }
  //std::cout << "done with init totalLogEmissionLookup setup" << std::endl;

  /*// set fixedParams with library sizes
  double totalAvgDiploidDepth = (*this->depths)[0]->getTotalDiploidDepth();
  this->setLibScalingFactor(0, (*this->depths)[0]->getTotalTumorDepth() / totalAvgDiploidDepth);
  //std::cout << "done with fixed ctor update cell 0" << std::endl;
  this->setLibScalingFactor(1, (*this->depths)[1]->getTotalTumorDepth() / totalAvgDiploidDepth);
  //std::cout << "done with fixed ctor update cell 1" << std::endl;*/
}

TwoCellFixLib3TrParam2DegPolyHMM::~TwoCellFixLib3TrParam2DegPolyHMM() {
  for(std::vector<gsl_matrix*>::iterator it = this->totalLogEmissionLookup->begin(); it != this->totalLogEmissionLookup->end(); ++it) {
    gsl_matrix_free(*it);
  }
  for(std::vector<gsl_matrix*>::iterator it = this->totalEmissionLookup->begin(); it != this->totalEmissionLookup->end(); ++it) {
    gsl_matrix_free(*it);
  }
  delete this->totalLogEmissionLookup;
  delete this->totalEmissionLookup;
}

void TwoCellFixLib3TrParam2DegPolyHMM::setLibScalingFactor(int cellNum, double libScalingFactor) {
  gsl_vector_set(this->fixedParams, this->LIB_SIZE_SCALING_FACTOR_START_IDX + cellNum, libScalingFactor);
  this->updateAllTotalEmissionLookup();
}
double TwoCellFixLib3TrParam2DegPolyHMM::getLibScalingFactor(int cellNum) const {
  return gsl_vector_get(this->fixedParams, this->LIB_SIZE_SCALING_FACTOR_START_IDX + cellNum);
}
/*
 * currChrDepthsVec is ignored. Assumes updateAllTotalEmissionLookup is called
 * every time lib sizes are updated
 */
double TwoCellFixLib3TrParam2DegPolyHMM::getTotalLogEmissionProb(int stateIdx, std::vector<std::vector<double>*>* currChrDepthsVec, int chrIdx, int depthIdx) {
  gsl_matrix* currChrEmissionMat = (*this->totalLogEmissionLookup)[chrIdx];
  //return gsl_matrix_get(currChrEmissionMat, stateIdx, depthIdx);
  return gsl_matrix_get(currChrEmissionMat, depthIdx, stateIdx); // transposed
}
double TwoCellFixLib3TrParam2DegPolyHMM::getTotalEmissionProb(int stateIdx, std::vector<std::vector<double>*>* currChrDepthsVec, int chrIdx, int depthIdx) {
  gsl_matrix* currChrEmissionMat = (*this->totalEmissionLookup)[chrIdx];
  return gsl_matrix_get(currChrEmissionMat, depthIdx, stateIdx); // transposed
}

/*
 * updates totalLogEmissionLookup, where each entry corresponds to one chr's matrix.
 * //Each matrix is rows: states, cols: depthIdx. Each entry is calculated by
 * Each matrix is rows: depthIdx, cols: states (trying out since  matrices are stored row major, and we iter through all states in a window before moving to the next window). Each entry is calculated by
 * TwoCell3Tr2DegPolyHMM's getTotalLogEmissionProb
 */
void TwoCellFixLib3TrParam2DegPolyHMM::updateAllTotalEmissionLookup() {
  std::vector<std::string>* chrVec = this->getChrVec();
  std::vector<std::vector<double>*>* currChrDepthsVec = new std::vector<std::vector<double>*>(this->NUM_CELLS + 1);
  DepthPair* firstDepthPair = (*this->depths)[0]; // for convenience
  DepthPair* currDepths = nullptr;
  int cellIdx = 0;
  std::string currChr;
  gsl_matrix* currChrLogEmissionMat = nullptr;
  gsl_matrix* currChrEmissionMat = nullptr;
  double currLogEmi = 0;

  // for each chr
  for(unsigned int chrIdx = 0; chrIdx < chrVec->size(); chrIdx++) {
    currChrLogEmissionMat = (*this->totalLogEmissionLookup)[chrIdx];
    currChrEmissionMat = (*this->totalEmissionLookup)[chrIdx];
    currChr = (*chrVec)[chrIdx];
    for(cellIdx = 0; cellIdx < this->NUM_CELLS; cellIdx++) {
      currDepths = (*this->depths)[cellIdx];
      (*currChrDepthsVec)[cellIdx] = (*currDepths->chrToTumorDepthMap)[currChr];
    }
    (*currChrDepthsVec)[this->NUM_CELLS] = (*firstDepthPair->chrToDiploidDepthMap)[currChr];

    //for(unsigned int stateIdx = 0; stateIdx < currChrEmissionMat->size1; stateIdx++) {
    //  for(unsigned int depthIdx = 0; depthIdx < currChrEmissionMat->size2; depthIdx++) {
    for(unsigned int stateIdx = 0; stateIdx < currChrLogEmissionMat->size2; stateIdx++) { // transposed
      for(unsigned int depthIdx = 0; depthIdx < currChrLogEmissionMat->size1; depthIdx++) {
        //gsl_matrix_set(currChrEmissionMat, stateIdx, depthIdx, TwoCell3TrParam2DegPolyHMM::getTotalLogEmissionProb(stateIdx, currChrDepthsVec, chrIdx, depthIdx));
        currLogEmi = TwoCell3TrParam2DegPolyHMM::getTotalLogEmissionProb(stateIdx, currChrDepthsVec, chrIdx, depthIdx);
        gsl_matrix_set(currChrLogEmissionMat, depthIdx, stateIdx, currLogEmi); // transposed
        gsl_matrix_set(currChrEmissionMat, depthIdx, stateIdx, exp(currLogEmi)); // transposed
      }
    }
  }
  //printMatrix(currChrEmissionMat);
}

void TwoCellFixLib3TrParam2DegPolyHMM::setMeanVarianceFn(gsl_vector* meanVarianceCoefVec) {
  HMM::setMeanVarianceFn(meanVarianceCoefVec);
  this->updateAllTotalEmissionLookup();
}

/*
 * same as TwoCell3TrParam2DegPolyHMM, but lib sizes are skipped
 */
void TwoCellFixLib3TrParam2DegPolyHMM::convertProbToParam(gsl_vector* dest, const gsl_vector* src) const {
  //double d = (double) (*this->depths)[0]->maxWindowSize;
  //double a = this->getAlpha();
  //double b = gsl_vector_get(src, this->TRANSITION_PROB_START_IDX + 0);
  //double g = gsl_vector_get(src, this->TRANSITION_PROB_START_IDX + 1);
  //double t1 = gsl_vector_get(src, this->TRANSITION_PROB_START_IDX + 2) / d;
  //double t2 = gsl_vector_get(src, this->TRANSITION_PROB_START_IDX + 3) / d;
  //double t3 = gsl_vector_get(src, this->TRANSITION_PROB_START_IDX + 4) / d;

  //double c = (b * this->getKploidy() + g - 1);
  //gsl_vector_set(dest, this->TRANSITION_PROB_START_IDX + 0, log(-b * (1-2*a) * (double) this->getKploidy() / c)); // set y
  //gsl_vector_set(dest, this->TRANSITION_PROB_START_IDX + 1, log(-g * (1-2*a) / c)); // set z

  //gsl_vector_set(dest, this->TRANSITION_PROB_START_IDX + 2, log(-(d * t1) / (d * (t1 + t2 + t3) - 1))); // set T1
  //gsl_vector_set(dest, this->TRANSITION_PROB_START_IDX + 3, log(-(d * t2) / (d * (t1 + t2 + t3) - 1))); // set T2
  //gsl_vector_set(dest, this->TRANSITION_PROB_START_IDX + 4, log(-(d * t3) / (d * (t1 + t2 + t3) - 1))); // set T3

  double beta = gsl_vector_get(src, this->TRANSITION_PROB_START_IDX + 0);
  double lambda = gsl_vector_get(src, this->TRANSITION_PROB_START_IDX + 1);
  double t1 = gsl_vector_get(src, this->TRANSITION_PROB_START_IDX + 2);
  double t2 = gsl_vector_get(src, this->TRANSITION_PROB_START_IDX + 3);
  double t3 = gsl_vector_get(src, this->TRANSITION_PROB_START_IDX + 4);

  gsl_vector_set(dest, this->TRANSITION_PROB_START_IDX + 0, log(beta)); // set y
  gsl_vector_set(dest, this->TRANSITION_PROB_START_IDX + 1, log(lambda)); // set z
  gsl_vector_set(dest, this->TRANSITION_PROB_START_IDX + 2, log(t1)); // set T1
  gsl_vector_set(dest, this->TRANSITION_PROB_START_IDX + 3, log(t2)); // set T2
  gsl_vector_set(dest, this->TRANSITION_PROB_START_IDX + 4, log(t3)); // set T3
}

/*
 * reverse of convertProbToParam (see above)
 */
void TwoCellFixLib3TrParam2DegPolyHMM::convertParamToProb(gsl_vector* dest, const gsl_vector* src) const {
  double y = gsl_vector_get(src, this->TRANSITION_PROB_START_IDX + 0);
  double z = gsl_vector_get(src, this->TRANSITION_PROB_START_IDX + 1);
  double T1 = gsl_vector_get(src, this->TRANSITION_PROB_START_IDX + 2);
  double T2 = gsl_vector_get(src, this->TRANSITION_PROB_START_IDX + 3);
  double T3 = gsl_vector_get(src, this->TRANSITION_PROB_START_IDX + 4);
  //double a = this->getAlpha();
  //double c = 1 - 2.0*a + exp(y) + exp(z);

  //gsl_vector_set(dest, this->TRANSITION_PROB_START_IDX + 0, exp(y) / ((double) this->getKploidy() * c)); // beta
  //gsl_vector_set(dest, this->TRANSITION_PROB_START_IDX + 1, exp(z) / c); // gamma
  gsl_vector_set(dest, this->TRANSITION_PROB_START_IDX + 0, exp(y)); // beta
  gsl_vector_set(dest, this->TRANSITION_PROB_START_IDX + 1, exp(z)); // lambda
  //c = 1.0 / (1 + exp(T1) + exp(T2) + exp(T3));

  //gsl_vector_set(dest, this->TRANSITION_PROB_START_IDX + 2, exp(T1) * c); // set t1
  //gsl_vector_set(dest, this->TRANSITION_PROB_START_IDX + 3, exp(T2) * c); // set t2
  //gsl_vector_set(dest, this->TRANSITION_PROB_START_IDX + 4, exp(T3) * c); // set t3
  gsl_vector_set(dest, this->TRANSITION_PROB_START_IDX + 2, exp(T1)); // set t1
  gsl_vector_set(dest, this->TRANSITION_PROB_START_IDX + 3, exp(T2)); // set t2
  gsl_vector_set(dest, this->TRANSITION_PROB_START_IDX + 4, exp(T3)); // set t3
}


/*
 * calls bfgs, returns a bestGuessHMM
 */
TwoCellFixLib3TrParam2DegPolyHMM* TwoCellFixLib3TrParam2DegPolyHMM::bfgs(gsl_vector* initGuess, int maxIters, bool verbose, bool debug) {
  // create new HMM with the best guess parameters and return it
  //TwoCellFixLib3TrParam2DegPolyHMM* bestGuessHMM = new TwoCellFixLib3TrParam2DegPolyHMM(*this);
  TwoCellFixLib3TrParam2DegPolyHMM* bestGuessHMM = this;//new TwoCellFixLib3TrParam2DegPolyHMM(*this);
  gsl_vector* initGuessAsParams = gsl_vector_alloc(initGuess->size);
  this->convertProbToParam(initGuessAsParams, initGuess);
  Optimizable::bfgs(initGuessAsParams, bestGuessHMM, maxIters, verbose, debug);
  gsl_vector_free(initGuessAsParams);
  return bestGuessHMM;
}

