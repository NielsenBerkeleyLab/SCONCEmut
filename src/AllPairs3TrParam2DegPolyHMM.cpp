#include "AllPairs3TrParam2DegPolyHMM.hpp"

// ctors and destructor
AllPairs3TrParam2DegPolyHMM::AllPairs3TrParam2DegPolyHMM(std::vector<DepthPair*>* depths, std::vector<std::string>* sampleList, int maxPloidy) : AllPairs3TrParam2DegPolyHMM(depths, sampleList, nullptr, maxPloidy, 2, 1, 0, boost::math::binomial_coefficient<double>(depths->size(), 2), boost::math::binomial_coefficient<double>(depths->size(), 2) * 3) { // 2 shared transition params (beta/lambda), 1 fixed transition param (alpha), 0 fixed libs, n choose 2 pairs, n choose 2 pairs * 3 branches
  // delegate everything into the other ctor. Separating them is purely to be able to set NUM_SHARED_TRANSITION_PARAMS_TO_EST and NUM_PAIRS to a custom value elsewhere
}
AllPairs3TrParam2DegPolyHMM* AllPairs3TrParam2DegPolyHMM::create(std::vector<DepthPair*>* depths, std::vector<std::string>* sampleList, int maxPloidy, gsl_vector* meanVarianceCoefVec, bool preallocIntermediates) {
  AllPairs3TrParam2DegPolyHMM* hmm = new AllPairs3TrParam2DegPolyHMM(depths, sampleList, maxPloidy);
  hmm->makeHMMPairs(meanVarianceCoefVec, preallocIntermediates);
  return hmm;
}

AllPairs3TrParam2DegPolyHMM::AllPairs3TrParam2DegPolyHMM(std::vector<DepthPair*>* depths, std::vector<std::string>* sampleList, gsl_vector* fixedParams,  int maxPloidy, int numSharedTrParamsToEst, int numFixedTrParams, int numFixedLibs, int numPairs, int numBranchesToEst) : AllCells3TrParam2DegPolyHMM(depths, sampleList, fixedParams, maxPloidy, numSharedTrParamsToEst, numFixedTrParams, numFixedLibs, numPairs, numBranchesToEst), // numPairs is the numHMMs
                                        NUM_PAIRS(numPairs) {
  this->margLikelihoodMatAcrossPairsVec = nullptr;
  this->summaryPathAcrossPairsVec = nullptr;
  this->summaryMethodUsed = -1;
  this->hmmsCellAppearsInVec = new std::vector<std::set<int>*>(this->NUM_CELLS);
  for(unsigned int i = 0; i < this->hmmsCellAppearsInVec->size(); i++) {
    (*this->hmmsCellAppearsInVec)[i] = new std::set<int>();
  }
}

/*
 * helper method to prep for HMM set up
 */
void AllPairs3TrParam2DegPolyHMM::makeHMMPairs(gsl_vector* meanVarianceCoefVec, bool preallocIntermediates) {
  // prep for HMM set up
  //gsl_vector* meanVarianceCoefVec = HMM::createMeanVarianceCoefVec();
  if(meanVarianceCoefVec == nullptr) {
    meanVarianceCoefVec = HMM::createMeanVarianceCoefVec();
  }
  if(this->meanVarianceCoefVec == nullptr) {
    this->meanVarianceCoefVec = gsl_vector_alloc(meanVarianceCoefVec->size);
  }
  gsl_vector_memcpy(this->meanVarianceCoefVec, meanVarianceCoefVec);

  // prep transition mat
  gsl_vector* transitionParams = gsl_vector_alloc(5);
  //gsl_vector_set(transitionParams, 0, 0.025); // alpha
  //gsl_vector_set(transitionParams, 1, 0.002); // beta
  //gsl_vector_set(transitionParams, 2, 0.04); // gamma
  //gsl_vector_set(transitionParams, 3, 0.05);  // t2
  //gsl_vector_set(transitionParams, 4, 0.05);  // t3
  //gsl_vector_set(transitionParams, 0, 0.025); // alpha
  gsl_vector_set(transitionParams, 0, 0.002); // beta
  gsl_vector_set(transitionParams, 1, 0.04); // lambda
  gsl_vector_set(transitionParams, 2, 0.05);  // t1
  gsl_vector_set(transitionParams, 3, 0.05);  // t2
  gsl_vector_set(transitionParams, 4, 0.05);  // t3

//make genSimulations && ./genSimluations -o input/simulations/simulated_k5_set2 -n 50 -k 5 -b 0.01 -g 0.02 --t1 0.1 --t2 0.2 --t3 0.3 --lib1 1 --lib2 1.5 -s 5
  /*gsl_vector_set(transitionParams, 0, 0.01); // beta
  gsl_vector_set(transitionParams, 1, 0.02); // gamma
  gsl_vector_set(transitionParams, 2, 0.1);  // t1
  gsl_vector_set(transitionParams, 3, 0.2);  // t2
  gsl_vector_set(transitionParams, 4, 0.3);  // t3*/

  // save into paramsToEst
  gsl_vector_set(this->paramsToEst, this->SHARED_TRANSITION_PROB_START_IDX + 0, gsl_vector_get(transitionParams, 0)); // beta
  gsl_vector_set(this->paramsToEst, this->SHARED_TRANSITION_PROB_START_IDX + 1, gsl_vector_get(transitionParams, 1)); // lambda
  for(int i = 0; i < this->NUM_PAIRS; i++) {
    gsl_vector_set(this->paramsToEst, this->BRANCH_LENGTH_START_IDX + 0 + i * 3, gsl_vector_get(transitionParams, 2)); // t1
    gsl_vector_set(this->paramsToEst, this->BRANCH_LENGTH_START_IDX + 1 + i * 3, gsl_vector_get(transitionParams, 3)); // t2
    gsl_vector_set(this->paramsToEst, this->BRANCH_LENGTH_START_IDX + 2 + i * 3, gsl_vector_get(transitionParams, 4)); // t3
  }

  // process all pairs (call subclass's makeHMMPairs)
  makeHMMPairs(meanVarianceCoefVec, transitionParams, preallocIntermediates);
  //printColVector(this->paramsToEst);
  //gsl_vector_free(meanVarianceCoefVec);
  gsl_vector_free(transitionParams);

  //std::cout << "IN AllPairs3TrParam2DegPolyHMM::makeHMMPairs, about to call miscFunctions" << std::endl;
  this->miscFunctions();
  //exit(0);

  this->setHMMNames();
}

void AllPairs3TrParam2DegPolyHMM::makeHMMPairs(gsl_vector* meanVarianceCoefVec, gsl_vector* transitionParams, bool preallocIntermediates) {
  //std::cout << "in AllPairs3TrParam2DegPolyHMM::makeHMMPairs" << std::endl;
  std::vector<DepthPair*>* currDepths = nullptr;
  gsl_vector* currMeanVarCoefVec = nullptr; // make them all have their own copies of this vector
  int hmmIdx = 0;
  for(unsigned int i = 0; i < this->depthsVec->size() && hmmIdx < this->NUM_PAIRS; i++) {
    for(unsigned int j = i+1; j < this->depthsVec->size() && hmmIdx < this->NUM_PAIRS; j++) {
      currDepths = new std::vector<DepthPair*>();
      currDepths->push_back((*this->depthsVec)[i]);
      currDepths->push_back((*this->depthsVec)[j]);
      TwoCell3TrParam2DegPolyHMM* hmm = new TwoCell3TrParam2DegPolyHMM(currDepths, this->getKploidy(), preallocIntermediates);
      (*this->hmmVec)[hmmIdx] = hmm;

      // do rest of HMM set up (usually happens in main.cpp)
      currMeanVarCoefVec = gsl_vector_alloc(meanVarianceCoefVec->size);
      gsl_vector_memcpy(currMeanVarCoefVec, meanVarianceCoefVec);
      hmm->setMeanVarianceFn(currMeanVarCoefVec);
      hmm->setTransition(transitionParams); // this vector isn't saved anywhere
      hmm->setLibScalingFactorsToTotalRatio();
      hmm->setAlpha(this->getAlpha());

      this->setLibScalingFactors(i, j, hmm->getLibScalingFactor(0), hmm->getLibScalingFactor(1));
      this->storeHMMIdxForCells(i, j, hmmIdx);
      hmmIdx++;
    }
  }
}

/*
 * using this->sampleList, set into this->hmmNames, depending on
 * subclass's getCell0IdxFromHMMIdx and getCell1IdxFromHMMIdx.
 * Needs to be its own method instead of in the ctor so that the subclass's
 * getCell0IdxFromHMMIdx and getCell1IdxFromHMMIdx methods can be accessed
 */
void AllPairs3TrParam2DegPolyHMM::setHMMNames() {
  this->hmmNames = new std::vector<std::string>(this->NUM_HMMS);
  int cell0Idx = -1;
  int cell1Idx = -1;
  for(int hmmIdx = 0; hmmIdx < this->NUM_HMMS; hmmIdx++) {
    cell0Idx = this->getCell0IdxFromHMMIdx(hmmIdx);
    cell1Idx = this->getCell1IdxFromHMMIdx(hmmIdx);
    std::string cell0Name = (*this->sampleList)[cell0Idx];
    std::string cell1Name = (*this->sampleList)[cell1Idx];
    (*this->hmmNames)[hmmIdx] = cell0Name + "," + cell1Name;
  }
}

AllPairs3TrParam2DegPolyHMM::~AllPairs3TrParam2DegPolyHMM() {
  // TODO
}

/*
 * Set of helper methods to convert between tumor cell 0 and cell 1 (DepthPair/lib size scaling) indices into the linear HMM vector index.
 * Because we're using all pairs of DepthPairs, they can be thought of as an upper triangular matrix
 * and the HMM vec is the linear representation of the pairs.
 * Calculations from https://stackoverflow.com/a/27088560
 * where k = linear index, n = number of cells, i = cell0, j = cell1
 *
 * from linear index to (i, j):
 *   i = n - 2 - floor(sqrt(-8*k + 4*n*(n-1)-7)/2.0 - 0.5)
 *   j = k + i + 1 - n*(n-1)/2 + (n-i)*((n-i)-1)/2
 *
 * from (i, j) to linear index:
 *   k = (n*(n-1)/2) - (n-i)*((n-i)-1)/2 + j - i - 1
 */
int AllPairs3TrParam2DegPolyHMM::getCell0IdxFromHMMIdx(int hmmIdx) {
  //std::cout <<" &&&&&&&&&&&&&&&&&& in AllPairs3TrParam2DegPolyHMM::getCell0IdxFromHMMIdx" << std::endl;
  return (int) (this->NUM_CELLS - 2 - floor(sqrt(-8 * hmmIdx + 4 * this->NUM_CELLS *(this->NUM_CELLS - 1) - 7) / 2.0 - 0.5));
}
int AllPairs3TrParam2DegPolyHMM::getCell1IdxFromHMMIdx(int hmmIdx) {
  int i = getCell0IdxFromHMMIdx(hmmIdx);
  return (int) (hmmIdx + i + 1 - this->NUM_CELLS * (this->NUM_CELLS-1) / 2 + (this->NUM_CELLS - i) * ((this->NUM_CELLS - i) - 1) / 2);
}
int AllPairs3TrParam2DegPolyHMM::getHMMIdxFromCellPair(int cell0Idx, int cell1Idx) {
  return (int) ((this->NUM_CELLS * (this->NUM_CELLS - 1) / 2) - (this->NUM_CELLS - cell0Idx) * ((this->NUM_CELLS - cell0Idx) - 1) / 2 + cell1Idx - cell0Idx - 1);
}

/*
 * function to return a vector of all the HMMs this particular cell appears in as cell 0 (left cell in pair)
 * note: does not return all possible HMMs for this cell, only the ones that this cell actually appears in (ie in case not all possible pairs were created)
 */
std::vector<HMM*>* AllPairs3TrParam2DegPolyHMM::getHMMsWithCellAsCell0(int cellIdx) {
  std::vector<HMM*>* hmmsWithCell = new std::vector<HMM*>();
  //std::cout << "cellIdx " << cellIdx << " is cell0 in: ";

  int cell0Idx = -1;
  int cell1Idx = -1;
  int hmmIdx = -1;
  HMM* currHMM = nullptr;

  // set cell0Idx to be cellIdx (ie cell of interest is cell0/left in pair)
  cell0Idx = cellIdx;

  // for cell1Idx in 0:numCells
  for(cell1Idx = 0; cell1Idx < this->NUM_CELLS; cell1Idx++) {
    // left cell should always be less than right cell
    if(cell1Idx <= cell0Idx) {
      continue;
    }

    // getHMMidx(cell0, cell1)
    hmmIdx = this->getHMMIdxFromCellPair(cell0Idx, cell1Idx);

    // if proposed pair doesn't yield a valid pair
    if(hmmIdx == -1) {
      continue;
    }

    //std::cout << hmmIdx << " ";
    // filter for which hmms this cell actually appears in (ie in the SelectPairs classes, cells don't necessarily appear in all possible hmms)
    if((*this->hmmsCellAppearsInVec)[cellIdx]->count(hmmIdx) == 0) {
      //std::cout << "but doesn't appear in hmmsCellAppearsInVec, skipping" << std::endl;
      continue;
    }

    // store pointer to hmmIdx
    currHMM = (*this->hmmVec)[hmmIdx];
    // nullptr guard
    if(currHMM != nullptr) {
      hmmsWithCell->push_back(currHMM);
    }
  }
  //std::cout << std::endl;
  return hmmsWithCell;
}

/*
 * function to return a vector of all the HMMs this particular cell appears in as cell 1 (right cell in pair)
 * note: does not return all possible HMMs for this cell, only the ones that this cell actually appears in (ie in case not all possible pairs were created)
 */
std::vector<HMM*>* AllPairs3TrParam2DegPolyHMM::getHMMsWithCellAsCell1(int cellIdx) {
  std::vector<HMM*>* hmmsWithCell = new std::vector<HMM*>();
  //std::cout << "cellIdx " << cellIdx << " is cell1 in: ";

  int cell0Idx = -1;
  int cell1Idx = -1;
  int hmmIdx = -1;
  HMM* currHMM = nullptr;

  // set cell1Idx to be cellIdx
  cell1Idx = cellIdx;
  // for cell0Idx in 1:numCells
  for(cell0Idx = 0; cell0Idx < this->NUM_CELLS; cell0Idx++) {
    if(cell1Idx <= cell0Idx) {
      continue;
    }
    hmmIdx = this->getHMMIdxFromCellPair(cell0Idx, cell1Idx);

    if(hmmIdx == -1) {
      continue;
    }
    //std::cout << hmmIdx << " ";

    if((*this->hmmsCellAppearsInVec)[cellIdx]->count(hmmIdx) == 0) {
      //std::cout << "but doesn't appear in hmmsCellAppearsInVec, skipping" << std::endl;
      continue;
    }

    currHMM = (*this->hmmVec)[hmmIdx];
    // nullptr guard
    if(currHMM != nullptr) {
      hmmsWithCell->push_back(currHMM);
    }
  }
  //std::cout << std::endl;
  return hmmsWithCell;
}

/*
 * function to add hmmIdx to this->hmmsCellAppearsInVec for both cell0Idx and cell1Idx
 */
void AllPairs3TrParam2DegPolyHMM::storeHMMIdxForCells(int cell0Idx, int cell1Idx, int hmmIdx) {
  //std::cout << "storing hmmIdx: " << hmmIdx << " for cell0: " << cell0Idx << ", and cell1: " << cell1Idx << std::endl;
  this->storeHMMIdxForCell(cell0Idx, hmmIdx);
  this->storeHMMIdxForCell(cell1Idx, hmmIdx);
}
/*
 * function to store hmmIdx for just one cell (ie to  remember directionality). useful for if you want to summarize over only some cells
 */
// only one cell
void AllPairs3TrParam2DegPolyHMM::storeHMMIdxForCell(int cellIdx, int hmmIdx) {
  //std::cout << "storing hmmIdx: " << hmmIdx << " for cell0: " << cell0Idx <<  std::endl;
  ((*this->hmmsCellAppearsInVec)[cellIdx])->insert(hmmIdx);
}

void AllPairs3TrParam2DegPolyHMM::setLibScalingFactors(int cell0Idx, int cell1Idx, double lib0, double lib1) {
  /*gsl_vector_set(this->paramsToEst, this->LIB_SIZE_SCALING_FACTOR_START_IDX + cellNum, libScalingFactor);

  // cell 0
  if(cellNum % 2 == 0) {
    int hmmIdx = getHMMIdxFromCellPair(cellNum, cellNum+1);
    (*this->hmmVec)[hmmIdx]->setLibScalingFactor(0, libScalingFactor);
  }
  // cell 1
  else {
    int hmmIdx = getHMMIdxFromCellPair(cellNum-1, cellNum);
    (*this->hmmVec)[hmmIdx]->setLibScalingFactor(1, libScalingFactor);
  }*/
  gsl_vector_set(this->paramsToEst, this->LIB_SIZE_SCALING_FACTOR_START_IDX + cell0Idx, lib0);
  gsl_vector_set(this->paramsToEst, this->LIB_SIZE_SCALING_FACTOR_START_IDX + cell1Idx, lib1);
  int hmmIdx = getHMMIdxFromCellPair(cell0Idx, cell1Idx);
  (*this->hmmVec)[hmmIdx]->setLibScalingFactor(0, lib0);
  (*this->hmmVec)[hmmIdx]->setLibScalingFactor(1, lib1);
}

/*
 * sets all cells in indicated pair to passed libScalingFactor. Sets in this class, as well as each individual HMM
 * cellNumInPair should be either 0 or 1 (assumed, no checks)
 */
void AllPairs3TrParam2DegPolyHMM::setAllLibScalingFactors(int cellNumInPair, double libScalingFactor) {
  for(unsigned int i = 0; i < this->hmmVec->size(); i++) {
    //this->setLibScalingFactor(i * 2 + cellNumInPair, libScalingFactor);
    gsl_vector_set(this->paramsToEst, this->LIB_SIZE_SCALING_FACTOR_START_IDX + i * 2 + cellNumInPair, libScalingFactor);
    (*this->hmmVec)[i]->setLibScalingFactor(cellNumInPair, libScalingFactor);
  }
}
double AllPairs3TrParam2DegPolyHMM::getLibScalingFactor(int cellNum) const {
  return gsl_vector_get(this->paramsToEst, this->LIB_SIZE_SCALING_FACTOR_START_IDX + cellNum);
}

/*
 * sets transition matrices for all HMMs to be these shared transitionParams.
 * Also saves into this->paramsToEst. Useful for initializing all HMMs at the same transition params
 * (setParamsToEst requires individual branch lengths, libs, etc).
 * //transitionParams should be [beta, gamma, t1, t2, t3] (ie shared for all HMMs)
 * transitionParams should be [beta, lambda, t1, t2, t3] (ie shared for all HMMs)
 */
double AllPairs3TrParam2DegPolyHMM::setAllTransition(gsl_vector* transitionParams) {
  // save beta and lambda
  double beta = gsl_vector_get(transitionParams, this->SHARED_TRANSITION_PROB_START_IDX - this->NUM_LIBS_TO_EST + 0); // idx should be 0
  //double gamma = gsl_vector_get(transitionParams, this->SHARED_TRANSITION_PROB_START_IDX - this->NUM_LIBS_TO_EST + 1);
  double lambda = gsl_vector_get(transitionParams, this->SHARED_TRANSITION_PROB_START_IDX - this->NUM_LIBS_TO_EST + 1);
  gsl_vector_set(this->paramsToEst, this->SHARED_TRANSITION_PROB_START_IDX + 0, beta);
  //gsl_vector_set(this->paramsToEst, this->SHARED_TRANSITION_PROB_START_IDX + 1, gamma);
  gsl_vector_set(this->paramsToEst, this->SHARED_TRANSITION_PROB_START_IDX + 1, lambda);

  double t1 = gsl_vector_get(transitionParams, this->BRANCH_LENGTH_START_IDX - this->NUM_LIBS_TO_EST + 0); // idx should be 2
  double t2 = gsl_vector_get(transitionParams, this->BRANCH_LENGTH_START_IDX - this->NUM_LIBS_TO_EST + 1);
  double t3 = gsl_vector_get(transitionParams, this->BRANCH_LENGTH_START_IDX - this->NUM_LIBS_TO_EST + 2);

  double status = 0;
  for(unsigned int hmmIdx = 0; hmmIdx < this->hmmVec->size(); hmmIdx++) {
    // nullptr guard
    if((*this->hmmVec)[hmmIdx] == nullptr) {
      continue;
    }
    gsl_vector_set(this->paramsToEst, this->BRANCH_LENGTH_START_IDX + hmmIdx * 3 + 0, t1);
    gsl_vector_set(this->paramsToEst, this->BRANCH_LENGTH_START_IDX + hmmIdx * 3 + 1, t2);
    gsl_vector_set(this->paramsToEst, this->BRANCH_LENGTH_START_IDX + hmmIdx * 3 + 2, t3);

    status += (*this->hmmVec)[hmmIdx]->setTransition(transitionParams);
  }
  return status;
}

/*
 * given params from hmmIdx'th HMM, save into this class's paramsToEst. Does not save into the hmm itself.
 * assumes params = [lib0, lib1, beta, lambda, t1, t2, t3]
 */
void AllPairs3TrParam2DegPolyHMM::setParamsToEstFromIthHMM(gsl_vector* params, int hmmIdx) {
  gsl_vector_set(this->paramsToEst, this->LIB_SIZE_SCALING_FACTOR_START_IDX + hmmIdx + 0, gsl_vector_get(params, 0)); // lib0
  gsl_vector_set(this->paramsToEst, this->LIB_SIZE_SCALING_FACTOR_START_IDX + hmmIdx + 1, gsl_vector_get(params, 1)); // lib1
  gsl_vector_set(this->paramsToEst, this->SHARED_TRANSITION_PROB_START_IDX + 0, gsl_vector_get(params, 2)); // beta
  gsl_vector_set(this->paramsToEst, this->SHARED_TRANSITION_PROB_START_IDX + 1, gsl_vector_get(params, 3)); // lambda
  gsl_vector_set(this->paramsToEst, this->BRANCH_LENGTH_START_IDX + hmmIdx + 0, gsl_vector_get(params, 4)); // t1
  gsl_vector_set(this->paramsToEst, this->BRANCH_LENGTH_START_IDX + hmmIdx + 1, gsl_vector_get(params, 5)); // t2
  gsl_vector_set(this->paramsToEst, this->BRANCH_LENGTH_START_IDX + hmmIdx + 2, gsl_vector_get(params, 6)); // t3
}
/*
 * given params from hmmIdx'th HMM, save into this class's fixedParams. Does not save into the hmm itself.
 * assumes params = [alpha]
 */
void AllPairs3TrParam2DegPolyHMM::setFixedParamsFromIthHMM(gsl_vector* params, int hmmIdx) {
  gsl_vector_set(this->fixedParams, this->FIXED_TRANSITION_PROB_START_IDX + 0, gsl_vector_get(params, 0)); // alpha
}



// Optimizable methods
/*
 * method to extract HMM specific variables and set those variables for each HMM,
 * using each HMM's own setParamsToEst method.
 * params is in probability space
 */
double AllPairs3TrParam2DegPolyHMM::setParamsToEst(gsl_vector* params) {
  gsl_vector_memcpy(this->paramsToEst, params); // Thu 30 Jan 2020 10:47:35 AM PST changed to be a memcpy instead of free/alloc cycle
  //std::cout << "AllPairs3TrParam2DegPolyHMM setParamsToEst params: " << std::endl;
  //printColVector(params);

  // nullptr guard
  HMM* hmm = this->getFirstNonNullHMM();
  int hmmLibIdx = hmm->LIB_SIZE_SCALING_FACTOR_START_IDX;
  int hmmTrIdx = hmm->TRANSITION_PROB_START_IDX;
  int hmmBranchIdx = hmm->BRANCH_LENGTH_START_IDX;
  int numHMMParams = hmm->getNumParamsToEst();

  gsl_vector* currHMMParams = gsl_vector_alloc(numHMMParams);
  //std::cout << hmmLibIdx << ", " << hmmTrIdx << ", " << hmmBranchIdx << ", " << currHMMParams->size << std::endl;
  //std::cout << this->LIB_SIZE_SCALING_FACTOR_START_IDX << ", " << this->SHARED_TRANSITION_PROB_START_IDX << ", " << this->BRANCH_LENGTH_START_IDX << std::endl;

  int cell0Idx = -1;
  int cell1Idx = -1;

  double currLib0BFGS = 0;
  double currLib1BFGS = 0;
  //double alphaBFGS = gsl_vector_get(params, this->SHARED_TRANSITION_PROB_START_IDX + 0);
  //double betaBFGS  = gsl_vector_get(params, this->SHARED_TRANSITION_PROB_START_IDX + 1);
  //double gammaBFGS = gsl_vector_get(params, this->SHARED_TRANSITION_PROB_START_IDX + 2);
  double betaBFGS  = gsl_vector_get(params, this->SHARED_TRANSITION_PROB_START_IDX + 0);
  //double gammaBFGS = gsl_vector_get(params, this->SHARED_TRANSITION_PROB_START_IDX + 1);
  double lambdaBFGS = gsl_vector_get(params, this->SHARED_TRANSITION_PROB_START_IDX + 1);
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
    //// alloc new currHMMParams for each HMM (old paramsToEst for each HMM is freed in HMM::setParamsToEst) Thu 30 Jan 2020 10:50:38 AM PST changed to be memcpy
    //currHMMParams = gsl_vector_alloc(numHMMParams);

    // get cell0 and cell1 indices for libs
    cell0Idx = this->getCell0IdxFromHMMIdx(hmmIdx);
    cell1Idx = this->getCell1IdxFromHMMIdx(hmmIdx);

    // get appropriate libs
    currLib0BFGS = gsl_vector_get(params, this->LIB_SIZE_SCALING_FACTOR_START_IDX + cell0Idx);
    currLib1BFGS = gsl_vector_get(params, this->LIB_SIZE_SCALING_FACTOR_START_IDX + cell1Idx);

    //// get apropriate branch lengths (there are 2 branch lengths stored per HMM
    //currT2BFGS = gsl_vector_get(params, this->BRANCH_LENGTH_START_IDX + 2 * hmmIdx + 0);
    //currT3BFGS = gsl_vector_get(params, this->BRANCH_LENGTH_START_IDX + 2 * hmmIdx + 1);
    // get apropriate branch lengths (there are 3 branch lengths stored per HMM
    currT1BFGS = gsl_vector_get(params, this->BRANCH_LENGTH_START_IDX + 3 * hmmIdx + 0);
    currT2BFGS = gsl_vector_get(params, this->BRANCH_LENGTH_START_IDX + 3 * hmmIdx + 1);
    currT3BFGS = gsl_vector_get(params, this->BRANCH_LENGTH_START_IDX + 3 * hmmIdx + 2);

    //std::cout << currLib0BFGS << ", " << currLib1BFGS << ", " << alphaBFGS << ", " << betaBFGS << ", " << gammaBFGS << ", " << currT2BFGS << ", " << currT3BFGS << std::endl;
    // set everything into currHMMParams
    gsl_vector_set(currHMMParams, hmmLibIdx + 0, currLib0BFGS);
    gsl_vector_set(currHMMParams, hmmLibIdx + 1, currLib1BFGS);
    //gsl_vector_set(currHMMParams, hmmTrIdx + 0, alphaBFGS);
    //gsl_vector_set(currHMMParams, hmmTrIdx + 1, betaBFGS);
    //gsl_vector_set(currHMMParams, hmmTrIdx + 2, gammaBFGS);
    gsl_vector_set(currHMMParams, hmmTrIdx + 0, betaBFGS);
    //gsl_vector_set(currHMMParams, hmmTrIdx + 1, gammaBFGS);
    gsl_vector_set(currHMMParams, hmmTrIdx + 1, lambdaBFGS);
    //gsl_vector_set(currHMMParams, hmmBranchIdx + 0, currT2BFGS);
    //gsl_vector_set(currHMMParams, hmmBranchIdx + 1, currT3BFGS);
    gsl_vector_set(currHMMParams, hmmBranchIdx + 0, currT1BFGS);
    gsl_vector_set(currHMMParams, hmmBranchIdx + 1, currT2BFGS);
    gsl_vector_set(currHMMParams, hmmBranchIdx + 2, currT3BFGS);

    // call setParamsToEst on subclass, summing return status for each one (GSL_SUCCESS = 0, as returned by HMM::findSteadyStateDist)
    //std::cout << "IN subclass::SETPARAMSTOEST, currHMMParams:" << std::endl;
    //printColVector(currHMMParams);
    status += (*this->hmmVec)[hmmIdx]->setParamsToEst(currHMMParams);
    //(*this->hmmVec)[hmmIdx]->print(stdout);
  }
  //gsl_vector_free(currHMMParams); // each HMM needs its own copy; don't free this here
  //std::cout << "end of subclass:setParamsToEst" << std::endl;
  //this->print(stdout);
  //std::cout << "#########" << std::endl;
  return status;

}

/*
 * convert probabilitiy space values in src to BFGS space values in dest.
 * See TwoCell3TrParam2DegPolyHMM versions and comments for constraints/calculations;
 * these are the same, just scaled up
 */
void AllPairs3TrParam2DegPolyHMM::convertProbToParam(gsl_vector* dest, const gsl_vector* src) const {
  //std::cout << "AllPairs3TrParam2DegPolyHMM::convertProbToParam" << std::endl;
  /*gsl_vector_memcpy(dest, src); // debugging penalty
  return;*/
  // lib scaling factors
  double r = 0;
  for(int cellIdx = 0; cellIdx < this->NUM_LIBS_TO_EST; cellIdx++) {
    r = gsl_vector_get(src, this->LIB_SIZE_SCALING_FACTOR_START_IDX + cellIdx);
    gsl_vector_set(dest, this->LIB_SIZE_SCALING_FACTOR_START_IDX + cellIdx, log(r));
  }

  //// shared transition parameters
  ////double a = gsl_vector_get(src, this->SHARED_TRANSITION_PROB_START_IDX + 0);
  ////double b = gsl_vector_get(src, this->SHARED_TRANSITION_PROB_START_IDX + 1);
  ////double g = gsl_vector_get(src, this->SHARED_TRANSITION_PROB_START_IDX + 2);
  ////double c = 1.0 / (1.0 - 2.0 * a - (double) this->getKploidy() * b - g);
  ////gsl_vector_set(dest, this->SHARED_TRANSITION_PROB_START_IDX + 0, log(a * 2.0 * c));
  ////gsl_vector_set(dest, this->SHARED_TRANSITION_PROB_START_IDX + 1, log(b * (double) this->getKploidy() * c));
  ////gsl_vector_set(dest, this->SHARED_TRANSITION_PROB_START_IDX + 2, log(g * c));
  //double a = this->getAlpha(); //this->alpha;
  //double b = gsl_vector_get(src, this->SHARED_TRANSITION_PROB_START_IDX + 0);
  //double g = gsl_vector_get(src, this->SHARED_TRANSITION_PROB_START_IDX + 1);
  //double c = (b * this->getKploidy() + g - 1);
  //gsl_vector_set(dest, this->SHARED_TRANSITION_PROB_START_IDX + 0, log(-b * (1-2*a) * (double) this->getKploidy() / c)); // set y
  //gsl_vector_set(dest, this->SHARED_TRANSITION_PROB_START_IDX + 1, log(-g * (1-2*a) / c)); // set z
  ////std::cout << "ALPHA: " << a << ", BETA: " << b << ", GAMMA: " << g << ", b transformed: " << gsl_vector_get(dest, this->SHARED_TRANSITION_PROB_START_IDX + 0) << ", g transformed: " << gsl_vector_get(dest, this->SHARED_TRANSITION_PROB_START_IDX + 1) << std::endl;

  //// pairwise branch lengths
  //double d = (double) (*this->depthsVec)[0]->maxWindowSize;
  //double t1 = 0;
  //double t2 = 0;
  //double t3 = 0;
  //for(int pairIdx = 0; pairIdx < this->NUM_PAIRS; pairIdx++) {
  //  //t2 = gsl_vector_get(src, this->BRANCH_LENGTH_START_IDX + 2 * pairIdx + 0) / d;
  //  //t3 = gsl_vector_get(src, this->BRANCH_LENGTH_START_IDX + 2 * pairIdx + 1) / d;
  //  //gsl_vector_set(dest, this->BRANCH_LENGTH_START_IDX + 2 * pairIdx + 0, log(-(d * t2) / (d * (t2 + t3) - 1))); // set T2
  //  //gsl_vector_set(dest, this->BRANCH_LENGTH_START_IDX + 2 * pairIdx + 1, log(-(d * t3) / (d * (t2 + t3) - 1))); // set T3

  //  t1 = gsl_vector_get(src, this->BRANCH_LENGTH_START_IDX + 3 * pairIdx + 0) / d;
  //  t2 = gsl_vector_get(src, this->BRANCH_LENGTH_START_IDX + 3 * pairIdx + 1) / d;
  //  t3 = gsl_vector_get(src, this->BRANCH_LENGTH_START_IDX + 3 * pairIdx + 2) / d;
  //  c = (d * (t1 + t2 + t3) - 1);
  //  gsl_vector_set(dest, this->BRANCH_LENGTH_START_IDX + 3 * pairIdx + 0, log(-(d * t1) / c)); // set T1
  //  gsl_vector_set(dest, this->BRANCH_LENGTH_START_IDX + 3 * pairIdx + 1, log(-(d * t2) / c)); // set T2
  //  gsl_vector_set(dest, this->BRANCH_LENGTH_START_IDX + 3 * pairIdx + 2, log(-(d * t3) / c)); // set T3
  //}

  // shared transition parameters
  double beta = gsl_vector_get(src, this->SHARED_TRANSITION_PROB_START_IDX + 0);
  double lambda = gsl_vector_get(src, this->SHARED_TRANSITION_PROB_START_IDX + 1);
  gsl_vector_set(dest, this->SHARED_TRANSITION_PROB_START_IDX + 0, log(beta)); // set y
  gsl_vector_set(dest, this->SHARED_TRANSITION_PROB_START_IDX + 1, log(lambda)); // set z

  // pairwise branch lengths
  double t1 = 0;
  double t2 = 0;
  double t3 = 0;
  for(int pairIdx = 0; pairIdx < this->NUM_PAIRS; pairIdx++) {
    t1 = gsl_vector_get(src, this->BRANCH_LENGTH_START_IDX + 3 * pairIdx + 0);
    t2 = gsl_vector_get(src, this->BRANCH_LENGTH_START_IDX + 3 * pairIdx + 1);
    t3 = gsl_vector_get(src, this->BRANCH_LENGTH_START_IDX + 3 * pairIdx + 2);
    gsl_vector_set(dest, this->BRANCH_LENGTH_START_IDX + 3 * pairIdx + 0, log(t1)); // set T1
    gsl_vector_set(dest, this->BRANCH_LENGTH_START_IDX + 3 * pairIdx + 1, log(t2)); // set T2
    gsl_vector_set(dest, this->BRANCH_LENGTH_START_IDX + 3 * pairIdx + 2, log(t3)); // set T3
  }
}

void AllPairs3TrParam2DegPolyHMM::convertParamToProb(gsl_vector* dest, const gsl_vector* src) const {
  //std::cout << "AllPairs3TrParam2DegPolyHMM::convertParamToProb" << std::endl;
  /*gsl_vector_memcpy(dest, src); // debugging penalty
  return;*/
  // lib scaling factors
  double w = 0;
  for(int cellIdx = 0; cellIdx < this->NUM_LIBS_TO_EST; cellIdx++) {
    w = gsl_vector_get(src, this->LIB_SIZE_SCALING_FACTOR_START_IDX + cellIdx);
    gsl_vector_set(dest, this->LIB_SIZE_SCALING_FACTOR_START_IDX + cellIdx, exp(w));
  }

  //// shared transition parameters
  ////double x = gsl_vector_get(src, this->SHARED_TRANSITION_PROB_START_IDX + 0);
  ////double y = gsl_vector_get(src, this->SHARED_TRANSITION_PROB_START_IDX + 1);
  ////double z = gsl_vector_get(src, this->SHARED_TRANSITION_PROB_START_IDX + 2);
  ////double c = 1 + exp(x) + exp(y) + exp(z);
  ////gsl_vector_set(dest, this->SHARED_TRANSITION_PROB_START_IDX + 0, exp(x) / (2.0 * c));
  ////gsl_vector_set(dest, this->SHARED_TRANSITION_PROB_START_IDX + 1, exp(y) / ((double) this->getKploidy() * c));
  ////gsl_vector_set(dest, this->SHARED_TRANSITION_PROB_START_IDX + 2, exp(z) / c);
  //double a = this->getAlpha(); //this->alpha;
  //double y = gsl_vector_get(src, this->SHARED_TRANSITION_PROB_START_IDX + 0);
  //double z = gsl_vector_get(src, this->SHARED_TRANSITION_PROB_START_IDX + 1);
  //double c = 1 - 2.0*a + exp(y) + exp(z);
  //gsl_vector_set(dest, this->SHARED_TRANSITION_PROB_START_IDX + 0, exp(y) / ((double) this->getKploidy() * c)); // beta
  //gsl_vector_set(dest, this->SHARED_TRANSITION_PROB_START_IDX + 1, exp(z) / c); // gamma

  //// pairwise branch lengths
  //double T1 = 0;
  //double T2 = 0;
  //double T3 = 0;
  //for(int pairIdx = 0; pairIdx < this->NUM_PAIRS; pairIdx++) {
  //  //T2 = gsl_vector_get(src, this->BRANCH_LENGTH_START_IDX + 2 * pairIdx + 0);
  //  //T3 = gsl_vector_get(src, this->BRANCH_LENGTH_START_IDX + 2 * pairIdx + 1);
  //  //c = 1.0 / (1 + exp(T2) + exp(T3));
  //  //gsl_vector_set(dest, this->BRANCH_LENGTH_START_IDX + 2 * pairIdx + 0, exp(T2) * c); // set t2
  //  //gsl_vector_set(dest, this->BRANCH_LENGTH_START_IDX + 2 * pairIdx + 1, exp(T3) * c); // set t3

  //  T1 = gsl_vector_get(src, this->BRANCH_LENGTH_START_IDX + 3 * pairIdx + 0);
  //  T2 = gsl_vector_get(src, this->BRANCH_LENGTH_START_IDX + 3 * pairIdx + 1);
  //  T3 = gsl_vector_get(src, this->BRANCH_LENGTH_START_IDX + 3 * pairIdx + 2);
  //  c = 1.0 / (1 + exp(T1) + exp(T2) + exp(T3));
  //  gsl_vector_set(dest, this->BRANCH_LENGTH_START_IDX + 3 * pairIdx + 0, exp(T1) * c); // set t1
  //  gsl_vector_set(dest, this->BRANCH_LENGTH_START_IDX + 3 * pairIdx + 1, exp(T2) * c); // set t2
  //  gsl_vector_set(dest, this->BRANCH_LENGTH_START_IDX + 3 * pairIdx + 2, exp(T3) * c); // set t3
  //}

  // shared transition parameters
  double y = gsl_vector_get(src, this->SHARED_TRANSITION_PROB_START_IDX + 0);
  double z = gsl_vector_get(src, this->SHARED_TRANSITION_PROB_START_IDX + 1);
  gsl_vector_set(dest, this->SHARED_TRANSITION_PROB_START_IDX + 0, exp(y)); // beta
  gsl_vector_set(dest, this->SHARED_TRANSITION_PROB_START_IDX + 1, exp(z)); // lambda

  // pairwise branch lengths
  double T1 = 0;
  double T2 = 0;
  double T3 = 0;
  for(int pairIdx = 0; pairIdx < this->NUM_PAIRS; pairIdx++) {
    T1 = gsl_vector_get(src, this->BRANCH_LENGTH_START_IDX + 3 * pairIdx + 0);
    T2 = gsl_vector_get(src, this->BRANCH_LENGTH_START_IDX + 3 * pairIdx + 1);
    T3 = gsl_vector_get(src, this->BRANCH_LENGTH_START_IDX + 3 * pairIdx + 2);
    gsl_vector_set(dest, this->BRANCH_LENGTH_START_IDX + 3 * pairIdx + 0, exp(T1)); // set t1
    gsl_vector_set(dest, this->BRANCH_LENGTH_START_IDX + 3 * pairIdx + 1, exp(T2)); // set t2
    gsl_vector_set(dest, this->BRANCH_LENGTH_START_IDX + 3 * pairIdx + 2, exp(T3)); // set t3
  }
}

/*
 * same as setParamsToEst, but using simParamsToEst
 */
void AllPairs3TrParam2DegPolyHMM::setSimParamsToEst(gsl_vector* params) {
  this->simParamsToEst = gsl_vector_alloc(params->size);
  gsl_vector_memcpy(this->simParamsToEst, params); // Thu 30 Jan 2020 10:47:35 AM PST changed to be a memcpy instead of free/alloc cycle

  // nullptr guard
  HMM* hmm = this->getFirstNonNullHMM();
  int hmmLibIdx = hmm->LIB_SIZE_SCALING_FACTOR_START_IDX;
  int hmmTrIdx = hmm->TRANSITION_PROB_START_IDX;
  int hmmBranchIdx = hmm->BRANCH_LENGTH_START_IDX;
  int numHMMParams = hmm->getNumParamsToEst();
  gsl_vector* currHMMParams = gsl_vector_alloc(numHMMParams);

  int cell0Idx = -1;
  int cell1Idx = -1;

  double currLib0BFGS = 0;
  double currLib1BFGS = 0;
  double betaBFGS  = gsl_vector_get(params, this->SHARED_TRANSITION_PROB_START_IDX + 0);
  //double gammaBFGS = gsl_vector_get(params, this->SHARED_TRANSITION_PROB_START_IDX + 1);
  double lambdaBFGS = gsl_vector_get(params, this->SHARED_TRANSITION_PROB_START_IDX + 1);
  double currT1BFGS = 0;
  double currT2BFGS = 0;
  double currT3BFGS = 0;

  // for each HMM, call that HMM's setSimParamsToEst method
  for(unsigned int hmmIdx = 0; hmmIdx < this->hmmVec->size(); hmmIdx++) {
    // nullptr guard
    if((*this->hmmVec)[hmmIdx] == nullptr) {
      continue;
    }

    // get cell0 and cell1 indices for libs
    cell0Idx = this->getCell0IdxFromHMMIdx(hmmIdx);
    cell1Idx = this->getCell1IdxFromHMMIdx(hmmIdx);

    // get appropriate libs
    if(this->NUM_LIBS_TO_EST > 0) {
      currLib0BFGS = gsl_vector_get(params, this->LIB_SIZE_SCALING_FACTOR_START_IDX + cell0Idx);
      currLib1BFGS = gsl_vector_get(params, this->LIB_SIZE_SCALING_FACTOR_START_IDX + cell1Idx);
      gsl_vector_set(currHMMParams, hmmLibIdx + 0, currLib0BFGS);
      gsl_vector_set(currHMMParams, hmmLibIdx + 1, currLib1BFGS);
    }

    // get apropriate branch lengths (there are 3 branch lengths stored per HMM
    currT1BFGS = gsl_vector_get(params, this->BRANCH_LENGTH_START_IDX + 3 * hmmIdx + 0);
    currT2BFGS = gsl_vector_get(params, this->BRANCH_LENGTH_START_IDX + 3 * hmmIdx + 1);
    currT3BFGS = gsl_vector_get(params, this->BRANCH_LENGTH_START_IDX + 3 * hmmIdx + 2);

    gsl_vector_set(currHMMParams, hmmTrIdx + 0, betaBFGS);
    //gsl_vector_set(currHMMParams, hmmTrIdx + 1, gammaBFGS);
    gsl_vector_set(currHMMParams, hmmTrIdx + 1, lambdaBFGS);
    gsl_vector_set(currHMMParams, hmmBranchIdx + 0, currT1BFGS);
    gsl_vector_set(currHMMParams, hmmBranchIdx + 1, currT2BFGS);
    gsl_vector_set(currHMMParams, hmmBranchIdx + 2, currT3BFGS);

    // call setParamsToEst on subclass, summing return status for each one (GSL_SUCCESS = 0, as returned by HMM::findSteadyStateDist)
    (*this->hmmVec)[hmmIdx]->setSimParamsToEst(currHMMParams);
  }
  gsl_vector_free(currHMMParams);
}
void AllPairs3TrParam2DegPolyHMM::setSimFixedParams(gsl_vector* params) {
  this->simFixedParams = gsl_vector_alloc(params->size);
  gsl_vector_memcpy(this->simFixedParams, params); // Thu 30 Jan 2020 10:47:35 AM PST changed to be a memcpy instead of free/alloc cycle

  // nullptr guard
  HMM* hmm = this->getFirstNonNullHMM();
  int hmmLibIdx = 0;//(*this->hmmVec)[0]->LIB_SIZE_SCALING_FACTOR_START_IDX;
  int hmmTrIdx = hmm->FIXED_TRANSITION_PROB_START_IDX;
  int numHMMParams = hmm->getNumFixedParams();
  gsl_vector* currHMMParams = gsl_vector_alloc(numHMMParams);

  int cell0Idx = -1;
  int cell1Idx = -1;

  double currLib0BFGS = 0;
  double currLib1BFGS = 0;
  double alpha = 0;

  for(unsigned int hmmIdx = 0; hmmIdx < this->hmmVec->size(); hmmIdx++) {
    // nullptr guard
    if((*this->hmmVec)[hmmIdx] == nullptr) {
      continue;
    }

    // get cell0 and cell1 indices for libs
    cell0Idx = this->getCell0IdxFromHMMIdx(hmmIdx);
    cell1Idx = this->getCell1IdxFromHMMIdx(hmmIdx);

    // get appropriate libs
    if(this->NUM_FIXED_LIBS > 0) {
      currLib0BFGS = gsl_vector_get(params, this->LIB_SIZE_SCALING_FACTOR_START_IDX + cell0Idx);
      currLib1BFGS = gsl_vector_get(params, this->LIB_SIZE_SCALING_FACTOR_START_IDX + cell1Idx);
      gsl_vector_set(currHMMParams, hmmLibIdx + 0, currLib0BFGS);
      gsl_vector_set(currHMMParams, hmmLibIdx + 1, currLib1BFGS);
    }
    alpha = gsl_vector_get(params, this->FIXED_TRANSITION_PROB_START_IDX);
    gsl_vector_set(currHMMParams, hmmTrIdx, alpha);
    (*this->hmmVec)[hmmIdx]->setSimFixedParams(currHMMParams);
  }
}

/*
 * this function copies all of the individual lib scaling factors into this's lib scaling factors
 */
void AllPairs3TrParam2DegPolyHMM::miscFunctions() {
  //std::cout << "In AllPairs3TrParam2DegPolyHMM::miscFunctions" << std::endl;
  int cell0Idx = 0;
  int cell1Idx = 0;
  for(unsigned int hmmIdx = 0; hmmIdx < this->hmmVec->size(); hmmIdx++) {
    // nullptr guard
    if((*this->hmmVec)[hmmIdx] == nullptr) {
      continue;
    }

    (*this->hmmVec)[hmmIdx]->miscFunctions();

    // save reestimated lib sizes into this object
    cell0Idx = getCell0IdxFromHMMIdx(hmmIdx);
    cell1Idx = getCell1IdxFromHMMIdx(hmmIdx);
    //this->setLibScalingFactor(cell0Idx, (*this->hmmVec)[hmmIdx]->getLibScalingFactor(0));
    //this->setLibScalingFactor(cell1Idx, (*this->hmmVec)[hmmIdx]->getLibScalingFactor(1));
    this->setLibScalingFactors(cell0Idx, cell1Idx, (*this->hmmVec)[hmmIdx]->getLibScalingFactor(0), (*this->hmmVec)[hmmIdx]->getLibScalingFactor(1));
  }

  // added for more granular printGradientPerIter, but this only prints lib scaling factors; tags these lines with an @ symbol
  fprintf(stderr, "@");
  for(std::vector<HMM*>::iterator itr = this->hmmVec->begin(); itr != this->hmmVec->end(); ++itr) {
    // nullptr guard
    if((*itr) == nullptr) {
      continue;
    }
    for(int cellNum = 0; cellNum < (*itr)->NUM_CELLS; cellNum++) {
      fprintf(stderr, "%.8f\t", (*itr)->getLibScalingFactor(cellNum));
    }
  }
  fprintf(stderr, "\n");
}

AllPairs3TrParam2DegPolyHMM* AllPairs3TrParam2DegPolyHMM::bfgs(gsl_vector* initGuess, int maxIters, bool verbose, bool debug) {
  // create new HMM with the best guess parameters and return it
  //std::cout << "subclass::bfgs this paramsToEst: ";
  //printColVector(this->getParamsToEst());
  //AllPairs3TrParam2DegPolyHMM* bestGuessOptim = new AllPairs3TrParam2DegPolyHMM(*this);
  //AllPairs3TrParam2DegPolyHMM* bestGuessOptim = AllPairs3TrParam2DegPolyHMM::create(*this);
  AllPairs3TrParam2DegPolyHMM* bestGuessOptim = this;//AllPairs3TrParam2DegPolyHMM::create(*this);
  //bestGuessOptim->print(stdout);
  //std::cout << "subclass::bfgs bestGuessOptim paramsToEst: ";
  //printColVector(bestGuessOptim->getParamsToEst());
  gsl_vector* initGuessAsParams = gsl_vector_alloc(initGuess->size);
  //std::cout << "initGuess: " << std::endl;;
  //printColVector(initGuess);
  this->convertProbToParam(initGuessAsParams, initGuess);
  //std::cout << "initGuessAsParams: " << std::endl;
  //printColVector(initGuessAsParams);
  //exit(0);
  Optimizable::bfgs(initGuessAsParams, bestGuessOptim, maxIters, verbose, debug);
  gsl_vector_free(initGuessAsParams);
  return bestGuessOptim;
}

void AllPairs3TrParam2DegPolyHMM::setInitGuessNthTime(gsl_vector* initGuess, int iter, int numTotalRuns) const {
  gsl_vector_set_zero(initGuess);
  if(iter == 0) {
    // lib scaling factors
    for(int cellIdx = 0; cellIdx < this->NUM_LIBS_TO_EST; cellIdx++) {
      gsl_vector_set(initGuess, this->LIB_SIZE_SCALING_FACTOR_START_IDX + cellIdx, 1);
    }

    // shared transition parameters
    //gsl_vector_set(initGuess, this->SHARED_TRANSITION_PROB_START_IDX + 0, 0.01); // alpha
    gsl_vector_set(initGuess, this->SHARED_TRANSITION_PROB_START_IDX + 0, 0.075); // beta
    gsl_vector_set(initGuess, this->SHARED_TRANSITION_PROB_START_IDX + 1, 0.075); // lambda

    // pairwise branch lengths
    for(int pairIdx = 0; pairIdx < this->NUM_PAIRS; pairIdx++) {
      gsl_vector_set(initGuess, this->BRANCH_LENGTH_START_IDX + 3 * pairIdx + 0, 0.225); // set t1
      gsl_vector_set(initGuess, this->BRANCH_LENGTH_START_IDX + 3 * pairIdx + 1, 0.225); // set t2
      gsl_vector_set(initGuess, this->BRANCH_LENGTH_START_IDX + 3 * pairIdx + 2, 0.225); // set t3
    }
  }
  else if(iter == 1) {
    // lib scaling factors
    for(int cellIdx = 0; cellIdx < this->NUM_LIBS_TO_EST; cellIdx++) {
      gsl_vector_set(initGuess, this->LIB_SIZE_SCALING_FACTOR_START_IDX + cellIdx, .75);
    }

    // shared transition parameters
    //gsl_vector_set(initGuess, this->SHARED_TRANSITION_PROB_START_IDX + 0, 0.05); // alpha
    gsl_vector_set(initGuess, this->SHARED_TRANSITION_PROB_START_IDX + 0, 0.025); // beta
    gsl_vector_set(initGuess, this->SHARED_TRANSITION_PROB_START_IDX + 1, 0.025); // lambda

    // pairwise branch lengths
    for(int pairIdx = 0; pairIdx < this->NUM_PAIRS; pairIdx++) {
      gsl_vector_set(initGuess, this->BRANCH_LENGTH_START_IDX + 3 * pairIdx + 0, 0.13); // set t1
      gsl_vector_set(initGuess, this->BRANCH_LENGTH_START_IDX + 3 * pairIdx + 1, 0.17); // set t2
      gsl_vector_set(initGuess, this->BRANCH_LENGTH_START_IDX + 3 * pairIdx + 2, 0.19); // set t3
    }
  }
  else if(iter == 2) {
    // lib scaling factors
    for(int cellIdx = 0; cellIdx < this->NUM_LIBS_TO_EST; cellIdx++) {
      gsl_vector_set(initGuess, this->LIB_SIZE_SCALING_FACTOR_START_IDX + cellIdx, 1.25);
    }

    // shared transition parameters
    //gsl_vector_set(initGuess, this->SHARED_TRANSITION_PROB_START_IDX + 0, 0.01); // alpha
    gsl_vector_set(initGuess, this->SHARED_TRANSITION_PROB_START_IDX + 0, 0.12); // beta
    gsl_vector_set(initGuess, this->SHARED_TRANSITION_PROB_START_IDX + 1, 0.09); // lambda

    // pairwise branch lengths
    for(int pairIdx = 0; pairIdx < this->NUM_PAIRS; pairIdx++) {
      gsl_vector_set(initGuess, this->BRANCH_LENGTH_START_IDX + 3 * pairIdx + 0, 0.02); // set t1
      gsl_vector_set(initGuess, this->BRANCH_LENGTH_START_IDX + 3 * pairIdx + 1, 0.02); // set t2
      gsl_vector_set(initGuess, this->BRANCH_LENGTH_START_IDX + 3 * pairIdx + 2, 0.02); // set t3
    }
  }
}

/*
 * internal method to actually calculate the sum of squared residuals for least squares. taken out of setBaumWelchInitGuess.
 * needs to be overridden by subclasses if there are different num and ordering of params
 */
double AllPairs3TrParam2DegPolyHMM::baumWelchLeastSquares_calcSumSqResid(const gsl_vector* v) {
  // convert v from BFGS space to prob space (like regular convertParamToProb)
  gsl_vector* probs = gsl_vector_alloc(v->size);
  this->baumWelchLeastSquares_convertParamToProb(probs, v);

  // nullptr guard
  HMM* hmm = this->getFirstNonNullHMM();

  // pull out relevant params (like setParamsToEst)
  int hmmTrIdx = hmm->TRANSITION_PROB_START_IDX - hmm->NUM_LIBS_TO_EST;
  int hmmBranchIdx = hmm->BRANCH_LENGTH_START_IDX - hmm->NUM_LIBS_TO_EST;
  int numHMMParams = hmm->getNumParamsToEst() - hmm->NUM_LIBS_TO_EST;
  gsl_vector* currHMMProbs = gsl_vector_alloc(numHMMParams);

  double beta  = gsl_vector_get(probs, this->SHARED_TRANSITION_PROB_START_IDX - this->NUM_LIBS_TO_EST + 0);
  double lambda = gsl_vector_get(probs, this->SHARED_TRANSITION_PROB_START_IDX - this->NUM_LIBS_TO_EST + 1);
  double currT1 = 0;
  double currT2 = 0;
  double currT3 = 0;

  double sumSqResid = 0;
  // for each pair HMM
  for(unsigned int hmmIdx = 0; hmmIdx < this->hmmVec->size(); hmmIdx++) {
    // nullptr guard
    if((*this->hmmVec)[hmmIdx] == nullptr) {
      continue;
    }

    // get apropriate branch lengths (there are 3 branch lengths stored per HMM
    currT1 = gsl_vector_get(probs, this->BRANCH_LENGTH_START_IDX - this->NUM_LIBS_TO_EST + 3 * hmmIdx + 0);
    currT2 = gsl_vector_get(probs, this->BRANCH_LENGTH_START_IDX - this->NUM_LIBS_TO_EST + 3 * hmmIdx + 1);
    currT3 = gsl_vector_get(probs, this->BRANCH_LENGTH_START_IDX - this->NUM_LIBS_TO_EST + 3 * hmmIdx + 2);

    //std::cout << currLib0BFGS << ", " << currLib1BFGS << ", " << alphaBFGS << ", " << betaBFGS << ", " << gammaBFGS << ", " << currT2BFGS << ", " << currT3BFGS << std::endl;
    // set everything into currHMMProbs
    gsl_vector_set(currHMMProbs, hmmTrIdx + 0, beta);
    gsl_vector_set(currHMMProbs, hmmTrIdx + 1, lambda);
    gsl_vector_set(currHMMProbs, hmmBranchIdx + 0, currT1);
    gsl_vector_set(currHMMProbs, hmmBranchIdx + 1, currT2);
    gsl_vector_set(currHMMProbs, hmmBranchIdx + 2, currT3);

    // call TwoCell3TrParam2DegPolyHMM::baumWelchLeastSquares_f(gsl_vector* probs)
    sumSqResid += (*this->hmmVec)[hmmIdx]->baumWelchLeastSquares_f(currHMMProbs);
  }
  gsl_vector_free(probs);
  gsl_vector_free(currHMMProbs);
  return sumSqResid;
}

/*
 * this is to be called during setBaumWelchInitGuess, in order to save params
 * estimated after least squares back into this->paramsToEst and initGuess.
 *
 * Each pairwise HMM already has params set correctly, just need to save them in this object. varsToEst_probSpace doesn't contain libs (varsToESt_probSpace is from least squares), so those need to be saved first (if any)
 *
 * //varsToEst_probSpace = [beta, gamma, t1_cell0_1, t2_cell0_1, t3_cell0_1, t1_cell0_2, t2_cell0_2, t3_cell0_2, ..., t3_cell(N-1)_N]
 * //initGuess = paramsToEst = [lib0, lib1, ..., libN, beta, gamma, t1_cell0_1, t2_cell0_1, t3_cell0_1, t1_cell0_2, t2_cell0_2, t3_cell0_2, ..., t3_cell(N-1)_N]
 * varsToEst_probSpace = [beta, lambda, t1_cell0_1, t2_cell0_1, t3_cell0_1, t1_cell0_2, t2_cell0_2, t3_cell0_2, ..., t3_cell(N-1)_N]
 * initGuess = paramsToEst = [lib0, lib1, ..., libN, beta, lambda, t1_cell0_1, t2_cell0_1, t3_cell0_1, t1_cell0_2, t2_cell0_2, t3_cell0_2, ..., t3_cell(N-1)_N]
 */
void AllPairs3TrParam2DegPolyHMM::saveBaumWelchEstsIntoParamsToEst(gsl_vector* varsToEst_probSpace, gsl_vector* initGuess) {
  gsl_vector* estsWithLibs = gsl_vector_alloc(this->getNumParamsToEst());
  gsl_vector_set_zero(estsWithLibs);
  int cell0Idx = 0;
  int cell1Idx = 0;
  double cell0Lib = 0;
  double cell1Lib = 0;
  int estsWithLibsIdx = this->LIB_SIZE_SCALING_FACTOR_START_IDX;
  for(unsigned int hmmIdx = 0; hmmIdx < this->hmmVec->size(); hmmIdx++) {
    // nullptr guard
    if((*this->hmmVec)[hmmIdx] == nullptr) {
      continue;
    }

    cell0Idx = getCell0IdxFromHMMIdx(hmmIdx);
    cell1Idx = getCell1IdxFromHMMIdx(hmmIdx);
    cell0Lib = (*this->hmmVec)[hmmIdx]->getLibScalingFactor(0);
    cell1Lib = (*this->hmmVec)[hmmIdx]->getLibScalingFactor(1);
    this->setLibScalingFactors(cell0Idx, cell1Idx, cell0Lib, cell1Lib);
    if(this->NUM_LIBS_TO_EST > 0) {
      gsl_vector_set(estsWithLibs, estsWithLibsIdx + 0, cell0Lib);
      gsl_vector_set(estsWithLibs, estsWithLibsIdx + 1, cell1Lib);
      estsWithLibsIdx += 2;
    }
  }
  // copy the rest of varsToEst_probSpace
  for(unsigned int varsToEstIdx = 0; varsToEstIdx < varsToEst_probSpace->size; varsToEstIdx++, estsWithLibsIdx++) {
    gsl_vector_set(estsWithLibs, estsWithLibsIdx, gsl_vector_get(varsToEst_probSpace, varsToEstIdx));
  }
  this->setParamsToEst(estsWithLibs);
  gsl_vector_memcpy(initGuess, estsWithLibs);
  //gsl_vector_free(varsToEst_probSpace);
  gsl_vector_free(estsWithLibs);
}

// same as convertProbToParam, just without library sizes. Fri 26 Jun 2020 05:36:29 PM PDT TODO depends on ordering of paramsToEst; there's probably a more elegant way of doing this
void AllPairs3TrParam2DegPolyHMM::baumWelchLeastSquares_convertProbToParam(gsl_vector* dest, const gsl_vector* src) const {
  ////std::cout << "AllPairs3TrParam2DegPolyHMM::baumWelchLeastSquares_convertProbToParam" << std::endl;
  ////std::cout << this->BRANCH_LENGTH_START_IDX - this->NUM_LIBS_TO_EST << std::endl;
  //// shared transition parameters
  //double a = this->getAlpha();
  //double b = gsl_vector_get(src, this->SHARED_TRANSITION_PROB_START_IDX - this->NUM_LIBS_TO_EST + 0);
  //double g = gsl_vector_get(src, this->SHARED_TRANSITION_PROB_START_IDX - this->NUM_LIBS_TO_EST + 1);
  //double c = (b * this->getKploidy() + g - 1);
  //gsl_vector_set(dest, this->SHARED_TRANSITION_PROB_START_IDX - this->NUM_LIBS_TO_EST + 0, log(-b * (1-2*a) * (double) this->getKploidy() / c)); // set y
  //gsl_vector_set(dest, this->SHARED_TRANSITION_PROB_START_IDX - this->NUM_LIBS_TO_EST + 1, log(-g * (1-2*a) / c)); // set z

  //// pairwise branch lengths
  //double d = (double) (*this->depthsVec)[0]->maxWindowSize;
  //double t1 = 0;
  //double t2 = 0;
  //double t3 = 0;
  //for(int pairIdx = 0; pairIdx < this->NUM_PAIRS; pairIdx++) {
  //  t1 = gsl_vector_get(src, this->BRANCH_LENGTH_START_IDX - this->NUM_LIBS_TO_EST + 3 * pairIdx + 0) / d;
  //  t2 = gsl_vector_get(src, this->BRANCH_LENGTH_START_IDX - this->NUM_LIBS_TO_EST + 3 * pairIdx + 1) / d;
  //  t3 = gsl_vector_get(src, this->BRANCH_LENGTH_START_IDX - this->NUM_LIBS_TO_EST + 3 * pairIdx + 2) / d;
  //  c = (d * (t1 + t2 + t3) - 1);
  //  gsl_vector_set(dest, this->BRANCH_LENGTH_START_IDX - this->NUM_LIBS_TO_EST + 3 * pairIdx + 0, log(-(d * t1) / c)); // set T1
  //  gsl_vector_set(dest, this->BRANCH_LENGTH_START_IDX - this->NUM_LIBS_TO_EST + 3 * pairIdx + 1, log(-(d * t2) / c)); // set T2
  //  gsl_vector_set(dest, this->BRANCH_LENGTH_START_IDX - this->NUM_LIBS_TO_EST + 3 * pairIdx + 2, log(-(d * t3) / c)); // set T3
  //}

  //std::cout << "AllPairs3TrParam2DegPolyHMM::baumWelchLeastSquares_convertProbToParam" << std::endl;
  // shared transition parameters
  double beta = gsl_vector_get(src, this->SHARED_TRANSITION_PROB_START_IDX - this->NUM_LIBS_TO_EST + 0);
  double lambda = gsl_vector_get(src, this->SHARED_TRANSITION_PROB_START_IDX - this->NUM_LIBS_TO_EST + 1);
  gsl_vector_set(dest, this->SHARED_TRANSITION_PROB_START_IDX - this->NUM_LIBS_TO_EST + 0, log(beta)); // set y
  gsl_vector_set(dest, this->SHARED_TRANSITION_PROB_START_IDX - this->NUM_LIBS_TO_EST + 1, log(lambda)); // set z

  // pairwise branch lengths
  double t1 = 0;
  double t2 = 0;
  double t3 = 0;
  for(int pairIdx = 0; pairIdx < this->NUM_PAIRS; pairIdx++) {
    t1 = gsl_vector_get(src, this->BRANCH_LENGTH_START_IDX - this->NUM_LIBS_TO_EST + 3 * pairIdx + 0);
    t2 = gsl_vector_get(src, this->BRANCH_LENGTH_START_IDX - this->NUM_LIBS_TO_EST + 3 * pairIdx + 1);
    t3 = gsl_vector_get(src, this->BRANCH_LENGTH_START_IDX - this->NUM_LIBS_TO_EST + 3 * pairIdx + 2);
    gsl_vector_set(dest, this->BRANCH_LENGTH_START_IDX - this->NUM_LIBS_TO_EST + 3 * pairIdx + 0, log(t1)); // set T1
    gsl_vector_set(dest, this->BRANCH_LENGTH_START_IDX - this->NUM_LIBS_TO_EST + 3 * pairIdx + 1, log(t2)); // set T2
    gsl_vector_set(dest, this->BRANCH_LENGTH_START_IDX - this->NUM_LIBS_TO_EST + 3 * pairIdx + 2, log(t3)); // set T3
  }
}

void AllPairs3TrParam2DegPolyHMM::baumWelchLeastSquares_convertParamToProb(gsl_vector* dest, const gsl_vector* src) const {
  ////std::cout << "AllPairs3TrParam2DegPolyHMM::baumWelchLeastSquares_convertParamToProb" << std::endl;
  //// shared transition parameters
  //double a = this->getAlpha(); //this->alpha;
  //double y = gsl_vector_get(src, this->SHARED_TRANSITION_PROB_START_IDX - this->NUM_LIBS_TO_EST + 0);
  //double z = gsl_vector_get(src, this->SHARED_TRANSITION_PROB_START_IDX - this->NUM_LIBS_TO_EST + 1);
  //double c = 1 - 2.0*a + exp(y) + exp(z);
  //gsl_vector_set(dest, this->SHARED_TRANSITION_PROB_START_IDX - this->NUM_LIBS_TO_EST + 0, exp(y) / ((double) this->getKploidy() * c)); // beta
  //gsl_vector_set(dest, this->SHARED_TRANSITION_PROB_START_IDX - this->NUM_LIBS_TO_EST + 1, exp(z) / c); // gamma

  //// pairwise branch lengths
  //double T1 = 0;
  //double T2 = 0;
  //double T3 = 0;
  //for(int pairIdx = 0; pairIdx < this->NUM_PAIRS; pairIdx++) {
  //  T1 = gsl_vector_get(src, this->BRANCH_LENGTH_START_IDX - this->NUM_LIBS_TO_EST + 3 * pairIdx + 0);
  //  T2 = gsl_vector_get(src, this->BRANCH_LENGTH_START_IDX - this->NUM_LIBS_TO_EST + 3 * pairIdx + 1);
  //  T3 = gsl_vector_get(src, this->BRANCH_LENGTH_START_IDX - this->NUM_LIBS_TO_EST + 3 * pairIdx + 2);
  //  c = 1.0 / (1 + exp(T1) + exp(T2) + exp(T3));
  //  gsl_vector_set(dest, this->BRANCH_LENGTH_START_IDX - this->NUM_LIBS_TO_EST + 3 * pairIdx + 0, exp(T1) * c); // set t1
  //  gsl_vector_set(dest, this->BRANCH_LENGTH_START_IDX - this->NUM_LIBS_TO_EST + 3 * pairIdx + 1, exp(T2) * c); // set t2
  //  gsl_vector_set(dest, this->BRANCH_LENGTH_START_IDX - this->NUM_LIBS_TO_EST + 3 * pairIdx + 2, exp(T3) * c); // set t3
  //}

  //std::cout << "AllPairs3TrParam2DegPolyHMM::baumWelchLeastSquares_convertParamToProb" << std::endl;
  //std::cout << "beta src idx: " << this->SHARED_TRANSITION_PROB_START_IDX - this->NUM_LIBS_TO_EST << ", branch idx start: " << this->BRANCH_LENGTH_START_IDX - this->NUM_LIBS_TO_EST << std::endl;
  // shared transition parameters
  double y = gsl_vector_get(src, this->SHARED_TRANSITION_PROB_START_IDX - this->NUM_LIBS_TO_EST + 0);
  double z = gsl_vector_get(src, this->SHARED_TRANSITION_PROB_START_IDX - this->NUM_LIBS_TO_EST + 1);
  gsl_vector_set(dest, this->SHARED_TRANSITION_PROB_START_IDX - this->NUM_LIBS_TO_EST + 0, exp(y)); // beta
  gsl_vector_set(dest, this->SHARED_TRANSITION_PROB_START_IDX - this->NUM_LIBS_TO_EST + 1, exp(z)); // lambda

  // pairwise branch lengths
  double T1 = 0;
  double T2 = 0;
  double T3 = 0;
  for(int pairIdx = 0; pairIdx < this->NUM_PAIRS; pairIdx++) {
    T1 = gsl_vector_get(src, this->BRANCH_LENGTH_START_IDX - this->NUM_LIBS_TO_EST + 3 * pairIdx + 0);
    T2 = gsl_vector_get(src, this->BRANCH_LENGTH_START_IDX - this->NUM_LIBS_TO_EST + 3 * pairIdx + 1);
    T3 = gsl_vector_get(src, this->BRANCH_LENGTH_START_IDX - this->NUM_LIBS_TO_EST + 3 * pairIdx + 2);
    gsl_vector_set(dest, this->BRANCH_LENGTH_START_IDX - this->NUM_LIBS_TO_EST + 3 * pairIdx + 0, exp(T1)); // set t1
    gsl_vector_set(dest, this->BRANCH_LENGTH_START_IDX - this->NUM_LIBS_TO_EST + 3 * pairIdx + 1, exp(T2)); // set t2
    gsl_vector_set(dest, this->BRANCH_LENGTH_START_IDX - this->NUM_LIBS_TO_EST + 3 * pairIdx + 2, exp(T3)); // set t3
  }
}

/*
 * function to save all CNAs to bed files
 */
void AllPairs3TrParam2DegPolyHMM::saveAllCNAToBed(std::string filename) {
  // for each HMM, save each cell to a bed file
  for(unsigned int hmmIdx = 0; hmmIdx < this->hmmVec->size(); hmmIdx++) {
    HMM* currHMM = (*this->hmmVec)[hmmIdx];
    // nullptr guard
    if(currHMM == nullptr) {
      continue;
    }
    std::string hmmName = (*this->hmmNames)[hmmIdx];
    std::vector<std::string> cellNames;
    boost::split(cellNames, hmmName, boost::is_any_of(","));
    boost::replace_all(hmmName, ",", "_");

    for(int cellIdx = 0; cellIdx < currHMM->NUM_CELLS && cellIdx < (int) cellNames.size(); cellIdx++) {
      currHMM->saveCNAToBed(filename + "__pair_" + hmmName + "__cell_" + cellNames[cellIdx] + "__k" + std::to_string(this->MAX_PLOIDY) + ".bed", cellIdx);
    }
  }
}

/*
 * helper function to return a string for the summaryMethod
 */
std::string AllPairs3TrParam2DegPolyHMM::getSummaryMethodName(int summaryMethod) {
  if(summaryMethod == AllPairs3TrParam2DegPolyHMM::SUMMARY_MEAN) {
    return "mean";
  } else if (summaryMethod == AllPairs3TrParam2DegPolyHMM::SUMMARY_MEDIAN) {
    return "median";
  } else if (summaryMethod == AllPairs3TrParam2DegPolyHMM::SUMMARY_MODE) {
    return "mode";
  } else if (summaryMethod == AllPairs3TrParam2DegPolyHMM::SUMMARY_MARGINAL) {
    return "marg";
  }
  return "unknown";
}

/*
 * wrapper function around summary decoding methods. summaryMethod should be one of the static consts in this class
 */
void AllPairs3TrParam2DegPolyHMM::summaryDecodeAllCellsAcrossPairs(int summaryMethod, bool runDecoding) {
  this->summaryMethodUsed = summaryMethod;

  // allocate intermediates
  delete this->summaryPathAcrossPairsVec;
  this->summaryPathAcrossPairsVec = new std::vector<std::unordered_map<std::string, std::vector<double>*>*>();
  for(int cellIdx = 0; cellIdx < this->NUM_CELLS; cellIdx++) {
    (this->summaryPathAcrossPairsVec)->push_back(new std::unordered_map<std::string, std::vector<double>*>());
  }

  if(summaryMethod == AllPairs3TrParam2DegPolyHMM::SUMMARY_MARGINAL) {
    this->decodeMargAllCellsAcrossPairs();
    return;
  }

  // if using mean, median, or mode, need the viterbi decoding first, unless have already run it previously (ie summarizing across multiple methods, no need to keep doing viterbi decoding since it won't change)
  if(runDecoding) {
    this->viterbiDecodeAll();
  }
  this->decodeStatSummaryAllCellsAcrossPairs(summaryMethod);
}

/*
 * function to decode marginally for one cell, across all pairs, by taking the maximum likelihood
 * from each pair's forBackMat
 * assigns margLikelihoodMatAcrossPairsVec and summaryPathAcrossPairsVec
 */
void AllPairs3TrParam2DegPolyHMM::decodeMargAllCellsAcrossPairs() {
  // call each HMM's calcMargLikelihoods()
  for(unsigned int hmmIdx = 0; hmmIdx < this->hmmVec->size(); hmmIdx++) {
    // nullptr guard
    if((*this->hmmVec)[hmmIdx] == nullptr) {
      continue;
    }
    (*this->hmmVec)[hmmIdx]->calcMargLikelihoods();
  }

  // allocate intermediates
  this->margLikelihoodMatAcrossPairsVec = new std::vector<std::unordered_map<std::string, gsl_matrix*>*>();

  // for each cell, call decodeMargOneCellAcrossPairs
  for(int cellIdx = 0; cellIdx < this->NUM_CELLS; cellIdx++) {
    // alloc intermediates
    (this->margLikelihoodMatAcrossPairsVec)->push_back(new std::unordered_map<std::string, gsl_matrix*>());
    this->decodeMargOneCellAcrossPairs(cellIdx);
  }
}

/*
 * function that will find the most likely copy number state of a particular cell, using the forBackMargMatMap field from each DepthPair.
 * sets the entries for cellNum in margLikelihoodMatAcrossPairsVec and summaryPathAcrossPairsVec
 */
void AllPairs3TrParam2DegPolyHMM::decodeMargOneCellAcrossPairs(int cellNum) {
  HMM* currHMM = this->getFirstNonNullHMM();
  std::vector<std::string>* chrVec = currHMM->getChrVec();

  gsl_matrix* summedLogLikMat = nullptr;
  gsl_matrix* currMatToAdd = nullptr;
  std::string currChr;
  int currNumWindows = -1;
  int currMaxIdx = -1;
  std::vector<double>* currPath = nullptr;
  std::unordered_map<std::string, gsl_matrix*>* currForBackMargMatMap = nullptr;

  // identify which hmms this cell is in
  std::vector<HMM*>* hmmsWithCellAsCell0 = this->getHMMsWithCellAsCell0(cellNum);
  std::vector<HMM*>* hmmsWithCellAsCell1 = this->getHMMsWithCellAsCell1(cellNum);

  // for each chr
  for(unsigned int chrIdx = 0; chrIdx < chrVec->size(); chrIdx++) {
    // allocate summedLogLikMat for this chr
    currChr = (*chrVec)[chrIdx];
    currNumWindows = (*(*this->depthsVec)[0]->regions)[currChr]->size();
    summedLogLikMat = gsl_matrix_alloc(this->MAX_PLOIDY+1, currNumWindows);
    gsl_matrix_set_zero(summedLogLikMat);

    // for each hmm this cell is in as the left in the pair (cell0)
    for(unsigned int hmmIdx = 0; hmmIdx < hmmsWithCellAsCell0->size(); hmmIdx++) {
      // go to hmmIdx, get cell0's forBackMargMatMap
      currHMM = (*hmmsWithCellAsCell0)[hmmIdx];
      //currMatToAdd = (*(*currHMM->getDepths())[0]->forBackMargMatMap)[currChr];
      currForBackMargMatMap = (*currHMM->getForBackMargMatMapVec())[0];
      currMatToAdd = (*currForBackMargMatMap)[currChr];

      // add to running total
      gsl_matrix_add(summedLogLikMat, currMatToAdd);
    }

    // for each hmm this cell is in as the right in the pair (cell1)
    for(unsigned int hmmIdx = 0; hmmIdx < hmmsWithCellAsCell1->size(); hmmIdx++) {
      currHMM = (*hmmsWithCellAsCell1)[hmmIdx];
      //currMatToAdd = (*(*currHMM->getDepths())[1]->forBackMargMatMap)[currChr];
      currForBackMargMatMap = (*currHMM->getForBackMargMatMapVec())[1];
      currMatToAdd = (*currForBackMargMatMap)[currChr];
      gsl_matrix_add(summedLogLikMat, currMatToAdd);
    }

    // store summedLogLikMat for this chr in this->margLikelihoodMatAcrossPairsVec[cellNum]
    (*(*this->margLikelihoodMatAcrossPairsVec)[cellNum])[currChr] = summedLogLikMat;

    // allocate a path vector
    currPath = new std::vector<double>();

    // for each position in chr
    for(unsigned int i = 0; i < summedLogLikMat->size2; i++) {
      // find argmax in that col
      gsl_vector_view currCol = gsl_matrix_column(summedLogLikMat, i);
      currMaxIdx = gsl_vector_max_index(&currCol.vector);

      // store in path vector
      currPath->push_back(currMaxIdx);
    }

    // store path vector in this->summaryPathAcrossPairsVec
    (*(*this->summaryPathAcrossPairsVec)[cellNum])[currChr] = currPath;
  }

  // clean up
  delete hmmsWithCellAsCell0;
  delete hmmsWithCellAsCell1;
}

void AllPairs3TrParam2DegPolyHMM::decodeStatSummaryAllCellsAcrossPairs(int summaryMethod) {
  for(int cellIdx = 0; cellIdx < this->NUM_CELLS; cellIdx++) {
    this->decodeStatSummaryOneCellAcrossPairs(summaryMethod, cellIdx);
  }
}

/*
 * for each cell, finds the hmm pairs it's part of and summarizies the viterbi decodings by the passed summaryMethod
 */
void AllPairs3TrParam2DegPolyHMM::decodeStatSummaryOneCellAcrossPairs(int summaryMethod, int cellNum) {
  // intermediates
  HMM* currHMM = this->getFirstNonNullHMM();
  std::vector<std::string>* chrVec = currHMM->getChrVec();
  DepthPair* firstDepthPair = (*this->depthsVec)[0];
  std::string currChr;
  std::unordered_map<std::string, std::vector<int>*>* currVitPathMap = nullptr;
  int ploidiesIdx = -1;
  int currPloidy = -1;
  int ploidyCount = 0;
  double currSummaryPloidy = -1;
  std::vector<double>* currPath = nullptr;

  // identify which hmms this cell is in
  std::vector<HMM*>* hmmsWithCellAsCell0 = this->getHMMsWithCellAsCell0(cellNum);
  std::vector<HMM*>* hmmsWithCellAsCell1 = this->getHMMsWithCellAsCell1(cellNum);
  //std::cout << "cell " << cellNum << "is 0 in: " << std::endl;
  //for(unsigned int i = 0; i < hmmsWithCellAsCell0->size(); i++) {
  //  std::cout << (*hmmsWithCellAsCell0)[i] << "\t";
  //}
  //std::cout << std::endl;

  //std::cout << "cell " << cellNum << "is 1 in: " << std::endl;
  //for(unsigned int i = 0; i < hmmsWithCellAsCell1->size(); i++) {
  //  std::cout << (*hmmsWithCellAsCell1)[i] << "\t";
  //}
  //std::cout << std::endl;

  // intermediate vector. stores all calc ploidies for a given genomic window
  gsl_vector* ploidies = nullptr;
  // if using mode, count occurences of each state so can calc max idx
  if(summaryMethod == AllPairs3TrParam2DegPolyHMM::SUMMARY_MODE) {
    ploidies = gsl_vector_alloc(this->MAX_PLOIDY + 1); // state:count
  }
  // otherwise, store a vec of every state so can calc mean/median
  else {
    //ploidies = gsl_vector_alloc(this->NUM_CELLS - 1); // each cell is in n-1 pairs. one entry for each pair
    ploidies = gsl_vector_alloc(hmmsWithCellAsCell0->size() + hmmsWithCellAsCell1->size()); // only alloc for the number of HMMs this cell appears in, else will introduce scaling errors
  }

  // for each chr
  for(unsigned int chrIdx = 0; chrIdx < chrVec->size(); chrIdx++) {
    currChr = (*chrVec)[chrIdx];

    // allocate a path vector
    currPath = new std::vector<double>();

    // for each window in currChr
    for(unsigned int regionIdx = 0; regionIdx < (*firstDepthPair->regions)[currChr]->size(); regionIdx++) {
      ploidiesIdx = 0;
      gsl_vector_set_zero(ploidies);

      // for each hmm this cell is in as the left in the pair (cell0)
      for(unsigned int hmmIdx = 0; hmmIdx < hmmsWithCellAsCell0->size(); hmmIdx++) {
        // store viterbi decoding for this idx
        currHMM = (*hmmsWithCellAsCell0)[hmmIdx];
        currVitPathMap = (*currHMM->getChrToViterbiPathMapVec())[0];
        currPloidy = (*(*currVitPathMap)[currChr])[regionIdx];
        if(summaryMethod == AllPairs3TrParam2DegPolyHMM::SUMMARY_MODE) {
          // store count + 1
          ploidyCount = gsl_vector_get(ploidies, currPloidy);
          gsl_vector_set(ploidies, currPloidy, ploidyCount + 1);
        }
        else {
          gsl_vector_set(ploidies, ploidiesIdx, currPloidy);
          ploidiesIdx++;
        }
      }
      // for each hmm this cell is in as the right in the pair (cell1)
      for(unsigned int hmmIdx = 0; hmmIdx < hmmsWithCellAsCell1->size(); hmmIdx++) {
        currHMM = (*hmmsWithCellAsCell1)[hmmIdx];
        currVitPathMap = (*currHMM->getChrToViterbiPathMapVec())[1];
        currPloidy = (*(*currVitPathMap)[currChr])[regionIdx];
        if(summaryMethod == AllPairs3TrParam2DegPolyHMM::SUMMARY_MODE) {
          ploidyCount = gsl_vector_get(ploidies, currPloidy);
          gsl_vector_set(ploidies, currPloidy, ploidyCount + 1);
        }
        else {
          gsl_vector_set(ploidies, ploidiesIdx, currPloidy);
          ploidiesIdx++;
        }
      }

      // calc summary ploidy for this window idx
      if(summaryMethod == AllPairs3TrParam2DegPolyHMM::SUMMARY_MEAN) {
        currSummaryPloidy = gsl_stats_mean(ploidies->data, 1, ploidies->size);
      }
      else if(summaryMethod == AllPairs3TrParam2DegPolyHMM::SUMMARY_MEDIAN) {
        gsl_sort(ploidies->data, 1, ploidies->size);
        currSummaryPloidy = gsl_stats_median_from_sorted_data(ploidies->data, 1, ploidies->size);
      }
      else if(summaryMethod == AllPairs3TrParam2DegPolyHMM::SUMMARY_MODE) {
        currSummaryPloidy = gsl_vector_max_index(ploidies);
      }

      // store summary ploidy in path
      currPath->push_back(currSummaryPloidy);
    }

    // store path vector in this->summaryPathAcrossPairsVec
    (*(*this->summaryPathAcrossPairsVec)[cellNum])[currChr] = currPath;
  }

  // clean up
  gsl_vector_free(ploidies);
  delete hmmsWithCellAsCell0;
  delete hmmsWithCellAsCell1;
}


/*
 * function to save each cell's path vector from the marginal across pairs decoding, similar to HMM::saveViterbiDecodedCNA
 * Assumes decodeMargAllCellsAcrossPairs was called first (to set summaryPathAcrossPairsVec), no safety checks are done
 * coord | diploid_mean | diploid_var | tumor0 | summaryDecoded0_0-MAX_PLOIDY
 */
void AllPairs3TrParam2DegPolyHMM::saveAllSummaryDecodedCNA(std::string filename) {
  HMM* currHMM = this->getFirstNonNullHMM();
  std::vector<std::string>* chrVec = currHMM->getChrVec();
  DepthPair* firstDepthPair = (*this->depthsVec)[0];

  std::unordered_map<std::string, std::vector<double>*>* currSummaryPathMap = nullptr;
  DepthPair* currDepths = nullptr;
  std::string sep = "\t";
  std::string currChr;
  std::string currMethod = this->getSummaryMethodName(this->summaryMethodUsed);

  for(int cellIdx = 0; cellIdx < this->NUM_CELLS; cellIdx++) {
    std::string currCellName = (*this->sampleList)[cellIdx];
    std::string currFilename = filename + "__" + currCellName + "__k" + std::to_string(this->MAX_PLOIDY) + "." + currMethod + "Decoded";
    std::ofstream outFile(currFilename);

    currDepths = (*this->depthsVec)[cellIdx];
    currSummaryPathMap = (*summaryPathAcrossPairsVec)[cellIdx];

    // write header
    outFile << "coord\tdiploid_mean\tdiploid_var";
    if(firstDepthPair->chrToDiploidSimStateMap->size() > 0) {
      outFile << "\tdiploid_simState";
    }
    outFile << "\t" << currCellName << "\t" << currMethod << "Decoded" << cellIdx << "_0-" << this->MAX_PLOIDY;
    if(currDepths->chrToTumorSimStateMap->size() > 0) {
      outFile << "\tsimulated" << cellIdx << "_0-" << this->MAX_PLOIDY;
    }
    outFile << std::endl;

    // for each chr
    for(unsigned int i = 0; i < chrVec->size(); i++) {
      currChr = (*chrVec)[i];

      // for each window in currChr
      for(unsigned int regionIdx = 0; regionIdx < (*firstDepthPair->regions)[currChr]->size(); regionIdx++) {
        // coord
        outFile << (*(*firstDepthPair->regions)[currChr])[regionIdx] << sep;

        // diploid mean
        outFile << (*(*firstDepthPair->chrToDiploidDepthMap)[currChr])[regionIdx] << sep;

        // diploid var
        outFile << (*(*firstDepthPair->chrToDiploidVarMap)[currChr])[regionIdx];

        // diploid simulated state (if set)
        if(firstDepthPair->chrToDiploidSimStateMap->size() > 0) {
          outFile << sep << (*(*firstDepthPair->chrToDiploidSimStateMap)[currChr])[regionIdx];
        }

        // tumor depth
        outFile << sep << (*(*currDepths->chrToTumorDepthMap)[currChr])[regionIdx];

        // decoded CNA
        outFile << sep << (*(*currSummaryPathMap)[currChr])[regionIdx];

        // tumor simulated state (if set)
        if(currDepths->chrToTumorSimStateMap->size() > 0) {
          outFile << sep << (*(*currDepths->chrToTumorSimStateMap)[currChr])[regionIdx];
        }

        // new line
        outFile << std::endl;
      }
    }
    outFile.close();
  }
}

/*
 * function to save each cell's path vector from the marginal across pairs decoding, similar to above but in bed format
 * Assumes summaryDecodeAllCellsAcrossPairs was called first (to set summaryPathAcrossPairsVec), no safety checks are done
 * chr | start | end | margDecoded
 */
void AllPairs3TrParam2DegPolyHMM::saveAllSummaryDecodedCNAToBed(std::string filename) {
  HMM* currHMM = this->getFirstNonNullHMM();
  std::vector<std::string>* chrVec = currHMM->getChrVec();
  DepthPair* firstDepthPair = (*this->depthsVec)[0];

  std::unordered_map<std::string, std::vector<double>*>* currSummaryPathMap = nullptr;
  DepthPair* currDepths = nullptr;
  std::string sep = "\t";
  std::string currChr;
  std::string currMethod = this->getSummaryMethodName(this->summaryMethodUsed);

  for(int cellIdx = 0; cellIdx < this->NUM_CELLS; cellIdx++) {
    std::string currCellName = (*this->sampleList)[cellIdx];
    std::string currFilename = filename + "__" + currCellName + "__k" + std::to_string(this->MAX_PLOIDY) + "__" + currMethod + ".bed";
    std::ofstream outFile(currFilename);

    currDepths = (*this->depthsVec)[cellIdx];
    currSummaryPathMap = (*summaryPathAcrossPairsVec)[cellIdx];

    // for each chr
    for(unsigned int i = 0; i < chrVec->size(); i++) {
      currChr = (*chrVec)[i];

      // for each window in currChr
      for(unsigned int regionIdx = 0; regionIdx < (*firstDepthPair->regions)[currChr]->size(); regionIdx++) {
        // transform stored coord "chr:start-end" ==> "chr\tstart\tend"
        std::string coord = (*(*currDepths->regions)[currChr])[regionIdx];
        boost::replace_all(coord, ":", sep);
        boost::replace_all(coord, "-", sep);
        outFile << coord << sep;

        // decoded CNA
        //outFile << (*(*currSummaryPathMap)[currChr])[regionIdx];
        outFile << std::setprecision(5) << (*(*currSummaryPathMap)[currChr])[regionIdx];

        // new line
        outFile << std::endl;
      }
    }
    outFile.close();
  }
}

/*
 * sets paramsToEst for each HMM file it reads, and saves into this->paramsToEst
 * assumes all HMM's have already been set up, now we're just overwriting preexisting this->paramsToEst values
 */
void AllPairs3TrParam2DegPolyHMM::getPairedOptimParamsToEstFromFiles(std::string pairedEstimatesPath, int numExpectedLinesPerFile) {
  this->pairedEstimatesPath = pairedEstimatesPath;
  HMM* currHMM = nullptr;
  for(unsigned int hmmIdx = 0; hmmIdx < this->hmmVec->size(); hmmIdx++) {
    // nullptr guard, if we don't care about this HMM
    currHMM = (*this->hmmVec)[hmmIdx];
    if(currHMM == nullptr) {
      continue;
    }
    this->getOnePairedOptimParamsToEstFromFile(hmmIdx, pairedEstimatesPath, numExpectedLinesPerFile);
  }
}

void AllPairs3TrParam2DegPolyHMM::getOnePairedOptimParamsToEstFromFile(int hmmIdx, std::string pairedEstimatesPath, int numExpectedLinesPerFile) {
  HMM* currHMM = (*this->hmmVec)[hmmIdx];

  std::string hmmName = (*this->hmmNames)[hmmIdx];
  boost::replace_all(hmmName, ",", "__");
  std::string currFile = pairedEstimatesPath + "__" + hmmName + "__k" + std::to_string(this->MAX_PLOIDY) + ".pairedParams";
  bool currFileExists = false;
  // check both orderings of cells in hmmName
  if(boost::filesystem::exists(currFile)) {
    currFileExists = true;
  }
  else {
    std::vector<std::string> cellNames;
    hmmName = (*this->hmmNames)[hmmIdx];
    boost::split(cellNames, hmmName, boost::is_any_of(","));
    hmmName = cellNames[1] + "__" + cellNames[0]; // swap ordering
    currFile = pairedEstimatesPath + "__" + hmmName + "__k" + std::to_string(this->MAX_PLOIDY) + ".pairedParams";
    if(boost::filesystem::exists(currFile)) {
      currFileExists = true;
    }
    else {
      currFileExists = false;
    }
  }
  // file exists guard, if we didn't get to this HMM in previous runs
  if(currFileExists) {
    double status = currHMM->getParamsToEstFromFile(currFile, numExpectedLinesPerFile);
    // if error occurred, then don't use the estimates read from file
    if(gsl_isnan(status)) {
      currHMM->hasBeenReadFromFile = false;
    }
    else {
      // if successfully read from file, then this HMM is finalized
      currHMM->hasBeenReadFromFile = true;
      currHMM->finalLl = currHMM->getLogLikelihood();
    }
  }
  else {
    currHMM->hasBeenReadFromFile = false;
  }

  // save into this->paramsToEst
  gsl_vector* currParamsToEst = currHMM->getParamsToEst();
  for(unsigned int i = 0; i < currParamsToEst->size; i++) {
    gsl_vector_set(this->paramsToEst, 3 * hmmIdx + i, gsl_vector_get(currParamsToEst, i));
  }
}

