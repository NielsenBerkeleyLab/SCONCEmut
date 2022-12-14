#include "AllCells3TrParam2DegPolyHMM.hpp"

// ctors and destructor
AllCells3TrParam2DegPolyHMM::AllCells3TrParam2DegPolyHMM(std::vector<DepthPair*>* depths, std::vector<std::string>* sampleList, gsl_vector* fixedParams,  int maxPloidy, int numSharedTrParamsToEst, int numFixedTrParams, int numFixedLibs, int numHMMs, int numBranchesToEst) :
                                        MAX_PLOIDY(maxPloidy),
                                        NUM_CELLS(depths->size()),
                                        NUM_HMMS(numHMMs),
                                        NUM_LIBS_TO_EST(this->NUM_CELLS - numFixedLibs),
                                        NUM_SHARED_TRANSITION_PARAMS_TO_EST(numSharedTrParamsToEst),
                                        NUM_BRANCH_LENGTHS_TO_EST(numBranchesToEst),
                                        NUM_FIXED_LIBS(numFixedLibs),
                                        NUM_FIXED_SHARED_TRANSITION_PARAMS(numFixedTrParams),
                                        LIB_SIZE_SCALING_FACTOR_START_IDX(0),
                                        SHARED_TRANSITION_PROB_START_IDX(this->NUM_CELLS - numFixedLibs), // libs go first, estimating one for each "unfixed lib" cell
                                        BRANCH_LENGTH_START_IDX(this->SHARED_TRANSITION_PROB_START_IDX + numSharedTrParamsToEst),
                                        FIXED_TRANSITION_PROB_START_IDX(numFixedLibs), // if no libs are fixed, then this starts at 0
                                        FIXED_MUTATION_PARAM_START_IDX(this->NUM_FIXED_LIBS + this->NUM_FIXED_SHARED_TRANSITION_PARAMS) // only used for classes *WithMuts
                                        {
  this->depthsVec = depths;
  this->sampleList = sampleList;
  this->hmmVec = new std::vector<HMM*>(this->NUM_HMMS);
  this->hmmNames = nullptr;
  this->baumWelchParamResults = nullptr;
  this->meanVarianceCoefVec = nullptr;
  this->paramsToEst = gsl_vector_alloc(this->NUM_LIBS_TO_EST + this->NUM_SHARED_TRANSITION_PARAMS_TO_EST + this->NUM_BRANCH_LENGTHS_TO_EST); // lib scaling for each cell + b/L + t1/t2/t3 for each pair
  gsl_vector_set_all(this->paramsToEst, 1.0);
  if(fixedParams == nullptr) {
    this->fixedParams = gsl_vector_alloc(numFixedLibs + numFixedTrParams);
    gsl_vector_set_zero(this->fixedParams);
  }
  else {
    this->fixedParams = fixedParams;
  }
  this->probParamConversionVec = gsl_vector_alloc(this->getNumParamsToEst());
  this->bwGradientIntermediateVec = gsl_vector_alloc(this->NUM_SHARED_TRANSITION_PARAMS_TO_EST + this->NUM_BRANCH_LENGTHS_TO_EST);
  this->maxNumBFGSStarts = 0;

  this->setAlpha(0.1); // arbitrary starting value
}

AllCells3TrParam2DegPolyHMM::~AllCells3TrParam2DegPolyHMM() {
  for(unsigned int hmmIdx = 0; hmmIdx < this->hmmVec->size(); hmmIdx++) {
    delete (*this->hmmVec)[hmmIdx];
  }
  delete this->hmmVec;
  if(this->baumWelchParamResults != nullptr) {
    for(std::vector<gsl_vector*>::iterator itr = this->baumWelchParamResults->begin(); itr != this->baumWelchParamResults->end(); ++itr) {
      gsl_vector_free(*itr);
    }
  }
  delete this->baumWelchParamResults;
  gsl_vector_free(this->meanVarianceCoefVec);
  delete this->bwGradientIntermediateVec;
}

int AllCells3TrParam2DegPolyHMM::getKploidy() const {
  return this->MAX_PLOIDY;
}

/*
 * function to return a chromosome vector from the first non null HMM.
 * the SelectPairs* classes have some nullptr entries in their HMM vectors so need to search for the first non null one (ie can't assume the 0th entry is non null)
 */
HMM* AllCells3TrParam2DegPolyHMM::getFirstNonNullHMM() const {
  HMM* currHMM = nullptr;
  for(unsigned int hmmIdx = 0; hmmIdx < this->hmmVec->size(); hmmIdx++) {
    currHMM = (*this->hmmVec)[hmmIdx];
    if(currHMM != nullptr) {
      return currHMM;
    }
  }
  return nullptr;
}

std::vector<HMM*>* AllCells3TrParam2DegPolyHMM::getHMMs() {
  return this->hmmVec;
}

void AllCells3TrParam2DegPolyHMM::print(FILE* stream) {
  // print each HMM
  for(unsigned int hmmIdx = 0; hmmIdx < this->hmmVec->size(); hmmIdx++) {
    // nullptr guard
    if((*this->hmmVec)[hmmIdx] == nullptr) {
      continue;
    }

    fprintf(stream, "HMM %i (%s):\n", hmmIdx, (*this->hmmNames)[hmmIdx].c_str());
    (*this->hmmVec)[hmmIdx]->print(stream);
    fprintf(stream, "\n");
  }

  fprintf(stream, "\n\n");

  // print paramsToEst
  fprintf(stream, "AllCells_TrParam2DegPolyHMM paramsToEst:\n");
  printColVector(stream, this->paramsToEst);

  // print fixedParams
  fprintf(stream, "AllCells_TrParam2DegPolyHMM fixedParams:\n");
  printColVector(stream, this->fixedParams);

  // print total loglikelihood
  // finalLl should only be set right before freeToMinimum (ie we clear everything else out)
  double ll = 0;
  if(this->finalLl != -std::numeric_limits<double>::infinity()) {
    ll = this->finalLl;
  }
  else {
    ll = this->getLogLikelihood();
  }
  fprintf(stream, "AllCells3TrParam2DegPolyHMM total loglikelihood: %.10f\n\n", ll);
}

std::vector<std::string>* AllCells3TrParam2DegPolyHMM::getSampleList() {
  return this->sampleList;
}

void AllCells3TrParam2DegPolyHMM::setSampleList(std::vector<std::string>* sampleList) {
  if(this->sampleList != nullptr && this->sampleList != sampleList) {
    delete this->sampleList;
  }
  this->sampleList = sampleList;
}

std::vector<gsl_vector*>* AllCells3TrParam2DegPolyHMM::getBaumWelchParamResults() const {
  return this->baumWelchParamResults;
}
void AllCells3TrParam2DegPolyHMM::setBaumWelchParamResults(std::vector<gsl_vector*>* paramResults) {
  this->baumWelchParamResults = paramResults;
}

double AllCells3TrParam2DegPolyHMM::getAlpha() const {
  return gsl_vector_get(this->fixedParams, this->FIXED_TRANSITION_PROB_START_IDX);
}

void AllCells3TrParam2DegPolyHMM::setAlpha(double alpha) {
  gsl_vector_set(this->fixedParams, this->FIXED_TRANSITION_PROB_START_IDX, alpha);
}

void AllCells3TrParam2DegPolyHMM::setAllMeanVarianceFn(gsl_vector* meanVarianceCoefVec) {
  // save in this object
  if(this->meanVarianceCoefVec == nullptr) {
    this->meanVarianceCoefVec = gsl_vector_alloc(meanVarianceCoefVec->size);
  }
  gsl_vector_memcpy(this->meanVarianceCoefVec, meanVarianceCoefVec);

  // save for each HMM
  gsl_vector* currMeanVarCoefVec = nullptr;
  HMM* hmm = nullptr;
  for(unsigned int hmmIdx = 0; hmmIdx < this->hmmVec->size(); hmmIdx++) {
    // nullptr guard
    if((*this->hmmVec)[hmmIdx] == nullptr) {
      continue;
    }

    hmm = (*this->hmmVec)[hmmIdx];
    currMeanVarCoefVec = gsl_vector_alloc(meanVarianceCoefVec->size);
    gsl_vector_memcpy(currMeanVarCoefVec, meanVarianceCoefVec);
    hmm->setMeanVarianceFn(currMeanVarCoefVec);
  }
}

// Optimizable methods
/*
 * same as HMM::checkOptimProbValidity
 */
double AllCells3TrParam2DegPolyHMM::checkOptimProbValidity(gsl_vector* probs) const {
  // shortcut for bad library sizes (too large or too small)
  double currLibSizeScalingFactor = -1;
  double storedLib = -1;
  for(int i = 0; i < this->NUM_LIBS_TO_EST; i++) {
    currLibSizeScalingFactor = gsl_vector_get(probs, this->LIB_SIZE_SCALING_FACTOR_START_IDX + i);
    storedLib = gsl_vector_get(this->initGuessCopyBeforeBFGS, this->LIB_SIZE_SCALING_FACTOR_START_IDX + i);
    if(gsl_isinf(currLibSizeScalingFactor) || gsl_isnan(currLibSizeScalingFactor) || currLibSizeScalingFactor < 1e-2 || currLibSizeScalingFactor > 1e2 || (storedLib - currLibSizeScalingFactor) / storedLib > 0.25) { // libguard can get arbitrarily large but can't get more than .25% smaller
      return GSL_NAN;
    }
  }

  // shortcut for any probabilities becoming too small
  double probMin = gsl_vector_min(probs);
  if(probMin < 1e-5 || gsl_isnan(probMin)) {
    return GSL_NAN;
  }

  // shortcut for any rates becoming too large
  double probMax = gsl_vector_max(probs);
  if(probMax > 10000.0 || gsl_isnan(probMax)) {
    return GSL_NAN;
  }
  return 0;
}

/*
 * function to check if the state of the HMMs is currently ok.
 * That is, returns nan if at least one HMM has a completely
 * uniform transition matrix, 0 o/w
 */
double AllCells3TrParam2DegPolyHMM::checkStateValidity(double epsilon) const {
  HMM* currHMM = nullptr;
  double status = 0;
  for(unsigned int hmmIdx = 0; hmmIdx < this->hmmVec->size(); hmmIdx++) {
    // nullptr guard
    if((*this->hmmVec)[hmmIdx] == nullptr) {
      continue;
    }

    currHMM = (*this->hmmVec)[hmmIdx];
    status = status + currHMM->checkStateValidity(epsilon);
  }
  return status;
}


/*
 * Function to check for transient states
 * (ex state 1 is only seen in 0->1->2 or 2->1->0 transitions)
 * Returns nan if transient states are found, 0 o/w
 */
double AllCells3TrParam2DegPolyHMM::checkForTransientStates() {
  HMM* currHMM = nullptr;
  double status = 0;
  for(unsigned int hmmIdx = 0; hmmIdx < this->hmmVec->size(); hmmIdx++) {
    // nullptr guard
    if((*this->hmmVec)[hmmIdx] == nullptr) {
      continue;
    }

    currHMM = (*this->hmmVec)[hmmIdx];
    status = status + currHMM->checkForTransientStates();
  }
  return status;
}

/*
 * Function to check if initProb has gaps
 * that is, check if the initProb vector is something like [0.15, 0, 0, 0.8, 0, 0, 0.05, 0]
 * ie check if initProb vector has 0 for idx 1 and 2
 */
double AllCells3TrParam2DegPolyHMM::checkForInitProbGaps() {
  HMM* currHMM = nullptr;
  double status = 0;
  for(unsigned int hmmIdx = 0; hmmIdx < this->hmmVec->size(); hmmIdx++) {
    // nullptr guard
    if((*this->hmmVec)[hmmIdx] == nullptr) {
      continue;
    }

    currHMM = (*this->hmmVec)[hmmIdx];
    status = status + currHMM->checkForInitProbGaps();
  }
  return status;
}

void AllCells3TrParam2DegPolyHMM::setCentralDiffFlag(bool flag) {
  this->centralDiff = flag;
  for(unsigned int hmmIdx = 0; hmmIdx < this->hmmVec->size(); hmmIdx++) {
    // nullptr guard
    if((*this->hmmVec)[hmmIdx] == nullptr) {
      continue;
    }
    (*this->hmmVec)[hmmIdx]->setCentralDiffFlag(flag);
  }
}

void AllCells3TrParam2DegPolyHMM::setBaumWelchInitGuess(gsl_vector* initGuess, int numBWIters, int numLibStarts, double libStartValue, bool verbose, bool debug) {
  std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
  std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
  double elapsedSec = 0;
  double totalTime = 0;
  double initTotalLik = 0;
  double bwTotalLik = 0;

  gsl_vector* libMultipliers = gsl_vector_alloc(numLibStarts);
  if(numLibStarts == 2) {
    gsl_vector_set(libMultipliers, 0, 1);
    gsl_vector_set(libMultipliers, 1, 2); // double orig
  }
  else if(numLibStarts == 3) {
    gsl_vector_set(libMultipliers, 0, 1);
    gsl_vector_set(libMultipliers, 1, 2); // double orig
    gsl_vector_set(libMultipliers, 2, 4);
  }
  else if(numLibStarts == 4) {
    gsl_vector_set(libMultipliers, 0, 1);
    gsl_vector_set(libMultipliers, 1, 2); // double orig
    gsl_vector_set(libMultipliers, 2, 4);
    gsl_vector_set(libMultipliers, 3, 2.0/3.0);
  }
  else {
    gsl_vector_set_all(libMultipliers, libStartValue);
  }
  this->baumWelchParamResults = new std::vector<gsl_vector*>();

  gsl_matrix* bwLibResults = gsl_matrix_alloc(numLibStarts + 1 + 1, this->NUM_LIBS_TO_EST); // record results from each bw run. rows: runs, cols: cells; add one for emergency anyInitProbGapsRescale case and one for emergency anySmallChangeBwLlRescale case
  gsl_matrix_set_zero(bwLibResults);
  std::vector<std::vector<gsl_matrix*>*>* bwTransitionMats = new std::vector<std::vector<gsl_matrix*>*>(); // one tr mat per hmm per libStart
  std::vector<std::vector<gsl_vector*>*>* bwInitProbs = new std::vector<std::vector<gsl_vector*>*>();
  std::vector<bool>* bwHasTransientStates = new std::vector<bool>();
  std::vector<bool>* bwHasInitProbGaps = new std::vector<bool>();
  std::vector<bool>* bwHasSmallChangeBwLl = new std::vector<bool>();
  std::vector<double>* bwTotalLiks = new std::vector<double>();

  gsl_vector* origParamsToEst = gsl_vector_alloc(this->getNumParamsToEst());
  gsl_vector_memcpy(origParamsToEst, this->paramsToEst);
  gsl_vector* origTransitionParams = gsl_vector_alloc(this->NUM_SHARED_TRANSITION_PARAMS_TO_EST + this->NUM_BRANCH_LENGTHS_TO_EST);
  for(unsigned origTrIdx = 0; origTrIdx < origTransitionParams->size; origTrIdx++) {
    gsl_vector_set(origTransitionParams, origTrIdx, gsl_vector_get(origParamsToEst, this->SHARED_TRANSITION_PROB_START_IDX + origTrIdx));
  }
  bool anyInitProbGapsRescale = false;
  bool anySmallChangeBwLlRescale = false;
  for(unsigned int libIdx = 0; libIdx < libMultipliers->size || anyInitProbGapsRescale || anySmallChangeBwLlRescale; libIdx++) {
    // reset to original transition params. libs set below to multiple of prev iter
    this->setAllTransition(origTransitionParams);

    // for each HMM, set the lib start for each cell to the current libStart value
    for(unsigned int hmmIdx = 0, i = 0; hmmIdx < this->hmmVec->size(); hmmIdx++) {
      HMM* currHMM = (*this->hmmVec)[hmmIdx];
      // nullptr guard
      if(currHMM == nullptr) {
        continue;
      }

      for(int cellIdx = 0; cellIdx < currHMM->NUM_CELLS; cellIdx++) {
        // if this is the first one, use the libRatio as the starting point. else, scale the prev iter's final lib
        if(libIdx == 0) {
          double libRatio = currHMM->calcLibScalingFactorsToTotalRatio(cellIdx);
          currHMM->setLibScalingFactor(cellIdx, gsl_vector_get(libMultipliers, libIdx) * libRatio); // scale libRatio by multiplier
          gsl_vector_set(this->paramsToEst, this->LIB_SIZE_SCALING_FACTOR_START_IDX + i, gsl_vector_get(libMultipliers, libIdx) * libRatio);
        }
        else {
          double prevLib = gsl_matrix_get(bwLibResults, 0, hmmIdx * currHMM->NUM_LIBS_TO_EST + cellIdx); // scaling relative to first run
          double currLibMult = 0;
          if(anyInitProbGapsRescale) {
            currLibMult = gsl_vector_get(libMultipliers, libMultipliers->size - 1) * 2; // double the last multiplier
          }
          else if(anySmallChangeBwLlRescale) {
            currLibMult = gsl_vector_get(libMultipliers, 0) * 2.0 / 3.0; // mult first multiplier by 2/3
          }
          else {
            currLibMult = gsl_vector_get(libMultipliers, libIdx);
          }
          currHMM->setLibScalingFactor(cellIdx, currLibMult * prevLib); // scale prevLib by multiplier
          gsl_vector_set(this->paramsToEst, this->LIB_SIZE_SCALING_FACTOR_START_IDX + i, currLibMult * prevLib);
        }
        i++;
      }
    }

    if(verbose) {
      std::lock_guard<std::mutex> lock(Optimizable::mtx); // from https://stackoverflow.com/a/18277334
      double currLibMult = 0;
      if(anyInitProbGapsRescale) {
        currLibMult = gsl_vector_get(libMultipliers, libMultipliers->size - 1) * 2; // double the last multiplier
        std::cout << "anyInitProbGapsRescale flag is true, doubling the last multiplier" << std::endl;
      }
      else if(anySmallChangeBwLlRescale) {
        currLibMult = gsl_vector_get(libMultipliers, 0) * 2.0 / 3.0; // mult first multiplier by 2/3
        std::cout << "anySmallChangeBwLlRescale flag is true, multiplying libRatio by 2/3" << std::endl;
      }
      else {
        currLibMult = gsl_vector_get(libMultipliers, libIdx);
      }
      std::cout << "Starting AllCells3TrParam2DegPolyHMM::setBaumWelchInitGuess, with libMultiplier " << currLibMult << " (RUN " << libIdx + 1 << " OF " << libMultipliers->size << ")" << std::endl;
      this->print(stdout);
    }

    // first call baum welch on each HMM
    initTotalLik = this->getLogLikelihood();

    begin = std::chrono::steady_clock::now();
    for(unsigned int hmmIdx = 0; hmmIdx < this->hmmVec->size(); hmmIdx++) {
      HMM* currHMM = (*this->hmmVec)[hmmIdx];
      // nullptr guard
      if(currHMM == nullptr) {
        continue;
      }

      if(verbose) {
        double currLibMult = 0;
        if(anyInitProbGapsRescale) {
          currLibMult = gsl_vector_get(libMultipliers, libMultipliers->size - 1) * 2; // double the last multiplier
        }
        else if(anySmallChangeBwLlRescale) {
          currLibMult = gsl_vector_get(libMultipliers, 0) * 2.0 / 3.0; // mult first multiplier by 2/3
        }
        else {
          currLibMult = gsl_vector_get(libMultipliers, libIdx);
        }
        std::lock_guard<std::mutex> lock(Optimizable::mtx); // from https://stackoverflow.com/a/18277334
        std::cout << "STARTING BAUM WELCH for HMM " << hmmIdx << ", WITH LIBMULTIPLIER " << currLibMult << " (RUN " << libIdx + 1 << " OF " << libMultipliers->size << ")" << std::endl;
      }
      if(debug) {
        std::lock_guard<std::mutex> lock(Optimizable::mtx); // from https://stackoverflow.com/a/18277334
        std::cerr << "##################################################################" << std::endl; // add a buffer line to the err file before each baum welch
      }
      currHMM->runBaumWelch(numBWIters, verbose, debug);
      currHMM->setUpBaumWelchLeastSquares();

      // set lib size from bw into this->paramsToEst
      for(int cellIdx = 0; cellIdx < currHMM->NUM_LIBS_TO_EST; cellIdx++) {
        double currLib = currHMM->getLibScalingFactor(cellIdx);
        gsl_vector_set(this->paramsToEst, this->LIB_SIZE_SCALING_FACTOR_START_IDX + (hmmIdx * currHMM->NUM_LIBS_TO_EST) + cellIdx, currLib);
        gsl_matrix_set(bwLibResults, libIdx, (hmmIdx * currHMM->NUM_LIBS_TO_EST) + cellIdx, currLib);
      }
    }
    bwTotalLik = this->getLogLikelihood();
    // rule out any runs with transient states
    double hasTransientStates = this->checkForTransientStates();
    double hasInitProbGaps = this->checkForInitProbGaps();

    bwTransitionMats->push_back(new std::vector<gsl_matrix*>());
    bwInitProbs->push_back(new std::vector<gsl_vector*>());
    for(unsigned int hmmIdx = 0; hmmIdx < this->hmmVec->size(); hmmIdx++) {
      HMM* currHMM = (*this->hmmVec)[hmmIdx];
      // nullptr guard
      if(currHMM == nullptr) {
        continue;
      }

      gsl_matrix* currTransMat = currHMM->getTransition();
      gsl_matrix* transitionMatToSave = gsl_matrix_alloc(currTransMat->size1, currTransMat->size2);
      gsl_matrix_memcpy(transitionMatToSave, currTransMat);
      (*bwTransitionMats)[libIdx]->push_back(transitionMatToSave);

      gsl_vector* currInitProb = currHMM->getInitProb();
      gsl_vector* initProbToSave = gsl_vector_alloc(currInitProb->size);
      gsl_vector_memcpy(initProbToSave, currInitProb);
      (*bwInitProbs)[libIdx]->push_back(initProbToSave);
    }
    if(!anyInitProbGapsRescale && anySmallChangeBwLlRescale) {
      // be lenient on special case of 2/3
      bwHasTransientStates->push_back(false);
      bwHasInitProbGaps->push_back(gsl_isnan(hasInitProbGaps));
    }
    else {
      bwHasTransientStates->push_back(gsl_isnan(hasTransientStates));
      bwHasInitProbGaps->push_back(gsl_isnan(hasInitProbGaps));
    }
    bwHasSmallChangeBwLl->push_back((bwTotalLik - initTotalLik) < 4000); // this cutoff is arbitrary, and based on empirical results from simulations
    bwTotalLiks->push_back(bwTotalLik);

    end = std::chrono::steady_clock::now();
    elapsedSec = std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() / 1000000.0;
    totalTime = elapsedSec;
    if(verbose) {
      std::lock_guard<std::mutex> lock(Optimizable::mtx); // from https://stackoverflow.com/a/18277334
      printf("BAUM WELCH TOTAL LOGLIKELIHOOD FOUND: %.40f\n", bwTotalLik);
      printf("CHANGE IN BAUM WELCH TOTAL LIKELIHOOD: %.40f\n", bwTotalLik - initTotalLik);
      printf("BAUM WELCH TOTAL TIME (sec): %.10f\n", elapsedSec);
      printf("DONE WITH BAUM WELCH TOTAL (AllCells3TrParam2DegPolyHMM):\n\n");
      this->print(stdout);
    }
    if(debug) {
      std::lock_guard<std::mutex> lock(Optimizable::mtx); // from https://stackoverflow.com/a/18277334
      std::cerr << "##################################################################" << std::endl; // add a buffer line to the err file before starting least squares bfgs
    }

    // check if all libs have transient states if not currently already doing an emergency rerun (ie only do a rescale once)
    if(anyInitProbGapsRescale) {
      anyInitProbGapsRescale = false;
    }
    else if(anySmallChangeBwLlRescale) {
      anySmallChangeBwLlRescale = false;
    }
    else if(libMultipliers->size > 1 && libIdx == libMultipliers->size - 1) { // and on the last planned iter
      anyInitProbGapsRescale = std::any_of(bwHasInitProbGaps->begin(), bwHasInitProbGaps->end(), [](bool v) { return v; }); // https://stackoverflow.com/a/26896132
      anySmallChangeBwLlRescale = (*bwHasSmallChangeBwLl)[0]; // limit to first run
    }
  } // end iterating over libMultipliers


    // figure out which libMultiplier result to send
    int bestLibIdx = -1;
    int secondBestLibIdx = -1;
    double bestBwLik = -std::numeric_limits<double>::max();
    double secondBestBwLik = -std::numeric_limits<double>::max(); // not strictly second best by ll, could be higher in ll but have transient states/be suspicious

    { // put one extra scope around this block for mutex https://stackoverflow.com/questions/54399719/including-stdlock-guard-in-extra-scope
    std::lock_guard<std::mutex> lock(Optimizable::mtx); // from https://stackoverflow.com/a/18277334
    if(numLibStarts > 1) {
      bool allHaveInitProbGaps = std::all_of(bwHasInitProbGaps->begin(), bwHasInitProbGaps->end(), [](bool v) { return v; }); // https://stackoverflow.com/a/26896132

      for(unsigned int bwIdx = 0; bwIdx < bwTotalLiks->size(); bwIdx++) {
        // having gaps in initProb disqualifies, unless all have gaps
        std::cout << "(*bwHasInitProbGaps)[bwIdx]: " << (*bwHasInitProbGaps)[bwIdx] << std::endl;
        bool hasInitProbGaps = (*bwHasInitProbGaps)[bwIdx];
        if(!allHaveInitProbGaps && hasInitProbGaps) {
          continue;
        }
        bwTotalLik = (*bwTotalLiks)[bwIdx];
        // bwTotalLik > bestBwLik > secondBestBwLik
        // bwTotalLik is better than both bests seen so far, this is now the best
        if(bwTotalLik > bestBwLik && bwTotalLik > secondBestBwLik) {
          secondBestBwLik = bestBwLik;
          secondBestLibIdx = bestLibIdx;
          bestBwLik = bwTotalLik;
          bestLibIdx = bwIdx;
        }
        // bestBwLik > bwTotalLik > secondBestBwLik
        // bwTotalLik is in the middle, new second best
        else if(bwTotalLik < bestBwLik && bwTotalLik > secondBestBwLik) {
          secondBestBwLik = bwTotalLik;
          secondBestLibIdx = bwIdx;
        }
        // bestBwLik > secondBestBwLik > bwTotalLik
        // bwTotalLik is worse, skip over
        else {
          continue;
        }
      }
      // if only one run has qualified due to all others having initProbGaps, just use the bestLibIdx
      if(secondBestLibIdx < 0) {
        if(verbose) {
          std::cout << "Run " << bestLibIdx + 1 << " is the only run qualified without any initProbGaps, using run " << bestLibIdx  + 1 << std::endl;
        }
      }
      else {
        bool bestHasTransientStates = (*bwHasTransientStates)[bestLibIdx];
        bool secondBestHasTransientStates = (*bwHasTransientStates)[secondBestLibIdx];

        if(verbose) {
          std::cout << "Run " << bestLibIdx + 1 << " has the best loglikelihood (" << bestBwLik << "), followed by run " << secondBestLibIdx + 1 << " (" << secondBestBwLik << "), diff (" << std::abs(secondBestBwLik - bestBwLik) << ")." << std::endl;
        }
        // if bestHasTransientStates, might get replaced
        int diffInLlThreshold = 100; // threshold for "close" ll
        if(bestHasTransientStates) {
          if(verbose) {
            std::cout << "Best loglik has transient states, ";
          }
          // if second best is nearby or is not suspicious, replace best
          if(std::abs(secondBestBwLik - bestBwLik) < diffInLlThreshold || !secondBestHasTransientStates) {
            if(std::abs(secondBestBwLik - bestBwLik) < diffInLlThreshold) {
              if(verbose) {
                std::cout << "diff is small, ";
              }
            }
            if(!secondBestHasTransientStates) {
              if(verbose) {
                std::cout << "second best does NOT have transient states, ";
              }
            }
            if(verbose) {
              std::cout << "so replacing with second best." << std::endl;
            }
            bestLibIdx = secondBestLibIdx;
            bestBwLik = secondBestBwLik;
            bestHasTransientStates = secondBestHasTransientStates;
          }
          // else, no compelling reason to replace
          else {
            if(verbose) {
              std::cout << "but second best ";
            }
            if(std::abs(secondBestBwLik - bestBwLik) >= diffInLlThreshold) {
              if(verbose) {
                std::cout << "is far, ";
              }
            }
            if(secondBestHasTransientStates) {
              if(verbose) {
                std::cout << "has transient states, ";
              }
            }
            if(verbose) {
              std::cout << "so NOT replacing with second best." << std::endl;
            }
          } // end no compelling reason to replace
        } // end bestHasTransientStates
        // else, just keep best
        else {
          if(verbose) {
            std::cout << "Best loglik does NOT have transient states, using best." << std::endl;
          }
        }
      } // end multiple qualifying runs (secondBestIdx >= 0)
    } // end more than one libStart
    else {
      if(verbose) {
        std::cout << "Only one run, using it." << std::endl;
      }
      bestLibIdx = 0;
      bestBwLik = (*bwTotalLiks)[bestLibIdx];
    }
    } // end extra mutex scope

    std::vector<gsl_matrix*>* bestBwTransitionMats = (*bwTransitionMats)[bestLibIdx];
    std::vector<gsl_vector*>* bestBwInitProbs = (*bwInitProbs)[bestLibIdx];

    // send the best ll transition mats and lib sizes to least squares
    for(unsigned int hmmIdx = 0, i = 0; hmmIdx < this->hmmVec->size(); hmmIdx++) {
      HMM* currHMM = (*this->hmmVec)[hmmIdx];
      // nullptr guard
      if(currHMM == nullptr) {
        continue;
      }

      gsl_matrix* savedTransMat = (*bestBwTransitionMats)[hmmIdx];
      currHMM->setTransition(savedTransMat);
      gsl_vector* savedInitProb = (*bestBwInitProbs)[hmmIdx];
      currHMM->setInitProb(savedInitProb);
      for(int cellIdx = 0; cellIdx < currHMM->NUM_CELLS; cellIdx++) {
        double storedLib = gsl_matrix_get(bwLibResults, bestLibIdx, hmmIdx * currHMM->NUM_LIBS_TO_EST + cellIdx);
        currHMM->setLibScalingFactor(cellIdx, storedLib);
        gsl_vector_set(this->paramsToEst, this->LIB_SIZE_SCALING_FACTOR_START_IDX + i, storedLib);
        i++;
      }
    }
    if(verbose) {
      std::lock_guard<std::mutex> lock(Optimizable::mtx); // from https://stackoverflow.com/a/18277334
      std::cout << "RESETTING TO BEST TRANSITION MATRICES AND LIB SIZE SCALING FACTORS, FROM RUN " << bestLibIdx + 1 << ":" << std::endl;
    }

    // then solve overdetermined system of equations for each HMM
    double leastSquaresStatus = this->doLeastSquares(initGuess, 3, verbose, debug);

    double lsTotalLik = this->getLogLikelihood();

    // if all least squares attempts failed, keep lib size and original rate params
    if(gsl_isnan(leastSquaresStatus) || compareDoubles(0, lsTotalLik - bestBwLik)) {
      if(verbose) {
        std::lock_guard<std::mutex> lock(Optimizable::mtx); // from https://stackoverflow.com/a/18277334
        std::cout << "WARNING: All least squares attempts failed, no change from Baum Welch loglik. Using best library scaling factor and original rate params" << std::endl;
        if(debug) {
          std::cerr << "##################################################################" << std::endl; // add a buffer line to the err file after ending least squares bfgs
        }
      }
      gsl_vector* origParamsToEstCopy = gsl_vector_alloc(origParamsToEst->size);
      gsl_vector_memcpy(origParamsToEstCopy, origParamsToEst);


      for(unsigned int hmmIdx = 0, i = 0; hmmIdx < this->hmmVec->size(); hmmIdx++) {
        HMM* currHMM = (*this->hmmVec)[hmmIdx];
        // nullptr guard
        if(currHMM == nullptr) {
          continue;
        }

        for(int cellIdx = 0; cellIdx < currHMM->NUM_CELLS; cellIdx++) {
          double bestLibFromBw = gsl_matrix_get(bwLibResults, bestLibIdx, hmmIdx * currHMM->NUM_LIBS_TO_EST + cellIdx);
          gsl_vector_set(origParamsToEstCopy, this->LIB_SIZE_SCALING_FACTOR_START_IDX + i, bestLibFromBw);
          i++;
        }
      }
      this->baumWelchParamResults->push_back(origParamsToEstCopy);
      this->setParamsToEst(origParamsToEstCopy);
    }

    end = std::chrono::steady_clock::now();
    elapsedSec = std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() / 1000000.0;
    totalTime += elapsedSec;
    if(verbose) {
      std::lock_guard<std::mutex> lock(Optimizable::mtx); // from https://stackoverflow.com/a/18277334
      printf("CHANGE IN LEAST SQUARES LIKELIHOOD FROM STARTING POINT: %.40f\n", lsTotalLik - initTotalLik);
      printf("CHANGE IN LEAST SQUARES LIKELIHOOD FROM BAUM WELCH: %.40f\n", lsTotalLik - bestBwLik);
      printf("BAUM WELCH WITH LEAST SQUARES TOTAL TIME (sec): %.10f\n", totalTime);
      printf("DONE WITH AllCells3TrParam2DegPolyHMM BAUM WELCH LEAST SQUARES\n");
    }

    // if optim failed (as determined by a uniform transition matrix), try again. else, save params
    double allValidTransitionMatrices = this->checkStateValidity();
    if(gsl_isnan(allValidTransitionMatrices)) {
      if(verbose) {
        std::lock_guard<std::mutex> lock(Optimizable::mtx); // from https://stackoverflow.com/a/18277334
        std::cout << "WARNING: at least one HMM had an invalid transition matrix from this baum welch + least squares parameter set:" << std::endl;
        this->print(stdout);
      }
    }
    gsl_vector* paramsToSave = gsl_vector_alloc(initGuess->size);
    gsl_vector_memcpy(paramsToSave, initGuess);
    this->baumWelchParamResults->push_back(paramsToSave);

    if(debug) {
      std::lock_guard<std::mutex> lock(Optimizable::mtx); // from https://stackoverflow.com/a/18277334
      std::cerr << "##################################################################" << std::endl; // add a buffer line to the err file after ending least squares bfgs
    }

  if(this->baumWelchParamResults->size() == 0) {
    if(verbose) {
      std::lock_guard<std::mutex> lock(Optimizable::mtx); // from https://stackoverflow.com/a/18277334
      std::cout << "WARNING: All baum welch + least squares runs failed. Using original params" << std::endl;
      if(debug) {
        std::cerr << "##################################################################" << std::endl; // add a buffer line to the err file after ending least squares bfgs
      }
    }
    gsl_vector* origParamsToEstCopy = gsl_vector_alloc(origParamsToEst->size);
    gsl_vector_memcpy(origParamsToEstCopy, origParamsToEst);
    this->baumWelchParamResults->push_back(origParamsToEstCopy);
    this->setParamsToEst(origParamsToEstCopy);
  }
  if(verbose) {
    std::lock_guard<std::mutex> lock(Optimizable::mtx); // from https://stackoverflow.com/a/18277334
    std::cout << "The baum welch parameter sets are:" << std::endl;
    for(unsigned int i = 0; i < this->baumWelchParamResults->size(); i++) {
      if((*this->baumWelchParamResults)[i] != nullptr) {
        printColVector((*this->baumWelchParamResults)[i]);
      }
    }
  }
  gsl_vector_memcpy(initGuess, (*this->baumWelchParamResults)[0]); // save first set of params into initGuess

  // clean up
  bestBwTransitionMats->clear();
  bestBwInitProbs->clear();
  gsl_vector_free(libMultipliers);
  gsl_vector_free(origParamsToEst);
  gsl_matrix_free(bwLibResults);
}
double AllCells3TrParam2DegPolyHMM::doLeastSquares(gsl_vector* initGuess, int maxNumAttempts, bool verbose, bool debug) {
  std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
  double overallStatus = GSL_NAN;
  double status = GSL_CONTINUE;
  double initSumSqResid = 0;
  double currSumSqResid = 0;
  double prevSumSqResid = 0;
  double deltaSumSqResid = 0;
  double totalSumSqResid = 0;
  double finalSumSqResid = 0;
  int iter = 0;
  int countTooClose = 0;
  std::chrono::steady_clock::time_point lsBegin;
  std::chrono::steady_clock::time_point lsEnd;
  double lsElapsedTime = 0;

  /*
   * structure of varsToEst
   * [b]
   * [L]
   * [pair0, t1]
   * [pair0, t2]
   * [pair0, t3]
   * ...
   * [pairN, t3]
   */
  gsl_vector* varsToEst = gsl_vector_alloc(this->NUM_SHARED_TRANSITION_PARAMS_TO_EST + this->NUM_BRANCH_LENGTHS_TO_EST);
  gsl_vector* varsToEst_probSpace = gsl_vector_alloc(varsToEst->size); // for printing in prob space

  for(int lsAttempt = 1; lsAttempt <= maxNumAttempts; lsAttempt++) {
    if(verbose) {
      std::lock_guard<std::mutex> lock(Optimizable::mtx); // from https://stackoverflow.com/a/18277334
      std::cout << "STARTING LEAST SQUARES, (ATTEMPT " << lsAttempt << " OF " << maxNumAttempts << ")" << std::endl;
    }

    if(lsAttempt == 1) {
      gsl_vector_set(varsToEst_probSpace, 0, 0.01); // beta
      gsl_vector_set(varsToEst_probSpace, 1, 500); // lambda
      if(this->NUM_BRANCH_LENGTHS_TO_EST == 1) {
        gsl_vector_set(varsToEst_probSpace, 2, 0.01); // t (for single independent cells)
      }
      else {
        for(int pairIdx = 0; pairIdx < this->NUM_BRANCH_LENGTHS_TO_EST / 3; pairIdx++) {
          gsl_vector_set(varsToEst_probSpace, 2 + 3 * pairIdx + 0, 0.01); // t1
          gsl_vector_set(varsToEst_probSpace, 2 + 3 * pairIdx + 1, 0.01); // t2
          gsl_vector_set(varsToEst_probSpace, 2 + 3 * pairIdx + 2, 0.01); // t3
        }
      }
    }
    else if(lsAttempt == 2) {
      gsl_vector_set(varsToEst_probSpace, 0, 0.1); // beta
      gsl_vector_set(varsToEst_probSpace, 1, 100); // lambda
      if(this->NUM_BRANCH_LENGTHS_TO_EST == 1) {
        gsl_vector_set(varsToEst_probSpace, 2, 0.001); // t (for single independent cells)
      }
      else {
        for(int pairIdx = 0; pairIdx < this->NUM_BRANCH_LENGTHS_TO_EST / 3; pairIdx++) {
          gsl_vector_set(varsToEst_probSpace, 2 + 3 * pairIdx + 0, 0.001); // t1
          gsl_vector_set(varsToEst_probSpace, 2 + 3 * pairIdx + 1, 0.001); // t2
          gsl_vector_set(varsToEst_probSpace, 2 + 3 * pairIdx + 2, 0.001); // t3
        }
      }
    }
    else {
      gsl_vector_set(varsToEst_probSpace, 0, 0.05); // beta
      gsl_vector_set(varsToEst_probSpace, 1, 1000); // lambda
      if(this->NUM_BRANCH_LENGTHS_TO_EST == 1) {
        gsl_vector_set(varsToEst_probSpace, 2, 0.05); // t (for single independent cells)
      }
      else {
        for(int pairIdx = 0; pairIdx < this->NUM_BRANCH_LENGTHS_TO_EST / 3; pairIdx++) {
          gsl_vector_set(varsToEst_probSpace, 2 + 3 * pairIdx + 0, 0.05); // t1
          gsl_vector_set(varsToEst_probSpace, 2 + 3 * pairIdx + 1, 0.05); // t2
          gsl_vector_set(varsToEst_probSpace, 2 + 3 * pairIdx + 2, 0.05); // t3
        }
      }
    }

    this->baumWelchLeastSquares_convertProbToParam(varsToEst, varsToEst_probSpace);

    gsl_multimin_function_fdf my_func; // which function to minimize
    my_func.df = &baumWelchLeastSquares_df; // gradient of function
    my_func.fdf = &baumWelchLeastSquares_fdf; // how to set both f and df

    my_func.f = &baumWelchLeastSquares_f; // function itself
    my_func.params = this;
    my_func.n = varsToEst->size;

    const gsl_multimin_fdfminimizer_type* T = gsl_multimin_fdfminimizer_vector_bfgs2; // which algorithm to use
    gsl_multimin_fdfminimizer* s = gsl_multimin_fdfminimizer_alloc(T, my_func.n);
    gsl_multimin_fdfminimizer_set(s, &my_func, varsToEst, 1e-1, .1); // minimizer s, function my_func, starting at x, first step_size, accuracy of line minimization specified by tol(erance

    status = GSL_CONTINUE;
    currSumSqResid = 0;
    deltaSumSqResid = 0;
    totalSumSqResid = 0;
    iter = 0;
    countTooClose = 0;
    lsElapsedTime = 0;

    initSumSqResid = baumWelchLeastSquares_f(varsToEst, this);
    prevSumSqResid = initSumSqResid;
    if(verbose) {
      std::lock_guard<std::mutex> lock(Optimizable::mtx); // from https://stackoverflow.com/a/18277334
      printf("BW LEAST SQUARES BFGS INITIAL currSumSqResid: %.20f\n", initSumSqResid);
    }

    while(status == GSL_CONTINUE) {
      lsBegin = std::chrono::steady_clock::now();

      status = gsl_multimin_fdfminimizer_iterate(s); // do one iter
      currSumSqResid = gsl_multimin_fdfminimizer_minimum(s);

      lsEnd = std::chrono::steady_clock::now();
      lsElapsedTime = std::chrono::duration_cast<std::chrono::microseconds>(lsEnd - lsBegin).count() / 1000000.0;
      deltaSumSqResid = std::abs(prevSumSqResid - currSumSqResid);
      totalSumSqResid += deltaSumSqResid;
      prevSumSqResid = currSumSqResid;

      if(verbose) {
        std::lock_guard<std::mutex> lock(Optimizable::mtx); // from https://stackoverflow.com/a/18277334
        printf("ON BW LEAST SQUARES BFGS ITER %d, currSumSqResid %.20f, time elapsed (sec) %.5f, total change in currSumSqResid %.5f\n", iter, currSumSqResid, lsElapsedTime, totalSumSqResid);
      }

      if(status) {
        if(verbose) {
          std::lock_guard<std::mutex> lock(Optimizable::mtx); // from https://stackoverflow.com/a/18277334
          std::cout << "STATUS IS: " << gsl_strerror(status) << std::endl;
        }
        this->baumWelchLeastSquares_convertParamToProb(varsToEst_probSpace, gsl_multimin_fdfminimizer_x(s));
        printRowVector(varsToEst_probSpace);
        break;
      }
      status = gsl_multimin_test_gradient(s->gradient, 1e-8);

      if(status == GSL_SUCCESS) {
        if(verbose) {
          std::lock_guard<std::mutex> lock(Optimizable::mtx); // from https://stackoverflow.com/a/18277334
          printf("SUCCESS: BW+LS minimum found at: ");
        }
        this->baumWelchLeastSquares_convertParamToProb(varsToEst_probSpace, gsl_multimin_fdfminimizer_x(s));
        printRowVector(varsToEst_probSpace);
      }
      if(iter >= 500) {
        this->baumWelchLeastSquares_convertParamToProb(varsToEst_probSpace, gsl_multimin_fdfminimizer_x(s)); // save anyway
        if(verbose) {
          std::lock_guard<std::mutex> lock(Optimizable::mtx); // from https://stackoverflow.com/a/18277334
          std::cout << "STATUS IS: minimum not found in 500 iterations" << std::endl;
          printRowVector(varsToEst_probSpace);
        }
        break;
      }

      if(deltaSumSqResid < 1e-5) {
        countTooClose++;
        if(countTooClose >= 3) {
          this->baumWelchLeastSquares_convertParamToProb(varsToEst_probSpace, gsl_multimin_fdfminimizer_x(s)); // save anyway
          if(verbose) {
            std::lock_guard<std::mutex> lock(Optimizable::mtx); // from https://stackoverflow.com/a/18277334
            std::cout << "STATUS IS: converged by consecutive small change in sumSqResid" << std::endl;
            printRowVector(varsToEst_probSpace);
          }
          break;
        }
      } else {
        countTooClose = 0;
      }

      iter++;
    }
    gsl_multimin_fdfminimizer_free(s);

    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    double elapsedSec = std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() / 1000000.0;
    if(verbose) {
      std::lock_guard<std::mutex> lock(Optimizable::mtx); // from https://stackoverflow.com/a/18277334
      printf("LEAST SQUARES TIME (sec): %.10f\n", elapsedSec);
    }

    // if ls succeeded (>0 change in currSumSqResid and no uniform matrix), then break. else, try again
    finalSumSqResid = currSumSqResid;
    deltaSumSqResid = std::abs(initSumSqResid - finalSumSqResid);
    double allValidTransitionMatrices = this->checkStateValidity();
    if(deltaSumSqResid > 0 && !gsl_isnan(allValidTransitionMatrices)) {
      if(verbose) {
        std::lock_guard<std::mutex> lock(Optimizable::mtx); // from https://stackoverflow.com/a/18277334
        std::cout << "LEAST SQUARES (ATTEMPT " << lsAttempt << " OF " << maxNumAttempts << ") SUCCEEDED, BREAKING\n" << std::endl;
      }
      overallStatus = 0;

      // save results. delegate to subclass since depends on number of cells to save libs
      this->saveBaumWelchEstsIntoParamsToEst(varsToEst_probSpace, initGuess);
      double lsTotalLik = this->getLogLikelihood();
      if(verbose) {
        std::lock_guard<std::mutex> lock(Optimizable::mtx); // from https://stackoverflow.com/a/18277334
        printf("LEAST SQUARES LOGLIKELIHOOD FOUND: %.40f\n", lsTotalLik);
      }

      break;
    }
    else {
      if(verbose) {
        std::lock_guard<std::mutex> lock(Optimizable::mtx); // from https://stackoverflow.com/a/18277334
        std::cout << "LEAST SQUARES (ATTEMPT " << lsAttempt << " OF " << maxNumAttempts << ") FAILED due to ";
        if(gsl_isnan(allValidTransitionMatrices)) {
          std::cout << "at least one invalid transition matrix";
        }
        else {
          std::cout << "no change in currSumSqResid";
        }
        std::cout << ". Trying again\n" << std::endl;
      }
    }
  }

  return overallStatus;
}
double AllCells3TrParam2DegPolyHMM::baumWelchLeastSquares_f(const gsl_vector* v, void* params) {
  AllCells3TrParam2DegPolyHMM* hmm = (AllCells3TrParam2DegPolyHMM*) params;
  double sumSqResid = hmm->baumWelchLeastSquares_calcSumSqResid(v); // delegate to subclass, since depends on param ordering
  return sumSqResid;
}
void AllCells3TrParam2DegPolyHMM::baumWelchLeastSquares_df(const gsl_vector* v, void* params, gsl_vector* df) {
  // approx using finite difference method:
  // f'(x) = [f(x + h) - f(x)] / h
  double h = 1e-4;
  double f_x_ph = 0; // f(x+h)
  double f_x_mh = 0; // f(x-h)
  double x = 0;
  double derivApprox = 0;
  gsl_vector_set_zero(df);
  gsl_vector* currVec = ((AllCells3TrParam2DegPolyHMM*)params)->bwGradientIntermediateVec;
  gsl_vector_set_zero(currVec);
  for(unsigned int i = 0; i < v->size; i++) {
    gsl_vector_memcpy(currVec, v);
    x = gsl_vector_get(v, i);
    gsl_vector_set(currVec, i, x + h);
    f_x_ph = AllCells3TrParam2DegPolyHMM::baumWelchLeastSquares_f(currVec, params);

    if(gsl_isnan(f_x_ph)) {
      derivApprox = GSL_NAN;
    } else {
      //derivApprox = (f_x_ph - f_x) / h; // forward difference

      gsl_vector_set(currVec, i, x - h);
      f_x_mh = AllCells3TrParam2DegPolyHMM::baumWelchLeastSquares_f(currVec, params);
      if(gsl_isnan(f_x_mh)) {
        derivApprox = GSL_NAN;
      } else {
        derivApprox = (f_x_ph - f_x_mh) / (2*h); // central difference*/
      }
    }
    gsl_vector_set(df, i, derivApprox);
  }

  // BFGS printing
  if(((AllCells3TrParam2DegPolyHMM*) params)->gradientDebug) {
    std::lock_guard<std::mutex> lock(Optimizable::mtx); // from https://stackoverflow.com/a/18277334
    gsl_vector* probs = gsl_vector_alloc(v->size);
    ((AllCells3TrParam2DegPolyHMM*) params)->baumWelchLeastSquares_convertParamToProb(probs, v);
    fprintf(stderr, "%.20f \t", AllCells3TrParam2DegPolyHMM::baumWelchLeastSquares_f(v, params));
    fprintf(stderr, "|\t");
    for(unsigned int idx = 0; idx < probs->size; idx++) {
      fprintf(stderr, "%.20f\t", gsl_vector_get(probs, idx));
    }
    fprintf(stderr, "|\t");
    for(unsigned int idx = 0; idx < v->size; idx++) {
      fprintf(stderr, "%.20f\t", gsl_vector_get(v, idx));
    }
    fprintf(stderr, "|\t");
    for(unsigned int idx = 0; idx < df->size; idx++) {
      fprintf(stderr, "%.20f\t", gsl_vector_get(df, idx));
    }
    fprintf(stderr, "\n");
    gsl_vector_free(probs);
  }
}
void AllCells3TrParam2DegPolyHMM::baumWelchLeastSquares_fdf(const gsl_vector* v, void* params, double* f, gsl_vector* df) {
  *f = baumWelchLeastSquares_f(v, params);
  baumWelchLeastSquares_df(v, params, df);
}

AllCells3TrParam2DegPolyHMM* AllCells3TrParam2DegPolyHMM::bfgs(gsl_vector* initGuess, int maxIters, bool verbose, bool debug) {
  AllCells3TrParam2DegPolyHMM* bestGuessOptim = this;
  gsl_vector* initGuessAsParams = gsl_vector_alloc(initGuess->size);
  this->convertProbToParam(initGuessAsParams, initGuess);
  Optimizable::bfgs(initGuessAsParams, bestGuessOptim, maxIters, verbose, debug);
  return bestGuessOptim;
}

double AllCells3TrParam2DegPolyHMM::getLogLikelihood() {
  double totalLogLikelihood = 0;
  double currLL = 0;
  for(std::vector<HMM*>::iterator itr = this->hmmVec->begin(); itr != this->hmmVec->end(); ++itr) {
    // nullptr guard
    if((*itr) == nullptr) {
      continue;
    }
    // if finalLl has been set, then we don't need to rerun forward alg
    if((*itr)->finalLl != -std::numeric_limits<double>::infinity()) {
      currLL = (*itr)->finalLl;
    }
    else {
      currLL = (*itr)->getLogLikelihood();
    }
    totalLogLikelihood += currLL;
  }
  return totalLogLikelihood;
}

double AllCells3TrParam2DegPolyHMM::getViterbiLogLikelihood() {
  double vitLL = 0;
  double currLL = 0;
  for(std::vector<HMM*>::iterator itr = this->hmmVec->begin(); itr != this->hmmVec->end(); ++itr) {
    // nullptr guard
    if((*itr) == nullptr) {
      continue;
    }
    currLL = (*itr)->getViterbiLogLikelihood();
    vitLL += currLL;
  }
  return vitLL;
}

/*
 * call each HMM's viterbiDecode
 */
void AllCells3TrParam2DegPolyHMM::viterbiDecodeAll() {
  for(std::vector<HMM*>::iterator itr = this->hmmVec->begin(); itr != this->hmmVec->end(); ++itr) {
    // nullptr guard
    if((*itr) == nullptr) {
      continue;
    }
    (*itr)->allocDecodingIntermediates(); // includes an internal check if it's already been called
    (*itr)->viterbiDecode();
  }
}

/*
 * call each HMM's saveViterbiDecodedCNA method. Each HMM's output is numbered from 0
 */
void AllCells3TrParam2DegPolyHMM::saveAllViterbiDecodedCNA(std::string filename) {
  for(unsigned int hmmIdx = 0; hmmIdx < this->hmmNames->size(); hmmIdx++) {
    std::string hmmName = (*this->hmmNames)[hmmIdx];
    boost::replace_all(hmmName, ",", "__");
    (*this->hmmVec)[hmmIdx]->saveViterbiDecodedCNA(filename + "__" + hmmName + "__k" + std::to_string(this->MAX_PLOIDY) + ".viterbiDecoded");
  }

  std::vector<gsl_vector*>* BFGSParamResults = this->getBFGSParamResults();
  if(BFGSParamResults != nullptr) {
    gsl_vector* bestParams = gsl_vector_alloc(this->getNumParamsToEst());
    gsl_vector_memcpy(bestParams, this->getParamsToEst());

    for(unsigned int i = 1; i <= BFGSParamResults->size(); i++) {
      this->setParamsToEst((*BFGSParamResults)[i - 1]);
      std::cout << "SAVING RESULTS FOR i: " << i << std::endl;
      this->print(stdout);
      for(unsigned int hmmIdx = 0; hmmIdx < this->hmmVec->size(); hmmIdx++) {
        std::string hmmName = (*this->hmmNames)[hmmIdx];
        boost::replace_all(hmmName, ",", "__");
        // redo viterbi decoding under these params before saving
        (*this->hmmVec)[hmmIdx]->viterbiDecode();
        (*this->hmmVec)[hmmIdx]->saveViterbiDecodedCNA(filename + "__" + hmmName + "__k" + std::to_string(this->MAX_PLOIDY) + "__run" + std::to_string(i) + ".viterbiDecoded");
      }
    }
    // reset to previously best params
    this->setParamsToEst(bestParams);
    for(unsigned int hmmIdx = 0; hmmIdx < this->hmmVec->size(); hmmIdx++) {
      // nullptr guard
      if((*this->hmmVec)[hmmIdx] == nullptr) {
        continue;
      }

      (*this->hmmVec)[hmmIdx]->viterbiDecode();
    }
    gsl_vector_free(bestParams);
  }
}

/*
 * returns if a given hmm has been read from file successfully
 */
bool AllCells3TrParam2DegPolyHMM::hasIthHMMBeenReadFromFile(int hmmIdx) {
  return (*this->hmmVec)[hmmIdx]->hasBeenReadFromFile;
}

gsl_vector* AllCells3TrParam2DegPolyHMM::getParamsToEstFromIthHMM(int hmmIdx) {
  return (*this->hmmVec)[hmmIdx]->getParamsToEst();
}

