#include "TwoCell0TrParam2DegPolyHMM.hpp"

/*
 ********
 * constructors and destructor
 ********
 */
TwoCell0TrParam2DegPolyHMM::TwoCell0TrParam2DegPolyHMM(std::vector<DepthPair*>* depths, gsl_vector* fixedParams, int maxPloidy, bool preallocIntermediates) : TwoCell0TrParam2DegPolyHMM(depths, fixedParams, maxPloidy, 0, 3, 0, 3, preallocIntermediates) { // 0 transition params to est, 3 fixedTrParams, 0 fixedLibs, 3 branches
}
TwoCell0TrParam2DegPolyHMM::TwoCell0TrParam2DegPolyHMM(std::vector<DepthPair*>* depths, gsl_vector* fixedParams, int maxPloidy, int numTrParamsToEst, int numFixedTrParams, int numFixedLibs, int numBranches, bool preallocIntermediates) : TwoCell3TrParam2DegPolyHMM(depths, fixedParams, maxPloidy, numTrParamsToEst, numFixedTrParams, numFixedLibs, numBranches, preallocIntermediates) {
  this->maxNumBFGSStarts = 3;
}

//TwoCell0TrParam2DegPolyHMM::TwoCell0TrParam2DegPolyHMM(const TwoCell0TrParam2DegPolyHMM& otherHMM) : TwoCell3TrParam2DegPolyHMM(otherHMM) {
//}
TwoCell0TrParam2DegPolyHMM::~TwoCell0TrParam2DegPolyHMM() {
}

//int TwoCell0TrParam2DegPolyHMM::getMaxNumBFGSStarts() const {
//  return this->maxNumBFGSStarts;
//}
/*
 ********
 * accessors and mutators
 ********
 */
/*
 ********
 * functions that depend on numbering and ordering of transition params
 ********
 */
//double TwoCell0TrParam2DegPolyHMM::setTransition(gsl_vector* transitionParams) {
double TwoCell0TrParam2DegPolyHMM::setTransition(gsl_matrix* dest, gsl_vector* transitionParams) {
  //double alpha = gsl_vector_get(this->fixedParams, 0);
  double beta  = gsl_vector_get(this->fixedParams, this->FIXED_TRANSITION_PROB_START_IDX + 1);
  //double gamma = gsl_vector_get(this->fixedParams, this->FIXED_TRANSITION_PROB_START_IDX + 2);
  double lambda = gsl_vector_get(this->fixedParams, this->FIXED_TRANSITION_PROB_START_IDX + 2);

  double t1 = 0;
  double t2 = 0;
  double t3 = 0;

  // if given new params to use, use those. Otherwise, we're changing values from fixedParams and want an easy way to call setTransition
  if(transitionParams != nullptr) {
    t1 = gsl_vector_get(transitionParams, 0);
    t2 = gsl_vector_get(transitionParams, 1);
    t3 = gsl_vector_get(transitionParams, 2);

    // save into paramsToEst
    int transitionParamsIdx = 0;
    for(int i = 0; i < this->NUM_BRANCH_LENGTHS_TO_EST; i++) {
      gsl_vector_set(this->paramsToEst, this->BRANCH_LENGTH_START_IDX + i, gsl_vector_get(transitionParams, transitionParamsIdx));
      transitionParamsIdx++;
    }
  }
  else {
    t1 = gsl_vector_get(this->paramsToEst, this->BRANCH_LENGTH_START_IDX + 0);
    t2 = gsl_vector_get(this->paramsToEst, this->BRANCH_LENGTH_START_IDX + 1);
    t3 = gsl_vector_get(this->paramsToEst, this->BRANCH_LENGTH_START_IDX + 2);
  }

  //return TwoCell3TrParam2DegPolyHMM::setTransition(dest, this->getAlpha(), beta, gamma, t1, t2, t3); // call parent class's setTransition helper method
  return TwoCell3TrParam2DegPolyHMM::setTransition(dest, this->getAlpha(), beta, lambda, t1, t2, t3); // call parent class's setTransition helper method
}

void TwoCell0TrParam2DegPolyHMM::convertProbToParam(gsl_vector* dest, const gsl_vector* src) const {
  //double d = (double) (*this->depths)[0]->maxWindowSize;
  //double t2 = gsl_vector_get(src, this->TRANSITION_PROB_START_IDX + 3) / d;
  //double t3 = gsl_vector_get(src, this->TRANSITION_PROB_START_IDX + 4) / d;
  //double t1 = gsl_vector_get(src, this->TRANSITION_PROB_START_IDX + 0) / d;
  //double t2 = gsl_vector_get(src, this->TRANSITION_PROB_START_IDX + 1) / d;
  //double t3 = gsl_vector_get(src, this->TRANSITION_PROB_START_IDX + 2) / d;
  double t1 = gsl_vector_get(src, this->TRANSITION_PROB_START_IDX + 0);
  double t2 = gsl_vector_get(src, this->TRANSITION_PROB_START_IDX + 1);
  double t3 = gsl_vector_get(src, this->TRANSITION_PROB_START_IDX + 2);
  double r = gsl_vector_get(src, this->LIB_SIZE_SCALING_FACTOR_START_IDX);
  double s = gsl_vector_get(src, this->LIB_SIZE_SCALING_FACTOR_START_IDX + 1);

  //gsl_vector_set(dest, this->TRANSITION_PROB_START_IDX + 3, log(-(d * t2) / (d * (t2 + t3) - 1))); // set T2
  //gsl_vector_set(dest, this->TRANSITION_PROB_START_IDX + 4, log(-(d * t3) / (d * (t2 + t3) - 1))); // set T3
  //gsl_vector_set(dest, this->TRANSITION_PROB_START_IDX + 0, log(-(d * t1) / (d * (t1 + t2 + t3) - 1))); // set T1
  //gsl_vector_set(dest, this->TRANSITION_PROB_START_IDX + 1, log(-(d * t2) / (d * (t1 + t2 + t3) - 1))); // set T2
  //gsl_vector_set(dest, this->TRANSITION_PROB_START_IDX + 2, log(-(d * t3) / (d * (t1 + t2 + t3) - 1))); // set T3
  gsl_vector_set(dest, this->TRANSITION_PROB_START_IDX + 0, log(t1)); // set T1
  gsl_vector_set(dest, this->TRANSITION_PROB_START_IDX + 1, log(t2)); // set T2
  gsl_vector_set(dest, this->TRANSITION_PROB_START_IDX + 2, log(t3)); // set T3

  gsl_vector_set(dest, this->LIB_SIZE_SCALING_FACTOR_START_IDX, log(r));
  gsl_vector_set(dest, this->LIB_SIZE_SCALING_FACTOR_START_IDX + 1, log(s));
}

/*
 * reverse of convertProbToParam (see above)
 */
void TwoCell0TrParam2DegPolyHMM::convertParamToProb(gsl_vector* dest, const gsl_vector* src) const {
  //double T2 = gsl_vector_get(src, this->TRANSITION_PROB_START_IDX + 3);
  //double T3 = gsl_vector_get(src, this->TRANSITION_PROB_START_IDX + 4);
  double T1 = gsl_vector_get(src, this->TRANSITION_PROB_START_IDX + 0);
  double T2 = gsl_vector_get(src, this->TRANSITION_PROB_START_IDX + 1);
  double T3 = gsl_vector_get(src, this->TRANSITION_PROB_START_IDX + 2);
  double w = gsl_vector_get(src, this->LIB_SIZE_SCALING_FACTOR_START_IDX);
  double v = gsl_vector_get(src, this->LIB_SIZE_SCALING_FACTOR_START_IDX + 1);
  //double c = 1.0 / (1 + exp(T2) + exp(T3));
  //double c = 1.0 / (1 + exp(T1) + exp(T2) + exp(T3));

  //gsl_vector_set(dest, this->TRANSITION_PROB_START_IDX + 3, exp(T2) * c); // set t2
  //gsl_vector_set(dest, this->TRANSITION_PROB_START_IDX + 4, exp(T3) * c); // set t3
  //gsl_vector_set(dest, this->TRANSITION_PROB_START_IDX + 0, exp(T1) * c); // set t1
  //gsl_vector_set(dest, this->TRANSITION_PROB_START_IDX + 1, exp(T2) * c); // set t2
  //gsl_vector_set(dest, this->TRANSITION_PROB_START_IDX + 2, exp(T3) * c); // set t3
  gsl_vector_set(dest, this->TRANSITION_PROB_START_IDX + 0, exp(T1)); // set t1
  gsl_vector_set(dest, this->TRANSITION_PROB_START_IDX + 1, exp(T2)); // set t2
  gsl_vector_set(dest, this->TRANSITION_PROB_START_IDX + 2, exp(T3)); // set t3

  gsl_vector_set(dest, this->LIB_SIZE_SCALING_FACTOR_START_IDX, exp(w));
  gsl_vector_set(dest, this->LIB_SIZE_SCALING_FACTOR_START_IDX + 1, exp(v));
}

/*
 ********
 * functions
 ********
 */

/*
 * calls bfgs, returns a bestGuessHMM
 */
TwoCell0TrParam2DegPolyHMM* TwoCell0TrParam2DegPolyHMM::bfgs(gsl_vector* initGuess, int maxIters, bool verbose, bool debug) {
  // create new HMM with the best guess parameters and return it
  //TwoCell0TrParam2DegPolyHMM* bestGuessHMM = new TwoCell0TrParam2DegPolyHMM(*this);
  TwoCell0TrParam2DegPolyHMM* bestGuessHMM = this;//new TwoCell0TrParam2DegPolyHMM(*this);
  gsl_vector* initGuessAsParams = gsl_vector_alloc(initGuess->size);
  this->convertProbToParam(initGuessAsParams, initGuess);
  Optimizable::bfgs(initGuessAsParams, bestGuessHMM, maxIters, verbose, debug);
  gsl_vector_free(initGuessAsParams);
  return bestGuessHMM;
}

/*
 * hard coded param sets for spread out starting points taken from AllPairs3TrParam2DegPolyHMM (Mon 02 Dec 2019 02:49:20 PM PST)
 */
void TwoCell0TrParam2DegPolyHMM::setInitGuessNthTime(gsl_vector* initGuess, int iter, int numTotalRuns) const {
  gsl_vector_set_zero(initGuess);
  // TODO write other cases
  if(iter == 0) {
    // lib scaling factors
    for(int cellIdx = 0; cellIdx < this->NUM_LIBS_TO_EST; cellIdx++) {
      gsl_vector_set(initGuess, this->LIB_SIZE_SCALING_FACTOR_START_IDX + cellIdx, 1);
    }

    // pairwise branch lengths
    gsl_vector_set(initGuess, this->BRANCH_LENGTH_START_IDX + 0, 0.2); // set t1
    gsl_vector_set(initGuess, this->BRANCH_LENGTH_START_IDX + 1, 0.2); // set t2
    gsl_vector_set(initGuess, this->BRANCH_LENGTH_START_IDX + 2, 0.2); // set t3
  }
  else if(iter == 1) {
    // lib scaling factors
    for(int cellIdx = 0; cellIdx < this->NUM_LIBS_TO_EST; cellIdx++) {
      gsl_vector_set(initGuess, this->LIB_SIZE_SCALING_FACTOR_START_IDX + cellIdx, .75);
    }

    // pairwise branch lengths
    gsl_vector_set(initGuess, this->BRANCH_LENGTH_START_IDX + 0, 0.1); // set t1
    gsl_vector_set(initGuess, this->BRANCH_LENGTH_START_IDX + 1, 0.1); // set t2
    gsl_vector_set(initGuess, this->BRANCH_LENGTH_START_IDX + 2, 0.1); // set t3
  }
  else if(iter == 2) {
    // lib scaling factors
    for(int cellIdx = 0; cellIdx < this->NUM_LIBS_TO_EST; cellIdx++) {
      gsl_vector_set(initGuess, this->LIB_SIZE_SCALING_FACTOR_START_IDX + cellIdx, 1.25);
    }

    // pairwise branch lengths
    gsl_vector_set(initGuess, this->BRANCH_LENGTH_START_IDX + 0, 0.02); // set t1
    gsl_vector_set(initGuess, this->BRANCH_LENGTH_START_IDX + 1, 0.02); // set t2
    gsl_vector_set(initGuess, this->BRANCH_LENGTH_START_IDX + 2, 0.02); // set t3
  }
}

