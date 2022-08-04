#include "OneCell3TrParam2DegPolyHMM.hpp"

/*
 ********
 * constructors and destructor
 ********
 */
OneCell3TrParam2DegPolyHMM::OneCell3TrParam2DegPolyHMM(std::vector<DepthPair*>* depths, int maxPloidy) : OneCell3TrParam2DegPolyHMM(depths, nullptr, maxPloidy, 2, 1, 0, 1) { // depths, nullptr fixed params, maxPloidy, 2 transition params to est (beta+gamma), 1 fixedTrParams (alpha), 0 numFixedLibs (est all 1 lib), 1 branch to est (t)
}
OneCell3TrParam2DegPolyHMM::OneCell3TrParam2DegPolyHMM(std::vector<DepthPair*>* depths, gsl_vector* fixedParams, int maxPloidy, int numTrParamsToEst, int numFixedTrParams, int numFixedLibs, int numBranches) : HMM (depths, fixedParams, numTrParamsToEst, 1, numBranches, maxPloidy, numFixedTrParams, numFixedLibs) {
  this->logFacKVec = nullptr;
  this->maxNumBFGSStarts = std::numeric_limits<int>::max();
  this->setAlpha(0.1); // arbitrary starting value
}

OneCell3TrParam2DegPolyHMM::~OneCell3TrParam2DegPolyHMM() {
  delete this->logFacKVec;
}

// ex hmm->print(stdout);
// ex hmm->print(stderr);
void OneCell3TrParam2DegPolyHMM::print(FILE* stream) {
  // print standard info first
  HMM::print(stream);

  // emission likelihood function
  fprintf(stream, "lambda_ij = ploidy_ij * %f * mu_i/2 + hat_mu/%i\n", this->getLibScalingFactor(0), this->DEPTH_ERROR_SCALING);
  fprintf(stream, "X_ij ~ NegBinom(mean=lambda_ij, variance=%f + %f*lambda_ij + %f*lambda_ij^2)\n", this->getMeanVarianceIntercept(), this->getMeanVarianceSlope(), this->getMeanVariancePoly2());
}

void OneCell3TrParam2DegPolyHMM::setLibScalingFactor(int cellNum, double libScalingFactor) {
  gsl_vector_set(this->paramsToEst, this->LIB_SIZE_SCALING_FACTOR_START_IDX + cellNum, libScalingFactor);
}
double OneCell3TrParam2DegPolyHMM::getLibScalingFactor(int cellNum) const {
  return gsl_vector_get(this->paramsToEst, this->LIB_SIZE_SCALING_FACTOR_START_IDX + cellNum);
}

/*
 ********
 * functions that depend on numbering and ordering of transition params
 ********
 */
/*
 * transitionParams = [beta, lambda, t]
 */
double OneCell3TrParam2DegPolyHMM::setTransition(gsl_matrix* dest, gsl_vector* transitionParams) {
  double beta   = gsl_vector_get(transitionParams, 0);
  double lambda = gsl_vector_get(transitionParams, 1);
  double t      = gsl_vector_get(transitionParams, 2);

  // save into paramsToEst if saving into this->transition
  if(dest == this->transition) {
    int transitionParamsIdx = 0;
    for(int i = 0; i < this->NUM_TRANSITION_PARAMS_TO_EST; i++) {
      gsl_vector_set(this->paramsToEst, this->TRANSITION_PROB_START_IDX + i, gsl_vector_get(transitionParams, transitionParamsIdx));
      transitionParamsIdx++;
    }
    for(int i = 0; i < this->NUM_BRANCH_LENGTHS_TO_EST; i++) {
      gsl_vector_set(this->paramsToEst, this->BRANCH_LENGTH_START_IDX + i, gsl_vector_get(transitionParams, transitionParamsIdx));
      transitionParamsIdx++;
    }
  }

  return this->setTransition(dest, this->getAlpha(), beta, lambda, t);
}

double OneCell3TrParam2DegPolyHMM::setTransition(gsl_matrix* dest, double alpha, double beta, double lambda, double t) {
  // first set rate matrix
  this->setRateMatrixQ(alpha, beta, lambda);

  // then set time dependent matrix P for time t
  double status = this->setTimeDepMatrixP(t);
  if(status != GSL_SUCCESS) {
    return status;
  }

  // then calculate matrix M(t) = {m_ij(t)}, where m_ij(t) = P_(2,2),(i,j)(t) / sum_v=0^k P_(2,2),(i,v)(t)
  int ancDiploidStateIdx = getStateIdxFromPloidyPair(2, 2); // index of ancestral diploid state (2,2)
  double currVal = 0;
  double rowsum = 0;
  for(unsigned int i = 0; i < this->states->size(); i++) {
    for(unsigned int j = 0; j < this->states->size(); j++) {
      currVal = gsl_matrix_get(this->timeDepMatrixP, ancDiploidStateIdx, getStateIdxFromPloidyPair(i, j));
      gsl_matrix_set(dest, i, j, currVal);
    }

    // rescale by rowsum
    gsl_vector_view currRow = gsl_matrix_row(dest, i);


    rowsum = gsl_blas_dasum(&currRow.vector); // Double Absolute SUM
    gsl_vector_scale(&currRow.vector, 1.0 / rowsum);

  }
  return status;
}

/*
 * Because BFGS optimizes unrestrained values, we need to transform
 * probabilities and parameters at different steps.
 *
 *   BFGS will optimize params x,y \in (-inf, inf)
 *   Probs a,b must be \in [0,1]
 *   and for the transition matrix, the constraint is
 *   //0 <= 1-2a-kb-g <= 1 <==> 2a+kb+g <= 1
 *   0 <= 1-t(2a+kb+g) <= 1 <==> 0 <= 2a+kb+g <= 1 && 0 <= t <= 1
 *
 * The following transformation converts params x,y,z to probabilities a,b,g
 *   2a = e^x / (1 + e^x + e^y + e^z) ==> a = e^x / [2 (1 + e^x + e^y + e^z)]
 *   kb = e^y / (1 + e^x + e^y + e^z) ==> b = e^y / [k (1 + e^x + e^y + e^z)]
 *   g = e^z / (1 + e^x + e^y + e^z) ==> g = e^z / [(1 + e^x + e^y + e^z)]
 * 
 * The following transformation converts probs a,b,g to params x,y,z
 *   Let c = 1+ e^x + e^y + e^z
 *   2a + kb + g + (1-2a-kb-g) = 1 ==> e^x / c + e^y / c + e^z / c + 1 / c = 1
 *     ==> (1-2a-kb-g) = 1 / c
 *     ==> c = 1 / (1-2a-kb-g)
 *   a = e^x / (2c) ==> a * 2c = e^x
 *   x = ln(a * 2c)
 *   y = ln(b * kc)
 *   z = ln(g * c)
 *
 * The following transformation converts branch length t (copied/reduced from TwoCell3TrParam2DegPolyHMM)
 *   d*t = e^(T) / (1 + exp(T))
 *   T = ln(-(d*t) / (d*t - 1))
 *
 * The library scaling factor (s) must be strictly positive. BFGS will optimize w \in (-inf, inf)
 *   s = exp(w) > 0
 *   w = ln(s)
 *
 * The following two functions convert the elements of src into the elements of dest
 * Also, copies over any other elements into dest (ie libSizeScalingFactor)
 */
void OneCell3TrParam2DegPolyHMM::convertProbToParam(gsl_vector* dest, const gsl_vector* src) const {
  // Wed 09 Sep 2020 04:09:01 PM PDT all vars must be > 0, but that's it. just exp/log everything for continuous time model
  double beta = gsl_vector_get(src, this->TRANSITION_PROB_START_IDX + 0);
  double lambda = gsl_vector_get(src, this->TRANSITION_PROB_START_IDX + 1);
  double t = gsl_vector_get(src, this->BRANCH_LENGTH_START_IDX);

  double s = gsl_vector_get(src, this->LIB_SIZE_SCALING_FACTOR_START_IDX);


  gsl_vector_set(dest, this->TRANSITION_PROB_START_IDX + 0, log(beta)); // set y
  gsl_vector_set(dest, this->TRANSITION_PROB_START_IDX + 1, log(lambda)); // set z
  gsl_vector_set(dest, this->BRANCH_LENGTH_START_IDX, log(t)); // set T
  gsl_vector_set(dest, this->LIB_SIZE_SCALING_FACTOR_START_IDX, log(s));
}

/*
 * reverse of convertProbToParam (see above)
 */
void OneCell3TrParam2DegPolyHMM::convertParamToProb(gsl_vector* dest, const gsl_vector* src) const {
  double y = gsl_vector_get(src, this->TRANSITION_PROB_START_IDX + 0);
  double z = gsl_vector_get(src, this->TRANSITION_PROB_START_IDX + 1);
  double T = gsl_vector_get(src, this->BRANCH_LENGTH_START_IDX);
  double w = gsl_vector_get(src, this->LIB_SIZE_SCALING_FACTOR_START_IDX);

  gsl_vector_set(dest, this->TRANSITION_PROB_START_IDX + 0, exp(y)); // set beta
  gsl_vector_set(dest, this->TRANSITION_PROB_START_IDX + 1, exp(z)); // set lambda
  gsl_vector_set(dest, this->BRANCH_LENGTH_START_IDX, exp(T)); // set t
  gsl_vector_set(dest, this->LIB_SIZE_SCALING_FACTOR_START_IDX, exp(w));
}

/*
 ********
 * functions
 ********
 */
double OneCell3TrParam2DegPolyHMM::getEmissionProb(double tumorDepth, double diploidDepth, int ploidy, int cellIdx) {
  return exp(this->getLogEmissionProb(tumorDepth, diploidDepth, ploidy, cellIdx));
}
double OneCell3TrParam2DegPolyHMM::getTotalEmissionProb(int stateIdx, std::vector<std::vector<double>*>* currChrDepthsVec, int chrIdx, int depthIdx) {
  return exp(this->getTotalLogEmissionProb(stateIdx, currChrDepthsVec, chrIdx, depthIdx)); // marginally faster
}

double OneCell3TrParam2DegPolyHMM::getLogEmissionProb(double tumorDepth, double diploidDepth, int ploidy, int cellIdx) {
  // if diploidDepth == 0, assume this is a hard to map region (ie any tumor reads mapping here are spurious)
  // P(any number reads | no data) = 1
  if(diploidDepth < 1e-4) {
    return 0; // return 0 since log(1) = 0
  }

  double currLibSizeScalingFactor = this->getLibScalingFactor(cellIdx);
  double lambda_i = (ploidy * diploidDepth/2.0) * currLibSizeScalingFactor + 272.5568 / this->DEPTH_ERROR_SCALING; // constant error

  // tumorDepth ~ NegBinom(lambda_i, var=intercept + slope * lambda_i + poly * lambda_i + poly2 * lambda_i * lambda_i
  double var = this->getMeanVarianceIntercept() + this->getMeanVarianceSlope() * lambda_i + this->getMeanVariancePoly2() * lambda_i * lambda_i;

  // in the boost library, p = P(success), r = # successes
  double p = lambda_i / var;
  double r = lambda_i * lambda_i / (var - lambda_i);

  // if diploidDepth is so low that p or r goes negative, assume this is a hard to map region
  if((p < 0 || r < 0) && diploidDepth < 1) {
    return 1;
  }

  if(r < 1) {
    r = 1;
    var = lambda_i * lambda_i + lambda_i;
    p = lambda_i / var;
  }

  double likelihood = -1;
  int k = (int) tumorDepth;
  likelihood = (lgamma(r+k) - lgamma(r) - this->getLogFacK(k) + r * log(p) + k * log(1-p));
  return likelihood;
}

// copied from TwoCell3TrParam2DegPolyHMM
double OneCell3TrParam2DegPolyHMM::getTotalLogEmissionProb(int stateIdx, std::vector<std::vector<double>*>* currChrDepthsVec, int chrIdx, int depthIdx) {
  double currDiploidDepth = (*(*currChrDepthsVec)[this->NUM_CELLS])[depthIdx];
  int currPloidy = getCellPloidyFromStateIdx(0, stateIdx);
  double currTumorDepth = (*(*currChrDepthsVec)[0])[depthIdx];
  return this->getLogEmissionProb(currTumorDepth, currDiploidDepth, currPloidy, 0);
}

/*
 * helper function to return log(k!) part of the negative binomial coefficient
 * given k (some read depth).
 * to be consistent with prior boost implementation:
 * f(k; r, p) = gamma(r+k)/(k! gamma(r)) * p^r * (1-p)^k
 * https://www.boost.org/doc/libs/1_70_0/libs/math/doc/html/math_toolkit/dist_ref/dists/negative_binomial_dist.html
 * The underlying member variable vector (this->logFacKVec)
 * has structure k:[log(k!) == lgamma(k+1)]
 */
double OneCell3TrParam2DegPolyHMM::getLogFacK(int k) {
  // if haven't called getLogFacK yet, alloc the lookup vector
  if(this->logFacKVec == nullptr) {
    this->setLogFacK();
  }
  return (*this->logFacKVec)[k];
}
/*
 * helper method to set this->logFacKVec
 */
void OneCell3TrParam2DegPolyHMM::setLogFacK() {
  // k is 0..maxReadDepth
  int maxDepth = this->maxObservedDepth;
  if(this->logFacKVec == nullptr) {
    this->logFacKVec = new std::vector<double>();
  }
  for(int i = this->logFacKVec->size(); i < maxDepth+1; i++) {
    this->logFacKVec->push_back(lgamma(i+1));
  }
}

/*
 * calls bfgs, returns a bestGuessHMM
 */
OneCell3TrParam2DegPolyHMM* OneCell3TrParam2DegPolyHMM::bfgs(gsl_vector* initGuess, int maxIters, bool verbose, bool debug) {
  OneCell3TrParam2DegPolyHMM* bestGuessHMM = this;
  gsl_vector* initGuessAsParams = gsl_vector_alloc(initGuess->size);
  this->convertProbToParam(initGuessAsParams, initGuess);
  Optimizable::bfgs(initGuessAsParams, bestGuessHMM, maxIters, verbose, debug);
  gsl_vector_free(initGuessAsParams);
  return bestGuessHMM;
}
  
/*
 * method to simulate a sequence of states, then a sequence of read depths
 * given those states.
 * Tue 11 Jun 2019 04:03:05 PM PDT seed is only used once (first time simulate is called)
 */
void OneCell3TrParam2DegPolyHMM::simulate() {
  this->simulate(43);
}
void OneCell3TrParam2DegPolyHMM::simulate(int seed) {
  this->simulate(seed, true); // by default, simulate diploid depths
}
void OneCell3TrParam2DegPolyHMM::simulate(int seed, bool simDiploid, double diploid_lambda_i, int numDiploidCells) {
  if(this->generator == nullptr) {
    this->seed = seed;
    this->generator = new base_generator_type(seed);
  }
  boost::uniform_real<> uni_dist(0,1);
  boost::variate_generator<base_generator_type&, boost::uniform_real<>> uni(*this->generator, uni_dist);

  // simulate states
  std::vector<std::string>* diploidStates = nullptr;
  std::vector<std::string>* tumorStates = nullptr;

  // clear out DepthPair variables. This is important if multiple simulations are run
  for(unsigned int i = 0; i < this->depths->size(); i++) {
    if(simDiploid) {
      (*this->depths)[i]->diploidLibrarySize = 0;
      (*this->depths)[i]->maxDiploidDepth = 0;
    }
    (*this->depths)[i]->tumorLibrarySize = 0;
    (*this->depths)[i]->maxTumorDepth = 0;
  }

  std::string currChr;
  int currNumWindows = -1;
  std::vector<std::string>* chrVec = this->getChrVec();
  int fromStateIdx = -1;
  int toStateIdx = -1;
  // for each chr
  for(unsigned int i = 0; i < chrVec->size(); i++) {
    // init this chr's vectors
    diploidStates = new std::vector<std::string>();
    tumorStates = new std::vector<std::string>();

    currChr = (*chrVec)[i];
    currNumWindows = (*(*this->depths)[0]->regions)[currChr]->size();

    // diploid always in state diploid
    diploidStates->push_back("2");

    // start tumor in steady state
    double stateProb = uni();
    fromStateIdx = getRandStateIdx(stateProb, this->initProb);
    tumorStates->push_back((*this->states)[fromStateIdx]);

    // for each window in this chr
    for(int winIdx = 1; winIdx < currNumWindows; winIdx++) {
      // diploid always in state diploid
      diploidStates->push_back("2");

      // tumor state changes according to transition matrix
      stateProb = uni();
      toStateIdx = getRandStateIdx(stateProb, fromStateIdx);
      tumorStates->push_back((*this->states)[toStateIdx]);
      fromStateIdx = toStateIdx;
    }

    // save
    (*(*this->depths)[0]->chrToDiploidSimStateMap)[currChr] = diploidStates;
    (*(*this->depths)[0]->chrToTumorSimStateMap)[currChr] = tumorStates;
  }

  double diploid_var = this->getMeanVarianceIntercept() + this->getMeanVarianceSlope() * diploid_lambda_i + this->getMeanVariancePoly2() * diploid_lambda_i * diploid_lambda_i; // Mon 20 Apr 2020 08:10:12 PM PDT now all diploid variance is calculated using the polynomial mean/var relationship
  double diploid_p = diploid_lambda_i / diploid_var;
  double diploid_r = diploid_lambda_i * diploid_lambda_i / (diploid_var - diploid_lambda_i);
  boost::random::negative_binomial_distribution<> diploid_negBinom_dist(diploid_r, diploid_p);
  boost::variate_generator<base_generator_type&, boost::random::negative_binomial_distribution<>> diploid_negBinom(*this->generator, diploid_negBinom_dist);

  // given states, simulate coverage according to neg binom model
  double tumor_lambda_i = -1;
  double tumor_var = -1;
  double tumor_p = -1;
  double tumor_r = -1;
  int currTumorPloidy = -1;
  double maxSim = -1;
  double currDiploidSim = -1;
  double currTumorSim = -1;
  double currLibSizeScalingFactor = -1;

  if(simDiploid) {
    // for each chr
    for(unsigned int j = 0; j < chrVec->size(); j++) {
      currChr = (*chrVec)[j];
      diploidStates = (*(*this->depths)[0]->chrToDiploidSimStateMap)[currChr];
      for(unsigned int i = 0; i < diploidStates->size(); i++) {
        double dipSimVec[numDiploidCells];
        for(int dipIdx = 0; dipIdx < numDiploidCells; dipIdx++) {
          dipSimVec[dipIdx] = diploid_negBinom();
          std::cerr << dipSimVec[dipIdx] << std::endl;
        }
        currDiploidSim = gsl_stats_mean(dipSimVec, 1, numDiploidCells);
        (*(*(*this->depths)[0]->chrToDiploidDepthMap)[currChr])[i] = currDiploidSim;
        if(currDiploidSim > maxSim) {
          maxSim = currDiploidSim;
          (*this->depths)[0]->maxDiploidDepth = currDiploidSim;
        }
        (*this->depths)[0]->diploidLibrarySize += currDiploidSim;
        (*this->depths)[1]->diploidLibrarySize += currDiploidSim;
        (*(*(*this->depths)[0]->chrToDiploidVarMap)[currChr])[i] = gsl_stats_variance(dipSimVec, 1, numDiploidCells); // Fri 17 Apr 2020 10:42:38 PM PDT debugging simulating multiple diploid cells, then averaging
      }
    }
  }

  // sim tumor
  // for each chr
  for(unsigned int j = 0; j < chrVec->size(); j++) {
    currChr = (*chrVec)[j];
    tumorStates = (*(*this->depths)[0]->chrToTumorSimStateMap)[currChr];

    // for each simulated state
    for(unsigned int i = 0; i < tumorStates->size(); i++) {
      diploid_lambda_i = (*(*(*this->depths)[0]->chrToDiploidDepthMap)[currChr])[i];

      // parse tumor states
      currTumorPloidy = atoi((*tumorStates)[i].c_str());
      currLibSizeScalingFactor = this->getLibScalingFactor(0);
      tumor_lambda_i = currLibSizeScalingFactor * (currTumorPloidy * diploid_lambda_i / 2 + diploid_lambda_i / this->DEPTH_ERROR_SCALING); // add error term withErr
      tumor_var = this->getMeanVarianceIntercept() + this->getMeanVarianceSlope() * tumor_lambda_i + this->getMeanVariancePoly2() * tumor_lambda_i * tumor_lambda_i;
      tumor_p = tumor_lambda_i / tumor_var;
      tumor_r = tumor_lambda_i * tumor_lambda_i / (tumor_var - tumor_lambda_i);
      if(tumor_r < 1) {
        // fix r=1, adj var, recalc p
        // r = lambda^2 / (var - lambda); v = lambda^2 / r + lambda
        tumor_r = 1;
        tumor_var = tumor_lambda_i * tumor_lambda_i + tumor_lambda_i;
        tumor_p = tumor_lambda_i / tumor_var;
      }

      boost::random::negative_binomial_distribution<> tumor_negBinom_dist(tumor_r,tumor_p);
      boost::variate_generator<base_generator_type&, boost::random::negative_binomial_distribution<>> tumor_negBinom(*this->generator, tumor_negBinom_dist);
      currTumorSim = tumor_negBinom();

      // save simulated values, as well as max for alphabet
      (*(*(*this->depths)[0]->chrToTumorDepthMap)[currChr])[i] = currTumorSim;
      if(currTumorSim > maxSim) {
        maxSim = currTumorSim;
        (*this->depths)[0]->maxTumorDepth = currTumorSim;
      }
      (*this->depths)[0]->tumorLibrarySize += currTumorSim;
    }
  }

  // save alphabet
  this->maxObservedDepth = (int)(maxSim) + 1;

  // save this->logFacKVec
  this->setLogFacK();

  // set library sizes according to final tumor size
  this->estLibScalingFactorsPosterior();
}

void OneCell3TrParam2DegPolyHMM::setUpBaumWelchLeastSquares() {
  this->baumWelchTransitionMat = gsl_matrix_alloc(this->transition->size1, this->transition->size2);
}

/*
 * param probs should be in prob space (ie already converted away from BFGS space). Order should be
 * [b]
 * [g]
 * [t]
 */
double OneCell3TrParam2DegPolyHMM::baumWelchLeastSquares_f(gsl_vector* probs) {
  double sumSqResid = 0;
  double resid = 0;

  // recalc what the transition matrix would be under these probs
  double status = this->setTransition(this->baumWelchTransitionMat, probs);
  if(gsl_isnan(status)) {
    return status;
  }

  // compare each entry in the hypothetical proposed transition matrix to the baum welch calculated transition matrix (ie sum up the squared difference)
  for(unsigned int row = 0; row < this->transition->size1; row++) {
    for(unsigned int col = 0; col < this->transition->size2; col++) {
      resid = (gsl_matrix_get(this->baumWelchTransitionMat, row, col) - gsl_matrix_get(this->transition, row, col));
      sumSqResid += resid * resid;
    }
  }
  return sumSqResid;
}

void OneCell3TrParam2DegPolyHMM::setInitGuessNthTime(gsl_vector* initGuess, int iter, int numTotalRuns) const {
  // TODO
}

/*
 * method to check validity of parameters proposed by BFGS. Assumes probs
 * is in probability space.
 * returns 0 if probs is valid, GSL_NAN otherwise
 */
double OneCell3TrParam2DegPolyHMM::checkOptimProbValidity(gsl_vector* probs) const {
  return 0; // Thu 27 Aug 2020 08:27:51 PM PDT debugging no validity check
  // shortcut for bad library sizes (too large or too small)
  double currLibSizeScalingFactor = -1;
  for(int i = 0; i < this->NUM_LIBS_TO_EST; i++) {
    currLibSizeScalingFactor = gsl_vector_get(probs, this->LIB_SIZE_SCALING_FACTOR_START_IDX + i);
    if(currLibSizeScalingFactor < 1e-2 || currLibSizeScalingFactor > 1e2) {
      return GSL_NAN;
    }
  }

  // shortcut for any probabilities becoming too small
  double probMin = gsl_vector_min(probs);
  if(probMin < 1e-5 || gsl_isnan(probMin)) {
    return GSL_NAN;
  }
  return 0;
}

