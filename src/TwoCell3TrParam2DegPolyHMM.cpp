#include "TwoCell3TrParam2DegPolyHMM.hpp"

/*
 ********
 * constructors and destructor
 ********
 */
TwoCell3TrParam2DegPolyHMM::TwoCell3TrParam2DegPolyHMM(std::vector<DepthPair*>* depths, int maxPloidy, bool preallocIntermediates) : TwoCell3TrParam2DegPolyHMM(depths, nullptr, maxPloidy, preallocIntermediates) {
}
TwoCell3TrParam2DegPolyHMM::TwoCell3TrParam2DegPolyHMM(std::vector<DepthPair*>* depths, gsl_vector* fixedParams, int maxPloidy, bool preallocIntermediates) : TwoCell3TrParam2DegPolyHMM(depths, fixedParams, maxPloidy, 2, 1, 0, 3, preallocIntermediates) { // 2 transition params to est (beta and gamma), 1 fixed param (alpha), 0 fixed libs (ie estimate all of them), 3 branches
  // delegate everything into the other ctor. Separating them is purely to be able to set NUM_TRANSITION_PARAMS_TO_EST to a custom value elsewhere
}
TwoCell3TrParam2DegPolyHMM::TwoCell3TrParam2DegPolyHMM(std::vector<DepthPair*>* depths, gsl_vector* fixedParams, int maxPloidy, int numTrParamsToEst, int numFixedTrParams, int numFixedLibs, int numBranches, bool preallocIntermediates) : HMM(depths, fixedParams, numTrParamsToEst, 2, numBranches, maxPloidy, numFixedTrParams, numFixedLibs, preallocIntermediates) { // 2 cells
  this->logFacKVec = nullptr;
  this->maxNumBFGSStarts = 3;
  this->setAlpha(0.1); // arbitrary starting value

  this->timeDepMatrixA  = gsl_matrix_alloc(this->timeDepMatrixP->size1, this->timeDepMatrixP->size2);
  this->timeDepMatrixP2  = gsl_matrix_alloc(this->timeDepMatrixP->size1, this->timeDepMatrixP->size2);
  this->timeDepMatrixP3  = gsl_matrix_alloc(this->timeDepMatrixP->size1, this->timeDepMatrixP->size2);
  gsl_matrix_set_zero(this->timeDepMatrixA);
  gsl_matrix_set_zero(this->timeDepMatrixP2);
  gsl_matrix_set_zero(this->timeDepMatrixP3);
}

TwoCell3TrParam2DegPolyHMM::~TwoCell3TrParam2DegPolyHMM() {
  delete this->logFacKVec;
  gsl_matrix_free(this->baumWelchTransitionMat);
  gsl_matrix_free(this->timeDepMatrixA);
  gsl_matrix_free(this->timeDepMatrixP2);
  gsl_matrix_free(this->timeDepMatrixP3);
}

// ex hmm->print(stdout);
// ex hmm->print(stderr);
void TwoCell3TrParam2DegPolyHMM::print(FILE* stream) {
  // print standard info first
  HMM::print(stream);

  // emission likelihood function
  fprintf(stream, "cell pair 0:\nlambda_ij = ploidy_ij * %f * mu_i/2 + hat_mu/%i\n", this->getLibScalingFactor(0), this->DEPTH_ERROR_SCALING);
  fprintf(stream, "X_ij ~ NegBinom(mean=lambda_ij, variance=%f + %f*lambda_ij + %f*lambda_ij^2)\n", this->getMeanVarianceIntercept(), this->getMeanVarianceSlope(), this->getMeanVariancePoly2());

  fprintf(stream, "\ncell pair 1:\nlambda_ij = ploidy_ij * %f * mu_i/2 + hat_mu/%i\n", this->getLibScalingFactor(1), this->DEPTH_ERROR_SCALING);
  fprintf(stream, "X_ij ~ NegBinom(mean=lambda_ij, variance=%f + %f*lambda_ij + %f*lambda_ij^2)\n", this->getMeanVarianceIntercept(), this->getMeanVarianceSlope(), this->getMeanVariancePoly2());
}

void TwoCell3TrParam2DegPolyHMM::setLibScalingFactor(int cellNum, double libScalingFactor) {
  gsl_vector_set(this->paramsToEst, this->LIB_SIZE_SCALING_FACTOR_START_IDX + cellNum, libScalingFactor);
}
double TwoCell3TrParam2DegPolyHMM::getLibScalingFactor(int cellNum) const {
  return gsl_vector_get(this->paramsToEst, this->LIB_SIZE_SCALING_FACTOR_START_IDX + cellNum);
}

/*
 ********
 * functions that depend on numbering and ordering of transition params
 ********
 */
double TwoCell3TrParam2DegPolyHMM::setTransition(gsl_matrix* dest, gsl_vector* transitionParams) {
  double beta  = gsl_vector_get(transitionParams, 0);
  double lambda = gsl_vector_get(transitionParams, 1);
  double t1 = gsl_vector_get(transitionParams, 2);
  double t2 = gsl_vector_get(transitionParams, 3);
  double t3 = gsl_vector_get(transitionParams, 4);

  if(dest == this->transition) {
    // save into paramsToEst
    int transitionParamsIdx = 0;
    for(int i = 0; i < this->NUM_TRANSITION_PARAMS_TO_EST; i++) {
      if((unsigned int) (this->TRANSITION_PROB_START_IDX + i) < this->paramsToEst->size) {
        gsl_vector_set(this->paramsToEst, this->TRANSITION_PROB_START_IDX + i, gsl_vector_get(transitionParams, transitionParamsIdx));
        transitionParamsIdx++;
      }
    }
    for(int i = 0; i < this->NUM_BRANCH_LENGTHS_TO_EST; i++) {
      gsl_vector_set(this->paramsToEst, this->BRANCH_LENGTH_START_IDX + i, gsl_vector_get(transitionParams, transitionParamsIdx));
      transitionParamsIdx++;
    }
  }

  return this->setTransition(dest, this->getAlpha(), beta, lambda, t1, t2, t3);
}

/*
 * helper method to set the transition matrix given these specific parameters:
 *   alpha: P(one CNA)
 *   beta: P(any CNA)
 *   gamma: P(return to diploid)
 *   t1: shared rate (fixed at 1)
 *   t2: cell A specific rate
 *   t3: cell B specific rate
*/
double TwoCell3TrParam2DegPolyHMM::setTransition(gsl_matrix* dest, double alpha, double beta, double lambda, double t1, double t2, double t3) {
  // first set rate matrix, if it's possible for beta and lambda to change
  this->setRateMatrixQ(alpha, beta, lambda);

  // then set time dependent matrix A for times t1,t2,t3
  double status = this->setTimeDepMatrixA(t1, t2, t3);

  // then calculate matrix M(t) = {m_(i,j),(i',j')(\bar t)}, where m_(i,j),(i',j')(\bar t) = A_(i,j),(i',j')(t) / sum_(i'',j'') A_(i,j),(i'',j'')
  // this is just normalizing matrix A by each rowsum
  double rowsum = 0;
  gsl_matrix_memcpy(dest, this->timeDepMatrixA);
  for(unsigned int row = 0; row < dest->size1; row++) {
    // rescale by rowsum
    gsl_vector_view currRow = gsl_matrix_row(dest, row);

    rowsum = gsl_blas_dasum(&currRow.vector); // Double Absolute SUM
    gsl_vector_scale(&currRow.vector, 1.0 / rowsum);
  }

  return status;
}

/*
 * helper function to set time dependent matrix A, which is the multiplication of time dependent matrix P.
 * that is, P depends on one time, and A depends on 3 times.
 * A_(i,j),(i',j')(t) = sum_(L,V) P_(2,2),(L,V)(t1) * P_(L,V),(i,i')(t2) * P_(L,V),(j,j')(t3)
 * a pair (i,j) or (i',j') or (L,V) is the equivalent of one state, where i is the ploidy for cell0, j is the ploidy for cell1
 */
double TwoCell3TrParam2DegPolyHMM::setTimeDepMatrixA(double t1, double t2, double t3) {
  // set timeDepMatrixP, timeDepMatrixP2, timeDepMatrixP3 TODO handle the status of each call
  double status = this->setTimeDepMatrixP(this->timeDepMatrixP, t1);
  status += this->setTimeDepMatrixP(this->timeDepMatrixP2, t2);
  status += this->setTimeDepMatrixP(this->timeDepMatrixP3, t3);

  // iterate over A
  int ancDiploidStateIdx = getStateIdxFromPloidyPair(2, 2); // index of ancestral diploid state (2,2)
  int i = 0;
  int j = 0;
  int iPrime = 0;
  int jPrime = 0;
  int lvStateIdx = 0;
  double currStateSum = 0;
  double currLVProd = 0;
  for(unsigned int row = 0; row < this->timeDepMatrixA->size1; row++) {
    i = getIndvPloidyFromStateIdx(0, row);
    j = getIndvPloidyFromStateIdx(1, row);
    for(unsigned int col = 0; col < this->timeDepMatrixA->size2; col++) {
      iPrime = getIndvPloidyFromStateIdx(0, col);
      jPrime = getIndvPloidyFromStateIdx(1, col);
      currStateSum = 0;
      for(int L = 0; L <= this->MAX_PLOIDY; L++) {
        for(int V = 0; V <= this->MAX_PLOIDY; V++) {
          lvStateIdx = getStateIdxFromPloidyPair(L, V);
          currLVProd = gsl_matrix_get(this->timeDepMatrixP, ancDiploidStateIdx, lvStateIdx);
          currLVProd *= gsl_matrix_get(this->timeDepMatrixP2, lvStateIdx, getStateIdxFromPloidyPair(i, iPrime));
          currLVProd *= gsl_matrix_get(this->timeDepMatrixP3, lvStateIdx, getStateIdxFromPloidyPair(j, jPrime));
          currStateSum += currLVProd;
        }
      }
      gsl_matrix_set(this->timeDepMatrixA, row, col, currStateSum);
    }
  }

  return status;
}

/*
 * Because BFGS optimizes unrestrained values, we need to transform
 * probabilities and parameters at different steps.
 *
 *   BFGS will optimize params dT1,dT2,dT3,v,w,x,y,z \in (-inf, inf)
 *   Probs a,b,g must be \in [0,1]
 *   branch lengths and library scaling factors must be > 0
 *
 * Wed 18 Dec 2019 04:19:36 PM PST
 * Now we're fixing alpha at some arbitrary value (ex. 0.1) and letting all branch lengths vary instead. Let c = -2a for wolfram alpha
 *   kb = exp(y) / (1 -2a + exp(y) + exp(z))
 *   g  = exp(z) / (1 -2a + exp(y) + exp(z))
 * 
 *   y = ln((-b(1-2a)k) / (bk + g - 1))
 *   z = ln((-g(1-2a))  / (bk + g - 1))
 * 
 * checked with https://www.wolframalpha.com/input/?i=Solve%5B+k*b+%3D%3D+Exp%5By%5D%2F%281%2B+c+%2B+Exp%5By%5D+%2B+Exp%5Bz%5D%29%2C+++g+%3D%3D+Exp%5Bz%5D%2F%281+%2Bc%2B+Exp%5By%5D+%2B+Exp%5Bz%5D%29%7D%2C+%7B+y%2C+z%7D%5D
 *
 * To allow t1 to vary as well:
 *   dt1 = e^(T1) / (1 + exp(T1) + exp(T2) + exp(T3))
 *   dt2 = e^(T2) / (1 + exp(T1) + exp(T2) + exp(T3))
 *   dt3 = e^(T3) / (1 + exp(T1) + exp(T2) + exp(T3))
 *   T1 = ln(-(dt1) / (d(t1 + t2 + t3) - 1))
 *   T2 = ln(-(dt2) / (d(t1 + t2 + t3) - 1))
 *   T3 = ln(-(dt3) / (d(t1 + t2 + t3) - 1))
 *
 * checked with https://www.wolframalpha.com/input/?i=Solve%5B%7B+d*s+%3D%3D+Exp%5BA%5D%2F%281+%2B+Exp%5BA%5D+%2B+Exp%5BB%5D+%2B+Exp%5BC%5D%29%2C++++d*t+%3D%3D+Exp%5BB%5D%2F%281+%2B+Exp%5BA%5D+%2B+Exp%5BB%5D+%2B+Exp%5BC%5D%29%2C+d*u+%3D%3D+Exp%5BC%5D%2F%281+%2B+Exp%5BA%5D+%2B+Exp%5BB%5D+%2B+Exp%5BC%5D%29%7D%2C+%7BA%2C+B%2C+C%7D%5D
 *
 * The library scaling factors (r,s) must be strictly positive. BFGS will optimize v,w \in (-inf, inf)
 *   r = exp(v) > 0
 *   s = exp(w) > 0
 *   v = ln(r)
 *   w = ln(s)
 *
 * The following two functions convert the elements of src into the elements of dest
 */
void TwoCell3TrParam2DegPolyHMM::convertProbToParam(gsl_vector* dest, const gsl_vector* src) const {
  double beta = gsl_vector_get(src, this->TRANSITION_PROB_START_IDX + 0);
  double lambda = gsl_vector_get(src, this->TRANSITION_PROB_START_IDX + 1);
  double t1 = gsl_vector_get(src, this->TRANSITION_PROB_START_IDX + 2);
  double t2 = gsl_vector_get(src, this->TRANSITION_PROB_START_IDX + 3);
  double t3 = gsl_vector_get(src, this->TRANSITION_PROB_START_IDX + 4);
  double r = gsl_vector_get(src, this->LIB_SIZE_SCALING_FACTOR_START_IDX);
  double s = gsl_vector_get(src, this->LIB_SIZE_SCALING_FACTOR_START_IDX + 1);

  gsl_vector_set(dest, this->TRANSITION_PROB_START_IDX + 0, log(beta)); // set y
  gsl_vector_set(dest, this->TRANSITION_PROB_START_IDX + 1, log(lambda)); // set z
  gsl_vector_set(dest, this->TRANSITION_PROB_START_IDX + 2, log(t1)); // set T1
  gsl_vector_set(dest, this->TRANSITION_PROB_START_IDX + 3, log(t2)); // set T2
  gsl_vector_set(dest, this->TRANSITION_PROB_START_IDX + 4, log(t3)); // set T3
  gsl_vector_set(dest, this->LIB_SIZE_SCALING_FACTOR_START_IDX, log(r));
  gsl_vector_set(dest, this->LIB_SIZE_SCALING_FACTOR_START_IDX + 1, log(s));

}

/*
 * reverse of convertProbToParam (see above)
 */
void TwoCell3TrParam2DegPolyHMM::convertParamToProb(gsl_vector* dest, const gsl_vector* src) const {
  double y = gsl_vector_get(src, this->TRANSITION_PROB_START_IDX + 0);
  double z = gsl_vector_get(src, this->TRANSITION_PROB_START_IDX + 1);
  double T1 = gsl_vector_get(src, this->TRANSITION_PROB_START_IDX + 2);
  double T2 = gsl_vector_get(src, this->TRANSITION_PROB_START_IDX + 3);
  double T3 = gsl_vector_get(src, this->TRANSITION_PROB_START_IDX + 4);
  double w = gsl_vector_get(src, this->LIB_SIZE_SCALING_FACTOR_START_IDX);
  double v = gsl_vector_get(src, this->LIB_SIZE_SCALING_FACTOR_START_IDX + 1);

  gsl_vector_set(dest, this->TRANSITION_PROB_START_IDX + 0, exp(y)); // beta
  gsl_vector_set(dest, this->TRANSITION_PROB_START_IDX + 1, exp(z)); // lambda
  gsl_vector_set(dest, this->TRANSITION_PROB_START_IDX + 2, exp(T1)); // set t1
  gsl_vector_set(dest, this->TRANSITION_PROB_START_IDX + 3, exp(T2)); // set t2
  gsl_vector_set(dest, this->TRANSITION_PROB_START_IDX + 4, exp(T3)); // set t3
  gsl_vector_set(dest, this->LIB_SIZE_SCALING_FACTOR_START_IDX, exp(w));
  gsl_vector_set(dest, this->LIB_SIZE_SCALING_FACTOR_START_IDX + 1, exp(v));
}


/*
 ********
 * functions
 ********
 */
/*
 * in this class, getEmissionProb is always calculated on the fly (with the exception
 * of the logFacK entries). So, chrIdx and depthIdx (which are used to look up stored
 * emission probs in FixLib subclasses) are ignored here)
 */
double TwoCell3TrParam2DegPolyHMM::getEmissionProb(double tumorDepth, double diploidDepth, int ploidy, int cellIdx) {
  return exp(this->getLogEmissionProb(tumorDepth, diploidDepth, ploidy, cellIdx));
}

/*
 * returns emission probability across all cells for a given state, chromosome, and depthIdx (position along chr).
 * emission prob is NOT in log space.
 * currChrDepthsVec is useful for caching and provides a speed up
 */
double TwoCell3TrParam2DegPolyHMM::getTotalEmissionProb(int stateIdx, std::vector<std::vector<double>*>* currChrDepthsVec, int chrIdx, int depthIdx) {
  double emissionProb = 0;
  int currPloidy = -1;
  double currTumorDepth = -1;
  double currDiploidDepth = (*(*currChrDepthsVec)[this->NUM_CELLS])[depthIdx];
  // for each cell, calc the emission prob and multiply by it
  for(int cellIdx = 0; cellIdx < this->NUM_CELLS; cellIdx++) {
    currPloidy = getCellPloidyFromStateIdx(cellIdx, stateIdx);
    currTumorDepth = (*(*currChrDepthsVec)[cellIdx])[depthIdx];
    emissionProb += this->getLogEmissionProb(currTumorDepth, currDiploidDepth, currPloidy, cellIdx);
  }
  return exp(emissionProb);
}

double TwoCell3TrParam2DegPolyHMM::getLogEmissionProb(double tumorDepth, double diploidDepth, int ploidy, int cellIdx) {
  // if diploidDepth == 0, assume this is a hard to map region (ie any tumor reads mapping here are spurious)
  // P(any number reads | no data) = 1
  if(diploidDepth < 1e-4) {
    return 0; // return 0 since log(1) = 0
  }

  // if scaling factor becomes super small or ridiculously large, return error code
  double currLibSizeScalingFactor = this->getLibScalingFactor(cellIdx);
  double lambda_i = (ploidy * diploidDepth/2.0) * currLibSizeScalingFactor + 272.5568 / this->DEPTH_ERROR_SCALING; // Thu 27 May 2021 10:50:36 AM PDT constant error

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

/*
 * returns emission probability across all cells for a given state, chromosome, and depthIdx (position along chr).
 * emission prob is in log space.
 * currChrDepthsVec is useful for caching and provides a speed up
 */
double TwoCell3TrParam2DegPolyHMM::getTotalLogEmissionProb(int stateIdx, std::vector<std::vector<double>*>* currChrDepthsVec, int chrIdx, int depthIdx) {
  double emissionProb = 0;
  int currPloidy = -1;
  double currTumorDepth = -1;
  double currDiploidDepth = (*(*currChrDepthsVec)[this->NUM_CELLS])[depthIdx];
  // for each cell, calc the emission prob and multiply by it
  for(int cellIdx = 0; cellIdx < this->NUM_CELLS; cellIdx++) {
    currPloidy = getCellPloidyFromStateIdx(cellIdx, stateIdx);
    currTumorDepth = (*(*currChrDepthsVec)[cellIdx])[depthIdx];
    emissionProb += this->getLogEmissionProb(currTumorDepth, currDiploidDepth, currPloidy, cellIdx);
  }
  return emissionProb;
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
double TwoCell3TrParam2DegPolyHMM::getLogFacK(int k) {
  // if haven't called getLogFacK yet, alloc the lookup vector
  if(this->logFacKVec == nullptr) {
    this->setLogFacK();
  }
  return (*this->logFacKVec)[k];
}
/*
 * helper method to set this->logFacKVec
 */
void TwoCell3TrParam2DegPolyHMM::setLogFacK() {
  // k is 0..maxReadDepth
  int maxDepth = this->maxObservedDepth;
  this->setLogFacK(maxDepth);
}
/*
 * helper method to set this->logFacKVec to a custom maxDepth, where maxDepth is included
 */
void TwoCell3TrParam2DegPolyHMM::setLogFacK(int maxDepth) {
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
TwoCell3TrParam2DegPolyHMM* TwoCell3TrParam2DegPolyHMM::bfgs(gsl_vector* initGuess, int maxIters, bool verbose, bool debug) {
  TwoCell3TrParam2DegPolyHMM* bestGuessHMM = this;//new TwoCell3TrParam2DegPolyHMM(*this);
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
void TwoCell3TrParam2DegPolyHMM::simulate() {
  this->simulate(43);
}
void TwoCell3TrParam2DegPolyHMM::simulate(int seed) {
  this->simulate(seed, true); // by default, simulate diploid depths
}

/*
 * function to simulate states and read depths. If(simDiploid) then diploid
 * depths are simulated. If not, then only tumor depths are simulated. This is
 * useful for multiple tumors simulated from the same parameters, since the average
 * diploid should not change between simulations. There no checks that diploid info has
 * been previously simulated
 *
 */
void TwoCell3TrParam2DegPolyHMM::simulate(int seed, bool simDiploid, double diploid_lambda_i, int numDiploidCells) {
  if(this->generator == nullptr) {
    this->seed = seed;
    this->generator = new base_generator_type(seed);
  }
  // set up random number generator
  // from https://www.boost.org/doc/libs/1_70_0/libs/random/example/random_demo.cpp
  boost::uniform_real<> uni_dist(0,1);
  boost::variate_generator<base_generator_type&, boost::uniform_real<>> uni(*this->generator, uni_dist);

  // simulate states
  std::vector<std::string>* diploidStates = nullptr;
  std::vector<std::string>* tumor0States = nullptr;
  std::vector<std::string>* tumor1States = nullptr;

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
    tumor0States = new std::vector<std::string>();
    tumor1States = new std::vector<std::string>();

    currChr = (*chrVec)[i];
    currNumWindows = (*(*this->depths)[0]->regions)[currChr]->size();

    // diploid always in state diploid
    diploidStates->push_back("2");

    // start tumor in steady state
    double stateProb = uni();
    //std::cout << stateProb << std::endl;
    fromStateIdx = getRandStateIdx(stateProb, this->initProb);

    // get tumor cell specific states
    tumor0States->push_back(std::to_string(getCellPloidyFromStateIdx(0, fromStateIdx)));
    tumor1States->push_back(std::to_string(getCellPloidyFromStateIdx(1, fromStateIdx)));

    // for each window in this chr
    for(int winIdx = 1; winIdx < currNumWindows; winIdx++) {
      // diploid always in state diploid
      diploidStates->push_back("2");

      // tumor state changes according to transition matrix
      stateProb = uni();
      //std::cout << stateProb << std::endl;
      toStateIdx = getRandStateIdx(stateProb, fromStateIdx);

      // if toStateIdx is different from where came from, then must have transitioned states
      fromStateIdx = toStateIdx;

      // get tumor cell specific states
      tumor0States->push_back(std::to_string(getCellPloidyFromStateIdx(0, fromStateIdx)));
      tumor1States->push_back(std::to_string(getCellPloidyFromStateIdx(1, fromStateIdx)));
    }

    // save
    (*(*this->depths)[0]->chrToDiploidSimStateMap)[currChr] = diploidStates;
    (*(*this->depths)[0]->chrToTumorSimStateMap)[currChr] = tumor0States;
    (*(*this->depths)[1]->chrToTumorSimStateMap)[currChr] = tumor1States;
  }

  // values from means of input/diploidBinsVarMeanNormLibSize_250kb.coordMeanVar
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
  int currTumor0Ploidy = -1;
  int currTumor1Ploidy = -1;
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
        //std::cout << currDiploidSim << std::endl;
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

  // sim tumor 0
  // for each chr
  for(unsigned int j = 0; j < chrVec->size(); j++) {
    currChr = (*chrVec)[j];
    tumor0States = (*(*this->depths)[0]->chrToTumorSimStateMap)[currChr];

    // for each simulated state
    for(unsigned int i = 0; i < tumor0States->size(); i++) {
      diploid_lambda_i = (*(*(*this->depths)[0]->chrToDiploidDepthMap)[currChr])[i];

      // parse tumor states
      currTumor0Ploidy = atoi((*tumor0States)[i].c_str());
      currLibSizeScalingFactor = this->getLibScalingFactor(0);
      tumor_lambda_i = currLibSizeScalingFactor * (currTumor0Ploidy * diploid_lambda_i / 2 + diploid_lambda_i / this->DEPTH_ERROR_SCALING); // add error term withErr
      tumor_var = this->getMeanVarianceIntercept() + this->getMeanVarianceSlope() * tumor_lambda_i + this->getMeanVariancePoly2() * tumor_lambda_i * tumor_lambda_i;
      tumor_p = tumor_lambda_i / tumor_var;
      tumor_r = tumor_lambda_i * tumor_lambda_i / (tumor_var - tumor_lambda_i);
      //std::cerr << currTumor0Ploidy << ", " << diploid_lambda_i << std::endl;
      if(tumor_r < 1) {
        // fix r=1, adj var, recalc p
        // r = lambda^2 / (var - lambda); v = lambda^2 / r + lambda
        tumor_r = 1;
        tumor_var = tumor_lambda_i * tumor_lambda_i + tumor_lambda_i;
        tumor_p = tumor_lambda_i / tumor_var;
      }

      boost::random::negative_binomial_distribution<> tumor0_negBinom_dist(tumor_r,tumor_p);
      boost::variate_generator<base_generator_type&, boost::random::negative_binomial_distribution<>> tumor0_negBinom(*this->generator, tumor0_negBinom_dist);
      currTumorSim = tumor0_negBinom();

      // save simulated values, as well as max for alphabet
      (*(*(*this->depths)[0]->chrToTumorDepthMap)[currChr])[i] = currTumorSim;
      if(currTumorSim > maxSim) {
        maxSim = currTumorSim;
        (*this->depths)[0]->maxTumorDepth = currTumorSim;
      }
      (*this->depths)[0]->tumorLibrarySize += currTumorSim;
    }
  }

  // simulate tumor 1
  for(unsigned int j = 0; j < chrVec->size(); j++) {
    currChr = (*chrVec)[j];
    tumor1States = (*(*this->depths)[1]->chrToTumorSimStateMap)[currChr];
    // for each simulated state
    for(unsigned int i = 0; i < tumor1States->size(); i++) {
      diploid_lambda_i = (*(*(*this->depths)[0]->chrToDiploidDepthMap)[currChr])[i];

      // parse tumor states
      currTumor1Ploidy = atoi((*tumor1States)[i].c_str());
      currLibSizeScalingFactor = this->getLibScalingFactor(1);
      tumor_lambda_i = currLibSizeScalingFactor * (currTumor1Ploidy * diploid_lambda_i / 2 + diploid_lambda_i / this->DEPTH_ERROR_SCALING); // add error term withErr
      tumor_var = this->getMeanVarianceIntercept() + this->getMeanVarianceSlope() * tumor_lambda_i + this->getMeanVariancePoly2() * tumor_lambda_i * tumor_lambda_i;
      tumor_p = tumor_lambda_i / tumor_var;
      tumor_r = tumor_lambda_i * tumor_lambda_i / (tumor_var - tumor_lambda_i);
      if(tumor_r < 1) {
        tumor_r = 1;
        tumor_var = tumor_lambda_i * tumor_lambda_i + tumor_lambda_i;
        tumor_p = tumor_lambda_i / tumor_var;
      }

      boost::random::negative_binomial_distribution<> tumor1_negBinom_dist(tumor_r,tumor_p);
      boost::variate_generator<base_generator_type&, boost::random::negative_binomial_distribution<>> tumor1_negBinom(*this->generator, tumor1_negBinom_dist);
      currTumorSim = tumor1_negBinom();

      // save simulated values, as well as max for alphabet
      (*(*(*this->depths)[1]->chrToTumorDepthMap)[currChr])[i] = currTumorSim;
      if(currTumorSim > maxSim) {
        maxSim = currTumorSim;
        (*this->depths)[1]->maxTumorDepth = currTumorSim;
      }
      (*this->depths)[1]->tumorLibrarySize += currTumorSim;
    }
  }

  // save alphabet
  this->maxObservedDepth = (int)(maxSim) + 1;

  // save this->logFacKVec
  this->setLogFacK();

  // set library sizes according to final tumor size
  this->estLibScalingFactorsPosterior();
}

/*
 * some param sets may not result in valid transition matrices (ie some entries will go
 * negative, depending on k). If this happens, all the likelihood methods will spit out nans, and there's
 * nothing done about it so far (Mon 16 Sep 2019 05:16:49 PM PDT)
 */
void TwoCell3TrParam2DegPolyHMM::setInitGuessNthTime(gsl_vector* initGuess, int iter, int numTotalRuns) const {
  switch(iter) {
    // beta/gamma all at 0.05, t1/t2/t3 at 0.1
    case 0:
      for(int cellIdx = 0; cellIdx < this->NUM_LIBS_TO_EST; cellIdx++) {
        gsl_vector_set(initGuess, this->LIB_SIZE_SCALING_FACTOR_START_IDX + cellIdx, 1);
      }
      gsl_vector_set(initGuess, this->TRANSITION_PROB_START_IDX + 0, 0.01); // beta
      gsl_vector_set(initGuess, this->TRANSITION_PROB_START_IDX + 1, 0.01); // gamma
      gsl_vector_set(initGuess, this->TRANSITION_PROB_START_IDX + 2, 0.2); // t1
      gsl_vector_set(initGuess, this->TRANSITION_PROB_START_IDX + 3, 0.2);
      gsl_vector_set(initGuess, this->TRANSITION_PROB_START_IDX + 4, 0.2);
      break;

    // same as case 0, but lib sizes start at 1
    case 1:
      for(int cellIdx = 0; cellIdx < this->NUM_LIBS_TO_EST; cellIdx++) {
        gsl_vector_set(initGuess, this->LIB_SIZE_SCALING_FACTOR_START_IDX + cellIdx, .75);
      }
      gsl_vector_set(initGuess, this->TRANSITION_PROB_START_IDX + 0, 0.05);
      gsl_vector_set(initGuess, this->TRANSITION_PROB_START_IDX + 1, 0.05);
      gsl_vector_set(initGuess, this->TRANSITION_PROB_START_IDX + 2, 0.1);
      gsl_vector_set(initGuess, this->TRANSITION_PROB_START_IDX + 3, 0.1);
      gsl_vector_set(initGuess, this->TRANSITION_PROB_START_IDX + 4, 0.1);
      break;

    // start lib sizes at right place, uniform start t 0.01
    case 2:
      for(int cellIdx = 0; cellIdx < this->NUM_LIBS_TO_EST; cellIdx++) {
        gsl_vector_set(initGuess, this->LIB_SIZE_SCALING_FACTOR_START_IDX + cellIdx, 1.25);
      }
      gsl_vector_set(initGuess, this->TRANSITION_PROB_START_IDX + 0, 0.02);
      gsl_vector_set(initGuess, this->TRANSITION_PROB_START_IDX + 1, 0.05);
      gsl_vector_set(initGuess, this->TRANSITION_PROB_START_IDX + 2, 0.02);
      gsl_vector_set(initGuess, this->TRANSITION_PROB_START_IDX + 3, 0.02);
      gsl_vector_set(initGuess, this->TRANSITION_PROB_START_IDX + 4, 0.02);
      break;
  }
}

/*
 * function to solve the overdetermined system of equations from the Baum Welch estimated transition matrix, a*
 * Assumes runBaumWelch() has already been called
 */
void TwoCell3TrParam2DegPolyHMM::setUpBaumWelchLeastSquares() {
  this->baumWelchTransitionMat = gsl_matrix_alloc(this->transition->size1, this->transition->size2);
}

/*
 * param probs should be in prob space (ie already converted away from BFGS space). Order should be
 * [b]
 * //[g]
 * [L]
 * [t1]
 * [t2]
 * [t3]
 */
double TwoCell3TrParam2DegPolyHMM::baumWelchLeastSquares_f(gsl_vector* probs) {
  // recalc what the transition matrix would be under these probs
  this->setTransition(this->baumWelchTransitionMat, probs);
  // compare each entry in the hypothetical proposed transition matrix to the baum welch calculated transition matrix (ie sum up the squared difference)
  double sumSqResid = 0;
  double resid = 0;
  for(unsigned int row = 0; row < this->transition->size1; row++) {
    for(unsigned int col = 0; col < this->transition->size2 - 1; col++) {
      double bwVal = gsl_matrix_get(this->baumWelchTransitionMat, row, col);
      double trVal = gsl_matrix_get(this->transition, row, col);
      resid = ((bwVal - trVal) * (bwVal - trVal));// / bwVal;
      sumSqResid += resid;
    }
  }
  return sumSqResid;
}

/*
 * method to check validity of parameters proposed by BFGS. Assumes probs
 * is in probability space.
 * returns 0 if probs is valid, GSL_NAN otherwise
 */
double TwoCell3TrParam2DegPolyHMM::checkOptimProbValidity(gsl_vector* probs) const {
  return 0; // Thu 27 Aug 2020 08:27:51 PM PDT debugging no validity check
  // shortcut for bad library sizes (too large or too small)
  double currLibSizeScalingFactor = -1;
  for(int i = 0; i < this->NUM_LIBS_TO_EST; i++) {
    currLibSizeScalingFactor = gsl_vector_get(probs, this->LIB_SIZE_SCALING_FACTOR_START_IDX + i);
    if(gsl_isinf(currLibSizeScalingFactor) || currLibSizeScalingFactor < 1e-2 || currLibSizeScalingFactor > 1e2) {
      std::cerr << "shortcut for bad lib sizes: ";// << std::endl;
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

