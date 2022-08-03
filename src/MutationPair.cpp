#include "MutationPair.hpp"

MutationPair::MutationPair(std::vector<MutationList*>* mutListVec, gsl_vector* initGuess, bool verbose, bool gradientDebug) : MutationListContainer(mutListVec, 3, 0, verbose, gradientDebug) { // 3 mutCounts, 0 overdispVars
  // get muts that have data in both (intersection) or either (union)
  MutationList* mutListA = (*this->mutListVec)[0];
  MutationList* mutListB = (*this->mutListVec)[1];
  std::vector<std::string> sortedA(mutListA->coordVec->begin(), mutListA->coordVec->end());
  std::vector<std::string> sortedB(mutListB->coordVec->begin(), mutListB->coordVec->end());
  std::sort(sortedA.begin(), sortedA.end());
  std::sort(sortedB.begin(), sortedB.end());
  //// intersection
  //this->sharedMutationVec = new std::vector<std::string>(sortedA.size());
  //std::vector<std::string>::iterator lastElement = std::set_intersection(sortedA.begin(), sortedA.end(), sortedB.begin(), sortedB.end(), this->sharedMutationVec->begin());
  // union
  this->sharedMutationVec = new std::vector<std::string>(sortedA.size() + sortedB.size());
  std::vector<std::string>::iterator lastElement = std::set_union(sortedA.begin(), sortedA.end(), sortedB.begin(), sortedB.end(), this->sharedMutationVec->begin());
  this->sharedMutationVec->resize(lastElement - this->sharedMutationVec->begin());

  // estimate x1/x2/x3
  this->internalInitGuess = gsl_vector_alloc(this->getNumParamsToEst());
  if(initGuess != nullptr) {
    gsl_vector_memcpy(this->internalInitGuess, initGuess);
    gsl_vector_memcpy(this->paramsToEst, initGuess);
  }
  else {
    gsl_vector_set_all(this->paramsToEst, 5);
    //gsl_vector_set_all(this->internalInitGuess, 10);
    gsl_vector_memcpy(this->internalInitGuess, this->paramsToEst);
  }
  //std::cout << "MutationPair ctor going to run bfgs" << std::endl;
  //this->print(stdout);
  //this->estMutCountsPerBranch(nullptr, this->internalInitGuess, 500, this->verbose); // mutation bfgs short circuit here
  //this->bfgs(this->internalInitGuess, 2, this->verbose, this->gradientDebug); // TODO change back to 500
  //this->bfgs(this->internalInitGuess, 500, this->verbose, this->gradientDebug);
  //gsl_vector_free(initGuess);
}

MutationPair::~MutationPair() {
  delete this->sharedMutationVec;
}

/*
 * calculates l(X1,X2,X3) = sum_sites [log(p(data_i | S_i=(1,1)) * X1/N + p(data_i | S_i=(1,0)) * X2/N + p(data_i | S_i=(0,1))*X3/N + p(data_i | S_i=(0,0))*(N-X1-X2-X3)/N)]
 */
double MutationPair::getLogLikelihood() {
  //return GSL_NAN;
  double totalLl = 0;
  double currLl = 0;
  double x1 = gsl_vector_get(this->paramsToEst, 0);
  double x2 = gsl_vector_get(this->paramsToEst, 1);
  double x3 = gsl_vector_get(this->paramsToEst, 2);
  double N = this->sharedMutationVec->size();
  double prob_s11 = 0; // P(D_i | S_i = (1,1))
  double prob_s10 = 0;
  double prob_s01 = 0;
  double prob_s00 = 0;
  std::string currSite;
  MutationList* mutListA = (*this->mutListVec)[0];
  MutationList* mutListB = (*this->mutListVec)[1];
  double prob_a1 = 0; // P(D_A | S_A = 1)
  double prob_a0 = 0; // P(D_A | S_A = 0)
  double prob_b1 = 0; // P(D_B | S_A = 1)
  double prob_b0 = 0; // P(D_B | S_A = 0)
  double totalProb = 0;

  // for each shared site
  for(unsigned int sharedMutIdx = 0; sharedMutIdx < this->sharedMutationVec->size(); sharedMutIdx++) {
    currSite = (*this->sharedMutationVec)[sharedMutIdx];
    prob_a1 = mutListA->getLikelihood(currSite, true);
    prob_a0 = mutListA->getLikelihood(currSite, false);
    prob_b1 = mutListB->getLikelihood(currSite, true);
    prob_b0 = mutListB->getLikelihood(currSite, false);

    // calc ll for each type of mutation using betabinomial dist
    //prob_s11 = mutListA->getLikelihood(currSite, true)  * mutListB->getLikelihood(currSite, true);
    //prob_s10 = mutListA->getLikelihood(currSite, true)  * mutListB->getLikelihood(currSite, false);
    //prob_s01 = mutListA->getLikelihood(currSite, false) * mutListB->getLikelihood(currSite, true);
    //prob_s00 = mutListA->getLikelihood(currSite, false) * mutListB->getLikelihood(currSite, false);

    prob_s11 = prob_a1 * prob_b1;
    prob_s10 = prob_a1 * prob_b0;
    prob_s01 = prob_a0 * prob_b1;
    prob_s00 = prob_a0 * prob_b0;
    // mutPairMissing
    // if A doesn't have any data, then prob of both being alt or A alone being alt is 0
    if(mutListA->coordNumRefReadsMap->count(currSite) == 0) {
      prob_s11 = 0;
      prob_s10 = 0;
    }
    if(mutListB->coordNumRefReadsMap->count(currSite) == 0) {
      prob_s11 = 0;
      prob_s01 = 0;
    }
    //std::cout << "mutPairProb: " << prob_a1 << ", " << prob_a0 << ", " << prob_b1 << ", " << prob_b0 << " ::: " << prob_s11 << ", " << prob_s10 << ", " << prob_s01 << ", " << prob_s00 << std::endl;

    //std::cout << "prob_s11: " << prob_s11 << ", " << prob_s10 << ", " << prob_s01 << ", " << prob_s00 << std::endl;
    // try normalizing by total prob
    totalProb = prob_s11 + prob_s10 + prob_s01 + prob_s00;
    if(compareDoubles(0, totalProb)) {
      currLl = 0;
    }
    else {
      prob_s11 /= totalProb;
      prob_s10 /= totalProb;
      prob_s01 /= totalProb;
      prob_s00 /= totalProb;
      //std::cout << "totalProb: " << totalProb << ", " << prob_s11 << ", " << prob_s10 << ", " << prob_s01 << ", " << prob_s00 << std::endl;
      //std::cout << "prob_s: " <<
      //  mutListA->getLikelihood(currSite, true) << " * " << mutListB->getLikelihood(currSite, true) << ", " <<
      //  mutListA->getLikelihood(currSite, true) << " * " << mutListB->getLikelihood(currSite, false) << ", " <<
      //  mutListA->getLikelihood(currSite, false) << " * " << mutListB->getLikelihood(currSite, true) << ", " <<
      //  mutListA->getLikelihood(currSite, false) << " * " << mutListB->getLikelihood(currSite, false) << ", " << std::endl;
      //std::cout << "prob_s: " << prob_s11 << ", " << prob_s10 << ", " << prob_s01 << ", " << prob_s00  << "; " << x1/N << ", " << x2/N << ", " << x3/N << ", " << (N-x1-x2-x3)/N << "; " << prob_s11 * x1/N + prob_s10 * x2/N + prob_s01 * x3/N + prob_s00 * (N-x1-x2-x3)/N << std::endl;

      // then sum up and store
      currLl = log(prob_s11 * x1/N + prob_s10 * x2/N + prob_s01 * x3/N + prob_s00 * (N-x1-x2-x3)/N); // P(D|S=(1,1)) * X1/N + P(D|S=(1,0)) * X2/N  + ...
      //std::cout << "currLl: " << currLl << std::endl;
      if(gsl_isinf(currLl) || gsl_isnan(currLl)) {
        //return GSL_NAN; // Sun 10 Jul 2022 01:38:03 PM PDT test?? ==> Tue 12 Jul 2022 03:31:07 PM PDT should probably be 0; inf comes from taking the log of 0 (from sites with no data)
        currLl = 0;
      }
    }
    totalLl += currLl;
  }
  return totalLl;
}

void MutationPair::print(FILE* stream) {
  MutationListContainer::print(stream);

  /*fprintf(stream, "\nAnalyzed mutation coords:\n");
  for(std::vector<std::string>::iterator itr = sharedMutationVec->begin(); itr != sharedMutationVec->end(); itr++) {
    fprintf(stream, "%s\n", (*itr).c_str());
  }*/
}

// check if x1+x2 > mutListA's length, or x1+x3 > mutListB's, or if x1+x2+x3 > N
double MutationPair::checkOptimProbValidity(gsl_vector* probs) const {
  // check if any mutation count is negative
  double curr = 0;
  double validMax = this->getValidOptimParamMax();
  for(unsigned int i = 0; i < probs->size; i++) {
    curr = gsl_vector_get(probs, i);
    //if(curr < 1e-2) {
    if(curr < 0) {
      //std::cout << "NAN: curr " << curr << " < 1e-2" << std::endl;
      return GSL_NAN;
    }
    if(curr > validMax) {
      //std::cout << "NAN: curr " << curr << " > " << validMax << std::endl;
      return GSL_NAN;
    }
    if(gsl_isnan(curr)) {
      //std::cout << "NAN: curr is nan" << std::endl;
      return GSL_NAN;
    }
  }
  double paramSum = gsl_blas_dasum(probs);
  if(paramSum > validMax) {
    //std::cout << "NAN: paramSum " << paramSum << " > " << validMax << std::endl;
    return GSL_NAN;
  }

  // check if x1+x2 > mutListA's size, or x1+x3 > mutListB's size
  MutationList* mutListA = (*this->mutListVec)[0];
  MutationList* mutListB = (*this->mutListVec)[1];
  double x1 = this->getMutCountEst(0);
  double x2 = this->getMutCountEst(1);
  double x3 = this->getMutCountEst(2);
  if(x1 + x2 > mutListA->coordVec->size() || x1 + x3 > mutListB->coordVec->size()) {
    return GSL_NAN;
  }
  return 0;
}




// ex estimated x1+x2+x3 shouldn't be greater than N
double MutationPair::getValidOptimParamMax() const {
  return sharedMutationVec->size();
  //return std::numeric_limits<int>::max();
}

