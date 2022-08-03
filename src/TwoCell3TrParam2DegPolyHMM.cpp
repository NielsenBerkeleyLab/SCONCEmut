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
//TwoCell3TrParam2DegPolyHMM::TwoCell3TrParam2DegPolyHMM(std::vector<DepthPair*>* depths, gsl_vector* fixedParams, int maxPloidy, int numTrParamsToEst, int numFixedTrParams, int numFixedLibs) : HMM(depths, fixedParams, numTrParamsToEst, 2, 3, maxPloidy, numFixedTrParams, numFixedLibs) { // 2 cells, 3 branches
TwoCell3TrParam2DegPolyHMM::TwoCell3TrParam2DegPolyHMM(std::vector<DepthPair*>* depths, gsl_vector* fixedParams, int maxPloidy, int numTrParamsToEst, int numFixedTrParams, int numFixedLibs, int numBranches, bool preallocIntermediates) : HMM(depths, fixedParams, numTrParamsToEst, 2, numBranches, maxPloidy, numFixedTrParams, numFixedLibs, preallocIntermediates) { // 2 cells
  this->logFacKVec = nullptr;
  this->maxNumBFGSStarts = 3;
  //this->initGuessCases = new std::vector<int>();
  this->setAlpha(0.1); // arbitrary starting value

  /*this->stateRNG = nullptr;
  this->diploidRNG = nullptr;
  this->tumor0RNG = nullptr;
  this->tumor1RNG = nullptr;*/
  this->timeDepMatrixA  = gsl_matrix_alloc(this->timeDepMatrixP->size1, this->timeDepMatrixP->size2);
  this->timeDepMatrixP2  = gsl_matrix_alloc(this->timeDepMatrixP->size1, this->timeDepMatrixP->size2);
  this->timeDepMatrixP3  = gsl_matrix_alloc(this->timeDepMatrixP->size1, this->timeDepMatrixP->size2);
  gsl_matrix_set_zero(this->timeDepMatrixA);
  gsl_matrix_set_zero(this->timeDepMatrixP2);
  gsl_matrix_set_zero(this->timeDepMatrixP3);
}

//TwoCell3TrParam2DegPolyHMM::TwoCell3TrParam2DegPolyHMM(const TwoCell3TrParam2DegPolyHMM& otherHMM) : HMM(otherHMM) {
//  if(otherHMM.logFacKVec != nullptr) {
//    this->logFacKVec = new std::vector<double>(*otherHMM.logFacKVec);
//  }
//  else {
//    this->logFacKVec = nullptr;
//  }
//  this->maxNumBFGSStarts = otherHMM.maxNumBFGSStarts;
//  //this->initGuessCases = new std::vector<int>();
//}
TwoCell3TrParam2DegPolyHMM::~TwoCell3TrParam2DegPolyHMM() {
  delete this->logFacKVec;
  //delete this->initGuessCases;
  /*delete this->stateRNG;
  delete this->diploidRNG;
  delete this->tumor0RNG;
  delete this->tumor1RNG;*/

  // least squares variables
  //gsl_matrix_free(this->coefs);
  //gsl_vector_free(this->residuals);
  //gsl_vector_free(this->a_ij);
  gsl_matrix_free(this->baumWelchTransitionMat);

  gsl_matrix_free(this->timeDepMatrixA);
  gsl_matrix_free(this->timeDepMatrixP2);
  gsl_matrix_free(this->timeDepMatrixP3);
}

//int TwoCell3TrParam2DegPolyHMM::getMaxNumBFGSStarts() const {
//  return this->maxNumBFGSStarts;
//}

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
//double TwoCell3TrParam2DegPolyHMM::getAlpha() const {
//  return gsl_vector_get(this->fixedParams, this->FIXED_TRANSITION_PROB_START_IDX);
//}
//void TwoCell3TrParam2DegPolyHMM::setAlpha(double alpha) {
//  gsl_vector_set(this->fixedParams, this->FIXED_TRANSITION_PROB_START_IDX, alpha);
//}

/*
 ********
 * functions that depend on numbering and ordering of transition params
 ********
 */
double TwoCell3TrParam2DegPolyHMM::setTransition(gsl_matrix* dest, gsl_vector* transitionParams) {
  //std::cout << "in TwoCell3TrParam2DegPolyHMM::setTransition" << std::endl;
  /*if(this->NUM_TRANSITION_PARAMS_TO_EST != 5) {
    std::cerr << "ERROR: incorrect number of transition parameters to estimate. Exiting" << std::endl;
    exit(EXIT_FAILURE);
  }*/
  //std::cout << "setting Transition mat" << std::endl;
  //printRowVector(transitionParams);
  //double alpha = gsl_vector_get(transitionParams, 0); // TODO consts here. does it make sense to have a variably sized vector of consts to define ordering of these params?
  //double beta  = gsl_vector_get(transitionParams, 1);
  //double gamma = gsl_vector_get(transitionParams, 2);
  //double t2 = gsl_vector_get(transitionParams, 3);
  //double t3 = gsl_vector_get(transitionParams, 4);
  //double t1 = 1; // fix t1 at 1

  // Thu 19 Dec 2019 01:10:41 PM PST fixing alpha instead of t1
  //double alpha = gsl_vector_get(transitionParams, 0); // TODO consts here. does it make sense to have a variably sized vector of consts to define ordering of these params?
  double beta  = gsl_vector_get(transitionParams, 0);
  //double gamma = gsl_vector_get(transitionParams, 1);
  double lambda = gsl_vector_get(transitionParams, 1);
  double t1 = gsl_vector_get(transitionParams, 2);
  double t2 = gsl_vector_get(transitionParams, 3);
  double t3 = gsl_vector_get(transitionParams, 4);

  if(dest == this->transition) {
  // save into paramsToEst
  int transitionParamsIdx = 0;
  for(int i = 0; i < this->NUM_TRANSITION_PARAMS_TO_EST; i++) {
    if((unsigned int) (this->TRANSITION_PROB_START_IDX + i) < this->paramsToEst->size) { // Tue 02 Jul 2019 04:32:25 PM PDT TODO fix this; was init added to make the fixAllButBeta thing work
      gsl_vector_set(this->paramsToEst, this->TRANSITION_PROB_START_IDX + i, gsl_vector_get(transitionParams, transitionParamsIdx));
      transitionParamsIdx++;
    }
  }
  //printColVector(transitionParams);
  //std::cout << "paramsToEst: " << std::endl;
  //printColVector(this->paramsToEst);
  //std::cout << "this->NUM_BRANCH_LENGTHS_TO_EST: " << this->NUM_BRANCH_LENGTHS_TO_EST << std::endl;
  for(int i = 0; i < this->NUM_BRANCH_LENGTHS_TO_EST; i++) {
    //std::cout << "i: " << i << ", of " << this->NUM_BRANCH_LENGTHS_TO_EST << ", transitionParamsIdx: " << transitionParamsIdx << std::endl;
    gsl_vector_set(this->paramsToEst, this->BRANCH_LENGTH_START_IDX + i, gsl_vector_get(transitionParams, transitionParamsIdx));
    transitionParamsIdx++;
  }
  }

  //return this->setTransition(this->alpha, beta, gamma, t1, t2, t3);
  //return this->setTransition(dest, this->getAlpha(), beta, gamma, t1, t2, t3);
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
//double TwoCell3TrParam2DegPolyHMM::setTransition(gsl_matrix* dest, double alpha, double beta, double gamma, double t1, double t2, double t3) {
double TwoCell3TrParam2DegPolyHMM::setTransition(gsl_matrix* dest, double alpha, double beta, double lambda, double t1, double t2, double t3) {
  // first set rate matrix, if it's possible for beta and lambda to change
  this->setRateMatrixQ(alpha, beta, lambda);
  //std::cout << "rate matrix q" << std::endl;
  //printMatrix(this->rateMatrixQ);

  // then set time dependent matrix A for times t1,t2,t3
  double status = this->setTimeDepMatrixA(t1, t2, t3);
  //std::cout << "time dep matrix A" << std::endl;
  //printMatrix(this->timeDepMatrixA);

  // then calculate matrix M(t) = {m_(i,j),(i',j')(\bar t)}, where m_(i,j),(i',j')(\bar t) = A_(i,j),(i',j')(t) / sum_(i'',j'') A_(i,j),(i'',j'')
  // this is just normalizing matrix A by each rowsum
  double rowsum = 0;
  gsl_matrix_memcpy(dest, this->timeDepMatrixA);
  for(unsigned int row = 0; row < dest->size1; row++) {
    // rescale by rowsum
    gsl_vector_view currRow = gsl_matrix_row(dest, row);

    //std::cout << "unscaled row " << row << " in transition mat" << std::endl;
    //printRowVector(&currRow.vector);

    rowsum = gsl_blas_dasum(&currRow.vector); // Double Absolute SUM
    gsl_vector_scale(&currRow.vector, 1.0 / rowsum);

    //std::cout << "scaled row " << row << " in transition mat" << std::endl;
    //printRowVector(&currRow.vector);
  }
  //std::cout << "scaled dest" << std::endl;
  //printMatrix(dest);

  return status;


  ///*// if transition matrix doesn't exist yet, create it
  //if(this->transition == nullptr) {
  //  this->transition = gsl_matrix_alloc(this->states->size(), this->states->size());
  //  this->allocIntermediates(); // intermediates probably weren't set if transition matrix is null
  //}*/ // Fri 26 Jun 2020 03:51:06 PM PDT now saving into dest instead of this->transition by default

  //// if haven't saved the string represetation of the transition matrix yet, save it now
  //bool saveStrs = false;
  //std::string* currTrStr = nullptr;
  //if(this->transitionStrings == nullptr) {
  //  saveStrs = true;
  //  this->transitionStrings = new std::vector<std::vector<std::string*>*>(this->states->size());
  //  for(unsigned int i = 0; i < this->states->size(); i++) {
  //    (*this->transitionStrings)[i] = new std::vector<std::string*>(this->states->size());
  //    for(unsigned int j = 0; j < this->states->size(); j++) {
  //      currTrStr = new std::string();
  //      *currTrStr += "(";
  //      (*(*this->transitionStrings)[i])[j] = currTrStr;
  //    }
  //  }
  //}

  //// zero out everything
  ////gsl_matrix_set_zero(this->transition);
  //gsl_matrix_set_zero(dest);

  ////////////////////////////////////////////////////
  ////// start iterating over the state space
  ////int frA = 0;
  ////int frB = 0;
  ////int toA = 0;
  ////int toB = 0;
  ////double currProb = 0;
  ////for(unsigned int row = 0; row < this->states->size(); row++) {
  ////  frA = getCellPloidyFromStateIdx(0, row);
  ////  frB = getCellPloidyFromStateIdx(1, row);
  ////  for(unsigned int col = 0; col < this->states->size(); col++) {
  ////    currTrStr = (*(*this->transitionStrings)[row])[col];
  ////    currProb = 0;
  ////    toA = getCellPloidyFromStateIdx(0, col);
  ////    toB = getCellPloidyFromStateIdx(1, col);

  ////    // if frA == toA (A does not move). These are the main diagonal blocks
  ////    if(frA == toA) {
  ////      // if frB == toB (B does not move)
  ////      if(frB == toB) {
  ////        // diagonal; continue, will set to normalizing constants later
  ////        if(saveStrs) {
  ////          *currTrStr += "r)";
  ////        }
  ////        continue;
  ////      }
  ////      currProb += beta;  // any CNA for B
  ////      if(saveStrs) {
  ////        *currTrStr += "b";
  ////      }
  ////      // if abs(frB - toB) == 1 (B adjacent move)
  ////      if(std::abs(frB - toB) == 1) {
  ////        currProb += alpha; // adj CNA
  ////        if(saveStrs) {
  ////          *currTrStr += "+a";
  ////        }
  ////      }
  ////      // if toB == 2
  ////      if(toB == 2) {
  ////        currProb += gamma; // return to diploid
  ////        if(saveStrs) {
  ////          *currTrStr += "+g";
  ////        }
  ////      }
  ////      // set transition[row][col] = t3 * currProb
  ////      currProb = t3 * currProb;
  ////      if(saveStrs) {
  ////        *currTrStr += ")*t3";
  ////      }
  ////    }
  ////    // else (A has moved). These are the main off diagonal blocks
  ////    else {
  ////      // if frB == toB (B does not move)
  ////      if(frB == toB) {
  ////        currProb += beta; // any CNA for A
  ////        if(saveStrs) {
  ////          *currTrStr += "b";
  ////        }
  ////        // if abs(frA - toA) == 1 (A adjacent move)
  ////        if(std::abs(frA - toA) == 1) {
  ////          currProb += alpha; // adj CNA
  ////          if(saveStrs) {
  ////            *currTrStr += "+a";
  ////          }
  ////        }
  ////        // if toA == 2
  ////        if(toA == 2) {
  ////          currProb += gamma; // return to diploid
  ////          if(saveStrs) {
  ////            *currTrStr += "+g";
  ////          }
  ////        }
  ////        // set transition[row][col] = t2 * currProb
  ////        currProb = t2 * currProb;
  ////        if(saveStrs) {
  ////          *currTrStr += ")*t2";
  ////        }
  ////      }
  ////      // else if frA == frB && toA == toB (move in parallel)
  ////      else if(frA == frB && toA == toB) {
  ////        currProb += beta;
  ////        if(saveStrs) {
  ////          *currTrStr += "b";
  ////        }
  ////        // if abs(frA - toA) == 1 (A and B adjacent move)
  ////        if(std::abs(frA - toA) == 1) {
  ////          currProb += alpha;
  ////          if(saveStrs) {
  ////            *currTrStr += "+a";
  ////          }
  ////        }
  ////        // if toA == 2
  ////        if(toA == 2) {
  ////          currProb += gamma; // return to diploid
  ////          if(saveStrs) {
  ////            *currTrStr += "+g";
  ////          }
  ////        }
  ////        // set transition[row][col] = t1 * currProb; continue
  ////        currProb = t1 * currProb;
  ////        if(saveStrs) {
  ////          *currTrStr += ")*t1";
  ////        }
  ////      }
  ////      // else, some impossible transition
  ////      else {
  ////        // set transition[row][col] = 0
  ////        currProb = 0;
  ////        if(saveStrs) {
  ////          *currTrStr += "0)";
  ////        }
  ////      }
  ////    }
  ////    //std::cerr << "built: " << *currTrStr << std::endl;
  ////    //std::cerr << "stored: " << (*(*this->transitionStrings)[row])[col] << std::endl;
  ////    //fprintf(stderr, "stored: %s\n",  (*(*this->transitionStrings)[row])[col]->c_str());
  ////    // store transition prob
  ////    //std::cerr << frA << ", " << frB << ", " << toA << ", " << toB << ", " << currProb << std::endl;
  ////    //gsl_matrix_set(this->transition, row, col, currProb);
  ////    gsl_matrix_set(dest, row, col, currProb);
  ////  }
  ////}  //////////////////////////////////////////////////

  //// staggered version (only +-1 transitions allowed), gamma*t1 terms all dropped
  //// start iterating over the state space
  //int frA = 0;
  //int frB = 0;
  //int toA = 0;
  //int toB = 0;
  //double currProb = 0;
  //for(unsigned int row = 0; row < this->states->size(); row++) {
  //  frA = getCellPloidyFromStateIdx(0, row);
  //  frB = getCellPloidyFromStateIdx(1, row);
  //  for(unsigned int col = 0; col < this->states->size(); col++) {
  //    currTrStr = (*(*this->transitionStrings)[row])[col];
  //    currProb = 0;
  //    toA = getCellPloidyFromStateIdx(0, col);
  //    toB = getCellPloidyFromStateIdx(1, col);

  //    // if frA == toA (A does not move). These are the main diagonal blocks
  //    if(frA == toA) {
  //      // if frB == toB (B does not move)
  //      if(frB == toB) {
  //        // diagonal; continue, will set to normalizing constants later
  //        if(saveStrs) {
  //          *currTrStr += "r)";
  //        }
  //        continue;
  //      }
  //      currProb += beta;  // any CNA for B
  //      if(saveStrs) {
  //        *currTrStr += "b";
  //      }
  //      // if abs(frB - toB) == 1 (B adjacent move)
  //      if(std::abs(frB - toB) == 1) {
  //        currProb += alpha; // adj CNA
  //        if(saveStrs) {
  //          *currTrStr += "+a";
  //        }
  //      }
  //      // if toB == 2
  //      if(toB == 2) {
  //        currProb += gamma; // return to diploid
  //        if(saveStrs) {
  //          *currTrStr += "+g";
  //        }
  //      }
  //      // special case for 0th and kth blocks, where losses can go to "0 or less" or "k or more". should be (b+a[+g])*(t1+t3)
  //      // "0 or less" is on subdiagonal; "k or more" is on superdiagonal
  //      if((toA == 0 && frB - toB == 1) || (toA == this->MAX_PLOIDY && toB - frB == 1)) {
  //        if(toB != 2) {
  //          currProb *= (t1+t3);
  //          if(saveStrs) {
  //            *currTrStr += ")*(t1+t3)";
  //          }
  //        }
  //        else { // 0,3 -> 0,2; MAX_PLOIDY,1 -> MAX_PLOIDY,2 cases, to get (b+a+g)*(t3)+(b+a)*t1
  //          currProb *= t3;
  //          currProb += (beta + alpha) * t1;
  //          if(saveStrs) {
  //            *currTrStr += ")*t3+(b+a)*t1";
  //          }
  //        }
  //      }
  //      // else, just multiply by t3
  //      // set transition[row][col] = t3 * currProb
  //      else {
  //        currProb *= t3;
  //        if(saveStrs) {
  //          *currTrStr += ")*t3";
  //        }
  //      }
  //    }
  //    // else (A has moved). These are the main off diagonal blocks
  //    else {
  //      // if frB == toB (B does not move); diagonal (b[+a+g]*t2 elements
  //      if(frB == toB) {
  //        currProb += beta; // any CNA for A
  //        if(saveStrs) {
  //          *currTrStr += "b";
  //        }
  //        // if abs(frA - toA) == 1 (A adjacent move)
  //        if(std::abs(frA - toA) == 1) {
  //          currProb += alpha; // adj CNA
  //          if(saveStrs) {
  //            *currTrStr += "+a";
  //          }
  //        }
  //        // if toA == 2
  //        if(toA == 2) {
  //          currProb += gamma; // return to diploid
  //          if(saveStrs) {
  //            *currTrStr += "+g";
  //          }
  //        }
  //        // 0 or less/k or more corners should be (b+a)*(t1+t2) or (b+a+g)*t2+(b+a)*t1
  //        // sub main diagonals are 0 or less
  //        // super main diagonals are k or more
  //        //if((toB == 0 || toB == this->MAX_PLOIDY) && std::abs(toA - frA) == 1) {
  //        if((toB == 0  && frA - toA == 1) || (toB == this->MAX_PLOIDY && toA - frA == 1)) {
  //          if(toA != 2) {
  //            currProb *= (t1+t2);
  //            if(saveStrs) {
  //              *currTrStr += ")*(t1+t2)";
  //            }
  //          }
  //          else {
  //            currProb *= t2;
  //            currProb += (beta + alpha) * t1;
  //            if(saveStrs) {
  //              *currTrStr += ")*t2+(b+a)*t1";
  //            }
  //          }
  //        }
  //        // set transition[row][col] = t2 * currProb
  //        else {
  //          currProb *= t2;
  //          if(saveStrs) {
  //            *currTrStr += ")*t2";
  //          }
  //        }
  //      }
  //      // t1 elements on sub/super diagonals on off diagonal blocks. (b+a)*t1 elements
  //      //else if(std::abs(frA - toA) == 1 && std::abs(frB - toB) == 1) {
  //      else if((frA - toA == 1 && frB - toB == 1) || (toA - frA == 1 && toB - frB == 1)) {
  //        currProb = (beta + alpha) * t1;
  //        if(saveStrs) {
  //          *currTrStr += "b+a)*t1";
  //        }
  //      }

  //      // else, some impossible transition
  //      else {
  //        // set transition[row][col] = 0
  //        currProb = 0;
  //        if(saveStrs) {
  //          *currTrStr += "0)";
  //        }
  //      }
  //    }
  //    //std::cerr << "built: " << *currTrStr << std::endl;
  //    //std::cerr << "stored: " << (*(*this->transitionStrings)[row])[col] << std::endl;
  //    //fprintf(stderr, "stored: %s\n",  (*(*this->transitionStrings)[row])[col]->c_str());
  //    // store transition prob
  //    //std::cerr << frA << ", " << frB << ", " << toA << ", " << toB << ", " << currProb << std::endl;
  //    //gsl_matrix_set(this->transition, row, col, currProb);
  //    gsl_matrix_set(dest, row, col, currProb);
  //  }
  //}//////////////////////////////////////////////////
  ////    // if frA == toA (A does not move). These are the main diagonal blocks
  ////    if(frA == toA) {
  ////      // if frB == toB (B does not move)
  ////      if(frB == toB) {
  ////        // diagonal; continue, will set to normalizing constants later
  ////        if(saveStrs) {
  ////          *currTrStr += "r)";
  ////        }
  ////        continue;
  ////      }
  ////      currProb += beta;  // any CNA for B
  ////      if(saveStrs) {
  ////        *currTrStr += "b";
  ////      }
  ////      // if abs(frB - toB) == 1 (B adjacent move)
  ////      if(std::abs(frB - toB) == 1) {
  ////        currProb += alpha; // adj CNA
  ////        if(saveStrs) {
  ////          *currTrStr += "+a";
  ////        }
  ////      }
  ////      // if toB == 2
  ////      if(toB == 2) {
  ////        currProb += gamma; // return to diploid
  ////        if(saveStrs) {
  ////          *currTrStr += "+g";
  ////        }
  ////      }
  ////      // set transition[row][col] = t3 * currProb
  ////      currProb = t3 * currProb;
  ////      if(saveStrs) {
  ////        *currTrStr += ")*t3";
  ////      }
  ////    }
  ////    // else (A has moved). These are the main off diagonal blocks
  ////    else {
  ////      // if frB == toB (B does not move)
  ////      if(frB == toB) {
  ////        currProb += beta; // any CNA for A
  ////        if(saveStrs) {
  ////          *currTrStr += "b";
  ////        }
  ////        // if abs(frA - toA) == 1 (A adjacent move)
  ////        if(std::abs(frA - toA) == 1) {
  ////          currProb += alpha; // adj CNA
  ////          if(saveStrs) {
  ////            *currTrStr += "+a";
  ////          }
  ////        }
  ////        // if toA == 2
  ////        if(toA == 2) {
  ////          currProb += gamma; // return to diploid
  ////          if(saveStrs) {
  ////            *currTrStr += "+g";
  ////          }
  ////        }
  ////        // set transition[row][col] = t2 * currProb
  ////        currProb = t2 * currProb;
  ////        if(saveStrs) {
  ////          *currTrStr += ")*t2";
  ////        }
  ////      }
  ////      // else if frA == frB && toA == toB (move in parallel)
  ////      else if(frA == frB && toA == toB) {
  ////        currProb += beta;
  ////        if(saveStrs) {
  ////          *currTrStr += "b";
  ////        }
  ////        // if abs(frA - toA) == 1 (A and B adjacent move)
  ////        if(std::abs(frA - toA) == 1) {
  ////          currProb += alpha;
  ////          if(saveStrs) {
  ////            *currTrStr += "+a";
  ////          }
  ////        }
  ////        // if toA == 2
  ////        if(toA == 2) {
  ////          currProb += gamma; // return to diploid
  ////          if(saveStrs) {
  ////            *currTrStr += "+g";
  ////          }
  ////        }
  ////        // set transition[row][col] = t1 * currProb; continue
  ////        currProb = t1 * currProb;
  ////        if(saveStrs) {
  ////          *currTrStr += ")*t1";
  ////        }
  ////      }

  ////      // else if frA == frB && toA == toB (move in parallel); ex 0,0 -> 3,3
  ////      else if(frA == frB && toA == toB) {
  ////        currProb += beta;
  ////        if(saveStrs) {
  ////          *currTrStr += "b";
  ////        }
  ////        // if abs(frA - toA) == 1 (A and B adjacent move)
  ////        if(std::abs(frA - toA) == 1) {
  ////          currProb += alpha;
  ////          if(saveStrs) {
  ////            *currTrStr += "+a";
  ////          }
  ////        }
  ////        // if toA == 2
  ////        if(toA == 2) {
  ////          currProb += gamma; // return to diploid
  ////          if(saveStrs) {
  ////            *currTrStr += "+g";
  ////          }
  ////        }
  ////        // set transition[row][col] = t1 * currProb; continue
  ////        currProb = t1 * currProb;
  ////        if(saveStrs) {
  ////          *currTrStr += ")*t1";
  ////        }
  ////      }
  ////      ////////////////////////////////////////////////

  ////      ////////////////////////////////////////////////
  ////      // same additive change in tandem
  ////      else if((frA - toA) == (frB - toB)) {
  ////        currProb += beta;
  ////        if(saveStrs) {
  ////          *currTrStr+= "b";
  ////        }
  ////        // if abs(frA - toA) == 1 (A and B adjacent move)
  ////        if(std::abs(frA - toA) == 1) {
  ////          currProb += alpha;
  ////          if(saveStrs) {
  ////            *currTrStr += "+a";
  ////          }
  ////        }
  ////        // if either goes to diploid
  ////        if(toA == 2 || toB == 2) {
  ////          currProb += gamma; // return to diploid
  ////          if(saveStrs) {
  ////            *currTrStr += "+g";
  ////          }
  ////        }
  ////        // set transition[row][col] = t1 * currProb; continue
  ////        currProb = t1 * currProb;
  ////        if(saveStrs) {
  ////          *currTrStr += ")*t1";
  ////        }
  ////      }
  ////      //////////////////////////////////////////////// 





  //// set diagonals to be normalizing constants
  //for(unsigned int row = 0; row < this->states->size(); row++) {
  //  double rowSum = 0;
  //  for(unsigned int col = 0; col < this->states->size(); col++) {
  //    //rowSum += gsl_matrix_get(this->transition, row, col);
  //    rowSum += gsl_matrix_get(dest, row, col);
  //  }
  //  //gsl_matrix_set(this->transition, row, row, 1 - rowSum);
  //  gsl_matrix_set(dest, row, row, 1 - rowSum);
  //}

  //// check if anything went negative
  //// if somehow the transition matrix went negative, return NAN
  ////double transitionMatMin = gsl_matrix_min(this->transition);
  //double transitionMatMin = gsl_matrix_min(dest);
  //if(transitionMatMin < 0 || gsl_isnan(transitionMatMin)) {
  //  return GSL_NAN;
  //}

  // otherwise, update initProb
  //return this->setInitProbSteadyState();
  //return 0;
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

  //std::cout << "timeDepMatrixP" << std::endl;
  //printMatrix(this->timeDepMatrixP);
  //std::cout << "timeDepMatrixP2" << std::endl;
  //printMatrix(this->timeDepMatrixP2);
  //std::cout << "timeDepMatrixP3" << std::endl;
  //printMatrix(this->timeDepMatrixP3);

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
      //std::cout << "Calculating for (row:i,j),(col:i',j'): (" << row << ":" << i << "," << j << "),(" << col << ":" << iPrime << "," << jPrime << ")" << std::endl;
      for(int L = 0; L <= this->MAX_PLOIDY; L++) {
        for(int V = 0; V <= this->MAX_PLOIDY; V++) {
          lvStateIdx = getStateIdxFromPloidyPair(L, V);
          currLVProd = gsl_matrix_get(this->timeDepMatrixP, ancDiploidStateIdx, lvStateIdx);
          currLVProd *= gsl_matrix_get(this->timeDepMatrixP2, lvStateIdx, getStateIdxFromPloidyPair(i, iPrime));
          currLVProd *= gsl_matrix_get(this->timeDepMatrixP3, lvStateIdx, getStateIdxFromPloidyPair(j, jPrime));
          currStateSum += currLVProd;
          //std::cout << "P_(" << ancDiploidStateIdx << ":2,2),(" << lvStateIdx << ":" << L << "," << V << ")(t1) * " <<
          //             "P_(" << lvStateIdx << ":" << L << "," << V << "),(" << getStateIdxFromPloidyPair(i, iPrime) << ":" << i << "," << iPrime << ")(t2) * " <<
          //             "P_(" << lvStateIdx << ":" << L << "," << V << "),(" << getStateIdxFromPloidyPair(j, jPrime) << ":" << j << "," << jPrime << ")(t3) = " << currLVProd << ", runningTotal: " << currStateSum << std::endl;
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
 *   For the transition matrix, the constraints are
 *   0 <= (dt1 + dt2 + dt3) * (2a + kb + g) <= 1
 *   (dt1 + dt2 + dt3) <= 1
 *   (2a + kb + g) <= 1
 *   where d = max_j (d_j) = max window size = (*this->depths)[0]->avgTumorDepth->maxWindowSize
 *   set T1 = 1, dT2 = dt2 / (1 + dt2 + dt3), dT3 = dt3 / (1 + dt2 + dt3) // Fri 05 Jul 2019 04:56:42 PM PDT questioning, may be backwards
 *
 * //The following transformation converts params x,y,z to probabilities a,b,g
 * //  2a = e^x / (1 + e^x + e^y + e^z) ==> a = e^x / [2 (1 + e^x + e^y + e^z)]
 * //  kb = e^y / (1 + e^x + e^y + e^z) ==> b = e^y / [k (1 + e^x + e^y + e^z)]
 * //  g = e^z / (1 + e^x + e^y + e^z) ==> g = e^z / [(1 + e^x + e^y + e^z)]
 *
 * //The following transformation converts probs a,b,g to params x,y,z
 * //  Let c = 1+ e^x + e^y + e^z
 * //  2a + kb + g + (1-2a-kb-g) = 1 ==> e^x / c + e^y / c + e^z / c + 1 / c = 1
 * //    ==> (1-2a-kb-g) = 1 / c
 * //    ==> c = 1 / (1-2a-kb-g)
 * //  a = e^x / (2c) ==> a * 2c = e^x
 * //  x = ln(a * 2c)
 * //  y = ln(b * kc)
 * //  z = ln(g * c)
 * //  checked with https://www.wolframalpha.com/input/?i=Solve%5B%7B2a+%3D%3D+Exp%5Bx%5D%2F(1+%2B+Exp%5Bx%5D+%2B+Exp%5By%5D+%2B+Exp%5Bz%5D),+++k*b+%3D%3D+Exp%5By%5D%2F(1+%2B+Exp%5Bx%5D+%2B+Exp%5By%5D+%2B+Exp%5Bz%5D),+++g+%3D%3D+Exp%5Bz%5D%2F(1+%2B+Exp%5Bx%5D+%2B+Exp%5By%5D+%2B+Exp%5Bz%5D)%7D,+%7Bx,+y,+z%7D%5D
 *
 * //The following transformation converts params T2,T3 to branch lengths t2,t3
 * //  dt1 = 1   / (1 + e^(T2) + e^(T3))
 * //  dt2 = e^(T2) / (1 + e^(T2) + e^(T3))
 * //  dt3 = e^(T3) / (1 + e^(T2) + e^(T3))
 *
 * //The following transformation converts branch lengths t2,t3 to params T2,T3
 * //  T2 = ln(-(d * t2) / (d * (t2 + t3) - 1)
 * //  T3 = ln(-(d * t3) / (d * (t2 + t3) - 1)
 * //  as determined by https://www.wolframalpha.com/input/?i=Solve%5B%7B+d*s+%3D%3D+Exp%5BA%5D%2F(1+%2B+Exp%5BA%5D+%2B+Exp%5BB%5D),++++d*t+%3D%3D+Exp%5BB%5D%2F(1+%2B+Exp%5BA%5D+%2B+Exp%5BB%5D)%7D,+%7BA,+B%7D%5D
 * //    (s == t2, t == t3, A == T2, B == T3 to make wolfram alpha happy)
 *

 * //Fri 05 Jul 2019 04:53:25 PM PDT
 * //if we're fixing t1 at 1 (true branch length), then are we missing a d term? lower t are equiv of a,b,g ==> should be the exp terms.
 * //trying 
 * //   t1 = 1      / (1 + 1/d + e^(T2) + e^(T3))
 * //  dt2 = e^(T2) / (1 + 1/d + e^(T2) + e^(T3))
 * //  dt3 = e^(T3) / (1 + 1/d + e^(T2) + e^(T3))
 * //  T2 = ln(-((d+1) * t2) / (d * (t2 + t3) - 1)
 * //  T3 = ln(-((d+1) * t3) / (d * (t2 + t3) - 1)
 * //https://www.wolframalpha.com/input/?i=Solve%5B%7Bd*+s+%3D%3D+Exp%5BA%5D%2F(1%2B+1%2Fd%2B+Exp%5BA%5D+%2B+Exp%5BB%5D),+d*t+%3D%3D+Exp%5BB%5D%2F(1%2B1%2Fd%2BExp%5BA%5D+%2B+Exp%5BB%5D)%7D,+%7BA,+B%7D%5D
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
//void TwoCell3TrParam2DegPolyHMM::convertProbToParam(gsl_vector* dest, const gsl_vector* src, int k) const {
void TwoCell3TrParam2DegPolyHMM::convertProbToParam(gsl_vector* dest, const gsl_vector* src) const {
  //std::cout << "TwoCell3TrParam2DegPolyHMM::convertProbToParam" << std::endl;
  /*//std::cerr << this->TRANSITION_PROB_START_IDX << ", " << this->LIB_SIZE_SCALING_FACTOR_START_IDX << std::endl;
  //printRowVector((gsl_vector*) src);
  double d = (double) (*this->depths)[0]->maxWindowSize;
  //std::cerr << d << std::endl;
  //double a = gsl_vector_get(src, this->TRANSITION_PROB_START_IDX + 0);
  //double b = gsl_vector_get(src, this->TRANSITION_PROB_START_IDX + 1);
  //double g = gsl_vector_get(src, this->TRANSITION_PROB_START_IDX + 2);
  //double t2 = gsl_vector_get(src, this->TRANSITION_PROB_START_IDX + 3) / d;
  //double t3 = gsl_vector_get(src, this->TRANSITION_PROB_START_IDX + 4) / d;
  //double r = gsl_vector_get(src, this->LIB_SIZE_SCALING_FACTOR_START_IDX);
  //double s = gsl_vector_get(src, this->LIB_SIZE_SCALING_FACTOR_START_IDX + 1);

  double a = this->getAlpha(); //this->alpha;//gsl_vector_get(src, this->TRANSITION_PROB_START_IDX + 0);
  double b = gsl_vector_get(src, this->TRANSITION_PROB_START_IDX + 0);
  double g = gsl_vector_get(src, this->TRANSITION_PROB_START_IDX + 1);
  double t1 = gsl_vector_get(src, this->TRANSITION_PROB_START_IDX + 2) / d;
  double t2 = gsl_vector_get(src, this->TRANSITION_PROB_START_IDX + 3) / d;
  double t3 = gsl_vector_get(src, this->TRANSITION_PROB_START_IDX + 4) / d;
  double r = gsl_vector_get(src, this->LIB_SIZE_SCALING_FACTOR_START_IDX);
  double s = gsl_vector_get(src, this->LIB_SIZE_SCALING_FACTOR_START_IDX + 1);

  //double c = 1.0 / (1.0 - 2.0 * a - (double) this->getKploidy() * b - g);
  double c = (b * this->getKploidy() + g - 1);
  //std::cout << "a, b, g: " << a << ", " << b << ", " << g << std::endl;
  //std::cout << "c: " << c << std::endl;
  //gsl_vector_set(dest, this->TRANSITION_PROB_START_IDX + 0, log(a * 2.0 * c));
  //gsl_vector_set(dest, this->TRANSITION_PROB_START_IDX + 1, log(b * (double) this->getKploidy() * c));
  //gsl_vector_set(dest, this->TRANSITION_PROB_START_IDX + 2, log(g * c));

  gsl_vector_set(dest, this->TRANSITION_PROB_START_IDX + 0, log(-b * (1-2*a) * (double) this->getKploidy() / c)); // set y
  gsl_vector_set(dest, this->TRANSITION_PROB_START_IDX + 1, log(-g * (1-2*a) / c)); // set z

  gsl_vector_set(dest, this->TRANSITION_PROB_START_IDX + 2, log(-(d * t1) / (d * (t1 + t2 + t3) - 1))); // set T1
  gsl_vector_set(dest, this->TRANSITION_PROB_START_IDX + 3, log(-(d * t2) / (d * (t1 + t2 + t3) - 1))); // set T2
  gsl_vector_set(dest, this->TRANSITION_PROB_START_IDX + 4, log(-(d * t3) / (d * (t1 + t2 + t3) - 1))); // set T3
  //gsl_vector_set(dest, this->TRANSITION_PROB_START_IDX + 3, log(-((d+1) * t2) / (d * (t2 + t3) - 1))); // set T2
  //gsl_vector_set(dest, this->TRANSITION_PROB_START_IDX + 4, log(-((d+1) * t3) / (d * (t2 + t3) - 1))); // set T3

  gsl_vector_set(dest, this->LIB_SIZE_SCALING_FACTOR_START_IDX, log(r));
  gsl_vector_set(dest, this->LIB_SIZE_SCALING_FACTOR_START_IDX + 1, log(s));*/


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
//void TwoCell3TrParam2DegPolyHMM::convertParamToProb(gsl_vector* dest, const gsl_vector* src, int k) const {
void TwoCell3TrParam2DegPolyHMM::convertParamToProb(gsl_vector* dest, const gsl_vector* src) const {
  //std::cout << "TwoCell3TrParam2DegPolyHMM::convertParamToProb" << std::endl;
  /*//double x = gsl_vector_get(src, this->TRANSITION_PROB_START_IDX + 0);
  //double y = gsl_vector_get(src, this->TRANSITION_PROB_START_IDX + 1);
  //double z = gsl_vector_get(src, this->TRANSITION_PROB_START_IDX + 2);
  double y = gsl_vector_get(src, this->TRANSITION_PROB_START_IDX + 0);
  double z = gsl_vector_get(src, this->TRANSITION_PROB_START_IDX + 1);
  double T1 = gsl_vector_get(src, this->TRANSITION_PROB_START_IDX + 2);
  double T2 = gsl_vector_get(src, this->TRANSITION_PROB_START_IDX + 3);
  double T3 = gsl_vector_get(src, this->TRANSITION_PROB_START_IDX + 4);
  double w = gsl_vector_get(src, this->LIB_SIZE_SCALING_FACTOR_START_IDX);
  double v = gsl_vector_get(src, this->LIB_SIZE_SCALING_FACTOR_START_IDX + 1);
  //double c = 1 + exp(x) + exp(y) + exp(z);
  double a = this->getAlpha(); //this->alpha;
  double c = 1 - 2.0*a + exp(y) + exp(z);
  //gsl_vector_set(dest, this->TRANSITION_PROB_START_IDX + 0, exp(x) / (2.0 * c));
  //gsl_vector_set(dest, this->TRANSITION_PROB_START_IDX + 1, exp(y) / ((double) this->getKploidy() * c));
  //gsl_vector_set(dest, this->TRANSITION_PROB_START_IDX + 2, exp(z) / c);
  gsl_vector_set(dest, this->TRANSITION_PROB_START_IDX + 0, exp(y) / ((double) this->getKploidy() * c)); // beta
  gsl_vector_set(dest, this->TRANSITION_PROB_START_IDX + 1, exp(z) / c); // gamma
  //c = 1.0 / (1 + exp(T2) + exp(T3));
  c = 1.0 / (1 + exp(T1) + exp(T2) + exp(T3));
  //double d = (double) (*this->depths)[0]->maxWindowSize;
  //c = 1.0 / (1 + 1/d + exp(T2) + exp(T3));

  gsl_vector_set(dest, this->TRANSITION_PROB_START_IDX + 2, exp(T1) * c); // set t1
  gsl_vector_set(dest, this->TRANSITION_PROB_START_IDX + 3, exp(T2) * c); // set t2
  gsl_vector_set(dest, this->TRANSITION_PROB_START_IDX + 4, exp(T3) * c); // set t3

  gsl_vector_set(dest, this->LIB_SIZE_SCALING_FACTOR_START_IDX, exp(w));
  gsl_vector_set(dest, this->LIB_SIZE_SCALING_FACTOR_START_IDX + 1, exp(v));*/


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
//double TwoCell3TrParam2DegPolyHMM::getEmissionProb(double tumorDepth, double diploidDepth, double ploidy, int cellIdx) {
//double TwoCell3TrParam2DegPolyHMM::getEmissionProb(double tumorDepth, double diploidDepth, int ploidy, int windowIdx, int cellIdx) {
//double TwoCell3TrParam2DegPolyHMM::getEmissionProb(double tumorDepth, double diploidDepth, int ploidy, int chrIdx, int depthIdx, int cellIdx) {
double TwoCell3TrParam2DegPolyHMM::getEmissionProb(double tumorDepth, double diploidDepth, int ploidy, int cellIdx) {
  return exp(this->getLogEmissionProb(tumorDepth, diploidDepth, ploidy, cellIdx));
}

/*
 * returns emission probability across all cells for a given state, chromosome, and depthIdx (position along chr).
 * emission prob is NOT in log space.
 * currChrDepthsVec is useful for caching and provides a speed up
 */
double TwoCell3TrParam2DegPolyHMM::getTotalEmissionProb(int stateIdx, std::vector<std::vector<double>*>* currChrDepthsVec, int chrIdx, int depthIdx) {
  //return exp(getTotalLogEmissionProb(stateIdx, currChrDepthsVec, chrIdx, depthIdx)); // marginally slower
  double emissionProb = 0;
  int currPloidy = -1;
  double currTumorDepth = -1;
  double currDiploidDepth = (*(*currChrDepthsVec)[this->NUM_CELLS])[depthIdx];
  // for each cell, calc the emission prob and multiply by it
  for(int cellIdx = 0; cellIdx < this->NUM_CELLS; cellIdx++) {
    currPloidy = getCellPloidyFromStateIdx(cellIdx, stateIdx);
    currTumorDepth = (*(*currChrDepthsVec)[cellIdx])[depthIdx];
    //emissionProb += log(this->getEmissionProb(currTumorDepth, currDiploidDepth, currPloidy, cellIdx)); // extra exp/log
    emissionProb += this->getLogEmissionProb(currTumorDepth, currDiploidDepth, currPloidy, cellIdx);
  }
  return exp(emissionProb);
}

double TwoCell3TrParam2DegPolyHMM::getLogEmissionProb(double tumorDepth, double diploidDepth, int ploidy, int cellIdx) {
  // if diploidDepth == 0, assume this is a hard to map region (ie any tumor reads mapping here are spurious)
  // P(any number reads | no data) = 1
  //if(compareDoubles(diploidDepth, 0.0)) {
  if(diploidDepth < 1e-4) {
    //std::cout << "diploidDepth 0, tumorDepth " << tumorDepth << std::endl;
    //return 1;
    return 0; // return 0 since log(1) = 0
  }

  // Wed 29 May 2019 12:05:41 PM PDT with new error term, shouldn't short circuit here any more
  //// if ploidy == 0, cannot have any reads mapping here (ie there isn't any DNA)
  //// P(0 reads) = 1
  //if(compareDoubles(ploidy, 0.0)) {
  //  if(compareDoubles(tumorDepth, 0.0)) {
  //    //std::cout << "ploidy 0, tumorDepth 0" << std::endl;
  //    return 1;
  //  }
  //  else {
  //    //std::cout << "ploidy 0, tumorDepth " << tumorDepth << std::endl;
  //    return 0;
  //  }
  //}

  // if scaling factor becomes super small or ridiculously large, return error code
  double currLibSizeScalingFactor = this->getLibScalingFactor(cellIdx);
  /*if(currLibSizeScalingFactor < 5e-3 || currLibSizeScalingFactor > 1e3) {
    return GSL_NAN;
  }*/ // removed Tue 05 Nov 2019 05:38:55 PM PST; redundant check since there's a call to checkOptimProbValidity in evalLikelihoodAtPoint, which is the only place params could go crazy

  //double lambda_i = ploidy * currLibSizeScalingFactor * diploidDepth / 2;
  //double lambda_i = ploidy * currLibSizeScalingFactor * diploidDepth / 2 + diploidDepth / this->DEPTH_ERROR_SCALING; // add error term


  //double diploidPloidyProd = (*this->diploidDepthPloidyPreCalc)[windowIdx * (this->MAX_PLOIDY + 1) + ploidy]; // vec of ploidy * diploidDepth / 2
  //double lambda_i = currLibSizeScalingFactor * diploidPloidyProd + diploidDepth / this->DEPTH_ERROR_SCALING; // add error term
  //double lambda_i = currLibSizeScalingFactor * (*this->diploidDepthPloidyPreCalc)[windowIdx * (this->MAX_PLOIDY + 1) + ploidy] + diploidDepth / this->DEPTH_ERROR_SCALING; // add error term
  //double lambda_i = currLibSizeScalingFactor * ploidy * diploidDepth / 2.0 + diploidDepth / this->DEPTH_ERROR_SCALING; // old, calc on the fly


  // Fri 07 Feb 2020 09:50:37 PM PST trying out new X_ij ~ NB(r_ik + epsilon)
  //double lambda_i = ploidy * diploidDepth * currLibSizeScalingFactor + 1.0 / this->DEPTH_ERROR_SCALING;
  //double lambda_i = ploidy * diploidDepth / 2.0 * currLibSizeScalingFactor + diploidDepth / this->DEPTH_ERROR_SCALING;
  //double lambda_i = ploidy * diploidDepth * currLibSizeScalingFactor + diploidDepth / this->DEPTH_ERROR_SCALING;
  //double lambda_i = (ploidy * diploidDepth/2.0 + 1.0 / this->DEPTH_ERROR_SCALING) * currLibSizeScalingFactor;
  //double lambda_i = (ploidy * diploidDepth/2.0 + 65.37034 / this->DEPTH_ERROR_SCALING) * currLibSizeScalingFactor; // withErr
  //double lambda_i = (ploidy * diploidDepth/2.0 + 272.5568 / this->DEPTH_ERROR_SCALING) * currLibSizeScalingFactor; // withErr, unnormDiploid
  double lambda_i = (ploidy * diploidDepth/2.0) * currLibSizeScalingFactor + 272.5568 / this->DEPTH_ERROR_SCALING; // Thu 27 May 2021 10:50:36 AM PDT constant error
  //double lambda_i = (ploidy * diploidDepth/2.0 /*+ 65.37034 / this->DEPTH_ERROR_SCALING*/) * currLibSizeScalingFactor; // Thu 26 Mar 2020 01:55:01 PM PDT incStateSep debugging noErr
  /*// with no errors, if we're in ploidy == 0, then tumorDepths must be 0 // Thu 26 Mar 2020 01:55:01 PM PDT incStateSep debugging noErr
  if(ploidy == 0) {
    if(tumorDepth < 1e-4) {
      return 1;
    }
    else {
      return 0;
    }
  }*/




  //if(std::abs(lambda_i - old) > 1e-12) {
  //  std::cerr << "diploidPloidyProd: " << diploidPloidyProd << ", windowIdx: " << windowIdx << ", " << cellIdx << ";\tdiploidDepth: " << diploidDepth << ", ploidy: " << ploidy << ";\t" << "ploidy * dipDepth/2: " << ploidy * diploidDepth / 2.0;// << ", " << cellIdx << std::endl;
  //  fprintf(stderr, "lookup: %.40f,\tcalc: %.40f,\tdiff: %.40f\n", lambda_i, old, lambda_i - old);
  //  //std::cout << "\tprecalc: " << lambda_i << ", old: " << old << ", diff: " << lambda_i - old << std::endl;
  //}
  // tumorDepth ~ NegBinom(lambda_i, var=intercept + slope*lambda_i)
  //double var = this->getMeanVarianceIntercept() + this->getMeanVarianceSlope() * lambda_i;

  // tumorDepth ~ NegBinom(lambda_i, var=exp(intercept + slope * lambda_i))
  //double var = gsl_sf_exp(this->getMeanVarianceIntercept() + this->getMeanVarianceSlope() * lambda_i);

  // tumorDepth ~ NegBinom(lambda_i, var=intercept + slope * lambda_i + poly * lambda_i + poly2 * lambda_i * lambda_i
  double var = this->getMeanVarianceIntercept() + this->getMeanVarianceSlope() * lambda_i + this->getMeanVariancePoly2() * lambda_i * lambda_i;
  //double var = this->getMeanVarianceIntercept() + this->getMeanVarianceSlope() * lambda_i + this->getMeanVariancePoly2() * lambda_i * lambda_i + 144.3447; // Wed 15 Apr 2020 03:56:35 PM PDT debugging try adding in diploid variance, just hardcoded here
  //double var = this->getMeanVarianceIntercept() + this->getMeanVarianceSlope() * lambda_i + this->getMeanVariancePoly2() * lambda_i * lambda_i + 144.3447 / 2; // Wed 15 Apr 2020 03:56:35 PM PDT debugging try adding in diploid variance, just hardcoded here

  // tumorDepth ~ NegBinom(lambda_i, var=intercept + slope * lambda_i + poly2 * lambda_i * lambda_i + poly2 * lambda_i * lambda_i + poly3 * lambda_i^3
  //double var = this->getMeanVarianceIntercept() + this->getMeanVarianceSlope() * lambda_i + this->getMeanVariancePoly2() * lambda_i * lambda_i + gsl_vector_get(this->meanVarianceCoefVec, 4) * lambda_i * lambda_i * lambda_i;

  // in the boost library, p = P(success), r = # successes
  double p = lambda_i / var;
  double r = lambda_i * lambda_i / (var - lambda_i);

  //std::cerr << "currLibSizeScalingFactor: " << currLibSizeScalingFactor << ", lambda_i: " << lambda_i << ", var: " << var << ", p: " << p << ", r: " << r << std::endl;
  // if diploidDepth is so low that p or r goes negative, assume this is a hard to map region
  if((p < 0 || r < 0) && diploidDepth < 1) {
    //std::cerr << "P OR R TOO SMALL AND DIPLOIDDEPTH < 1; RETURNING 1" << std::endl;
    return 1;
  }

  // Mon 30 Mar 2020 03:55:48 PM PDT try matching the simulation fix for small r
  //std::cerr << r << std::endl;
  if(r < 1) {
    //std::cerr << "R TOO SMALL; ADJUSTING" << std::endl;
    r = 1;
    var = lambda_i * lambda_i + lambda_i;
    p = lambda_i / var;
  }

  double likelihood = -1;
  int k = (int) tumorDepth;
  // if underdispersed, push back to poisson
  /*
  if(var < lambda_i) {
    //std::cout << "POISSON lambda_i: " << lambda_i;
    likelihood = pdf(boost::math::poisson(lambda_i), k);
    //std::cout << ", likelihood: " << likelihood << std::endl;
  }
  // otherwise, use neg binom
  else {*/
    //double boostMean = mean(boost::math::negative_binomial(r,p)); // check r, p calcs are correct
    //double boostVar = variance(boost::math::negative_binomial(r,p));
    //std::cout << "mean: " << boostMean << ", var: " << boostVar << std::endl;
    //std::cout << "p: " << p << ", r: " << r << ", tumorDepth: " << tumorDepth  << ", diploidDepth: " << diploidDepth << ", ploidy: " << ploidy << ", scalingFactor: " << this->getLibScalingFactor(cellIdx);
    //likelihood = pdf(boost::math::negative_binomial(r, p), k); // this is the boost library version
    //std::cout << "r: " << r << ", p: " << p << ", k: " << k << std::endl;
    //std::cout << "boost: " << pdf(boost::math::negative_binomial(r, p), k) << std::endl; // this is the boost library version

    //double coef = exp(lgamma(r + k) - lgamma(r) - lgamma(k+1));
    //likelihood = exp(lgamma(r+k) - lgamma(r) - this->getLogFacK(k) + r * log(p) + k * log(1-p)); // extra exp/log
    likelihood = (lgamma(r+k) - lgamma(r) - this->getLogFacK(k) + r * log(p) + k * log(1-p));
    //std::cerr << ", likelihood: " << likelihood << std::endl;

    //double manualLikelihood = exp(lgamma(r + k) - lgamma(r) - lgamma(k+1)) * pow(p, r) * pow(1-p, k);
    //std::cout << "manualLikelihood: " << manualLikelihood << std::endl;

    //std::cout << "log coef: " << lgamma(r+k) - lgamma(r) - this->getLogFacK(k) << std::endl;
    //std::cout << "log gam r+k: " << lgamma(r+k) << std::endl;
    //std::cout << "log gam r: " << lgamma(r) << std::endl;
    //std::cout << "log k: " << this->getLogFacK(k) << std::endl;
    //std::cout << "log probs: " << r * log(p) + k * log(1-p) << std::endl;
    //std::cout << "sum: " << lgamma(r+k) - lgamma(r) - this->getLogFacK(k) + r * log(p) + k * log(1-p) << std::endl;
  //}
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
    //emissionProb += log(this->getEmissionProb(currTumorDepth, currDiploidDepth, currPloidy, cellIdx)); // extra exp/log
    emissionProb += this->getLogEmissionProb(currTumorDepth, currDiploidDepth, currPloidy, cellIdx);
  }
  return emissionProb;
  //return exp(emissionProb);
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
  //int maxDepth = *(this->alphabet->rbegin()); // sets are stored in sorted order
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
  // create new HMM with the best guess parameters and return it
  //TwoCell3TrParam2DegPolyHMM* bestGuessHMM = new TwoCell3TrParam2DegPolyHMM(*this);
  TwoCell3TrParam2DegPolyHMM* bestGuessHMM = this;//new TwoCell3TrParam2DegPolyHMM(*this);
  //std::cerr << "HERE" << std::endl;
  //bestGuessHMM->print(stdout);
  //std::cerr << "HERE2" << std::endl;
  gsl_vector* initGuessAsParams = gsl_vector_alloc(initGuess->size);
  //bestGuessHMM->setParamsToEst(initGuess);
  //this->convertProbToParam(initGuessAsParams, initGuess, this->getKploidy());
  this->convertProbToParam(initGuessAsParams, initGuess);
  //HMM::bfgs(initGuess, bestGuessHMM, verbose);
  //HMM::bfgs(initGuessAsParams, bestGuessHMM, verbose);
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
///*
// * function to simulate states and read depths. If(simDiploid) then diploid
// * depths are simulated. If not, then only tumor depths are simulated. This is
// * useful for multiple tumors simulated from the same parameters, since the average
// * diploid should not change between simulations. There no checks that diploid info has
// * been previously simulated
// *
// */
//void TwoCell3TrParam2DegPolyHMM::simulate(int seed, bool simDiploid, double diploid_lambda_i, int numDiploidCells) {
//  if(this->generator == nullptr) {
//    this->generator = new base_generator_type(seed);
//  }
//  if(this->stateRNG == nullptr) {
//    this->stateRNG = new base_generator_type(seed * 333);
//    //this->stateRNG = new base_generator_type(seed << 3);
//  }
//  if(this->diploidRNG == nullptr) {
//    this->diploidRNG = new base_generator_type(seed * 555555);
//    //this->diploidRNG = new base_generator_type(seed << 5);
//  }
//  if(this->tumor0RNG == nullptr) {
//    this->tumor0RNG = new base_generator_type(seed * 7777777);
//    //this->tumor0RNG = new base_generator_type(seed << 7);
//  }
//  if(this->tumor1RNG == nullptr) {
//    this->tumor1RNG = new base_generator_type(seed * 11111111111);
//    //this->tumor1RNG = new base_generator_type(seed << 11);
//  }
//  //std::cout << this->generator << std::endl;
//  // set up random number generator
//  // from https://www.boost.org/doc/libs/1_70_0/libs/random/example/random_demo.cpp
//  //base_generator_type generator(seed);
//  boost::uniform_real<> uni_dist(0,1);
//  //boost::variate_generator<base_generator_type&, boost::uniform_real<>> uni(*this->generator, uni_dist);
//  boost::variate_generator<base_generator_type&, boost::uniform_real<>> uni(*this->stateRNG, uni_dist);
//
//  // simulate states
//  std::vector<std::string>* diploidStates = nullptr;
//  std::vector<std::string>* tumor0States = nullptr;
//  std::vector<std::string>* tumor1States = nullptr;
//  //int numStateTransitions = 0;
//  //gsl_matrix* numTransitionsMat = gsl_matrix_alloc(this->transition->size1, this->transition->size2);
//  gsl_matrix_set_zero(this->numTransitionsMat);
//
//  // clear out DepthPair variables. This is important if multiple simulations are run
//  for(unsigned int i = 0; i < this->depths->size(); i++) {
//    if(simDiploid) {
//      (*this->depths)[i]->diploidLibrarySize = 0;
//      (*this->depths)[i]->maxDiploidDepth = 0;
//    }
//    (*this->depths)[i]->tumorLibrarySize = 0;
//    (*this->depths)[i]->maxTumorDepth = 0;
//  }
//
//  std::string currChr;
//  int currNumWindows = -1;
//  std::vector<std::string>* chrVec = this->getChrVec();
//  int fromStateIdx = -1;
//  int toStateIdx = -1;
//  // for each chr
//  for(unsigned int i = 0; i < chrVec->size(); i++) {
//    // init this chr's vectors
//    diploidStates = new std::vector<std::string>();
//    tumor0States = new std::vector<std::string>();
//    tumor1States = new std::vector<std::string>();
//
//    currChr = (*chrVec)[i];
//    currNumWindows = (*(*this->depths)[0]->regions)[currChr]->size();
//
//    // diploid always in state diploid
//    diploidStates->push_back("2");
//
//    // start tumor in steady state
//    double stateProb = uni();
//    //std::cout << stateProb << std::endl;
//    fromStateIdx = getRandStateIdx(stateProb, this->initProb);
//
//    // get tumor cell specific states
//    tumor0States->push_back(std::to_string(getCellPloidyFromStateIdx(0, fromStateIdx)));
//    tumor1States->push_back(std::to_string(getCellPloidyFromStateIdx(1, fromStateIdx)));
//
//    // for each window in this chr
//    for(int winIdx = 1; winIdx < currNumWindows; winIdx++) {
//      // diploid always in state diploid
//      diploidStates->push_back("2");
//
//      // tumor state changes according to transition matrix
//      stateProb = uni();
//      //std::cout << stateProb << std::endl;
//      toStateIdx = getRandStateIdx(stateProb, fromStateIdx);
//
//      // if toStateIdx is different from where came from, then must have transitioned states
//      //if(toStateIdx != fromStateIdx) {
//      //  //numStateTransitions++;
//      //  // add one to count of these transitions
//      //}
//      gsl_matrix_set(this->numTransitionsMat, fromStateIdx, toStateIdx, gsl_matrix_get(this->numTransitionsMat, fromStateIdx, toStateIdx) + 1);
//      fromStateIdx = toStateIdx;
//
//      // get tumor cell specific states
//      tumor0States->push_back(std::to_string(getCellPloidyFromStateIdx(0, fromStateIdx)));
//      tumor1States->push_back(std::to_string(getCellPloidyFromStateIdx(1, fromStateIdx)));
//    }
//
//    // save
//    (*(*this->depths)[0]->chrToDiploidSimStateMap)[currChr] = diploidStates;
//    (*(*this->depths)[0]->chrToTumorSimStateMap)[currChr] = tumor0States;
//    (*(*this->depths)[1]->chrToTumorSimStateMap)[currChr] = tumor1States;
//  }
//
//  // values from means of input/diploidBinsVarMeanNormLibSize_250kb.coordMeanVar
//  //double diploid_lambda_i = 65.37034;
//  //double diploid_var = 144.3447;
//  double diploid_var = this->getMeanVarianceIntercept() + this->getMeanVarianceSlope() * diploid_lambda_i + this->getMeanVariancePoly2() * diploid_lambda_i * diploid_lambda_i; // Mon 20 Apr 2020 08:10:12 PM PDT now all diploid variance is calculated using the polynomial mean/var relationship
//  double diploid_p = diploid_lambda_i / diploid_var;
//  double diploid_r = diploid_lambda_i * diploid_lambda_i / (diploid_var - diploid_lambda_i);
//  //std::cout << diploid_r << ", " << diploid_p  << ", " << std::endl;
//  boost::random::negative_binomial_distribution<> diploid_negBinom_dist(diploid_r, diploid_p);
//  //boost::variate_generator<base_generator_type&, boost::random::negative_binomial_distribution<>> diploid_negBinom(*this->generator, diploid_negBinom_dist);
//  //boost::variate_generator<base_generator_type&, boost::random::negative_binomial_distribution<>> diploid_negBinom(*dipGen, diploid_negBinom_dist);
//  boost::variate_generator<base_generator_type&, boost::random::negative_binomial_distribution<>> diploid_negBinom(*this->diploidRNG, diploid_negBinom_dist);
//  //boost::random::poisson_distribution<> diploid_pois_dist(diploid_lambda_i); // set1 Mon 06 Apr 2020 05:28:58 PM PDT debugging poisDiploid
//  //boost::random::poisson_distribution<> diploid_pois_dist(diploid_lambda_i * 10); // set2
//  //boost::random::poisson_distribution<> diploid_pois_dist(diploid_lambda_i * 100); // set3
//  //boost::variate_generator<base_generator_type&, boost::random::poisson_distribution<>> diploid_pois(*this->diploidRNG, diploid_pois_dist);
//  //boost::random::normal_distribution<> diploid_norm_dist(diploid_lambda_i, sqrt(diploid_var)); // set1 Thu 09 Apr 2020 07:00:50 PM PDT debugging normDiploid; params are mean and sd
//  //boost::variate_generator<base_generator_type&, boost::random::normal_distribution<>> diploid_norm(*this->diploidRNG, diploid_norm_dist);
//
//  // given states, simulate coverage according to neg binom model
//  double tumor_lambda_i = -1;
//  double tumor_var = -1;
//  double tumor_p = -1;
//  double tumor_r = -1;
//  int currTumor0Ploidy = -1;
//  int currTumor1Ploidy = -1;
//  double maxSim = -1;
//  double currDiploidSim = -1;
//  double currTumorSim = -1;
//  double currLibSizeScalingFactor = -1;
//  // for each chr
//  for(unsigned int j = 0; j < chrVec->size(); j++) {
//    currChr = (*chrVec)[j];
//    diploidStates = (*(*this->depths)[0]->chrToDiploidSimStateMap)[currChr];
//    tumor0States = (*(*this->depths)[0]->chrToTumorSimStateMap)[currChr];
//    tumor1States = (*(*this->depths)[1]->chrToTumorSimStateMap)[currChr];
//    // for each simulated state
//    for(unsigned int i = 0; i < diploidStates->size(); i++) {
//      // simulate diploid
//      if(simDiploid) {
//        //currDiploidSim = diploid_negBinom();
//        //currDiploidSim = diploid_pois();
//        //currDiploidSim = diploid_norm();
//        //if(currDiploidSim < 0) {
//        //  currDiploidSim = 0;
//        //}
//        //currDiploidSim = diploid_lambda_i;
//        // Fri 17 Apr 2020 10:42:38 PM PDT debugging simulating multiple diploid cells, then averaging
//        //int numDiploidCells = 35;
//        double dipSimVec[numDiploidCells];
//        for(int dipIdx = 0; dipIdx < numDiploidCells; dipIdx++) {
//          dipSimVec[dipIdx] = diploid_negBinom();
//          std::cerr << dipSimVec[dipIdx] << std::endl;
//        }
//        currDiploidSim = gsl_stats_mean(dipSimVec, 1, numDiploidCells);
//        //std::cout << currDiploidSim << std::endl;
//        (*(*(*this->depths)[0]->chrToDiploidDepthMap)[currChr])[i] = currDiploidSim;
//        if(currDiploidSim > maxSim) {
//          maxSim = currDiploidSim;
//          (*this->depths)[0]->maxDiploidDepth = currDiploidSim;
//        }
//        (*this->depths)[0]->diploidLibrarySize += currDiploidSim;
//        (*this->depths)[1]->diploidLibrarySize += currDiploidSim;
//        //(*(*(*this->depths)[0]->chrToDiploidVarMap)[currChr])[i] = diploid_var;
//        (*(*(*this->depths)[0]->chrToDiploidVarMap)[currChr])[i] = gsl_stats_variance(dipSimVec, 1, numDiploidCells); // Fri 17 Apr 2020 10:42:38 PM PDT debugging simulating multiple diploid cells, then averaging
//        //diploid_lambda_i = currDiploidSim;
//      }
//      /*else {
//        diploid_lambda_i = (*(*(*this->depths)[0]->chrToDiploidDepthMap)[currChr])[i];
//      }*/ // Thu 20 Feb 2020 02:41:19 PM PST use true diploid_lambda_i while simulating tumors
//
//      // parse tumor states
//      currTumor0Ploidy = atoi((*tumor0States)[i].c_str());
//      currTumor1Ploidy = atoi((*tumor1States)[i].c_str());
//
//      // simulate tumor 0
//      currLibSizeScalingFactor = this->getLibScalingFactor(0);
//      //tumor_lambda_i = currTumor0Ploidy * currLibSizeScalingFactor * diploid_lambda_i / 2 + diploid_lambda_i / this->DEPTH_ERROR_SCALING; // add error term
//      tumor_lambda_i = currLibSizeScalingFactor * (currTumor0Ploidy * diploid_lambda_i / 2 + diploid_lambda_i / this->DEPTH_ERROR_SCALING); // add error term withErr
//      //tumor_lambda_i = currLibSizeScalingFactor * (currTumor0Ploidy * diploid_lambda_i / 2 /*+ diploid_lambda_i / this->DEPTH_ERROR_SCALING*/); // add error term  // Thu 26 Mar 2020 01:55:01 PM PDT incStateSep debugging noErr
//      tumor_var = this->getMeanVarianceIntercept() + this->getMeanVarianceSlope() * tumor_lambda_i + this->getMeanVariancePoly2() * tumor_lambda_i * tumor_lambda_i;
//      tumor_p = tumor_lambda_i / tumor_var;
//      tumor_r = tumor_lambda_i * tumor_lambda_i / (tumor_var - tumor_lambda_i);
//      //std::cerr << "currLibSizeScalingFactor: " << currLibSizeScalingFactor << ", tumor_lambda_i: " << tumor_lambda_i << ", tumor_var: " << tumor_var << ", tumor_p: " << tumor_p << ", tumor_r: " << tumor_r << std::endl;
//      //std::cerr << currTumor0Ploidy << ", " << currLibSizeScalingFactor << ", " << diploid_lambda_i << ", " << diploid_lambda_i / this->DEPTH_ERROR_SCALING << ", const:" << this->DEPTH_ERROR_SCALING << ", " << tumor_lambda_i << ", " << tumor_var << ", " << tumor_r << ", " << tumor_p << std::endl;
//      //std::cerr << currTumor0Ploidy << ", " << diploid_lambda_i << std::endl;
//      if(tumor_r < 1) {
//        //likelihood = pdf(boost::math::poisson(lambda_i), (int)tumorDepth);
//
//        /*
//        boost::random::poisson_distribution<> tumor0_pois_dist(diploid_lambda_i / this->DEPTH_ERROR_SCALING);
//        boost::variate_generator<base_generator_type&, boost::random::poisson_distribution<>> tumor0_pois(generator, tumor0_pois_dist);
//        currTumorSim = tumor0_pois();
//        */
//
//        // fix r=1, adj var, recalc p
//        // r = lambda^2 / (var - lambda); v = lambda^2 / r + lambda
//        tumor_r = 1;
//        tumor_var = tumor_lambda_i * tumor_lambda_i + tumor_lambda_i;
//        tumor_p = tumor_lambda_i / tumor_var;
//
//        //currTumorSim = 0; // TODO can incorporate more complicated error model here. currently, just emit 0 if ploidy is 0
//        //std::cerr << tumor_r << ", ";
//        //tumor_r = -(tumor_r-1);
//        //tumor_r = tumor_r + 1;
//        //std::cerr << tumor_r << std::endl;
//      }
//
//
//      /*// Thu 26 Mar 2020 01:55:01 PM PDT incStateSep debugging noErr
//      if(currTumor0Ploidy == 0) {
//        currTumorSim = 0;
//      } else {*/
//
//      //else {
//        boost::random::negative_binomial_distribution<> tumor0_negBinom_dist(tumor_r,tumor_p);
//        //boost::variate_generator<base_generator_type&, boost::random::negative_binomial_distribution<>> tumor0_negBinom(*this->generator, tumor0_negBinom_dist);
//        boost::variate_generator<base_generator_type&, boost::random::negative_binomial_distribution<>> tumor0_negBinom(*this->tumor0RNG, tumor0_negBinom_dist);
//        currTumorSim = tumor0_negBinom();
//        //currTumorSim = tumor_lambda_i;
//        //std::cout << currTumorSim << std::endl;
//        //std::cout << tumor0_negBinom_dist.p() << ", " << tumor0_negBinom_dist.k() << std::endl;
//        //exit(0);
//        //for(int tmp = 0; tmp < 10000; tmp++) {
//        //  fprintf(stderr, "%i\n", tumor0_negBinom());
//        //}
//        //exit(0);
//      //}
//      //} // Thu 26 Mar 2020 01:55:01 PM PDT incStateSep debugging noErr
//
//      // save simulated values, as well as max for alphabet
//      (*(*(*this->depths)[0]->chrToTumorDepthMap)[currChr])[i] = currTumorSim;
//      if(currTumorSim > maxSim) {
//        maxSim = currTumorSim;
//        (*this->depths)[0]->maxTumorDepth = currTumorSim;
//      }
//      (*this->depths)[0]->tumorLibrarySize += currTumorSim;
//
//      // simulate tumor 1
//      currLibSizeScalingFactor = this->getLibScalingFactor(1);
//      //tumor_lambda_i = currTumor1Ploidy * currLibSizeScalingFactor * diploid_lambda_i / 2 + diploid_lambda_i / this->DEPTH_ERROR_SCALING; // add error term
//      tumor_lambda_i = currLibSizeScalingFactor * (currTumor1Ploidy * diploid_lambda_i / 2 + diploid_lambda_i / this->DEPTH_ERROR_SCALING); // add error term withErr
//      //tumor_lambda_i = currLibSizeScalingFactor * (currTumor1Ploidy * diploid_lambda_i / 2 /*+ diploid_lambda_i / this->DEPTH_ERROR_SCALING*/); // add error term  // Thu 26 Mar 2020 01:55:01 PM PDT incStateSep debugging noErr
//      tumor_var = this->getMeanVarianceIntercept() + this->getMeanVarianceSlope() * tumor_lambda_i + this->getMeanVariancePoly2() * tumor_lambda_i * tumor_lambda_i;
//      tumor_p = tumor_lambda_i / tumor_var;
//      tumor_r = tumor_lambda_i * tumor_lambda_i / (tumor_var - tumor_lambda_i);
//      //std::cerr << currTumorPloidy << ", " << diploid_lambda_i << ", " << diploid_lambda_i / this->DEPTH_ERROR_SCALING << ", const:" << this->DEPTH_ERROR_SCALING << ", " << tumor_lambda_i << ", " << tumor_var << ", " << tumor_r << ", " << tumor_p << std::endl;
//      //std::cout << "before: " << tumor_p << ", " << tumor_r << std::endl;
//      if(tumor_r < 1) {
//        //currTumorSim = 0; // TODO can incorporate more complicated error model here. currently, just emit 0 if ploidy is 0
//        //tumor_r = -(tumor_r-1);
//        //tumor_r = tumor_r + 1;
//        /*
//        boost::random::poisson_distribution<> tumor1_pois_dist(diploid_lambda_i / this->DEPTH_ERROR_SCALING);
//        boost::variate_generator<base_generator_type&, boost::random::poisson_distribution<>> tumor1_pois(generator, tumor1_pois_dist);
//        currTumorSim = tumor1_pois();
//        */
//        tumor_r = 1;
//        tumor_var = tumor_lambda_i * tumor_lambda_i + tumor_lambda_i;
//        tumor_p = tumor_lambda_i / tumor_var;
//        //std::cout << "after: " << tumor_p << ", " << tumor_r << std::endl;
//      }
//
//
//      /*// Thu 26 Mar 2020 01:55:01 PM PDT incStateSep debugging noErr
//      if(currTumor1Ploidy == 0) {
//        currTumorSim = 0;
//      } else {*/
//
//      //else {
//        boost::random::negative_binomial_distribution<> tumor1_negBinom_dist(tumor_r,tumor_p);
//        //boost::variate_generator<base_generator_type&, boost::random::negative_binomial_distribution<>> tumor1_negBinom(*this->generator, tumor1_negBinom_dist);
//        boost::variate_generator<base_generator_type&, boost::random::negative_binomial_distribution<>> tumor1_negBinom(*this->tumor1RNG, tumor1_negBinom_dist);
//        currTumorSim = tumor1_negBinom();
//        //currTumorSim = tumor_lambda_i;
//        //std::cout << currTumorSim << std::endl;
//      //}
//      //} // Thu 26 Mar 2020 01:55:01 PM PDT incStateSep debugging noErr
//
//      // save simulated values, as well as max for alphabet
//      (*(*(*this->depths)[1]->chrToTumorDepthMap)[currChr])[i] = currTumorSim;
//      if(currTumorSim > maxSim) {
//        maxSim = currTumorSim;
//        (*this->depths)[1]->maxTumorDepth = currTumorSim;
//      }
//      (*this->depths)[1]->tumorLibrarySize += currTumorSim;
//
//    }
//  }
//
//  // save alphabet
//  std::vector<int> fullAlphabet((int)maxSim + 1); // add one to make inclusive
//  std::iota(fullAlphabet.begin(), fullAlphabet.end(), 0);
//  this->setAlphabet(new std::set<int>(fullAlphabet.begin(), fullAlphabet.end()));
//
//  // save this->logFacKVec
//  this->setLogFacK();
//
//  //// precalc the diploidDepth * ploidy / 2 vector for getEmissionProb
//  //this->fillDiploidDepthPloidyPreCalc();
//
//  //return numStateTransitions;
//  //return this->numTransitionsMat;
//
//  // set library sizes according to final tumor size
//  this->estLibScalingFactorsPosterior();
//}

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
  //std::cout << this->generator << std::endl;
  // set up random number generator
  // from https://www.boost.org/doc/libs/1_70_0/libs/random/example/random_demo.cpp
  //base_generator_type generator(seed);
  boost::uniform_real<> uni_dist(0,1);
  boost::variate_generator<base_generator_type&, boost::uniform_real<>> uni(*this->generator, uni_dist);
  //boost::variate_generator<base_generator_type&, boost::uniform_real<>> uni(*this->stateRNG, uni_dist);

  // simulate states
  std::vector<std::string>* diploidStates = nullptr;
  std::vector<std::string>* tumor0States = nullptr;
  std::vector<std::string>* tumor1States = nullptr;
  //int numStateTransitions = 0;
  //gsl_matrix* numTransitionsMat = gsl_matrix_alloc(this->transition->size1, this->transition->size2);
  //gsl_matrix_set_zero(this->numTransitionsMat);

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
      //if(toStateIdx != fromStateIdx) {
      //  //numStateTransitions++;
      //  // add one to count of these transitions
      //}
      //gsl_matrix_set(this->numTransitionsMat, fromStateIdx, toStateIdx, gsl_matrix_get(this->numTransitionsMat, fromStateIdx, toStateIdx) + 1);
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
  //double diploid_lambda_i = 65.37034;
  //double diploid_var = 144.3447;
  double diploid_var = this->getMeanVarianceIntercept() + this->getMeanVarianceSlope() * diploid_lambda_i + this->getMeanVariancePoly2() * diploid_lambda_i * diploid_lambda_i; // Mon 20 Apr 2020 08:10:12 PM PDT now all diploid variance is calculated using the polynomial mean/var relationship
  double diploid_p = diploid_lambda_i / diploid_var;
  double diploid_r = diploid_lambda_i * diploid_lambda_i / (diploid_var - diploid_lambda_i);
  //std::cout << diploid_r << ", " << diploid_p  << ", " << std::endl;
  boost::random::negative_binomial_distribution<> diploid_negBinom_dist(diploid_r, diploid_p);
  boost::variate_generator<base_generator_type&, boost::random::negative_binomial_distribution<>> diploid_negBinom(*this->generator, diploid_negBinom_dist);
  //boost::variate_generator<base_generator_type&, boost::random::negative_binomial_distribution<>> diploid_negBinom(*this->diploidRNG, diploid_negBinom_dist);

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
        //(*(*(*this->depths)[0]->chrToDiploidVarMap)[currChr])[i] = diploid_var;
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
      //tumor_lambda_i = currTumor0Ploidy * currLibSizeScalingFactor * diploid_lambda_i / 2 + diploid_lambda_i / this->DEPTH_ERROR_SCALING; // add error term
      tumor_lambda_i = currLibSizeScalingFactor * (currTumor0Ploidy * diploid_lambda_i / 2 + diploid_lambda_i / this->DEPTH_ERROR_SCALING); // add error term withErr
      //tumor_lambda_i = currLibSizeScalingFactor * (currTumor0Ploidy * diploid_lambda_i / 2 /*+ diploid_lambda_i / this->DEPTH_ERROR_SCALING*/); // add error term  // Thu 26 Mar 2020 01:55:01 PM PDT incStateSep debugging noErr
      tumor_var = this->getMeanVarianceIntercept() + this->getMeanVarianceSlope() * tumor_lambda_i + this->getMeanVariancePoly2() * tumor_lambda_i * tumor_lambda_i;
      tumor_p = tumor_lambda_i / tumor_var;
      tumor_r = tumor_lambda_i * tumor_lambda_i / (tumor_var - tumor_lambda_i);
      //std::cerr << "currLibSizeScalingFactor: " << currLibSizeScalingFactor << ", tumor_lambda_i: " << tumor_lambda_i << ", tumor_var: " << tumor_var << ", tumor_p: " << tumor_p << ", tumor_r: " << tumor_r << std::endl;
      //std::cerr << currTumor0Ploidy << ", " << currLibSizeScalingFactor << ", " << diploid_lambda_i << ", " << diploid_lambda_i / this->DEPTH_ERROR_SCALING << ", const:" << this->DEPTH_ERROR_SCALING << ", " << tumor_lambda_i << ", " << tumor_var << ", " << tumor_r << ", " << tumor_p << std::endl;
      //std::cerr << currTumor0Ploidy << ", " << diploid_lambda_i << std::endl;
      if(tumor_r < 1) {
        // fix r=1, adj var, recalc p
        // r = lambda^2 / (var - lambda); v = lambda^2 / r + lambda
        tumor_r = 1;
        tumor_var = tumor_lambda_i * tumor_lambda_i + tumor_lambda_i;
        tumor_p = tumor_lambda_i / tumor_var;

        //currTumorSim = 0; // TODO can incorporate more complicated error model here. currently, just emit 0 if ploidy is 0
        //std::cerr << tumor_r << ", ";
        //tumor_r = -(tumor_r-1);
        //tumor_r = tumor_r + 1;
        //std::cerr << tumor_r << std::endl;
      }

      /*// Thu 26 Mar 2020 01:55:01 PM PDT incStateSep debugging noErr
      if(currTumor0Ploidy == 0) {
        currTumorSim = 0;
      } else {*/

      //else {
        boost::random::negative_binomial_distribution<> tumor0_negBinom_dist(tumor_r,tumor_p);
        boost::variate_generator<base_generator_type&, boost::random::negative_binomial_distribution<>> tumor0_negBinom(*this->generator, tumor0_negBinom_dist);
        //boost::variate_generator<base_generator_type&, boost::random::negative_binomial_distribution<>> tumor0_negBinom(*this->tumor0RNG, tumor0_negBinom_dist);
        currTumorSim = tumor0_negBinom();
        //currTumorSim = tumor_lambda_i;
        //std::cout << currTumorSim << std::endl;
        //std::cout << tumor0_negBinom_dist.p() << ", " << tumor0_negBinom_dist.k() << std::endl;
        //exit(0);
        //for(int tmp = 0; tmp < 10000; tmp++) {
        //  fprintf(stderr, "%i\n", tumor0_negBinom());
        //}
        //exit(0);
      //}
      //} // Thu 26 Mar 2020 01:55:01 PM PDT incStateSep debugging noErr

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
      //tumor_lambda_i = currTumor1Ploidy * currLibSizeScalingFactor * diploid_lambda_i / 2 + diploid_lambda_i / this->DEPTH_ERROR_SCALING; // add error term
      tumor_lambda_i = currLibSizeScalingFactor * (currTumor1Ploidy * diploid_lambda_i / 2 + diploid_lambda_i / this->DEPTH_ERROR_SCALING); // add error term withErr
      //tumor_lambda_i = currLibSizeScalingFactor * (currTumor1Ploidy * diploid_lambda_i / 2 /*+ diploid_lambda_i / this->DEPTH_ERROR_SCALING*/); // add error term  // Thu 26 Mar 2020 01:55:01 PM PDT incStateSep debugging noErr
      tumor_var = this->getMeanVarianceIntercept() + this->getMeanVarianceSlope() * tumor_lambda_i + this->getMeanVariancePoly2() * tumor_lambda_i * tumor_lambda_i;
      tumor_p = tumor_lambda_i / tumor_var;
      tumor_r = tumor_lambda_i * tumor_lambda_i / (tumor_var - tumor_lambda_i);
      //std::cerr << currTumorPloidy << ", " << diploid_lambda_i << ", " << diploid_lambda_i / this->DEPTH_ERROR_SCALING << ", const:" << this->DEPTH_ERROR_SCALING << ", " << tumor_lambda_i << ", " << tumor_var << ", " << tumor_r << ", " << tumor_p << std::endl;
      //std::cout << "before: " << tumor_p << ", " << tumor_r << std::endl;
      if(tumor_r < 1) {
        //currTumorSim = 0; // TODO can incorporate more complicated error model here. currently, just emit 0 if ploidy is 0
        //tumor_r = -(tumor_r-1);
        //tumor_r = tumor_r + 1;
        tumor_r = 1;
        tumor_var = tumor_lambda_i * tumor_lambda_i + tumor_lambda_i;
        tumor_p = tumor_lambda_i / tumor_var;
        //std::cout << "after: " << tumor_p << ", " << tumor_r << std::endl;
      }

      /*// Thu 26 Mar 2020 01:55:01 PM PDT incStateSep debugging noErr
      if(currTumor1Ploidy == 0) {
        currTumorSim = 0;
      } else {*/

      //else {
        boost::random::negative_binomial_distribution<> tumor1_negBinom_dist(tumor_r,tumor_p);
        boost::variate_generator<base_generator_type&, boost::random::negative_binomial_distribution<>> tumor1_negBinom(*this->generator, tumor1_negBinom_dist);
        //boost::variate_generator<base_generator_type&, boost::random::negative_binomial_distribution<>> tumor1_negBinom(*this->tumor1RNG, tumor1_negBinom_dist);
        currTumorSim = tumor1_negBinom();
        //currTumorSim = tumor_lambda_i;
        //std::cout << currTumorSim << std::endl;
      //}
      //} // Thu 26 Mar 2020 01:55:01 PM PDT incStateSep debugging noErr

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
  //std::vector<int> fullAlphabet((int)maxSim + 1); // add one to make inclusive
  //std::iota(fullAlphabet.begin(), fullAlphabet.end(), 0);
  //this->setAlphabet(new std::set<int>(fullAlphabet.begin(), fullAlphabet.end()));
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
//void TwoCell3TrParam2DegPolyHMM::setInitGuessNthTime(gsl_vector* initGuess, int iter, int n, double max) const {
void TwoCell3TrParam2DegPolyHMM::setInitGuessNthTime(gsl_vector* initGuess, int iter, int numTotalRuns) const {
  /*// fill initGuessCases if this is the first iter
  if(iter == 0) {
    //this->initGuessCases->clear();
    for(int i = 0; i < this->maxNumBFGSStarts; i++) {
      this->initGuessCases->push_back(i);
    }
    std::shuffle(this->initGuessCases->begin(), this->initGuessCases->end(), *(this->generator));
  }
  int currInitGuessCase = (*this->initGuessCases)[iter];*/

  /*boost::random::uniform_int_distribution<> initGuess_unif_dist(0,18);
  boost::variate_generator<base_generator_type&, boost::random::uniform_int_distribution<>> initGuess_unif(*this->generator, initGuess_unif_dist);
  int currInitGuessCase = -1;
  // pick cases as long as they haven't been picked before
  while(true) {
    currInitGuessCase = initGuess_unif();
    std::cout << "&&&&& " << currInitGuessCase << std::endl;

    // if haven't seen before, pick it
    if(this->initGuessCases->find(currInitGuessCase) == this->initGuessCases->end()) {
      this->initGuessCases->insert(currInitGuessCase);
      break;
    }
  }*/

  //std::cout << "PARAM SET: " << currInitGuessCase << std::endl;
  //switch(currInitGuessCase) {
  switch(iter) {
    // beta/gamma all at 0.05, t1/t2/t3 at 0.1
    case 0:
      for(int cellIdx = 0; cellIdx < this->NUM_LIBS_TO_EST; cellIdx++) {
        gsl_vector_set(initGuess, this->LIB_SIZE_SCALING_FACTOR_START_IDX + cellIdx, 1);
      }
      //gsl_vector_set(initGuess, this->LIB_SIZE_SCALING_FACTOR_START_IDX + 0, .75);
      //gsl_vector_set(initGuess, this->LIB_SIZE_SCALING_FACTOR_START_IDX + 1, .75);
      //gsl_vector_set(initGuess, this->LIB_SIZE_SCALING_FACTOR_START_IDX + 0, this->getLibScalingFactor(0));
      //gsl_vector_set(initGuess, this->LIB_SIZE_SCALING_FACTOR_START_IDX + 1, this->getLibScalingFactor(1));
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
      //gsl_vector_set(initGuess, this->LIB_SIZE_SCALING_FACTOR_START_IDX + 0, 1);
      //gsl_vector_set(initGuess, this->LIB_SIZE_SCALING_FACTOR_START_IDX + 1, 1);
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
      //gsl_vector_set(initGuess, this->LIB_SIZE_SCALING_FACTOR_START_IDX + 0, 1);
      //gsl_vector_set(initGuess, this->LIB_SIZE_SCALING_FACTOR_START_IDX + 1, 1);
      gsl_vector_set(initGuess, this->TRANSITION_PROB_START_IDX + 0, 0.02);
      gsl_vector_set(initGuess, this->TRANSITION_PROB_START_IDX + 1, 0.05);
      gsl_vector_set(initGuess, this->TRANSITION_PROB_START_IDX + 2, 0.02);
      gsl_vector_set(initGuess, this->TRANSITION_PROB_START_IDX + 3, 0.02);
      gsl_vector_set(initGuess, this->TRANSITION_PROB_START_IDX + 4, 0.02);
      break;

    /*// start lib sizes at right place, start branches higher up
    case 3:
      for(int cellIdx = 0; cellIdx < this->NUM_LIBS_TO_EST; cellIdx++) {
        gsl_vector_set(initGuess, this->LIB_SIZE_SCALING_FACTOR_START_IDX + cellIdx, 1);
      }
      //gsl_vector_set(initGuess, this->LIB_SIZE_SCALING_FACTOR_START_IDX + 0, 1);
      //gsl_vector_set(initGuess, this->LIB_SIZE_SCALING_FACTOR_START_IDX + 1, 1);
      gsl_vector_set(initGuess, this->TRANSITION_PROB_START_IDX + 0, 0.01);
      gsl_vector_set(initGuess, this->TRANSITION_PROB_START_IDX + 1, 0.01);
      gsl_vector_set(initGuess, this->TRANSITION_PROB_START_IDX + 2, 0.2);
      gsl_vector_set(initGuess, this->TRANSITION_PROB_START_IDX + 3, 0.2);
      gsl_vector_set(initGuess, this->TRANSITION_PROB_START_IDX + 4, 0.2);
      break;

    // start lib sizes at right place, start branches even higher up
    case 4:
      for(int cellIdx = 0; cellIdx < this->NUM_LIBS_TO_EST; cellIdx++) {
        gsl_vector_set(initGuess, this->LIB_SIZE_SCALING_FACTOR_START_IDX + cellIdx, 1);
      }
      //gsl_vector_set(initGuess, this->LIB_SIZE_SCALING_FACTOR_START_IDX + 0, 1);
      //gsl_vector_set(initGuess, this->LIB_SIZE_SCALING_FACTOR_START_IDX + 1, 1);
      gsl_vector_set(initGuess, this->TRANSITION_PROB_START_IDX + 0, 0.01);
      gsl_vector_set(initGuess, this->TRANSITION_PROB_START_IDX + 1, 0.01);
      gsl_vector_set(initGuess, this->TRANSITION_PROB_START_IDX + 2, 0.4);
      gsl_vector_set(initGuess, this->TRANSITION_PROB_START_IDX + 3, 0.4);
      gsl_vector_set(initGuess, this->TRANSITION_PROB_START_IDX + 4, 0.4);
      break;

    // start lib sizes at right place, start branches and gamma higher up
    case 5:
      for(int cellIdx = 0; cellIdx < this->NUM_LIBS_TO_EST; cellIdx++) {
        gsl_vector_set(initGuess, this->LIB_SIZE_SCALING_FACTOR_START_IDX + cellIdx, 1);
      }
      //gsl_vector_set(initGuess, this->LIB_SIZE_SCALING_FACTOR_START_IDX + 0, 1);
      //gsl_vector_set(initGuess, this->LIB_SIZE_SCALING_FACTOR_START_IDX + 1, 1);
      gsl_vector_set(initGuess, this->TRANSITION_PROB_START_IDX + 0, 0.01);
      gsl_vector_set(initGuess, this->TRANSITION_PROB_START_IDX + 1, 0.01);
      gsl_vector_set(initGuess, this->TRANSITION_PROB_START_IDX + 2, 0.2);
      gsl_vector_set(initGuess, this->TRANSITION_PROB_START_IDX + 3, 0.2);
      gsl_vector_set(initGuess, this->TRANSITION_PROB_START_IDX + 4, 0.2);
      break;

    // start lib sizes at right place, start branches and gamma even higher up
    case 6:
      for(int cellIdx = 0; cellIdx < this->NUM_LIBS_TO_EST; cellIdx++) {
        gsl_vector_set(initGuess, this->LIB_SIZE_SCALING_FACTOR_START_IDX + cellIdx, 1);
      }
      //gsl_vector_set(initGuess, this->LIB_SIZE_SCALING_FACTOR_START_IDX + 0, 1);
      //gsl_vector_set(initGuess, this->LIB_SIZE_SCALING_FACTOR_START_IDX + 1, 1);
      gsl_vector_set(initGuess, this->TRANSITION_PROB_START_IDX + 0, 0.01);
      gsl_vector_set(initGuess, this->TRANSITION_PROB_START_IDX + 1, 0.1);
      gsl_vector_set(initGuess, this->TRANSITION_PROB_START_IDX + 2, 0.4);
      gsl_vector_set(initGuess, this->TRANSITION_PROB_START_IDX + 3, 0.4);
      gsl_vector_set(initGuess, this->TRANSITION_PROB_START_IDX + 4, 0.4);
      break;

    // start lib sizes at right place, start everything a little higher
    case 7:
      for(int cellIdx = 0; cellIdx < this->NUM_LIBS_TO_EST; cellIdx++) {
        gsl_vector_set(initGuess, this->LIB_SIZE_SCALING_FACTOR_START_IDX + cellIdx, 1);
      }
      //gsl_vector_set(initGuess, this->LIB_SIZE_SCALING_FACTOR_START_IDX + 0, 1);
      //gsl_vector_set(initGuess, this->LIB_SIZE_SCALING_FACTOR_START_IDX + 1, 1);
      gsl_vector_set(initGuess, this->TRANSITION_PROB_START_IDX + 0, 0.025);
      gsl_vector_set(initGuess, this->TRANSITION_PROB_START_IDX + 1, 0.1);
      gsl_vector_set(initGuess, this->TRANSITION_PROB_START_IDX + 2, 0.2);
      gsl_vector_set(initGuess, this->TRANSITION_PROB_START_IDX + 3, 0.2);
      gsl_vector_set(initGuess, this->TRANSITION_PROB_START_IDX + 4, 0.2);
      break;

    // start lib sizes at right place, start everything a little higher
    case 8:
      for(int cellIdx = 0; cellIdx < this->NUM_LIBS_TO_EST; cellIdx++) {
        gsl_vector_set(initGuess, this->LIB_SIZE_SCALING_FACTOR_START_IDX + cellIdx, 1);
      }
      //gsl_vector_set(initGuess, this->LIB_SIZE_SCALING_FACTOR_START_IDX + 0, 1);
      //gsl_vector_set(initGuess, this->LIB_SIZE_SCALING_FACTOR_START_IDX + 1, 1);
      gsl_vector_set(initGuess, this->TRANSITION_PROB_START_IDX + 0, 0.025);
      gsl_vector_set(initGuess, this->TRANSITION_PROB_START_IDX + 1, 0.1);
      gsl_vector_set(initGuess, this->TRANSITION_PROB_START_IDX + 2, 0.4);
      gsl_vector_set(initGuess, this->TRANSITION_PROB_START_IDX + 3, 0.4);
      gsl_vector_set(initGuess, this->TRANSITION_PROB_START_IDX + 4, 0.4);
      break;*/

    /*/////////////////////////////////////////////////////////////////
    // these ones are the true params picked in gridSimulation.sh
    // 13.42394947; time ./gridSimulation 100 443 0.001 0.001 0.001 0.001 0.001 2> gridSimulations/grid_443_a0.001_b0.001_g0.001_t0.001.out 1> gridSimulations/grid_443_a0.001_b0.001_g0.001_t0.001.log &
    case 9:
      gsl_vector_set(initGuess, this->LIB_SIZE_SCALING_FACTOR_START_IDX + 0, 1);
      gsl_vector_set(initGuess, this->LIB_SIZE_SCALING_FACTOR_START_IDX + 1, 1);
      gsl_vector_set(initGuess, this->TRANSITION_PROB_START_IDX + 0, 0.001);
      gsl_vector_set(initGuess, this->TRANSITION_PROB_START_IDX + 1, 0.001);
      gsl_vector_set(initGuess, this->TRANSITION_PROB_START_IDX + 2, 0.001);
      gsl_vector_set(initGuess, this->TRANSITION_PROB_START_IDX + 3, 0.001);
      gsl_vector_set(initGuess, this->TRANSITION_PROB_START_IDX + 4, 0.001);
      break;

    // 135.8353296; time ./gridSimulation 100 443 0.025 0.002 0.01 0.01 0.01 2> gridSimulations/grid_443_a0.025_b0.002_g0.01_t0.01.out 1> gridSimulations/grid_443_a0.025_b0.002_g0.01_t0.01.log &
    case 10:
      gsl_vector_set(initGuess, this->LIB_SIZE_SCALING_FACTOR_START_IDX + 0, 1);
      gsl_vector_set(initGuess, this->LIB_SIZE_SCALING_FACTOR_START_IDX + 1, 1);
      gsl_vector_set(initGuess, this->TRANSITION_PROB_START_IDX + 0, 0.025);
      gsl_vector_set(initGuess, this->TRANSITION_PROB_START_IDX + 1, 0.002);
      gsl_vector_set(initGuess, this->TRANSITION_PROB_START_IDX + 2, 0.01);
      gsl_vector_set(initGuess, this->TRANSITION_PROB_START_IDX + 3, 0.01);
      gsl_vector_set(initGuess, this->TRANSITION_PROB_START_IDX + 4, 0.01);
      break;

    // 375.8504695; time ./gridSimulation 100 443 0.025 0.002 0.1 0.01 0.01 2> gridSimulations/grid_443_a0.025_b0.002_g0.1_t0.01.out 1> gridSimulations/grid_443_a0.025_b0.002_g0.1_t0.01.log &
    case 11:
      gsl_vector_set(initGuess, this->LIB_SIZE_SCALING_FACTOR_START_IDX + 0, 1);
      gsl_vector_set(initGuess, this->LIB_SIZE_SCALING_FACTOR_START_IDX + 1, 1);
      gsl_vector_set(initGuess, this->TRANSITION_PROB_START_IDX + 0, 0.025);
      gsl_vector_set(initGuess, this->TRANSITION_PROB_START_IDX + 1, 0.002);
      gsl_vector_set(initGuess, this->TRANSITION_PROB_START_IDX + 2, 0.1);
      gsl_vector_set(initGuess, this->TRANSITION_PROB_START_IDX + 3, 0.01);
      gsl_vector_set(initGuess, this->TRANSITION_PROB_START_IDX + 4, 0.01);
      break;

    // 554.1692459; time ./gridSimulation 100 443 0.025 0.002 0.4 0.1 0.1 2> gridSimulations/grid_443_a0.025_b0.002_g0.4_t0.1.out 1> gridSimulations/grid_443_a0.025_b0.002_g0.4_t0.1.log &
    case 12:
      gsl_vector_set(initGuess, this->LIB_SIZE_SCALING_FACTOR_START_IDX + 0, 1);
      gsl_vector_set(initGuess, this->LIB_SIZE_SCALING_FACTOR_START_IDX + 1, 1);
      gsl_vector_set(initGuess, this->TRANSITION_PROB_START_IDX + 0, 0.025);
      gsl_vector_set(initGuess, this->TRANSITION_PROB_START_IDX + 1, 0.002);
      gsl_vector_set(initGuess, this->TRANSITION_PROB_START_IDX + 2, 0.4);
      gsl_vector_set(initGuess, this->TRANSITION_PROB_START_IDX + 3, 0.1);
      gsl_vector_set(initGuess, this->TRANSITION_PROB_START_IDX + 4, 0.1);
      break;

    // 844.0242963; time ./gridSimulation 100 443 0.025 0.002 0.4 0.05 0.05 2> gridSimulations/grid_443_a0.025_b0.002_g0.4_t0.05.out 1> gridSimulations/grid_443_a0.025_b0.002_g0.4_t0.05.log &
    case 13:
      gsl_vector_set(initGuess, this->LIB_SIZE_SCALING_FACTOR_START_IDX + 0, 1);
      gsl_vector_set(initGuess, this->LIB_SIZE_SCALING_FACTOR_START_IDX + 1, 1);
      gsl_vector_set(initGuess, this->TRANSITION_PROB_START_IDX + 0, 0.025);
      gsl_vector_set(initGuess, this->TRANSITION_PROB_START_IDX + 1, 0.002);
      gsl_vector_set(initGuess, this->TRANSITION_PROB_START_IDX + 2, 0.4);
      gsl_vector_set(initGuess, this->TRANSITION_PROB_START_IDX + 3, 0.05);
      gsl_vector_set(initGuess, this->TRANSITION_PROB_START_IDX + 4, 0.05);
      break;

    // 1315.07599; time ./gridSimulation 100 443 0.025 0.002 0.1 0.5 0.5 2> gridSimulations/grid_443_a0.025_b0.002_g0.1_t0.5.out 1> gridSimulations/grid_443_a0.025_b0.002_g0.1_t0.5.log &
    case 14:
      gsl_vector_set(initGuess, this->LIB_SIZE_SCALING_FACTOR_START_IDX + 0, 1);
      gsl_vector_set(initGuess, this->LIB_SIZE_SCALING_FACTOR_START_IDX + 1, 1);
      gsl_vector_set(initGuess, this->TRANSITION_PROB_START_IDX + 0, 0.025);
      gsl_vector_set(initGuess, this->TRANSITION_PROB_START_IDX + 1, 0.002);
      gsl_vector_set(initGuess, this->TRANSITION_PROB_START_IDX + 2, 0.1);
      gsl_vector_set(initGuess, this->TRANSITION_PROB_START_IDX + 3, 0.49); // TODO setting this to 0.5 is causing nan issues
      gsl_vector_set(initGuess, this->TRANSITION_PROB_START_IDX + 4, 0.49);
      break;

    // 507.1099853; time ./gridSimulation 100 443 0.025 0.002 0.025 0.25 0.25 2> gridSimulations/grid_443_a0.025_b0.002_g0.025_t0.25.out 1> gridSimulations/grid_443_a0.025_b0.002_g0.025_t0.25.log &
    case 15:
      gsl_vector_set(initGuess, this->LIB_SIZE_SCALING_FACTOR_START_IDX + 0, 1);
      gsl_vector_set(initGuess, this->LIB_SIZE_SCALING_FACTOR_START_IDX + 1, 1);
      gsl_vector_set(initGuess, this->TRANSITION_PROB_START_IDX + 0, 0.025);
      gsl_vector_set(initGuess, this->TRANSITION_PROB_START_IDX + 1, 0.002);
      gsl_vector_set(initGuess, this->TRANSITION_PROB_START_IDX + 2, 0.025);
      gsl_vector_set(initGuess, this->TRANSITION_PROB_START_IDX + 3, 0.25);
      gsl_vector_set(initGuess, this->TRANSITION_PROB_START_IDX + 4, 0.25);
      break;

    // 844.1918028; time ./gridSimulation 100 443 0.025 0.002 0.1 0.25 0.25 2> gridSimulations/grid_443_a0.025_b0.002_g0.1_t0.25.out 1> gridSimulations/grid_443_a0.025_b0.002_g0.1_t0.25.log &
    case 16:
      gsl_vector_set(initGuess, this->LIB_SIZE_SCALING_FACTOR_START_IDX + 0, 1);
      gsl_vector_set(initGuess, this->LIB_SIZE_SCALING_FACTOR_START_IDX + 1, 1);
      gsl_vector_set(initGuess, this->TRANSITION_PROB_START_IDX + 0, 0.025);
      gsl_vector_set(initGuess, this->TRANSITION_PROB_START_IDX + 1, 0.002);
      gsl_vector_set(initGuess, this->TRANSITION_PROB_START_IDX + 2, 0.1);
      gsl_vector_set(initGuess, this->TRANSITION_PROB_START_IDX + 3, 0.25);
      gsl_vector_set(initGuess, this->TRANSITION_PROB_START_IDX + 4, 0.25);
      break;

    // 758.3894924; time ./gridSimulation 100 443 0.025 0.002 0.2 0.1 0.1 2> gridSimulations/grid_443_a0.025_b0.002_g0.2_t0.1.out 1> gridSimulations/grid_443_a0.025_b0.002_g0.2_t0.1.log &
    case 17:
      gsl_vector_set(initGuess, this->LIB_SIZE_SCALING_FACTOR_START_IDX + 0, 1);
      gsl_vector_set(initGuess, this->LIB_SIZE_SCALING_FACTOR_START_IDX + 1, 1);
      gsl_vector_set(initGuess, this->TRANSITION_PROB_START_IDX + 0, 0.025);
      gsl_vector_set(initGuess, this->TRANSITION_PROB_START_IDX + 1, 0.002);
      gsl_vector_set(initGuess, this->TRANSITION_PROB_START_IDX + 2, 0.2);
      gsl_vector_set(initGuess, this->TRANSITION_PROB_START_IDX + 3, 0.1);
      gsl_vector_set(initGuess, this->TRANSITION_PROB_START_IDX + 4, 0.1);
      break;

    // 888.4588808; time ./gridSimulation 100 443 0.025 0.002 0.3 0.1 0.1 2> gridSimulations/grid_443_a0.025_b0.002_g0.3_t0.1.out 1> gridSimulations/grid_443_a0.025_b0.002_g0.3_t0.1.log &
    case 18:
      gsl_vector_set(initGuess, this->LIB_SIZE_SCALING_FACTOR_START_IDX + 0, 1);
      gsl_vector_set(initGuess, this->LIB_SIZE_SCALING_FACTOR_START_IDX + 1, 1);
      gsl_vector_set(initGuess, this->TRANSITION_PROB_START_IDX + 0, 0.025);
      gsl_vector_set(initGuess, this->TRANSITION_PROB_START_IDX + 1, 0.002);
      gsl_vector_set(initGuess, this->TRANSITION_PROB_START_IDX + 2, 0.3);
      gsl_vector_set(initGuess, this->TRANSITION_PROB_START_IDX + 3, 0.1);
      gsl_vector_set(initGuess, this->TRANSITION_PROB_START_IDX + 4, 0.1);
      break;
      */
  }
}

/*
 * function to solve the overdetermined system of equations from the Baum Welch estimated transition matrix, a*
 * Assumes runBaumWelch() has already been called
 */
//void TwoCell3TrParam2DegPolyHMM::solveBaumWelchEstimates() {
void TwoCell3TrParam2DegPolyHMM::setUpBaumWelchLeastSquares() {
  this->baumWelchTransitionMat = gsl_matrix_alloc(this->transition->size1, this->transition->size2);
//  //std::cout << "TwoCell3TrParam2DegPolyHMM::setUpBaumWelchLeastSquares" << std::endl;
//  /*
//   * transition matrix entries are of the form:
//   * (b+a+g)*t1
//   * (b+a+g)*t2
//   * (b+a+g)*t3
//   * (b+a)*t1
//   * (b+a)*t2
//   * (b+a)*t3
//   * (b+g)*t1
//   * (b+g)*t2
//   * (b+g)*t3
//   * (b)*t1
//   * (b)*t2
//   * (b)*t3
//   * which can be rewritten as f(.) = (xa + yb + zg) * (ut1 + vt2 + wt3)
//   * where u,v,w,x,y,z are all indicator variables.
//   * 
//   * We want to solve the overdetermined system by BFGS on least squares
//   * S = sum_i r_i^2
//   * r_i = a_ij* - f(.)
//   */
//
//  // one residual per coef entry, which has one entry per transition matrix element
//  this->residuals = gsl_vector_alloc(this->states->size() * this->states->size());
//  gsl_vector_set_zero(this->residuals);
//
//  /*
//   * structure of coefs:
//   *   x  y  z  u  v  w
//   * [                  ] // a*[0,0]
//   * [                  ] // a*[0,1]
//   *         ...
//   * [                  ] // a*[k,k]
//   */
//  int xIdx = 0;
//  int yIdx = 1;
//  int zIdx = 2;
//  int uIdx = 3;
//  int vIdx = 4;
//  int wIdx = 5;
//  this->coefs = gsl_matrix_alloc(this->states->size() * this->states->size(), 1 + this->NUM_TRANSITION_PARAMS_TO_EST + this->NUM_BRANCH_LENGTHS_TO_EST); // +1 for alpha coef
//  gsl_matrix_set_zero(coefs);
//  //printMatrix(coefs);
//
//  // fill in coefs. logic taken from setTransition
//  int frA = 0;
//  int frB = 0;
//  int toA = 0;
//  int toB = 0;
//  int coefsIdx = 0;
//  for(unsigned int row = 0; row < this->states->size(); row++) {
//    frA = getCellPloidyFromStateIdx(0, row);
//    frB = getCellPloidyFromStateIdx(1, row);
//    for(unsigned int col = 0; col < this->states->size(); col++, coefsIdx++) {
//      toA = getCellPloidyFromStateIdx(0, col);
//      toB = getCellPloidyFromStateIdx(1, col);
//
//      // if frA == toA (A does not move)
//      if(frA == toA) {
//        // if frB == toB (B does not move)
//        if(frB == toB) {
//          // diagonal; do nothing and continue to next
//          continue;
//        }
//        // any CNA for B
//        gsl_matrix_set(coefs, coefsIdx, yIdx, 1);
//        // if abs(frB - toB) == 1 (B adjacent move)
//        if(std::abs(frB - toB) == 1) { // adj CNA (alpha)
//          gsl_matrix_set(coefs, coefsIdx, xIdx, 1);
//        }
//        // if toB == 2 (return to diploid (gamma))
//        if(toB == 2) {
//          gsl_matrix_set(coefs, coefsIdx, zIdx, 1);
//        }
//        // A stays but B moves (t3 branch)
//        gsl_matrix_set(coefs, coefsIdx, wIdx, 1);
//      }
//      // else (A has moved)
//      else {
//        // if frB == toB (B does not move)
//        if(frB == toB) {
//          gsl_matrix_set(coefs, coefsIdx, yIdx, 1); // any CNA (beta)
//          // if abs(frA - toA) == 1 (A adjacent move)
//          if(std::abs(frA - toA) == 1) { // adj CNA (alpha)
//            gsl_matrix_set(coefs, coefsIdx, xIdx, 1);
//          }
//          // if toA == 2 (return to diploid (gamma))
//          if(toA == 2) {
//            gsl_matrix_set(coefs, coefsIdx, zIdx, 1);
//          }
//          // A moves but B stays (t2 branch)
//          gsl_matrix_set(coefs, coefsIdx, vIdx, 1);
//        }
//        // else if frA == frB && toA == toB (move in parallel)
//        else if(frA == frB && toA == toB) {
//          gsl_matrix_set(coefs, coefsIdx, yIdx, 1); // any CNA (beta)
//          // if abs(frA - toA) == 1 (A and B adjacent move)
//          if(std::abs(frA - toA) == 1) { // adj CNA (alpha)
//            gsl_matrix_set(coefs, coefsIdx, xIdx, 1);
//          }
//          // if toA == 2 (return to diploid (gamma))
//          if(toA == 2) {
//            gsl_matrix_set(coefs, coefsIdx, zIdx, 1);
//          }
//          // both A and B move (t1 branch)
//          gsl_matrix_set(coefs, coefsIdx, uIdx, 1);
//        }
//        // else, some impossible transition
//        else {
//          // do nothing
//        }
//        //gsl_vector_view currRow = gsl_matrix_row(coefs, coefsIdx);
//        //printRowVector(&currRow.vector);
//      }
//    }
//  }
//
//  this->a_ij = gsl_vector_alloc(this->states->size() * this->states->size());
//  int a_ijIdx = 0;
//  for(unsigned int row = 0; row < this->states->size(); row++) {
//    for(unsigned int col = 0; col < this->states->size(); col++, a_ijIdx++) {
//      gsl_vector_set(a_ij, a_ijIdx, gsl_matrix_get(this->transition, row, col));
//    }
//  }
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
  /*if(gsl_vector_max(probs) > 200) {
    return GSL_NAN; // makes matrices singular
    //return 1e20;
  }*/

  //std::cout << "probs" << std::endl;
  //printColVector(probs);
  //std::cout << "this->baumWelchTransitionMat" << std::endl;
  //printMatrix(this->baumWelchTransitionMat);
  //std::cout << "this->transition" << std::endl;
  //printMatrix(this->transition);
  // compare each entry in the hypothetical proposed transition matrix to the baum welch calculated transition matrix (ie sum up the squared difference)
  double sumSqResid = 0;
  double resid = 0;
  //std::cout << "this->baumWelchTransitionMat" << std::endl;
  //printMatrix(this->baumWelchTransitionMat);
  //std::cout << "this->transition" << std::endl;
  //printMatrix(this->transition);
  for(unsigned int row = 0; row < this->transition->size1; row++) {
    //gsl_vector_view bwRow = gsl_matrix_subrow(this->baumWelchTransitionMat, row, 0, this->baumWelchTransitionMat->size2 - 1);
    //gsl_vector_view trRow = gsl_matrix_subrow(this->transition, row, 0, this->transition->size2 - 1);
    //resid = gsl_stats_correlation( (double*) bwRow.vector.data, 1, (double*) trRow.vector.data, 1, this->baumWelchTransitionMat->size2 - 1);
    //gsl_vector_view bwRow = gsl_matrix_subrow(this->baumWelchTransitionMat, row, 0, this->baumWelchTransitionMat->size2);
    //gsl_vector_view trRow = gsl_matrix_subrow(this->transition, row, 0, this->transition->size2);
    //resid = gsl_stats_correlation(bwRow.vector.data, 1, trRow.vector.data, 1, this->baumWelchTransitionMat->size2);
    //sumSqResid += resid;
    //double bwRowSum = 0;
    //double trRowSum = 0;
    //for(unsigned int col = 0; col < this->transition->size2; col++) {
    //std::cout << "row: " << row << std::endl;
    for(unsigned int col = 0; col < this->transition->size2 - 1; col++) {
      //if(row == col) continue;
      //resid = (gsl_matrix_get(this->baumWelchTransitionMat, row, col) - gsl_matrix_get(this->transition, row, col));
      //sumSqResid += resid * resid;
      double bwVal = gsl_matrix_get(this->baumWelchTransitionMat, row, col);
      double trVal = gsl_matrix_get(this->transition, row, col);
      //bwRowSum += bwVal;
      //trRowSum += trVal;
      resid = ((bwVal - trVal) * (bwVal - trVal));// / bwVal;
      //resid = std::sqrt((bwVal - trVal) * (bwVal - trVal));// / bwVal;
      sumSqResid += resid;
      //std::cout << "bwVal: " << bwVal << ", trVal: " << trVal << ", diff: " << bwVal - trVal << ", resid: " << resid << ", sumSqResid: " << sumSqResid << std::endl;
    }
    //bwRowSum /= (double) this->transition->size2;
    //trRowSum /= (double) this->transition->size2;
    //for(unsigned int col = 0; col < this->transition->size2; col++) {
    //  if(row == col) continue;
    //  double bwVal = gsl_matrix_get(this->baumWelchTransitionMat, row, col);
    //  double trVal = gsl_matrix_get(this->transition, row, col);
    //  resid = (bwVal - bwRowSum) * (trVal - trRowSum);
    //  sumSqResid += resid;
    //}
  }
  //sumSqResid = sqrt(sumSqResid);
  //sumSqResid = calcChiSqOfMatrix(this->transition, this->baumWelchTransitionMat); // (exp, obs)
  //sumSqResid = calcChiSqOfMatrix(this->baumWelchTransitionMat, this->transition); // (exp, obs)
  //std::cout << "sumSqResid: " << sumSqResid << std::endl;
  return sumSqResid;


  //double transitionProbs = 0;
  //double branchLengths = 0;
  //double coefVarProd = 0;
  //for(unsigned int residIdx = 0; residIdx < residuals->size; residIdx++) {
  //  // skip over the diagonals
  //  if(residIdx % (this->states->size() + 1) == 0) {
  //    //std::cout << "DIAG: " << residIdx << std::endl;
  //    resid = 0;
  //  }
  //  else {
  //    transitionProbs  = gsl_matrix_get(coefs, residIdx, 0) * this->getAlpha(); // alpha
  //    transitionProbs += gsl_matrix_get(coefs, residIdx, 1) * gsl_vector_get(probs, 0); // beta
  //    transitionProbs += gsl_matrix_get(coefs, residIdx, 2) * gsl_vector_get(probs, 1); // gamma

  //    branchLengths  = gsl_matrix_get(coefs, residIdx, 3) * gsl_vector_get(probs, 2); // t1
  //    branchLengths += gsl_matrix_get(coefs, residIdx, 4) * gsl_vector_get(probs, 3); // t2
  //    branchLengths += gsl_matrix_get(coefs, residIdx, 5) * gsl_vector_get(probs, 4); // t3
  //    /*for(unsigned int coefIdx = 1; coefIdx < coefs->size2; coefIdx++) {
  //      std::cout << "coef[" <<coefIdx << "]: " << gsl_matrix_get(coefs, residIdx, coefIdx) << ", probs[" << coefIdx-1 << "]: " << gsl_vector_get(probs, coefIdx-1) << std::endl;
  //      //std::cout << "coefVarProd: " << coefVarProd << ", residIdx " << residIdx << ", coefIdx " << coefIdx << ", coefs[,] " << gsl_matrix_get(coefs, residIdx, coefIdx) << ", v[] " << gsl_vector_get(v, coefIdx - 1) << std::endl;
  //      //coefVarProd += gsl_matrix_get(coefs, residIdx, coefIdx) * gsl_vector_get(probs, coefIdx - 1);
  //    }*/
  //    //std::cout << "transitionProbs: " << transitionProbs << ", branchLengths: " << branchLengths << std::endl;
  //    coefVarProd = transitionProbs * branchLengths;
  //    resid = gsl_vector_get(a_ij, residIdx) - coefVarProd;
  //  }
  //  gsl_vector_set(residuals, residIdx, resid);
  //  //std::cout << "resid: "  << resid << std::endl;
  //  sumSqResid += resid * resid;
  //  //std::cout << "sumSqResid: "  << sumSqResid << std::endl;
  //}
  /*std::cout << "coefs: " << std::endl;
  printMatrix(coefs);
  std::cout << "residuals: " << std::endl;
  printColVector(residuals);
  std::cout << "a_ij: " << std::endl;
  printColVector(a_ij);
  std::cout << "probs: " << std::endl;
  printColVector((gsl_vector*)probs);*/

  //return sumSqResid;
}
/*
 * function to save all parameters to file for a specific cell.
 * Always saves in order of [lib0, lib1, alpha, beta, lambda, t1, t2, t3]
 */
//void TwoCell3TrParam2DegPolyHMM::saveParamEstimates(std::string filename) const {
//  // TODO how to generate filename for two cells to also make it easily findable?
//}

/*
 * method to check validity of parameters proposed by BFGS. Assumes probs
 * is in probability space.
 * returns 0 if probs is valid, GSL_NAN otherwise
 */
double TwoCell3TrParam2DegPolyHMM::checkOptimProbValidity(gsl_vector* probs) const {
  return 0; // Thu 27 Aug 2020 08:27:51 PM PDT debugging no validity check
  //std::cout << " TwoCell3TrParam2DegPolyHMM::checkOptimProbValidity" << std::endl;
  // shortcut for bad library sizes (too large or too small)
  double currLibSizeScalingFactor = -1;
  for(int i = 0; i < this->NUM_LIBS_TO_EST; i++) {
    currLibSizeScalingFactor = gsl_vector_get(probs, this->LIB_SIZE_SCALING_FACTOR_START_IDX + i);
    if(gsl_isinf(currLibSizeScalingFactor) || currLibSizeScalingFactor < 1e-2 || currLibSizeScalingFactor > 1e2) {
      std::cerr << "shortcut for bad lib sizes: ";// << std::endl;
      //printRowVector(stderr, probs);
      return GSL_NAN;
    }
  }

  // shortcut for any probabilities becoming too small
  double probMin = gsl_vector_min(probs);
  if(probMin < 1e-5 || gsl_isnan(probMin)) {
    //std::cerr << "shortcut for any probabilities becoming too small: " ;//<< probMin << std::endl;
    //printRowVector(stderr, probs);
    return GSL_NAN;
  }
  return 0;
}


