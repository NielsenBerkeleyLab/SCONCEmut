#include "SelectPairsFixLib0TrParam2DegPolyHMM.hpp"

// constructors and destructors
SelectPairsFixLib0TrParam2DegPolyHMM::SelectPairsFixLib0TrParam2DegPolyHMM(std::vector<DepthPair*>* depths, std::vector<std::string>* sampleList, gsl_vector* fixedParams, int maxPloidy, int numPairs, int numPairsToSummarize, int ordering, unsigned int seed, gsl_matrix* stage1NearestCellIdxMat, int numBranchesToEst) : AllPairsFixLib0TrParam2DegPolyHMM(depths, sampleList, fixedParams, maxPloidy, numPairs, numBranchesToEst) {
  this->numPairsToSummarize = numPairsToSummarize;
  this->ordering = ordering;
  this->seed = seed;
  this->stage1NearestCellIdxMat = stage1NearestCellIdxMat;
}

SelectPairsFixLib0TrParam2DegPolyHMM* SelectPairsFixLib0TrParam2DegPolyHMM::create(std::vector<DepthPair*>* depths, std::vector<std::string>* sampleList, gsl_vector* fixedParams, int maxPloidy, int numPairs, int numPairsToSummarize, int ordering, unsigned int seed, gsl_matrix* stage1NearestCellIdxMat, gsl_vector* meanVarianceCoefVec, bool preallocIntermediates) {
  SelectPairsFixLib0TrParam2DegPolyHMM* hmm = new SelectPairsFixLib0TrParam2DegPolyHMM(depths, sampleList, fixedParams, maxPloidy, numPairs, numPairsToSummarize, ordering, seed, stage1NearestCellIdxMat, numPairs * 3); // numPairs*3 numBranchesToEst
  hmm->makeHMMPairs(meanVarianceCoefVec, preallocIntermediates);
  return hmm;
}

void SelectPairsFixLib0TrParam2DegPolyHMM::makeHMMPairs(gsl_vector* meanVarianceCoefVec, gsl_vector* transitionParams, bool preallocIntermediates) {
  // first init all entries of hmmVec to nullptr
  for(unsigned int hmmIdx = 0; hmmIdx < this->hmmVec->size(); hmmIdx++) {
    (*this->hmmVec)[hmmIdx] = nullptr;
    (*this->shouldCallBFGSOnHmmIdx)[hmmIdx] = false;
  }

  //std::vector<DepthPair*>* currDepths = nullptr;
  //gsl_vector* currMeanVarCoefVec = nullptr; // make them all have their own copies of this vector
  //gsl_vector* currFixedParams = nullptr;
  int hmmIdx = 0;
  //int refCell = 0; // refCell and nearestCell do not change
  int nearestCell = 0; // nearestCell does not change
  int leftCell = 0; // leftCell and rightCell might swap
  int rightCell = 0;
  int swap = 0;

  std::vector<int>* randomOrdering = nullptr;
  if(this->ordering == IndThenPairs2Stages3TrParam2DegPolyHMM::ORDERING_RANDOM) {
    if(this->generator == nullptr) {
      this->generator = new base_generator_type(this->seed);
    }
  }

  // for each cell
  //std::cout << "depths size: " << this->depthsVec->size() << ", (int) depths size: " << (int) this->depthsVec->size() << std::endl;
  //std::cout << "nearestCellMat size2: " << this->stage1NearestCellIdxMat->size2 << ", (int) nearestCellMat size2: " << (int) this->stage1NearestCellIdxMat->size2 << std::endl;
  bool atLeastOneHmmCreated = false;
  for(int i = 0; i < (int) this->depthsVec->size(); i++) {
    // regen a random ordering if necessary
    if(this->ordering == IndThenPairs2Stages3TrParam2DegPolyHMM::ORDERING_RANDOM) {
      randomOrdering = genUniqueRandomSeqInRange(0, this->stage1NearestCellIdxMat->size2, i, this->generator);
    }
    // init hmm's for its numPairsToSummarize nearest cells (unless already alloc)
    for(int nearestCellIdx = 0; nearestCellIdx < this->numPairsToSummarize && nearestCellIdx < (int) this->stage1NearestCellIdxMat->size2; nearestCellIdx++) {
      //refCell = i;
      if(this->ordering == IndThenPairs2Stages3TrParam2DegPolyHMM::ORDERING_RANDOM) {
        nearestCell = (*randomOrdering)[nearestCellIdx];
      }
      else if(this->ordering == IndThenPairs2Stages3TrParam2DegPolyHMM::ORDERING_FURTHEST) {
        nearestCell = gsl_matrix_get(this->stage1NearestCellIdxMat, i, this->stage1NearestCellIdxMat->size2 - nearestCellIdx - 1);
      }
      else { // IndThenPairs2Stages3TrParam2DegPolyHMM::ORDERING_NEAREST is default
        nearestCell = gsl_matrix_get(this->stage1NearestCellIdxMat, i, nearestCellIdx);
      }

      // first cellidx should be smaller than second since this calculates an *upper* triangular matrix
      // nearestCell will never contain self, since that mat is constructed to exclude self
      leftCell = i;
      rightCell = nearestCell;
      if(leftCell > rightCell) {
        swap = leftCell;
        leftCell = rightCell;
        rightCell = swap;
      }
      //if(refCell > nearestCell) {
      //  swap = refCell;
      //  refCell = nearestCell;
      //  nearestCell = swap;
      //}
      //hmmIdx = this->getHMMIdxFromCellPair(refCell, nearestCell);
      hmmIdx = this->getHMMIdxFromCellPair(leftCell, rightCell);
      //std::cout << "row: " << i << ", refCell: " << refCell << ", nearestCellIdx: " << nearestCellIdx << ", nearestCell: " << nearestCell << ", hmmIdx: " << hmmIdx << std::endl;
      this->storeHMMIdxForCell(i, hmmIdx);
      (*this->shouldCallBFGSOnHmmIdx)[hmmIdx] = true;

      // always create at least one hmm (needed for indexing constants)
      if(!atLeastOneHmmCreated) {
        this->makeOneHMMPair(leftCell, rightCell, preallocIntermediates);
        atLeastOneHmmCreated = true;
      }
      // if preallocating intermediates, set it up
      else if(preallocIntermediates) {
        if((*this->hmmVec)[hmmIdx] != nullptr) {
          //std::cout << "already exists, continue" << std::endl;
          continue;
        }
        this->makeOneHMMPair(leftCell, rightCell, preallocIntermediates);
      }

      //currDepths = new std::vector<DepthPair*>();
      ////currDepths->push_back((*this->depthsVec)[refCell]);
      ////currDepths->push_back((*this->depthsVec)[nearestCell]);
      //currDepths->push_back((*this->depthsVec)[leftCell]);
      //currDepths->push_back((*this->depthsVec)[rightCell]);

      //currFixedParams = gsl_vector_alloc(2 + 3); // 2 libs, alpha/beta/lambda
      //// set libs
      ////gsl_vector_set(currFixedParams, 0, gsl_vector_get(this->fixedParams, this->LIB_SIZE_SCALING_FACTOR_START_IDX + i));
      ////gsl_vector_set(currFixedParams, 1, gsl_vector_get(this->fixedParams, this->LIB_SIZE_SCALING_FACTOR_START_IDX + nearestCell));
      //gsl_vector_set(currFixedParams, 0, gsl_vector_get(this->fixedParams, this->LIB_SIZE_SCALING_FACTOR_START_IDX + leftCell));
      //gsl_vector_set(currFixedParams, 1, gsl_vector_get(this->fixedParams, this->LIB_SIZE_SCALING_FACTOR_START_IDX + rightCell));

      //// set shared transition params
      //gsl_vector_set(currFixedParams, 2, gsl_vector_get(this->fixedParams, this->FIXED_TRANSITION_PROB_START_IDX + 0)); // alpha
      //gsl_vector_set(currFixedParams, 3, gsl_vector_get(this->fixedParams, this->FIXED_TRANSITION_PROB_START_IDX + 1)); // beta
      //gsl_vector_set(currFixedParams, 4, gsl_vector_get(this->fixedParams, this->FIXED_TRANSITION_PROB_START_IDX + 2)); // lambda

      //TwoCellFixLib0TrParam2DegPolyHMM* hmm = new TwoCellFixLib0TrParam2DegPolyHMM(currDepths, currFixedParams, this->getKploidy(), preallocIntermediates);
      //(*this->hmmVec)[hmmIdx] = hmm;

      //// do rest of HMM set up (usually happens in main.cpp)
      //currMeanVarCoefVec = gsl_vector_alloc(meanVarianceCoefVec->size);
      //gsl_vector_memcpy(currMeanVarCoefVec, meanVarianceCoefVec);
      //hmm->setMeanVarianceFn(currMeanVarCoefVec);
      //hmm->setTransition(transitionParams); // this vector isn't saved anywhere
      //hmm->setAlpha(this->getAlpha());

      ////this->setLibScalingFactors(refCell, nearestCell, hmm->getLibScalingFactor(0), hmm->getLibScalingFactor(1));
      //this->setLibScalingFactors(leftCell, rightCell, hmm->getLibScalingFactor(0), hmm->getLibScalingFactor(1));
      ////this->storeHMMIdxForCells(refCell, nearestCell, hmmIdx);
    }
  }
  delete randomOrdering;
}

SelectPairsFixLib0TrParam2DegPolyHMM::~SelectPairsFixLib0TrParam2DegPolyHMM() {
  // TODO
}

void SelectPairsFixLib0TrParam2DegPolyHMM::print(FILE* stream) {
  AllCells3TrParam2DegPolyHMM::print(stream);
  fprintf(stream, "SelectPairsFixLib0TrParam2DegPolyHMM numPairsToSummarize: %i\n\n", this->numPairsToSummarize);
  if(this->stage1NearestCellIdxMat != nullptr) {
    fprintf(stream, "SelectPairsFixLib0TrParam2DegPolyHMM stage1NearestCellIdxMat:\n\n");
    printMatrix(stream, this->stage1NearestCellIdxMat);
  }
}

// override specifically so can call resetSkippedParams
SelectPairsFixLib0TrParam2DegPolyHMM* SelectPairsFixLib0TrParam2DegPolyHMM::bfgs(gsl_vector* initGuess, std::string filename, int maxIters, bool verbose, bool debug) {
  //std::cout << "SelectPairsFixLib0TrParam2DegPolyHMM::bfgs" << std::endl;
  SelectPairsFixLib0TrParam2DegPolyHMM* bestGuess = (SelectPairsFixLib0TrParam2DegPolyHMM*) AllPairsFixLib0TrParam2DegPolyHMM::bfgs(initGuess, filename, maxIters, verbose, debug);
  bestGuess->resetSkippedParams();
  return bestGuess;
}

/*
 * function to reset values in paramsToEst that weren't actually estimated. that is, there wasn't an HMM corresponding to those branch lengths, so setting to -1 to be clear they weren't set by BFGS estimation. TODO perhaps don't reset, in case want to keep bw+ls ests?
 */
void SelectPairsFixLib0TrParam2DegPolyHMM::resetSkippedParams() {
  //std::cout << "SelectPairsFixLib0TrParam2DegPolyHMM::resetSkippedParams" << std::endl;
  for(unsigned int hmmIdx = 0; hmmIdx < this->hmmVec->size(); hmmIdx++) {
    if((*this->hmmVec)[hmmIdx] == nullptr) {
      gsl_vector_set(this->paramsToEst, this->BRANCH_LENGTH_START_IDX + 3 * hmmIdx + 0, -1); // t1
      gsl_vector_set(this->paramsToEst, this->BRANCH_LENGTH_START_IDX + 3 * hmmIdx + 1, -1); // t2
      gsl_vector_set(this->paramsToEst, this->BRANCH_LENGTH_START_IDX + 3 * hmmIdx + 2, -1); // t3
    }
  }
}
