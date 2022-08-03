#include "SelectPairs0TrParam2DegPolyHMM.hpp"

// ctors and destructor
SelectPairs0TrParam2DegPolyHMM::SelectPairs0TrParam2DegPolyHMM(std::vector<DepthPair*>* depths, std::vector<std::string>* sampleList, gsl_vector* fixedParams, int maxPloidy, int numFixedLibs, int numPairs, int numPairsToSummarize, int ordering, unsigned int seed, gsl_matrix* stage1NearestCellIdxMat, int numBranchesToEst) : AllPairs0TrParam2DegPolyHMM(depths, sampleList, fixedParams, maxPloidy, numFixedLibs, numPairs, numBranchesToEst) {
  this->numPairsToSummarize = numPairsToSummarize;
  this->ordering = ordering;
  this->seed = seed;
  this->stage1NearestCellIdxMat = stage1NearestCellIdxMat;
}

SelectPairs0TrParam2DegPolyHMM* SelectPairs0TrParam2DegPolyHMM::create(std::vector<DepthPair*>* depths, std::vector<std::string>* sampleList, gsl_vector* fixedParams, int maxPloidy, int numPairs, int numPairsToSummarize, int ordering, unsigned int seed, gsl_matrix* stage1NearestCellIdxMat, gsl_vector* meanVarianceCoefVec, bool preallocIntermediates) {
  SelectPairs0TrParam2DegPolyHMM* hmm = new SelectPairs0TrParam2DegPolyHMM(depths, sampleList, fixedParams, maxPloidy, 0, numPairs, numPairsToSummarize, ordering, seed, stage1NearestCellIdxMat, numPairs * 3); // 0 numFixedLibs, numPairs*3 numBranchesToEst
  hmm->makeHMMPairs(meanVarianceCoefVec, preallocIntermediates);
  return hmm;
}

SelectPairs0TrParam2DegPolyHMM::~SelectPairs0TrParam2DegPolyHMM() {
  // TODO
}

void SelectPairs0TrParam2DegPolyHMM::makeHMMPairs(gsl_vector* meanVarianceCoefVec, gsl_vector* transitionParams, bool preallocIntermediates) {
  // first init all entries of hmmVec to nullptr
  for(unsigned int hmmIdx = 0; hmmIdx < this->hmmVec->size(); hmmIdx++) {
    (*this->hmmVec)[hmmIdx] = nullptr;
  }

  //std::vector<DepthPair*>* currDepths = nullptr;
  //gsl_vector* currMeanVarCoefVec = nullptr; // make them all have their own copies of this vector
  //gsl_vector* currFixedParams = nullptr;
  int hmmIdx = 0;
  //int refCell = 0;
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
  for(int i = 0; i < (int) this->depthsVec->size() && hmmIdx < this->NUM_PAIRS; i++) {
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
      this->storeHMMIdxForCell(i, hmmIdx);
      if((*this->hmmVec)[hmmIdx] != nullptr) {
        continue;
      }

      this->makeOneHMMPair(leftCell, rightCell, preallocIntermediates);
      //currDepths = new std::vector<DepthPair*>();
      ////currDepths->push_back((*this->depthsVec)[refCell]);
      ////currDepths->push_back((*this->depthsVec)[nearestCell]);
      //currDepths->push_back((*this->depthsVec)[leftCell]);
      //currDepths->push_back((*this->depthsVec)[rightCell]);

      //currFixedParams = gsl_vector_alloc(this->fixedParams->size);
      //gsl_vector_memcpy(currFixedParams, this->fixedParams);
      //TwoCell0TrParam2DegPolyHMM* hmm = new TwoCell0TrParam2DegPolyHMM(currDepths, currFixedParams, this->getKploidy(), preallocIntermediates);
      //(*this->hmmVec)[hmmIdx] = hmm;

      //// do rest of HMM set up (usually happens in main.cpp)
      //currMeanVarCoefVec = gsl_vector_alloc(meanVarianceCoefVec->size);
      //gsl_vector_memcpy(currMeanVarCoefVec, meanVarianceCoefVec);
      //hmm->setMeanVarianceFn(currMeanVarCoefVec);
      //hmm->setTransition(transitionParams); // this vector isn't saved anywhere
      //hmm->setLibScalingFactorsToTotalRatio();
      //hmm->setAlpha(this->getAlpha());

      ////this->setLibScalingFactors(refCell, nearestCell, hmm->getLibScalingFactor(0), hmm->getLibScalingFactor(1));
      //this->setLibScalingFactors(leftCell, rightCell, hmm->getLibScalingFactor(0), hmm->getLibScalingFactor(1));
      hmmIdx++;
    }
  }
  delete randomOrdering;
}

void SelectPairs0TrParam2DegPolyHMM::print(FILE* stream) {
  AllCells3TrParam2DegPolyHMM::print(stream);
  fprintf(stream, "SelectPairs0TrParam2DegPolyHMM numPairsToSummarize: %i\n\n", this->numPairsToSummarize);
  if(this->stage1NearestCellIdxMat != nullptr) {
    fprintf(stream, "SelectPairsFixLib0TrParam2DegPolyHMM stage1NearestCellIdxMat:\n\n");
    printMatrix(stream, this->stage1NearestCellIdxMat);
  }
}

//void SelectPairs0TrParam2DegPolyHMM::decodeStatSummaryOneCellAcrossPairs(int summaryMethod, int cellNum) {
//  // TODO HERE
//}

// override specifically so can call resetSkippedParams
SelectPairs0TrParam2DegPolyHMM* SelectPairs0TrParam2DegPolyHMM::bfgs(gsl_vector* initGuess, int maxIters, bool verbose, bool debug) {
  SelectPairs0TrParam2DegPolyHMM* bestGuess = (SelectPairs0TrParam2DegPolyHMM*) AllPairs0TrParam2DegPolyHMM::bfgs(initGuess, maxIters, verbose, debug);
  bestGuess->resetSkippedParams();
  return bestGuess;
}

/*
 * function to reset values in paramsToEst that weren't actually estimated. that is, there wasn't an HMM corresponding to those branch lengths, so setting to -1 to be clear they weren't set
 * copied from SelectPairsFixLib0TrParam2DegPolyHMM
 */
void SelectPairs0TrParam2DegPolyHMM::resetSkippedParams() {
  for(unsigned int hmmIdx = 0; hmmIdx < this->hmmVec->size(); hmmIdx++) {
    if((*this->hmmVec)[hmmIdx] == nullptr) {
      gsl_vector_set(this->paramsToEst, this->BRANCH_LENGTH_START_IDX + 3 * hmmIdx + 0, -1); // t1
      gsl_vector_set(this->paramsToEst, this->BRANCH_LENGTH_START_IDX + 3 * hmmIdx + 1, -1); // t2
      gsl_vector_set(this->paramsToEst, this->BRANCH_LENGTH_START_IDX + 3 * hmmIdx + 2, -1); // t3
    }
  }
}

