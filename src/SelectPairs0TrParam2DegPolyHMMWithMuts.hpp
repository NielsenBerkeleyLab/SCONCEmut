#ifndef SELECTPAIRS0TRPARAM2DEGPOLYHMMWITHMUTS_HPP
#define SELECTPAIRS0TRPARAM2DEGPOLYHMMWITHMUTS_HPP

#include "SelectPairs0TrParam2DegPolyHMM.hpp"
#include "TwoCell0TrParam2DegPolyHMMWithMuts.hpp"
#include "MutationList.hpp"
#include "MutationPair.hpp"

/*
 * This class is for making an HMM that has fixed alpha/beta/lambda values and estimates branches
 * and library scaling sizes only across pairs of tumor cells. This is similar to AllPairs3TrParam2DegPolyHMM,
 * but all the transition params are fixed.
 * This class also incorporates mutations, by TwoCell0TrParam2DegPolyHMMWithMuts (and MutationList and MutationPair internally)
 *
 * Specifically, only entries that are selected to be optimized have non null values in this->hmmVec. this->hmmVec still has the full size expected (ie numPairs).
 * this->paramsToEst also has the full paramset, even if they're never used in an HMM (this will be helpful if we want to alloc HMMs on the fly, as all param values will already be stored)
 * Selected pairs are determined by the numPairsToSummarize'th nearest pairs for each cell, as determined by stage1NearestCellIdxMat.
 * The selection process is done in makeHMMPairs.
 *
 * this->paramsToEst = [lib0, lib1, ..., libN, t1_cell0_1, t2_cell0_1, t3_cell0_1, t1_cell0_2, t2_cell0_2, t3_cell0_2, ..., t3_cell(N-1)_N]
 * this->fixedParams = [alpha, beta, lambda, cnaToMutRateMu]
 */
class SelectPairs0TrParam2DegPolyHMMWithMuts : public SelectPairs0TrParam2DegPolyHMM {
  private:
  protected:
    // member variables
    std::vector<MutationList*>* mutListVec;
    std::vector<MutationPair*>* mutPairVec;
    gsl_vector* allIndEstMutCounts;

    SelectPairs0TrParam2DegPolyHMMWithMuts(std::vector<DepthPair*>* depths, std::vector<std::string>* sampleList, std::vector<MutationList*>* mutListVec, gsl_vector* allIndEstMutCounts, gsl_vector* fixedParams, int maxPloidy, int numFixedLibs, int numPairs, int numPairsToSummarize, int ordering, unsigned int seed, gsl_matrix* stage1NearestCellIdxMat, int numBranchesToEst);

    virtual void makeOneHMMPair(int i, int j, bool preallocIntermediates = true);
    using AllPairs0TrParam2DegPolyHMM::makeHMMPairs; // unhide parent method of same name https://stackoverflow.com/a/18100999

  public:
    // constructors and destructor
    virtual ~SelectPairs0TrParam2DegPolyHMMWithMuts();
    static SelectPairs0TrParam2DegPolyHMMWithMuts* create(std::vector<DepthPair*>* depths, std::vector<std::string>* sampleList, std::vector<MutationList*>* mutListVec, gsl_vector* allIndEstMutCounts, gsl_vector* fixedParams, int maxPloidy, int numPairs, int numPairsToSummarize, int ordering, unsigned int seed, gsl_matrix* stage1NearestCellIdxMat, gsl_vector* meanVarianceCoefVec, bool preallocIntermediates = true);

    virtual gsl_vector* getAllPairedEstMutCounts() override;
};

#endif

