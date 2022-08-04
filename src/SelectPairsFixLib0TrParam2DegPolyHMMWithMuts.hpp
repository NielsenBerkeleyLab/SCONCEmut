#ifndef SELECTPAIRSFIXLIB0TRPARAM2DEGPOLYHMMWITHMUTS_HPP
#define SELECTPAIRSFIXLIB0TRPARAM2DEGPOLYHMMWITHMUTS_HPP

#include "SelectPairsFixLib0TrParam2DegPolyHMM.hpp"
#include "TwoCellFixLib0TrParam2DegPolyHMMWithMuts.hpp"
#include "MutationList.hpp"
#include "MutationPair.hpp"

/*
 * This class is for making an HMM that has fixed alpha/beta/lambda and library size scaling factor values and
 * estimates branches only across selected pairs of tumor cells. This is similar to AllPairsFixLib0TrParam2DegPolyHMM,
 * but only selected pairs are computed.
 * This class also incorporates mutations, by TwoCellFixLib0TrParam2DegPolyHMMWithMuts (and MutationList and MutationPair internally)
 *
 * Specifically, only entries that are selected to be optimized have non null values in this->hmmVec. this->hmmVec still has the full size expected (ie numPairs).
 * this->paramsToEst also has the full paramset, even if they're never used in an HMM (this will be helpful if we want to alloc HMMs on the fly, as all param values will already be stored)
 * Selected pairs are determined by the numPairsToSummarize'th nearest pairs for each cell, as determined by stage1NearestCellIdxMat.
 * The selection process is done in makeHMMPairs.
 *
 * This class is intended to only be used as a second stage BFGS class. Each HMM's bfgs is called sequentially
 * (ie can be parallelized). No parameters are shared between HMMs (ie no joint estimations), as each HMM
 * estimates its own branch lengths
 *
 * this->paramsToEst = [t1_cell0_1, t2_cell0_1, t3_cell0_1, t1_cell0_2, t2_cell0_2, t3_cell0_2, ..., t3_cell(N-1)_N]
 * this->fixedParams = [lib0, lib1, ..., libN, alpha, beta, lambda, cnaToMutRateMu]
 */
class SelectPairsFixLib0TrParam2DegPolyHMMWithMuts : public SelectPairsFixLib0TrParam2DegPolyHMM {
  private:

  protected:
    // member variables
    std::vector<MutationList*>* mutListVec;
    std::vector<MutationPair*>* mutPairVec;
    gsl_vector* allIndEstMutCounts;

    SelectPairsFixLib0TrParam2DegPolyHMMWithMuts(std::vector<DepthPair*>* depths, std::vector<std::string>* sampleList, std::vector<MutationList*>* mutListVec, gsl_vector* allIndEstMutCounts, gsl_vector* fixedParams, int maxPloidy, int numPairs, int numPairsToSummarize, int ordering, unsigned int seed, gsl_matrix* stage1NearestCellIdxMat, int numBranchesToEst);
    using AllPairs0TrParam2DegPolyHMM::makeHMMPairs; // unhide parent method of same name https://stackoverflow.com/a/18100999
    virtual void makeOneHMMPair(int i, int j, bool preallocIntermediates = true) override;

  public:
    // constructors and destructor
    virtual ~SelectPairsFixLib0TrParam2DegPolyHMMWithMuts();
    static SelectPairsFixLib0TrParam2DegPolyHMMWithMuts* create(std::vector<DepthPair*>* depths, std::vector<std::string>* sampleList, std::vector<MutationList*>* mutListVec, gsl_vector* allIndEstMutCounts, gsl_vector* fixedParams, int maxPloidy, int numPairs, int numPairsToSummarize, int ordering, unsigned int seed, gsl_matrix* stage1NearestCellIdxMat, gsl_vector* meanVarianceCoefVec, bool preallocIntermediates = true);

    virtual gsl_vector* getAllPairedEstMutCounts() override;
};

#endif

