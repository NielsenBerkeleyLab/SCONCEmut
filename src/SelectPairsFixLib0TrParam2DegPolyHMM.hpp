#ifndef SELECTPAIRSFIXLIB0TRPARAM2DEGPOLYHMM_HPP
#define SELECTPAIRSFIXLIB0TRPARAM2DEGPOLYHMM_HPP

#include "IndThenPairs2Stages3TrParam2DegPolyHMM.hpp"
#include "AllPairsFixLib0TrParam2DegPolyHMM.hpp"
#include "TwoCellFixLib0TrParam2DegPolyHMM.hpp"

/*
 * This class is for making an HMM that has fixed alpha/beta/lambda and library size scaling factor values and
 * estimates branches only across selected pairs of tumor cells. This is similar to AllPairsFixLib0TrParam2DegPolyHMM,
 * but only selected pairs are computed.
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
 * this->fixedParams = [lib0, lib1, ..., libN, alpha, beta, lambda]
 */
class SelectPairsFixLib0TrParam2DegPolyHMM : public AllPairsFixLib0TrParam2DegPolyHMM {
  private:

  protected:
    int numPairsToSummarize;
    int ordering;
    unsigned int seed;
    gsl_matrix* stage1NearestCellIdxMat;

    SelectPairsFixLib0TrParam2DegPolyHMM(std::vector<DepthPair*>* depths, std::vector<std::string>* sampleList, gsl_vector* fixedParams, int maxPloidy, int numPairs, int numPairsToSummarize, int ordering, unsigned int seed, gsl_matrix* stage1NearestCellIdxMat, int numBranchesToEst);
    virtual void makeHMMPairs(gsl_vector* meanVarianceCoefVec, gsl_vector* transitionParams, bool preallocIntermediates = true) override;
    using AllPairs0TrParam2DegPolyHMM::makeHMMPairs; // unhide parent method of same name https://stackoverflow.com/a/18100999

  public:
    // constructors and destructor
    virtual ~SelectPairsFixLib0TrParam2DegPolyHMM();
    static SelectPairsFixLib0TrParam2DegPolyHMM* create(std::vector<DepthPair*>* depths, std::vector<std::string>* sampleList, gsl_vector* fixedParams, int maxPloidy, int numPairs, int numPairsToSummarize, int ordering, unsigned int seed, gsl_matrix* stage1NearestCellIdxMat, gsl_vector* meanVarianceCoefVec, bool preallocIntermediates = true);

    virtual void print(FILE* stream) override;
    virtual SelectPairsFixLib0TrParam2DegPolyHMM* bfgs(gsl_vector* initGuess, std::string filename, int maxIters, bool verbose = true, bool debug = false) override;
    void resetSkippedParams();
};

#endif

