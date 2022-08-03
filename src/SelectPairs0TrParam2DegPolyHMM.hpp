#ifndef SELECTPAIRS0TRPARAM2DEGPOLYHMM_HPP
#define SELECTPAIRS0TRPARAM2DEGPOLYHMM_HPP

#include "IndThenPairs2Stages3TrParam2DegPolyHMM.hpp"
#include "AllPairs0TrParam2DegPolyHMM.hpp"
#include "TwoCell0TrParam2DegPolyHMM.hpp"

/*
 * This class is for making an HMM that has fixed alpha/beta/lambda values and estimates branches
 * and library scaling sizes only across pairs of tumor cells. This is similar to AllPairs3TrParam2DegPolyHMM,
 * but all the transition params are fixed.
 *
 * Specifically, only entries that are selected to be optimized have non null values in this->hmmVec. this->hmmVec still has the full size expected (ie numPairs).
 * this->paramsToEst also has the full paramset, even if they're never used in an HMM (this will be helpful if we want to alloc HMMs on the fly, as all param values will already be stored)
 * Selected pairs are determined by the numPairsToSummarize'th nearest pairs for each cell, as determined by stage1NearestCellIdxMat.
 * The selection process is done in makeHMMPairs.
 *
 * this->paramsToEst = [lib0, lib1, ..., libN, t1_cell0_1, t2_cell0_1, t3_cell0_1, t1_cell0_2, t2_cell0_2, t3_cell0_2, ..., t3_cell(N-1)_N]
 * this->fixedParams = [alpha, beta, lambda]
 */
class SelectPairs0TrParam2DegPolyHMM : public AllPairs0TrParam2DegPolyHMM {
  private:
  protected:
    int numPairsToSummarize;
    int ordering;
    unsigned int seed;
    gsl_matrix* stage1NearestCellIdxMat;

    SelectPairs0TrParam2DegPolyHMM(std::vector<DepthPair*>* depths, std::vector<std::string>* sampleList, gsl_vector* fixedParams, int maxPloidy, int numFixedLibs, int numPairs, int numPairsToSummarize, int ordering, unsigned int seed, gsl_matrix* stage1NearestCellIdxMat, int numBranchesToEst);

    virtual void makeHMMPairs(gsl_vector* meanVarianceCoefVec, gsl_vector* transitionParams, bool preallocIntermediates = true) override;
    using AllPairs0TrParam2DegPolyHMM::makeHMMPairs; // unhide parent method of same name https://stackoverflow.com/a/18100999

  public:
    // constructors and destructor
    virtual ~SelectPairs0TrParam2DegPolyHMM();
    static SelectPairs0TrParam2DegPolyHMM* create(std::vector<DepthPair*>* depths, std::vector<std::string>* sampleList, gsl_vector* fixedParams, int maxPloidy, int numPairs, int numPairsToSummarize, int ordering, unsigned int seed, gsl_matrix* stage1NearestCellIdxMat, gsl_vector* meanVarianceCoefVec, bool preallocIntermediates = true);

    virtual void print(FILE* stream) override;
    //virtual void decodeStatSummaryOneCellAcrossPairs(int summaryMethod, int cellNum) override;
    virtual SelectPairs0TrParam2DegPolyHMM* bfgs(gsl_vector* initGuess, int maxIters, bool verbose = true, bool debug = false) override;
    void resetSkippedParams();
};

#endif


