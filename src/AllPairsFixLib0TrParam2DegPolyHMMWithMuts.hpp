#ifndef ALLPAIRSFIXLIB0TRPARAM2DEGPOLYHMMWITHMUTS_HPP
#define ALLPAIRSFIXLIB0TRPARAM2DEGPOLYHMMWITHMUTS_HPP

#include "AllPairsFixLib0TrParam2DegPolyHMM.hpp"
#include "TwoCellFixLib0TrParam2DegPolyHMMWithMuts.hpp"
#include "MutationList.hpp"
#include "MutationPair.hpp"

/*
 * This class is for making an HMM that has fixed alpha/beta/lambda and library size scaling factor values and
 * estimates branches only across pairs of tumor cells. This is similar to AllPairs0TrParam2DegPolyHMM,
 * but all the transition params and library size scaling factors are fixed.
 * This class also incorporates mutations, by TwoCellFixLib0TrParam2DegPolyHMMWithMuts (and MutationList and MutationPair internally)
 *
 * This class is intended to only be used as a second stage BFGS class. Each HMM's bfgs is called sequentially
 * (ie can be parallelized). No parameters are shared between HMMs (ie no joint estimations), as each HMM
 * estimates its own branch lengths
 *
 * this->paramsToEst = [t1_cell0_1, t2_cell0_1, t3_cell0_1, t1_cell0_2, t2_cell0_2, t3_cell0_2, ..., t3_cell(N-1)_N]
 * this->fixedParams = [lib0, lib1, ..., libN, alpha, beta, lambda, cnaToMutRateMu]
 */
class AllPairsFixLib0TrParam2DegPolyHMMWithMuts : public AllPairsFixLib0TrParam2DegPolyHMM {
  private:

  protected:
    // member variables
    std::vector<MutationList*>* mutListVec;
    std::vector<MutationPair*>* mutPairVec;
    gsl_vector* allIndEstMutCounts;

    AllPairsFixLib0TrParam2DegPolyHMMWithMuts(std::vector<DepthPair*>* depths, std::vector<std::string>* sampleList, std::vector<MutationList*>* mutListVec, gsl_vector* allIndEstMutCounts, gsl_vector* fixedParams, int maxPloidy, int numPairs, int numBranchesToEst);
    using AllPairsFixLib0TrParam2DegPolyHMM::makeHMMPairs; // unhide parent method of same name https://stackoverflow.com/a/18100999
    virtual void makeOneHMMPair(int i, int j, bool preallocIntermediates = true) override;

  public:
    // constructors and destructor
    virtual ~AllPairsFixLib0TrParam2DegPolyHMMWithMuts();
    static AllPairsFixLib0TrParam2DegPolyHMMWithMuts* create(std::vector<DepthPair*>* depths, std::vector<std::string>* sampleList, std::vector<MutationList*>* mutListVec, gsl_vector* allIndEstMutCounts, gsl_vector* fixedParams, int maxPloidy, int numPairs, gsl_vector* meanVarianceCoefVec, bool preallocIntermediates = true);

    // accessors and mutators
    virtual gsl_vector* getAllPairedEstMutCounts() override;

};

#endif

