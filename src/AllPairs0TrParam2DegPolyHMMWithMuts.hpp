#ifndef ALLPAIRS0TRPARAM2DEGPOLYHMMWITHMUTS_HPP
#define ALLPAIRS0TRPARAM2DEGPOLYHMMWITHMUTS_HPP

#include "AllPairs0TrParam2DegPolyHMM.hpp"
#include "TwoCell0TrParam2DegPolyHMMWithMuts.hpp"
#include "MutationList.hpp"
#include "MutationPair.hpp"

/*
 * This class is for making an HMM that has fixed alpha/beta/lambda values and estimates branches
 * and library scaling sizes only across pairs of tumor cells. This is similar to AllPairs3TrParam2DegPolyHMM,
 * but all the transition params are fixed.
 * This class also incorporates mutations, by TwoCell0TrParam2DegPolyHMMWithMuts (and MutationList and MutationPair internally)
 *
 * this->paramsToEst = [lib0, lib1, ..., libN, t1_cell0_1, t2_cell0_1, t3_cell0_1, t1_cell0_2, t2_cell0_2, t3_cell0_2, ..., t3_cell(N-1)_N]
 * this->fixedParams = [alpha, beta, lambda, cnaToMutRateMu]
 */
class AllPairs0TrParam2DegPolyHMMWithMuts : public AllPairs0TrParam2DegPolyHMM {
  private:
  protected:
    std::vector<MutationList*>* mutListVec;
    std::vector<MutationPair*>* mutPairVec;
    gsl_vector* allIndEstMutCounts;

    AllPairs0TrParam2DegPolyHMMWithMuts(std::vector<DepthPair*>* depths, std::vector<std::string>* sampleList, std::vector<MutationList*>* mutListVec, gsl_vector* allIndEstMutCounts, gsl_vector* fixedParams, int maxPloidy, int numFixedLibs, int numPairs, int numBranchesToEst);

    virtual void makeOneHMMPair(int i, int j, bool preallocIntermediates = true) override;
    using AllPairs0TrParam2DegPolyHMM::makeHMMPairs; // unhide parent method of same name https://stackoverflow.com/a/18100999

  public:
    // constructors and destructor
    virtual ~AllPairs0TrParam2DegPolyHMMWithMuts();
    static AllPairs0TrParam2DegPolyHMMWithMuts* create(std::vector<DepthPair*>* depths, std::vector<std::string>* sampleList, std::vector<MutationList*>* mutListVec, gsl_vector* allIndEstMutCounts, gsl_vector* fixedParams, int maxPloidy, int numPairs, gsl_vector* meanVarianceCoefVec, bool preallocIntermediates = true);

    // accessors and mutators
    virtual gsl_vector* getAllPairedEstMutCounts() override;
};

#endif

