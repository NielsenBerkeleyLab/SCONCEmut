#ifndef ALLPAIRS0TRPARAM2DEGPOLYHMM_HPP
#define ALLPAIRS0TRPARAM2DEGPOLYHMM_HPP

#include "AllPairs3TrParam2DegPolyHMM.hpp"
#include "TwoCell0TrParam2DegPolyHMM.hpp"

/*
 * This class is for making an HMM that has fixed alpha/beta/lambda values and estimates branches
 * and library scaling sizes only across pairs of tumor cells. This is similar to AllPairs3TrParam2DegPolyHMM,
 * but all the transition params are fixed.
 *
 * this->paramsToEst = [lib0, lib1, ..., libN, t1_cell0_1, t2_cell0_1, t3_cell0_1, t1_cell0_2, t2_cell0_2, t3_cell0_2, ..., t3_cell(N-1)_N]
 * this->fixedParams = [alpha, beta, lambda]
 */
class AllPairs0TrParam2DegPolyHMM : public AllPairs3TrParam2DegPolyHMM {
  private:
  protected:
    AllPairs0TrParam2DegPolyHMM(std::vector<DepthPair*>* depths, std::vector<std::string>* sampleList, gsl_vector* fixedParams, int maxPloidy, int numFixedLibs, int numPairs, int numBranchesToEst);

    virtual void makeHMMPairs(gsl_vector* meanVarianceCoefVec, bool preallocIntermediates = true) override;
    virtual void makeHMMPairs(gsl_vector* meanVarianceCoefVec, gsl_vector* transitionParams, bool preallocIntermediates = true) override;
    virtual void makeOneHMMPair(int i, int j, bool preallocIntermediates = true);
    using AllPairs3TrParam2DegPolyHMM::makeHMMPairs; // unhide parent method of same name https://stackoverflow.com/a/18100999

  public:
    // constructors and destructor
    virtual ~AllPairs0TrParam2DegPolyHMM();
    static AllPairs0TrParam2DegPolyHMM* create(std::vector<DepthPair*>* depths, std::vector<std::string>* sampleList, gsl_vector* fixedParams, int maxPloidy, int numPairs, gsl_vector* meanVarianceCoefVec, bool preallocIntermediates = true);
    static AllPairs0TrParam2DegPolyHMM* create(const AllPairs0TrParam2DegPolyHMM& otherAllPairsHMM);

    // accessors and mutators
    virtual void setFixedParams(gsl_vector* params);
    gsl_vector* getFixedParams();
    void setFixedTrParams(gsl_vector* params);
    double setAllBranches(gsl_vector* branches);
    virtual void setParamsToEstFromIthHMM(gsl_vector* params, int hmmIdx) override;
    virtual void setFixedParamsFromIthHMM(gsl_vector* params, int hmmIdx) override;
    virtual gsl_vector* getAllPairedEstMutCounts();

    // override Optimizable methods
    virtual double setParamsToEst(gsl_vector* params) override;
    virtual void convertProbToParam(gsl_vector* dest, const gsl_vector* src) const override;
    virtual void convertParamToProb(gsl_vector* dest, const gsl_vector* src) const override;
    virtual double checkOptimProbValidity(gsl_vector* probs) const override;

    // BFGS
    virtual void setInitGuessNthTime(gsl_vector* initGuess, int iter, int numTotalRuns) const override;
    virtual AllPairs0TrParam2DegPolyHMM* bfgs(gsl_vector* initGuess, int maxIters, bool verbose = true, bool debug = false) override;
    virtual AllPairs0TrParam2DegPolyHMM* bfgs(gsl_vector* initGuess, std::string filename, int maxIters, bool verbose = true, bool debug = false);

};

#endif

