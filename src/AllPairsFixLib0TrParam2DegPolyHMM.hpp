#ifndef ALLPAIRSFIXLIB0TRPARAM2DEGPOLYHMM_HPP
#define ALLPAIRSFIXLIB0TRPARAM2DEGPOLYHMM_HPP

#include <boost/asio/thread_pool.hpp>
#include <boost/asio.hpp>
#include "AllPairs0TrParam2DegPolyHMM.hpp"
#include "TwoCellFixLib0TrParam2DegPolyHMM.hpp"

/*
 * This class is for making an HMM that has fixed alpha/beta/lambda and library size scaling factor values and
 * estimates branches only across pairs of tumor cells. This is similar to AllPairs0TrParam2DegPolyHMM,
 * but all the transition params and library size scaling factors are fixed.
 *
 * This class is intended to only be used as a second stage BFGS class. Each HMM's bfgs is called sequentially
 * (ie can be parallelized). No parameters are shared between HMMs (ie no joint estimations), as each HMM
 * estimates its own branch lengths
 *
 * this->paramsToEst = [t1_cell0_1, t2_cell0_1, t3_cell0_1, t1_cell0_2, t2_cell0_2, t3_cell0_2, ..., t3_cell(N-1)_N]
 * this->fixedParams = [lib0, lib1, ..., libN, alpha, beta, lambda]
 */
class AllPairsFixLib0TrParam2DegPolyHMM : public AllPairs0TrParam2DegPolyHMM {
  private:

  protected:
    std::vector<bool>* shouldCallBFGSOnHmmIdx;

    AllPairsFixLib0TrParam2DegPolyHMM(std::vector<DepthPair*>* depths, std::vector<std::string>* sampleList, gsl_vector* fixedParams, int maxPloidy, int numPairs, int numBranchesToEst);
    virtual void makeHMMPairs(gsl_vector* meanVarianceCoefVec, gsl_vector* transitionParams, bool preallocIntermediates = true) override;
    using AllPairs0TrParam2DegPolyHMM::makeHMMPairs; // unhide parent method of same name https://stackoverflow.com/a/18100999
    virtual void makeOneHMMPair(int i, int j, bool preallocIntermediates = true) override;
    virtual void callIndvHMMBFGS(int hmmIdx, gsl_vector* initGuess, std::string filename, int maxIters, bool verbose, bool debug);

  public:
    // constructors and destructor
    virtual ~AllPairsFixLib0TrParam2DegPolyHMM();
    static AllPairsFixLib0TrParam2DegPolyHMM* create(std::vector<DepthPair*>* depths, std::vector<std::string>* sampleList, gsl_vector* fixedParams, int maxPloidy, int numPairs, gsl_vector* meanVarianceCoefVec, bool preallocIntermediates = true);
    static AllPairsFixLib0TrParam2DegPolyHMM* create(const AllPairsFixLib0TrParam2DegPolyHMM& otherAllPairsHMM);

    // accessors and mutators
    virtual void setLibScalingFactors(int cell0Idx, int cell1Idx, double lib0, double lib1);
    virtual double getLibScalingFactor(int cellNum) const override;
    virtual void setAllLibScalingFactors(int cellNumInPair, double libScalingFactor) override;
    virtual void setParamsToEstFromIthHMM(gsl_vector* params, int hmmIdx) override;
    virtual void setFixedParamsFromIthHMM(gsl_vector* params, int hmmIdx) override;

    // override Optimizable methods
    virtual double setParamsToEst(gsl_vector* params) override;
    virtual AllPairsFixLib0TrParam2DegPolyHMM* bfgs(gsl_vector* initGuess, std::string filename, int maxIters, bool verbose = true, bool debug = false) override;
};

#endif

