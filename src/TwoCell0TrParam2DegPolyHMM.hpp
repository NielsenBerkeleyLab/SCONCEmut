#ifndef TWOCELL0TRPARAM2DEGPOLYHMM_HPP
#define TWOCELL0TRPARAM2DEGPOLYHMM_HPP

#include "TwoCell3TrParam2DegPolyHMM.hpp"

/*
 * This class is uesd in the second stage of AllPairs2Stages3TrParam2DegPolyHMM
 * by the AllPairs0TrParam2DegPolyHMM class
 *
 * this->paramsToEst = [lib0, lib1, t1, t2, t3]
 * //this->fixedParams = [alpha, beta, gamma]
 * this->fixedParams = [alpha, beta, lambda]
 */
class TwoCell0TrParam2DegPolyHMM : public TwoCell3TrParam2DegPolyHMM {
  protected:
    TwoCell0TrParam2DegPolyHMM(std::vector<DepthPair*>* depths, gsl_vector* fixedParams, int maxPloidy, int numTrParamsToEst, int numFixedTrParams, int numFixedLibs, int numBranches, bool preallocIntermediates = true);

  public:
    // constructors and destructor
    TwoCell0TrParam2DegPolyHMM(std::vector<DepthPair*>* depths, gsl_vector* fixedParams, int maxPloidy, bool preallocIntermediates = true);
    //TwoCell0TrParam2DegPolyHMM(const TwoCell0TrParam2DegPolyHMM& otherHMM);
    virtual ~TwoCell0TrParam2DegPolyHMM();

    // accessors and mutators

    // functions that depend on numbering and ordering of transition params
    //virtual double setTransition(gsl_vector* transitionParams) override;
    virtual double setTransition(gsl_matrix* dest, gsl_vector* transitionParams) override;
    using HMM::setTransition;
    virtual void convertProbToParam(gsl_vector* dest, const gsl_vector* src) const override;
    virtual void convertParamToProb(gsl_vector* dest, const gsl_vector* src) const override;

    // functions that depend on model
    //virtual int getMaxNumBFGSStarts() const override;
    virtual TwoCell0TrParam2DegPolyHMM* bfgs(gsl_vector* initGuess, int maxIters, bool verbose = true, bool debug = false) override;
    //virtual void simulate() override;
    //virtual void simulate(int seed) override;
    virtual void setInitGuessNthTime(gsl_vector* initGuess, int iter, int numTotalRuns) const override;

};

#endif

