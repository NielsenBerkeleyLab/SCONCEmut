#ifndef TWOCELLFIXLIB3TRPARAM2DEGPOLYHMM_HPP
#define TWOCELLFIXLIB3TRPARAM2DEGPOLYHMM_HPP

#include "TwoCell3TrParam2DegPolyHMM.hpp"

/*
 * This class is uesd in the second stage of AllPairs2Stages3TrParam2DegPolyHMM
 * by the AllPairs0TrParam2DegPolyHMM class.
 *
 * It's the same as TwoCell3TrParam2DegPolyHMM but all the library sizes are fixed
 *
 * this->paramsToEst = [beta, lambda, t1, t2, t3]
 * this->fixedParams = [lib0, lib1, alpha]
 */
class TwoCellFixLib3TrParam2DegPolyHMM : public TwoCell3TrParam2DegPolyHMM {
  protected:
    TwoCellFixLib3TrParam2DegPolyHMM(std::vector<DepthPair*>* depths, gsl_vector* fixedParams, int maxPloidy, int numTrParamsToEst, int numFixedTrParams, int numFixedLibs, int numBranches, bool preallocIntermediates = true);
    std::vector<gsl_matrix*>* totalLogEmissionLookup; // chrIdx:[rows: depthIdx, cols: stateIdx]
    std::vector<gsl_matrix*>* totalEmissionLookup; // chrIdx:[rows: depthIdx, cols: stateIdx]
    void updateAllTotalEmissionLookup();

  public:
    // constructors and destructor
    TwoCellFixLib3TrParam2DegPolyHMM(std::vector<DepthPair*>* depths, gsl_vector* fixedParams, int maxPloidy, bool preallocIntermediates = true);
    virtual ~TwoCellFixLib3TrParam2DegPolyHMM();

    // accessors and mutators
    virtual void setLibScalingFactor(int cellNum, double libScalingFactor) override;
    virtual double getLibScalingFactor(int cellNum) const override;
    virtual double getTotalLogEmissionProb(int stateIdx, std::vector<std::vector<double>*>* currChrDepthsVec, int chrIdx, int depthIdx) override;
    virtual double getTotalEmissionProb(int stateIdx, std::vector<std::vector<double>*>* currChrDepthsVec, int chrIdx, int depthIdx) override;
    virtual void setMeanVarianceFn(gsl_vector* meanVarianceCoefVec) override;

    // functions that depend on numbering and ordering of transition params
    virtual void convertProbToParam(gsl_vector* dest, const gsl_vector* src) const override;
    virtual void convertParamToProb(gsl_vector* dest, const gsl_vector* src) const override;

    // functions that depend on model
    virtual TwoCellFixLib3TrParam2DegPolyHMM* bfgs(gsl_vector* initGuess, int maxIters, bool verbose = true, bool debug = false) override;
};

#endif

