#ifndef TWOCELLFIXLIB0TRPARAM2DEGPOLYHMM_HPP
#define TWOCELLFIXLIB0TRPARAM2DEGPOLYHMM_HPP

#include "TwoCell0TrParam2DegPolyHMM.hpp"

/*
 * This class is uesd in the second stage of AllPairs2Stages3TrParam2DegPolyHMM
 * by the AllPairs0TrParam2DegPolyHMM class, when library sizes are fixed
 * alpha, beta, gamma are also all fixed
 *
 * this->paramsToEst = [t1, t2, t3]
 * //this->fixedParams = [lib0, lib1, alpha, beta, gamma]
 * this->fixedParams = [lib0, lib1, alpha, beta, lambda]
 */
class TwoCellFixLib0TrParam2DegPolyHMM : public TwoCell0TrParam2DegPolyHMM {
  protected:
    std::vector<gsl_matrix*>* totalLogEmissionLookup; // chrIdx:[rows: depthIdx, cols: stateIdx]. Copied directly from TwoCellFixLib3TrParam2DegPolyHMM
    std::vector<gsl_matrix*>* totalEmissionLookup; // chrIdx:[rows: depthIdx, cols: stateIdx]. Copied directly from TwoCellFixLib3TrParam2DegPolyHMM
    void updateAllTotalEmissionLookup();

    TwoCellFixLib0TrParam2DegPolyHMM(std::vector<DepthPair*>* depths, gsl_vector* fixedParams, int maxPloidy, int numTrParamsToEst, int numFixedTrParams, int numFixedLibs, int numBranches, bool preallocIntermediates = true);

  public:
    // constructors and destructor
    TwoCellFixLib0TrParam2DegPolyHMM(std::vector<DepthPair*>* depths, gsl_vector* fixedParams, int maxPloidy, bool preallocIntermediates = true);
    //TwoCellFixLib0TrParam2DegPolyHMM(const TwoCellFixLib0TrParam2DegPolyHMM& otherHMM);
    virtual ~TwoCellFixLib0TrParam2DegPolyHMM();

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
    virtual TwoCellFixLib0TrParam2DegPolyHMM* bfgs(gsl_vector* initGuess, int maxIters, bool verbose = true, bool debug = false) override;
    //virtual void simulate() override;
    //virtual void simulate(int seed) override;

    virtual void miscFunctions() override;

    //virtual double setParamsToEst(gsl_vector* params) override;
};

#endif

