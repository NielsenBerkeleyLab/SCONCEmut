#ifndef TWOCELL3TRPARAM2DEGPOLYHMM_HPP
#define TWOCELL3TRPARAM2DEGPOLYHMM_HPP

#include <gsl/gsl_statistics_double.h>
#include "HMM.hpp"

/*
 * HMM for 2 cells, where alpha is fixed, beta/gamma are estimated, as are lib sizes and branches
 *
 * this->paramsToEst = [lib0, lib1, beta, lambda, t1, t2, t3]
 * this->fixedParams = [alpha]
 */
class TwoCell3TrParam2DegPolyHMM : public HMM {
  private:
    // member variables
    std::vector<double>* logFacKVec; // k:[log(k!)]
    double getLogFacK(int k);
    void setLogFacK();

  protected:
    gsl_matrix* timeDepMatrixA; // matrix A for time dependent transitions along lineage, for set of times [t1,t2,t3]
    gsl_matrix* timeDepMatrixP2; // matrix P for time dependent transitions along t2 lineage. HMM::timeDepMatrixP is for shared lineage
    gsl_matrix* timeDepMatrixP3; // matrix P for time dependent transitions along t3 lineage

    TwoCell3TrParam2DegPolyHMM(std::vector<DepthPair*>* depths, gsl_vector* fixedParams, int maxPloidy, int numTrParamsToEst, int numFixedTrParams, int numFixedLibs, int numBranches, bool preallocIntermediates = true);
    double setTransition(gsl_matrix* dest, double alpha, double beta, double lambda, double t1, double t2, double t3);
    double setTimeDepMatrixA(double t1, double t2, double t3);

  public:
    // constructors and destructor
    TwoCell3TrParam2DegPolyHMM(std::vector<DepthPair*>* depths, gsl_vector* fixedParams, int maxPloidy, bool preallocIntermediates = true);
    TwoCell3TrParam2DegPolyHMM(std::vector<DepthPair*>* depths, int maxPloidy, bool preallocIntermediates = true);
    virtual ~TwoCell3TrParam2DegPolyHMM();

    // accessors and mutators
    virtual void print(FILE* stream) override;
    virtual void setLibScalingFactor(int cellNum, double libScalingFactor) override;
    virtual double getLibScalingFactor(int cellNum) const override;
    void setLogFacK(int maxDepth);

    // functions that depend on numbering and ordering of transition params
    virtual double setTransition(gsl_matrix* dest, gsl_vector* transitionParams) override;
    using HMM::setTransition;
    virtual void convertProbToParam(gsl_vector* dest, const gsl_vector* src) const override;
    virtual void convertParamToProb(gsl_vector* dest, const gsl_vector* src) const override;

    // functions that depend on model
    virtual double getEmissionProb(double tumorDepth, double diploidDepth, int ploidy, int cellIdx) override;
    virtual double getTotalEmissionProb(int stateIdx, std::vector<std::vector<double>*>* currChrDepthsVec, int chrIdx, int depthIdx) override;
    virtual double getLogEmissionProb(double tumorDepth, double diploidDepth, int ploidy, int cellIdx) override;
    virtual double getTotalLogEmissionProb(int stateIdx, std::vector<std::vector<double>*>* currChrDepthsVec, int chrIdx, int depthIdx) override;
    virtual TwoCell3TrParam2DegPolyHMM* bfgs(gsl_vector* initGuess, int maxIters, bool verbose = true, bool debug = false) override;
    virtual void simulate() override;
    virtual void simulate(int seed) override;
    virtual void simulate(int seed, bool simDiploid, double diploid_lambda_i = 272.5568, int numDiploidCells = 35);
    virtual void setInitGuessNthTime(gsl_vector* initGuess, int iter, int numTotalRuns) const override;
    virtual void setUpBaumWelchLeastSquares() override;
    virtual double baumWelchLeastSquares_f(gsl_vector* probs) override;
    virtual double checkOptimProbValidity(gsl_vector* probs) const override;
};

#endif

