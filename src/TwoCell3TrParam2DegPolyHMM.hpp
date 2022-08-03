#ifndef TWOCELL3TRPARAM2DEGPOLYHMM_HPP
#define TWOCELL3TRPARAM2DEGPOLYHMM_HPP

#include <gsl/gsl_statistics_double.h>
#include "HMM.hpp"

/*
 * HMM for 2 cells, where alpha is fixed, beta/gamma are estimated, as are lib sizes and branches
 *
 * //this->paramsToEst = [lib0, lib1, beta, gamma, t1, t2, t3]
 * this->paramsToEst = [lib0, lib1, beta, lambda, t1, t2, t3]
 * this->fixedParams = [alpha]
 */
class TwoCell3TrParam2DegPolyHMM : public HMM {
  private:
    // member variables
    std::vector<double>* logFacKVec; // k:[log(k!)]
    double getLogFacK(int k);
    void setLogFacK();

    // variables and fns that have to do with least squares
    //gsl_matrix* coefs = nullptr;
    //gsl_vector* residuals = nullptr;
    //gsl_vector* a_ij = nullptr; // residuals from baum welch, not very useful
    //gsl_matrix* baumWelchTransitionMat = nullptr;
    //static double baumWelchLeastSquares_f(const gsl_vector* v, void* params);
    //static void baumWelchLeastSquares_df(const gsl_vector* v, void* params, gsl_vector* df);
    //static void baumWelchLeastSquares_fdf(const gsl_vector* v, void* params, double* f, gsl_vector* df);
    //void baumWelchLeastSquares_convertProbToParam(gsl_vector* dest, const gsl_vector* src) const;
    //void baumWelchLeastSquares_convertParamToProb(gsl_vector* dest, const gsl_vector* src) const;

    //std::vector<int>* initGuessCases;

  protected:
    //double alpha = 0.1; // fixed transition param, added Wed 18 Dec 2019 05:50:53 PM PST. should prob make setable to match AllPairs3TrParam2DegPolyHMM
    gsl_matrix* timeDepMatrixA; // matrix A for time dependent transitions along lineage, for set of times [t1,t2,t3]
    gsl_matrix* timeDepMatrixP2; // matrix P for time dependent transitions along t2 lineage. HMM::timeDepMatrixP is for shared lineage
    gsl_matrix* timeDepMatrixP3; // matrix P for time dependent transitions along t3 lineage

    TwoCell3TrParam2DegPolyHMM(std::vector<DepthPair*>* depths, gsl_vector* fixedParams, int maxPloidy, int numTrParamsToEst, int numFixedTrParams, int numFixedLibs, int numBranches, bool preallocIntermediates = true);
    //double setTransition(gsl_matrix* dest, double alpha, double beta, double gamma, double t1, double t2, double t3);
    double setTransition(gsl_matrix* dest, double alpha, double beta, double lambda, double t1, double t2, double t3);
    double setTimeDepMatrixA(double t1, double t2, double t3);

  public:
    //boost::minstd_rand* stateRNG;
    //boost::minstd_rand* diploidRNG;
    //boost::minstd_rand* tumor0RNG;
    //boost::minstd_rand* tumor1RNG;

    // constructors and destructor
    TwoCell3TrParam2DegPolyHMM(std::vector<DepthPair*>* depths, gsl_vector* fixedParams, int maxPloidy, bool preallocIntermediates = true);
    TwoCell3TrParam2DegPolyHMM(std::vector<DepthPair*>* depths, int maxPloidy, bool preallocIntermediates = true);
    //TwoCell3TrParam2DegPolyHMM(const TwoCell3TrParam2DegPolyHMM& otherHMM);
    virtual ~TwoCell3TrParam2DegPolyHMM();

    // accessors and mutators
    virtual void print(FILE* stream) override;
    virtual void setLibScalingFactor(int cellNum, double libScalingFactor) override;
    virtual double getLibScalingFactor(int cellNum) const override;
    //virtual double getAlpha() const;
    //virtual void setAlpha(double alpha);
    void setLogFacK(int maxDepth); // TODO debugging; should be private eventually

    // functions that depend on numbering and ordering of transition params
    virtual double setTransition(gsl_matrix* dest, gsl_vector* transitionParams) override;
    using HMM::setTransition;
    virtual void convertProbToParam(gsl_vector* dest, const gsl_vector* src) const override;
    virtual void convertParamToProb(gsl_vector* dest, const gsl_vector* src) const override;

    // functions that depend on model
    //virtual int getMaxNumBFGSStarts() const override;
    //virtual double getEmissionProb(double tumorDepth, double diploidDepth, int ploidy, int windowIdx, int cellIdx) override;
    virtual double getEmissionProb(double tumorDepth, double diploidDepth, int ploidy, int cellIdx) override;
    virtual double getTotalEmissionProb(int stateIdx, std::vector<std::vector<double>*>* currChrDepthsVec, int chrIdx, int depthIdx) override;
    virtual double getLogEmissionProb(double tumorDepth, double diploidDepth, int ploidy, int cellIdx) override;
    virtual double getTotalLogEmissionProb(int stateIdx, std::vector<std::vector<double>*>* currChrDepthsVec, int chrIdx, int depthIdx) override;
    virtual TwoCell3TrParam2DegPolyHMM* bfgs(gsl_vector* initGuess, int maxIters, bool verbose = true, bool debug = false) override;
    virtual void simulate() override;
    virtual void simulate(int seed) override;
    //virtual void simulate(int seed, bool simDiploid, double diploid_lambda_i = 65.37034, double diploid_var = 144.3447);
    //virtual void simulate(int seed, bool simDiploid, double diploid_lambda_i = 65.37034, int numDiploidCells = 35);
    virtual void simulate(int seed, bool simDiploid, double diploid_lambda_i = 272.5568, int numDiploidCells = 35);
    virtual void setInitGuessNthTime(gsl_vector* initGuess, int iter, int numTotalRuns) const override;
    //virtual void solveBaumWelchEstimates() override;
    virtual void setUpBaumWelchLeastSquares() override;
    virtual double baumWelchLeastSquares_f(gsl_vector* probs) override;
    //virtual void saveParamEstimates(std::string filename) const override;

    virtual double checkOptimProbValidity(gsl_vector* probs) const override;

};

#endif

