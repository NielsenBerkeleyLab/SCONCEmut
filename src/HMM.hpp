#ifndef HMM_HPP
#define HMM_HPP

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_permutation.h>
#include <boost/math/distributions.hpp>

#include <iostream>
#include <string>
#include <numeric>
#include <algorithm>
#include <chrono>

#include <set>
#include <vector>
#include <forward_list>

#include "Optimizable.hpp"
#include "DepthPair.hpp"
#include "util.hpp"

/*
 * This class is the main abstract class for analyzing one pair of cells or one cell for CNVs.
 * Parameters are
 *   - 3 transition parameters (alpha=P(adjacent CNA); beta=P(any CNA); gamma=P(back to diploid))
 *   - a library size scaling factor for each cell
 *   - branch lengths (3 branches between 2 cells, 1 scaling "branch" for 1 cell, or split time between 2 cells (given total branch lengths))
 * These parameters are stored in this->paramsToEst (parameters estimated by one of GSL's optimization functions) or this->fixedParams,
 * depending on the subclass.
 *
 * Other things that can be set:
 *   - mean/variance relationship for the negative binomial relationship
 */
class HMM : public Optimizable {
  protected:
    // member variables
    std::vector<DepthPair*>* depths;
    std::vector<std::string>* states; // ploidy, from 0[,0] to maxPloidy+1[,maxPloidy+1]
    int maxObservedDepth;
    gsl_matrix* transition;
    gsl_vector* initProb;
    gsl_vector* meanVarianceCoefVec; // intercept, slope (first order), second order coef, ...
    int maxNumBFGSStarts;

    std::vector<std::unordered_map<std::string, std::vector<int>*>*>* chrToViterbiPathMapVec; // cellIdx:[chr:vector of viterbi decoded paths]; assigned after viterbiDecode() is called in the HMM
    std::vector<std::unordered_map<std::string, gsl_matrix*>*>* forBackMargMatMapVec; // cellIdx:[chr:marginalized likelihood matrix from forBackAlg]

    // for continuous time
    gsl_matrix* rateMatrixQ; // matrix Q for mutations in time along lineage
    gsl_matrix* timeDepMatrixP; // matrix P for time dependent transitions along [shared] lineage
    std::vector<std::vector<std::string*>*>* rateMatrixStrings;

    // intermediates
    std::vector<gsl_matrix*>* forwardMatVec;
    std::vector<gsl_matrix*>* backwardMatVec;
    std::vector<gsl_matrix*>* forBackMatVec;
    std::vector<gsl_matrix*>* backTraceVec;
    std::vector<gsl_vector*>* scalingVecVec;
    gsl_vector* prevForwardCol;
    gsl_vector* currForwardCol;
    gsl_vector* nextBackwardCol;
    gsl_vector* currBackwardCol;
    gsl_matrix* transitionTranspose;

    // used for eigenvalue/eigenvector calculations for finding steady state dist
    gsl_vector_complex* ssEval;
    gsl_matrix_complex* ssEvec;
    gsl_eigen_nonsymmv_workspace* ssEigenWorkspace;

    // used for eigendecomposition for matrixExponential of continuous time rate matrix. These are constantly overwritten, no guarantees are made about their contents
    gsl_vector_complex* rateEval;
    gsl_matrix_complex* rateEvec;
    gsl_eigen_nonsymmv_workspace* rateEigenWorkspace;
    gsl_matrix* rateRealEvecMat;
    gsl_matrix* rateDiagMat;
    gsl_matrix* rateLUdecompMat;
    gsl_permutation* ratePerm;
    gsl_matrix* rateInverseMat;

    // constructors (protected so only derived classes can call)
    HMM(std::vector<DepthPair*>* depths, gsl_vector* fixedParams, int numTransitionParamsToEst, int numCells, int numBranches, int maxPloidy, int numFixedTrParams, int numFixedLibs, bool preallocIntermediates = true);

    // protected helper methods
    int getRandStateIdx(double p, int fromStateIdx) const; // helper method to get a state idx according to prob p
    int getRandStateIdx(double p, gsl_vector* probVec) const;
    double setTimeDepMatrixP(gsl_matrix* destMat, double time);

    // for baum welch + least squares step
    gsl_matrix* baumWelchTransitionMat = nullptr;

  public:
    // constants
    // paramsToEst has up to 3 sections: library size scaling factors, transition probabilities, and branch lengths. Each section must be contiguous within itself
    // ie: [lib1, lib2, ..., libN, trProb1, trProb2, trProb3, ..., trProbM, t1, ..., t(NUM_BRANCH_LENGTHS_TO_EST)]
    const int MAX_PLOIDY;
    const int NUM_CELLS;
    const int NUM_LIBS_TO_EST;
    const int NUM_TRANSITION_PARAMS_TO_EST; // ex alpha, beta, gamma. Does not include branch lengths
    const int NUM_BRANCH_LENGTHS_TO_EST;
    const int NUM_FIXED_LIBS;
    const int NUM_FIXED_TRANSITION_PARAMS;
    const int LIB_SIZE_SCALING_FACTOR_START_IDX;
    const int TRANSITION_PROB_START_IDX;
    const int BRANCH_LENGTH_START_IDX;
    const int FIXED_TRANSITION_PROB_START_IDX;
    const int FIXED_MUTATION_PARAM_START_IDX; // only used for classes *WithMuts. these params must be passed in fixedParams

    static const int DEPTH_ERROR_SCALING = 50; // for emission prob. lambda_i = ___ + errorTerm, where errorTerm = diploidDepth_i / DEPTH_ERROR_SCALING. Super hacky to set here, should prob move to an external file of constants
    //static const int DEPTH_ERROR_SCALING = 100; // for emission prob. lambda_i = ___ + errorTerm, where errorTerm = diploidDepth_i / DEPTH_ERROR_SCALING. Super hacky to set here, should prob move to an external file of constants

    bool hasBeenReadFromFile = false;

    // accessors and mutators
    void allocMinimalIntermediates();
    void allocBWIntermediates();
    void allocDecodingIntermediates();
    void freeToMinimum();
    std::vector<DepthPair*>* getDepths();
    void setStates(std::vector<std::string>* states);
    std::vector<std::string>* getStates() const;
    void setTransition(gsl_matrix* transition);
    virtual double setTransition();
    void setRateMatrixQ(double alpha, double beta, double lambda);
    double setTimeDepMatrixP(double time);
    gsl_matrix* getTransition() const;
    gsl_matrix* getBaumWelchTransitionMat() const;
    int getKploidy() const;
    void setInitProb(gsl_vector* initProb);
    gsl_vector* getInitProb() const;
    virtual void setMeanVarianceFn(gsl_vector* meanVarianceCoefVec);
    static gsl_vector* createMeanVarianceCoefVec();
    double getMeanVarianceIntercept() const;
    double getMeanVarianceSlope() const;
    double getMeanVariancePoly2() const;
    gsl_vector* getMeanVarianceCoefVec() const;
    std::vector<std::unordered_map<std::string, std::vector<int>*>*>* getChrToViterbiPathMapVec() const;
    std::vector<std::unordered_map<std::string, gsl_matrix*>*>* getForBackMargMatMapVec() const;
    virtual double setParamsToEst(gsl_vector* params) override;
    double setFixedParams(gsl_vector* fixedParams);
    virtual void setSimParamsToEst(gsl_vector* params) override;
    virtual void setSimFixedParams(gsl_vector* params) override;
    std::vector<std::string>* getChrVec() const; // list of all chromosomes from DepthPair
    double setInitProbSteadyState();
    virtual void setLibScalingFactorsToTotalRatio();
    double calcLibScalingFactorsToTotalRatio(int cellIdx) const;
    virtual void setLibScalingFactors(gsl_vector* libScalingFactors);
    virtual void setLibScalingFactor(int cellNum, double libScalingFactor) = 0;
    virtual double getLibScalingFactor(int cellNum) const = 0;
    virtual double getAlpha() const;
    virtual void setAlpha(double alpha);
    virtual void estLibScalingFactorsPosterior();
    virtual double estLibScalingFactorsPosterior(int cellNum);
    int getCellPloidyFromStateIdx(int cellIdx, int stateIdx) const; // convert cellIdx and stateIdx to a ploidy. Works for any number of cells
    int getIndvPloidyFromStateIdx(int idxInPair, int stateIdx) const; // given stateIdx that corresponds *to a pair only*, get the idxInPair's ploidy
    int getStateIdxFromPloidyPair(int ploidy0, int ploidy1) const; // convert a pair of ploidies into a state idx
    std::unordered_map<std::string, std::vector<int>*>* getIthChrToViterbiPathMapVec(int cellIdx) const; // get specified chr:vector of viterbi decoded paths

    // functions that depend on numbering and ordering of transition params and change based on model
    virtual double setTransition(gsl_matrix* dest, gsl_vector* transitionParams) = 0;
    virtual double setTransition(gsl_vector* transitionParams); // always saves into this->transition

    // functions that change based on model
    // note: the log emission functions are the fastest, as the non log ones call the log ones and exp them
    virtual double getEmissionProb(double tumorDepth, double diploidDepth, int ploidy, int cellIdx) = 0;
    virtual double getTotalEmissionProb(int stateIdx, std::vector<std::vector<double>*>* currChrDepthsVec, int chrIdx, int depthIdx) = 0;
    virtual double getLogEmissionProb(double tumorDepth, double diploidDepth, int ploidy, int cellIdx) = 0;
    virtual double getTotalLogEmissionProb(int stateIdx, std::vector<std::vector<double>*>* currChrDepthsVec, int chrIdx, int depthIdx) = 0;
    virtual void simulate() = 0;
    virtual void simulate(int seed) = 0;
    virtual void setUpBaumWelchLeastSquares() = 0;
    virtual double baumWelchLeastSquares_f(gsl_vector* probs) = 0;

    // shared functions
    virtual double getLogLikelihood() override;
    virtual double getViterbiLogLikelihood();
    double runForwardAlg();
    double runBackwardAlg();
    int runForBackAlg();
    void runBaumWelch(int numBWIters = 20, bool verbose = true, bool debug = false);
    void viterbiDecode();
    gsl_vector* getAveragePloidy();
    void calcMargLikelihoods();
    double findSteadyStateDist(gsl_vector* steadyStateVec) const;
    virtual void print(FILE* stream) override;
    void saveViterbiDecodedCNA(std::string filename);
    void saveCNAToBed(std::string filename, int cellIdx);
    virtual void saveHMMToFile(std::string filename); // internally calls print
    virtual void saveParamEstimates(std::string filename) const;
    virtual double getParamsToEstFromFile(std::string filename, int numExpectedLinesPerFile);
    virtual double checkStateValidity(double epsilon = 1e-8) const override; // should check if transition matrices are ok
    virtual double checkStateValidity(gsl_matrix* mat, double epsilon = 1e-8) const; // should check if passed matrix is ok
    virtual double checkForTransientStates(); // check if states only show up transiently
    virtual double checkForInitProbGaps(); // check if initProb has big gaps in between states

    // destructor
    virtual ~HMM();

    virtual void miscFunctions() override;

};

#endif

