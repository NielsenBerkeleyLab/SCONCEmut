#ifndef DISJOINTPAIRS3TRPARAM2DEGPOLYHMM_HPP
#define DISJOINTPAIRS3TRPARAM2DEGPOLYHMM_HPP

#include <boost/random.hpp>

#include "Optimizable.hpp"
#include "TwoCell3TrParam2DegPolyHMM.hpp"
#include "AllPairs3TrParam2DegPolyHMM.hpp"

/*
 * This class is for making an HMM that jointly estimates alpha, beta, and lambda across disjoint pairs of
 * tumor cells, where each pair is independent of the other pairs (except for the shared parameters
 * and cell specific library size scaling factor (this is used every time a cell is used). Branch
 * lengths are estimated for each pair of cells
 *
 * This is the same as AllPairs3TrParam2DegPolyHMM, but with a different constructor so that
 * only disjoint/non overlapping pairs are formed. This is to reduce the overall number of calculations
 * needed in AllPairs2Stages3TrParam2DegPolyHMM.
 *
 * //this->paramsToEst = [lib0, lib1, ..., libN, beta, gamma, t1_cell0_1, t2_cell0_1, t3_cell0_1, t1_cell2_3, t2_cell2_3, t3_cell2_3, ..., t3_cell(N-1)_N]
 * this->paramsToEst = [lib0, lib1, ..., libN, beta, lambda, t1_cell0_1, t2_cell0_1, t3_cell0_1, t1_cell2_3, t2_cell2_3, t3_cell2_3, ..., t3_cell(N-1)_N]
 * this->fixedParams = [alpha]
 */
class DisjointPairs3TrParam2DegPolyHMM : public AllPairs3TrParam2DegPolyHMM {
  private:
    /*// member variables
    std::vector<HMM*>* hmmVec;
    std::vector<DepthPair*>* depthsVec;
    std::vector<std::string>* sampleList; // list of cell/sample names as determined by filename (corresponds to depthsVec)
    */

  protected:
    DisjointPairs3TrParam2DegPolyHMM(std::vector<DepthPair*>* depths, std::vector<std::string>* sampleList, gsl_vector* fixedParams, int maxPloidy, int numFixedLibs, int numPairs, int numBranchesToEst);
    //DisjointPairs3TrParam2DegPolyHMM(const DisjointPairs3TrParam2DegPolyHMM& otherDisjointPairs3TrParam2DegPolyHMM);

    virtual int getCell0IdxFromHMMIdx(int hmmIdx) override;
    virtual int getCell1IdxFromHMMIdx(int hmmIdx) override;
    virtual int getHMMIdxFromCellPair(int cell0Idx, int cell1Idx) override;
    virtual void makeHMMPairs(gsl_vector* meanVarianceCoefVec, gsl_vector* transitionParams, bool preallocIntermediates = true) override;
    using AllPairs3TrParam2DegPolyHMM::makeHMMPairs; // unhide parent method of same name https://stackoverflow.com/a/18100999

    /*AllPairs3TrParam2DegPolyHMM(std::vector<DepthPair*>* depths, std::vector<std::string>* sampleList, int maxPloidy, int numSharedTrParamsToEst);
    virtual HMM* makeNewHMM(std::vector<DepthPair*>* depths, int maxPloidy) const;
    virtual HMM* copyHMM(HMM* hmm) const;
    */

  public:
    /*// constants
    const int LIB_SIZE_SCALING_FACTOR_START_IDX;
    const int NUM_CELLS;
    const int NUM_PAIRS;
    const int NUM_SHARED_TRANSITION_PARAMS_TO_EST; // ie. alpha, beta, gamma. Does not count branch lengths or lib size scaling factors
    const int SHARED_TRANSITION_PROB_START_IDX;
    const int BRANCH_LENGTH_START_IDX;
    const int NUM_BRANCH_LENGTHS_TO_EST;
    const int MAX_PLOIDY;
    */

    // constructors and destructor
    //AllPairs3TrParam2DegPolyHMM();
    virtual ~DisjointPairs3TrParam2DegPolyHMM();
    static DisjointPairs3TrParam2DegPolyHMM* create(std::vector<DepthPair*>* depths, std::vector<std::string>* sampleList, int maxPloidy, int numPairs, gsl_vector* meanVarianceCoefVec, bool preallocIntermediates = true);
    static DisjointPairs3TrParam2DegPolyHMM* create(const DisjointPairs3TrParam2DegPolyHMM& otherDisjointPairs3TrParam2DegPolyHMM);

    /*// accessors and mutators
    int getKploidy() const;
    //std::vector<TwoCell3TrParam2DegPolyHMM*>* getHMMs();
    std::vector<HMM*>* getHMMs();
    virtual void print(FILE* stream) override;
    std::vector<std::string>* getSampleList();
    void setSampleList(std::vector<std::string>* sampleList);
    //virtual double getLibScalingFactor(int cellNum) const override;
    //virtual double setTransition() override;

    // override Optimizable methods
    virtual int getMaxNumBFGSStarts() const override;
    virtual double setParamsToEst(gsl_vector* params) override;
    virtual void convertProbToParam(gsl_vector* dest, const gsl_vector* src) const override;
    virtual void convertParamToProb(gsl_vector* dest, const gsl_vector* src) const override;
    virtual double getLogLikelihood() override;
    virtual double checkOptimProbValidity(gsl_vector* probs) override;

    // BFGS
    //AllPairs3TrParam2DegPolyHMM* callBFGSNTimes(int numRuns, bool verbose = true, int seed = 43);
    virtual void setInitGuessNthTime(gsl_vector* initGuess, int iter, int numTotalRuns) const override;*/
    DisjointPairs3TrParam2DegPolyHMM* bfgs(gsl_vector* initGuess, int maxIters, bool verbose = true, bool debug = false) override;
    /*double runForwardAlg();

    // methods to save results
    void viterbiDecodeAll();
    void saveAllViterbiDecodedCNA(std::string filename);
    //// functions that depend on numbering and ordering of transition params
    //virtual double setTransition(gsl_vector* transitionParams) override;

    //// functions that depend on model
    //virtual double getEmissionProb(double tumorDepth, double diploidDepth, double ploidy, int cellIdx) override;
    //virtual void simulate() override;
    //virtual void simulate(int seed) override;*/

};

#endif

