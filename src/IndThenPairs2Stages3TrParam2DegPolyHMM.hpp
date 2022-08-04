#ifndef INDTHENPAIRS2STAGES3TRPARAM2DEGPOLYHMM_HPP
#define INDTHENPAIRS2STAGES3TRPARAM2DEGPOLYHMM_HPP

#include <boost/random.hpp>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_sort_vector.h>

#include "Optimizable.hpp"
#include "DepthPair.hpp"

#include "AllCells3TrParam2DegPolyHMM.hpp"

#include "AllInd2Stages3TrParam2DegPolyHMM.hpp"
#include "AllPairs0TrParam2DegPolyHMM.hpp"
#include "AllPairsFixLib0TrParam2DegPolyHMM.hpp"

#include "SelectPairs0TrParam2DegPolyHMM.hpp"
#include "SelectPairsFixLib0TrParam2DegPolyHMM.hpp"


/*
 * This class is for:
 * SCONCE as stage 0/1 to est [lib, beta, lambda, t] for each cell independently, including bw+ls and bfgs
 *
 * Because all parameters are relative, alpha is fixed in all cases.
 *
 * Then, stage 2 is all pairs of cells. libs = as is, beta = median(beta_A, beta_B,...), lambda = median(lambda_A, lambda_B,...).
 * BFGS is used to estimate [t1, t2, t3]_AB, [t1, t2, t3]_AC,..., where l(t1,t2,t3) = l_CNA(t1,t2,t3)
 * Starting points for t's, for pair AB, is:
 *   t1 = min(t_A, t_B) / 2
 *   t2 = t_A - t1
 *   t3 = t_B - t1
 * Then should get good ests for these params, can output some sort of summary CNA calls
 *
 * Note this class doesn't store any parameter values directly
 */
class IndThenPairs2Stages3TrParam2DegPolyHMM {
  protected:
    // member variables
    std::vector<AllInd2Stages3TrParam2DegPolyHMM*>* ind2StagesVec;
    AllPairs0TrParam2DegPolyHMM* stage2Pairs;
    std::vector<DepthPair*>* depthsVec;
    std::vector<std::string>* sampleList; // list of cell/sample names as determined by filename (corresponds to depthsVec)
    bool runIndBFGS;
    bool verbose;
    bool gradientDebug;
    bool centralDiff;

    gsl_matrix* stage1PairwiseDistMat;
    gsl_matrix* stage1NearestCellIdxMat;
    double calcPairwiseDist(int cell0Idx, int cell1Idx);

    gsl_vector* getIndOptimAllLibEsts() const;
    gsl_vector* getIndOptimIthTrParams(int trParamIdx) const;
    gsl_vector* getIndOptimSummaryTrParams() const;
    gsl_vector* getIndOptimAllBranchEsts() const;
    std::vector<std::unordered_map<std::string, std::vector<int>*>*>* getIndOptimAllChrToViterbiPathMapVec() const;

    // parallelizing
    void callIndvBaumWelch(int hmmIdx, gsl_vector* initGuess, int numBWIters, int numLibStarts, double libStartVal, bool verbose, bool debug); // helper method for parallelizing
    void callIndvBFGS(int hmmIdx, gsl_vector* initGuess, std::string filename, int maxIters, bool verbose, bool debug); // helper method for parallelizing

  public:
    // constants
    const int NUM_CELLS;
    const int NUM_SHARED_TRANSITION_PARAMS_TO_EST; // ie. beta, lambda. Does not count branch lengths or lib size scaling factors
    const int MAX_PLOIDY;
    const int NUM_PAIRS_STAGE_2;

    static const int ORDERING_NEAREST  = 0;
    static const int ORDERING_RANDOM   = 1;
    static const int ORDERING_FURTHEST = 2;

    // parallelizing
    int numThreads;

    // constructors and destructor
    IndThenPairs2Stages3TrParam2DegPolyHMM(std::vector<DepthPair*>* depths, std::vector<std::string>* sampleList, int maxPloidy, bool runIndBFGS, int numPairsStage2, bool verbose, bool gradientDebug, bool centralDiff);
    virtual ~IndThenPairs2Stages3TrParam2DegPolyHMM();
    void cleanUpIndOptim();

    // accessors and mutators
    std::vector<AllInd2Stages3TrParam2DegPolyHMM*>* getInd2StagesVec();
    AllPairs0TrParam2DegPolyHMM* getStage2HMM();
    virtual void print(FILE* stream);
    std::vector<std::string>* getSampleList();
    void setSampleList(std::vector<std::string>* sampleList);
    int getIndInitGuessSize() const;
    int getPairsInitGuessSize() const;
    gsl_vector* getIndOptimParamsToEstSummary() const;
    void setAllIndMeanVarianceFn(gsl_vector* meanVarianceCoefVec);
    void setAllPairsMeanVarianceFn(gsl_vector* meanVarianceCoefVec);
    void setNumThreads(int newNumThreads);
    int getNumThreads() const;

    double getTotalIndLogLikelihood();
    double getStage2LogLikelihood();

    // optimization methods
    virtual void setUpInd2Stages();
    virtual void copyStage1EstsIntoStage2FixedParams(gsl_vector* stage2FixedParams, gsl_vector* stage2InitGuess, gsl_vector* indEstParams, int numStage2LibsToFix);
    virtual void setUpAllPairsBFGS(gsl_vector* fixedParams, gsl_vector* initGuess, gsl_vector* meanVarianceCoefVec, bool fixLib = true, bool preallocIntermediates = false);
    virtual void setUpSelectPairsBFGS(gsl_vector* fixedParams, gsl_vector* initGuess, int numPairsToSummarize, int ordering, int seed, gsl_vector* meanVarianceCoefVec, bool fixLib = true, bool preallocIntermediates = false);
    virtual void optimIndCells(gsl_vector* initGuess, std::string filename, int numBWIters, int numLibStarts, int libStartVal, int maxBFGSIters, bool verbose = true, bool debug = false);
    virtual AllPairs0TrParam2DegPolyHMM* bfgsStage2(gsl_vector* initGuess, int maxIters, std::string filename, bool verbose = true, bool debug = false);

    // methods to compare and order stage 1 profiles
    void calcStage1PairwiseDistMat();
    void calcStage1NearestCellIdxMat();

    // methods to save results
    void viterbiDecodeStage1();
    void viterbiDecodeStage2();
    void saveStage1ViterbiDecodedCNA(std::string filename);
    void saveStage2ViterbiDecodedCNA(std::string filename);
    void saveStage1CNAToBed(std::string filename);
    void saveStage2CNAToBed(std::string filename);
    void saveStage1ParamEstimates(std::string filename);
    void summaryDecodeAllCellsAcrossPairs(int summaryMethod, bool runDecoding);
    void saveAllSummaryDecodedCNAToBed(std::string filename);

    // methods to read intermediate results
    void getIndOptimParamsToEstFromIndvFiles(std::string tumorFileList, std::string sconceEstimatesPath, int numExpectedLines, gsl_vector* meanVarianceCoefVec) const;
    gsl_vector* getIndOptimParamsToEstSummaryFromJointFile(std::string sconceEstimatesJointFile, int numExpectedLines) const;
    void getPairedOptimParamsToEstFromFiles(std::string pairedEstimatesPath, int numExpectedLinesPerFile);

};

#endif

