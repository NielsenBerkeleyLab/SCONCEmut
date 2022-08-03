#ifndef ALLPAIRS3TRPARAM2DEGPOLYHMM_HPP
#define ALLPAIRS3TRPARAM2DEGPOLYHMM_HPP

#include <boost/random.hpp>

#include <vector>
#include <unordered_map>
#include <gsl/gsl_sort.h>
#include <iomanip>

#include "Optimizable.hpp"
#include "TwoCell3TrParam2DegPolyHMM.hpp"
#include "AllCells3TrParam2DegPolyHMM.hpp"

/*
 * This class is for making a collection of HMMs that jointly estimates beta and lambda across a given number of pairs of
 * tumor cells, where each pair is independent of the other pairs (except for the shared parameters
 * and cell specific library size scaling factor (this is used every time a cell is used). Branch
 * lengths (t1/t2/t3) are estimated for each pair of cells
 *
 * Because there are so many pairs of cells, it makes sense to make the likelihood calc (at each iter
 * in BFGS) parallelizable. So, this class will hold a collection of TwoCell* HMMs where at each iteration
 * of BFGS, the appropriate parameters are set
 *
 * //this->paramsToEst = [lib0, lib1, ..., libN, beta, gamma, t1_cell0_1, t2_cell0_1, t3_cell0_1, t1_cell0_2, t2_cell0_2, t3_cell0_2, ..., t3_cell(N-1)_N]
 * this->paramsToEst = [lib0, lib1, ..., libN, beta, lambda, t1_cell0_1, t2_cell0_1, t3_cell0_1, t1_cell0_2, t2_cell0_2, t3_cell0_2, ..., t3_cell(N-1)_N]
 * this->fixedParams = [alpha]
 */
class AllPairs3TrParam2DegPolyHMM : public AllCells3TrParam2DegPolyHMM {
  private:
    int summaryMethodUsed;

  protected:
    // protected ctors
    AllPairs3TrParam2DegPolyHMM(std::vector<DepthPair*>* depths, std::vector<std::string>* sampleList, int maxPloidy);

    // member variables
    std::vector<std::unordered_map<std::string, gsl_matrix*>*>* margLikelihoodMatAcrossPairsVec; // cell num:[chr:summed forBackMargMat]
    std::vector<std::unordered_map<std::string, std::vector<double>*>*>* summaryPathAcrossPairsVec; // cell num:[chr:vector of decoded paths]
    std::vector<std::set<int>*>* hmmsCellAppearsInVec; // cellIdx:[vector of hmmIdx's the cell appears in]

    AllPairs3TrParam2DegPolyHMM(std::vector<DepthPair*>* depths, std::vector<std::string>* sampleList, gsl_vector* fixedParams, int maxPloidy, int numSharedTrParamsToEst, int numFixedTrParams, int numFixedLibs, int numPairs, int numBranchesToEst);
    virtual void makeHMMPairs(gsl_vector* meanVarianceCoefVec, bool preallocIntermediates = true);
    virtual void makeHMMPairs(gsl_vector* meanVarianceCoefVec, gsl_vector* transitionParams, bool preallocIntermediates = true);
    virtual void setHMMNames();

    // protected accessors and mutators
    virtual int getCell0IdxFromHMMIdx(int hmmIdx);
    virtual int getCell1IdxFromHMMIdx(int hmmIdx);
    virtual int getHMMIdxFromCellPair(int cell0Idx, int cell1Idx);
    virtual std::vector<HMM*>* getHMMsWithCellAsCell0(int cellIdx);
    virtual std::vector<HMM*>* getHMMsWithCellAsCell1(int cellIdx);
    virtual void storeHMMIdxForCells(int cell0Idx, int cell1Idx, int hmmIdx);
    virtual void storeHMMIdxForCell(int cellIdx, int hmmIdx);

    void decodeMargAllCellsAcrossPairs();
    void decodeMargOneCellAcrossPairs(int cellNum);

    void decodeStatSummaryAllCellsAcrossPairs(int summaryMethod);
    virtual void decodeStatSummaryOneCellAcrossPairs(int summaryMethod, int cellNum);

    std::string getSummaryMethodName(int summaryMethod);

    // methods for least squares for solving for Baum Welch estimates
    virtual void baumWelchLeastSquares_convertProbToParam(gsl_vector* dest, const gsl_vector* src) const override;
    virtual void baumWelchLeastSquares_convertParamToProb(gsl_vector* dest, const gsl_vector* src) const override;
    virtual double baumWelchLeastSquares_calcSumSqResid(const gsl_vector* v) override;
    virtual void saveBaumWelchEstsIntoParamsToEst(gsl_vector* varsToEst_probSpace, gsl_vector* initGuess) override;

  public:
    std::string pairedEstimatesPath; // path to read in .pairedParams files. needed for creating hmm's on the fly and reading in their params

    // constants
    const int NUM_PAIRS;

    static const int SUMMARY_MEAN     = 0;
    static const int SUMMARY_MEDIAN   = 1;
    static const int SUMMARY_MODE     = 2;
    static const int SUMMARY_MARGINAL = 3;

    // constructors and destructor
    virtual ~AllPairs3TrParam2DegPolyHMM();
    static AllPairs3TrParam2DegPolyHMM* create(std::vector<DepthPair*>* depths, std::vector<std::string>* sampleList, int maxPloidy, gsl_vector* meanVarianceCoefVec, bool preallocIntermediates = true);
    static AllPairs3TrParam2DegPolyHMM* create(const AllPairs3TrParam2DegPolyHMM& otherAllPairs3TrParam2DegPolyHMM);

    // accessors and mutators
    virtual void setLibScalingFactors(int cell0Idx, int cell1Idx, double lib0, double lib1);
    virtual void setAllLibScalingFactors(int cellNumInPair, double libScalingFactor);
    virtual double getLibScalingFactor(int cellNum) const;
    virtual double setAllTransition(gsl_vector* transitionParams) override;
    virtual void setParamsToEstFromIthHMM(gsl_vector* params, int hmmIdx) override;
    virtual void setFixedParamsFromIthHMM(gsl_vector* params, int hmmIdx) override;

    // override Optimizable methods
    virtual double setParamsToEst(gsl_vector* params) override;
    virtual void convertProbToParam(gsl_vector* dest, const gsl_vector* src) const override;
    virtual void convertParamToProb(gsl_vector* dest, const gsl_vector* src) const override;
    virtual void setSimParamsToEst(gsl_vector* params) override;
    virtual void setSimFixedParams(gsl_vector* params) override;
    virtual void miscFunctions() override;

    // BFGS
    virtual void setInitGuessNthTime(gsl_vector* initGuess, int iter, int numTotalRuns) const override;
    AllPairs3TrParam2DegPolyHMM* bfgs(gsl_vector* initGuess, int maxIters, bool verbose = true, bool debug = false) override;

    // methods to save results
    void saveAllCNAToBed(std::string filename) override;
    void summaryDecodeAllCellsAcrossPairs(int summaryMethod, bool runDecoding);
    void saveAllSummaryDecodedCNA(std::string filename);
    void saveAllSummaryDecodedCNAToBed(std::string filename);

    // methods to read intermediate results
    void getPairedOptimParamsToEstFromFiles(std::string pairedEstimatesPath, int numExpectedLinesPerFile);
    void getOnePairedOptimParamsToEstFromFile(int hmmIdx, std::string pairedEstimatesPath, int numExpectedLinesPerFile);
};

#endif

