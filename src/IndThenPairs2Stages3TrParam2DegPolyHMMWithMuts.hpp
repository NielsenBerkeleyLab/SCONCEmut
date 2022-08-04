#ifndef INDTHENPAIRS2STAGES3TRPARAM2DEGPOLYHMMWITHMUTS_HPP
#define INDTHENPAIRS2STAGES3TRPARAM2DEGPOLYHMMWITHMUTS_HPP

#include <boost/random.hpp>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_sort_vector.h>

#include "Optimizable.hpp"
#include "DepthPair.hpp"
#include "MutationListContainer.hpp"
#include "MutationList.hpp"
#include "MutationInd.hpp"
#include "MutationPair.hpp"
#include "MutationJointOverdisp.hpp"
#include "MutationJointMutRateOverdisp.hpp"

#include "IndThenPairs2Stages3TrParam2DegPolyHMM.hpp"

#include "AllCells3TrParam2DegPolyHMM.hpp"
#include "AllInd2Stages3TrParam2DegPolyHMM.hpp"
#include "AllPairs0TrParam2DegPolyHMM.hpp"
#include "AllPairsFixLib0TrParam2DegPolyHMM.hpp"
#include "SelectPairs0TrParam2DegPolyHMM.hpp"
#include "SelectPairsFixLib0TrParam2DegPolyHMM.hpp"

#include "AllPairs0TrParam2DegPolyHMMWithMuts.hpp"
#include "AllPairsFixLib0TrParam2DegPolyHMMWithMuts.hpp"
#include "SelectPairs0TrParam2DegPolyHMMWithMuts.hpp"
#include "SelectPairsFixLib0TrParam2DegPolyHMMWithMuts.hpp"


/*
 * This class is for adding mutations to IndThenPairs2Stages3TrParam2DegPolyHMM.
 *
 * SCONCE as stage 0/1 to est [lib, beta, lambda, t] for each cell independently, including bw+ls and bfgs
 * Then mutations are added to indv cells, and [mu] is est for each cell. mu is the connecting rate between CNAs and mutations
 *
 * Because all parameters are relative, alpha is fixed in all cases.
 *
 * Then, stage 2 is all pairs of cells. libs = as is, beta = median(beta_A, beta_B,...), lambda = median(lambda_A, lambda_B,...), mu = median(mu_A, mu_B,...).
 * BFGS is used to estimate [t1, t2, t3]_AB, [t1, t2, t3]_AC,..., where l(t1,t2,t3) = l_CNA(t1,t2,t3) + l_SNP(t1,t2,t3)
 * Starting points for t's, for pair AB, is:
 *   t1 = min(t_A, t_B) / 2
 *   t2 = t_A - t1
 *   t3 = t_B - t1
 * Then should get good ests for these params, can output some sort of summary CNA calls
 *
 */
class IndThenPairs2Stages3TrParam2DegPolyHMMWithMuts : public IndThenPairs2Stages3TrParam2DegPolyHMM {
  protected:
    // member variables
    std::vector<MutationList*>* mutListVec; // just a list of coord:counts, etc
    std::vector<MutationInd*>* mutIndVec; // Optimizable objects, used to estimate cnaToMutRateMu
    MutationJointOverdisp* mutJointOverdisp; // Optimizable object, used to jointly estimate omega (overdispersion factor for beta binomial) across indv cells
    MutationJointMutRateOverdisp* mutJointMutRateOverdisp; // Optimizable object, used to jointly estimate global mu (mutation to cna rate) and omega (overdispersion factor for beta binomial) across indv cells
    gsl_vector* allIndEstMutCounts;
    gsl_vector* allPairedEstMutCounts;
    double cnaToMutRateMu; // this scales branch lengths to mutation count estimates. also in MutationInd (as member var) and All/SelectPair/TwoCell*WithMuts (in fixedParams)
    double mutOverdispOmega; // in the beta binomial distribution, this is the overdispersion parameter. also stored directly in MutationList (as member var). consider making a vector to store diff values of omega

  public:
    // constructors and destructor
    IndThenPairs2Stages3TrParam2DegPolyHMMWithMuts(std::vector<DepthPair*>* depths, std::vector<std::string>* sampleList, std::vector<MutationList*>* mutListVec, int maxPloidy, bool runIndBFGS, int numPairsStage2, bool verbose, bool gradientDebug, bool centralDiff);
    virtual ~IndThenPairs2Stages3TrParam2DegPolyHMMWithMuts();

    // accessors and mutators
    virtual void print(FILE* stream) override;
    double getCnaToMutRateMu();
    double getMutOverdispOmega();
    gsl_vector* getAllIndEstMutCounts();
    gsl_vector* getAllPairedEstMutCounts();

    // optimization methods
    virtual void estCnaToMutRateMuIndependentlyForAllIndCells();
    virtual void estMutOverdispOmegaJointlyAcrossAllIndCells();
    virtual void estCnaToMutRateMuAndMutOverdispOmegaJointlyAcrossAllIndCells(std::string sconceEstimatesPath);
    virtual void estimateMutParamsFromIndCells(std::string sconceEstimatesPath, std::string filename);
    virtual void reestIndBranchLengths();
    virtual void copyStage1EstsIntoStage2FixedParams(gsl_vector* stage2FixedParams, gsl_vector* stage2InitGuess, gsl_vector* indEstParams, int numStage2LibsToFix) override;
    virtual void setUpAllPairsBFGS(gsl_vector* fixedParams, gsl_vector* initGuess, gsl_vector* meanVarianceCoefVec, bool fixLib = true, bool preallocIntermediates = false) override;
    virtual void setUpSelectPairsBFGS(gsl_vector* fixedParams, gsl_vector* initGuess, int numPairsToSummarize, int ordering, int seed, gsl_vector* meanVarianceCoefVec, bool fixLib = true, bool preallocIntermediates = false) override;
    virtual AllPairs0TrParam2DegPolyHMM* bfgsStage2(gsl_vector* initGuess, int maxIters, std::string filename, bool verbose = true, bool debug = false) override;

};

#endif

