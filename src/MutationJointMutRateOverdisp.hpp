#ifndef MUTATIONJOINTMUTRATEOVERDISP_HPP
#define MUTATIONJOINTMUTRATEOVERDISP_HPP

#include <vector>

#include "MutationListContainer.hpp"
#include "MutationList.hpp"
#include "MutationInd.hpp"

/*
 * This class estimates both cnaToMutRateMu (MutationInd param) and mutOverdispOmega (MutationList param)  across all cells jointly, given branch length and copy number estimates from SCONCE
 * In each optimization step, omega is set for each MutationInd object, then the number of mutations (x) is estimated for that MutationInd obejct. Then, the probability of observing x mutations on branch length t, with cnaToMutRateMu, is calculated (Poisson like dist). This log probability is summed across all cells
 *
 * getLogLikelihood() returns the joint (across all Mutation objects) probability of mu and omega
 *
 * paramsToEst = [mu, omega]
 */
class MutationJointMutRateOverdisp : public MutationListContainer {
  public:
    // member variables
    std::vector<MutationInd*>* mutIndVec;
    gsl_vector* branchLengths;

    // constructors and destructor
    MutationJointMutRateOverdisp(std::vector<MutationList*>* mutListVec, std::vector<std::unordered_map<std::string, std::vector<int>*>*>* chrToViterbiPathMapVec, gsl_vector* branchLengths, bool verbose = true, bool gradientDebug = false);
    virtual ~MutationJointMutRateOverdisp();

    // functions
    virtual double getCnaToMutRateMu();
    virtual double getMutOverdispOmega();
    virtual gsl_vector* getAllIndEstMutCounts();

    // override Optimizable functions
    virtual double getLogLikelihood() override;
    virtual void print(FILE* stream) override;
    virtual Optimizable* bfgs(gsl_vector* initGuess, int maxIters, bool verbose = true, bool debug = false) override; // this one should construct a subclass Optimizable, convert initGuess into BFGS space, and call Optimizable::bfgs

    void callIndvMutCountEst(int mutIndIdx);
    virtual double getValidOptimParamMax() const override;

    double getMutParamsFromFile(std::vector<std::string>* sampleList, int maxPloidy, std::string sconceEstimatesPath);
};

#endif

