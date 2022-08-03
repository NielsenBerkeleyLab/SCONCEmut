#ifndef MUTATIONIND_HPP
#define MUTATIONIND_HPP

#include <vector>

#include "MutationListContainer.hpp"
#include "MutationList.hpp"

/*
 * This class stores one individual MutationList, where the mutation count is estimated (independently) for one cell. This is useful for estimating cnaToMutRateMu
 *
 * getLikelihood() returns the probability of mutation count x, given read data (MutationLists with fixed omega and epsilon))
 *
 * paramsToEst = [x]
 */
class MutationInd : public MutationListContainer {
  public:
    // member variables
    double cnaToMutRateMu;
    bool hasMutEsts;

    // constructors and destructor
    MutationInd(std::vector<MutationList*>* mutListVec, bool verbose = true, bool gradientDebug = false);
    virtual ~MutationInd();

    // functions
    virtual MutationListContainer* estMutCountsPerBranch(std::vector<std::unordered_map<std::string, std::vector<int>*>*>* chrToViterbiPathMapVec, gsl_vector* initGuess, int maxIters, bool verbose) override;
    using MutationListContainer::estMutCountsPerBranch;
    double getMutCountEst();
    using MutationListContainer::getMutCountEst;
    void setCnaToMutRateMu(double mu);
    double getCnaToMutRateMu();
    double getMutOverdispOmega(int mutListIdx = 0);
    virtual void saveMutParamsToFile(int mutListIdx, std::string filename) override;
    virtual double getMutParamsFromFile(std::string filename, int numExpectedLinesPerFile) override;

    // override Optimizable functions
    virtual double getLogLikelihood() override;

    virtual double getValidOptimParamMax() const override; // ex estimated x shouldn't be greater than N
};

#endif

