#ifndef MUTATIONLISTCONTAINER_HPP
#define MUTATIONLISTCONTAINER_HPP

#include <unordered_map>
#include <vector>
#include <set>

#include "Optimizable.hpp"
#include "MutationList.hpp"

/*
 * This class is the super class for storing MutationLists. It extends Optimizable, and the function it maximizes depends on the subclass's implementation of getLoglikelihood()
 * Assumes all params should be non negative
 *
 * Note MutationJointOverdisp stores MutationInd objects, not MutationList objects directly
 *
 * paramsToEst = [x] (MutationInd), [x1,x2,x3] (MutationPair), [w] (MutationJointOverdisp), [mu, omega] (MutationJointMutRateOverdisp)
 */
class MutationListContainer : public Optimizable {
  public:
    // constants
    const int NUM_MUTATION_COUNTS_TO_EST; // x, x1/x2/x3
    const int NUM_MUTATION_VARS_TO_EST; // mu, omega
    const int MUTATION_COUNT_START_IDX;
    const int MUTATION_VARS_START_IDX;

    // member variables
    std::vector<MutationList*>* mutListVec; // list of MutationLists, one per cell
    gsl_vector* internalInitGuess;

    // constructors and destructor
    MutationListContainer(std::vector<MutationList*>* mutListVec, int numMutCountsToEst, int numOverdispVarsToEst, bool verbose = true, bool gradientDebug = false);
    virtual ~MutationListContainer();

    // functions
    virtual MutationListContainer* estMutCountsPerBranch(std::vector<std::unordered_map<std::string, std::vector<int>*>*>* chrToViterbiPathMapVec, gsl_vector* initGuess, int maxIters, bool verbose);
    virtual void setAllCoordSconceCNMap(std::vector<std::unordered_map<std::string, std::vector<int>*>*>* chrToViterbiPathMapVec);
    virtual double getMutCountEst(int branchIdx) const;
    virtual double getMutationVarEst(int varIdx);
    virtual void setAllMutOverdispOmega(double omega);
    virtual void saveMutParamsToFile(int mutListIdx, std::string filename);
    virtual void saveMutParamsToFile(int mutListIdx, std::string filename, const char* mode);
    virtual double getMutParamsFromFile(std::string filename, int numExpectedLinesPerFile);

    // override Optimizable functions
    virtual double setParamsToEst(gsl_vector* params) override;
    virtual void convertProbToParam(gsl_vector* dest, const gsl_vector* src) const override;
    virtual void convertParamToProb(gsl_vector* dest, const gsl_vector* src) const override;
    virtual void setInitGuessNthTime(gsl_vector* initGuess, int iter, int numTotalRuns) const override;
    virtual Optimizable* bfgs(gsl_vector* initGuess, int maxIters, bool verbose = true, bool debug = false) override; // this one should construct a subclass Optimizable, convert initGuess into BFGS space, and call Optimizable::bfgs
    virtual void print(FILE* stream) override;
    virtual double checkOptimProbValidity(gsl_vector* probs) const override;
    virtual double getValidOptimParamMax() const = 0; // ex estimated x1+x2+x3 shouldn't be greater than N
    virtual double checkStateValidity(double epsilon = 1e-8) const override; // should check if transition matrices are ok
    virtual void setSimParamsToEst(gsl_vector* params) override;
    virtual void setSimFixedParams(gsl_vector* params) override;
    virtual void miscFunctions() override;

};

#endif

