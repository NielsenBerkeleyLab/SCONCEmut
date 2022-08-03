#ifndef MUTATIONPAIR_HPP
#define MUTATIONPAIR_HPP

#include <unordered_map>
#include <vector>
#include <set>

#include "MutationListContainer.hpp"
#include "MutationList.hpp"

/*
 * This class stores pairs of individual MutationLists (one for each cell). It keeps track of which mutations are shared between these cells
 *
 * getLikelihood() returns the probability of mutation counts x1,x2,x3, given MutationList objects (fixed mu and epsilon values)
 * calculates l(X1,X2,X3) = sum_sites [log(p(data_i | S_i=(1,1)) * X1/N + p(data_i | S_i=(1,0)) * X2/N + p(data_i | S_i=(0,1))*X3/N + p(data_i | S_i=(0,0))*(N-X1-X2-X3)/N)]
 *
 * paramsToEst = [x1,x2,x3]
 */
class MutationPair : public MutationListContainer {
  public:
    // member variables
    std::vector<std::string>* sharedMutationVec; // chr+coord of shared mutations Wed 06 Apr 2022 04:28:03 PM PDT all mutations in either cell, not necessarily both

    // constructors and destructor
    MutationPair(std::vector<MutationList*>* mutListVec, gsl_vector* initGuess, bool verbose = true, bool gradientDebug = false);
    virtual ~MutationPair();

    // override Optimizable functions
    virtual double getLogLikelihood() override;
    virtual void print(FILE* stream) override;
    virtual double checkOptimProbValidity(gsl_vector* probs) const override;

    virtual double getValidOptimParamMax() const override; // ex estimated x1+x2+x3 shouldn't be greater than N
};

#endif

