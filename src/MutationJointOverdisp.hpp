#ifndef MUTATIONJOINTOVERDISP_HPP
#define MUTATIONJOINTOVERDISP_HPP

#include <vector>

#include "MutationListContainer.hpp"
#include "MutationList.hpp"
#include "MutationInd.hpp"

/*
 * This class stores a vector of MutationInd objects (ie *not* MutationList objects), where the overdispersion parameter omega is estimated jointly across all cell
 *
 * getLikelihood() returns the joint (across all MutationInd objects) probability of omega, given MutationInd objects (with fixed mu estimates)
 *
 * paramsToEst = [omega]
 */
class MutationJointOverdisp : public MutationListContainer {
  public:
    // member variables
    std::vector<MutationInd*>* mutIndVec;

    // constructors and destructor
    MutationJointOverdisp(std::vector<MutationInd*>* mutIndVec, bool verbose = true, bool gradientDebug = false);
    virtual ~MutationJointOverdisp();

    // override Optimizable functions
    virtual double getLogLikelihood() override;

    virtual double getValidOptimParamMax() const override;
};

#endif

