#ifndef TWOCELL0TRPARAM2DEGPOLYHMMWITHMUTS_HPP
#define TWOCELL0TRPARAM2DEGPOLYHMMWITHMUTS_HPP

#include "TwoCell0TrParam2DegPolyHMM.hpp"
#include "MutationPair.hpp"
#include "MutationList.hpp"

/*
 * This class is uesd in the second stage of AllPairs2Stages3TrParam2DegPolyHMM
 * by the AllPairs0TrParam2DegPolyHMM class
 * This class also has a MutaionPair
 *
 * this->paramsToEst = [lib0, lib1, t1, t2, t3]
 * this->fixedParams = [alpha, beta, lambda, cnaToMutRateMu]
 */
class TwoCell0TrParam2DegPolyHMMWithMuts : public TwoCell0TrParam2DegPolyHMM {
  protected:
    TwoCell0TrParam2DegPolyHMMWithMuts(std::vector<DepthPair*>* depths, MutationPair* mutPair, gsl_vector* fixedParams, int maxPloidy, int numTrParamsToEst, int numFixedTrParams, int numFixedLibs, int numBranches, bool preallocIntermediates = true);
    //double cnaToMutRateMu; // consider moving to fixedParams (would need to add idx to HMM)
    MutationPair* mutPair;

  public:
    // constructors and destructor
    TwoCell0TrParam2DegPolyHMMWithMuts(std::vector<DepthPair*>* depths, MutationPair* mutPair, gsl_vector* fixedParams, int maxPloidy, bool preallocIntermediates = true);
    virtual ~TwoCell0TrParam2DegPolyHMMWithMuts();

    // accessors and mutators
    double getCnaToMutRateMu();
    void setCnaToMutRateMu(double rate);

    // functions that depend on model
    virtual TwoCell0TrParam2DegPolyHMMWithMuts* bfgs(gsl_vector* initGuess, int maxIters, bool verbose = true, bool debug = false) override;
    virtual void print(FILE* stream) override;
    virtual double getLogLikelihood() override;
    virtual double getMutLogLikelihood();
    virtual void estMutCountsPerBranch();
};

#endif

