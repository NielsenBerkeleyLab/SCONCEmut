#ifndef TWOCELLFIXLIB0TRPARAM2DEGPOLYHMMWITHMUTS_HPP
#define TWOCELLFIXLIB0TRPARAM2DEGPOLYHMMWITHMUTS_HPP

#include "TwoCellFixLib0TrParam2DegPolyHMM.hpp"
#include "MutationPair.hpp"

/*
 * This class is uesd in the second stage of sconceAllPairs
 * by the AllPairs0TrParam2DegPolyHMM class, when library sizes are fixed
 * alpha, beta, gamma are also all fixed
 *
 * In this class, we also include the loglik from somatic mutations.
 * That is, l(t1,t2,t3) = l(t1,t2,t3)_cna + l(t1,t2,t3)_mut
 * l(t1,t2,t3)_cna = forward loglikelihood
 * l(t1,t2,t3)_mut = log(P(cnaToMutRateMu*t1*X1) + P(cnaToMutRateMu*t2*X2) + P(cnaToMutRateMu*t3*X3))
 * where x1,x2,x3 are first estimated using bfgs on MutationPair (l(X1,X2,X3)), and then fixed
 *
 * this->paramsToEst = [t1, t2, t3]
 * this->fixedParams = [lib0, lib1, alpha, beta, lambda, cnaToMutRateMu]
 */
class TwoCellFixLib0TrParam2DegPolyHMMWithMuts : public TwoCellFixLib0TrParam2DegPolyHMM {
  protected:
    TwoCellFixLib0TrParam2DegPolyHMMWithMuts(std::vector<DepthPair*>* depths, MutationPair* mutPair, gsl_vector* fixedParams, int maxPloidy, int numTrParamsToEst, int numFixedTrParams, int numFixedLibs, int numBranches, bool preallocIntermediates = true);
    //double cnaToMutRateMu; // consider moving to fixedParams (would need to add idx to HMM)
    MutationPair* mutPair;

  public:
    // constructors and destructor
    TwoCellFixLib0TrParam2DegPolyHMMWithMuts(std::vector<DepthPair*>* depths, MutationPair* mutPair, gsl_vector* fixedParams, int maxPloidy, bool preallocIntermediates = true);
    virtual ~TwoCellFixLib0TrParam2DegPolyHMMWithMuts();

    // accessors and mutators
    double getCnaToMutRateMu();
    void setCnaToMutRateMu(double rate);

    // functions that depend on model
    virtual TwoCellFixLib0TrParam2DegPolyHMMWithMuts* bfgs(gsl_vector* initGuess, int maxIters, bool verbose = true, bool debug = false) override;
    virtual void print(FILE* stream) override;
    virtual double getLogLikelihood() override;

    virtual double getMutLogLikelihood();
    virtual void estMutCountsPerBranch();
    virtual void saveParamEstimates(std::string filename) const override;
};

#endif

