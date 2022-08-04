#ifndef MUTATIONLIST_HPP
#define MUTATIONLIST_HPP

#include <unordered_map>
#include <vector>
#include <sstream>
#include <iostream>
#include <boost/filesystem.hpp> // for checking if files exist
#include <boost/algorithm/string.hpp> // for string splitting
#include <boost/math/special_functions/binomial.hpp> // for binomial coefficient, n choose k
#include <boost/math/special_functions/beta.hpp> // for beta function
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_nan.h>


/*
 * This class stores a list of mutations for one cell. Mutations are stored in a three maps (chr+coord: #ref reads; chr+coord: #alt reads; chr+coord: windowIdx)
 *
 * getLikelihood() returns the probability of observing read data at site i, given the allelic state, sconce CN, omega, and epsilon
 */
class MutationList {
  public:
    // constants
    static constexpr double SEQUENCING_ERROR_RATE = 0.005;

    // member variables
    double mutOverdispOmega; // in the beta binomial distribution, this is the overdispersion parameter, w. consider making a vector to store diff values of omega

    std::vector<std::string>* coordVec; // "chr:coord" list of all mutations
    std::unordered_map<std::string, int>* coordNumRefReadsMap; // "chr:coord": num ref reads
    std::unordered_map<std::string, int>* coordNumAltReadsMap; // "chr:coord": num alt reads
    std::unordered_map<std::string, int>* coordWindowIdxMap; // "chr:coord": window index, relative to chromosome, 0 based (ie is offset from start of chr)
    std::unordered_map<std::string, int>* coordWindowLineNumMap; // "chr:coord": line number relative to genome/entire depth file, 1 based (ie matches bed file exactly); for debugging
    std::unordered_map<std::string, int>* coordSconceCNMap; // "chr:coord": copy number estimate from sconce
    std::unordered_map<std::string, double>* coordBinomCoefMap; // "chr:coord": n choose numAlt binomial coefficient
    std::vector<double>* logCNvec; // cn: log(1/CN)

    // constructors and destructor
    MutationList(std::string mutationFilename, std::pair<std::unordered_map<std::string, std::vector<long long int>*>*, std::unordered_map<std::string, std::vector<int>*>*>* chrWindowIdxLineNumMaps); // chr:windowStartVal
    ~MutationList();

    // functions
    void setMutOverdispOmega(double omega);
    double getMutOverdispOmega();
    int findWindowIdx(std::string chr, long long int coord, std::unordered_map<std::string, std::vector<long long int>*>* chrWindowStartMap);
    void print(FILE* stream);
    void setCoordSconceCNMap(std::unordered_map<std::string, std::vector<int>*>* chrToViterbiPathMap); // chrToViterbiPathMap is chr:vector of viterbi decoded paths
    double getLikelihood(std::string site, bool isAltAllele);
    double getBinomCoef(std::string site);
    void setupBinomCoefs();
    double getLogCN(int cn);
    void setupLogCN(int maxCN);
};

#endif

