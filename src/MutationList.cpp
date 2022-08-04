#include "MutationList.hpp"

/*
 * chrWindowStartMap should be chr:vec of windowStartVal, needed to find windows idxs for each mutation
 */
MutationList::MutationList(std::string mutationFilename, std::pair<std::unordered_map<std::string, std::vector<long long int>*>*, std::unordered_map<std::string, std::vector<int>*>*>* chrWindowIdxLineNumMaps) {
  // set up member vars
  this->coordVec = new std::vector<std::string>();
  this->coordNumRefReadsMap = new std::unordered_map<std::string, int>();
  this->coordNumAltReadsMap = new std::unordered_map<std::string, int>();
  this->coordWindowIdxMap = new std::unordered_map<std::string, int>();
  this->coordWindowLineNumMap = new std::unordered_map<std::string, int>();
  this->coordSconceCNMap = new std::unordered_map<std::string, int>();
  this->setMutOverdispOmega(10);

  if(!boost::filesystem::exists(mutationFilename)) {
    std::cerr << "Error: " << mutationFilename << " does not exist. Exiting." << std::endl;
    exit(EXIT_FAILURE);
  }

  std::ifstream mutation(mutationFilename);
  std::string mutationLine;
  std::string key;
  std::string chr;
  long long int coord = -1;
  long long int end = -1;
  int numRefReads = -1;
  int numAltReads = -1;
  int windowIdx = -1;
  int windowLineNum = -1;
  std::unordered_map<std::string, std::vector<long long int>*>* chrWindowStartMap = chrWindowIdxLineNumMaps->first;
  std::unordered_map<std::string, std::vector<int>*>* chrWindowLineNumMap = chrWindowIdxLineNumMaps->second;
  // read in mutation file
  while(getline(mutation, mutationLine)) {
    // parse mutation line
    std::istringstream iss(mutationLine);
    iss >> chr >> coord >> end >> numRefReads >> numAltReads;
    key = chr + ":" + std::to_string(coord);

    // find window this mutation belongs to
    windowIdx = this->findWindowIdx(chr, coord, chrWindowStartMap);
    windowLineNum = (*(*chrWindowLineNumMap)[chr])[windowIdx];

    // store into appropriate maps
    this->coordVec->push_back(key);
    (*this->coordNumRefReadsMap)[key] = numRefReads;
    (*this->coordNumAltReadsMap)[key] = numAltReads;
    (*this->coordWindowIdxMap)[key] = windowIdx;
    (*this->coordWindowLineNumMap)[key] = windowLineNum;
  }
  mutation.close();
  this->coordBinomCoefMap = new std::unordered_map<std::string, double>();
  this->setupBinomCoefs();
  this->logCNvec = nullptr;
}

MutationList::~MutationList() {
  delete this->coordVec;
  delete this->coordNumRefReadsMap;
  delete this->coordNumAltReadsMap;
  delete this->coordWindowIdxMap;
  delete this->coordWindowLineNumMap;
  delete this->coordSconceCNMap;
  delete this->coordBinomCoefMap;
}

void MutationList::setMutOverdispOmega(double omega) {
  this->mutOverdispOmega = omega;
}
double MutationList::getMutOverdispOmega() {
  return this->mutOverdispOmega;
}

int MutationList::findWindowIdx(std::string chr, long long int coord, std::unordered_map<std::string, std::vector<long long int>*>* chrWindowStartMap) {
  std::vector<long long int>* currWindowStartVec = (*chrWindowStartMap)[chr];
  std::vector<long long int>::iterator bound = std::upper_bound(currWindowStartVec->begin(), currWindowStartVec->end(), coord);

  int windowIdx = -1;
  if(bound == currWindowStartVec->end()) {
    windowIdx = currWindowStartVec->size() - 1;
  }
  else {
    windowIdx = std::distance(currWindowStartVec->begin(), bound) - 1;
  }
  return windowIdx; // subtract 1 since windows are [start, end) and upper_bound returns the first element that is greater than coord https://stackoverflow.com/a/61971601
}

/*
 * prints tab delimited:
 * chr:coord	numRef	numAlt	max(sconceCN,-1)	windowIdx(relative to chr)	windowLineNum(relative to genome)
 */
void MutationList::print(FILE* stream) {
  /*std::string currCoord;
  int sconceCN = -1;
  fprintf(stream, "MutationList\n");
  for(std::vector<std::string>::iterator itr = this->coordVec->begin(); itr != this->coordVec->end(); ++itr) {
    currCoord = *itr;
    if(this->coordSconceCNMap->size() > 0) {
      sconceCN = (*this->coordSconceCNMap)[currCoord];
    }
    //fprintf(stream, "%s\t%i\t%i\t%i\t%i\t%i\n", currCoord.c_str(), (*this->coordNumRefReadsMap)[currCoord], (*this->coordNumAltReadsMap)[currCoord], sconceCN, (*this->coordWindowIdxMap)[currCoord], (*this->coordWindowLineNumMap)[currCoord]);
    double altLl = this->getLikelihood(currCoord, true);
    double refLl = this->getLikelihood(currCoord, false);
    fprintf(stream, "%s\t%i\t%i\t%i\t%i\t%i\t%0.5f\t%0.5f\t%0.5f\t%0.5f\n", currCoord.c_str(), (*this->coordNumRefReadsMap)[currCoord], (*this->coordNumAltReadsMap)[currCoord], sconceCN, (*this->coordWindowIdxMap)[currCoord], (*this->coordWindowLineNumMap)[currCoord], altLl, refLl, altLl / (altLl + refLl), refLl / (altLl + refLl));
  }*/
  fprintf(stream, "MutationList mutOverdispOmega:\n%0.40f\n\n", this->mutOverdispOmega);
}

/*
 * function to store per window sconce copy number estimates for each stored mutation.
 * chrToViterbiPathMap is chr:vector of viterbi decoded paths
 */
void MutationList::setCoordSconceCNMap(std::unordered_map<std::string, std::vector<int>*>* chrToViterbiPathMap) {
  int currWindowIdx = -1;
  int currCN = -1;
  int maxCN = 0;
  std::string currChr;
  std::string currCoord;
  std::vector<std::string> coordParts;
  for(std::vector<std::string>::iterator itr = coordVec->begin(); itr != coordVec->end(); ++itr) {
    currCoord = *itr;
    boost::split(coordParts, currCoord, boost::is_any_of(":"));
    currChr = coordParts[0];
    currWindowIdx = (*this->coordWindowIdxMap)[currCoord];
    currCN = (*(*chrToViterbiPathMap)[currChr])[currWindowIdx];
    (*this->coordSconceCNMap)[currCoord] = currCN;
    if(currCN > maxCN) {
      maxCN = currCN;
    }
  }
  this->setupLogCN(maxCN);
}

/*
 * function to calculate the likelihood of site, given if we should see the alt allele or not
 * uses the beta binomial distribution
 *
 * Assumes this->coordSconceCNMap has been set; returns nan if not (because default CN values are -1)
 *
 * manually calculates the probability, using boost's binomial coefficient and beta function (boost does not contain a beta binomial dist)
 * https://www.boost.org/doc/libs/1_78_0/libs/math/doc/html/math_toolkit/factorials/sf_binomial.html
 * https://www.boost.org/doc/libs/1_78_0/libs/math/doc/html/math_toolkit/sf_beta/beta_function.html
 *
 * p(D_ij | S_ij = 1) (S_ij == 1 if isAltAllele)
 * where D_ij = (a_ij = # alt reads, n_ij = # total reads) is the collection of reads in the i'th site, j'th cell
 * and S_ij is the true underlying allelic state
 *
 * A_ij is the number of alt reads in the i'th site, j'th cell
 * A_ij ~ BetaBinom(n=n_ij, alpha=wf, beta=w(1-f)), following notation from https://en.wikipedia.org/wiki/Beta-binomial_distribution
 * w=mutOverdispOmega, f=alleleFreq
 *
 * P(a_ij | n_ij, f, w) = (n_ij choose a_ij) * Beta(a_ij + wf, n_ij - a_ij + w(1-f)) / Beta(wf, w(1-f))
 * ie, the BetaBinomial distribution, where alpha=wf, beta=w(1-f)
 * alternatively, A_ij ~ Binom(n_ij, p); p ~ Beta(wf, w(1-f))
 * this gives E(A_ij) = n_ij * f
 *
 * let k=refCN, l=altCN; k+l = CN_ij
 *   f = l/(k+l)*(1-eps) + k/(k+l)*eps
 *
 * returns
 * if S_ij == 1 (isAltAlelle) {
 *   sum_(0<=k<CN) P(D_ij | f, w, k, CN_ij) * P(k | CN_ij)
 *   P(k | CN_ij) = 1/CN_ij
 * else
 *   P(D_ij | f, w, k=CN_ij, CN_ij)
 */
double MutationList::getLikelihood(std::string site, bool isAltAllele) {
  // if missing data (ie this cell does not have data for this site), return 0 (ie no contribution to the likelihood function). if n_ij = A_ij = 0; P(D_ij | f, w, k, CN_ij) = 0
  if(this->coordNumRefReadsMap->count(site) == 0) {
    return 1;
    return 0;
  }
  double likelihood = 0;
  int altCN = 0; // l
  int refCN = 0; // k
  double altAlleleFreq = 0; // f
  double alpha = 0; // for beta function
  double beta = 0; // for beta function

  double numRefReads = (*this->coordNumRefReadsMap)[site]; // R_ij = n_ij - A_ij, where n_ij is the total read depth at this site and cell
  double numAltReads = (*this->coordNumAltReadsMap)[site]; // A_ij
  double numTotalReads = numRefReads + numAltReads;
  double nChooseA = this->getBinomCoef(site); // (n_ij choose A_ij)

  // if infer copy number 0, return 0 (no contribution to the likelihood function)
  int sconceCN = (*this->coordSconceCNMap)[site]; // CN_ij
  if(sconceCN == 0) {
    return 0;
  }
  double logSconceCN = this->getLogCN(sconceCN);
  double currLl = 0;

  // if true underlying state is the alternate allele, then we sum P(A_ij | D_ij, f, mutOverdispOmega) over 0<=k<CN,0<l<=CN, k+l=CN, with f = l/(k+l)*(1-eps) + k/(k+l)*eps
  if(isAltAllele) {
    for(refCN = 0; refCN < sconceCN; refCN++) {
      altCN = sconceCN - refCN;
      altAlleleFreq = (double) altCN / (double) (sconceCN) * (1.0 - this->SEQUENCING_ERROR_RATE) + (double) refCN / (double) (sconceCN) * this->SEQUENCING_ERROR_RATE;
      alpha = this->mutOverdispOmega * altAlleleFreq; // alpha = wf
      beta = this->mutOverdispOmega * (1.0 - altAlleleFreq); // beta = w(1-f)
      if(alpha < 0 || beta < 0 || (numAltReads + alpha) < 0 || (numTotalReads - numAltReads + beta) < 0) {
        return GSL_NAN;
      }
      currLl = exp(nChooseA + gsl_sf_lnbeta(numAltReads + alpha, numTotalReads - numAltReads + beta) - gsl_sf_lnbeta(alpha, beta) + logSconceCN); // all allele freqs equally likely
      likelihood += currLl;
    }
  }
  // else is ref allele, we set k=CN,l=0
  else {
    altCN = 0;
    refCN = sconceCN;
    altAlleleFreq = (double) altCN / (double) (sconceCN) * (1.0 - this->SEQUENCING_ERROR_RATE) + (double) refCN / (double) (sconceCN) * this->SEQUENCING_ERROR_RATE;
    alpha = this->mutOverdispOmega * altAlleleFreq; // alpha = wf
    beta = this->mutOverdispOmega * (1.0 - altAlleleFreq); // beta = w(1-f)
    if(alpha < 0 || beta < 0 || (numAltReads + alpha) < 0 || (numTotalReads - numAltReads + beta) < 0) {
      return GSL_NAN;
    }
    likelihood = exp(nChooseA + gsl_sf_lnbeta(numAltReads + alpha, numTotalReads - numAltReads + beta) - gsl_sf_lnbeta(alpha, beta));
  }
  return likelihood;
}

/*
 * function to get precalculated binomial coefficient for a specified site
 */
double MutationList::getBinomCoef(std::string site) {
  return (*this->coordBinomCoefMap)[site];
}

/*
 * function to precalculate binomial coefficients for all sites. Stores in log space
 */
void MutationList::setupBinomCoefs() {
  std::string site;
  double numRefReads = 0;
  double numAltReads = 0;
  double numTotalReads = 0;
  for(unsigned int i = 0; i < this->coordVec->size(); i++) {
    site = (*this->coordVec)[i];
    numRefReads = (*this->coordNumRefReadsMap)[site]; // R_ij = n_ij - A_ij, where n_ij is the total read depth at this site and cell
    numAltReads = (*this->coordNumAltReadsMap)[site]; // A_ij
    numTotalReads = numRefReads + numAltReads; // n_ij
    (*this->coordBinomCoefMap)[site] = gsl_sf_lnchoose(numTotalReads, numAltReads);
  }
}

/*
 * function to get log(1/cn) from lookup
 */
double MutationList::getLogCN(int cn) {
  return (*this->logCNvec)[cn];
}
/*
 * function to setup logCNvec, which stores log(1/cn) for up to observed maxCN
 */
void MutationList::setupLogCN(int maxCN) {
  if(this->logCNvec == nullptr) {
    this->logCNvec = new std::vector<double>(maxCN + 1);
  }
  else if((int) this->logCNvec->size() != maxCN + 1) {
    this->logCNvec->resize(maxCN);
  }
  for(int i = 0; i <= maxCN; i++) {
    (*this->logCNvec)[i] = log(1.0 / i);
  }
}

