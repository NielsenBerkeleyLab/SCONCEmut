#include <iostream>
#include <cmath>
#include <string>
#include <sstream>
#include <vector>
#include <algorithm>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_errno.h>


/*
 * Sun 17 Jul 2022 12:36:29 PM PDT
 * program to calculate LLR for tumor sites, given omega (overdispersion) value
 * H0: site is not variable, so allele freq of minor/derived allele  = 0
 * H1: site is variable, so allele freq of minor/derived allele != 0
 *
 * g++ fitBetaBinomTumor.cpp `gsl-config --cflags --libs` -o fitBetaBinomTumor
 *
 * input: tab separated
 * <chr> <start> <end> <dbsnp> <cell0_totalNumReads> <cell0_numMajorAlleleReads> <cell0_llr> <...>
 *
 * prints to stdout
 * <chr> <start> <end> <dbsnp> <llr>
 *
 * sample usage:
 *   ./fitBetaBinomTumor dbsnpIntersection_tumor/pooledDiploidTumorMajorAllele.llr.dbsnp.countSummary.minCells5_minAvgReads3.5.readsMat.bed dbsnpIntersection_tumor/pooledDiploidTumorMajorAllele.llr.dbsnp.header > dbsnpIntersection_tumor/pooledDiploidTumorMajorAllele.llr.dbsnp.countSummary.minCells5_minAvgReads3.5.readsMat.pooledLlr.bed
 *   cat dbsnpIntersection_tumor/pooledDiploidTumorMajorAllele.llr.dbsnp.countSummary.minCells5_minAvgReads3.5.readsMat.bed | ./fitBetaBinomTumor > dbsnpIntersection_tumor/pooledDiploidTumorMajorAllele.llr.dbsnp.countSummary.minCells5_minAvgReads3.5.readsMat.pooledLlr.bed
 *
 */

//double calcBetaBinomLl(double f, void* params) {
double calcBetaBinomLl(double f, std::pair<int, int>* pair) {
  //std::pair<int, int>* pair = (std::pair<int, int>*) params;
  int numTotalReads = pair->first;
  int numPooledMajorReads = pair->second; // number of reads mapping to the major allele (relative to the pooled diploid samples)
  int numDerivedReads = numTotalReads - numPooledMajorReads;
  double omega = 8.68008260176282497866; // omega, from fitBetaBinomNull.R, for dbsnpIntersection/pooledDiploidMajorAllele.llr.dbsnp.llr32.notInDbsnp.minCells5_minAvgReads3.5.bed

  double alpha = f * omega;
  double beta = (1-f) * omega;
  //std::cout << "alpha: " << alpha << ", beta: " << beta << std::endl;
  if(alpha < 0 || beta < 0 || numDerivedReads + alpha < 0 || numTotalReads - numDerivedReads + beta < 0) {
    return GSL_NAN;
  }

  double ll = gsl_sf_lnchoose(numTotalReads, numDerivedReads) + gsl_sf_lnbeta(numDerivedReads + alpha, numTotalReads - numDerivedReads + beta) - gsl_sf_lnbeta(alpha, beta);
  //std::cout << numDerivedReads << ", " << numTotalReads << ", " <<ll << std::endl;
  return ll;
}
//double calcBetaBinomLlOnVec(double f, std::vector<std::pair<int, int>>* readCountsVec) {
double calcBetaBinomLlOnVec(double f, void* params) {
  std::vector<std::pair<int, int>>* readCountsVec = (std::vector<std::pair<int, int>>*) params;
  double totalLl = 0;
  for(std::vector<std::pair<int, int>>::iterator itr = readCountsVec->begin(); itr != readCountsVec->end(); ++itr) {
    totalLl += calcBetaBinomLl(f, &(*itr));
  }
  return -totalLl;
}

double optimAlleleFreq(std::vector<std::pair<int, int>>* readCountsVec) {
  int status;
  int iter = 0;
  int max_iter = 100;
  double lowerBound = 1e-8; // allele frequencies can range form 0 to 1
  double upperBound = 1-1e-8;
  double f = 0.5;

  gsl_function my_func;
  my_func.function = &calcBetaBinomLlOnVec;
  my_func.params = readCountsVec;

  const gsl_min_fminimizer_type* T = gsl_min_fminimizer_brent;
  gsl_min_fminimizer* s = gsl_min_fminimizer_alloc(T);
  gsl_min_fminimizer_set(s, &my_func, f, lowerBound, upperBound);
  //printf ("%5d [%.7f, %.7f] "
  //    "%.7f %.7f\n",
  //    iter, lowerBound, upperBound,
  //    f, upperBound - lowerBound);

  gsl_set_error_handler_off();

  do {
    iter++;
    status = gsl_min_fminimizer_iterate(s);

    f = gsl_min_fminimizer_x_minimum(s);
    lowerBound = gsl_min_fminimizer_x_lower(s);
    upperBound = gsl_min_fminimizer_x_upper(s);

    status = gsl_min_test_interval(lowerBound, upperBound, 0.001, 0.0);

    //if(status == GSL_SUCCESS) {
    //  printf("Converged:\n");
    //}

    //printf ("%5d [%.7f, %.7f] "
    //    "%.7f %.7f\n",
    //    iter, lowerBound, upperBound,
    //    f, upperBound - lowerBound);
  }
  while(status == GSL_CONTINUE && iter < max_iter);
  if(status == GSL_EINVAL) {
    return GSL_NAN;
  }
  f = gsl_min_fminimizer_x_minimum(s);

  gsl_min_fminimizer_free (s);
  return f;
}


/*
 * check if this site is variable. If all pairs have the same num total and pooled major reads, then not variable
 */
bool isVariable(std::vector<std::pair<int, int>>* readCountsVec) {
  return !std::all_of(readCountsVec->begin(), readCountsVec->end(),
      [&](std::pair<int, int> pair) -> bool {
      return (pair.first == pair.second);
      });
}

//int main(char* argv, int argc) {
int main() {
  // constants
  double errorRate = 0.005;
  //double logErr = log(errorRate);
  //double log1Err = log(1.0 - errorRate);
  //double denomErr = 1.0 - 2.0 * errorRate;

  //// read in header to get number of fields
  //for(std::string line; std::getline(std::cin, line);) {
  //  std::stringstream ss(line);
  //  while(ss >> buf) {
  //    tokens.push_back(buf)
  //  }

  std::string chr;
  std::string start;
  std::string end;
  std::string dbsnp;
  std::string numTotalReads;
  std::string numPooledMajorReads; // number of reads mapping to the major allele (relative to the pooled samples)
  std::string indvLlr; // ignored
  std::vector<std::pair<int, int>>* readCountsVec = new std::vector<std::pair<int, int>>();

  double f = 0;
  double llr = 0;
  double altHypLogLik = 0;
  double nullHypLogLik = 0;

  // loop reading from https://stackoverflow.com/a/10464355
  for(std::string line; std::getline(std::cin, line);) {
    // parse input
    std::istringstream iss(line);
    iss >> chr >> start >> end >> dbsnp;

    //std::cout << line << std::endl;
    readCountsVec->clear();
    do {
      iss >> numTotalReads >> numPooledMajorReads >> indvLlr;
      //std::cout << numPooledMajorReads << ", " << numPooledMajorReads << ", " << indvLlr << std::endl;
      if(numTotalReads != "NA") {
        //std::cout << "inserting (" << numTotalReads << ", " << numPooledMajorReads << ")" << std::endl;
        readCountsVec->push_back(std::make_pair(stoi(numTotalReads), stoi(numPooledMajorReads)));
      }
    } while(iss);

    // skip over non variable sites
    if(!isVariable(readCountsVec)) {
      continue;
    }
    //for(std::vector<std::pair<int, int>>::iterator itr = readCountsVec->begin(); itr != readCountsVec->end(); ++itr) {
    //  //for(int i = 0; i < readCountsVec->size(); i++) {
    //  //std::cout << (*itr).first << ", " << (*itr).second << std::endl;
    //  std::cout << itr->first << ", " << itr->second << std::endl;
    //  //std::pair<int, int> pair = (*readCountsVec)[i];
    //  //std::cout << pair.first << ", " << pair.second << std::endl;
    //  //std::cout << (*itr) << std::endl;
    //}


    // estimate allele frequency by maximizing beta binom loglikelihood
    f = optimAlleleFreq(readCountsVec);
    // if optim failed, the site was probably not variable (probably all pairs equal except for one)
    if(gsl_isnan(f)) {
      continue;
    }
    altHypLogLik = -calcBetaBinomLlOnVec(f, readCountsVec); // need to (un)negate; sign flipped for minimization

    // calculate loglikelihood of allele frequency under the null model: site is not variable, so minor allele freq is 0
    nullHypLogLik = -calcBetaBinomLlOnVec(errorRate, readCountsVec);

    //std::cout << altHypLogLik << ", " << nullHypLogLik << std::endl;
    llr = altHypLogLik - nullHypLogLik;
    // print
    // <chr> <start> <end> <dbsnp> <llr>
    printf("%s\t%s\t%s\t%s\t%0.10f\t%0.10f\n", chr.c_str(), start.c_str(), end.c_str(), dbsnp.c_str(), f, llr);

    //// if pooled and cell specific major alleles match, then #maj = #total - #min
    //if(pooledMajorAllele == cellMajorAllele) {
    //  numPooledMajorReads = numTotalReads - numCellMinorReads;
    //}
    //// else, don't match, so cell's minor is assumed to be pooled major allele, so read counts match
    //else {
    //  numPooledMajorReads = numCellMinorReads;
    //}

    //// f_a = min{ (numPooledMajorReads - errorRate * numTotalReads) / ((1-2*errorRate) * numTotalReads), 1}
    //// LLR = l(f_a) - l(f_a = 1)
    ////     = numPooledMajorReads * (log(errorRate + f_a - 2*errorRate*f_a) - log(1-errorRate)) - (numTotalReads - numPooledMajorReads) * (log(errorRate) - log(1-errorRate - f_a + 2*errorRate*f_a))

    //f_a = (numPooledMajorReads - errorRate * numTotalReads) / (denomErr * numTotalReads);
    //if(f_a > 1) {
    //  f_a = 1;
    //}
    //else if(f_a < 0) {
    //  f_a = 0;
    //}

    //llr = numPooledMajorReads * (log(errorRate + f_a - 2.0 * errorRate * f_a) - log1Err) - (numTotalReads - numPooledMajorReads) * (logErr - log(1.0 - errorRate - f_a + 2.0 * errorRate * f_a));


    //// print output
    //// <chr> <start> <end> <pooledMajorAllele> <cellMajorAllele> <numTotalReads> <numPooledMajorReads> <LLR>
    //printf("%s\t%s\t%s\t%s\t%s\t%i\t%i\t%0.5f\n", chr.c_str(), start.c_str(), end.c_str(), pooledMajorAllele.c_str(), cellMajorAllele.c_str(), numTotalReads, numPooledMajorReads, llr);
    }
    return 0;
  }


