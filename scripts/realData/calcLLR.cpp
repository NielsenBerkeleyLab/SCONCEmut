#include <iostream>
#include <cmath>
#include <string>
#include <sstream>

/*
 * Thu 05 May 2022 12:38:59 PM PDT
 * script to calculate the log likelihood ratio testing if the allele frequency != 0 (that is, high LLR ==> variable site)
 *
 * g++ calcLLR.cpp -o calcLLR
 *
 * reads from stdin, prints to stdout
 *
 * expected line format:
 * <chr> <start> <end> <pooledMajorAllele> <cellMajorAllele> <numTotalReads> <numCellMinorReads>
 *
 * appends LLR to end of line and swaps numCellMinorReads to numPooledMajorReads for consistency
 * <chr> <start> <end> <pooledMajorAllele> <cellMajorAllele> <numTotalReads> <numPooledMajorReads> <LLR>
 */

int main() {
  // constants
  double errorRate = 0.005;
  double logErr = log(errorRate);
  double log1Err = log(1.0 - errorRate);
  double denomErr = 1.0 - 2.0 * errorRate;

  std::string chr;
  std::string start;
  std::string end;
  std::string pooledMajorAllele;
  std::string cellMajorAllele;
  int numTotalReads = 0;
  int numCellMinorReads = 0; // number of reads mapping to the minor allele (relative to one cell)
  int numPooledMajorReads = 0; // number of reads mapping to the major allele (relative to the pooled samples)
  double f_a = 0; // allele frequency for pooled major allele
  double llr = 0;

  // loop reading from https://stackoverflow.com/a/10464355
  for(std::string line; std::getline(std::cin, line);) {
    // parse input
    std::istringstream iss(line);
    iss >> chr >> start >> end >> pooledMajorAllele >> cellMajorAllele >> numTotalReads >> numCellMinorReads;

    // if pooled and cell specific major alleles match, then #maj = #total - #min
    if(pooledMajorAllele == cellMajorAllele) {
      numPooledMajorReads = numTotalReads - numCellMinorReads;
    }
    // else, don't match, so cell's minor is assumed to be pooled major allele, so read counts match
    else {
      numPooledMajorReads = numCellMinorReads;
    }

    // f_a = min{ (numPooledMajorReads - errorRate * numTotalReads) / ((1-2*errorRate) * numTotalReads), 1}
    // LLR = l(f_a) - l(f_a = 1)
    //     = numPooledMajorReads * (log(errorRate + f_a - 2*errorRate*f_a) - log(1-errorRate)) - (numTotalReads - numPooledMajorReads) * (log(errorRate) - log(1-errorRate - f_a + 2*errorRate*f_a))

    f_a = (numPooledMajorReads - errorRate * numTotalReads) / (denomErr * numTotalReads);
    if(f_a > 1) {
      f_a = 1;
    }
    else if(f_a < 0) {
      f_a = 0;
    }

    llr = numPooledMajorReads * (log(errorRate + f_a - 2.0 * errorRate * f_a) - log1Err) - (numTotalReads - numPooledMajorReads) * (logErr - log(1.0 - errorRate - f_a + 2.0 * errorRate * f_a));

    // print output
    // <chr> <start> <end> <pooledMajorAllele> <cellMajorAllele> <numTotalReads> <numPooledMajorReads> <LLR>
    printf("%s\t%s\t%s\t%s\t%s\t%i\t%i\t%0.5f\n", chr.c_str(), start.c_str(), end.c_str(), pooledMajorAllele.c_str(), cellMajorAllele.c_str(), numTotalReads, numPooledMajorReads, llr);
  }
  return 0;
}

