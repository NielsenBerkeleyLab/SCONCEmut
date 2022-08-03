#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <vector>
#include <boost/program_options.hpp>
#include "HMM.hpp"
#include "DepthPair.hpp"
#include "util.hpp"
#include "Optimizable.hpp"
#include "MutationList.hpp"

#include "IndThenPairs2Stages3TrParam2DegPolyHMM.hpp"
#include "IndThenPairs2Stages3TrParam2DegPolyHMMWithMuts.hpp"

namespace po = boost::program_options; // see https://www.boost.org/doc/libs/1_71_0/doc/html/program_options/tutorial.html and https://www.boost.org/doc/libs/1_71_0/libs/program_options/example/options_description.cpp

/*
 * Tue 29 Mar 2022 06:13:12 PM PDT
 * This file is a main function to run sconce to estimate library sizes, beta, lambda, and tree branch lengths independently for all cells.
 * Then, lib sizes are fixed, beta and lambda are fixed at the median of all of their estimates. Also, we estimate mu (rate connecting CNAs and somatic point mutations) from t and X (est of mutation count, from read data).
 * Then, cells are paired up, and we get ests of X1/X2/X3 from read data, which are then fixed.
 * Then, BFGS is used to estimate tree branch lengths across all pairs of cells. The likelihood function is a function of t1/t2/t3, and includes CNAs and somatic point mutations
 *
 * For a cell pair (A,B), starting points for branch lengths are:
 *   t1 = min(t_A, t_B) / 2
 *   t2 = t_A
 *   t3 = t_B
 *
 * As of Thu 19 Aug 2021 05:18:20 PM PDT, used allPairs2Stages.cpp as a base for this file
 */

int main(int argc, char** argv) {
  // DEBUGGING: TURN OFF BUFFERING. from https://stackoverflow.com/a/1716621
  setbuf(stdout, NULL);

  // parse command line options
  int maxKploid = 0;
  std::string diploidFile;
  std::string tumorFileList;
  std::string mutationFileList;
  int numPairsStage2 = -1;
  bool estStage2Libs = false;
  bool disableIndBFGS = false;
  bool runIndBFGS = true;
  bool forwardDiff = false;
  int numThreads = 1;

  // output options
  std::string outputBase;
  bool saveSconce = false;
  bool saveVitDec = false;
  bool verbose = false;
  bool debug = false;

  // summarizing options
  int numPairsToSummarize = -1;
  int summarizeOrdering = -1;
  unsigned int seed;
  bool summarizeNearest = true;
  bool disableSummarizeNearest = false;
  bool summarizeRandom = false;
  bool summarizeFurthest = false;
  bool summarizeMean = true;
  bool disableSummarizeMean = false;
  bool summarizeMedian = false;
  bool summarizeMode = false;
  bool summarizeAll = false;

  // optimization start points
  int numBWIters = 0;
  int numLibStarts = 0;
  double libStartVal = 0;
  int maxIndBFGSIters = 0;
  int maxPairsBFGSIters = 0;

  // run options
  bool preallocIntermediates = false; // preallocate all intermediates ahead of time if true; allocates as needed if false
  bool readSconceEstsFromFile = false; // should we read sconce estimates from file? if false, will estimate them from scratch
  bool readSummarizedSconceEstsFromFile = false; // should we read summarized sconce estimates from file? precludes estimating nearest cell, since that depends on having indv estimates
  bool readPairedEstsFromFile = false; // should we read paired estimates from files?
  std::string sconceEstimatesPath; // key used to find sconce estimates per indv cell (ie if previously ran sconce or just some cells here, and want to just run paired stuff). uses sample names (from tumorFileList)
  std::string sconceEstimatesJointFile; // path to file containing col vector of all sconce params (ie if rerunning just paired stuff, use these params to skip over one cell sconce)
  std::string pairedEstimatesPath;

  // config file options if want to specify mean/var coefs for negative binomial. Otherwise, uses default from HMM::createMeanVarianceCoefVec()
  std::string meanVarCoefFile;
  double intercept = 0;
  double slope = 0;
  double poly2 = 0;
  gsl_vector* meanVarianceCoefVec = nullptr;
  try {
    po::options_description cmdLineOps("Command line options");
    cmdLineOps.add_options()
      ("diploid,d", po::value<std::string>(&diploidFile)->required(), "path to diploid depth file")
      ("tumorFileList,t", po::value<std::string>(&tumorFileList)->required(), "path to file listing tumor depth files")
      ("mutationFileList,m", po::value<std::string>(&mutationFileList)->required(), "path to file listing mutation files, assumed to be in the same order as tumorFileList")
      ("meanVarCoefFile", po::value<std::string>(&meanVarCoefFile), "path to negative binomial mean/variance coefficients file")
      ("outputBase,o", po::value<std::string>(&outputBase)->required(), "path to output files")
      ("maxKploid,k", po::value<int>(&maxKploid)->required(), "maximum allowed ploidy")
      ("numThreads,j", po::value<int>(&numThreads)->default_value(1), "number of threads to use")
      ("forwardDiff", po::bool_switch(&forwardDiff)->default_value(false), "use one-sided forward finite difference to approximate the gradient (uses two-sided central difference by default)")
      ("verbose,v", po::bool_switch(&verbose)->default_value(false), "enable verbose progress statements throughout the program")
      ("debug", po::bool_switch(&debug)->default_value(false), "enable debugging statements for gradient calculations in BFGS")
      ("saveSconce", po::bool_switch(&saveSconce)->default_value(false), "enable saving CNA calls after sconce has been run independently for each cell")
      ("saveViterbiDecoded", po::bool_switch(&saveVitDec)->default_value(false), "enable saving CNA calls in more verbose viterbiDecoded format")
      ("disableSummarizeMean", po::bool_switch(&disableSummarizeMean)->default_value(false), "disable summarizing CNA calls using mean across pairs")
      ("summarizeMedian", po::bool_switch(&summarizeMedian)->default_value(false), "summarize CNA calls using median across pairs")
      ("summarizeMode", po::bool_switch(&summarizeMode)->default_value(false), "summarize CNA calls using mode across pairs")
      ("summarizeAll", po::bool_switch(&summarizeAll)->default_value(false), "summarize CNA calls using all summary methods across pairs")
      ("numPairsStage2", po::value<int>(&numPairsStage2), "number of sequential pairs to run in stage 2. n choose 2 if not specified")
      ("numPairsToSummarize", po::value<int>(&numPairsToSummarize)->default_value(0), "number of pairs to summarize across per cell in stage 2. 0 for all pairs")
      ("disableSummarizeNearest", po::bool_switch(&disableSummarizeNearest)->default_value(false), "disable summarizing across pairs that are nearest to this cell")
      ("summarizeRandom", po::bool_switch(&summarizeRandom)->default_value(false), "summarize across pairs in a random ordering")
      ("summarizeFurthest", po::bool_switch(&summarizeFurthest)->default_value(false), "summarize across pairs that are furthest from this cell")
      ("seed", po::value<unsigned int>(&seed)->default_value(time(NULL)), "random seed, used if summarizing pairs randomly. uses current time if not passed")
      ("disableIndBFGS", po::bool_switch(&disableIndBFGS)->default_value(false), "disable BFGS on independent cells (in addition to baum welch) before optimizing over pairs of cells")

      ("estStage2Libs", po::bool_switch(&estStage2Libs)->default_value(false), "enable library size re-estimation (stage 2: all pairs). Note, this means parallelization is impossible")
      ("bwIters", po::value<int>(&numBWIters)->default_value(20), "number of baum welch iterations")
      ("maxIndBFGSIters", po::value<int>(&maxIndBFGSIters)->default_value(500), "maximum number of BFGS iterations for independent cell optim")
      ("maxPairsBFGSIters", po::value<int>(&maxPairsBFGSIters)->default_value(500), "maximum number of BFGS iterations for paired cell optim")
      ("numLibStarts", po::value<int>(&numLibStarts)->default_value(3), "number of library starting points for baum welch (should be 1 for [libStartVal], 2 for [1, 2], 3 for [1, 2, 4], or 4 for [1, 2, 4, 2/3]). These specify the multipliers for the lib estimate at the end of the first round of BW")
      ("libStartVal", po::value<double>(&libStartVal)->default_value(1.0), "if numLibStarts == 1, the value to start the library size scaling factor at")
      ("disableLowMem", po::bool_switch(&preallocIntermediates)->default_value(false), "disable low memory usage (will preallocate every intermediate ahead of time if this option is specified)")
      ("sconceEstimatesPath", po::value<std::string>(&sconceEstimatesPath)->default_value(""), "file key to construct filenames in tumorFileList to find SCONCE estimates when stored per indv cell")
      ("sconceEstimatesJointFile", po::value<std::string>(&sconceEstimatesJointFile)->default_value(""), "path to find SCONCE estimates when stored in one joint file")
      ("pairedEstimatesPath", po::value<std::string>(&pairedEstimatesPath)->default_value(""), "file key to construct filenames in tumorFileList to find paired parameter estimates")
      ("help,h", "help message")
      ;

    po::options_description meanVarCoefFileOps("Negative Binomial Mean and Variance Coefficient configuration file options");
    meanVarCoefFileOps.add_options()
      ("intercept", po::value<double>(&intercept)->required(), "var = {intercept} + slope * mean + poly2 * mean^2")
      ("slope", po::value<double>(&slope)->required(), "var = intercept + {slope} * mean + poly2 * mean^2")
      ("poly2", po::value<double>(&poly2)->required(), "var = intercept + slope * mean + {poly2} * mean^2");

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, cmdLineOps), vm);

    if(vm.count("help") || argc == 1) {
      std::cout << cmdLineOps << std::endl;;
      std::cout << meanVarCoefFileOps << std::endl;;
      return 0;
    }

    po::notify(vm); // deal with any command line errors (ie any missing required values)

    // check for conflicting options
    conflicting_options(vm, "summarizeFurthest", "summarizeRandom");
    conflicting_options(vm, "saveSconce", "sconceEstimatesJointFile");
    conflicting_options(vm, "sconceEstimatesPath", "sconceEstimatesJointFile");
    conflicting_options(vm, "numPairsToSummarize", "sconceEstimatesJointFile");

    if(vm.count("meanVarCoefFile")) {
      std::ifstream cfg(meanVarCoefFile.c_str());
      if(!cfg) {
        std::cerr << "Error: cannot open " << meanVarCoefFile << std::endl;
        exit(EXIT_FAILURE);
      }
      po::store(po::parse_config_file(cfg, meanVarCoefFileOps, true), vm); // true is to allow unknown options (ie for alpha to be in the masterParams file without messing up param parsing; see https://stackoverflow.com/a/31922646
      cfg.close();
      meanVarianceCoefVec = gsl_vector_alloc(3);
    }
    po::notify(vm); // store params from file into local vars, then store into meanVarianceCoefVec
    if(meanVarianceCoefVec != nullptr) {
      gsl_vector_set(meanVarianceCoefVec, 0, intercept);
      gsl_vector_set(meanVarianceCoefVec, 1, slope);
      gsl_vector_set(meanVarianceCoefVec, 2, poly2);
    }
  }
  catch(std::exception& e)
  {
    std::cerr << "Error: " << e.what() << "\n";
    exit(EXIT_FAILURE);
  }
  if(disableIndBFGS) {
    runIndBFGS = false;
  }
  if(disableSummarizeMean) {
    summarizeMean = false;
  }
  if(disableSummarizeNearest || summarizeRandom || summarizeFurthest) {
    summarizeNearest = false;
  }
  if(summarizeNearest) {
    summarizeOrdering = IndThenPairs2Stages3TrParam2DegPolyHMM::ORDERING_NEAREST;
  }
  else if(summarizeRandom) {
    summarizeOrdering = IndThenPairs2Stages3TrParam2DegPolyHMM::ORDERING_RANDOM;
  }
  else if(summarizeFurthest) {
    summarizeOrdering = IndThenPairs2Stages3TrParam2DegPolyHMM::ORDERING_FURTHEST;
  }

  if(sconceEstimatesPath.length() > 0) {
    readSconceEstsFromFile = true;
  }
  if(sconceEstimatesJointFile.length() > 0) {
    if(!boost::filesystem::exists(sconceEstimatesJointFile)) {
      std::cerr << "Error: --sconceEstimatesJointFile " << sconceEstimatesJointFile << " does not exist. Exiting" << std::endl;
      exit(EXIT_FAILURE);
    }
    readSummarizedSconceEstsFromFile = true;
  }
  if(pairedEstimatesPath.length() > 0) {
    readPairedEstsFromFile = true;
    if(!readSconceEstsFromFile) {
      std::cerr << "Error: --pairedEstimatesPath passed but missing --sconceEstimatesPath. Exiting" << std::endl;
      exit(EXIT_FAILURE);
    }
  }

  if(!boost::filesystem::exists(diploidFile)) {
    std::cerr << "Error: --diploid " << diploidFile << " does not exist. Exiting" << std::endl;
    exit(EXIT_FAILURE);
  }
  if(!boost::filesystem::exists(tumorFileList)) {
    std::cerr << "Error: --tumorFileList " << tumorFileList << " does not exist. Exiting" << std::endl;
    exit(EXIT_FAILURE);
  }
  if(!boost::filesystem::exists(mutationFileList)) {
    std::cerr << "Error: --mutationFileList " << mutationFileList << " does not exist. Exiting" << std::endl;
    exit(EXIT_FAILURE);
  }

  // ##### set up AllPairs2Stages obj #####
  std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

  // read through tumorFileList, making DepthPairs
  std::vector<DepthPair*>* depthsVec = new std::vector<DepthPair*>();
  DepthPair* firstDepthPair = nullptr;
  std::ifstream tumorListStream(tumorFileList);
  std::string currTumorFile;
  std::vector<std::string>* sampleList = new std::vector<std::string>();
  while(getline(tumorListStream, currTumorFile)) {
    // if on the first one, set firstDepthPair
    if(firstDepthPair == nullptr) {
      firstDepthPair = new DepthPair(diploidFile, currTumorFile);
      depthsVec->push_back(firstDepthPair);
    }
    else {
      depthsVec->push_back(new DepthPair(firstDepthPair, currTumorFile));
    }
    sampleList->push_back(parseSampleName(currTumorFile));
  }
  tumorListStream.close();

  // read through mutationFileList
  std::pair<std::unordered_map<std::string, std::vector<long long int>*>*, std::unordered_map<std::string, std::vector<int>*>*>* chrWindowIdxLineNumMaps = firstDepthPair->createChrWindowIdxLineNumMaps();
  std::vector<MutationList*>* mutListVec = new std::vector<MutationList*>();
  std::ifstream mutListStream(mutationFileList);
  std::string currMutFile;
  while(getline(mutListStream, currMutFile)) {
    mutListVec->push_back(new MutationList(currMutFile, chrWindowIdxLineNumMaps));
  }
  mutListStream.close();
  // ##### end file reading #####

  // if wasn't passed a param for numPairsStage2 (used for debugging early out), then calculate n choose 2
  if(numPairsStage2 == -1) {
    numPairsStage2 = boost::math::binomial_coefficient<double>(depthsVec->size(), 2);
  }

  // create all paired HMMs of tumor cells
  IndThenPairs2Stages3TrParam2DegPolyHMMWithMuts* indThenPairsHMM = new IndThenPairs2Stages3TrParam2DegPolyHMMWithMuts(depthsVec, sampleList, mutListVec, maxKploid, runIndBFGS, numPairsStage2, verbose, debug, !forwardDiff); // flip forwardDiff flag to centralDiff
  indThenPairsHMM->setNumThreads(numThreads);
  gsl_vector* optimizedIndParamsToEst = nullptr;

  // if not reading paired estimates from file, then we have to start from the beginning
  FILE* oFile = nullptr;
  gsl_vector* stage2FixedParams = nullptr;
  gsl_vector* stage2InitGuess = nullptr;
  // ##### run independent cell estimation #####
  gsl_vector* indInitGuess = nullptr;
  std::cout << "######## STARTING INDEPENDENT CELL SET UP ########" << std::endl;
  indThenPairsHMM->setUpInd2Stages();
  if(meanVarianceCoefVec != nullptr) {
    indThenPairsHMM->setAllIndMeanVarianceFn(meanVarianceCoefVec);
  }

  // if reading sconce estimates from file
  if(readSconceEstsFromFile) {
    // read from indv cell files
    std::cout << "######## READING INDV CELL .sconceParams FILES AND COMBINING INTO optimizedIndParamsToEst FROM FILES WITH PATHS FROM " << tumorFileList << " WITH PATH " << sconceEstimatesPath << " ######## " << std::endl;
    indThenPairsHMM->getIndOptimParamsToEstFromIndvFiles(tumorFileList, sconceEstimatesPath, 5, meanVarianceCoefVec); // lib + alpha/beta/gamma + branch; this is per cell
  }

  // estimate lib size, beta, lambda, t for each cell individually, using only copy numbers
  indInitGuess = gsl_vector_alloc(indThenPairsHMM->getIndInitGuessSize());
  indThenPairsHMM->optimIndCells(indInitGuess, outputBase, numBWIters, numLibStarts, libStartVal, maxIndBFGSIters, verbose, debug);

  indThenPairsHMM->viterbiDecodeStage1();
  if(saveSconce) {
    if(saveVitDec) {
      indThenPairsHMM->saveStage1ViterbiDecodedCNA(outputBase);
    }
    indThenPairsHMM->saveStage1CNAToBed(outputBase);
    indThenPairsHMM->saveStage1ParamEstimates(outputBase);
  }

  // INSERT indInitGuess SHORTCUTS HERE
/*gsl_vector_set(indInitGuess,   0   , 0.9736357438666501940005559845303650945425);
gsl_vector_set(indInitGuess,   101 , 0.1226746711770437536781486187464906834066);*/

  printf("TOTAL LOGLIKLIHOOD AFTER INDEPENDENT CELL ESTIMATION: %.5f\n\n", indThenPairsHMM->getTotalIndLogLikelihood());
  //exit(0);

  // estimate mutation parameters on individual cells
  // mu: cnaToMutRateMu for l(t1,t2,t3) ~ P(cnaToMutRateMu * t1 * X1)...
  // omega: mutOverdispOmega for l(D_ij | S_ij) ~ BetaBinom(D_ij, f, mutOverdispOmega))
  std::chrono::steady_clock::time_point jointMutEstBegin = std::chrono::steady_clock::now();
  std::cout << "######## STARTING JOINT MUTATION PARAMETERS ESTIMATION ########" << std::endl;
  indThenPairsHMM->estimateMutParamsFromIndCells(sconceEstimatesPath, outputBase);
  std::chrono::steady_clock::time_point jointMutEstEnd = std::chrono::steady_clock::now();
  double jointMutEstElapsedSec = std::chrono::duration_cast<std::chrono::microseconds>(jointMutEstEnd - jointMutEstBegin).count() / 1000000.0;
  printf("JOINT MUTATION PARAMETERS ESTIMATION TOTAL TIME (sec): %.5f\n\n", jointMutEstElapsedSec);
  std::cout << "######## DONE WITH JOINT MUTATION PARAMETERS ESTIMATION ########" << std::endl;

  // write results of ind stage to file
  oFile = fopen((outputBase + ".hmm").c_str(), "w");
  indThenPairsHMM->print(oFile);
  fclose(oFile);

  // if passed a joint filename, read from joint file
  if(readSummarizedSconceEstsFromFile) {
    std::cout << "######## READING JOINT optimizedIndParamsToEst FROM " << sconceEstimatesJointFile << " ######## " << std::endl;
    int numExpectedLines = depthsVec->size() + 3 + depthsVec->size(); // libs + alpha/beta/lambda + branch lengths
    optimizedIndParamsToEst = indThenPairsHMM->getIndOptimParamsToEstSummaryFromJointFile(sconceEstimatesJointFile, numExpectedLines);
  }
  // else calc summary
  else {
    optimizedIndParamsToEst = indThenPairsHMM->getIndOptimParamsToEstSummary();
  }

  // ##### end independent cell estimation #####
  std::cout << "######## DONE WITH INDEPENDENT CELL ESTIMATION ########" << std::endl;
  if(verbose) {
    std::cout << "optimizedIndParamsToEst:" << std::endl;
    printColVector(optimizedIndParamsToEst);
    fprintf(stdout, "cnaToMutRateMu, mutOverdispOmega:\n%.40f\n%.40f\n\n", indThenPairsHMM->getCnaToMutRateMu(), indThenPairsHMM->getMutOverdispOmega());
    gsl_vector* stage1MutCounts = indThenPairsHMM->getAllIndEstMutCounts();
    std::cout << "stage1MutCounts:" << std::endl;
    printColVector(stage1MutCounts);
  }

  std::chrono::steady_clock::time_point stage1End = std::chrono::steady_clock::now();
  double stage1ElapsedSec = std::chrono::duration_cast<std::chrono::microseconds>(stage1End - begin).count() / 1000000.0;
  printf("STAGE 1 TOTAL TIME (sec): %.5f\n", stage1ElapsedSec);

  // ##### set up stage 2 #####
  std::chrono::steady_clock::time_point stage2setupBegin = std::chrono::steady_clock::now();
  std::cout << "######## STARTING STAGE 2 SET UP ########" << std::endl;
  // by default, do not estimate stage 2 libs, as this would remove the possibility of parallelizing. But estimate them if the user insists
  if(!estStage2Libs) {
    stage2FixedParams = gsl_vector_alloc(depthsVec->size() + 3 + 1); // all libs + alpha/beta/lambda are fixed + cnaToMutRateMu
    stage2InitGuess = gsl_vector_alloc(3 * numPairsStage2); // t1/t2/t3 for each pair are estimated
    indThenPairsHMM->copyStage1EstsIntoStage2FixedParams(stage2FixedParams, stage2InitGuess, optimizedIndParamsToEst, depthsVec->size());
  }
  else {
    stage2FixedParams = gsl_vector_alloc(3 + 1); // alpha/beta/lambda are fixed + cnaToMutRateMu
    stage2InitGuess = gsl_vector_alloc(depthsVec->size() + 3 * numPairsStage2); // all libs and t1/t2/t3 for each pair are estimated
    indThenPairsHMM->copyStage1EstsIntoStage2FixedParams(stage2FixedParams, stage2InitGuess, optimizedIndParamsToEst, 0); // save 0 libs into fixedParams
  }

  // INSERT optimizedIndParamsToEst AND stage2FixedParams SHORTCUTS HERE

  if(verbose) {
    std::cout << "optimizedIndParamsToEst" << std::endl;
    printColVector(optimizedIndParamsToEst);
    std::cout << "stage2FixedParams:" << std::endl;
    printColVector(stage2FixedParams);
    std::cout << "stage2InitGuess:" << std::endl;
    printColVector(stage2InitGuess);
  }

  if(numPairsToSummarize > 0 && !readSummarizedSconceEstsFromFile) { // readSummarizedSconceEstsFromFile and numPairsToSummarize conflict done earlier, extra safety check
    indThenPairsHMM->calcStage1PairwiseDistMat();
    indThenPairsHMM->calcStage1NearestCellIdxMat();
    indThenPairsHMM->setUpSelectPairsBFGS(stage2FixedParams, stage2InitGuess, numPairsToSummarize, summarizeOrdering, seed, meanVarianceCoefVec, !estStage2Libs, preallocIntermediates);
  }
  else {
    indThenPairsHMM->setUpAllPairsBFGS(stage2FixedParams, stage2InitGuess, meanVarianceCoefVec, !estStage2Libs, preallocIntermediates);
  }
  if(readPairedEstsFromFile) {
    std::cout << "######## READING STAGE 2 PAIRED ESTIMATES FROM FILES WITH PATH " << pairedEstimatesPath << " ########" << std::endl;
    indThenPairsHMM->getPairedOptimParamsToEstFromFiles(pairedEstimatesPath, 2 + 3 + 3); // 2 libs + alpha/beta/lambda + 3 branches per pair
    if(verbose) {
      std::cout << "stage2->paramsToEst" << std::endl;
      printColVector(indThenPairsHMM->getStage2HMM()->getParamsToEst());
      std::cout << std::endl;
      std::cout << "stage2->fixedParams" << std::endl;
      printColVector(indThenPairsHMM->getStage2HMM()->getFixedParams());
      std::cout << std::endl;
    }
    //indThenPairsHMM->getStage2HMM()->print(stdout);
    //exit(0);
  }

  // clean up stage 1 intermediates
  indThenPairsHMM->cleanUpIndOptim();

  std::chrono::steady_clock::time_point stage2setupEnd = std::chrono::steady_clock::now();
  double stage2setupElapsedSec = std::chrono::duration_cast<std::chrono::nanoseconds>(stage2setupEnd - stage2setupBegin).count() / 1e9;
  fprintf(stdout, "STAGE 2 SETUP TIME (sec): %.20f\n", stage2setupElapsedSec);

  // mutation counts are estimated when MutationPairs are created to save memory
  //// ##### estimate mutation counts for each MutationPair #####
  //std::chrono::steady_clock::time_point mutPairMutCountBegin = std::chrono::steady_clock::now();
  //std::cout << "######## STARTING MUTATIONPAIR:MUTATION COUNT BFGS  ########" << std::endl;
  //indThenPairsHMM->bfgsStage2(stage2InitGuess, maxPairsBFGSIters, outputBase, verbose, debug);
  //indThenPairsHMM->est
  //std::cout << "######## DONE WITH MUTATIONPAIR:MUTATION COUNT BFGS ########" << std::endl;
  //std::chrono::steady_clock::time_point mutPairMutCountEnd = std::chrono::steady_clock::now();
  //double mutPairMutCountElapsedSec = std::chrono::duration_cast<std::chrono::nanoseconds>(mutPairMutCountEnd - mutPairMutCountBegin).count() / 1e9;
  //fprintf(stdout, "MUTATIONPAIR:MUTATION COUNT BFGS TIME (sec): %.20f\n", mutPairMutCountElapsedSec);

  // ##### run bfgs on stage 2 #####
  std::cout << "######## STARTING STAGE 2 BFGS ########" << std::endl;
  indThenPairsHMM->bfgsStage2(stage2InitGuess, maxPairsBFGSIters, outputBase, verbose, debug);
  std::cout << "######## DONE WITH STAGE 2 BFGS ########" << std::endl;
  std::cout << "######## STARTING DECODING ########" << std::endl;

  // do viterbi decoding of all pairs
  std::chrono::steady_clock::time_point decodeStart = std::chrono::steady_clock::now();
  indThenPairsHMM->viterbiDecodeStage2();
  if(saveVitDec) {
    indThenPairsHMM->saveStage2ViterbiDecodedCNA(outputBase);
  }
  std::cout << "######## SAVING STAGE 2 PAIRED CNA TO BED (" << outputBase << ") ########" << std::endl;
  indThenPairsHMM->saveStage2CNAToBed(outputBase);
  std::cout << "######## DONE SAVING STAGE 2 PAIRED CNA TO BED ########" << std::endl;

  // get summary calls across pairs
  bool runDecoding = true; // don't need to run decoding again if using multiple summarizing methods, since viterbi decoding won't change
  if(summarizeAll || summarizeMean) {
    indThenPairsHMM->summaryDecodeAllCellsAcrossPairs(AllPairs3TrParam2DegPolyHMM::SUMMARY_MEAN, runDecoding);
    indThenPairsHMM->saveAllSummaryDecodedCNAToBed(outputBase);
    runDecoding = false;
  }
  if(summarizeAll || summarizeMedian) {
    indThenPairsHMM->summaryDecodeAllCellsAcrossPairs(AllPairs3TrParam2DegPolyHMM::SUMMARY_MEDIAN, runDecoding);
    indThenPairsHMM->saveAllSummaryDecodedCNAToBed(outputBase);
    runDecoding = false;
  }
  if(summarizeAll || summarizeMode) {
    indThenPairsHMM->summaryDecodeAllCellsAcrossPairs(AllPairs3TrParam2DegPolyHMM::SUMMARY_MODE, runDecoding);
    indThenPairsHMM->saveAllSummaryDecodedCNAToBed(outputBase);
    runDecoding = false;
  }

  std::chrono::steady_clock::time_point decodeEnd = std::chrono::steady_clock::now();
  double decodeElapsedSec = std::chrono::duration_cast<std::chrono::microseconds>(decodeEnd - decodeStart).count() / 1000000.0;
  printf("TOTAL DECODING TIME ELAPSED (sec): %.5f\n", decodeElapsedSec);
  printf("AVERAGE DECODING TIME ELAPSED (sec): %.5f\n", decodeElapsedSec / (double) numPairsStage2);
  std::cout << "######## DONE WITH DECODING ########" << std::endl;

  // #### save the optimized HMM ####
  std::cout << "######## SAVING HMM ########" << std::endl;
  oFile = fopen((outputBase + ".hmm").c_str(), "a");
  fprintf(oFile, "\n\n################################################################################\n");
  fprintf(oFile, "#### FINAL HMM ####\n");
  fprintf(oFile, "################################################################################\n\n");
  indThenPairsHMM->print(oFile);
  fclose(oFile);
  std::chrono::steady_clock::time_point saveHMMEnd = std::chrono::steady_clock::now();
  double saveHMMElapsedSec = std::chrono::duration_cast<std::chrono::microseconds>(saveHMMEnd - decodeEnd).count() / 1000000.0;
  printf("TOTAL SAVING FINAL HMM TIME ELAPSED (sec): %.5f\n", saveHMMElapsedSec);
  std::cout << "######## DONE WITH SAVING HMM ########" << std::endl;

  std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
  double elapsedSec = std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() / 1000000.0;
  printf("TOTAL TIME ELAPSED (sec): %.5f\n", elapsedSec);

  // clean up
  // TODO
  deleteChrWindowIdxLineNumMaps(chrWindowIdxLineNumMaps); // includes a call to delete for pair pointer

  return 0;
}


