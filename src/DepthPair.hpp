#ifndef DEPTHPAIR_HPP
#define DEPTHPAIR_HPP

#include <string>
#include <numeric>
#include <vector>
#include <unordered_map>
#include <fstream>
#include <sstream>
#include <iostream>

#include <gsl/gsl_matrix.h>

#include <boost/filesystem.hpp> // for checking if files exist
#include <boost/algorithm/string.hpp> // for string splitting

#include "util.hpp" // for printMatrix

/*
 * This class stores a pair of tumor + mean diploid (windowed) read depth data
 */
class DepthPair {
  public:
   // member variables
    double diploidLibrarySize;
    double tumorLibrarySize;
    unsigned int numWindows;
    int maxWindowSize;
    double maxDiploidDepth;
    double maxTumorDepth;
    std::vector<std::string>* chrVec; // list of all chromosomes
    std::vector<std::string>* allRegions; // list of all regions
    std::unordered_map<std::string, std::vector<std::string>*>* regions; // chr:vector of region
    std::unordered_map<std::string, std::vector<double>*>* chrToDiploidDepthMap; // chr:vector of mean diploid depths
    std::unordered_map<std::string, std::vector<double>*>* chrToDiploidVarMap; // chr:vector of diploid variances
    std::unordered_map<std::string, std::vector<double>*>* chrToTumorDepthMap; // chr:vector of tumor depths
    std::unordered_map<std::string, std::vector<std::string>*>* chrToDiploidSimStateMap; // chr:vector of simulated states
    std::unordered_map<std::string, std::vector<std::string>*>* chrToTumorSimStateMap; // chr:vector of simulated states
    std::vector<double>* avgDiploidDepth;
    std::vector<double>* avgTumorDepth;

    // constructors and destructor
    DepthPair(std::string diploidFilename, std::string tumorFilename);
    DepthPair(DepthPair* otherDepths, std::string tumorFilename);
    DepthPair(int numWindows, int numChr = 1, int windowSize = 250000); // ctor for simulation
    DepthPair(bool shareDiploid, DepthPair* otherDepths); // ctor for simulation for subsequent cells
    DepthPair(std::string windowsFilename); // ctor for simulation with a windows file (ie if want to match a reference genome)
    DepthPair(const DepthPair& other);
    ~DepthPair();

     // functions
    void print(FILE* stream);
    double getDiploidDepthAt(std::string chr, int i);
    double getTumorDepthAt(std::string chr, int i);
    std::vector<double>* getAverageDepth(bool getDiploid);
    double getTotalDiploidDepth();
    double getTotalTumorDepth();
    void saveSimDiploidAsDepthFile(std::string filename);
    void saveSimTumorAsDepthFile(std::string filename);
    std::pair<std::unordered_map<std::string, std::vector<long long int>*>*, std::unordered_map<std::string, std::vector<int>*>*>* createChrWindowIdxLineNumMaps();
};

#endif

