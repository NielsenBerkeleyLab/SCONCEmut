# SCONCEmut

This is the repository for SCONCEmut.

## Dependencies
SCONCEmut is written in C++11 and requires the following libraries
- GNU make (tested on v4.1 and v4.2.1)
- g++ (tested on 7.5.0 and v9.4.0)
- Boost C++ Libraries (tested on v1.66 and v1.71; requires >= v1.66)
- GNU Scientific Library (tested on v2.4 and v2.5)

Additional [scripts](scripts/) require
- R (tested on v4.2.1)
- ape (tested on v5.6-2)
- cowplot (tested on v1.1.1)
- ggplot2 (tested on v3.3.6)
- ggtree (tested on v3.4.1)
- grid (tested on v0.5-1)
- gtools (tested on v3.9.2)
- phangorn (tested on v2.9.0)
- plyr (tested on v1.8.7)
- reshape2 (tested on v1.4.4)
- scales (tested on v1.2.0)
- stringr (tested on v1.4.0)
- python (tested on v3.10.5)

SCONCEmut was developed and tested on Ubuntu 20.04.4.

## Installation instructions
1. Clone this repo:
```
git clone git@github.com:NielsenBerkeleyLab/SCONCEmut.git
```
2. Run `make`. This will build intermediates into the `build/` directory and create an executable named `SCONCEmut`.

## Brief parameter descriptions for SCONCE2
- `--diploid` This file should be the averaged read depth across all diploid cells. It can be generated using [scripts/avgDiploid.R](scripts/avgDiploid.R).
- `--tumorFileList` This file should list filenames of the per window read depth of tumor cells to analyze.
- `--mutationFileList` This file should list filenames of the cell specific mutation read data, in the same order as `tumorFileList`
- `--meanVarCoefFile` This file should define the coefficients for the relationship between the mean and variance of the negative binomial distribution used for emission probabilities. It should be generated using [scripts/fitMeanVarRlnshp.R](scripts/fitMeanVarRlnshp.R).
- `--outputBase` This gives the basename for all [output files](#output-files).
- `--maxKploid` This gives the maximum allowed ploidy (recommended `k=10`).
- `-j` This gives how many threads sconce2 should use (default 1)
- `--sconceEstimatesPath` This is useful if one needs to stop and resume an analysis, or if one has previously run sconce and would like to skip the first sconce step. This should be the `--outputBase` value of the previous run. If previous runs cannot be found, sconce2 will be run.
- `--pairedEstimatesPath` as above, but for the paired estimates. Should be the `--outputBase` value of the previous run.
- Run `./SCONCEmut -h` for the full list of options


## Input files
Averaged diploid read depth files should be tab separated, with columns `<chr>\t<start>\t<end>\t<meanReadDepth>\t<varianceOfReadDepth>`. They should be generated using [scripts/avgDiploid.R](scripts/avgDiploid.R), given a file providing a list of paths to the observed diploid read depths. For example:
```
Rscript scripts/avgDiploid.R test/diploidFileList test/test_healthy_avg.bed
```
produces the following output:
```
$ head test/ref_healthy_avg.bed
chr1	0	250000	314.96	2413.43272727273
chr1	250000	500000	318.9	2269.84848484848
chr1	500000	750000	322.03	2715.32232323232
chr1	750000	1000000	321.58	1781.64
chr1	1000000	1250000	318.57	2367.21727272727
chr1	1250000	1500000	318.14	2574.24282828283
chr1	1500000	1750000	327.14	2439.11151515152
chr1	1750000	2000000	313.13	1945.95262626263
chr1	2000000	2250000	329.12	2260.00565656566
chr1	2250000	2500000	321.15	2663.78535353535
```

`tumorFileList` should list tumor read depth files to analyze, one per line. Tumor read depth files should be tab separated, with columns `<chr>\t<start>\t<end>\t<readDepth>`. See [simulations/README.md](simulations/README.md) for how to generate simulations with this format. For real data, a tool like [bedtools coverage](https://bedtools.readthedocs.io/en/latest/content/tools/coverage.html) can be used to create this file from a bam file. For example:
```
$ cat test/tumorFileList
test/simu_readDepth_cancer_cell_0.hg19_lite.bed
test/simu_readDepth_cancer_cell_1.hg19_lite.bed
test/simu_readDepth_cancer_cell_2.hg19_lite.bed
test/simu_readDepth_cancer_cell_3.hg19_lite.bed

$ head test/simu_cancer_cell_0.hg19_lite.bed
chr1	0	250000	348
chr1	250000	500000	401
chr1	500000	750000	307
chr1	750000	1000000	349
chr1	1000000	1250000	259
chr1	1250000	1500000	337
chr1	1500000	1750000	300
chr1	1750000	2000000	295
chr1	2000000	2250000	362
chr1	2250000	2500000	350
```

`mutationFileList` should list tumor mutation files to analyze, one per line, in the same order as `tumorFileList`. Tumor mutation files should be tab separated, with columns `<chr>\t<start>\t<end>\t<numReadsMappingToAncestralAllele>\t<numReadsMappingToDerivedAllele>`. See [simulations/setupMutSims.sh](simulations/setupMutSims.sh) for how to generate simulations with this format. See [scripts/realData/fullPipeline.sh](scripts/realData/fullPipeline.sh) for an example of analyzing real data. For example:
```
$ cat test/tumorMutList
test/simu_snpAlleles_cell_0.hg19_lite.bed
test/simu_snpAlleles_cell_1.hg19_lite.bed
test/simu_snpAlleles_cell_2.hg19_lite.bed
test/simu_snpAlleles_cell_3.hg19_lite.bed

$ head test/simu_cancer_cell_0.hg19_lite.bed
chr1	250001	250002	96	0
chr1	250002	250003	117	0
chr1	250003	250004	128	0
chr1	250004	250005	127	0
chr1	250005	250006	107	0
chr1	250006	250007	113	0
chr1	250007	250008	132	0
chr1	250008	250009	98	0
chr1	250010	250011	118	0
chr1	1250048	1250049	80	71
chr1	2250073	2250074	66	61
```


Mean and variance coefficient files should have one parameter and value pair per line. They should be generated using [scripts/fitMeanVarRlnshp.R](scripts/fitMeanVarRlnshp.R). For example:
```
Rscript scripts/fitMeanVarRlnshp.R test/diploidFileList test/ref.meanVar
```
produces the following output:
```
$ cat test/ref.meanVar
intercept=10.79950512968
slope= 1.18499408815
poly2= 0.01910756218
```

## Output files
SCONCE2 will create the following files automatically:
- `<output>.hmm` This file contains the state of the HMM after the Baum Welch step and the state of the HMM after the BFGS step.
- `<output>__<cellName>__k<maxKploid>.sconceParams` This file contains the model estimates from running sconce on this individual cell. Useful for resuming an analysis.
- `<output>__<cellName>__k<maxKploid>.sconceMutParams` This file contains the model estimates for the overdispersion (omega) estimate, the CNA to point mutation rate scaler (varphi), and the estimated number of somatic mutations for this cell. Useful for resuming an analysis.
- `<output>__<cell0Name>__<cell1Name>__k<maxKploid>.pairedParams` This file contains the model estimates from running cell0 and cell1 as a pair. Useful for resuming an analysis.
- `<output>__sconce__<cellName>__k<maxKploid>.bed` This file contains the copy number calls in tab separated bed format for this individual cell from running sconce. Only saved if the `--saveSconce` option is passed.
- `<output>__pair_<cell0Name>_<cell1Name>__<cellXName>__k<maxKploid>.bed` This file contains the copy number calls in tab separated bed format for cellX (out of the pair cell0 and cell1), created from the joint/paired analysis of cell0 and cell1.
- `<output>__<cellName>__k<maxKploid>__<summaryMethod>.bed` This file contains the summary copy number calls in tab separated bed format for this indivdual cell, created using <summaryMethod>. Only mean is saved by default, median and mode can be enabled by passing `--summarizeAll` or `--summarizeMedian`/`--summarizeMode`, respectively.


SCONCE2 also prints log messages to stdout and error messages to stderr.
If the `--debug` flag is used, debugging statements will be printed to stderr.

## Simulations
To compile and run the simulation program, see [simulations/README.md](simulations/README.md).


