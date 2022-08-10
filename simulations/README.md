# sconce_sim
This file describes the simulation program for SCONCEmut.

## Compiling
To compile the program, run `make`. All dependencies are included here.

## Generating simulated datasets
To run it, run `./SCONCEmut_sim <infile> <paramfile> [optional random seed]`. Or, first create an `infile1.txt` and `paramfile1.txt`, then run `./setupMutSims.sh infile1.txt paramfile1.txt 0`.

This will create the following files in the current working directory:
- `simu_readDepth_cancer_cell_*` (observed read depth for cancer cells)
- `simu_snpAlleles_cancer_cell_*` (observed number of reads mapping to each allele for cancer cells)
- `simu_readDepth_healthy_cell_*` (observed read depth for healthy cells)
- `true_copyNumber_cancer_cell_*` (true copy number for cancer cells)
- `true_snpAlleles_cancer_cell_*` (true number of copies for each allele for cancer cells)
- `true_copyNumber_healthy_cell_*` (true copy number for healthy cells)


## Parameter file formats
Any lines in the infile and paramfile after the arguments are ignored (ie good for comments).


The `<infile>` format is:
```
<flag for linesegment (0) model> <genome length> <# tumor cells> <flag for coalescence simulations (0)> <# healthy cells>
<rate of deletion> <rate of amplification> <mean deletion length> <mean amplification length> <point mutation rate>
<length of edge leading to root of tree = length of simulation time for single cell>
<tree in Newick format>
```

All published simulation sets use the same `<paramfile>` format. The `<paramfile>` format is:
```
<# bins to divide the genome into> <parameter r of the negative binomial (must be integer, converges to Poisson as r goes to infinity> <total expected number of reads (note that the observed number is random)> <length of a read, as fraction of the length of the genome> <sequencing error rate> <flag for even coverage in expectation before accounting for CNAs>
```

## Example parameter files
Parameter files for SCONCEmut are included in the inputFiles directory for reproducibility.

Newick tree strings were generated using createNewickTrees.R.

