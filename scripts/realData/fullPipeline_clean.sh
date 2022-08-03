#!/bin/bash

# clean version of fullPipeline.sh

export ref="/space/s2/sandra/hg19/hg19_lite.fa"

# prep
bamlist=/space/s2/sandra/alleleFreqHmm/Navin_Nature2011_hg19/diploidBamList
pileup=/space/s2/sandra/alleleFreqHmm/Navin_Nature2011_hg19/monovarOutput/diploid.mpileup.gz


####################
# find major allele pooled across all diploid cells
####################
# create pooled diploid pileup file
samtools mpileup -Q 0 -d 10000 -f "$ref" -q 40 -b "$bamlist" | gzip > "$pileup" # pooled mpileup file

# split pooled diploid pileup from full genome into pooled diploid chrs
parallel -j 10 zgrep -E -w chr{} "$pileup" > pooledPileupFiles/diploid_chr_{}.mpileup ::: `seq 22`

# find pooled diploid major allele (./calc_reads_parallel.sh) TODO missing chrX,Y since sometimes labelled chr22. just throw out for now
for i in $(seq 1 22) ; do
  echo "python3 find_major_allele.py -f pooledPileupFiles/diploid_chr_"$i".mpileup -t diploid -i "$i" > pooledMajorAlleleReadCounts/diploid_chr_"$i".err 2>&1"
done | parallel -j 10

# reformat calc_reads_parallel.sh output into bed file (./create_chromosome_bedfiles.sh)
for chr in $(seq 1 22) ; do
  tail -n +2 pooledMajorAlleleReadCounts/diploid_chr_"$chr"_read_allele_counts.tsv | awk '{print $1"\t"$2"\t"$2+1"\t"$3}' > pooledMajorAlleleBedFiles/diploid_chr"$chr".bed
done
find pooledMajorAlleleBedFiles -maxdepth 1 -name "diploid_chr[0-9]*bed" -exec cat {} \+ | sort -k1,1 -k2,2n > pooledMajorAlleleBedFiles/diploid_sorted.bed


####################
# for each individual diploid cell, count up the total number of reads and the number of reads that map to the pooled diploid major allele, and then calculate the LLR (checks if a site is variable, ie high LLR ==> variable), and put this all into a big matrix with a flag for dbsnp membership
####################
# create individual diploid pileup files
pileupDir=indvDiploidPileups
mkdir -p "$pileupDir"
while read bamFile ; do
  shortName=$(basename "$bamFile" .bam)
  samtools mpileup -Q 0 -d 10000 -f "$ref" -q 40 "$bamFile" > "$pileupDir""/""$shortName"".mpileup"
done < "$bamlist"

# for each indvidual diploid cell,
# - add pooled major allele to indv diploid cell read counts
# - calc LLR for each site (relative to pooled major allele) (calls ./calcLLR); later discarded
# - add "inDbsnp" "flag"
find ./indvDiploidPileups/ -maxdepth 1 -name "*mpileup" -exec basename {} .q20.adaptTr.qualTr.dusted.sorted.noDups.mpileup \; | sort | parallel -j 10 --delay 0.5 ./createDiploidLLRdbsnpBedFiles.sh {} "indvDiploidPileups/" "dbsnpIntersection"

# combine cell specific bed files into one big matrix (dbsnpIntersection/pooledDiploidMajorAllele.llr.dbsnp.bed)
./createPooledDiploidLlrMatrix.sh

# summarize llr across diploid cells
awk '{cellCount=0; pooledReadCount=0; pooledMajorReadCount=0; summedLlr=0;
  for(i=5; i <= NF; i+=3) {
    if($i != "NA") {
      cellCount+=1
      pooledReadCount+=$i
      pooledMajorReadCount+=$(i+1)
      summedLlr+=$(i+2)
    }
  } print $1"\t"$2"\t"$3"\t"$4"\t"cellCount"\t"pooledReadCount"\t"pooledReadCount/cellCount"\t"pooledMajorReadCount"\t"summedLlr}' dbsnpIntersection/pooledDiploidMajorAllele.llr.dbsnp.bed > dbsnpIntersection/pooledDiploidMajorAllele.llr.dbsnp.countLlrSummary.bed

# find 95% of sites in dbSNP llr cutoff
# g++ countLLRcutoff.cpp -o countLLRcutoff
tail -n +2 dbsnpIntersection/pooledDiploidMajorAllele.llr.dbsnp.countLlrSummary.bed | ./countLLRcutoff > dbsnpIntersection/pooledDiploidMajorAllele.llr.dbsnp.countLlrSummary.cutoffCounts.txt
Rscript plotCutoffCounts.R dbsnpIntersection/pooledDiploidMajorAllele.llr.dbsnp.countLlrSummary.cutoffCounts.txt
# [1] "95% llr cutoff is: 33"
# ==> 32 or lower is probably error if not in dbsnp

awk '{if($9 <= 32 && $4 == ".") print}' dbsnpIntersection/pooledDiploidMajorAllele.llr.dbsnp.countLlrSummary.bed > dbsnpIntersection/pooledDiploidMajorAllele.llr.dbsnp.countLlrSummary.llr32.notInDbsnp.bed
# ==> coords of sites that are probably errors

# calculate avg num cells/reads for these sites to set minimum cell and read average cutoffs
awk '{cellSum += $5; readSum += $7} END{print cellSum / NR"\t"readSum / NR}' dbsnpIntersection/pooledDiploidMajorAllele.llr.dbsnp.countLlrSummary.llr32.notInDbsnp.bed
# 2.47439 avg num cells; 2*2.47 ~= 5
# 1.75195 avg num reads; 2*1.75 ~= 3.5
# 869878372 total sites
minCells=5
minAvgReads=3.5
awk -v minCells=$minCells minAvgReads=$minAvgReads'{if($5 >= minCells && $7 >= minAvgReads) {print}}' dbsnpIntersection/pooledDiploidMajorAllele.llr.dbsnp.countLlrSummary.llr32.notInDbsnp.bed > dbsnpIntersection/pooledDiploidMajorAllele.llr.dbsnp.countLlrSummary.llr32.notInDbsnp.minCells"$minCells"_minAvgReads"$minAvgReads".bed # think this might end up being essentially the same as dbsnpIntersection/pooledDiploidMajorAllele.llr.dbsnp.notInDbsnp.countSummary.minCells"$minCells"_minAvgReads"$minAvgReads".bed (ie don't think llr 32 filter did much) ==> 3959722 sites
bedtools intersect -b dbsnpIntersection/pooledDiploidMajorAllele.llr.dbsnp.countLlrSummary.llr32.notInDbsnp.minCells"$minCells"_minAvgReads"$minAvgReads".bed -a dbsnpIntersection/pooledDiploidMajorAllele.llr.dbsnp.bed -sorted > dbsnpIntersection/pooledDiploidMajorAllele.llr.dbsnp.llr32.notInDbsnp.minCells"$minCells"_minAvgReads"$minAvgReads".bed


####################
# estimate overdispersion parameter, omega, from filtered diploid data
####################
Rscript fitBetaBinomNull.R dbsnpIntersection/pooledDiploidMajorAllele.llr.dbsnp.llr32.notInDbsnp.minCells"$minCells"_minAvgReads"$minAvgReads".bed dbsnpIntersection/pooledDiploidMajorAllele.llr.dbsnp.header
# 8.68008260176282497866


####################
# preprocess tumor data: generate mpileups, then convert into a matrix of coords, total read counts, and major allele read counts
####################
bamlist=/space/s2/sandra/alleleFreqHmm/Navin_Nature2011_hg19/tumorBamList

# create individual diploid pileup files
pileupDir=indvDiploidPileups
mkdir -p "$pileupDir"
while read bamFile ; do
  shortName=$(basename "$bamFile" .bam)
  samtools mpileup -Q 0 -d 10000 -f "$ref" -q 40 "$bamFile" > "$pileupDir""/""$shortName"".mpileup"
done < "$bamlist"

find indvTumorPileups/ -maxdepth 1 -name "*mpileup" -exec basename {} .q20.adaptTr.qualTr.dusted.sorted.noDups.mpileup \; | sort | parallel -j 8 ./createDiploidLLRdbsnpBedFiles.sh {} "indvTumorPileups/" "dbsnpIntersection_tumor"
./createPooledTumorLLRmatrix.sh # creates dbsnpIntersection_tumor/pooledDiploidTumorMajorAllele.llr.dbsnp.bed


####################
# using omega from diploid data, calc maximum likelihood allele freqs for tumor sites
####################
cat dbsnpIntersection_tumor/pooledDiploidTumorMajorAllele.llr.dbsnp.countSummary.minCells5_minAvgReads3.5.readsMat.bed | ./fitBetaBinomTumor > dbsnpIntersection_tumor/pooledDiploidTumorMajorAllele.llr.dbsnp.countSummary.minCells5_minAvgReads3.5.readsMat.pooledLlr.bed
# $ wc -l pooledDiploidTumorMajorAllele.llr.dbsnp.countSummary.minCells5_minAvgReads3.5.readsMat.bed
# 8966211 pooledDiploidTumorMajorAllele.llr.dbsnp.countSummary.minCells5_minAvgReads3.5.readsMat.bed
# $ wc -l pooledDiploidTumorMajorAllele.llr.dbsnp.countSummary.minCells5_minAvgReads3.5.readsMat.pooledLlr.bed # any non variable (ie only see all major alleles) sites are removed
# 433947 pooledDiploidTumorMajorAllele.llr.dbsnp.countSummary.minCells5_minAvgReads3.5.readsMat.pooledLlr.bed

# then count how many tumor sites survive each llr cutoff
for llr in $(seq 1 15) ; do
  awk -v llr=$llr '{if($6 > llr) print}' dbsnpIntersection_tumor/pooledDiploidTumorMajorAllele.llr.dbsnp.countSummary.minCells5_minAvgReads3.5.readsMat.pooledLlr.bed > dbsnpIntersection_tumor/pooledDiploidTumorMajorAllele.llr.dbsnp.countSummary.minCells5_minAvgReads3.5.readsMat.pooledLlr.llr"$llr".bed
done
#$ wc -l dbsnpIntersection_tumor/pooledDiploidTumorMajorAllele.llr.dbsnp.countSummary.minCells5_minAvgReads3.5.readsMat.pooledLlr.llr*.bed
#  416847 dbsnpIntersection_tumor/pooledDiploidTumorMajorAllele.llr.dbsnp.countSummary.minCells5_minAvgReads3.5.readsMat.pooledLlr.llr1.bed
#   52853 dbsnpIntersection_tumor/pooledDiploidTumorMajorAllele.llr.dbsnp.countSummary.minCells5_minAvgReads3.5.readsMat.pooledLlr.llr2.bed
#   39180 dbsnpIntersection_tumor/pooledDiploidTumorMajorAllele.llr.dbsnp.countSummary.minCells5_minAvgReads3.5.readsMat.pooledLlr.llr3.bed
#   15410 dbsnpIntersection_tumor/pooledDiploidTumorMajorAllele.llr.dbsnp.countSummary.minCells5_minAvgReads3.5.readsMat.pooledLlr.llr4.bed ==> use 4 as cutoff? results in 4.15254 somatic sites on avg
#    4800 dbsnpIntersection_tumor/pooledDiploidTumorMajorAllele.llr.dbsnp.countSummary.minCells5_minAvgReads3.5.readsMat.pooledLlr.llr5.bed
#    3121 dbsnpIntersection_tumor/pooledDiploidTumorMajorAllele.llr.dbsnp.countSummary.minCells5_minAvgReads3.5.readsMat.pooledLlr.llr6.bed
#    1840 dbsnpIntersection_tumor/pooledDiploidTumorMajorAllele.llr.dbsnp.countSummary.minCells5_minAvgReads3.5.readsMat.pooledLlr.llr7.bed
#    1116 dbsnpIntersection_tumor/pooledDiploidTumorMajorAllele.llr.dbsnp.countSummary.minCells5_minAvgReads3.5.readsMat.pooledLlr.llr8.bed
#     780 dbsnpIntersection_tumor/pooledDiploidTumorMajorAllele.llr.dbsnp.countSummary.minCells5_minAvgReads3.5.readsMat.pooledLlr.llr9.bed
#     558 dbsnpIntersection_tumor/pooledDiploidTumorMajorAllele.llr.dbsnp.countSummary.minCells5_minAvgReads3.5.readsMat.pooledLlr.llr10.bed
#     426 dbsnpIntersection_tumor/pooledDiploidTumorMajorAllele.llr.dbsnp.countSummary.minCells5_minAvgReads3.5.readsMat.pooledLlr.llr11.bed
#     322 dbsnpIntersection_tumor/pooledDiploidTumorMajorAllele.llr.dbsnp.countSummary.minCells5_minAvgReads3.5.readsMat.pooledLlr.llr12.bed
#     247 dbsnpIntersection_tumor/pooledDiploidTumorMajorAllele.llr.dbsnp.countSummary.minCells5_minAvgReads3.5.readsMat.pooledLlr.llr13.bed
#     194 dbsnpIntersection_tumor/pooledDiploidTumorMajorAllele.llr.dbsnp.countSummary.minCells5_minAvgReads3.5.readsMat.pooledLlr.llr14.bed
#     155 dbsnpIntersection_tumor/pooledDiploidTumorMajorAllele.llr.dbsnp.countSummary.minCells5_minAvgReads3.5.readsMat.pooledLlr.llr15.bed
#  537849 total

for bed in $(find dbsnpIntersection_tumor -maxdepth 1 -name "SRR*_withPooledMajorAllele.llr.dbsnp.bed.gz" | sort -V ) ; do
  #for llr in "3" "4" ; do
  for llr in "2" ; do
    tumorCoordsOutfile=$(echo $bed | sed "s/.bed.gz/.minCells5_minAvgReads3.5.readsMat.llr_t"$llr".bed/")
    # bed files have: chr start end pooledMajorAllele cellMajorAllele totalNumReads numPooledMajorAlleleReads LLR 
    # want to get: chr start end numAnc numDer
    zcat $bed | bedtools intersect -a - -b dbsnpIntersection_tumor/pooledDiploidTumorMajorAllele.llr.dbsnp.countSummary.minCells5_minAvgReads3.5.readsMat.pooledLlr.llr"$llr".bed -wa | awk '{print $1"\t"$2"\t"$3"\t"$7"\t"($6-$7)}' > $tumorCoordsOutfile
  done
done

# avg num derived sites vs avg num total sites; both approx 10/3136 ~= 4/1215 ~= 0.003 frac derived sites
#$ find . -name "SRR*_withPooledMajorAllele.llr.dbsnp.minCells5_minAvgReads3.5.readsMat.llr_t3.bed" -exec bash -c 'grep -v "0$" {} | wc -l ' \; | awk '{sum += $0} END{print sum /NR}'
# 10.0847
#$ find . -name "SRR*_withPooledMajorAllele.llr.dbsnp.minCells5_minAvgReads3.5.readsMat.llr_t3.bed" -exec bash -c ' wc -l {}' \; | awk '{sum += $0} END{print sum /NR}'
# 3136.97
#$ find . -name "SRR*_withPooledMajorAllele.llr.dbsnp.minCells5_minAvgReads3.5.readsMat.llr_t4.bed" -exec bash -c 'grep -v "0$" {} | wc -l ' \; | awk '{sum += $0} END{print sum /NR}'
# 4.15254
#$ find . -name "SRR*_withPooledMajorAllele.llr.dbsnp.minCells5_minAvgReads3.5.readsMat.llr_t4.bed" -exec bash -c ' wc -l {}' \; | awk '{sum += $0} END{print sum /NR}'
# 1215.31

