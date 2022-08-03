#!/bin/bash

# Fri 13 May 2022 02:46:30 PM PDT
# script to create a big matrix of tumor LLR values, ***relative to pooled diploid major allele***. That is, every "pooled major allele" referenced is from the pooled diploid
# based on createPooledDiploidLLRmatrix.sh, just different file names

indvBedList=$(find dbsnpIntersection_tumor -maxdepth 1 -name "*_withPooledMajorAllele.llr.dbsnp.bed" | sort -V)

# first get a `chr start end` coord file of all coords from all tumor cells
bedops -m $indvBedList | bedops --chop 1 - > dbsnpIntersection_tumor/pooledTumorCoords.bed

# filter for positions that have diploid data
#bedtools intersect -sorted <(sort -k1,1 -k2,2n dbsnpIntersection/pooledDiploidCoords.bed) <(sort -k1,1 -k2,2n dbsnpIntersection_tumor/pooledTumorCoords.bed) > dbsnpIntersection_tumor/pooledDiploidTumorCoords.bed
bedtools intersect -sorted -a dbsnpIntersection/pooledDiploidCoords.bed -b dbsnpIntersection_tumor/pooledTumorCoords.bed > dbsnpIntersection_tumor/pooledDiploidTumorCoords.bed

# and add dbsnp status to it
bedtools intersect -sorted -loj -a dbsnpIntersection_tumor/pooledDiploidTumorCoords.bed -b dbsnpFiles/dbsnp_sorted.bed | cut -f 1-4 > dbsnpIntersection_tumor/pooledDiploidTumorCoords.dbsnp.bed

# for each cell, fill in NAs for lines in pooled coord file that aren't in cell file
for f in $indvBedList ; do
  out=$(echo "$f" | sed 's/.bed/.diploidTumorCoords.NAs.bed/')
  # chr start end pooledMajorAllele totalNumReads numPooledMajorAlleleReads LLR dbsnpChr
  bedops -n 1 dbsnpIntersection_tumor/pooledDiploidTumorCoords.bed "$f" | awk '{print $0"\tNA\tNA\tNA\tNA\tNA"}' | bedops -u - <(cut -f1-4,6-9 "$f") | cut -f5-7> "$out"
done

# combine into a big matrix
naBedList=$(find dbsnpIntersection_tumor -maxdepth 1 -name "*_withPooledMajorAllele.llr.dbsnp.diploidTumorCoords.NAs.bed" | sort -V )
colNames=$(echo "$naBedList" | xargs -I '{}' basename '{}' | sed -e 's/_withPooledMajorAllele.llr.dbsnp.diploidTumorCoords.NAs.bed//' | awk '{print $0"_totalNumReads\t"$0"_numPooledMajorAlleleReads\t"$0"_LLR"}' | tr '\n' '\t')

echo -e "chr\tstart\tend\tdbsnp\t""$colNames" > dbsnpIntersection_tumor/pooledDiploidTumorMajorAllele.llr.dbsnp.bed
paste dbsnpIntersection_tumor/pooledDiploidTumorCoords.dbsnp.bed $naBedList >> dbsnpIntersection_tumor/pooledDiploidTumorMajorAllele.llr.dbsnp.bed

