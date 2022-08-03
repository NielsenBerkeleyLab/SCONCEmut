#!/bin/bash

# script to combine individual diploid llr bed files into a big matrix for plotting

# input is indv bed files from ./createDiploidLLRdbsnpBedFiles.sh
# chr start end pooledMajorAllele cellMajorAllele totalNumReads numPooledMajorAlleleReads LLR dbsnpChr dbsnpStart dbsnpEnd


# following https://www.biostars.org/p/320196/#320204

indvBedList=$(find dbsnpIntersection -maxdepth 1 -name "*_withPooledMajorAllele.llr.dbsnp.bed" | sort -V)

# first get a `chr start end` coord file of all coords from all diploid cells
bedops -m $indvBedList | bedops --chop 1 - > dbsnpIntersection/pooledDiploidCoords.bed
# and add dbsnp status to it
bedtools intersect -sorted -loj -a dbsnpIntersection/pooledDiploidCoords.bed -b dbsnpFiles/dbsnp_sorted.bed | cut -f 1-4 > dbsnpIntersection/pooledDiploidCoords.dbsnp.bed

# for each cell, fill in NAs for lines in pooled coord file that aren't in cell file
for f in $indvBedList ; do
  out=$(echo "$f" | sed 's/.bed/.NAs.bed/')
  # chr start end pooledMajorAllele totalNumReads numPooledMajorAlleleReads LLR dbsnpChr
  bedops -n 1 dbsnpIntersection/pooledDiploidCoords.bed "$f" | awk '{print $0"\tNA\tNA\tNA\tNA\tNA"}' | bedops -u - <(cut -f1-4,6-9 "$f") | cut -f5-7> "$out"
done

# combine into a big matrix
naBedList=$(find dbsnpIntersection -maxdepth 1 -name "*_withPooledMajorAllele.llr.dbsnp.NAs.bed" | sort -V )
colNames=$(echo "$naBedList" | xargs -I '{}' basename '{}' | sed -e 's/_withPooledMajorAllele.llr.dbsnp.NAs.bed//' | awk '{print $0"_totalNumReads\t"$0"_numPooledMajorAlleleReads\t"$0"_LLR"}' | tr '\n' '\t')

echo -e "chr\tstart\tend\tdbsnp\t""$colNames" > dbsnpIntersection/pooledDiploidMajorAllele.llr.dbsnp.bed
paste dbsnpIntersection/pooledDiploidCoords.dbsnp.bed $naBedList >> dbsnpIntersection/pooledDiploidMajorAllele.llr.dbsnp.bed

