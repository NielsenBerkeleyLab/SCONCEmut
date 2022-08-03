#!/bin/bash

# Thu 05 May 2022 03:43:33 PM PDT
# script to
# - add pooled major allele to indv diploid cell read counts
# - calc LLR for each site (relative to pooled major allele)
# - add "inDbsnp" "flag"

# output:
# chr start end pooledMajorAllele cellMajorAllele totalNumReads numPooledMajorAlleleReads LLR dbsnpChr dbsnpStart dbsnpEnd
# dbsnpChr, dbsnpStart, dbsnpEnd are ".	-1	-1" if not in dbsnp

# based on dbsnpIntersection/README and pipeline_overview.txt code

cellID="$1"
pileupDir="$2"
dbsnpOutputDir="$3"
mkdir -p "$dbsnpOutputDir"

# first create all the $pileupDir/*_minor_allele_reads.bed files
# 'chr', 'start', 'end', 'total_read_count', 'cell_major_allele', 'minor_allele_read_count'
#find "$pileupDir" -name *.mpileup | parallel -j 2 python3 calculate_cell_read_counts.py -f {} '>' {}".err" '2>&1'
cellPileup=$(find "$pileupDir" -maxdepth 1 -name "$cellID""*mpileup")
python3 calculate_cell_read_counts.py -f "$cellPileup" -o "$pileupDir" > "$cellPileup"".err" 2>&1


#for cellID in $(find "$pileupDir" -name "*mpileup" -exec basename {} .q20.adaptTr.qualTr.dusted.sorted.noDups.mpileup \; | sort) ; do
  #echo $cellID

  # 1:chr 2:start 3:end 4:total_num_reads      5:cell_major_allele 6:minor_allele_read_counts
  # 7:chr 8:start 9:end 10:pooled_major_allele
  sort -k1,1 -k2,2n "$pileupDir""/""$cellID""_minor_allele_reads.bed" | bedtools intersect -sorted -wa -wb -a stdin -b "pooledMajorAlleleBedFiles/diploid_sorted.bed" | awk '{print $1"\t"$2"\t"$3"\t"$10"\t"$5"\t"$4"\t"$6}' | tee "$pileupDir""/""$cellID""_withPooledMajorAllele.bed" | ./calcLLR > "$pileupDir""/""$cellID""_withPooledMajorAllele.llr.bed" 2> "$pileupDir""/""$cellID""_withPooledMajorAllele.llr.err"
  # chr start end pooledMajorAllele cellMajorAllele totalNumReads numPooledMajorAlleleReads LLR 
  #cat "$pileupDir""/""$cellID""_withPooledMajorAllele.bed" | ./calcLLR > "$pileupDir""/""$cellID""_withPooledMajorAllele.llr.bed" 2> "$pileupDir""/""$cellID""_withPooledMajorAllele.llr.err"
  
  bedtools intersect -sorted -loj -a "$pileupDir""/""$cellID""_withPooledMajorAllele.llr.bed" -b "dbsnpFiles/dbsnp_sorted.bed" > "$dbsnpOutputDir""/""$cellID""_withPooledMajorAllele.llr.dbsnp.bed"
#done # end cell loop

