#!/bin/bash

# Fri 06 May 2022 06:30:23 PM PDT
# script to set up simulations using known newick trees and somatic mutations
#for j in 8 9 10 11 ; do
#echo "./setupMutSims.sh muts/params"$j" > muts/params"$j"sim.log"
#done | parallel -j 4 --delay 1

p="$1" # something like segments/params13
i=$(echo "$p" | grep -E "[0-9]+" -o) # pull out the param set number (assumes a number exists in the passed param name)
modType=$(dirname $p)
echo "setting up "$p""

export scripts="/space/s2/sandra/alleleFreqHmm/scripts"
export src="/space/s2/sandra/alleleFreqHmm/src"
export ref="/space/s2/sandra/hg19/hg19_lite.fa"
export windows="/space/s2/sandra/hg19/hg19_lite.250kb_windows"

mkdir -p "$p"
muts/CNV_muts "$modType"/infile"$i".txt "$modType"/paramfile"$i".txt 0 2>&1 > "$p"sim.log
simRet=$?
if [ $simRet -ne 0 ] ; then 
  echo "simulation failed"
  exit -1
fi
mv simu_* "$p"
mv true_* "$p"
#mv "$p".log "$p"
cp "$modType"/infile"$i".txt "$p"/
cp "$modType"/paramfile"$i".txt "$p"/

find $(pwd)/$p -maxdepth 1 -name "simu_readDepth_cancer*[0-9]" | sort -V > $p/tumor_base
find $(pwd)/$p -maxdepth 1 -name "simu_readDepth_healthy*[0-9]" -not -name "*bed" | sort -V > $p/healthy_base

# "convert" to hg19
echo "creating bed files"
while read simu ; do
  cell=$(basename $simu)
  depth=$(pwd)/$p/$cell".hg19_lite.depth"
  paste "$windows" <(awk '{print $4}' $simu) | awk 'BEGIN{OFS="\t"} {if($4 == "") {$4=0} print}' > $depth
done < $p/tumor_base
while read simu ; do
  cell=$(basename $simu)
  depth=$(pwd)/$p/$cell".hg19_lite.depth"
  paste "$windows" <(awk '{print $4}' $simu) | awk 'BEGIN{OFS="\t"} {if($4 == "") {$4=0} print}' > $depth
done < $p/healthy_base
while read simu ; do
  trueFile=$(echo "$simu" | sed 's/simu_readDepth/true_copyNumber/')
  trueBed="$trueFile"".hg19_lite.bed"
  paste "$windows" <(awk '{print $4}' $trueFile) | awk 'BEGIN{OFS="\t"} {if($4 == "") {$4=0} print}' > $trueBed
done < $p/tumor_base
for f in $(find $(pwd)/$p -maxdepth 1 -name "simu_snpAlleles_cancer_cell_*[0-9]" | sort -V) ; do
  $src/convert_simu_snpAlleles_toBed $windows $f > $f".full.hg19_lite.bed"
  # filter out 0,0 mutation lines
  grep -v "0\s0$" $f".full.hg19_lite.bed" > $f".hg19_lite.bed"
  # filter for at least 3 ancestral reads and at least 1 derived read
  awk '{if($4 >= 3 || $5 >= 1) print}' $f".full.hg19_lite.bed" > $f".ancr3_derr1.hg19_lite.bed"
  awk '{if($4 >= 2 || $5 >= 1) print}' $f".full.hg19_lite.bed" > $f".ancr2_derr1.hg19_lite.bed"
  awk '{if($4 >= 1 || $5 >= 1) print}' $f".full.hg19_lite.bed" > $f".ancr1_derr1.hg19_lite.bed"
done

echo "creating sample lists"
find $(pwd)/$p -maxdepth 1 -name "simu_readDepth_cancer*hg19_lite.depth" | sort -V > $p/tumor_depths
shuf $p/tumor_depths > $p/tumor_depths_shuf
head -n 4 $p/tumor_depths_shuf > $p/tumor_depths_4
head -n 20 $p/tumor_depths_shuf > $p/tumor_depths_20
head -n 40 $p/tumor_depths_shuf | tail -n 20 > $p/tumor_depths_40
head -n 60 $p/tumor_depths_shuf | tail -n 20 > $p/tumor_depths_60
head -n 80 $p/tumor_depths_shuf | tail -n 20 > $p/tumor_depths_80
head -n 100 $p/tumor_depths_shuf | tail -n 20 > $p/tumor_depths_100
head -n 120 $p/tumor_depths_shuf | tail -n 20 > $p/tumor_depths_120


for r in $(seq 1 3); do
#for c in $(seq 0 127) ; do
#  mv $p/simu_snpAlleles_cancer_cell_"$c".anc"$r"_der1.hg19_lite.bed $p/simu_snpAlleles_cancer_cell_"$c".ancr"$r"_derr1.hg19_lite.bed
#done
cat $p/tumor_depths     | sed -e 's/readDepth/snpAlleles/' -e 's/.depth/.bed/' -e "s/hg19/ancr"$r"_derr1.hg19/" > $p/tumor_mutList_ancr"$r"_derr1
cat $p/tumor_depths_4   | sed -e 's/readDepth/snpAlleles/' -e 's/.depth/.bed/' -e "s/hg19/ancr"$r"_derr1.hg19/" > $p/tumor_mutList_4_ancr"$r"_derr1
cat $p/tumor_depths_20  | sed -e 's/readDepth/snpAlleles/' -e 's/.depth/.bed/' -e "s/hg19/ancr"$r"_derr1.hg19/" > $p/tumor_mutList_20_ancr"$r"_derr1
cat $p/tumor_depths_40  | sed -e 's/readDepth/snpAlleles/' -e 's/.depth/.bed/' -e "s/hg19/ancr"$r"_derr1.hg19/" > $p/tumor_mutList_40_ancr"$r"_derr1
cat $p/tumor_depths_60  | sed -e 's/readDepth/snpAlleles/' -e 's/.depth/.bed/' -e "s/hg19/ancr"$r"_derr1.hg19/" > $p/tumor_mutList_60_ancr"$r"_derr1
cat $p/tumor_depths_80  | sed -e 's/readDepth/snpAlleles/' -e 's/.depth/.bed/' -e "s/hg19/ancr"$r"_derr1.hg19/" > $p/tumor_mutList_80_ancr"$r"_derr1
cat $p/tumor_depths_100 | sed -e 's/readDepth/snpAlleles/' -e 's/.depth/.bed/' -e "s/hg19/ancr"$r"_derr1.hg19/" > $p/tumor_mutList_100_ancr"$r"_derr1
cat $p/tumor_depths_120 | sed -e 's/readDepth/snpAlleles/' -e 's/.depth/.bed/' -e "s/hg19/ancr"$r"_derr1.hg19/" > $p/tumor_mutList_120_ancr"$r"_derr1

done

#exit

uniqCells=$(find "$p" -name "true_copyNumber_cancer_cell_*hg19_lite.bed" -exec cksum {} + | awk '{ ck[$1$2] = ck[$1$2] ? ck[$1$2] OFS $3 : $3 } END { for (i in ck) print ck[i] }')
echo -e "uniqueness\n$uniqCells" | tee -a "$p"sim.log
numUniqCells=$(echo "$uniqCells" | wc -l)
echo "numUniqCells: "$numUniqCells"" | tee -a "$p"sim.log

avgPloidy=$(for f in $p/true_copyNumber_cancer_cell_*.hg19_lite.bed ; do awk '{sum+=$4} END{print sum/NR}' "$f" ; done | awk '{sum+=$0} END{print sum/NR}')
echo -e "avgPloidy: $avgPloidy" | tee -a "$p"sim.log

echo "ploidy > 10"
numCellsLargePloidy=0
for f in $p/true_copyNumber_cancer_cell_*hg19_lite.bed ; do
  numBinsLargePloidy=$(awk '{if($4 > 10) {print $4}}' $f | sort | uniq -c | awk '{sum+=$1} END{print sum}')
  if [[ ! -z $numBinsLargePloidy ]] ; then
    echo -e "$f""\t""$numBinsLargePloidy" | tee -a "$p"sim.log
    numCellsLargePloidy=$((numCellsLargePloidy + 1))
  fi
done
echo "numCellsLargePloidy: ""$numCellsLargePloidy" | tee -a "$p"sim.log

avgNumDerSites=$(for f in $(find $p -name "simu_snpAlleles_cancer_cell_[0-9]*full.hg19_lite.bed" ) ; do
  grep -v "0$" $f | wc -l
done | awk '{sum+=$1} END{print sum/NR}')
echo -e "avg num derived alleles: $avgNumDerSites" | tee -a "$p"sim.log

Rscript plotTreeSim.R $p

#exit 0

# filter out snps
# create a matrix of tumor ancestral (col 4) and derived (col 5) read counts and coords
for f in $(find $p -name "simu_snpAlleles_cancer_cell*.full.hg19_lite.bed" | sort -V) ; do
  cut -f 4-5 $f > $f".readCounts"
done
readCountList=$(find $p -name "simu_snpAlleles_cancer_cell*.full.hg19_lite.bed.readCounts" | sort -V)
colNames=$(echo "$readCountList" | grep -o -E "cancer_cell_[0-9]+" | awk '{print $0"_ancReads\t"$0"_derReads"}' | sed 's/ /\t/g')
echo -e "chr\tstart\tend\t"$colNames > $p/readCountMatrix.txt
paste <(cut -f 1-3 $p/simu_snpAlleles_cancer_cell_0.full.hg19_lite.bed) $readCountList  >> $p/readCountMatrix.txt

# get coords of sites that appear in enough cells (ie allele must have reads in >= cutoff num cells)
#derPresenceCutoff=3
derPresenceCutoff=0
#ancPresenceCutoff=60
ancPresenceCutoff=0

awk -v cutoff="$derPresenceCutoff" '{count=0; for(i=5; i <= NF; i+=2) {count+=($i!=0)} if(count >= cutoff) {print $0}}' $p/readCountMatrix.txt | cut -f 1-3 |tail -n +2 > $p/mutCoords_"$derPresenceCutoff"cells_der.bed
awk -v cutoff="$ancPresenceCutoff" '{count=0; for(i=4; i <= NF; i+=2) {count+=($i!=0)} if(count >= cutoff) {print $0}}' $p/readCountMatrix.txt | cut -f 1-3 |tail -n +2 > $p/mutCoords_"$ancPresenceCutoff"cells_anc.bed


#downsampleNum="$2"
#downsampleNum=10000
#downsampleNum=7500
#downsampleNum=5000
#downsampleNum=2000
#downsampleNum=1500
#downsampleNum=1000
#downsampleNum=500
for downsampleNum in 500 1000 1500 2000 5000 7500 10000; do
downsampleStr="_downSamp""$downsampleNum"

echo -e "ancPresenceCutoff: ""$ancPresenceCutoff"", derPresenceCutoff: ""$derPresenceCutoff""$downsampleStr"  | tee -a "$p"sim.log
#
#avgNumCellsDer=$(awk '{count=0; for(i=5; i <= NF; i+=2) {count+=($i!=0)} {print count}}' $p/readCountMatrix.txt | tail -n +2 | awk '{if($1 > 0) {sum += $1; nz+=1}} END{print sum/nz}')
#avgNumCellsAnc=$(awk '{count=0; for(i=4; i <= NF; i+=2) {count+=($i!=0)} {print count}}' $p/readCountMatrix.txt | tail -n +2 | awk '{if($1 > 0) {sum += $1; nz+=1}} END{print sum/nz}')
#echo -e "avgNumCellsDer: ""$avgNumCellsDer" | tee -a "$p"sim.log
#echo -e "avgNumCellsAnc: ""$avgNumCellsAnc" | tee -a "$p"sim.log

# make seeding possible for shuf, from https://www.gnu.org/software/coreutils/manual/html_node/Random-sources.html#Random-sources
get_seeded_random()
{
  seed="$1"
  openssl enc -aes-256-ctr -pass pass:"$seed" -nosalt \
    </dev/zero 2>/dev/null
}

# select coords above cutoff, and filter out sites without any reads in this cell
avgNumKeptSites=$(for f in $(find $p -name "simu_snpAlleles_cancer_cell*.full.hg19_lite.bed" | sort -V) ; do
  #out=$(echo $f | sed "s/.full./.anc"$ancPresenceCutoff"_der"$derPresenceCutoff"./")
  #cat $p/mutCoords_"$derPresenceCutoff"cells_der.bed $p/mutCoords_"$ancPresenceCutoff"cells_anc.bed | bedtools intersect -a $f -b - | grep -v "0\s0$" > "$out"
  out=$(echo $f | sed "s/.full./.anc"$ancPresenceCutoff"_der"$derPresenceCutoff""$downsampleStr"./")
  cat $p/mutCoords_"$derPresenceCutoff"cells_der.bed $p/mutCoords_"$ancPresenceCutoff"cells_anc.bed | sort | uniq | bedtools intersect -a $f -b - | grep -v "0\s0$" | shuf -n $downsampleNum --random-source=<(get_seeded_random "$out") | sort -k1,1 -k2,2n > "$out" # downsample to a fixed number of sites with reads
  wc -l $out
done | awk '{sum+=$1} END{print sum/NR}')
echo -e "avgNumKeptSites: ""$avgNumKeptSites" | tee -a "$p"sim.log

avgNumDerSites=$(for f in $(find $p -name "simu_snpAlleles_cancer_cell*.anc"$ancPresenceCutoff"_der"$derPresenceCutoff""$downsampleStr".hg19_lite.bed" | sort -V) ; do
  grep -v "0$" $f | wc -l
done | awk '{sum+=$1} END{print sum/NR}')
echo -e "avgNumDerSites: ""$avgNumDerSites" | tee -a "$p"sim.log


cat $p/tumor_depths     | sed -e 's/readDepth/snpAlleles/' -e 's/.depth/.bed/' > $p/tumor_mutList
cat $p/tumor_depths_4   | sed -e 's/readDepth/snpAlleles/' -e 's/.depth/.bed/' > $p/tumor_mutList_4
cat $p/tumor_depths_20  | sed -e 's/readDepth/snpAlleles/' -e 's/.depth/.bed/' > $p/tumor_mutList_20
cat $p/tumor_depths_40  | sed -e 's/readDepth/snpAlleles/' -e 's/.depth/.bed/' > $p/tumor_mutList_40
cat $p/tumor_depths_60  | sed -e 's/readDepth/snpAlleles/' -e 's/.depth/.bed/' > $p/tumor_mutList_60
cat $p/tumor_depths_80  | sed -e 's/readDepth/snpAlleles/' -e 's/.depth/.bed/' > $p/tumor_mutList_80
cat $p/tumor_depths_100 | sed -e 's/readDepth/snpAlleles/' -e 's/.depth/.bed/' > $p/tumor_mutList_100
cat $p/tumor_depths_120 | sed -e 's/readDepth/snpAlleles/' -e 's/.depth/.bed/' > $p/tumor_mutList_120
#exit
cat $p/tumor_depths     | sed -e 's/readDepth/snpAlleles/' -e 's/.depth/.bed/' -e "s/hg19/anc"$ancPresenceCutoff"_der"$derPresenceCutoff""$downsampleStr".hg19/" > $p/tumor_mutList_anc"$ancPresenceCutoff"_der"$derPresenceCutoff""$downsampleStr"
cat $p/tumor_depths_4   | sed -e 's/readDepth/snpAlleles/' -e 's/.depth/.bed/' -e "s/hg19/anc"$ancPresenceCutoff"_der"$derPresenceCutoff""$downsampleStr".hg19/" > $p/tumor_mutList_4_anc"$ancPresenceCutoff"_der"$derPresenceCutoff""$downsampleStr"
cat $p/tumor_depths_20  | sed -e 's/readDepth/snpAlleles/' -e 's/.depth/.bed/' -e "s/hg19/anc"$ancPresenceCutoff"_der"$derPresenceCutoff""$downsampleStr".hg19/" > $p/tumor_mutList_20_anc"$ancPresenceCutoff"_der"$derPresenceCutoff""$downsampleStr"
cat $p/tumor_depths_40  | sed -e 's/readDepth/snpAlleles/' -e 's/.depth/.bed/' -e "s/hg19/anc"$ancPresenceCutoff"_der"$derPresenceCutoff""$downsampleStr".hg19/" > $p/tumor_mutList_40_anc"$ancPresenceCutoff"_der"$derPresenceCutoff""$downsampleStr"
cat $p/tumor_depths_60  | sed -e 's/readDepth/snpAlleles/' -e 's/.depth/.bed/' -e "s/hg19/anc"$ancPresenceCutoff"_der"$derPresenceCutoff""$downsampleStr".hg19/" > $p/tumor_mutList_60_anc"$ancPresenceCutoff"_der"$derPresenceCutoff""$downsampleStr"
cat $p/tumor_depths_80  | sed -e 's/readDepth/snpAlleles/' -e 's/.depth/.bed/' -e "s/hg19/anc"$ancPresenceCutoff"_der"$derPresenceCutoff""$downsampleStr".hg19/" > $p/tumor_mutList_80_anc"$ancPresenceCutoff"_der"$derPresenceCutoff""$downsampleStr"
cat $p/tumor_depths_100 | sed -e 's/readDepth/snpAlleles/' -e 's/.depth/.bed/' -e "s/hg19/anc"$ancPresenceCutoff"_der"$derPresenceCutoff""$downsampleStr".hg19/" > $p/tumor_mutList_100_anc"$ancPresenceCutoff"_der"$derPresenceCutoff""$downsampleStr"
cat $p/tumor_depths_120 | sed -e 's/readDepth/snpAlleles/' -e 's/.depth/.bed/' -e "s/hg19/anc"$ancPresenceCutoff"_der"$derPresenceCutoff""$downsampleStr".hg19/" > $p/tumor_mutList_120_anc"$ancPresenceCutoff"_der"$derPresenceCutoff""$downsampleStr"

echo -n > $p/mutFilt_anc"$ancPresenceCutoff"_der"$derPresenceCutoff""$downsampleStr".log
while read ln ; do
  summary=$(awk '{if($5 != 0) {count +=1}} END{print count"\t"NR"\t"count/NR}' "$ln")
  echo -e "$ln""\t""$summary" >> $p/mutFilt_anc"$ancPresenceCutoff"_der"$derPresenceCutoff""$downsampleStr".log
done < $p/tumor_mutList_anc"$ancPresenceCutoff"_der"$derPresenceCutoff""$downsampleStr"

done

#exit 0

find $(pwd)/$p -maxdepth 1 -name "simu_readDepth_healthy*hg19_lite.depth" -not -name "*bed" | sort -V > $p/healthy_depths
find $(pwd)/$p -maxdepth 1 -name "simu_readDepth_cancer*hg19_lite.depth" | sort -V | sed 's/.depth//g' > $p/tumor_base_hg19
cat $p/tumor_base_hg19 | parallel -j 1 basename {} > $p/tumor_ids_hg19

# create snpSummaryCN_* files
numTumorCells=$(awk '{print $3; exit}' "$modType"/infile"$i".txt)
for j in $(seq 0 $(($numTumorCells -1))) ; do
  #paste <(cut -f 2,4,5 $p/simu_snpAlleles_cancer_cell_"$j" ) $p/true_snpAlleles_cancer_cell_"$j" | awk '{print $4"\t"$1"\t"$5"\t"$6"\t"$2"\t"$3}' > $p/snpSummary_cancer_cell_"$j"
  paste <(cut -f 2,4,5 $p/simu_snpAlleles_cancer_cell_"$j" ) $p/true_snpAlleles_cancer_cell_"$j" | sed 's/^\t/0\t0\t0\t/' | awk '{print $4"\t"$1"\t"$5"\t"$6"\t"$2"\t"$3}' > $p/snpSummary_cancer_cell_"$j"
  # snpSummary_cancer ==> snpIdx(0based) binNum(1based) ancAll derAll ancRead derRead
  awk '{print NR"\t"$4}' $p/true_copyNumber_cancer_cell_"$j" > $p/true_copyNumber_binNum_cancer_cell_"$j"
  # true_copyNumber_cancer_cell_binNum_ ==> binNum(1based) CN
  #join <(sort -k2,2b $p/snpSummary_cancer_cell_"$j") <(sort -k1,1b $p/true_copyNumber_binNum_cancer_cell_"$j") -1 2 -2 1 -t $'\t' | sort -k1,1n > $p/snpSummaryCN_cancer_cell_"$j"
  #export LANG=en_EN
  LC_ALL=C join <(LC_ALL=C sort -k1 <(awk '{print $2"\t"$1"\t"$3"\t"$4"\t"$5"\t"$6}' $p/snpSummary_cancer_cell_"$j") ) <(LC_ALL=C sort -k1 $p/true_copyNumber_binNum_cancer_cell_"$j") -1 1 -2 1 -t $'\t' | sort -k1,1n > $p/snpSummaryCN_cancer_cell_"$j" 2> $p/snpSummaryCN_cancer_cell_"$j".err
  #LC_ALL=C join <(LC_ALL=C sort -k1 $p/snpSummary_cancer_cell_"$j") <(LC_ALL=C sort -k1 $p/true_copyNumber_binNum_cancer_cell_"$j") -1 2 -2 1 -t $'\t' | sort -k1,1n > $p/snpSummaryCN_cancer_cell_"$j" 2> $p/snpSummaryCN_cancer_cell_"$j".err
  if [[ $(cat $p/snpSummaryCN_cancer_cell_"$j".err | wc -l) != 0 ]] ; then
    echo "join <(awk '{print $2"\t"$1"\t"$3"\t"$4"\t"$5"\t"$6}' $p/snpSummary_cancer_cell_"$j" | sort) <(sort $p/true_copyNumber_binNum_cancer_cell_"$j") -1 1 -2 1 -t $'\t' | sort -k1,1n > $p/snpSummaryCN_cancer_cell_"$j""
    exit
  fi
  # snpSummaryCN_ ==> binNum(1based) snpIdx(0based) ancAll derAll ancRead derRead CN
done 

exit
echo "Rscript "$scripts"/avgDiploid.R $p/healthy_depths $p/simu_healthy_avg.bed"
Rscript "$scripts"/avgDiploid.R $p/healthy_depths $p/simu_healthy_avg.bed
if [ ! -f $p/params"$i".meanVar ] ; then
  echo "Rscript "$scripts"/fitMeanVarRlnshp.R $p/healthy_depths $p/params"$i".meanVar # ~ 30 mins"
  Rscript "$scripts"/fitMeanVarRlnshp.R $p/healthy_depths $p/params"$i".meanVar # ~ 30 mins
fi



