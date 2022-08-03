#!/bin/bash

# Mon 14 Feb 2022 09:16:12 PM PST
# wrapper to call cnaBedToFasta.sh on all paramsets

#p="$1" # "segments/params21"
#k="$2" # "10"
#c="$3" # "20"
#bedType="$4" # "mean", "mode", "median". Note cnp2cnp expects integers (which mean cannot guarantee)

#scAllPFileKey="scAllP_v21"
k=10
#for p in "segments/params21" "segments/params22" "segments/params23" "segments/params24" "segments/params25" "segments/params26" "segments/params27" "segments/params28" "segments/params29" "segments/params30" ; do
#for p in "segments/params21" "segments/params22" "segments/params23" "segments/params24" "segments/params25" "segments/params26" "segments/params27" "segments/params28" ; do
#for p in "segments/params31" ; do
#for p in "segments/params21" "segments/params22" "segments/params23" "segments/params24" "segments/params25" "segments/params26" ; do
#for p in "segments/params21" "segments/params22" "segments/params25" "segments/params26" ; do
#for p in "muts/params14" "muts/params15" "muts/params16" "muts/params17" ; do #"muts/params20" "muts/params32" "muts/params33" "muts/params34" ; do
#for p in "muts/params33" "muts/params34" ; do
for p in "muts/params32" "muts/params33" "muts/params34" ; do
#for p in "segments/params22" "segments/params24" ; do
  dataset=$(echo "$p" | sed 's/\//_/')
  for c in "20" ; do
  #for c in "40" "60" "80" "100" "120" ; do
  #for c in "20" "40" "60" "80" "100" "120" "128"; do
  #for c in "20" "40" "60" "80" "100" "120" ; do
    # first do true bed files if they don't already exist
    #trueBase="$p""/tumor_depths_""$c""_true"
    #if [[ ! -f "$trueBase""Rev.zzs.cnp2cnp" ]] ; then
    ##if [[ ! -f "$trueBase""Rev.cnp2cnp" ]] ; then
    #  #echo -n > "$trueBase"".fasta"
    #  #while read simu ; do
    #  #  bed=$(echo $simu | sed -e 's/simu/true/' -e 's/depth/bed/' -e 's/readDepth/copyNumber/')
    #  #  line1=">""$bed"
    #  #  line2=$(awk '{printf("%.0f\n", $4)}' "$bed" | tr '\n' ','i | sed 's/,$//')
    #  #  echo -e "$line1""\n""$line2" >> "$trueBase"".fasta"
    #  #done < "$p""/tumor_depths_""$c"

    #  #echo -n > "$trueBase""Rev.fasta"
    #  #tac "$p""/tumor_depths_""$c" | while read simu ; do
    #  #  bed=$(echo $simu | sed -e 's/simu/true/' -e 's/depth/bed/' -e 's/readDepth/copyNumber/')
    #  #  line1=">""$bed"
    #  #  line2=$(awk '{printf("%.0f\n", $4)}' "$bed" | tr '\n' ','i | sed 's/,$//')
    #  #  echo -e "$line1""\n""$line2" >> "$trueBase""Rev.fasta"
    #  #done 

    #  #python /space/s2/sandra/methodComparisons/cnp2cnp/cnp2cnp.py -m "matrix" -d euclidean -i "$trueBase"".fasta" -o "$trueBase"".euc.cnp2cnp"
    #  #python /space/s2/sandra/methodComparisons/cnp2cnp/cnp2cnp.py -m "matrix" -d euclidean -i "$trueBase"".fasta" -o "$trueBase""Rev.euc.cnp2cnp"
    #  echo "~/Python-3.10.5/debug/python /space/s1/sandra/cnp2cnp/cnp2cnp.py -m "matrix" -d any -i "$trueBase"".fasta" -o "$trueBase"".cnp2cnp""
    #  echo "~/Python-3.10.5/debug/python /space/s1/sandra/cnp2cnp/cnp2cnp.py -m "matrix" -d zzs -i "$trueBase"".fasta" -o "$trueBase"".zzs.cnp2cnp""
    #  echo "~/Python-3.10.5/debug/python /space/s1/sandra/cnp2cnp/cnp2cnp.py -m "matrix" -d any -i "$trueBase""Rev.fasta" -o "$trueBase""Rev.cnp2cnp""
    #  echo "~/Python-3.10.5/debug/python /space/s1/sandra/cnp2cnp/cnp2cnp.py -m "matrix" -d zzs -i "$trueBase""Rev.fasta" -o "$trueBase""Rev.zzs.cnp2cnp""
    #fi

   ##for key in "scAllP_v21" ; do
    #for scAllPFileKey in "scAllP_v24" ; do
    #for scAllPFileKey in "scAllP_v26" "scAllPMuts_v14" ; do
    #for scAllPFileKey in "scAllP_v26" ; do #"scAllPMuts_v14" ; do
    for scAllPFileKey in "scAllPMuts_v16" ; do
    #for scAllPFileKey in "scAllPMuts_v16_reuseMutEsts_shortcut" ; do
    #for scAllPFileKey in "scAllP_v24_nearest10" ; do
      # check if last cell in this tumor_depths_$c set has a mode file (ie is this run done)
      lastCell=$(basename "$(tail -n 1 "$p""/tumor_depths_""$c")")
      lastCellMode="$p""/output_""$scAllPFileKey""_""$dataset""_k""$k""_c""$c""__""$lastCell""__k""$k""__mode.bed"
      #if [[ ! -f "$lastCellMode" ]] ; then
      #  echo "could not find last cell's mode bed file: $lastCellMode. skipping"
      #  continue
      #fi
       #for bedType in "sconce" "mean" "median" "mode" ; do
       for bedType in "sconce" "mean" ; do
      #for bedType in "sconce" ; do
        outBase="$p""/output_""$scAllPFileKey""_""$dataset""_k""$k""_c""$c""_rounded""$bedType"
        #echo "going to run ""$outBase"

        #if [[ -f "$outBase""Rev.zzs.cnp2cnp" ]] ; then
        #  continue
        #fi

        ## create fastas
        #echo "scripts/cnaBedToFasta.sh "$p" "$k" "$c" "$bedType" "$scAllPFileKey" "F" > "$outBase"".fasta""
        #echo "scripts/cnaBedToFasta.sh "$p" "$k" "$c" "$bedType" "$scAllPFileKey" "T" > "$outBase""Rev.fasta""

        # forward
        #python /space/s2/sandra/methodComparisons/cnp2cnp/cnp2cnp.py -m "matrix" -d euclidean -i "$outBase"".fasta" -o "$outBase"".euc.cnp2cnp"
        #python /space/s2/sandra/methodComparisons/cnp2cnp/cnp2cnp.py -m "matrix" -d euclidean -i "$outBase""Rev.fasta" -o "$outBase""Rev.euc.cnp2cnp"
        # reversed. calc both since cnp2cnp is not symmetric
        echo "~/Python-3.10.5/debug/python /space/s1/sandra/cnp2cnp/cnp2cnp.py -m "matrix" -d any -i "$outBase"".fasta" -o "$outBase"".cnp2cnp""
        echo "~/Python-3.10.5/debug/python /space/s1/sandra/cnp2cnp/cnp2cnp.py -m "matrix" -d zzs -i "$outBase"".fasta" -o "$outBase"".zzs.cnp2cnp""
        echo "~/Python-3.10.5/debug/python /space/s1/sandra/cnp2cnp/cnp2cnp.py -m "matrix" -d any -i "$outBase""Rev.fasta" -o "$outBase""Rev.cnp2cnp""
        echo "~/Python-3.10.5/debug/python /space/s1/sandra/cnp2cnp/cnp2cnp.py -m "matrix" -d zzs -i "$outBase""Rev.fasta" -o "$outBase""Rev.zzs.cnp2cnp""

        ##wait

      done
    done
  done
done

#for p in "segments/params27" "segments/params28" ; do
#  dataset=$(echo "$p" | sed 's/\//_/')
#  for c in "extremes" ; do
#    for key in "scAllP_v21" ; do
#      for bedType in "mean" "median" "mode" ; do
#        outBase="$p""/output_""$scAllPFileKey""_""$dataset""_k""$k""_c""$c""_rounded""$bedType"
#        nearestOutBase="$p""/output_""$scAllPFileKey""_nearest10_""$dataset""_k""$k""_c""$c""_rounded""$bedType"
#
#        # forward
#        scripts/cnaBedToFasta.sh "$p" "$k" "$c" "$bedType" "$scAllPFileKey" "F" > "$outBase"".fasta"
#        retVal=$?
#        if [ $retVal -ne 0 ]; then
#          continue
#        fi
#        python /space/s2/sandra/methodComparisons/cnp2cnp/cnp2cnp.py -m "matrix" -i "$outBase"".fasta" -o "$outBase"".cnp2cnp"
#
#        scripts/cnaBedToFasta.sh "$p" "$k" "$c" "$bedType" "$scAllPFileKey" "F" "nearest10_" > "$nearestOutBase"".fasta"
#        python /space/s2/sandra/methodComparisons/cnp2cnp/cnp2cnp.py -m "matrix" -i "$nearestOutBase"".fasta" -o "$nearestOutBase"".cnp2cnp"
#
#        # reversed. calc both since cnp2cnp is not symmetric
#        scripts/cnaBedToFasta.sh "$p" "$k" "$c" "$bedType" "$scAllPFileKey" "T" > "$outBase""Rev.fasta"
#        python /space/s2/sandra/methodComparisons/cnp2cnp/cnp2cnp.py -m "matrix" -i "$outBase""Rev.fasta" -o "$outBase""Rev.cnp2cnp"
#
#        scripts/cnaBedToFasta.sh "$p" "$k" "$c" "$bedType" "$scAllPFileKey" "T" "nearest10_" > "$nearestOutBase""Rev.fasta"
#        python /space/s2/sandra/methodComparisons/cnp2cnp/cnp2cnp.py -m "matrix" -i "$nearestOutBase""Rev.fasta" -o "$nearestOutBase""Rev.cnp2cnp"
#
#      done
#    done
#  done
#done
# 

