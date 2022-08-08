#!/bin/bash

# Mon 14 Feb 2022 09:16:12 PM PST
# wrapper to call cnaBedToFasta.sh on all paramsets

cnp2cnpDir="."
k=10
for p in "muts/params14" "muts/params15" "muts/params16" "muts/params17" "muts/params32" "muts/params33" "muts/params34" ; do
  dataset=$(echo "$p" | sed 's/\//_/')
  for c in "20" ; do
  #for c in "20" "40" "60" "80" "100" "120" ; do
    # first do true bed files if they don't already exist
    trueBase="$p""/tumor_depths_""$c""_true"
    if [[ ! -f "$trueBase""Rev.zzs.cnp2cnp" ]] ; then
      echo -n > "$trueBase"".fasta"
      while read simu ; do
        bed=$(echo $simu | sed -e 's/simu/true/' -e 's/depth/bed/' -e 's/readDepth/copyNumber/')
        line1=">""$bed"
        line2=$(awk '{printf("%.0f\n", $4)}' "$bed" | tr '\n' ','i | sed 's/,$//')
        echo -e "$line1""\n""$line2" >> "$trueBase"".fasta"
      done < "$p""/tumor_depths_""$c"

      echo -n > "$trueBase""Rev.fasta"
      tac "$p""/tumor_depths_""$c" | while read simu ; do
        bed=$(echo $simu | sed -e 's/simu/true/' -e 's/depth/bed/' -e 's/readDepth/copyNumber/')
        line1=">""$bed"
        line2=$(awk '{printf("%.0f\n", $4)}' "$bed" | tr '\n' ','i | sed 's/,$//')
        echo -e "$line1""\n""$line2" >> "$trueBase""Rev.fasta"
      done

      python "$cnp2cnpDir"/cnp2cnp.py -m "matrix" -d any -i "$trueBase"".fasta" -o "$trueBase".cnp2cnp
      python "$cnp2cnpDir"/cnp2cnp.py -m "matrix" -d zzs -i "$trueBase"".fasta" -o "$trueBase".zzs.cnp2cnp
      python "$cnp2cnpDir"/cnp2cnp.py -m "matrix" -d any -i "$trueBase""Rev.fasta" -o "$trueBase"Rev.cnp2cnp
      python "$cnp2cnpDir"/cnp2cnp.py -m "matrix" -d zzs -i "$trueBase""Rev.fasta" -o "$trueBase"Rev.zzs.cnp2cnp
    fi

    for scAllPFileKey in "scAllPMuts_v16" "scAllP_v26" ; do
      #for bedType in "sconce" "mean" "median" "mode" ; do
      for bedType in "sconce" "mean" ; do
        outBase="$p""/output_""$scAllPFileKey""_""$dataset""_k""$k""_c""$c""_rounded""$bedType"

        # create fastas
        scripts/cnaBedToFasta.sh "$p" "$k" "$c" "$bedType" "$scAllPFileKey" "F" > "$outBase".fasta
        scripts/cnaBedToFasta.sh "$p" "$k" "$c" "$bedType" "$scAllPFileKey" "T" > "$outBase"Rev.fasta

        # forward and reversed. calc both since cnp2cnp is not symmetric
        python "$cnp2cnpDir"/cnp2cnp/cnp2cnp.py -m "matrix" -d any -i "$outBase"".fasta" -o "$outBase".cnp2cnp
        python "$cnp2cnpDir"/cnp2cnp/cnp2cnp.py -m "matrix" -d zzs -i "$outBase"".fasta" -o "$outBase".zzs.cnp2cnp
        python "$cnp2cnpDir"/cnp2cnp/cnp2cnp.py -m "matrix" -d any -i "$outBase""Rev.fasta" -o "$outBase"Rev.cnp2cnp
        python "$cnp2cnpDir"/cnp2cnp/cnp2cnp.py -m "matrix" -d zzs -i "$outBase""Rev.fasta" -o "$outBase"Rev.zzs.cnp2cnp

      done
    done
  done
done

