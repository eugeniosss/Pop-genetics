#!/bin/bash

set -e 

# Check if correct number of arguments are provided
if [ "$#" -ne 4 ]; then
    echo "Usage: $0 <bed_file> <vcf_file> <tmp_file> <final_file>"
    exit 1
fi

# Assign arguments to variables
BED_FILE="$1"
VCF_FILE="$2"
TMP_FILE="$3"
FINAL="$4"

echo NOTE!! ${BED_FILE} NEEDS TO BE SORTED!!

#NUMBER_OF_INDS=$(awk '{if ($1 == "#CHROM"){print NF-9; exit}}' $VCF_FILE)

echo "Sample $(grep -m 1 "CHR" $VCF_FILE | cut -f 10- | tr '\t' ' ')" >> $TMP_FILE

paste -d " " <(awk '{print $4}' ${BED_FILE}) <(bedtools intersect -a $VCF_FILE -b ${BED_FILE} -wa | awk '{for(i=10; i<=NF; i++) {split($i, a, ":"); printf "%s ", a[1]} printf "\n"}') >> $TMP_FILE

awk '{ for (i=1; i<=NF; i++) RtoC[i]= (i in RtoC?RtoC[i] OFS :"") $i; } 
        END{ for (i=1; i<=NF; i++) print RtoC[i] }' $TMP_FILE > $FINAL

rm $TMP_FILE
