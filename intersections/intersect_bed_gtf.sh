#!/bin/bash

set -e

# Default values
BED=""
GTF=""
OUTPUT="overlaps.bed"
OPTION=""

# Function to display script usage
usage() {
    echo "Usage: $0 -b <BED_FILE> -g <GTF_FILE> [-s <OPTION>] [-o <OUTPUT_FILE>]"
    echo "Options:"
    echo "  -s <OPTION>: Specify 'upper' or 'lower' to preserve only one entry per gene based on maximum or minimum value in the 4th column, respectively."
    exit 1
}

# Parse command-line options
while getopts ":b:g:s:o:" opt; do
    case $opt in
        b)
            BED="$OPTARG"
            ;;
        g)
            GTF="$OPTARG"
            ;;
        s)
            if [[ "$OPTARG" == "upper" || "$OPTARG" == "lower" ]]; then
                OPTION="$OPTARG"
            else
                echo "Invalid option for -s: $OPTARG. Must be 'upper' or 'lower'."
                usage
            fi
            ;;
        o)
            OUTPUT="$OPTARG"
            ;;
        \?)
            echo "Invalid option: -$OPTARG" >&2
            usage
            ;;
        :)
            echo "Option -$OPTARG requires an argument." >&2
            usage
            ;;
    esac
done

# Check if required options are provided
if [ -z "$BED" ] || [ -z "$GTF" ]; then
    echo "BED and GTF files are required options."
    usage
fi


# Intersect BED file with filtered GTF file and extract relevant information
bedtools intersect -a $BED -b <(grep "gene_name" $GTF) -wa -wb | awk -v OFS='\t' '{if($3-$2 > 0) print $1, $2, $3, $18}' | sed -e 's,;,,g; s,",,g'  > ${OUTPUT}_tmp

awk 'NR==FNR{a[$1"_"$2"_"$3]=$4;next}{if (!($1"_"$2"_"$3 in b)) {print $0,a[$1"_"$2"_"$3]; b[$1"_"$2"_"$3]}}' $BED  ${OUTPUT}_tmp > ${OUTPUT}_tmp2

# Preserve only one entry per gene based on maximum or minimum value in the 4th column
if [ "$OPTION" == "upper" ]; then
    awk '!seen[$4] || $5 > max[$4] {max[$4] = $5; line[$4] = $0} {seen[$4] = 1} END {for (key in line) print line[key]}' "${OUTPUT}_tmp2" > "$OUTPUT"
elif [ "$OPTION" == "lower" ]; then
    awk '!seen[$4] || $5 < min[$4] {min[$4] = $5; line[$4] = $0} {seen[$4] = 1} END {for (key in line) print line[key]}' "${OUTPUT}_tmp2" > "$OUTPUT"
else
    cp "${OUTPUT}_tmp2" "$OUTPUT"
fi

# Clean up temporary files
rm -f "${OUTPUT}_tmp" "${OUTPUT}_tmp2"

echo "Overlaps saved to $OUTPUT"
