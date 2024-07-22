#!/bin/bash

# Check if input VCF file and output directory are provided
if [ $# -ne 2 ]; then
    echo "Usage: $0 <input_vcf_file> <output_directory>"
    exit 1
fi

# Input VCF file
input_vcf="$1"

# Output directory
output_dir="$2"
mkdir -p "$output_dir"

# Subset the data set using 1% of the variants
#bcftools view "$input_vcf" | vcfrandomsample -r 0.01 > "${output_dir}/subset.vcf"

# Compress subset VCF
#gzip "${output_dir}/subset.vcf"

# Calculate allele frequency
#vcftools --gzvcf "${output_dir}/subset.vcf.gz" --freq2 --out "${output_dir}/allele_freq2" --min-alleles 2 --max-alleles 2

# Mean depth per individual
#vcftools --gzvcf "${output_dir}/subset.vcf.gz" --depth --out "${output_dir}/mean_depth_per_individual"

# Mean depth per site
#vcftools --gzvcf "${output_dir}/subset.vcf.gz" --site-mean-depth --out "${output_dir}/mean_depth_per_site"

# Calculating site quality
#vcftools --gzvcf "${output_dir}/subset.vcf.gz" --site-quality --out "${output_dir}/site_quality"

# Proportion of Missing Data per indv
#vcftools --gzvcf "${output_dir}/subset.vcf.gz" --missing-indv --out "${output_dir}/missing_data_per_indv"

# Missing data per site
#vcftools --gzvcf "${output_dir}/subset.vcf.gz" --missing-site --out "${output_dir}/missing_data_per_site"

# Calculating heterozygosity and inbreeding coefficient per individual
vcftools --gzvcf "${output_dir}/subset.vcf.gz" --het --out "${output_dir}/heterozygosity_inbreeding"

