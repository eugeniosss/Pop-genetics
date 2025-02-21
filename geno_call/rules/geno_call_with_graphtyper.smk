rule geno_call_with_graphtyper:
    input:
        reference=config["ref_genome"],  # Reference genome
        bed_file="split_bed/{chr}.bed"   # Per-chromosome BED file
    output:
        "graphtyper/{chr}_variants.vcf"   # Output VCF per chromosome
    log:
        "graphtyper/{chr}.log"
    conda:
        config["dir"] + "envs/geno_callers_env.yml"
    threads: 1
    params:
        bams=config["bams"],  # BAM file list (define in config or Snakefile)
    shell:
        """
        mkdir -p graphtyper
        
        awk -F "\t" '{{print $1 ":" $2 "-" $3}}' {input.bed_file} > split_bed/{wildcards.chr}_graphtyper.bed


        graphtyper genotype {input.reference} --sams={params.bams} --region_file=split_bed/{wildcards.chr}_graphtyper.bed --threads=1 --output={wildcards.chr} > {log} 2>&1
        
        bcftools concat \
            -o {output} \
            {wildcards.chr}/*/*.vcf.gz
        
        rm -r {wildcards.chr} split_bed/{wildcards.chr}_graphtyper.bed
        """   