rule prepare_bed_file_for_graphtyper:
    input:
        "split_bed/{chr}.bed"   # Per-chromosome BED file
    output:
        temp("split_bed_graph/{chr}.bed")
    log:
        "split_bed_graph/{chr}.log"
    shell:
        """
        (awk -F "\t" '{{print $1 ":" $2 "-" $3}}' {input} > {output}) > {log} 2>&1
        """

rule geno_call_with_graphtyper:
    input:
        reference=config["ref_genome"],  # Reference genome
        bed="split_bed_graph/{chr}.bed"
    output:
        temp("graphtyper_{dataset}/{chr}_variants.vcf")   # Output VCF per chromosome
    log:
        "graphtyper_{dataset}/{chr}.log"
    conda:
        config["dir"] + "envs/geno_callers_env.yml"
    threads: 1
    params:
        bams=lambda wildcards: f"{wildcards.dataset}.txt",  # Ensure dataset is a wildcard
        options=config["options"]["graphtyper"]
    shell:
        """    
        (graphtyper genotype \
            {input.reference} \
            --sams={params.bams} \
            --region_file={input.bed} \
            --threads=1 \
            {params.options} \
            --output={wildcards.chr}_{wildcards.dataset}) > {log} 2>&1
        
        bcftools concat \
            -o {output} \
            {wildcards.chr}_{wildcards.dataset}/*/*.vcf.gz
        
        rm -r {wildcards.chr}_{wildcards.dataset}
        """   