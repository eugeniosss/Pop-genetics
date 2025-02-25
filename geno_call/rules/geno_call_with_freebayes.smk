
rule geno_call_with_freebayes:
    input:
        reference=config["ref_genome"],  # Reference genome
        bed_file="split_bed/{chr}.bed"   # Per-chromosome BED file
    output:
        temp("freebayes_{dataset}/{chr}_variants.vcf")   # Output VCF per chromosome
    log:
        "freebayes_{dataset}/{chr}.log"
    conda:
        config["dir"] + "envs/geno_callers_env.yml"
    threads: 1
    params:
        bams=lambda wildcards: f"{wildcards.dataset}.txt",  # Ensure dataset is a wildcard
        options=config["options"]["freebayes"]
    shell:
        """
        mkdir -p freebayes
        (freebayes \
            --fasta {input.reference} \
            -L {params.bams} \
            {params.options} \
            -t {input.bed_file} > {output}) > {log} 2>&1
        """
