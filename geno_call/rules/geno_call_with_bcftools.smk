
rule geno_call_with_bcftools:
    input:
        reference=config["ref_genome"],  # Reference genome
        bed_file="split_bed/{chr}.bed"   # Per-chromosome BED file
    output:
        temp("bcftools_{dataset}/{chr}_variants.vcf")   # Output VCF per chromosome
    log:
        "bcftools_{dataset}/{chr}.log"
    conda:
        config["dir"] + "envs/geno_callers_env.yml"
    threads: 1
    params:
        bams=lambda wildcards: f"{wildcards.dataset}.txt",  # Ensure dataset is a wildcard
        options=config["options"]["bcftools"]
    shell:
        """
        mkdir -p bcftools 
        (bcftools mpileup \
            -f {input.reference} \
            -b {params.bams} \
            {params.options} \
            -T {input.bed_file} | bcftools call -mv -Ov -o {output}) > {log} 2>&1
        """
