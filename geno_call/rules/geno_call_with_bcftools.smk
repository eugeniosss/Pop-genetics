
rule geno_call_with_bcftools:
    input:
        reference=config["ref_genome"],  # Reference genome
        bed_file="split_bed/{chr}.bed"   # Per-chromosome BED file
    output:
        "bcftools/{chr}_variants.vcf"   # Output VCF per chromosome
    log:
        "bcftools/{chr}.log"
    conda:
        config["dir"] + "envs/geno_callers_env.yml"
    threads: 1
    params:
        bams=config["bams"],  # BAM file list (define in config or Snakefile)
    shell:
        "(bcftools mpileup -f {input.reference} -b {params.bams} -T {input.bed_file} | bcftools call -mv -Ov -o {output}) > {log} 2>&1"
